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

/// \file EventSelectionModule.h
/// \brief
/// \author ALICE

#ifndef COMMON_TOOLS_EVENTSELECTIONMODULE_H_
#define COMMON_TOOLS_EVENTSELECTIONMODULE_H_

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/RCTSelectionFlags.h"
#include "Common/CCDB/TriggerAliases.h"
#include "Common/Core/TableHelper.h"
#include "Common/DataModel/EventSelection.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/LHCConstants.h>
#include <DataFormatsCTP/Configuration.h>
#include <DataFormatsCTP/Scalers.h>
#include <DataFormatsFT0/Digit.h>
#include <DataFormatsITSMFT/TimeDeadMap.h>
#include <DataFormatsParameters/AggregatedRunInfo.h>
#include <DataFormatsParameters/GRPLHCIFData.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/Logger.h>
#include <ITSMFTBase/DPLAlpideParam.h>
#include <ITSMFTReconstruction/ChipMappingITS.h>

#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <TString.h>

#include <Rtypes.h>
#include <RtypesCore.h>

#include <algorithm>
#include <bitset>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iterator>
#include <limits>
#include <map>
#include <string>
#include <utility>
#include <vector>

#define bitcheck(var, nbit) ((var) & (static_cast<uint32_t>(1) << (nbit)))
#define bitcheck64(var, nbit) ((var) & (static_cast<uint64_t>(1) << (nbit)))

//__________________________________________
// MultModule

namespace o2
{
namespace common
{
namespace eventselection
{
static const double bcNS = o2::constants::lhc::LHCBunchSpacingNS;
static const int32_t nBCsPerOrbit = o2::constants::lhc::LHCMaxBunches;

// for providing temporary buffer
// FIXME ideally cursors could be readable
// to avoid duplicate memory allocation but ok
struct bcselEntry {
  uint32_t alias{0};
  uint64_t selection{0};
  uint32_t rct{0};
  int foundFT0Id = -1;
  int foundFV0Id = -1;
  int foundFDDId = -1;
  int foundZDCId = -1;
};

// bc selection configurables
struct bcselConfigurables : o2::framework::ConfigurableGroup {
  std::string prefix = "bcselOpts";
  o2::framework::Configurable<int> amIneeded{"amIneeded", -1, "run BC selection or not. -1: automatic; 0: no; 1: yes"};                                                           // o2-linter: disable=name/configurable (temporary fix)
  o2::framework::Configurable<int> confTriggerBcShift{"triggerBcShift", 0, "set either custom shift or 999 for apass2/apass3 in LHC22o-t"};                                       // o2-linter: disable=name/configurable (temporary fix)
  o2::framework::Configurable<int> confITSROFrameStartBorderMargin{"ITSROFrameStartBorderMargin", -1, "Number of bcs at the start of ITS RO Frame border. Take from CCDB if -1"}; // o2-linter: disable=name/configurable (temporary fix)
  o2::framework::Configurable<int> confITSROFrameEndBorderMargin{"ITSROFrameEndBorderMargin", -1, "Number of bcs at the end of ITS RO Frame border. Take from CCDB if -1"};       // o2-linter: disable=name/configurable (temporary fix)
  o2::framework::Configurable<int> confTimeFrameStartBorderMargin{"TimeFrameStartBorderMargin", -1, "Number of bcs to cut at the start of the Time Frame. Take from CCDB if -1"}; // o2-linter: disable=name/configurable (temporary fix)
  o2::framework::Configurable<int> confTimeFrameEndBorderMargin{"TimeFrameEndBorderMargin", -1, "Number of bcs to cut at the end of the Time Frame. Take from CCDB if -1"};       // o2-linter: disable=name/configurable (temporary fix)
  o2::framework::Configurable<bool> confCheckRunDurationLimits{"checkRunDurationLimits", false, "Check if the BCs are within the run duration limits"};                           // o2-linter: disable=name/configurable (temporary fix)
  o2::framework::Configurable<std::vector<int>> maxInactiveChipsPerLayer{"maxInactiveChipsPerLayer", {8, 8, 8, 111, 111, 195, 195}, "Maximum allowed number of inactive ITS chips per layer"};
  o2::framework::Configurable<int> confNumberOfOrbitsPerTF{"NumberOfOrbitsPerTF", -1, "Number of orbits per Time Frame. Take from CCDB if -1"}; // o2-linter: disable=name/configurable (temporary fix)
};

// event selection configurables
struct evselConfigurables : o2::framework::ConfigurableGroup {
  std::string prefix = "evselOpts";
  bool isMC_metadata = false;
  o2::framework::Configurable<int> amIneeded{"amIneeded", -1, "run event selection or not. -1: automatic; 0: no; 1: yes"}; // o2-linter: disable=name/configurable (temporary fix)
  o2::framework::Configurable<int> muonSelection{"muonSelection", 0, "0 - barrel, 1 - muon selection with pileup cuts, 2 - muon selection without pileup cuts"};
  o2::framework::Configurable<float> maxDiffZvtxFT0vsPV{"maxDiffZvtxFT0vsPV", 1., "maximum difference (in cm) between z-vertex from FT0 and PV"};
  o2::framework::Configurable<int> isMC{"isMC", -1, "-1 - autoset, 0 - data, 1 - MC"};
  o2::framework::Configurable<int> confSigmaBCforHighPtTracks{"confSigmaBCforHighPtTracks", 4, "Custom sigma (in bcs) for collisions with high-pt tracks"};

  // configurables for occupancy-based event selection
  o2::framework::Configurable<float> confTimeIntervalForOccupancyCalculationMin{"TimeIntervalForOccupancyCalculationMin", -40, "Min time diff window for TPC occupancy calculation, us"};                        // o2-linter: disable=name/configurable (temporary fix)
  o2::framework::Configurable<float> confTimeIntervalForOccupancyCalculationMax{"TimeIntervalForOccupancyCalculationMax", 100, "Max time diff window for TPC occupancy calculation, us"};                        // o2-linter: disable=name/configurable (temporary fix)
  o2::framework::Configurable<float> confTimeRangeVetoOnCollStrict{"TimeRangeVetoOnCollStrict", 10.0, "Exclusion of a collision if there are other collisions nearby, +/- us"};                                  // o2-linter: disable=name/configurable (temporary fix)
  o2::framework::Configurable<float> confTimeRangeVetoOnCollNarrow{"TimeRangeVetoOnCollNarrow", 0.25, "Exclusion of a collision if other collisions nearby, to suppress bc-collision mis-associations, +/- us"}; // o2-linter: disable=name/configurable (temporary fix)
  o2::framework::Configurable<int> confFT0CamplCutVetoOnCollInTimeRange{"FT0CamplPerCollCutVetoOnCollInTimeRange", 8000, "Max allowed FT0C amplitude for each nearby collision in +/- time range"};              // o2-linter: disable=name/configurable (temporary fix)
  o2::framework::Configurable<float> confFT0CamplCutVetoOnCollInROF{"FT0CamplPerCollCutVetoOnCollInROF", 5000, "Max allowed FT0C amplitude for each nearby collision inside this ITS ROF"};                      // o2-linter: disable=name/configurable (temporary fix)
  o2::framework::Configurable<float> confEpsilonVzDiffVetoInROF{"EpsilonVzDiffVetoInROF", 0.3, "Minumum distance to nearby collisions along z inside this ITS ROF, cm"};                                         // o2-linter: disable=name/configurable (temporary fix)
  o2::framework::Configurable<bool> confUseWeightsForOccupancyVariable{"UseWeightsForOccupancyEstimator", 1, "Use or not the delta-time weights for the occupancy estimator"};                                   // o2-linter: disable=name/configurable (temporary fix)
  o2::framework::Configurable<int> confNumberOfOrbitsPerTF{"NumberOfOrbitsPerTF", -1, "Number of orbits per Time Frame. Take from CCDB if -1"};                                                                  // o2-linter: disable=name/configurable (temporary fix)

  // configurables for light-ion event selection
  o2::framework::Configurable<float> confLightIonsNsigmaOnVzDiff{"VzDiffNsigma", 3.0, "+/- nSigma on vZ difference by FT0 and by tracks"};                  // o2-linter: disable=name/configurable (temporary fix)
  o2::framework::Configurable<float> confLightIonsMarginVzDiff{"VzDiffMargin", 0.2, "margin for +/- nSigma cut on vZ difference by FT0 and by tracks, cm"}; // o2-linter: disable=name/configurable (temporary fix)
};

// luminosity configurables
struct lumiConfigurables : o2::framework::ConfigurableGroup {
  std::string prefix = "lumiOpts";
  o2::framework::Configurable<int> amIneeded{"amIneeded", -1, "run BC selection or not. -1: automatic; 0: no; 1: yes"}; // o2-linter: disable=name/configurable (temporary fix)
};

class BcSelectionModule
{
 public:
  BcSelectionModule()
  {
    // constructor
  }
  // declaration of structs here
  // (N.B.: will be invisible to the outside, create your own copies)
  o2::common::eventselection::bcselConfigurables bcselOpts;

  int lastRun = -1;
  int64_t lastTF = -1;
  uint32_t lastRCT = 0;
  uint64_t sorTimestamp = 0;             // default SOR timestamp
  uint64_t eorTimestamp = 1;             // default EOR timestamp
  int64_t bcSOR = -1;                    // global bc of the start of run
  int64_t nBCsPerTF = -1;                // duration of TF in bcs, should be 128*3564 or 32*3564
  int rofOffset = -1;                    // ITS ROF offset, in bc
  int rofLength = -1;                    // ITS ROF length, in bc
  int mITSROFrameStartBorderMargin = 10; // default value
  int mITSROFrameEndBorderMargin = 20;   // default value
  int mTimeFrameStartBorderMargin = 300; // default value
  int mTimeFrameEndBorderMargin = 4000;  // default value
  std::string strLPMProductionTag = "";  // MC production tag to be retrieved from AO2D metadata

  TriggerAliases* aliases = nullptr;
  EventSelectionParams* par = nullptr;
  std::map<uint64_t, uint32_t>* mapRCT = nullptr;
  std::map<int64_t, std::vector<int16_t>> mapInactiveChips; // number of inactive chips vs orbit per layer
  int64_t prevOrbitForInactiveChips = 0;                    // cached next stored orbit in the inactive chip map
  int64_t nextOrbitForInactiveChips = 0;                    // cached previous stored orbit in the inactive chip map
  bool isGoodITSLayer3 = true;                              // default value
  bool isGoodITSLayer0123 = true;                           // default value
  bool isGoodITSLayersAll = true;                           // default value

  template <typename TContext, typename TBcSelOpts, typename THistoRegistry, typename TMetadataInfo>
  void init(TContext& context, TBcSelOpts const& external_bcselopts, THistoRegistry& histos, TMetadataInfo const& metadataInfo)
  {
    // read in configurations from the task where it's used
    bcselOpts = external_bcselopts;

    if (bcselOpts.amIneeded.value < 0) {
      int bcSelNeeded = -1, evSelNeeded = -1;
      bcselOpts.amIneeded.value = 0;
      enableFlagIfTableRequired(context, "BcSels", bcSelNeeded);
      enableFlagIfTableRequired(context, "EvSels", evSelNeeded);
      if (bcSelNeeded == 1) {
        bcselOpts.amIneeded.value = 1;
        LOGF(info, "BC Selection / Autodetection for aod::BcSels: subscription present, will generate.");
      }
      if (evSelNeeded == 1 && bcSelNeeded == 0) {
        bcselOpts.amIneeded.value = 1;
        LOGF(info, "BC Selection / Autodetection for aod::BcSels: not there, but EvSel needed. Will generate.");
      }
      if (bcSelNeeded == 0 && evSelNeeded == 0) {
        LOGF(info, "BC Selection / Autodetection for aod::BcSels: not required. Skipping generation.");
        return;
      }
    }
    strLPMProductionTag = metadataInfo.get("LPMProductionTag"); // to extract info from ccdb by the tag

    // add counter
    histos.add("bcselection/hCounterInvalidBCTimestamp", "", o2::framework::kTH1D, {{1, 0., 1.}});
  }

  //__________________________________________________
  template <typename TCCDB, typename TBCs>
  bool configure(TCCDB& ccdb, TBCs const& bcs)
  {
    if (bcs.size() == 0) {
      return false;
    }
    int run = bcs.iteratorAt(0).runNumber();
    if (run != lastRun) {
      lastRun = run;
      int run3min = 500000;
      if (run < run3min) {                                  // unanchored Run3 MC
        auto runDuration = ccdb->getRunDuration(run, true); // fatalise if timestamps are not found
        // SOR and EOR timestamps
        sorTimestamp = runDuration.first;  // timestamp of the SOR/SOX/STF in ms
        eorTimestamp = runDuration.second; // timestamp of the EOR/EOX/ETF in ms
        auto ctp = ccdb->template getForTimeStamp<std::vector<int64_t>>("CTP/Calib/OrbitReset", sorTimestamp / 2 + eorTimestamp / 2);
        auto orbitResetMUS = (*ctp)[0];
        // first bc of the first orbit
        bcSOR = static_cast<int64_t>((sorTimestamp * 1000 - orbitResetMUS) / o2::constants::lhc::LHCOrbitMUS) * nBCsPerOrbit;
        // duration of TF in bcs
        nBCsPerTF = 32; // hard-coded for Run3 MC (no info from ccdb at the moment)
      } else {
        auto runInfo = o2::parameters::AggregatedRunInfo::buildAggregatedRunInfo(o2::ccdb::BasicCCDBManager::instance(), run, strLPMProductionTag);
        // SOR and EOR timestamps
        sorTimestamp = runInfo.sor;
        eorTimestamp = runInfo.eor;
        // first bc of the first orbit
        bcSOR = runInfo.orbitSOR * nBCsPerOrbit;
        // duration of TF in bcs
        nBCsPerTF = bcselOpts.confNumberOfOrbitsPerTF < 0 ? runInfo.orbitsPerTF * nBCsPerOrbit : bcselOpts.confNumberOfOrbitsPerTF * nBCsPerOrbit;
      }

      // timestamp of the middle of the run used to access run-wise CCDB entries
      int64_t ts = sorTimestamp / 2 + eorTimestamp / 2;
      // access ITSROF and TF border margins
      par = ccdb->template getForTimeStamp<EventSelectionParams>("EventSelection/EventSelectionParams", ts);
      mITSROFrameStartBorderMargin = bcselOpts.confITSROFrameStartBorderMargin < 0 ? par->fITSROFrameStartBorderMargin : bcselOpts.confITSROFrameStartBorderMargin;
      mITSROFrameEndBorderMargin = bcselOpts.confITSROFrameEndBorderMargin < 0 ? par->fITSROFrameEndBorderMargin : bcselOpts.confITSROFrameEndBorderMargin;
      mTimeFrameStartBorderMargin = bcselOpts.confTimeFrameStartBorderMargin < 0 ? par->fTimeFrameStartBorderMargin : bcselOpts.confTimeFrameStartBorderMargin;
      mTimeFrameEndBorderMargin = bcselOpts.confTimeFrameEndBorderMargin < 0 ? par->fTimeFrameEndBorderMargin : bcselOpts.confTimeFrameEndBorderMargin;
      // ITSROF parameters
      auto alppar = ccdb->template getForTimeStamp<o2::itsmft::DPLAlpideParam<0>>("ITS/Config/AlpideParam", ts);
      rofOffset = alppar->roFrameBiasInBC;
      rofLength = alppar->roFrameLengthInBC;
      // Trigger aliases
      aliases = ccdb->template getForTimeStamp<TriggerAliases>("EventSelection/TriggerAliases", ts);

      // prepare map of inactive chips
      auto itsDeadMap = ccdb->template getForTimeStamp<o2::itsmft::TimeDeadMap>("ITS/Calib/TimeDeadMap", ts);
      auto itsDeadMapOrbits = itsDeadMap->getEvolvingMapKeys(); // roughly every second, ~350 TFs = 350x32 orbits
      std::vector<uint16_t> vClosest;                           // temporary vector of inactive chip ids for the current orbit range
      for (const auto& orbit : itsDeadMapOrbits) {
        itsDeadMap->getMapAtOrbit(orbit, vClosest);
        // insert initial (orbit,vector) pair for each layer
        mapInactiveChips[orbit].resize(o2::itsmft::ChipMappingITS::NLayers, 0);

        // fill map of inactive chips
        for (size_t iel = 0; iel < vClosest.size(); iel++) {
          uint16_t w1 = vClosest[iel];
          bool isLastInSequence = (w1 & 0x8000) == 0;
          uint16_t w2 = isLastInSequence ? w1 + 1 : vClosest[iel + 1];
          uint16_t chipId1 = w1 & 0x7FFF;
          uint16_t chipId2 = w2 & 0x7FFF;
          for (int chipId = chipId1; chipId < chipId2; chipId++) {
            auto layer = o2::itsmft::ChipMappingITS::getLayer(chipId);
            mapInactiveChips[orbit][layer]++;
          }
        } // loop over vector of inactive chip ids
      } // loop over orbits

      // QC info
      std::map<std::string, std::string> metadata;
      metadata["run"] = Form("%d", run);
      ccdb->setFatalWhenNull(0);
      mapRCT = ccdb->template getSpecific<std::map<uint64_t, uint32_t>>("RCT/Flags/RunFlags", ts, metadata);
      ccdb->setFatalWhenNull(1);
      if (mapRCT == nullptr) {
        LOGP(info, "rct object missing... inserting dummy rct flags");
        mapRCT = new std::map<uint64_t, uint32_t>;
        uint32_t dummyValue = 1u << 31; // setting bit 31 to indicate that rct object is missing
        mapRCT->insert(std::pair<uint64_t, uint32_t>(sorTimestamp, dummyValue));
      }
    }
    return true;
  }

  //__________________________________________________
  template <typename TCCDB, typename TBCs, typename TTimestamps, typename TBcSelBuffer, typename TBcSelCursor>
  void processRun2(TCCDB const& ccdb, TBCs const& bcs, TTimestamps const& timestamps, TBcSelBuffer& bcselbuffer, TBcSelCursor& bcsel)
  {
    if (bcselOpts.amIneeded.value == 0) {
      bcselbuffer.clear();
      return;
    }
    bcselbuffer.clear();
    for (const auto& bc : bcs) {
      uint64_t timestamp = timestamps[bc.globalIndex()];
      par = ccdb->template getForTimeStamp<EventSelectionParams>("EventSelection/EventSelectionParams", timestamp);
      aliases = ccdb->template getForTimeStamp<TriggerAliases>("EventSelection/TriggerAliases", timestamp);
      // fill fired aliases
      uint32_t alias{0};
      uint64_t triggerMask = bc.triggerMask();
      for (const auto& al : aliases->GetAliasToTriggerMaskMap()) {
        if (triggerMask & al.second) {
          alias |= BIT(al.first);
        }
      }
      uint64_t triggerMaskNext50 = bc.triggerMaskNext50();
      for (const auto& al : aliases->GetAliasToTriggerMaskNext50Map()) {
        if (triggerMaskNext50 & al.second) {
          alias |= BIT(al.first);
        }
      }
      alias |= BIT(kALL);

      // get timing info from ZDC, FV0, FT0 and FDD
      float timeZNA = bc.has_zdc() ? bc.zdc().timeZNA() : -999.f;
      float timeZNC = bc.has_zdc() ? bc.zdc().timeZNC() : -999.f;
      float timeV0A = bc.has_fv0a() ? bc.fv0a().time() : -999.f;
      float timeV0C = bc.has_fv0c() ? bc.fv0c().time() : -999.f;
      float timeT0A = bc.has_ft0() ? bc.ft0().timeA() : -999.f;
      float timeT0C = bc.has_ft0() ? bc.ft0().timeC() : -999.f;
      float timeFDA = bc.has_fdd() ? bc.fdd().timeA() : -999.f;
      float timeFDC = bc.has_fdd() ? bc.fdd().timeC() : -999.f;

      LOGF(debug, "timeZNA=%f timeZNC=%f", timeZNA, timeZNC);
      LOGF(debug, "timeV0A=%f timeV0C=%f", timeV0A, timeV0C);
      LOGF(debug, "timeFDA=%f timeFDC=%f", timeFDA, timeFDC);
      LOGF(debug, "timeT0A=%f timeT0C=%f", timeT0A, timeT0C);

      // fill time-based selection criteria
      uint64_t selection{0};
      selection |= timeV0A > par->fV0ABBlower && timeV0A < par->fV0ABBupper ? BIT(aod::evsel::kIsBBV0A) : 0;
      selection |= timeV0C > par->fV0CBBlower && timeV0C < par->fV0CBBupper ? BIT(aod::evsel::kIsBBV0C) : 0;
      selection |= timeFDA > par->fFDABBlower && timeFDA < par->fFDABBupper ? BIT(aod::evsel::kIsBBFDA) : 0;
      selection |= timeFDC > par->fFDCBBlower && timeFDC < par->fFDCBBupper ? BIT(aod::evsel::kIsBBFDC) : 0;
      selection |= !(timeV0A > par->fV0ABGlower && timeV0A < par->fV0ABGupper) ? BIT(aod::evsel::kNoBGV0A) : 0;
      selection |= !(timeV0C > par->fV0CBGlower && timeV0C < par->fV0CBGupper) ? BIT(aod::evsel::kNoBGV0C) : 0;
      selection |= !(timeFDA > par->fFDABGlower && timeFDA < par->fFDABGupper) ? BIT(aod::evsel::kNoBGFDA) : 0;
      selection |= !(timeFDC > par->fFDCBGlower && timeFDC < par->fFDCBGupper) ? BIT(aod::evsel::kNoBGFDC) : 0;
      selection |= (timeT0A > par->fT0ABBlower && timeT0A < par->fT0ABBupper) ? BIT(aod::evsel::kIsBBT0A) : 0;
      selection |= (timeT0C > par->fT0CBBlower && timeT0C < par->fT0CBBupper) ? BIT(aod::evsel::kIsBBT0C) : 0;
      selection |= (timeZNA > par->fZNABBlower && timeZNA < par->fZNABBupper) ? BIT(aod::evsel::kIsBBZNA) : 0;
      selection |= (timeZNC > par->fZNCBBlower && timeZNC < par->fZNCBBupper) ? BIT(aod::evsel::kIsBBZNC) : 0;
      selection |= !(std::fabs(timeZNA) > par->fZNABGlower && std::fabs(timeZNA) < par->fZNABGupper) ? BIT(aod::evsel::kNoBGZNA) : 0;
      selection |= !(std::fabs(timeZNC) > par->fZNCBGlower && std::fabs(timeZNC) < par->fZNCBGupper) ? BIT(aod::evsel::kNoBGZNC) : 0;
      selection |= (std::pow((timeZNA + timeZNC - par->fZNSumMean) / par->fZNSumSigma, 2) + std::pow((timeZNA - timeZNC - par->fZNDifMean) / par->fZNDifSigma, 2) < 1) ? BIT(aod::evsel::kIsBBZAC) : 0;

      // Calculate V0 multiplicity per ring
      float multRingV0A[5] = {0.};
      float multRingV0C[4] = {0.};
      float multFV0A = 0;
      float multFV0C = 0;
      if (bc.has_fv0a()) {
        for (unsigned int i = 0; i < bc.fv0a().amplitude().size(); ++i) {
          int ring = bc.fv0a().channel()[i] / 8;
          multRingV0A[ring] += bc.fv0a().amplitude()[i];
          multFV0A += bc.fv0a().amplitude()[i];
        }
      }

      if (bc.has_fv0c()) {
        for (unsigned int i = 0; i < bc.fv0c().amplitude().size(); ++i) {
          int ring = bc.fv0c().channel()[i] / 8;
          multRingV0C[ring] += bc.fv0c().amplitude()[i];
          multFV0C += bc.fv0c().amplitude()[i];
        }
      }

      // Calculate pileup and background related selection flags
      // V0A0 excluded from online V0A charge sum => excluding also from offline sum for consistency
      float ofV0M = multFV0A + multFV0C - multRingV0A[0];
      float onV0M = bc.v0TriggerChargeA() + bc.v0TriggerChargeC();
      float ofSPD = bc.spdFiredChipsL0() + bc.spdFiredChipsL1();
      float onSPD = bc.spdFiredFastOrL0() + bc.spdFiredFastOrL1();
      float multV0C012 = multRingV0C[0] + multRingV0C[1] + multRingV0C[2];

      selection |= (onV0M > par->fV0MOnVsOfA + par->fV0MOnVsOfB * ofV0M) ? BIT(aod::evsel::kNoV0MOnVsOfPileup) : 0;
      selection |= (onSPD > par->fSPDOnVsOfA + par->fSPDOnVsOfB * ofSPD) ? BIT(aod::evsel::kNoSPDOnVsOfPileup) : 0;
      selection |= (multRingV0C[3] > par->fV0CasymA + par->fV0CasymB * multV0C012) ? BIT(aod::evsel::kNoV0Casymmetry) : 0;
      selection |= (TESTBIT(selection, aod::evsel::kIsBBV0A) || TESTBIT(selection, aod::evsel::kIsBBV0C) || ofSPD) ? BIT(aod::evsel::kIsINT1) : 0;
      selection |= (bc.has_ft0() ? TESTBIT(bc.ft0().triggerMask(), o2::ft0::Triggers::bitVertex) : 0) ? BIT(aod::evsel::kIsTriggerTVX) : 0;

      // copy remaining selection decisions from eventCuts
      uint32_t eventCuts = bc.eventCuts();
      selection |= (eventCuts & 1 << aod::kTimeRangeCut) ? BIT(aod::evsel::kIsGoodTimeRange) : 0;
      selection |= (eventCuts & 1 << aod::kIncompleteDAQ) ? BIT(aod::evsel::kNoIncompleteDAQ) : 0;
      selection |= !(eventCuts & 1 << aod::kIsTPCLaserWarmUp) ? BIT(aod::evsel::kNoTPCLaserWarmUp) : 0;
      selection |= !(eventCuts & 1 << aod::kIsTPCHVdip) ? BIT(aod::evsel::kNoTPCHVdip) : 0;
      selection |= !(eventCuts & 1 << aod::kIsPileupFromSPD) ? BIT(aod::evsel::kNoPileupFromSPD) : 0;
      selection |= !(eventCuts & 1 << aod::kIsV0PFPileup) ? BIT(aod::evsel::kNoV0PFPileup) : 0;
      selection |= (eventCuts & 1 << aod::kConsistencySPDandTrackVertices) ? BIT(aod::evsel::kNoInconsistentVtx) : 0;
      selection |= (eventCuts & 1 << aod::kPileupInMultBins) ? BIT(aod::evsel::kNoPileupInMultBins) : 0;
      selection |= (eventCuts & 1 << aod::kPileUpMV) ? BIT(aod::evsel::kNoPileupMV) : 0;
      selection |= (eventCuts & 1 << aod::kTPCPileUp) ? BIT(aod::evsel::kNoPileupTPC) : 0;

      int32_t foundFT0 = bc.has_ft0() ? bc.ft0().globalIndex() : -1;
      int32_t foundFV0 = bc.has_fv0a() ? bc.fv0a().globalIndex() : -1;
      int32_t foundFDD = bc.has_fdd() ? bc.fdd().globalIndex() : -1;
      int32_t foundZDC = bc.has_zdc() ? bc.zdc().globalIndex() : -1;

      // // Fill TVX (T0 vertex) counters FIXME - this is a bug, there's no hCounterTVX in this task
      // if (TESTBIT(selection, aod::evsel::kIsTriggerTVX)) {
      //   histos.get<TH1>(HIST("hCounterTVX"))->Fill(Form("%d", bc.runNumber()), 1);
      // }

      uint32_t rct = 0;

      // initialize properties
      o2::common::eventselection::bcselEntry entry;
      entry.alias = alias;
      entry.selection = selection;
      entry.rct = rct;
      entry.foundFT0Id = foundFT0;
      entry.foundFV0Id = foundFV0;
      entry.foundFDDId = foundFDD;
      entry.foundZDCId = foundZDC;
      bcselbuffer.push_back(entry);

      // Fill bc selection columns
      bcsel(alias, selection, rct, foundFT0, foundFV0, foundFDD, foundZDC);
    } // end bc loop
  } // end processRun2

  //__________________________________________________
  template <typename TCCDB, typename THistoRegistry, typename TBCs, typename TTimestamps, typename TBcSelBuffer, typename TBcSelCursor>
  void processRun3(TCCDB const& ccdb, THistoRegistry& histos, TBCs const& bcs, TTimestamps const& timestamps, TBcSelBuffer& bcselbuffer, TBcSelCursor& bcsel)
  {
    if (bcselOpts.amIneeded.value == 0) {
      bcselbuffer.clear();
      return;
    }
    bcselbuffer.clear();
    if (!configure(ccdb, bcs))
      return; // don't do anything in case configuration reported not ok

    int run = bcs.iteratorAt(0).runNumber();
    // map from GlobalBC to BcId needed to find triggerBc
    std::map<uint64_t, int32_t> mapGlobalBCtoBcId;
    for (const auto& bc : bcs) {
      mapGlobalBCtoBcId[bc.globalBC()] = bc.globalIndex();
    }

    int triggerBcShift = bcselOpts.confTriggerBcShift;
    if (bcselOpts.confTriggerBcShift == 999) {                                                                                                                               // o2-linter: disable=magic-number (special shift for early 2022 data)
      triggerBcShift = (run <= 526766 || (run >= 526886 && run <= 527237) || (run >= 527259 && run <= 527518) || run == 527523 || run == 527734 || run >= 534091) ? 0 : 294; // o2-linter: disable=magic-number (magic list of runs)
    }

    // bc loop
    for (auto bc : bcs) { // o2-linter: disable=const-ref-in-for-loop (use bc as nonconst iterator)
      uint64_t timestamp = timestamps[bc.globalIndex()];
      // store rct flags
      uint32_t rct = lastRCT;
      int64_t thisTF = (bc.globalBC() - bcSOR) / nBCsPerTF;
      if (mapRCT != nullptr && thisTF != lastTF) { // skip for unanchored runs; do it once per TF
        auto itrct = mapRCT->upper_bound(timestamp);
        if (itrct != mapRCT->begin())
          itrct--;
        rct = itrct->second;
        LOGP(debug, "sor={} eor={} ts={} rct={}", sorTimestamp, eorTimestamp, timestamp, rct);
        lastRCT = rct;
        lastTF = thisTF;
      }

      uint32_t alias{0};
      // workaround for pp2022 (trigger info is shifted by -294 bcs)
      int32_t triggerBcId = mapGlobalBCtoBcId[bc.globalBC() + triggerBcShift];
      if (triggerBcId && aliases) {
        auto triggerBc = bcs.iteratorAt(triggerBcId);
        uint64_t triggerMask = triggerBc.triggerMask();
        for (const auto& al : aliases->GetAliasToTriggerMaskMap()) {
          if (triggerMask & al.second) {
            alias |= BIT(al.first);
          }
        }
      }
      alias |= BIT(kALL);

      // get timing info from ZDC, FV0, FT0 and FDD
      float timeZNA = bc.has_zdc() ? bc.zdc().timeZNA() : -999.f;
      float timeZNC = bc.has_zdc() ? bc.zdc().timeZNC() : -999.f;
      float timeV0A = bc.has_fv0a() ? bc.fv0a().time() : -999.f;
      float timeT0A = bc.has_ft0() ? bc.ft0().timeA() : -999.f;
      float timeT0C = bc.has_ft0() ? bc.ft0().timeC() : -999.f;
      float timeFDA = bc.has_fdd() ? bc.fdd().timeA() : -999.f;
      float timeFDC = bc.has_fdd() ? bc.fdd().timeC() : -999.f;
      float timeV0ABG = -999.f;
      float timeT0ABG = -999.f;
      float timeT0CBG = -999.f;
      float timeFDABG = -999.f;
      float timeFDCBG = -999.f;

      uint64_t globalBC = bc.globalBC();
      // move to previous bcs to check beam-gas in FT0, FV0 and FDD
      int64_t backwardMoveCount = 0;
      int64_t deltaBC = 6; // up to 6 bcs back
      while (bc.globalBC() + deltaBC >= globalBC) {
        if (bc == bcs.begin()) {
          break;
        }
        --bc;
        backwardMoveCount++;
        int bcDistanceToBeamGasForFT0 = 1;
        int bcDistanceToBeamGasForFDD = 5;
        if (bc.globalBC() + bcDistanceToBeamGasForFT0 == globalBC) {
          timeV0ABG = bc.has_fv0a() ? bc.fv0a().time() : -999.f;
          timeT0ABG = bc.has_ft0() ? bc.ft0().timeA() : -999.f;
          timeT0CBG = bc.has_ft0() ? bc.ft0().timeC() : -999.f;
        }
        if (bc.globalBC() + bcDistanceToBeamGasForFDD == globalBC) {
          timeFDABG = bc.has_fdd() ? bc.fdd().timeA() : -999.f;
          timeFDCBG = bc.has_fdd() ? bc.fdd().timeC() : -999.f;
        }
      }
      // move back to initial position
      bc.moveByIndex(backwardMoveCount);

      // fill time-based selection criteria
      uint64_t selection{0};
      selection |= timeV0A > par->fV0ABBlower && timeV0A < par->fV0ABBupper ? BIT(aod::evsel::kIsBBV0A) : 0;
      selection |= timeFDA > par->fFDABBlower && timeFDA < par->fFDABBupper ? BIT(aod::evsel::kIsBBFDA) : 0;
      selection |= timeFDC > par->fFDCBBlower && timeFDC < par->fFDCBBupper ? BIT(aod::evsel::kIsBBFDC) : 0;
      selection |= !(timeV0ABG > par->fV0ABGlower && timeV0ABG < par->fV0ABGupper) ? BIT(aod::evsel::kNoBGV0A) : 0;
      selection |= !(timeFDABG > par->fFDABGlower && timeFDABG < par->fFDABGupper) ? BIT(aod::evsel::kNoBGFDA) : 0;
      selection |= !(timeFDCBG > par->fFDCBGlower && timeFDCBG < par->fFDCBGupper) ? BIT(aod::evsel::kNoBGFDC) : 0;
      selection |= !(timeT0ABG > par->fT0ABGlower && timeT0ABG < par->fT0ABGupper) ? BIT(aod::evsel::kNoBGT0A) : 0;
      selection |= !(timeT0CBG > par->fT0CBGlower && timeT0CBG < par->fT0CBGupper) ? BIT(aod::evsel::kNoBGT0C) : 0;
      selection |= (timeT0A > par->fT0ABBlower && timeT0A < par->fT0ABBupper) ? BIT(aod::evsel::kIsBBT0A) : 0;
      selection |= (timeT0C > par->fT0CBBlower && timeT0C < par->fT0CBBupper) ? BIT(aod::evsel::kIsBBT0C) : 0;
      selection |= (timeZNA > par->fZNABBlower && timeZNA < par->fZNABBupper) ? BIT(aod::evsel::kIsBBZNA) : 0;
      selection |= (timeZNC > par->fZNCBBlower && timeZNC < par->fZNCBBupper) ? BIT(aod::evsel::kIsBBZNC) : 0;
      selection |= (std::pow((timeZNA + timeZNC - par->fZNSumMean) / par->fZNSumSigma, 2) + std::pow((timeZNA - timeZNC - par->fZNDifMean) / par->fZNDifSigma, 2) < 1) ? BIT(aod::evsel::kIsBBZAC) : 0;
      selection |= !(std::fabs(timeZNA) > par->fZNABGlower && std::fabs(timeZNA) < par->fZNABGupper) ? BIT(aod::evsel::kNoBGZNA) : 0;
      selection |= !(std::fabs(timeZNC) > par->fZNCBGlower && std::fabs(timeZNC) < par->fZNCBGupper) ? BIT(aod::evsel::kNoBGZNC) : 0;
      selection |= (bc.has_ft0() ? (bc.ft0().triggerMask() & BIT(o2::ft0::Triggers::bitVertex)) > 0 : 0) ? BIT(aod::evsel::kIsTriggerTVX) : 0;

      // check if bc is far from start and end of the ITS RO Frame border
      uint16_t bcInITSROF = (globalBC + nBCsPerOrbit - rofOffset) % rofLength;
      LOGP(debug, "bcInITSROF={}", bcInITSROF);
      selection |= bcInITSROF > mITSROFrameStartBorderMargin && bcInITSROF < rofLength - mITSROFrameEndBorderMargin ? BIT(aod::evsel::kNoITSROFrameBorder) : 0;

      // check if bc is far from the Time Frame borders
      int64_t bcInTF = (globalBC - bcSOR) % nBCsPerTF;
      LOGP(debug, "bcInTF={}", bcInTF);
      selection |= bcInTF > mTimeFrameStartBorderMargin && bcInTF < nBCsPerTF - mTimeFrameEndBorderMargin ? BIT(aod::evsel::kNoTimeFrameBorder) : 0;

      // check number of inactive chips and set kIsGoodITSLayer3, kIsGoodITSLayer0123, kIsGoodITSLayersAll flags
      int64_t orbit = globalBC / nBCsPerOrbit;
      if (mapInactiveChips.size() > 0 && (orbit < prevOrbitForInactiveChips || orbit > nextOrbitForInactiveChips)) {
        auto it = mapInactiveChips.upper_bound(orbit);
        bool isEnd = (it == mapInactiveChips.end());
        if (isEnd)
          it--;
        nextOrbitForInactiveChips = isEnd ? orbit : it->first; // setting current orbit in case we reached the end of mapInactiveChips
        auto vNextInactiveChips = it->second;
        if (it != mapInactiveChips.begin() && !isEnd)
          it--;
        prevOrbitForInactiveChips = it->first;
        auto vPrevInactiveChips = it->second;
        LOGP(debug, "orbit: {}, previous orbit: {}, next orbit: {} ", orbit, prevOrbitForInactiveChips, nextOrbitForInactiveChips);
        LOGP(debug, "next inactive chips: {} {} {} {} {} {} {}", vNextInactiveChips[0], vNextInactiveChips[1], vNextInactiveChips[2], vNextInactiveChips[3], vNextInactiveChips[4], vNextInactiveChips[5], vNextInactiveChips[6]);
        LOGP(debug, "prev inactive chips: {} {} {} {} {} {} {}", vPrevInactiveChips[0], vPrevInactiveChips[1], vPrevInactiveChips[2], vPrevInactiveChips[3], vPrevInactiveChips[4], vPrevInactiveChips[5], vPrevInactiveChips[6]);
        isGoodITSLayer3 = vPrevInactiveChips[3] <= bcselOpts.maxInactiveChipsPerLayer->at(3) && vNextInactiveChips[3] <= bcselOpts.maxInactiveChipsPerLayer->at(3);
        isGoodITSLayer0123 = true;
        for (int i = 0; i < 4; i++) { // o2-linter: disable=magic-number (counting first 4 ITS layers)
          isGoodITSLayer0123 &= vPrevInactiveChips[i] <= bcselOpts.maxInactiveChipsPerLayer->at(i) && vNextInactiveChips[i] <= bcselOpts.maxInactiveChipsPerLayer->at(i);
        }
        isGoodITSLayersAll = true;
        for (int i = 0; i < o2::itsmft::ChipMappingITS::NLayers; i++) {
          isGoodITSLayersAll &= vPrevInactiveChips[i] <= bcselOpts.maxInactiveChipsPerLayer->at(i) && vNextInactiveChips[i] <= bcselOpts.maxInactiveChipsPerLayer->at(i);
        }
      }

      selection |= isGoodITSLayer3 ? BIT(aod::evsel::kIsGoodITSLayer3) : 0;
      selection |= isGoodITSLayer0123 ? BIT(aod::evsel::kIsGoodITSLayer0123) : 0;
      selection |= isGoodITSLayersAll ? BIT(aod::evsel::kIsGoodITSLayersAll) : 0;

      // fill found indices
      int32_t foundFT0 = bc.has_ft0() ? bc.ft0().globalIndex() : -1;
      int32_t foundFV0 = bc.has_fv0a() ? bc.fv0a().globalIndex() : -1;
      int32_t foundFDD = bc.has_fdd() ? bc.fdd().globalIndex() : -1;
      int32_t foundZDC = bc.has_zdc() ? bc.zdc().globalIndex() : -1;
      LOGP(debug, "foundFT0={}", foundFT0);

      const char* srun = Form("%d", run);
      if (timestamp < sorTimestamp || timestamp > eorTimestamp) {
        histos.template get<TH1>(HIST("bcselection/hCounterInvalidBCTimestamp"))->Fill(srun, 1);
        if (bcselOpts.confCheckRunDurationLimits.value) {
          LOGF(warn, "Invalid BC timestamp: %d, run: %d, sor: %d, eor: %d", timestamp, run, sorTimestamp, eorTimestamp);
          alias = 0u;
          selection = 0u;
        }
      }

      // initialize properties
      o2::common::eventselection::bcselEntry entry;
      entry.alias = alias;
      entry.selection = selection;
      entry.rct = rct;
      entry.foundFT0Id = foundFT0;
      entry.foundFV0Id = foundFV0;
      entry.foundFDDId = foundFDD;
      entry.foundZDCId = foundZDC;
      bcselbuffer.push_back(entry);

      // Fill bc selection columns
      bcsel(alias, selection, rct, foundFT0, foundFV0, foundFDD, foundZDC);
    } // end bc loop
  } // end processRun3
}; // end BcSelectionModule

class EventSelectionModule
{
 public:
  EventSelectionModule()
  {
    // constructor
  }

  int run3min = 500000;
  int lastRun = -1;                     // last run number (needed to access ccdb only if run!=lastRun)
  std::bitset<nBCsPerOrbit> bcPatternB; // bc pattern of colliding bunches
  std::vector<int> bcsPattern;          // pattern of colliding BCs

  int64_t bcSOR = -1;                   // global bc of the start of the first orbit
  int64_t nBCsPerTF = -1;               // duration of TF in bcs, should be 128*3564 or 32*3564
  int rofOffset = -1;                   // ITS ROF offset, in bc
  int rofLength = -1;                   // ITS ROF length, in bc
  std::string strLPMProductionTag = ""; // MC production tag to be retrieved from AO2D metadata

  // temporary (?) parameterizations for light ion runs
  int runLightIons = -1;
  int runListLightIons[11] = {564356, 564359, 564373, 564374, 564387, 564400, 564414, 564430, 564445, 564468, 564472};
  std::vector<float> diffVzParMean;  // parameterization for mean of diff vZ by FT0 vs by tracks
  std::vector<float> diffVzParSigma; // parameterization for stddev of diff vZ by FT0 vs by tracks

  int32_t findClosest(const int64_t globalBC, const std::map<int64_t, int32_t>& bcs)
  {
    auto it = bcs.lower_bound(globalBC);
    int64_t bc1 = it->first;
    int32_t index1 = it->second;
    if (it != bcs.begin())
      --it;
    int64_t bc2 = it->first;
    int32_t index2 = it->second;
    int64_t dbc1 = std::abs(bc1 - globalBC);
    int64_t dbc2 = std::abs(bc2 - globalBC);
    return (dbc1 <= dbc2) ? index1 : index2;
  }

  // helper function to find median time in the vector of TOF or TRD-track times
  float getMedian(std::vector<float> v)
  {
    int medianIndex = v.size() / 2;
    std::nth_element(v.begin(), v.begin() + medianIndex, v.end());
    return v[medianIndex];
  }

  // helper function to find closest TVX signal in time and in zVtx
  int64_t findBestGlobalBC(int64_t meanBC, int64_t sigmaBC, int32_t nContrib, float zVtxCol, std::map<int64_t, float>& mapGlobalBcVtxZ)
  {
    // protection against
    if (sigmaBC < 1)
      sigmaBC = 1;

    int64_t minBC = meanBC - 3 * sigmaBC;
    int64_t maxBC = meanBC + 3 * sigmaBC;
    // TODO: use ITS ROF bounds to reduce the search range?

    float zVtxSigma = 2.7 * std::pow(nContrib, -0.466) + 0.024;
    zVtxSigma += 1.0; // additional uncertainty due to imperfectections of FT0 time calibration

    auto itMin = mapGlobalBcVtxZ.lower_bound(minBC);
    auto itMax = mapGlobalBcVtxZ.upper_bound(maxBC);

    float bestChi2 = 1e+10;
    int64_t bestGlobalBC = 0;
    for (std::map<int64_t, float>::iterator it = itMin; it != itMax; ++it) {
      float chi2 = std::pow((it->second - zVtxCol) / zVtxSigma, 2) + std::pow(static_cast<float>(it->first - meanBC) / sigmaBC, 2.);
      if (chi2 < bestChi2) {
        bestChi2 = chi2;
        bestGlobalBC = it->first;
      }
    }

    return bestGlobalBC;
  }

  float calcWeightForOccupancy(float dt)
  {
    float wOccup = 0;
    if (dt >= -40 && dt < -5)                     // collisions in the past                    // o2-linter: disable=magic-number
      wOccup = 1. / 1225 * (dt + 40) * (dt + 40); // o2-linter: disable=magic-number
    else if (dt >= -5 && dt < 15)                 // collisions near a given one           // o2-linter: disable=magic-number
      wOccup = 1;
    else if (dt >= 15 && dt < 40)              // collisions from the future            // o2-linter: disable=magic-number
      wOccup = -0.4 / 25 * dt + 1.24;          // o2-linter: disable=magic-number
    else if (dt >= 40 && dt < 100)             // collisions from the distant future   // o2-linter: disable=magic-number
      wOccup = -0.4 / 60 * dt + 0.6 + 0.8 / 3; // o2-linter: disable=magic-number
    return wOccup;
  }

  // declaration of structs here
  // (N.B.: will be invisible to the outside, create your own copies)
  o2::common::eventselection::evselConfigurables evselOpts;

  template <typename TContext, typename TEvSelOpts, typename THistoRegistry, typename TMetadataInfo>
  void init(TContext& context, TEvSelOpts const& external_evselopts, THistoRegistry& histos, TMetadataInfo const& metadataInfo)
  {
    // read in configurations from the task where it's used
    evselOpts = external_evselopts;

    if (evselOpts.amIneeded.value < 0) {
      enableFlagIfTableRequired(context, "EvSels", evselOpts.amIneeded.value);
      if (evselOpts.amIneeded.value == 0) {
        LOGF(info, "Event Selection / Autodetecting for aod::EvSels: not required, won't generate.");
        return;
      } else {
        LOGF(info, "Event Selection / Autodetecting for aod::EvSels: subscription present, will generate.");
      }
    }

    if (metadataInfo.isFullyDefined()) { // Check if the metadata is initialized (only if not forced from the workflow configuration)
      if (evselOpts.isMC == -1) {
        LOGF(info, "Autosetting the MC mode based on metadata (isMC? %i)", metadataInfo.isMC());
        if (metadataInfo.isMC()) {
          evselOpts.isMC.value = 1;
        } else {
          evselOpts.isMC.value = 0;
        }
      }
    }
    strLPMProductionTag = metadataInfo.get("LPMProductionTag"); // to extract info from ccdb by the tag

    histos.add("eventselection/hColCounterAll", "", framework::kTH1D, {{1, 0., 1.}});
    histos.add("eventselection/hColCounterTVX", "", framework::kTH1D, {{1, 0., 1.}});
    histos.add("eventselection/hColCounterAcc", "", framework::kTH1D, {{1, 0., 1.}});
  }

  //__________________________________________________
  template <typename TCCDB, typename TTimestamps, typename TBCs>
  bool configure(TCCDB& ccdb, TTimestamps const& timestamps, TBCs const& bcs)
  {
    if (bcs.size() == 0) {
      return false;
    }
    int run = bcs.iteratorAt(0).runNumber();
    // extract bc pattern from CCDB for data or anchored MC only
    if (run != lastRun && run >= run3min) {
      lastRun = run;
      auto runInfo = o2::parameters::AggregatedRunInfo::buildAggregatedRunInfo(o2::ccdb::BasicCCDBManager::instance(), run, strLPMProductionTag);
      // first bc of the first orbit
      bcSOR = runInfo.orbitSOR * nBCsPerOrbit;
      // duration of TF in bcs
      nBCsPerTF = evselOpts.confNumberOfOrbitsPerTF < 0 ? runInfo.orbitsPerTF * nBCsPerOrbit : evselOpts.confNumberOfOrbitsPerTF * nBCsPerOrbit;
      // colliding bc pattern
      int64_t ts = timestamps[0];

      // getForTimeStamp replaced with getSpecific to set metadata to zero
      // avoids crash related to specific run number
      auto grplhcif = ccdb->template getSpecific<o2::parameters::GRPLHCIFData>("GLO/Config/GRPLHCIF", ts);
      bcPatternB = grplhcif->getBunchFilling().getBCPattern();
      bcsPattern = grplhcif->getBunchFilling().getFilledBCs();
      if (runLightIons >= 0) {
        for (uint32_t i = 0; i < bcsPattern.size(); i++)
          LOGP(debug, "bcsPattern: i={} bc={}", i, bcsPattern.at(i));
      }

      // extract ITS ROF parameters
      auto alppar = ccdb->template getForTimeStamp<o2::itsmft::DPLAlpideParam<0>>("ITS/Config/AlpideParam", ts);
      rofOffset = alppar->roFrameBiasInBC;
      rofLength = alppar->roFrameLengthInBC;
      LOGP(debug, "ITS ROF Offset={} ITS ROF Length={}", rofOffset, rofLength);

      // special treatment of light ion runs
      if (lastRun >= 564356 && lastRun <= 564472) {
        for (uint32_t i = 0; i < sizeof(runListLightIons) / sizeof(*runListLightIons); i++) {
          if (runListLightIons[i] == lastRun) {
            runLightIons = lastRun;
            // extract parameterization for diff of vZ by FT0 vs by tracks
            auto parMeans = ccdb->template getForTimeStamp<std::vector<float>>("Users/a/altsybee/diffVzCollVsFTOmeanPar", ts);
            auto parSigmas = ccdb->template getForTimeStamp<std::vector<float>>("Users/a/altsybee/diffVzCollVsFTOsigmaPar", ts);
            diffVzParMean = *parMeans;
            diffVzParSigma = *parSigmas;
            LOGP(info, ">>> special treatment for diffVz for light ion run {}", runLightIons);
            for (int i = 0; i < 5; i++)
              LOGP(info, " mean par {} = {}", i, diffVzParMean[i]);
            for (int i = 0; i < 5; i++)
              LOGP(info, " sigma par {} = {}", i, diffVzParSigma[i]);
            break;
          }
        }
      }

    } // if run != lastRun
    return true;
  }

  //__________________________________________________
  template <typename TCCDB, typename THistoRegistry, typename TCollisions, typename TTracklets, typename TSlicecache, typename TTimestamps, typename TBcSelBuffer, typename TEvselCursor>
  void processRun2(TCCDB const& ccdb, THistoRegistry& histos, TCollisions const& collisions, TTracklets const& tracklets, TSlicecache& cache, TTimestamps const& timestamps, TBcSelBuffer const& bcselbuffer, TEvselCursor& evsel)
  {
    if (evselOpts.amIneeded.value == 0) {
      return; // dummy process
    }
    for (const auto& col : collisions) {
      auto bc = col.template bc_as<soa::Join<aod::BCs, aod::Run2BCInfos, aod::Run2MatchedToBCSparse>>();
      uint64_t timestamp = timestamps[bc.globalIndex()];
      EventSelectionParams* par = ccdb->template getForTimeStamp<EventSelectionParams>("EventSelection/EventSelectionParams", timestamp);
      bool* applySelection = par->getSelection(evselOpts.muonSelection);
      if (evselOpts.isMC == 1) {
        applySelection[aod::evsel::kIsBBZAC] = 0;
        applySelection[aod::evsel::kNoV0MOnVsOfPileup] = 0;
        applySelection[aod::evsel::kNoSPDOnVsOfPileup] = 0;
        applySelection[aod::evsel::kNoV0Casymmetry] = 0;
        applySelection[aod::evsel::kNoV0PFPileup] = 0;
      }

      int32_t foundBC = bc.globalIndex();
      int32_t foundFT0 = bcselbuffer[foundBC].foundFT0Id;
      int32_t foundFV0 = bcselbuffer[foundBC].foundFV0Id;
      int32_t foundFDD = bcselbuffer[foundBC].foundFDDId;
      int32_t foundZDC = bcselbuffer[foundBC].foundZDCId;

      // copy alias decisions from bcsel table
      uint32_t alias = bcselbuffer[foundBC].alias;

      // copy selection decisions from bcsel table
      uint64_t selection = bcselbuffer[foundBC].selection;

      // copy rct flags from bcsel table
      uint32_t rct = bcselbuffer[foundBC].rct;

      // calculate V0C012 multiplicity
      float multRingV0C[4] = {0.};
      if (bc.has_fv0c()) {
        for (unsigned int i = 0; i < bc.fv0c().amplitude().size(); ++i) {
          int ring = bc.fv0c().channel()[i] / 8;
          multRingV0C[ring] += bc.fv0c().amplitude()[i];
        }
      }
      float multV0C012 = multRingV0C[0] + multRingV0C[1] + multRingV0C[2];

      // applying selections depending on the number of tracklets
      auto trackletsGrouped = tracklets.sliceByCached(aod::track::collisionId, col.globalIndex(), cache);
      int nTkl = trackletsGrouped.size();
      int spdClusters = bc.spdClustersL0() + bc.spdClustersL1();

      selection |= (spdClusters < par->fSPDClsVsTklA + nTkl * par->fSPDClsVsTklB) ? BIT(aod::evsel::kNoSPDClsVsTklBG) : 0;
      selection |= !(nTkl < 6 && multV0C012 > par->fV0C012vsTklA + nTkl * par->fV0C012vsTklB) ? BIT(aod::evsel::kNoV0C012vsTklBG) : 0; // o2-linter: disable=magic-number (nTkl dependent parameterization)

      // apply int7-like selections
      bool sel7 = 1;
      for (int i = 0; i < aod::evsel::kNsel; i++) {
        sel7 = sel7 && (applySelection[i] ? TESTBIT(selection, i) : 1);
      }

      // TODO introduce array of sel[0]... sel[8] or similar?
      bool sel8 = bitcheck64(selection, aod::evsel::kIsBBT0A) && bitcheck64(selection, aod::evsel::kIsBBT0C); // TODO apply other cuts for sel8
      bool sel1 = bitcheck64(selection, aod::evsel::kIsINT1);
      sel1 = sel1 && bitcheck64(selection, aod::evsel::kNoBGV0A);
      sel1 = sel1 && bitcheck64(selection, aod::evsel::kNoBGV0C);
      sel1 = sel1 && bitcheck64(selection, aod::evsel::kNoTPCLaserWarmUp);
      sel1 = sel1 && bitcheck64(selection, aod::evsel::kNoTPCHVdip);

      // INT1 (SPDFO>0 | V0A | V0C) minimum bias trigger logic used in pp2010 and pp2011
      bool isINT1period = bc.runNumber() <= 136377 || (bc.runNumber() >= 144871 && bc.runNumber() <= 159582); // o2-linter: disable=magic-number (magic run numbers)

      // fill counters
      if (evselOpts.isMC == 1 || (!isINT1period && bitcheck(alias, kINT7)) || (isINT1period && bitcheck(alias, kINT1))) {
        histos.template get<TH1>(HIST("eventselection/hColCounterAll"))->Fill(Form("%d", bc.runNumber()), 1);
        if ((!isINT1period && sel7) || (isINT1period && sel1)) {
          histos.template get<TH1>(HIST("eventselection/hColCounterAcc"))->Fill(Form("%d", bc.runNumber()), 1);
        }
      }

      evsel(alias, selection, rct, sel7, sel8, foundBC, foundFT0, foundFV0, foundFDD, foundZDC, 0, 0);
    }
  } // end processRun2

  //__________________________________________________
  template <typename TCCDB, typename THistoRegistry, typename TBCs, typename TCollisions, typename TPVTracks, typename TFT0s, typename TSlicecache, typename TTimestamps, typename TBcSelBuffer, typename TEvselCursor>
  void processRun3(TCCDB const& ccdb, THistoRegistry& histos, TBCs const& bcs, TCollisions const& cols, TPVTracks const& pvTracks, TFT0s const& ft0s, TSlicecache& cache, TTimestamps const& timestamps, TBcSelBuffer const& bcselbuffer, TEvselCursor& evsel)
  {
    if (evselOpts.amIneeded.value == 0) {
      return; // dummy process
    }
    if (!configure(ccdb, timestamps, bcs))
      return; // don't do anything in case configuration reported not ok

    int run = bcs.iteratorAt(0).runNumber();
    // create maps from globalBC to bc index for TVX-fired bcs
    // to be used for closest TVX searches
    std::map<int64_t, int32_t> mapGlobalBcWithTVX;
    std::map<int64_t, int32_t> mapGlobalBcWithOrInFT0;
    std::map<int64_t, float> mapGlobalBcVtxZ;
    for (const auto& bc : bcs) {
      int64_t globalBC = bc.globalBC();
      // skip non-colliding bcs for data and anchored runs
      if (run >= run3min && bcPatternB[globalBC % nBCsPerOrbit] == 0) {
        continue;
      }

      if (bc.has_ft0()) {
        mapGlobalBcWithOrInFT0[globalBC] = bc.globalIndex();
      }

      auto selection = bcselbuffer[bc.globalIndex()].selection;
      if (bitcheck64(selection, aod::evsel::kIsTriggerTVX)) {
        mapGlobalBcWithTVX[globalBC] = bc.globalIndex();
        mapGlobalBcVtxZ[globalBC] = bc.has_ft0() ? bc.ft0().posZ() : 0;
      }
    }

    // protection against empty FT0 maps
    if (mapGlobalBcWithTVX.size() == 0) {
      LOGP(error, "FT0 table is empty or corrupted. Filling evsel table with dummy values");
      for (const auto& col : cols) {
        auto bc = col.template bc_as<soa::Join<aod::BCs, aod::Run3MatchedToBCSparse>>();
        int32_t foundBC = bc.globalIndex();
        int32_t foundFT0 = bcselbuffer[bc.globalIndex()].foundFT0Id;
        int32_t foundFV0 = bcselbuffer[bc.globalIndex()].foundFV0Id;
        int32_t foundFDD = bcselbuffer[bc.globalIndex()].foundFDDId;
        int32_t foundZDC = bcselbuffer[bc.globalIndex()].foundZDCId;
        uint32_t rct = 0;
        evsel(bcselbuffer[bc.globalIndex()].alias, bcselbuffer[bc.globalIndex()].selection, rct, kFALSE, kFALSE, foundBC, foundFT0, foundFV0, foundFDD, foundZDC, -1, -1);
      }
      return;
    }
    std::vector<int> vTracksITS567perColl(cols.size(), 0); // counter of tracks per collision for occupancy studies
    std::vector<float> vAmpFT0CperColl(cols.size(), 0);    // amplitude FT0C per collision
    std::vector<float> vCollVz(cols.size(), 0);            // vector with vZ positions for each collision
    std::vector<bool> vIsVertexITSTPC(cols.size(), 0);     // at least one of vertex contributors is ITS-TPC track
    std::vector<bool> vIsVertexTOFmatched(cols.size(), 0); // at least one of vertex contributors is matched to TOF
    std::vector<bool> vIsVertexTRDmatched(cols.size(), 0); // at least one of vertex contributors is matched to TRD

    std::vector<int> vCollisionsPerBc(bcs.size(), 0);          // counter of collisions per found bc for pileup checks
    std::vector<int> vCollisionsPileupPerColl(cols.size(), 0); // counter of pileup in the same bc as a given collision
    std::vector<int64_t> vBCinPatternPerColl(cols.size(), 0);  // found nominal BCs for collisions
    std::vector<int> vFoundBCindex(cols.size(), -1);           // indices of found bcs
    std::vector<int64_t> vFoundGlobalBC(cols.size(), 0);       // global BCs for collisions

    std::vector<bool> vIsVertexTOF(cols.size(), 0);
    std::vector<bool> vIsVertexTRD(cols.size(), 0);
    std::vector<bool> vIsVertexTPC(cols.size(), 0);
    std::vector<bool> vIsVertexHighPtTPC(cols.size(), 0);
    std::vector<int> vNcontributors(cols.size(), 0);
    std::vector<float> vWeightedTimesTPCnoTOFnoTRD(cols.size(), 0);
    std::vector<float> vWeightedSigmaTPCnoTOFnoTRD(cols.size(), 0);

    // temporary vectors to find tracks with median time
    std::vector<float> vTrackTimesTOF;
    std::vector<float> vTrackTimesTRDnoTOF;

    // first loop to match collisions to TVX, also extract other per-collision information for further use
    for (const auto& col : cols) {
      int32_t colIndex = col.globalIndex();
      vCollVz[colIndex] = col.posZ();

      vTrackTimesTOF.clear();
      vTrackTimesTRDnoTOF.clear();
      int nPvTracksTPCnoTOFnoTRD = 0;
      int nPvTracksHighPtTPCnoTOFnoTRD = 0;
      const auto& colPvTracks = pvTracks.sliceByCached(aod::track::collisionId, col.globalIndex(), cache);
      float sumTime = 0, sumW = 0, sumHighPtTime = 0, sumHighPtW = 0;
      for (const auto& track : colPvTracks) {
        float trackTime = track.trackTime();
        if (track.itsNCls() >= 5) // o2-linter: disable=magic-number (indeed counting layers 5 6 7)
          vTracksITS567perColl[colIndex]++;
        if (track.hasTRD())
          vIsVertexTRDmatched[colIndex] = 1;
        if (track.hasTPC())
          vIsVertexITSTPC[colIndex] = 1;
        if (track.hasTOF()) {
          vTrackTimesTOF.push_back(trackTime);
          vIsVertexTOFmatched[colIndex] = 1;
        } else if (track.hasTRD()) {
          vTrackTimesTRDnoTOF.push_back(trackTime);
        } else if (track.hasTPC()) {
          float trackTimeRes = track.trackTimeRes();
          float trackPt = track.pt();
          float w = 1. / (trackTimeRes * trackTimeRes);
          sumTime += trackTime * w;
          sumW += w;
          nPvTracksTPCnoTOFnoTRD++;
          if (trackPt > 1) {
            sumHighPtTime += trackTime * w;
            sumHighPtW += w;
            nPvTracksHighPtTPCnoTOFnoTRD++;
          }
        }
      }
      vWeightedTimesTPCnoTOFnoTRD[colIndex] = sumW > 0 ? sumTime / sumW : 0;
      vWeightedSigmaTPCnoTOFnoTRD[colIndex] = sumW > 0 ? std::sqrt(1. / sumW) : 0;
      vNcontributors[colIndex] = colPvTracks.size();
      int nPvTracksTOF = vTrackTimesTOF.size();
      int nPvTracksTRDnoTOF = vTrackTimesTRDnoTOF.size();
      // collision type
      vIsVertexTOF[colIndex] = nPvTracksTOF > 0;
      vIsVertexTRD[colIndex] = nPvTracksTRDnoTOF > 0;
      vIsVertexTPC[colIndex] = nPvTracksTPCnoTOFnoTRD > 0;
      vIsVertexHighPtTPC[colIndex] = nPvTracksHighPtTPCnoTOFnoTRD > 0;

      // collision-bc association, other bc-related routine
      auto bc = col.template bc_as<soa::Join<aod::BCs, aod::Run3MatchedToBCSparse>>();
      int64_t globalBC = bc.globalBC();
      int64_t foundGlobalBC = 0;
      int32_t foundBCindex = -1;

      // alternative collision-BC matching (currently: test mode, the aim is to improve pileup rejection)
      if (runLightIons >= 0) {
        foundGlobalBC = globalBC;
        // find closest nominal bc in pattern
        for (uint32_t i = 0; i < bcsPattern.size(); i++) {
          int32_t localBC = globalBC % nBCsPerOrbit;
          int32_t bcFromPattern = bcsPattern.at(i);
          int64_t bcDiff = bcFromPattern - localBC;
          if (std::abs(bcDiff) <= 20) {
            foundGlobalBC = (globalBC / nBCsPerOrbit) * nBCsPerOrbit + bcFromPattern;
            break; // the bc in pattern is found
          }
        }

        // matched with TOF --> precise time, match to TVX, but keep the nominal foundGlobalBC from pattern
        if (vIsVertexTOFmatched[colIndex]) {
          std::map<int64_t, int32_t>::iterator it = mapGlobalBcWithTVX.find(foundGlobalBC);
          if (it != mapGlobalBcWithTVX.end()) {
            foundBCindex = it->second;                       // TVX at foundGlobalBC is found
          } else {                                           // check if TVX is in nearby bcs
            it = mapGlobalBcWithTVX.find(foundGlobalBC + 1); // next bc
            if (it != mapGlobalBcWithTVX.end()) {
              // foundGlobalBC += 1;
              foundBCindex = it->second;
            } else {
              it = mapGlobalBcWithTVX.find(foundGlobalBC - 1); // previous bc
              if (it != mapGlobalBcWithTVX.end()) {
                // foundGlobalBC -= 1;
                foundBCindex = it->second;
              } else {
                foundBCindex = bc.globalIndex(); // keep original BC index
              }
            }
          }
        } // end of if TOF-matched vertex
        else { // for non-TOF and low-mult vertices, consider nearby nominal bcs
          int64_t meanBC = globalBC + TMath::Nint(sumHighPtTime / sumHighPtW / bcNS);
          int64_t bestGlobalBC = findBestGlobalBC(meanBC, evselOpts.confSigmaBCforHighPtTracks, vNcontributors[colIndex], col.posZ(), mapGlobalBcVtxZ);
          if (bestGlobalBC > 0) {
            foundGlobalBC = bestGlobalBC;
            // find closest nominal bc in pattern
            for (uint32_t j = 0; j < bcsPattern.size(); j++) {
              int32_t bcFromPatternBest = bcsPattern.at(j);
              int64_t bcDiff = bcFromPatternBest - (bestGlobalBC % nBCsPerOrbit);
              if (std::abs(bcDiff) <= 20) {
                foundGlobalBC = (bestGlobalBC / nBCsPerOrbit) * nBCsPerOrbit + bcFromPatternBest;
                break; // the bc in pattern is found
              }
            }
            foundBCindex = mapGlobalBcWithTVX[bestGlobalBC];
          } else {                           // failed to find a proper TVX with small vZ difference
            foundBCindex = bc.globalIndex(); // keep original BC index
          }
        } // end of non-TOF matched vertices
        //  sanitity check: if BC was not found
        if (foundBCindex == -1) {
          foundBCindex = bc.globalIndex();
        }
        vBCinPatternPerColl[colIndex] = foundGlobalBC;
        // end of alternative coll-BC matching (test)
      } else if (nPvTracksTOF > 0) { // "standard matching":
        // for collisions with TOF tracks:
        // take bc corresponding to TOF track with median time
        int64_t tofGlobalBC = globalBC + TMath::Nint(getMedian(vTrackTimesTOF) / bcNS);
        std::map<int64_t, int32_t>::iterator it = mapGlobalBcWithTVX.find(tofGlobalBC);
        if (it != mapGlobalBcWithTVX.end()) {
          foundGlobalBC = it->first;
          foundBCindex = it->second;
        }
      } else if (nPvTracksTPCnoTOFnoTRD == 0 && nPvTracksTRDnoTOF > 0) {
        // for collisions with TRD tracks but without TOF or ITSTPC-only tracks:
        // take bc corresponding to TRD track with median time
        int64_t trdGlobalBC = globalBC + TMath::Nint(getMedian(vTrackTimesTRDnoTOF) / bcNS);
        std::map<int64_t, int32_t>::iterator it = mapGlobalBcWithTVX.find(trdGlobalBC);
        if (it != mapGlobalBcWithTVX.end()) {
          foundGlobalBC = it->first;
          foundBCindex = it->second;
        }
      } else if (nPvTracksHighPtTPCnoTOFnoTRD > 0) {
        // for collisions with high-pt ITSTPC-nonTOF-nonTRD tracks
        // search in 3*confSigmaBCforHighPtTracks range (3*4 bcs by default)
        int64_t meanBC = globalBC + TMath::Nint(sumHighPtTime / sumHighPtW / bcNS);
        int64_t bestGlobalBC = findBestGlobalBC(meanBC, evselOpts.confSigmaBCforHighPtTracks, vNcontributors[colIndex], col.posZ(), mapGlobalBcVtxZ);
        if (bestGlobalBC > 0) {
          foundGlobalBC = bestGlobalBC;
          foundBCindex = mapGlobalBcWithTVX[bestGlobalBC];
        }
      }

      // fill foundBC indices and global BCs
      // keep current bc if TVX matching failed at this step
      vFoundBCindex[colIndex] = foundBCindex >= 0 ? foundBCindex : bc.globalIndex();
      vFoundGlobalBC[colIndex] = foundGlobalBC > 0 ? foundGlobalBC : globalBC;

      // erase found global BC with TVX from the pool of bcs for the next loop over low-pt TPCnoTOFnoTRD collisions
      if (foundBCindex >= 0)
        mapGlobalBcVtxZ.erase(foundGlobalBC);
    }
    // alternative matching: looking for collisions with the same nominal BC
    if (runLightIons >= 0) {
      for (uint32_t iCol = 0; iCol < vBCinPatternPerColl.size(); iCol++) {
        int64_t foundNominalBC = vBCinPatternPerColl[iCol];
        for (uint32_t jCol = 0; jCol < vBCinPatternPerColl.size(); jCol++) {
          int64_t foundNominalBC2 = vBCinPatternPerColl[jCol];
          if (foundNominalBC2 == foundNominalBC) {
            vCollisionsPileupPerColl[iCol]++;
          }
        }
      }
    } else { // continue standard matching: second loop to match remaining low-pt TPCnoTOFnoTRD collisions
      for (const auto& col : cols) {
        int32_t colIndex = col.globalIndex();
        if (vIsVertexTPC[colIndex] > 0 && vIsVertexTOF[colIndex] == 0 && vIsVertexHighPtTPC[colIndex] == 0) {
          float weightedTime = vWeightedTimesTPCnoTOFnoTRD[colIndex];
          float weightedSigma = vWeightedSigmaTPCnoTOFnoTRD[colIndex];
          auto bc = col.template bc_as<soa::Join<aod::BCs, aod::Run3MatchedToBCSparse>>();
          int64_t globalBC = bc.globalBC();
          int64_t meanBC = globalBC + TMath::Nint(weightedTime / bcNS);
          int64_t sigmaBC = TMath::CeilNint(weightedSigma / bcNS);
          int64_t bestGlobalBC = findBestGlobalBC(meanBC, sigmaBC, vNcontributors[colIndex], col.posZ(), mapGlobalBcVtxZ);
          vFoundGlobalBC[colIndex] = bestGlobalBC > 0 ? bestGlobalBC : globalBC;
          vFoundBCindex[colIndex] = bestGlobalBC > 0 ? mapGlobalBcWithTVX[bestGlobalBC] : bc.globalIndex();
        }
        // fill pileup counter
        vCollisionsPerBc[vFoundBCindex[colIndex]]++;
      }
    }

    // pre-loop for occupancy calculation
    std::vector<bool> vIsFullInfoForOccupancy(cols.size(), 0);               // info for occupancy in +/- windows is available (i.e. a given coll is not too close to the TF borders)
    std::vector<bool> vIsCollAtROFborder(cols.size(), 0);                    // collision is close to ITS ROF border
    std::vector<bool> vIsCollRejectedByTFborderCut(cols.size(), 0);          // helper vector with
    std::vector<bool> vCanHaveAssocCollsWithinLastDriftTime(cols.size(), 0); // to see if for some collisions in the occupancy calc (that are close to TF border) we will switch to FT0C based occupancy estimation

    const float timeWinOccupancyCalcMinNS = evselOpts.confTimeIntervalForOccupancyCalculationMin * 1e3; // ns
    const float timeWinOccupancyCalcMaxNS = evselOpts.confTimeIntervalForOccupancyCalculationMax * 1e3; // ns

    for (const auto& col : cols) {
      int32_t colIndex = col.globalIndex();
      int64_t foundGlobalBC = vFoundGlobalBC[colIndex];
      int bcInTF = (foundGlobalBC - bcSOR) % nBCsPerTF;
      vIsFullInfoForOccupancy[colIndex] = ((bcInTF - 10) * bcNS > -timeWinOccupancyCalcMinNS) && ((nBCsPerTF - 10 - bcInTF) * bcNS > timeWinOccupancyCalcMaxNS) ? true : false; // 10 BCs is margin

      int32_t foundBC = vFoundBCindex[colIndex];
      auto bcselEntry = bcselbuffer[foundBC];
      // check if we are close to ROF or TF borders => N tracks are not reliable, but FT0 can be used for occupancy estimation
      if (!bitcheck64(bcselEntry.selection, aod::evsel::kNoITSROFrameBorder)) {
        vIsCollAtROFborder[colIndex] = true;
      }

      if (!bitcheck64(bcselEntry.selection, aod::evsel::kNoTimeFrameBorder)) {
        vIsCollRejectedByTFborderCut[colIndex] = true;
      }
      if (nBCsPerTF - bcInTF < 4000 * 2) {
        vCanHaveAssocCollsWithinLastDriftTime[colIndex] = true;
      }
    }

    // save indices of collisions for occupancy calculation (both in ROF and in time range)
    std::vector<std::vector<int>> vCollsInSameITSROF;
    std::vector<std::vector<int>> vCollsInPrevITSROF;
    std::vector<std::vector<int>> vCollsInTimeWin;
    std::vector<std::vector<float>> vTimeDeltaForColls; // delta time wrt a given collision
    for (const auto& col : cols) {
      int32_t colIndex = col.globalIndex();
      int64_t foundGlobalBC = vFoundGlobalBC[colIndex];
      auto bcselEntr = bcselbuffer[vFoundBCindex[colIndex]];
      if (bcselEntr.foundFT0Id > -1) {
        // required: explicit ft0s table
        auto foundFT0 = ft0s.rawIteratorAt(bcselEntr.foundFT0Id);
        vAmpFT0CperColl[colIndex] = foundFT0.sumAmpC();
      }

      int64_t tfId = (foundGlobalBC - bcSOR) / nBCsPerTF;
      int64_t rofId = (foundGlobalBC + nBCsPerOrbit - rofOffset) / rofLength;

      // ### for in-ROF occupancy
      std::vector<int> vAssocCollInSameROF;
      // find all collisions in the same ROF before a given collision
      int32_t minColIndex = colIndex - 1;
      while (minColIndex >= 0) {
        int64_t thisBC = vFoundGlobalBC[minColIndex];
        // check if this is still the same TF
        int64_t thisTFid = (thisBC - bcSOR) / nBCsPerTF;
        if (thisTFid != tfId)
          break;
        // int thisRofIdInTF = (thisBC - rofOffset) / rofLength;
        int64_t thisRofId = (thisBC + nBCsPerOrbit - rofOffset) / rofLength;

        // check if we are within the same ROF
        if (thisRofId != rofId)
          break;
        vAssocCollInSameROF.push_back(minColIndex);
        minColIndex--;
      }
      // find all collisions in the same ROF after the current one
      int32_t maxColIndex = colIndex + 1;
      while (maxColIndex < cols.size()) {
        int64_t thisBC = vFoundGlobalBC[maxColIndex];
        int64_t thisTFid = (thisBC - bcSOR) / nBCsPerTF;
        if (thisTFid != tfId)
          break;
        int64_t thisRofId = (thisBC + nBCsPerOrbit - rofOffset) / rofLength;
        if (thisRofId != rofId)
          break;
        vAssocCollInSameROF.push_back(maxColIndex);
        maxColIndex++;
      }
      vCollsInSameITSROF.push_back(vAssocCollInSameROF);

      // ### bookkeep collisions in previous ROF
      std::vector<int> vAssocCollInPrevROF;
      minColIndex = colIndex - 1;
      while (minColIndex >= 0) {
        int64_t thisBC = vFoundGlobalBC[minColIndex];
        // check if this is still the same TF
        int64_t thisTFid = (thisBC - bcSOR) / nBCsPerTF;
        if (thisTFid != tfId)
          break;
        int64_t thisRofId = (thisBC + nBCsPerOrbit - rofOffset) / rofLength;
        if (thisRofId == rofId - 1)
          vAssocCollInPrevROF.push_back(minColIndex);
        else if (thisRofId < rofId - 1)
          break;
        minColIndex--;
      }
      vCollsInPrevITSROF.push_back(vAssocCollInPrevROF);

      // ### for occupancy in time windows
      std::vector<int> vAssocToThisCol;
      std::vector<float> vCollsTimeDeltaWrtGivenColl;
      // find all collisions in time window before the current one
      minColIndex = colIndex - 1;
      while (minColIndex >= 0) {
        int64_t thisBC = vFoundGlobalBC[minColIndex];
        // check if this is still the same TF
        int64_t thisTFid = (thisBC - bcSOR) / nBCsPerTF;
        if (thisTFid != tfId)
          break;
        float dt = (thisBC - foundGlobalBC) * bcNS; // ns
        // check if we are within the chosen time range
        if (dt < timeWinOccupancyCalcMinNS)
          break;
        vAssocToThisCol.push_back(minColIndex);
        vCollsTimeDeltaWrtGivenColl.push_back(dt);
        minColIndex--;
      }
      // find all collisions in time window after the current one
      maxColIndex = colIndex + 1;
      while (maxColIndex < cols.size()) {
        int64_t thisBC = vFoundGlobalBC[maxColIndex];
        int64_t thisTFid = (thisBC - bcSOR) / nBCsPerTF;
        if (thisTFid != tfId)
          break;
        float dt = (thisBC - foundGlobalBC) * bcNS; // ns
        if (dt > timeWinOccupancyCalcMaxNS)
          break;
        vAssocToThisCol.push_back(maxColIndex);
        vCollsTimeDeltaWrtGivenColl.push_back(dt);
        maxColIndex++;
      }
      vCollsInTimeWin.push_back(vAssocToThisCol);
      vTimeDeltaForColls.push_back(vCollsTimeDeltaWrtGivenColl);
    }

    // perform the occupancy calculation per ITS ROF and also in the pre-defined time window
    std::vector<int> vNumTracksITS567inFullTimeWin(cols.size(), 0); // counter of tracks in full time window for occupancy studies (excluding given event)
    std::vector<float> vSumAmpFT0CinFullTimeWin(cols.size(), 0);    // sum of FT0C of tracks in full time window for occupancy studies (excluding given event)

    std::vector<bool> vNoCollInTimeRangeStrict(cols.size(), 0);   // no collisions in a specified time range
    std::vector<bool> vNoCollInTimeRangeNarrow(cols.size(), 0);   // no collisions in a specified time range (narrow)
    std::vector<bool> vNoHighMultCollInTimeRange(cols.size(), 0); // no high-mult collisions in a specified time range

    std::vector<bool> vNoCollInSameRofStrict(cols.size(), 0);      // to veto events with other collisions in the same ITS ROF
    std::vector<bool> vNoCollInSameRofStandard(cols.size(), 0);    // to veto events with other collisions in the same ITS ROF, with per-collision multiplicity above threshold
    std::vector<bool> vNoCollInSameRofWithCloseVz(cols.size(), 0); // to veto events with nearby collisions with close vZ
    std::vector<bool> vNoHighMultCollInPrevRof(cols.size(), 0);    // veto events if FT0C amplitude in previous ITS ROF is above threshold

    for (const auto& col : cols) {
      int32_t colIndex = col.globalIndex();
      float vZ = col.posZ();

      // ### in-ROF occupancy
      std::vector<int> vAssocCollInSameROF = vCollsInSameITSROF[colIndex];
      int nITS567tracksForSameRofVetoStrict = 0;    // to veto events with other collisions in the same ITS ROF
      int nCollsInRofWithFT0CAboveVetoStandard = 0; // to veto events with other collisions in the same ITS ROF, with per-collision multiplicity above threshold
      int nITS567tracksForRofVetoOnCloseVz = 0;     // to veto events with nearby collisions with close vZ
      for (uint32_t iCol = 0; iCol < vAssocCollInSameROF.size(); iCol++) {
        int thisColIndex = vAssocCollInSameROF[iCol];
        nITS567tracksForSameRofVetoStrict += vTracksITS567perColl[thisColIndex];
        if (vAmpFT0CperColl[thisColIndex] > evselOpts.confFT0CamplCutVetoOnCollInROF)
          nCollsInRofWithFT0CAboveVetoStandard++;
        if (std::fabs(vCollVz[thisColIndex] - vZ) < evselOpts.confEpsilonVzDiffVetoInROF)
          nITS567tracksForRofVetoOnCloseVz += vTracksITS567perColl[thisColIndex];
      }
      // in-ROF occupancy flags
      vNoCollInSameRofStrict[colIndex] = (nITS567tracksForSameRofVetoStrict == 0);
      vNoCollInSameRofStandard[colIndex] = (nCollsInRofWithFT0CAboveVetoStandard == 0);
      vNoCollInSameRofWithCloseVz[colIndex] = (nITS567tracksForRofVetoOnCloseVz == 0);

      // ### occupancy in previous ROF
      std::vector<int> vAssocCollInPrevROF = vCollsInPrevITSROF[colIndex];
      float totalFT0amplInPrevROF = 0;
      for (uint32_t iCol = 0; iCol < vAssocCollInPrevROF.size(); iCol++) {
        int thisColIndex = vAssocCollInPrevROF[iCol];
        totalFT0amplInPrevROF += vAmpFT0CperColl[thisColIndex];
      }
      // veto events if FT0C amplitude in previous ITS ROF is above threshold
      vNoHighMultCollInPrevRof[colIndex] = (totalFT0amplInPrevROF < evselOpts.confFT0CamplCutVetoOnCollInROF);

      // ### occupancy in time windows
      std::vector<int> vAssocToThisCol = vCollsInTimeWin[colIndex];
      std::vector<float> vCollsTimeDeltaWrtGivenColl = vTimeDeltaForColls[colIndex];
      int nITS567tracksInFullTimeWindow = 0;
      float sumAmpFT0CInFullTimeWindow = 0;
      int nITS567tracksForVetoNarrow = 0;      // to veto events with nearby collisions (narrow range) with per-collision multiplicity above threshold
      int nITS567tracksForVetoStrict = 0;      // to veto events with nearby collisions
      int nCollsWithFT0CAboveVetoStandard = 0; // to veto events with nearby collisions that have per-collision multiplicity above threshold
      int colIndexFirstRejectedByTFborderCut = -1;
      for (uint32_t iCol = 0; iCol < vAssocToThisCol.size(); iCol++) {
        int thisColIndex = vAssocToThisCol[iCol];
        float dt = vCollsTimeDeltaWrtGivenColl[iCol] / 1e3; // ns -> us

        // check if we are close to ITS ROF borders => N ITS tracks is not reliable, and FT0C ampl can be used for occupancy estimation
        // denominator for vAmpFT0CperColl is the approximate conversion factor b/n FT0C ampl and number of PV tracks after cuts
        int nItsTracksAssocColl = !vIsCollAtROFborder[thisColIndex] ? vTracksITS567perColl[thisColIndex] : vAmpFT0CperColl[thisColIndex] / 10.;
        // counting tracks from other collisions in fixed time windows
        if (std::fabs(dt) < evselOpts.confTimeRangeVetoOnCollNarrow)
          nITS567tracksForVetoNarrow += nItsTracksAssocColl;
        if (std::fabs(dt) < evselOpts.confTimeRangeVetoOnCollStrict)
          nITS567tracksForVetoStrict += nItsTracksAssocColl;

        // veto on high-mult collisions nearby, where artificial structures in the dt-occupancy plots are observed
        if (dt > -4.0 && dt < 2.0 && vAmpFT0CperColl[thisColIndex] > evselOpts.confFT0CamplCutVetoOnCollInTimeRange) { // dt in us // o2-linter: disable=magic-number
          nCollsWithFT0CAboveVetoStandard++;
        }

        // check if we are close to TF borders => N ITS tracks is not reliable, and FT0C ampl will be used for occupancy estimation (a loop below)
        if (vIsCollRejectedByTFborderCut[thisColIndex]) {
          if (colIndexFirstRejectedByTFborderCut == -1)
            colIndexFirstRejectedByTFborderCut = thisColIndex;
          continue;
        }

        // weighted occupancy calc:
        if (vIsFullInfoForOccupancy[colIndex]) {
          float wOccup = 1.;
          if (evselOpts.confUseWeightsForOccupancyVariable) {
            // weighted occupancy
            wOccup = calcWeightForOccupancy(dt);
          }
          nITS567tracksInFullTimeWindow += wOccup * nItsTracksAssocColl;
          sumAmpFT0CInFullTimeWindow += wOccup * vAmpFT0CperColl[thisColIndex];
        }
      }

      // if some associated collisions are close to TF border - take FT0C amplitude instead of nTracks, using BC table
      if (vIsFullInfoForOccupancy[colIndex] && vCanHaveAssocCollsWithinLastDriftTime[colIndex] && colIndexFirstRejectedByTFborderCut >= 0) {
        int64_t foundGlobalBC = vFoundGlobalBC[colIndex];
        int64_t tfId = (foundGlobalBC - bcSOR) / nBCsPerTF;
        std::map<int64_t, int32_t>::iterator it = mapGlobalBcWithTVX.find(vFoundGlobalBC[colIndexFirstRejectedByTFborderCut]);
        while (it != mapGlobalBcWithTVX.end()) {
          int64_t thisFoundGlobalBC = it->first;
          int32_t thisFoundBCindex = it->second;
          auto bc = bcs.iteratorAt(thisFoundBCindex);
          int64_t thisTFid = (bc.globalBC() - bcSOR) / nBCsPerTF;
          if (thisTFid != tfId)
            break;

          float dt = (thisFoundGlobalBC - foundGlobalBC) * bcNS; // ns
          if (dt > timeWinOccupancyCalcMaxNS)
            break;

          float multT0C = -1;
          if (bc.has_ft0()) {
            multT0C = bc.ft0().sumAmpC();
            float wOccup = 1.;
            if (evselOpts.confUseWeightsForOccupancyVariable) {
              wOccup = calcWeightForOccupancy(dt / 1e3); // ns -> us
            }
            if (multT0C > 50.) // multiplicity in TVX is non-negligible, take it into occupancy calc
            {
              nITS567tracksInFullTimeWindow += wOccup * multT0C / 10.;
              sumAmpFT0CInFullTimeWindow += wOccup * multT0C;
            }
          }
          it++;
        }
      }

      // protection against TF borders
      if (!vIsFullInfoForOccupancy[colIndex]) { // occupancy in undefined (too close to TF borders)
        nITS567tracksInFullTimeWindow = -1;
        sumAmpFT0CInFullTimeWindow = -1;
      }

      vNumTracksITS567inFullTimeWin[colIndex] = nITS567tracksInFullTimeWindow; // occupancy by a sum of number of ITS tracks (without a current collision)
      vSumAmpFT0CinFullTimeWin[colIndex] = sumAmpFT0CInFullTimeWindow;         // occupancy by a sum of FT0C amplitudes (without a current collision)
      // occupancy flags based on nearby collisions
      vNoCollInTimeRangeNarrow[colIndex] = (nITS567tracksForVetoNarrow == 0);
      vNoCollInTimeRangeStrict[colIndex] = (nITS567tracksForVetoStrict == 0);
      vNoHighMultCollInTimeRange[colIndex] = (nCollsWithFT0CAboveVetoStandard == 0 && nITS567tracksForVetoNarrow == 0);
    }

    for (const auto& col : cols) {
      int32_t colIndex = col.globalIndex();
      int32_t foundBC = vFoundBCindex[colIndex];
      auto bc = bcs.iteratorAt(foundBC);
      auto bcselEntry = bcselbuffer[foundBC];
      int32_t foundFT0 = bcselEntry.foundFT0Id;
      int32_t foundFV0 = bcselEntry.foundFV0Id;
      int32_t foundFDD = bcselEntry.foundFDDId;
      int32_t foundZDC = bcselEntry.foundZDCId;

      // compare zVtx from FT0 and from PV
      bool isGoodZvtxFT0vsPV = 0;
      if (bcselEntry.foundFT0Id > -1) {
        auto foundFT0Inner = ft0s.rawIteratorAt(bcselEntry.foundFT0Id);
        float diffVz = foundFT0Inner.posZ() - col.posZ();
        if (runLightIons == -1) {
          isGoodZvtxFT0vsPV = std::fabs(diffVz) < evselOpts.maxDiffZvtxFT0vsPV;
        } else { // special treatment of light ion runs
          float multT0A = bc.ft0().sumAmpA();
          float multT0C = bc.ft0().sumAmpC();
          float T0M = multT0A + multT0C;
          // calc mean at this T0 ampl.
          float x = (T0M < 50 ? 50 : T0M);
          double diffMean = diffVzParMean[0] + diffVzParMean[1] * pow(x, diffVzParMean[2]) + diffVzParMean[3] * pow(x, diffVzParMean[4]);
          // calc sigma at this T0 ampl.
          x = (T0M < 20 ? 20 : (T0M > 1.2e4 ? 1.2e4 : T0M));
          double diffSigma = diffVzParSigma[0] + diffVzParSigma[1] * pow(x, diffVzParSigma[2]) + diffVzParSigma[3] * pow(x, diffVzParSigma[4]);
          float nSigma = evselOpts.confLightIonsNsigmaOnVzDiff;
          float margin = evselOpts.confLightIonsMarginVzDiff;
          isGoodZvtxFT0vsPV = (diffVz > diffMean - nSigma * diffSigma - margin && diffVz < diffMean + nSigma * diffSigma + margin);
        }
      }

      // copy alias decisions from bcsel table
      uint32_t alias = bcselEntry.alias;

      // copy selection decisions from bcsel table
      uint64_t selection = bcselbuffer[bc.globalIndex()].selection;
      if (runLightIons >= 0) // for light ions, apply different condition to assign pileup flags
        selection |= vCollisionsPileupPerColl[colIndex] <= 1 ? BIT(aod::evsel::kNoSameBunchPileup) : 0;
      else
        selection |= vCollisionsPerBc[foundBC] <= 1 ? BIT(aod::evsel::kNoSameBunchPileup) : 0;
      selection |= vIsVertexITSTPC[colIndex] ? BIT(aod::evsel::kIsVertexITSTPC) : 0;
      selection |= vIsVertexTOFmatched[colIndex] ? BIT(aod::evsel::kIsVertexTOFmatched) : 0;
      selection |= vIsVertexTRDmatched[colIndex] ? BIT(aod::evsel::kIsVertexTRDmatched) : 0;
      selection |= isGoodZvtxFT0vsPV ? BIT(aod::evsel::kIsGoodZvtxFT0vsPV) : 0;

      // selection bits based on occupancy time pattern
      selection |= vNoCollInTimeRangeNarrow[colIndex] ? BIT(aod::evsel::kNoCollInTimeRangeNarrow) : 0;
      selection |= vNoCollInTimeRangeStrict[colIndex] ? BIT(aod::evsel::kNoCollInTimeRangeStrict) : 0;
      selection |= vNoHighMultCollInTimeRange[colIndex] ? BIT(aod::evsel::kNoCollInTimeRangeStandard) : 0;

      // selection bits based on ITS in-ROF occupancy
      selection |= vNoCollInSameRofStrict[colIndex] ? BIT(aod::evsel::kNoCollInRofStrict) : 0;
      selection |= (vNoCollInSameRofStandard[colIndex] && vNoCollInSameRofWithCloseVz[colIndex]) ? BIT(aod::evsel::kNoCollInRofStandard) : 0;
      selection |= vNoHighMultCollInPrevRof[colIndex] ? BIT(aod::evsel::kNoHighMultCollInPrevRof) : 0;

      // copy rct flags from bcsel table
      uint32_t rct = bcselEntry.rct;

      // apply int7-like selections
      bool sel7 = 0;

      // TODO apply other cuts for sel8
      // TODO introduce sel1 etc?
      // TODO introduce array of sel[0]... sel[8] or similar?
      bool sel8 = bitcheck64(bcselEntry.selection, aod::evsel::kIsTriggerTVX) && bitcheck64(bcselEntry.selection, aod::evsel::kNoTimeFrameBorder) && bitcheck64(bcselEntry.selection, aod::evsel::kNoITSROFrameBorder);

      // fill counters
      histos.template get<TH1>(HIST("eventselection/hColCounterAll"))->Fill(Form("%d", bc.runNumber()), 1);
      if (bitcheck64(bcselEntry.selection, aod::evsel::kIsTriggerTVX)) {
        histos.template get<TH1>(HIST("eventselection/hColCounterTVX"))->Fill(Form("%d", bc.runNumber()), 1);
      }
      if (sel8) {
        histos.template get<TH1>(HIST("eventselection/hColCounterAcc"))->Fill(Form("%d", bc.runNumber()), 1);
      }

      evsel(alias, selection, rct, sel7, sel8, foundBC, foundFT0, foundFV0, foundFDD, foundZDC,
            vNumTracksITS567inFullTimeWin[colIndex], vSumAmpFT0CinFullTimeWin[colIndex]);
    }
  } // end processRun3
}; // end EventSelectionModule

class LumiModule
{
 public:
  LumiModule()
  {
    // constructor
  }

  int lastRun = -1; // last run number (needed to access ccdb only if run!=lastRun)
  float csTVX = -1; // dummy -1 for the visible TVX cross section (in ub) used in lumi accounting
  float csTCE = -1; // dummy -1 for the visible TCE cross section (in ub) used in lumi accounting
  float csZEM = -1; // dummy -1 for the visible ZEM cross section (in ub) used in lumi accounting
  float csZNC = -1; // dummy -1 for the visible ZNC cross section (in ub) used in lumi accounting

  std::vector<int64_t> mOrbits;
  std::vector<double> mPileupCorrectionTVX;
  std::vector<double> mPileupCorrectionTCE;
  std::vector<double> mPileupCorrectionZEM;
  std::vector<double> mPileupCorrectionZNC;

  int64_t minOrbitInRange = std::numeric_limits<int64_t>::max();
  int64_t maxOrbitInRange = 0;
  uint32_t currentOrbitIndex = 0;
  std::bitset<nBCsPerOrbit> bcPatternB; // bc pattern of colliding bunches
  std::vector<o2::aod::rctsel::RCTFlagsChecker> mRCTFlagsCheckers;

  // declaration of structs here
  // (N.B.: will be invisible to the outside, create your own copies)
  o2::common::eventselection::lumiConfigurables lumiOpts;

  template <typename TContext, typename TLumiOpts, typename THistoRegistry>
  void init(TContext& context, TLumiOpts const& external_lumiopts, THistoRegistry& histos)
  {
    lumiOpts = external_lumiopts;

    if (lumiOpts.amIneeded.value < 0) {
      int bcSelNeeded = -1, evSelNeeded = -1;
      lumiOpts.amIneeded.value = 0;
      enableFlagIfTableRequired(context, "BcSels", bcSelNeeded);
      enableFlagIfTableRequired(context, "EvSels", evSelNeeded);
      if (bcSelNeeded == 1) {
        lumiOpts.amIneeded.value = 1;
        LOGF(info, "Luminosity / Autodetection for aod::BcSels: subscription present, will generate.");
      }
      if (evSelNeeded == 1 && bcSelNeeded == 0) {
        lumiOpts.amIneeded.value = 1;
        LOGF(info, "Luminosity / Autodetection for aod::BcSels: not there, but EvSel needed. Will generate.");
      }
      if (bcSelNeeded == 0 && evSelNeeded == 0) {
        LOGF(info, "Luminosity / Autodetection for aod::BcSels: not required. Skipping generation.");
        return;
      }
    }

    histos.add("luminosity/hCounterTVX", "", framework::kTH1D, {{1, 0., 1.}});
    histos.add("luminosity/hCounterTCE", "", framework::kTH1D, {{1, 0., 1.}});
    histos.add("luminosity/hCounterZEM", "", framework::kTH1D, {{1, 0., 1.}});
    histos.add("luminosity/hCounterZNC", "", framework::kTH1D, {{1, 0., 1.}});
    histos.add("luminosity/hCounterTVXafterBCcuts", "", framework::kTH1D, {{1, 0., 1.}});
    histos.add("luminosity/hCounterTCEafterBCcuts", "", framework::kTH1D, {{1, 0., 1.}});
    histos.add("luminosity/hCounterZEMafterBCcuts", "", framework::kTH1D, {{1, 0., 1.}});
    histos.add("luminosity/hCounterZNCafterBCcuts", "", framework::kTH1D, {{1, 0., 1.}});
    histos.add("luminosity/hLumiTVX", ";;Luminosity, 1/#mub", framework::kTH1D, {{1, 0., 1.}});
    histos.add("luminosity/hLumiTCE", ";;Luminosity, 1/#mub", framework::kTH1D, {{1, 0., 1.}});
    histos.add("luminosity/hLumiZEM", ";;Luminosity, 1/#mub", framework::kTH1D, {{1, 0., 1.}});
    histos.add("luminosity/hLumiZNC", ";;Luminosity, 1/#mub", framework::kTH1D, {{1, 0., 1.}});
    histos.add("luminosity/hLumiTVXafterBCcuts", ";;Luminosity, 1/#mub", framework::kTH1D, {{1, 0., 1.}});
    histos.add("luminosity/hLumiTCEafterBCcuts", ";;Luminosity, 1/#mub", framework::kTH1D, {{1, 0., 1.}});
    histos.add("luminosity/hLumiZEMafterBCcuts", ";;Luminosity, 1/#mub", framework::kTH1D, {{1, 0., 1.}});
    histos.add("luminosity/hLumiZNCafterBCcuts", ";;Luminosity, 1/#mub", framework::kTH1D, {{1, 0., 1.}});

    const int nLists = 6;
    TString rctListNames[] = {"CBT", "CBT_hadronPID", "CBT_electronPID", "CBT_calo", "CBT_muon", "CBT_muon_glo"};
    histos.add("luminosity/hLumiTVXafterBCcutsRCT", ";;Luminosity, 1/#mub", framework::kTH2D, {{1, 0., 1.}, {4 * nLists, -0.5, 4. * nLists - 0.5}});
    histos.add("luminosity/hLumiTCEafterBCcutsRCT", ";;Luminosity, 1/#mub", framework::kTH2D, {{1, 0., 1.}, {4 * nLists, -0.5, 4. * nLists - 0.5}});
    histos.add("luminosity/hLumiZEMafterBCcutsRCT", ";;Luminosity, 1/#mub", framework::kTH2D, {{1, 0., 1.}, {4 * nLists, -0.5, 4. * nLists - 0.5}});
    histos.add("luminosity/hLumiZNCafterBCcutsRCT", ";;Luminosity, 1/#mub", framework::kTH2D, {{1, 0., 1.}, {4 * nLists, -0.5, 4. * nLists - 0.5}});

    for (int i = 0; i < nLists; i++) {
      const auto& rctListName = rctListNames[i];
      mRCTFlagsCheckers.emplace_back(rctListName.Data(), false, false); // disable zdc check, disable lim. acc. check
      mRCTFlagsCheckers.emplace_back(rctListName.Data(), false, true);  // disable zdc check, enable lim. acc. check
      mRCTFlagsCheckers.emplace_back(rctListName.Data(), true, false);  // enable zdc check, disable lim. acc. check
      mRCTFlagsCheckers.emplace_back(rctListName.Data(), true, true);   // enable zdc check, enable lim. acc. check
      histos.template get<TH2>(HIST("luminosity/hLumiTVXafterBCcutsRCT"))->GetYaxis()->SetBinLabel(4 * i + 1, rctListName.Data());
      histos.template get<TH2>(HIST("luminosity/hLumiTCEafterBCcutsRCT"))->GetYaxis()->SetBinLabel(4 * i + 1, rctListName.Data());
      histos.template get<TH2>(HIST("luminosity/hLumiZEMafterBCcutsRCT"))->GetYaxis()->SetBinLabel(4 * i + 1, rctListName.Data());
      histos.template get<TH2>(HIST("luminosity/hLumiZNCafterBCcutsRCT"))->GetYaxis()->SetBinLabel(4 * i + 1, rctListName.Data());
      histos.template get<TH2>(HIST("luminosity/hLumiTVXafterBCcutsRCT"))->GetYaxis()->SetBinLabel(4 * i + 2, (rctListName + "_fullacc").Data());
      histos.template get<TH2>(HIST("luminosity/hLumiTCEafterBCcutsRCT"))->GetYaxis()->SetBinLabel(4 * i + 2, (rctListName + "_fullacc").Data());
      histos.template get<TH2>(HIST("luminosity/hLumiZEMafterBCcutsRCT"))->GetYaxis()->SetBinLabel(4 * i + 2, (rctListName + "_fullacc").Data());
      histos.template get<TH2>(HIST("luminosity/hLumiZNCafterBCcutsRCT"))->GetYaxis()->SetBinLabel(4 * i + 2, (rctListName + "_fullacc").Data());
      histos.template get<TH2>(HIST("luminosity/hLumiTVXafterBCcutsRCT"))->GetYaxis()->SetBinLabel(4 * i + 3, (rctListName + "_zdc").Data());
      histos.template get<TH2>(HIST("luminosity/hLumiTCEafterBCcutsRCT"))->GetYaxis()->SetBinLabel(4 * i + 3, (rctListName + "_zdc").Data());
      histos.template get<TH2>(HIST("luminosity/hLumiZEMafterBCcutsRCT"))->GetYaxis()->SetBinLabel(4 * i + 3, (rctListName + "_zdc").Data());
      histos.template get<TH2>(HIST("luminosity/hLumiZNCafterBCcutsRCT"))->GetYaxis()->SetBinLabel(4 * i + 3, (rctListName + "_zdc").Data());
      histos.template get<TH2>(HIST("luminosity/hLumiTVXafterBCcutsRCT"))->GetYaxis()->SetBinLabel(4 * i + 4, (rctListName + "_zdc" + "_fullacc").Data());
      histos.template get<TH2>(HIST("luminosity/hLumiTCEafterBCcutsRCT"))->GetYaxis()->SetBinLabel(4 * i + 4, (rctListName + "_zdc" + "_fullacc").Data());
      histos.template get<TH2>(HIST("luminosity/hLumiZEMafterBCcutsRCT"))->GetYaxis()->SetBinLabel(4 * i + 4, (rctListName + "_zdc" + "_fullacc").Data());
      histos.template get<TH2>(HIST("luminosity/hLumiZNCafterBCcutsRCT"))->GetYaxis()->SetBinLabel(4 * i + 4, (rctListName + "_zdc" + "_fullacc").Data());
    }
  }

  template <typename TCCDB, typename TTimestamps, typename TBCs>
  bool configure(TCCDB& ccdb, TTimestamps const& timestamps, TBCs const& bcs)
  {
    if (bcs.size() == 0) {
      return false;
    }
    int run = bcs.iteratorAt(0).runNumber();
    if (run < 500000) { // o2-linter: disable=magic-number (skip for unanchored MCs)
      return false;
    }
    if (run != lastRun && run >= 520259) { // o2-linter: disable=magic-number (scalers available for runs above 520120)
      lastRun = run;
      int64_t ts = timestamps[0];

      // getting GRP LHCIF object to extract colliding system, energy and colliding bc pattern
      auto grplhcif = ccdb->template getForTimeStamp<parameters::GRPLHCIFData>("GLO/Config/GRPLHCIF", ts);
      int beamZ1 = grplhcif->getBeamZ(constants::lhc::BeamA);
      int beamZ2 = grplhcif->getBeamZ(constants::lhc::BeamC);
      float sqrts = grplhcif->getSqrtS();
      int nCollidingBCs = grplhcif->getBunchFilling().getNBunches();
      bcPatternB = grplhcif->getBunchFilling().getBCPattern();
      LOGP(info, "beamZ1={} beamZ2={} sqrts={}", beamZ1, beamZ2, sqrts);
      // visible cross sections in ub. Using dummy -1 if lumi estimator is not reliable for this colliding system
      csTVX = -1;
      csTCE = -1;
      csZEM = -1;
      csZNC = -1;
      // Temporary workaround to get visible cross section. TODO: store run-by-run visible cross sections in CCDB
      if (beamZ1 == 8 && beamZ2 == 1) {
        csTVX = 0.3874e6; // eff(TVX) = 0.807 (based on LHC25e6f); sigma(INEL)=0.48b; arxiv:2507.05853
      } else if (beamZ1 == 8 && beamZ2 == 8) {
        csTVX = 1.2050e6; // eff(TVX) = 0.886 (based on LHC25e6b); sigma(INEL)=1.36b; arxiv:2507.05853
      } else if (beamZ1 == 10 && beamZ2 == 10) {
        csTVX = 1.5411e6; // eff(TVX) = 0.896 (based on LHC25e6g); sigma(INEL)=1.72b; arxiv:2507.05853
      } else if (beamZ1 == 1 && beamZ2 == 1) {
        if (std::fabs(sqrts - 900.) < 100.) {          // o2-linter: disable=magic-number (TODO store and extract cross sections from ccdb)
          csTVX = 0.0357e6;                            // ub
        } else if (std::fabs(sqrts - 5360.) < 100.) {  // pp-ref     // o2-linter: disable=magic-number (TODO store and extract cross sections from ccdb)
          csTVX = 0.0503e6;                            // ub
        } else if (std::fabs(sqrts - 13600.) < 300.) { // o2-linter: disable=magic-number (TODO store and extract cross sections from ccdb)
          csTVX = 0.0594e6;                            // ub
        } else {
          LOGP(warn, "Cross section for pp @ {} GeV is not defined", sqrts);
        }
      } else if (beamZ1 == 82 && beamZ2 == 82) { // o2-linter: disable=magic-number (PbPb colliding system)
        // see AN: https://alice-notes.web.cern.ch/node/1515
        if (std::fabs(sqrts - 5360) < 20) {           // o2-linter: disable=magic-number (TODO store and extract cross sections from ccdb)
          csZNC = 214.5e6;                            // ub
          csZEM = 415.2e6;                            // ub
          csTCE = 10.36e6;                            // ub
          if (run > 543437 && run < 543514) {         // o2-linter: disable=magic-number (TODO store and extract cross sections from ccdb)
            csTCE = 8.3e6;                            // ub
          } else if (run >= 543514 && run < 545367) { // o2-linter: disable=magic-number (TODO store and extract cross sections from ccdb)
            csTCE = 4.10e6;                           // ub
          } else if (run >= 559544) {                 // o2-linter: disable=magic-number (TODO store and extract cross sections from ccdb)
            csTCE = 3.86e6;                           // ub
          }
        } else {
          LOGP(warn, "Cross section for PbPb @ {} GeV is not defined", sqrts);
        }
      } else {
        LOGP(warn, "Cross section for z={} + z={} @ {} GeV is not defined", beamZ1, beamZ2, sqrts);
      }
      // getting CTP config to extract lumi class indices (used for rate fetching and pileup correction)
      std::map<std::string, std::string> metadata;
      metadata["runNumber"] = std::to_string(run);
      auto config = ccdb->template getSpecific<o2::ctp::CTPConfiguration>("CTP/Config/Config", ts, metadata);
      auto classes = config->getCTPClasses();
      TString lumiClassNameZNC = "C1ZNC-B-NOPF-CRU";
      TString lumiClassNameTCE = "CMTVXTCE-B-NOPF-CRU";
      TString lumiClassNameTVX1 = "MINBIAS_TVX";         // run >= 534467
      TString lumiClassNameTVX2 = "MINBIAS_TVX_NOMASK";  // run >= 534468
      TString lumiClassNameTVX3 = "CMTVX-NONE-NOPF-CRU"; // run >= 534996
      TString lumiClassNameTVX4 = "CMTVX-B-NOPF-CRU";    // run >= 543437

      // find class indices
      int classIdZNC = -1;
      int classIdTCE = -1;
      int classIdTVX = -1;
      for (unsigned int i = 0; i < classes.size(); i++) {
        TString clname = classes[i].name;
        clname.ToUpper();
        // using position (i) in the vector of classes instead of classes[i].getIndex()
        // due to bug or inconsistencies in scaler record and class indices
        if (clname == lumiClassNameZNC)
          classIdZNC = i;
        if (clname == lumiClassNameTCE)
          classIdTCE = i;
        if (clname == lumiClassNameTVX4 || clname == lumiClassNameTVX3 || clname == lumiClassNameTVX2 || clname == lumiClassNameTVX1)
          classIdTVX = i;
      }

      // extract trigger counts from CTP scalers
      auto scalers = ccdb->template getSpecific<ctp::CTPRunScalers>("CTP/Calib/Scalers", ts, metadata);
      scalers->convertRawToO2();
      std::vector<int64_t> mCounterTVX;
      std::vector<int64_t> mCounterTCE;
      std::vector<int64_t> mCounterZNC;
      std::vector<int64_t> mCounterZEM;
      mOrbits.clear();
      for (const auto& record : scalers->getScalerRecordO2()) {
        mOrbits.push_back(record.intRecord.orbit);
        mCounterTVX.push_back(classIdTVX >= 0 ? record.scalers[classIdTVX].lmBefore : 0);
        mCounterTCE.push_back(classIdTCE >= 0 ? record.scalers[classIdTCE].lmBefore : 0);
        if (run >= 543437 && run < 544448 && record.scalersInps.size() >= 26) { // o2-linter: disable=magic-number (ZNC class not defined for this run range)
          mCounterZNC.push_back(record.scalersInps[25]);                        // see ZNC=1ZNC input index in https://indico.cern.ch/event/1153630/contributions/4844362/
        } else {
          mCounterZNC.push_back(classIdZNC >= 0 ? record.scalers[classIdZNC].l1Before : 0);
        }
        // ZEM class not defined, using inputs instead
        uint32_t indexZEM = 24; // see ZEM=1ZED input index in https://indico.cern.ch/event/1153630/contributions/4844362/
        mCounterZEM.push_back(record.scalersInps.size() >= indexZEM + 1 ? record.scalersInps[indexZEM] : 0);
      }

      // calculate pileup corrections
      mPileupCorrectionTVX.clear();
      mPileupCorrectionTCE.clear();
      mPileupCorrectionZEM.clear();
      mPileupCorrectionZNC.clear();
      for (uint32_t i = 0; i < mOrbits.size() - 1; i++) {
        int64_t nOrbits = mOrbits[i + 1] - mOrbits[i];
        if (nOrbits <= 0 || nCollidingBCs == 0)
          continue;
        double perBcRateTVX = static_cast<double>(mCounterTVX[i + 1] - mCounterTVX[i]) / nOrbits / nCollidingBCs;
        double perBcRateTCE = static_cast<double>(mCounterTCE[i + 1] - mCounterTCE[i]) / nOrbits / nCollidingBCs;
        double perBcRateZNC = static_cast<double>(mCounterZNC[i + 1] - mCounterZNC[i]) / nOrbits / nCollidingBCs;
        double perBcRateZEM = static_cast<double>(mCounterZEM[i + 1] - mCounterZEM[i]) / nOrbits / nCollidingBCs;
        double muTVX = (perBcRateTVX < 1 && perBcRateTVX > 1e-10) ? -std::log(1 - perBcRateTVX) : 0;
        double muTCE = (perBcRateTCE < 1 && perBcRateTCE > 1e-10) ? -std::log(1 - perBcRateTCE) : 0;
        double muZNC = (perBcRateZNC < 1 && perBcRateZNC > 1e-10) ? -std::log(1 - perBcRateZNC) : 0;
        double muZEM = (perBcRateZEM < 1 && perBcRateZEM > 1e-10) ? -std::log(1 - perBcRateZEM) : 0;
        LOGP(debug, "orbit={} muTVX={} muTCE={} muZNC={} muZEM={}", mOrbits[i], muTVX, muTCE, muZNC, muZEM);
        mPileupCorrectionTVX.push_back(muTVX > 1e-10 ? muTVX / (1 - std::exp(-muTVX)) : 1);
        mPileupCorrectionTCE.push_back(muTCE > 1e-10 ? muTCE / (1 - std::exp(-muTCE)) : 1);
        mPileupCorrectionZNC.push_back(muZNC > 1e-10 ? muZNC / (1 - std::exp(-muZNC)) : 1);
        mPileupCorrectionZEM.push_back(muZEM > 1e-10 ? muZEM / (1 - std::exp(-muZEM)) : 1);
      }
      // filling last orbit range using previous orbit range
      mPileupCorrectionTVX.push_back(mPileupCorrectionTVX.back());
      mPileupCorrectionTCE.push_back(mPileupCorrectionTCE.back());
      mPileupCorrectionZNC.push_back(mPileupCorrectionZNC.back());
      mPileupCorrectionZEM.push_back(mPileupCorrectionZEM.back());
    } // access ccdb once per run
    return true; // carry on, please
  }

  //__________________________________________________
  template <typename TCCDB, typename THistoRegistry, typename TBCs, typename TTimestamps, typename TBcSelBuffer>
  void process(TCCDB& ccdb, THistoRegistry& histos, TBCs const& bcs, TTimestamps const& timestamps, TBcSelBuffer const& bcselBuffer)
  {
    if (lumiOpts.amIneeded.value == 0) {
      return;
    }

    if (!configure(ccdb, timestamps, bcs))
      return; // don't do anything in case configuration reported not ok

    int run = bcs.iteratorAt(0).runNumber();
    const char* srun = Form("%d", run);

    // processing loop
    for (const auto& bc : bcs) {
      auto selection = bcselBuffer[bc.globalIndex()].selection;
      if (bcPatternB[bc.globalBC() % nBCsPerOrbit] == 0) // skip non-colliding bcs
        continue;

      bool noBorder = TESTBIT(selection, aod::evsel::kNoTimeFrameBorder) && TESTBIT(selection, aod::evsel::kNoITSROFrameBorder);
      bool isTriggerTVX = TESTBIT(selection, aod::evsel::kIsTriggerTVX);
      bool isTriggerTCE = bc.has_ft0() ? (TESTBIT(selection, aod::evsel::kIsTriggerTVX) && TESTBIT(bc.ft0().triggerMask(), o2::ft0::Triggers::bitCen)) : 0;
      bool isTriggerZNA = TESTBIT(selection, aod::evsel::kIsBBZNA);
      bool isTriggerZNC = TESTBIT(selection, aod::evsel::kIsBBZNC);
      bool isTriggerZEM = isTriggerZNA || isTriggerZNC;

      // determine pileup correction
      int64_t orbit = bc.globalBC() / nBCsPerOrbit;
      if ((orbit < minOrbitInRange || orbit > maxOrbitInRange) && mOrbits.size() > 1) {
        auto it = std::lower_bound(mOrbits.begin(), mOrbits.end(), orbit);
        uint32_t nextOrbitIndex = std::distance(mOrbits.begin(), it);
        if (nextOrbitIndex == 0) // if orbit is below stored scaler orbits
          nextOrbitIndex = 1;
        else if (nextOrbitIndex == mOrbits.size()) // if orbit is above stored scaler orbits
          nextOrbitIndex = mOrbits.size() - 1;
        currentOrbitIndex = nextOrbitIndex - 1;
        minOrbitInRange = mOrbits[currentOrbitIndex];
        maxOrbitInRange = mOrbits[nextOrbitIndex];
      }
      double pileupCorrectionTVX = currentOrbitIndex < mPileupCorrectionTVX.size() ? mPileupCorrectionTVX[currentOrbitIndex] : 1.;
      double pileupCorrectionTCE = currentOrbitIndex < mPileupCorrectionTCE.size() ? mPileupCorrectionTCE[currentOrbitIndex] : 1.;
      double pileupCorrectionZNC = currentOrbitIndex < mPileupCorrectionZNC.size() ? mPileupCorrectionZNC[currentOrbitIndex] : 1.;
      double pileupCorrectionZEM = currentOrbitIndex < mPileupCorrectionZEM.size() ? mPileupCorrectionZEM[currentOrbitIndex] : 1.;

      double lumiTVX = 1. / csTVX * pileupCorrectionTVX;
      double lumiTCE = 1. / csTCE * pileupCorrectionTCE;
      double lumiZNC = 1. / csZNC * pileupCorrectionZNC;
      double lumiZEM = 1. / csZEM * pileupCorrectionZEM;

      auto rct = bcselBuffer[bc.globalIndex()].rct;

      if (isTriggerTVX) {
        histos.template get<TH1>(HIST("luminosity/hCounterTVX"))->Fill(srun, 1);
        histos.template get<TH1>(HIST("luminosity/hLumiTVX"))->Fill(srun, lumiTVX);
        if (noBorder) {
          histos.template get<TH1>(HIST("luminosity/hCounterTVXafterBCcuts"))->Fill(srun, 1);
          histos.template get<TH1>(HIST("luminosity/hLumiTVXafterBCcuts"))->Fill(srun, lumiTVX);
          for (size_t i = 0; i < mRCTFlagsCheckers.size(); i++) {
            if ((rct & mRCTFlagsCheckers[i].value()) == 0)
              histos.template get<TH2>(HIST("luminosity/hLumiTVXafterBCcutsRCT"))->Fill(srun, i, lumiTVX);
          }
        }
      }

      if (isTriggerTCE) {
        histos.template get<TH1>(HIST("luminosity/hCounterTCE"))->Fill(srun, 1);
        histos.template get<TH1>(HIST("luminosity/hLumiTCE"))->Fill(srun, lumiTCE);
        if (noBorder) {
          histos.template get<TH1>(HIST("luminosity/hCounterTCEafterBCcuts"))->Fill(srun, 1);
          histos.template get<TH1>(HIST("luminosity/hLumiTCEafterBCcuts"))->Fill(srun, lumiTCE);
          for (size_t i = 0; i < mRCTFlagsCheckers.size(); i++) {
            if ((rct & mRCTFlagsCheckers[i].value()) == 0)
              histos.template get<TH2>(HIST("luminosity/hLumiTCEafterBCcutsRCT"))->Fill(srun, i, lumiTCE);
          }
        }
      }

      if (isTriggerZEM) {
        histos.template get<TH1>(HIST("luminosity/hCounterZEM"))->Fill(srun, 1);
        histos.template get<TH1>(HIST("luminosity/hLumiZEM"))->Fill(srun, lumiZEM);
        if (noBorder) {
          histos.template get<TH1>(HIST("luminosity/hCounterZEMafterBCcuts"))->Fill(srun, 1);
          histos.template get<TH1>(HIST("luminosity/hLumiZEMafterBCcuts"))->Fill(srun, lumiZEM);
          for (size_t i = 0; i < mRCTFlagsCheckers.size(); i++) {
            if ((rct & mRCTFlagsCheckers[i].value()) == 0)
              histos.template get<TH2>(HIST("luminosity/hLumiZEMafterBCcutsRCT"))->Fill(srun, i, lumiZEM);
          }
        }
      }

      if (isTriggerZNC) {
        histos.template get<TH1>(HIST("luminosity/hCounterZNC"))->Fill(srun, 1);
        histos.template get<TH1>(HIST("luminosity/hLumiZNC"))->Fill(srun, lumiZNC);
        if (noBorder) {
          histos.template get<TH1>(HIST("luminosity/hCounterZNCafterBCcuts"))->Fill(srun, 1);
          histos.template get<TH1>(HIST("luminosity/hLumiZNCafterBCcuts"))->Fill(srun, lumiZNC);
          for (size_t i = 0; i < mRCTFlagsCheckers.size(); i++) {
            if ((rct & mRCTFlagsCheckers[i].value()) == 0)
              histos.template get<TH2>(HIST("luminosity/hLumiZNCafterBCcutsRCT"))->Fill(srun, i, lumiZNC);
          }
        }
      }
    } // bcs
  } // process
}; // end LumiModule

} // namespace eventselection
} // namespace common
} // namespace o2

#endif // COMMON_TOOLS_EVENTSELECTIONMODULE_H_
