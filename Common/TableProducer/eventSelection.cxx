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

/// \file eventSelection.cxx
/// \brief Event selection task
///
/// \author Evgeny Kryshen <evgeny.kryshen@cern.ch> and Igor Altsybeev <Igor.Altsybeev@cern.ch>

#include "Common/DataModel/EventSelection.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/RCTSelectionFlags.h"
#include "Common/CCDB/TriggerAliases.h"
#include "Common/Core/MetadataHelper.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/LHCConstants.h>
#include <DataFormatsCTP/Configuration.h>
#include <DataFormatsCTP/Scalers.h>
#include <DataFormatsFT0/Digit.h>
#include <DataFormatsITSMFT/TimeDeadMap.h>
#include <DataFormatsParameters/AggregatedRunInfo.h>
#include <DataFormatsParameters/GRPLHCIFData.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/DataTypes.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>
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
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iterator>
#include <limits>
#include <map>
#include <string>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::aod::evsel;

o2::common::core::MetadataHelper metadataInfo; // Metadata helper

using BCsWithRun2InfosTimestampsAndMatches = soa::Join<aod::BCs, aod::Run2BCInfos, aod::Timestamps, aod::Run2MatchedToBCSparse>;
using BCsWithRun3Matchings = soa::Join<aod::BCs, aod::Timestamps, aod::Run3MatchedToBCSparse>;
using BCsWithBcSelsRun2 = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels, aod::Run2BCInfos, aod::Run2MatchedToBCSparse>;
using BCsWithBcSelsRun3 = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;
using FullTracks = soa::Join<aod::Tracks, aod::TracksExtra>;
using FullTracksIU = soa::Join<aod::TracksIU, aod::TracksExtra>;
static const double bcNS = o2::constants::lhc::LHCBunchSpacingNS;
static const int32_t nBCsPerOrbit = o2::constants::lhc::LHCMaxBunches;

struct BcSelectionTask {
  Produces<aod::BcSels> bcsel;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Configurable<int> confTriggerBcShift{"triggerBcShift", 0, "set either custom shift or 999 for apass2/apass3 in LHC22o-t"};                                       // o2-linter: disable=name/configurable (temporary fix)
  Configurable<int> confITSROFrameStartBorderMargin{"ITSROFrameStartBorderMargin", -1, "Number of bcs at the start of ITS RO Frame border. Take from CCDB if -1"}; // o2-linter: disable=name/configurable (temporary fix)
  Configurable<int> confITSROFrameEndBorderMargin{"ITSROFrameEndBorderMargin", -1, "Number of bcs at the end of ITS RO Frame border. Take from CCDB if -1"};       // o2-linter: disable=name/configurable (temporary fix)
  Configurable<int> confTimeFrameStartBorderMargin{"TimeFrameStartBorderMargin", -1, "Number of bcs to cut at the start of the Time Frame. Take from CCDB if -1"}; // o2-linter: disable=name/configurable (temporary fix)
  Configurable<int> confTimeFrameEndBorderMargin{"TimeFrameEndBorderMargin", -1, "Number of bcs to cut at the end of the Time Frame. Take from CCDB if -1"};       // o2-linter: disable=name/configurable (temporary fix)
  Configurable<bool> confCheckRunDurationLimits{"checkRunDurationLimits", false, "Check if the BCs are within the run duration limits"};                           // o2-linter: disable=name/configurable (temporary fix)
  Configurable<std::vector<int>> maxInactiveChipsPerLayer{"maxInactiveChipsPerLayer", {8, 8, 8, 111, 111, 195, 195}, "Maximum allowed number of inactive ITS chips per layer"};
  Configurable<int> confNumberOfOrbitsPerTF{"NumberOfOrbitsPerTF", -1, "Number of orbits per Time Frame. Take from CCDB if -1"}; // o2-linter: disable=name/configurable (temporary fix)

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
  void init(InitContext&)
  {
    // ccdb->setURL("http://ccdb-test.cern.ch:8080");
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    strLPMProductionTag = metadataInfo.get("LPMProductionTag"); // to extract info from ccdb by the tag

    histos.add("hCounterInvalidBCTimestamp", "", kTH1D, {{1, 0., 1.}});
  }

  void processRun2(
    BCsWithRun2InfosTimestampsAndMatches const& bcs,
    aod::Zdcs const&,
    aod::FV0As const&,
    aod::FV0Cs const&,
    aod::FT0s const&,
    aod::FDDs const&)
  {
    bcsel.reserve(bcs.size());

    for (const auto& bc : bcs) {
      par = ccdb->getForTimeStamp<EventSelectionParams>("EventSelection/EventSelectionParams", bc.timestamp());
      aliases = ccdb->getForTimeStamp<TriggerAliases>("EventSelection/TriggerAliases", bc.timestamp());
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
      selection |= timeV0A > par->fV0ABBlower && timeV0A < par->fV0ABBupper ? BIT(kIsBBV0A) : 0;
      selection |= timeV0C > par->fV0CBBlower && timeV0C < par->fV0CBBupper ? BIT(kIsBBV0C) : 0;
      selection |= timeFDA > par->fFDABBlower && timeFDA < par->fFDABBupper ? BIT(kIsBBFDA) : 0;
      selection |= timeFDC > par->fFDCBBlower && timeFDC < par->fFDCBBupper ? BIT(kIsBBFDC) : 0;
      selection |= !(timeV0A > par->fV0ABGlower && timeV0A < par->fV0ABGupper) ? BIT(kNoBGV0A) : 0;
      selection |= !(timeV0C > par->fV0CBGlower && timeV0C < par->fV0CBGupper) ? BIT(kNoBGV0C) : 0;
      selection |= !(timeFDA > par->fFDABGlower && timeFDA < par->fFDABGupper) ? BIT(kNoBGFDA) : 0;
      selection |= !(timeFDC > par->fFDCBGlower && timeFDC < par->fFDCBGupper) ? BIT(kNoBGFDC) : 0;
      selection |= (timeT0A > par->fT0ABBlower && timeT0A < par->fT0ABBupper) ? BIT(kIsBBT0A) : 0;
      selection |= (timeT0C > par->fT0CBBlower && timeT0C < par->fT0CBBupper) ? BIT(kIsBBT0C) : 0;
      selection |= (timeZNA > par->fZNABBlower && timeZNA < par->fZNABBupper) ? BIT(kIsBBZNA) : 0;
      selection |= (timeZNC > par->fZNCBBlower && timeZNC < par->fZNCBBupper) ? BIT(kIsBBZNC) : 0;
      selection |= !(std::fabs(timeZNA) > par->fZNABGlower && std::fabs(timeZNA) < par->fZNABGupper) ? BIT(kNoBGZNA) : 0;
      selection |= !(std::fabs(timeZNC) > par->fZNCBGlower && std::fabs(timeZNC) < par->fZNCBGupper) ? BIT(kNoBGZNC) : 0;
      selection |= (std::pow((timeZNA + timeZNC - par->fZNSumMean) / par->fZNSumSigma, 2) + std::pow((timeZNA - timeZNC - par->fZNDifMean) / par->fZNDifSigma, 2) < 1) ? BIT(kIsBBZAC) : 0;

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

      selection |= (onV0M > par->fV0MOnVsOfA + par->fV0MOnVsOfB * ofV0M) ? BIT(kNoV0MOnVsOfPileup) : 0;
      selection |= (onSPD > par->fSPDOnVsOfA + par->fSPDOnVsOfB * ofSPD) ? BIT(kNoSPDOnVsOfPileup) : 0;
      selection |= (multRingV0C[3] > par->fV0CasymA + par->fV0CasymB * multV0C012) ? BIT(kNoV0Casymmetry) : 0;
      selection |= (TESTBIT(selection, kIsBBV0A) || TESTBIT(selection, kIsBBV0C) || ofSPD) ? BIT(kIsINT1) : 0;
      selection |= (bc.has_ft0() ? TESTBIT(bc.ft0().triggerMask(), o2::ft0::Triggers::bitVertex) : 0) ? BIT(kIsTriggerTVX) : 0;

      // copy remaining selection decisions from eventCuts
      uint32_t eventCuts = bc.eventCuts();
      selection |= (eventCuts & 1 << aod::kTimeRangeCut) ? BIT(kIsGoodTimeRange) : 0;
      selection |= (eventCuts & 1 << aod::kIncompleteDAQ) ? BIT(kNoIncompleteDAQ) : 0;
      selection |= !(eventCuts & 1 << aod::kIsTPCLaserWarmUp) ? BIT(kNoTPCLaserWarmUp) : 0;
      selection |= !(eventCuts & 1 << aod::kIsTPCHVdip) ? BIT(kNoTPCHVdip) : 0;
      selection |= !(eventCuts & 1 << aod::kIsPileupFromSPD) ? BIT(kNoPileupFromSPD) : 0;
      selection |= !(eventCuts & 1 << aod::kIsV0PFPileup) ? BIT(kNoV0PFPileup) : 0;
      selection |= (eventCuts & 1 << aod::kConsistencySPDandTrackVertices) ? BIT(kNoInconsistentVtx) : 0;
      selection |= (eventCuts & 1 << aod::kPileupInMultBins) ? BIT(kNoPileupInMultBins) : 0;
      selection |= (eventCuts & 1 << aod::kPileUpMV) ? BIT(kNoPileupMV) : 0;
      selection |= (eventCuts & 1 << aod::kTPCPileUp) ? BIT(kNoPileupTPC) : 0;

      int32_t foundFT0 = bc.has_ft0() ? bc.ft0().globalIndex() : -1;
      int32_t foundFV0 = bc.has_fv0a() ? bc.fv0a().globalIndex() : -1;
      int32_t foundFDD = bc.has_fdd() ? bc.fdd().globalIndex() : -1;
      int32_t foundZDC = bc.has_zdc() ? bc.zdc().globalIndex() : -1;

      // Fill TVX (T0 vertex) counters
      if (TESTBIT(selection, kIsTriggerTVX)) {
        histos.get<TH1>(HIST("hCounterTVX"))->Fill(Form("%d", bc.runNumber()), 1);
      }

      uint32_t rct = 0;
      // Fill bc selection columns
      bcsel(alias, selection, rct, foundFT0, foundFV0, foundFDD, foundZDC);
    }
  }
  PROCESS_SWITCH(BcSelectionTask, processRun2, "Process Run2 event selection", true);

  void processRun3(BCsWithRun3Matchings const& bcs,
                   aod::Zdcs const&,
                   aod::FV0As const&,
                   aod::FT0s const&,
                   aod::FDDs const&)
  {
    if (bcs.size() == 0)
      return;
    bcsel.reserve(bcs.size());

    int run = bcs.iteratorAt(0).runNumber();

    if (run != lastRun) {
      lastRun = run;
      int run3min = 500000;
      if (run < run3min) {                                  // unanchored Run3 MC
        auto runDuration = ccdb->getRunDuration(run, true); // fatalise if timestamps are not found
        // SOR and EOR timestamps
        sorTimestamp = runDuration.first;  // timestamp of the SOR/SOX/STF in ms
        eorTimestamp = runDuration.second; // timestamp of the EOR/EOX/ETF in ms
        auto ctp = ccdb->getForTimeStamp<std::vector<int64_t>>("CTP/Calib/OrbitReset", sorTimestamp / 2 + eorTimestamp / 2);
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
        nBCsPerTF = confNumberOfOrbitsPerTF < 0 ? runInfo.orbitsPerTF * nBCsPerOrbit : static_cast<int64_t>(confNumberOfOrbitsPerTF) * nBCsPerOrbit;
        if (strLPMProductionTag == "LHC25f3") // temporary workaround for MC production LHC25f3 anchored to Pb-Pb 2023 apass5 (to be removed once the info is in ccdb)
          nBCsPerTF = static_cast<int64_t>(8) * nBCsPerOrbit;
      }

      // timestamp of the middle of the run used to access run-wise CCDB entries
      int64_t ts = sorTimestamp / 2 + eorTimestamp / 2;
      // access ITSROF and TF border margins
      par = ccdb->getForTimeStamp<EventSelectionParams>("EventSelection/EventSelectionParams", ts);
      mITSROFrameStartBorderMargin = confITSROFrameStartBorderMargin < 0 ? par->fITSROFrameStartBorderMargin : confITSROFrameStartBorderMargin;
      mITSROFrameEndBorderMargin = confITSROFrameEndBorderMargin < 0 ? par->fITSROFrameEndBorderMargin : confITSROFrameEndBorderMargin;
      mTimeFrameStartBorderMargin = confTimeFrameStartBorderMargin < 0 ? par->fTimeFrameStartBorderMargin : confTimeFrameStartBorderMargin;
      mTimeFrameEndBorderMargin = confTimeFrameEndBorderMargin < 0 ? par->fTimeFrameEndBorderMargin : confTimeFrameEndBorderMargin;
      // ITSROF parameters
      auto alppar = ccdb->getForTimeStamp<o2::itsmft::DPLAlpideParam<0>>("ITS/Config/AlpideParam", ts);
      rofOffset = alppar->roFrameBiasInBC;
      rofLength = alppar->roFrameLengthInBC;
      // Trigger aliases
      aliases = ccdb->getForTimeStamp<TriggerAliases>("EventSelection/TriggerAliases", ts);

      // prepare map of inactive chips
      auto itsDeadMap = ccdb->getForTimeStamp<o2::itsmft::TimeDeadMap>("ITS/Calib/TimeDeadMap", ts);
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
      mapRCT = ccdb->getSpecific<std::map<uint64_t, uint32_t>>("RCT/Flags/RunFlags", ts, metadata);
      ccdb->setFatalWhenNull(1);
      if (mapRCT == nullptr) {
        LOGP(info, "rct object missing... inserting dummy rct flags");
        mapRCT = new std::map<uint64_t, uint32_t>;
        uint32_t dummyValue = 1u << 31; // setting bit 31 to indicate that rct object is missing
        mapRCT->insert(std::pair<uint64_t, uint32_t>(sorTimestamp, dummyValue));
      }
    }

    // map from GlobalBC to BcId needed to find triggerBc
    std::map<uint64_t, int32_t> mapGlobalBCtoBcId;
    for (const auto& bc : bcs) {
      mapGlobalBCtoBcId[bc.globalBC()] = bc.globalIndex();
    }

    int triggerBcShift = confTriggerBcShift;
    if (confTriggerBcShift == 999) {                                                                                                                                         // o2-linter: disable=magic-number (special shift for early 2022 data)
      triggerBcShift = (run <= 526766 || (run >= 526886 && run <= 527237) || (run >= 527259 && run <= 527518) || run == 527523 || run == 527734 || run >= 534091) ? 0 : 294; // o2-linter: disable=magic-number (magic list of runs)
    }

    // bc loop
    for (auto bc : bcs) { // o2-linter: disable=const-ref-in-for-loop (use bc as nonconst iterator)
      // store rct flags
      uint32_t rct = lastRCT;
      int64_t thisTF = (bc.globalBC() - bcSOR) / nBCsPerTF;
      if (mapRCT != nullptr && thisTF != lastTF) { // skip for unanchored runs; do it once per TF
        auto itrct = mapRCT->upper_bound(bc.timestamp());
        if (itrct != mapRCT->begin())
          itrct--;
        rct = itrct->second;
        LOGP(debug, "sor={} eor={} ts={} rct={}", sorTimestamp, eorTimestamp, bc.timestamp(), rct);
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
      selection |= timeV0A > par->fV0ABBlower && timeV0A < par->fV0ABBupper ? BIT(kIsBBV0A) : 0;
      selection |= timeFDA > par->fFDABBlower && timeFDA < par->fFDABBupper ? BIT(kIsBBFDA) : 0;
      selection |= timeFDC > par->fFDCBBlower && timeFDC < par->fFDCBBupper ? BIT(kIsBBFDC) : 0;
      selection |= !(timeV0ABG > par->fV0ABGlower && timeV0ABG < par->fV0ABGupper) ? BIT(kNoBGV0A) : 0;
      selection |= !(timeFDABG > par->fFDABGlower && timeFDABG < par->fFDABGupper) ? BIT(kNoBGFDA) : 0;
      selection |= !(timeFDCBG > par->fFDCBGlower && timeFDCBG < par->fFDCBGupper) ? BIT(kNoBGFDC) : 0;
      selection |= !(timeT0ABG > par->fT0ABGlower && timeT0ABG < par->fT0ABGupper) ? BIT(kNoBGT0A) : 0;
      selection |= !(timeT0CBG > par->fT0CBGlower && timeT0CBG < par->fT0CBGupper) ? BIT(kNoBGT0C) : 0;
      selection |= (timeT0A > par->fT0ABBlower && timeT0A < par->fT0ABBupper) ? BIT(kIsBBT0A) : 0;
      selection |= (timeT0C > par->fT0CBBlower && timeT0C < par->fT0CBBupper) ? BIT(kIsBBT0C) : 0;
      selection |= (timeZNA > par->fZNABBlower && timeZNA < par->fZNABBupper) ? BIT(kIsBBZNA) : 0;
      selection |= (timeZNC > par->fZNCBBlower && timeZNC < par->fZNCBBupper) ? BIT(kIsBBZNC) : 0;
      selection |= (std::pow((timeZNA + timeZNC - par->fZNSumMean) / par->fZNSumSigma, 2) + std::pow((timeZNA - timeZNC - par->fZNDifMean) / par->fZNDifSigma, 2) < 1) ? BIT(kIsBBZAC) : 0;
      selection |= !(std::fabs(timeZNA) > par->fZNABGlower && std::fabs(timeZNA) < par->fZNABGupper) ? BIT(kNoBGZNA) : 0;
      selection |= !(std::fabs(timeZNC) > par->fZNCBGlower && std::fabs(timeZNC) < par->fZNCBGupper) ? BIT(kNoBGZNC) : 0;
      selection |= (bc.has_ft0() ? (bc.ft0().triggerMask() & BIT(o2::ft0::Triggers::bitVertex)) > 0 : 0) ? BIT(kIsTriggerTVX) : 0;

      // check if bc is far from start and end of the ITS RO Frame border
      uint16_t bcInITSROF = (globalBC + nBCsPerOrbit - rofOffset) % rofLength;
      LOGP(debug, "bcInITSROF={}", bcInITSROF);
      selection |= bcInITSROF > mITSROFrameStartBorderMargin && bcInITSROF < rofLength - mITSROFrameEndBorderMargin ? BIT(kNoITSROFrameBorder) : 0;

      // check if bc is far from the Time Frame borders
      int64_t bcInTF = (globalBC - bcSOR) % nBCsPerTF;
      LOGP(debug, "bcInTF={}", bcInTF);
      selection |= bcInTF > mTimeFrameStartBorderMargin && bcInTF < nBCsPerTF - mTimeFrameEndBorderMargin ? BIT(kNoTimeFrameBorder) : 0;

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
        isGoodITSLayer3 = vPrevInactiveChips[3] <= maxInactiveChipsPerLayer->at(3) && vNextInactiveChips[3] <= maxInactiveChipsPerLayer->at(3);
        isGoodITSLayer0123 = true;
        for (int i = 0; i < 4; i++) { // o2-linter: disable=magic-number (counting first 4 ITS layers)
          isGoodITSLayer0123 &= vPrevInactiveChips[i] <= maxInactiveChipsPerLayer->at(i) && vNextInactiveChips[i] <= maxInactiveChipsPerLayer->at(i);
        }
        isGoodITSLayersAll = true;
        for (int i = 0; i < o2::itsmft::ChipMappingITS::NLayers; i++) {
          isGoodITSLayersAll &= vPrevInactiveChips[i] <= maxInactiveChipsPerLayer->at(i) && vNextInactiveChips[i] <= maxInactiveChipsPerLayer->at(i);
        }
      }

      selection |= isGoodITSLayer3 ? BIT(kIsGoodITSLayer3) : 0;
      selection |= isGoodITSLayer0123 ? BIT(kIsGoodITSLayer0123) : 0;
      selection |= isGoodITSLayersAll ? BIT(kIsGoodITSLayersAll) : 0;

      // fill found indices
      int32_t foundFT0 = bc.has_ft0() ? bc.ft0().globalIndex() : -1;
      int32_t foundFV0 = bc.has_fv0a() ? bc.fv0a().globalIndex() : -1;
      int32_t foundFDD = bc.has_fdd() ? bc.fdd().globalIndex() : -1;
      int32_t foundZDC = bc.has_zdc() ? bc.zdc().globalIndex() : -1;
      LOGP(debug, "foundFT0={}", foundFT0);

      const char* srun = Form("%d", run);
      if (bc.timestamp() < sorTimestamp || bc.timestamp() > eorTimestamp) {
        histos.get<TH1>(HIST("hCounterInvalidBCTimestamp"))->Fill(srun, 1);
        if (confCheckRunDurationLimits.value) {
          LOGF(warn, "Invalid BC timestamp: %d, run: %d, sor: %d, eor: %d", bc.timestamp(), run, sorTimestamp, eorTimestamp);
          alias = 0u;
          selection = 0u;
        }
      }

      // Fill bc selection columns
      bcsel(alias, selection, rct, foundFT0, foundFV0, foundFDD, foundZDC);
    }
  }
  PROCESS_SWITCH(BcSelectionTask, processRun3, "Process Run3 event selection", false);
};

struct EventSelectionTask {
  SliceCache cache;
  Produces<aod::EvSels> evsel;
  Configurable<int> muonSelection{"muonSelection", 0, "0 - barrel, 1 - muon selection with pileup cuts, 2 - muon selection without pileup cuts"};
  Configurable<float> maxDiffZvtxFT0vsPV{"maxDiffZvtxFT0vsPV", 1., "maximum difference (in cm) between z-vertex from FT0 and PV"};
  Configurable<int> isMC{"isMC", 0, "-1 - autoset, 0 - data, 1 - MC"};
  Configurable<int> confSigmaBCforHighPtTracks{"confSigmaBCforHighPtTracks", 4, "Custom sigma (in bcs) for collisions with high-pt tracks"};

  // configurables for occupancy-based event selection
  Configurable<float> confTimeIntervalForOccupancyCalculationMin{"TimeIntervalForOccupancyCalculationMin", -40, "Min time diff window for TPC occupancy calculation, us"};           // o2-linter: disable=name/configurable (temporary fix)
  Configurable<float> confTimeIntervalForOccupancyCalculationMax{"TimeIntervalForOccupancyCalculationMax", 100, "Max time diff window for TPC occupancy calculation, us"};           // o2-linter: disable=name/configurable (temporary fix)
  Configurable<float> confTimeRangeVetoOnCollStandard{"TimeRangeVetoOnCollStandard", 10.0, "Exclusion of a collision if there are other collisions nearby, +/- us"};                 // o2-linter: disable=name/configurable (temporary fix)
  Configurable<float> confTimeRangeVetoOnCollNarrow{"TimeRangeVetoOnCollNarrow", 2.0, "Exclusion of a collision if there are other collisions nearby, +/- us"};                      // o2-linter: disable=name/configurable (temporary fix)
  Configurable<int> confFT0CamplCutVetoOnCollInTimeRange{"FT0CamplPerCollCutVetoOnCollInTimeRange", 8000, "Max allowed FT0C amplitude for each nearby collision in +/- time range"}; // o2-linter: disable=name/configurable (temporary fix)
  Configurable<float> confFT0CamplCutVetoOnCollInROF{"FT0CamplPerCollCutVetoOnCollInROF", 5000, "Max allowed FT0C amplitude for each nearby collision inside this ITS ROF"};         // o2-linter: disable=name/configurable (temporary fix)
  Configurable<float> confEpsilonVzDiffVetoInROF{"EpsilonVzDiffVetoInROF", 0.3, "Minumum distance to nearby collisions along z inside this ITS ROF, cm"};                            // o2-linter: disable=name/configurable (temporary fix)
  Configurable<bool> confUseWeightsForOccupancyVariable{"UseWeightsForOccupancyEstimator", 1, "Use or not the delta-time weights for the occupancy estimator"};                      // o2-linter: disable=name/configurable (temporary fix)
  Configurable<int> confNumberOfOrbitsPerTF{"NumberOfOrbitsPerTF", -1, "Number of orbits per Time Frame. Take from CCDB if -1"};                                                     // o2-linter: disable=name/configurable (temporary fix)

  Partition<FullTracks> tracklets = (aod::track::trackType == static_cast<uint8_t>(o2::aod::track::TrackTypeEnum::Run2Tracklet));

  Preslice<FullTracks> perCollision = aod::track::collisionId;
  Preslice<FullTracksIU> perCollisionIU = aod::track::collisionId;

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  int lastRun = -1;                     // last run number (needed to access ccdb only if run!=lastRun)
  std::bitset<nBCsPerOrbit> bcPatternB; // bc pattern of colliding bunches

  int64_t bcSOR = -1;                   // global bc of the start of the first orbit
  int64_t nBCsPerTF = -1;               // duration of TF in bcs, should be 128*3564 or 32*3564
  int rofOffset = -1;                   // ITS ROF offset, in bc
  int rofLength = -1;                   // ITS ROF length, in bc
  std::string strLPMProductionTag = ""; // MC production tag to be retrieved from AO2D metadata

  int32_t findClosest(int64_t globalBC, const std::map<int64_t, int32_t>& bcs)
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

  void init(InitContext&)
  {
    if (metadataInfo.isFullyDefined()) { // Check if the metadata is initialized (only if not forced from the workflow configuration)
      if (isMC == -1) {
        LOG(info) << "Autosetting the MC mode based on metadata";
        if (metadataInfo.isMC()) {
          isMC.value = 1;
        } else {
          isMC.value = 0;
        }
      }
    }
    strLPMProductionTag = metadataInfo.get("LPMProductionTag"); // to extract info from ccdb by the tag

    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    histos.add("hColCounterAll", "", kTH1D, {{1, 0., 1.}});
    histos.add("hColCounterTVX", "", kTH1D, {{1, 0., 1.}});
    histos.add("hColCounterAcc", "", kTH1D, {{1, 0., 1.}});
  }

  void process(aod::Collisions const& collisions)
  {
    evsel.reserve(collisions.size());
  }

  void processRun2(aod::Collision const& col, BCsWithBcSelsRun2 const&, FullTracks const&, aod::FV0Cs const&)
  {
    auto bc = col.bc_as<BCsWithBcSelsRun2>();
    EventSelectionParams* par = ccdb->getForTimeStamp<EventSelectionParams>("EventSelection/EventSelectionParams", bc.timestamp());
    bool* applySelection = par->getSelection(muonSelection);
    if (isMC == 1) {
      applySelection[kIsBBZAC] = 0;
      applySelection[kNoV0MOnVsOfPileup] = 0;
      applySelection[kNoSPDOnVsOfPileup] = 0;
      applySelection[kNoV0Casymmetry] = 0;
      applySelection[kNoV0PFPileup] = 0;
    }

    int32_t foundBC = bc.globalIndex();
    int32_t foundFT0 = bc.foundFT0Id();
    int32_t foundFV0 = bc.foundFV0Id();
    int32_t foundFDD = bc.foundFDDId();
    int32_t foundZDC = bc.foundZDCId();

    // copy alias decisions from bcsel table
    uint32_t alias = bc.alias_raw();

    // copy selection decisions from bcsel table
    uint64_t selection = bc.selection_raw();

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
    auto trackletsGrouped = tracklets->sliceByCached(aod::track::collisionId, col.globalIndex(), cache);
    int nTkl = trackletsGrouped.size();
    int spdClusters = bc.spdClustersL0() + bc.spdClustersL1();

    selection |= (spdClusters < par->fSPDClsVsTklA + nTkl * par->fSPDClsVsTklB) ? BIT(kNoSPDClsVsTklBG) : 0;
    selection |= !(nTkl < 6 && multV0C012 > par->fV0C012vsTklA + nTkl * par->fV0C012vsTklB) ? BIT(kNoV0C012vsTklBG) : 0; // o2-linter: disable=magic-number (nTkl dependent parameterization)

    // copy rct flags from bcsel table
    uint32_t rct = bc.rct_raw();

    // apply int7-like selections
    bool sel7 = 1;
    for (int i = 0; i < kNsel; i++) {
      sel7 = sel7 && (applySelection[i] ? TESTBIT(selection, i) : 1);
    }

    // TODO introduce array of sel[0]... sel[8] or similar?
    bool sel8 = bc.selection_bit(kIsBBT0A) && bc.selection_bit(kIsBBT0C); // TODO apply other cuts for sel8
    bool sel1 = bc.selection_bit(kIsINT1);
    sel1 = sel1 && bc.selection_bit(kNoBGV0A);
    sel1 = sel1 && bc.selection_bit(kNoBGV0C);
    sel1 = sel1 && bc.selection_bit(kNoTPCLaserWarmUp);
    sel1 = sel1 && bc.selection_bit(kNoTPCHVdip);

    // INT1 (SPDFO>0 | V0A | V0C) minimum bias trigger logic used in pp2010 and pp2011
    bool isINT1period = bc.runNumber() <= 136377 || (bc.runNumber() >= 144871 && bc.runNumber() <= 159582); // o2-linter: disable=magic-number (magic run numbers)

    // fill counters
    if (isMC == 1 || (!isINT1period && bc.alias_bit(kINT7)) || (isINT1period && bc.alias_bit(kINT1))) {
      histos.get<TH1>(HIST("hColCounterAll"))->Fill(Form("%d", bc.runNumber()), 1);
      if ((!isINT1period && sel7) || (isINT1period && sel1)) {
        histos.get<TH1>(HIST("hColCounterAcc"))->Fill(Form("%d", bc.runNumber()), 1);
      }
    }

    evsel(alias, selection, rct, sel7, sel8, foundBC, foundFT0, foundFV0, foundFDD, foundZDC, 0, 0, 0);
  }
  PROCESS_SWITCH(EventSelectionTask, processRun2, "Process Run2 event selection", true);

  Partition<FullTracksIU> pvTracks = ((aod::track::flags & static_cast<uint32_t>(o2::aod::track::PVContributor)) == static_cast<uint32_t>(o2::aod::track::PVContributor));
  void processRun3(aod::Collisions const& cols, FullTracksIU const&, BCsWithBcSelsRun3 const& bcs, aod::FT0s const&)
  {
    int run = bcs.iteratorAt(0).runNumber();
    // extract bc pattern from CCDB for data or anchored MC only
    int run3min = 500000;
    if (run != lastRun && run >= run3min) {
      lastRun = run;
      auto runInfo = o2::parameters::AggregatedRunInfo::buildAggregatedRunInfo(o2::ccdb::BasicCCDBManager::instance(), run, strLPMProductionTag);
      // first bc of the first orbit
      bcSOR = runInfo.orbitSOR * nBCsPerOrbit;
      // duration of TF in bcs
      nBCsPerTF = confNumberOfOrbitsPerTF < 0 ? runInfo.orbitsPerTF * nBCsPerOrbit : confNumberOfOrbitsPerTF * nBCsPerOrbit;
      // colliding bc pattern
      int64_t ts = bcs.iteratorAt(0).timestamp();
      auto grplhcif = ccdb->getForTimeStamp<o2::parameters::GRPLHCIFData>("GLO/Config/GRPLHCIF", ts);
      bcPatternB = grplhcif->getBunchFilling().getBCPattern();

      // extract ITS ROF parameters
      auto alppar = ccdb->getForTimeStamp<o2::itsmft::DPLAlpideParam<0>>("ITS/Config/AlpideParam", ts);
      rofOffset = alppar->roFrameBiasInBC;
      rofLength = alppar->roFrameLengthInBC;
      LOGP(debug, "ITS ROF Offset={} ITS ROF Length={}", rofOffset, rofLength);
    } // if run != lastRun

    // create maps from globalBC to bc index for TVX-fired bcs
    // to be used for closest TVX searches
    std::map<int64_t, int32_t> mapGlobalBcWithTVX;
    std::map<int64_t, float> mapGlobalBcVtxZ;
    for (const auto& bc : bcs) {
      int64_t globalBC = bc.globalBC();
      // skip non-colliding bcs for data and anchored runs
      if (run >= run3min && bcPatternB[globalBC % nBCsPerOrbit] == 0) {
        continue;
      }
      if (bc.selection_bit(kIsTriggerTVX)) {
        mapGlobalBcWithTVX[globalBC] = bc.globalIndex();
        mapGlobalBcVtxZ[globalBC] = bc.has_ft0() ? bc.ft0().posZ() : 0;
      }
    }

    // protection against empty FT0 maps
    if (mapGlobalBcWithTVX.size() == 0) {
      LOGP(error, "FT0 table is empty or corrupted. Filling evsel table with dummy values");
      for (const auto& col : cols) {
        auto bc = col.bc_as<BCsWithBcSelsRun3>();
        int32_t foundBC = bc.globalIndex();
        int32_t foundFT0 = bc.foundFT0Id();
        int32_t foundFV0 = bc.foundFV0Id();
        int32_t foundFDD = bc.foundFDDId();
        int32_t foundZDC = bc.foundZDCId();
        uint32_t rct = 0;
        evsel(bc.alias_raw(), bc.selection_raw(), rct, kFALSE, kFALSE, foundBC, foundFT0, foundFV0, foundFDD, foundZDC, -1, -1, -1);
      }
      return;
    }
    std::vector<int> vTracksITS567perColl(cols.size(), 0);                                    // counter of tracks per collision for occupancy studies
    std::vector<float> vAmpFT0CperColl(cols.size(), 0);                                       // amplitude FT0C per collision
    std::vector<float> vCollVz(cols.size(), 0);                                               // vector with vZ positions for each collision
    std::vector<bool> vIsFullInfoForOccupancy(cols.size(), 0);                                // info for occupancy in +/- windows is available (i.e. a given coll is not too close to the TF borders)
    const float timeWinOccupancyCalcMinNS = confTimeIntervalForOccupancyCalculationMin * 1e3; // ns
    const float timeWinOccupancyCalcMaxNS = confTimeIntervalForOccupancyCalculationMax * 1e3; // ns
    std::vector<bool> vIsVertexITSTPC(cols.size(), 0);                                        // at least one of vertex contributors is ITS-TPC track
    std::vector<bool> vIsVertexTOFmatched(cols.size(), 0);                                    // at least one of vertex contributors is matched to TOF
    std::vector<bool> vIsVertexTRDmatched(cols.size(), 0);                                    // at least one of vertex contributors is matched to TRD

    std::vector<int> vCollisionsPerBc(bcs.size(), 0);    // counter of collisions per found bc for pileup checks
    std::vector<int> vFoundBCindex(cols.size(), -1);     // indices of found bcs
    std::vector<int64_t> vFoundGlobalBC(cols.size(), 0); // global BCs for collisions

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
      auto bc = col.bc_as<BCsWithBcSelsRun3>();

      vCollVz[colIndex] = col.posZ();

      int64_t globalBC = bc.globalBC();
      int bcInTF = (bc.globalBC() - bcSOR) % nBCsPerTF;
      vIsFullInfoForOccupancy[colIndex] = ((bcInTF - 300) * bcNS > -timeWinOccupancyCalcMinNS) && ((nBCsPerTF - 4000 - bcInTF) * bcNS > timeWinOccupancyCalcMaxNS) ? true : false;

      const auto& colPvTracks = pvTracks.sliceByCached(aod::track::collisionId, col.globalIndex(), cache);
      vTrackTimesTOF.clear();
      vTrackTimesTRDnoTOF.clear();
      int nPvTracksTPCnoTOFnoTRD = 0;
      int nPvTracksHighPtTPCnoTOFnoTRD = 0;
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

      int64_t foundGlobalBC = 0;
      int32_t foundBCindex = -1;

      if (nPvTracksTOF > 0) {
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
        int64_t bestGlobalBC = findBestGlobalBC(meanBC, confSigmaBCforHighPtTracks, vNcontributors[colIndex], col.posZ(), mapGlobalBcVtxZ);
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

    // second loop to match remaining low-pt TPCnoTOFnoTRD collisions
    for (const auto& col : cols) {
      int32_t colIndex = col.globalIndex();
      if (vIsVertexTPC[colIndex] > 0 && vIsVertexTOF[colIndex] == 0 && vIsVertexHighPtTPC[colIndex] == 0) {
        float weightedTime = vWeightedTimesTPCnoTOFnoTRD[colIndex];
        float weightedSigma = vWeightedSigmaTPCnoTOFnoTRD[colIndex];
        auto bc = col.bc_as<BCsWithBcSelsRun3>();
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

    // save indices of collisions for occupancy calculation (both in ROF and in time range)
    std::vector<std::vector<int>> vCollsInSameITSROF;
    std::vector<std::vector<int>> vCollsInPrevITSROF;
    std::vector<std::vector<int>> vCollsInTimeWin;
    std::vector<std::vector<float>> vTimeDeltaForColls; // delta time wrt a given collision
    for (const auto& col : cols) {
      int32_t colIndex = col.globalIndex();
      int64_t foundGlobalBC = vFoundGlobalBC[colIndex];
      auto bc = bcs.iteratorAt(vFoundBCindex[colIndex]);
      if (bc.has_foundFT0())
        vAmpFT0CperColl[colIndex] = bc.foundFT0().sumAmpC();

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
      // protection against the TF borders
      if (!vIsFullInfoForOccupancy[colIndex]) {
        vCollsInTimeWin.push_back(vAssocToThisCol);
        vTimeDeltaForColls.push_back(vCollsTimeDeltaWrtGivenColl);
        continue;
      }
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
        if (vAmpFT0CperColl[thisColIndex] > confFT0CamplCutVetoOnCollInROF)
          nCollsInRofWithFT0CAboveVetoStandard++;
        if (std::fabs(vCollVz[thisColIndex] - vZ) < confEpsilonVzDiffVetoInROF)
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
      vNoHighMultCollInPrevRof[colIndex] = (totalFT0amplInPrevROF < confFT0CamplCutVetoOnCollInROF);

      // ### occupancy in time windows
      // protection against TF borders
      if (!vIsFullInfoForOccupancy[colIndex]) { // occupancy in undefined (too close to TF borders)
        vNumTracksITS567inFullTimeWin[colIndex] = -1;
        vSumAmpFT0CinFullTimeWin[colIndex] = -1;
        continue;
      }
      std::vector<int> vAssocToThisCol = vCollsInTimeWin[colIndex];
      std::vector<float> vCollsTimeDeltaWrtGivenColl = vTimeDeltaForColls[colIndex];
      int nITS567tracksInFullTimeWindow = 0;
      float sumAmpFT0CInFullTimeWindow = 0;
      int nITS567tracksForVetoNarrow = 0;      // to veto events with nearby collisions (narrower range)
      int nITS567tracksForVetoStrict = 0;      // to veto events with nearby collisions
      int nCollsWithFT0CAboveVetoStandard = 0; // to veto events with per-collision multiplicity above threshold
      for (uint32_t iCol = 0; iCol < vAssocToThisCol.size(); iCol++) {
        int thisColIndex = vAssocToThisCol[iCol];
        float dt = vCollsTimeDeltaWrtGivenColl[iCol] / 1e3; // ns -> us
        float wOccup = 1.;
        if (confUseWeightsForOccupancyVariable) {
          // weighted occupancy
          wOccup = 0;
          if (dt >= -40 && dt < -5)                     // collisions in the past                    // o2-linter: disable=magic-number (to be checked by Igor)
            wOccup = 1. / 1225 * (dt + 40) * (dt + 40); // o2-linter: disable=magic-number (to be checked by Igor)
          else if (dt >= -5 && dt < 15)                 // collisions near a given one           // o2-linter: disable=magic-number (to be checked by Igor)
            wOccup = 1;
          // else if (dt >= 15 && dt < 100) // collisions from the future
          //   wOccup = -1. / 85 * dt + 20. / 17;
          else if (dt >= 15 && dt < 40)              // collisions from the future            // o2-linter: disable=magic-number (to be checked by Igor)
            wOccup = -0.4 / 25 * dt + 1.24;          // o2-linter: disable=magic-number (to be checked by Igor)
          else if (dt >= 40 && dt < 100)             // collisions from the distant future   // o2-linter: disable=magic-number (to be checked by Igor)
            wOccup = -0.4 / 60 * dt + 0.6 + 0.8 / 3; // o2-linter: disable=magic-number (to be checked by Igor)
        }
        nITS567tracksInFullTimeWindow += wOccup * vTracksITS567perColl[thisColIndex];
        sumAmpFT0CInFullTimeWindow += wOccup * vAmpFT0CperColl[thisColIndex];

        // counting tracks from other collisions in fixed time windows
        if (std::fabs(dt) < confTimeRangeVetoOnCollNarrow)
          nITS567tracksForVetoNarrow += vTracksITS567perColl[thisColIndex];
        if (std::fabs(dt) < confTimeRangeVetoOnCollStandard)
          nITS567tracksForVetoStrict += vTracksITS567perColl[thisColIndex];

        // standard cut on other collisions vs delta-times
        const float driftV = 2.5;  // drift velocity in cm/us, TPC drift_length / drift_time = 250 cm / 100 us
        if (std::fabs(dt) < 2.0) { // us, complete veto on other collisions  // o2-linter: disable=magic-number (to be checked by Igor)
          nCollsWithFT0CAboveVetoStandard++;
        } else if (dt > -4.0 && dt <= -2.0) { // us, strict veto to suppress fake ITS-TPC matches more  // o2-linter: disable=magic-number (to be checked by Igor)
          if (vAmpFT0CperColl[thisColIndex] > confFT0CamplCutVetoOnCollInTimeRange / 5)
            nCollsWithFT0CAboveVetoStandard++;
        } else if (std::fabs(dt) < 8 + std::fabs(vZ) / driftV) { // loose veto, 8 us corresponds to maximum possible |vZ|, which is ~20 cm  // o2-linter: disable=magic-number (to be checked by Igor)
          // counting number of other collisions with multiplicity above threshold
          if (vAmpFT0CperColl[thisColIndex] > confFT0CamplCutVetoOnCollInTimeRange)
            nCollsWithFT0CAboveVetoStandard++;
        }
      }
      vNumTracksITS567inFullTimeWin[colIndex] = nITS567tracksInFullTimeWindow; // occupancy by a sum of number of ITS tracks (without a current collision)
      vSumAmpFT0CinFullTimeWin[colIndex] = sumAmpFT0CInFullTimeWindow;         // occupancy by a sum of FT0C amplitudes (without a current collision)
      // occupancy flags based on nearby collisions
      vNoCollInTimeRangeNarrow[colIndex] = (nITS567tracksForVetoNarrow == 0);
      vNoCollInTimeRangeStrict[colIndex] = (nITS567tracksForVetoStrict == 0);
      vNoHighMultCollInTimeRange[colIndex] = (nCollsWithFT0CAboveVetoStandard == 0);
    }

    for (const auto& col : cols) {
      int32_t colIndex = col.globalIndex();
      int32_t foundBC = vFoundBCindex[colIndex];
      auto bc = bcs.iteratorAt(foundBC);
      int32_t foundFT0 = bc.foundFT0Id();
      int32_t foundFV0 = bc.foundFV0Id();
      int32_t foundFDD = bc.foundFDDId();
      int32_t foundZDC = bc.foundZDCId();

      // compare zVtx from FT0 and from PV
      bool isGoodZvtxFT0vsPV = bc.has_foundFT0() ? std::fabs(bc.foundFT0().posZ() - col.posZ()) < maxDiffZvtxFT0vsPV : 0;

      // copy alias decisions from bcsel table
      uint32_t alias = bc.alias_raw();

      // copy selection decisions from bcsel table
      uint64_t selection = bc.selection_raw();
      selection |= vCollisionsPerBc[foundBC] <= 1 ? BIT(kNoSameBunchPileup) : 0;
      selection |= vIsVertexITSTPC[colIndex] ? BIT(kIsVertexITSTPC) : 0;
      selection |= vIsVertexTOFmatched[colIndex] ? BIT(kIsVertexTOFmatched) : 0;
      selection |= vIsVertexTRDmatched[colIndex] ? BIT(kIsVertexTRDmatched) : 0;
      selection |= isGoodZvtxFT0vsPV ? BIT(kIsGoodZvtxFT0vsPV) : 0;

      // selection bits based on occupancy time pattern
      selection |= vNoCollInTimeRangeNarrow[colIndex] ? BIT(kNoCollInTimeRangeNarrow) : 0;
      selection |= vNoCollInTimeRangeStrict[colIndex] ? BIT(kNoCollInTimeRangeStrict) : 0;
      selection |= vNoHighMultCollInTimeRange[colIndex] ? BIT(kNoCollInTimeRangeStandard) : 0;

      // selection bits based on ITS in-ROF occupancy
      selection |= vNoCollInSameRofStrict[colIndex] ? BIT(kNoCollInRofStrict) : 0;
      selection |= (vNoCollInSameRofStandard[colIndex] && vNoCollInSameRofWithCloseVz[colIndex]) ? BIT(kNoCollInRofStandard) : 0;
      selection |= vNoHighMultCollInPrevRof[colIndex] ? BIT(kNoHighMultCollInPrevRof) : 0;

      // copy rct flags from bcsel table
      uint32_t rct = bc.rct_raw();

      // apply int7-like selections
      bool sel7 = 0;

      // TODO apply other cuts for sel8
      // TODO introduce sel1 etc?
      // TODO introduce array of sel[0]... sel[8] or similar?
      bool sel8 = bc.selection_bit(kIsTriggerTVX) && bc.selection_bit(kNoTimeFrameBorder) && bc.selection_bit(kNoITSROFrameBorder);

      // fill counters
      histos.get<TH1>(HIST("hColCounterAll"))->Fill(Form("%d", bc.runNumber()), 1);
      if (bc.selection_bit(kIsTriggerTVX)) {
        histos.get<TH1>(HIST("hColCounterTVX"))->Fill(Form("%d", bc.runNumber()), 1);
      }
      if (sel8) {
        histos.get<TH1>(HIST("hColCounterAcc"))->Fill(Form("%d", bc.runNumber()), 1);
      }

      evsel(alias, selection, rct, sel7, sel8, foundBC, foundFT0, foundFV0, foundFDD, foundZDC,
            vNumTracksITS567inFullTimeWin[colIndex], vSumAmpFT0CinFullTimeWin[colIndex], 0);
    }
  }

  PROCESS_SWITCH(EventSelectionTask, processRun3, "Process Run3 event selection", false);
};

struct LumiTask {
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

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

  void init(InitContext&)
  {
    histos.add("hCounterTVX", "", kTH1D, {{1, 0., 1.}});
    histos.add("hCounterTCE", "", kTH1D, {{1, 0., 1.}});
    histos.add("hCounterZEM", "", kTH1D, {{1, 0., 1.}});
    histos.add("hCounterZNC", "", kTH1D, {{1, 0., 1.}});
    histos.add("hCounterTVXZDC", "", kTH1D, {{1, 0., 1.}});
    histos.add("hCounterTVXafterBCcuts", "", kTH1D, {{1, 0., 1.}});
    histos.add("hCounterTCEafterBCcuts", "", kTH1D, {{1, 0., 1.}});
    histos.add("hCounterZEMafterBCcuts", "", kTH1D, {{1, 0., 1.}});
    histos.add("hCounterZNCafterBCcuts", "", kTH1D, {{1, 0., 1.}});
    histos.add("hCounterTVXZDCafterBCcuts", "", kTH1D, {{1, 0., 1.}});
    histos.add("hLumiTVX", ";;Luminosity, 1/#mub", kTH1D, {{1, 0., 1.}});
    histos.add("hLumiTCE", ";;Luminosity, 1/#mub", kTH1D, {{1, 0., 1.}});
    histos.add("hLumiZEM", ";;Luminosity, 1/#mub", kTH1D, {{1, 0., 1.}});
    histos.add("hLumiZNC", ";;Luminosity, 1/#mub", kTH1D, {{1, 0., 1.}});
    histos.add("hLumiTVXafterBCcuts", ";;Luminosity, 1/#mub", kTH1D, {{1, 0., 1.}});
    histos.add("hLumiTCEafterBCcuts", ";;Luminosity, 1/#mub", kTH1D, {{1, 0., 1.}});
    histos.add("hLumiZEMafterBCcuts", ";;Luminosity, 1/#mub", kTH1D, {{1, 0., 1.}});
    histos.add("hLumiZNCafterBCcuts", ";;Luminosity, 1/#mub", kTH1D, {{1, 0., 1.}});

    const int nLists = 6;
    TString rctListNames[] = {"CBT", "CBT_hadronPID", "CBT_electronPID", "CBT_calo", "CBT_muon", "CBT_muon_glo"};
    histos.add("hLumiTVXafterBCcutsRCT", ";;Luminosity, 1/#mub", kTH2D, {{1, 0., 1.}, {4 * nLists, -0.5, 4. * nLists - 0.5}});
    histos.add("hLumiTCEafterBCcutsRCT", ";;Luminosity, 1/#mub", kTH2D, {{1, 0., 1.}, {4 * nLists, -0.5, 4. * nLists - 0.5}});
    histos.add("hLumiZEMafterBCcutsRCT", ";;Luminosity, 1/#mub", kTH2D, {{1, 0., 1.}, {4 * nLists, -0.5, 4. * nLists - 0.5}});
    histos.add("hLumiZNCafterBCcutsRCT", ";;Luminosity, 1/#mub", kTH2D, {{1, 0., 1.}, {4 * nLists, -0.5, 4. * nLists - 0.5}});

    for (int i = 0; i < nLists; i++) {
      const auto& rctListName = rctListNames[i];
      mRCTFlagsCheckers.emplace_back(rctListName.Data(), false, false); // disable zdc check, disable lim. acc. check
      mRCTFlagsCheckers.emplace_back(rctListName.Data(), false, true);  // disable zdc check, enable lim. acc. check
      mRCTFlagsCheckers.emplace_back(rctListName.Data(), true, false);  // enable zdc check, disable lim. acc. check
      mRCTFlagsCheckers.emplace_back(rctListName.Data(), true, true);   // enable zdc check, enable lim. acc. check
      histos.get<TH2>(HIST("hLumiTVXafterBCcutsRCT"))->GetYaxis()->SetBinLabel(4 * i + 1, rctListName.Data());
      histos.get<TH2>(HIST("hLumiTCEafterBCcutsRCT"))->GetYaxis()->SetBinLabel(4 * i + 1, rctListName.Data());
      histos.get<TH2>(HIST("hLumiZEMafterBCcutsRCT"))->GetYaxis()->SetBinLabel(4 * i + 1, rctListName.Data());
      histos.get<TH2>(HIST("hLumiZNCafterBCcutsRCT"))->GetYaxis()->SetBinLabel(4 * i + 1, rctListName.Data());
      histos.get<TH2>(HIST("hLumiTVXafterBCcutsRCT"))->GetYaxis()->SetBinLabel(4 * i + 2, (rctListName + "_fullacc").Data());
      histos.get<TH2>(HIST("hLumiTCEafterBCcutsRCT"))->GetYaxis()->SetBinLabel(4 * i + 2, (rctListName + "_fullacc").Data());
      histos.get<TH2>(HIST("hLumiZEMafterBCcutsRCT"))->GetYaxis()->SetBinLabel(4 * i + 2, (rctListName + "_fullacc").Data());
      histos.get<TH2>(HIST("hLumiZNCafterBCcutsRCT"))->GetYaxis()->SetBinLabel(4 * i + 2, (rctListName + "_fullacc").Data());
      histos.get<TH2>(HIST("hLumiTVXafterBCcutsRCT"))->GetYaxis()->SetBinLabel(4 * i + 3, (rctListName + "_zdc").Data());
      histos.get<TH2>(HIST("hLumiTCEafterBCcutsRCT"))->GetYaxis()->SetBinLabel(4 * i + 3, (rctListName + "_zdc").Data());
      histos.get<TH2>(HIST("hLumiZEMafterBCcutsRCT"))->GetYaxis()->SetBinLabel(4 * i + 3, (rctListName + "_zdc").Data());
      histos.get<TH2>(HIST("hLumiZNCafterBCcutsRCT"))->GetYaxis()->SetBinLabel(4 * i + 3, (rctListName + "_zdc").Data());
      histos.get<TH2>(HIST("hLumiTVXafterBCcutsRCT"))->GetYaxis()->SetBinLabel(4 * i + 4, (rctListName + "_zdc" + "_fullacc").Data());
      histos.get<TH2>(HIST("hLumiTCEafterBCcutsRCT"))->GetYaxis()->SetBinLabel(4 * i + 4, (rctListName + "_zdc" + "_fullacc").Data());
      histos.get<TH2>(HIST("hLumiZEMafterBCcutsRCT"))->GetYaxis()->SetBinLabel(4 * i + 4, (rctListName + "_zdc" + "_fullacc").Data());
      histos.get<TH2>(HIST("hLumiZNCafterBCcutsRCT"))->GetYaxis()->SetBinLabel(4 * i + 4, (rctListName + "_zdc" + "_fullacc").Data());
    }
  }

  void processRun2(aod::BCs const&)
  {
    LOGP(debug, "Dummy process function for Run 2");
  }

  PROCESS_SWITCH(LumiTask, processRun2, "Process Run2 lumi task", true);

  void processRun3(BCsWithBcSelsRun3 const& bcs, aod::FT0s const&)
  {
    if (bcs.size() == 0)
      return;
    int run = bcs.iteratorAt(0).runNumber();
    if (run < 500000) // o2-linter: disable=magic-number (skip for unanchored MCs)
      return;
    if (run != lastRun && run >= 520259) { // o2-linter: disable=magic-number (scalers available for runs above 520120)
      lastRun = run;
      int64_t ts = bcs.iteratorAt(0).timestamp();

      // getting GRP LHCIF object to extract colliding system, energy and colliding bc pattern
      auto grplhcif = ccdb->getForTimeStamp<parameters::GRPLHCIFData>("GLO/Config/GRPLHCIF", ts);
      int beamZ1 = grplhcif->getBeamZ(constants::lhc::BeamA);
      int beamZ2 = grplhcif->getBeamZ(constants::lhc::BeamC);
      float sqrts = grplhcif->getSqrtS();
      int nCollidingBCs = grplhcif->getBunchFilling().getNBunches();
      bcPatternB = grplhcif->getBunchFilling().getBCPattern();

      // visible cross sections in ub. Using dummy -1 if lumi estimator is not reliable for this colliding system
      csTVX = -1;
      csTCE = -1;
      csZEM = -1;
      csZNC = -1;
      // Temporary workaround to get visible cross section. TODO: store run-by-run visible cross sections in CCDB
      if (beamZ1 == 1 && beamZ2 == 1) {
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
      auto config = ccdb->getSpecific<o2::ctp::CTPConfiguration>("CTP/Config/Config", ts, metadata);
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
      auto scalers = ccdb->getSpecific<ctp::CTPRunScalers>("CTP/Calib/Scalers", ts, metadata);
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

    const char* srun = Form("%d", run);

    for (const auto& bc : bcs) {
      auto& selection = bc.selection_raw();
      if (bcPatternB[bc.globalBC() % nBCsPerOrbit] == 0) // skip non-colliding bcs
        continue;

      bool noBorder = TESTBIT(selection, kNoTimeFrameBorder) && TESTBIT(selection, kNoITSROFrameBorder);
      bool isTriggerTVX = TESTBIT(selection, kIsTriggerTVX);
      bool isTriggerTCE = bc.has_ft0() ? (TESTBIT(selection, kIsTriggerTVX) && TESTBIT(bc.ft0().triggerMask(), o2::ft0::Triggers::bitCen)) : 0;
      bool isTriggerZNA = TESTBIT(selection, kIsBBZNA);
      bool isTriggerZNC = TESTBIT(selection, kIsBBZNC);
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

      if (isTriggerTVX) {
        histos.get<TH1>(HIST("hCounterTVX"))->Fill(srun, 1);
        histos.get<TH1>(HIST("hLumiTVX"))->Fill(srun, lumiTVX);
        if (isTriggerZNA && isTriggerZNC) {
          histos.get<TH1>(HIST("hCounterTVXZDC"))->Fill(srun, 1);
        }
        if (noBorder) {
          histos.get<TH1>(HIST("hCounterTVXafterBCcuts"))->Fill(srun, 1);
          histos.get<TH1>(HIST("hLumiTVXafterBCcuts"))->Fill(srun, lumiTVX);
          if (isTriggerZNA && isTriggerZNC) {
            histos.get<TH1>(HIST("hCounterTVXZDCafterBCcuts"))->Fill(srun, 1);
          }
          for (size_t i = 0; i < mRCTFlagsCheckers.size(); i++) {
            if (mRCTFlagsCheckers[i](bc))
              histos.get<TH2>(HIST("hLumiTVXafterBCcutsRCT"))->Fill(srun, i, lumiTVX);
          }
        }
      }

      if (isTriggerTCE) {
        histos.get<TH1>(HIST("hCounterTCE"))->Fill(srun, 1);
        histos.get<TH1>(HIST("hLumiTCE"))->Fill(srun, lumiTCE);
        if (noBorder) {
          histos.get<TH1>(HIST("hCounterTCEafterBCcuts"))->Fill(srun, 1);
          histos.get<TH1>(HIST("hLumiTCEafterBCcuts"))->Fill(srun, lumiTCE);
          for (size_t i = 0; i < mRCTFlagsCheckers.size(); i++) {
            if (mRCTFlagsCheckers[i](bc))
              histos.get<TH2>(HIST("hLumiTCEafterBCcutsRCT"))->Fill(srun, i, lumiTCE);
          }
        }
      }

      if (isTriggerZEM) {
        histos.get<TH1>(HIST("hCounterZEM"))->Fill(srun, 1);
        histos.get<TH1>(HIST("hLumiZEM"))->Fill(srun, lumiZEM);
        if (noBorder) {
          histos.get<TH1>(HIST("hCounterZEMafterBCcuts"))->Fill(srun, 1);
          histos.get<TH1>(HIST("hLumiZEMafterBCcuts"))->Fill(srun, lumiZEM);
          for (size_t i = 0; i < mRCTFlagsCheckers.size(); i++) {
            if (mRCTFlagsCheckers[i](bc))
              histos.get<TH2>(HIST("hLumiZEMafterBCcutsRCT"))->Fill(srun, i, lumiZEM);
          }
        }
      }

      if (isTriggerZNC) {
        histos.get<TH1>(HIST("hCounterZNC"))->Fill(srun, 1);
        histos.get<TH1>(HIST("hLumiZNC"))->Fill(srun, lumiZNC);
        if (noBorder) {
          histos.get<TH1>(HIST("hCounterZNCafterBCcuts"))->Fill(srun, 1);
          histos.get<TH1>(HIST("hLumiZNCafterBCcuts"))->Fill(srun, lumiZNC);
          for (size_t i = 0; i < mRCTFlagsCheckers.size(); i++) {
            if (mRCTFlagsCheckers[i](bc))
              histos.get<TH2>(HIST("hLumiZNCafterBCcutsRCT"))->Fill(srun, i, lumiZNC);
          }
        }
      }

    } // bcs
  } // process
  PROCESS_SWITCH(LumiTask, processRun3, "Process Run3 lumi task", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  // Parse the metadata
  metadataInfo.initMetadata(cfgc);

  return WorkflowSpec{
    adaptAnalysisTask<BcSelectionTask>(cfgc),
    adaptAnalysisTask<EventSelectionTask>(cfgc),
    adaptAnalysisTask<LumiTask>(cfgc)};
}
