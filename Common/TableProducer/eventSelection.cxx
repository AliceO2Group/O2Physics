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
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/TriggerAliases.h"
#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/LHCConstants.h"
#include "Framework/HistogramRegistry.h"
#include "DataFormatsFT0/Digit.h"
#include "DataFormatsParameters/GRPLHCIFData.h"
#include "DataFormatsParameters/GRPECSObject.h"
#include "ITSMFTBase/DPLAlpideParam.h"
#include "MetadataHelper.h"

#include "TH1D.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::aod::evsel;

MetadataHelper metadataInfo; // Metadata helper

using BCsWithRun2InfosTimestampsAndMatches = soa::Join<aod::BCs, aod::Run2BCInfos, aod::Timestamps, aod::Run2MatchedToBCSparse>;
using BCsWithRun3Matchings = soa::Join<aod::BCs, aod::Timestamps, aod::Run3MatchedToBCSparse>;
using BCsWithBcSelsRun2 = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels, aod::Run2BCInfos, aod::Run2MatchedToBCSparse>;
using BCsWithBcSelsRun3 = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels>;
using FullTracksIU = soa::Join<aod::TracksIU, aod::TracksExtra>;

struct BcSelectionTask {
  Produces<aod::BcSels> bcsel;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Configurable<int> confTriggerBcShift{"triggerBcShift", 999, "set to 294 for apass2/apass3 in LHC22o-t"};
  Configurable<int> confITSROFrameStartBorderMargin{"ITSROFrameStartBorderMargin", -1, "Number of bcs at the start of ITS RO Frame border. Take from CCDB if -1"};
  Configurable<int> confITSROFrameEndBorderMargin{"ITSROFrameEndBorderMargin", -1, "Number of bcs at the end of ITS RO Frame border. Take from CCDB if -1"};
  Configurable<int> confTimeFrameStartBorderMargin{"TimeFrameStartBorderMargin", -1, "Number of bcs to cut at the start of the Time Frame. Take from CCDB if -1"};
  Configurable<int> confTimeFrameEndBorderMargin{"TimeFrameEndBorderMargin", -1, "Number of bcs to cut at the end of the Time Frame. Take from CCDB if -1"};
  Configurable<bool> confCheckRunDurationLimits{"checkRunDurationLimits", false, "Check if the BCs are within the run duration limits"};

  int lastRunNumber = -1;
  int64_t bcSOR = -1;                    // global bc of the start of the first orbit
  int64_t nBCsPerTF = -1;                // duration of TF in bcs, should be 128*3564 or 32*3564
  int mITSROFrameStartBorderMargin = 10; // default value
  int mITSROFrameEndBorderMargin = 20;   // default value
  int mTimeFrameStartBorderMargin = 300; // default value
  int mTimeFrameEndBorderMargin = 4000;  // default value

  void init(InitContext&)
  {
    if (metadataInfo.isFullyDefined() && !doprocessRun2 && !doprocessRun3) { // Check if the metadata is initialized (only if not forced from the workflow configuration)
      LOG(info) << "Autosetting the processing mode (Run2 or Run3) based on metadata";
      if (metadataInfo.isRun3()) {
        doprocessRun3.value = true;
      } else {
        doprocessRun2.value = false;
      }
    }

    // ccdb->setURL("http://ccdb-test.cern.ch:8080");
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    histos.add("hCounterTVX", "", kTH1D, {{1, 0., 1.}});
    histos.add("hCounterTCE", "", kTH1D, {{1, 0., 1.}});
    histos.add("hCounterZEM", "", kTH1D, {{1, 0., 1.}});
    histos.add("hCounterZNC", "", kTH1D, {{1, 0., 1.}});
    histos.add("hCounterTVXafterBCcuts", "", kTH1D, {{1, 0., 1.}});
    histos.add("hCounterTCEafterBCcuts", "", kTH1D, {{1, 0., 1.}});
    histos.add("hCounterZEMafterBCcuts", "", kTH1D, {{1, 0., 1.}});
    histos.add("hCounterZNCafterBCcuts", "", kTH1D, {{1, 0., 1.}});
    histos.add("hCounterInvalidBCTimestamp", "", kTH1D, {{1, 0., 1.}});
    histos.add("hLumiTVX", ";;Luminosity, 1/#mub", kTH1D, {{1, 0., 1.}});
    histos.add("hLumiTCE", ";;Luminosity, 1/#mub", kTH1D, {{1, 0., 1.}});
    histos.add("hLumiZEM", ";;Luminosity, 1/#mub", kTH1D, {{1, 0., 1.}});
    histos.add("hLumiZNC", ";;Luminosity, 1/#mub", kTH1D, {{1, 0., 1.}});
    histos.add("hLumiTVXafterBCcuts", ";;Luminosity, 1/#mub", kTH1D, {{1, 0., 1.}});
    histos.add("hLumiTCEafterBCcuts", ";;Luminosity, 1/#mub", kTH1D, {{1, 0., 1.}});
    histos.add("hLumiZEMafterBCcuts", ";;Luminosity, 1/#mub", kTH1D, {{1, 0., 1.}});
    histos.add("hLumiZNCafterBCcuts", ";;Luminosity, 1/#mub", kTH1D, {{1, 0., 1.}});
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

    for (auto& bc : bcs) {
      EventSelectionParams* par = ccdb->getForTimeStamp<EventSelectionParams>("EventSelection/EventSelectionParams", bc.timestamp());
      TriggerAliases* aliases = ccdb->getForTimeStamp<TriggerAliases>("EventSelection/TriggerAliases", bc.timestamp());
      // fill fired aliases
      uint32_t alias{0};
      uint64_t triggerMask = bc.triggerMask();
      for (auto& al : aliases->GetAliasToTriggerMaskMap()) {
        if (triggerMask & al.second) {
          alias |= BIT(al.first);
        }
      }
      uint64_t triggerMaskNext50 = bc.triggerMaskNext50();
      for (auto& al : aliases->GetAliasToTriggerMaskNext50Map()) {
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
      selection |= !(fabs(timeZNA) > par->fZNABGlower && fabs(timeZNA) < par->fZNABGupper) ? BIT(kNoBGZNA) : 0;
      selection |= !(fabs(timeZNC) > par->fZNCBGlower && fabs(timeZNC) < par->fZNCBGupper) ? BIT(kNoBGZNC) : 0;
      selection |= (pow((timeZNA + timeZNC - par->fZNSumMean) / par->fZNSumSigma, 2) + pow((timeZNA - timeZNC - par->fZNDifMean) / par->fZNDifSigma, 2) < 1) ? BIT(kIsBBZAC) : 0;

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

      // Fill bc selection columns
      bcsel(alias, selection, foundFT0, foundFV0, foundFDD, foundZDC);
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
    // extract ITS time frame parameters
    int run = bcs.iteratorAt(0).runNumber();
    auto timestamps = ccdb->getRunDuration(run, true); /// fatalise if timestamps are not found
    int64_t sorTimestamp = timestamps.first;           // timestamp of the SOR/SOX/STF in ms
    int64_t eorTimestamp = timestamps.second;          // timestamp of the EOR/EOX/ETF in ms
    int64_t ts = eorTimestamp / 2 + sorTimestamp / 2;  // timestamp of the middle of the run
    auto alppar = ccdb->getForTimeStamp<o2::itsmft::DPLAlpideParam<0>>("ITS/Config/AlpideParam", ts);
    EventSelectionParams* par = ccdb->getForTimeStamp<EventSelectionParams>("EventSelection/EventSelectionParams", ts);
    TriggerAliases* aliases = ccdb->getForTimeStamp<TriggerAliases>("EventSelection/TriggerAliases", ts);
    // map from GlobalBC to BcId needed to find triggerBc
    std::map<uint64_t, int32_t> mapGlobalBCtoBcId;
    for (auto& bc : bcs) {
      mapGlobalBCtoBcId[bc.globalBC()] = bc.globalIndex();
    }
    int triggerBcShift = confTriggerBcShift;
    if (confTriggerBcShift == 999) {
      triggerBcShift = (run <= 526766 || (run >= 526886 && run <= 527237) || (run >= 527259 && run <= 527518) || run == 527523 || run == 527734 || run >= 534091) ? 0 : 294;
    }

    // extract run number and related information
    if (run != lastRunNumber) {
      lastRunNumber = run; // do it only once
      if (run >= 500000) { // access CCDB for data or anchored MC only
        int64_t ts = bcs.iteratorAt(0).timestamp();
        // access orbitShift, ITSROF and TF border margins
        EventSelectionParams* par = ccdb->getForTimeStamp<EventSelectionParams>("EventSelection/EventSelectionParams", ts);
        mITSROFrameStartBorderMargin = confITSROFrameStartBorderMargin < 0 ? par->fITSROFrameStartBorderMargin : confITSROFrameStartBorderMargin;
        mITSROFrameEndBorderMargin = confITSROFrameEndBorderMargin < 0 ? par->fITSROFrameEndBorderMargin : confITSROFrameEndBorderMargin;
        mTimeFrameStartBorderMargin = confTimeFrameStartBorderMargin < 0 ? par->fTimeFrameStartBorderMargin : confTimeFrameStartBorderMargin;
        mTimeFrameEndBorderMargin = confTimeFrameEndBorderMargin < 0 ? par->fTimeFrameEndBorderMargin : confTimeFrameEndBorderMargin;
        // access orbit-reset timestamp
        auto ctpx = ccdb->getForTimeStamp<std::vector<Long64_t>>("CTP/Calib/OrbitReset", ts);
        int64_t tsOrbitReset = (*ctpx)[0]; // us
        // access TF duration, start-of-run and end-of-run timestamps from ECS GRP
        std::map<std::string, std::string> metadata;
        metadata["runNumber"] = Form("%d", run);
        auto grpecs = ccdb->getSpecific<o2::parameters::GRPECSObject>("GLO/Config/GRPECS", ts, metadata);
        uint32_t nOrbitsPerTF = grpecs->getNHBFPerTF(); // assuming 1 orbit = 1 HBF;  nOrbitsPerTF=128 in 2022, 32 in 2023
        int64_t tsSOR = grpecs->getTimeStart();         // ms
        // calculate SOR orbit
        int64_t orbitSOR = (tsSOR * 1000 - tsOrbitReset) / o2::constants::lhc::LHCOrbitMUS;
        // adjust to the nearest TF edge
        orbitSOR = orbitSOR / nOrbitsPerTF * nOrbitsPerTF + par->fTimeFrameOrbitShift;
        // first bc of the first orbit (should coincide with TF start)
        bcSOR = orbitSOR * o2::constants::lhc::LHCMaxBunches;
        // duration of TF in bcs
        nBCsPerTF = nOrbitsPerTF * o2::constants::lhc::LHCMaxBunches;
        LOGP(info, "tsOrbitReset={} us, SOR = {} ms, orbitSOR = {}, nBCsPerTF = {}", tsOrbitReset, tsSOR, orbitSOR, nBCsPerTF);
      }
    }

    // bc loop
    for (auto bc : bcs) {
      uint32_t alias{0};
      // workaround for pp2022 (trigger info is shifted by -294 bcs)
      int32_t triggerBcId = mapGlobalBCtoBcId[bc.globalBC() + triggerBcShift];
      if (triggerBcId) {
        auto triggerBc = bcs.iteratorAt(triggerBcId);
        uint64_t triggerMask = triggerBc.triggerMask();
        for (auto& al : aliases->GetAliasToTriggerMaskMap()) {
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
        if (bc.globalBC() + 1 == globalBC) {
          timeV0ABG = bc.has_fv0a() ? bc.fv0a().time() : -999.f;
          timeT0ABG = bc.has_ft0() ? bc.ft0().timeA() : -999.f;
          timeT0CBG = bc.has_ft0() ? bc.ft0().timeC() : -999.f;
        }
        if (bc.globalBC() + 5 == globalBC) {
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
      selection |= (pow((timeZNA + timeZNC - par->fZNSumMean) / par->fZNSumSigma, 2) + pow((timeZNA - timeZNC - par->fZNDifMean) / par->fZNDifSigma, 2) < 1) ? BIT(kIsBBZAC) : 0;
      selection |= !(fabs(timeZNA) > par->fZNABGlower && fabs(timeZNA) < par->fZNABGupper) ? BIT(kNoBGZNA) : 0;
      selection |= !(fabs(timeZNC) > par->fZNCBGlower && fabs(timeZNC) < par->fZNCBGupper) ? BIT(kNoBGZNC) : 0;
      selection |= (bc.has_ft0() ? (bc.ft0().triggerMask() & BIT(o2::ft0::Triggers::bitVertex)) > 0 : 0) ? BIT(kIsTriggerTVX) : 0;

      // check if bc is far (at least confITSROFrameBorderMargin) from the end of ITS RO Frame border
      // 2bc margin is also introduced at ehe beginning of ITS RO Frame to account for the uncertainty of the roFrameBiasInBC
      uint16_t bcInITSROF = (globalBC + 3564 - alppar->roFrameBiasInBC) % alppar->roFrameLengthInBC;
      LOGP(debug, "bcInITSROF={}", bcInITSROF);
      selection |= bcInITSROF > mITSROFrameStartBorderMargin && bcInITSROF < alppar->roFrameLengthInBC - mITSROFrameEndBorderMargin ? BIT(kNoITSROFrameBorder) : 0;

      // check if bc is far from the Time Frame borders
      int64_t bcInTF = (globalBC - bcSOR) % nBCsPerTF;
      LOGP(debug, "bcInTF={}", bcInTF);
      selection |= bcInTF > mTimeFrameStartBorderMargin && bcInTF < nBCsPerTF - mTimeFrameEndBorderMargin ? BIT(kNoTimeFrameBorder) : 0;

      int32_t foundFT0 = bc.has_ft0() ? bc.ft0().globalIndex() : -1;
      int32_t foundFV0 = bc.has_fv0a() ? bc.fv0a().globalIndex() : -1;
      int32_t foundFDD = bc.has_fdd() ? bc.fdd().globalIndex() : -1;
      int32_t foundZDC = bc.has_zdc() ? bc.zdc().globalIndex() : -1;
      LOGP(debug, "foundFT0={}", foundFT0);

      // Temporary workaround to get visible cross section. TODO: store run-by-run visible cross sections in CCDB
      const char* srun = Form("%d", run);
      auto grplhcif = ccdb->getForTimeStamp<o2::parameters::GRPLHCIFData>("GLO/Config/GRPLHCIF", bc.timestamp());
      int beamZ1 = grplhcif->getBeamZ(o2::constants::lhc::BeamA);
      int beamZ2 = grplhcif->getBeamZ(o2::constants::lhc::BeamC);
      bool isPP = beamZ1 == 1 && beamZ2 == 1;
      bool injectionEnergy = (run >= 500000 && run <= 520099) || (run >= 534133 && run <= 534468);
      // Cross sections in ub. Using dummy -1 if lumi estimator is not reliable
      float csTVX = isPP ? (injectionEnergy ? 0.0355e6 : 0.0594e6) : -1.;
      float csTCE = isPP ? -1. : 10.36e6;
      float csZEM = isPP ? -1. : 415.2e6; // see AN: https://alice-notes.web.cern.ch/node/1515
      float csZNC = isPP ? -1. : 214.5e6; // see AN: https://alice-notes.web.cern.ch/node/1515
      if (run > 543437 && run < 543514) {
        csTCE = 8.3e6;
      }
      if (run >= 543514) {
        csTCE = 4.10e6; // see AN: https://alice-notes.web.cern.ch/node/1515
      }

      // Fill TVX (T0 vertex) counters
      if (TESTBIT(selection, kIsTriggerTVX)) {
        histos.get<TH1>(HIST("hCounterTVX"))->Fill(srun, 1);
        histos.get<TH1>(HIST("hLumiTVX"))->Fill(srun, 1. / csTVX);
        if (TESTBIT(selection, kNoITSROFrameBorder) && TESTBIT(selection, kNoTimeFrameBorder)) {
          histos.get<TH1>(HIST("hCounterTVXafterBCcuts"))->Fill(srun, 1);
          histos.get<TH1>(HIST("hLumiTVXafterBCcuts"))->Fill(srun, 1. / csTVX);
        }
      }
      // Fill counters and lumi histograms for Pb-Pb lumi monitoring
      // TODO: introduce pileup correction
      if (bc.has_ft0() ? (TESTBIT(selection, kIsTriggerTVX) && TESTBIT(bc.ft0().triggerMask(), o2::ft0::Triggers::bitCen)) : 0) {
        histos.get<TH1>(HIST("hCounterTCE"))->Fill(srun, 1);
        histos.get<TH1>(HIST("hLumiTCE"))->Fill(srun, 1. / csTCE);
        if (TESTBIT(selection, kNoITSROFrameBorder) && TESTBIT(selection, kNoTimeFrameBorder)) {
          histos.get<TH1>(HIST("hCounterTCEafterBCcuts"))->Fill(srun, 1);
          histos.get<TH1>(HIST("hLumiTCEafterBCcuts"))->Fill(srun, 1. / csTCE);
        }
      }
      if (TESTBIT(selection, kIsBBZNA) || TESTBIT(selection, kIsBBZNC)) {
        histos.get<TH1>(HIST("hCounterZEM"))->Fill(srun, 1);
        histos.get<TH1>(HIST("hLumiZEM"))->Fill(srun, 1. / csZEM);
        if (TESTBIT(selection, kNoITSROFrameBorder) && TESTBIT(selection, kNoTimeFrameBorder)) {
          histos.get<TH1>(HIST("hCounterZEMafterBCcuts"))->Fill(srun, 1);
          histos.get<TH1>(HIST("hLumiZEMafterBCcuts"))->Fill(srun, 1. / csZEM);
        }
      }
      if (TESTBIT(selection, kIsBBZNC)) {
        histos.get<TH1>(HIST("hCounterZNC"))->Fill(srun, 1);
        histos.get<TH1>(HIST("hLumiZNC"))->Fill(srun, 1. / csZNC);
        if (TESTBIT(selection, kNoITSROFrameBorder) && TESTBIT(selection, kNoTimeFrameBorder)) {
          histos.get<TH1>(HIST("hCounterZNCafterBCcuts"))->Fill(srun, 1);
          histos.get<TH1>(HIST("hLumiZNCafterBCcuts"))->Fill(srun, 1. / csZNC);
        }
      }

      if (bc.timestamp() < static_cast<uint64_t>(sorTimestamp) || bc.timestamp() > static_cast<uint64_t>(eorTimestamp)) {
        histos.get<TH1>(HIST("hCounterInvalidBCTimestamp"))->Fill(srun, 1);
        if (confCheckRunDurationLimits.value) {
          LOGF(warn, "Invalid BC timestamp: %d, run: %d, sor: %d, eor: %d", bc.timestamp(), run, sorTimestamp, eorTimestamp);
          alias = 0u;
          selection = 0u;
        }
      }

      // Fill bc selection columns
      bcsel(alias, selection, foundFT0, foundFV0, foundFDD, foundZDC);
    }
  }
  PROCESS_SWITCH(BcSelectionTask, processRun3, "Process Run3 event selection", false);
};

struct EventSelectionTask {
  SliceCache cache;
  Produces<aod::EvSels> evsel;
  Configurable<std::string> syst{"syst", "PbPb", "pp, pPb, Pbp, PbPb, XeXe"}; // TODO determine from AOD metadata or from CCDB
  Configurable<int> muonSelection{"muonSelection", 0, "0 - barrel, 1 - muon selection with pileup cuts, 2 - muon selection without pileup cuts"};
  Configurable<float> maxDiffZvtxFT0vsPV{"maxDiffZvtxFT0vsPV", 1., "maximum difference (in cm) between z-vertex from FT0 and PV"};
  Configurable<int> isMC{"isMC", 0, "-1 - autoset, 0 - data, 1 - MC"};
  // configurables for occupancy-based event selection
  Configurable<float> confTimeIntervalForOccupancyCalculationMin{"TimeIntervalForOccupancyCalculationMin", -40, "Min time diff window for TPC occupancy calculation, us"};
  Configurable<float> confTimeIntervalForOccupancyCalculationMax{"TimeIntervalForOccupancyCalculationMax", 100, "Max time diff window for TPC occupancy calculation, us"};
  Configurable<std::vector<float>> confTimeBinsForOccupancyCalculation{"TimeBinsForOccupancyCalculation", {-40, -20, 0, 25, 50, 75, 100}, "Time bins for occupancy calculation and corresponding cuts (us)"};
  Configurable<std::vector<float>> confReferenceOccupanciesInTimeBins{"ReferenceOccupanciesInTimeBins", {3000, 1400, 750, 1000, 1750, 4000}, "Occupancy cuts in time bins (n tracks)"};
  Configurable<float> confTimeRangeVetoOnCollStandard{"TimeRangeVetoOnCollStandard", 10, "Exclusion of a collision if there are other collisions nearby, +/- us"};
  Configurable<float> confTimeRangeVetoOnCollNarrow{"TimeRangeVetoOnCollNarrow", 4, "Exclusion of a collision if there are other collisions nearby, +/- us"};
  Configurable<bool> confUseWeightsForOccupancyVariable{"UseWeightsForOccupancyEstimator", 1, "Use or not the delta-time weights for the occupancy estimator"};

  Partition<aod::Tracks> tracklets = (aod::track::trackType == static_cast<uint8_t>(o2::aod::track::TrackTypeEnum::Run2Tracklet));

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  int lastRun = -1;                                          // last run number (needed to access ccdb only if run!=lastRun)
  std::bitset<o2::constants::lhc::LHCMaxBunches> bcPatternB; // bc pattern of colliding bunches

  int64_t bcSOR = -1;     // global bc of the start of the first orbit
  int64_t nBCsPerTF = -1; // duration of TF in bcs, should be 128*3564 or 32*3564

  int32_t findClosest(int64_t globalBC, std::map<int64_t, int32_t>& bcs)
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

  void init(InitContext&)
  {
    if (metadataInfo.isFullyDefined()) { // Check if the metadata is initialized (only if not forced from the workflow configuration)
      if (!doprocessRun2 && !doprocessRun3) {
        LOG(info) << "Autosetting the processing mode (Run2 or Run3) based on metadata";
        if (metadataInfo.isRun3()) {
          doprocessRun3.value = true;
        } else {
          doprocessRun2.value = false;
        }
      }
      if (isMC == -1) {
        LOG(info) << "Autosetting the MC mode based on metadata";
        if (metadataInfo.isMC()) {
          isMC.value = 1;
        } else {
          isMC.value = 0;
        }
      }
    }

    // ccdb->setURL("http://ccdb-test.cern.ch:8080");
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    histos.add("hColCounterAll", "", kTH1D, {{1, 0., 1.}});
    histos.add("hColCounterTVX", "", kTH1D, {{1, 0., 1.}});
    histos.add("hColCounterAcc", "", kTH1D, {{1, 0., 1.}});
    // histos.add("hOccupancy", "", kTH1D, {{200, 0., 10000}});
  }

  void process(aod::Collisions const& collisions)
  {
    evsel.reserve(collisions.size());
  }

  void processRun2(aod::Collision const& col, BCsWithBcSelsRun2 const&, aod::Tracks const&, aod::FV0Cs const&)
  {
    auto bc = col.bc_as<BCsWithBcSelsRun2>();
    EventSelectionParams* par = ccdb->getForTimeStamp<EventSelectionParams>("EventSelection/EventSelectionParams", bc.timestamp());
    bool* applySelection = par->GetSelection(muonSelection);
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
    selection |= !(nTkl < 6 && multV0C012 > par->fV0C012vsTklA + nTkl * par->fV0C012vsTklB) ? BIT(kNoV0C012vsTklBG) : 0;

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
    bool isINT1period = bc.runNumber() <= 136377 || (bc.runNumber() >= 144871 && bc.runNumber() <= 159582);

    // fill counters
    if (isMC == 1 || (!isINT1period && bc.alias_bit(kINT7)) || (isINT1period && bc.alias_bit(kINT1))) {
      histos.get<TH1>(HIST("hColCounterAll"))->Fill(Form("%d", bc.runNumber()), 1);
      if ((!isINT1period && sel7) || (isINT1period && sel1)) {
        histos.get<TH1>(HIST("hColCounterAcc"))->Fill(Form("%d", bc.runNumber()), 1);
      }
    }

    evsel(alias, selection, sel7, sel8, foundBC, foundFT0, foundFV0, foundFDD, foundZDC, 0);
  }
  PROCESS_SWITCH(EventSelectionTask, processRun2, "Process Run2 event selection", true);

  Preslice<FullTracksIU> perCollision = aod::track::collisionId;
  void processRun3(aod::Collisions const& cols, FullTracksIU const& tracks, BCsWithBcSelsRun3 const& bcs, aod::FT0s const&)
  {
    int run = bcs.iteratorAt(0).runNumber();
    // extract bc pattern from CCDB for data or anchored MC only
    if (run != lastRun && run >= 500000) {
      lastRun = run;
      int64_t ts = bcs.iteratorAt(0).timestamp();
      auto grplhcif = ccdb->getForTimeStamp<o2::parameters::GRPLHCIFData>("GLO/Config/GRPLHCIF", ts);
      bcPatternB = grplhcif->getBunchFilling().getBCPattern();

      //
      EventSelectionParams* par = ccdb->getForTimeStamp<EventSelectionParams>("EventSelection/EventSelectionParams", ts);
      // access orbit-reset timestamp
      auto ctpx = ccdb->getForTimeStamp<std::vector<Long64_t>>("CTP/Calib/OrbitReset", ts);
      int64_t tsOrbitReset = (*ctpx)[0]; // us
      // access TF duration, start-of-run timestamp from ECS GRP
      std::map<std::string, std::string> metadata;
      metadata["runNumber"] = Form("%d", run);
      auto grpecs = ccdb->getSpecific<o2::parameters::GRPECSObject>("GLO/Config/GRPECS", ts, metadata);
      uint32_t nOrbitsPerTF = grpecs->getNHBFPerTF(); // assuming 1 orbit = 1 HBF;  nOrbitsPerTF=128 in 2022, 32 in 2023
      int64_t tsSOR = grpecs->getTimeStart();         // ms
      // calculate SOR orbit
      int64_t orbitSOR = (tsSOR * 1000 - tsOrbitReset) / o2::constants::lhc::LHCOrbitMUS;
      // adjust to the nearest TF edge
      orbitSOR = orbitSOR / nOrbitsPerTF * nOrbitsPerTF + par->fTimeFrameOrbitShift;
      // first bc of the first orbit (should coincide with TF start)
      bcSOR = orbitSOR * o2::constants::lhc::LHCMaxBunches;
      // duration of TF in bcs
      nBCsPerTF = nOrbitsPerTF * o2::constants::lhc::LHCMaxBunches;
    }

    // create maps from globalBC to bc index for TVX or FT0-OR fired bcs
    // to be used for closest TVX (FT0-OR) searches
    std::map<int64_t, int32_t> mapGlobalBcWithTVX;
    std::map<int64_t, int32_t> mapGlobalBcWithTOR;
    for (auto& bc : bcs) {
      int64_t globalBC = bc.globalBC();
      // skip non-colliding bcs for data and anchored runs
      if (run >= 500000 && bcPatternB[globalBC % o2::constants::lhc::LHCMaxBunches] == 0) {
        continue;
      }
      if (bc.selection_bit(kIsBBT0A) || bc.selection_bit(kIsBBT0C)) {
        mapGlobalBcWithTOR[globalBC] = bc.globalIndex();
      }
      if (bc.selection_bit(kIsTriggerTVX)) {
        mapGlobalBcWithTVX[globalBC] = bc.globalIndex();
      }
    }

    // protection against empty FT0 maps
    if (mapGlobalBcWithTOR.size() == 0 || mapGlobalBcWithTVX.size() == 0) {
      LOGP(error, "FT0 table is empty or corrupted. Filling evsel table with dummy values");
      for (auto& col : cols) {
        auto bc = col.bc_as<BCsWithBcSelsRun3>();
        int32_t foundBC = bc.globalIndex();
        int32_t foundFT0 = bc.foundFT0Id();
        int32_t foundFV0 = bc.foundFV0Id();
        int32_t foundFDD = bc.foundFDDId();
        int32_t foundZDC = bc.foundZDCId();
        evsel(bc.alias_raw(), bc.selection_raw(), kFALSE, kFALSE, foundBC, foundFT0, foundFV0, foundFDD, foundZDC, -1);
      }
      return;
    }

    std::vector<int> vFoundBCindex(cols.size(), -1);       // indices of found bcs
    std::vector<bool> vIsVertexITSTPC(cols.size(), 0);     // at least one of vertex contributors is ITS-TPC track
    std::vector<bool> vIsVertexTOFmatched(cols.size(), 0); // at least one of vertex contributors is matched to TOF
    std::vector<bool> vIsVertexTRDmatched(cols.size(), 0); // at least one of vertex contributors is matched to TRD
    std::vector<int> vCollisionsPerBc(bcs.size(), 0);      // counter of collisions per found bc for pileup checks

    // for the occupancy study
    std::vector<int64_t> vFoundGlobalBC(cols.size(), 0);                                      // global BCs for collisions
    std::vector<int> vTracksITS567perColl(cols.size(), 0);                                    // counter of tracks per found bc for occupancy studies
    std::vector<bool> vIsFullInfoForOccupancy(cols.size(), 0);                                // info for occupancy in +/- windows is available (i.e. a given coll is not too close to the TF borders)
    const float timeWinOccupancyCalcMinNS = confTimeIntervalForOccupancyCalculationMin * 1e3; // ns
    const float timeWinOccupancyCalcMaxNS = confTimeIntervalForOccupancyCalculationMax * 1e3; // ns
    // const double timeWinOccupancyExclusionRangeNS = confExclusionIntervalForOccupancyCalculation * 1e3; // ns
    const double bcNS = o2::constants::lhc::LHCBunchSpacingNS;

    // loop to find nearest bc with FT0 entry -> foundBC index
    for (auto& col : cols) {
      auto bc = col.bc_as<BCsWithBcSelsRun3>();
      int64_t meanBC = bc.globalBC();
      const double bcNS = o2::constants::lhc::LHCBunchSpacingNS;
      int64_t deltaBC = std::ceil(col.collisionTimeRes() / bcNS * 4);

      // count tracks of different types
      int nITS567cls = 0;
      int nITSTPCtracks = 0;
      int nTOFtracks = 0;
      int nTRDtracks = 0;
      double timeFromTOFtracks = 0;
      auto tracksGrouped = tracks.sliceBy(perCollision, col.globalIndex());
      for (auto& track : tracksGrouped) {
        if (!track.isPVContributor()) {
          continue;
        }
        nITSTPCtracks += track.hasITS() && track.hasTPC();
        nTOFtracks += track.hasTOF();
        nTRDtracks += track.hasTRD();
        // calculate average time using TOF tracks
        if (track.hasTOF()) {
          timeFromTOFtracks += track.trackTime();
        }

        if (track.itsNCls() >= 5)
          nITS567cls++;
      }
      LOGP(debug, "nContrib={} nITSTPCtracks={} nTOFtracks={} nTRDtracks={}", col.numContrib(), nITSTPCtracks, nTOFtracks, nTRDtracks);

      if (nTOFtracks > 0) {
        meanBC += TMath::FloorNint(timeFromTOFtracks / nTOFtracks / bcNS); // assign collision bc using TOF-matched tracks
        deltaBC = 4;                                                       // use precise bc from TOF tracks with +/-4 bc margin
      } else if (nITSTPCtracks > 0) {
        deltaBC += 30; // extend deltaBC for collisions built with ITS-TPC tracks only
      }

      int64_t minBC = meanBC - deltaBC;
      int64_t maxBC = meanBC + deltaBC;

      int32_t indexClosestTVX = findClosest(meanBC, mapGlobalBcWithTVX);
      int64_t tvxBC = bcs.iteratorAt(indexClosestTVX).globalBC();
      if (tvxBC >= minBC && tvxBC <= maxBC) { // closest TVX within search region
        bc.setCursor(indexClosestTVX);
      } else { // no TVX within search region, searching for TOR = T0A | T0C
        int32_t indexClosestTOR = findClosest(meanBC, mapGlobalBcWithTOR);
        int64_t torBC = bcs.iteratorAt(indexClosestTOR).globalBC();
        if (torBC >= minBC && torBC <= maxBC) {
          bc.setCursor(indexClosestTOR);
        }
      }
      int32_t foundBC = bc.globalIndex();
      int32_t colIndex = col.globalIndex();
      LOGP(debug, "foundBC = {} globalBC = {}", foundBC, bc.globalBC());
      vFoundBCindex[colIndex] = foundBC;
      vIsVertexITSTPC[colIndex] = nITSTPCtracks > 0;
      vIsVertexTOFmatched[colIndex] = nTOFtracks > 0;
      vIsVertexTRDmatched[colIndex] = nTRDtracks > 0;
      vCollisionsPerBc[foundBC]++;
      vTracksITS567perColl[colIndex] = nITS567cls;
      vFoundGlobalBC[colIndex] = bc.globalBC();

      // check that this collision has full information inside the time window (taking into account TF borders)
      int64_t bcInTF = (bc.globalBC() - bcSOR) % nBCsPerTF;
      vIsFullInfoForOccupancy[colIndex] = ((bcInTF - 300) * bcNS > -timeWinOccupancyCalcMinNS) && ((nBCsPerTF - 4000 - bcInTF) * bcNS > timeWinOccupancyCalcMaxNS) ? true : false;
    }

    // save indices of collisions in time range for occupancy calculation
    std::vector<std::vector<int>> vCollsInTimeWin;
    std::vector<std::vector<float>> vTimeDeltaForColls; // delta time wrt a given collision
    for (auto& col : cols) {
      int32_t colIndex = col.globalIndex();
      std::vector<int> vAssocToThisCol;
      std::vector<float> vCollsTimeDeltaWrtGivenColl;

      // protection against the TF borders
      if (!vIsFullInfoForOccupancy[colIndex]) {
        vCollsInTimeWin.push_back(vAssocToThisCol);
        vTimeDeltaForColls.push_back(vCollsTimeDeltaWrtGivenColl);
        continue;
      }

      int64_t foundGlobalBC = vFoundGlobalBC[colIndex];
      int64_t TFid = (foundGlobalBC - bcSOR) / nBCsPerTF;

      // find all collisions in time window before the current one (start with the current collision)
      int32_t minColIndex = colIndex;
      while (minColIndex >= 0) {
        int64_t thisBC = vFoundGlobalBC[minColIndex];
        // check if this is still the same TF
        int64_t thisTFid = (thisBC - bcSOR) / nBCsPerTF;
        if (thisTFid != TFid)
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
      int32_t maxColIndex = colIndex + 1;
      while (maxColIndex < cols.size()) {
        int64_t thisBC = vFoundGlobalBC[maxColIndex];
        int64_t thisTFid = (thisBC - bcSOR) / nBCsPerTF;
        if (thisTFid != TFid)
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

    // perform the occupancy calculation in the pre-defined time window
    std::vector<int> vNumTracksITS567inFullTimeWin(cols.size(), 0); // counter of tracks in full time window for occupancy studies
    std::vector<bool> vNoOccupAggressiveCuts(cols.size(), 0);       // no occupancy according to the agressive cuts
    std::vector<bool> vNoOccupStrictCuts(cols.size(), 0);           // no occupancy according to the strict cuts
    std::vector<bool> vNoOccupMediumCuts(cols.size(), 0);           // no occupancy according to the medium cuts
    std::vector<bool> vNoOccupRelaxedCuts(cols.size(), 0);          // no occupancy according to the relaxed cuts
    std::vector<bool> vNoOccupGentleCuts(cols.size(), 0);           // no occupancy according to the gentle cuts
    std::vector<bool> vNoCollInTimeRangeStandard(cols.size(), 0);   // no collisions in a specified time range
    std::vector<bool> vNoCollInTimeRangeNarrow(cols.size(), 0);     // no collisions in a specified time range

    // time ranges for occupancy calculation
    const int nTimeIntervals = 6;
    const float coeffOccupInTimeBins[] = {0.2, 0.4, 0.6, 1.};

    for (auto& col : cols) {
      int32_t colIndex = col.globalIndex();
      // protection against TF borders
      if (!vIsFullInfoForOccupancy[colIndex]) {
        vNumTracksITS567inFullTimeWin[colIndex] = -1; // occupancy in undefined (too close to TF borders)
        continue;
      }
      std::vector<int> vAssocToThisCol = vCollsInTimeWin[colIndex];
      std::vector<float> vCollsTimeDeltaWrtGivenColl = vTimeDeltaForColls[colIndex];
      int nITS567tracksInFullTimeWindow = 0;
      int nITS567tracksInTimeBins[nTimeIntervals] = {};
      int nITS567tracksForVetoStandard = 0; // to veto events with nearby collisions
      int nITS567tracksForVetoNarrow = 0;   // to veto events with nearby collisions (narrower range)
      for (int iCol = 0; iCol < vAssocToThisCol.size(); iCol++) {
        int thisColIndex = vAssocToThisCol[iCol];
        if (thisColIndex == colIndex) // skip the same collision
          continue;
        float dt = vCollsTimeDeltaWrtGivenColl[iCol] / 1e3; // ns -> us

        if (!confUseWeightsForOccupancyVariable) {
          nITS567tracksInFullTimeWindow += vTracksITS567perColl[thisColIndex];
        } else {
          // weighted occupancy
          float wOccup = 0;
          if (dt >= -40 && dt < -5) // collisions in the past
            wOccup = 1. / 1225 * (dt + 40) * (dt + 40);
          else if (dt >= -5 && dt < 15) // collisions near a given one
            wOccup = 1;
          else if (dt >= 15 && dt < 100) // collisions from the future
            wOccup = -1. / 85 * dt + 20. / 17;
          if (wOccup > 0)
            nITS567tracksInFullTimeWindow += wOccup * vTracksITS567perColl[thisColIndex];
        }

        for (int iTime = 0; iTime < nTimeIntervals; iTime++) {
          if (confTimeBinsForOccupancyCalculation->at(iTime) < dt && dt <= confTimeBinsForOccupancyCalculation->at(iTime + 1))
            nITS567tracksInTimeBins[iTime] += vTracksITS567perColl[thisColIndex];
          if (fabs(dt) < confTimeRangeVetoOnCollStandard)
            nITS567tracksForVetoStandard += vTracksITS567perColl[thisColIndex];
          if (fabs(dt) < confTimeRangeVetoOnCollNarrow)
            nITS567tracksForVetoNarrow += vTracksITS567perColl[thisColIndex];
        }
      }
      vNumTracksITS567inFullTimeWin[colIndex] = nITS567tracksInFullTimeWindow; // occupancy (without a current collision)

      // decisions based on occupancies in time bins
      bool decisions[4];
      for (int iCut = 0; iCut < 4; iCut++) {
        decisions[iCut] = true;
        for (int iTime = 0; iTime < nTimeIntervals; iTime++) {
          if (nITS567tracksInTimeBins[iTime] >= coeffOccupInTimeBins[iCut] * confReferenceOccupanciesInTimeBins->at(iTime)) {
            decisions[iCut] = false;
            break;
          }
        }
      }
      vNoOccupStrictCuts[colIndex] = decisions[0];
      vNoOccupMediumCuts[colIndex] = decisions[1];
      vNoOccupRelaxedCuts[colIndex] = decisions[2];
      vNoOccupGentleCuts[colIndex] = decisions[3];
      vNoOccupAggressiveCuts[colIndex] = ((nITS567tracksInTimeBins[0] < 300) && (nITS567tracksInTimeBins[1] == 0) && (nITS567tracksInTimeBins[2] == 0) && (nITS567tracksInTimeBins[3] == 0) && (nITS567tracksInTimeBins[4] < 200) && (nITS567tracksInTimeBins[5] < 400));
      vNoCollInTimeRangeStandard[colIndex] = (nITS567tracksForVetoStandard == 0);
      vNoCollInTimeRangeNarrow[colIndex] = (nITS567tracksForVetoNarrow == 0);
    }

    for (auto& col : cols) {
      int32_t colIndex = col.globalIndex();
      int32_t foundBC = vFoundBCindex[colIndex];
      auto bc = bcs.iteratorAt(foundBC);
      int32_t foundFT0 = bc.foundFT0Id();
      int32_t foundFV0 = bc.foundFV0Id();
      int32_t foundFDD = bc.foundFDDId();
      int32_t foundZDC = bc.foundZDCId();

      // compare zVtx from FT0 and from PV
      bool isGoodZvtxFT0vsPV = bc.has_foundFT0() ? fabs(bc.foundFT0().posZ() - col.posZ()) < maxDiffZvtxFT0vsPV : 0;

      // copy alias decisions from bcsel table
      uint32_t alias = bc.alias_raw();

      // copy selection decisions from bcsel table
      uint64_t selection = bc.selection_raw();
      selection |= vCollisionsPerBc[foundBC] <= 1 ? BIT(kNoSameBunchPileup) : 0;
      selection |= vIsVertexITSTPC[colIndex] ? BIT(kIsVertexITSTPC) : 0;
      selection |= vIsVertexTOFmatched[colIndex] ? BIT(kIsVertexTOFmatched) : 0;
      selection |= vIsVertexTRDmatched[colIndex] ? BIT(kIsVertexTRDmatched) : 0;
      selection |= isGoodZvtxFT0vsPV ? BIT(kIsGoodZvtxFT0vsPV) : 0;

      // selection bits based on occupancy pattern
      selection |= vNoOccupAggressiveCuts[colIndex] ? BIT(kNoHighOccupancyAgressive) : 0;
      selection |= vNoOccupStrictCuts[colIndex] && vNoCollInTimeRangeStandard[colIndex] ? BIT(kNoHighOccupancyStrict) : 0;
      selection |= vNoOccupMediumCuts[colIndex] && vNoCollInTimeRangeStandard[colIndex] ? BIT(kNoHighOccupancyMedium) : 0;
      selection |= vNoOccupRelaxedCuts[colIndex] && vNoCollInTimeRangeStandard[colIndex] ? BIT(kNoHighOccupancyRelaxed) : 0;
      selection |= vNoOccupGentleCuts[colIndex] && vNoCollInTimeRangeNarrow[colIndex] ? BIT(kNoHighOccupancyGentle) : 0;
      selection |= vNoCollInTimeRangeStandard[colIndex] ? BIT(kNoCollInTimeRangeStandard) : 0;
      selection |= vNoCollInTimeRangeNarrow[colIndex] ? BIT(kNoCollInTimeRangeNarrow) : 0;

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

      int nTracksITS567inFullTimeWin = vNumTracksITS567inFullTimeWin[colIndex];
      // histos.get<TH1>(HIST("hOccupancy"))->Fill(nTracksITS567inFullTimeWin);

      evsel(alias, selection, sel7, sel8, foundBC, foundFT0, foundFV0, foundFDD, foundZDC, nTracksITS567inFullTimeWin);
    }
  }

  PROCESS_SWITCH(EventSelectionTask, processRun3, "Process Run3 event selection", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  // Parse the metadata
  metadataInfo.initMetadata(cfgc);

  return WorkflowSpec{
    adaptAnalysisTask<BcSelectionTask>(cfgc),
    adaptAnalysisTask<EventSelectionTask>(cfgc)};
}
