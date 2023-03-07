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
using namespace o2;
using namespace o2::framework;

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
#include "TH1F.h"
using namespace evsel;

using BCsWithRun2InfosTimestampsAndMatches = soa::Join<aod::BCs, aod::Run2BCInfos, aod::Timestamps, aod::Run2MatchedToBCSparse>;
using BCsWithRun3Matchings = soa::Join<aod::BCs, aod::Timestamps, aod::Run3MatchedToBCSparse>;
using BCsWithBcSels = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels>;

struct BcSelectionTask {
  Produces<aod::BcSels> bcsel;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {
    // ccdb->setURL("http://ccdb-test.cern.ch:8080");
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    histos.add("hCounterTVX", "", kTH1F, {{1, 0., 1.}});
  }

  void processRun2(
    BCsWithRun2InfosTimestampsAndMatches const& bcs,
    aod::Zdcs const&,
    aod::FV0As const&,
    aod::FV0Cs const&,
    aod::FT0s const&,
    aod::FDDs const&)
  {

    for (auto& bc : bcs) {
      EventSelectionParams* par = ccdb->getForTimeStamp<EventSelectionParams>("EventSelection/EventSelectionParams", bc.timestamp());
      TriggerAliases* aliases = ccdb->getForTimeStamp<TriggerAliases>("EventSelection/TriggerAliases", bc.timestamp());
      // fill fired aliases
      int32_t alias[kNaliases] = {0};
      uint64_t triggerMask = bc.triggerMask();
      for (auto& al : aliases->GetAliasToTriggerMaskMap()) {
        alias[al.first] |= (triggerMask & al.second) > 0;
      }
      uint64_t triggerMaskNext50 = bc.triggerMaskNext50();
      for (auto& al : aliases->GetAliasToTriggerMaskNext50Map()) {
        alias[al.first] |= (triggerMaskNext50 & al.second) > 0;
      }
      alias[kALL] = 1;

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

      // applying timing selections
      bool bbV0A = timeV0A > par->fV0ABBlower && timeV0A < par->fV0ABBupper;
      bool bbV0C = timeV0C > par->fV0CBBlower && timeV0C < par->fV0CBBupper;
      bool bbFDA = timeFDA > par->fFDABBlower && timeFDA < par->fFDABBupper;
      bool bbFDC = timeFDC > par->fFDCBBlower && timeFDC < par->fFDCBBupper;
      bool bgV0A = timeV0A > par->fV0ABGlower && timeV0A < par->fV0ABGupper;
      bool bgV0C = timeV0C > par->fV0CBGlower && timeV0C < par->fV0CBGupper;
      bool bgFDA = timeFDA > par->fFDABGlower && timeFDA < par->fFDABGupper;
      bool bgFDC = timeFDC > par->fFDCBGlower && timeFDC < par->fFDCBGupper;
      float znSum = timeZNA + timeZNC;
      float znDif = timeZNA - timeZNC;

      // fill time-based selection criteria
      int32_t selection[kNsel] = {0}; // TODO switch to bool array
      selection[kIsBBV0A] = bbV0A;
      selection[kIsBBV0C] = bbV0C;
      selection[kIsBBFDA] = bbFDA;
      selection[kIsBBFDC] = bbFDC;
      selection[kNoBGV0A] = !bgV0A;
      selection[kNoBGV0C] = !bgV0C;
      selection[kNoBGFDA] = !bgFDA;
      selection[kNoBGFDC] = !bgFDC;
      selection[kIsBBT0A] = timeT0A > par->fT0ABBlower && timeT0A < par->fT0ABBupper;
      selection[kIsBBT0C] = timeT0C > par->fT0CBBlower && timeT0C < par->fT0CBBupper;
      selection[kIsBBZNA] = timeZNA > par->fZNABBlower && timeZNA < par->fZNABBupper;
      selection[kIsBBZNC] = timeZNC > par->fZNCBBlower && timeZNC < par->fZNCBBupper;
      selection[kNoBGZNA] = !(fabs(timeZNA) > par->fZNABGlower && fabs(timeZNA) < par->fZNABGupper);
      selection[kNoBGZNC] = !(fabs(timeZNC) > par->fZNCBGlower && fabs(timeZNC) < par->fZNCBGupper);
      selection[kIsBBZAC] = pow((znSum - par->fZNSumMean) / par->fZNSumSigma, 2) + pow((znDif - par->fZNDifMean) / par->fZNDifSigma, 2) < 1;

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
      uint32_t spdClusters = bc.spdClustersL0() + bc.spdClustersL1();

      // Calculate pileup and background related selection flags
      // V0A0 excluded from online V0A charge sum => excluding also from offline sum for consistency
      float ofV0M = multFV0A + multFV0C - multRingV0A[0];
      float onV0M = bc.v0TriggerChargeA() + bc.v0TriggerChargeC();
      float ofSPD = bc.spdFiredChipsL0() + bc.spdFiredChipsL1();
      float onSPD = bc.spdFiredFastOrL0() + bc.spdFiredFastOrL1();
      float multV0C012 = multRingV0C[0] + multRingV0C[1] + multRingV0C[2];

      selection[kNoV0MOnVsOfPileup] = onV0M > par->fV0MOnVsOfA + par->fV0MOnVsOfB * ofV0M;
      selection[kNoSPDOnVsOfPileup] = onSPD > par->fSPDOnVsOfA + par->fSPDOnVsOfB * ofSPD;
      selection[kNoV0Casymmetry] = multRingV0C[3] > par->fV0CasymA + par->fV0CasymB * multV0C012;

      // copy remaining selection decisions from eventCuts
      uint32_t eventCuts = bc.eventCuts();

      selection[kIsGoodTimeRange] = (eventCuts & 1 << aod::kTimeRangeCut) > 0;
      selection[kNoIncompleteDAQ] = (eventCuts & 1 << aod::kIncompleteDAQ) > 0;
      selection[kNoTPCLaserWarmUp] = (eventCuts & 1 << aod::kIsTPCLaserWarmUp) == 0;
      selection[kNoTPCHVdip] = (eventCuts & 1 << aod::kIsTPCHVdip) == 0;
      selection[kNoPileupFromSPD] = (eventCuts & 1 << aod::kIsPileupFromSPD) == 0;
      selection[kNoV0PFPileup] = (eventCuts & 1 << aod::kIsV0PFPileup) == 0;
      selection[kNoInconsistentVtx] = (eventCuts & 1 << aod::kConsistencySPDandTrackVertices) > 0;
      selection[kNoPileupInMultBins] = (eventCuts & 1 << aod::kPileupInMultBins) > 0;
      selection[kNoPileupMV] = (eventCuts & 1 << aod::kPileUpMV) > 0;
      selection[kNoPileupTPC] = (eventCuts & 1 << aod::kTPCPileUp) > 0;
      selection[kIsTriggerTVX] = bc.has_ft0() ? (bc.ft0().triggerMask() & BIT(o2::ft0::Triggers::bitVertex)) > 0 : 0;
      selection[kIsINT1] = bbV0A || bbV0C || ofSPD > 0;

      int32_t foundFT0 = bc.has_ft0() ? bc.ft0().globalIndex() : -1;
      int32_t foundFV0 = bc.has_fv0a() ? bc.fv0a().globalIndex() : -1;
      int32_t foundFDD = bc.has_fdd() ? bc.fdd().globalIndex() : -1;
      int32_t foundZDC = bc.has_zdc() ? bc.zdc().globalIndex() : -1;

      // Fill TVX (T0 vertex) counters
      if (selection[kIsTriggerTVX]) {
        histos.get<TH1>(HIST("hCounterTVX"))->Fill(Form("%d", bc.runNumber()), 1);
      }

      // Fill bc selection columns
      bcsel(alias, selection,
            bbV0A, bbV0C, bgV0A, bgV0C,
            bbFDA, bbFDC, bgFDA, bgFDC,
            multRingV0A, multRingV0C, spdClusters, foundFT0, foundFV0, foundFDD, foundZDC);
    }
  }
  PROCESS_SWITCH(BcSelectionTask, processRun2, "Process Run2 event selection", true);

  void processRun3(BCsWithRun3Matchings const& bcs,
                   aod::Zdcs const& zdcs,
                   aod::FV0As const&,
                   aod::FT0s const&,
                   aod::FDDs const&)
  {
    // map from GlobalBC to BcId needed to find triggerBc
    std::map<uint64_t, int32_t> mapGlobalBCtoBcId;
    for (auto& bc : bcs) {
      mapGlobalBCtoBcId[bc.globalBC()] = bc.globalIndex();
    }

    for (auto bc : bcs) {
      EventSelectionParams* par = ccdb->getForTimeStamp<EventSelectionParams>("EventSelection/EventSelectionParams", bc.timestamp());
      TriggerAliases* aliases = ccdb->getForTimeStamp<TriggerAliases>("EventSelection/TriggerAliases", bc.timestamp());
      int32_t alias[kNaliases] = {0};

      // workaround for pp2022 apass2-apass3 (trigger info is shifted by -294 bcs)
      int32_t triggerBcId = mapGlobalBCtoBcId[bc.globalBC() - 294];
      if (triggerBcId) {
        auto triggerBc = bcs.iteratorAt(triggerBcId);
        uint64_t triggerMask = triggerBc.triggerMask();
        for (auto& al : aliases->GetAliasToTriggerMaskMap()) {
          alias[al.first] |= (triggerMask & al.second) > 0;
        }
      }
      alias[kALL] = 1;

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
      float znSum = timeZNA + timeZNC;
      float znDif = timeZNA - timeZNC;

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

      // applying timing selections
      bool bbV0A = timeV0A > par->fV0ABBlower && timeV0A < par->fV0ABBupper;
      bool bbFDA = timeFDA > par->fFDABBlower && timeFDA < par->fFDABBupper;
      bool bbFDC = timeFDC > par->fFDCBBlower && timeFDC < par->fFDCBBupper;
      bool bgV0A = timeV0ABG > par->fV0ABGlower && timeV0ABG < par->fV0ABGupper;
      bool bgFDA = timeFDABG > par->fFDABGlower && timeFDABG < par->fFDABGupper;
      bool bgFDC = timeFDCBG > par->fFDCBGlower && timeFDCBG < par->fFDCBGupper;
      bool bgT0A = timeT0ABG > par->fT0ABGlower && timeT0ABG < par->fT0ABGupper;
      bool bgT0C = timeT0CBG > par->fT0CBGlower && timeT0CBG < par->fT0CBGupper;
      bool bbV0C = 0;
      bool bgV0C = 0;

      // fill time-based selection criteria
      int32_t selection[kNsel] = {0}; // TODO switch to bool array
      selection[kIsBBV0A] = bbV0A;
      selection[kIsBBFDA] = bbFDA;
      selection[kIsBBFDC] = bbFDC;
      selection[kNoBGV0A] = !bgV0A;
      selection[kNoBGFDA] = !bgFDA;
      selection[kNoBGFDC] = !bgFDC;
      selection[kNoBGT0A] = !bgT0A;
      selection[kNoBGT0C] = !bgT0C;
      selection[kIsBBT0A] = timeT0A > par->fT0ABBlower && timeT0A < par->fT0ABBupper;
      selection[kIsBBT0C] = timeT0C > par->fT0CBBlower && timeT0C < par->fT0CBBupper;
      selection[kIsBBZNA] = timeZNA > par->fZNABBlower && timeZNA < par->fZNABBupper;
      selection[kIsBBZNC] = timeZNC > par->fZNCBBlower && timeZNC < par->fZNCBBupper;
      selection[kIsBBZAC] = pow((znSum - par->fZNSumMean) / par->fZNSumSigma, 2) + pow((znDif - par->fZNDifMean) / par->fZNDifSigma, 2) < 1;
      selection[kNoBGZNA] = !(fabs(timeZNA) > par->fZNABGlower && fabs(timeZNA) < par->fZNABGupper);
      selection[kNoBGZNC] = !(fabs(timeZNC) > par->fZNCBGlower && fabs(timeZNC) < par->fZNCBGupper);
      selection[kIsTriggerTVX] = bc.has_ft0() ? (bc.ft0().triggerMask() & BIT(o2::ft0::Triggers::bitVertex)) > 0 : 0;

      // Calculate V0 multiplicity per ring
      float multRingV0A[5] = {0.};
      float multRingV0C[4] = {0.};
      if (bc.has_fv0a()) {
        for (unsigned int i = 0; i < bc.fv0a().amplitude().size(); ++i) {
          int ring = bc.fv0a().channel()[i] / 8;
          if (ring == 5) {
            ring = 4; // Outermost ring has 16 channels
          }
          multRingV0A[ring] += bc.fv0a().amplitude()[i];
        }
      }

      uint32_t spdClusters = 0;

      int32_t foundFT0 = bc.has_ft0() ? bc.ft0().globalIndex() : -1;
      int32_t foundFV0 = bc.has_fv0a() ? bc.fv0a().globalIndex() : -1;
      int32_t foundFDD = bc.has_fdd() ? bc.fdd().globalIndex() : -1;
      int32_t foundZDC = bc.has_zdc() ? bc.zdc().globalIndex() : -1;
      LOGP(debug, "foundFT0={}\n", foundFT0);

      // Fill TVX (T0 vertex) counters
      if (selection[kIsTriggerTVX]) {
        histos.get<TH1>(HIST("hCounterTVX"))->Fill(Form("%d", bc.runNumber()), 1);
      }

      // Fill bc selection columns
      bcsel(alias, selection,
            bbV0A, bbV0C, bgV0A, bgV0C,
            bbFDA, bbFDC, bgFDA, bgFDC,
            multRingV0A, multRingV0C, spdClusters, foundFT0, foundFV0, foundFDD, foundZDC);
    }
  }
  PROCESS_SWITCH(BcSelectionTask, processRun3, "Process Run3 event selection", false);
};

struct EventSelectionTask {
  Produces<aod::EvSels> evsel;
  Configurable<std::string> syst{"syst", "PbPb", "pp, pPb, Pbp, PbPb, XeXe"}; // TODO determine from AOD metadata or from CCDB
  Configurable<int> muonSelection{"muonSelection", 0, "0 - barrel, 1 - muon selection with pileup cuts, 2 - muon selection without pileup cuts"};
  Configurable<int> customDeltaBC{"customDeltaBC", 300, "custom BC delta for FIT-collision matching"};
  Configurable<bool> isMC{"isMC", 0, "0 - data, 1 - MC"};
  Partition<aod::Tracks> tracklets = (aod::track::trackType == static_cast<uint8_t>(o2::aod::track::TrackTypeEnum::Run2Tracklet));

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {
    // ccdb->setURL("http://ccdb-test.cern.ch:8080");
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    histos.add("hColCounterAll", "", kTH1F, {{1, 0., 1.}});
    histos.add("hColCounterAcc", "", kTH1F, {{1, 0., 1.}});
  }

  void processRun2(aod::Collision const& col, BCsWithBcSels const& bcs, aod::Tracks const& tracks)
  {
    auto bc = col.bc_as<BCsWithBcSels>();
    EventSelectionParams* par = ccdb->getForTimeStamp<EventSelectionParams>("EventSelection/EventSelectionParams", bc.timestamp());
    bool* applySelection = par->GetSelection(muonSelection);
    if (isMC) {
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
    int32_t alias[kNaliases];
    for (int i = 0; i < kNaliases; i++) {
      alias[i] = bc.alias()[i];
    }

    // copy selection decisions from bcsel table
    int32_t selection[kNsel] = {0};
    for (int i = 0; i < kNsel; i++) {
      selection[i] = bc.selection()[i];
    }

    // copy multiplicity per ring and calculate V0C012 multiplicity
    float multRingV0A[5] = {0.};
    float multRingV0C[4] = {0.};
    for (int i = 0; i < 5; i++) {
      multRingV0A[i] = bc.multRingV0A()[i];
    }
    for (int i = 0; i < 4; i++) {
      multRingV0C[i] = bc.multRingV0C()[i];
    }
    float multV0C012 = bc.multRingV0C()[0] + bc.multRingV0C()[1] + bc.multRingV0C()[2];

    // applying selections depending on the number of tracklets
    auto trackletsGrouped = tracklets->sliceByCached(aod::track::collisionId, col.globalIndex());
    int nTkl = trackletsGrouped.size();

    uint32_t spdClusters = bc.spdClusters();
    selection[kNoSPDClsVsTklBG] = spdClusters < par->fSPDClsVsTklA + nTkl * par->fSPDClsVsTklB;
    selection[kNoV0C012vsTklBG] = !(nTkl < 6 && multV0C012 > par->fV0C012vsTklA + nTkl * par->fV0C012vsTklB);

    // copy beam-beam and beam-gas flags from bcsel table
    bool bbV0A = bc.bbV0A();
    bool bbV0C = bc.bbV0C();
    bool bgV0A = bc.bgV0A();
    bool bgV0C = bc.bgV0C();
    bool bbFDA = bc.bbFDA();
    bool bbFDC = bc.bbFDC();
    bool bgFDA = bc.bgFDA();
    bool bgFDC = bc.bgFDC();

    // apply int7-like selections
    bool sel7 = 1;
    for (int i = 0; i < kNsel; i++) {
      sel7 &= applySelection[i] ? selection[i] : 1;
    }

    // TODO introduce array of sel[0]... sel[8] or similar?
    bool sel8 = selection[kIsBBT0A] & selection[kIsBBT0C]; // TODO apply other cuts for sel8
    bool sel1 = selection[kIsINT1] & selection[kNoBGV0A] & selection[kNoBGV0C] & selection[kNoTPCLaserWarmUp] & selection[kNoTPCHVdip];

    // INT1 (SPDFO>0 | V0A | V0C) minimum bias trigger logic used in pp2010 and pp2011
    bool isINT1period = bc.runNumber() <= 136377 || (bc.runNumber() >= 144871 && bc.runNumber() <= 159582);

    // fill counters
    if (isMC || (!isINT1period && alias[kINT7]) || (isINT1period && alias[kINT1])) {
      histos.get<TH1>(HIST("hColCounterAll"))->Fill(Form("%d", bc.runNumber()), 1);
      if ((!isINT1period && sel7) || (isINT1period && sel1)) {
        histos.get<TH1>(HIST("hColCounterAcc"))->Fill(Form("%d", bc.runNumber()), 1);
      }
    }

    evsel(alias, selection,
          bbV0A, bbV0C, bgV0A, bgV0C,
          bbFDA, bbFDC, bgFDA, bgFDC,
          multRingV0A, multRingV0C, spdClusters, nTkl, sel7, sel8,
          foundBC, foundFT0, foundFV0, foundFDD, foundZDC);
  }
  PROCESS_SWITCH(EventSelectionTask, processRun2, "Process Run2 event selection", true);

  void processRun3(aod::Collision const& col, soa::Join<aod::TracksIU, aod::TracksExtra> const& tracks, BCsWithBcSels const& bcs)
  {
    // count tracks of different types
    int nITStracks = 0;
    int nTPCtracks = 0;
    int nTOFtracks = 0;
    int nTRDtracks = 0;
    for (auto& track : tracks) {
      if (!track.isPVContributor()) {
        continue;
      }
      nITStracks += track.hasITS();
      nTPCtracks += track.hasTPC();
      nTOFtracks += track.hasTOF();
      nTRDtracks += track.hasTRD();
    }

    LOGP(debug, "nContrib={} nITStracks={} nTPCtracks={} nTOFtracks={} nTRDtracks={}", col.numContrib(), nITStracks, nTPCtracks, nTOFtracks, nTRDtracks);

    auto bc = col.bc_as<BCsWithBcSels>();
    int run = bc.runNumber();

    int64_t meanBC = bc.globalBC();
    int64_t deltaBC = std::ceil(col.collisionTimeRes() / o2::constants::lhc::LHCBunchSpacingNS * 4);

    // pilot runs with isolated bunches: using custom delta bc to improve FT0-collision matching
    if (run <= 520297) {
      int min_bunch_spacing = (run >= 520259 && run <= 520297) || (run >= 519041 && run <= 519507) ? 100 : 1200;
      deltaBC = deltaBC < min_bunch_spacing ? min_bunch_spacing / 2 : deltaBC;
    }

    // use custom delta
    if (customDeltaBC > 0) {
      deltaBC = customDeltaBC;
    }

    uint64_t minBC = meanBC - deltaBC;
    uint64_t maxBC = meanBC + deltaBC;

    // temporary workaround for runs at high rate (>22m)
    // significant eta-dependent biases up to 70 bcs for single ITS-TPC track times
    // extend deltaBC for collisions built with ITS-TPC tracks only
    if (run >= 523141 && nTRDtracks == 0 && nTOFtracks == 0 && nTPCtracks > 0) {
      minBC -= 100;
      maxBC += 100;
    }

    // quick fix to account for reduction of TVX efficiency at high n contibutors in LHC22s
    if (run >= 529397 && run <= 529418 && nTRDtracks == 0 && nTOFtracks > 0) {
      minBC -= 100;
      maxBC += 100;
    }

    // temporary workaround for runs without proper TOF calibration
    if (run > 520297 && run < 523306 && nTOFtracks > 0) {
      minBC = meanBC - deltaBC - 2;
      maxBC = meanBC + deltaBC - 2;
    }

    // precise timing for collisions with TRD-matched tracks
    if (nTRDtracks > 0) {
      minBC = meanBC;
      maxBC = meanBC;
      // collisions with TRD tracks shifted by 15 bcs in LHC22cdef
      if (run > 520297 && run <= 521326) {
        minBC = meanBC - 15;
        maxBC = meanBC - 15;
      }
      // collisions with TRD tracks shifted by -1 bc in LHC22s
      if (run >= 529397 && run <= 529418) {
        minBC = meanBC - 1;
        maxBC = meanBC + 1;
      }
    }

    int forwardMoveCount = 0, backwardMoveCount = 0;
    int forwardMoveCountTvx = 0, backwardMoveCountTvx = 0;
    uint64_t backwardBC = minBC - 1;
    uint64_t forwardBC = maxBC + 1;
    uint64_t backwardTvxBC = minBC - 1;
    uint64_t forwardTvxBC = maxBC + 1;

    LOGP(debug, "meanBC={} minBC={} maxBC={} collisionTimeRes={}", meanBC, minBC, maxBC, col.collisionTimeRes());

    // search TVX in forward direction starting from the current bc
    while (bc != bcs.end() && bc.globalBC() <= maxBC && bc.globalBC() >= minBC) {
      if (bc.selection()[kIsTriggerTVX]) {
        forwardTvxBC = bc.globalBC();
        break;
      }
      bc++;
      forwardMoveCountTvx++;
    }
    bc.moveByIndex(-forwardMoveCountTvx);

    // search TVX in backward direction
    while (bc.globalIndex() > 0 && bc.globalBC() >= minBC) {
      bc--;
      backwardMoveCountTvx--;
      if (bc.globalBC() > maxBC || bc.globalBC() < minBC) {
        continue;
      }
      if (bc.selection()[kIsTriggerTVX]) {
        backwardTvxBC = bc.globalBC();
        break;
      }
    }
    bc.moveByIndex(-backwardMoveCountTvx);

    // search FT0-OR in forward direction starting from the current bc
    while (bc != bcs.end() && bc.globalBC() <= maxBC && bc.globalBC() >= minBC) {
      if (bc.selection()[kIsBBT0A] || bc.selection()[kIsBBT0C]) {
        forwardBC = bc.globalBC();
        break;
      }
      bc++;
      forwardMoveCount++;
    }
    bc.moveByIndex(-forwardMoveCount);

    // search FT0-OR in backward direction
    while (bc.globalIndex() > 0 && bc.globalBC() >= minBC) {
      bc--;
      backwardMoveCount--;
      if (bc.globalBC() > maxBC || bc.globalBC() < minBC) {
        continue;
      }
      if (bc.selection()[kIsBBT0A] || bc.selection()[kIsBBT0C]) {
        backwardBC = bc.globalBC();
        break;
      }
    }
    bc.moveByIndex(-backwardMoveCount);

    // first check for found TVX signal. If TVX is not found, search for FT0-OR
    if (forwardTvxBC <= maxBC && backwardTvxBC >= minBC) {
      // if TVX is found on both sides from meanBC, move to closest one
      if (labs(int64_t(forwardTvxBC) - meanBC) < labs(int64_t(backwardTvxBC) - meanBC)) {
        bc.moveByIndex(forwardMoveCountTvx);
      } else {
        bc.moveByIndex(backwardMoveCountTvx);
      }
    } else if (forwardTvxBC <= maxBC) {
      // if TVX is found only in forward, move forward
      bc.moveByIndex(forwardMoveCountTvx);
    } else if (backwardTvxBC >= minBC) {
      // if TVX is found only in backward, move backward
      bc.moveByIndex(backwardMoveCountTvx);
    } else if (forwardBC <= maxBC && backwardBC >= minBC) {
      // if FT0-OR is found on both sides from meanBC, move to closest one
      if (labs(int64_t(forwardBC) - meanBC) < labs(int64_t(backwardBC) - meanBC)) {
        bc.moveByIndex(forwardMoveCount);
      } else {
        bc.moveByIndex(backwardMoveCount);
      }
    } else if (forwardBC <= maxBC) {
      // if FT0-OR is found only in forward, move forward
      bc.moveByIndex(forwardMoveCount);
    } else if (backwardBC >= minBC) {
      // if FT0-OR is found only in backward, move backward
      bc.moveByIndex(backwardMoveCount);
    }

    int32_t foundBC = bc.globalIndex();
    int32_t foundFT0 = bc.foundFT0Id();
    int32_t foundFV0 = bc.foundFV0Id();
    int32_t foundFDD = bc.foundFDDId();
    int32_t foundZDC = bc.foundZDCId();

    LOGP(debug, "foundFT0 = {} globalBC = {}", foundFT0, bc.globalBC());

    // copy alias decisions from bcsel table
    int32_t alias[kNaliases];
    for (int i = 0; i < kNaliases; i++) {
      alias[i] = bc.alias()[i];
    }

    // copy selection decisions from bcsel table
    int32_t selection[kNsel] = {0};
    for (int i = 0; i < kNsel; i++) {
      selection[i] = bc.selection()[i];
    }

    // copy multiplicity per ring
    float multRingV0A[5] = {0.};
    float multRingV0C[4] = {0.};
    for (int i = 0; i < 5; i++) {
      multRingV0A[i] = bc.multRingV0A()[i];
    }
    for (int i = 0; i < 4; i++) {
      multRingV0C[i] = bc.multRingV0C()[i];
    }

    int nTkl = 0;
    uint32_t spdClusters = 0;

    // copy beam-beam and beam-gas flags from bcsel table
    bool bbV0A = bc.bbV0A();
    bool bbV0C = bc.bbV0C();
    bool bgV0A = bc.bgV0A();
    bool bgV0C = bc.bgV0C();
    bool bbFDA = bc.bbFDA();
    bool bbFDC = bc.bbFDC();
    bool bgFDA = bc.bgFDA();
    bool bgFDC = bc.bgFDC();

    // apply int7-like selections
    bool sel7 = 0;

    // TODO apply other cuts for sel8
    // TODO introduce sel1 etc?
    // TODO introduce array of sel[0]... sel[8] or similar?
    bool sel8 = selection[kIsBBT0A] & selection[kIsBBT0C];

    // fill counters
    histos.get<TH1>(HIST("hColCounterAll"))->Fill(Form("%d", bc.runNumber()), 1);
    if (sel8) {
      histos.get<TH1>(HIST("hColCounterAcc"))->Fill(Form("%d", bc.runNumber()), 1);
    }

    evsel(alias, selection,
          bbV0A, bbV0C, bgV0A, bgV0C,
          bbFDA, bbFDC, bgFDA, bgFDC,
          multRingV0A, multRingV0C, spdClusters, nTkl, sel7, sel8,
          foundBC, foundFT0, foundFV0, foundFDD, foundZDC);
  }
  PROCESS_SWITCH(EventSelectionTask, processRun3, "Process Run3 event selection", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<BcSelectionTask>(cfgc),
    adaptAnalysisTask<EventSelectionTask>(cfgc)};
}
