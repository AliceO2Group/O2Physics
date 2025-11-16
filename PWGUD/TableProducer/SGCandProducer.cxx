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
/// \file SGCandProducer.cxx
/// \brief Produces PWGUD derived table from standard tables
///
/// \author Alexander Bylinkin <roman.lavicka@cern.ch>, Uniersity of Bergen
/// \since  23.11.2023
/// \author Adam Matyja <adam.tomasz.matyja@cern.ch>, INP PAN Krakow, Poland
/// \since  May 2025
//

#include "PWGUD/Core/SGSelector.h"
#include "PWGUD/Core/UPCHelpers.h"
#include "PWGUD/DataModel/UDTables.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/ctpRateFetcher.h"
#include "Common/DataModel/EventSelection.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/LHCConstants.h"
#include "DataFormatsFIT/Triggers.h"
#include "DataFormatsParameters/AggregatedRunInfo.h"
#include "DataFormatsParameters/GRPLHCIFData.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "Framework/ASoA.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Vertex.h"

#include <cmath>
#include <map>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::dataformats;
using namespace o2::aod::rctsel;

#define getHist(type, name) std::get<std::shared_ptr<type>>(histPointers[name])

struct SGCandProducer {
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  // data inputs
  using CCs = soa::Join<aod::Collisions, aod::EvSels>;
  using CC = CCs::iterator;
  using BCs = soa::Join<aod::BCsWithTimestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;
  using BC = BCs::iterator;
  using TCs = soa::Join<aod::Tracks, /*aod::TracksCov,*/ aod::TracksExtra, aod::TracksDCA, aod::TrackSelection,
                        aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr,
                        aod::pidTPCFullDe, aod::pidTPCFullTr, aod::pidTPCFullHe, aod::pidTPCFullAl,
                        aod::TOFSignal, aod::pidTOFbeta,
                        aod::pidTOFFullDe, aod::pidTOFFullTr, aod::pidTOFFullHe, aod::pidTOFFullAl,
                        aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>;
  using FWs = aod::FwdTracks;

  using MCCCs = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>;
  using MCCC = MCCCs::iterator;

  // get an SGCutparHolder
  SGCutParHolder sameCuts = SGCutParHolder(); // SGCutparHolder
  Configurable<SGCutParHolder> SGCuts{"SGCuts", {}, "SG event cuts"};
  Configurable<bool> verboseInfo{"verboseInfo", false, "Print general info to terminal; default it false."};
  Configurable<bool> saveAllTracks{"saveAllTracks", true, "save only PV contributors or all tracks associated to a collision"};
  Configurable<bool> savenonPVCITSOnlyTracks{"savenonPVCITSOnlyTracks", false, "save non PV contributors with ITS only information"};
  Configurable<bool> rejectAtTFBoundary{"rejectAtTFBoundary", true, "reject collisions at a TF boundary"};
  Configurable<bool> noITSROFrameBorder{"noITSROFrameBorder", true, "reject ITS RO Frame Border"};
  Configurable<bool> noSameBunchPileUp{"noSameBunchPileUp", true, "reject SameBunchPileUp"};
  Configurable<bool> IsGoodVertex{"IsGoodVertex", false, "Select FT0 PV vertex matching"};
  Configurable<bool> ITSTPCVertex{"ITSTPCVertex", true, "reject ITS-only vertex"}; // if one wants to look at Single Gap pp events
  Configurable<std::vector<int>> generatorIds{"generatorIds", std::vector<int>{-1}, "MC generatorIds to process"};
  Configurable<bool> storeSG{"storeSG", true, "Store SG events in the output"};
  Configurable<bool> storeDG{"storeDG", true, "Store DG events in the output"};

  Configurable<bool> isGoodRCTCollision{"isGoodRCTCollision", true, "Check RCT flags for FT0,ITS,TPC and tracking"};
  Configurable<bool> isGoodRCTZdc{"isGoodRCTZdc", false, "Check RCT flags for ZDC if present in run"};

  // Configurables to decide which tables are filled
  Configurable<bool> fillTrackTables{"fillTrackTables", true, "Fill track tables"};
  Configurable<bool> fillFwdTrackTables{"fillFwdTrackTables", true, "Fill forward track tables"};

  //  SG selector
  SGSelector sgSelector;
  ctpRateFetcher mRateFetcher;

  // initialize RCT flag checker
  RCTFlagsChecker myRCTChecker{"CBT"};

  // data tables
  Produces<aod::SGCollisions> outputSGCollisions;
  Produces<aod::UDCollisions> outputCollisions;
  Produces<aod::UDCollisionsSels> outputCollisionsSels;
  Produces<aod::UDCollisionSelExtras> outputCollisionSelExtras;
  Produces<aod::UDCollsLabels> outputCollsLabels;
  Produces<aod::UDZdcs> outputZdcs;
  Produces<aod::UDZdcsReduced> udZdcsReduced;
  Produces<aod::UDTracks> outputTracks;
  Produces<aod::UDTracksCov> outputTracksCov;
  Produces<aod::UDTracksDCA> outputTracksDCA;
  Produces<aod::UDTracksPID> outputTracksPID;
  Produces<aod::UDTracksPIDExtra> outputTracksPIDExtra;
  Produces<aod::UDTracksExtra> outputTracksExtra;
  Produces<aod::UDTracksFlags> outputTracksFlag;
  Produces<aod::UDFwdTracks> outputFwdTracks;
  Produces<aod::UDFwdTracksExtra> outputFwdTracksExtra;
  Produces<aod::UDTracksLabels> outputTracksLabel;

  // initialize histogram registry
  HistogramRegistry registry{
    "registry",
    {}};
  std::map<std::string, HistPtr> histPointers;

  int runNumber = -1;

  // function to update UDFwdTracks, UDFwdTracksExtra
  template <typename TFwdTrack>
  void updateUDFwdTrackTables(TFwdTrack const& fwdtrack, uint64_t const& bcnum)
  {
    outputFwdTracks(outputCollisions.lastIndex(),
                    fwdtrack.px(), fwdtrack.py(), fwdtrack.pz(), fwdtrack.sign(),
                    bcnum, fwdtrack.trackTime(), fwdtrack.trackTimeRes());
    outputFwdTracksExtra(fwdtrack.trackType(),
                         fwdtrack.nClusters(),
                         fwdtrack.pDca(),
                         fwdtrack.rAtAbsorberEnd(),
                         fwdtrack.chi2(),
                         fwdtrack.chi2MatchMCHMID(),
                         fwdtrack.chi2MatchMCHMFT(),
                         fwdtrack.mchBitMap(),
                         fwdtrack.midBitMap(),
                         fwdtrack.midBoards());
  }

  // function to update UDTracks, UDTracksCov, UDTracksDCA, UDTracksPID, UDTracksExtra, UDTracksFlag,
  // and UDTrackCollisionIDs
  template <typename TTrack>
  void updateUDTrackTables(int64_t lastIndex, TTrack const& track, uint64_t const& bcnum)
  {
    outputTracks(lastIndex,
                 track.px(), track.py(), track.pz(), track.sign(),
                 bcnum, track.trackTime(), track.trackTimeRes());

    // float sigmaY = track.sigmaY();
    // float sigmaZ = track.sigmaZ();
    float sigmaY = -1.;
    float sigmaZ = -1.;
    outputTracksCov(track.x(), track.y(), track.z(), sigmaY, sigmaZ);

    outputTracksDCA(track.dcaZ(), track.dcaXY());
    outputTracksPID(track.tpcNSigmaEl(),
                    track.tpcNSigmaMu(),
                    track.tpcNSigmaPi(),
                    track.tpcNSigmaKa(),
                    track.tpcNSigmaPr(),
                    track.beta(),
                    track.betaerror(),
                    track.tofNSigmaEl(),
                    track.tofNSigmaMu(),
                    track.tofNSigmaPi(),
                    track.tofNSigmaKa(),
                    track.tofNSigmaPr());
    outputTracksPIDExtra(track.tpcNSigmaDe(),
                         track.tpcNSigmaTr(),
                         track.tpcNSigmaHe(),
                         track.tpcNSigmaAl(),
                         track.tofNSigmaDe(),
                         track.tofNSigmaTr(),
                         track.tofNSigmaHe(),
                         track.tofNSigmaAl());
    outputTracksExtra(track.tpcInnerParam(),
                      track.itsClusterSizes(),
                      track.tpcNClsFindable(),
                      track.tpcNClsFindableMinusFound(),
                      track.tpcNClsFindableMinusCrossedRows(),
                      track.tpcNClsShared(),
                      track.trdPattern(),
                      track.itsChi2NCl(),
                      track.tpcChi2NCl(),
                      track.trdChi2(),
                      track.tofChi2(),
                      track.tpcSignal(),
                      track.tofSignal(),
                      track.trdSignal(),
                      track.length(),
                      track.tofExpMom(),
                      track.detectorMap());
    outputTracksFlag(track.has_collision(),
                     track.isPVContributor());
    outputTracksLabel(track.globalIndex());
  }

  // function to process trigger counters, accounting for BC selection bits
  void processCountersTrg(BCs const& bcs, aod::FT0s const&, aod::Zdcs const&)
  {
    const auto& firstBc = bcs.iteratorAt(0);
    if (runNumber != firstBc.runNumber())
      runNumber = firstBc.runNumber();

    auto hCountersTrg = getHist(TH1, "reco/hCountersTrg");
    auto hCountersTrgBcSel = getHist(TH1, "reco/hCountersTrgBcSel");
    auto hLumi = getHist(TH1, "reco/hLumi");
    auto hLumiBcSel = getHist(TH1, "reco/hLumiBcSel");

    // Cross sections in ub. Using dummy -1 if lumi estimator is not reliable
    float csTCE = 10.36e6;
    const float csZEM = 415.2e6; // see AN: https://alice-notes.web.cern.ch/node/1515
    const float csZNC = 214.5e6; // see AN: https://alice-notes.web.cern.ch/node/1515
    if (runNumber > 543437 && runNumber < 543514) {
      csTCE = 8.3e6;
    }
    if (runNumber >= 543514) {
      csTCE = 4.10e6; // see AN: https://alice-notes.web.cern.ch/node/1515
    }

    for (const auto& bc : bcs) {
      bool hasFT0 = bc.has_foundFT0();
      bool hasZDC = bc.has_foundZDC();
      if (!hasFT0 && !hasZDC)
        continue;
      bool isSelectedBc = true;
      if (rejectAtTFBoundary && !bc.selection_bit(aod::evsel::kNoTimeFrameBorder))
        isSelectedBc = false;
      if (noITSROFrameBorder && !bc.selection_bit(aod::evsel::kNoITSROFrameBorder))
        isSelectedBc = false;
      if (hasFT0) {
        auto ft0TrgMask = bc.ft0().triggerMask();
        if (TESTBIT(ft0TrgMask, o2::fit::Triggers::bitVertex)) {
          hCountersTrg->Fill("TVX", 1);
          if (isSelectedBc)
            hCountersTrgBcSel->Fill("TVX", 1);
        }
        if (TESTBIT(ft0TrgMask, o2::fit::Triggers::bitVertex) && TESTBIT(ft0TrgMask, o2::fit::Triggers::bitCen)) {
          hCountersTrg->Fill("TCE", 1);
          hLumi->Fill("TCE", 1. / csTCE);
          if (isSelectedBc) {
            hCountersTrgBcSel->Fill("TCE", 1);
            hLumiBcSel->Fill("TCE", 1. / csTCE);
          }
        }
      }
      if (hasZDC) {
        if (bc.selection_bit(aod::evsel::kIsBBZNA) || bc.selection_bit(aod::evsel::kIsBBZNC)) {
          hCountersTrg->Fill("ZEM", 1);
          hLumi->Fill("ZEM", 1. / csZEM);
          if (isSelectedBc) {
            hCountersTrgBcSel->Fill("ZEM", 1);
            hLumiBcSel->Fill("ZEM", 1. / csZEM);
          }
        }
        if (bc.selection_bit(aod::evsel::kIsBBZNC)) {
          hCountersTrg->Fill("ZNC", 1);
          hLumi->Fill("ZNC", 1. / csZNC);
          if (isSelectedBc) {
            hCountersTrgBcSel->Fill("ZNC", 1);
            hLumiBcSel->Fill("ZNC", 1. / csZNC);
          }
        }
      }
    }
  }

  PROCESS_SWITCH(SGCandProducer, processCountersTrg, "Produce trigger counters and luminosity histograms", true);

  // function to process reconstructed data
  template <typename TCol>
  void processReco(std::string histdir, TCol const& collision, BCs const& bcs,
                   TCs const& tracks, FWs const& fwdtracks,
                   aod::FV0As const& fv0as, aod::FT0s const& ft0s, aod::FDDs const& fdds)
  {
    if (verboseInfo)
      LOGF(debug, "<SGCandProducer>  collision %d", collision.globalIndex());
    getHist(TH1, histdir + "/Stat")->Fill(0., 1.);
    // reject collisions at TF boundaries
    if (rejectAtTFBoundary && !collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
      return;
    }
    getHist(TH1, histdir + "/Stat")->Fill(1., 1.);
    // reject collisions at ITS RO TF boundaries
    if (noITSROFrameBorder && !collision.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
      return;
    }
    getHist(TH1, histdir + "/Stat")->Fill(2., 1.);
    // reject Same Bunch PileUp
    if (noSameBunchPileUp && !collision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      return;
    }
    getHist(TH1, histdir + "/Stat")->Fill(3., 1.);
    // check vertex matching to FT0
    if (IsGoodVertex && !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return;
    }
    getHist(TH1, histdir + "/Stat")->Fill(4., 1.);
    // reject ITS Only vertices
    if (ITSTPCVertex && !collision.selection_bit(aod::evsel::kIsVertexITSTPC)) {
      return;
    }
    getHist(TH1, histdir + "/Stat")->Fill(5., 1.);
    // nominal BC
    if (!collision.has_foundBC()) {
      return;
    }
    getHist(TH1, histdir + "/Stat")->Fill(6., 1.);
    // RCT CBT for collision check
    if (isGoodRCTCollision && !myRCTChecker(collision)) {
      return;
    }
    getHist(TH1, histdir + "/Stat")->Fill(7., 1.);
    // RCT CBT+ZDC for collision check
    if (isGoodRCTZdc && !myRCTChecker(collision)) {
      return;
    }
    getHist(TH1, histdir + "/Stat")->Fill(8., 1.);

    //
    const int trs = collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard) ? 1 : 0;
    const int trofs = collision.selection_bit(o2::aod::evsel::kNoCollInRofStandard) ? 1 : 0;
    const int hmpr = collision.selection_bit(o2::aod::evsel::kNoHighMultCollInPrevRof) ? 1 : 0;
    const int tfb = collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder) ? 1 : 0;
    const int itsROFb = collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder) ? 1 : 0;
    const int sbp = collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup) ? 1 : 0;
    const int zVtxFT0vPv = collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV) ? 1 : 0;
    const int vtxITSTPC = collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC) ? 1 : 0;
    auto bc = collision.template foundBC_as<BCs>();
    double ir = 0.;
    const uint64_t ts = bc.timestamp();
    const int runnumber = bc.runNumber();
    if (bc.has_zdc()) {
      ir = mRateFetcher.fetch(ccdb.service, ts, runnumber, "ZNC hadronic") * 1.e-3;
    }
    auto newbc = bc;

    // obtain slice of compatible BCs
    auto bcRange = udhelpers::compatibleBCs(collision, sameCuts.NDtcoll(), bcs, sameCuts.minNBCs());
    auto isSGEvent = sgSelector.IsSelected(sameCuts, collision, bcRange, bc);
    // auto isSGEvent = sgSelector.IsSelected(sameCuts, collision, bcRange, tracks);
    int issgevent = isSGEvent.value;
    if (isSGEvent.bc && issgevent < 2) {
      newbc = *(isSGEvent.bc);
    } else {
      if (verboseInfo)
        LOGF(info, "No Newbc %i", bc.globalBC());
    }
    getHist(TH1, histdir + "/Stat")->Fill(issgevent + 10, 1.);
    if ((storeDG && issgevent == o2::aod::sgselector::DoubleGap) || (storeSG && (issgevent == o2::aod::sgselector::SingleGapA || issgevent == o2::aod::sgselector::SingleGapC))) {
      if (verboseInfo)
        LOGF(info, "Current BC: %i, %i, %i", bc.globalBC(), newbc.globalBC(), issgevent);
      if (sameCuts.minRgtrwTOF()) {
        if (udhelpers::rPVtrwTOF<true>(tracks, collision.numContrib()) < sameCuts.minRgtrwTOF())
          return;
      }
      upchelpers::FITInfo fitInfo{};
      const uint8_t chFT0A = 0;
      const uint8_t chFT0C = 0;
      const uint8_t chFDDA = 0;
      const uint8_t chFDDC = 0;
      const uint8_t chFV0A = 0;
      const int occ = collision.trackOccupancyInTimeRange();
      udhelpers::getFITinfo(fitInfo, newbc, bcs, ft0s, fv0as, fdds);
      const int upc_flag = (collision.flags() & dataformats::Vertex<o2::dataformats::TimeStamp<int>>::Flags::UPCMode) ? 1 : 0;
      // update SG candidates tables
      outputCollisions(bc.globalBC(), bc.runNumber(),
                       collision.posX(), collision.posY(), collision.posZ(), upc_flag,
                       collision.numContrib(), udhelpers::netCharge<true>(tracks),
                       1.); // rtrwTOF); //omit the calculation to speed up the things while skimming

      outputSGCollisions(issgevent);
      outputCollisionsSels(fitInfo.ampFT0A, fitInfo.ampFT0C, fitInfo.timeFT0A, fitInfo.timeFT0C,
                           fitInfo.triggerMaskFT0,
                           fitInfo.ampFDDA, fitInfo.ampFDDC, fitInfo.timeFDDA, fitInfo.timeFDDC,
                           fitInfo.triggerMaskFDD,
                           fitInfo.ampFV0A, fitInfo.timeFV0A, fitInfo.triggerMaskFV0A,
                           fitInfo.BBFT0Apf, fitInfo.BBFT0Cpf, fitInfo.BGFT0Apf, fitInfo.BGFT0Cpf,
                           fitInfo.BBFV0Apf, fitInfo.BGFV0Apf,
                           fitInfo.BBFDDApf, fitInfo.BBFDDCpf, fitInfo.BGFDDApf, fitInfo.BGFDDCpf);
      outputCollisionSelExtras(chFT0A, chFT0C, chFDDA, chFDDC, chFV0A, occ, ir, trs, trofs, hmpr, tfb, itsROFb, sbp, zVtxFT0vPv, vtxITSTPC, collision.rct_raw());
      outputCollsLabels(collision.globalIndex());
      if (newbc.has_zdc()) {
        auto zdc = newbc.zdc();
        udZdcsReduced(outputCollisions.lastIndex(), zdc.timeZNA(), zdc.timeZNC(), zdc.energyCommonZNA(), zdc.energyCommonZNC());
      } else {
        udZdcsReduced(outputCollisions.lastIndex(), -999, -999, -999, -999);
      }
      // update SGTracks tables
      if (fillTrackTables) {
        for (const auto& track : tracks) {
          if (track.pt() > sameCuts.minPt() && track.eta() > sameCuts.minEta() && track.eta() < sameCuts.maxEta()) {
            if (track.isPVContributor()) {
              updateUDTrackTables(outputCollisions.lastIndex(), track, bc.globalBC());
            } else if (saveAllTracks) {
              if (track.itsClusterSizes() && track.itsChi2NCl() > 0 && ((track.tpcNClsFindable() == 0 && savenonPVCITSOnlyTracks) || track.tpcNClsFindable() > 50))
                updateUDTrackTables(outputCollisions.lastIndex(), track, bc.globalBC());
              // if (track.isPVContributor())  updateUDTrackTables(outputCollisions.lastIndex(), track, bc.globalBC());
            }
          }
        }
      }
      // update SGFwdTracks tables
      if (fillFwdTrackTables) {
        if (sameCuts.withFwdTracks()) {
          for (const auto& fwdtrack : fwdtracks) {
            if (!sgSelector.FwdTrkSelector(fwdtrack))
              updateUDFwdTrackTables(fwdtrack, bc.globalBC());
          }
        }
      }
    }
  }

  void init(InitContext& context)
  {
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setFatalWhenNull(false);
    sameCuts = (SGCutParHolder)SGCuts;

    // add histograms for the different process functions
    histPointers.clear();
    if (context.mOptions.get<bool>("processData")) {
      histPointers.insert({"reco/Stat", registry.add("reco/Stat", "Cut statistics; Selection criterion; Collisions", {HistType::kTH1F, {{14, -0.5, 13.5}}})});

      const AxisSpec axisCountersTrg{10, 0.5, 10.5, ""};
      histPointers.insert({"reco/hCountersTrg", registry.add("reco/hCountersTrg", "Trigger counts before selections; Trigger; Counts", {HistType::kTH1F, {axisCountersTrg}})});
      histPointers.insert({"reco/hCountersTrgBcSel", registry.add("reco/hCountersTrgSel", "Trigger counts after BC selections; Trigger; Counts", {HistType::kTH1F, {axisCountersTrg}})});
      histPointers.insert({"reco/hLumi", registry.add("reco/hLumi", "Integrated luminosity before selections; Trigger; Luminosity, 1/#mub", {HistType::kTH1F, {axisCountersTrg}})});
      histPointers.insert({"reco/hLumiBcSel", registry.add("reco/hLumiBcSel", "Integrated luminosity before selections; Trigger; Luminosity, 1/#mub", {HistType::kTH1F, {axisCountersTrg}})});
      auto hCountersTrg = getHist(TH1, "reco/hCountersTrg");
      auto hCountersTrgBcSel = getHist(TH1, "reco/hCountersTrgBcSel");
      auto hLumi = getHist(TH1, "reco/hLumi");
      auto hLumiBcSel = getHist(TH1, "reco/hLumiBcSel");
      for (const auto& h : {hCountersTrg, hCountersTrgBcSel, hLumi, hLumiBcSel}) {
        h->GetXaxis()->SetBinLabel(1, "TVX");
        h->GetXaxis()->SetBinLabel(2, "TCE");
        h->GetXaxis()->SetBinLabel(3, "ZEM");
        h->GetXaxis()->SetBinLabel(4, "ZNC");
      }
    }
    if (context.mOptions.get<bool>("processMcData")) {
      histPointers.insert({"MCreco/Stat", registry.add("MCreco/Stat", "Cut statistics; Selection criterion; Collisions", {HistType::kTH1F, {{14, -0.5, 13.5}}})});
    }

    if (isGoodRCTZdc) {
      myRCTChecker.init("CBT", true);
    }
  }

  // process function for reconstructed data
  void processData(CC const& collision, BCs const& bcs, TCs const& tracks, FWs const& fwdtracks,
                   aod::Zdcs const& /*zdcs*/, aod::FV0As const& fv0as, aod::FT0s const& ft0s, aod::FDDs const& fdds)
  {
    processReco(std::string("reco"), collision, bcs, tracks, fwdtracks, fv0as, ft0s, fdds);
  }
  PROCESS_SWITCH(SGCandProducer, processData, "Produce UD table with data", true);

  // process function for reconstructed MC data
  void processMcData(MCCC const& collision, aod::McCollisions const& /*mccollisions*/, BCs const& bcs,
                     TCs const& tracks, FWs const& fwdtracks, aod::Zdcs const& /*zdcs*/, aod::FV0As const& fv0as,
                     aod::FT0s const& ft0s, aod::FDDs const& fdds)
  {
    // select specific processes with the GeneratorID
    if (!collision.has_mcCollision())
      return;
    auto mccol = collision.mcCollision();
    if (verboseInfo)
      LOGF(info, "GeneratorId %d (%d)", mccol.getGeneratorId(), generatorIds->size());

    if (std::find(generatorIds->begin(), generatorIds->end(), mccol.getGeneratorId()) != generatorIds->end()) {
      if (verboseInfo)
        LOGF(info, "Event with good generatorId");
      processReco(std::string("MCreco"), collision, bcs, tracks, fwdtracks, fv0as, ft0s, fdds);
    }
  }
  PROCESS_SWITCH(SGCandProducer, processMcData, "Produce UD tables with MC data", false);
};

struct McSGCandProducer {
  // MC tables
  Configurable<bool> verboseInfoMC{"verboseInfoMC", false, "Print general info to terminal; default it false."};
  Produces<aod::UDMcCollisions> outputMcCollisions;
  Produces<aod::UDMcParticles> outputMcParticles;
  Produces<aod::UDMcCollsLabels> outputMcCollsLabels;
  Produces<aod::UDMcTrackLabels> outputMcTrackLabels;

  // save all McTruth, even if the collisions is not reconstructed
  Configurable<std::vector<int>> generatorIds{"generatorIds", std::vector<int>{-1}, "MC generatorIds to process"};
  Configurable<bool> saveAllMcCollisions{"saveAllMcCollisions", true, "save all McCollisions"};

  using CCs = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>;
  using BCs = soa::Join<aod::BCsWithTimestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;
  using TCs = soa::Join<aod::Tracks, aod::McTrackLabels>;
  using UDCCs = soa::Join<aod::UDCollisions, aod::UDCollsLabels, aod::SGCollisions>;
  using UDTCs = soa::Join<aod::UDTracks, aod::UDTracksLabels>;

  // prepare slices
  SliceCache cache;
  PresliceUnsorted<aod::McParticles> mcPartsPerMcCollision = aod::mcparticle::mcCollisionId;
  Preslice<UDTCs> udtracksPerUDCollision = aod::udtrack::udCollisionId;

  // initialize histogram registry
  HistogramRegistry registry{
    "registry",
    {}};

  template <typename TMcCollision>
  void updateUDMcCollisions(TMcCollision const& mccol, uint64_t globBC)
  {
    // save mccol
    outputMcCollisions(globBC,
                       mccol.generatorsID(),
                       mccol.posX(),
                       mccol.posY(),
                       mccol.posZ(),
                       mccol.t(),
                       mccol.weight(),
                       mccol.impactParameter());
  }

  template <typename TMcParticle>
  void updateUDMcParticle(TMcParticle const& McPart, int64_t McCollisionId, std::map<int64_t, int64_t>& mcPartIsSaved)
  {
    // save McPart
    // mother and daughter indices are set to -1
    // ATTENTION: this can be improved to also include mother and daughter indices
    std::vector<int32_t> newmids;
    int32_t newdids[2] = {-1, -1};

    // update UDMcParticles
    if (mcPartIsSaved.find(McPart.globalIndex()) == mcPartIsSaved.end()) {
      outputMcParticles(McCollisionId,
                        McPart.pdgCode(),
                        McPart.statusCode(),
                        McPart.flags(),
                        newmids,
                        newdids,
                        McPart.weight(),
                        McPart.px(),
                        McPart.py(),
                        McPart.pz(),
                        McPart.e());
      mcPartIsSaved[McPart.globalIndex()] = outputMcParticles.lastIndex();
    }
  }

  template <typename TMcParticles>
  void updateUDMcParticles(TMcParticles const& McParts, int64_t McCollisionId, std::map<int64_t, int64_t>& mcPartIsSaved)
  {
    // save McParts
    // new mother and daughter ids
    std::vector<int32_t> newmids;
    int32_t newdids[2] = {-1, -1};
    int64_t newval = -1;

    // Determine the particle indices within the UDMcParticles table
    // before filling the table
    // This is needed to be able to assign the new daughter indices
    std::map<int64_t, int64_t> oldnew;
    auto lastId = outputMcParticles.lastIndex();
    for (const auto& mcpart : McParts) {
      auto oldId = mcpart.globalIndex();
      if (mcPartIsSaved.find(oldId) != mcPartIsSaved.end()) {
        oldnew[oldId] = mcPartIsSaved[oldId];
      } else {
        lastId++;
        oldnew[oldId] = lastId;
      }
    }

    // all particles of the McCollision are saved
    for (const auto& mcpart : McParts) {
      if (mcPartIsSaved.find(mcpart.globalIndex()) == mcPartIsSaved.end()) {
        // mothers
        newmids.clear();
        auto oldmids = mcpart.mothersIds();
        for (const auto& oldmid : oldmids) {
          auto m = McParts.rawIteratorAt(oldmid);
          if (verboseInfoMC)
            LOGF(debug, "    m %d", m.globalIndex());
          if (mcPartIsSaved.find(oldmid) != mcPartIsSaved.end()) {
            newval = mcPartIsSaved[oldmid];
          } else {
            newval = -1;
          }
          newmids.push_back(newval);
        }
        // daughters
        auto olddids = mcpart.daughtersIds();
        for (uint ii = 0; ii < olddids.size(); ii++) {
          if (oldnew.find(olddids[ii]) != oldnew.end()) {
            newval = oldnew[olddids[ii]];
          } else {
            newval = -1;
          }
          newdids[ii] = newval;
        }
        if (verboseInfoMC)
          LOGF(debug, " ms %i ds %i", oldmids.size(), olddids.size());

        // update UDMcParticles
        outputMcParticles(McCollisionId,
                          mcpart.pdgCode(),
                          mcpart.statusCode(),
                          mcpart.flags(),
                          newmids,
                          newdids,
                          mcpart.weight(),
                          mcpart.px(),
                          mcpart.py(),
                          mcpart.pz(),
                          mcpart.e());
        mcPartIsSaved[mcpart.globalIndex()] = outputMcParticles.lastIndex();
      }
    }
  }

  template <typename TTrack>
  void updateUDMcTrackLabel(TTrack const& udtrack, std::map<int64_t, int64_t>& mcPartIsSaved)
  {
    // udtrack (UDTCs) -> track (TCs) -> mcTrack (McParticles) -> udMcTrack (UDMcParticles)
    auto trackId = udtrack.trackId();
    if (trackId >= 0) {
      auto track = udtrack.template track_as<TCs>();
      auto mcTrackId = track.mcParticleId();
      if (mcTrackId >= 0) {
        if (mcPartIsSaved.find(mcTrackId) != mcPartIsSaved.end()) {
          outputMcTrackLabels(mcPartIsSaved[mcTrackId], track.mcMask());
        } else {
          outputMcTrackLabels(-1, track.mcMask());
        }
      } else {
        outputMcTrackLabels(-1, track.mcMask());
      }
    } else {
      outputMcTrackLabels(-1, -1);
    }
  }

  template <typename TTrack>
  void updateUDMcTrackLabels(TTrack const& udtracks, std::map<int64_t, int64_t>& mcPartIsSaved)
  {
    // loop over all tracks
    for (const auto& udtrack : udtracks) {
      // udtrack (UDTCs) -> track (TCs) -> mcTrack (McParticles) -> udMcTrack (UDMcParticles)
      auto trackId = udtrack.trackId();
      if (trackId >= 0) {
        auto track = udtrack.template track_as<TCs>();
        auto mcTrackId = track.mcParticleId();
        if (mcTrackId >= 0) {
          if (mcPartIsSaved.find(mcTrackId) != mcPartIsSaved.end()) {
            outputMcTrackLabels(mcPartIsSaved[mcTrackId], track.mcMask());
          } else {
            outputMcTrackLabels(-1, track.mcMask());
          }
        } else {
          outputMcTrackLabels(-1, track.mcMask());
        }
      } else {
        outputMcTrackLabels(-1, -1);
      }
    }
  }

  // updating McTruth data and links to reconstructed data
  void procWithSgCand(aod::McCollisions const& mccols, aod::McParticles const& mcparts,
                      UDCCs const& sgcands, UDTCs const& udtracks)
  {
    // use a hash table to keep track of the McCollisions which have been added to the UDMcCollision table
    // {McCollisionId : udMcCollisionId}
    // similar for the McParticles which have been added to the UDMcParticle table
    // {McParticleId : udMcParticleId}
    std::map<int64_t, int64_t> mcColIsSaved;
    std::map<int64_t, int64_t> mcPartIsSaved;

    // loop over McCollisions and UDCCs simultaneously
    auto mccol = mccols.iteratorAt(0);
    auto mcOfInterest = std::find(generatorIds->begin(), generatorIds->end(), mccol.getGeneratorId()) != generatorIds->end();
    auto lastmccol = mccols.iteratorAt(mccols.size() - 1);
    auto mccolAtEnd = false;

    auto sgcand = sgcands.iteratorAt(0);
    auto lastsgcand = sgcands.iteratorAt(sgcands.size() - 1);
    auto sgcandAtEnd = false;

    // advance dgcand and mccol until both are AtEnd
    int64_t mccolId = mccol.globalIndex();
    int64_t mcsgId = -1;
    bool goon = true;
    while (goon) {
      auto globBC = mccol.bc_as<BCs>().globalBC();
      // check if dgcand has an associated McCollision
      if (sgcand.has_collision()) {
        auto sgcandCol = sgcand.collision_as<CCs>();
        // colId = sgcandCol.globalIndex();
        if (sgcandCol.has_mcCollision()) {
          mcsgId = sgcandCol.mcCollision().globalIndex();
        } else {
          mcsgId = -1;
        }
      } else {
        mcsgId = -1;
      }
      if (verboseInfoMC)
        LOGF(info, "\nStart of loop mcsgId %d mccolId %d", mcsgId, mccolId);

      // two cases to consider
      // 1. mcdgId <= mccolId: the event to process is a dgcand. In this case the Mc tables as well as the McLabel tables are updated
      // 2. mccolId < mcdgId: the event to process is an MC event of interest without reconstructed dgcand. In this case only the Mc tables are updated
      if ((!sgcandAtEnd && !mccolAtEnd && (mcsgId <= mccolId)) || mccolAtEnd) {
        // this is case 1.
        if (verboseInfoMC)
          LOGF(info, "Doing case 1 with mcsgId %d", mcsgId);

        // update UDMcCollisions and UDMcColsLabels (for each UDCollision -> UDMcCollisions)
        // update UDMcParticles and UDMcTrackLabels (for each UDTrack -> UDMcParticles)
        // get dgcand tracks
        auto sgTracks = udtracks.sliceByCached(aod::udtrack::udCollisionId, sgcand.globalIndex(), cache);

        // If the sgcand has an associated McCollision then the McCollision and all associated
        // McParticles are saved
        // but only consider generated events of interest
        if (mcsgId >= 0 && mcOfInterest) {
          if (mcColIsSaved.find(mcsgId) == mcColIsSaved.end()) {
            if (verboseInfoMC)
              LOGF(info, "  Saving McCollision %d", mcsgId);
            // update UDMcCollisions
            auto sgcandMcCol = sgcand.collision_as<CCs>().mcCollision();
            updateUDMcCollisions(sgcandMcCol, globBC);
            mcColIsSaved[mcsgId] = outputMcCollisions.lastIndex();
          }

          // update UDMcColsLabels (for each UDCollision -> UDMcCollisions)
          outputMcCollsLabels(mcColIsSaved[mcsgId]);

          // update UDMcParticles
          auto mcPartsSlice = mcparts.sliceBy(mcPartsPerMcCollision, mcsgId);
          updateUDMcParticles(mcPartsSlice, mcColIsSaved[mcsgId], mcPartIsSaved);

          // update UDMcTrackLabels (for each UDTrack -> UDMcParticles)
          updateUDMcTrackLabels(sgTracks, mcPartIsSaved);

        } else {
          // If the sgcand has no associated McCollision then only the McParticles which are associated
          // with the tracks of the sgcand are saved
          if (verboseInfoMC)
            LOGF(info, "  Saving McCollision %d", -1);

          // update UDMcColsLabels (for each UDCollision -> UDMcCollisions)
          outputMcCollsLabels(-1);

          // update UDMcParticles and UDMcTrackLabels (for each UDTrack -> UDMcParticles)
          // loop over tracks of dgcand
          for (const auto& sgtrack : sgTracks) {
            if (sgtrack.has_track()) {
              auto track = sgtrack.track_as<TCs>();
              if (track.has_mcParticle()) {
                auto mcPart = track.mcParticle();
                auto mcCol = mcPart.mcCollision();
                if (mcColIsSaved.find(mcCol.globalIndex()) == mcColIsSaved.end()) {
                  updateUDMcCollisions(mcCol, globBC);
                  mcColIsSaved[mcCol.globalIndex()] = outputMcCollisions.lastIndex();
                }
                updateUDMcParticle(mcPart, mcColIsSaved[mcCol.globalIndex()], mcPartIsSaved);
                updateUDMcTrackLabel(sgtrack, mcPartIsSaved);
              } else {
                outputMcTrackLabels(-1, track.mcMask());
              }
            } else {
              outputMcTrackLabels(-1, -1);
            }
          }
        }
        // advance sgcand
        if (sgcand != lastsgcand) {
          sgcand++;
        } else {
          sgcandAtEnd = true;
        }
      } else {
        // this is case 2.
        if (verboseInfoMC)
          LOGF(info, "Doing case 2");

        // update UDMcCollisions and UDMcParticles
        // but only consider generated events of interest
        if (mcOfInterest && mcColIsSaved.find(mccolId) == mcColIsSaved.end()) {
          if (verboseInfoMC)
            LOGF(info, "  Saving McCollision %d", mccolId);
          // update UDMcCollisions
          updateUDMcCollisions(mccol, globBC);
          mcColIsSaved[mccolId] = outputMcCollisions.lastIndex();

          // update UDMcParticles
          auto mcPartsSlice = mcparts.sliceBy(mcPartsPerMcCollision, mccolId);
          updateUDMcParticles(mcPartsSlice, mcColIsSaved[mccolId], mcPartIsSaved);
        }

        // advance mccol
        if (mccol != lastmccol) {
          mccol++;
          mcOfInterest = std::find(generatorIds->begin(), generatorIds->end(), mccol.getGeneratorId()) != generatorIds->end();
          mccolId = mccol.globalIndex();
        } else {
          mccolAtEnd = true;
        }
      }

      goon = !sgcandAtEnd || !mccolAtEnd;
      if (verboseInfoMC)
        LOGF(info, "End of loop mcsgId %d mccolId %d", mcsgId, mccolId);
    }
  }

  // updating McTruth data only
  void procWithoutSgCand(aod::McCollisions const& mccols, aod::McParticles const& mcparts)
  {
    // use a hash table to keep track of the McCollisions which have been added to the UDMcCollision table
    // {McCollisionId : udMcCollisionId}
    // similar for the McParticles which have been added to the UDMcParticle table
    // {McParticleId : udMcParticleId}
    std::map<int64_t, int64_t> mcColIsSaved;
    std::map<int64_t, int64_t> mcPartIsSaved;

    // loop over McCollisions
    for (auto const& mccol : mccols) {
      int64_t mccolId = mccol.globalIndex();
      uint64_t globBC = mccol.bc_as<BCs>().globalBC();

      // update UDMcCollisions and UDMcParticles
      if (mcColIsSaved.find(mccolId) == mcColIsSaved.end()) {
        if (verboseInfoMC)
          LOGF(info, "  Saving McCollision %d", mccolId);

        // update UDMcCollisions
        updateUDMcCollisions(mccol, globBC);
        mcColIsSaved[mccolId] = outputMcCollisions.lastIndex();

        // update UDMcParticles
        auto mcPartsSlice = mcparts.sliceBy(mcPartsPerMcCollision, mccolId);
        updateUDMcParticles(mcPartsSlice, mcColIsSaved[mccolId], mcPartIsSaved);
      }
    }
  }

  void init(InitContext& context)
  {
    // add histograms for the different process functions
    if (context.mOptions.get<bool>("processMC")) {
      registry.add("mcTruth/collisions", "Number of associated collisions", {HistType::kTH1F, {{11, -0.5, 10.5}}});
      registry.add("mcTruth/collType", "Collision type", {HistType::kTH1F, {{5, -0.5, 4.5}}});
      registry.add("mcTruth/IVMpt", "Invariant mass versus p_{T}", {HistType::kTH2F, {{150, 0.0, 3.0}, {150, 0.0, 3.0}}});
    }
  }

  // process function for MC data
  // save the MC truth of all events of interest and of the DG events
  void processMC(aod::McCollisions const& mccols, aod::McParticles const& mcparts,
                 UDCCs const& sgcands, UDTCs const& udtracks,
                 CCs const& /*collisions*/, BCs const& /*bcs*/, TCs const& /*tracks*/)
  {
    if (verboseInfoMC) {
      LOGF(info, "Number of McCollisions %d", mccols.size());
      LOGF(info, "Number of SG candidates %d", sgcands.size());
      LOGF(info, "Number of UD tracks %d", udtracks.size());
    }
    if (mccols.size() > 0) {
      if (sgcands.size() > 0) {
        procWithSgCand(mccols, mcparts, sgcands, udtracks);
      } else {
        if (saveAllMcCollisions) {
          procWithoutSgCand(mccols, mcparts);
        }
      }
    }
  }
  PROCESS_SWITCH(McSGCandProducer, processMC, "Produce MC tables", false);

  void processDummy(aod::Collisions const& /*collisions*/)
  {
    // do nothing
    if (verboseInfoMC)
      LOGF(info, "Running dummy process function!");
  }
  PROCESS_SWITCH(McSGCandProducer, processDummy, "Dummy function", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<SGCandProducer>(cfgc, TaskName{"sgcandproducer"}),
    adaptAnalysisTask<McSGCandProducer>(cfgc, TaskName{"mcsgcandproducer"})};
  return workflow;
}
