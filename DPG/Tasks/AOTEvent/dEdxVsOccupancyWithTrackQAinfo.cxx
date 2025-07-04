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

/// \file dEdxVsOccupancyWithTrackQAinfo.cxx
/// \brief dE/dx vs occupancy QA task with more detailed checks
///
/// \author Igor Altsybeev <Igor.Altsybeev@cern.ch>

#include <vector>
#include <map>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/ctpRateFetcher.h"
#include "CCDB/BasicCCDBManager.h"
#include "Framework/HistogramRegistry.h"
#include "CommonDataFormat/BunchFilling.h"
#include "DataFormatsParameters/GRPLHCIFData.h"
#include "DataFormatsParameters/GRPECSObject.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "DataFormatsParameters/AggregatedRunInfo.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TH3.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::aod::evsel;

using BCsRun3 = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;
using ColEvSels = soa::Join<aod::Collisions, aod::EvSels, aod::Mults>;
// using ColEvSels = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Cs>;
using FullTracksIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TrackSelection, aod::TrackSelectionExtension, aod::TracksDCA>;

struct dEdxVsOccupancyWithTrackQAinfoTask {
  // configurables for study of occupancy in time windows
  // Configurable<bool> confAddBasicQAhistos{"AddBasicQAhistos", true, "0 - add basic histograms, 1 - skip"}; // o2-linter: disable=name/configurable (temporary fix)
  // Configurable<float> confTimeIntervalForOccupancyCalculation{"TimeIntervalForOccupancyCalculation", 100, "Time interval for TPC occupancy calculation, us"}; // o2-linter: disable=name/configurable (temporary fix)
  // Configurable<bool> confFlagCentralityIsAvailable{"FlagCentralityIsAvailable", true, "Fill centrality-related historams"};                                   // o2-linter: disable=name/configurable (temporary fix)
  // Configurable<bool> confFlagManyHeavyHistos{"FlagManyHeavyHistos", true, "Fill more TH2, TH3, THn historams"};                                               // o2-linter: disable=name/configurable (temporary fix)

  // event and track cuts for given event
  Configurable<float> confCutVertZMinThisEvent{"VzMinThisEvent", -10, "vZ cut for a current event"};                           // o2-linter: disable=name/configurable (temporary fix)
  Configurable<float> confCutVertZMaxThisEvent{"VzMaxThisEvent", 10, "vZ cut for a current event"};                            // o2-linter: disable=name/configurable (temporary fix)
  Configurable<float> confCutPtMinThisEvent{"PtMinThisEvent", 0.2, "pt cut for particles in a current event"};                 // o2-linter: disable=name/configurable (temporary fix)
  Configurable<float> confCutPtMaxThisEvent{"PtMaxThisEvent", 100., "pt cut for particles in a current event"};                // o2-linter: disable=name/configurable (temporary fix)
  Configurable<float> confCutEtaMinTracksThisEvent{"EtaMinTracksThisEvent", -0.8, "eta cut for particles in a current event"}; // o2-linter: disable=name/configurable (temporary fix)
  Configurable<float> confCutEtaMaxTracksThisEvent{"EtaMaxTracksThisEvent", 0.8, "eta cut for particles in a current event"};  // o2-linter: disable=name/configurable (temporary fix)
  Configurable<int> confCutMinTPCcls{"MinNumTPCcls", 70, "min number of TPC clusters for a current event"};                    // o2-linter: disable=name/configurable (temporary fix)

  uint64_t minGlobalBC = 0;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  bool* applySelection = NULL;
  int nBCsPerOrbit = 3564;
  int lastRunNumber = -1;
  int nOrbits;
  double minOrbit;
  int64_t bcSOR = 0;                      // global bc of the start of the first orbit, setting 0 by default for unanchored MC
  int64_t nBCsPerTF = 128 * nBCsPerOrbit; // duration of TF in bcs, should be 128*3564 or 32*3564, setting 128 orbits by default sfor unanchored MC
  ctpRateFetcher mRateFetcher;

  // save time "slices" for several collisions for QA
  bool flagFillQAtimeOccupHist = false;
  int nCollisionsForTimeBinQA = 40;
  int counterQAtimeOccupHistos = 0;

  void init(InitContext&)
  {
    // ccdb->setURL("http://ccdb-test.cern.ch:8080");
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    // dE/dx
    AxisSpec axisDeDx{800, 0.0, 800.0, "dE/dx (a. u.)"};
    histos.add("dEdx_vs_Momentum_CORRECTED", "dE/dx", kTH2F, {{1000, -5.0, 5.0, "#it{p}/Z (GeV/c)"}, axisDeDx});
    histos.add("dEdx_vs_Momentum", "dE/dx", kTH2F, {{1000, -5.0, 5.0, "#it{p}/Z (GeV/c)"}, axisDeDx});
    // histos.add("dEdx_vs_Momentum_occupBelow200", "dE/dx", kTH2F, {{1000, -5.0, 5.0, "#it{p}/Z (GeV/c)"}, axisDeDx});
    // histos.add("dEdx_vs_Momentum_occupBelow200_kNoCollStd", "dE/dx", kTH2F, {{1000, -5.0, 5.0, "#it{p}/Z (GeV/c)"}, axisDeDx});
    // histos.add("dEdx_vs_Momentum_occupAbove4000", "dE/dx", kTH2F, {{1000, -5.0, 5.0, "#it{p}/Z (GeV/c)"}, axisDeDx});
    // histos.add("dEdx_vs_Momentum_NegativeFractionNclsPID", "dE/dx", kTH2F, {{1000, -5.0, 5.0, "#it{p}/Z (GeV/c)"}, axisDeDx});
    // histos.add("dEdx_vs_Momentum_HighFractionNclsNonPID", "dE/dx", kTH2F, {{1000, -5.0, 5.0, "#it{p}/Z (GeV/c)"}, axisDeDx});
    AxisSpec axisBinsOccupStudydEdx{{0., 500, 1000, 2000, 4000, 6000, 8000, 15000}, "p_{T}"};
    // histos.add("dEdx_vs_Momentum_vs_occup", "dE/dx", kTH3F, {{1000, -5.0, 5.0, "#it{p}/Z (GeV/c)"}, axisDeDx, axisBinsOccupStudydEdx});
    // if (confFlagManyHeavyHistos) {
    // histos.add("dEdx_vs_Momentum_vs_occup_eta_02_04", "dE/dx", kTH3F, {{1000, -5.0, 5.0, "#it{p}/Z (GeV/c)"}, axisDeDx, axisBinsOccupStudydEdx});
    // histos.add("dEdx_vs_Momentum_vs_occup_eta_04_02", "dE/dx", kTH3F, {{1000, -5.0, 5.0, "#it{p}/Z (GeV/c)"}, axisDeDx, axisBinsOccupStudydEdx});
    // }
    histos.add("dEdx_3OROC_tot_vs_Momentum_vs_occup_eta_02_04", "dE/dx", kTH3F, {{1000, -5.0, 5.0, "#it{p}/Z (GeV/c)"}, axisDeDx, axisBinsOccupStudydEdx});
    histos.add("dEdx_3OROC_tot_vs_Momentum_vs_occup_eta_04_02", "dE/dx", kTH3F, {{1000, -5.0, 5.0, "#it{p}/Z (GeV/c)"}, axisDeDx, axisBinsOccupStudydEdx});
    histos.add("dEdx_3OROC_max_vs_Momentum_vs_occup_eta_02_04", "dE/dx", kTH3F, {{1000, -5.0, 5.0, "#it{p}/Z (GeV/c)"}, axisDeDx, axisBinsOccupStudydEdx});
    histos.add("dEdx_3OROC_max_vs_Momentum_vs_occup_eta_04_02", "dE/dx", kTH3F, {{1000, -5.0, 5.0, "#it{p}/Z (GeV/c)"}, axisDeDx, axisBinsOccupStudydEdx});

    // track QA info
    histos.add("tpcdEdxMax0R_vs_Momentum", "dE/dx", kTH2F, {{1000, -5.0, 5.0, "#it{p}/Z (GeV/c)"}, axisDeDx});
    histos.add("tpcdEdxMax1R_vs_Momentum", "dE/dx", kTH2F, {{1000, -5.0, 5.0, "#it{p}/Z (GeV/c)"}, axisDeDx});
    histos.add("tpcdEdxMax2R_vs_Momentum", "dE/dx", kTH2F, {{1000, -5.0, 5.0, "#it{p}/Z (GeV/c)"}, axisDeDx});
    histos.add("tpcdEdxMax3R_vs_Momentum", "dE/dx", kTH2F, {{1000, -5.0, 5.0, "#it{p}/Z (GeV/c)"}, axisDeDx});

    histos.add("tpcdEdxTot0R_vs_Momentum", "dE/dx", kTH2F, {{1000, -5.0, 5.0, "#it{p}/Z (GeV/c)"}, axisDeDx});
    histos.add("tpcdEdxTot1R_vs_Momentum", "dE/dx", kTH2F, {{1000, -5.0, 5.0, "#it{p}/Z (GeV/c)"}, axisDeDx});
    histos.add("tpcdEdxTot2R_vs_Momentum", "dE/dx", kTH2F, {{1000, -5.0, 5.0, "#it{p}/Z (GeV/c)"}, axisDeDx});
    histos.add("tpcdEdxTot3R_vs_Momentum", "dE/dx", kTH2F, {{1000, -5.0, 5.0, "#it{p}/Z (GeV/c)"}, axisDeDx});

    histos.add("tpcdEdxTot3R_vs_dEdxFromTracks", "dE/dx", kTH2F, {axisDeDx, axisDeDx});
    histos.add("tpcdEdxTotSUM_vs_dEdxFromTracks", "dE/dx", kTH2F, {axisDeDx, axisDeDx});
    histos.add("tpcdEdxMaxSUM_vs_dEdxFromTracks", "dE/dx", kTH2F, {axisDeDx, axisDeDx});

    histos.add("tpcdEdxAverageMax_vs_dEdxFromTracks", "dE/dx", kTH2F, {axisDeDx, axisDeDx});
    histos.add("tpcdEdxAverageTot_vs_dEdxFromTracks", "dE/dx", kTH2F, {axisDeDx, axisDeDx});

    histos.add("tpcdEdxAverageMax_3OROC_vs_dEdxFromTracks", "dE/dx", kTH2F, {axisDeDx, axisDeDx});
    histos.add("tpcdEdxAverageTot_3OROC_vs_dEdxFromTracks", "dE/dx", kTH2F, {axisDeDx, axisDeDx});

    histos.add("tpcdEdxCORRECTED_vs_dEdxFromTracks", "dE/dx", kTH2F, {axisDeDx, axisDeDx});

    const AxisSpec axisDcaZ{1000, -5., 5., "DCA_{z}, cm"};
    histos.add("dcaXY_vs_dcaXYqa", "dE/dx", kTH2F, {axisDcaZ, {601, -300.5, 300.5, "DCA_{z}, cm"}});
    histos.add("dcaZ_vs_dcaZqa", "dE/dx", kTH2F, {axisDcaZ, {601, -300.5, 300.5, "DCA_{z}, cm"}});

    AxisSpec axisOccupancyForDeDxStudies{60, 0, 15000, "occupancy"};
    histos.add("dEdx_vs_centr_vs_occup_narrow_p_win", "dE/dx", kTH3F, {{20, 0, 4000, "nITStrk cls567"}, axisOccupancyForDeDxStudies, axisDeDx});
    // histos.add("dEdx_vs_centr_vs_occup_narrow_p_win_pos", "dE/dx", kTH3F, {{20, 0, 4000, "nITStrk cls567"}, axisOccupancyForDeDxStudies, axisDeDx});
    // histos.add("dEdx_vs_centr_vs_occup_narrow_p_win_neg", "dE/dx", kTH3F, {{20, 0, 4000, "nITStrk cls567"}, axisOccupancyForDeDxStudies, axisDeDx});
    // histos.add("dEdx_vs_centr_vs_occup_narrow_p_win_pos_FractionPIDclsInRange", "dE/dx", kTH3F, {{20, 0, 4000, "nITStrk cls567"}, axisOccupancyForDeDxStudies, axisDeDx});
    // histos.add("dEdx_vs_centr_vs_occup_narrow_p_win_neg_FractionPIDclsInRange", "dE/dx", kTH3F, {{20, 0, 4000, "nITStrk cls567"}, axisOccupancyForDeDxStudies, axisDeDx});

    histos.add("dEdx_3OROC_max_vs_centr_vs_occup_narrow_p_win", "dE/dx", kTH3F, {{20, 0, 4000, "nITStrk cls567"}, axisOccupancyForDeDxStudies, axisDeDx});
    histos.add("dEdx_3OROC_tot_vs_centr_vs_occup_narrow_p_win", "dE/dx", kTH3F, {{20, 0, 4000, "nITStrk cls567"}, axisOccupancyForDeDxStudies, axisDeDx});

    histos.add("dEdx_vs_centr_vs_occup_narrow_p_win_CORRECTED", "dE/dx", kTH3F, {{20, 0, 4000, "nITStrk cls567"}, axisOccupancyForDeDxStudies, axisDeDx});

    // AxisSpec axisFractionNclsFindableMinusPID{110, -1.1, 1.1, "TPC nClsFindableMinusPID / nClsFindable"};
    // histos.add("fraction_tpcNClsFindableMinusPID_vs_occup", "", kTH2D, {axisOccupancyForDeDxStudies, axisFractionNclsFindableMinusPID});
    // histos.add("fraction_tpcNClsFindableMinusPID_vs_occup_peripheralByV0A", "", kTH2D, {axisOccupancyForDeDxStudies, axisFractionNclsFindableMinusPID});
    // histos.add("fraction_tpcNClsFindableMinusPID_vs_occup_centralByV0A", "", kTH2D, {axisOccupancyForDeDxStudies, axisFractionNclsFindableMinusPID});
    // histos.add("fraction_tpcNClsFindableMinusPID_vs_occup_eta02", "", kTH2D, {axisOccupancyForDeDxStudies, axisFractionNclsFindableMinusPID});
    // histos.add("fraction_tpcNClsFindableMinusPID_vs_occup_pos", "", kTH2D, {axisOccupancyForDeDxStudies, axisFractionNclsFindableMinusPID});
    // histos.add("fraction_tpcNClsFindableMinusPID_vs_occup_neg", "", kTH2D, {axisOccupancyForDeDxStudies, axisFractionNclsFindableMinusPID});
    // histos.add("fraction_tpcNClsFindableMinusPID_vs_occup_lowPt", "", kTH2D, {axisOccupancyForDeDxStudies, axisFractionNclsFindableMinusPID});
    // histos.add("fraction_tpcNClsFindableMinusPID_vs_occup_highPt", "", kTH2D, {axisOccupancyForDeDxStudies, axisFractionNclsFindableMinusPID});

    // QA dEdx correction coeff
    histos.add("dEdx_CORRECTION_COEFF", "coeff", kTH1F, {{1000, -2.0, 2.0, "correction coeff"}});
  }

  Preslice<FullTracksIU> perCollision = aod::track::collisionId;

  float fReal_fTPCSignalN(float mbb0R1, float a1pt, float atgl, float atglmbb0R1, float a1ptmbb0R1, float side, float a1pt2, float fTrackOccN, float fOccTPCN, float fTrackOccMeanN)
  {
    return ((0.017012 * mbb0R1) + (-0.0018469 * a1pt) + (-0.0052177 * atgl) + (-0.0035655 * atglmbb0R1) + (0.0017846 * a1ptmbb0R1) + (0.0019127 * side) + (-0.00012964 * a1pt2) + (0.013066)) * fTrackOccN + ((0.0055592 * mbb0R1) + (-0.0010618 * a1pt) + (-0.0016134 * atgl) + (-0.0059098 * atglmbb0R1) + (0.0013335 * a1ptmbb0R1) + (0.00052133 * side) + (3.1119e-05 * a1pt2) + (0.0049428)) * fOccTPCN + ((0.00077317 * mbb0R1) + (-0.0013827 * a1pt) + (0.003249 * atgl) + (-0.00063689 * atglmbb0R1) + (0.0016218 * a1ptmbb0R1) + (-0.00045215 * side) + (-1.5815e-05 * a1pt2) + (-0.004882)) * fTrackOccMeanN + ((-0.015053 * mbb0R1) + (0.0018912 * a1pt) + (-0.012305 * atgl) + (0.081387 * atglmbb0R1) + (0.003205 * a1ptmbb0R1) + (-0.0087404 * side) + (-0.0028608 * a1pt2) + (0.99091));
  };
  void processRun3(
    ColEvSels const& cols,
    FullTracksIU const& tracks,
    BCsRun3 const& bcs,
    aod::TracksQA_002 const& tracksQA,
    aod::FT0s const&)
  {
    int runNumber = bcs.iteratorAt(0).runNumber();
    if (runNumber != lastRunNumber) {
      lastRunNumber = runNumber; // do it only once

      if (runNumber >= 500000) {
        auto runInfo = o2::parameters::AggregatedRunInfo::buildAggregatedRunInfo(o2::ccdb::BasicCCDBManager::instance(), runNumber);
        // first bc of the first orbit
        bcSOR = runInfo.orbitSOR * o2::constants::lhc::LHCMaxBunches;
        // duration of TF in bcs
        nBCsPerTF = runInfo.orbitsPerTF * o2::constants::lhc::LHCMaxBunches;

        LOGP(info, "bcSOR = {}, nBCsPerTF = {}", bcSOR, nBCsPerTF);
      }
    }

    // track QA table
    std::vector<int64_t> labelTrack2TrackQA;
    labelTrack2TrackQA.clear();
    labelTrack2TrackQA.resize(tracks.size(), -1);
    for (const auto& trackQA : tracksQA) {
      int64_t trackId = trackQA.trackId();
      int64_t trackQAIndex = trackQA.globalIndex();
      labelTrack2TrackQA[trackId] = trackQAIndex;
    }

    for (const auto& col : cols) {
      if (!col.sel8())
        continue;

      // check hadronic rate
      auto bc = col.foundBC_as<BCsRun3>();
      int64_t ts = bc.timestamp();
      double hadronicRate = mRateFetcher.fetch(ccdb.service, ts, runNumber, "ZNC hadronic") * 1.e-3; // kHz
      const int multTPC = col.multTPC();

      int occupancy = col.trackOccupancyInTimeRange();

      auto tracksGrouped = tracks.sliceBy(perCollision, col.globalIndex());

      // pre-calc nPV
      int nPV = 0;
      for (const auto& track : tracksGrouped) {
        if (!track.isPVContributor())
          continue;
        if (track.pt() < confCutPtMinThisEvent || track.pt() > confCutPtMaxThisEvent)
          continue;
        if (track.eta() < confCutEtaMinTracksThisEvent || track.eta() > confCutEtaMaxTracksThisEvent)
          continue;
        if (track.itsNCls() < 5)
          continue;
        nPV++;
      }

      // main track loop
      for (const auto& track : tracksGrouped) {
        if (!track.isPVContributor())
          continue;
        if (track.pt() < confCutPtMinThisEvent || track.pt() > confCutPtMaxThisEvent)
          continue;
        if (track.eta() < confCutEtaMinTracksThisEvent || track.eta() > confCutEtaMaxTracksThisEvent)
          continue;
        if (track.itsNCls() < 5)
          continue;
        if (!track.isGlobalTrack())
          continue;
        if (track.tpcNClsFound() < confCutMinTPCcls)
          continue;

        float signedP = track.sign() * track.tpcInnerParam();

        aod::TracksQA_002::iterator trackQA;
        // bool existPosTrkQA;
        if (labelTrack2TrackQA[track.globalIndex()] != -1) {
          trackQA = tracksQA.iteratorAt(labelTrack2TrackQA[track.globalIndex()]);
          // existPosTrkQA = true;

          float signedP = track.sign() * track.tpcInnerParam();
          float dEdx = track.tpcSignal();

          float tpcdEdxMax0Rabs = trackQA.tpcdEdxMax0R() * dEdx / 100;
          float tpcdEdxMax1Rabs = trackQA.tpcdEdxMax1R() * dEdx / 100;
          float tpcdEdxMax2Rabs = trackQA.tpcdEdxMax2R() * dEdx / 100;
          float tpcdEdxMax3Rabs = trackQA.tpcdEdxMax3R() * dEdx / 100;

          float tpcdEdxTot0Rabs = trackQA.tpcdEdxTot0R() * dEdx / 100;
          float tpcdEdxTot1Rabs = trackQA.tpcdEdxTot1R() * dEdx / 100;
          float tpcdEdxTot2Rabs = trackQA.tpcdEdxTot2R() * dEdx / 100;
          float tpcdEdxTot3Rabs = trackQA.tpcdEdxTot3R() * dEdx / 100;

          histos.fill(HIST("tpcdEdxMax0R_vs_Momentum"), signedP, tpcdEdxMax0Rabs);
          histos.fill(HIST("tpcdEdxMax1R_vs_Momentum"), signedP, tpcdEdxMax1Rabs);
          histos.fill(HIST("tpcdEdxMax2R_vs_Momentum"), signedP, tpcdEdxMax2Rabs);
          histos.fill(HIST("tpcdEdxMax3R_vs_Momentum"), signedP, tpcdEdxMax3Rabs);

          histos.fill(HIST("tpcdEdxTot0R_vs_Momentum"), signedP, tpcdEdxTot0Rabs);
          histos.fill(HIST("tpcdEdxTot1R_vs_Momentum"), signedP, tpcdEdxTot1Rabs);
          histos.fill(HIST("tpcdEdxTot2R_vs_Momentum"), signedP, tpcdEdxTot2Rabs);
          histos.fill(HIST("tpcdEdxTot3R_vs_Momentum"), signedP, tpcdEdxTot3Rabs);

          // FROM: https://github.com/AliceO2Group/AliceO2/blob/d4afff4276fae2d31f6c3c79d9ec4246deff95f8/Detectors/AOD/src/AODProducerWorkflowSpec.cxx#L2628C1-L2629C84
          // const float dEdxNorm = (tpcOrig.getdEdx().dEdxTotTPC > 0) ? 100. / tpcOrig.getdEdx().dEdxTotTPC : 0;
          // trackQAHolder.tpcdEdxMax0R = uint8_t(tpcOrig.getdEdx().dEdxMaxIROC * dEdxNorm);
          histos.fill(HIST("tpcdEdxTot3R_vs_dEdxFromTracks"), track.tpcSignal(), trackQA.tpcdEdxTot3R() * dEdx / 100);
          histos.fill(HIST("tpcdEdxTotSUM_vs_dEdxFromTracks"), track.tpcSignal(), tpcdEdxTot0Rabs + tpcdEdxTot1Rabs + tpcdEdxTot2Rabs + tpcdEdxTot3Rabs);
          histos.fill(HIST("tpcdEdxMaxSUM_vs_dEdxFromTracks"), track.tpcSignal(), tpcdEdxMax0Rabs + tpcdEdxMax1Rabs + tpcdEdxMax2Rabs + tpcdEdxMax3Rabs);

          // ### dEdx MAX
          if (1) {
            float sum_dEdx_max = 0;
            int counter_has_dEdx_max = 0;
            if (tpcdEdxMax1Rabs > 0) {
              sum_dEdx_max += tpcdEdxMax1Rabs;
              counter_has_dEdx_max++;
            }
            if (tpcdEdxMax2Rabs > 0) {
              sum_dEdx_max += tpcdEdxMax2Rabs;
              counter_has_dEdx_max++;
            }
            if (tpcdEdxMax3Rabs > 0) {
              sum_dEdx_max += tpcdEdxMax3Rabs;
              counter_has_dEdx_max++;
            }
            // only 3 OROC:
            float sum_3OROC_dEdx_max = sum_dEdx_max;
            int counter_3OROC_has_dEdx_max = counter_has_dEdx_max;
            if (counter_3OROC_has_dEdx_max > 0) {
              sum_3OROC_dEdx_max /= counter_3OROC_has_dEdx_max;
              histos.fill(HIST("tpcdEdxAverageMax_3OROC_vs_dEdxFromTracks"), track.tpcSignal(), sum_3OROC_dEdx_max);
            }
            // now IROC:
            if (tpcdEdxMax0Rabs > 0) {
              sum_dEdx_max += tpcdEdxMax0Rabs;
              counter_has_dEdx_max++;
            }
            // average and fill histos
            if (counter_has_dEdx_max > 0) {
              sum_dEdx_max /= counter_has_dEdx_max;
              histos.fill(HIST("tpcdEdxAverageMax_vs_dEdxFromTracks"), track.tpcSignal(), sum_dEdx_max);
            }
            if (occupancy >= 0) {
              if (std::fabs(signedP) > 0.38 && std::fabs(signedP) < 0.4)
                histos.fill(HIST("dEdx_3OROC_max_vs_centr_vs_occup_narrow_p_win"), nPV, occupancy, sum_dEdx_max);

              if (track.eta() > 0.2 && track.eta() < 0.4)
                histos.fill(HIST("dEdx_3OROC_max_vs_Momentum_vs_occup_eta_02_04"), signedP, sum_dEdx_max, occupancy);
              if (track.eta() > -0.4 && track.eta() < -0.2)
                histos.fill(HIST("dEdx_3OROC_max_vs_Momentum_vs_occup_eta_04_02"), signedP, sum_dEdx_max, occupancy);
            }
          }
          // ### dEdx TOT
          if (1) {
            float sum_dEdx_tot = 0;
            int counter_has_dEdx_tot = 0;
            if (tpcdEdxTot1Rabs > 0) {
              sum_dEdx_tot += tpcdEdxTot1Rabs;
              counter_has_dEdx_tot++;
            }
            if (tpcdEdxTot2Rabs > 0) {
              sum_dEdx_tot += tpcdEdxTot2Rabs;
              counter_has_dEdx_tot++;
            }
            if (tpcdEdxTot3Rabs > 0) {
              sum_dEdx_tot += tpcdEdxTot3Rabs;
              counter_has_dEdx_tot++;
            }
            // only 3 OROC:
            float sum_3OROC_dEdx_tot = sum_dEdx_tot;
            int counter_3OROC_has_dEdx_tot = counter_has_dEdx_tot;
            if (counter_3OROC_has_dEdx_tot > 0) {
              sum_3OROC_dEdx_tot /= counter_3OROC_has_dEdx_tot;
              histos.fill(HIST("tpcdEdxAverageTot_3OROC_vs_dEdxFromTracks"), track.tpcSignal(), sum_3OROC_dEdx_tot);
            }
            // now IROC:
            if (tpcdEdxTot0Rabs > 0) {
              sum_dEdx_tot += tpcdEdxTot0Rabs;
              counter_has_dEdx_tot++;
            }
            // average and fill histos
            if (counter_has_dEdx_tot > 0) {
              sum_dEdx_tot /= counter_has_dEdx_tot;
              histos.fill(HIST("tpcdEdxAverageTot_vs_dEdxFromTracks"), track.tpcSignal(), sum_dEdx_tot);
            }

            if (occupancy >= 0) {
              if (std::fabs(signedP) > 0.38 && std::fabs(signedP) < 0.4)
                histos.fill(HIST("dEdx_3OROC_tot_vs_centr_vs_occup_narrow_p_win"), nPV, occupancy, sum_dEdx_tot);

              if (track.eta() > 0.2 && track.eta() < 0.4)
                histos.fill(HIST("dEdx_3OROC_tot_vs_Momentum_vs_occup_eta_02_04"), signedP, sum_dEdx_tot, occupancy);
              if (track.eta() > -0.4 && track.eta() < -0.2)
                histos.fill(HIST("dEdx_3OROC_tot_vs_Momentum_vs_occup_eta_04_02"), signedP, sum_dEdx_tot, occupancy);
            }
          }

          histos.fill(HIST("dcaXY_vs_dcaXYqa"), track.dcaXY(), trackQA.tpcdcaR());
          histos.fill(HIST("dcaZ_vs_dcaZqa"), track.dcaZ(), trackQA.tpcdcaZ());
        }
        // else {
        //   existPosTrkQA = false;
        // }

        // ### dE/dx by Marian:
        float fTPCSignal = track.tpcSignal();
        float fNormMultTPC = multTPC / 11000.; // IA: my guess: it's https://github.com/AliceO2Group/O2Physics/blob/f681d9cc71214c4eb5613a3f473cbea41e48a61f/DPG/Tasks/TPC/tpcSkimsTableCreator.cxx#L575C30-L575C47

        // df["mdEdx"]=(50/df["fTPCSignal"]).clip(0.05,1.1)
        // df["fTPCSignalN"]=(df["fTPCSignal"]/df["bb0"]/50.).clip(0.5,1.5)
        // df["fTrackOccN"]=df.eval("fTrackOcc/1000.")
        // df["mdEdxExp"]=df.eval("1./bb0")
        // df["fFt0OccN"]=df["fFt0Occ"]*df.eval("fFt0Occ/fTrackOcc").median()
        // df["mdEdxExpOcc"]=df.eval("mdEdxExp*fTrackOccN")
        // df["fTrackOccMeanN"]=(df["fHadronicRate"]/5)                 # normalization 5 - 10 bins
        // df["fTrackOccN2"]=df.eval("fTrackOccN*fTrackOccN")
        // df["fOccTPCN"]=(df["fNormMultTPC"]*10).clip(0,12)           # normalization 10 - 12 bins
        // df["mdEdxOccTPCN"]=df.eval("mdEdx*fOccTPCN")
        // df["mdEdxMeanOccTPCN"]=df.eval("mdEdx*fTrackOccMeanN")

        float fTrackOccN = occupancy / 1000.;
        float fOccTPCN = fNormMultTPC * 10; //(fNormMultTPC*10).clip(0,12)
        if (fOccTPCN > 12)
          fOccTPCN = 12;
        else if (fOccTPCN < 0)
          fOccTPCN = 0;

        float fTrackOccMeanN = hadronicRate / 5;

        float side = track.tgl() > 0 ? 1 : 0;
        float a1pt = std::abs(track.signed1Pt());
        float a1pt2 = a1pt * a1pt;
        float atgl = std::abs(track.tgl());
        float mbb0R = 50 / fTPCSignal;
        if (mbb0R > 1.05)
          mbb0R = 1.05;
        else if (mbb0R < 0.05)
          mbb0R = 0.05;
        // float mbb0R =  max(0.05,  min(50 / fTPCSignal, 1.05));
        float a1ptmbb0R = a1pt * mbb0R;
        float atglmbb0R = atgl * mbb0R;

        // tree->SetAlias("side","fTgl>0");
        // tree->SetAlias("a1pt","abs(fSigned1Pt)");
        // tree->SetAlias("a1pt2","abs(fSigned1Pt**2)");
        // tree->SetAlias("atgl","abs(fTgl)");
        // tree->SetAlias("mbb0R","max(0.05,min(50/fTPCSignal,1.05))");
        // tree->SetAlias("a1ptmbb0R","a1pt*mbb0R");
        // tree->SetAlias("atglmbb0R","atgl*mbb0R");

        // ### iteration 1 correction
        // float fTPCSignalN_CBB = fReal_fTPCSignalN(mbb0,a1pt,atgl,atglmbb0,a1ptmbb0,side,a1pt2,fTrackOccN,fOccTPCN,fTrackOccMeanN+0); // atglmbb0 is != atglmbb0R!!! etc.
        float fTPCSignalN_CR0 = fReal_fTPCSignalN(mbb0R, a1pt, atgl, atglmbb0R, a1ptmbb0R, side, a1pt2, fTrackOccN, fOccTPCN, fTrackOccMeanN + 0);

        // tree->SetAlias("fTPCSignalN_CBB","fReal_fTPCSignalN(mbb0,a1pt,atgl,atglmbb0,a1ptmbb0,side,a1pt2,fTrackOccN,fOccTPCN,fTrackOccMeanN+0)");
        // tree->SetAlias("fTPCSignalN_CR0","fReal_fTPCSignalN(mbb0R,a1pt,atgl,atglmbb0R,a1ptmbb0R,side,a1pt2,fTrackOccN,fOccTPCN,fTrackOccMeanN+0)");

        float mbb0R1 = 50 / (fTPCSignal / fTPCSignalN_CR0);
        if (mbb0R1 > 1.05)
          mbb0R1 = 1.05;
        else if (mbb0R1 < 0.05)
          mbb0R1 = 0.05;
        // float mbb0R1 = max(0.05, min(50 / (fTPCSignal / fTPCSignalN_CR0), 1.05 + 0));
        // tree->SetAlias("mbb0R1","max(0.05,min(50/(fTPCSignal/fTPCSignalN_CR0),1.05+0))");
        float fTPCSignalN_CR1 = fReal_fTPCSignalN(mbb0R1, a1pt, atgl, atgl * mbb0R1, a1pt * mbb0R1, side, a1pt2, fTrackOccN, fOccTPCN, fTrackOccMeanN + 0);
        // tree->SetAlias("fTPCSignalN_CR1","fReal_fTPCSignalN(mbb0R1,a1pt,atgl,atgl*mbb0R1,a1pt*mbb0R1,side,a1pt2,fTrackOccN,fOccTPCN,fTrackOccMeanN+0)");
        //
        // tree->SetAlias("fTPCSignalN_mad_BB","fReal_fTPCSignalN_mad(mbb0,a1pt,atgl,atglmbb0,a1ptmbb0,side,a1pt2,fTrackOccN,fOccTPCN,fTrackOccMeanN+0)");
        // tree->SetAlias("fTPCSignalN_mad_R0","fReal_fTPCSignalN_mad(mbb0R1,a1pt,atgl,atgl*mbb0R1,a1pt*mbb0R1,side,a1pt2,fTrackOccN,fOccTPCN,fTrackOccMeanN+0)");
        //
        // tree->SetAlias("fTPCSignal_CorrR1","fTPCSignal/fTPCSignalN_CR1");
        // tree->SetAlias("fTPCSignal_CorrBB","fTPCSignal/fTPCSignalN_CBB");

        histos.fill(HIST("dEdx_vs_Momentum"), signedP, fTPCSignal);

        float corrected_dEdx = fTPCSignal / fTPCSignalN_CR1;
        histos.fill(HIST("dEdx_CORRECTION_COEFF"), fTPCSignalN_CR1);
        histos.fill(HIST("dEdx_vs_Momentum_CORRECTED"), signedP, corrected_dEdx);
        histos.fill(HIST("tpcdEdxCORRECTED_vs_dEdxFromTracks"), fTPCSignal, corrected_dEdx);

        if (occupancy >= 0) {
          if (std::fabs(signedP) > 0.38 && std::fabs(signedP) < 0.4) {
            histos.fill(HIST("dEdx_vs_centr_vs_occup_narrow_p_win"), nPV, occupancy, fTPCSignal);
            histos.fill(HIST("dEdx_vs_centr_vs_occup_narrow_p_win_CORRECTED"), nPV, occupancy, corrected_dEdx);
          }
        }

      } // end of track loop
    } // end of collision loop
  }
  PROCESS_SWITCH(dEdxVsOccupancyWithTrackQAinfoTask, processRun3, "Process Run3 tracking vs detector occupancy QA", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<dEdxVsOccupancyWithTrackQAinfoTask>(cfgc)};
}
