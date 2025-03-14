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

/// \file timeDependentQa.cxx
/// \brief Time-dependent QA for a number of observables
///
/// \author Evgeny Kryshen <evgeny.kryshen@cern.ch> and Igor Altsybeev <Igor.Altsybeev@cern.ch>

#include <map>
#include <vector>
#include <string>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"
#include "CCDB/BasicCCDBManager.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/CCDB/EventSelectionParams.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/CCDB/ctpRateFetcher.h"
#include "TPCCalibration/TPCMShapeCorrection.h"
#include "DataFormatsParameters/AggregatedRunInfo.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "ReconstructionDataFormats/Vertex.h"
#include "Common/DataModel/Multiplicity.h"

#include "TTree.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::aod::evsel;

using ColEvSels = soa::Join<aod::Collisions, aod::EvSels, aod::Mults>;
using BCsRun3 = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels>;
using BarrelTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>;

const AxisSpec axisQoverPt{100, -1., 1., "q/p_{T}, 1/GeV"};
const AxisSpec axisDcaR{1000, -5., 5., "DCA_{r}, cm"};
const AxisSpec axisDcaZ{1000, -5., 5., "DCA_{z}, cm"};
const AxisSpec axisSparseQoverPt{20, -1., 1., "q/p_{T}, 1/GeV"};
const AxisSpec axisSparseDcaR{100, -1., 1., "DCA_{r}, cm"};
const AxisSpec axisSparseDcaZ{100, -1., 1., "DCA_{z}, cm"};

struct TimeDependentQaTask {
  Configurable<float> confTimeBinWidthInSec{"TimeBinWidthInSec", 0.25, "Width of time bins in seconds"};                                                                         // o2-linter: disable=name/configurable (temporary fix)
  Configurable<int> confTakeVerticesWithUPCsettings{"ConsiderVerticesWithUPCsettings", 0, "Take vertices: 0 - all , 1 - only without UPC settings, 2 - only with UPC settings"}; // o2-linter: disable=name/configurable (temporary fix)
  Configurable<int> confFillPhiVsTimeHist{"FlagFillPhiVsTimeHist", 2, "0 - don't fill , 1 - fill only for global/7cls/TRD/TOF tracks, 2 - fill also layer-by-layer"};            // o2-linter: disable=name/configurable (temporary fix)
  Configurable<int> confFillEtaPhiVsTimeHist{"FlagFillEtaPhiVsTimeHist", 0, "0 - don't fill , 1 - fill"};                                                                        // o2-linter: disable=name/configurable (temporary fix)
  Configurable<float> confCutOnNtpcClsForSharedFractAndDeDxCalc{"CutOnNtpcClsForSharedFractAndDeDxCalc", 70, ""};                                                                // o2-linter: disable=name/configurable (temporary fix)

  enum EvSelBitsToMonitor {
    enCollisionsAll,
    enIsTriggerTVX,
    enNoTimeFrameBorder,
    enNoITSROFrameBorder,
    enCollisionsSel8,
    enNoSameBunchPileup,
    enIsGoodZvtxFT0vsPV,
    enIsVertexITSTPC,
    enIsVertexTOFmatched,
    enIsVertexTRDmatched,
    enNoCollInTimeRangeNarrow,
    enNoCollInTimeRangeStrict,
    enNoCollInTimeRangeStandard,
    enNoCollInRofStrict,
    enNoCollInRofStandard,
    enNoHighMultCollInPrevRof,
    enIsGoodITSLayer3,
    enIsGoodITSLayer0123,
    enIsGoodITSLayersAll,
    enIsLowOccupStd,
    enIsLowOccupStdAlsoInPrevRof,
    enIsLowOccupStdCut500,
    enIsLowOccupStdCut2000,
    enIsLowOccupStdCut4000,
    enIsLowOccupStdAlsoInPrevRofCut2000noDeadStaves,
    enNumEvSelBits, // counter
  };

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  o2::tpc::TPCMShapeCorrection mshape; // object for simple access
  int lastRunNumber = -1;
  double maxSec = 1;
  double minSec = 0;
  static const int32_t nBCsPerOrbit = o2::constants::lhc::LHCMaxBunches;
  int64_t bcSOR = 0;      // global bc of the start of the first orbit, setting 0 for unanchored MC
  int64_t nBCsPerTF = -1; // duration of TF in bcs
  ctpRateFetcher mRateFetcher;

  void init(InitContext&)
  {
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    histos.add("allTracks/hQoverPt", "", kTH1F, {axisQoverPt});
    histos.add("allTracks/hDcaR", "", kTH1F, {axisDcaR});
    histos.add("allTracks/hDcaZ", "", kTH1F, {axisDcaZ});
    histos.add("allTracks/hDcaRafterCuts", "", kTH1F, {axisDcaR});
    histos.add("allTracks/hDcaZafterCuts", "", kTH1F, {axisDcaZ});
    histos.add("allTracks/hQoverPtDcaR", "", kTH2F, {axisSparseQoverPt, axisSparseDcaR});
    histos.add("allTracks/hQoverPtDcaZ", "", kTH2F, {axisSparseQoverPt, axisSparseDcaZ});

    histos.add("PVcontrib/hDcaRafterCuts", "", kTH1F, {axisDcaR});
    histos.add("PVcontrib/hDcaZafterCuts", "", kTH1F, {axisDcaZ});

    histos.add("A/global/hDcaRafterCuts", "", kTH1F, {axisDcaR});
    histos.add("A/global/hDcaZafterCuts", "", kTH1F, {axisDcaZ});
    histos.add("A/globalPV/hDcaRafterCuts", "", kTH1F, {axisDcaR});
    histos.add("A/globalPV/hDcaZafterCuts", "", kTH1F, {axisDcaZ});

    histos.add("C/global/hDcaRafterCuts", "", kTH1F, {axisDcaR});
    histos.add("C/global/hDcaZafterCuts", "", kTH1F, {axisDcaZ});
    histos.add("C/globalPV/hDcaRafterCuts", "", kTH1F, {axisDcaR});
    histos.add("C/globalPV/hDcaZafterCuts", "", kTH1F, {axisDcaZ});

    const AxisSpec axisBCinTF{150000, 0, 150000, "bc in TF"};
    histos.add("hNcolVsBcInTF", ";bc in TF; n collisions", kTH1F, {axisBCinTF});
    histos.add("hNcolVsBcInTFantiBorderCut", ";bc in TF; n collisions", kTH1F, {axisBCinTF});
  }

  void processRun3(
    ColEvSels const& cols,
    BarrelTracks const& tracks,
    BCsRun3 const& bcs,
    aod::FT0s const&)
  {
    int runNumber = bcs.iteratorAt(0).runNumber();
    if (runNumber != lastRunNumber) {
      LOGP(debug, "  >> QA: run number = {}", runNumber);
      lastRunNumber = runNumber;

      int64_t tsSOR = 0; // dummy start-of-run timestamp
      int64_t tsEOR = 1; // dummy end-of-run timestamp
      if (runNumber >= 500000) {
        auto runInfo = o2::parameters::AggregatedRunInfo::buildAggregatedRunInfo(o2::ccdb::BasicCCDBManager::instance(), runNumber);
        // first bc of the first orbit
        bcSOR = runInfo.orbitSOR * nBCsPerOrbit;
        // duration of TF in bcs
        nBCsPerTF = runInfo.orbitsPerTF * nBCsPerOrbit;
        // start-of-run timestamp
        tsSOR = runInfo.sor;
        // end-of-run timestamp
        tsEOR = runInfo.eor;
      }

      minSec = floor(tsSOR / 1000.);
      maxSec = ceil(tsEOR / 1000.);
      int nTimeBins = static_cast<int>((maxSec - minSec) / confTimeBinWidthInSec);
      double timeInterval = nTimeBins * confTimeBinWidthInSec;

      const AxisSpec axisSeconds{nTimeBins, 0, timeInterval, "seconds"};
      histos.add("hSecondsBCsTVX", "", kTH1D, {axisSeconds});
      // histos.add("hSecondsBCsTFborder", "", kTH1D, {axisSeconds});
      histos.add("hSecondsBCsTVXandTFborder", "", kTH1D, {axisSeconds});

      histos.add("hSecondsCollisionsBeforeAllCuts", "", kTH1D, {axisSeconds});
      histos.add("hSecondsCollisionsNoVzInTVX", "", kTH1D, {axisSeconds});
      histos.add("hSecondsCollisionsNoVzNoTFborder", "", kTH1D, {axisSeconds});
      histos.add("hSecondsCollisionsNoVzInTVXandNoTFborder", "", kTH1D, {axisSeconds});

      histos.add("hSecondsCollisions", "", kTH1D, {axisSeconds});
      histos.add("hSecondsIR", "", kTH1D, {axisSeconds});
      histos.add("hSecondsVz", "", kTH1D, {axisSeconds});
      histos.add("hSecondsFT0Camlp", "", kTH1D, {axisSeconds});
      histos.add("hSecondsFT0CamlpByColMult", "", kTH1D, {axisSeconds});
      histos.add("hSecondsFT0AamlpByColMult", "", kTH1D, {axisSeconds});
      histos.add("hSecondsV0Aamlp", "", kTH1D, {axisSeconds});
      histos.add("hSecondsOccupancyByTracks", "", kTH1D, {axisSeconds});
      histos.add("hSecondsOccupancyByFT0C", "", kTH1D, {axisSeconds});

      // QA for UPC settings
      histos.add("hSecondsUPCverticesBeforeSel8", "", kTH2F, {axisSeconds, {2, -0.5, 1.5, "Is vertex with UPC settings"}});
      histos.add("hSecondsUPCvertices", "", kTH2F, {axisSeconds, {2, -0.5, 1.5, "Is vertex with UPC settings after sel8"}});

      // ### QA event selection bits
      int nEvSelBits = enNumEvSelBits;
      histos.add("hSecondsEventSelBits", "", kTH2F, {axisSeconds, {nEvSelBits, -0.5, nEvSelBits - 0.5, "Monitoring of event selection bits"}});
      TAxis* axSelBits = reinterpret_cast<TAxis*>(histos.get<TH2>(HIST("hSecondsEventSelBits"))->GetYaxis());
      axSelBits->SetBinLabel(1 + enCollisionsAll, "collisionsAll");
      axSelBits->SetBinLabel(1 + enIsTriggerTVX, "IsTriggerTVX");
      axSelBits->SetBinLabel(1 + enNoTimeFrameBorder, "NoTimeFrameBorder");
      axSelBits->SetBinLabel(1 + enNoITSROFrameBorder, "NoITSROFrameBorder");

      // bits after sel8
      axSelBits->SetBinLabel(1 + enCollisionsSel8, "collisionsSel8");
      axSelBits->SetBinLabel(1 + enNoSameBunchPileup, "NoSameBunchPileup");
      axSelBits->SetBinLabel(1 + enIsGoodZvtxFT0vsPV, "IsGoodZvtxFT0vsPV");
      axSelBits->SetBinLabel(1 + enIsVertexITSTPC, "IsVertexITSTPC");
      axSelBits->SetBinLabel(1 + enIsVertexTOFmatched, "IsVertexTOFmatched");
      axSelBits->SetBinLabel(1 + enIsVertexTRDmatched, "IsVertexTRDmatched");

      axSelBits->SetBinLabel(1 + enNoCollInTimeRangeNarrow, "NoCollInTimeRangeNarrow");
      axSelBits->SetBinLabel(1 + enNoCollInTimeRangeStrict, "NoCollInTimeRangeStrict");
      axSelBits->SetBinLabel(1 + enNoCollInTimeRangeStandard, "NoCollInTimeRangeStandard");
      axSelBits->SetBinLabel(1 + enNoCollInRofStrict, "NoCollInRofStrict");
      axSelBits->SetBinLabel(1 + enNoCollInRofStandard, "NoCollInRofStandard");
      axSelBits->SetBinLabel(1 + enNoHighMultCollInPrevRof, "NoHighMultCollInPrevRof");

      axSelBits->SetBinLabel(1 + enIsGoodITSLayer3, "IsGoodITSLayer3");
      axSelBits->SetBinLabel(1 + enIsGoodITSLayer0123, "IsGoodITSLayer0123");
      axSelBits->SetBinLabel(1 + enIsGoodITSLayersAll, "IsGoodITSLayersAll");

      // combined conditions on occupancy
      axSelBits->SetBinLabel(1 + enIsLowOccupStd, "isLowOccupStd");
      axSelBits->SetBinLabel(1 + enIsLowOccupStdAlsoInPrevRof, "isLowOccupStdAlsoInPrevRof");
      axSelBits->SetBinLabel(1 + enIsLowOccupStdCut500, "isLowOccupStdCut500");
      axSelBits->SetBinLabel(1 + enIsLowOccupStdCut2000, "isLowOccupStdCut2000");
      axSelBits->SetBinLabel(1 + enIsLowOccupStdCut4000, "isLowOccupStdCut4000");
      axSelBits->SetBinLabel(1 + enIsLowOccupStdAlsoInPrevRofCut2000noDeadStaves, "isLowOccupStdAlsoInPrevRofCut2000noDeadStaves");

      // const AxisSpec axisChi2ITS{40, 0., 20., "chi2/ndof"};
      // const AxisSpec axisChi2TPC{40, 0., 20., "chi2/ndof"};
      const AxisSpec axisNclsITS{5, 3.5, 8.5, "n ITS cls"};
      const AxisSpec axisNclsTPC{40, -0.5, 159.5, "n TPC cls"};
      const AxisSpec axisFraction{20, 0, 1., "Fraction shared cls Tpc"};

      // QA for all tracks
      histos.add("allTracks/hSecondsTracks", "", kTH1D, {axisSeconds});
      histos.add("allTracks/hSecondsQoverPtSumDcaR", "", kTH2D, {axisSeconds, axisSparseQoverPt});
      histos.add("allTracks/hSecondsQoverPtSumDcaZ", "", kTH2D, {axisSeconds, axisSparseQoverPt});
      histos.add("allTracks/hSecondsSumDcaR", "", kTH1D, {axisSeconds});
      histos.add("allTracks/hSecondsSumDcaZ", "", kTH1D, {axisSeconds});
      histos.add("allTracks/hSecondsSumPt", "", kTH1D, {axisSeconds});
      histos.add("allTracks/hSecondsNumClsIts", "", kTH2D, {axisSeconds, axisNclsITS});
      histos.add("allTracks/hSecondsChi2NClIts", "", kTH1D, {axisSeconds});
      histos.add("allTracks/hSecondsTracksMshape", "", kTH1D, {axisSeconds});

      // QA for PV contributors
      histos.add("PVcontrib/hSecondsTracks", "", kTH1D, {axisSeconds});
      // histos.add("PVcontrib/hSecondsQoverPtSumDcaR", "", kTH2D, {axisSeconds, axisSparseQoverPt});
      // histos.add("PVcontrib/hSecondsQoverPtSumDcaZ", "", kTH2D, {axisSeconds, axisSparseQoverPt});
      histos.add("PVcontrib/hSecondsSumDcaR", "", kTH1D, {axisSeconds});
      histos.add("PVcontrib/hSecondsSumDcaZ", "", kTH1D, {axisSeconds});
      histos.add("PVcontrib/hSecondsSumPt", "", kTH1D, {axisSeconds});
      histos.add("PVcontrib/hSecondsNumClsIts", "", kTH2D, {axisSeconds, axisNclsITS});
      histos.add("PVcontrib/hSecondsChi2NClIts", "", kTH1D, {axisSeconds});

      // QA for global tracks
      // ### A side
      // global tracks
      histos.add("A/global/hSecondsNumTracks", "", kTH1D, {axisSeconds});
      histos.add("A/global/hSecondsQoverPtSumDcaR", "", kTH2D, {axisSeconds, axisSparseQoverPt});
      histos.add("A/global/hSecondsQoverPtSumDcaZ", "", kTH2D, {axisSeconds, axisSparseQoverPt});
      histos.add("A/global/hSecondsSumDcaR", "", kTH1D, {axisSeconds});
      histos.add("A/global/hSecondsSumDcaZ", "", kTH1D, {axisSeconds});
      histos.add("A/global/hSecondsSumPt", "", kTH1D, {axisSeconds});
      histos.add("A/global/hSecondsNumClsIts", "", kTH2D, {axisSeconds, axisNclsITS});
      histos.add("A/global/hSecondsChi2NClIts", "", kTH1D, {axisSeconds});
      histos.add("A/global/hSecondsNumClsTpc", "", kTH2D, {axisSeconds, axisNclsTPC});
      histos.add("A/global/hSecondsChi2NClTpc", "", kTH1D, {axisSeconds});
      histos.add("A/global/hSecondsTpcFractionSharedCls", "", kTH2D, {axisSeconds, axisFraction});
      histos.add("A/global/hSecondsDeDx", "", kTH1D, {axisSeconds});

      // global && PV tracks
      histos.add("A/globalPV/hSecondsNumPVcontributors", "", kTH1D, {axisSeconds});
      histos.add("A/globalPV/hSecondsSumDcaR", "", kTH1D, {axisSeconds});
      histos.add("A/globalPV/hSecondsSumDcaZ", "", kTH1D, {axisSeconds});
      histos.add("A/globalPV/hSecondsSumPt", "", kTH1D, {axisSeconds});
      histos.add("A/globalPV/hSecondsNumClsIts", "", kTH2D, {axisSeconds, axisNclsITS});
      histos.add("A/globalPV/hSecondsChi2NClIts", "", kTH1D, {axisSeconds});
      histos.add("A/globalPV/hSecondsNumClsTpc", "", kTH2D, {axisSeconds, axisNclsTPC});
      histos.add("A/globalPV/hSecondsChi2NClTpc", "", kTH1D, {axisSeconds});
      histos.add("A/globalPV/hSecondsTpcFractionSharedCls", "", kTH2D, {axisSeconds, axisFraction});
      histos.add("A/globalPV/hSecondsDeDx", "", kTH1D, {axisSeconds});

      // ### C side
      // global tracks
      histos.add("C/global/hSecondsNumTracks", "", kTH1D, {axisSeconds});
      histos.add("C/global/hSecondsQoverPtSumDcaR", "", kTH2D, {axisSeconds, axisSparseQoverPt});
      histos.add("C/global/hSecondsQoverPtSumDcaZ", "", kTH2D, {axisSeconds, axisSparseQoverPt});
      histos.add("C/global/hSecondsSumDcaR", "", kTH1D, {axisSeconds});
      histos.add("C/global/hSecondsSumDcaZ", "", kTH1D, {axisSeconds});
      histos.add("C/global/hSecondsSumPt", "", kTH1D, {axisSeconds});
      histos.add("C/global/hSecondsNumClsIts", "", kTH2D, {axisSeconds, axisNclsITS});
      histos.add("C/global/hSecondsChi2NClIts", "", kTH1D, {axisSeconds});
      histos.add("C/global/hSecondsNumClsTpc", "", kTH2D, {axisSeconds, axisNclsTPC});
      histos.add("C/global/hSecondsChi2NClTpc", "", kTH1D, {axisSeconds});
      histos.add("C/global/hSecondsTpcFractionSharedCls", "", kTH2D, {axisSeconds, axisFraction});
      histos.add("C/global/hSecondsDeDx", "", kTH1D, {axisSeconds});

      // global && PV tracks
      histos.add("C/globalPV/hSecondsNumPVcontributors", "", kTH1D, {axisSeconds});
      histos.add("C/globalPV/hSecondsSumDcaR", "", kTH1D, {axisSeconds});
      histos.add("C/globalPV/hSecondsSumDcaZ", "", kTH1D, {axisSeconds});
      histos.add("C/globalPV/hSecondsSumPt", "", kTH1D, {axisSeconds});
      histos.add("C/globalPV/hSecondsNumClsIts", "", kTH2D, {axisSeconds, axisNclsITS});
      histos.add("C/globalPV/hSecondsChi2NClIts", "", kTH1D, {axisSeconds});
      histos.add("C/globalPV/hSecondsNumClsTpc", "", kTH2D, {axisSeconds, axisNclsTPC});
      histos.add("C/globalPV/hSecondsChi2NClTpc", "", kTH1D, {axisSeconds});
      histos.add("C/globalPV/hSecondsTpcFractionSharedCls", "", kTH2D, {axisSeconds, axisFraction});
      histos.add("C/globalPV/hSecondsDeDx", "", kTH1D, {axisSeconds});

      // phi holes vs time
      const AxisSpec axisPhi{64, 0, TMath::TwoPi(), "#varphi"}; // o2-linter: disable=external-pi (temporary fix)
      const AxisSpec axisEta{10, -0.8, 0.8, "#eta"};
      if (confFillPhiVsTimeHist == 2) {
        histos.add("hSecondsITSlayer0vsPhi", "", kTH2F, {axisSeconds, axisPhi});
        histos.add("hSecondsITSlayer1vsPhi", "", kTH2F, {axisSeconds, axisPhi});
        histos.add("hSecondsITSlayer2vsPhi", "", kTH2F, {axisSeconds, axisPhi});
        histos.add("hSecondsITSlayer3vsPhi", "", kTH2F, {axisSeconds, axisPhi});
        histos.add("hSecondsITSlayer4vsPhi", "", kTH2F, {axisSeconds, axisPhi});
        histos.add("hSecondsITSlayer5vsPhi", "", kTH2F, {axisSeconds, axisPhi});
        histos.add("hSecondsITSlayer6vsPhi", "", kTH2F, {axisSeconds, axisPhi});
      }
      if (confFillPhiVsTimeHist > 0) {
        histos.add("hSecondsITS7clsVsPhi", "", kTH2F, {axisSeconds, axisPhi});
        histos.add("hSecondsITSglobalVsPhi", "", kTH2F, {axisSeconds, axisPhi});
        histos.add("hSecondsITSTRDVsPhi", "", kTH2F, {axisSeconds, axisPhi});
        histos.add("hSecondsITSTOFVsPhi", "", kTH2F, {axisSeconds, axisPhi});
      }
      if (confFillEtaPhiVsTimeHist)
        histos.add("hSecondsITSglobalVsEtaPhi", "", kTH3F, {axisSeconds, axisEta, axisPhi});
    }

    // count TVX triggers per DF
    for (const auto& bc : bcs) {
      // auto bc = col.foundBC_as<BCsRun3>();
      int64_t ts = bc.timestamp();
      double secFromSOR = ts / 1000. - minSec;
      if (bc.selection_bit(kIsTriggerTVX)) {
        histos.fill(HIST("hSecondsBCsTVX"), secFromSOR);
      }
      // if (bc.selection_bit(kNoTimeFrameBorder)) {
      //   histos.fill(HIST("hSecondsBCsTFborder"), secFromSOR);
      // }
      if (bc.selection_bit(kIsTriggerTVX) && bc.selection_bit(kNoTimeFrameBorder)) {
        histos.fill(HIST("hSecondsBCsTVXandTFborder"), secFromSOR);
      }
    }

    // ### collision loop
    for (const auto& col : cols) {
      // check if a vertex is found in the UPC mode ITS ROF
      // flags from: https://github.com/AliceO2Group/AliceO2/blob/dev/DataFormats/Reconstruction/include/ReconstructionDataFormats/Vertex.h
      ushort flags = col.flags();
      bool isVertexUPC = flags & dataformats::Vertex<o2::dataformats::TimeStamp<int>>::Flags::UPCMode; // is vertex with UPC settings
      if (confTakeVerticesWithUPCsettings > 0) {                                                       // otherwise analyse all collisions
        if (confTakeVerticesWithUPCsettings == 1 && isVertexUPC)                                       // reject vertices with UPC settings
          continue;
        if (confTakeVerticesWithUPCsettings == 2 && !isVertexUPC) // we want to select vertices with UPC settings --> reject vertices reconstructed with "normal" settings
          continue;
        // LOGP(info, "flags={} nTracks = {}", flags, tracks.size());
      }

      auto bc = col.foundBC_as<BCsRun3>();
      int64_t ts = bc.timestamp();
      double secFromSOR = ts / 1000. - minSec;

      histos.fill(HIST("hSecondsCollisionsBeforeAllCuts"), secFromSOR);
      if (col.selection_bit(kIsTriggerTVX))
        histos.fill(HIST("hSecondsCollisionsNoVzInTVX"), secFromSOR);
      if (col.selection_bit(kNoTimeFrameBorder))
        histos.fill(HIST("hSecondsCollisionsNoVzNoTFborder"), secFromSOR);
      if (col.selection_bit(kIsTriggerTVX) && col.selection_bit(kNoTimeFrameBorder))
        histos.fill(HIST("hSecondsCollisionsNoVzInTVXandNoTFborder"), secFromSOR);

      if (std::fabs(col.posZ()) > 10)
        continue;

      histos.fill(HIST("hSecondsUPCverticesBeforeSel8"), secFromSOR, isVertexUPC ? 1 : 0);

      histos.fill(HIST("hSecondsEventSelBits"), secFromSOR, enCollisionsAll);
      histos.fill(HIST("hSecondsEventSelBits"), secFromSOR, enIsTriggerTVX, col.selection_bit(kIsTriggerTVX));
      histos.fill(HIST("hSecondsEventSelBits"), secFromSOR, enNoTimeFrameBorder, col.selection_bit(kNoTimeFrameBorder));
      histos.fill(HIST("hSecondsEventSelBits"), secFromSOR, enNoITSROFrameBorder, col.selection_bit(kNoITSROFrameBorder));

      // for QA:
      uint64_t globalBC = bc.globalBC();
      int64_t bcInTF = (globalBC - bcSOR) % nBCsPerTF;

      histos.fill(HIST("hNcolVsBcInTF"), bcInTF);
      if (!col.selection_bit(kNoTimeFrameBorder))
        histos.fill(HIST("hNcolVsBcInTFantiBorderCut"), bcInTF);

      // sel8 selection:
      if (!col.sel8())
        continue;

      histos.fill(HIST("hSecondsUPCvertices"), secFromSOR, isVertexUPC ? 1 : 0);
      histos.fill(HIST("hSecondsCollisions"), secFromSOR);
      histos.fill(HIST("hSecondsVz"), secFromSOR, col.posZ());
      histos.fill(HIST("hSecondsFT0Camlp"), secFromSOR, bc.foundFT0().sumAmpC());
      histos.fill(HIST("hSecondsFT0CamlpByColMult"), secFromSOR, col.multFT0C());
      histos.fill(HIST("hSecondsFT0AamlpByColMult"), secFromSOR, col.multFT0A());
      histos.fill(HIST("hSecondsV0Aamlp"), secFromSOR, col.multFV0A());

      histos.fill(HIST("hSecondsOccupancyByTracks"), secFromSOR, col.trackOccupancyInTimeRange());
      histos.fill(HIST("hSecondsOccupancyByFT0C"), secFromSOR, col.ft0cOccupancyInTimeRange());

      histos.fill(HIST("hSecondsEventSelBits"), secFromSOR, enCollisionsSel8);
      histos.fill(HIST("hSecondsEventSelBits"), secFromSOR, enNoSameBunchPileup, col.selection_bit(kNoSameBunchPileup));
      histos.fill(HIST("hSecondsEventSelBits"), secFromSOR, enIsGoodZvtxFT0vsPV, col.selection_bit(kIsGoodZvtxFT0vsPV));
      histos.fill(HIST("hSecondsEventSelBits"), secFromSOR, enIsVertexITSTPC, col.selection_bit(kIsVertexITSTPC));
      histos.fill(HIST("hSecondsEventSelBits"), secFromSOR, enIsVertexTOFmatched, col.selection_bit(kIsVertexTOFmatched));
      histos.fill(HIST("hSecondsEventSelBits"), secFromSOR, enIsVertexTRDmatched, col.selection_bit(kIsVertexTRDmatched));
      histos.fill(HIST("hSecondsEventSelBits"), secFromSOR, enNoCollInTimeRangeNarrow, col.selection_bit(kNoCollInTimeRangeNarrow));
      histos.fill(HIST("hSecondsEventSelBits"), secFromSOR, enNoCollInTimeRangeStrict, col.selection_bit(kNoCollInTimeRangeStrict));
      histos.fill(HIST("hSecondsEventSelBits"), secFromSOR, enNoCollInTimeRangeStandard, col.selection_bit(kNoCollInTimeRangeStandard));
      histos.fill(HIST("hSecondsEventSelBits"), secFromSOR, enNoCollInRofStrict, col.selection_bit(kNoCollInRofStrict));
      histos.fill(HIST("hSecondsEventSelBits"), secFromSOR, enNoCollInRofStandard, col.selection_bit(kNoCollInRofStandard));
      histos.fill(HIST("hSecondsEventSelBits"), secFromSOR, enNoHighMultCollInPrevRof, col.selection_bit(kNoHighMultCollInPrevRof));
      histos.fill(HIST("hSecondsEventSelBits"), secFromSOR, enIsGoodITSLayer3, col.selection_bit(kIsGoodITSLayer3));
      histos.fill(HIST("hSecondsEventSelBits"), secFromSOR, enIsGoodITSLayer0123, col.selection_bit(kIsGoodITSLayer0123));
      histos.fill(HIST("hSecondsEventSelBits"), secFromSOR, enIsGoodITSLayersAll, col.selection_bit(kIsGoodITSLayersAll));

      // occupancy selection combinations
      float occupByTracks = col.trackOccupancyInTimeRange();

      bool isLowOccupStd = col.selection_bit(kNoCollInTimeRangeStandard) && col.selection_bit(kNoCollInRofStandard);
      histos.fill(HIST("hSecondsEventSelBits"), secFromSOR, enIsLowOccupStd, isLowOccupStd);

      bool isLowOccupStdAlsoInPrevRof = isLowOccupStd && col.selection_bit(kNoHighMultCollInPrevRof);
      histos.fill(HIST("hSecondsEventSelBits"), secFromSOR, enIsLowOccupStdAlsoInPrevRof, isLowOccupStdAlsoInPrevRof);

      bool isLowOccupStdCut500 = isLowOccupStd && occupByTracks >= 0 && occupByTracks < 500;
      histos.fill(HIST("hSecondsEventSelBits"), secFromSOR, enIsLowOccupStdCut500, isLowOccupStdCut500);

      bool isLowOccupStdCut2000 = isLowOccupStd && occupByTracks >= 0 && occupByTracks < 2000;
      histos.fill(HIST("hSecondsEventSelBits"), secFromSOR, enIsLowOccupStdCut2000, isLowOccupStdCut2000);

      bool isLowOccupStdCut4000 = isLowOccupStd && occupByTracks >= 0 && occupByTracks < 4000;
      histos.fill(HIST("hSecondsEventSelBits"), secFromSOR, enIsLowOccupStdCut4000, isLowOccupStdCut4000);

      bool isLowOccupStdAlsoInPrevRofCut2000noDeadStaves = isLowOccupStdCut2000 && col.selection_bit(kNoHighMultCollInPrevRof) && col.selection_bit(kIsGoodITSLayersAll);
      histos.fill(HIST("hSecondsEventSelBits"), secFromSOR, enIsLowOccupStdAlsoInPrevRofCut2000noDeadStaves, isLowOccupStdAlsoInPrevRofCut2000noDeadStaves);

      double hadronicRate = mRateFetcher.fetch(ccdb.service, ts, runNumber, "ZNC hadronic") * 1.e-3; // kHz
      histos.fill(HIST("hSecondsIR"), secFromSOR, hadronicRate);

      // checking mShape flags in time:
      auto mShapeTree = ccdb->getForTimeStamp<TTree>("TPC/Calib/MShapePotential", ts);
      mshape.setFromTree(*mShapeTree);
      bool isMshape = !mshape.getBoundaryPotential(ts).mPotential.empty();

      // ##### track loop
      for (const auto& track : tracks) {
        // if (!track.hasTPC() || !track.hasITS())
        //   continue;
        if (std::fabs(track.eta()) > 0.8 || std::fabs(track.pt()) < 0.2)
          continue;

        double dcaR = track.dcaXY();
        double dcaZ = track.dcaZ();
        LOGP(debug, "dcaR = {} dcaZ = {}", dcaR, dcaZ);
        histos.fill(HIST("allTracks/hDcaR"), dcaR);
        histos.fill(HIST("allTracks/hDcaZ"), dcaZ);

        // now DCA cuts:
        if (std::fabs(dcaR) > 1. || std::fabs(dcaZ) > 1.)
          continue;

        histos.fill(HIST("allTracks/hSecondsTracks"), secFromSOR);

        histos.fill(HIST("allTracks/hDcaRafterCuts"), dcaR);
        histos.fill(HIST("allTracks/hDcaZafterCuts"), dcaZ);

        double qpt = track.signed1Pt();
        histos.fill(HIST("allTracks/hQoverPt"), qpt);
        histos.fill(HIST("allTracks/hQoverPtDcaR"), qpt, dcaR);
        histos.fill(HIST("allTracks/hQoverPtDcaZ"), qpt, dcaZ);

        // now consider only abs values for DCAs:
        double dcaRabs = std::fabs(dcaR);
        double dcaZabs = std::fabs(dcaZ);

        histos.fill(HIST("allTracks/hSecondsSumDcaR"), secFromSOR, dcaRabs);
        histos.fill(HIST("allTracks/hSecondsSumDcaZ"), secFromSOR, dcaZabs);
        histos.fill(HIST("allTracks/hSecondsSumPt"), secFromSOR, track.pt());
        histos.fill(HIST("allTracks/hSecondsQoverPtSumDcaR"), secFromSOR, qpt, dcaRabs);
        histos.fill(HIST("allTracks/hSecondsQoverPtSumDcaZ"), secFromSOR, qpt, dcaZabs);
        histos.fill(HIST("allTracks/hSecondsNumClsIts"), secFromSOR, track.itsNCls());
        histos.fill(HIST("allTracks/hSecondsChi2NClIts"), secFromSOR, track.itsChi2NCl());
        if (isMshape) {
          histos.fill(HIST("allTracks/hSecondsTracksMshape"), secFromSOR);
        }

        // ### PV contributors
        if (track.isPVContributor()) {
          histos.fill(HIST("PVcontrib/hDcaRafterCuts"), dcaR);
          histos.fill(HIST("PVcontrib/hDcaZafterCuts"), dcaZ);

          histos.fill(HIST("PVcontrib/hSecondsTracks"), secFromSOR);
          histos.fill(HIST("PVcontrib/hSecondsSumDcaR"), secFromSOR, dcaRabs);
          histos.fill(HIST("PVcontrib/hSecondsSumDcaZ"), secFromSOR, dcaZabs);
          histos.fill(HIST("PVcontrib/hSecondsSumPt"), secFromSOR, track.pt());
          // histos.fill(HIST("PVcontrib/hSecondsQoverPtSumDcaR"), secFromSOR, qpt, dcaRabs);
          // histos.fill(HIST("PVcontrib/hSecondsQoverPtSumDcaZ"), secFromSOR, qpt, dcaZabs);
          histos.fill(HIST("PVcontrib/hSecondsNumClsIts"), secFromSOR, track.itsNCls());
          histos.fill(HIST("PVcontrib/hSecondsChi2NClIts"), secFromSOR, track.itsChi2NCl());
        }

        // ### global tracks
        float dedx = track.tpcSignal();
        if (track.isGlobalTrack()) { // A side
          if (track.tgl() > 0.) {
            histos.fill(HIST("A/global/hDcaRafterCuts"), dcaR);
            histos.fill(HIST("A/global/hDcaZafterCuts"), dcaZ);

            histos.fill(HIST("A/global/hSecondsNumTracks"), secFromSOR);
            histos.fill(HIST("A/global/hSecondsQoverPtSumDcaR"), secFromSOR, qpt, dcaRabs);
            histos.fill(HIST("A/global/hSecondsQoverPtSumDcaZ"), secFromSOR, qpt, dcaZabs);
            histos.fill(HIST("A/global/hSecondsSumDcaR"), secFromSOR, dcaRabs);
            histos.fill(HIST("A/global/hSecondsSumDcaZ"), secFromSOR, dcaZabs);
            histos.fill(HIST("A/global/hSecondsSumPt"), secFromSOR, track.pt());
            histos.fill(HIST("A/global/hSecondsNumClsIts"), secFromSOR, track.itsNCls());
            histos.fill(HIST("A/global/hSecondsChi2NClIts"), secFromSOR, track.itsChi2NCl());
            histos.fill(HIST("A/global/hSecondsNumClsTpc"), secFromSOR, track.tpcNClsFound());
            histos.fill(HIST("A/global/hSecondsChi2NClTpc"), secFromSOR, track.tpcChi2NCl());
            if (track.tpcNClsFound() >= confCutOnNtpcClsForSharedFractAndDeDxCalc) {
              histos.fill(HIST("A/global/hSecondsTpcFractionSharedCls"), secFromSOR, track.tpcFractionSharedCls());
              if (dedx < 1.e4) // protection from weird values
                histos.fill(HIST("A/global/hSecondsDeDx"), secFromSOR, dedx);
            }

            if (track.isPVContributor()) {
              histos.fill(HIST("A/globalPV/hDcaRafterCuts"), dcaR);
              histos.fill(HIST("A/globalPV/hDcaZafterCuts"), dcaZ);

              histos.fill(HIST("A/globalPV/hSecondsNumPVcontributors"), secFromSOR);
              histos.fill(HIST("A/globalPV/hSecondsSumDcaR"), secFromSOR, dcaRabs);
              histos.fill(HIST("A/globalPV/hSecondsSumDcaZ"), secFromSOR, dcaZabs);
              histos.fill(HIST("A/globalPV/hSecondsSumPt"), secFromSOR, track.pt());
              histos.fill(HIST("A/globalPV/hSecondsNumClsIts"), secFromSOR, track.itsNCls());
              histos.fill(HIST("A/globalPV/hSecondsChi2NClIts"), secFromSOR, track.itsChi2NCl());
              histos.fill(HIST("A/globalPV/hSecondsNumClsTpc"), secFromSOR, track.tpcNClsFound());
              histos.fill(HIST("A/globalPV/hSecondsChi2NClTpc"), secFromSOR, track.tpcChi2NCl());
              if (track.tpcNClsFound() >= confCutOnNtpcClsForSharedFractAndDeDxCalc) {
                histos.fill(HIST("A/globalPV/hSecondsTpcFractionSharedCls"), secFromSOR, track.tpcFractionSharedCls());
                if (dedx < 1.e4) // protection from weird values
                  histos.fill(HIST("A/globalPV/hSecondsDeDx"), secFromSOR, dedx);
              }
            }
          } else { // C side
            histos.fill(HIST("C/global/hDcaRafterCuts"), dcaR);
            histos.fill(HIST("C/global/hDcaZafterCuts"), dcaZ);

            histos.fill(HIST("C/global/hSecondsNumTracks"), secFromSOR);
            histos.fill(HIST("C/global/hSecondsQoverPtSumDcaR"), secFromSOR, qpt, dcaRabs);
            histos.fill(HIST("C/global/hSecondsQoverPtSumDcaZ"), secFromSOR, qpt, dcaZabs);
            histos.fill(HIST("C/global/hSecondsSumDcaR"), secFromSOR, dcaRabs);
            histos.fill(HIST("C/global/hSecondsSumDcaZ"), secFromSOR, dcaZabs);
            histos.fill(HIST("C/global/hSecondsSumPt"), secFromSOR, track.pt());
            histos.fill(HIST("C/global/hSecondsNumClsIts"), secFromSOR, track.itsNCls());
            histos.fill(HIST("C/global/hSecondsChi2NClIts"), secFromSOR, track.itsChi2NCl());
            histos.fill(HIST("C/global/hSecondsNumClsTpc"), secFromSOR, track.tpcNClsFound());
            histos.fill(HIST("C/global/hSecondsChi2NClTpc"), secFromSOR, track.tpcChi2NCl());
            if (track.tpcNClsFound() >= confCutOnNtpcClsForSharedFractAndDeDxCalc) {
              histos.fill(HIST("C/global/hSecondsTpcFractionSharedCls"), secFromSOR, track.tpcFractionSharedCls());
              if (dedx < 1.e4) // protection from weird values
                histos.fill(HIST("C/global/hSecondsDeDx"), secFromSOR, dedx);
            }

            if (track.isPVContributor()) {
              histos.fill(HIST("C/globalPV/hDcaRafterCuts"), dcaR);
              histos.fill(HIST("C/globalPV/hDcaZafterCuts"), dcaZ);

              histos.fill(HIST("C/globalPV/hSecondsNumPVcontributors"), secFromSOR);
              histos.fill(HIST("C/globalPV/hSecondsSumDcaR"), secFromSOR, dcaRabs);
              histos.fill(HIST("C/globalPV/hSecondsSumDcaZ"), secFromSOR, dcaZabs);
              histos.fill(HIST("C/globalPV/hSecondsSumPt"), secFromSOR, track.pt());
              histos.fill(HIST("C/globalPV/hSecondsNumClsIts"), secFromSOR, track.itsNCls());
              histos.fill(HIST("C/globalPV/hSecondsChi2NClIts"), secFromSOR, track.itsChi2NCl());
              histos.fill(HIST("C/globalPV/hSecondsNumClsTpc"), secFromSOR, track.tpcNClsFound());
              histos.fill(HIST("C/globalPV/hSecondsChi2NClTpc"), secFromSOR, track.tpcChi2NCl());
              if (track.tpcNClsFound() >= confCutOnNtpcClsForSharedFractAndDeDxCalc) {
                histos.fill(HIST("C/globalPV/hSecondsTpcFractionSharedCls"), secFromSOR, track.tpcFractionSharedCls());
                if (dedx < 1.e4) // protection from weird values
                  histos.fill(HIST("C/globalPV/hSecondsDeDx"), secFromSOR, dedx);
              }
            }
          }
        } // end of global tracks

        // study ITS cluster pattern vs phi vs time (pt>1 GeV/c cut selects straight tracks)
        if (track.isPVContributor() && track.pt() > 1) {
          // layer-by-layer check
          if (confFillPhiVsTimeHist == 2) {
            if (track.itsClusterMap() & (1 << 0))
              histos.fill(HIST("hSecondsITSlayer0vsPhi"), secFromSOR, track.phi());
            if (track.itsClusterMap() & (1 << 1))
              histos.fill(HIST("hSecondsITSlayer1vsPhi"), secFromSOR, track.phi());
            if (track.itsClusterMap() & (1 << 2))
              histos.fill(HIST("hSecondsITSlayer2vsPhi"), secFromSOR, track.phi());
            if (track.itsClusterMap() & (1 << 3))
              histos.fill(HIST("hSecondsITSlayer3vsPhi"), secFromSOR, track.phi());
            if (track.itsClusterMap() & (1 << 4))
              histos.fill(HIST("hSecondsITSlayer4vsPhi"), secFromSOR, track.phi());
            if (track.itsClusterMap() & (1 << 5))
              histos.fill(HIST("hSecondsITSlayer5vsPhi"), secFromSOR, track.phi());
            if (track.itsClusterMap() & (1 << 6))
              histos.fill(HIST("hSecondsITSlayer6vsPhi"), secFromSOR, track.phi());
          }
          // tracks with conditions
          if (confFillPhiVsTimeHist > 0) {
            if (track.itsNCls() == 7)
              histos.fill(HIST("hSecondsITS7clsVsPhi"), secFromSOR, track.phi());
            if (track.isGlobalTrack())
              histos.fill(HIST("hSecondsITSglobalVsPhi"), secFromSOR, track.phi());
            if (track.hasTRD())
              histos.fill(HIST("hSecondsITSTRDVsPhi"), secFromSOR, track.phi());
            if (track.hasTOF())
              histos.fill(HIST("hSecondsITSTOFVsPhi"), secFromSOR, track.phi());
          }
          // eta-phi histogram for global tracks
          if (confFillEtaPhiVsTimeHist && track.isGlobalTrack()) {
            histos.fill(HIST("hSecondsITSglobalVsEtaPhi"), secFromSOR, track.eta(), track.phi());
          }
        }
      }
    }
  } // end of collision loop
  PROCESS_SWITCH(TimeDependentQaTask, processRun3, "Process Run3 QA vs time", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<TimeDependentQaTask>(cfgc)};
}
