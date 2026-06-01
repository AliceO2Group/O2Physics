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

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/RCTSelectionFlags.h"
#include "Common/CCDB/ctpRateFetcher.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/LHCConstants.h>
#include <CommonDataFormat/TimeStamp.h>
#include <DataFormatsParameters/AggregatedRunInfo.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/Vertex.h>
#include <TPCCalibration/TPCMShapeCorrection.h>

#include <TAxis.h>
#include <TH2.h>
#include <TMath.h>
#include <TTree.h>

#include <sys/types.h>

#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::aod::evsel;
using namespace o2::aod::rctsel;

using ColEvSels = soa::Join<aod::Collisions, aod::EvSels, aod::Mults>;
using BCsRun3 = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;
using BarrelTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>;

const AxisSpec axisQoverPt{100, -1., 1., "q/p_{T}, 1/GeV"};
const AxisSpec axisDcaR{1000, -5., 5., "DCA_{r}, cm"};
const AxisSpec axisDcaZ{1000, -5., 5., "DCA_{z}, cm"};
const AxisSpec axisSparseQoverPt{20, -1., 1., "q/p_{T}, 1/GeV"};
const AxisSpec axisSparseDcaR{100, -1., 1., "DCA_{r}, cm"};
const AxisSpec axisSparseDcaZ{100, -1., 1., "DCA_{z}, cm"};

struct TimeDependentQaTask {
  Configurable<float> confTimeBinWidthInSec{"TimeBinWidthInSec", 0.5, "Width of time bins in seconds"};                                                                          // o2-linter: disable=name/configurable (temporary fix)
  Configurable<float> confTimeWiderBinFactor{"TimeWideBinFactor", 4, "Factor for wider time bins for some 2D histograms"};                                                       // o2-linter: disable=name/configurable (temporary fix)
  Configurable<float> confTimeMuchWiderBinFactor{"TimeMuchWiderBinFactor", 20, "Factor for even wider time bins for some 2D histograms"};                                        // o2-linter: disable=name/configurable (temporary fix)
  Configurable<float> confTimeMuchMuchWiderBinFactor{"TimeMuchMuchWiderBinFactor", 120, "Factor for super wide time bins for some 2D histograms"};                               // o2-linter: disable=name/configurable (temporary fix)
  Configurable<int> confTakeVerticesWithUPCsettings{"ConsiderVerticesWithUPCsettings", 0, "Take vertices: 0 - all , 1 - only without UPC settings, 2 - only with UPC settings"}; // o2-linter: disable=name/configurable (temporary fix)
  Configurable<int> confFlagFillPhiVsTimeHist{"FlagFillPhiVsTimeHist", 2, "0 - don't fill , 1 - fill only for global/7cls/TRD/TOF tracks, 2 - fill also layer-by-layer"};        // o2-linter: disable=name/configurable (temporary fix)
  Configurable<int> confFlagFillEtaPhiVsTimeHist{"FlagFillEtaPhiVsTimeHist", 0, "0 - don't fill , 1 - fill"};                                                                    // o2-linter: disable=name/configurable (temporary fix)
  Configurable<float> confCutOnNtpcClsForSharedFractAndDeDxCalc{"CutOnNtpcClsForSharedFractAndDeDxCalc", 70, ""};                                                                // o2-linter: disable=name/configurable (temporary fix)
  Configurable<int> confFlagCheckMshape{"FlagCheckMshape", 0, "0 - don't check , 1 - check"};                                                                                    // o2-linter: disable=name/configurable (temporary fix)
  Configurable<int> confFlagCheckQoverPtHist{"FlagCheckQoverPtHist", 1, "0 - don't check , 1 - check"};                                                                          // o2-linter: disable=name/configurable (temporary fix)

  // for O-O and Ne-Ne run
  Configurable<int> confIncludeMultDistrVsTimeHistos{"IncludeMultDistrVsTimeHistos", 0, ""};                                                                    // o2-linter: disable=name/configurable (temporary fix)
  Configurable<float> confMaxNtracksForTimeDepDistributions{"MaxNtracksForTimeDepDistributions", 800, ""};                                                      // o2-linter: disable=name/configurable (temporary fix)
  Configurable<float> confMaxZNACenergyForTimeDepDistributions{"MaxZNACenergyForTimeDepDistributions", 4000, ""};                                               // o2-linter: disable=name/configurable (temporary fix)
  Configurable<float> confMaxT0ACamplForTimeDepDistributions{"MaxT0ACamplForTimeDepDistributions", 25000, ""};                                                  // o2-linter: disable=name/configurable (temporary fix)
  Configurable<float> confMaxV0AamplForTimeDepDistributions{"MaxV0AamplForTimeDepDistributions", 40000, ""};                                                    // o2-linter: disable=name/configurable (temporary fix)
  Configurable<int> confFlagUseGlobalTracksForTimeDepDistributions{"FlagUseGlobalTracksForTimeDepDistributions", 0, "0 - PV contributors , 1 - global tracks"}; // o2-linter: disable=name/configurable (temporary fix)
  Configurable<int> confFlagUseGoodZvtxFT0vsPVForTimeDepDistributions{"FlagUseGoodZvtxFT0vsPVForTimeDepDistributions", 0, "0 - no , 1 - yes"};                  // o2-linter: disable=name/configurable (temporary fix)

  enum EvSelBitsToMonitor {
    enCollisionsAll = 0,
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

  enum RctCombFlagsToMonitor {
    enCBT = kNRCTSelectionFlags,
    enCBT_hadronPID,
    enCBT_electronPID,
    enCBT_calo,
    enCBT_muon,
    enCBT_muon_glo,
    enNumRctFlagsTotal, // counter
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

  // RCT flag combinations: checkers (based on presentation https://indico.cern.ch/event/1513866/#18-how-to-use-the-rct-flags-at)
  RCTFlagsChecker rctCheckerCBT{"CBT"};                         // o2-linter: disable=name/function-variable (temporary fix)
  RCTFlagsChecker rctCheckerCBT_hadronPID{"CBT_hadronPID"};     // o2-linter: disable=name/function-variable (temporary fix)
  RCTFlagsChecker rctCheckerCBT_electronPID{"CBT_electronPID"}; // o2-linter: disable=name/function-variable (temporary fix)
  RCTFlagsChecker rctCheckerCBT_calo{"CBT_calo"};               // o2-linter: disable=name/function-variable (temporary fix)
  RCTFlagsChecker rctCheckerCBT_muon{"CBT_muon"};               // o2-linter: disable=name/function-variable (temporary fix)
  RCTFlagsChecker rctCheckerCBT_muon_glo{"CBT_muon_glo"};       // o2-linter: disable=name/function-variable (temporary fix)

  TAxis* axRctFlags;

  void init(InitContext&)
  {
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    histos.add("allTracks/hQoverPt", "", kTH1F, {axisQoverPt});
    if (confFlagCheckQoverPtHist) {
      histos.add("allTracks/hQoverPtDcaR", "", kTH2F, {axisSparseQoverPt, axisSparseDcaR});
      histos.add("allTracks/hQoverPtDcaZ", "", kTH2F, {axisSparseQoverPt, axisSparseDcaZ});
    }
    histos.add("allTracks/hDcaR", "", kTH1F, {axisDcaR});
    histos.add("allTracks/hDcaZ", "", kTH1F, {axisDcaZ});
    histos.add("allTracks/hDcaRafterCuts", "", kTH1F, {axisDcaR});
    histos.add("allTracks/hDcaZafterCuts", "", kTH1F, {axisDcaZ});

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
  }

  Preslice<BarrelTracks> perCollision = aod::track::collisionId;

  void processRun3(
    ColEvSels const& cols,
    BarrelTracks const& tracks,
    BCsRun3 const& bcs,
    aod::Zdcs const&,
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
      int nTimeWideBins = static_cast<int>((maxSec - minSec) / confTimeBinWidthInSec / confTimeWiderBinFactor);
      int nTimeVeryWideBins = static_cast<int>((maxSec - minSec) / confTimeBinWidthInSec / confTimeMuchWiderBinFactor);
      int nTimeSuperWideBins = static_cast<int>((maxSec - minSec) / confTimeBinWidthInSec / confTimeMuchMuchWiderBinFactor);
      double timeInterval = nTimeBins * confTimeBinWidthInSec;

      const AxisSpec axisSeconds{nTimeBins, 0, timeInterval, "seconds"};
      const AxisSpec axisSecondsWideBins{nTimeWideBins, 0, timeInterval, "seconds"};
      const AxisSpec axisSecondsVeryWideBins{nTimeVeryWideBins, 0, timeInterval, "seconds"};
      const AxisSpec axisSecondsSuperWideBins{nTimeSuperWideBins, 0, timeInterval, "seconds"};

      histos.add("hSecondsBCsTVX", "", kTH1D, {axisSeconds});
      histos.add("hSecondsBCsTVXandTFborderCuts", "", kTH1D, {axisSeconds});

      histos.add("hSecondsCollisionsBeforeAllCuts", "", kTH1D, {axisSeconds});
      histos.add("hSecondsCollisionsTVXNoVzCut", "", kTH1D, {axisSeconds});
      histos.add("hSecondsCollisionsTFborderCutNoVzCut", "", kTH1D, {axisSeconds});
      histos.add("hSecondsCollisionsTVXTFborderCutNoVzCut", "", kTH1D, {axisSeconds});

      histos.add("hSecondsCollisions", "", kTH1D, {axisSeconds});
      histos.add("hSecondsCollisionsNoPileup", "", kTH1D, {axisSeconds});
      histos.add("hSecondsIR", "", kTH1D, {axisSeconds});
      histos.add("hSecondsVz", "", kTH1D, {axisSeconds});
      histos.add("hSecondsFT0Camlp", "", kTH1D, {axisSeconds});
      histos.add("hSecondsFT0CamlpByColMult", "", kTH1D, {axisSeconds});
      histos.add("hSecondsFT0AamlpByColMult", "", kTH1D, {axisSeconds});
      histos.add("hSecondsV0Aamlp", "", kTH1D, {axisSeconds});
      histos.add("hSecondsOccupancyByTracks", "", kTH1D, {axisSeconds});
      histos.add("hSecondsOccupancyByFT0C", "", kTH1D, {axisSeconds});

      // QA for UPC settings
      histos.add("hSecondsUPCverticesBeforeAllCuts", "", kTH2F, {axisSeconds, {2, -0.5, 1.5, "Is vertex with UPC settings"}});
      histos.add("hSecondsUPCverticesBeforeSel8", "", kTH2F, {axisSeconds, {2, -0.5, 1.5, "Is vertex with UPC settings after |vZ|<10 cut"}});
      histos.add("hSecondsUPCvertices", "", kTH2F, {axisSeconds, {2, -0.5, 1.5, "Is vertex with UPC settings after |vZ|<10 and sel8 cuts"}});

      const int32_t nBCsPerOrbit = o2::constants::lhc::LHCMaxBunches;
      const AxisSpec axisBCs{nBCsPerOrbit, 0., static_cast<double>(nBCsPerOrbit), ""};
      histos.add("hSecondsBCsMap", "", kTH2F, {axisSecondsSuperWideBins, axisBCs});

      // shapes of distributions (added for the O-O run monitoring)
      if (confIncludeMultDistrVsTimeHistos) {
        int maxNtracks = confMaxNtracksForTimeDepDistributions;
        float maxZNACenergyForTimeDepDistributions = confMaxZNACenergyForTimeDepDistributions;
        float maxT0ACamplForTimeDepDistributions = confMaxT0ACamplForTimeDepDistributions;
        float maxV0AamplForTimeDepDistributions = confMaxV0AamplForTimeDepDistributions;
        histos.add("multDistributions/hSecondsDistrPVtracks", "", kTH2D, {axisSecondsVeryWideBins, {maxNtracks, -0.5, maxNtracks - 0.5, "n PV tracks"}});
        histos.add("multDistributions/hSecondsDistrT0A", "", kTH2D, {axisSecondsVeryWideBins, {250, 0, maxT0ACamplForTimeDepDistributions, "T0A ampl"}});
        histos.add("multDistributions/hSecondsDistrT0C", "", kTH2D, {axisSecondsVeryWideBins, {250, 0, maxT0ACamplForTimeDepDistributions, "T0C ampl"}});
        histos.add("multDistributions/hSecondsDistrV0A", "", kTH2D, {axisSecondsVeryWideBins, {400, 0, maxV0AamplForTimeDepDistributions, "V0A ampl"}});
        histos.add("multDistributions/hSecondsDistrZNA", "", kTH2D, {axisSecondsVeryWideBins, {320, 0, maxZNACenergyForTimeDepDistributions, "ZNA ampl"}});
        histos.add("multDistributions/hSecondsDistrZNC", "", kTH2D, {axisSecondsVeryWideBins, {320, 0, maxZNACenergyForTimeDepDistributions, "ZNC ampl"}});
        histos.add("multDistributions/hSecondsDistrZNACdiff", "", kTH2D, {axisSecondsVeryWideBins, {600, -maxZNACenergyForTimeDepDistributions, maxZNACenergyForTimeDepDistributions, "ZN A-C diff"}});
        histos.add("multDistributions/hSecondsDistrZNACdiffNorm", "", kTH2D, {axisSecondsVeryWideBins, {200, -1., 1., "ZN A-C diff"}});
        histos.add("multDistributions/hSecondsDistrZNAampl", "", kTH2D, {axisSecondsVeryWideBins, {320, 0, maxZNACenergyForTimeDepDistributions, "ZNA ampl"}});
        histos.add("multDistributions/hSecondsDistrZNCampl", "", kTH2D, {axisSecondsVeryWideBins, {320, 0, maxZNACenergyForTimeDepDistributions, "ZNC ampl"}});
        histos.add("multDistributions/hSecondsDistrZNACdiffAmpl", "", kTH2D, {axisSecondsVeryWideBins, {200, -1., 1., "ZN A-C diff"}});
        histos.add("multDistributions/hSecondsDistrZNACdiffNormAmpl", "", kTH2D, {axisSecondsVeryWideBins, {200, -1., 1., "ZN A-C diff"}});

        // 2D vs time
        histos.add("multDistributions/hSecondsNPVvsZNAampl", "", kTH3F, {axisSecondsSuperWideBins, {200, 0, 500, "n PV tracks"}, {200, 0, maxZNACenergyForTimeDepDistributions, "ZNA ampl"}});
        histos.add("multDistributions/hSecondsNPVvsZNCampl", "", kTH3F, {axisSecondsSuperWideBins, {200, 0, 500, "n PV tracks"}, {200, 0, maxZNACenergyForTimeDepDistributions, "ZNC ampl"}});

        // now histos after kNoSameBunchPileup cut:
        histos.add("multDistributionsNoPileup/hSecondsDistrPVtracks", "", kTH2D, {axisSecondsVeryWideBins, {maxNtracks, -0.5, maxNtracks - 0.5, "n PV tracks"}});
        histos.add("multDistributionsNoPileup/hSecondsDistrT0A", "", kTH2D, {axisSecondsVeryWideBins, {250, 0, maxT0ACamplForTimeDepDistributions, "T0A ampl"}});
        histos.add("multDistributionsNoPileup/hSecondsDistrT0C", "", kTH2D, {axisSecondsVeryWideBins, {250, 0, maxT0ACamplForTimeDepDistributions, "T0C ampl"}});
        histos.add("multDistributionsNoPileup/hSecondsDistrV0A", "", kTH2D, {axisSecondsVeryWideBins, {400, 0, maxV0AamplForTimeDepDistributions, "V0A ampl"}});
        histos.add("multDistributionsNoPileup/hSecondsDistrZNA", "", kTH2D, {axisSecondsVeryWideBins, {320, 0, maxZNACenergyForTimeDepDistributions, "ZNA ampl"}});
        histos.add("multDistributionsNoPileup/hSecondsDistrZNC", "", kTH2D, {axisSecondsVeryWideBins, {320, 0, maxZNACenergyForTimeDepDistributions, "ZNC ampl"}});
        histos.add("multDistributionsNoPileup/hSecondsDistrZNACdiff", "", kTH2D, {axisSecondsVeryWideBins, {600, -maxZNACenergyForTimeDepDistributions, maxZNACenergyForTimeDepDistributions, "ZN A-C diff"}});
        histos.add("multDistributionsNoPileup/hSecondsDistrZNACdiffNorm", "", kTH2D, {axisSecondsVeryWideBins, {200, -1., 1., "ZN A-C diff"}});
        histos.add("multDistributionsNoPileup/hSecondsDistrZNAampl", "", kTH2D, {axisSecondsVeryWideBins, {320, 0, maxZNACenergyForTimeDepDistributions, "ZNA ampl"}});
        histos.add("multDistributionsNoPileup/hSecondsDistrZNCampl", "", kTH2D, {axisSecondsVeryWideBins, {320, 0, maxZNACenergyForTimeDepDistributions, "ZNC ampl"}});
        histos.add("multDistributionsNoPileup/hSecondsDistrZNACdiffAmpl", "", kTH2D, {axisSecondsVeryWideBins, {200, -1., 1., "ZN A-C diff"}});
        histos.add("multDistributionsNoPileup/hSecondsDistrZNACdiffNormAmpl", "", kTH2D, {axisSecondsVeryWideBins, {200, -1., 1., "ZN A-C diff"}});

        // 2D vs time
        histos.add("multDistributionsNoPileup/hSecondsNPVvsZNAampl", "", kTH3F, {axisSecondsSuperWideBins, {200, 0, 500, "n PV tracks"}, {200, 0, maxZNACenergyForTimeDepDistributions, "ZNA ampl"}});
        histos.add("multDistributionsNoPileup/hSecondsNPVvsZNCampl", "", kTH3F, {axisSecondsSuperWideBins, {200, 0, 500, "n PV tracks"}, {200, 0, maxZNACenergyForTimeDepDistributions, "ZNC ampl"}});
      }

      // ### QA event selection bits
      int nEvSelBits = enNumEvSelBits;
      histos.add("hSecondsEventSelBits", "", kTH2F, {axisSecondsWideBins, {nEvSelBits, -0.5, nEvSelBits - 0.5, "Monitoring of event selection bits"}});
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

      // ### QA RCT flags
      int nRctFlagsTotal = enNumRctFlagsTotal;
      histos.add("hSecondsRCTflags", "", kTH2F, {axisSecondsWideBins, {nRctFlagsTotal + 2, -0.5, nRctFlagsTotal + 2 - 0.5, "Monitoring of RCT flags"}});
      axRctFlags = reinterpret_cast<TAxis*>(histos.get<TH2>(HIST("hSecondsRCTflags"))->GetYaxis());
      axRctFlags->SetBinLabel(1, "NcollisionsSel8");
      axRctFlags->SetBinLabel(2, "CcdbNotFound");
      axRctFlags->SetBinLabel(3 + kCPVBad, "CPVBad");
      axRctFlags->SetBinLabel(3 + kEMCBad, "EMCBad");
      axRctFlags->SetBinLabel(3 + kEMCLimAccMCRepr, "EMCLimAccMCRepr");
      axRctFlags->SetBinLabel(3 + kFDDBad, "FDDBad");
      axRctFlags->SetBinLabel(3 + kFT0Bad, "FT0Bad");
      axRctFlags->SetBinLabel(3 + kFV0Bad, "FV0Bad");
      axRctFlags->SetBinLabel(3 + kHMPBad, "HMPBad");
      axRctFlags->SetBinLabel(3 + kITSBad, "ITSBad");
      axRctFlags->SetBinLabel(3 + kITSLimAccMCRepr, "ITSLimAccMCRepr");
      axRctFlags->SetBinLabel(3 + kMCHBad, "MCHBad");
      axRctFlags->SetBinLabel(3 + kMCHLimAccMCRepr, "MCHLimAccMCRepr");
      axRctFlags->SetBinLabel(3 + kMFTBad, "MFTBad");
      axRctFlags->SetBinLabel(3 + kMFTLimAccMCRepr, "MFTLimAccMCRepr");
      axRctFlags->SetBinLabel(3 + kMIDBad, "MIDBad");
      axRctFlags->SetBinLabel(3 + kMIDLimAccMCRepr, "MIDLimAccMCRepr");
      axRctFlags->SetBinLabel(3 + kPHSBad, "PHSBad");
      axRctFlags->SetBinLabel(3 + kTOFBad, "TOFBad");
      axRctFlags->SetBinLabel(3 + kTOFLimAccMCRepr, "TOFLimAccMCRepr");
      axRctFlags->SetBinLabel(3 + kTPCBadTracking, "TPCBadTracking");
      axRctFlags->SetBinLabel(3 + kTPCBadPID, "TPCBadPID");
      axRctFlags->SetBinLabel(3 + kTPCLimAccMCRepr, "TPCLimAccMCRepr");
      axRctFlags->SetBinLabel(3 + kTRDBad, "TRDBad");
      axRctFlags->SetBinLabel(3 + kZDCBad, "ZDCBad");
      // combined flags
      axRctFlags->SetBinLabel(3 + enCBT, "CBT");
      axRctFlags->SetBinLabel(3 + enCBT_hadronPID, "CBT_hadronPID");
      axRctFlags->SetBinLabel(3 + enCBT_electronPID, "CBT_electronPID");
      axRctFlags->SetBinLabel(3 + enCBT_calo, "CBT_calo");
      axRctFlags->SetBinLabel(3 + enCBT_muon, "CBT_muon");
      axRctFlags->SetBinLabel(3 + enCBT_muon_glo, "CBT_muon_glo");

      // QA for all tracks
      // const AxisSpec axisChi2ITS{40, 0., 20., "chi2/ndof"};
      // const AxisSpec axisChi2TPC{40, 0., 20., "chi2/ndof"};
      const AxisSpec axisNclsITS{5, 3.5, 8.5, "n ITS cls"};
      const AxisSpec axisNclsTPC{40, -0.5, 159.5, "n TPC cls"};
      const AxisSpec axisFraction{20, 0, 1., "Fraction shared cls Tpc"};
      histos.add("allTracks/hSecondsTracks", "", kTH1D, {axisSeconds});
      if (confFlagCheckQoverPtHist) {
        histos.add("allTracks/hSecondsQoverPtSumDcaR", "", kTH2D, {axisSecondsWideBins, axisSparseQoverPt});
        histos.add("allTracks/hSecondsQoverPtSumDcaZ", "", kTH2D, {axisSecondsWideBins, axisSparseQoverPt});
      }
      histos.add("allTracks/hSecondsSumDcaR", "", kTH1D, {axisSeconds});
      histos.add("allTracks/hSecondsSumDcaZ", "", kTH1D, {axisSeconds});
      histos.add("allTracks/hSecondsSumPt", "", kTH1D, {axisSeconds});
      histos.add("allTracks/hSecondsNumClsIts", "", kTH1D, {axisSeconds});
      histos.add("allTracks/hSeconds2DNumClsIts", "", kTH2D, {axisSecondsWideBins, axisNclsITS});
      histos.add("allTracks/hSecondsChi2NClIts", "", kTH1D, {axisSeconds});
      if (confFlagCheckMshape)
        histos.add("allTracks/hSecondsTracksMshape", "", kTH1D, {axisSeconds});

      // QA for PV contributors
      histos.add("PVcontrib/hSecondsTracks", "", kTH1D, {axisSeconds});
      histos.add("PVcontrib/hSecondsSumDcaR", "", kTH1D, {axisSeconds});
      histos.add("PVcontrib/hSecondsSumDcaZ", "", kTH1D, {axisSeconds});
      histos.add("PVcontrib/hSecondsSumPt", "", kTH1D, {axisSeconds});
      histos.add("PVcontrib/hSecondsNumClsIts", "", kTH1D, {axisSeconds});
      histos.add("PVcontrib/hSeconds2DNumClsIts", "", kTH2D, {axisSecondsWideBins, axisNclsITS});
      histos.add("PVcontrib/hSecondsChi2NClIts", "", kTH1D, {axisSeconds});

      // QA for global tracks
      // ### A side
      // global tracks
      histos.add("A/global/hSecondsNumTracks", "", kTH1D, {axisSeconds});
      if (confFlagCheckQoverPtHist) {
        histos.add("A/global/hSecondsQoverPtSumDcaR", "", kTH2D, {axisSecondsWideBins, axisSparseQoverPt});
        histos.add("A/global/hSecondsQoverPtSumDcaZ", "", kTH2D, {axisSecondsWideBins, axisSparseQoverPt});
      }
      histos.add("A/global/hSecondsSumDcaR", "", kTH1D, {axisSeconds});
      histos.add("A/global/hSecondsSumDcaZ", "", kTH1D, {axisSeconds});
      histos.add("A/global/hSecondsSumPt", "", kTH1D, {axisSeconds});
      histos.add("A/global/hSecondsNumClsIts", "", kTH1D, {axisSeconds});
      histos.add("A/global/hSeconds2DNumClsIts", "", kTH2D, {axisSecondsWideBins, axisNclsITS});
      histos.add("A/global/hSecondsChi2NClIts", "", kTH1D, {axisSeconds});
      histos.add("A/global/hSecondsNumClsTpc", "", kTH1D, {axisSeconds});
      histos.add("A/global/hSeconds2DNumClsTpc", "", kTH2D, {axisSecondsWideBins, axisNclsTPC});
      histos.add("A/global/hSecondsChi2NClTpc", "", kTH1D, {axisSeconds});
      histos.add("A/global/hSecondsTpcFractionSharedCls", "", kTH1D, {axisSeconds});
      histos.add("A/global/hSeconds2DTpcFractionSharedCls", "", kTH2D, {axisSecondsWideBins, axisFraction});
      histos.add("A/global/hSecondsDeDx", "", kTH1D, {axisSeconds});

      // global && PV tracks
      histos.add("A/globalPV/hSecondsNumPVcontributors", "", kTH1D, {axisSeconds});
      histos.add("A/globalPV/hSecondsSumDcaR", "", kTH1D, {axisSeconds});
      histos.add("A/globalPV/hSecondsSumDcaZ", "", kTH1D, {axisSeconds});
      histos.add("A/globalPV/hSecondsSumPt", "", kTH1D, {axisSeconds});
      histos.add("A/globalPV/hSecondsNumClsIts", "", kTH1D, {axisSeconds});
      histos.add("A/globalPV/hSeconds2DNumClsIts", "", kTH2D, {axisSecondsWideBins, axisNclsITS});
      histos.add("A/globalPV/hSecondsChi2NClIts", "", kTH1D, {axisSeconds});
      histos.add("A/globalPV/hSecondsNumClsTpc", "", kTH1D, {axisSeconds});
      histos.add("A/globalPV/hSeconds2DNumClsTpc", "", kTH2D, {axisSecondsWideBins, axisNclsTPC});
      histos.add("A/globalPV/hSecondsChi2NClTpc", "", kTH1D, {axisSeconds});
      histos.add("A/globalPV/hSecondsTpcFractionSharedCls", "", kTH1D, {axisSeconds});
      histos.add("A/globalPV/hSeconds2DTpcFractionSharedCls", "", kTH2D, {axisSecondsWideBins, axisFraction});
      histos.add("A/globalPV/hSecondsDeDx", "", kTH1D, {axisSeconds});

      // ### C side
      // global tracks
      histos.add("C/global/hSecondsNumTracks", "", kTH1D, {axisSeconds});
      if (confFlagCheckQoverPtHist) {
        histos.add("C/global/hSecondsQoverPtSumDcaR", "", kTH2D, {axisSecondsWideBins, axisSparseQoverPt});
        histos.add("C/global/hSecondsQoverPtSumDcaZ", "", kTH2D, {axisSecondsWideBins, axisSparseQoverPt});
      }
      histos.add("C/global/hSecondsSumDcaR", "", kTH1D, {axisSeconds});
      histos.add("C/global/hSecondsSumDcaZ", "", kTH1D, {axisSeconds});
      histos.add("C/global/hSecondsSumPt", "", kTH1D, {axisSeconds});
      histos.add("C/global/hSecondsNumClsIts", "", kTH1D, {axisSeconds});
      histos.add("C/global/hSeconds2DNumClsIts", "", kTH2D, {axisSecondsWideBins, axisNclsITS});
      histos.add("C/global/hSecondsChi2NClIts", "", kTH1D, {axisSeconds});
      histos.add("C/global/hSecondsNumClsTpc", "", kTH1D, {axisSeconds});
      histos.add("C/global/hSeconds2DNumClsTpc", "", kTH2D, {axisSecondsWideBins, axisNclsTPC});
      histos.add("C/global/hSecondsChi2NClTpc", "", kTH1D, {axisSeconds});
      histos.add("C/global/hSecondsTpcFractionSharedCls", "", kTH1D, {axisSeconds});
      histos.add("C/global/hSeconds2DTpcFractionSharedCls", "", kTH2D, {axisSecondsWideBins, axisFraction});
      histos.add("C/global/hSecondsDeDx", "", kTH1D, {axisSeconds});

      // global && PV tracks
      histos.add("C/globalPV/hSecondsNumPVcontributors", "", kTH1D, {axisSeconds});
      histos.add("C/globalPV/hSecondsSumDcaR", "", kTH1D, {axisSeconds});
      histos.add("C/globalPV/hSecondsSumDcaZ", "", kTH1D, {axisSeconds});
      histos.add("C/globalPV/hSecondsSumPt", "", kTH1D, {axisSeconds});
      histos.add("C/globalPV/hSecondsNumClsIts", "", kTH1D, {axisSeconds});
      histos.add("C/globalPV/hSeconds2DNumClsIts", "", kTH2D, {axisSecondsWideBins, axisNclsITS});
      histos.add("C/globalPV/hSecondsChi2NClIts", "", kTH1D, {axisSeconds});
      histos.add("C/globalPV/hSecondsNumClsTpc", "", kTH1D, {axisSeconds});
      histos.add("C/globalPV/hSeconds2DNumClsTpc", "", kTH2D, {axisSecondsWideBins, axisNclsTPC});
      histos.add("C/globalPV/hSecondsChi2NClTpc", "", kTH1D, {axisSeconds});
      histos.add("C/globalPV/hSecondsTpcFractionSharedCls", "", kTH1D, {axisSeconds});
      histos.add("C/globalPV/hSeconds2DTpcFractionSharedCls", "", kTH2D, {axisSecondsWideBins, axisFraction});
      histos.add("C/globalPV/hSecondsDeDx", "", kTH1D, {axisSeconds});

      // phi holes vs time
      const AxisSpec axisPhi{64, 0, TMath::TwoPi(), "#varphi"}; // o2-linter: disable=external-pi (temporary fix)
      const AxisSpec axisEta{10, -0.8, 0.8, "#eta"};
      if (confFlagFillPhiVsTimeHist == 2) {
        histos.add("hSecondsITSlayer0vsPhi", "", kTH2F, {axisSeconds, axisPhi});
        histos.add("hSecondsITSlayer1vsPhi", "", kTH2F, {axisSeconds, axisPhi});
        histos.add("hSecondsITSlayer2vsPhi", "", kTH2F, {axisSeconds, axisPhi});
        histos.add("hSecondsITSlayer3vsPhi", "", kTH2F, {axisSeconds, axisPhi});
        histos.add("hSecondsITSlayer4vsPhi", "", kTH2F, {axisSeconds, axisPhi});
        histos.add("hSecondsITSlayer5vsPhi", "", kTH2F, {axisSeconds, axisPhi});
        histos.add("hSecondsITSlayer6vsPhi", "", kTH2F, {axisSeconds, axisPhi});
      }
      if (confFlagFillPhiVsTimeHist > 0) {
        histos.add("hSecondsITS7clsVsPhi", "", kTH2F, {axisSeconds, axisPhi});
        histos.add("hSecondsITSglobalVsPhi", "", kTH2F, {axisSeconds, axisPhi});
        histos.add("hSecondsITSTRDVsPhi", "", kTH2F, {axisSeconds, axisPhi});
        histos.add("hSecondsITSTOFVsPhi", "", kTH2F, {axisSeconds, axisPhi});
      }
      if (confFlagFillEtaPhiVsTimeHist)
        histos.add("hSecondsITSglobalVsEtaPhi", "", kTH3F, {axisSeconds, axisEta, axisPhi});
    }

    // count TVX triggers per DF
    for (const auto& bc : bcs) {
      // auto bc = col.foundBC_as<BCsRun3>();
      int64_t ts = bc.timestamp();
      double secFromSOR = ts / 1000. - minSec;
      if (bc.selection_bit(kIsTriggerTVX)) {
        histos.fill(HIST("hSecondsBCsTVX"), secFromSOR);

        uint64_t globalBC = bc.globalBC();
        int localBC = globalBC % nBCsPerOrbit;
        histos.fill(HIST("hSecondsBCsMap"), secFromSOR, localBC);

        if (bc.selection_bit(kNoTimeFrameBorder)) {
          histos.fill(HIST("hSecondsBCsTVXandTFborderCuts"), secFromSOR);
        }
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
        histos.fill(HIST("hSecondsCollisionsTVXNoVzCut"), secFromSOR);
      if (col.selection_bit(kNoTimeFrameBorder))
        histos.fill(HIST("hSecondsCollisionsTFborderCutNoVzCut"), secFromSOR);
      if (col.selection_bit(kIsTriggerTVX) && col.selection_bit(kNoTimeFrameBorder))
        histos.fill(HIST("hSecondsCollisionsTVXTFborderCutNoVzCut"), secFromSOR);

      histos.fill(HIST("hSecondsUPCverticesBeforeAllCuts"), secFromSOR, isVertexUPC ? 1 : 0);

      if (std::fabs(col.posZ()) > 10)
        continue;

      histos.fill(HIST("hSecondsUPCverticesBeforeSel8"), secFromSOR, isVertexUPC ? 1 : 0);

      histos.fill(HIST("hSecondsEventSelBits"), secFromSOR, enCollisionsAll);
      histos.fill(HIST("hSecondsEventSelBits"), secFromSOR, enIsTriggerTVX, col.selection_bit(kIsTriggerTVX));
      histos.fill(HIST("hSecondsEventSelBits"), secFromSOR, enNoTimeFrameBorder, col.selection_bit(kNoTimeFrameBorder));
      histos.fill(HIST("hSecondsEventSelBits"), secFromSOR, enNoITSROFrameBorder, col.selection_bit(kNoITSROFrameBorder));

      // sel8 selection:
      if (!col.sel8())
        continue;

      histos.fill(HIST("hSecondsUPCvertices"), secFromSOR, isVertexUPC ? 1 : 0);
      histos.fill(HIST("hSecondsCollisions"), secFromSOR);
      if (col.selection_bit(kNoSameBunchPileup))
        histos.fill(HIST("hSecondsCollisionsNoPileup"), secFromSOR);
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

      // check RCT flags
      histos.fill(HIST("hSecondsRCTflags"), secFromSOR, 0);                                 // n collisions sel8
      histos.fill(HIST("hSecondsRCTflags"), secFromSOR, 1, col.rct_bit(kCcdbObjectLoaded)); // CCDB object not loaded
      LOGP(debug, "i = 1, bitValue = {}, binLabel={}, binCenter={}", col.rct_bit(kCcdbObjectLoaded), axRctFlags->GetBinLabel(2), axRctFlags->GetBinCenter(2));
      for (int iFlag = 0; iFlag < kNRCTSelectionFlags; iFlag++) {
        histos.fill(HIST("hSecondsRCTflags"), secFromSOR, 2 + iFlag, col.rct_bit(iFlag));
        LOGP(debug, "i = {}, bitValue = {}, binLabel={}, binCenter={}", iFlag, col.rct_bit(iFlag), axRctFlags->GetBinLabel(3 + iFlag), axRctFlags->GetBinCenter(3 + iFlag));
      }
      LOGP(debug, "CBT_hadronPID = {}, kFT0Bad = {}, kITSBad = {}, kTPCBadTracking = {}, kTPCBadPID = {}, kTOFBad = {}, 1 + enCBT_hadronPID = {}, binLabel={}, binCenter={}", rctCheckerCBT_hadronPID(col),
           col.rct_bit(kFT0Bad), col.rct_bit(kITSBad), col.rct_bit(kTPCBadTracking), col.rct_bit(kTPCBadPID), col.rct_bit(kTOFBad), 2 + enCBT_hadronPID, axRctFlags->GetBinLabel(3 + enCBT_hadronPID), axRctFlags->GetBinCenter(3 + enCBT_hadronPID));
      histos.fill(HIST("hSecondsRCTflags"), secFromSOR, 2 + enCBT, rctCheckerCBT(col));
      histos.fill(HIST("hSecondsRCTflags"), secFromSOR, 2 + enCBT_hadronPID, rctCheckerCBT_hadronPID(col));
      histos.fill(HIST("hSecondsRCTflags"), secFromSOR, 2 + enCBT_electronPID, rctCheckerCBT_electronPID(col));
      histos.fill(HIST("hSecondsRCTflags"), secFromSOR, 2 + enCBT_calo, rctCheckerCBT_calo(col));
      histos.fill(HIST("hSecondsRCTflags"), secFromSOR, 2 + enCBT_muon, rctCheckerCBT_muon(col));
      histos.fill(HIST("hSecondsRCTflags"), secFromSOR, 2 + enCBT_muon_glo, rctCheckerCBT_muon_glo(col));

      // check hadronic rate
      double hadronicRate = mRateFetcher.fetch(ccdb.service, ts, runNumber, "ZNC hadronic") * 1.e-3; // kHz
      histos.fill(HIST("hSecondsIR"), secFromSOR, hadronicRate);

      // checking mShape flags in time:
      bool isMshape = false;
      if (confFlagCheckMshape) {
        auto mShapeTree = ccdb->getForTimeStamp<TTree>("TPC/Calib/MShapePotential", ts);
        mshape.setFromTree(*mShapeTree);
        isMshape = !mshape.getBoundaryPotential(ts).mPotential.empty();
      }

      // ##### track loop
      auto tracksGrouped = tracks.sliceBy(perCollision, col.globalIndex());
      int nPVtracks = 0;
      int nGlobalTracks = 0;
      for (const auto& track : tracksGrouped) {
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
        if (confFlagCheckQoverPtHist) {
          histos.fill(HIST("allTracks/hQoverPtDcaR"), qpt, dcaR);
          histos.fill(HIST("allTracks/hQoverPtDcaZ"), qpt, dcaZ);
        }
        // now consider only abs values for DCAs:
        double dcaRabs = std::fabs(dcaR);
        double dcaZabs = std::fabs(dcaZ);

        histos.fill(HIST("allTracks/hSecondsSumDcaR"), secFromSOR, dcaRabs);
        histos.fill(HIST("allTracks/hSecondsSumDcaZ"), secFromSOR, dcaZabs);
        histos.fill(HIST("allTracks/hSecondsSumPt"), secFromSOR, track.pt());
        if (confFlagCheckQoverPtHist) {
          histos.fill(HIST("allTracks/hSecondsQoverPtSumDcaR"), secFromSOR, qpt, dcaRabs);
          histos.fill(HIST("allTracks/hSecondsQoverPtSumDcaZ"), secFromSOR, qpt, dcaZabs);
        }
        histos.fill(HIST("allTracks/hSecondsNumClsIts"), secFromSOR, track.itsNCls());
        histos.fill(HIST("allTracks/hSeconds2DNumClsIts"), secFromSOR, track.itsNCls());
        histos.fill(HIST("allTracks/hSecondsChi2NClIts"), secFromSOR, track.itsChi2NCl());
        if (confFlagCheckMshape && isMshape) {
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
          histos.fill(HIST("PVcontrib/hSeconds2DNumClsIts"), secFromSOR, track.itsNCls());
          histos.fill(HIST("PVcontrib/hSecondsChi2NClIts"), secFromSOR, track.itsChi2NCl());

          nPVtracks++;
        }

        // ### global tracks
        float dedx = track.tpcSignal();
        if (track.isGlobalTrack()) {
          nGlobalTracks++;
          if (track.tgl() > 0.) { // A side
            histos.fill(HIST("A/global/hDcaRafterCuts"), dcaR);
            histos.fill(HIST("A/global/hDcaZafterCuts"), dcaZ);

            histos.fill(HIST("A/global/hSecondsNumTracks"), secFromSOR);
            if (confFlagCheckQoverPtHist) {
              histos.fill(HIST("A/global/hSecondsQoverPtSumDcaR"), secFromSOR, qpt, dcaRabs);
              histos.fill(HIST("A/global/hSecondsQoverPtSumDcaZ"), secFromSOR, qpt, dcaZabs);
            }
            histos.fill(HIST("A/global/hSecondsSumDcaR"), secFromSOR, dcaRabs);
            histos.fill(HIST("A/global/hSecondsSumDcaZ"), secFromSOR, dcaZabs);
            histos.fill(HIST("A/global/hSecondsSumPt"), secFromSOR, track.pt());
            histos.fill(HIST("A/global/hSecondsNumClsIts"), secFromSOR, track.itsNCls());
            histos.fill(HIST("A/global/hSeconds2DNumClsIts"), secFromSOR, track.itsNCls());
            histos.fill(HIST("A/global/hSecondsChi2NClIts"), secFromSOR, track.itsChi2NCl());
            histos.fill(HIST("A/global/hSecondsNumClsTpc"), secFromSOR, track.tpcNClsFound());
            histos.fill(HIST("A/global/hSeconds2DNumClsTpc"), secFromSOR, track.tpcNClsFound());
            histos.fill(HIST("A/global/hSecondsChi2NClTpc"), secFromSOR, track.tpcChi2NCl());
            if (track.tpcNClsFound() >= confCutOnNtpcClsForSharedFractAndDeDxCalc) {
              histos.fill(HIST("A/global/hSecondsTpcFractionSharedCls"), secFromSOR, track.tpcFractionSharedCls());
              histos.fill(HIST("A/global/hSeconds2DTpcFractionSharedCls"), secFromSOR, track.tpcFractionSharedCls());
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
              histos.fill(HIST("A/globalPV/hSeconds2DNumClsIts"), secFromSOR, track.itsNCls());
              histos.fill(HIST("A/globalPV/hSecondsChi2NClIts"), secFromSOR, track.itsChi2NCl());
              histos.fill(HIST("A/globalPV/hSecondsNumClsTpc"), secFromSOR, track.tpcNClsFound());
              histos.fill(HIST("A/globalPV/hSeconds2DNumClsTpc"), secFromSOR, track.tpcNClsFound());
              histos.fill(HIST("A/globalPV/hSecondsChi2NClTpc"), secFromSOR, track.tpcChi2NCl());
              if (track.tpcNClsFound() >= confCutOnNtpcClsForSharedFractAndDeDxCalc) {
                histos.fill(HIST("A/globalPV/hSecondsTpcFractionSharedCls"), secFromSOR, track.tpcFractionSharedCls());
                histos.fill(HIST("A/globalPV/hSeconds2DTpcFractionSharedCls"), secFromSOR, track.tpcFractionSharedCls());
                if (dedx < 1.e4) // protection from weird values
                  histos.fill(HIST("A/globalPV/hSecondsDeDx"), secFromSOR, dedx);
              }
            }
          } else { // C side
            histos.fill(HIST("C/global/hDcaRafterCuts"), dcaR);
            histos.fill(HIST("C/global/hDcaZafterCuts"), dcaZ);

            histos.fill(HIST("C/global/hSecondsNumTracks"), secFromSOR);
            if (confFlagCheckQoverPtHist) {
              histos.fill(HIST("C/global/hSecondsQoverPtSumDcaR"), secFromSOR, qpt, dcaRabs);
              histos.fill(HIST("C/global/hSecondsQoverPtSumDcaZ"), secFromSOR, qpt, dcaZabs);
            }
            histos.fill(HIST("C/global/hSecondsSumDcaR"), secFromSOR, dcaRabs);
            histos.fill(HIST("C/global/hSecondsSumDcaZ"), secFromSOR, dcaZabs);
            histos.fill(HIST("C/global/hSecondsSumPt"), secFromSOR, track.pt());
            histos.fill(HIST("C/global/hSecondsNumClsIts"), secFromSOR, track.itsNCls());
            histos.fill(HIST("C/global/hSeconds2DNumClsIts"), secFromSOR, track.itsNCls());
            histos.fill(HIST("C/global/hSecondsChi2NClIts"), secFromSOR, track.itsChi2NCl());
            histos.fill(HIST("C/global/hSecondsNumClsTpc"), secFromSOR, track.tpcNClsFound());
            histos.fill(HIST("C/global/hSeconds2DNumClsTpc"), secFromSOR, track.tpcNClsFound());
            histos.fill(HIST("C/global/hSecondsChi2NClTpc"), secFromSOR, track.tpcChi2NCl());
            if (track.tpcNClsFound() >= confCutOnNtpcClsForSharedFractAndDeDxCalc) {
              histos.fill(HIST("C/global/hSecondsTpcFractionSharedCls"), secFromSOR, track.tpcFractionSharedCls());
              histos.fill(HIST("C/global/hSeconds2DTpcFractionSharedCls"), secFromSOR, track.tpcFractionSharedCls());
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
              histos.fill(HIST("C/globalPV/hSeconds2DNumClsIts"), secFromSOR, track.itsNCls());
              histos.fill(HIST("C/globalPV/hSecondsChi2NClIts"), secFromSOR, track.itsChi2NCl());
              histos.fill(HIST("C/globalPV/hSecondsNumClsTpc"), secFromSOR, track.tpcNClsFound());
              histos.fill(HIST("C/globalPV/hSeconds2DNumClsTpc"), secFromSOR, track.tpcNClsFound());
              histos.fill(HIST("C/globalPV/hSecondsChi2NClTpc"), secFromSOR, track.tpcChi2NCl());
              if (track.tpcNClsFound() >= confCutOnNtpcClsForSharedFractAndDeDxCalc) {
                histos.fill(HIST("C/globalPV/hSecondsTpcFractionSharedCls"), secFromSOR, track.tpcFractionSharedCls());
                histos.fill(HIST("C/globalPV/hSeconds2DTpcFractionSharedCls"), secFromSOR, track.tpcFractionSharedCls());
                if (dedx < 1.e4) // protection from weird values
                  histos.fill(HIST("C/globalPV/hSecondsDeDx"), secFromSOR, dedx);
              }
            }
          }
        } // end of global tracks

        // study ITS cluster pattern vs phi vs time (pt>1 GeV/c cut selects straight tracks)
        if (track.isPVContributor() && track.pt() > 1) {
          // layer-by-layer check
          if (confFlagFillPhiVsTimeHist == 2) {
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
          if (confFlagFillPhiVsTimeHist > 0) {
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
          if (confFlagFillEtaPhiVsTimeHist && track.isGlobalTrack()) {
            histos.fill(HIST("hSecondsITSglobalVsEtaPhi"), secFromSOR, track.eta(), track.phi());
          }
        }
      } // end of track loop

      // fill mult distributions vs time
      if (confIncludeMultDistrVsTimeHistos && (!confFlagUseGoodZvtxFT0vsPVForTimeDepDistributions ? true : col.selection_bit(kIsGoodZvtxFT0vsPV))) {
        bool noPileup = col.selection_bit(kNoSameBunchPileup);

        int nTracksForTimeDep = confFlagUseGlobalTracksForTimeDepDistributions ? nGlobalTracks : nPVtracks;

        histos.fill(HIST("multDistributions/hSecondsDistrPVtracks"), secFromSOR, nTracksForTimeDep);

        // ZNA,C
        // float multZNA = bc.has_zdc() ? bc.zdc().energyCommonZNA() : -999.f;
        // float multZNC = bc.has_zdc() ? bc.zdc().energyCommonZNC() : -999.f;
        histos.fill(HIST("multDistributions/hSecondsDistrZNA"), secFromSOR, col.multZNA());
        histos.fill(HIST("multDistributions/hSecondsDistrZNC"), secFromSOR, col.multZNC());
        float ZNdiff = col.multZNA() - col.multZNC();
        float ZNsum = col.multZNA() + col.multZNC();
        histos.fill(HIST("multDistributions/hSecondsDistrZNACdiff"), secFromSOR, ZNdiff);
        if (ZNsum > 0)
          histos.fill(HIST("multDistributions/hSecondsDistrZNACdiffNorm"), secFromSOR, ZNdiff / ZNsum);

        // ZNA,C by amplitudes (suggested by Chiara O.)
        float ZNAampl = bc.has_zdc() ? bc.zdc().amplitudeZNA() : 0;
        float ZNCampl = bc.has_zdc() ? bc.zdc().amplitudeZNC() : 0;
        histos.fill(HIST("multDistributions/hSecondsDistrZNAampl"), secFromSOR, ZNAampl);
        histos.fill(HIST("multDistributions/hSecondsDistrZNCampl"), secFromSOR, ZNCampl);
        float ZNdiffAmpl = ZNAampl - ZNCampl;
        float ZNsumAmpl = ZNAampl + ZNCampl;
        histos.fill(HIST("multDistributions/hSecondsDistrZNACdiffAmpl"), secFromSOR, ZNdiffAmpl);
        if (ZNsumAmpl > 0)
          histos.fill(HIST("multDistributions/hSecondsDistrZNACdiffNormAmpl"), secFromSOR, ZNdiffAmpl / ZNsumAmpl);

        histos.fill(HIST("multDistributions/hSecondsNPVvsZNAampl"), secFromSOR, nTracksForTimeDep, ZNAampl);
        histos.fill(HIST("multDistributions/hSecondsNPVvsZNCampl"), secFromSOR, nTracksForTimeDep, ZNCampl);

        // FT0A,C, V0A
        // float multT0A = bc.has_ft0() ? bc.ft0().sumAmpA() : -999.f;
        // float multT0C = bc.has_ft0() ? fbcundBC.ft0().sumAmpC() : -999.f;
        histos.fill(HIST("multDistributions/hSecondsDistrT0A"), secFromSOR, col.multFT0A());
        histos.fill(HIST("multDistributions/hSecondsDistrT0C"), secFromSOR, col.multFT0C());
        histos.fill(HIST("multDistributions/hSecondsDistrV0A"), secFromSOR, col.multFV0A());

        if (noPileup) {
          histos.fill(HIST("multDistributionsNoPileup/hSecondsDistrPVtracks"), secFromSOR, nTracksForTimeDep);

          histos.fill(HIST("multDistributionsNoPileup/hSecondsDistrZNA"), secFromSOR, col.multZNA());
          histos.fill(HIST("multDistributionsNoPileup/hSecondsDistrZNC"), secFromSOR, col.multZNC());
          histos.fill(HIST("multDistributionsNoPileup/hSecondsDistrZNACdiff"), secFromSOR, ZNdiff);
          if (ZNsum > 0)
            histos.fill(HIST("multDistributionsNoPileup/hSecondsDistrZNACdiffNorm"), secFromSOR, ZNdiff / ZNsum);

          histos.fill(HIST("multDistributionsNoPileup/hSecondsDistrZNAampl"), secFromSOR, ZNAampl);
          histos.fill(HIST("multDistributionsNoPileup/hSecondsDistrZNCampl"), secFromSOR, ZNCampl);
          histos.fill(HIST("multDistributionsNoPileup/hSecondsDistrZNACdiffAmpl"), secFromSOR, ZNdiffAmpl);
          if (ZNsumAmpl > 0)
            histos.fill(HIST("multDistributionsNoPileup/hSecondsDistrZNACdiffNormAmpl"), secFromSOR, ZNdiffAmpl / ZNsumAmpl);

          histos.fill(HIST("multDistributionsNoPileup/hSecondsNPVvsZNAampl"), secFromSOR, nTracksForTimeDep, ZNAampl);
          histos.fill(HIST("multDistributionsNoPileup/hSecondsNPVvsZNCampl"), secFromSOR, nTracksForTimeDep, ZNCampl);

          //
          histos.fill(HIST("multDistributionsNoPileup/hSecondsDistrT0A"), secFromSOR, col.multFT0A());
          histos.fill(HIST("multDistributionsNoPileup/hSecondsDistrT0C"), secFromSOR, col.multFT0C());
          histos.fill(HIST("multDistributionsNoPileup/hSecondsDistrV0A"), secFromSOR, col.multFV0A());
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
