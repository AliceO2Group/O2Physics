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
///
/// \file neutronProtonCorrZdc.cxx
/// \brief Correlations between protons and neutrons in the ZDC
/// \author Olaf Massen <olaf.massen@cern.ch>

#include <string>
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"
#include "Framework/StaticFor.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

enum EventCounter { kNoSelection = 0,
                    kSel8 = 1,
                    kNoSameBunchPileUp = 2,
                    kIsGoodZvtxFT0vsPV = 3,
                    kNoCollInTimeRangeStandard = 4,
                    kMaxCentralitySelection = 5,
                    kZDCSelection = 6,
                    kTimeDifferenceZDC = 7 };

struct NeutronProtonCorrZdc {
  // Histogram registry: an object to hold your histograms
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<float> cfgZVertexCut{"cfgZVertexCut", 10., "Cut on Z vertex position"};
  Configurable<bool> cfgNoSameBunchPileupCut{"cfgNoSameBunchPileupCut", true, "kNoSameBunchPileUp Cut"};
  Configurable<bool> cfgIsGoodZvtxFT0vsPV{"cfgIsGoodZvtxFT0vsPV", true, "kIsGoodZvtxFT0vsPV Cut"};
  Configurable<bool> cfgNoCollInTimeRangeStandard{"cfgNoCollInTimeRangeStandard", true, "kNoCollInTimeRangeStandard Cut"};
  Configurable<double> cfgMaxCentrality{"cfgMaxCentrality", 80, "Maximum collision centrality"};
  Configurable<bool> cfgZDCTimingInformationCut{"cfgZDCTimingInformationCut", true, "Use timing information in ZDC event selection"};
  Configurable<int> cfgTimingBins{"cfgTimingBins", 200, "N bins for timing histograms"};
  Configurable<float> cfgTDCZNmincut{"cfgTDCZNmincut", -3.0, "Min ZN TDC cut"};
  Configurable<float> cfgTDCZNmaxcut{"cfgTDCZNmaxcut", 3.0, "Max ZN TDC cut"};
  Configurable<float> cfgTDCZPmincut{"cfgTDCZPmincut", -3.0, "Min ZP TDC cut"};
  Configurable<float> cfgTDCZPmaxcut{"cfgTDCZPmaxcut", 3.0, "Max ZP TDC cut"};
  Configurable<int> cfgNBinsZN{"cfgNBinsZN", 100, "N bins for ZNA and ZNC"};
  Configurable<int> cfgNBinsZP{"cfgNBinsZP", 100, "N bins for ZPA and ZPC"};
  Configurable<double> cfgZNmin{"cfgZNmin", -10, "Minimum value for ZN signal"};
  Configurable<double> cfgZNmax{"cfgZNmax", 350, "Maximum value for ZN signal"};
  Configurable<double> cfgZPmin{"cfgZPmin", -10, "Minimum value for ZP signal"};
  Configurable<double> cfgZPmax{"cfgZPmax", 150, "Maximum value for ZP signal"};
  Configurable<double> cfgAmplMin{"cfgAmplMin", -50, "Minimum value for sector amplitude"};
  Configurable<double> cfgAmplMax{"cfgAmplMax", 300, "Maximum value for sector amplitude"};
  Configurable<int> cfgNBinsAmpl{"cfgNBinsAmpl", 100, "Number of bins for sector amplitude histograms"};
  Configurable<double> cfgDiffZmin{"cfgDiffZmin", -30, "Minimum value for the diffZ signal"};
  Configurable<double> cfgDiffZmax{"cfgDiffZmax", 50, "Maximum value for the diffZ signal"};
  Configurable<int> cfgNBinsAlpha{"cfgNBinsAlpha", 100, "Number of bins for ZDC asymmetry"};
  Configurable<double> cfgAlphaZmin{"cfgAlphaZmin", -1, "Minimum value for ZDC asymmetry"};
  Configurable<double> cfgAlphaZmax{"cfgAlphaZmax", 1, "Maximum value for ZDC asymmetry"};
  Configurable<int> cfgCentralityEstimator{"cfgCentralityEstimator", 0, "Choice of centrality estimator"};
  Configurable<bool> cfgFillMultiplicityQAHistograms{"cfgFillMultiplicityQAHistograms", true, "Fill multiplicity QA plots"};
  Configurable<int> cfgNBinsMultiplicity{"cfgNBinsMultiplicity", 500, "N bins for multiplicity histograms"};

  ConfigurableAxis cfgAxisCent{"cfgAxisCent", {VARIABLE_WIDTH, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0, 31.0, 32.0, 33.0, 34.0, 35.0, 36.0, 37.0, 38.0, 39.0, 40.0, 41.0, 42.0, 43.0, 44.0, 45.0, 46.0, 47.0, 48.0, 49.0, 50.0, 51.0, 52.0, 53.0, 54.0, 55.0, 56.0, 57.0, 58.0, 59.0, 60.0, 61.0, 62.0, 63.0, 64.0, 65.0, 66.0, 67.0, 68.0, 69.0, 70.0, 71.0, 72.0, 73.0, 74.0, 75.0, 76.0, 77.0, 78.0, 79.0, 80.0, 81.0, 82.0, 83.0, 84.0, 85.0, 86.0, 87.0, 88.0, 89.0, 90.0, 91.0, 92.0, 93.0, 94.0, 95.0, 96.0, 97.0, 98.0, 99.0, 100.0}, "Centrality [%]"};

  Filter collisionVtxZ = nabs(aod::collision::posZ) < cfgZVertexCut;

  using CentralitiesRun3 = soa::Join<aod::CentFT0Cs, aod::CentFT0Ms, aod::CentFV0As, aod::CentFT0CVariant1s, aod::CentNGlobals>;
  using CentralitiesRun2 = aod::CentRun2V0Ms;
  using BCsRun3 = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;

  void init(InitContext const&)
  {
    // define axes you want to use
    const AxisSpec axisCounter{8, -0.5, 7.5, ""};
    const AxisSpec axisZDCTiming{cfgTimingBins, -10, 10};
    const AxisSpec axisZNSectorSignal{cfgNBinsZN, cfgZNmin, cfgZNmax / 3.};
    const AxisSpec axisZPSectorSignal{cfgNBinsZP, cfgZPmin, cfgZPmax / 3.};
    const AxisSpec axisZNASignal{cfgNBinsZN, cfgZNmin, cfgZNmax, "ZNA (a.u.)"};
    const AxisSpec axisZNCSignal{cfgNBinsZN, cfgZNmin, cfgZNmax, "ZNC (a.u.)"};
    const AxisSpec axisZPASignal{cfgNBinsZP, cfgZPmin, cfgZPmax, "ZPA (a.u.)"};
    const AxisSpec axisZPCSignal{cfgNBinsZP, cfgZPmin, cfgZPmax, "ZPC (a.u.)"};
    const AxisSpec axisZNSignal{2 * cfgNBinsZN, cfgZNmin, 1.5 * cfgZNmax, "ZN (a.u.)"};
    const AxisSpec axisZPSignal{2 * cfgNBinsZP, cfgZPmin, 1.5 * cfgZPmax, "ZP (a.u.)"};
    const AxisSpec axisZNAmplitude{cfgNBinsAmpl, cfgAmplMin, cfgAmplMax, "Ampl (a.u.)"};
    const AxisSpec axisZPAmplitude{cfgNBinsAmpl, cfgAmplMin * 0.3, cfgAmplMax * 0.3, "Ampl (a.u.)"};
    const AxisSpec axisTotalAmplitude{cfgNBinsAmpl, 4.0 * cfgAmplMin, 4.0 * cfgAmplMax, "Ampl (a.u.)"};
    const AxisSpec axisAlphaZ{cfgNBinsAlpha, cfgAlphaZmin, cfgAlphaZmax, "#alpha_{spec}"};
    const AxisSpec axisZDiffSignal{cfgNBinsZN, cfgDiffZmin, cfgDiffZmax, "#Delta E"};
    const AxisSpec axisMultiplicityF0A{cfgNBinsMultiplicity, 0, 200000, "F0A"};
    const AxisSpec axisMultiplicityF0C{cfgNBinsMultiplicity, 0, 100000, "F0C"};
    const AxisSpec axisMultiplicityF0M{cfgNBinsMultiplicity, 0, 300000, "F0M"};
    const AxisSpec axisMultiplicityFDD{cfgNBinsMultiplicity, 0, 50000, "FDD"};
    const AxisSpec axisMultiplicityTPC{cfgNBinsMultiplicity, 0, 100000, "TPC"};
    const AxisSpec axisMultiplicityMultNGlobal{cfgNBinsMultiplicity, 0, 3500, "MultsNGlobal"};

    HistogramConfigSpec defaultTimingHistogram({HistType::kTH2F, {cfgAxisCent, axisZDCTiming}});
    HistogramConfigSpec defaultZNSectorHist({HistType::kTH2F, {cfgAxisCent, axisZNSectorSignal}});
    HistogramConfigSpec defaultZPSectorHist({HistType::kTH2F, {cfgAxisCent, axisZPSectorSignal}});
    HistogramConfigSpec defaultZNAmplSectorHist({HistType::kTH2F, {cfgAxisCent, axisZNAmplitude}});
    HistogramConfigSpec defaultZPAmplSectorHist({HistType::kTH2F, {cfgAxisCent, axisZPAmplitude}});
    HistogramConfigSpec defaultEnergyZNAmplSectorHist({HistType::kTH2F, {axisZNSignal, axisZNAmplitude}});
    HistogramConfigSpec defaultEnergyZPAmplSectorHist({HistType::kTH2F, {axisZPSignal, axisZPAmplitude}});
    HistogramConfigSpec defaultZDCDiffHist({HistType::kTH2F, {cfgAxisCent, axisZDiffSignal}});

    // create histograms
    histos.add("eventCounter", "eventCounter", kTH1F, {axisCounter});
    histos.get<TH1>(HIST("eventCounter"))->GetXaxis()->SetBinLabel(EventCounter::kSel8 + 1, "Sel8");
    histos.get<TH1>(HIST("eventCounter"))->GetXaxis()->SetBinLabel(EventCounter::kNoSameBunchPileUp + 1, "kNoSameBunchPileup");
    histos.get<TH1>(HIST("eventCounter"))->GetXaxis()->SetBinLabel(EventCounter::kIsGoodZvtxFT0vsPV + 1, "kIsGoodZvtxFT0vsPV");
    histos.get<TH1>(HIST("eventCounter"))->GetXaxis()->SetBinLabel(EventCounter::kNoCollInTimeRangeStandard + 1, "kNoCollInTimeRangeStandard");
    histos.get<TH1>(HIST("eventCounter"))->GetXaxis()->SetBinLabel(EventCounter::kMaxCentralitySelection + 1, "Cenrality range");
    histos.get<TH1>(HIST("eventCounter"))->GetXaxis()->SetBinLabel(EventCounter::kZDCSelection + 1, "isSelectedZDC");
    histos.get<TH1>(HIST("eventCounter"))->GetXaxis()->SetBinLabel(EventCounter::kTimeDifferenceZDC + 1, "ZDC time difference");

    histos.add("CentralityPercentile", "CentralityPercentile", kTH1F, {cfgAxisCent});
    histos.add("TimingZNAvsCent", "TimingZNAvsCent", defaultTimingHistogram);
    histos.add("TimingZNCvsCent", "TimingZNCvsCent", defaultTimingHistogram);
    histos.add("TimingZPAvsCent", "TimingZPAvsCent", defaultTimingHistogram);
    histos.add("TimingZPCvsCent", "TimingZPCvsCent", defaultTimingHistogram);

    histos.add("ASide/CentvsZNSector0Signal", "CentvsZNASector0Signal", defaultZNSectorHist);
    histos.add("ASide/CentvsZNSector1Signal", "CentvsZNASector1Signal", defaultZNSectorHist);
    histos.add("ASide/CentvsZNSector2Signal", "CentvsZNASector2Signal", defaultZNSectorHist);
    histos.add("ASide/CentvsZNSector3Signal", "CentvsZNASector3Signal", defaultZNSectorHist);
    histos.add("ASide/CentvsZPSector0Signal", "CentvsZPASector0Signal", defaultZPSectorHist);
    histos.add("ASide/CentvsZPSector1Signal", "CentvsZPASector1Signal", defaultZPSectorHist);
    histos.add("ASide/CentvsZPSector2Signal", "CentvsZPASector2Signal", defaultZPSectorHist);
    histos.add("ASide/CentvsZPSector3Signal", "CentvsZPASector3Signal", defaultZPSectorHist);
    histos.add("ASide/CentvsZNSignalSum", "CentvsZNASignalSum", kTH2F, {cfgAxisCent, axisZNASignal});
    histos.add("ASide/CentvsZNSignalCommon", "CentvsZNASignalCommon", kTH2F, {cfgAxisCent, axisZNASignal});
    histos.add("ASide/CentvsZPSignalSum", "CentvsZNASignalSum", kTH2F, {cfgAxisCent, axisZPASignal});
    histos.add("ASide/CentvsZPSignalCommon", "CentvsZNASignalCommon", kTH2F, {cfgAxisCent, axisZPASignal});
    histos.add("ASide/CentvsdiffZNSignal", "CentvsdiffZNSignal", defaultZDCDiffHist);
    histos.add("ASide/CentvsdiffZPSignal", "CentvsdiffZPSignal", defaultZDCDiffHist);

    histos.add("ASide/CentvsZNAmplitude", "CentvsZNAmplitude", defaultZNAmplSectorHist);
    histos.add("ASide/CentvsZPAmplitude", "CentvsZPAmplitude", defaultZPAmplSectorHist);
    histos.add("ASide/ZNSignalvsZNAmplitudeSum", "ZNSignalvsZNAmplitudeSum", defaultEnergyZNAmplSectorHist);
    histos.add("ASide/ZNSignalvsZNAmplitudeCommon", "ZNSignalvsZNAmplitudeCommon", defaultEnergyZNAmplSectorHist);
    histos.add("ASide/ZPSignalvsZPAmplitudeSum", "ZPSignalvsZPAmplitudeSum", defaultEnergyZPAmplSectorHist);
    histos.add("ASide/ZPSignalvsZPAmplitudeCommon", "ZPSignalvsZPAmplitudeCommon", defaultEnergyZPAmplSectorHist);

    // Cloning the folder
    histos.addClone("ASide/", "CSide/");

    histos.add("CentvsZNSignalCommon", "CentvsZNSignalCommon", kTH2F, {cfgAxisCent, axisZNSignal});
    histos.add("CentvsZNSignalSum", "CentvsZNSignalSum", kTH2F, {cfgAxisCent, axisZNSignal});
    histos.add("CentvsZPSignalCommon", "CentvsZPSignalCommon", kTH2F, {cfgAxisCent, axisZPSignal});
    histos.add("CentvsZPSignalSum", "CentvsZPSignalSum", kTH2F, {cfgAxisCent, axisZPSignal});
    histos.add("CentvsAlphaZN", "CentvsAlphaZN", kTH2F, {cfgAxisCent, axisAlphaZ});
    histos.add("CentvsAlphaZP", "CentvsAlphaZP", kTH2F, {cfgAxisCent, axisAlphaZ});
    histos.add("CentvsAlphaZNcommon", "CentvsAlphaZNcommon", kTH2F, {cfgAxisCent, axisAlphaZ});
    histos.add("CentvsAlphaZPcommon", "CentvsAlphaZPcommon", kTH2F, {cfgAxisCent, axisAlphaZ});
    histos.add("CentvsAlphaZNAmplitude", "CentvsAlphaZNAmplitude", kTH2F, {cfgAxisCent, axisAlphaZ});
    histos.add("CentvsAlphaZPAmplitude", "CentvsAlphaZPAmplitude", kTH2F, {cfgAxisCent, axisAlphaZ});
    histos.add("CentvsDiffZNSignal", "CentvsDiffZNSignal", defaultZDCDiffHist);
    histos.add("CentvsDiffZPSignal", "CentvsDiffZPSignal", defaultZDCDiffHist);
    histos.add("CentvsZNAvsZNC", "CentvsZNAvsZNC", kTH3F, {cfgAxisCent, axisZNASignal, axisZNCSignal});
    histos.add("CentvsZNAvsZPA", "CentvsZNAvsZPA", kTH3F, {cfgAxisCent, axisZNASignal, axisZPASignal});
    histos.add("CentvsZNAvsZPC", "CentvsZNAvsZPC", kTH3F, {cfgAxisCent, axisZNASignal, axisZPCSignal});
    histos.add("CentvsZPAvsZNC", "CentvsZPAvsZNC", kTH3F, {cfgAxisCent, axisZPASignal, axisZNCSignal});
    histos.add("CentvsZPAvsZPC", "CentvsZNAvsZPC", kTH3F, {cfgAxisCent, axisZPASignal, axisZPCSignal});
    histos.add("CentvsZNCvsZPC", "CentvsZNCvsZPC", kTH3F, {cfgAxisCent, axisZNCSignal, axisZPCSignal});
    histos.add("CentvsZNvsZP", "CentvsZNvsZP", kTH3F, {cfgAxisCent, axisZNSignal, axisZPSignal});

    if (cfgFillMultiplicityQAHistograms) {
      histos.add("MultiplicityHistograms/FV0A", "FV0A", kTH1F, {axisMultiplicityF0A});
      histos.add("MultiplicityHistograms/FT0A", "FT0A", kTH1F, {axisMultiplicityF0A});
      histos.add("MultiplicityHistograms/FT0C", "FT0C", kTH1F, {axisMultiplicityF0C});
      histos.add("MultiplicityHistograms/FDDA", "FDDA", kTH1F, {axisMultiplicityFDD});
      histos.add("MultiplicityHistograms/FDDC", "FDDC", kTH1F, {axisMultiplicityFDD});
      histos.add("MultiplicityHistograms/TPC", "TPC", kTH1F, {axisMultiplicityTPC});
      histos.add("MultiplicityHistograms/NGlobal", "NGlobal", kTH1F, {axisMultiplicityMultNGlobal});
      histos.add("MultiplicityHistograms/CentvsFT0C", "CentvsFT0C", kTH2F, {cfgAxisCent, axisMultiplicityF0C});
      histos.add("MultiplicityHistograms/CentvsFT0CVar1", "CentvsFT0CVar1", kTH2F, {cfgAxisCent, axisMultiplicityF0C});
      histos.add("MultiplicityHistograms/CentvsFT0M", "CentvsFT0M", kTH2F, {cfgAxisCent, axisMultiplicityF0M});
      histos.add("MultiplicityHistograms/CentvsFV0A", "CentvsFV0A", kTH2F, {cfgAxisCent, axisMultiplicityF0A});
      histos.add("MultiplicityHistograms/CentvsNGlobal", "CentvsNGlobal", kTH2F, {cfgAxisCent, axisMultiplicityMultNGlobal});
    }
  }

  template <typename TCollision>
  bool eventSelected(TCollision coll, const float centrality)
  {
    if (!coll.sel8())
      return 0;
    histos.fill(HIST("eventCounter"), kSel8);

    if (cfgNoSameBunchPileupCut) {
      if (!coll.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
        // rejects collisions which are associated with the same "found-by-T0" bunch crossing
        // https://indico.cern.ch/event/1396220/#1-event-selection-with-its-rof
        return 0;
      }
      histos.fill(HIST("eventCounter"), EventCounter::kNoSameBunchPileUp);
    }

    if (cfgIsGoodZvtxFT0vsPV) {
      if (!coll.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
        // removes collisions with large differences between z of PV by tracks and z of PV from FT0 A-C time difference
        // use this cut at low multiplicities with caution
        return 0;
      }
      histos.fill(HIST("eventCounter"), EventCounter::kIsGoodZvtxFT0vsPV);
    }

    if (cfgNoCollInTimeRangeStandard) {
      if (!coll.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
        //  Rejection of the collisions which have other events nearby
        return 0;
      }
      histos.fill(HIST("eventCounter"), EventCounter::kNoCollInTimeRangeStandard);
    }

    if (centrality > cfgMaxCentrality) {
      return 0;
    }
    histos.fill(HIST("eventCounter"), EventCounter::kMaxCentralitySelection);

    return 1;
  }

  template <int mult, typename C>
  void fillMultHistosRun3(const C& col)
  {
    static constexpr std::string_view MultLabels[] = {"FT0C", "FT0A", "FV0A", "FDDC", "FDDA", "TPC", "NGlobal"};
    std::array<float, 7> multarray = {col.multFT0C(), col.multFT0A(), col.multFV0A(), col.multFDDC(), col.multFDDA(), static_cast<float>(col.multTPC()), static_cast<float>(col.multNTracksGlobal())};

    histos.fill(HIST("MultiplicityHistograms/") + HIST(MultLabels[mult]), multarray[mult]);
  }

  template <int cent, typename C>
  void fillCentHistosRun3(const C& col)
  {
    static constexpr std::string_view CentLabels[] = {"CentvsFT0C", "CentvsFT0CVar1", "CentvsFT0M", "CentvsFV0A", "CentvsNGlobal"};
    std::array<float, 5> centarray = {col.centFT0C(), col.centFT0CVariant1(), col.centFT0M(), col.centFV0A(), col.centNGlobal()};
    std::array<float, 5> multarray = {col.multFT0C(), col.multFT0C(), col.multFT0C() + col.multFT0A(), col.multFV0A(), static_cast<float>(col.multNTracksGlobal())};

    histos.fill(HIST("MultiplicityHistograms/") + HIST(CentLabels[cent]), centarray[cent], multarray[cent]);
  }

  template <int side, typename Z>
  void fillZDCSideCommonHistos(const float centr, const Z& zdc)
  {
    static constexpr std::string_view SubDir[] = {"ASide/", "CSide/"};

    std::array<std::array<float, 4>, 2> znEnergyResponse = {zdc.energySectorZNA(), zdc.energySectorZNC()};
    std::array<std::array<float, 4>, 2> zpEnergyResponse = {zdc.energySectorZPA(), zdc.energySectorZPC()};
    std::array<float, 2> znEnergyResponseCommon = {zdc.energyCommonZNA(), zdc.energyCommonZNC()};
    std::array<float, 2> zpEnergyResponseCommon = {zdc.energyCommonZPA(), zdc.energyCommonZPC()};
    std::array<float, 2> znAmplitudeResponse = {zdc.amplitudeZNA(), zdc.amplitudeZNC()};
    std::array<float, 2> zpAmplitudeResponse = {zdc.amplitudeZPA(), zdc.amplitudeZPC()};

    float sumZN = znEnergyResponse[side][0] + znEnergyResponse[side][1] + znEnergyResponse[side][2] + znEnergyResponse[side][3];
    float sumZP = zpEnergyResponse[side][0] + zpEnergyResponse[side][1] + zpEnergyResponse[side][2] + zpEnergyResponse[side][3];

    histos.fill(HIST(SubDir[side]) + HIST("CentvsZNSignalSum"), centr, sumZN);
    histos.fill(HIST(SubDir[side]) + HIST("CentvsZNSignalCommon"), centr, znEnergyResponseCommon[side]);
    histos.fill(HIST(SubDir[side]) + HIST("CentvsdiffZNSignal"), centr, sumZN - znEnergyResponseCommon[side]);
    histos.fill(HIST(SubDir[side]) + HIST("CentvsZPSignalSum"), centr, sumZP);
    histos.fill(HIST(SubDir[side]) + HIST("CentvsZPSignalCommon"), centr, zpEnergyResponseCommon[side]);
    histos.fill(HIST(SubDir[side]) + HIST("CentvsdiffZPSignal"), centr, sumZP - zpEnergyResponseCommon[side]);
    histos.fill(HIST(SubDir[side]) + HIST("CentvsZNAmplitude"), centr, znAmplitudeResponse[side]);
    histos.fill(HIST(SubDir[side]) + HIST("CentvsZPAmplitude"), centr, zpAmplitudeResponse[side]);

    histos.fill(HIST(SubDir[side]) + HIST("ZNSignalvsZNAmplitudeSum"), sumZN, znAmplitudeResponse[side]);
    histos.fill(HIST(SubDir[side]) + HIST("ZNSignalvsZNAmplitudeCommon"), znEnergyResponseCommon[side], znAmplitudeResponse[side]);
    histos.fill(HIST(SubDir[side]) + HIST("ZPSignalvsZPAmplitudeSum"), sumZP, zpAmplitudeResponse[side]);
    histos.fill(HIST(SubDir[side]) + HIST("ZPSignalvsZPAmplitudeCommon"), zpEnergyResponseCommon[side], zpAmplitudeResponse[side]);
  }

  template <int side, int sector, typename Z>
  void fillZDCSideSectorHistos(const float centr, const Z& zdc)
  {
    static constexpr std::string_view SubDir[] = {"ASide/", "CSide/"};
    static constexpr std::string_view ZNSector[] = {"CentvsZNSector0Signal", "CentvsZNSector1Signal", "CentvsZNSector2Signal", "CentvsZNSector3Signal"};
    static constexpr std::string_view ZPSector[] = {"CentvsZPSector0Signal", "CentvsZPSector1Signal", "CentvsZPSector2Signal", "CentvsZPSector3Signal"};

    std::array<std::array<float, 4>, 2> znEnergyResponse = {zdc.energySectorZNA(), zdc.energySectorZNC()};
    std::array<std::array<float, 4>, 2> zpEnergyResponse = {zdc.energySectorZPA(), zdc.energySectorZPC()};

    histos.fill(HIST(SubDir[side]) + HIST(ZNSector[sector]), centr, znEnergyResponse[side][sector]);
    histos.fill(HIST(SubDir[side]) + HIST(ZPSector[sector]), centr, zpEnergyResponse[side][sector]);
  }

  void processRun3(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::MultsGlobal, CentralitiesRun3>>::iterator const& collision, BCsRun3 const&, aod::Zdcs const&)
  {
    histos.fill(HIST("eventCounter"), EventCounter::kNoSelection);

    const float centArray[] = {collision.centFT0C(), collision.centFT0CVariant1(), collision.centFT0M(), collision.centFV0A(), collision.centNGlobal()};
    const auto cent = centArray[cfgCentralityEstimator];

    if (!eventSelected(collision, cent))
      return;

    const auto& foundBC = collision.foundBC_as<BCsRun3>();
    if (foundBC.has_zdc()) {
      const auto& zdcread = foundBC.zdc();
      histos.fill(HIST("eventCounter"), EventCounter::kZDCSelection);

      auto tZNA = zdcread.timeZNA();
      auto tZNC = zdcread.timeZNC();
      auto tZPA = zdcread.timeZPA();
      auto tZPC = zdcread.timeZPC();

      histos.fill(HIST("TimingZNAvsCent"), cent, tZNA);
      histos.fill(HIST("TimingZNCvsCent"), cent, tZNC);
      histos.fill(HIST("TimingZPAvsCent"), cent, tZPA);
      histos.fill(HIST("TimingZPCvsCent"), cent, tZPC);

      // Selection on timing for the ZDC
      if (cfgZDCTimingInformationCut) {
        if (tZNA <= cfgTDCZNmincut || tZNA >= cfgTDCZNmaxcut) {
          return;
        }
        if (tZNC <= cfgTDCZNmincut || tZNC >= cfgTDCZNmaxcut) {
          return;
        }
        if (tZPA <= cfgTDCZPmincut || tZPA >= cfgTDCZPmaxcut) {
          return;
        }
        if (tZPC <= cfgTDCZPmincut || tZPC >= cfgTDCZPmaxcut) {
          return;
        }
      }
      histos.fill(HIST("eventCounter"), EventCounter::kTimeDifferenceZDC);
      histos.fill(HIST("CentralityPercentile"), cent);

      if (cfgFillMultiplicityQAHistograms) {
        static_for<0, 6>([&](auto i) {
          fillMultHistosRun3<i>(collision); // Fill multiplicity histograms
        });

        static_for<0, 4>([&](auto i) {
          fillCentHistosRun3<i>(collision); // Fill centrality vs multiplicity histograms
        });
      }

      static_for<0, 1>([&](auto i) {
        fillZDCSideCommonHistos<i>(cent, zdcread); // Fill i-side common histograms
        static_for<0, 3>([&](auto j) {
          fillZDCSideSectorHistos<i, j>(cent, zdcread); // Fill i-side sector j histograms
        });
      });

      float sumZNC = (zdcread.energySectorZNC())[0] + (zdcread.energySectorZNC())[1] + (zdcread.energySectorZNC())[2] + (zdcread.energySectorZNC())[3];
      float sumZNA = (zdcread.energySectorZNA())[0] + (zdcread.energySectorZNA())[1] + (zdcread.energySectorZNA())[2] + (zdcread.energySectorZNA())[3];
      float sumZPC = (zdcread.energySectorZPC())[0] + (zdcread.energySectorZPC())[1] + (zdcread.energySectorZPC())[2] + (zdcread.energySectorZPC())[3];
      float sumZPA = (zdcread.energySectorZPA())[0] + (zdcread.energySectorZPA())[1] + (zdcread.energySectorZPA())[2] + (zdcread.energySectorZPA())[3];

      float alphaZN = (sumZNA - sumZNC) / (sumZNA + sumZNC);
      float alphaZP = (sumZPA - sumZPC) / (sumZPA + sumZPC);

      float alphaZNAmplitude = (zdcread.amplitudeZNA() - zdcread.amplitudeZNC()) / (zdcread.amplitudeZNA() + zdcread.amplitudeZNC());
      float alphaZPAmplitude = (zdcread.amplitudeZPA() - zdcread.amplitudeZPC()) / (zdcread.amplitudeZPA() + zdcread.amplitudeZPC());

      histos.fill(HIST("CentvsDiffZNSignal"), cent, (sumZNA + sumZNC) - (zdcread.energyCommonZNA() + zdcread.energyCommonZNC()));
      histos.fill(HIST("CentvsDiffZPSignal"), cent, (sumZPA + sumZPC) - (zdcread.energyCommonZPA() + zdcread.energyCommonZPC()));
      histos.fill(HIST("CentvsZNSignalSum"), cent, sumZNA + sumZNC);
      histos.fill(HIST("CentvsZPSignalSum"), cent, sumZPA + sumZPC);
      histos.fill(HIST("CentvsZNSignalCommon"), cent, (zdcread.energyCommonZNA() + zdcread.energyCommonZNC()));
      histos.fill(HIST("CentvsZPSignalCommon"), cent, (zdcread.energyCommonZPA() + zdcread.energyCommonZPC()));
      histos.fill(HIST("CentvsAlphaZN"), cent, alphaZN);
      histos.fill(HIST("CentvsAlphaZP"), cent, alphaZP);
      histos.fill(HIST("CentvsAlphaZNcommon"), cent, (zdcread.energyCommonZNA() - zdcread.energyCommonZNC()) / (zdcread.energyCommonZNA() + zdcread.energyCommonZNC()));
      histos.fill(HIST("CentvsAlphaZPcommon"), cent, (zdcread.energyCommonZPA() - zdcread.energyCommonZPC()) / (zdcread.energyCommonZPA() + zdcread.energyCommonZPC()));
      histos.fill(HIST("CentvsAlphaZNAmplitude"), cent, alphaZNAmplitude);
      histos.fill(HIST("CentvsAlphaZPAmplitude"), cent, alphaZPAmplitude);

      histos.fill(HIST("CentvsZNAvsZNC"), cent, sumZNA, sumZNC);
      histos.fill(HIST("CentvsZNAvsZPA"), cent, sumZNA, sumZPA);
      histos.fill(HIST("CentvsZNAvsZPC"), cent, sumZNA, sumZPC);
      histos.fill(HIST("CentvsZPAvsZNC"), cent, sumZPA, sumZNC);
      histos.fill(HIST("CentvsZPAvsZPC"), cent, sumZPA, sumZPC);
      histos.fill(HIST("CentvsZNCvsZPC"), cent, sumZNC, sumZPC);
      histos.fill(HIST("CentvsZNvsZP"), cent, sumZNA + sumZNC, sumZPA + sumZPC);
    }
  }
  PROCESS_SWITCH(NeutronProtonCorrZdc, processRun3, "Process analysis for Run 3 data", true);

  void processRun2(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Run2MatchedSparse, CentralitiesRun2>>::iterator const& collision, aod::BCsWithTimestamps const&, aod::Zdcs const&)
  {
    histos.fill(HIST("eventCounter"), EventCounter::kNoSelection);
    if (!collision.alias_bit(kINT7)) {
      return;
    }
    histos.fill(HIST("eventCounter"), EventCounter::kSel8);
    if (collision.centRun2V0M() > cfgMaxCentrality) {
      return;
    }
    histos.fill(HIST("eventCounter"), EventCounter::kMaxCentralitySelection);

    if (collision.has_zdc()) {
      const auto& zdcread = collision.zdc();
      const auto cent = collision.centRun2V0M();

      histos.fill(HIST("eventCounter"), EventCounter::kZDCSelection);
      histos.fill(HIST("CentralityPercentile"), cent);

      static_for<0, 1>([&](auto i) {
        fillZDCSideCommonHistos<i>(cent, zdcread); // Fill i-side common channels
        static_for<0, 3>([&](auto j) {
          fillZDCSideSectorHistos<i, j>(cent, zdcread); // Fill i-side sector j
        });
      });

      float sumZNC = (zdcread.energySectorZNC())[0] + (zdcread.energySectorZNC())[1] + (zdcread.energySectorZNC())[2] + (zdcread.energySectorZNC())[3];
      float sumZNA = (zdcread.energySectorZNA())[0] + (zdcread.energySectorZNA())[1] + (zdcread.energySectorZNA())[2] + (zdcread.energySectorZNA())[3];
      float sumZPC = (zdcread.energySectorZPC())[0] + (zdcread.energySectorZPC())[1] + (zdcread.energySectorZPC())[2] + (zdcread.energySectorZPC())[3];
      float sumZPA = (zdcread.energySectorZPA())[0] + (zdcread.energySectorZPA())[1] + (zdcread.energySectorZPA())[2] + (zdcread.energySectorZPA())[3];

      float alphaZN = (sumZNA - sumZNC) / (sumZNA + sumZNC);
      float alphaZP = (sumZPA - sumZPC) / (sumZPA + sumZPC);

      histos.fill(HIST("CentvsDiffZNSignal"), cent, (sumZNA + sumZNC) - (zdcread.energyCommonZNA() + zdcread.energyCommonZNC()));
      histos.fill(HIST("CentvsDiffZPSignal"), cent, (sumZPA + sumZPC) - (zdcread.energyCommonZPA() + zdcread.energyCommonZPC()));
      histos.fill(HIST("CentvsZNSignalSum"), cent, sumZNA + sumZNC);
      histos.fill(HIST("CentvsZPSignalSum"), cent, sumZPA + sumZPC);
      histos.fill(HIST("CentvsZNSignalCommon"), cent, (zdcread.energyCommonZNA() + zdcread.energyCommonZNC()));
      histos.fill(HIST("CentvsZPSignalCommon"), cent, (zdcread.energyCommonZPA() + zdcread.energyCommonZPC()));
      histos.fill(HIST("CentvsAlphaZN"), cent, alphaZN);
      histos.fill(HIST("CentvsAlphaZP"), cent, alphaZP);
      histos.fill(HIST("CentvsAlphaZNcommon"), cent, (zdcread.energyCommonZNA() - zdcread.energyCommonZNC()) / (zdcread.energyCommonZNA() + zdcread.energyCommonZNC()));
      histos.fill(HIST("CentvsAlphaZPcommon"), cent, (zdcread.energyCommonZPA() - zdcread.energyCommonZPC()) / (zdcread.energyCommonZPA() + zdcread.energyCommonZPC()));

      histos.fill(HIST("CentvsZNAvsZNC"), cent, sumZNA, sumZNC);
      histos.fill(HIST("CentvsZNAvsZPA"), cent, sumZNA, sumZPA);
      histos.fill(HIST("CentvsZNAvsZPC"), cent, sumZNA, sumZPC);
      histos.fill(HIST("CentvsZPAvsZNC"), cent, sumZPA, sumZNC);
      histos.fill(HIST("CentvsZPAvsZPC"), cent, sumZPA, sumZPC);
      histos.fill(HIST("CentvsZNCvsZPC"), cent, sumZNC, sumZPC);
      histos.fill(HIST("CentvsZNvsZP"), cent, sumZNA + sumZNC, sumZPA + sumZPC);
    }
  }
  PROCESS_SWITCH(NeutronProtonCorrZdc, processRun2, "Process analysis for Run 2 converted data", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<NeutronProtonCorrZdc>(cfgc)};
}
