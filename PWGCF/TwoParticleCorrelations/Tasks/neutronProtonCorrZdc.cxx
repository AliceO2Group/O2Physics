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
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"
#include "Framework/StaticFor.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

enum EventCounter { kNoSelection = 0,
                    kQualitySelection = 1,
                    kMaxCentralitySelection = 2,
                    kZDCSelection = 3 };

struct NeutronProtonCorrZdc {
  // Histogram registry: an object to hold your histograms
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<int> cfgNBinsZN{"cfgNBinsZN", 100, "N bins for ZNA and ZNC"};
  Configurable<int> cfgNBinsZP{"cfgNBinsZP", 100, "N bins for ZPA and ZPC"};
  Configurable<double> cfgZNmin{"cfgZNmin", -10, "Minimum value for ZN signal"};
  Configurable<double> cfgZNmax{"cfgZNmax", 350, "Maximum value for ZN signal"};
  Configurable<double> cfgZPmin{"cfgZPmin", -10, "Minimum value for ZP signal"};
  Configurable<double> cfgZPmax{"cfgZPmax", 200, "Maximum value for ZP signal"};
  Configurable<double> cfgDiffZmin{"cfgDiffZmin", -30, "Minimum value for the diffZ signal"};
  Configurable<double> cfgDiffZmax{"cfgDiffZmax", 50, "Maximum value for the diffZ signal"};
  Configurable<int> cfgNBinsAlpha{"cfgNBinsAlpha", 100, "Number of bins for ZDC asymmetry"};
  Configurable<double> cfgAlphaZmin{"cfgAlphaZmin", -1, "Minimum value for ZDC asymmetry"};
  Configurable<double> cfgAlphaZmax{"cfgAlphaZmax", 1, "Maximum value for ZDC asymmetry"};
  Configurable<double> cfgMaxCentrality{"cfgMaxCentrality", 80, "Maximum collision centrality"};

  ConfigurableAxis cfgAxisCent{"cfgAxisCent", {VARIABLE_WIDTH, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0, 31.0, 32.0, 33.0, 34.0, 35.0, 36.0, 37.0, 38.0, 39.0, 40.0, 41.0, 42.0, 43.0, 44.0, 45.0, 46.0, 47.0, 48.0, 49.0, 50.0, 51.0, 52.0, 53.0, 54.0, 55.0, 56.0, 57.0, 58.0, 59.0, 60.0, 61.0, 62.0, 63.0, 64.0, 65.0, 66.0, 67.0, 68.0, 69.0, 70.0, 71.0, 72.0, 73.0, 74.0, 75.0, 76.0, 77.0, 78.0, 79.0, 80.0, 81.0, 82.0, 83.0, 84.0, 85.0, 86.0, 87.0, 88.0, 89.0, 90.0, 91.0, 92.0, 93.0, 94.0, 95.0, 96.0, 97.0, 98.0, 99.0, 100.0}, "Centrality [%]"};

  Filter collisionVtxZ = nabs(aod::collision::posZ) < 10.f;

  using CentralitiesRun3 = aod::CentFT0Cs;
  using CentralitiesRun2 = aod::CentRun2V0Ms;
  using BCsRun3 = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;

  void init(InitContext const&)
  {
    // define axes you want to use
    const AxisSpec axisCounter{4, -0.5, 3.5, ""};
    const AxisSpec axisZNSectorSignal{cfgNBinsZN, cfgZNmin, cfgZNmax / 3.};
    const AxisSpec axisZPSectorSignal{cfgNBinsZP, cfgZPmin, cfgZPmax / 3.};
    const AxisSpec axisZNASignal{cfgNBinsZN, cfgZNmin, cfgZNmax, "ZNA (a.u.)"};
    const AxisSpec axisZNCSignal{cfgNBinsZN, cfgZNmin, cfgZNmax, "ZNC (a.u.)"};
    const AxisSpec axisZPASignal{cfgNBinsZP, cfgZPmin, cfgZPmax, "ZPA (a.u.)"};
    const AxisSpec axisZPCSignal{cfgNBinsZP, cfgZPmin, cfgZPmax, "ZPC (a.u.)"};
    const AxisSpec axisZNSignal{2 * cfgNBinsZN, cfgZNmin, 1.5 * cfgZNmax, "ZN (a.u.)"};
    const AxisSpec axisZPSignal{2 * cfgNBinsZP, cfgZPmin, 1.5 * cfgZPmax, "ZP (a.u.)"};
    const AxisSpec axisAlphaZ{cfgNBinsAlpha, cfgAlphaZmin, cfgAlphaZmax, "#alpha_{spec}"};
    const AxisSpec axisZDiffSignal{cfgNBinsZN, cfgDiffZmin, cfgDiffZmax, "#Delta E"};

    HistogramConfigSpec defaultZNSectorHist({HistType::kTH2F, {cfgAxisCent, axisZNSectorSignal}});
    HistogramConfigSpec defaultZPSectorHist({HistType::kTH2F, {cfgAxisCent, axisZPSectorSignal}});
    HistogramConfigSpec defaultZDCDiffHist({HistType::kTH2F, {cfgAxisCent, axisZDiffSignal}});

    // create histograms
    histos.add("eventCounter", "eventCounter", kTH1F, {axisCounter});
    histos.add("CentralityPercentile", "CentralityPercentile", kTH1F, {cfgAxisCent});

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

    // Cloning the folder
    histos.addClone("ASide/", "CSide/");

    histos.add("CentvsZNSignalCommon", "CentvsZNSignalCommon", kTH2F, {cfgAxisCent, axisZNSignal});
    histos.add("CentvsZNSignalSum", "CentvsZNSignalSum", kTH2F, {cfgAxisCent, axisZNSignal});
    histos.add("CentvsZPSignalCommon", "CentvsZPSignalCommon", kTH2F, {cfgAxisCent, axisZPSignal});
    histos.add("CentvsZPSignalSum", "CentvsZPSignalSum", kTH2F, {cfgAxisCent, axisZPSignal});
    histos.add("CentvsAlphaZN", "CentvsAlphaZN", kTH2F, {cfgAxisCent, axisAlphaZ});
    histos.add("CentvsAlphaZP", "CentvsAlphaZP", kTH2F, {cfgAxisCent, axisAlphaZ});
    histos.add("CentvsDiffZNSignal", "CentvsDiffZNSignal", defaultZDCDiffHist);
    histos.add("CentvsDiffZPSignal", "CentvsDiffZPSignal", defaultZDCDiffHist);
  }
  template <int side, typename Z>
  void fillZDCHistos(const float centr, const Z& zdc)
  {
    static constexpr std::string_view SubDir[] = {"ASide/", "CSide/"};

    std::array<std::array<float, 4>, 2> znEnergyResponse = {zdc.energySectorZNA(), zdc.energySectorZNC()};
    std::array<std::array<float, 4>, 2> zpEnergyResponse = {zdc.energySectorZPA(), zdc.energySectorZPC()};
    std::array<float, 2> znEnergyResponseCommon = {zdc.energyCommonZNA(), zdc.energyCommonZNC()};
    std::array<float, 2> zpEnergyResponseCommon = {zdc.energyCommonZPA(), zdc.energyCommonZPC()};

    // Fill Neutron ZDC historgrams
    histos.fill(HIST(SubDir[side]) + HIST("CentvsZNSector0Signal"), centr, znEnergyResponse[side][0]);
    histos.fill(HIST(SubDir[side]) + HIST("CentvsZNSector1Signal"), centr, znEnergyResponse[side][1]);
    histos.fill(HIST(SubDir[side]) + HIST("CentvsZNSector2Signal"), centr, znEnergyResponse[side][2]);
    histos.fill(HIST(SubDir[side]) + HIST("CentvsZNSector3Signal"), centr, znEnergyResponse[side][3]);

    float sumZN = znEnergyResponse[side][0] + znEnergyResponse[side][1] + znEnergyResponse[side][2] + znEnergyResponse[side][3];

    histos.fill(HIST(SubDir[side]) + HIST("CentvsZNSignalSum"), centr, sumZN);
    histos.fill(HIST(SubDir[side]) + HIST("CentvsZNSignalCommon"), centr, znEnergyResponseCommon[side]);
    histos.fill(HIST(SubDir[side]) + HIST("CentvsdiffZNSignal"), centr, sumZN - znEnergyResponseCommon[side]);

    // Fill Proton ZDC histograms
    histos.fill(HIST(SubDir[side]) + HIST("CentvsZPSector0Signal"), centr, zpEnergyResponse[side][0]);
    histos.fill(HIST(SubDir[side]) + HIST("CentvsZPSector1Signal"), centr, zpEnergyResponse[side][1]);
    histos.fill(HIST(SubDir[side]) + HIST("CentvsZPSector2Signal"), centr, zpEnergyResponse[side][2]);
    histos.fill(HIST(SubDir[side]) + HIST("CentvsZPSector3Signal"), centr, zpEnergyResponse[side][3]);

    float sumZP = zpEnergyResponse[side][0] + zpEnergyResponse[side][1] + zpEnergyResponse[side][2] + zpEnergyResponse[side][3];

    histos.fill(HIST(SubDir[side]) + HIST("CentvsZPSignalSum"), centr, sumZP);
    histos.fill(HIST(SubDir[side]) + HIST("CentvsZPSignalCommon"), centr, zpEnergyResponseCommon[side]);
    histos.fill(HIST(SubDir[side]) + HIST("CentvsdiffZPSignal"), centr, sumZP - zpEnergyResponseCommon[side]);
  }

  void processRun3(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, CentralitiesRun3>>::iterator const& collision, BCsRun3 const&, aod::Zdcs const&)
  {
    histos.fill(HIST("eventCounter"), EventCounter::kNoSelection);
    if (!collision.sel8()) {
      return;
    }
    histos.fill(HIST("eventCounter"), EventCounter::kQualitySelection);
    if (collision.centFT0C() > cfgMaxCentrality) {
      return;
    }
    histos.fill(HIST("eventCounter"), EventCounter::kMaxCentralitySelection);
    const auto& foundBC = collision.foundBC_as<BCsRun3>();
    if (foundBC.has_zdc()) {
      const auto& zdcread = foundBC.zdc();
      const auto cent = collision.centFT0C();

      histos.fill(HIST("eventCounter"), EventCounter::kZDCSelection);
      histos.fill(HIST("CentralityPercentile"), cent);

      fillZDCHistos<0>(cent, zdcread); // Fill A-side
      fillZDCHistos<1>(cent, zdcread); // Fill C-side

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
    }
  }
  PROCESS_SWITCH(NeutronProtonCorrZdc, processRun3, "Process analysis for Run 3 data", true);

  void processRun2(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Run2MatchedSparse, CentralitiesRun2>>::iterator const& collision, aod::BCsWithTimestamps const&, aod::Zdcs const&)
  {
    histos.fill(HIST("eventCounter"), EventCounter::kNoSelection);
    if (!collision.alias_bit(kINT7)) {
      return;
    }
    histos.fill(HIST("eventCounter"), EventCounter::kQualitySelection);
    if (collision.centRun2V0M() > cfgMaxCentrality) {
      return;
    }
    histos.fill(HIST("eventCounter"), EventCounter::kMaxCentralitySelection);

    if (collision.has_zdc()) {
      const auto& zdcread = collision.zdc();
      const auto cent = collision.centRun2V0M();

      histos.fill(HIST("eventCounter"), EventCounter::kZDCSelection);
      histos.fill(HIST("CentralityPercentile"), cent);

      fillZDCHistos<0>(cent, zdcread); // Fill A-side
      fillZDCHistos<1>(cent, zdcread); // Fill C-side

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
    }
  }
  PROCESS_SWITCH(NeutronProtonCorrZdc, processRun2, "Process analysis for Run 2 converted data", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<NeutronProtonCorrZdc>(cfgc)};
}
