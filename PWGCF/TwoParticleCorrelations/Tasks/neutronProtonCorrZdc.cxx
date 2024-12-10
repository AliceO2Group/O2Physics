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

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"

#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct NeutronProtonCorrZdc {
  // Histogram registry: an object to hold your histograms
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<int> cfgNBinsZN{"cfgNBinsZN", 100, "N bins for ZNA and ZNC"};
  Configurable<int> cfgNBinsZP{"cfgNBinsZP", 100, "N bins for ZPA and ZPC"};
  Configurable<double> cfgZNmin{"cfgZNmin", -10, "Minimum value for ZN signal"};
  Configurable<double> cfgZNmax{"cfgZNmax", 350, "Maximum value for ZN signal"};
  Configurable<double> cfgZPmin{"cfgZPmin", -10, "Minimum value for ZP signal"};
  Configurable<double> cfgZPmax{"cfgZPmax", 200, "Maximum value for ZP signal"};
  Configurable<double> cfgDiffZmin{"cfgDiffZmin", -50, "Minimum value for the diffZ signal"};
  Configurable<double> cfgDiffZmax{"cfgDiffZmax", 50, "Maximum value for the diffZ signal"};
  Configurable<int> cfgNBinsAlpha{"cfgNBinsAlpha", 100, "Number of bins for ZDC asymmetry"};
  Configurable<double> cfgAlphaZmin{"cfgAlphaZmin", -1, "Minimum value for ZDC asymmetry"};
  Configurable<double> cfgAlphaZmax{"cfgAlphaZmax", 1, "Maximum value for ZDC asymmetry"};
  Configurable<bool> cfgProcessRun2{"cfgProcessRun2", false, "Analyse Run 2 converted data"};
  // Configurable<int> cfgCentralityEstimator{"cfgCentralityEstimator", 0, "Choice of centrality estimator"};//0 for FTOC, 1 for FTOA, 2 for FTOM, 3 for FVOA. //To be included at a later stage

  ConfigurableAxis cfgAxisCent{"cfgAxisCent", {VARIABLE_WIDTH, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0, 31.0, 32.0, 33.0, 34.0, 35.0, 36.0, 37.0, 38.0, 39.0, 40.0, 41.0, 42.0, 43.0, 44.0, 45.0, 46.0, 47.0, 48.0, 49.0, 50.0, 51.0, 52.0, 53.0, 54.0, 55.0, 56.0, 57.0, 58.0, 59.0, 60.0, 61.0, 62.0, 63.0, 64.0, 65.0, 66.0, 67.0, 68.0, 69.0, 70.0, 71.0, 72.0, 73.0, 74.0, 75.0, 76.0, 77.0, 78.0, 79.0, 80.0, 81.0, 82.0, 83.0, 84.0, 85.0, 86.0, 87.0, 88.0, 89.0, 90.0, 91.0, 92.0, 93.0, 94.0, 95.0, 96.0, 97.0, 98.0, 99.0, 100.0}, "Centrality [%]"};

  Filter collisionVtxZ = nabs(aod::collision::posZ) < 10.f;

  using CentralitiesRun3 = aod::CentFT0Cs;
  using CentralitiesRun2 = aod::CentRun2V0Ms;
  using BCsRun3 = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;

  void init(InitContext const&)
  {
    // define axes you want to use
    const AxisSpec axisCounter{4, 0, +2, ""};
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

    // create histograms
    histos.add("eventCounter", "eventCounter", kTH1F, {axisCounter});

    histos.add("ZNASector0Signal", "ZNASector0Signal", kTH1F, {axisZNSectorSignal});
    histos.add("ZNASector1Signal", "ZNASector1Signal", kTH1F, {axisZNSectorSignal});
    histos.add("ZNASector2Signal", "ZNASector2Signal", kTH1F, {axisZNSectorSignal});
    histos.add("ZNASector3Signal", "ZNASector3Signal", kTH1F, {axisZNSectorSignal});

    histos.add("ZNCSector0Signal", "ZNCSector0Signal", kTH1F, {axisZNSectorSignal});
    histos.add("ZNCSector1Signal", "ZNCSector1Signal", kTH1F, {axisZNSectorSignal});
    histos.add("ZNCSector2Signal", "ZNCSector2Signal", kTH1F, {axisZNSectorSignal});
    histos.add("ZNCSector3Signal", "ZNCSector3Signal", kTH1F, {axisZNSectorSignal});

    histos.add("ZPASector0Signal", "ZPASector0Signal", kTH1F, {axisZPSectorSignal});
    histos.add("ZPASector1Signal", "ZPASector1Signal", kTH1F, {axisZPSectorSignal});
    histos.add("ZPASector2Signal", "ZPASector2Signal", kTH1F, {axisZPSectorSignal});
    histos.add("ZPASector3Signal", "ZPASector3Signal", kTH1F, {axisZPSectorSignal});

    histos.add("ZPCSector0Signal", "ZPCSector0Signal", kTH1F, {axisZPSectorSignal});
    histos.add("ZPCSector1Signal", "ZPCSector1Signal", kTH1F, {axisZPSectorSignal});
    histos.add("ZPCSector2Signal", "ZPCSector2Signal", kTH1F, {axisZPSectorSignal});
    histos.add("ZPCSector3Signal", "ZPCSector3Signal", kTH1F, {axisZPSectorSignal});

    histos.add("ZNASignal", "ZNASignal", kTH1F, {axisZNASignal});
    histos.add("ZNCSignal", "ZNCSignal", kTH1F, {axisZNCSignal});
    histos.add("ZNSignal", "ZNSignal", kTH1F, {axisZNSignal});
    histos.add("ZPASignal", "ZPASignal", kTH1F, {axisZPASignal});
    histos.add("ZPCSignal", "ZPCSignal", kTH1F, {axisZPCSignal});
    histos.add("ZPSignal", "ZPSignal", kTH1F, {axisZPSignal});

    histos.add("alphaZN", "alphaZN", kTH1F, {axisAlphaZ});
    histos.add("alphaZP", "alphaZP", kTH1F, {axisAlphaZ});

    histos.add("diffZNASignal", "diffZNASignal", kTH1F, {axisZDiffSignal});
    histos.add("diffZNCSignal", "diffZNCSignal", kTH1F, {axisZDiffSignal});
    histos.add("diffZPASignal", "diffZPASignal", kTH1F, {axisZDiffSignal});
    histos.add("diffZPCSignal", "diffZPCSignal", kTH1F, {axisZDiffSignal});
    histos.add("diffZNSignal", "diffZNSignal", kTH1F, {axisZDiffSignal});
    histos.add("diffZPSignal", "diffZPSignal", kTH1F, {axisZDiffSignal});

    histos.add("CentralityPercentile", "CentralityPercentile", kTH1F, {cfgAxisCent});

    histos.add("CentvsZNASignal", "FT0CvsZNASignal", kTH2F, {cfgAxisCent, axisZNASignal});
    histos.add("CentvsZNCSignal", "FT0CvsZNCSignal", kTH2F, {cfgAxisCent, axisZNCSignal});
    histos.add("CentvsZPASignal", "FT0CvsZPASignal", kTH2F, {cfgAxisCent, axisZPASignal});
    histos.add("CentvsZPCSignal", "FT0CvsZPCSignal", kTH2F, {cfgAxisCent, axisZPCSignal});
    histos.add("CentvsZNSignal", "FT0CvsZNSignal", kTH2F, {cfgAxisCent, axisZNSignal});
    histos.add("CentvsZPSignal", "FT0CvsZPSignal", kTH2F, {cfgAxisCent, axisZPSignal});
  }

  void processRun3(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, CentralitiesRun3>>::iterator const& collision, BCsRun3 const&, aod::Zdcs const&)
  {
    if (!collision.sel8()) {
      histos.fill(HIST("eventCounter"), 0.25);
      return;
    }
    const auto& foundBC = collision.foundBC_as<BCsRun3>();
    if (foundBC.has_zdc()) {
      const auto& zdcread = foundBC.zdc();
      const auto cent = collision.centFT0C();

      histos.fill(HIST("eventCounter"), 1.25);
      histos.fill(HIST("ZNASignal"), zdcread.energyCommonZNA());
      histos.fill(HIST("ZNCSignal"), zdcread.energyCommonZNC());
      histos.fill(HIST("ZPASignal"), zdcread.energyCommonZPA());
      histos.fill(HIST("ZPCSignal"), zdcread.energyCommonZPC());
      histos.fill(HIST("ZNASector0Signal"), (zdcread.energySectorZNA())[0]);
      histos.fill(HIST("ZNASector1Signal"), (zdcread.energySectorZNA())[1]);
      histos.fill(HIST("ZNASector2Signal"), (zdcread.energySectorZNA())[2]);
      histos.fill(HIST("ZNASector3Signal"), (zdcread.energySectorZNA())[3]);
      histos.fill(HIST("ZNCSector0Signal"), (zdcread.energySectorZNC())[0]);
      histos.fill(HIST("ZNCSector1Signal"), (zdcread.energySectorZNC())[1]);
      histos.fill(HIST("ZNCSector2Signal"), (zdcread.energySectorZNC())[2]);
      histos.fill(HIST("ZNCSector3Signal"), (zdcread.energySectorZNC())[3]);
      histos.fill(HIST("ZPASector0Signal"), (zdcread.energySectorZPA())[0]);
      histos.fill(HIST("ZPASector1Signal"), (zdcread.energySectorZPA())[1]);
      histos.fill(HIST("ZPASector2Signal"), (zdcread.energySectorZPA())[2]);
      histos.fill(HIST("ZPASector3Signal"), (zdcread.energySectorZPA())[3]);
      histos.fill(HIST("ZPCSector0Signal"), (zdcread.energySectorZPC())[0]);
      histos.fill(HIST("ZPCSector1Signal"), (zdcread.energySectorZPC())[1]);
      histos.fill(HIST("ZPCSector2Signal"), (zdcread.energySectorZPC())[2]);
      histos.fill(HIST("ZPCSector3Signal"), (zdcread.energySectorZPC())[3]);

      float sumZNC = (zdcread.energySectorZNC())[0] + (zdcread.energySectorZNC())[1] + (zdcread.energySectorZNC())[2] + (zdcread.energySectorZNC())[3];
      float sumZNA = (zdcread.energySectorZNA())[0] + (zdcread.energySectorZNA())[1] + (zdcread.energySectorZNA())[2] + (zdcread.energySectorZNA())[3];
      float sumZPC = (zdcread.energySectorZPC())[0] + (zdcread.energySectorZPC())[1] + (zdcread.energySectorZPC())[2] + (zdcread.energySectorZPC())[3];
      float sumZPA = (zdcread.energySectorZPA())[0] + (zdcread.energySectorZPA())[1] + (zdcread.energySectorZPA())[2] + (zdcread.energySectorZPA())[3];
      float alphaZN = (sumZNA - sumZNC) / (sumZNA + sumZNC);
      float alphaZP = (sumZPA - sumZPC) / (sumZPA + sumZPC);

      histos.fill(HIST("alphaZN"), alphaZN);
      histos.fill(HIST("alphaZP"), alphaZP);

      histos.fill(HIST("diffZNASignal"), sumZNA - zdcread.energyCommonZNA());
      histos.fill(HIST("diffZNCSignal"), sumZNC - zdcread.energyCommonZNC());
      histos.fill(HIST("diffZPASignal"), sumZPA - zdcread.energyCommonZPA());
      histos.fill(HIST("diffZPCSignal"), sumZPC - zdcread.energyCommonZPC());
      histos.fill(HIST("ZNSignal"), sumZNA + sumZNC);
      histos.fill(HIST("ZPSignal"), sumZPA + sumZPC);
      histos.fill(HIST("diffZNSignal"), (sumZNA + sumZNC) - (zdcread.energyCommonZNA() + zdcread.energyCommonZNC()));
      histos.fill(HIST("diffZPSignal"), (sumZPA + sumZPC) - (zdcread.energyCommonZPA() + zdcread.energyCommonZPC()));
      histos.fill(HIST("CentralityPercentile"), cent);
      histos.fill(HIST("CentvsZNASignal"), cent, sumZNA);
      histos.fill(HIST("CentvsZNCSignal"), cent, sumZNC);
      histos.fill(HIST("CentvsZPASignal"), cent, sumZPA);
      histos.fill(HIST("CentvsZPCSignal"), cent, sumZPC);
      histos.fill(HIST("CentvsZNSignal"), cent, sumZNA + sumZNC);
      histos.fill(HIST("CentvsZPSignal"), cent, sumZPA + sumZPC);
    }
  }
  PROCESS_SWITCH(NeutronProtonCorrZdc, processRun3, "Process analysis for Run 3 data", true);

  void processRun2(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Run2MatchedSparse, CentralitiesRun2>>::iterator const& collision, aod::BCsWithTimestamps const&, aod::Zdcs const&)
  {
    if (!collision.alias_bit(kINT7)) {
      histos.fill(HIST("eventCounter"), 0.25);
      return;
    }

    if (collision.has_zdc()) {
      const auto& zdcread = collision.zdc();
      const auto cent = collision.centRun2V0M();
      histos.fill(HIST("eventCounter"), 1.25);
      histos.fill(HIST("ZNASignal"), zdcread.energyCommonZNA());
      histos.fill(HIST("ZNCSignal"), zdcread.energyCommonZNC());
      histos.fill(HIST("ZPASignal"), zdcread.energyCommonZPA());
      histos.fill(HIST("ZPCSignal"), zdcread.energyCommonZPC());
      histos.fill(HIST("ZNASector0Signal"), (zdcread.energySectorZNA())[0]);
      histos.fill(HIST("ZNASector1Signal"), (zdcread.energySectorZNA())[1]);
      histos.fill(HIST("ZNASector2Signal"), (zdcread.energySectorZNA())[2]);
      histos.fill(HIST("ZNASector3Signal"), (zdcread.energySectorZNA())[3]);
      histos.fill(HIST("ZNCSector0Signal"), (zdcread.energySectorZNC())[0]);
      histos.fill(HIST("ZNCSector1Signal"), (zdcread.energySectorZNC())[1]);
      histos.fill(HIST("ZNCSector2Signal"), (zdcread.energySectorZNC())[2]);
      histos.fill(HIST("ZNCSector3Signal"), (zdcread.energySectorZNC())[3]);
      histos.fill(HIST("ZPASector0Signal"), (zdcread.energySectorZPA())[0]);
      histos.fill(HIST("ZPASector1Signal"), (zdcread.energySectorZPA())[1]);
      histos.fill(HIST("ZPASector2Signal"), (zdcread.energySectorZPA())[2]);
      histos.fill(HIST("ZPASector3Signal"), (zdcread.energySectorZPA())[3]);
      histos.fill(HIST("ZPCSector0Signal"), (zdcread.energySectorZPC())[0]);
      histos.fill(HIST("ZPCSector1Signal"), (zdcread.energySectorZPC())[1]);
      histos.fill(HIST("ZPCSector2Signal"), (zdcread.energySectorZPC())[2]);
      histos.fill(HIST("ZPCSector3Signal"), (zdcread.energySectorZPC())[3]);
      float sumZNC = (zdcread.energySectorZNC())[0] + (zdcread.energySectorZNC())[1] + (zdcread.energySectorZNC())[2] + (zdcread.energySectorZNC())[3];
      float sumZNA = (zdcread.energySectorZNA())[0] + (zdcread.energySectorZNA())[1] + (zdcread.energySectorZNA())[2] + (zdcread.energySectorZNA())[3];
      float sumZPC = (zdcread.energySectorZPC())[0] + (zdcread.energySectorZPC())[1] + (zdcread.energySectorZPC())[2] + (zdcread.energySectorZPC())[3];
      float sumZPA = (zdcread.energySectorZPA())[0] + (zdcread.energySectorZPA())[1] + (zdcread.energySectorZPA())[2] + (zdcread.energySectorZPA())[3];
      float alphaZN = (sumZNA - sumZNC) / (sumZNA + sumZNC);
      float alphaZP = (sumZPA - sumZPC) / (sumZPA + sumZPC);
      histos.fill(HIST("alphaZN"), alphaZN);
      histos.fill(HIST("alphaZP"), alphaZP);

      histos.fill(HIST("diffZNASignal"), sumZNA - zdcread.energyCommonZNA());
      histos.fill(HIST("diffZNCSignal"), sumZNC - zdcread.energyCommonZNC());
      histos.fill(HIST("diffZPASignal"), sumZPA - zdcread.energyCommonZPA());
      histos.fill(HIST("diffZPCSignal"), sumZPC - zdcread.energyCommonZPC());
      histos.fill(HIST("ZNSignal"), sumZNA + sumZNC);
      histos.fill(HIST("ZPSignal"), sumZPA + sumZPC);
      histos.fill(HIST("diffZNSignal"), (sumZNA + sumZNC) - (zdcread.energyCommonZNA() + zdcread.energyCommonZNC()));
      histos.fill(HIST("diffZPSignal"), (sumZPA + sumZPC) - (zdcread.energyCommonZPA() + zdcread.energyCommonZPC()));
      histos.fill(HIST("CentralityPercentile"), cent);
      histos.fill(HIST("CentvsZNASignal"), cent, sumZNA);
      histos.fill(HIST("CentvsZNCSignal"), cent, sumZNC);
      histos.fill(HIST("CentvsZPASignal"), cent, sumZPA);
      histos.fill(HIST("CentvsZPCSignal"), cent, sumZPC);
      histos.fill(HIST("CentvsZNSignal"), cent, sumZNA + sumZNC);
      histos.fill(HIST("CentvsZPSignal"), cent, sumZPA + sumZPC);
    }
  }
  PROCESS_SWITCH(NeutronProtonCorrZdc, processRun2, "Process analysis for Run 2 converted data", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<NeutronProtonCorrZdc>(cfgc)};
}
