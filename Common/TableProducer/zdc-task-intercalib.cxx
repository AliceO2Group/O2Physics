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
/// \brief Task for ZDC tower inter-calibration
/// \author chiara.oppedisano@cern.ch
// Minimal example to run this task:
// export OPTIONS="-b --configuration json://config.json --aod-file AO2D.root"
//  o2-analysis-timestamp ${OPTIONS} |
//  o2-analysis-event-selection ${OPTIONS} |
//  o2-analysis-track-propagation ${OPTIONS} |
//  o2-analysis-mm-zdc-task-intercalib ${OPTIONS}

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/TriggerAliases.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/ZDCInterCalib.h"

#include "TH1F.h"
#include "TH2F.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::evsel;

constexpr double kVeryNegative = -1.e12;

using BCsRun3 = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;
using ColEvSels = soa::Join<aod::Collisions, aod::EvSels>;

struct zdcInterCalib {

  Produces<aod::ZDCInterCalib> zTab;

  // Configurable parameters
  //
  Configurable<int> nBins{"nBins", 400, "n bins"};
  Configurable<float> MaxZN{"MaxZN", 399.5, "Max ZN signal"};
  //
  HistogramRegistry registry{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&)
  {
    registry.add("ZNApmc", "ZNApmc; ZNA PMC; Entries", {HistType::kTH1F, {{nBins, -0.5, MaxZN}}});
    registry.add("ZNCpmc", "ZNCpmc; ZNC PMC; Entries", {HistType::kTH1F, {{nBins, -0.5, MaxZN}}});
    registry.add("ZNApm1", "ZNApm1; ZNA PM1; Entries", {HistType::kTH1F, {{nBins, -0.5, MaxZN}}});
    registry.add("ZNApm2", "ZNApm2; ZNA PM2; Entries", {HistType::kTH1F, {{nBins, -0.5, MaxZN}}});
    registry.add("ZNApm3", "ZNApm3; ZNA PM3; Entries", {HistType::kTH1F, {{nBins, -0.5, MaxZN}}});
    registry.add("ZNApm4", "ZNApm4; ZNA PM4; Entries", {HistType::kTH1F, {{nBins, -0.5, MaxZN}}});
    registry.add("ZNCpm1", "ZNCpm1; ZNC PM1; Entries", {HistType::kTH1F, {{nBins, -0.5, MaxZN}}});
    registry.add("ZNCpm2", "ZNCpm2; ZNC PM2; Entries", {HistType::kTH1F, {{nBins, -0.5, MaxZN}}});
    registry.add("ZNCpm3", "ZNCpm3; ZNC PM3; Entries", {HistType::kTH1F, {{nBins, -0.5, MaxZN}}});
    registry.add("ZNCpm4", "ZNCpm4; ZNC PM4; Entries", {HistType::kTH1F, {{nBins, -0.5, MaxZN}}});
    registry.add("ZNAsumq", "ZNAsumq; ZNA uncalib. sum PMQ; Entries", {HistType::kTH1F, {{nBins, -0.5, MaxZN}}});
    registry.add("ZNCsumq", "ZNCsumq; ZNC uncalib. sum PMQ; Entries", {HistType::kTH1F, {{nBins, -0.5, MaxZN}}});
  }

  void process(ColEvSels const& cols, BCsRun3 const& /*bcs*/, aod::Zdcs const& /*zdcs*/)
  {
    // collision-based event selection
    for (auto& collision : cols) {
      const auto& foundBC = collision.foundBC_as<BCsRun3>();
      if (foundBC.has_zdc()) {
        const auto& zdc = foundBC.zdc();

        // To assure that ZN have a genuine signal (tagged by the relative TDC)
        // we can check that the amplitude is >0 or that ADC is NOT very negative (-inf)
        // If this is not the case, signals are set to kVeryNegative values
        double pmcZNC = zdc.energyCommonZNC();
        double pmcZNA = zdc.energyCommonZNA();
        //
        // double tdcZNC = zdc.zdc.timeZNC();
        // double tdcZNA = zdc.zdc.timeZNA();
        //
        bool isZNChit = true, isZNAhit = true;
        if (pmcZNC < kVeryNegative) {
          pmcZNC = kVeryNegative;
          isZNChit = false;
        }
        if (pmcZNA < kVeryNegative) {
          pmcZNA = kVeryNegative;
          isZNAhit = false;
        }
        //
        double sumZNC = 0;
        double sumZNA = 0;
        double pmqZNC[4] = {
          0,
          0,
          0,
          0,
        };
        double pmqZNA[4] = {
          0,
          0,
          0,
          0,
        };
        //
        if (isZNChit) {
          for (int it = 0; it < 4; it++) {
            pmqZNC[it] = (zdc.energySectorZNC())[it];
            sumZNC += pmqZNC[it];
          }
          registry.get<TH1>(HIST("ZNCpmc"))->Fill(pmcZNC);
          registry.get<TH1>(HIST("ZNCpm1"))->Fill(pmqZNC[0]);
          registry.get<TH1>(HIST("ZNCpm2"))->Fill(pmqZNC[1]);
          registry.get<TH1>(HIST("ZNCpm3"))->Fill(pmqZNC[2]);
          registry.get<TH1>(HIST("ZNCpm4"))->Fill(pmqZNC[3]);
          registry.get<TH1>(HIST("ZNCsumq"))->Fill(sumZNC);
        }
        if (isZNAhit) {
          for (int it = 0; it < 4; it++) {
            pmqZNA[it] = (zdc.energySectorZNA())[it];
            sumZNA += pmqZNA[it];
          }
          //
          registry.get<TH1>(HIST("ZNApmc"))->Fill(pmcZNA);
          registry.get<TH1>(HIST("ZNApm1"))->Fill(pmqZNA[0]);
          registry.get<TH1>(HIST("ZNApm2"))->Fill(pmqZNA[1]);
          registry.get<TH1>(HIST("ZNApm3"))->Fill(pmqZNA[2]);
          registry.get<TH1>(HIST("ZNApm4"))->Fill(pmqZNA[3]);
          registry.get<TH1>(HIST("ZNAsumq"))->Fill(sumZNA);
        }
        if (isZNAhit || isZNChit)
          zTab(pmcZNA, pmqZNA[0], pmqZNA[1], pmqZNA[2], pmqZNA[3], pmcZNC, pmqZNC[0], pmqZNC[1], pmqZNC[2], pmqZNC[3]);
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<zdcInterCalib>(cfgc)};
}
