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
/// \brief Task for ZDC
/// \author
/// \since
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "Common/DataModel/EventSelection.h"
//#include "Common/CCDB/EventSelectionParams.h"
//#include "Common/CCDB/TriggerAliases.h"
#include "CCDB/BasicCCDBManager.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TObjArray.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// using BCsWithRun3Matchings = soa::Join<aod::BCs, aod::Timestamps, aod::Run3MatchedToBCSparse>;
// using BCsWithBcSels = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels>;
using BCsRun3 = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;

struct ZDCAnalysis {

  // Configurable for number of bins
  Configurable<int> nBins1{"nBins1", 400, "Nbins 1d ZDC histos"};
  Configurable<int> nBins2{"nBins2", 800, "Nbins 2d ZDC histos"};
  // CInfigurable maximum limit
  Configurable<float> MaxZN{"MaxZN", 8000, "Max ZN signal"};
  Configurable<float> MaxZP{"MaxZP", 6000, "Max ZP signal"};
  Configurable<float> MaxMultFV0{"MaxMultFV0", 20000, "Max FV0 signal"};
  Configurable<float> MaxMultFT0{"MaxMultFT0", 10000, "Max FT0 signal"};
  Configurable<float> MaxMultFDD{"MaxMultFDD", 10000, "Max FDD signal"};
  Configurable<float> MaxMultNTracks{"MaxMultNTracks", 1000, "Max Ntracks"};

  HistogramRegistry registry;

  void init(InitContext const&)
  {
    if (doprocessZdcAuto) { // Check if the process function for ZDCAuto is enabled
      registry.add("ZNApmc", "ZNApmc", {HistType::kTH1F, {{nBins1, -10., MaxZN}}});
      registry.add("ZPApmc", "ZPApmc", {HistType::kTH1F, {{nBins1, -10., MaxZN}}});
      registry.add("ZNCpmc", "ZNCpmc", {HistType::kTH1F, {{nBins1, -10., MaxZP}}});
      registry.add("ZPCpmc", "ZPCpmc", {HistType::kTH1F, {{nBins1, -10., MaxZP}}});
      registry.add("ZEM1", "ZEM1", {HistType::kTH1F, {{nBins1, -10., MaxZP}}});
      registry.add("ZEM2", "ZEM2", {HistType::kTH1F, {{nBins1, -10., MaxZP}}});
      registry.add("ZNvsZEM", "ZNvsZEM", {HistType::kTH2F, {{{nBins2, -10., 10000.5}, {nBins2, -10., 2. * MaxZN}}}});
      registry.add("ZNAvsZNC", "ZNAvsZNC", {HistType::kTH2F, {{{nBins2, -10., MaxZN}, {nBins2, -10., MaxZN}}}});
      registry.add("ZPAvsZPC", "ZPAvsZPC", {HistType::kTH2F, {{{nBins1, -10., MaxZP}, {nBins2, -10., MaxZP}}}});
      registry.add("ZNAvsZPA", "ZNAvsZPA", {HistType::kTH2F, {{{nBins2, -10., MaxZP}, {nBins2, -10., MaxZN}}}});
      registry.add("ZNCvsZPC", "ZNCvsZPC", {HistType::kTH2F, {{{nBins2, -10., MaxZP}, {nBins2, -10., MaxZN}}}});
    }
    if (doprocessZdcBcAss) { // Check if the process function for ZDCBcAss is enabled
      registry.add("ZNAbc", "ZNAbc", {HistType::kTH1F, {{nBins1, -10., MaxZN}}});
      registry.add("ZPAbc", "ZPAbc", {HistType::kTH1F, {{nBins1, -10., MaxZN}}});
      registry.add("ZNCbc", "ZNCbc", {HistType::kTH1F, {{nBins1, -10., MaxZP}}});
      registry.add("ZPCbc", "ZPCbc", {HistType::kTH1F, {{nBins1, -10., MaxZP}}});
      registry.add("ZEM1bc", "ZEM1bc", {HistType::kTH1F, {{nBins1, -10., MaxZP}}});
      registry.add("ZEM2bc", "ZEM2bc", {HistType::kTH1F, {{nBins1, -10., MaxZP}}});
      registry.add("ZNvsZEMbc", "ZNvsZEMbc", {HistType::kTH2F, {{{nBins2, -10., 10000.5}, {nBins2, -10., 2. * MaxZN}}}});
      registry.add("ZNAvsZNCbc", "ZNAvsZNCbc", {HistType::kTH2F, {{{nBins2, -10., MaxZN}, {nBins2, -10., MaxZN}}}});
      registry.add("ZPAvsZPCbc", "ZPAvsZPCbc", {HistType::kTH2F, {{{nBins1, -10., MaxZP}, {nBins2, -10., MaxZP}}}});
      registry.add("ZNAvsZPAbc", "ZNAvsZPAbc", {HistType::kTH2F, {{{nBins2, -10., MaxZP}, {nBins2, -10., MaxZN}}}});
      registry.add("ZNCvsZPCbc", "ZNCvsZPCbc", {HistType::kTH2F, {{{nBins2, -10., MaxZP}, {nBins2, -10., MaxZN}}}});
    }
    if (doprocessZdcCollAss) { // Check if the process function for ZDCCollAss is enabled
      registry.add("ZNAcoll", "ZNAcoll", {HistType::kTH1F, {{nBins1, -10., MaxZN}}});
      registry.add("ZPAcoll", "ZPAcoll", {HistType::kTH1F, {{nBins1, -10., MaxZN}}});
      registry.add("ZNCcoll", "ZNCcoll", {HistType::kTH1F, {{nBins1, -10., MaxZP}}});
      registry.add("ZPCcoll", "ZPCcoll", {HistType::kTH1F, {{nBins1, -10., MaxZP}}});
      registry.add("ZEM1coll", "ZEM1coll", {HistType::kTH1F, {{nBins1, -10., MaxZP}}});
      registry.add("ZEM2coll", "ZEM2coll", {HistType::kTH1F, {{nBins1, -10., MaxZP}}});
      registry.add("ZNvsZEMcoll", "ZNvsZEMcoll", {HistType::kTH2F, {{{nBins2, -10., 10000.5}, {nBins2, -10., 2. * MaxZN}}}});
      registry.add("ZNAvsZNCcoll", "ZNAvsZNCcoll", {HistType::kTH2F, {{{nBins2, -10., MaxZN}, {nBins2, -10., MaxZN}}}});
      registry.add("ZPAvsZPCcoll", "ZPAvsZPCcoll", {HistType::kTH2F, {{{nBins1, -10., MaxZP}, {nBins2, -10., MaxZP}}}});
      registry.add("ZNAvsZPAcoll", "ZNAvsZPAcoll", {HistType::kTH2F, {{{nBins2, -10., MaxZP}, {nBins2, -10., MaxZN}}}});
      registry.add("ZNCvsZPCcoll", "ZNCvsZPCcoll", {HistType::kTH2F, {{{nBins2, -10., MaxZP}, {nBins2, -10., MaxZN}}}});
    }
    if (doprocessZdcCorrela) { // Check if the process function for ZDCCollCorrela is enabled
      registry.add("ZNvsFV0Acorrel", "ZNvsFV0Acorrel", {HistType::kTH2F, {{{nBins2, 0., MaxMultFV0}, {nBins2, -10., 2. * MaxZN}}}});
      registry.add("ZNvsFT0correl", "ZNvsFT0correl", {HistType::kTH2F, {{{nBins2, 0., MaxMultFT0}, {nBins2, -10., 2. * MaxZN}}}});
      registry.add("ZNvsFDDcorrel", "ZNvsFDDcorrel", {HistType::kTH2F, {{{nBins2, 0., MaxMultFDD}, {nBins2, -10., 2. * MaxZN}}}});
    }
  }

  void processZdcAuto(aod::Zdc const& zdc)
  {
    registry.get<TH1>(HIST("ZNApmc"))->Fill(zdc.energyCommonZNA());
    registry.get<TH1>(HIST("ZNCpmc"))->Fill(zdc.energyCommonZNC());
    registry.get<TH1>(HIST("ZPApmc"))->Fill(zdc.energyCommonZPA());
    registry.get<TH1>(HIST("ZPCpmc"))->Fill(zdc.energyCommonZPC());
    registry.get<TH1>(HIST("ZEM1"))->Fill(zdc.energyZEM1());
    registry.get<TH1>(HIST("ZEM2"))->Fill(zdc.energyZEM2());
    registry.get<TH2>(HIST("ZNvsZEM"))->Fill(zdc.energyZEM1() + zdc.energyZEM2(), zdc.energyCommonZNA() + zdc.energyCommonZNC());
    registry.get<TH2>(HIST("ZNAvsZNC"))->Fill(zdc.energyCommonZNC(), zdc.energyCommonZNA());
    registry.get<TH2>(HIST("ZPAvsZPC"))->Fill(zdc.energyCommonZPC(), zdc.energyCommonZPA());
    registry.get<TH2>(HIST("ZNAvsZPA"))->Fill(zdc.energyCommonZPA(), zdc.energyCommonZNA());
    registry.get<TH2>(HIST("ZNCvsZPC"))->Fill(zdc.energyCommonZPC(), zdc.energyCommonZNC());
  }
  /// name, description, function pointer, default value
  /// note that it has to be declared after the function, so that the pointer is known
  PROCESS_SWITCH(ZDCAnalysis, processZdcAuto, "Processing ZDC autotriggered events", true);

  void processZdcBcAss(
    // soa::Join<aod::BCs, aod::Timestamps> const& bcs,
    BCsRun3 const& bcs,
    // BCsWithBcSels const& bcs,
    // BCsWithRun3Matchings const& bcs,
    aod::Zdcs const& zdcs)
  {
    for (const auto& bc : bcs) {
      if (bc.has_zdc()) {

        registry.get<TH1>(HIST("ZNAbc"))->Fill(bc.zdc().energyCommonZNA());
        registry.get<TH1>(HIST("ZNCbc"))->Fill(bc.zdc().energyCommonZNC());
        registry.get<TH1>(HIST("ZPAbc"))->Fill(bc.zdc().energyCommonZPA());
        registry.get<TH1>(HIST("ZPCbc"))->Fill(bc.zdc().energyCommonZPC());
        registry.get<TH1>(HIST("ZEM1bc"))->Fill(bc.zdc().energyZEM1());
        registry.get<TH1>(HIST("ZEM2bc"))->Fill(bc.zdc().energyZEM2());
        registry.get<TH2>(HIST("ZNvsZEMbc"))->Fill(bc.zdc().energyZEM1() + bc.zdc().energyZEM2(), bc.zdc().energyCommonZNA() + bc.zdc().energyCommonZNC());
        registry.get<TH2>(HIST("ZNAvsZNCbc"))->Fill(bc.zdc().energyCommonZNC(), bc.zdc().energyCommonZNA());
        registry.get<TH2>(HIST("ZPAvsZPCbc"))->Fill(bc.zdc().energyCommonZPC(), bc.zdc().energyCommonZPA());
        registry.get<TH2>(HIST("ZNAvsZPAbc"))->Fill(bc.zdc().energyCommonZPA(), bc.zdc().energyCommonZNA());
        registry.get<TH2>(HIST("ZNCvsZPCbc"))->Fill(bc.zdc().energyCommonZPC(), bc.zdc().energyCommonZNC());
      }
    }
  }
  PROCESS_SWITCH(ZDCAnalysis, processZdcBcAss, "Processing ZDC w. BC association", true);

  void processZdcCollAss(
    soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision,
    aod::Zdcs const& zdcs)
  {
    if (collision.foundZDCId() >= 0) {
      registry.get<TH1>(HIST("ZNAcoll"))->Fill(collision.foundZDC().energyCommonZNA());
      registry.get<TH1>(HIST("ZNCcoll"))->Fill(collision.foundZDC().energyCommonZNC());
      registry.get<TH1>(HIST("ZPAcoll"))->Fill(collision.foundZDC().energyCommonZPA());
      registry.get<TH1>(HIST("ZPCcoll"))->Fill(collision.foundZDC().energyCommonZPC());
      registry.get<TH1>(HIST("ZEM1coll"))->Fill(collision.foundZDC().energyZEM1());
      registry.get<TH1>(HIST("ZEM2coll"))->Fill(collision.foundZDC().energyZEM2());
      registry.get<TH2>(HIST("ZNvsZEMcoll"))->Fill(collision.foundZDC().energyZEM1() + collision.foundZDC().energyZEM2(), collision.foundZDC().energyCommonZNA() + collision.foundZDC().energyCommonZNC());
      registry.get<TH2>(HIST("ZNAvsZNCcoll"))->Fill(collision.foundZDC().energyCommonZNC(), collision.foundZDC().energyCommonZNA());
      registry.get<TH2>(HIST("ZPAvsZPCcoll"))->Fill(collision.foundZDC().energyCommonZPC(), collision.foundZDC().energyCommonZPA());
      registry.get<TH2>(HIST("ZNAvsZPAcoll"))->Fill(collision.foundZDC().energyCommonZPA(), collision.foundZDC().energyCommonZNA());
      registry.get<TH2>(HIST("ZNCvsZPCcoll"))->Fill(collision.foundZDC().energyCommonZPC(), collision.foundZDC().energyCommonZNC());
    }
  }
  PROCESS_SWITCH(ZDCAnalysis, processZdcCollAss, "Processing ZDC w. collision association", true);

  void processZdcCorrela(
    soa::Join<aod::Collisions, aod::EvSels>::iterator const& coll,
    BCsRun3 const& bcs,
    aod::Zdcs const& zdcs,
    aod::FV0As const& fv0as,
    aod::FT0s const& ft0s,
    aod::FDDs const& fdds)
  {
    const auto& foundBC = coll.foundBC_as<BCsRun3>();

    // FT0
    float multT0A = 0;
    float multT0C = 0;
    if (foundBC.has_ft0()) {
      for (auto amplitude : foundBC.ft0().amplitudeA()) {
        multT0A += amplitude;
      }
      for (auto amplitude : foundBC.ft0().amplitudeC()) {
        multT0C += amplitude;
      }
    } else {
      multT0A = multT0C = -999;
    }
    // FV0
    float multV0A = 0;
    if (foundBC.has_fv0a()) {
      for (auto amplitude : foundBC.fv0a().amplitude()) {
        multV0A += amplitude;
      }
    } else {
      multV0A = -999;
    }
    // FDD
    float multFDA = 0;
    float multFDC = 0;
    if (foundBC.has_fdd()) {
      for (auto amplitude : foundBC.fdd().chargeA()) {
        multFDA += amplitude;
      }
      for (auto amplitude : foundBC.fdd().chargeC()) {
        multFDC += amplitude;
      }
    } else {
      multFDA = multFDC = -999;
    }

    if (foundBC.has_zdc()) {
      registry.get<TH2>(HIST("ZNvsFV0Acorrel"))->Fill(multV0A, foundBC.zdc().energyCommonZNA() + foundBC.zdc().energyCommonZNC());
      registry.get<TH2>(HIST("ZNvsFT0correl"))->Fill(multT0A + multT0C, foundBC.zdc().energyCommonZNC() + foundBC.zdc().energyCommonZNA());
      registry.get<TH2>(HIST("ZNvsFDDcorrel"))->Fill(multFDC + multFDA, foundBC.zdc().energyCommonZNC() + foundBC.zdc().energyCommonZNA());
    }
  }
  PROCESS_SWITCH(ZDCAnalysis, processZdcCorrela, "Processing ZDC vs. mult. w. association", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<ZDCAnalysis>(cfgc) //
  };
}