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

using BCsWithBcSels = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels>;
// using BCsRun3 = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;
using BCsWithRun3Matchings = soa::Join<aod::BCs, aod::Timestamps, aod::Run3MatchedToBCSparse>;

struct ZDCAnalysis {

  // Configurable for number of bins
  Configurable<int> nBins1{"nBins1", 400, "Nbins1"};
  Configurable<int> nBins2{"nBins2", 800, "Nbins2"};
  // CInfigurable maximum limit
  Configurable<float> MaxZN{"MaxZN", 4000, "Max ZN signal"};
  Configurable<float> MaxZP{"MaxZP", 3000, "Max ZP signal"};
  Configurable<float> MaxMultFV0{"MaxMultFV0", 3000, "Max FV0 signal"};
  Configurable<float> MaxMultFT0{"MaxMultFT0", 3000, "Max FT0 signal"};
  Configurable<float> MaxMultFDD{"MaxMultFDD", 80000, "Max FDD signal"};
  Configurable<float> MaxMultNTracks{"MaxMultNTracks", 1000, "Max Ntracks"};

  HistogramRegistry registry;

  void init(InitContext const&)
  {
    if (doprocessZdcAuto) { // Check if the process function for ZDCAuto is enabled
      registry.add("ZNApmc", "ZNApmc", {HistType::kTH1F, {{nBins1, -10., MaxZN}}});
      registry.add("ZPApmc", "ZPApmc", {HistType::kTH1F, {{nBins1, -10., MaxZP}}});
      registry.add("ZNCpmc", "ZNCpmc", {HistType::kTH1F, {{nBins1, -10., MaxZN}}});
      registry.add("ZPCpmc", "ZPCpmc", {HistType::kTH1F, {{nBins1, -10., MaxZP}}});
      registry.add("ZEM1", "ZEM1", {HistType::kTH1F, {{nBins1, -10., MaxZP}}});
      registry.add("ZEM2", "ZEM2", {HistType::kTH1F, {{nBins1, -10., MaxZP}}});
      registry.add("ZNvsZEM", "ZNvsZEM", {HistType::kTH2F, {{{nBins1, -10., MaxZP}, {nBins1, -10., MaxZN}}}});
      registry.add("ZNAvsZNC", "ZNAvsZNC", {HistType::kTH2F, {{{nBins1, -10., MaxZN}, {nBins1, -10., MaxZN}}}});
      registry.add("ZPAvsZPC", "ZPAvsZPC", {HistType::kTH2F, {{{nBins1, -10., MaxZP}, {nBins1, -10., MaxZP}}}});
      registry.add("ZNAvsZPA", "ZNAvsZPA", {HistType::kTH2F, {{{nBins1, -10., MaxZP}, {nBins1, -10., MaxZN}}}});
      registry.add("ZNCvsZPC", "ZNCvsZPC", {HistType::kTH2F, {{{nBins1, -10., MaxZP}, {nBins1, -10., MaxZN}}}});
      //
      registry.add("ZNCcvsZNCsum", "ZNCcvsZNCsum", {HistType::kTH2F, {{{nBins1, -10., 3. * MaxZN}, {nBins1, -10., 3. * MaxZN}}}});
      registry.add("ZNAcvsZNAsum", "ZNAcvsZNAsum", {HistType::kTH2F, {{{nBins1, -10., 3. * MaxZN}, {nBins1, -10., 3. * MaxZN}}}});
      registry.add("ZPCcvsZPCsum", "ZPCcvsZPCsum", {HistType::kTH2F, {{{nBins1, -10., 3. * MaxZP}, {nBins1, -10., 3. * MaxZP}}}});
      registry.add("ZPAcvsZPAsum", "ZPAcvsZPAsum", {HistType::kTH2F, {{{nBins1, -10., 3. * MaxZP}, {nBins1, -10., 3. * MaxZP}}}});
      //
      registry.add("ZNCadcvstdc", "ZNCadcvstdc", {HistType::kTH2F, {{{400, -50., 50.}, {nBins1, -10., MaxZN}}}});
      registry.add("ZNAadcvstdc", "ZNAadcvstdc", {HistType::kTH2F, {{{400, -50., 50.}, {nBins1, -10., MaxZN}}}});
      registry.add("ZPCadcvstdc", "ZPCadcvstdc", {HistType::kTH2F, {{{400, -50., 50.}, {nBins1, -10., MaxZP}}}});
      registry.add("ZPAadcvstdc", "ZPAadcvstdc", {HistType::kTH2F, {{{400, -50., 50.}, {nBins1, -10., MaxZP}}}});
      registry.add("ZEM1adcvstdc", "ZEM1adcvstdc", {HistType::kTH2F, {{{400, -50., 50.}, {nBins1, -10., MaxZP}}}});
      registry.add("ZEM2adcvstdc", "ZEM2adcvstdc", {HistType::kTH2F, {{{400, -50., 50.}, {nBins1, -10., MaxZP}}}});
    }
    if (doprocessZdcBcAss) { // Check if the process function for ZDCBcAss is enabled
      registry.add("ZNAbc", "ZNAbc", {HistType::kTH1F, {{nBins1, -10., MaxZN}}});
      registry.add("ZNCbc", "ZNCbc", {HistType::kTH1F, {{nBins1, -10., MaxZN}}});
      registry.add("ZPAbc", "ZPAbc", {HistType::kTH1F, {{nBins1, -10., MaxZP}}});
      registry.add("ZPCbc", "ZPCbc", {HistType::kTH1F, {{nBins1, -10., MaxZP}}});
      registry.add("ZEM1bc", "ZEM1bc", {HistType::kTH1F, {{nBins1, -10., MaxZP}}});
      registry.add("ZEM2bc", "ZEM2bc", {HistType::kTH1F, {{nBins1, -10., MaxZP}}});
      registry.add("ZNvsZEMbc", "ZNvsZEMbc", {HistType::kTH2F, {{{nBins1, -10., MaxZP}, {nBins1, -10., MaxZN}}}});
      registry.add("ZNAvsZNCbc", "ZNAvsZNCbc", {HistType::kTH2F, {{{nBins1, -10., MaxZN}, {nBins1, -10., MaxZN}}}});
      registry.add("ZPAvsZPCbc", "ZPAvsZPCbc", {HistType::kTH2F, {{{nBins1, -10., MaxZP}, {nBins1, -10., MaxZP}}}});
      registry.add("ZNAvsZPAbc", "ZNAvsZPAbc", {HistType::kTH2F, {{{nBins1, -10., MaxZP}, {nBins1, -10., MaxZN}}}});
      registry.add("ZNCvsZPCbc", "ZNCvsZPCbc", {HistType::kTH2F, {{{nBins1, -10., MaxZP}, {nBins1, -10., MaxZN}}}});
    }
    if (doprocessZdcCollAss) { // Check if the process function for ZDCCollAss is enabled
      registry.add("ZNAcoll", "ZNAcoll", {HistType::kTH1F, {{nBins1, -10., MaxZN}}});
      registry.add("ZPAcoll", "ZPAcoll", {HistType::kTH1F, {{nBins1, -10., MaxZP}}});
      registry.add("ZNCcoll", "ZNCcoll", {HistType::kTH1F, {{nBins1, -10., MaxZN}}});
      registry.add("ZPCcoll", "ZPCcoll", {HistType::kTH1F, {{nBins1, -10., MaxZP}}});
      registry.add("ZEM1coll", "ZEM1coll", {HistType::kTH1F, {{nBins1, -10., MaxZP}}});
      registry.add("ZEM2coll", "ZEM2coll", {HistType::kTH1F, {{nBins1, -10., MaxZP}}});
      registry.add("ZNvsZEMcoll", "ZNvsZEMcoll", {HistType::kTH2F, {{{nBins1, -10., MaxZP}, {nBins1, -10., MaxZN}}}});
      registry.add("ZNAvsZNCcoll", "ZNAvsZNCcoll", {HistType::kTH2F, {{{nBins1, -10., MaxZN}, {nBins1, -10., MaxZN}}}});
      registry.add("ZPAvsZPCcoll", "ZPAvsZPCcoll", {HistType::kTH2F, {{{nBins1, -10., MaxZP}, {nBins1, -10., MaxZP}}}});
      registry.add("ZNAvsZPAcoll", "ZNAvsZPAcoll", {HistType::kTH2F, {{{nBins1, -10., MaxZP}, {nBins1, -10., MaxZN}}}});
      registry.add("ZNCvsZPCcoll", "ZNCvsZPCcoll", {HistType::kTH2F, {{{nBins1, -10., MaxZP}, {nBins1, -10., MaxZN}}}});
      //
      registry.add("ZNCadcvstdccoll", "ZNCadcvstdccoll", {HistType::kTH2F, {{{400, -50., 50.}, {nBins1, -10., MaxZN}}}});
      registry.add("ZNAadcvstdccoll", "ZNAadcvstdccoll", {HistType::kTH2F, {{{400, -50., 50.}, {nBins1, -10., MaxZN}}}});
      registry.add("ZPCadcvstdccoll", "ZPCadcvstdccoll", {HistType::kTH2F, {{{400, -50., 50.}, {nBins1, -10., MaxZP}}}});
      registry.add("ZPAadcvstdccoll", "ZPAadcvstdccoll", {HistType::kTH2F, {{{400, -50., 50.}, {nBins1, -10., MaxZP}}}});
      registry.add("ZEM1adcvstdccoll", "ZEM1adcvstdccoll", {HistType::kTH2F, {{{400, -50., 50.}, {nBins1, -10., MaxZP}}}});
      registry.add("ZEM2adcvstdccoll", "ZEM2adcvstdccoll", {HistType::kTH2F, {{{400, -50., 50.}, {nBins1, -10., MaxZP}}}});
    }
    if (doprocessZdcCorrela) { // Check if the process function for ZDCCollCorrela is enabled
      registry.add("ZNvsFV0Acorrel", "ZNvsFV0Acorrel", {HistType::kTH2F, {{{nBins2, 0., MaxMultFV0}, {nBins2, -10., MaxZN}}}});
      registry.add("ZNvsFT0correl", "ZNvsFT0correl", {HistType::kTH2F, {{{nBins2, 0., MaxMultFT0}, {nBins2, -10., MaxZN}}}});
      registry.add("ZNvsFDDcorrel", "ZNvsFDDcorrel", {HistType::kTH2F, {{{nBins2, 0., MaxMultFDD}, {nBins2, -10., MaxZN}}}});
    }
  }

  void processZdcAuto(aod::Zdc const& zdc)
  {
    registry.get<TH1>(HIST("ZNApmc"))->Fill(zdc.amplitudeZNA());
    registry.get<TH1>(HIST("ZNCpmc"))->Fill(zdc.amplitudeZNC());
    registry.get<TH1>(HIST("ZPApmc"))->Fill(zdc.amplitudeZPA());
    registry.get<TH1>(HIST("ZPCpmc"))->Fill(zdc.amplitudeZPC());
    registry.get<TH1>(HIST("ZEM1"))->Fill(zdc.amplitudeZEM1());
    registry.get<TH1>(HIST("ZEM2"))->Fill(zdc.amplitudeZEM2());
    registry.get<TH2>(HIST("ZNvsZEM"))->Fill(zdc.amplitudeZEM1() + zdc.amplitudeZEM2(), zdc.amplitudeZNA() + zdc.amplitudeZNC());
    registry.get<TH2>(HIST("ZNAvsZNC"))->Fill(zdc.amplitudeZNC(), zdc.amplitudeZNA());
    registry.get<TH2>(HIST("ZPAvsZPC"))->Fill(zdc.amplitudeZPC(), zdc.amplitudeZPA());
    registry.get<TH2>(HIST("ZNAvsZPA"))->Fill(zdc.amplitudeZPA(), zdc.amplitudeZNA());
    registry.get<TH2>(HIST("ZNCvsZPC"))->Fill(zdc.amplitudeZPC(), zdc.amplitudeZNC());
    //
    float sumZNC = (zdc.energySectorZNC())[0] + (zdc.energySectorZNC())[1] + (zdc.energySectorZNC())[2] + (zdc.energySectorZNC())[3];
    float sumZNA = (zdc.energySectorZNA())[0] + (zdc.energySectorZNA())[1] + (zdc.energySectorZNA())[2] + (zdc.energySectorZNA())[3];
    float sumZPC = (zdc.energySectorZPC())[0] + (zdc.energySectorZPC())[1] + (zdc.energySectorZPC())[2] + (zdc.energySectorZPC())[3];
    float sumZPA = (zdc.energySectorZPA())[0] + (zdc.energySectorZPA())[1] + (zdc.energySectorZPA())[2] + (zdc.energySectorZPA())[3];
    registry.get<TH2>(HIST("ZNCcvsZNCsum"))->Fill(sumZNC, zdc.energyCommonZNC());
    registry.get<TH2>(HIST("ZNAcvsZNAsum"))->Fill(sumZNA, zdc.energyCommonZNA());
    registry.get<TH2>(HIST("ZPCcvsZPCsum"))->Fill(sumZPC, zdc.energyCommonZPC());
    registry.get<TH2>(HIST("ZPAcvsZPAsum"))->Fill(sumZPA, zdc.energyCommonZPA());
    //
    registry.get<TH2>(HIST("ZNCadcvstdc"))->Fill(zdc.timeZNC(), zdc.amplitudeZNC());
    registry.get<TH2>(HIST("ZNAadcvstdc"))->Fill(zdc.timeZNA(), zdc.amplitudeZNA());
    registry.get<TH2>(HIST("ZPCadcvstdc"))->Fill(zdc.timeZPC(), zdc.amplitudeZPC());
    registry.get<TH2>(HIST("ZPAadcvstdc"))->Fill(zdc.timeZPA(), zdc.amplitudeZPA());
    registry.get<TH2>(HIST("ZEM1adcvstdc"))->Fill(zdc.timeZEM1(), zdc.amplitudeZEM1());
    registry.get<TH2>(HIST("ZEM2adcvstdc"))->Fill(zdc.timeZEM2(), zdc.amplitudeZEM2());
  }
  /// name, description, function pointer, default value
  /// note that it has to be declared after the function, so that the pointer is known
  PROCESS_SWITCH(ZDCAnalysis, processZdcAuto, "Processing ZDC autotriggered events", true);

  void processZdcBcAss(
    // soa::Join<aod::BCs, aod::Timestamps> const& bcs,
    // BCsRun3 const& bcs,
    // BCsWithBcSels const& bcs,
    BCsWithRun3Matchings const& bcs,
    aod::Zdcs const& zdcs)
  {
    for (const auto& bc : bcs) {
      if (bc.has_zdc()) {

        registry.get<TH1>(HIST("ZNAbc"))->Fill(bc.zdc().amplitudeZNA());
        registry.get<TH1>(HIST("ZNCbc"))->Fill(bc.zdc().amplitudeZNC());
        registry.get<TH1>(HIST("ZPAbc"))->Fill(bc.zdc().amplitudeZPA());
        registry.get<TH1>(HIST("ZPCbc"))->Fill(bc.zdc().amplitudeZPC());
        registry.get<TH1>(HIST("ZEM1bc"))->Fill(bc.zdc().amplitudeZEM1());
        registry.get<TH1>(HIST("ZEM2bc"))->Fill(bc.zdc().amplitudeZEM2());
        registry.get<TH2>(HIST("ZNvsZEMbc"))->Fill(bc.zdc().amplitudeZEM1() + bc.zdc().amplitudeZEM2(), bc.zdc().amplitudeZNA() + bc.zdc().amplitudeZNC());
        registry.get<TH2>(HIST("ZNAvsZNCbc"))->Fill(bc.zdc().amplitudeZNC(), bc.zdc().amplitudeZNA());
        registry.get<TH2>(HIST("ZPAvsZPCbc"))->Fill(bc.zdc().amplitudeZPC(), bc.zdc().amplitudeZPA());
        registry.get<TH2>(HIST("ZNAvsZPAbc"))->Fill(bc.zdc().amplitudeZPA(), bc.zdc().amplitudeZNA());
        registry.get<TH2>(HIST("ZNCvsZPCbc"))->Fill(bc.zdc().amplitudeZPC(), bc.zdc().amplitudeZNC());
      }
    }
  }
  PROCESS_SWITCH(ZDCAnalysis, processZdcBcAss, "Processing ZDC w. BC association", true);

  void processZdcCollAss(
    soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision,
    aod::Zdcs const& zdcs)
  {
    if (collision.foundZDCId() >= 0) {
      registry.get<TH1>(HIST("ZNAcoll"))->Fill(collision.foundZDC().amplitudeZNA());
      registry.get<TH1>(HIST("ZNCcoll"))->Fill(collision.foundZDC().amplitudeZNC());
      registry.get<TH1>(HIST("ZPAcoll"))->Fill(collision.foundZDC().amplitudeZPA());
      registry.get<TH1>(HIST("ZPCcoll"))->Fill(collision.foundZDC().amplitudeZPC());
      registry.get<TH1>(HIST("ZEM1coll"))->Fill(collision.foundZDC().amplitudeZEM1());
      registry.get<TH1>(HIST("ZEM2coll"))->Fill(collision.foundZDC().amplitudeZEM2());
      registry.get<TH2>(HIST("ZNvsZEMcoll"))->Fill(collision.foundZDC().amplitudeZEM1() + collision.foundZDC().amplitudeZEM2(), collision.foundZDC().amplitudeZNA() + collision.foundZDC().amplitudeZNC());
      registry.get<TH2>(HIST("ZNAvsZNCcoll"))->Fill(collision.foundZDC().amplitudeZNC(), collision.foundZDC().amplitudeZNA());
      registry.get<TH2>(HIST("ZPAvsZPCcoll"))->Fill(collision.foundZDC().amplitudeZPC(), collision.foundZDC().amplitudeZPA());
      registry.get<TH2>(HIST("ZNAvsZPAcoll"))->Fill(collision.foundZDC().amplitudeZPA(), collision.foundZDC().amplitudeZNA());
      registry.get<TH2>(HIST("ZNCvsZPCcoll"))->Fill(collision.foundZDC().amplitudeZPC(), collision.foundZDC().amplitudeZNC());
      //
      registry.get<TH2>(HIST("ZNCadcvstdccoll"))->Fill(collision.foundZDC().timeZNC(), collision.foundZDC().amplitudeZNC());
      registry.get<TH2>(HIST("ZNAadcvstdccoll"))->Fill(collision.foundZDC().timeZNA(), collision.foundZDC().amplitudeZNA());
      registry.get<TH2>(HIST("ZPCadcvstdccoll"))->Fill(collision.foundZDC().timeZPC(), collision.foundZDC().amplitudeZPC());
      registry.get<TH2>(HIST("ZPAadcvstdccoll"))->Fill(collision.foundZDC().timeZPA(), collision.foundZDC().amplitudeZPA());
      registry.get<TH2>(HIST("ZEM1adcvstdccoll"))->Fill(collision.foundZDC().timeZEM1(), collision.foundZDC().amplitudeZEM1());
      registry.get<TH2>(HIST("ZEM2adcvstdccoll"))->Fill(collision.foundZDC().timeZEM2(), collision.foundZDC().amplitudeZEM2());
    }
  }
  PROCESS_SWITCH(ZDCAnalysis, processZdcCollAss, "Processing ZDC w. collision association", true);

  void processZdcCorrela(
    soa::Join<aod::Collisions, aod::EvSels>::iterator const& coll,
    BCsWithRun3Matchings const& bcs,
    aod::Zdcs const& zdcs,
    aod::FV0As const& fv0as,
    aod::FT0s const& ft0s,
    aod::FDDs const& fdds)
  {
    const auto& foundBC = coll.foundBC_as<BCsWithRun3Matchings>();

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
      registry.get<TH2>(HIST("ZNvsFV0Acorrel"))->Fill(multV0A / 100., foundBC.zdc().amplitudeZNA() + foundBC.zdc().amplitudeZNC());
      registry.get<TH2>(HIST("ZNvsFT0correl"))->Fill((multT0A + multT0C) / 100., foundBC.zdc().amplitudeZNC() + foundBC.zdc().amplitudeZNA());
      registry.get<TH2>(HIST("ZNvsFDDcorrel"))->Fill(multFDC + multFDA, foundBC.zdc().amplitudeZNC() + foundBC.zdc().amplitudeZNA());
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
