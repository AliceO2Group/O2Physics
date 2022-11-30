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

void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  std::vector<ConfigParamSpec> options{
    {"zdc-trig", VariantType::Int, 1, {"ZDC triggered events"}},
    {"zdc-bc", VariantType::Int, 1, {"ZDC in BC association"}},
    {"zdc-coll", VariantType::Int, 1, {"ZDC in collision association"}},
    {"zdc-correl", VariantType::Int, 1, {"ZDC correlated with other det.s"}}};
  std::swap(workflowOptions, options);
}

#include "Framework/runDataProcessing.h"

// using BCsWithRun3Matchings = soa::Join<aod::BCs, aod::Timestamps, aod::Run3MatchedToBCSparse>;
// using BCsWithBcSels = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels>;
using BCsRun3 = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;

struct ZDCTriggeredEvents {

  // Configurable for number of bins
  Configurable<int> nBins1{"nBins1", 400, "Nbins 1d ZDC histos"};
  Configurable<int> nBins2{"nBins2", 800, "Nbins 2d ZDC histos"};
  // CInfigurable maximum limit
  Configurable<float> MaxZN{"MaxZN", 8000, "Max ZN signal"};
  Configurable<float> MaxZP{"MaxZP", 6000, "Max ZP signal"};

  HistogramRegistry registry{
    "registry",
    {
      //{"hVertexZ", "hVertexZ", {HistType::kTH1F, {{100, -15., 15.}}}},
      {"ZNApmc", "ZNApmc", {HistType::kTH1F, {{nBins1, -10., MaxZN}}}},
      {"ZPApmc", "ZPApmc", {HistType::kTH1F, {{nBins1, -10., MaxZN}}}},
      {"ZNCpmc", "ZNCpmc", {HistType::kTH1F, {{nBins1, -10., MaxZP}}}},
      {"ZPCpmc", "ZPCpmc", {HistType::kTH1F, {{nBins1, -10., MaxZP}}}},
      {"ZEM1", "ZEM1", {HistType::kTH1F, {{nBins1, -10., MaxZP}}}},
      {"ZEM2", "ZEM2", {HistType::kTH1F, {{nBins1, -10., MaxZP}}}},
      {"ZNvsZEM", "ZNvsZEM", {HistType::kTH2F, {{{nBins2, -10., 10000.5}, {nBins2, -10., 2. * MaxZN}}}}},
      {"ZNAvsZNC", "ZNAvsZNC", {HistType::kTH2F, {{{nBins2, -10., MaxZN}, {nBins2, -10., MaxZN}}}}},
      {"ZPAvsZPC", "ZPAvsZPC", {HistType::kTH2F, {{{nBins1, -10., MaxZP}, {nBins2, -10., MaxZP}}}}},
      {"ZNAvsZPA", "ZNAvsZPA", {HistType::kTH2F, {{{nBins2, -10., MaxZP}, {nBins2, -10., MaxZN}}}}},
      {"ZNCvsZPC", "ZNCvsZPC", {HistType::kTH2F, {{{nBins2, -10., MaxZP}, {nBins2, -10., MaxZN}}}}},
    }};

  void process(aod::Zdc const& zdc)
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
};

struct ZDCBCEvents {

  // Configurable for number of bins
  Configurable<int> nBins1{"nBins1", 400, "Nbins 1d ZDC histos"};
  Configurable<int> nBins2{"nBins2", 800, "Nbins 2d ZDC histos"};
  // CInfigurable maximum limit
  Configurable<float> MaxZN{"MaxZN", 8000, "Max ZN signal"};
  Configurable<float> MaxZP{"MaxZP", 6000, "Max ZP signal"};

  HistogramRegistry registry{
    "registry",
    {
      {"ZNAbc", "ZNAbc", {HistType::kTH1F, {{nBins1, -10., MaxZN}}}},
      {"ZPAbc", "ZPAbc", {HistType::kTH1F, {{nBins1, -10., MaxZN}}}},
      {"ZNCbc", "ZNCbc", {HistType::kTH1F, {{nBins1, -10., MaxZP}}}},
      {"ZPCbc", "ZPCbc", {HistType::kTH1F, {{nBins1, -10., MaxZP}}}},
      {"ZEM1bc", "ZEM1bc", {HistType::kTH1F, {{nBins1, -10., MaxZP}}}},
      {"ZEM2bc", "ZEM2bc", {HistType::kTH1F, {{nBins1, -10., MaxZP}}}},
      {"ZNvsZEMbc", "ZNvsZEMbc", {HistType::kTH2F, {{{nBins2, -10., 10000.5}, {nBins2, -10., 2. * MaxZN}}}}},
      {"ZNAvsZNCbc", "ZNAvsZNCbc", {HistType::kTH2F, {{{nBins2, -10., MaxZN}, {nBins2, -10., MaxZN}}}}},
      {"ZPAvsZPCbc", "ZPAvsZPCbc", {HistType::kTH2F, {{{nBins1, -10., MaxZP}, {nBins2, -10., MaxZP}}}}},
      {"ZNAvsZPAbc", "ZNAvsZPAbc", {HistType::kTH2F, {{{nBins2, -10., MaxZP}, {nBins2, -10., MaxZN}}}}},
      {"ZNCvsZPCbc", "ZNCvsZPCbc", {HistType::kTH2F, {{{nBins2, -10., MaxZP}, {nBins2, -10., MaxZN}}}}},
    }};

  void process(
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
};

struct ZDCCollEvents {

  // Configurable for number of bins
  Configurable<int> nBins1{"nBins1", 400, "Nbins 1d ZDC histos"};
  Configurable<int> nBins2{"nBins2", 800, "Nbins 2d ZDC histos"};
  // CInfigurable maximum limit
  Configurable<float> MaxZN{"MaxZN", 8000, "Max ZN signal"};
  Configurable<float> MaxZP{"MaxZP", 6000, "Max ZP signal"};

  HistogramRegistry registry{
    "registry",
    {
      {"ZNAcoll", "ZNAcoll", {HistType::kTH1F, {{nBins1, -10., MaxZN}}}},
      {"ZPAcoll", "ZPAcoll", {HistType::kTH1F, {{nBins1, -10., MaxZN}}}},
      {"ZNCcoll", "ZNCcoll", {HistType::kTH1F, {{nBins1, -10., MaxZP}}}},
      {"ZPCcoll", "ZPCcoll", {HistType::kTH1F, {{nBins1, -10., MaxZP}}}},
      {"ZEM1coll", "ZEM1coll", {HistType::kTH1F, {{nBins1, -10., MaxZP}}}},
      {"ZEM2coll", "ZEM2coll", {HistType::kTH1F, {{nBins1, -10., MaxZP}}}},
      {"ZNvsZEMcoll", "ZNvsZEMcoll", {HistType::kTH2F, {{{nBins2, -10., 10000.5}, {nBins2, -10., 2. * MaxZN}}}}},
      {"ZNAvsZNCcoll", "ZNAvsZNCcoll", {HistType::kTH2F, {{{nBins2, -10., MaxZN}, {nBins2, -10., MaxZN}}}}},
      {"ZPAvsZPCcoll", "ZPAvsZPCcoll", {HistType::kTH2F, {{{nBins1, -10., MaxZP}, {nBins2, -10., MaxZP}}}}},
      {"ZNAvsZPAcoll", "ZNAvsZPAcoll", {HistType::kTH2F, {{{nBins2, -10., MaxZP}, {nBins2, -10., MaxZN}}}}},
      {"ZNCvsZPCcoll", "ZNCvsZPCcoll", {HistType::kTH2F, {{{nBins2, -10., MaxZP}, {nBins2, -10., MaxZN}}}}},
    }};

  void process(
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
};

struct ZDCCorrelate {

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

  HistogramRegistry registry{
    "registry",
    {
      {"ZNvsFV0Acorrel", "ZNvsFV0Acorrel", {HistType::kTH2F, {{{nBins2, 0., MaxMultFV0}, {nBins2, -10., 2. * MaxZN}}}}},
      {"ZNvsFT0correl", "ZNvsFT0correl", {HistType::kTH2F, {{{nBins2, 0., MaxMultFT0}, {nBins2, -10., 2. * MaxZN}}}}},
      {"ZNvsFDDcorrel", "ZNvsFDDcorrel", {HistType::kTH2F, {{{nBins2, 0., MaxMultFDD}, {nBins2, -10., 2. * MaxZN}}}}},
    }};

  void process(
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
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  const bool nosel = cfgc.options().get<int>("zdc-trig");
  const bool bcsel = cfgc.options().get<int>("zdc-bc");
  const bool colsel = cfgc.options().get<int>("zdc-coll");
  const bool correla = cfgc.options().get<int>("zdc-correl");

  WorkflowSpec workflow{};

  if (nosel) {
    workflow.push_back(adaptAnalysisTask<ZDCTriggeredEvents>(cfgc, TaskName{"zdc-analysis-triggered"}));
  }
  if (bcsel) {
    workflow.push_back(adaptAnalysisTask<ZDCBCEvents>(cfgc, TaskName{"zdc-analysis-bc"}));
  }
  if (colsel) {
    workflow.push_back(adaptAnalysisTask<ZDCCollEvents>(cfgc, TaskName{"zdc-analysis-coll"}));
  }
  if (correla) {
    workflow.push_back(adaptAnalysisTask<ZDCCorrelate>(cfgc, TaskName{"zdc-analysis-corr"}));
  }
  //
  return workflow;
}