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
#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/TriggerAliases.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/Core/TrackSelection.h"
#include "ReconstructionDataFormats/GlobalTrackID.h"
#include "ReconstructionDataFormats/Track.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TObjArray.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::track;
using namespace o2::aod::evsel;

using BCsRun3 = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;
// using BCsWithRun3Matchings = soa::Join<aod::BCs, aod::Timestamps, aod::Run3MatchedToBCSparse>;
using ColEvSels = soa::Join<aod::Collisions, aod::EvSels>;

struct ZDCAnalysis {

  // Configurable number of bins
  Configurable<int> nBinsADC{"nBinsADC", 1000, "nbinsADC"};
  Configurable<int> nBinsAmp{"nBinsAmp", 1025, "nbinsAmp"};
  Configurable<int> nBinsTDC{"nBinsTDC", 480, "nbinsTDC"};
  Configurable<int> nBinsFit{"nBinsFit", 1000, "nbinsFit"};
  // Configurable flags
  Configurable<bool> TDCcut{"TDCcut", false, "Flag for TDC cut"};
  // Configurable limits
  Configurable<float> MaxZN{"MaxZN", 4099.5, "Max ZN signal"};
  Configurable<float> MaxZP{"MaxZP", 3099.5, "Max ZP signal"};
  Configurable<float> MaxZEM{"MaxZEM", 3099.5, "Max ZEM signal"};
  Configurable<float> tdcZNmincut{"tdcZNmincut", -4.0, "Min ZN TDC cut"};
  Configurable<float> tdcZNmaxcut{"tdcZNmaxcut", -4.0, "Max ZN TDC cut"};
  Configurable<float> tdcZPmincut{"tdcZPmincut", -4.0, "Min ZP TDC cut"};
  Configurable<float> tdcZPmaxcut{"tdcZPmaxcut", -4.0, "Max ZP TDC cut"};
  //
  Configurable<float> MaxMultFV0{"MaxMultFV0", 3000, "Max FV0 signal"};
  Configurable<float> MaxMultFT0{"MaxMultFT0", 3000, "Max FT0 signal"};
  Configurable<float> MaxMultFDD{"MaxMultFDD", 80000, "Max FDD signal"};
  //
  HistogramRegistry registry{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&)
  {
    if (doprocessZdcAuto) { // Check if the process function for ZDCAuto is enabled
      registry.add("ZNApmc", "ZNApmc; ZNA amplitude; Entries", {HistType::kTH1F, {{nBinsAmp, -0.5, MaxZN}}});
      registry.add("ZPApmc", "ZPApmc; ZPA amplitude; Entries", {HistType::kTH1F, {{nBinsAmp, -0.5, MaxZP}}});
      registry.add("ZNCpmc", "ZNCpmc; ZNC amplitude; Entries", {HistType::kTH1F, {{nBinsAmp, -0.5, MaxZN}}});
      registry.add("ZPCpmc", "ZPCpmc; ZPC amplitude; Entries", {HistType::kTH1F, {{nBinsAmp, -0.5, MaxZP}}});
      registry.add("ZEM1", "ZEM1; ZEM1 amplitude; Entries", {HistType::kTH1F, {{nBinsAmp, -0.5, MaxZEM}}});
      registry.add("ZEM2", "ZEM2; ZEM2 amplitude; Entries", {HistType::kTH1F, {{nBinsAmp, -0.5, MaxZEM}}});
      registry.add("ZNvsZEM", "ZNvsZEM; ZEM; ZNA+ZNC", {HistType::kTH2F, {{{nBinsAmp, -0.5, MaxZEM}, {nBinsAmp, -0.5, 2. * MaxZN}}}});
      registry.add("ZNAvsZNC", "ZNAvsZNC; ZNC; ZNA", {HistType::kTH2F, {{{nBinsAmp, -0.5, MaxZN}, {nBinsAmp, -0.5, MaxZN}}}});
      registry.add("ZPAvsZPC", "ZPAvsZPC; ZPC; ZPA", {HistType::kTH2F, {{{nBinsAmp, -0.5, MaxZP}, {nBinsAmp, -0.5, MaxZP}}}});
      registry.add("ZNAvsZPA", "ZNAvsZPA; ZPA; ZNA", {HistType::kTH2F, {{{nBinsAmp, -0.5, MaxZP}, {nBinsAmp, -0.5, MaxZN}}}});
      registry.add("ZNCvsZPC", "ZNCvsZPC; ZPC; ZNC", {HistType::kTH2F, {{{nBinsAmp, -0.5, MaxZP}, {nBinsAmp, -0.5, MaxZN}}}});
      //
      registry.add("ZNCcvsZNCsum", "ZNCcvsZNCsum; ZNCC ADC; ZNCsum", {HistType::kTH2F, {{{nBinsADC, -0.5, 3. * MaxZN}, {nBinsADC, -0.5, 3. * MaxZN}}}});
      registry.add("ZNAcvsZNAsum", "ZNAcvsZNAsum ZNAC ADC; ZNAsum", {HistType::kTH2F, {{{nBinsADC, -0.5, 3. * MaxZN}, {nBinsADC, -0.5, 3. * MaxZN}}}});
      registry.add("ZPCcvsZPCsum", "ZPCcvsZPCsum ZPCC ADC; ZPCsum", {HistType::kTH2F, {{{nBinsADC, -0.5, 3. * MaxZP}, {nBinsADC, -0.5, 3. * MaxZP}}}});
      registry.add("ZPAcvsZPAsum", "ZPAcvsZPAsum ZPAC ADC; ZPAsum", {HistType::kTH2F, {{{nBinsADC, -0.5, 3. * MaxZP}, {nBinsADC, -0.5, 3. * MaxZP}}}});
      //
      registry.add("ZNCvstdc", "ZNCvstdc; ZNC amplitude; ZNC TDC", {HistType::kTH2F, {{{480, -13.5, 11.45}, {nBinsAmp, -0.5, MaxZN}}}});
      registry.add("ZNAvstdc", "ZNAvstdc; ZNA amplitude; ZNA TDC", {HistType::kTH2F, {{{480, -13.5, 11.45}, {nBinsAmp, -0.5, MaxZN}}}});
      registry.add("ZPCvstdc", "ZPCvstdc; ZPC amplitude; ZPC TDC", {HistType::kTH2F, {{{480, -13.5, 11.45}, {nBinsAmp, -0.5, MaxZP}}}});
      registry.add("ZPAvstdc", "ZPAvstdc; ZPA amplitude; ZPA TDC", {HistType::kTH2F, {{{480, -13.5, 11.45}, {nBinsAmp, -0.5, MaxZP}}}});
      registry.add("ZEM1vstdc", "ZEM1vstdc; ZEM1 amplitude; ZEM1 TDC", {HistType::kTH2F, {{{480, -13.5, 11.45}, {nBinsAmp, -0.5, MaxZEM}}}});
      registry.add("ZEM2vstdc", "ZEM2vstdc; ZEM2 amplitude; ZEM2 TDC", {HistType::kTH2F, {{{480, -13.5, 11.45}, {nBinsAmp, -0.5, MaxZEM}}}});
    }
    if (doprocessZdcBcAss) { // Check if the process function for ZDCBcAss is enabled
      registry.add("ZNAbc", "ZNAbc; ZNA amplitude; Entries", {HistType::kTH1F, {{nBinsAmp, -0.5, MaxZN}}});
      registry.add("ZNCbc", "ZNCbc; ZNC amplitude; Entries", {HistType::kTH1F, {{nBinsAmp, -0.5, MaxZN}}});
      registry.add("ZPAbc", "ZPAbc; ZPA amplitude; Entries", {HistType::kTH1F, {{nBinsAmp, -0.5, MaxZP}}});
      registry.add("ZPCbc", "ZPCbc; ZPC amplitude; Entries", {HistType::kTH1F, {{nBinsAmp, -0.5, MaxZP}}});
      registry.add("ZEM1bc", "ZEM1bc; ZEM1 amplitude; Entries", {HistType::kTH1F, {{nBinsAmp, -0.5, MaxZEM}}});
      registry.add("ZEM2bc", "ZEM2bc; ZEM2 amplitude; Entries", {HistType::kTH1F, {{nBinsAmp, -0.5, MaxZEM}}});
      registry.add("ZNvsZEMbc", "ZNvsZEMbc; ZEM; ZNA+ZNC", {HistType::kTH2F, {{{nBinsAmp, -0.5, MaxZEM}, {nBinsAmp, -0.5, 2.0 * MaxZN}}}});
      registry.add("ZNAvsZNCbc", "ZNAvsZNCbc; ZNC; ZNA", {HistType::kTH2F, {{{nBinsAmp, -0.5, MaxZN}, {nBinsAmp, -0.5, MaxZN}}}});
      registry.add("ZPAvsZPCbc", "ZPAvsZPCbc; ZPC; ZPA", {HistType::kTH2F, {{{nBinsAmp, -0.5, MaxZP}, {nBinsAmp, -0.5, MaxZP}}}});
      registry.add("ZNAvsZPAbc", "ZNAvsZPAbc; ZPA; ZNA", {HistType::kTH2F, {{{nBinsAmp, -0.5, MaxZP}, {nBinsAmp, -0.5, MaxZN}}}});
      registry.add("ZNCvsZPCbc", "ZNCvsZPCbc; ZPC; ZNC", {HistType::kTH2F, {{{nBinsAmp, -0.5, MaxZP}, {nBinsAmp, -0.5, MaxZN}}}});
    }
    if (doprocessZdcCollAss) { // Check if the process function for ZDCCollAss is enabled
      registry.add("ZNAcoll", "ZNAcoll; ZNA amplitude; Entries", {HistType::kTH1F, {{nBinsAmp, -0.5, MaxZN}}});
      registry.add("ZPAcoll", "ZPAcoll; ZPA amplitude; Entries", {HistType::kTH1F, {{nBinsAmp, -0.5, MaxZP}}});
      registry.add("ZNCcoll", "ZNCcoll; ZNC amplitude; Entries", {HistType::kTH1F, {{nBinsAmp, -0.5, MaxZN}}});
      registry.add("ZPCcoll", "ZPCcoll; ZPC amplitude; Entries", {HistType::kTH1F, {{nBinsAmp, -0.5, MaxZP}}});
      registry.add("ZEM1coll", "ZEM1coll; ZEM1 amplitude; Entries", {HistType::kTH1F, {{nBinsAmp, -0.5, MaxZEM}}});
      registry.add("ZEM2coll", "ZEM2coll; ZEM2 amplitude; Entries", {HistType::kTH1F, {{nBinsAmp, -0.5, MaxZEM}}});
      registry.add("ZNvsZEMcoll", "ZNvsZEMcoll; ZEM; ZNA+ZNC", {HistType::kTH2F, {{{nBinsAmp, -0.5, MaxZEM}, {nBinsAmp, -0.5, 2. * MaxZN}}}});
      registry.add("ZNAvsZNCcoll", "ZNAvsZNCcoll; ZNC; ZNA", {HistType::kTH2F, {{{nBinsAmp, -0.5, MaxZN}, {nBinsAmp, -0.5, MaxZN}}}});
      registry.add("ZPAvsZPCcoll", "ZPAvsZPCcoll; ZPA; ZPC", {HistType::kTH2F, {{{nBinsAmp, -0.5, MaxZP}, {nBinsAmp, -0.5, MaxZP}}}});
      registry.add("ZNAvsZPAcoll", "ZNAvsZPAcoll; ZPA; ZNA", {HistType::kTH2F, {{{nBinsAmp, -0.5, MaxZP}, {nBinsAmp, -0.5, MaxZN}}}});
      registry.add("ZNCvsZPCcoll", "ZNCvsZPCcoll; ZPC; ZNC", {HistType::kTH2F, {{{nBinsAmp, -0.5, MaxZP}, {nBinsAmp, -0.5, MaxZN}}}});
      //
      registry.add("ZNCvstdccoll", "ZNCvstdccoll", {HistType::kTH2F, {{{nBinsTDC, -13.5, 11.45}, {nBinsAmp, -0.5, MaxZN}}}});
      registry.add("ZNAvstdccoll", "ZNAvstdccoll", {HistType::kTH2F, {{{nBinsTDC, -13.5, 11.45}, {nBinsAmp, -0.5, MaxZN}}}});
      registry.add("ZPCvstdccoll", "ZPCvstdccoll", {HistType::kTH2F, {{{nBinsTDC, -13.5, 11.45}, {nBinsAmp, -0.5, MaxZP}}}});
      registry.add("ZPAvstdccoll", "ZPAvstdccoll", {HistType::kTH2F, {{{nBinsTDC, -13.5, 11.45}, {nBinsAmp, -0.5, MaxZP}}}});
      registry.add("ZEM1vstdccoll", "ZEM1vstdccoll", {HistType::kTH2F, {{{nBinsTDC, -13.5, 11.45}, {nBinsAmp, -0.5, MaxZEM}}}});
      registry.add("ZEM2vstdccoll", "ZEM2vstdccoll", {HistType::kTH2F, {{{nBinsTDC, -13.5, 11.45}, {nBinsAmp, -0.5, MaxZEM}}}});
      registry.add("debunch", "ZN sum vs. ZN diff.", {HistType::kTH2F, {{{240, -12., 12.}, {240, -12., 12.}}}});
      registry.add("centroidZNA", "ZNA centroid", {HistType::kTH2F, {{{350, -1.75, 1.75}, {350, -1.75, 1.75}}}});
      registry.add("centroidZNC", "ZNC centroid", {HistType::kTH2F, {{{350, -1.75, 1.75}, {350, -1.75, 1.75}}}});
    }
    if (doprocessZdcCorrela) { // Check if the process function for ZDCCollCorrela is enabled
      registry.add("ZNvsFV0Acorrel", "ZNvsFV0Acorrel", {HistType::kTH2F, {{{nBinsFit, 0., MaxMultFV0}, {nBinsAmp, -0.5, 2. * MaxZN}}}});
      registry.add("ZNvsFT0correl", "ZNvsFT0correl", {HistType::kTH2F, {{{nBinsFit, 0., MaxMultFT0}, {nBinsAmp, -0.5, 2. * MaxZN}}}});
      registry.add("ZNvsFDDcorrel", "ZNvsFDDcorrel", {HistType::kTH2F, {{{nBinsFit, 0., MaxMultFDD}, {nBinsAmp, -0.5, 2. * MaxZN}}}});
    }
  }

  void processZdcAuto(aod::Zdc const& zdc)
  {
    auto aZNA = zdc.amplitudeZNA();
    auto aZNC = zdc.amplitudeZNC();
    auto aZPA = zdc.amplitudeZPA();
    auto aZPC = zdc.amplitudeZPC();
    auto aZEM1 = zdc.amplitudeZEM1();
    auto aZEM2 = zdc.amplitudeZEM2();
    //
    registry.get<TH1>(HIST("ZNApmc"))->Fill(aZNA);
    registry.get<TH1>(HIST("ZNCpmc"))->Fill(aZNC);
    registry.get<TH1>(HIST("ZPApmc"))->Fill(aZPA);
    registry.get<TH1>(HIST("ZPCpmc"))->Fill(aZPC);
    registry.get<TH1>(HIST("ZEM1"))->Fill(aZEM1);
    registry.get<TH1>(HIST("ZEM2"))->Fill(aZEM2);
    registry.get<TH2>(HIST("ZNvsZEM"))->Fill(aZEM1 + aZEM2, aZNA + aZNC);
    registry.get<TH2>(HIST("ZNAvsZNC"))->Fill(aZNC, aZNA);
    registry.get<TH2>(HIST("ZPAvsZPC"))->Fill(aZPC, aZPA);
    registry.get<TH2>(HIST("ZNAvsZPA"))->Fill(aZPA, aZNA);
    registry.get<TH2>(HIST("ZNCvsZPC"))->Fill(aZPC, aZNC);
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
    registry.get<TH2>(HIST("ZNCvstdc"))->Fill(zdc.timeZNC(), aZNC);
    registry.get<TH2>(HIST("ZNAvstdc"))->Fill(zdc.timeZNA(), aZNA);
    registry.get<TH2>(HIST("ZPCvstdc"))->Fill(zdc.timeZPC(), aZPC);
    registry.get<TH2>(HIST("ZPAvstdc"))->Fill(zdc.timeZPA(), aZPA);
    registry.get<TH2>(HIST("ZEM1vstdc"))->Fill(zdc.timeZEM1(), aZEM1);
    registry.get<TH2>(HIST("ZEM2vstdc"))->Fill(zdc.timeZEM2(), aZEM2);
  }
  /// name, description, function pointer, default value
  /// note that it has to be declared after the function, so that the pointer is known
  PROCESS_SWITCH(ZDCAnalysis, processZdcAuto, "Processing ZDC auto-triggered events", true);

  void processZdcBcAss(
    BCsRun3 const& bcs,
    aod::Zdcs const& /*zdcs*/)
  {
    for (const auto& bc : bcs) {
      if (bc.has_zdc()) {

        auto aZNA = bc.zdc().amplitudeZNA();
        auto aZNC = bc.zdc().amplitudeZNC();
        auto aZPA = bc.zdc().amplitudeZPA();
        auto aZPC = bc.zdc().amplitudeZPC();
        auto aZEM1 = bc.zdc().amplitudeZEM1();
        auto aZEM2 = bc.zdc().amplitudeZEM2();

        registry.get<TH1>(HIST("ZNAbc"))->Fill(aZNA);
        registry.get<TH1>(HIST("ZNCbc"))->Fill(aZNC);
        registry.get<TH1>(HIST("ZPAbc"))->Fill(aZPA);
        registry.get<TH1>(HIST("ZPCbc"))->Fill(aZPC);
        registry.get<TH1>(HIST("ZEM1bc"))->Fill(aZEM1);
        registry.get<TH1>(HIST("ZEM2bc"))->Fill(aZEM2);
        registry.get<TH2>(HIST("ZNvsZEMbc"))->Fill(aZEM1 + aZEM2, aZNA + aZNC);
        registry.get<TH2>(HIST("ZNAvsZNCbc"))->Fill(aZNC, aZNA);
        registry.get<TH2>(HIST("ZPAvsZPCbc"))->Fill(aZPC, aZPA);
        registry.get<TH2>(HIST("ZNAvsZPAbc"))->Fill(aZPA, aZNA);
        registry.get<TH2>(HIST("ZNCvsZPCbc"))->Fill(aZPC, aZNC);
      }
    }
  }
  PROCESS_SWITCH(ZDCAnalysis, processZdcBcAss, "Processing ZDC w. BC association", true);

  void processZdcCollAss(
    ColEvSels const& cols,
    BCsRun3 const& /*bcs*/,
    aod::Zdcs const& /*zdcs*/)
  {
    // collision-based event selection
    for (auto& collision : cols) {
      const auto& foundBC = collision.foundBC_as<BCsRun3>();
      if (foundBC.has_zdc()) {
        const auto& zdcread = foundBC.zdc();
        //
        auto aZNA = zdcread.amplitudeZNA();
        auto aZNC = zdcread.amplitudeZNC();
        auto aZPA = zdcread.amplitudeZPA();
        auto aZPC = zdcread.amplitudeZPC();
        auto aZEM1 = zdcread.amplitudeZEM1();
        auto aZEM2 = zdcread.amplitudeZEM2();
        auto tZNA = zdcread.timeZNA();
        auto tZNC = zdcread.timeZNC();
        auto tZPA = zdcread.timeZPA();
        auto tZPC = zdcread.timeZPC();
        //
        if (TDCcut) {
          if ((tZNA >= tdcZNmincut) && (tZNA <= tdcZNmaxcut))
            registry.get<TH1>(HIST("ZNAcoll"))->Fill(aZNA);
          if ((tZNC >= tdcZNmincut) && (tZNC <= tdcZNmaxcut))
            registry.get<TH1>(HIST("ZNCcoll"))->Fill(aZNC);
          if ((tZPA >= tdcZPmincut) && (tZPA <= tdcZPmaxcut))
            registry.get<TH1>(HIST("ZPAcoll"))->Fill(aZPA);
          if ((tZPC >= tdcZPmincut) && (tZPC <= tdcZPmaxcut))
            registry.get<TH1>(HIST("ZPCcoll"))->Fill(aZPC);
          if (((tZNC >= tdcZNmincut) && (tZNC <= tdcZNmaxcut)) && ((tZNA >= tdcZNmincut) && (tZNA <= tdcZNmaxcut)))
            registry.get<TH2>(HIST("ZNvsZEMcoll"))->Fill(aZEM1 + aZEM2, aZNA + aZNC);
          if (((tZNC >= tdcZNmincut) && (tZNC <= tdcZNmaxcut)) && ((tZNA >= tdcZNmincut) && (tZNA <= tdcZNmaxcut)))
            registry.get<TH2>(HIST("ZNAvsZNCcoll"))->Fill(aZNC, aZNA);
          if (((tZPC >= tdcZPmincut) && (tZPC <= tdcZPmaxcut)) && ((tZPA >= tdcZPmincut) && (tZPA <= tdcZPmaxcut)))
            registry.get<TH2>(HIST("ZPAvsZPCcoll"))->Fill(aZPC, aZPA);
          if ((tZNA >= tdcZNmincut) && (tZNA <= tdcZNmaxcut))
            registry.get<TH2>(HIST("ZNAvsZPAcoll"))->Fill(aZPA, aZNA);
          if ((tZNC >= tdcZNmincut) && (tZNC <= tdcZNmaxcut))
            registry.get<TH2>(HIST("ZNCvsZPCcoll"))->Fill(aZPC, aZNC);
          //
        } else {
          registry.get<TH1>(HIST("ZNAcoll"))->Fill(aZNA);
          registry.get<TH1>(HIST("ZNCcoll"))->Fill(aZNC);
          registry.get<TH1>(HIST("ZPAcoll"))->Fill(aZPA);
          registry.get<TH1>(HIST("ZPCcoll"))->Fill(aZPC);
          registry.get<TH2>(HIST("ZNvsZEMcoll"))->Fill(aZEM1 + aZEM2, aZNA + aZNC);
          registry.get<TH2>(HIST("ZNAvsZNCcoll"))->Fill(aZNC, aZNA);
          registry.get<TH2>(HIST("ZPAvsZPCcoll"))->Fill(aZPC, aZPA);
          registry.get<TH2>(HIST("ZNAvsZPAcoll"))->Fill(aZPA, aZNA);
          registry.get<TH2>(HIST("ZNCvsZPCcoll"))->Fill(aZPC, aZNC);
        }
        registry.get<TH1>(HIST("ZEM1coll"))->Fill(aZEM1);
        registry.get<TH1>(HIST("ZEM2coll"))->Fill(aZEM2);
        //
        registry.get<TH2>(HIST("ZNCvstdccoll"))->Fill(tZNC, aZNC);
        registry.get<TH2>(HIST("ZNAvstdccoll"))->Fill(tZNA, aZNA);
        registry.get<TH2>(HIST("ZPCvstdccoll"))->Fill(tZPC, aZPC);
        registry.get<TH2>(HIST("ZPAvstdccoll"))->Fill(tZPA, aZPA);
        registry.get<TH2>(HIST("ZEM1vstdccoll"))->Fill(zdcread.timeZEM1(), aZEM1);
        registry.get<TH2>(HIST("ZEM2vstdccoll"))->Fill(zdcread.timeZEM2(), aZEM2);
        //
        registry.get<TH2>(HIST("debunch"))->Fill(tZNA - tZNC, tZNA + tZNC);
        //
        // Calculating centroid over ZNA and ZNC
        // const float beamEne = 5.36/2.;
        const float x[4] = {-1.75, 1.75, -1.75, 1.75};
        const float y[4] = {-1.75, -1.75, 1.75, 1.75};
        const float alpha = 0.395;
        double numXZNC = 0., numYZNC = 0., denZNC = 0.;
        double numXZNA = 0., numYZNA = 0., denZNA = 0.;
        double centrZNC[2] = {0., 0}, centrZNA[2] = {0., 0.};
        //
        auto zncTower = zdcread.energySectorZNC();
        auto znaTower = zdcread.energySectorZNA();
        auto zncC = zdcread.energyCommonZNC();
        auto znaC = zdcread.energyCommonZNA();
        //
        for (int i = 0; i < 4; i++) {
          if (zncTower[i] > 0.) {
            float wZNC = std::pow(zncTower[i], alpha);
            numXZNC += x[i] * wZNC;
            numYZNC += y[i] * wZNC;
            denZNC += wZNC;
          }
          if (znaTower[i] > 0.) {
            float wZNA = std::pow(znaTower[i], alpha);
            numXZNA += x[i] * wZNA;
            numYZNA += y[i] * wZNA;
            denZNA += wZNA;
          }
        }
        //
        if (denZNC != 0) {
          centrZNC[0] = numXZNC / denZNC;
          centrZNC[1] = numYZNC / denZNC;
        } else {
          centrZNC[0] = centrZNC[1] = 999.;
        }
        //
        if (denZNA != 0) {
          centrZNA[0] = -numXZNA / denZNA;
          centrZNA[1] = numYZNA / denZNA;
        } else {
          centrZNA[0] = centrZNA[1] = 999.;
        }
        //
        registry.get<TH2>(HIST("centroidZNA"))->Fill(centrZNA[0], centrZNA[1], znaC);
        registry.get<TH2>(HIST("centroidZNC"))->Fill(centrZNC[0], centrZNC[1], zncC);
      }
    }
  }
  PROCESS_SWITCH(ZDCAnalysis, processZdcCollAss, "Processing ZDC w. collision association", true);

  void processZdcCorrela(
    soa::Join<aod::Collisions, aod::EvSels>::iterator const& coll,
    BCsRun3 const& /*bcs*/,
    aod::Zdcs const& /*zdcs*/,
    aod::FV0As const& /*fv0as*/,
    aod::FT0s const& /*ft0s*/,
    aod::FDDs const& /*fdds*/)
  {
    const auto& foundBC = coll.foundBC_as<BCsRun3>();

    // FT0
    float multT0A = 0.;
    float multT0C = 0.;
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
      auto const& fdd = foundBC.fdd();
      for (auto const& amplitude : fdd.chargeA()) {
        multFDA += amplitude;
      }
      for (auto const& amplitude : fdd.chargeC()) {
        multFDC += amplitude;
      }
    } else {
      multFDA = multFDC = -999;
    }

    if (foundBC.has_zdc()) {
      const auto& zdcread = foundBC.zdc();
      auto aZNA = zdcread.amplitudeZNA();
      auto aZNC = zdcread.amplitudeZNC();
      registry.get<TH2>(HIST("ZNvsFV0Acorrel"))->Fill(multV0A / 100., aZNA + aZNC);
      registry.get<TH2>(HIST("ZNvsFT0correl"))->Fill((multT0A + multT0C) / 100., aZNC + aZNA);
      registry.get<TH2>(HIST("ZNvsFDDcorrel"))->Fill(multFDC + multFDA, aZNC + aZNA);
    }
  }
  PROCESS_SWITCH(ZDCAnalysis, processZdcCorrela, "Processing ZDC vs. mult. w. collision association", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<ZDCAnalysis>(cfgc) //
  };
}
