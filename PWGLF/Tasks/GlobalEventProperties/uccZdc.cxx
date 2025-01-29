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
/// \file uccZdc.cxx
///
/// \brief task for analysis of UCC with the ZDC
/// \author Omar Vazquez (omar.vazquez.rueda@cern.ch)
/// \since January 29, 2025

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/TriggerAliases.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/GlobalTrackID.h"
#include "ReconstructionDataFormats/Track.h"

using namespace std;
using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::evsel;

using ColEvSels = soa::Join<aod::Collisions, aod::EvSels>;
using BCsRun3 = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels,
                          aod::Run3MatchedToBCSparse>;

struct UccZdc {
  // Configurables, binning
  Configurable<int> nBinsAmp{"nBinsAmp", 1025, "nbinsAmp"};
  Configurable<float> maxZN{"maxZN", 4099.5, "Max ZN signal"};
  Configurable<float> maxZP{"maxZP", 3099.5, "Max ZP signal"};
  Configurable<float> maxZEM{"maxZEM", 3099.5, "Max ZEM signal"};
  Configurable<int> nBinsTDC{"nBinsTDC", 480, "nbinsTDC"};
  Configurable<int> nBinsFit{"nBinsFit", 1000, "nbinsFit"};
  Configurable<float> maxMultFV0{"maxMultFV0", 3000, "Max FV0 signal"};
  Configurable<float> maxMultFT0{"maxMultFT0", 3000, "Max FT0 signal"};

  // Analysis Histograms: Data
  HistogramRegistry histos{
    "histos",
    {},
    OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&)
  {
    // define axes you want to use
    const AxisSpec axisEvent{3, 0., +3.0, ""};
    const AxisSpec axisEta{30, -1.5, +1.5, "#eta"};

    // Histograms
    histos.add("etaHistogram", "etaHistogram", kTH1F, {axisEta});
    histos.add("hEventCounter", "Event counter", kTH1F, {axisEvent});
    histos.add("ZNAcomm", "; ZNA common energy; Entries",
               {HistType::kTH1F, {{nBinsAmp, -0.5, maxZN}}});
    histos.add("ZNCcomm", "; ZNC common energy; Entries",
               {HistType::kTH1F, {{nBinsAmp, -0.5, maxZN}}});
    histos.add("ZNAcoll", "ZNAcoll; ZNA amplitude; Entries",
               {HistType::kTH1F, {{nBinsAmp, -0.5, maxZN}}});
    histos.add("ZPAcoll", "ZPAcoll; ZPA amplitude; Entries",
               {HistType::kTH1F, {{nBinsAmp, -0.5, maxZP}}});
    histos.add("ZNCcoll", "ZNCcoll; ZNC amplitude; Entries",
               {HistType::kTH1F, {{nBinsAmp, -0.5, maxZN}}});
    histos.add("ZPCcoll", "ZPCcoll; ZPC amplitude; Entries",
               {HistType::kTH1F, {{nBinsAmp, -0.5, maxZP}}});
    histos.add("ZEM1coll", "ZEM1coll; ZEM1 amplitude; Entries",
               {HistType::kTH1F, {{nBinsAmp, -0.5, maxZEM}}});
    histos.add("ZEM2coll", "ZEM2coll; ZEM2 amplitude; Entries",
               {HistType::kTH1F, {{nBinsAmp, -0.5, maxZEM}}});
    histos.add("ZNvsZEMcoll", "ZNvsZEMcoll; ZEM; ZNA+ZNC",
               {HistType::kTH2F,
                {{{nBinsAmp, -0.5, maxZEM}, {nBinsAmp, -0.5, 2. * maxZN}}}});
    histos.add("ZNAvsZNCcoll", "ZNAvsZNCcoll; ZNC; ZNA",
               {HistType::kTH2F,
                {{{nBinsAmp, -0.5, maxZN}, {nBinsAmp, -0.5, maxZN}}}});
    histos.add("ZPAvsZPCcoll", "ZPAvsZPCcoll; ZPA; ZPC",
               {HistType::kTH2F,
                {{{nBinsAmp, -0.5, maxZP}, {nBinsAmp, -0.5, maxZP}}}});
    histos.add("ZNAvsZPAcoll", "ZNAvsZPAcoll; ZPA; ZNA",
               {HistType::kTH2F,
                {{{nBinsAmp, -0.5, maxZP}, {nBinsAmp, -0.5, maxZN}}}});
    histos.add("ZNCvsZPCcoll", "ZNCvsZPCcoll; ZPC; ZNC",
               {HistType::kTH2F,
                {{{nBinsAmp, -0.5, maxZP}, {nBinsAmp, -0.5, maxZN}}}});
    //
    histos.add("ZNCvstdccoll", "ZNCvstdccoll; time ZNC; ZNC",
               {HistType::kTH2F,
                {{{nBinsTDC, -13.5, 11.45}, {nBinsAmp, -0.5, maxZN}}}});
    histos.add("ZNAvstdccoll", "ZNAvstdccoll; time ZNA; ZNA",
               {HistType::kTH2F,
                {{{nBinsTDC, -13.5, 11.45}, {nBinsAmp, -0.5, maxZN}}}});
    histos.add("ZPCvstdccoll", "ZPCvstdccoll; time ZPC; ZPC",
               {HistType::kTH2F,
                {{{nBinsTDC, -13.5, 11.45}, {nBinsAmp, -0.5, maxZP}}}});
    histos.add("ZPAvstdccoll", "ZPAvstdccoll; time ZPA; ZPA",
               {HistType::kTH2F,
                {{{nBinsTDC, -13.5, 11.45}, {nBinsAmp, -0.5, maxZP}}}});
    histos.add("ZEM1vstdccoll", "ZEM1vstdccoll; time ZEM1; ZEM1",
               {HistType::kTH2F,
                {{{nBinsTDC, -13.5, 11.45}, {nBinsAmp, -0.5, maxZEM}}}});
    histos.add("ZEM2vstdccoll", "ZEM2vstdccoll; time ZEM2; ZEM2",
               {HistType::kTH2F,
                {{{nBinsTDC, -13.5, 11.45}, {nBinsAmp, -0.5, maxZEM}}}});
    histos.add("debunch", "ZN sum vs. ZN diff.",
               {HistType::kTH2F, {{{240, -12., 12.}, {240, -12., 12.}}}});
    histos.add("ZNvsFV0Acorrel", "ZNvsFV0Acorrel",
               {HistType::kTH2F,
                {{{nBinsFit, 0., maxMultFV0}, {nBinsAmp, -0.5, 2. * maxZN}}}});
    histos.add("ZNvsFT0correl", "ZNvsFT0correl",
               {HistType::kTH2F,
                {{{nBinsFit, 0., maxMultFT0}, {nBinsAmp, -0.5, 2. * maxZN}}}});
  }

  void processZdcCollAss(
    ColEvSels const& cols, BCsRun3 const& /*bcs*/, aod::Zdcs const& /*zdcs*/,
    aod::FV0As const& /*fv0as*/, aod::FT0s const& /*ft0s*/,
    soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA> const& tracks)
  {
    // Collision loop
    for (const auto& collision : cols) {
      histos.fill(HIST("hEventCounter"), 0.5);
      if (!collision.sel8()) {
        return;
      }

      histos.fill(HIST("hEventCounter"), 1.5);

      const auto& foundBC = collision.foundBC_as<BCsRun3>();
      if (foundBC.has_zdc()) {
        histos.fill(HIST("hEventCounter"), 2.5);

        const auto& zdcread = foundBC.zdc();
        auto aZNA = zdcread.amplitudeZNA();
        auto aZNC = zdcread.amplitudeZNC();
        auto aZPA = zdcread.amplitudeZPA();
        auto aZPC = zdcread.amplitudeZPC();
        auto aZEM1 = zdcread.amplitudeZEM1();
        auto aZEM2 = zdcread.amplitudeZEM2();
        auto tZEM1 = zdcread.timeZEM1();
        auto tZEM2 = zdcread.timeZEM2();
        auto tZNA = zdcread.timeZNA();
        auto tZNC = zdcread.timeZNC();
        auto tZPA = zdcread.timeZPA();
        auto tZPC = zdcread.timeZPC();
        float multT0A{0.};
        float multT0C{0.};
        float multV0A{0.};

        if (foundBC.has_ft0()) {
          for (const auto& amplitude : foundBC.ft0().amplitudeA()) {
            multT0A += amplitude;
          }
          for (const auto& amplitude : foundBC.ft0().amplitudeC()) {
            multT0C += amplitude;
          }
        } else {
          multT0A = multT0C = -999;
        }

        if (foundBC.has_fv0a()) {
          for (const auto& amplitude : foundBC.fv0a().amplitude()) {
            multV0A += amplitude;
          }
        } else {
          multV0A = -999;
        }

        histos.get<TH1>(HIST("ZNAcomm"))->Fill(zdcread.energyCommonZNA());
        histos.get<TH1>(HIST("ZNCcomm"))->Fill(zdcread.energyCommonZNC());
        histos.get<TH1>(HIST("ZNAcoll"))->Fill(aZNA);
        histos.get<TH1>(HIST("ZNCcoll"))->Fill(aZNC);
        histos.get<TH1>(HIST("ZPAcoll"))->Fill(aZPA);
        histos.get<TH1>(HIST("ZPCcoll"))->Fill(aZPC);
        histos.get<TH2>(HIST("ZNvsZEMcoll"))->Fill(aZEM1 + aZEM2, aZNA + aZNC);
        histos.get<TH2>(HIST("ZNAvsZNCcoll"))->Fill(aZNC, aZNA);
        histos.get<TH2>(HIST("ZPAvsZPCcoll"))->Fill(aZPC, aZPA);
        histos.get<TH2>(HIST("ZNAvsZPAcoll"))->Fill(aZPA, aZNA);
        histos.get<TH2>(HIST("ZNCvsZPCcoll"))->Fill(aZPC, aZNC);

        histos.get<TH1>(HIST("ZEM1coll"))->Fill(aZEM1);
        histos.get<TH1>(HIST("ZEM2coll"))->Fill(aZEM2);
        histos.get<TH2>(HIST("ZNCvstdccoll"))->Fill(tZNC, aZNC);
        histos.get<TH2>(HIST("ZNAvstdccoll"))->Fill(tZNA, aZNA);
        histos.get<TH2>(HIST("ZPCvstdccoll"))->Fill(tZPC, aZPC);
        histos.get<TH2>(HIST("ZPAvstdccoll"))->Fill(tZPA, aZPA);
        histos.get<TH2>(HIST("ZEM1vstdccoll"))->Fill(tZEM1, aZEM1);
        histos.get<TH2>(HIST("ZEM2vstdccoll"))->Fill(tZEM2, aZEM2);
        histos.get<TH2>(HIST("debunch"))->Fill(tZNA - tZNC, tZNA + tZNC);

        histos.get<TH2>(HIST("ZNvsFV0Acorrel"))
          ->Fill(multV0A / 100., aZNA + aZNC);
        histos.get<TH2>(HIST("ZNvsFT0correl"))
          ->Fill((multT0A + multT0C) / 100., aZNC + aZNA);
      } // foundBC.has_zdc()

      for (const auto& track : tracks) {
        if (track.tpcNClsCrossedRows() < 70)
          continue;
        if (std::fabs(track.dcaXY()) > 0.2)
          continue;

        histos.fill(HIST("etaHistogram"), track.eta());
      }
    }
  }
  PROCESS_SWITCH(UccZdc, processZdcCollAss,
                 "Processing ZDC w. collision association", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<UccZdc>(cfgc)};
}
