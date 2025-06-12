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
//
// This task converts tiny PID tables into full PID tables.
// It is meant to be used with Run 2 converted data to maintain
// full compatibility with any task that may subscribe to the Full
// tables (at the cost of some memory consumption).
// It is also able to produce very simple QA plots on the stored
// quantities (optionally disabled for simplicity)
//
// Warning: expected resolution is NOT provided.

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/DataModel/PIDResponse.h"
#include "TableHelper.h"

using namespace o2;
using namespace o2::framework;

using tinyPidTracks = soa::Join<aod::Tracks, aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTPCFullDe, aod::pidTPCFullTr, aod::pidTPCFullHe, aod::pidTPCFullAl>;

static constexpr int nParameters = 1;
static const std::vector<std::string> tableNames{"Electron", // 0
                                                 "Muon",     // 1
                                                 "Pion",     // 2
                                                 "Kaon",     // 3
                                                 "Proton",   // 4
                                                 "Deuteron", // 5
                                                 "Triton",   // 6
                                                 "Helium",   // 7
                                                 "Alpha"};   // 8
static const std::vector<std::string> parameterNames{"enable"};
static const int defaultParameters[9][nParameters]{{0}, {0}, {1}, {1}, {1}, {0}, {0}, {0}, {0}};

static constexpr int kPidEl = 0;
static constexpr int kPidMu = 1;
static constexpr int kPidPi = 2;
static constexpr int kPidKa = 3;
static constexpr int kPidPr = 4;
static constexpr int kPidDe = 5;
static constexpr int kPidTr = 6;
static constexpr int kPidHe = 7;
static constexpr int kPidAl = 8;
static constexpr int nTables = 9;

struct TpcPidQa {
  Configurable<LabeledArray<int>> enabledTables{"enabledTables",
                                                {defaultParameters[0], nTables, nParameters, tableNames, parameterNames},
                                                "Produce QA for this species: 0 - no, 1 - yes"};
  std::vector<int> mEnabledTables; // Vector of enabled tables

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  ConfigurableAxis axisMomentum{"axisMomentum", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "momentum"};

  ConfigurableAxis axisNSigma{"axisNSigma", {48, -6.0f, 6.0f}, "axisNSigma"};

  void init(InitContext&)
  {
    mEnabledTables.resize(9, 0);

    for (int i = 0; i < nTables; i++) {
      int f = enabledTables->get(tableNames[i].c_str(), "enable");
      if (f == 1) {
        mEnabledTables[i] = 1;
        histos.add(fmt::format("hNSigmaVsPTot{}", tableNames[i]).c_str(), "", kTH2F, {axisMomentum, axisNSigma});
      }
    }
  }

  void process(tinyPidTracks const& tracks)
  {
    for (const auto& track : tracks) {
      if (mEnabledTables[kPidEl]) {
        histos.fill(HIST("hNSigmaVsPTotElectron"), track.p(), track.tpcNSigmaEl());
      }
      if (mEnabledTables[kPidMu]) {
        histos.fill(HIST("hNSigmaVsPTotMuon"), track.p(), track.tpcNSigmaMu());
      }
      if (mEnabledTables[kPidPi]) {
        histos.fill(HIST("hNSigmaVsPTotPion"), track.p(), track.tpcNSigmaPi());
      }
      if (mEnabledTables[kPidKa]) {
        histos.fill(HIST("hNSigmaVsPTotKaon"), track.p(), track.tpcNSigmaKa());
      }
      if (mEnabledTables[kPidPr]) {
        histos.fill(HIST("hNSigmaVsPTotProton"), track.p(), track.tpcNSigmaPr());
      }
      if (mEnabledTables[kPidDe]) {
        histos.fill(HIST("hNSigmaVsPTotDeuteron"), track.p(), track.tpcNSigmaDe());
      }
      if (mEnabledTables[kPidTr]) {
        histos.fill(HIST("hNSigmaVsPTotTriton"), track.p(), track.tpcNSigmaTr());
      }
      if (mEnabledTables[kPidHe]) {
        histos.fill(HIST("hNSigmaVsPTotHelium"), track.p(), track.tpcNSigmaHe());
      }
      if (mEnabledTables[kPidAl]) {
        histos.fill(HIST("hNSigmaVsPTotAlpha"), track.p(), track.tpcNSigmaAl());
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<TpcPidQa>(cfgc)};
}
