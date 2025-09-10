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

#include "Common/Core/TableHelper.h"
#include "Common/DataModel/PIDResponseTPC.h"

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Array2D.h>
#include <Framework/Configurable.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;

using tinyPidTracks = soa::Join<aod::Tracks, aod::pidTPCEl, aod::pidTPCMu, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr, aod::pidTPCDe, aod::pidTPCTr, aod::pidTPCHe, aod::pidTPCAl>;

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

static constexpr int nParameters = 1;
static const std::vector<std::string> tableNames{"pidTPCFullEl",  // 0
                                                 "pidTPCFullMu",  // 1
                                                 "pidTPCFullPi",  // 2
                                                 "pidTPCFullKa",  // 3
                                                 "pidTPCFullPr",  // 4
                                                 "pidTPCFullDe",  // 5
                                                 "pidTPCFullTr",  // 6
                                                 "pidTPCFullHe",  // 7
                                                 "pidTPCFullAl"}; // 8
static const std::vector<std::string> parameterNames{"enable"};
static const int defaultParameters[nTables][nParameters]{{-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}};

struct Run2TinyToFullPID {
  Produces<aod::pidTPCFullEl> pidTPCFullEl;
  Produces<aod::pidTPCFullMu> pidTPCFullMu;
  Produces<aod::pidTPCFullPi> pidTPCFullPi;
  Produces<aod::pidTPCFullKa> pidTPCFullKa;
  Produces<aod::pidTPCFullPr> pidTPCFullPr;
  Produces<aod::pidTPCFullDe> pidTPCFullDe;
  Produces<aod::pidTPCFullTr> pidTPCFullTr;
  Produces<aod::pidTPCFullHe> pidTPCFullHe;
  Produces<aod::pidTPCFullAl> pidTPCFullAl;

  Configurable<LabeledArray<int>> enabledTables{"enabledTables",
                                                {defaultParameters[0], nTables, nParameters, tableNames, parameterNames},
                                                "Produce full PID tables depending on needs. Autodetect is -1, Force no is 0 and force yes is 1."};

  std::vector<int> mEnabledTables; // Vector of enabled tables

  void init(InitContext& context)
  {
    for (int i = 0; i < nTables; i++) {
      LOGF(info, "test %i", i);
      int f = enabledTables->get(tableNames[i].c_str(), "enable");
      enableFlagIfTableRequired(context, tableNames[i], f);
      if (f == 1) {
        mEnabledTables.push_back(i);
      }
    }
  }

  void process(tinyPidTracks const& tracks)
  {
    // reserve memory
    for (auto i : mEnabledTables) {
      switch (i) {
        case kPidEl:
          pidTPCFullEl.reserve(tracks.size());
          break;
        case kPidMu:
          pidTPCFullMu.reserve(tracks.size());
          break;
        case kPidPi:
          pidTPCFullPi.reserve(tracks.size());
          break;
        case kPidKa:
          pidTPCFullKa.reserve(tracks.size());
          break;
        case kPidPr:
          pidTPCFullPr.reserve(tracks.size());
          break;
        case kPidDe:
          pidTPCFullDe.reserve(tracks.size());
          break;
        case kPidTr:
          pidTPCFullTr.reserve(tracks.size());
          break;
        case kPidHe:
          pidTPCFullHe.reserve(tracks.size());
          break;
        case kPidAl:
          pidTPCFullAl.reserve(tracks.size());
          break;
        default:
          LOG(fatal) << "Unknown table requested: " << i;
          break;
      }
    }

    for (const auto& track : tracks) {
      for (auto i : mEnabledTables) {
        switch (i) {
          case kPidEl:
            pidTPCFullEl(0.0f, track.tpcNSigmaEl());
            break;
          case kPidMu:
            pidTPCFullMu(0.0f, track.tpcNSigmaMu());
            break;
          case kPidPi:
            pidTPCFullPi(0.0f, track.tpcNSigmaPi());
            break;
          case kPidKa:
            pidTPCFullKa(0.0f, track.tpcNSigmaKa());
            break;
          case kPidPr:
            pidTPCFullPr(0.0f, track.tpcNSigmaPr());
            break;
          case kPidDe:
            pidTPCFullDe(0.0f, track.tpcNSigmaDe());
            break;
          case kPidTr:
            pidTPCFullTr(0.0f, track.tpcNSigmaTr());
            break;
          case kPidHe:
            pidTPCFullHe(0.0f, track.tpcNSigmaHe());
            break;
          case kPidAl:
            pidTPCFullAl(0.0f, track.tpcNSigmaAl());
            break;
          default:
            LOG(fatal) << "Unknown table requested: " << i;
            break;
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<Run2TinyToFullPID>(cfgc)};
}
