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
/// \file singleTrackSelectorPIDMaker.cxx
/// \brief creates dummy tables for PID columns that are not in the derived data
/// \author Sofia Tomassini, Gleb Romanenko, Nicolò Jacazio
/// \since 22 January 2025

#include <fairlogger/Logger.h>
#include <Framework/AnalysisDataModel.h>

#include <vector>
#include <string>
#include <unordered_map>

#include "PWGCF/Femto3D/DataModel/singletrackselector.h"

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::track;
using namespace o2::aod;
//::singletrackselector; // the namespace defined in .h

struct StPidEl {
  Produces<o2::aod::SinglePIDEls> table;
  void process(o2::aod::SingleTrackSels const& tracks)
  {
    table.reserve(tracks.size());
    for (int i = 0; i < tracks.size(); i++) {
      table(singletrackselector::binning::nsigma::underflowBin,
            singletrackselector::binning::nsigma::underflowBin);
    }
  }
};
struct StPidPi {
  Produces<o2::aod::SinglePIDPis> table;
  void process(o2::aod::SingleTrackSels const& tracks)
  {
    table.reserve(tracks.size());
    for (int i = 0; i < tracks.size(); i++) {
      table(singletrackselector::binning::nsigma::underflowBin,
            singletrackselector::binning::nsigma::underflowBin);
    }
  }
};
struct StPidKa {
  Produces<o2::aod::SinglePIDKas> table;
  void process(o2::aod::SingleTrackSels const& tracks)
  {
    table.reserve(tracks.size());
    for (int i = 0; i < tracks.size(); i++) {
      table(singletrackselector::binning::nsigma::underflowBin,
            singletrackselector::binning::nsigma::underflowBin);
    }
  }
};
struct StPidPr {
  Produces<o2::aod::SinglePIDPrs> table;
  void process(o2::aod::SingleTrackSels const& tracks)
  {
    table.reserve(tracks.size());
    for (int i = 0; i < tracks.size(); i++) {
      table(singletrackselector::binning::nsigma::underflowBin,
            singletrackselector::binning::nsigma::underflowBin);
    }
  }
};
struct StPidDe {
  Produces<o2::aod::SinglePIDDes> table;
  void process(o2::aod::SingleTrackSels const& tracks)
  {
    table.reserve(tracks.size());
    for (int i = 0; i < tracks.size(); i++) {
      table(singletrackselector::binning::nsigma::underflowBin,
            singletrackselector::binning::nsigma::underflowBin);
    }
  }
};
struct StPidTr {
  Produces<o2::aod::SinglePIDTrs> table;
  void process(o2::aod::SingleTrackSels const& tracks)
  {
    table.reserve(tracks.size());
    for (int i = 0; i < tracks.size(); i++) {
      table(singletrackselector::binning::nsigma::underflowBin,
            singletrackselector::binning::nsigma::underflowBin);
    }
  }
};
struct StPidHe {
  Produces<o2::aod::SinglePIDHes> table;
  void process(o2::aod::SingleTrackSels const& tracks)
  {
    table.reserve(tracks.size());
    for (int i = 0; i < tracks.size(); i++) {
      table(singletrackselector::binning::nsigma::underflowBin,
            singletrackselector::binning::nsigma::underflowBin);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  auto workflow = WorkflowSpec{};

  // Check if 'aod-metadata-tables' option is available in the config context
  if (cfgc.options().hasOption("aod-metadata-tables")) {
    const std::vector<std::string> tables = cfgc.options().get<std::vector<std::string>>("aod-metadata-tables");

    // Map of table names to their corresponding converter task functions
    std::unordered_map<std::string, std::vector<std::function<void()>>> tableToTasks = {
      {"O2singlepidel", {[&]() { workflow.push_back(adaptAnalysisTask<StPidEl>(cfgc)); }}},
      {"O2singlepidpi", {[&]() { workflow.push_back(adaptAnalysisTask<StPidPi>(cfgc)); }}},
      {"O2singlepidka", {[&]() { workflow.push_back(adaptAnalysisTask<StPidKa>(cfgc)); }}},
      {"O2singlepidpr", {[&]() { workflow.push_back(adaptAnalysisTask<StPidPr>(cfgc)); }}},
      {"O2singlepidde", {[&]() { workflow.push_back(adaptAnalysisTask<StPidDe>(cfgc)); }}},
      {"O2singlepidtr", {[&]() { workflow.push_back(adaptAnalysisTask<StPidTr>(cfgc)); }}},
      {"O2singlepidhe", {[&]() { workflow.push_back(adaptAnalysisTask<StPidHe>(cfgc)); }}}

    };

    for (auto const& tableInWorkflow : tables) {
      LOG(info) << tableInWorkflow;
    }

    // Iterate through the tables and process based on the mapping
    for (auto const& table : tableToTasks) {
      bool foundIt = false;
      for (auto const& tableInWorkflow : tables) {
        if (tableInWorkflow == table.first) {
          foundIt = true;
          break;
        }
      }
      if (foundIt)
        continue;
      for (auto const& task : table.second) {
        LOG(info) << "Adding task " << table.first;
        task();
      }
    }
  } else {
    LOG(warning) << "AOD converter: No tables found in the meta data";
  }
  return workflow;
}
