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
// O2 includes

#include <iostream>
#include <cstdio>
#include <random>
#include <string>
#include <string_view>
#include <vector>

#include <fmt/format.h>
#include <rapidjson/document.h>
#include <rapidjson/filereadstream.h>

#include "filterTables.h"

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"
#include "Common/DataModel/TrackSelectionTables.h"

// we need to add workflow options before including Framework/runDataProcessing
void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  // option allowing to set parameters
  std::vector<o2::framework::ConfigParamSpec> options{o2::framework::ConfigParamSpec{"train_config", o2::framework::VariantType::String, "full_config.json", {"Configuration of the filtering train"}}};

  std::swap(workflowOptions, options);
}

#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace rapidjson;

namespace
{
bool readJsonFile(std::string& config, Document& d)
{
  FILE* fp = fopen(config.data(), "rb");
  if (!fp) {
    LOG(warning) << "Missing configuration json file: " << config;
    return false;
  }

  char readBuffer[65536];
  FileReadStream is(fp, readBuffer, sizeof(readBuffer));

  d.ParseStream(is);
  fclose(fp);
  return true;
}

std::unordered_map<std::string, std::unordered_map<std::string, float>> mDownscaling;
static const std::vector<std::string> downscalingName{"Downscaling"};
static const float defaultDownscaling[128][1]{
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f},
  {1.f}}; /// Max number of columns for triggers is 128 (extendible)

#define FILTER_CONFIGURABLE(_TYPE_) \
  Configurable<LabeledArray<float>> cfg##_TYPE_ { #_TYPE_, {defaultDownscaling[0], NumberOfColumns < _TYPE_>(), 1, ColumnsNames(typename _TYPE_::iterator::persistent_columns_t{}), downscalingName }, #_TYPE_ " downscalings" }

} // namespace

struct centralEventFilterTask {

  HistogramRegistry scalers{"scalers", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  Produces<aod::CefpDecisions> tags;
  Configurable<float> cfgTimingCut{"cfgTimingCut", 1.f, "nsigma timing cut associating BC and collisions"};

  FILTER_CONFIGURABLE(NucleiFilters);
  FILTER_CONFIGURABLE(DiffractionFilters);
  FILTER_CONFIGURABLE(DqFilters);
  FILTER_CONFIGURABLE(HfFilters);
  FILTER_CONFIGURABLE(CFFiltersTwoN);
  FILTER_CONFIGURABLE(CFFilters);
  FILTER_CONFIGURABLE(JetFilters);
  FILTER_CONFIGURABLE(StrangenessFilters);
  FILTER_CONFIGURABLE(MultFilters);

  void init(o2::framework::InitContext& initc)
  {
    LOG(info) << "Start init";
    int nCols{0};
    for (auto& table : mDownscaling) {
      nCols += table.second.size();
    }
    LOG(info) << "Middle init, total number of columns " << nCols;

    auto mScalers = std::get<std::shared_ptr<TH1>>(scalers.add("mScalers", ";;Number of events", HistType::kTH1F, {{nCols + 1, -0.5, 0.5 + nCols}}));
    auto mFiltered = std::get<std::shared_ptr<TH1>>(scalers.add("mFiltered", ";;Number of filtered events", HistType::kTH1F, {{nCols + 1, -0.5, 0.5 + nCols}}));

    mScalers->GetXaxis()->SetBinLabel(1, "Total number of events");
    mFiltered->GetXaxis()->SetBinLabel(1, "Total number of events");
    int bin{2};

    // for (auto& spec : reinterpret_cast<std::unique_ptr<ConfigParamStore>*>(&(initc.mOptions))->get()->specs()) {
    //   std::cout << "Configuration available: " << spec.name << "\t" << int(spec.type) << std::endl;
    //   auto filterOpt = initc.mOptions.get<LabeledArray<float>>(spec.name.data());
    //   std::cout << " -- row labs: ";
    //   for (auto& lab : filterOpt.labels_rows) {
    //     std::cout << lab << "\t";
    //   }
    //   std::cout << "\n -- col labs: ";
    //   for (auto& lab : filterOpt.labels_cols) {
    //     std::cout << lab << "\t";
    //   }
    //   std::cout << std::endl;
    // }

    for (auto& table : mDownscaling) {
      LOG(info) << "Setting downscalings for table " << table.first;
      for (auto& column : table.second) {
        mScalers->GetXaxis()->SetBinLabel(bin, column.first.data());
        mFiltered->GetXaxis()->SetBinLabel(bin++, column.first.data());
      }
      if (!initc.options().isSet(table.first.data())) {
        continue;
      }
      auto filterOpt = initc.mOptions.get<LabeledArray<float>>(table.first.data());
      for (auto& col : table.second) {
        LOG(info) << "- Channel " << col.first << ": " << filterOpt.get(col.first.data(), 0u);
        col.second = filterOpt.get(col.first.data(), 0u);
      }
    }
  }

  void run(ProcessingContext& pc)
  {
    auto mScalers{scalers.get<TH1>(HIST("mScalers"))};
    auto mFiltered{scalers.get<TH1>(HIST("mFiltered"))};

    int64_t nEvents{-1};
    std::vector<bool> outDecision;
    for (auto& tableName : mDownscaling) {
      if (!pc.inputs().isValid(tableName.first)) {
        LOG(fatal) << tableName.first << " table is not valid.";
      }
      auto tableConsumer = pc.inputs().get<TableConsumer>(tableName.first);
      auto tablePtr{tableConsumer->asArrowTable()};
      int64_t nRows{tablePtr->num_rows()};
      nEvents = nEvents < 0 ? nRows : nEvents;
      if (nEvents != nRows) {
        LOG(fatal) << "Inconsistent number of rows across trigger tables.";
      }

      if (outDecision.size() == 0)
        outDecision.resize(nEvents, false);

      auto schema{tablePtr->schema()};
      for (auto& colName : tableName.second) {
        int bin{mScalers->GetXaxis()->FindBin(colName.first.data())};
        double binCenter{mScalers->GetXaxis()->GetBinCenter(bin)};
        auto column{tablePtr->GetColumnByName(colName.first)};
        double downscaling{colName.second};
        if (column) {
          int entry = 0;
          for (int64_t iC{0}; iC < column->num_chunks(); ++iC) {
            auto chunk{column->chunk(iC)};
            auto boolArray = std::static_pointer_cast<arrow::BooleanArray>(chunk);
            for (int64_t iS{0}; iS < chunk->length(); ++iS) {
              if (boolArray->Value(iS)) {
                mScalers->Fill(binCenter);
                if (mUniformGenerator(mGeneratorEngine) < downscaling) {
                  mFiltered->Fill(binCenter);
                  outDecision[entry] = true;
                }
              }
              entry++;
            }
          }
        }
      }
    }
    mScalers->SetBinContent(1, mScalers->GetBinContent(1) + nEvents);
    mFiltered->SetBinContent(1, mFiltered->GetBinContent(1) + nEvents);

    // Filling output table
    auto collTabConsumer = pc.inputs().get<TableConsumer>("Collisions");
    auto collTabPtr{collTabConsumer->asArrowTable()};
    if (outDecision.size() != static_cast<uint64_t>(collTabPtr->num_rows())) {
      LOG(fatal) << "Inconsistent number of rows across Collision table and CEFP decision vector.";
    }
    auto columnBCId{collTabPtr->GetColumnByName("fIndexBCs")};
    auto columnCollTime{collTabPtr->GetColumnByName("fCollisionTime")};
    auto columnCollTimeRes{collTabPtr->GetColumnByName("fCollisionTimeRes")};

    std::unordered_map<int32_t, int64_t> decisions;

    for (int64_t iC{0}; iC < columnBCId->num_chunks(); ++iC) {
      auto chunkBC{columnBCId->chunk(iC)};
      auto chunkCollTime{columnCollTime->chunk(iC)};
      auto chunkCollTimeRes{columnCollTimeRes->chunk(iC)};

      auto BCArray = std::static_pointer_cast<arrow::NumericArray<arrow::Int32Type>>(chunkBC);
      auto CollTimeArray = std::static_pointer_cast<arrow::NumericArray<arrow::DoubleType>>(chunkCollTime);
      auto CollTimeResArray = std::static_pointer_cast<arrow::NumericArray<arrow::DoubleType>>(chunkCollTimeRes);
      for (int64_t iD{0}; iD < chunkBC->length(); ++iD) {
        auto collTime = CollTimeArray->Value(iD);
        auto collTimeRes = CollTimeResArray->Value(iD);
        int32_t startBC{BCArray->Value(iD) - static_cast<int>(std::floor(collTime - cfgTimingCut * collTimeRes))};
        int32_t endBC{BCArray->Value(iD) + static_cast<int>(std::ceil((collTime + cfgTimingCut * collTimeRes) / 25.f))};
        for (int32_t iB{startBC}; iB < endBC; ++iB) {
          if (decisions.find(iB) == decisions.end()) {
            decisions[iB] = static_cast<int64_t>(outDecision[iD]);
          } else {
            decisions[iB] |= static_cast<int64_t>(outDecision[iD]);
          }
        }
      }
    }

    for (auto& decision : decisions) {
      tags(decision.first, decision.second);
    }
  }

  std::mt19937_64 mGeneratorEngine;
  std::uniform_real_distribution<double> mUniformGenerator = std::uniform_real_distribution<double>(0., 1.);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{
  std::vector<InputSpec> inputs;
  inputs.emplace_back("Collisions", "AOD", "COLLISION", 0, Lifetime::Timeframe);

  auto config = cfg.options().get<std::string>("train_config");
  Document d;
  std::unordered_map<std::string, std::unordered_map<std::string, float>> downscalings;
  FillFiltersMap(FiltersPack, downscalings);

  std::array<bool, NumberOfFilters> enabledFilters = {false};
  if (readJsonFile(config, d)) {
    LOG(info) << "Reading workflow json... ";
    for (auto& workflow : d["workflows"].GetArray()) {
      LOG(info) << "- Reading workflow list... " << workflow["workflow_name"].GetString();
      for (uint32_t iFilter{0}; iFilter < NumberOfFilters; ++iFilter) {
        if (std::string_view(workflow["workflow_name"].GetString()) == std::string_view(FilteringTaskNames[iFilter])) {
          inputs.emplace_back(std::string(AvailableFilters[iFilter]), "AOD", FilterDescriptions[iFilter], 0, Lifetime::Timeframe);
          enabledFilters[iFilter] = true;
          LOG(info) << "    * Adding inputs from " << AvailableFilters[iFilter] << " to workflow";
          break;
        }
      }
    }
  }

  for (uint32_t iFilter{0}; iFilter < NumberOfFilters; ++iFilter) {
    if (!enabledFilters[iFilter]) {
      LOG(info) << std::string_view(AvailableFilters[iFilter]) << " not present in the configuration, removing it.";
      downscalings.erase(std::string(AvailableFilters[iFilter]));
    }
  }

  DataProcessorSpec spec{adaptAnalysisTask<centralEventFilterTask>(cfg)};
  for (auto& input : inputs) {
    spec.inputs.emplace_back(input);
  }
  mDownscaling.swap(downscalings);

  return WorkflowSpec{spec};
}
