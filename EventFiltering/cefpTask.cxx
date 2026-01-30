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

#include <fmt/format.h>
#include <rapidjson/document.h>
#include <rapidjson/filereadstream.h>

#include <cstdio>
#include <random>
#include <string>
#include <string_view>
#include <vector>
#include <unordered_map>
#include <utility>

#include "filterTables.h"

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "CommonConstants/LHCConstants.h"
#include "CommonDataFormat/InteractionRecord.h"
#include "DataFormatsCTP/Scalers.h"

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

#define FILTER_CONFIGURABLE(_TYPE_)                                                                                                                                                                                  \
  Configurable<LabeledArray<float>> cfg##_TYPE_                                                                                                                                                                      \
  {                                                                                                                                                                                                                  \
    #_TYPE_, {defaultDownscaling[0], NumberOfColumns(typename _TYPE_::table_t::persistent_columns_t{}), 1, ColumnsNames(typename _TYPE_::table_t::persistent_columns_t{}), downscalingName}, #_TYPE_ " downscalings" \
  }
} // namespace

struct centralEventFilterTask {

  HistogramRegistry scalers{"scalers", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  Produces<aod::CefpDecisions> tags;

  Configurable<bool> cfgDisableDownscalings{"cfgDisableDownscalings", false, "Disable downscalings"};
  Configurable<bool> cfgSkipUntriggeredEvents{"cfgSkipUntriggeredEvents", false, "Skip untriggered events"};

  FILTER_CONFIGURABLE(F1ProtonFilters);
  FILTER_CONFIGURABLE(DoublePhiFilters);
  FILTER_CONFIGURABLE(GlueballFilters);
  FILTER_CONFIGURABLE(NucleiFilters);
  FILTER_CONFIGURABLE(DiffractionFilters);
  FILTER_CONFIGURABLE(DqFilters);
  FILTER_CONFIGURABLE(HfFilters);
  FILTER_CONFIGURABLE(CFFilters);
  FILTER_CONFIGURABLE(JetFilters);
  FILTER_CONFIGURABLE(StrangenessFilters);
  FILTER_CONFIGURABLE(MultFilters);
  FILTER_CONFIGURABLE(FullJetFilters);
  FILTER_CONFIGURABLE(PhotonFilters);
  FILTER_CONFIGURABLE(HeavyNeutralMesonFilters);

  void init(o2::framework::InitContext& initc)
  {
    LOG(debug) << "Start init";
    int nCols{0};
    for (auto& table : mDownscaling) {
      nCols += table.second.size();
    }
    LOG(debug) << "Middle init, total number of columns " << nCols;

    auto mScalers = std::get<std::shared_ptr<TH1>>(scalers.add("mScalers", ";;Number of events", HistType::kTH1D, {{nCols + 2, -0.5, 1.5 + nCols}}));
    auto mFiltered = std::get<std::shared_ptr<TH1>>(scalers.add("mFiltered", ";;Number of filtered events", HistType::kTH1D, {{nCols + 2, -0.5, 1.5 + nCols}}));
    auto mCovariance = std::get<std::shared_ptr<TH2>>(scalers.add("mCovariance", "Selection covariance", HistType::kTH2D, {{nCols, -0.5, nCols - 0.5}, {nCols, -0.5, nCols - 0.5}}));

    mScalers->GetXaxis()->SetBinLabel(1, "Total number of events");
    mFiltered->GetXaxis()->SetBinLabel(1, "Total number of events");
    mScalers->GetXaxis()->SetBinLabel(nCols + 2, "Triggered events");
    mFiltered->GetXaxis()->SetBinLabel(nCols + 2, "Filtered events");
    int bin{2};

    for (auto& table : mDownscaling) {
      LOG(info) << "Setting downscalings for table " << table.first;
      for (auto& column : table.second) {
        mCovariance->GetXaxis()->SetBinLabel(bin - 1, column.first.data());
        mCovariance->GetYaxis()->SetBinLabel(bin - 1, column.first.data());
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
    if (cfgDisableDownscalings.value) {
      LOG(info) << "Downscalings are disabled for all channels.";
    }
  }

  void run(ProcessingContext& pc)
  {

    // Filling output table
    auto bcTabConsumer = pc.inputs().get<TableConsumer>(o2::soa::getTableLabel<aod::BCs>());
    auto bcTabPtr{bcTabConsumer->asArrowTable()};
    auto collTabConsumer = pc.inputs().get<TableConsumer>(o2::soa::getTableLabel<aod::Collisions>());
    auto collTabPtr{collTabConsumer->asArrowTable()};
    auto evSelConsumer = pc.inputs().get<TableConsumer>(o2::soa::getTableLabel<aod::EvSels>());
    auto evSelTabPtr{evSelConsumer->asArrowTable()};

    auto columnGloBCId{bcTabPtr->GetColumnByName(aod::BC::GlobalBC::mLabel)};
    auto columnCollBCId{collTabPtr->GetColumnByName(aod::Collision::BCId::mLabel)};
    auto columnCollTime{collTabPtr->GetColumnByName(aod::Collision::CollisionTime::mLabel)};
    auto columnCollTimeRes{collTabPtr->GetColumnByName(aod::Collision::CollisionTimeRes::mLabel)};
    auto columnFoundBC{evSelTabPtr->GetColumnByName(o2::aod::evsel::FoundBCId::mLabel)};

    auto chunkCollBCid{columnCollBCId->chunk(0)};
    auto chunkCollTime{columnCollTime->chunk(0)};
    auto chunkCollTimeRes{columnCollTimeRes->chunk(0)};
    auto chunkGloBC{columnGloBCId->chunk(0)};
    auto chunkFoundBC{columnFoundBC->chunk(0)};

    auto CollBCIdArray = std::static_pointer_cast<arrow::NumericArray<arrow::Int32Type>>(chunkCollBCid);
    auto CollTimeArray = std::static_pointer_cast<arrow::NumericArray<arrow::FloatType>>(chunkCollTime);
    auto CollTimeResArray = std::static_pointer_cast<arrow::NumericArray<arrow::FloatType>>(chunkCollTimeRes);
    auto GloBCArray = std::static_pointer_cast<arrow::NumericArray<arrow::UInt64Type>>(chunkGloBC);
    auto FoundBCArray = std::static_pointer_cast<arrow::NumericArray<arrow::Int32Type>>(chunkFoundBC);

    if (CollTimeArray->length() == 0) {
      LOG(warn) << "The collision table in the current folder is empty.";
    }

    int startCollision{0};

    auto mScalers{scalers.get<TH1>(HIST("mScalers"))};
    auto mFiltered{scalers.get<TH1>(HIST("mFiltered"))};
    auto mCovariance{scalers.get<TH2>(HIST("mCovariance"))};

    int64_t nEvents{collTabPtr->num_rows()};
    std::vector<std::array<uint64_t, 2>> outTrigger, outDecision;
    for (auto& tableName : mDownscaling) {
      if (!pc.inputs().isValid(tableName.first)) {
        LOG(fatal) << tableName.first << " table is not valid.";
      }
      auto tableConsumer = pc.inputs().get<TableConsumer>(tableName.first);
      auto tablePtr{tableConsumer->asArrowTable()};
      int64_t nRows{tablePtr->num_rows()};
      if (nEvents != nRows) {
        LOGF(fatal, "Inconsistent number of rows in the trigger table %s: %lld but it should be %lld", tableName.first.data(), nRows, nEvents);
      }

      if (outDecision.size() == 0) {
        outDecision.resize(nEvents, {0ull, 0ull});
        outTrigger.resize(nEvents, {0ull, 0ull});
      }

      auto schema{tablePtr->schema()};
      for (auto& colName : tableName.second) {
        uint64_t bin{static_cast<uint64_t>(mScalers->GetXaxis()->FindBin(colName.first.data()))};
        double binCenter{mScalers->GetXaxis()->GetBinCenter(bin)};
        uint64_t decisionBin{(bin - 2) / 64};
        uint64_t triggerBit{BIT((bin - 2) % 64)};
        auto column{tablePtr->GetColumnByName(colName.first)};
        double downscaling{cfgDisableDownscalings.value ? 1. : colName.second};
        if (column) {
          int entry = 0;
          for (int64_t iC{0}; iC < column->num_chunks(); ++iC) {
            auto chunk{column->chunk(iC)};
            auto boolArray = std::static_pointer_cast<arrow::BooleanArray>(chunk);
            for (int64_t iS{startCollision}; iS < chunk->length(); ++iS) {
              if (boolArray->Value(iS)) {
                mScalers->Fill(binCenter);
                outTrigger[entry][decisionBin] |= triggerBit;
                if (mUniformGenerator(mGeneratorEngine) < downscaling) {
                  mFiltered->Fill(binCenter);
                  outDecision[entry][decisionBin] |= triggerBit;
                }
              }
              entry++;
            }
          }
        }
      }
    }
    mScalers->SetBinContent(1, mScalers->GetBinContent(1) + nEvents - startCollision);
    mFiltered->SetBinContent(1, mFiltered->GetBinContent(1) + nEvents - startCollision);

    for (uint64_t iE{0}; iE < outTrigger.size(); ++iE) {
      const auto& triggerWord{outTrigger[iE]};
      bool triggered{false}, selected{false};
      for (uint64_t iD{0}; iD < triggerWord.size(); ++iD) {
        for (int iB{0}; iB < 64; ++iB) {
          if (!(triggerWord[iD] & BIT(iB))) {
            continue;
          }
          uint64_t xIndex{iD * 64 + iB};
          for (uint64_t jD{0}; jD < triggerWord.size(); ++jD) {
            for (int jB{0}; jB < 64; ++jB) {
              uint64_t yIndex{jD * 64 + jB};
              if (xIndex <= yIndex && triggerWord[jD] & BIT(jB)) {
                mCovariance->Fill(iD * 64 + iB, jD * 64 + jB);
              }
            }
          }
        }
        triggered = triggered || triggerWord[iD];
        selected = selected || outDecision[iE][iD];
      }
      if (triggered) {
        mScalers->Fill(mScalers->GetNbinsX() - 1);
      }
      if (selected) {
        mFiltered->Fill(mFiltered->GetNbinsX() - 1);
      }
    }

    if (outDecision.size() != static_cast<uint64_t>(nEvents)) {
      LOGF(fatal, "Inconsistent number of rows across Collision table and CEFP decision vector.");
    }
    if (outDecision.size() != static_cast<uint64_t>(evSelTabPtr->num_rows())) {
      LOGF(fatal, "Inconsistent number of rows across EvSel table and CEFP decision vector %ull vs %ll", outDecision.size(), evSelTabPtr->num_rows());
    }
    for (uint64_t iD{0}; iD < outDecision.size(); ++iD) {
      uint64_t foundBC = FoundBCArray->Value(iD) >= 0 && FoundBCArray->Value(iD) < GloBCArray->length() ? GloBCArray->Value(FoundBCArray->Value(iD)) : -1;
      if (cfgSkipUntriggeredEvents.value && !outDecision[iD][0] && !outDecision[iD][1]) {
        continue;
      }
      tags(CollBCIdArray->Value(iD), GloBCArray->Value(CollBCIdArray->Value(iD)), foundBC, CollTimeArray->Value(iD), CollTimeResArray->Value(iD), outTrigger[iD][0], outTrigger[iD][1], outDecision[iD][0], outDecision[iD][1]);
    }
  }

  using BCs = soa::Join<aod::BCs, aod::BcSels, aod::Run3MatchedToBCSparse>;
  using BC = BCs::iterator;
  using CCs = soa::Join<aod::Collisions, aod::EvSels>;
  using CC = CCs::iterator;
  void process(CCs const&, BCs const&)
  {
  }

  std::mt19937_64 mGeneratorEngine;
  std::uniform_real_distribution<double> mUniformGenerator = std::uniform_real_distribution<double>(0., 1.);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{
  std::vector<InputSpec> inputs;
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
