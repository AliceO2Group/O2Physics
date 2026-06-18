// Copyright 2019-2025 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file mcCollisionHistManager.h
/// \brief histogram manager for mc particle histograms
/// \author Anton Riedel, TU München, anton.riedel@cern.ch

#ifndef PWGCF_FEMTO_CORE_MCCOLLISIONHISTMANAGER_H_
#define PWGCF_FEMTO_CORE_MCCOLLISIONHISTMANAGER_H_

#include "PWGCF/Femto/Core/histManager.h"
#include "PWGCF/Femto/Core/modes.h"

#include <CommonConstants/MathConstants.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/Logger.h>

#include <TH1.h>

#include <array>
#include <map>
#include <string>
#include <string_view>
#include <vector>

namespace o2::analysis::femto
{
namespace mccollisionhistmanager
{

// enum for mc particle histograms
enum McCollisionHist {
  kPosZ,
  kMult,
  kCent,
  kPoszVsMult,
  kPoszVsCent,
  kCentVsMult,
  kMcCollisionHistLast
};

struct ConfMcCollisionBinning : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("McCollisionBinning");
  o2::framework::ConfigurableAxis vtxZ{"vtxZ", {200, -10, 10}, "Vertex Z binning"};
  o2::framework::ConfigurableAxis mult{"mult", {200, 0, 200}, "Multiplicity binning"};
  o2::framework::ConfigurableAxis cent{"cent", {100, 0.0f, 100.0f}, "Centrality (multiplicity percentile) binning"};
};

// the enum gives the correct index in the array
constexpr std::array<histmanager::HistInfo<McCollisionHist>, kMcCollisionHistLast>
  HistTable = {
    {
      {kPosZ, o2::framework::HistType::kTH1F, "hPosZ", "Vertex Z; V_{Z} (cm); Entries"},
      {kMult, o2::framework::HistType::kTH1F, "hMult", "Multiplicity; Multiplicity; Entries"},
      {kCent, o2::framework::HistType::kTH1F, "hCent", "Centrality; Centrality (%); Entries"},
      {kPoszVsMult, o2::framework::HistType::kTH2F, "hPoszVsMult", "Vertex Z vs Multiplicity; V_{Z} (cm); Multiplicity"},
      {kPoszVsCent, o2::framework::HistType::kTH2F, "hPoszVsCent", "Vertex Z vs Centrality; V_{Z} (cm); Centrality (%)"},
      {kCentVsMult, o2::framework::HistType::kTH2F, "hCentVsMult", "Centrality vs Multiplicity; Centrality (%); Multiplicity"},
    }};

template <typename T>
auto makeMcCollisionHistSpecMap(const T& confBinning)
{
  return std::map<McCollisionHist, std::vector<o2::framework::AxisSpec>>{
    {kPosZ, {confBinning.vtxZ}},
    {kMult, {confBinning.mult}},
    {kCent, {confBinning.cent}},
    {kPoszVsCent, {confBinning.vtxZ, confBinning.cent}},
    {kPoszVsMult, {confBinning.vtxZ, confBinning.mult}},
    {kCentVsMult, {confBinning.cent, confBinning.mult}},
  };
}

constexpr char PrefixMcCollsion[] = "McCollision/";

constexpr std::string_view McDir = "MC/";

template <const char* prefix>
class McCollisionHistManager
{
 public:
  McCollisionHistManager() = default;
  ~McCollisionHistManager() = default;

  // init for analysis and mc
  void init(o2::framework::HistogramRegistry* registry,
            std::map<McCollisionHist, std::vector<o2::framework::AxisSpec>> const& Specs)

  {
    mHistogramRegistry = registry;

    std::string mcDir = std::string(prefix) + std::string(McDir);
    mHistogramRegistry->add(mcDir + getHistNameV2(kPosZ, HistTable), getHistDesc(kPosZ, HistTable), getHistType(kPosZ, HistTable), {Specs.at(kPosZ)});
    mHistogramRegistry->add(mcDir + getHistNameV2(kMult, HistTable), getHistDesc(kMult, HistTable), getHistType(kMult, HistTable), {Specs.at(kMult)});
    mHistogramRegistry->add(mcDir + getHistNameV2(kCent, HistTable), getHistDesc(kCent, HistTable), getHistType(kCent, HistTable), {Specs.at(kCent)});

    mHistogramRegistry->add(mcDir + getHistNameV2(kPoszVsMult, HistTable), getHistDesc(kPoszVsMult, HistTable), getHistType(kPoszVsMult, HistTable), {Specs.at(kPoszVsMult)});
    mHistogramRegistry->add(mcDir + getHistNameV2(kPoszVsCent, HistTable), getHistDesc(kPoszVsCent, HistTable), getHistType(kPoszVsCent, HistTable), {Specs.at(kPoszVsCent)});
    mHistogramRegistry->add(mcDir + getHistNameV2(kCentVsMult, HistTable), getHistDesc(kCentVsMult, HistTable), getHistType(kCentVsMult, HistTable), {Specs.at(kCentVsMult)});
  }

  template <typename T1>
  void fill(T1 const& col)
  {
    mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kPosZ, HistTable)), col.posZ());
    mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kMult, HistTable)), col.mult());
    mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kCent, HistTable)), col.cent());

    mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kPoszVsMult, HistTable)), col.posZ(), col.mult());
    mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kPoszVsCent, HistTable)), col.posZ(), col.cent());
    mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kCentVsMult, HistTable)), col.cent(), col.mult());
  }

 private:
  o2::framework::HistogramRegistry* mHistogramRegistry = nullptr;
};
}; // namespace mccollisionhistmanager
}; // namespace o2::analysis::femto
#endif // PWGCF_FEMTO_CORE_MCCOLLISIONHISTMANAGER_H_
