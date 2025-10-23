// Copyright 2019-2022 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file collisionHistManager.h
/// \brief collision histogram manager
/// \author anton.riedel@tum.de, TU München, anton.riedel@tum.de

#ifndef PWGCF_FEMTO_CORE_COLLISIONHISTMANAGER_H_
#define PWGCF_FEMTO_CORE_COLLISIONHISTMANAGER_H_

#include "PWGCF/Femto/Core/histManager.h"
#include "PWGCF/Femto/Core/modes.h"

#include "Framework/Configurable.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/HistogramSpec.h"

#include <array>
#include <map>
#include <string>
#include <string_view>
#include <vector>

namespace o2::analysis::femto
{
namespace colhistmanager
{

enum ColHist {
  kPosZ,
  kMult,
  kCent,
  kMagField,
  // qa
  kPosX,
  kPosY,
  kPos,
  kOccupancy,
  kSphericity,
  // 2d
  kPoszVsMult,
  kPoszVsCent,
  kCentVsMult,
  kCentVsSphericity,
  kMultVsSphericity,
  kColHistLast
};

constexpr std::string_view ColAnalysisDir = "Collisions/Analysis/";
constexpr std::string_view ColQaDir = "Collisions/QA/";

constexpr std::array<histmanager::HistInfo<ColHist>, kColHistLast> HistTable = {
  {
    {kPosZ, o2::framework::kTH1F, "hPosZ", "Vertex Z; V_{Z} (cm); Entries"},
    {kMult, o2::framework::kTH1F, "hMult", "Multiplicity; Multiplicity; Entries"},
    {kCent, o2::framework::kTH1F, "hCent", "Centrality; Centrality (%); Entries"},
    {kMagField, o2::framework::kTH1F, "hMagField", "Magnetic Field; B (kG); Entries"},
    {kPosX, o2::framework::kTH1F, "hPosX", "Vertex X; V_{X} (cm); Entries"},
    {kPosY, o2::framework::kTH1F, "hPosY", "Vertex Z; V_{Y} (cm); Entries"},
    {kPos, o2::framework::kTH1F, "hPos", "Primary vertex; V_{pos} (cm); Entries"},
    {kSphericity, o2::framework::kTH1F, "hSphericity", "Sphericity; Sphericity; Entries"},
    {kOccupancy, o2::framework::kTH1F, "hOccupancy", "Occupancy; Occupancy; Entries"},
    {kPoszVsMult, o2::framework::kTH2F, "hPoszVsMult", "Vertex Z vs Multiplicity; V_{Z} (cm); Multiplicity"},
    {kPoszVsCent, o2::framework::kTH2F, "hPoszVsCent", "Vertex Z vs Centrality; V_{Z} (cm); Centrality (%)"},
    {kCentVsMult, o2::framework::kTH2F, "hCentVsMult", "Centrality vs Multiplicity; Centrality (%); Multiplicity"},
    {kMultVsSphericity, o2::framework::kTH2F, "hMultVsSphericity", "Multiplicity vs Sphericity; Multiplicity; Sphericity"},
    {kCentVsSphericity, o2::framework::kTH2F, "hCentVsSphericity", "Centrality vs Sphericity; Centrality (%); Sphericity"},
  }};

template <typename T>
auto makeColHistSpecMap(const T& confBinningAnalysis)
{
  return std::map<ColHist, std::vector<framework::AxisSpec>>{
    {kPosZ, {confBinningAnalysis.vtxZ}},
    {kMult, {confBinningAnalysis.mult}},
    {kCent, {confBinningAnalysis.cent}},
    {kMagField, {confBinningAnalysis.magField}}};
}

template <typename T1, typename T2>
auto makeColQaHistSpecMap(const T1& confBinningAnalysis, const T2& confBinningQa)
{
  return std::map<ColHist, std::vector<framework::AxisSpec>>{
    {kPosZ, {confBinningAnalysis.vtxZ}},
    {kMult, {confBinningAnalysis.mult}},
    {kCent, {confBinningAnalysis.cent}},
    {kMagField, {confBinningAnalysis.magField}},
    {kPosX, {confBinningQa.vtxXY}},
    {kPosY, {confBinningQa.vtxXY}},
    {kPos, {confBinningQa.vtx}},
    {kSphericity, {confBinningQa.sphericity}},
    {kOccupancy, {confBinningQa.occupancy}},
    {kPoszVsMult, {confBinningAnalysis.vtxZ, confBinningAnalysis.mult}},
    {kPoszVsCent, {confBinningAnalysis.vtxZ, confBinningAnalysis.cent}},
    {kCentVsMult, {confBinningAnalysis.cent, confBinningAnalysis.mult}},
    {kMultVsSphericity, {confBinningAnalysis.mult, confBinningQa.sphericity}},
    {kCentVsSphericity, {confBinningAnalysis.cent, confBinningQa.sphericity}}};
}

struct ConfCollisionBinning : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("CollisionBinning");
  o2::framework::ConfigurableAxis vtxZ{"vtxZ", {200, -10, 10}, "Vertex Z binning"};
  o2::framework::ConfigurableAxis mult{"mult", {200, 0, 200}, "Multiplicity binning"};
  o2::framework::ConfigurableAxis cent{"cent", {100, 0.0f, 100.0f}, "Centrality (multiplicity percentile) binning"};
  o2::framework::ConfigurableAxis magField{"magField", {11, -5.5, 5.5}, "Magnetic field binning"};
};

struct ConfCollisionQaBinning : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("CollisionQaBinning");
  o2::framework::ConfigurableAxis vtx{"vtx", {120, 0.f, 12.f}, "Vertex position binning"};
  o2::framework::ConfigurableAxis vtxXY{"vtxXY", {100, -1.f, 1.f}, "Vertex X/Y binning"};
  o2::framework::ConfigurableAxis sphericity{"sphericity", {100, 0.f, 1.f}, "Spericity Binning"};
  o2::framework::ConfigurableAxis occupancy{"occupancy", {500, 0.f, 5000.f}, "Spericity Binning"};
};

template <modes::Mode mode>
class CollisionHistManager
{
 public:
  CollisionHistManager() = default;
  virtual ~CollisionHistManager() = default;
  /// Initializes histograms for the task
  /// \param registry Histogram registry to be passed
  void init(o2::framework::HistogramRegistry* registry, std::map<ColHist, std::vector<o2::framework::AxisSpec>> Specs)
  {
    mHistogramRegistry = registry;
    if constexpr (isFlagSet(mode, modes::Mode::kAnalysis)) {
      std::string analysisDir = std::string(ColAnalysisDir);
      mHistogramRegistry->add(analysisDir + getHistNameV2(kPosZ, HistTable), getHistDesc(kPosZ, HistTable), getHistType(kPosZ, HistTable), {Specs[kPosZ]});
      mHistogramRegistry->add(analysisDir + getHistNameV2(kMult, HistTable), getHistDesc(kMult, HistTable), getHistType(kMult, HistTable), {Specs[kMult]});
      mHistogramRegistry->add(analysisDir + getHistNameV2(kCent, HistTable), getHistDesc(kCent, HistTable), getHistType(kCent, HistTable), {Specs[kCent]});
      mHistogramRegistry->add(analysisDir + getHistNameV2(kMagField, HistTable), getHistDesc(kMagField, HistTable), getHistType(kMagField, HistTable), {Specs[kMagField]});
    }

    if constexpr (isFlagSet(mode, modes::Mode::kQa)) {
      std::string qaDir = std::string(ColQaDir);

      mHistogramRegistry->add(qaDir + getHistNameV2(kPosX, HistTable), getHistDesc(kPosX, HistTable), getHistType(kPosX, HistTable), {Specs[kPosX]});
      mHistogramRegistry->add(qaDir + getHistNameV2(kPosY, HistTable), getHistDesc(kPosY, HistTable), getHistType(kPosY, HistTable), {Specs[kPosY]});
      mHistogramRegistry->add(qaDir + getHistNameV2(kPos, HistTable), getHistDesc(kPos, HistTable), getHistType(kPos, HistTable), {Specs[kPos]});
      mHistogramRegistry->add(qaDir + getHistNameV2(kSphericity, HistTable), getHistDesc(kSphericity, HistTable), getHistType(kSphericity, HistTable), {Specs[kSphericity]});
      mHistogramRegistry->add(qaDir + getHistNameV2(kOccupancy, HistTable), getHistDesc(kOccupancy, HistTable), getHistType(kOccupancy, HistTable), {Specs[kOccupancy]});
      mHistogramRegistry->add(qaDir + getHistNameV2(kPoszVsMult, HistTable), getHistDesc(kPoszVsMult, HistTable), getHistType(kPoszVsMult, HistTable), {Specs[kPoszVsMult]});
      mHistogramRegistry->add(qaDir + getHistNameV2(kPoszVsCent, HistTable), getHistDesc(kPoszVsCent, HistTable), getHistType(kPoszVsCent, HistTable), {Specs[kPoszVsCent]});
      mHistogramRegistry->add(qaDir + getHistNameV2(kCentVsMult, HistTable), getHistDesc(kCentVsMult, HistTable), getHistType(kCentVsMult, HistTable), {Specs[kCentVsMult]});
      mHistogramRegistry->add(qaDir + getHistNameV2(kMultVsSphericity, HistTable), getHistDesc(kMultVsSphericity, HistTable), getHistType(kMultVsSphericity, HistTable), {Specs[kMultVsSphericity]});
      mHistogramRegistry->add(qaDir + getHistNameV2(kCentVsSphericity, HistTable), getHistDesc(kCentVsSphericity, HistTable), getHistType(kCentVsSphericity, HistTable), {Specs[kCentVsSphericity]});
    }
  } // namespace o2::analysis::femtounited

  template <typename T>
  void fill(T const& col)
  {
    if constexpr (isFlagSet(mode, modes::Mode::kAnalysis)) {
      mHistogramRegistry->fill(HIST(ColAnalysisDir) + HIST(getHistName(kPosZ, HistTable)), col.posZ());
      mHistogramRegistry->fill(HIST(ColAnalysisDir) + HIST(getHistName(kMult, HistTable)), col.mult());
      mHistogramRegistry->fill(HIST(ColAnalysisDir) + HIST(getHistName(kCent, HistTable)), col.cent());
      mHistogramRegistry->fill(HIST(ColAnalysisDir) + HIST(getHistName(kMagField, HistTable)), col.magField());
    }

    if constexpr (isFlagSet(mode, modes::Mode::kQa)) {
      mHistogramRegistry->fill(HIST(ColQaDir) + HIST(getHistName(kPosX, HistTable)), col.posX());
      mHistogramRegistry->fill(HIST(ColQaDir) + HIST(getHistName(kPosY, HistTable)), col.posY());
      mHistogramRegistry->fill(HIST(ColQaDir) + HIST(getHistName(kPos, HistTable)), std::hypot(col.posX(), col.posY(), col.posZ()));
      mHistogramRegistry->fill(HIST(ColQaDir) + HIST(getHistName(kSphericity, HistTable)), col.sphericity());
      mHistogramRegistry->fill(HIST(ColQaDir) + HIST(getHistName(kOccupancy, HistTable)), col.trackOccupancyInTimeRange());
      mHistogramRegistry->fill(HIST(ColQaDir) + HIST(getHistName(kPoszVsMult, HistTable)), col.posZ(), col.mult());
      mHistogramRegistry->fill(HIST(ColQaDir) + HIST(getHistName(kPoszVsCent, HistTable)), col.posZ(), col.cent());
      mHistogramRegistry->fill(HIST(ColQaDir) + HIST(getHistName(kCentVsMult, HistTable)), col.cent(), col.mult());
      mHistogramRegistry->fill(HIST(ColQaDir) + HIST(getHistName(kMultVsSphericity, HistTable)), col.mult(), col.sphericity());
      mHistogramRegistry->fill(HIST(ColQaDir) + HIST(getHistName(kCentVsSphericity, HistTable)), col.cent(), col.sphericity());
    }
  }

 private:
  o2::framework::HistogramRegistry* mHistogramRegistry = nullptr;
}; // namespace femtounitedcolhistmanager
}; // namespace colhistmanager
}; // namespace o2::analysis::femto
#endif // PWGCF_FEMTO_CORE_COLLISIONHISTMANAGER_H_
