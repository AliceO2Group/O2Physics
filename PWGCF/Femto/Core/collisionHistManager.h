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
  kSphericity,
  // qa
  kPosX,
  kPosY,
  kPos,
  kOccupancy,
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
    {kMagField, o2::framework::kTH1F, "hMagField", "Magnetic Field; B (T); Entries"},
    {kSphericity, o2::framework::kTH1F, "hSphericity", "Sphericity; Sphericity; Entries"},
    {kPosX, o2::framework::kTH1F, "hPosX", "Vertex X; V_{X} (cm); Entries"},
    {kPosY, o2::framework::kTH1F, "hPosY", "Vertex Z; V_{Y} (cm); Entries"},
    {kPos, o2::framework::kTH1F, "hPos", "Primary vertex; V_{pos} (cm); Entries"},
    {kOccupancy, o2::framework::kTH1F, "hOccupancy", "Occupancy; Occupancy; Entries"},
    {kPoszVsMult, o2::framework::kTH2F, "hPoszVsMult", "Vertex Z vs Multiplicity; V_{Z} (cm); Multiplicity"},
    {kPoszVsCent, o2::framework::kTH2F, "hPoszVsCent", "Vertex Z vs Centrality; V_{Z} (cm); Centrality (%)"},
    {kCentVsMult, o2::framework::kTH2F, "hCentVsMult", "Centrality vs Multiplicity; Centrality (%); Multiplicity"},
    {kMultVsSphericity, o2::framework::kTH2F, "hMultVsSphericity", "Multiplicity vs Sphericity; Multiplicity; Sphericity"},
    {kCentVsSphericity, o2::framework::kTH2F, "hCentVsSphericity", "Centrality vs Sphericity; Centrality (%); Sphericity"},
  }};

template <typename BinningStruct>
auto makeColHistSpecMap(const BinningStruct& binning)
{
  return std::map<ColHist, std::vector<framework::AxisSpec>>{
    {kPosZ, {binning.vtZ}},
    {kMult, {binning.mult}},
    {kCent, {binning.cent}},
    {kSphericity, {binning.spher}},
    {kMagField, {binning.magField}},
    {kPoszVsMult, {binning.vtZ, binning.mult}},
    {kPoszVsCent, {binning.vtZ, binning.cent}},
    {kCentVsMult, {binning.cent, binning.mult}},
    {kMultVsSphericity, {binning.mult, binning.spher}},
    {kCentVsSphericity, {binning.cent, binning.spher}}};
}

struct ConfCollisionBinning : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("CollisionBinning");
  o2::framework::ConfigurableAxis vtZ{"vtZ", {200, -10, 10}, "Vertex Z binning"};
  o2::framework::ConfigurableAxis mult{"mult", {200, 0, 200}, "Multiplicity binning"};
  o2::framework::ConfigurableAxis cent{"cent", {100, 0.0f, 100.0f}, "Centrality (multiplicity percentile) binning"};
  o2::framework::ConfigurableAxis spher{"spher", {200, 0.0f, 2.0f}, "Sphericity binning"};
  o2::framework::ConfigurableAxis magField{"magField", {2, -1, 1}, "Magnetic field binning"};
};

template <modes::Mode mode>
class CollisionHistManager
{
 public:
  /// Destructor
  virtual ~CollisionHistManager() = default;
  /// Initializes histograms for the task
  /// \param registry Histogram registry to be passed
  void init(o2::framework::HistogramRegistry* registry, std::map<ColHist, std::vector<o2::framework::AxisSpec>> Specs)
  {
    mHistogramRegistry = registry;
    if constexpr (isFlagSet(mode, modes::Mode::kAnalysis)) {
      std::string analysisDir = std::string(ColAnalysisDir);
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kPosZ, HistTable), GetHistDesc(kPosZ, HistTable), GetHistType(kPosZ, HistTable), {Specs[kPosZ]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kMult, HistTable), GetHistDesc(kMult, HistTable), GetHistType(kMult, HistTable), {Specs[kMult]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kCent, HistTable), GetHistDesc(kCent, HistTable), GetHistType(kCent, HistTable), {Specs[kCent]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kSphericity, HistTable), GetHistDesc(kSphericity, HistTable), GetHistType(kSphericity, HistTable), {Specs[kSphericity]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kMagField, HistTable), GetHistDesc(kMagField, HistTable), GetHistType(kMagField, HistTable), {Specs[kMagField]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kPoszVsMult, HistTable), GetHistDesc(kPoszVsMult, HistTable), GetHistType(kPoszVsMult, HistTable), {Specs[kPoszVsMult]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kPoszVsCent, HistTable), GetHistDesc(kPoszVsCent, HistTable), GetHistType(kPoszVsCent, HistTable), {Specs[kPoszVsCent]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kCentVsMult, HistTable), GetHistDesc(kCentVsMult, HistTable), GetHistType(kCentVsMult, HistTable), {Specs[kCentVsMult]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kMultVsSphericity, HistTable), GetHistDesc(kMultVsSphericity, HistTable), GetHistType(kMultVsSphericity, HistTable), {Specs[kMultVsSphericity]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kCentVsSphericity, HistTable), GetHistDesc(kCentVsSphericity, HistTable), GetHistType(kCentVsSphericity, HistTable), {Specs[kCentVsSphericity]});
    }

    if constexpr (isFlagSet(mode, modes::Mode::kQa)) {
      std::string qaDir = std::string(ColQaDir);
      // to be implemented
    }
  } // namespace o2::analysis::femtounited

  template <typename T>
  void fill(T const& col)
  {
    if constexpr (isFlagSet(mode, modes::Mode::kAnalysis)) {
      mHistogramRegistry->fill(HIST(ColAnalysisDir) + HIST(GetHistName(kPosZ, HistTable)), col.posZ());
      mHistogramRegistry->fill(HIST(ColAnalysisDir) + HIST(GetHistName(kMult, HistTable)), col.mult());
      mHistogramRegistry->fill(HIST(ColAnalysisDir) + HIST(GetHistName(kCent, HistTable)), col.cent());
      mHistogramRegistry->fill(HIST(ColAnalysisDir) + HIST(GetHistName(kSphericity, HistTable)), col.sphericity());
      mHistogramRegistry->fill(HIST(ColAnalysisDir) + HIST(GetHistName(kMagField, HistTable)), col.magField());
      mHistogramRegistry->fill(HIST(ColAnalysisDir) + HIST(GetHistName(kPoszVsMult, HistTable)), col.posZ(), col.mult());
      mHistogramRegistry->fill(HIST(ColAnalysisDir) + HIST(GetHistName(kPoszVsCent, HistTable)), col.posZ(), col.cent());
      mHistogramRegistry->fill(HIST(ColAnalysisDir) + HIST(GetHistName(kCentVsMult, HistTable)), col.cent(), col.mult());
      mHistogramRegistry->fill(HIST(ColAnalysisDir) + HIST(GetHistName(kMultVsSphericity, HistTable)), col.mult(), col.sphericity());
      mHistogramRegistry->fill(HIST(ColAnalysisDir) + HIST(GetHistName(kCentVsSphericity, HistTable)), col.cent(), col.sphericity());
    }

    if constexpr (isFlagSet(mode, modes::Mode::kQa)) {
      // to be implemented
    }
  }

 private:
  o2::framework::HistogramRegistry* mHistogramRegistry;
}; // namespace femtounitedcolhistmanager
}; // namespace colhistmanager
}; // namespace o2::analysis::femto
#endif // PWGCF_FEMTO_CORE_COLLISIONHISTMANAGER_H_
