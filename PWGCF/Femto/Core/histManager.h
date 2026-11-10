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

/// \file histManager.h
/// \brief common structs for histogram managers
/// \author anton.riedel@tum.de, TU München, anton.riedel@tum.de

#ifndef PWGCF_FEMTO_CORE_HISTMANAGER_H_
#define PWGCF_FEMTO_CORE_HISTMANAGER_H_

#include <Framework/HistogramSpec.h>

#include <string>
#include <string_view>

namespace o2::analysis::femto::histmanager
{
// plot level for secondaries during mc processing
enum SecondaryPlotLevel {
  kSecondaryPlotLevel1 = 1,
  kSecondaryPlotLevel2 = 2,
  kSecondaryPlotLevel3 = 3
};

template <typename Hist>
struct HistInfo {
  constexpr HistInfo(Hist h,
                     o2::framework::HistType t,
                     std::string_view n,
                     std::string_view d)
    : hist(h), histtype(t), histname(n), histdesc(d)
  {
  }

  Hist hist{};
  o2::framework::HistType histtype{o2::framework::kUndefinedHist};
  std::string_view histname;
  std::string_view histdesc;
};

template <typename EnumType, typename ArrayType>
constexpr o2::framework::HistType getHistType(EnumType variable, const ArrayType& array)
{
  const auto it = std::find_if(array.begin(), array.end(), [=](const auto& entry) {
    return entry.hist == variable;
  });

  return it != array.end() ? it->histtype : o2::framework::kUndefinedHist;
}

template <typename EnumType, typename ArrayType>
constexpr std::string_view getHistName(EnumType variable, const ArrayType& array)
{
  auto it = std::find_if(array.begin(), array.end(), [=](const auto& entry) {
    return entry.hist == variable;
  });

  return (it != array.end()) ? it->histname : std::string_view{};
}

template <typename EnumType, typename ArrayType>
std::string getHistNameV2(EnumType variable, const ArrayType& array)
{
  auto it = std::find_if(array.begin(), array.end(), [=](const auto& entry) {
    return entry.hist == variable;
  });

  return (it != array.end()) ? std::string(it->histname) : std::string{};
}

template <typename EnumType, typename ArrayType>
constexpr const char* getHistDesc(EnumType variable, const ArrayType& array)
{
  auto it = std::find_if(array.begin(), array.end(), [=](const auto& entry) {
    return entry.hist == variable;
  });

  return it != array.end() ? it->histdesc.data() : "";
}

// Replace the x-axis segment of a "title; x axis; y axis[; z axis...]" histogram description
// with newTitle. Used to patch a static HistTable description at runtime once the actual axis
// quantity (e.g. which momentum type, which kinematic variable) is resolved from a Configurable.
inline std::string replaceFirstAxisTitle(std::string_view desc, std::string_view newTitle)
{
  std::string result(desc);
  const auto first = result.find(';');
  const auto second = (first == std::string::npos) ? std::string::npos : result.find(';', first + 1);
  if (second != std::string::npos) {
    result.replace(first + 1, second - first - 1, " " + std::string(newTitle) + " ");
  }
  return result;
}
} // namespace o2::analysis::femto::histmanager
#endif // PWGCF_FEMTO_CORE_HISTMANAGER_H_
