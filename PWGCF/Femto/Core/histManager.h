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

#include "Framework/HistogramSpec.h"

#include <string>
#include <string_view>

namespace o2::analysis::femto
{
namespace histmanager
{

template <typename Hist>
struct HistInfo {
  Hist hist;
  o2::framework::HistType histtype;
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
} // namespace histmanager
} // namespace o2::analysis::femto
#endif // PWGCF_FEMTO_CORE_HISTMANAGER_H_
