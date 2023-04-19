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
/// \author Deependra Sharma, IITB, deependra.sharma@cern.ch

#ifndef PWGCF_FEMTOWORLD_CORE_FEMTOWORLDPAIRWITHCENTRALITY_H_
#define PWGCF_FEMTOWORLD_CORE_FEMTOWORLDPAIRWITHCENTRALITY_H_

#include <string>
#include <iostream>
#include <vector>
#include "Framework/HistogramRegistry.h"

using namespace o2;
using namespace o2::framework;

namespace o2::analysis::femtoWorld
{
class PairWithCentrality
{
 public:
  virtual ~PairWithCentrality() = default;
  /// @brief
  /// @tparam t1
  /// @param registry
  /// @param kstarbins
  /// @param centbins
  template <typename t1>
  void init(HistogramRegistry* registry, t1& kstarbins, t1& centbins)
  {
    PairWithCentralityRegistry = registry;
    AxisSpec kstarAxis = {kstarbins, "#it{k*} (GeV/#it{c})"};
    CentBins = centbins;
    CentBins.erase(CentBins.begin());

    for (int i = 0; i < static_cast<int>(CentBins.size() - 1); i++) {
      std::string HistTitle = "kstar_cent_" + std::to_string(CentBins[i]) + "_" + std::to_string(CentBins[i + 1]);
      std::string HistSuffix1 = static_cast<std::string>(HistSuffix[i]);
      std::string HistSuffix2 = static_cast<std::string>(HistSuffix[i + 1]);
      std::string HistName = "kstar_cent_" + HistSuffix1 + "_" + HistSuffix2;
      PairWithCentralityRegistry->add(HistName.c_str(), HistTitle.c_str(), HistType::kTH1F, {kstarAxis});
    }
    PairWithCentralityRegistry->add("Beyond_Max_Cent", "Beyond_Max_Cent", HistType::kTH1F, {kstarAxis});
  }

  /// @brief
  /// @tparam t1
  /// @param kstar_value
  /// @param centV0M_value
  template <typename t1>
  void fill(t1 kstar_value, t1 centV0M_value)
  {
    if (centV0M_value > CentBins[CentBins.size() - 1] || centV0M_value < CentBins[0]) {
      PairWithCentralityRegistry->fill(HIST("Beyond_Max_Cent"), kstar_value);
    } else if (centV0M_value <= CentBins[1]) {
      PairWithCentralityRegistry->fill(HIST("kstar_cent_0_1"), kstar_value);
    } else if (centV0M_value <= CentBins[2]) {
      PairWithCentralityRegistry->fill(HIST("kstar_cent_1_2"), kstar_value);
    } else if (centV0M_value <= CentBins[3]) {
      PairWithCentralityRegistry->fill(HIST("kstar_cent_2_3"), kstar_value);
    } else if (centV0M_value <= CentBins[4]) {
      PairWithCentralityRegistry->fill(HIST("kstar_cent_3_4"), kstar_value);
    } else if (centV0M_value <= CentBins[5]) {
      PairWithCentralityRegistry->fill(HIST("kstar_cent_4_5"), kstar_value);
    } else if (centV0M_value <= CentBins[6]) {
      PairWithCentralityRegistry->fill(HIST("kstar_cent_5_6"), kstar_value);
    } else if (centV0M_value <= CentBins[7]) {
      PairWithCentralityRegistry->fill(HIST("kstar_cent_6_7"), kstar_value);
    } else if (centV0M_value <= CentBins[8]) {
      PairWithCentralityRegistry->fill(HIST("kstar_cent_7_8"), kstar_value);
    } else if (centV0M_value <= CentBins[9]) {
      PairWithCentralityRegistry->fill(HIST("kstar_cent_8_9"), kstar_value);
    }
  }

 protected:
  HistogramRegistry* PairWithCentralityRegistry = nullptr;
  std::vector<double> CentBins;
  static constexpr std::string_view HistSuffix[10] = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9"};
};
} // namespace o2::analysis::femtoWorld

#endif // PWGCF_FEMTOWORLD_CORE_FEMTOWORLDPAIRWITHCENTRALITY_H_
