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
/// \author Alicja Płachta, WUT Warsaw, alicja.plachta.stud@pw.edu.pl

#ifndef PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSEPAIRWITHCENTMULTKT_H_
#define PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSEPAIRWITHCENTMULTKT_H_

#include <string>
#include <iostream>
#include <vector>
#include "Framework/HistogramRegistry.h"

using namespace o2;
using namespace o2::framework;

namespace o2::analysis::femtoUniverse
{
class PairWithCentMultKt
{
 public:
  virtual ~PairWithCentMultKt() = default;
  /// @brief
  /// @tparam t1
  /// @param registry
  /// @param kstarbins
  /// @param centmultbins
  template <typename t1>
  void init(HistogramRegistry* registry, t1& kstarbins, t1& centmultbins, t1& ktbins, bool processKT)
  {
    PairWithCentMultKtRegistry = registry;
    AxisSpec kstarAxis = {kstarbins, "#it{k*} (GeV/#it{c})"};
    CentMultBins = centmultbins;
    KtBins = ktbins;
    KtBins.erase(KtBins.begin());
    CentMultBins.erase(CentMultBins.begin());
    UseKt = processKT;

    for (int i = 0; i < static_cast<int>(CentMultBins.size() - 1); i++) {
      int lowBin = static_cast<int>((CentMultBins[i]));
      int highBin = static_cast<int>((CentMultBins[i + 1]));
      std::string HistTitle = "mult_" + std::to_string(lowBin) + "-" + std::to_string(highBin);
      std::string HistSuffix1 = static_cast<std::string>(HistSuffix[i]);
      std::string HistSuffix2 = static_cast<std::string>(HistSuffix[i + 1]);
      std::string HistFolderMult = "mult_" + HistSuffix1 + "_" + HistSuffix2;
      std::string HistName = HistFolderMult + "/kstar";
      PairWithCentMultKtRegistry->add(HistName.c_str(), HistTitle.c_str(), HistType::kTH1F, {kstarAxis});
      if (UseKt) {
        for (int i = 0; i < static_cast<int>(KtBins.size() - 1); i++) {
          std::string kt_bin1_string = std::to_string(KtBins[i]);
          std::replace(kt_bin1_string.begin(), kt_bin1_string.end(), '.', '_');
          std::string kt_bin2_string = std::to_string(KtBins[i + 1]);
          std::replace(kt_bin2_string.begin(), kt_bin2_string.end(), '.', '_');
          kt_bin1_string.resize(4);
          kt_bin2_string.resize(4);
          std::string HistTitleKt = "kt_" + kt_bin1_string + "-" + kt_bin2_string;
          std::string HistSuffix1Kt = static_cast<std::string>(HistSuffix[i]);
          std::string HistSuffix2Kt = static_cast<std::string>(HistSuffix[i + 1]);
          std::string HistNameKt = HistFolderMult + "/kstar_kt_" + HistSuffix1Kt + "_" + HistSuffix2Kt;
          PairWithCentMultKtRegistry->add(HistNameKt.c_str(), HistTitleKt.c_str(), HistType::kTH1F, {kstarAxis});
        }
      }
    }
    PairWithCentMultKtRegistry->add("Beyond_Max", "Beyond_Max", HistType::kTH1F, {kstarAxis});
  }

  /// @brief
  /// @tparam t1
  /// @param kstar_value
  /// @param cent_mult_value
  template <typename t1>
  void fill(t1 kstar_value, t1 cent_mult_value, t1 kt_value)
  {

    if (cent_mult_value > CentMultBins[CentMultBins.size() - 1] || cent_mult_value < CentMultBins[0]) {
      PairWithCentMultKtRegistry->fill(HIST("Beyond_Max"), kstar_value);
    } else if (cent_mult_value <= CentMultBins[1]) {
      PairWithCentMultKtRegistry->fill(HIST("mult_0_1/kstar"), kstar_value);
      if (UseKt) {
        auto histMultFolder = HIST("mult_0_1/");
        fill_kT(kstar_value, kt_value, histMultFolder);
      }
    } else if (cent_mult_value <= CentMultBins[2]) {
      PairWithCentMultKtRegistry->fill(HIST("mult_1_2/kstar"), kstar_value);
      if (UseKt) {
        auto histMultFolder = HIST("mult_1_2/");
        fill_kT(kstar_value, kt_value, histMultFolder);
      }
    } else if (cent_mult_value <= CentMultBins[3]) {
      PairWithCentMultKtRegistry->fill(HIST("mult_2_3/kstar"), kstar_value);
      if (UseKt) {
        auto histMultFolder = HIST("mult_2_3/");
        fill_kT(kstar_value, kt_value, histMultFolder);
      }
    } else if (cent_mult_value <= CentMultBins[4]) {
      PairWithCentMultKtRegistry->fill(HIST("mult_3_4/kstar"), kstar_value);
      if (UseKt) {
        auto histMultFolder = HIST("mult_3_4/");
        fill_kT(kstar_value, kt_value, histMultFolder);
      }
    } else if (cent_mult_value <= CentMultBins[5]) {
      PairWithCentMultKtRegistry->fill(HIST("mult_4_5/kstar"), kstar_value);
      if (UseKt) {
        auto histMultFolder = HIST("mult_4_5/");
        fill_kT(kstar_value, kt_value, histMultFolder);
      }
    } else if (cent_mult_value <= CentMultBins[6]) {
      PairWithCentMultKtRegistry->fill(HIST("mult_5_6/kstar"), kstar_value);
      if (UseKt) {
        auto histMultFolder = HIST("mult_5_6/");
        fill_kT(kstar_value, kt_value, histMultFolder);
      }
    } else if (cent_mult_value <= CentMultBins[7]) {
      PairWithCentMultKtRegistry->fill(HIST("mult_6_7/kstar"), kstar_value);
      if (UseKt) {
        auto histMultFolder = HIST("mult_6_7/");
        fill_kT(kstar_value, kt_value, histMultFolder);
      }
    } else if (cent_mult_value <= CentMultBins[8]) {
      PairWithCentMultKtRegistry->fill(HIST("mult_7_8/kstar"), kstar_value);
      if (UseKt) {
        auto histMultFolder = HIST("mult_7_8/");
        fill_kT(kstar_value, kt_value, histMultFolder);
      }
    } else if (cent_mult_value <= CentMultBins[9]) {
      PairWithCentMultKtRegistry->fill(HIST("mult_8_9/kstar"), kstar_value);
      if (UseKt) {
        auto histMultFolder = HIST("mult_8_9/");
        fill_kT(kstar_value, kt_value, histMultFolder);
      }
    }
  }

  /// @brief
  /// @tparam t1
  /// @tparam t2
  /// @param kstar_value
  /// @param kt_value
  /// @param folder
  template <typename t1, typename t2>
  void fill_kT(t1 kstar_value, t1 kt_value, t2 folder)
  {
    if (kt_value <= KtBins[1]) {
      PairWithCentMultKtRegistry->fill(folder + HIST("kstar_kt_0_1"), kstar_value);
    } else if (kt_value <= KtBins[2]) {
      PairWithCentMultKtRegistry->fill(folder + HIST("kstar_kt_1_2"), kstar_value);
    } else if (kt_value <= KtBins[3]) {
      PairWithCentMultKtRegistry->fill(folder + HIST("kstar_kt_2_3"), kstar_value);
    } else if (kt_value <= KtBins[4]) {
      PairWithCentMultKtRegistry->fill(folder + HIST("kstar_kt_3_4"), kstar_value);
    } else if (kt_value <= KtBins[5]) {
      PairWithCentMultKtRegistry->fill(folder + HIST("kstar_kt_4_5"), kstar_value);
    } else if (kt_value <= KtBins[6]) {
      PairWithCentMultKtRegistry->fill(folder + HIST("kstar_kt_5_6"), kstar_value);
    } else if (kt_value <= KtBins[7]) {
      PairWithCentMultKtRegistry->fill(folder + HIST("kstar_kt_6_7"), kstar_value);
    } else if (kt_value <= KtBins[8]) {
      PairWithCentMultKtRegistry->fill(folder + HIST("kstar_kt_7_8"), kstar_value);
    } else if (kt_value <= KtBins[9]) {
      PairWithCentMultKtRegistry->fill(folder + HIST("kstar_kt_8_9"), kstar_value);
    }
  }

 protected:
  HistogramRegistry* PairWithCentMultKtRegistry = nullptr;
  std::vector<double> CentMultBins;
  std::vector<double> KtBins;
  bool UseKt = false;
  static constexpr std::string_view HistSuffix[10] = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9"};
};
} // namespace o2::analysis::femtoUniverse

#endif // PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSEPAIRWITHCENTMULTKT_H_
