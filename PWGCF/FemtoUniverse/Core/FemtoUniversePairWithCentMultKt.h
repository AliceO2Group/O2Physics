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
/// \file FemtoUniversePairWithCentMultKt.h
/// \brief FemtoUniversePairWithCentMultKt - Histogram class for tracks with centrality and multiplicity
/// \author Deependra Sharma, IITB, deependra.sharma@cern.ch
/// \author Alicja Płachta, WUT Warsaw, alicja.plachta.stud@pw.edu.pl
/// \author Zuzanna Chochulska, WUT Warsaw & CTU Prague, zchochul@cern.ch

#ifndef PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSEPAIRWITHCENTMULTKT_H_
#define PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSEPAIRWITHCENTMULTKT_H_

#include <string>
#include <vector>
#include "Framework/HistogramRegistry.h"

namespace o2::analysis::femto_universe
{
class FemtoUniversePairWithCentMultKt
{
 public:
  virtual ~FemtoUniversePairWithCentMultKt() = default;
  /// @brief
  /// @tparam t1
  /// @param registry
  /// @param kstarbins
  /// @param centmultbins
  template <typename t1>
  void init(HistogramRegistry* registry, t1& kstarbins, t1& centmultbins, t1& ktbins, bool processKT, bool process3D)
  {
    pairWithCentMultKtRegistry = registry;
    AxisSpec kstarAxis = {kstarbins, "#it{k*} (GeV/#it{c})"};
    AxisSpec kOutAxis = {kstarbins, "#it{q}_{out} (GeV/#it{c})"};
    AxisSpec kSideAxis = {kstarbins, "#it{q}_{side} (GeV/#it{c})"};
    AxisSpec kLongAxis = {kstarbins, "#it{q}_{long} (GeV/#it{c})"};
    centMultBins = centmultbins;
    ktBins = ktbins;
    ktBins.erase(ktBins.begin());
    centMultBins.erase(centMultBins.begin());
    useKt = processKT;
    use3D = process3D;

    for (int i = 0; i < static_cast<int>(centMultBins.size() - 1); i++) {
      int lowBin = static_cast<int>((centMultBins[i]));
      int highBin = static_cast<int>((centMultBins[i + 1]));
      std::string histTitle = "mult_" + std::to_string(lowBin) + "-" + std::to_string(highBin);
      std::string histSuffix1 = static_cast<std::string>(HistSuffix[i]);
      std::string histSuffix2 = static_cast<std::string>(HistSuffix[i + 1]);
      LOGF(info, "histSuffix1 %s histSuffix2", histSuffix1, histSuffix2);
      std::string histFolderMult = "mult_" + histSuffix1 + "_" + histSuffix2;
      std::string histName = histFolderMult + "/kstar";
      std::string histName3D = histFolderMult + "/q3D";
      pairWithCentMultKtRegistry->add(histName.c_str(), histTitle.c_str(), HistType::kTH1F, {kstarAxis});
      pairWithCentMultKtRegistry->add(histName3D.c_str(), histTitle.c_str(), HistType::kTH3F, {kOutAxis, kSideAxis, kLongAxis});
      if (useKt) {
        for (int i = 0; i < static_cast<int>(ktBins.size() - 1); i++) {
          std::string ktBin1String = std::to_string(ktBins[i]);
          std::replace(ktBin1String.begin(), ktBin1String.end(), '.', '_');
          std::string ktBin2String = std::to_string(ktBins[i + 1]);
          std::replace(ktBin2String.begin(), ktBin2String.end(), '.', '_');
          ktBin1String.resize(4);
          ktBin2String.resize(4);
          std::string histTitleKt = "kt_" + ktBin1String + "-" + ktBin2String;
          std::string histSuffix1Kt = static_cast<std::string>(HistSuffix[i]);
          std::string histSuffix2Kt = static_cast<std::string>(HistSuffix[i + 1]);
          std::string histNameKt = histFolderMult + "/kstar_kt_" + histSuffix1Kt + "_" + histSuffix2Kt;
          LOGF(info, "histNameKt %s", histNameKt);
          pairWithCentMultKtRegistry->add(histNameKt.c_str(), histTitleKt.c_str(), HistType::kTH1F, {kstarAxis});
        }
      }
      if (use3D) {
        for (int i = 0; i < static_cast<int>(ktBins.size() - 1); i++) {
          std::string ktBin1String = std::to_string(ktBins[i]);
          std::replace(ktBin1String.begin(), ktBin1String.end(), '.', '_');
          std::string ktBin2String = std::to_string(ktBins[i + 1]);
          std::replace(ktBin2String.begin(), ktBin2String.end(), '.', '_');
          ktBin1String.resize(4);
          ktBin2String.resize(4);
          std::string histTitleKt = "kt_" + ktBin1String + "-" + ktBin2String;
          std::string histSuffix1Kt = static_cast<std::string>(HistSuffix[i]);
          std::string histSuffix2Kt = static_cast<std::string>(HistSuffix[i + 1]);
          std::string histNameKt = histFolderMult + "/q3D_kt_" + histSuffix1Kt + "_" + histSuffix2Kt;
          LOGF(info, "histNameKt %s", histNameKt);
          pairWithCentMultKtRegistry->add(histNameKt.c_str(), histTitleKt.c_str(), HistType::kTH3F, {kOutAxis, kSideAxis, kLongAxis});
        }
      }
    }
    pairWithCentMultKtRegistry->add("Beyond_Max", "Beyond_Max", HistType::kTH1F, {kstarAxis});
    pairWithCentMultKtRegistry->add("Beyond_Max_3D", "Beyond_Max_3D", HistType::kTH3F, {kOutAxis, kSideAxis, kLongAxis});
  }

  /// @brief
  /// @tparam t1
  /// @param kstar_value
  /// @param cent_mult_value
  template <typename t1>
  void fill(t1 kstar_value, t1 cent_mult_value, t1 kt_value)
  {

    if (cent_mult_value >= centMultBins[centMultBins.size() - 1] || cent_mult_value < centMultBins[0]) {
      pairWithCentMultKtRegistry->fill(HIST("Beyond_Max"), kstar_value);
    } else if (cent_mult_value >= centMultBins[0] && cent_mult_value < centMultBins[1]) {
      pairWithCentMultKtRegistry->fill(HIST("mult_0_1/kstar"), kstar_value);
      if (useKt) {
        auto histMultFolder = HIST("mult_0_1/");
        fillkT(kstar_value, kt_value, histMultFolder);
      }
    } else if (cent_mult_value >= centMultBins[1] && cent_mult_value < centMultBins[2]) {
      pairWithCentMultKtRegistry->fill(HIST("mult_1_2/kstar"), kstar_value);
      if (useKt) {
        auto histMultFolder = HIST("mult_1_2/");
        fillkT(kstar_value, kt_value, histMultFolder);
      }
    } else if (cent_mult_value >= centMultBins[2] && cent_mult_value < centMultBins[3]) {
      pairWithCentMultKtRegistry->fill(HIST("mult_2_3/kstar"), kstar_value);
      if (useKt) {
        auto histMultFolder = HIST("mult_2_3/");
        fillkT(kstar_value, kt_value, histMultFolder);
      }
    } else if (cent_mult_value >= centMultBins[3] && cent_mult_value < centMultBins[4]) {
      pairWithCentMultKtRegistry->fill(HIST("mult_3_4/kstar"), kstar_value);
      if (useKt) {
        auto histMultFolder = HIST("mult_3_4/");
        fillkT(kstar_value, kt_value, histMultFolder);
      }
    } else if (cent_mult_value >= centMultBins[4] && cent_mult_value < centMultBins[5]) {
      pairWithCentMultKtRegistry->fill(HIST("mult_4_5/kstar"), kstar_value);
      if (useKt) {
        auto histMultFolder = HIST("mult_4_5/");
        fillkT(kstar_value, kt_value, histMultFolder);
      }
    } else if (cent_mult_value >= centMultBins[5] && cent_mult_value < centMultBins[6]) {
      pairWithCentMultKtRegistry->fill(HIST("mult_5_6/kstar"), kstar_value);
      if (useKt) {
        auto histMultFolder = HIST("mult_5_6/");
        fillkT(kstar_value, kt_value, histMultFolder);
      }
    } else if (cent_mult_value >= centMultBins[6] && cent_mult_value < centMultBins[7]) {
      pairWithCentMultKtRegistry->fill(HIST("mult_6_7/kstar"), kstar_value);
      if (useKt) {
        auto histMultFolder = HIST("mult_6_7/");
        fillkT(kstar_value, kt_value, histMultFolder);
      }
    } else if (cent_mult_value >= centMultBins[7] && cent_mult_value < centMultBins[8]) {
      pairWithCentMultKtRegistry->fill(HIST("mult_7_8/kstar"), kstar_value);
      if (useKt) {
        auto histMultFolder = HIST("mult_7_8/");
        fillkT(kstar_value, kt_value, histMultFolder);
      }
    } else if (cent_mult_value >= centMultBins[8] && cent_mult_value < centMultBins[9]) {
      pairWithCentMultKtRegistry->fill(HIST("mult_8_9/kstar"), kstar_value);
      if (useKt) {
        auto histMultFolder = HIST("mult_8_9/");
        fillkT(kstar_value, kt_value, histMultFolder);
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
  void fillkT(t1 kstar_value, t1 kt_value, t2 folder)
  {
    if (kt_value >= ktBins[0] && kt_value < ktBins[1]) {
      pairWithCentMultKtRegistry->fill(folder + HIST("kstar_kt_0_1"), kstar_value);
    } else if (kt_value >= ktBins[1] && kt_value < ktBins[2]) {
      pairWithCentMultKtRegistry->fill(folder + HIST("kstar_kt_1_2"), kstar_value);
    } else if (kt_value >= ktBins[2] && kt_value < ktBins[3]) {
      pairWithCentMultKtRegistry->fill(folder + HIST("kstar_kt_2_3"), kstar_value);
    } else if (kt_value >= ktBins[3] && kt_value < ktBins[4]) {
      pairWithCentMultKtRegistry->fill(folder + HIST("kstar_kt_3_4"), kstar_value);
    } else if (kt_value >= ktBins[4] && kt_value < ktBins[5]) {
      pairWithCentMultKtRegistry->fill(folder + HIST("kstar_kt_4_5"), kstar_value);
    } else if (kt_value >= ktBins[5] && kt_value < ktBins[6]) {
      pairWithCentMultKtRegistry->fill(folder + HIST("kstar_kt_5_6"), kstar_value);
    } else if (kt_value >= ktBins[6] && kt_value < ktBins[7]) {
      pairWithCentMultKtRegistry->fill(folder + HIST("kstar_kt_6_7"), kstar_value);
    } else if (kt_value >= ktBins[7] && kt_value < ktBins[8]) {
      pairWithCentMultKtRegistry->fill(folder + HIST("kstar_kt_7_8"), kstar_value);
    } else if (kt_value >= ktBins[8] && kt_value < ktBins[9]) {
      pairWithCentMultKtRegistry->fill(folder + HIST("kstar_kt_8_9"), kstar_value);
    }
  }

  /// @brief
  /// @tparam t1
  /// @param qout_value
  /// @param qside_value
  /// @param qlong_value
  /// @param cent_mult_value
  template <typename t1>
  void fill3D(t1 qout_value, t1 qside_value, t1 qlong_value, t1 cent_mult_value, t1 kt_value)
  {

    if (cent_mult_value >= centMultBins[centMultBins.size() - 1] || cent_mult_value < centMultBins[0]) {
      pairWithCentMultKtRegistry->fill(HIST("Beyond_Max_3D"), qout_value, qside_value, qlong_value);
    } else if (cent_mult_value >= centMultBins[0] && cent_mult_value < centMultBins[1]) {
      pairWithCentMultKtRegistry->fill(HIST("mult_0_1/q3D"), qout_value, qside_value, qlong_value);
      if (use3D) {
        auto histMultFolder = HIST("mult_0_1/");
        fillkT3D(qout_value, qside_value, qlong_value, kt_value, histMultFolder);
      }
    } else if (cent_mult_value >= centMultBins[1] && cent_mult_value < centMultBins[2]) {
      // pairWithCentMultKtRegistry->fill(HIST("mult_1_2/q3D"), qout_value, qside_value, qlong_value);
      if (use3D) {
        auto histMultFolder = HIST("mult_1_2/");
        fillkT3D(qout_value, qside_value, qlong_value, kt_value, histMultFolder);
      }
    } else if (cent_mult_value >= centMultBins[2] && cent_mult_value < centMultBins[3]) {
      // pairWithCentMultKtRegistry->fill(HIST("mult_2_3/q3D"), qout_value, qside_value, qlong_value);
      if (use3D) {
        auto histMultFolder = HIST("mult_2_3/");
        fillkT3D(qout_value, qside_value, qlong_value, kt_value, histMultFolder);
      }
    } else if (cent_mult_value >= centMultBins[3] && cent_mult_value < centMultBins[4]) {
      // pairWithCentMultKtRegistry->fill(HIST("mult_3_4/q3D"), qout_value, qside_value, qlong_value);
      if (use3D) {
        auto histMultFolder = HIST("mult_3_4/");
        fillkT3D(qout_value, qside_value, qlong_value, kt_value, histMultFolder);
      }
    } else if (cent_mult_value >= centMultBins[4] && cent_mult_value < centMultBins[5]) {
      // pairWithCentMultKtRegistry->fill(HIST("mult_4_5/q3D"), qout_value, qside_value, qlong_value);
      if (use3D) {
        auto histMultFolder = HIST("mult_4_5/");
        fillkT3D(qout_value, qside_value, qlong_value, kt_value, histMultFolder);
      }
    } else if (cent_mult_value >= centMultBins[5] && cent_mult_value < centMultBins[6]) {
      // pairWithCentMultKtRegistry->fill(HIST("mult_5_6/q3D"), qout_value, qside_value, qlong_value);
      if (use3D) {
        auto histMultFolder = HIST("mult_5_6/");
        fillkT3D(qout_value, qside_value, qlong_value, kt_value, histMultFolder);
      }
    } else if (cent_mult_value >= centMultBins[6] && cent_mult_value < centMultBins[7]) {
      // pairWithCentMultKtRegistry->fill(HIST("mult_6_7/q3D"), qout_value, qside_value, qlong_value);
      if (use3D) {
        auto histMultFolder = HIST("mult_6_7/");
        fillkT3D(qout_value, qside_value, qlong_value, kt_value, histMultFolder);
      }
    } else if (cent_mult_value >= centMultBins[7] && cent_mult_value < centMultBins[8]) {
      // pairWithCentMultKtRegistry->fill(HIST("mult_7_8/q3D"), qout_value, qside_value, qlong_value);
      if (use3D) {
        auto histMultFolder = HIST("mult_7_8/");
        fillkT3D(qout_value, qside_value, qlong_value, kt_value, histMultFolder);
      }
    } else if (cent_mult_value >= centMultBins[8] && cent_mult_value < centMultBins[9]) {
      // pairWithCentMultKtRegistry->fill(HIST("mult_8_9/q3D"), qout_value, qside_value, qlong_value);
      if (use3D) {
        auto histMultFolder = HIST("mult_8_9/");
        fillkT3D(qout_value, qside_value, qlong_value, kt_value, histMultFolder);
      }
    }
  }

  /// @brief
  /// @tparam t1
  /// @tparam t2
  /// @param qout_value
  /// @param qside_value
  /// @param qlong_value
  /// @param folder
  template <typename t1, typename t2>
  void fillkT3D(t1 qout_value, t1 qside_value, t1 qlong_value, t1 kt_value, t2 folder)
  {
    if (kt_value >= ktBins[0] && kt_value < ktBins[1]) {
      pairWithCentMultKtRegistry->fill(folder + HIST("q3D_kt_0_1"), qout_value, qside_value, qlong_value);
    } else if (kt_value >= ktBins[1] && kt_value < ktBins[2]) {
      pairWithCentMultKtRegistry->fill(folder + HIST("q3D_kt_1_2"), qout_value, qside_value, qlong_value);
    } else if (kt_value >= ktBins[2] && kt_value < ktBins[3]) {
      pairWithCentMultKtRegistry->fill(folder + HIST("q3D_kt_2_3"), qout_value, qside_value, qlong_value);
    } else if (kt_value >= ktBins[3] && kt_value < ktBins[4]) {
      pairWithCentMultKtRegistry->fill(folder + HIST("q3D_kt_3_4"), qout_value, qside_value, qlong_value);
    } else if (kt_value >= ktBins[4] && kt_value < ktBins[5]) {
      pairWithCentMultKtRegistry->fill(folder + HIST("q3D_kt_4_5"), qout_value, qside_value, qlong_value);
    } else if (kt_value >= ktBins[5] && kt_value < ktBins[6]) {
      pairWithCentMultKtRegistry->fill(folder + HIST("q3D_kt_5_6"), qout_value, qside_value, qlong_value);
    }
  }

 protected:
  HistogramRegistry* pairWithCentMultKtRegistry = nullptr;
  std::vector<double> centMultBins;
  std::vector<double> ktBins;
  bool useKt = false;
  bool use3D = false;
  static constexpr std::string_view HistSuffix[10] = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9"};
};
} // namespace o2::analysis::femto_universe

#endif // PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSEPAIRWITHCENTMULTKT_H_
