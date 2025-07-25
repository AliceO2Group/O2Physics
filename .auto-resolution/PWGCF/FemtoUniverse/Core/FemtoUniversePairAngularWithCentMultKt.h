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
/// \file FemtoUniversePairAngularWithCentMultKt.h
/// \brief FemtoUniversePairAngularWithCentMultKt - Histogram class for angular pair tracks with centrality and multiplicity
/// \author Deependra Sharma, IITB, deependra.sharma@cern.ch
/// \author Alicja Płachta, WUT Warsaw, alicja.plachta.stud@pw.edu.pl
/// \author Zuzanna Chochulska, WUT Warsaw & CTU Prague, zchochul@cern.ch

#ifndef PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSEPAIRANGULARWITHCENTMULTKT_H_
#define PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSEPAIRANGULARWITHCENTMULTKT_H_

#include "Framework/HistogramRegistry.h"

#include <string>
#include <vector>

namespace o2::analysis::femto_universe
{
class FemtoUniversePairAngularWithCentMultKt
{
 public:
  virtual ~FemtoUniversePairAngularWithCentMultKt() = default;
  /// @brief
  /// @tparam t1
  /// @param registry
  /// @param kstarbins
  /// @param centmultbins
  template <typename t1, typename t2>
  void init(HistogramRegistry* registry, t1& /*kstarbins*/, t1& centmultbins, t2& phiBins, t2& etaBins, bool processKT)
  {
    pairWithCentMultKtRegistry = registry;
    // AxisSpec kstarAxis = {kstarbins, "#it{k*} (GeV/#it{c})"};

    kPhiLow = (-static_cast<int>(phiBins / 4) + 0.5) * o2::constants::math::TwoPI / phiBins;
    kPhiHigh = o2::constants::math::TwoPI + (-static_cast<int>(phiBins / 4) + 0.5) * o2::constants::math::TwoPI / phiBins;
    framework::AxisSpec phiAxis = {phiBins, kPhiLow, kPhiHigh};
    framework::AxisSpec etaAxis = {etaBins, -2.0, 2.0};

    kCentMultBins = centmultbins;
    ktBins = {0.0f, 0.2f, 0.4f, 0.6f, 0.8f, 1.0f, 2.0f, 99999.f}; // temporary
    ktBins.erase(ktBins.begin());
    kCentMultBins.erase(kCentMultBins.begin());
    useKt = processKT;

    for (int i = 0; i < static_cast<int>(kCentMultBins.size() - 1); i++) {
      int lowBin = static_cast<int>((kCentMultBins[i]));
      int highBin = static_cast<int>((kCentMultBins[i + 1]));
      std::string kHistTitle = "mult_" + std::to_string(lowBin) + "-" + std::to_string(highBin);
      std::string kHistSuffix1 = static_cast<std::string>(HistSuffix[i]);
      std::string kHistSuffix2 = static_cast<std::string>(HistSuffix[i + 1]);
      std::string kHistFolderMult = "mult_" + kHistSuffix1 + "_" + kHistSuffix2;
      std::string kHistName = kHistFolderMult + "/DeltaEtaDeltaPhi";
      pairWithCentMultKtRegistry->add(kHistName.c_str(), kHistTitle.c_str(), HistType::kTH2F, {phiAxis, etaAxis});
      if (useKt) {
        for (int i = 0; i < static_cast<int>(ktBins.size() - 1); i++) {
          std::string ktBin1String = std::to_string(ktBins[i]);
          std::replace(ktBin1String.begin(), ktBin1String.end(), '.', '_');
          std::string ktBin2String = std::to_string(ktBins[i + 1]);
          std::replace(ktBin2String.begin(), ktBin2String.end(), '.', '_');
          ktBin1String.resize(4);
          ktBin2String.resize(4);
          std::string kHistTitleKt = "kt_" + ktBin1String + "-" + ktBin2String;
          std::string kHistSuffix1Kt = static_cast<std::string>(HistSuffix[i]);
          std::string kHistSuffix2Kt = static_cast<std::string>(HistSuffix[i + 1]);
          std::string kHistNameKt = kHistFolderMult + "/DeltaEtaDeltaPhi" + kHistSuffix1Kt + "_" + kHistSuffix2Kt;
          pairWithCentMultKtRegistry->add(kHistNameKt.c_str(), kHistTitleKt.c_str(), HistType::kTH2F, {phiAxis, etaAxis});
        }
      }
    }
    pairWithCentMultKtRegistry->add("Beyond_Max", "Beyond_Max", HistType::kTH2F, {phiAxis, etaAxis});
  }

  /// @brief
  /// @tparam t1
  /// @param kstar_value
  /// @param cent_mult_value
  template <typename t1>
  void fill(t1 kstar_value, t1 cent_mult_value, t1 d_phi_value, t1 d_eta_value)
  {

    if (cent_mult_value > kCentMultBins[kCentMultBins.size() - 1] || cent_mult_value < kCentMultBins[0]) {
      pairWithCentMultKtRegistry->fill(HIST("Beyond_Max"), kstar_value);
    } else if (cent_mult_value <= kCentMultBins[1]) {
      pairWithCentMultKtRegistry->fill(HIST("mult_0_1/DeltaEtaDeltaPhi"), d_phi_value, d_eta_value);
      if (useKt) {
        // auto histMultFolder = HIST("mult_0_1/");
        // fill_kT(kstar_value, kt_value, histMultFolder);
      }
    } else if (cent_mult_value <= kCentMultBins[2]) {
      pairWithCentMultKtRegistry->fill(HIST("mult_1_2/DeltaEtaDeltaPhi"), d_phi_value, d_eta_value);
      if (useKt) {
        // auto histMultFolder = HIST("mult_1_2/");
        // fill_kT(kstar_value, kt_value, histMultFolder);
      }
    } else if (cent_mult_value <= kCentMultBins[3]) {
      pairWithCentMultKtRegistry->fill(HIST("mult_2_3/DeltaEtaDeltaPhi"), d_phi_value, d_eta_value);
      if (useKt) {
        // auto histMultFolder = HIST("mult_2_3/");
        // fill_kT(kstar_value, kt_value, histMultFolder);
      }
    } else if (cent_mult_value <= kCentMultBins[4]) {
      pairWithCentMultKtRegistry->fill(HIST("mult_3_4/DeltaEtaDeltaPhi"), d_phi_value, d_eta_value);
      if (useKt) {
        // auto histMultFolder = HIST("mult_3_4/");
        // fill_kT(kstar_value, kt_value, histMultFolder);
      }
    } else if (cent_mult_value <= kCentMultBins[5]) {
      pairWithCentMultKtRegistry->fill(HIST("mult_4_5/DeltaEtaDeltaPhi"), d_phi_value, d_eta_value);
      if (useKt) {
        // auto histMultFolder = HIST("mult_4_5/");
        // fill_kT(kstar_value, kt_value, histMultFolder);
      }
    } else if (cent_mult_value <= kCentMultBins[6]) {
      pairWithCentMultKtRegistry->fill(HIST("mult_5_6/DeltaEtaDeltaPhi"), d_phi_value, d_eta_value);
      if (useKt) {
        // auto histMultFolder = HIST("mult_5_6/");
        // fill_kT(kstar_value, kt_value, histMultFolder);
      }
    } else if (cent_mult_value <= kCentMultBins[7]) {
      pairWithCentMultKtRegistry->fill(HIST("mult_6_7/DeltaEtaDeltaPhi"), d_phi_value, d_eta_value);
      if (useKt) {
        // auto histMultFolder = HIST("mult_6_7/");
        // fill_kT(kstar_value, kt_value, histMultFolder);
      }
    } else if (cent_mult_value <= kCentMultBins[8]) {
      pairWithCentMultKtRegistry->fill(HIST("mult_7_8/DeltaEtaDeltaPhi"), d_phi_value, d_eta_value);
      if (useKt) {
        // auto histMultFolder = HIST("mult_7_8/");
        // fill_kT(kstar_value, kt_value, histMultFolder);
      }
    } else if (cent_mult_value <= kCentMultBins[9]) {
      pairWithCentMultKtRegistry->fill(HIST("mult_8_9/kstar"), d_phi_value, d_eta_value);
      if (useKt) {
        // auto histMultFolder = HIST("mult_8_9/");
        // fill_kT(kstar_value, kt_value, histMultFolder);
      }
    }
  }

  //   /// @brief
  //   /// @tparam t1
  //   /// @tparam t2
  //   /// @param kstar_value
  //   /// @param kt_value
  //   /// @param folder
  //   template <typename t1, typename t2>
  //   void fill_kT(t1 kstar_value, t1 kt_value, t2 folder)
  //   {
  //     if (kt_value <= ktBins[1]) {
  //       pairWithCentMultKtRegistry->fill(folder + HIST("kstar_kt_0_1"), kstar_value);
  //     } else if (kt_value <= ktBins[2]) {
  //       pairWithCentMultKtRegistry->fill(folder + HIST("kstar_kt_1_2"), kstar_value);
  //     } else if (kt_value <= ktBins[3]) {
  //       pairWithCentMultKtRegistry->fill(folder + HIST("kstar_kt_2_3"), kstar_value);
  //     } else if (kt_value <= ktBins[4]) {
  //       pairWithCentMultKtRegistry->fill(folder + HIST("kstar_kt_3_4"), kstar_value);
  //     } else if (kt_value <= ktBins[5]) {
  //       pairWithCentMultKtRegistry->fill(folder + HIST("kstar_kt_4_5"), kstar_value);
  //     } else if (kt_value <= ktBins[6]) {
  //       pairWithCentMultKtRegistry->fill(folder + HIST("kstar_kt_5_6"), kstar_value);
  //     } else if (kt_value <= ktBins[7]) {
  //       pairWithCentMultKtRegistry->fill(folder + HIST("kstar_kt_6_7"), kstar_value);
  //     } else if (kt_value <= ktBins[8]) {
  //       pairWithCentMultKtRegistry->fill(folder + HIST("kstar_kt_7_8"), kstar_value);
  //     } else if (kt_value <= ktBins[9]) {
  //       pairWithCentMultKtRegistry->fill(folder + HIST("kstar_kt_8_9"), kstar_value);
  //     }
  //   }

 protected:
  HistogramRegistry* pairWithCentMultKtRegistry = nullptr;
  std::vector<double> kCentMultBins;
  std::vector<double> ktBins;
  bool useKt = false;
  double kPhiLow;
  double kPhiHigh;
  double deltaEta;
  double deltaPhi;
  static constexpr std::string_view HistSuffix[10] = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9"};
};
} // namespace o2::analysis::femto_universe

#endif // PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSEPAIRANGULARWITHCENTMULTKT_H_
