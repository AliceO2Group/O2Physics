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

///
/// \file   PIDResponseTPC.h
/// \since  2024-11-15
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \brief  Set of tables, tasks and utilities to provide the interface between
///         the analysis data model and the PID response of the TPC
///

#ifndef COMMON_DATAMODEL_PIDRESPONSETPC_H_
#define COMMON_DATAMODEL_PIDRESPONSETPC_H_

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/Logger.h>
#include <ReconstructionDataFormats/PID.h>

#include <cstdint>
#include <experimental/type_traits>

namespace o2::aod
{
namespace pidutils
{

// Checkers for TPC PID hypothesis availability (runtime)
template <class T>
using hasTPCEl = decltype(std::declval<T&>().tpcNSigmaEl());
template <class T>
using hasTPCMu = decltype(std::declval<T&>().tpcNSigmaMu());
template <class T>
using hasTPCPi = decltype(std::declval<T&>().tpcNSigmaPi());
template <class T>
using hasTPCKa = decltype(std::declval<T&>().tpcNSigmaKa());
template <class T>
using hasTPCPr = decltype(std::declval<T&>().tpcNSigmaPr());
template <class T>
using hasTPCDe = decltype(std::declval<T&>().tpcNSigmaDe());
template <class T>
using hasTPCTr = decltype(std::declval<T&>().tpcNSigmaTr());
template <class T>
using hasTPCHe = decltype(std::declval<T&>().tpcNSigmaHe());
template <class T>
using hasTPCAl = decltype(std::declval<T&>().tpcNSigmaAl());

// PID index as template argument
#define perSpeciesWrapper(functionName)                       \
  template <o2::track::PID::ID index, typename TrackType>     \
  auto functionName(const TrackType& track)                   \
  {                                                           \
    if constexpr (index == o2::track::PID::Electron) {        \
      return track.functionName##El();                        \
    } else if constexpr (index == o2::track::PID::Muon) {     \
      return track.functionName##Mu();                        \
    } else if constexpr (index == o2::track::PID::Pion) {     \
      return track.functionName##Pi();                        \
    } else if constexpr (index == o2::track::PID::Kaon) {     \
      return track.functionName##Ka();                        \
    } else if constexpr (index == o2::track::PID::Proton) {   \
      return track.functionName##Pr();                        \
    } else if constexpr (index == o2::track::PID::Deuteron) { \
      return track.functionName##De();                        \
    } else if constexpr (index == o2::track::PID::Triton) {   \
      return track.functionName##Tr();                        \
    } else if constexpr (index == o2::track::PID::Helium3) {  \
      return track.functionName##He();                        \
    } else if constexpr (index == o2::track::PID::Alpha) {    \
      return track.functionName##Al();                        \
    }                                                         \
  }

perSpeciesWrapper(tpcNSigma);
perSpeciesWrapper(tpcExpSigma);
template <o2::track::PID::ID index, typename TrackType>
auto tpcExpSignal(const TrackType& track)
{
  if constexpr (index == o2::track::PID::Electron) {
    return track.tpcExpSignalEl(track.tpcSignal());
  } else if constexpr (index == o2::track::PID::Muon) {
    return track.tpcExpSignalMu(track.tpcSignal());
  } else if constexpr (index == o2::track::PID::Pion) {
    return track.tpcExpSignalPi(track.tpcSignal());
  } else if constexpr (index == o2::track::PID::Kaon) {
    return track.tpcExpSignalKa(track.tpcSignal());
  } else if constexpr (index == o2::track::PID::Proton) {
    return track.tpcExpSignalPr(track.tpcSignal());
  } else if constexpr (index == o2::track::PID::Deuteron) {
    return track.tpcExpSignalDe(track.tpcSignal());
  } else if constexpr (index == o2::track::PID::Triton) {
    return track.tpcExpSignalTr(track.tpcSignal());
  } else if constexpr (index == o2::track::PID::Helium3) {
    return track.tpcExpSignalHe(track.tpcSignal());
  } else if constexpr (index == o2::track::PID::Alpha) {
    return track.tpcExpSignalAl(track.tpcSignal());
  }
}
perSpeciesWrapper(tpcExpSignalDiff);

#undef perSpeciesWrapper

// PID index as function argument for TPC
#define perSpeciesWrapper(functionName)                                                                             \
  template <typename TrackType>                                                                                     \
  auto functionName(const o2::track::PID::ID index, const TrackType& track)                                         \
  {                                                                                                                 \
    switch (index) {                                                                                                \
      case o2::track::PID::Electron:                                                                                \
        if constexpr (std::experimental::is_detected<hasTPCEl, TrackType>::value) {                                 \
          return track.functionName##El();                                                                          \
        }                                                                                                           \
      case o2::track::PID::Muon:                                                                                    \
        if constexpr (std::experimental::is_detected<hasTPCMu, TrackType>::value) {                                 \
          return track.functionName##Mu();                                                                          \
        }                                                                                                           \
      case o2::track::PID::Pion:                                                                                    \
        if constexpr (std::experimental::is_detected<hasTPCPi, TrackType>::value) {                                 \
          return track.functionName##Pi();                                                                          \
        }                                                                                                           \
      case o2::track::PID::Kaon:                                                                                    \
        if constexpr (std::experimental::is_detected<hasTPCKa, TrackType>::value) {                                 \
          return track.functionName##Ka();                                                                          \
        }                                                                                                           \
      case o2::track::PID::Proton:                                                                                  \
        if constexpr (std::experimental::is_detected<hasTPCPr, TrackType>::value) {                                 \
          return track.functionName##Pr();                                                                          \
        }                                                                                                           \
      case o2::track::PID::Deuteron:                                                                                \
        if constexpr (std::experimental::is_detected<hasTPCDe, TrackType>::value) {                                 \
          return track.functionName##De();                                                                          \
        }                                                                                                           \
      case o2::track::PID::Triton:                                                                                  \
        if constexpr (std::experimental::is_detected<hasTPCTr, TrackType>::value) {                                 \
          return track.functionName##Tr();                                                                          \
        }                                                                                                           \
      case o2::track::PID::Helium3:                                                                                 \
        if constexpr (std::experimental::is_detected<hasTPCHe, TrackType>::value) {                                 \
          return track.functionName##He();                                                                          \
        }                                                                                                           \
      case o2::track::PID::Alpha:                                                                                   \
        if constexpr (std::experimental::is_detected<hasTPCAl, TrackType>::value) {                                 \
          return track.functionName##Al();                                                                          \
        }                                                                                                           \
      default:                                                                                                      \
        LOGF(fatal, "TPC PID table for PID index %i (%s) is not available", index, o2::track::PID::getName(index)); \
        return 0.f;                                                                                                 \
    }                                                                                                               \
  }

perSpeciesWrapper(tpcNSigma);
perSpeciesWrapper(tpcExpSigma);
template <typename TrackType>
auto tpcExpSignal(const o2::track::PID::ID index, const TrackType& track)
{
  switch (index) {
    case o2::track::PID::Electron:
      if constexpr (std::experimental::is_detected<hasTPCEl, TrackType>::value) {
        return track.tpcExpSignalEl(track.tpcSignal());
      }
    case o2::track::PID::Muon:
      if constexpr (std::experimental::is_detected<hasTPCMu, TrackType>::value) {
        return track.tpcExpSignalMu(track.tpcSignal());
      }
    case o2::track::PID::Pion:
      if constexpr (std::experimental::is_detected<hasTPCPi, TrackType>::value) {
        return track.tpcExpSignalPi(track.tpcSignal());
      }
    case o2::track::PID::Kaon:
      if constexpr (std::experimental::is_detected<hasTPCKa, TrackType>::value) {
        return track.tpcExpSignalKa(track.tpcSignal());
      }
    case o2::track::PID::Proton:
      if constexpr (std::experimental::is_detected<hasTPCPr, TrackType>::value) {
        return track.tpcExpSignalPr(track.tpcSignal());
      }
    case o2::track::PID::Deuteron:
      if constexpr (std::experimental::is_detected<hasTPCDe, TrackType>::value) {
        return track.tpcExpSignalDe(track.tpcSignal());
      }
    case o2::track::PID::Triton:
      if constexpr (std::experimental::is_detected<hasTPCTr, TrackType>::value) {
        return track.tpcExpSignalTr(track.tpcSignal());
      }
    case o2::track::PID::Helium3:
      if constexpr (std::experimental::is_detected<hasTPCHe, TrackType>::value) {
        return track.tpcExpSignalHe(track.tpcSignal());
      }
    case o2::track::PID::Alpha:
      if constexpr (std::experimental::is_detected<hasTPCAl, TrackType>::value) {
        return track.tpcExpSignalAl(track.tpcSignal());
      }
    default:
      LOGF(fatal, "TPC PID table for PID index %i (%s) is not available", index, o2::track::PID::getName(index));
      return 0.f;
  }
}
perSpeciesWrapper(tpcExpSignalDiff);

#undef perSpeciesWrapper

} // namespace pidutils

namespace pidtpc
{
// Expected signals
DECLARE_SOA_DYNAMIC_COLUMN(TPCExpSignalEl, tpcExpSignalEl, //! Expected signal with the TPC detector for electron
                           [](float nsigma, float sigma, float tpcsignal) -> float { return tpcsignal - nsigma * sigma; });
DECLARE_SOA_DYNAMIC_COLUMN(TPCExpSignalMu, tpcExpSignalMu, //! Expected signal with the TPC detector for muon
                           [](float nsigma, float sigma, float tpcsignal) -> float { return tpcsignal - nsigma * sigma; });
DECLARE_SOA_DYNAMIC_COLUMN(TPCExpSignalPi, tpcExpSignalPi, //! Expected signal with the TPC detector for pion
                           [](float nsigma, float sigma, float tpcsignal) -> float { return tpcsignal - nsigma * sigma; });
DECLARE_SOA_DYNAMIC_COLUMN(TPCExpSignalKa, tpcExpSignalKa, //! Expected signal with the TPC detector for kaon
                           [](float nsigma, float sigma, float tpcsignal) -> float { return tpcsignal - nsigma * sigma; });
DECLARE_SOA_DYNAMIC_COLUMN(TPCExpSignalPr, tpcExpSignalPr, //! Expected signal with the TPC detector for proton
                           [](float nsigma, float sigma, float tpcsignal) -> float { return tpcsignal - nsigma * sigma; });
DECLARE_SOA_DYNAMIC_COLUMN(TPCExpSignalDe, tpcExpSignalDe, //! Expected signal with the TPC detector for deuteron
                           [](float nsigma, float sigma, float tpcsignal) -> float { return tpcsignal - nsigma * sigma; });
DECLARE_SOA_DYNAMIC_COLUMN(TPCExpSignalTr, tpcExpSignalTr, //! Expected signal with the TPC detector for triton
                           [](float nsigma, float sigma, float tpcsignal) -> float { return tpcsignal - nsigma * sigma; });
DECLARE_SOA_DYNAMIC_COLUMN(TPCExpSignalHe, tpcExpSignalHe, //! Expected signal with the TPC detector for helium3
                           [](float nsigma, float sigma, float tpcsignal) -> float { return tpcsignal - nsigma * sigma; });
DECLARE_SOA_DYNAMIC_COLUMN(TPCExpSignalAl, tpcExpSignalAl, //! Expected signal with the TPC detector for alpha
                           [](float nsigma, float sigma, float tpcsignal) -> float { return tpcsignal - nsigma * sigma; });
// Delta with respect to signal
DECLARE_SOA_DYNAMIC_COLUMN(TPCExpSignalDiffEl, tpcExpSignalDiffEl, //! Difference between signal and expected for electron
                           [](float nsigma, float sigma) -> float { return nsigma * sigma; });
DECLARE_SOA_DYNAMIC_COLUMN(TPCExpSignalDiffMu, tpcExpSignalDiffMu, //! Difference between signal and expected for muon
                           [](float nsigma, float sigma) -> float { return nsigma * sigma; });
DECLARE_SOA_DYNAMIC_COLUMN(TPCExpSignalDiffPi, tpcExpSignalDiffPi, //! Difference between signal and expected for pion
                           [](float nsigma, float sigma) -> float { return nsigma * sigma; });
DECLARE_SOA_DYNAMIC_COLUMN(TPCExpSignalDiffKa, tpcExpSignalDiffKa, //! Difference between signal and expected for kaon
                           [](float nsigma, float sigma) -> float { return nsigma * sigma; });
DECLARE_SOA_DYNAMIC_COLUMN(TPCExpSignalDiffPr, tpcExpSignalDiffPr, //! Difference between signal and expected for proton
                           [](float nsigma, float sigma) -> float { return nsigma * sigma; });
DECLARE_SOA_DYNAMIC_COLUMN(TPCExpSignalDiffDe, tpcExpSignalDiffDe, //! Difference between signal and expected for deuteron
                           [](float nsigma, float sigma) -> float { return nsigma * sigma; });
DECLARE_SOA_DYNAMIC_COLUMN(TPCExpSignalDiffTr, tpcExpSignalDiffTr, //! Difference between signal and expected for triton
                           [](float nsigma, float sigma) -> float { return nsigma * sigma; });
DECLARE_SOA_DYNAMIC_COLUMN(TPCExpSignalDiffHe, tpcExpSignalDiffHe, //! Difference between signal and expected for helium3
                           [](float nsigma, float sigma) -> float { return nsigma * sigma; });
DECLARE_SOA_DYNAMIC_COLUMN(TPCExpSignalDiffAl, tpcExpSignalDiffAl, //! Difference between signal and expected for alpha
                           [](float nsigma, float sigma) -> float { return nsigma * sigma; });
// Expected sigma
DECLARE_SOA_COLUMN(TPCExpSigmaEl, tpcExpSigmaEl, float); //! Expected resolution with the TPC detector for electron
DECLARE_SOA_COLUMN(TPCExpSigmaMu, tpcExpSigmaMu, float); //! Expected resolution with the TPC detector for muon
DECLARE_SOA_COLUMN(TPCExpSigmaPi, tpcExpSigmaPi, float); //! Expected resolution with the TPC detector for pion
DECLARE_SOA_COLUMN(TPCExpSigmaKa, tpcExpSigmaKa, float); //! Expected resolution with the TPC detector for kaon
DECLARE_SOA_COLUMN(TPCExpSigmaPr, tpcExpSigmaPr, float); //! Expected resolution with the TPC detector for proton
DECLARE_SOA_COLUMN(TPCExpSigmaDe, tpcExpSigmaDe, float); //! Expected resolution with the TPC detector for deuteron
DECLARE_SOA_COLUMN(TPCExpSigmaTr, tpcExpSigmaTr, float); //! Expected resolution with the TPC detector for triton
DECLARE_SOA_COLUMN(TPCExpSigmaHe, tpcExpSigmaHe, float); //! Expected resolution with the TPC detector for helium3
DECLARE_SOA_COLUMN(TPCExpSigmaAl, tpcExpSigmaAl, float); //! Expected resolution with the TPC detector for alpha
// NSigma
DECLARE_SOA_COLUMN(TPCNSigmaEl, tpcNSigmaEl, float); //! Nsigma separation with the TPC detector for electron
DECLARE_SOA_COLUMN(TPCNSigmaMu, tpcNSigmaMu, float); //! Nsigma separation with the TPC detector for muon
DECLARE_SOA_COLUMN(TPCNSigmaPi, tpcNSigmaPi, float); //! Nsigma separation with the TPC detector for pion
DECLARE_SOA_COLUMN(TPCNSigmaKa, tpcNSigmaKa, float); //! Nsigma separation with the TPC detector for kaon
DECLARE_SOA_COLUMN(TPCNSigmaPr, tpcNSigmaPr, float); //! Nsigma separation with the TPC detector for proton
DECLARE_SOA_COLUMN(TPCNSigmaDe, tpcNSigmaDe, float); //! Nsigma separation with the TPC detector for deuteron
DECLARE_SOA_COLUMN(TPCNSigmaTr, tpcNSigmaTr, float); //! Nsigma separation with the TPC detector for triton
DECLARE_SOA_COLUMN(TPCNSigmaHe, tpcNSigmaHe, float); //! Nsigma separation with the TPC detector for helium3
DECLARE_SOA_COLUMN(TPCNSigmaAl, tpcNSigmaAl, float); //! Nsigma separation with the TPC detector for alpha

} // namespace pidtpc

namespace pidtpc_tiny
{
struct binning {
 public:
  typedef int8_t binned_t;
  static constexpr int nbins = (1 << 8 * sizeof(binned_t)) - 2;
  static constexpr binned_t overflowBin = nbins >> 1;
  static constexpr binned_t underflowBin = -(nbins >> 1);
  static constexpr float binned_max = 6.35;
  static constexpr float binned_min = -6.35;
  static constexpr float bin_width = (binned_max - binned_min) / nbins;

  // Function to pack a float into a binned value in table
  template <typename T>
  static void packInTable(const float& valueToBin, T& table)
  {
    if (valueToBin <= binned_min) {
      table(underflowBin);
    } else if (valueToBin >= binned_max) {
      table(overflowBin);
    } else if (valueToBin >= 0) {
      table(static_cast<binned_t>((valueToBin / bin_width) + 0.5f));
    } else {
      table(static_cast<binned_t>((valueToBin / bin_width) - 0.5f));
    }
  }

  // Function to unpack a binned value into a float
  static float unPackInTable(const binned_t& valueToUnpack)
  {
    return bin_width * static_cast<float>(valueToUnpack);
  }
};

// NSigma with reduced size 8 bit
DECLARE_SOA_COLUMN(TPCNSigmaStoreEl, tpcNSigmaStoreEl, binning::binned_t); //! Stored binned nsigma with the TPC detector for electron
DECLARE_SOA_COLUMN(TPCNSigmaStoreMu, tpcNSigmaStoreMu, binning::binned_t); //! Stored binned nsigma with the TPC detector for muon
DECLARE_SOA_COLUMN(TPCNSigmaStorePi, tpcNSigmaStorePi, binning::binned_t); //! Stored binned nsigma with the TPC detector for pion
DECLARE_SOA_COLUMN(TPCNSigmaStoreKa, tpcNSigmaStoreKa, binning::binned_t); //! Stored binned nsigma with the TPC detector for kaon
DECLARE_SOA_COLUMN(TPCNSigmaStorePr, tpcNSigmaStorePr, binning::binned_t); //! Stored binned nsigma with the TPC detector for proton
DECLARE_SOA_COLUMN(TPCNSigmaStoreDe, tpcNSigmaStoreDe, binning::binned_t); //! Stored binned nsigma with the TPC detector for deuteron
DECLARE_SOA_COLUMN(TPCNSigmaStoreTr, tpcNSigmaStoreTr, binning::binned_t); //! Stored binned nsigma with the TPC detector for triton
DECLARE_SOA_COLUMN(TPCNSigmaStoreHe, tpcNSigmaStoreHe, binning::binned_t); //! Stored binned nsigma with the TPC detector for helium3
DECLARE_SOA_COLUMN(TPCNSigmaStoreAl, tpcNSigmaStoreAl, binning::binned_t); //! Stored binned nsigma with the TPC detector for alpha

// NSigma with reduced size in [binned_min, binned_max] bin size bin_width
DECLARE_SOA_DYNAMIC_COLUMN(TPCNSigmaEl, tpcNSigmaEl, //! Unwrapped (float) nsigma with the TPC detector for electron
                           [](binning::binned_t nsigma_binned) -> float { return binning::unPackInTable(nsigma_binned); });
DECLARE_SOA_DYNAMIC_COLUMN(TPCNSigmaMu, tpcNSigmaMu, //! Unwrapped (float) nsigma with the TPC detector for muon
                           [](binning::binned_t nsigma_binned) -> float { return binning::unPackInTable(nsigma_binned); });
DECLARE_SOA_DYNAMIC_COLUMN(TPCNSigmaPi, tpcNSigmaPi, //! Unwrapped (float) nsigma with the TPC detector for pion
                           [](binning::binned_t nsigma_binned) -> float { return binning::unPackInTable(nsigma_binned); });
DECLARE_SOA_DYNAMIC_COLUMN(TPCNSigmaKa, tpcNSigmaKa, //! Unwrapped (float) nsigma with the TPC detector for kaon
                           [](binning::binned_t nsigma_binned) -> float { return binning::unPackInTable(nsigma_binned); });
DECLARE_SOA_DYNAMIC_COLUMN(TPCNSigmaPr, tpcNSigmaPr, //! Unwrapped (float) nsigma with the TPC detector for proton
                           [](binning::binned_t nsigma_binned) -> float { return binning::unPackInTable(nsigma_binned); });
DECLARE_SOA_DYNAMIC_COLUMN(TPCNSigmaDe, tpcNSigmaDe, //! Unwrapped (float) nsigma with the TPC detector for deuteron
                           [](binning::binned_t nsigma_binned) -> float { return binning::unPackInTable(nsigma_binned); });
DECLARE_SOA_DYNAMIC_COLUMN(TPCNSigmaTr, tpcNSigmaTr, //! Unwrapped (float) nsigma with the TPC detector for triton
                           [](binning::binned_t nsigma_binned) -> float { return binning::unPackInTable(nsigma_binned); });
DECLARE_SOA_DYNAMIC_COLUMN(TPCNSigmaHe, tpcNSigmaHe, //! Unwrapped (float) nsigma with the TPC detector for helium3
                           [](binning::binned_t nsigma_binned) -> float { return binning::unPackInTable(nsigma_binned); });
DECLARE_SOA_DYNAMIC_COLUMN(TPCNSigmaAl, tpcNSigmaAl, //! Unwrapped (float) nsigma with the TPC detector for alpha
                           [](binning::binned_t nsigma_binned) -> float { return binning::unPackInTable(nsigma_binned); });

} // namespace pidtpc_tiny

// Per particle tables
DECLARE_SOA_TABLE(pidTPCFullEl, "AOD", "pidTPCFullEl", //! Table of the TPC (full) response with expected signal, expected resolution and Nsigma for electron
                  pidtpc::TPCExpSignalEl<pidtpc::TPCNSigmaEl, pidtpc::TPCExpSigmaEl>,
                  pidtpc::TPCExpSignalDiffEl<pidtpc::TPCNSigmaEl, pidtpc::TPCExpSigmaEl>,
                  pidtpc::TPCExpSigmaEl, pidtpc::TPCNSigmaEl);
DECLARE_SOA_TABLE(pidTPCFullMu, "AOD", "pidTPCFullMu", //! Table of the TPC (full) response with expected signal, expected resolution and Nsigma for muon
                  pidtpc::TPCExpSignalMu<pidtpc::TPCNSigmaMu, pidtpc::TPCExpSigmaMu>,
                  pidtpc::TPCExpSignalDiffMu<pidtpc::TPCNSigmaMu, pidtpc::TPCExpSigmaMu>,
                  pidtpc::TPCExpSigmaMu, pidtpc::TPCNSigmaMu);
DECLARE_SOA_TABLE(pidTPCFullPi, "AOD", "pidTPCFullPi", //! Table of the TPC (full) response with expected signal, expected resolution and Nsigma for pion
                  pidtpc::TPCExpSignalPi<pidtpc::TPCNSigmaPi, pidtpc::TPCExpSigmaPi>,
                  pidtpc::TPCExpSignalDiffPi<pidtpc::TPCNSigmaPi, pidtpc::TPCExpSigmaPi>,
                  pidtpc::TPCExpSigmaPi, pidtpc::TPCNSigmaPi);
DECLARE_SOA_TABLE(pidTPCFullKa, "AOD", "pidTPCFullKa", //! Table of the TPC (full) response with expected signal, expected resolution and Nsigma for kaon
                  pidtpc::TPCExpSignalKa<pidtpc::TPCNSigmaKa, pidtpc::TPCExpSigmaKa>,
                  pidtpc::TPCExpSignalDiffKa<pidtpc::TPCNSigmaKa, pidtpc::TPCExpSigmaKa>,
                  pidtpc::TPCExpSigmaKa, pidtpc::TPCNSigmaKa);
DECLARE_SOA_TABLE(pidTPCFullPr, "AOD", "pidTPCFullPr", //! Table of the TPC (full) response with expected signal, expected resolution and Nsigma for proton
                  pidtpc::TPCExpSignalPr<pidtpc::TPCNSigmaPr, pidtpc::TPCExpSigmaPr>,
                  pidtpc::TPCExpSignalDiffPr<pidtpc::TPCNSigmaPr, pidtpc::TPCExpSigmaPr>,
                  pidtpc::TPCExpSigmaPr, pidtpc::TPCNSigmaPr);
DECLARE_SOA_TABLE(pidTPCFullDe, "AOD", "pidTPCFullDe", //! Table of the TPC (full) response with expected signal, expected resolution and Nsigma for deuteron
                  pidtpc::TPCExpSignalDe<pidtpc::TPCNSigmaDe, pidtpc::TPCExpSigmaDe>,
                  pidtpc::TPCExpSignalDiffDe<pidtpc::TPCNSigmaDe, pidtpc::TPCExpSigmaDe>,
                  pidtpc::TPCExpSigmaDe, pidtpc::TPCNSigmaDe);
DECLARE_SOA_TABLE(pidTPCFullTr, "AOD", "pidTPCFullTr", //! Table of the TPC (full) response with expected signal, expected resolution and Nsigma for triton
                  pidtpc::TPCExpSignalTr<pidtpc::TPCNSigmaTr, pidtpc::TPCExpSigmaTr>,
                  pidtpc::TPCExpSignalDiffTr<pidtpc::TPCNSigmaTr, pidtpc::TPCExpSigmaTr>,
                  pidtpc::TPCExpSigmaTr, pidtpc::TPCNSigmaTr);
DECLARE_SOA_TABLE(pidTPCFullHe, "AOD", "pidTPCFullHe", //! Table of the TPC (full) response with expected signal, expected resolution and Nsigma for helium3
                  pidtpc::TPCExpSignalHe<pidtpc::TPCNSigmaHe, pidtpc::TPCExpSigmaHe>,
                  pidtpc::TPCExpSignalDiffHe<pidtpc::TPCNSigmaHe, pidtpc::TPCExpSigmaHe>,
                  pidtpc::TPCExpSigmaHe, pidtpc::TPCNSigmaHe);
DECLARE_SOA_TABLE(pidTPCFullAl, "AOD", "pidTPCFullAl", //! Table of the TPC (full) response with expected signal, expected resolution and Nsigma for alpha
                  pidtpc::TPCExpSignalAl<pidtpc::TPCNSigmaAl, pidtpc::TPCExpSigmaAl>,
                  pidtpc::TPCExpSignalDiffAl<pidtpc::TPCNSigmaAl, pidtpc::TPCExpSigmaAl>,
                  pidtpc::TPCExpSigmaAl, pidtpc::TPCNSigmaAl);

// Tiny size tables
DECLARE_SOA_TABLE(pidTPCEl, "AOD", "pidTPCEl", //! Table of the TPC response with binned Nsigma for electron
                  pidtpc_tiny::TPCNSigmaStoreEl, pidtpc_tiny::TPCNSigmaEl<pidtpc_tiny::TPCNSigmaStoreEl>);
DECLARE_SOA_TABLE(pidTPCMu, "AOD", "pidTPCMu", //! Table of the TPC response with binned Nsigma for muon
                  pidtpc_tiny::TPCNSigmaStoreMu, pidtpc_tiny::TPCNSigmaMu<pidtpc_tiny::TPCNSigmaStoreMu>);
DECLARE_SOA_TABLE(pidTPCPi, "AOD", "pidTPCPi", //! Table of the TPC response with binned Nsigma for pion
                  pidtpc_tiny::TPCNSigmaStorePi, pidtpc_tiny::TPCNSigmaPi<pidtpc_tiny::TPCNSigmaStorePi>);
DECLARE_SOA_TABLE(pidTPCKa, "AOD", "pidTPCKa", //! Table of the TPC response with binned Nsigma for kaon
                  pidtpc_tiny::TPCNSigmaStoreKa, pidtpc_tiny::TPCNSigmaKa<pidtpc_tiny::TPCNSigmaStoreKa>);
DECLARE_SOA_TABLE(pidTPCPr, "AOD", "pidTPCPr", //! Table of the TPC response with binned Nsigma for proton
                  pidtpc_tiny::TPCNSigmaStorePr, pidtpc_tiny::TPCNSigmaPr<pidtpc_tiny::TPCNSigmaStorePr>);
DECLARE_SOA_TABLE(pidTPCDe, "AOD", "pidTPCDe", //! Table of the TPC response with binned Nsigma for deuteron
                  pidtpc_tiny::TPCNSigmaStoreDe, pidtpc_tiny::TPCNSigmaDe<pidtpc_tiny::TPCNSigmaStoreDe>);
DECLARE_SOA_TABLE(pidTPCTr, "AOD", "pidTPCTr", //! Table of the TPC response with binned Nsigma for triton
                  pidtpc_tiny::TPCNSigmaStoreTr, pidtpc_tiny::TPCNSigmaTr<pidtpc_tiny::TPCNSigmaStoreTr>);
DECLARE_SOA_TABLE(pidTPCHe, "AOD", "pidTPCHe", //! Table of the TPC response with binned Nsigma for helium3
                  pidtpc_tiny::TPCNSigmaStoreHe, pidtpc_tiny::TPCNSigmaHe<pidtpc_tiny::TPCNSigmaStoreHe>);
DECLARE_SOA_TABLE(pidTPCAl, "AOD", "pidTPCAl", //! Table of the TPC response with binned Nsigma for alpha
                  pidtpc_tiny::TPCNSigmaStoreAl, pidtpc_tiny::TPCNSigmaAl<pidtpc_tiny::TPCNSigmaStoreAl>);

// Extra tables
namespace mcpidtpc
{
// Tuned MC on data
DECLARE_SOA_COLUMN(DeDxTunedMc, mcTunedTPCSignal, float); //! TPC signal after TuneOnData application for MC
} // namespace mcpidtpc

DECLARE_SOA_TABLE(mcTPCTuneOnData, "AOD", "MCTPCTUNEONDATA", mcpidtpc::DeDxTunedMc);

} // namespace o2::aod

#endif // COMMON_DATAMODEL_PIDRESPONSETPC_H_
