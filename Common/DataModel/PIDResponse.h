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
/// \file   PIDResponse.h
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \brief  Set of tables, tasks and utilities to provide the interface between
///         the analysis data model and the PID response
///

#ifndef COMMON_DATAMODEL_PIDRESPONSE_H_
#define COMMON_DATAMODEL_PIDRESPONSE_H_

#include <experimental/type_traits>

// O2 includes
#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "ReconstructionDataFormats/PID.h"
#include "Framework/Logger.h"

namespace o2::aod
{
namespace pidutils
{
// Function to pack a float into a binned value in table
template <typename binningType, typename T>
void packInTable(const float& valueToBin, T& table)
{
  if (valueToBin <= binningType::binned_min) {
    table(binningType::underflowBin);
  } else if (valueToBin >= binningType::binned_max) {
    table(binningType::overflowBin);
  } else if (valueToBin >= 0) {
    table(static_cast<typename binningType::binned_t>((valueToBin / binningType::bin_width) + 0.5f));
  } else {
    table(static_cast<typename binningType::binned_t>((valueToBin / binningType::bin_width) - 0.5f));
  }
}

// Function to unpack a binned value into a float
template <typename binningType>
float unPackInTable(const typename binningType::binned_t& valueToUnpack)
{
  return binningType::bin_width * static_cast<float>(valueToUnpack);
}

// Checkers for TOF PID hypothesis availability (runtime)
template <class T>
using hasTOFEl = decltype(std::declval<T&>().tofNSigmaEl());
template <class T>
using hasTOFMu = decltype(std::declval<T&>().tofNSigmaMu());
template <class T>
using hasTOFPi = decltype(std::declval<T&>().tofNSigmaPi());
template <class T>
using hasTOFKa = decltype(std::declval<T&>().tofNSigmaKa());
template <class T>
using hasTOFPr = decltype(std::declval<T&>().tofNSigmaPr());
template <class T>
using hasTOFDe = decltype(std::declval<T&>().tofNSigmaDe());
template <class T>
using hasTOFTr = decltype(std::declval<T&>().tofNSigmaTr());
template <class T>
using hasTOFHe = decltype(std::declval<T&>().tofNSigmaHe());
template <class T>
using hasTOFAl = decltype(std::declval<T&>().tofNSigmaAl());

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
#define PER_SPECIES_WRAPPER(functionName)                     \
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

PER_SPECIES_WRAPPER(tofNSigma);
PER_SPECIES_WRAPPER(tofExpSigma);
template <o2::track::PID::ID index, typename TrackType>
auto tofExpSignal(const TrackType& track)
{
  if constexpr (index == o2::track::PID::Electron) {
    return track.tofExpSignalEl(track.tofSignal());
  } else if constexpr (index == o2::track::PID::Muon) {
    return track.tofExpSignalMu(track.tofSignal());
  } else if constexpr (index == o2::track::PID::Pion) {
    return track.tofExpSignalPi(track.tofSignal());
  } else if constexpr (index == o2::track::PID::Kaon) {
    return track.tofExpSignalKa(track.tofSignal());
  } else if constexpr (index == o2::track::PID::Proton) {
    return track.tofExpSignalPr(track.tofSignal());
  } else if constexpr (index == o2::track::PID::Deuteron) {
    return track.tofExpSignalDe(track.tofSignal());
  } else if constexpr (index == o2::track::PID::Triton) {
    return track.tofExpSignalTr(track.tofSignal());
  } else if constexpr (index == o2::track::PID::Helium3) {
    return track.tofExpSignalHe(track.tofSignal());
  } else if constexpr (index == o2::track::PID::Alpha) {
    return track.tofExpSignalAl(track.tofSignal());
  }
}
PER_SPECIES_WRAPPER(tofExpSignalDiff);

PER_SPECIES_WRAPPER(tpcNSigma);
PER_SPECIES_WRAPPER(tpcExpSigma);
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
PER_SPECIES_WRAPPER(tpcExpSignalDiff);

#undef PER_SPECIES_WRAPPER

// PID index as function argument for TOF
#define PER_SPECIES_WRAPPER(functionName)                                                                           \
  template <typename TrackType>                                                                                     \
  auto functionName(const o2::track::PID::ID index, const TrackType& track)                                         \
  {                                                                                                                 \
    switch (index) {                                                                                                \
      case o2::track::PID::Electron:                                                                                \
        if constexpr (std::experimental::is_detected<hasTOFEl, TrackType>::value) {                                 \
          return track.functionName##El();                                                                          \
        }                                                                                                           \
      case o2::track::PID::Muon:                                                                                    \
        if constexpr (std::experimental::is_detected<hasTOFMu, TrackType>::value) {                                 \
          return track.functionName##Mu();                                                                          \
        }                                                                                                           \
      case o2::track::PID::Pion:                                                                                    \
        if constexpr (std::experimental::is_detected<hasTOFPi, TrackType>::value) {                                 \
          return track.functionName##Pi();                                                                          \
        }                                                                                                           \
      case o2::track::PID::Kaon:                                                                                    \
        if constexpr (std::experimental::is_detected<hasTOFKa, TrackType>::value) {                                 \
          return track.functionName##Ka();                                                                          \
        }                                                                                                           \
      case o2::track::PID::Proton:                                                                                  \
        if constexpr (std::experimental::is_detected<hasTOFPr, TrackType>::value) {                                 \
          return track.functionName##Pr();                                                                          \
        }                                                                                                           \
      case o2::track::PID::Deuteron:                                                                                \
        if constexpr (std::experimental::is_detected<hasTOFDe, TrackType>::value) {                                 \
          return track.functionName##De();                                                                          \
        }                                                                                                           \
      case o2::track::PID::Triton:                                                                                  \
        if constexpr (std::experimental::is_detected<hasTOFTr, TrackType>::value) {                                 \
          return track.functionName##Tr();                                                                          \
        }                                                                                                           \
      case o2::track::PID::Helium3:                                                                                 \
        if constexpr (std::experimental::is_detected<hasTOFHe, TrackType>::value) {                                 \
          return track.functionName##He();                                                                          \
        }                                                                                                           \
      case o2::track::PID::Alpha:                                                                                   \
        if constexpr (std::experimental::is_detected<hasTOFAl, TrackType>::value) {                                 \
          return track.functionName##Al();                                                                          \
        }                                                                                                           \
      default:                                                                                                      \
        LOGF(fatal, "TOF PID table for PID index %i (%s) is not available", index, o2::track::PID::getName(index)); \
        return 0.f;                                                                                                 \
    }                                                                                                               \
  }

PER_SPECIES_WRAPPER(tofNSigma);
PER_SPECIES_WRAPPER(tofExpSigma);
template <typename TrackType>
auto tofExpSignal(const o2::track::PID::ID index, const TrackType& track)
{
  switch (index) {
    case o2::track::PID::Electron:
      if constexpr (std::experimental::is_detected<hasTOFEl, TrackType>::value) {
        return track.tofExpSignalEl(track.tofSignal());
      }
    case o2::track::PID::Muon:
      if constexpr (std::experimental::is_detected<hasTOFMu, TrackType>::value) {
        return track.tofExpSignalMu(track.tofSignal());
      }
    case o2::track::PID::Pion:
      if constexpr (std::experimental::is_detected<hasTOFPi, TrackType>::value) {
        return track.tofExpSignalPi(track.tofSignal());
      }
    case o2::track::PID::Kaon:
      if constexpr (std::experimental::is_detected<hasTOFKa, TrackType>::value) {
        return track.tofExpSignalKa(track.tofSignal());
      }
    case o2::track::PID::Proton:
      if constexpr (std::experimental::is_detected<hasTOFPr, TrackType>::value) {
        return track.tofExpSignalPr(track.tofSignal());
      }
    case o2::track::PID::Deuteron:
      if constexpr (std::experimental::is_detected<hasTOFDe, TrackType>::value) {
        return track.tofExpSignalDe(track.tofSignal());
      }
    case o2::track::PID::Triton:
      if constexpr (std::experimental::is_detected<hasTOFTr, TrackType>::value) {
        return track.tofExpSignalTr(track.tofSignal());
      }
    case o2::track::PID::Helium3:
      if constexpr (std::experimental::is_detected<hasTOFHe, TrackType>::value) {
        return track.tofExpSignalHe(track.tofSignal());
      }
    case o2::track::PID::Alpha:
      if constexpr (std::experimental::is_detected<hasTOFAl, TrackType>::value) {
        return track.tofExpSignalAl(track.tofSignal());
      }
    default:
      LOGF(fatal, "TOF PID table for PID index %i (%s) is not available", index, o2::track::PID::getName(index));
      return 0.f;
  }
}
PER_SPECIES_WRAPPER(tofExpSignalDiff);

#undef PER_SPECIES_WRAPPER

// PID index as function argument for TPC
#define PER_SPECIES_WRAPPER(functionName)                                                                           \
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

PER_SPECIES_WRAPPER(tpcNSigma);
PER_SPECIES_WRAPPER(tpcExpSigma);
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
PER_SPECIES_WRAPPER(tpcExpSignalDiff);

#undef PER_SPECIES_WRAPPER

} // namespace pidutils

namespace pidflags
{

namespace enums
{
enum PIDFlags : uint8_t {
  EvTimeUndef = 0x0,  // Event collision not set, corresponding to the LHC Fill event time
  EvTimeTOF = 0x1,    // Event collision time from TOF
  EvTimeT0AC = 0x2,   // Event collision time from the FT0AC
  EvTimeTOFT0AC = 0x4 // Event collision time from the TOF and FT0AC
};
}

DECLARE_SOA_COLUMN(GoodTOFMatch, goodTOFMatch, bool);        //! Bool for the TOF PID information on the single track information
DECLARE_SOA_COLUMN(TOFFlags, tofFlags, uint8_t);             //! Flag for the complementary TOF PID information for the event time
DECLARE_SOA_DYNAMIC_COLUMN(IsEvTimeDefined, isEvTimeDefined, //! True if the Event Time was computed with any method i.e. there is a usable event time
                           [](uint8_t flags) -> bool { return (flags > 0); });
DECLARE_SOA_DYNAMIC_COLUMN(IsEvTimeTOF, isEvTimeTOF, //! True if the Event Time was computed with the TOF
                           [](uint8_t flags) -> bool { return (flags & enums::PIDFlags::EvTimeTOF) == enums::PIDFlags::EvTimeTOF; });
DECLARE_SOA_DYNAMIC_COLUMN(IsEvTimeT0AC, isEvTimeT0AC, //! True if the Event Time was computed with the T0AC
                           [](uint8_t flags) -> bool { return (flags & enums::PIDFlags::EvTimeT0AC) == enums::PIDFlags::EvTimeT0AC; });
DECLARE_SOA_DYNAMIC_COLUMN(IsEvTimeTOFT0AC, isEvTimeTOFT0AC, //! True if the Event Time was computed with the TOF and T0AC
                           [](uint8_t flags) -> bool { return (flags & enums::PIDFlags::EvTimeTOFT0AC) == enums::PIDFlags::EvTimeTOFT0AC; });

} // namespace pidflags

namespace pidtofsignal
{
DECLARE_SOA_COLUMN(TOFSignal, tofSignal, float);                   //! TOF signal from track time
DECLARE_SOA_DYNAMIC_COLUMN(EventCollisionTime, eventCollisionTime, //! Event collision time used for the track. Needs the TOF
                           [](float signal, float tMinusTexp, float texp) -> float { return texp + tMinusTexp - signal; });

} // namespace pidtofsignal

namespace pidtofevtime
{
DECLARE_SOA_COLUMN(TOFEvTime, tofEvTime, float);       //! event time for TOF signal. Can be obtained via a combination of detectors e.g. TOF, FT0A, FT0C
DECLARE_SOA_COLUMN(TOFEvTimeErr, tofEvTimeErr, float); //! event time error for TOF. Can be obtained via a combination of detectors e.g. TOF, FT0A, FT0C
} // namespace pidtofevtime

namespace pidtofbeta
{
DECLARE_SOA_COLUMN(Beta, beta, float);           //! TOF beta
DECLARE_SOA_COLUMN(BetaError, betaerror, float); //! Uncertainty on the TOF beta
//
DECLARE_SOA_COLUMN(ExpBetaEl, expbetael, float);           //! Expected beta of electron
DECLARE_SOA_COLUMN(ExpBetaElError, expbetaelerror, float); //! Expected uncertainty on the beta of electron
//
DECLARE_SOA_COLUMN(SeparationBetaEl, separationbetael, float); //! Separation computed with the expected beta for electrons
DECLARE_SOA_DYNAMIC_COLUMN(DiffBetaEl, diffbetael,             //! Difference between the measured and the expected beta for electrons
                           [](float beta, float expbetael) -> float { return beta - expbetael; });
} // namespace pidtofbeta

namespace pidtofmass
{
DECLARE_SOA_COLUMN(TOFMass, mass, float); //! TOF mass
} // namespace pidtofmass

namespace pidtof
{
// Expected times
DECLARE_SOA_DYNAMIC_COLUMN(TOFExpSignalEl, tofExpSignalEl, //! Expected time for electron
                           [](float nsigma, float sigma, float tofsignal) -> float { return tofsignal - nsigma * sigma; });
DECLARE_SOA_DYNAMIC_COLUMN(TOFExpSignalMu, tofExpSignalMu, //! Expected time for muon
                           [](float nsigma, float sigma, float tofsignal) -> float { return tofsignal - nsigma * sigma; });
DECLARE_SOA_DYNAMIC_COLUMN(TOFExpSignalPi, tofExpSignalPi, //! Expected time for pion
                           [](float nsigma, float sigma, float tofsignal) -> float { return tofsignal - nsigma * sigma; });
DECLARE_SOA_DYNAMIC_COLUMN(TOFExpSignalKa, tofExpSignalKa, //! Expected time for kaon
                           [](float nsigma, float sigma, float tofsignal) -> float { return tofsignal - nsigma * sigma; });
DECLARE_SOA_DYNAMIC_COLUMN(TOFExpSignalPr, tofExpSignalPr, //! Expected time for proton
                           [](float nsigma, float sigma, float tofsignal) -> float { return tofsignal - nsigma * sigma; });
DECLARE_SOA_DYNAMIC_COLUMN(TOFExpSignalDe, tofExpSignalDe, //! Expected time for deuteron
                           [](float nsigma, float sigma, float tofsignal) -> float { return tofsignal - nsigma * sigma; });
DECLARE_SOA_DYNAMIC_COLUMN(TOFExpSignalTr, tofExpSignalTr, //! Expected time for triton
                           [](float nsigma, float sigma, float tofsignal) -> float { return tofsignal - nsigma * sigma; });
DECLARE_SOA_DYNAMIC_COLUMN(TOFExpSignalHe, tofExpSignalHe, //! Expected time for helium3
                           [](float nsigma, float sigma, float tofsignal) -> float { return tofsignal - nsigma * sigma; });
DECLARE_SOA_DYNAMIC_COLUMN(TOFExpSignalAl, tofExpSignalAl, //! Expected time for alpha
                           [](float nsigma, float sigma, float tofsignal) -> float { return tofsignal - nsigma * sigma; });
// Delta with respect to signal
DECLARE_SOA_DYNAMIC_COLUMN(TOFExpSignalDiffEl, tofExpSignalDiffEl, //! Difference between signal and expected for electron
                           [](float nsigma, float sigma) -> float { return nsigma * sigma; });
DECLARE_SOA_DYNAMIC_COLUMN(TOFExpSignalDiffMu, tofExpSignalDiffMu, //! Difference between signal and expected for muon
                           [](float nsigma, float sigma) -> float { return nsigma * sigma; });
DECLARE_SOA_DYNAMIC_COLUMN(TOFExpSignalDiffPi, tofExpSignalDiffPi, //! Difference between signal and expected for pion
                           [](float nsigma, float sigma) -> float { return nsigma * sigma; });
DECLARE_SOA_DYNAMIC_COLUMN(TOFExpSignalDiffKa, tofExpSignalDiffKa, //! Difference between signal and expected for kaon
                           [](float nsigma, float sigma) -> float { return nsigma * sigma; });
DECLARE_SOA_DYNAMIC_COLUMN(TOFExpSignalDiffPr, tofExpSignalDiffPr, //! Difference between signal and expected for proton
                           [](float nsigma, float sigma) -> float { return nsigma * sigma; });
DECLARE_SOA_DYNAMIC_COLUMN(TOFExpSignalDiffDe, tofExpSignalDiffDe, //! Difference between signal and expected for deuteron
                           [](float nsigma, float sigma) -> float { return nsigma * sigma; });
DECLARE_SOA_DYNAMIC_COLUMN(TOFExpSignalDiffTr, tofExpSignalDiffTr, //! Difference between signal and expected for triton
                           [](float nsigma, float sigma) -> float { return nsigma * sigma; });
DECLARE_SOA_DYNAMIC_COLUMN(TOFExpSignalDiffHe, tofExpSignalDiffHe, //! Difference between signal and expected for helium3
                           [](float nsigma, float sigma) -> float { return nsigma * sigma; });
DECLARE_SOA_DYNAMIC_COLUMN(TOFExpSignalDiffAl, tofExpSignalDiffAl, //! Difference between signal and expected for alpha
                           [](float nsigma, float sigma) -> float { return nsigma * sigma; });
// Expected sigma
DECLARE_SOA_COLUMN(TOFExpSigmaEl, tofExpSigmaEl, float); //! Expected resolution with the TOF detector for electron
DECLARE_SOA_COLUMN(TOFExpSigmaMu, tofExpSigmaMu, float); //! Expected resolution with the TOF detector for muon
DECLARE_SOA_COLUMN(TOFExpSigmaPi, tofExpSigmaPi, float); //! Expected resolution with the TOF detector for pion
DECLARE_SOA_COLUMN(TOFExpSigmaKa, tofExpSigmaKa, float); //! Expected resolution with the TOF detector for kaon
DECLARE_SOA_COLUMN(TOFExpSigmaPr, tofExpSigmaPr, float); //! Expected resolution with the TOF detector for proton
DECLARE_SOA_COLUMN(TOFExpSigmaDe, tofExpSigmaDe, float); //! Expected resolution with the TOF detector for deuteron
DECLARE_SOA_COLUMN(TOFExpSigmaTr, tofExpSigmaTr, float); //! Expected resolution with the TOF detector for triton
DECLARE_SOA_COLUMN(TOFExpSigmaHe, tofExpSigmaHe, float); //! Expected resolution with the TOF detector for helium3
DECLARE_SOA_COLUMN(TOFExpSigmaAl, tofExpSigmaAl, float); //! Expected resolution with the TOF detector for alpha
// NSigma
DECLARE_SOA_COLUMN(TOFNSigmaEl, tofNSigmaEl, float); //! Nsigma separation with the TOF detector for electron
DECLARE_SOA_COLUMN(TOFNSigmaMu, tofNSigmaMu, float); //! Nsigma separation with the TOF detector for muon
DECLARE_SOA_COLUMN(TOFNSigmaPi, tofNSigmaPi, float); //! Nsigma separation with the TOF detector for pion
DECLARE_SOA_COLUMN(TOFNSigmaKa, tofNSigmaKa, float); //! Nsigma separation with the TOF detector for kaon
DECLARE_SOA_COLUMN(TOFNSigmaPr, tofNSigmaPr, float); //! Nsigma separation with the TOF detector for proton
DECLARE_SOA_COLUMN(TOFNSigmaDe, tofNSigmaDe, float); //! Nsigma separation with the TOF detector for deuteron
DECLARE_SOA_COLUMN(TOFNSigmaTr, tofNSigmaTr, float); //! Nsigma separation with the TOF detector for triton
DECLARE_SOA_COLUMN(TOFNSigmaHe, tofNSigmaHe, float); //! Nsigma separation with the TOF detector for helium3
DECLARE_SOA_COLUMN(TOFNSigmaAl, tofNSigmaAl, float); //! Nsigma separation with the TOF detector for alpha
} // namespace pidtof

// Macro to convert the stored Nsigmas to floats
#define DEFINE_UNWRAP_NSIGMA_COLUMN(COLUMN, COLUMN_NAME) \
  DECLARE_SOA_DYNAMIC_COLUMN(COLUMN, COLUMN_NAME,        \
                             [](binning::binned_t nsigma_binned) -> float { return o2::aod::pidutils::unPackInTable<binning>(nsigma_binned); });

namespace pidtof_tiny
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
};

// NSigma with reduced size 8 bit
DECLARE_SOA_COLUMN(TOFNSigmaStoreEl, tofNSigmaStoreEl, binning::binned_t); //! Stored binned nsigma with the TOF detector for electron
DECLARE_SOA_COLUMN(TOFNSigmaStoreMu, tofNSigmaStoreMu, binning::binned_t); //! Stored binned nsigma with the TOF detector for muon
DECLARE_SOA_COLUMN(TOFNSigmaStorePi, tofNSigmaStorePi, binning::binned_t); //! Stored binned nsigma with the TOF detector for pion
DECLARE_SOA_COLUMN(TOFNSigmaStoreKa, tofNSigmaStoreKa, binning::binned_t); //! Stored binned nsigma with the TOF detector for kaon
DECLARE_SOA_COLUMN(TOFNSigmaStorePr, tofNSigmaStorePr, binning::binned_t); //! Stored binned nsigma with the TOF detector for proton
DECLARE_SOA_COLUMN(TOFNSigmaStoreDe, tofNSigmaStoreDe, binning::binned_t); //! Stored binned nsigma with the TOF detector for deuteron
DECLARE_SOA_COLUMN(TOFNSigmaStoreTr, tofNSigmaStoreTr, binning::binned_t); //! Stored binned nsigma with the TOF detector for triton
DECLARE_SOA_COLUMN(TOFNSigmaStoreHe, tofNSigmaStoreHe, binning::binned_t); //! Stored binned nsigma with the TOF detector for helium3
DECLARE_SOA_COLUMN(TOFNSigmaStoreAl, tofNSigmaStoreAl, binning::binned_t); //! Stored binned nsigma with the TOF detector for alpha
// NSigma with reduced size in [binned_min, binned_max] bin size bin_width
DEFINE_UNWRAP_NSIGMA_COLUMN(TOFNSigmaEl, tofNSigmaEl); //! Unwrapped (float) nsigma with the TOF detector for electron
DEFINE_UNWRAP_NSIGMA_COLUMN(TOFNSigmaMu, tofNSigmaMu); //! Unwrapped (float) nsigma with the TOF detector for muon
DEFINE_UNWRAP_NSIGMA_COLUMN(TOFNSigmaPi, tofNSigmaPi); //! Unwrapped (float) nsigma with the TOF detector for pion
DEFINE_UNWRAP_NSIGMA_COLUMN(TOFNSigmaKa, tofNSigmaKa); //! Unwrapped (float) nsigma with the TOF detector for kaon
DEFINE_UNWRAP_NSIGMA_COLUMN(TOFNSigmaPr, tofNSigmaPr); //! Unwrapped (float) nsigma with the TOF detector for proton
DEFINE_UNWRAP_NSIGMA_COLUMN(TOFNSigmaDe, tofNSigmaDe); //! Unwrapped (float) nsigma with the TOF detector for deuteron
DEFINE_UNWRAP_NSIGMA_COLUMN(TOFNSigmaTr, tofNSigmaTr); //! Unwrapped (float) nsigma with the TOF detector for triton
DEFINE_UNWRAP_NSIGMA_COLUMN(TOFNSigmaHe, tofNSigmaHe); //! Unwrapped (float) nsigma with the TOF detector for helium3
DEFINE_UNWRAP_NSIGMA_COLUMN(TOFNSigmaAl, tofNSigmaAl); //! Unwrapped (float) nsigma with the TOF detector for alpha

} // namespace pidtof_tiny

DECLARE_SOA_TABLE(TOFSignal, "AOD", "TOFSignal", //! Table of the TOF signal
                  pidtofsignal::TOFSignal,
                  pidtofsignal::EventCollisionTime<pidtofsignal::TOFSignal>);

DECLARE_SOA_TABLE(TOFEvTime, "AOD", "TOFEvTime", //! Table of the TOF event time. One entry per track.
                  pidtofevtime::TOFEvTime,
                  pidtofevtime::TOFEvTimeErr);

DECLARE_SOA_TABLE(pidTOFFlags, "AOD", "pidTOFFlags", //! Table of the flags for TOF signal quality on the track level
                  pidflags::GoodTOFMatch);

DECLARE_SOA_TABLE(pidTOFbeta, "AOD", "pidTOFbeta", //! Table of the TOF beta
                  pidtofbeta::Beta, pidtofbeta::BetaError);

DECLARE_SOA_TABLE(pidTOFmass, "AOD", "pidTOFmass", //! Table of the TOF mass
                  pidtofmass::TOFMass);

DECLARE_SOA_TABLE(pidEvTimeFlags, "AOD", "pidEvTimeFlags", //! Table of the PID flags for the event time tables
                  pidflags::TOFFlags,
                  pidflags::IsEvTimeDefined<pidflags::TOFFlags>,
                  pidflags::IsEvTimeTOF<pidflags::TOFFlags>,
                  pidflags::IsEvTimeT0AC<pidflags::TOFFlags>,
                  pidflags::IsEvTimeTOFT0AC<pidflags::TOFFlags>);

// Per particle tables
DECLARE_SOA_TABLE(pidTOFFullEl, "AOD", "pidTOFFullEl", //! Table of the TOF (full) response with expected signal, expected resolution and Nsigma for electron
                  pidtof::TOFExpSignalDiffEl<pidtof::TOFNSigmaEl, pidtof::TOFExpSigmaEl>,
                  pidtof::TOFExpSignalEl<pidtof::TOFNSigmaEl, pidtof::TOFExpSigmaEl>,
                  pidtof::TOFExpSigmaEl, pidtof::TOFNSigmaEl);
DECLARE_SOA_TABLE(pidTOFFullMu, "AOD", "pidTOFFullMu", //! Table of the TOF (full) response with expected signal, expected resolution and Nsigma for muon
                  pidtof::TOFExpSignalDiffMu<pidtof::TOFNSigmaMu, pidtof::TOFExpSigmaMu>,
                  pidtof::TOFExpSignalMu<pidtof::TOFNSigmaMu, pidtof::TOFExpSigmaMu>,
                  pidtof::TOFExpSigmaMu, pidtof::TOFNSigmaMu);
DECLARE_SOA_TABLE(pidTOFFullPi, "AOD", "pidTOFFullPi", //! Table of the TOF (full) response with expected signal, expected resolution and Nsigma for pion
                  pidtof::TOFExpSignalDiffPi<pidtof::TOFNSigmaPi, pidtof::TOFExpSigmaPi>,
                  pidtof::TOFExpSignalPi<pidtof::TOFNSigmaPi, pidtof::TOFExpSigmaPi>,
                  pidtof::TOFExpSigmaPi, pidtof::TOFNSigmaPi);
DECLARE_SOA_TABLE(pidTOFFullKa, "AOD", "pidTOFFullKa", //! Table of the TOF (full) response with expected signal, expected resolution and Nsigma for kaon
                  pidtof::TOFExpSignalDiffKa<pidtof::TOFNSigmaKa, pidtof::TOFExpSigmaKa>,
                  pidtof::TOFExpSignalKa<pidtof::TOFNSigmaKa, pidtof::TOFExpSigmaKa>,
                  pidtof::TOFExpSigmaKa, pidtof::TOFNSigmaKa);
DECLARE_SOA_TABLE(pidTOFFullPr, "AOD", "pidTOFFullPr", //! Table of the TOF (full) response with expected signal, expected resolution and Nsigma for proton
                  pidtof::TOFExpSignalDiffPr<pidtof::TOFNSigmaPr, pidtof::TOFExpSigmaPr>,
                  pidtof::TOFExpSignalPr<pidtof::TOFNSigmaPr, pidtof::TOFExpSigmaPr>,
                  pidtof::TOFExpSigmaPr, pidtof::TOFNSigmaPr);
DECLARE_SOA_TABLE(pidTOFFullDe, "AOD", "pidTOFFullDe", //! Table of the TOF (full) response with expected signal, expected resolution and Nsigma for deuteron
                  pidtof::TOFExpSignalDiffDe<pidtof::TOFNSigmaDe, pidtof::TOFExpSigmaDe>,
                  pidtof::TOFExpSignalDe<pidtof::TOFNSigmaDe, pidtof::TOFExpSigmaDe>,
                  pidtof::TOFExpSigmaDe, pidtof::TOFNSigmaDe);
DECLARE_SOA_TABLE(pidTOFFullTr, "AOD", "pidTOFFullTr", //! Table of the TOF (full) response with expected signal, expected resolution and Nsigma for triton
                  pidtof::TOFExpSignalDiffTr<pidtof::TOFNSigmaTr, pidtof::TOFExpSigmaTr>,
                  pidtof::TOFExpSignalTr<pidtof::TOFNSigmaTr, pidtof::TOFExpSigmaTr>,
                  pidtof::TOFExpSigmaTr, pidtof::TOFNSigmaTr);
DECLARE_SOA_TABLE(pidTOFFullHe, "AOD", "pidTOFFullHe", //! Table of the TOF (full) response with expected signal, expected resolution and Nsigma for helium3
                  pidtof::TOFExpSignalDiffHe<pidtof::TOFNSigmaHe, pidtof::TOFExpSigmaHe>,
                  pidtof::TOFExpSignalHe<pidtof::TOFNSigmaHe, pidtof::TOFExpSigmaHe>,
                  pidtof::TOFExpSigmaHe, pidtof::TOFNSigmaHe);
DECLARE_SOA_TABLE(pidTOFFullAl, "AOD", "pidTOFFullAl", //! Table of the TOF (full) response with expected signal, expected resolution and Nsigma for alpha
                  pidtof::TOFExpSignalDiffAl<pidtof::TOFNSigmaAl, pidtof::TOFExpSigmaAl>,
                  pidtof::TOFExpSignalAl<pidtof::TOFNSigmaAl, pidtof::TOFExpSigmaAl>,
                  pidtof::TOFExpSigmaAl, pidtof::TOFNSigmaAl);

// Tiny size tables
DECLARE_SOA_TABLE(pidTOFEl, "AOD", "pidTOFEl", //! Table of the TOF response with binned Nsigma for electron
                  pidtof_tiny::TOFNSigmaStoreEl, pidtof_tiny::TOFNSigmaEl<pidtof_tiny::TOFNSigmaStoreEl>);
DECLARE_SOA_TABLE(pidTOFMu, "AOD", "pidTOFMu", //! Table of the TOF response with binned Nsigma for muon
                  pidtof_tiny::TOFNSigmaStoreMu, pidtof_tiny::TOFNSigmaMu<pidtof_tiny::TOFNSigmaStoreMu>);
DECLARE_SOA_TABLE(pidTOFPi, "AOD", "pidTOFPi", //! Table of the TOF response with binned Nsigma for pion
                  pidtof_tiny::TOFNSigmaStorePi, pidtof_tiny::TOFNSigmaPi<pidtof_tiny::TOFNSigmaStorePi>);
DECLARE_SOA_TABLE(pidTOFKa, "AOD", "pidTOFKa", //! Table of the TOF response with binned Nsigma for kaon
                  pidtof_tiny::TOFNSigmaStoreKa, pidtof_tiny::TOFNSigmaKa<pidtof_tiny::TOFNSigmaStoreKa>);
DECLARE_SOA_TABLE(pidTOFPr, "AOD", "pidTOFPr", //! Table of the TOF response with binned Nsigma for proton
                  pidtof_tiny::TOFNSigmaStorePr, pidtof_tiny::TOFNSigmaPr<pidtof_tiny::TOFNSigmaStorePr>);
DECLARE_SOA_TABLE(pidTOFDe, "AOD", "pidTOFDe", //! Table of the TOF response with binned Nsigma for deuteron
                  pidtof_tiny::TOFNSigmaStoreDe, pidtof_tiny::TOFNSigmaDe<pidtof_tiny::TOFNSigmaStoreDe>);
DECLARE_SOA_TABLE(pidTOFTr, "AOD", "pidTOFTr", //! Table of the TOF response with binned Nsigma for triton
                  pidtof_tiny::TOFNSigmaStoreTr, pidtof_tiny::TOFNSigmaTr<pidtof_tiny::TOFNSigmaStoreTr>);
DECLARE_SOA_TABLE(pidTOFHe, "AOD", "pidTOFHe", //! Table of the TOF response with binned Nsigma for helium3
                  pidtof_tiny::TOFNSigmaStoreHe, pidtof_tiny::TOFNSigmaHe<pidtof_tiny::TOFNSigmaStoreHe>);
DECLARE_SOA_TABLE(pidTOFAl, "AOD", "pidTOFAl", //! Table of the TOF response with binned Nsigma for alpha
                  pidtof_tiny::TOFNSigmaStoreAl, pidtof_tiny::TOFNSigmaAl<pidtof_tiny::TOFNSigmaStoreAl>);

namespace mcpidtpc
{
// Tuned MC on data
DECLARE_SOA_COLUMN(DeDxTunedMc, mcTunedTPCSignal, float); //! TPC signal after TuneOnData application for MC
} // namespace mcpidtpc

DECLARE_SOA_TABLE(mcTPCTuneOnData, "AOD", "MCTPCTUNEONDATA", mcpidtpc::DeDxTunedMc);

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
// Expected signals difference
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
};

// NSigma with reduced size
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
DEFINE_UNWRAP_NSIGMA_COLUMN(TPCNSigmaEl, tpcNSigmaEl); //! Unwrapped (float) nsigma with the TPC detector for electron
DEFINE_UNWRAP_NSIGMA_COLUMN(TPCNSigmaMu, tpcNSigmaMu); //! Unwrapped (float) nsigma with the TPC detector for muon
DEFINE_UNWRAP_NSIGMA_COLUMN(TPCNSigmaPi, tpcNSigmaPi); //! Unwrapped (float) nsigma with the TPC detector for pion
DEFINE_UNWRAP_NSIGMA_COLUMN(TPCNSigmaKa, tpcNSigmaKa); //! Unwrapped (float) nsigma with the TPC detector for kaon
DEFINE_UNWRAP_NSIGMA_COLUMN(TPCNSigmaPr, tpcNSigmaPr); //! Unwrapped (float) nsigma with the TPC detector for proton
DEFINE_UNWRAP_NSIGMA_COLUMN(TPCNSigmaDe, tpcNSigmaDe); //! Unwrapped (float) nsigma with the TPC detector for deuteron
DEFINE_UNWRAP_NSIGMA_COLUMN(TPCNSigmaTr, tpcNSigmaTr); //! Unwrapped (float) nsigma with the TPC detector for triton
DEFINE_UNWRAP_NSIGMA_COLUMN(TPCNSigmaHe, tpcNSigmaHe); //! Unwrapped (float) nsigma with the TPC detector for helium3
DEFINE_UNWRAP_NSIGMA_COLUMN(TPCNSigmaAl, tpcNSigmaAl); //! Unwrapped (float) nsigma with the TPC detector for alpha

} // namespace pidtpc_tiny

// Per particle tables
DECLARE_SOA_TABLE(pidTPCFullEl, "AOD", "pidTPCFullEl", //! Table of the TPC (full) response with expected signal, expected resolution and Nsigma for electron
                  pidtpc::TPCExpSignalEl<pidtpc::TPCNSigmaEl, pidtpc::TPCExpSigmaEl>, pidtpc::TPCExpSignalDiffEl<pidtpc::TPCNSigmaEl, pidtpc::TPCExpSigmaEl>, pidtpc::TPCExpSigmaEl, pidtpc::TPCNSigmaEl);
DECLARE_SOA_TABLE(pidTPCFullMu, "AOD", "pidTPCFullMu", //! Table of the TPC (full) response with expected signal, expected resolution and Nsigma for muon
                  pidtpc::TPCExpSignalMu<pidtpc::TPCNSigmaMu, pidtpc::TPCExpSigmaMu>, pidtpc::TPCExpSignalDiffMu<pidtpc::TPCNSigmaMu, pidtpc::TPCExpSigmaMu>, pidtpc::TPCExpSigmaMu, pidtpc::TPCNSigmaMu);
DECLARE_SOA_TABLE(pidTPCFullPi, "AOD", "pidTPCFullPi", //! Table of the TPC (full) response with expected signal, expected resolution and Nsigma for pion
                  pidtpc::TPCExpSignalPi<pidtpc::TPCNSigmaPi, pidtpc::TPCExpSigmaPi>, pidtpc::TPCExpSignalDiffPi<pidtpc::TPCNSigmaPi, pidtpc::TPCExpSigmaPi>, pidtpc::TPCExpSigmaPi, pidtpc::TPCNSigmaPi);
DECLARE_SOA_TABLE(pidTPCFullKa, "AOD", "pidTPCFullKa", //! Table of the TPC (full) response with expected signal, expected resolution and Nsigma for kaon
                  pidtpc::TPCExpSignalKa<pidtpc::TPCNSigmaKa, pidtpc::TPCExpSigmaKa>, pidtpc::TPCExpSignalDiffKa<pidtpc::TPCNSigmaKa, pidtpc::TPCExpSigmaKa>, pidtpc::TPCExpSigmaKa, pidtpc::TPCNSigmaKa);
DECLARE_SOA_TABLE(pidTPCFullPr, "AOD", "pidTPCFullPr", //! Table of the TPC (full) response with expected signal, expected resolution and Nsigma for proton
                  pidtpc::TPCExpSignalPr<pidtpc::TPCNSigmaPr, pidtpc::TPCExpSigmaPr>, pidtpc::TPCExpSignalDiffPr<pidtpc::TPCNSigmaPr, pidtpc::TPCExpSigmaPr>, pidtpc::TPCExpSigmaPr, pidtpc::TPCNSigmaPr);
DECLARE_SOA_TABLE(pidTPCFullDe, "AOD", "pidTPCFullDe", //! Table of the TPC (full) response with expected signal, expected resolution and Nsigma for deuteron
                  pidtpc::TPCExpSignalDe<pidtpc::TPCNSigmaDe, pidtpc::TPCExpSigmaDe>, pidtpc::TPCExpSignalDiffDe<pidtpc::TPCNSigmaDe, pidtpc::TPCExpSigmaDe>, pidtpc::TPCExpSigmaDe, pidtpc::TPCNSigmaDe);
DECLARE_SOA_TABLE(pidTPCFullTr, "AOD", "pidTPCFullTr", //! Table of the TPC (full) response with expected signal, expected resolution and Nsigma for triton
                  pidtpc::TPCExpSignalTr<pidtpc::TPCNSigmaTr, pidtpc::TPCExpSigmaTr>, pidtpc::TPCExpSignalDiffTr<pidtpc::TPCNSigmaTr, pidtpc::TPCExpSigmaTr>, pidtpc::TPCExpSigmaTr, pidtpc::TPCNSigmaTr);
DECLARE_SOA_TABLE(pidTPCFullHe, "AOD", "pidTPCFullHe", //! Table of the TPC (full) response with expected signal, expected resolution and Nsigma for helium3
                  pidtpc::TPCExpSignalHe<pidtpc::TPCNSigmaHe, pidtpc::TPCExpSigmaHe>, pidtpc::TPCExpSignalDiffHe<pidtpc::TPCNSigmaHe, pidtpc::TPCExpSigmaHe>, pidtpc::TPCExpSigmaHe, pidtpc::TPCNSigmaHe);
DECLARE_SOA_TABLE(pidTPCFullAl, "AOD", "pidTPCFullAl", //! Table of the TPC (full) response with expected signal, expected resolution and Nsigma for alpha
                  pidtpc::TPCExpSignalAl<pidtpc::TPCNSigmaAl, pidtpc::TPCExpSigmaAl>, pidtpc::TPCExpSignalDiffAl<pidtpc::TPCNSigmaAl, pidtpc::TPCExpSigmaAl>, pidtpc::TPCExpSigmaAl, pidtpc::TPCNSigmaAl);

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

#undef DEFINE_UNWRAP_NSIGMA_COLUMN

namespace pidbayes
{
typedef int8_t binned_prob_t;
// Bayesian probabilities with reduced size
DECLARE_SOA_COLUMN(BayesEl, bayesEl, binned_prob_t);                //! Bayesian probability for electron expressed in %
DECLARE_SOA_COLUMN(BayesMu, bayesMu, binned_prob_t);                //! Bayesian probability for muon expressed in %
DECLARE_SOA_COLUMN(BayesPi, bayesPi, binned_prob_t);                //! Bayesian probability for pion expressed in %
DECLARE_SOA_COLUMN(BayesKa, bayesKa, binned_prob_t);                //! Bayesian probability for kaon expressed in %
DECLARE_SOA_COLUMN(BayesPr, bayesPr, binned_prob_t);                //! Bayesian probability for proton expressed in %
DECLARE_SOA_COLUMN(BayesDe, bayesDe, binned_prob_t);                //! Bayesian probability for deuteron expressed in %
DECLARE_SOA_COLUMN(BayesTr, bayesTr, binned_prob_t);                //! Bayesian probability for triton expressed in %
DECLARE_SOA_COLUMN(BayesHe, bayesHe, binned_prob_t);                //! Bayesian probability for helium3 expressed in %
DECLARE_SOA_COLUMN(BayesAl, bayesAl, binned_prob_t);                //! Bayesian probability for alpha expressed in %
DECLARE_SOA_COLUMN(BayesProb, bayesProb, binned_prob_t);            //! Bayesian probability of the most probable ID
DECLARE_SOA_COLUMN(BayesID, bayesID, o2::track::pid_constants::ID); //! Most probable ID

} // namespace pidbayes

// Table for each particle hypothesis
DECLARE_SOA_TABLE(pidBayesEl, "AOD", "pidBayesEl", //! Binned (in percentage) Bayesian probability of having a Electron
                  pidbayes::BayesEl);
DECLARE_SOA_TABLE(pidBayesMu, "AOD", "pidBayesMu", //! Binned (in percentage) Bayesian probability of having a Muon
                  pidbayes::BayesMu);
DECLARE_SOA_TABLE(pidBayesPi, "AOD", "pidBayesPi", //! Binned (in percentage) Bayesian probability of having a Pion
                  pidbayes::BayesPi);
DECLARE_SOA_TABLE(pidBayesKa, "AOD", "pidBayesKa", //! Binned (in percentage) Bayesian probability of having a Kaon
                  pidbayes::BayesKa);
DECLARE_SOA_TABLE(pidBayesPr, "AOD", "pidBayesPr", //! Binned (in percentage) Bayesian probability of having a Proton
                  pidbayes::BayesPr);
DECLARE_SOA_TABLE(pidBayesDe, "AOD", "pidBayesDe", //! Binned (in percentage) Bayesian probability of having a Deuteron
                  pidbayes::BayesDe);
DECLARE_SOA_TABLE(pidBayesTr, "AOD", "pidBayesTr", //! Binned (in percentage) Bayesian probability of having a Triton
                  pidbayes::BayesTr);
DECLARE_SOA_TABLE(pidBayesHe, "AOD", "pidBayesHe", //! Binned (in percentage) Bayesian probability of having a Helium3
                  pidbayes::BayesHe);
DECLARE_SOA_TABLE(pidBayesAl, "AOD", "pidBayesAl", //! Binned (in percentage) Bayesian probability of having a Alpha
                  pidbayes::BayesAl);

// Table for the most probable particle
DECLARE_SOA_TABLE(pidBayes, "AOD", "pidBayes", pidbayes::BayesProb, pidbayes::BayesID); //! Index of the most probable ID and its bayesian probability

} // namespace o2::aod

#endif // COMMON_DATAMODEL_PIDRESPONSE_H_
