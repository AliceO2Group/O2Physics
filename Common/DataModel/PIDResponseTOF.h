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
/// \file   PIDResponseTOF.h
/// \since  2024-11-15
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \brief  Set of tables, tasks and utilities to provide the interface between
///         the analysis data model and the PID response of the TOF
///

#ifndef COMMON_DATAMODEL_PIDRESPONSETOF_H_
#define COMMON_DATAMODEL_PIDRESPONSETOF_H_

#include "Common/Core/PID/PIDTOF.h"
#include "Common/Core/PID/PIDTOFParamService.h"

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

perSpeciesWrapper(tofNSigma);
perSpeciesWrapper(tofExpSigma);
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
perSpeciesWrapper(tofExpSignalDiff);

#undef perSpeciesWrapper

// PID index as function argument for TOF
#define perSpeciesWrapper(functionName)                                                                             \
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

perSpeciesWrapper(tofNSigma);
perSpeciesWrapper(tofExpSigma);
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
perSpeciesWrapper(tofExpSignalDiff);

#undef perSpeciesWrapper

} // namespace pidutils

// Extra tables
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

DECLARE_SOA_TABLE(pidTOFFlags, "AOD", "pidTOFFlags", //! Table of the flags for TOF signal quality on the track level
                  pidflags::GoodTOFMatch);

DECLARE_SOA_TABLE(pidEvTimeFlags, "AOD", "pidEvTimeFlags", //! Table of the PID flags for the event time tables
                  pidflags::TOFFlags,
                  pidflags::IsEvTimeDefined<pidflags::TOFFlags>,
                  pidflags::IsEvTimeTOF<pidflags::TOFFlags>,
                  pidflags::IsEvTimeT0AC<pidflags::TOFFlags>,
                  pidflags::IsEvTimeTOFT0AC<pidflags::TOFFlags>);

namespace pidtofsignal
{
DECLARE_SOA_COLUMN(TOFSignal, tofSignal, float);                   //! TOF signal from track time
DECLARE_SOA_DYNAMIC_COLUMN(EventCollisionTime, eventCollisionTime, //! Event collision time used for the track. Needs the TOF
                           [](float signal, float tMinusTexp, float texp) -> float { return texp + tMinusTexp - signal; });

} // namespace pidtofsignal

DECLARE_SOA_TABLE(TOFSignal, "AOD", "TOFSignal", //! Table of the TOF signal
                  pidtofsignal::TOFSignal,
                  pidtofsignal::EventCollisionTime<pidtofsignal::TOFSignal>);

namespace pidtofevtime
{
DECLARE_SOA_COLUMN(TOFEvTime, tofEvTime, float);       //! event time for TOF signal. Can be obtained via a combination of detectors e.g. TOF, FT0A, FT0C
DECLARE_SOA_COLUMN(TOFEvTimeErr, tofEvTimeErr, float); //! event time error for TOF. Can be obtained via a combination of detectors e.g. TOF, FT0A, FT0C
} // namespace pidtofevtime

DECLARE_SOA_TABLE(TOFEvTime, "AOD", "TOFEvTime", //! Table of the TOF event time. One entry per track.
                  pidtofevtime::TOFEvTime,
                  pidtofevtime::TOFEvTimeErr);

namespace pidtof
{
// Expected signals
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

//! Expected resolution with the TOF detector for electron (computed on the fly)
#define PERSPECIES_TOF_SIGMA_COLUMN(name, id)                                                         \
  DECLARE_SOA_DYNAMIC_COLUMN(TOFExpSigma##name##Imp, tofExpSigmaDyn##name,                            \
                             [](float tofSignal,                                                      \
                                float tofExpMom,                                                      \
                                float momentum,                                                       \
                                float eta,                                                            \
                                float tofEvTimeErr) -> float {                                        \
                               return o2::pid::tof::TOFResponseImpl::expectedSigma<id>(tofSignal,     \
                                                                                       tofExpMom,     \
                                                                                       momentum,      \
                                                                                       eta,           \
                                                                                       tofEvTimeErr); \
                             });

PERSPECIES_TOF_SIGMA_COLUMN(El, o2::track::PID::Electron);
PERSPECIES_TOF_SIGMA_COLUMN(Mu, o2::track::PID::Muon);
PERSPECIES_TOF_SIGMA_COLUMN(Pi, o2::track::PID::Pion);
PERSPECIES_TOF_SIGMA_COLUMN(Ka, o2::track::PID::Kaon);
PERSPECIES_TOF_SIGMA_COLUMN(Pr, o2::track::PID::Proton);
PERSPECIES_TOF_SIGMA_COLUMN(De, o2::track::PID::Deuteron);
PERSPECIES_TOF_SIGMA_COLUMN(Tr, o2::track::PID::Triton);
PERSPECIES_TOF_SIGMA_COLUMN(He, o2::track::PID::Helium3);
PERSPECIES_TOF_SIGMA_COLUMN(Al, o2::track::PID::Alpha);
#undef PERSPECIES_TOF_SIGMA_COLUMN

#define PERSPECIES_TOF_SEPARATION_COLUMN(name, id)                                             \
  DECLARE_SOA_DYNAMIC_COLUMN(TOFNSigma##name##Imp, tofNSigmaDyn##name,                         \
                             [](const float tofSignal,                                         \
                                const float tofExpMom,                                         \
                                const float length,                                            \
                                const float momentum,                                          \
                                const float eta,                                               \
                                const float tofEvTime,                                         \
                                const float tofEvTimeErr) -> float {                           \
                               return o2::pid::tof::TOFResponseImpl::nSigma<id>(tofSignal,     \
                                                                                tofExpMom,     \
                                                                                length,        \
                                                                                momentum,      \
                                                                                eta,           \
                                                                                tofEvTime,     \
                                                                                tofEvTimeErr); \
                             });

PERSPECIES_TOF_SEPARATION_COLUMN(El, o2::track::PID::Electron);
PERSPECIES_TOF_SEPARATION_COLUMN(Mu, o2::track::PID::Muon);
PERSPECIES_TOF_SEPARATION_COLUMN(Pi, o2::track::PID::Pion);
PERSPECIES_TOF_SEPARATION_COLUMN(Ka, o2::track::PID::Kaon);
PERSPECIES_TOF_SEPARATION_COLUMN(Pr, o2::track::PID::Proton);
PERSPECIES_TOF_SEPARATION_COLUMN(De, o2::track::PID::Deuteron);
PERSPECIES_TOF_SEPARATION_COLUMN(Tr, o2::track::PID::Triton);
PERSPECIES_TOF_SEPARATION_COLUMN(He, o2::track::PID::Helium3);
PERSPECIES_TOF_SEPARATION_COLUMN(Al, o2::track::PID::Alpha);
#undef PERSPECIES_TOF_SEPARATION_COLUMN

} // namespace pidtof

using TOFExpSigmaDynEl = pidtof::TOFExpSigmaElImp<pidtofsignal::TOFSignal, track::TOFExpMom, track::P, track::Eta, pidtofevtime::TOFEvTimeErr>;
using TOFExpSigmaDynMu = pidtof::TOFExpSigmaMuImp<pidtofsignal::TOFSignal, track::TOFExpMom, track::P, track::Eta, pidtofevtime::TOFEvTimeErr>;
using TOFExpSigmaDynPi = pidtof::TOFExpSigmaPiImp<pidtofsignal::TOFSignal, track::TOFExpMom, track::P, track::Eta, pidtofevtime::TOFEvTimeErr>;
using TOFExpSigmaDynKa = pidtof::TOFExpSigmaKaImp<pidtofsignal::TOFSignal, track::TOFExpMom, track::P, track::Eta, pidtofevtime::TOFEvTimeErr>;
using TOFExpSigmaDynPr = pidtof::TOFExpSigmaPrImp<pidtofsignal::TOFSignal, track::TOFExpMom, track::P, track::Eta, pidtofevtime::TOFEvTimeErr>;
using TOFExpSigmaDynDe = pidtof::TOFExpSigmaDeImp<pidtofsignal::TOFSignal, track::TOFExpMom, track::P, track::Eta, pidtofevtime::TOFEvTimeErr>;
using TOFExpSigmaDynTr = pidtof::TOFExpSigmaTrImp<pidtofsignal::TOFSignal, track::TOFExpMom, track::P, track::Eta, pidtofevtime::TOFEvTimeErr>;
using TOFExpSigmaDynHe = pidtof::TOFExpSigmaHeImp<pidtofsignal::TOFSignal, track::TOFExpMom, track::P, track::Eta, pidtofevtime::TOFEvTimeErr>;
using TOFExpSigmaDynAl = pidtof::TOFExpSigmaAlImp<pidtofsignal::TOFSignal, track::TOFExpMom, track::P, track::Eta, pidtofevtime::TOFEvTimeErr>;

using TOFNSigmaDynEl = pidtof::TOFNSigmaElImp<pidtofsignal::TOFSignal, track::TOFExpMom, track::Length, track::P, track::Eta, pidtofevtime::TOFEvTime, pidtofevtime::TOFEvTimeErr>;
using TOFNSigmaDynMu = pidtof::TOFNSigmaMuImp<pidtofsignal::TOFSignal, track::TOFExpMom, track::Length, track::P, track::Eta, pidtofevtime::TOFEvTime, pidtofevtime::TOFEvTimeErr>;
using TOFNSigmaDynPi = pidtof::TOFNSigmaPiImp<pidtofsignal::TOFSignal, track::TOFExpMom, track::Length, track::P, track::Eta, pidtofevtime::TOFEvTime, pidtofevtime::TOFEvTimeErr>;
using TOFNSigmaDynKa = pidtof::TOFNSigmaKaImp<pidtofsignal::TOFSignal, track::TOFExpMom, track::Length, track::P, track::Eta, pidtofevtime::TOFEvTime, pidtofevtime::TOFEvTimeErr>;
using TOFNSigmaDynPr = pidtof::TOFNSigmaPrImp<pidtofsignal::TOFSignal, track::TOFExpMom, track::Length, track::P, track::Eta, pidtofevtime::TOFEvTime, pidtofevtime::TOFEvTimeErr>;
using TOFNSigmaDynDe = pidtof::TOFNSigmaDeImp<pidtofsignal::TOFSignal, track::TOFExpMom, track::Length, track::P, track::Eta, pidtofevtime::TOFEvTime, pidtofevtime::TOFEvTimeErr>;
using TOFNSigmaDynTr = pidtof::TOFNSigmaTrImp<pidtofsignal::TOFSignal, track::TOFExpMom, track::Length, track::P, track::Eta, pidtofevtime::TOFEvTime, pidtofevtime::TOFEvTimeErr>;
using TOFNSigmaDynHe = pidtof::TOFNSigmaHeImp<pidtofsignal::TOFSignal, track::TOFExpMom, track::Length, track::P, track::Eta, pidtofevtime::TOFEvTime, pidtofevtime::TOFEvTimeErr>;
using TOFNSigmaDynAl = pidtof::TOFNSigmaAlImp<pidtofsignal::TOFSignal, track::TOFExpMom, track::Length, track::P, track::Eta, pidtofevtime::TOFEvTime, pidtofevtime::TOFEvTimeErr>;

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
DECLARE_SOA_DYNAMIC_COLUMN(TOFNSigmaEl, tofNSigmaEl, //! Unwrapped (float) nsigma with the TOF detector for electron
                           [](binning::binned_t nsigma_binned) -> float { return binning::unPackInTable(nsigma_binned); });
DECLARE_SOA_DYNAMIC_COLUMN(TOFNSigmaMu, tofNSigmaMu, //! Unwrapped (float) nsigma with the TOF detector for muon
                           [](binning::binned_t nsigma_binned) -> float { return binning::unPackInTable(nsigma_binned); });
DECLARE_SOA_DYNAMIC_COLUMN(TOFNSigmaPi, tofNSigmaPi, //! Unwrapped (float) nsigma with the TOF detector for pion
                           [](binning::binned_t nsigma_binned) -> float { return binning::unPackInTable(nsigma_binned); });
DECLARE_SOA_DYNAMIC_COLUMN(TOFNSigmaKa, tofNSigmaKa, //! Unwrapped (float) nsigma with the TOF detector for kaon
                           [](binning::binned_t nsigma_binned) -> float { return binning::unPackInTable(nsigma_binned); });
DECLARE_SOA_DYNAMIC_COLUMN(TOFNSigmaPr, tofNSigmaPr, //! Unwrapped (float) nsigma with the TOF detector for proton
                           [](binning::binned_t nsigma_binned) -> float { return binning::unPackInTable(nsigma_binned); });
DECLARE_SOA_DYNAMIC_COLUMN(TOFNSigmaDe, tofNSigmaDe, //! Unwrapped (float) nsigma with the TOF detector for deuteron
                           [](binning::binned_t nsigma_binned) -> float { return binning::unPackInTable(nsigma_binned); });
DECLARE_SOA_DYNAMIC_COLUMN(TOFNSigmaTr, tofNSigmaTr, //! Unwrapped (float) nsigma with the TOF detector for triton
                           [](binning::binned_t nsigma_binned) -> float { return binning::unPackInTable(nsigma_binned); });
DECLARE_SOA_DYNAMIC_COLUMN(TOFNSigmaHe, tofNSigmaHe, //! Unwrapped (float) nsigma with the TOF detector for helium3
                           [](binning::binned_t nsigma_binned) -> float { return binning::unPackInTable(nsigma_binned); });
DECLARE_SOA_DYNAMIC_COLUMN(TOFNSigmaAl, tofNSigmaAl, //! Unwrapped (float) nsigma with the TOF detector for alpha
                           [](binning::binned_t nsigma_binned) -> float { return binning::unPackInTable(nsigma_binned); });

} // namespace pidtof_tiny

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

namespace pidtofbeta
{
DECLARE_SOA_COLUMN(Beta, beta, float);           //! TOF beta
DECLARE_SOA_COLUMN(BetaError, betaerror, float); //! Uncertainty on the TOF beta
// Dynamic column, i.e. the future
DECLARE_SOA_DYNAMIC_COLUMN(TOFBetaImp, tofBeta, //! TOF Beta value
                           [](const float length, const float tofSignal, const float collisionTime) -> float {
                             return o2::pid::tof::Beta::GetBeta(length, tofSignal, collisionTime);
                           });
//
DECLARE_SOA_COLUMN(ExpBetaEl, expbetael, float);           //! Expected beta of electron
DECLARE_SOA_COLUMN(ExpBetaElError, expbetaelerror, float); //! Expected uncertainty on the beta of electron
//
DECLARE_SOA_COLUMN(SeparationBetaEl, separationbetael, float); //! Separation computed with the expected beta for electrons
DECLARE_SOA_DYNAMIC_COLUMN(DiffBetaEl, diffbetael,             //! Difference between the measured and the expected beta for electrons
                           [](float beta, float expbetael) -> float { return beta - expbetael; });
} // namespace pidtofbeta

using TOFBeta = pidtofbeta::TOFBetaImp<o2::aod::track::Length, o2::aod::pidtofsignal::TOFSignal, o2::aod::pidtofevtime::TOFEvTime>;

DECLARE_SOA_TABLE(pidTOFbeta, "AOD", "pidTOFbeta", //! Table of the TOF beta
                  pidtofbeta::Beta, pidtofbeta::BetaError);

namespace pidtofmass
{
DECLARE_SOA_COLUMN(TOFMass, mass, float); //! TOF mass
// Dynamic column, i.e. the future
DECLARE_SOA_DYNAMIC_COLUMN(TOFMassImp, tofMass, //! TOF Mass value
                           [](const float length, const float tofSignal, const float collisionTime, const float momentum) -> float {
                             const float beta = o2::pid::tof::Beta::GetBeta(length, tofSignal, collisionTime);
                             return o2::pid::tof::TOFMass::GetTOFMass(momentum, beta);
                           });
} // namespace pidtofmass

using TOFMass = pidtofmass::TOFMassImp<o2::aod::track::Length, o2::aod::pidtofsignal::TOFSignal, o2::aod::pidtofevtime::TOFEvTime, o2::aod::track::TOFExpMom>;

DECLARE_SOA_TABLE(pidTOFmass, "AOD", "pidTOFmass", //! Table of the TOF mass
                  pidtofmass::TOFMass);

} // namespace o2::aod

#endif // COMMON_DATAMODEL_PIDRESPONSETOF_H_
