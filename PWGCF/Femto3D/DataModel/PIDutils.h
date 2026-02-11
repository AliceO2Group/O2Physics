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
/// \file   PIDutils.h
/// \author Gleb Romanenko gleb.romanenko@cern.ch
/// \brief  SFINAE checks for existance of PID info
/// \since  24/03/2025
///

#ifndef PWGCF_FEMTO3D_DATAMODEL_PIDUTILS_H_
#define PWGCF_FEMTO3D_DATAMODEL_PIDUTILS_H_

#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"

#include <type_traits>
#include <utility>
#include <vector>

namespace o2::aod::singletrackselector
{
namespace pidutils
{

//========================================== SFINAE checks ==========================================

template <typename T, typename = void>
struct hasTPCPi : std::false_type {
};
template <typename T>
struct hasTPCPi<T, std::void_t<o2::aod::pidutils::hasTPCPi<T>>> : std::true_type {
};

template <typename T, typename = void>
struct hasTPCKa : std::false_type {
};
template <typename T>
struct hasTPCKa<T, std::void_t<o2::aod::pidutils::hasTPCKa<T>>> : std::true_type {
};

template <typename T, typename = void>
struct hasTPCPr : std::false_type {
};
template <typename T>
struct hasTPCPr<T, std::void_t<o2::aod::pidutils::hasTPCPr<T>>> : std::true_type {
};

template <typename T, typename = void>
struct hasTPCDe : std::false_type {
};
template <typename T>
struct hasTPCDe<T, std::void_t<o2::aod::pidutils::hasTPCDe<T>>> : std::true_type {
};

template <typename T, typename = void>
struct hasTPCTr : std::false_type {
};
template <typename T>
struct hasTPCTr<T, std::void_t<o2::aod::pidutils::hasTPCTr<T>>> : std::true_type {
};

template <typename T, typename = void>
struct hasTPCHe : std::false_type {
};
template <typename T>
struct hasTPCHe<T, std::void_t<o2::aod::pidutils::hasTPCHe<T>>> : std::true_type {
};

template <typename T, typename = void>
struct hasTOFPi : std::false_type {
};
template <typename T>
struct hasTOFPi<T, std::void_t<o2::aod::pidutils::hasTOFPi<T>>> : std::true_type {
};

template <typename T, typename = void>
struct hasTOFKa : std::false_type {
};
template <typename T>
struct hasTOFKa<T, std::void_t<o2::aod::pidutils::hasTOFKa<T>>> : std::true_type {
};

template <typename T, typename = void>
struct hasTOFPr : std::false_type {
};
template <typename T>
struct hasTOFPr<T, std::void_t<o2::aod::pidutils::hasTOFPr<T>>> : std::true_type {
};

template <typename T, typename = void>
struct hasTOFDe : std::false_type {
};
template <typename T>
struct hasTOFDe<T, std::void_t<o2::aod::pidutils::hasTOFDe<T>>> : std::true_type {
};

template <typename T, typename = void>
struct hasTOFTr : std::false_type {
};
template <typename T>
struct hasTOFTr<T, std::void_t<o2::aod::pidutils::hasTOFTr<T>>> : std::true_type {
};

template <typename T, typename = void>
struct hasTOFHe : std::false_type {
};
template <typename T>
struct hasTOFHe<T, std::void_t<o2::aod::pidutils::hasTOFHe<T>>> : std::true_type {
};

} // namespace pidutils

//========================================== ITS PID ==========================================

template <typename TrackType>
inline float getITSNsigma(TrackType const& track, int const& PDG)
{
  switch (PDG) {
    case 211:
      return track.itsNSigmaPi();
    case 321:
      return track.itsNSigmaKa();
    case 2212:
      return track.itsNSigmaPr();
    case 1000010020:
      return track.itsNSigmaDe();
    case 1000020030:
      return track.itsNSigmaHe();
    case 1000010030:
      return track.itsNSigmaTr();
    case 0:
      return -1000.0;
    default:
      LOG(fatal) << "Cannot interpret PDG for ITS selection: " << PDG;
      return -1000.0;
  }
}

template <typename TrackType>
inline bool ITSselection(TrackType const& track, std::pair<int, std::vector<float>> const& PIDcuts)
{
  float Nsigma = getITSNsigma(track, PIDcuts.first);

  if (Nsigma > PIDcuts.second[0] && Nsigma < PIDcuts.second[1]) {
    return true;
  }
  return false;
}

//========================================== TPC PID ==========================================

template <typename TrackType>
inline float getTPCNsigma(TrackType const& track, int const& PDG)
{
  switch (PDG) {
    case 211:
      if constexpr (o2::aod::singletrackselector::pidutils::hasTPCPi<TrackType>::value)
        return track.tpcNSigmaPi();
    case 321:
      if constexpr (o2::aod::singletrackselector::pidutils::hasTPCKa<TrackType>::value)
        return track.tpcNSigmaKa();
    case 2212:
      if constexpr (o2::aod::singletrackselector::pidutils::hasTPCPr<TrackType>::value)
        return track.tpcNSigmaPr();
    case 1000010020:
      if constexpr (o2::aod::singletrackselector::pidutils::hasTPCDe<TrackType>::value)
        return track.tpcNSigmaDe();
    case 1000020030:
      if constexpr (o2::aod::singletrackselector::pidutils::hasTPCHe<TrackType>::value)
        return track.tpcNSigmaHe();
    case 1000010030:
      if constexpr (o2::aod::singletrackselector::pidutils::hasTPCTr<TrackType>::value)
        return track.tpcNSigmaTr();
    case 0:
      return -1000.0;
    default:
      LOG(fatal) << "Cannot interpret PDG for TPC selection: " << PDG;
      return -1000.0;
  }
}

template <bool useITS, typename TrackType>
inline bool TPCselection(TrackType const& track, std::pair<int, std::vector<float>> const& PIDcuts, std::vector<float> const& ITSCut = std::vector<float>{})
{
  int PDG = PIDcuts.first;

  if constexpr (useITS) {
    if (ITSCut.size() != 0 && !ITSselection(track, std::make_pair(PDG, ITSCut)))
      return false;
  }

  float Nsigma = getTPCNsigma(track, PDG);

  if (Nsigma > PIDcuts.second[0] && Nsigma < PIDcuts.second[1]) {
    return true;
  }
  return false;
}

//========================================== TOF PID ==========================================

template <typename TrackType>
inline float getTOFNsigma(TrackType const& track, int const& PDG)
{
  switch (PDG) {
    case 211:
      if constexpr (o2::aod::singletrackselector::pidutils::hasTOFPi<TrackType>::value)
        return track.tofNSigmaPi();
    case 321:
      if constexpr (o2::aod::singletrackselector::pidutils::hasTOFKa<TrackType>::value)
        return track.tofNSigmaKa();
    case 2212:
      if constexpr (o2::aod::singletrackselector::pidutils::hasTOFPr<TrackType>::value)
        return track.tofNSigmaPr();
    case 1000010020:
      if constexpr (o2::aod::singletrackselector::pidutils::hasTOFDe<TrackType>::value)
        return track.tofNSigmaDe();
    case 1000020030:
      if constexpr (o2::aod::singletrackselector::pidutils::hasTOFHe<TrackType>::value)
        return track.tofNSigmaHe();
    case 1000010030:
      if constexpr (o2::aod::singletrackselector::pidutils::hasTOFTr<TrackType>::value)
        return track.tofNSigmaTr();
    case 0:
      return -1000.0;
    default:
      LOG(fatal) << "Cannot interpret PDG for TOF selection: " << PDG;
      return -1000.0;
  }
}

template <typename TrackType>
inline bool TOFselection(TrackType const& track, std::pair<int, std::vector<float>> const& PIDcuts, std::vector<float> const& TPCresidualCut = std::vector<float>{-5.0f, 5.0f})
{
  int PDG = PIDcuts.first;
  if (!TPCselection<false>(track, std::make_pair(PDG, TPCresidualCut)))
    return false;

  float Nsigma = getTOFNsigma(track, PDG);

  if (Nsigma > PIDcuts.second[0] && Nsigma < PIDcuts.second[1]) {
    return true;
  }
  return false;
}
} // namespace o2::aod::singletrackselector

#endif // PWGCF_FEMTO3D_DATAMODEL_PIDUTILS_H_
