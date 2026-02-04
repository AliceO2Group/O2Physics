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

/// \file TrackSelection.h
/// \brief helper functions for pair track selection
/// \author felix.schlepper@cern.ch

#ifndef PWGEM_PHOTONMESON_UTILS_TRACKSELECTION_H_
#define PWGEM_PHOTONMESON_UTILS_TRACKSELECTION_H_

#include <Framework/ASoA.h>

#include <TPDGCode.h>

#include <cmath>

namespace o2::pwgem::photonmeson
{

/**
 * @brief Track has ITS and TPC
 *
 * @tparam TTrack track
 * @param track track
 * @return true if has both
 */
template <o2::soa::is_iterator TTrack>
inline bool isITSTPCTrack(TTrack const& track)
{
  return track.hasITS() && track.hasTPC();
}

/**
 * @brief Track has TPC and TRD
 *
 * @tparam TTrack track
 * @param track track
 * @return true if has both
 */
template <o2::soa::is_iterator TTrack>
inline bool isTPCTRDTrack(TTrack const& track)
{
  return !track.hasITS() && track.hasTPC() && track.hasTRD() && !track.hasTOF();
}

/**
 * @brief Track has ITS,TPC and TRD
 *
 * @tparam TTrack track
 * @param track track
 * @return true if has all
 */
template <o2::soa::is_iterator TTrack>
inline bool isITSTPCTRDTrack(TTrack const& track)
{
  return track.hasITS() && track.hasTPC() && track.hasTRD() && !track.hasTOF();
}

/**
 * @brief Track has TPC and TOF
 * @tparam TTrack track
 * @param track track
 * @return true if has both
 */
template <o2::soa::is_iterator TTrack>
inline bool isTPCTOFTrack(TTrack const& track)
{
  return !track.hasITS() && track.hasTPC() && !track.hasTRD() && track.hasTOF();
}

/**
 * @brief Track has TPC,TRD and TOF
 * @tparam TTrack track
 * @param track track
 * @return true if has all
 */
template <o2::soa::is_iterator TTrack>
inline bool isTPCTRDTOFTrack(TTrack const& track)
{
  return !track.hasITS() && track.hasTPC() && track.hasTRD() && track.hasTOF();
}

/**
 * @brief Track has ITS,TPC,TRD and TOF
 * @tparam TTrack track
 * @param track track
 * @return true if has all
 */
template <o2::soa::is_iterator TTrack>
inline bool isITSTPCTRDTOFTrack(TTrack const& track)
{
  return track.hasITS() && track.hasTPC() && track.hasTRD() && track.hasTOF();
}

/**
 * @brief Track is TPC-only
 *
 * @tparam TTrack track
 * @param track track
 * @return true if tracks is TPC-only
 */
template <o2::soa::is_iterator TTrack>
inline bool isTPConlyTrack(TTrack const& track)
{
  return !track.hasITS() && track.hasTPC() && !track.hasTRD() && !track.hasTOF();
}

/**
 * @brief Track is ITS-only
 *
 * @tparam TTrack track
 * @param track track
 * @return true if tracks is ITS-only
 */
template <o2::soa::is_iterator TTrack>
inline bool isITSonlyTrack(TTrack const& track)
{
  return track.hasITS() && !track.hasTPC() && !track.hasTRD() && !track.hasTOF();
}

/**
 * @brief If both V0 pairs are ITSTPC-tracks
 *
 * @tparam TTrack track
 * @param track0 track from daughter 0
 * @param track1 track from daughter 1
 * @return true if V0 pairs are ITSTPC-tracks
 */
template <o2::soa::is_iterator TTrack>
inline bool isITSTPC_ITSTPC(TTrack const& track0, TTrack const& track1)
{
  return isITSTPCTrack(track0) && isITSTPCTrack(track1);
}

/**
 * @brief If one track is TPC-only the other ITSTPC
 *
 * @tparam TTrack track
 * @param track0 track from daughter 0
 * @param track1 track from daughter 1
 * @return true if one is TPC-only and the other ITSTPC
 */
template <o2::soa::is_iterator TTrack>
inline bool isITSTPC_TPConly(TTrack const& track0, TTrack const& track1)
{
  return (isITSTPCTrack(track0) && isTPConlyTrack(track1)) || (isITSTPCTrack(track1) && isTPConlyTrack(track0));
}

/**
 * @brief If one track is ITS-only the other ITSTPC
 *
 * @tparam TTrack track
 * @param track0 track from daughter 0
 * @param track1 track from daughter 1
 * @return true if one is ITS-only and the other ITSTPC
 */
template <o2::soa::is_iterator TTrack>
inline bool isITSTPC_ITSonly(TTrack const& track0, TTrack const& track1)
{
  return (isITSTPCTrack(track0) && isITSonlyTrack(track1)) || (isITSTPCTrack(track1) && isITSonlyTrack(track0));
}

/**
 * @brief If V0 pairs are TPC-only tracks
 *
 * @tparam TTrack track
 * @param track0 track from daughter 0
 * @param track1 track from daughter 1
 * @return true if both are TPC-only tracks
 */
template <o2::soa::is_iterator TTrack>
inline bool isTPConly_TPConly(TTrack const& track0, TTrack const& track1)
{
  return isTPConlyTrack(track0) && isTPConlyTrack(track1);
}

/**
 * @brief If V0 pairs are ITS-only tracks
 *
 * @tparam TTrack track
 * @param track0 track from daughter 0
 * @param track1 track from daughter 1
 * @return true if both are ITS-only tracks
 */
template <o2::soa::is_iterator TTrack>
inline bool isITSonly_ITSonly(TTrack const& track0, TTrack const& track1)
{
  return isITSonlyTrack(track0) && isITSonlyTrack(track1);
}

/**
 * @brief If one V0 pair is ITS-only and the other TPC-only
 *
 * @tparam TTrack track
 * @param track0 track from daughter 0
 * @param track1 track from daughter 1
 * @return true if either one is ITS-only while the other one is TPC-only
 */
template <o2::soa::is_iterator TTrack>
inline bool isTPConly_ITSonly(TTrack const& track0, TTrack const& track1)
{
  return (isTPConlyTrack(track0) && isITSonlyTrack(track1)) || (isTPConlyTrack(track1) && isITSonlyTrack(track0));
}

/**
 * @brief Check if MC particles have the expected&same mother particle
 *
 * @param mc1 MCParticle 0
 * @param mc2 MCParticle 1
 * @return true if the mother particle is the expected type and the same for both
 */
template <PDG_t motherType, o2::soa::is_iterator T>
inline bool checkMCParticles(T const& mc1, T const& mc2)
{
  if (std::abs(mc1.pdgCode()) != kElectron || std::abs(mc2.pdgCode()) != kElectron) {
    return false;
  }
  if (!mc1.has_mothers() || !mc2.has_mothers()) {
    return false;
  }
  if (mc1.mothersIds()[0] != mc2.mothersIds()[0]) {
    return false;
  }
  if (const auto& mothers = mc1.template mothers_as<typename T::parent_t>(); mothers.empty() || mothers.size() > 1) {
    return false;
  }
  if (mc1.template mothers_first_as<typename T::parent_t>().pdgCode() != motherType) {
    return false;
  }
  return true;
}
} // namespace o2::pwgem::photonmeson

#endif // PWGEM_PHOTONMESON_UTILS_TRACKSELECTION_H_
