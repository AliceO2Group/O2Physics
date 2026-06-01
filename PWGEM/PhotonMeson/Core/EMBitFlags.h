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

/// \file EMBitFlags.h
/// \brief Header of bit flag class for particle selection.
/// \author M. Hemmer, marvin.hemmer@cern.ch

#ifndef PWGEM_PHOTONMESON_CORE_EMBITFLAGS_H_
#define PWGEM_PHOTONMESON_CORE_EMBITFLAGS_H_

#include <cstddef> // size_t
#include <cstdint> // uint64_t
#include <vector>  // vector

/// \class EMBitFlags
/// \brief Dynamically-sized bit container with bit-level storage.
///
/// Bits can be set beyond the current size, in which case the container
/// grows automatically. Access via test() and reset() requires the index
/// to be within the current size. Bits are all on by default and will be
/// switched off if particle fails a cut
class EMBitFlags
{
 public:
  explicit EMBitFlags(std::size_t nBits = 0);

  /// \brief get number of stored bits
  std::size_t size() const;

  /// \brief check bit i
  /// \param i index of bit that should be checked
  bool test(std::size_t i) const;

  /// \brief set bit i
  /// \param i index of bit which value should be set
  void set(std::size_t i);

  /// \brief reset bit i
  /// \param i index of bit which value should be reset
  void reset(std::size_t i);

  /// \brief resetting all flags to false
  void clear();

  /// \brief reserve space in the underlying storage for nBits bits
  /// \param nBits number of bits to reserve capacity for
  void reserve(std::size_t nBits);

  /// \brief resize the container to hold nBits bits
  /// \param nBits new number of bits (new bits are initialized to false)
  void resize(std::size_t nBits);

 private:
  /// \brief ensure that the container can hold at least nBits bits
  /// \param nBits required number of bits
  void ensureSize(std::size_t nBits);

  static constexpr std::size_t word(std::size_t i) { return i >> 6; }
  static constexpr std::uint64_t mask(std::size_t i) { return 1ULL << (i & 63); }

  std::vector<std::uint64_t> mBits;
  std::size_t mSize = 0;
};

#endif // PWGEM_PHOTONMESON_CORE_EMBITFLAGS_H_
