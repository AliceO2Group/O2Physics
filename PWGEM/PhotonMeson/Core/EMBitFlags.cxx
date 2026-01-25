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

/// \file EMBitFlags.cxx
/// \brief Source of bit flag class for particle selection.
/// \author M. Hemmer, marvin.hemmer@cern.ch

#include "EMBitFlags.h"

#include <algorithm> // fill
#include <cassert>   // assert
#include <cstddef>   // size_t

EMBitFlags::EMBitFlags(std::size_t nBits)
  : mBits((nBits + 63) / 64, ~0ULL),
    mSize(nBits)
{
}

std::size_t EMBitFlags::size() const
{
  return mSize;
}

bool EMBitFlags::test(std::size_t i) const
{
  assert(i < mSize);
  return mBits[word(i)] & mask(i);
}

void EMBitFlags::set(std::size_t i)
{
  ensureSize(i + 1);
  mBits[word(i)] &= ~mask(i);
}

void EMBitFlags::reset(std::size_t i)
{
  assert(i < mSize);
  mBits[word(i)] |= mask(i);
}

void EMBitFlags::clear()
{
  std::fill(mBits.begin(), mBits.end(), ~0ULL);
}

void EMBitFlags::reserve(std::size_t nBits)
{
  mBits.reserve((nBits + 63) / 64);
}

void EMBitFlags::resize(std::size_t nBits)
{
  mBits.resize((nBits + 63) / 64, ~0ULL);
  mSize = nBits;
}

void EMBitFlags::ensureSize(std::size_t nBits)
{
  if (nBits > mSize) {
    resize(nBits);
  }
}
