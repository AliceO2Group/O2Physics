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

#include "FlatLutEntry.h"

#include <Framework/Logger.h>

#include <cstring>

namespace o2::delphes
{

float map_t::fracPositionWithinBin(float val) const
{
  float width = (max - min) / nbins;
  int bin;
  float returnVal = 0.5f;
  if (log) {
    bin = static_cast<int>((std::log10(val) - min) / width);
    returnVal = ((std::log10(val) - min) / width) - bin;
  } else {
    bin = static_cast<int>((val - min) / width);
    returnVal = val / width - bin;
  }
  return returnVal;
}

int map_t::find(float val) const
{
  float width = (max - min) / nbins;
  int bin;
  if (log) {
    bin = static_cast<int>((std::log10(val) - min) / width);
  } else {
    bin = static_cast<int>((val - min) / width);
  }
  if (bin < 0) {
    return 0;
  }
  if (bin > nbins - 1) {
    return nbins - 1;
  }
  return bin;
}

void map_t::print() const
{
  LOGF(info, "nbins = %d, min = %f, max = %f, log = %s \n", nbins, min, max, log ? "on" : "off");
}

bool lutHeader_t::check_version() const
{
  return (version == LUTCOVM_VERSION);
}

void lutHeader_t::print() const
{
  LOGF(info, " version: %d \n", version);
  LOGF(info, "     pdg: %d \n", pdg);
  LOGF(info, "   field: %f \n", field);
  LOGF(info, "  nchmap: ");
  nchmap.print();
  LOGF(info, "  radmap: ");
  radmap.print();
  LOGF(info, "  etamap: ");
  etamap.print();
  LOGF(info, "   ptmap: ");
  ptmap.print();
}

void FlatLutData::initialize(const lutHeader_t& header)
{
  mNchBins = header.nchmap.nbins;
  mRadBins = header.radmap.nbins;
  mEtaBins = header.etamap.nbins;
  mPtBins = header.ptmap.nbins;

  size_t headerSize = sizeof(lutHeader_t);
  size_t numEntries = static_cast<size_t>(mNchBins) * mRadBins * mEtaBins * mPtBins;
  size_t entriesSize = numEntries * sizeof(lutEntry_t);
  size_t totalSize = headerSize + entriesSize;

  mData.resize(totalSize);

  // Write header at the beginning
  std::memcpy(mData.data(), &header, headerSize);
}

size_t FlatLutData::getEntryOffset(int nch_bin, int rad_bin, int eta_bin, int pt_bin) const
{
  size_t headerSize = sizeof(lutHeader_t);

  // Linear index: nch varies slowest, pt varies fastest
  // idx = nch * (rad*eta*pt) + rad * (eta*pt) + eta * pt + pt
  size_t linearIdx = static_cast<size_t>(nch_bin) * (mRadBins * mEtaBins * mPtBins) + static_cast<size_t>(rad_bin) * (mEtaBins * mPtBins) + static_cast<size_t>(eta_bin) * mPtBins + static_cast<size_t>(pt_bin);

  return headerSize + linearIdx * sizeof(lutEntry_t);
}

lutEntry_t* FlatLutData::getEntry(int nch_bin, int rad_bin, int eta_bin, int pt_bin)
{
  size_t offset = getEntryOffset(nch_bin, rad_bin, eta_bin, pt_bin);
  return reinterpret_cast<lutEntry_t*>(mData.data() + offset);
}

const lutEntry_t* FlatLutData::getEntry(int nch_bin, int rad_bin, int eta_bin, int pt_bin) const
{
  size_t offset = getEntryOffset(nch_bin, rad_bin, eta_bin, pt_bin);
  return reinterpret_cast<const lutEntry_t*>(mData.data() + offset);
}

FlatLutData FlatLutData::fromBuffer(const uint8_t* buffer, size_t size)
{
  FlatLutData data;
  // Validate buffer
  if (size < sizeof(lutHeader_t)) {
    LOG(fatal) << "Buffer too small for LUT header";
  }

  const auto* header = reinterpret_cast<const lutHeader_t*>(buffer);
  data.mNchBins = header->nchmap.nbins;
  data.mRadBins = header->radmap.nbins;
  data.mEtaBins = header->etamap.nbins;
  data.mPtBins = header->ptmap.nbins;

  size_t expectedSize = sizeof(lutHeader_t) + static_cast<size_t>(data.mNchBins) * data.mRadBins * data.mEtaBins * data.mPtBins * sizeof(lutEntry_t);

  if (size < expectedSize) {
    LOG(fatal) << "Buffer size mismatch: expected " << expectedSize << ", got " << size;
  }

  // Copy buffer
  data.mData.resize(size);
  std::memcpy(data.mData.data(), buffer, size);

  return data;
}

FlatLutData FlatLutData::fromExternalBuffer(uint8_t* buffer, size_t size)
{
  FlatLutData data;
  // Validate buffer
  if (size < sizeof(lutHeader_t)) {
    LOG(fatal) << "Buffer too small for LUT header";
  }

  const auto* header = reinterpret_cast<const lutHeader_t*>(buffer);
  data.mNchBins = header->nchmap.nbins;
  data.mRadBins = header->radmap.nbins;
  data.mEtaBins = header->etamap.nbins;
  data.mPtBins = header->ptmap.nbins;

  size_t expectedSize = sizeof(lutHeader_t) + static_cast<size_t>(data.mNchBins) * data.mRadBins * data.mEtaBins * data.mPtBins * sizeof(lutEntry_t);

  if (size < expectedSize) {
    LOG(fatal) << "Buffer size mismatch: expected " << expectedSize << " got " << size;
  }

  // Store reference to external buffer (no copy)
  // WARNING: Caller must ensure buffer lifetime exceeds FlatLutData usage
  data.mData.assign(buffer, buffer + size);

  return data;
}

} // namespace o2::delphes
