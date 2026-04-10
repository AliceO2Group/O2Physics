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
#include <Framework/RuntimeError.h>

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
  updateRef();
}

size_t FlatLutData::getEntryOffset(int nch_bin, int rad_bin, int eta_bin, int pt_bin) const
{
  size_t headerSize = sizeof(lutHeader_t);

  // Linear index: nch varies slowest, pt varies fastest
  // idx = nch * (rad*eta*pt) + rad * (eta*pt) + eta * pt + pt
  size_t linearIdx = static_cast<size_t>(nch_bin) * (mRadBins * mEtaBins * mPtBins) + static_cast<size_t>(rad_bin) * (mEtaBins * mPtBins) + static_cast<size_t>(eta_bin) * mPtBins + static_cast<size_t>(pt_bin);

  return headerSize + linearIdx * sizeof(lutEntry_t);
}

const lutEntry_t* FlatLutData::getEntryRef(int nch_bin, int rad_bin, int eta_bin, int pt_bin) const
{
  size_t offset = getEntryOffset(nch_bin, rad_bin, eta_bin, pt_bin);
  return reinterpret_cast<const lutEntry_t*>(mDataRef.data() + offset);
}

lutEntry_t* FlatLutData::getEntry(int nch_bin, int rad_bin, int eta_bin, int pt_bin)
{
  size_t offset = getEntryOffset(nch_bin, rad_bin, eta_bin, pt_bin);
  return reinterpret_cast<lutEntry_t*>(mData.data() + offset);
}

const lutHeader_t& FlatLutData::getHeaderRef() const
{
  return *reinterpret_cast<const lutHeader_t*>(mDataRef.data());
}

lutHeader_t& FlatLutData::getHeader()
{
  return *reinterpret_cast<lutHeader_t*>(mData.data());
}

void FlatLutData::updateRef()
{
  mDataRef = std::span{mData.data(), mData.size()};
}

void FlatLutData::cacheDimensions()
{
  auto const& header = getHeaderRef();
  mNchBins = header.nchmap.nbins;
  mRadBins = header.radmap.nbins;
  mEtaBins = header.etamap.nbins;
  mPtBins = header.ptmap.nbins;
}

void FlatLutData::resetDimensions()
{
  mNchBins = 0;
  mRadBins = 0;
  mEtaBins = 0;
  mPtBins = 0;
}

void FlatLutData::adopt(const uint8_t* buffer, size_t size)
{
  mData.resize(size);
  std::memcpy(mData.data(), buffer, size);
  updateRef();
  cacheDimensions();
}

void FlatLutData::view(const uint8_t* buffer, size_t size)
{
  mData.clear();
  mDataRef = std::span{buffer, size};
  cacheDimensions();
}

void FlatLutData::validateBuffer(const uint8_t* buffer, size_t size)
{
  auto header = PreviewHeader(buffer, size);
  if (!header.check_version()) {
    throw framework::runtime_error_f("LUT header version mismatch: expected %d, got %d", LUTCOVM_VERSION, header.version);
  }
  auto mNchBins = header.nchmap.nbins;
  auto mRadBins = header.radmap.nbins;
  auto mEtaBins = header.etamap.nbins;
  auto mPtBins = header.ptmap.nbins;

  size_t expectedSize = sizeof(lutHeader_t) + static_cast<size_t>(mNchBins) * mRadBins * mEtaBins * mPtBins * sizeof(lutEntry_t);

  if (size < expectedSize) {
    throw framework::runtime_error_f("Buffer size mismatch: expected %zu, got %zu", expectedSize, size);
  }
}

lutHeader_t FlatLutData::PreviewHeader(const uint8_t* buffer, size_t size)
{
  if (size < sizeof(lutHeader_t)) {
    throw framework::runtime_error_f("Buffer too small for LUT header: expected at least %zu, got %zu", sizeof(lutHeader_t), size);
  }
  const auto* header = reinterpret_cast<const lutHeader_t*>(buffer);
  if (!header->check_version()) {
    throw framework::runtime_error_f("LUT header version mismatch: expected %d, got %d", LUTCOVM_VERSION, header->version);
  }
  return *header;
}

FlatLutData FlatLutData::AdoptFromBuffer(const uint8_t* buffer, size_t size)
{
  validateBuffer(buffer, size);
  FlatLutData data;

  // Copy buffer
  data.adopt(buffer, size);
  return data;
}

FlatLutData FlatLutData::ViewFromBuffer(const uint8_t* buffer, size_t size)
{
  validateBuffer(buffer, size);
  FlatLutData data;

  // Store reference to external buffer
  // WARNING: Caller must ensure buffer lifetime exceeds FlatLutData usage
  data.view(buffer, size);
  return data;
}

bool FlatLutData::isLoaded() const
{
  return ((!mData.empty()) || (!mDataRef.empty()));
}

lutHeader_t FlatLutData::PreviewHeader(std::ifstream& file, const char* filename)
{
  lutHeader_t tempHeader;
  file.read(reinterpret_cast<char*>(&tempHeader), sizeof(lutHeader_t));
  if (file.gcount() != static_cast<std::streamsize>(sizeof(lutHeader_t))) {
    throw framework::runtime_error_f("Failed to read LUT header from %s", filename);
  }
  if (!tempHeader.check_version()) {
    throw framework::runtime_error_f("LUT header version mismatch: expected %d, got %d", LUTCOVM_VERSION, tempHeader.version);
  }
  return tempHeader;
}

FlatLutData FlatLutData::loadFromFile(std::ifstream& file, const char* filename)
{
  // Read header first
  lutHeader_t tempHeader = PreviewHeader(file, filename);

  FlatLutData data;

  // Initialize flat data structure
  data.initialize(tempHeader);

  // Read all entries sequentially into flat buffer
  size_t headerSize = sizeof(lutHeader_t);
  size_t numEntries = static_cast<size_t>(data.mNchBins) * data.mRadBins * data.mEtaBins * data.mPtBins;
  size_t entriesSize = numEntries * sizeof(lutEntry_t);

  file.read(reinterpret_cast<char*>(data.data() + headerSize), entriesSize);
  if (file.gcount() != static_cast<std::streamsize>(entriesSize)) {
    throw framework::runtime_error_f("Failed to read LUT entries from %s: expected %zu bytes, got %zu", filename, entriesSize, static_cast<size_t>(file.gcount()));
  }

  LOGF(info, "Successfully loaded LUT from %s: %zu entries", filename, numEntries);
  return data;
}

void FlatLutData::reset()
{
  mData.clear();
  updateRef();
  resetDimensions();
}

} // namespace o2::delphes
