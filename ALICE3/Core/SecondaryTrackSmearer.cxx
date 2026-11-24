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

#include "ALICE3/Core/SecondaryTrackSmearer.h"

#include "ALICE3/Core/ConfigurationParser.h"
#include "ALICE3/Core/GeometryContainer.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/Logger.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>

namespace o2::fastsim
{

namespace
{
int particleKey(int pdg)
{
  const auto magnitude = std::abs(static_cast<std::int64_t>(pdg));
  return magnitude <= std::numeric_limits<int>::max() ? static_cast<int>(magnitude) : 0;
}

float clampToAxis(const SecondaryLutAxis& axis, float value)
{
  return std::clamp(value, axis.lowEdge(0), axis.upEdge(axis.nbins - 1));
}
} // namespace

bool SecondaryTrackSmearer::loadTable(int pdg, const char* filename, bool forceReload, const char* objectName)
{
  if (!filename || !filename[0] || !objectName || !objectName[0]) {
    LOGF(error, "A filename and object name are required for secondary LUT PDG %d", pdg);
    return false;
  }
  if (particleKey(pdg) == 0) {
    LOGF(error, "Invalid secondary LUT PDG %d", pdg);
    return false;
  }
  if (hasTable(pdg) && !forceReload) {
    LOGF(info, "Secondary LUT for PDG %d already loaded", pdg);
    return false;
  }
  const std::string path(filename);
  if (path.rfind("ccdb:", 0) == 0 && !mCcdbManager) {
    LOGF(error, "A CCDB manager is required to load secondary LUT %s", filename);
    return false;
  }
  const auto localFile = ConfigurationParser::accessFile(path, "./.ALICE3/LUTs/", mCcdbManager,
                                                         GeometryContainer::cleanLutWhenLoaded() ? 10 : 0);
  std::unique_ptr<SecondaryLut> table(SecondaryLut::readFromFile(localFile.c_str(), objectName));
  if (!table) {
    LOGF(error, "Cannot read secondary LUT object '%s' from '%s'", objectName, localFile.c_str());
    return false;
  }
  return installTable(pdg, std::move(table), forceReload);
}

bool SecondaryTrackSmearer::adoptTable(int pdg, const SecondaryLut& table, bool forceReload)
{
  if (hasTable(pdg) && !forceReload) {
    LOGF(info, "Secondary LUT for PDG %d already loaded", pdg);
    return false;
  }
  return installTable(pdg, std::make_unique<SecondaryLut>(table), forceReload);
}

bool SecondaryTrackSmearer::installTable(int pdg, std::unique_ptr<SecondaryLut> table, bool forceReload)
{
  const int key = particleKey(pdg);
  const int tableKey = particleKey(table->getPdg());
  if (key == 0 || (key != tableKey && !(key == o2::constants::physics::kAlpha && tableKey == o2::constants::physics::kHelium3))) {
    LOGF(error, "Secondary LUT PDG mismatch: requested %d, table has %d", pdg, table->getPdg());
    return false;
  }
  if (hasTable(pdg) && !forceReload) {
    LOGF(info, "Secondary LUT for PDG %d already loaded", pdg);
    return false;
  }
  std::size_t expectedSize = 1;
  for (int d = 0; d < SecondaryLut::kNDims; ++d) {
    const auto& axis = table->getAxis(static_cast<SecondaryLut::Dim>(d));
    if (axis.nbins <= 0 || !std::isfinite(axis.min) || !std::isfinite(axis.max) ||
        !std::isfinite(axis.max - axis.min) || axis.max <= axis.min ||
        !std::isfinite(axis.lowEdge(0)) || !std::isfinite(axis.upEdge(axis.nbins - 1)) ||
        (axis.log && axis.lowEdge(0) <= 0.f) ||
        expectedSize > std::numeric_limits<std::size_t>::max() / static_cast<std::size_t>(axis.nbins)) {
      LOGF(error, "Invalid axis %d in secondary LUT for PDG %d", d, pdg);
      return false;
    }
    expectedSize *= static_cast<std::size_t>(axis.nbins);
  }
  if (!table->isInitialized() || table->size() != expectedSize) {
    LOGF(error, "Invalid entry count in secondary LUT for PDG %d", pdg);
    return false;
  }
  mTables[key] = std::move(table);
  LOGF(info, "Loaded secondary LUT for PDG %d", pdg);
  return true;
}

bool SecondaryTrackSmearer::hasTable(int pdg) const
{
  return getLUT(pdg) != nullptr;
}

const SecondaryLut* SecondaryTrackSmearer::getLUT(int pdg) const
{
  const auto it = mTables.find(particleKey(pdg));
  return it == mTables.end() ? nullptr : it->second.get();
}

void SecondaryTrackSmearer::setWhatEfficiency(EfficiencyType value)
{
  if (value != kWhatEfficiencyReco && value != kWhatEfficiencyRecoAndTOF) {
    throw std::invalid_argument("Unknown secondary LUT efficiency type");
  }
  mWhatEfficiency = value;
  mSmearer.setWhatEfficiency(value);
}

float SecondaryTrackSmearer::efficiency(const SecondaryLutEntry& entry) const
{
  return mWhatEfficiency == kWhatEfficiencyReco ? entry.eff : entry.eff2;
}

const SecondaryLutEntry* SecondaryTrackSmearer::getLUTEntry(int pdg, float nch, float radius, float eta, float pt, int nHits, float& interpolatedEff) const
{
  interpolatedEff = std::numeric_limits<float>::quiet_NaN();
  const auto* table = getLUT(pdg);
  if (!table) {
    LOGF(warn, "No secondary LUT loaded for PDG %d", pdg);
    return nullptr;
  }
  if (nHits < 0 || !std::isfinite(nch) || !std::isfinite(radius) || !std::isfinite(eta) || !std::isfinite(pt)) {
    LOGF(error, "Invalid secondary LUT lookup values for PDG %d", pdg);
    return nullptr;
  }
  nch = clampToAxis(table->getAxis(SecondaryLut::kNch), nch);
  radius = clampToAxis(table->getAxis(SecondaryLut::kRad), radius);
  eta = clampToAxis(table->getAxis(SecondaryLut::kEta), eta);
  pt = clampToAxis(table->getAxis(SecondaryLut::kPt), pt);
  const float hits = clampToAxis(table->getAxis(SecondaryLut::kNHits), static_cast<float>(nHits));
  const auto* entry = table->lookup(nch, radius, eta, pt, hits);
  interpolatedEff = efficiency(*entry);
  if (!mInterpolateEfficiency || !entry->valid) {
    return entry;
  }

  const auto& axis = table->getAxis(SecondaryLut::kNch);
  const float value = axis.log ? (nch > 0.f ? std::log10(nch) : axis.min) : nch;
  const float position = std::clamp((value - axis.min) / (axis.max - axis.min) * axis.nbins - 0.5f, 0.f, static_cast<float>(axis.nbins - 1));
  const int lower = static_cast<int>(position);
  const int upper = std::min(lower + 1, axis.nbins - 1);
  const int irad = table->getAxis(SecondaryLut::kRad).find(radius);
  const int ieta = table->getAxis(SecondaryLut::kEta).find(eta);
  const int ipt = table->getAxis(SecondaryLut::kPt).find(pt);
  const int ihits = table->getAxis(SecondaryLut::kNHits).find(hits);
  const auto& low = table->at(lower, irad, ieta, ipt, ihits);
  const auto& high = table->at(upper, irad, ieta, ipt, ihits);
  if (!low.valid || !high.valid) {
    LOGF(warn, "Cannot interpolate secondary LUT efficiency across an invalid bin for PDG %d", pdg);
    return entry;
  }
  const float fraction = position - lower;
  interpolatedEff = (1.f - fraction) * efficiency(low) + fraction * efficiency(high);
  return entry;
}

bool SecondaryTrackSmearer::smearTrack(O2Track& track, const SecondaryLutEntry* entry, float interpolatedEff, bool* wasSmeared)
{
  if (wasSmeared) {
    *wasSmeared = false;
  }
  if (!entry || !entry->valid) {
    LOG(warn) << "Cannot smear a track with a missing or invalid secondary LUT entry";
    return false;
  }
  const float selectedEfficiency = mInterpolateEfficiency ? interpolatedEff : efficiency(*entry);
  if (mUseEfficiency && (!std::isfinite(selectedEfficiency) || selectedEfficiency < 0.f || selectedEfficiency > 1.f)) {
    LOG(error) << "Cannot smear a track with an invalid efficiency";
    return false;
  }
  const auto finite = [](float value) { return std::isfinite(value); };
  if (!std::all_of(entry->covm.begin(), entry->covm.end(), finite) ||
      !std::all_of(entry->eigval.begin(), entry->eigval.end(), finite) ||
      !std::all_of(entry->eigvec.begin(), entry->eigvec.end(), finite) ||
      !std::all_of(entry->eiginv.begin(), entry->eiginv.end(), finite)) {
    LOG(error) << "Cannot smear a track with non-finite secondary LUT parameters";
    return false;
  }
  lutEntry_t flatEntry;
  flatEntry.valid = entry->valid;
  flatEntry.eff = entry->eff;
  flatEntry.eff2 = entry->eff2;
  flatEntry.covm = entry->covm;
  flatEntry.eigval = entry->eigval;
  for (int i = 0; i < kNumEigenModes; ++i) {
    for (int j = 0; j < kNumEigenModes; ++j) {
      flatEntry.eigvec[i][j] = entry->eigvec[i * kNumEigenModes + j];
      flatEntry.eiginv[i][j] = entry->eiginv[i * kNumEigenModes + j];
    }
  }
  const bool reconstructed = mSmearer.smearTrack(track, &flatEntry, interpolatedEff);
  if (wasSmeared) {
    *wasSmeared = reconstructed || !mSkipUnreconstructed;
  }
  return reconstructed;
}

bool SecondaryTrackSmearer::smearTrack(O2Track& track, int pdg, float nch, float radius, int nHits, bool* wasSmeared)
{
  float pt = track.getPt();
  const auto key = particleKey(pdg);
  if (key == o2::constants::physics::kHelium3 || key == o2::constants::physics::kAlpha) {
    pt *= 2.f;
  }
  float interpolatedEff;
  const auto* entry = getLUTEntry(pdg, nch, radius, track.getEta(), pt, nHits, interpolatedEff);
  return smearTrack(track, entry, interpolatedEff, wasSmeared);
}

double SecondaryTrackSmearer::getEfficiency(int pdg, float nch, float radius, float eta, float pt, int nHits) const
{
  float result;
  const auto* entry = getLUTEntry(pdg, nch, radius, eta, pt, nHits, result);
  if (!entry || !entry->valid || !std::isfinite(result) || result < 0.f || result > 1.f) {
    LOGF(warn, "No valid secondary LUT efficiency for PDG %d", pdg);
    return std::numeric_limits<double>::quiet_NaN();
  }
  return result;
}

double SecondaryTrackSmearer::getPtRes(int pdg, float nch, float radius, float eta, float pt, int nHits) const
{
  float ignored;
  const auto* entry = getLUTEntry(pdg, nch, radius, eta, pt, nHits, ignored);
  if (!entry || !entry->valid) {
    LOGF(warn, "No valid secondary LUT resolution for PDG %d", pdg);
    return std::numeric_limits<double>::quiet_NaN();
  }
  if (!std::isfinite(entry->covm[14]) || entry->covm[14] < 0.f) {
    LOGF(error, "Invalid inverse-pT variance in secondary LUT for PDG %d", pdg);
    return std::numeric_limits<double>::quiet_NaN();
  }
  const auto& axis = getLUT(pdg)->getAxis(SecondaryLut::kPt);
  return std::sqrt(entry->covm[14]) * axis.center(axis.find(clampToAxis(axis, pt)));
}

double SecondaryTrackSmearer::getAbsPtRes(int pdg, float nch, float radius, float eta, float pt, int nHits) const
{
  const double relative = getPtRes(pdg, nch, radius, eta, pt, nHits);
  if (!std::isfinite(relative)) {
    return relative;
  }
  const auto& axis = getLUT(pdg)->getAxis(SecondaryLut::kPt);
  return relative * axis.center(axis.find(clampToAxis(axis, pt)));
}

double SecondaryTrackSmearer::getAbsEtaRes(int pdg, float nch, float radius, float eta, float pt, int nHits) const
{
  float ignored;
  const auto* entry = getLUTEntry(pdg, nch, radius, eta, pt, nHits, ignored);
  if (!entry || !entry->valid) {
    LOGF(warn, "No valid secondary LUT resolution for PDG %d", pdg);
    return std::numeric_limits<double>::quiet_NaN();
  }
  if (!std::isfinite(entry->covm[9]) || entry->covm[9] < 0.f) {
    LOGF(error, "Invalid tan(lambda) variance in secondary LUT for PDG %d", pdg);
    return std::numeric_limits<double>::quiet_NaN();
  }
  return std::fabs(std::sin(2.0 * std::atan(std::exp(-eta)))) * std::sqrt(entry->covm[9]);
}

double SecondaryTrackSmearer::getEtaRes(int pdg, float nch, float radius, float eta, float pt, int nHits) const
{
  const double absolute = getAbsEtaRes(pdg, nch, radius, eta, pt, nHits);
  if (!std::isfinite(absolute)) {
    return absolute;
  }
  const auto& axis = getLUT(pdg)->getAxis(SecondaryLut::kEta);
  const float binEta = axis.center(axis.find(clampToAxis(axis, eta)));
  if (binEta == 0.f) {
    LOGF(warn, "Cannot compute relative eta resolution for zero eta bin center (PDG %d)", pdg);
    return std::numeric_limits<double>::quiet_NaN();
  }
  return absolute / binEta;
}

} // namespace o2::fastsim
