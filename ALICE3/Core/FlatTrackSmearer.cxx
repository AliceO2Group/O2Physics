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

#include "FlatTrackSmearer.h"

#include "FlatLutEntry.h"

#include "ALICE3/Core/GeometryContainer.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/RuntimeError.h>

#include <TRandom.h>

#include <fairlogger/Logger.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <span>
#include <string>

namespace o2::delphes
{
int TrackSmearer::getIndexPDG(int pdg)
{
  switch (std::abs(pdg)) {
    case 11:
      return 0; // Electron
    case 13:
      return 1; // Muon
    case 211:
      return 2; // Pion
    case 321:
      return 3; // Kaon
    case 2212:
      return 4; // Proton
    case 1000010020:
      return 5; // Deuteron
    case 1000010030:
      return 6; // Triton
    case 1000020030:
      return 7; // Helium3
    case 1000020040:
      return 8; // Alphas
    default:
      return 2; // Default: pion
  }
}

const char* TrackSmearer::getParticleName(int pdg)
{
  switch (std::abs(pdg)) {
    case 11:
      return "electron";
    case 13:
      return "muon";
    case 211:
      return "pion";
    case 321:
      return "kaon";
    case 2212:
      return "proton";
    case 1000010020:
      return "deuteron";
    case 1000010030:
      return "triton";
    case 1000020030:
      return "helium3";
    case 1000020040:
      return "alpha";
    default:
      return "pion"; // Default: pion
  }
}

void TrackSmearer::setWhatEfficiency(int val)
{
  // FIXME: this really should be an enum
  if (val > 2) {
    throw framework::runtime_error_f("getLUTEntry: unknown efficiency type %d", mWhatEfficiency);
  }
  mWhatEfficiency = val;
}

bool TrackSmearer::loadTable(int pdg, const char* filename, bool forceReload)
{
  if (!filename || filename[0] == '\0') {
    LOGF(info, "No LUT file provided for PDG %d. Skipping load.", pdg);
    return false;
  }

  const auto ipdg = getIndexPDG(pdg);
  if (mLUTData[ipdg].isLoaded() && !forceReload) {
    LOGF(info, "LUT table for PDG %d already loaded (index %d)", pdg, ipdg);
    return false;
  }

  LOGF(info, "Loading %s LUT file: '%s'", getParticleName(pdg), filename);
  const std::string localFilename = o2::fastsim::GeometryEntry::accessFile(filename, "./.ALICE3/LUTs/", mCcdbManager, 10);

  std::ifstream lutFile(localFilename, std::ifstream::binary);
  if (!lutFile.is_open()) {
    throw framework::runtime_error_f("Cannot open LUT file: %s", localFilename.c_str());
  }

  try {
    mLUTData[ipdg] = FlatLutData::loadFromFile(lutFile, localFilename.c_str());

    // Validate header
    auto header = mLUTData[ipdg].getHeader();
    if (header.pdg != pdg && !checkSpecialCase(pdg, header)) {
      LOGF(error, "LUT header PDG mismatch: expected %d, got %d; not loading", pdg, header.pdg);
      return false;
    }
  } catch (framework::RuntimeErrorRef ref) {
    LOGF(error, "%s", framework::error_from_ref(ref).what);
    return false;
  }

  LOGF(info, "Successfully read LUT for PDG %d: %s", pdg, localFilename.c_str());
  mLUTData[ipdg].getHeaderRef().print();
  return true;
}

bool TrackSmearer::adoptTable(int pdg, const uint8_t* buffer, size_t size, bool forceReload)
{
  const auto ipdg = getIndexPDG(pdg);
  if (mLUTData[ipdg].isLoaded() && !forceReload) {
    LOGF(info, "LUT table for PDG %d already loaded (index %d)", pdg, ipdg);
    return false;
  }
  try {
    auto header = FlatLutData::PreviewHeader(buffer, size);
    if (header.pdg != pdg && !checkSpecialCase(pdg, header)) {
      LOGF(error, "LUT header PDG mismatch: expected %d, got %d", pdg, header.pdg);
      return false;
    }
    mLUTData[ipdg] = FlatLutData::AdoptFromBuffer(buffer, size);
  } catch (framework::RuntimeErrorRef ref) {
    LOGF(error, "%s", framework::error_from_ref(ref).what);
  }

  LOGF(info, "Successfully adopted LUT for PDG %d", pdg);
  mLUTData[ipdg].getHeaderRef().print();
  return true;
}

bool TrackSmearer::viewTable(int pdg, const uint8_t* buffer, size_t size, bool forceReload)
{
  const auto ipdg = getIndexPDG(pdg);
  if (mLUTData[ipdg].isLoaded() && !forceReload) {
    LOGF(info, "LUT table for PDG %d already loaded (index %d)", pdg, ipdg);
    return false;
  }
  try {
    auto header = FlatLutData::PreviewHeader(buffer, size);
    if (header.pdg != pdg && !checkSpecialCase(pdg, header)) {
      LOGF(error, "LUT header PDG mismatch: expected %d, got %d", pdg, header.pdg);
      return false;
    }
    mLUTData[ipdg] = FlatLutData::ViewFromBuffer(buffer, size);
  } catch (framework::RuntimeErrorRef ref) {
    LOGF(error, "%s", framework::error_from_ref(ref).what);
  }

  LOGF(info, "Successfully viewing LUT for PDG %d", pdg);
  mLUTData[ipdg].getHeaderRef().print();
  return true;
}

bool TrackSmearer::viewTable(int pdg, std::span<std::byte> const& span, bool forceReload)
{
  return viewTable(pdg, reinterpret_cast<const uint8_t*>(span.data()), span.size_bytes(), forceReload);
}

bool TrackSmearer::hasTable(int pdg) const
{
  const int ipdg = getIndexPDG(pdg);
  return mLUTData[ipdg].isLoaded();
}

bool TrackSmearer::checkSpecialCase(int pdg, lutHeader_t const& header)
{
  // Validate header
  bool specialPdgCase = false;
  switch (pdg) {
    case o2::constants::physics::kAlpha:
      // Special case: Allow Alpha particles to use He3 LUT
      specialPdgCase = (header.pdg == o2::constants::physics::kHelium3);
      if (specialPdgCase) {
        LOGF(info, "Alpha particles (PDG %d) will use He3 LUT data (PDG %d)", pdg, header.pdg);
      }
      break;
    default:
      break;
  }
  return specialPdgCase;
}

const lutHeader_t* TrackSmearer::getLUTHeader(int pdg) const
{
  const int ipdg = getIndexPDG(pdg);
  if (!mLUTData[ipdg].isLoaded()) {
    return nullptr;
  }
  return &mLUTData[ipdg].getHeaderRef();
}

const lutEntry_t* TrackSmearer::getLUTEntry(const int pdg, const float nch, const float radius, const float eta, const float pt, float& interpolatedEff) const
{
  const int ipdg = getIndexPDG(pdg);
  if (!mLUTData[ipdg].isLoaded()) {
    return nullptr;
  }

  const auto& header = mLUTData[ipdg].getHeaderRef();

  auto inch = header.nchmap.find(nch);
  auto irad = header.radmap.find(radius);
  auto ieta = header.etamap.find(eta);
  auto ipt = header.ptmap.find(pt);

  // Interpolate efficiency if requested
  if (mInterpolateEfficiency) {
    auto fraction = header.nchmap.fracPositionWithinBin(nch);
    static constexpr float kFractionThreshold = 0.5f;
    if (fraction > kFractionThreshold) {
      switch (mWhatEfficiency) {
        case 1: {
          const auto* entry_curr = mLUTData[ipdg].getEntryRef(inch, irad, ieta, ipt);
          if (inch < header.nchmap.nbins - 1) {
            const auto* entry_next = mLUTData[ipdg].getEntryRef(inch + 1, irad, ieta, ipt);
            interpolatedEff = (1.5f - fraction) * entry_curr->eff + (-0.5f + fraction) * entry_next->eff;
          } else {
            interpolatedEff = entry_curr->eff;
          }
          break;
        }
        case 2: {
          const auto* entry_curr = mLUTData[ipdg].getEntryRef(inch, irad, ieta, ipt);
          if (inch < header.nchmap.nbins - 1) {
            const auto* entry_next = mLUTData[ipdg].getEntryRef(inch + 1, irad, ieta, ipt);
            interpolatedEff = (1.5f - fraction) * entry_curr->eff2 + (-0.5f + fraction) * entry_next->eff2;
          } else {
            interpolatedEff = entry_curr->eff2;
          }
          break;
        }
      }
    } else {
      float comparisonValue = header.nchmap.log ? std::log10(nch) : nch;
      switch (mWhatEfficiency) {
        case 1: {
          const auto* entry_curr = mLUTData[ipdg].getEntryRef(inch, irad, ieta, ipt);
          if (inch > 0 && comparisonValue < header.nchmap.max) {
            const auto* entry_prev = mLUTData[ipdg].getEntryRef(inch - 1, irad, ieta, ipt);
            interpolatedEff = (0.5f + fraction) * entry_curr->eff + (0.5f - fraction) * entry_prev->eff;
          } else {
            interpolatedEff = entry_curr->eff;
          }
          break;
        }
        case 2: {
          const auto* entry_curr = mLUTData[ipdg].getEntryRef(inch, irad, ieta, ipt);
          if (inch > 0 && comparisonValue < header.nchmap.max) {
            const auto* entry_prev = mLUTData[ipdg].getEntryRef(inch - 1, irad, ieta, ipt);
            interpolatedEff = (0.5f + fraction) * entry_curr->eff2 + (0.5f - fraction) * entry_prev->eff2;
          } else {
            interpolatedEff = entry_curr->eff2;
          }
          break;
        }
      }
    }
  } else {
    const auto* entry = mLUTData[ipdg].getEntryRef(inch, irad, ieta, ipt);
    if (entry) {
      switch (mWhatEfficiency) {
        case 1:
          interpolatedEff = entry->eff;
          break;
        case 2:
          interpolatedEff = entry->eff2;
          break;
      }
    }
  }

  return mLUTData[ipdg].getEntryRef(inch, irad, ieta, ipt);
}

bool TrackSmearer::smearTrack(O2Track& o2track, const lutEntry_t* lutEntry, float interpolatedEff)
{
  bool isReconstructed = true;

  // Generate efficiency
  if (mUseEfficiency) {
    auto eff = 0.f;
    switch (mWhatEfficiency) {
      case 1:
        eff = lutEntry->eff;
        break;
      case 2:
        eff = lutEntry->eff2;
        break;
    }
    if (mInterpolateEfficiency) {
      eff = interpolatedEff;
    }
    if (gRandom->Uniform() > eff) { // FIXME: use a fixed RNG instead of whatever ROOT has as a default
      isReconstructed = false;
    }
  }

  // Return false already now in case not reco'ed
  if (!isReconstructed && mSkipUnreconstructed) {
    return false;
  }

  // Transform params vector and smear
  static constexpr int kParSize = 5;
  double params[kParSize];
  for (int i = 0; i < kParSize; ++i) {
    double val = 0.;
    for (int j = 0; j < kParSize; ++j) {
      val += lutEntry->eigvec[j][i] * o2track.getParam(j);
    }
    params[i] = gRandom->Gaus(val, std::sqrt(lutEntry->eigval[i]));
  }

  // Transform back params vector
  for (int i = 0; i < kParSize; ++i) {
    double val = 0.;
    for (int j = 0; j < kParSize; ++j) {
      val += lutEntry->eiginv[j][i] * params[j];
    }
    o2track.setParam(val, i);
  }

  // Sanity check that par[2] sin(phi) is in [-1, 1]
  if (std::fabs(o2track.getParam(2)) > 1.) {
    LOGF(warn, "smearTrack failed sin(phi) sanity check: %f", o2track.getParam(2));
  }

  // Set covariance matrix
  static constexpr int kCovMatSize = 15;
  for (int i = 0; i < kCovMatSize; ++i) {
    o2track.setCov(lutEntry->covm[i], i);
  }

  return isReconstructed;
}

bool TrackSmearer::smearTrack(O2Track& o2track, int pdg, float nch)
{
  auto pt = o2track.getPt();
  switch (pdg) {
    case o2::constants::physics::kHelium3:
    case -o2::constants::physics::kHelium3:
      pt *= 2.f;
      break;
  }

  auto eta = o2track.getEta();
  float interpolatedEff = 0.0f;
  const lutEntry_t* lutEntry = getLUTEntry(pdg, nch, 0.f, eta, pt, interpolatedEff);

  if (!lutEntry || !lutEntry->valid) {
    return false;
  }

  return smearTrack(o2track, lutEntry, interpolatedEff);
}

double TrackSmearer::getPtRes(const int pdg, const float nch, const float eta, const float pt) const
{
  float dummy = 0.0f;
  const lutEntry_t* lutEntry = getLUTEntry(pdg, nch, 0.f, eta, pt, dummy);
  auto val = std::sqrt(lutEntry->covm[14]) * lutEntry->pt;
  return val;
}

double TrackSmearer::getEtaRes(const int pdg, const float nch, const float eta, const float pt) const
{
  float dummy = 0.0f;
  const lutEntry_t* lutEntry = getLUTEntry(pdg, nch, 0.f, eta, pt, dummy);
  auto sigmatgl = std::sqrt(lutEntry->covm[9]);                                  // sigmatgl2
  auto etaRes = std::fabs(std::sin(2.0 * std::atan(std::exp(-eta)))) * sigmatgl; // propagate tgl to eta uncertainty
  etaRes /= lutEntry->eta;                                                       // relative uncertainty
  return etaRes;
}

double TrackSmearer::getAbsPtRes(const int pdg, const float nch, const float eta, const float pt) const
{
  float dummy = 0.0f;
  const lutEntry_t* lutEntry = getLUTEntry(pdg, nch, 0.f, eta, pt, dummy);
  auto val = std::sqrt(lutEntry->covm[14]) * lutEntry->pt * lutEntry->pt;
  return val;
}

double TrackSmearer::getAbsEtaRes(const int pdg, const float nch, const float eta, const float pt) const
{
  float dummy = 0.0f;
  const lutEntry_t* lutEntry = getLUTEntry(pdg, nch, 0.f, eta, pt, dummy);
  auto sigmatgl = std::sqrt(lutEntry->covm[9]);                                  // sigmatgl2
  auto etaRes = std::fabs(std::sin(2.0 * std::atan(std::exp(-eta)))) * sigmatgl; // propagate tgl to eta uncertainty
  return etaRes;
}

double TrackSmearer::getEfficiency(const int pdg, const float nch, const float eta, const float pt) const
{
  float efficiency = 0.0f;
  (void)getLUTEntry(pdg, nch, 0.f, eta, pt, efficiency);
  return efficiency;
}

} // namespace o2::delphes
