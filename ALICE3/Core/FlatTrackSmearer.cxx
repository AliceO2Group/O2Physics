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

#include "ALICE3/Core/GeometryContainer.h"

#include <Framework/RuntimeError.h>

namespace o2::delphes
{
int TrackSmearer::getIndexPDG(int pdg) const
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

const char* TrackSmearer::getParticleName(int pdg) const
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

bool TrackSmearer::loadTable(int pdg, const char* filename, bool forceReload)
{
  if (!filename || filename[0] == '\0') {
    LOGF(info, "No LUT file provided for PDG %d. Skipping load.", pdg);
    return false;
  }

  const auto ipdg = getIndexPDG(pdg);
  LOGF(info, "Loading %s LUT file: '%s'", getParticleName(pdg), filename);

  if (mLUTData[ipdg].isLoaded() && !forceReload) {
    LOGF(info, "LUT table for PDG %d already loaded (index %d)", pdg, ipdg);
    return false;
  }

  const std::string localFilename = o2::fastsim::GeometryEntry::accessFile(filename, "./.ALICE3/LUTs/", mCcdbManager, 10);

  try {
    mLUTData[ipdg] = FlatLutData::loadFromFile(localFilename.c_str());
  } catch (framework::RuntimeErrorRef ref) {
    LOGF(error, "%s", framework::error_from_ref(ref).what);
    return false;
  }

  // Validate header
  const auto& header = mLUTData[ipdg].getHeaderRef();

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

  if (header.pdg != pdg && !specialPdgCase) {
    LOGF(error, "LUT header PDG mismatch: expected %d, got %d", pdg, header.pdg);
    mLUTData[ipdg].reset();
    return false;
  }

  LOGF(info, "Successfully read LUT for PDG %d: %s", pdg, localFilename.c_str());
  header.print();

  return true;
}

bool TrackSmearer::hasTable(int pdg) const
{
  const int ipdg = getIndexPDG(pdg);
  return mLUTData[ipdg].isLoaded();
}

} // namespace o2::delphes
