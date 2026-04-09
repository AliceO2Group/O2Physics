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

namespace o2::delphes {
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
} // namespace o2::delphes
