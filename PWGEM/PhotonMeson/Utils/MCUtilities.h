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

/// \commonly used for MC analysis.
/// \author daiki.sekihata@cern.ch

#ifndef PWGEM_PHOTONMESON_UTILS_MCUTILITIES_H_
#define PWGEM_PHOTONMESON_UTILS_MCUTILITIES_H_

#include "Framework/AnalysisTask.h"

//_______________________________________________________________________
template <typename T, typename TMCs>
bool IsPhysicalPrimary(T const& mctrack, TMCs const& mcTracks)
{
  // This is to check mctrack is ALICE physical primary.
  // https://inspirehep.net/files/4c26ef5fb432df99bdc1ff847653502f

  if (!mctrack.producedByGenerator())
    return false;
  float r3D = sqrt(pow(mctrack.vx(), 2) + pow(mctrack.vy(), 2) + pow(mctrack.vz(), 2)); // cm
  if (r3D > 1.0)
    return false;

  // exclude weak decay. K0S is the most relevant strange particle for neutral mesons.
  if (mctrack.has_mothers()) {
    for (auto& m : mctrack.mothersIds()) {
      if (m < mcTracks.size()) { // protect against bad mother indices
        auto mp = mcTracks.iteratorAt(m);
        int pdg_mother = mp.pdgCode();
        if (pdg_mother == 310 || pdg_mother == 3122) {
          return false;
        }
      }
    }
  }
  return true;
}
//_______________________________________________________________________
template <typename T, typename TMCs>
int IsEleFromPC(T const& mctrack, TMCs const& mcTracks)
{
  // is election from photon conversion? returns index of mother photon
  if (abs(mctrack.pdgCode()) != 11)
    return -1;
  if (mctrack.producedByGenerator())
    return -1;
  if (mctrack.has_mothers()) {
    int motherid = mctrack.mothersIds()[0]; // first mother
    auto mp = mcTracks.iteratorAt(motherid);
    int pdg_mother = mp.pdgCode();
    if (pdg_mother == 22) {
      return motherid;
    }
  } else {
    return -1;
  }
  return -1;
}
//_______________________________________________________________________
//_______________________________________________________________________
#endif // PWGEM_PHOTONMESON_UTILS_MCUTILITIES_H_
