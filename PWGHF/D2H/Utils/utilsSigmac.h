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

/// \file utilsSigmac.h
/// \brief Utilities for Sigmac analysis
/// \author Mattia Faggin <mattia.faggin@cern.ch>, INFN Padova

#ifndef PWGHF_D2H_UTILS_UTILSSIGMAC_H_
#define PWGHF_D2H_UTILS_UTILSSIGMAC_H_

#include "PWGHF/DataModel/CandidateReconstructionTables.h"

#include <Rtypes.h>

#include <cstdint>

namespace o2::hf_sigmac_utils
{
/// @brief Function to determine if the reconstructed candidate Σc0,++ decays into Λc+ → pK-π+, Λc+ → π+K-p or both
/// @tparam L template for Lc daughter of Sc candidate
/// @tparam S template for Sc candidate
/// @param candidateLc Lc daughter of Sc candidate
/// @param candSc Sc candidate
/// @return 0: none; 1: only Λc+ → pK-π+ possible; 2: Λc+ → π+K-p possible; 3: both possible
template <typename L, typename S>
int8_t isDecayToPKPiToPiKP(L& candidateLc, S& candSc)
{
  int8_t channel = 0;
  if ((candidateLc.isSelLcToPKPi() >= 1) && candSc.statusSpreadLcMinvPKPiFromPDG()) {
    // Λc+ → pK-π+ and within the requested mass to build the Σc0,++
    SETBIT(channel, o2::aod::hf_cand_sigmac::Decays::PKPi);
  }
  if ((candidateLc.isSelLcToPiKP() >= 1) && candSc.statusSpreadLcMinvPiKPFromPDG()) {
    // Λc+ → π+K-p and within the requested mass to build the Σc0,++
    SETBIT(channel, o2::aod::hf_cand_sigmac::Decays::PiKP);
  }
  return channel; /// 0: none; 1: pK-π+ only; 2: π+K-p only; 3: both possible
}
} // namespace o2::hf_sigmac_utils

#endif // PWGHF_D2H_UTILS_UTILSSIGMAC_H_
