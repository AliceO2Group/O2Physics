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
//
/// \file DalitzBinning.h
/// \brief Binning of the generated pi0/eta Dalitz decays (pT, |y|, m_ee) shared by the producer and the analysis. Both sides must use these bins.
/// \author Joshua Koenig (joshua.konig@cern.ch)

#ifndef PWGEM_PHOTONMESON_UTILS_DALITZBINNING_H_
#define PWGEM_PHOTONMESON_UTILS_DALITZBINNING_H_

#include <vector>

namespace o2::aod::pwgem::photonmeson::utils::dalitzbinning
{
inline std::vector<double> ptBins()
{
  std::vector<double> bins;
  bins.reserve(34);
  for (int i = 0; i < 16; i++) { // o2-linter: disable=magic-number (just numbers for binning)
    bins.emplace_back(0.25 * i); // 0 to 3.75 GeV/c, every 0.25 GeV/c
  }
  for (int i = 0; i < 12; i++) {      // o2-linter: disable=magic-number (just numbers for binning)
    bins.emplace_back(4.0 + 0.5 * i); // 4 to 9.5 GeV/c, every 0.5 GeV/c
  }
  for (int i = 0; i < 6; i++) {        // o2-linter: disable=magic-number (just numbers for binning)
    bins.emplace_back(10.0 + 2.0 * i); // 10 to 20 GeV/c, every 2 GeV/c
  }
  return bins;
}

inline std::vector<double> yBins()
{
  return {0.0, 0.8, 0.9};
}

inline std::vector<double> meeBins()
{
  std::vector<double> bins;
  bins.reserve(46);
  for (int i = 0; i < 20; i++) {  // o2-linter: disable=magic-number (just numbers for binning)
    bins.emplace_back(0.005 * i); // 0 to 0.095 GeV/c^2, every 0.005 GeV/c^2
  }
  for (int i = 0; i < 26; i++) {       // o2-linter: disable=magic-number (just numbers for binning)
    bins.emplace_back(0.1 + 0.02 * i); // 0.1 to 0.6 GeV/c^2, every 0.02 GeV/c^2
  }
  return bins;
}
} // namespace o2::aod::pwgem::photonmeson::utils::dalitzbinning

#endif // PWGEM_PHOTONMESON_UTILS_DALITZBINNING_H_
