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

#ifndef O2_ANALYSIS_PARAMETER_HOLDER_H_
#define O2_ANALYSIS_PARAMETER_HOLDER_H_

#include <Rtypes.h>

// object to hold customizable analysis parameters
class anaparHolder
{
 public:
  // constructor
  anaparHolder(int NCombine = 0,
               std::vector<float> TPCnSigmas = {0.}) : mNCombine{NCombine}, mTPCnSigmas{TPCnSigmas}
  {
    // definition of maximum 10 particles
    // for each particle 12 parameters
    //     1: PID
    //     2: sign
    //  3, 4: min/max nsigma for e
    //  5, 6: min/max nsigma for pi
    //  7, 8: min/max nsigma for mu
    //  9,10: min/max nsigma for Ka
    // 11,12: min/max nsigma for Pr
    //
    // the limits (min/max) of particle PID are inclusive, for all other particles the limits are exclusive
    //
    mTPCnSigmas.resize(120, 0.);
  }

  // getter
  int nCombine() const;
  std::vector<float> TPCnSigmas() const;

 private:
  // number of tracks to combine
  int mNCombine;

  // PID information
  std::vector<float> mTPCnSigmas;

  ClassDefNV(anaparHolder, 1);
};

#endif // O2_ANALYSIS_PARAMETER_HOLDER_H_
