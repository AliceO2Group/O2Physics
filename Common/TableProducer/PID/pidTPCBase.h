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

///
/// \file   pidTPCBase.h
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \brief  Base to build tasks for TPC PID tasks.
///

#ifndef COMMON_TABLEPRODUCER_PID_PIDTPCBASE_H_
#define COMMON_TABLEPRODUCER_PID_PIDTPCBASE_H_

#include "Common/DataModel/Multiplicity.h"

#include <Framework/ASoA.h>
#include <ReconstructionDataFormats/PID.h>

#include <TMatrixD.h> // IWYU pragma: keep (do not replace with TMatrixDfwd.h)
#include <TMatrixDfwd.h>

#include <cstdlib>
#include <vector>

namespace o2::aod
{

namespace pid
{
DECLARE_SOA_COLUMN(TpcSignalCorrected, tpcSignalCorrected, float); //!
}; // namespace pid

DECLARE_SOA_TABLE(PIDMults, "AOD", "PIDMults", //! TPC auxiliary table for the PID
                  o2::soa::Marker<1>,
                  mult::MultTPC);
DECLARE_SOA_TABLE_FULL(DEdxsCorrected, "DEdxsCorrected", "AOD", "DEDXCORR", pid::TpcSignalCorrected); //!
using PIDMult = PIDMults::iterator;
using DEdxCorrected = DEdxsCorrected::iterator; //!

} // namespace o2::aod

int getPIDIndex(const int pdgCode) // Get O2 PID index corresponding to MC PDG code
{
  switch (abs(pdgCode)) {
    case 11:
      return o2::track::PID::Electron;
    case 13:
      return o2::track::PID::Muon;
    case 211:
      return o2::track::PID::Pion;
    case 321:
      return o2::track::PID::Kaon;
    case 2212:
      return o2::track::PID::Proton;
    case 1000010020:
      return o2::track::PID::Deuteron;
    case 1000010030:
      return o2::track::PID::Triton;
    case 1000020030:
      return o2::track::PID::Helium3;
    case 1000020040:
      return o2::track::PID::Alpha;
    default: // treat as pion if not any of the above
      return o2::track::PID::Pion;
  }
}

typedef struct Str_dEdx_correction {
  TMatrixD fMatrix;
  bool warning = true;

  // void init(std::vector<double>& params)
  void init()
  {
    double elements[32] = {0.99091, -0.015053, 0.0018912, -0.012305,
                           0.081387, 0.003205, -0.0087404, -0.0028608,
                           0.013066, 0.017012, -0.0018469, -0.0052177,
                           -0.0035655, 0.0017846, 0.0019127, -0.00012964,
                           0.0049428, 0.0055592, -0.0010618, -0.0016134,
                           -0.0059098, 0.0013335, 0.00052133, 3.1119e-05,
                           -0.004882, 0.00077317, -0.0013827, 0.003249,
                           -0.00063689, 0.0016218, -0.00045215, -1.5815e-05};
    fMatrix.ResizeTo(4, 8);
    fMatrix.SetMatrixArray(elements);
  }

  float fReal_fTPCSignalN(std::vector<float> vec1, std::vector<float> vec2)
  {
    float result = 0.f;
    // push 1.
    vec1.insert(vec1.begin(), 1.0);
    vec2.insert(vec2.begin(), 1.0);
    for (int i = 0; i < fMatrix.GetNrows(); i++) {
      for (int j = 0; j < fMatrix.GetNcols(); j++) {
        double param = fMatrix(i, j);
        double value1 = i > static_cast<int>(vec1.size()) ? 0 : vec1[i];
        double value2 = j > static_cast<int>(vec2.size()) ? 0 : vec2[j];
        result += param * value1 * value2;
      }
    }
    return result;
  }
} Str_dEdx_correction;

#endif // COMMON_TABLEPRODUCER_PID_PIDTPCBASE_H_
