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
#ifndef O2_ANALYSIS_MULTIPLICITY_H_
#define O2_ANALYSIS_MULTIPLICITY_H_

#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
namespace mult
{
DECLARE_SOA_COLUMN(MultFV0A, multFV0A, float); //!
DECLARE_SOA_COLUMN(MultFV0C, multFV0C, float); //!
DECLARE_SOA_COLUMN(MultFT0A, multFT0A, float); //!
DECLARE_SOA_COLUMN(MultFT0C, multFT0C, float); //!
DECLARE_SOA_COLUMN(MultFDDA, multFDDA, float); //!
DECLARE_SOA_COLUMN(MultFDDC, multFDDC, float); //!
DECLARE_SOA_COLUMN(MultZNA, multZNA, float);   //!
DECLARE_SOA_COLUMN(MultZNC, multZNC, float);   //!
DECLARE_SOA_DYNAMIC_COLUMN(MultFV0M, multFV0M, //!
                           [](float multFV0A, float multFV0C) -> float { return multFV0A + multFV0C; });
DECLARE_SOA_DYNAMIC_COLUMN(MultFT0M, multFT0M, //!
                           [](float multFT0A, float multFT0C) -> float { return multFT0A + multFT0C; });
DECLARE_SOA_DYNAMIC_COLUMN(MultFDDM, multFDDM, //!
                           [](float multFDDA, float multFDDC) -> float { return multFDDA + multFDDC; });
DECLARE_SOA_COLUMN(MultTracklets, multTracklets, int);
DECLARE_SOA_COLUMN(MultTPC, multTPC, int);
DECLARE_SOA_COLUMN(MultNTracksPV, multNTracksPV, int);
DECLARE_SOA_COLUMN(MultNTracksPVeta1, multNTracksPVeta1, int);
DECLARE_SOA_COLUMN(MultNTracksPVetaHalf, multNTracksPVetaHalf, int);

} // namespace mult
DECLARE_SOA_TABLE(Mults, "AOD", "MULT", //!
                  mult::MultFV0A, mult::MultFV0C,
                  mult::MultFT0A, mult::MultFT0C,
                  mult::MultFDDA, mult::MultFDDC,
                  mult::MultZNA, mult::MultZNC,
                  mult::MultFV0M<mult::MultFV0A, mult::MultFV0C>,
                  mult::MultFT0M<mult::MultFT0A, mult::MultFT0C>,
                  mult::MultFDDM<mult::MultFDDA, mult::MultFDDC>,
                  mult::MultTracklets,
                  mult::MultTPC,
                  mult::MultNTracksPV,
                  mult::MultNTracksPVeta1,
                  mult::MultNTracksPVetaHalf);
using Mult = Mults::iterator;

namespace multZeq
{
DECLARE_SOA_COLUMN(MultZeqFV0A, multZeqFV0A, float);           //!
DECLARE_SOA_COLUMN(MultZeqFT0A, multZeqFT0A, float);           //!
DECLARE_SOA_COLUMN(MultZeqFT0C, multZeqFT0C, float);           //!
DECLARE_SOA_COLUMN(MultZeqFDDA, multZeqFDDA, float);           //!
DECLARE_SOA_COLUMN(MultZeqFDDC, multZeqFDDC, float);           //!
DECLARE_SOA_COLUMN(MultZeqNTracksPV, multZeqNTracksPV, float); //!
} // namespace multZeq
DECLARE_SOA_TABLE(MultZeqs, "AOD", "MULTZEQ", //!
                  multZeq::MultZeqFV0A,
                  multZeq::MultZeqFT0A, multZeq::MultZeqFT0C,
                  multZeq::MultZeqFDDA, multZeq::MultZeqFDDC,
                  multZeq::MultZeqNTracksPV);
using MultZeq = MultZeqs::iterator;
} // namespace o2::aod

#endif // O2_ANALYSIS_MULTIPLICITY_H_
