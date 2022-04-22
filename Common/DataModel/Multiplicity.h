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
DECLARE_SOA_COLUMN(MultV0A, multV0A, float);   //!
DECLARE_SOA_COLUMN(MultV0C, multV0C, float);   //!
DECLARE_SOA_COLUMN(MultT0A, multT0A, float);   //!
DECLARE_SOA_COLUMN(MultT0C, multT0C, float);   //!
DECLARE_SOA_COLUMN(MultFDDA, multFDDA, float); //!
DECLARE_SOA_COLUMN(MultFDDC, multFDDC, float); //!
DECLARE_SOA_COLUMN(MultZNA, multZNA, float);   //!
DECLARE_SOA_COLUMN(MultZNC, multZNC, float);   //!
DECLARE_SOA_DYNAMIC_COLUMN(MultV0M, multV0M,   //!
                           [](float multV0A, float multV0C) -> float { return multV0A + multV0C; });
DECLARE_SOA_DYNAMIC_COLUMN(MultT0M, multT0M, //!
                           [](float multT0A, float multT0C) -> float { return multT0A + multT0C; });
DECLARE_SOA_DYNAMIC_COLUMN(MultFDD, multFDD, //!
                           [](float multFDDA, float multFDDC) -> float { return multFDDA + multFDDC; });
DECLARE_SOA_COLUMN(MultTracklets, multTracklets, int);
DECLARE_SOA_COLUMN(MultTPC, multTPC, int);
DECLARE_SOA_COLUMN(MultNContribs, multNContribs, int);
DECLARE_SOA_COLUMN(MultZeqV0A, multZeqV0A, float);   //!
DECLARE_SOA_COLUMN(MultZeqT0A, multZeqT0A, float);   //!
DECLARE_SOA_COLUMN(MultZeqT0C, multZeqT0C, float);   //!
DECLARE_SOA_COLUMN(MultZeqFDDA, multZeqFDDA, float); //!
DECLARE_SOA_COLUMN(MultZeqFDDC, multZeqFDDC, float); //!
DECLARE_SOA_COLUMN(MultZeqNContribs, multZeqNContribs, float); //!

} // namespace mult
DECLARE_SOA_TABLE(Mults, "AOD", "MULT", //!
                  mult::MultV0A, mult::MultV0C,
                  mult::MultT0A, mult::MultT0C,
                  mult::MultFDDA, mult::MultFDDC,
                  mult::MultZNA, mult::MultZNC,
                  mult::MultV0M<mult::MultV0A, mult::MultV0C>,
                  mult::MultT0M<mult::MultT0A, mult::MultT0C>,
                  mult::MultFDD<mult::MultFDDA, mult::MultFDDC>,
                  mult::MultTracklets,
                  mult::MultTPC,
                  mult::MultNContribs,
                  mult::MultZeqV0A,
                  mult::MultZeqT0A,  mult::MultZeqT0C,
                  mult::MultZeqFDDA, mult::MultZeqFDDC,
                  mult::MultZeqNContribs);
using Mult = Mults::iterator;
} // namespace o2::aod

#endif // O2_ANALYSIS_MULTIPLICITY_H_
