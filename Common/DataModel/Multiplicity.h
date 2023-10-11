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
DECLARE_SOA_DYNAMIC_COLUMN(IsInelGt0, isInelGt0, //! is INEL > 0
                           [](int multPveta1) -> bool { return multPveta1 > 0; });
DECLARE_SOA_DYNAMIC_COLUMN(IsInelGt1, isInelGt1, //! is INEL > 1
                           [](int multPveta1) -> bool { return multPveta1 > 1; });

// complementary / MultsExtra table
DECLARE_SOA_COLUMN(MultPVTotalContributors, multPVTotalContributors, float); //!
DECLARE_SOA_COLUMN(MultPVChi2, multPVChi2, float);                           //!
DECLARE_SOA_COLUMN(MultCollisionTimeRes, multCollisionTimeRes, float);       //!
DECLARE_SOA_COLUMN(MultRunNumber, multRunNumber, int);                       //!
DECLARE_SOA_COLUMN(MultPVz, multPVz, float);                                 //!
DECLARE_SOA_COLUMN(MultSel8, multSel8, bool);                                //!

DECLARE_SOA_COLUMN(MultNTracksHasITS, multNTracksHasITS, int); //!
DECLARE_SOA_COLUMN(MultNTracksHasTPC, multNTracksHasTPC, int); //!
DECLARE_SOA_COLUMN(MultNTracksHasTOF, multNTracksHasTOF, int); //!
DECLARE_SOA_COLUMN(MultNTracksHasTRD, multNTracksHasTRD, int); //!

// further QA
DECLARE_SOA_COLUMN(MultNTracksITSOnly, multNTracksITSOnly, int); //!
DECLARE_SOA_COLUMN(MultNTracksTPCOnly, multNTracksTPCOnly, int); //!
DECLARE_SOA_COLUMN(MultNTracksITSTPC, multNTracksITSTPC, int);   //!

DECLARE_SOA_COLUMN(BCNumber, bcNumber, int); //!

} // namespace mult
DECLARE_SOA_TABLE(FV0Mults, "AOD", "FV0MULT", //! Multiplicity with the FV0 detector
                  mult::MultFV0A, mult::MultFV0C,
                  mult::MultFV0M<mult::MultFV0A, mult::MultFV0C>);
DECLARE_SOA_TABLE(FT0Mults, "AOD", "FT0MULT", //! Multiplicity with the FT0 detector
                  mult::MultFT0A, mult::MultFT0C,
                  mult::MultFT0M<mult::MultFT0A, mult::MultFT0C>);
DECLARE_SOA_TABLE(FDDMults, "AOD", "FDDMULT", //! Multiplicity with the FDD detector
                  mult::MultFDDA, mult::MultFDDC,
                  mult::MultFDDM<mult::MultFDDA, mult::MultFDDC>);
DECLARE_SOA_TABLE(ZDCMults, "AOD", "ZDCMULT", //! Multiplicity with the ZDC detector
                  mult::MultZNA, mult::MultZNC);
DECLARE_SOA_TABLE(BarrelMults, "AOD", "BARRELMULT", //! Multiplicity in the barrel
                  mult::MultTracklets,
                  mult::MultTPC,
                  mult::MultNTracksPV,
                  mult::MultNTracksPVeta1,
                  mult::MultNTracksPVetaHalf,
                  mult::IsInelGt0<mult::MultNTracksPVeta1>,
                  mult::IsInelGt1<mult::MultNTracksPVeta1>);
using Mults = soa::Join<BarrelMults, FV0Mults, FT0Mults, FDDMults, ZDCMults>;
using Mult = Mults::iterator;

// for QA purposes
DECLARE_SOA_TABLE(MultsExtra, "AOD", "MULTEXTRA", //!
                  mult::MultPVTotalContributors, mult::MultPVChi2, mult::MultCollisionTimeRes, mult::MultRunNumber, mult::MultPVz, mult::MultSel8,
                  mult::MultNTracksHasITS, mult::MultNTracksHasTPC, mult::MultNTracksHasTOF, mult::MultNTracksHasTRD,
                  mult::MultNTracksITSOnly, mult::MultNTracksTPCOnly, mult::MultNTracksITSTPC, mult::BCNumber);
using MultExtra = MultsExtra::iterator;

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

namespace multBC
{
DECLARE_SOA_COLUMN(MultBCFT0A, multBCFT0A, float);    //!
DECLARE_SOA_COLUMN(MultBCFT0C, multBCFT0C, float);    //!
DECLARE_SOA_COLUMN(MultBCFV0A, multBCFV0A, float);    //!
DECLARE_SOA_COLUMN(MultBCTVX, multBCTVX, bool);       //!
DECLARE_SOA_COLUMN(MultBCFV0OrA, multBCFV0OrA, bool); //!
} // namespace multDebug
DECLARE_SOA_TABLE(MultsBC, "AOD", "MULTBC", //!
                  multBC::MultBCFT0A,
                  multBC::MultBCFT0C,
                  multBC::MultBCFV0A,
                  multBC::MultBCTVX,
                  multBC::MultBCFV0OrA);
using MultBC = MultsBC::iterator;

} // namespace o2::aod

#endif // O2_ANALYSIS_MULTIPLICITY_H_
