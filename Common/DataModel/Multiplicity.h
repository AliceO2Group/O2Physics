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
#ifndef COMMON_DATAMODEL_MULTIPLICITY_H_
#define COMMON_DATAMODEL_MULTIPLICITY_H_

#include "Framework/AnalysisDataModel.h"
#include "Common/DataModel/EventSelection.h"

namespace o2::aod
{
namespace mult
{
DECLARE_SOA_COLUMN(MultFV0A, multFV0A, float);           //!
DECLARE_SOA_COLUMN(MultFV0AOuter, multFV0AOuter, float); //!
DECLARE_SOA_COLUMN(MultFV0C, multFV0C, float);           //!
DECLARE_SOA_COLUMN(MultFT0A, multFT0A, float);           //!
DECLARE_SOA_COLUMN(MultFT0C, multFT0C, float);           //!
DECLARE_SOA_COLUMN(MultFDDA, multFDDA, float);           //!
DECLARE_SOA_COLUMN(MultFDDC, multFDDC, float);           //!
DECLARE_SOA_COLUMN(MultZNA, multZNA, float);             //!
DECLARE_SOA_COLUMN(MultZNC, multZNC, float);             //!
DECLARE_SOA_COLUMN(MultZEM1, multZEM1, float);           //!
DECLARE_SOA_COLUMN(MultZEM2, multZEM2, float);           //!
DECLARE_SOA_COLUMN(MultZPA, multZPA, float);             //!
DECLARE_SOA_COLUMN(MultZPC, multZPC, float);             //!
DECLARE_SOA_DYNAMIC_COLUMN(MultFV0M, multFV0M,           //!
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

// forward track counters
DECLARE_SOA_COLUMN(MFTNalltracks, mftNalltracks, int); //! overall counter, uses AO2D coll assoc
DECLARE_SOA_COLUMN(MFTNtracks, mftNtracks, int);       //! reassigned, uses mult group software

// MC
DECLARE_SOA_COLUMN(MultMCFT0A, multMCFT0A, int);                       //!
DECLARE_SOA_COLUMN(MultMCFT0C, multMCFT0C, int);                       //!
DECLARE_SOA_COLUMN(MultMCFV0A, multMCFV0A, int);                       //!
DECLARE_SOA_COLUMN(MultMCFDDA, multMCFDDA, int);                       //!
DECLARE_SOA_COLUMN(MultMCFDDC, multMCFDDC, int);                       //!
DECLARE_SOA_COLUMN(MultMCNParticlesEta10, multMCNParticlesEta10, int); //!
DECLARE_SOA_COLUMN(MultMCNParticlesEta08, multMCNParticlesEta08, int); //!
DECLARE_SOA_COLUMN(MultMCNParticlesEta05, multMCNParticlesEta05, int); //!
DECLARE_SOA_COLUMN(MultMCPVz, multMCPVz, float);                       //!

// complementary / MultsExtra table
DECLARE_SOA_COLUMN(MultPVTotalContributors, multPVTotalContributors, int); //!
DECLARE_SOA_COLUMN(MultPVChi2, multPVChi2, float);                         //!
DECLARE_SOA_COLUMN(MultCollisionTimeRes, multCollisionTimeRes, float);     //!
DECLARE_SOA_COLUMN(MultRunNumber, multRunNumber, int);                     //!
DECLARE_SOA_COLUMN(MultPVz, multPVz, float);                               //!
DECLARE_SOA_COLUMN(MultSel8, multSel8, bool);                              //!

DECLARE_SOA_COLUMN(MultNTracksHasITS, multNTracksHasITS, int); //!
DECLARE_SOA_COLUMN(MultNTracksHasTPC, multNTracksHasTPC, int); //!
DECLARE_SOA_COLUMN(MultNTracksHasTOF, multNTracksHasTOF, int); //!
DECLARE_SOA_COLUMN(MultNTracksHasTRD, multNTracksHasTRD, int); //!

// further QA
DECLARE_SOA_COLUMN(MultNTracksITSOnly, multNTracksITSOnly, int);     //!
DECLARE_SOA_COLUMN(MultNTracksTPCOnly, multNTracksTPCOnly, int);     //!
DECLARE_SOA_COLUMN(MultNTracksITSTPC, multNTracksITSTPC, int);       //!
DECLARE_SOA_COLUMN(MultAllTracksTPCOnly, multAllTracksTPCOnly, int); //!
DECLARE_SOA_COLUMN(MultAllTracksITSTPC, multAllTracksITSTPC, int);   //!
DECLARE_SOA_COLUMN(MultNTracksGlobal, multNTracksGlobal, int);       //!
DECLARE_SOA_COLUMN(MultNGlobalTracksPV, multNGlobalTracksPV, int);
DECLARE_SOA_COLUMN(MultNGlobalTracksPVeta1, multNGlobalTracksPVeta1, int);
DECLARE_SOA_COLUMN(MultNGlobalTracksPVetaHalf, multNGlobalTracksPVetaHalf, int);

// even further QA: timing information for neighboring events
DECLARE_SOA_COLUMN(TimeToPrePrevious, timeToPrePrevious, float); //!
DECLARE_SOA_COLUMN(TimeToPrevious, timeToPrevious, float);       //!
DECLARE_SOA_COLUMN(TimeToNext, timeToNext, float);               //!
DECLARE_SOA_COLUMN(TimeToNeNext, timeToNeNext, float);           //!

} // namespace mult
DECLARE_SOA_TABLE(FV0Mults, "AOD", "FV0MULT", //! Multiplicity with the FV0 detector
                  mult::MultFV0A, mult::MultFV0C,
                  mult::MultFV0M<mult::MultFV0A, mult::MultFV0C>);
DECLARE_SOA_TABLE(FV0AOuterMults, "AOD", "FVOAOUTERMULT", //! FV0 without innermost ring
                  mult::MultFV0AOuter);
DECLARE_SOA_TABLE(FT0Mults, "AOD", "FT0MULT", //! Multiplicity with the FT0 detector
                  mult::MultFT0A, mult::MultFT0C,
                  mult::MultFT0M<mult::MultFT0A, mult::MultFT0C>);
DECLARE_SOA_TABLE(FDDMults, "AOD", "FDDMULT", //! Multiplicity with the FDD detector
                  mult::MultFDDA, mult::MultFDDC,
                  mult::MultFDDM<mult::MultFDDA, mult::MultFDDC>);
DECLARE_SOA_TABLE(ZDCMults, "AOD", "ZDCMULT", //! Multiplicity with the ZDC detector
                  mult::MultZNA, mult::MultZNC, mult::MultZEM1, mult::MultZEM2, mult::MultZPA, mult::MultZPC);
DECLARE_SOA_TABLE(TrackletMults, "AOD", "TRKLTMULT", //! Multiplicity with tracklets (only Run2)
                  mult::MultTracklets);
DECLARE_SOA_TABLE(TPCMults, "AOD", "TPCMULT", //! Multiplicity with TPC
                  mult::MultTPC);
DECLARE_SOA_TABLE(PVMults, "AOD", "PVMULT", //! Multiplicity from the PV contributors
                  mult::MultNTracksPV,
                  mult::MultNTracksPVeta1,
                  mult::MultNTracksPVetaHalf,
                  mult::IsInelGt0<mult::MultNTracksPVeta1>,
                  mult::IsInelGt1<mult::MultNTracksPVeta1>);
DECLARE_SOA_TABLE(MFTMults, "AOD", "MFTMULT", //! Multiplicity with MFT
                  mult::MFTNalltracks, mult::MFTNtracks);
using BarrelMults = soa::Join<TrackletMults, TPCMults, PVMults>;
using Mults = soa::Join<BarrelMults, FV0Mults, FT0Mults, FDDMults, ZDCMults>;
using FT0Mult = FT0Mults::iterator;
using MFTMult = MFTMults::iterator;
using Mult = Mults::iterator;

DECLARE_SOA_TABLE(MultsExtra_000, "AOD", "MULTEXTRA", //!
                  mult::MultPVTotalContributors, mult::MultPVChi2, mult::MultCollisionTimeRes, mult::MultRunNumber, mult::MultPVz, mult::MultSel8,
                  mult::MultNTracksHasITS, mult::MultNTracksHasTPC, mult::MultNTracksHasTOF, mult::MultNTracksHasTRD,
                  mult::MultNTracksITSOnly, mult::MultNTracksTPCOnly, mult::MultNTracksITSTPC,
                  mult::MultAllTracksTPCOnly, mult::MultAllTracksITSTPC,
                  evsel::NumTracksInTimeRange,
                  collision::Flags);

DECLARE_SOA_TABLE_VERSIONED(MultsExtra_001, "AOD", "MULTEXTRA", 1, //! debug information
                            mult::MultPVTotalContributors, mult::MultPVChi2, mult::MultCollisionTimeRes, mult::MultRunNumber, mult::MultPVz, mult::MultSel8,
                            mult::MultNTracksHasITS, mult::MultNTracksHasTPC, mult::MultNTracksHasTOF, mult::MultNTracksHasTRD,
                            mult::MultNTracksITSOnly, mult::MultNTracksTPCOnly, mult::MultNTracksITSTPC,
                            mult::MultAllTracksTPCOnly, mult::MultAllTracksITSTPC,
                            evsel::NumTracksInTimeRange,
                            evsel::SumAmpFT0CInTimeRange,
                            collision::Flags);

using MultsExtra = MultsExtra_001;

DECLARE_SOA_TABLE(MultNeighs, "AOD", "MULTNEIGH", //!
                  mult::TimeToPrePrevious, mult::TimeToPrevious,
                  mult::TimeToNext, mult::TimeToNeNext);

// for QA purposes
DECLARE_SOA_TABLE(MultsGlobal, "AOD", "MULTGLOBAL", //! counters that use Track Selection (optional)
                  mult::MultNTracksGlobal,
                  mult::MultNGlobalTracksPV,
                  mult::MultNGlobalTracksPVeta1,
                  mult::MultNGlobalTracksPVetaHalf);
DECLARE_SOA_TABLE(MultSelections, "AOD", "MULTSELECTIONS", //!
                  evsel::Selection);                       // for derived data / QA studies
using MultExtra = MultsExtra::iterator;

// mc collisions table - indexed to Mult
DECLARE_SOA_TABLE(MultMCExtras_000, "AOD", "MULTMCEXTRA", //! Table for MC information
                  mult::MultMCFT0A,
                  mult::MultMCFT0C,
                  mult::MultMCNParticlesEta05,
                  mult::MultMCNParticlesEta08,
                  mult::MultMCNParticlesEta10,
                  mult::MultMCPVz,
                  mult::IsInelGt0<mult::MultMCNParticlesEta10>,
                  mult::IsInelGt1<mult::MultMCNParticlesEta10>,
                  o2::soa::Marker<1>);

// mc collisions table - indexed to Mult
DECLARE_SOA_TABLE_VERSIONED(MultMCExtras_001, "AOD", "MULTMCEXTRA", 1, //! Table for MC information
                            mult::MultMCFT0A,
                            mult::MultMCFT0C,
                            mult::MultMCFV0A,
                            mult::MultMCFDDA,
                            mult::MultMCFDDC,
                            mult::MultMCNParticlesEta05,
                            mult::MultMCNParticlesEta08,
                            mult::MultMCNParticlesEta10,
                            mult::MultMCPVz,
                            mult::IsInelGt0<mult::MultMCNParticlesEta10>,
                            mult::IsInelGt1<mult::MultMCNParticlesEta10>,
                            o2::soa::Marker<1>);

using MultMCExtras = MultMCExtras_001;
using MultMCExtra = MultMCExtras::iterator;
using MultsExtraMC = MultMCExtras; // for backwards compatibility with previous naming scheme

// crosslinks
namespace mult
{
DECLARE_SOA_INDEX_COLUMN(MultMCExtra, multMCExtra);
}

DECLARE_SOA_TABLE(Mult2MCExtras, "AOD", "Mult2MCEXTRA", //! Relate reco mult entry to MC extras entry
                  o2::soa::Index<>, mult::MultMCExtraId);

namespace multZeq
{
DECLARE_SOA_COLUMN(MultZeqFV0A, multZeqFV0A, float);           //! Multiplicity equalized for the vertex position with the FV0A detector
DECLARE_SOA_COLUMN(MultZeqFT0A, multZeqFT0A, float);           //! Multiplicity equalized for the vertex position with the FT0A detector
DECLARE_SOA_COLUMN(MultZeqFT0C, multZeqFT0C, float);           //! Multiplicity equalized for the vertex position with the FT0C detector
DECLARE_SOA_COLUMN(MultZeqFDDA, multZeqFDDA, float);           //! Multiplicity equalized for the vertex position with the FDDA detector
DECLARE_SOA_COLUMN(MultZeqFDDC, multZeqFDDC, float);           //! Multiplicity equalized for the vertex position with the FDDC detector
DECLARE_SOA_COLUMN(MultZeqNTracksPV, multZeqNTracksPV, float); //! Multiplicity equalized for the vertex position from the PV contributors
} // namespace multZeq
DECLARE_SOA_TABLE(FV0MultZeqs, "AOD", "FV0MULTZEQ", //! Multiplicity equalized for the vertex position with the FV0 detector
                  multZeq::MultZeqFV0A);
DECLARE_SOA_TABLE(FT0MultZeqs, "AOD", "FT0MULTZEQ", //! Multiplicity equalized for the vertex position with the FT0 detector
                  multZeq::MultZeqFT0A, multZeq::MultZeqFT0C);
DECLARE_SOA_TABLE(FDDMultZeqs, "AOD", "FDDMULTZEQ", //! Multiplicity equalized for the vertex position with the FDD detector
                  multZeq::MultZeqFDDA, multZeq::MultZeqFDDC);
DECLARE_SOA_TABLE(PVMultZeqs, "AOD", "PVMULTZEQ", //! Multiplicity equalized for the vertex position from the PV contributors
                  multZeq::MultZeqNTracksPV);
using MultZeqs = soa::Join<FV0MultZeqs, FT0MultZeqs, FDDMultZeqs, PVMultZeqs>;
using MultZeq = MultZeqs::iterator;

namespace mult
{
// extra BC information
DECLARE_SOA_COLUMN(MultTVX, multTVX, bool);                          //!
DECLARE_SOA_COLUMN(MultFV0OrA, multFV0OrA, bool);                    //!
DECLARE_SOA_COLUMN(MultV0triggerBits, multV0triggerBits, uint8_t);   //!
DECLARE_SOA_COLUMN(MultT0triggerBits, multT0triggerBits, uint8_t);   //!
DECLARE_SOA_COLUMN(MultFDDtriggerBits, multFDDtriggerBits, uint8_t); //!
DECLARE_SOA_COLUMN(MultTriggerMask, multTriggerMask, uint64_t);      //! CTP trigger mask
DECLARE_SOA_COLUMN(MultCollidingBC, multCollidingBC, bool);          //! CTP trigger mask

DECLARE_SOA_COLUMN(MultFT0PosZ, multFT0PosZ, float);          //! Position along Z computed with the FT0 information within the BC
DECLARE_SOA_COLUMN(MultFT0PosZValid, multFT0PosZValid, bool); //! Validity of the position along Z computed with the FT0 information
} // namespace mult
DECLARE_SOA_TABLE(MultBCs, "AOD", "MULTBC", //!
                  mult::MultFT0A,
                  mult::MultFT0C,
                  mult::MultFT0PosZ,
                  mult::MultFT0PosZValid,
                  mult::MultFV0A,
                  mult::MultFDDA,
                  mult::MultFDDC,
                  mult::MultZNA,
                  mult::MultZNC,
                  mult::MultZEM1,
                  mult::MultZEM2,
                  mult::MultZPA,
                  mult::MultZPC,
                  mult::MultTVX,
                  mult::MultFV0OrA,
                  mult::MultV0triggerBits,
                  mult::MultT0triggerBits,
                  mult::MultFDDtriggerBits,
                  mult::MultTriggerMask,
                  mult::MultCollidingBC,
                  timestamp::Timestamp,
                  bc::Flags);
using MultBC = MultBCs::iterator;

// crosslinks
namespace mult
{
DECLARE_SOA_INDEX_COLUMN(MultBC, multBC);
}
namespace multBC
{
DECLARE_SOA_INDEX_COLUMN(FT0Mult, ft0Mult);
}

// for QA purposes
DECLARE_SOA_TABLE(Mults2BC, "AOD", "MULTS2BC", //! Relate mult -> BC
                  o2::soa::Index<>, mult::MultBCId);
DECLARE_SOA_TABLE(BC2Mults, "AOD", "BC2MULTS", //! Relate BC -> mult
                  o2::soa::Index<>, multBC::FT0MultId);

} // namespace o2::aod

#endif // COMMON_DATAMODEL_MULTIPLICITY_H_
