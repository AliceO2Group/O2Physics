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

/// \file LFEventTopologyTables.h
/// \brief Per-collision forward/backward sub-event tables, reconstructed and generated level.
/// \author Cristian Andrei <Cristian.Andrei@cern.ch>

#ifndef PWGLF_DATAMODEL_LFEVENTTOPOLOGYTABLES_H_
#define PWGLF_DATAMODEL_LFEVENTTOPOLOGYTABLES_H_

#include <Framework/AnalysisDataModel.h>

#include <cstdint>

namespace o2::aod
{

namespace evshapecoex
{
DECLARE_SOA_COLUMN(CollisionId, collisionId, int);    //! source collision global index (DF-local; not an index column)
DECLARE_SOA_COLUMN(RunNumber, runNumber, int);        //! run number
DECLARE_SOA_COLUMN(GlobalBC, globalBC, uint64_t);     //! global bunch crossing
DECLARE_SOA_COLUMN(PosZ, posZ, float);                //! primary-vertex z (cm)
DECLARE_SOA_COLUMN(MultFT0A, multFT0A, float);        //! FT0-A amplitude
DECLARE_SOA_COLUMN(MultFT0C, multFT0C, float);        //! FT0-C amplitude
DECLARE_SOA_COLUMN(NF, nF, uint16_t);                 //! track count, forward sub-event
DECLARE_SOA_COLUMN(NB, nB, uint16_t);                 //! track count, backward sub-event
DECLARE_SOA_COLUMN(NGap, nGap, uint16_t);             //! track count, central gap
DECLARE_SOA_COLUMN(AF, aF, float);                    //! mean pT, forward (GeV/c); NaN if empty
DECLARE_SOA_COLUMN(AB, aB, float);                    //! mean pT, backward (GeV/c); NaN if empty
DECLARE_SOA_COLUMN(SumPt2F, sumPt2F, float);          //! sum of pT^2, forward (GeV^2/c^2)
DECLARE_SOA_COLUMN(SumPt2B, sumPt2B, float);          //! sum of pT^2, backward (GeV^2/c^2)
DECLARE_SOA_COLUMN(QxF, qxF, float);                  //! sum of cos(2 phi), forward
DECLARE_SOA_COLUMN(QyF, qyF, float);                  //! sum of sin(2 phi), forward
DECLARE_SOA_COLUMN(QxB, qxB, float);                  //! sum of cos(2 phi), backward
DECLARE_SOA_COLUMN(QyB, qyB, float);                  //! sum of sin(2 phi), backward
DECLARE_SOA_COLUMN(Qx4F, qx4F, float);                //! sum of cos(4 phi), forward
DECLARE_SOA_COLUMN(Qy4F, qy4F, float);                //! sum of sin(4 phi), forward
DECLARE_SOA_COLUMN(Qx4B, qx4B, float);                //! sum of cos(4 phi), backward
DECLARE_SOA_COLUMN(Qy4B, qy4B, float);                //! sum of sin(4 phi), backward
DECLARE_SOA_COLUMN(NPlus, nPlus, uint16_t);           //! positive tracks, forward + backward
DECLARE_SOA_COLUMN(NMinus, nMinus, uint16_t);         //! negative tracks, forward + backward
DECLARE_SOA_COLUMN(QaBits, qaBits, uint16_t);         //! event-selection bits, order as QaBitOrder in eventShapeCoex.cxx
DECLARE_SOA_COLUMN(OccTracks, occTracks, int);        //! track occupancy in time range
DECLARE_SOA_COLUMN(OccFT0C, occFT0C, float);          //! FT0C occupancy in time range
DECLARE_SOA_COLUMN(NumContrib, numContrib, uint16_t); //! number of PV contributors
DECLARE_SOA_COLUMN(CollTimeRes, collTimeRes, float);  //! collision time resolution (ns)
DECLARE_SOA_COLUMN(NPVC, nPVC, uint16_t);             //! PV contributors among the selected tracks
} // namespace evshapecoex

DECLARE_SOA_TABLE(EvShapeCoex, "AOD", "EVSHAPECOEX", //! per-collision forward/backward sub-events, reconstructed level
                  o2::soa::Index<>,
                  evshapecoex::CollisionId, evshapecoex::RunNumber, evshapecoex::GlobalBC, evshapecoex::PosZ,
                  evshapecoex::MultFT0A, evshapecoex::MultFT0C,
                  evshapecoex::NF, evshapecoex::NB, evshapecoex::NGap,
                  evshapecoex::AF, evshapecoex::AB, evshapecoex::SumPt2F, evshapecoex::SumPt2B,
                  evshapecoex::QxF, evshapecoex::QyF, evshapecoex::QxB, evshapecoex::QyB,
                  evshapecoex::Qx4F, evshapecoex::Qy4F, evshapecoex::Qx4B, evshapecoex::Qy4B,
                  evshapecoex::NPlus, evshapecoex::NMinus,
                  evshapecoex::QaBits, evshapecoex::OccTracks, evshapecoex::OccFT0C,
                  evshapecoex::NumContrib, evshapecoex::CollTimeRes, evshapecoex::NPVC);
using EvShapeCoexRow = EvShapeCoex::iterator;

// Generator level: charged physical primaries, one row per McCollision (no event selection).
namespace evshapecoexgen
{
DECLARE_SOA_COLUMN(McCollisionId, mcCollisionId, int); //! source McCollision global index (DF-local; not an index column)
DECLARE_SOA_COLUMN(PosZ, posZ, float);                 //! generated primary-vertex z (cm)
DECLARE_SOA_COLUMN(NFwdA, nFwdA, uint16_t);            //! charged primaries in the FT0-A acceptance
DECLARE_SOA_COLUMN(NFwdC, nFwdC, uint16_t);            //! charged primaries in the FT0-C acceptance
DECLARE_SOA_COLUMN(NF, nF, uint16_t);                  //! particle count, forward sub-event
DECLARE_SOA_COLUMN(NB, nB, uint16_t);                  //! particle count, backward sub-event
DECLARE_SOA_COLUMN(NGap, nGap, uint16_t);              //! particle count, central gap
DECLARE_SOA_COLUMN(AF, aF, float);                     //! mean pT, forward (GeV/c); NaN if empty
DECLARE_SOA_COLUMN(AB, aB, float);                     //! mean pT, backward (GeV/c); NaN if empty
DECLARE_SOA_COLUMN(SumPt2F, sumPt2F, float);           //! sum of pT^2, forward (GeV^2/c^2)
DECLARE_SOA_COLUMN(SumPt2B, sumPt2B, float);           //! sum of pT^2, backward (GeV^2/c^2)
DECLARE_SOA_COLUMN(QxF, qxF, float);                   //! sum of cos(2 phi), forward
DECLARE_SOA_COLUMN(QyF, qyF, float);                   //! sum of sin(2 phi), forward
DECLARE_SOA_COLUMN(QxB, qxB, float);                   //! sum of cos(2 phi), backward
DECLARE_SOA_COLUMN(QyB, qyB, float);                   //! sum of sin(2 phi), backward
DECLARE_SOA_COLUMN(Qx4F, qx4F, float);                 //! sum of cos(4 phi), forward
DECLARE_SOA_COLUMN(Qy4F, qy4F, float);                 //! sum of sin(4 phi), forward
DECLARE_SOA_COLUMN(Qx4B, qx4B, float);                 //! sum of cos(4 phi), backward
DECLARE_SOA_COLUMN(Qy4B, qy4B, float);                 //! sum of sin(4 phi), backward
DECLARE_SOA_COLUMN(NPlus, nPlus, uint16_t);            //! positive particles, forward + backward
DECLARE_SOA_COLUMN(NMinus, nMinus, uint16_t);          //! negative particles, forward + backward
} // namespace evshapecoexgen

DECLARE_SOA_TABLE(EvShapeCoexGen, "AOD", "EVSHAPECOEXGEN", //! per-McCollision forward/backward sub-events, generator level
                  o2::soa::Index<>,
                  evshapecoexgen::McCollisionId, evshapecoexgen::PosZ,
                  evshapecoexgen::NFwdA, evshapecoexgen::NFwdC,
                  evshapecoexgen::NF, evshapecoexgen::NB, evshapecoexgen::NGap,
                  evshapecoexgen::AF, evshapecoexgen::AB, evshapecoexgen::SumPt2F, evshapecoexgen::SumPt2B,
                  evshapecoexgen::QxF, evshapecoexgen::QyF, evshapecoexgen::QxB, evshapecoexgen::QyB,
                  evshapecoexgen::Qx4F, evshapecoexgen::Qy4F, evshapecoexgen::Qx4B, evshapecoexgen::Qy4B,
                  evshapecoexgen::NPlus, evshapecoexgen::NMinus);
using EvShapeCoexGenRow = EvShapeCoexGen::iterator;

namespace evshapecoexmclabel
{
DECLARE_SOA_INDEX_COLUMN_FULL(EvShapeCoexGen, evShapeCoexGen, int, EvShapeCoexGen, ""); //! matching EvShapeCoexGen row; negative if unlabelled
} // namespace evshapecoexmclabel

DECLARE_SOA_TABLE(EvShapeCoexMcLabels, "AOD", "EVSHAPECOEXLBL", //! MC label, joinable with EvShapeCoex
                  o2::soa::Index<>, evshapecoexmclabel::EvShapeCoexGenId);
using EvShapeCoexMcLabel = EvShapeCoexMcLabels::iterator;

} // namespace o2::aod

#endif // PWGLF_DATAMODEL_LFEVENTTOPOLOGYTABLES_H_
