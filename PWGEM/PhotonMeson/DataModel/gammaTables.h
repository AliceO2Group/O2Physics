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

#include "Framework/AnalysisDataModel.h"
#include "Common/Core/PID/PIDResponse.h"
#include "Common/DataModel/StrangenessTables.h"
#include "Common/DataModel/TrackSelectionTables.h"

// todo: declare more columns in this file dynamic or expression, atm I save a lot of redundant information
namespace o2::aod
{
namespace gammatrackreco
{
DECLARE_SOA_COLUMN(Eta, eta, float);                                                     //! Pseudorapidity
DECLARE_SOA_COLUMN(P, p, float);                                                         //! Total momentum in GeV/c
DECLARE_SOA_COLUMN(Phi, phi, float);                                                     //! Azimuthal angle
DECLARE_SOA_COLUMN(Pt, pt, float);                                                       //! Transversal momentum in GeV/c
DECLARE_SOA_COLUMN(PositivelyCharged, positivelyCharged, bool);                          //! True for positively charged track
DECLARE_SOA_COLUMN(TpcCrossedRowsOverFindableCls, tpcCrossedRowsOverFindableCls, float); //! Ratio  crossed rows over findable clusters
DECLARE_SOA_COLUMN(TpcFoundOverFindableCls, tpcFoundOverFindableCls, float);             //! Ratio of found over findable clusters
DECLARE_SOA_COLUMN(TpcNClsCrossedRows, tpcNClsCrossedRows, float);                       //! Number of crossed TPC Rows
} // namespace gammatrackreco

DECLARE_SOA_TABLE(V0DaughterTracks, "AOD", "V0TRACKS",
                  o2::soa::Index<>,
                  v0data::V0Id,
                  track::DcaXY,
                  gammatrackreco::Eta,
                  gammatrackreco::P,
                  gammatrackreco::Phi,
                  gammatrackreco::Pt,
                  gammatrackreco::PositivelyCharged,
                  gammatrackreco::TpcCrossedRowsOverFindableCls,
                  gammatrackreco::TpcFoundOverFindableCls,
                  gammatrackreco::TpcNClsCrossedRows,
                  pidtpc::TPCNSigmaEl,
                  pidtpc::TPCNSigmaPi,
                  track::TPCSignal);

namespace gammamctrue
{
DECLARE_SOA_COLUMN(Gamma, gamma, int64_t);       //! Used as reference for the daughters
DECLARE_SOA_COLUMN(NDaughters, nDaughters, int); //! Number of daughters
DECLARE_SOA_COLUMN(Eta, eta, float);             //! Pseudorapidity
DECLARE_SOA_COLUMN(Phi, phi, float);             //! Angle phi in rad
DECLARE_SOA_COLUMN(P, p, float);                 //! Absolute momentum in GeV/c
DECLARE_SOA_COLUMN(Pt, pt, float);               //! Transversal momentum in GeV/c
DECLARE_SOA_COLUMN(Y, y, float);                 //! Rapidity

DECLARE_SOA_COLUMN(ConversionX, conversionX, float); //! x of conversion point in cm
DECLARE_SOA_COLUMN(ConversionY, conversionY, float); //! y of conversion point in cm
DECLARE_SOA_COLUMN(ConversionZ, conversionZ, float); //! z of conversion point in cm
DECLARE_SOA_COLUMN(V0Radius, v0Radius, float);       //! 2d radius of conversion point
} // namespace gammamctrue

DECLARE_SOA_TABLE(McGammasTrue, "AOD", "MCGATRUE",
                  o2::soa::Index<>,
                  mcparticle::McCollisionId,
                  gammamctrue::Gamma,
                  v0data::V0Id, // reference to reconstructed v0 (if its a task with reconstucted info)
                  mcparticle::PdgCode, mcparticle::StatusCode, mcparticle::Flags,
                  mcparticle::Px, mcparticle::Py, mcparticle::Pz,
                  mcparticle::Vx, mcparticle::Vy, mcparticle::Vz, mcparticle::Vt,
                  gammamctrue::NDaughters,
                  // todo: make those expression columns in an extended table
                  gammamctrue::Eta, gammamctrue::Phi, gammamctrue::P, gammamctrue::Pt, gammamctrue::Y,
                  gammamctrue::ConversionX, gammamctrue::ConversionY, gammamctrue::ConversionZ,
                  gammamctrue::V0Radius,

                  // Dynamic columns
                  mcparticle::ProducedByGenerator<mcparticle::Flags>,
                  mcparticle::FromBackgroundEvent<mcparticle::Flags>,
                  mcparticle::GetGenStatusCode<mcparticle::Flags, mcparticle::StatusCode>,
                  mcparticle::GetProcess<mcparticle::Flags, mcparticle::StatusCode>,
                  mcparticle::IsPhysicalPrimary<mcparticle::Flags>);

namespace gammadaughtermctrue
{
DECLARE_SOA_INDEX_COLUMN_FULL(Mother0, mother0, int64_t, McGammasTrue, ""); //! index of first mother
DECLARE_SOA_COLUMN(NMothers, nMothers, int);                                //! the number of mothers
} // namespace gammadaughtermctrue

// table to hold mc truth information of daughter particles of MC gammas
DECLARE_SOA_TABLE(McGammaDaughtersTrue, "AOD", "MCGADAUGHTRUE",
                  o2::soa::Index<>,
                  mcparticle::McCollisionId, // SFS maybe drop since there is already a pointer to MCGammas which point to mccollision. But maybe still good to have for correct automatic grouping - not sure about that
                  gammadaughtermctrue::Mother0Id,
                  gammadaughtermctrue::NMothers,

                  mcparticle::PdgCode, mcparticle::StatusCode, mcparticle::Flags,
                  mcparticle::Px, mcparticle::Py, mcparticle::Pz, mcparticle::E,
                  mcparticle::Vx, mcparticle::Vy, mcparticle::Vz, mcparticle::Vt,

                  // dynamic columns
                  mcparticle::ProducedByGenerator<mcparticle::Flags>,
                  mcparticle::FromBackgroundEvent<mcparticle::Flags>,
                  mcparticle::GetGenStatusCode<mcparticle::Flags, mcparticle::StatusCode>,
                  mcparticle::GetProcess<mcparticle::Flags, mcparticle::StatusCode>,
                  mcparticle::IsPhysicalPrimary<mcparticle::Flags>);
} // namespace o2::aod
