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

#ifndef PWGCF_FEMTOWORLD_DATAMODEL_FEMTOWORLDDERIVED_H_
#define PWGCF_FEMTOWORLD_DATAMODEL_FEMTOWORLDDERIVED_H_

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/DataTypes.h"
#include "Framework/Expressions.h"
#include "MathUtils/Utils.h"

#include <cmath>

namespace o2::aod
{
/// FemtoWorldCollision
namespace femtoworldcollision
{
DECLARE_SOA_COLUMN(MultV0M, multV0M, float);       //! V0M multiplicity
DECLARE_SOA_COLUMN(Sphericity, sphericity, float); //! Sphericity of the event
DECLARE_SOA_COLUMN(MagField, magField, float);     //! Sphericity of the event
DECLARE_SOA_COLUMN(RunCentrality, runCent, float); //! Sphericity of the event

} // namespace femtoworldcollision

DECLARE_SOA_TABLE(FemtoWorldCollisions, "AOD", "FEMTOWORLDCOLS",
                  o2::soa::Index<>,
                  o2::aod::collision::PosZ,
                  femtoworldcollision::MultV0M,
                  femtoworldcollision::Sphericity,
                  femtoworldcollision::MagField,
                  femtoworldcollision::RunCentrality);
using FemtoWorldCollision = FemtoWorldCollisions::iterator;

/// FemtoWorldTrack
namespace femtoworldparticle
{
/// Distinuishes the different particle types
enum ParticleType {
  kTrack,           //! Track
  kV0,              //! V0
  kV0Child,         //! Child track of a V0
  kPhi,             // PhiCandidate type
  kPhiChild,        // PhiCandidate type
  kD0D0bar,         //! D0/D0bar meson
  kD0D0barChild,    //! Child track of a D0/D0bar meson
  kCascade,         //! Cascade
  kCascadeBachelor, //! Bachelor track of a cascade
  kNParticleTypes   //! Number of particle types

};

static constexpr std::string_view ParticleTypeName[kNParticleTypes] = {"Tracks", "V0", "V0Child", "Phi", "PhiChild", "D0D0bar", "D0D0barChild", "Cascade", "CascadeBachelor"}; //! Naming of the different particle types

using cutContainerType = uint32_t; //! Definition of the data type for the bit-wise container for the different selection criteria

enum TrackType {
  kNoChild,    //! Not a V0 child
  kPosChild,   //! Positive V0 child
  kNegChild,   //! Negative V0 child
  kNTrackTypes //! Number of child types
};

static constexpr std::string_view TrackTypeName[kNTrackTypes] = {"Trk", "Pos", "Neg"}; //! Naming of the different particle types

DECLARE_SOA_INDEX_COLUMN(FemtoWorldCollision, femtoWorldCollision);
DECLARE_SOA_COLUMN(Pt, pt, float);                    //! p_T (GeV/c)
DECLARE_SOA_COLUMN(Eta, eta, float);                  //! Eta
DECLARE_SOA_COLUMN(Phi, phi, float);                  //! Phi
DECLARE_SOA_COLUMN(Mass, mass, float);                //! Mass of the particle
DECLARE_SOA_COLUMN(MassD0, massD0, float);            //! Mass of the D0 meson
DECLARE_SOA_COLUMN(MassD0bar, massD0bar, float);      //! Mass of the D0bar meson
DECLARE_SOA_COLUMN(FlagD0, flagD0, uint8_t);          //! HF-Flag of the D0 meson
DECLARE_SOA_COLUMN(FlagD0bar, flagD0bar, uint8_t);    //! HF-Flag of the D0bar meson
DECLARE_SOA_COLUMN(PartType, partType, uint8_t);      //! Type of the particle, according to femtoworldparticle::ParticleType
DECLARE_SOA_COLUMN(Cut, cut, cutContainerType);       //! Bit-wise container for the different selection criteria
DECLARE_SOA_COLUMN(PIDCut, pidcut, cutContainerType); //! Bit-wise container for the different PID selection criteria \todo since bit-masking cannot be done yet with filters we use a second field for the PID
DECLARE_SOA_COLUMN(TempFitVar, tempFitVar, float);    //! Observable for the template fitting (Track: DCA_xy, V0: CPA)
DECLARE_SOA_COLUMN(Indices, indices, int[2]);         //! Field for the track indices to remove auto-correlations
DECLARE_SOA_COLUMN(MLambda, mLambda, float);          //! The invariant mass of V0 candidate, assuming lambda
DECLARE_SOA_COLUMN(MAntiLambda, mAntiLambda, float);  //! The invariant mass of V0 candidate, assuming antilambda

DECLARE_SOA_DYNAMIC_COLUMN(Theta, theta, //! Compute the theta of the track
                           [](float eta) -> float {
                             return 2.f * std::atan(std::exp(-eta));
                           });
DECLARE_SOA_DYNAMIC_COLUMN(Px, px, //! Compute the momentum in x in GeV/c
                           [](float pt, float phi) -> float {
                             return pt * std::sin(phi);
                           });
DECLARE_SOA_DYNAMIC_COLUMN(Py, py, //! Compute the momentum in y in GeV/c
                           [](float pt, float phi) -> float {
                             return pt * std::cos(phi);
                           });
DECLARE_SOA_DYNAMIC_COLUMN(Pz, pz, //! Compute the momentum in z in GeV/c
                           [](float pt, float eta) -> float {
                             return pt * std::sinh(eta);
                           });
DECLARE_SOA_COLUMN(Sign, sign, int8_t);                                                  //! Sign of the track charge
DECLARE_SOA_COLUMN(Beta, beta, float);                                                   //! TOF Beta signal
DECLARE_SOA_COLUMN(ITSChi2NCl, itsChi2NCl, float);                                       //! Chi2 / cluster for the ITS track segment
DECLARE_SOA_COLUMN(TPCChi2NCl, tpcChi2NCl, float);                                       //! Chi2 / cluster for the TPC track segment
DECLARE_SOA_COLUMN(TPCNClsFound, tpcNClsFound, uint8_t);                                 //! Number of TPC clusters
DECLARE_SOA_COLUMN(TPCNClsCrossedRows, tpcNClsCrossedRows, uint8_t);                     //! Number of TPC crossed rows
DECLARE_SOA_COLUMN(TPCNSigmaKaon, tpcNSigmaKaon, float);                                 //! TPCNSigmaKaon
DECLARE_SOA_COLUMN(TOFNSigmaKaon, tofNSigmaKaon, float);                                 //! TOFNSigmaKaon
DECLARE_SOA_COLUMN(ITSNCls, itsNCls, uint8_t);                                           //! Number of ITS clusters
DECLARE_SOA_COLUMN(ITSNClsInnerBarrel, itsNClsInnerBarrel, uint8_t);                     //! Number of ITS clusters in the inner barrel                             //! TPC signal
DECLARE_SOA_DYNAMIC_COLUMN(TPCCrossedRowsOverFindableCls, tpcCrossedRowsOverFindableCls, //! Compute the number of crossed rows over findable TPC clusters
                           [](uint8_t tpcNClsFindable, uint8_t tpcNClsCrossedRows) -> float {
                             return (float)tpcNClsCrossedRows / (float)tpcNClsFindable;
                           });
// DECLARE_SOA_COLUMN(TPCNSigmaStoreKa, tpcNSigmaStoreKa, binning::binned_t);
DECLARE_SOA_COLUMN(DaughDCA, daughDCA, float);       //! DCA between daughters
DECLARE_SOA_COLUMN(DecayLength, decayLength, float); //! Decay length of the D0/D0bar meson
DECLARE_SOA_COLUMN(TransRadius, transRadius, float); //! Transverse radius of the decay vertex
DECLARE_SOA_COLUMN(DecayVtxX, decayVtxX, float);     //! X position of the decay vertex
DECLARE_SOA_COLUMN(DecayVtxY, decayVtxY, float);     //! Y position of the decay vertex
DECLARE_SOA_COLUMN(DecayVtxZ, decayVtxZ, float);     //! Z position of the decay vertex
DECLARE_SOA_COLUMN(MKaon, mKaon, float);             //! The invariant mass of V0 candidate, assuming kaon

} // namespace femtoworldparticle

DECLARE_SOA_TABLE(FemtoWorldParticles, "AOD", "FEMTOWORLDPARTS",
                  o2::soa::Index<>,
                  femtoworldparticle::FemtoWorldCollisionId,
                  femtoworldparticle::Pt,
                  femtoworldparticle::Eta,
                  femtoworldparticle::Phi,
                  track::P,
                  femtoworldparticle::Mass,
                  femtoworldparticle::MassD0,
                  femtoworldparticle::MassD0bar,
                  femtoworldparticle::FlagD0,
                  femtoworldparticle::FlagD0bar,
                  femtoworldparticle::PartType,
                  femtoworldparticle::Cut,
                  femtoworldparticle::PIDCut,
                  femtoworldparticle::TempFitVar,
                  femtoworldparticle::Indices,
                  femtoworldparticle::MLambda,
                  femtoworldparticle::MAntiLambda,
                  femtoworldparticle::Theta<femtoworldparticle::Eta>,
                  femtoworldparticle::Px<femtoworldparticle::Pt, femtoworldparticle::Phi>,
                  femtoworldparticle::Py<femtoworldparticle::Pt, femtoworldparticle::Phi>,
                  femtoworldparticle::Pz<femtoworldparticle::Pt, femtoworldparticle::Eta>,
                  femtoworldparticle::Sign,
                  pidtofbeta::Beta,
                  femtoworldparticle::ITSChi2NCl,
                  femtoworldparticle::TPCChi2NCl,
                  femtoworldparticle::TPCNSigmaKaon,
                  femtoworldparticle::TOFNSigmaKaon,
                  femtoworldparticle::TPCNClsFound,
                  track::TPCNClsFindable,
                  femtoworldparticle::TPCNClsCrossedRows,
                  track::TPCNClsShared,
                  track::TPCInnerParam,
                  femtoworldparticle::ITSNCls,
                  femtoworldparticle::ITSNClsInnerBarrel,
                  track::DcaXY,
                  track::DcaZ,
                  track::TPCSignal,
                  pidtpc_tiny::TPCNSigmaStoreEl,
                  pidtpc_tiny::TPCNSigmaStorePi,
                  pidtpc_tiny::TPCNSigmaStoreKa,
                  pidtpc_tiny::TPCNSigmaStorePr,
                  pidtpc_tiny::TPCNSigmaStoreDe,
                  pidtof_tiny::TOFNSigmaStoreEl,
                  pidtof_tiny::TOFNSigmaStorePi,
                  pidtof_tiny::TOFNSigmaStoreKa,
                  pidtof_tiny::TOFNSigmaStorePr,
                  pidtof_tiny::TOFNSigmaStoreDe,
                  femtoworldparticle::DaughDCA,
                  femtoworldparticle::DecayLength,
                  femtoworldparticle::TransRadius,
                  femtoworldparticle::DecayVtxX,
                  femtoworldparticle::DecayVtxY,
                  femtoworldparticle::DecayVtxZ,
                  femtoworldparticle::MKaon,
                  femtoworldparticle::TPCCrossedRowsOverFindableCls<track::TPCNClsFindable, femtoworldparticle::TPCNClsCrossedRows>,
                  pidtpc_tiny::TPCNSigmaEl<pidtpc_tiny::TPCNSigmaStoreEl>,
                  pidtpc_tiny::TPCNSigmaPi<pidtpc_tiny::TPCNSigmaStorePi>,
                  pidtpc_tiny::TPCNSigmaKa<pidtpc_tiny::TPCNSigmaStoreKa>,
                  pidtpc_tiny::TPCNSigmaPr<pidtpc_tiny::TPCNSigmaStorePr>,
                  pidtpc_tiny::TPCNSigmaDe<pidtpc_tiny::TPCNSigmaStoreDe>,
                  pidtof_tiny::TOFNSigmaEl<pidtof_tiny::TOFNSigmaStoreEl>,
                  pidtof_tiny::TOFNSigmaPi<pidtof_tiny::TOFNSigmaStorePi>,
                  pidtof_tiny::TOFNSigmaKa<pidtof_tiny::TOFNSigmaStoreKa>,
                  pidtof_tiny::TOFNSigmaPr<pidtof_tiny::TOFNSigmaStorePr>,
                  pidtof_tiny::TOFNSigmaDe<pidtof_tiny::TOFNSigmaStoreDe>);
using FemtoWorldParticle = FemtoWorldParticles::iterator;

namespace femtohash
{
DECLARE_SOA_COLUMN(FemtoBin, bin, int); //! Hash for the event mixing
} // namespace femtohash
DECLARE_SOA_TABLE(FemtoHashes, "AOD", "HASH", femtohash::FemtoBin);
using FemtoHash = FemtoHashes::iterator;

} // namespace o2::aod
#endif // PWGCF_FEMTOWORLD_DATAMODEL_FEMTOWORLDDERIVED_H_
