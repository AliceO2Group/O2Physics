// Copyright 2019-2025 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file FemtoDerived.h
/// \brief Declaration of FemtoFlow tables
/// \author Zuzanna Chochulska, WUT Warsaw & CTU Prague, zchochul@cern.ch

#ifndef PWGCF_FEMTO_DATAMODEL_FEMTODERIVED_H_
#define PWGCF_FEMTO_DATAMODEL_FEMTODERIVED_H_

#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/DataTypes.h"
#include "Framework/Expressions.h"
#include "MathUtils/Utils.h"

#include <cmath>

namespace o2::aod
{
/// FemtoFlowCollision
namespace femtoflowcollision
{
DECLARE_SOA_COLUMN(MultV0M, multV0M, float);                 //! V0M multiplicity
DECLARE_SOA_COLUMN(MultNtr, multNtr, int);                   //! multiplicity of charged tracks as defined in the producer
DECLARE_SOA_COLUMN(Sphericity, sphericity, float);           //! Sphericity of the event
DECLARE_SOA_COLUMN(qnBin, qnbin, int);                       //! Bin of qn-vector of the event
DECLARE_SOA_COLUMN(MagField, magField, float);               //! Magnetic field of the event
DECLARE_SOA_COLUMN(InteractionRate, interactionRate, float); //! Interaction rate
DECLARE_SOA_COLUMN(Occupancy, occupancy, int);               //! TPC occupancy

} // namespace femtoflowcollision

DECLARE_SOA_TABLE(FDCollisions, "AOD", "FDCOLLISION",
                  o2::soa::Index<>,
                  o2::aod::collision::PosZ,
                  femtoflowcollision::MultV0M,
                  femtoflowcollision::MultNtr,
                  femtoflowcollision::Sphericity,
                  femtoflowcollision::qnBin,
                  femtoflowcollision::MagField);
using FDCollision = FDCollisions::iterator;

DECLARE_SOA_TABLE(FDExtCollisions, "AOD", "FDEXTCOLLISION",
                  femtoflowcollision::InteractionRate,
                  femtoflowcollision::Occupancy);
using FDExtCollision = FDExtCollisions::iterator;

/// FemtoFlowTrack
namespace femtoflowparticle
{
/// Distinuishes the different particle types
enum ParticleType {
  kTrack,           //! Track
  kMCTruthTrack,    //! MC Truth Track
  kV0,              //! V0
  kV0Child,         //! Child track of a V0
  kCascade,         //! Cascade
  kCascadeBachelor, //! Bachelor track of a cascade
  kPhi,             //! Phi meson
  kPhiChild,        //! Child track of a Phi meson
  kD0,              //! D0/D0bar meson
  kD0Child,         //! Child track of a D0/D0bar meson
  kNParticleTypes   //! Number of particle types
};

static constexpr std::string_view ParticleTypeName[kNParticleTypes] = {"Tracks", "MCTruthTracks", "V0", "V0Child", "Cascade", "CascadeBachelor", "Phi", "PhiChild", "D0", "D0Child"}; //! Naming of the different particle types
static constexpr std::string_view TempFitVarName[kNParticleTypes] = {"/hDCAxy", "/hPDGvspT", "/hCPA", "/hDCAxy", "/hCPA", "/hDCAxy", "/hInvMass", "/hDCAxy", "/hInvMass", "/hDCAxy"};

using CutContainerType = uint32_t; //! Definition of the data type for the bit-wise container for the different selection criteria

enum TrackType {
  kNoChild,    //! Not a V0 child
  kPosChild,   //! Positive V0 child
  kNegChild,   //! Negative V0 child
  kBachelor,   //! Cascade bachelor
  kNTrackTypes //! Number of child types
};

static constexpr std::string_view TrackTypeName[kNTrackTypes] = {"Trk", "Pos", "Neg", "Bach"}; //! Naming of the different particle types

DECLARE_SOA_INDEX_COLUMN(FDCollision, fdCollision);
DECLARE_SOA_COLUMN(Pt, pt, float);                       //! p_T (GeV/c)
DECLARE_SOA_COLUMN(Eta, eta, float);                     //! Eta
DECLARE_SOA_COLUMN(Phi, phi, float);                     //! Phi
DECLARE_SOA_COLUMN(PartType, partType, uint8_t);         //! Type of the particle, according to femtoflowparticle::ParticleType
DECLARE_SOA_COLUMN(Cut, cut, CutContainerType);          //! Bit-wise container for the different selection criteria
DECLARE_SOA_COLUMN(PidCut, pidCut, CutContainerType);    //! Bit-wise container for the different PID selection criteria \todo since bit-masking cannot be done yet with filters we use a second field for the PID
DECLARE_SOA_COLUMN(TempFitVar, tempFitVar, float);       //! Observable for the template fitting (Track: DCA_xy, V0: CPA)
DECLARE_SOA_SELF_ARRAY_INDEX_COLUMN(Children, children); //! Field for the track indices to remove auto-correlations
DECLARE_SOA_COLUMN(MLambda, mLambda, float);             //! The invariant mass of V0 candidate, assuming lambda
DECLARE_SOA_COLUMN(MAntiLambda, mAntiLambda, float);     //! The invariant mass of V0 candidate, assuming antilambda

DECLARE_SOA_DYNAMIC_COLUMN(Theta, theta, //! Compute the theta of the track
                           [](float eta) -> float {
                             return 2.f * std::atan(std::exp(-eta));
                           });
DECLARE_SOA_DYNAMIC_COLUMN(Px, px, //! Compute the momentum in x in GeV/c
                           [](float pt, float phi) -> float {
                             return pt * std::cos(phi);
                           });
DECLARE_SOA_DYNAMIC_COLUMN(Py, py, //! Compute the momentum in y in GeV/c
                           [](float pt, float phi) -> float {
                             return pt * std::sin(phi);
                           });
DECLARE_SOA_DYNAMIC_COLUMN(Pz, pz, //! Compute the momentum in z in GeV/c
                           [](float pt, float eta) -> float {
                             return pt * std::sinh(eta);
                           });
DECLARE_SOA_DYNAMIC_COLUMN(P, p, //! Compute the overall momentum in GeV/c
                           [](float pt, float eta) -> float {
                             return pt * std::cosh(eta);
                           });
// debug variables
DECLARE_SOA_COLUMN(Sign, sign, int8_t);                                                  //! Sign of the track charge
DECLARE_SOA_COLUMN(TpcNClsFound, tpcNClsFound, uint8_t);                                 //! Number of TPC clusters
DECLARE_SOA_COLUMN(TpcNClsCrossedRows, tpcNClsCrossedRows, uint8_t);                     //! Number of TPC crossed rows
DECLARE_SOA_COLUMN(TpcFractionSharedCls, tpcFractionSharedCls, float);                   //! Number of TPC crossed rows
DECLARE_SOA_COLUMN(ItsNCls, itsNCls, uint8_t);                                           //! Number of ITS clusters
DECLARE_SOA_COLUMN(ItsNClsInnerBarrel, itsNClsInnerBarrel, uint8_t);                     //! Number of ITS clusters in the inner barrel                             //! TPC signal
DECLARE_SOA_DYNAMIC_COLUMN(TpcCrossedRowsOverFindableCls, tpcCrossedRowsOverFindableCls, //! Compute the number of crossed rows over findable TPC clusters
                           [](uint8_t tpcNClsFindable, uint8_t tpcNClsCrossedRows) -> float {
                             return (float)tpcNClsCrossedRows / (float)tpcNClsFindable;
                           });
DECLARE_SOA_COLUMN(DaughDCA, daughDCA, float);       //! DCA between daughters
DECLARE_SOA_COLUMN(TransRadius, transRadius, float); //! Transverse radius of the decay vertex
DECLARE_SOA_COLUMN(DecayVtxX, decayVtxX, float);     //! X position of the decay vertex
DECLARE_SOA_COLUMN(DecayVtxY, decayVtxY, float);     //! Y position of the decay vertex
DECLARE_SOA_COLUMN(DecayVtxZ, decayVtxZ, float);     //! Z position of the decay vertex
DECLARE_SOA_COLUMN(MKaon, mKaon, float);             //! The invariant mass of V0 candidate, assuming kaon

} // namespace femtoflowparticle

DECLARE_SOA_TABLE(FDParticles, "AOD", "FDPARTICLE",
                  o2::soa::Index<>,
                  femtoflowparticle::FDCollisionId,
                  femtoflowparticle::Pt,
                  femtoflowparticle::Eta,
                  femtoflowparticle::Phi,
                  femtoflowparticle::PartType,
                  femtoflowparticle::Cut,
                  femtoflowparticle::PidCut,
                  femtoflowparticle::TempFitVar,
                  femtoflowparticle::ChildrenIds,
                  femtoflowparticle::MLambda,
                  femtoflowparticle::MAntiLambda,
                  femtoflowparticle::Theta<femtoflowparticle::Eta>,
                  femtoflowparticle::Px<femtoflowparticle::Pt, femtoflowparticle::Phi>,
                  femtoflowparticle::Py<femtoflowparticle::Pt, femtoflowparticle::Phi>,
                  femtoflowparticle::Pz<femtoflowparticle::Pt, femtoflowparticle::Eta>,
                  femtoflowparticle::P<femtoflowparticle::Pt, femtoflowparticle::Eta>);
using FDParticle = FDParticles::iterator;

/// FemtoFlowCascadeTrack
namespace femtoflowcascparticle
{
DECLARE_SOA_INDEX_COLUMN(FDParticle, fdParticle);
DECLARE_SOA_COLUMN(DcaV0daughters, dcaV0daughters, float);     //! DCA between V0 daughters
DECLARE_SOA_COLUMN(Cpav0, cpav0, float);                       //! V0 cos of pointing angle
DECLARE_SOA_COLUMN(V0radius, v0radius, float);                 //! V0 transverse radius*/
DECLARE_SOA_COLUMN(CpaCasc, cpaCasc, float);                   //! cascade cosinus of pointing angle
DECLARE_SOA_COLUMN(Dcacascdaughters, dcacascdaughters, float); //! DCA between cascade daughters
DECLARE_SOA_COLUMN(Cascradius, cascradius, float);             //! cascade transverse radius
DECLARE_SOA_COLUMN(Dcapostopv, dcapostopv, float);             //! DCA of positive daughter to PV
DECLARE_SOA_COLUMN(Dcanegtopv, dcanegtopv, float);             //! DCA of negative daughter to PV
DECLARE_SOA_COLUMN(Dcabachtopv, dcabachtopv, float);           //! DCA of bachelor track to PV
DECLARE_SOA_COLUMN(Dcav0topv, dcav0topv, float);               //! DCA of V0 to PV

} // namespace femtoflowcascparticle

DECLARE_SOA_TABLE(FDExtParticles, "AOD", "FDEXTPARTICLE",
                  femtoflowparticle::Sign,
                  femtoflowparticle::TpcNClsFound,
                  track::TPCNClsFindable,
                  femtoflowparticle::TpcNClsCrossedRows,
                  track::TPCNClsShared,
                  femtoflowparticle::TpcFractionSharedCls,
                  track::TPCInnerParam,
                  femtoflowparticle::ItsNCls,
                  femtoflowparticle::ItsNClsInnerBarrel,
                  track::DcaXY,
                  track::DcaZ,
                  track::TPCSignal,
                  pidtofbeta::Beta,
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
                  femtoflowparticle::DaughDCA,
                  femtoflowparticle::TransRadius,
                  femtoflowparticle::DecayVtxX,
                  femtoflowparticle::DecayVtxY,
                  femtoflowparticle::DecayVtxZ,
                  femtoflowparticle::MKaon,
                  femtoflowparticle::TpcCrossedRowsOverFindableCls<track::TPCNClsFindable, femtoflowparticle::TpcNClsCrossedRows>,
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
using FDFullParticle = FDExtParticles::iterator;

DECLARE_SOA_TABLE(FDCascParticles, "AOD", "FDCASCPARTICLE",
                  o2::soa::Index<>,
                  femtoflowparticle::FDCollisionId,
                  femtoflowcascparticle::FDParticleId,
                  femtoflowparticle::Theta<femtoflowparticle::Eta>,
                  femtoflowparticle::Px<femtoflowparticle::Pt, femtoflowparticle::Phi>,
                  femtoflowparticle::Py<femtoflowparticle::Pt, femtoflowparticle::Phi>,
                  femtoflowparticle::Pz<femtoflowparticle::Pt, femtoflowparticle::Eta>,
                  femtoflowparticle::P<femtoflowparticle::Pt, femtoflowparticle::Eta>,
                  femtoflowcascparticle::DcaV0daughters,
                  femtoflowcascparticle::Cpav0,
                  femtoflowcascparticle::V0radius,
                  femtoflowcascparticle::CpaCasc,
                  femtoflowcascparticle::Dcacascdaughters,
                  femtoflowcascparticle::Cascradius,
                  femtoflowcascparticle::Dcapostopv,
                  femtoflowcascparticle::Dcanegtopv,
                  femtoflowcascparticle::Dcabachtopv,
                  femtoflowcascparticle::Dcav0topv);
using FDCascParticle = FDCascParticles::iterator;

/// FemtoFlowTrackMC
namespace femtoflow_mc_particle
{
/// Distinuishes the different particle origins
enum ParticleOriginMCTruth {
  kPrimary,           //! Primary track or V0
  kDaughter,          //! Particle from a decay
  kMaterial,          //! Particle from a material
  kNotPrimary,        //! Not primary particles (kept for compatibility reasons with the FullProducer task. will be removed, since we look at "non primaries" more differentially now)
  kFake,              //! Particle, that has NOT the PDG code of the current analysed particle
  kDaughterLambda,    //! Daughter from a Lambda decay
  kDaughterSigmaplus, //! Daughter from a Sigma^plus decay
  kPrompt,            //! Origin for D0/D0bar mesons
  kNonPrompt,         //! Origin for D0/D0bar mesons
  kNOriginMCTruthTypes,
  kElse,
  kWrongCollision //! Origin for the wrong collision
};

//! Naming of the different OriginMCTruth types
static constexpr std::string_view ParticleOriginMCTruthName[kNOriginMCTruthTypes] = {
  "_Primary",
  "_Daughter",
  "_Material",
  "_NotPrimary",
  "_Fake",
  "_DaughterLambda",
  "DaughterSigmaPlus",
  "_Prompt",
  "_NonPrompt"};

/// Distinguished between reconstructed and truth
enum MCType {
  kRecon, //! Reconstructed in case of MC and used as default in case of data
  kTruth, //! MC truth
  kNMCTypes
};

static constexpr std::string_view MCTypeName[kNMCTypes] = {"", "_MC"};

DECLARE_SOA_COLUMN(PartOriginMCTruth, partOriginMCTruth, uint8_t); //! Origin of the particle, according to femtoflowparticle::ParticleOriginMCTruth
DECLARE_SOA_COLUMN(PdgMCTruth, pdgMCTruth, int);                   //! Particle PDG

// debug variables
DECLARE_SOA_COLUMN(MotherPDG, motherPDG, int); //! Checks mother PDG, where mother is the primary particle for that decay chain
} // namespace femtoflow_mc_particle

DECLARE_SOA_TABLE(FDMCParticles, "AOD", "FDMCPARTICLE",
                  o2::soa::Index<>,
                  femtoflow_mc_particle::PartOriginMCTruth,
                  femtoflow_mc_particle::PdgMCTruth,
                  femtoflowparticle::Pt,
                  femtoflowparticle::Eta,
                  femtoflowparticle::Phi);
using FDMCParticle = FDMCParticles::iterator;

DECLARE_SOA_TABLE(FDExtMCParticles, "AOD", "FDEXTMCPARTICLE",
                  femtoflow_mc_particle::MotherPDG);
using FDExtMCParticle = FDExtMCParticles::iterator;

namespace mcfdlabel
{
DECLARE_SOA_INDEX_COLUMN(FDMCParticle, fdMCParticle); //! MC particle for femtoflowparticle
} // namespace mcfdlabel
DECLARE_SOA_TABLE(FDMCLabels, "AOD", "FDMCLabel", //! Table joinable to FemtoFlowParticle containing the MC labels
                  mcfdlabel::FDMCParticleId);

/// Hash
namespace hash
{
DECLARE_SOA_COLUMN(Bin, bin, int); //! Hash for the event mixing
} // namespace hash
DECLARE_SOA_TABLE(MixingHashes, "AOD", "HASH", hash::Bin);
using MixingHash = MixingHashes::iterator;

} // namespace o2::aod

#endif // PWGCF_FEMTO_DATAMODEL_FEMTODERIVED_H_
