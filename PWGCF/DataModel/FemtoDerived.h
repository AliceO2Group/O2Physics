// Copyright 2019-2022 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef PWGCF_DATAMODEL_FEMTODERIVED_H_
#define PWGCF_DATAMODEL_FEMTODERIVED_H_

#include <cmath>
#include "Framework/ASoA.h"
#include "MathUtils/Utils.h"
#include "Framework/DataTypes.h"
#include "Common/DataModel/Multiplicity.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/Expressions.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/PIDResponse.h"

namespace o2::aod
{
/// FemtoDreamCollision
namespace femtodreamcollision
{
DECLARE_SOA_COLUMN(MultV0M, multV0M, float);       //! V0M multiplicity
DECLARE_SOA_COLUMN(MultNtr, multNtr, int);         //! multiplicity of charged tracks as defined in the producer
DECLARE_SOA_COLUMN(Sphericity, sphericity, float); //! Sphericity of the event
DECLARE_SOA_COLUMN(MagField, magField, float);     //! Magnetic field of the event

} // namespace femtodreamcollision

DECLARE_SOA_TABLE(FDCollisions, "AOD", "FDCOLLISION",
                  o2::soa::Index<>,
                  o2::aod::collision::PosZ,
                  femtodreamcollision::MultV0M,
                  femtodreamcollision::MultNtr,
                  femtodreamcollision::Sphericity,
                  femtodreamcollision::MagField);
using FDCollision = FDCollisions::iterator;

/// FemtoDreamTrack
namespace femtodreamparticle
{
/// Distinuishes the different particle types
enum ParticleType {
  kTrack,           //! Track
  kV0,              //! V0
  kV0Child,         //! Child track of a V0
  kCascade,         //! Cascade
  kCascadeBachelor, //! Bachelor track of a cascade
  kNParticleTypes   //! Number of particle types
};

enum MomentumType {
  kPt,    //! transverse momentum
  kPreco, //! reconstructed/propagated momentum at the vertex
  kPtpc   //! momentum at the inner wall of the TPC (useful for PID plots)
};

static constexpr std::string_view ParticleTypeName[kNParticleTypes] = {"Tracks", "V0", "V0Child", "Cascade", "CascadeBachelor"}; //! Naming of the different particle types
static constexpr std::string_view TempFitVarName[kNParticleTypes] = {"/hDCAxy", "/hCPA", "/hDCAxy", "/hCPA", "/hDCAxy"};

using cutContainerType = uint32_t; //! Definition of the data type for the bit-wise container for the different selection criteria

enum TrackType {
  kNoChild,    //! Not a V0 child
  kPosChild,   //! Positive V0 child
  kNegChild,   //! Negative V0 child
  kNTrackTypes //! Number of child types
};

static constexpr std::string_view TrackTypeName[kNTrackTypes] = {"Trk", "Pos", "Neg"}; //! Naming of the different particle types

DECLARE_SOA_INDEX_COLUMN(FDCollision, fdCollision);
DECLARE_SOA_COLUMN(Pt, pt, float);                       //! p_T (GeV/c)
DECLARE_SOA_COLUMN(Eta, eta, float);                     //! Eta
DECLARE_SOA_COLUMN(Phi, phi, float);                     //! Phi
DECLARE_SOA_COLUMN(PartType, partType, uint8_t);         //! Type of the particle, according to femtodreamparticle::ParticleType
DECLARE_SOA_COLUMN(Cut, cut, cutContainerType);          //! Bit-wise container for the different selection criteria
DECLARE_SOA_COLUMN(PIDCut, pidcut, cutContainerType);    //! Bit-wise container for the different PID selection criteria \todo since bit-masking cannot be done yet with filters we use a second field for the PID
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
DECLARE_SOA_DYNAMIC_COLUMN(P, p, //! Compute the overall momentum in GeV/c
                           [](float pt, float eta) -> float {
                             return pt * std::cosh(eta);
                           });
// debug variables
DECLARE_SOA_COLUMN(Sign, sign, int8_t);                                                  //! Sign of the track charge
DECLARE_SOA_COLUMN(TPCNClsFound, tpcNClsFound, uint8_t);                                 //! Number of TPC clusters
DECLARE_SOA_COLUMN(TPCNClsCrossedRows, tpcNClsCrossedRows, uint8_t);                     //! Number of TPC crossed rows
DECLARE_SOA_COLUMN(ITSNCls, itsNCls, uint8_t);                                           //! Number of ITS clusters
DECLARE_SOA_COLUMN(ITSNClsInnerBarrel, itsNClsInnerBarrel, uint8_t);                     //! Number of ITS clusters in the inner barrel                             //! TPC signal
DECLARE_SOA_DYNAMIC_COLUMN(TPCCrossedRowsOverFindableCls, tpcCrossedRowsOverFindableCls, //! Compute the number of crossed rows over findable TPC clusters
                           [](uint8_t tpcNClsFindable, uint8_t tpcNClsCrossedRows) -> float {
                             return (float)tpcNClsCrossedRows / (float)tpcNClsFindable;
                           });
DECLARE_SOA_COLUMN(TPCNSigmaEl, tpcNSigmaEl, float); //! Nsigma separation with the TPC detector for electron
DECLARE_SOA_COLUMN(TPCNSigmaPi, tpcNSigmaPi, float); //! Nsigma separation with the TPC detector for pion
DECLARE_SOA_COLUMN(TPCNSigmaKa, tpcNSigmaKa, float); //! Nsigma separation with the TPC detector for kaon
DECLARE_SOA_COLUMN(TPCNSigmaPr, tpcNSigmaPr, float); //! Nsigma separation with the TPC detector for proton
DECLARE_SOA_COLUMN(TPCNSigmaDe, tpcNSigmaDe, float); //! Nsigma separation with the TPC detector for deuteron
DECLARE_SOA_COLUMN(TOFNSigmaEl, tofNSigmaEl, float); //! Nsigma separation with the TPC detector for electron
DECLARE_SOA_COLUMN(TOFNSigmaPi, tofNSigmaPi, float); //! Nsigma separation with the TPC detector for pion
DECLARE_SOA_COLUMN(TOFNSigmaKa, tofNSigmaKa, float); //! Nsigma separation with the TPC detector for kaon
DECLARE_SOA_COLUMN(TOFNSigmaPr, tofNSigmaPr, float); //! Nsigma separation with the TPC detector for proton
DECLARE_SOA_COLUMN(TOFNSigmaDe, tofNSigmaDe, float); //! Nsigma separation with the TPC detector for deuteron
DECLARE_SOA_COLUMN(DaughDCA, daughDCA, float);       //! DCA between daughters
DECLARE_SOA_COLUMN(TransRadius, transRadius, float); //! Transverse radius of the decay vertex
DECLARE_SOA_COLUMN(DecayVtxX, decayVtxX, float);     //! X position of the decay vertex
DECLARE_SOA_COLUMN(DecayVtxY, decayVtxY, float);     //! Y position of the decay vertex
DECLARE_SOA_COLUMN(DecayVtxZ, decayVtxZ, float);     //! Z position of the decay vertex
DECLARE_SOA_COLUMN(MKaon, mKaon, float);             //! The invariant mass of V0 candidate, assuming kaon
} // namespace femtodreamparticle

DECLARE_SOA_TABLE(FDParticles, "AOD", "FDPARTICLE",
                  o2::soa::Index<>,
                  femtodreamparticle::FDCollisionId,
                  femtodreamparticle::Pt,
                  femtodreamparticle::Eta,
                  femtodreamparticle::Phi,
                  femtodreamparticle::PartType,
                  femtodreamparticle::Cut,
                  femtodreamparticle::PIDCut,
                  femtodreamparticle::TempFitVar,
                  femtodreamparticle::ChildrenIds,
                  femtodreamparticle::MLambda,
                  femtodreamparticle::MAntiLambda,
                  femtodreamparticle::Theta<femtodreamparticle::Eta>,
                  femtodreamparticle::Px<femtodreamparticle::Pt, femtodreamparticle::Phi>,
                  femtodreamparticle::Py<femtodreamparticle::Pt, femtodreamparticle::Phi>,
                  femtodreamparticle::Pz<femtodreamparticle::Pt, femtodreamparticle::Eta>,
                  femtodreamparticle::P<femtodreamparticle::Pt, femtodreamparticle::Eta>);
using FDParticle = FDParticles::iterator;

DECLARE_SOA_TABLE(FDExtParticles, "AOD", "FDEXTPARTICLE",
                  femtodreamparticle::Sign,
                  femtodreamparticle::TPCNClsFound,
                  track::TPCNClsFindable,
                  femtodreamparticle::TPCNClsCrossedRows,
                  track::TPCNClsShared,
                  track::TPCInnerParam,
                  femtodreamparticle::ITSNCls,
                  femtodreamparticle::ITSNClsInnerBarrel,
                  track::DcaXY,
                  track::DcaZ,
                  track::TPCSignal,
                  femtodreamparticle::TPCNSigmaEl,
                  femtodreamparticle::TPCNSigmaPi,
                  femtodreamparticle::TPCNSigmaKa,
                  femtodreamparticle::TPCNSigmaPr,
                  femtodreamparticle::TPCNSigmaDe,
                  femtodreamparticle::TOFNSigmaEl,
                  femtodreamparticle::TOFNSigmaPi,
                  femtodreamparticle::TOFNSigmaKa,
                  femtodreamparticle::TOFNSigmaPr,
                  femtodreamparticle::TOFNSigmaDe,
                  femtodreamparticle::DaughDCA,
                  femtodreamparticle::TransRadius,
                  femtodreamparticle::DecayVtxX,
                  femtodreamparticle::DecayVtxY,
                  femtodreamparticle::DecayVtxZ,
                  femtodreamparticle::MKaon,
                  femtodreamparticle::TPCCrossedRowsOverFindableCls<track::TPCNClsFindable, femtodreamparticle::TPCNClsCrossedRows>)
using FDFullParticle = FDExtParticles::iterator;

/// FemtoDreamTrackMC
namespace femtodreamMCparticle
{
/// Distinuishes the different particle origins
enum ParticleOriginMCTruth {
  kPrimary,                    //! Primary track or V0
  kSecondary,                  //! Particle from a decay
  kMaterial,                   //! Particle from a material
  kNotPrimary,                 //! Not primary particles (kept for compatibility reasons with the FullProducer task. will be removed, since we look at "non primaries" more differentially now)
  kFake,                       //! particle, that has NOT the PDG code of the current analysed particle
  kSecondaryDaughterLambda,    //! Daughter from a Lambda decay
  kSecondaryDaughterSigmaplus, //! Daughter from a Sigma^plus decay
  kNOriginMCTruthTypes
};

//! Naming of the different OriginMCTruth types
static constexpr std::string_view ParticleOriginMCTruthName[kNOriginMCTruthTypes] = {
  "_Primary",
  "_Secondary",
  "_Material",
  "_NotPrimary",
  "_Fake",
  "_SecondaryDaughterLambda",
  "_SecondaryDaughterSigmaPlus"};

/// Distinguished between reconstructed and truth
enum MCType {
  kRecon, //! Reconstructed in case of MC and used as default in case of data
  kTruth, //! MC truth
  kNMCTypes
};

static constexpr std::string_view MCTypeName[kNMCTypes] = {"", "_MC"};

DECLARE_SOA_COLUMN(PartOriginMCTruth, partOriginMCTruth, uint8_t); //! Origin of the particle, according to femtodreamparticle::ParticleOriginMCTruth
DECLARE_SOA_COLUMN(PDGMCTruth, pdgMCTruth, int);                   //! Particle PDG

// debug variables
DECLARE_SOA_COLUMN(MotherPDG, motherPDG, int); //! Checks mother PDG, where mother is the primary particle for that decay chain
} // namespace femtodreamMCparticle

DECLARE_SOA_TABLE(FDMCParticles, "AOD", "FDMCPARTICLE",
                  o2::soa::Index<>,
                  femtodreamMCparticle::PartOriginMCTruth,
                  femtodreamMCparticle::PDGMCTruth,
                  femtodreamparticle::Pt,
                  femtodreamparticle::Eta,
                  femtodreamparticle::Phi);
using FDMCParticle = FDMCParticles::iterator;

DECLARE_SOA_TABLE(FDExtMCParticles, "AOD", "FDEXTMCPARTICLE",
                  femtodreamMCparticle::MotherPDG);
using FDExtMCParticle = FDExtMCParticles::iterator;

namespace mcfdlabel
{
DECLARE_SOA_INDEX_COLUMN(FDMCParticle, fdMCParticle); //! MC particle for femtodreamparticle
} // namespace mcfdlabel
DECLARE_SOA_TABLE(FDMCLabels, "AOD", "FDMCLabel", //! Table joinable to FemtoDreamParticle containing the MC labels
                  mcfdlabel::FDMCParticleId);

/// Hash
namespace hash
{
DECLARE_SOA_COLUMN(Bin, bin, int); //! Hash for the event mixing
} // namespace hash
DECLARE_SOA_TABLE(Hashes, "AOD", "HASH", hash::Bin);
using Hash = Hashes::iterator;

} // namespace o2::aod

#endif // PWGCF_DATAMODEL_FEMTODERIVED_H_
