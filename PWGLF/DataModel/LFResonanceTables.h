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

/// \file LFResonanceTables.h
/// \brief Definitions of tables of resonance decay candidates
///
/// Inspired by StrangenessTables.h, FemtoDerived.h
///
/// \author Bong-Hwi Lim <bong-hwi.lim@cern.ch>

#ifndef O2_ANALYSIS_LFRESONANCETABLES_H_
#define O2_ANALYSIS_LFRESONANCETABLES_H_

#include <cmath>

#include "Common/DataModel/PIDResponse.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
/// Resonance Collisions
namespace resocollision
{
DECLARE_SOA_COLUMN(MultV0M, multV0M, float);       //! V0M multiplicity
DECLARE_SOA_COLUMN(Sphericity, sphericity, float); //! Sphericity of the event
} // namespace resocollision
DECLARE_SOA_TABLE(ResoCollisions, "AOD", "RESOCOL",
                  o2::soa::Index<>,
                  o2::aod::collision::PosX,
                  o2::aod::collision::PosY,
                  o2::aod::collision::PosZ,
                  resocollision::MultV0M,
                  resocollision::Sphericity,
                  timestamp::Timestamp);
using ResoCollision = ResoCollisions::iterator;

// Resonance Daughters
// inspired from PWGCF/DataModel/FemtoDerived.h
namespace resodaughter
{
/// Distinuishes the different daughter types
enum DaughterType {
  kTrack,         //! Track
  kV0,            //! V0
  kCascade,       //! Cascade
  kNDaughterTypes //! Number of Daughter types
};
static constexpr std::string_view DaughterTypeName[kNDaughterTypes] = {"Tracks", "V0", "Cascade"}; //! Naming of the different particle types

enum PDGtype {
  kPion = BIT(0),
  kKaon = BIT(1),
  kProton = BIT(2),
  kHasTOF = BIT(6) // Save hasTOF info for TOF selection
};

#define requireTPCPIDCutInFilter(mask) ((aod::resodaughter::tpcPIDselectionFlag & (uint8_t)aod::resodaughter::mask) == (uint8_t)aod::resodaughter::mask)
#define requireTOFPIDCutInFilter(mask) (((aod::resodaughter::tofPIDselectionFlag & (uint8_t)aod::resodaughter::mask) == (uint8_t)aod::resodaughter::mask) && ((aod::resodaughter::tofPIDselectionFlag & (uint8_t)aod::resodaughter::kHasTOF) == (uint8_t)aod::resodaughter::kHasTOF))
#define requireTPCPIDPionCutInFilter() requireTPCPIDCutInFilter(PDGtype::kPion)
#define requireTPCPIDKaonCutInFilter() requireTPCPIDCutInFilter(PDGtype::kKaon)
#define requireTPCPIDProtonCutInFilter() requireTPCPIDCutInFilter(PDGtype::kProton)
#define requireTOFPIDPionCutInFilter() requireTOFPIDCutInFilter(PDGtype::kPion)
#define requireTOFPIDKaonCutInFilter() requireTOFPIDCutInFilter(PDGtype::kKaon)
#define requireTOFPIDProtonCutInFilter() requireTOFPIDCutInFilter(PDGtype::kProton)

DECLARE_SOA_INDEX_COLUMN(ResoCollision, resoCollision);
DECLARE_SOA_COLUMN(Pt, pt, float);                                     //! p_T (GeV/c)
DECLARE_SOA_COLUMN(Px, px, float);                                     //! p_x (GeV/c)
DECLARE_SOA_COLUMN(Py, py, float);                                     //! p_y (GeV/c)
DECLARE_SOA_COLUMN(Pz, pz, float);                                     //! p_z (GeV/c)
DECLARE_SOA_COLUMN(Eta, eta, float);                                   //! Eta
DECLARE_SOA_COLUMN(Phi, phi, float);                                   //! Phi
DECLARE_SOA_COLUMN(PartType, partType, uint8_t);                       //! Type of the particle, according to resodaughter::ParticleType
DECLARE_SOA_COLUMN(TempFitVar, tempFitVar, float);                     //! Observable for the template fitting (Track: DCA_xy, V0: CPA)
DECLARE_SOA_COLUMN(Indices, indices, int[2]);                          //! Field for the track indices to remove auto-correlations
DECLARE_SOA_COLUMN(Sign, sign, int8_t);                                //! Sign of the track charge
DECLARE_SOA_COLUMN(TPCNClsCrossedRows, tpcNClsCrossedRows, uint8_t);   //! Number of TPC crossed rows
DECLARE_SOA_COLUMN(TPCPIDselectionFlag, tpcPIDselectionFlag, uint8_t); //! TPC PID selection
DECLARE_SOA_COLUMN(TOFPIDselectionFlag, tofPIDselectionFlag, uint8_t); //! TOF PID selection
DECLARE_SOA_COLUMN(DaughDCA, daughDCA, float);                         //! DCA between daughters
DECLARE_SOA_COLUMN(MLambda, mLambda, float);                           //! The invariant mass of V0 candidate, assuming lambda
DECLARE_SOA_COLUMN(MAntiLambda, mAntiLambda, float);                   //! The invariant mass of V0 candidate, assuming antilambda
DECLARE_SOA_COLUMN(TransRadius, transRadius, float);                   //! Transverse radius of the decay vertex
DECLARE_SOA_COLUMN(DecayVtxX, decayVtxX, float);                       //! X position of the decay vertex
DECLARE_SOA_COLUMN(DecayVtxY, decayVtxY, float);                       //! Y position of the decay vertex
DECLARE_SOA_COLUMN(DecayVtxZ, decayVtxZ, float);                       //! Z position of the decay vertex
// For MC
DECLARE_SOA_COLUMN(IsPhysicalPrimary, isPhysicalPrimary, bool);
DECLARE_SOA_COLUMN(ProducedByGenerator, producedByGenerator, bool);
} // namespace resodaughter
DECLARE_SOA_TABLE(ResoDaughters, "AOD", "RESODAUGHTERS",
                  o2::soa::Index<>,
                  resodaughter::ResoCollisionId,
                  resodaughter::Pt,
                  resodaughter::Px,
                  resodaughter::Py,
                  resodaughter::Pz,
                  resodaughter::Eta,
                  resodaughter::Phi,
                  resodaughter::PartType,
                  resodaughter::TempFitVar,
                  resodaughter::Indices,
                  resodaughter::Sign,
                  resodaughter::TPCNClsCrossedRows,
                  o2::aod::track::DcaXY,
                  o2::aod::track::DcaZ,
                  o2::aod::track::X,
                  o2::aod::track::Alpha,
                  resodaughter::TPCPIDselectionFlag,
                  resodaughter::TOFPIDselectionFlag,
                  o2::aod::pidtpc::TPCNSigmaPi,
                  o2::aod::pidtpc::TPCNSigmaKa,
                  o2::aod::pidtpc::TPCNSigmaPr,
                  o2::aod::pidtof::TOFNSigmaPi,
                  o2::aod::pidtof::TOFNSigmaKa,
                  o2::aod::pidtof::TOFNSigmaPr,
                  resodaughter::DaughDCA,
                  resodaughter::MLambda,
                  resodaughter::MAntiLambda,
                  resodaughter::TransRadius,
                  resodaughter::DecayVtxX,
                  resodaughter::DecayVtxY,
                  resodaughter::DecayVtxZ);
using ResoDaughter = ResoDaughters::iterator;

DECLARE_SOA_TABLE(ResoDaughtersMC, "AOD", "RESODAUGHTERSMC",
                  mcparticle::PdgCode,
                  mcparticle::MothersIds,
                  mcparticle::DaughtersIdSlice,
                  resodaughter::IsPhysicalPrimary,
                  resodaughter::ProducedByGenerator);
using Reso2TracksExt = soa::Join<aod::FullTracks, aod::TracksDCA>; // without Extra
using Reso2TracksMC = soa::Join<aod::FullTracks, McTrackLabels>;
using Reso2TracksPID = soa::Join<aod::FullTracks, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr, aod::pidTOFPi, aod::pidTOFKa, aod::pidTOFPr>;
using Reso2TracksPIDExt = soa::Join<Reso2TracksPID, aod::TracksDCA>; // Without Extra

} // namespace o2::aod
#endif // O2_ANALYSIS_LFRESONANCETABLES_H_