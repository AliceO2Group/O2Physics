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
/// \author Nasir Mehdi Malik <nasir.mehdi.malik@cern.ch>
/// \author Minjae Kim <minjae.kim@cern.ch>
///

#ifndef PWGLF_DATAMODEL_LFRESONANCETABLES_H_
#define PWGLF_DATAMODEL_LFRESONANCETABLES_H_

#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/AnalysisDataModel.h"

#include <algorithm>
#include <cmath>

namespace o2::aod
{
/// Resonance Collisions
namespace resocollision
{
enum {
  kECbegin = 0,
  kINEL = 1,
  kINEL10,
  kINELg0,
  kINELg010,
  kTrig,
  kTrig10,
  kTrigINELg0,
  kTrigINELg010,
  kSel8,
  kSel810,
  kSel8INELg0,
  kSel8INELg010,
  kAllCuts,
  kAllCuts10,
  kAllCutsINELg0,
  kAllCutsINELg010,
  kECend,
};
DECLARE_SOA_INDEX_COLUMN_FULL(Collision, collision, int, Collisions, "_Col"); //!
DECLARE_SOA_COLUMN(Cent, cent, float);                                        //! Centrality (Multiplicity) percentile (Default: FT0M)
DECLARE_SOA_COLUMN(Spherocity, spherocity, float);                            //! Spherocity of the event
DECLARE_SOA_COLUMN(EvtPl, evtPl, float);                                      //! Second harmonic event plane
DECLARE_SOA_COLUMN(EvtPlResAB, evtPlResAB, float);                            //! Second harmonic event plane resolution of A-B sub events
DECLARE_SOA_COLUMN(EvtPlResAC, evtPlResAC, float);                            //! Second harmonic event plane resolution of A-C sub events
DECLARE_SOA_COLUMN(EvtPlResBC, evtPlResBC, float);                            //! Second harmonic event plane resolution of B-C sub events
DECLARE_SOA_COLUMN(BMagField, bMagField, float);                              //! Magnetic field
DECLARE_SOA_COLUMN(IsRecINELgt0, isRecINELgt0, bool);                         //! Is reconstructed INEL>0 event
// MC
DECLARE_SOA_COLUMN(IsVtxIn10, isVtxIn10, bool);               //! Vtx10
DECLARE_SOA_COLUMN(IsINELgt0, isINELgt0, bool);               //! INEL>0
DECLARE_SOA_COLUMN(IsTriggerTVX, isTriggerTVX, bool);         //! TriggerTVX
DECLARE_SOA_COLUMN(IsInSel8, isInSel8, bool);                 //! InSel8
DECLARE_SOA_COLUMN(IsInAfterAllCuts, isInAfterAllCuts, bool); //! InAfterAllCuts
DECLARE_SOA_COLUMN(ImpactParameter, impactParameter, float);  //! ImpactParameter
DECLARE_SOA_COLUMN(MCMultiplicity, mcMultiplicity, float);    //! MC Multiplicity

} // namespace resocollision
DECLARE_SOA_TABLE(ResoCollisions, "AOD", "RESOCOLLISION",
                  o2::soa::Index<>,
                  o2::aod::mult::MultNTracksPV,
                  o2::aod::mult::MultNTracksPVeta1,
                  o2::aod::mult::MultNTracksPVetaHalf,
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ,
                  resocollision::Cent,
                  resocollision::BMagField,
                  resocollision::IsRecINELgt0);
using ResoCollision = ResoCollisions::iterator;

DECLARE_SOA_TABLE(ResoCollisionColls, "AOD", "RESOCOLLISIONCOL",
                  resocollision::CollisionId);
using ResoCollisionColl = ResoCollisionColls::iterator;

DECLARE_SOA_TABLE(ResoMCCollisions, "AOD", "RESOMCCOLLISION",
                  o2::soa::Index<>,
                  resocollision::IsVtxIn10,
                  resocollision::IsINELgt0,
                  resocollision::IsTriggerTVX,
                  resocollision::IsInSel8,
                  resocollision::IsInAfterAllCuts,
                  resocollision::ImpactParameter,
                  resocollision::MCMultiplicity);
using ResoMCCollision = ResoMCCollisions::iterator;

DECLARE_SOA_TABLE(ResoSpheroCollisions, "AOD", "RESOSPHEROCOLLISION",
                  o2::soa::Index<>,
                  resocollision::Spherocity);
using ResoSpheroCollision = ResoSpheroCollisions::iterator;

DECLARE_SOA_TABLE(ResoEvtPlCollisions, "AOD", "RESOEVTPLCOLLISION",
                  o2::soa::Index<>,
                  resocollision::EvtPl,
                  resocollision::EvtPlResAB,
                  resocollision::EvtPlResAC,
                  resocollision::EvtPlResBC);
using ResoEvtPlCollision = ResoEvtPlCollisions::iterator;

// For DF mixing study
DECLARE_SOA_TABLE(ResoCollisionDFs, "AOD", "RESOCOLLISIONDF",
                  o2::soa::Index<>,
                  o2::aod::mult::MultNTracksPV,
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ,
                  resocollision::Cent,
                  resocollision::Spherocity,
                  resocollision::EvtPl,
                  resocollision::EvtPlResAB,
                  resocollision::EvtPlResAC,
                  resocollision::EvtPlResBC,
                  resocollision::BMagField,
                  timestamp::Timestamp,
                  evsel::NumTracksInTimeRange);
using ResoCollisionDF = ResoCollisionDFs::iterator;

// Resonance Daughters
// inspired from PWGCF/DataModel/FemtoDerived.h
namespace resodaughter
{
struct ResoTrackFlags {
 public:
  typedef uint8_t flagtype;
  static constexpr flagtype kPassedITSRefit = 1 << 0;
  static constexpr flagtype kPassedTPCRefit = 1 << 1;
  static constexpr flagtype kIsGlobalTrackWoDCA = 1 << 2;
  static constexpr flagtype kIsGlobalTrack = 1 << 3;
  static constexpr flagtype kIsPrimaryTrack = 1 << 4;
  static constexpr flagtype kIsPVContributor = 1 << 5;
  static constexpr flagtype kHasTOF = 1 << 6;
  static constexpr flagtype kSign = 1 << 7;
  /// @brief check if the flag is set
  static bool checkFlag(const flagtype flags, const flagtype mask)
  {
    return (flags & mask) == mask;
  }
};
#define requireTrackFlag(mask) ((o2::aod::resodaughter::trackFlags & o2::aod::resodaughter::mask) == o2::aod::resodaughter::mask)

#define requirePassedITSRefit() requireTrackFlag(ResoTrackFlags::kPassedITSRefit)
#define requirePassedTPCRefit() requireTrackFlag(ResoTrackFlags::kPassedTPCRefit)
#define requireGlobalTrack() requireTrackFlag(ResoTrackFlags::kIsGlobalTrack)
#define requireGlobalTrackWoDCA() requireTrackFlag(ResoTrackFlags::kIsGlobalTrackWoDCA)
#define requirePrimaryTrack() requireTrackFlag(ResoTrackFlags::kIsPrimaryTrack)
#define requirePVContributor() requireTrackFlag(ResoTrackFlags::kIsPVContributor)
#define requireHasTOF() requireTrackFlag(ResoTrackFlags::kHasTOF)
#define requireSign() requireTrackFlag(ResoTrackFlags::kSign)

#define DECLARE_DYN_TRKSEL_COLUMN(name, getter, mask) \
  DECLARE_SOA_DYNAMIC_COLUMN(name, getter, [](ResoTrackFlags::flagtype flags) -> bool { return ResoTrackFlags::checkFlag(flags, mask); });

DECLARE_SOA_INDEX_COLUMN(ResoCollision, resoCollision);
DECLARE_SOA_INDEX_COLUMN(ResoCollisionDF, resoCollisionDF);
DECLARE_SOA_INDEX_COLUMN_FULL(Track, track, int, Tracks, "_Trk");                 //!
DECLARE_SOA_INDEX_COLUMN_FULL(V0, v0, int, V0s, "_V0");                           //!
DECLARE_SOA_INDEX_COLUMN_FULL(Cascade, cascade, int, Cascades, "_Cas");           //!
DECLARE_SOA_COLUMN(Pt, pt, float);                                                //! p_t (GeV/c)
DECLARE_SOA_COLUMN(Px, px, float);                                                //! p_x (GeV/c)
DECLARE_SOA_COLUMN(Py, py, float);                                                //! p_y (GeV/c)
DECLARE_SOA_COLUMN(Pz, pz, float);                                                //! p_z (GeV/c)
DECLARE_SOA_COLUMN(PartType, partType, uint8_t);                                  //! Type of the particle, according to resodaughter::ParticleType
DECLARE_SOA_COLUMN(TempFitVar, tempFitVar, float);                                //! Observable for the template fitting (Track: DCA_xy, V0: CPA)
DECLARE_SOA_COLUMN(Indices, indices, int[2]);                                     //! Field for the track indices to remove auto-correlations
DECLARE_SOA_COLUMN(CascadeIndices, cascadeIndices, int[3]);                       //! Field for the track indices to remove auto-correlations (ordered: positive, negative, bachelor)
DECLARE_SOA_COLUMN(TpcNClsCrossedRows, tpcNClsCrossedRows, uint8_t);              //! Number of TPC crossed rows
DECLARE_SOA_COLUMN(TpcNClsFound, tpcNClsFound, uint8_t);                          //! Number of TPC clusters found
DECLARE_SOA_COLUMN(DcaXY10000, dcaXY10000, int16_t);                              //! DCA_xy x10,000 in int16_t, resolution 10 um
DECLARE_SOA_COLUMN(DcaZ10000, dcaZ10000, int16_t);                                //! DCA_z x10,000 in int16_t, resolution 10 um
DECLARE_SOA_COLUMN(TrackFlags, trackFlags, uint8_t);                              //! Track flags
DECLARE_SOA_COLUMN(TpcNSigmaPi10, tpcNSigmaPi10, int8_t);                         //! TPC PID x10 of the track as Pion
DECLARE_SOA_COLUMN(TpcNSigmaKa10, tpcNSigmaKa10, int8_t);                         //! TPC PID x10 of the track as Kaon
DECLARE_SOA_COLUMN(TpcNSigmaPr10, tpcNSigmaPr10, int8_t);                         //! TPC PID x10 of the track as Proton
DECLARE_SOA_COLUMN(TofNSigmaPi10, tofNSigmaPi10, int8_t);                         //! TOF PID x10 of the track as Pion
DECLARE_SOA_COLUMN(TofNSigmaKa10, tofNSigmaKa10, int8_t);                         //! TOF PID x10 of the track as Kaon
DECLARE_SOA_COLUMN(TofNSigmaPr10, tofNSigmaPr10, int8_t);                         //! TOF PID x10 of the track as Proton
DECLARE_SOA_COLUMN(DaughDCA, daughDCA, float);                                    //! DCA between daughters
DECLARE_SOA_COLUMN(CascDaughDCA, cascDaughDCA, float);                            //! DCA between daughters from cascade
DECLARE_SOA_COLUMN(V0CosPA, v0CosPA, float);                                      //! V0 Cosine of Pointing Angle
DECLARE_SOA_COLUMN(CascCosPA, cascCosPA, float);                                  //! Cascade Cosine of Pointing Angle
DECLARE_SOA_COLUMN(MLambda, mLambda, float);                                      //! The invariant mass of V0 candidate, assuming lambda
DECLARE_SOA_COLUMN(MAntiLambda, mAntiLambda, float);                              //! The invariant mass of V0 candidate, assuming antilambda
DECLARE_SOA_COLUMN(MK0Short, mK0Short, float);                                    //! The invariant mass of V0 candidate, assuming k0s
DECLARE_SOA_COLUMN(MXi, mXi, float);                                              //! The invariant mass of Xi candidate
DECLARE_SOA_COLUMN(TransRadius, transRadius, float);                              //! Transverse radius of the decay vertex
DECLARE_SOA_COLUMN(CascTransRadius, cascTransRadius, float);                      //! Transverse radius of the decay vertex from cascade
DECLARE_SOA_COLUMN(DecayVtxX, decayVtxX, float);                                  //! X position of the decay vertex
DECLARE_SOA_COLUMN(DecayVtxY, decayVtxY, float);                                  //! Y position of the decay vertex
DECLARE_SOA_COLUMN(DecayVtxZ, decayVtxZ, float);                                  //! Z position of the decay vertex
DECLARE_SOA_COLUMN(Alpha, alpha, float);                                          //! Alpha of the decay vertex
DECLARE_SOA_COLUMN(QtArm, qtarm, float);                                          //! Armenteros Qt of the decay vertex
DECLARE_SOA_COLUMN(TpcSignal10, tpcSignal10, int8_t);                             //! TPC signal of the track x10
DECLARE_SOA_COLUMN(DaughterTPCNSigmaPosPi10, daughterTPCNSigmaPosPi10, int8_t);   //! TPC PID x10 of the positive daughter as Pion
DECLARE_SOA_COLUMN(DaughterTPCNSigmaPosKa10, daughterTPCNSigmaPosKa10, int8_t);   //! TPC PID x10 of the positive daughter as Kaon
DECLARE_SOA_COLUMN(DaughterTPCNSigmaPosPr10, daughterTPCNSigmaPosPr10, int8_t);   //! TPC PID x10 of the positive daughter as Proton
DECLARE_SOA_COLUMN(DaughterTPCNSigmaNegPi10, daughterTPCNSigmaNegPi10, int8_t);   //! TPC PID x10 of the negative daughter as Pion
DECLARE_SOA_COLUMN(DaughterTPCNSigmaNegKa10, daughterTPCNSigmaNegKa10, int8_t);   //! TPC PID x10 of the negative daughter as Kaon
DECLARE_SOA_COLUMN(DaughterTPCNSigmaNegPr10, daughterTPCNSigmaNegPr10, int8_t);   //! TPC PID x10 of the negative daughter as Proton
DECLARE_SOA_COLUMN(DaughterTPCNSigmaBachPi10, daughterTPCNSigmaBachPi10, int8_t); //! TPC PID x10 of the bachelor daughter as Pion
DECLARE_SOA_COLUMN(DaughterTPCNSigmaBachKa10, daughterTPCNSigmaBachKa10, int8_t); //! TPC PID x10 of the bachelor daughter as Kaon
DECLARE_SOA_COLUMN(DaughterTPCNSigmaBachPr10, daughterTPCNSigmaBachPr10, int8_t); //! TPC PID x10 of the bachelor daughter as Proton
DECLARE_SOA_COLUMN(DaughterTOFNSigmaPosPi10, daughterTOFNSigmaPosPi10, int8_t);   //! TOF PID x10 of the positive daughter as Pion
DECLARE_SOA_COLUMN(DaughterTOFNSigmaPosKa10, daughterTOFNSigmaPosKa10, int8_t);   //! TOF PID x10 of the positive daughter as Kaon
DECLARE_SOA_COLUMN(DaughterTOFNSigmaPosPr10, daughterTOFNSigmaPosPr10, int8_t);   //! TOF PID x10 of the positive daughter as Proton
DECLARE_SOA_COLUMN(DaughterTOFNSigmaNegPi10, daughterTOFNSigmaNegPi10, int8_t);   //! TOF PID x10 of the negative daughter as Pion
DECLARE_SOA_COLUMN(DaughterTOFNSigmaNegKa10, daughterTOFNSigmaNegKa10, int8_t);   //! TOF PID x10 of the negative daughter as Kaon
DECLARE_SOA_COLUMN(DaughterTOFNSigmaNegPr10, daughterTOFNSigmaNegPr10, int8_t);   //! TOF PID x10 of the negative daughter as Proton
DECLARE_SOA_COLUMN(DaughterTOFNSigmaBachPi10, daughterTOFNSigmaBachPi10, int8_t); //! TOF PID x10 of the bachelor daughter as Pion
DECLARE_SOA_COLUMN(DaughterTOFNSigmaBachKa10, daughterTOFNSigmaBachKa10, int8_t); //! TOF PID x10 of the bachelor daughter as Kaon
DECLARE_SOA_COLUMN(DaughterTOFNSigmaBachPr10, daughterTOFNSigmaBachPr10, int8_t); //! TOF PID x10 of the bachelor daughter as Proton
// For MC
DECLARE_SOA_INDEX_COLUMN(McParticle, mcParticle); //! Index of the corresponding MC particle
DECLARE_SOA_COLUMN(IsPhysicalPrimary, isPhysicalPrimary, bool);
DECLARE_SOA_COLUMN(ProducedByGenerator, producedByGenerator, bool);
DECLARE_SOA_COLUMN(MotherId, motherId, int);         //! Id of the mother particle
DECLARE_SOA_COLUMN(MotherPDG, motherPDG, int);       //! PDG code of the mother particle
DECLARE_SOA_COLUMN(MotherPt, motherPt, float);       //! Pt of the mother particle
DECLARE_SOA_COLUMN(MotherRap, motherRap, float);     //! Rapidity of the mother particle
DECLARE_SOA_COLUMN(DaughterPDG1, daughterPDG1, int); //! PDG code of the first Daughter particle
DECLARE_SOA_COLUMN(DaughterPDG2, daughterPDG2, int); //! PDG code of the second Daughter particle
DECLARE_SOA_COLUMN(DaughterID1, daughterID1, int);   //! Id of the first Daughter particle
DECLARE_SOA_COLUMN(DaughterID2, daughterID2, int);   //! Id of the second Daughter particle
DECLARE_SOA_COLUMN(SiblingIds, siblingIds, int[2]);  //! Index of the particles with the same mother
DECLARE_SOA_COLUMN(BachTrkID, bachTrkID, int);       //! Id of the bach track from cascade
DECLARE_SOA_COLUMN(V0ID, v0ID, int);                 //! Id of the V0 from cascade
// Dynamic columns
// DCA_xy x10,000
DECLARE_SOA_DYNAMIC_COLUMN(DcaXY, dcaXY,
                           [](int16_t dcaXY10000) { return (float)dcaXY10000 / 10000.f; });
// DCA_z x10,000
DECLARE_SOA_DYNAMIC_COLUMN(DcaZ, dcaZ,
                           [](int16_t dcaZ10000) { return (float)dcaZ10000 / 10000.f; });
// TPC PID return value/10
DECLARE_SOA_DYNAMIC_COLUMN(TpcNSigmaPi, tpcNSigmaPi,
                           [](int8_t tpcNSigmaPi10) { return (float)tpcNSigmaPi10 / 10.f; });
DECLARE_SOA_DYNAMIC_COLUMN(TpcNSigmaKa, tpcNSigmaKa,
                           [](int8_t tpcNSigmaKa10) { return (float)tpcNSigmaKa10 / 10.f; });
DECLARE_SOA_DYNAMIC_COLUMN(TpcNSigmaPr, tpcNSigmaPr,
                           [](int8_t tpcNSigmaPr10) { return (float)tpcNSigmaPr10 / 10.f; });
DECLARE_SOA_DYNAMIC_COLUMN(TofNSigmaPi, tofNSigmaPi,
                           [](int8_t tofNSigmaPi10) { return (float)tofNSigmaPi10 / 10.f; });
DECLARE_SOA_DYNAMIC_COLUMN(TofNSigmaKa, tofNSigmaKa,
                           [](int8_t tofNSigmaKa10) { return (float)tofNSigmaKa10 / 10.f; });
DECLARE_SOA_DYNAMIC_COLUMN(TofNSigmaPr, tofNSigmaPr,
                           [](int8_t tofNSigmaPr10) { return (float)tofNSigmaPr10 / 10.f; });
DECLARE_SOA_DYNAMIC_COLUMN(DaughterTPCNSigmaPosPi, daughterTPCNSigmaPosPi,
                           [](int8_t daughterTPCNSigmaPosPi10) { return (float)daughterTPCNSigmaPosPi10 / 10.f; });
DECLARE_SOA_DYNAMIC_COLUMN(DaughterTPCNSigmaPosKa, daughterTPCNSigmaPosKa,
                           [](int8_t daughterTPCNSigmaPosKa10) { return (float)daughterTPCNSigmaPosKa10 / 10.f; });
DECLARE_SOA_DYNAMIC_COLUMN(DaughterTPCNSigmaPosPr, daughterTPCNSigmaPosPr,
                           [](int8_t daughterTPCNSigmaPosPr10) { return (float)daughterTPCNSigmaPosPr10 / 10.f; });
DECLARE_SOA_DYNAMIC_COLUMN(DaughterTPCNSigmaNegPi, daughterTPCNSigmaNegPi,
                           [](int8_t daughterTPCNSigmaNegPi10) { return (float)daughterTPCNSigmaNegPi10 / 10.f; });
DECLARE_SOA_DYNAMIC_COLUMN(DaughterTPCNSigmaNegKa, daughterTPCNSigmaNegKa,
                           [](int8_t daughterTPCNSigmaNegKa10) { return (float)daughterTPCNSigmaNegKa10 / 10.f; });
DECLARE_SOA_DYNAMIC_COLUMN(DaughterTPCNSigmaNegPr, daughterTPCNSigmaNegPr,
                           [](int8_t daughterTPCNSigmaNegPr10) { return (float)daughterTPCNSigmaNegPr10 / 10.f; });
DECLARE_SOA_DYNAMIC_COLUMN(DaughterTPCNSigmaBachPi, daughterTPCNSigmaBachPi,
                           [](int8_t daughterTPCNSigmaBachPi10) { return (float)daughterTPCNSigmaBachPi10 / 10.f; });
DECLARE_SOA_DYNAMIC_COLUMN(DaughterTPCNSigmaBachKa, daughterTPCNSigmaBachKa,
                           [](int8_t daughterTPCNSigmaBachKa10) { return (float)daughterTPCNSigmaBachKa10 / 10.f; });
DECLARE_SOA_DYNAMIC_COLUMN(DaughterTPCNSigmaBachPr, daughterTPCNSigmaBachPr,
                           [](int8_t daughterTPCNSigmaBachPr10) { return (float)daughterTPCNSigmaBachPr10 / 10.f; });
DECLARE_SOA_DYNAMIC_COLUMN(DaughterTOFNSigmaPosPi, daughterTOFNSigmaPosPi,
                           [](int8_t daughterTOFNSigmaPosPi10) { return (float)daughterTOFNSigmaPosPi10 / 10.f; });
DECLARE_SOA_DYNAMIC_COLUMN(DaughterTOFNSigmaPosKa, daughterTOFNSigmaPosKa,
                           [](int8_t daughterTOFNSigmaPosKa10) { return (float)daughterTOFNSigmaPosKa10 / 10.f; });
DECLARE_SOA_DYNAMIC_COLUMN(DaughterTOFNSigmaPosPr, daughterTOFNSigmaPosPr,
                           [](int8_t daughterTOFNSigmaPosPr10) { return (float)daughterTOFNSigmaPosPr10 / 10.f; });
DECLARE_SOA_DYNAMIC_COLUMN(DaughterTOFNSigmaNegPi, daughterTOFNSigmaNegPi,
                           [](int8_t daughterTOFNSigmaNegPi10) { return (float)daughterTOFNSigmaNegPi10 / 10.f; });
DECLARE_SOA_DYNAMIC_COLUMN(DaughterTOFNSigmaNegKa, daughterTOFNSigmaNegKa,
                           [](int8_t daughterTOFNSigmaNegKa10) { return (float)daughterTOFNSigmaNegKa10 / 10.f; });
DECLARE_SOA_DYNAMIC_COLUMN(DaughterTOFNSigmaNegPr, daughterTOFNSigmaNegPr,
                           [](int8_t daughterTOFNSigmaNegPr10) { return (float)daughterTOFNSigmaNegPr10 / 10.f; });
DECLARE_SOA_DYNAMIC_COLUMN(DaughterTOFNSigmaBachPi, daughterTOFNSigmaBachPi,
                           [](int8_t daughterTOFNSigmaBachPi10) { return (float)daughterTOFNSigmaBachPi10 / 10.f; });
DECLARE_SOA_DYNAMIC_COLUMN(DaughterTOFNSigmaBachKa, daughterTOFNSigmaBachKa,
                           [](int8_t daughterTOFNSigmaBachKa10) { return (float)daughterTOFNSigmaBachKa10 / 10.f; });
DECLARE_SOA_DYNAMIC_COLUMN(DaughterTOFNSigmaBachPr, daughterTOFNSigmaBachPr,
                           [](int8_t daughterTOFNSigmaBachPr10) { return (float)daughterTOFNSigmaBachPr10 / 10.f; });
// TPC signal x10
DECLARE_SOA_DYNAMIC_COLUMN(TpcSignal, tpcSignal,
                           [](int8_t tpcSignal10) { return (float)tpcSignal10 / 10.f; });
// pT, Eta, Phi
// DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt, [](float px, float py) -> float { return RecoDecay::sqrtSumOfSquares(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(Eta, eta, [](float px, float py, float pz) -> float { return RecoDecay::eta(std::array{px, py, pz}); });
DECLARE_SOA_DYNAMIC_COLUMN(Phi, phi, [](float px, float py) -> float { return RecoDecay::phi(px, py); });
// Track flags
DECLARE_SOA_DYNAMIC_COLUMN(PassedITSRefit, passedITSRefit,
                           [](ResoTrackFlags::flagtype trackFlags) -> bool {
                             return ResoTrackFlags::checkFlag(trackFlags, ResoTrackFlags::kPassedITSRefit);
                           });
DECLARE_SOA_DYNAMIC_COLUMN(PassedTPCRefit, passedTPCRefit,
                           [](ResoTrackFlags::flagtype trackFlags) -> bool {
                             return ResoTrackFlags::checkFlag(trackFlags, ResoTrackFlags::kPassedTPCRefit);
                           });
DECLARE_SOA_DYNAMIC_COLUMN(IsGlobalTrackWoDCA, isGlobalTrackWoDCA,
                           [](ResoTrackFlags::flagtype trackFlags) -> bool {
                             return ResoTrackFlags::checkFlag(trackFlags, ResoTrackFlags::kIsGlobalTrackWoDCA);
                           });
DECLARE_SOA_DYNAMIC_COLUMN(IsGlobalTrack, isGlobalTrack,
                           [](ResoTrackFlags::flagtype trackFlags) -> bool {
                             return ResoTrackFlags::checkFlag(trackFlags, ResoTrackFlags::kIsGlobalTrack);
                           });
DECLARE_SOA_DYNAMIC_COLUMN(IsPrimaryTrack, isPrimaryTrack,
                           [](ResoTrackFlags::flagtype trackFlags) -> bool {
                             return ResoTrackFlags::checkFlag(trackFlags, ResoTrackFlags::kIsPrimaryTrack);
                           });
DECLARE_SOA_DYNAMIC_COLUMN(IsPVContributor, isPVContributor,
                           [](ResoTrackFlags::flagtype trackFlags) -> bool {
                             return ResoTrackFlags::checkFlag(trackFlags, ResoTrackFlags::kIsPVContributor);
                           });
DECLARE_SOA_DYNAMIC_COLUMN(HasTOF, hasTOF,
                           [](ResoTrackFlags::flagtype trackFlags) -> bool {
                             return ResoTrackFlags::checkFlag(trackFlags, ResoTrackFlags::kHasTOF);
                           });
DECLARE_SOA_DYNAMIC_COLUMN(Sign, sign,
                           [](ResoTrackFlags::flagtype trackFlags) -> int8_t {
                             return (trackFlags & ResoTrackFlags::kSign) ? 1 : -1;
                           });

} // namespace resodaughter

namespace resomicrodaughter
{
// micro track for primary pion

/// @brief Save TPC & TOF nSigma info with 8-bit variable
struct PidNSigma {
  uint8_t flag;

  /// @brief Constructor: Convert TPC & TOF values and save
  PidNSigma(float TPCnSigma, float TOFnSigma, bool hasTOF)
  {
    uint8_t TPCencoded = encodeNSigma(TPCnSigma);
    uint8_t TOFencoded = hasTOF ? encodeNSigma(TOFnSigma) : 0x0F; // If TOF is not available, set all 4 bits to 1
    flag = (TPCencoded << 4) | TOFencoded;                        // Upper 4 bits = TPC, Lower 4 bits = TOF
  }

  /// @brief Encode 0.2 sigma interval to 0~10 range
  static uint8_t encodeNSigma(float nSigma)
  {
    const float x = std::abs(nSigma);
    if (x <= 1.5)
      return 0; // Return 0 when absolute nSigma is smaller than 1.5
    float t = (x - 1.5) / 0.2;
    int encoded = static_cast<int>(std::ceil(t)); // (1.5,1.7]->1, ..., (3.3,3.5]->10
    if (encoded < 1)
      encoded = 1;
    if (encoded > 10)
      encoded = 10;
    return static_cast<uint8_t>(encoded);
  }

  /// @brief Decode 0~10 value to original 1.5~3.5 sigma range
  static float decodeNSigma(uint8_t encoded)
  {
    if (encoded == 0)
      return 1.5;
    if (encoded > 10)
      encoded = 10;
    return 1.5 + static_cast<float>(encoded) * 0.2;
  }

  /// @brief Check if TOF info is available
  bool hasTOF() const
  {
    return (flag & 0x0F) != 0x0F; // Check if lower 4 bits are not all 1
  }

  /// @brief Restore TPC nSigma value
  static float getTPCnSigma(uint8_t encoded)
  {
    return decodeNSigma((encoded >> 4) & 0x0F); // Extract upper 4 bits
  }

  /// @brief Restore TOF nSigma value (if not available, return NAN)
  static float getTOFnSigma(uint8_t encoded)
  {
    uint8_t TOFencoded = encoded & 0x0F; // Extract lower 4 bits
    return (TOFencoded == 0x0F) ? NAN : decodeNSigma(TOFencoded);
  }

  /// @brief Operator to convert to uint8_t (automatic conversion support)
  operator uint8_t() const
  {
    return flag;
  }
};

DECLARE_SOA_COLUMN(PidNSigmaPiFlag, pidNSigmaPiFlag, uint8_t);        //! Pid flag for the track as Pion
DECLARE_SOA_COLUMN(PidNSigmaKaFlag, pidNSigmaKaFlag, uint8_t);        //! Pid flag for the track as Kaon
DECLARE_SOA_COLUMN(PidNSigmaPrFlag, pidNSigmaPrFlag, uint8_t);        //! Pid flag for the track as Proton
DECLARE_SOA_COLUMN(TrackSelectionFlags, trackSelectionFlags, int8_t); //! Track selection flags
DECLARE_SOA_DYNAMIC_COLUMN(HasTOF, hasTOF,
                           [](uint8_t pidNSigmaFlags) -> bool {
                             return (pidNSigmaFlags & 0x0F) != 0x0F;
                           });

/// @brief DCAxy & DCAz selection flag
struct ResoMicroTrackSelFlag {
  uint8_t flag; // Flag for DCAxy & DCAz selection (8-bit variable)

  /// @brief Default constructor
  ResoMicroTrackSelFlag()
  {
    flag = 0x00;
  }

  /// @brief Constructor: Convert DCAxy/DCAz and save (default 1~15 values)
  ResoMicroTrackSelFlag(float DCAxy, float DCAz)
  {
    uint8_t DCAxyEncoded = encodeDCA(DCAxy);
    uint8_t DCAzEncoded = encodeDCA(DCAz);
    flag = (DCAxyEncoded << 4) | DCAzEncoded; // Upper 4 bits = DCAxy, Lower 4 bits = DCAz
  }

  /// @brief Convert DCA to 1~15 steps (|DCA|<0.1 is saved in 0)
  static uint8_t encodeDCA(float DCA)
  {
    float x = std::fabs(DCA);
    if (x < 0.1)
      return 0;
    int encoded = static_cast<int>(std::ceil((x - 0.1) / 0.1)); // (0.1, 0.2] -> 1, ..., (1.4, 1.5] -> 14
    if (encoded < 1)
      encoded = 1;
    if (encoded > 14)
      encoded = 15;
    return static_cast<uint8_t>(encoded);
  }

  /// @brief Operator to convert to `uint8_t` (for SOA storage)
  operator uint8_t() const
  {
    return flag;
  }

  /// @brief Get DCAxy value
  uint8_t getDCAxyFlag() const
  {
    return (flag >> 4) & 0x0F; // Extract upper 4 bits
  }

  /// @brief Get DCAz value
  uint8_t getDCAzFlag() const
  {
    return flag & 0x0F; // Extract lower 4 bits
  }

  /// @brief Apply DCAxy tight cut (0 value)
  void setDCAxy0()
  {
    flag &= 0x0F; // Set DCAxy to 0 (delete upper 4 bits)
  }

  /// @brief Apply DCAz tight cut (0 value)
  void setDCAz0()
  {
    flag &= 0xF0; // Set DCAz to 0 (delete lower 4 bits)
  }
  /// @brief Decode DCAxy
  static float decodeDCAxy(uint8_t flag_saved)
  {
    uint8_t DCAxyFlag = (flag_saved >> 4) & 0x0F;      // Extract upper 4 bits
    return (DCAxyFlag == 0) ? 0.0f : DCAxyFlag * 0.1f; // Tight cut(0) is 0.0, otherwise flag * 0.1 cm
  }

  /// @brief Decode DCAz
  static float decodeDCAz(uint8_t flag_saved)
  {
    uint8_t DCAzFlag = flag_saved & 0x0F;            // Extract lower 4 bits
    return (DCAzFlag == 0) ? 0.0f : DCAzFlag * 0.1f; // Tight cut(0) is 0.0, otherwise flag * 0.1 cm
  }
};

DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt, [](float px, float py) -> float { return RecoDecay::sqrtSumOfSquares(px, py); });
} // namespace resomicrodaughter

DECLARE_SOA_TABLE(ResoTracks, "AOD", "RESOTRACK",
                  o2::soa::Index<>,
                  resodaughter::ResoCollisionId,
                  resodaughter::Pt,
                  resodaughter::Px,
                  resodaughter::Py,
                  resodaughter::Pz,
                  resodaughter::TpcNClsCrossedRows,
                  resodaughter::TpcNClsFound,
                  resodaughter::DcaXY10000,
                  resodaughter::DcaZ10000,
                  resodaughter::TpcNSigmaPi10,
                  resodaughter::TpcNSigmaKa10,
                  resodaughter::TpcNSigmaPr10,
                  resodaughter::TofNSigmaPi10,
                  resodaughter::TofNSigmaKa10,
                  resodaughter::TofNSigmaPr10,
                  resodaughter::TpcSignal10,
                  resodaughter::TrackFlags,
                  // Dynamic columns
                  resodaughter::TpcNSigmaPi<resodaughter::TpcNSigmaPi10>,
                  resodaughter::TpcNSigmaKa<resodaughter::TpcNSigmaKa10>,
                  resodaughter::TpcNSigmaPr<resodaughter::TpcNSigmaPr10>,
                  resodaughter::TofNSigmaPi<resodaughter::TofNSigmaPi10>,
                  resodaughter::TofNSigmaKa<resodaughter::TofNSigmaKa10>,
                  resodaughter::TofNSigmaPr<resodaughter::TofNSigmaPr10>,
                  resodaughter::TpcSignal<resodaughter::TpcSignal10>,
                  // resodaughter::Pt<resodaughter::Px, resodaughter::Py>,
                  resodaughter::DcaXY<resodaughter::DcaXY10000>,
                  resodaughter::DcaZ<resodaughter::DcaZ10000>,
                  resodaughter::Eta<resodaughter::Px, resodaughter::Py, resodaughter::Pz>,
                  resodaughter::Phi<resodaughter::Px, resodaughter::Py>,
                  resodaughter::PassedITSRefit<resodaughter::TrackFlags>,
                  resodaughter::PassedTPCRefit<resodaughter::TrackFlags>,
                  resodaughter::IsGlobalTrackWoDCA<resodaughter::TrackFlags>,
                  resodaughter::IsGlobalTrack<resodaughter::TrackFlags>,
                  resodaughter::IsPrimaryTrack<resodaughter::TrackFlags>,
                  resodaughter::IsPVContributor<resodaughter::TrackFlags>,
                  resodaughter::HasTOF<resodaughter::TrackFlags>,
                  resodaughter::Sign<resodaughter::TrackFlags>);
using ResoTrack = ResoTracks::iterator;

DECLARE_SOA_TABLE(ResoTrackTracks, "AOD", "RESOTRACKTRACK",
                  resodaughter::TrackId);
using ResoTrackTrack = ResoTrackTracks::iterator;

DECLARE_SOA_TABLE(ResoMicroTracks, "AOD", "RESOMICROTRACK",
                  o2::soa::Index<>,
                  resodaughter::ResoCollisionId,
                  resodaughter::Px,
                  resodaughter::Py,
                  resodaughter::Pz,
                  resomicrodaughter::PidNSigmaPiFlag,
                  resomicrodaughter::PidNSigmaKaFlag,
                  resomicrodaughter::PidNSigmaPrFlag,
                  resomicrodaughter::TrackSelectionFlags,
                  resodaughter::TrackFlags,
                  // Dynamic columns
                  resomicrodaughter::Pt<resodaughter::Px, resodaughter::Py>,
                  resodaughter::Eta<resodaughter::Px, resodaughter::Py, resodaughter::Pz>,
                  resodaughter::Phi<resodaughter::Px, resodaughter::Py>,
                  resodaughter::PassedITSRefit<resodaughter::TrackFlags>,
                  resodaughter::PassedTPCRefit<resodaughter::TrackFlags>,
                  resodaughter::IsGlobalTrackWoDCA<resodaughter::TrackFlags>,
                  resodaughter::IsGlobalTrack<resodaughter::TrackFlags>,
                  resodaughter::IsPrimaryTrack<resodaughter::TrackFlags>,
                  resodaughter::IsPVContributor<resodaughter::TrackFlags>,
                  resomicrodaughter::HasTOF<resomicrodaughter::PidNSigmaPiFlag>,
                  resodaughter::Sign<resodaughter::TrackFlags>);
using ResoMicroTrack = ResoMicroTracks::iterator;

DECLARE_SOA_TABLE(ResoMicroTrackTracks, "AOD", "RESOMICROTRACKTRACK",
                  resodaughter::TrackId);
using ResoMicroTrackTrack = ResoMicroTrackTracks::iterator;

// For DF mixing study
DECLARE_SOA_TABLE(ResoTrackDFs, "AOD", "RESOTRACKDF",
                  o2::soa::Index<>,
                  resodaughter::ResoCollisionDFId,
                  resodaughter::Pt,
                  resodaughter::Px,
                  resodaughter::Py,
                  resodaughter::Pz,
                  resodaughter::TpcNClsCrossedRows,
                  resodaughter::TpcNClsFound,
                  resodaughter::DcaXY10000,
                  resodaughter::DcaZ10000,
                  resodaughter::TpcNSigmaPi10,
                  resodaughter::TpcNSigmaKa10,
                  resodaughter::TpcNSigmaPr10,
                  resodaughter::TofNSigmaPi10,
                  resodaughter::TofNSigmaKa10,
                  resodaughter::TofNSigmaPr10,
                  resodaughter::TpcSignal10,
                  resodaughter::TrackFlags,
                  // Dynamic columns
                  resodaughter::TpcNSigmaPi<resodaughter::TpcNSigmaPi10>,
                  resodaughter::TpcNSigmaKa<resodaughter::TpcNSigmaKa10>,
                  resodaughter::TpcNSigmaPr<resodaughter::TpcNSigmaPr10>,
                  resodaughter::TofNSigmaPi<resodaughter::TofNSigmaPi10>,
                  resodaughter::TofNSigmaKa<resodaughter::TofNSigmaKa10>,
                  resodaughter::TofNSigmaPr<resodaughter::TofNSigmaPr10>,
                  resodaughter::TpcSignal<resodaughter::TpcSignal10>,
                  // resodaughter::Pt<resodaughter::Px, resodaughter::Py>,
                  resodaughter::DcaXY<resodaughter::DcaXY10000>,
                  resodaughter::DcaZ<resodaughter::DcaZ10000>,
                  resodaughter::Eta<resodaughter::Px, resodaughter::Py, resodaughter::Pz>,
                  resodaughter::Phi<resodaughter::Px, resodaughter::Py>,
                  resodaughter::PassedITSRefit<resodaughter::TrackFlags>,
                  resodaughter::PassedTPCRefit<resodaughter::TrackFlags>,
                  resodaughter::IsGlobalTrackWoDCA<resodaughter::TrackFlags>,
                  resodaughter::IsGlobalTrack<resodaughter::TrackFlags>,
                  resodaughter::IsPrimaryTrack<resodaughter::TrackFlags>,
                  resodaughter::IsPVContributor<resodaughter::TrackFlags>,
                  resodaughter::HasTOF<resodaughter::TrackFlags>,
                  resodaughter::Sign<resodaughter::TrackFlags>);
using ResoTrackDF = ResoTrackDFs::iterator;

DECLARE_SOA_TABLE(ResoV0s, "AOD", "RESOV0",
                  o2::soa::Index<>,
                  resodaughter::ResoCollisionId,
                  resodaughter::Pt,
                  resodaughter::Px,
                  resodaughter::Py,
                  resodaughter::Pz,
                  resodaughter::Indices,
                  resodaughter::DaughterTPCNSigmaPosPi10,
                  resodaughter::DaughterTPCNSigmaPosKa10,
                  resodaughter::DaughterTPCNSigmaPosPr10,
                  resodaughter::DaughterTPCNSigmaNegPi10,
                  resodaughter::DaughterTPCNSigmaNegKa10,
                  resodaughter::DaughterTPCNSigmaNegPr10,
                  resodaughter::DaughterTOFNSigmaPosPi10,
                  resodaughter::DaughterTOFNSigmaPosKa10,
                  resodaughter::DaughterTOFNSigmaPosPr10,
                  resodaughter::DaughterTOFNSigmaNegPi10,
                  resodaughter::DaughterTOFNSigmaNegKa10,
                  resodaughter::DaughterTOFNSigmaNegPr10,
                  resodaughter::V0CosPA,
                  resodaughter::DaughDCA,
                  v0data::DCAPosToPV,
                  v0data::DCANegToPV,
                  v0data::DCAV0ToPV,
                  resodaughter::MLambda,
                  resodaughter::MAntiLambda,
                  resodaughter::MK0Short,
                  resodaughter::TransRadius,
                  resodaughter::DecayVtxX,
                  resodaughter::DecayVtxY,
                  resodaughter::DecayVtxZ,
                  resodaughter::Alpha,
                  resodaughter::QtArm,
                  // resodaughter::Pt<resodaughter::Px, resodaughter::Py>,
                  resodaughter::Eta<resodaughter::Px, resodaughter::Py, resodaughter::Pz>,
                  resodaughter::Phi<resodaughter::Px, resodaughter::Py>,
                  resodaughter::DaughterTPCNSigmaPosPi<resodaughter::DaughterTPCNSigmaPosPi10>,
                  resodaughter::DaughterTPCNSigmaPosKa<resodaughter::DaughterTPCNSigmaPosKa10>,
                  resodaughter::DaughterTPCNSigmaPosPr<resodaughter::DaughterTPCNSigmaPosPr10>,
                  resodaughter::DaughterTPCNSigmaNegPi<resodaughter::DaughterTPCNSigmaNegPi10>,
                  resodaughter::DaughterTPCNSigmaNegKa<resodaughter::DaughterTPCNSigmaNegKa10>,
                  resodaughter::DaughterTPCNSigmaNegPr<resodaughter::DaughterTPCNSigmaNegPr10>,
                  resodaughter::DaughterTOFNSigmaPosPi<resodaughter::DaughterTOFNSigmaPosPi10>,
                  resodaughter::DaughterTOFNSigmaPosKa<resodaughter::DaughterTOFNSigmaPosKa10>,
                  resodaughter::DaughterTOFNSigmaPosPr<resodaughter::DaughterTOFNSigmaPosPr10>,
                  resodaughter::DaughterTOFNSigmaNegPi<resodaughter::DaughterTOFNSigmaNegPi10>,
                  resodaughter::DaughterTOFNSigmaNegKa<resodaughter::DaughterTOFNSigmaNegKa10>,
                  resodaughter::DaughterTOFNSigmaNegPr<resodaughter::DaughterTOFNSigmaNegPr10>);
using ResoV0 = ResoV0s::iterator;

DECLARE_SOA_TABLE(ResoV0V0s, "AOD", "RESOV0V0",
                  resodaughter::V0Id);
using ResoV0V0 = ResoV0V0s::iterator;

DECLARE_SOA_TABLE(ResoCascades, "AOD", "RESOCASCADE",
                  o2::soa::Index<>,
                  resodaughter::ResoCollisionId,
                  resodaughter::Pt,
                  resodaughter::Px,
                  resodaughter::Py,
                  resodaughter::Pz,
                  resodaughter::CascadeIndices,
                  resodaughter::DaughterTPCNSigmaPosPi10,
                  resodaughter::DaughterTPCNSigmaPosKa10,
                  resodaughter::DaughterTPCNSigmaPosPr10,
                  resodaughter::DaughterTPCNSigmaNegPi10,
                  resodaughter::DaughterTPCNSigmaNegKa10,
                  resodaughter::DaughterTPCNSigmaNegPr10,
                  resodaughter::DaughterTPCNSigmaBachPi10,
                  resodaughter::DaughterTPCNSigmaBachKa10,
                  resodaughter::DaughterTPCNSigmaBachPr10,
                  resodaughter::DaughterTOFNSigmaPosPi10,
                  resodaughter::DaughterTOFNSigmaPosKa10,
                  resodaughter::DaughterTOFNSigmaPosPr10,
                  resodaughter::DaughterTOFNSigmaNegPi10,
                  resodaughter::DaughterTOFNSigmaNegKa10,
                  resodaughter::DaughterTOFNSigmaNegPr10,
                  resodaughter::DaughterTOFNSigmaBachPi10,
                  resodaughter::DaughterTOFNSigmaBachKa10,
                  resodaughter::DaughterTOFNSigmaBachPr10,
                  resodaughter::V0CosPA,
                  resodaughter::CascCosPA,
                  resodaughter::DaughDCA,
                  resodaughter::CascDaughDCA,
                  cascdata::DCAPosToPV,
                  cascdata::DCANegToPV,
                  cascdata::DCABachToPV,
                  v0data::DCAV0ToPV,
                  cascdata::DCAXYCascToPV,
                  cascdata::DCAZCascToPV,
                  cascdata::Sign,
                  resodaughter::MLambda,
                  resodaughter::MXi,
                  resodaughter::TransRadius,
                  resodaughter::CascTransRadius,
                  resodaughter::DecayVtxX,
                  resodaughter::DecayVtxY,
                  resodaughter::DecayVtxZ,
                  // resodaughter::Pt<resodaughter::Px, resodaughter::Py>,
                  resodaughter::Eta<resodaughter::Px, resodaughter::Py, resodaughter::Pz>,
                  resodaughter::Phi<resodaughter::Px, resodaughter::Py>,
                  resodaughter::DaughterTPCNSigmaPosPi<resodaughter::DaughterTPCNSigmaPosPi10>,
                  resodaughter::DaughterTPCNSigmaPosKa<resodaughter::DaughterTPCNSigmaPosKa10>,
                  resodaughter::DaughterTPCNSigmaPosPr<resodaughter::DaughterTPCNSigmaPosPr10>,
                  resodaughter::DaughterTPCNSigmaNegPi<resodaughter::DaughterTPCNSigmaNegPi10>,
                  resodaughter::DaughterTPCNSigmaNegKa<resodaughter::DaughterTPCNSigmaNegKa10>,
                  resodaughter::DaughterTPCNSigmaNegPr<resodaughter::DaughterTPCNSigmaNegPr10>,
                  resodaughter::DaughterTPCNSigmaBachPi<resodaughter::DaughterTPCNSigmaBachPi10>,
                  resodaughter::DaughterTPCNSigmaBachKa<resodaughter::DaughterTPCNSigmaBachKa10>,
                  resodaughter::DaughterTPCNSigmaBachPr<resodaughter::DaughterTPCNSigmaBachPr10>,
                  resodaughter::DaughterTOFNSigmaPosPi<resodaughter::DaughterTOFNSigmaPosPi10>,
                  resodaughter::DaughterTOFNSigmaPosKa<resodaughter::DaughterTOFNSigmaPosKa10>,
                  resodaughter::DaughterTOFNSigmaPosPr<resodaughter::DaughterTOFNSigmaPosPr10>,
                  resodaughter::DaughterTOFNSigmaNegPi<resodaughter::DaughterTOFNSigmaNegPi10>,
                  resodaughter::DaughterTOFNSigmaNegKa<resodaughter::DaughterTOFNSigmaNegKa10>,
                  resodaughter::DaughterTOFNSigmaNegPr<resodaughter::DaughterTOFNSigmaNegPr10>,
                  resodaughter::DaughterTOFNSigmaBachPi<resodaughter::DaughterTOFNSigmaBachPi10>,
                  resodaughter::DaughterTOFNSigmaBachKa<resodaughter::DaughterTOFNSigmaBachKa10>,
                  resodaughter::DaughterTOFNSigmaBachPr<resodaughter::DaughterTOFNSigmaBachPr10>);
using ResoCascade = ResoCascades::iterator;

DECLARE_SOA_TABLE(ResoCascadeCascades, "AOD", "RESOCASCADECASCADE",
                  resodaughter::CascadeId);
using ResoCascadeCascade = ResoCascadeCascades::iterator;

DECLARE_SOA_TABLE(ResoCascadeDFs, "AOD", "RESOCASCADEDF",
                  o2::soa::Index<>,
                  resodaughter::ResoCollisionDFId,
                  resodaughter::Pt,
                  resodaughter::Px,
                  resodaughter::Py,
                  resodaughter::Pz,
                  resodaughter::CascadeIndices,
                  resodaughter::DaughterTPCNSigmaPosPi10,
                  resodaughter::DaughterTPCNSigmaPosKa10,
                  resodaughter::DaughterTPCNSigmaPosPr10,
                  resodaughter::DaughterTPCNSigmaNegPi10,
                  resodaughter::DaughterTPCNSigmaNegKa10,
                  resodaughter::DaughterTPCNSigmaNegPr10,
                  resodaughter::DaughterTPCNSigmaBachPi10,
                  resodaughter::DaughterTPCNSigmaBachKa10,
                  resodaughter::DaughterTPCNSigmaBachPr10,
                  resodaughter::DaughterTOFNSigmaPosPi10,
                  resodaughter::DaughterTOFNSigmaPosKa10,
                  resodaughter::DaughterTOFNSigmaPosPr10,
                  resodaughter::DaughterTOFNSigmaNegPi10,
                  resodaughter::DaughterTOFNSigmaNegKa10,
                  resodaughter::DaughterTOFNSigmaNegPr10,
                  resodaughter::DaughterTOFNSigmaBachPi10,
                  resodaughter::DaughterTOFNSigmaBachKa10,
                  resodaughter::DaughterTOFNSigmaBachPr10,
                  resodaughter::V0CosPA,
                  resodaughter::CascCosPA,
                  resodaughter::DaughDCA,
                  resodaughter::CascDaughDCA,
                  cascdata::DCAPosToPV,
                  cascdata::DCANegToPV,
                  cascdata::DCABachToPV,
                  v0data::DCAV0ToPV,
                  cascdata::DCAXYCascToPV,
                  cascdata::DCAZCascToPV,
                  cascdata::Sign,
                  resodaughter::MLambda,
                  resodaughter::MXi,
                  resodaughter::TransRadius,
                  resodaughter::CascTransRadius,
                  resodaughter::DecayVtxX,
                  resodaughter::DecayVtxY,
                  resodaughter::DecayVtxZ,
                  // resodaughter::Pt<resodaughter::Px, resodaughter::Py>,
                  resodaughter::Eta<resodaughter::Px, resodaughter::Py, resodaughter::Pz>,
                  resodaughter::Phi<resodaughter::Px, resodaughter::Py>,
                  resodaughter::DaughterTPCNSigmaPosPi<resodaughter::DaughterTPCNSigmaPosPi10>,
                  resodaughter::DaughterTPCNSigmaPosKa<resodaughter::DaughterTPCNSigmaPosKa10>,
                  resodaughter::DaughterTPCNSigmaPosPr<resodaughter::DaughterTPCNSigmaPosPr10>,
                  resodaughter::DaughterTPCNSigmaNegPi<resodaughter::DaughterTPCNSigmaNegPi10>,
                  resodaughter::DaughterTPCNSigmaNegKa<resodaughter::DaughterTPCNSigmaNegKa10>,
                  resodaughter::DaughterTPCNSigmaNegPr<resodaughter::DaughterTPCNSigmaNegPr10>,
                  resodaughter::DaughterTPCNSigmaBachPi<resodaughter::DaughterTPCNSigmaBachPi10>,
                  resodaughter::DaughterTPCNSigmaBachKa<resodaughter::DaughterTPCNSigmaBachKa10>,
                  resodaughter::DaughterTPCNSigmaBachPr<resodaughter::DaughterTPCNSigmaBachPr10>,
                  resodaughter::DaughterTOFNSigmaPosPi<resodaughter::DaughterTOFNSigmaPosPi10>,
                  resodaughter::DaughterTOFNSigmaPosKa<resodaughter::DaughterTOFNSigmaPosKa10>,
                  resodaughter::DaughterTOFNSigmaPosPr<resodaughter::DaughterTOFNSigmaPosPr10>,
                  resodaughter::DaughterTOFNSigmaNegPi<resodaughter::DaughterTOFNSigmaNegPi10>,
                  resodaughter::DaughterTOFNSigmaNegKa<resodaughter::DaughterTOFNSigmaNegKa10>,
                  resodaughter::DaughterTOFNSigmaNegPr<resodaughter::DaughterTOFNSigmaNegPr10>,
                  resodaughter::DaughterTOFNSigmaBachPi<resodaughter::DaughterTOFNSigmaBachPi10>,
                  resodaughter::DaughterTOFNSigmaBachKa<resodaughter::DaughterTOFNSigmaBachKa10>,
                  resodaughter::DaughterTOFNSigmaBachPr<resodaughter::DaughterTOFNSigmaBachPr10>);
using ResoCascadeDF = ResoCascadeDFs::iterator;

DECLARE_SOA_TABLE(ResoMCTracks, "AOD", "RESOMCTRACK",
                  mcparticle::PdgCode,
                  resodaughter::MotherId,
                  resodaughter::MotherPDG,
                  resodaughter::SiblingIds,
                  resodaughter::IsPhysicalPrimary,
                  resodaughter::ProducedByGenerator);
using ResoMCTrack = ResoMCTracks::iterator;

DECLARE_SOA_TABLE(ResoMCV0s, "AOD", "RESOMCV0",
                  mcparticle::PdgCode,
                  resodaughter::MotherId,
                  resodaughter::MotherPDG,
                  resodaughter::MotherPt,
                  resodaughter::MotherRap,
                  resodaughter::DaughterID1,
                  resodaughter::DaughterID2,
                  resodaughter::DaughterPDG1,
                  resodaughter::DaughterPDG2,
                  resodaughter::IsPhysicalPrimary,
                  resodaughter::ProducedByGenerator);
using ResoMCV0 = ResoMCV0s::iterator;

DECLARE_SOA_TABLE(ResoMCCascades, "AOD", "RESOMCCASCADE",
                  mcparticle::PdgCode,
                  resodaughter::MotherId,
                  resodaughter::MotherPDG,
                  resodaughter::MotherPt,
                  resodaughter::MotherRap,
                  resodaughter::BachTrkID,
                  resodaughter::V0ID,
                  resodaughter::DaughterPDG1,
                  resodaughter::DaughterPDG2,
                  resodaughter::IsPhysicalPrimary,
                  resodaughter::ProducedByGenerator);
using ResoMCCascade = ResoMCCascades::iterator;

DECLARE_SOA_TABLE(ResoMCParents, "AOD", "RESOMCPARENT",
                  o2::soa::Index<>,
                  resodaughter::ResoCollisionId,
                  resodaughter::McParticleId,
                  mcparticle::PdgCode,
                  resodaughter::DaughterPDG1,
                  resodaughter::DaughterPDG2,
                  resodaughter::IsPhysicalPrimary,
                  resodaughter::ProducedByGenerator,
                  resodaughter::Pt,
                  resodaughter::Px,
                  resodaughter::Py,
                  resodaughter::Pz,
                  mcparticle::Y,
                  mcparticle::E,
                  mcparticle::StatusCode,
                  // resodaughter::Pt<resodaughter::Px, resodaughter::Py>,
                  resodaughter::Eta<resodaughter::Px, resodaughter::Py, resodaughter::Pz>,
                  resodaughter::Phi<resodaughter::Px, resodaughter::Py>);
using ResoMCParent = ResoMCParents::iterator;

using Reso2TracksExt = soa::Join<aod::FullTracks, aod::TracksDCA>; // without Extra
using Reso2TracksMC = soa::Join<aod::FullTracks, McTrackLabels>;
using Reso2TracksPID = soa::Join<aod::FullTracks, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr, aod::pidTOFPi, aod::pidTOFKa, aod::pidTOFPr>;
using Reso2TracksPIDExt = soa::Join<Reso2TracksPID, aod::TracksDCA, aod::TrackSelection, aod::TrackSelectionExtension>; // Without Extra

using ResoCollisionCandidates = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0As, aod::Mults>;
using ResoRun2CollisionCandidates = soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0Ms>;
using ResoCollisionCandidatesMC = soa::Join<ResoCollisionCandidates, aod::McCollisionLabels>;
using ResoRun2CollisionCandidatesMC = soa::Join<ResoRun2CollisionCandidates, aod::McCollisionLabels>;
using ResoTrackCandidates = aod::Reso2TracksPIDExt;
using ResoTrackCandidatesMC = soa::Join<ResoTrackCandidates, aod::McTrackLabels>;
using ResoV0Candidates = aod::V0Datas;
using ResoV0CandidatesMC = soa::Join<ResoV0Candidates, aod::McV0Labels>;
using ResoCascadesCandidates = aod::CascDatas;
using ResoCascadesCandidatesMC = soa::Join<ResoCascadesCandidates, aod::McCascLabels>;
using BCsWithRun2Info = soa::Join<aod::BCs, aod::Run2BCInfos, aod::Timestamps>;

} // namespace o2::aod
#endif // PWGLF_DATAMODEL_LFRESONANCETABLES_H_
