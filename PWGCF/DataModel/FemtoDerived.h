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
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
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
// Define different methods for the event mixing
enum CollisionBinning {
  kMult,               //! Bin collision in number of charged tracks for mixing
  kMultPercentile,     //! Bin collision in multiplicity percentile for mixing
  kMultMultPercentile, //! Bin collision in number of charged tracks and multiplicity percentile for mixing
  kNCollisionBinning
};

DECLARE_SOA_COLUMN(MultV0M, multV0M, float);       //! V0M multiplicity
DECLARE_SOA_COLUMN(MultNtr, multNtr, int);         //! multiplicity of charged tracks as defined in the producer
DECLARE_SOA_COLUMN(Sphericity, sphericity, float); //! Sphericity of the event
DECLARE_SOA_COLUMN(MagField, magField, float);     //! Magnetic field of the event

using BitMaskType = uint32_t; //! Definition of the data type for the collision masks

DECLARE_SOA_COLUMN(BitMaskTrackOne, bitmaskTrackOne, BitMaskType);     //! Bit for track one
DECLARE_SOA_COLUMN(BitMaskTrackTwo, bitmaskTrackTwo, BitMaskType);     //! Bit for track two
DECLARE_SOA_COLUMN(BitMaskTrackThree, bitmaskTrackThree, BitMaskType); //! Bit for track three

DECLARE_SOA_COLUMN(Downsample, downsample, bool); //! Flag for downsampling
} // namespace femtodreamcollision

DECLARE_SOA_TABLE_STAGED(FDCollisions, "FDCOLLISION",
                         o2::soa::Index<>,
                         o2::aod::collision::PosZ,
                         femtodreamcollision::MultV0M,
                         femtodreamcollision::MultNtr,
                         femtodreamcollision::Sphericity,
                         femtodreamcollision::MagField);
using FDCollision = FDCollisions::iterator;

DECLARE_SOA_TABLE(FDColMasks, "AOD", "FDCOLMASK",
                  femtodreamcollision::BitMaskTrackOne,
                  femtodreamcollision::BitMaskTrackTwo,
                  femtodreamcollision::BitMaskTrackThree);

DECLARE_SOA_TABLE(FDDownSample, "AOD", "FDDOWNSAMPLE",
                  femtodreamcollision::Downsample);

namespace femtodreamMCcollision
{
DECLARE_SOA_COLUMN(MultMCgenPartEta08, multMCgenPartEta08, int); //! Multiplicity of the event as given by the generator in |eta|<0.8
}

DECLARE_SOA_TABLE_STAGED(FDMCCollisions, "FDMCCOLLISION",
                         o2::soa::Index<>,
                         femtodreamMCcollision::MultMCgenPartEta08);
using FDMCCollision = FDMCCollisions::iterator;

namespace mcfdcolllabel
{
DECLARE_SOA_INDEX_COLUMN(FDMCCollision, fdMCCollision); //! MC collision for femtodreamcollision
}
DECLARE_SOA_TABLE_STAGED(FDMCCollLabels, "FDMCCollLabel", mcfdcolllabel::FDMCCollisionId);

/// FemtoDreamTrack
namespace femtodreamparticle
{
/// Distinuishes the different particle types
enum ParticleType {
  kTrack,   //! Track
  kV0,      //! V0
  kV0Child, //! Child track of a V0
  kCascade, //! Cascade
  kCascadeV0,
  kCascadeV0Child,
  kCascadeBachelor, //! Bachelor track of a cascade
  kCharmHadron,     //! Bachelor track of a cascade
  kNParticleTypes   //! Number of particle types
};

enum MomentumType {
  kPt,    //! transverse momentum
  kPreco, //! reconstructed/propagated momentum at the vertex
  kPtpc   //! momentum at the inner wall of the TPC (useful for PID plots)
};

static constexpr std::string_view ParticleTypeName[kNParticleTypes] = {"Tracks", "V0", "V0Child", "Cascade", "CascadeV0", "CascadeV0Child", "CascadeBachelor", "CharmHadron"}; //! Naming of the different particle types
static constexpr std::string_view TempFitVarName[kNParticleTypes] = {"/hDCAxy", "/hCPA", "/hDCAxy", "/hCPA", "/hCPA", "/hDCAxy", "/hDCAxy", "/hCPA"};

using cutContainerType = uint32_t; //! Definition of the data type for the bit-wise container for the different selection criteria

enum TrackType {
  kNoChild,    //! Not any child
  kPosChild,   //! Positive V0 child
  kNegChild,   //! Negative V0 child
  kBachelor,   //! Bachelor Cascade child
  kNTrackTypes //! Number of child types
};

static constexpr std::string_view TrackTypeName[kNTrackTypes] = {"Trk", "Pos", "Neg", "Bach"}; //! Naming of the different particle types

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
DECLARE_SOA_COLUMN(TPCNSigmaTr, tpcNSigmaTr, float); //! Nsigma separation with the TPC detector for triton
DECLARE_SOA_COLUMN(TPCNSigmaHe, tpcNSigmaHe, float); //! Nsigma separation with the TPC detector for helium3
DECLARE_SOA_COLUMN(TOFNSigmaEl, tofNSigmaEl, float); //! Nsigma separation with the TOF detector for electron
DECLARE_SOA_COLUMN(TOFNSigmaPi, tofNSigmaPi, float); //! Nsigma separation with the TOF detector for pion
DECLARE_SOA_COLUMN(TOFNSigmaKa, tofNSigmaKa, float); //! Nsigma separation with the TOF detector for kaon
DECLARE_SOA_COLUMN(TOFNSigmaPr, tofNSigmaPr, float); //! Nsigma separation with the TOF detector for proton
DECLARE_SOA_COLUMN(TOFNSigmaDe, tofNSigmaDe, float); //! Nsigma separation with the TOF detector for deuteron
DECLARE_SOA_COLUMN(TOFNSigmaTr, tofNSigmaTr, float); //! Nsigma separation with the TOF detector for triton
DECLARE_SOA_COLUMN(TOFNSigmaHe, tofNSigmaHe, float); //! Nsigma separation with the TOF detector for helium3
DECLARE_SOA_COLUMN(ITSSignal, itsSignal, float);
DECLARE_SOA_COLUMN(ITSNSigmaEl, itsNSigmaEl, float); //! Nsigma separation with the Its detector for electron
DECLARE_SOA_COLUMN(ITSNSigmaPi, itsNSigmaPi, float); //! Nsigma separation with the Its detector for pion
DECLARE_SOA_COLUMN(ITSNSigmaKa, itsNSigmaKa, float); //! Nsigma separation with the Its detector for kaon
DECLARE_SOA_COLUMN(ITSNSigmaPr, itsNSigmaPr, float); //! Nsigma separation with the Its detector for proton
DECLARE_SOA_COLUMN(ITSNSigmaDe, itsNSigmaDe, float); //! Nsigma separation with the Its detector for deuteron
DECLARE_SOA_COLUMN(ITSNSigmaTr, itsNSigmaTr, float); //! Nsigma separation with the Its detector for triton
DECLARE_SOA_COLUMN(ITSNSigmaHe, itsNSigmaHe, float); //! Nsigma separation with the Its detector for helium3
DECLARE_SOA_COLUMN(DaughDCA, daughDCA, float);       //! DCA between daughters
DECLARE_SOA_COLUMN(TransRadius, transRadius, float); //! Transverse radius of the decay vertex
DECLARE_SOA_COLUMN(DecayVtxX, decayVtxX, float);     //! X position of the decay vertex
DECLARE_SOA_COLUMN(DecayVtxY, decayVtxY, float);     //! Y position of the decay vertex
DECLARE_SOA_COLUMN(DecayVtxZ, decayVtxZ, float);     //! Z position of the decay vertex
DECLARE_SOA_COLUMN(MKaon, mKaon, float);             //! The invariant mass of V0 candidate, assuming kaon
} // namespace femtodreamparticle

namespace fdhf
{

enum CharmHadronMassHypo {
  wrongParticle = 0,
  lcToPKPi = 1,
  lcToPiKP = 2
};

DECLARE_SOA_COLUMN(TrackId, trackId, int);                          //! track id to match associate particle with charm hadron prongs
DECLARE_SOA_COLUMN(Charge, charge, int8_t);                         //! Charge of charm hadron
DECLARE_SOA_COLUMN(Prong0Id, prong0Id, int);                        //! Track id of charm hadron prong0
DECLARE_SOA_COLUMN(Prong1Id, prong1Id, int);                        //! Track id of charm hadron prong1
DECLARE_SOA_COLUMN(Prong2Id, prong2Id, int);                        //! Track id of charm hadron prong2
DECLARE_SOA_COLUMN(Prong0Pt, prong0Pt, float);                      //! Track pT of charm hadron prong0
DECLARE_SOA_COLUMN(Prong1Pt, prong1Pt, float);                      //! Track pT of charm hadron prong1
DECLARE_SOA_COLUMN(Prong2Pt, prong2Pt, float);                      //! Track pT of charm hadron prong2
DECLARE_SOA_COLUMN(Prong0Eta, prong0Eta, float);                    //! Track eta of charm hadron prong0
DECLARE_SOA_COLUMN(Prong1Eta, prong1Eta, float);                    //! Track eta of charm hadron prong1
DECLARE_SOA_COLUMN(Prong2Eta, prong2Eta, float);                    //! Track eta of charm hadron prong2
DECLARE_SOA_COLUMN(Prong0Phi, prong0Phi, float);                    //! Track phi of charm hadron prong0
DECLARE_SOA_COLUMN(Prong1Phi, prong1Phi, float);                    //! Track phi of charm hadron prong1
DECLARE_SOA_COLUMN(Prong2Phi, prong2Phi, float);                    //! Track phi of charm hadron prong2
DECLARE_SOA_COLUMN(CandidateSelFlag, candidateSelFlag, int8_t);     //! Selection of mass hypothesis for charm hadron (1 for Lc -> pkpi, 2 for Lc -> pikp)
DECLARE_SOA_COLUMN(BDTBkg, bdtBkg, float);                          //! Background score using Boosted Decision Tree for charm hadron
DECLARE_SOA_COLUMN(BDTPrompt, bdtPrompt, float);                    //! Prompt signal score using Boosted Decision Tree for charm hadron
DECLARE_SOA_COLUMN(BDTFD, bdtFD, float);                            //! Feed-down score using Boosted Decision Tree for charm hadron
DECLARE_SOA_COLUMN(FlagMc, flagMc, int8_t);                         //! To select MC particle among charm hadrons, { DplusToPiKPi = 1, LcToPKPi = 2, DsToKKPi = 4, XicToPKP = 8, N3ProngD = 2ecays };
DECLARE_SOA_COLUMN(OriginMcRec, originMcRec, int8_t);               //! flag for reconstruction level matching (1 for prompt, 2 for non-prompt)
DECLARE_SOA_COLUMN(OriginMcGen, originMcGen, int8_t);               //! flag for generator level matching (1 for prompt, 2 for non-prompt)
DECLARE_SOA_COLUMN(IsCandidateSwapped, isCandidateSwapped, int8_t); //! swapping of the prongs order (0 for Lc -> pkpi, 1 for Lc -> pikp)
DECLARE_SOA_COLUMN(PtAssoc, ptAssoc, float);                        //! Transverse momentum of associate femto particle
DECLARE_SOA_COLUMN(Kstar, kstar, float);                            //! Relative momentum in particles pair frame
DECLARE_SOA_COLUMN(KT, kT, float);                                  //! kT distribution of particle pairs
DECLARE_SOA_COLUMN(MT, mT, float);                                  //! Transverse mass distribution
DECLARE_SOA_COLUMN(CharmM, charmM, float);                          //! Charm hadron mass
DECLARE_SOA_COLUMN(CharmPt, charmPt, float);                        //! Transverse momentum of charm hadron for result task
DECLARE_SOA_COLUMN(Mult, mult, int);                                //! Charge particle multiplicity
DECLARE_SOA_COLUMN(MultPercentile, multPercentile, float);          //! Multiplicity precentile
DECLARE_SOA_COLUMN(PairSign, pairSign, int8_t);                     //! Selection between like sign (1) and unlike sign pair (2)
DECLARE_SOA_COLUMN(ProcessType, processType, int64_t);              //! Selection between same-event (1), and mixed-event (2)
DECLARE_SOA_DYNAMIC_COLUMN(M, m,                                    //!
                           [](float pt0, float phi0, float eta0, float pt1, float phi1, float eta1, float pt2, float phi2, float eta2, const std::array<double, 3>& m) -> float { return RecoDecay::m(std::array{
                                                                                                                                                                                                        RecoDecayPtEtaPhi::pVector(pt0, eta0, phi0),
                                                                                                                                                                                                        RecoDecayPtEtaPhi::pVector(pt1, eta1, phi1),
                                                                                                                                                                                                        RecoDecayPtEtaPhi::pVector(pt2, eta2, phi2)},
                                                                                                                                                                                                      m); }); //! Charm hadron mass
DECLARE_SOA_DYNAMIC_COLUMN(P, p,                                                                                                                                                                                                                                                                                                                                 //!
                           [](float pt0, float phi0, float eta0, float pt1, float phi1, float eta1, float pt2, float phi2, float eta2) -> float { return RecoDecay::p(RecoDecay::pVec(
                                                                                                                                                    RecoDecayPtEtaPhi::pVector(pt0, eta0, phi0),
                                                                                                                                                    RecoDecayPtEtaPhi::pVector(pt1, eta1, phi1),
                                                                                                                                                    RecoDecayPtEtaPhi::pVector(pt2, eta2, phi2))); }); //! Momentum of charm hadron
DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt,                                                                                                                                                                                                                                                                                                 //!
                           [](float pt0, float phi0, float eta0, float pt1, float phi1, float eta1, float pt2, float phi2, float eta2) -> float { return RecoDecay::pt(RecoDecay::pVec(
                                                                                                                                                    RecoDecayPtEtaPhi::pVector(pt0, eta0, phi0),
                                                                                                                                                    RecoDecayPtEtaPhi::pVector(pt1, eta1, phi1),
                                                                                                                                                    RecoDecayPtEtaPhi::pVector(pt2, eta2, phi2))); }); //! Transverse momentum of charm hadron
DECLARE_SOA_DYNAMIC_COLUMN(Phi, phi,                                                                                                                                                                                                                                                                                                //!
                           [](float pt0, float phi0, float eta0, float pt1, float phi1, float eta1, float pt2, float phi2, float eta2) -> float { return RecoDecay::phi(RecoDecay::pVec(
                                                                                                                                                    RecoDecayPtEtaPhi::pVector(pt0, eta0, phi0),
                                                                                                                                                    RecoDecayPtEtaPhi::pVector(pt1, eta1, phi1),
                                                                                                                                                    RecoDecayPtEtaPhi::pVector(pt2, eta2, phi2))); }); //! Phi distribution of charm hadron
DECLARE_SOA_DYNAMIC_COLUMN(Y, y,                                                                                                                                                                                                                                                                                                     //!
                           [](float pt0, float phi0, float eta0, float pt1, float phi1, float eta1, float pt2, float phi2, float eta2) -> float { return RecoDecay::y(RecoDecay::pVec(
                                                                                                                                                                        RecoDecayPtEtaPhi::pVector(pt0, eta0, phi0),
                                                                                                                                                                        RecoDecayPtEtaPhi::pVector(pt1, eta1, phi1),
                                                                                                                                                                        RecoDecayPtEtaPhi::pVector(pt2, eta2, phi2)),
                                                                                                                                                                      o2::constants::physics::MassLambdaCPlus); }); //! Rapidity distribution of charm hadron
DECLARE_SOA_DYNAMIC_COLUMN(Eta, eta,                                                                                                                                                                                                                                                                                                                                        //!
                           [](float pt0, float phi0, float eta0, float pt1, float phi1, float eta1, float pt2, float phi2, float eta2) -> float { return RecoDecay::eta(RecoDecay::pVec(
                                                                                                                                                    RecoDecayPtEtaPhi::pVector(pt0, eta0, phi0),
                                                                                                                                                    RecoDecayPtEtaPhi::pVector(pt1, eta1, phi1),
                                                                                                                                                    RecoDecayPtEtaPhi::pVector(pt2, eta2, phi2))); }); //! Eta distribution of charm hadron

} // namespace fdhf

DECLARE_SOA_TABLE(FDHfCand, "AOD", "FDHFCAND", //! Table to store the derived data for charm hadron candidates
                  o2::soa::Index<>,
                  femtodreamparticle::FDCollisionId,
                  fdhf::Charge,
                  fdhf::Prong0Id,
                  fdhf::Prong1Id,
                  fdhf::Prong2Id,
                  fdhf::Prong0Pt,
                  fdhf::Prong1Pt,
                  fdhf::Prong2Pt,
                  fdhf::Prong0Eta,
                  fdhf::Prong1Eta,
                  fdhf::Prong2Eta,
                  fdhf::Prong0Phi,
                  fdhf::Prong1Phi,
                  fdhf::Prong2Phi,
                  fdhf::CandidateSelFlag,
                  fdhf::BDTBkg,
                  fdhf::BDTPrompt,
                  fdhf::BDTFD,
                  fdhf::M<fdhf::Prong0Pt, fdhf::Prong0Phi, fdhf::Prong0Eta, fdhf::Prong1Pt, fdhf::Prong1Phi, fdhf::Prong1Eta, fdhf::Prong2Pt, fdhf::Prong2Phi, fdhf::Prong2Eta>,
                  fdhf::P<fdhf::Prong0Pt, fdhf::Prong0Phi, fdhf::Prong0Eta, fdhf::Prong1Pt, fdhf::Prong1Phi, fdhf::Prong1Eta, fdhf::Prong2Pt, fdhf::Prong2Phi, fdhf::Prong2Eta>,
                  fdhf::Y<fdhf::Prong0Pt, fdhf::Prong0Phi, fdhf::Prong0Eta, fdhf::Prong1Pt, fdhf::Prong1Phi, fdhf::Prong1Eta, fdhf::Prong2Pt, fdhf::Prong2Phi, fdhf::Prong2Eta>,
                  fdhf::Eta<fdhf::Prong0Pt, fdhf::Prong0Phi, fdhf::Prong0Eta, fdhf::Prong1Pt, fdhf::Prong1Phi, fdhf::Prong1Eta, fdhf::Prong2Pt, fdhf::Prong2Phi, fdhf::Prong2Eta>,
                  fdhf::Phi<fdhf::Prong0Pt, fdhf::Prong0Phi, fdhf::Prong0Eta, fdhf::Prong1Pt, fdhf::Prong1Phi, fdhf::Prong1Eta, fdhf::Prong2Pt, fdhf::Prong2Phi, fdhf::Prong2Eta>,
                  fdhf::Pt<fdhf::Prong0Pt, fdhf::Prong0Phi, fdhf::Prong0Eta, fdhf::Prong1Pt, fdhf::Prong1Phi, fdhf::Prong1Eta, fdhf::Prong2Pt, fdhf::Prong2Phi, fdhf::Prong2Eta>);

DECLARE_SOA_TABLE(FDResultsHF, "AOD", "FDRESULTSHF", //! table to store results for HF femtoscopy
                  fdhf::CharmM,
                  fdhf::CharmPt,
                  fdhf::PtAssoc,
                  fdhf::BDTBkg,
                  fdhf::BDTPrompt,
                  fdhf::BDTFD,
                  fdhf::Kstar,
                  fdhf::KT,
                  fdhf::MT,
                  fdhf::Mult,
                  fdhf::MultPercentile,
                  fdhf::Charge,
                  fdhf::PairSign,
                  fdhf::ProcessType,
                  fdhf::FlagMc,
                  fdhf::OriginMcRec);

DECLARE_SOA_TABLE(FDHfCandMC, "AOD", "FDHFCANDMC", //! Table for reconstructed MC charm hadron candidates
                  o2::soa::Index<>,
                  fdhf::FlagMc,
                  fdhf::OriginMcRec);

DECLARE_SOA_TABLE(FDParticlesIndex, "AOD", "FDPARTICLEINDEX", //! Table track index to match associate particle with charm hadron prongs
                  o2::soa::Index<>,
                  fdhf::TrackId);

DECLARE_SOA_TABLE_STAGED(FDParticles, "FDPARTICLE",
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

DECLARE_SOA_TABLE_STAGED(FDExtParticles, "FDEXTPARTICLE",
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
                         femtodreamparticle::TPCNSigmaTr,
                         femtodreamparticle::TPCNSigmaHe,
                         femtodreamparticle::TOFNSigmaEl,
                         femtodreamparticle::TOFNSigmaPi,
                         femtodreamparticle::TOFNSigmaKa,
                         femtodreamparticle::TOFNSigmaPr,
                         femtodreamparticle::TOFNSigmaDe,
                         femtodreamparticle::TOFNSigmaTr,
                         femtodreamparticle::TOFNSigmaHe,
                         femtodreamparticle::ITSSignal,
                         femtodreamparticle::ITSNSigmaEl,
                         femtodreamparticle::ITSNSigmaPi,
                         femtodreamparticle::ITSNSigmaKa,
                         femtodreamparticle::ITSNSigmaPr,
                         femtodreamparticle::ITSNSigmaDe,
                         femtodreamparticle::ITSNSigmaTr,
                         femtodreamparticle::ITSNSigmaHe,
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
  kWrongCollision,             //! particle, that was associated wrongly to the collision
  kSecondaryDaughterLambda,    //! Daughter from a Lambda decay
  kSecondaryDaughterSigmaplus, //! Daughter from a Sigma^plus decay
  kSecondaryDaughterSigma0,    //! Daughter from a Sigma^0 decay
  kSecondaryDaughterXiMinus,   //! Daughter from a Xi^- decay
  kSecondaryDaughterXi0,       //! Daughter from a Xi^0 decay
  kElse,                       //! none of the above; (NOTE: used to catch bugs. will be removed once MC usage is properly validated)
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

DECLARE_SOA_TABLE_STAGED(FDMCParticles, "FDMCPARTICLE",
                         o2::soa::Index<>,
                         femtodreamMCparticle::PartOriginMCTruth,
                         femtodreamMCparticle::PDGMCTruth,
                         femtodreamparticle::Pt,
                         femtodreamparticle::Eta,
                         femtodreamparticle::Phi);
using FDMCParticle = FDMCParticles::iterator;

DECLARE_SOA_TABLE_STAGED(FDExtMCParticles, "FDEXTMCPARTICLE",
                         femtodreamMCparticle::MotherPDG);
using FDExtMCParticle = FDExtMCParticles::iterator;

namespace mcfdlabel
{
DECLARE_SOA_INDEX_COLUMN(FDMCParticle, fdMCParticle); //! MC particle for femtodreamparticle
} // namespace mcfdlabel
DECLARE_SOA_TABLE_STAGED(FDMCLabels, "FDMCLabel", //! Table joinable to FemtoDreamParticle containing the MC labels
                         mcfdlabel::FDMCParticleId);
namespace mcfdextlabel
{
DECLARE_SOA_INDEX_COLUMN(FDExtMCParticle, fdExtMCParticle); //! MC particle for femtodreamparticle
} // namespace mcfdextlabel
DECLARE_SOA_TABLE(FDExtMCLabels, "AOD", "FDExtMCLabel", //! Table joinable to FemtoDreamParticle containing the MC labels
                  mcfdextlabel::FDExtMCParticleId);

DECLARE_SOA_TABLE(FDHfCandMCGen, "AOD", "FDHFCANDMCGEN", //! Table for generated MC charm hadron
                  mcfdlabel::FDMCParticleId,
                  fdhf::Pt<fdhf::Prong0Pt, fdhf::Prong0Phi, fdhf::Prong1Pt, fdhf::Prong1Phi, fdhf::Prong2Pt, fdhf::Prong2Phi>,
                  fdhf::Eta<fdhf::Prong0Pt, fdhf::Prong0Phi, fdhf::Prong0Eta, fdhf::Prong1Pt, fdhf::Prong1Phi, fdhf::Prong1Eta, fdhf::Prong2Pt, fdhf::Prong2Phi, fdhf::Prong2Eta>,
                  fdhf::Phi<fdhf::Prong0Pt, fdhf::Prong0Phi, fdhf::Prong0Eta, fdhf::Prong1Pt, fdhf::Prong1Phi, fdhf::Prong1Eta, fdhf::Prong2Pt, fdhf::Prong2Phi, fdhf::Prong2Eta>,
                  fdhf::Y<fdhf::Prong0Pt, fdhf::Prong0Phi, fdhf::Prong0Eta, fdhf::Prong1Pt, fdhf::Prong1Phi, fdhf::Prong1Eta, fdhf::Prong2Pt, fdhf::Prong2Phi, fdhf::Prong2Eta>,
                  fdhf::FlagMc,
                  fdhf::OriginMcGen);

/// Hash
namespace hash
{
DECLARE_SOA_COLUMN(Bin, bin, int); //! Hash for the event mixing
} // namespace hash
DECLARE_SOA_TABLE(MixingHashes, "AOD", "HASH", hash::Bin);
using MixingHash = MixingHashes::iterator;

} // namespace o2::aod

#endif // PWGCF_DATAMODEL_FEMTODERIVED_H_
       //
