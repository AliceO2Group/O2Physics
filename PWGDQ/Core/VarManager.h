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
//
// Contact: iarsene@cern.ch, i.c.arsene@fys.uio.no
//
// Class to handle analysis variables
//

#ifndef PWGDQ_CORE_VARMANAGER_H_
#define PWGDQ_CORE_VARMANAGER_H_

#include <Math/Vector4Dfwd.h>
#include <cstdint>
#ifndef HomogeneousField
#define HomogeneousField
#endif

#include <vector>
#include <map>
#include <cmath>
#include <iostream>
#include <utility>
#include <complex>

#include <TObject.h>
#include <TString.h>
#include "TRandom.h"
#include "TH3F.h"
#include "Math/Vector4D.h"
#include "Math/Vector3D.h"
#include "Math/GenVector/Boost.h"

#include "Framework/DataTypes.h"
#include "TGeoGlobalMagField.h"
#include "Field/MagneticField.h"
#include "ReconstructionDataFormats/Track.h"
#include "ReconstructionDataFormats/Vertex.h"
#include "DCAFitter/DCAFitterN.h"
#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/TriggerAliases.h"
#include "ReconstructionDataFormats/DCA.h"
#include "DetectorsBase/Propagator.h"
#include "Common/Core/trackUtilities.h"

#include "Math/SMatrix.h"
#include "ReconstructionDataFormats/TrackFwd.h"
#include "DCAFitter/FwdDCAFitterN.h"
#include "GlobalTracking/MatchGlobalFwd.h"
#include "CommonConstants/PhysicsConstants.h"
#include "CommonConstants/LHCConstants.h"

#include "KFParticle.h"
#include "KFPTrack.h"
#include "KFPVertex.h"
#include "KFParticleBase.h"
#include "KFVertex.h"

#include "Common/Core/EventPlaneHelper.h"

using std::complex;
using std::cout;
using std::endl;

using SMatrix55 = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
using SMatrix5 = ROOT::Math::SVector<double, 5>;
using Vec3D = ROOT::Math::SVector<double, 3>;

//_________________________________________________________________________
class VarManager : public TObject
{
 public:
  // map the information contained in the objects passed to the Fill functions
  enum ObjTypes {
    // NOTE: Elements containing "Reduced" in their name refer strictly to skimmed data tables
    //       and the ones that don't refer to tables from the Framework data model or both models
    BC = BIT(0),
    Collision = BIT(1),
    CollisionCentRun2 = BIT(2),
    CollisionTimestamp = BIT(3),
    ReducedEvent = BIT(4),
    ReducedEventExtended = BIT(5),
    ReducedEventVtxCov = BIT(6),
    CollisionMC = BIT(7),
    ReducedEventMC = BIT(8),
    ReducedEventQvector = BIT(9),
    CollisionCent = BIT(10),
    CollisionMult = BIT(11),
    EventFilter = BIT(12),
    CollisionQvect = BIT(13),
    ReducedEventQvectorExtra = BIT(14),
    ReducedEventRefFlow = BIT(15),
    Zdc = BIT(16),
    ReducedZdc = BIT(17),
    CollisionMultExtra = BIT(18),
    ReducedEventMultExtra = BIT(19),
    CollisionQvectCentr = BIT(20),
    Track = BIT(0),
    TrackCov = BIT(1),
    TrackExtra = BIT(2),
    TrackPID = BIT(3),      // used for basic PID properties (needed such that we can subscribe to a minimal set of PID tables): e,pi,K,p for TPC and TOF
    TrackPIDExtra = BIT(4), // extra PID information
    TrackDCA = BIT(5),
    TrackSelection = BIT(6),
    TrackV0Bits = BIT(7),
    ReducedTrack = BIT(8),
    ReducedTrackBarrel = BIT(9),
    ReducedTrackBarrelCov = BIT(10),
    ReducedTrackBarrelPID = BIT(11),
    Muon = BIT(12),
    MuonCov = BIT(13),
    ReducedMuon = BIT(14),
    ReducedMuonExtra = BIT(15),
    ReducedMuonCov = BIT(16),
    ParticleMC = BIT(17),
    Pair = BIT(18), // TODO: check whether we really need the Pair member here
    AmbiTrack = BIT(19),
    AmbiMuon = BIT(20),
    DalitzBits = BIT(21),
    TrackTPCPID = BIT(22),
    TrackMFT = BIT(23),
    ReducedTrackCollInfo = BIT(24), // TODO: remove it once new reduced data tables are produced for dielectron with ReducedTracksBarrelInfo
    ReducedMuonCollInfo = BIT(25)   // TODO: remove it once new reduced data tables are produced for dielectron with ReducedTracksBarrelInfo
  };

  enum PairCandidateType {
    // TODO: need to agree on a scheme to incorporate all various hypotheses (e.g. e - mu, jpsi - K+, Jpsi - pipi,...)
    kDecayToEE = 0, // e.g. J/psi        -> e+ e-
    kDecayToMuMu,   // e.g. J/psi        -> mu+ mu-
    kDecayToPiPi,
    kElectronMuon,              // e.g. Electron - muon correlations
    kBcToThreeMuons,            // e.g. Bc           -> mu+ mu- mu+
    kBtoJpsiEEK,                // e.g. B+           -> e+ e- K+
    kXtoJpsiPiPi,               // e.g. X(3872)      -> J/psi pi+ pi-
    kChictoJpsiEE,              // e.g. Chi_c1      -> J/psi e+ e-
    kDstarToD0KPiPi,            // e.g. D*+ -> D0 pi+ -> K- pi+ pi+
    kTripleCandidateToEEPhoton, // e.g. chi_c   -> e+ e- photon or pi0 -> e+ e- photon
    kDecayToKPi,                // e.g. D0           -> K+ pi- or cc.
    kTripleCandidateToKPiPi,    // e.g. D+ -> K- pi+ pi+
    kTripleCandidateToPKPi,     // e.g. Lambda_c -> p K- pi+
    kNMaxCandidateTypes
  };

  enum BarrelTrackFilteringBits {
    kIsConversionLeg = 0,     // electron from conversions
    kIsK0sLeg,                // pion from K0s
    kIsLambdaLeg,             // proton or pion from Lambda
    kIsALambdaLeg,            // proton or pion from anti-Lambda
    kIsOmegaLeg,              // kaon from Omega baryon decay
    kDalitzBits = 5,          // first bit for Dalitz tagged tracks
    kBarrelUserCutsBits = 13, // first bit for user track cuts
    kIsTPCPostcalibrated = 63 // tracks were postcalibrated for the TPC PID
  };

  enum MuonTrackFilteringBits {
    kMuonUserCutsBits = 0, // first bit for user muon cuts
    kMuonIsPropagated = 7  // whether the muon was propagated already
  };

 public:
  enum Variables {
    kNothing = -1,
    // Run wise variables
    kRunNo = 0,
    kRunId,
    kRunIndex,
    kNRunWiseVariables,

    // Event wise variables
    kTimestamp,
    kTimeFromSOR, // Time since Start of Run (SOR) in minutes
    kCollisionTime,
    kCollisionTimeRes,
    kBC,
    kBCOrbit,
    kIsPhysicsSelection,
    kIsNoTFBorder,               // No time frame border
    kIsNoITSROFBorder,           // No ITS read out frame border (from event selection)
    kIsNoITSROFBorderRecomputed, // No ITS read out frame border, computed here
    kIsNoSameBunch,              // No collisions with same T0 BC
    kIsGoodZvtxFT0vsPV,          // No collisions w/ difference between z_ {PV, tracks} and z_{PV FT0A-C}
    kIsVertexITSTPC,             // At least one ITS-TPC track
    kIsVertexTOFmatched,         // At least one TOF-matched track
    kIsSel8,                     // TVX in Run3
    kIsINT7,
    kIsEMC7,
    kIsINT7inMUON,
    kIsMuonSingleLowPt7,
    kIsMuonSingleHighPt7,
    kIsMuonUnlikeLowPt7,
    kIsMuonLikeLowPt7,
    kIsCUP8,
    kIsCUP9,
    kIsMUP10,
    kIsMUP11,
    kVtxX,
    kVtxY,
    kVtxZ,
    kVtxNcontrib,
    kVtxNcontribReal,
    kVtxCovXX,
    kVtxCovXY,
    kVtxCovXZ,
    kVtxCovYY,
    kVtxCovYZ,
    kVtxCovZZ,
    kVtxChi2,
    kCentVZERO,
    kCentFT0C,
    kMultTPC,
    kMultFV0A,
    kMultFV0C,
    kMultFT0A,
    kMultFT0C,
    kMultFDDA,
    kMultFDDC,
    kMultZNA,
    kMultZNC,
    kMultTracklets,
    kMultDimuons,
    kMultNTracksHasITS,
    kMultNTracksHasTPC,
    kMultNTracksHasTOF,
    kMultNTracksHasTRD,
    kMultNTracksITSOnly,
    kMultNTracksTPCOnly,
    kMultNTracksITSTPC,
    kTrackOccupancyInTimeRange,
    kMultAllTracksTPCOnly,
    kMultAllTracksITSTPC,
    kNTPCpileupContribA,
    kNTPCpileupContribC,
    kNTPCpileupZA,
    kNTPCpileupZC,
    kNTPCtracksInPast,
    kNTPCtracksInFuture,
    kMCEventGeneratorId,
    kMCEventSubGeneratorId,
    kMCVtxX,
    kMCVtxY,
    kMCVtxZ,
    kMCEventTime,
    kMCEventWeight,
    kMCEventImpParam,
    kQ1ZNAX,
    kQ1ZNAY,
    kQ1ZNCX,
    kQ1ZNCY,
    KIntercalibZNA,
    KIntercalibZNC,
    kQ1ZNACXX,
    kQ1ZNACYY,
    kQ1ZNACYX,
    kQ1ZNACXY,
    kQ1X0A, // q-vector (e.g. from TPC) with x component (harmonic 1 and power 0), sub-event A
    kQ1Y0A, // q-vector (e.g. from TPC) with y component (harmonic 1 and power 0), sub-event A
    kQ1X0B,
    kQ1Y0B,
    kQ1X0C,
    kQ1Y0C,
    kQ2X0A,    // q-vector (e.g. from TPC) with x component (harmonic 2 and power 0), sub-event A
    kQ2Y0A,    // q-vector (e.g. from TPC) with y component (harmonic 2 and power 0), sub-event A
    kQ2X0APOS, // q-vector (e.g. from TPC) with x component (harmonic 2 and power 1), Pos. TPC
    kQ2Y0APOS, // q-vector (e.g. from TPC) with y component (harmonic 2 and power 1), Pos. TPC
    kQ2X0ANEG, // q-vector (e.g. from TPC) with x component (harmonic 2 and power 1), Neg. TPC
    kQ2Y0ANEG, // q-vector (e.g. from TPC) with y component (harmonic 2 and power 1), Neg. TPC
    kQ2X0B,
    kQ2Y0B,
    kQ2X0C,
    kQ2Y0C,
    kQ2YYAB,
    kQ2XXAB,
    kQ2XYAB,
    kQ2YXAB,
    kQ2YYAC,
    kQ2XXAC,
    kQ2XYAC,
    kQ2YXAC,
    kQ2YYBC,
    kQ2XXBC,
    kQ2XYBC,
    kQ2YXBC,
    kMultA,    // Multiplicity of the sub-event A
    kMultAPOS, // Multiplicity of the sub-event A
    kMultANEG, // Multiplicity of the sub-event A
    kMultB,
    kMultC,
    kQ3X0A, // q-vector (e.g. from TPC) with x component (harmonic 3 and power 0), sub-event A
    kQ3Y0A, // q-vector (e.g. from TPC) with y component (harmonic 3 and power 0), sub-event A
    kQ3X0B,
    kQ3Y0B,
    kQ3X0C,
    kQ3Y0C,
    kQ4X0A, // q-vector (e.g. from TPC) with x component (harmonic 4 and power 0), sub-event A
    kQ4Y0A, // q-vector (e.g. from TPC) with y component (harmonic 4 and power 0), sub-event A
    kQ4X0B,
    kQ4Y0B,
    kQ4X0C,
    kQ4Y0C,
    kR2SP_AB,
    kR2SP_AC,
    kR2SP_BC,
    kR2SP_FT0CTPCPOS,
    kR2SP_FT0CTPCNEG,
    kR2SP_FT0ATPCPOS,
    kR2SP_FT0ATPCNEG,
    kR2SP_FT0MTPCPOS,
    kR2SP_FT0MTPCNEG,
    kR2SP_FV0ATPCPOS,
    kR2SP_FV0ATPCNEG,
    kR3SP,
    kR2EP_AB,
    kR2EP_AC,
    kR2EP_BC,
    kR2EP_FT0CTPCPOS,
    kR2EP_FT0CTPCNEG,
    kR2EP_FT0ATPCPOS,
    kR2EP_FT0ATPCNEG,
    kR2EP_FT0MTPCPOS,
    kR2EP_FT0MTPCNEG,
    kR2EP_FV0ATPCPOS,
    kR2EP_FV0ATPCNEG,
    kR3EP,
    kIsDoubleGap,  // Double rapidity gap
    kIsSingleGapA, // Rapidity gap on side A
    kIsSingleGapC, // Rapidity gap on side C
    kIsSingleGap,  // Rapidity gap on either side
    kTwoEvPosZ1,   // vtx-z for collision 1 in two events correlations
    kTwoEvPosZ2,   // vtx-z for collision 2 in two events correlations
    kTwoEvPosR1,   // vtx-R for collision 1 in two events correlations
    kTwoEvPosR2,
    kTwoEvPVcontrib1, // n-contributors for collision 1 in two events correlations
    kTwoEvPVcontrib2,
    kTwoEvDeltaZ, // distance in z between collisions
    kTwoEvDeltaX, // distance in x between collisions
    kTwoEvDeltaY, // distance in y between collisions
    kTwoEvDeltaR, // distance in (x,y) plane between collisions
    kEnergyCommonZNA,
    kEnergyCommonZNC,
    kEnergyCommonZPA,
    kEnergyCommonZPC,
    kEnergyZNA1,
    kEnergyZNA2,
    kEnergyZNA3,
    kEnergyZNA4,
    kEnergyZNC1,
    kEnergyZNC2,
    kEnergyZNC3,
    kEnergyZNC4,
    kTimeZNA,
    kTimeZNC,
    kTimeZPA,
    kTimeZPC,
    kNEventWiseVariables,
    kQ2X0A1,
    kQ2X0A2,
    kQ2Y0A1,
    kQ2Y0A2,
    kU2Q2Ev1,
    kU2Q2Ev2,
    kTwoR2SP1, // Scalar product resolution of event1 for ME technique
    kTwoR2SP2, // Scalar product resolution of event2 for ME technique
    kTwoR2EP1, // Event plane resolution of event2 for ME technique
    kTwoR2EP2, // Event plane resolution of event2 for ME technique

    // Basic track/muon/pair wise variables
    kX,
    kY,
    kZ,
    kPt,
    kSignedPt,
    kInvPt,
    kEta,
    kTgl,
    kPhi,
    kP,
    kPx,
    kPy,
    kPz,
    kRap,
    kMass,
    kCharge,
    kNBasicTrackVariables,
    kUsedKF,
    kKFMass,
    kKFMassGeoTop,

    kPt1,
    kEta1,
    kPhi1,
    kCharge1,
    kPt2,
    kEta2,
    kPhi2,
    kCharge2,

    // Barrel track variables
    kPin,
    kSignedPin,
    kTOFExpMom,
    kTrackTime,
    kTrackTimeRes,
    kTrackTimeResRelative,
    kDetectorMap,
    kHasITS,
    kHasTRD,
    kHasTOF,
    kHasTPC,
    kIsGlobalTrack,
    kIsGlobalTrackSDD,
    kIsITSrefit,
    kIsSPDany,
    kIsSPDfirst,
    kIsSPDboth,
    kIsITSibAny,
    kIsITSibFirst,
    kIsITSibAll,
    kITSncls,
    kITSchi2,
    kITSlayerHit,
    kITSmeanClsSize,
    kIsTPCrefit,
    kTPCncls,
    kITSClusterMap,
    kTPCnclsCR,
    kTPCchi2,
    kTPCsignal,
    kTPCsignalRandomized,
    kTPCsignalRandomizedDelta,
    kPhiTPCOuter,
    kTrackIsInsideTPCModule,
    kTRDsignal,
    kTRDPattern,
    kTOFbeta,
    kTrackLength,
    kTrackDCAxy,
    kTrackDCAxyProng1,
    kTrackDCAxyProng2,
    kTrackDCAz,
    kTrackDCAzProng1,
    kTrackDCAzProng2,
    kTrackDCAsigXY,
    kTrackDCAsigZ,
    kTrackDCAresXY,
    kTrackDCAresZ,
    kIsGoldenChi2,
    kTrackCYY,
    kTrackCZZ,
    kTrackCSnpSnp,
    kTrackCTglTgl,
    kTrackC1Pt21Pt2,
    kTPCnSigmaEl,
    kTPCnSigmaElRandomized,
    kTPCnSigmaElRandomizedDelta,
    kTPCnSigmaMu,
    kTPCnSigmaPi,
    kTPCnSigmaPiRandomized,
    kTPCnSigmaPiRandomizedDelta,
    kTPCnSigmaKa,
    kTPCnSigmaPr,
    kTPCnSigmaEl_Corr,
    kTPCnSigmaPi_Corr,
    kTPCnSigmaKa_Corr,
    kTPCnSigmaPr_Corr,
    kTPCnSigmaPrRandomized,
    kTPCnSigmaPrRandomizedDelta,
    kTOFnSigmaEl,
    kTOFnSigmaMu,
    kTOFnSigmaPi,
    kTOFnSigmaKa,
    kTOFnSigmaPr,
    kTrackTimeResIsRange, // Gaussian or range (see Framework/DataTypes)
    kPVContributor,       // This track has contributed to the collision vertex fit (see Framework/DataTypes)
    kOrphanTrack,         // Track has no association with any collision vertex (see Framework/DataTypes)
    kIsAmbiguous,
    kIsLegFromGamma,
    kIsLegFromK0S,
    kIsLegFromLambda,
    kIsLegFromAntiLambda,
    kIsLegFromOmega,
    kIsProtonFromLambdaAndAntiLambda,
    kIsDalitzLeg,                             // Up to 8 dalitz selections
    kBarrelNAssocsInBunch = kIsDalitzLeg + 8, // number of in bunch collision associations
    kBarrelNAssocsOutOfBunch,                 // number of out of bunch collision associations
    kNBarrelTrackVariables,

    // Muon track variables
    kMuonNClusters,
    kMuonPDca,
    kMuonRAtAbsorberEnd,
    kMCHBitMap,
    kMuonChi2,
    kMuonChi2MatchMCHMID,
    kMuonChi2MatchMCHMFT,
    kMuonMatchScoreMCHMFT,
    kMuonCXX,
    kMuonCXY,
    kMuonCYY,
    kMuonCPhiX,
    kMuonCPhiY,
    kMuonCPhiPhi,
    kMuonCTglX,
    kMuonCTglY,
    kMuonCTglPhi,
    kMuonCTglTgl,
    kMuonC1Pt2X,
    kMuonC1Pt2Y,
    kMuonC1Pt2Phi,
    kMuonC1Pt2Tgl,
    kMuonC1Pt21Pt2,
    kMuonTrackType,
    kMuonDCAx,
    kMuonDCAy,
    kMuonTime,
    kMuonTimeRes,
    kMftNClusters,
    kMftClusterSize,
    kMftMeanClusterSize,
    kMuonNAssocsInBunch,
    kMuonNAssocsOutOfBunch,
    kNMuonTrackVariables,

    // MC particle variables
    kMCPdgCode,
    kMCParticleWeight,
    kMCPx,
    kMCPy,
    kMCPz,
    kMCE,
    kMCVx,
    kMCVy,
    kMCVz,
    kMCPt,
    kMCPhi,
    kMCEta,
    kMCY,
    kMCParticleGeneratorId,
    kNMCParticleVariables,

    // MC mother particle variables
    kMCMotherPdgCode,

    // Pair variables
    kCandidateId,
    kPairType,
    kVertexingLxy,
    kVertexingLxyErr,
    kVertexingPseudoCTau,
    kVertexingLxyz,
    kVertexingLxyzErr,
    kVertexingLz,
    kVertexingLzErr,
    kVertexingTauxy,
    kVertexingTauxyErr,
    kVertexingLzProjected,
    kVertexingLxyProjected,
    kVertexingLxyzProjected,
    kVertexingTauzProjected,
    kVertexingTauxyProjected,
    kVertexingTauxyProjectedNs,
    kVertexingTauxyzProjected,
    kVertexingTauz,
    kVertexingTauzErr,
    kVertexingPz,
    kVertexingSV,
    kVertexingProcCode,
    kVertexingChi2PCA,
    kCosThetaHE,
    kCosThetaCS,
    kPhiHE,
    kPhiCS,
    kDeltaPhiPair2,
    kPsiPair,
    kDeltaPhiPair,
    kOpeningAngle,
    kQuadDCAabsXY,
    kQuadDCAsigXY,
    kQuadDCAabsZ,
    kQuadDCAsigZ,
    kQuadDCAsigXYZ,
    kSignQuadDCAsigXY,
    kCosPointingAngle,
    kImpParXYJpsi,
    kImpParXYK,
    kDCATrackProd,
    kDCATrackVtxProd,
    kU2Q2,
    kU3Q3,
    kQ42XA,
    kQ42YA,
    kQ23XA,
    kQ23YA,
    kS11A,
    kS12A,
    kS13A,
    kS31A,
    kM11REF,
    kM01POI,
    kM1111REF,
    kM0111POI,
    kCORR2REF,
    kCORR2POI,
    kCORR4REF,
    kCORR4POI,
    kM11REFoverMp,
    kM01POIoverMp,
    kM1111REFoverMp,
    kM0111POIoverMp,
    kCORR2POIMp,
    kCORR4POIMp,
    kPsi2A,
    kPsi2APOS,
    kPsi2ANEG,
    kPsi2B,
    kPsi2C,
    kCos2DeltaPhi,
    kCos2DeltaPhiMu1, // cos(phi - phi1) for muon1
    kCos2DeltaPhiMu2, ////cos(phi - phi2) for muon2
    kCos3DeltaPhi,
    kDeltaPtotTracks,
    kVertexingLxyOverErr,
    kVertexingLzOverErr,
    kVertexingLxyzOverErr,
    kKFTrack0DCAxyz,
    kKFTrack1DCAxyz,
    kKFTracksDCAxyzMax,
    kKFDCAxyzBetweenProngs,
    kKFTrack0DCAxy,
    kKFTrack1DCAxy,
    kKFTracksDCAxyMax,
    kKFDCAxyBetweenProngs,
    kKFTrack0DeviationFromPV,
    kKFTrack1DeviationFromPV,
    kKFTrack0DeviationxyFromPV,
    kKFTrack1DeviationxyFromPV,
    kKFChi2OverNDFGeo,
    kKFNContributorsPV,
    kKFCosPA,
    kKFChi2OverNDFGeoTop,
    kKFJpsiDCAxyz,
    kKFJpsiDCAxy,
    kKFPairDeviationFromPV,
    kKFPairDeviationxyFromPV,
    kNPairVariables,

    // Candidate-track correlation variables
    kPairMass,
    kPairMassDau,
    kMassDau,
    kPairPt,
    kPairPtDau,
    kPairEta,
    kPairPhi,
    kPairPhiv,
    kDeltaEta,
    kDeltaPhi,
    kDeltaPhiSym,
    kNCorrelationVariables,

    // Dilepton-track-track variables
    kQuadMass,
    kQuadDefaultDileptonMass,
    kQuadPt,
    kQuadEta,
    kQuadPhi,
    kCosthetaDileptonDitrack,
    kDitrackMass,
    kDitrackPt,
    kQ,
    kDeltaR1,
    kDeltaR2,

    // DQ-HF correlation variables
    kMassCharmHadron,
    kPtCharmHadron,
    kRapCharmHadron,
    kPhiCharmHadron,
    kBdtCharmHadron,

    // Index used to scan bit maps
    kBitMapIndex,

    // deltaMass = kPairMass - kPairMassDau
    kDeltaMass,
    // deltaMass_jpsi = kPairMass - kPairMassDau +3.096900
    kDeltaMass_jpsi,

    kNVars
  }; // end of Variables enumeration

  enum CalibObjects {
    kTPCElectronMean = 0,
    kTPCElectronSigma,
    kTPCPionMean,
    kTPCPionSigma,
    kTPCKaonMean,
    kTPCKaonSigma,
    kTPCProtonMean,
    kTPCProtonSigma,
    kNCalibObjects
  };

  enum DileptonCharmHadronTypes {
    kJPsi = 0,
    kD0ToPiK,
    kD0barToKPi
  };

  enum EventFilters {
    kDoubleGap = 0,
    kSingleGapA,
    kSingleGapC
  };

  enum MuonExtrapolation {
    // Index used to set different options for Muon propagation
    kToVertex = 0, // propagtion to vertex by default
    kToDCA,
    kToRabs
  };

  static TString fgVariableNames[kNVars]; // variable names
  static TString fgVariableUnits[kNVars]; // variable units
  static void SetDefaultVarNames();

  static void SetUseVariable(int var)
  {
    if (var >= 0 && var < kNVars) {
      fgUsedVars[var] = kTRUE;
    }
    SetVariableDependencies();
  }
  static void SetUseVars(const bool* usedVars)
  {
    for (int i = 0; i < kNVars; ++i) {
      if (usedVars[i]) {
        fgUsedVars[i] = true; // overwrite only the variables that are being used since there are more channels to modify the used variables array, independently
      }
    }
    SetVariableDependencies();
  }
  static void SetUseVars(const std::vector<int> usedVars)
  {
    for (auto& var : usedVars) {
      fgUsedVars[var] = true;
    }
  }
  static bool GetUsedVar(int var)
  {
    if (var >= 0 && var < kNVars) {
      return fgUsedVars[var];
    }
    return false;
  }

  static void SetRunNumbers(int n, int* runs);
  static void SetRunNumbers(std::vector<int> runs);
  static float GetRunIndex(double);
  static void SetDummyRunlist(int InitRunnumber);
  static int GetDummyFirst();
  static int GetDummyLast();
  static int GetDummyNRuns()
  {
    if (fgRunMap.size() == 0) {
      return 101;
    } else {
      return fgRunMap.size();
    }
  }
  static int GetNRuns()
  {
    return fgRunMap.size();
  }
  static TString GetRunStr()
  {
    return fgRunStr;
  }

  // Setup the collision system
  static void SetCollisionSystem(TString system, float energy);

  static void SetMagneticField(float magField)
  {
    fgMagField = magField;
  }

  // Setup the 2 prong KFParticle
  static void SetupTwoProngKFParticle(float magField)
  {
    KFParticle::SetField(magField);
    fgUsedKF = true;
  }
  // Setup magnetic field for muon propagation
  static void SetupMuonMagField()
  {
    o2::mch::TrackExtrap::setField();
  }

  // Setup the 2 prong DCAFitterN
  static void SetupTwoProngDCAFitter(float magField, bool propagateToPCA, float maxR, float maxDZIni, float minParamChange, float minRelChi2Change, bool useAbsDCA)
  {
    fgFitterTwoProngBarrel.setBz(magField);
    fgFitterTwoProngBarrel.setPropagateToPCA(propagateToPCA);
    fgFitterTwoProngBarrel.setMaxR(maxR);
    fgFitterTwoProngBarrel.setMaxDZIni(maxDZIni);
    fgFitterTwoProngBarrel.setMinParamChange(minParamChange);
    fgFitterTwoProngBarrel.setMinRelChi2Change(minRelChi2Change);
    fgFitterTwoProngBarrel.setUseAbsDCA(useAbsDCA);
    fgUsedKF = false;
  }

  // Setup the 2 prong FwdDCAFitterN
  static void SetupTwoProngFwdDCAFitter(float magField, bool propagateToPCA, float maxR, float minParamChange, float minRelChi2Change, bool useAbsDCA)
  {
    fgFitterTwoProngFwd.setBz(magField);
    fgFitterTwoProngFwd.setPropagateToPCA(propagateToPCA);
    fgFitterTwoProngFwd.setMaxR(maxR);
    fgFitterTwoProngFwd.setMinParamChange(minParamChange);
    fgFitterTwoProngFwd.setMinRelChi2Change(minRelChi2Change);
    fgFitterTwoProngFwd.setUseAbsDCA(useAbsDCA);
    fgUsedKF = false;
  }
  // Use MatLayerCylSet to correct MCS in fwdtrack propagation
  static void SetupMatLUTFwdDCAFitter(o2::base::MatLayerCylSet* m)
  {
    fgFitterTwoProngFwd.setTGeoMat(false);
    fgFitterTwoProngFwd.setMatLUT(m);
  }
  // Use GeometryManager to correct MCS in fwdtrack propagation
  static void SetupTGeoFwdDCAFitter()
  {
    fgFitterTwoProngFwd.setTGeoMat(true);
  }
  // No material budget in fwdtrack propagation
  static void SetupFwdDCAFitterNoCorr()
  {
    fgFitterTwoProngFwd.setTGeoMat(false);
  }
  // Setup the 3 prong KFParticle
  static void SetupThreeProngKFParticle(float magField)
  {
    KFParticle::SetField(magField);
    fgUsedKF = true;
  }

  // Setup the 3 prong DCAFitterN
  static void SetupThreeProngDCAFitter(float magField, bool propagateToPCA, float maxR, float /*maxDZIni*/, float minParamChange, float minRelChi2Change, bool useAbsDCA)
  {
    fgFitterThreeProngBarrel.setBz(magField);
    fgFitterThreeProngBarrel.setPropagateToPCA(propagateToPCA);
    fgFitterThreeProngBarrel.setMaxR(maxR);
    fgFitterThreeProngBarrel.setMinParamChange(minParamChange);
    fgFitterThreeProngBarrel.setMinRelChi2Change(minRelChi2Change);
    fgFitterThreeProngBarrel.setUseAbsDCA(useAbsDCA);
    fgUsedKF = false;
  }

  static auto getEventPlane(int harm, float qnxa, float qnya)
  {
    // Compute event plane angle from qn vector components for the sub-event A
    return (1.0 / harm) * TMath::ATan2(qnya, qnxa);
  };

  static float getDeltaPsiInRange(float psi1, float psi2, float harmonic)
  {
    float deltaPsi = psi1 - psi2;
    if (std::abs(deltaPsi) > o2::constants::math::PI / harmonic) {
      if (deltaPsi > 0.) {
        deltaPsi -= o2::constants::math::TwoPI / harmonic;
      } else {
        deltaPsi += o2::constants::math::TwoPI / harmonic;
      }
    }
    return deltaPsi;
  }

  template <typename T, typename C>
  static o2::dataformats::GlobalFwdTrack PropagateMuon(const T& muon, const C& collision, int endPoint = kToVertex);
  template <uint32_t fillMap, typename T, typename C>
  static void FillMuonPDca(const T& muon, const C& collision, float* values = nullptr);
  template <uint32_t fillMap, typename T, typename C>
  static void FillPropagateMuon(const T& muon, const C& collision, float* values = nullptr);
  template <typename T>
  static void FillBC(T const& bc, float* values = nullptr);
  template <uint32_t fillMap, typename T>
  static void FillEvent(T const& event, float* values = nullptr);
  template <uint32_t fillMap, typename TEvent, typename TAssoc, typename TTracks>
  static void FillEventTrackEstimators(TEvent const& collision, TAssoc const& groupedTrackIndices, TTracks const& tracks, float* values = nullptr);
  template <typename T>
  static void FillTwoEvents(T const& event1, T const& event2, float* values = nullptr);
  template <uint32_t fillMap, typename T1, typename T2>
  static void FillTwoMixEvents(T1 const& event1, T1 const& event2, T2 const& tracks1, T2 const& tracks2, float* values = nullptr);
  template <uint32_t fillMap, typename T>
  static void FillTrack(T const& track, float* values = nullptr);
  template <uint32_t fillMap, typename T>
  static void FillPhoton(T const& photon, float* values = nullptr);
  template <uint32_t fillMap, typename T, typename C>
  static void FillTrackCollision(T const& track, C const& collision, float* values = nullptr);
  template <uint32_t fillMap, typename T, typename C, typename M, typename P>
  static void FillTrackCollisionMatCorr(T const& track, C const& collision, M const& materialCorr, P const& propagator, float* values = nullptr);
  template <typename U, typename T>
  static void FillTrackMC(const U& mcStack, T const& track, float* values = nullptr);
  template <uint32_t fillMap, typename T1, typename T2, typename C>
  static void FillPairPropagateMuon(T1 const& muon1, T2 const& muon2, const C& collision, float* values = nullptr);
  template <uint32_t fillMap, typename T1, typename T2, typename C>
  static void FillGlobalMuonRefit(T1 const& muontrack, T2 const& mfttrack, const C& collision, float* values = nullptr);
  template <int pairType, uint32_t fillMap, typename T1, typename T2>
  static void FillPair(T1 const& t1, T2 const& t2, float* values = nullptr);
  template <int pairType, uint32_t fillMap, typename C, typename T1, typename T2>
  static void FillPairCollision(C const& collision, T1 const& t1, T2 const& t2, float* values = nullptr);
  template <int pairType, uint32_t fillMap, typename C, typename T1, typename T2, typename M, typename P>
  static void FillPairCollisionMatCorr(C const& collision, T1 const& t1, T2 const& t2, M const& materialCorr, P const& propagator, float* values = nullptr);
  template <typename T1, typename T2, typename T3>
  static void FillTriple(T1 const& t1, T2 const& t2, T3 const& t3, float* values = nullptr, PairCandidateType pairType = kTripleCandidateToEEPhoton);
  template <int pairType, typename T1, typename T2>
  static void FillPairME(T1 const& t1, T2 const& t2, float* values = nullptr);
  template <typename T1, typename T2>
  static void FillPairMC(T1 const& t1, T2 const& t2, float* values = nullptr, PairCandidateType pairType = kDecayToEE);
  template <typename T1, typename T2, typename T3>
  static void FillTripleMC(T1 const& t1, T2 const& t2, T3 const& t3, float* values = nullptr, PairCandidateType pairType = kTripleCandidateToEEPhoton);
  template <int pairType, uint32_t collFillMap, uint32_t fillMap, typename C, typename T>
  static void FillPairVertexing(C const& collision, T const& t1, T const& t2, bool propToSV = false, float* values = nullptr);
  template <uint32_t collFillMap, uint32_t fillMap, typename C, typename T>
  static void FillTripletVertexing(C const& collision, T const& t1, T const& t2, T const& t3, PairCandidateType tripletType, float* values = nullptr);
  template <int candidateType, uint32_t collFillMap, uint32_t fillMap, typename C, typename T1>
  static void FillDileptonTrackVertexing(C const& collision, T1 const& lepton1, T1 const& lepton2, T1 const& track, float* values);
  template <typename T1, typename T2>
  static void FillDileptonHadron(T1 const& dilepton, T2 const& hadron, float* values = nullptr, float hadronMass = 0.0f);
  template <typename T1, typename T2>
  static void FillDileptonPhoton(T1 const& dilepton, T2 const& photon, float* values = nullptr);
  template <typename T>
  static void FillHadron(T const& hadron, float* values = nullptr, float hadronMass = 0.0f);
  template <int partType, typename Cand, typename H, typename T>
  static void FillSingleDileptonCharmHadron(Cand const& candidate, H hfHelper, T& bdtScoreCharmHad, float* values = nullptr);
  template <int partTypeCharmHad, typename DQ, typename HF, typename H, typename T>
  static void FillDileptonCharmHadron(DQ const& dilepton, HF const& charmHadron, H hfHelper, T& bdtScoreCharmHad, float* values = nullptr);
  template <typename C, typename A>
  static void FillQVectorFromGFW(C const& collision, A const& compA11, A const& compB11, A const& compC11, A const& compA21, A const& compB21, A const& compC21, A const& compA31, A const& compB31, A const& compC31, A const& compA41, A const& compB41, A const& compC41, A const& compA23, A const& compA42, float S10A = 1.0, float S10B = 1.0, float S10C = 1.0, float S11A = 1.0, float S11B = 1.0, float S11C = 1.0, float S12A = 1.0, float S13A = 1.0, float S14A = 1.0, float S21A = 1.0, float S22A = 1.0, float S31A = 1.0, float S41A = 1.0, float* values = nullptr);
  template <typename C>
  static void FillQVectorFromCentralFW(C const& collision, float* values = nullptr);
  template <typename C>
  static void FillNewQVectorFromCentralFW(C const& collision, float* values = nullptr);
  template <typename C>
  static void FillSpectatorPlane(C const& collision, float* values = nullptr);
  template <uint32_t fillMap, int pairType, typename T1, typename T2>
  static void FillPairVn(T1 const& t1, T2 const& t2, float* values = nullptr);
  template <int candidateType, typename T1, typename T2, typename T3>
  static void FillDileptonTrackTrack(T1 const& dilepton, T2 const& hadron1, T3 const& hadron2, float* values = nullptr);
  template <typename T>
  static void FillZDC(const T& zdc, float* values = nullptr);

  static void SetCalibrationObject(CalibObjects calib, TObject* obj)
  {
    fgCalibs[calib] = obj;
    // Check whether all the needed objects for TPC postcalibration are available
    if (fgCalibs.find(kTPCElectronMean) != fgCalibs.end() && fgCalibs.find(kTPCElectronSigma) != fgCalibs.end()) {
      fgRunTPCPostCalibration[0] = true;
    }
    if (fgCalibs.find(kTPCPionMean) != fgCalibs.end() && fgCalibs.find(kTPCPionSigma) != fgCalibs.end()) {
      fgRunTPCPostCalibration[1] = true;
    }
    if (fgCalibs.find(kTPCKaonMean) != fgCalibs.end() && fgCalibs.find(kTPCKaonSigma) != fgCalibs.end()) {
      fgRunTPCPostCalibration[2] = true;
    }
    if (fgCalibs.find(kTPCProtonMean) != fgCalibs.end() && fgCalibs.find(kTPCProtonSigma) != fgCalibs.end()) {
      fgRunTPCPostCalibration[3] = true;
    }
  }
  static TObject* GetCalibrationObject(CalibObjects calib)
  {
    auto obj = fgCalibs.find(calib);
    if (obj == fgCalibs.end()) {
      return 0x0;
    } else {
      return obj->second;
    }
  }
  static void SetTPCInterSectorBoundary(float boundarySize)
  {
    fgTPCInterSectorBoundary = boundarySize;
  }
  static void SetITSROFBorderselection(int bias, int length, int marginLow, int marginHigh)
  {
    fgITSROFbias = bias;
    fgITSROFlength = length;
    fgITSROFBorderMarginLow = marginLow;
    fgITSROFBorderMarginHigh = marginHigh;
  }

  static void SetSORandEOR(uint64_t sor, uint64_t eor)
  {
    fgSOR = sor;
    fgEOR = eor;
  }

 public:
  VarManager();
  ~VarManager() override;

  static float fgValues[kNVars]; // array holding all variables computed during analysis
  static void ResetValues(int startValue = 0, int endValue = kNVars, float* values = nullptr);

 private:
  static bool fgUsedVars[kNVars]; // holds flags for when the corresponding variable is needed (e.g., in the histogram manager, in cuts, mixing handler, etc.)
  static bool fgUsedKF;
  static void SetVariableDependencies(); // toggle those variables on which other used variables might depend

  static float fgMagField;
  static std::map<int, int> fgRunMap;     // map of runs to be used in histogram axes
  static TString fgRunStr;                // semi-colon separated list of runs, to be used for histogram axis labels
  static std::vector<int> fgRunList;      // vector of runs, to be used for histogram axis
  static float fgCenterOfMassEnergy;      // collision energy
  static float fgMassofCollidingParticle; // mass of the colliding particle
  static float fgTPCInterSectorBoundary;  // TPC inter-sector border size at the TPC outer radius, in cm
  static int fgITSROFbias;                // ITS ROF bias (from ALPIDE parameters)
  static int fgITSROFlength;              // ITS ROF length (from ALPIDE parameters)
  static int fgITSROFBorderMarginLow;     // ITS ROF border low margin
  static int fgITSROFBorderMarginHigh;    // ITS ROF border high margin
  static uint64_t fgSOR;                  // Timestamp for start of run
  static uint64_t fgEOR;                  // Timestamp for end of run

  static void FillEventDerived(float* values = nullptr);
  static void FillTrackDerived(float* values = nullptr);
  template <typename T, typename U, typename V>
  static auto getRotatedCovMatrixXX(const T& matrix, U phi, V theta);
  template <typename T>
  static KFPTrack createKFPTrackFromTrack(const T& track);
  template <typename T>
  static KFPTrack createKFPFwdTrackFromFwdTrack(const T& muon);
  template <typename T>
  static KFPVertex createKFPVertexFromCollision(const T& collision);
  static float calculateCosPA(KFParticle kfp, KFParticle PV);
  template <int pairType, typename T1, typename T2>
  static float calculatePhiV(const T1& t1, const T2& t2);

  static o2::vertexing::DCAFitterN<2> fgFitterTwoProngBarrel;
  static o2::vertexing::DCAFitterN<3> fgFitterThreeProngBarrel;
  static o2::vertexing::FwdDCAFitterN<2> fgFitterTwoProngFwd;
  static o2::vertexing::FwdDCAFitterN<3> fgFitterThreeProngFwd;
  static o2::globaltracking::MatchGlobalFwd mMatching;

  static std::map<CalibObjects, TObject*> fgCalibs; // map of calibration histograms
  static bool fgRunTPCPostCalibration[4];           // 0-electron, 1-pion, 2-kaon, 3-proton

  VarManager& operator=(const VarManager& c);
  VarManager(const VarManager& c);

  ClassDef(VarManager, 3);
};

template <typename T, typename U, typename V>
auto VarManager::getRotatedCovMatrixXX(const T& matrix, U phi, V theta)
{
  //
  auto cp = std::cos(phi);
  auto sp = std::sin(phi);
  auto ct = std::cos(theta);
  auto st = std::sin(theta);
  return matrix[0] * cp * cp * ct * ct        // covXX
         + matrix[1] * 2. * cp * sp * ct * ct // covXY
         + matrix[2] * sp * sp * ct * ct      // covYY
         + matrix[3] * 2. * cp * ct * st      // covXZ
         + matrix[4] * 2. * sp * ct * st      // covYZ
         + matrix[5] * st * st;               // covZZ
}

template <typename T>
KFPTrack VarManager::createKFPTrackFromTrack(const T& track)
{
  std::array<float, 5> trackpars = {track.y(), track.z(), track.snp(), track.tgl(), track.signed1Pt()};
  std::array<float, 15> trackcovs = {track.cYY(), track.cZY(), track.cZZ(), track.cSnpY(), track.cSnpZ(),
                                     track.cSnpSnp(), track.cTglY(), track.cTglZ(), track.cTglSnp(), track.cTglTgl(),
                                     track.c1PtY(), track.c1PtZ(), track.c1PtSnp(), track.c1PtTgl(), track.c1Pt21Pt2()};
  o2::track::TrackParametrizationWithError trackparCov{track.x(), track.alpha(), std::move(trackpars), std::move(trackcovs)};
  std::array<float, 3> trkpos_par;
  std::array<float, 3> trkmom_par;
  std::array<float, 21> trk_cov;
  trackparCov.getXYZGlo(trkpos_par);
  trackparCov.getPxPyPzGlo(trkmom_par);
  trackparCov.getCovXYZPxPyPzGlo(trk_cov);
  float trkpar_KF[6] = {trkpos_par[0], trkpos_par[1], trkpos_par[2],
                        trkmom_par[0], trkmom_par[1], trkmom_par[2]};
  float trkcov_KF[21];
  for (int i = 0; i < 21; i++) {
    trkcov_KF[i] = trk_cov[i];
  }
  KFPTrack kfpTrack;
  kfpTrack.SetParameters(trkpar_KF);
  kfpTrack.SetCovarianceMatrix(trkcov_KF);
  kfpTrack.SetCharge(track.sign());
  kfpTrack.SetNDF(track.tpcNClsFound() - 5);
  kfpTrack.SetChi2(track.tpcChi2NCl() * track.tpcNClsFound());
  return kfpTrack;
}

template <typename T>
KFPTrack VarManager::createKFPFwdTrackFromFwdTrack(const T& muon)
{
  double chi2 = muon.chi2();
  SMatrix5 tpars(muon.x(), muon.y(), muon.phi(), muon.tgl(), muon.signed1Pt());
  std::vector<double> v1{muon.cXX(), muon.cXY(), muon.cYY(), muon.cPhiX(), muon.cPhiY(),
                         muon.cPhiPhi(), muon.cTglX(), muon.cTglY(), muon.cTglPhi(), muon.cTglTgl(),
                         muon.c1PtX(), muon.c1PtY(), muon.c1PtPhi(), muon.c1PtTgl(), muon.c1Pt21Pt2()};
  SMatrix55 tcovs(v1.begin(), v1.end());
  o2::track::TrackParCovFwd trackparCov{muon.z(), tpars, tcovs, chi2};

  std::array<float, 21> trk_cov;
  trackparCov.getCovXYZPxPyPzGlo(trk_cov);
  double trkpar_KF[6] = {trackparCov.getX(), trackparCov.getY(), trackparCov.getZ(),
                         trackparCov.getPx(), trackparCov.getPy(), trackparCov.getPz()};
  float trkcov_KF[21];
  for (int i = 0; i < 21; i++) {
    trkcov_KF[i] = trk_cov[i];
  }
  KFPTrack kfpTrack;
  kfpTrack.SetParameters(trkpar_KF);
  kfpTrack.SetCovarianceMatrix(trkcov_KF);
  kfpTrack.SetCharge(muon.sign());
  kfpTrack.SetNDF(muon.nClusters() - 5);
  kfpTrack.SetChi2(chi2);
  return kfpTrack;
}

template <typename T>
KFPVertex VarManager::createKFPVertexFromCollision(const T& collision)
{
  KFPVertex kfpVertex;
  kfpVertex.SetXYZ(collision.posX(), collision.posY(), collision.posZ());
  kfpVertex.SetCovarianceMatrix(collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ());
  kfpVertex.SetChi2(collision.chi2());
  kfpVertex.SetNDF(2 * collision.numContrib() - 3);
  kfpVertex.SetNContributors(collision.numContrib());
  return kfpVertex;
}

template <typename T, typename C>
o2::dataformats::GlobalFwdTrack VarManager::PropagateMuon(const T& muon, const C& collision, const int endPoint)
{
  double chi2 = muon.chi2();
  SMatrix5 tpars(muon.x(), muon.y(), muon.phi(), muon.tgl(), muon.signed1Pt());
  std::vector<double> v1{muon.cXX(), muon.cXY(), muon.cYY(), muon.cPhiX(), muon.cPhiY(),
                         muon.cPhiPhi(), muon.cTglX(), muon.cTglY(), muon.cTglPhi(), muon.cTglTgl(),
                         muon.c1PtX(), muon.c1PtY(), muon.c1PtPhi(), muon.c1PtTgl(), muon.c1Pt21Pt2()};
  SMatrix55 tcovs(v1.begin(), v1.end());
  o2::track::TrackParCovFwd fwdtrack{muon.z(), tpars, tcovs, chi2};
  o2::dataformats::GlobalFwdTrack propmuon;
  if (static_cast<int>(muon.trackType()) > 2) {
    o2::dataformats::GlobalFwdTrack track;
    track.setParameters(tpars);
    track.setZ(fwdtrack.getZ());
    track.setCovariances(tcovs);
    auto mchTrack = mMatching.FwdtoMCH(track);

    if (endPoint == kToVertex) {
      o2::mch::TrackExtrap::extrapToVertex(mchTrack, collision.posX(), collision.posY(), collision.posZ(), collision.covXX(), collision.covYY());
    }
    if (endPoint == kToDCA) {
      o2::mch::TrackExtrap::extrapToVertexWithoutBranson(mchTrack, collision.posZ());
    }
    if (endPoint == kToRabs) {
      o2::mch::TrackExtrap::extrapToZ(mchTrack, -505.);
    }

    auto proptrack = mMatching.MCHtoFwd(mchTrack);
    propmuon.setParameters(proptrack.getParameters());
    propmuon.setZ(proptrack.getZ());
    propmuon.setCovariances(proptrack.getCovariances());

  } else if (static_cast<int>(muon.trackType()) < 2) {
    double centerMFT[3] = {0, 0, -61.4};
    o2::field::MagneticField* field = static_cast<o2::field::MagneticField*>(TGeoGlobalMagField::Instance()->GetField());
    auto Bz = field->getBz(centerMFT); // Get field at centre of MFT
    auto geoMan = o2::base::GeometryManager::meanMaterialBudget(muon.x(), muon.y(), muon.z(), collision.posX(), collision.posY(), collision.posZ());
    auto x2x0 = static_cast<float>(geoMan.meanX2X0);
    fwdtrack.propagateToVtxhelixWithMCS(collision.posZ(), {collision.posX(), collision.posY()}, {collision.covXX(), collision.covYY()}, Bz, x2x0);
    propmuon.setParameters(fwdtrack.getParameters());
    propmuon.setZ(fwdtrack.getZ());
    propmuon.setCovariances(fwdtrack.getCovariances());
  }
  return propmuon;
}

template <uint32_t fillMap, typename T, typename C>
void VarManager::FillMuonPDca(const T& muon, const C& collision, float* values)
{
  if (!values) {
    values = fgValues;
  }

  if constexpr ((fillMap & MuonCov) > 0 || (fillMap & ReducedMuonCov) > 0) {

    o2::dataformats::GlobalFwdTrack propmuon = PropagateMuon(muon, collision);
    o2::dataformats::GlobalFwdTrack propmuonAtDCA = PropagateMuon(muon, collision, kToDCA);

    float dcaX = (propmuonAtDCA.getX() - collision.posX());
    float dcaY = (propmuonAtDCA.getY() - collision.posY());
    float dcaXY = std::sqrt(dcaX * dcaX + dcaY * dcaY);
    values[kMuonPDca] = muon.p() * dcaXY;
  }
}

template <uint32_t fillMap, typename T, typename C>
void VarManager::FillPropagateMuon(const T& muon, const C& collision, float* values)
{
  if (!values) {
    values = fgValues;
  }

  if constexpr ((fillMap & ReducedMuonCov) > 0) {
    if (muon.filteringFlags() & (uint8_t(1) << VarManager::kMuonIsPropagated)) { // the muon is already propagated, so nothing to do
      return;
    }
  }

  if constexpr ((fillMap & MuonCov) > 0 || (fillMap & ReducedMuonCov) > 0) {
    o2::dataformats::GlobalFwdTrack propmuon = PropagateMuon(muon, collision);
    values[kPt] = propmuon.getPt();
    values[kX] = propmuon.getX();
    values[kY] = propmuon.getY();
    values[kZ] = propmuon.getZ();
    values[kEta] = propmuon.getEta();
    values[kTgl] = propmuon.getTgl();
    values[kPhi] = propmuon.getPhi();

    // Redo propagation only for muon tracks
    // propagation of MFT tracks alredy done in fwdtrack-extention task
    if (static_cast<int>(muon.trackType()) > 2) {
      o2::dataformats::GlobalFwdTrack propmuonAtDCA = PropagateMuon(muon, collision, kToDCA);
      o2::dataformats::GlobalFwdTrack propmuonAtRabs = PropagateMuon(muon, collision, kToRabs);
      float dcaX = (propmuonAtDCA.getX() - collision.posX());
      float dcaY = (propmuonAtDCA.getY() - collision.posY());
      values[kMuonDCAx] = dcaX;
      values[kMuonDCAy] = dcaY;
      double xAbs = propmuonAtRabs.getX();
      double yAbs = propmuonAtRabs.getY();
      values[kMuonRAtAbsorberEnd] = std::sqrt(xAbs * xAbs + yAbs * yAbs);
    }

    SMatrix55 cov = propmuon.getCovariances();
    values[kMuonCXX] = cov(0, 0);
    values[kMuonCXY] = cov(1, 0);
    values[kMuonCYY] = cov(1, 1);
    values[kMuonCPhiX] = cov(2, 0);
    values[kMuonCPhiY] = cov(2, 1);
    values[kMuonCPhiPhi] = cov(2, 2);
    values[kMuonCTglX] = cov(3, 0);
    values[kMuonCTglY] = cov(3, 1);
    values[kMuonCTglPhi] = cov(3, 2);
    values[kMuonCTglTgl] = cov(3, 3);
    values[kMuonC1Pt2X] = cov(4, 0);
    values[kMuonC1Pt2Y] = cov(4, 1);
    values[kMuonC1Pt2Phi] = cov(4, 2);
    values[kMuonC1Pt2Tgl] = cov(4, 3);
    values[kMuonC1Pt21Pt2] = cov(4, 4);
  }
}

template <uint32_t fillMap, typename T1, typename T2, typename C>
void VarManager::FillGlobalMuonRefit(T1 const& muontrack, T2 const& mfttrack, const C& collision, float* values)
{
  if (!values) {
    values = fgValues;
  }
  if constexpr ((fillMap & MuonCov) > 0 || (fillMap & ReducedMuonCov) > 0) {
    o2::dataformats::GlobalFwdTrack propmuon = PropagateMuon(muontrack, collision);
    double px = propmuon.getP() * sin(M_PI / 2 - atan(mfttrack.tgl())) * cos(mfttrack.phi());
    double py = propmuon.getP() * sin(M_PI / 2 - atan(mfttrack.tgl())) * sin(mfttrack.phi());
    double pz = propmuon.getP() * cos(M_PI / 2 - atan(mfttrack.tgl()));
    double pt = std::sqrt(std::pow(px, 2) + std::pow(py, 2));
    values[kX] = mfttrack.x();
    values[kY] = mfttrack.y();
    values[kZ] = mfttrack.z();
    values[kTgl] = mfttrack.tgl();
    values[kPt] = pt;
    values[kPz] = pz;
    values[kEta] = mfttrack.eta();
    values[kPhi] = mfttrack.phi();
  }
}

template <typename T>
void VarManager::FillBC(T const& bc, float* values)
{
  if (!values) {
    values = fgValues;
  }
  values[kRunNo] = bc.runNumber();
  values[kBC] = bc.globalBC();
  values[kBCOrbit] = bc.globalBC() % o2::constants::lhc::LHCMaxBunches;
  values[kTimestamp] = bc.timestamp();
  values[kTimeFromSOR] = (fgSOR > 0 ? (bc.timestamp() - fgSOR) / 60000. : -1.0);
  values[kRunIndex] = GetRunIndex(bc.runNumber());
}

template <uint32_t fillMap, typename T>
void VarManager::FillEvent(T const& event, float* values)
{
  if (!values) {
    values = fgValues;
  }

  if constexpr ((fillMap & CollisionTimestamp) > 0) {
    values[kTimestamp] = event.timestamp();
  }

  if constexpr ((fillMap & Collision) > 0) {
    // TODO: trigger info from the event selection requires a separate flag
    //       so that it can be switched off independently of the rest of Collision variables (e.g. if event selection is not available)

    if (fgUsedVars[kIsNoITSROFBorder]) {
      values[kIsNoITSROFBorder] = event.selection_bit(o2::aod::evsel::kNoITSROFrameBorder);
    }
    if (fgUsedVars[kIsNoTFBorder]) {
      values[kIsNoTFBorder] = event.selection_bit(o2::aod::evsel::kNoTimeFrameBorder);
    }
    if (fgUsedVars[kIsNoSameBunch]) {
      values[kIsNoSameBunch] = event.selection_bit(o2::aod::evsel::kNoSameBunchPileup);
    }
    if (fgUsedVars[kIsGoodZvtxFT0vsPV]) {
      values[kIsGoodZvtxFT0vsPV] = event.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV);
    }
    if (fgUsedVars[kIsVertexITSTPC]) {
      values[kIsVertexITSTPC] = event.selection_bit(o2::aod::evsel::kIsVertexITSTPC);
    }
    if (fgUsedVars[kIsVertexTOFmatched]) {
      values[kIsVertexTOFmatched] = event.selection_bit(o2::aod::evsel::kIsVertexTOFmatched);
    }
    if (fgUsedVars[kIsSel8]) {
      values[kIsSel8] = event.selection_bit(o2::aod::evsel::kIsTriggerTVX) && event.selection_bit(o2::aod::evsel::kNoITSROFrameBorder) && event.selection_bit(o2::aod::evsel::kNoTimeFrameBorder);
    }
    if (fgUsedVars[kIsINT7]) {
      values[kIsINT7] = (event.alias_bit(kINT7) > 0);
    }
    if (fgUsedVars[kIsEMC7]) {
      values[kIsEMC7] = (event.alias_bit(kEMC7) > 0);
    }
    if (fgUsedVars[kIsINT7inMUON]) {
      values[kIsINT7inMUON] = (event.alias_bit(kINT7inMUON) > 0);
    }
    if (fgUsedVars[kIsMuonSingleLowPt7]) {
      values[kIsMuonSingleLowPt7] = (event.alias_bit(kMuonSingleLowPt7) > 0);
    }
    if (fgUsedVars[kIsMuonSingleHighPt7]) {
      values[kIsMuonSingleHighPt7] = (event.alias_bit(kMuonSingleHighPt7) > 0);
    }
    if (fgUsedVars[kIsMuonUnlikeLowPt7]) {
      values[kIsMuonUnlikeLowPt7] = (event.alias_bit(kMuonUnlikeLowPt7) > 0);
    }
    if (fgUsedVars[kIsMuonLikeLowPt7]) {
      values[kIsMuonLikeLowPt7] = (event.alias_bit(kMuonLikeLowPt7) > 0);
    }
    if (fgUsedVars[kIsCUP8]) {
      values[kIsCUP8] = (event.alias_bit(kCUP8) > 0);
    }
    if (fgUsedVars[kIsCUP9]) {
      values[kIsCUP9] = (event.alias_bit(kCUP9) > 0);
    }
    if (fgUsedVars[kIsMUP10]) {
      values[kIsMUP10] = (event.alias_bit(kMUP10) > 0);
    }
    if (fgUsedVars[kIsMUP11]) {
      values[kIsMUP11] = (event.alias_bit(kMUP11) > 0);
    }
    values[kVtxX] = event.posX();
    values[kVtxY] = event.posY();
    values[kVtxZ] = event.posZ();
    values[kVtxNcontrib] = event.numContrib();
    values[kVtxCovXX] = event.covXX();
    values[kVtxCovXY] = event.covXY();
    values[kVtxCovXZ] = event.covXZ();
    values[kVtxCovYY] = event.covYY();
    values[kVtxCovYZ] = event.covYZ();
    values[kVtxCovZZ] = event.covZZ();
    values[kVtxChi2] = event.chi2();
    values[kCollisionTime] = event.collisionTime();
    values[kCollisionTimeRes] = event.collisionTimeRes();
  }

  if constexpr ((fillMap & CollisionCentRun2) > 0) {
    values[kCentVZERO] = event.centRun2V0M();
  }

  if constexpr ((fillMap & CollisionCent) > 0 || (fillMap & ReducedEventExtended) > 0) {
    values[kCentFT0C] = event.centFT0C();
  }

  if constexpr ((fillMap & CollisionMult) > 0 || (fillMap & ReducedEventExtended) > 0) {
    values[kMultTPC] = event.multTPC();
    values[kMultFV0A] = event.multFV0A();
    values[kMultFV0C] = event.multFV0C();
    values[kMultFT0A] = event.multFT0A();
    values[kMultFT0C] = event.multFT0C();
    values[kMultFDDA] = event.multFDDA();
    values[kMultFDDC] = event.multFDDC();
    values[kMultZNA] = event.multZNA();
    values[kMultZNC] = event.multZNC();
    values[kMultTracklets] = event.multTracklets();
    values[kVtxNcontribReal] = event.multNTracksPV();
  }

  if constexpr ((fillMap & CollisionMultExtra) > 0 || (fillMap & ReducedEventMultExtra) > 0) {
    values[kMultNTracksHasITS] = event.multNTracksHasITS();
    values[kMultNTracksHasTPC] = event.multNTracksHasTPC();
    values[kMultNTracksHasTOF] = event.multNTracksHasTOF();
    values[kMultNTracksHasTRD] = event.multNTracksHasTRD();
    values[kMultNTracksITSOnly] = event.multNTracksITSOnly();
    values[kMultNTracksTPCOnly] = event.multNTracksTPCOnly();
    values[kMultNTracksITSTPC] = event.multNTracksITSTPC();
    values[kTrackOccupancyInTimeRange] = event.trackOccupancyInTimeRange();
    values[kMultAllTracksTPCOnly] = event.multAllTracksTPCOnly();
    values[kMultAllTracksITSTPC] = event.multAllTracksITSTPC();
    if constexpr ((fillMap & ReducedEventMultExtra) > 0) {
      values[kNTPCpileupContribA] = event.nTPCpileupContribA();
      values[kNTPCpileupContribC] = event.nTPCpileupContribC();
      values[kNTPCpileupZA] = event.nTPCpileupZA();
      values[kNTPCpileupZC] = event.nTPCpileupZC();
      values[kNTPCtracksInPast] = event.nTPCtracksInPast();
      values[kNTPCtracksInFuture] = event.nTPCtracksInFuture();
    }
  }
  // TODO: need to add EvSels and Cents tables, etc. in case of the central data model

  if constexpr ((fillMap & ReducedEvent) > 0) {
    values[kRunNo] = event.runNumber();
    values[kRunIndex] = GetRunIndex(event.runNumber());
    values[kVtxX] = event.posX();
    values[kVtxY] = event.posY();
    values[kVtxZ] = event.posZ();
    values[kVtxNcontrib] = event.numContrib();
    if (fgUsedVars[kIsNoITSROFBorder]) {
      values[kIsNoITSROFBorder] = (event.selection_bit(o2::aod::evsel::kNoITSROFrameBorder) > 0);
    }
    if (fgUsedVars[kIsNoTFBorder]) {
      values[kIsNoTFBorder] = (event.selection_bit(o2::aod::evsel::kNoTimeFrameBorder) > 0);
    }
    if (fgUsedVars[kIsNoSameBunch]) {
      values[kIsNoSameBunch] = (event.selection_bit(o2::aod::evsel::kNoSameBunchPileup) > 0);
    }
    if (fgUsedVars[kIsGoodZvtxFT0vsPV]) {
      values[kIsGoodZvtxFT0vsPV] = (event.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV) > 0);
    }
    if (fgUsedVars[kIsVertexITSTPC]) {
      values[kIsVertexITSTPC] = (event.selection_bit(o2::aod::evsel::kIsVertexITSTPC) > 0);
    }
    if (fgUsedVars[kIsVertexTOFmatched]) {
      values[kIsVertexTOFmatched] = (event.selection_bit(o2::aod::evsel::kIsVertexTOFmatched) > 0);
    }
    if (fgUsedVars[kIsSel8]) {
      values[kIsSel8] = event.selection_bit(o2::aod::evsel::kIsTriggerTVX) && event.selection_bit(o2::aod::evsel::kNoTimeFrameBorder) && event.selection_bit(o2::aod::evsel::kNoITSROFrameBorder);
    }
    if (fgUsedVars[kIsDoubleGap]) {
      values[kIsDoubleGap] = (event.tag_bit(56 + kDoubleGap) > 0);
    }
    if (fgUsedVars[kIsSingleGap] || fgUsedVars[kIsSingleGapA] || fgUsedVars[kIsSingleGapC]) {
      values[kIsSingleGapA] = (event.tag_bit(56 + kSingleGapA) > 0);
      values[kIsSingleGapC] = (event.tag_bit(56 + kSingleGapC) > 0);
      values[kIsSingleGap] = values[kIsSingleGapA] || values[kIsSingleGapC];
    }
    values[kCollisionTime] = event.collisionTime();
    values[kCollisionTimeRes] = event.collisionTimeRes();
  }

  if constexpr ((fillMap & ReducedEventExtended) > 0) {
    values[kBC] = event.globalBC();
    values[kBCOrbit] = event.globalBC() % o2::constants::lhc::LHCMaxBunches;
    values[kTimestamp] = event.timestamp();
    values[kTimeFromSOR] = (fgSOR > 0 ? (event.timestamp() - fgSOR) / 60000. : -1.0);
    values[kCentVZERO] = event.centRun2V0M();
    values[kCentFT0C] = event.centFT0C();
    if (fgUsedVars[kIsNoITSROFBorderRecomputed]) {
      uint16_t bcInITSROF = (event.globalBC() + 3564 - fgITSROFbias) % fgITSROFlength;
      values[kIsNoITSROFBorderRecomputed] = bcInITSROF > fgITSROFBorderMarginLow && bcInITSROF < fgITSROFlength - fgITSROFBorderMarginHigh ? 1.0 : 0.0;
    }
    if (fgUsedVars[kIsINT7]) {
      values[kIsINT7] = (event.alias_bit(kINT7) > 0);
    }
    if (fgUsedVars[kIsEMC7]) {
      values[kIsEMC7] = (event.alias_bit(kEMC7) > 0);
    }
    if (fgUsedVars[kIsINT7inMUON]) {
      values[kIsINT7inMUON] = (event.alias_bit(kINT7inMUON) > 0);
    }
    if (fgUsedVars[kIsMuonSingleLowPt7]) {
      values[kIsMuonSingleLowPt7] = (event.alias_bit(kMuonSingleLowPt7) > 0);
    }
    if (fgUsedVars[kIsMuonSingleHighPt7]) {
      values[kIsMuonSingleHighPt7] = (event.alias_bit(kMuonSingleHighPt7) > 0);
    }
    if (fgUsedVars[kIsMuonUnlikeLowPt7]) {
      values[kIsMuonUnlikeLowPt7] = (event.alias_bit(kMuonUnlikeLowPt7) > 0);
    }
    if (fgUsedVars[kIsMuonLikeLowPt7]) {
      values[kIsMuonLikeLowPt7] = (event.alias_bit(kMuonLikeLowPt7) > 0);
    }
    if (fgUsedVars[kIsCUP8]) {
      values[kIsCUP8] = (event.alias_bit(kCUP8) > 0);
    }
    if (fgUsedVars[kIsCUP9]) {
      values[kIsCUP9] = (event.alias_bit(kCUP9) > 0);
    }
    if (fgUsedVars[kIsMUP10]) {
      values[kIsMUP10] = (event.alias_bit(kMUP10) > 0);
    }
    if (fgUsedVars[kIsMUP11]) {
      values[kIsMUP11] = (event.alias_bit(kMUP11) > 0);
    }
  }

  if constexpr ((fillMap & ReducedEventVtxCov) > 0) {
    values[kVtxCovXX] = event.covXX();
    values[kVtxCovXY] = event.covXY();
    values[kVtxCovXZ] = event.covXZ();
    values[kVtxCovYY] = event.covYY();
    values[kVtxCovYZ] = event.covYZ();
    values[kVtxCovZZ] = event.covZZ();
    values[kVtxChi2] = event.chi2();
  }

  if constexpr ((fillMap & ReducedEventQvector) > 0) {
    values[kQ1X0A] = event.q1x0a();
    values[kQ1Y0A] = event.q1y0a();
    values[kQ1X0B] = event.q1x0b();
    values[kQ1Y0B] = event.q1y0b();
    values[kQ1X0C] = event.q1x0c();
    values[kQ1Y0C] = event.q1y0c();
    values[kQ2X0A] = event.q2x0a();
    values[kQ2Y0A] = event.q2y0a();
    values[kQ2X0B] = event.q2x0b();
    values[kQ2Y0B] = event.q2y0b();
    values[kQ2X0C] = event.q2x0c();
    values[kQ2Y0C] = event.q2y0c();
    values[kMultA] = event.multa();
    values[kMultB] = event.multb();
    values[kMultC] = event.multc();
    values[kQ3X0A] = event.q3x0a();
    values[kQ3Y0A] = event.q3y0a();
    values[kQ3X0B] = event.q3x0b();
    values[kQ3Y0B] = event.q3y0b();
    values[kQ3X0C] = event.q3x0c();
    values[kQ3Y0C] = event.q3y0c();
    values[kQ4X0A] = event.q4x0a();
    values[kQ4Y0A] = event.q4y0a();
    values[kQ4X0B] = event.q4x0b();
    values[kQ4Y0B] = event.q4y0b();
    values[kQ4X0C] = event.q4x0c();
    values[kQ4Y0C] = event.q4y0c();

    EventPlaneHelper epHelper;
    float Psi2A = epHelper.GetEventPlane(values[kQ2X0A], values[kQ2Y0A], 2);
    float Psi2B = epHelper.GetEventPlane(values[kQ2X0B], values[kQ2Y0B], 2);
    float Psi2C = epHelper.GetEventPlane(values[kQ2X0C], values[kQ2Y0C], 2);

    values[VarManager::kPsi2A] = Psi2A;
    values[VarManager::kPsi2B] = Psi2B;
    values[VarManager::kPsi2C] = Psi2C;

    if constexpr ((fillMap & ReducedEventQvectorExtra) > 0) {
      values[kQ42XA] = event.q42xa();
      values[kQ42YA] = event.q42ya();
      values[kQ23XA] = event.q23xa();
      values[kQ23YA] = event.q23ya();
      values[kS11A] = event.s11a();
      values[kS12A] = event.s12a();
      values[kS13A] = event.s13a();
      values[kS31A] = event.s31a();
    }

    if constexpr ((fillMap & ReducedEventRefFlow) > 0) {
      values[kM1111REF] = event.m1111ref();
      values[kM11REF] = event.m11ref();
      values[kCORR4REF] = event.corr4ref();
      values[kCORR2REF] = event.corr2ref();
      values[kMultA] = event.multa();
    }
  }

  if constexpr ((fillMap & CollisionQvect) > 0) {
    values[kQ1X0A] = -999;
    values[kQ1Y0A] = -999;
    values[kQ1X0B] = -999;
    values[kQ1Y0B] = -999;
    values[kQ1X0C] = -999;
    values[kQ1Y0C] = -999;
    values[kQ2X0A] = (event.qvecBPosRe() * event.nTrkBPos() + event.qvecBNegRe() * event.nTrkBNeg()) / (event.nTrkBPos() + event.nTrkBNeg());
    values[kQ2Y0A] = (event.qvecBPosIm() * event.nTrkBPos() + event.qvecBNegIm() * event.nTrkBNeg()) / (event.nTrkBPos() + event.nTrkBNeg());
    values[kQ2X0APOS] = event.qvecBPosRe();
    values[kQ2Y0APOS] = event.qvecBPosIm();
    values[kQ2X0ANEG] = event.qvecBNegRe();
    values[kQ2Y0ANEG] = event.qvecBNegIm();
    values[kQ2X0B] = event.qvecFT0ARe();
    values[kQ2Y0B] = event.qvecFT0AIm();
    values[kQ2X0C] = event.qvecFT0CRe();
    values[kQ2Y0C] = event.qvecFT0CIm();
    values[kMultA] = event.nTrkBPos() + event.nTrkBNeg();
    values[kMultAPOS] = event.nTrkBPos();
    values[kMultANEG] = event.nTrkBNeg();
    values[kMultB] = event.sumAmplFT0A();
    values[kMultC] = event.sumAmplFT0C();
    values[kQ3X0A] = -999;
    values[kQ3Y0A] = -999;
    values[kQ3X0B] = -999;
    values[kQ3Y0B] = -999;
    values[kQ3X0C] = -999;
    values[kQ3Y0C] = -999;
    values[kQ4X0A] = -999;
    values[kQ4Y0A] = -999;
    values[kQ4X0B] = -999;
    values[kQ4Y0B] = -999;
    values[kQ4X0C] = -999;
    values[kQ4Y0C] = -999;

    EventPlaneHelper epHelper;
    float Psi2A = epHelper.GetEventPlane(values[kQ2X0A], values[kQ2Y0A], 2);
    float Psi2APOS = epHelper.GetEventPlane(values[kQ2X0APOS], values[kQ2Y0APOS], 2);
    float Psi2ANEG = epHelper.GetEventPlane(values[kQ2X0ANEG], values[kQ2Y0ANEG], 2);
    float Psi2B = epHelper.GetEventPlane(values[kQ2X0B], values[kQ2Y0B], 2);
    float Psi2C = epHelper.GetEventPlane(values[kQ2X0C], values[kQ2Y0C], 2);

    values[kPsi2A] = Psi2A;
    values[kPsi2APOS] = Psi2APOS;
    values[kPsi2ANEG] = Psi2ANEG;
    values[kPsi2B] = Psi2B;
    values[kPsi2C] = Psi2C;

    values[kR2SP_AB] = (values[kQ2X0A] * values[kQ2X0B] + values[kQ2Y0A] * values[kQ2Y0B]);
    values[kR2SP_AC] = (values[kQ2X0A] * values[kQ2X0C] + values[kQ2Y0A] * values[kQ2Y0C]);
    values[kR2SP_BC] = (values[kQ2X0B] * values[kQ2X0C] + values[kQ2Y0B] * values[kQ2Y0C]);

    values[kR2EP_AB] = TMath::Cos(2 * (Psi2A - Psi2B));
    values[kR2EP_AC] = TMath::Cos(2 * (Psi2A - Psi2C));
    values[kR2EP_BC] = TMath::Cos(2 * (Psi2B - Psi2C));
  }

  if constexpr ((fillMap & CollisionMC) > 0) {
    values[kMCEventGeneratorId] = event.generatorsID();
    values[kMCEventSubGeneratorId] = event.getSubGeneratorId();
    values[kMCVtxX] = event.posX();
    values[kMCVtxY] = event.posY();
    values[kMCVtxZ] = event.posZ();
    values[kMCEventTime] = event.t();
    values[kMCEventWeight] = event.weight();
    values[kMCEventImpParam] = event.impactParameter();
  }

  if constexpr ((fillMap & ReducedEventMC) > 0) {
    values[kMCEventGeneratorId] = event.generatorsID();
    values[kMCEventGeneratorId] = -999; // to be added in reduced events
    values[kMCVtxX] = event.mcPosX();
    values[kMCVtxY] = event.mcPosY();
    values[kMCVtxZ] = event.mcPosZ();
    values[kMCEventTime] = event.t();
    values[kMCEventWeight] = event.weight();
    values[kMCEventImpParam] = event.impactParameter();
  }

  if constexpr ((fillMap & EventFilter) > 0) {
    values[kIsDoubleGap] = (event.eventFilter() & (uint64_t(1) << kDoubleGap)) > 0;
    values[kIsSingleGapA] = (event.eventFilter() & (uint64_t(1) << kSingleGapA)) > 0;
    values[kIsSingleGapC] = (event.eventFilter() & (uint64_t(1) << kSingleGapC)) > 0;
    values[kIsSingleGap] = values[kIsSingleGapA] || values[kIsSingleGapC];
  }

  if constexpr ((fillMap & ReducedZdc) > 0) {
    FillZDC(event, values);
  }

  FillEventDerived(values);
}

template <uint32_t fillMap, typename TEvent, typename TAssoc, typename TTracks>
void VarManager::FillEventTrackEstimators(TEvent const& collision, TAssoc const& assocs, TTracks const& /*tracks*/, float* values)
{
  // Compute median Z for the large dcaZ tracks in the TPC
  // This is for studies of the pileup impact on the TPC

  if (!values) {
    values = fgValues;
  }

  if constexpr ((fillMap & Track) > 0 && (fillMap & TrackDCA) > 0) {

    std::vector<float> tracksP;
    std::vector<float> tracksM;

    for (const auto& assoc : assocs) {
      auto track = assoc.template track_as<TTracks>();
      // compute the dca of this track wrt the collision
      auto trackPar = getTrackPar(track);
      std::array<float, 2> dca{1e10f, 1e10f};
      trackPar.propagateParamToDCA({collision.posX(), collision.posY(), collision.posZ()}, fgMagField, &dca);

      // if it is a displaced track longitudinally, add it to the track vector
      if (abs(dca[0]) < 3.0 && abs(dca[1]) > 4.0) {
        if (track.tgl() > 0.1) {
          tracksP.push_back(track.z());
        }
        if (track.tgl() < -0.1) {
          tracksM.push_back(track.z());
        }
      }
    } // end loop over associations

    // compute the number of pileup contributors and the median z for pileup
    if (tracksP.size() > 0) {
      std::sort(tracksP.begin(), tracksP.end());
      auto midP = tracksP.size() / 2;
      values[kNTPCpileupContribA] = tracksP.size();
      values[kNTPCpileupZA] = (tracksP.size() % 2 ? (tracksP[midP] + tracksP[midP - 1]) / 2 : tracksP[midP]);
    }

    if (tracksM.size() > 0) {
      std::sort(tracksM.begin(), tracksM.end());
      values[kNTPCpileupContribC] = tracksM.size();
      auto midM = tracksM.size() / 2;
      values[kNTPCpileupZC] = (tracksM.size() % 2 ? (tracksM[midM] + tracksM[midM - 1]) / 2 : tracksM[midM]);
    }
  }
}

template <typename T>
void VarManager::FillTwoEvents(T const& ev1, T const& ev2, float* values)
{
  if (!values) {
    values = fgValues;
  }

  values[kTwoEvPosZ1] = ev1.posZ();
  values[kTwoEvPosZ2] = ev2.posZ();
  values[kTwoEvPosR1] = std::sqrt(ev1.posX() * ev1.posX() + ev1.posY() * ev1.posY());
  values[kTwoEvPosR2] = std::sqrt(ev2.posX() * ev2.posX() + ev2.posY() * ev2.posY());
  values[kTwoEvPVcontrib1] = ev1.numContrib();
  values[kTwoEvPVcontrib2] = ev2.numContrib();
  if (ev1.numContrib() < ev2.numContrib()) {
    values[kTwoEvPosZ1] = ev2.posZ();
    values[kTwoEvPosZ2] = ev1.posZ();
    values[kTwoEvPVcontrib1] = ev2.numContrib();
    values[kTwoEvPVcontrib2] = ev1.numContrib();
    values[kTwoEvPosR1] = std::sqrt(ev2.posX() * ev2.posX() + ev2.posY() * ev2.posY());
    ;
    values[kTwoEvPosR2] = std::sqrt(ev1.posX() * ev1.posX() + ev1.posY() * ev1.posY());
  }
  values[kTwoEvDeltaZ] = ev1.posZ() - ev2.posZ();
  values[kTwoEvDeltaX] = ev1.posX() - ev2.posX();
  values[kTwoEvDeltaY] = ev1.posY() - ev2.posY();
  values[kTwoEvDeltaR] = std::sqrt(values[kTwoEvDeltaX] * values[kTwoEvDeltaX] + values[kTwoEvDeltaY] * values[kTwoEvDeltaY]);
}

template <uint32_t fillMap, typename T1, typename T2>
void VarManager::FillTwoMixEvents(T1 const& ev1, T1 const& ev2, T2 const& /*tracks1*/, T2 const& /*tracks2*/, float* values)
{
  if (!values) {
    values = fgValues;
  }
  values[kTwoEvPosZ1] = ev1.posZ();
  values[kTwoEvPosZ2] = ev2.posZ();
  values[kTwoEvPosR1] = std::sqrt(ev1.posX() * ev1.posX() + ev1.posY() * ev1.posY());
  values[kTwoEvPosR2] = std::sqrt(ev2.posX() * ev2.posX() + ev2.posY() * ev2.posY());
  values[kTwoEvPVcontrib1] = ev1.numContrib();
  values[kTwoEvPVcontrib2] = ev2.numContrib();
  /*
   uint32_t Track1Filter = 0;
   uint32_t Track2Filter = 0;
    for (auto& track1 : tracks1) { Track1Filter = uint32_t(track1.isMuonSelected());}
    for (auto& track2 : tracks2) { Track2Filter = uint32_t(track2.isMuonSelected());}
   */
  if constexpr ((fillMap & ReducedEventQvector) > 0) {
    values[kTwoR2SP1] = (ev1.q2x0b() * ev1.q2x0c() + ev1.q2y0b() * ev1.q2y0c());
    values[kTwoR2SP2] = (ev2.q2x0b() * ev2.q2x0c() + ev2.q2y0b() * ev2.q2y0c());

    if (ev1.q2y0b() * ev1.q2y0c() != 0.0) {
      values[kTwoR2EP1] = TMath::Cos(2 * (getEventPlane(2, ev1.q2x0b(), ev1.q2y0b()) - getEventPlane(2, ev1.q2x0c(), ev1.q2y0c())));
    }

    if (ev2.q2y0b() * ev2.q2y0c() != 0.0) {
      values[kTwoR2EP2] = TMath::Cos(2 * (getEventPlane(2, ev2.q2x0b(), ev2.q2y0b()) - getEventPlane(2, ev2.q2x0c(), ev2.q2y0c())));
    }
    // Tobe used for the calculation of u1q1 and u2q2
    values[kQ2X0A1] = ev1.q2x0a();
    values[kQ2X0A2] = ev2.q2x0a();
    values[kQ2Y0A1] = ev1.q2y0a();
    values[kQ2Y0A2] = ev2.q2y0a();
  }

  if (std::isnan(VarManager::fgValues[VarManager::kTwoR2SP1]) == true || std::isnan(VarManager::fgValues[VarManager::kTwoR2EP1]) == true) {
    values[kTwoR2SP1] = -999.;
    values[kTwoR2SP2] = -999.;
    values[kTwoR2EP1] = -999.;
    values[kTwoR2EP2] = -999.;
  }
}

template <uint32_t fillMap, typename T>
void VarManager::FillTrack(T const& track, float* values)
{
  if (!values) {
    values = fgValues;
  }

  if constexpr ((fillMap & TrackMFT) > 0) {
    values[kPt] = track.pt();
    values[kEta] = track.eta();
    values[kPhi] = track.phi();
    values[kMftNClusters] = track.nClusters();

    uint64_t mftClsAndFlgs = track.mftClusterSizesAndTrackFlags();
    double meanClusterSize = 0;
    for (int i = 0; i < 10; ++i) {
      double size = (mftClsAndFlgs >> (i * 6)) & 0x3fULL;
      values[kMftClusterSize + i] = (mftClsAndFlgs >> (i * 6)) & 0x3fULL;
      if (size > 0) {
        meanClusterSize += size;
      }
    }
    meanClusterSize /= track.nClusters();
    values[kMftMeanClusterSize] = meanClusterSize;
  }

  // Quantities based on the basic table (contains just kine information and filter bits)
  if constexpr ((fillMap & Track) > 0 || (fillMap & Muon) > 0 || (fillMap & ReducedTrack) > 0 || (fillMap & ReducedMuon) > 0) {
    values[kPt] = track.pt();
    values[kSignedPt] = track.pt() * track.sign();
    if (fgUsedVars[kP]) {
      values[kP] = track.p();
    }
    if (fgUsedVars[kPx]) {
      values[kPx] = track.px();
    }
    if (fgUsedVars[kPy]) {
      values[kPy] = track.py();
    }
    if (fgUsedVars[kPz]) {
      values[kPz] = track.pz();
    }
    if (fgUsedVars[kInvPt]) {
      values[kInvPt] = 1. / track.pt();
    }
    values[kEta] = track.eta();
    values[kPhi] = track.phi();
    values[kCharge] = track.sign();
    if (fgUsedVars[kPhiTPCOuter]) {
      values[kPhiTPCOuter] = track.phi() - (track.sign() > 0 ? 1.0 : -1.0) * (TMath::PiOver2() - TMath::ACos(0.22 * fgMagField / track.pt()));
      if (values[kPhiTPCOuter] > TMath::TwoPi()) {
        values[kPhiTPCOuter] -= TMath::TwoPi();
      }
      if (values[kPhiTPCOuter] < 0.0) {
        values[kPhiTPCOuter] += TMath::TwoPi();
      }
    }
    if (fgUsedVars[kTrackIsInsideTPCModule]) {
      float localSectorPhi = values[kPhiTPCOuter] - TMath::Floor(18.0 * values[kPhiTPCOuter] / TMath::TwoPi()) * (TMath::TwoPi() / 18.0);
      float edge = fgTPCInterSectorBoundary / 2.0 / 246.6; // minimal inter-sector boundary as angle
      float curvature = 3.0 * 3.33 * track.pt() / fgMagField * (1.0 - TMath::Sin(TMath::ACos(0.22 * fgMagField / track.pt())));
      if (curvature / 2.466 > edge) {
        edge = curvature / 2.466;
      }
      double min = edge;
      double max = TMath::TwoPi() / 18.0 - edge;
      values[kTrackIsInsideTPCModule] = (localSectorPhi > min && localSectorPhi < max ? 1.0 : 0.0);
    }

    if constexpr ((fillMap & ReducedTrack) > 0 && !((fillMap & Pair) > 0)) {
      // values[kIsGlobalTrack] = track.filteringFlags_bit(0);
      // values[kIsGlobalTrackSDD] = track.filteringFlags_bit(1);
      values[kIsAmbiguous] = track.isAmbiguous();

      values[kIsLegFromGamma] = track.filteringFlags_bit(VarManager::kIsConversionLeg);
      values[kIsLegFromK0S] = track.filteringFlags_bit(VarManager::kIsK0sLeg);
      values[kIsLegFromLambda] = track.filteringFlags_bit(VarManager::kIsLambdaLeg);
      values[kIsLegFromAntiLambda] = track.filteringFlags_bit(VarManager::kIsALambdaLeg);
      values[kIsLegFromOmega] = track.filteringFlags_bit(VarManager::kIsOmegaLeg);

      values[kIsProtonFromLambdaAndAntiLambda] = static_cast<bool>((values[kIsLegFromLambda] * track.sign() > 0) || (values[kIsLegFromAntiLambda] * (-track.sign()) > 0));

      for (int i = 0; i < 8; i++) {
        values[kIsDalitzLeg + i] = track.filteringFlags_bit(VarManager::kDalitzBits + i);
      }
    }
  }

  // Quantities based on the barrel tables
  if constexpr ((fillMap & TrackExtra) > 0 || (fillMap & ReducedTrackBarrel) > 0) {
    values[kPin] = track.tpcInnerParam();
    values[kSignedPin] = track.tpcInnerParam() * track.sign();
    if (fgUsedVars[kIsITSrefit]) {
      values[kIsITSrefit] = (track.flags() & o2::aod::track::ITSrefit) > 0; // NOTE: This is just for Run-2
    }
    if (fgUsedVars[kTrackTimeResIsRange]) {
      values[kTrackTimeResIsRange] = (track.flags() & o2::aod::track::TrackTimeResIsRange) > 0; // NOTE: This is NOT for Run-2
    }
    if (fgUsedVars[kIsTPCrefit]) {
      values[kIsTPCrefit] = (track.flags() & o2::aod::track::TPCrefit) > 0; // NOTE: This is just for Run-2
    }
    if (fgUsedVars[kPVContributor]) {
      values[kPVContributor] = (track.flags() & o2::aod::track::PVContributor) > 0; // NOTE: This is NOT for Run-2
    }
    if (fgUsedVars[kIsGoldenChi2]) {
      values[kIsGoldenChi2] = (track.flags() & o2::aod::track::GoldenChi2) > 0; // NOTE: This is just for Run-2
    }
    if (fgUsedVars[kOrphanTrack]) {
      values[kOrphanTrack] = (track.flags() & o2::aod::track::OrphanTrack) > 0; // NOTE: This is NOT for Run-2
    }
    if (fgUsedVars[kIsSPDfirst]) {
      values[kIsSPDfirst] = (track.itsClusterMap() & uint8_t(1)) > 0;
    }
    if (fgUsedVars[kIsSPDboth]) {
      values[kIsSPDboth] = (track.itsClusterMap() & uint8_t(3)) > 0;
    }
    if (fgUsedVars[kIsSPDany]) {
      values[kIsSPDany] = (track.itsClusterMap() & uint8_t(1)) || (track.itsClusterMap() & uint8_t(2));
    }
    if (fgUsedVars[kITSClusterMap]) {
      values[kITSClusterMap] = track.itsClusterMap();
    }

    if (fgUsedVars[kIsITSibFirst]) {
      values[kIsITSibFirst] = (track.itsClusterMap() & uint8_t(1)) > 0;
    }
    if (fgUsedVars[kIsITSibAny]) {
      values[kIsITSibAny] = (track.itsClusterMap() & (1 << uint8_t(0))) > 0 || (track.itsClusterMap() & (1 << uint8_t(1))) > 0 || (track.itsClusterMap() & (1 << uint8_t(2))) > 0;
    }
    if (fgUsedVars[kIsITSibAll]) {
      values[kIsITSibAll] = (track.itsClusterMap() & (1 << uint8_t(0))) > 0 && (track.itsClusterMap() & (1 << uint8_t(1))) > 0 && (track.itsClusterMap() & (1 << uint8_t(2))) > 0;
    }

    values[kTrackTime] = track.trackTime();
    values[kTrackTimeRes] = track.trackTimeRes();
    values[kTrackTimeResRelative] = track.trackTimeRes() / track.trackTime();
    values[kTOFExpMom] = track.tofExpMom();
    values[kITSchi2] = track.itsChi2NCl();
    values[kTPCncls] = track.tpcNClsFound();
    values[kTPCchi2] = track.tpcChi2NCl();
    values[kTrackLength] = track.length();
    values[kTPCnclsCR] = track.tpcNClsCrossedRows();
    values[kTRDPattern] = track.trdPattern();

    values[kTPCsignal] = track.tpcSignal();
    values[kTRDsignal] = track.trdSignal();

    values[kDetectorMap] = track.detectorMap();
    values[kHasITS] = track.hasITS();
    values[kHasTRD] = track.hasTRD();
    values[kHasTOF] = track.hasTOF();
    values[kHasTPC] = track.hasTPC();

    if constexpr ((fillMap & TrackExtra) > 0) {
      if (fgUsedVars[kITSncls]) {
        values[kITSncls] = track.itsNCls(); // dynamic column
      }
      if (fgUsedVars[kITSmeanClsSize]) {
        values[kITSmeanClsSize] = 0.0;
        uint32_t clsizeflag = track.itsClusterSizes();
        float mcls = 0.;
        for (unsigned int layer = 0; layer < 7; layer++) {
          mcls += (clsizeflag >> (layer * 4)) & 0xF;
        }
        if (track.itsNCls() > 0) {
          values[kITSmeanClsSize] = mcls / track.itsNCls();
        }
      }
    }
    if constexpr ((fillMap & ReducedTrackBarrel) > 0) {
      if (fgUsedVars[kITSncls]) {
        values[kITSncls] = 0.0;
        for (int i = 0; i < 7; ++i) {
          values[kITSncls] += ((track.itsClusterMap() & (1 << i)) ? 1 : 0);
        }
      }
      values[kTrackDCAxy] = track.dcaXY();
      values[kTrackDCAz] = track.dcaZ();
      if constexpr ((fillMap & ReducedTrackBarrelCov) > 0) {
        if (fgUsedVars[kTrackDCAsigXY]) {
          values[kTrackDCAsigXY] = track.dcaXY() / std::sqrt(track.cYY());
        }
        if (fgUsedVars[kTrackDCAsigZ]) {
          values[kTrackDCAsigZ] = track.dcaZ() / std::sqrt(track.cZZ());
        }
        if (fgUsedVars[kTrackDCAresXY]) {
          values[kTrackDCAresXY] = std::sqrt(track.cYY());
        }
        if (fgUsedVars[kTrackDCAresZ]) {
          values[kTrackDCAresZ] = std::sqrt(track.cZZ());
        }
      }
    }
  }

  // Quantities based on the barrel track selection table
  if constexpr ((fillMap & TrackDCA) > 0) {
    values[kTrackDCAxy] = track.dcaXY();
    values[kTrackDCAz] = track.dcaZ();
    if constexpr ((fillMap & TrackCov) > 0) {
      if (fgUsedVars[kTrackDCAsigXY]) {
        values[kTrackDCAsigXY] = track.dcaXY() / std::sqrt(track.cYY());
      }
      if (fgUsedVars[kTrackDCAsigZ]) {
        values[kTrackDCAsigZ] = track.dcaZ() / std::sqrt(track.cZZ());
      }
      if (fgUsedVars[kTrackDCAresXY]) {
        values[kTrackDCAresXY] = std::sqrt(track.cYY());
      }
      if (fgUsedVars[kTrackDCAresZ]) {
        values[kTrackDCAresZ] = std::sqrt(track.cZZ());
      }
    }
  }

  // Quantities based on the barrel track selection table
  if constexpr ((fillMap & TrackSelection) > 0) {
    values[kIsGlobalTrack] = track.isGlobalTrack();
    values[kIsGlobalTrackSDD] = track.isGlobalTrackSDD();
  }

  // Quantities based on the barrel covariance tables
  if constexpr ((fillMap & TrackCov) > 0 || (fillMap & ReducedTrackBarrelCov) > 0) {
    values[kTrackCYY] = track.cYY();
    values[kTrackCZZ] = track.cZZ();
    values[kTrackCSnpSnp] = track.cSnpSnp();
    values[kTrackCTglTgl] = track.cTglTgl();
    values[kTrackC1Pt21Pt2] = track.c1Pt21Pt2();
  }

  // Quantities based on the dalitz selections
  if constexpr ((fillMap & DalitzBits) > 0) {
    for (int i = 0; i < 8; i++) {
      values[kIsDalitzLeg + i] = static_cast<bool>(track.dalitzBits() & (uint8_t(1) << i));
    }
  }

  // Quantities based on the V0 selections
  if constexpr ((fillMap & TrackV0Bits) > 0) {
    values[kIsLegFromGamma] = static_cast<bool>(track.pidbit() & (uint8_t(1) << VarManager::kIsConversionLeg));
    values[kIsLegFromK0S] = static_cast<bool>(track.pidbit() & (uint8_t(1) << VarManager::kIsK0sLeg));
    values[kIsLegFromLambda] = static_cast<bool>(track.pidbit() & (uint8_t(1) << VarManager::kIsLambdaLeg));
    values[kIsLegFromAntiLambda] = static_cast<bool>(track.pidbit() & (uint8_t(1) << VarManager::kIsALambdaLeg));
    values[kIsLegFromOmega] = static_cast<bool>(track.pidbit() & (uint8_t(1) << VarManager::kIsOmegaLeg));
    values[kIsProtonFromLambdaAndAntiLambda] = static_cast<bool>((values[kIsLegFromLambda] * track.sign() > 0) || (values[kIsLegFromAntiLambda] * (-track.sign()) > 0));
  }

  // Quantities based on the barrel PID tables
  if constexpr ((fillMap & TrackPID) > 0 || (fillMap & ReducedTrackBarrelPID) > 0) {
    values[kTPCnSigmaEl] = track.tpcNSigmaEl();
    values[kTPCnSigmaPi] = track.tpcNSigmaPi();
    values[kTPCnSigmaKa] = track.tpcNSigmaKa();
    values[kTPCnSigmaPr] = track.tpcNSigmaPr();

    bool isTPCCalibrated = false;
    if constexpr ((fillMap & ReducedTrackBarrelPID) > 0) {
      if (track.filteringFlags_bit(kIsTPCPostcalibrated)) {
        isTPCCalibrated = true;
      }
    }
    // compute TPC postcalibrated electron nsigma based on calibration histograms from CCDB
    if (fgUsedVars[kTPCnSigmaEl_Corr] && fgRunTPCPostCalibration[0]) {
      TH3F* calibMean = reinterpret_cast<TH3F*>(fgCalibs[kTPCElectronMean]);
      TH3F* calibSigma = reinterpret_cast<TH3F*>(fgCalibs[kTPCElectronSigma]);

      int binTPCncls = calibMean->GetXaxis()->FindBin(values[kTPCncls]);
      binTPCncls = (binTPCncls == 0 ? 1 : binTPCncls);
      binTPCncls = (binTPCncls > calibMean->GetXaxis()->GetNbins() ? calibMean->GetXaxis()->GetNbins() : binTPCncls);
      int binPin = calibMean->GetYaxis()->FindBin(values[kPin]);
      binPin = (binPin == 0 ? 1 : binPin);
      binPin = (binPin > calibMean->GetYaxis()->GetNbins() ? calibMean->GetYaxis()->GetNbins() : binPin);
      int binEta = calibMean->GetZaxis()->FindBin(values[kEta]);
      binEta = (binEta == 0 ? 1 : binEta);
      binEta = (binEta > calibMean->GetZaxis()->GetNbins() ? calibMean->GetZaxis()->GetNbins() : binEta);

      double mean = calibMean->GetBinContent(binTPCncls, binPin, binEta);
      double width = calibSigma->GetBinContent(binTPCncls, binPin, binEta);
      if (!isTPCCalibrated) {
        values[kTPCnSigmaEl_Corr] = (values[kTPCnSigmaEl] - mean) / width;
      } else {
        values[kTPCnSigmaEl_Corr] = track.tpcNSigmaEl();
      }
    }
    // compute TPC postcalibrated pion nsigma if required
    if (fgUsedVars[kTPCnSigmaPi_Corr] && fgRunTPCPostCalibration[1]) {
      TH3F* calibMean = reinterpret_cast<TH3F*>(fgCalibs[kTPCPionMean]);
      TH3F* calibSigma = reinterpret_cast<TH3F*>(fgCalibs[kTPCPionSigma]);

      int binTPCncls = calibMean->GetXaxis()->FindBin(values[kTPCncls]);
      binTPCncls = (binTPCncls == 0 ? 1 : binTPCncls);
      binTPCncls = (binTPCncls > calibMean->GetXaxis()->GetNbins() ? calibMean->GetXaxis()->GetNbins() : binTPCncls);
      int binPin = calibMean->GetYaxis()->FindBin(values[kPin]);
      binPin = (binPin == 0 ? 1 : binPin);
      binPin = (binPin > calibMean->GetYaxis()->GetNbins() ? calibMean->GetYaxis()->GetNbins() : binPin);
      int binEta = calibMean->GetZaxis()->FindBin(values[kEta]);
      binEta = (binEta == 0 ? 1 : binEta);
      binEta = (binEta > calibMean->GetZaxis()->GetNbins() ? calibMean->GetZaxis()->GetNbins() : binEta);

      double mean = calibMean->GetBinContent(binTPCncls, binPin, binEta);
      double width = calibSigma->GetBinContent(binTPCncls, binPin, binEta);
      if (!isTPCCalibrated) {
        values[kTPCnSigmaPi_Corr] = (values[kTPCnSigmaPi] - mean) / width;
      } else {
        values[kTPCnSigmaPi_Corr] = track.tpcNSigmaPi();
      }
    }
    if (fgUsedVars[kTPCnSigmaKa_Corr] && fgRunTPCPostCalibration[2]) {
      TH3F* calibMean = reinterpret_cast<TH3F*>(fgCalibs[kTPCKaonMean]);
      TH3F* calibSigma = reinterpret_cast<TH3F*>(fgCalibs[kTPCKaonSigma]);

      int binTPCncls = calibMean->GetXaxis()->FindBin(values[kTPCncls]);
      binTPCncls = (binTPCncls == 0 ? 1 : binTPCncls);
      binTPCncls = (binTPCncls > calibMean->GetXaxis()->GetNbins() ? calibMean->GetXaxis()->GetNbins() : binTPCncls);
      int binPin = calibMean->GetYaxis()->FindBin(values[kPin]);
      binPin = (binPin == 0 ? 1 : binPin);
      binPin = (binPin > calibMean->GetYaxis()->GetNbins() ? calibMean->GetYaxis()->GetNbins() : binPin);
      int binEta = calibMean->GetZaxis()->FindBin(values[kEta]);
      binEta = (binEta == 0 ? 1 : binEta);
      binEta = (binEta > calibMean->GetZaxis()->GetNbins() ? calibMean->GetZaxis()->GetNbins() : binEta);

      double mean = calibMean->GetBinContent(binTPCncls, binPin, binEta);
      double width = calibSigma->GetBinContent(binTPCncls, binPin, binEta);
      if (!isTPCCalibrated) {
        values[kTPCnSigmaKa_Corr] = (values[kTPCnSigmaKa] - mean) / width;
      } else {
        values[kTPCnSigmaKa_Corr] = track.tpcNSigmaKa();
      }
    }
    // compute TPC postcalibrated proton nsigma if required
    if (fgUsedVars[kTPCnSigmaPr_Corr] && fgRunTPCPostCalibration[3]) {
      TH3F* calibMean = reinterpret_cast<TH3F*>(fgCalibs[kTPCProtonMean]);
      TH3F* calibSigma = reinterpret_cast<TH3F*>(fgCalibs[kTPCProtonSigma]);

      int binTPCncls = calibMean->GetXaxis()->FindBin(values[kTPCncls]);
      binTPCncls = (binTPCncls == 0 ? 1 : binTPCncls);
      binTPCncls = (binTPCncls > calibMean->GetXaxis()->GetNbins() ? calibMean->GetXaxis()->GetNbins() : binTPCncls);
      int binPin = calibMean->GetYaxis()->FindBin(values[kPin]);
      binPin = (binPin == 0 ? 1 : binPin);
      binPin = (binPin > calibMean->GetYaxis()->GetNbins() ? calibMean->GetYaxis()->GetNbins() : binPin);
      int binEta = calibMean->GetZaxis()->FindBin(values[kEta]);
      binEta = (binEta == 0 ? 1 : binEta);
      binEta = (binEta > calibMean->GetZaxis()->GetNbins() ? calibMean->GetZaxis()->GetNbins() : binEta);

      double mean = calibMean->GetBinContent(binTPCncls, binPin, binEta);
      double width = calibSigma->GetBinContent(binTPCncls, binPin, binEta);
      if (!isTPCCalibrated) {
        values[kTPCnSigmaPr_Corr] = (values[kTPCnSigmaPr] - mean) / width;
      } else {
        values[kTPCnSigmaPr_Corr] = track.tpcNSigmaPr();
      }
    }
    values[kTOFnSigmaEl] = track.tofNSigmaEl();
    values[kTOFnSigmaPi] = track.tofNSigmaPi();
    values[kTOFnSigmaKa] = track.tofNSigmaKa();
    values[kTOFnSigmaPr] = track.tofNSigmaPr();

    if (fgUsedVars[kTPCsignalRandomized] || fgUsedVars[kTPCnSigmaElRandomized] || fgUsedVars[kTPCnSigmaPiRandomized] || fgUsedVars[kTPCnSigmaPrRandomized]) {
      // NOTE: this is needed temporarily for the study of the impact of TPC pid degradation on the quarkonium triggers in high lumi pp
      //     This study involves a degradation from a dE/dx resolution of 5% to one of 6% (20% worsening)
      //     For this we smear the dE/dx and n-sigmas using a gaus distribution with a width of 3.3%
      //         which is approx the needed amount to get dE/dx to a resolution of 6%
      double randomX = gRandom->Gaus(0.0, 0.033);
      values[kTPCsignalRandomized] = values[kTPCsignal] * (1.0 + randomX);
      values[kTPCsignalRandomizedDelta] = values[kTPCsignal] * randomX;
      values[kTPCnSigmaElRandomized] = values[kTPCnSigmaEl] * (1.0 + randomX);
      values[kTPCnSigmaElRandomizedDelta] = values[kTPCnSigmaEl] * randomX;
      values[kTPCnSigmaPiRandomized] = values[kTPCnSigmaPi] * (1.0 + randomX);
      values[kTPCnSigmaPiRandomizedDelta] = values[kTPCnSigmaPi] * randomX;
      values[kTPCnSigmaPrRandomized] = values[kTPCnSigmaPr] * (1.0 + randomX);
      values[kTPCnSigmaPrRandomizedDelta] = values[kTPCnSigmaPr] * randomX;
    }

    if constexpr ((fillMap & ReducedTrackBarrelPID) > 0) {
      values[kTPCnSigmaMu] = track.tpcNSigmaMu();
      values[kTOFnSigmaMu] = track.tofNSigmaMu();
      values[kTPCsignal] = track.tpcSignal();
      values[kTRDsignal] = track.trdSignal();
      values[kTOFbeta] = track.beta();
    }
  }
  if constexpr ((fillMap & TrackTPCPID) > 0) {
    values[kTPCnSigmaEl] = track.tpcNSigmaEl();
    values[kTPCnSigmaPi] = track.tpcNSigmaPi();
    values[kTPCnSigmaKa] = track.tpcNSigmaKa();
    values[kTPCnSigmaPr] = track.tpcNSigmaPr();
  }
  if constexpr ((fillMap & TrackPIDExtra) > 0) {
    values[kTPCnSigmaMu] = track.tpcNSigmaMu();
    values[kTOFnSigmaMu] = track.tofNSigmaMu();
    values[kTOFbeta] = track.beta();
  }

  // Quantities based on the muon extra table
  if constexpr ((fillMap & ReducedMuonExtra) > 0 || (fillMap & Muon) > 0) {
    values[kMuonNClusters] = track.nClusters();
    values[kMuonPDca] = track.pDca();
    values[kMCHBitMap] = track.mchBitMap();
    values[kMuonRAtAbsorberEnd] = track.rAtAbsorberEnd();
    values[kMuonChi2] = track.chi2();
    values[kMuonChi2MatchMCHMID] = track.chi2MatchMCHMID();
    values[kMuonChi2MatchMCHMFT] = track.chi2MatchMCHMFT();
    values[kMuonMatchScoreMCHMFT] = track.matchScoreMCHMFT();
    values[kMuonTrackType] = track.trackType();
    values[kMuonDCAx] = track.fwdDcaX();
    values[kMuonDCAy] = track.fwdDcaY();
    values[kMuonTime] = track.trackTime();
    values[kMuonTimeRes] = track.trackTimeRes();
  }
  // Quantities based on the muon covariance table
  if constexpr ((fillMap & ReducedMuonCov) > 0 || (fillMap & MuonCov) > 0) {
    values[kX] = track.x();
    values[kY] = track.y();
    values[kZ] = track.z();
    values[kTgl] = track.tgl();
    values[kMuonCXX] = track.cXX();
    values[kMuonCXY] = track.cXY();
    values[kMuonCYY] = track.cYY();
    values[kMuonCPhiX] = track.cPhiX();
    values[kMuonCPhiY] = track.cPhiY();
    values[kMuonCPhiPhi] = track.cPhiPhi();
    values[kMuonCTglX] = track.cTglX();
    values[kMuonCTglY] = track.cTglY();
    values[kMuonCTglPhi] = track.cTglPhi();
    values[kMuonCTglTgl] = track.cTglTgl();
    values[kMuonC1Pt2X] = track.c1PtX();
    values[kMuonC1Pt2Y] = track.c1PtY();
    values[kMuonC1Pt2Phi] = track.c1PtPhi();
    values[kMuonC1Pt2Tgl] = track.c1PtTgl();
    values[kMuonC1Pt21Pt2] = track.c1Pt21Pt2();
  }

  // Quantities based on the pair table(s)
  if constexpr ((fillMap & Pair) > 0) {
    values[kMass] = track.mass();
    ROOT::Math::PtEtaPhiMVector vpair(track.pt(), track.eta(), track.phi(), track.mass());
    values[kRap] = vpair.Rapidity();
  }

  // Derived quantities which can be computed based on already filled variables
  FillTrackDerived(values);
}

template <uint32_t fillMap, typename T, typename C>
void VarManager::FillTrackCollision(T const& track, C const& collision, float* values)
{

  if (!values) {
    values = fgValues;
  }
  if constexpr ((fillMap & ReducedTrackBarrel) > 0 || (fillMap & TrackDCA) > 0) {
    auto trackPar = getTrackPar(track);
    std::array<float, 2> dca{1e10f, 1e10f};
    trackPar.propagateParamToDCA({collision.posX(), collision.posY(), collision.posZ()}, fgMagField, &dca);

    values[kTrackDCAxy] = dca[0];
    values[kTrackDCAz] = dca[1];

    if constexpr ((fillMap & ReducedTrackBarrelCov) > 0 || (fillMap & TrackCov) > 0) {
      if (fgUsedVars[kTrackDCAsigXY]) {
        values[kTrackDCAsigXY] = dca[0] / std::sqrt(track.cYY());
      }
      if (fgUsedVars[kTrackDCAsigZ]) {
        values[kTrackDCAsigZ] = dca[1] / std::sqrt(track.cZZ());
      }
    }
  }
  if constexpr ((fillMap & MuonCov) > 0 || (fillMap & ReducedMuonCov) > 0) {

    o2::dataformats::GlobalFwdTrack propmuonAtDCA = PropagateMuon(track, collision, kToDCA);

    float dcaX = (propmuonAtDCA.getX() - collision.posX());
    float dcaY = (propmuonAtDCA.getY() - collision.posY());
    float dcaXY = std::sqrt(dcaX * dcaX + dcaY * dcaY);
    values[kMuonPDca] = track.p() * dcaXY;
  }
}

template <uint32_t fillMap, typename T, typename C, typename M, typename P>
void VarManager::FillTrackCollisionMatCorr(T const& track, C const& collision, M const& materialCorr, P const& propagator, float* values)
{
  if (!values) {
    values = fgValues;
  }
  if constexpr ((fillMap & ReducedTrackBarrel) > 0 || (fillMap & TrackDCA) > 0) {
    auto trackPar = getTrackPar(track);
    std::array<float, 2> dca{1e10f, 1e10f};
    std::array<float, 3> pVec = {track.px(), track.py(), track.pz()};
    // trackPar.propagateParamToDCA({collision.posX(), collision.posY(), collision.posZ()}, fgMagField, &dca);
    propagator->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackPar, 2.f, materialCorr, &dca);
    getPxPyPz(trackPar, pVec);

    values[kTrackDCAxy] = dca[0];
    values[kTrackDCAz] = dca[1];

    if constexpr ((fillMap & ReducedTrackBarrelCov) > 0 || (fillMap & TrackCov) > 0) {
      if (fgUsedVars[kTrackDCAsigXY]) {
        values[kTrackDCAsigXY] = dca[0] / std::sqrt(track.cYY());
      }
      if (fgUsedVars[kTrackDCAsigZ]) {
        values[kTrackDCAsigZ] = dca[1] / std::sqrt(track.cZZ());
      }
    }
  }
}

template <uint32_t fillMap, typename T>
void VarManager::FillPhoton(T const& track, float* values)
{
  if (!values) {
    values = fgValues;
  }

  // Quantities based on the basic table (contains just kine information and filter bits)
  if constexpr ((fillMap & Track) > 0 || (fillMap & ReducedTrack) > 0) {
    values[kPt] = track.pt();
    if (fgUsedVars[kP]) {
      values[kP] = track.p();
    }
    if (fgUsedVars[kPx]) {
      values[kPx] = track.px();
    }
    if (fgUsedVars[kPy]) {
      values[kPy] = track.py();
    }
    if (fgUsedVars[kPz]) {
      values[kPz] = track.pz();
    }
    if (fgUsedVars[kInvPt]) {
      values[kInvPt] = 1. / track.pt();
    }
    values[kEta] = track.eta();
    values[kPhi] = track.phi();
    values[kRap] = track.eta(); // photon does not know rapidity .y()
    values[kMassDau] = track.mGamma();
  }
}

template <typename U, typename T>
void VarManager::FillTrackMC(const U& mcStack, T const& track, float* values)
{
  if (!values) {
    values = fgValues;
  }

  // Quantities based on the mc particle table
  values[kMCPdgCode] = track.pdgCode();
  values[kMCParticleWeight] = track.weight();
  values[kMCPx] = track.px();
  values[kMCPy] = track.py();
  values[kMCPz] = track.pz();
  values[kMCE] = track.e();
  values[kMCVx] = track.vx();
  values[kMCVy] = track.vy();
  values[kMCVz] = track.vz();
  values[kMCPt] = track.pt();
  values[kMCPhi] = track.phi();
  values[kMCEta] = track.eta();
  values[kMCY] = -track.y();
  values[kMCParticleGeneratorId] = track.producedByGenerator();
  if (fgUsedVars[kMCMotherPdgCode]) {
    if (track.has_mothers()) {
      auto motherId = track.mothersIds()[0];
      auto mother = mcStack.rawIteratorAt(motherId);
      values[kMCMotherPdgCode] = mother.pdgCode();
    }
  }

  FillTrackDerived(values);
}

template <uint32_t fillMap, typename T1, typename T2, typename C>
void VarManager::FillPairPropagateMuon(T1 const& muon1, T2 const& muon2, const C& collision, float* values)
{
  if (!values) {
    values = fgValues;
  }
  o2::dataformats::GlobalFwdTrack propmuon1 = PropagateMuon(muon1, collision);
  o2::dataformats::GlobalFwdTrack propmuon2 = PropagateMuon(muon2, collision);

  float m = o2::constants::physics::MassMuon;

  ROOT::Math::PtEtaPhiMVector v1(propmuon1.getPt(), propmuon1.getEta(), propmuon1.getPhi(), m);
  ROOT::Math::PtEtaPhiMVector v2(propmuon2.getPt(), propmuon2.getEta(), propmuon2.getPhi(), m);
  ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
  values[kMass] = v12.M();
  values[kPt] = v12.Pt();
  values[kEta] = v12.Eta();
  values[kPhi] = v12.Phi();
  values[kRap] = -v12.Rapidity();

  double Ptot1 = TMath::Sqrt(v1.Px() * v1.Px() + v1.Py() * v1.Py() + v1.Pz() * v1.Pz());
  double Ptot2 = TMath::Sqrt(v2.Px() * v2.Px() + v2.Py() * v2.Py() + v2.Pz() * v2.Pz());
  values[kDeltaPtotTracks] = Ptot1 - Ptot2;
}

template <int pairType, uint32_t fillMap, typename T1, typename T2>
void VarManager::FillPair(T1 const& t1, T2 const& t2, float* values)
{
  if (!values) {
    values = fgValues;
  }

  float m1 = o2::constants::physics::MassElectron;
  float m2 = o2::constants::physics::MassElectron;
  if constexpr (pairType == kDecayToMuMu) {
    m1 = o2::constants::physics::MassMuon;
    m2 = o2::constants::physics::MassMuon;
  }

  if constexpr (pairType == kDecayToPiPi) {
    m1 = o2::constants::physics::MassPionCharged;
    m2 = o2::constants::physics::MassPionCharged;
  }

  if constexpr (pairType == kDecayToKPi) {
    m1 = o2::constants::physics::MassKaonCharged;
    m2 = o2::constants::physics::MassPionCharged;
  }

  if constexpr (pairType == kElectronMuon) {
    m2 = o2::constants::physics::MassMuon;
  }

  ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), m1);
  ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), m2);
  ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
  values[kMass] = v12.M();
  values[kPt] = v12.Pt();
  values[kEta] = v12.Eta();
  // values[kPhi] = v12.Phi();
  values[kPhi] = v12.Phi() > 0 ? v12.Phi() : v12.Phi() + 2. * M_PI;
  values[kRap] = -v12.Rapidity();
  double Ptot1 = TMath::Sqrt(v1.Px() * v1.Px() + v1.Py() * v1.Py() + v1.Pz() * v1.Pz());
  double Ptot2 = TMath::Sqrt(v2.Px() * v2.Px() + v2.Py() * v2.Py() + v2.Pz() * v2.Pz());
  values[kDeltaPtotTracks] = Ptot1 - Ptot2;

  if (fgUsedVars[kDeltaPhiPair2]) {
    double phipair2 = v1.Phi() - v2.Phi();
    if (phipair2 > 3 * TMath::Pi() / 2) {
      values[kDeltaPhiPair2] = phipair2 - 2 * TMath::Pi();
    } else if (phipair2 < -TMath::Pi() / 2) {
      values[kDeltaPhiPair2] = phipair2 + 2 * TMath::Pi();
    } else {
      values[kDeltaPhiPair2] = phipair2;
    }
  }

  if (fgUsedVars[kPsiPair]) {
    values[kDeltaPhiPair] = (t1.sign() * fgMagField > 0.) ? (v1.Phi() - v2.Phi()) : (v2.Phi() - v1.Phi());
    double xipair = TMath::ACos((v1.Px() * v2.Px() + v1.Py() * v2.Py() + v1.Pz() * v2.Pz()) / v1.P() / v2.P());
    values[kPsiPair] = (t1.sign() * fgMagField > 0.) ? TMath::ASin((v1.Theta() - v2.Theta()) / xipair) : TMath::ASin((v2.Theta() - v1.Theta()) / xipair);
  }

  if (fgUsedVars[kOpeningAngle]) {
    double scalar = v1.Px() * v2.Px() + v1.Py() * v2.Py() + v1.Pz() * v2.Pz();
    double Ptot12 = Ptot1 * Ptot2;
    if (Ptot12 <= 0) {
      values[kOpeningAngle] = 0.;
    } else {
      double arg = scalar / Ptot12;
      if (arg > 1.)
        arg = 1.;
      if (arg < -1)
        arg = -1;
      values[kOpeningAngle] = TMath::ACos(arg);
    }
  }

  // TO DO: get the correct values from CCDB
  double BeamMomentum = TMath::Sqrt(fgCenterOfMassEnergy * fgCenterOfMassEnergy / 4 - fgMassofCollidingParticle * fgMassofCollidingParticle); // GeV
  ROOT::Math::PxPyPzEVector Beam1(0., 0., -BeamMomentum, fgCenterOfMassEnergy / 2);
  ROOT::Math::PxPyPzEVector Beam2(0., 0., BeamMomentum, fgCenterOfMassEnergy / 2);

  // Boost to center of mass frame
  ROOT::Math::Boost boostv12{v12.BoostToCM()};
  ROOT::Math::XYZVectorF v1_CM{(boostv12(v1).Vect()).Unit()};
  ROOT::Math::XYZVectorF v2_CM{(boostv12(v2).Vect()).Unit()};
  ROOT::Math::XYZVectorF Beam1_CM{(boostv12(Beam1).Vect()).Unit()};
  ROOT::Math::XYZVectorF Beam2_CM{(boostv12(Beam2).Vect()).Unit()};

  // Helicity frame
  ROOT::Math::XYZVectorF zaxis_HE{(v12.Vect()).Unit()};
  ROOT::Math::XYZVectorF yaxis_HE{(Beam1_CM.Cross(Beam2_CM)).Unit()};
  ROOT::Math::XYZVectorF xaxis_HE{(yaxis_HE.Cross(zaxis_HE)).Unit()};

  // Collins-Soper frame
  ROOT::Math::XYZVectorF zaxis_CS{((Beam1_CM.Unit() - Beam2_CM.Unit()).Unit())};
  ROOT::Math::XYZVectorF yaxis_CS{(Beam1_CM.Cross(Beam2_CM)).Unit()};
  ROOT::Math::XYZVectorF xaxis_CS{(yaxis_CS.Cross(zaxis_CS)).Unit()};

  if (fgUsedVars[kCosThetaHE]) {
    values[kCosThetaHE] = (t1.sign() > 0 ? zaxis_HE.Dot(v1_CM) : zaxis_HE.Dot(v2_CM));
  }

  if (fgUsedVars[kPhiHE]) {
    values[kPhiHE] = (t1.sign() > 0 ? TMath::ATan2(yaxis_HE.Dot(v1_CM), xaxis_HE.Dot(v1_CM)) : TMath::ATan2(yaxis_HE.Dot(v2_CM), xaxis_HE.Dot(v2_CM)));
  }

  if (fgUsedVars[kCosThetaCS]) {
    values[kCosThetaCS] = (t1.sign() > 0 ? zaxis_CS.Dot(v1_CM) : zaxis_CS.Dot(v2_CM));
  }

  if (fgUsedVars[kPhiCS]) {
    values[kPhiCS] = (t1.sign() > 0 ? TMath::ATan2(yaxis_CS.Dot(v1_CM), xaxis_CS.Dot(v1_CM)) : TMath::ATan2(yaxis_CS.Dot(v2_CM), xaxis_CS.Dot(v2_CM)));
  }

  if constexpr ((pairType == kDecayToEE) && ((fillMap & TrackCov) > 0 || (fillMap & ReducedTrackBarrelCov) > 0)) {

    if (fgUsedVars[kQuadDCAabsXY] || fgUsedVars[kQuadDCAsigXY] || fgUsedVars[kQuadDCAabsZ] || fgUsedVars[kQuadDCAsigZ] || fgUsedVars[kQuadDCAsigXYZ] || fgUsedVars[kSignQuadDCAsigXY]) {
      // Quantities based on the barrel tables
      double dca1XY = t1.dcaXY();
      double dca2XY = t2.dcaXY();
      double dca1Z = t1.dcaZ();
      double dca2Z = t2.dcaZ();
      double dca1sigXY = dca1XY / std::sqrt(t1.cYY());
      double dca2sigXY = dca2XY / std::sqrt(t2.cYY());
      double dca1sigZ = dca1Z / std::sqrt(t1.cZZ());
      double dca2sigZ = dca2Z / std::sqrt(t2.cZZ());

      values[kQuadDCAabsXY] = std::sqrt((dca1XY * dca1XY + dca2XY * dca2XY) / 2);
      values[kQuadDCAsigXY] = std::sqrt((dca1sigXY * dca1sigXY + dca2sigXY * dca2sigXY) / 2);
      values[kQuadDCAabsZ] = std::sqrt((dca1Z * dca1Z + dca2Z * dca2Z) / 2);
      values[kQuadDCAsigZ] = std::sqrt((dca1sigZ * dca1sigZ + dca2sigZ * dca2sigZ) / 2);
      values[kSignQuadDCAsigXY] = t1.sign() * t2.sign() * TMath::Sign(1., dca1sigXY) * TMath::Sign(1., dca2sigXY) * std::sqrt((dca1sigXY * dca1sigXY + dca2sigXY * dca2sigXY) / 2);

      double det1 = t1.cYY() * t1.cZZ() - t1.cZY() * t1.cZY();
      double det2 = t2.cYY() * t2.cZZ() - t2.cZY() * t2.cZY();
      if ((det1 < 0) || (det2 < 0)) {
        values[kQuadDCAsigXYZ] = -999;
      } else {
        double chi2t1 = (dca1XY * dca1XY * t1.cZZ() + dca1Z * dca1Z * t1.cYY() - 2. * dca1XY * dca1Z * t1.cZY()) / det1;
        double chi2t2 = (dca2XY * dca2XY * t2.cZZ() + dca2Z * dca2Z * t2.cYY() - 2. * dca2XY * dca2Z * t2.cZY()) / det2;

        double dca1sigXYZ = std::sqrt(std::abs(chi2t1) / 2.);
        double dca2sigXYZ = std::sqrt(std::abs(chi2t2) / 2.);

        values[kQuadDCAsigXYZ] = std::sqrt((dca1sigXYZ * dca1sigXYZ + dca2sigXYZ * dca2sigXYZ) / 2);
      }
    }
  }
  if constexpr ((pairType == kDecayToMuMu) && ((fillMap & Muon) > 0 || (fillMap & ReducedMuon) > 0)) {
    if (fgUsedVars[kQuadDCAabsXY]) {
      double dca1X = t1.fwdDcaX();
      double dca1Y = t1.fwdDcaY();
      double dca1XY = std::sqrt(dca1X * dca1X + dca1Y * dca1Y);
      double dca2X = t2.fwdDcaX();
      double dca2Y = t2.fwdDcaY();
      double dca2XY = std::sqrt(dca2X * dca2X + dca2Y * dca2Y);
      values[kQuadDCAabsXY] = std::sqrt((dca1XY * dca1XY + dca2XY * dca2XY) / 2.);
    }
  }
  if (fgUsedVars[kPairPhiv]) {
    values[kPairPhiv] = calculatePhiV<pairType>(t1, t2);
  }
}

template <int pairType, uint32_t fillMap, typename C, typename T1, typename T2>
void VarManager::FillPairCollision(const C& collision, T1 const& t1, T2 const& t2, float* values)
{
  if (!values) {
    values = fgValues;
  }

  if constexpr ((pairType == kDecayToEE) && ((fillMap & TrackCov) > 0 || (fillMap & ReducedTrackBarrelCov) > 0)) {

    if (fgUsedVars[kQuadDCAabsXY] || fgUsedVars[kQuadDCAsigXY] || fgUsedVars[kQuadDCAabsZ] || fgUsedVars[kQuadDCAsigZ] || fgUsedVars[kQuadDCAsigXYZ] || fgUsedVars[kSignQuadDCAsigXY]) {

      auto trackPart1 = getTrackPar(t1);
      std::array<float, 2> dca1{1e10f, 1e10f};
      trackPart1.propagateParamToDCA({collision.posX(), collision.posY(), collision.posZ()}, fgMagField, &dca1);

      auto trackPart2 = getTrackPar(t2);
      std::array<float, 2> dca2{1e10f, 1e10f};
      trackPart2.propagateParamToDCA({collision.posX(), collision.posY(), collision.posZ()}, fgMagField, &dca2);

      // Recalculated quantities
      double dca1XY = dca1[0];
      double dca2XY = dca2[0];
      double dca1Z = dca1[1];
      double dca2Z = dca2[1];
      double dca1sigXY = dca1XY / std::sqrt(t1.cYY());
      double dca2sigXY = dca2XY / std::sqrt(t2.cYY());
      double dca1sigZ = dca1Z / std::sqrt(t1.cZZ());
      double dca2sigZ = dca2Z / std::sqrt(t2.cZZ());

      values[kQuadDCAabsXY] = std::sqrt((dca1XY * dca1XY + dca2XY * dca2XY) / 2);
      values[kQuadDCAsigXY] = std::sqrt((dca1sigXY * dca1sigXY + dca2sigXY * dca2sigXY) / 2);
      values[kQuadDCAabsZ] = std::sqrt((dca1Z * dca1Z + dca2Z * dca2Z) / 2);
      values[kQuadDCAsigZ] = std::sqrt((dca1sigZ * dca1sigZ + dca2sigZ * dca2sigZ) / 2);
      values[kSignQuadDCAsigXY] = t1.sign() * t2.sign() * TMath::Sign(1., dca1sigXY) * TMath::Sign(1., dca2sigXY) * std::sqrt((dca1sigXY * dca1sigXY + dca2sigXY * dca2sigXY) / 2);

      double det1 = t1.cYY() * t1.cZZ() - t1.cZY() * t1.cZY();
      double det2 = t2.cYY() * t2.cZZ() - t2.cZY() * t2.cZY();
      if ((det1 < 0) || (det2 < 0)) {
        values[kQuadDCAsigXYZ] = -999;
      } else {
        double chi2t1 = (dca1XY * dca1XY * t1.cZZ() + dca1Z * dca1Z * t1.cYY() - 2. * dca1XY * dca1Z * t1.cZY()) / det1;
        double chi2t2 = (dca2XY * dca2XY * t2.cZZ() + dca2Z * dca2Z * t2.cYY() - 2. * dca2XY * dca2Z * t2.cZY()) / det2;

        double dca1sigXYZ = std::sqrt(std::abs(chi2t1) / 2.);
        double dca2sigXYZ = std::sqrt(std::abs(chi2t2) / 2.);

        values[kQuadDCAsigXYZ] = std::sqrt((dca1sigXYZ * dca1sigXYZ + dca2sigXYZ * dca2sigXYZ) / 2);
      }
    }
  }
}

template <int pairType, uint32_t fillMap, typename C, typename T1, typename T2, typename M, typename P>
void VarManager::FillPairCollisionMatCorr(C const& collision, T1 const& t1, T2 const& t2, M const& materialCorr, P const& propagator, float* values)
{
  if (!values) {
    values = fgValues;
  }

  if constexpr ((pairType == kDecayToEE) && ((fillMap & TrackCov) > 0 || (fillMap & ReducedTrackBarrelCov) > 0)) {

    if (fgUsedVars[kQuadDCAabsXY] || fgUsedVars[kQuadDCAsigXY] || fgUsedVars[kQuadDCAabsZ] || fgUsedVars[kQuadDCAsigZ] || fgUsedVars[kQuadDCAsigXYZ] || fgUsedVars[kSignQuadDCAsigXY]) {

      auto trackPart1 = getTrackPar(t1);
      std::array<float, 2> dca1{1e10f, 1e10f};
      std::array<float, 3> pVect1 = {t1.px(), t1.py(), t1.pz()};
      // trackPar.propagateParamToDCA({collision.posX(), collision.posY(), collision.posZ()}, fgMagField, &dca);
      propagator->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackPart1, 2.f, materialCorr, &dca1);
      getPxPyPz(trackPart1, pVect1);

      auto trackPart2 = getTrackPar(t2);
      std::array<float, 2> dca2{1e10f, 1e10f};
      std::array<float, 3> pVect2 = {t2.px(), t2.py(), t2.pz()};
      // trackPar.propagateParamToDCA({collision.posX(), collision.posY(), collision.posZ()}, fgMagField, &dca);
      propagator->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackPart2, 2.f, materialCorr, &dca2);
      getPxPyPz(trackPart2, pVect2);

      // Recalculated quantities
      double dca1XY = dca1[0];
      double dca2XY = dca2[0];
      double dca1Z = dca1[1];
      double dca2Z = dca2[1];
      double dca1sigXY = dca1XY / std::sqrt(t1.cYY());
      double dca2sigXY = dca2XY / std::sqrt(t2.cYY());
      double dca1sigZ = dca1Z / std::sqrt(t1.cZZ());
      double dca2sigZ = dca2Z / std::sqrt(t2.cZZ());

      values[kQuadDCAabsXY] = std::sqrt((dca1XY * dca1XY + dca2XY * dca2XY) / 2);
      values[kQuadDCAsigXY] = std::sqrt((dca1sigXY * dca1sigXY + dca2sigXY * dca2sigXY) / 2);
      values[kQuadDCAabsZ] = std::sqrt((dca1Z * dca1Z + dca2Z * dca2Z) / 2);
      values[kQuadDCAsigZ] = std::sqrt((dca1sigZ * dca1sigZ + dca2sigZ * dca2sigZ) / 2);
      values[kSignQuadDCAsigXY] = t1.sign() * t2.sign() * TMath::Sign(1., dca1sigXY) * TMath::Sign(1., dca2sigXY) * std::sqrt((dca1sigXY * dca1sigXY + dca2sigXY * dca2sigXY) / 2);

      double det1 = t1.cYY() * t1.cZZ() - t1.cZY() * t1.cZY();
      double det2 = t2.cYY() * t2.cZZ() - t2.cZY() * t2.cZY();
      if ((det1 < 0) || (det2 < 0)) {
        values[kQuadDCAsigXYZ] = -999;
      } else {
        double chi2t1 = (dca1XY * dca1XY * t1.cZZ() + dca1Z * dca1Z * t1.cYY() - 2. * dca1XY * dca1Z * t1.cZY()) / det1;
        double chi2t2 = (dca2XY * dca2XY * t2.cZZ() + dca2Z * dca2Z * t2.cYY() - 2. * dca2XY * dca2Z * t2.cZY()) / det2;

        double dca1sigXYZ = std::sqrt(std::abs(chi2t1) / 2.);
        double dca2sigXYZ = std::sqrt(std::abs(chi2t2) / 2.);

        values[kQuadDCAsigXYZ] = std::sqrt((dca1sigXYZ * dca1sigXYZ + dca2sigXYZ * dca2sigXYZ) / 2);
      }
    }
  }
}

template <typename T1, typename T2, typename T3>
void VarManager::FillTriple(T1 const& t1, T2 const& t2, T3 const& t3, float* values, PairCandidateType pairType)
{

  if (!values) {
    values = fgValues;
  }
  if (pairType == kTripleCandidateToEEPhoton) {
    float m1 = o2::constants::physics::MassElectron;
    float m3 = o2::constants::physics::MassPhoton;
    float m4 = o2::constants::physics::MassJPsi;

    ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), m1);
    ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), m1);
    ROOT::Math::PtEtaPhiMVector v3(t3.pt(), t3.eta(), t3.phi(), m3);
    ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
    ROOT::Math::PtEtaPhiMVector v123 = v12 + v3;
    values[kPairMass] = v123.M();
    values[kPairPt] = v123.Pt();
    values[kPairEta] = v123.Eta();
    values[kPairPhi] = v123.Phi();
    values[kPairMassDau] = v12.M();
    values[kMassDau] = m3;
    values[kPairPtDau] = v12.Pt();
    values[kPt] = t3.pt();
    values[kEta] = t3.eta();
    values[kEta1] = t1.eta();
    values[kEta2] = t2.eta();
    values[kDeltaEta] = v12.Eta();
    values[VarManager::kDeltaMass] = v123.M() - v12.M();
    values[VarManager::kDeltaMass_jpsi] = v123.M() - v12.M() + m4;
    values[kRap] = v123.Rapidity();
    values[kPt1] = t1.pt();
    values[kPt2] = t2.pt();
  }

  if (pairType == kTripleCandidateToKPiPi) {
    float m1 = o2::constants::physics::MassKaonCharged;
    float m2 = o2::constants::physics::MassPionCharged;

    ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), m1);
    ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), m2);
    ROOT::Math::PtEtaPhiMVector v3(t3.pt(), t3.eta(), t3.phi(), m2);
    ROOT::Math::PtEtaPhiMVector v123 = v1 + v2 + v3;
    values[kMass] = v123.M();
    values[kPt] = v123.Pt();
    values[kEta] = v123.Eta();
    values[kPhi] = v123.Phi();
    values[kRap] = -v123.Rapidity();
  }

  if (pairType == kTripleCandidateToPKPi) {
    float m1 = o2::constants::physics::MassProton;
    float m2 = o2::constants::physics::MassKaonCharged;
    float m3 = o2::constants::physics::MassPionCharged;

    ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), m1);
    ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), m2);
    ROOT::Math::PtEtaPhiMVector v3(t3.pt(), t3.eta(), t3.phi(), m3);
    ROOT::Math::PtEtaPhiMVector v123 = v1 + v2 + v3;
    values[kMass] = v123.M();
    values[kPt] = v123.Pt();
    values[kEta] = v123.Eta();
    values[kPhi] = v123.Phi();
    values[kRap] = -v123.Rapidity();
  }
}

template <int pairType, typename T1, typename T2>
void VarManager::FillPairME(T1 const& t1, T2 const& t2, float* values)
{
  //
  // Lightweight fill function called from the innermost event mixing loop
  //
  if (!values) {
    values = fgValues;
  }

  float m1 = o2::constants::physics::MassElectron;
  float m2 = o2::constants::physics::MassElectron;
  if constexpr (pairType == kDecayToMuMu) {
    m1 = o2::constants::physics::MassMuon;
    m2 = o2::constants::physics::MassMuon;
  }

  if constexpr (pairType == kDecayToPiPi) {
    m1 = o2::constants::physics::MassPionCharged;
    m2 = o2::constants::physics::MassPionCharged;
  }

  if constexpr (pairType == kElectronMuon) {
    m2 = o2::constants::physics::MassMuon;
  }

  ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), m1);
  ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), m2);
  ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
  values[kMass] = v12.M();
  values[kPt] = v12.Pt();
  values[kEta] = v12.Eta();
  // values[kPhi] = v12.Phi();
  values[kPhi] = v12.Phi() > 0 ? v12.Phi() : v12.Phi() + 2. * M_PI;
  values[kRap] = -v12.Rapidity();

  // TODO: provide different computations for vn
  // Compute the scalar product UQ for two muon from different event using Q-vector from A, for second and third harmonic
  values[kU2Q2Ev1] = values[kQ2X0A1] * std::cos(2 * v1.Phi()) + values[kQ2Y0A1] * std::sin(2 * v1.Phi());
  values[kU2Q2Ev2] = values[kQ2X0A2] * std::cos(2 * v2.Phi()) + values[kQ2Y0A2] * std::sin(2 * v2.Phi());
  values[kCos2DeltaPhiMu1] = std::cos(2 * (v1.Phi() - v12.Phi()));
  values[kCos2DeltaPhiMu2] = std::cos(2 * (v2.Phi() - v12.Phi()));

  if constexpr (pairType == kDecayToMuMu) {
    if (fgUsedVars[kQuadDCAabsXY]) {
      double dca1X = t1.fwdDcaX();
      double dca1Y = t1.fwdDcaY();
      double dca1XY = std::sqrt(dca1X * dca1X + dca1Y * dca1Y);
      double dca2X = t2.fwdDcaX();
      double dca2Y = t2.fwdDcaY();
      double dca2XY = std::sqrt(dca2X * dca2X + dca2Y * dca2Y);
      values[kQuadDCAabsXY] = std::sqrt((dca1XY * dca1XY + dca2XY * dca2XY) / 2.);
    }
  }
  if (fgUsedVars[kPairPhiv]) {
    values[kPairPhiv] = calculatePhiV<pairType>(t1, t2);
  }
}

template <typename T1, typename T2>
void VarManager::FillPairMC(T1 const& t1, T2 const& t2, float* values, PairCandidateType pairType)
{
  if (!values) {
    values = fgValues;
  }

  float m1 = o2::constants::physics::MassElectron;
  float m2 = o2::constants::physics::MassElectron;
  if (pairType == kDecayToMuMu) {
    m1 = o2::constants::physics::MassMuon;
    m2 = o2::constants::physics::MassMuon;
  }

  if (pairType == kDecayToPiPi) {
    m1 = o2::constants::physics::MassPionCharged;
    m2 = o2::constants::physics::MassPionCharged;
  }

  if (pairType == kElectronMuon) {
    m2 = o2::constants::physics::MassMuon;
  }

  // TODO : implement resolution smearing.
  ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), m1);
  ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), m2);
  ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
  values[kMass] = v12.M();
  values[kPt] = v12.Pt();
  values[kEta] = v12.Eta();
  values[kPhi] = v12.Phi();
  values[kRap] = -v12.Rapidity();
}

template <typename T1, typename T2, typename T3>
void VarManager::FillTripleMC(T1 const& t1, T2 const& t2, T3 const& t3, float* values, PairCandidateType pairType)
{
  if (!values) {
    values = fgValues;
  }

  if (pairType == kTripleCandidateToEEPhoton) {
    float m1 = o2::constants::physics::MassElectron;
    float m2 = o2::constants::physics::MassElectron;
    float m3 = o2::constants::physics::MassPhoton;
    float m4 = o2::constants::physics::MassJPsi;
    ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), m1);
    ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), m2);
    ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
    ROOT::Math::PtEtaPhiMVector v3(t3.pt(), t3.eta(), t3.phi(), m3);
    ROOT::Math::PtEtaPhiMVector v123 = v12 + v3;
    values[kPairMass] = v123.M();
    values[kPairPt] = v123.Pt();
    values[kPairEta] = v123.Eta();
    values[kPhi] = v123.Phi();
    values[kPairMassDau] = v12.M();
    values[kMassDau] = m3;
    values[kPairPtDau] = v12.Pt();
    values[kPt] = t3.pt();
    values[kEta] = t3.eta();
    values[kEta1] = t1.eta();
    values[kEta2] = t2.eta();
    values[kDeltaEta] = v12.Eta();
    values[VarManager::kDeltaMass] = v123.M() - v12.M();
    values[VarManager::kDeltaMass_jpsi] = v123.M() - v12.M() + m4;
    values[kRap] = -v123.Rapidity();
    values[kPt1] = t1.pt();
    values[kPt2] = t2.pt();
  }
}

template <int pairType, uint32_t collFillMap, uint32_t fillMap, typename C, typename T>
void VarManager::FillPairVertexing(C const& collision, T const& t1, T const& t2, bool propToSV, float* values)
{
  // check at compile time that the event and cov matrix have the cov matrix
  constexpr bool eventHasVtxCov = ((collFillMap & Collision) > 0 || (collFillMap & ReducedEventVtxCov) > 0);
  constexpr bool trackHasCov = ((fillMap & TrackCov) > 0 || (fillMap & ReducedTrackBarrelCov) > 0);
  constexpr bool muonHasCov = ((fillMap & MuonCov) > 0 || (fillMap & ReducedMuonCov) > 0);

  if (!values) {
    values = fgValues;
  }
  float m1 = o2::constants::physics::MassElectron;
  float m2 = o2::constants::physics::MassElectron;
  if constexpr (pairType == kDecayToKPi) {
    m1 = o2::constants::physics::MassKaonCharged;
    m2 = o2::constants::physics::MassPionCharged;
  }
  if constexpr (pairType == kDecayToMuMu && muonHasCov) {
    m1 = o2::constants::physics::MassMuon;
    m2 = o2::constants::physics::MassMuon;
  }
  ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), m1);
  ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), m2);
  ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;

  values[kUsedKF] = fgUsedKF;
  if (!fgUsedKF) {
    int procCode = 0;

    // TODO: use trackUtilities functions to initialize the various matrices to avoid code duplication
    // auto pars1 = getTrackParCov(t1);
    // auto pars2 = getTrackParCov(t2);
    // We need to hide the cov data members from the cases when no cov table is provided
    if constexpr ((pairType == kDecayToEE || pairType == kDecayToKPi) && trackHasCov) {
      std::array<float, 5> t1pars = {t1.y(), t1.z(), t1.snp(), t1.tgl(), t1.signed1Pt()};
      std::array<float, 15> t1covs = {t1.cYY(), t1.cZY(), t1.cZZ(), t1.cSnpY(), t1.cSnpZ(),
                                      t1.cSnpSnp(), t1.cTglY(), t1.cTglZ(), t1.cTglSnp(), t1.cTglTgl(),
                                      t1.c1PtY(), t1.c1PtZ(), t1.c1PtSnp(), t1.c1PtTgl(), t1.c1Pt21Pt2()};
      o2::track::TrackParCov pars1{t1.x(), t1.alpha(), t1pars, t1covs};
      std::array<float, 5> t2pars = {t2.y(), t2.z(), t2.snp(), t2.tgl(), t2.signed1Pt()};
      std::array<float, 15> t2covs = {t2.cYY(), t2.cZY(), t2.cZZ(), t2.cSnpY(), t2.cSnpZ(),
                                      t2.cSnpSnp(), t2.cTglY(), t2.cTglZ(), t2.cTglSnp(), t2.cTglTgl(),
                                      t2.c1PtY(), t2.c1PtZ(), t2.c1PtSnp(), t2.c1PtTgl(), t2.c1Pt21Pt2()};
      o2::track::TrackParCov pars2{t2.x(), t2.alpha(), t2pars, t2covs};
      procCode = fgFitterTwoProngBarrel.process(pars1, pars2);
    } else if constexpr ((pairType == kDecayToMuMu) && muonHasCov) {
      // Initialize track parameters for forward
      double chi21 = t1.chi2();
      double chi22 = t2.chi2();
      SMatrix5 t1pars(t1.x(), t1.y(), t1.phi(), t1.tgl(), t1.signed1Pt());
      std::vector<double> v1{t1.cXX(), t1.cXY(), t1.cYY(), t1.cPhiX(), t1.cPhiY(),
                             t1.cPhiPhi(), t1.cTglX(), t1.cTglY(), t1.cTglPhi(), t1.cTglTgl(),
                             t1.c1PtX(), t1.c1PtY(), t1.c1PtPhi(), t1.c1PtTgl(), t1.c1Pt21Pt2()};
      SMatrix55 t1covs(v1.begin(), v1.end());
      o2::track::TrackParCovFwd pars1{t1.z(), t1pars, t1covs, chi21};
      SMatrix5 t2pars(t2.x(), t2.y(), t2.phi(), t2.tgl(), t2.signed1Pt());
      std::vector<double> v2{t2.cXX(), t2.cXY(), t2.cYY(), t2.cPhiX(), t2.cPhiY(),
                             t2.cPhiPhi(), t2.cTglX(), t2.cTglY(), t2.cTglPhi(), t2.cTglTgl(),
                             t2.c1PtX(), t2.c1PtY(), t2.c1PtPhi(), t2.c1PtTgl(), t2.c1Pt21Pt2()};
      SMatrix55 t2covs(v2.begin(), v2.end());
      o2::track::TrackParCovFwd pars2{t2.z(), t2pars, t2covs, chi22};
      procCode = fgFitterTwoProngFwd.process(pars1, pars2);
    } else {
      return;
    }

    values[kVertexingProcCode] = procCode;
    if (procCode == 0) {
      // TODO: set the other variables to appropriate values and return
      values[kVertexingChi2PCA] = -999.;
      values[kVertexingLxy] = -999.;
      values[kVertexingLxyz] = -999.;
      values[kVertexingLz] = -999.;
      values[kVertexingLxyErr] = -999.;
      values[kVertexingLxyzErr] = -999.;
      values[kVertexingLzErr] = -999.;

      values[kVertexingTauxy] = -999.;
      values[kVertexingTauz] = -999.;
      values[kVertexingTauxyErr] = -999.;
      values[kVertexingTauzErr] = -999.;
      values[kVertexingPz] = -999.;
      values[kVertexingSV] = -999.;
      return;
    }

    Vec3D secondaryVertex;

    if constexpr (eventHasVtxCov) {

      std::array<float, 6> covMatrixPCA;
      // get track impact parameters
      // This modifies track momenta!
      o2::math_utils::Point3D<float> vtxXYZ(collision.posX(), collision.posY(), collision.posZ());
      std::array<float, 6> vtxCov{collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ()};
      o2::dataformats::VertexBase primaryVertex = {std::move(vtxXYZ), std::move(vtxCov)};
      // auto primaryVertex = getPrimaryVertex(collision);
      auto covMatrixPV = primaryVertex.getCov();

      if constexpr ((pairType == kDecayToEE || pairType == kDecayToKPi) && trackHasCov) {
        secondaryVertex = fgFitterTwoProngBarrel.getPCACandidate();
        covMatrixPCA = fgFitterTwoProngBarrel.calcPCACovMatrixFlat();
        auto chi2PCA = fgFitterTwoProngBarrel.getChi2AtPCACandidate();
        auto trackParVar0 = fgFitterTwoProngBarrel.getTrack(0);
        auto trackParVar1 = fgFitterTwoProngBarrel.getTrack(1);
        values[kVertexingChi2PCA] = chi2PCA;
        v1 = {trackParVar0.getPt(), trackParVar0.getEta(), trackParVar0.getPhi(), m1};
        v2 = {trackParVar1.getPt(), trackParVar1.getEta(), trackParVar1.getPhi(), m2};
        v12 = v1 + v2;

      } else if constexpr (pairType == kDecayToMuMu && muonHasCov) {
        // Get pca candidate from forward DCA fitter
        secondaryVertex = fgFitterTwoProngFwd.getPCACandidate();
        covMatrixPCA = fgFitterTwoProngFwd.calcPCACovMatrixFlat();
        auto chi2PCA = fgFitterTwoProngFwd.getChi2AtPCACandidate();
        auto trackParVar0 = fgFitterTwoProngFwd.getTrack(0);
        auto trackParVar1 = fgFitterTwoProngFwd.getTrack(1);
        values[kVertexingChi2PCA] = chi2PCA;
        v1 = {trackParVar0.getPt(), trackParVar0.getEta(), trackParVar0.getPhi(), m1};
        v2 = {trackParVar1.getPt(), trackParVar1.getEta(), trackParVar1.getPhi(), m2};
        v12 = v1 + v2;

        values[kPt1] = trackParVar0.getPt();
        values[kEta1] = trackParVar0.getEta();
        values[kPhi1] = trackParVar0.getPhi();

        values[kPt2] = trackParVar1.getPt();
        values[kEta2] = trackParVar1.getEta();
        values[kPhi2] = trackParVar1.getPhi();
      }
      double phi = std::atan2(secondaryVertex[1] - collision.posY(), secondaryVertex[0] - collision.posX());
      double theta = std::atan2(secondaryVertex[2] - collision.posZ(),
                                std::sqrt((secondaryVertex[0] - collision.posX()) * (secondaryVertex[0] - collision.posX()) +
                                          (secondaryVertex[1] - collision.posY()) * (secondaryVertex[1] - collision.posY())));

      values[kVertexingLxyzErr] = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, theta) + getRotatedCovMatrixXX(covMatrixPCA, phi, theta));
      values[kVertexingLxyErr] = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, 0.) + getRotatedCovMatrixXX(covMatrixPCA, phi, 0.));
      values[kVertexingLzErr] = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, 0, theta) + getRotatedCovMatrixXX(covMatrixPCA, 0, theta));

      values[kVertexingLxy] = (collision.posX() - secondaryVertex[0]) * (collision.posX() - secondaryVertex[0]) +
                              (collision.posY() - secondaryVertex[1]) * (collision.posY() - secondaryVertex[1]);
      values[kVertexingLz] = (collision.posZ() - secondaryVertex[2]) * (collision.posZ() - secondaryVertex[2]);
      values[kVertexingLxyz] = values[kVertexingLxy] + values[kVertexingLz];
      values[kVertexingLxy] = std::sqrt(values[kVertexingLxy]);
      values[kVertexingLz] = std::sqrt(values[kVertexingLz]);
      values[kVertexingLxyz] = std::sqrt(values[kVertexingLxyz]);

      values[kVertexingTauz] = (collision.posZ() - secondaryVertex[2]) * v12.M() / (TMath::Abs(v12.Pz()) * o2::constants::physics::LightSpeedCm2NS);
      values[kVertexingTauxy] = values[kVertexingLxy] * v12.M() / (v12.Pt() * o2::constants::physics::LightSpeedCm2NS);

      values[kVertexingPz] = TMath::Abs(v12.Pz());
      values[kVertexingSV] = secondaryVertex[2];

      values[kVertexingTauzErr] = values[kVertexingLzErr] * v12.M() / (TMath::Abs(v12.Pz()) * o2::constants::physics::LightSpeedCm2NS);
      values[kVertexingTauxyErr] = values[kVertexingLxyErr] * v12.M() / (v12.Pt() * o2::constants::physics::LightSpeedCm2NS);

      values[kCosPointingAngle] = ((collision.posX() - secondaryVertex[0]) * v12.Px() +
                                   (collision.posY() - secondaryVertex[1]) * v12.Py() +
                                   (collision.posZ() - secondaryVertex[2]) * v12.Pz()) /
                                  (v12.P() * values[VarManager::kVertexingLxyz]);
      // Decay length defined as in Run 2
      values[kVertexingLzProjected] = ((secondaryVertex[2] - collision.posZ()) * v12.Pz()) / TMath::Sqrt(v12.Pz() * v12.Pz());
      values[kVertexingLxyProjected] = ((secondaryVertex[0] - collision.posX()) * v12.Px()) + ((secondaryVertex[1] - collision.posY()) * v12.Py());
      values[kVertexingLxyProjected] = values[kVertexingLxyProjected] / TMath::Sqrt((v12.Px() * v12.Px()) + (v12.Py() * v12.Py()));
      values[kVertexingLxyzProjected] = ((secondaryVertex[0] - collision.posX()) * v12.Px()) + ((secondaryVertex[1] - collision.posY()) * v12.Py()) + ((secondaryVertex[2] - collision.posZ()) * v12.Pz());
      values[kVertexingLxyzProjected] = values[kVertexingLxyzProjected] / TMath::Sqrt((v12.Px() * v12.Px()) + (v12.Py() * v12.Py()) + (v12.Pz() * v12.Pz()));
      values[kVertexingTauxyProjected] = values[kVertexingLxyProjected] * v12.M() / (v12.Pt());
      values[kVertexingTauxyProjectedNs] = values[kVertexingTauxyProjected] / o2::constants::physics::LightSpeedCm2NS;
      values[kVertexingTauzProjected] = values[kVertexingLzProjected] * v12.M() / TMath::Abs(v12.Pz());
      values[kVertexingTauxyzProjected] = values[kVertexingLxyzProjected] * v12.M() / (v12.P());
    }
  } else {
    KFParticle trk0KF;
    KFParticle trk1KF;
    KFParticle KFGeoTwoProng;
    if constexpr ((pairType == kDecayToEE) && trackHasCov) {
      KFPTrack kfpTrack0 = createKFPTrackFromTrack(t1);
      trk0KF = KFParticle(kfpTrack0, -11 * t1.sign());
      KFPTrack kfpTrack1 = createKFPTrackFromTrack(t2);
      trk1KF = KFParticle(kfpTrack1, -11 * t2.sign());

      KFGeoTwoProng.SetConstructMethod(2);
      KFGeoTwoProng.AddDaughter(trk0KF);
      KFGeoTwoProng.AddDaughter(trk1KF);

    } else if constexpr ((pairType == kDecayToMuMu) && muonHasCov) {
      KFPTrack kfpTrack0 = createKFPFwdTrackFromFwdTrack(t1);
      trk0KF = KFParticle(kfpTrack0, -13 * t1.sign());
      KFPTrack kfpTrack1 = createKFPFwdTrackFromFwdTrack(t2);
      trk1KF = KFParticle(kfpTrack1, -13 * t2.sign());

      KFGeoTwoProng.SetConstructMethod(2);
      KFGeoTwoProng.AddDaughter(trk0KF);
      KFGeoTwoProng.AddDaughter(trk1KF);

    } else if constexpr ((pairType == kDecayToKPi) && trackHasCov) {
      KFPTrack kfpTrack0 = createKFPTrackFromTrack(t1);
      trk0KF = KFParticle(kfpTrack0, 321 * t1.sign());
      KFPTrack kfpTrack1 = createKFPTrackFromTrack(t2);
      trk1KF = KFParticle(kfpTrack1, 211 * t2.sign());

      KFGeoTwoProng.SetConstructMethod(2);
      KFGeoTwoProng.AddDaughter(trk0KF);
      KFGeoTwoProng.AddDaughter(trk1KF);
    }
    if (fgUsedVars[kKFMass]) {
      float mass = 0., massErr = 0.;
      if (!KFGeoTwoProng.GetMass(mass, massErr))
        values[kKFMass] = mass;
      else
        values[kKFMass] = -999.;
    }

    if constexpr (eventHasVtxCov) {
      KFPVertex kfpVertex = createKFPVertexFromCollision(collision);
      values[kKFNContributorsPV] = kfpVertex.GetNContributors();
      KFParticle KFPV(kfpVertex);
      double dxPair2PV = KFGeoTwoProng.GetX() - KFPV.GetX();
      double dyPair2PV = KFGeoTwoProng.GetY() - KFPV.GetY();
      double dzPair2PV = KFGeoTwoProng.GetZ() - KFPV.GetZ();
      if (fgUsedVars[kVertexingLxy] || fgUsedVars[kVertexingLz] || fgUsedVars[kVertexingLxyz] || fgUsedVars[kVertexingLxyErr] || fgUsedVars[kVertexingLzErr] || fgUsedVars[kVertexingTauxy] || fgUsedVars[kVertexingLxyOverErr] || fgUsedVars[kVertexingLzOverErr] || fgUsedVars[kVertexingLxyzOverErr] || fgUsedVars[kCosPointingAngle]) {
        values[kVertexingLxy] = std::sqrt(dxPair2PV * dxPair2PV + dyPair2PV * dyPair2PV);
        values[kVertexingLz] = std::sqrt(dzPair2PV * dzPair2PV);
        values[kVertexingLxyz] = std::sqrt(dxPair2PV * dxPair2PV + dyPair2PV * dyPair2PV + dzPair2PV * dzPair2PV);
        values[kVertexingLxyErr] = (KFPV.GetCovariance(0) + KFGeoTwoProng.GetCovariance(0)) * dxPair2PV * dxPair2PV + (KFPV.GetCovariance(2) + KFGeoTwoProng.GetCovariance(2)) * dyPair2PV * dyPair2PV + 2 * ((KFPV.GetCovariance(1) + KFGeoTwoProng.GetCovariance(1)) * dxPair2PV * dyPair2PV);
        values[kVertexingLzErr] = (KFPV.GetCovariance(5) + KFGeoTwoProng.GetCovariance(5)) * dzPair2PV * dzPair2PV;
        values[kVertexingLxyzErr] = (KFPV.GetCovariance(0) + KFGeoTwoProng.GetCovariance(0)) * dxPair2PV * dxPair2PV + (KFPV.GetCovariance(2) + KFGeoTwoProng.GetCovariance(2)) * dyPair2PV * dyPair2PV + (KFPV.GetCovariance(5) + KFGeoTwoProng.GetCovariance(5)) * dzPair2PV * dzPair2PV + 2 * ((KFPV.GetCovariance(1) + KFGeoTwoProng.GetCovariance(1)) * dxPair2PV * dyPair2PV + (KFPV.GetCovariance(3) + KFGeoTwoProng.GetCovariance(3)) * dxPair2PV * dzPair2PV + (KFPV.GetCovariance(4) + KFGeoTwoProng.GetCovariance(4)) * dyPair2PV * dzPair2PV);
        if (fabs(values[kVertexingLxy]) < 1.e-8f)
          values[kVertexingLxy] = 1.e-8f;
        values[kVertexingLxyErr] = values[kVertexingLxyErr] < 0. ? 1.e8f : std::sqrt(values[kVertexingLxyErr]) / values[kVertexingLxy];
        if (fabs(values[kVertexingLz]) < 1.e-8f)
          values[kVertexingLz] = 1.e-8f;
        values[kVertexingLzErr] = values[kVertexingLzErr] < 0. ? 1.e8f : std::sqrt(values[kVertexingLzErr]) / values[kVertexingLz];
        if (fabs(values[kVertexingLxyz]) < 1.e-8f)
          values[kVertexingLxyz] = 1.e-8f;
        values[kVertexingLxyzErr] = values[kVertexingLxyzErr] < 0. ? 1.e8f : std::sqrt(values[kVertexingLxyzErr]) / values[kVertexingLxyz];
        values[kVertexingTauxy] = KFGeoTwoProng.GetPseudoProperDecayTime(KFPV, KFGeoTwoProng.GetMass()) / (o2::constants::physics::LightSpeedCm2NS);
        values[kVertexingTauz] = -1 * dzPair2PV * KFGeoTwoProng.GetMass() / (TMath::Abs(KFGeoTwoProng.GetPz()) * o2::constants::physics::LightSpeedCm2NS);
        values[kVertexingPz] = TMath::Abs(KFGeoTwoProng.GetPz());
        values[kVertexingSV] = KFGeoTwoProng.GetZ();
        values[kVertexingTauxyErr] = values[kVertexingLxyErr] * KFGeoTwoProng.GetMass() / (KFGeoTwoProng.GetPt() * o2::constants::physics::LightSpeedCm2NS);
        values[kVertexingTauzErr] = values[kVertexingLzErr] * KFGeoTwoProng.GetMass() / (TMath::Abs(KFGeoTwoProng.GetPz()) * o2::constants::physics::LightSpeedCm2NS);
        values[kCosPointingAngle] = (std::sqrt(dxPair2PV * dxPair2PV) * v12.Px() +
                                     std::sqrt(dyPair2PV * dyPair2PV) * v12.Py() +
                                     std::sqrt(dzPair2PV * dzPair2PV) * v12.Pz()) /
                                    (v12.P() * values[VarManager::kVertexingLxyz]);
      }
      // As defined in Run 2 (projected onto momentum)
      if (fgUsedVars[kVertexingLxyProjected] || fgUsedVars[kVertexingLxyzProjected] || fgUsedVars[kVertexingLzProjected]) {
        values[kVertexingLzProjected] = (dzPair2PV * KFGeoTwoProng.GetPz()) / TMath::Sqrt(KFGeoTwoProng.GetPz() * KFGeoTwoProng.GetPz());
        values[kVertexingLxyProjected] = (dxPair2PV * KFGeoTwoProng.GetPx()) + (dyPair2PV * KFGeoTwoProng.GetPy());
        values[kVertexingLxyProjected] = values[kVertexingLxyProjected] / TMath::Sqrt((KFGeoTwoProng.GetPx() * KFGeoTwoProng.GetPx()) + (KFGeoTwoProng.GetPy() * KFGeoTwoProng.GetPy()));
        values[kVertexingLxyzProjected] = (dxPair2PV * KFGeoTwoProng.GetPx()) + (dyPair2PV * KFGeoTwoProng.GetPy()) + (dzPair2PV * KFGeoTwoProng.GetPz());
        values[kVertexingLxyzProjected] = values[kVertexingLxyzProjected] / TMath::Sqrt((KFGeoTwoProng.GetPx() * KFGeoTwoProng.GetPx()) + (KFGeoTwoProng.GetPy() * KFGeoTwoProng.GetPy()) + (KFGeoTwoProng.GetPz() * KFGeoTwoProng.GetPz()));
        values[kVertexingTauxyProjected] = values[kVertexingLxyProjected] * KFGeoTwoProng.GetMass() / (KFGeoTwoProng.GetPt());
        values[kVertexingTauxyProjectedNs] = values[kVertexingTauxyProjected] / o2::constants::physics::LightSpeedCm2NS;
        values[kVertexingTauzProjected] = values[kVertexingLzProjected] * KFGeoTwoProng.GetMass() / TMath::Abs(KFGeoTwoProng.GetPz());
      }

      if (fgUsedVars[kVertexingLxyOverErr] || fgUsedVars[kVertexingLzOverErr] || fgUsedVars[kVertexingLxyzOverErr]) {
        values[kVertexingLxyOverErr] = values[kVertexingLxy] / values[kVertexingLxyErr];
        values[kVertexingLzOverErr] = values[kVertexingLz] / values[kVertexingLzErr];
        values[kVertexingLxyzOverErr] = values[kVertexingLxyz] / values[kVertexingLxyzErr];
      }

      if (fgUsedVars[kKFChi2OverNDFGeo])
        values[kKFChi2OverNDFGeo] = KFGeoTwoProng.GetChi2() / KFGeoTwoProng.GetNDF();
      if (fgUsedVars[kKFCosPA])
        values[kKFCosPA] = calculateCosPA(KFGeoTwoProng, KFPV);

      // in principle, they should be in FillTrack
      if (fgUsedVars[kKFTrack0DCAxyz] || fgUsedVars[kKFTrack1DCAxyz]) {
        values[kKFTrack0DCAxyz] = trk0KF.GetDistanceFromVertex(KFPV);
        values[kKFTrack1DCAxyz] = trk1KF.GetDistanceFromVertex(KFPV);
      }
      if (fgUsedVars[kKFTrack0DCAxy] || fgUsedVars[kKFTrack1DCAxy]) {
        values[kKFTrack0DCAxy] = trk0KF.GetDistanceFromVertexXY(KFPV);
        values[kKFTrack1DCAxy] = trk1KF.GetDistanceFromVertexXY(KFPV);
      }
      if (fgUsedVars[kKFDCAxyzBetweenProngs])
        values[kKFDCAxyzBetweenProngs] = trk0KF.GetDistanceFromParticle(trk1KF);
      if (fgUsedVars[kKFDCAxyBetweenProngs])
        values[kKFDCAxyBetweenProngs] = trk0KF.GetDistanceFromParticleXY(trk1KF);

      if (fgUsedVars[kKFTracksDCAxyzMax]) {
        values[kKFTracksDCAxyzMax] = values[kKFTrack0DCAxyz] > values[kKFTrack1DCAxyz] ? values[kKFTrack0DCAxyz] : values[kKFTrack1DCAxyz];
      }
      if (fgUsedVars[kKFTracksDCAxyMax]) {
        values[kKFTracksDCAxyMax] = TMath::Abs(values[kKFTrack0DCAxy]) > TMath::Abs(values[kKFTrack1DCAxy]) ? values[kKFTrack0DCAxy] : values[kKFTrack1DCAxy];
      }
      if (fgUsedVars[kKFTrack0DeviationFromPV] || fgUsedVars[kKFTrack1DeviationFromPV]) {
        values[kKFTrack0DeviationFromPV] = trk0KF.GetDeviationFromVertex(KFPV);
        values[kKFTrack1DeviationFromPV] = trk1KF.GetDeviationFromVertex(KFPV);
      }
      if (fgUsedVars[kKFTrack0DeviationxyFromPV] || fgUsedVars[kKFTrack1DeviationxyFromPV]) {
        values[kKFTrack0DeviationxyFromPV] = trk0KF.GetDeviationFromVertexXY(KFPV);
        values[kKFTrack1DeviationxyFromPV] = trk1KF.GetDeviationFromVertexXY(KFPV);
      }
      if (fgUsedVars[kKFJpsiDCAxyz]) {
        values[kKFJpsiDCAxyz] = KFGeoTwoProng.GetDistanceFromVertex(KFPV);
      }
      if (fgUsedVars[kKFJpsiDCAxy]) {
        values[kKFJpsiDCAxy] = KFGeoTwoProng.GetDistanceFromVertexXY(KFPV);
      }
      if (fgUsedVars[kKFPairDeviationFromPV] || fgUsedVars[kKFPairDeviationxyFromPV]) {
        values[kKFPairDeviationFromPV] = KFGeoTwoProng.GetDeviationFromVertex(KFPV);
        values[kKFPairDeviationxyFromPV] = KFGeoTwoProng.GetDeviationFromVertexXY(KFPV);
      }
      if (fgUsedVars[kKFChi2OverNDFGeoTop] || fgUsedVars[kKFMassGeoTop]) {
        KFParticle KFGeoTopTwoProngBarrel = KFGeoTwoProng;
        KFGeoTopTwoProngBarrel.SetProductionVertex(KFPV);
        values[kKFChi2OverNDFGeoTop] = KFGeoTopTwoProngBarrel.GetChi2() / KFGeoTopTwoProngBarrel.GetNDF();
        float mass = 0., massErr = 0.;
        if (!KFGeoTopTwoProngBarrel.GetMass(mass, massErr))
          values[kKFMassGeoTop] = mass;
        else
          values[kKFMassGeoTop] = -999.;
      }
      if (propToSV) {
        if constexpr ((pairType == kDecayToMuMu) && muonHasCov) {
          double chi21 = t1.chi2();
          double chi22 = t2.chi2();
          SMatrix5 t1pars(t1.x(), t1.y(), t1.phi(), t1.tgl(), t1.signed1Pt());
          std::vector<double> c1{t1.cXX(), t1.cXY(), t1.cYY(), t1.cPhiX(), t1.cPhiY(),
                                 t1.cPhiPhi(), t1.cTglX(), t1.cTglY(), t1.cTglPhi(), t1.cTglTgl(),
                                 t1.c1PtX(), t1.c1PtY(), t1.c1PtPhi(), t1.c1PtTgl(), t1.c1Pt21Pt2()};
          SMatrix55 t1covs(c1.begin(), c1.end());
          o2::track::TrackParCovFwd pars1{t1.z(), t1pars, t1covs, chi21};
          SMatrix5 t2pars(t2.x(), t2.y(), t2.phi(), t2.tgl(), t2.signed1Pt());
          std::vector<double> c2{t2.cXX(), t2.cXY(), t2.cYY(), t2.cPhiX(), t2.cPhiY(),
                                 t2.cPhiPhi(), t2.cTglX(), t2.cTglY(), t2.cTglPhi(), t2.cTglTgl(),
                                 t2.c1PtX(), t2.c1PtY(), t2.c1PtPhi(), t2.c1PtTgl(), t2.c1Pt21Pt2()};
          SMatrix55 t2covs(c2.begin(), c2.end());
          o2::track::TrackParCovFwd pars2{t2.z(), t2pars, t2covs, chi22};

          auto geoMan1 = o2::base::GeometryManager::meanMaterialBudget(t1.x(), t1.y(), t1.z(), KFGeoTwoProng.GetX(), KFGeoTwoProng.GetY(), KFGeoTwoProng.GetZ());
          auto geoMan2 = o2::base::GeometryManager::meanMaterialBudget(t2.x(), t2.y(), t2.z(), KFGeoTwoProng.GetX(), KFGeoTwoProng.GetY(), KFGeoTwoProng.GetZ());
          auto x2x01 = static_cast<float>(geoMan1.meanX2X0);
          auto x2x02 = static_cast<float>(geoMan2.meanX2X0);
          float B[3];
          float xyz[3] = {0, 0, 0};
          KFGeoTwoProng.GetFieldValue(xyz, B);
          // TODO: find better soluton to handle cases where KF outputs negative variances
          /*float covXX = 0.1;
          float covYY = 0.1;
          if (KFGeoTwoProng.GetCovariance(0, 0) > 0) {
            covXX = KFGeoTwoProng.GetCovariance(0, 0);
          }
          if (KFGeoTwoProng.GetCovariance(1, 1) > 0) {
            covYY = KFGeoTwoProng.GetCovariance(0, 0);
          }*/
          pars1.propagateToVtxhelixWithMCS(KFGeoTwoProng.GetZ(), {KFGeoTwoProng.GetX(), KFGeoTwoProng.GetY()}, {KFGeoTwoProng.GetCovariance(0, 0), KFGeoTwoProng.GetCovariance(1, 1)}, B[2], x2x01);
          pars2.propagateToVtxhelixWithMCS(KFGeoTwoProng.GetZ(), {KFGeoTwoProng.GetX(), KFGeoTwoProng.GetY()}, {KFGeoTwoProng.GetCovariance(0, 0), KFGeoTwoProng.GetCovariance(1, 1)}, B[2], x2x02);
          v1 = {pars1.getPt(), pars1.getEta(), pars1.getPhi(), m1};
          v2 = {pars2.getPt(), pars2.getEta(), pars2.getPhi(), m2};
          v12 = v1 + v2;
          values[kMass] = v12.M();
          values[kPt] = v12.Pt();
          values[kEta] = v12.Eta();
          values[kPhi] = v12.Phi();
          values[kRap] = -v12.Rapidity();
          values[kVertexingTauxy] = KFGeoTwoProng.GetPseudoProperDecayTime(KFPV, v12.M()) / (o2::constants::physics::LightSpeedCm2NS);
          values[kVertexingTauz] = -1 * dzPair2PV * v12.M() / (TMath::Abs(v12.Pz()) * o2::constants::physics::LightSpeedCm2NS);
          values[kVertexingTauxyErr] = values[kVertexingLxyErr] * v12.M() / (v12.Pt() * o2::constants::physics::LightSpeedCm2NS);
          values[kVertexingTauzErr] = values[kVertexingLzErr] * v12.M() / (TMath::Abs(v12.Pz()) * o2::constants::physics::LightSpeedCm2NS);
          values[kVertexingPz] = TMath::Abs(v12.Pz());
          values[kVertexingSV] = KFGeoTwoProng.GetZ();

          values[kPt1] = pars1.getPt();
          values[kEta1] = pars1.getEta();
          values[kPhi1] = pars1.getPhi();

          values[kPt2] = pars2.getPt();
          values[kEta2] = pars2.getEta();
          values[kPhi2] = pars2.getPhi();
        }
      }
    }
  }
  if (propToSV) {
    values[kMass] = v12.M();
    values[kPt] = v12.Pt();
    values[kEta] = v12.Eta();
    // values[kPhi] = v12.Phi();
    values[kPhi] = v12.Phi() > 0 ? v12.Phi() : v12.Phi() + 2. * M_PI;
  } else {
    values[kPt1] = t1.pt();
    values[kEta1] = t1.eta();
    values[kPhi1] = t1.phi();

    values[kPt2] = t2.pt();
    values[kEta2] = t2.eta();
    values[kPhi2] = t2.phi();
  }
}

template <uint32_t collFillMap, uint32_t fillMap, typename C, typename T>
void VarManager::FillTripletVertexing(C const& collision, T const& t1, T const& t2, T const& t3, VarManager::PairCandidateType tripletType, float* values)
{
  // TODO: Vertexing error variables
  constexpr bool eventHasVtxCov = ((collFillMap & Collision) > 0 || (collFillMap & ReducedEventVtxCov) > 0);
  bool trackHasCov = ((fillMap & ReducedTrackBarrelCov) > 0);

  if (!values) {
    values = fgValues;
  }

  float m1, m2, m3;

  if (tripletType == kTripleCandidateToKPiPi) {
    m1 = o2::constants::physics::MassKaonCharged;
    m2 = o2::constants::physics::MassPionCharged;
    m3 = o2::constants::physics::MassPionCharged;
  }
  if (tripletType == kTripleCandidateToPKPi) {
    m1 = o2::constants::physics::MassProton;
    m2 = o2::constants::physics::MassKaonCharged;
    m3 = o2::constants::physics::MassPionCharged;
  }
  ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), m1);
  ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), m2);
  ROOT::Math::PtEtaPhiMVector v3(t3.pt(), t3.eta(), t3.phi(), m3);
  ROOT::Math::PtEtaPhiMVector v123 = v1 + v2 + v3;

  values[kUsedKF] = fgUsedKF;
  if (!fgUsedKF) {
    int procCode = 0;

    if (trackHasCov) {
      std::array<float, 5> t1pars = {t1.y(), t1.z(), t1.snp(), t1.tgl(), t1.signed1Pt()};
      std::array<float, 15> t1covs = {t1.cYY(), t1.cZY(), t1.cZZ(), t1.cSnpY(), t1.cSnpZ(),
                                      t1.cSnpSnp(), t1.cTglY(), t1.cTglZ(), t1.cTglSnp(), t1.cTglTgl(),
                                      t1.c1PtY(), t1.c1PtZ(), t1.c1PtSnp(), t1.c1PtTgl(), t1.c1Pt21Pt2()};
      o2::track::TrackParCov pars1{t1.x(), t1.alpha(), t1pars, t1covs};
      std::array<float, 5> t2pars = {t2.y(), t2.z(), t2.snp(), t2.tgl(), t2.signed1Pt()};
      std::array<float, 15> t2covs = {t2.cYY(), t2.cZY(), t2.cZZ(), t2.cSnpY(), t2.cSnpZ(),
                                      t2.cSnpSnp(), t2.cTglY(), t2.cTglZ(), t2.cTglSnp(), t2.cTglTgl(),
                                      t2.c1PtY(), t2.c1PtZ(), t2.c1PtSnp(), t2.c1PtTgl(), t2.c1Pt21Pt2()};
      o2::track::TrackParCov pars2{t2.x(), t2.alpha(), t2pars, t2covs};
      std::array<float, 5> t3pars = {t3.y(), t3.z(), t3.snp(), t3.tgl(), t3.signed1Pt()};
      std::array<float, 15> t3covs = {t3.cYY(), t3.cZY(), t3.cZZ(), t3.cSnpY(), t3.cSnpZ(),
                                      t3.cSnpSnp(), t3.cTglY(), t3.cTglZ(), t3.cTglSnp(), t3.cTglTgl(),
                                      t3.c1PtY(), t3.c1PtZ(), t3.c1PtSnp(), t3.c1PtTgl(), t3.c1Pt21Pt2()};
      o2::track::TrackParCov pars3{t3.x(), t3.alpha(), t3pars, t3covs};
      procCode = VarManager::fgFitterThreeProngBarrel.process(pars1, pars2, pars3);
    } else {
      return;
    }

    values[VarManager::kVertexingProcCode] = procCode;
    if (procCode == 0) {
      // TODO: set the other variables to appropriate values and return
      values[kVertexingLxy] = -999.;
      values[kVertexingLxyz] = -999.;
      values[kVertexingLz] = -999.;
      values[kVertexingLxyErr] = -999.;
      values[kVertexingLxyzErr] = -999.;
      values[kVertexingLzErr] = -999.;

      values[kVertexingTauxy] = -999.;
      values[kVertexingTauz] = -999.;
      values[kVertexingTauxyErr] = -999.;
      values[kVertexingTauzErr] = -999.;

      values[kVertexingLzProjected] = -999.;
      values[kVertexingLxyProjected] = -999.;
      values[kVertexingLxyzProjected] = -999.;
      values[kVertexingTauzProjected] = -999.;
      values[kVertexingTauxyProjected] = -999.;
      values[kVertexingTauxyzProjected] = -999.;

      return;
    }

    Vec3D secondaryVertex;

    if constexpr (eventHasVtxCov) {
      secondaryVertex = fgFitterThreeProngBarrel.getPCACandidate();

      std::array<float, 6> covMatrixPCA = fgFitterThreeProngBarrel.calcPCACovMatrixFlat();

      o2::math_utils::Point3D<float> vtxXYZ(collision.posX(), collision.posY(), collision.posZ());
      std::array<float, 6> vtxCov{collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ()};
      o2::dataformats::VertexBase primaryVertex = {std::move(vtxXYZ), std::move(vtxCov)};
      auto covMatrixPV = primaryVertex.getCov();

      double phi = std::atan2(secondaryVertex[1] - collision.posY(), secondaryVertex[0] - collision.posX());
      double theta = std::atan2(secondaryVertex[2] - collision.posZ(),
                                std::sqrt((secondaryVertex[0] - collision.posX()) * (secondaryVertex[0] - collision.posX()) +
                                          (secondaryVertex[1] - collision.posY()) * (secondaryVertex[1] - collision.posY())));

      values[kVertexingLxy] = (collision.posX() - secondaryVertex[0]) * (collision.posX() - secondaryVertex[0]) +
                              (collision.posY() - secondaryVertex[1]) * (collision.posY() - secondaryVertex[1]);
      values[kVertexingLz] = (collision.posZ() - secondaryVertex[2]) * (collision.posZ() - secondaryVertex[2]);
      values[kVertexingLxyz] = values[kVertexingLxy] + values[kVertexingLz];
      values[kVertexingLxy] = std::sqrt(values[kVertexingLxy]);
      values[kVertexingLz] = std::sqrt(values[kVertexingLz]);
      values[kVertexingLxyz] = std::sqrt(values[kVertexingLxyz]);

      values[kVertexingLxyzErr] = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, theta) + getRotatedCovMatrixXX(covMatrixPCA, phi, theta));
      values[kVertexingLxyErr] = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, 0.) + getRotatedCovMatrixXX(covMatrixPCA, phi, 0.));
      values[kVertexingLzErr] = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, 0, theta) + getRotatedCovMatrixXX(covMatrixPCA, 0, theta));

      values[kVertexingTauz] = (collision.posZ() - secondaryVertex[2]) * v123.M() / (TMath::Abs(v123.Pz()) * o2::constants::physics::LightSpeedCm2NS);
      values[kVertexingTauxy] = values[kVertexingLxy] * v123.M() / (v123.Pt() * o2::constants::physics::LightSpeedCm2NS);

      values[kVertexingTauzErr] = values[kVertexingLzErr] * v123.M() / (TMath::Abs(v123.Pz()) * o2::constants::physics::LightSpeedCm2NS);
      values[kVertexingTauxyErr] = values[kVertexingLxyErr] * v123.M() / (v123.Pt() * o2::constants::physics::LightSpeedCm2NS);

      // run 2 definitions: Decay length projected onto the momentum vector of the candidate
      values[kVertexingLzProjected] = (secondaryVertex[2] - collision.posZ()) * v123.Pz();
      values[kVertexingLzProjected] = values[kVertexingLzProjected] / TMath::Sqrt(v123.Pz() * v123.Pz());
      values[kVertexingLxyProjected] = ((secondaryVertex[0] - collision.posX()) * v123.Px()) + ((secondaryVertex[1] - collision.posY()) * v123.Py());
      values[kVertexingLxyProjected] = values[kVertexingLxyProjected] / TMath::Sqrt((v123.Px() * v123.Px()) + (v123.Py() * v123.Py()));
      values[kVertexingLxyzProjected] = ((secondaryVertex[0] - collision.posX()) * v123.Px()) + ((secondaryVertex[1] - collision.posY()) * v123.Py()) + ((secondaryVertex[2] - collision.posZ()) * v123.Pz());
      values[kVertexingLxyzProjected] = values[kVertexingLxyzProjected] / TMath::Sqrt((v123.Px() * v123.Px()) + (v123.Py() * v123.Py()) + (v123.Pz() * v123.Pz()));

      values[kVertexingTauzProjected] = values[kVertexingLzProjected] * v123.M() / TMath::Abs(v123.Pz());
      values[kVertexingTauxyProjected] = values[kVertexingLxyProjected] * v123.M() / (v123.Pt());
      values[kVertexingTauxyzProjected] = values[kVertexingLxyzProjected] * v123.M() / (v123.P());
    }
  } else {
    KFParticle trk0KF;
    KFParticle trk1KF;
    KFParticle trk2KF;
    KFParticle KFGeoThreeProng;
    if ((tripletType == kTripleCandidateToKPiPi) && trackHasCov) {
      KFPTrack kfpTrack0 = createKFPTrackFromTrack(t1);
      trk0KF = KFParticle(kfpTrack0, 321 * t1.sign());
      KFPTrack kfpTrack1 = createKFPTrackFromTrack(t2);
      trk1KF = KFParticle(kfpTrack1, 211 * t2.sign());
      KFPTrack kfpTrack2 = createKFPTrackFromTrack(t3);
      trk2KF = KFParticle(kfpTrack2, 211 * t3.sign());

      KFGeoThreeProng.SetConstructMethod(3);
      KFGeoThreeProng.AddDaughter(trk0KF);
      KFGeoThreeProng.AddDaughter(trk1KF);
      KFGeoThreeProng.AddDaughter(trk2KF);

    } else if ((tripletType == kTripleCandidateToPKPi) && trackHasCov) {
      KFPTrack kfpTrack0 = createKFPTrackFromTrack(t1);
      trk0KF = KFParticle(kfpTrack0, 2212 * t1.sign());
      KFPTrack kfpTrack1 = createKFPTrackFromTrack(t2);
      trk1KF = KFParticle(kfpTrack1, 321 * t2.sign());
      KFPTrack kfpTrack2 = createKFPTrackFromTrack(t3);
      trk2KF = KFParticle(kfpTrack2, 211 * t3.sign());

      KFGeoThreeProng.SetConstructMethod(3);
      KFGeoThreeProng.AddDaughter(trk0KF);
      KFGeoThreeProng.AddDaughter(trk1KF);
      KFGeoThreeProng.AddDaughter(trk2KF);
    }
    if (fgUsedVars[kKFMass]) {
      float mass = 0., massErr = 0.;
      if (!KFGeoThreeProng.GetMass(mass, massErr))
        values[kKFMass] = mass;
      else
        values[kKFMass] = -999.;
    }

    if constexpr (eventHasVtxCov) {
      KFPVertex kfpVertex = createKFPVertexFromCollision(collision);
      values[kKFNContributorsPV] = kfpVertex.GetNContributors();
      KFParticle KFPV(kfpVertex);
      double dxTriplet2PV = KFGeoThreeProng.GetX() - KFPV.GetX();
      double dyTriplet2PV = KFGeoThreeProng.GetY() - KFPV.GetY();
      double dzTriplet2PV = KFGeoThreeProng.GetZ() - KFPV.GetZ();

      values[kVertexingLxy] = std::sqrt(dxTriplet2PV * dxTriplet2PV + dyTriplet2PV * dyTriplet2PV);
      values[kVertexingLz] = std::sqrt(dzTriplet2PV * dzTriplet2PV);
      values[kVertexingLxyz] = std::sqrt(dxTriplet2PV * dxTriplet2PV + dyTriplet2PV * dyTriplet2PV + dzTriplet2PV * dzTriplet2PV);

      values[kVertexingLxyErr] = (KFPV.GetCovariance(0) + KFGeoThreeProng.GetCovariance(0)) * dxTriplet2PV * dxTriplet2PV + (KFPV.GetCovariance(2) + KFGeoThreeProng.GetCovariance(2)) * dyTriplet2PV * dyTriplet2PV + 2 * ((KFPV.GetCovariance(1) + KFGeoThreeProng.GetCovariance(1)) * dxTriplet2PV * dyTriplet2PV);
      values[kVertexingLzErr] = (KFPV.GetCovariance(5) + KFGeoThreeProng.GetCovariance(5)) * dzTriplet2PV * dzTriplet2PV;
      values[kVertexingLxyzErr] = (KFPV.GetCovariance(0) + KFGeoThreeProng.GetCovariance(0)) * dxTriplet2PV * dxTriplet2PV + (KFPV.GetCovariance(2) + KFGeoThreeProng.GetCovariance(2)) * dyTriplet2PV * dyTriplet2PV + (KFPV.GetCovariance(5) + KFGeoThreeProng.GetCovariance(5)) * dzTriplet2PV * dzTriplet2PV + 2 * ((KFPV.GetCovariance(1) + KFGeoThreeProng.GetCovariance(1)) * dxTriplet2PV * dyTriplet2PV + (KFPV.GetCovariance(3) + KFGeoThreeProng.GetCovariance(3)) * dxTriplet2PV * dzTriplet2PV + (KFPV.GetCovariance(4) + KFGeoThreeProng.GetCovariance(4)) * dyTriplet2PV * dzTriplet2PV);
      if (fabs(values[kVertexingLxy]) < 1.e-8f)
        values[kVertexingLxy] = 1.e-8f;
      values[kVertexingLxyErr] = values[kVertexingLxyErr] < 0. ? 1.e8f : std::sqrt(values[kVertexingLxyErr]) / values[kVertexingLxy];
      if (fabs(values[kVertexingLz]) < 1.e-8f)
        values[kVertexingLz] = 1.e-8f;
      values[kVertexingLzErr] = values[kVertexingLzErr] < 0. ? 1.e8f : std::sqrt(values[kVertexingLzErr]) / values[kVertexingLz];
      if (fabs(values[kVertexingLxyz]) < 1.e-8f)
        values[kVertexingLxyz] = 1.e-8f;
      values[kVertexingLxyzErr] = values[kVertexingLxyzErr] < 0. ? 1.e8f : std::sqrt(values[kVertexingLxyzErr]) / values[kVertexingLxyz];

      values[kVertexingTauxy] = KFGeoThreeProng.GetPseudoProperDecayTime(KFPV, KFGeoThreeProng.GetMass()) / (o2::constants::physics::LightSpeedCm2NS);
      values[kVertexingTauz] = -1 * dzTriplet2PV * KFGeoThreeProng.GetMass() / (TMath::Abs(KFGeoThreeProng.GetPz()) * o2::constants::physics::LightSpeedCm2NS);
      values[kVertexingPz] = TMath::Abs(KFGeoThreeProng.GetPz());
      values[kVertexingSV] = KFGeoThreeProng.GetZ();
      values[kVertexingTauxyErr] = values[kVertexingLxyErr] * KFGeoThreeProng.GetMass() / (KFGeoThreeProng.GetPt() * o2::constants::physics::LightSpeedCm2NS);
      values[kVertexingTauzErr] = values[kVertexingLzErr] * KFGeoThreeProng.GetMass() / (TMath::Abs(KFGeoThreeProng.GetPz()) * o2::constants::physics::LightSpeedCm2NS);
      values[kCosPointingAngle] = (std::sqrt(dxTriplet2PV * dxTriplet2PV) * v123.Px() +
                                   std::sqrt(dyTriplet2PV * dyTriplet2PV) * v123.Py() +
                                   std::sqrt(dzTriplet2PV * dzTriplet2PV) * v123.Pz()) /
                                  (v123.P() * values[VarManager::kVertexingLxyz]);

      values[kVertexingLzProjected] = (dzTriplet2PV * KFGeoThreeProng.GetPz()) / TMath::Sqrt(KFGeoThreeProng.GetPz() * KFGeoThreeProng.GetPz());
      values[kVertexingLxyProjected] = (dxTriplet2PV * KFGeoThreeProng.GetPx()) + (dyTriplet2PV * KFGeoThreeProng.GetPy());
      values[kVertexingLxyProjected] = values[kVertexingLxyProjected] / TMath::Sqrt((KFGeoThreeProng.GetPx() * KFGeoThreeProng.GetPx()) + (KFGeoThreeProng.GetPy() * KFGeoThreeProng.GetPy()));
      values[kVertexingLxyzProjected] = (dxTriplet2PV * KFGeoThreeProng.GetPx()) + (dyTriplet2PV * KFGeoThreeProng.GetPy()) + (dzTriplet2PV * KFGeoThreeProng.GetPz());
      values[kVertexingLxyzProjected] = values[kVertexingLxyzProjected] / TMath::Sqrt((KFGeoThreeProng.GetPx() * KFGeoThreeProng.GetPx()) + (KFGeoThreeProng.GetPy() * KFGeoThreeProng.GetPy()) + (KFGeoThreeProng.GetPz() * KFGeoThreeProng.GetPz()));
      values[kVertexingTauxyProjected] = values[kVertexingLxyProjected] * KFGeoThreeProng.GetMass() / (KFGeoThreeProng.GetPt());
      values[kVertexingTauxyProjectedNs] = values[kVertexingTauxyProjected] / o2::constants::physics::LightSpeedCm2NS;
      values[kVertexingTauzProjected] = values[kVertexingLzProjected] * KFGeoThreeProng.GetMass() / TMath::Abs(KFGeoThreeProng.GetPz());
    }
  }
}

template <int candidateType, uint32_t collFillMap, uint32_t fillMap, typename C, typename T1>
void VarManager::FillDileptonTrackVertexing(C const& collision, T1 const& lepton1, T1 const& lepton2, T1 const& track, float* values)
{

  constexpr bool eventHasVtxCov = ((collFillMap & Collision) > 0 || (collFillMap & ReducedEventVtxCov) > 0);
  constexpr bool trackHasCov = ((fillMap & TrackCov) > 0 || (fillMap & ReducedTrackBarrelCov) > 0);
  constexpr bool muonHasCov = ((fillMap & MuonCov) > 0 || (fillMap & ReducedMuonCov) > 0);
  if (!values) {
    values = fgValues;
  }

  float mtrack;
  float mlepton1, mlepton2;

  int procCode = 0;
  int procCodeJpsi = 0;

  values[kUsedKF] = fgUsedKF;
  if (!fgUsedKF) {
    if constexpr ((candidateType == kBcToThreeMuons) && muonHasCov) {
      mlepton1 = o2::constants::physics::MassMuon;
      mlepton2 = o2::constants::physics::MassMuon;
      mtrack = o2::constants::physics::MassMuon;

      double chi21 = lepton1.chi2();
      double chi22 = lepton2.chi2();
      double chi23 = track.chi2();
      SMatrix5 t1pars(lepton1.x(), lepton1.y(), lepton1.phi(), lepton1.tgl(), lepton1.signed1Pt());
      std::vector<double> v1{lepton1.cXX(), lepton1.cXY(), lepton1.cYY(), lepton1.cPhiX(), lepton1.cPhiY(),
                             lepton1.cPhiPhi(), lepton1.cTglX(), lepton1.cTglY(), lepton1.cTglPhi(), lepton1.cTglTgl(),
                             lepton1.c1PtX(), lepton1.c1PtY(), lepton1.c1PtPhi(), lepton1.c1PtTgl(), lepton1.c1Pt21Pt2()};
      SMatrix55 t1covs(v1.begin(), v1.end());
      o2::track::TrackParCovFwd pars1{lepton1.z(), t1pars, t1covs, chi21};

      SMatrix5 t2pars(lepton2.x(), lepton2.y(), lepton2.phi(), lepton2.tgl(), lepton2.signed1Pt());
      std::vector<double> v2{lepton2.cXX(), lepton2.cXY(), lepton2.cYY(), lepton2.cPhiX(), lepton2.cPhiY(),
                             lepton2.cPhiPhi(), lepton2.cTglX(), lepton2.cTglY(), lepton2.cTglPhi(), lepton2.cTglTgl(),
                             lepton2.c1PtX(), lepton2.c1PtY(), lepton2.c1PtPhi(), lepton2.c1PtTgl(), lepton2.c1Pt21Pt2()};
      SMatrix55 t2covs(v2.begin(), v2.end());
      o2::track::TrackParCovFwd pars2{lepton2.z(), t2pars, t2covs, chi22};

      SMatrix5 t3pars(track.x(), track.y(), track.phi(), track.tgl(), track.signed1Pt());
      std::vector<double> v3{track.cXX(), track.cXY(), track.cYY(), track.cPhiX(), track.cPhiY(),
                             track.cPhiPhi(), track.cTglX(), track.cTglY(), track.cTglPhi(), track.cTglTgl(),
                             track.c1PtX(), track.c1PtY(), track.c1PtPhi(), track.c1PtTgl(), track.c1Pt21Pt2()};
      SMatrix55 t3covs(v3.begin(), v3.end());
      o2::track::TrackParCovFwd pars3{track.z(), t3pars, t3covs, chi23};
      procCode = VarManager::fgFitterThreeProngFwd.process(pars1, pars2, pars3);
      procCodeJpsi = VarManager::fgFitterTwoProngFwd.process(pars1, pars2);
    } else if constexpr ((candidateType == kBtoJpsiEEK || candidateType == kDstarToD0KPiPi) && trackHasCov) {
      if constexpr ((candidateType == kBtoJpsiEEK) && trackHasCov) {
        mlepton1 = o2::constants::physics::MassElectron;
        mlepton2 = o2::constants::physics::MassElectron;
        mtrack = o2::constants::physics::MassKaonCharged;
      } else if constexpr ((candidateType == kDstarToD0KPiPi) && trackHasCov) {
        mlepton1 = o2::constants::physics::MassKaonCharged;
        mlepton2 = o2::constants::physics::MassPionCharged;
        mtrack = o2::constants::physics::MassPionCharged;
      }
      std::array<float, 5> lepton1pars = {lepton1.y(), lepton1.z(), lepton1.snp(), lepton1.tgl(), lepton1.signed1Pt()};
      std::array<float, 15> lepton1covs = {lepton1.cYY(), lepton1.cZY(), lepton1.cZZ(), lepton1.cSnpY(), lepton1.cSnpZ(),
                                           lepton1.cSnpSnp(), lepton1.cTglY(), lepton1.cTglZ(), lepton1.cTglSnp(), lepton1.cTglTgl(),
                                           lepton1.c1PtY(), lepton1.c1PtZ(), lepton1.c1PtSnp(), lepton1.c1PtTgl(), lepton1.c1Pt21Pt2()};
      o2::track::TrackParCov pars1{lepton1.x(), lepton1.alpha(), lepton1pars, lepton1covs};
      std::array<float, 5> lepton2pars = {lepton2.y(), lepton2.z(), lepton2.snp(), lepton2.tgl(), lepton2.signed1Pt()};
      std::array<float, 15> lepton2covs = {lepton2.cYY(), lepton2.cZY(), lepton2.cZZ(), lepton2.cSnpY(), lepton2.cSnpZ(),
                                           lepton2.cSnpSnp(), lepton2.cTglY(), lepton2.cTglZ(), lepton2.cTglSnp(), lepton2.cTglTgl(),
                                           lepton2.c1PtY(), lepton2.c1PtZ(), lepton2.c1PtSnp(), lepton2.c1PtTgl(), lepton2.c1Pt21Pt2()};
      o2::track::TrackParCov pars2{lepton2.x(), lepton2.alpha(), lepton2pars, lepton2covs};
      std::array<float, 5> lepton3pars = {track.y(), track.z(), track.snp(), track.tgl(), track.signed1Pt()};
      std::array<float, 15> lepton3covs = {track.cYY(), track.cZY(), track.cZZ(), track.cSnpY(), track.cSnpZ(),
                                           track.cSnpSnp(), track.cTglY(), track.cTglZ(), track.cTglSnp(), track.cTglTgl(),
                                           track.c1PtY(), track.c1PtZ(), track.c1PtSnp(), track.c1PtTgl(), track.c1Pt21Pt2()};
      o2::track::TrackParCov pars3{track.x(), track.alpha(), lepton3pars, lepton3covs};
      procCode = VarManager::fgFitterThreeProngBarrel.process(pars1, pars2, pars3);
      procCodeJpsi = VarManager::fgFitterTwoProngBarrel.process(pars1, pars2);
    } else {
      return;
    }

    ROOT::Math::PtEtaPhiMVector v1(lepton1.pt(), lepton1.eta(), lepton1.phi(), mlepton1);
    ROOT::Math::PtEtaPhiMVector v2(lepton2.pt(), lepton2.eta(), lepton2.phi(), mlepton2);
    ROOT::Math::PtEtaPhiMVector v3(track.pt(), track.eta(), track.phi(), mtrack);
    ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
    ROOT::Math::PtEtaPhiMVector vdilepton(v12.pt(), v12.eta(), v12.phi(), v12.M());
    ROOT::Math::PtEtaPhiMVector v123 = vdilepton + v3;
    values[VarManager::kPairMass] = v123.M();
    values[VarManager::kMassDau] = mtrack;
    values[VarManager::kDeltaMass] = v123.M() - v12.M();
    values[VarManager::kPairPt] = v123.Pt();
    values[VarManager::kPairEta] = v123.Eta();
    if (fgUsedVars[kPairMassDau] || fgUsedVars[kPairPtDau]) {
      values[VarManager::kPairMassDau] = v12.M();
      values[VarManager::kPairPtDau] = v12.Pt();
    }
    values[VarManager::kPt] = track.pt();

    values[VarManager::kVertexingProcCode] = procCode;
    if (procCode == 0 || procCodeJpsi == 0) {
      // TODO: set the other variables to appropriate values and return
      values[VarManager::kVertexingChi2PCA] = -999.;
      values[VarManager::kVertexingLxy] = -999.;
      values[VarManager::kVertexingLxyz] = -999.;
      values[VarManager::kVertexingLz] = -999.;
      values[VarManager::kVertexingLxyErr] = -999.;
      values[VarManager::kVertexingLxyzErr] = -999.;
      values[VarManager::kVertexingLzErr] = -999.;

      values[VarManager::kVertexingTauxy] = -999.;
      values[VarManager::kVertexingTauz] = -999.;
      values[VarManager::kVertexingTauxyErr] = -999.;
      values[VarManager::kVertexingTauzErr] = -999.;
      values[VarManager::kVertexingPz] = -999.;
      values[VarManager::kVertexingSV] = -999.;
      return;
    }

    Vec3D secondaryVertex;

    if constexpr (eventHasVtxCov) {
      std::array<float, 6> covMatrixPCA;
      o2::dataformats::DCA impactParameter0;
      o2::dataformats::DCA impactParameter1;

      o2::math_utils::Point3D<float> vtxXYZ(collision.posX(), collision.posY(), collision.posZ());
      std::array<float, 6> vtxCov{collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ()};
      o2::dataformats::VertexBase primaryVertex = {std::move(vtxXYZ), std::move(vtxCov)};
      auto covMatrixPV = primaryVertex.getCov();

      if constexpr ((candidateType == kBtoJpsiEEK || candidateType == kDstarToD0KPiPi) && trackHasCov) {
        secondaryVertex = fgFitterThreeProngBarrel.getPCACandidate();
        covMatrixPCA = fgFitterThreeProngBarrel.calcPCACovMatrixFlat();
      } else if constexpr (candidateType == kBcToThreeMuons && muonHasCov) {
        secondaryVertex = fgFitterThreeProngFwd.getPCACandidate();
        covMatrixPCA = fgFitterThreeProngFwd.calcPCACovMatrixFlat();
      }

      auto chi2PCA = fgFitterThreeProngBarrel.getChi2AtPCACandidate();
      if (fgUsedVars[kVertexingChi2PCA])
        values[VarManager::kVertexingChi2PCA] = chi2PCA;

      double phi = std::atan2(secondaryVertex[1] - collision.posY(), secondaryVertex[0] - collision.posX());
      double theta = std::atan2(secondaryVertex[2] - collision.posZ(),
                                std::sqrt((secondaryVertex[0] - collision.posX()) * (secondaryVertex[0] - collision.posX()) +
                                          (secondaryVertex[1] - collision.posY()) * (secondaryVertex[1] - collision.posY())));
      if (fgUsedVars[kVertexingLxy] || fgUsedVars[kVertexingLz] || fgUsedVars[kVertexingLxyz]) {

        values[VarManager::kVertexingLxy] = (collision.posX() - secondaryVertex[0]) * (collision.posX() - secondaryVertex[0]) +
                                            (collision.posY() - secondaryVertex[1]) * (collision.posY() - secondaryVertex[1]);
        values[VarManager::kVertexingLxy] = std::sqrt(values[VarManager::kVertexingLxy]);
        values[VarManager::kVertexingLz] = (collision.posZ() - secondaryVertex[2]) * (collision.posZ() - secondaryVertex[2]);
        values[VarManager::kVertexingLz] = std::sqrt(values[VarManager::kVertexingLz]);
        values[VarManager::kVertexingLxyz] = values[VarManager::kVertexingLxy] + values[VarManager::kVertexingLz];
        values[VarManager::kVertexingLxyz] = std::sqrt(values[VarManager::kVertexingLxyz]);
      }

      if (fgUsedVars[kVertexingLxyzErr] || fgUsedVars[kVertexingLxyErr] || fgUsedVars[kVertexingLzErr]) {
        values[kVertexingLxyzErr] = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, theta) + getRotatedCovMatrixXX(covMatrixPCA, phi, theta));
        values[kVertexingLxyErr] = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, 0.) + getRotatedCovMatrixXX(covMatrixPCA, phi, 0.));
        values[kVertexingLzErr] = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, 0, theta) + getRotatedCovMatrixXX(covMatrixPCA, 0, theta));
      }

      values[kVertexingTauz] = (collision.posZ() - secondaryVertex[2]) * v123.M() / (TMath::Abs(v123.Pz()) * o2::constants::physics::LightSpeedCm2NS);
      values[kVertexingTauxy] = values[kVertexingLxy] * v123.M() / (v123.Pt() * o2::constants::physics::LightSpeedCm2NS);

      values[kVertexingPz] = TMath::Abs(v123.Pz());
      values[kVertexingSV] = secondaryVertex[2];

      values[kVertexingTauzErr] = values[kVertexingLzErr] * v123.M() / (TMath::Abs(v123.Pz()) * o2::constants::physics::LightSpeedCm2NS);
      values[kVertexingTauxyErr] = values[kVertexingLxyErr] * v123.M() / (v123.Pt() * o2::constants::physics::LightSpeedCm2NS);

      if (fgUsedVars[kCosPointingAngle] && fgUsedVars[kVertexingLxyz]) {
        values[VarManager::kCosPointingAngle] = ((collision.posX() - secondaryVertex[0]) * v123.Px() +
                                                 (collision.posY() - secondaryVertex[1]) * v123.Py() +
                                                 (collision.posZ() - secondaryVertex[2]) * v123.Pz()) /
                                                (v123.P() * values[VarManager::kVertexingLxyz]);
      }
      // run 2 definitions: Lxy projected onto the momentum vector of the candidate
      if (fgUsedVars[kVertexingLxyProjected] || fgUsedVars[kVertexingLxyzProjected] || values[kVertexingTauxyProjected]) {
        values[kVertexingLzProjected] = (secondaryVertex[2] - collision.posZ()) * v123.Pz();
        values[kVertexingLzProjected] = values[kVertexingLzProjected] / TMath::Sqrt(v123.Pz() * v123.Pz());
        values[kVertexingLxyProjected] = ((secondaryVertex[0] - collision.posX()) * v123.Px()) + ((secondaryVertex[1] - collision.posY()) * v123.Py());
        values[kVertexingLxyProjected] = values[kVertexingLxyProjected] / TMath::Sqrt((v123.Px() * v123.Px()) + (v123.Py() * v123.Py()));
        values[kVertexingLxyzProjected] = ((secondaryVertex[0] - collision.posX()) * v123.Px()) + ((secondaryVertex[1] - collision.posY()) * v123.Py()) + ((secondaryVertex[2] - collision.posZ()) * v123.Pz());
        values[kVertexingLxyzProjected] = values[kVertexingLxyzProjected] / TMath::Sqrt((v123.Px() * v123.Px()) + (v123.Py() * v123.Py()) + (v123.Pz() * v123.Pz()));
        values[kVertexingTauzProjected] = values[kVertexingLzProjected] * v123.M() / TMath::Abs(v123.Pz());
        values[kVertexingTauxyProjected] = values[kVertexingLxyProjected] * v123.M() / v123.Pt();
      }
    }
  } else {
    KFParticle lepton1KF;
    KFParticle lepton2KF;
    KFParticle hadronKF;
    KFParticle KFGeoTwoLeptons;
    KFParticle KFGeoThreeProng;

    if constexpr ((candidateType == kBtoJpsiEEK) && trackHasCov) {
      KFPTrack kfpTrack0 = createKFPTrackFromTrack(lepton1);
      lepton1KF = KFParticle(kfpTrack0, -11 * lepton1.sign());
      KFPTrack kfpTrack1 = createKFPTrackFromTrack(lepton2);
      lepton2KF = KFParticle(kfpTrack1, -11 * lepton2.sign());
      KFPTrack kfpTrack2 = createKFPTrackFromTrack(track);
      hadronKF = KFParticle(kfpTrack2, 321 * track.sign()); // kaon mass

      KFGeoTwoLeptons.SetConstructMethod(2);
      KFGeoTwoLeptons.AddDaughter(lepton1KF);
      KFGeoTwoLeptons.AddDaughter(lepton2KF);

      if (fgUsedVars[kPairMassDau] || fgUsedVars[kPairPtDau]) {
        values[VarManager::kPairMassDau] = KFGeoTwoLeptons.GetMass();
        values[VarManager::kPairPtDau] = KFGeoTwoLeptons.GetPt();
      }

      // Quantities between 3rd prong and candidate
      if (fgUsedVars[kKFDCAxyzBetweenProngs])
        values[kKFDCAxyzBetweenProngs] = KFGeoTwoLeptons.GetDistanceFromParticle(hadronKF);

      KFGeoThreeProng.SetConstructMethod(2);
      KFGeoThreeProng.AddDaughter(KFGeoTwoLeptons);
      KFGeoThreeProng.AddDaughter(hadronKF);

      if (fgUsedVars[kKFMass])
        values[kKFMass] = KFGeoThreeProng.GetMass();

      if constexpr (eventHasVtxCov) {
        KFPVertex kfpVertex = createKFPVertexFromCollision(collision);
        KFParticle KFPV(kfpVertex);
        double dxTriplet3PV = KFGeoThreeProng.GetX() - KFPV.GetX();
        double dyTriplet3PV = KFGeoThreeProng.GetY() - KFPV.GetY();
        double dzTriplet3PV = KFGeoThreeProng.GetZ() - KFPV.GetZ();

        if (fgUsedVars[kVertexingLxy] || fgUsedVars[kVertexingLz] || fgUsedVars[kVertexingLxyz] || fgUsedVars[kVertexingLxyErr] || fgUsedVars[kVertexingLzErr] || fgUsedVars[kVertexingTauxy] || fgUsedVars[kVertexingLxyOverErr] || fgUsedVars[kVertexingLzOverErr] || fgUsedVars[kVertexingLxyzOverErr] || fgUsedVars[kCosPointingAngle]) {
          values[kVertexingLxy] = std::sqrt(dxTriplet3PV * dxTriplet3PV + dyTriplet3PV * dyTriplet3PV);
          values[kVertexingLz] = std::sqrt(dzTriplet3PV * dzTriplet3PV);
          values[kVertexingLxyz] = std::sqrt(dxTriplet3PV * dxTriplet3PV + dyTriplet3PV * dyTriplet3PV + dzTriplet3PV * dzTriplet3PV);
          values[kVertexingLxyErr] = (KFPV.GetCovariance(0) + KFGeoThreeProng.GetCovariance(0)) * dxTriplet3PV * dxTriplet3PV + (KFPV.GetCovariance(2) + KFGeoThreeProng.GetCovariance(2)) * dyTriplet3PV * dyTriplet3PV + 2 * ((KFPV.GetCovariance(1) + KFGeoThreeProng.GetCovariance(1)) * dxTriplet3PV * dyTriplet3PV);
          values[kVertexingLzErr] = (KFPV.GetCovariance(5) + KFGeoThreeProng.GetCovariance(5)) * dzTriplet3PV * dzTriplet3PV;
          values[kVertexingLxyzErr] = (KFPV.GetCovariance(0) + KFGeoThreeProng.GetCovariance(0)) * dxTriplet3PV * dxTriplet3PV + (KFPV.GetCovariance(2) + KFGeoThreeProng.GetCovariance(2)) * dyTriplet3PV * dyTriplet3PV + (KFPV.GetCovariance(5) + KFGeoThreeProng.GetCovariance(5)) * dzTriplet3PV * dzTriplet3PV + 2 * ((KFPV.GetCovariance(1) + KFGeoThreeProng.GetCovariance(1)) * dxTriplet3PV * dyTriplet3PV + (KFPV.GetCovariance(3) + KFGeoThreeProng.GetCovariance(3)) * dxTriplet3PV * dzTriplet3PV + (KFPV.GetCovariance(4) + KFGeoThreeProng.GetCovariance(4)) * dyTriplet3PV * dzTriplet3PV);
          if (fabs(values[kVertexingLxy]) < 1.e-8f)
            values[kVertexingLxy] = 1.e-8f;
          values[kVertexingLxyErr] = values[kVertexingLxyErr] < 0. ? 1.e8f : std::sqrt(values[kVertexingLxyErr]) / values[kVertexingLxy];
          if (fabs(values[kVertexingLz]) < 1.e-8f)
            values[kVertexingLz] = 1.e-8f;
          values[kVertexingLzErr] = values[kVertexingLzErr] < 0. ? 1.e8f : std::sqrt(values[kVertexingLzErr]) / values[kVertexingLz];
          if (fabs(values[kVertexingLxyz]) < 1.e-8f)
            values[kVertexingLxyz] = 1.e-8f;
          values[kVertexingLxyzErr] = values[kVertexingLxyzErr] < 0. ? 1.e8f : std::sqrt(values[kVertexingLxyzErr]) / values[kVertexingLxyz];

          if (fgUsedVars[kVertexingTauxy])
            values[kVertexingTauxy] = KFGeoThreeProng.GetPseudoProperDecayTime(KFPV, KFGeoThreeProng.GetMass()) / (o2::constants::physics::LightSpeedCm2NS);
          if (fgUsedVars[kVertexingTauxyErr])
            values[kVertexingTauxyErr] = values[kVertexingLxyErr] * KFGeoThreeProng.GetMass() / (KFGeoThreeProng.GetPt() * o2::constants::physics::LightSpeedCm2NS);

          if (fgUsedVars[kCosPointingAngle])
            values[VarManager::kCosPointingAngle] = (dxTriplet3PV * KFGeoThreeProng.GetPx() +
                                                     dyTriplet3PV * KFGeoThreeProng.GetPy() +
                                                     dzTriplet3PV * KFGeoThreeProng.GetPz()) /
                                                    (KFGeoThreeProng.GetP() * values[VarManager::kVertexingLxyz]);
        } // end calculate vertex variables

        // As defined in Run 2 (projected onto momentum)
        if (fgUsedVars[kVertexingLxyProjected] || fgUsedVars[kVertexingLxyzProjected] || fgUsedVars[kVertexingLzProjected]) {
          values[kVertexingLzProjected] = (dzTriplet3PV * KFGeoThreeProng.GetPz()) / TMath::Sqrt(KFGeoThreeProng.GetPz() * KFGeoThreeProng.GetPz());
          values[kVertexingLxyProjected] = (dxTriplet3PV * KFGeoThreeProng.GetPx()) + (dyTriplet3PV * KFGeoThreeProng.GetPy());
          values[kVertexingLxyProjected] = values[kVertexingLxyProjected] / TMath::Sqrt((KFGeoThreeProng.GetPx() * KFGeoThreeProng.GetPx()) + (KFGeoThreeProng.GetPy() * KFGeoThreeProng.GetPy()));
          values[kVertexingLxyzProjected] = (dxTriplet3PV * KFGeoThreeProng.GetPx()) + (dyTriplet3PV * KFGeoThreeProng.GetPy()) + (dzTriplet3PV * KFGeoThreeProng.GetPz());
          values[kVertexingLxyzProjected] = values[kVertexingLxyzProjected] / TMath::Sqrt((KFGeoThreeProng.GetPx() * KFGeoThreeProng.GetPx()) + (KFGeoThreeProng.GetPy() * KFGeoThreeProng.GetPy()) + (KFGeoThreeProng.GetPz() * KFGeoThreeProng.GetPz()));
          values[kVertexingTauxyProjected] = (values[kVertexingLxyProjected] * KFGeoThreeProng.GetMass()) / (KFGeoThreeProng.GetPt());
          values[kVertexingTauxyProjectedNs] = values[kVertexingTauxyProjected] / o2::constants::physics::LightSpeedCm2NS;
          values[kVertexingTauzProjected] = (values[kVertexingLzProjected] * KFGeoThreeProng.GetMass()) / TMath::Abs(KFGeoThreeProng.GetPz());
        } // end Run 2 quantities
      }   // end eventHasVtxCov
    }     // end (candidateType == kBtoJpsiEEK) && trackHasCov
  }       // end KF
}

template <typename C, typename A>
void VarManager::FillQVectorFromGFW(C const& /*collision*/, A const& compA11, A const& compB11, A const& compC11, A const& compA21, A const& compB21, A const& compC21, A const& compA31, A const& compB31, A const& compC31, A const& compA41, A const& compB41, A const& compC41, A const& compA23, A const& compA42, float S10A, float S10B, float S10C, float S11A, float S11B, float S11C, float S12A, float S13A, float S14A, float S21A, float S22A, float S31A, float S41A, float* values)
{
  if (!values) {
    values = fgValues;
  }

  // Fill Qn vectors from generic flow framework for different eta gap A, B, C (n=1,2,3,4) with proper normalisation
  // Use normalized Q-vectors for SP and EP
  values[kQ1X0A] = compA11.real() / S11A;
  values[kQ1Y0A] = compA11.imag() / S11A;
  values[kQ1X0B] = compB11.real() / S11B;
  values[kQ1Y0B] = compB11.imag() / S11B;
  values[kQ1X0C] = compC11.real() / S11C;
  values[kQ1Y0C] = compC11.imag() / S11C;
  values[kQ2X0A] = compA21.real() / S11A;
  values[kQ2Y0A] = compA21.imag() / S11A;
  values[kQ2X0B] = compB21.real() / S11B;
  values[kQ2Y0B] = compB21.imag() / S11B;
  values[kQ2X0C] = compC21.real() / S11C;
  values[kQ2Y0C] = compC21.imag() / S11C;
  values[kQ3X0A] = compA31.real() / S11A;
  values[kQ3Y0A] = compA31.imag() / S11A;
  values[kQ3X0B] = compB31.real() / S11B;
  values[kQ3Y0B] = compB31.imag() / S11B;
  values[kQ3X0C] = compC31.real() / S11C;
  values[kQ3Y0C] = compC31.imag() / S11C;
  values[kQ4X0A] = compA41.real() / S11A;
  values[kQ4Y0A] = compA41.imag() / S11A;
  values[kQ4X0B] = compB41.real() / S11B;
  values[kQ4Y0B] = compB41.imag() / S11B;
  values[kQ4X0C] = compC41.real() / S11C;
  values[kQ4Y0C] = compC41.imag() / S11C;
  values[kQ42XA] = compA42.real(); // Only being used by cumulants, no need for normalization
  values[kQ42YA] = compA42.imag(); // Only being used by cumulants, no need for normalization
  values[kQ23XA] = compA23.real(); // Only being used by cumulants, no need for normalization
  values[kQ23YA] = compA23.imag(); // Only being used by cumulants, no need for normalization
  values[kS11A] = S11A;
  values[kS12A] = S12A;
  values[kS13A] = S13A;
  values[kS31A] = S31A;

  // Q-vectors components correlation (A, B, C)
  values[kQ2YYAB] = values[kQ2Y0A] * values[kQ2Y0B];
  values[kQ2XXAB] = values[kQ2X0A] * values[kQ2X0B];
  values[kQ2XYAB] = values[kQ2X0A] * values[kQ2Y0B];
  values[kQ2YXAB] = values[kQ2Y0A] * values[kQ2X0B];
  values[kQ2YYAC] = values[kQ2Y0A] * values[kQ2Y0C];
  values[kQ2XXAC] = values[kQ2X0A] * values[kQ2X0C];
  values[kQ2XYAC] = values[kQ2X0A] * values[kQ2Y0C];
  values[kQ2YXAC] = values[kQ2Y0A] * values[kQ2X0C];
  values[kQ2YYBC] = values[kQ2Y0B] * values[kQ2Y0C];
  values[kQ2XXBC] = values[kQ2X0B] * values[kQ2X0C];
  values[kQ2XYBC] = values[kQ2X0B] * values[kQ2Y0C];
  values[kQ2YXBC] = values[kQ2Y0B] * values[kQ2X0C];

  // Fill event multiplicities
  values[kMultA] = S10A;
  values[kMultB] = S10B;
  values[kMultC] = S10C;

  // Fill necessary quantities for cumulant calculations with weighted Q-vectors
  values[kM11REF] = S21A - S12A;
  values[kM1111REF] = S41A - 6. * S12A * S21A + 8. * S13A * S11A + 3. * S22A - 6. * S14A;
  values[kCORR2REF] = (norm(compA21) - S12A) / values[kM11REF];
  values[kCORR4REF] = (pow(norm(compA21), 2) + norm(compA42) - 2. * (compA42 * conj(compA21) * conj(compA21)).real() + 8. * (compA23 * conj(compA21)).real() - 4. * S12A * norm(compA21) - 6. * S14A - 2. * S22A) / values[kM1111REF];

  // TODO: provide different computations for R
  // Compute the R factor using the 2 sub-events technique for second and third harmonic
  // Compute event planes
  auto Psi2A = getEventPlane(2, values[kQ2X0A], values[kQ2Y0A]);
  auto Psi2B = getEventPlane(2, values[kQ2X0B], values[kQ2Y0B]);
  auto Psi3B = getEventPlane(3, values[kQ3X0B], values[kQ3Y0B]);
  auto Psi2C = getEventPlane(2, values[kQ2X0C], values[kQ2Y0C]);
  auto Psi3C = getEventPlane(3, values[kQ3X0C], values[kQ3Y0C]);
  values[kPsi2A] = Psi2A;
  values[kPsi2B] = Psi2B;
  values[kPsi2C] = Psi2C;

  values[kR2SP_AB] = (values[kQ2X0A] * values[kQ2X0B] + values[kQ2Y0A] * values[kQ2Y0B]);
  values[kR2SP_AC] = (values[kQ2X0A] * values[kQ2X0C] + values[kQ2Y0A] * values[kQ2Y0C]);
  values[kR2SP_BC] = (values[kQ2X0B] * values[kQ2X0C] + values[kQ2Y0B] * values[kQ2Y0C]);
  values[kR3SP] = (values[kQ3X0B] * values[kQ3X0C] + values[kQ3Y0B] * values[kQ3Y0C]);

  values[kR2EP_AB] = TMath::Cos(2 * (Psi2A - Psi2B));
  values[kR2EP_AC] = TMath::Cos(2 * (Psi2A - Psi2C));
  values[kR2EP_BC] = TMath::Cos(2 * (Psi2B - Psi2C));
  values[kR3EP] = TMath::Cos(3 * (Psi3B - Psi3C));
}

template <typename C>
void VarManager::FillQVectorFromCentralFW(C const& collision, float* values)
{
  if (!values) {
    values = fgValues;
  }

  float xQVecFT0a = collision.qvecFT0ARe(); // already normalised
  float yQVecFT0a = collision.qvecFT0AIm(); // already normalised
  float xQVecFT0c = collision.qvecFT0CRe(); // already normalised
  float yQVecFT0c = collision.qvecFT0CIm(); // already normalised
  float xQVecFT0m = collision.qvecFT0MRe(); // already normalised
  float yQVecFT0m = collision.qvecFT0MIm(); // already normalised
  float xQVecFV0a = collision.qvecFV0ARe(); // already normalised
  float yQVecFV0a = collision.qvecFV0AIm(); // already normalised
  float xQVecBPos = collision.qvecTPCposRe(); // already normalised
  float yQVecBPos = collision.qvecTPCposIm(); // already normalised
  float xQVecBNeg = collision.qvecTPCnegRe(); // already normalised
  float yQVecBNeg = collision.qvecTPCnegIm(); // already normalised

  values[kQ2X0A] = collision.qvecTPCallRe();
  values[kQ2Y0A] = collision.qvecTPCallIm();
  values[kQ2X0APOS] = xQVecBPos;
  values[kQ2Y0APOS] = yQVecBPos;
  values[kQ2X0ANEG] = xQVecBNeg;
  values[kQ2Y0ANEG] = yQVecBNeg;
  values[kQ2X0B] = xQVecFT0a;
  values[kQ2Y0B] = yQVecFT0a;
  values[kQ2X0C] = xQVecFT0c;
  values[kQ2Y0C] = yQVecFT0c;
  values[kMultA] = collision.nTrkTPCpos() + collision.nTrkTPCneg();
  values[kMultAPOS] = collision.nTrkTPCpos();
  values[kMultANEG] = collision.nTrkTPCneg();
  values[kMultB] = collision.sumAmplFT0A(); // Be careful, this is weighted sum of multiplicity
  values[kMultC] = collision.sumAmplFT0C(); // Be careful, this is weighted sum of multiplicity

  // Q-vectors components correlation (A, B, C)
  values[kQ2YYAB] = values[kQ2Y0A] * values[kQ2Y0B];
  values[kQ2XXAB] = values[kQ2X0A] * values[kQ2X0B];
  values[kQ2XYAB] = values[kQ2X0A] * values[kQ2Y0B];
  values[kQ2YXAB] = values[kQ2Y0A] * values[kQ2X0B];
  values[kQ2YYAC] = values[kQ2Y0A] * values[kQ2Y0C];
  values[kQ2XXAC] = values[kQ2X0A] * values[kQ2X0C];
  values[kQ2XYAC] = values[kQ2X0A] * values[kQ2Y0C];
  values[kQ2YXAC] = values[kQ2Y0A] * values[kQ2X0C];
  values[kQ2YYBC] = values[kQ2Y0B] * values[kQ2Y0C];
  values[kQ2XXBC] = values[kQ2X0B] * values[kQ2X0C];
  values[kQ2XYBC] = values[kQ2X0B] * values[kQ2Y0C];
  values[kQ2YXBC] = values[kQ2Y0B] * values[kQ2X0C];

  EventPlaneHelper epHelper;
  float Psi2A = epHelper.GetEventPlane(values[kQ2X0A], values[kQ2Y0A], 2);
  float Psi2APOS = epHelper.GetEventPlane(values[kQ2X0APOS], values[kQ2Y0APOS], 2);
  float Psi2ANEG = epHelper.GetEventPlane(values[kQ2X0ANEG], values[kQ2Y0ANEG], 2);
  float Psi2B = epHelper.GetEventPlane(values[kQ2X0B], values[kQ2Y0B], 2);
  float Psi2C = epHelper.GetEventPlane(values[kQ2X0C], values[kQ2Y0C], 2);

  values[kPsi2A] = Psi2A;
  values[kPsi2APOS] = Psi2APOS;
  values[kPsi2ANEG] = Psi2ANEG;
  values[kPsi2B] = Psi2B;
  values[kPsi2C] = Psi2C;

  values[kR2SP_AB] = (values[kQ2X0A] * values[kQ2X0B] + values[kQ2Y0A] * values[kQ2Y0B]);
  values[kR2SP_AC] = (values[kQ2X0A] * values[kQ2X0C] + values[kQ2Y0A] * values[kQ2Y0C]);
  values[kR2SP_BC] = (values[kQ2X0B] * values[kQ2X0C] + values[kQ2Y0B] * values[kQ2Y0C]);
  values[kR2SP_FT0CTPCPOS] = (xQVecFT0c * xQVecBPos + yQVecFT0c * yQVecBPos);
  values[kR2SP_FT0CTPCNEG] = (xQVecFT0c * xQVecBNeg + yQVecFT0c * yQVecBNeg);
  values[kR2SP_FT0ATPCPOS] = (xQVecFT0a * xQVecBPos + yQVecFT0a * yQVecBPos);
  values[kR2SP_FT0ATPCNEG] = (xQVecFT0a * xQVecBNeg + yQVecFT0a * yQVecBNeg);
  values[kR2SP_FT0MTPCPOS] = (xQVecFT0m * xQVecBPos + yQVecFT0m * yQVecBPos);
  values[kR2SP_FT0MTPCNEG] = (xQVecFT0m * xQVecBNeg + yQVecFT0m * yQVecBNeg);
  values[kR2SP_FV0ATPCPOS] = (xQVecFV0a * xQVecBPos + yQVecFV0a * yQVecBPos);
  values[kR2SP_FV0ATPCNEG] = (xQVecFV0a * xQVecBNeg + yQVecFV0a * yQVecBNeg);

  float epFT0a = epHelper.GetEventPlane(xQVecFT0a, yQVecFT0a, 2);
  float epFT0c = epHelper.GetEventPlane(xQVecFT0c, yQVecFT0c, 2);
  float epFT0m = epHelper.GetEventPlane(xQVecFT0m, yQVecFT0m, 2);
  float epFV0a = epHelper.GetEventPlane(xQVecFV0a, yQVecFV0a, 2);
  float epBPoss = epHelper.GetEventPlane(xQVecBPos, yQVecBPos, 2);
  float epBNegs = epHelper.GetEventPlane(xQVecBNeg, yQVecBNeg, 2);
  float epTPCFull = epHelper.GetEventPlane(values[kQ2X0A], values[kQ2Y0A], 2);

  values[kR2EP_AB] = std::cos(2 * getDeltaPsiInRange(epTPCFull, epFT0a, 2));
  values[kR2EP_AC] = std::cos(2 * getDeltaPsiInRange(epTPCFull, epFT0c, 2));
  values[kR2EP_BC] = std::cos(2 * getDeltaPsiInRange(epFT0a, epFT0c, 2));
  values[kR2EP_FT0CTPCPOS] = std::cos(2 * getDeltaPsiInRange(epFT0c, epBPoss, 2));
  values[kR2EP_FT0CTPCNEG] = std::cos(2 * getDeltaPsiInRange(epFT0c, epBNegs, 2));
  values[kR2EP_FT0ATPCPOS] = std::cos(2 * getDeltaPsiInRange(epFT0a, epBPoss, 2));
  values[kR2EP_FT0ATPCNEG] = std::cos(2 * getDeltaPsiInRange(epFT0a, epBNegs, 2));
  values[kR2EP_FT0MTPCPOS] = std::cos(2 * getDeltaPsiInRange(epFT0m, epBPoss, 2));
  values[kR2EP_FT0MTPCNEG] = std::cos(2 * getDeltaPsiInRange(epFT0m, epBNegs, 2));
  values[kR2EP_FV0ATPCPOS] = std::cos(2 * getDeltaPsiInRange(epFV0a, epBPoss, 2));
  values[kR2EP_FV0ATPCNEG] = std::cos(2 * getDeltaPsiInRange(epFV0a, epBNegs, 2));
}

template <typename C>
void VarManager::FillSpectatorPlane(C const& collision, float* values)
{
  if (!values) {
    values = fgValues;
  }

  auto zncEnergy = collision.energySectorZNC();
  auto znaEnergy = collision.energySectorZNA();
  for (int i = 0; i < 4; i++) { // avoid std::numeric_limits<float>::infinity() in the table
    if (zncEnergy[i] < -1.e12) {
      zncEnergy[i] = -1.f;
    }
    if (znaEnergy[i] < -1.e12) {
      znaEnergy[i] = -1.f;
    }
  }
  float znaCommon = collision.energyCommonZNA() < 0 ? -1.f : collision.energyCommonZNA();
  float zncCommon = collision.energyCommonZNC() < 0 ? -1.f : collision.energyCommonZNC();
  float zpaCommon = collision.energyCommonZPA() < 0 ? -1.f : collision.energyCommonZPA();
  float zpcCommon = collision.energyCommonZPC() < 0 ? -1.f : collision.energyCommonZPC();

  // Store ZNA and ZNC energies for calibrations
  values[kEnergyCommonZNA] = znaCommon;
  values[kEnergyCommonZNC] = zncCommon;
  values[kEnergyCommonZPA] = zpaCommon;
  values[kEnergyCommonZPC] = zpcCommon;
  values[kEnergyZNA1] = znaEnergy[0];
  values[kEnergyZNA2] = znaEnergy[1];
  values[kEnergyZNA3] = znaEnergy[2];
  values[kEnergyZNA4] = znaEnergy[3];
  values[kEnergyZNC1] = zncEnergy[0];
  values[kEnergyZNC2] = zncEnergy[1];
  values[kEnergyZNC3] = zncEnergy[2];
  values[kEnergyZNC4] = zncEnergy[3];
  values[kTimeZNA] = collision.timeZNA();
  values[kTimeZNC] = collision.timeZNC();
  values[kTimeZPA] = collision.timeZPA();
  values[kTimeZPC] = collision.timeZPC();

  constexpr float beamEne = 5.36 * 0.5;
  constexpr float x[4] = {-1.75, 1.75, -1.75, 1.75};
  constexpr float y[4] = {-1.75, -1.75, 1.75, 1.75};
  // constexpr float intcalibZNA[4] = {0.7997028, 0.8453715, 0.7879917, 0.7695486};
  // constexpr float intcalibZNC[4] = {0.7631577, 0.8408003, 0.7083920, 0.7731769};
  // constexpr float alpha = 0.395; // WARNING: Run 2 coorection, to be checked
  constexpr float alpha = 1.;
  float numXZNC = 0., numYZNC = 0., denZNC = 0.;
  float numXZNA = 0., numYZNA = 0., denZNA = 0.;

  float sumZNA = 0;
  float sumZNC = 0;

  for (int i = 0; i < 4; i++) {
    if (zncEnergy[i] > 0.) {
      float wZNC = std::pow(zncEnergy[i], alpha);
      // sumZNC += intcalibZNC[i] * wZNC;
      sumZNC += wZNC;
      numXZNC -= x[i] * wZNC;
      numYZNC += y[i] * wZNC;
      denZNC += wZNC;
    }
    if (znaEnergy[i] > 0.) {
      float wZNA = std::pow(znaEnergy[i], alpha);
      // sumZNA += intcalibZNA[i] * wZNA;
      sumZNA += wZNA;
      numXZNA += x[i] * wZNA;
      numYZNA += y[i] * wZNA;
      denZNA += wZNA;
    }
  }

  if (denZNC != 0.) {
    float nSpecnC = zncCommon / beamEne;                  // WARNING: Run 2 coorection, to be checked
    float cZNC = 1.89358 - 0.71262 / (nSpecnC + 0.71789); // WARNING: Run 2 coorection, to be checked
    cZNC = 1.;
    values[kQ1ZNCX] = cZNC * numXZNC / denZNC;
    values[kQ1ZNCY] = cZNC * numYZNC / denZNC;
  } else {
    values[kQ1ZNCX] = values[kQ1ZNCY] = 999.;
  }

  if (denZNA != 0.) {
    float nSpecnA = znaCommon / beamEne;                  // WARNING: Run 2 coorection, to be checked
    float cZNA = 1.89358 - 0.71262 / (nSpecnA + 0.71789); // WARNING: Run 2 coorection, to be checked
    cZNA = 1.;
    values[kQ1ZNAX] = cZNA * numXZNA / denZNA;
    values[kQ1ZNAY] = cZNA * numYZNA / denZNA;
  } else {
    values[kQ1ZNAX] = values[kQ1ZNAY] = 999.;
  }

  if (denZNA != 0. && denZNC != 0.) {
    values[kQ1ZNACXX] = values[kQ1ZNAX] * values[kQ1ZNCX];
    values[kQ1ZNACYY] = values[kQ1ZNAY] * values[kQ1ZNCY];
    values[kQ1ZNACYX] = values[kQ1ZNAY] * values[kQ1ZNCX];
    values[kQ1ZNACXY] = values[kQ1ZNAX] * values[kQ1ZNCY];
  } else {
    values[kQ1ZNACXX] = values[kQ1ZNACYY] = values[kQ1ZNACYX] = values[kQ1ZNACXY] = 999.;
  }

  if (znaCommon != 0 && sumZNA != 0 && zncCommon != 0 && sumZNC) {
    values[KIntercalibZNA] = znaCommon - sumZNA;
    values[KIntercalibZNC] = zncCommon - sumZNC;
  }
}

template <uint32_t fillMap, int pairType, typename T1, typename T2>
void VarManager::FillPairVn(T1 const& t1, T2 const& t2, float* values)
{

  if (!values) {
    values = fgValues;
  }

  float m1 = o2::constants::physics::MassElectron;
  float m2 = o2::constants::physics::MassElectron;
  if constexpr (pairType == kDecayToMuMu) {
    m1 = o2::constants::physics::MassMuon;
    m2 = o2::constants::physics::MassMuon;
  }

  if constexpr (pairType == kDecayToPiPi) {
    m1 = o2::constants::physics::MassPionCharged;
    m2 = o2::constants::physics::MassPionCharged;
  }

  if constexpr (pairType == kElectronMuon) {
    m2 = o2::constants::physics::MassMuon;
  }

  // Fill dilepton information
  ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), m1);
  ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), m2);
  ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;

  // TODO: provide different computations for vn
  // Compute the scalar product UQ using Q-vector from A, for second and third harmonic
  // Dilepton vn could be accessible after dividing this product with the R factor
  values[kU2Q2] = values[kQ2X0A] * std::cos(2 * v12.Phi()) + values[kQ2Y0A] * std::sin(2 * v12.Phi());
  values[kU3Q3] = values[kQ3X0A] * std::cos(3 * v12.Phi()) + values[kQ3Y0A] * std::sin(3 * v12.Phi());
  values[kR2SP_AB] = (values[kQ2X0A] * values[kQ2X0B] + values[kQ2Y0A] * values[kQ2Y0B]);
  values[kR2SP_AC] = (values[kQ2X0A] * values[kQ2X0C] + values[kQ2Y0A] * values[kQ2Y0C]);
  values[kR2SP_BC] = (values[kQ2X0B] * values[kQ2X0C] + values[kQ2Y0B] * values[kQ2Y0C]);
  values[kR3SP] = (values[kQ3X0B] * values[kQ3X0C] + values[kQ3Y0B] * values[kQ3Y0C]);

  float Psi2A = getEventPlane(2, values[kQ2X0A], values[kQ2Y0A]);
  float Psi3A = getEventPlane(3, values[kQ3X0A], values[kQ3Y0A]);
  float Psi2B = getEventPlane(2, values[kQ2X0B], values[kQ2Y0B]);
  float Psi3B = getEventPlane(3, values[kQ3X0B], values[kQ3Y0B]);
  float Psi2C = getEventPlane(2, values[kQ2X0C], values[kQ2Y0C]);
  float Psi3C = getEventPlane(3, values[kQ3X0C], values[kQ3Y0C]);
  values[kCos2DeltaPhi] = std::cos(2 * (v12.Phi() - Psi2A));
  values[kCos3DeltaPhi] = std::cos(3 * (v12.Phi() - Psi3A));
  values[kR2EP_AB] = TMath::Cos(2 * (Psi2A - Psi2B));
  values[kR2EP_AC] = TMath::Cos(2 * (Psi2A - Psi2C));
  values[kR2EP_BC] = TMath::Cos(2 * (Psi2B - Psi2C));
  values[kR3EP] = TMath::Cos(3 * (Psi3B - Psi3C));

  values[kCos2DeltaPhiMu1] = std::cos(2 * (v1.Phi() - v12.Phi()));
  values[kCos2DeltaPhiMu2] = std::cos(2 * (v2.Phi() - v12.Phi()));

  if (std::isnan(VarManager::fgValues[VarManager::kU2Q2]) == true) {
    values[kU2Q2] = -999.;
    values[kR2SP_AB] = -999.;
    values[kR2SP_AC] = -999.;
    values[kR2SP_BC] = -999.;
  }
  if (std::isnan(VarManager::fgValues[VarManager::kU3Q3]) == true) {
    values[kU3Q3] = -999.;
    values[kR3SP] = -999.;
  }
  if (std::isnan(VarManager::fgValues[VarManager::kCos2DeltaPhi]) == true) {
    values[kCos2DeltaPhi] = -999.;
    values[kR2EP_AB] = -999.;
    values[kR2EP_AC] = -999.;
    values[kR2EP_BC] = -999.;
  }
  if (std::isnan(VarManager::fgValues[VarManager::kCos3DeltaPhi]) == true) {
    values[kCos3DeltaPhi] = -999.;
    values[kR3EP] = -999.;
  }

  //  kV4, kC4POI, kC4REF etc.
  if constexpr ((fillMap & ReducedEventQvectorExtra) > 0) {
    complex<double> Q21(values[kQ2X0A] * values[kS11A], values[kQ2Y0A] * values[kS11A]);
    complex<double> Q42(values[kQ42XA], values[kQ42YA]);
    complex<double> Q23(values[kQ23XA], values[kQ23YA]);
    complex<double> P2(std::cos(2 * v12.Phi()), std::sin(2 * v12.Phi()));
    values[kM01POI] = values[kMultDimuons] * values[kS11A];
    values[kM0111POI] = values[kMultDimuons] * (values[kS31A] - 3. * values[kS11A] * values[kS12A] + 2. * values[kS13A]);
    values[kCORR2POI] = (P2 * conj(Q21)).real() / values[kM01POI];
    values[kCORR4POI] = (P2 * Q21 * conj(Q21) * conj(Q21) - P2 * Q21 * conj(Q42) - 2. * values[kS12A] * P2 * conj(Q21) + 2. * P2 * conj(Q23)).real() / values[kM0111POI];
    values[kM01POIoverMp] = values[kMultDimuons] > 0 && !(std::isnan(values[kM01POI]) || std::isinf(values[kM01POI])) ? values[kM01POI] / values[kMultDimuons] : 0;
    values[kM0111POIoverMp] = values[kMultDimuons] > 0 && !(std::isnan(values[kM0111POI]) || std::isinf(values[kM0111POI])) ? values[kM0111POI] / values[kMultDimuons] : 0;
    values[kM11REFoverMp] = values[kMultDimuons] > 0 && !(std::isnan(values[kM11REF]) || std::isinf(values[kM11REF])) ? values[kM11REF] / values[kMultDimuons] : 0;
    values[kM1111REFoverMp] = values[kMultDimuons] > 0 && !(std::isnan(values[kM1111REF]) || std::isinf(values[kM1111REF])) ? values[kM1111REF] / values[kMultDimuons] : 0;
    values[kCORR2POIMp] = std::isnan(values[kCORR2POI]) || std::isinf(values[kCORR2POI]) ? 0 : values[kCORR2POI] * values[kMultDimuons];
    values[kCORR4POIMp] = std::isnan(values[kCORR4POI]) || std::isinf(values[kCORR4POI]) ? 0 : values[kCORR4POI] * values[kMultDimuons];
    values[kCORR2REF] = std::isnan(values[kCORR2REF]) || std::isinf(values[kCORR2REF]) ? 0 : values[kCORR2REF];
    values[kCORR4REF] = std::isnan(values[kCORR4REF]) || std::isinf(values[kCORR4REF]) ? 0 : values[kCORR4REF];
  }
}

template <typename T>
void VarManager::FillZDC(T const& zdc, float* values)
{
  if (!values) {
    values = fgValues;
  }

  values[kEnergyCommonZNA] = (zdc.energyCommonZNA() > 0) ? zdc.energyCommonZNA() : -1.;
  values[kEnergyCommonZNC] = (zdc.energyCommonZNC() > 0) ? zdc.energyCommonZNC() : -1.;
  values[kEnergyCommonZPA] = (zdc.energyCommonZPA() > 0) ? zdc.energyCommonZPA() : -1.;
  values[kEnergyCommonZPC] = (zdc.energyCommonZPC() > 0) ? zdc.energyCommonZPC() : -1.;
  values[kTimeZNA] = zdc.timeZNA();
  values[kTimeZNC] = zdc.timeZNC();
  values[kTimeZPA] = zdc.timeZPA();
  values[kTimeZPC] = zdc.timeZPC();
}

template <typename T1, typename T2>
void VarManager::FillDileptonHadron(T1 const& dilepton, T2 const& hadron, float* values, float hadronMass)
{
  if (!values) {
    values = fgValues;
  }

  if (fgUsedVars[kPairMass] || fgUsedVars[kPairPt] || fgUsedVars[kPairEta] || fgUsedVars[kPairPhi] || fgUsedVars[kPairMassDau] || fgUsedVars[kPairPtDau]) {
    ROOT::Math::PtEtaPhiMVector v1(dilepton.pt(), dilepton.eta(), dilepton.phi(), dilepton.mass());
    ROOT::Math::PtEtaPhiMVector v2(hadron.pt(), hadron.eta(), hadron.phi(), hadronMass);
    ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
    values[kPairMass] = v12.M();
    values[kPairPt] = v12.Pt();
    values[kPairEta] = v12.Eta();
    values[kPairPhi] = v12.Phi();
    values[kPairMassDau] = dilepton.mass();
    values[kPairPtDau] = dilepton.pt();
    values[kMassDau] = hadronMass;
    values[kDeltaMass] = v12.M() - dilepton.mass();
  }
  if (fgUsedVars[kDeltaPhi]) {
    double delta = dilepton.phi() - hadron.phi();
    if (delta > 3.0 / 2.0 * M_PI) {
      delta -= 2.0 * M_PI;
    }
    if (delta < -0.5 * M_PI) {
      delta += 2.0 * M_PI;
    }
    values[kDeltaPhi] = delta;
  }
  if (fgUsedVars[kDeltaPhiSym]) {
    double delta = std::abs(dilepton.phi() - hadron.phi());
    if (delta > M_PI) {
      delta = 2 * M_PI - delta;
    }
    values[kDeltaPhiSym] = delta;
  }
  if (fgUsedVars[kDeltaEta]) {
    values[kDeltaEta] = dilepton.eta() - hadron.eta();
  }
}

template <typename T1, typename T2>
void VarManager::FillDileptonPhoton(T1 const& dilepton, T2 const& photon, float* values)
{
  if (!values) {
    values = fgValues;
  }
  if (fgUsedVars[kPairMass] || fgUsedVars[kPairPt] || fgUsedVars[kPairEta] || fgUsedVars[kPairPhi]) {
    ROOT::Math::PtEtaPhiMVector v1(dilepton.pt(), dilepton.eta(), dilepton.phi(), dilepton.mass());
    ROOT::Math::PtEtaPhiMVector v2(photon.pt(), photon.eta(), photon.phi(), photon.mGamma());
    ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
    values[kPairMass] = v12.M();
    values[kPairPt] = v12.Pt();
    values[kPairEta] = v12.Eta();
    values[kPairPhi] = v12.Phi();
    values[kPairMassDau] = dilepton.mass();
    values[kMassDau] = photon.mGamma();
    values[kPairPtDau] = dilepton.pt();
    values[kPt] = photon.pt();
    values[kDeltaEta] = dilepton.eta();
    values[kEta] = photon.eta();
    values[VarManager::kDeltaMass] = v12.M() - dilepton.mass();
    float m4 = o2::constants::physics::MassJPsi;
    values[VarManager::kDeltaMass_jpsi] = v12.M() - dilepton.mass() + m4;
    values[kRap] = v12.Rapidity();
  }
}

template <typename T>
void VarManager::FillHadron(T const& hadron, float* values, float hadronMass)
{
  if (!values) {
    values = fgValues;
  }

  ROOT::Math::PtEtaPhiMVector vhadron(hadron.pt(), hadron.eta(), hadron.phi(), hadronMass);
  values[kMass] = hadronMass;
  values[kPt] = hadron.pt();
  values[kEta] = hadron.eta();
  values[kPhi] = hadron.phi();
  values[kRap] = vhadron.Rapidity();
}

template <int partType, typename Cand, typename H, typename T>
void VarManager::FillSingleDileptonCharmHadron(Cand const& candidate, H hfHelper, T& bdtScoreCharmHad, float* values)
{
  if (!values) {
    values = fgValues;
  }

  if constexpr (partType == kJPsi) {
    values[kMass] = candidate.mass();
    values[kPt] = candidate.pt();
    values[kPhi] = candidate.phi();
    values[kRap] = candidate.rap();
  }
  if constexpr (partType == kD0ToPiK) {
    values[kMassCharmHadron] = hfHelper.invMassD0ToPiK(candidate);
    values[kPtCharmHadron] = candidate.pt();
    values[kPhiCharmHadron] = candidate.phi();
    values[kRapCharmHadron] = hfHelper.yD0(candidate);
    values[kBdtCharmHadron] = static_cast<float>(bdtScoreCharmHad);
  }
  if constexpr (partType == kD0barToKPi) {
    values[kMassCharmHadron] = hfHelper.invMassD0barToKPi(candidate);
    values[kPtCharmHadron] = candidate.pt();
    values[kPhiCharmHadron] = candidate.phi();
    values[kRapCharmHadron] = hfHelper.yD0(candidate);
    values[kBdtCharmHadron] = static_cast<float>(bdtScoreCharmHad);
  }
}

template <int partTypeCharmHad, typename DQ, typename HF, typename H, typename T>
void VarManager::FillDileptonCharmHadron(DQ const& dilepton, HF const& charmHadron, H hfHelper, T& bdtScoreCharmHad, float* values)
{
  FillSingleDileptonCharmHadron<kJPsi>(dilepton, hfHelper, bdtScoreCharmHad, values);
  FillSingleDileptonCharmHadron<partTypeCharmHad>(charmHadron, hfHelper, bdtScoreCharmHad, values);
}

template <int candidateType, typename T1, typename T2, typename T3>
void VarManager::FillDileptonTrackTrack(T1 const& dilepton, T2 const& hadron1, T3 const& hadron2, float* values)
{
  if (!values) {
    values = fgValues;
  }

  double defaultDileptonMass = 3.096;
  double hadronMass1 = o2::constants::physics::MassPionCharged;
  double hadronMass2 = o2::constants::physics::MassPionCharged;
  if (candidateType == kXtoJpsiPiPi) {
    defaultDileptonMass = 3.096;
    hadronMass1 = o2::constants::physics::MassPionCharged;
    hadronMass2 = o2::constants::physics::MassPionCharged;
  }
  if (candidateType == kChictoJpsiEE) {
    defaultDileptonMass = 3.096;
    hadronMass1 = o2::constants::physics::MassElectron;
    hadronMass2 = o2::constants::physics::MassElectron;
  }

  ROOT::Math::PtEtaPhiMVector v1(dilepton.pt(), dilepton.eta(), dilepton.phi(), dilepton.mass());
  ROOT::Math::PtEtaPhiMVector v2(hadron1.pt(), hadron1.eta(), hadron1.phi(), hadronMass1);
  ROOT::Math::PtEtaPhiMVector v3(hadron2.pt(), hadron2.eta(), hadron2.phi(), hadronMass2);
  ROOT::Math::PtEtaPhiMVector v123 = v1 + v2 + v3;
  values[kQuadMass] = v123.M();
  values[kQuadDefaultDileptonMass] = v123.M() - v1.M() + defaultDileptonMass;
  values[kQuadPt] = v123.Pt();
  values[kQuadEta] = v123.Eta();
  values[kQuadPhi] = v123.Phi();

  values[kTrackDCAxyProng1] = hadron1.dcaXY();
  values[kTrackDCAzProng1] = hadron1.dcaZ();
  values[kPt1] = hadron1.pt();

  values[kTrackDCAxyProng2] = hadron2.dcaXY();
  values[kTrackDCAzProng2] = hadron2.dcaZ();
  values[kPt2] = hadron2.pt();

  if (fgUsedVars[kCosthetaDileptonDitrack] || fgUsedVars[kPairMass] || fgUsedVars[kPairPt] || fgUsedVars[kDitrackPt] || fgUsedVars[kDitrackMass] || fgUsedVars[kQ] || fgUsedVars[kDeltaR1] || fgUsedVars[kDeltaR2] || fgUsedVars[kRap]) {
    ROOT::Math::PtEtaPhiMVector v23 = v2 + v3;
    values[kPairMass] = v1.M();
    values[kPairPt] = v1.Pt();
    values[kDitrackMass] = v23.M();
    values[kDitrackPt] = v23.Pt();
    values[kCosthetaDileptonDitrack] = (v1.Px() * v123.Px() + v1.Py() * v123.Py() + v1.Pz() * v123.Pz()) / (v1.P() * v123.P());
    values[kQ] = v123.M() - defaultDileptonMass - v23.M();
    values[kDeltaR1] = sqrt(pow(v1.Eta() - v2.Eta(), 2) + pow(v1.Phi() - v2.Phi(), 2));
    values[kDeltaR2] = sqrt(pow(v1.Eta() - v3.Eta(), 2) + pow(v1.Phi() - v3.Phi(), 2));
    values[kRap] = v123.Rapidity();
  }
}

//__________________________________________________________________
template <int pairType, typename T1, typename T2>
float VarManager::calculatePhiV(T1 const& t1, T2 const& t2)
{
  // cos(phiv) = w*a /|w||a|
  // with w = u x v
  // and  a = u x z / |u x z|   , unit vector perpendicular to v12 and z-direction (magnetic field)
  // u = v12 / |v12|            , the unit vector of v12
  // v = v1 x v2 / |v1 x v2|    , unit vector perpendicular to v1 and v2

  float m1 = o2::constants::physics::MassElectron;
  float m2 = o2::constants::physics::MassElectron;
  if constexpr (pairType == kDecayToMuMu) {
    m1 = o2::constants::physics::MassMuon;
    m2 = o2::constants::physics::MassMuon;
  }

  if constexpr (pairType == kDecayToPiPi) {
    m1 = o2::constants::physics::MassPionCharged;
    m2 = o2::constants::physics::MassPionCharged;
  }

  if constexpr (pairType == kElectronMuon) {
    m2 = o2::constants::physics::MassMuon;
  }

  ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), m1);
  ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), m2);
  ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;

  float pairPhiV = -999;
  float bz = fgMagField;

  bool swapTracks = false;
  if (v1.Pt() < v2.Pt()) { // ordering of track, pt1 > pt2
    ROOT::Math::PtEtaPhiMVector v3 = v1;
    v1 = v2;
    v2 = v3;
    swapTracks = true;
  }

  // momentum of e+ and e- in (ax,ay,az) axis. Note that az=0 by definition.
  // vector product of pep X pem
  float vpx = 0, vpy = 0, vpz = 0;
  if (t1.sign() * t2.sign() > 0) { // Like Sign
    if (!swapTracks) {
      if (bz * t1.sign() < 0) {
        vpx = v1.Py() * v2.Pz() - v1.Pz() * v2.Py();
        vpy = v1.Pz() * v2.Px() - v1.Px() * v2.Pz();
        vpz = v1.Px() * v2.Py() - v1.Py() * v2.Px();
      } else {
        vpx = v2.Py() * v1.Pz() - v2.Pz() * v1.Py();
        vpy = v2.Pz() * v1.Px() - v2.Px() * v1.Pz();
        vpz = v2.Px() * v1.Py() - v2.Py() * v1.Px();
      }
    } else { // swaped tracks
      if (bz * t2.sign() < 0) {
        vpx = v1.Py() * v2.Pz() - v1.Pz() * v2.Py();
        vpy = v1.Pz() * v2.Px() - v1.Px() * v2.Pz();
        vpz = v1.Px() * v2.Py() - v1.Py() * v2.Px();
      } else {
        vpx = v2.Py() * v1.Pz() - v2.Pz() * v1.Py();
        vpy = v2.Pz() * v1.Px() - v2.Px() * v1.Pz();
        vpz = v2.Px() * v1.Py() - v2.Py() * v1.Px();
      }
    }
  } else { // Unlike Sign
    if (!swapTracks) {
      if (bz * t1.sign() > 0) {
        vpx = v1.Py() * v2.Pz() - v1.Pz() * v2.Py();
        vpy = v1.Pz() * v2.Px() - v1.Px() * v2.Pz();
        vpz = v1.Px() * v2.Py() - v1.Py() * v2.Px();
      } else {
        vpx = v2.Py() * v1.Pz() - v2.Pz() * v1.Py();
        vpy = v2.Pz() * v1.Px() - v2.Px() * v1.Pz();
        vpz = v2.Px() * v1.Py() - v2.Py() * v1.Px();
      }
    } else { // swaped tracks
      if (bz * t2.sign() > 0) {
        vpx = v1.Py() * v2.Pz() - v1.Pz() * v2.Py();
        vpy = v1.Pz() * v2.Px() - v1.Px() * v2.Pz();
        vpz = v1.Px() * v2.Py() - v1.Py() * v2.Px();
      } else {
        vpx = v2.Py() * v1.Pz() - v2.Pz() * v1.Py();
        vpy = v2.Pz() * v1.Px() - v2.Px() * v1.Pz();
        vpz = v2.Px() * v1.Py() - v2.Py() * v1.Px();
      }
    }
  }

  // unit vector of pep X pem
  float vx = vpx / TMath::Sqrt(vpx * vpx + vpy * vpy + vpz * vpz);
  float vy = vpy / TMath::Sqrt(vpx * vpx + vpy * vpy + vpz * vpz);
  float vz = vpz / TMath::Sqrt(vpx * vpx + vpy * vpy + vpz * vpz);

  float px = v12.Px();
  float py = v12.Py();
  float pz = v12.Pz();

  // unit vector of (pep+pem)
  float ux = px / TMath::Sqrt(px * px + py * py + pz * pz);
  float uy = py / TMath::Sqrt(px * px + py * py + pz * pz);
  float uz = pz / TMath::Sqrt(px * px + py * py + pz * pz);
  float ax = uy / TMath::Sqrt(ux * ux + uy * uy);
  float ay = -ux / TMath::Sqrt(ux * ux + uy * uy);

  // The third axis defined by vector product (ux,uy,uz)X(vx,vy,vz)
  float wx = uy * vz - uz * vy;
  float wy = uz * vx - ux * vz;
  // by construction, (wx,wy,wz) must be a unit vector. Measure angle between (wx,wy,wz) and (ax,ay,0).
  // The angle between them should be small if the pair is conversion. This function then returns values close to pi!
  pairPhiV = TMath::ACos(wx * ax + wy * ay); // phiv in [0,pi] //cosPhiV = wx * ax + wy * ay;
  return pairPhiV;
}

#endif // PWGDQ_CORE_VARMANAGER_H_
