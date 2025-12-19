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

#ifndef HomogeneousField
#define HomogeneousField
#endif

#include "PWGUD/Core/UDHelpers.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/TriggerAliases.h"
#include "Common/Core/CollisionTypeHelper.h"
#include "Common/Core/EventPlaneHelper.h"
#include "Common/Core/PID/PIDTOFParamService.h"
#include "Common/Core/fwdtrackUtilities.h"
#include "Common/Core/trackUtilities.h"

#include <CommonConstants/LHCConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <DCAFitter/DCAFitterN.h>
#include <DCAFitter/FwdDCAFitterN.h>
#include <DetectorsBase/Propagator.h>
#include <Field/MagneticField.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/DataTypes.h>
#include <GlobalTracking/MatchGlobalFwd.h>
#include <ReconstructionDataFormats/DCA.h>
#include <ReconstructionDataFormats/Track.h>
#include <ReconstructionDataFormats/TrackFwd.h>
#include <ReconstructionDataFormats/Vertex.h>

#include <Math/GenVector/Boost.h>
#include <Math/SMatrix.h>
#include <Math/Vector3D.h>
#include <Math/Vector4D.h>
#include <Math/Vector4Dfwd.h>
#include <Math/VectorUtil.h>
#include <TGeoGlobalMagField.h>
#include <TH3F.h>
#include <THn.h>
#include <TObject.h>
#include <TRandom.h>
#include <TString.h>

#include <KFPTrack.h>
#include <KFPVertex.h>
#include <KFParticle.h>
#include <KFParticleBase.h>
#include <KFVertex.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <iostream>
#include <map>
#include <utility>
#include <vector>

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
    RapidityGapFilter = BIT(21),
    Fit = BIT(22),
    ReducedFit = BIT(23),
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
    Pair = BIT(17), // TODO: check whether we really need the Pair member here
    AmbiTrack = BIT(18),
    AmbiMuon = BIT(19),
    DalitzBits = BIT(20),
    TrackTPCPID = BIT(21),
    TrackMFT = BIT(22),
    ReducedTrackCollInfo = BIT(23), // TODO: remove it once new reduced data tables are produced for dielectron with ReducedTracksBarrelInfo
    ReducedMuonCollInfo = BIT(24),  // TODO: remove it once new reduced data tables are produced for dielectron with ReducedTracksBarrelInfo
    MuonRealign = BIT(25),
    MuonCovRealign = BIT(26),
    MFTCov = BIT(27),
    TrackTOFService = BIT(28),
    ParticleMC = BIT(29)
  };

  enum PairCandidateType {
    // TODO: need to agree on a scheme to incorporate all various hypotheses (e.g. e - mu, jpsi - K+, Jpsi - pipi,...)
    kDecayToEE = 0, // e.g. J/psi        -> e+ e-
    kDecayToMuMu,   // e.g. J/psi        -> mu+ mu-
    kDecayToPiPi,
    kElectronMuon,              // e.g. Electron - muon correlations
    kBcToThreeMuons,            // e.g. Bc           -> mu+ mu- mu+
    kBtoJpsiEEK,                // e.g. B+           -> e+ e- K+
    kJpsiEEProton,              // e.g. Jpsi-proton correlation, Jpsi to e+e-
    kXtoJpsiPiPi,               // e.g. X(3872)      -> J/psi pi+ pi-
    kPsi2StoJpsiPiPi,           // e.g. Psi(2S)      -> J/psi pi+ pi-
    kChictoJpsiEE,              // e.g. Chi_c1      -> J/psi e+ e-
    kDstarToD0KPiPi,            // e.g. D*+ -> D0 pi+ -> K- pi+ pi+
    kTripleCandidateToEEPhoton, // e.g. chi_c   -> e+ e- photon or pi0 -> e+ e- photon
    kDecayToKPi,                // e.g. D0           -> K+ pi- or cc.
    kTripleCandidateToKPiPi,    // e.g. D+ -> K- pi+ pi+
    kTripleCandidateToPKPi,     // e.g. Lambda_c -> p K- pi+
    kJpsiHadronMass,            // using the real hadron mass
    kJpsiPionMass,              // treat the hadron as pion
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
    kNRunWiseVariables,

    // Timeframe wise variables
    kTFNBCs = kNRunWiseVariables,
    kTFNCollisions,
    kTFNMCCollisions,
    kTFNTracks,
    kTFNMuons,
    kTFNMFTs,
    kNTFWiseVariables,

    // Event wise variables
    kTimestamp = kNTFWiseVariables,
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
    kIsSel8,                     // TVX in Run3 && No time frame border && No ITS read out frame border (from event selection)
    kIsGoodITSLayer3,            // number of inactive chips on ITS layer 3 is below maximum allowed value
    kIsGoodITSLayer0123,         // numbers of inactive chips on ITS layers 0-3 are below maximum allowed values
    kIsGoodITSLayersAll,         // numbers of inactive chips on all ITS layers are below maximum allowed values
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
    kCentFT0A,
    kCentFT0M,
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
    kMultAntiMuons,
    kMultMuons,
    kMultSingleMuons,
    kMultNTracksHasITS,
    kMultNTracksHasTPC,
    kMultNTracksHasTOF,
    kMultNTracksHasTRD,
    kMultNTracksITSOnly,
    kMultNTracksTPCOnly,
    kMultNTracksITSTPC,
    kMultNTracksPVeta1,
    kMultNTracksPVetaHalf,
    kTrackOccupancyInTimeRange,
    kFT0COccupancyInTimeRange,
    kNoCollInTimeRangeStandard,
    kMultAllTracksTPCOnly,
    kMultAllTracksITSTPC,
    kNTPCpileupContribA,
    kNTPCpileupContribC,
    kNTPCpileupZA,
    kNTPCpileupZC,
    kNTPCtracksInPast,
    kNTPCtracksInFuture,
    kNTPCcontribLongA,
    kNTPCcontribLongC,
    kNTPCmeanTimeLongA,
    kNTPCmeanTimeLongC,
    kNTPCmedianTimeLongA,
    kNTPCmedianTimeLongC,
    kNTPCcontribShortA,
    kNTPCcontribShortC,
    kNTPCmeanTimeShortA,
    kNTPCmeanTimeShortC,
    kNTPCmedianTimeShortA,
    kNTPCmedianTimeShortC,
    kMCEventGeneratorId,
    kMCEventSubGeneratorId,
    kMCVtxX,
    kMCVtxY,
    kMCVtxZ,
    kMCEventTime,
    kMCEventWeight,
    kMCEventImpParam,
    kMCEventCentrFT0C,
    kMultMCNParticlesEta10,
    kMultMCNParticlesEta08,
    kMultMCNParticlesEta05,
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
    kWR2SP_AB,
    kWR2SP_AC,
    kWR2SP_BC,
    kR2SP_AB_Im,
    kR2SP_AC_Im,
    kR2SP_BC_Im,
    kWR2SP_AB_Im,
    kWR2SP_AC_Im,
    kWR2SP_BC_Im,
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
    kWR2EP_AB,
    kWR2EP_AC,
    kWR2EP_BC,
    kR2EP_AB_Im,
    kR2EP_AC_Im,
    kR2EP_BC_Im,
    kWR2EP_AB_Im,
    kWR2EP_AC_Im,
    kWR2EP_BC_Im,
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
    kIsITSUPCMode, // UPC mode used for event
    kTwoEvPosZ1,   // vtx-z for collision 1 in two events correlations
    kTwoEvPosZ2,   // vtx-z for collision 2 in two events correlations
    kTwoEvPosR1,   // vtx-R for collision 1 in two events correlations
    kTwoEvPosR2,
    kTwoEvCentFT0C1,
    kTwoEvCentFT0C2,
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
    kQ2X0A1,
    kQ2X0A2,
    kQ2Y0A1,
    kQ2Y0A2,
    kU2Q2Ev1,
    kU2Q2Ev2,
    kCos2DeltaPhiEv1,
    kCos2DeltaPhiEv2,
    kV2SP1,
    kV2SP2,
    kV2EP1,
    kV2EP2,
    kV2ME_SP,
    kV2ME_EP,
    kWV2ME_SP,
    kWV2ME_EP,
    kTwoR2SP1, // Scalar product resolution of event1 for ME technique
    kTwoR2SP2, // Scalar product resolution of event2 for ME technique
    kTwoR2EP1, // Event plane resolution of event2 for ME technique
    kTwoR2EP2, // Event plane resolution of event2 for ME technique

    // Variables for event mixing with cumulant
    kV22m,
    kV24m,
    kV22p,
    kV24p,
    kV22ME,
    kV24ME,
    kWV22ME,
    kWV24ME,
    kNEventWiseVariables,

    // Basic track/muon/pair wise variables
    kX = kNEventWiseVariables,
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
    kPin_leg1,
    kTPCnSigmaKa_leg1,
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
    kTPCnCRoverFindCls,
    kTPCchi2,
    kTPCsignal,
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
    kTPCnSigmaMu,
    kTPCnSigmaPi,
    kTPCnSigmaKa,
    kTPCnSigmaPr,
    kTPCnSigmaEl_Corr,
    kTPCnSigmaPi_Corr,
    kTPCnSigmaKa_Corr,
    kTPCnSigmaPr_Corr,
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
    kMuonNClusters = kNBarrelTrackVariables,
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
    kMCMass,
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
    kMCHadronPdgCode,
    kMCCosTheta,
    kMCJpsiPt,
    kMCCosChi,
    kMCdeltaphi,
    kMCdeltaeta,
    kMCHadronPt,
    kMCHadronEta,
    kMCHadronPhi,
    kMCWeight,
    kMCCosChi_randomPhi_toward,
    kMCWeight_randomPhi_toward,
    kMCCosChi_randomPhi_away,
    kMCWeight_randomPhi_away,
    kMCCosChi_randomPhi_trans,
    kMCWeight_randomPhi_trans,
    kMCdeltaphi_randomPhi_toward,
    kMCdeltaphi_randomPhi_away,
    kMCdeltaphi_randomPhi_trans,
    kMCWeight_before,

    // MC mother particle variables
    kMCMotherPdgCode,

    // MC pair variables
    kMCCosThetaHE,
    kMCPhiHE,
    kMCPhiTildeHE,
    kMCCosThetaCS,
    kMCPhiCS,
    kMCPhiTildeCS,
    kMCCosThetaPP,
    kMCPhiPP,
    kMCPhiTildePP,
    kMCCosThetaRM,

    // Pair variables
    kCandidateId,
    kPairType,
    kVertexingLxy,
    kVertexingLxyErr,
    kVertexingPseudoCTau,
    kVertexingLxyz,
    kVertexingLxyzErr,
    kMCVertexingLxy,
    kMCVertexingLxyz,
    kMCLxyExpected,
    kMCLxyzExpected,
    kVertexingLz,
    kVertexingLzErr,
    kMCVertexingLz,
    kVertexingTauxy,
    kVertexingTauxyErr,
    kMCVertexingTauxy,
    kVertexingLzProjected,
    kVertexingLxyProjected,
    kVertexingLxyzProjected,
    kMCVertexingLzProjected,
    kMCVertexingLxyProjected,
    kMCVertexingLxyzProjected,
    kVertexingTauzProjected,
    kVertexingTauxyProjected,
    kVertexingTauxyProjectedPoleJPsiMass,
    kVertexingTauxyProjectedNs,
    kVertexingTauxyzProjected,
    kMCVertexingTauzProjected,
    kMCVertexingTauxyProjected,
    kMCVertexingTauxyProjectedNs,
    kMCVertexingTauxyzProjected,
    kMCCosPointingAngle,
    kVertexingTauz,
    kMCVertexingTauz,
    kVertexingTauzErr,
    kVertexingPz,
    kVertexingSV,
    kVertexingProcCode,
    kVertexingChi2PCA,
    kCosThetaHE,
    kPhiHE,
    kPhiTildeHE,
    kCosThetaCS,
    kPhiCS,
    kPhiTildeCS,
    kCosThetaPP,
    kPhiPP,
    kPhiTildePP,
    kCosThetaRM,
    kCosThetaStarTPC,
    kCosThetaStarFT0A,
    kCosThetaStarFT0C,
    kCosPhiVP,
    kPhiVP,
    kDeltaPhiPair2,
    kDeltaEtaPair2,
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
    kV2SP,
    kV2EP,
    kWV2SP,
    kWV2EP,
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
    kM11REFetagap,
    kM01POI,
    kM1111REF,
    kM11M1111REF,
    kM11M1111REFoverMp,
    kM01M0111POIoverMp,
    kM0111POI,
    kCORR2REF,
    kCORR2REFbydimuons,
    kCORR2REFbysinglemu,
    kCORR2REFetagap,
    kCORR2POI,
    kCORR2POICORR4POI,
    kCORR2REFCORR4POI,
    kCORR2REFCORR2POI,
    kM01M0111overMp,
    kM11M0111overMp,
    kM11M01overMp,
    kCORR2CORR4REF,
    kCORR4REF,
    kCORR4REFbydimuons,
    kCORR4REFbysinglemu,
    kCORR4POI,
    kM11REFoverMp,
    kM01POIoverMp,
    kM1111REFoverMp,
    kM0111POIoverMp,
    kCORR2POIMp,
    kCORR4POIMp,
    kM01POIplus,
    kM0111POIplus,
    kM01POIminus,
    kM0111POIminus,
    kM01POIsingle,
    kM0111POIsingle,
    kM01POIoverMpminus,
    kM01POIoverMpplus,
    kM01POIoverMpsingle,
    kM01POIoverMpmoins,
    kM0111POIoverMpminus,
    kM0111POIoverMpplus,
    kM0111POIoverMpsingle,
    kCORR2POIplus,
    kCORR2POIminus,
    kCORR2POIsingle,
    kCORR4POIplus,
    kCORR4POIminus,
    kCORR4POIsingle,
    kM11REFoverMpplus,
    kM1111REFoverMpplus,
    kM11REFoverMpminus,
    kM1111REFoverMpminus,
    kM11REFoverMpsingle,
    kM1111REFoverMpsingle,
    kM01POIME,
    kMultDimuonsME,
    kM0111POIME,
    kCORR2POIME,
    kCORR4POIME,
    kM01POIoverMpME,
    kM0111POIoverMpME,
    kM11REFoverMpME,
    kM1111REFoverMpME,
    kCORR2REFbydimuonsME,
    kCORR4REFbydimuonsME,
    kR2SP,
    kR2EP,
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
    kS12,
    kS13,
    kS23,
    kNPairVariables,

    // Candidate-track correlation variables
    kPairMass,
    kPairMassDau,
    kMassDau,
    kPairPt,
    kPairPtDau,
    kPairEta,
    kPairRap,
    kPairPhi,
    kPairPhiv,
    kDeltaEta,
    kDeltaPhi,
    kDeltaPhiSym,
    kNCorrelationVariables,
    kDileptonHadronKstar,
    kCosChi,
    kEtaDau,
    kPhiDau,
    kECWeight,
    kPtDau,
    kCosTheta,
    kEWeight_before,
    kCosChi_randomPhi_trans,
    kCosChi_randomPhi_toward,
    kCosChi_randomPhi_away,
    kWeight_randomPhi_trans,
    kWeight_randomPhi_toward,
    kWeight_randomPhi_away,
    kdeltaphi_randomPhi_trans,
    kdeltaphi_randomPhi_toward,
    kdeltaphi_randomPhi_away,

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
    kDeltaR,

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

    // BDT score
    kBdtBackground,
    kBdtPrompt,
    kBdtNonprompt,

    // FIT detector variables
    kAmplitudeFT0A,
    kAmplitudeFT0C,
    kTimeFT0A,
    kTimeFT0C,
    kTriggerMaskFT0,
    kAmplitudeFDDA,
    kAmplitudeFDDC,
    kTimeFDDA,
    kTimeFDDC,
    kTriggerMaskFDD,
    kAmplitudeFV0A,
    kTimeFV0A,
    kTriggerMaskFV0A,
    kBBFT0Apf,
    kBGFT0Apf,
    kBBFT0Cpf,
    kBGFT0Cpf,
    kBBFV0Apf,
    kBGFV0Apf,
    kBBFDDApf,
    kBGFDDApf,
    kBBFDDCpf,
    kBGFDDCpf,
    kNFiredChannelsFT0A,
    kNFiredChannelsFT0C,
    kNFiredChannelsFV0A,

    // ALICE 3 Variables
    kMultDensity,
    kMultMCNParticlesEta40,
    kMultMCNParticlesEta20,
    kIsReconstructed,
    kNSiliconHits,
    kNTPCHits,
    kOTTOTSignal,
    kOTnSigmaEl,
    kOTnSigmaMu,
    kOTnSigmaPi,
    kOTnSigmaKa,
    kOTnSigmaPr,
    kOTnSigmaDe,
    kOTnSigmaTr,
    kOTnSigmaHe3,
    kOTnSigmaAl,
    kHasRICHSig,
    kHasRICHSigInGas,
    kHasRICHSigEl,
    kHasRICHSigMu,
    kHasRICHSigPi,
    kHasRICHSigKa,
    kHasRICHSigPr,
    kHasRICHSigDe,
    kHasRICHSigTr,
    kHasRICHSigHe3,
    kHasRICHSigAl,
    kRICHnSigmaEl,
    kRICHnSigmaMu,
    kRICHnSigmaPi,
    kRICHnSigmaKa,
    kRICHnSigmaPr,
    kRICHnSigmaDe,
    kRICHnSigmaTr,
    kRICHnSigmaHe3,
    kRICHnSigmaAl,
    kTOFEventTime,
    kTOFEventTimeErr,
    kOuterTOFnSigmaEl,
    kOuterTOFnSigmaMu,
    kOuterTOFnSigmaPi,
    kOuterTOFnSigmaKa,
    kOuterTOFnSigmaPr,
    kOuterTOFnSigmaDe,
    kOuterTOFnSigmaTr,
    kOuterTOFnSigmaHe3,
    kOuterTOFnSigmaAl,
    kInnerTOFnSigmaEl,
    kInnerTOFnSigmaMu,
    kInnerTOFnSigmaPi,
    kInnerTOFnSigmaKa,
    kInnerTOFnSigmaPr,
    kInnerTOFnSigmaDe,
    kInnerTOFnSigmaTr,
    kInnerTOFnSigmaHe3,
    kInnerTOFnSigmaAl,
    kA3Variables,

    kNVars
  }; // end of Variables enumeration

  enum CalibObjects {
    kTPCElectronMean = 0,
    kTPCElectronSigma,
    kTPCElectronStatus,
    kTPCPionMean,
    kTPCPionSigma,
    kTPCPionStatus,
    kTPCKaonMean,
    kTPCKaonSigma,
    kTPCKaonStatus,
    kTPCProtonMean,
    kTPCProtonSigma,
    kTPCProtonStatus,
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
    kSingleGapC,
    kITSUPCMode
  };

  enum MuonExtrapolation {
    // Index used to set different options for Muon propagation
    kToVertex = 0, // propagtion to vertex by default
    kToDCA,
    kToRabs,
    kToMatching
  };

  static TString fgVariableNames[kNVars];      // variable names
  static TString fgVariableUnits[kNVars];      // variable units
  static std::map<TString, int> fgVarNamesMap; // key: variables short name, value: order in the Variables enum
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

  // Setup the collision system
  static void SetCollisionSystem(TString system, float energy);
  static void SetCollisionSystem(o2::parameters::GRPLHCIFData* grplhcif);

  static void SetMagneticField(float magField)
  {
    fgMagField = magField;
  }

  // Setup plane position for MFT-MCH matching
  static void SetMatchingPlane(float z)
  {
    fgzMatching = z;
  }

  static float GetMatchingPlane()
  {
    return fgzMatching;
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

  // Setup the 4 prong KFParticle
  static void SetupFourProngKFParticle(float magField)
  {
    KFParticle::SetField(magField);
    fgUsedKF = true;
  }

  // Setup the 4 prong DCAFitterN
  static void SetupFourProngDCAFitter(float magField, bool propagateToPCA, float maxR, float /*maxDZIni*/, float minParamChange, float minRelChi2Change, bool useAbsDCA)
  {
    fgFitterFourProngBarrel.setBz(magField);
    fgFitterFourProngBarrel.setPropagateToPCA(propagateToPCA);
    fgFitterFourProngBarrel.setMaxR(maxR);
    fgFitterFourProngBarrel.setMinParamChange(minParamChange);
    fgFitterFourProngBarrel.setMinRelChi2Change(minRelChi2Change);
    fgFitterFourProngBarrel.setUseAbsDCA(useAbsDCA);
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
  static o2::track::TrackParCovFwd FwdToTrackPar(const T& track, const C& cov);
  template <typename T, typename C>
  static o2::dataformats::GlobalFwdTrack PropagateMuon(const T& muon, const C& collision, int endPoint = kToVertex);
  template <typename T, typename C>
  static o2::track::TrackParCovFwd PropagateFwd(const T& track, const C& cov, float z);
  template <uint32_t fillMap, typename T, typename C>
  static void FillMuonPDca(const T& muon, const C& collision, float* values = nullptr);
  template <uint32_t fillMap, typename T, typename C>
  static void FillPropagateMuon(const T& muon, const C& collision, float* values = nullptr);
  template <typename T>
  static void FillBC(T const& bc, float* values = nullptr);
  template <uint32_t fillMap, typename T>
  static void FillEvent(T const& event, float* values = nullptr);
  template <typename T>
  static void FillTimeFrame(T const& tfTable, float* values = nullptr);
  template <typename T>
  static void FillEventFlowResoFactor(T const& hs_sp, T const& hs_ep, float* values = nullptr);
  template <typename T>
  static void FillTwoEvents(T const& event1, T const& event2, float* values = nullptr);
  template <uint32_t fillMap, typename T1, typename T2>
  static void FillTwoMixEvents(T1 const& event1, T1 const& event2, T2 const& tracks1, T2 const& tracks2, float* values = nullptr);
  template <typename T>
  static void FillTwoMixEventsFlowResoFactor(T const& hs_sp, T const& hs_ep, float* values = nullptr);
  template <typename T, typename T1, typename T2>
  static void FillTwoMixEventsCumulants(T const& h_v22m, T const& h_v24m, T const& h_v22p, T const& h_v24p, T1 const& t1, T2 const& t2, float* values = nullptr);
  template <uint32_t fillMap, typename T>
  static void FillTrack(T const& track, float* values = nullptr);
  template <uint32_t fillMap, typename T>
  static void FillPhoton(T const& photon, float* values = nullptr);
  template <uint32_t fillMap, typename T, typename C>
  static void FillTrackCollision(T const& track, C const& collision, float* values = nullptr);
  template <int candidateType, uint32_t fillMap, typename T1, typename T2, typename C>
  static void FillTrackCollisionMC(T1 const& track, T2 const& MotherTrack, C const& collision, float* values = nullptr);
  template <uint32_t fillMap, typename T, typename C, typename M, typename P>
  static void FillTrackCollisionMatCorr(T const& track, C const& collision, M const& materialCorr, P const& propagator, float* values = nullptr);
  template <typename U, typename T>
  static void FillTrackMC(const U& mcStack, T const& track, float* values = nullptr);
  template <int pairType, typename T, typename T1>
  static void FillEnergyCorrelatorsMC(T const& track, T1 const& t1, float* values = nullptr);
  template <uint32_t fillMap, typename T1, typename T2, typename C>
  static void FillPairPropagateMuon(T1 const& muon1, T2 const& muon2, const C& collision, float* values = nullptr);
  template <uint32_t fillMap, typename T1, typename T2, typename C>
  static void FillGlobalMuonRefit(T1 const& muontrack, T2 const& mfttrack, const C& collision, float* values = nullptr);
  template <uint32_t MuonfillMap, uint32_t MFTfillMap, typename T1, typename T2, typename C, typename C2>
  static void FillGlobalMuonRefitCov(T1 const& muontrack, T2 const& mfttrack, const C& collision, C2 const& mftcov, float* values = nullptr);
  template <int pairType, uint32_t fillMap, typename T1, typename T2>
  static void FillPair(T1 const& t1, T2 const& t2, float* values = nullptr);
  template <int pairType, uint32_t fillMap, typename C, typename T1, typename T2>
  static void FillPairCollision(C const& collision, T1 const& t1, T2 const& t2, float* values = nullptr);
  template <int pairType, uint32_t fillMap, typename C, typename T1, typename T2, typename M, typename P>
  static void FillPairCollisionMatCorr(C const& collision, T1 const& t1, T2 const& t2, M const& materialCorr, P const& propagator, float* values = nullptr);
  template <typename T1, typename T2, typename T3>
  static void FillTriple(T1 const& t1, T2 const& t2, T3 const& t3, float* values = nullptr, PairCandidateType pairType = kTripleCandidateToEEPhoton);
  template <uint32_t fillMap, int pairType, typename T1, typename T2>
  static void FillPairME(T1 const& t1, T2 const& t2, float* values = nullptr);
  template <int pairType, typename T1, typename T2>
  static void FillPairMC(T1 const& t1, T2 const& t2, float* values = nullptr);
  template <int candidateType, typename T1, typename T2, typename T3>
  static void FillTripleMC(T1 const& t1, T2 const& t2, T3 const& t3, float* values = nullptr);
  template <int candidateType, typename T1, typename T2>
  static void FillQuadMC(T1 const& t1, T2 const& t2, T2 const& t3, float* values = nullptr);
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
  template <int candidateType, uint32_t collFillMap, uint32_t fillMap, typename C, typename T1>
  static void FillDileptonTrackTrackVertexing(C const& collision, T1 const& lepton1, T1 const& lepton2, T1 const& track1, T1 const& track2, float* values);
  template <typename T>
  static void FillZDC(const T& zdc, float* values = nullptr);
  template <typename T>
  static void FillBdtScore(const T& bdtScore, float* values = nullptr);
  template <typename T1, typename T2, typename T3, typename T4, typename T5>
  static void FillFIT(const T1& bc, const T2& bcs, const T3& ft0s, const T4& fv0as, const T5& fdds, float* values = nullptr);
  template <int pairType, uint32_t fillMap, typename T1, typename T2>
  static void FillPairAlice3(T1 const& t1, T2 const& t2, float* values = nullptr);
  template <uint32_t fillMap, typename T>
  static void FillEventAlice3(T const& event, float* values = nullptr);
  template <uint32_t fillMap, typename T>
  static void FillTrackAlice3(T const& track, float* values  = nullptr);
  
  static void SetCalibrationObject(CalibObjects calib, TObject* obj)
  {
    fgCalibs[calib] = obj;
    // Check whether all the needed objects for TPC postcalibration are available
    if (fgCalibs.find(kTPCElectronMean) != fgCalibs.end() && fgCalibs.find(kTPCElectronSigma) != fgCalibs.end()) {
      fgRunTPCPostCalibration[0] = true;
      fgUsedVars[kTPCnSigmaEl_Corr] = true;
    }
    if (fgCalibs.find(kTPCPionMean) != fgCalibs.end() && fgCalibs.find(kTPCPionSigma) != fgCalibs.end()) {
      fgRunTPCPostCalibration[1] = true;
      fgUsedVars[kTPCnSigmaPi_Corr] = true;
    }
    if (fgCalibs.find(kTPCKaonMean) != fgCalibs.end() && fgCalibs.find(kTPCKaonSigma) != fgCalibs.end()) {
      fgRunTPCPostCalibration[2] = true;
      fgUsedVars[kTPCnSigmaKa_Corr] = true;
    }
    if (fgCalibs.find(kTPCProtonMean) != fgCalibs.end() && fgCalibs.find(kTPCProtonSigma) != fgCalibs.end()) {
      fgRunTPCPostCalibration[3] = true;
      fgUsedVars[kTPCnSigmaPr_Corr] = true;
    }
  }

  static void SetCalibrationType(int type, bool useInterpolation = true)
  {
    if (type < 0 || type > 2) {
      LOG(fatal) << "Invalid calibration type. Must be 0, 1, or 2.";
    }
    fgCalibrationType = type;
    fgUseInterpolatedCalibration = useInterpolation;
  }
  static double ComputePIDcalibration(int species, double nSigmaValue);

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
  static float fgzMatching;
  static float fgCenterOfMassEnergy;        // collision energy
  static float fgMassofCollidingParticle;   // mass of the colliding particle
  static float fgTPCInterSectorBoundary;    // TPC inter-sector border size at the TPC outer radius, in cm
  static int fgITSROFbias;                  // ITS ROF bias (from ALPIDE parameters)
  static int fgITSROFlength;                // ITS ROF length (from ALPIDE parameters)
  static int fgITSROFBorderMarginLow;       // ITS ROF border low margin
  static int fgITSROFBorderMarginHigh;      // ITS ROF border high margin
  static uint64_t fgSOR;                    // Timestamp for start of run
  static uint64_t fgEOR;                    // Timestamp for end of run
  static ROOT::Math::PxPyPzEVector fgBeamA; // beam from A-side 4-momentum vector
  static ROOT::Math::PxPyPzEVector fgBeamC; // beam from C-side 4-momentum vector

  // static void FillEventDerived(float* values = nullptr);
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
  template <typename T1, typename T2>
  static float LorentzTransformJpsihadroncosChi(TString Option, const T1& v1, const T2& v2);

  static o2::vertexing::DCAFitterN<2> fgFitterTwoProngBarrel;
  static o2::vertexing::DCAFitterN<3> fgFitterThreeProngBarrel;
  static o2::vertexing::DCAFitterN<4> fgFitterFourProngBarrel;
  static o2::vertexing::FwdDCAFitterN<2> fgFitterTwoProngFwd;
  static o2::vertexing::FwdDCAFitterN<3> fgFitterThreeProngFwd;
  static o2::globaltracking::MatchGlobalFwd mMatching;

  static std::map<CalibObjects, TObject*> fgCalibs; // map of calibration histograms
  static bool fgRunTPCPostCalibration[4];           // 0-electron, 1-pion, 2-kaon, 3-proton
  static int fgCalibrationType;                     // 0 - no calibration, 1 - calibration vs (TPCncls,pIN,eta) typically for pp, 2 - calibration vs (eta,nPV,nLong,tLong) typically for PbPb
  static bool fgUseInterpolatedCalibration;         // use interpolated calibration histograms (default: true)

  VarManager& operator=(const VarManager& c);
  VarManager(const VarManager& c);

  ClassDef(VarManager, 5);
};

template <typename T, typename C>
o2::track::TrackParCovFwd VarManager::FwdToTrackPar(const T& track, const C& cov)
{
  double chi2 = track.chi2();
  SMatrix5 tpars(track.x(), track.y(), track.phi(), track.tgl(), track.signed1Pt());
  std::vector<double> v1{cov.cXX(), cov.cXY(), cov.cYY(), cov.cPhiX(), cov.cPhiY(),
                         cov.cPhiPhi(), cov.cTglX(), cov.cTglY(), cov.cTglPhi(), cov.cTglTgl(),
                         cov.c1PtX(), cov.c1PtY(), cov.c1PtPhi(), cov.c1PtTgl(), cov.c1Pt21Pt2()};
  SMatrix55 tcovs(v1.begin(), v1.end());
  o2::track::TrackParCovFwd trackparCov{track.z(), tpars, tcovs, chi2};
  return trackparCov;
}

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
  o2::track::TrackParCovFwd trackparCov = FwdToTrackPar(muon, muon);

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
  kfpTrack.SetChi2(muon.chi2());
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
  o2::track::TrackParCovFwd fwdtrack = FwdToTrackPar(muon, muon);
  o2::dataformats::GlobalFwdTrack propmuon;
  if (static_cast<int>(muon.trackType()) > 2) {
    o2::dataformats::GlobalFwdTrack track;
    track.setParameters(fwdtrack.getParameters());
    track.setZ(fwdtrack.getZ());
    track.setCovariances(fwdtrack.getCovariances());
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
    if (endPoint == kToMatching) {
      o2::mch::TrackExtrap::extrapToVertexWithoutBranson(mchTrack, fgzMatching);
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

template <typename T, typename C>
o2::track::TrackParCovFwd VarManager::PropagateFwd(const T& track, const C& cov, float z)
{
  o2::track::TrackParCovFwd fwdtrack = FwdToTrackPar(track, cov);
  fwdtrack.propagateToZhelix(z, fgMagField);
  return fwdtrack;
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

  if constexpr ((fillMap & MuonCov) > 0 || (fillMap & ReducedMuonCov) > 0 || (fillMap & MuonCovRealign) > 0) {
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

template <uint32_t MuonfillMap, uint32_t MFTfillMap, typename T1, typename T2, typename C, typename C2>
void VarManager::FillGlobalMuonRefitCov(T1 const& muontrack, T2 const& mfttrack, const C& collision, C2 const& mftcov, float* values)
{
  if (!values) {
    values = fgValues;
  }
  if constexpr ((MuonfillMap & MuonCov) > 0) {
    if constexpr ((MFTfillMap & MFTCov) > 0) {
      o2::dataformats::GlobalFwdTrack propmuon = PropagateMuon(muontrack, collision);
      o2::track::TrackParCovFwd mft = FwdToTrackPar(mfttrack, mftcov);

      o2::dataformats::GlobalFwdTrack globalRefit = o2::aod::fwdtrackutils::refitGlobalMuonCov(propmuon, mft);
      values[kX] = globalRefit.getX();
      values[kY] = globalRefit.getY();
      values[kZ] = globalRefit.getZ();
      values[kTgl] = globalRefit.getTgl();
      values[kPt] = globalRefit.getPt();
      values[kPz] = globalRefit.getPz();
      values[kEta] = globalRefit.getEta();
      values[kPhi] = globalRefit.getPhi();
    }
  }
}

template <typename T>
void VarManager::FillTimeFrame(T const& tf, float* values)
{
  if (!values) {
    values = fgValues;
  }
  if constexpr (T::template contains<o2::aod::BCs>()) {
    values[kTFNBCs] = tf.size();
  }
  if constexpr (T::template contains<o2::aod::Collisions>()) {
    values[kTFNCollisions] = tf.size();
  }
  if constexpr (T::template contains<o2::aod::McCollisions>()) {
    values[kTFNMCCollisions] = tf.size();
  }
  if constexpr (T::template contains<o2::aod::Tracks>()) {
    values[kTFNTracks] = tf.size();
  }
  if constexpr (T::template contains<o2::aod::FwdTracks>()) {
    values[kTFNMuons] = tf.size();
  }
  if constexpr (T::template contains<o2::aod::MFTTracks>()) {
    values[kTFNMFTs] = tf.size();
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
    if (fgUsedVars[kTrackOccupancyInTimeRange]) {
      values[kTrackOccupancyInTimeRange] = event.trackOccupancyInTimeRange();
    }
    if (fgUsedVars[kFT0COccupancyInTimeRange]) {
      values[kFT0COccupancyInTimeRange] = event.ft0cOccupancyInTimeRange();
    }
    if (fgUsedVars[kNoCollInTimeRangeStandard]) {
      values[kNoCollInTimeRangeStandard] = event.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard);
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
    if (fgUsedVars[kIsGoodITSLayer3]) {
      values[kIsGoodITSLayer3] = event.selection_bit(o2::aod::evsel::kIsGoodITSLayer3);
    }
    if (fgUsedVars[kIsGoodITSLayer0123]) {
      values[kIsGoodITSLayer0123] = event.selection_bit(o2::aod::evsel::kIsGoodITSLayer0123);
    }
    if (fgUsedVars[kIsGoodITSLayersAll]) {
      values[kIsGoodITSLayersAll] = event.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll);
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
    if constexpr ((fillMap & CollisionMC) == 0) {
      values[kCentFT0C] = event.centFT0C();
      values[kCentFT0A] = event.centFT0A();
      values[kCentFT0M] = event.centFT0M();
    }
  }

  if constexpr ((fillMap & CollisionMult) > 0 || (fillMap & ReducedEventExtended) > 0) {
    values[kMultFV0A] = event.multFV0A();
    values[kMultFT0A] = event.multFT0A();
    values[kMultFT0C] = event.multFT0C();
    values[kMultFDDA] = event.multFDDA();
    values[kMultFDDC] = event.multFDDC();
    values[kMultTPC] = event.multTPC();
    values[kMultFV0C] = event.multFV0C();
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
    values[kMultNTracksPVeta1] = event.multNTracksPVeta1();
    values[kMultNTracksPVetaHalf] = event.multNTracksPVetaHalf();
    values[kMultAllTracksTPCOnly] = event.multAllTracksTPCOnly();
    values[kMultAllTracksITSTPC] = event.multAllTracksITSTPC();
    values[kTrackOccupancyInTimeRange] = event.trackOccupancyInTimeRange();
    values[kFT0COccupancyInTimeRange] = event.ft0cOccupancyInTimeRange();
    if constexpr ((fillMap & ReducedEventMultExtra) > 0) {
      values[kNTPCcontribLongA] = event.nTPCoccupContribLongA();
      values[kNTPCcontribLongC] = event.nTPCoccupContribLongC();
      values[kNTPCcontribShortA] = event.nTPCoccupContribShortA();
      values[kNTPCcontribShortC] = event.nTPCoccupContribShortC();
      values[kNTPCmeanTimeLongA] = event.nTPCoccupMeanTimeLongA();
      values[kNTPCmeanTimeLongC] = event.nTPCoccupMeanTimeLongC();
      values[kNTPCmeanTimeShortA] = event.nTPCoccupMeanTimeShortA();
      values[kNTPCmeanTimeShortC] = event.nTPCoccupMedianTimeShortC();
      values[kNTPCmedianTimeLongA] = event.nTPCoccupMedianTimeLongA();
      values[kNTPCmedianTimeLongC] = event.nTPCoccupMedianTimeLongC();
      values[kNTPCmedianTimeShortA] = event.nTPCoccupMedianTimeShortA();
      values[kNTPCmedianTimeShortC] = event.nTPCoccupMedianTimeShortC();
    }
  }
  // TODO: need to add EvSels and Cents tables, etc. in case of the central data model

  if constexpr ((fillMap & ReducedEvent) > 0) {
    values[kRunNo] = event.runNumber();
    values[kVtxX] = event.posX();
    values[kVtxY] = event.posY();
    values[kVtxZ] = event.posZ();
    values[kVtxNcontrib] = event.numContrib();
    if (fgUsedVars[kIsDoubleGap]) {
      values[kIsDoubleGap] = (event.tag_bit(56 + kDoubleGap) > 0);
    }
    if (fgUsedVars[kIsSingleGap] || fgUsedVars[kIsSingleGapA] || fgUsedVars[kIsSingleGapC]) {
      values[kIsSingleGapA] = (event.tag_bit(56 + kSingleGapA) > 0);
      values[kIsSingleGapC] = (event.tag_bit(56 + kSingleGapC) > 0);
      values[kIsSingleGap] = values[kIsSingleGapA] || values[kIsSingleGapC];
    }
    if (fgUsedVars[kIsITSUPCMode]) {
      values[kIsITSUPCMode] = (event.tag_bit(56 + kITSUPCMode) > 0);
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
    if (fgUsedVars[kIsNoITSROFBorder]) {
      values[kIsNoITSROFBorder] = (event.selection_bit(o2::aod::evsel::kNoITSROFrameBorder) > 0);
    }
    if (fgUsedVars[kIsNoTFBorder]) {
      values[kIsNoTFBorder] = (event.selection_bit(o2::aod::evsel::kNoTimeFrameBorder) > 0);
    }
    if (fgUsedVars[kNoCollInTimeRangeStandard]) {
      values[kNoCollInTimeRangeStandard] = (event.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard) > 0);
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
    if (fgUsedVars[kIsGoodITSLayer3]) {
      values[kIsGoodITSLayer3] = event.selection_bit(o2::aod::evsel::kIsGoodITSLayer3);
    }
    if (fgUsedVars[kIsGoodITSLayer0123]) {
      values[kIsGoodITSLayer0123] = event.selection_bit(o2::aod::evsel::kIsGoodITSLayer0123);
    }
    if (fgUsedVars[kIsGoodITSLayersAll]) {
      values[kIsGoodITSLayersAll] = event.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll);
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
      values[kMultA] = event.multa();
      values[kCORR2REFetagap] = event.corr2refetagap();
      values[kM11REFetagap] = event.m11refetagap();
      values[kCORR2REF] = std::isnan(values[kM11REF]) || std::isinf(values[kM11REF]) || std::isnan(values[kCORR2REF]) || std::isinf(values[kCORR2REF]) || std::isnan(values[kM1111REF]) || std::isinf(values[kM1111REF]) || std::isnan(values[kCORR4REF]) || std::isinf(values[kCORR4REF]) ? 0 : event.corr2ref();
      values[kCORR4REF] = std::isnan(values[kM1111REF]) || std::isinf(values[kM1111REF]) || std::isnan(values[kCORR4REF]) || std::isinf(values[kCORR4REF]) || std::isnan(values[kM11REF]) || std::isinf(values[kM11REF]) || std::isnan(values[kCORR2REF]) || std::isinf(values[kCORR2REF]) ? 0 : event.corr4ref();
      values[kCORR2CORR4REF] = std::isnan(values[kM11M1111REF]) || std::isinf(values[kM11M1111REF]) || std::isnan(values[kCORR2CORR4REF]) || std::isinf(values[kCORR2CORR4REF]) ? 0 : event.corr2ref() * event.corr4ref();
      values[kM11REF] = !(std::isnan(values[kM11REF]) || std::isinf(values[kM11REF]) || std::isnan(values[kCORR2REF]) || std::isinf(values[kCORR2REF]) || std::isnan(values[kM1111REF]) || std::isinf(values[kM1111REF]) || std::isnan(values[kCORR4REF]) || std::isinf(values[kCORR4REF])) ? event.m11ref() : 0;
      values[kM1111REF] = !(std::isnan(values[kM1111REF]) || std::isinf(values[kM1111REF]) || std::isnan(values[kCORR4REF]) || std::isinf(values[kCORR4REF]) || std::isnan(values[kM11REF]) || std::isinf(values[kM11REF]) || std::isnan(values[kCORR2REF]) || std::isinf(values[kCORR2REF])) ? event.m1111ref() : 0;
      values[kM11M1111REF] = !(std::isnan(values[kM11M1111REF]) || std::isinf(values[kM11M1111REF]) || std::isnan(values[kCORR2CORR4REF]) || std::isinf(values[kCORR2CORR4REF])) ? event.m11ref() * event.m1111ref() : 0;
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

    float R2SP_AB = (values[kQ2X0A] * values[kQ2X0B] + values[kQ2Y0A] * values[kQ2Y0B]);
    float R2SP_AC = (values[kQ2X0A] * values[kQ2X0C] + values[kQ2Y0A] * values[kQ2Y0C]);
    float R2SP_BC = (values[kQ2X0B] * values[kQ2X0C] + values[kQ2Y0B] * values[kQ2Y0C]);
    float R2SP_AB_Im = (values[kQ2Y0A] * values[kQ2X0B] - values[kQ2X0A] * values[kQ2Y0B]);
    float R2SP_AC_Im = (values[kQ2Y0A] * values[kQ2X0C] - values[kQ2X0A] * values[kQ2Y0C]);
    float R2SP_BC_Im = (values[kQ2Y0B] * values[kQ2X0C] - values[kQ2X0B] * values[kQ2Y0C]);
    values[kR2SP_AB] = std::isnan(R2SP_AB) || std::isinf(R2SP_AB) ? 0. : R2SP_AB;
    values[kWR2SP_AB] = std::isnan(R2SP_AB) || std::isinf(R2SP_AB) ? 0. : 1.0;
    values[kR2SP_AC] = std::isnan(R2SP_AC) || std::isinf(R2SP_AC) ? 0. : R2SP_AC;
    values[kWR2SP_AC] = std::isnan(R2SP_AC) || std::isinf(R2SP_AC) ? 0. : 1.0;
    values[kR2SP_BC] = std::isnan(R2SP_BC) || std::isinf(R2SP_BC) ? 0. : R2SP_BC;
    values[kWR2SP_BC] = std::isnan(R2SP_BC) || std::isinf(R2SP_BC) ? 0. : 1.0;
    values[kR2SP_AB_Im] = std::isnan(R2SP_AB_Im) || std::isinf(R2SP_AB_Im) ? 0. : R2SP_AB_Im;
    values[kWR2SP_AB_Im] = std::isnan(R2SP_AB_Im) || std::isinf(R2SP_AB_Im) ? 0. : 1.0;
    values[kR2SP_AC_Im] = std::isnan(R2SP_AC_Im) || std::isinf(R2SP_AC_Im) ? 0. : R2SP_AC_Im;
    values[kWR2SP_AC_Im] = std::isnan(R2SP_AC_Im) || std::isinf(R2SP_AC_Im) ? 0. : 1.0;
    values[kR2SP_BC_Im] = std::isnan(R2SP_BC_Im) || std::isinf(R2SP_BC_Im) ? 0. : R2SP_BC_Im;
    values[kWR2SP_BC_Im] = std::isnan(R2SP_BC_Im) || std::isinf(R2SP_BC_Im) ? 0. : 1.0;

    float R2EP_AB = std::isnan(Psi2A) || std::isinf(Psi2A) || std::isnan(Psi2B) || std::isinf(Psi2B) ? 0. : TMath::Cos(2 * (Psi2A - Psi2B));
    float R2EP_AC = std::isnan(Psi2A) || std::isinf(Psi2A) || std::isnan(Psi2C) || std::isinf(Psi2C) ? 0. : TMath::Cos(2 * (Psi2A - Psi2C));
    float R2EP_BC = std::isnan(Psi2B) || std::isinf(Psi2B) || std::isnan(Psi2C) || std::isinf(Psi2C) ? 0. : TMath::Cos(2 * (Psi2B - Psi2C));
    float R2EP_AB_Im = std::isnan(Psi2A) || std::isinf(Psi2A) || std::isnan(Psi2B) || std::isinf(Psi2B) ? 0. : TMath::Sin(2 * (Psi2A - Psi2B));
    float R2EP_AC_Im = std::isnan(Psi2A) || std::isinf(Psi2A) || std::isnan(Psi2C) || std::isinf(Psi2C) ? 0. : TMath::Sin(2 * (Psi2A - Psi2C));
    float R2EP_BC_Im = std::isnan(Psi2B) || std::isinf(Psi2B) || std::isnan(Psi2C) || std::isinf(Psi2C) ? 0. : TMath::Sin(2 * (Psi2B - Psi2C));
    values[kR2EP_AB] = std::isnan(R2EP_AB) || std::isinf(R2EP_AB) ? 0. : R2EP_AB;
    values[kWR2EP_AB] = std::isnan(R2EP_AB) || std::isinf(R2EP_AB) ? 0. : 1.0;
    values[kR2EP_AC] = std::isnan(R2EP_AC) || std::isinf(R2EP_AC) ? 0. : R2EP_AC;
    values[kWR2EP_AC] = std::isnan(R2EP_AC) || std::isinf(R2EP_AC) ? 0. : 1.0;
    values[kR2EP_BC] = std::isnan(R2EP_BC) || std::isinf(R2EP_BC) ? 0. : R2EP_BC;
    values[kWR2EP_BC] = std::isnan(R2EP_BC) || std::isinf(R2EP_BC) ? 0. : 1.0;
    values[kR2EP_AB_Im] = std::isnan(R2EP_AB_Im) || std::isinf(R2EP_AB_Im) ? 0. : R2EP_AB_Im;
    values[kWR2EP_AB_Im] = std::isnan(R2EP_AB_Im) || std::isinf(R2EP_AB_Im) ? 0. : 1.0;
    values[kR2EP_AC_Im] = std::isnan(R2EP_AC_Im) || std::isinf(R2EP_AC_Im) ? 0. : R2EP_AC_Im;
    values[kWR2EP_AC_Im] = std::isnan(R2EP_AC_Im) || std::isinf(R2EP_AC_Im) ? 0. : 1.0;
    values[kR2EP_BC_Im] = std::isnan(R2EP_BC_Im) || std::isinf(R2EP_BC_Im) ? 0. : R2EP_BC_Im;
    values[kWR2EP_BC_Im] = std::isnan(R2EP_BC_Im) || std::isinf(R2EP_BC_Im) ? 0. : 1.0;
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
    if constexpr ((fillMap & CollisionCent) > 0) {
      // WARNING: temporary solution, ongoing work to provide proper MC gen. centrality
      values[kMCEventCentrFT0C] = event.bestCollisionCentFT0C();
      values[kMultMCNParticlesEta05] = event.multMCNParticlesEta05();
      values[kMultMCNParticlesEta08] = event.multMCNParticlesEta08();
      values[kMultMCNParticlesEta10] = event.multMCNParticlesEta10();
    }
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
    values[kMCEventCentrFT0C] = event.centFT0C();
  }

  if constexpr ((fillMap & EventFilter) > 0 || (fillMap & RapidityGapFilter) > 0) {
    values[kIsDoubleGap] = (event.eventFilter() & (static_cast<uint64_t>(1) << kDoubleGap)) > 0;
    values[kIsSingleGapA] = (event.eventFilter() & (static_cast<uint64_t>(1) << kSingleGapA)) > 0;
    values[kIsSingleGapC] = (event.eventFilter() & (static_cast<uint64_t>(1) << kSingleGapC)) > 0;
    values[kIsSingleGap] = values[kIsSingleGapA] || values[kIsSingleGapC];
    values[kIsITSUPCMode] = (event.eventFilter() & (static_cast<uint64_t>(1) << kITSUPCMode)) > 0;
  }

  if constexpr ((fillMap & ReducedZdc) > 0) {
    FillZDC(event, values);
  }

  if constexpr ((fillMap & ReducedFit) > 0) {
    values[kAmplitudeFT0A] = event.amplitudeFT0A();
    values[kAmplitudeFT0C] = event.amplitudeFT0C();
    values[kTimeFT0A] = event.timeFT0A();
    values[kTimeFT0C] = event.timeFT0C();
    values[kTriggerMaskFT0] = event.triggerMaskFT0();
    values[kNFiredChannelsFT0A] = event.nFiredChannelsFT0A();
    values[kNFiredChannelsFT0C] = event.nFiredChannelsFT0C();
    values[kAmplitudeFDDA] = event.amplitudeFDDA();
    values[kAmplitudeFDDC] = event.amplitudeFDDC();
    values[kTimeFDDA] = event.timeFDDA();
    values[kTimeFDDC] = event.timeFDDC();
    values[kTriggerMaskFDD] = event.triggerMaskFDD();
    values[kAmplitudeFV0A] = event.amplitudeFV0A();
    values[kTimeFV0A] = event.timeFV0A();
    values[kTriggerMaskFV0A] = event.triggerMaskFV0A();
    values[kNFiredChannelsFV0A] = event.nFiredChannelsFV0A();
    values[kBBFT0Apf] = event.bbFT0Apf();
    values[kBGFT0Apf] = event.bgFT0Apf();
    values[kBBFT0Cpf] = event.bbFT0Cpf();
    values[kBGFT0Cpf] = event.bgFT0Cpf();
    values[kBBFV0Apf] = event.bbFV0Apf();
    values[kBGFV0Apf] = event.bgFV0Apf();
    values[kBBFDDApf] = event.bbFDDApf();
    values[kBGFDDApf] = event.bgFDDApf();
    values[kBBFDDCpf] = event.bbFDDCpf();
    values[kBGFDDCpf] = event.bgFDDCpf();
  }

  // FillEventDerived(values);
}

template <typename T>
void VarManager::FillEventFlowResoFactor(T const& hs_sp, T const& hs_ep, float* values)
{
  if (!values) {
    values = fgValues;
  }

  if (values[kCentFT0C] >= 0.) {
    int idx_sp = hs_sp->FindBin(values[kCentFT0C]);
    int idx_ep = hs_ep->FindBin(values[kCentFT0C]);

    values[kR2SP] = hs_sp->GetBinContent(idx_sp);
    values[kR2EP] = hs_ep->GetBinContent(idx_ep);
  }
}

template <typename T>
void VarManager::FillTwoMixEventsFlowResoFactor(T const& hs_sp, T const& hs_ep, float* values)
{
  if (!values) {
    values = fgValues;
  }

  if (values[kTwoEvCentFT0C1] >= 0.) {
    int idx_sp1 = hs_sp->FindBin(values[kTwoEvCentFT0C1]);
    int idx_ep1 = hs_ep->FindBin(values[kTwoEvCentFT0C1]);

    values[kTwoR2SP1] = hs_sp->GetBinContent(idx_sp1);
    values[kTwoR2EP1] = hs_ep->GetBinContent(idx_ep1);
  }

  if (values[kTwoEvCentFT0C2] >= 0.) {
    int idx_sp2 = hs_sp->FindBin(values[kTwoEvCentFT0C2]);
    int idx_ep2 = hs_ep->FindBin(values[kTwoEvCentFT0C2]);

    values[kTwoR2SP2] = hs_sp->GetBinContent(idx_sp2);
    values[kTwoR2EP2] = hs_ep->GetBinContent(idx_ep2);
  }
}

template <typename T, typename T1, typename T2>
void VarManager::FillTwoMixEventsCumulants(T const& h_v22ev1, T const& h_v24ev1, T const& h_v22ev2, T const& h_v24ev2, T1 const& t1, T2 const& t2, float* values)
{
  if (!values) {
    values = fgValues;
  }

  int idx_v22ev1;
  int idx_v24ev1;
  int idx_v22ev2;
  int idx_v24ev2;

  if (values[kTwoEvCentFT0C1] >= 0.) {
    if (t1.sign() < 0) {

      idx_v22ev1 = h_v22ev1->FindBin(values[kTwoEvCentFT0C1], t1.pt());
      idx_v24ev1 = h_v24ev1->FindBin(values[kTwoEvCentFT0C1], t1.pt());
      values[kV22m] = h_v22ev1->GetBinContent(idx_v22ev1);
      values[kV24m] = h_v24ev1->GetBinContent(idx_v24ev1);

    } else {

      idx_v22ev1 = h_v22ev2->FindBin(values[kTwoEvCentFT0C1], t1.pt());
      idx_v24ev1 = h_v24ev2->FindBin(values[kTwoEvCentFT0C1], t1.pt());
      values[kV22m] = h_v22ev2->GetBinContent(idx_v22ev1);
      values[kV24m] = h_v24ev2->GetBinContent(idx_v24ev1);
    }
  }
  if (values[kTwoEvCentFT0C2] >= 0.) {
    if (t2.sign() < 0) {

      idx_v22ev2 = h_v22ev1->FindBin(values[kTwoEvCentFT0C2], t2.pt());
      idx_v24ev2 = h_v24ev1->FindBin(values[kTwoEvCentFT0C2], t2.pt());
      values[kV22p] = h_v22ev1->GetBinContent(idx_v22ev2);
      values[kV24p] = h_v24ev1->GetBinContent(idx_v24ev2);

    } else {

      idx_v22ev2 = h_v22ev2->FindBin(values[kTwoEvCentFT0C2], t2.pt());
      idx_v24ev2 = h_v24ev2->FindBin(values[kTwoEvCentFT0C2], t2.pt());
      values[kV22p] = h_v22ev2->GetBinContent(idx_v22ev2);
      values[kV24p] = h_v24ev2->GetBinContent(idx_v24ev2);
    }
  }
}

template <typename T>
void VarManager::FillTwoEvents(T const& ev1, T const& ev2, float* values)
{
  if (!values) {
    values = fgValues;
  }
  // if constexpr (T::template contains<o2::aod::Collision>()) {
  values[kTwoEvPosZ1] = ev1.posZ();
  values[kTwoEvPosZ2] = ev2.posZ();
  values[kTwoEvPosR1] = std::sqrt(ev1.posX() * ev1.posX() + ev1.posY() * ev1.posY());
  values[kTwoEvPosR2] = std::sqrt(ev2.posX() * ev2.posX() + ev2.posY() * ev2.posY());
  values[kTwoEvDeltaZ] = ev1.posZ() - ev2.posZ();
  values[kTwoEvDeltaX] = ev1.posX() - ev2.posX();
  values[kTwoEvDeltaY] = ev1.posY() - ev2.posY();
  //}
  values[kTwoEvPVcontrib1] = ev1.numContrib();
  values[kTwoEvPVcontrib2] = ev2.numContrib();
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
  if constexpr ((fillMap & CollisionCent) > 0 || (fillMap & ReducedEventExtended) > 0) {
    values[kTwoEvCentFT0C1] = ev1.centFT0C();
    values[kTwoEvCentFT0C2] = ev2.centFT0C();
  }
  if constexpr ((fillMap & ReducedEventQvector) > 0) {
    // Tobe used for the calculation of u1q1 and u2q2
    values[kQ2X0A1] = ev1.q2x0a();
    values[kQ2X0A2] = ev2.q2x0a();
    values[kQ2Y0A1] = ev1.q2y0a();
    values[kQ2Y0A2] = ev2.q2y0a();
  }
  if constexpr ((fillMap & CollisionQvect) > 0) {
    // Tobe used for the calculation of u1q1 and u2q2
    values[kQ2X0A1] = (ev1.qvecBPosRe() * ev1.nTrkBPos() + ev1.qvecBNegRe() * ev1.nTrkBNeg()) / (ev1.nTrkBPos() + ev1.nTrkBNeg());
    values[kQ2X0A2] = (ev2.qvecBPosRe() * ev2.nTrkBPos() + ev2.qvecBNegRe() * ev2.nTrkBNeg()) / (ev2.nTrkBPos() + ev2.nTrkBNeg());
    values[kQ2Y0A1] = (ev1.qvecBPosIm() * ev1.nTrkBPos() + ev1.qvecBNegIm() * ev1.nTrkBNeg()) / (ev1.nTrkBPos() + ev1.nTrkBNeg());
    values[kQ2Y0A2] = (ev2.qvecBPosIm() * ev2.nTrkBPos() + ev2.qvecBNegIm() * ev2.nTrkBNeg()) / (ev2.nTrkBPos() + ev2.nTrkBNeg());
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
  if constexpr ((fillMap & Track) > 0 || (fillMap & Muon) > 0 || (fillMap & MuonRealign) > 0 || (fillMap & ReducedTrack) > 0 || (fillMap & ReducedMuon) > 0) {
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

    if (fgUsedVars[kM11REFoverMpsingle]) {
      float m = o2::constants::physics::MassMuon;
      ROOT::Math::PtEtaPhiMVector v(track.pt(), track.eta(), track.phi(), m);
      complex<double> Q21(values[kQ2X0A] * values[kS11A], values[kQ2Y0A] * values[kS11A]);
      complex<double> Q42(values[kQ42XA], values[kQ42YA]);
      complex<double> Q23(values[kQ23XA], values[kQ23YA]);
      complex<double> P2(TMath::Cos(2 * v.Phi()), TMath::Sin(2 * v.Phi()));
      values[kM11REFoverMpsingle] = values[kMultSingleMuons] > 0 && !(std::isnan(values[kM11REF]) || std::isinf(values[kM11REF]) || std::isnan(values[kCORR2REF]) || std::isinf(values[kCORR2REF]) || std::isnan(values[kM1111REF]) || std::isinf(values[kM1111REF]) || std::isnan(values[kCORR4REF]) || std::isinf(values[kCORR4REF])) ? values[kM11REF] / values[kMultSingleMuons] : 0;
      values[kM1111REFoverMpsingle] = values[kMultSingleMuons] > 0 && !(std::isnan(values[kM1111REF]) || std::isinf(values[kM1111REF]) || std::isnan(values[kCORR4REF]) || std::isinf(values[kCORR4REF]) || std::isnan(values[kM11REF]) || std::isinf(values[kM11REF]) || std::isnan(values[kCORR2REF]) || std::isinf(values[kCORR2REF])) ? values[kM1111REF] / values[kMultSingleMuons] : 0;
      values[kCORR2REFbysinglemu] = std::isnan(values[kM11REFoverMpsingle]) || std::isinf(values[kM11REFoverMpsingle]) || std::isnan(values[kCORR2REF]) || std::isinf(values[kCORR2REF]) || std::isnan(values[kM1111REFoverMpsingle]) || std::isinf(values[kM1111REFoverMpsingle]) || std::isnan(values[kCORR4REF]) || std::isinf(values[kCORR4REF]) ? 0 : values[kCORR2REF];
      values[kCORR4REFbysinglemu] = std::isnan(values[kM1111REFoverMpsingle]) || std::isinf(values[kM1111REFoverMpsingle]) || std::isnan(values[kCORR4REF]) || std::isinf(values[kCORR4REF]) || std::isnan(values[kM11REFoverMpsingle]) || std::isinf(values[kM11REFoverMpsingle]) || std::isnan(values[kCORR2REF]) || std::isinf(values[kCORR2REF]) ? 0 : values[kCORR4REF];
      values[kCORR2POIsingle] = (P2 * conj(Q21)).real() / values[kM01POI];
      values[kM01POIsingle] = values[kMultSingleMuons] * values[kS11A];
      values[kM0111POIsingle] = values[kMultSingleMuons] * (values[kS31A] - 3. * values[kS11A] * values[kS12A] + 2. * values[kS13A]);
      values[kCORR2POIsingle] = (P2 * conj(Q21)).real() / values[kM01POIsingle];
      values[kCORR4POIsingle] = (P2 * Q21 * conj(Q21) * conj(Q21) - P2 * Q21 * conj(Q42) - 2. * values[kS12A] * P2 * conj(Q21) + 2. * P2 * conj(Q23)).real() / values[kM0111POIsingle];
      values[kM01POIsingle] = std::isnan(values[kM01POIsingle]) || std::isinf(values[kM01POIsingle]) || std::isnan(values[kM0111POIsingle]) || std::isinf(values[kM0111POIsingle]) || std::isnan(values[kCORR2POIsingle]) || std::isinf(values[kCORR2POIsingle]) || std::isnan(values[kCORR4POIsingle]) || std::isinf(values[kCORR4POIsingle]) ? 0 : values[kM01POIsingle];
      values[kM0111POIsingle] = std::isnan(values[kM0111POIsingle]) || std::isinf(values[kM0111POIsingle]) || std::isnan(values[kCORR2POIsingle]) || std::isinf(values[kCORR2POIsingle]) || std::isnan(values[kCORR4POIsingle]) || std::isinf(values[kCORR4POIsingle]) ? 0 : values[kM0111POIsingle];
      values[kCORR2POIsingle] = std::isnan(values[kM01POIsingle]) || std::isinf(values[kM01POIsingle]) || std::isnan(values[kM0111POIsingle]) || std::isinf(values[kM0111POIsingle]) || std::isnan(values[kCORR2POIsingle]) || std::isinf(values[kCORR2POIsingle]) || std::isnan(values[kCORR4POIsingle]) || std::isinf(values[kCORR4POIsingle]) ? 0 : values[kCORR2POIsingle];
      values[kCORR4POIsingle] = std::isnan(values[kM01POIsingle]) || std::isinf(values[kM01POIsingle]) || std::isnan(values[kM0111POIsingle]) || std::isinf(values[kM0111POIsingle]) || std::isnan(values[kCORR2POIsingle]) || std::isinf(values[kCORR2POIsingle]) || std::isnan(values[kCORR4POIsingle]) || std::isinf(values[kCORR4POIsingle]) ? 0 : values[kCORR4POIsingle];
      values[kM01POIoverMpsingle] = values[kMultSingleMuons] > 0 && !(std::isnan(values[kM0111POIsingle]) || std::isinf(values[kM0111POIsingle]) || std::isnan(values[kCORR4POIsingle]) || std::isinf(values[kCORR4POIsingle]) || std::isnan(values[kM01POIsingle]) || std::isinf(values[kM01POIsingle]) || std::isnan(values[kCORR2POIsingle]) || std::isinf(values[kCORR2POIsingle])) ? values[kM01POIsingle] / values[kMultSingleMuons] : 0;
      values[kM0111POIoverMpsingle] = values[kMultSingleMuons] > 0 && !(std::isnan(values[kM0111POIsingle]) || std::isinf(values[kM0111POIsingle]) || std::isnan(values[kCORR4POIsingle]) || std::isinf(values[kCORR4POIsingle]) || std::isnan(values[kM01POIsingle]) || std::isinf(values[kM01POIsingle]) || std::isnan(values[kCORR2POIsingle]) || std::isinf(values[kCORR2POIsingle])) ? values[kM0111POIsingle] / values[kMultSingleMuons] : 0;
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
      if (fgUsedVars[kTPCnCRoverFindCls]) {
        values[kTPCnCRoverFindCls] = track.tpcCrossedRowsOverFindableCls();
      }
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
  if constexpr ((fillMap & TrackPID) > 0 || (fillMap & TrackTPCPID) > 0 || (fillMap & ReducedTrackBarrelPID) > 0) {
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
      if (!isTPCCalibrated) {
        values[kTPCnSigmaEl_Corr] = ComputePIDcalibration(0, values[kTPCnSigmaEl]);
      } else {
        LOG(fatal) << "TPC PID postcalibration is configured but the tracks are already postcalibrated. This is not allowed. Please check your configuration.";
        values[kTPCnSigmaEl_Corr] = track.tpcNSigmaEl();
      }
    }

    // compute TPC postcalibrated pion nsigma if required
    if (fgUsedVars[kTPCnSigmaPi_Corr] && fgRunTPCPostCalibration[1]) {
      if (!isTPCCalibrated) {
        values[kTPCnSigmaPi_Corr] = ComputePIDcalibration(1, values[kTPCnSigmaPi]);
      } else {
        LOG(fatal) << "TPC PID postcalibration is configured but the tracks are already postcalibrated. This is not allowed. Please check your configuration.";
        values[kTPCnSigmaPi_Corr] = track.tpcNSigmaPi();
      }
    }
    if (fgUsedVars[kTPCnSigmaKa_Corr] && fgRunTPCPostCalibration[2]) {
      // compute TPC postcalibrated kaon nsigma if required
      if (!isTPCCalibrated) {
        values[kTPCnSigmaKa_Corr] = ComputePIDcalibration(2, values[kTPCnSigmaKa]);
      } else {
        LOG(fatal) << "TPC PID postcalibration is configured but the tracks are already postcalibrated. This is not allowed. Please check your configuration.";
        values[kTPCnSigmaKa_Corr] = track.tpcNSigmaKa();
      }
    }
    // compute TPC postcalibrated proton nsigma if required
    if (fgUsedVars[kTPCnSigmaPr_Corr] && fgRunTPCPostCalibration[3]) {
      if (!isTPCCalibrated) {
        values[kTPCnSigmaPr_Corr] = ComputePIDcalibration(3, values[kTPCnSigmaPr]);
      } else {
        LOG(fatal) << "TPC PID postcalibration is configured but the tracks are already postcalibrated. This is not allowed. Please check your configuration.";
        values[kTPCnSigmaPr_Corr] = track.tpcNSigmaPr();
      }
    }

    if constexpr ((fillMap & TrackPID) > 0 || (fillMap & ReducedTrackBarrelPID) > 0) {
      values[kTOFnSigmaEl] = track.tofNSigmaEl();
      values[kTOFnSigmaPi] = track.tofNSigmaPi();
      values[kTOFnSigmaKa] = track.tofNSigmaKa();
      values[kTOFnSigmaPr] = track.tofNSigmaPr();
    }

    if constexpr ((fillMap & ReducedTrackBarrelPID) > 0) {
      values[kTPCnSigmaMu] = track.tpcNSigmaMu();
      values[kTOFnSigmaMu] = track.tofNSigmaMu();
      values[kTPCsignal] = track.tpcSignal();
      values[kTRDsignal] = track.trdSignal();
      values[kTOFbeta] = track.beta();
    }
  }
  if constexpr ((fillMap & TrackTOFService) > 0) {
    values[kTOFnSigmaEl] = track.tofNSigmaDynEl();
    values[kTOFnSigmaEl] = track.tofNSigmaDynPi();
    values[kTOFnSigmaEl] = track.tofNSigmaDynKa();
    values[kTOFnSigmaEl] = track.tofNSigmaDynPr();
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
  if constexpr ((fillMap & ReducedMuonExtra) > 0 || (fillMap & Muon) > 0 || (fillMap & MuonRealign) > 0) {
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
  if constexpr ((fillMap & ReducedMuonCov) > 0 || (fillMap & MuonCov) > 0 || (fillMap & MuonCovRealign) > 0) {
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
  if constexpr ((fillMap & MuonCov) > 0 || (fillMap & MuonCovRealign) > 0 || (fillMap & ReducedMuonCov) > 0) {

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

template <int candidateType, uint32_t fillMap, typename T1, typename T2, typename C>
void VarManager::FillTrackCollisionMC(T1 const& track, T2 const& MotherTrack, C const& collision, float* values)
{

  if (!values) {
    values = fgValues;
  }

  float m = o2::constants::physics::MassBPlus;
  float pdgLifetime = 1.638e-12; // s
  if (std::abs(MotherTrack.pdgCode()) == 511) {
    m = o2::constants::physics::MassB0;
    pdgLifetime = 1.517e-12; // s
  }

  // Extract the collision primary vertex position using constexpr, since the collision type may be CollisionMC or ReducedMCEvent
  double collPos[3] = {0.0, 0.0, 0.0};
  if constexpr (fillMap & ObjTypes::CollisionMC) {
    collPos[0] = collision.posX();
    collPos[1] = collision.posY();
    collPos[2] = collision.posZ();
  } else if constexpr (fillMap & ObjTypes::ReducedEventMC) {
    collPos[0] = collision.mcPosX();
    collPos[1] = collision.mcPosY();
    collPos[2] = collision.mcPosZ();
  }

  // displaced vertex is compued with decay product (track) and momentum of mother particle (MotherTrack)
  values[kMCVertexingLxy] = (collPos[0] - track.vx()) * (collPos[0] - track.vx()) +
                            (collPos[1] - track.vy()) * (collPos[1] - track.vy());
  values[kMCVertexingLz] = (collPos[2] - track.vz()) * (collPos[2] - track.vz());
  values[kMCVertexingLxyz] = values[kMCVertexingLxy] + values[kMCVertexingLz];
  values[kMCVertexingLxy] = std::sqrt(values[kMCVertexingLxy]);
  values[kMCVertexingLz] = std::sqrt(values[kMCVertexingLz]);
  values[kMCVertexingLxyz] = std::sqrt(values[kMCVertexingLxyz]);
  values[kMCVertexingTauz] = (collPos[2] - track.vz()) * m / (TMath::Abs(MotherTrack.pz()) * o2::constants::physics::LightSpeedCm2NS);
  values[kMCVertexingTauxy] = values[kMCVertexingLxy] * m / (MotherTrack.pt() * o2::constants::physics::LightSpeedCm2NS);

  values[kMCCosPointingAngle] = ((collPos[0] - track.vx()) * MotherTrack.px() +
                                 (collPos[1] - track.vy()) * MotherTrack.py() +
                                 (collPos[2] - track.vz()) * MotherTrack.pz()) /
                                (MotherTrack.p() * values[VarManager::kMCVertexingLxyz]);

  values[kMCLxyExpected] = (MotherTrack.pt() / m) * (pdgLifetime * o2::constants::physics::LightSpeedCm2S);
  values[kMCLxyzExpected] = (MotherTrack.p() / m) * (pdgLifetime * o2::constants::physics::LightSpeedCm2S);

  values[kMCVertexingLzProjected] = ((track.vz() - collPos[2]) * MotherTrack.pz()) / TMath::Abs(MotherTrack.pz());
  values[kMCVertexingLxyProjected] = (((track.vx() - collPos[0]) * MotherTrack.px()) + ((track.vy() - collPos[1]) * MotherTrack.py())) / TMath::Abs(MotherTrack.pt());
  values[kMCVertexingLxyzProjected] = (((track.vx() - collPos[1]) * MotherTrack.px()) + ((track.vy() - collPos[1]) * MotherTrack.py()) + ((track.vz() - collPos[2]) * MotherTrack.pz())) / MotherTrack.p();
  values[kMCVertexingTauxyProjected] = values[kMCVertexingLxyProjected] * m / (MotherTrack.pt());
  values[kMCVertexingTauzProjected] = values[kMCVertexingLzProjected] * m / TMath::Abs(MotherTrack.pz());
  values[kMCVertexingTauxyzProjected] = values[kMCVertexingLxyzProjected] * m / (MotherTrack.p());
}

template <int pairType, typename T, typename T1>
void VarManager::FillEnergyCorrelatorsMC(T const& track, T1 const& t1, float* values)
{
  // energy correlators
  float MassHadron;
  if constexpr (pairType == kJpsiHadronMass) {
    MassHadron = TMath::Sqrt(t1.e() * t1.e() - t1.p() * t1.p());
  }
  if constexpr (pairType == kJpsiPionMass) {
    MassHadron = o2::constants::physics::MassPionCharged;
  }
  ROOT::Math::PtEtaPhiMVector v1(track.pt(), track.eta(), track.phi(), o2::constants::physics::MassJPsi);
  float deltaphi = RecoDecay::constrainAngle(track.phi() - t1.phi(), -o2::constants::math::PIHalf);
  float deltaeta = t1.eta() - track.eta();
  ROOT::Math::PtEtaPhiMVector v2(t1.pt(), t1.eta(), t1.phi(), MassHadron);
  float E_boost = LorentzTransformJpsihadroncosChi("weight_boost", v1, v2);
  float CosChi = LorentzTransformJpsihadroncosChi("coschi", v1, v2);
  float CosTheta = LorentzTransformJpsihadroncosChi("costheta", v1, v2);
  values[kMCCosChi] = CosChi;
  values[kMCWeight_before] = t1.pt() / o2::constants::physics::MassJPsi;
  values[kMCCosTheta] = CosTheta;
  values[kMCdeltaphi] = deltaphi;
  values[kMCdeltaeta] = deltaeta;
  values[kMCHadronPt] = t1.pt();
  values[kMCHadronEta] = t1.eta();
  values[kMCHadronPhi] = RecoDecay::constrainAngle(t1.phi(), -o2::constants::math::PIHalf);
  values[kMCHadronPdgCode] = t1.pdgCode();
  values[kMCWeight] = E_boost / o2::constants::physics::MassJPsi;

  values[kMCCosChi_randomPhi_trans] = -999.9f;
  values[kMCCosChi_randomPhi_toward] = -999.9f;
  values[kMCCosChi_randomPhi_away] = -999.9f;

  values[kMCdeltaphi_randomPhi_trans] = -999.9f;
  values[kMCdeltaphi_randomPhi_toward] = -999.9f;
  values[kMCdeltaphi_randomPhi_away] = -999.9f;

  float randomPhi_trans = -o2::constants::math::PIHalf;
  float randomPhi_toward = -o2::constants::math::PIHalf;
  float randomPhi_away = -o2::constants::math::PIHalf;

  if ((deltaphi > -0.5 * TMath::Pi() && deltaphi < -1. / 3 * TMath::Pi()) || (deltaphi > 4. / 3 * TMath::Pi() && deltaphi < 1.5 * TMath::Pi()) || (deltaphi > 1. / 3 * TMath::Pi() && deltaphi < 2. / 3 * TMath::Pi())) {
    randomPhi_trans = gRandom->Uniform(-o2::constants::math::PIHalf, 3. * o2::constants::math::PIHalf);
    randomPhi_toward = gRandom->Uniform(-o2::constants::math::PIHalf, 3. * o2::constants::math::PIHalf);
    randomPhi_away = gRandom->Uniform(-o2::constants::math::PIHalf, 3. * o2::constants::math::PIHalf);

    ROOT::Math::PtEtaPhiMVector v2_randomPhi_trans(v2.pt(), v2.eta(), randomPhi_trans, o2::constants::physics::MassPionCharged);
    values[kMCCosChi_randomPhi_trans] = LorentzTransformJpsihadroncosChi("coschi", v1, v2_randomPhi_trans);
    values[kMCWeight_randomPhi_trans] = LorentzTransformJpsihadroncosChi("weight_boost", v1, v2_randomPhi_trans) / v1.M();

    ROOT::Math::PtEtaPhiMVector v2_randomPhi_toward(v2.pt(), v2.eta(), randomPhi_toward, o2::constants::physics::MassPionCharged);
    values[kMCCosChi_randomPhi_toward] = LorentzTransformJpsihadroncosChi("coschi", v1, v2_randomPhi_toward);
    values[kMCWeight_randomPhi_toward] = LorentzTransformJpsihadroncosChi("weight_boost", v1, v2_randomPhi_toward) / v1.M();

    ROOT::Math::PtEtaPhiMVector v2_randomPhi_away(v2.pt(), v2.eta(), randomPhi_away, o2::constants::physics::MassPionCharged);
    values[kMCCosChi_randomPhi_away] = LorentzTransformJpsihadroncosChi("coschi", v1, v2_randomPhi_away);
    values[kMCWeight_randomPhi_away] = LorentzTransformJpsihadroncosChi("weight_boost", v1, v2_randomPhi_away) / v1.M();

    values[kMCdeltaphi_randomPhi_trans] = RecoDecay::constrainAngle(v1.phi() - randomPhi_trans, -o2::constants::math::PIHalf);
    values[kMCdeltaphi_randomPhi_toward] = RecoDecay::constrainAngle(v1.phi() - randomPhi_toward, -o2::constants::math::PIHalf);
    values[kMCdeltaphi_randomPhi_away] = RecoDecay::constrainAngle(v1.phi() - randomPhi_away, -o2::constants::math::PIHalf);
  }
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
    // Make the TPC information of the kaon available for pair histograms
    values[kPin_leg1] = t1.tpcInnerParam();
    values[kTPCnSigmaKa_leg1] = t1.tpcNSigmaKa();
  }

  if constexpr (pairType == kElectronMuon) {
    m2 = o2::constants::physics::MassMuon;
  }

  values[kCharge] = t1.sign() + t2.sign();
  values[kCharge1] = t1.sign();
  values[kCharge2] = t2.sign();
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

  values[kPt1] = t1.pt();
  values[kEta1] = t1.eta();
  values[kPhi1] = t1.phi();
  values[kPt2] = t2.pt();
  values[kEta2] = t2.eta();
  values[kPhi2] = t2.phi();

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

  if (fgUsedVars[kDeltaEtaPair2]) {
    values[kDeltaEtaPair2] = v1.Eta() - v2.Eta();
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

  // polarization parameters
  bool useHE = fgUsedVars[kCosThetaHE] || fgUsedVars[kPhiHE]; // helicity frame
  bool useCS = fgUsedVars[kCosThetaCS] || fgUsedVars[kPhiCS]; // Collins-Soper frame
  bool usePP = fgUsedVars[kCosThetaPP];                       // production plane frame
  bool useRM = fgUsedVars[kCosThetaRM];                       // Random frame

  if (useHE || useCS || usePP || useRM) {
    ROOT::Math::Boost boostv12{v12.BoostToCM()};
    ROOT::Math::XYZVectorF v1_CM{(boostv12(v1).Vect()).Unit()};
    ROOT::Math::XYZVectorF v2_CM{(boostv12(v2).Vect()).Unit()};
    ROOT::Math::XYZVectorF Beam1_CM{(boostv12(fgBeamA).Vect()).Unit()};
    ROOT::Math::XYZVectorF Beam2_CM{(boostv12(fgBeamC).Vect()).Unit()};

    // using positive sign convention for the first track
    ROOT::Math::XYZVectorF v_CM = (t1.sign() > 0 ? v1_CM : v2_CM);

    if (useHE) {
      ROOT::Math::XYZVectorF zaxis_HE{(v12.Vect()).Unit()};
      ROOT::Math::XYZVectorF yaxis_HE{(Beam1_CM.Cross(Beam2_CM)).Unit()};
      ROOT::Math::XYZVectorF xaxis_HE{(yaxis_HE.Cross(zaxis_HE)).Unit()};
      if (fgUsedVars[kCosThetaHE])
        values[kCosThetaHE] = zaxis_HE.Dot(v_CM);
      if (fgUsedVars[kPhiHE]) {
        values[kPhiHE] = TMath::ATan2(yaxis_HE.Dot(v_CM), xaxis_HE.Dot(v_CM));
        if (values[kPhiHE] < 0) {
          values[kPhiHE] += 2 * TMath::Pi(); // ensure phi is in [0, 2pi]
        }
      }
      if (fgUsedVars[kPhiTildeHE]) {
        if (fgUsedVars[kCosThetaHE] && fgUsedVars[kPhiHE]) {
          if (values[kCosThetaHE] > 0) {
            values[kPhiTildeHE] = values[kPhiHE] - 0.25 * TMath::Pi(); // phi_tilde = phi - pi/4
            if (values[kPhiTildeHE] < 0) {
              values[kPhiTildeHE] += 2 * TMath::Pi(); // ensure phi_tilde is in [0, 2pi]
            }
          } else {
            values[kPhiTildeHE] = values[kPhiHE] - 0.75 * TMath::Pi(); // phi_tilde = phi - 3pi/4
            if (values[kPhiTildeHE] < 0) {
              values[kPhiTildeHE] += 2 * TMath::Pi(); // ensure phi_tilde is in [0, 2pi]
            }
          }
        } else {
          values[kPhiTildeHE] = -999; // not computable
        }
      }
    }

    if (useCS) {
      ROOT::Math::XYZVectorF zaxis_CS{(Beam1_CM - Beam2_CM).Unit()};
      ROOT::Math::XYZVectorF yaxis_CS{(Beam1_CM.Cross(Beam2_CM)).Unit()};
      ROOT::Math::XYZVectorF xaxis_CS{(yaxis_CS.Cross(zaxis_CS)).Unit()};
      if (fgUsedVars[kCosThetaCS])
        values[kCosThetaCS] = zaxis_CS.Dot(v_CM);
      if (fgUsedVars[kPhiCS]) {
        values[kPhiCS] = TMath::ATan2(yaxis_CS.Dot(v_CM), xaxis_CS.Dot(v_CM));
        if (values[kPhiCS] < 0) {
          values[kPhiCS] += 2 * TMath::Pi(); // ensure phi is in [0, 2pi]
        }
      }
      if (fgUsedVars[kPhiTildeCS]) {
        if (fgUsedVars[kCosThetaCS] && fgUsedVars[kPhiCS]) {
          if (values[kCosThetaCS] > 0) {
            values[kPhiTildeCS] = values[kPhiCS] - 0.25 * TMath::Pi(); // phi_tilde = phi - pi/4
            if (values[kPhiTildeCS] < 0) {
              values[kPhiTildeCS] += 2 * TMath::Pi(); // ensure phi_tilde is in [0, 2pi]
            }
          } else {
            values[kPhiTildeCS] = values[kPhiCS] - 0.75 * TMath::Pi(); // phi_tilde = phi - 3pi/4
            if (values[kPhiTildeCS] < 0) {
              values[kPhiTildeCS] += 2 * TMath::Pi(); // ensure phi_tilde is in [0, 2pi]
            }
          }
        } else {
          values[kPhiTildeCS] = -999; // not computable
        }
      }
    }

    if (usePP) {
      ROOT::Math::XYZVector zaxis_PP = ROOT::Math::XYZVector(v12.Py(), -v12.Px(), 0.f);
      ROOT::Math::XYZVector yaxis_PP{(v12.Vect()).Unit()};
      ROOT::Math::XYZVector xaxis_PP{(yaxis_PP.Cross(zaxis_PP)).Unit()};
      if (fgUsedVars[kCosThetaPP]) {
        values[kCosThetaPP] = zaxis_PP.Dot(v_CM) / std::sqrt(zaxis_PP.Mag2());
      }
      if (fgUsedVars[kPhiPP]) {
        values[kPhiPP] = TMath::ATan2(yaxis_PP.Dot(v_CM), xaxis_PP.Dot(v_CM));
        if (values[kPhiPP] < 0) {
          values[kPhiPP] += 2 * TMath::Pi(); // ensure phi is in [0, 2pi]
        }
      }
      if (fgUsedVars[kPhiTildePP]) {
        if (fgUsedVars[kCosThetaPP] && fgUsedVars[kPhiPP]) {
          if (values[kCosThetaPP] > 0) {
            values[kPhiTildePP] = values[kPhiPP] - 0.25 * TMath::Pi(); // phi_tilde = phi - pi/4
            if (values[kPhiTildePP] < 0) {
              values[kPhiTildePP] += 2 * TMath::Pi(); // ensure phi_tilde is in [0, 2pi]
            }
          } else {
            values[kPhiTildePP] = values[kPhiPP] - 0.75 * TMath::Pi(); // phi_tilde = phi - 3pi/4
            if (values[kPhiTildePP] < 0) {
              values[kPhiTildePP] += 2 * TMath::Pi(); // ensure phi_tilde is in [0, 2pi]
            }
          }
        } else {
          values[kPhiTildePP] = -999; // not computable
        }
      }
    }

    if (useRM) {
      double randomCostheta = gRandom->Uniform(-1., 1.);
      double randomPhi = gRandom->Uniform(0., 2. * TMath::Pi());
      ROOT::Math::XYZVectorF zaxis_RM(randomCostheta, std::sqrt(1 - randomCostheta * randomCostheta) * std::cos(randomPhi), std::sqrt(1 - randomCostheta * randomCostheta) * std::sin(randomPhi));
      if (fgUsedVars[kCosThetaRM])
        values[kCosThetaRM] = zaxis_RM.Dot(v_CM);
    }
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
    values[kCharge] = t1.sign() + t2.sign() + t3.sign();
    values[kS12] = (v1 + v2).M2();
    values[kS13] = (v1 + v3).M2();
    values[kS23] = (v2 + v3).M2();
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

template <uint32_t fillMap, int pairType, typename T1, typename T2>
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

  if (fgUsedVars[kDeltaPhiPair2]) {
    double phipair2ME = v1.Phi() - v2.Phi();
    if (phipair2ME > 3 * TMath::Pi() / 2) {
      values[kDeltaPhiPair2] = phipair2ME - 2 * TMath::Pi();
    } else if (phipair2ME < -TMath::Pi() / 2) {
      values[kDeltaPhiPair2] = phipair2ME + 2 * TMath::Pi();
    } else {
      values[kDeltaPhiPair2] = phipair2ME;
    }
  }

  if (fgUsedVars[kDeltaEtaPair2]) {
    values[kDeltaEtaPair2] = v1.Eta() - v2.Eta();
  }

  // polarization parameters
  bool useHE = fgUsedVars[kCosThetaHE] || fgUsedVars[kPhiHE]; // helicity frame
  bool useCS = fgUsedVars[kCosThetaCS] || fgUsedVars[kPhiCS]; // Collins-Soper frame
  bool usePP = fgUsedVars[kCosThetaPP];                       // production plane frame
  bool useRM = fgUsedVars[kCosThetaRM];                       // Random frame

  if (useHE || useCS || usePP || useRM) {
    ROOT::Math::Boost boostv12{v12.BoostToCM()};
    ROOT::Math::XYZVectorF v1_CM{(boostv12(v1).Vect()).Unit()};
    ROOT::Math::XYZVectorF v2_CM{(boostv12(v2).Vect()).Unit()};
    ROOT::Math::XYZVectorF Beam1_CM{(boostv12(fgBeamA).Vect()).Unit()};
    ROOT::Math::XYZVectorF Beam2_CM{(boostv12(fgBeamC).Vect()).Unit()};

    // using positive sign convention for the first track
    ROOT::Math::XYZVectorF v_CM = (t1.sign() > 0 ? v1_CM : v2_CM);

    if (useHE) {
      ROOT::Math::XYZVectorF zaxis_HE{(v12.Vect()).Unit()};
      ROOT::Math::XYZVectorF yaxis_HE{(Beam1_CM.Cross(Beam2_CM)).Unit()};
      ROOT::Math::XYZVectorF xaxis_HE{(yaxis_HE.Cross(zaxis_HE)).Unit()};
      if (fgUsedVars[kCosThetaHE])
        values[kCosThetaHE] = zaxis_HE.Dot(v_CM);
      if (fgUsedVars[kPhiHE]) {
        values[kPhiHE] = TMath::ATan2(yaxis_HE.Dot(v_CM), xaxis_HE.Dot(v_CM));
        if (values[kPhiHE] < 0) {
          values[kPhiHE] += 2 * TMath::Pi(); // ensure phi is in [0, 2pi]
        }
      }
      if (fgUsedVars[kPhiTildeHE]) {
        if (fgUsedVars[kCosThetaHE] && fgUsedVars[kPhiHE]) {
          if (values[kCosThetaHE] > 0) {
            values[kPhiTildeHE] = values[kPhiHE] - 0.25 * TMath::Pi(); // phi_tilde = phi - pi/4
            if (values[kPhiTildeHE] < 0) {
              values[kPhiTildeHE] += 2 * TMath::Pi(); // ensure phi_tilde is in [0, 2pi]
            }
          } else {
            values[kPhiTildeHE] = values[kPhiHE] - 0.75 * TMath::Pi(); // phi_tilde = phi - 3pi/4
            if (values[kPhiTildeHE] < 0) {
              values[kPhiTildeHE] += 2 * TMath::Pi(); // ensure phi_tilde is in [0, 2pi]
            }
          }
        } else {
          values[kPhiTildeHE] = -999; // not computable
        }
      }
    }

    if (useCS) {
      ROOT::Math::XYZVectorF zaxis_CS{(Beam1_CM - Beam2_CM).Unit()};
      ROOT::Math::XYZVectorF yaxis_CS{(Beam1_CM.Cross(Beam2_CM)).Unit()};
      ROOT::Math::XYZVectorF xaxis_CS{(yaxis_CS.Cross(zaxis_CS)).Unit()};
      if (fgUsedVars[kCosThetaCS])
        values[kCosThetaCS] = zaxis_CS.Dot(v_CM);
      if (fgUsedVars[kPhiCS]) {
        values[kPhiCS] = TMath::ATan2(yaxis_CS.Dot(v_CM), xaxis_CS.Dot(v_CM));
        if (values[kPhiCS] < 0) {
          values[kPhiCS] += 2 * TMath::Pi(); // ensure phi is in [0, 2pi]
        }
      }
      if (fgUsedVars[kPhiTildeCS]) {
        if (fgUsedVars[kCosThetaCS] && fgUsedVars[kPhiCS]) {
          if (values[kCosThetaCS] > 0) {
            values[kPhiTildeCS] = values[kPhiCS] - 0.25 * TMath::Pi(); // phi_tilde = phi - pi/4
            if (values[kPhiTildeCS] < 0) {
              values[kPhiTildeCS] += 2 * TMath::Pi(); // ensure phi_tilde is in [0, 2pi]
            }
          } else {
            values[kPhiTildeCS] = values[kPhiCS] - 0.75 * TMath::Pi(); // phi_tilde = phi - 3pi/4
            if (values[kPhiTildeCS] < 0) {
              values[kPhiTildeCS] += 2 * TMath::Pi(); // ensure phi_tilde is in [0, 2pi]
            }
          }
        } else {
          values[kPhiTildeCS] = -999; // not computable
        }
      }
    }

    if (usePP) {
      ROOT::Math::XYZVector zaxis_PP = ROOT::Math::XYZVector(v12.Py(), -v12.Px(), 0.f);
      ROOT::Math::XYZVector yaxis_PP{(v12.Vect()).Unit()};
      ROOT::Math::XYZVector xaxis_PP{(yaxis_PP.Cross(zaxis_PP)).Unit()};
      if (fgUsedVars[kCosThetaPP]) {
        values[kCosThetaPP] = zaxis_PP.Dot(v_CM) / std::sqrt(zaxis_PP.Mag2());
      }
      if (fgUsedVars[kPhiPP]) {
        values[kPhiPP] = TMath::ATan2(yaxis_PP.Dot(v_CM), xaxis_PP.Dot(v_CM));
        if (values[kPhiPP] < 0) {
          values[kPhiPP] += 2 * TMath::Pi(); // ensure phi is in [0, 2pi]
        }
      }
      if (fgUsedVars[kPhiTildePP]) {
        if (fgUsedVars[kCosThetaPP] && fgUsedVars[kPhiPP]) {
          if (values[kCosThetaPP] > 0) {
            values[kPhiTildePP] = values[kPhiPP] - 0.25 * TMath::Pi(); // phi_tilde = phi - pi/4
            if (values[kPhiTildePP] < 0) {
              values[kPhiTildePP] += 2 * TMath::Pi(); // ensure phi_tilde is in [0, 2pi]
            }
          } else {
            values[kPhiTildePP] = values[kPhiPP] - 0.75 * TMath::Pi(); // phi_tilde = phi - 3pi/4
            if (values[kPhiTildePP] < 0) {
              values[kPhiTildePP] += 2 * TMath::Pi(); // ensure phi_tilde is in [0, 2pi]
            }
          }
        } else {
          values[kPhiTildePP] = -999; // not computable
        }
      }
    }

    if (useRM) {
      double randomCostheta = gRandom->Uniform(-1., 1.);
      double randomPhi = gRandom->Uniform(0., 2. * TMath::Pi());
      ROOT::Math::XYZVectorF zaxis_RM(randomCostheta, std::sqrt(1 - randomCostheta * randomCostheta) * std::cos(randomPhi), std::sqrt(1 - randomCostheta * randomCostheta) * std::sin(randomPhi));
      if (fgUsedVars[kCosThetaRM])
        values[kCosThetaRM] = zaxis_RM.Dot(v_CM);
    }
  }

  if constexpr ((fillMap & ReducedEventQvector) > 0 || (fillMap & CollisionQvect) > 0) {
    // TODO: provide different computations for vn
    // Compute the scalar product UQ for two muon from different event using Q-vector from A, for second and third harmonic
    float Psi2A1 = getEventPlane(2, values[kQ2X0A1], values[kQ2Y0A1]);
    float Psi2A2 = getEventPlane(2, values[kQ2X0A2], values[kQ2Y0A2]);
    values[kCos2DeltaPhi] = TMath::Cos(2 * (v12.Phi() - Psi2A1)); // WARNING: using the first event EP
    values[kCos2DeltaPhiEv1] = TMath::Cos(2 * (v1.Phi() - Psi2A1));
    values[kCos2DeltaPhiEv2] = TMath::Cos(2 * (v2.Phi() - Psi2A2));
    values[kU2Q2] = values[kQ2X0A1] * TMath::Cos(2 * v12.Phi()) + values[kQ2Y0A1] * TMath::Sin(2 * v12.Phi()); // WARNING: using the first event EP
    values[kU2Q2Ev1] = values[kQ2X0A1] * TMath::Cos(2 * v1.Phi()) + values[kQ2Y0A1] * TMath::Sin(2 * v1.Phi());
    values[kU2Q2Ev2] = values[kQ2X0A2] * TMath::Cos(2 * v2.Phi()) + values[kQ2Y0A2] * TMath::Sin(2 * v2.Phi());

    values[kCos2DeltaPhiMu1] = TMath::Cos(2 * (v1.Phi() - v12.Phi()));
    values[kCos2DeltaPhiMu2] = TMath::Cos(2 * (v2.Phi() - v12.Phi()));

    values[kV2SP1] = values[kU2Q2Ev1] / values[kTwoR2SP1];
    values[kV2SP2] = values[kU2Q2Ev2] / values[kTwoR2SP2];
    values[kV2EP1] = values[kCos2DeltaPhiEv1] / values[kTwoR2EP1];
    values[kV2EP2] = values[kCos2DeltaPhiEv2] / values[kTwoR2EP2];

    float V2ME_SP = values[kV2SP1] * values[kCos2DeltaPhiMu1] + values[kV2SP2] * values[kCos2DeltaPhiMu2];
    float V2ME_EP = values[kV2EP1] * values[kCos2DeltaPhiMu1] + values[kV2EP2] * values[kCos2DeltaPhiMu2];
    values[kV2ME_SP] = std::isnan(V2ME_SP) || std::isinf(V2ME_SP) ? 0. : V2ME_SP;
    values[kWV2ME_SP] = std::isnan(V2ME_SP) || std::isinf(V2ME_SP) ? 0. : 1.0;
    values[kV2ME_EP] = std::isnan(V2ME_EP) || std::isinf(V2ME_EP) ? 0. : V2ME_EP;
    values[kWV2ME_EP] = std::isnan(V2ME_EP) || std::isinf(V2ME_EP) ? 0. : 1.0;

    // Cumulant part
    float V22ME = values[kV22m] * values[kCos2DeltaPhiMu1] + values[kV22p] * values[kCos2DeltaPhiMu2];
    float V24ME = values[kV24m] * values[kCos2DeltaPhiMu1] + values[kV24p] * values[kCos2DeltaPhiMu2];
    values[kV22ME] = (std::isnan(V22ME) || std::isinf(V22ME) || std::isnan(V24ME) || std::isinf(V24ME)) ? 0. : V22ME;
    values[kWV22ME] = (std::isnan(V22ME) || std::isinf(V22ME) || std::isnan(V24ME) || std::isinf(V24ME)) ? 0. : 1.0;
    values[kV24ME] = (std::isnan(V22ME) || std::isinf(V22ME) || std::isnan(V24ME) || std::isinf(V24ME)) ? 0. : V24ME;
    values[kWV24ME] = (std::isnan(V22ME) || std::isinf(V22ME) || std::isnan(V24ME) || std::isinf(V24ME)) ? 0. : 1.0;

    if constexpr ((fillMap & ReducedEventQvectorExtra) > 0) {
      complex<double> Q21(values[kQ2X0A] * values[kS11A], values[kQ2Y0A] * values[kS11A]);
      complex<double> Q42(values[kQ42XA], values[kQ42YA]);
      complex<double> Q23(values[kQ23XA], values[kQ23YA]);
      complex<double> P2(TMath::Cos(2 * v12.Phi()), TMath::Sin(2 * v12.Phi()));
      values[kM01POIME] = values[kMultDimuonsME] * values[kS11A];
      values[kM0111POIME] = values[kMultDimuonsME] * (values[kS31A] - 3. * values[kS11A] * values[kS12A] + 2. * values[kS13A]);
      values[kCORR2POIME] = (P2 * conj(Q21)).real() / values[kM01POIME];
      values[kCORR4POIME] = (P2 * Q21 * conj(Q21) * conj(Q21) - P2 * Q21 * conj(Q42) - 2. * values[kS12A] * P2 * conj(Q21) + 2. * P2 * conj(Q23)).real() / values[kM0111POIME];
      values[kM01POIoverMpME] = values[kMultDimuonsME] > 0 && !(std::isnan(values[kM01POIME]) || std::isinf(values[kM01POIME]) || std::isnan(values[kCORR2POIME]) || std::isinf(values[kCORR2POIME]) || std::isnan(values[kM0111POIME]) || std::isinf(values[kM0111POIME]) || std::isnan(values[kCORR4POIME]) || std::isinf(values[kCORR4POIME])) ? values[kM01POIME] / values[kMultDimuonsME] : 0;
      values[kM0111POIoverMpME] = values[kMultDimuonsME] > 0 && !(std::isnan(values[kM0111POIME]) || std::isinf(values[kM0111POIME]) || std::isnan(values[kCORR4POIME]) || std::isinf(values[kCORR4POIME]) || std::isnan(values[kM01POIME]) || std::isinf(values[kM01POIME]) || std::isnan(values[kCORR2POIME]) || std::isinf(values[kCORR2POIME])) ? values[kM0111POIME] / values[kMultDimuonsME] : 0;
      values[kM11REFoverMpME] = values[kMultDimuonsME] > 0 && !(std::isnan(values[kM11REF]) || std::isinf(values[kM11REF]) || std::isnan(values[kCORR2REF]) || std::isinf(values[kCORR2REF]) || std::isnan(values[kM1111REF]) || std::isinf(values[kM1111REF]) || std::isnan(values[kCORR4REF]) || std::isinf(values[kCORR4REF])) ? values[kM11REF] / values[kMultDimuonsME] : 0;
      values[kM1111REFoverMpME] = values[kMultDimuonsME] > 0 && !(std::isnan(values[kM1111REF]) || std::isinf(values[kM1111REF]) || std::isnan(values[kCORR4REF]) || std::isinf(values[kCORR4REF]) || std::isnan(values[kM11REF]) || std::isinf(values[kM11REF]) || std::isnan(values[kCORR2REF]) || std::isinf(values[kCORR2REF])) ? values[kM1111REF] / values[kMultDimuonsME] : 0;
      values[kCORR2REFbydimuonsME] = std::isnan(values[kM11REFoverMpME]) || std::isinf(values[kM11REFoverMpME]) || std::isnan(values[kCORR2REF]) || std::isinf(values[kCORR2REF]) || std::isnan(values[kM1111REFoverMpME]) || std::isinf(values[kM1111REFoverMpME]) || std::isnan(values[kCORR4REF]) || std::isinf(values[kCORR4REF]) ? 0 : values[kCORR2REF];
      values[kCORR4REFbydimuonsME] = std::isnan(values[kM1111REFoverMpME]) || std::isinf(values[kM1111REFoverMpME]) || std::isnan(values[kCORR4REF]) || std::isinf(values[kCORR4REF]) || std::isnan(values[kM11REFoverMpME]) || std::isinf(values[kM11REFoverMpME]) || std::isnan(values[kCORR2REF]) || std::isinf(values[kCORR2REF]) ? 0 : values[kCORR4REF];
      values[kCORR2POIME] = std::isnan(values[kCORR2POIME]) || std::isinf(values[kCORR2POIME]) || std::isnan(values[kM01POIME]) || std::isinf(values[kM01POIME]) || std::isnan(values[kCORR4POIME]) || std::isinf(values[kCORR4POIME]) || std::isnan(values[kM0111POIME]) || std::isinf(values[kM0111POIME]) ? 0 : values[kCORR2POIME];
      values[kCORR4POIME] = std::isnan(values[kCORR4POIME]) || std::isinf(values[kCORR4POIME]) || std::isnan(values[kM0111POIME]) || std::isinf(values[kM0111POIME]) || std::isnan(values[kCORR2POIME]) || std::isinf(values[kCORR2POIME]) || std::isnan(values[kM01POIME]) || std::isinf(values[kM01POIME]) ? 0 : values[kCORR4POIME];
    }
  }
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

template <int pairType, typename T1, typename T2>
void VarManager::FillPairMC(T1 const& t1, T2 const& t2, float* values)
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

  if (pairType == kDecayToKPi) {
    m1 = o2::constants::physics::MassKaonCharged;
    m2 = o2::constants::physics::MassPionCharged;
  }

  if (pairType == kElectronMuon) {
    m2 = o2::constants::physics::MassMuon;
  }

  // TODO : implement resolution smearing.
  ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), m1);
  ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), m2);
  ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
  values[kMCMass] = v12.M();
  values[kMCPt] = v12.Pt();
  values[kMCEta] = v12.Eta();
  values[kMCPhi] = v12.Phi();
  values[kMCY] = -v12.Rapidity();

  // polarization parameters
  bool useHE = fgUsedVars[kMCCosThetaHE] || fgUsedVars[kMCPhiHE]; // helicity frame
  bool useCS = fgUsedVars[kMCCosThetaCS] || fgUsedVars[kMCPhiCS]; // Collins-Soper frame
  bool usePP = fgUsedVars[kMCCosThetaPP];                         // production plane frame
  bool useRM = fgUsedVars[kMCCosThetaRM];                         // Random frame

  if (useHE || useCS || usePP || useRM) {
    ROOT::Math::Boost boostv12{v12.BoostToCM()};
    ROOT::Math::XYZVectorF v1_CM{(boostv12(v1).Vect()).Unit()};
    ROOT::Math::XYZVectorF v2_CM{(boostv12(v2).Vect()).Unit()};
    ROOT::Math::XYZVectorF Beam1_CM{(boostv12(fgBeamA).Vect()).Unit()};
    ROOT::Math::XYZVectorF Beam2_CM{(boostv12(fgBeamC).Vect()).Unit()};

    // using positive sign convention for the first track
    ROOT::Math::XYZVectorF v_CM = (t1.pdgCode() > 0 ? v1_CM : v2_CM);

    if (useHE) {
      ROOT::Math::XYZVectorF zaxis_HE{(v12.Vect()).Unit()};
      ROOT::Math::XYZVectorF yaxis_HE{(Beam1_CM.Cross(Beam2_CM)).Unit()};
      ROOT::Math::XYZVectorF xaxis_HE{(yaxis_HE.Cross(zaxis_HE)).Unit()};
      if (fgUsedVars[kMCCosThetaHE])
        values[kMCCosThetaHE] = zaxis_HE.Dot(v_CM);
      if (fgUsedVars[kMCPhiHE]) {
        values[kMCPhiHE] = TMath::ATan2(yaxis_HE.Dot(v_CM), xaxis_HE.Dot(v_CM));
        if (values[kMCPhiHE] < 0) {
          values[kMCPhiHE] += 2 * TMath::Pi(); // ensure phi is in [0, 2pi]
        }
      }
      if (fgUsedVars[kMCPhiTildeHE]) {
        if (fgUsedVars[kMCCosThetaHE] && fgUsedVars[kMCPhiHE]) {
          if (values[kMCCosThetaHE] > 0) {
            values[kMCPhiTildeHE] = values[kMCPhiHE] - 0.25 * TMath::Pi(); // phi_tilde = phi - pi/4
            if (values[kMCPhiTildeHE] < 0) {
              values[kMCPhiTildeHE] += 2 * TMath::Pi(); // ensure phi_tilde is in [0, 2pi]
            }
          } else {
            values[kMCPhiTildeHE] = values[kMCPhiHE] - 0.75 * TMath::Pi(); // phi_tilde = phi - 3pi/4
            if (values[kMCPhiTildeHE] < 0) {
              values[kMCPhiTildeHE] += 2 * TMath::Pi(); // ensure phi_tilde is in [0, 2pi]
            }
          }
        } else {
          values[kMCPhiTildeHE] = -999; // not computable
        }
      }
    }

    if (useCS) {
      ROOT::Math::XYZVectorF zaxis_CS{(Beam1_CM - Beam2_CM).Unit()};
      ROOT::Math::XYZVectorF yaxis_CS{(Beam1_CM.Cross(Beam2_CM)).Unit()};
      ROOT::Math::XYZVectorF xaxis_CS{(yaxis_CS.Cross(zaxis_CS)).Unit()};
      if (fgUsedVars[kMCCosThetaCS])
        values[kMCCosThetaCS] = zaxis_CS.Dot(v_CM);
      if (fgUsedVars[kMCPhiCS]) {
        values[kMCPhiCS] = TMath::ATan2(yaxis_CS.Dot(v_CM), xaxis_CS.Dot(v_CM));
        if (values[kMCPhiCS] < 0) {
          values[kMCPhiCS] += 2 * TMath::Pi(); // ensure phi is in [0, 2pi]
        }
      }
      if (fgUsedVars[kMCPhiTildeCS]) {
        if (fgUsedVars[kMCCosThetaCS] && fgUsedVars[kMCPhiCS]) {
          if (values[kMCCosThetaCS] > 0) {
            values[kMCPhiTildeCS] = values[kMCPhiCS] - 0.25 * TMath::Pi(); // phi_tilde = phi - pi/4
            if (values[kMCPhiTildeCS] < 0) {
              values[kMCPhiTildeCS] += 2 * TMath::Pi(); // ensure phi_tilde is in [0, 2pi]
            }
          } else {
            values[kMCPhiTildeCS] = values[kMCPhiCS] - 0.75 * TMath::Pi(); // phi_tilde = phi - 3pi/4
            if (values[kMCPhiTildeCS] < 0) {
              values[kMCPhiTildeCS] += 2 * TMath::Pi(); // ensure phi_tilde is in [0, 2pi]
            }
          }
        } else {
          values[kMCPhiTildeCS] = -999; // not computable
        }
      }
    }

    if (usePP) {
      ROOT::Math::XYZVector zaxis_PP = ROOT::Math::XYZVector(v12.Py(), -v12.Px(), 0.f);
      ROOT::Math::XYZVector yaxis_PP{v12.Vect().Unit()};
      ROOT::Math::XYZVector xaxis_PP{(yaxis_PP.Cross(zaxis_PP)).Unit()};
      if (fgUsedVars[kMCCosThetaPP]) {
        values[kMCCosThetaPP] = zaxis_PP.Dot(v_CM);
      }
      if (fgUsedVars[kMCPhiPP]) {
        values[kMCPhiPP] = TMath::ATan2(yaxis_PP.Dot(v_CM), xaxis_PP.Dot(v_CM));
        if (values[kMCPhiPP] < 0) {
          values[kMCPhiPP] += 2 * TMath::Pi(); // ensure phi is in [0, 2pi]
        }
      }
      if (fgUsedVars[kMCPhiTildePP]) {
        if (fgUsedVars[kMCCosThetaPP] && fgUsedVars[kMCPhiPP]) {
          if (values[kMCCosThetaPP] > 0) {
            values[kMCPhiTildePP] = values[kMCPhiPP] - 0.25 * TMath::Pi(); // phi_tilde = phi - pi/4
            if (values[kMCPhiTildePP] < 0) {
              values[kMCPhiTildePP] += 2 * TMath::Pi(); // ensure phi_tilde is in [0, 2pi]
            }
          } else {
            values[kMCPhiTildePP] = values[kMCPhiPP] - 0.75 * TMath::Pi(); // phi_tilde = phi - 3pi/4
            if (values[kMCPhiTildePP] < 0) {
              values[kMCPhiTildePP] += 2 * TMath::Pi(); // ensure phi_tilde is in [0, 2pi]
            }
          }
        } else {
          values[kMCPhiTildePP] = -999; // not computable
        }
      }
    }

    if (useRM) {
      double randomCostheta = gRandom->Uniform(-1., 1.);
      double randomPhi = gRandom->Uniform(0., 2. * TMath::Pi());
      ROOT::Math::XYZVectorF zaxis_RM(randomCostheta, std::sqrt(1 - randomCostheta * randomCostheta) * std::cos(randomPhi), std::sqrt(1 - randomCostheta * randomCostheta) * std::sin(randomPhi));
      if (fgUsedVars[kMCCosThetaRM])
        values[kMCCosThetaRM] = zaxis_RM.Dot(v_CM);
    }
  }
}

template <int candidateType, typename T1, typename T2, typename T3>
void VarManager::FillTripleMC(T1 const& t1, T2 const& t2, T3 const& t3, float* values)
{
  if (!values) {
    values = fgValues;
  }

  if constexpr (candidateType == kTripleCandidateToEEPhoton) {
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

  if constexpr (candidateType == kTripleCandidateToKPiPi) {
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
    values[kCharge] = t1.sign() + t2.sign() + t3.sign();
    values[kS12] = (v1 + v2).M2();
    values[kS13] = (v1 + v3).M2();
    values[kS23] = (v2 + v3).M2();
  }
  if constexpr (candidateType == kBtoJpsiEEK) {
    float m1 = o2::constants::physics::MassElectron;
    float m2 = o2::constants::physics::MassElectron;
    float m3 = o2::constants::physics::MassKaonCharged;
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
    values[kMCY] = -v123.Rapidity();
    values[kPairMassDau] = v12.M();
    values[kPairPtDau] = v12.Pt();
    values[kRap] = -v123.Rapidity();
    values[kMassDau] = m3;
    values[VarManager::kDeltaMass] = v123.M() - v12.M();
    values[VarManager::kDeltaMass_jpsi] = v123.M() - v12.M() + m4;
    values[kPt] = t3.pt();
    values[kEta] = t3.eta();
    values[kEta1] = t1.eta();
    values[kEta2] = t2.eta();
    values[kDeltaEta] = v12.Eta();
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
      o2::track::TrackParCovFwd pars1 = FwdToTrackPar(t1, t1);
      o2::track::TrackParCovFwd pars2 = FwdToTrackPar(t2, t2);
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

      std::array<float, 6> covMatrixPCA{};
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
      values[kVertexingTauxyProjectedPoleJPsiMass] = values[kVertexingLxyProjected] * o2::constants::physics::MassJPsi / (v12.Pt());
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
        values[kVertexingTauxyProjectedPoleJPsiMass] = values[kVertexingLxyProjected] * o2::constants::physics::MassJPsi / (KFGeoTwoProng.GetPt());
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
          o2::track::TrackParCovFwd pars1 = FwdToTrackPar(t1, t1);
          o2::track::TrackParCovFwd pars2 = FwdToTrackPar(t2, t2);

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

      if (fgUsedVars[kVertexingChi2PCA]) {
        auto chi2PCA = fgFitterThreeProngBarrel.getChi2AtPCACandidate();
        values[VarManager::kVertexingChi2PCA] = chi2PCA;
      }

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

      values[kCosPointingAngle] = ((collision.posX() - secondaryVertex[0]) * v123.Px() +
                                   (collision.posY() - secondaryVertex[1]) * v123.Py() +
                                   (collision.posZ() - secondaryVertex[2]) * v123.Pz()) /
                                  (v123.P() * values[VarManager::kVertexingLxyz]);
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

      o2::track::TrackParCovFwd pars1 = FwdToTrackPar(lepton1, lepton1);
      o2::track::TrackParCovFwd pars2 = FwdToTrackPar(lepton2, lepton2);
      o2::track::TrackParCovFwd pars3 = FwdToTrackPar(track, track);

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
    values[VarManager::kPairRap] = -v123.Rapidity();
    values[VarManager::kPairEta] = v123.Eta();
    if (fgUsedVars[kPairMassDau] || fgUsedVars[kPairPtDau]) {
      values[VarManager::kPairMassDau] = v12.M();
      values[VarManager::kPairPtDau] = v12.Pt();
    }
    values[VarManager::kPt] = track.pt();
    values[kS12] = (v1 + v2).M2();
    values[kS13] = (v1 + v3).M2();
    values[kS23] = (v2 + v3).M2();

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

      if (fgUsedVars[kVertexingChi2PCA]) {
        auto chi2PCA = fgFitterThreeProngBarrel.getChi2AtPCACandidate();
        values[VarManager::kVertexingChi2PCA] = chi2PCA;
      }

      double phi = std::atan2(secondaryVertex[1] - collision.posY(), secondaryVertex[0] - collision.posX());
      double theta = std::atan2(secondaryVertex[2] - collision.posZ(),
                                std::sqrt((secondaryVertex[0] - collision.posX()) * (secondaryVertex[0] - collision.posX()) +
                                          (secondaryVertex[1] - collision.posY()) * (secondaryVertex[1] - collision.posY())));
      if (fgUsedVars[kVertexingLxy] || fgUsedVars[kVertexingLz] || fgUsedVars[kVertexingLxyz]) {

        values[VarManager::kVertexingLxy] = (collision.posX() - secondaryVertex[0]) * (collision.posX() - secondaryVertex[0]) +
                                            (collision.posY() - secondaryVertex[1]) * (collision.posY() - secondaryVertex[1]);
        values[VarManager::kVertexingLz] = (collision.posZ() - secondaryVertex[2]) * (collision.posZ() - secondaryVertex[2]);
        values[VarManager::kVertexingLxyz] = values[VarManager::kVertexingLxy] + values[VarManager::kVertexingLz];
        values[VarManager::kVertexingLxy] = std::sqrt(values[VarManager::kVertexingLxy]);
        values[VarManager::kVertexingLz] = std::sqrt(values[VarManager::kVertexingLz]);
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
      } // end eventHasVtxCov
    } // end (candidateType == kBtoJpsiEEK) && trackHasCov
  } // end KF
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
  values[kCORR4REF] = (pow(norm(compA21), 2) + norm(compA42) - 2. * (compA42 * conj(compA21) * conj(compA21)).real() + 8. * (compA23 * conj(compA21)).real() - 4. * S12A * norm(compA21) - 6. * S14A + 2. * S22A) / values[kM1111REF];
  values[kCORR2REF] = std::isnan(values[kM11REF]) || std::isinf(values[kM11REF]) || std::isnan(values[kCORR2REF]) || std::isinf(values[kCORR2REF]) ? 0 : values[kCORR2REF];
  values[kM11REF] = std::isnan(values[kM11REF]) || std::isinf(values[kM11REF]) || std::isnan(values[kCORR2REF]) || std::isinf(values[kCORR2REF]) ? 0 : values[kM11REF];
  values[kCORR4REF] = std::isnan(values[kM1111REF]) || std::isinf(values[kM1111REF]) || std::isnan(values[kCORR4REF]) || std::isinf(values[kCORR4REF]) ? 0 : values[kCORR4REF];
  values[kM1111REF] = std::isnan(values[kM1111REF]) || std::isinf(values[kM1111REF]) || std::isnan(values[kCORR4REF]) || std::isinf(values[kCORR4REF]) ? 0 : values[kM1111REF];
  values[kCORR2CORR4REF] = values[kCORR2REF] * values[kCORR4REF];
  values[kM11M1111REF] = values[kM11REF] * values[kM1111REF];

  // For cumulants: A = Full TPC, B = Negative TPC, C = Positive TPC
  complex<double> QA(values[kQ2X0A] * values[kS11A], values[kQ2Y0A] * values[kS11A]);
  complex<double> QB(values[kQ2X0B] * S11B, values[kQ2Y0B] * S11B);
  complex<double> QC(values[kQ2X0C] * S11C, values[kQ2Y0C] * S11C);
  values[kM11REFetagap] = S11B * S11C;
  values[kCORR2REFetagap] = ((QB * conj(QC)).real()) / values[kM11REFetagap];
  values[kCORR2REFetagap] = std::isnan(values[kM11REFetagap]) || std::isinf(values[kM11REFetagap]) || std::isnan(values[kCORR2REFetagap]) || std::isinf(values[kCORR2REFetagap]) ? 0 : values[kCORR2REFetagap];
  values[kM11REFetagap] = std::isnan(values[kM11REFetagap]) || std::isinf(values[kM11REFetagap]) || std::isnan(values[kCORR2REFetagap]) || std::isinf(values[kCORR2REFetagap]) ? 0 : values[kM11REFetagap];

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

  float xQVecFT0a = collision.qvecFT0ARe();   // already normalised
  float yQVecFT0a = collision.qvecFT0AIm();   // already normalised
  float xQVecFT0c = collision.qvecFT0CRe();   // already normalised
  float yQVecFT0c = collision.qvecFT0CIm();   // already normalised
  float xQVecFT0m = collision.qvecFT0MRe();   // already normalised
  float yQVecFT0m = collision.qvecFT0MIm();   // already normalised
  float xQVecFV0a = collision.qvecFV0ARe();   // already normalised
  float yQVecFV0a = collision.qvecFV0AIm();   // already normalised
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

  values[kR2EP_AB] = TMath::Cos(2 * getDeltaPsiInRange(epTPCFull, epFT0a, 2));
  values[kR2EP_AC] = TMath::Cos(2 * getDeltaPsiInRange(epTPCFull, epFT0c, 2));
  values[kR2EP_BC] = TMath::Cos(2 * getDeltaPsiInRange(epFT0a, epFT0c, 2));
  values[kR2EP_FT0CTPCPOS] = TMath::Cos(2 * getDeltaPsiInRange(epFT0c, epBPoss, 2));
  values[kR2EP_FT0CTPCNEG] = TMath::Cos(2 * getDeltaPsiInRange(epFT0c, epBNegs, 2));
  values[kR2EP_FT0ATPCPOS] = TMath::Cos(2 * getDeltaPsiInRange(epFT0a, epBPoss, 2));
  values[kR2EP_FT0ATPCNEG] = TMath::Cos(2 * getDeltaPsiInRange(epFT0a, epBNegs, 2));
  values[kR2EP_FT0MTPCPOS] = TMath::Cos(2 * getDeltaPsiInRange(epFT0m, epBPoss, 2));
  values[kR2EP_FT0MTPCNEG] = TMath::Cos(2 * getDeltaPsiInRange(epFT0m, epBNegs, 2));
  values[kR2EP_FV0ATPCPOS] = TMath::Cos(2 * getDeltaPsiInRange(epFV0a, epBPoss, 2));
  values[kR2EP_FV0ATPCNEG] = TMath::Cos(2 * getDeltaPsiInRange(epFV0a, epBNegs, 2));
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
  values[kPt1] = t1.pt();
  values[kPt2] = t2.pt();

  // TODO: provide different computations for vn
  // Compute the scalar product UQ using Q-vector from A, for second and third harmonic
  // Dilepton vn could be accessible after dividing this product with the R factor
  values[kU2Q2] = values[kQ2X0A] * TMath::Cos(2 * v12.Phi()) + values[kQ2Y0A] * TMath::Sin(2 * v12.Phi());
  values[kU3Q3] = values[kQ3X0A] * TMath::Cos(3 * v12.Phi()) + values[kQ3Y0A] * TMath::Sin(3 * v12.Phi());
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
  values[kCos2DeltaPhi] = TMath::Cos(2 * (v12.Phi() - Psi2A));
  values[kCos3DeltaPhi] = TMath::Cos(3 * (v12.Phi() - Psi3A));
  values[kR2EP_AB] = TMath::Cos(2 * (Psi2A - Psi2B));
  values[kR2EP_AC] = TMath::Cos(2 * (Psi2A - Psi2C));
  values[kR2EP_BC] = TMath::Cos(2 * (Psi2B - Psi2C));
  values[kR3EP] = TMath::Cos(3 * (Psi3B - Psi3C));

  float V2SP = values[kU2Q2] / values[kR2SP];
  float V2EP = values[kCos2DeltaPhi] / values[kR2EP];
  values[kV2SP] = std::isnan(V2SP) || std::isinf(V2SP) ? 0. : V2SP;
  values[kWV2SP] = std::isnan(V2SP) || std::isinf(V2SP) ? 0. : 1.0;
  values[kV2EP] = std::isnan(V2EP) || std::isinf(V2EP) ? 0. : V2EP;
  values[kWV2EP] = std::isnan(V2EP) || std::isinf(V2EP) ? 0. : 1.0;

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

  // global polarization parameters
  bool useGlobalPolarizatiobSpinOne = fgUsedVars[kCosThetaStarTPC] || fgUsedVars[kCosThetaStarFT0A] || fgUsedVars[kCosThetaStarFT0C];
  if (useGlobalPolarizatiobSpinOne) {
    ROOT::Math::Boost boostv12{v12.BoostToCM()};
    ROOT::Math::XYZVectorF v1_CM{(boostv12(v1).Vect()).Unit()};
    ROOT::Math::XYZVectorF v2_CM{(boostv12(v2).Vect()).Unit()};

    // using positive sign convention for the first track
    ROOT::Math::XYZVectorF v_CM = (t1.sign() > 0 ? v1_CM : v2_CM);

    ROOT::Math::XYZVector zaxisTPC = ROOT::Math::XYZVector(TMath::Cos(Psi2A), TMath::Sin(Psi2A), 0).Unit();
    values[kCosThetaStarTPC] = v_CM.Dot(zaxisTPC);

    ROOT::Math::XYZVector zaxisFT0A = ROOT::Math::XYZVector(TMath::Cos(Psi2B), TMath::Sin(Psi2B), 0).Unit();
    values[kCosThetaStarFT0A] = v_CM.Dot(zaxisFT0A);

    ROOT::Math::XYZVector zaxisFT0C = ROOT::Math::XYZVector(TMath::Cos(Psi2C), TMath::Sin(Psi2C), 0).Unit();
    values[kCosThetaStarFT0C] = v_CM.Dot(zaxisFT0C);
  }

  //  kV4, kC4POI, kC4REF etc.
  if constexpr ((fillMap & ReducedEventQvectorExtra) > 0) {
    complex<double> Q21(values[kQ2X0A] * values[kS11A], values[kQ2Y0A] * values[kS11A]);
    complex<double> Q42(values[kQ42XA], values[kQ42YA]);
    complex<double> Q23(values[kQ23XA], values[kQ23YA]);
    complex<double> P2(TMath::Cos(2 * v12.Phi()), TMath::Sin(2 * v12.Phi()));
    values[kM01POI] = values[kMultDimuons] * values[kS11A];
    values[kM0111POI] = values[kMultDimuons] * (values[kS31A] - 3. * values[kS11A] * values[kS12A] + 2. * values[kS13A]);
    values[kCORR2POI] = (P2 * conj(Q21)).real() / values[kM01POI];
    values[kCORR4POI] = (P2 * Q21 * conj(Q21) * conj(Q21) - P2 * Q21 * conj(Q42) - 2. * values[kS12A] * P2 * conj(Q21) + 2. * P2 * conj(Q23)).real() / values[kM0111POI];
    values[kM01POIoverMp] = values[kMultDimuons] > 0 && !(std::isnan(values[kM01POI]) || std::isinf(values[kM01POI]) || std::isnan(values[kCORR2POI]) || std::isinf(values[kCORR2POI]) || std::isnan(values[kM0111POI]) || std::isinf(values[kM0111POI]) || std::isnan(values[kCORR4POI]) || std::isinf(values[kCORR4POI])) ? values[kM01POI] / values[kMultDimuons] : 0;
    values[kM0111POIoverMp] = values[kMultDimuons] > 0 && !(std::isnan(values[kM0111POI]) || std::isinf(values[kM0111POI]) || std::isnan(values[kCORR4POI]) || std::isinf(values[kCORR4POI]) || std::isnan(values[kM01POI]) || std::isinf(values[kM01POI]) || std::isnan(values[kCORR2POI]) || std::isinf(values[kCORR2POI])) ? values[kM0111POI] / values[kMultDimuons] : 0;
    values[kM11REFoverMp] = values[kMultDimuons] > 0 && !(std::isnan(values[kM11REF]) || std::isinf(values[kM11REF]) || std::isnan(values[kCORR2REF]) || std::isinf(values[kCORR2REF]) || std::isnan(values[kM1111REF]) || std::isinf(values[kM1111REF]) || std::isnan(values[kCORR4REF]) || std::isinf(values[kCORR4REF])) ? values[kM11REF] / values[kMultDimuons] : 0;
    values[kM1111REFoverMp] = values[kMultDimuons] > 0 && !(std::isnan(values[kM1111REF]) || std::isinf(values[kM1111REF]) || std::isnan(values[kCORR4REF]) || std::isinf(values[kCORR4REF]) || std::isnan(values[kM11REF]) || std::isinf(values[kM11REF]) || std::isnan(values[kCORR2REF]) || std::isinf(values[kCORR2REF])) ? values[kM1111REF] / values[kMultDimuons] : 0;
    values[kCORR2REFbydimuons] = std::isnan(values[kM11REFoverMp]) || std::isinf(values[kM11REFoverMp]) || std::isnan(values[kCORR2REF]) || std::isinf(values[kCORR2REF]) || std::isnan(values[kM1111REFoverMp]) || std::isinf(values[kM1111REFoverMp]) || std::isnan(values[kCORR4REF]) || std::isinf(values[kCORR4REF]) ? 0 : values[kCORR2REF];
    values[kCORR4REFbydimuons] = std::isnan(values[kM1111REFoverMp]) || std::isinf(values[kM1111REFoverMp]) || std::isnan(values[kCORR4REF]) || std::isinf(values[kCORR4REF]) || std::isnan(values[kM11REFoverMp]) || std::isinf(values[kM11REFoverMp]) || std::isnan(values[kCORR2REF]) || std::isinf(values[kCORR2REF]) ? 0 : values[kCORR4REF];
    values[kCORR2POI] = std::isnan(values[kCORR2POI]) || std::isinf(values[kCORR2POI]) || std::isnan(values[kM01POI]) || std::isinf(values[kM01POI]) || std::isnan(values[kCORR4POI]) || std::isinf(values[kCORR4POI]) || std::isnan(values[kM0111POI]) || std::isinf(values[kM0111POI]) ? 0 : values[kCORR2POI];
    values[kCORR4POI] = std::isnan(values[kCORR4POI]) || std::isinf(values[kCORR4POI]) || std::isnan(values[kM0111POI]) || std::isinf(values[kM0111POI]) || std::isnan(values[kCORR2POI]) || std::isinf(values[kCORR2POI]) || std::isnan(values[kM01POI]) || std::isinf(values[kM01POI]) ? 0 : values[kCORR4POI];
    values[kCORR2CORR4REF] = std::isnan(values[kM11M1111REFoverMp]) || std::isinf(values[kM11M1111REFoverMp]) || std::isnan(values[kCORR2CORR4REF]) || std::isinf(values[kCORR2CORR4REF]) ? 0 : values[kCORR2CORR4REF];
    values[kCORR2POICORR4POI] = std::isnan(values[kCORR2POI]) || std::isinf(values[kCORR2POI]) || std::isnan(values[kCORR4POI]) || std::isinf(values[kCORR4POI]) || std::isnan(values[kM01POI]) || std::isinf(values[kM01POI]) || std::isnan(values[kM0111POI]) || std::isinf(values[kM0111POI]) ? 0 : values[kCORR2POI] * values[kCORR4POI];
    values[kCORR2REFCORR4POI] = std::isnan(values[kCORR2REF]) || std::isinf(values[kCORR2REF]) || std::isnan(values[kCORR4POI]) || std::isinf(values[kCORR4POI]) || std::isnan(values[kM11REF]) || std::isinf(values[kM11REF]) || std::isnan(values[kM0111POI]) || std::isinf(values[kM0111POI]) ? 0 : values[kCORR2REF] * values[kCORR4POI];
    values[kCORR2REFCORR2POI] = std::isnan(values[kCORR2REF]) || std::isinf(values[kCORR2REF]) || std::isnan(values[kCORR2POI]) || std::isinf(values[kCORR2POI]) || std::isnan(values[kM11REF]) || std::isinf(values[kM11REF]) || std::isnan(values[kM01POI]) || std::isinf(values[kM01POI]) ? 0 : values[kCORR2REF] * values[kCORR2POI];
    values[kM11M1111REFoverMp] = values[kMultDimuons] > 0 && !(std::isnan(values[kM11M1111REF]) || std::isinf(values[kM11M1111REF]) || std::isnan(values[kCORR2CORR4REF]) || std::isinf(values[kCORR2CORR4REF])) ? values[kM11M1111REF] / values[kMultDimuons] : 0;
    values[kM01M0111overMp] = values[kMultDimuons] > 0 && !(std::isnan(values[kM01POI]) || std::isinf(values[kM01POI]) || std::isnan(values[kM0111POI]) || std::isinf(values[kM0111POI]) || std::isnan(values[kCORR2POI]) || std::isinf(values[kCORR2POI]) || std::isnan(values[kCORR4POI]) || std::isinf(values[kCORR4POI])) ? (values[kM01POI] * values[kM0111POI]) / values[kMultDimuons] : 0;
    values[kM11M0111overMp] = values[kMultDimuons] > 0 && !(std::isnan(values[kM11REF]) || std::isinf(values[kM11REF]) || std::isnan(values[kM0111POI]) || std::isinf(values[kM0111POI]) || std::isnan(values[kCORR2REF]) || std::isinf(values[kCORR2REF]) || std::isnan(values[kCORR4POI]) || std::isinf(values[kCORR4POI])) ? (values[kM11REF] * values[kM0111POI]) / values[kMultDimuons] : 0;
    values[kM11M01overMp] = values[kMultDimuons] > 0 && !(std::isnan(values[kM11REF]) || std::isinf(values[kM11REF]) || std::isnan(values[kM01POI]) || std::isinf(values[kM01POI]) || std::isnan(values[kCORR2REF]) || std::isinf(values[kCORR2REF]) || std::isnan(values[kCORR2POI]) || std::isinf(values[kCORR2POI])) ? (values[kM11REF] * values[kM01POI]) / values[kMultDimuons] : 0;

    complex<double> P2plus(TMath::Cos(2 * v1.Phi()), TMath::Sin(2 * v1.Phi()));
    complex<double> P2minus(TMath::Cos(2 * v2.Phi()), TMath::Sin(2 * v2.Phi()));
    values[kM11REFoverMpplus] = values[kMultAntiMuons] > 0 && !(std::isnan(values[kM11REF]) || std::isinf(values[kM11REF]) || std::isnan(values[kCORR2REF]) || std::isinf(values[kCORR2REF]) || std::isnan(values[kM1111REF]) || std::isinf(values[kM1111REF]) || std::isnan(values[kCORR4REF]) || std::isinf(values[kCORR4REF])) ? values[kM11REF] / values[kMultAntiMuons] : 0;
    values[kM1111REFoverMpplus] = values[kMultAntiMuons] > 0 && !(std::isnan(values[kM11REF]) || std::isinf(values[kM11REF]) || std::isnan(values[kCORR2REF]) || std::isinf(values[kCORR2REF]) || std::isnan(values[kM1111REF]) || std::isinf(values[kM1111REF]) || std::isnan(values[kCORR4REF]) || std::isinf(values[kCORR4REF])) ? values[kM1111REF] / values[kMultAntiMuons] : 0;
    values[kM11REFoverMpminus] = values[kMultMuons] > 0 && !(std::isnan(values[kM11REF]) || std::isinf(values[kM11REF]) || std::isnan(values[kCORR2REF]) || std::isinf(values[kCORR2REF]) || std::isnan(values[kM1111REF]) || std::isinf(values[kM1111REF]) || std::isnan(values[kCORR4REF]) || std::isinf(values[kCORR4REF])) ? values[kM11REF] / values[kMultMuons] : 0;
    values[kM1111REFoverMpminus] = values[kMultMuons] > 0 && !(std::isnan(values[kM11REF]) || std::isinf(values[kM11REF]) || std::isnan(values[kCORR2REF]) || std::isinf(values[kCORR2REF]) || std::isnan(values[kM1111REF]) || std::isinf(values[kM1111REF]) || std::isnan(values[kCORR4REF]) || std::isinf(values[kCORR4REF])) ? values[kM1111REF] / values[kMultMuons] : 0;
    values[kCORR2POIplus] = (P2plus * conj(Q21)).real() / values[kM01POI];
    values[kCORR2POIminus] = (P2minus * conj(Q21)).real() / values[kM01POI];
    values[kM01POIplus] = values[kMultAntiMuons] * values[kS11A];
    values[kM0111POIplus] = values[kMultAntiMuons] * (values[kS31A] - 3. * values[kS11A] * values[kS12A] + 2. * values[kS13A]);
    values[kCORR2POIplus] = (P2plus * conj(Q21)).real() / values[kM01POIplus];
    values[kCORR4POIplus] = (P2plus * Q21 * conj(Q21) * conj(Q21) - P2plus * Q21 * conj(Q42) - 2. * values[kS12A] * P2plus * conj(Q21) + 2. * P2plus * conj(Q23)).real() / values[kM0111POIplus];
    values[kM01POIminus] = values[kMultMuons] * values[kS11A];
    values[kM0111POIminus] = values[kMultMuons] * (values[kS31A] - 3. * values[kS11A] * values[kS12A] + 2. * values[kS13A]);
    values[kCORR2POIminus] = (P2minus * conj(Q21)).real() / values[kM01POIminus];
    values[kCORR4POIminus] = (P2minus * Q21 * conj(Q21) * conj(Q21) - P2minus * Q21 * conj(Q42) - 2. * values[kS12A] * P2minus * conj(Q21) + 2. * P2minus * conj(Q23)).real() / values[kM0111POIminus];
    values[kM01POIplus] = std::isnan(values[kM01POIplus]) || std::isinf(values[kM01POIplus]) || std::isnan(values[kM0111POIplus]) || std::isinf(values[kM0111POIplus]) || std::isnan(values[kCORR2POIplus]) || std::isinf(values[kCORR2POIplus]) || std::isnan(values[kCORR4POIplus]) || std::isinf(values[kCORR4POIplus]) ? 0 : values[kM01POIplus];
    values[kM0111POIplus] = std::isnan(values[kM0111POIplus]) || std::isinf(values[kM0111POIplus]) || std::isnan(values[kCORR2POIplus]) || std::isinf(values[kCORR2POIplus]) || std::isnan(values[kCORR4POIplus]) || std::isinf(values[kCORR4POIplus]) ? 0 : values[kM0111POIplus];
    values[kCORR2POIplus] = std::isnan(values[kM01POIplus]) || std::isinf(values[kM01POIplus]) || std::isnan(values[kM0111POIplus]) || std::isinf(values[kM0111POIplus]) || std::isnan(values[kCORR2POIplus]) || std::isinf(values[kCORR2POIplus]) || std::isnan(values[kCORR4POIplus]) || std::isinf(values[kCORR4POIplus]) ? 0 : values[kCORR2POIplus];
    values[kCORR4POIplus] = std::isnan(values[kM01POIplus]) || std::isinf(values[kM01POIplus]) || std::isnan(values[kM0111POIplus]) || std::isinf(values[kM0111POIplus]) || std::isnan(values[kCORR2POIplus]) || std::isinf(values[kCORR2POIplus]) || std::isnan(values[kCORR4POIplus]) || std::isinf(values[kCORR4POIplus]) ? 0 : values[kCORR4POIplus];
    values[kM01POIminus] = std::isnan(values[kM01POIminus]) || std::isinf(values[kM01POIminus]) || std::isnan(values[kM0111POIminus]) || std::isinf(values[kM0111POIminus]) || std::isnan(values[kCORR2POIminus]) || std::isinf(values[kCORR2POIminus]) || std::isnan(values[kCORR4POIminus]) || std::isinf(values[kCORR4POIminus]) ? 0 : values[kM01POIminus];
    values[kM0111POIminus] = std::isnan(values[kM01POIminus]) || std::isinf(values[kM01POIminus]) || std::isnan(values[kM0111POIminus]) || std::isinf(values[kM0111POIminus]) || std::isnan(values[kCORR2POIminus]) || std::isinf(values[kCORR2POIminus]) || std::isnan(values[kCORR4POIminus]) || std::isinf(values[kCORR4POIminus]) ? 0 : values[kM0111POIminus];
    values[kCORR2POIminus] = std::isnan(values[kM01POIminus]) || std::isinf(values[kM01POIminus]) || std::isnan(values[kM0111POIminus]) || std::isinf(values[kM0111POIminus]) || std::isnan(values[kCORR2POIminus]) || std::isinf(values[kCORR2POIminus]) || std::isnan(values[kCORR4POIminus]) || std::isinf(values[kCORR4POIminus]) ? 0 : values[kCORR2POIminus];
    values[kCORR4POIminus] = std::isnan(values[kM01POIminus]) || std::isinf(values[kM01POIminus]) || std::isnan(values[kM0111POIminus]) || std::isinf(values[kM0111POIminus]) || std::isnan(values[kCORR2POIminus]) || std::isinf(values[kCORR2POIminus]) || std::isnan(values[kCORR4POIminus]) || std::isinf(values[kCORR4POIminus]) ? 0 : values[kCORR4POIminus];
    values[kM01POIoverMpminus] = values[kMultMuons] > 0 && !(std::isnan(values[kM0111POIminus]) || std::isinf(values[kM0111POIminus]) || std::isnan(values[kCORR4POIminus]) || std::isinf(values[kCORR4POIminus]) || std::isnan(values[kM01POIminus]) || std::isinf(values[kM01POIminus]) || std::isnan(values[kCORR2POIminus]) || std::isinf(values[kCORR2POIminus])) ? values[kM01POIminus] / values[kMultMuons] : 0;
    values[kM0111POIoverMpminus] = values[kMultMuons] > 0 && !(std::isnan(values[kM0111POIminus]) || std::isinf(values[kM0111POIminus]) || std::isnan(values[kCORR4POIminus]) || std::isinf(values[kCORR4POIminus]) || std::isnan(values[kM01POIminus]) || std::isinf(values[kM01POIminus]) || std::isnan(values[kCORR2POIminus]) || std::isinf(values[kCORR2POIminus])) ? values[kM0111POIminus] / values[kMultMuons] : 0;
    values[kM01POIoverMpplus] = values[kMultAntiMuons] > 0 && !(std::isnan(values[kM0111POIplus]) || std::isinf(values[kM0111POIplus]) || std::isnan(values[kCORR4POIplus]) || std::isinf(values[kCORR4POIplus]) || std::isnan(values[kM01POIplus]) || std::isinf(values[kM01POIplus]) || std::isnan(values[kCORR2POIplus]) || std::isinf(values[kCORR2POIplus])) ? values[kM01POIplus] / values[kMultAntiMuons] : 0;
    values[kM0111POIoverMpplus] = values[kMultMuons] > 0 && !(std::isnan(values[kM0111POIplus]) || std::isinf(values[kM0111POIplus]) || std::isnan(values[kCORR4POIplus]) || std::isinf(values[kCORR4POIplus]) || std::isnan(values[kM01POIplus]) || std::isinf(values[kM01POIplus]) || std::isnan(values[kCORR2POIplus]) || std::isinf(values[kCORR2POIplus])) ? values[kM0111POIplus] / values[kMultAntiMuons] : 0;
  }

  ROOT::Math::PtEtaPhiMVector v1_vp(v1.Pt(), v1.Eta(), v1.Phi() - Psi2B, v1.M());
  ROOT::Math::PtEtaPhiMVector v2_vp(v2.Pt(), v2.Eta(), v2.Phi() - Psi2B, v2.M());
  ROOT::Math::PtEtaPhiMVector v12_vp = v1_vp + v2_vp;
  auto p12_vp = ROOT::Math::XYZVectorF(v12_vp.Px(), v12_vp.Py(), v12_vp.Pz());
  auto p12_vp_projXZ = ROOT::Math::XYZVectorF(p12_vp.X(), 0, p12_vp.Z());
  auto vDimu = (t1.sign() > 0 ? ROOT::Math::XYZVectorF(v1_vp.Px(), v1_vp.Py(), v1_vp.Pz()).Cross(ROOT::Math::XYZVectorF(v2_vp.Px(), v2_vp.Py(), v2_vp.Pz()))
                              : ROOT::Math::XYZVectorF(v2_vp.Px(), v2_vp.Py(), v2_vp.Pz()).Cross(ROOT::Math::XYZVectorF(v1_vp.Px(), v1_vp.Py(), v1_vp.Pz())));
  auto vRef = p12_vp.Cross(p12_vp_projXZ);
  values[kCosPhiVP] = vDimu.Dot(vRef) / (vRef.R() * vDimu.R());
  values[kPhiVP] = std::acos(vDimu.Dot(vRef) / (vRef.R() * vDimu.R()));
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

  if (fgUsedVars[kPairMass] || fgUsedVars[kPairPt] || fgUsedVars[kPairEta] || fgUsedVars[kPairPhi] || fgUsedVars[kPairMassDau] || fgUsedVars[kPairPtDau] || fgUsedVars[kDileptonHadronKstar]) {
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
    // Calculate kstar of Dilepton and hadron pair
    ROOT::Math::PtEtaPhiMVector v12_Qvect = v1 - v2;
    double Pinv = v12.M();
    double Q1 = (dilepton.mass() * dilepton.mass() - hadronMass * hadronMass) / Pinv;
    values[kDileptonHadronKstar] = sqrt(Q1 * Q1 - v12_Qvect.M2()) / 2.0;
  }
  if (fgUsedVars[kCosChi] || fgUsedVars[kECWeight] || fgUsedVars[kCosTheta] || fgUsedVars[kEWeight_before] || fgUsedVars[kPtDau] || fgUsedVars[kEtaDau] || fgUsedVars[kPhiDau] || fgUsedVars[kCosChi_randomPhi_trans] || fgUsedVars[kCosChi_randomPhi_toward] || fgUsedVars[kCosChi_randomPhi_away]) {
    ROOT::Math::PtEtaPhiMVector v1(dilepton.pt(), dilepton.eta(), dilepton.phi(), dilepton.mass());
    ROOT::Math::PtEtaPhiMVector v2(hadron.pt(), hadron.eta(), hadron.phi(), o2::constants::physics::MassPionCharged);
    values[kCosChi] = LorentzTransformJpsihadroncosChi("coschi", v1, v2);
    float E_boost = LorentzTransformJpsihadroncosChi("weight_boost", v1, v2);
    values[kECWeight] = E_boost / v1.M();
    values[kCosTheta] = LorentzTransformJpsihadroncosChi("costheta", v1, v2);
    values[kEWeight_before] = v2.Pt() / v1.M();
    values[kPtDau] = v2.pt();
    values[kEtaDau] = v2.eta();
    values[kPhiDau] = RecoDecay::constrainAngle(v2.phi(), -o2::constants::math::PIHalf);

    float deltaphi = RecoDecay::constrainAngle(v1.phi() - v2.phi(), -o2::constants::math::PIHalf);
    values[kCosChi_randomPhi_trans] = -999.9f;
    values[kCosChi_randomPhi_toward] = -999.9f;
    values[kCosChi_randomPhi_away] = -999.9f;

    values[kdeltaphi_randomPhi_trans] = -999.9f;
    values[kdeltaphi_randomPhi_toward] = -999.9f;
    values[kdeltaphi_randomPhi_away] = -999.9f;

    float randomPhi_trans = -o2::constants::math::PIHalf;
    float randomPhi_toward = -o2::constants::math::PIHalf;
    float randomPhi_away = -o2::constants::math::PIHalf;

    if ((deltaphi > -0.5 * TMath::Pi() && deltaphi < -1. / 3 * TMath::Pi()) || (deltaphi > 4. / 3 * TMath::Pi() && deltaphi < 1.5 * TMath::Pi()) || (deltaphi > 1. / 3 * TMath::Pi() && deltaphi < 2. / 3 * TMath::Pi())) {
      randomPhi_trans = gRandom->Uniform(-o2::constants::math::PIHalf, 3. * o2::constants::math::PIHalf);
      randomPhi_toward = gRandom->Uniform(-o2::constants::math::PIHalf, 3. * o2::constants::math::PIHalf);
      randomPhi_away = gRandom->Uniform(-o2::constants::math::PIHalf, 3. * o2::constants::math::PIHalf);

      ROOT::Math::PtEtaPhiMVector v2_randomPhi_trans(v2.pt(), v2.eta(), randomPhi_trans, o2::constants::physics::MassPionCharged);
      values[kCosChi_randomPhi_trans] = LorentzTransformJpsihadroncosChi("coschi", v1, v2_randomPhi_trans);
      values[kWeight_randomPhi_trans] = LorentzTransformJpsihadroncosChi("weight_boost", v1, v2_randomPhi_trans) / v1.M();

      ROOT::Math::PtEtaPhiMVector v2_randomPhi_toward(v2.pt(), v2.eta(), randomPhi_toward, o2::constants::physics::MassPionCharged);
      values[kCosChi_randomPhi_toward] = LorentzTransformJpsihadroncosChi("coschi", v1, v2_randomPhi_toward);
      values[kWeight_randomPhi_toward] = LorentzTransformJpsihadroncosChi("weight_boost", v1, v2_randomPhi_toward) / v1.M();

      ROOT::Math::PtEtaPhiMVector v2_randomPhi_away(v2.pt(), v2.eta(), randomPhi_away, o2::constants::physics::MassPionCharged);
      values[kCosChi_randomPhi_away] = LorentzTransformJpsihadroncosChi("coschi", v1, v2_randomPhi_away);
      values[kWeight_randomPhi_away] = LorentzTransformJpsihadroncosChi("weight_boost", v1, v2_randomPhi_away) / v1.M();

      values[kdeltaphi_randomPhi_trans] = RecoDecay::constrainAngle(v1.phi() - randomPhi_trans, -o2::constants::math::PIHalf);
      values[kdeltaphi_randomPhi_toward] = RecoDecay::constrainAngle(v1.phi() - randomPhi_toward, -o2::constants::math::PIHalf);
      values[kdeltaphi_randomPhi_away] = RecoDecay::constrainAngle(v1.phi() - randomPhi_away, -o2::constants::math::PIHalf);
    }
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
    values[kDeltaR1] = ROOT::Math::VectorUtil::DeltaR(v1, v2);
    values[kDeltaR2] = ROOT::Math::VectorUtil::DeltaR(v1, v3);
    values[kDeltaR] = sqrt(pow(values[kDeltaR1], 2) + pow(values[kDeltaR2], 2));
    values[kRap] = v123.Rapidity();
  }
}

//__________________________________________________________________
template <int candidateType, uint32_t collFillMap, uint32_t fillMap, typename C, typename T1>
void VarManager::FillDileptonTrackTrackVertexing(C const& collision, T1 const& lepton1, T1 const& lepton2, T1 const& track1, T1 const& track2, float* values)
{
  constexpr bool eventHasVtxCov = ((collFillMap & Collision) > 0 || (collFillMap & ReducedEventVtxCov) > 0);
  constexpr bool trackHasCov = ((fillMap & TrackCov) > 0 || (fillMap & ReducedTrackBarrelCov) > 0);

  if (!eventHasVtxCov || !trackHasCov) {
    return;
  }

  if (!values) {
    values = fgValues;
  }

  float mtrack1, mtrack2;
  float mlepton1, mlepton2;

  if constexpr (candidateType == kXtoJpsiPiPi || candidateType == kPsi2StoJpsiPiPi) {
    mlepton1 = o2::constants::physics::MassElectron;
    mlepton2 = o2::constants::physics::MassElectron;
    mtrack1 = o2::constants::physics::MassPionCharged;
    mtrack2 = o2::constants::physics::MassPionCharged;
  }

  ROOT::Math::PtEtaPhiMVector v1(lepton1.pt(), lepton1.eta(), lepton1.phi(), mlepton1);
  ROOT::Math::PtEtaPhiMVector v2(lepton2.pt(), lepton2.eta(), lepton2.phi(), mlepton2);
  ROOT::Math::PtEtaPhiMVector v3(track1.pt(), track1.eta(), track1.phi(), mtrack1);
  ROOT::Math::PtEtaPhiMVector v4(track2.pt(), track2.eta(), track2.phi(), mtrack2);
  ROOT::Math::PtEtaPhiMVector v1234 = v1 + v2 + v3 + v4;

  int procCodeDilepton = 0;
  int procCodeDileptonTrackTrack = 0;

  values[kUsedKF] = fgUsedKF;
  if (!fgUsedKF) {
    // create covariance matrix
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
    std::array<float, 5> track1pars = {track1.y(), track1.z(), track1.snp(), track1.tgl(), track1.signed1Pt()};
    std::array<float, 15> track1covs = {track1.cYY(), track1.cZY(), track1.cZZ(), track1.cSnpY(), track1.cSnpZ(),
                                        track1.cSnpSnp(), track1.cTglY(), track1.cTglZ(), track1.cTglSnp(), track1.cTglTgl(),
                                        track1.c1PtY(), track1.c1PtZ(), track1.c1PtSnp(), track1.c1PtTgl(), track1.c1Pt21Pt2()};
    o2::track::TrackParCov pars3{track1.x(), track1.alpha(), track1pars, track1covs};
    std::array<float, 5> track2pars = {track2.y(), track2.z(), track2.snp(), track2.tgl(), track2.signed1Pt()};
    std::array<float, 15> track2covs = {track2.cYY(), track2.cZY(), track2.cZZ(), track2.cSnpY(), track2.cSnpZ(),
                                        track2.cSnpSnp(), track2.cTglY(), track2.cTglZ(), track2.cTglSnp(), track2.cTglTgl(),
                                        track2.c1PtY(), track2.c1PtZ(), track2.c1PtSnp(), track2.c1PtTgl(), track2.c1Pt21Pt2()};
    o2::track::TrackParCov pars4{track2.x(), track2.alpha(), track2pars, track2covs};

    procCodeDilepton = VarManager::fgFitterTwoProngBarrel.process(pars1, pars2);
    // create dilepton track
    // o2::track::TrackParCov parsDilepton = VarManager::fgFitterTwoProngBarrel.createParentTrackParCov(0);
    // procCodeDileptonTrackTrack = VarManager::fgFitterThreeProngBarrel.process(parsDilepton, pars3, pars4);
    procCodeDileptonTrackTrack = VarManager::fgFitterFourProngBarrel.process(pars1, pars2, pars3, pars4);

    // fill values
    if (procCodeDilepton == 0 && procCodeDileptonTrackTrack == 0) {
      // TODO: set the other variables to appropriate values and return
      values[kVertexingLxy] = -999.;
      values[kVertexingLxyz] = -999.;
      values[kVertexingLz] = -999.;
      values[kVertexingLxyErr] = -999.;
      values[kVertexingLxyzErr] = -999.;
      values[kVertexingLzErr] = -999.;
      values[kVertexingTauxy] = -999.;
      values[kVertexingTauxyErr] = -999.;
      values[kVertexingTauz] = -999.;
      values[kVertexingTauzErr] = -999.;
      values[kVertexingLzProjected] = -999.;
      values[kVertexingLxyProjected] = -999.;
      values[kVertexingLxyzProjected] = -999.;
      values[kVertexingTauzProjected] = -999.;
      values[kVertexingTauxyProjected] = -999.;
      values[kVertexingTauxyzProjected] = -999.;
      return;
    } else {
      Vec3D secondaryVertex;
      std::array<float, 6> covMatrixPCA;
      secondaryVertex = fgFitterFourProngBarrel.getPCACandidate();
      covMatrixPCA = fgFitterFourProngBarrel.calcPCACovMatrixFlat();

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

      values[kVertexingTauz] = (collision.posZ() - secondaryVertex[2]) * v1234.M() / (TMath::Abs(v1234.Pz()) * o2::constants::physics::LightSpeedCm2NS);
      values[kVertexingTauxy] = values[kVertexingLxy] * v1234.M() / (v1234.Pt() * o2::constants::physics::LightSpeedCm2NS);

      values[kVertexingTauzErr] = values[kVertexingLzErr] * v1234.M() / (TMath::Abs(v1234.Pz()) * o2::constants::physics::LightSpeedCm2NS);
      values[kVertexingTauxyErr] = values[kVertexingLxyErr] * v1234.M() / (v1234.Pt() * o2::constants::physics::LightSpeedCm2NS);

      values[kCosPointingAngle] = ((collision.posX() - secondaryVertex[0]) * v1234.Px() +
                                   (collision.posY() - secondaryVertex[1]) * v1234.Py() +
                                   (collision.posZ() - secondaryVertex[2]) * v1234.Pz()) /
                                  (v1234.P() * values[VarManager::kVertexingLxyz]);
      // // run 2 definitions: Decay length projected onto the momentum vector of the candidate
      values[kVertexingLzProjected] = (secondaryVertex[2] - collision.posZ()) * v1234.Pz();
      values[kVertexingLzProjected] = values[kVertexingLzProjected] / TMath::Sqrt(v1234.Pz() * v1234.Pz());
      values[kVertexingLxyProjected] = ((secondaryVertex[0] - collision.posX()) * v1234.Px()) + ((secondaryVertex[1] - collision.posY()) * v1234.Py());
      values[kVertexingLxyProjected] = values[kVertexingLxyProjected] / TMath::Sqrt((v1234.Px() * v1234.Px()) + (v1234.Py() * v1234.Py()));
      values[kVertexingLxyzProjected] = ((secondaryVertex[0] - collision.posX()) * v1234.Px()) + ((secondaryVertex[1] - collision.posY()) * v1234.Py()) + ((secondaryVertex[2] - collision.posZ()) * v1234.Pz());
      values[kVertexingLxyzProjected] = values[kVertexingLxyzProjected] / TMath::Sqrt((v1234.Px() * v1234.Px()) + (v1234.Py() * v1234.Py()) + (v1234.Pz() * v1234.Pz()));

      values[kVertexingTauzProjected] = values[kVertexingLzProjected] * v1234.M() / TMath::Abs(v1234.Pz());
      values[kVertexingTauxyProjected] = values[kVertexingLxyProjected] * v1234.M() / (v1234.Pt());
      values[kVertexingTauxyzProjected] = values[kVertexingLxyzProjected] * v1234.M() / (v1234.P());
    }
  } else if (fgUsedKF) {
    KFParticle lepton1KF; // lepton1
    KFParticle lepton2KF; // lepton2
    KFParticle KFGeoTwoLeptons;
    KFParticle trk1KF; // track1
    KFParticle trk2KF; // track2
    KFParticle KFGeoTwoTracks;
    KFParticle KFGeoFourProng;
    if constexpr (candidateType == kXtoJpsiPiPi) {
      KFPTrack kfpTrack0 = createKFPTrackFromTrack(lepton1);
      lepton1KF = KFParticle(kfpTrack0, -11 * lepton1.sign());
      KFPTrack kfpTrack1 = createKFPTrackFromTrack(lepton2);
      lepton2KF = KFParticle(kfpTrack1, -11 * lepton2.sign());
      KFPTrack kfpTrack2 = createKFPTrackFromTrack(track1);
      trk1KF = KFParticle(kfpTrack2, 211 * track1.sign());
      KFPTrack kfpTrack3 = createKFPTrackFromTrack(track2);
      trk2KF = KFParticle(kfpTrack3, 211 * track2.sign());

      KFGeoTwoLeptons.SetConstructMethod(2);
      KFGeoTwoLeptons.AddDaughter(lepton1KF);
      KFGeoTwoLeptons.AddDaughter(lepton2KF);

      if (fgUsedVars[kPairMass] || fgUsedVars[kPairPt]) {
        values[VarManager::kPairMass] = KFGeoTwoLeptons.GetMass();
        values[VarManager::kPairPt] = KFGeoTwoLeptons.GetPt();
      }

      KFGeoTwoTracks.SetConstructMethod(2);
      KFGeoTwoTracks.AddDaughter(trk1KF);
      KFGeoTwoTracks.AddDaughter(trk2KF);

      if (fgUsedVars[kDitrackMass] || fgUsedVars[kDitrackPt]) {
        values[VarManager::kDitrackMass] = KFGeoTwoTracks.GetMass();
        values[VarManager::kDitrackPt] = KFGeoTwoTracks.GetPt();
      }

      KFGeoFourProng.SetConstructMethod(2);
      KFGeoFourProng.AddDaughter(KFGeoTwoLeptons);
      KFGeoFourProng.AddDaughter(KFGeoTwoTracks);
    }

    if constexpr (candidateType == kPsi2StoJpsiPiPi) {
      KFPTrack kfpTrack0 = createKFPTrackFromTrack(lepton1);
      lepton1KF = KFParticle(kfpTrack0, -11 * lepton1.sign());
      KFPTrack kfpTrack1 = createKFPTrackFromTrack(lepton2);
      lepton2KF = KFParticle(kfpTrack1, -11 * lepton2.sign());
      KFPTrack kfpTrack2 = createKFPTrackFromTrack(track1);
      trk1KF = KFParticle(kfpTrack2, 211 * track1.sign());
      KFPTrack kfpTrack3 = createKFPTrackFromTrack(track2);
      trk2KF = KFParticle(kfpTrack3, 211 * track2.sign());

      KFGeoTwoLeptons.SetConstructMethod(2);
      KFGeoTwoLeptons.AddDaughter(lepton1KF);
      KFGeoTwoLeptons.AddDaughter(lepton2KF);

      if (fgUsedVars[kPairMass] || fgUsedVars[kPairPt]) {
        values[VarManager::kPairMass] = KFGeoTwoLeptons.GetMass();
        values[VarManager::kPairPt] = KFGeoTwoLeptons.GetPt();
      }

      KFGeoFourProng.SetConstructMethod(2);
      KFGeoFourProng.AddDaughter(KFGeoTwoLeptons);
      KFGeoFourProng.AddDaughter(trk1KF);
      KFGeoFourProng.AddDaughter(trk2KF);
    }

    if (fgUsedVars[kKFMass]) {
      float mass = 0., massErr = 0.;
      if (!KFGeoFourProng.GetMass(mass, massErr))
        values[kKFMass] = mass;
      else
        values[kKFMass] = -999.;
    }

    KFPVertex kfpVertex = createKFPVertexFromCollision(collision);
    values[kKFNContributorsPV] = kfpVertex.GetNContributors();
    KFParticle KFPV(kfpVertex);
    double dxQuadlet2PV = KFGeoFourProng.GetX() - KFPV.GetX();
    double dyQuadlet2PV = KFGeoFourProng.GetY() - KFPV.GetY();
    double dzQuadlet2PV = KFGeoFourProng.GetZ() - KFPV.GetZ();

    values[kVertexingLxy] = std::sqrt(dxQuadlet2PV * dxQuadlet2PV + dyQuadlet2PV * dyQuadlet2PV);
    values[kVertexingLz] = std::sqrt(dzQuadlet2PV * dzQuadlet2PV);
    values[kVertexingLxyz] = std::sqrt(dxQuadlet2PV * dxQuadlet2PV + dyQuadlet2PV * dyQuadlet2PV + dzQuadlet2PV * dzQuadlet2PV);

    values[kVertexingLxyErr] = (KFPV.GetCovariance(0) + KFGeoFourProng.GetCovariance(0)) * dxQuadlet2PV * dxQuadlet2PV + (KFPV.GetCovariance(2) + KFGeoFourProng.GetCovariance(2)) * dyQuadlet2PV * dyQuadlet2PV + 2 * ((KFPV.GetCovariance(1) + KFGeoFourProng.GetCovariance(1)) * dxQuadlet2PV * dyQuadlet2PV);
    values[kVertexingLzErr] = (KFPV.GetCovariance(5) + KFGeoFourProng.GetCovariance(5)) * dzQuadlet2PV * dzQuadlet2PV;
    values[kVertexingLxyzErr] = (KFPV.GetCovariance(0) + KFGeoFourProng.GetCovariance(0)) * dxQuadlet2PV * dxQuadlet2PV + (KFPV.GetCovariance(2) + KFGeoFourProng.GetCovariance(2)) * dyQuadlet2PV * dyQuadlet2PV + (KFPV.GetCovariance(5) + KFGeoFourProng.GetCovariance(5)) * dzQuadlet2PV * dzQuadlet2PV + 2 * ((KFPV.GetCovariance(1) + KFGeoFourProng.GetCovariance(1)) * dxQuadlet2PV * dyQuadlet2PV + (KFPV.GetCovariance(3) + KFGeoFourProng.GetCovariance(3)) * dxQuadlet2PV * dzQuadlet2PV + (KFPV.GetCovariance(4) + KFGeoFourProng.GetCovariance(4)) * dyQuadlet2PV * dzQuadlet2PV);

    if (fabs(values[kVertexingLxy]) < 1.e-8f)
      values[kVertexingLxy] = 1.e-8f;
    values[kVertexingLxyErr] = values[kVertexingLxyErr] < 0. ? 1.e8f : std::sqrt(values[kVertexingLxyErr]) / values[kVertexingLxy];
    if (fabs(values[kVertexingLz]) < 1.e-8f)
      values[kVertexingLz] = 1.e-8f;
    values[kVertexingLzErr] = values[kVertexingLzErr] < 0. ? 1.e8f : std::sqrt(values[kVertexingLzErr]) / values[kVertexingLz];
    if (fabs(values[kVertexingLxyz]) < 1.e-8f)
      values[kVertexingLxyz] = 1.e-8f;
    values[kVertexingLxyzErr] = values[kVertexingLxyzErr] < 0. ? 1.e8f : std::sqrt(values[kVertexingLxyzErr]) / values[kVertexingLxyz];

    values[kVertexingTauxy] = KFGeoFourProng.GetPseudoProperDecayTime(KFPV, KFGeoFourProng.GetMass()) / (o2::constants::physics::LightSpeedCm2NS);
    values[kVertexingTauz] = -1 * dzQuadlet2PV * KFGeoFourProng.GetMass() / (TMath::Abs(KFGeoFourProng.GetPz()) * o2::constants::physics::LightSpeedCm2NS);
    values[kVertexingPz] = TMath::Abs(KFGeoFourProng.GetPz());
    values[kVertexingSV] = KFGeoFourProng.GetZ();

    values[kVertexingTauxyErr] = values[kVertexingLxyErr] * KFGeoFourProng.GetMass() / (KFGeoFourProng.GetPt() * o2::constants::physics::LightSpeedCm2NS);
    values[kVertexingTauzErr] = values[kVertexingLzErr] * KFGeoFourProng.GetMass() / (TMath::Abs(KFGeoFourProng.GetPz()) * o2::constants::physics::LightSpeedCm2NS);
    values[kVertexingChi2PCA] = KFGeoFourProng.GetChi2();
    values[kCosPointingAngle] = (std::sqrt(dxQuadlet2PV * dxQuadlet2PV) * v1234.Px() +
                                 std::sqrt(dyQuadlet2PV * dyQuadlet2PV) * v1234.Py() +
                                 std::sqrt(dzQuadlet2PV * dzQuadlet2PV) * v1234.Pz()) /
                                (v1234.P() * values[VarManager::kVertexingLxyz]);
    // // run 2 definitions: Decay length projected onto the momentum vector of the candidate
    values[kVertexingLzProjected] = (dzQuadlet2PV * KFGeoFourProng.GetPz()) / TMath::Sqrt(KFGeoFourProng.GetPz() * KFGeoFourProng.GetPz());
    values[kVertexingLxyProjected] = (dxQuadlet2PV * KFGeoFourProng.GetPx()) + (dyQuadlet2PV * KFGeoFourProng.GetPy());
    values[kVertexingLxyProjected] = values[kVertexingLxyProjected] / TMath::Sqrt((KFGeoFourProng.GetPx() * KFGeoFourProng.GetPx()) + (KFGeoFourProng.GetPy() * KFGeoFourProng.GetPy()));
    values[kVertexingLxyzProjected] = (dxQuadlet2PV * KFGeoFourProng.GetPx()) + (dyQuadlet2PV * KFGeoFourProng.GetPy()) + (dzQuadlet2PV * KFGeoFourProng.GetPz());
    values[kVertexingLxyzProjected] = values[kVertexingLxyzProjected] / TMath::Sqrt((KFGeoFourProng.GetPx() * KFGeoFourProng.GetPx()) + (KFGeoFourProng.GetPy() * KFGeoFourProng.GetPy()) + (KFGeoFourProng.GetPz() * KFGeoFourProng.GetPz()));
    values[kVertexingTauxyProjected] = values[kVertexingLxyProjected] * KFGeoFourProng.GetMass() / (KFGeoFourProng.GetPt());
    values[kVertexingTauxyProjectedNs] = values[kVertexingTauxyProjected] / o2::constants::physics::LightSpeedCm2NS;
    values[kVertexingTauzProjected] = values[kVertexingLzProjected] * KFGeoFourProng.GetMass() / TMath::Abs(KFGeoFourProng.GetPz());
    values[kKFChi2OverNDFGeo] = KFGeoFourProng.GetChi2() / KFGeoFourProng.GetNDF();
  } else {
    return;
  }
}

//__________________________________________________________________
template <int candidateType, typename T1, typename T2>
void VarManager::FillQuadMC(T1 const& dilepton, T2 const& track1, T2 const& track2, float* values)
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

  ROOT::Math::PtEtaPhiMVector v1(dilepton.pt(), dilepton.eta(), dilepton.phi(), defaultDileptonMass);
  ROOT::Math::PtEtaPhiMVector v2(track1.pt(), track1.eta(), track1.phi(), hadronMass1);
  ROOT::Math::PtEtaPhiMVector v3(track2.pt(), track2.eta(), track2.phi(), hadronMass2);
  ROOT::Math::PtEtaPhiMVector v123 = v1 + v2 + v3;
  ROOT::Math::PtEtaPhiMVector v23 = v2 + v3;
  values[kQuadMass] = v123.M();
  values[kQuadDefaultDileptonMass] = v123.M();
  values[kQuadPt] = v123.Pt();
  values[kQuadEta] = v123.Eta();
  values[kQuadPhi] = v123.Phi();
  values[kQ] = v123.M() - defaultDileptonMass - v23.M();
  values[kDeltaR1] = ROOT::Math::VectorUtil::DeltaR(v1, v2);
  values[kDeltaR2] = ROOT::Math::VectorUtil::DeltaR(v1, v3);
  values[kDeltaR] = sqrt(pow(values[kDeltaR1], 2) + pow(values[kDeltaR2], 2));
  values[kDitrackMass] = v23.M();
  values[kDitrackPt] = v23.Pt();
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

/// Fill BDT score values.
/// Supports binary (1 output) and multiclass (3 outputs) models.
template <typename T1>
void VarManager::FillBdtScore(T1 const& bdtScore, float* values)
{
  if (!values) {
    values = fgValues;
  }

  if (bdtScore.size() == 1) {
    values[kBdtBackground] = bdtScore[0];
  } else if (bdtScore.size() == 3) {
    values[kBdtBackground] = bdtScore[0];
    values[kBdtPrompt] = bdtScore[1];
    values[kBdtNonprompt] = bdtScore[2];
  } else {
    LOG(warning) << "Unexpected number of BDT outputs: " << bdtScore.size();
  }
}
//__________________________________________________________________
template <typename T1, typename T2>
float VarManager::LorentzTransformJpsihadroncosChi(TString Option, T1 const& v1, T2 const& v2)
{
  float value = -999.0f;
  auto beta_v1 = v1.BoostToCM();
  ROOT::Math::Boost boostv1{beta_v1};
  auto v2_boost = boostv1(v2);
  auto p_v1_lab = v1.Vect();
  float p1_lab = p_v1_lab.R();
  if (Option == "coschi") {
    auto p_v2_boost = v2_boost.Vect();
    float p_boost = p_v2_boost.R();
    value = p_v2_boost.Dot(p_v1_lab) / (p1_lab * p_boost);
  } else if (Option == "costheta") {
    auto p_v2_lab = v2.Vect();
    float p2_lab = p_v2_lab.R();
    value = p_v1_lab.Dot(p_v2_lab) / (p1_lab * p2_lab);
  } else if (Option == "weight_boost") {
    value = v2_boost.E();
  }
  return value;
}
template <typename T1, typename T2, typename T3, typename T4, typename T5>
void VarManager::FillFIT(T1 const& bc, T2 const& bcs, T3 const& ft0s, T4 const& fv0as, T5 const& fdds, float* values)
{
  if (!values) {
    values = fgValues;
  }

  // Initialize FIT info structure
  upchelpers::FITInfo fitInfo{};
  udhelpers::getFITinfo(fitInfo, bc, bcs, ft0s, fv0as, fdds);

  // Fill FT0 information
  values[kAmplitudeFT0A] = fitInfo.ampFT0A;
  values[kAmplitudeFT0C] = fitInfo.ampFT0C;
  values[kTimeFT0A] = fitInfo.timeFT0A;
  values[kTimeFT0C] = fitInfo.timeFT0C;
  values[kTriggerMaskFT0] = static_cast<float>(fitInfo.triggerMaskFT0);
  const auto ft0Index = bc.ft0Id();
  if (ft0Index < 0 || ft0Index >= ft0s.size()) {
    values[kNFiredChannelsFT0A] = -1;
    values[kNFiredChannelsFT0C] = -1;
  } else {
    const auto ft0 = ft0s.iteratorAt(ft0Index);
    values[kNFiredChannelsFT0A] = ft0.channelA().size();
    values[kNFiredChannelsFT0C] = ft0.channelC().size();
  }
  // Fill FDD information
  values[kAmplitudeFDDA] = fitInfo.ampFDDA;
  values[kAmplitudeFDDC] = fitInfo.ampFDDC;
  values[kTimeFDDA] = fitInfo.timeFDDA;
  values[kTimeFDDC] = fitInfo.timeFDDC;
  values[kTriggerMaskFDD] = static_cast<float>(fitInfo.triggerMaskFDD);

  // Fill FV0A information
  values[kAmplitudeFV0A] = fitInfo.ampFV0A;
  values[kTimeFV0A] = fitInfo.timeFV0A;
  values[kTriggerMaskFV0A] = static_cast<float>(fitInfo.triggerMaskFV0A);
  const auto fv0aIndex = bc.fv0aId();
  if (fv0aIndex < 0 || fv0aIndex >= fv0as.size()) {
    values[kNFiredChannelsFV0A] = -1;
  } else {
    const auto fv0a = fv0as.iteratorAt(fv0aIndex);
    values[kNFiredChannelsFV0A] = fv0a.channel().size();
  }
  // Fill pileup flags
  values[kBBFT0Apf] = static_cast<float>(fitInfo.BBFT0Apf);
  values[kBGFT0Apf] = static_cast<float>(fitInfo.BGFT0Apf);
  values[kBBFT0Cpf] = static_cast<float>(fitInfo.BBFT0Cpf);
  values[kBGFT0Cpf] = static_cast<float>(fitInfo.BGFT0Cpf);
  values[kBBFV0Apf] = static_cast<float>(fitInfo.BBFV0Apf);
  values[kBGFV0Apf] = static_cast<float>(fitInfo.BGFV0Apf);
  values[kBBFDDApf] = static_cast<float>(fitInfo.BBFDDApf);
  values[kBGFDDApf] = static_cast<float>(fitInfo.BGFDDApf);
  values[kBBFDDCpf] = static_cast<float>(fitInfo.BBFDDCpf);
  values[kBGFDDCpf] = static_cast<float>(fitInfo.BGFDDCpf);
}


template <uint32_t fillMap, typename T>
void VarManager::FillEventAlice3(T const& event, float* values)
{
  if (!values) {
    values = fgValues;
  }
  if constexpr ((fillMap & CollisionTimestamp) > 0) {
    values[kTimestamp] = event.timestamp();
  }

  if constexpr ((fillMap & Collision) > 0) {
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

  if constexpr ((fillMap & ReducedEvent) > 0) {
    values[kRunNo] = -1;
    values[kVtxX] = event.posX();
    values[kVtxY] = event.posY();
    values[kVtxZ] = event.posZ();
    values[kVtxNcontrib] = event.numContrib();
    values[kCollisionTime] = event.collisionTime();
    values[kCollisionTimeRes] = event.collisionTimeRes();
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
}

template <uint32_t fillMap, typename T>
void VarManager::FillTrackAlice3(T const& track, float* values)
{
  if (!values) {
    values = fgValues;
  }

  if constexpr ((fillMap & Track) > 0 || (fillMap & ReducedTrack) > 0) {
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

    if (fgUsedVars[kPVContributor]) {
      values[kPVContributor] = (track.flags() & o2::aod::track::PVContributor) > 0;
    }

    if (fgUsedVars[kITSClusterMap]) {
      values[kITSClusterMap] = track.itsClusterMap();
    }

    values[kITSchi2] = track.itsChi2NCl();
  }

  // Quantities based on the barrel tables
  if constexpr ((fillMap & TrackExtra) > 0 || (fillMap & ReducedTrackBarrel) > 0) {

    values[kTrackLength] = track.length();

    if constexpr ((fillMap & ReducedTrackBarrel) > 0) {
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

  // Quantities based on the barrel covariance tables
  if constexpr ((fillMap & TrackCov) > 0 || (fillMap & ReducedTrackBarrelCov) > 0) {
    values[kTrackCYY] = track.cYY();
    values[kTrackCZZ] = track.cZZ();
    values[kTrackCSnpSnp] = track.cSnpSnp();
    values[kTrackCTglTgl] = track.cTglTgl();
    values[kTrackC1Pt21Pt2] = track.c1Pt21Pt2();
  }

  if constexpr ((fillMap & TrackPID) > 0 || (fillMap & ReducedTrackBarrelPID) > 0) {
    
    values[kOTTOTSignal] = track.timeOverThresholdBarrel();
    values[kOTnSigmaEl]  = track.nSigmaTrkEl();
    values[kOTnSigmaMu]  = track.nSigmaTrkMu();
    values[kOTnSigmaPi]  = track.nSigmaTrkPi();
    values[kOTnSigmaKa]  = track.nSigmaTrkKa();
    values[kOTnSigmaPr]  = track.nSigmaTrkPr();
    values[kOTnSigmaDe]  = track.nSigmaTrkDe();
    values[kOTnSigmaTr]  = track.nSigmaTrkTr();
    values[kOTnSigmaHe3] = track.nSigmaTrkHe();
    values[kOTnSigmaAl]  = track.nSigmaTrkAl();
    values[kHasRICHSig]  = track.hasSig();
    values[kHasRICHSigInGas] = track.hasSigInGas();
    values[kHasRICHSigEl] = track.hasSigEl();
    values[kHasRICHSigMu] = track.hasSigMu();
    values[kHasRICHSigPi] = track.hasSigPi();
    values[kHasRICHSigKa] = track.hasSigKa();
    values[kHasRICHSigPr] = track.hasSigPr();
    values[kHasRICHSigDe] = track.hasSigDe();
    values[kHasRICHSigTr] = track.hasSigTr();
    values[kHasRICHSigHe3] = track.hasSigHe3();
    values[kHasRICHSigAl] = track.hasSigAl();
    values[kRICHnSigmaEl]  = track.nSigmaElectronRich();
    values[kRICHnSigmaMu]  = track.nSigmaMuonRich();
    values[kRICHnSigmaPi]  = track.nSigmaPionRich();
    values[kRICHnSigmaKa]  = track.nSigmaKaonRich();
    values[kRICHnSigmaPr]  = track.nSigmaProtonRich();
    values[kRICHnSigmaDe]  = track.nSigmaDeuteronRich();
    values[kRICHnSigmaTr]  = track.nSigmaTritonRich();
    values[kRICHnSigmaHe3] = track.nSigmaHelium3Rich();
    values[kRICHnSigmaAl]  = track.nSigmaAlphaRich();
    values[kTOFEventTime]  = track.tofEventTime();
    values[kTOFEventTimeErr] = track.tofEventTimeErr();
    values[kOuterTOFnSigmaEl]  = track.nSigmaElectronOuterTOF();
    values[kOuterTOFnSigmaMu]  = track.nSigmaMuonOuterTOF();
    values[kOuterTOFnSigmaPi]  = track.nSigmaPionOuterTOF();
    values[kOuterTOFnSigmaKa]  = track.nSigmaKaonOuterTOF();
    values[kOuterTOFnSigmaPr]  = track.nSigmaProtonOuterTOF();
    values[kOuterTOFnSigmaDe]  = track.nSigmaDeuteronOuterTOF();
    values[kOuterTOFnSigmaTr]  = track.nSigmaTritonOuterTOF();
    values[kOuterTOFnSigmaHe3] = track.nSigmaHelium3OuterTOF();
    values[kOuterTOFnSigmaAl]  = track.nSigmaAlphaOuterTOF();
    values[kInnerTOFnSigmaEl]  = track.nSigmaElectronInnerTOF();
    values[kInnerTOFnSigmaMu]  = track.nSigmaMuonInnerTOF();
    values[kInnerTOFnSigmaPi]  = track.nSigmaPionInnerTOF();
    values[kInnerTOFnSigmaKa]  = track.nSigmaKaonInnerTOF();
    values[kInnerTOFnSigmaPr]  = track.nSigmaProtonInnerTOF();
    values[kInnerTOFnSigmaDe]  = track.nSigmaDeuteronInnerTOF();
    values[kInnerTOFnSigmaTr]  = track.nSigmaTritonInnerTOF();
    values[kInnerTOFnSigmaHe3] = track.nSigmaHelium3InnerTOF();
    values[kInnerTOFnSigmaAl]  = track.nSigmaAlphaInnerTOF();
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


template <int pairType, uint32_t fillMap, typename T1, typename T2>
void VarManager::FillPairAlice3(T1 const& t1, T2 const& t2, float* values)
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

  values[kCharge] = t1.sign() + t2.sign();
  values[kCharge1] = t1.sign();
  values[kCharge2] = t2.sign();
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

  values[kPt1] = t1.pt();
  values[kEta1] = t1.eta();
  values[kPhi1] = t1.phi();
  values[kPt2] = t2.pt();
  values[kEta2] = t2.eta();
  values[kPhi2] = t2.phi();

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

  if (fgUsedVars[kDeltaEtaPair2]) {
    values[kDeltaEtaPair2] = v1.Eta() - v2.Eta();
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

  // polarization parameters
  bool useHE = fgUsedVars[kCosThetaHE] || fgUsedVars[kPhiHE]; // helicity frame
  bool useCS = fgUsedVars[kCosThetaCS] || fgUsedVars[kPhiCS]; // Collins-Soper frame
  bool usePP = fgUsedVars[kCosThetaPP];                       // production plane frame
  bool useRM = fgUsedVars[kCosThetaRM];                       // Random frame

  if (useHE || useCS || usePP || useRM) {
    ROOT::Math::Boost boostv12{v12.BoostToCM()};
    ROOT::Math::XYZVectorF v1_CM{(boostv12(v1).Vect()).Unit()};
    ROOT::Math::XYZVectorF v2_CM{(boostv12(v2).Vect()).Unit()};
    ROOT::Math::XYZVectorF Beam1_CM{(boostv12(fgBeamA).Vect()).Unit()};
    ROOT::Math::XYZVectorF Beam2_CM{(boostv12(fgBeamC).Vect()).Unit()};

    // using positive sign convention for the first track
    ROOT::Math::XYZVectorF v_CM = (t1.sign() > 0 ? v1_CM : v2_CM);

    if (useHE) {
      ROOT::Math::XYZVectorF zaxis_HE{(v12.Vect()).Unit()};
      ROOT::Math::XYZVectorF yaxis_HE{(Beam1_CM.Cross(Beam2_CM)).Unit()};
      ROOT::Math::XYZVectorF xaxis_HE{(yaxis_HE.Cross(zaxis_HE)).Unit()};
      if (fgUsedVars[kCosThetaHE])
        values[kCosThetaHE] = zaxis_HE.Dot(v_CM);
      if (fgUsedVars[kPhiHE]) {
        values[kPhiHE] = TMath::ATan2(yaxis_HE.Dot(v_CM), xaxis_HE.Dot(v_CM));
        if (values[kPhiHE] < 0) {
          values[kPhiHE] += 2 * TMath::Pi(); // ensure phi is in [0, 2pi]
        }
      }
      if (fgUsedVars[kPhiTildeHE]) {
        if (fgUsedVars[kCosThetaHE] && fgUsedVars[kPhiHE]) {
          if (values[kCosThetaHE] > 0) {
            values[kPhiTildeHE] = values[kPhiHE] - 0.25 * TMath::Pi(); // phi_tilde = phi - pi/4
            if (values[kPhiTildeHE] < 0) {
              values[kPhiTildeHE] += 2 * TMath::Pi(); // ensure phi_tilde is in [0, 2pi]
            }
          } else {
            values[kPhiTildeHE] = values[kPhiHE] - 0.75 * TMath::Pi(); // phi_tilde = phi - 3pi/4
            if (values[kPhiTildeHE] < 0) {
              values[kPhiTildeHE] += 2 * TMath::Pi(); // ensure phi_tilde is in [0, 2pi]
            }
          }
        } else {
          values[kPhiTildeHE] = -999; // not computable
        }
      }
    }

    if (useCS) {
      ROOT::Math::XYZVectorF zaxis_CS{(Beam1_CM - Beam2_CM).Unit()};
      ROOT::Math::XYZVectorF yaxis_CS{(Beam1_CM.Cross(Beam2_CM)).Unit()};
      ROOT::Math::XYZVectorF xaxis_CS{(yaxis_CS.Cross(zaxis_CS)).Unit()};
      if (fgUsedVars[kCosThetaCS])
        values[kCosThetaCS] = zaxis_CS.Dot(v_CM);
      if (fgUsedVars[kPhiCS]) {
        values[kPhiCS] = TMath::ATan2(yaxis_CS.Dot(v_CM), xaxis_CS.Dot(v_CM));
        if (values[kPhiCS] < 0) {
          values[kPhiCS] += 2 * TMath::Pi(); // ensure phi is in [0, 2pi]
        }
      }
      if (fgUsedVars[kPhiTildeCS]) {
        if (fgUsedVars[kCosThetaCS] && fgUsedVars[kPhiCS]) {
          if (values[kCosThetaCS] > 0) {
            values[kPhiTildeCS] = values[kPhiCS] - 0.25 * TMath::Pi(); // phi_tilde = phi - pi/4
            if (values[kPhiTildeCS] < 0) {
              values[kPhiTildeCS] += 2 * TMath::Pi(); // ensure phi_tilde is in [0, 2pi]
            }
          } else {
            values[kPhiTildeCS] = values[kPhiCS] - 0.75 * TMath::Pi(); // phi_tilde = phi - 3pi/4
            if (values[kPhiTildeCS] < 0) {
              values[kPhiTildeCS] += 2 * TMath::Pi(); // ensure phi_tilde is in [0, 2pi]
            }
          }
        } else {
          values[kPhiTildeCS] = -999; // not computable
        }
      }
    }

    if (usePP) {
      ROOT::Math::XYZVector zaxis_PP = ROOT::Math::XYZVector(v12.Py(), -v12.Px(), 0.f);
      ROOT::Math::XYZVector yaxis_PP{(v12.Vect()).Unit()};
      ROOT::Math::XYZVector xaxis_PP{(yaxis_PP.Cross(zaxis_PP)).Unit()};
      if (fgUsedVars[kCosThetaPP]) {
        values[kCosThetaPP] = zaxis_PP.Dot(v_CM) / std::sqrt(zaxis_PP.Mag2());
      }
      if (fgUsedVars[kPhiPP]) {
        values[kPhiPP] = TMath::ATan2(yaxis_PP.Dot(v_CM), xaxis_PP.Dot(v_CM));
        if (values[kPhiPP] < 0) {
          values[kPhiPP] += 2 * TMath::Pi(); // ensure phi is in [0, 2pi]
        }
      }
      if (fgUsedVars[kPhiTildePP]) {
        if (fgUsedVars[kCosThetaPP] && fgUsedVars[kPhiPP]) {
          if (values[kCosThetaPP] > 0) {
            values[kPhiTildePP] = values[kPhiPP] - 0.25 * TMath::Pi(); // phi_tilde = phi - pi/4
            if (values[kPhiTildePP] < 0) {
              values[kPhiTildePP] += 2 * TMath::Pi(); // ensure phi_tilde is in [0, 2pi]
            }
          } else {
            values[kPhiTildePP] = values[kPhiPP] - 0.75 * TMath::Pi(); // phi_tilde = phi - 3pi/4
            if (values[kPhiTildePP] < 0) {
              values[kPhiTildePP] += 2 * TMath::Pi(); // ensure phi_tilde is in [0, 2pi]
            }
          }
        } else {
          values[kPhiTildePP] = -999; // not computable
        }
      }
    }

    if (useRM) {
      double randomCostheta = gRandom->Uniform(-1., 1.);
      double randomPhi = gRandom->Uniform(0., 2. * TMath::Pi());
      ROOT::Math::XYZVectorF zaxis_RM(randomCostheta, std::sqrt(1 - randomCostheta * randomCostheta) * std::cos(randomPhi), std::sqrt(1 - randomCostheta * randomCostheta) * std::sin(randomPhi));
      if (fgUsedVars[kCosThetaRM])
        values[kCosThetaRM] = zaxis_RM.Dot(v_CM);
    }
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

#endif // PWGDQ_CORE_VARMANAGER_H_
