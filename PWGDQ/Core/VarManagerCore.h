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

#ifndef PWGDQ_CORE_VARMANAGERCORE_H_
#define PWGDQ_CORE_VARMANAGERCORE_H_

#ifndef HomogeneousField
#define HomogeneousField
#endif

#include <TObject.h>

#include <DCAFitter/DCAFitterN.h>
#include <DCAFitter/FwdDCAFitterN.h>
#include <DetectorsBase/MatLayerCylSet.h>
#include <GlobalTracking/MatchGlobalFwd.h>
#include <ReconstructionDataFormats/GlobalFwdTrack.h>
#include <ReconstructionDataFormats/TrackFwd.h>
#include <ReconstructionDataFormats/Vertex.h>
#include <DataFormatsParameters/GRPLHCIFData.h>

#include <KFPTrack.h>
#include <KFPVertex.h>
#include <KFParticle.h>

#include <Math/Vector4Dfwd.h>
#include <Math/MatrixRepresentationsStatic.h>
#include <Math/SVector.h>

#include <map>

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
    ParticleMC = BIT(29),
    MuonDca = BIT(30)
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
    kCollisionRandom, // random number generated per collision (if required, can be used to perform random selections at the collision level)
    kIsPhysicsSelection,
    kIsTVXTriggered,             // Is trigger TVX
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
    kIsTriggerZNAZNC,            // trigger ZNA && ZNC
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
    kDCAzBimodalityCoefficient,
    kDCAzMean,
    kDCAzRMS,
    kDCAzSkewness,
    kDCAzKurtosis,
    kDCAzBimodalityCoefficientBinned,
    kDCAzBimodalityCoefficientBinnedTrimmed1,
    kDCAzBimodalityCoefficientBinnedTrimmed2,
    kDCAzBimodalityCoefficientBinnedTrimmed3,
    kDCAzMeanBinnedTrimmed1,
    kDCAzMeanBinnedTrimmed2,
    kDCAzMeanBinnedTrimmed3,
    kDCAzRMSBinnedTrimmed1,
    kDCAzRMSBinnedTrimmed2,
    kDCAzRMSBinnedTrimmed3,
    kDCAzFracAbove100um,
    kDCAzFracAbove200um,
    kDCAzFracAbove500um,
    kDCAzFracAbove1mm,
    kDCAzFracAbove2mm,
    kDCAzFracAbove5mm,
    kDCAzFracAbove10mm,
    kDCAzNPeaks,
    kDCAzNPeaksTrimmed1,
    kDCAzNPeaksTrimmed2,
    kDCAzNPeaksTrimmed3,
    kMCEventGeneratorId,
    kMCEventSubGeneratorId,
    kMCVtxX,
    kMCVtxY,
    kMCVtxZ,
    kMCEventTime,
    kMCEventWeight,
    kMCEventImpParam,
    kMCEventCentrFT0C,
    kMCEventPlaneAngle,
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
    kIsNoGap,      // No rapidity gap
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
    kTwoR2SP1,       // Scalar product resolution of event1 for ME technique
    kTwoR2SP2,       // Scalar product resolution of event2 for ME technique
    kTwoR2EP1,       // Event plane resolution of event2 for ME technique
    kTwoR2EP2,       // Event plane resolution of event2 for ME technique
    kNPairsPerEvent, // number of pairs per event in same-event or mixed-event pairing
    kInteractionRate,

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
    kDCAxy1,
    kDCAz1,
    kITSclusterMap1,
    kTPCnSigmaEl1,
    kPin_leg1,
    kTPCnSigmaKa_leg1,
    kPt2,
    kEta2,
    kPhi2,
    kCharge2,
    kDCAxy2,
    kDCAz2,
    kITSclusterMap2,
    kTPCnSigmaEl2,

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
    kMCAccweight,
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
    kMCHadronpt_randomPhi_trans,
    kMCWeight_before,
    kMCEWeight_before,
    kMCCosChi_gen,
    kMCWeight_gen,
    kMCdeltaeta_gen,
    kMCCosChi_rec,
    kMCWeight_rec,
    kMCdeltaeta_rec,
    kMCCosChi_randomPhi_trans_rec,
    kMCWeight_randomPhi_trans_rec,
    kMCCosChi_randomPhi_trans_gen,
    kMCWeight_randomPhi_trans_gen,

    // MC mother particle variables
    kMCMotherPdgCode,

    // MC pair variables
    kMCPt1,
    kMCEta1,
    kMCP1,
    kMCPt2,
    kMCEta2,
    kMCP2,
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
    kMCCosThetaStar,

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
    kVertexingLxyProjectedRecalculatePV,
    kVertexingLxyzProjected,
    kMCVertexingLzProjected,
    kMCVertexingLxyProjected,
    kMCVertexingLxyzProjected,
    kVertexingTauzProjected,
    kVertexingTauxyProjected,
    kVertexingTauxyProjectedPoleJPsiMass,
    kVertexingTauxyProjectedPoleJPsiMassRecalculatePV,
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
    kAbsCosThetaStarFT0C,
    kCos2ThetaStarFT0C,
    kCosThetaStarRandom,
    kCos2ThetaStarRandom,
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
    kU2Q2POS,
    kU2Q2NEG,
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
    kRandomPsi2,
    kCos2DeltaPhi,
    kCos2DeltaPhiPOS,
    kCos2DeltaPhiNEG,
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
    kPairEfficiency,
    kPairWeight,
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
    kWeight,
    kECWeight,
    kPtDau,
    kCosTheta,
    kEWeight_before,
    kWeight_before,
    kCosChi_randomPhi_trans,
    kCosChi_randomPhi_toward,
    kCosChi_randomPhi_away,
    kWeight_randomPhi_trans,
    kWeight_randomPhi_toward,
    kWeight_randomPhi_away,
    kdeltaphi_randomPhi_trans,
    kdeltaphi_randomPhi_toward,
    kdeltaphi_randomPhi_away,
    kdileptonmass,
    kPtDau_randomPhi_trans,

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

    // Resolution variables
    kDeltaPt,
    kPtResolution,
    kEtaResolution,

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
    kAmplitudeFT0M,
    kTimeFT0A,
    kTimeFT0C,
    kTriggerMaskFT0,
    kFT0OrA,
    kFT0OrC,
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
    kiTOFBeta,
    koTOFBeta,
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

  enum EfficiencyType {
    kNone = 0,
    kPairPtCentFT0cCosThetaStarFT0c,
    kPairPtCentFT0cCosThetaStarRandom,
    // Add more efficiency types as needed
    kNEfficiencyTypes
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

  static void SetUseVariable(int var);
  static void SetUseVars(const bool* usedVars);
  static void SetUseVars(const std::vector<int> usedVars);
  static bool GetUsedVar(int var);

         // Flag to  set PV recalculation via KF
  static void SetPVrecalculationKF(const bool pvRecalKF);

         // Setup the collision system
  static void SetCollisionSystem(TString system, float energy);
  static void SetCollisionSystem(o2::parameters::GRPLHCIFData* grplhcif);

  static void SetMagneticField(float magField);

         // Setup plane position for MFT-MCH matching
  static void SetMatchingPlane(float z);

  static float GetMatchingPlane();

         // Set z shift for forward tracks
  static void SetZShift(float z);

         // Setup the 2 prong KFParticle
  static void SetupTwoProngKFParticle(float magField);
  // Setup magnetic field for muon propagation
  static void SetupMuonMagField();

         // Setup the 2 prong DCAFitterN
  static void SetupTwoProngDCAFitter(float magField, bool propagateToPCA, float maxR, float maxDZIni, float minParamChange, float minRelChi2Change, bool useAbsDCA);

         // Setup the 2 prong FwdDCAFitterN
  static void SetupTwoProngFwdDCAFitter(float magField, bool propagateToPCA, float maxR, float minParamChange, float minRelChi2Change, bool useAbsDCA);
  // Use MatLayerCylSet to correct MCS in fwdtrack propagation
  static void SetupMatLUTFwdDCAFitter(o2::base::MatLayerCylSet* m);
  // Use GeometryManager to correct MCS in fwdtrack propagation
  static void SetupTGeoFwdDCAFitter();
  // No material budget in fwdtrack propagation
  static void SetupFwdDCAFitterNoCorr();
  // Setup the 3 prong KFParticle
  static void SetupThreeProngKFParticle(float magField);

         // Setup the 3 prong DCAFitterN
  static void SetupThreeProngDCAFitter(float magField, bool propagateToPCA, float maxR, float /*maxDZIni*/, float minParamChange, float minRelChi2Change, bool useAbsDCA);

         // Setup the 4 prong KFParticle
  static void SetupFourProngKFParticle(float magField);

         // Setup the 4 prong DCAFitterN
  static void SetupFourProngDCAFitter(float magField, bool propagateToPCA, float maxR, float /*maxDZIni*/, float minParamChange, float minRelChi2Change, bool useAbsDCA);

  static auto getEventPlane(int harm, float qnxa, float qnya);

  static float getDeltaPsiInRange(float psi1, float psi2, float harmonic);
  template <typename T, typename T1>
  static o2::dataformats::VertexBase RecalculatePrimaryVertex(T const& track0, T const& track1, const T1& collision);

  static std::tuple<float, float, float, float, float> BimodalityCoefficientUnbinned(const std::vector<float>& data);
  static std::tuple<float, float, float, float, float, int> BimodalityCoefficientAndNPeaks(const std::vector<float>& data, float binWidth, int trim = 0, float min = -15.0, float max = 15.0);

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
  static void FillEventTracks(T const& tracks, float* values = nullptr);
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
  template <int candidateType, typename T1>
  static void FillTrackCollisionMC(T1 const& track, const std::array<double, 3>& collPos, float massHyp = -1., float* values = nullptr);
  template <uint32_t fillMap, typename T, typename C, typename M, typename P>
  static void FillTrackCollisionMatCorr(T const& track, C const& collision, M const& materialCorr, P const& propagator, float* values = nullptr);
  template <typename U, typename T>
  static void FillTrackMC(const U& mcStack, T const& track, float* values = nullptr);
  template <int pairType, typename T, typename T1>
  static void FillEnergyCorrelatorsMC(T const& track, T1 const& t1, float* values = nullptr, float Translow = 1. / 3, float Transhigh = 2. / 3, float Accweight = 1.0f);
  template <uint32_t fillMap, typename T1, typename T2, typename C>
  static void FillPairPropagateMuon(T1 const& muon1, T2 const& muon2, const C& collision, float* values = nullptr);
  template <uint32_t fillMap, typename T1, typename T2, typename C>
  static void FillGlobalMuonRefit(T1 const& muontrack, T2 const& mfttrack, const C& collision, float* values = nullptr);
  template <uint32_t MuonfillMap, uint32_t MFTfillMap, typename T1, typename T2, typename C, typename C2>
  static void FillGlobalMuonRefitCov(T1 const& muontrack, T2 const& mfttrack, const C& collision, C2 const& mftcov, float* values = nullptr);
  template <int pairType, uint32_t fillMap, typename T1, typename T2>
  static void FillPair(T1 const& t1, T2 const& t2, float* values = nullptr);
  template <int pairType, uint32_t fillMap, typename T1, typename T2>
  static void FillPairRotation(T1 const& t1, T2 const& t2, float* values = nullptr);
  template <int pairType, uint32_t fillMap, typename C, typename T1, typename T2>
  static void FillPairCollision(C const& collision, T1 const& t1, T2 const& t2, float* values = nullptr);
  template <int pairType, uint32_t fillMap, typename C, typename T1, typename T2, typename M, typename P>
  static void FillPairCollisionMatCorr(C const& collision, T1 const& t1, T2 const& t2, M const& materialCorr, P const& propagator, float* values = nullptr);
  template <typename T1, typename T2, typename T3>
  static void FillTriple(T1 const& t1, T2 const& t2, T3 const& t3, float* values = nullptr, PairCandidateType pairType = kTripleCandidateToEEPhoton);
  template <uint32_t fillMap, int pairType, typename T1, typename T2>
  static void FillPairME(T1 const& t1, T2 const& t2, float* values = nullptr);
  template <typename T>
  static void FillPairMEAcrossTFs(T const& t1, T const& t2, float* values = nullptr);
  template <int pairType, typename T1, typename T2>
  static void FillPairMC(T1 const& t1, T2 const& t2, float* values = nullptr);
  template <int candidateType, typename T1, typename T2, typename T3>
  static void FillTripleMC(T1 const& t1, T2 const& t2, T3 const& t3, float* values = nullptr);
  template <int candidateType, typename T1, typename T2>
  static void FillQuadMC(T1 const& t1, T2 const& t2, T2 const& t3, float* values = nullptr);
  template <int pairType, uint32_t collFillMap, uint32_t fillMap, typename C, typename T>
  static void FillPairVertexing(C const& collision, T const& t1, T const& t2, bool propToSV = false, float* values = nullptr);
  template <int pairType, uint32_t collFillMap, uint32_t fillMap, typename C, typename T>
  static void FillPairVertexingRecomputePV(C const& /*collision*/, T const& t1, T const& t2, o2::dataformats::VertexBase pvRefitted, float* values = nullptr);
  template <uint32_t collFillMap, uint32_t fillMap, typename C, typename T>
  static void FillTripletVertexing(C const& collision, T const& t1, T const& t2, T const& t3, PairCandidateType tripletType, float* values = nullptr);
  template <int candidateType, uint32_t collFillMap, uint32_t fillMap, typename C, typename T1>
  static void FillDileptonTrackVertexing(C const& collision, T1 const& lepton1, T1 const& lepton2, T1 const& track, float* values);
  template <typename T1, typename T2>
  static void FillDileptonHadron(T1 const& dilepton, T2 const& hadron, float* values = nullptr, float hadronMass = 0.0f);
  template <typename T1, typename T2, typename T3>
  static void FillEnergyCorrelatorTriple(T1 const& lepton1, T2 const& lepton2, T3 const& hadron, float* values = nullptr, float Translow = 1. / 3, float Transhigh = 2. / 3, bool applyFitMass = false, float sidebandMass = 0.0f, float weight = 1.0f);
  template <int pairType, typename T1, typename T2, typename T3, typename T4, typename T5>
  static void FillEnergyCorrelatorsUnfoldingTriple(T1 const& lepton1, T2 const& lepton2, T3 const& hadron, T4 const& track, T5 const& t1, float* values = nullptr, bool applyFitMass = false, float Effweight_rec = 1.f, float Accweight_gen = 1.f, float Translow = 1. / 3, float Transhigh = 2. / 3);
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
  static void FillTrackAlice3(T const& track, float* values = nullptr);
  template <typename M, typename T>
  static void FillResolutions(M const& mcTrack, T const& track, float* values = nullptr);

  static void SetCalibrationObject(CalibObjects calib, TObject* obj);

  static void SetCalibrationType(int type, bool useInterpolation = true);
  static double ComputePIDcalibration(int species, double nSigmaValue);

  static void SetEfficiencyObject(int type, TObject* obj);
  static void FillEfficiency(float* values = nullptr);
  static TObject* GetCalibrationObject(CalibObjects calib);
  static void SetTPCInterSectorBoundary(float boundarySize);
  static void SetITSROFBorderselection(int bias, int length, int marginLow, int marginHigh);

  static void SetSORandEOR(uint64_t sor, uint64_t eor);

 public:
  VarManager();
  ~VarManager() override;

  static float fgValues[kNVars]; // array holding all variables computed during analysis
  static void ResetValues(int startValue = 0, int endValue = kNVars, float* values = nullptr);

 private:
  static bool fgUsedVars[kNVars]; // holds flags for when the corresponding variable is needed (e.g., in the histogram manager, in cuts, mixing handler, etc.)
  static bool fgUsedKF;
  static bool fgPVrecalKF;
  static void SetVariableDependencies(); // toggle those variables on which other used variables might depend

  static float fgMagField;
  static float fgzMatching;
  static float fgzShiftFwd;
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

  static int fgEfficiencyType;      // type of efficiency correction to apply
  static TObject* fgEfficiencyHist; // histogram for efficiency correction

  VarManager& operator=(const VarManager& c);
  VarManager(const VarManager& c);

  ClassDef(VarManager, 6);
};

#endif // PWGDQ_CORE_VARMANAGERCORE_H_
