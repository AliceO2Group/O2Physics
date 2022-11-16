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

#include <vector>
#include <map>
#include <cmath>
#include <iostream>
#include <utility>

#include <TObject.h>
#include <TString.h>
#include "TRandom.h"
#include "Math/Vector4D.h"
#include "Math/Vector3D.h"
#include "Math/GenVector/Boost.h"

#include "Framework/DataTypes.h"
#include "ReconstructionDataFormats/Track.h"
#include "ReconstructionDataFormats/Vertex.h"
#include "DetectorsVertexing/DCAFitterN.h"
#include "Common/CCDB/TriggerAliases.h"
#include "ReconstructionDataFormats/DCA.h"

#include "Math/SMatrix.h"
#include "ReconstructionDataFormats/TrackFwd.h"
#include "DetectorsVertexing/FwdDCAFitterN.h"

using std::cout;
using std::endl;
using SMatrix55 = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
using SMatrix5 = ROOT::Math::SVector<double, 5>;
using Vec3D = ROOT::Math::SVector<double, 3>;

// TODO: create an array holding these constants for all needed particles or check for a place where these are already defined
static const float fgkElectronMass = 0.000511; // GeV
static const float fgkMuonMass = 0.105658;     // GeV
static const float fgkKaonMass = 0.493677;     // GeV

//_________________________________________________________________________
class VarManager : public TObject
{
 public:
  // map the information contained in the objects passed to the Fill functions
  enum ObjTypes {
    // NOTE: Elements containing "Reduced" in their name refer to skimmed data tables
    //       and the ones that don't refer to tables from the Framework data model
    BC = BIT(0),
    Collision = BIT(1),
    CollisionCent = BIT(2),
    CollisionTimestamp = BIT(3),
    ReducedEvent = BIT(4),
    ReducedEventExtended = BIT(5),
    ReducedEventVtxCov = BIT(6),
    CollisionMC = BIT(7),
    ReducedEventMC = BIT(8),
    ReducedEventQvector = BIT(9),
    Track = BIT(0),
    TrackCov = BIT(1),
    TrackExtra = BIT(2),
    TrackPID = BIT(3),
    TrackDCA = BIT(4),
    TrackSelection = BIT(5),
    TrackV0Bits = BIT(6),
    ReducedTrack = BIT(7),
    ReducedTrackBarrel = BIT(8),
    ReducedTrackBarrelCov = BIT(9),
    ReducedTrackBarrelPID = BIT(10),
    Muon = BIT(11),
    MuonCov = BIT(12),
    ReducedMuon = BIT(13),
    ReducedMuonExtra = BIT(14),
    ReducedMuonCov = BIT(15),
    ParticleMC = BIT(16),
    Pair = BIT(17), // TODO: check whether we really need the Pair member here
    AmbiTrack = BIT(18),
    AmbiMuon = BIT(19),
    DalitzBits = BIT(20)
  };

  enum PairCandidateType {
    // TODO: need to agree on a scheme to incorporate all various hypotheses (e.g. e - mu, jpsi - K+, Jpsi - pipi,...)
    kJpsiToEE = 0,   // J/psi        -> e+ e-
    kJpsiToMuMu,     // J/psi        -> mu+ mu-
    kElectronMuon,   // Electron - muon correlations
    kBcToThreeMuons, // Bc         -> mu+ mu- mu+
    kBtoJpsiEEK,     // B+             -> e+ e- K+
    kNMaxCandidateTypes
  };

 public:
  enum Variables {
    kNothing = -1,
    // Run wise variables
    kRunNo = 0,
    kRunId,
    kNRunWiseVariables,

    // Event wise variables
    kTimestamp,
    kBC,
    kIsPhysicsSelection,
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
    kVtxCovXX,
    kVtxCovXY,
    kVtxCovXZ,
    kVtxCovYY,
    kVtxCovYZ,
    kVtxCovZZ,
    kVtxChi2,
    kCentVZERO,
    kMCEventGeneratorId,
    kMCVtxX,
    kMCVtxY,
    kMCVtxZ,
    kMCEventTime,
    kMCEventWeight,
    kMCEventImpParam,
    kQ2X0A, // q-vector (e.g. from TPC) with x component (harmonic 2 and power 0), sub-event A
    kQ2Y0A, // q-vector (e.g. from TPC) with y component (harmonic 2 and power 0), sub-event A
    kQ2X0B,
    kQ2Y0B,
    kQ2X0C,
    kQ2Y0C,
    kMultA, // Multiplicity of the sub-event A
    kMultB,
    kMultC,
    kQ3X0A, // q-vector (e.g. from TPC) with x component (harmonic 2 and power 0), sub-event A
    kQ3Y0A, // q-vector (e.g. from TPC) with y component (harmonic 2 and power 0), sub-event A
    kQ3X0B,
    kQ3Y0B,
    kQ3X0C,
    kQ3Y0C,
    kR2SP,
    kR3SP,
    kR2EP,
    kR3EP,
    kNEventWiseVariables,

    // Basic track/muon/pair wise variables
    kPt,
    kEta,
    kPhi,
    kP,
    kPx,
    kPy,
    kPz,
    kRap,
    kMass,
    kCharge,
    kNBasicTrackVariables,

    // Barrel track variables
    kPin,
    kIsGlobalTrack,
    kIsGlobalTrackSDD,
    kIsITSrefit,
    kIsSPDany,
    kIsSPDfirst,
    kIsSPDboth,
    kITSncls,
    kITSchi2,
    kITSlayerHit,
    kIsTPCrefit,
    kTPCncls,
    kITSClusterMap,
    kTPCnclsCR,
    kTPCchi2,
    kTPCsignal,
    kTPCsignalRandomized,
    kTPCsignalRandomizedDelta,
    kTRDsignal,
    kTRDPattern,
    kTOFbeta,
    kTrackLength,
    kTrackDCAxy,
    kTrackDCAz,
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
    kIsLegFromGamma,
    kIsLegFromK0S,
    kIsLegFromLambda,
    kIsLegFromAntiLambda,
    kIsLegFromOmega,
    kIsProtonFromLambdaAndAntiLambda,
    kIsDalitzLeg,  // Up to 8 dalitz selections
    kNBarrelTrackVariables = kIsDalitzLeg + 8,

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
    kMuonCYY,
    kMuonCPhiPhi,
    kMuonCTglTgl,
    kMuonC1Pt21Pt2,
    kNMuonTrackVariables,
    kMuonTrackType,
    kMuonDCAx,
    kMuonDCAy,

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
    kVertexingTauz,
    kVertexingTauzErr,
    kVertexingProcCode,
    kVertexingChi2PCA,
    kCosThetaHE,
    kPsiPair,
    kDeltaPhiPair,
    kQuadDCAabsXY,
    kQuadDCAsigXY,
    kCosPointingAngle,
    kImpParXYJpsi,
    kImpParXYK,
    kDCATrackProd,
    kDCATrackVtxProd,
    kU2Q2,
    kU3Q3,
    kCos2DeltaPhi,
    kCos3DeltaPhi,
    kNPairVariables,

    // Candidate-track correlation variables
    kPairMass,
    kPairMassDau,
    kPairPt,
    kPairPtDau,
    kPairEta,
    kPairPhi,
    kPairPhiv,
    kDeltaEta,
    kDeltaPhi,
    kDeltaPhiSym,
    kNCorrelationVariables,

    // Index used to scan bit maps
    kBitMapIndex,

    kNVars
  }; // end of Variables enumeration

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
  static int GetNRuns()
  {
    return fgRunMap.size();
  }
  static TString GetRunStr()
  {
    return fgRunStr;
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
  }
  static auto getEventPlane(int harm, float qnxa, float qnya)
  {
    // Compute event plane angle from qn vector components for the sub-event A
    return (1.0 / harm) * TMath::ATan(qnya / qnxa);
  };

  template <uint32_t fillMap, typename T>
  static void FillEvent(T const& event, float* values = nullptr);
  template <uint32_t fillMap, typename T>
  static void FillTrack(T const& track, float* values = nullptr);
  template <int pairType, uint32_t fillMap, typename T1, typename T2>
  static void FillPair(T1 const& t1, T2 const& t2, float* values = nullptr);
  template <int pairType, typename T1, typename T2>
  static void FillPairME(T1 const& t1, T2 const& t2, float* values = nullptr);
  template <typename T1, typename T2>
  static void FillPairMC(T1 const& t1, T2 const& t2, float* values = nullptr, PairCandidateType pairType = kJpsiToEE);
  template <int pairType, uint32_t collFillMap, uint32_t fillMap, typename C, typename T>
  static void FillPairVertexing(C const& collision, T const& t1, T const& t2, float* values = nullptr);
  template <int candidateType, uint32_t collFillMap, uint32_t fillMap, typename C, typename T1>
  static void FillDileptonTrackVertexing(C const& collision, T1 const& lepton1, T1 const& lepton2, T1 const& track, float* values);
  template <typename T1, typename T2>
  static void FillDileptonHadron(T1 const& dilepton, T2 const& hadron, float* values = nullptr, float hadronMass = 0.0f);
  template <typename C, typename A>
  static void FillQVectorFromGFW(C const& collision, A const& compA2, A const& compB2, A const& compC2, A const& compA3, A const& compB3, A const& compC3, float normA = 1.0, float normB = 1.0, float normC = 1.0, float* values = nullptr);
  template <int pairType, typename T1, typename T2>
  static void FillPairVn(T1 const& t1, T2 const& t2, float* values = nullptr);

 public:
  VarManager();
  ~VarManager() override;

  static float fgValues[kNVars]; // array holding all variables computed during analysis
  static void ResetValues(int startValue = 0, int endValue = kNVars, float* values = nullptr);

 private:
  static bool fgUsedVars[kNVars];        // holds flags for when the corresponding variable is needed (e.g., in the histogram manager, in cuts, mixing handler, etc.)
  static void SetVariableDependencies(); // toggle those variables on which other used variables might depend

  static std::map<int, int> fgRunMap; // map of runs to be used in histogram axes
  static TString fgRunStr;            // semi-colon separated list of runs, to be used for histogram axis labels

  static void FillEventDerived(float* values = nullptr);
  static void FillTrackDerived(float* values = nullptr);
  static float GetTPCPostCalibMap(float pin, float eta, int particle_type, TString period);
  static TString GetRunPeriod(float runNumber);
  template <typename T, typename U, typename V>
  static auto getRotatedCovMatrixXX(const T& matrix, U phi, V theta);

  static o2::vertexing::DCAFitterN<2> fgFitterTwoProngBarrel;
  static o2::vertexing::DCAFitterN<3> fgFitterThreeProngBarrel;
  static o2::vertexing::FwdDCAFitterN<2> fgFitterTwoProngFwd;
  static o2::vertexing::FwdDCAFitterN<3> fgFitterThreeProngFwd;

  VarManager& operator=(const VarManager& c);
  VarManager(const VarManager& c);

  ClassDef(VarManager, 1)
};

template <typename T, typename U, typename V>
auto VarManager::getRotatedCovMatrixXX(const T& matrix, U phi, V theta)
{
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

template <uint32_t fillMap, typename T>
void VarManager::FillEvent(T const& event, float* values)
{
  if (!values) {
    values = fgValues;
  }

  if constexpr ((fillMap & BC) > 0) {
    values[kRunNo] = event.bc().runNumber(); // accessed via Collisions table
    values[kBC] = event.bc().globalBC();
  }

  if constexpr ((fillMap & CollisionTimestamp) > 0) {
    values[kTimestamp] = event.timestamp();
  }

  if constexpr ((fillMap & Collision) > 0) {
    // TODO: trigger info from the event selection requires a separate flag
    //       so that it can be switched off independently of the rest of Collision variables (e.g. if event selection is not available)
    if (fgUsedVars[kIsINT7]) {
      values[kIsINT7] = (event.alias()[kINT7] > 0);
    }
    if (fgUsedVars[kIsEMC7]) {
      values[kIsEMC7] = (event.alias()[kEMC7] > 0);
    }
    if (fgUsedVars[kIsINT7inMUON]) {
      values[kIsINT7inMUON] = (event.alias()[kINT7inMUON] > 0);
    }
    if (fgUsedVars[kIsMuonSingleLowPt7]) {
      values[kIsMuonSingleLowPt7] = (event.alias()[kMuonSingleLowPt7] > 0);
    }
    if (fgUsedVars[kIsMuonSingleHighPt7]) {
      values[kIsMuonSingleHighPt7] = (event.alias()[kMuonSingleHighPt7] > 0);
    }
    if (fgUsedVars[kIsMuonUnlikeLowPt7]) {
      values[kIsMuonUnlikeLowPt7] = (event.alias()[kMuonUnlikeLowPt7] > 0);
    }
    if (fgUsedVars[kIsMuonLikeLowPt7]) {
      values[kIsMuonLikeLowPt7] = (event.alias()[kMuonLikeLowPt7] > 0);
    }
    if (fgUsedVars[kIsCUP8]) {
      values[kIsCUP8] = (event.alias()[kCUP8] > 0);
    }
    if (fgUsedVars[kIsCUP9]) {
      values[kIsCUP9] = (event.alias()[kCUP9] > 0);
    }
    if (fgUsedVars[kIsMUP10]) {
      values[kIsMUP10] = (event.alias()[kMUP10] > 0);
    }
    if (fgUsedVars[kIsMUP11]) {
      values[kIsMUP11] = (event.alias()[kMUP11] > 0);
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
  }

  if constexpr ((fillMap & CollisionCent) > 0) {
    values[kCentVZERO] = event.centRun2V0M();
  }

  // TODO: need to add EvSels and Cents tables, etc. in case of the central data model

  if constexpr ((fillMap & ReducedEvent) > 0) {
    values[kRunNo] = event.runNumber();
    values[kVtxX] = event.posX();
    values[kVtxY] = event.posY();
    values[kVtxZ] = event.posZ();
    values[kVtxNcontrib] = event.numContrib();
  }

  if constexpr ((fillMap & ReducedEventExtended) > 0) {
    values[kBC] = event.globalBC();
    values[kTimestamp] = event.timestamp();
    values[kCentVZERO] = event.centRun2V0M();
    if (fgUsedVars[kIsINT7]) {
      values[kIsINT7] = (event.triggerAlias() & (uint32_t(1) << kINT7)) > 0;
    }
    if (fgUsedVars[kIsEMC7]) {
      values[kIsEMC7] = (event.triggerAlias() & (uint32_t(1) << kEMC7)) > 0;
    }
    if (fgUsedVars[kIsINT7inMUON]) {
      values[kIsINT7inMUON] = (event.triggerAlias() & (uint32_t(1) << kINT7inMUON)) > 0;
    }
    if (fgUsedVars[kIsMuonSingleLowPt7]) {
      values[kIsMuonSingleLowPt7] = (event.triggerAlias() & (uint32_t(1) << kMuonSingleLowPt7)) > 0;
    }
    if (fgUsedVars[kIsMuonSingleHighPt7]) {
      values[kIsMuonSingleHighPt7] = (event.triggerAlias() & (uint32_t(1) << kMuonSingleHighPt7)) > 0;
    }
    if (fgUsedVars[kIsMuonUnlikeLowPt7]) {
      values[kIsMuonUnlikeLowPt7] = (event.triggerAlias() & (uint32_t(1) << kMuonUnlikeLowPt7)) > 0;
    }
    if (fgUsedVars[kIsMuonLikeLowPt7]) {
      values[kIsMuonLikeLowPt7] = (event.triggerAlias() & (uint32_t(1) << kMuonLikeLowPt7)) > 0;
    }
    if (fgUsedVars[kIsCUP8]) {
      values[kIsCUP8] = (event.triggerAlias() & (uint32_t(1) << kCUP8)) > 0;
    }
    if (fgUsedVars[kIsCUP9]) {
      values[kIsCUP9] = (event.triggerAlias() & (uint32_t(1) << kCUP9)) > 0;
    }
    if (fgUsedVars[kIsMUP10]) {
      values[kIsMUP10] = (event.triggerAlias() & (uint32_t(1) << kMUP10)) > 0;
    }
    if (fgUsedVars[kIsMUP11]) {
      values[kIsMUP11] = (event.triggerAlias() & (uint32_t(1) << kMUP11)) > 0;
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
    values[kR2SP] = (event.q2x0b() * event.q2x0c() + event.q2y0b() * event.q2y0c());
    values[kR3SP] = (event.q3x0b() * event.q3x0c() + event.q3y0b() * event.q3y0c());
    if (event.q2y0b() * event.q2y0c() != 0.0) {
      values[kR2EP] = TMath::Cos(2 * (getEventPlane(2, event.q2x0b(), event.q2y0b()) - getEventPlane(2, event.q2x0c(), event.q2y0c())));
    }
    if (event.q3y0b() * event.q3y0c() != 0.0) {
      values[kR3EP] = TMath::Cos(3 * (getEventPlane(3, event.q3x0b(), event.q3y0b()) - getEventPlane(3, event.q3x0c(), event.q3y0c())));
    }
  }

  if constexpr ((fillMap & CollisionMC) > 0) {
    values[kMCEventGeneratorId] = event.generatorsID();
    values[kMCVtxX] = event.posX();
    values[kMCVtxY] = event.posY();
    values[kMCVtxZ] = event.posZ();
    values[kMCEventTime] = event.t();
    values[kMCEventWeight] = event.weight();
    values[kMCEventImpParam] = event.impactParameter();
  }

  if constexpr ((fillMap & ReducedEventMC) > 0) {
    values[kMCEventGeneratorId] = event.generatorsID();
    values[kMCVtxX] = event.mcPosX();
    values[kMCVtxY] = event.mcPosY();
    values[kMCVtxZ] = event.mcPosZ();
    values[kMCEventTime] = event.t();
    values[kMCEventWeight] = event.weight();
    values[kMCEventImpParam] = event.impactParameter();
  }

  FillEventDerived(values);
}

template <uint32_t fillMap, typename T>
void VarManager::FillTrack(T const& track, float* values)
{
  if (!values) {
    values = fgValues;
  }

  // Quantities based on the basic table (contains just kine information and filter bits)
  if constexpr ((fillMap & Track) > 0 || (fillMap & Muon) > 0 || (fillMap & ReducedTrack) > 0 || (fillMap & ReducedMuon) > 0) {
    values[kPt] = track.pt();
    if (fgUsedVars[kPx]) {
      values[kPx] = track.px();
    }
    if (fgUsedVars[kPy]) {
      values[kPy] = track.py();
    }
    if (fgUsedVars[kPz]) {
      values[kPz] = track.pz();
    }
    values[kEta] = track.eta();
    values[kPhi] = track.phi();
    values[kCharge] = track.sign();

    if constexpr ((fillMap & ReducedTrack) > 0 && !((fillMap & Pair) > 0)) {
      values[kIsGlobalTrack] = track.filteringFlags() & (uint64_t(1) << 0);
      values[kIsGlobalTrackSDD] = track.filteringFlags() & (uint64_t(1) << 1);

      values[kIsLegFromGamma] = static_cast<bool>(track.filteringFlags() & (uint64_t(1) << 2));
      values[kIsLegFromK0S] = static_cast<bool>(track.filteringFlags() & (uint64_t(1) << 3));
      values[kIsLegFromLambda] = static_cast<bool>(track.filteringFlags() & (uint64_t(1) << 4));
      values[kIsLegFromAntiLambda] = static_cast<bool>(track.filteringFlags() & (uint64_t(1) << 5));
      values[kIsLegFromOmega] = static_cast<bool>(track.filteringFlags() & (uint64_t(1) << 6));

      values[kIsProtonFromLambdaAndAntiLambda] = static_cast<bool>((values[kIsLegFromLambda] * track.sign() > 0) || (values[kIsLegFromAntiLambda] * (-track.sign()) > 0));
      
      for (int i = 0; i < 8; i++) {
        values[kIsDalitzLeg + i] = static_cast<bool>(track.filteringFlags() & (uint64_t(1) << (15 + i)));
      }

    }
  }

  // Quantities based on the barrel tables
  if constexpr ((fillMap & TrackExtra) > 0 || (fillMap & ReducedTrackBarrel) > 0) {
    values[kPin] = track.tpcInnerParam();
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
    values[kITSchi2] = track.itsChi2NCl();
    values[kTPCncls] = track.tpcNClsFound();
    values[kTPCchi2] = track.tpcChi2NCl();
    values[kTrackLength] = track.length();
    values[kTPCnclsCR] = track.tpcNClsCrossedRows();
    values[kTRDPattern] = track.trdPattern();

    if constexpr ((fillMap & TrackExtra) > 0) {
      if (fgUsedVars[kITSncls]) {
        values[kITSncls] = track.itsNCls(); // dynamic column
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
  if constexpr((fillMap & DalitzBits) > 0) {
    for (int i = 0; i < 8; i++) {
      values[kIsDalitzLeg + i] = bool(track.dalitzBits() & (uint8_t(1) << i));
    }
  }

  // Quantities based on the barrel PID tables
  if constexpr ((fillMap & TrackPID) > 0 || (fillMap & ReducedTrackBarrelPID) > 0) {
    values[kTPCnSigmaEl] = track.tpcNSigmaEl();
    values[kTPCnSigmaMu] = track.tpcNSigmaMu();
    values[kTPCnSigmaPi] = track.tpcNSigmaPi();
    values[kTPCnSigmaKa] = track.tpcNSigmaKa();
    values[kTPCnSigmaPr] = track.tpcNSigmaPr();
    values[kTOFnSigmaEl] = track.tofNSigmaEl();
    values[kTOFnSigmaMu] = track.tofNSigmaMu();
    values[kTOFnSigmaPi] = track.tofNSigmaPi();
    values[kTOFnSigmaKa] = track.tofNSigmaKa();
    values[kTOFnSigmaPr] = track.tofNSigmaPr();
    values[kTPCsignal] = track.tpcSignal();
    values[kTRDsignal] = track.trdSignal();
    values[kTOFbeta] = track.beta();
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
    if (fgUsedVars[kTPCnSigmaEl_Corr] || fgUsedVars[kTPCnSigmaPi_Corr] || fgUsedVars[kTPCnSigmaPr_Corr]) {
      values[kTPCnSigmaEl_Corr] = values[kTPCnSigmaEl] - GetTPCPostCalibMap(values[kPin], values[kEta], 0, GetRunPeriod(values[kRunNo]));
      values[kTPCnSigmaPi_Corr] = values[kTPCnSigmaPi] - GetTPCPostCalibMap(values[kPin], values[kEta], 1, GetRunPeriod(values[kRunNo]));
      values[kTPCnSigmaPr_Corr] = values[kTPCnSigmaPr] - GetTPCPostCalibMap(values[kPin], values[kEta], 2, GetRunPeriod(values[kRunNo]));
    }
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
  }
  // Quantities based on the muon covariance table
  if constexpr ((fillMap & ReducedMuonCov) > 0 || (fillMap & MuonCov) > 0) {
    values[kMuonCXX] = track.cXX();
    values[kMuonCYY] = track.cYY();
    values[kMuonCPhiPhi] = track.cPhiPhi();
    values[kMuonCTglTgl] = track.cTglTgl();
    values[kMuonC1Pt21Pt2] = track.c1Pt21Pt2();
  }

  // Quantities based on the pair table(s)
  if constexpr ((fillMap & Pair) > 0) {
    values[kMass] = track.mass();
  }

  if constexpr ((fillMap & ParticleMC) > 0) {
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
    values[kMCY] = track.y();
    values[kMCParticleGeneratorId] = track.producedByGenerator();
  }

  // Derived quantities which can be computed based on already filled variables
  FillTrackDerived(values);
}

template <int pairType, uint32_t fillMap, typename T1, typename T2>
void VarManager::FillPair(T1 const& t1, T2 const& t2, float* values)
{
  if (!values) {
    values = fgValues;
  }

  float m1 = fgkElectronMass;
  float m2 = fgkElectronMass;
  if constexpr (pairType == kJpsiToMuMu) {
    m1 = fgkMuonMass;
    m2 = fgkMuonMass;
  }

  if constexpr (pairType == kElectronMuon) {
    m2 = fgkMuonMass;
  }

  ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), m1);
  ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), m2);
  ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
  values[kMass] = v12.M();
  values[kPt] = v12.Pt();
  values[kEta] = v12.Eta();
  values[kPhi] = v12.Phi();
  values[kRap] = -v12.Rapidity();

  if (fgUsedVars[kPsiPair]) {
    values[kDeltaPhiPair] = v1.Phi() - v2.Phi();
    double xipair = TMath::ACos((v1.Px() * v2.Px() + v1.Py() * v2.Py() + v1.Pz() * v2.Pz()) / v1.P() / v2.P());
    values[kPsiPair] = TMath::ASin((v1.Theta() - v2.Theta()) / xipair);
  }

  // CosTheta Helicity calculation
  ROOT::Math::Boost boostv12{v12.BoostToCM()};
  ROOT::Math::XYZVectorF v1_CM{(boostv12(v1).Vect()).Unit()};
  ROOT::Math::XYZVectorF v2_CM{(boostv12(v2).Vect()).Unit()};
  ROOT::Math::XYZVectorF zaxis{(v12.Vect()).Unit()};

  if (fgUsedVars[kCosThetaHE]) {
    values[kCosThetaHE] = (t1.sign() > 0 ? zaxis.Dot(v1_CM) : zaxis.Dot(v2_CM));
  }

  if constexpr ((pairType == kJpsiToEE) && ((fillMap & TrackCov) > 0 || (fillMap & ReducedTrackBarrelCov) > 0)) {

    if (fgUsedVars[kQuadDCAabsXY] || fgUsedVars[kQuadDCAsigXY]) {
      // Quantities based on the barrel tables
      double dca1 = t1.dcaXY();
      double dca2 = t2.dcaXY();
      double dca1sig = dca1 / std::sqrt(t1.cYY());
      double dca2sig = dca2 / std::sqrt(t2.cYY());

      values[kQuadDCAabsXY] = std::sqrt((dca1 * dca1 + dca2 * dca2) / 2);
      values[kQuadDCAsigXY] = std::sqrt((dca1sig * dca1sig + dca2sig * dca2sig) / 2);
    }
  }
  if (fgUsedVars[kPairPhiv]) {
    // cos(phiv) = w*a /|w||a|
    // with w = u x v
    // and  a = u x z / |u x z|   , unit vector perpendicular to v12 and z-direction (magnetic field)
    // u = v12 / |v12|            , the unit vector of v12
    // v = v1 x v2 / |v1 x v2|    , unit vector perpendicular to v1 and v2

    float bz = fgFitterTwoProngBarrel.getBz();

    // ordering of tracks, so v1 has larger momentum
    if (v1.P() < v2.P()) {
      ROOT::Math::PtEtaPhiMVector v3 = v1;
      v1 = v2;
      v2 = v3;
    }

    // momentum of e+ and e- in (ax,ay,az) axis. Note that az=0 by definition.
    // vector product of pep X pem
    float vpx = 0, vpy = 0, vpz = 0;
    if (t1.sign() * t2.sign() > 0) { // Like Sign
      if (bz * t1.sign() < 0) {
        vpx = v1.Py() * v2.Pz() - v1.Pz() * v2.Py();
        vpy = v1.Pz() * v2.Px() - v1.Px() * v2.Pz();
        vpz = v1.Px() * v2.Py() - v1.Py() * v2.Px();
      } else {
        vpx = v2.Py() * v1.Pz() - v2.Pz() * v1.Py();
        vpy = v2.Pz() * v1.Px() - v2.Px() * v1.Pz();
        vpz = v2.Px() * v1.Py() - v2.Py() * v1.Px();
      }
    } else { // Unlike Sign
      if (bz * t1.sign() > 0) {
        vpx = v1.Py() * v2.Pz() - v1.Pz() * v2.Py();
        vpy = v1.Pz() * v2.Px() - v1.Px() * v2.Pz();
        vpz = v1.Px() * v2.Py() - v1.Py() * v2.Px();
      } else {
        vpx = v2.Py() * v1.Pz() - v2.Pz() * v1.Py();
        vpy = v2.Pz() * v1.Px() - v2.Px() * v1.Pz();
        vpz = v2.Px() * v1.Py() - v2.Py() * v1.Px();
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
    values[kPairPhiv] = TMath::ACos(wx * ax + wy * ay); // phiv in [0,pi] //cosPhiV = wx * ax + wy * ay;
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

  float m1 = fgkElectronMass;
  float m2 = fgkElectronMass;
  if constexpr (pairType == kJpsiToMuMu) {
    m1 = fgkMuonMass;
    m2 = fgkMuonMass;
  }

  if constexpr (pairType == kElectronMuon) {
    m2 = fgkMuonMass;
  }

  ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), m1);
  ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), m2);
  ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
  values[kMass] = v12.M();
  values[kPt] = v12.Pt();
  values[kEta] = v12.Eta();
  values[kPhi] = v12.Phi();
  values[kRap] = -v12.Rapidity();
}

template <typename T1, typename T2>
void VarManager::FillPairMC(T1 const& t1, T2 const& t2, float* values, PairCandidateType pairType)
{
  if (!values) {
    values = fgValues;
  }

  float m1 = fgkElectronMass;
  float m2 = fgkElectronMass;
  if (pairType == kJpsiToMuMu) {
    m1 = fgkMuonMass;
    m2 = fgkMuonMass;
  }

  if (pairType == kElectronMuon) {
    m2 = fgkMuonMass;
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

template <int pairType, uint32_t collFillMap, uint32_t fillMap, typename C, typename T>
void VarManager::FillPairVertexing(C const& collision, T const& t1, T const& t2, float* values)
{
  // check at compile time that the event and cov matrix have the cov matrix
  constexpr bool eventHasVtxCov = ((collFillMap & Collision) > 0 || (collFillMap & ReducedEventVtxCov) > 0);
  constexpr bool trackHasCov = ((fillMap & TrackCov) > 0 || (fillMap & ReducedTrackBarrelCov) > 0);
  constexpr bool muonHasCov = ((fillMap & MuonCov) > 0 || (fillMap & ReducedMuonCov) > 0);

  if (!values) {
    values = fgValues;
  }

  int procCode = 0;

  // TODO: use trackUtilities functions to initialize the various matrices to avoid code duplication
  // auto pars1 = getTrackParCov(t1);
  // auto pars2 = getTrackParCov(t2);
  // We need to hide the cov data members from the cases when no cov table is provided
  if constexpr ((pairType == kJpsiToEE) && trackHasCov) {
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
  } else if constexpr ((pairType == kJpsiToMuMu) && muonHasCov) {
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
    return;
  }

  float m1 = fgkElectronMass;
  float m2 = fgkElectronMass;
  Vec3D secondaryVertex;
  float bz = 0;
  std::array<float, 3> pvec0;
  std::array<float, 3> pvec1;

  if constexpr (eventHasVtxCov) {
    std::array<float, 6> covMatrixPCA;
    o2::dataformats::DCA impactParameter0;
    o2::dataformats::DCA impactParameter1;
    // get track impact parameters
    // This modifies track momenta!
    o2::math_utils::Point3D<float> vtxXYZ(collision.posX(), collision.posY(), collision.posZ());
    std::array<float, 6> vtxCov{collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ()};
    o2::dataformats::VertexBase primaryVertex = {std::move(vtxXYZ), std::move(vtxCov)};
    // auto primaryVertex = getPrimaryVertex(collision);
    auto covMatrixPV = primaryVertex.getCov();

    if constexpr (pairType == kJpsiToEE && trackHasCov) {
      secondaryVertex = fgFitterTwoProngBarrel.getPCACandidate();
      bz = fgFitterTwoProngBarrel.getBz();
      covMatrixPCA = fgFitterTwoProngBarrel.calcPCACovMatrixFlat();
      auto chi2PCA = fgFitterTwoProngBarrel.getChi2AtPCACandidate();
      auto trackParVar0 = fgFitterTwoProngBarrel.getTrack(0);
      auto trackParVar1 = fgFitterTwoProngBarrel.getTrack(1);
      values[kVertexingChi2PCA] = chi2PCA;
      trackParVar0.getPxPyPzGlo(pvec0);
      trackParVar1.getPxPyPzGlo(pvec1);
      trackParVar0.propagateToDCA(primaryVertex, bz, &impactParameter0);
      trackParVar1.propagateToDCA(primaryVertex, bz, &impactParameter1);
    } else if constexpr (pairType == kJpsiToMuMu && muonHasCov) {
      // Get pca candidate from forward DCA fitter
      m1 = fgkMuonMass;
      m2 = fgkMuonMass;

      secondaryVertex = fgFitterTwoProngFwd.getPCACandidate();
      bz = fgFitterTwoProngFwd.getBz();
      covMatrixPCA = fgFitterTwoProngFwd.calcPCACovMatrixFlat();
      auto chi2PCA = fgFitterTwoProngFwd.getChi2AtPCACandidate();
      auto trackParVar0 = fgFitterTwoProngFwd.getTrack(0);
      auto trackParVar1 = fgFitterTwoProngFwd.getTrack(1);
      values[kVertexingChi2PCA] = chi2PCA;
      pvec0[0] = trackParVar0.getPx();
      pvec0[1] = trackParVar0.getPy();
      pvec0[2] = trackParVar0.getPz();
      pvec1[0] = trackParVar1.getPx();
      pvec1[1] = trackParVar1.getPy();
      pvec1[2] = trackParVar1.getPz();
      trackParVar0.propagateToZlinear(primaryVertex.getZ());
      trackParVar1.propagateToZlinear(primaryVertex.getZ());
    }
    // get uncertainty of the decay length
    // double phi, theta;
    // getPointDirection(array{collision.posX(), collision.posY(), collision.posZ()}, secondaryVertex, phi, theta);
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

    ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), m1);
    ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), m2);
    ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;

    values[kVertexingTauz] = (collision.posZ() - secondaryVertex[2]) * v12.M() / (TMath::Abs(v12.Pz()) * o2::constants::physics::LightSpeedCm2NS);
    values[kVertexingTauxy] = values[kVertexingLxy] * v12.M() / (v12.P() * o2::constants::physics::LightSpeedCm2NS);

    values[kVertexingTauzErr] = values[kVertexingLzErr] * v12.M() / (TMath::Abs(v12.Pz()) * o2::constants::physics::LightSpeedCm2NS);
    values[kVertexingTauxyErr] = values[kVertexingLxyErr] * v12.M() / (v12.P() * o2::constants::physics::LightSpeedCm2NS);
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
  float mlepton;

  int procCode = 0;
  int procCodeJpsi = 0;

  if constexpr ((candidateType == kBcToThreeMuons) && muonHasCov) {
    mlepton = fgkMuonMass;
    mtrack = fgkMuonMass;

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
  } else if constexpr ((candidateType == kBtoJpsiEEK) && trackHasCov) {
    mlepton = fgkElectronMass;
    mtrack = fgkKaonMass;
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

  ROOT::Math::PtEtaPhiMVector v1(lepton1.pt(), lepton1.eta(), lepton1.phi(), mlepton);
  ROOT::Math::PtEtaPhiMVector v2(lepton2.pt(), lepton2.eta(), lepton2.phi(), mlepton);
  ROOT::Math::PtEtaPhiMVector v3(track.pt(), track.eta(), track.phi(), mtrack);
  ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
  ROOT::Math::PtEtaPhiMVector vdilepton(v12.pt(), v12.eta(), v12.phi(), v12.M());
  ROOT::Math::PtEtaPhiMVector v123 = vdilepton + v3;
  values[VarManager::kPairMass] = v123.M();
  values[VarManager::kPairPt] = v123.Pt();
  values[VarManager::kPairEta] = v123.Eta();

  values[VarManager::kPairMassDau] = v12.M();
  values[VarManager::kPairPtDau] = v12.Pt();
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

    if constexpr (candidateType == kBtoJpsiEEK && trackHasCov) {
      secondaryVertex = fgFitterThreeProngBarrel.getPCACandidate();
      covMatrixPCA = fgFitterThreeProngBarrel.calcPCACovMatrixFlat();
    } else if constexpr (candidateType == kBcToThreeMuons && muonHasCov) {
      secondaryVertex = fgFitterThreeProngFwd.getPCACandidate();
      covMatrixPCA = fgFitterThreeProngFwd.calcPCACovMatrixFlat();
    }

    double phi = std::atan2(secondaryVertex[1] - collision.posY(), secondaryVertex[0] - collision.posX());
    double theta = std::atan2(secondaryVertex[2] - collision.posZ(),
                              std::sqrt((secondaryVertex[0] - collision.posX()) * (secondaryVertex[0] - collision.posX()) +
                                        (secondaryVertex[1] - collision.posY()) * (secondaryVertex[1] - collision.posY())));

    values[kVertexingLxyzErr] = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, theta) + getRotatedCovMatrixXX(covMatrixPCA, phi, theta));
    values[kVertexingLxyErr] = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, 0.) + getRotatedCovMatrixXX(covMatrixPCA, phi, 0.));
    values[kVertexingLzErr] = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, 0, theta) + getRotatedCovMatrixXX(covMatrixPCA, 0, theta));

    values[VarManager::kVertexingLxy] = (collision.posX() - secondaryVertex[0]) * (collision.posX() - secondaryVertex[0]) +
                                        (collision.posY() - secondaryVertex[1]) * (collision.posY() - secondaryVertex[1]);
    values[VarManager::kVertexingLz] = (collision.posZ() - secondaryVertex[2]) * (collision.posZ() - secondaryVertex[2]);
    values[VarManager::kVertexingLxyz] = values[VarManager::kVertexingLxy] + values[VarManager::kVertexingLz];
    values[VarManager::kVertexingLxy] = std::sqrt(values[VarManager::kVertexingLxy]);
    values[VarManager::kVertexingLz] = std::sqrt(values[VarManager::kVertexingLz]);
    values[VarManager::kVertexingLxyz] = std::sqrt(values[VarManager::kVertexingLxyz]);

    values[kVertexingTauz] = (collision.posZ() - secondaryVertex[2]) * v123.M() / (TMath::Abs(v123.Pz()) * o2::constants::physics::LightSpeedCm2NS);
    values[kVertexingTauxy] = values[kVertexingLxy] * v123.M() / (v123.P() * o2::constants::physics::LightSpeedCm2NS);

    values[kVertexingTauzErr] = values[kVertexingLzErr] * v123.M() / (TMath::Abs(v123.Pz()) * o2::constants::physics::LightSpeedCm2NS);
    values[kVertexingTauxyErr] = values[kVertexingLxyErr] * v123.M() / (v123.P() * o2::constants::physics::LightSpeedCm2NS);

    values[VarManager::kCosPointingAngle] = ((collision.posX() - secondaryVertex[0]) * v123.Px() +
                                             (collision.posY() - secondaryVertex[1]) * v123.Py() +
                                             (collision.posZ() - secondaryVertex[2]) * v123.Pz()) /
                                            (v123.P() * values[VarManager::kVertexingLxyz]);
  }
}

template <typename C, typename A>
void VarManager::FillQVectorFromGFW(C const& collision, A const& compA2, A const& compB2, A const& compC2, A const& compA3, A const& compB3, A const& compC3, float normA, float normB, float normC, float* values)
{
  if (!values) {
    values = fgValues;
  }

  // Fill Qn vectors from generic flow framework for different eta gap A, B, C (n=2,3)
  values[kQ2X0A] = compA2.Re() / normA;
  values[kQ2Y0A] = compA2.Im() / normA;
  values[kQ2X0B] = compB2.Re() / normB;
  values[kQ2Y0B] = compB2.Im() / normB;
  values[kQ2X0C] = compC2.Re() / normC;
  values[kQ2Y0C] = compC2.Im() / normC;
  values[kQ3X0A] = compA3.Re() / normA;
  values[kQ3Y0A] = compA3.Im() / normA;
  values[kQ3X0B] = compB3.Re() / normB;
  values[kQ3Y0B] = compB3.Im() / normB;
  values[kQ3X0C] = compC3.Re() / normC;
  values[kQ3Y0C] = compC3.Im() / normC;
  values[kMultA] = normA;
  values[kMultB] = normB;
  values[kMultC] = normC;

  // TODO: provide different computations for R
  // Compute the R factor using the 2 sub-events technique for second and third harmonic
  auto Psi2B = getEventPlane(2, values[kQ2X0B], values[kQ2Y0B]);
  auto Psi3B = getEventPlane(3, values[kQ3X0B], values[kQ3Y0B]);
  auto Psi2C = getEventPlane(2, values[kQ2X0C], values[kQ2Y0C]);
  auto Psi3C = getEventPlane(3, values[kQ3X0C], values[kQ3Y0C]);
  values[kR2SP] = (values[kQ2X0B] * values[kQ2X0C] + values[kQ2Y0B] * values[kQ2Y0C]);
  values[kR3SP] = (values[kQ3X0B] * values[kQ3X0C] + values[kQ3Y0B] * values[kQ3Y0C]);
  if (values[kQ2Y0B] * values[kQ2Y0C] != 0.0) {
    values[kR2EP] = TMath::Cos(2 * (Psi2B - Psi2C));
  }
  if (values[kQ3Y0B] * values[kQ3Y0C] != 0.0) {
    values[kR3EP] = TMath::Cos(3 * (Psi3B - Psi3C));
  }
}

template <int pairType, typename T1, typename T2>
void VarManager::FillPairVn(T1 const& t1, T2 const& t2, float* values)
{

  if (!values) {
    values = fgValues;
  }

  float m1 = fgkElectronMass;
  float m2 = fgkElectronMass;
  if constexpr (pairType == kJpsiToMuMu) {
    m1 = fgkMuonMass;
    m2 = fgkMuonMass;
  }

  if constexpr (pairType == kElectronMuon) {
    m2 = fgkMuonMass;
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
  values[kCos2DeltaPhi] = std::cos(2 * (v12.Phi() - getEventPlane(2, values[kQ2X0A], values[kQ2Y0A])));
  values[kCos3DeltaPhi] = std::cos(3 * (v12.Phi() - getEventPlane(3, values[kQ3X0A], values[kQ3Y0A])));
}

template <typename T1, typename T2>
void VarManager::FillDileptonHadron(T1 const& dilepton, T2 const& hadron, float* values, float hadronMass)
{
  if (!values) {
    values = fgValues;
  }

  if (fgUsedVars[kPairMass] || fgUsedVars[kPairPt] || fgUsedVars[kPairEta] || fgUsedVars[kPairPhi]) {
    ROOT::Math::PtEtaPhiMVector v1(dilepton.pt(), dilepton.eta(), dilepton.phi(), dilepton.mass());
    ROOT::Math::PtEtaPhiMVector v2(hadron.pt(), hadron.eta(), hadron.phi(), hadronMass);
    ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
    values[kPairMass] = v12.M();
    values[kPairPt] = v12.Pt();
    values[kPairEta] = v12.Eta();
    values[kPairPhi] = v12.Phi();
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
#endif // PWGDQ_CORE_VARMANAGER_H_
