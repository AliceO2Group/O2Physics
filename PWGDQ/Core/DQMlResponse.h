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

/// \file DQMlResponse.h
/// \brief Class to compute the ML response
/// \author Daniel Samitz <daniel.samitz@cern.ch>, SMI Vienna
///         Elisa Meninno, <elisa.meninno@cern.ch>, SMI Vienna

#ifndef PWGDQ_CORE_DQMLRESPONSE_H_
#define PWGDQ_CORE_DQMLRESPONSE_H_

#include <map>
#include <string>
#include <vector>

#include "Tools/ML/MlResponse.h"
#include "PWGDQ/Core/VarManager.h"

namespace dq_ml_cuts
{
// direction of the cut
enum CutDirection {
  CutGreater = 0, // require score < cut value
  CutSmaller,     // require score > cut value
  CutNot          // do not cut on score
};

static constexpr int nBins = 1;
static constexpr int nCutScores = 2;
// default values for the pT bin edges), offset by 1 from the bin numbers in cuts array
constexpr double bins[nBins + 1] = {
  0.,
  9999.};
auto vecBins = std::vector<double>{bins, bins + nBins + 1};

// default values for the cut directions
constexpr int cutDir[nCutScores] = {CutSmaller, CutGreater};
auto vecCutDir = std::vector<int>{cutDir, cutDir + nCutScores};

// default values for the cuts
constexpr double cuts[nBins][nCutScores] = {
  {0.5, 0.5}};

// row labels
static const std::vector<std::string> labels = {
  "bin 0"};

// column labels
static const std::vector<std::string> labelsCutScore = {"Signal", "Background"};
} // namespace dq_ml_cuts

#define FILL_MAP_VARMANAGER(FEATURE)                                 \
  {                                                                  \
#FEATURE, static_cast < uint8_t>(VarManager::Variables::FEATURE) \
  }

namespace o2::analysis
{

template <typename TypeOutputScore = float>
class DQMlResponse : public MlResponse<TypeOutputScore>
{
 public:
  /// Default constructor
  DQMlResponse() = default;
  /// Default destructor
  virtual ~DQMlResponse() = default;

  std::vector<float> getInputFeatures()
  {
    std::vector<float> inputFeatures;
    for (const auto& idx : MlResponse<TypeOutputScore>::mCachedIndices) {
      inputFeatures.emplace_back(VarManager::fgValues[idx]);
    }
    return inputFeatures;
  }

  void cacheBinVariableIndex(std::string const& cfgBinVariable)
  {
    mCachedBinVariableIndex = MlResponse<TypeOutputScore>::mAvailableInputFeatures[cfgBinVariable];
  }

  float getBinVariable()
  {
    return VarManager::fgValues[mCachedBinVariableIndex];
  }

  uint32_t getCachedBinVariableIndex()
  {
    return mCachedBinVariableIndex;
  }

  std::vector<int> getCachedIndices()
  {
    std::vector<int> convertedVec = {};
    for (int v : MlResponse<TypeOutputScore>::mCachedIndices) {
      convertedVec.push_back(v);
    }
    return convertedVec;
  }

 protected:
  /// Method to fill the map of available input features
  void setAvailableInputFeatures()
  {
    // TODO: give an error for variables that are not set
    MlResponse<TypeOutputScore>::mAvailableInputFeatures = {
      // Event wise variables
      FILL_MAP_VARMANAGER(kTimestamp),
      FILL_MAP_VARMANAGER(kBC),
      FILL_MAP_VARMANAGER(kIsPhysicsSelection),
      FILL_MAP_VARMANAGER(kIsSel8), // TVX in Run3
      FILL_MAP_VARMANAGER(kIsINT7),
      FILL_MAP_VARMANAGER(kIsEMC7),
      FILL_MAP_VARMANAGER(kIsINT7inMUON),
      FILL_MAP_VARMANAGER(kIsMuonSingleLowPt7),
      FILL_MAP_VARMANAGER(kIsMuonSingleHighPt7),
      FILL_MAP_VARMANAGER(kIsMuonUnlikeLowPt7),
      FILL_MAP_VARMANAGER(kIsMuonLikeLowPt7),
      FILL_MAP_VARMANAGER(kIsCUP8),
      FILL_MAP_VARMANAGER(kIsCUP9),
      FILL_MAP_VARMANAGER(kIsMUP10),
      FILL_MAP_VARMANAGER(kIsMUP11),
      FILL_MAP_VARMANAGER(kVtxX),
      FILL_MAP_VARMANAGER(kVtxY),
      FILL_MAP_VARMANAGER(kVtxZ),
      FILL_MAP_VARMANAGER(kVtxNcontrib),
      FILL_MAP_VARMANAGER(kVtxNcontribReal),
      FILL_MAP_VARMANAGER(kVtxCovXX),
      FILL_MAP_VARMANAGER(kVtxCovXY),
      FILL_MAP_VARMANAGER(kVtxCovXZ),
      FILL_MAP_VARMANAGER(kVtxCovYY),
      FILL_MAP_VARMANAGER(kVtxCovYZ),
      FILL_MAP_VARMANAGER(kVtxCovZZ),
      FILL_MAP_VARMANAGER(kVtxChi2),
      FILL_MAP_VARMANAGER(kCentVZERO),
      FILL_MAP_VARMANAGER(kCentFT0C),
      FILL_MAP_VARMANAGER(kMultTPC),
      FILL_MAP_VARMANAGER(kMultFV0A),
      FILL_MAP_VARMANAGER(kMultFV0C),
      FILL_MAP_VARMANAGER(kMultFT0A),
      FILL_MAP_VARMANAGER(kMultFT0C),
      FILL_MAP_VARMANAGER(kMultFDDA),
      FILL_MAP_VARMANAGER(kMultFDDC),
      FILL_MAP_VARMANAGER(kMultZNA),
      FILL_MAP_VARMANAGER(kMultZNC),
      FILL_MAP_VARMANAGER(kMultTracklets),
      FILL_MAP_VARMANAGER(kMultDimuons),
      FILL_MAP_VARMANAGER(kMCEventGeneratorId),
      FILL_MAP_VARMANAGER(kMCVtxX),
      FILL_MAP_VARMANAGER(kMCVtxY),
      FILL_MAP_VARMANAGER(kMCVtxZ),
      FILL_MAP_VARMANAGER(kMCEventTime),
      FILL_MAP_VARMANAGER(kMCEventWeight),
      FILL_MAP_VARMANAGER(kMCEventImpParam),
      FILL_MAP_VARMANAGER(kQ1X0A), // q-vector (e.g. from TPC) with x component (harmonic 1 and power 0)), sub-event A
      FILL_MAP_VARMANAGER(kQ1Y0A), // q-vector (e.g. from TPC) with y component (harmonic 1 and power 0)), sub-event A
      FILL_MAP_VARMANAGER(kQ1X0B),
      FILL_MAP_VARMANAGER(kQ1Y0B),
      FILL_MAP_VARMANAGER(kQ1X0C),
      FILL_MAP_VARMANAGER(kQ1Y0C),
      FILL_MAP_VARMANAGER(kQ2X0A), // q-vector (e.g. from TPC) with x component (harmonic 2 and power 0)), sub-event A
      FILL_MAP_VARMANAGER(kQ2Y0A), // q-vector (e.g. from TPC) with y component (harmonic 2 and power 0)), sub-event A
      FILL_MAP_VARMANAGER(kQ2X0B),
      FILL_MAP_VARMANAGER(kQ2Y0B),
      FILL_MAP_VARMANAGER(kQ2X0C),
      FILL_MAP_VARMANAGER(kQ2Y0C),
      FILL_MAP_VARMANAGER(kMultA), // Multiplicity of the sub-event A
      FILL_MAP_VARMANAGER(kMultB),
      FILL_MAP_VARMANAGER(kMultC),
      FILL_MAP_VARMANAGER(kQ3X0A), // q-vector (e.g. from TPC) with x component (harmonic 3 and power 0)), sub-event A
      FILL_MAP_VARMANAGER(kQ3Y0A), // q-vector (e.g. from TPC) with y component (harmonic 3 and power 0)), sub-event A
      FILL_MAP_VARMANAGER(kQ3X0B),
      FILL_MAP_VARMANAGER(kQ3Y0B),
      FILL_MAP_VARMANAGER(kQ3X0C),
      FILL_MAP_VARMANAGER(kQ3Y0C),
      FILL_MAP_VARMANAGER(kQ4X0A), // q-vector (e.g. from TPC) with x component (harmonic 4 and power 0)), sub-event A
      FILL_MAP_VARMANAGER(kQ4Y0A), // q-vector (e.g. from TPC) with y component (harmonic 4 and power 0)), sub-event A
      FILL_MAP_VARMANAGER(kQ4X0B),
      FILL_MAP_VARMANAGER(kQ4Y0B),
      FILL_MAP_VARMANAGER(kQ4X0C),
      FILL_MAP_VARMANAGER(kQ4Y0C),
      FILL_MAP_VARMANAGER(kR2SP),
      FILL_MAP_VARMANAGER(kR3SP),
      FILL_MAP_VARMANAGER(kR2EP),
      FILL_MAP_VARMANAGER(kR3EP),
      FILL_MAP_VARMANAGER(kIsDoubleGap),  // Double rapidity gap
      FILL_MAP_VARMANAGER(kIsSingleGapA), // Rapidity gap on side A
      FILL_MAP_VARMANAGER(kIsSingleGapC), // Rapidity gap on side C
      FILL_MAP_VARMANAGER(kIsSingleGap),  // Rapidity gap on either side
      FILL_MAP_VARMANAGER(kTwoEvPosZ1),   // vtx-z for collision 1 in two events correlations
      FILL_MAP_VARMANAGER(kTwoEvPosZ2),   // vtx-z for collision 2 in two events correlations
      FILL_MAP_VARMANAGER(kTwoEvPosR1),   // vtx-R for collision 1 in two events correlations
      FILL_MAP_VARMANAGER(kTwoEvPosR2),
      FILL_MAP_VARMANAGER(kTwoEvPVcontrib1), // n-contributors for collision 1 in two events correlations
      FILL_MAP_VARMANAGER(kTwoEvPVcontrib2),
      FILL_MAP_VARMANAGER(kTwoEvDeltaZ), // distance in z between collisions
      FILL_MAP_VARMANAGER(kTwoEvDeltaX), // distance in x between collisions
      FILL_MAP_VARMANAGER(kTwoEvDeltaY), // distance in y between collisions
      FILL_MAP_VARMANAGER(kTwoEvDeltaR), // distance in (x),y) plane between collisions
      // Basic track/muon/pair wise variables
      FILL_MAP_VARMANAGER(kX),
      FILL_MAP_VARMANAGER(kY),
      FILL_MAP_VARMANAGER(kZ),
      FILL_MAP_VARMANAGER(kPt),
      FILL_MAP_VARMANAGER(kSignedPt),
      FILL_MAP_VARMANAGER(kInvPt),
      FILL_MAP_VARMANAGER(kEta),
      FILL_MAP_VARMANAGER(kTgl),
      FILL_MAP_VARMANAGER(kPhi),
      FILL_MAP_VARMANAGER(kP),
      FILL_MAP_VARMANAGER(kPx),
      FILL_MAP_VARMANAGER(kPy),
      FILL_MAP_VARMANAGER(kPz),
      FILL_MAP_VARMANAGER(kRap),
      FILL_MAP_VARMANAGER(kMass),
      FILL_MAP_VARMANAGER(kCharge),
      FILL_MAP_VARMANAGER(kNBasicTrackVariables),
      FILL_MAP_VARMANAGER(kUsedKF),
      FILL_MAP_VARMANAGER(kKFMass),
      FILL_MAP_VARMANAGER(kPt1),
      FILL_MAP_VARMANAGER(kEta1),
      FILL_MAP_VARMANAGER(kPhi1),
      FILL_MAP_VARMANAGER(kCharge1),
      FILL_MAP_VARMANAGER(kPt2),
      FILL_MAP_VARMANAGER(kEta2),
      FILL_MAP_VARMANAGER(kPhi2),
      FILL_MAP_VARMANAGER(kCharge2),
      // Barrel track variables
      FILL_MAP_VARMANAGER(kPin),
      FILL_MAP_VARMANAGER(kSignedPin),
      FILL_MAP_VARMANAGER(kTOFExpMom),
      FILL_MAP_VARMANAGER(kTrackTime),
      FILL_MAP_VARMANAGER(kTrackTimeRes),
      FILL_MAP_VARMANAGER(kTrackTimeResRelative),
      FILL_MAP_VARMANAGER(kDetectorMap),
      FILL_MAP_VARMANAGER(kHasITS),
      FILL_MAP_VARMANAGER(kHasTRD),
      FILL_MAP_VARMANAGER(kHasTOF),
      FILL_MAP_VARMANAGER(kHasTPC),
      FILL_MAP_VARMANAGER(kIsGlobalTrack),
      FILL_MAP_VARMANAGER(kIsGlobalTrackSDD),
      FILL_MAP_VARMANAGER(kIsITSrefit),
      FILL_MAP_VARMANAGER(kIsSPDany),
      FILL_MAP_VARMANAGER(kIsSPDfirst),
      FILL_MAP_VARMANAGER(kIsSPDboth),
      FILL_MAP_VARMANAGER(kIsITSibAny),
      FILL_MAP_VARMANAGER(kIsITSibFirst),
      FILL_MAP_VARMANAGER(kIsITSibAll),
      FILL_MAP_VARMANAGER(kITSncls),
      FILL_MAP_VARMANAGER(kITSchi2),
      FILL_MAP_VARMANAGER(kITSlayerHit),
      FILL_MAP_VARMANAGER(kITSmeanClsSize),
      FILL_MAP_VARMANAGER(kIsTPCrefit),
      FILL_MAP_VARMANAGER(kTPCncls),
      FILL_MAP_VARMANAGER(kITSClusterMap),
      FILL_MAP_VARMANAGER(kTPCnclsCR),
      FILL_MAP_VARMANAGER(kTPCchi2),
      FILL_MAP_VARMANAGER(kTPCsignal),
      FILL_MAP_VARMANAGER(kTPCsignalRandomized),
      FILL_MAP_VARMANAGER(kTPCsignalRandomizedDelta),
      FILL_MAP_VARMANAGER(kTRDsignal),
      FILL_MAP_VARMANAGER(kTRDPattern),
      FILL_MAP_VARMANAGER(kTOFbeta),
      FILL_MAP_VARMANAGER(kTrackLength),
      FILL_MAP_VARMANAGER(kTrackDCAxy),
      FILL_MAP_VARMANAGER(kTrackDCAz),
      FILL_MAP_VARMANAGER(kTrackDCAsigXY),
      FILL_MAP_VARMANAGER(kTrackDCAsigZ),
      FILL_MAP_VARMANAGER(kTrackDCAresXY),
      FILL_MAP_VARMANAGER(kTrackDCAresZ),
      FILL_MAP_VARMANAGER(kIsGoldenChi2),
      FILL_MAP_VARMANAGER(kTrackCYY),
      FILL_MAP_VARMANAGER(kTrackCZZ),
      FILL_MAP_VARMANAGER(kTrackCSnpSnp),
      FILL_MAP_VARMANAGER(kTrackCTglTgl),
      FILL_MAP_VARMANAGER(kTrackC1Pt21Pt2),
      FILL_MAP_VARMANAGER(kTPCnSigmaEl),
      FILL_MAP_VARMANAGER(kTPCnSigmaElRandomized),
      FILL_MAP_VARMANAGER(kTPCnSigmaElRandomizedDelta),
      FILL_MAP_VARMANAGER(kTPCnSigmaMu),
      FILL_MAP_VARMANAGER(kTPCnSigmaPi),
      FILL_MAP_VARMANAGER(kTPCnSigmaPiRandomized),
      FILL_MAP_VARMANAGER(kTPCnSigmaPiRandomizedDelta),
      FILL_MAP_VARMANAGER(kTPCnSigmaKa),
      FILL_MAP_VARMANAGER(kTPCnSigmaPr),
      FILL_MAP_VARMANAGER(kTPCnSigmaEl_Corr),
      FILL_MAP_VARMANAGER(kTPCnSigmaPi_Corr),
      FILL_MAP_VARMANAGER(kTPCnSigmaKa_Corr),
      FILL_MAP_VARMANAGER(kTPCnSigmaPr_Corr),
      FILL_MAP_VARMANAGER(kTPCnSigmaPrRandomized),
      FILL_MAP_VARMANAGER(kTPCnSigmaPrRandomizedDelta),
      FILL_MAP_VARMANAGER(kTOFnSigmaEl),
      FILL_MAP_VARMANAGER(kTOFnSigmaMu),
      FILL_MAP_VARMANAGER(kTOFnSigmaPi),
      FILL_MAP_VARMANAGER(kTOFnSigmaKa),
      FILL_MAP_VARMANAGER(kTOFnSigmaPr),
      FILL_MAP_VARMANAGER(kTrackTimeResIsRange), // Gaussian or range (see Framework/DataTypes)
      FILL_MAP_VARMANAGER(kPVContributor),       // This track has contributed to the collision vertex fit (see Framework/DataTypes)
      FILL_MAP_VARMANAGER(kOrphanTrack),         // Track has no association with any collision vertex (see Framework/DataTypes)
      FILL_MAP_VARMANAGER(kIsAmbiguous),
      FILL_MAP_VARMANAGER(kIsLegFromGamma),
      FILL_MAP_VARMANAGER(kIsLegFromK0S),
      FILL_MAP_VARMANAGER(kIsLegFromLambda),
      FILL_MAP_VARMANAGER(kIsLegFromAntiLambda),
      FILL_MAP_VARMANAGER(kIsLegFromOmega),
      FILL_MAP_VARMANAGER(kIsProtonFromLambdaAndAntiLambda),
      FILL_MAP_VARMANAGER(kIsDalitzLeg),             // Up to 8 dalitz selections
      FILL_MAP_VARMANAGER(kBarrelNAssocsInBunch),    // number of in bunch collision associations
      FILL_MAP_VARMANAGER(kBarrelNAssocsOutOfBunch), // number of out of bunch collision associations
      // Muon track variables
      FILL_MAP_VARMANAGER(kMuonNClusters),
      FILL_MAP_VARMANAGER(kMuonPDca),
      FILL_MAP_VARMANAGER(kMuonRAtAbsorberEnd),
      FILL_MAP_VARMANAGER(kMCHBitMap),
      FILL_MAP_VARMANAGER(kMuonChi2),
      FILL_MAP_VARMANAGER(kMuonChi2MatchMCHMID),
      FILL_MAP_VARMANAGER(kMuonChi2MatchMCHMFT),
      FILL_MAP_VARMANAGER(kMuonMatchScoreMCHMFT),
      FILL_MAP_VARMANAGER(kMuonCXX),
      FILL_MAP_VARMANAGER(kMuonCXY),
      FILL_MAP_VARMANAGER(kMuonCYY),
      FILL_MAP_VARMANAGER(kMuonCPhiX),
      FILL_MAP_VARMANAGER(kMuonCPhiY),
      FILL_MAP_VARMANAGER(kMuonCPhiPhi),
      FILL_MAP_VARMANAGER(kMuonCTglX),
      FILL_MAP_VARMANAGER(kMuonCTglY),
      FILL_MAP_VARMANAGER(kMuonCTglPhi),
      FILL_MAP_VARMANAGER(kMuonCTglTgl),
      FILL_MAP_VARMANAGER(kMuonC1Pt2X),
      FILL_MAP_VARMANAGER(kMuonC1Pt2Y),
      FILL_MAP_VARMANAGER(kMuonC1Pt2Phi),
      FILL_MAP_VARMANAGER(kMuonC1Pt2Tgl),
      FILL_MAP_VARMANAGER(kMuonC1Pt21Pt2),
      FILL_MAP_VARMANAGER(kMuonTrackType),
      FILL_MAP_VARMANAGER(kMuonDCAx),
      FILL_MAP_VARMANAGER(kMuonDCAy),
      FILL_MAP_VARMANAGER(kMuonTime),
      FILL_MAP_VARMANAGER(kMuonTimeRes),
      FILL_MAP_VARMANAGER(kMftNClusters),
      FILL_MAP_VARMANAGER(kMftClusterSize),
      FILL_MAP_VARMANAGER(kMftMeanClusterSize),
      // Pair variables
      FILL_MAP_VARMANAGER(kCandidateId),
      FILL_MAP_VARMANAGER(kPairType),
      FILL_MAP_VARMANAGER(kVertexingLxy),
      FILL_MAP_VARMANAGER(kVertexingLxyErr),
      FILL_MAP_VARMANAGER(kVertexingPseudoCTau),
      FILL_MAP_VARMANAGER(kVertexingLxyz),
      FILL_MAP_VARMANAGER(kVertexingLxyzErr),
      FILL_MAP_VARMANAGER(kVertexingLz),
      FILL_MAP_VARMANAGER(kVertexingLzErr),
      FILL_MAP_VARMANAGER(kVertexingTauxy),
      FILL_MAP_VARMANAGER(kVertexingTauxyErr),
      FILL_MAP_VARMANAGER(kVertexingLzProjected),
      FILL_MAP_VARMANAGER(kVertexingLxyProjected),
      FILL_MAP_VARMANAGER(kVertexingLxyzProjected),
      FILL_MAP_VARMANAGER(kVertexingTauzProjected),
      FILL_MAP_VARMANAGER(kVertexingTauxyProjected),
      FILL_MAP_VARMANAGER(kVertexingTauxyProjectedNs),
      FILL_MAP_VARMANAGER(kVertexingTauz),
      FILL_MAP_VARMANAGER(kVertexingTauzErr),
      FILL_MAP_VARMANAGER(kVertexingPz),
      FILL_MAP_VARMANAGER(kVertexingSV),
      FILL_MAP_VARMANAGER(kVertexingProcCode),
      FILL_MAP_VARMANAGER(kVertexingChi2PCA),
      FILL_MAP_VARMANAGER(kCosThetaHE),
      FILL_MAP_VARMANAGER(kCosThetaCS),
      FILL_MAP_VARMANAGER(kPhiHE),
      FILL_MAP_VARMANAGER(kPhiCS),
      FILL_MAP_VARMANAGER(kPsiPair),
      FILL_MAP_VARMANAGER(kDeltaPhiPair),
      FILL_MAP_VARMANAGER(kOpeningAngle),
      FILL_MAP_VARMANAGER(kQuadDCAabsXY),
      FILL_MAP_VARMANAGER(kQuadDCAsigXY),
      FILL_MAP_VARMANAGER(kQuadDCAabsZ),
      FILL_MAP_VARMANAGER(kQuadDCAsigZ),
      FILL_MAP_VARMANAGER(kQuadDCAsigXYZ),
      FILL_MAP_VARMANAGER(kCosPointingAngle),
      FILL_MAP_VARMANAGER(kImpParXYJpsi),
      FILL_MAP_VARMANAGER(kImpParXYK),
      FILL_MAP_VARMANAGER(kDCATrackProd),
      FILL_MAP_VARMANAGER(kDCATrackVtxProd),
      FILL_MAP_VARMANAGER(kU2Q2),
      FILL_MAP_VARMANAGER(kU3Q3),
      FILL_MAP_VARMANAGER(kCORR2REF),
      FILL_MAP_VARMANAGER(kCORR2POI),
      FILL_MAP_VARMANAGER(kCORR4REF),
      FILL_MAP_VARMANAGER(kCORR4POI),
      FILL_MAP_VARMANAGER(kC4REF),
      FILL_MAP_VARMANAGER(kC4POI),
      FILL_MAP_VARMANAGER(kV4),
      FILL_MAP_VARMANAGER(kPsi2A),
      FILL_MAP_VARMANAGER(kPsi2B),
      FILL_MAP_VARMANAGER(kPsi2C),
      FILL_MAP_VARMANAGER(kCos2DeltaPhi),
      FILL_MAP_VARMANAGER(kCos3DeltaPhi),
      FILL_MAP_VARMANAGER(kDeltaPtotTracks),
      FILL_MAP_VARMANAGER(kVertexingLxyOverErr),
      FILL_MAP_VARMANAGER(kVertexingLzOverErr),
      FILL_MAP_VARMANAGER(kVertexingLxyzOverErr),
      FILL_MAP_VARMANAGER(kKFTrack0DCAxyz),
      FILL_MAP_VARMANAGER(kKFTrack1DCAxyz),
      FILL_MAP_VARMANAGER(kKFTracksDCAxyzMax),
      FILL_MAP_VARMANAGER(kKFDCAxyzBetweenProngs),
      FILL_MAP_VARMANAGER(kKFTrack0DCAxy),
      FILL_MAP_VARMANAGER(kKFTrack1DCAxy),
      FILL_MAP_VARMANAGER(kKFTracksDCAxyMax),
      FILL_MAP_VARMANAGER(kKFDCAxyBetweenProngs),
      FILL_MAP_VARMANAGER(kKFChi2OverNDFGeo),
      FILL_MAP_VARMANAGER(kKFNContributorsPV),
      FILL_MAP_VARMANAGER(kKFCosPA),
      // Candidate-track correlation variables
      FILL_MAP_VARMANAGER(kPairMass),
      FILL_MAP_VARMANAGER(kPairMassDau),
      FILL_MAP_VARMANAGER(kMassDau),
      FILL_MAP_VARMANAGER(kPairPt),
      FILL_MAP_VARMANAGER(kPairPtDau),
      FILL_MAP_VARMANAGER(kPairEta),
      FILL_MAP_VARMANAGER(kPairPhi),
      FILL_MAP_VARMANAGER(kPairPhiv),
      FILL_MAP_VARMANAGER(kDeltaEta),
      FILL_MAP_VARMANAGER(kDeltaPhi),
      FILL_MAP_VARMANAGER(kDeltaPhiSym),
      // DQ-HF correlation variables
      FILL_MAP_VARMANAGER(kMassCharmHadron),
      FILL_MAP_VARMANAGER(kPtCharmHadron),
      FILL_MAP_VARMANAGER(kRapCharmHadron),
      FILL_MAP_VARMANAGER(kPhiCharmHadron),
      FILL_MAP_VARMANAGER(kBdtCharmHadron)};
  }

  uint8_t mCachedBinVariableIndex;
};

} // namespace o2::analysis

#undef FILL_MAP_VARMANAGER

#endif // PWGDQ_CORE_DQMLRESPONSE_H_
