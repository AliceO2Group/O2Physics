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

#include "PWGDQ/Core/VarManager.h"

#include <cmath>

ClassImp(VarManager);

TString VarManager::fgVariableNames[VarManager::kNVars] = {""};
TString VarManager::fgVariableUnits[VarManager::kNVars] = {""};
bool VarManager::fgUsedVars[VarManager::kNVars] = {kFALSE};
float VarManager::fgValues[VarManager::kNVars] = {0.0f};
std::map<int, int> VarManager::fgRunMap;
TString VarManager::fgRunStr = "";
o2::vertexing::DCAFitterN<2> VarManager::fgFitterTwoProngBarrel;
o2::vertexing::DCAFitterN<3> VarManager::fgFitterThreeProngBarrel;
o2::vertexing::FwdDCAFitterN<2> VarManager::fgFitterTwoProngFwd;
o2::vertexing::FwdDCAFitterN<3> VarManager::fgFitterThreeProngFwd;

//__________________________________________________________________
VarManager::VarManager() : TObject()
{
  //
  // constructor
  //
  SetDefaultVarNames();
}

//__________________________________________________________________
VarManager::~VarManager() = default;

//__________________________________________________________________
void VarManager::SetVariableDependencies()
{
  //
  // Set as used variables on which other variables calculation depends
  //
  if (fgUsedVars[kP]) {
    fgUsedVars[kPt] = kTRUE;
    fgUsedVars[kEta] = kTRUE;
  }
}

//__________________________________________________________________
void VarManager::ResetValues(int startValue, int endValue, float* values)
{
  //
  // reset all variables to an "innocent" value
  // NOTE: here we use -9999.0 as a neutral value, but depending on situation, this may not be the case
  if (!values) {
    values = fgValues;
  }
  for (Int_t i = startValue; i < endValue; ++i) {
    values[i] = -9999.;
  }
}

//__________________________________________________________________
void VarManager::SetRunNumbers(int n, int* runs)
{
  //
  // maps the list of runs such that one can plot the list of runs nicely in a histogram axis
  //
  for (int i = 0; i < n; ++i) {
    fgRunMap[runs[i]] = i + 1;
    fgRunStr += Form("%d;", runs[i]);
  }
}

//__________________________________________________________________
void VarManager::SetRunNumbers(std::vector<int> runs)
{
  //
  // maps the list of runs such that one can plot the list of runs nicely in a histogram axis
  //
  int i = 0;
  for (auto run = runs.begin(); run != runs.end(); run++, i++) {
    fgRunMap[*run] = i + 1;
    fgRunStr += Form("%d;", *run);
  }
}

//__________________________________________________________________
void VarManager::FillEventDerived(float* values)
{
  //
  // Fill event-wise derived quantities (these are all quantities which can be computed just based on the values already filled in the FillEvent() function)
  //
  if (fgUsedVars[kRunId]) {
    values[kRunId] = (fgRunMap.size() > 0 ? fgRunMap[static_cast<int>(values[kRunNo])] : 0);
  }
}

//__________________________________________________________________
void VarManager::FillTrackDerived(float* values)
{
  //
  // Fill track-wise derived quantities (these are all quantities which can be computed just based on the values already filled in the FillTrack() function)
  //
  if (fgUsedVars[kP]) {
    values[kP] = values[kPt] * std::cosh(values[kEta]);
  }
}
//_________________________________________________________________________________________________________________________________________________________________________________
float VarManager::GetTPCPostCalibMap(float pin, float eta, int particle_type, TString period)
{
  if (period.Contains("LHC22m_pass1_subset")) {
    float El_mean_curve_pin = (pin < 0.3) ? 1.74338 : ((pin < 3.5) ? 1 / (0.694318 - 5.66879 * pin) + 2.73696 + 0.000483342 * pin : 2.70321);
    float Pi_mean_curve_pin = (pin < 0.3) ? 1.95561 : ((pin < 3) ? -0.0606029 / pin + 2.18796 - 0.101135 * pin : 1.86435);
    float Pr_mean_curve_pin = (pin < 0.4) ? 1.94664 : ((pin < 5) ? 1 / (-0.495484 - 1.03622 * pin) + 3.0513 - 0.0143036 * pin : 2.80362);
    float El_mean_curve_eta = (std::abs(eta) < 0.9) ? -1.17942e-01 + -1.51932e-01 * std::cos(4.32572e+00 * eta) : 0.0;
    float Pi_mean_curve_eta = (std::abs(eta) < 0.9) ? -5.21188e-02 + -3.35413e-01 * std::cos(4.16710e+00 * eta) : 0.0;
    float Pr_mean_curve_eta = (std::abs(eta) < 0.9) ? -1.33413e-02 + -3.83269e-01 * std::cos(3.95128e+00 * eta) : 0.0;

    float pin_map = (particle_type == 0) ? El_mean_curve_pin : ((particle_type == 1) ? Pi_mean_curve_pin : Pr_mean_curve_pin);
    float eta_map = (particle_type == 0) ? El_mean_curve_eta : ((particle_type == 1) ? Pi_mean_curve_eta : Pr_mean_curve_eta);
    float map = pin_map + eta_map;
    return map;
  } else if (period.Contains("LHC22f_pass1")) {
    float El_mean_curve_pin = (pin < 0.3) ? 0.24335236 : ((pin < 3.5) ? 1 / (0.0113621 - 2.44516 * pin) + 1.63907 - 0.0367754 * pin : 1.3933518);
    float Pi_mean_curve_pin = (pin < 0.3) ? -0.30059726 : ((pin < 3.5) ? 1 / (-4.51007e+06 - 5.52635e+06 * pin) - 0.349193 + 0.171139 * pin - 0.0305089 * pin * pin : -0.12394057);
    float Pr_mean_curve_pin = (pin < 0.3) ? 0.1 : ((pin < 3.5) ? 1 / (0.482973 - 3.55557 * pin) + 1.38574 - 0.627066 * pin + 0.103612 * pin * pin : 0.37665460);
    float El_mean_curve_eta = (std::abs(eta) < 0.9) ? -0.25877 + -0.136459 * eta : 0.0;
    float Pi_mean_curve_eta = (std::abs(eta) < 0.9) ? -0.250769 + -0.321296 * eta + 0.509874 * eta * eta + 0.445708 * eta * eta * eta : 0.0;
    float Pr_mean_curve_eta = (std::abs(eta) < 0.9) ? -0.330935 + -0.395157 * eta + 0.582457 * eta * eta + 0.501215 * eta * eta * eta : 0.0;

    float pin_map = (particle_type == 0) ? El_mean_curve_pin : ((particle_type == 1) ? Pi_mean_curve_pin : Pr_mean_curve_pin);
    float eta_map = (particle_type == 0) ? El_mean_curve_eta : ((particle_type == 1) ? Pi_mean_curve_eta : Pr_mean_curve_eta);

    float map = pin_map + eta_map;
    return map;
  } else {
    float map = 0.0;
    return map;
  }
}
//__________________________________________________________________
TString VarManager::GetRunPeriod(float runNumber)
{
  int runlist_22f[2] = {520259, 520473};
  int runlist_22m[2] = {523393, 523397};

  if (runNumber >= runlist_22f[0] && runNumber <= runlist_22f[1]) {
    TString runperiod = "LHC22f_pass1";
    return runperiod;
  } else if (runNumber >= runlist_22m[0] && runNumber <= runlist_22m[1]) {
    TString runperiod = "LHC22m_pass1_subset";
    return runperiod;
  } else {
    TString runperiod = "none";
    // LOGF(info, "can't find run period for run %.0d", runNumber);
    return runperiod;
  }
};
//__________________________________________________________________
void VarManager::SetDefaultVarNames()
{
  //
  // Set default variable names
  //
  for (Int_t ivar = 0; ivar < kNVars; ++ivar) {
    fgVariableNames[ivar] = "DEFAULT NOT DEFINED";
    fgVariableUnits[ivar] = "n/a";
  }

  fgVariableNames[kRunNo] = "Run number";
  fgVariableUnits[kRunNo] = "";
  fgVariableNames[kRunId] = "Run number";
  fgVariableUnits[kRunId] = "";
  fgVariableNames[kBC] = "Bunch crossing";
  fgVariableUnits[kBC] = "";
  fgVariableNames[kIsPhysicsSelection] = "Physics selection";
  fgVariableUnits[kIsPhysicsSelection] = "";
  fgVariableNames[kVtxX] = "Vtx X ";
  fgVariableUnits[kVtxX] = "cm";
  fgVariableNames[kVtxY] = "Vtx Y ";
  fgVariableUnits[kVtxY] = "cm";
  fgVariableNames[kVtxZ] = "Vtx Z ";
  fgVariableUnits[kVtxZ] = "cm";
  fgVariableNames[kVtxNcontrib] = "Vtx contrib.";
  fgVariableUnits[kVtxNcontrib] = "";
  fgVariableNames[kVtxCovXX] = "Vtx covXX";
  fgVariableUnits[kVtxCovXX] = "cm";
  fgVariableNames[kVtxCovXY] = "Vtx covXY";
  fgVariableUnits[kVtxCovXY] = "cm";
  fgVariableNames[kVtxCovXZ] = "Vtx covXZ";
  fgVariableUnits[kVtxCovXZ] = "cm";
  fgVariableNames[kVtxCovYY] = "Vtx covYY";
  fgVariableUnits[kVtxCovYY] = "cm";
  fgVariableNames[kVtxCovYZ] = "Vtx covYZ";
  fgVariableUnits[kVtxCovYZ] = "cm";
  fgVariableNames[kVtxCovZZ] = "Vtx covZZ";
  fgVariableUnits[kVtxCovZZ] = "cm";
  fgVariableNames[kVtxChi2] = "Vtx chi2";
  fgVariableUnits[kVtxChi2] = "";
  fgVariableNames[kCentVZERO] = "Centrality VZERO";
  fgVariableUnits[kCentVZERO] = "%";
  fgVariableNames[kMCEventGeneratorId] = "MC Generator ID";
  fgVariableNames[kMCVtxX] = "MC Vtx X";
  fgVariableNames[kMCVtxY] = "MC Vtx Y";
  fgVariableNames[kMCVtxZ] = "MC Vtx Z";
  fgVariableNames[kMCEventTime] = "MC event time";
  fgVariableNames[kMCEventWeight] = "MC event weight";
  fgVariableNames[kMCEventImpParam] = "MC impact parameter";
  fgVariableUnits[kMCEventGeneratorId] = "";
  fgVariableUnits[kMCVtxX] = "cm";
  fgVariableUnits[kMCVtxY] = "cm";
  fgVariableUnits[kMCVtxZ] = "cm";
  fgVariableUnits[kMCEventTime] = ""; // TODO: add proper unit
  fgVariableUnits[kMCEventWeight] = "";
  fgVariableUnits[kMCEventImpParam] = "b";
  fgVariableNames[kPt] = "p_{T}";
  fgVariableUnits[kPt] = "GeV/c";
  fgVariableNames[kP] = "p";
  fgVariableUnits[kP] = "GeV/c";
  fgVariableNames[kPx] = "p_{x}";
  fgVariableUnits[kPy] = "GeV/c";
  fgVariableNames[kPy] = "p_{y}";
  fgVariableUnits[kPz] = "GeV/c";
  fgVariableNames[kPz] = "p_{z}";
  fgVariableUnits[kPx] = "GeV/c";
  fgVariableNames[kEta] = "#eta";
  fgVariableUnits[kEta] = "";
  fgVariableNames[kPhi] = "#varphi";
  fgVariableUnits[kPhi] = "rad.";
  fgVariableNames[kRap] = "y";
  fgVariableUnits[kRap] = "";
  fgVariableNames[kMass] = "mass";
  fgVariableUnits[kMass] = "GeV/c2";
  fgVariableNames[kCharge] = "charge";
  fgVariableUnits[kCharge] = "";
  fgVariableNames[kPin] = "p_{IN}";
  fgVariableUnits[kPin] = "GeV/c";
  fgVariableNames[kITSncls] = "ITS #cls";
  fgVariableUnits[kITSncls] = "";
  fgVariableUnits[kITSClusterMap] = "ITSClusterMap";
  fgVariableNames[kITSchi2] = "ITS chi2";
  fgVariableUnits[kITSchi2] = "";
  fgVariableNames[kITSlayerHit] = "ITS layer";
  fgVariableUnits[kITSlayerHit] = "";
  fgVariableNames[kTPCncls] = "TPC #cls";
  fgVariableUnits[kTPCncls] = "";
  fgVariableNames[kTPCnclsCR] = "TPC #cls crossed rows";
  fgVariableUnits[kTPCnclsCR] = "";
  fgVariableNames[kTPCchi2] = "TPC chi2";
  fgVariableUnits[kTPCchi2] = "";
  fgVariableNames[kTPCsignal] = "TPC dE/dx";
  fgVariableUnits[kTPCsignal] = "";
  fgVariableNames[kTRDsignal] = "TRD dE/dx";
  fgVariableUnits[kTRDsignal] = "";
  fgVariableNames[kTOFbeta] = "TOF #beta";
  fgVariableUnits[kTOFbeta] = "";
  fgVariableNames[kTrackLength] = "track length";
  fgVariableUnits[kTrackLength] = "cm";
  fgVariableNames[kTrackDCAxy] = "DCA_{xy}";
  fgVariableUnits[kTrackDCAxy] = "cm";
  fgVariableNames[kTrackDCAz] = "DCA_{z}";
  fgVariableUnits[kTrackDCAz] = "cm";
  fgVariableNames[kTPCnSigmaEl] = "n #sigma_{e}^{TPC}";
  fgVariableUnits[kTPCnSigmaEl] = "";
  fgVariableNames[kTPCnSigmaEl_Corr] = "n #sigma_{e}^{TPC} Corr.";
  fgVariableUnits[kTPCnSigmaEl_Corr] = "";
  fgVariableNames[kTPCnSigmaMu] = "n #sigma_{#mu}^{TPC}";
  fgVariableUnits[kTPCnSigmaMu] = "";
  fgVariableNames[kTPCnSigmaPi] = "n #sigma_{#pi}^{TPC}";
  fgVariableUnits[kTPCnSigmaPi] = "";
  fgVariableNames[kTPCnSigmaPi_Corr] = "n #sigma_{#pi}^{TPC} Corr.";
  fgVariableUnits[kTPCnSigmaPi_Corr] = "";
  fgVariableNames[kTPCnSigmaKa] = "n #sigma_{K}^{TPC}";
  fgVariableUnits[kTPCnSigmaKa] = "";
  fgVariableNames[kTPCnSigmaPr] = "n #sigma_{p}^{TPC}";
  fgVariableUnits[kTPCnSigmaPr] = "";
  fgVariableNames[kTPCnSigmaPr_Corr] = "n #sigma_{p}^{TPC} Corr.";
  fgVariableUnits[kTPCnSigmaPr_Corr] = "";
  fgVariableNames[kTOFnSigmaEl] = "n #sigma_{e}^{TOF}";
  fgVariableUnits[kTOFnSigmaEl] = "";
  fgVariableNames[kTOFnSigmaMu] = "n #sigma_{#mu}^{TOF}";
  fgVariableUnits[kTOFnSigmaMu] = "";
  fgVariableNames[kTOFnSigmaPi] = "n #sigma_{#pi}^{TOF}";
  fgVariableUnits[kTOFnSigmaPi] = "";
  fgVariableNames[kTOFnSigmaKa] = "n #sigma_{K}^{TOF}";
  fgVariableUnits[kTOFnSigmaKa] = "";
  fgVariableNames[kTOFnSigmaPr] = "n #sigma_{p}^{TOF}";
  fgVariableUnits[kTOFnSigmaPr] = "";
  fgVariableNames[kIsLegFromGamma] = "is leg from #gamma #rightarror e^{+}e^{-}";
  fgVariableUnits[kIsLegFromGamma] = "";
  fgVariableNames[kIsLegFromK0S] = "is leg from K_{S}^{0} #rightarror #pi^{+}#pi^{-}";
  fgVariableUnits[kIsLegFromK0S] = "";
  fgVariableNames[kIsLegFromLambda] = "is leg from #Lambda #rightarror p#pi^{-}";
  fgVariableUnits[kIsLegFromLambda] = "";
  fgVariableNames[kIsLegFromAntiLambda] = "is leg from #bar{#Lambda} #rightarrow #bar{p}#pi^{+}";
  fgVariableUnits[kIsLegFromAntiLambda] = "";
  fgVariableNames[kIsLegFromOmega] = "is leg from #Omega^{#mp} #rightarrow #LambdaKi^{#pm}";
  fgVariableUnits[kIsLegFromOmega] = "";
  fgVariableNames[kMuonNClusters] = "muon n-clusters";
  fgVariableUnits[kMuonNClusters] = "";
  fgVariableNames[kMuonRAtAbsorberEnd] = "R at the end of the absorber";
  fgVariableUnits[kMuonRAtAbsorberEnd] = "cm";
  fgVariableNames[kMuonPDca] = "p x dca";
  fgVariableUnits[kMuonPDca] = "cm x GeV/c";
  fgVariableNames[kMCHBitMap] = "MCH bitmap";
  fgVariableUnits[kMCHBitMap] = "";
  fgVariableNames[kMuonChi2] = "#chi^{2}";
  fgVariableUnits[kMuonChi2] = "";
  fgVariableNames[kMuonChi2MatchMCHMID] = "#chi^{2} MCH-MID";
  fgVariableUnits[kMuonChi2MatchMCHMID] = "";
  fgVariableNames[kMuonChi2MatchMCHMFT] = "#chi^{2} MCH-MFT";
  fgVariableUnits[kMuonChi2MatchMCHMFT] = "";
  fgVariableNames[kMuonMatchScoreMCHMFT] = "match score MCH-MFT";
  fgVariableUnits[kMuonMatchScoreMCHMFT] = "";
  fgVariableNames[kMuonDCAx] = "dca_X";
  fgVariableUnits[kMuonDCAx] = "cm";
  fgVariableNames[kMuonDCAy] = "dca_Y";
  fgVariableUnits[kMuonDCAy] = "cm";
  fgVariableNames[kMuonCXX] = "cov XX";
  fgVariableUnits[kMuonCXX] = "";
  fgVariableNames[kMuonCYY] = "cov YY";
  fgVariableUnits[kMuonCYY] = "";
  fgVariableNames[kMuonCPhiPhi] = "cov PhiPhi";
  fgVariableUnits[kMuonCPhiPhi] = "";
  fgVariableNames[kMuonCTglTgl] = "cov TglTgl";
  fgVariableUnits[kMuonCTglTgl] = "";
  fgVariableNames[kMuonC1Pt21Pt2] = "cov 1Pt1Pt";
  fgVariableUnits[kMuonC1Pt21Pt2] = "";
  fgVariableNames[kMCPdgCode] = "MC PDG code";
  fgVariableUnits[kMCPdgCode] = "";
  fgVariableNames[kMCParticleWeight] = "MC particle weight";
  fgVariableUnits[kMCParticleWeight] = "";
  fgVariableNames[kMCPx] = "MC px";
  fgVariableUnits[kMCPx] = "GeV/c";
  fgVariableNames[kMCPy] = "MC py";
  fgVariableUnits[kMCPy] = "GeV/c";
  fgVariableNames[kMCPz] = "MC pz";
  fgVariableUnits[kMCPz] = "GeV/c";
  fgVariableNames[kMCPt] = "MC pt";
  fgVariableUnits[kMCPt] = "GeV/c";
  fgVariableNames[kMCPhi] = "#varphi";
  fgVariableUnits[kMCPhi] = "rad";
  fgVariableNames[kMCEta] = "MC #eta";
  fgVariableUnits[kMCEta] = "";
  fgVariableNames[kMCY] = "MC y";
  fgVariableUnits[kMCY] = "";
  fgVariableNames[kMCE] = "MC Energy";
  fgVariableUnits[kMCE] = "GeV";
  fgVariableNames[kMCVx] = "MC vx";
  fgVariableUnits[kMCVx] = "cm"; // TODO: check the unit
  fgVariableNames[kMCVy] = "MC vy";
  fgVariableUnits[kMCVy] = "cm"; // TODO: check the unit
  fgVariableNames[kMCVz] = "MC vz";
  fgVariableUnits[kMCVz] = "cm"; // TODO: check the unit
  fgVariableNames[kCandidateId] = "";
  fgVariableUnits[kCandidateId] = "";
  fgVariableNames[kPairType] = "Pair type";
  fgVariableUnits[kPairType] = "";
  fgVariableNames[kVertexingLxy] = "Pair Lxy";
  fgVariableUnits[kVertexingLxy] = "cm";
  fgVariableNames[kVertexingLz] = "Pair Lz";
  fgVariableUnits[kVertexingLz] = "cm";
  fgVariableNames[kVertexingLxyz] = "Pair Lxyz";
  fgVariableUnits[kVertexingLxyz] = "cm";
  fgVariableNames[kVertexingLxyErr] = "Pair Lxy err.";
  fgVariableUnits[kVertexingLxyErr] = "cm";
  fgVariableNames[kVertexingLzErr] = "Pair Lz err.";
  fgVariableUnits[kVertexingLzErr] = "cm";
  fgVariableNames[kVertexingLxyzErr] = "Pair Lxyz err.";
  fgVariableUnits[kVertexingLxyzErr] = "cm";
  fgVariableNames[kVertexingTauz] = "Pair pseudo-proper Tauz";
  fgVariableUnits[kVertexingTauz] = "ns";
  fgVariableNames[kVertexingTauxy] = "Pair pseudo-proper Tauxy";
  fgVariableUnits[kVertexingTauxy] = "ns";
  fgVariableNames[kVertexingTauzErr] = "Pair pseudo-proper Tauz err.";
  fgVariableUnits[kVertexingTauzErr] = "ns";
  fgVariableNames[kVertexingTauxyErr] = "Pair pseudo-proper Tauxy err.";
  fgVariableUnits[kVertexingTauxyErr] = "ns";
  fgVariableNames[kVertexingProcCode] = "DCAFitterN<2> processing code";
  fgVariableUnits[kVertexingProcCode] = "";
  fgVariableNames[kVertexingChi2PCA] = "Pair #chi^{2} at PCA";
  fgVariableUnits[kVertexingChi2PCA] = "";
  fgVariableNames[kQ2X0A] = "Q_{2,x}^{A} ";
  fgVariableUnits[kQ2X0A] = "";
  fgVariableNames[kQ2Y0A] = "Q_{2,y}^{A} ";
  fgVariableUnits[kQ2Y0A] = "";
  fgVariableNames[kQ2X0B] = "Q_{2,x}^{B} ";
  fgVariableUnits[kQ2X0B] = "";
  fgVariableNames[kQ2Y0B] = "Q_{2,y}^{B} ";
  fgVariableUnits[kQ2Y0B] = "";
  fgVariableNames[kQ2X0C] = "Q_{2,x}^{C} ";
  fgVariableUnits[kQ2X0C] = "";
  fgVariableNames[kQ2Y0C] = "Q_{2,y}^{C} ";
  fgVariableUnits[kQ2Y0C] = "";
  fgVariableNames[kMultA] = "N_{ch}^{A} ";
  fgVariableUnits[kMultA] = "";
  fgVariableNames[kMultB] = "N_{ch}^{B} ";
  fgVariableUnits[kMultB] = "";
  fgVariableNames[kMultC] = "N_{ch}^{C} ";
  fgVariableUnits[kMultC] = "";
  fgVariableNames[kQ3X0A] = "Q_{3,x}^{A} ";
  fgVariableUnits[kQ3X0A] = "";
  fgVariableNames[kQ3Y0A] = "Q_{3,y}^{A} ";
  fgVariableUnits[kQ3Y0A] = "";
  fgVariableNames[kQ3X0B] = "Q_{3,x}^{B} ";
  fgVariableUnits[kQ3X0B] = "";
  fgVariableNames[kQ3Y0B] = "Q_{3,y}^{B} ";
  fgVariableUnits[kQ3Y0B] = "";
  fgVariableNames[kQ3X0C] = "Q_{3,x}^{C} ";
  fgVariableUnits[kQ3X0C] = "";
  fgVariableNames[kQ3Y0C] = "Q_{3,y}^{C} ";
  fgVariableUnits[kQ3Y0C] = "";
  fgVariableNames[kU2Q2] = "u_{2}Q_{2}^{A} ";
  fgVariableUnits[kU2Q2] = "";
  fgVariableNames[kU3Q3] = "u_{3}Q_{3}^{A} ";
  fgVariableUnits[kU3Q3] = "";
  fgVariableNames[kCos2DeltaPhi] = "cos 2(#varphi-#Psi_{2}^{A}) ";
  fgVariableUnits[kCos2DeltaPhi] = "";
  fgVariableNames[kCos3DeltaPhi] = "cos 3(#varphi-#Psi_{3}^{A}) ";
  fgVariableUnits[kCos3DeltaPhi] = "";
  fgVariableNames[kR2SP] = "R_{2}^{SP} ";
  fgVariableUnits[kR2SP] = "";
  fgVariableNames[kR3SP] = "R_{3}^{SP} ";
  fgVariableUnits[kR3SP] = "";
  fgVariableNames[kR2EP] = "R_{2}^{EP} ";
  fgVariableUnits[kR2EP] = "";
  fgVariableNames[kR3EP] = "R_{3}^{EP} ";
  fgVariableUnits[kR3EP] = "";
  fgVariableNames[kPairMass] = "mass";
  fgVariableUnits[kPairMass] = "GeV/c2";
  fgVariableNames[kPairPt] = "p_{T}";
  fgVariableUnits[kPairPt] = "GeV/c";
  fgVariableNames[kPairEta] = "#eta";
  fgVariableUnits[kPairEta] = "";
  fgVariableNames[kPairPhi] = "#varphi";
  fgVariableUnits[kPairPhi] = "rad.";
  fgVariableNames[kPairPhiv] = "#varphi_{V}";
  fgVariableUnits[kPairPhiv] = "rad.";
  fgVariableNames[kDeltaEta] = "#Delta#eta";
  fgVariableUnits[kDeltaEta] = "";
  fgVariableNames[kDeltaPhi] = "#Delta#phi";
  fgVariableUnits[kDeltaPhi] = "rad.";
  fgVariableNames[kDeltaPhiSym] = "#Delta#phi";
  fgVariableUnits[kDeltaPhiSym] = "rad.";
  fgVariableNames[kCosThetaHE] = "cos#it{#theta}";
  fgVariableUnits[kCosThetaHE] = "";
  fgVariableNames[kPsiPair] = "#Psi_{pair}";
  fgVariableUnits[kPsiPair] = "rad.";
  fgVariableNames[kDeltaPhiPair] = "#Delta#phi";
  fgVariableUnits[kDeltaPhiPair] = "rad.";
  fgVariableNames[kQuadDCAabsXY] = "DCA_{xy}^{quad}";
  fgVariableUnits[kQuadDCAabsXY] = "cm";
  fgVariableNames[kQuadDCAsigXY] = "DCA_{xy}^{quad}";
  fgVariableUnits[kQuadDCAsigXY] = "#sigma";
  fgVariableNames[kTrackDCAsigXY] = "DCA_{xy}";
  fgVariableUnits[kTrackDCAsigXY] = "#sigma";
  fgVariableNames[kTrackDCAsigZ] = "DCA_{z}";
  fgVariableUnits[kTrackDCAsigZ] = "#sigma";
  fgVariableNames[kTrackDCAresXY] = "#Delta DCA_{xy}";
  fgVariableUnits[kTrackDCAresXY] = "cm";
  fgVariableNames[kTrackDCAresZ] = "#Delta DCA_{z}";
  fgVariableUnits[kTrackDCAresZ] = "cm";
  fgVariableNames[kBitMapIndex] = " ";
  fgVariableUnits[kBitMapIndex] = "";
}
