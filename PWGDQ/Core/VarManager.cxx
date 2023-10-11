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
#include <cmath>
#include "PWGDQ/Core/VarManager.h"
#include "Tools/KFparticle/KFUtilities.h"

using std::cout;
using std::endl;
using namespace o2::constants::physics;

ClassImp(VarManager);

TString VarManager::fgVariableNames[VarManager::kNVars] = {""};
TString VarManager::fgVariableUnits[VarManager::kNVars] = {""};
bool VarManager::fgUsedVars[VarManager::kNVars] = {false};
bool VarManager::fgUsedKF = false;
float VarManager::fgMagField = 0.5;
float VarManager::fgValues[VarManager::kNVars] = {0.0f};
std::map<int, int> VarManager::fgRunMap;
TString VarManager::fgRunStr = "";
std::vector<int> VarManager::fgRunList = {0};
float VarManager::fgCenterOfMassEnergy = 13600;         // GeV
float VarManager::fgMassofCollidingParticle = 9.382720; // GeV
o2::vertexing::DCAFitterN<2> VarManager::fgFitterTwoProngBarrel;
o2::vertexing::DCAFitterN<3> VarManager::fgFitterThreeProngBarrel;
o2::vertexing::FwdDCAFitterN<2> VarManager::fgFitterTwoProngFwd;
o2::vertexing::FwdDCAFitterN<3> VarManager::fgFitterThreeProngFwd;
std::map<VarManager::CalibObjects, TObject*> VarManager::fgCalibs;
bool VarManager::fgRunTPCPostCalibration[4] = {false, false, false, false};

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

  if (fgUsedVars[kVertexingLxyOverErr]) {
    fgUsedVars[kVertexingLxy] = kTRUE;
    fgUsedVars[kVertexingLxyErr] = kTRUE;
  }
  if (fgUsedVars[kVertexingLzOverErr]) {
    fgUsedVars[kVertexingLz] = kTRUE;
    fgUsedVars[kVertexingLzErr] = kTRUE;
  }
  if (fgUsedVars[kVertexingLxyzOverErr]) {
    fgUsedVars[kVertexingLxyz] = kTRUE;
    fgUsedVars[kVertexingLxyzErr] = kTRUE;
  }
  if (fgUsedVars[kKFTracksDCAxyzMax]) {
    fgUsedVars[kKFTrack0DCAxyz] = kTRUE;
    fgUsedVars[kKFTrack1DCAxyz] = kTRUE;
  }
  if (fgUsedVars[kKFTracksDCAxyMax]) {
    fgUsedVars[kKFTrack0DCAxy] = kTRUE;
    fgUsedVars[kKFTrack1DCAxy] = kTRUE;
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
  fgRunList = runs;
}

//__________________________________________________________________
void VarManager::SetRunlist(TString period)
{
  //
  // runlist for the different periods
  // TODO: add more periods
  std::vector<int> LHC22f = {520259, 520294, 520471, 520472, 520473};
  std::vector<int> LHC22m = {523141, 523142, 523148, 523182, 523186, 523298, 523306, 523308, 523309, 523397, 523399, 523401, 523441, 523541, 523559, 523669, 523671, 523677, 523728, 523731, 523779, 523783, 523786, 523788, 523789, 523792, 523797, 523821, 523897};
  std::vector<int> LHC22o = {526463, 526465, 526466, 526467, 526468, 526486, 526505, 526508, 526510, 526512, 526525, 526526, 526528, 526534, 526559, 526596, 526606, 526612, 526638, 526639, 526641, 526643, 526647, 526649, 526689, 526712, 526713, 526714, 526715, 526776, 526865, 526886, 526926, 526927, 526928, 526929, 526934, 526935, 526937, 526966, 526968, 527015, 527109, 527228, 527237, 527240, 527259, 527260, 527261, 527262, 527345, 527347, 527349, 527446, 527518, 527522, 527523, 527799, 527821, 527825, 527826, 527828, 527848, 527850, 527852, 527863, 527864, 527865, 527869, 527871, 527895, 527898, 527899, 527902, 527963, 527976, 527978, 527979, 528021, 528026, 528036, 528093, 528094, 528097, 528105, 528107, 528109, 528110, 528231, 528232, 528233, 528263, 528266, 528292, 528294, 528316, 528319, 528328, 528329, 528330, 528332, 528336, 528347, 528359, 528379, 528381, 528386, 528448, 528451, 528461, 528463, 528529, 528530, 528531, 528534, 528537, 528543};
  std::vector<int> LHC22p = {528602, 528604, 528617, 528781, 528782, 528783, 528784, 528798, 528801};
  std::vector<int> LHC22q = {528991, 528997, 529003, 529005, 529006, 529009, 529015, 529035, 529037, 529038, 529039, 529043};
  std::vector<int> LHC22r = {529066, 529067, 529077, 529078, 529084, 529088, 529115, 529116, 529117, 529128, 529129, 529208, 529209, 529210, 529211, 529235, 529237, 529242, 529248, 529252, 529270, 529306, 529310, 529317, 529320, 529324, 529337, 529338, 529341};
  std::vector<int> LHC22s = {529397, 529399, 529403, 529414, 529418};
  std::vector<int> LHC22t = {529450, 529452, 529454, 529458, 529460, 529461, 529462, 529542, 529552, 529554, 529610, 529662, 529663, 529664, 529674, 529675, 529690, 529691};
  if (period.Contains("LHC22f")) {
    SetRunNumbers(LHC22f);
  }
  if (period.Contains("LHC22q")) {
    SetRunNumbers(LHC22q);
  }
  if (period.Contains("LHC22m")) {
    SetRunNumbers(LHC22m);
  }
  if (period.Contains("LHC22o")) {
    SetRunNumbers(LHC22o);
  }
  if (period.Contains("LHC22p")) {
    SetRunNumbers(LHC22p);
  }
  if (period.Contains("LHC22r")) {
    SetRunNumbers(LHC22r);
  }
  if (period.Contains("LHC22s")) {
    SetRunNumbers(LHC22s);
  }
  if (period.Contains("LHC22t")) {
    SetRunNumbers(LHC22t);
  }
}
//__________________________________________________________________
void VarManager::SetDummyRunlist(int InitRunnumber)
{
  //
  // runlist for the different periods
  fgRunList.clear();
  fgRunList.push_back(InitRunnumber);
  fgRunList.push_back(InitRunnumber + 100);
}

//__________________________________________________________________
int VarManager::GetDummyFirst()
{
  //
  // Get the fist index of the vector of run numbers
  //
  return fgRunList[0];
}
//__________________________________________________________________
int VarManager::GetDummyLast()
{
  //
  // Get the last index of the vector of run numbers
  //
  return fgRunList[fgRunList.size() - 1];
}
//_________________________________________________________________
float VarManager::GetRunIndex(double Runnumber)
{
  //
  // Get the index of RunNumber in it's runlist
  //
  int runNumber = static_cast<int>(Runnumber);
  auto runIndex = std::find(fgRunList.begin(), fgRunList.end(), runNumber);
  float index = std::distance(fgRunList.begin(), runIndex);
  return index;
}
//__________________________________________________________________
void VarManager::SetCollisionSystem(TString system, float energy)
{
  //
  // Set the collision system and the center of mass energy
  //
  fgCenterOfMassEnergy = energy;

  if (system.Contains("PbPb")) {
    fgMassofCollidingParticle = MassProton * 208;
  }
  if (system.Contains("pp")) {
    fgMassofCollidingParticle = MassProton;
  }
  // TO Do: add more systems
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

//__________________________________________________________________
float VarManager::calculateCosPA(KFParticle kfp, KFParticle PV)
{
  return cpaFromKF(kfp, PV);
}
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
  fgVariableNames[kVtxNcontribReal] = "Real Vtx contrib.";
  fgVariableUnits[kVtxNcontribReal] = "";
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
  fgVariableNames[kCentFT0C] = "Centrality FT0C";
  fgVariableUnits[kCentFT0C] = "%";
  fgVariableNames[kMultTPC] = "Multiplicity TPC";
  fgVariableUnits[kMultTPC] = "";
  fgVariableNames[kMultFV0A] = "Multiplicity FV0A";
  fgVariableUnits[kMultFV0A] = "";
  fgVariableNames[kMultFV0C] = "Multiplicity FV0C";
  fgVariableUnits[kMultFV0C] = "";
  fgVariableNames[kMultFT0A] = "Multiplicity FT0A";
  fgVariableUnits[kMultFT0A] = "";
  fgVariableNames[kMultFT0C] = "Multiplicity FT0C";
  fgVariableUnits[kMultFT0C] = "";
  fgVariableNames[kMultFDDA] = "Multiplicity FDDA";
  fgVariableUnits[kMultFDDA] = "";
  fgVariableNames[kMultFDDC] = "Multiplicity FDDC";
  fgVariableUnits[kMultFDDC] = "";
  fgVariableNames[kMultZNA] = "Multiplicity ZNA";
  fgVariableUnits[kMultZNA] = "";
  fgVariableNames[kMultZNC] = "Multiplicity ZNC";
  fgVariableUnits[kMultZNC] = "";
  fgVariableNames[kMultTracklets] = "Multiplicity Tracklets";
  fgVariableUnits[kMultTracklets] = "";
  fgVariableNames[kCentFT0C] = "Centrality FT0C";
  fgVariableUnits[kCentFT0C] = "%";
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
  fgVariableNames[kInvPt] = "1/p_{T}";
  fgVariableUnits[kInvPt] = "1/(GeV/c)";
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
  fgVariableNames[kDeltaPtotTracks] = "#it{p}_{Tot}^{#mu+} - #it{p}_{Tot}^{#mu-}";
  fgVariableUnits[kDeltaPtotTracks] = "GeV/c";
  fgVariableNames[kCharge] = "charge";
  fgVariableUnits[kCharge] = "";
  fgVariableNames[kPin] = "p_{IN}";
  fgVariableUnits[kPin] = "GeV/c";
  fgVariableNames[kSignedPin] = "p_{IN} x charge";
  fgVariableUnits[kSignedPin] = "GeV/c";
  fgVariableNames[kTOFExpMom] = "TOF expected momentum";
  fgVariableUnits[kTOFExpMom] = "GeV/c";
  fgVariableNames[kTrackTime] = "Track time wrt collision().bc()";
  fgVariableUnits[kTrackTime] = "ns";
  fgVariableNames[kTrackTimeRes] = "Resolution of the track time";
  fgVariableUnits[kTrackTimeRes] = "ns";
  fgVariableNames[kTrackTimeResRelative] = "Relative resolution of the track time";
  fgVariableUnits[kTrackTimeResRelative] = "";
  fgVariableNames[kDetectorMap] = "DetectorMap";
  fgVariableUnits[kDetectorMap] = "";
  fgVariableNames[kHasITS] = "HasITS";
  fgVariableUnits[kHasITS] = "";
  fgVariableNames[kHasTRD] = "HasTRD";
  fgVariableUnits[kHasTRD] = "";
  fgVariableNames[kHasTOF] = "HasTOF";
  fgVariableUnits[kHasTOF] = "";
  fgVariableNames[kHasTPC] = "HasTPC";
  fgVariableUnits[kHasTPC] = "";
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
  fgVariableNames[kTPCnSigmaKa_Corr] = "n #sigma_{K}^{TPC} Corr.";
  fgVariableUnits[kTPCnSigmaKa_Corr] = "";
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
  fgVariableNames[kIsAmbiguous] = "is ambiguous track";
  fgVariableUnits[kIsAmbiguous] = "";
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
  fgVariableNames[kMuonTime] = "Track time wrt collision().bc()";
  fgVariableUnits[kMuonTime] = "ns";
  fgVariableNames[kMuonTimeRes] = "Resolution of the track time";
  fgVariableUnits[kMuonTimeRes] = "ns";
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
  fgVariableNames[kVertexingLxyOverErr] = "Pair Lxy/DLxy";
  fgVariableUnits[kVertexingLxyOverErr] = "";
  fgVariableNames[kVertexingLzOverErr] = "Pair Lz/DLz";
  fgVariableUnits[kVertexingLzOverErr] = "";
  fgVariableNames[kVertexingLxyzOverErr] = "Pair Lxyz/DLxyz";
  fgVariableUnits[kVertexingLxyzOverErr] = "";
  fgVariableNames[kKFTrack0DCAxyz] = "Daughter0 DCAxyz";
  fgVariableUnits[kKFTrack0DCAxyz] = "cm";
  fgVariableNames[kKFTrack1DCAxyz] = "Daughter1 DCAxyz";
  fgVariableUnits[kKFTrack1DCAxyz] = "cm";
  fgVariableNames[kKFTracksDCAxyzMax] = "Maximum DCAxyz of two daughters";
  fgVariableUnits[kKFTracksDCAxyzMax] = "cm";
  fgVariableNames[kKFDCAxyzBetweenProngs] = "DCAxyz between two daughters";
  fgVariableUnits[kKFDCAxyzBetweenProngs] = "cm";
  fgVariableNames[kKFTrack0DCAxy] = "Daughter0 DCAxy";
  fgVariableUnits[kKFTrack0DCAxy] = "cm";
  fgVariableNames[kKFTrack1DCAxy] = "Daughter1 DCAxy";
  fgVariableUnits[kKFTrack1DCAxy] = "cm";
  fgVariableNames[kKFTracksDCAxyMax] = "Maximum DCAxy of two daughters";
  fgVariableUnits[kKFTracksDCAxyMax] = "cm";
  fgVariableNames[kKFDCAxyBetweenProngs] = "DCAxy between two daughters";
  fgVariableUnits[kKFDCAxyBetweenProngs] = "cm";
  fgVariableNames[kKFChi2OverNDFGeo] = "Pair geometrical #chi^{2}/ndf";
  fgVariableUnits[kKFChi2OverNDFGeo] = "";
  fgVariableNames[kKFCosPA] = "cosPA";
  fgVariableUnits[kKFCosPA] = "";
  fgVariableNames[kKFNContributorsPV] = "Real Number of Trks to PV";
  fgVariableUnits[kKFNContributorsPV] = "";
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
  fgVariableNames[kPairMassDau] = "mass dilepton";
  fgVariableUnits[kPairMassDau] = "GeV/c2";
  fgVariableNames[kMassDau] = "mass HF";
  fgVariableUnits[kMassDau] = "GeV/c2";
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
  fgVariableNames[kPhiHE] = "#varphi_{HE}";
  fgVariableUnits[kPhiHE] = "rad.";
  fgVariableNames[kCosThetaCS] = "cos#it{#theta}_{CS}";
  fgVariableUnits[kCosThetaCS] = "";
  fgVariableNames[kPhiCS] = "#varphi_{CS}";
  fgVariableUnits[kPhiCS] = "rad.";
  fgVariableNames[kPsiPair] = "#Psi_{pair}";
  fgVariableUnits[kPsiPair] = "rad.";
  fgVariableNames[kDeltaPhiPair] = "#Delta#phi";
  fgVariableUnits[kDeltaPhiPair] = "rad.";
  fgVariableNames[kOpeningAngle] = "Opening angle";
  fgVariableUnits[kOpeningAngle] = "rad.";
  fgVariableNames[kQuadDCAabsXY] = "DCA_{xy}^{quad}";
  fgVariableUnits[kQuadDCAabsXY] = "cm";
  fgVariableNames[kQuadDCAsigXY] = "DCA_{xy}^{quad}";
  fgVariableUnits[kQuadDCAsigXY] = "#sigma";
  fgVariableNames[kQuadDCAabsZ] = "DCA_{z}^{quad}";
  fgVariableUnits[kQuadDCAabsZ] = "cm";
  fgVariableNames[kQuadDCAsigZ] = "DCA_{z}^{quad}";
  fgVariableUnits[kQuadDCAsigZ] = "#sigma";
  fgVariableNames[kQuadDCAsigXYZ] = "DCA_{xyz}^{quad}";
  fgVariableUnits[kQuadDCAsigXYZ] = "#sigma";
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
  fgVariableNames[kMassCharmHadron] = "mass (charm hadron)";
  fgVariableUnits[kMassCharmHadron] = "GeV/c2";
  fgVariableNames[kPtCharmHadron] = "p_{T} (charm hadron)";
  fgVariableUnits[kPtCharmHadron] = "GeV/c";
  fgVariableNames[kRapCharmHadron] = "y (charm hadron)";
  fgVariableUnits[kRapCharmHadron] = " ";
  fgVariableNames[kPhiCharmHadron] = "#varphi (charm hadron)";
  fgVariableUnits[kPhiCharmHadron] = "rad.";
}
