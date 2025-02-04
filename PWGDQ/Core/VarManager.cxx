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
#include <iostream>
#include <vector>
#include <map>
#include "PWGDQ/Core/VarManager.h"
#include "Tools/KFparticle/KFUtilities.h"

using std::cout;
using std::endl;
using namespace o2::constants::physics;

ClassImp(VarManager);

TString VarManager::fgVariableNames[VarManager::kNVars] = {""};
TString VarManager::fgVariableUnits[VarManager::kNVars] = {""};
std::map<TString, int> VarManager::fgVarNamesMap;
bool VarManager::fgUsedVars[VarManager::kNVars] = {false};
bool VarManager::fgUsedKF = false;
float VarManager::fgMagField = 0.5;
float VarManager::fgValues[VarManager::kNVars] = {0.0f};
std::map<int, int> VarManager::fgRunMap;
TString VarManager::fgRunStr = "";
std::vector<int> VarManager::fgRunList = {0};
float VarManager::fgCenterOfMassEnergy = 13600;         // GeV
float VarManager::fgMassofCollidingParticle = 9.382720; // GeV
float VarManager::fgTPCInterSectorBoundary = 1.0;       // cm
int VarManager::fgITSROFbias = 0;
int VarManager::fgITSROFlength = 100;
int VarManager::fgITSROFBorderMarginLow = 0;
int VarManager::fgITSROFBorderMarginHigh = 0;
uint64_t VarManager::fgSOR = 0;
uint64_t VarManager::fgEOR = 0;
o2::vertexing::DCAFitterN<2> VarManager::fgFitterTwoProngBarrel;
o2::vertexing::DCAFitterN<3> VarManager::fgFitterThreeProngBarrel;
o2::vertexing::FwdDCAFitterN<2> VarManager::fgFitterTwoProngFwd;
o2::vertexing::FwdDCAFitterN<3> VarManager::fgFitterThreeProngFwd;
o2::globaltracking::MatchGlobalFwd VarManager::mMatching;
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
    fgUsedVars[kPt] = true;
    fgUsedVars[kEta] = true;
  }

  if (fgUsedVars[kVertexingLxyOverErr]) {
    fgUsedVars[kVertexingLxy] = true;
    fgUsedVars[kVertexingLxyErr] = true;
  }
  if (fgUsedVars[kVertexingLzOverErr]) {
    fgUsedVars[kVertexingLz] = true;
    fgUsedVars[kVertexingLzErr] = true;
  }
  if (fgUsedVars[kVertexingLxyzOverErr]) {
    fgUsedVars[kVertexingLxyz] = true;
    fgUsedVars[kVertexingLxyzErr] = true;
  }
  if (fgUsedVars[kKFTracksDCAxyzMax]) {
    fgUsedVars[kKFTrack0DCAxyz] = true;
    fgUsedVars[kKFTrack1DCAxyz] = true;
  }
  if (fgUsedVars[kKFTracksDCAxyMax]) {
    fgUsedVars[kKFTrack0DCAxy] = true;
    fgUsedVars[kKFTrack1DCAxy] = true;
  }
  if (fgUsedVars[kTrackIsInsideTPCModule]) {
    fgUsedVars[kPhiTPCOuter] = true;
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
  fgVariableNames[kTimeFromSOR] = "time since SOR";
  fgVariableUnits[kTimeFromSOR] = "min.";
  fgVariableNames[kBCOrbit] = "Bunch crossing";
  fgVariableUnits[kBCOrbit] = "";
  fgVariableNames[kIsPhysicsSelection] = "Physics selection";
  fgVariableUnits[kIsPhysicsSelection] = "";
  fgVariableNames[kVtxX] = "Vtx X ";
  fgVariableUnits[kVtxX] = "cm";
  fgVariableNames[kVtxY] = "Vtx Y ";
  fgVariableUnits[kVtxY] = "cm";
  fgVariableNames[kVtxZ] = "Vtx Z ";
  fgVariableUnits[kVtxZ] = "cm";
  fgVariableNames[kCollisionTime] = "collision time wrt BC";
  fgVariableUnits[kCollisionTime] = "ns";
  fgVariableNames[kCollisionTimeRes] = "collision time resolution";
  fgVariableUnits[kCollisionTimeRes] = "ns";
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
  fgVariableNames[kMultDimuons] = "Multiplicity Dimuons Unlike Sign";
  fgVariableUnits[kMultDimuons] = "";
  fgVariableNames[kCentFT0C] = "Centrality FT0C";
  fgVariableUnits[kCentFT0C] = "%";
  fgVariableNames[kMCEventGeneratorId] = "MC Generator ID";
  fgVariableNames[kMCEventSubGeneratorId] = "MC SubGenerator ID";
  fgVariableNames[kMCVtxX] = "MC Vtx X";
  fgVariableNames[kMCVtxY] = "MC Vtx Y";
  fgVariableNames[kMCVtxZ] = "MC Vtx Z";
  fgVariableNames[kMCEventTime] = "MC event time";
  fgVariableNames[kMCEventWeight] = "MC event weight";
  fgVariableNames[kMCEventImpParam] = "MC impact parameter";
  fgVariableUnits[kMCEventGeneratorId] = "";
  fgVariableUnits[kMCEventSubGeneratorId] = "";
  fgVariableUnits[kMCVtxX] = "cm";
  fgVariableUnits[kMCVtxY] = "cm";
  fgVariableUnits[kMCVtxZ] = "cm";
  fgVariableUnits[kMCEventTime] = ""; // TODO: add proper unit
  fgVariableUnits[kMCEventWeight] = "";
  fgVariableUnits[kMCEventImpParam] = "b";
  fgVariableNames[kTwoEvPosZ1] = "vtx-z_{1}";
  fgVariableUnits[kTwoEvPosZ1] = "cm";
  fgVariableNames[kTwoEvPosZ2] = "vtx-z_{2}";
  fgVariableUnits[kTwoEvPosZ2] = "cm";
  fgVariableNames[kTwoEvPosR1] = "vtx-R_{1}";
  fgVariableUnits[kTwoEvPosR1] = "cm";
  fgVariableNames[kTwoEvPosR2] = "vtx-R_{2}";
  fgVariableUnits[kTwoEvPosR2] = "cm";
  fgVariableNames[kTwoEvDeltaZ] = "#Delta_{z}";
  fgVariableUnits[kTwoEvDeltaZ] = "cm";
  fgVariableNames[kTwoEvDeltaR] = "#Delta_{R}";
  fgVariableUnits[kTwoEvDeltaR] = "cm";
  fgVariableNames[kTwoEvDeltaX] = "#Delta_{x}";
  fgVariableUnits[kTwoEvDeltaX] = "cm";
  fgVariableNames[kTwoEvDeltaY] = "#Delta_{y}";
  fgVariableUnits[kTwoEvDeltaY] = "cm";
  fgVariableNames[kTwoEvPVcontrib1] = "n.contrib 1";
  fgVariableUnits[kTwoEvPVcontrib1] = "";
  fgVariableNames[kTwoEvPVcontrib2] = "n.contrib 2";
  fgVariableUnits[kTwoEvPVcontrib2] = "";
  fgVariableNames[kEnergyCommonZNA] = "ZNA common energy";
  fgVariableUnits[kEnergyCommonZNA] = "";
  fgVariableNames[kEnergyCommonZNC] = "ZNC common energy";
  fgVariableUnits[kEnergyCommonZNC] = "";
  fgVariableNames[kEnergyCommonZPA] = "ZPA common energy";
  fgVariableUnits[kEnergyCommonZPA] = "";
  fgVariableNames[kEnergyCommonZPC] = "ZPC common energy";
  fgVariableUnits[kEnergyCommonZPC] = "";
  fgVariableNames[kTimeZNA] = "ZNA time";
  fgVariableUnits[kTimeZNA] = "";
  fgVariableNames[kTimeZNC] = "ZNC time";
  fgVariableUnits[kTimeZNC] = "";
  fgVariableNames[kTimeZPA] = "ZPA time";
  fgVariableUnits[kTimeZPA] = "";
  fgVariableNames[kTimeZPC] = "ZPC time";
  fgVariableUnits[kTimeZPC] = "";
  fgVariableNames[kMultNTracksHasITS] = "#tracks in PV with ITS";
  fgVariableUnits[kMultNTracksHasITS] = "";
  fgVariableNames[kMultNTracksHasTPC] = "#tracks in PV with TPC";
  fgVariableUnits[kMultNTracksHasTPC] = "";
  fgVariableNames[kMultNTracksHasTOF] = "#tracks in PV with TOF";
  fgVariableUnits[kMultNTracksHasTOF] = "";
  fgVariableNames[kMultNTracksHasTRD] = "#tracks in PV with TRD";
  fgVariableUnits[kMultNTracksHasTRD] = "";
  fgVariableNames[kMultNTracksITSOnly] = "# ITS only tracks in PV";
  fgVariableUnits[kMultNTracksITSOnly] = "";
  fgVariableNames[kMultNTracksTPCOnly] = "# TPC only tracks in PV";
  fgVariableUnits[kMultNTracksTPCOnly] = "";
  fgVariableNames[kMultNTracksITSTPC] = "# ITS-TPC tracks in PV";
  fgVariableUnits[kMultNTracksITSTPC] = "";
  fgVariableNames[kTrackOccupancyInTimeRange] = "track occupancy in TPC drift time (PV tracks)";
  fgVariableUnits[kTrackOccupancyInTimeRange] = "";
  fgVariableNames[kNoCollInTimeRangeStandard] = "track occupancy in TPC drift standart time";
  fgVariableUnits[kNoCollInTimeRangeStandard] = "";
  fgVariableNames[kMultAllTracksITSTPC] = "# ITS-TPC tracks";
  fgVariableUnits[kMultAllTracksITSTPC] = "";
  fgVariableNames[kMultAllTracksTPCOnly] = "# TPC only tracks";
  fgVariableUnits[kMultAllTracksTPCOnly] = "";
  fgVariableNames[kNTPCpileupContribA] = "# TPC pileup contributors on A side";
  fgVariableUnits[kNTPCpileupContribA] = "";
  fgVariableNames[kNTPCpileupContribC] = "# TPC pileup contributors on C side";
  fgVariableUnits[kNTPCpileupContribC] = "";
  fgVariableNames[kNTPCpileupZA] = "# TPC pileup mean-Z on A side";
  fgVariableUnits[kNTPCpileupZA] = "";
  fgVariableNames[kNTPCpileupZC] = "# TPC pileup mean-Z on C side";
  fgVariableUnits[kNTPCpileupZC] = "";
  fgVariableNames[kNTPCtracksInPast] = "# TPC tracks in past";
  fgVariableUnits[kNTPCtracksInPast] = "";
  fgVariableNames[kNTPCtracksInFuture] = "# TPC tracks in future";
  fgVariableUnits[kNTPCtracksInFuture] = "";
  fgVariableNames[kNTPCcontribLongA] = "# TPC-A pileup, long time range";
  fgVariableUnits[kNTPCcontribLongA] = "";
  fgVariableNames[kNTPCcontribLongC] = "# TPC-C pileup, long time range";
  fgVariableUnits[kNTPCcontribLongC] = "";
  fgVariableNames[kNTPCmeanTimeLongA] = "# TPC-A pileup mean time, long time range";
  fgVariableUnits[kNTPCmeanTimeLongA] = "#mu s";
  fgVariableNames[kNTPCmeanTimeLongC] = "# TPC-C pileup mean time, long time range";
  fgVariableUnits[kNTPCmeanTimeLongC] = "#mu s";
  fgVariableNames[kNTPCmedianTimeLongA] = "# TPC-A pileup median time, long time range";
  fgVariableUnits[kNTPCmedianTimeLongA] = "#mu s";
  fgVariableNames[kNTPCmedianTimeLongC] = "# TPC-C pileup median time, long time range";
  fgVariableUnits[kNTPCmedianTimeLongC] = "#mu s";
  fgVariableNames[kNTPCcontribShortA] = "# TPC-A pileup, short time range";
  fgVariableUnits[kNTPCcontribShortA] = "";
  fgVariableNames[kNTPCcontribShortC] = "# TPC-C pileup, short time range";
  fgVariableUnits[kNTPCcontribShortC] = "";
  fgVariableNames[kNTPCmeanTimeShortA] = "# TPC-A pileup mean time, short time range";
  fgVariableUnits[kNTPCmeanTimeShortA] = "#mu s";
  fgVariableNames[kNTPCmeanTimeShortC] = "# TPC-C pileup mean time, short time range";
  fgVariableUnits[kNTPCmeanTimeShortC] = "#mu s";
  fgVariableNames[kNTPCmedianTimeShortA] = "# TPC-A pileup median time, short time range";
  fgVariableUnits[kNTPCmedianTimeShortA] = "#mu s";
  fgVariableNames[kNTPCmedianTimeShortC] = "# TPC-C pileup median time, short time range";
  fgVariableUnits[kNTPCmedianTimeShortC] = "#mu s";
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
  fgVariableNames[kPin_leg1] = "p_{IN}";
  fgVariableUnits[kPin_leg1] = "GeV/c";
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
  fgVariableNames[kITSmeanClsSize] = "ITS mean Cls Size";
  fgVariableUnits[kITSmeanClsSize] = "";
  fgVariableNames[kTPCncls] = "TPC #cls";
  fgVariableUnits[kTPCncls] = "";
  fgVariableNames[kTPCnclsCR] = "TPC #cls crossed rows";
  fgVariableUnits[kTPCnclsCR] = "";
  fgVariableNames[kTPCnCRoverFindCls] = "TPC crossed rows over findable cls";
  fgVariableUnits[kTPCnCRoverFindCls] = "";
  fgVariableNames[kTPCchi2] = "TPC chi2";
  fgVariableUnits[kTPCchi2] = "";
  fgVariableNames[kTPCsignal] = "TPC dE/dx";
  fgVariableUnits[kTPCsignal] = "";
  fgVariableNames[kPhiTPCOuter] = "#varphi_{TPCout}";
  fgVariableUnits[kPhiTPCOuter] = "rad.";
  fgVariableNames[kTrackIsInsideTPCModule] = "Track is in TPC module";
  fgVariableUnits[kTrackIsInsideTPCModule] = "";
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
  fgVariableNames[kTPCnSigmaKa_leg1] = "n #sigma_{K}^{TPC}";
  fgVariableUnits[kTPCnSigmaKa_leg1] = "";
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
  fgVariableNames[kMftNClusters] = "MFT n-clusters";
  fgVariableUnits[kMftNClusters] = "";
  fgVariableNames[kMftClusterSize] = "MFT cluster size";
  fgVariableUnits[kMftClusterSize] = "";
  fgVariableNames[kMftMeanClusterSize] = "<MFT cluster size>";
  fgVariableUnits[kMftMeanClusterSize] = "";
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
  fgVariableNames[kVertexingLxyProjected] = "Pair Lxy";
  fgVariableUnits[kVertexingLxyProjected] = "cm";
  fgVariableNames[kVertexingLzProjected] = "Pair Lz";
  fgVariableUnits[kVertexingLzProjected] = "cm";
  fgVariableNames[kVertexingLxyzProjected] = "Pair Lxyz";
  fgVariableUnits[kVertexingLxyzProjected] = "cm";
  fgVariableNames[kVertexingTauzProjected] = "Pair pseudo-proper Tauz";
  fgVariableUnits[kVertexingTauzProjected] = "ns";
  fgVariableNames[kVertexingTauxyProjected] = "Pair pseudo-proper Tauxy";
  fgVariableUnits[kVertexingTauxyProjected] = "ns";
  fgVariableNames[kVertexingTauxyzProjected] = "Pair pseudo-proper Tauxyz";
  fgVariableUnits[kVertexingTauxyzProjected] = "ns";
  fgVariableNames[kCosPointingAngle] = "cos(#theta_{pointing})";
  fgVariableUnits[kCosPointingAngle] = "";
  fgVariableNames[kVertexingPz] = "Pz Pair";
  fgVariableUnits[kVertexingPz] = "GeV/c";
  fgVariableNames[kVertexingSV] = "Secondary Vertexing z";
  fgVariableUnits[kVertexingSV] = "cm";
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
  fgVariableNames[kCosPointingAngle] = "Cos #theta_{pointing}";
  fgVariableUnits[kCosPointingAngle] = "";
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
  fgVariableNames[kQ1ZNAX] = "Q_{1,x}^{ZNA} ";
  fgVariableUnits[kQ1ZNAX] = "";
  fgVariableNames[kQ1ZNAY] = "Q_{1,y}^{ZNA} ";
  fgVariableUnits[kQ1ZNAY] = "";
  fgVariableNames[kQ1ZNCX] = "Q_{1,x}^{ZNC} ";
  fgVariableUnits[kQ1ZNCX] = "";
  fgVariableNames[kQ1ZNCY] = "Q_{1,y}^{ZNC} ";
  fgVariableUnits[kQ1ZNCY] = "";
  fgVariableNames[KIntercalibZNA] = "ZNA^{common} - (ZNA1 + ZNA2 + ZNA3 + ZNA4)";
  fgVariableUnits[KIntercalibZNA] = "";
  fgVariableNames[KIntercalibZNC] = "ZNC^{common} - (ZNC1 + ZNC2 + ZNC3 + ZNC4)";
  fgVariableUnits[KIntercalibZNC] = "";
  fgVariableNames[kQ1ZNACXX] = "Q_{1,x}^{ZNC} #dot Q_{1,x}^{ZNA} ";
  fgVariableUnits[kQ1ZNACXX] = "";
  fgVariableNames[kQ1ZNACYY] = "Q_{1,y}^{ZNC} #dot Q_{1,y}^{ZNA} ";
  fgVariableUnits[kQ1ZNACYY] = "";
  fgVariableNames[kQ1ZNACYX] = "Q_{1,y}^{ZNC} #dot Q_{1,x}^{ZNA} ";
  fgVariableUnits[kQ1ZNACYX] = "";
  fgVariableNames[kQ1ZNACXY] = "Q_{1,x}^{ZNC} #dot Q_{1,y}^{ZNA} ";
  fgVariableUnits[kQ1ZNACXY] = "";
  fgVariableNames[kQ1X0A] = "Q_{1,x}^{A} ";
  fgVariableUnits[kQ1X0A] = "";
  fgVariableNames[kQ1Y0A] = "Q_{1,y}^{A} ";
  fgVariableUnits[kQ1Y0A] = "";
  fgVariableNames[kQ1X0B] = "Q_{1,x}^{B} ";
  fgVariableUnits[kQ1X0B] = "";
  fgVariableNames[kQ1Y0B] = "Q_{1,y}^{B} ";
  fgVariableUnits[kQ1Y0B] = "";
  fgVariableNames[kQ1X0C] = "Q_{1,x}^{C} ";
  fgVariableUnits[kQ1X0C] = "";
  fgVariableNames[kQ1Y0C] = "Q_{1,y}^{C} ";
  fgVariableUnits[kQ1Y0C] = "";
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
  fgVariableNames[kQ2YYAB] = "<Q_{2,y}^{A}*Q_{2,y}^{B}> ";
  fgVariableUnits[kQ2YYAB] = "";
  fgVariableNames[kQ2XXAB] = "<Q_{2,x}^{A}*Q_{2,x}^{B}> ";
  fgVariableUnits[kQ2XXAB] = "";
  fgVariableNames[kQ2XYAB] = "<Q_{2,x}^{A}*Q_{2,y}^{B}> ";
  fgVariableUnits[kQ2XYAB] = "";
  fgVariableNames[kQ2YXAB] = "<Q_{2,y}^{A}*Q_{2,x}^{B}> ";
  fgVariableUnits[kQ2YXAB] = "";
  fgVariableNames[kQ2YYAC] = "<Q_{2,y}^{A}*Q_{2,y}^{C}> ";
  fgVariableUnits[kQ2YYAC] = "";
  fgVariableNames[kQ2XXAC] = "<Q_{2,x}^{A}*Q_{2,x}^{C}> ";
  fgVariableUnits[kQ2XXAC] = "";
  fgVariableNames[kQ2XYAC] = "<Q_{2,x}^{A}*Q_{2,y}^{C}> ";
  fgVariableUnits[kQ2XYAC] = "";
  fgVariableNames[kQ2YXAC] = "<Q_{2,y}^{A}*Q_{2,x}^{C}> ";
  fgVariableUnits[kQ2YXAC] = "";
  fgVariableNames[kQ2YYBC] = "<Q_{2,y}^{B}*Q_{2,y}^{C}> ";
  fgVariableUnits[kQ2YYBC] = "";
  fgVariableNames[kQ2XXBC] = "<Q_{2,x}^{B}*Q_{2,x}^{C}> ";
  fgVariableUnits[kQ2XXBC] = "";
  fgVariableNames[kQ2XYBC] = "<Q_{2,x}^{B}*Q_{2,y}^{C}> ";
  fgVariableUnits[kQ2XYBC] = "";
  fgVariableNames[kQ2YXBC] = "<Q_{2,y}^{B}*Q_{2,x}^{C}> ";
  fgVariableUnits[kQ2YXBC] = "";
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
  fgVariableNames[kQ4X0A] = "Q_{4,x}^{A} ";
  fgVariableUnits[kQ4X0A] = "";
  fgVariableNames[kQ4Y0A] = "Q_{4,y}^{A} ";
  fgVariableUnits[kQ4Y0A] = "";
  fgVariableNames[kQ4X0B] = "Q_{4,x}^{B} ";
  fgVariableUnits[kQ4X0B] = "";
  fgVariableNames[kQ4Y0B] = "Q_{4,y}^{B} ";
  fgVariableUnits[kQ4Y0B] = "";
  fgVariableNames[kQ4X0C] = "Q_{4,x}^{C} ";
  fgVariableUnits[kQ4X0C] = "";
  fgVariableNames[kQ4Y0C] = "Q_{4,y}^{C} ";
  fgVariableUnits[kQ4Y0C] = "";
  fgVariableNames[kU2Q2] = "u_{2}Q_{2}^{A} ";
  fgVariableUnits[kU2Q2] = "";
  fgVariableNames[kU3Q3] = "u_{3}Q_{3}^{A} ";
  fgVariableUnits[kU3Q3] = "";
  fgVariableNames[kQ42XA] = "Q_{42,x}^{A} ";
  fgVariableUnits[kQ42XA] = "";
  fgVariableNames[kQ42YA] = "Q_{42,y}^{A} ";
  fgVariableUnits[kQ42YA] = "";
  fgVariableNames[kQ23XA] = "Q_{23,x}^{A} ";
  fgVariableUnits[kQ23XA] = "";
  fgVariableNames[kQ23YA] = "Q_{23,y}^{A} ";
  fgVariableUnits[kQ23YA] = "";
  fgVariableNames[kS11A] = "S_{11}^{A} ";
  fgVariableUnits[kS11A] = "";
  fgVariableNames[kS12A] = "S_{12}^{A} ";
  fgVariableUnits[kS12A] = "";
  fgVariableNames[kS13A] = "S_{13}^{A} ";
  fgVariableUnits[kS13A] = "";
  fgVariableNames[kS31A] = "S_{31}^{A} ";
  fgVariableUnits[kS31A] = "";
  fgVariableNames[kM11REF] = "M_{11}^{REF} ";
  fgVariableUnits[kM11REF] = "";
  fgVariableNames[kM11REFetagap] = "M_{11}^{REF}-etagap ";
  fgVariableUnits[kM11REFetagap] = "";
  fgVariableNames[kM01POI] = "M^{'}_{01}^{POI} ";
  fgVariableUnits[kM01POI] = "";
  fgVariableNames[kM1111REF] = "M_{1111}^{REF} ";
  fgVariableUnits[kM1111REF] = "";
  fgVariableNames[kM11M1111REF] = "M_{11}_{REF}M_{1111}_{REF}";
  fgVariableUnits[kM11M1111REF] = "";
  fgVariableNames[kM11M1111REFoverMp] = "M_{11}_{REF}M_{1111}_{REF} / M_{p}  ";
  fgVariableUnits[kM11M1111REFoverMp] = "";
  fgVariableNames[kM01M0111POIoverMp] = "M_{01}_{POI}M_{0111}_{POI} / M_{p}";
  fgVariableUnits[kM01M0111POIoverMp] = "";
  fgVariableNames[kCORR2CORR4REF] = "<2><4>";
  fgVariableUnits[kCORR2CORR4REF] = "";
  fgVariableNames[kCORR2POICORR4POI] = "<2'><4'>";
  fgVariableUnits[kCORR2POICORR4POI] = "";
  fgVariableNames[kCORR2REFCORR4POI] = "<2><4'>";
  fgVariableUnits[kCORR2REFCORR4POI] = "";
  fgVariableNames[kCORR2REFCORR2POI] = "<2><2'>";
  fgVariableUnits[kCORR2REFCORR2POI] = "";
  fgVariableNames[kM01M0111overMp] = "M_{01}_{POI} M_{0111}_{POI} / M_{p}  ";
  fgVariableUnits[kM01M0111overMp] = "";
  fgVariableNames[kM11M0111overMp] = "M_{11}_{REF}M_{0111}_{POI} / M_{p}  ";
  fgVariableUnits[kM11M0111overMp] = "";
  fgVariableNames[kM11M01overMp] = "M_{11}_{REF}M_{01}_{POI} / M_{p}  ";
  fgVariableUnits[kM11M01overMp] = "";
  fgVariableNames[kM0111POI] = "M^{'}_{0111}^{POI} ";
  fgVariableUnits[kM0111POI] = "";
  fgVariableNames[kCORR2REF] = "<2> ";
  fgVariableUnits[kCORR2REF] = "";
  fgVariableNames[kM01POIplus] = "M_{01}_{POI}^{+} ";
  fgVariableUnits[kM01POIplus] = "";
  fgVariableNames[kM0111POIplus] = "M^{'}_{0111}^{POI+} ";
  fgVariableUnits[kM0111POIplus] = "";
  fgVariableNames[kM01POIminus] = "M_{01}_{POI}^{-} ";
  fgVariableUnits[kM01POIminus] = "";
  fgVariableNames[kM0111POIminus] = "M^{'}_{0111}^{POI-} ";
  fgVariableUnits[kM0111POIminus] = "";
  fgVariableNames[kM01POIoverMpminus] = "M_{01}_{POI}^{-} / M_{p} ";
  fgVariableUnits[kM01POIoverMpminus] = "";
  fgVariableNames[kM01POIoverMpplus] = "M_{01}_{POI}^{+} / M_{p} ";
  fgVariableUnits[kM01POIoverMpplus] = "";
  fgVariableNames[kM01POIoverMpmoins] = "M_{01}_{POI}^{-} / M_{p} ";
  fgVariableUnits[kM01POIoverMpmoins] = "";
  fgVariableNames[kM01POIoverMpplus] = "M_{01}_{POI}^{+} / M_{p} ";
  fgVariableUnits[kM01POIoverMpplus] = "";
  fgVariableNames[kM01POIoverMpmoins] = "M_{01}_{POI}^{-} / M_{p} ";
  fgVariableUnits[kM01POIoverMpmoins] = "";
  fgVariableNames[kM0111POIoverMpminus] = "M^{'}_{0111}^{POI-} / M_{p} ";
  fgVariableUnits[kM0111POIoverMpminus] = "";
  fgVariableNames[kM0111POIoverMpplus] = "M^{'}_{0111}^{POI+} / M_{p} ";
  fgVariableUnits[kM0111POIoverMpplus] = "";
  fgVariableNames[kCORR2POIplus] = "<2>_{POI}^{+} ";
  fgVariableUnits[kCORR2POIplus] = "";
  fgVariableNames[kCORR2POIminus] = "<2>_{POI}^{-} ";
  fgVariableUnits[kCORR2POIminus] = "";
  fgVariableNames[kCORR4POIplus] = "<4>_{POI}^{+} ";
  fgVariableUnits[kCORR4POIplus] = "";
  fgVariableNames[kCORR4POIminus] = "<4>_{POI}^{-} ";
  fgVariableUnits[kCORR4POIminus] = "";
  fgVariableNames[kCORR2REFbydimuons] = "<2> only for events with dimuons";
  fgVariableUnits[kCORR2REFbydimuons] = "";
  fgVariableNames[kCORR2REFetagap] = "<2-etagap> ";
  fgVariableUnits[kCORR2REFetagap] = "";
  fgVariableNames[kCORR2POI] = "<2'> ";
  fgVariableUnits[kCORR2POI] = "";
  fgVariableNames[kCORR4REF] = "<4> ";
  fgVariableUnits[kCORR4REF] = "";
  fgVariableNames[kCORR4REFbydimuons] = "<4> only for events with dimuons";
  fgVariableUnits[kCORR4REFbydimuons] = "";
  fgVariableNames[kCORR4POI] = "<4'> ";
  fgVariableUnits[kCORR4POI] = "";
  fgVariableNames[kM11REFoverMp] = "M_{11}^{REF}/M_{p} ";
  fgVariableUnits[kM11REFoverMp] = "";
  fgVariableNames[kM01POIoverMp] = "M^{'}_{01}^{POI}/M_{p} ";
  fgVariableUnits[kM01POIoverMp] = "";
  fgVariableNames[kM1111REFoverMp] = "M_{1111}^{REF}/M_{p} ";
  fgVariableUnits[kM1111REFoverMp] = "";
  fgVariableNames[kM0111POIoverMp] = "M^{'}_{0111}^{POI}/M_{p} ";
  fgVariableUnits[kM0111POIoverMp] = "";
  fgVariableNames[kCORR2POIMp] = "<2'> M_{p} ";
  fgVariableUnits[kCORR2POIMp] = "";
  fgVariableNames[kCORR4POIMp] = "<4'> M_{p} ";
  fgVariableUnits[kCORR4POIMp] = "";
  fgVariableNames[kCos2DeltaPhi] = "cos 2(#varphi-#Psi_{2}^{A}) ";
  fgVariableUnits[kCos2DeltaPhi] = "";
  fgVariableNames[kCos3DeltaPhi] = "cos 3(#varphi-#Psi_{3}^{A}) ";
  fgVariableUnits[kCos3DeltaPhi] = "";
  fgVariableNames[kPsi2A] = "#Psi_{2}^{A} ";
  fgVariableUnits[kPsi2A] = "";
  fgVariableNames[kPsi2B] = "#Psi_{2}^{B} ";
  fgVariableUnits[kPsi2B] = "";
  fgVariableNames[kPsi2C] = "#Psi_{2}^{C} ";
  fgVariableUnits[kPsi2C] = "";
  fgVariableNames[kR2SP_AB] = "R_{2}^{SP} (AB) ";
  fgVariableUnits[kR2SP_AB] = "";
  fgVariableNames[kR2SP_AC] = "R_{2}^{SP} (AC) ";
  fgVariableUnits[kR2SP_AC] = "";
  fgVariableNames[kR2SP_BC] = "R_{2}^{SP} (BC) ";
  fgVariableUnits[kR2SP_BC] = "";
  fgVariableNames[kR2SP_FT0CTPCPOS] = "R_{2}^{SP} (FT0C-TPCpos) ";
  fgVariableUnits[kR2SP_FT0CTPCPOS] = "";
  fgVariableNames[kR2SP_FT0CTPCNEG] = "R_{2}^{SP}  (FT0C-TPCneg) ";
  fgVariableUnits[kR2SP_FT0CTPCNEG] = "";
  fgVariableNames[kR2SP_FT0ATPCPOS] = "R_{2}^{SP} (FT0A-TPCpos) ";
  fgVariableUnits[kR2SP_FT0ATPCPOS] = "";
  fgVariableNames[kR2SP_FT0ATPCNEG] = "R_{2}^{SP} (FT0A-TPCneg) ";
  fgVariableUnits[kR2SP_FT0ATPCNEG] = "";
  fgVariableNames[kR2EP_AB] = "R_{2}^{EP} (TPC-FT0A) ";
  fgVariableUnits[kR2EP_AB] = "";
  fgVariableNames[kR2EP_AC] = "R_{2}^{EP} (TPC-FT0C) ";
  fgVariableUnits[kR2EP_AC] = "";
  fgVariableNames[kR2EP_BC] = "R_{2}^{EP} (FT0C-FT0A) ";
  fgVariableUnits[kR2EP_BC] = "";
  fgVariableNames[kR2EP_FT0CTPCPOS] = "R_{2}^{EP} (FT0C-TPCpos) ";
  fgVariableUnits[kR2EP_FT0CTPCPOS] = "";
  fgVariableNames[kR2EP_FT0CTPCNEG] = "R_{2}^{EP}  (FT0C-TPCneg) ";
  fgVariableUnits[kR2EP_FT0CTPCNEG] = "";
  fgVariableNames[kR2EP_FT0ATPCPOS] = "R_{2}^{EP} (FT0A-TPCpos) ";
  fgVariableUnits[kR2EP_FT0ATPCPOS] = "";
  fgVariableNames[kR2EP_FT0ATPCNEG] = "R_{2}^{EP} (FT0A-TPCneg) ";
  fgVariableUnits[kR2EP_FT0ATPCNEG] = "";
  fgVariableNames[kR3SP] = "R_{3}^{SP} ";
  fgVariableUnits[kR3SP] = "";
  fgVariableNames[kR3EP] = "R_{3}^{EP} ";
  fgVariableUnits[kR3EP] = "";
  fgVariableNames[kPairMass] = "mass";
  fgVariableUnits[kPairMass] = "GeV/c2";
  fgVariableNames[kPairMassDau] = "mass dilepton";
  fgVariableUnits[kPairMassDau] = "GeV/c2";
  fgVariableNames[kDeltaMass] = "mass - dilepton mass";
  fgVariableUnits[kDeltaMass] = "GeV/c2";
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
  fgVariableNames[kCosPhiVP] = "cos#it{#varphi}_{VP}";
  fgVariableUnits[kCosPhiVP] = "";
  fgVariableNames[kPhiVP] = "#varphi_{VP} - #Psi_{2}";
  fgVariableUnits[kPhiVP] = "rad.";
  fgVariableNames[kDeltaPhiPair2] = "#Delta#phi";
  fgVariableUnits[kDeltaPhiPair2] = "rad.";
  fgVariableNames[kDeltaEtaPair2] = "#Delta#eta";
  fgVariableUnits[kDeltaEtaPair2] = "";
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
  fgVariableNames[kSignQuadDCAsigXY] = "signDCA_{xy}^{quad}";
  fgVariableUnits[kSignQuadDCAsigXY] = "#sigma";
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
  fgVariableNames[kBdtCharmHadron] = "BDT score (charm hadron)";
  fgVariableUnits[kBdtCharmHadron] = " ";
  fgVariableNames[kIsDoubleGap] = "is double gap event";
  fgVariableUnits[kIsDoubleGap] = "";
  fgVariableNames[kIsSingleGapA] = "is single gap event side A";
  fgVariableUnits[kIsSingleGapA] = "";
  fgVariableNames[kIsSingleGapC] = "is single gap event side C";
  fgVariableUnits[kIsSingleGapC] = "";
  fgVariableNames[kIsSingleGap] = "is single gap event";
  fgVariableUnits[kIsSingleGap] = "";
  fgVariableNames[kIsITSUPCMode] = "UPC settings used";
  fgVariableUnits[kIsITSUPCMode] = "";
  fgVariableNames[kQuadMass] = "mass quadruplet";
  fgVariableUnits[kQuadMass] = "GeV/c2";
  fgVariableNames[kQuadPt] = "p_{T}";
  fgVariableUnits[kQuadPt] = "GeV/c";
  fgVariableNames[kQuadEta] = "#eta";
  fgVariableUnits[kQuadEta] = "";
  fgVariableNames[kQuadPhi] = "#varphi";
  fgVariableUnits[kQuadPhi] = "rad.";
  fgVariableNames[kCosthetaDileptonDitrack] = "cos#it{#theta}_{dilepton-ditrack}";
  fgVariableUnits[kCosthetaDileptonDitrack] = "";
  fgVariableNames[kDitrackMass] = "mass di-track";
  fgVariableUnits[kDitrackMass] = "GeV/c2";
  fgVariableNames[kDitrackPt] = "p_{T}";
  fgVariableUnits[kDitrackPt] = "GeV/c";
  fgVariableNames[kQ] = "mass difference";
  fgVariableUnits[kQ] = "GeV/c2";
  fgVariableNames[kDeltaR1] = "angular distance prong 1";
  fgVariableUnits[kDeltaR1] = "";
  fgVariableNames[kDeltaR2] = "angular distance prong 2";
  fgVariableUnits[kDeltaR2] = "";

  // Set the variables short names map. This is needed for dynamic configuration via JSON files
  fgVarNamesMap["kNothing"] = kNothing;
  fgVarNamesMap["kRunNo"] = kRunNo;
  fgVarNamesMap["kRunId"] = kRunId;
  fgVarNamesMap["kRunIndex"] = kRunIndex;
  fgVarNamesMap["kNRunWiseVariables"] = kNRunWiseVariables;
  fgVarNamesMap["kTimestamp"] = kTimestamp;
  fgVarNamesMap["kTimeFromSOR"] = kTimeFromSOR;
  fgVarNamesMap["kCollisionTime"] = kCollisionTime;
  fgVarNamesMap["kCollisionTimeRes"] = kCollisionTimeRes;
  fgVarNamesMap["kBC"] = kBC;
  fgVarNamesMap["kBCOrbit"] = kBCOrbit;
  fgVarNamesMap["kIsPhysicsSelection"] = kIsPhysicsSelection;
  fgVarNamesMap["kIsNoTFBorder"] = kIsNoTFBorder;
  fgVarNamesMap["kIsNoITSROFBorder"] = kIsNoITSROFBorder;
  fgVarNamesMap["kIsNoITSROFBorderRecomputed"] = kIsNoITSROFBorderRecomputed;
  fgVarNamesMap["kIsNoSameBunch"] = kIsNoSameBunch;
  fgVarNamesMap["kIsGoodZvtxFT0vsPV"] = kIsGoodZvtxFT0vsPV;
  fgVarNamesMap["kIsVertexITSTPC"] = kIsVertexITSTPC;
  fgVarNamesMap["kIsVertexTOFmatched"] = kIsVertexTOFmatched;
  fgVarNamesMap["kIsSel8"] = kIsSel8;
  fgVarNamesMap["kIsGoodITSLayer3"] = kIsGoodITSLayer3;
  fgVarNamesMap["kIsGoodITSLayer0123"] = kIsGoodITSLayer0123;
  fgVarNamesMap["kIsGoodITSLayersAll"] = kIsGoodITSLayersAll;
  fgVarNamesMap["kIsINT7"] = kIsINT7;
  fgVarNamesMap["kIsEMC7"] = kIsEMC7;
  fgVarNamesMap["kIsINT7inMUON"] = kIsINT7inMUON;
  fgVarNamesMap["kIsMuonSingleLowPt7"] = kIsMuonSingleLowPt7;
  fgVarNamesMap["kIsMuonSingleHighPt7"] = kIsMuonSingleHighPt7;
  fgVarNamesMap["kIsMuonUnlikeLowPt7"] = kIsMuonUnlikeLowPt7;
  fgVarNamesMap["kIsMuonLikeLowPt7"] = kIsMuonLikeLowPt7;
  fgVarNamesMap["kIsCUP8"] = kIsCUP8;
  fgVarNamesMap["kIsCUP9"] = kIsCUP9;
  fgVarNamesMap["kIsMUP10"] = kIsMUP10;
  fgVarNamesMap["kIsMUP11"] = kIsMUP11;
  fgVarNamesMap["kVtxX"] = kVtxX;
  fgVarNamesMap["kVtxY"] = kVtxY;
  fgVarNamesMap["kVtxZ"] = kVtxZ;
  fgVarNamesMap["kVtxNcontrib"] = kVtxNcontrib;
  fgVarNamesMap["kVtxNcontribReal"] = kVtxNcontribReal;
  fgVarNamesMap["kVtxCovXX"] = kVtxCovXX;
  fgVarNamesMap["kVtxCovXY"] = kVtxCovXY;
  fgVarNamesMap["kVtxCovXZ"] = kVtxCovXZ;
  fgVarNamesMap["kVtxCovYY"] = kVtxCovYY;
  fgVarNamesMap["kVtxCovYZ"] = kVtxCovYZ;
  fgVarNamesMap["kVtxCovZZ"] = kVtxCovZZ;
  fgVarNamesMap["kVtxChi2"] = kVtxChi2;
  fgVarNamesMap["kCentVZERO"] = kCentVZERO;
  fgVarNamesMap["kCentFT0C"] = kCentFT0C;
  fgVarNamesMap["kMultTPC"] = kMultTPC;
  fgVarNamesMap["kMultFV0A"] = kMultFV0A;
  fgVarNamesMap["kMultFV0C"] = kMultFV0C;
  fgVarNamesMap["kMultFT0A"] = kMultFT0A;
  fgVarNamesMap["kMultFT0C"] = kMultFT0C;
  fgVarNamesMap["kMultFDDA"] = kMultFDDA;
  fgVarNamesMap["kMultFDDC"] = kMultFDDC;
  fgVarNamesMap["kMultZNA"] = kMultZNA;
  fgVarNamesMap["kMultZNC"] = kMultZNC;
  fgVarNamesMap["kMultTracklets"] = kMultTracklets;
  fgVarNamesMap["kMultDimuons"] = kMultDimuons;
  fgVarNamesMap["kMultNTracksHasITS"] = kMultNTracksHasITS;
  fgVarNamesMap["kMultNTracksHasTPC"] = kMultNTracksHasTPC;
  fgVarNamesMap["kMultNTracksHasTOF"] = kMultNTracksHasTOF;
  fgVarNamesMap["kMultNTracksHasTRD"] = kMultNTracksHasTRD;
  fgVarNamesMap["kMultNTracksITSOnly"] = kMultNTracksITSOnly;
  fgVarNamesMap["kMultNTracksTPCOnly"] = kMultNTracksTPCOnly;
  fgVarNamesMap["kMultNTracksITSTPC"] = kMultNTracksITSTPC;
  fgVarNamesMap["kTrackOccupancyInTimeRange"] = kTrackOccupancyInTimeRange;
  fgVarNamesMap["kNoCollInTimeRangeStandard"] = kNoCollInTimeRangeStandard;
  fgVarNamesMap["kMultAllTracksTPCOnly"] = kMultAllTracksTPCOnly;
  fgVarNamesMap["kMultAllTracksITSTPC"] = kMultAllTracksITSTPC;
  fgVarNamesMap["kNTPCpileupContribA"] = kNTPCpileupContribA;
  fgVarNamesMap["kNTPCpileupContribC"] = kNTPCpileupContribC;
  fgVarNamesMap["kNTPCpileupZA"] = kNTPCpileupZA;
  fgVarNamesMap["kNTPCpileupZC"] = kNTPCpileupZC;
  fgVarNamesMap["kNTPCtracksInPast"] = kNTPCtracksInPast;
  fgVarNamesMap["kNTPCtracksInFuture"] = kNTPCtracksInFuture;
  fgVarNamesMap["kNTPCcontribLongA"] = kNTPCcontribLongA;
  fgVarNamesMap["kNTPCcontribLongC"] = kNTPCcontribLongC;
  fgVarNamesMap["kNTPCmeanTimeLongA"] = kNTPCmeanTimeLongA;
  fgVarNamesMap["kNTPCmeanTimeLongC"] = kNTPCmeanTimeLongC;
  fgVarNamesMap["kNTPCmedianTimeLongA"] = kNTPCmedianTimeLongA;
  fgVarNamesMap["kNTPCmedianTimeLongC"] = kNTPCmedianTimeLongC;
  fgVarNamesMap["kNTPCcontribShortA"] = kNTPCcontribShortA;
  fgVarNamesMap["kNTPCcontribShortC"] = kNTPCcontribShortC;
  fgVarNamesMap["kNTPCmeanTimeShortA"] = kNTPCmeanTimeShortA;
  fgVarNamesMap["kNTPCmeanTimeShortC"] = kNTPCmeanTimeShortC;
  fgVarNamesMap["kNTPCmedianTimeShortA"] = kNTPCmedianTimeShortA;
  fgVarNamesMap["kNTPCmedianTimeShortC"] = kNTPCmedianTimeShortC;
  fgVarNamesMap["kMCEventGeneratorId"] = kMCEventGeneratorId;
  fgVarNamesMap["kMCEventSubGeneratorId"] = kMCEventSubGeneratorId;
  fgVarNamesMap["kMCVtxX"] = kMCVtxX;
  fgVarNamesMap["kMCVtxY"] = kMCVtxY;
  fgVarNamesMap["kMCVtxZ"] = kMCVtxZ;
  fgVarNamesMap["kMCEventTime"] = kMCEventTime;
  fgVarNamesMap["kMCEventWeight"] = kMCEventWeight;
  fgVarNamesMap["kMCEventImpParam"] = kMCEventImpParam;
  fgVarNamesMap["kQ1ZNAX"] = kQ1ZNAX;
  fgVarNamesMap["kQ1ZNAY"] = kQ1ZNAY;
  fgVarNamesMap["kQ1ZNCX"] = kQ1ZNCX;
  fgVarNamesMap["kQ1ZNCY"] = kQ1ZNCY;
  fgVarNamesMap["KIntercalibZNA"] = KIntercalibZNA;
  fgVarNamesMap["KIntercalibZNC"] = KIntercalibZNC;
  fgVarNamesMap["kQ1ZNACXX"] = kQ1ZNACXX;
  fgVarNamesMap["kQ1ZNACYY"] = kQ1ZNACYY;
  fgVarNamesMap["kQ1ZNACYX"] = kQ1ZNACYX;
  fgVarNamesMap["kQ1ZNACXY"] = kQ1ZNACXY;
  fgVarNamesMap["kQ1X0A"] = kQ1X0A;
  fgVarNamesMap["kQ1Y0A"] = kQ1Y0A;
  fgVarNamesMap["kQ1X0B"] = kQ1X0B;
  fgVarNamesMap["kQ1Y0B"] = kQ1Y0B;
  fgVarNamesMap["kQ1X0C"] = kQ1X0C;
  fgVarNamesMap["kQ1Y0C"] = kQ1Y0C;
  fgVarNamesMap["kQ2X0A"] = kQ2X0A;
  fgVarNamesMap["kQ2Y0A"] = kQ2Y0A;
  fgVarNamesMap["kQ2X0APOS"] = kQ2X0APOS;
  fgVarNamesMap["kQ2Y0APOS"] = kQ2Y0APOS;
  fgVarNamesMap["kQ2X0ANEG"] = kQ2X0ANEG;
  fgVarNamesMap["kQ2Y0ANEG"] = kQ2Y0ANEG;
  fgVarNamesMap["kQ2X0B"] = kQ2X0B;
  fgVarNamesMap["kQ2Y0B"] = kQ2Y0B;
  fgVarNamesMap["kQ2X0C"] = kQ2X0C;
  fgVarNamesMap["kQ2Y0C"] = kQ2Y0C;
  fgVarNamesMap["kQ2YYAB"] = kQ2YYAB;
  fgVarNamesMap["kQ2XXAB"] = kQ2XXAB;
  fgVarNamesMap["kQ2XYAB"] = kQ2XYAB;
  fgVarNamesMap["kQ2YXAB"] = kQ2YXAB;
  fgVarNamesMap["kQ2YYAC"] = kQ2YYAC;
  fgVarNamesMap["kQ2XXAC"] = kQ2XXAC;
  fgVarNamesMap["kQ2XYAC"] = kQ2XYAC;
  fgVarNamesMap["kQ2YXAC"] = kQ2YXAC;
  fgVarNamesMap["kQ2YYBC"] = kQ2YYBC;
  fgVarNamesMap["kQ2XXBC"] = kQ2XXBC;
  fgVarNamesMap["kQ2XYBC"] = kQ2XYBC;
  fgVarNamesMap["kQ2YXBC"] = kQ2YXBC;
  fgVarNamesMap["kMultA"] = kMultA;
  fgVarNamesMap["kMultAPOS"] = kMultAPOS;
  fgVarNamesMap["kMultANEG"] = kMultANEG;
  fgVarNamesMap["kMultB"] = kMultB;
  fgVarNamesMap["kMultC"] = kMultC;
  fgVarNamesMap["kQ3X0A"] = kQ3X0A;
  fgVarNamesMap["kQ3Y0A"] = kQ3Y0A;
  fgVarNamesMap["kQ3X0B"] = kQ3X0B;
  fgVarNamesMap["kQ3Y0B"] = kQ3Y0B;
  fgVarNamesMap["kQ3X0C"] = kQ3X0C;
  fgVarNamesMap["kQ3Y0C"] = kQ3Y0C;
  fgVarNamesMap["kQ4X0A"] = kQ4X0A;
  fgVarNamesMap["kQ4Y0A"] = kQ4Y0A;
  fgVarNamesMap["kQ4X0B"] = kQ4X0B;
  fgVarNamesMap["kQ4Y0B"] = kQ4Y0B;
  fgVarNamesMap["kQ4X0C"] = kQ4X0C;
  fgVarNamesMap["kQ4Y0C"] = kQ4Y0C;
  fgVarNamesMap["kR2SP_AB"] = kR2SP_AB;
  fgVarNamesMap["kR2SP_AC"] = kR2SP_AC;
  fgVarNamesMap["kR2SP_BC"] = kR2SP_BC;
  fgVarNamesMap["kWR2SP_AB"] = kWR2SP_AB;
  fgVarNamesMap["kWR2SP_AC"] = kWR2SP_AC;
  fgVarNamesMap["kWR2SP_BC"] = kWR2SP_BC;
  fgVarNamesMap["kR2SP_AB_Im"] = kR2SP_AB_Im;
  fgVarNamesMap["kR2SP_AC_Im"] = kR2SP_AC_Im;
  fgVarNamesMap["kR2SP_BC_Im"] = kR2SP_BC_Im;
  fgVarNamesMap["kWR2SP_AB_Im"] = kWR2SP_AB_Im;
  fgVarNamesMap["kWR2SP_AC_Im"] = kWR2SP_AC_Im;
  fgVarNamesMap["kWR2SP_BC_Im"] = kWR2SP_BC_Im;
  fgVarNamesMap["kR2SP_FT0CTPCPOS"] = kR2SP_FT0CTPCPOS;
  fgVarNamesMap["kR2SP_FT0CTPCNEG"] = kR2SP_FT0CTPCNEG;
  fgVarNamesMap["kR2SP_FT0ATPCPOS"] = kR2SP_FT0ATPCPOS;
  fgVarNamesMap["kR2SP_FT0ATPCNEG"] = kR2SP_FT0ATPCNEG;
  fgVarNamesMap["kR2SP_FT0MTPCPOS"] = kR2SP_FT0MTPCPOS;
  fgVarNamesMap["kR2SP_FT0MTPCNEG"] = kR2SP_FT0MTPCNEG;
  fgVarNamesMap["kR2SP_FV0ATPCPOS"] = kR2SP_FV0ATPCPOS;
  fgVarNamesMap["kR2SP_FV0ATPCNEG"] = kR2SP_FV0ATPCNEG;
  fgVarNamesMap["kR3SP"] = kR3SP;
  fgVarNamesMap["kR2EP_AB"] = kR2EP_AB;
  fgVarNamesMap["kR2EP_AC"] = kR2EP_AC;
  fgVarNamesMap["kR2EP_BC"] = kR2EP_BC;
  fgVarNamesMap["kWR2EP_AB"] = kWR2EP_AB;
  fgVarNamesMap["kWR2EP_AC"] = kWR2EP_AC;
  fgVarNamesMap["kWR2EP_BC"] = kWR2EP_BC;
  fgVarNamesMap["kR2EP_AB_Im"] = kR2EP_AB_Im;
  fgVarNamesMap["kR2EP_AC_Im"] = kR2EP_AC_Im;
  fgVarNamesMap["kR2EP_BC_Im"] = kR2EP_BC_Im;
  fgVarNamesMap["kWR2EP_AB_Im"] = kWR2EP_AB_Im;
  fgVarNamesMap["kWR2EP_AC_Im"] = kWR2EP_AC_Im;
  fgVarNamesMap["kWR2EP_BC_Im"] = kWR2EP_BC_Im;
  fgVarNamesMap["kR2EP_FT0CTPCPOS"] = kR2EP_FT0CTPCPOS;
  fgVarNamesMap["kR2EP_FT0CTPCNEG"] = kR2EP_FT0CTPCNEG;
  fgVarNamesMap["kR2EP_FT0ATPCPOS"] = kR2EP_FT0ATPCPOS;
  fgVarNamesMap["kR2EP_FT0ATPCNEG"] = kR2EP_FT0ATPCNEG;
  fgVarNamesMap["kR2EP_FT0MTPCPOS"] = kR2EP_FT0MTPCPOS;
  fgVarNamesMap["kR2EP_FT0MTPCNEG"] = kR2EP_FT0MTPCNEG;
  fgVarNamesMap["kR2EP_FV0ATPCPOS"] = kR2EP_FV0ATPCPOS;
  fgVarNamesMap["kR2EP_FV0ATPCNEG"] = kR2EP_FV0ATPCNEG;
  fgVarNamesMap["kR3EP"] = kR3EP;
  fgVarNamesMap["kIsDoubleGap"] = kIsDoubleGap;
  fgVarNamesMap["kIsSingleGapA"] = kIsSingleGapA;
  fgVarNamesMap["kIsSingleGapC"] = kIsSingleGapC;
  fgVarNamesMap["kIsSingleGap"] = kIsSingleGap;
  fgVarNamesMap["kIsITSUPCMode"] = kIsITSUPCMode;
  fgVarNamesMap["kTwoEvPosZ1"] = kTwoEvPosZ1;
  fgVarNamesMap["kTwoEvPosZ2"] = kTwoEvPosZ2;
  fgVarNamesMap["kTwoEvPosR1"] = kTwoEvPosR1;
  fgVarNamesMap["kTwoEvPosR2"] = kTwoEvPosR2;
  fgVarNamesMap["kTwoEvCentFT0C1"] = kTwoEvCentFT0C1;
  fgVarNamesMap["kTwoEvCentFT0C2"] = kTwoEvCentFT0C2;
  fgVarNamesMap["kTwoEvPVcontrib1"] = kTwoEvPVcontrib1;
  fgVarNamesMap["kTwoEvPVcontrib2"] = kTwoEvPVcontrib2;
  fgVarNamesMap["kTwoEvDeltaZ"] = kTwoEvDeltaZ;
  fgVarNamesMap["kTwoEvDeltaX"] = kTwoEvDeltaX;
  fgVarNamesMap["kTwoEvDeltaY"] = kTwoEvDeltaY;
  fgVarNamesMap["kTwoEvDeltaR"] = kTwoEvDeltaR;
  fgVarNamesMap["kEnergyCommonZNA"] = kEnergyCommonZNA;
  fgVarNamesMap["kEnergyCommonZNC"] = kEnergyCommonZNC;
  fgVarNamesMap["kEnergyCommonZPA"] = kEnergyCommonZPA;
  fgVarNamesMap["kEnergyCommonZPC"] = kEnergyCommonZPC;
  fgVarNamesMap["kEnergyZNA1"] = kEnergyZNA1;
  fgVarNamesMap["kEnergyZNA2"] = kEnergyZNA2;
  fgVarNamesMap["kEnergyZNA3"] = kEnergyZNA3;
  fgVarNamesMap["kEnergyZNA4"] = kEnergyZNA4;
  fgVarNamesMap["kEnergyZNC1"] = kEnergyZNC1;
  fgVarNamesMap["kEnergyZNC2"] = kEnergyZNC2;
  fgVarNamesMap["kEnergyZNC3"] = kEnergyZNC3;
  fgVarNamesMap["kEnergyZNC4"] = kEnergyZNC4;
  fgVarNamesMap["kTimeZNA"] = kTimeZNA;
  fgVarNamesMap["kTimeZNC"] = kTimeZNC;
  fgVarNamesMap["kTimeZPA"] = kTimeZPA;
  fgVarNamesMap["kTimeZPC"] = kTimeZPC;
  fgVarNamesMap["kQ2X0A1"] = kQ2X0A1;
  fgVarNamesMap["kQ2X0A2"] = kQ2X0A2;
  fgVarNamesMap["kQ2Y0A1"] = kQ2Y0A1;
  fgVarNamesMap["kQ2Y0A2"] = kQ2Y0A2;
  fgVarNamesMap["kU2Q2Ev1"] = kU2Q2Ev1;
  fgVarNamesMap["kU2Q2Ev2"] = kU2Q2Ev2;
  fgVarNamesMap["kCos2DeltaPhiEv1"] = kCos2DeltaPhiEv1;
  fgVarNamesMap["kCos2DeltaPhiEv2"] = kCos2DeltaPhiEv2;
  fgVarNamesMap["kV2SP1"] = kV2SP1;
  fgVarNamesMap["kV2SP2"] = kV2SP2;
  fgVarNamesMap["kV2EP1"] = kV2EP1;
  fgVarNamesMap["kV2EP2"] = kV2EP2;
  fgVarNamesMap["kV2ME_SP"] = kV2ME_SP;
  fgVarNamesMap["kV2ME_EP"] = kV2ME_EP;
  fgVarNamesMap["kWV2ME_SP"] = kWV2ME_SP;
  fgVarNamesMap["kWV2ME_EP"] = kWV2ME_EP;
  fgVarNamesMap["kTwoR2SP1"] = kTwoR2SP1;
  fgVarNamesMap["kTwoR2SP2"] = kTwoR2SP2;
  fgVarNamesMap["kTwoR2EP1"] = kTwoR2EP1;
  fgVarNamesMap["kTwoR2EP2"] = kTwoR2EP2;
  fgVarNamesMap["kNEventWiseVariables"] = kNEventWiseVariables;
  fgVarNamesMap["kX"] = kX;
  fgVarNamesMap["kY"] = kY;
  fgVarNamesMap["kZ"] = kZ;
  fgVarNamesMap["kPt"] = kPt;
  fgVarNamesMap["kSignedPt"] = kSignedPt;
  fgVarNamesMap["kInvPt"] = kInvPt;
  fgVarNamesMap["kEta"] = kEta;
  fgVarNamesMap["kTgl"] = kTgl;
  fgVarNamesMap["kPhi"] = kPhi;
  fgVarNamesMap["kP"] = kP;
  fgVarNamesMap["kPx"] = kPx;
  fgVarNamesMap["kPy"] = kPy;
  fgVarNamesMap["kPz"] = kPz;
  fgVarNamesMap["kRap"] = kRap;
  fgVarNamesMap["kMass"] = kMass;
  fgVarNamesMap["kCharge"] = kCharge;
  fgVarNamesMap["kNBasicTrackVariables"] = kNBasicTrackVariables;
  fgVarNamesMap["kUsedKF"] = kUsedKF;
  fgVarNamesMap["kKFMass"] = kKFMass;
  fgVarNamesMap["kKFMassGeoTop"] = kKFMassGeoTop;
  fgVarNamesMap["kPt1"] = kPt1;
  fgVarNamesMap["kEta1"] = kEta1;
  fgVarNamesMap["kPhi1"] = kPhi1;
  fgVarNamesMap["kCharge1"] = kCharge1;
  fgVarNamesMap["kPin_leg1"] = kPin_leg1;
  fgVarNamesMap["kTPCnSigmaKa_leg1"] = kTPCnSigmaKa_leg1;
  fgVarNamesMap["kPt2"] = kPt2;
  fgVarNamesMap["kEta2"] = kEta2;
  fgVarNamesMap["kPhi2"] = kPhi2;
  fgVarNamesMap["kCharge2"] = kCharge2;
  fgVarNamesMap["kPin"] = kPin;
  fgVarNamesMap["kSignedPin"] = kSignedPin;
  fgVarNamesMap["kTOFExpMom"] = kTOFExpMom;
  fgVarNamesMap["kTrackTime"] = kTrackTime;
  fgVarNamesMap["kTrackTimeRes"] = kTrackTimeRes;
  fgVarNamesMap["kTrackTimeResRelative"] = kTrackTimeResRelative;
  fgVarNamesMap["kDetectorMap"] = kDetectorMap;
  fgVarNamesMap["kHasITS"] = kHasITS;
  fgVarNamesMap["kHasTRD"] = kHasTRD;
  fgVarNamesMap["kHasTOF"] = kHasTOF;
  fgVarNamesMap["kHasTPC"] = kHasTPC;
  fgVarNamesMap["kIsGlobalTrack"] = kIsGlobalTrack;
  fgVarNamesMap["kIsGlobalTrackSDD"] = kIsGlobalTrackSDD;
  fgVarNamesMap["kIsITSrefit"] = kIsITSrefit;
  fgVarNamesMap["kIsSPDany"] = kIsSPDany;
  fgVarNamesMap["kIsSPDfirst"] = kIsSPDfirst;
  fgVarNamesMap["kIsSPDboth"] = kIsSPDboth;
  fgVarNamesMap["kIsITSibAny"] = kIsITSibAny;
  fgVarNamesMap["kIsITSibFirst"] = kIsITSibFirst;
  fgVarNamesMap["kIsITSibAll"] = kIsITSibAll;
  fgVarNamesMap["kITSncls"] = kITSncls;
  fgVarNamesMap["kITSchi2"] = kITSchi2;
  fgVarNamesMap["kITSlayerHit"] = kITSlayerHit;
  fgVarNamesMap["kITSmeanClsSize"] = kITSmeanClsSize;
  fgVarNamesMap["kIsTPCrefit"] = kIsTPCrefit;
  fgVarNamesMap["kTPCncls"] = kTPCncls;
  fgVarNamesMap["kITSClusterMap"] = kITSClusterMap;
  fgVarNamesMap["kTPCnclsCR"] = kTPCnclsCR;
  fgVarNamesMap["kTPCnCRoverFindCls"] = kTPCnCRoverFindCls;
  fgVarNamesMap["kTPCchi2"] = kTPCchi2;
  fgVarNamesMap["kTPCsignal"] = kTPCsignal;
  fgVarNamesMap["kTPCsignalRandomized"] = kTPCsignalRandomized;
  fgVarNamesMap["kTPCsignalRandomizedDelta"] = kTPCsignalRandomizedDelta;
  fgVarNamesMap["kPhiTPCOuter"] = kPhiTPCOuter;
  fgVarNamesMap["kTrackIsInsideTPCModule"] = kTrackIsInsideTPCModule;
  fgVarNamesMap["kTRDsignal"] = kTRDsignal;
  fgVarNamesMap["kTRDPattern"] = kTRDPattern;
  fgVarNamesMap["kTOFbeta"] = kTOFbeta;
  fgVarNamesMap["kTrackLength"] = kTrackLength;
  fgVarNamesMap["kTrackDCAxy"] = kTrackDCAxy;
  fgVarNamesMap["kTrackDCAxyProng1"] = kTrackDCAxyProng1;
  fgVarNamesMap["kTrackDCAxyProng2"] = kTrackDCAxyProng2;
  fgVarNamesMap["kTrackDCAz"] = kTrackDCAz;
  fgVarNamesMap["kTrackDCAzProng1"] = kTrackDCAzProng1;
  fgVarNamesMap["kTrackDCAzProng2"] = kTrackDCAzProng2;
  fgVarNamesMap["kTrackDCAsigXY"] = kTrackDCAsigXY;
  fgVarNamesMap["kTrackDCAsigZ"] = kTrackDCAsigZ;
  fgVarNamesMap["kTrackDCAresXY"] = kTrackDCAresXY;
  fgVarNamesMap["kTrackDCAresZ"] = kTrackDCAresZ;
  fgVarNamesMap["kIsGoldenChi2"] = kIsGoldenChi2;
  fgVarNamesMap["kTrackCYY"] = kTrackCYY;
  fgVarNamesMap["kTrackCZZ"] = kTrackCZZ;
  fgVarNamesMap["kTrackCSnpSnp"] = kTrackCSnpSnp;
  fgVarNamesMap["kTrackCTglTgl"] = kTrackCTglTgl;
  fgVarNamesMap["kTrackC1Pt21Pt2"] = kTrackC1Pt21Pt2;
  fgVarNamesMap["kTPCnSigmaEl"] = kTPCnSigmaEl;
  fgVarNamesMap["kTPCnSigmaElRandomized"] = kTPCnSigmaElRandomized;
  fgVarNamesMap["kTPCnSigmaElRandomizedDelta"] = kTPCnSigmaElRandomizedDelta;
  fgVarNamesMap["kTPCnSigmaMu"] = kTPCnSigmaMu;
  fgVarNamesMap["kTPCnSigmaPi"] = kTPCnSigmaPi;
  fgVarNamesMap["kTPCnSigmaPiRandomized"] = kTPCnSigmaPiRandomized;
  fgVarNamesMap["kTPCnSigmaPiRandomizedDelta"] = kTPCnSigmaPiRandomizedDelta;
  fgVarNamesMap["kTPCnSigmaKa"] = kTPCnSigmaKa;
  fgVarNamesMap["kTPCnSigmaPr"] = kTPCnSigmaPr;
  fgVarNamesMap["kTPCnSigmaEl_Corr"] = kTPCnSigmaEl_Corr;
  fgVarNamesMap["kTPCnSigmaPi_Corr"] = kTPCnSigmaPi_Corr;
  fgVarNamesMap["kTPCnSigmaKa_Corr"] = kTPCnSigmaKa_Corr;
  fgVarNamesMap["kTPCnSigmaPr_Corr"] = kTPCnSigmaPr_Corr;
  fgVarNamesMap["kTPCnSigmaPrRandomized"] = kTPCnSigmaPrRandomized;
  fgVarNamesMap["kTPCnSigmaPrRandomizedDelta"] = kTPCnSigmaPrRandomizedDelta;
  fgVarNamesMap["kTOFnSigmaEl"] = kTOFnSigmaEl;
  fgVarNamesMap["kTOFnSigmaMu"] = kTOFnSigmaMu;
  fgVarNamesMap["kTOFnSigmaPi"] = kTOFnSigmaPi;
  fgVarNamesMap["kTOFnSigmaKa"] = kTOFnSigmaKa;
  fgVarNamesMap["kTOFnSigmaPr"] = kTOFnSigmaPr;
  fgVarNamesMap["kTrackTimeResIsRange"] = kTrackTimeResIsRange;
  fgVarNamesMap["kPVContributor"] = kPVContributor;
  fgVarNamesMap["kOrphanTrack"] = kOrphanTrack;
  fgVarNamesMap["kIsAmbiguous"] = kIsAmbiguous;
  fgVarNamesMap["kIsLegFromGamma"] = kIsLegFromGamma;
  fgVarNamesMap["kIsLegFromK0S"] = kIsLegFromK0S;
  fgVarNamesMap["kIsLegFromLambda"] = kIsLegFromLambda;
  fgVarNamesMap["kIsLegFromAntiLambda"] = kIsLegFromAntiLambda;
  fgVarNamesMap["kIsLegFromOmega"] = kIsLegFromOmega;
  fgVarNamesMap["kIsProtonFromLambdaAndAntiLambda"] = kIsProtonFromLambdaAndAntiLambda;
  fgVarNamesMap["kIsDalitzLeg"] = kIsDalitzLeg;
  fgVarNamesMap["kBarrelNAssocsInBunch"] = kBarrelNAssocsInBunch;
  fgVarNamesMap["kBarrelNAssocsOutOfBunch"] = kBarrelNAssocsOutOfBunch;
  fgVarNamesMap["kNBarrelTrackVariables"] = kNBarrelTrackVariables;
  fgVarNamesMap["kMuonNClusters"] = kMuonNClusters;
  fgVarNamesMap["kMuonPDca"] = kMuonPDca;
  fgVarNamesMap["kMuonRAtAbsorberEnd"] = kMuonRAtAbsorberEnd;
  fgVarNamesMap["kMCHBitMap"] = kMCHBitMap;
  fgVarNamesMap["kMuonChi2"] = kMuonChi2;
  fgVarNamesMap["kMuonChi2MatchMCHMID"] = kMuonChi2MatchMCHMID;
  fgVarNamesMap["kMuonChi2MatchMCHMFT"] = kMuonChi2MatchMCHMFT;
  fgVarNamesMap["kMuonMatchScoreMCHMFT"] = kMuonMatchScoreMCHMFT;
  fgVarNamesMap["kMuonCXX"] = kMuonCXX;
  fgVarNamesMap["kMuonCXY"] = kMuonCXY;
  fgVarNamesMap["kMuonCYY"] = kMuonCYY;
  fgVarNamesMap["kMuonCPhiX"] = kMuonCPhiX;
  fgVarNamesMap["kMuonCPhiY"] = kMuonCPhiY;
  fgVarNamesMap["kMuonCPhiPhi"] = kMuonCPhiPhi;
  fgVarNamesMap["kMuonCTglX"] = kMuonCTglX;
  fgVarNamesMap["kMuonCTglY"] = kMuonCTglY;
  fgVarNamesMap["kMuonCTglPhi"] = kMuonCTglPhi;
  fgVarNamesMap["kMuonCTglTgl"] = kMuonCTglTgl;
  fgVarNamesMap["kMuonC1Pt2X"] = kMuonC1Pt2X;
  fgVarNamesMap["kMuonC1Pt2Y"] = kMuonC1Pt2Y;
  fgVarNamesMap["kMuonC1Pt2Phi"] = kMuonC1Pt2Phi;
  fgVarNamesMap["kMuonC1Pt2Tgl"] = kMuonC1Pt2Tgl;
  fgVarNamesMap["kMuonC1Pt21Pt2"] = kMuonC1Pt21Pt2;
  fgVarNamesMap["kMuonTrackType"] = kMuonTrackType;
  fgVarNamesMap["kMuonDCAx"] = kMuonDCAx;
  fgVarNamesMap["kMuonDCAy"] = kMuonDCAy;
  fgVarNamesMap["kMuonTime"] = kMuonTime;
  fgVarNamesMap["kMuonTimeRes"] = kMuonTimeRes;
  fgVarNamesMap["kMftNClusters"] = kMftNClusters;
  fgVarNamesMap["kMftClusterSize"] = kMftClusterSize;
  fgVarNamesMap["kMftMeanClusterSize"] = kMftMeanClusterSize;
  fgVarNamesMap["kMuonNAssocsInBunch"] = kMuonNAssocsInBunch;
  fgVarNamesMap["kMuonNAssocsOutOfBunch"] = kMuonNAssocsOutOfBunch;
  fgVarNamesMap["kNMuonTrackVariables"] = kNMuonTrackVariables;
  fgVarNamesMap["kMCPdgCode"] = kMCPdgCode;
  fgVarNamesMap["kMCParticleWeight"] = kMCParticleWeight;
  fgVarNamesMap["kMCPx"] = kMCPx;
  fgVarNamesMap["kMCPy"] = kMCPy;
  fgVarNamesMap["kMCPz"] = kMCPz;
  fgVarNamesMap["kMCE"] = kMCE;
  fgVarNamesMap["kMCVx"] = kMCVx;
  fgVarNamesMap["kMCVy"] = kMCVy;
  fgVarNamesMap["kMCVz"] = kMCVz;
  fgVarNamesMap["kMCPt"] = kMCPt;
  fgVarNamesMap["kMCPhi"] = kMCPhi;
  fgVarNamesMap["kMCEta"] = kMCEta;
  fgVarNamesMap["kMCY"] = kMCY;
  fgVarNamesMap["kMCParticleGeneratorId"] = kMCParticleGeneratorId;
  fgVarNamesMap["kNMCParticleVariables"] = kNMCParticleVariables;
  fgVarNamesMap["kMCMotherPdgCode"] = kMCMotherPdgCode;
  fgVarNamesMap["kCandidateId"] = kCandidateId;
  fgVarNamesMap["kPairType"] = kPairType;
  fgVarNamesMap["kVertexingLxy"] = kVertexingLxy;
  fgVarNamesMap["kVertexingLxyErr"] = kVertexingLxyErr;
  fgVarNamesMap["kVertexingPseudoCTau"] = kVertexingPseudoCTau;
  fgVarNamesMap["kVertexingLxyz"] = kVertexingLxyz;
  fgVarNamesMap["kVertexingLxyzErr"] = kVertexingLxyzErr;
  fgVarNamesMap["kVertexingLz"] = kVertexingLz;
  fgVarNamesMap["kVertexingLzErr"] = kVertexingLzErr;
  fgVarNamesMap["kVertexingTauxy"] = kVertexingTauxy;
  fgVarNamesMap["kVertexingTauxyErr"] = kVertexingTauxyErr;
  fgVarNamesMap["kVertexingLzProjected"] = kVertexingLzProjected;
  fgVarNamesMap["kVertexingLxyProjected"] = kVertexingLxyProjected;
  fgVarNamesMap["kVertexingLxyzProjected"] = kVertexingLxyzProjected;
  fgVarNamesMap["kVertexingTauzProjected"] = kVertexingTauzProjected;
  fgVarNamesMap["kVertexingTauxyProjected"] = kVertexingTauxyProjected;
  fgVarNamesMap["kVertexingTauxyProjectedNs"] = kVertexingTauxyProjectedNs;
  fgVarNamesMap["kVertexingTauxyzProjected"] = kVertexingTauxyzProjected;
  fgVarNamesMap["kVertexingTauz"] = kVertexingTauz;
  fgVarNamesMap["kVertexingTauzErr"] = kVertexingTauzErr;
  fgVarNamesMap["kVertexingPz"] = kVertexingPz;
  fgVarNamesMap["kVertexingSV"] = kVertexingSV;
  fgVarNamesMap["kVertexingProcCode"] = kVertexingProcCode;
  fgVarNamesMap["kVertexingChi2PCA"] = kVertexingChi2PCA;
  fgVarNamesMap["kCosThetaHE"] = kCosThetaHE;
  fgVarNamesMap["kCosThetaCS"] = kCosThetaCS;
  fgVarNamesMap["kPhiHE"] = kPhiHE;
  fgVarNamesMap["kPhiCS"] = kPhiCS;
  fgVarNamesMap["kCosPhiVP"] = kCosPhiVP;
  fgVarNamesMap["kPhiVP"] = kPhiVP;
  fgVarNamesMap["kDeltaPhiPair2"] = kDeltaPhiPair2;
  fgVarNamesMap["kDeltaEtaPair2"] = kDeltaEtaPair2;
  fgVarNamesMap["kPsiPair"] = kPsiPair;
  fgVarNamesMap["kDeltaPhiPair"] = kDeltaPhiPair;
  fgVarNamesMap["kOpeningAngle"] = kOpeningAngle;
  fgVarNamesMap["kQuadDCAabsXY"] = kQuadDCAabsXY;
  fgVarNamesMap["kQuadDCAsigXY"] = kQuadDCAsigXY;
  fgVarNamesMap["kQuadDCAabsZ"] = kQuadDCAabsZ;
  fgVarNamesMap["kQuadDCAsigZ"] = kQuadDCAsigZ;
  fgVarNamesMap["kQuadDCAsigXYZ"] = kQuadDCAsigXYZ;
  fgVarNamesMap["kSignQuadDCAsigXY"] = kSignQuadDCAsigXY;
  fgVarNamesMap["kCosPointingAngle"] = kCosPointingAngle;
  fgVarNamesMap["kImpParXYJpsi"] = kImpParXYJpsi;
  fgVarNamesMap["kImpParXYK"] = kImpParXYK;
  fgVarNamesMap["kDCATrackProd"] = kDCATrackProd;
  fgVarNamesMap["kDCATrackVtxProd"] = kDCATrackVtxProd;
  fgVarNamesMap["kV2SP"] = kV2SP;
  fgVarNamesMap["kV2EP"] = kV2EP;
  fgVarNamesMap["kWV2SP"] = kWV2SP;
  fgVarNamesMap["kWV2EP"] = kWV2EP;
  fgVarNamesMap["kU2Q2"] = kU2Q2;
  fgVarNamesMap["kU3Q3"] = kU3Q3;
  fgVarNamesMap["kQ42XA"] = kQ42XA;
  fgVarNamesMap["kQ42YA"] = kQ42YA;
  fgVarNamesMap["kQ23XA"] = kQ23XA;
  fgVarNamesMap["kQ23YA"] = kQ23YA;
  fgVarNamesMap["kS11A"] = kS11A;
  fgVarNamesMap["kS12A"] = kS12A;
  fgVarNamesMap["kS13A"] = kS13A;
  fgVarNamesMap["kS31A"] = kS31A;
  fgVarNamesMap["kM11REF"] = kM11REF;
  fgVarNamesMap["kM11REFetagap"] = kM11REFetagap;
  fgVarNamesMap["kM01POI"] = kM01POI;
  fgVarNamesMap["kM1111REF"] = kM1111REF;
  fgVarNamesMap["kM11M1111REF"] = kM11M1111REF;
  fgVarNamesMap["kM11M1111REFoverMp"] = kM11M1111REFoverMp;
  fgVarNamesMap["kM01M0111POIoverMp"] = kM01M0111POIoverMp;
  fgVarNamesMap["kM0111POI"] = kM0111POI;
  fgVarNamesMap["kCORR2REF"] = kCORR2REF;
  fgVarNamesMap["kCORR2REFbydimuons"] = kCORR2REFbydimuons;
  fgVarNamesMap["kCORR2REFetagap"] = kCORR2REFetagap;
  fgVarNamesMap["kCORR2POI"] = kCORR2POI;
  fgVarNamesMap["kCORR2POICORR4POI"] = kCORR2POICORR4POI;
  fgVarNamesMap["kCORR2REFCORR4POI"] = kCORR2REFCORR4POI;
  fgVarNamesMap["kCORR2REFCORR2POI"] = kCORR2REFCORR2POI;
  fgVarNamesMap["kM01M0111overMp"] = kM01M0111overMp;
  fgVarNamesMap["kM11M0111overMp"] = kM11M0111overMp;
  fgVarNamesMap["kM11M01overMp"] = kM11M01overMp;
  fgVarNamesMap["kCORR2CORR4REF"] = kCORR2CORR4REF;
  fgVarNamesMap["kCORR4REF"] = kCORR4REF;
  fgVarNamesMap["kCORR4REFbydimuons"] = kCORR4REFbydimuons;
  fgVarNamesMap["kCORR4POI"] = kCORR4POI;
  fgVarNamesMap["kM11REFoverMp"] = kM11REFoverMp;
  fgVarNamesMap["kM01POIoverMp"] = kM01POIoverMp;
  fgVarNamesMap["kM1111REFoverMp"] = kM1111REFoverMp;
  fgVarNamesMap["kM0111POIoverMp"] = kM0111POIoverMp;
  fgVarNamesMap["kCORR2POIMp"] = kCORR2POIMp;
  fgVarNamesMap["kCORR4POIMp"] = kCORR4POIMp;
  fgVarNamesMap["kR2SP"] = kR2SP;
  fgVarNamesMap["kR2EP"] = kR2EP;
  fgVarNamesMap["kPsi2A"] = kPsi2A;
  fgVarNamesMap["kPsi2APOS"] = kPsi2APOS;
  fgVarNamesMap["kPsi2ANEG"] = kPsi2ANEG;
  fgVarNamesMap["kPsi2B"] = kPsi2B;
  fgVarNamesMap["kPsi2C"] = kPsi2C;
  fgVarNamesMap["kCos2DeltaPhi"] = kCos2DeltaPhi;
  fgVarNamesMap["kCos2DeltaPhiMu1"] = kCos2DeltaPhiMu1;
  fgVarNamesMap["kCos2DeltaPhiMu2"] = kCos2DeltaPhiMu2;
  fgVarNamesMap["kCos3DeltaPhi"] = kCos3DeltaPhi;
  fgVarNamesMap["kDeltaPtotTracks"] = kDeltaPtotTracks;
  fgVarNamesMap["kVertexingLxyOverErr"] = kVertexingLxyOverErr;
  fgVarNamesMap["kVertexingLzOverErr"] = kVertexingLzOverErr;
  fgVarNamesMap["kVertexingLxyzOverErr"] = kVertexingLxyzOverErr;
  fgVarNamesMap["kKFTrack0DCAxyz"] = kKFTrack0DCAxyz;
  fgVarNamesMap["kKFTrack1DCAxyz"] = kKFTrack1DCAxyz;
  fgVarNamesMap["kKFTracksDCAxyzMax"] = kKFTracksDCAxyzMax;
  fgVarNamesMap["kKFDCAxyzBetweenProngs"] = kKFDCAxyzBetweenProngs;
  fgVarNamesMap["kKFTrack0DCAxy"] = kKFTrack0DCAxy;
  fgVarNamesMap["kKFTrack1DCAxy"] = kKFTrack1DCAxy;
  fgVarNamesMap["kKFTracksDCAxyMax"] = kKFTracksDCAxyMax;
  fgVarNamesMap["kKFDCAxyBetweenProngs"] = kKFDCAxyBetweenProngs;
  fgVarNamesMap["kKFTrack0DeviationFromPV"] = kKFTrack0DeviationFromPV;
  fgVarNamesMap["kKFTrack1DeviationFromPV"] = kKFTrack1DeviationFromPV;
  fgVarNamesMap["kKFTrack0DeviationxyFromPV"] = kKFTrack0DeviationxyFromPV;
  fgVarNamesMap["kKFTrack1DeviationxyFromPV"] = kKFTrack1DeviationxyFromPV;
  fgVarNamesMap["kKFChi2OverNDFGeo"] = kKFChi2OverNDFGeo;
  fgVarNamesMap["kKFNContributorsPV"] = kKFNContributorsPV;
  fgVarNamesMap["kKFCosPA"] = kKFCosPA;
  fgVarNamesMap["kKFChi2OverNDFGeoTop"] = kKFChi2OverNDFGeoTop;
  fgVarNamesMap["kKFJpsiDCAxyz"] = kKFJpsiDCAxyz;
  fgVarNamesMap["kKFJpsiDCAxy"] = kKFJpsiDCAxy;
  fgVarNamesMap["kKFPairDeviationFromPV"] = kKFPairDeviationFromPV;
  fgVarNamesMap["kKFPairDeviationxyFromPV"] = kKFPairDeviationxyFromPV;
  fgVarNamesMap["kNPairVariables"] = kNPairVariables;
  fgVarNamesMap["kPairMass"] = kPairMass;
  fgVarNamesMap["kPairMassDau"] = kPairMassDau;
  fgVarNamesMap["kMassDau"] = kMassDau;
  fgVarNamesMap["kPairPt"] = kPairPt;
  fgVarNamesMap["kPairPtDau"] = kPairPtDau;
  fgVarNamesMap["kPairEta"] = kPairEta;
  fgVarNamesMap["kPairPhi"] = kPairPhi;
  fgVarNamesMap["kPairPhiv"] = kPairPhiv;
  fgVarNamesMap["kDeltaEta"] = kDeltaEta;
  fgVarNamesMap["kDeltaPhi"] = kDeltaPhi;
  fgVarNamesMap["kDeltaPhiSym"] = kDeltaPhiSym;
  fgVarNamesMap["kNCorrelationVariables"] = kNCorrelationVariables;
  fgVarNamesMap["kQuadMass"] = kQuadMass;
  fgVarNamesMap["kQuadDefaultDileptonMass"] = kQuadDefaultDileptonMass;
  fgVarNamesMap["kQuadPt"] = kQuadPt;
  fgVarNamesMap["kQuadEta"] = kQuadEta;
  fgVarNamesMap["kQuadPhi"] = kQuadPhi;
  fgVarNamesMap["kCosthetaDileptonDitrack"] = kCosthetaDileptonDitrack;
  fgVarNamesMap["kDitrackMass"] = kDitrackMass;
  fgVarNamesMap["kDitrackPt"] = kDitrackPt;
  fgVarNamesMap["kQ"] = kQ;
  fgVarNamesMap["kDeltaR1"] = kDeltaR1;
  fgVarNamesMap["kDeltaR2"] = kDeltaR2;
  fgVarNamesMap["kMassCharmHadron"] = kMassCharmHadron;
  fgVarNamesMap["kPtCharmHadron"] = kPtCharmHadron;
  fgVarNamesMap["kRapCharmHadron"] = kRapCharmHadron;
  fgVarNamesMap["kPhiCharmHadron"] = kPhiCharmHadron;
  fgVarNamesMap["kBdtCharmHadron"] = kBdtCharmHadron;
  fgVarNamesMap["kBitMapIndex"] = kBitMapIndex;
  fgVarNamesMap["kDeltaMass"] = kDeltaMass;
  fgVarNamesMap["kDeltaMass_jpsi"] = kDeltaMass_jpsi;
}
