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

#include "Tools/KFparticle/KFUtilities.h"

#include <cmath>
#include <iostream>
#include <map>
#include <vector>

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
float VarManager::fgzMatching = -77.5;
float VarManager::fgValues[VarManager::kNVars] = {0.0f};
float VarManager::fgTPCInterSectorBoundary = 1.0; // cm
int VarManager::fgITSROFbias = 0;
int VarManager::fgITSROFlength = 100;
int VarManager::fgITSROFBorderMarginLow = 0;
int VarManager::fgITSROFBorderMarginHigh = 0;
uint64_t VarManager::fgSOR = 0;
uint64_t VarManager::fgEOR = 0;
ROOT::Math::PxPyPzEVector VarManager::fgBeamA(0, 0, 6799.99, 6800);  // GeV, beam from A-side 4-momentum vector
ROOT::Math::PxPyPzEVector VarManager::fgBeamC(0, 0, -6799.99, 6800); // GeV, beam from C-side 4-momentum vector
o2::vertexing::DCAFitterN<2> VarManager::fgFitterTwoProngBarrel;
o2::vertexing::DCAFitterN<3> VarManager::fgFitterThreeProngBarrel;
o2::vertexing::DCAFitterN<4> VarManager::fgFitterFourProngBarrel;
o2::vertexing::FwdDCAFitterN<2> VarManager::fgFitterTwoProngFwd;
o2::vertexing::FwdDCAFitterN<3> VarManager::fgFitterThreeProngFwd;
o2::globaltracking::MatchGlobalFwd VarManager::mMatching;
std::map<VarManager::CalibObjects, TObject*> VarManager::fgCalibs;
bool VarManager::fgRunTPCPostCalibration[4] = {false, false, false, false};
int VarManager::fgCalibrationType = 0;                // 0 - no calibration, 1 - calibration vs (TPCncls,pIN,eta) typically for pp, 2 - calibration vs (eta,nPV,nLong,tLong) typically for PbPb
bool VarManager::fgUseInterpolatedCalibration = true; // use interpolated calibration histograms (default: true)

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
void VarManager::SetCollisionSystem(TString system, float energy)
{
  //
  // Set the collision system and the center of mass energy
  //
  int NumberOfNucleonsA = 1; // default value for pp collisions
  int NumberOfNucleonsC = 1; // default value for pp collisions
  int NumberOfProtonsA = 1;  // default value for pp collisions
  int NumberOfProtonsC = 1;  // default value for pp collisions
  if (system.EqualTo("PbPb")) {
    NumberOfNucleonsA = 208;
    NumberOfNucleonsC = 208;
    NumberOfProtonsA = 82; // Pb has 82 protons
    NumberOfProtonsC = 82; // Pb has 82 protons
  } else if (system.EqualTo("pp")) {
    NumberOfNucleonsA = 1;
    NumberOfNucleonsC = 1;
    NumberOfProtonsA = 1; // proton has 1 proton
    NumberOfProtonsC = 1; // proton has 1 proton
  } else if (system.EqualTo("XeXe")) {
    NumberOfNucleonsA = 129;
    NumberOfNucleonsC = 129;
    NumberOfProtonsA = 54; // Xe has 54 protons
    NumberOfProtonsC = 54; // Xe has 54 protons
  } else if (system.EqualTo("pPb")) {
    NumberOfNucleonsA = 1;
    NumberOfNucleonsC = 208;
    NumberOfProtonsA = 1;  // proton has 1 proton
    NumberOfProtonsC = 82; // Pb has 82 protons
  } else if (system.EqualTo("Pbp")) {
    NumberOfNucleonsA = 208;
    NumberOfNucleonsC = 1;
    NumberOfProtonsA = 82; // Pb has 82 protons
    NumberOfProtonsC = 1;  // proton has 1 proton
  } else if (system.EqualTo("OO")) {
    NumberOfNucleonsA = 16;
    NumberOfNucleonsC = 16;
    NumberOfProtonsA = 8; // O has 8 protons
    NumberOfProtonsC = 8; // O has 8 protons
  } else if (system.EqualTo("pO")) {
    NumberOfNucleonsA = 1;
    NumberOfNucleonsC = 16;
    NumberOfProtonsA = 1; // proton has 1 proton
    NumberOfProtonsC = 8; // O has 8 protons
  } else if (system.EqualTo("NeNe")) {
    NumberOfNucleonsA = 20;
    NumberOfNucleonsC = 20;
    NumberOfProtonsA = 10; // Ne has 5 protons
    NumberOfProtonsC = 10; // Ne has 5 protons
  }
  // TO Do: add more systems

  // set the beam 4-momentum vectors
  float beamAEnergy = energy / 2.0 * sqrt(NumberOfProtonsA * NumberOfProtonsC / NumberOfProtonsC / NumberOfProtonsA); // GeV
  float beamCEnergy = energy / 2.0 * sqrt(NumberOfProtonsC * NumberOfProtonsA / NumberOfProtonsA / NumberOfProtonsC); // GeV
  float beamAMomentum = std::sqrt(beamAEnergy * beamAEnergy - NumberOfNucleonsA * NumberOfNucleonsA * MassProton * MassProton);
  float beamCMomentum = std::sqrt(beamCEnergy * beamCEnergy - NumberOfNucleonsC * NumberOfNucleonsC * MassProton * MassProton);
  fgBeamA.SetPxPyPzE(0, 0, beamAMomentum, beamAEnergy);
  fgBeamC.SetPxPyPzE(0, 0, -beamCMomentum, beamCEnergy);
}

//__________________________________________________________________
void VarManager::SetCollisionSystem(o2::parameters::GRPLHCIFData* grplhcif)
{
  //
  // Set the collision system and the center of mass energy from the GRP information
  double beamAEnergy = grplhcif->getBeamEnergyPerNucleonInGeV(o2::constants::lhc::BeamDirection::BeamA);
  double beamCEnergy = grplhcif->getBeamEnergyPerNucleonInGeV(o2::constants::lhc::BeamDirection::BeamC);
  double beamANucleons = grplhcif->getBeamA(o2::constants::lhc::BeamDirection::BeamA);
  double beamCNucleons = grplhcif->getBeamA(o2::constants::lhc::BeamDirection::BeamC);
  double beamAMomentum = std::sqrt(beamAEnergy * beamAEnergy - beamANucleons * beamANucleons * MassProton * MassProton);
  double beamCMomentum = std::sqrt(beamCEnergy * beamCEnergy - beamCNucleons * beamCNucleons * MassProton * MassProton);
  fgBeamA.SetPxPyPzE(0, 0, beamAMomentum, beamAEnergy);
  fgBeamC.SetPxPyPzE(0, 0, -beamCMomentum, beamCEnergy);
}

//__________________________________________________________________
// void VarManager::FillEventDerived(float* values)
// {
//   //
//   // Fill event-wise derived quantities (these are all quantities which can be computed just based on the values already filled in the FillEvent() function)
//   //
// }

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
double VarManager::ComputePIDcalibration(int species, double nSigmaValue)
{
  // species: 0 - electron, 1 - pion, 2 - kaon, 3 - proton
  // Depending on the PID calibration type, we use different types of calibration histograms

  if (fgCalibrationType == 1) {
    // get the calibration histograms
    CalibObjects calibMean, calibSigma;
    switch (species) {
      case 0:
        calibMean = kTPCElectronMean;
        calibSigma = kTPCElectronSigma;
        break;
      case 1:
        calibMean = kTPCPionMean;
        calibSigma = kTPCPionSigma;
        break;
      case 2:
        calibMean = kTPCKaonMean;
        calibSigma = kTPCKaonSigma;
        break;
      case 3:
        calibMean = kTPCProtonMean;
        calibSigma = kTPCProtonSigma;
        break;
      default:
        LOG(fatal) << "Invalid species for PID calibration: " << species;
        return -999.0; // Return zero if species is invalid
    }

    TH3F* calibMeanHist = reinterpret_cast<TH3F*>(fgCalibs[calibMean]);
    TH3F* calibSigmaHist = reinterpret_cast<TH3F*>(fgCalibs[calibSigma]);
    if (!calibMeanHist || !calibSigmaHist) {
      LOG(fatal) << "Calibration histograms not found for species: " << species;
      return -999.0; // Return zero if histograms are not found
    }

    // Get the bin indices for the calibration histograms
    int binTPCncls = calibMeanHist->GetXaxis()->FindBin(fgValues[kTPCncls]);
    binTPCncls = (binTPCncls == 0 ? 1 : binTPCncls);
    binTPCncls = (binTPCncls > calibMeanHist->GetXaxis()->GetNbins() ? calibMeanHist->GetXaxis()->GetNbins() : binTPCncls);
    int binPin = calibMeanHist->GetYaxis()->FindBin(fgValues[kPin]);
    binPin = (binPin == 0 ? 1 : binPin);
    binPin = (binPin > calibMeanHist->GetYaxis()->GetNbins() ? calibMeanHist->GetYaxis()->GetNbins() : binPin);
    int binEta = calibMeanHist->GetZaxis()->FindBin(fgValues[kEta]);
    binEta = (binEta == 0 ? 1 : binEta);
    binEta = (binEta > calibMeanHist->GetZaxis()->GetNbins() ? calibMeanHist->GetZaxis()->GetNbins() : binEta);

    double mean = calibMeanHist->GetBinContent(binTPCncls, binPin, binEta);
    double sigma = calibSigmaHist->GetBinContent(binTPCncls, binPin, binEta);
    return (nSigmaValue - mean) / sigma; // Return the calibrated nSigma value
  } else if (fgCalibrationType == 2) {
    // get the calibration histograms
    CalibObjects calibMean, calibSigma, calibStatus;
    switch (species) {
      case 0:
        calibMean = kTPCElectronMean;
        calibSigma = kTPCElectronSigma;
        calibStatus = kTPCElectronStatus;
        break;
      case 1:
        calibMean = kTPCPionMean;
        calibSigma = kTPCPionSigma;
        calibStatus = kTPCPionStatus;
        break;
      case 2:
        calibMean = kTPCKaonMean;
        calibSigma = kTPCKaonSigma;
        calibStatus = kTPCKaonStatus;
        break;
      case 3:
        calibMean = kTPCProtonMean;
        calibSigma = kTPCProtonSigma;
        calibStatus = kTPCProtonStatus;
        break;
      default:
        LOG(fatal) << "Invalid species for PID calibration: " << species;
        return -999.0; // Return zero if species is invalid
    }

    THnF* calibMeanHist = reinterpret_cast<THnF*>(fgCalibs[calibMean]);
    THnF* calibSigmaHist = reinterpret_cast<THnF*>(fgCalibs[calibSigma]);
    THnF* calibStatusHist = reinterpret_cast<THnF*>(fgCalibs[calibStatus]);
    if (!calibMeanHist || !calibSigmaHist || !calibStatusHist) {
      LOG(fatal) << "Calibration histograms not found for species: " << species;
      return -999.0; // Return zero if histograms are not found
    }

    // Get the bin indices for the calibration histograms
    int binEta = calibMeanHist->GetAxis(0)->FindBin(fgValues[kEta]);
    binEta = (binEta == 0 ? 1 : binEta);
    binEta = (binEta > calibMeanHist->GetAxis(0)->GetNbins() ? calibMeanHist->GetAxis(0)->GetNbins() : binEta);
    int binNpv = calibMeanHist->GetAxis(1)->FindBin(fgValues[kVtxNcontribReal]);
    binNpv = (binNpv == 0 ? 1 : binNpv);
    binNpv = (binNpv > calibMeanHist->GetAxis(1)->GetNbins() ? calibMeanHist->GetAxis(1)->GetNbins() : binNpv);
    int binNlong = calibMeanHist->GetAxis(2)->FindBin(fgValues[kNTPCcontribLongA]);
    binNlong = (binNlong == 0 ? 1 : binNlong);
    binNlong = (binNlong > calibMeanHist->GetAxis(2)->GetNbins() ? calibMeanHist->GetAxis(2)->GetNbins() : binNlong);
    int binTlong = calibMeanHist->GetAxis(3)->FindBin(fgValues[kNTPCmedianTimeLongA]);
    binTlong = (binTlong == 0 ? 1 : binTlong);
    binTlong = (binTlong > calibMeanHist->GetAxis(3)->GetNbins() ? calibMeanHist->GetAxis(3)->GetNbins() : binTlong);

    int bin[4] = {binEta, binNpv, binNlong, binTlong};
    int status = static_cast<int>(calibStatusHist->GetBinContent(bin));
    double mean = calibMeanHist->GetBinContent(bin);
    double sigma = calibSigmaHist->GetBinContent(bin);
    switch (status) {
      case 0:
        // good calibration, return the calibrated nSigma value
        return (nSigmaValue - mean) / sigma;
        break;
      case 1:
        // calibration not valid, return the original nSigma value
        return nSigmaValue;
        break;
      case 2: // calibration constant has poor stat uncertainty, consider the user option for what to do
      case 3:
        // calibration constants have been interpolated
        if (fgUseInterpolatedCalibration) {
          return (nSigmaValue - mean) / sigma;
        } else {
          // return the original nSigma value
          return nSigmaValue;
        }
        break;
      case 4:
        // calibration constants interpolation failed, return the original nSigma value
        return nSigmaValue;
        break;
      default:
        return nSigmaValue; // unknown status, return the original nSigma value
        break;
    }
  } else {
    // unknown calibration type, return the original nSigma value
    LOG(fatal) << "Unknown calibration type: " << fgCalibrationType;
    return nSigmaValue; // Return the original nSigma value
  }
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
  fgVariableNames[kCentFT0A] = "Centrality FT0A";
  fgVariableUnits[kCentFT0A] = "%";
  fgVariableNames[kCentFT0M] = "Centrality FT0M";
  fgVariableUnits[kCentFT0M] = "%";
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
  fgVariableNames[kMultDimuonsME] = "Multiplicity Dimuons Unlike Sign Mixed Events";
  fgVariableUnits[kMultDimuonsME] = "";
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
  fgVariableNames[kMCEventCentrFT0C] = "MC Centrality FT0C";
  fgVariableNames[kMultMCNParticlesEta05] = "MC Multiplicity Central Barrel for |eta| < 0.5";
  fgVariableNames[kMultMCNParticlesEta08] = "MC Multiplicity Central Barrel for |eta| < 0.8";
  fgVariableNames[kMultMCNParticlesEta10] = "MC Multiplicity Central Barrel for |eta| < 1.0";
  fgVariableUnits[kMCEventGeneratorId] = "";
  fgVariableUnits[kMCEventSubGeneratorId] = "";
  fgVariableUnits[kMCVtxX] = "cm";
  fgVariableUnits[kMCVtxY] = "cm";
  fgVariableUnits[kMCVtxZ] = "cm";
  fgVariableUnits[kMCEventTime] = ""; // TODO: add proper unit
  fgVariableUnits[kMCEventWeight] = "";
  fgVariableUnits[kMCEventImpParam] = "b";
  fgVariableUnits[kMCEventCentrFT0C] = "%";
  fgVariableUnits[kMultMCNParticlesEta05] = "Multiplicity_eta05";
  fgVariableUnits[kMultMCNParticlesEta08] = "Multiplicity_eta08";
  fgVariableUnits[kMultMCNParticlesEta10] = "Multiplicity_eta10";
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
  fgVariableNames[kMultNTracksPVeta1] = "# Mult Tracks PV |#eta| < 1";
  fgVariableUnits[kMultNTracksPVeta1] = "";
  fgVariableNames[kMultNTracksPVetaHalf] = "# Mult Tracks PV |#eta| < 0.5";
  fgVariableUnits[kMultNTracksPVetaHalf] = "";
  fgVariableNames[kTrackOccupancyInTimeRange] = "track occupancy in TPC drift time (PV tracks)";
  fgVariableUnits[kTrackOccupancyInTimeRange] = "";
  fgVariableNames[kFT0COccupancyInTimeRange] = "FT0C occupancy";
  fgVariableUnits[kFT0COccupancyInTimeRange] = "";
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
  fgVariableNames[kPt1] = "p_{T1}";
  fgVariableUnits[kPt1] = "GeV/c";
  fgVariableNames[kPt2] = "p_{T2}";
  fgVariableUnits[kPt2] = "GeV/c";
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
  fgVariableNames[kCharge1] = "charge track 1";
  fgVariableUnits[kCharge1] = "";
  fgVariableNames[kCharge2] = "charge track 2";
  fgVariableUnits[kCharge2] = "";
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
  fgVariableNames[kMCCosTheta] = "Cos#theta";
  fgVariableUnits[kMCCosTheta] = "";
  fgVariableNames[kMCHadronPdgCode] = "HadronPdgCode";
  fgVariableUnits[kMCHadronPdgCode] = "";
  fgVariableNames[kMCCosChi] = "Cos#chi";
  fgVariableUnits[kMCCosChi] = "";
  fgVariableNames[kMCJpsiPt] = "Jpsi p_{T}";
  fgVariableUnits[kMCJpsiPt] = "GeV/c";
  fgVariableNames[kMCHadronPt] = "Hadron p_{T}";
  fgVariableUnits[kMCHadronPt] = "GeV/c";
  fgVariableNames[kMCHadronEta] = "Hadron #eta";
  fgVariableUnits[kMCHadronEta] = "";
  fgVariableNames[kMCdeltaphi] = "#Delta#phi";
  fgVariableUnits[kMCdeltaphi] = "";
  fgVariableNames[kMCdeltaeta] = "#Delta#eta";
  fgVariableUnits[kMCdeltaeta] = "";
  fgVariableNames[kMCParticleWeight] = "MC particle weight";
  fgVariableUnits[kMCParticleWeight] = "";
  fgVariableNames[kMCPx] = "MC px";
  fgVariableUnits[kMCPx] = "GeV/c";
  fgVariableNames[kMCPy] = "MC py";
  fgVariableUnits[kMCPy] = "GeV/c";
  fgVariableNames[kMCPz] = "MC pz";
  fgVariableUnits[kMCPz] = "GeV/c";
  fgVariableNames[kMCPt] = "MC p_{T}";
  fgVariableUnits[kMCPt] = "GeV/c";
  fgVariableNames[kMCPhi] = "#varphi";
  fgVariableUnits[kMCPhi] = "rad";
  fgVariableNames[kMCEta] = "MC #eta";
  fgVariableUnits[kMCEta] = "";
  fgVariableNames[kMCY] = "MC y";
  fgVariableUnits[kMCY] = "";
  fgVariableNames[kMCE] = "MC Energy";
  fgVariableUnits[kMCE] = "GeV";
  fgVariableNames[kMCMass] = "MC Mass";
  fgVariableUnits[kMCMass] = "GeV/c2";
  fgVariableNames[kMCVx] = "MC vx";
  fgVariableUnits[kMCVx] = "cm"; // TODO: check the unit
  fgVariableNames[kMCVy] = "MC vy";
  fgVariableUnits[kMCVy] = "cm"; // TODO: check the unit
  fgVariableNames[kMCVz] = "MC vz";
  fgVariableUnits[kMCVz] = "cm"; // TODO: check the unit
  fgVariableNames[kMCCosThetaHE] = "MC cos(#theta_{HE})";
  fgVariableUnits[kMCCosThetaHE] = "";
  fgVariableNames[kMCPhiHE] = "MC #varphi_{HE}";
  fgVariableUnits[kMCPhiHE] = "rad";
  fgVariableNames[kMCPhiTildeHE] = "MC #tilde{#varphi}_{HE}";
  fgVariableUnits[kMCPhiTildeHE] = "rad";
  fgVariableNames[kMCCosThetaCS] = "MC cos(#theta_{CS})";
  fgVariableUnits[kMCCosThetaCS] = "";
  fgVariableNames[kMCPhiCS] = "MC #varphi_{CS}";
  fgVariableUnits[kMCPhiCS] = "rad";
  fgVariableNames[kMCPhiTildeCS] = "MC #tilde{#varphi}_{CS}";
  fgVariableUnits[kMCPhiTildeCS] = "rad";
  fgVariableNames[kMCCosThetaPP] = "MC cos(#theta_{PP})";
  fgVariableUnits[kMCCosThetaPP] = "";
  fgVariableNames[kMCPhiPP] = "MC #varphi_{PP}";
  fgVariableUnits[kMCPhiPP] = "rad";
  fgVariableNames[kMCPhiTildePP] = "MC #tilde{#varphi}_{PP}";
  fgVariableUnits[kMCPhiTildePP] = "rad";
  fgVariableNames[kMCCosThetaRM] = "MC cos(#theta_{RM})";
  fgVariableUnits[kMCCosThetaRM] = "";
  fgVariableNames[kCandidateId] = "";
  fgVariableUnits[kCandidateId] = "";
  fgVariableNames[kPairType] = "Pair type";
  fgVariableUnits[kPairType] = "";
  fgVariableNames[kVertexingLxy] = "Pair Lxy";
  fgVariableUnits[kVertexingLxy] = "cm";
  fgVariableNames[kMCVertexingLxy] = "MC Lxy";
  fgVariableUnits[kMCVertexingLxy] = "cm";
  fgVariableNames[kVertexingLz] = "Pair Lz";
  fgVariableUnits[kVertexingLz] = "cm";
  fgVariableNames[kMCVertexingLz] = "MC Lz";
  fgVariableUnits[kMCVertexingLz] = "cm";
  fgVariableNames[kVertexingLxyz] = "Pair Lxyz";
  fgVariableUnits[kVertexingLxyz] = "cm";
  fgVariableNames[kMCVertexingLxyz] = "MC Lxyz";
  fgVariableUnits[kMCVertexingLxyz] = "cm";
  fgVariableNames[kMCLxyExpected] = "MC Expected Lxy";
  fgVariableUnits[kMCLxyExpected] = "cm";
  fgVariableNames[kMCLxyzExpected] = "MC Expected Lxyz";
  fgVariableUnits[kMCLxyzExpected] = "cm";
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
  fgVariableNames[kMCVertexingTauz] = "MC pseudo-proper Tauz";
  fgVariableUnits[kMCVertexingTauz] = "ns";
  fgVariableNames[kMCVertexingTauxy] = "MC pseudo-proper Tauxy";
  fgVariableUnits[kMCVertexingTauxy] = "ns";
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
  fgVariableNames[kVertexingTauxyProjectedPoleJPsiMass] = "Pair pseudo-proper Tauxy (with pole JPsi mass)";
  fgVariableUnits[kVertexingTauxyProjectedPoleJPsiMass] = "ns";
  fgVariableNames[kVertexingTauxyzProjected] = "Pair pseudo-proper Tauxyz";
  fgVariableUnits[kVertexingTauxyzProjected] = "ns";
  fgVariableNames[kMCVertexingLxyProjected] = "MC Lxy_{proj}";
  fgVariableUnits[kMCVertexingLxyProjected] = "cm";
  fgVariableNames[kMCVertexingLzProjected] = "MC Lz_{proj}";
  fgVariableUnits[kMCVertexingLzProjected] = "cm";
  fgVariableNames[kMCVertexingLxyzProjected] = "MC Lxyz_{proj}";
  fgVariableUnits[kMCVertexingLxyzProjected] = "cm";
  fgVariableNames[kMCVertexingTauzProjected] = "MC Tauz_{proj}";
  fgVariableUnits[kMCVertexingTauzProjected] = "ns";
  fgVariableNames[kMCVertexingTauxyProjected] = "MC Tauxy_{proj}";
  fgVariableUnits[kMCVertexingTauxyProjected] = "ns";
  fgVariableNames[kMCVertexingTauxyzProjected] = "MC Tauxyz_{proj}";
  fgVariableUnits[kMCVertexingTauxyzProjected] = "ns";
  fgVariableNames[kCosPointingAngle] = "cos(#theta_{pointing})";
  fgVariableUnits[kCosPointingAngle] = "";
  fgVariableNames[kMCCosPointingAngle] = "MC cos(#theta_{pointing})";
  fgVariableUnits[kMCCosPointingAngle] = "";
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
  fgVariableNames[kMultMuons] = "Multiplicity muons";
  fgVariableUnits[kMultMuons] = "";
  fgVariableNames[kMultAntiMuons] = "Multiplicity anti-muons";
  fgVariableUnits[kMultAntiMuons] = "";
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
  fgVariableNames[kM11REFoverMpminus] = "M^{-}_{11}^{REF}/M^{-}_{p} ";
  fgVariableUnits[kM11REFoverMpminus] = "";
  fgVariableNames[kM11REFoverMpplus] = "M^{+}_{11}^{REF}/M^{+}_{p} ";
  fgVariableUnits[kM11REFoverMpplus] = "";
  fgVariableNames[kM1111REFoverMpplus] = "M^{+}_{1111}^{REF}/M^{+}_{p} ";
  fgVariableUnits[kM1111REFoverMpplus] = "";
  fgVariableNames[kM1111REFoverMpminus] = "M^{-}_{1111}^{REF}/M^{-}_{p} ";
  fgVariableUnits[kM1111REFoverMpminus] = "";
  fgVariableNames[kM01POIME] = "M_{01}^{POI, ME}";
  fgVariableUnits[kM01POIME] = "";
  fgVariableNames[kM0111POIME] = "M_{0111}^{POI, ME}";
  fgVariableUnits[kM0111POIME] = "";
  fgVariableNames[kCORR2POIME] = "CORR2^{POI, ME}";
  fgVariableUnits[kCORR2POIME] = "";
  fgVariableNames[kCORR4POIME] = "CORR4^{POI, ME}";
  fgVariableUnits[kCORR4POIME] = "";
  fgVariableNames[kM01POIoverMpME] = "M_{01}^{POI, ME} / M_p";
  fgVariableUnits[kM01POIoverMpME] = "";
  fgVariableNames[kM0111POIoverMpME] = "M_{0111}^{POI, ME} / M_p";
  fgVariableUnits[kM0111POIoverMpME] = "";
  fgVariableNames[kM11REFoverMpME] = "M_{11}^{REF} / M_p";
  fgVariableUnits[kM11REFoverMpME] = "";
  fgVariableNames[kM1111REFoverMpME] = "M_{1111}^{REF} / M_p";
  fgVariableUnits[kM1111REFoverMpME] = "";
  fgVariableNames[kCORR2REFbydimuonsME] = "CORR2^{REF} / dimuons ME";
  fgVariableUnits[kCORR2REFbydimuonsME] = "";
  fgVariableNames[kCORR4REFbydimuonsME] = "CORR4^{REF} / dimuons ME";
  fgVariableUnits[kCORR4REFbydimuonsME] = "";
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
  fgVariableNames[kPairRap] = "#rap";
  fgVariableUnits[kPairRap] = "";
  fgVariableNames[kPairPhi] = "#varphi";
  fgVariableUnits[kPairPhi] = "rad.";
  fgVariableNames[kPairPhiv] = "#varphi_{V}";
  fgVariableUnits[kPairPhiv] = "rad.";
  fgVariableNames[kDileptonHadronKstar] = "Dilepton-hadron k^{*}";
  fgVariableUnits[kDileptonHadronKstar] = "GeV/c^{2}";
  fgVariableNames[kDeltaEta] = "#Delta#eta";
  fgVariableUnits[kDeltaEta] = "";
  fgVariableNames[kDeltaPhi] = "#Delta#phi";
  fgVariableUnits[kDeltaPhi] = "rad.";
  fgVariableNames[kDeltaPhiSym] = "#Delta#phi";
  fgVariableUnits[kDeltaPhiSym] = "rad.";
  fgVariableNames[kCosChi] = "Cos#chi";
  fgVariableUnits[kCosChi] = "";
  fgVariableNames[kCosTheta] = "Cos#theta";
  fgVariableUnits[kCosTheta] = "";
  fgVariableNames[kPtDau] = "hadron P_{T}";
  fgVariableUnits[kPtDau] = "GeV/c";
  fgVariableNames[kEtaDau] = "hadron #eta";
  fgVariableUnits[kEtaDau] = "";
  fgVariableNames[kPhiDau] = "hadron #phi";
  fgVariableUnits[kPhiDau] = "";
  fgVariableNames[kCosThetaHE] = "cos#it{#theta}";
  fgVariableUnits[kCosThetaHE] = "";
  fgVariableNames[kPhiHE] = "#varphi_{HE}";
  fgVariableUnits[kPhiHE] = "rad.";
  fgVariableNames[kPhiTildeHE] = "#tilde{#varphi}_{HE}";
  fgVariableUnits[kPhiTildeHE] = "rad.";
  fgVariableNames[kCosThetaCS] = "cos#it{#theta}_{CS}";
  fgVariableUnits[kCosThetaCS] = "";
  fgVariableNames[kPhiCS] = "#varphi_{CS}";
  fgVariableUnits[kPhiCS] = "rad.";
  fgVariableNames[kPhiTildeCS] = "#tilde{#varphi}_{CS}";
  fgVariableUnits[kPhiTildeCS] = "rad.";
  fgVariableNames[kCosThetaPP] = "cos#it{#theta}_{PP}";
  fgVariableUnits[kCosThetaPP] = "";
  fgVariableNames[kPhiPP] = "#varphi_{PP}";
  fgVariableUnits[kPhiPP] = "rad.";
  fgVariableNames[kPhiTildePP] = "#tilde{#varphi}_{PP}";
  fgVariableUnits[kPhiTildePP] = "rad.";
  fgVariableNames[kCosThetaRM] = "cos#it{#theta}_{RM}";
  fgVariableUnits[kCosThetaRM] = "";
  fgVariableNames[kCosThetaStarTPC] = "cos#it{#theta}^{*}_{TPC}";
  fgVariableUnits[kCosThetaStarTPC] = "";
  fgVariableNames[kCosThetaStarFT0A] = "cos#it{#theta}^{*}_{FT0A}";
  fgVariableUnits[kCosThetaStarFT0A] = "";
  fgVariableNames[kCosThetaStarFT0C] = "cos#it{#theta}^{*}_{FT0C}";
  fgVariableUnits[kCosThetaStarFT0C] = "";
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
  fgVariableNames[kV22m] = "v_{2}(2)_{#mu^{-}}";
  fgVariableUnits[kV22m] = "";
  fgVariableNames[kV24m] = "v_{2}(4)_{#mu^{-}}";
  fgVariableUnits[kV24m] = "";
  fgVariableNames[kV22p] = "v_{2}(2)_{#mu^{+}}";
  fgVariableUnits[kV22p] = "";
  fgVariableNames[kV24p] = "v_{2}(4)_{#mu^{+}}";
  fgVariableUnits[kV24p] = "";
  fgVariableNames[kV22ME] = "v_{2}(2)_{ME}";
  fgVariableUnits[kV22ME] = "";
  fgVariableNames[kV24ME] = "v_{2}(4)_{ME}";
  fgVariableUnits[kV24ME] = "";
  fgVariableNames[kWV22ME] = "W_{2}(2)_{ME}";
  fgVariableUnits[kWV22ME] = "";
  fgVariableNames[kWV24ME] = "W_{2}(4)_{ME}";
  fgVariableUnits[kWV24ME] = "";
  fgVariableNames[kS12] = "m_{12}^{2}";
  fgVariableUnits[kS12] = "GeV^{2}/c^{4}";
  fgVariableNames[kS13] = "m_{13}^{2}";
  fgVariableUnits[kS13] = "GeV^{2}/c^{4}";
  fgVariableNames[kS23] = "m_{23}^{2}";
  fgVariableUnits[kS23] = "GeV^{2}/c^{4}";
  fgVariableNames[kBdtBackground] = "kBdtBackground";
  fgVariableUnits[kBdtBackground] = " ";
  fgVariableNames[kBdtPrompt] = "kBdtPrompt";
  fgVariableUnits[kBdtPrompt] = " ";
  fgVariableNames[kBdtNonprompt] = "kBdtNonprompt";
  fgVariableUnits[kBdtNonprompt] = " ";

  // Set the variables short names map. This is needed for dynamic configuration via JSON files
  fgVarNamesMap["kNothing"] = kNothing;
  fgVarNamesMap["kRunNo"] = kRunNo;
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
  fgVarNamesMap["kCentFT0A"] = kCentFT0A;
  fgVarNamesMap["kCentFT0M"] = kCentFT0M;
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
  fgVarNamesMap["kMultDimuonsME"] = kMultDimuonsME;
  fgVarNamesMap["kMultNTracksHasITS"] = kMultNTracksHasITS;
  fgVarNamesMap["kMultNTracksHasTPC"] = kMultNTracksHasTPC;
  fgVarNamesMap["kMultNTracksHasTOF"] = kMultNTracksHasTOF;
  fgVarNamesMap["kMultNTracksHasTRD"] = kMultNTracksHasTRD;
  fgVarNamesMap["kMultNTracksITSOnly"] = kMultNTracksITSOnly;
  fgVarNamesMap["kMultNTracksTPCOnly"] = kMultNTracksTPCOnly;
  fgVarNamesMap["kMultNTracksITSTPC"] = kMultNTracksITSTPC;
  fgVarNamesMap["kMultNTracksPVeta1"] = kMultNTracksPVeta1;
  fgVarNamesMap["kMultNTracksPVetaHalf"] = kMultNTracksPVetaHalf;
  fgVarNamesMap["kTrackOccupancyInTimeRange"] = kTrackOccupancyInTimeRange;
  fgVarNamesMap["kFT0COccupancyInTimeRange"] = kFT0COccupancyInTimeRange;
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
  fgVarNamesMap["kTPCnSigmaMu"] = kTPCnSigmaMu;
  fgVarNamesMap["kTPCnSigmaPi"] = kTPCnSigmaPi;
  fgVarNamesMap["kTPCnSigmaKa"] = kTPCnSigmaKa;
  fgVarNamesMap["kTPCnSigmaPr"] = kTPCnSigmaPr;
  fgVarNamesMap["kTPCnSigmaEl_Corr"] = kTPCnSigmaEl_Corr;
  fgVarNamesMap["kTPCnSigmaPi_Corr"] = kTPCnSigmaPi_Corr;
  fgVarNamesMap["kTPCnSigmaKa_Corr"] = kTPCnSigmaKa_Corr;
  fgVarNamesMap["kTPCnSigmaPr_Corr"] = kTPCnSigmaPr_Corr;
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
  fgVarNamesMap["kMCCosTheta"] = kMCCosTheta;
  fgVarNamesMap["kMCHadronPdgCode"] = kMCHadronPdgCode;
  fgVarNamesMap["kMCCosChi"] = kMCCosChi;
  fgVarNamesMap["kMCHadronPt"] = kMCHadronPt;
  fgVarNamesMap["kMCWeight_before"] = kMCWeight_before;
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
  fgVarNamesMap["kMCCosThetaHE"] = kMCCosThetaHE;
  fgVarNamesMap["kMCPhiHE"] = kMCPhiHE;
  fgVarNamesMap["kMCPhiTildeHE"] = kMCPhiTildeHE;
  fgVarNamesMap["kMCCosThetaCS"] = kMCCosThetaCS;
  fgVarNamesMap["kMCPhiCS"] = kMCPhiCS;
  fgVarNamesMap["kMCPhiTildeCS"] = kMCPhiTildeCS;
  fgVarNamesMap["kMCCosThetaPP"] = kMCCosThetaPP;
  fgVarNamesMap["kMCPhiPP"] = kMCPhiPP;
  fgVarNamesMap["kMCPhiTildePP"] = kMCPhiTildePP;
  fgVarNamesMap["kMCCosThetaRM"] = kMCCosThetaRM;
  fgVarNamesMap["kMCParticleGeneratorId"] = kMCParticleGeneratorId;
  fgVarNamesMap["kNMCParticleVariables"] = kNMCParticleVariables;
  fgVarNamesMap["kMCMotherPdgCode"] = kMCMotherPdgCode;
  fgVarNamesMap["kCandidateId"] = kCandidateId;
  fgVarNamesMap["kPairType"] = kPairType;
  fgVarNamesMap["kVertexingLxy"] = kVertexingLxy;
  fgVarNamesMap["kVertexingLxyErr"] = kVertexingLxyErr;
  fgVarNamesMap["kMCVertexingLxy"] = kMCVertexingLxy;
  fgVarNamesMap["kVertexingPseudoCTau"] = kVertexingPseudoCTau;
  fgVarNamesMap["kVertexingLxyz"] = kVertexingLxyz;
  fgVarNamesMap["kVertexingLxyzErr"] = kVertexingLxyzErr;
  fgVarNamesMap["kMCVertexingLxyz"] = kMCVertexingLxyz;
  fgVarNamesMap["kVertexingLz"] = kVertexingLz;
  fgVarNamesMap["kVertexingLzErr"] = kVertexingLzErr;
  fgVarNamesMap["kMCVertexingLz"] = kMCVertexingLz;
  fgVarNamesMap["kVertexingTauxy"] = kVertexingTauxy;
  fgVarNamesMap["kVertexingTauxyErr"] = kVertexingTauxyErr;
  fgVarNamesMap["kMCVertexingTauxy"] = kMCVertexingTauxy;
  fgVarNamesMap["kVertexingLzProjected"] = kVertexingLzProjected;
  fgVarNamesMap["kVertexingLxyProjected"] = kVertexingLxyProjected;
  fgVarNamesMap["kVertexingLxyzProjected"] = kVertexingLxyzProjected;
  fgVarNamesMap["kVertexingTauzProjected"] = kVertexingTauzProjected;
  fgVarNamesMap["kVertexingTauxyProjected"] = kVertexingTauxyProjected;
  fgVarNamesMap["kVertexingTauxyProjectedPoleJPsiMass"] = kVertexingTauxyProjectedPoleJPsiMass;
  fgVarNamesMap["kVertexingTauxyProjectedNs"] = kVertexingTauxyProjectedNs;
  fgVarNamesMap["kVertexingTauxyzProjected"] = kVertexingTauxyzProjected;
  fgVarNamesMap["kMCVertexingTauzProjected"] = kVertexingTauzProjected;
  fgVarNamesMap["kMCVertexingTauxyProjected"] = kVertexingTauxyProjected;
  fgVarNamesMap["kMCVertexingTauxyzProjected"] = kVertexingTauxyzProjected;
  fgVarNamesMap["kVertexingTauz"] = kVertexingTauz;
  fgVarNamesMap["kVertexingTauzErr"] = kVertexingTauzErr;
  fgVarNamesMap["kMCVertexingTauz"] = kMCVertexingTauz;
  fgVarNamesMap["kVertexingPz"] = kVertexingPz;
  fgVarNamesMap["kVertexingSV"] = kVertexingSV;
  fgVarNamesMap["kVertexingProcCode"] = kVertexingProcCode;
  fgVarNamesMap["kVertexingChi2PCA"] = kVertexingChi2PCA;
  fgVarNamesMap["kCosThetaHE"] = kCosThetaHE;
  fgVarNamesMap["kPhiHE"] = kPhiHE;
  fgVarNamesMap["kPhiTildeHE"] = kPhiTildeHE;
  fgVarNamesMap["kCosThetaCS"] = kCosThetaCS;
  fgVarNamesMap["kPhiCS"] = kPhiCS;
  fgVarNamesMap["kPhiTildeCS"] = kPhiTildeCS;
  fgVarNamesMap["kCosThetaPP"] = kCosThetaPP;
  fgVarNamesMap["kPhiPP"] = kPhiPP;
  fgVarNamesMap["kPhiTildePP"] = kPhiTildePP;
  fgVarNamesMap["kCosThetaRM"] = kCosThetaRM;
  fgVarNamesMap["kCosThetaStarTPC"] = kCosThetaStarTPC;
  fgVarNamesMap["kCosThetaStarFT0A"] = kCosThetaStarFT0A;
  fgVarNamesMap["kCosThetaStarFT0C"] = kCosThetaStarFT0C;
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
  fgVarNamesMap["kMCCosPointingAngle"] = kMCCosPointingAngle;
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
  fgVarNamesMap["kMultAntiMuons"] = kMultAntiMuons;
  fgVarNamesMap["kMultMuons"] = kMultMuons;
  fgVarNamesMap["kM01POIplus"] = kM01POIplus;
  fgVarNamesMap["kM0111POIplus"] = kM0111POIplus;
  fgVarNamesMap["kM01POIminus"] = kM01POIminus;
  fgVarNamesMap["kM0111POIminus"] = kM0111POIminus;
  fgVarNamesMap["kM01POIoverMpminus"] = kM01POIoverMpminus;
  fgVarNamesMap["kM01POIoverMpplus"] = kM01POIoverMpplus;
  fgVarNamesMap["kM01POIoverMpmoins"] = kM01POIoverMpmoins;
  fgVarNamesMap["kM0111POIoverMpminus"] = kM0111POIoverMpminus;
  fgVarNamesMap["kM0111POIoverMpplus"] = kM0111POIoverMpplus;
  fgVarNamesMap["kCORR2POIplus"] = kCORR2POIplus;
  fgVarNamesMap["kCORR2POIminus"] = kCORR2POIminus;
  fgVarNamesMap["kCORR4POIplus"] = kCORR4POIplus;
  fgVarNamesMap["kCORR4POIminus"] = kCORR4POIminus;
  fgVarNamesMap["kM11REFoverMpminus"] = kM11REFoverMpminus;
  fgVarNamesMap["kM11REFoverMpplus"] = kM11REFoverMpplus;
  fgVarNamesMap["kM1111REFoverMpplus"] = kM1111REFoverMpplus;
  fgVarNamesMap["kM1111REFoverMpminus"] = kM1111REFoverMpminus;
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
  fgVarNamesMap["kM01POIME"] = kM01POIME;
  fgVarNamesMap["kMultDimuonsME"] = kMultDimuonsME;
  fgVarNamesMap["kM0111POIME"] = kM0111POIME;
  fgVarNamesMap["kCORR2POIME"] = kCORR2POIME;
  fgVarNamesMap["kCORR4POIME"] = kCORR4POIME;
  fgVarNamesMap["kM01POIoverMpME"] = kM01POIoverMpME;
  fgVarNamesMap["kM0111POIoverMpME"] = kM0111POIoverMpME;
  fgVarNamesMap["kM11REFoverMpME"] = kM11REFoverMpME;
  fgVarNamesMap["kM1111REFoverMpME"] = kM1111REFoverMpME;
  fgVarNamesMap["kCORR2REFbydimuonsME"] = kCORR2REFbydimuonsME;
  fgVarNamesMap["kCORR4REFbydimuonsME"] = kCORR4REFbydimuonsME;
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
  fgVarNamesMap["kS12"] = kS12,
  fgVarNamesMap["kS13"] = kS13,
  fgVarNamesMap["kS23"] = kS23,
  fgVarNamesMap["kNPairVariables"] = kNPairVariables;
  fgVarNamesMap["kPairMass"] = kPairMass;
  fgVarNamesMap["kPairMassDau"] = kPairMassDau;
  fgVarNamesMap["kMassDau"] = kMassDau;
  fgVarNamesMap["kPairPt"] = kPairPt;
  fgVarNamesMap["kPairPtDau"] = kPairPtDau;
  fgVarNamesMap["kPairEta"] = kPairEta;
  fgVarNamesMap["kPairRap"] = kPairRap;
  fgVarNamesMap["kPairPhi"] = kPairPhi;
  fgVarNamesMap["kPairPhiv"] = kPairPhiv;
  fgVarNamesMap["kDileptonHadronKstar"] = kDileptonHadronKstar;
  fgVarNamesMap["kDeltaEta"] = kDeltaEta;
  fgVarNamesMap["kDeltaPhi"] = kDeltaPhi;
  fgVarNamesMap["kDeltaPhiSym"] = kDeltaPhiSym;
  fgVarNamesMap["kCosTheta"] = kCosTheta;
  fgVarNamesMap["kCosChi"] = kCosChi;
  fgVarNamesMap["kECWeight"] = kECWeight;
  fgVarNamesMap["kEWeight_before"] = kEWeight_before;
  fgVarNamesMap["kPtDau"] = kPtDau;
  fgVarNamesMap["kEtaDau"] = kEtaDau;
  fgVarNamesMap["kPhiDau"] = kPhiDau;
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
  fgVarNamesMap["kV22m"] = kV22m;
  fgVarNamesMap["kV24m"] = kV24m;
  fgVarNamesMap["kV22p"] = kV22p;
  fgVarNamesMap["kV24p"] = kV24p;
  fgVarNamesMap["kV22ME"] = kV22ME;
  fgVarNamesMap["kV24ME"] = kV24ME;
  fgVarNamesMap["kWV22ME"] = kWV22ME;
  fgVarNamesMap["kWV24ME"] = kWV24ME;
  fgVarNamesMap["kBdtBackground"] = kBdtBackground;
  fgVarNamesMap["kBdtPrompt"] = kBdtPrompt;
  fgVarNamesMap["kBdtNonprompt"] = kBdtNonprompt;
  fgVariableNames[kAmplitudeFT0A] = "FT0A amplitude";
  fgVariableUnits[kAmplitudeFT0A] = "a.u.";
  fgVariableNames[kAmplitudeFT0C] = "FT0C amplitude";
  fgVariableUnits[kAmplitudeFT0C] = "a.u.";
  fgVariableNames[kTimeFT0A] = "FT0A time";
  fgVariableUnits[kTimeFT0A] = "ns";
  fgVariableNames[kTimeFT0C] = "FT0C time";
  fgVariableUnits[kTimeFT0C] = "ns";
  fgVariableNames[kTriggerMaskFT0] = "FT0 trigger mask";
  fgVariableUnits[kTriggerMaskFT0] = "";
  fgVariableNames[kNFiredChannelsFT0A] = "FT0A fired channels";
  fgVariableUnits[kNFiredChannelsFT0A] = "";
  fgVariableNames[kNFiredChannelsFT0C] = "FT0C fired channels";
  fgVariableUnits[kNFiredChannelsFT0C] = "";
  fgVariableNames[kAmplitudeFDDA] = "FDDA amplitude";
  fgVariableUnits[kAmplitudeFDDA] = "a.u.";
  fgVariableNames[kAmplitudeFDDC] = "FDDC amplitude";
  fgVariableUnits[kAmplitudeFDDC] = "a.u.";
  fgVariableNames[kTimeFDDA] = "FDDA time";
  fgVariableUnits[kTimeFDDA] = "ns";
  fgVariableNames[kTimeFDDC] = "FDDC time";
  fgVariableUnits[kTimeFDDC] = "ns";
  fgVariableNames[kTriggerMaskFDD] = "FDD trigger mask";
  fgVariableUnits[kTriggerMaskFDD] = "";
  fgVariableNames[kAmplitudeFV0A] = "FV0A amplitude";
  fgVariableUnits[kAmplitudeFV0A] = "a.u.";
  fgVariableNames[kTimeFV0A] = "FV0A time";
  fgVariableUnits[kTimeFV0A] = "ns";
  fgVariableNames[kTriggerMaskFV0A] = "FV0A trigger mask";
  fgVariableUnits[kTriggerMaskFV0A] = "";
  fgVariableNames[kNFiredChannelsFV0A] = "FV0A fired channels";
  fgVariableUnits[kNFiredChannelsFV0A] = "";
  fgVariableNames[kBBFT0Apf] = "FT0A BB pileup flag";
  fgVariableUnits[kBBFT0Apf] = "";
  fgVariableNames[kBGFT0Apf] = "FT0A BG pileup flag";
  fgVariableUnits[kBGFT0Apf] = "";
  fgVariableNames[kBBFT0Cpf] = "FT0C BB pileup flag";
  fgVariableUnits[kBBFT0Cpf] = "";
  fgVariableNames[kBGFT0Cpf] = "FT0C BG pileup flag";
  fgVariableUnits[kBGFT0Cpf] = "";
  fgVariableNames[kBBFV0Apf] = "FV0A BB pileup flag";
  fgVariableUnits[kBBFV0Apf] = "";
  fgVariableNames[kBGFV0Apf] = "FV0A BG pileup flag";
  fgVariableUnits[kBGFV0Apf] = "";
  fgVariableNames[kBBFDDApf] = "FDDA BB pileup flag";
  fgVariableUnits[kBBFDDApf] = "";
  fgVariableNames[kBGFDDApf] = "FDDA BG pileup flag";
  fgVariableUnits[kBGFDDApf] = "";
  fgVariableNames[kBBFDDCpf] = "FDDC BB pileup flag";
  fgVariableUnits[kBBFDDCpf] = "";
  fgVariableNames[kBGFDDCpf] = "FDDC BG pileup flag";
  fgVariableUnits[kBGFDDCpf] = "";
}
