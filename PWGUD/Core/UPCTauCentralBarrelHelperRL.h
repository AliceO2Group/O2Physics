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
///
/// \brief
/// \author Roman Lavicka, roman.lavicka@cern.ch
/// \since  27.10.2022

#ifndef PWGUD_CORE_UPCTAUCENTRALBARRELHELPERRL_H_
#define PWGUD_CORE_UPCTAUCENTRALBARRELHELPERRL_H_

#include <string>
#include <algorithm>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

enum MyParticle {
  P_ELECTRON = 0,
  P_MUON = 1,
  P_PION = 2,
  P_KAON = 3,
  P_PROTON = 4
};

void printLargeMessage(std::string info)
// Helper to printf info message to terminal
{
  LOGF(info, "################################### %s ###################################", info);
}

void printMediumMessage(std::string info)
// Helper to printf info message to terminal
{
  LOGF(info, "+++++++++++++ %s +++++++++++++", info);
}

template <typename T>
int testPIDhypothesis(T trackPIDinfo, float nSigmaShift = 0., bool isMC = false)
// Choose, which particle it is according to PID
{
  float nSigmaTPC[5];
  nSigmaTPC[P_ELECTRON] = std::abs(trackPIDinfo.tpcNSigmaEl());
  nSigmaTPC[P_MUON] = std::abs(trackPIDinfo.tpcNSigmaMu());
  nSigmaTPC[P_PION] = std::abs(trackPIDinfo.tpcNSigmaPi());
  nSigmaTPC[P_KAON] = std::abs(trackPIDinfo.tpcNSigmaKa());
  nSigmaTPC[P_PROTON] = std::abs(trackPIDinfo.tpcNSigmaPr());
  // Correction if TPC tuneOnData is wrong
  if (isMC) {
    for (int i = 0; i < 5; i++)
      nSigmaTPC[i] -= nSigmaShift;
  }
  int enumChoiceTPC = std::distance(std::begin(nSigmaTPC),
                                    std::min_element(std::begin(nSigmaTPC), std::end(nSigmaTPC)));

  float nSigmaTOF[5];
  nSigmaTOF[P_ELECTRON] = std::abs(trackPIDinfo.tofNSigmaEl());
  nSigmaTOF[P_MUON] = std::abs(trackPIDinfo.tofNSigmaMu());
  nSigmaTOF[P_PION] = std::abs(trackPIDinfo.tofNSigmaPi());
  nSigmaTOF[P_KAON] = std::abs(trackPIDinfo.tofNSigmaKa());
  nSigmaTOF[P_PROTON] = std::abs(trackPIDinfo.tofNSigmaPr());
  int enumChoiceTOF = std::distance(std::begin(nSigmaTOF),
                                    std::min_element(std::begin(nSigmaTOF), std::end(nSigmaTOF)));

  if (trackPIDinfo.hasTPC() || trackPIDinfo.hasTOF()) {
    if (trackPIDinfo.hasTOF()) {
      return enumChoiceTOF;
    } else {
      return enumChoiceTPC;
    }
  } else {
    LOGF(debug, "testPIDhypothesis failed - track did not leave information in TPC or TOF");
    return -1;
  }
}

template <typename T>
int trackPDG(T trackPIDinfo)
// using testPIDhypothesis, reads enumMyParticle and return pdg value
{
  if (testPIDhypothesis(trackPIDinfo) == P_ELECTRON) {
    return 11;
  } else if (testPIDhypothesis(trackPIDinfo) == P_MUON) {
    return 13;
  } else if (testPIDhypothesis(trackPIDinfo) == P_PION) {
    return 211;
  } else if (testPIDhypothesis(trackPIDinfo) == P_KAON) {
    return 321;
  } else if (testPIDhypothesis(trackPIDinfo) == P_PROTON) {
    return 2212;
  } else {
    printMediumMessage("Something is wrong with track PDG selector");
    return -1.;
  }
}

int enumMyParticle(int valuePDG)
// reads pdg value and returns particle number as in enumMyParticle
{
  if (std::abs(valuePDG) == 11) {
    return P_ELECTRON;
  } else if (std::abs(valuePDG) == 13) {
    return P_MUON;
  } else if (std::abs(valuePDG) == 211) {
    return P_PION;
  } else if (std::abs(valuePDG) == 321) {
    return P_KAON;
  } else if (std::abs(valuePDG) == 2212) {
    return P_PROTON;
  } else {
    printMediumMessage("PDG value not found in enumMyParticle. Returning -1.");
    return -1.;
  }
}

float momentum(float px, float py, float pz)
// Just a simple function to return momentum
{
  return std::sqrt(px * px + py * py + pz * pz);
}

float invariantMass(float E, float px, float py, float pz)
// Just a simple function to return invariant mass
{
  return std::sqrt(E * E - px * px - py * py - pz * pz);
}

float phi(float px, float py)
// Just a simple function to return azimuthal angle
{
  if (px != 0)
    return std::atan(py / px);
  return -999.;
}

float eta(float px, float py, float pz)
// Just a simple function to return pseudorapidity
{
  float eta = -999.;
  float mom = momentum(px, py, pz);
  if (mom != 0)
    eta = std::atanh(pz / mom);
  if (-1. < eta && eta < 1.)
    return eta;
  return -999.;
}

float energy(float mass, float px, float py, float pz)
// Just a simple function to return track energy
{
  return std::sqrt(mass * mass + px * px + py * py + pz * pz);
}

float rapidity(float mass, float px, float py, float pz)
// Just a simple function to return track rapidity
{
  return 0.5 * std::log((energy(mass, px, py, pz) + pz) / (energy(mass, px, py, pz) - pz));
}

double calculateAcoplanarity(double phi_trk1, double phi_trk2)
// Function to calculate acoplanarity of two tracks based on phi of both tracks, which is in interval (0,2*pi)
{
  double aco = std::abs(phi_trk1 - phi_trk2);
  if (aco <= o2::constants::math::PI)
    return aco;
  else
    return (o2::constants::math::TwoPI - aco);
}

template <typename T>
float getAvgITSClSize(T const& track)
{
  float sum = 0.0;
  for (int iL = 0; iL < 6; iL++) {
    sum += (track.itsClusterSizes() >> (iL * 4)) & 0xf;
  }
  return sum / track.itsNCls();
}

template <typename T>
float getCosLambda(T const& track)
{
  // lambda is track inclination
  // tan(lambda) = track.tgl()
  // track.pz() = track.pt() * track.tgl
  float lambda = std::atan(track.pz() / track.pt());
  return std::cos(lambda);
}

template <typename T>
bool passITSAvgClsSizesLowMomCut(T const& track, double itscut = 3.5, double ptcut = 0.7)
{
  if (getAvgITSClSize(track) * getCosLambda(track) < itscut && track.pt() < ptcut) {
    return false;
  } else {
    return true;
  }
}

#endif // PWGUD_CORE_UPCTAUCENTRALBARRELHELPERRL_H_
