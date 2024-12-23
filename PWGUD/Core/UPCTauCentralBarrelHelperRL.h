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

enum MyParticle {
  P_ELECTRON = 0,
  P_MUON = 1,
  P_PION = 2,
  P_KAON = 3,
  P_PROTON = 4,
  P_ENUM_COUNTER = 5
};

enum MyTauChannel {
  CH_EE = 0,
  CH_MUMU = 1,
  CH_EMU = 2,
  CH_PIPI = 3,
  CH_EPI = 4,
  CH_MUPI = 5,
  CH_FOURPI = 6,
  CH_ETHREEPI = 7,
  CH_MUTHREEPI = 8,
  CH_SIXPI = 9,
  CH_EMUPI = 10,
  CH_ENUM_COUNTER = 11
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

void printDebugMessage(std::string info)
// Helper to printf info message to terminal
{
  LOGF(debug, "X!X!X!X!X!X!X!X!X %s X!X!X!X!X!X!X!X!X", info);
}

template <typename T>
int testPIDhypothesis(T trackPIDinfo, float maxNsigmaTPC = 5.0, float maxNsigmaTOF = 5.0, bool useTOF = true, bool useTOFsigmaAfterTPC = true, float nSigmaShift = 0., bool isMC = false)
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

  if (trackPIDinfo.hasTPC()) {
    if (trackPIDinfo.hasTOF() && useTOF) {
      if (nSigmaTOF[enumChoiceTOF] < maxNsigmaTOF) {
        return enumChoiceTOF;
      } else {
        printDebugMessage(Form("testPIDhypothesis cut - the lowest nSigmaTOF is higher than %f", maxNsigmaTPC));
        return -1;
      }
    } else if (trackPIDinfo.hasTOF() && useTOFsigmaAfterTPC) {
      if (nSigmaTPC[enumChoiceTPC] < maxNsigmaTPC && nSigmaTOF[enumChoiceTPC] < maxNsigmaTOF) {
        return enumChoiceTPC;
      } else {
        printDebugMessage(Form("testPIDhypothesis cut - the lowest nSigmaTPC is higher than %f or the lowest nSigmaTOF is higher than %f", maxNsigmaTPC, maxNsigmaTOF));
        return -1;
      }
    } else {
      if (nSigmaTPC[enumChoiceTPC] < maxNsigmaTPC) {
        return enumChoiceTPC;
      } else {
        printDebugMessage(Form("testPIDhypothesis cut - the lowest nSigmaTPC is higher than %f", maxNsigmaTPC));
        return -1;
      }
    }
  } else {
    printDebugMessage("testPIDhypothesis failed - track did not leave information in TPC");
    return -1;
  }
}

template <typename T>
int trackPDG(T trackPIDinfo, float maxNsigmaTPC = 5.0, float maxNsigmaTOF = 5.0, bool useTOF = true, bool useTOFsigmaAfterTPC = true, float nSigmaShift = 0., bool isMC = false)
// using testPIDhypothesis, reads enumMyParticle and return pdg value
{
  if (testPIDhypothesis(trackPIDinfo, maxNsigmaTPC, maxNsigmaTOF, useTOF, useTOFsigmaAfterTPC, nSigmaShift, isMC) == P_ELECTRON) {
    return 11;
  } else if (testPIDhypothesis(trackPIDinfo, maxNsigmaTPC, maxNsigmaTOF, useTOF, useTOFsigmaAfterTPC, nSigmaShift, isMC) == P_MUON) {
    return 13;
  } else if (testPIDhypothesis(trackPIDinfo, maxNsigmaTPC, maxNsigmaTOF, useTOF, useTOFsigmaAfterTPC, nSigmaShift, isMC) == P_PION) {
    return 211;
  } else if (testPIDhypothesis(trackPIDinfo, maxNsigmaTPC, maxNsigmaTOF, useTOF, useTOFsigmaAfterTPC, nSigmaShift, isMC) == P_KAON) {
    return 321;
  } else if (testPIDhypothesis(trackPIDinfo, maxNsigmaTPC, maxNsigmaTOF, useTOF, useTOFsigmaAfterTPC, nSigmaShift, isMC) == P_PROTON) {
    return 2212;
  } else {
    printDebugMessage("Something is wrong with track PDG selector");
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
    printDebugMessage("PDG value not found in enumMyParticle. Returning -1.");
    return -1.;
  }
}

float pt(float px, float py)
// Just a simple function to return pt
{
  return std::sqrt(px * px + py * py);
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
// Just a simple function to return azimuthal angle from 0 to 2pi
{
  if (px != 0)
    return (std::atan2(py, px) + o2::constants::math::PI);
  return -999.;
}

float eta(float px, float py, float pz)
// Just a simple function to return pseudorapidity
{
  float arg = -2.; // outside valid range for std::atanh
  float mom = momentum(px, py, pz);
  if (mom != 0)
    arg = pz / mom;
  if (-1. < arg && arg < 1.)
    return std::atanh(arg); // definition of eta
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

double calculateAcoplanarity(double phiTrk1, double phiTrk2)
// Function to calculate acoplanarity of two tracks based on phi of both tracks, which is in interval (0,2*pi)
{
  double aco = std::abs(phiTrk1 - phiTrk2);
  if (aco <= o2::constants::math::PI)
    return aco;
  else
    return (o2::constants::math::TwoPI - aco);
}

double calculateCollinearity(double etaTrk1, double etaTrk2, double phiTrk1, double phiTrk2)
// Function to calculate deltaR(trk1,trk2) = sqrt(deltaEta^2+deltaPhi^2)
{
  double deltaEta = etaTrk1 - etaTrk2;
  double deltaPhi = phiTrk1 - phiTrk2;
  return std::sqrt(deltaEta * deltaEta + deltaPhi * deltaPhi);
}

template <typename Ps>
int countPhysicalPrimary(Ps particles)
// Function to loop over particles associated to a mcCollision and return total of physical primary particles
{
  int nTotal = 0;
  for (const auto& particle : particles) {
    if (!particle.isPhysicalPrimary())
      continue;
    nTotal++;
  }
  return nTotal;
}

template <typename Ps>
int countParticlesWithoutMother(Ps particles)
// Function to loop over particles associated to a mcCollision and return total of particles without mothers (hopely alternative to isPhysicalPrimary)
{
  int nTotal = 0;
  for (const auto& particle : particles) {
    if (particle.has_mothers())
      continue;
    nTotal++;
  }
  return nTotal;
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
