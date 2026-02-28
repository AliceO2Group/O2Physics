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
/// \file DelphesO2LutWriter.h
/// \brief Porting to O2Physics of DelphesO2 code.
///        Minimal changes have been made to the original code for adaptation purposes, formatting and commented parts have been considered.
///        Relevant sources:
///                 DelphesO2/src/lutWrite.cc https://github.com/AliceO2Group/DelphesO2/blob/master/src/lutWrite.cc
/// \author: Roberto Preghenella
/// \email: preghenella@bo.infn.it
///

#ifndef ALICE3_CORE_DELPHESO2LUTWRITER_H_
#define ALICE3_CORE_DELPHESO2LUTWRITER_H_

#include "ALICE3/Core/DelphesO2TrackSmearer.h"
#include "ALICE3/Core/FastTracker.h"

#include "ReconstructionDataFormats/PID.h"

#include "TGraph.h"

#include <string>

namespace o2::fastsim
{
class DelphesO2LutWriter
{
 public:
  DelphesO2LutWriter() = default;
  virtual ~DelphesO2LutWriter() = default;

  // Setters
  void setBinningNch(bool log, int nbins, float min, float max) { mNchBinning = {log, nbins, min, max}; }
  void setBinningRadius(bool log, int nbins, float min, float max) { mRadiusBinning = {log, nbins, min, max}; }
  void setBinningEta(bool log, int nbins, float min, float max) { mEtaBinning = {log, nbins, min, max}; }
  void setBinningPt(bool log, int nbins, float min, float max) { mPtBinning = {log, nbins, min, max}; }
  void setEtaMaxBarrel(float eta) { etaMaxBarrel = eta; }
  void setAtLeastHits(int n) { mAtLeastHits = n; }
  void setAtLeastCorr(int n) { mAtLeastCorr = n; }
  void setAtLeastFake(int n) { mAtLeastFake = n; }
  bool fatSolve(o2::delphes::DelphesO2TrackSmearer::lutEntry_t& lutEntry,
                float pt = 0.1,
                float eta = 0.0,
                const float mass = o2::track::pid_constants::sMasses[o2::track::PID::Pion],
                size_t itof = 0,
                size_t otof = 0,
                int q = 1,
                const float nch = 1);

  void print() const;
  bool fwdSolve(float* covm, float pt = 0.1, float eta = 0.0, float mass = o2::track::pid_constants::sMasses[o2::track::PID::Pion]);
  bool fwdPara(o2::delphes::DelphesO2TrackSmearer::lutEntry_t& lutEntry, float pt = 0.1, float eta = 0.0, float mass = o2::track::pid_constants::sMasses[o2::track::PID::Pion], float Bfield = 0.5);
  void lutWrite(const char* filename = "lutCovm.dat", int pdg = 211, float field = 0.2, size_t itof = 0, size_t otof = 0);
  TGraph* lutRead(const char* filename, int pdg, int what, int vs, float nch = 0., float radius = 0., float eta = 0., float pt = 0.);

  o2::fastsim::FastTracker fat;

 private:
  void diagonalise(o2::delphes::DelphesO2TrackSmearer::lutEntry_t& lutEntry);
  float etaMaxBarrel = 1.75f;
  bool usePara = true;        // use fwd parameterisation
  bool useDipole = false;     // use dipole i.e. flat parametrization for efficiency and momentum resolution
  bool useFlatDipole = false; // use dipole i.e. flat parametrization outside of the barrel

  int mAtLeastHits = 4;
  int mAtLeastCorr = 4;
  int mAtLeastFake = 0;

  // Binning of the LUT to make
  struct LutBinning {
    bool log;
    int nbins;
    float min;
    float max;
    std::string toString() const;
  };
  LutBinning mNchBinning = {true, 20, 0.5f, 3.5f};
  LutBinning mRadiusBinning = {false, 1, 0.0f, 100.0f};
  LutBinning mEtaBinning = {false, 80, -4.0f, 4.0f};
  LutBinning mPtBinning = {true, 200, -2.0f, 2.0f};

  ClassDef(DelphesO2LutWriter, 1);
};
} // namespace o2::fastsim

#endif // ALICE3_CORE_DELPHESO2LUTWRITER_H_
