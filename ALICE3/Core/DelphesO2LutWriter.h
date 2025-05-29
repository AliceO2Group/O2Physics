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
/// @file DelphesO2LutWriter.h
/// @brief Porting to O2Physics of DelphesO2 code.
///        Minimal changes have been made to the original code for adaptation purposes, formatting and commented parts have been considered.
///        Relevant sources:
///                 DelphesO2/src/lutWrite.cc https://github.com/AliceO2Group/DelphesO2/blob/master/src/lutWrite.cc
/// @author: Roberto Preghenella
/// @email: preghenella@bo.infn.it
///

#ifndef ALICE3_CORE_DELPHESO2LUTWRITER_H_
#define ALICE3_CORE_DELPHESO2LUTWRITER_H_

#include "ALICE3/Core/DelphesO2TrackSmearer.h"
#include "ALICE3/Core/FastTracker.h"
#include "TGraph.h"

namespace o2::fastsim
{
class DelphesO2LutWriter
{
 public:
  DelphesO2LutWriter() = default;
  virtual ~DelphesO2LutWriter() = default;

  o2::fastsim::FastTracker fat;
  void diagonalise(lutEntry_t& lutEntry);
  float etaMaxBarrel = 1.75f;
  bool usePara = true;        // use fwd parameterisation
  bool useDipole = false;     // use dipole i.e. flat parametrization for efficiency and momentum resolution
  bool useFlatDipole = false; // use dipole i.e. flat parametrization outside of the barrel

  int mAtLeastHits = 4;
  int mAtLeastCorr = 4;
  int mAtLeastFake = 0;
  void SetAtLeastHits(int n) { mAtLeastHits = n; }
  void SetAtLeastCorr(int n) { mAtLeastCorr = n; }
  void SetAtLeastFake(int n) { mAtLeastFake = n; }

  void printLutWriterConfiguration();
  bool fatSolve(lutEntry_t& lutEntry,
                float pt = 0.1,
                float eta = 0.0,
                const float mass = 0.13957000,
                int itof = 0,
                int otof = 0,
                int q = 1,
                const float nch = 1);
  bool fwdSolve(float* covm, float pt = 0.1, float eta = 0.0, float mass = 0.13957000);
  bool fwdPara(lutEntry_t& lutEntry, float pt = 0.1, float eta = 0.0, float mass = 0.13957000, float Bfield = 0.5);
  void lutWrite(const char* filename = "lutCovm.dat", int pdg = 211, float field = 0.2, int itof = 0, int otof = 0);
  TGraph* lutRead(const char* filename, int pdg, int what, int vs, float nch = 0., float radius = 0., float eta = 0., float pt = 0.);

  ClassDef(DelphesO2LutWriter, 1);
};
} // namespace o2::fastsim

#endif // ALICE3_CORE_DELPHESO2LUTWRITER_H_
