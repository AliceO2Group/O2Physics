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

#ifndef ALICE3_CORE_FLATLUTWRITER_H_
#define ALICE3_CORE_FLATLUTWRITER_H_

#include "ALICE3/Core/FastTracker.h"
#include "ALICE3/Core/FlatLutEntry.h"

#include <ReconstructionDataFormats/PID.h>

#include <TGraph.h>

#include <cstddef>
#include <string>
namespace o2::fastsim
{
using lutEntry_t = o2::delphes::lutEntry_t;
/**
 * @brief LUT writer using flat binary format
 *
 * Generates lookup tables in the optimized flat format compatible with FlatLutData.
 * The output format is:
 *   [lutHeader_t][lutEntry_t_0][lutEntry_t_1]...[lutEntry_t_N]
 *
 * This enables zero-copy loading into the TrackSmearer and direct shared memory mapping.
 */
class FlatLutWriter
{
 public:
  FlatLutWriter() = default;

  // Setters for binning configuration
  void setBinningNch(bool log, int nbins, float min, float max) { mNchBinning = {log, nbins, min, max}; }
  void setBinningRadius(bool log, int nbins, float min, float max) { mRadiusBinning = {log, nbins, min, max}; }
  void setBinningEta(bool log, int nbins, float min, float max) { mEtaBinning = {log, nbins, min, max}; }
  void setBinningPt(bool log, int nbins, float min, float max) { mPtBinning = {log, nbins, min, max}; }

  void setEtaMaxBarrel(float eta) { etaMaxBarrel = eta; }
  void setAtLeastHits(int n) { mAtLeastHits = n; }
  void setAtLeastCorr(int n) { mAtLeastCorr = n; }
  void setAtLeastFake(int n) { mAtLeastFake = n; }

  bool fatSolve(lutEntry_t& lutEntry,
                float pt = 0.1f,
                float eta = 0.0f,
                const float mass = o2::track::pid_constants::sMasses[o2::track::PID::Pion],
                size_t itof = 0,
                size_t otof = 0,
                int q = 1,
                const float nch = 1.0f);

  void print() const;
  bool fwdSolve(float* covm, float pt = 0.1f, float eta = 0.0f, float mass = o2::track::pid_constants::sMasses[o2::track::PID::Pion]);
  bool fwdPara(lutEntry_t& lutEntry, float pt = 0.1f, float eta = 0.0f, float mass = o2::track::pid_constants::sMasses[o2::track::PID::Pion], float Bfield = 0.5f);
  void lutWrite(const char* filename = "lutCovm.dat", int pdg = 211, float field = 0.2f, size_t itof = 0, size_t otof = 0);
  TGraph* lutRead(const char* filename, int pdg, int what, int vs, float nch = 0.f, float radius = 0.f, float eta = 0.f, float pt = 0.f);

  o2::fastsim::FastTracker fat;

 private:
  void diagonalise(lutEntry_t& lutEntry);

  float etaMaxBarrel = 1.75f;
  bool usePara = true;        // use fwd parametrisation
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
};

} // namespace o2::fastsim

#endif // ALICE3_CORE_FLATLUTWRITER_H_
