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

#ifndef ALICE3_CORE_FASTTRACKER_H_
#define ALICE3_CORE_FASTTRACKER_H_

#include <fairlogger/Logger.h> // not a system header but megalinter thinks so
#include <vector>
#include "DetLayer.h"
#include "ReconstructionDataFormats/Track.h"

namespace o2
{
namespace fastsim
{

// +-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+

// this class implements a synthetic smearer that allows
// for on-demand smearing of TrackParCovs in a certain flexible t
// detector layout.
class FastTracker
{
 public:
  // Constructor/destructor
  FastTracker();
  virtual ~FastTracker() {}

  void AddLayer(TString name, float r, float z, float x0, float xrho, float resRPhi = 0.0f, float resZ = 0.0f, float eff = 0.0f, int type = 0);
  DetLayer GetLayer(const int layer, bool ignoreBarrelLayers = true);

  void AddSiliconALICE3v4(std::vector<float> pixelResolution);
  void AddSiliconALICE3v2(std::vector<float> pixelResolution);
  void AddTPC(float phiResMean, float zResMean);

  void Print();
  int FastTrack(o2::track::TrackParCov inputTrack, o2::track::TrackParCov& outputTrack, float nch);

  // For efficiency calculation
  float Dist(float z, float radius);
  float OneEventHitDensity(float multiplicity, float radius);
  float IntegratedHitDensity(float multiplicity, float radius);
  float UpcHitDensity(float radius);
  float HitDensity(float radius);
  float ProbGoodChiSqHit(float radius, float searchRadiusRPhi, float searchRadiusZ);

  // Definition of detector layers
  std::vector<DetLayer> layers;
  std::vector<std::vector<float>> hits; // bookkeep last added hits

  // operational
  bool applyZacceptance;     // check z acceptance or not
  bool applyMSCorrection;    // Apply correction for multiple scattering
  bool applyElossCorrection; // Apply correction for eloss (requires MS correction)
  bool applyEffCorrection;   // Apply correction for hit efficiency
  int verboseLevel;          // 0: not verbose, >0 more verbose
  int crossSectionMinB;
  int dNdEtaCent;
  int dNdEtaMinB;
  float integrationTime;
  float magneticField; // in kiloGauss (5 = 0.5T, etc)
  float covMatFactor;  // covmat off-diagonal factor to use for covmat fix (negative: no factor)
  float sigmaD;
  float luminosity;
  float otherBackground;
  float maxRadiusSlowDet;
  float avgRapidity;
  float lhcUPCScale;
  float upcBackgroundMultiplier;

  uint64_t covMatOK;    // cov mat has negative eigenvals
  uint64_t covMatNotOK; // cov mat has negative eigenvals

  // last track information
  int nIntercepts;    // found in first outward propagation
  int nSiliconPoints; // silicon-based space points added to track
  int nGasPoints;     // tpc-based space points added to track
  std::vector<float> goodHitProbability;

  ClassDef(FastTracker, 1);
};

// +-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+

} // namespace fastsim
} // namespace o2

#endif // ALICE3_CORE_FASTTRACKER_H_
