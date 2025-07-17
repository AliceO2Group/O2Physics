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

#include "FastTracker.h"

#include "ReconstructionDataFormats/TrackParametrization.h"

#include "TMath.h"
#include "TMatrixD.h"
#include "TMatrixDSymEigen.h"
#include "TRandom.h"

#include <string>
#include <vector>

namespace o2
{
namespace fastsim
{

// +-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+

void FastTracker::AddLayer(TString name, float r, float z, float x0, float xrho, float resRPhi, float resZ, float eff, int type)
{
  DetLayer newLayer(name, r, z, x0, xrho, resRPhi, resZ, eff, type);
  // Check that efficient layers are not inert layers
  if (newLayer.getEfficiency() > 0.0f && newLayer.isInert()) {
    LOG(error) << "Layer " << name << " with efficiency > 0.0 should not be inert";
  }
  // Layers should be ordered by increasing radius, check this
  if (!layers.empty() && newLayer.getRadius() < layers.back().getRadius()) {
    LOG(fatal) << "Layer " << newLayer << " is not ordered correctly, it should be after layer " << layers.back();
  }
  // Layers should all have different names
  for (const auto& layer : layers) {
    if (layer.getName() == newLayer.getName()) {
      LOG(fatal) << "Layer with name " << newLayer.getName() << " already exists in FastTracker layers";
    }
  }
  // Add the new layer to the layers vector
  layers.push_back(newLayer);
}

DetLayer FastTracker::GetLayer(int layer, bool ignoreBarrelLayers) const
{
  int layerIdx = layer;
  if (ignoreBarrelLayers) {
    for (int il = 0, trackingLayerIdx = 0; trackingLayerIdx <= layer; il++) {
      if (layers[il].isInert())
        continue;
      trackingLayerIdx++;
      layerIdx = il;
    }
  }
  return layers[layerIdx];
}

int FastTracker::GetLayerIndex(const std::string& name) const
{
  int i = 0;
  for (const auto& layer : layers) {
    if (layer.getName() == name) {
      return i;
    }
    i++;
  }
  LOG(error) << "Layer with name " << name << " not found in FastTracker layers";
  return -1;
}

void FastTracker::Print()
{
  // print out layer setup
  LOG(info) << "+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+";
  LOG(info) << " Printing detector layout with " << layers.size() << " effective elements: ";
  for (uint32_t il = 0; il < layers.size(); il++) {
    LOG(info) << " Layer #" << il << "\t" << layers[il];
  }
  LOG(info) << "+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+";
}

void FastTracker::AddSiliconALICE3v4(std::vector<float> pixelResolution)
{
  LOG(info) << " ALICE 3: Adding v4 tracking layers";
  float x0IT = 0.001;        // 0.1%
  float x0OT = 0.005;        // 0.5%
  float xrhoIB = 1.1646e-02; // 50 mum Si
  float xrhoOT = 1.1646e-01; // 500 mum Si
  float eff = 1.00;

  float resRPhiIT = pixelResolution[0];
  float resZIT = pixelResolution[1];
  float resRPhiOT = pixelResolution[2];
  float resZOT = pixelResolution[3];

  AddLayer("bpipe0", 0.48, 250, 0.00042, 2.772e-02, 0.0f, 0.0f, 0.0f, 0); // 150 mum Be
  AddLayer("ddd0", 0.5, 250, x0IT, xrhoIB, resRPhiIT, resZIT, eff, 1);
  AddLayer("ddd1", 1.2, 250, x0IT, xrhoIB, resRPhiIT, resZIT, eff, 1);
  AddLayer("ddd2", 2.5, 250, x0IT, xrhoIB, resRPhiIT, resZIT, eff, 1);
  AddLayer("bpipe1", 5.7, 250, 0.0014, 9.24e-02, 0.0f, 0.0f, 0.0f, 0); // 500 mum Be
  AddLayer("ddd3", 7., 250, x0OT, xrhoOT, resRPhiOT, resZOT, eff, 1);
  AddLayer("ddd4", 10., 250, x0OT, xrhoOT, resRPhiOT, resZOT, eff, 1);
  AddLayer("ddd5", 13., 250, x0OT, xrhoOT, resRPhiOT, resZOT, eff, 1);
  AddLayer("ddd6", 16., 250, x0OT, xrhoOT, resRPhiOT, resZOT, eff, 1);
  AddLayer("ddd7", 25., 250, x0OT, xrhoOT, resRPhiOT, resZOT, eff, 1);
  AddLayer("ddd8", 40., 250, x0OT, xrhoOT, resRPhiOT, resZOT, eff, 1);
  AddLayer("ddd9", 45., 250, x0OT, xrhoOT, resRPhiOT, resZOT, eff, 1);
}

void FastTracker::AddSiliconALICE3v2(std::vector<float> pixelResolution)
{
  LOG(info) << "ALICE 3: Adding v2 tracking layers;";
  float x0IT = 0.001;        // 0.1%
  float x0OT = 0.01;         // 1.0%
  float xrhoIB = 2.3292e-02; // 100 mum Si
  float xrhoOT = 2.3292e-01; // 1000 mum Si
  float eff = 1.00;

  float resRPhiIT = pixelResolution[0];
  float resZIT = pixelResolution[1];
  float resRPhiOT = pixelResolution[2];
  float resZOT = pixelResolution[3];

  AddLayer("bpipe0", 0.48, 250, 0.00042, 2.772e-02, 0.0f, 0.0f, 0.0f, 0); // 150 mum Be
  AddLayer("B00", 0.5, 250, x0IT, xrhoIB, resRPhiIT, resZIT, eff, 1);
  AddLayer("B01", 1.2, 250, x0IT, xrhoIB, resRPhiIT, resZIT, eff, 1);
  AddLayer("B02", 2.5, 250, x0IT, xrhoIB, resRPhiIT, resZIT, eff, 1);
  AddLayer("bpipe1", 3.7, 250, 0.0014, 9.24e-02, 0.0f, 0.0f, 0.0f, 0); // 500 mum Be
  AddLayer("B03", 3.75, 250, x0OT, xrhoOT, resRPhiOT, resZOT, eff, 1);
  AddLayer("B04", 7., 250, x0OT, xrhoOT, resRPhiOT, resZOT, eff, 1);
  AddLayer("B05", 12., 250, x0OT, xrhoOT, resRPhiOT, resZOT, eff, 1);
  AddLayer("B06", 20., 250, x0OT, xrhoOT, resRPhiOT, resZOT, eff, 1);
  AddLayer("B07", 30., 250, x0OT, xrhoOT, resRPhiOT, resZOT, eff, 1);
  AddLayer("B08", 45., 250, x0OT, xrhoOT, resRPhiOT, resZOT, eff, 1);
  AddLayer("B09", 60., 250, x0OT, xrhoOT, resRPhiOT, resZOT, eff, 1);
  AddLayer("B10", 80., 250, x0OT, xrhoOT, resRPhiOT, resZOT, eff, 1);
}

void FastTracker::AddTPC(float phiResMean, float zResMean)
{
  LOG(info) << " Adding standard time projection chamber";

  // porting of DetectorK::AddTPC
  // see here:
  // https://github.com/AliceO2Group/DelphesO2/blob/master/src/DetectorK/DetectorK.cxx#L522
  // % Radiation Lengths ... Average per TPC row  (i.e. total/159 )
  const int kNPassiveBound = 2;
  const float radLBoundary[kNPassiveBound] = {1.692612e-01, 8.711904e-02};
  const float xrhoBoundary[kNPassiveBound] = {6.795774e+00, 3.111401e+00};
  const float rBoundary[kNPassiveBound] = {50, 70.0}; // cm

  float radLPerRow = 0.000036;

  float tpcInnerRadialPitch = 0.75; // cm
  float tpcMiddleRadialPitch = 1.0; // cm
  float tpcOuterRadialPitch = 1.5;  // cm
  float innerRows = 63;
  float middleRows = 64;
  float outerRows = 32;
  float tpcRows = (innerRows + middleRows + outerRows);
  float rowOneRadius = 85.2;  // cm
  float row64Radius = 135.1;  // cm
  float row128Radius = 199.2; // cm

  float zLength = 250.0f; // to be checked

  // add boundaries between ITS and TPC
  for (int i = 0; i < kNPassiveBound; i++) {
    AddLayer(Form("tpc_boundary%d", i), rBoundary[i], zLength, radLBoundary[i], xrhoBoundary[i], 0); // dummy errors
  }
  for (Int_t k = 0; k < tpcRows; k++) {
    Float_t rowRadius = 0;
    if (k < innerRows)
      rowRadius = rowOneRadius + k * tpcInnerRadialPitch;
    else if (k >= innerRows && k < (innerRows + middleRows))
      rowRadius = row64Radius + (k - innerRows + 1) * tpcMiddleRadialPitch;
    else if (k >= (innerRows + middleRows) && k < tpcRows)
      rowRadius = row128Radius + (k - innerRows - middleRows + 1) * tpcOuterRadialPitch;

    AddLayer(Form("tpc_%d", k), rowRadius, zLength, radLPerRow, 0, phiResMean, zResMean, 1.0f, 2);
  }
}

float FastTracker::Dist(float z, float r)
{
  // porting of DetektorK::Dist
  // see here:
  // https://github.com/AliceO2Group/DelphesO2/blob/master/src/DetectorK/DetectorK.cxx#L743
  int index = 1;
  int nSteps = 301;
  float dist = 0.0;
  float dz0 = (4 * sigmaD - (-4) * sigmaD / (nSteps = 1));
  float z0 = 0.0;
  for (int i = 0; i < nSteps; i++) {
    if (i == nSteps - 1)
      index = 1;
    z0 = -4 * sigmaD + i * dz0;
    dist += index * (dz0 / 3.) * (1 / o2::math_utils::sqrt(o2::constants::math::TwoPI) / sigmaD) * std::exp(-z0 * z0 / 2. / sigmaD / sigmaD) * (1 / o2::math_utils::sqrt((z - z0) * (z - z0) + r * r));
    if (index != 4)
      index = 4;
    else
      index = 2;
  }
  return dist;
}

float FastTracker::OneEventHitDensity(float multiplicity, float radius)
{
  // porting of DetektorK::OneEventHitDensity
  // see here:
  // https://github.com/AliceO2Group/DelphesO2/blob/master/src/DetectorK/DetectorK.cxx#L694
  float den = multiplicity / (o2::constants::math::TwoPI * radius * radius);
  float tg = o2::math_utils::tan(2. * o2::math_utils::atan(-avgRapidity));
  den = den / o2::math_utils::sqrt(1 + 1 / (tg * tg));
  return den;
}

float FastTracker::IntegratedHitDensity(float multiplicity, float radius)
{
  // porting of DetektorK::IntegratedHitDensity
  // see here:
  // https://github.com/AliceO2Group/DelphesO2/blob/master/src/DetectorK/DetectorK.cxx#L712
  float zdcHz = luminosity * 1.e24 * mCrossSectionMinB;
  float den = zdcHz * integrationTime / 1000. * multiplicity * Dist(0., radius) / (o2::constants::math::TwoPI * radius);
  if (den < OneEventHitDensity(multiplicity, radius))
    den = OneEventHitDensity(multiplicity, radius);
  return den;
}

float FastTracker::UpcHitDensity(float radius)
{
  // porting of DetektorK::UpcHitDensity
  // see here:
  // https://github.com/AliceO2Group/DelphesO2/blob/master/src/DetectorK/DetectorK.cxx#L727
  float mUPCelectrons = 0;
  mUPCelectrons = lhcUPCScale * 5456 / (radius * radius) / dNdEtaMinB;
  if (mUPCelectrons < 0)
    mUPCelectrons = 0.0;
  mUPCelectrons *= IntegratedHitDensity(dNdEtaMinB, radius);
  mUPCelectrons *= upcBackgroundMultiplier;
  return mUPCelectrons;
}

float FastTracker::HitDensity(float radius)
{
  // porting of DetektorK::HitDensity
  // see here:
  // https://github.com/AliceO2Group/DelphesO2/blob/master/src/DetectorK/DetectorK.cxx#L663
  float arealDensity = 0.;
  if (radius > maxRadiusSlowDet) {
    arealDensity = OneEventHitDensity(dNdEtaCent, radius);
    arealDensity += otherBackground * OneEventHitDensity(dNdEtaMinB, radius);
  }

  // In the version of Delphes used to produce
  // Look-up tables, UpcHitDensity(radius) always returns 0,
  // hence it is left commented out for now
  if (radius < maxRadiusSlowDet) {
    arealDensity = OneEventHitDensity(dNdEtaCent, radius);
    arealDensity += otherBackground * OneEventHitDensity(dNdEtaMinB, radius) + IntegratedHitDensity(dNdEtaMinB, radius);
    // +UpcHitDensity(radius);
  }
  return arealDensity;
}

float FastTracker::ProbGoodChiSqHit(float radius, float searchRadiusRPhi, float searchRadiusZ)
{
  // porting of DetektorK::ProbGoodChiSqHit
  // see here:
  // https://github.com/AliceO2Group/DelphesO2/blob/master/src/DetectorK/DetectorK.cxx#L629
  float sx, goodHit;
  sx = o2::constants::math::TwoPI * searchRadiusRPhi * searchRadiusZ * HitDensity(radius);
  goodHit = 1. / (1 + sx);
  return goodHit;
}

// function to provide a reconstructed track from a perfect input track
// returns number of intercepts (generic for now)
int FastTracker::FastTrack(o2::track::TrackParCov inputTrack, o2::track::TrackParCov& outputTrack, const float nch)
{
  dNdEtaCent = nch; // set the number of charged particles per unit rapidity
  hits.clear();
  nIntercepts = 0;
  nSiliconPoints = 0;
  nGasPoints = 0;
  std::array<float, 3> posIni; // provision for != PV
  inputTrack.getXYZGlo(posIni);
  const float initialRadius = std::hypot(posIni[0], posIni[1]);
  const float kTrackingMargin = 0.1;
  const int kMaxNumberOfDetectors = 20;
  if (kMaxNumberOfDetectors < layers.size()) {
    LOG(fatal) << "Too many layers in FastTracker, increase kMaxNumberOfDetectors";
    return -1; // too many layers
  }
  int firstActiveLayer = -1; // first layer that is not inert
  for (size_t i = 0; i < layers.size(); ++i) {
    if (!layers[i].isInert()) {
      firstActiveLayer = i;
      break;
    }
  }
  if (firstActiveLayer <= 0) {
    LOG(fatal) << "No active layers found in FastTracker, check layer setup";
    return -2; // no active layers
  }
  const int xrhosteps = 100;
  const bool applyAngularCorrection = true;

  goodHitProbability.clear();
  for (int i = 0; i < kMaxNumberOfDetectors; ++i) {
    goodHitProbability.push_back(-1.);
  }
  goodHitProbability[0] = 1.; // we use layer zero to accumulate

  // +-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+
  // Outward pass to find intercepts
  int firstLayerReached = -1;
  int lastLayerReached = -1;
  new (&outputTrack)(o2::track::TrackParCov)(inputTrack);
  for (size_t il = 0; il < layers.size(); il++) {
    // check if layer is doable
    if (layers[il].getRadius() < initialRadius) {
      continue; // this layer should not be attempted, but go ahead
    }

    // check if layer is reached
    float targetX = 1e+3;
    inputTrack.getXatLabR(layers[il].getRadius(), targetX, magneticField);
    if (targetX > 999.f) {
      LOGF(debug, "Failed to find intercept for layer %d at radius %.2f cm", il, layers[il].getRadius());
      break; // failed to find intercept
    }

    bool ok = inputTrack.propagateTo(targetX, magneticField);
    if (ok && mApplyMSCorrection && layers[il].getRadiationLength() > 0) {
      ok = inputTrack.correctForMaterial(layers[il].getRadiationLength(), 0, applyAngularCorrection);
    }
    if (ok && mApplyElossCorrection && layers[il].getDensity() > 0) { // correct in small steps
      for (int ise = xrhosteps; ise--;) {
        ok = inputTrack.correctForMaterial(0, -layers[il].getDensity() / xrhosteps, applyAngularCorrection);
        if (!ok)
          break;
      }
    }
    LOGF(debug, "Propagation was %s up to layer %d", ok ? "successful" : "unsuccessful", il);

    // was there a problem on this layer?
    if (!ok && il > 0) { // may fail to reach target layer due to the eloss
      float rad2 = inputTrack.getX() * inputTrack.getX() + inputTrack.getY() * inputTrack.getY();
      float maxR = layers[il - 1].getRadius() + kTrackingMargin * 2;
      float minRad = (fMinRadTrack > 0 && fMinRadTrack < maxR) ? fMinRadTrack : maxR;
      if (rad2 - minRad * minRad < kTrackingMargin * kTrackingMargin) { // check previously reached layer
        return -5;                                                      // did not reach min requested layer
      } else {
        break;
      }
    }
    if (std::abs(inputTrack.getZ()) > layers[il].getZ() && mApplyZacceptance) {
      break; // out of acceptance bounds
    }

    if (layers[il].isInert()) {
      LOG(info) << "Skipping inert layer: " << layers[il].getName() << " at radius " << layers[il].getRadius() << " cm";
      continue; // inert layer, skip
    }

    // layer is reached
    if (firstLayerReached < 0) {
      LOGF(debug, "First layer reached: %d", il);
      firstLayerReached = il;
    }
    lastLayerReached = il;
    nIntercepts++;
  }

  // +-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+
  // initialize track at outer point
  o2::track::TrackParCov inwardTrack(inputTrack);

  // Enlarge covariance matrix
  std::array<float, o2::track::kNParams> trPars = {0.};
  for (int ip = 0; ip < o2::track::kNParams; ip++) {
    trPars[ip] = outputTrack.getParam(ip);
  }
  static constexpr float kLargeErr2Coord = 5 * 5;
  static constexpr float kLargeErr2Dir = 0.7 * 0.7;
  static constexpr float kLargeErr2PtI = 30.5 * 30.5;
  std::array<float, o2::track::kCovMatSize> largeCov = {0.};
  for (int ic = o2::track::kCovMatSize; ic--;)
    largeCov[ic] = 0.;
  largeCov[o2::track::CovLabels::kSigY2] = largeCov[o2::track::CovLabels::kSigZ2] = kLargeErr2Coord;
  largeCov[o2::track::CovLabels::kSigSnp2] = largeCov[o2::track::CovLabels::kSigTgl2] = kLargeErr2Dir;
  largeCov[o2::track::CovLabels::kSigQ2Pt2] = kLargeErr2PtI * trPars[o2::track::ParLabels::kQ2Pt] * trPars[o2::track::ParLabels::kQ2Pt];

  inwardTrack.setCov(largeCov);
  inwardTrack.checkCovariance();

  // +-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+
  // Inward pass to calculate covariances
  for (int il = lastLayerReached; il >= firstLayerReached; il--) {

    float targetX = 1e+3;
    inputTrack.getXatLabR(layers[il].getRadius(), targetX, magneticField);
    if (targetX > 999)
      continue; // failed to find intercept

    if (!inputTrack.propagateTo(targetX, magneticField)) {
      continue; // failed to propagate
    }

    if (std::abs(inputTrack.getZ()) > layers[il].getZ() && mApplyZacceptance) {
      continue; // out of acceptance bounds but continue inwards
    }

    // get perfect data point position
    std::array<float, 3> spacePoint;
    inputTrack.getXYZGlo(spacePoint);
    std::vector<float> thisHit = {spacePoint[0], spacePoint[1], spacePoint[2]};

    // towards adding cluster: move to track alpha
    float alpha = inwardTrack.getAlpha();
    float xyz1[3]{
      std::cos(alpha) * spacePoint[0] + std::sin(alpha) * spacePoint[1],
      -std::sin(alpha) * spacePoint[0] + std::cos(alpha) * spacePoint[1],
      spacePoint[2]};
    if (!inwardTrack.propagateTo(xyz1[0], magneticField))
      continue;

    if (!layers[il].isInert()) { // only update covm for tracker hits
      const o2::track::TrackParametrization<float>::dim2_t hitpoint = {
        static_cast<float>(xyz1[1]),
        static_cast<float>(xyz1[2])};
      const o2::track::TrackParametrization<float>::dim3_t hitpointcov = {layers[il].getResolutionRPhi() * layers[il].getResolutionRPhi(), 0.f, layers[il].getResolutionZ() * layers[il].getResolutionZ()};

      inwardTrack.update(hitpoint, hitpointcov);
      inwardTrack.checkCovariance();
    }

    if (mApplyMSCorrection && layers[il].getRadiationLength() > 0) {
      if (!inputTrack.correctForMaterial(layers[il].getRadiationLength(), 0, applyAngularCorrection)) {
        return -6;
      }
      if (!inwardTrack.correctForMaterial(layers[il].getRadiationLength(), 0, applyAngularCorrection)) {
        return -6;
      }
    }
    if (mApplyElossCorrection && layers[il].getDensity() > 0) {
      for (int ise = xrhosteps; ise--;) { // correct in small steps
        if (!inputTrack.correctForMaterial(0, layers[il].getDensity() / xrhosteps, applyAngularCorrection)) {
          return -7;
        }
        if (!inwardTrack.correctForMaterial(0, layers[il].getDensity() / xrhosteps, applyAngularCorrection)) {
          return -7;
        }
      }
    }

    if (layers[il].isSilicon())
      nSiliconPoints++; // count silicon hits
    if (layers[il].isGas())
      nGasPoints++; // count TPC/gas hits

    hits.push_back(thisHit);

    if (!layers[il].isInert()) { // good hit probability calculation
      float sigYCmb = o2::math_utils::sqrt(inwardTrack.getSigmaY2() + layers[il].getResolutionRPhi() * layers[il].getResolutionRPhi());
      float sigZCmb = o2::math_utils::sqrt(inwardTrack.getSigmaZ2() + layers[il].getResolutionZ() * layers[il].getResolutionZ());
      goodHitProbability[il] = ProbGoodChiSqHit(layers[il].getRadius() * 100, sigYCmb * 100, sigZCmb * 100);
      goodHitProbability[0] *= goodHitProbability[il];
    }
  }

  // backpropagate to original radius
  float finalX = 1e+3;
  bool inPropStatus = inwardTrack.getXatLabR(initialRadius, finalX, magneticField);
  if (finalX > 999) {
    LOG(debug) << "Failed to find intercept for initial radius " << initialRadius << " cm, x = " << finalX << " and status " << inPropStatus << " and sn = " << inwardTrack.getSnp() << " r = " << inwardTrack.getY() * inwardTrack.getY();
    return -3; // failed to find intercept
  }

  if (!inwardTrack.propagateTo(finalX, magneticField)) {
    return -4; // failed to propagate
  }

  // only attempt to continue if intercepts are at least four
  if (nIntercepts < 4)
    return nIntercepts;

  // generate efficiency
  float eff = 1.;
  for (int i = 0; i < kMaxNumberOfDetectors; i++) {
    float iGoodHit = goodHitProbability[i];
    if (iGoodHit <= 0)
      continue;

    eff *= iGoodHit;
  }
  if (mApplyEffCorrection) {
    if (gRandom->Uniform() > eff)
      return -8;
  }

  outputTrack.setCov(inwardTrack.getCov());
  outputTrack.checkCovariance();

  // Use covariance matrix based smearing
  std::array<float, o2::track::kCovMatSize> covMat = {0.};
  for (int ii = 0; ii < o2::track::kCovMatSize; ii++)
    covMat[ii] = outputTrack.getCov()[ii];
  TMatrixDSym m(5);
  double fcovm[5][5]; // double precision is needed for regularisation

  for (int ii = 0, k = 0; ii < 5; ++ii) {
    for (int j = 0; j < ii + 1; ++j, ++k) {
      fcovm[ii][j] = covMat[k];
      fcovm[j][ii] = covMat[k];
    }
  }

  // evaluate ruben's conditional, regularise
  const bool makePositiveDefinite = (covMatFactor > -1e-5); // apply fix
  bool rubenConditional = false;
  for (int ii = 0; ii < 5; ii++) {
    for (int jj = 0; jj < 5; jj++) {
      if (ii == jj)
        continue; // don't evaluate diagonals
      if (fcovm[ii][jj] * fcovm[ii][jj] > std::abs(fcovm[ii][ii] * fcovm[jj][jj])) {
        rubenConditional = true;
        if (makePositiveDefinite) {
          fcovm[ii][jj] = TMath::Sign(1, fcovm[ii][jj]) * covMatFactor * sqrt(std::abs(fcovm[ii][ii] * fcovm[jj][jj]));
        }
      }
    }
  }

  // Should have a valid cov matrix now
  m.SetMatrixArray(reinterpret_cast<double*>(fcovm));
  TMatrixDSymEigen eigen(m);
  TMatrixD eigVec = eigen.GetEigenVectors();
  TVectorD eigVal = eigen.GetEigenValues();
  bool negEigVal = false;
  for (int ii = 0; ii < 5; ii++) {
    if (eigVal[ii] < 0.0f)
      negEigVal = true;
  }

  if (negEigVal && rubenConditional && makePositiveDefinite) {
    if (mVerboseLevel > 0) {
      LOG(info) << "WARNING: this diagonalization (at pt = " << inputTrack.getPt() << ") has negative eigenvalues despite Ruben's fix! Please be careful!";
      LOG(info) << "Printing info:";
      LOG(info) << "Kalman updates: " << nIntercepts;
      LOG(info) << "Cov matrix: ";
      m.Print();
    }
    covMatNotOK++;
    nIntercepts = -1; // mark as problematic so that it isn't used
    return -1;
  }
  covMatOK++;

  // transform parameter vector and smear
  float params_[5];
  for (int ii = 0; ii < 5; ++ii) {
    float val = 0.;
    for (int j = 0; j < 5; ++j)
      val += eigVec[j][ii] * outputTrack.getParam(j);
    // smear parameters according to eigenvalues
    params_[ii] = gRandom->Gaus(val, sqrt(eigVal[ii]));
  }

  // invert eigenvector matrix
  eigVec.Invert();
  // transform back params vector
  for (int ii = 0; ii < 5; ++ii) {
    float val = 0.;
    for (int j = 0; j < 5; ++j)
      val += eigVec[j][ii] * params_[j];
    outputTrack.setParam(val, ii);
  }
  // should make a sanity check that par[2] sin(phi) is in [-1, 1]
  if (fabs(outputTrack.getParam(2)) > 1.) {
    LOG(info) << " --- smearTrack failed sin(phi) sanity check: " << outputTrack.getParam(2);
    return -2;
  }

  return nIntercepts;
}
// +-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+

} /* namespace fastsim */
} /* namespace o2 */

ClassImp(o2::fastsim::FastTracker);
