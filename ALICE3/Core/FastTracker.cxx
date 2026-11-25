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

/// \file FastTracker.cxx
/// \brief On the fly implementation of DelphesO2 solveTrack
/// \author David Dobrigkeit Chinellato
/// \author Nicolò Jacazio
/// \author Jesper Karlsson Gumprecht

#include "FastTracker.h"

#include "DetLayer.h"
#include "GeometryContainer.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <Framework/Logger.h>
#include <MathUtils/Utils.h>
#include <ReconstructionDataFormats/Track.h>
#include <ReconstructionDataFormats/TrackParametrization.h>

#include <TFile.h>
#include <TGraph.h>
#include <TMathBase.h>
#include <TMatrixD.h> // IWYU pragma: keep (do not replace with TMatrixDfwd.h)
#include <TMatrixDSymEigen.h>
#include <TMatrixDSymfwd.h>
#include <TMatrixDfwd.h>
#include <TObject.h>
#include <TRandom.h>
#include <TString.h>
#include <TVectorDfwd.h>

#include <Rtypes.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <string>
#include <vector>

namespace o2::fastsim
{

FastTracker::SiliconLayerHits FastTracker::countSiliconLayerHits(const o2::track::TrackParCov& track) const
{
  SiliconLayerHits result;
  std::array<float, 3> position{};
  track.getXYZGlo(position);
  const double x = position[0], y = position[1], z = position[2];
  const double phi = track.getPhi(), tgl = track.getTgl(), curvature = track.getCurvature(magneticField);
  if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z) ||
      !std::isfinite(phi) || !std::isfinite(tgl) || !std::isfinite(curvature)) {
    LOG(error) << "Cannot count silicon layer hits for non-finite track parameters";
    result.count = -1;
    return result;
  }
  const double ux = std::cos(phi), uy = std::sin(phi);
  constexpr double MinCurvature = 1.e-8;
  constexpr double TwoPi = 6.2831853071795864769;

  for (const auto& layer : layers) {
    if (!layer.isSilicon() && !layer.isSiliconDisk()) {
      continue;
    }
    if (layer.isSiliconDisk()) {
      const double diskZ = layer.getDiskZ();
      const double rMin = layer.getDiskInnerRadius();
      const double rMax = layer.getDiskOuterRadius();
      if (!std::isfinite(diskZ) || !std::isfinite(rMin) || !std::isfinite(rMax) ||
          rMin < 0. || rMax < rMin) {
        LOG(fatal) << "Invalid silicon disk geometry: " << layer;
        continue;
      }
      if (std::abs(tgl) < 1.e-12) {
        continue;
      }
      const double distance = (diskZ - z) / tgl;
      if (distance < 0.) {
        continue;
      }
      if (std::abs(curvature) >= MinCurvature) {
        if (distance > TwoPi / std::abs(curvature)) {
          continue;
        }
      }

      const double cx = std::abs(curvature) < MinCurvature ? x : x - uy / curvature;
      const double cy = std::abs(curvature) < MinCurvature ? y : y + ux / curvature;
      const double startTheta = std::atan2(y - cy, x - cx);
      const double diskX = std::abs(curvature) < MinCurvature ? x + ux * distance : cx + std::cos(startTheta + curvature * distance) / std::abs(curvature);
      const double diskY = std::abs(curvature) < MinCurvature ? y + uy * distance : cy + std::sin(startTheta + curvature * distance) / std::abs(curvature);
      const double hitRadius = std::hypot(diskX, diskY);
      if (hitRadius >= rMin && hitRadius <= rMax) {
        ++result.count;
        const float rHit = static_cast<float>(hitRadius);
        result.innerRadius = result.innerRadius < 0.f ? rHit : std::min(result.innerRadius, rHit);
        result.outerRadius = std::max(result.outerRadius, rHit);
      }
      continue;
    }
    const double r = layer.getRadius();
    const double halfLength = layer.getZ();
    if (!std::isfinite(r) || !std::isfinite(halfLength) || r <= 0. || halfLength < 0.) {
      LOG(fatal) << "Invalid silicon layer geometry: " << layer;
      continue;
    }
    bool reached = false;
    if (std::abs(curvature) < MinCurvature) {
      const double projection = x * ux + y * uy;
      const double discriminant = projection * projection + r * r - x * x - y * y;
      if (discriminant >= 0.) {
        const double root = std::sqrt(discriminant);
        for (const double distance : {-projection - root, -projection + root}) {
          if (distance >= 0. && std::abs(z + distance * tgl) <= halfLength) {
            reached = true;
            break;
          }
        }
      }
    } else {
      const double cx = x - uy / curvature, cy = y + ux / curvature;
      const double trackRadius = 1. / std::abs(curvature);
      const double centerDistance = std::hypot(cx, cy);
      if (centerDistance < 1.e-10 && std::abs(trackRadius - r) < 1.e-5) {
        const double endZ = z + TwoPi * trackRadius * tgl;
        reached = std::min(z, endZ) <= halfLength && std::max(z, endZ) >= -halfLength;
      } else if (centerDistance > 0. && centerDistance <= trackRadius + r &&
                 centerDistance >= std::abs(trackRadius - r)) {
        const double along = (r * r + centerDistance * centerDistance - trackRadius * trackRadius) / (2. * centerDistance);
        const double perpendicular = std::sqrt(std::max(0., r * r - along * along));
        const double nx = cx / centerDistance, ny = cy / centerDistance;
        const double startX = x - cx, startY = y - cy;
        for (const double sign : {-1., 1.}) {
          const double ix = along * nx - sign * perpendicular * ny;
          const double iy = along * ny + sign * perpendicular * nx;
          const double vx = ix - cx, vy = iy - cy;
          double angle = std::atan2(startX * vy - startY * vx, startX * vx + startY * vy);
          if (std::hypot(ix - x, iy - y) < 1.e-5) {
            angle = 0.;
          }
          if (curvature < 0.) {
            angle = -angle;
          }
          if (angle < 0.) {
            angle += TwoPi;
          }
          const double distance = angle / std::abs(curvature);
          if (std::abs(z + distance * tgl) <= halfLength) {
            reached = true;
            break;
          }
        }
      }
    }
    if (reached) {
      ++result.count;
      if (result.innerRadius < 0.f || r < result.innerRadius) {
        result.innerRadius = r;
      }
      result.outerRadius = std::max(result.outerRadius, static_cast<float>(r));
    }
  }
  return result;
}

// +-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+

DetLayer* FastTracker::addLayer(const TString& name, float r, float z, float x0, float xrho, float resRPhi, float resZ, float eff, int type)
{
  LOG(debug) << "Adding layer " << name << " r=" << r << " z=" << z << " x0=" << x0 << " xrho=" << xrho << " resRPhi=" << resRPhi << " resZ=" << resZ << " eff=" << eff << " type=" << type;
  DetLayer newLayer(name, r, z, x0, xrho, resRPhi, resZ, eff, type);
  if (newLayer.isSiliconDisk()) {
    newLayer.setSiliconDiskAcceptance(r, r, z);
  }
  // Check that efficient layers are not inert layers
  if (newLayer.getEfficiency() > 0.0f && newLayer.isInert()) {
    LOG(error) << "Layer " << name << " with efficiency > 0.0 should not be inert";
  }
  // Layers should be ordered by increasing radius, check this
  if (!newLayer.isSiliconDisk()) {
    const auto previous = std::find_if(layers.rbegin(), layers.rend(),
                                       [](const auto& layer) {
                                         return !layer.isSiliconDisk();
                                       });

    if (previous != layers.rend() && newLayer.getRadius() < previous->getRadius()) {
      LOG(fatal) << "Layer " << newLayer << " is not ordered correctly, it should be after layer " << previous->getName();
    }
  }
  // Layers should all have different names
  for (const auto& layer : layers) {
    if (layer.getName() == newLayer.getName()) {
      LOG(fatal) << "Layer with name " << newLayer.getName() << " already exists in FastTracker layers";
    }
  }
  // Add the new layer to the layers vector
  layers.push_back(newLayer);
  // Return the last added layer
  return &layers.back();
}

void FastTracker::addDeadPhiRegionInLayer(const std::string& layerName, float phiStart, float phiEnd)
{
  const int layerIdx = getLayerIndex(layerName);
  if (layerIdx < 0) {
    LOG(fatal) << "Cannot add dead phi region to non-existing layer " << layerName;
    return;
  }
  layers[layerIdx].addDeadPhiRegion(phiStart, phiEnd);
}

int FastTracker::getLayerIndex(const std::string& name) const
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

void FastTracker::print()
{
  // print out layer setup
  LOG(info) << "+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+";
  LOG(info) << " Printing detector layout with " << layers.size() << " effective elements: ";
  for (uint32_t il = 0; il < layers.size(); il++) {
    LOG(info) << " Layer #" << il << "\t" << layers[il];
  }
  LOG(info) << "+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+";
}

void FastTracker::addTPC(float phiResMean, float zResMean)
{
  LOG(info) << " Adding standard time projection chamber";

  // porting of DetectorK::AddTPC
  // see here:
  // https://github.com/AliceO2Group/DelphesO2/blob/master/src/DetectorK/DetectorK.cxx#L522
  // % Radiation Lengths ... Average per TPC row  (i.e. total/159 )
  const int kNPassiveBound = 2;
  const std::array<float, kNPassiveBound> radLBoundary = {1.692612e-01, 8.711904e-02};
  const std::array<float, kNPassiveBound> xrhoBoundary = {6.795774e+00, 3.111401e+00};
  const std::array<float, kNPassiveBound> rBoundary = {50, 70.0}; // cm

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
    addLayer(Form("tpc_boundary%d", i), rBoundary[i], zLength, radLBoundary[i], xrhoBoundary[i], 0); // dummy errors
  }
  for (int k = 0; k < tpcRows; k++) {
    float rowRadius = 0;
    if (k < innerRows) {
      rowRadius = rowOneRadius + k * tpcInnerRadialPitch;
    } else if (k >= innerRows && k < (innerRows + middleRows)) {
      rowRadius = row64Radius + (k - innerRows + 1) * tpcMiddleRadialPitch;
    } else if (k >= (innerRows + middleRows) && k < tpcRows) {
      rowRadius = row128Radius + (k - innerRows - middleRows + 1) * tpcOuterRadialPitch;
    }

    addLayer(Form("tpc_%d", k), rowRadius, zLength, radLPerRow, 0, phiResMean, zResMean, 1.0f, 2);
  }
}

void FastTracker::addGenericDetector(const o2::fastsim::GeometryEntry& configMap, o2::ccdb::BasicCCDBManager* ccdbManager)
{
  // Layers
  for (const auto& layer : configMap.getLayerNames()) {
    if (layer.find("global") != std::string::npos) { // Layers with global tag are skipped
      LOG(info) << " Skipping global configuration entry " << layer;
      continue;
    }

    LOG(info) << " Reading layer " << layer;
    const float r = configMap.getFloatValue(layer, "r");
    LOG(info) << " Layer " << layer << " has radius " << r;
    const float z = configMap.getFloatValue(layer, "z");
    const float x0 = configMap.getFloatValue(layer, "x0");
    const float xrho = configMap.getFloatValue(layer, "xrho");
    const float resRPhi = configMap.getFloatValue(layer, "resRPhi");
    const float resZ = configMap.getFloatValue(layer, "resZ");
    const float eff = configMap.getFloatValue(layer, "eff");
    const int type = configMap.getIntValue(layer, "type");
    const std::string deadPhiRegions = configMap.getValue(layer, "deadPhiRegions", false);

    LOG(info) << " Adding layer " << layer << " r=" << r << " z=" << z << " x0=" << x0 << " xrho=" << xrho << " resRPhi=" << resRPhi << " resZ=" << resZ << " eff=" << eff << " type=" << type << " deadPhiRegions=" << deadPhiRegions;
    DetLayer* addedLayer = addLayer(layer.c_str(), r, z, x0, xrho, resRPhi, resZ, eff, type);
    if (type == DetLayer::kLayerSiliconDisk) {
      const float rMax = configMap.getFloatValue(layer, "diskOuterRadius");
      addedLayer->setSiliconDiskAcceptance(r, rMax, z);
    }
    if (!deadPhiRegions.empty()) { // Taking it as ccdb path or local file
                                   // Check if it begins with ccdb:
      if (deadPhiRegions.starts_with("ccdb:")) {
        std::string ccdbPath = deadPhiRegions.substr(5); // remove "ccdb:" prefix
        if (ccdbManager == nullptr) {
          LOG(fatal) << "CCDB manager is null, cannot retrieve file " << ccdbPath;
          return;
        }
        auto g = ccdbManager->getForTimeStamp<TGraph>(ccdbPath, 1);
        addedLayer->setDeadPhiRegions(g);
      } else {
        // Taking it as local file
        TFile infile(deadPhiRegions.c_str(), "READ");
        if (!infile.IsOpen()) {
          LOG(fatal) << "Cannot open dead phi regions file " << deadPhiRegions;
          return;
        }
        auto g = dynamic_cast<TGraph*>(infile.Get(infile.GetListOfKeys()->At(0)->GetName()));
        infile.Close();
        addedLayer->setDeadPhiRegions(g);
      }
    } else {
      LOG(debug) << " No dead phi regions for layer " << layer;
    }
  }
}

float FastTracker::dist(float z, float r)
{
  // porting of DetektorK::Dist
  // see here:
  // https://github.com/AliceO2Group/DelphesO2/blob/master/src/DetectorK/DetectorK.cxx#L743
  int index = 1;
  int nSteps = 301;
  const int nSigma = 4;
  float dist = 0.0;
  float dz0 = (nSigma * sigmaD - (-nSigma) * sigmaD) / (nSteps - 1);
  float z0 = 0.0;
  for (int i = 0; i < nSteps; i++) {
    if (i == nSteps - 1) {
      index = 1;
    }
    z0 = -4 * sigmaD + i * dz0;
    dist += index * (dz0 / 3.) * (1 / o2::math_utils::sqrt(o2::constants::math::TwoPI) / sigmaD) * std::exp(-z0 * z0 / 2. / sigmaD / sigmaD) * (1 / o2::math_utils::sqrt((z - z0) * (z - z0) + r * r));
    if (index != nSigma) {
      index = nSigma;
    } else {
      index = 2;
    }
  }
  return dist;
}

float FastTracker::oneEventHitDensity(float multiplicity, float radius)
{
  // porting of DetektorK::oneEventHitDensity
  // see here:
  // https://github.com/AliceO2Group/DelphesO2/blob/master/src/DetectorK/DetectorK.cxx#L694
  float den = multiplicity / (o2::constants::math::TwoPI * radius * radius);
  float tg = o2::math_utils::tan(2. * o2::math_utils::atan(std::exp(-avgRapidity))); // theta = 2 atan(exp(-eta))
  den = den / o2::math_utils::sqrt(1 + 1 / (tg * tg));
  return den;
}

float FastTracker::integratedHitDensity(float multiplicity, float radius)
{
  // porting of DetektorK::IntegratedHitDensity
  // see here:
  // https://github.com/AliceO2Group/DelphesO2/blob/master/src/DetectorK/DetectorK.cxx#L712
  float zdcHz = luminosity * 1.e-24 * mCrossSectionMinB; // cm^-2 s^-1 * barn (1 b = 1e-24 cm^2)
  float den = zdcHz * integrationTime / 1000. * multiplicity * dist(0., radius) / (o2::constants::math::TwoPI * radius);
  if (den < oneEventHitDensity(multiplicity, radius)) {
    den = oneEventHitDensity(multiplicity, radius);
  }
  return den;
}

float FastTracker::upcHitDensity(float radius)
{
  // porting of DetektorK::UpcHitDensity
  // see here:
  // https://github.com/AliceO2Group/DelphesO2/blob/master/src/DetectorK/DetectorK.cxx#L727
  float mUPCelectrons = lhcUPCScale * 5456 / (radius * radius) / dNdEtaMinB;
  if (mUPCelectrons < 0) {
    mUPCelectrons = 0.0;
  }
  mUPCelectrons *= integratedHitDensity(dNdEtaMinB, radius);
  mUPCelectrons *= upcBackgroundMultiplier;
  return mUPCelectrons;
}

float FastTracker::hitDensity(float radius)
{
  // porting of DetektorK::HitDensity
  // see here:
  // https://github.com/AliceO2Group/DelphesO2/blob/master/src/DetectorK/DetectorK.cxx#L663
  float arealDensity = 0.;
  if (radius > maxRadiusSlowDet) {
    arealDensity = oneEventHitDensity(dNdEtaCent, radius);
    arealDensity += otherBackground * oneEventHitDensity(dNdEtaMinB, radius);
  }

  // In the version of Delphes used to produce
  // Look-up tables, UpcHitDensity(radius) always returns 0,
  // hence it is left commented out for now
  if (radius < maxRadiusSlowDet) {
    arealDensity = oneEventHitDensity(dNdEtaCent, radius) + integratedHitDensity(dNdEtaMinB, radius);
    // +UpcHitDensity(radius);
    arealDensity += otherBackground * integratedHitDensity(dNdEtaMinB, radius);
  }
  return arealDensity;
}

float FastTracker::probGoodChiSqHit(float radius, float searchRadiusRPhi, float searchRadiusZ)
{
  // porting of DetektorK::ProbGoodChiSqHit
  // see here:
  // https://github.com/AliceO2Group/DelphesO2/blob/master/src/DetectorK/DetectorK.cxx#L629
  const float sx = o2::constants::math::TwoPI * searchRadiusRPhi * searchRadiusZ * hitDensity(radius);
  const float goodHit = 1. / (1 + sx);
  return goodHit;
}

// function to provide a reconstructed track from a perfect input track
// returns number of intercepts (generic for now)
int FastTracker::fastTrack(o2::track::TrackParCov inputTrack, o2::track::TrackParCov& outputTrack, const float nch, const float maxRadius)
{
  dNdEtaCent = nch; // set the number of charged particles per unit rapidity
  hits.clear();
  nIntercepts = 0;
  nSiliconPoints = 0;
  nGasPoints = 0;
  trackIU = -1.f;
  trackOU = -1.f;
  std::array<float, 3> posIni{}; // provision for != PV
  inputTrack.getXYZGlo(posIni);
  const float initialRadius = std::hypot(posIni[0], posIni[1]);
  const float initialX = inputTrack.getX();         // local X of the input track, used for back-propagation
  const float initialAlpha = inputTrack.getAlpha(); // frame of the input track
  const float kTrackingMargin = 0.1;

  int firstActiveLayer = -1; // first supported layer that is not inert
  for (size_t i = 0; i < layers.size(); ++i) {
    if (!layers[i].isInert() && !layers[i].isSiliconDisk()) {
      firstActiveLayer = i;
      break;
    }
  }
  if (firstActiveLayer < 0) {
    LOG(fatal) << "No active non-disk layers found in FastTracker, check layer setup";
    return -2; // no active layers
  }
  const int xrhosteps = 100;
  const bool applyAngularCorrection = true;
  static constexpr float InterceptFailed = 999.f;

  // Delphes sets this to 20 instead of the number of layers,
  // but does not count all points in the tpc as layers which we do here
  // Loop over all the added layers to prevent crash when adding the tpc
  // Should not affect efficiency calculation
  goodHitProbability.clear();
  for (size_t i = 0; i < layers.size(); ++i) {
    goodHitProbability.push_back(-1.);
  }
  goodHitProbability[0] = 1.; // we use layer zero to accumulate

  // +-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+
  // Outward pass to find intercepts
  int firstLayerReached = -1;
  int lastLayerReached = -1;
  int previousLayer = -1;
  new (&outputTrack)(o2::track::TrackParCov)(inputTrack);
  for (size_t il = 0; il < layers.size(); il++) {
    if (layers[il].isSiliconDisk()) {
      continue;
    }
    const int previousNonDiskLayer = previousLayer;
    previousLayer = il;
    // check if layer is doable
    if (layers[il].getRadius() < initialRadius) {
      continue; // this layer should not be attempted, but go ahead
    }

    if (layers[il].getRadius() > maxRadius) {
      if (lastLayerReached == -1) {
        // This means that we didn't reach the first layer
        return -9;
      }
      break; // could not reach
    }

    // check if layer is reached
    float targetX = 1e+3;
    // propagate outward: for a track starting at the vertex both crossings of the cylinder are at about the same
    // distance and DirAuto may pick the one behind the vertex, i.e. propagate the track backwards
    inputTrack.getXatLabR(layers[il].getRadius(), targetX, magneticField, o2::track::DirOutward);
    if (targetX > InterceptFailed) {
      LOGF(debug, "Failed to find intercept for layer %d at radius %.2f cm", il, layers[il].getRadius());
      break; // failed to find intercept
    }

    bool ok = inputTrack.propagateTo(targetX, magneticField);
    if (ok && mApplyMSCorrection && layers[il].getRadiationLength() > 0) {
      ok = inputTrack.correctForMaterial(layers[il].getRadiationLength(), 0, applyAngularCorrection);
    }
    if (ok && mApplyElossCorrection && layers[il].getDensity() > 0) { // correct in small steps
      for (int ise = xrhosteps; ise > 0; --ise) {
        ok = inputTrack.correctForMaterial(0, -layers[il].getDensity() / xrhosteps, applyAngularCorrection);
        if (!ok) {
          break;
        }
      }
    }
    // as in DetectorK::SolveTrack: in the frame of an active layer the track must have |sin(phi)| < mMaxSnp,
    // otherwise the layer is considered not reached (avoids grazing crossings)
    if (ok && !layers[il].isInert()) {
      std::array<float, 3> posLayer{};
      inputTrack.getXYZGlo(posLayer);
      o2::track::TrackParCov trackInLayerFrame(inputTrack);
      ok = trackInLayerFrame.rotate(std::atan2(posLayer[1], posLayer[0])) && std::abs(trackInLayerFrame.getSnp()) < mMaxSnp;
    }
    LOGF(debug, "Propagation was %s up to layer %d", ok ? "successful" : "unsuccessful", il);

    // was there a problem on this layer?
    if (!ok && previousNonDiskLayer >= 0) { // may fail to reach target layer due to the eloss
      const float rad2 = inputTrack.getX() * inputTrack.getX() + inputTrack.getY() * inputTrack.getY();
      const float maxR = layers[previousNonDiskLayer].getRadius() + kTrackingMargin * 2;
      const float minRad = (fMinRadTrack > 0 && fMinRadTrack < maxR) ? fMinRadTrack : maxR;
      if (rad2 - minRad * minRad < kTrackingMargin * kTrackingMargin) { // check previously reached layer
        return -5;                                                      // did not reach min requested layer
      }
      break;
    }

    if (std::abs(inputTrack.getZ()) > layers[il].getZ() && mApplyZacceptance) {
      break; // out of acceptance bounds
    }

    if (layers[il].isInert()) {
      if (mVerboseLevel > 0) {
        LOG(info) << "Skipping inert layer: " << layers[il].getName() << " at radius " << layers[il].getRadius() << " cm";
      }
      continue; // inert layer, skip
    }

    if (layers[il].isInDeadPhiRegion(inputTrack.getPhi())) {
      LOGF(debug, "Track is in dead region of layer %d", il);
      continue; // dead region, skip
    }

    // layer is reached
    if (firstLayerReached < 0) {
      LOGF(debug, "First layer reached: %d", il);
      firstLayerReached = il;
    }
    lastLayerReached = il;
    if (layers[il].isTOF()) {
      continue; // TOF layers are not used for tracking, but are kept in range for the innermost/outermost update
    }
    nIntercepts++;
  }

  // +-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+
  // initialize track at outer point
  o2::track::TrackParCov inwardTrack(inputTrack);

  // Enlarge covariance matrix
  static constexpr float LargeErr2Coord = 5 * 5;
  static constexpr float LargeErr2Dir = 0.7 * 0.7;
  static constexpr float LargeErr2PtI = 30.5 * 30.5;
  auto resetToLargeCovariance = [](o2::track::TrackParCov& trk) {
    std::array<float, o2::track::kCovMatSize> largeCov{}; // zero-initialised
    largeCov[o2::track::CovLabels::kSigY2] = largeCov[o2::track::CovLabels::kSigZ2] = LargeErr2Coord;
    largeCov[o2::track::CovLabels::kSigSnp2] = largeCov[o2::track::CovLabels::kSigTgl2] = LargeErr2Dir;
    largeCov[o2::track::CovLabels::kSigQ2Pt2] = LargeErr2PtI * trk.getQ2Pt() * trk.getQ2Pt();
    trk.setCov(largeCov);
    trk.checkCovariance();
  };
  resetToLargeCovariance(inwardTrack);

  // as in DetectorK, radius and search radii are scaled by 100 in the good-hit probability: the scaling
  // cancels in the probability, but it means that all layers are treated as fast detectors (r * 100 > maxRadiusSlowDet)
  static constexpr float DelphesScale = 100.f;

  // bookkeeping of the layer crossings of the inward pass, re-used by the second outward pass
  struct LayerCrossing {
    bool crossed = false;    // the layer was crossed in the inward pass
    bool trackerHit = false; // a tracker measurement was added on this layer
    bool tofHit = false;     // TOF layer crossed (probability computed, but no measurement)
    float phi = 0.f;         // azimuth of the crossing point (frame of the layer)
    float r = 0.f;           // radius of the crossing point
    float z = 0.f;           // z of the crossing point
    float sigY2Inw = 0.f;    // inward prediction before the update
    float sigZ2Inw = 0.f;
  };
  std::vector<LayerCrossing> crossings(layers.size());

  // +-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+
  // Inward pass to calculate covariances
  // As in DetectorK::SolveTrack, the track is rotated to the frame of each layer (X axis normal to the layer)
  for (int il = lastLayerReached; il >= firstLayerReached; il--) {
    if (layers[il].isSiliconDisk()) {
      continue;
    }

    float targetX = 1e+3;
    inputTrack.getXatLabR(layers[il].getRadius(), targetX, magneticField);
    if (targetX > InterceptFailed) {
      continue; // failed to find intercept
    }

    if (!inputTrack.propagateTo(targetX, magneticField)) {
      continue; // failed to propagate
    }

    if (std::abs(inputTrack.getZ()) > layers[il].getZ() && mApplyZacceptance) {
      continue; // out of acceptance bounds but continue inwards
    }

    // get perfect data point position
    std::array<float, 3> spacePoint{};
    inputTrack.getXYZGlo(spacePoint);
    std::vector<float> thisHit = {spacePoint[0], spacePoint[1], spacePoint[2]};
    const float phiHit = std::atan2(spacePoint[1], spacePoint[0]);
    const float rHit = std::hypot(spacePoint[0], spacePoint[1]);

    // move to the frame of the layer at the crossing point, where the hit is at local (rHit, 0, z)
    if (!inwardTrack.rotate(phiHit) || !inwardTrack.propagateTo(rHit, magneticField)) {
      continue;
    }

    // dead phi regions produce no hit (as in the outward pass), but their material is still traversed
    const bool inDeadRegion = layers[il].isInDeadPhiRegion(inputTrack.getPhi());
    const bool isTrackerHit = !layers[il].isInert() && !layers[il].isTOF() && !inDeadRegion;

    LayerCrossing& crossing = crossings[il];
    crossing.crossed = true;
    crossing.trackerHit = isTrackerHit;
    crossing.tofHit = layers[il].isTOF() && !inDeadRegion;
    crossing.phi = phiHit;
    crossing.r = rHit;
    crossing.z = spacePoint[2];
    // save the inward prediction before the update (DetectorK: saveParInward)
    crossing.sigY2Inw = inwardTrack.getSigmaY2();
    crossing.sigZ2Inw = inwardTrack.getSigmaZ2();

    if (isTrackerHit) { // only update covm for tracker hits
      const o2::track::TrackParametrization<float>::dim2_t hitpoint = {0.f, spacePoint[2]};
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
      for (int ise = xrhosteps; ise > 0; --ise) { // correct in small steps
        if (!inputTrack.correctForMaterial(0, layers[il].getDensity() / xrhosteps, applyAngularCorrection)) {
          return -7;
        }
        if (!inwardTrack.correctForMaterial(0, layers[il].getDensity() / xrhosteps, applyAngularCorrection)) {
          return -7;
        }
      }
    }

    if (inDeadRegion) {
      continue; // no hit in dead region, material already accounted for
    }

    if (layers[il].isSilicon()) {
      nSiliconPoints++; // count silicon hits
    }
    if (layers[il].isGas()) {
      nGasPoints++; // count TPC/gas hits
    }

    hits.push_back(thisHit);
    if (!layers[il].isInert()) {
      trackIU = std::hypot(spacePoint[0], spacePoint[1]);
      if (trackOU < 0.f) {
        trackOU = trackIU;
      }
    }
  }

  // +-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+
  // Second outward pass for the good hit probabilities (DetectorK "backward tracking"):
  // the search window at each layer is the combination of the inward and outward predictions,
  // both taken before the update with the hit on that layer
  o2::track::TrackParCov outwardTrack(inwardTrack);
  resetToLargeCovariance(outwardTrack);
  for (int il = firstLayerReached; il <= lastLayerReached; il++) {
    const LayerCrossing& crossing = crossings[il];
    if (!crossing.crossed) {
      continue;
    }
    if (!outwardTrack.rotate(crossing.phi) || !outwardTrack.propagateTo(crossing.r, magneticField)) {
      LOGF(debug, "Second outward pass: failed to reach layer %d", il);
      break; // the remaining layers keep goodHitProbability = -1 (not used)
    }

    if (crossing.trackerHit || crossing.tofHit) {
      const float sigY2Out = outwardTrack.getSigmaY2();
      const float sigZ2Out = outwardTrack.getSigmaZ2();
      const float sigY2Cmb = crossing.sigY2Inw * sigY2Out / (crossing.sigY2Inw + sigY2Out);
      const float sigZ2Cmb = crossing.sigZ2Inw * sigZ2Out / (crossing.sigZ2Inw + sigZ2Out);
      const float sigYCmb = o2::math_utils::sqrt(sigY2Cmb + layers[il].getResolutionRPhi() * layers[il].getResolutionRPhi());
      const float sigZCmb = o2::math_utils::sqrt(sigZ2Cmb + layers[il].getResolutionZ() * layers[il].getResolutionZ());
      const float layerGoodHitProbability = probGoodChiSqHit(layers[il].getRadius() * DelphesScale, sigYCmb * DelphesScale, sigZCmb * DelphesScale);
      if (il > 0) { // index 0 holds the accumulated product, do not overwrite it
        goodHitProbability[il] = layerGoodHitProbability;
      }
      if (crossing.trackerHit) { // TOF layers are not part of the tracking efficiency
        goodHitProbability[0] *= layerGoodHitProbability;
      }
    }

    if (crossing.trackerHit) {
      const o2::track::TrackParametrization<float>::dim2_t hitpoint = {0.f, crossing.z};
      const o2::track::TrackParametrization<float>::dim3_t hitpointcov = {layers[il].getResolutionRPhi() * layers[il].getResolutionRPhi(), 0.f, layers[il].getResolutionZ() * layers[il].getResolutionZ()};
      outwardTrack.update(hitpoint, hitpointcov);
      outwardTrack.checkCovariance();
    }

    bool ok = true;
    if (mApplyMSCorrection && layers[il].getRadiationLength() > 0) {
      ok = outwardTrack.correctForMaterial(layers[il].getRadiationLength(), 0, applyAngularCorrection);
    }
    if (ok && mApplyElossCorrection && layers[il].getDensity() > 0) {
      for (int ise = xrhosteps; ise > 0; --ise) { // correct in small steps, energy loss going outward
        ok = outwardTrack.correctForMaterial(0, -layers[il].getDensity() / xrhosteps, applyAngularCorrection);
        if (!ok) {
          break;
        }
      }
    }
    if (!ok) {
      LOGF(debug, "Second outward pass: material correction failed at layer %d", il);
      break;
    }
  }

  // backpropagate to the original position: rotate back to the frame of the input track and
  // propagate to its initial local X. Searching for an intercept at initialRadius fails when the
  // refitted track's DCA exceeds it (e.g. initialRadius ~ 0 for tracks at the PV)
  if (!inwardTrack.rotate(initialAlpha) || !inwardTrack.propagateTo(initialX, magneticField)) {
    return -4; // failed to propagate
  }

  // only attempt to continue if intercepts are at least four
  static constexpr int MinIntercepts = 4;
  if (nIntercepts < MinIntercepts) {
    return nIntercepts;
  }

  // generate efficiency: product of the good-hit probabilities, accumulated in index 0
  const float eff = goodHitProbability[0];
  if (mApplyEffCorrection) {
    if (gRandom->Uniform() > eff) {
      return -8;
    }
  }

  outputTrack.setCov(inwardTrack.getCov());
  outputTrack.checkCovariance();

  // Use covariance matrix based smearing
  std::array<float, o2::track::kCovMatSize> covMat = {0.};
  for (int ii = 0; ii < o2::track::kCovMatSize; ii++) {
    covMat[ii] = outputTrack.getCov()[ii];
  }

  std::array<std::array<double, o2::track::kNParams>, o2::track::kNParams> fcovm{}; // double precision is needed for regularisation
  TMatrixDSym m(o2::track::kNParams);
  for (int ii = 0, k = 0; ii < o2::track::kNParams; ++ii) {
    for (int jj = 0; jj < ii + 1; ++jj, ++k) {
      fcovm[ii][jj] = covMat[k];
      fcovm[jj][ii] = covMat[k];
    }
  }

  // evaluate ruben's conditional, regularise
  static constexpr float CovMatFactorThreshold = -1e-5;
  const bool makePositiveDefinite = (covMatFactor > CovMatFactorThreshold); // apply fix
  bool rubenConditional = false;
  for (int ii = 0; ii < o2::track::kNParams; ii++) {
    for (int jj = 0; jj < o2::track::kNParams; jj++) {
      if (ii == jj) {
        continue; // don't evaluate diagonals
      }
      if (fcovm[ii][jj] * fcovm[ii][jj] > std::abs(fcovm[ii][ii] * fcovm[jj][jj])) {
        rubenConditional = true;
        if (makePositiveDefinite) {
          fcovm[ii][jj] = TMath::Sign(1, fcovm[ii][jj]) * covMatFactor * std::sqrt(std::abs(fcovm[ii][ii] * fcovm[jj][jj]));
        }
      }
    }
  }

  // Should have a valid cov matrix now
  m.SetMatrixArray(fcovm[0].data());
  TMatrixDSymEigen eigen(m);
  TMatrixD eigVec = eigen.GetEigenVectors();
  const TVectorD& eigVal = eigen.GetEigenValues();
  bool negEigVal = false;
  for (int ii = 0; ii < o2::track::kNParams; ii++) {
    if (eigVal[ii] < 0.0f) {
      negEigVal = true;
    }
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
  std::array<float, o2::track::kNParams> transformedParams{};
  for (int ii = 0; ii < o2::track::kNParams; ++ii) {
    float val = 0.;
    for (int jj = 0; jj < o2::track::kNParams; ++jj) {
      val += eigVec[jj][ii] * outputTrack.getParam(jj);
    }
    // smear parameters according to eigenvalues
    transformedParams[ii] = gRandom->Gaus(val, std::sqrt(eigVal[ii]));
  }

  // invert eigenvector matrix
  eigVec.Invert();
  // transform back params vector
  for (int ii = 0; ii < o2::track::kNParams; ++ii) {
    float val = 0.;
    for (int jj = 0; jj < o2::track::kNParams; ++jj) {
      val += eigVec[jj][ii] * transformedParams[jj];
    }
    outputTrack.setParam(val, ii);
  }
  // should make a sanity check that par[2] sin(phi) is in [-1, 1]
  if (std::abs(outputTrack.getParam(2)) > 1.) {
    LOG(info) << " --- smearTrack failed sin(phi) sanity check: " << outputTrack.getParam(2);
    return -2;
  }

  return nIntercepts;
}
// +-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+

} // namespace o2::fastsim

ClassImp(o2::fastsim::FastTracker);
