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

#include <vector>
#include "TMath.h"
#include "TMatrixD.h"
#include "TRandom.h"
#include "TMatrixDSymEigen.h"
#include "FastTracker.h"

namespace o2
{
namespace fastsim
{

// +-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+

FastTracker::FastTracker()
{
  // base constructor
  magneticField = 20; // in kiloGauss
  applyZacceptance = false;
  covMatFactor = 0.99f;
  verboseLevel = 0;

  // last fast-tracked track properties
  covMatOK = 0;
  covMatNotOK = 0;
  nIntercepts = 0;
  nSiliconPoints = 0;
  nGasPoints = 0;
}

void FastTracker::AddLayer(TString name, float r, float z, float x0, float xrho, float resRPhi, float resZ, float eff, int type)
{
  DetLayer newLayer{name.Data(), r, z, x0, xrho, resRPhi, resZ, eff, type};
  layers.push_back(newLayer);
}

DetLayer FastTracker::GetLayer(int layer, bool ignoreBarrelLayers)
{
  int layerIdx = layer;
  if (ignoreBarrelLayers) {
    for (int il = 0, trackingLayerIdx = 0; trackingLayerIdx <= layer; il++) {
      if (layers[il].type == 0)
        continue;
      trackingLayerIdx++;
      layerIdx = il;
    }
  }
  return layers[layerIdx];
}

void FastTracker::Print()
{
  // print out layer setup
  LOG(info) << "+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+";
  LOG(info) << " Printing detector layout with " << layers.size() << " effective elements: ";
  for (uint32_t il = 0; il < layers.size(); il++) {
    LOG(info) << " Layer #" << il << "\t" << layers[il].name.Data() << "\tr = " << Form("%.2f", layers[il].r) << "cm\tz = " << layers[il].z << "\t"
              << "x0 = " << layers[il].x0 << "\txrho = " << layers[il].xrho << "\tresRPhi = " << layers[il].resRPhi << "\tresZ = " << layers[il].resZ << "\teff = " << layers[il].eff;
  }
  LOG(info) << "+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+";
}

void FastTracker::AddSiliconALICE3v4()
{
  LOG(info) << " Adding ALICE 3 v4 ITS layers";
  float x0IT = 0.001;        // 0.1%
  float x0OT = 0.005;        // 0.5%
  float xrhoIB = 1.1646e-02; // 50 mum Si
  float xrhoOB = 1.1646e-01; // 500 mum Si

  float resRPhiIT = 0.00025; //  2.5 mum
  float resZIT = 0.00025;    //  2.5 mum
  float resRPhiOT = 0.0005;  // 5 mum
  float resZOT = 0.0005;     // 5 mum
  float eff = 1.00;

  layers.push_back(DetLayer{"bpipe0", 0.48, 250, 0.00042, 2.772e-02, 0.0f, 0.0f, 0.0f, 0}); // 150 mum Be
  layers.push_back(DetLayer{"ddd0", 0.5, 250, x0IT, xrhoIB, resRPhiIT, resZIT, eff, 1});
  layers.push_back(DetLayer{"ddd1", 1.2, 250, x0IT, xrhoIB, resRPhiIT, resZIT, eff, 1});
  layers.push_back(DetLayer{"ddd2", 2.5, 250, x0IT, xrhoIB, resRPhiIT, resZIT, eff, 1});
  layers.push_back(DetLayer{"bpipe1", 5.7, 250, 0.0014, 9.24e-02, 0.0f, 0.0f, 0.0f, 0}); // 500 mum Be
  layers.push_back(DetLayer{"ddd3", 7., 250, x0IT, xrhoIB, resRPhiIT, resZIT, eff, 1});
  layers.push_back(DetLayer{"ddd4", 10., 250, x0OT, xrhoOB, resRPhiOT, resZOT, eff, 1});
  layers.push_back(DetLayer{"ddd5", 13., 250, x0OT, xrhoOB, resRPhiOT, resZOT, eff, 1});
  layers.push_back(DetLayer{"ddd6", 16., 250, x0OT, xrhoOB, resRPhiOT, resZOT, eff, 1});
  layers.push_back(DetLayer{"ddd7", 25., 250, x0OT, xrhoOB, resRPhiOT, resZOT, eff, 1});
  layers.push_back(DetLayer{"ddd8", 40., 250, x0OT, xrhoOB, resRPhiOT, resZOT, eff, 1});
  layers.push_back(DetLayer{"ddd9", 45., 250, x0OT, xrhoOB, resRPhiOT, resZOT, eff, 1});
}

void FastTracker::AddSiliconALICE3v1()
{
  LOG(info) << " Adding ALICE 3 v1 ITS layers";
  float x0IT = 0.001;        // 0.1%
  float x0OT = 0.005;        // 0.5%
  float xrhoIB = 2.3292e-02; // 100 mum Si
  float xrhoOB = 2.3292e-01; // 1000 mum Si

  float resRPhiIT = 0.00025; //  2.5 mum
  float resZIT = 0.00025;    //  2.5 mum
  float resRPhiOT = 0.00100; // 5 mum
  float resZOT = 0.00100;    // 5 mum
  float eff = 1.00;

  layers.push_back(DetLayer{"bpipe0", 0.48, 250, 0.00042, 2.772e-02, 0.0f, 0.0f, 0.0f, 0}); // 150 mum Be
  layers.push_back(DetLayer{"B00", 0.5, 250, x0IT, xrhoIB, resRPhiIT, resZIT, eff, 1});
  layers.push_back(DetLayer{"B01", 1.2, 250, x0IT, xrhoIB, resRPhiIT, resZIT, eff, 1});
  layers.push_back(DetLayer{"B02", 2.5, 250, x0IT, xrhoIB, resRPhiIT, resZIT, eff, 1});
  layers.push_back(DetLayer{"bpipe1", 3.7, 250, 0.0014, 9.24e-02, 0.0f, 0.0f, 0.0f, 0}); // 500 mum Be
  layers.push_back(DetLayer{"B03", 3.75, 250, x0IT, xrhoIB, resRPhiIT, resZIT, eff, 1});
  layers.push_back(DetLayer{"B04", 7., 250, x0OT, xrhoOB, resRPhiOT, resZOT, eff, 1});
  layers.push_back(DetLayer{"B05", 12., 250, x0OT, xrhoOB, resRPhiOT, resZOT, eff, 1});
  layers.push_back(DetLayer{"B06", 20., 250, x0OT, xrhoOB, resRPhiOT, resZOT, eff, 1});
  layers.push_back(DetLayer{"B07", 30., 250, x0OT, xrhoOB, resRPhiOT, resZOT, eff, 1});
  layers.push_back(DetLayer{"B08", 45., 250, x0OT, xrhoOB, resRPhiOT, resZOT, eff, 1});
  layers.push_back(DetLayer{"B09", 60., 250, x0OT, xrhoOB, resRPhiOT, resZOT, eff, 1});
  layers.push_back(DetLayer{"B10", 80., 250, x0OT, xrhoOB, resRPhiOT, resZOT, eff, 1});
}

void FastTracker::AddTPC(float phiResMean, float zResMean)
{
  LOG(info) << " Adding standard time projection chamber";

  // porting of DetectorK::AddTPC
  // see here:
  // https://github.com/AliceO2Group/DelphesO2/blob/master/src/DetectorK/DetectorK.cxx
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

// function to provide a reconstructed track from a perfect input track
// returns number of intercepts (generic for now)
int FastTracker::FastTrack(o2::track::TrackParCov inputTrack, o2::track::TrackParCov& outputTrack)
{
  hits.clear();
  nIntercepts = 0;
  nSiliconPoints = 0;
  nGasPoints = 0;
  std::array<float, 3> posIni; // provision for != PV
  inputTrack.getXYZGlo(posIni);
  float initialRadius = std::hypot(posIni[0], posIni[1]);
  float kTrackingMargin = 0.1;
  int xrhosteps = 100;
  bool applyAngularCorrection = true;

  // +-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+
  // Outward pass to find intercepts
  int firstLayerReached = -1;
  int lastLayerReached = -1;
  new (&outputTrack)(o2::track::TrackParCov)(inputTrack);
  for (uint32_t il = 0; il < layers.size(); il++) {
    // check if layer is doable
    if (layers[il].r < initialRadius)
      continue; // this layer should not be attempted, but go ahead

    // check if layer is reached
    float targetX = 1e+3;
    bool ok = true;
    inputTrack.getXatLabR(layers[il].r, targetX, magneticField);
    if (targetX > 999)
      break; // failed to find intercept

    ok = inputTrack.propagateTo(targetX, magneticField);
    if (ok && applyMSCorrection && layers[il].x0 > 0) {
      ok = inputTrack.correctForMaterial(layers[il].x0, 0, applyAngularCorrection);
    }
    if (ok && applyElossCorrection && layers[il].xrho > 0) { // correct in small steps
      for (int ise = xrhosteps; ise--;) {
        ok = inputTrack.correctForMaterial(0, -layers[il].xrho / xrhosteps, applyAngularCorrection);
        if (!ok)
          break;
      }
    }

    // was there a problem on this layer?
    if (!ok && il > 0) { // may fail to reach target layer due to the eloss
      float rad2 = inputTrack.getX() * inputTrack.getX() + inputTrack.getY() * inputTrack.getY();
      float fMinRadTrack = 132.;
      float maxR = layers[il - 1].r + kTrackingMargin * 2;
      float minRad = (fMinRadTrack > 0 && fMinRadTrack < maxR) ? fMinRadTrack : maxR;
      if (rad2 - minRad * minRad < kTrackingMargin * kTrackingMargin) { // check previously reached layer
        return -5;                                                      // did not reach min requested layer
      } else {
        break;
      }
    }
    if (std::abs(inputTrack.getZ()) > layers[il].z && applyZacceptance) {
      break; // out of acceptance bounds
    }

    if (layers[il].type == 0)
      continue; // inert layer, skip

    // layer is reached
    if (firstLayerReached < 0)
      firstLayerReached = il;
    lastLayerReached = il;
    nIntercepts++;
  }

  // +-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+
  // initialize track at outer point
  o2::track::TrackParCov inwardTrack(inputTrack);

  // Enlarge covariance matrix
  std::array<float, 5> trPars = {0.};
  for (int ip = 0; ip < 5; ip++) {
    trPars[ip] = outputTrack.getParam(ip);
  }
  std::array<float, 15> largeCov = {0.};
  enum { kY,
         kZ,
         kSnp,
         kTgl,
         kPtI }; // track parameter aliases
  enum { kY2,
         kYZ,
         kZ2,
         kYSnp,
         kZSnp,
         kSnp2,
         kYTgl,
         kZTgl,
         kSnpTgl,
         kTgl2,
         kYPtI,
         kZPtI,
         kSnpPtI,
         kTglPtI,
         kPtI2 }; // cov.matrix aliases
  const double kLargeErr2Coord = 5 * 5;
  const double kLargeErr2Dir = 0.7 * 0.7;
  const double kLargeErr2PtI = 30.5 * 30.5;
  for (int ic = 15; ic--;)
    largeCov[ic] = 0.;
  largeCov[kY2] = largeCov[kZ2] = kLargeErr2Coord;
  largeCov[kSnp2] = largeCov[kTgl2] = kLargeErr2Dir;
  largeCov[kPtI2] = kLargeErr2PtI * trPars[kPtI] * trPars[kPtI];

  inwardTrack.setCov(largeCov);
  inwardTrack.checkCovariance();

  // +-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+
  // Inward pass to calculate covariances
  for (int il = lastLayerReached; il >= firstLayerReached; il--) {

    float targetX = 1e+3;
    inputTrack.getXatLabR(layers[il].r, targetX, magneticField);
    if (targetX > 999)
      continue; // failed to find intercept

    if (!inputTrack.propagateTo(targetX, magneticField)) {
      continue; // failed to propagate
    }

    if (std::abs(inputTrack.getZ()) > layers[il].z && applyZacceptance) {
      continue; // out of acceptance bounds but continue inwards
    }

    // get perfect data point position
    std::array<float, 3> spacePoint;
    inputTrack.getXYZGlo(spacePoint);
    std::vector<float> thisHit = {spacePoint[0], spacePoint[1], spacePoint[2]};

    // towards adding cluster: move to track alpha
    double alpha = inwardTrack.getAlpha();
    double xyz1[3]{
      TMath::Cos(alpha) * spacePoint[0] + TMath::Sin(alpha) * spacePoint[1],
      -TMath::Sin(alpha) * spacePoint[0] + TMath::Cos(alpha) * spacePoint[1],
      spacePoint[2]};
    if (!inwardTrack.propagateTo(xyz1[0], magneticField))
      continue;

    if (layers[il].type != 0) { // only update covm for tracker hits
      const o2::track::TrackParametrization<float>::dim2_t hitpoint = {
        static_cast<float>(xyz1[1]),
        static_cast<float>(xyz1[2])};
      const o2::track::TrackParametrization<float>::dim3_t hitpointcov = {layers[il].resRPhi * layers[il].resRPhi, 0.f, layers[il].resZ * layers[il].resZ};

      inwardTrack.update(hitpoint, hitpointcov);
      inwardTrack.checkCovariance();
    }

    if (applyMSCorrection && layers[il].x0 > 0) {
      if (!inputTrack.correctForMaterial(layers[il].x0, 0, applyAngularCorrection)) {
        return -6;
      }
      if (!inwardTrack.correctForMaterial(layers[il].x0, 0, applyAngularCorrection)) {
        return -6;
      }
    }
    if (applyElossCorrection && layers[il].xrho > 0) {
      for (int ise = xrhosteps; ise--;) { // correct in small steps
        if (!inputTrack.correctForMaterial(0, layers[il].xrho / xrhosteps, applyAngularCorrection)) {
          return -7;
        }
        if (!inwardTrack.correctForMaterial(0, layers[il].xrho / xrhosteps, applyAngularCorrection)) {
          return -7;
        }
      }
    }

    if (layers[il].type == 1)
      nSiliconPoints++; // count silicon hits
    if (layers[il].type == 2)
      nGasPoints++; // count TPC/gas hits

    hits.push_back(thisHit);
  }

  // backpropagate to original radius
  float finalX = 1e+3;
  inwardTrack.getXatLabR(initialRadius, finalX, magneticField);
  if (finalX > 999)
    return -3; // failed to find intercept

  if (!inwardTrack.propagateTo(finalX, magneticField)) {
    return -4; // failed to propagate
  }

  // only attempt to continue if intercepts are at least four
  if (nIntercepts < 4)
    return nIntercepts;

  outputTrack.setCov(inwardTrack.getCov());
  outputTrack.checkCovariance();

  // Use covariance matrix based smearing
  std::array<double, 15> covMat = {0.};
  for (int ii = 0; ii < 15; ii++)
    covMat[ii] = outputTrack.getCov()[ii];
  TMatrixDSym m(5);
  double fcovm[5][5];

  for (int ii = 0, k = 0; ii < 5; ++ii) {
    for (int j = 0; j < ii + 1; ++j, ++k) {
      fcovm[ii][j] = covMat[k];
      fcovm[j][ii] = covMat[k];
    }
  }

  // evaluate ruben's conditional, regularise
  bool makePositiveDefinite = (covMatFactor > -1e-5); // apply fix
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
    if (verboseLevel > 0) {
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
  double params_[5];
  for (int ii = 0; ii < 5; ++ii) {
    double val = 0.;
    for (int j = 0; j < 5; ++j)
      val += eigVec[j][ii] * outputTrack.getParam(j);
    // smear parameters according to eigenvalues
    params_[ii] = gRandom->Gaus(val, sqrt(eigVal[ii]));
  }

  // invert eigenvector matrix
  eigVec.Invert();
  // transform back params vector
  for (int ii = 0; ii < 5; ++ii) {
    double val = 0.;
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
