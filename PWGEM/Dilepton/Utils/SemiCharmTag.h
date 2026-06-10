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

/// \event mixing handler
/// \author daiki.sekihata@cern.ch
/// see https://arxiv.org/abs/2604.11574 invented by LHCb

#ifndef PWGEM_DILEPTON_UTILS_SEMICHARMTAG_H_
#define PWGEM_DILEPTON_UTILS_SEMICHARMTAG_H_

#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/Logger.h>
#include <ReconstructionDataFormats/DCA.h>
#include <ReconstructionDataFormats/PID.h>
#include <ReconstructionDataFormats/V0.h>

#include <Math/Vector4D.h> // IWYU pragma: keep (do not replace with Math/Vector4Dfwd.h)
#include <Math/Vector4Dfwd.h>

#include <array>

namespace o2::aod::pwgem::dilepton::utils
{

struct LHPair { // struct to store electron-hadron pair information
  float mass{-999.f};
  float pt{-999.f};
  float p{-999.f};
  float deta{-999.f};
  float dphi{-999.f};
  float dca2legs{-999.f};
  float chi2PCA{-999.f};
  float cospa{-999.f};
  float cospaXY{-999.f};
  float cospaRZ{-999.f};
  float lxy{-999.f};
  float lz{-999.f};
  float lxyz{-999.f};
  float lxyErr{-999.f};
  float lzErr{-999.f};
  float lxyzErr{-999.f};
  float impParXY{-999.f};
  float impParZ{-999.f};
  float impParCYY{-999.f};
  float impParCZY{-999.f};
  float impParCZZ{-999.f};

  // float ptSVL{-999.f};
  // float plSVL{-999.f};
  // float ptSVH{-999.f};
  // float plSVH{-999.f};

  // float ptFDL{-999.f};
  // float plFDL{-999.f};
  // float ptFDH{-999.f};
  // float plFDH{-999.f};

  // float ptFD{-999.f};
  // float plFD{-999.f};

  bool isOK{false};
};

template <typename TFitter, typename TCollision, typename TLepton, typename TTrack>
LHPair makePairLeptonTrack(TFitter& fitter, TCollision const& collision, TLepton const& lepton, TTrack const& track, o2::track::PID::ID leptonId, o2::track::PID::ID strHadId)
{
  LHPair pair;
  auto leptonParCov = getTrackParCov(lepton);
  leptonParCov.setPID(leptonId);

  auto trackParCov = getTrackParCov(track);
  trackParCov.setPID(strHadId);

  const std::array<float, 3> vertex = {collision.posX(), collision.posY(), collision.posZ()};
  std::array<float, 3> svpos = {0.}; // secondary vertex position
  std::array<float, 3> pvec0 = {0.};
  std::array<float, 3> pvec1 = {0.};

  int nCand = 0;
  try {
    nCand = fitter.process(leptonParCov, trackParCov);
  } catch (...) {
    LOG(error) << "Exception caught in DCA fitter process call!";
    pair.isOK = false;
    return pair;
  }
  if (nCand == 0) {
    pair.isOK = false;
    return pair;
  }

  // fitter.propagateTracksToVertex();
  const auto& vtx = fitter.getPCACandidate();
  for (int i = 0; i < 3; i++) {
    svpos[i] = vtx[i];
  }
  // const std::array<float, 3> coordVtxLK = df.getPCACandidatePos();

  fitter.getTrack(0).getPxPyPzGlo(pvec0); // lepton
  fitter.getTrack(1).getPxPyPzGlo(pvec1); // track
  std::array<float, 3> pvecSum = {pvec0[0] + pvec1[0], pvec0[1] + pvec1[1], pvec0[2] + pvec1[2]};

  // pair.ptSVL = RecoDecay::sqrtSumOfSquares(pvec0[0], pvec0[1]);
  // pair.ptSVH = RecoDecay::sqrtSumOfSquares(pvec1[0], pvec1[1]);
  // pair.plSVL = pvec0[2];
  // pair.plSVH = pvec1[2];

  pair.cospa = RecoDecay::cpa(vertex, svpos, pvecSum);
  pair.cospaXY = RecoDecay::cpaXY(vertex, svpos, pvecSum);
  pair.cospaRZ = RecoDecay::cpaRZ(vertex, svpos, pvecSum);
  pair.dca2legs = std::sqrt(fitter.getChi2AtPCACandidate());
  pair.chi2PCA = fitter.getChi2AtPCACandidate();
  pair.lxy = std::sqrt(std::pow(svpos[0] - collision.posX(), 2) + std::pow(svpos[1] - collision.posY(), 2));
  pair.lz = svpos[2] - collision.posZ();
  pair.lxyz = std::sqrt(std::pow(svpos[0] - collision.posX(), 2) + std::pow(svpos[1] - collision.posY(), 2) + std::pow(svpos[2] - collision.posZ(), 2));

  auto primaryVertex = getPrimaryVertex(collision);
  std::array<float, 6> covVtxLH = fitter.calcPCACovMatrixFlat();
  double phiLH{}, thetaLH{};
  getPointDirection(std::array{primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ()}, svpos, phiLH, thetaLH);
  pair.lxyzErr = std::sqrt(getRotatedCovMatrixXX(primaryVertex.getCov(), phiLH, thetaLH) + getRotatedCovMatrixXX(covVtxLH, phiLH, thetaLH));
  pair.lxyErr = std::sqrt(getRotatedCovMatrixXX(primaryVertex.getCov(), phiLH, 0.) + getRotatedCovMatrixXX(covVtxLH, phiLH, 0.));
  pair.lzErr = std::sqrt(getRotatedCovMatrixXX(primaryVertex.getCov(), 0, thetaLH) + getRotatedCovMatrixXX(covVtxLH, 0, thetaLH));

  // std::array<float, 3> uvFD = {(svpos[0] - collision.posX()) / pair.lxyz, (svpos[1] - collision.posY()) / pair.lxyz, (svpos[2] - collision.posZ()) / pair.lxyz}; // unit vector of flight direction
  // pair.plFD = RecoDecay::dotProd(pvecSum, uvFD);
  // pair.ptFD = RecoDecay::sqrtSumOfSquares(pvecSum[0] - pair.plFD * uvFD[0], pvecSum[1] - pair.plFD * uvFD[1], pvecSum[2] - pair.plFD * uvFD[2]);

  // pair.plFDL = RecoDecay::dotProd(pvec0, uvFD);
  // pair.plFDH = RecoDecay::dotProd(pvec1, uvFD);
  // pair.ptFDL = RecoDecay::sqrtSumOfSquares(pvec0[0] - pair.plFDL * uvFD[0], pvec0[1] - pair.plFDL * uvFD[1], pvec0[2] - pair.plFDL * uvFD[2]);
  // pair.ptFDH = RecoDecay::sqrtSumOfSquares(pvec1[0] - pair.plFDH * uvFD[0], pvec1[1] - pair.plFDH * uvFD[1], pvec1[2] - pair.plFDH * uvFD[2]);

  // // propagate the 2 prongs to the secondary vertex
  // leptonParCov.propagateTo(vtx[0], fitter.getBz());
  // trackParCov.propagateTo(vtx[0], fitter.getBz());

  // // calculate impact parameter
  // o2::dataformats::DCA dcaLH;
  // auto trackParCovLH = o2::dataformats::V0(fitter.getPCACandidatePos(), pvecSum, fitter.calcPCACovMatrixFlat(), leptonParCov, trackParCov);
  // trackParCovLH.propagateToDCA(primaryVertex, fitter.getBz(), &dcaLH);

  // pair.impParXY = dcaLH.getY();
  // pair.impParZ = dcaLH.getZ();
  // pair.impParCYY = dcaLH.getSigmaY2();
  // pair.impParCZY = dcaLH.getSigmaYZ();
  // pair.impParCZZ = dcaLH.getSigmaZ2();

  // LOGF(info, "fitter.getBz() = %f, dcaLH.getY() = %f, dcaLH.getZ() = %f", fitter.getBz(), dcaLH.getY(), dcaLH.getZ());

  ROOT::Math::PxPyPzMVector v1(pvec0[0], pvec0[1], pvec0[2], o2::constants::physics::MassElectron);
  if (leptonId == o2::track::PID::Electron) {
    v1.SetM(o2::constants::physics::MassElectron);
  } else if (leptonId == o2::track::PID::Muon) {
    v1.SetM(o2::constants::physics::MassMuon);
  } else {
    LOGF(info, "leptonId supports only Electron or Muon.");
    pair.isOK = false;
    return pair;
  }

  ROOT::Math::PxPyPzMVector v2(pvec1[0], pvec1[1], pvec1[2], o2::constants::physics::MassKaonCharged);
  ROOT::Math::PxPyPzMVector v12 = v1 + v2;
  pair.mass = v12.M();
  pair.pt = v12.Pt();
  pair.p = v12.P();
  pair.deta = v1.Eta() - v2.Eta();                                                                                                           // lepton - hadron
  pair.dphi = RecoDecay::constrainAngle(RecoDecay::constrainAngle(v1.Phi(), 0, 1U) - RecoDecay::constrainAngle(v2.Phi(), 0, 1U), -M_PI, 1U); // lepton - hadron
  pair.isOK = true;
  return pair;
}

template <typename TFitter, typename TCollision, typename TLepton, typename TV0>
LHPair makePairLeptonV0(TFitter& fitter, TCollision const& collision, TLepton const& lepton, TV0 const& v0, o2::track::PID::ID leptonId, o2::track::PID::ID strHadId)
{
  LHPair pair;
  auto leptonParCov = getTrackParCov(lepton);
  leptonParCov.setPID(leptonId);

  const std::array<float, 3> vertex = {collision.posX(), collision.posY(), collision.posZ()};
  const std::array<float, 3> vertexV0 = {v0.x(), v0.y(), v0.z()};
  const std::array<float, 3> momV0 = {v0.px(), v0.py(), v0.pz()};
  std::array<float, 21> covV0 = {0.f};

  constexpr int MomInd[6] = {9, 13, 14, 18, 19, 20}; // cov matrix elements for momentum component
  for (int i = 0; i < 6; i++) {
    covV0[MomInd[i]] = v0.momentumCovMat()[i];
    covV0[i] = v0.positionCovMat()[i];
  }

  auto v0ParCov = o2::track::TrackParCov(vertexV0, momV0, covV0, 0, true);
  v0ParCov.setAbsCharge(0);
  // v0ParCov.setPID(o2::track::PID::Lambda);
  v0ParCov.setPID(strHadId);

  std::array<float, 3> svpos = {0.}; // secondary vertex position
  std::array<float, 3> pvec0 = {0.};
  std::array<float, 3> pvec1 = {0.};

  int nCand = 0;
  try {
    nCand = fitter.process(leptonParCov, v0ParCov);
  } catch (...) {
    LOG(error) << "Exception caught in DCA fitter process call!";
    pair.isOK = false;
    return pair;
  }
  if (nCand == 0) {
    pair.isOK = false;
    return pair;
  }

  fitter.propagateTracksToVertex();
  const auto& vtx = fitter.getPCACandidate();
  for (int i = 0; i < 3; i++) {
    svpos[i] = vtx[i];
  }
  fitter.getTrack(0).getPxPyPzGlo(pvec0); // lepton
  fitter.getTrack(1).getPxPyPzGlo(pvec1); // v0
  std::array<float, 3> pvecSum = {pvec0[0] + pvec1[0], pvec0[1] + pvec1[1], pvec0[2] + pvec1[2]};

  // pair.ptSVL = RecoDecay::sqrtSumOfSquares(pvec0[0], pvec0[1]);
  // pair.ptSVH = RecoDecay::sqrtSumOfSquares(pvec1[0], pvec1[1]);
  // pair.plSVL = pvec0[2];
  // pair.plSVH = pvec1[2];

  pair.cospa = RecoDecay::cpa(vertex, svpos, pvecSum);
  pair.cospaXY = RecoDecay::cpaXY(vertex, svpos, pvecSum);
  pair.cospaRZ = RecoDecay::cpaRZ(vertex, svpos, pvecSum);
  pair.dca2legs = std::sqrt(fitter.getChi2AtPCACandidate());
  pair.chi2PCA = fitter.getChi2AtPCACandidate();
  pair.lxy = std::sqrt(std::pow(svpos[0] - collision.posX(), 2) + std::pow(svpos[1] - collision.posY(), 2));
  pair.lz = svpos[2] - collision.posZ();
  pair.lxyz = std::sqrt(std::pow(svpos[0] - collision.posX(), 2) + std::pow(svpos[1] - collision.posY(), 2) + std::pow(svpos[2] - collision.posZ(), 2));

  auto primaryVertex = getPrimaryVertex(collision);
  std::array<float, 6> covVtxLV0 = fitter.calcPCACovMatrixFlat();
  double phiLV0{}, thetaLV0{};
  getPointDirection(std::array{primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ()}, svpos, phiLV0, thetaLV0);
  pair.lxyzErr = std::sqrt(getRotatedCovMatrixXX(primaryVertex.getCov(), phiLV0, thetaLV0) + getRotatedCovMatrixXX(covVtxLV0, phiLV0, thetaLV0));
  pair.lxyErr = std::sqrt(getRotatedCovMatrixXX(primaryVertex.getCov(), phiLV0, 0.) + getRotatedCovMatrixXX(covVtxLV0, phiLV0, 0.));
  pair.lzErr = std::sqrt(getRotatedCovMatrixXX(primaryVertex.getCov(), 0, thetaLV0) + getRotatedCovMatrixXX(covVtxLV0, 0, thetaLV0));

  // std::array<float, 3> uvFD = {(svpos[0] - collision.posX()) / pair.lxyz, (svpos[1] - collision.posY()) / pair.lxyz, (svpos[2] - collision.posZ()) / pair.lxyz}; // unit vector of flight direction
  // pair.plFD = RecoDecay::dotProd(pvecSum, uvFD);
  // pair.ptFD = RecoDecay::sqrtSumOfSquares(pvecSum[0] - pair.plFD * uvFD[0], pvecSum[1] - pair.plFD * uvFD[1], pvecSum[2] - pair.plFD * uvFD[2]);

  // pair.plFDL = RecoDecay::dotProd(pvec0, uvFD);
  // pair.plFDH = RecoDecay::dotProd(pvec1, uvFD);
  // pair.ptFDL = RecoDecay::sqrtSumOfSquares(pvec0[0] - pair.plFDL * uvFD[0], pvec0[1] - pair.plFDL * uvFD[1], pvec0[2] - pair.plFDL * uvFD[2]);
  // pair.ptFDH = RecoDecay::sqrtSumOfSquares(pvec1[0] - pair.plFDH * uvFD[0], pvec1[1] - pair.plFDH * uvFD[1], pvec1[2] - pair.plFDH * uvFD[2]);

  // // propagate the 2 prongs to the secondary vertex
  // leptonParCov.propagateTo(vtx[0], fitter.getBz());
  // v0ParCov.propagateTo(vtx[0], fitter.getBz());

  // // calculate impact parameter
  // o2::dataformats::DCA dcaLH;
  // auto trackParCovLH = o2::dataformats::V0(fitter.getPCACandidatePos(), pvecSum, fitter.calcPCACovMatrixFlat(), leptonParCov, v0ParCov);
  // trackParCovLH.propagateToDCA(primaryVertex, fitter.getBz(), &dcaLH);

  // pair.impParXY = dcaLH.getY();
  // pair.impParZ = dcaLH.getZ();
  // pair.impParCYY = dcaLH.getSigmaY2();
  // pair.impParCZY = dcaLH.getSigmaYZ();
  // pair.impParCZZ = dcaLH.getSigmaZ2();

  ROOT::Math::PxPyPzMVector v1(pvec0[0], pvec0[1], pvec0[2], o2::constants::physics::MassElectron);
  if (leptonId == o2::track::PID::Electron) {
    v1.SetM(o2::constants::physics::MassElectron);
  } else if (leptonId == o2::track::PID::Muon) {
    v1.SetM(o2::constants::physics::MassMuon);
  } else {
    LOGF(info, "leptonId supports only Electron or Muon.");
    pair.isOK = false;
    return pair;
  }

  ROOT::Math::PxPyPzMVector v2(pvec1[0], pvec1[1], pvec1[2], o2::constants::physics::MassLambda);
  if (strHadId == o2::track::PID::Lambda) {
    v2.SetM(o2::constants::physics::MassLambda);
  } else if (strHadId == o2::track::PID::K0) {
    v2.SetM(o2::constants::physics::MassK0Short);
  } else {
    LOGF(info, "strHadId supports only K0 and Lambda.");
    pair.isOK = false;
    return pair;
  }

  ROOT::Math::PxPyPzMVector v12 = v1 + v2;
  pair.mass = v12.M();
  pair.pt = v12.Pt();
  pair.p = v12.P();
  pair.deta = v1.Eta() - v2.Eta();                                                                                                           // lepton - hadron
  pair.dphi = RecoDecay::constrainAngle(RecoDecay::constrainAngle(v1.Phi(), 0, 1U) - RecoDecay::constrainAngle(v2.Phi(), 0, 1U), -M_PI, 1U); // lepton - hadron
  pair.isOK = true;

  return pair;
}

template <typename TFitter, typename TCollision, typename TLepton, typename TCascade>
LHPair makePairLeptonCascade(TFitter& fitter, TCollision const& collision, TLepton const& lepton, TCascade const& cascade, o2::track::PID::ID leptonId, o2::track::PID::ID strHadId)
{
  LHPair pair;
  auto leptonParCov = getTrackParCov(lepton);
  leptonParCov.setPID(leptonId);

  const std::array<float, 3> vertex = {collision.posX(), collision.posY(), collision.posZ()};
  const std::array<float, 3> vertexCasc = {cascade.x(), cascade.y(), cascade.z()};
  const std::array<float, 3> momCasc = {cascade.px(), cascade.py(), cascade.pz()};

  std::array<float, 21> covCasc = {0.};
  constexpr int MomInd[6] = {9, 13, 14, 18, 19, 20}; // cov matrix elements for momentum component
  for (int i = 0; i < 6; i++) {
    covCasc[MomInd[i]] = cascade.momentumCovMat()[i];
    covCasc[i] = cascade.positionCovMat()[i];
  }

  auto cascParCov = o2::track::TrackParCov(vertexCasc, momCasc, covCasc, cascade.sign(), true);
  cascParCov.setAbsCharge(1);
  cascParCov.setPID(strHadId);
  // if constexpr (cascType == 0) {
  //   cascParCov.setPID(o2::track::PID::XiMinus);
  // } else if constexpr (cascType == 1) {
  //   cascParCov.setPID(o2::track::PID::OmegaMinus);
  // }
  std::array<float, 3> svpos = {0.}; // secondary vertex position
  std::array<float, 3> pvec0 = {0.};
  std::array<float, 3> pvec1 = {0.};

  int nCand = 0;
  try {
    nCand = fitter.process(leptonParCov, cascParCov);
  } catch (...) {
    LOG(error) << "Exception caught in DCA fitter process call!";
    pair.isOK = false;
    return pair;
  }

  if (nCand == 0) {
    pair.isOK = false;
    return pair;
  }

  fitter.propagateTracksToVertex(); // propagate e and Xi/Omega to decay vertex of charm baryon
  const auto& vtx = fitter.getPCACandidate();
  for (int i = 0; i < 3; i++) {
    svpos[i] = vtx[i];
  }
  fitter.getTrack(0).getPxPyPzGlo(pvec0); // lepton
  fitter.getTrack(1).getPxPyPzGlo(pvec1); // cascade
  std::array<float, 3> pvecSum = {pvec0[0] + pvec1[0], pvec0[1] + pvec1[1], pvec0[2] + pvec1[2]};

  // pair.ptSVL = RecoDecay::sqrtSumOfSquares(pvec0[0], pvec0[1]);
  // pair.ptSVH = RecoDecay::sqrtSumOfSquares(pvec1[0], pvec1[1]);
  // pair.plSVL = pvec0[2];
  // pair.plSVH = pvec1[2];

  pair.cospa = RecoDecay::cpa(vertex, svpos, pvecSum);
  pair.cospaXY = RecoDecay::cpaXY(vertex, svpos, pvecSum);
  pair.cospaRZ = RecoDecay::cpaRZ(vertex, svpos, pvecSum);
  pair.dca2legs = std::sqrt(fitter.getChi2AtPCACandidate());
  pair.chi2PCA = fitter.getChi2AtPCACandidate();
  pair.lxy = std::sqrt(std::pow(svpos[0] - collision.posX(), 2) + std::pow(svpos[1] - collision.posY(), 2));
  pair.lz = svpos[2] - collision.posZ();
  pair.lxyz = std::sqrt(std::pow(svpos[0] - collision.posX(), 2) + std::pow(svpos[1] - collision.posY(), 2) + std::pow(svpos[2] - collision.posZ(), 2));

  auto primaryVertex = getPrimaryVertex(collision);
  std::array<float, 6> covVtxLC = fitter.calcPCACovMatrixFlat();
  double phiLC{}, thetaLC{};
  getPointDirection(std::array{primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ()}, svpos, phiLC, thetaLC);
  pair.lxyErr = std::sqrt(getRotatedCovMatrixXX(primaryVertex.getCov(), phiLC, 0.) + getRotatedCovMatrixXX(covVtxLC, phiLC, 0.));
  pair.lzErr = std::sqrt(getRotatedCovMatrixXX(primaryVertex.getCov(), 0, thetaLC) + getRotatedCovMatrixXX(covVtxLC, 0, thetaLC));
  pair.lxyzErr = std::sqrt(getRotatedCovMatrixXX(primaryVertex.getCov(), phiLC, thetaLC) + getRotatedCovMatrixXX(covVtxLC, phiLC, thetaLC));

  // std::array<float, 3> uvFD = {(svpos[0] - collision.posX()) / pair.lxyz, (svpos[1] - collision.posY()) / pair.lxyz, (svpos[2] - collision.posZ()) / pair.lxyz}; // unit vector of flight direction
  // pair.plFD = RecoDecay::dotProd(pvecSum, uvFD);
  // pair.ptFD = RecoDecay::sqrtSumOfSquares(pvecSum[0] - pair.plFD * uvFD[0], pvecSum[1] - pair.plFD * uvFD[1], pvecSum[2] - pair.plFD * uvFD[2]);

  // pair.plFDL = RecoDecay::dotProd(pvec0, uvFD);
  // pair.plFDH = RecoDecay::dotProd(pvec1, uvFD);
  // pair.ptFDL = RecoDecay::sqrtSumOfSquares(pvec0[0] - pair.plFDL * uvFD[0], pvec0[1] - pair.plFDL * uvFD[1], pvec0[2] - pair.plFDL * uvFD[2]);
  // pair.ptFDH = RecoDecay::sqrtSumOfSquares(pvec1[0] - pair.plFDH * uvFD[0], pvec1[1] - pair.plFDH * uvFD[1], pvec1[2] - pair.plFDH * uvFD[2]);

  // // propagate the 2 prongs to the secondary vertex
  // leptonParCov.propagateTo(vtx[0], fitter.getBz());
  // cascParCov.propagateTo(vtx[0], fitter.getBz());

  // // calculate impact parameter
  // o2::dataformats::DCA dcaLH;
  // auto trackParCovLH = o2::dataformats::V0(fitter.getPCACandidatePos(), pvecSum, fitter.calcPCACovMatrixFlat(), leptonParCov, cascParCov);
  // trackParCovLH.propagateToDCA(primaryVertex, fitter.getBz(), &dcaLH);

  // pair.impParXY = dcaLH.getY();
  // pair.impParZ = dcaLH.getZ();
  // pair.impParCYY = dcaLH.getSigmaY2();
  // pair.impParCZY = dcaLH.getSigmaYZ();
  // pair.impParCZZ = dcaLH.getSigmaZ2();

  ROOT::Math::PxPyPzMVector v1(pvec0[0], pvec0[1], pvec0[2], o2::constants::physics::MassElectron);
  if (leptonId == o2::track::PID::Electron) {
    v1.SetM(o2::constants::physics::MassElectron);
  } else if (leptonId == o2::track::PID::Muon) {
    v1.SetM(o2::constants::physics::MassMuon);
  } else {
    LOGF(info, "leptonId supports only Electron or Muon.");
    pair.isOK = false;
    return pair;
  }

  ROOT::Math::PxPyPzMVector v2(pvec1[0], pvec1[1], pvec1[2], o2::constants::physics::MassXiMinus);
  if (strHadId == o2::track::PID::XiMinus) {
    v2.SetM(o2::constants::physics::MassXiMinus);
  } else if (strHadId == o2::track::PID::OmegaMinus) {
    v2.SetM(o2::constants::physics::MassOmegaMinus);
  } else {
    LOGF(info, "strHadId supports only Xi and Omega.");
    pair.isOK = false;
    return pair;
  }

  ROOT::Math::PxPyPzMVector v12 = v1 + v2;
  pair.mass = v12.M();
  pair.pt = v12.Pt();
  pair.p = v12.P();
  pair.deta = v1.Eta() - v2.Eta();                                                                                                           // lepton - hadron
  pair.dphi = RecoDecay::constrainAngle(RecoDecay::constrainAngle(v1.Phi(), 0, 1U) - RecoDecay::constrainAngle(v2.Phi(), 0, 1U), -M_PI, 1U); // lepton - hadron
  pair.isOK = true;

  return pair;
}

} // namespace o2::aod::pwgem::dilepton::utils
#endif // PWGEM_DILEPTON_UTILS_SEMICHARMTAG_H_
