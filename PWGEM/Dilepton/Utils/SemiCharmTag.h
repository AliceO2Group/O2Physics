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
#include <ReconstructionDataFormats/PID.h>
#include <ReconstructionDataFormats/Track.h>
#include <ReconstructionDataFormats/TrackParametrizationWithError.h>

#include <Math/Vector4D.h> // IWYU pragma: keep (do not replace with Math/Vector4Dfwd.h)
#include <Math/Vector4Dfwd.h>

#include <array>

namespace o2::aod::pwgem::dilepton::utils
{

struct LHPair { // struct to store electron-hadron pair information
  float mass{-999.f};
  float dca2legs{-999.f};
  float cospa{-999.f};
  float lxy{-999.f};
  float lz{-999.f};
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
  fitter.getTrack(0).getPxPyPzGlo(pvec0); // lepton
  fitter.getTrack(1).getPxPyPzGlo(pvec1); // track
  std::array<float, 3> pvecSum = {pvec0[0] + pvec1[0], pvec0[1] + pvec1[1], pvec0[2] + pvec1[2]};

  float cospa = RecoDecay::cpa(vertex, svpos, pvecSum);
  float dca2legs = std::sqrt(fitter.getChi2AtPCACandidate());
  float lxy = std::sqrt(std::pow(svpos[0] - collision.posX(), 2) + std::pow(svpos[1] - collision.posY(), 2));
  float lz = std::fabs(svpos[2] - collision.posZ());

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
  if (strHadId == o2::track::PID::Kaon) {
    v2.SetM(o2::constants::physics::MassKaonCharged);
  } else {
    LOGF(info, "strHadId supports only Kaon.");
    pair.isOK = false;
    return pair;
  }

  ROOT::Math::PxPyPzMVector v12 = v1 + v2;

  pair.mass = v12.M();
  pair.dca2legs = dca2legs;
  pair.cospa = cospa;
  pair.lxy = lxy;
  pair.lz = lz;
  pair.isOK = true;

  return pair;
}

template <typename TFitter, typename TCollision, typename TLepton, typename TV0>
LHPair makePairLeptonV0(TFitter& fitter, TCollision const& collision, TLepton const& lepton, TV0 const& v0, o2::track::PID::ID leptonId, o2::track::PID::ID strHadId)
{
  LHPair pair;
  auto trackParCov = getTrackParCov(lepton);
  trackParCov.setPID(leptonId);

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
    nCand = fitter.process(trackParCov, v0ParCov);
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
  fitter.getTrack(0).getPxPyPzGlo(pvec0); // electron
  fitter.getTrack(1).getPxPyPzGlo(pvec1); // v0
  std::array<float, 3> pvecSum = {pvec0[0] + pvec1[0], pvec0[1] + pvec1[1], pvec0[2] + pvec1[2]};

  float cospa = RecoDecay::cpa(vertex, svpos, pvecSum);
  float dca2legs = std::sqrt(fitter.getChi2AtPCACandidate());
  float lxy = std::sqrt(std::pow(svpos[0] - collision.posX(), 2) + std::pow(svpos[1] - collision.posY(), 2));
  float lz = std::fabs(svpos[2] - collision.posZ());

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
  pair.dca2legs = dca2legs;
  pair.cospa = cospa;
  pair.lxy = lxy;
  pair.lz = lz;
  pair.isOK = true;

  return pair;
}

template <typename TFitter, typename TCollision, typename TLepton, typename TCascade>
LHPair makePairLeptonCascade(TFitter& fitter, TCollision const& collision, TLepton const& lepton, TCascade const& cascade, o2::track::PID::ID leptonId, o2::track::PID::ID strHadId)
{
  LHPair pair;
  auto trackParCov = getTrackParCov(lepton);
  trackParCov.setPID(leptonId);

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
    nCand = fitter.process(trackParCov, cascParCov);
  } catch (...) {
    LOG(error) << "Exception caught in DCA fitter process call!";
    return pair;
  }
  if (nCand == 0) {
    return pair;
  }

  fitter.propagateTracksToVertex(); // propagate e and Xi/Omega to decay vertex of charm baryon
  const auto& vtx = fitter.getPCACandidate();
  for (int i = 0; i < 3; i++) {
    svpos[i] = vtx[i];
  }
  fitter.getTrack(0).getPxPyPzGlo(pvec0); // electron
  fitter.getTrack(1).getPxPyPzGlo(pvec1); // v0
  std::array<float, 3> pvecSum = {pvec0[0] + pvec1[0], pvec0[1] + pvec1[1], pvec0[2] + pvec1[2]};

  float cospa = RecoDecay::cpa(vertex, svpos, pvecSum);
  float dca2legs = std::sqrt(fitter.getChi2AtPCACandidate());
  float lxy = std::sqrt(std::pow(svpos[0] - collision.posX(), 2) + std::pow(svpos[1] - collision.posY(), 2));
  float lz = std::fabs(svpos[2] - collision.posZ());
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
  pair.dca2legs = dca2legs;
  pair.cospa = cospa;
  pair.lxy = lxy;
  pair.lz = lz;
  pair.isOK = true;

  return pair;
}

} // namespace o2::aod::pwgem::dilepton::utils
#endif // PWGEM_DILEPTON_UTILS_SEMICHARMTAG_H_
