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
/// \file  strangenessBuilderHelper.cxx
/// \since 07/01/2025
/// \brief Utilities to build strange decays in a modular way. Used by strangeness builder task
///

#include "strangenessBuilderHelper.h"
#include "ReconstructionDataFormats/Track.h"
#include "DetectorsBase/GeometryManager.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/Core/trackUtilities.h"

namespace o2
{
namespace pwglf
{
//_______________________________________________________________
strangenessBuilderHelper::strangenessBuilderHelper() {
  // standards hardcoded in builder ...
  // ...but can be changed easily since fitter is public 
  fitter.setPropagateToPCA(true);
  fitter.setMaxR(200.);
  fitter.setMinParamChange(1e-3);
  fitter.setMinRelChi2Change(0.9);
  fitter.setMaxDZIni(1e9);
  fitter.setMaxDXYIni(4.0f);
  fitter.setMaxChi2(1e9);
  fitter.setUseAbsDCA(true);
  fitter.setWeightedFinalPCA(false);

  // LUT has to be loaded later
  lut = nullptr;
  fitter.setMatCorrType(o2::base::Propagator::MatCorrType::USEMatCorrLUT);

  // mag field has to be set later
  fitter.setBz(-999.9f); // will NOT make sense if not changed
}

//_______________________________________________________________
// builds V0 from two tracks. Does not check any conditionals
// except for DCA fitter convergence (should be minimal overhead)
// Resulting properties can be checked with strangenessBuilderHelper::v0
bool strangenessBuilderHelper::buildV0Candidate(
                                                o2::aod::Collision const& collision,
                                                soa::Join<o2::aod::TracksIU, o2::aod::TracksExtra, o2::aod::TracksCovIU>::iterator const& positiveTrack, 
                                                soa::Join<o2::aod::TracksIU, o2::aod::TracksExtra, o2::aod::TracksCovIU>::iterator const& negativeTrack, 
                                                bool useCollinearFit){
  // Calculate DCA with respect to the collision associated to the V0, not individual tracks
  gpu::gpustd::array<float, 2> dcaInfo;

  auto posTrackPar = getTrackPar(positiveTrack);
  o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, posTrackPar, 2.f, fitter.getMatCorrType(), &dcaInfo);
  v0.positiveDCAxy = dcaInfo[0];

  auto negTrackPar = getTrackPar(negativeTrack);
  o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, negTrackPar, 2.f, fitter.getMatCorrType(), &dcaInfo);
  v0.negativeDCAxy = dcaInfo[0];

  o2::track::TrackParCov positiveTrackParam = getTrackParCov(positiveTrack);
  o2::track::TrackParCov negativeTrackParam = getTrackParCov(negativeTrack);

  // Perform DCA fit
  int nCand = 0;
  fitter.setCollinear(useCollinearFit);
  try {
    nCand = fitter.process(positiveTrackParam, negativeTrackParam);
  } catch (...) {
    return false;
  }
  if (nCand == 0) {
    return false;
  }

  v0.positiveTrackX = fitter.getTrack(0).getX();
  v0.negativeTrackX = fitter.getTrack(1).getX();
  positiveTrackParam = fitter.getTrack(0);
  negativeTrackParam = fitter.getTrack(1);
  positiveTrackParam.getPxPyPzGlo(v0.positiveMomentum);
  negativeTrackParam.getPxPyPzGlo(v0.negativeMomentum);
  positiveTrackParam.getXYZGlo(v0.positivePosition);
  negativeTrackParam.getXYZGlo(v0.negativePosition);

  // get decay vertex coordinates
  const auto& vtx = fitter.getPCACandidate();
  for (int i = 0; i < 3; i++) {
    v0.position[i] = vtx[i];
  }

  v0.daughterDCA = TMath::Sqrt(fitter.getChi2AtPCACandidate());
  v0.pointingAngle = TMath::ACos(RecoDecay::cpa(
    std::array{collision.posX(), collision.posY(), collision.posZ()}, 
    std::array{v0.position[0], v0.position[1], v0.position[2]}, 
    std::array{v0.positiveMomentum[0] + v0.negativeMomentum[0], v0.positiveMomentum[1] + v0.negativeMomentum[1], v0.positiveMomentum[2] + v0.negativeMomentum[2]}
  ));

  v0.dcaXY = CalculateDCAStraightToPV(
    v0.position[0], v0.position[1], v0.position[2],
    v0.positiveMomentum[0] + v0.negativeMomentum[0],
    v0.positiveMomentum[1] + v0.negativeMomentum[1],
    v0.positiveMomentum[2] + v0.negativeMomentum[2],
    collision.posX(), collision.posY(), collision.posZ());

  // Calculate masses
  v0.massGamma = RecoDecay::m(std::array{
    std::array{v0.positiveMomentum[0], v0.positiveMomentum[1], v0.positiveMomentum[2]}, 
    std::array{v0.negativeMomentum[0], v0.negativeMomentum[1], v0.negativeMomentum[2]}}, 
    std::array{o2::constants::physics::MassElectron, o2::constants::physics::MassElectron});
  v0.massK0Short = RecoDecay::m(std::array{
    std::array{v0.positiveMomentum[0], v0.positiveMomentum[1], v0.positiveMomentum[2]}, 
    std::array{v0.negativeMomentum[0], v0.negativeMomentum[1], v0.negativeMomentum[2]}}, 
    std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassPionCharged});
  v0.massLambda = RecoDecay::m(std::array{
    std::array{v0.positiveMomentum[0], v0.positiveMomentum[1], v0.positiveMomentum[2]}, 
    std::array{v0.negativeMomentum[0], v0.negativeMomentum[1], v0.negativeMomentum[2]}}, 
    std::array{o2::constants::physics::MassProton, o2::constants::physics::MassPionCharged});
  v0.massAntiLambda = RecoDecay::m(std::array{
    std::array{v0.positiveMomentum[0], v0.positiveMomentum[1], v0.positiveMomentum[2]}, 
    std::array{v0.negativeMomentum[0], v0.negativeMomentum[1], v0.negativeMomentum[2]}}, 
    std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassProton});

  // information validated, V0 built successfully. Signal OK
  return true;   
  }

//_______________________________________________________________
// internal helper to calculate DCAxy of a straight line to a given PV analytically
float strangenessBuilderHelper::CalculateDCAStraightToPV(float X, float Y, float Z, float Px, float Py, float Pz, float pvX, float pvY, float pvZ)
{
  return std::sqrt((std::pow((pvY - Y) * Pz - (pvZ - Z) * Py, 2) + std::pow((pvX - X) * Pz - (pvZ - Z) * Px, 2) + std::pow((pvX - X) * Py - (pvY - Y) * Px, 2)) / (Px * Px + Py * Py + Pz * Pz));
}

//_______________________________________________________________
} // namespace pwglf
} // namespace o2