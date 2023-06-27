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

//
// Task to add a table of track parameters propagated to the primary vertex
//

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/ASoAHelpers.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/trackUtilities.h"
#include "ReconstructionDataFormats/DCA.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "CommonUtils/NameConf.h"
#include "CCDB/CcdbApi.h"
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsCalibration/MeanVertexObject.h"
#include "CommonConstants/GeomConstants.h"
#include "CommonConstants/PhysicsConstants.h"
#include "TRandom3.h"
#include "ALICE3/DataModel/OTFTOF.h"
#include "DetectorsVertexing/HelixHelper.h"

/// \file onTheFlyTOFPID.cxx
///
/// \brief This task goes straight from a combination of track table and mcParticles
/// and a custom TOF configuration to a table of TOF NSigmas for the particles
/// being analysed. It currently contemplates 5 particle types:
/// electrons, pions, kaons, protons and muons
///
/// More particles could be added but would have to be added to the LUT
/// being used in the onTheFly tracker task.
///
/// \author David Dobrigkeit Chinellato, UNICAMP

using namespace o2;
using namespace o2::framework;

struct OnTheFlyTOFPID {
  Produces<aod::UpgradeTof> upgradeTof;

  // necessary for particle charges
  Service<O2DatabasePDG> pdg;

  // these are the settings governing the TOF layers to be used
  // note that there are two layers foreseen for now: inner and outer TOF
  // more could be added (especially a disk TOF at a certain z?)
  // in the evolution of this effort
  Configurable<float> dBz{"dBz", 20, "magnetic field (kilogauss)"};
  Configurable<float> innerTOFRadius{"innerTOFRadius", 20, "barrel inner TOF radius (cm)"};
  Configurable<float> outerTOFRadius{"outerTOFRadius", 80, "barrel outer TOF radius (cm)"};
  Configurable<float> innerTOFTimeReso{"innerTOFTimeReso", 20, "barrel inner TOF time error (ps)"};
  Configurable<float> outerTOFTimeReso{"outerTOFTimeReso", 20, "barrel outer TOF time error (ps)"};
  Configurable<int> nStepsLIntegrator{"nStepsLIntegrator", 200, "number of steps in length integrator"};
  Configurable<bool> doQAplots{"doQAplots", true, "do basic velocity plot qa"};
  Configurable<int> nBinsBeta{"nBinsBeta", 2200, "number of bins in beta"};
  Configurable<int> nBinsP{"nBinsP", 80, "number of bins in momentum"};
  Configurable<int> nBinsTrackLengthInner{"nBinsTrackLengthInner", 300, "number of bins in track length"};
  Configurable<int> nBinsTrackLengthOuter{"nBinsTrackLengthOuter", 300, "number of bins in track length"};
  Configurable<int> nBinsTrackDeltaLength{"nBinsTrackDeltaLength", 100, "number of bins in delta track length"};

  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;

  // needed: random number generator for smearing
  TRandom3 pRandomNumberGenerator;

  // for handling basic QA histograms if requested
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext& initContext)
  {
    pRandomNumberGenerator.SetSeed(0); // fully randomize

    if (doQAplots) {
      const AxisSpec axisMomentum{static_cast<int>(nBinsP), 0.0f, +4.0f, "#it{p} (GeV/#it{c})"};
      const AxisSpec axisMomentumSmall{static_cast<int>(nBinsP), 0.0f, +1.0f, "#it{p} (GeV/#it{c})"};
      const AxisSpec axisVelocity{static_cast<int>(nBinsBeta), 0.0f, +1.1f, "Measured #beta"};
      const AxisSpec axisTrackLengthInner{static_cast<int>(nBinsTrackLengthInner), 0.0f, 60.0f, "Track length (cm)"};
      const AxisSpec axisTrackLengthOuter{static_cast<int>(nBinsTrackLengthOuter), 0.0f, 300.0f, "Track length (cm)"};
      const AxisSpec axisTrackDeltaLength{static_cast<int>(nBinsTrackDeltaLength), 0.0f, 30.0f, "Delta Track length (cm)"};
      histos.add("h2dVelocityVsMomentumInner", "h2dVelocityVsMomentumInner", kTH2F, {axisMomentum, axisVelocity});
      histos.add("h2dVelocityVsMomentumOuter", "h2dVelocityVsMomentumOuter", kTH2F, {axisMomentum, axisVelocity});
      histos.add("h2dTrackLengthInnerVsPt", "h2dTrackLengthInnerVsPt", kTH2F, {axisMomentumSmall, axisTrackLengthInner});
      histos.add("h2dTrackLengthOuterVsPt", "h2dTrackLengthOuterVsPt", kTH2F, {axisMomentumSmall, axisTrackLengthOuter});

      histos.add("h2dTrackLengthInnerRecoVsPt", "h2dTrackLengthInnerRecoVsPt", kTH2F, {axisMomentumSmall, axisTrackLengthInner});
      histos.add("h2dTrackLengthOuterRecoVsPt", "h2dTrackLengthOuterRecoVsPt", kTH2F, {axisMomentumSmall, axisTrackLengthOuter});

      histos.add("h2dDeltaTrackLengthInnerVsPt", "h2dDeltaTrackLengthInnerVsPt", kTH2F, {axisMomentumSmall, axisTrackDeltaLength});
      histos.add("h2dDeltaTrackLengthOuterVsPt", "h2dDeltaTrackLengthOuterVsPt", kTH2F, {axisMomentumSmall, axisTrackDeltaLength});
    }
  }

  /// Function to convert a McParticle into a perfect Track
  /// \param particle the particle to convert (mcParticle)
  /// \param o2track the address of the resulting TrackParCov
  template <typename McParticleType>
  void convertMCParticleToO2Track(McParticleType& particle, o2::track::TrackParCov& o2track)
  {
    // FIXME: this is a fundamentally important piece of code.
    // It could be placed in a utility file instead of here.
    auto pdgInfo = pdg->GetParticle(particle.pdgCode());
    int charge = 0;
    if (pdgInfo != nullptr) {
      charge = pdgInfo->Charge() / 3;
    }
    std::array<float, 5> params;
    std::array<float, 15> covm = {0.};
    float s, c, x;
    o2::math_utils::sincos(particle.phi(), s, c);
    o2::math_utils::rotateZInv(particle.vx(), particle.vy(), x, params[0], s, c);
    params[1] = particle.vz();
    params[2] = 0.; // since alpha = phi
    auto theta = 2. * std::atan(std::exp(-particle.eta()));
    params[3] = 1. / std::tan(theta);
    params[4] = charge / particle.pt();

    // Initialize TrackParCov in-place
    new (&o2track)(o2::track::TrackParCov)(x, particle.phi(), params, covm);
  }

  /// function to calculate track length of this track up to a certain radius
  /// \param track the input track
  /// \param radius the radius of the layer you're calculating the length to
  /// \param magneticField the magnetic field to use when propagating
  float trackLength(o2::track::TrackParCov track, float radius, float magneticField)
  {
    // don't make use of the track parametrization
    float length = -100;

    o2::math_utils::CircleXYf_t trcCircle;
    float sna, csa;
    track.getCircleParams(magneticField, trcCircle, sna, csa);

    // distance between circle centers (one circle is at origin -> easy)
    float centerDistance = std::hypot(trcCircle.xC, trcCircle.yC);

    // condition of circles touching - if not satisfied returned length will be -100
    if (centerDistance < trcCircle.rC + radius && centerDistance > fabs(trcCircle.rC - radius)) {
      length = 0.0f;

      // base radical direction
      float ux = trcCircle.xC / centerDistance;
      float uy = trcCircle.yC / centerDistance;
      // calculate perpendicular vector (normalized) for +/- displacement
      float vx = -uy;
      float vy = +ux;
      // calculate coordinate for radical line
      float radical = (centerDistance * centerDistance - trcCircle.rC * trcCircle.rC + radius * radius) / (2.0f * centerDistance);
      // calculate absolute displacement from center-to-center axis
      float displace = (0.5f / centerDistance) * TMath::Sqrt(
                                                   (-centerDistance + trcCircle.rC - radius) *
                                                   (-centerDistance - trcCircle.rC + radius) *
                                                   (-centerDistance + trcCircle.rC + radius) *
                                                   (centerDistance + trcCircle.rC + radius));

      // possible intercept points of track and TOF layer in 2D plane
      float point1[2] = {radical * ux + displace * vx, radical * uy + displace * vy};
      float point2[2] = {radical * ux - displace * vx, radical * uy - displace * vy};

      // decide on correct intercept point
      std::array<float, 3> mom;
      track.getPxPyPzGlo(mom);
      float scalarProduct1 = point1[0] * mom[0] + point1[1] * mom[1];
      float scalarProduct2 = point2[0] * mom[0] + point2[1] * mom[1];

      // get start point
      std::array<float, 3> startPoint;
      track.getXYZGlo(startPoint);

      float cosAngle = -1000, modulus = -1000;

      if (scalarProduct1 > scalarProduct2) {
        modulus = std::hypot(point1[0] - trcCircle.xC, point1[1] - trcCircle.yC) * std::hypot(startPoint[0] - trcCircle.xC, startPoint[1] - trcCircle.yC);
        cosAngle = (point1[0] - trcCircle.xC) * (startPoint[0] - trcCircle.xC) + (point1[1] - trcCircle.yC) * (startPoint[0] - trcCircle.yC);
      } else {
        modulus = std::hypot(point2[0] - trcCircle.xC, point2[1] - trcCircle.yC) * std::hypot(startPoint[0] - trcCircle.xC, startPoint[1] - trcCircle.yC);
        cosAngle = (point2[0] - trcCircle.xC) * (startPoint[0] - trcCircle.xC) + (point2[1] - trcCircle.yC) * (startPoint[0] - trcCircle.yC);
      }
      cosAngle /= modulus;
      length = trcCircle.rC * TMath::ACos(cosAngle);
      length *= sqrt(1.0f + track.getTgl() * track.getTgl());
    }
    return length;
  }

  /// returns velocity in centimeters per picoseconds
  /// \param momentum the momentum of the tarck
  /// \param mass the mass of the particle
  float velocity(float momentum, float mass)
  {
    float a = std::pow(momentum / mass, 2);
    // uses light speed in cm/ps so output is in those units
    return (o2::constants::physics::LightSpeedCm2NS / 1e+3) * std::sqrt(a / (1 + a));
  }

  void process(soa::Join<aod::Collisions, aod::McCollisionLabels>::iterator const& collision, soa::Join<aod::Tracks, aod::TracksCov, aod::McTrackLabels> const& tracks, aod::McParticles const&, aod::McCollisions const&)
  {
    o2::dataformats::VertexBase pvVtx({collision.posX(), collision.posY(), collision.posZ()},
                                      {collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ()});

    std::array<float, 6> mcPvCov = {0.};
    o2::dataformats::VertexBase mcPvVtx({0.0f, 0.0f, 0.0f}, mcPvCov);
    if (collision.has_mcCollision()) {
      auto mcCollision = collision.mcCollision();
      mcPvVtx.setX(mcCollision.posX());
      mcPvVtx.setY(mcCollision.posY());
      mcPvVtx.setZ(mcCollision.posZ());
    } // else remains untreated for now

    for (const auto& track : tracks) {
      // first step: find precise arrival time (if any)
      // --- convert track into perfect track
      if (!track.has_mcParticle()) // should always be OK but check please
        continue;

      o2::track::TrackParCov o2track;
      auto mcParticle = track.mcParticle();
      convertMCParticleToO2Track(mcParticle, o2track);

      float xPv = -100, trackLengthInnerTOF = -1, trackLengthOuterTOF = -1;
      if (o2track.propagateToDCA(mcPvVtx, dBz))
        xPv = o2track.getX();
      if (xPv > -99.) {
        trackLengthInnerTOF = trackLength(o2track, innerTOFRadius, dBz);
        trackLengthOuterTOF = trackLength(o2track, outerTOFRadius, dBz);
      }

      // get mass to calculate velocity
      auto pdgInfo = pdg->GetParticle(mcParticle.pdgCode());
      if (pdgInfo == nullptr) {
        continue;
      }
      float expectedTimeInnerTOF = trackLengthInnerTOF / velocity(o2track.getP(), pdgInfo->Mass());
      float expectedTimeOuterTOF = trackLengthOuterTOF / velocity(o2track.getP(), pdgInfo->Mass());

      // Smear with expected resolutions
      float measuredTimeInnerTOF = pRandomNumberGenerator.Gaus(expectedTimeInnerTOF, innerTOFTimeReso);
      float measuredTimeOuterTOF = pRandomNumberGenerator.Gaus(expectedTimeOuterTOF, outerTOFTimeReso);

      // Now we calculate the expected arrival time following certain mass hypotheses
      // and the (imperfect!) reconstructed track parametrizations
      float trackLengthRecoInnerTOF = -1, trackLengthRecoOuterTOF = -1;
      auto recoTrack = getTrackParCov(track);
      if (recoTrack.propagateToDCA(pvVtx, dBz))
        xPv = recoTrack.getX();
      if (xPv > -99.) {
        trackLengthRecoInnerTOF = trackLength(recoTrack, innerTOFRadius, dBz);
        trackLengthRecoOuterTOF = trackLength(recoTrack, outerTOFRadius, dBz);
      }

      // Straight to Nsigma
      float deltaTimeInnerTOF[5], nSigmaInnerTOF[5];
      float deltaTimeOuterTOF[5], nSigmaOuterTOF[5];
      int lpdg_array[5] = {kElectron, kMuonMinus, kPiPlus, kKPlus, kProton};
      float masses[5];

      if (doQAplots) {
        float momentum = recoTrack.getP();
        // unit conversion: length in cm, time in ps
        float innerBeta = 1e+3 * (trackLengthInnerTOF / measuredTimeInnerTOF) / o2::constants::physics::LightSpeedCm2NS;
        float outerBeta = 1e+3 * (trackLengthOuterTOF / measuredTimeOuterTOF) / o2::constants::physics::LightSpeedCm2NS;
        if (trackLengthRecoInnerTOF > 0) {
          histos.fill(HIST("h2dVelocityVsMomentumInner"), momentum, innerBeta);
          histos.fill(HIST("h2dTrackLengthInnerVsPt"), o2track.getPt(), trackLengthInnerTOF);
          histos.fill(HIST("h2dTrackLengthInnerRecoVsPt"), o2track.getPt(), trackLengthRecoInnerTOF);
        }
        if (trackLengthRecoOuterTOF > 0) {
          histos.fill(HIST("h2dVelocityVsMomentumOuter"), momentum, outerBeta);
          histos.fill(HIST("h2dTrackLengthOuterVsPt"), o2track.getPt(), trackLengthOuterTOF);
          histos.fill(HIST("h2dTrackLengthOuterRecoVsPt"), o2track.getPt(), trackLengthRecoOuterTOF);
        }
      }

      for (int ii = 0; ii < 5; ii++) {
        nSigmaInnerTOF[ii] = -100;
        nSigmaOuterTOF[ii] = -100;

        auto pdgInfoThis = pdg->GetParticle(lpdg_array[ii]);
        masses[ii] = pdgInfoThis->Mass();
        deltaTimeInnerTOF[ii] = trackLengthRecoInnerTOF / velocity(recoTrack.getP(), masses[ii]) - measuredTimeInnerTOF;
        deltaTimeOuterTOF[ii] = trackLengthRecoOuterTOF / velocity(recoTrack.getP(), masses[ii]) - measuredTimeOuterTOF;

        // Fixme: assumes dominant resolution effect is the TOF resolution
        // and not the tracking itself. It's *probably* a fair assumption
        // but it should be tested further!
        if (trackLengthInnerTOF > 0 && trackLengthRecoInnerTOF > 0)
          nSigmaInnerTOF[ii] = deltaTimeInnerTOF[ii] / innerTOFTimeReso;
        if (trackLengthOuterTOF > 0 && trackLengthRecoOuterTOF > 0)
          nSigmaOuterTOF[ii] = deltaTimeOuterTOF[ii] / outerTOFTimeReso;
      }

      float deltaTrackLengthInnerTOF = abs(trackLengthInnerTOF - trackLengthRecoInnerTOF);
      if (trackLengthInnerTOF > 0 && trackLengthRecoInnerTOF > 0) {
        histos.fill(HIST("h2dDeltaTrackLengthInnerVsPt"), o2track.getPt(), deltaTrackLengthInnerTOF);
      }
      float deltaTrackLengthOuterTOF = abs(trackLengthOuterTOF - trackLengthRecoOuterTOF);
      if (trackLengthOuterTOF > 0 && trackLengthRecoOuterTOF > 0) {
        histos.fill(HIST("h2dDeltaTrackLengthOuterVsPt"), o2track.getPt(), deltaTrackLengthInnerTOF);
      }

      // Sigmas have been fully calculated. Please populate the NSigma helper table (once per track)
      upgradeTof(nSigmaInnerTOF[0], nSigmaInnerTOF[1], nSigmaInnerTOF[2], nSigmaInnerTOF[3], nSigmaInnerTOF[4], trackLengthInnerTOF, trackLengthRecoInnerTOF, deltaTrackLengthInnerTOF,
                 nSigmaOuterTOF[0], nSigmaOuterTOF[1], nSigmaOuterTOF[2], nSigmaOuterTOF[3], nSigmaOuterTOF[4], trackLengthOuterTOF, trackLengthRecoOuterTOF, deltaTrackLengthOuterTOF);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<OnTheFlyTOFPID>(cfgc)};
}
