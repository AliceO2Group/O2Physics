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

namespace o2::aod
{
namespace upgrade_tof
{
DECLARE_SOA_COLUMN(NSigmaElectronInner, nSigmaElectronInner, float); //! NSigma electron InnerTOF
DECLARE_SOA_COLUMN(NSigmaMuonInner, nSigmaMuonInner, float);         //! NSigma muon InnerTOF
DECLARE_SOA_COLUMN(NSigmaPionInner, nSigmaPionInner, float);         //! NSigma pion InnerTOF
DECLARE_SOA_COLUMN(NSigmaKaonInner, nSigmaKaonInner, float);         //! NSigma kaon InnerTOF
DECLARE_SOA_COLUMN(NSigmaProtonInner, nSigmaProtonInner, float);     //! NSigma proton InnerTOF
DECLARE_SOA_COLUMN(NSigmaElectronOuter, nSigmaElectronOuter, float); //! NSigma electron OuterTOF
DECLARE_SOA_COLUMN(NSigmaMuonOuter, nSigmaMuonOuter, float);         //! NSigma muon OuterTOF
DECLARE_SOA_COLUMN(NSigmaPionOuter, nSigmaPionOuter, float);         //! NSigma pion OuterTOF
DECLARE_SOA_COLUMN(NSigmaKaonOuter, nSigmaKaonOuter, float);         //! NSigma kaon OuterTOF
DECLARE_SOA_COLUMN(NSigmaProtonOuter, nSigmaProtonOuter, float);     //! NSigma proton OuterTOF
} // namespace upgrade_tof
DECLARE_SOA_TABLE(UpgradeTof, "AOD", "UPGRADETOF",
                  upgrade_tof::NSigmaElectronInner,
                  upgrade_tof::NSigmaMuonInner,
                  upgrade_tof::NSigmaPionInner,
                  upgrade_tof::NSigmaKaonInner,
                  upgrade_tof::NSigmaProtonInner,
                  upgrade_tof::NSigmaElectronOuter,
                  upgrade_tof::NSigmaMuonOuter,
                  upgrade_tof::NSigmaPionOuter,
                  upgrade_tof::NSigmaKaonOuter,
                  upgrade_tof::NSigmaProtonOuter);
} // namespace o2::aod

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
      const AxisSpec axisVelocity{static_cast<int>(nBinsBeta), 0.0f, +1.1f, "Measured #beta"};
      histos.add("h2dVelocityVsMomentumInner", "h2dVelocityVsMomentumInner", kTH2F, {axisMomentum, axisVelocity});
      histos.add("h2dVelocityVsMomentumOuter", "h2dVelocityVsMomentumOuter", kTH2F, {axisMomentum, axisVelocity});
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
      charge = pdgInfo->Charge();
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

  /// function to calculate track length
  /// \param track the input track
  /// \param x0 the initial position
  /// \param x1 the final position
  /// \param magneticField the magnetic field to use when propagating
  float trackLength(o2::track::TrackParCov track, float x0, float x1, float magneticField)
  {
    std::array<float, 3> pointN;
    std::array<float, 3> pointNplus;
    float length = 0.0;
    track.propagateTo(x0, magneticField);
    for (int iStep = 1; iStep < nStepsLIntegrator; iStep++) {
      track.getXYZGlo(pointN);
      float position = x0 + (x1 - x0) * (static_cast<float>(iStep)) / (static_cast<float>(nStepsLIntegrator - 1));
      track.propagateTo(position, magneticField);
      track.getXYZGlo(pointNplus);
      length += std::hypot(pointNplus[0] - pointN[0], pointNplus[1] - pointN[1], pointNplus[2] - pointN[2]);
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

      float xPv = -100, xInnerTOF = -100, xOuterTOF = -100, trackLengthInnerTOF = -1, trackLengthOuterTOF = -1;
      if (o2track.propagateToDCA(mcPvVtx, dBz))
        xPv = o2track.getX();
      if (!o2track.getXatLabR(innerTOFRadius, xInnerTOF, dBz, o2::track::DirOutward))
        xInnerTOF = -100;
      if (!o2track.getXatLabR(outerTOFRadius, xOuterTOF, dBz, o2::track::DirOutward))
        xOuterTOF = -100;
      if (xPv > -99. && xInnerTOF > -99.)
        trackLengthInnerTOF = trackLength(o2track, xPv, xInnerTOF, dBz);
      if (xPv > -99. && xOuterTOF > -99.)
        trackLengthOuterTOF = trackLength(o2track, xPv, xOuterTOF, dBz);

      // get mass to calculate velocity
      auto pdgInfo = pdg->GetParticle(mcParticle.pdgCode());
      if (pdgInfo == nullptr) {
        continue;
      }
      float expectedTimeInnerTOF = trackLengthInnerTOF / velocity(o2track.getP(), pdgInfo->Mass());
      float expectedTimeOuterTOF = trackLengthOuterTOF / velocity(o2track.getP(), pdgInfo->Mass());

      // Smear with expected resolutions
      float measuredTimeInnerTOF = pRandomNumberGenerator.Gaus(expectedTimeInnerTOF, innerTOFTimeReso);
      float measuredTimeOuterTOF = pRandomNumberGenerator.Gaus(expectedTimeOuterTOF, innerTOFTimeReso);

      // Now we calculate the expected arrival time following certain mass hypotheses
      // and the (imperfect!) reconstructed track parametrizations
      float trackLengthRecoInnerTOF = -1, trackLengthRecoOuterTOF = -1;
      auto recoTrack = getTrackParCov(track);
      if (recoTrack.propagateToDCA(pvVtx, dBz))
        xPv = recoTrack.getX();
      if (!recoTrack.getXatLabR(innerTOFRadius, xInnerTOF, dBz, o2::track::DirOutward))
        xInnerTOF = -100;
      if (!recoTrack.getXatLabR(outerTOFRadius, xOuterTOF, dBz, o2::track::DirOutward))
        xOuterTOF = -100;
      if (xPv > -99. && xInnerTOF > -99.)
        trackLengthRecoInnerTOF = trackLength(recoTrack, xPv, xInnerTOF, dBz);
      if (xPv > -99. && xOuterTOF > -99.)
        trackLengthRecoOuterTOF = trackLength(recoTrack, xPv, xOuterTOF, dBz);

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
        if (trackLengthRecoInnerTOF > 0)
          histos.fill(HIST("h2dVelocityVsMomentumInner"), momentum, innerBeta);
        if (trackLengthRecoOuterTOF > 0)
          histos.fill(HIST("h2dVelocityVsMomentumOuter"), momentum, outerBeta);
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
        if (trackLengthInnerTOF > 0 && trackLengthRecoInnerTOF)
          nSigmaInnerTOF[ii] = deltaTimeInnerTOF[ii] / innerTOFTimeReso;
        if (trackLengthOuterTOF > 0 && trackLengthRecoOuterTOF)
          nSigmaOuterTOF[ii] = deltaTimeOuterTOF[ii] / outerTOFTimeReso;
      }

      // Sigmas have been fully calculated. Please populate the NSigma helper table (once per track)
      upgradeTof(nSigmaInnerTOF[0], nSigmaInnerTOF[1], nSigmaInnerTOF[2], nSigmaInnerTOF[3], nSigmaInnerTOF[4],
                 nSigmaOuterTOF[0], nSigmaOuterTOF[1], nSigmaOuterTOF[2], nSigmaOuterTOF[3], nSigmaOuterTOF[4]);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<OnTheFlyTOFPID>(cfgc)};
}
