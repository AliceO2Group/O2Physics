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


// This task goes straight from a combination of track table and mcParticles
// and a custom TOF configuration to a table of TOF NSigmas for the particles
// being analysed. It currently contemplates 5 particle types:
// electrons, pions, kaons, protons and muons
//
// More particles could be added but would have to be added to the LUT
// being used in the onTheFly tracker task.

using namespace o2;
using namespace o2::framework;

namespace o2::aod
{
namespace upgrade_tof{
DECLARE_SOA_COLUMN(NSigmaElectronInner, nSigmaElectronInner, float); //! NSigma electron iTOF
DECLARE_SOA_COLUMN(NSigmaMuonInner, nSigmaMuonInner, float);         //! NSigma muon iTOF
DECLARE_SOA_COLUMN(NSigmaPionInner, nSigmaPionInner, float);         //! NSigma pion iTOF
DECLARE_SOA_COLUMN(NSigmaKaonInner, nSigmaKaonInner, float);         //! NSigma kaon iTOF
DECLARE_SOA_COLUMN(NSigmaProtonInner, nSigmaProtonInner, float);     //! NSigma proton iTOF
DECLARE_SOA_COLUMN(NSigmaElectronOuter, nSigmaElectronOuter, float); //! NSigma electron oTOF
DECLARE_SOA_COLUMN(NSigmaMuonOuter, nSigmaMuonOuter, float);         //! NSigma muon oTOF
DECLARE_SOA_COLUMN(NSigmaPionOuter, nSigmaPionOuter, float);         //! NSigma pion oTOF
DECLARE_SOA_COLUMN(NSigmaKaonOuter, nSigmaKaonOuter, float);         //! NSigma kaon oTOF
DECLARE_SOA_COLUMN(NSigmaProtonOuter, nSigmaProtonOuter, float);     //! NSigma proton oTOF
} // namespace upgradetof
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
                  upgrade_tof::NSigmaProtonOuter
                  );
} // namespace o2::aod

struct OnTheFlyTOFPID {
  Produces<aod::UpgradeTof> upgradetof;

  // necessary for particle charges
  Service<O2DatabasePDG> pdg;

  // these are the settings governing the TOF layers to be used
  // note that there are two layers foreseen for now: inner and outer TOF
  // more could be added (especially a disk TOF at a certain z?)
  // in the evolution of this effort
  Configurable<float> dBz{"dBz", 20, "magnetic field (kilogauss)"};
  Configurable<float> iTOFRadius{"iTOFRadius", 20, "barrel inner TOF radius (cm)"};
  Configurable<float> oTOFRadius{"oTOFRadius", 80, "barrel outer TOF radius (cm)"};
  Configurable<float> iTOFTimeReso{"iTOFTimeReso", 20, "barrel inner TOF time error (ps)"};
  Configurable<float> oTOFTimeReso{"oTOFTimeReso", 20, "barrel outer TOF time error (ps)"};
  Configurable<int> lNStepsLIntegrator{"lNStepsLIntegrator", 200, "number of steps in length integrator"};
  Configurable<bool> doQAplots{"doQAplots", true, "do basic velocity plot qa"};

  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;

  // needed: random number generator for smearing
  TRandom3 lPRNG;

  // for handling basic QA histograms if requested
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext& initContext)
  {
    lPRNG.SetSeed(0); // fully randomize

    if(doQAplots){
      const AxisSpec axisMomentum{static_cast<int>(80), 0.0f, +4.0f, "#it{p} (GeV/#it{c})"};
      const AxisSpec axisVelocity{static_cast<int>(110), 0.0f, +1.1f, "Measured #beta"};
      histos.add("h2dVelocityVsMomentumInner", "h2dVelocityVsMomentumInner", kTH2F, {axisMomentum, axisVelocity});
      histos.add("h2dVelocityVsMomentumOuter", "h2dVelocityVsMomentumOuter", kTH2F, {axisMomentum, axisVelocity});
    }
  }

  template <typename mcParticleType>
  void convertMCParticleToO2Track(mcParticleType& particle, o2::track::TrackParCov& o2track)
  {
    // FIXME: this is a fundamentally important piece of code.
    // It could be placed in a utility file instead of here.
    std::array<float, 3> xyz = {static_cast<float>(particle.vx()), static_cast<float>(particle.vy()), static_cast<float>(particle.vz())};
    std::array<float, 3> ptetaphi = {static_cast<float>(particle.pt()), static_cast<float>(particle.eta()), static_cast<float>(particle.phi())};
    auto pdgInfo = pdg->GetParticle(particle.pdgCode());
    int charge = 0;
    if (pdgInfo != nullptr)
      charge = pdgInfo->Charge();
    std::array<float, 5> params;
    std::array<float, 15> covm = {0.};
    float s, c, x;
    o2::math_utils::sincos(ptetaphi[2], s, c);
    o2::math_utils::rotateZInv(xyz[0], xyz[1], x, params[0], s, c);
    params[1] = xyz[2];
    params[2] = 0.; // since alpha = phi
    auto theta = 2. * std::atan(std::exp(-ptetaphi[1]));
    params[3] = 1. / std::tan(theta);
    params[4] = charge / ptetaphi[0];

    // Initialize TrackParCov in-place
    new (&o2track)(o2::track::TrackParCov)(x, ptetaphi[2], params, covm);
  }

  float trackLength(o2::track::TrackParCov track, float lX0, float lX1, float lMagneticField)
  {
    // utility class to calculate track length from track position X0 to track position X1 
    // with a given magnetic field
    std::array<float, 3> lPointN;
    std::array<float, 3> lPointNplus;
    float lLength = 0.0;
    track.propagateTo(lX0, lMagneticField);
    for (int iStep = 1; iStep < lNStepsLIntegrator; iStep++) {
      track.getXYZGlo(lPointN);
      float lPosition = lX0 + (lX1 - lX0) * ((float)(iStep)) / ((float)(lNStepsLIntegrator - 1));
      track.propagateTo(lPosition, lMagneticField);
      track.getXYZGlo(lPointNplus);
      lLength += std::hypot(lPointNplus[0] - lPointN[0], lPointNplus[1] - lPointN[1], lPointNplus[2] - lPointN[2]);
    }
    return lLength;
  }

  float Velocity(float lMomentum, float lMass){
    //Momentum p and mass m -> returns speed in centimeters per picosecond
    //Useful for TOF calculations
    float lA = std::pow(lMomentum / lMass, 2);
    return (o2::constants::physics::LightSpeedCm2NS / 1e+3)*TMath::Sqrt(lA/(1+lA));
  }

  void process(soa::Join<aod::Collisions, aod::McCollisionLabels>::iterator const& collision, soa::Join<aod::Tracks, aod::TracksCov, aod::McTrackLabels> const& tracks, aod::McParticles const&)
  {
    o2::math_utils::Point3D<float> pvPos{collision.posX(),collision.posY(),collision.posZ()};
    std::array<float, 6> pvCov;
    pvCov[0] = collision.covXX();
    pvCov[1] = collision.covXY();
    pvCov[2] = collision.covYY();
    pvCov[3] = collision.covXZ();
    pvCov[4] = collision.covXY();
    pvCov[5] = collision.covXZ();
    o2::dataformats::VertexBase pvvtx(pvPos, pvCov);

    auto mcCollision = collision.mcCollision();
    o2::math_utils::Point3D<float> mcPvPos{mcCollision.posX(),mcCollision.posY(),mcCollision.posZ()};
    std::array<float, 6> mcPvCov; // dummy! 
    for(int ii=0; ii<6; ii++) mcPvCov[ii] = 1e-6;
    o2::dataformats::VertexBase mcpvvtx(mcPvPos, mcPvCov);

    for (const auto& track : tracks) {
      // first step: find precise arrival time (if any)
      // --- convert track into perfect track
      if (!track.has_mcParticle()) // should always be OK but check please
        LOG(error) << "Oh no! No mcParticle label for this track! This shouldn't happen!";

      o2::track::TrackParCov o2track;
      auto mcParticle = track.mcParticle();
      convertMCParticleToO2Track(mcParticle, o2track);

      float lX_PV = -100, lX_iTOF = -100, lX_oTOF = -100, lThisTrackLength_iTOF = -1, lThisTrackLength_oTOF = -1;
      if (o2track.propagateToDCA(mcpvvtx, dBz))
        lX_PV = o2track.getX();
      if (!o2track.getXatLabR(iTOFRadius, lX_iTOF, dBz, o2::track::DirOutward))
        lX_iTOF = -100;
      if (!o2track.getXatLabR(oTOFRadius, lX_oTOF, dBz, o2::track::DirOutward))
        lX_oTOF = -100;
      if (lX_PV > -99. && lX_iTOF > -99.)
        lThisTrackLength_iTOF = trackLength(o2track, lX_PV, lX_iTOF, dBz);
      if (lX_PV > -99. && lX_oTOF > -99.)
        lThisTrackLength_oTOF = trackLength(o2track, lX_PV, lX_oTOF, dBz);

      // get mass to calculate velocity
      auto pdgInfo = pdg->GetParticle(mcParticle.pdgCode());
      float lExpectedTime_iTOF = lThisTrackLength_iTOF / Velocity(o2track.getP(), pdgInfo->Mass());
      float lExpectedTime_oTOF = lThisTrackLength_oTOF / Velocity(o2track.getP(), pdgInfo->Mass());

      // Smear with expected resolutions
      float lMeasuredTime_iTOF = lPRNG.Gaus(lExpectedTime_iTOF, iTOFTimeReso);
      float lMeasuredTime_oTOF = lPRNG.Gaus(lExpectedTime_oTOF, iTOFTimeReso);

      // Now we calculate the expected arrival time following certain mass hypotheses
      // and the (imperfect!) reconstructed track parametrizations
      float lThisTrackLengthReco_iTOF = -1, lThisTrackLengthReco_oTOF = -1;
      auto recoTrack = getTrackParCov(track);
      if (recoTrack.propagateToDCA(pvvtx, dBz))
        lX_PV = recoTrack.getX();
      if (!recoTrack.getXatLabR(iTOFRadius, lX_iTOF, dBz, o2::track::DirOutward))
        lX_iTOF = -100;
      if (!recoTrack.getXatLabR(oTOFRadius, lX_oTOF, dBz, o2::track::DirOutward))
        lX_oTOF = -100;
      if (lX_PV > -99. && lX_iTOF > -99.)
        lThisTrackLengthReco_iTOF = TrackLength(recoTrack, lX_PV, lX_iTOF, dBz);
      if (lX_PV > -99. && lX_oTOF > -99.)
        lThisTrackLengthReco_oTOF = TrackLength(recoTrack, lX_PV, lX_oTOF, dBz);

      // Straight to Nsigma
      float lDeltaTime_iTOF[5], lNSigma_iTOF[5];
      float lDeltaTime_oTOF[5], lNSigma_oTOF[5];
      int lpdg_array[5] = {11, 13, 211, 321, 2212};
      float lMasses[5];

      if (doQAplots) {
        float momentum = recoTrack.getP();
        float innerbeta = (lThisTrackLength_iTOF / (1e+3 * lMeasuredTime_iTOF)) / o2::constants::physics::LightSpeedCm2NS;
        float outerbeta = (lThisTrackLength_oTOF / (1e+3 * lMeasuredTime_oTOF)) / o2::constants::physics::LightSpeedCm2NS;

        histos.fill(HIST("h2dVelocityVsMomentumInner"), momentum, innerbeta);
        histos.fill(HIST("h2dVelocityVsMomentumOuter"), momentum, outerbeta);
      }

      for (int ii = 0; ii < 5; ii++) {
        lNSigma_iTOF[ii] = -100;
        lNSigma_oTOF[ii] = -100;

        auto pdgInfoThis = pdg->GetParticle(lpdg_array[ii]);
        lMasses[ii] = pdgInfoThis->Mass();
        lDeltaTime_iTOF[ii] = lThisTrackLengthReco_iTOF / Velocity(recoTrack.getP(), lMasses[ii]) - lMeasuredTime_iTOF;
        lDeltaTime_oTOF[ii] = lThisTrackLengthReco_oTOF / Velocity(recoTrack.getP(), lMasses[ii]) - lMeasuredTime_oTOF;

        // Fixme: assumes dominant resolution effect is the TOF resolution
        // and not the tracking itself. It's *probably* a fair assumption
        // but it should be tested further!
        if (lThisTrackLength_iTOF > 0 && lThisTrackLengthReco_iTOF)
          lNSigma_iTOF[ii] = lDeltaTime_iTOF[ii] / iTOFTimeReso;
        if (lThisTrackLength_oTOF > 0 && lThisTrackLengthReco_oTOF)
          lNSigma_oTOF[ii] = lDeltaTime_oTOF[ii] / oTOFTimeReso;
      }

      // Sigmas have been fully calculated. Please populate the NSigma helper table (once per track)
      upgradetof(lNSigma_iTOF[0], lNSigma_iTOF[1], lNSigma_iTOF[2], lNSigma_iTOF[3], lNSigma_iTOF[4],
                 lNSigma_oTOF[0], lNSigma_oTOF[1], lNSigma_oTOF[2], lNSigma_oTOF[3], lNSigma_oTOF[4]);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<OnTheFlyTOFPID>(cfgc)};
}
