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
//  *+-+*+-+*+-+*+-+*+-+*+-+*
//     Lambdakzero PID
//  *+-+*+-+*+-+*+-+*+-+*+-+*
//
/// \author Nicolò Jacazio
/// \author David Dobrigkeit Chinellato
/// \since  11/05/2023
/// \brief  Table producer for V0 daughter PID info
//
// This task produces daughter PID information for strange daughters
// taking into account the (candidate-by-candidate) time spent as a heavier
// (strange, weakly-decaying) particle. This task is meant to be a test, as
// it hasn't been fully tested yet! Use at your own peril for now :-)

#include <cmath>
#include <array>
#include <cstdlib>
#include <map>
#include <iterator>
#include <utility>

#include "Framework/runDataProcessing.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "DCAFitter/DCAFitterN.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "PWGLF/DataModel/LFParticleIdentification.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Common/TableProducer/PID/pidTOFBase.h"
#include "Common/DataModel/PIDResponse.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

// use parameters + cov mat non-propagated, aux info + (extension propagated)
using FullTracksExt = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov>;
using FullTracksExtIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TOFEvTime, aod::TOFSignal>;
using TracksWithExtra = soa::Join<aod::Tracks, aod::TracksExtra>;

// For dE/dx association in pre-selection
using TracksExtraWithPID = soa::Join<aod::TracksExtra, aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullPr, aod::pidTPCFullHe>;

// For MC and dE/dx association
using TracksExtraWithPIDandLabels = soa::Join<aod::TracksExtra, aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullPr, aod::pidTPCFullHe, aod::McTrackLabels>;

// Pre-selected V0s
using TaggedV0s = soa::Join<aod::V0s, aod::V0Tags>;

// For MC association in pre-selection
using LabeledTracksExtra = soa::Join<aod::TracksExtra, aod::McTrackLabels>;

struct lambdakzeropid {
  // TOF pid for strangeness (recalculated with topology)
  Produces<aod::V0TOF> v0tof;       // raw table for checks
  Produces<aod::V0TOFPID> v0tofpid; // table with Nsigmas

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  // For manual sliceBy
  Preslice<aod::V0Datas> perCollision = o2::aod::v0data::collisionId;

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Operation and minimisation criteria
  Configurable<double> d_bz_input{"d_bz", -999, "bz field, -999 is automatic"};
  Configurable<float> tofPosition{"tofPosition", 370, "TOF position for tests"};
  Configurable<bool> checkTPCCompatibility{"checkTPCCompatibility", true, "check compatibility with dE/dx in QA plots"};
  Configurable<bool> fillRawPID{"fillRawPID", true, "fill raw PID tables for debug/x-check"};

  // CCDB options
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};

  ConfigurableAxis axisPtQA{"axisPtQA", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "pt axis for QA histograms"};
  ConfigurableAxis axisDeltaTime{"axisDeltaTime", {2000, -1000.0f, +1000.0f}, "delta-time (ps)"};
  ConfigurableAxis axisK0ShortMass{"axisK0ShortMass", {200, 0.400f, 0.600f}, "Inv. Mass (GeV/c^{2})"};
  ConfigurableAxis axisLambdaMass{"axisLambdaMass", {200, 1.01f, 1.21f}, "Inv. Mass (GeV/c^{2})"};

  int mRunNumber;
  float d_bz;
  float maxSnp;  // max sine phi for propagation
  float maxStep; // max step size (cm) for propagation

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

  void init(InitContext& context)
  {
    mRunNumber = 0;
    d_bz = 0;
    maxSnp = 0.85f;  // could be changed later
    maxStep = 2.00f; // could be changed later

    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    histos.add("h3dMassK0ShortPositive", "h3dMassK0ShortPositive", kTH3F, {axisPtQA, axisDeltaTime, axisK0ShortMass});
    histos.add("h3dMassLambdaPositive", "h3dMassLambdaPositive", kTH3F, {axisPtQA, axisDeltaTime, axisLambdaMass});
    histos.add("h3dMassAntiLambdaPositive", "h3dMassAntiLambdaPositive", kTH3F, {axisPtQA, axisDeltaTime, axisLambdaMass});
    histos.add("h3dMassK0ShortNegative", "h3dMassK0ShortNegative", kTH3F, {axisPtQA, axisDeltaTime, axisK0ShortMass});
    histos.add("h3dMassLambdaNegative", "h3dMassLambdaNegative", kTH3F, {axisPtQA, axisDeltaTime, axisLambdaMass});
    histos.add("h3dMassAntiLambdaNegative", "h3dMassAntiLambdaNegative", kTH3F, {axisPtQA, axisDeltaTime, axisLambdaMass});
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }

    // In case override, don't proceed, please - no CCDB access required
    if (d_bz_input > -990) {
      d_bz = d_bz_input;
      o2::parameters::GRPMagField grpmag;
      if (fabs(d_bz) > 1e-5) {
        grpmag.setL3Current(30000.f / (d_bz / 5.0f));
      }
      o2::base::Propagator::initFieldFromGRP(&grpmag);
      mRunNumber = bc.runNumber();
      return;
    }

    auto run3grp_timestamp = bc.timestamp();
    o2::parameters::GRPObject* grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(grpPath, run3grp_timestamp);
    o2::parameters::GRPMagField* grpmag = 0x0;
    if (grpo) {
      o2::base::Propagator::initFieldFromGRP(grpo);
      // Fetch magnetic field from ccdb for current collision
      d_bz = grpo->getNominalL3Field();
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
    } else {
      grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, run3grp_timestamp);
      if (!grpmag) {
        LOG(fatal) << "Got nullptr from CCDB for path " << grpmagPath << " of object GRPMagField and " << grpPath << " of object GRPObject for timestamp " << run3grp_timestamp;
      }
      o2::base::Propagator::initFieldFromGRP(grpmag);
      // Fetch magnetic field from ccdb for current collision
      d_bz = std::lround(5.f * grpmag->getL3Current() / 30000.f);
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
    }
    mRunNumber = bc.runNumber();
  }

  float velocity(float lMomentum, float lMass)
  {
    // Momentum p and mass m -> returns speed in centimeters per picosecond
    // Useful for TOF calculations
    float lA = (lMomentum / lMass) * (lMomentum / lMass);
    return 0.0299792458 * TMath::Sqrt(lA / (1 + lA));
  }

  void process(aod::Collisions const& collisions, aod::V0Datas const& V0s, FullTracksExtIU const&, aod::BCsWithTimestamps const&, TaggedV0s const& allV0s)
  {
    for (const auto& collision : collisions) {
      // Fire up CCDB
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      // Do analysis with collision-grouped V0s, retain full collision information
      const uint64_t collIdx = collision.globalIndex();
      auto V0Table_thisCollision = V0s.sliceBy(perCollision, collIdx);
      // V0 table sliced
      for (auto const& v0 : V0Table_thisCollision) {
        // time of V0 segment
        float lengthV0 = std::hypot(v0.x() - collision.posX(), v0.y() - collision.posY(), v0.z() - collision.posZ());
        float velocityK0Short = velocity(v0.p(), o2::constants::physics::MassKaonNeutral);
        float velocityLambda = velocity(v0.p(), o2::constants::physics::MassLambda);
        float timeK0Short = lengthV0 / velocityK0Short; // in picoseconds
        float timeLambda = lengthV0 / velocityLambda;   // in picoseconds

        auto const& posTrackRow = v0.posTrack_as<FullTracksExtIU>();
        auto const& negTrackRow = v0.negTrack_as<FullTracksExtIU>();

        auto posTrack = getTrackParCov(posTrackRow);
        auto negTrack = getTrackParCov(negTrackRow);

        float deltaTimePositiveLambdaPi = -1e+6;
        float deltaTimeNegativeLambdaPi = -1e+6;
        float deltaTimePositiveLambdaPr = -1e+6;
        float deltaTimeNegativeLambdaPr = -1e+6;
        float deltaTimePositiveK0ShortPi = -1e+6;
        float deltaTimeNegativeK0ShortPi = -1e+6;

        float velocityPositivePr = velocity(posTrack.getP(), o2::constants::physics::MassProton);
        float velocityPositivePi = velocity(posTrack.getP(), o2::constants::physics::MassPionCharged);
        float velocityNegativePr = velocity(negTrack.getP(), o2::constants::physics::MassProton);
        float velocityNegativePi = velocity(negTrack.getP(), o2::constants::physics::MassPionCharged);

        // propagate to V0 decay vertex
        posTrack.propagateTo(v0.posX(), d_bz);
        negTrack.propagateTo(v0.negX(), d_bz);

        float lengthPositive = trackLength(posTrack, tofPosition, d_bz); // FIXME: tofPosition ok? adjust?
        float lengthNegative = trackLength(negTrack, tofPosition, d_bz); // FIXME: tofPosition ok? adjust?
        float timePositivePr = lengthPositive / velocityPositivePr;
        float timePositivePi = lengthPositive / velocityPositivePi;
        float timeNegativePr = lengthNegative / velocityNegativePr;
        float timeNegativePi = lengthNegative / velocityNegativePi;

        deltaTimePositiveLambdaPr = (posTrackRow.tofSignal() - posTrackRow.tofEvTime()) - (timeLambda + timePositivePr);
        deltaTimePositiveLambdaPi = (posTrackRow.tofSignal() - posTrackRow.tofEvTime()) - (timeLambda + timePositivePi);
        deltaTimeNegativeLambdaPr = (negTrackRow.tofSignal() - negTrackRow.tofEvTime()) - (timeLambda + timeNegativePr);
        deltaTimeNegativeLambdaPi = (negTrackRow.tofSignal() - negTrackRow.tofEvTime()) - (timeLambda + timeNegativePi);
        deltaTimePositiveK0ShortPi = (posTrackRow.tofSignal() - posTrackRow.tofEvTime()) - (timeK0Short + timeNegativePi);
        deltaTimeNegativeK0ShortPi = (negTrackRow.tofSignal() - negTrackRow.tofEvTime()) - (timeK0Short + timeNegativePi);

        if (fillRawPID) {
          v0tof(posTrackRow.length(), negTrackRow.length(),
                deltaTimePositiveLambdaPi, deltaTimePositiveLambdaPr,
                deltaTimeNegativeLambdaPi, deltaTimeNegativeLambdaPr,
                deltaTimePositiveK0ShortPi, deltaTimeNegativeK0ShortPi);
        }

        auto originalV0 = v0.v0_as<TaggedV0s>(); // this could look confusing, so:
        // the first v0 is the v0data row; the getter de-references the v0 (stored indices) row
        // the v0 (stored indices) contain the tags of the lambdakzero preselector

        if (originalV0.isdEdxK0Short() || !checkTPCCompatibility) {
          histos.fill(HIST("h3dMassK0ShortPositive"), v0.pt(), deltaTimePositiveK0ShortPi, v0.mK0Short());
          histos.fill(HIST("h3dMassK0ShortNegative"), v0.pt(), deltaTimePositiveK0ShortPi, v0.mK0Short());
        }
        if (originalV0.isdEdxLambda() || !checkTPCCompatibility) {
          histos.fill(HIST("h3dMassLambdaPositive"), v0.pt(), deltaTimePositiveLambdaPr, v0.mLambda());
          histos.fill(HIST("h3dMassLambdaNegative"), v0.pt(), deltaTimeNegativeLambdaPi, v0.mLambda());
        }
        if (originalV0.isdEdxAntiLambda() || !checkTPCCompatibility) {
          histos.fill(HIST("h3dMassAntiLambdaPositive"), v0.pt(), deltaTimePositiveK0ShortPi, v0.mAntiLambda());
          histos.fill(HIST("h3dMassAntiLambdaNegative"), v0.pt(), deltaTimeNegativeK0ShortPi, v0.mAntiLambda());
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<lambdakzeropid>(cfgc)};
}
