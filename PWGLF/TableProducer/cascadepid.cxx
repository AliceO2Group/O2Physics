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
//     Cascade PID tables
//  *+-+*+-+*+-+*+-+*+-+*+-+*
//
/// \author Nicolò Jacazio
/// \author David Dobrigkeit Chinellato
/// \since  22/11/2023
/// \brief  Table producer for Casc daughter PID info
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
using TracksExtraWithPID = soa::Join<aod::TracksExtra, aod::pidTPCFullPi, aod::pidTPCFullPr, aod::pidTPCFullKa>;

// For MC and dE/dx association
using TracksExtraWithPIDandLabels = soa::Join<aod::TracksExtra, aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullPr, aod::pidTPCFullHe, aod::McTrackLabels>;

// Pre-selected V0s
using TaggedCascades = soa::Join<aod::Cascades, aod::CascTags>;

// For MC association in pre-selection
using LabeledTracksExtra = soa::Join<aod::TracksExtra, aod::McTrackLabels>;

struct cascadepid {
  // TOF pid for strangeness (recalculated with topology)
  Produces<aod::CascTOF> casctof;       // raw table for checks
  Produces<aod::CascTOFPID> casctofpid; // table with Nsigmas

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  // For manual sliceBy
  Preslice<aod::CascDatas> perCollision = o2::aod::cascdata::collisionId;

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Operation and minimisation criteria
  Configurable<double> d_bz_input{"d_bz", -999, "bz field, -999 is automatic"};
  Configurable<float> tofPosition{"tofPosition", 370, "TOF position for tests"};
  Configurable<bool> fillRawPID{"fillRawPID", true, "fill raw PID tables for debug/x-check"};

  // CCDB options
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};

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

  void process(aod::Collisions const& collisions, aod::CascDatas const& Cascades, FullTracksExtIU const&, aod::BCsWithTimestamps const&, aod::V0s const&)
  {
    for (const auto& collision : collisions) {
      // Fire up CCDB
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      // Do analysis with collision-grouped cascades, retain full collision information
      const uint64_t collIdx = collision.globalIndex();
      auto CascTable_thisCollision = Cascades.sliceBy(perCollision, collIdx);
      // cascade table sliced
      for (auto const& cascade : CascTable_thisCollision) {
        // Track casting
        auto bachTrack = cascade.bachelor_as<FullTracksExtIU>();
        auto v0 = cascade.v0();
        auto posTrack = v0.posTrack_as<FullTracksExtIU>();
        auto negTrack = v0.negTrack_as<FullTracksExtIU>();

        // FIXME: TOF calculation: under construction, to follow

        if (fillRawPID) {
          casctof(posTrack.length(), negTrack.length(), bachTrack.length(),
                  0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
                  0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<cascadepid>(cfgc)};
}
