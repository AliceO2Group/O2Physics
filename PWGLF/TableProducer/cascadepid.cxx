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

// Cores with references and TOF pid
using dauTracks = soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>;
using CascFullCores = soa::Join<aod::CascCores, aod::CascTOFs, aod::CascExtras, aod::CascCollRefs>;

struct cascadepid {
  // TOF pid for strangeness (recalculated with topology)
  Produces<aod::CascTOFPIDs> casctofpids; // table with Nsigmas

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  // For manual sliceBy
  Preslice<CascFullCores> perCollision = o2::aod::v0data::straCollisionId;

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Operation and minimisation criteria
  Configurable<double> d_bz_input{"d_bz", -999, "bz field, -999 is automatic"};
  Configurable<float> tofPosition{"tofPosition", 377.934f, "TOF effective (inscribed) radius"};
  Configurable<bool> doQA{"doQA", true, "create QA histos"};
  Configurable<float> qaV0DCADau{"qaV0DCADau", 0.5, "DCA daughters (cm) for QA plots"};
  Configurable<float> qaCascDCADau{"qaCascDCADau", 0.5, "DCA daughters (cm) for QA plots"};
  Configurable<float> qaV0CosPA{"qaV0CosPA", 0.995, "CosPA for QA plots"};
  Configurable<float> qaCascCosPA{"qaCascCosPA", 0.995, "CosPA for QA plots"};
  Configurable<float> qaMassWindow{"qaMassWindow", 0.005, "Mass window around expected (in GeV/c2) for QA plots"};
  Configurable<float> qaTPCNSigma{"qaTPCNSigma", 5, "TPC N-sigma to apply for qa plots"};

  // CCDB options
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};

  ConfigurableAxis axisEta{"axisEta", {20, -1.0f, +1.0f}, "#eta"};
  ConfigurableAxis axisDeltaTime{"axisDeltaTime", {2000, -1000.0f, +1000.0f}, "delta-time (ps)"};
  ConfigurableAxis axisTime{"axisTime", {200, 0.0f, +20000.0f}, "T (ps)"};
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "p_{T} (GeV/c)"};

  int mRunNumber;
  float d_bz;
  float maxSnp;  // max sine phi for propagation
  float maxStep; // max step size (cm) for propagation

  /// function to calculate track length of this track up to a certain segment of a detector
  /// to be used internally in another funcrtion that calculates length until it finds the proper one
  /// warning: this could be optimised further for speed
  /// \param track the input track
  /// \param x1 x of the first point of the detector segment
  /// \param y1 y of the first point of the detector segment
  /// \param x2 x of the first point of the detector segment
  /// \param y2 y of the first point of the detector segment
  /// \param magneticField the magnetic field to use when propagating
  float trackLengthToSegment(o2::track::TrackPar track, float x1, float y1, float x2, float y2, float magneticField)
  {
    // don't make use of the track parametrization
    float length = -104;

    // causality protection
    std::array<float, 3> mom;
    track.getPxPyPzGlo(mom);
    // get start point
    std::array<float, 3> startPoint;
    track.getXYZGlo(startPoint);

    // better replaced with scalar momentum check later
    // if (((x1 + x2) * mom[0] + (y1 + y2) * mom[1]) < 0.0f)
    //   return -101;

    // get circle X, Y please
    o2::math_utils::CircleXYf_t trcCircle;
    float sna, csa;
    track.getCircleParams(magneticField, trcCircle, sna, csa);

    // Calculate necessary inner product
    float segmentModulus = std::hypot(x2 - x1, y2 - y1);
    float alongSegment = ((trcCircle.xC - x1) * (x2 - x1) + (trcCircle.yC - y1) * (y2 - y1)) / segmentModulus;

    // find point of closest approach between segment and circle center
    float pcaX = (x2 - x1) * alongSegment / segmentModulus + x1;
    float pcaY = (y2 - y1) * alongSegment / segmentModulus + y1;

    float centerDistToPC = std::hypot(pcaX - trcCircle.xC, pcaY - trcCircle.yC);

    // distance pca-to-intercept in multiples of segment modulus (for convenience)
    if (centerDistToPC > trcCircle.rC)
      return -103;

    float pcaToIntercept = TMath::Sqrt(TMath::Abs(trcCircle.rC * trcCircle.rC - centerDistToPC * centerDistToPC));

    float interceptX1 = pcaX + (x2 - x1) / segmentModulus * pcaToIntercept;
    float interceptY1 = pcaY + (y2 - y1) / segmentModulus * pcaToIntercept;
    float interceptX2 = pcaX - (x2 - x1) / segmentModulus * pcaToIntercept;
    float interceptY2 = pcaY - (y2 - y1) / segmentModulus * pcaToIntercept;

    float scalarCheck1 = ((x2 - x1) * (interceptX1 - x1) + (y2 - y1) * (interceptY1 - y1)) / segmentModulus;
    float scalarCheck2 = ((x2 - x1) * (interceptX2 - x1) + (y2 - y1) * (interceptY2 - y1)) / segmentModulus;

    float cosAngle1 = -1000, sinAngle1 = -1000, modulus1 = -1000;
    float cosAngle2 = -1000, sinAngle2 = -1000, modulus2 = -1000;
    float length1 = -1000, length2 = -1000;

    modulus1 = std::hypot(interceptX1 - trcCircle.xC, interceptY1 - trcCircle.yC) * std::hypot(startPoint[0] - trcCircle.xC, startPoint[1] - trcCircle.yC);
    cosAngle1 = (interceptX1 - trcCircle.xC) * (startPoint[0] - trcCircle.xC) + (interceptY1 - trcCircle.yC) * (startPoint[1] - trcCircle.yC);
    sinAngle1 = (interceptX1 - trcCircle.xC) * (startPoint[1] - trcCircle.yC) - (interceptY1 - trcCircle.yC) * (startPoint[0] - trcCircle.xC);
    cosAngle1 /= modulus1;
    sinAngle1 /= modulus1;
    length1 = trcCircle.rC * TMath::ACos(cosAngle1);
    length1 *= sqrt(1.0f + track.getTgl() * track.getTgl());

    modulus2 = std::hypot(interceptX2 - trcCircle.xC, interceptY2 - trcCircle.yC) * std::hypot(startPoint[0] - trcCircle.xC, startPoint[1] - trcCircle.yC);
    cosAngle2 = (interceptX2 - trcCircle.xC) * (startPoint[0] - trcCircle.xC) + (interceptY2 - trcCircle.yC) * (startPoint[1] - trcCircle.yC);
    sinAngle2 = (interceptX2 - trcCircle.xC) * (startPoint[1] - trcCircle.yC) - (interceptY2 - trcCircle.yC) * (startPoint[0] - trcCircle.xC);
    cosAngle2 /= modulus2;
    sinAngle2 /= modulus2;
    length2 = trcCircle.rC * TMath::ACos(cosAngle2);
    length2 *= sqrt(1.0f + track.getTgl() * track.getTgl());

    // rotate transverse momentum vector such that it is at intercepts
    float angle1 = TMath::ACos(cosAngle1);
    if (sinAngle1 < 0)
      angle1 *= -1.0f;
    float px1 = +TMath::Cos(angle1) * mom[0] + TMath::Sin(angle1) * mom[1];
    float py1 = -TMath::Sin(angle1) * mom[0] + TMath::Cos(angle1) * mom[1];

    float angle2 = TMath::ACos(cosAngle2);
    if (sinAngle2 < 0)
      angle2 *= -1.0f;
    float px2 = +TMath::Cos(angle2) * mom[0] + TMath::Sin(angle2) * mom[1];
    float py2 = -TMath::Sin(angle2) * mom[0] + TMath::Cos(angle2) * mom[1];

    float midSegX = 0.5f * (x2 + x1);
    float midSegY = 0.5f * (y2 + y1);

    float scalarMomentumCheck1 = px1 * midSegX + py1 * midSegY;
    float scalarMomentumCheck2 = px2 * midSegX + py2 * midSegY;

    if (scalarCheck1 > 0.0f && scalarCheck1 < segmentModulus && scalarMomentumCheck1 > 0.0f) {
      length = length1;
      // X = interceptX1; Y = interceptY1;
    }
    if (scalarCheck2 > 0.0f && scalarCheck2 < segmentModulus && scalarMomentumCheck2 > 0.0f) {
      length = length2;
      // X = interceptX2; Y = interceptY2;
    }
    return length;
  }

  /// function to calculate track length of this track up to a certain segmented detector
  /// \param track the input track
  /// \param magneticField the magnetic field to use when propagating
  float findInterceptLength(o2::track::TrackPar track, float magneticField)
  {
    float length = 1e+6;
    for (int iSeg = 0; iSeg < 18; iSeg++) {
      // Detector segmentation loop
      float segmentAngle = 20.0f / 180.0f * TMath::Pi();
      float theta = static_cast<float>(iSeg) * 20.0f / 180.0f * TMath::Pi();
      float halfWidth = tofPosition * TMath::Tan(0.5f * segmentAngle);
      float x1 = TMath::Cos(theta) * (-halfWidth) + TMath::Sin(theta) * tofPosition;
      float y1 = -TMath::Sin(theta) * (-halfWidth) + TMath::Cos(theta) * tofPosition;
      float x2 = TMath::Cos(theta) * (+halfWidth) + TMath::Sin(theta) * tofPosition;
      float y2 = -TMath::Sin(theta) * (+halfWidth) + TMath::Cos(theta) * tofPosition;
      float thisLength = trackLengthToSegment(track, x1, y1, x2, y2, magneticField);
      if (thisLength < length && thisLength > 0)
        length = thisLength;
    }
    if (length > 1e+5)
      length = -100; // force negative to avoid misunderstandings
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

    // measured vs expected total time QA
    if (doQA) {
      // standard deltaTime values
      histos.add("hArcDebug", "hArcDebug", kTH2F, {axisPt, {500, -5.0f, 10.0f}});
      histos.add("h2dposDeltaTimeAsXiPi", "h2dposDeltaTimeAsXiPi", {HistType::kTH3F, {axisPt, axisEta, axisDeltaTime}});
      histos.add("h2dposDeltaTimeAsXiPr", "h2dposDeltaTimeAsXiPr", {HistType::kTH3F, {axisPt, axisEta, axisDeltaTime}});
      histos.add("h2dnegDeltaTimeAsXiPi", "h2dnegDeltaTimeAsXiPi", {HistType::kTH3F, {axisPt, axisEta, axisDeltaTime}});
      histos.add("h2dnegDeltaTimeAsXiPr", "h2dnegDeltaTimeAsXiPr", {HistType::kTH3F, {axisPt, axisEta, axisDeltaTime}});
      histos.add("h2dbachDeltaTimeAsXiPi", "h2dbachDeltaTimeAsXiPi", {HistType::kTH3F, {axisPt, axisEta, axisDeltaTime}});

      histos.add("h2dposDeltaTimeAsOmPi", "h2dposDeltaTimeAsOmPi", {HistType::kTH3F, {axisPt, axisEta, axisDeltaTime}});
      histos.add("h2dposDeltaTimeAsOmPr", "h2dposDeltaTimeAsOmPr", {HistType::kTH3F, {axisPt, axisEta, axisDeltaTime}});
      histos.add("h2dnegDeltaTimeAsOmPi", "h2dnegDeltaTimeAsOmPi", {HistType::kTH3F, {axisPt, axisEta, axisDeltaTime}});
      histos.add("h2dnegDeltaTimeAsOmPr", "h2dnegDeltaTimeAsOmPr", {HistType::kTH3F, {axisPt, axisEta, axisDeltaTime}});
      histos.add("h2dbachDeltaTimeAsOmKa", "h2dbachDeltaTimeAsOmKa", {HistType::kTH3F, {axisPt, axisEta, axisDeltaTime}});
    }
  }

  void initCCDB(soa::Join<aod::StraCollisions, aod::StraStamps>::iterator const& collision)
  {
    if (mRunNumber == collision.runNumber()) {
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
      mRunNumber = collision.runNumber();
      return;
    }

    auto run3grp_timestamp = collision.timestamp();
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
    mRunNumber = collision.runNumber();
  }

  float velocity(float lMomentum, float lMass)
  {
    // Momentum p and mass m -> returns speed in centimeters per picosecond
    // Useful for TOF calculations
    float lA = (lMomentum / lMass) * (lMomentum / lMass);
    return 0.0299792458 * TMath::Sqrt(lA / (1 + lA));
  }

  void process(soa::Join<aod::StraCollisions, aod::StraStamps> const& collisions, CascFullCores const& Cascades, dauTracks const&)
  {
    for (const auto& collision : collisions) {
      // Fire up CCDB - based on StraCollisions for derived analysis
      initCCDB(collision);
      // Do analysis with collision-grouped V0s, retain full collision information
      const uint64_t collIdx = collision.globalIndex();
      auto CascTable_thisCollision = Cascades.sliceBy(perCollision, collIdx);
      // cascade table sliced
      for (auto const& cascade : CascTable_thisCollision) {
        // initialize from positions and momenta as needed
        o2::track::TrackPar posTrack = o2::track::TrackPar({cascade.xlambda(), cascade.ylambda(), cascade.zlambda()}, {cascade.pxpos(), cascade.pypos(), cascade.pzpos()}, +1);
        o2::track::TrackPar negTrack = o2::track::TrackPar({cascade.xlambda(), cascade.ylambda(), cascade.zlambda()}, {cascade.pxneg(), cascade.pyneg(), cascade.pzneg()}, -1);
        o2::track::TrackPar bachTrack = o2::track::TrackPar({cascade.x(), cascade.y(), cascade.z()}, {cascade.pxbach(), cascade.pybach(), cascade.pzbach()}, cascade.sign());
        o2::track::TrackPar cascTrack = o2::track::TrackPar({cascade.x(), cascade.y(), cascade.z()}, {cascade.px(), cascade.py(), cascade.pz()}, cascade.sign());

        // start calculation: calculate velocities
        float velocityPositivePr = velocity(posTrack.getP(), o2::constants::physics::MassProton);
        float velocityPositivePi = velocity(posTrack.getP(), o2::constants::physics::MassPionCharged);
        float velocityNegativePr = velocity(negTrack.getP(), o2::constants::physics::MassProton);
        float velocityNegativePi = velocity(negTrack.getP(), o2::constants::physics::MassPionCharged);
        float velocityBachelorPi = velocity(bachTrack.getP(), o2::constants::physics::MassPionCharged);
        float velocityBachelorKa = velocity(bachTrack.getP(), o2::constants::physics::MassKaonCharged);
        float velocityXi = velocity(cascTrack.getP(), o2::constants::physics::MassXiMinus);
        float velocityOm = velocity(cascTrack.getP(), o2::constants::physics::MassOmegaMinus);
        float velocityLa = velocity(std::hypot(cascade.pxlambda(), cascade.pylambda(), cascade.pzlambda()), o2::constants::physics::MassLambda);

        // calculate daughter length to TOF intercept
        float lengthPositive = findInterceptLength(posTrack, d_bz);  // FIXME: tofPosition ok? adjust?
        float lengthNegative = findInterceptLength(negTrack, d_bz);  // FIXME: tofPosition ok? adjust?
        float lengthBachelor = findInterceptLength(bachTrack, d_bz); // FIXME: tofPosition ok? adjust?

        // calculate mother lengths
        float lengthV0 = std::hypot(cascade.xlambda() - cascade.x(), cascade.ylambda() - cascade.y(), cascade.zlambda() - cascade.z());
        float lengthCascade = -1e+6;
        const o2::math_utils::Point3D<float> collVtx{collision.posX(), collision.posY(), collision.posZ()};
        bool successPropag = o2::base::Propagator::Instance()->propagateToDCA(collVtx, cascTrack, d_bz, 2.f, o2::base::Propagator::MatCorrType::USEMatCorrNONE);
        float d = -1.0f, d3d = 0.0f;
        float linearToPV = std::hypot(cascade.x() - collision.posX(), cascade.y() - collision.posY(), cascade.z() - collision.posZ());
        if (successPropag) {
          std::array<float, 3> cascCloseToPVPosition;
          cascTrack.getXYZGlo(cascCloseToPVPosition);
          o2::math_utils::CircleXYf_t trcCircleCascade;
          float sna, csa;
          cascTrack.getCircleParams(d_bz, trcCircleCascade, sna, csa);

          // calculate 2D distance between two points
          d = std::hypot(cascade.x() - cascCloseToPVPosition[0], cascade.y() - cascCloseToPVPosition[1]);
          d3d = std::hypot(cascade.x() - cascCloseToPVPosition[0], cascade.y() - cascCloseToPVPosition[1], cascade.z() - cascCloseToPVPosition[2]); // cross-check variable
          float sinThetaOverTwo = d / (2.0f * trcCircleCascade.rC);
          lengthCascade = 2.0f * trcCircleCascade.rC * TMath::ASin(sinThetaOverTwo);
          lengthCascade *= sqrt(1.0f + cascTrack.getTgl() * cascTrack.getTgl());
        }

        if (!successPropag) {
          lengthCascade = linearToPV; // if propagation failed, use linear estimate (optional: actually do not define?)
        }

        // lambda, xi and omega flight time is always defined
        float lambdaFlight = lengthV0 / velocityLa;
        float xiFlight = lengthCascade / velocityXi;
        float omFlight = lengthCascade / velocityOm;
        float posFlightPi = lengthPositive / velocityPositivePi;
        float posFlightPr = lengthPositive / velocityPositivePr;
        float negFlightPi = lengthNegative / velocityNegativePi;
        float negFlightPr = lengthNegative / velocityNegativePr;
        float bachFlightPi = lengthBachelor / velocityBachelorPi;
        float bachFlightKa = lengthBachelor / velocityBachelorKa;

        // initialize delta-times (actual PID variables)
        float posDeltaTimeAsXiPi = -1e+6, posDeltaTimeAsXiPr = -1e+6;
        float negDeltaTimeAsXiPi = -1e+6, negDeltaTimeAsXiPr = -1e+6;
        float bachDeltaTimeAsXiPi = -1e+6;
        float posDeltaTimeAsOmPi = -1e+6, posDeltaTimeAsOmPr = -1e+6;
        float negDeltaTimeAsOmPi = -1e+6, negDeltaTimeAsOmPr = -1e+6;
        float bachDeltaTimeAsOmKa = -1e+6;

        if (cascade.posTOFSignal() > 0 && cascade.posTOFEventTime() > 0) {
          posDeltaTimeAsXiPi = (cascade.posTOFSignal() - cascade.posTOFEventTime()) - (xiFlight + lambdaFlight + posFlightPi);
          posDeltaTimeAsXiPr = (cascade.posTOFSignal() - cascade.posTOFEventTime()) - (xiFlight + lambdaFlight + posFlightPr);
          posDeltaTimeAsOmPi = (cascade.posTOFSignal() - cascade.posTOFEventTime()) - (omFlight + lambdaFlight + posFlightPi);
          posDeltaTimeAsOmPr = (cascade.posTOFSignal() - cascade.posTOFEventTime()) - (omFlight + lambdaFlight + posFlightPr);
        }
        if (cascade.negTOFSignal() > 0 && cascade.negTOFEventTime() > 0) {
          negDeltaTimeAsXiPi = (cascade.negTOFSignal() - cascade.negTOFEventTime()) - (xiFlight + lambdaFlight + negFlightPi);
          negDeltaTimeAsXiPr = (cascade.negTOFSignal() - cascade.negTOFEventTime()) - (xiFlight + lambdaFlight + negFlightPr);
          negDeltaTimeAsOmPi = (cascade.negTOFSignal() - cascade.negTOFEventTime()) - (omFlight + lambdaFlight + negFlightPi);
          negDeltaTimeAsOmPr = (cascade.negTOFSignal() - cascade.negTOFEventTime()) - (omFlight + lambdaFlight + negFlightPr);
        }
        if (cascade.bachTOFSignal() > 0 && cascade.bachTOFEventTime() > 0) {
          bachDeltaTimeAsXiPi = (cascade.bachTOFSignal() - cascade.bachTOFEventTime()) - (xiFlight + bachFlightPi);
          bachDeltaTimeAsOmKa = (cascade.bachTOFSignal() - cascade.bachTOFEventTime()) - (omFlight + bachFlightKa);
        }

        casctofpids(
          posDeltaTimeAsXiPi, posDeltaTimeAsXiPr, negDeltaTimeAsXiPi, negDeltaTimeAsXiPr, bachDeltaTimeAsXiPi,
          posDeltaTimeAsOmPi, posDeltaTimeAsOmPr, negDeltaTimeAsOmPi, negDeltaTimeAsOmPr, bachDeltaTimeAsOmKa,
          0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f // no Nsigmas yet, note: could be fewer
        );

        if (doQA) {
          auto pTra = cascade.posTrackExtra_as<dauTracks>();
          auto nTra = cascade.negTrackExtra_as<dauTracks>();
          auto bTra = cascade.negTrackExtra_as<dauTracks>();

          // fill QA histograms for cross-checking
          histos.fill(HIST("hArcDebug"), cascade.pt(), lengthCascade - d3d); // for debugging purposes

          if (cascade.dcaV0daughters() < qaV0DCADau && cascade.dcacascdaughters() < qaCascDCADau && cascade.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) > qaV0CosPA && cascade.casccosPA(collision.posX(), collision.posY(), collision.posZ()) > qaCascCosPA) {
            if (cascade.sign() < 0) {
              if (std::abs(cascade.mXi() - 1.32171) < qaMassWindow && fabs(pTra.tpcNSigmaPr()) < qaTPCNSigma && fabs(nTra.tpcNSigmaPi()) < qaTPCNSigma && fabs(bTra.tpcNSigmaPi()) < qaTPCNSigma) {
                histos.fill(HIST("h2dposDeltaTimeAsXiPr"), cascade.pt(), cascade.eta(), posDeltaTimeAsXiPr);
                histos.fill(HIST("h2dnegDeltaTimeAsXiPi"), cascade.pt(), cascade.eta(), negDeltaTimeAsXiPi);
                histos.fill(HIST("h2dbachDeltaTimeAsXiPi"), cascade.pt(), cascade.eta(), bachDeltaTimeAsXiPi);
              }
              if (std::abs(cascade.mOmega() - 1.67245) < qaMassWindow && fabs(pTra.tpcNSigmaPr()) < qaTPCNSigma && fabs(nTra.tpcNSigmaPi()) < qaTPCNSigma && fabs(bTra.tpcNSigmaKa()) < qaTPCNSigma) {
                histos.fill(HIST("h2dposDeltaTimeAsOmPr"), cascade.pt(), cascade.eta(), posDeltaTimeAsOmPr);
                histos.fill(HIST("h2dnegDeltaTimeAsOmPi"), cascade.pt(), cascade.eta(), negDeltaTimeAsOmPi);
                histos.fill(HIST("h2dbachDeltaTimeAsOmKa"), cascade.pt(), cascade.eta(), bachDeltaTimeAsOmKa);
              }
            } else {
              if (std::abs(cascade.mXi() - 1.32171) < qaMassWindow && fabs(pTra.tpcNSigmaPi()) < qaTPCNSigma && fabs(nTra.tpcNSigmaPr()) < qaTPCNSigma && fabs(bTra.tpcNSigmaPi()) < qaTPCNSigma) {
                histos.fill(HIST("h2dposDeltaTimeAsXiPi"), cascade.pt(), cascade.eta(), posDeltaTimeAsXiPi);
                histos.fill(HIST("h2dnegDeltaTimeAsXiPr"), cascade.pt(), cascade.eta(), negDeltaTimeAsXiPr);
                histos.fill(HIST("h2dbachDeltaTimeAsXiPi"), cascade.pt(), cascade.eta(), bachDeltaTimeAsXiPi);
              }
              if (std::abs(cascade.mOmega() - 1.67245) < qaMassWindow && fabs(pTra.tpcNSigmaPi()) < qaTPCNSigma && fabs(nTra.tpcNSigmaPr()) < qaTPCNSigma && fabs(bTra.tpcNSigmaKa()) < qaTPCNSigma) {
                histos.fill(HIST("h2dposDeltaTimeAsOmPi"), cascade.pt(), cascade.eta(), posDeltaTimeAsOmPi);
                histos.fill(HIST("h2dnegDeltaTimeAsOmPr"), cascade.pt(), cascade.eta(), negDeltaTimeAsOmPr);
                histos.fill(HIST("h2dbachDeltaTimeAsOmKa"), cascade.pt(), cascade.eta(), bachDeltaTimeAsOmKa);
              }
            }
          }
        }
      }
    }
  }
};

Configurable<float> qaV0DCADau{"qaV0DCADau", 0.5, "DCA daughters (cm) for QA plots"};
Configurable<float> qaCascDCADau{"qaCascDCADau", 0.5, "DCA daughters (cm) for QA plots"};
Configurable<float> qaV0CosPA{"qaV0CosPA", 0.995, "CosPA for QA plots"};
Configurable<float> qaCascCosPA{"qaCascCosPA", 0.995, "CosPA for QA plots"};
Configurable<float> qaMassWindow{"qaMassWindow", 0.005, "Mass window around expected (in GeV/c2) for QA plots"};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<cascadepid>(cfgc)};
}
