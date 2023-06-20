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
//     Lambdakzero TOF PID 
//  *+-+*+-+*+-+*+-+*+-+*+-+*
//
/// \author Nicolò Jacazio
/// \author David Dobrigkeit Chinellato
/// \since  11/05/2023
/// \brief  Table producer for V0 daughter TOF info
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

struct lambdakzerotofpid {
  Produces<aod::V0DeltaTimeTOF> v0DeltaTimeTOF;
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  // For manual sliceBy
  Preslice<aod::V0Datas> perCollision = o2::aod::v0data::collisionId;

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Operation and minimisation criteria
  Configurable<double> d_bz_input{"d_bz", -999, "bz field, -999 is automatic"};
  Configurable<int> nStepsLIntegrator{"nStepsLIntegrator", 200, "number of steps in length integrator"};

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

    const AxisSpec axisDeltaTime{(int)2000, -1000.0f, +1000.0f, "p_{T} (GeV/c)"};

    histos.add("hDeltaTimePositive", "hDeltaTimePositive", kTH1F, {axisDeltaTime});
    histos.add("hDeltaTimeNegative", "hDeltaTimeNegative", kTH1F, {axisDeltaTime});
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

  /// function to calculate track length
  /// \param track the input track
  /// \param x0 the initial position
  /// \param x1 the final position
  float trackLength(o2::track::TrackParCov track, float x0, float x1)
  {
    std::array<float, 3> pointN;
    std::array<float, 3> pointNplus;
    float length = 0.0;
    track.propagateTo(x0, d_bz);
    for (int iStep = 1; iStep < nStepsLIntegrator; iStep++) {
      track.getXYZGlo(pointN);
      float position = x0 + (x1 - x0) * (static_cast<float>(iStep)) / (static_cast<float>(nStepsLIntegrator - 1));
      track.propagateTo(position, d_bz);
      track.getXYZGlo(pointNplus);
      length += std::hypot(pointNplus[0] - pointN[0], pointNplus[1] - pointN[1], pointNplus[2] - pointN[2]);
    }
    return length;
  }

  float velocity(float lMomentum, float lMass){
    //Momentum p and mass m -> returns speed in centimeters per picosecond
    //Useful for TOF calculations
    float lA = (lMomentum / lMass)*(lMomentum / lMass);
    return 0.0299792458*TMath::Sqrt(lA/(1+lA));
  }

  void process(aod::Collisions const& collisions, aod::V0Datas const& V0s, FullTracksExtIU const&, aod::BCsWithTimestamps const&)
  {
    for (const auto& collision : collisions) {
      // Fire up CCDB
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      // Do analysis with collision-grouped V0s, retain full collision information
      const uint64_t collIdx = collision.globalIndex();
      auto V0Table_thisCollision = V0s.sliceBy(perCollision, collIdx);
      // V0 table sliced
      for(auto const& v0 : V0Table_thisCollision){
        // time of V0 segment
        float lengthV0 = std::hypot( v0.x() - collision.posX(), v0.y() - collision.posY(), v0.z() - collision.posZ() ); 
        float velocityV0 = velocity( v0.p(), o2::constants::physics::MassLambda ); 
        //float velocityV0K0short = velocity( v0.p(), o2::constants::physics::MassLambda ); 
        float timeV0 = lengthV0/velocityV0; //in picoseconds

        auto const& posTrackRow = v0.posTrack_as<FullTracksExtIU>();
        auto const& negTrackRow = v0.negTrack_as<FullTracksExtIU>();

        auto posTrack = getTrackParCov(posTrackRow);
        auto negTrack = getTrackParCov(negTrackRow);

        float posTofX = -1; 
        float negTofX = -1; 
        if (!posTrack.getXatLabR(370.0, posTofX, d_bz, o2::track::DirOutward)){
          posTofX = -1;
        }
        if (!negTrack.getXatLabR(370.0, negTofX, d_bz, o2::track::DirOutward)){
          negTofX = -1;
        }
        float deltaTimePositive = -1e+6;
        float deltaTimeNegative = -1e+6;
        if( posTofX > 10 ){
          float velocityPositive = velocity( posTrack.getP(), o2::constants::physics::MassProton );
          float lengthPositive = trackLength( posTrack, v0.posX(),  posTofX );
          float timePositive = lengthPositive/velocityPositive;
          deltaTimePositive = (posTrackRow.tofSignal() - posTrackRow.tofEvTime()) - (timeV0+timePositive) ; 
        }
        if( negTofX > 10 ){
          float velocityNegative = velocity( negTrack.getP(), o2::constants::physics::MassPionCharged );
          float lengthNegative = trackLength( negTrack, v0.negX(),  negTofX );
          float timeNegative = lengthNegative/velocityNegative;
          deltaTimeNegative = (negTrackRow.tofSignal() - negTrackRow.tofEvTime()) - (timeV0+timeNegative) ; 
        }
        v0DeltaTimeTOF(-1, deltaTimePositive, deltaTimeNegative, -1); 

        histos.fill(HIST("hDeltaTimePositive"), deltaTimePositive);
        histos.fill(HIST("hDeltaTimeNegative"), deltaTimeNegative);
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<lambdakzerotofpid>(cfgc)};
}
