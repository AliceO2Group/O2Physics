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

/// \file extended3bodyCreator.cxx

#include <cmath>
#include <array>
#include <cstdlib>
#include <string>
#include <algorithm>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "PWGLF/DataModel/pidTOFGeneric.h"
#include "PWGLF/DataModel/Reduced3BodyTables.h"
#include "PWGLF/DataModel/Reduced3BodyTablesLocal.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/Core/PID/PIDTOF.h"
#include "TableHelper.h"
#include "Tools/KFparticle/KFUtilities.h"

#include "EventFiltering/Zorro.h"
#include "EventFiltering/ZorroSummary.h"

#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"

#ifndef HomogeneousField
#define HomogeneousField
#endif

// includes KFParticle
#include "KFParticle.h"
#include "KFPTrack.h"
#include "KFPVertex.h"
#include "KFParticleBase.h"
#include "KFVertex.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using ReducedCollisionsMultsCents = soa::Join<aod::RedCollisions, aod::RedPVMults, aod::RedCentFT0Cs>;

struct extended3bodyCreator {

  Produces<aod::RedLocCollisions> reducedCollisions;
  Produces<aod::RedLocPVMults> reducedPVMults;
  Produces<aod::RedLocCentFT0Cs> reducedCentFT0Cs;
  Produces<aod::RedLocDecay3Bodys> reducedDecay3Bodys;
  Produces<aod::RedLoc3BodyInfo> reduced3BodyInfo;
  Produces<aod::StoredRedLocIUTracks> reducedFullTracksPIDIU;

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  // CCDB options
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};

  int mRunNumber;
  float d_bz;

  HistogramRegistry registry{"registry", {}};

  void init(InitContext&)
  {
    mRunNumber = 0;

    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);
  }

  void initCCDBfromRunNumber(int runNumber)
  {
    // Check if the CCDB data for this run is already cached
    if (ccdbCache.find(runNumber) != ccdbCache.end()) {
      LOG(debug) << "CCDB data already cached for run " << runNumber;

      // get magnetic field info from cache
      float d_bz = ccdbCache[runNumber];

      // Set magnetic field for KF vertexing
#ifdef HomogeneousField
      KFParticle::SetField(d_bz);
#endif
      // Set field for DCAfitter
      fitter3body.setBz(d_bz);

      if (useMatCorrType == 2) {
        // setMatLUT only after magfield has been initalized
        o2::base::Propagator::Instance()->setMatLUT(lut);
      }
    } else { // fetch data from CCDB and cache it
      o2::parameters::GRPMagField* grpmag = ccdb->getForRun<o2::parameters::GRPMagField>(grpmagPath, runNumber);
      if (!grpmag) {
        LOG(fatal) << "Got nullptr from CCDB for path " << grpmagPath << " of object GRPMagField and " << grpPath << " of object GRPObject for run number " << runNumber;
      }
      o2::base::Propagator::initFieldFromGRP(grpmag);
      // Fetch magnetic field from ccdb for current collision
      d_bz = o2::base::Propagator::Instance()->getNominalBz();
      LOG(info) << "Retrieved GRP for run number " << runNumber << " with magnetic field of " << d_bz << " kZG";

      // Set magnetic field for KF vertexing
#ifdef HomogeneousField
      KFParticle::SetField(d_bz);
#endif
      // Set field for DCAfitter
      fitter3body.setBz(d_bz);

      if (useMatCorrType == 2) {
        // setMatLUT only after magfield has been initalized
        o2::base::Propagator::Instance()->setMatLUT(lut);
      }

      // cache magnetic field info
      ccdbCache[runNumber] = d_bz;
    }
  }

  //------------------------------------------------------------------
  // function to fit KFParticle 3body vertex
  template <typename TKFParticle>
  void fit3bodyVertex(TKFParticle& kfpProton, TKFParticle& kfpPion, TKFParticle& kfpDeuteron, TKFParticle& KFHt)
  {
    // Construct 3body vertex
    int nDaughters3body = 3;
    const KFParticle* Daughters3body[3] = {&kfpProton, &kfpPion, &kfpDeuteron};
    KFHt.SetConstructMethod(2);
    try {
      KFHt.Construct(Daughters3body, nDaughters3body);
    } catch (std::runtime_error& e) {
      LOG(debug) << "Failed to create Hyper triton 3-body vertex." << e.what();
      return;
    }
    LOG(debug) << "Hypertriton vertex constructed.";
  }

  void processCollisions(ReducedCollisionsMultsCents const& collisions)
  {

    for (const auto& collision : collisions) {
      reducedCollisions(
        collision.posX(), collision.posY(), collision.posZ(),
        collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ(),
        collision.flags(), collision.chi2(), collision.numContrib(),
        collision.collisionTime(), collision.collisionTimeRes(),
        collision.runNumber());
      reducedPVMults(collision.multNTracksPV());
      reducedCentFT0Cs(collision.centFT0C());
    }
  }

  void processTracks(aod::RedIUTracks const& tracks)
  {
    for (const auto& track : tracks) {
      reducedFullTracksPIDIU(
        // TrackIU
        track.collisionId(),
        track.x(), track.alpha(),
        track.y(), track.z(), track.snp(), track.tgl(),
        track.signed1Pt(),
        // TracksCovIU
        track.sigmaY(), track.sigmaZ(), track.sigmaSnp(), track.sigmaTgl(), track.sigma1Pt(),
        track.rhoZY(), track.rhoSnpY(), track.rhoSnpZ(), track.rhoTglY(), track.rhoTglZ(),
        track.rhoTglSnp(), track.rho1PtY(), track.rho1PtZ(), track.rho1PtSnp(), track.rho1PtTgl(),
        // TracksExtra
        track.tpcInnerParam(), track.flags(), track.itsClusterSizes(),
        track.tpcNClsFindable(), track.tpcNClsFindableMinusFound(), track.tpcNClsFindableMinusCrossedRows(),
        track.trdPattern(), track.tpcChi2NCl(), track.tofChi2(),
        track.tpcSignal(), track.tofExpMom(),
        // PID
        track.tpcNSigmaPr(), track.tpcNSigmaPi(), track.tpcNSigmaDe(),
        track.tofNSigmaDe());
    }
  }

  void process3Bodys(ReducedCollisionsMultsCents const&, aod::RedIUTracks const&, aod::RedDecay3Bodys const& decay3bodys)
  {
    for (const auto& d3body : decay3bodys) {

      // save reduced decay3body table
      reducedDecay3Bodys(d3body.collisionId(), d3body.track0Id(), d3body.track1Id(), d3body.track2Id());

      // get collision
      auto collision = d3body.template collision_as<ReducedCollisionsMultsCents>();
      initCCDBfromRunNumber(collision.runNumber());

      // Get daughter tracks
      const auto daughter0 = d3body.template track0_as<RedIUTracks>();
      const auto daughter1 = d3body.template track1_as<RedIUTracks>();
      const auto daughter2 = d3body.template track2_as<RedIUTracks>();

      // get trackParCov daughters
      auto trackParCovPos = getTrackParCov(daughter0);
      auto trackParCovNeg = getTrackParCov(daughter1);
      auto trackParCovBach = getTrackParCov(daughter2);

      // create KFParticle daughters
      KFParticle kfpProton, kfpPion, kfpDeuteron;
      if (daughter2.sign() > 0) {
        kfpProton = createKFParticleFromTrackParCov(trackParCovPos, daughter0.sign(), constants::physics::MassProton);
        kfpPion = createKFParticleFromTrackParCov(trackParCovNeg, daughter1.sign(), constants::physics::MassPionCharged);
      } else if (!(daughter2.sign() > 0)) {
        kfpProton = createKFParticleFromTrackParCov(trackParCovNeg, daughter1.sign(), constants::physics::MassProton);
        kfpPion = createKFParticleFromTrackParCov(trackParCovPos, daughter0.sign(), constants::physics::MassPionCharged);
      }
      kfpDeuteron = createKFParticleFromTrackParCov(trackParCovBach, daughter2.sign(), constants::physics::MassDeuteron);

      // fit 3body vertex
      KFParticle KFHt;
      fit3bodyVertex(kfpProton, kfpPion, kfpDeuteron, KFHt);

      // calculate radius and phi
      auto radius = std::sqrt(KFHt.GetX() * KFHt.GetX() + KFHt.GetY() * KFHt.GetY());
      float phi, sigma;
      KFHt.GetPhi(phi, sigma);

      // fill 3body info table
      reduced3BodyInfo(radius, phi, KFHt.GetZ());
    } // end decay3body loop
  }
};

struct extended3bodyInitializer {
  Spawns<aod::RedLocIUTracks> reducedLocTracksIU;
  void init(InitContext const&) {}
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<extended3bodyInitializer>(cfgc),
    adaptAnalysisTask<extended3bodyCreator>(cfgc),
  };
}
