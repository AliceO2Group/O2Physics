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

/// \file reduced3bodyCreator.cxx
/// \brief Task to produce reduced AO2Ds for use in the hypertriton 3body reconstruction with the decay3bodybuilder.cxx
/// \author Yuanzhe Wang <yuanzhe.wang@cern.ch>
/// \author Carolina Reetz <c.reetz@cern.ch>

#include "TableHelper.h"

#include "PWGLF/DataModel/LFPIDTOFGenericTables.h"
#include "PWGLF/DataModel/Reduced3BodyTables.h"
#include "PWGLF/Utils/pidTOFGeneric.h"

#include "Common/Core/PID/PIDTOF.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/Zorro.h"
#include "Common/Core/ZorroSummary.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Tools/KFparticle/KFUtilities.h"

#include "CCDB/BasicCCDBManager.h"
#include "DCAFitter/DCAFitterN.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>

#ifndef HomogeneousField
#define HomogeneousField
#endif

// includes KFParticle
#include "KFPTrack.h"
#include "KFPVertex.h"
#include "KFParticle.h"
#include "KFParticleBase.h"
#include "KFVertex.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

o2::common::core::MetadataHelper metadataInfo;

using FullTracksExtIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU>;
using FullTracksExtPIDIU = soa::Join<FullTracksExtIU, aod::pidTPCFullPr, aod::pidTPCFullPi, aod::pidTPCFullDe>;

using ColwithEvTimes = o2::soa::Join<aod::Collisions, aod::EvSels, aod::EvTimeTOFFT0>;
using ColwithEvTimesMultsCents = o2::soa::Join<ColwithEvTimes, aod::PVMults, aod::CentFT0Cs>;
using TrackExtIUwithEvTimes = soa::Join<FullTracksExtIU, aod::EvTimeTOFFT0ForTrack>;
using TrackExtPIDIUwithEvTimes = soa::Join<FullTracksExtPIDIU, aod::EvTimeTOFFT0ForTrack>;

struct reduced3bodyCreator {

  Produces<aod::RedCollisions> reducedCollisions;
  Produces<aod::RedPVMults> reducedPVMults;
  Produces<aod::RedCentFT0Cs> reducedCentFT0Cs;
  Produces<aod::RedDecay3Bodys> reducedDecay3Bodys;
  Produces<aod::Red3BodyInfo> reduced3BodyInfo;
  Produces<aod::StoredRedIUTracks> reducedFullTracksPIDIU;

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Zorro zorro;
  OutputObj<ZorroSummary> zorroSummary{"zorroSummary"};

  o2::vertexing::DCAFitterN<3> fitter3body;
  o2::aod::pidtofgeneric::TofPidNewCollision<TrackExtPIDIUwithEvTimes::iterator> bachelorTOFPID;

  Configurable<bool> disableITSROFCut{"disableITSROFCut", false, "Disable ITS ROF border cut"};
  // CCDB options
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  // CCDB TOF PID paras
  Configurable<int64_t> timestamp{"ccdb-timestamp", -1, "timestamp of the object"};
  Configurable<std::string> paramFileName{"paramFileName", "", "Path to the parametrization object. If empty the parametrization is not taken from file"};
  Configurable<std::string> parametrizationPath{"parametrizationPath", "TOF/Calib/Params", "Path of the TOF parametrization on the CCDB or in the file, if the paramFileName is not empty"};
  Configurable<std::string> passName{"passName", "", "Name of the pass inside of the CCDB parameter collection. If empty, the automatically deceted from metadata (to be implemented!!!)"};
  Configurable<std::string> timeShiftCCDBPath{"timeShiftCCDBPath", "", "Path of the TOF time shift vs eta. If empty none is taken"};
  Configurable<bool> loadResponseFromCCDB{"loadResponseFromCCDB", false, "Flag to load the response from the CCDB"};
  Configurable<bool> fatalOnPassNotAvailable{"fatalOnPassNotAvailable", true, "Flag to throw a fatal if the pass is not available in the retrieved CCDB object"};
  // Zorro counting
  Configurable<bool> cfgSkimmedProcessing{"cfgSkimmedProcessing", false, "Skimmed dataset processing"};
  // Flag for trigger
  Configurable<bool> cfgOnlyKeepInterestedTrigger{"cfgOnlyKeepInterestedTrigger", false, "Flag to keep only interested trigger"};
  Configurable<std::string> triggerList{"triggerList", "fH3L3Body", "List of triggers used to select events"};

  Configurable<int> cfgMaterialCorrection{"cfgMaterialCorrection", static_cast<int>(o2::base::Propagator::MatCorrType::USEMatCorrNONE), "Type of material correction for DCAFitter"};

  int mRunNumber;
  float mBz;
  // TOF response and input parameters
  o2::pid::tof::TOFResoParamsV3 mRespParamsV3;
  o2::aod::pidtofgeneric::TOFCalibConfig mTOFCalibConfig; // TOF Calib configuration

  // tracked cluster size
  std::vector<int> fTrackedClSizeVector;

  HistogramRegistry registry{"registry", {}};

  void init(InitContext& initContext)
  {
    mRunNumber = 0;
    zorroSummary.setObject(zorro.getZorroSummary());
    bachelorTOFPID.SetPidType(o2::track::PID::Deuteron);

    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    // TOF PID parameters initialization
    mTOFCalibConfig.metadataInfo = metadataInfo;
    mTOFCalibConfig.inheritFromBaseTask(initContext);
    mTOFCalibConfig.initSetup(mRespParamsV3, ccdb);

    fitter3body.setPropagateToPCA(true);
    fitter3body.setMaxR(200.); //->maxRIni3body
    fitter3body.setMinParamChange(1e-3);
    fitter3body.setMinRelChi2Change(0.9);
    fitter3body.setMaxDZIni(1e9);
    fitter3body.setMaxChi2(1e9);
    fitter3body.setUseAbsDCA(true);
    int mat{static_cast<int>(cfgMaterialCorrection)};
    fitter3body.setMatCorrType(static_cast<o2::base::Propagator::MatCorrType>(mat));

    registry.add("hAllSelEventsVtxZ", "hAllSelEventsVtxZ", HistType::kTH1F, {{500, -15.0f, 15.0f, "PV Z (cm)"}});

    auto hEventCounter = registry.add<TH1>("hEventCounter", "hEventCounter", HistType::kTH1D, {{3, 0.0f, 3.0f}});
    hEventCounter->GetXaxis()->SetBinLabel(1, "all");
    hEventCounter->GetXaxis()->SetBinLabel(2, "selected");
    hEventCounter->GetXaxis()->SetBinLabel(3, "reduced");
    hEventCounter->LabelsOption("v");

    auto hEventCounterZorro = registry.add<TH1>("hEventCounterZorro", "hEventCounterZorro", HistType::kTH1D, {{2, 0, 2}});
    hEventCounterZorro->GetXaxis()->SetBinLabel(1, "Zorro before evsel");
    hEventCounterZorro->GetXaxis()->SetBinLabel(2, "Zorro after evsel");
  }

  void initZorroBC(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (cfgSkimmedProcessing) {
      zorro.initCCDB(ccdb.service, bc.runNumber(), bc.timestamp(), triggerList);
      zorro.populateHistRegistry(registry, bc.runNumber());
    }
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    // In case override, don't proceed, please - no CCDB access required
    if (mRunNumber == bc.runNumber()) {
      return;
    }
    mRunNumber = bc.runNumber();

    // In case override, don't proceed, please - no CCDB access required
    auto run3grp_timestamp = bc.timestamp();
    ccdb->clearCache(grpmagPath);
    o2::parameters::GRPMagField* grpmag = ccdb->getSpecific<o2::parameters::GRPMagField>(grpmagPath, run3grp_timestamp);
    if (!grpmag) {
      LOG(fatal) << "Got nullptr from CCDB for path " << grpmagPath << " of object GRPMagField and " << grpPath << " of object GRPObject for timestamp " << run3grp_timestamp;
    }
    o2::base::Propagator::initFieldFromGRP(grpmag);
    // Fetch magnetic field from ccdb for current collision
    // mBz = std::lround(5.f * grpmag->getL3Current() / 30000.f);
    mBz = o2::base::Propagator::Instance()->getNominalBz();
    LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << mBz << " kZG";
// Set magnetic field for KF vertexing
#ifdef HomogeneousField
    KFParticle::SetField(mBz);
#endif

    fitter3body.setBz(mBz);
    mTOFCalibConfig.processSetup(mRespParamsV3, ccdb, bc);
  }

  //------------------------------------------------------------------
  // function to fill reduced track table
  template <typename TTrack>
  void fillTrackTable(TTrack const& daughter, double tofNSigmaTrack, auto collisionIndex)
  {
    reducedFullTracksPIDIU(
      // TrackIU
      collisionIndex,
      daughter.x(), daughter.alpha(),
      daughter.y(), daughter.z(), daughter.snp(), daughter.tgl(),
      daughter.signed1Pt(),
      // TracksCovIU
      daughter.sigmaY(), daughter.sigmaZ(), daughter.sigmaSnp(), daughter.sigmaTgl(), daughter.sigma1Pt(),
      daughter.rhoZY(), daughter.rhoSnpY(), daughter.rhoSnpZ(), daughter.rhoTglY(), daughter.rhoTglZ(),
      daughter.rhoTglSnp(), daughter.rho1PtY(), daughter.rho1PtZ(), daughter.rho1PtSnp(), daughter.rho1PtTgl(),
      // TracksExtra
      daughter.tpcInnerParam(), daughter.flags(), daughter.itsClusterSizes(),
      daughter.tpcNClsFindable(), daughter.tpcNClsFindableMinusFound(), daughter.tpcNClsFindableMinusCrossedRows(),
      daughter.trdPattern(), daughter.tpcChi2NCl(), daughter.tofChi2(),
      daughter.tpcSignal(), daughter.tofExpMom(),
      // PID
      daughter.tpcNSigmaPr(), daughter.tpcNSigmaPi(), daughter.tpcNSigmaDe(),
      tofNSigmaTrack);
  }

  //------------------------------------------------------------------
  // function to fit KFParticle 3body vertex
  template <typename TKFParticle>
  bool fit3bodyVertex(TKFParticle& kfpProton, TKFParticle& kfpPion, TKFParticle& kfpDeuteron, TKFParticle& KFHt)
  {
    // Construct 3body vertex
    int nDaughters3body = 3;
    const KFParticle* Daughters3body[3] = {&kfpProton, &kfpPion, &kfpDeuteron};
    KFHt.SetConstructMethod(2);
    try {
      KFHt.Construct(Daughters3body, nDaughters3body);
    } catch (std::runtime_error& e) {
      LOG(debug) << "Failed to create Hyper triton 3-body vertex." << e.what();
      return false;
    }
    LOG(debug) << "Hypertriton vertex constructed.";
    return true;
  }

  void process(ColwithEvTimesMultsCents const& collisions, TrackExtPIDIUwithEvTimes const&, aod::Decay3Bodys const& decay3bodys, aod::Tracked3Bodys const& tracked3bodys, aod::BCsWithTimestamps const&)
  {
    std::vector<bool> isTriggeredCollision(collisions.size(), false);

    int lastRunNumber = -1; // RunNumber of last collision, used for zorro counting
    // Event counting
    for (const auto& collision : collisions) {

      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      if (bc.runNumber() != lastRunNumber) {
        initZorroBC(bc);
        lastRunNumber = bc.runNumber(); // Update the last run number
      }

      // all events
      registry.fill(HIST("hEventCounter"), 0.5);

      // ITS ROF boarder cut if not disabled
      if (!collision.selection_bit(aod::evsel::kNoITSROFrameBorder) && !disableITSROFCut) {
        continue;
      }

      // Zorro event counting
      bool isZorroSelected = false;
      if (cfgSkimmedProcessing) {
        isZorroSelected = zorro.isSelected(bc.globalBC());
        if (isZorroSelected) {
          registry.fill(HIST("hEventCounterZorro"), 0.5);
          isTriggeredCollision[collision.globalIndex()] = true;
        }
      }

      // event selection
      if (!collision.selection_bit(aod::evsel::kIsTriggerTVX) || !collision.selection_bit(aod::evsel::kNoTimeFrameBorder) || (collision.posZ() >= 10.0f || collision.posZ() <= -10.0f)) {
        continue;
      }

      // selected events
      registry.fill(HIST("hEventCounter"), 1.5);
      registry.fill(HIST("hAllSelEventsVtxZ"), collision.posZ());

      if (cfgSkimmedProcessing && isZorroSelected) {
        registry.fill(HIST("hEventCounterZorro"), 1.5);
      }
    }

    int lastCollisionID = -1; // collisionId of last analysed decay3body. Table is sorted.

    // get tracked cluster size info
    fTrackedClSizeVector.clear();
    fTrackedClSizeVector.resize(decay3bodys.size(), 0);
    for (const auto& tvtx3body : tracked3bodys) {
      fTrackedClSizeVector[tvtx3body.decay3BodyId()] = tvtx3body.itsClsSize();
    }

    // Create reduced table
    for (const auto& d3body : decay3bodys) {

      auto collision = d3body.template collision_as<ColwithEvTimesMultsCents>();

      // event selection
      if (!collision.selection_bit(aod::evsel::kNoITSROFrameBorder) && !disableITSROFCut) { // ITS ROF boarder cut if not disabled
        continue;
      }
      if (!collision.selection_bit(aod::evsel::kIsTriggerTVX) || !collision.selection_bit(aod::evsel::kNoTimeFrameBorder) || (collision.posZ() >= 10.0f || collision.posZ() <= -10.0f)) {
        continue;
      }

      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      if (cfgSkimmedProcessing && cfgOnlyKeepInterestedTrigger && !isTriggeredCollision[collision.globalIndex()]) {
        continue;
      }

      // Save the collision
      if (collision.globalIndex() != lastCollisionID) {
        int runNumber = bc.runNumber();
        reducedCollisions(
          collision.posX(), collision.posY(), collision.posZ(),
          collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ(),
          collision.flags(), collision.chi2(), collision.numContrib(),
          collision.collisionTime(), collision.collisionTimeRes(),
          runNumber);
        reducedPVMults(collision.multNTracksPV());
        reducedCentFT0Cs(collision.centFT0C());

        lastCollisionID = collision.globalIndex();
      }

      // Precompute collision index
      const auto collisionIndex = reducedCollisions.lastIndex();

      // Save daughter tracks
      const auto daughter0 = d3body.template track0_as<TrackExtPIDIUwithEvTimes>();
      const auto daughter1 = d3body.template track1_as<TrackExtPIDIUwithEvTimes>();
      const auto daughter2 = d3body.template track2_as<TrackExtPIDIUwithEvTimes>();

      // TOF PID of bachelor must be calcualted here
      // ----------------------------------------------
      auto originalcol = daughter2.template collision_as<ColwithEvTimesMultsCents>();
      double tofNSigmaBach = bachelorTOFPID.GetTOFNSigma(mRespParamsV3, daughter2, originalcol, collision);
      // ----------------------------------------------

      // -------- save reduced track table with decay3body daughters ----------
      fillTrackTable(daughter0, -999, collisionIndex);
      fillTrackTable(daughter1, -999, collisionIndex);
      fillTrackTable(daughter2, tofNSigmaBach, collisionIndex);

      // -------- save reduced decay3body table --------
      const auto trackStartIndex = reducedFullTracksPIDIU.lastIndex();
      reducedDecay3Bodys(collisionIndex, trackStartIndex - 2, trackStartIndex - 1, trackStartIndex);

      // -------- get decay3body info with KF --------
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
      // fit 3body vertex and caclulate radius, phi, z position
      float radius, phi, posZ;
      KFParticle KFHt;
      if (fit3bodyVertex(kfpProton, kfpPion, kfpDeuteron, KFHt)) {
        radius = std::hypot(KFHt.GetX(), KFHt.GetY());
        phi = std::atan2(KFHt.GetPx(), KFHt.GetPy());
        posZ = KFHt.GetZ();
      } else {
        radius = -999.;
        phi = -999.;
        posZ = -999.;
      }

      // -------- get decay3body info with DCA fitter --------
      auto Track0 = getTrackParCov(daughter0);
      auto Track1 = getTrackParCov(daughter1);
      auto Track2 = getTrackParCov(daughter2);
      int n3bodyVtx = fitter3body.process(Track0, Track1, Track2);
      float phiVtx, rVtx, zVtx;
      if (n3bodyVtx == 0) { // discard this pair
        phiVtx = -999.;
        rVtx = -999.;
        zVtx = -999.;
      } else {
        const auto& vtxXYZ = fitter3body.getPCACandidate();

        std::array<float, 3> p0 = {0.}, p1 = {0.}, p2{0.};
        const auto& propagatedTrack0 = fitter3body.getTrack(0);
        const auto& propagatedTrack1 = fitter3body.getTrack(1);
        const auto& propagatedTrack2 = fitter3body.getTrack(2);
        propagatedTrack0.getPxPyPzGlo(p0);
        propagatedTrack1.getPxPyPzGlo(p1);
        propagatedTrack2.getPxPyPzGlo(p2);
        std::array<float, 3> p3B = {p0[0] + p1[0] + p2[0], p0[1] + p1[1] + p2[1], p0[2] + p1[2] + p2[2]};
        phiVtx = std::atan2(p3B[1], p3B[0]);
        rVtx = std::hypot(vtxXYZ[0], vtxXYZ[1]);
        zVtx = vtxXYZ[2];
      }

      // fill 3body info table (KF and DCA fitter info)
      reduced3BodyInfo(radius, phi, posZ, rVtx, phiVtx, zVtx, fTrackedClSizeVector[d3body.globalIndex()]);
    } // end decay3body loop

    registry.fill(HIST("hEventCounter"), 2.5, reducedCollisions.lastIndex() + 1);
  }
};

struct reduced3bodyInitializer {
  Spawns<aod::RedIUTracks> reducedTracksIU;
  void init(InitContext const&) {}
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  metadataInfo.initMetadata(cfgc);
  return WorkflowSpec{
    adaptAnalysisTask<reduced3bodyInitializer>(cfgc),
    adaptAnalysisTask<reduced3bodyCreator>(cfgc),
  };
}
