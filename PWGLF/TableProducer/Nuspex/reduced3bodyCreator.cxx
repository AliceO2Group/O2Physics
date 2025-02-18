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

/// \brief Task to produce reduced AO2Ds for use in the hypertriton 3body reconstruction with the decay3bodybuilder.cxx
/// \author Yuanzhe Wang <yuanzhe.wang@cern.ch>
/// \author Carolina Reetz <c.reetz@cern.ch>

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
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "PWGLF/DataModel/pidTOFGeneric.h"
#include "PWGLF/DataModel/Reduced3BodyTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/Core/PID/PIDTOF.h"
#include "TableHelper.h"

#include "EventFiltering/Zorro.h"
#include "EventFiltering/ZorroSummary.h"

#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using FullTracksExtIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU>;
using FullTracksExtPIDIU = soa::Join<FullTracksExtIU, aod::pidTPCFullPr, aod::pidTPCFullPi, aod::pidTPCFullDe>;

using ColwithEvTimes = o2::soa::Join<aod::Collisions, aod::EvSels, aod::EvTimeTOFFT0>;
using ColwithEvTimesMultsCents = o2::soa::Join<ColwithEvTimes, aod::PVMults, aod::CentFT0Cs>;
using TrackExtIUwithEvTimes = soa::Join<FullTracksExtIU, aod::EvTimeTOFFT0ForTrack>;
using TrackExtPIDIUwithEvTimes = soa::Join<FullTracksExtPIDIU, aod::EvTimeTOFFT0ForTrack>;

struct reduced3bodyCreator {

  Produces<aod::ReducedCollisions> reducedCollisions;
  Produces<aod::ReducedPVMults> reducedPVMults;
  Produces<aod::ReducedCentFT0Cs> reducedCentFTOCs;
  Produces<aod::ReducedDecay3Bodys> reducedDecay3Bodys;
  Produces<aod::StoredReducedTracksIU> reducedFullTracksPIDIU;

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Zorro zorro;
  OutputObj<ZorroSummary> zorroSummary{"zorroSummary"};

  o2::aod::pidtofgeneric::TofPidNewCollision<TrackExtPIDIUwithEvTimes::iterator> bachelorTOFPID;

  std::vector<TrackExtPIDIUwithEvTimes::iterator> daughterTracks;

  Configurable<bool> event_sel8_selection{"event_sel8_selection", true, "event selection count post sel8 cut"};
  Configurable<bool> mc_event_selection{"mc_event_selection", true, "mc event selection count post kIsTriggerTVX and kNoTimeFrameBorder"};
  Configurable<bool> event_posZ_selection{"event_posZ_selection", true, "event selection count post poZ cut"};
  // CCDB options
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
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

  int mRunNumber;
  o2::pid::tof::TOFResoParamsV2 mRespParamsV2;

  HistogramRegistry registry{"registry", {}};

  void init(InitContext&)
  {
    mRunNumber = 0;
    zorroSummary.setObject(zorro.getZorroSummary());
    bachelorTOFPID.SetPidType(o2::track::PID::Deuteron);

    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    registry.add("hAllSelEventsVtxZ", "hAllSelEventsVtxZ", HistType::kTH1F, {{500, -15.0f, 15.0f, "PV Z (cm)"}});

    auto hEventCounter = registry.add<TH1>("hEventCounter", "hEventCounter", HistType::kTH1F, {{4, 0.0f, 4.0f}});
    hEventCounter->GetXaxis()->SetBinLabel(1, "total");
    hEventCounter->GetXaxis()->SetBinLabel(2, "sel8");
    hEventCounter->GetXaxis()->SetBinLabel(3, "vertexZ");
    hEventCounter->GetXaxis()->SetBinLabel(4, "reduced");
    hEventCounter->LabelsOption("v");

    auto hEventCounterZorro = registry.add<TH1>("hEventCounterZorro", "hEventCounterZorro", HistType::kTH1D, {{2, -0.5, 1.5}});
    hEventCounterZorro->GetXaxis()->SetBinLabel(1, "Zorro before evsel");
    hEventCounterZorro->GetXaxis()->SetBinLabel(2, "Zorro after evsel");
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    // In case override, don't proceed, please - no CCDB access required
    if (mRunNumber == bc.runNumber()) {
      return;
    }

    mRunNumber = bc.runNumber();
    if (cfgSkimmedProcessing) {
      zorro.initCCDB(ccdb.service, bc.runNumber(), bc.timestamp(), "fH3L3Body");
      zorro.populateHistRegistry(registry, bc.runNumber());
    }

    // Initial TOF PID Paras, copied from PIDTOF.h
    timestamp.value = bc.timestamp();
    ccdb->setTimestamp(timestamp.value);
    // Not later than now objects
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    // TODO: implement the automatic pass name detection from metadata
    if (passName.value == "") {
      passName.value = "unanchored"; // temporary default
      LOG(warning) << "Passed autodetect mode for pass, not implemented yet, waiting for metadata. Taking '" << passName.value << "'";
    }
    LOG(info) << "Using parameter collection, starting from pass '" << passName.value << "'";

    const std::string fname = paramFileName.value;
    if (!fname.empty()) { // Loading the parametrization from file
      LOG(info) << "Loading exp. sigma parametrization from file " << fname << ", using param: " << parametrizationPath.value;
      if (1) {
        o2::tof::ParameterCollection paramCollection;
        paramCollection.loadParamFromFile(fname, parametrizationPath.value);
        LOG(info) << "+++ Loaded parameter collection from file +++";
        if (!paramCollection.retrieveParameters(mRespParamsV2, passName.value)) {
          if (fatalOnPassNotAvailable) {
            LOGF(fatal, "Pass '%s' not available in the retrieved CCDB object", passName.value.data());
          } else {
            LOGF(warning, "Pass '%s' not available in the retrieved CCDB object", passName.value.data());
          }
        } else {
          mRespParamsV2.setShiftParameters(paramCollection.getPars(passName.value));
          mRespParamsV2.printShiftParameters();
        }
      } else {
        mRespParamsV2.loadParamFromFile(fname.data(), parametrizationPath.value);
      }
    } else if (loadResponseFromCCDB) { // Loading it from CCDB
      LOG(info) << "Loading exp. sigma parametrization from CCDB, using path: " << parametrizationPath.value << " for timestamp " << timestamp.value;
      o2::tof::ParameterCollection* paramCollection = ccdb->getForTimeStamp<o2::tof::ParameterCollection>(parametrizationPath.value, timestamp.value);
      paramCollection->print();
      if (!paramCollection->retrieveParameters(mRespParamsV2, passName.value)) { // Attempt at loading the parameters with the pass defined
        if (fatalOnPassNotAvailable) {
          LOGF(fatal, "Pass '%s' not available in the retrieved CCDB object", passName.value.data());
        } else {
          LOGF(warning, "Pass '%s' not available in the retrieved CCDB object", passName.value.data());
        }
      } else { // Pass is available, load non standard parameters
        mRespParamsV2.setShiftParameters(paramCollection->getPars(passName.value));
        mRespParamsV2.printShiftParameters();
      }
    }
    mRespParamsV2.print();
    if (timeShiftCCDBPath.value != "") {
      if (timeShiftCCDBPath.value.find(".root") != std::string::npos) {
        mRespParamsV2.setTimeShiftParameters(timeShiftCCDBPath.value, "gmean_Pos", true);
        mRespParamsV2.setTimeShiftParameters(timeShiftCCDBPath.value, "gmean_Neg", false);
      } else {
        mRespParamsV2.setTimeShiftParameters(ccdb->getForTimeStamp<TGraph>(Form("%s/pos", timeShiftCCDBPath.value.c_str()), timestamp.value), true);
        mRespParamsV2.setTimeShiftParameters(ccdb->getForTimeStamp<TGraph>(Form("%s/neg", timeShiftCCDBPath.value.c_str()), timestamp.value), false);
      }
    }

    bachelorTOFPID.SetParams(mRespParamsV2);
  }

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

  void process(ColwithEvTimesMultsCents const& collisions, TrackExtPIDIUwithEvTimes const&, aod::Decay3Bodys const& decay3bodys, aod::BCsWithTimestamps const&)
  {

    int lastCollisionID = -1; // collisionId of last analysed decay3body. Table is sorted.

    // Event counting
    for (const auto& collision : collisions) {
      // Zorro event counting
      bool isZorroSelected = false;
      if (cfgSkimmedProcessing) {
        isZorroSelected = zorro.isSelected(collision.bc_as<aod::BCsWithTimestamps>().globalBC());
        if (isZorroSelected) {
          registry.fill(HIST("hEventCounterZorro"), 0.5);
        }
      }

      // Event selection
      registry.fill(HIST("hEventCounter"), 0.5);
      if (event_sel8_selection && !collision.sel8()) {
        continue;
      }
      registry.fill(HIST("hEventCounter"), 1.5);
      if (event_posZ_selection && (collision.posZ() >= 10.0f || collision.posZ() <= -10.0f)) { // 10cm
        continue;
      }
      registry.fill(HIST("hEventCounter"), 2.5);
      registry.fill(HIST("hAllSelEventsVtxZ"), collision.posZ());

      if (cfgSkimmedProcessing && isZorroSelected) {
        registry.fill(HIST("hEventCounterZorro"), 1.5);
      }
    }

    // Creat reduced table
    for (const auto& d3body : decay3bodys) {

      daughterTracks.clear();

      auto collision = d3body.template collision_as<ColwithEvTimesMultsCents>();

      if (event_sel8_selection && !collision.sel8()) {
        continue;
      }
      if (event_posZ_selection && (collision.posZ() >= 10.0f || collision.posZ() <= -10.0f)) { // 10cm
        continue;
      }

      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

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
        reducedCentFTOCs(collision.centFT0C());

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
      double tofNSigmaBach = bachelorTOFPID.GetTOFNSigma(daughter2, originalcol, collision);
      // ----------------------------------------------

      // save reduced track table with decay3body daughters
      fillTrackTable(daughter0, -999, collisionIndex);
      fillTrackTable(daughter1, -999, collisionIndex);
      fillTrackTable(daughter2, tofNSigmaBach, collisionIndex);

      // save reduced decay3body table
      const auto trackStartIndex = reducedFullTracksPIDIU.lastIndex();
      reducedDecay3Bodys(collisionIndex, trackStartIndex - 2, trackStartIndex - 1, trackStartIndex);
    }

    registry.fill(HIST("hEventCounter"), 3.5, reducedCollisions.lastIndex() + 1);
  }
};

struct reduced3bodyInitializer {
  Spawns<aod::ReducedTracksIU> reducedTracksIU;
  void init(InitContext const&) {}
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<reduced3bodyInitializer>(cfgc),
    adaptAnalysisTask<reduced3bodyCreator>(cfgc),
  };
}
