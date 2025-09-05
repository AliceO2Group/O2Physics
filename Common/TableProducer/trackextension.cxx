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
// Task performing basic track selection.
//

#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <DataFormatsParameters/GRPObject.h>
#include <DetectorsBase/MatLayerCylSet.h>
#include <DetectorsBase/Propagator.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/DataTypes.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <array>
#include <cstdlib>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

//****************************************************************************************
/**
 * Produce the more complicated derived track quantities needed for track selection.
 * FIXME: we shall run this only if all other selections are passed to avoid
 * FIXME: computing overhead and errors in calculations
 */
//****************************************************************************************
namespace o2
{
namespace analysis
{
namespace trackextension
{
const char* ccdbpath_lut = "GLO/Param/MatLUT";
const char* ccdbpath_grp = "GLO/GRP/GRP";
const char* ccdburl = "http://alice-ccdb.cern.ch"; /* test  "http://alice-ccdb.cern.ch:8080"; */
} // namespace trackextension
} // namespace analysis
} // namespace o2

struct TrackExtension {
  Produces<aod::TracksDCA> extendedTrackQuantities;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Configurable<bool> compatibilityIU{"compatibilityIU", false, "compatibility option to allow the processing of tracks before the introduction of IU tracks"};

  o2::base::MatLayerCylSet* lut;
  int mRunNumber;
  float mMagField;

  void init(InitContext&)
  {
    using namespace analysis::trackextension;

    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(ccdbpath_lut));

    mRunNumber = 0;
    mMagField = 0.0;
  }

  void processRun2(aod::FullTracks const& tracks, aod::Collisions const&, aod::BCsWithTimestamps const&)
  {
    using namespace analysis::trackextension;

    for (auto& track : tracks) {
      std::array<float, 2> dca{1e10f, 1e10f};
      if (track.has_collision()) {
        if (track.trackType() == o2::aod::track::TrackTypeEnum::Run2Track && track.itsChi2NCl() != 0.f && track.tpcChi2NCl() != 0.f && std::abs(track.x()) < 10.f) {
          auto bc = track.collision_as<aod::Collisions>().bc_as<aod::BCsWithTimestamps>();
          if (mRunNumber != bc.runNumber()) {
            o2::parameters::GRPObject* grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(ccdbpath_grp, bc.timestamp());
            if (grpo != nullptr) {
              mMagField = grpo->getNominalL3Field();
              LOGF(info, "Setting magnetic field to %f kG for run %d", mMagField, bc.runNumber());
            } else {
              LOGF(fatal, "GRP object is not available in CCDB for run=%d at timestamp=%llu", bc.runNumber(), bc.timestamp());
            }
            mRunNumber = bc.runNumber();
          }
          auto trackPar = getTrackPar(track);
          auto const& collision = track.collision();
          trackPar.propagateParamToDCA({collision.posX(), collision.posY(), collision.posZ()}, mMagField, &dca);
        }
      }
      extendedTrackQuantities(dca[0], dca[1]);
    }
  }
  PROCESS_SWITCH(TrackExtension, processRun2, "Process Run2 track extension task", true);

  void processRun3(aod::Tracks const& tracks, aod::Collisions const&, aod::BCsWithTimestamps const&)
  {
    using namespace analysis::trackextension;

    /* it is not clear yet if we will the GRP object per run number */
    /* if that is the case something similar to what has been       */
    /* done for Run2 needs to be implemented                        */
    /* but incorporating here only the GRP object makes appear      */
    /* the double peak in the DCAxy distribution again              */
    /* so probably the whole initialization sequence is needed        */
    /* when a new run (Run3) is processed                           */
    o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;

    for (auto& track : tracks) {
      std::array<float, 2> dca{1e10f, 1e10f};
      if (track.has_collision()) {
        if (((compatibilityIU.value) && track.trackType() == o2::aod::track::TrackTypeEnum::TrackIU) ||
            ((!compatibilityIU.value) && track.trackType() == o2::aod::track::TrackTypeEnum::Track)) {
          auto bc = track.collision_as<aod::Collisions>().bc_as<aod::BCsWithTimestamps>();
          if (mRunNumber != bc.runNumber()) {
            auto grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(ccdbpath_grp, bc.timestamp());
            if (grpo != nullptr) {
              o2::base::Propagator::initFieldFromGRP(grpo);
              o2::base::Propagator::Instance()->setMatLUT(lut);
              LOGF(info, "Setting magnetic field to %d kG for run %d from its GRP CCDB object", grpo->getNominalL3Field(), bc.runNumber());
            } else {
              LOGF(fatal, "GRP object is not available in CCDB for run=%d at timestamp=%llu", bc.runNumber(), bc.timestamp());
            }
            mRunNumber = bc.runNumber();
          }
          auto trackPar = getTrackPar(track);
          auto const& collision = track.collision();
          std::array<float, 2> dcaInfo;
          if (o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackPar, 2.f, matCorr, &dcaInfo)) {
            dca[0] = dcaInfo[0];
            dca[1] = dcaInfo[1];
          }
        }
      }
      extendedTrackQuantities(dca[0], dca[1]);
    }
  }
  PROCESS_SWITCH(TrackExtension, processRun3, "Process Run3 track extension task", false);
};

//****************************************************************************************
/**
 * Workflow definition.
 */
//****************************************************************************************
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<TrackExtension>(cfgc)};
  return workflow;
}
