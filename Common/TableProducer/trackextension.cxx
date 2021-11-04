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

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/Core/trackUtilities.h"
#include "ReconstructionDataFormats/DCA.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsCommonDataFormats/NameConf.h"
#include "DataFormatsParameters/GRPObject.h"
#include <CCDB/BasicCCDBManager.h>

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
constexpr long run3grp_timestamp = (1619781650000 + 1619781529000) / 2;
const char* ccdbpath_lut = "GLO/Param/MatLUT";
const char* ccdbpath_geo = "GLO/Config/Geometry";
const char* ccdbpath_grp = "GLO/GRP/GRP";
const char* ccdburl = "https://alice-ccdb.cern.ch"; /* test  "http://alice-ccdb.cern.ch:8080"; */
} // namespace trackextension
} // namespace analysis
} // namespace o2

struct TrackExtension {
  Produces<aod::TracksExtended> extendedTrackQuantities;
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  int mRunNumber;
  float mMagField;

  void init(InitContext& context)
  {
    using namespace analysis::trackextension;

    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    auto lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(ccdbpath_lut));

    if (!o2::base::GeometryManager::isGeometryLoaded()) {
      auto* gm = ccdb->get<TGeoManager>(ccdbpath_geo);
      /* it seems this is needed at this level for the material LUT to work properly */
      /* but what happens if the run changes while doing the processing?             */
      o2::parameters::GRPObject* grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(ccdbpath_grp, analysis::trackextension::run3grp_timestamp);
      o2::base::Propagator::initFieldFromGRP(grpo);
      o2::base::Propagator::Instance()->setMatLUT(lut);
    }
    mRunNumber = 0;
    mMagField = 0.0f;
  }

  void processRun2(aod::FullTracks const& tracks, aod::Collisions const&, aod::BCsWithTimestamps const&)
  {
    using namespace analysis::trackextension;

    int lastCollId = -1;
    for (auto& track : tracks) {
      std::array<float, 2> dca{1e10f, 1e10f};
      if (track.has_collision()) {
        if (track.trackType() == o2::aod::track::TrackTypeEnum::Run2Track && track.itsChi2NCl() != 0.f && track.tpcChi2NCl() != 0.f && std::abs(track.x()) < 10.f) {
          if (lastCollId != track.collisionId()) {
            auto bc = track.collision_as<aod::Collisions>().bc_as<aod::BCsWithTimestamps>();
            if (mRunNumber != bc.runNumber()) {
              o2::parameters::GRPObject* grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(ccdbpath_grp, bc.timestamp());
              if (grpo != nullptr) {
                float l3current = grpo->getL3Current();
                mMagField = l3current / 30000.0f * 5.0f;
                LOGF(info, "Setting magnetic field to %f from an L3 current %f for run %d", mMagField, l3current, bc.runNumber());
              } else {
                LOGF(fatal, "GRP object is not available in CCDB for run=%d at timestamp=%llu", bc.runNumber(), bc.timestamp());
              }
              mRunNumber = bc.runNumber();
            }
            lastCollId = track.collisionId();
          }
          auto trackPar = getTrackPar(track);
          auto const& collision = track.collision();
          trackPar.propagateParamToDCA({collision.posX(), collision.posY(), collision.posZ()}, mMagField, &dca);
        }
      }
      extendedTrackQuantities(dca[0], dca[1]);

      // TODO: add realtive pt resolution sigma(pt)/pt \approx pt * sigma(1/pt)
      // TODO: add geometrical length / fiducial volume
    }
  }
  PROCESS_SWITCH(TrackExtension, processRun2, "Process Run2 track extension task", true);

  void processRun3(aod::FullTracks const& tracks, aod::Collisions const&)
  {
    using namespace analysis::trackextension;

    /* it is not clear yet if we will the GRP object per run number */
    /* if that is the case something similar to what has been       */
    /* done for Run2 needs to be implemented                        */
    /* but incorporating here only the GRP object makes appear      */
    /* the double peak in the DCAxy distribution again              */
    /* so probaly the whole initalization sequence is needed        */
    /* when a new run (Run3) is processed                           */
    o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;

    for (auto& track : tracks) {
      std::array<float, 2> dca{1e10f, 1e10f};
      if (track.has_collision()) {
        if (track.trackType() == o2::aod::track::TrackTypeEnum::Track) {
          auto trackPar = getTrackPar(track);
          auto const& collision = track.collision();
          gpu::gpustd::array<float, 2> dcaInfo;
          if (o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackPar, 2.f, matCorr, &dcaInfo)) {
            dca[0] = dcaInfo[0];
            dca[1] = dcaInfo[1];
          }
        }
      }
      extendedTrackQuantities(dca[0], dca[1]);

      // TODO: add realtive pt resolution sigma(pt)/pt \approx pt * sigma(1/pt)
      // TODO: add geometrical length / fiducial volume
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
