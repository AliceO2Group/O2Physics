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
#include "Common/Core/trackUtilities.h"
#include "ReconstructionDataFormats/DCA.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsCommonDataFormats/NameConf.h"

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
struct TrackExtensionTask {
  Configurable<int> cfgDcaMethod{"dcamethod", 1, "Method to estimate the track DCA: 0 = crude, 1 = minimum, 2 = rigorous. Default minimum"};

  Produces<aod::TracksExtended> extendedTrackQuantities;

  void process(aod::FullTracks const& tracks, aod::Collisions const&)
  {
    o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
    if ((cfgDcaMethod == 1) or (cfgDcaMethod == 2)) {
      if (!o2::base::GeometryManager::isGeometryLoaded()) {
        o2::base::GeometryManager::loadGeometry();
        o2::base::Propagator::initFieldFromGRP();
        auto matLUTFile = o2::base::NameConf::getMatLUTFileName();
        if (o2::utils::Str::pathExists(matLUTFile)) {
          auto* lut = o2::base::MatLayerCylSet::loadFromFile(matLUTFile);
          o2::base::Propagator::Instance()->setMatLUT(lut);
        }
      }
    }
    for (auto& track : tracks) {

      std::array<float, 2> dca{1e10f, 1e10f};
      if (track.has_collision()) {
        if ((track.trackType() == o2::aod::track::TrackTypeEnum::Track) ||
            (track.trackType() == o2::aod::track::TrackTypeEnum::Run2Track && track.itsChi2NCl() != 0.f && track.tpcChi2NCl() != 0.f && std::abs(track.x()) < 10.f)) {
          auto trackPar = getTrackPar(track);
          auto const& collision = track.collision();
          if (cfgDcaMethod == 1) {
            trackPar.propagateParamToDCA({collision.posX(), collision.posY(), collision.posZ()}, o2::base::Propagator::Instance()->getNominalBz(), &dca);
          } else if (cfgDcaMethod == 2) {
            gpu::gpustd::array<float, 2> dcaInfo;
            if (o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackPar, 2.f, matCorr, &dcaInfo)) {
              dca[0] = dcaInfo[0];
              dca[1] = dcaInfo[1];
            }
          } else {
            float magField = 5.0; // in kG (FIXME: get this from CCDB)
            trackPar.propagateParamToDCA({collision.posX(), collision.posY(), collision.posZ()}, magField, &dca);
          }
        }
      }
      extendedTrackQuantities(dca[0], dca[1]);

      // TODO: add realtive pt resolution sigma(pt)/pt \approx pt * sigma(1/pt)
      // TODO: add geometrical length / fiducial volume
    }
  }
};

//****************************************************************************************
/**
 * Workflow definition.
 */
//****************************************************************************************
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<TrackExtensionTask>(cfgc)};
  return workflow;
}