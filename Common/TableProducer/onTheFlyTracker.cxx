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
// Task to add a table of track parameters propagated to the primary vertex
//

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/trackUtilities.h"
#include "ReconstructionDataFormats/DCA.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "CommonUtils/NameConf.h"
#include "CCDB/CcdbApi.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "DataFormatsCalibration/MeanVertexObject.h"
#include "CommonConstants/GeomConstants.h"

// This task allows for the calculation of aod::collisions and aod::Tracks in a synthetic mannter, 
// smearing MC particles with very configurable settings

using namespace o2;
using namespace o2::framework;

struct onTheFlyTracker {
  Produces<aod::Collisions> collisions;
  Produces<aod::McCollisionLabels> collLabels;
  Produces<aod::StoredTracks> tracksPar;
  Produces<aod::TracksExtension> tracksParExtension;
  Produces<aod::StoredTracksCov> tracksParCov;
  Produces<aod::TracksCovExtension> tracksParCovExtension;
  Produces<aod::McTrackLabels> tracksLabels;
  Produces<aod::TracksDCA> tracksDCA;

  bool fillTracksDCA = false;

  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;
  
  void init(o2::framework::InitContext& initContext)
  {
    // Checking if the tables are requested in the workflow and enabling them
    auto& workflows = initContext.services().get<RunningWorkflowInfo const>();
    for (DeviceSpec const& device : workflows.devices) {
      for (auto const& input : device.inputs) {
        if (input.matcher.binding == "TracksDCA") {
          fillTracksDCA = true;
        }
      }
    }
  }

  template <typename CollType, typename TTrackPar>
  void FillTracksPar(CollType& coll, aod::track::TrackTypeEnum trackType, TTrackPar& trackPar)
  {
    tracksPar(coll.globalIndex(), trackType, trackPar.getX(), trackPar.getAlpha(), trackPar.getY(), trackPar.getZ(), trackPar.getSnp(), trackPar.getTgl(), trackPar.getQ2Pt());
    tracksParExtension(trackPar.getPt(), trackPar.getP(), trackPar.getEta(), trackPar.getPhi());
  }

  void process(aod::McCollision const& mcCollision, aod::McParticles const& mcParticles)
  {
    o2::dataformats::DCA dcaInfoCov;
    o2::dataformats::VertexBase vtx;

    for (auto& mcParticle : mcParticles) {
      // Ideally: o2::delphes::TrackSmearer
      // For now: demonstration

      auto pdg = std::abs(mcParticle.pdgCode());
      if (pdg != 11 && pdg != 13 && pdg != 211 && pdg != 321 && pdg != 2212) continue;

      std::array<float, 21> covV = {0.};
      constexpr int MomInd[6] = {9, 13, 14, 18, 19, 20}; // cov matrix elements for momentum component
      for (int i = 0; i < 6; i++) {
        covV[MomInd[i]] = 1e-3;
        covV[i] = 1e-3;
      }
      o2::track::TrackParCov trackParCov = o2::track::TrackParCov(
        {mcParticle.vx(), mcParticle.vy(), mcParticle.vz()},
        {mcParticle.px(), mcParticle.py(), mcParticle.pz()},
        covV, 0, true);

      aod::track::TrackTypeEnum trackType = aod::track::Track;

      //Fixme: collision index could be changeable
      FillTracksPar(mcCollision, trackType, trackParCov);
      if (fillTracksDCA) {
        tracksDCA(1e-3, 1e-3);
      }
      // TODO do we keep the rho as 0? Also the sigma's are duplicated information
      tracksParCov(std::sqrt(trackParCov.getSigmaY2()), std::sqrt(trackParCov.getSigmaZ2()), std::sqrt(trackParCov.getSigmaSnp2()),
                             std::sqrt(trackParCov.getSigmaTgl2()), std::sqrt(trackParCov.getSigma1Pt2()), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
      tracksParCovExtension(trackParCov.getSigmaY2(), trackParCov.getSigmaZY(), trackParCov.getSigmaZ2(), trackParCov.getSigmaSnpY(),
                            trackParCov.getSigmaSnpZ(), trackParCov.getSigmaSnp2(), trackParCov.getSigmaTglY(), trackParCov.getSigmaTglZ(), trackParCov.getSigmaTglSnp(),
                            trackParCov.getSigmaTgl2(), trackParCov.getSigma1PtY(), trackParCov.getSigma1PtZ(), trackParCov.getSigma1PtSnp(), trackParCov.getSigma1PtTgl(),
                            trackParCov.getSigma1Pt2());
      tracksLabels( mcParticle.globalIndex(), 0 );
    }
    collisions(-1,
              mcCollision.posX(), mcCollision.posY(), mcCollision.posZ(), 
              1e-3, 0.0, 1e-3, 0.0, 0.0, 1e-3,
              0, 1e-3, mcParticles.size(),
              0, 0);
    collLabels(mcCollision.globalIndex(), 0);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<onTheFlyTracker>(cfgc)};
  return workflow;
}
