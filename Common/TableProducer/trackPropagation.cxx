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

#include "Common/Core/TableHelper.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Tools/TrackTuner.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/GeomConstants.h>
#include <DataFormatsCalibration/MeanVertexObject.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DetectorsBase/MatLayerCylSet.h>
#include <DetectorsBase/Propagator.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/ConfigContext.h>
#include <Framework/Configurable.h>
#include <Framework/DataTypes.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/WorkflowSpec.h>
#include <ReconstructionDataFormats/DCA.h>
#include <ReconstructionDataFormats/TrackParametrization.h>
#include <ReconstructionDataFormats/TrackParametrizationWithError.h>
#include <ReconstructionDataFormats/Vertex.h>

#include <TH1.h>

#include <array>
#include <cmath>
#include <string>

// The Run 3 AO2D stores the tracks at the point of innermost update. For a track with ITS this is the innermost (or second innermost)
// ITS layer. For a track without ITS, this is the TPC inner wall or for loopers in the TPC even a radius beyond that.
// In order to use the track parameters, the tracks have to be propagated to the collision vertex which is done by this task.
// The task consumes the TracksIU and TracksCovIU tables and produces Tracks and TracksCov to which then the user analysis can subscribe.
//
// This task is not needed for Run 2 converted data.
// There are two versions of the task (see process flags), one producing also the covariance matrix and the other only the tracks table.

using namespace o2;
using namespace o2::framework;
// using namespace o2::framework::expressions;

struct TrackPropagation {
  Produces<aod::StoredTracks> tracksParPropagated;
  Produces<aod::TracksExtension> tracksParExtensionPropagated;

  Produces<aod::StoredTracksCov> tracksParCovPropagated;
  Produces<aod::TracksCovExtension> tracksParCovExtensionPropagated;

  Produces<aod::TracksDCA> tracksDCA;
  Produces<aod::TracksDCACov> tracksDCACov;

  Produces<aod::TrackTunerTable> tunertable;

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  bool fillTracksDCA = false;
  bool fillTracksDCACov = false;
  int runNumber = -1;

  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;

  const o2::dataformats::MeanVertexObject* mMeanVtx = nullptr;
  o2::parameters::GRPMagField* grpmag = nullptr;
  o2::base::MatLayerCylSet* lut = nullptr;
  TrackTuner trackTunerObj;

  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> mVtxPath{"mVtxPath", "GLO/Calib/MeanVertex", "Path of the mean vertex file"};
  Configurable<float> minPropagationRadius{"minPropagationDistance", o2::constants::geom::XTPCInnerRef + 0.1, "Only tracks which are at a smaller radius will be propagated, defaults to TPC inner wall"};
  // for TrackTuner only (MC smearing)
  Configurable<bool> useTrackTuner{"useTrackTuner", false, "Apply track tuner corrections to MC"};
  Configurable<bool> fillTrackTunerTable{"fillTrackTunerTable", false, "flag to fill track tuner table"};
  Configurable<int> trackTunerConfigSource{"trackTunerConfigSource", aod::track_tuner::InputString, "1: input string; 2: TrackTuner Configurables"};
  Configurable<std::string> trackTunerParams{"trackTunerParams", "debugInfo=0|updateTrackDCAs=1|updateTrackCovMat=1|updateCurvature=0|updateCurvatureIU=0|updatePulls=0|isInputFileFromCCDB=1|pathInputFile=Users/m/mfaggin/test/inputsTrackTuner/PbPb2022|nameInputFile=trackTuner_DataLHC22sPass5_McLHC22l1b2_run529397.root|pathFileQoverPt=Users/h/hsharma/qOverPtGraphs|nameFileQoverPt=D0sigma_Data_removal_itstps_MC_LHC22b1b.root|usePvRefitCorrections=0|qOverPtMC=-1.|qOverPtData=-1.", "TrackTuner parameter initialization (format: <name>=<value>|<name>=<value>)"};
  ConfigurableAxis axisPtQA{"axisPtQA", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "pt axis for QA histograms"};
  OutputObj<TH1D> trackTunedTracks{TH1D("trackTunedTracks", "", 1, 0.5, 1.5), OutputObjHandlingPolicy::AnalysisObject};

  // OutputObj<TH2F> hDCAxyVsPtRec{TH2F("hDCAxyVsPtRec", ";DCAxy;PtRec", 600, -0.15, 0.15, axisPtQA)};
  // OutputObj<TH2F> hDCAxyVsPtMC{TH2F("hDCAxyVsPtMC", ";DCAxy;PtMC", 600, -0.15, 0.15, axisPtQA)};

  using TracksIUWithMc = soa::Join<aod::StoredTracksIU, aod::McTrackLabels, aod::TracksCovIU>;

  HistogramRegistry registry{"registry"};

  void init(o2::framework::InitContext& initContext)
  {
    int nEnabledProcesses = 0;
    if (doprocessStandard) {
      LOG(info) << "Enabling processStandard";
      nEnabledProcesses++;
    }
    if (doprocessCovarianceMc) {
      LOG(info) << "Enabling processCovarianceMc";
      nEnabledProcesses++;
    }

    if (doprocessCovariance) {
      LOG(info) << "Enabling processCovariance";
      nEnabledProcesses++;
    }

    if (doprocessStandardWithPID) {
      LOG(info) << "Enabling processStandardWithPID";
      nEnabledProcesses++;
    }
    if (doprocessCovarianceWithPID) {
      LOG(info) << "Enabling processCovarianceWithPID";
      nEnabledProcesses++;
    }
    if (nEnabledProcesses != 1) {
      LOG(fatal) << "Exactly one process flag must be set to true. Please choose one.";
    }
    // Checking if the tables are requested in the workflow and enabling them
    fillTracksDCA = isTableRequiredInWorkflow(initContext, "TracksDCA");
    fillTracksDCACov = isTableRequiredInWorkflow(initContext, "TracksDCACov");

    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    // Histograms for track tuner
    AxisSpec axisBinsDCA = {600, -0.15f, 0.15f, "#it{dca}_{xy} (cm)"};
    registry.add("hDCAxyVsPtRec", "hDCAxyVsPtRec", kTH2F, {axisBinsDCA, axisPtQA});
    registry.add("hDCAxyVsPtMC", "hDCAxyVsPtMC", kTH2F, {axisBinsDCA, axisPtQA});
    registry.add("hDCAzVsPtRec", "hDCAzVsPtRec", kTH2F, {axisBinsDCA, axisPtQA});
    registry.add("hDCAzVsPtMC", "hDCAzVsPtMC", kTH2F, {axisBinsDCA, axisPtQA});

    /// TrackTuner initialization
    if (useTrackTuner) {
      std::string outputStringParams = "";
      switch (trackTunerConfigSource) {
        case aod::track_tuner::InputString:
          outputStringParams = trackTunerObj.configParams(trackTunerParams);
          break;
        case aod::track_tuner::Configurables:
          outputStringParams = trackTunerObj.configParams();
          break;

        default:
          LOG(fatal) << "TrackTuner configuration source not defined. Fix it! (Supported options: input string (1); Configurables (2))";
          break;
      }

      trackTunerObj.getDcaGraphs();
      trackTunedTracks->SetTitle(outputStringParams.c_str());
      trackTunedTracks->GetXaxis()->SetBinLabel(1, "all tracks");
    }
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (runNumber == bc.runNumber()) {
      return;
    }

    // load matLUT for this timestamp
    if (!lut) {
      LOG(info) << "Loading material look-up table for timestamp: " << bc.timestamp();
      lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->getForTimeStamp<o2::base::MatLayerCylSet>(lutPath, bc.timestamp()));
    } else {
      LOG(info) << "Material look-up table already in place. Not reloading.";
    }

    grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, bc.timestamp());
    LOG(info) << "Setting magnetic field to current " << grpmag->getL3Current() << " A for run " << bc.runNumber() << " from its GRPMagField CCDB object";
    o2::base::Propagator::initFieldFromGRP(grpmag);
    o2::base::Propagator::Instance()->setMatLUT(lut);
    mMeanVtx = ccdb->getForTimeStamp<o2::dataformats::MeanVertexObject>(mVtxPath, bc.timestamp());
    runNumber = bc.runNumber();
  }

  // Running variables
  std::array<float, 2> mDcaInfo;
  o2::dataformats::DCA mDcaInfoCov;
  o2::dataformats::VertexBase mVtx;
  o2::track::TrackParametrization<float> mTrackPar;
  o2::track::TrackParametrizationWithError<float> mTrackParCov;

  template <typename TTrack, typename TParticle, bool isMc, bool fillCovMat = false, bool useTrkPid = false>
  void fillTrackTables(TTrack const& tracks,
                       TParticle const&,
                       aod::Collisions const&,
                       aod::BCsWithTimestamps const& bcs)
  {
    if (bcs.size() == 0) {
      return;
    }
    initCCDB(bcs.begin());

    if constexpr (fillCovMat) {
      tracksParCovPropagated.reserve(tracks.size());
      tracksParCovExtensionPropagated.reserve(tracks.size());
      if (fillTracksDCACov) {
        tracksDCACov.reserve(tracks.size());
      }
    } else {
      tracksParPropagated.reserve(tracks.size());
      tracksParExtensionPropagated.reserve(tracks.size());
      if (fillTracksDCA) {
        tracksDCA.reserve(tracks.size());
      }
    }

    for (auto& track : tracks) {
      if constexpr (fillCovMat) {
        if (fillTracksDCA || fillTracksDCACov) {
          mDcaInfoCov.set(999, 999, 999, 999, 999);
        }
        setTrackParCov(track, mTrackParCov);
        if constexpr (useTrkPid) {
          mTrackParCov.setPID(track.pidForTracking());
        }
      } else {
        if (fillTracksDCA) {
          mDcaInfo[0] = 999;
          mDcaInfo[1] = 999;
        }
        setTrackPar(track, mTrackPar);
        if constexpr (useTrkPid) {
          mTrackPar.setPID(track.pidForTracking());
        }
      }
      // auto trackParCov = getTrackParCov(track);
      aod::track::TrackTypeEnum trackType = (aod::track::TrackTypeEnum)track.trackType();
      // std::array<float, 3> trackPxPyPz;
      // std::array<float, 3> trackPxPyPzTuned = {0.0, 0.0, 0.0};
      double q2OverPtNew = -9999.;
      // Only propagate tracks which have passed the innermost wall of the TPC (e.g. skipping loopers etc). Others fill unpropagated.
      if (track.trackType() == aod::track::TrackIU && track.x() < minPropagationRadius) {
        if constexpr (isMc && fillCovMat) { // checking MC and fillCovMat block begins
          // bool hasMcParticle = track.has_mcParticle();
          if (useTrackTuner) {
            trackTunedTracks->Fill(1); // all tracks
            bool hasMcParticle = track.has_mcParticle();
            if (hasMcParticle) {
              auto mcParticle = track.mcParticle();
              trackTunerObj.tuneTrackParams(mcParticle, mTrackParCov, matCorr, &mDcaInfoCov, trackTunedTracks);
              q2OverPtNew = mTrackParCov.getQ2Pt();
            }
          }
        } // MC and fillCovMat block ends
        bool isPropagationOK = true;

        if (track.has_collision()) {
          auto const& collision = track.collision();
          if constexpr (fillCovMat) {
            mVtx.setPos({collision.posX(), collision.posY(), collision.posZ()});
            mVtx.setCov(collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ());
            isPropagationOK = o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, mTrackParCov, 2.f, matCorr, &mDcaInfoCov);
          } else {
            isPropagationOK = o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, mTrackPar, 2.f, matCorr, &mDcaInfo);
          }
        } else {
          if constexpr (fillCovMat) {
            mVtx.setPos({mMeanVtx->getX(), mMeanVtx->getY(), mMeanVtx->getZ()});
            mVtx.setCov(mMeanVtx->getSigmaX() * mMeanVtx->getSigmaX(), 0.0f, mMeanVtx->getSigmaY() * mMeanVtx->getSigmaY(), 0.0f, 0.0f, mMeanVtx->getSigmaZ() * mMeanVtx->getSigmaZ());
            isPropagationOK = o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, mTrackParCov, 2.f, matCorr, &mDcaInfoCov);
          } else {
            isPropagationOK = o2::base::Propagator::Instance()->propagateToDCABxByBz({mMeanVtx->getX(), mMeanVtx->getY(), mMeanVtx->getZ()}, mTrackPar, 2.f, matCorr, &mDcaInfo);
          }
        }
        if (isPropagationOK) {
          trackType = aod::track::Track;
        }
        // filling some QA histograms for track tuner test purpose
        if constexpr (isMc && fillCovMat) { // checking MC and fillCovMat block begins
          if (track.has_mcParticle() && isPropagationOK) {
            auto mcParticle1 = track.mcParticle();
            // && abs(mcParticle1.pdgCode())==211
            if (mcParticle1.isPhysicalPrimary()) {
              registry.fill(HIST("hDCAxyVsPtRec"), mDcaInfoCov.getY(), mTrackParCov.getPt());
              registry.fill(HIST("hDCAxyVsPtMC"), mDcaInfoCov.getY(), mcParticle1.pt());
              registry.fill(HIST("hDCAzVsPtRec"), mDcaInfoCov.getZ(), mTrackParCov.getPt());
              registry.fill(HIST("hDCAzVsPtMC"), mDcaInfoCov.getZ(), mcParticle1.pt());
            }
          }
        } // MC and fillCovMat block ends
      }
      // Filling modified Q/Pt values at IU/production point by track tuner in track tuner table
      if (useTrackTuner && fillTrackTunerTable) {
        tunertable(q2OverPtNew);
      }
      // LOG(info) <<  " trackPropagation (this value filled in tuner table)--> "  << q2OverPtNew;
      if constexpr (fillCovMat) {
        tracksParPropagated(track.collisionId(), trackType, mTrackParCov.getX(), mTrackParCov.getAlpha(), mTrackParCov.getY(), mTrackParCov.getZ(), mTrackParCov.getSnp(), mTrackParCov.getTgl(), mTrackParCov.getQ2Pt());
        tracksParExtensionPropagated(mTrackParCov.getPt(), mTrackParCov.getP(), mTrackParCov.getEta(), mTrackParCov.getPhi());
        // TODO do we keep the rho as 0? Also the sigma's are duplicated information
        tracksParCovPropagated(std::sqrt(mTrackParCov.getSigmaY2()), std::sqrt(mTrackParCov.getSigmaZ2()), std::sqrt(mTrackParCov.getSigmaSnp2()),
                               std::sqrt(mTrackParCov.getSigmaTgl2()), std::sqrt(mTrackParCov.getSigma1Pt2()), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
        tracksParCovExtensionPropagated(mTrackParCov.getSigmaY2(), mTrackParCov.getSigmaZY(), mTrackParCov.getSigmaZ2(), mTrackParCov.getSigmaSnpY(),
                                        mTrackParCov.getSigmaSnpZ(), mTrackParCov.getSigmaSnp2(), mTrackParCov.getSigmaTglY(), mTrackParCov.getSigmaTglZ(), mTrackParCov.getSigmaTglSnp(),
                                        mTrackParCov.getSigmaTgl2(), mTrackParCov.getSigma1PtY(), mTrackParCov.getSigma1PtZ(), mTrackParCov.getSigma1PtSnp(), mTrackParCov.getSigma1PtTgl(),
                                        mTrackParCov.getSigma1Pt2());
        if (fillTracksDCA) {
          tracksDCA(mDcaInfoCov.getY(), mDcaInfoCov.getZ());
        }
        if (fillTracksDCACov) {
          tracksDCACov(mDcaInfoCov.getSigmaY2(), mDcaInfoCov.getSigmaZ2());
        }
      } else {
        tracksParPropagated(track.collisionId(), trackType, mTrackPar.getX(), mTrackPar.getAlpha(), mTrackPar.getY(), mTrackPar.getZ(), mTrackPar.getSnp(), mTrackPar.getTgl(), mTrackPar.getQ2Pt());
        tracksParExtensionPropagated(mTrackPar.getPt(), mTrackPar.getP(), mTrackPar.getEta(), mTrackPar.getPhi());
        if (fillTracksDCA) {
          tracksDCA(mDcaInfo[0], mDcaInfo[1]);
        }
      }
    }
  }

  void processStandard(aod::StoredTracksIU const& tracks, aod::Collisions const& collisions, aod::BCsWithTimestamps const& bcs)
  {
    fillTrackTables</*TTrack*/ aod::StoredTracksIU, /*Particle*/ aod::StoredTracksIU, /*isMc = */ false, /*fillCovMat =*/false, /*useTrkPid =*/false>(tracks, tracks, collisions, bcs);
  }
  PROCESS_SWITCH(TrackPropagation, processStandard, "Process without covariance", true);

  void processStandardWithPID(soa::Join<aod::StoredTracksIU, aod::TracksExtra> const& tracks, aod::Collisions const& collisions, aod::BCsWithTimestamps const& bcs)
  {
    fillTrackTables</*TTrack*/ soa::Join<aod::StoredTracksIU, aod::TracksExtra>, /*Particle*/ soa::Join<aod::StoredTracksIU, aod::TracksExtra>, /*isMc = */ false, /*fillCovMat =*/false, /*useTrkPid =*/true>(tracks, tracks, collisions, bcs);
  }
  PROCESS_SWITCH(TrackPropagation, processStandardWithPID, "Process without covariance and with PID in tracking", false);

  // -----------------------
  void processCovarianceMc(TracksIUWithMc const& tracks, aod::McParticles const& mcParticles, aod::Collisions const& collisions, aod::BCsWithTimestamps const& bcs)
  {
    // auto table_extension = soa::Extend<TracksIUWithMc, aod::extension::MomX>(tracks);
    fillTrackTables</*TTrack*/ TracksIUWithMc, /*Particle*/ aod::McParticles, /*isMc = */ true, /*fillCovMat =*/true, /*useTrkPid =*/false>(tracks, mcParticles, collisions, bcs);
  }
  PROCESS_SWITCH(TrackPropagation, processCovarianceMc, "Process with covariance on MC", false);

  void processCovariance(soa::Join<aod::StoredTracksIU, aod::TracksCovIU> const& tracks, aod::Collisions const& collisions, aod::BCsWithTimestamps const& bcs)
  {
    fillTrackTables</*TTrack*/ soa::Join<aod::StoredTracksIU, aod::TracksCovIU>, /*Particle*/ soa::Join<aod::StoredTracksIU, aod::TracksCovIU>, /*isMc = */ false, /*fillCovMat =*/true, /*useTrkPid =*/false>(tracks, tracks, collisions, bcs);
  }
  PROCESS_SWITCH(TrackPropagation, processCovariance, "Process with covariance", false);
  // ------------------------

  void processCovarianceWithPID(soa::Join<aod::StoredTracksIU, aod::TracksCovIU, aod::TracksExtra> const& tracks, aod::Collisions const& collisions, aod::BCsWithTimestamps const& bcs)
  {
    fillTrackTables</*TTrack*/ soa::Join<aod::StoredTracksIU, aod::TracksCovIU, aod::TracksExtra>, /*Particle*/ soa::Join<aod::StoredTracksIU, aod::TracksCovIU, aod::TracksExtra>, /*isMc = */ false, /*fillCovMat =*/true, /*useTrkPid =*/false>(tracks, tracks, collisions, bcs);
  }
  PROCESS_SWITCH(TrackPropagation, processCovarianceWithPID, "Process with covariance and with PID in tracking", false);
};

//****************************************************************************************
/**
 * Workflow definition.
 */
//****************************************************************************************
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<TrackPropagation>(cfgc)};
  return workflow;
}
