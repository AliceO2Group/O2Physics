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

#ifndef COMMON_TOOLS_TRACKPROPAGATIONMODULE_H_
#define COMMON_TOOLS_TRACKPROPAGATIONMODULE_H_

#include <cstdlib>
#include <cmath>
#include <array>
#include "Framework/AnalysisDataModel.h"
#include "Framework/Configurable.h"
#include "Framework/HistogramSpec.h"
#include "Common/Tools/TrackTuner.h"

//__________________________________________
// track propagation module
//
// this class is capable of performing the usual track propagation
// and table creation it is a demonstration of core service
// plug-in functionality that could be used to reduce the number of
// heavyweight (e.g. mat-LUT-using, propagating) core services to
// reduce overhead and make it easier to pipeline / parallelize
// bottlenecks in core services

class TrackPropagationModule
{
 public:
  TrackPropagationModule()
  {
    // constructor: set defaults
    minPropagationRadius = o2::constants::geom::XTPCInnerRef + 0.1;
    useTrackTuner = false;
    useTrkPid = false;
    fillTrackTunerTable = false;
    trackTunerConfigSource = o2::aod::track_tuner::InputString;
    std::string trackTunerParams = "debugInfo=0|updateTrackDCAs=1|updateTrackCovMat=1|updateCurvature=0|updateCurvatureIU=0|updatePulls=0|isInputFileFromCCDB=1|pathInputFile=Users/m/mfaggin/test/inputsTrackTuner/PbPb2022|nameInputFile=trackTuner_DataLHC22sPass5_McLHC22l1b2_run529397.root|pathFileQoverPt=Users/h/hsharma/qOverPtGraphs|nameFileQoverPt=D0sigma_Data_removal_itstps_MC_LHC22b1b.root|usePvRefitCorrections=0|qOverPtMC=-1.|qOverPtData=-1.";
  };

  float minPropagationRadius;
  bool useTrackTuner;
  bool useTrkPid;
  bool fillTrackTunerTable;
  int trackTunerConfigSource;
  std::string trackTunerParams;

  // controls behaviour
  bool fillTracksCov = false;
  bool fillTracksDCA = false;
  bool fillTracksDCACov = false;
  int runNumber = -1;

  // pointers to objs needed for operation
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
  const o2::dataformats::MeanVertexObject* mMeanVtx = nullptr;
  o2::parameters::GRPMagField* grpmag = nullptr;
  o2::base::MatLayerCylSet* lut = nullptr;
  std::shared_ptr<TH1> trackTunedTracks;
  TrackTuner trackTunerObj;

  // Running variables
  o2::gpu::gpustd::array<float, 2> mDcaInfo;
  o2::dataformats::DCA mDcaInfoCov;
  o2::dataformats::VertexBase mVtx;
  o2::track::TrackParametrization<float> mTrackPar;
  o2::track::TrackParametrizationWithError<float> mTrackParCov;

  // takes a configurableGroup and reads in necessary configs
  template <typename TConfigurableGroup>
  void readConfiguration(TConfigurableGroup const& cGroup)
  {
    minPropagationRadius = cGroup.minPropagationRadius.value;
    useTrackTuner = cGroup.minPropagationRadius.value;
    useTrkPid = cGroup.useTrkPid.value;
    fillTrackTunerTable = cGroup.fillTrackTunerTable.value;
    trackTunerConfigSource = cGroup.trackTunerConfigSource.value;
    trackTunerParams = cGroup.trackTunerParams.value;
    LOGF(info, "Track propagation module configuration done.");
  }
  template <typename TCCDBLoader>
  void getFromCCDBLoader(TCCDBLoader const& loader)
  {
    lut = loader.lut;
    mMeanVtx = loader.mMeanVtx;
    grpmag = loader.grpmag;
  }

  template <typename TInitContext>
  void init(TInitContext& initContext)
  {
    // Checking if the tables are requested in the workflow and enabling them
    fillTracksCov = isTableRequiredInWorkflow(initContext, "TracksCov");
    fillTracksDCA = isTableRequiredInWorkflow(initContext, "TracksDCA");
    fillTracksDCACov = isTableRequiredInWorkflow(initContext, "TracksDCACov");

    /// TrackTuner initialization
    if (useTrackTuner) {
      std::string outputStringParams = "";
      switch (trackTunerConfigSource) {
        case o2::aod::track_tuner::InputString:
          outputStringParams = trackTunerObj.configParams(trackTunerParams);
          break;
        case o2::aod::track_tuner::Configurables:
          outputStringParams = trackTunerObj.configParams();
          break;

        default:
          LOG(fatal) << "TrackTuner configuration source not defined. Fix it! (Supported options: input string (1); Configurables (2))";
          break;
      }

      trackTunerObj.getDcaGraphs();
    }
  }

  template <typename THistoRegistry, typename TConfigurableGroup>
  void initHistograms(THistoRegistry& registry, TConfigurableGroup& cGroup)
  {
    trackTunedTracks = registry.template add<TH1>("trackTunedTracks", "trackTunedTracks", o2::framework::kTH1D, {{1, 0.5f, 1.5f}});

    // Histograms for track tuner
    o2::framework::AxisSpec axisBinsDCA = {600, -0.15f, 0.15f, "#it{dca}_{xy} (cm)"};
    registry.template add("hDCAxyVsPtRec", "hDCAxyVsPtRec", o2::framework::kTH2F, {axisBinsDCA, cGroup.axisPtQA});
    registry.template add("hDCAxyVsPtMC", "hDCAxyVsPtMC", o2::framework::kTH2F, {axisBinsDCA, cGroup.axisPtQA});
    registry.template add("hDCAzVsPtRec", "hDCAzVsPtRec", o2::framework::kTH2F, {axisBinsDCA, cGroup.axisPtQA});
    registry.template add("hDCAzVsPtMC", "hDCAzVsPtMC", o2::framework::kTH2F, {axisBinsDCA, cGroup.axisPtQA});
  }

  template <bool isMc, typename TTracks, typename TOutputGroup, typename THistoRegistry>
  void fillTrackTables(TTracks const& tracks, TOutputGroup& cursors, THistoRegistry& registry)
  {
    if (fillTracksCov) {
      cursors.tracksParCovPropagated.reserve(tracks.size());
      cursors.tracksParCovExtensionPropagated.reserve(tracks.size());
      if (fillTracksDCACov) {
        cursors.tracksDCACov.reserve(tracks.size());
      }
    } else {
      cursors.tracksParPropagated.reserve(tracks.size());
      cursors.tracksParExtensionPropagated.reserve(tracks.size());
      if (fillTracksDCA) {
        cursors.tracksDCA.reserve(tracks.size());
      }
    }

    for (auto& track : tracks) {
      if (fillTracksCov) {
        if (fillTracksDCA || fillTracksDCACov) {
          mDcaInfoCov.set(999, 999, 999, 999, 999);
        }
        setTrackParCov(track, mTrackParCov);
        if (useTrkPid) {
          mTrackParCov.setPID(track.pidForTracking());
        }
      } else {
        if (fillTracksDCA) {
          mDcaInfo[0] = 999;
          mDcaInfo[1] = 999;
        }
        setTrackPar(track, mTrackPar);
        if (useTrkPid) {
          mTrackPar.setPID(track.pidForTracking());
        }
      }
      // auto trackParCov = getTrackParCov(track);
      o2::aod::track::TrackTypeEnum trackType = (o2::aod::track::TrackTypeEnum)track.trackType();
      // std::array<float, 3> trackPxPyPz;
      // std::array<float, 3> trackPxPyPzTuned = {0.0, 0.0, 0.0};
      double q2OverPtNew = -9999.;
      // Only propagate tracks which have passed the innermost wall of the TPC (e.g. skipping loopers etc). Others fill unpropagated.
      if (track.trackType() == o2::aod::track::TrackIU && track.x() < minPropagationRadius) {
        if (fillTracksCov) {
          if constexpr (isMc) { // checking MC and fillCovMat block begins
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
        }
        bool isPropagationOK = true;

        if (track.has_collision()) {
          auto const& collision = track.collision();
          if (fillTracksCov) {
            mVtx.setPos({collision.posX(), collision.posY(), collision.posZ()});
            mVtx.setCov(collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ());
            isPropagationOK = o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, mTrackParCov, 2.f, matCorr, &mDcaInfoCov);
          } else {
            isPropagationOK = o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, mTrackPar, 2.f, matCorr, &mDcaInfo);
          }
        } else {
          if (fillTracksCov) {
            mVtx.setPos({mMeanVtx->getX(), mMeanVtx->getY(), mMeanVtx->getZ()});
            mVtx.setCov(mMeanVtx->getSigmaX() * mMeanVtx->getSigmaX(), 0.0f, mMeanVtx->getSigmaY() * mMeanVtx->getSigmaY(), 0.0f, 0.0f, mMeanVtx->getSigmaZ() * mMeanVtx->getSigmaZ());
            isPropagationOK = o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, mTrackParCov, 2.f, matCorr, &mDcaInfoCov);
          } else {
            isPropagationOK = o2::base::Propagator::Instance()->propagateToDCABxByBz({mMeanVtx->getX(), mMeanVtx->getY(), mMeanVtx->getZ()}, mTrackPar, 2.f, matCorr, &mDcaInfo);
          }
        }
        if (isPropagationOK) {
          trackType = o2::aod::track::Track;
        }
        // filling some QA histograms for track tuner test purpose
        if (fillTracksCov) {
          if constexpr (isMc) { // checking MC and fillCovMat block begins
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
      }
      // Filling modified Q/Pt values at IU/production point by track tuner in track tuner table
      if (useTrackTuner && fillTrackTunerTable) {
        cursors.tunertable(q2OverPtNew);
      }
      // LOG(info) <<  " trackPropagation (this value filled in tuner table)--> "  << q2OverPtNew;
      if (fillTracksCov) {
        cursors.tracksParPropagated(track.collisionId(), trackType, mTrackParCov.getX(), mTrackParCov.getAlpha(), mTrackParCov.getY(), mTrackParCov.getZ(), mTrackParCov.getSnp(), mTrackParCov.getTgl(), mTrackParCov.getQ2Pt());
        cursors.tracksParExtensionPropagated(mTrackParCov.getPt(), mTrackParCov.getP(), mTrackParCov.getEta(), mTrackParCov.getPhi());
        // TODO do we keep the rho as 0? Also the sigma's are duplicated information
        cursors.tracksParCovPropagated(std::sqrt(mTrackParCov.getSigmaY2()), std::sqrt(mTrackParCov.getSigmaZ2()), std::sqrt(mTrackParCov.getSigmaSnp2()),
                                       std::sqrt(mTrackParCov.getSigmaTgl2()), std::sqrt(mTrackParCov.getSigma1Pt2()), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
        cursors.tracksParCovExtensionPropagated(mTrackParCov.getSigmaY2(), mTrackParCov.getSigmaZY(), mTrackParCov.getSigmaZ2(), mTrackParCov.getSigmaSnpY(),
                                                mTrackParCov.getSigmaSnpZ(), mTrackParCov.getSigmaSnp2(), mTrackParCov.getSigmaTglY(), mTrackParCov.getSigmaTglZ(), mTrackParCov.getSigmaTglSnp(),
                                                mTrackParCov.getSigmaTgl2(), mTrackParCov.getSigma1PtY(), mTrackParCov.getSigma1PtZ(), mTrackParCov.getSigma1PtSnp(), mTrackParCov.getSigma1PtTgl(),
                                                mTrackParCov.getSigma1Pt2());
        if (fillTracksDCA) {
          cursors.tracksDCA(mDcaInfoCov.getY(), mDcaInfoCov.getZ());
        }
        if (fillTracksDCACov) {
          cursors.tracksDCACov(mDcaInfoCov.getSigmaY2(), mDcaInfoCov.getSigmaZ2());
        }
      } else {
        cursors.tracksParPropagated(track.collisionId(), trackType, mTrackPar.getX(), mTrackPar.getAlpha(), mTrackPar.getY(), mTrackPar.getZ(), mTrackPar.getSnp(), mTrackPar.getTgl(), mTrackPar.getQ2Pt());
        cursors.tracksParExtensionPropagated(mTrackPar.getPt(), mTrackPar.getP(), mTrackPar.getEta(), mTrackPar.getPhi());
        if (fillTracksDCA) {
          cursors.tracksDCA(mDcaInfo[0], mDcaInfo[1]);
        }
      }
    }
  }
};

#endif // COMMON_TOOLS_TRACKPROPAGATIONMODULE_H_
