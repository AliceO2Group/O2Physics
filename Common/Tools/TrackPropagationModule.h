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

/// \file TrackPropagationModule.h
/// \brief track propagation module functionality to be used in core services
/// \author ALICE

#ifndef COMMON_TOOLS_TRACKPROPAGATIONMODULE_H_
#define COMMON_TOOLS_TRACKPROPAGATIONMODULE_H_

#include "Common/Core/TableHelper.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Tools/TrackTuner.h"

#include <CommonConstants/GeomConstants.h>
#include <DetectorsBase/Propagator.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/Configurable.h>
#include <Framework/DataTypes.h>
#include <Framework/DeviceSpec.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/Logger.h>
#include <ReconstructionDataFormats/DCA.h>
#include <ReconstructionDataFormats/TrackParametrization.h>
#include <ReconstructionDataFormats/TrackParametrizationWithError.h>

#include <TH1.h>
#include <TH2.h>

#include <array>
#include <cmath>
#include <cstdlib>
#include <memory>
#include <string>

//__________________________________________
// track propagation module
//
// this class is capable of performing the usual track propagation
// and table creation it is a demonstration of core service
// plug-in functionality that could be used to reduce the number of
// heavyweight (e.g. mat-LUT-using, propagating) core services to
// reduce overhead and make it easier to pipeline / parallelize
// bottlenecks in core services

namespace o2
{
namespace common
{

struct TrackPropagationProducts : o2::framework::ProducesGroup {
  o2::framework::Produces<aod::StoredTracks> tracksParPropagated;
  o2::framework::Produces<aod::TracksExtension> tracksParExtensionPropagated;
  o2::framework::Produces<aod::StoredTracksCov> tracksParCovPropagated;
  o2::framework::Produces<aod::TracksCovExtension> tracksParCovExtensionPropagated;
  o2::framework::Produces<aod::TracksDCA> tracksDCA;
  o2::framework::Produces<aod::TracksDCACov> tracksDCACov;
  o2::framework::Produces<aod::TrackTunerTable> tunertable;
};

struct TrackPropagationConfigurables : o2::framework::ConfigurableGroup {
  std::string prefix = "trackPropagation";
  o2::framework::Configurable<float> minPropagationRadius{"minPropagationDistance", o2::constants::geom::XTPCInnerRef + 0.1, "Only tracks which are at a smaller radius will be propagated, defaults to TPC inner wall"};
  // for TrackTuner only (MC smearing)
  o2::framework::Configurable<bool> useTrackTuner{"useTrackTuner", false, "Apply track tuner corrections to MC"};
  o2::framework::Configurable<bool> useTrkPid{"useTrkPid", false, "use pid in tracking"};
  o2::framework::Configurable<bool> fillTrackTunerTable{"fillTrackTunerTable", false, "flag to fill track tuner table"};
  o2::framework::Configurable<int> trackTunerConfigSource{"trackTunerConfigSource", aod::track_tuner::InputString, "1: input string; 2: TrackTuner Configurables"};
  o2::framework::Configurable<std::string> trackTunerParams{"trackTunerParams", "debugInfo=0|updateTrackDCAs=1|updateTrackCovMat=1|updateCurvature=0|updateCurvatureIU=0|updatePulls=0|isInputFileFromCCDB=1|pathInputFile=Users/m/mfaggin/test/inputsTrackTuner/PbPb2022|nameInputFile=trackTuner_DataLHC22sPass5_McLHC22l1b2_run529397.root|pathFileQoverPt=Users/h/hsharma/qOverPtGraphs|nameFileQoverPt=D0sigma_Data_removal_itstps_MC_LHC22b1b.root|usePvRefitCorrections=0|qOverPtMC=-1.|qOverPtData=-1.", "TrackTuner parameter initialization (format: <name>=<value>|<name>=<value>)"};
  o2::framework::ConfigurableAxis axisPtQA{"axisPtQA", {o2::framework::VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "pt axis for QA histograms"};
};

class TrackPropagationModule
{
 public:
  TrackPropagationModule()
  {
    // constructor
  }

  // controls behaviour
  bool fillTracks = false;
  bool fillTracksCov = false;
  bool fillTracksDCA = false;
  bool fillTracksDCACov = false;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;

  // pointers to objs needed for operation
  std::shared_ptr<TH1> trackTunedTracks;

  // Running variables
  std::array<float, 2> mDcaInfo{};
  o2::dataformats::DCA mDcaInfoCov;
  o2::dataformats::VertexBase mVtx;
  o2::track::TrackParametrization<float> mTrackPar;
  o2::track::TrackParametrizationWithError<float> mTrackParCov;
  bool autoDetectDcaCalib = false; // track tuner setting

  template <typename TConfigurableGroup, typename TInitContext, typename THistoRegistry>
  void init(TConfigurableGroup const& cGroup, TrackTuner& trackTunerObj, THistoRegistry& registry, TInitContext& initContext)
  {
    // Checking if the tables are requested in the workflow and enabling them
    fillTracks = isTableRequiredInWorkflow(initContext, "Tracks");
    fillTracksCov = isTableRequiredInWorkflow(initContext, "TracksCov");
    fillTracksDCA = isTableRequiredInWorkflow(initContext, "TracksDCA");
    fillTracksDCACov = isTableRequiredInWorkflow(initContext, "TracksDCACov");

    // enable Tracks in case Tracks have been requested
    if (fillTracksDCA && !fillTracks) {
      LOGF(info, "******************************************************************");
      LOGF(info, " There is no task subscribed to Tracks, but I have detected a");
      LOGF(info, " subscription to TracksDCA. Now enabling tracks as algorithmic");
      LOGF(info, " dependency. Note: please be sure this is intentional! For");
      LOGF(info, " secondary analyses, the proper DCA to test against is the DCA");
      LOGF(info, " that the V0 or Cascade is assigned to and not necessarily the");
      LOGF(info, " the one that the Track is assigned to (if any). ");
      LOGF(info, "******************************************************************");
      fillTracks = true;
    }

    if (!fillTracks) {
      LOGF(info, "Track propagation to PV not required. Suppressing all further processing and logs.");
    }

    LOGF(info, " Track propagation table detection results:");
    if (fillTracks) {
      LOGF(info, " ---> Will generate Tracks table.");
    }
    if (fillTracksCov) {
      LOGF(info, " ---> Will generate TracksCov table.");
    }
    if (fillTracksDCA) {
      LOGF(info, " ---> Will generate TracksDCA table.");
    }
    if (fillTracksDCACov) {
      LOGF(info, " ---> Will generate TracksDCACov table.");
    }
    if (fillTracksCov) {
      LOGF(info, "**************************************************************");
      LOGF(info, " Warning: TracksCov has been requested due to a subscription!");
      LOGF(info, " Please be mindful that generating track covariances requires");
      LOGF(info, " a significant extra amount of CPU and memory. If not strictly");
      LOGF(info, " necessary, requesting TracksCov should be avoided to save");
      LOGF(info, " these additional resouces.");
      LOGF(info, "**************************************************************");
    }

    /// TrackTuner initialization
    std::string outputStringParams = "";
    if (cGroup.useTrackTuner.value) {
      switch (cGroup.trackTunerConfigSource.value) {
        case o2::aod::track_tuner::InputString:
          outputStringParams = trackTunerObj.configParams(cGroup.trackTunerParams.value);
          break;
        case o2::aod::track_tuner::Configurables:
          outputStringParams = trackTunerObj.configParams();
          break;

        default:
          LOG(fatal) << "TrackTuner configuration source not defined. Fix it! (Supported options: input string (1); Configurables (2))";
          break;
      }

      /// read the track tuner instance configurations,
      /// to understand whether the TrackTuner::getDcaGraphs function can be called here (input path from string/configurables)
      /// or inside the process function, to "auto-detect" the input file based on the run number
      const auto& workflows = initContext.services().template get<o2::framework::RunningWorkflowInfo const>();
      for (const o2::framework::DeviceSpec& device : workflows.devices) { /// loop over devices
        if (device.name == "propagation-service") {
          // loop over the options
          // to find the value of TrackTuner::autoDetectDcaCalib
          for (const auto& option : device.options) { /// loop over options
            if (option.name == "trackTuner.autoDetectDcaCalib") {
              // found it!
              autoDetectDcaCalib = option.defaultValue.get<bool>();
              break;
            }
          } /// end loop over options
          break;
        }
      } /// end loop over devices
      LOG(info) << "[TrackPropagationModule]  trackTuner.autoDetectDcaCalib it's equal to " << autoDetectDcaCalib;
      if (!autoDetectDcaCalib) {
        LOG(info) << "[TrackPropagationModule]  retrieve the graphs already (we are in propagationService::Init() function)";
        trackTunerObj.getDcaGraphs();
      } else {
        LOG(info) << "[TrackPropagationModule]  trackTunerObj.getDcaGraphs() function to be called later, in the process function!";
      }
    }

    trackTunedTracks = registry.template add<TH1>("trackTunedTracks", outputStringParams.c_str(), o2::framework::kTH1D, {{1, 0.5f, 1.5f}});

    // Histograms for track tuner
    o2::framework::AxisSpec axisBinsDCA = {600, -0.15f, 0.15f, "#it{dca}_{xy} (cm)"};
    registry.template add<TH2>("hDCAxyVsPtRec", "hDCAxyVsPtRec", o2::framework::kTH2F, {axisBinsDCA, cGroup.axisPtQA});
    registry.template add<TH2>("hDCAxyVsPtMC", "hDCAxyVsPtMC", o2::framework::kTH2F, {axisBinsDCA, cGroup.axisPtQA});
    registry.template add<TH2>("hDCAzVsPtRec", "hDCAzVsPtRec", o2::framework::kTH2F, {axisBinsDCA, cGroup.axisPtQA});
    registry.template add<TH2>("hDCAzVsPtMC", "hDCAzVsPtMC", o2::framework::kTH2F, {axisBinsDCA, cGroup.axisPtQA});
  }

  template <bool isMc, typename TConfigurableGroup, typename TCCDBLoader, typename TCollisions, typename TTracks, typename TOutputGroup, typename THistoRegistry>
  void fillTrackTables(TConfigurableGroup const& cGroup, TrackTuner& trackTunerObj, TCCDBLoader const& ccdbLoader, TCollisions const& collisions, TTracks const& tracks, TOutputGroup& cursors, THistoRegistry& registry)
  {

    /// retrieve the TrackTuner calibration graphs *if not done yet*
    /// i.e. if autodetect is required
    if (cGroup.useTrackTuner.value && autoDetectDcaCalib && !trackTunerObj.areGraphsConfigured) {

      /// get the run number from the ccdb loader, already initialized
      const int runNumber = ccdbLoader.runNumber;
      trackTunerObj.setRunNumber(runNumber);

      /// setup the "auto-detected" path based on the run number
      trackTunerObj.getPathInputFileAutomaticFromCCDB();
      trackTunedTracks->SetTitle(trackTunerObj.outputString.c_str());

      /// now that the path is ok, retrieve the graphs
      trackTunerObj.getDcaGraphs();
    }

    if (!fillTracks) {
      return; // suppress everything
    }

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

    for (const auto& track : tracks) {
      if (fillTracksCov) {
        if (fillTracksDCA || fillTracksDCACov) {
          mDcaInfoCov.set(999, 999, 999, 999, 999);
        }
        setTrackParCov(track, mTrackParCov);
        if (cGroup.useTrkPid.value) {
          mTrackParCov.setPID(track.pidForTracking());
        }
      } else {
        if (fillTracksDCA) {
          mDcaInfo[0] = 999;
          mDcaInfo[1] = 999;
        }
        setTrackPar(track, mTrackPar);
        if (cGroup.useTrkPid.value) {
          mTrackPar.setPID(track.pidForTracking());
        }
      }
      // auto trackParCov = getTrackParCov(track);
      o2::aod::track::TrackTypeEnum trackType = (o2::aod::track::TrackTypeEnum)track.trackType();
      // std::array<float, 3> trackPxPyPz;
      // std::array<float, 3> trackPxPyPzTuned = {0.0, 0.0, 0.0};
      double q2OverPtNew = -9999.;
      bool isPropagationOK = true;
      // Only propagate tracks which have passed the innermost wall of the TPC (e.g. skipping loopers etc). Others fill unpropagated.
      if (track.x() < cGroup.minPropagationRadius.value) {
        if (fillTracksCov) {
          if constexpr (isMc) { // checking MC and fillCovMat block begins
            // bool hasMcParticle = track.has_mcParticle();
            if (cGroup.useTrackTuner.value) {
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

        if (track.trackType() == o2::aod::track::TrackIU) {
          if (track.has_collision()) {
            auto const& collision = collisions.rawIteratorAt(track.collisionId());
            if (fillTracksCov) {
              mVtx.setPos({collision.posX(), collision.posY(), collision.posZ()});
              mVtx.setCov(collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ());
              isPropagationOK = o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, mTrackParCov, 2.f, matCorr, &mDcaInfoCov);
            } else {
              isPropagationOK = o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, mTrackPar, 2.f, matCorr, &mDcaInfo);
            }
          } else {
            if (fillTracksCov) {
              mVtx.setPos({ccdbLoader.mMeanVtx->getX(), ccdbLoader.mMeanVtx->getY(), ccdbLoader.mMeanVtx->getZ()});
              mVtx.setCov(ccdbLoader.mMeanVtx->getSigmaX() * ccdbLoader.mMeanVtx->getSigmaX(), 0.0f, ccdbLoader.mMeanVtx->getSigmaY() * ccdbLoader.mMeanVtx->getSigmaY(), 0.0f, 0.0f, ccdbLoader.mMeanVtx->getSigmaZ() * ccdbLoader.mMeanVtx->getSigmaZ());
              isPropagationOK = o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, mTrackParCov, 2.f, matCorr, &mDcaInfoCov);
            } else {
              isPropagationOK = o2::base::Propagator::Instance()->propagateToDCABxByBz({ccdbLoader.mMeanVtx->getX(), ccdbLoader.mMeanVtx->getY(), ccdbLoader.mMeanVtx->getZ()}, mTrackPar, 2.f, matCorr, &mDcaInfo);
            }
          }
          if (isPropagationOK) {
            trackType = o2::aod::track::Track;
          }
        } else {
          if (fillTracksDCA || fillTracksDCACov) {
            if (track.has_collision()) {
              auto const& collision = collisions.rawIteratorAt(track.collisionId());
              mVtx.setPos({collision.posX(), collision.posY(), collision.posZ()});
              if (fillTracksDCACov) {
                mVtx.setCov(collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ());
              }
            } else {
              mVtx.setPos({ccdbLoader.mMeanVtx->getX(), ccdbLoader.mMeanVtx->getY(), ccdbLoader.mMeanVtx->getZ()});
              if (fillTracksDCACov) {
                mVtx.setCov(ccdbLoader.mMeanVtx->getSigmaX() * ccdbLoader.mMeanVtx->getSigmaX(), 0.0f, ccdbLoader.mMeanVtx->getSigmaY() * ccdbLoader.mMeanVtx->getSigmaY(), 0.0f, 0.0f, ccdbLoader.mMeanVtx->getSigmaZ() * ccdbLoader.mMeanVtx->getSigmaZ());
              }
            }

            if (fillTracksDCACov) {
              calculateDCA(mTrackParCov, mVtx, o2::base::Propagator::Instance()->getNominalBz(), &mDcaInfoCov, 999.f);
            } else {
              calculateDCA(mTrackPar, mVtx, o2::base::Propagator::Instance()->getNominalBz(), &mDcaInfo, 999.f);
            }
          }
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
      if (cGroup.useTrackTuner.value && cGroup.fillTrackTunerTable.value) {
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

  //_______________________________________________________________________
  // TTrackPar type: either aod::TrackPar or aod::TrackParCov
  // TVertex type: either math_utils::Point3D<value_t> or o2::dataformats::VertexBase
  // TDCA type: either dim2_t or o2::dataformats::DCA
  template <typename TTrackPar, typename TVertex, typename TDCA>
  bool calculateDCA(TTrackPar trackPar, const TVertex& vtx, double b, TDCA* dca, double maxD)
  {
    // propagate track to DCA to the vertex
    double sn, cs, alp = trackPar.getAlpha();
    o2::math_utils::detail::sincos(alp, sn, cs);
    double x = trackPar.getX(), y = trackPar.getY(), snp = trackPar.getSnp(), csp = gpu::CAMath::Sqrt((1.f - snp) * (1.f + snp));
    double xv = vtx.getX() * cs + vtx.getY() * sn, yv = -vtx.getX() * sn + vtx.getY() * cs, zv = vtx.getZ();
    x -= xv;
    y -= yv;
    // Estimate the impact parameter neglecting the track curvature
    double d = gpu::CAMath::Abs(x * snp - y * csp);
    if (d > maxD) {
      // provide default DCA for failed propag
      if constexpr (requires { trackPar.getSigmaY2(); vtx.getSigmaX2(); }) {
        dca->set(o2::track::DefaultDCA, o2::track::DefaultDCA,
                 o2::track::DefaultDCACov, o2::track::DefaultDCACov, o2::track::DefaultDCACov);
      } else {
        (*dca)[0] = o2::track::DefaultDCA;
        (*dca)[1] = o2::track::DefaultDCA;
      }
      return false;
    }
    double crv = trackPar.getCurvature(b);
    double tgfv = -(crv * x - snp) / (crv * y + csp);
    sn = tgfv / gpu::CAMath::Sqrt(1.f + tgfv * tgfv);
    cs = gpu::CAMath::Sqrt((1.f - sn) * (1.f + sn));
    cs = (gpu::CAMath::Abs(tgfv) > constants::math::Almost0) ? sn / tgfv : constants::math::Almost1;

    x = xv * cs + yv * sn;
    yv = -xv * sn + yv * cs;
    xv = x;

    alp += gpu::CAMath::ASin(sn);
    if (!trackPar.rotate(alp)) {
      // provide default DCA for failed propag
      if constexpr (requires { trackPar.getSigmaY2(); vtx.getSigmaX2(); }) {
        dca->set(o2::track::DefaultDCA, o2::track::DefaultDCA,
                 o2::track::DefaultDCACov, o2::track::DefaultDCACov, o2::track::DefaultDCACov);
      } else {
        (*dca)[0] = o2::track::DefaultDCA;
        (*dca)[1] = o2::track::DefaultDCA;
      }

      return false;
    }
    if constexpr (requires { trackPar.getSigmaY2(); vtx.getSigmaX2(); }) {
      o2::math_utils::detail::sincos(alp, sn, cs);
      auto s2ylocvtx = vtx.getSigmaX2() * sn * sn + vtx.getSigmaY2() * cs * cs - 2. * vtx.getSigmaXY() * cs * sn;
      dca->set(trackPar.getY() - yv, trackPar.getZ() - zv, trackPar.getSigmaY2() + s2ylocvtx, trackPar.getSigmaZY(), trackPar.getSigmaZ2() + vtx.getSigmaZ2());
    } else {
      (*dca)[0] = trackPar.getY() - yv;
      (*dca)[1] = trackPar.getZ() - zv;
    }
    return true;
  }
};

} // namespace common
} // namespace o2

#endif // COMMON_TOOLS_TRACKPROPAGATIONMODULE_H_
