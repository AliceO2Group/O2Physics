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
// Task producing QA histograms to study track (pre-)propagation
// Work in progress! More to follow, use at your own peril
//

#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Tools/trackSelectionRequest.h"

#include <CCDB/BasicCCDBManager.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DataFormatsParameters/GRPObject.h>
#include <DetectorsBase/GeometryManager.h>
#include <DetectorsBase/MatLayerCylSet.h>
#include <DetectorsBase/Propagator.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/DataTypes.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>
#include <MathUtils/Primitive2D.h>
#include <ReconstructionDataFormats/Track.h>

#include <TMath.h>
#include <TMathBase.h>

#include <array>
#include <cmath>
#include <string>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct propagatorQa {
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  int mRunNumber;
  float d_bz;

  o2::base::MatLayerCylSet* lut = nullptr;
  o2::base::Propagator::MatCorrType matCorr;

  // Configurable based on a struct
  Configurable<trackSelectionRequest> trackSels{"trackSels", {}, "track selections"};

  Configurable<float> windowDCA{"windowDCA", 50, "windowDCA"};
  Configurable<int> NbinsX{"NbinsX", 500, "NbinsX"};
  Configurable<int> NbinsDCA{"NbinsDCA", 2000, "NbinsDCA"};
  Configurable<int> NbinsPt{"NbinsPt", 100, "NbinsPt"};
  Configurable<int> NbinsPtCoarse{"NbinsPtCoarse", 100, "NbinsPtCoarse"};
  Configurable<float> maxXtoConsider{"maxXtoConsider", 10000, "max X to consider"};
  Configurable<float> maxPropagStep{"maxPropagStep", 2.0, "max propag step"};
  // Operation and minimisation criteria
  Configurable<double> d_bz_input{"d_bz", -999, "bz field, -999 is automatic"};
  Configurable<int> dQANBinsRadius{"dQANBinsRadius", 100, "binning for radius x itsmap histo"};

  Configurable<int> NbinsTanLambda{"NbinsTanLambda", 100, "binning for tan(lambda)"};
  Configurable<float> TanLambdaLimit{"TanLambdaLimit", 1, "limit for tan(lambda)"};

  Configurable<int> NbinsDeltaPt{"NbinsDeltaPt", 100, "binning for delta-pt"};
  Configurable<float> DeltaPtLimit{"DeltaPtLimit", 1, "limit for delta-pt"};
  Configurable<int> minTPCClustersRequired{"minTPCClustersRequired", -1, "minimum number of TPC clusters"};

  // CCDB options
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};

  Configurable<int> useMatCorrType{"useMatCorrType", 0, "0: none, 1: TGeo, 2: LUT"};

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  o2::track::TrackPar lTrackParametrization;

  void init(InitContext&)
  {
    if (useMatCorrType == 1) {
      LOGF(info, "TGeo correction requested, loading geometry");
      if (!o2::base::GeometryManager::isGeometryLoaded()) {
        ccdb->get<TGeoManager>(geoPath);
      }
    }
    if (useMatCorrType == 2) {
      LOGF(info, "LUT correction requested, loading LUT");
      lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(lutPath));
    }

    mRunNumber = 0;
    d_bz = 0;

    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    // output objects
    const AxisSpec axisX{NbinsX, 0.0f, +250.0f, "X value"};
    const AxisSpec axisDCAxy{NbinsDCA, -windowDCA, windowDCA, "DCA_{xy} (cm)"};
    const AxisSpec axisPt{NbinsPt, 0.0f, 10.0f, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec axisPtCoarse{NbinsPtCoarse, 0.0f, 10.0f, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec axisTanLambda{NbinsTanLambda, -TanLambdaLimit, +TanLambdaLimit, "tan(#lambda)"};
    const AxisSpec axisDeltaPt{NbinsDeltaPt, -DeltaPtLimit, +DeltaPtLimit, "#it{p}_{T} (GeV/#it{c})"};

    // All tracks
    histos.add("hTrackX", "hTrackX", kTH1F, {axisX});
    histos.add("hUpdateRadii", "hUpdateRadii", kTH1F, {axisX});
    histos.add("hdcaXYall", "hdcaXYall", kTH1F, {axisDCAxy});
    histos.add("hCircleDCA", "hCircleDCA", kTH1F, {axisDCAxy});
    histos.add("hTrackXVsDCA", "hTrackXVsDCA", kTH2F, {axisDCAxy, axisX});
    histos.add("hLastUpdateRadiusVsDCA", "hLastUpdateRadiusVsDCA", kTH2F, {axisDCAxy, axisX});
    histos.add("hTrackXVsCircleDCA", "hTrackXVsCircleDCA", kTH2F, {axisDCAxy, axisX});
    histos.add("hLastUpdateRadiusVsCircleDCA", "hLastUpdateRadiusVsCircleDCA", kTH2F, {axisDCAxy, axisX});
    histos.add("hCircleDCAVsDCA", "hCircleDCAVsDCA", kTH2F, {axisDCAxy, axisDCAxy});
    histos.add("hDeltaDCAs", "hDeltaDCAs", kTH1F, {axisDCAxy});
    histos.add("hDeltaDCAsVsPt", "hDeltaDCAsVsPt", kTH2F, {axisPt, axisDCAxy});
    histos.add("hRecalculatedDeltaDCAsVsPt", "hRecalculatedDeltaDCAsVsPt", kTH2F, {axisPt, axisDCAxy});

    // TPC PID checks: difference in tan(lambda) and q/pT between propagated and non propagated
    histos.add("hDeltaTanLambdaVsPt", "hDeltaTanLambdaVsPt", kTH2F, {axisPt, axisTanLambda});
    histos.add("hDeltaPtVsPt", "hDeltaPtVsPt", kTH2F, {axisPt, axisDeltaPt});
    histos.add("hPrimaryDeltaTanLambdaVsPt", "hPrimaryDeltaTanLambdaVsPt", kTH2F, {axisPt, axisTanLambda});
    histos.add("hPrimaryDeltaPtVsPt", "hPrimaryDeltaPtVsPt", kTH2F, {axisPt, axisDeltaPt});

    // Primaries
    histos.add("hPrimaryTrackX", "hPrimaryTrackX", kTH1F, {axisX});
    histos.add("hPrimaryUpdateRadii", "hPrimaryUpdateRadii", kTH1F, {axisX});
    histos.add("hPrimarydcaXYall", "hPrimarydcaXYall", kTH1F, {axisDCAxy});
    histos.add("hPrimaryCircleDCA", "hPrimaryCircleDCA", kTH1F, {axisDCAxy});
    histos.add("hPrimaryTrackXVsDCA", "hPrimaryTrackXVsDCA", kTH2F, {axisDCAxy, axisX});
    histos.add("hPrimaryLastUpdateRadiusVsDCA", "hPrimaryLastUpdateRadiusVsDCA", kTH2F, {axisDCAxy, axisX});
    histos.add("hPrimaryTrackXVsCircleDCA", "hPrimaryTrackXVsCircleDCA", kTH2F, {axisDCAxy, axisX});
    histos.add("hPrimaryLastUpdateRadiusVsCircleDCA", "hPrimaryLastUpdateRadiusVsCircleDCA", kTH2F, {axisDCAxy, axisX});
    histos.add("hPrimaryCircleDCAVsDCA", "hPrimaryCircleDCAVsDCA", kTH2F, {axisDCAxy, axisDCAxy});
    histos.add("hPrimaryDeltaDCAs", "hPrimaryDeltaDCAs", kTH1F, {axisDCAxy});
    histos.add("hPrimaryDeltaDCAsVsPt", "hPrimaryDeltaDCAsVsPt", kTH2F, {axisPt, axisDCAxy});
    histos.add("hPrimaryRecalculatedDeltaDCAsVsPt", "hPrimaryRecalculatedDeltaDCAsVsPt", kTH2F, {axisPt, axisDCAxy});

    // Used in vertexer
    histos.add("hdcaXYusedInSVertexer", "hdcaXYusedInSVertexer", kTH1F, {axisDCAxy});
    histos.add("hUpdateRadiiusedInSVertexer", "hUpdateRadiiusedInSVertexer", kTH1F, {axisX});
    // bit packed ITS cluster map
    const AxisSpec axisITSCluMap{128, -0.5f, +127.5f, "Packed ITS map"};
    const AxisSpec axisRadius{dQANBinsRadius, 0.0f, +50.0f, "Radius (cm)"};

    // Histogram to bookkeep cluster maps
    histos.add("h2dITSCluMap", "h2dITSCluMap", kTH3D, {axisITSCluMap, axisRadius, axisPtCoarse});
    histos.add("h2dITSCluMapPrimaries", "h2dITSCluMapPrimaries", kTH3D, {axisITSCluMap, axisRadius, axisPtCoarse});

    // Material correction
    matCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;
    if (useMatCorrType == 1)
      matCorr = o2::base::Propagator::MatCorrType::USEMatCorrTGeo;
    if (useMatCorrType == 2)
      matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }

    // In case override, don't proceed, please - no CCDB access required
    if (d_bz_input > -990) {
      d_bz = d_bz_input;
      o2::parameters::GRPMagField grpmag;
      if (std::fabs(d_bz) > 1e-5) {
        grpmag.setL3Current(30000.f / (d_bz / 5.0f));
      }
      o2::base::Propagator::initFieldFromGRP(&grpmag);
      mRunNumber = bc.runNumber();
      return;
    }

    auto run3grp_timestamp = bc.timestamp();
    o2::parameters::GRPObject* grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(grpPath, run3grp_timestamp);
    o2::parameters::GRPMagField* grpmag = 0x0;
    if (grpo) {
      o2::base::Propagator::initFieldFromGRP(grpo);
      // Fetch magnetic field from ccdb for current collision
      d_bz = grpo->getNominalL3Field();
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
    } else {
      grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, run3grp_timestamp);
      if (!grpmag) {
        LOG(fatal) << "Got nullptr from CCDB for path " << grpmagPath << " of object GRPMagField and " << grpPath << " of object GRPObject for timestamp " << run3grp_timestamp;
      }
      o2::base::Propagator::initFieldFromGRP(grpmag);
      // Fetch magnetic field from ccdb for current collision
      d_bz = std::lround(5.f * grpmag->getL3Current() / 30000.f);
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
    }
    mRunNumber = bc.runNumber();

    if (useMatCorrType == 2) {
      o2::base::Propagator::Instance()->setMatLUT(lut);
    }
  }

  void processMC(aod::Collision const& collision, aod::V0s const& V0s, aod::Cascades const& cascades, soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksDCA, aod::McTrackLabels> const& tracks, aod::BCsWithTimestamps const&, aod::McParticles const&)
  {
    /* check the previous run number */
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);
    std::array<float, 2> dcaInfo;

    for (const auto& track : tracks) {
      if (track.tpcNClsFound() < minTPCClustersRequired)
        continue;

      if (!track.has_mcParticle())
        continue;
      auto mctrack = track.mcParticle();
      bool lIsPrimary = mctrack.isPhysicalPrimary();

      if (track.trackType() != aod::track::TrackIU && track.x() > maxXtoConsider)
        continue;

      std::array<float, 3> pos;
      lTrackParametrization = getTrackPar(track);
      lTrackParametrization.getXYZGlo(pos);
      float lRadiusOfLastUpdate = TMath::Sqrt(pos[0] * pos[0] + pos[1] * pos[1]);
      float lDCA = track.dcaXY();

      //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
      // Simple snippet for analytical DCA (no e-loss)
      // |<----L---->|
      // |<-R->|<-d->|
      // *-----)     X
      // |           ^ primary vertex
      // ^ circle center

      o2::math_utils::CircleXYf_t lCircle;
      float sna, csa;
      lTrackParametrization.getCircleParams(d_bz, lCircle, sna, csa);
      float lR = lCircle.rC;
      float lL = TMath::Sqrt(
        TMath::Power(lCircle.xC - collision.posX(), 2) +
        TMath::Power(lCircle.yC - collision.posY(), 2));
      float lCircleDCA = lTrackParametrization.getSign() * (lL - lR); // signed dca
      //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*

      dcaInfo[0] = 999;
      dcaInfo[1] = 999;

      //*+-+*
      // Recalculate the propagation
      o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, lTrackParametrization, maxPropagStep, matCorr, &dcaInfo);
      float lRecalculatedDCA = dcaInfo[0];
      //*+-+*
      histos.fill(HIST("hDeltaTanLambdaVsPt"), track.tgl(), track.tgl() - lTrackParametrization.getTgl());
      histos.fill(HIST("hDeltaPtVsPt"), track.pt(), track.pt() - lTrackParametrization.getPt());

      histos.fill(HIST("hUpdateRadii"), lRadiusOfLastUpdate);
      histos.fill(HIST("hTrackX"), lTrackParametrization.getX());
      histos.fill(HIST("hdcaXYall"), lDCA);
      histos.fill(HIST("hCircleDCA"), lCircleDCA);
      histos.fill(HIST("hLastUpdateRadiusVsDCA"), lDCA, lRadiusOfLastUpdate);
      histos.fill(HIST("hTrackXVsDCA"), lDCA, lTrackParametrization.getX());
      histos.fill(HIST("hLastUpdateRadiusVsCircleDCA"), lCircleDCA, lRadiusOfLastUpdate);
      histos.fill(HIST("hTrackXVsCircleDCA"), lCircleDCA, lTrackParametrization.getX());
      histos.fill(HIST("hCircleDCAVsDCA"), lDCA, lCircleDCA);
      histos.fill(HIST("hDeltaDCAs"), lCircleDCA - lDCA);
      histos.fill(HIST("hDeltaDCAsVsPt"), track.pt(), lCircleDCA - lDCA);
      histos.fill(HIST("hRecalculatedDeltaDCAsVsPt"), track.pt(), lRecalculatedDCA - lDCA);

      // ITS cluster map
      float lMCCreation = TMath::Sqrt(mctrack.vx() * mctrack.vx() + mctrack.vy() * mctrack.vy());

      histos.fill(HIST("h2dITSCluMap"), static_cast<float>(track.itsClusterMap()), lMCCreation, track.pt());

      if (lIsPrimary) {
        histos.fill(HIST("hPrimaryDeltaTanLambdaVsPt"), track.tgl(), track.tgl() - lTrackParametrization.getTgl());
        histos.fill(HIST("hPrimaryDeltaPtVsPt"), track.pt(), track.pt() - lTrackParametrization.getPt());

        histos.fill(HIST("hPrimaryUpdateRadii"), lRadiusOfLastUpdate);
        histos.fill(HIST("hPrimaryTrackX"), lTrackParametrization.getX());
        histos.fill(HIST("hPrimarydcaXYall"), lDCA);
        histos.fill(HIST("hPrimaryCircleDCA"), lCircleDCA);
        histos.fill(HIST("hPrimaryLastUpdateRadiusVsDCA"), lDCA, lRadiusOfLastUpdate);
        histos.fill(HIST("hPrimaryTrackXVsDCA"), lDCA, lTrackParametrization.getX());
        histos.fill(HIST("hPrimaryLastUpdateRadiusVsCircleDCA"), lCircleDCA, lRadiusOfLastUpdate);
        histos.fill(HIST("hPrimaryTrackXVsCircleDCA"), lCircleDCA, lTrackParametrization.getX());
        histos.fill(HIST("hPrimaryCircleDCAVsDCA"), lDCA, lCircleDCA);
        histos.fill(HIST("hPrimaryDeltaDCAs"), lCircleDCA - lDCA);
        histos.fill(HIST("hPrimaryDeltaDCAsVsPt"), track.pt(), lCircleDCA - lDCA);
        histos.fill(HIST("hPrimaryRecalculatedDeltaDCAsVsPt"), track.pt(), lRecalculatedDCA - lDCA);
        histos.fill(HIST("h2dITSCluMapPrimaries"), static_cast<float>(track.itsClusterMap()), lMCCreation, track.pt());
      }
      // determine if track was used in svertexer
      bool usedInSVertexer = false;
      bool lUsedByV0 = false, lUsedByCascade = false;
      for (const auto& V0 : V0s) {
        if (V0.posTrackId() == track.globalIndex()) {
          lUsedByV0 = true;
          break;
        }
        if (V0.negTrackId() == track.globalIndex()) {
          lUsedByV0 = true;
          break;
        }
      }
      for (const auto& cascade : cascades) {
        if (cascade.bachelorId() == track.globalIndex()) {
          lUsedByCascade = true;
          break;
        }
      }
      if (lUsedByV0 || lUsedByCascade)
        usedInSVertexer = true;

      if (usedInSVertexer) {
        histos.fill(HIST("hUpdateRadiiusedInSVertexer"), lRadiusOfLastUpdate);
        histos.fill(HIST("hdcaXYusedInSVertexer"), lDCA);
      }
    }
  }
  PROCESS_SWITCH(propagatorQa, processMC, "process MC", true);

  void processData(aod::Collision const& collision, aod::V0s const& V0s, aod::Cascades const& cascades, soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksDCA> const& tracks, aod::BCsWithTimestamps const&)
  {
    /* check the previous run number */
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);
    std::array<float, 2> dcaInfo;

    for (const auto& track : tracks) {
      if (track.tpcNClsFound() < minTPCClustersRequired)
        continue;

      if (track.trackType() != aod::track::TrackIU && track.x() > maxXtoConsider)
        continue;

      std::array<float, 3> pos;
      lTrackParametrization = getTrackPar(track);
      lTrackParametrization.getXYZGlo(pos);
      float lRadiusOfLastUpdate = TMath::Sqrt(pos[0] * pos[0] + pos[1] * pos[1]);
      float lDCA = track.dcaXY();

      //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
      // Simple snippet for analytical DCA (no e-loss)
      // |<----L---->|
      // |<-R->|<-d->|
      // *-----)     X
      // |           ^ primary vertex
      // ^ circle center

      o2::math_utils::CircleXYf_t lCircle;
      float sna, csa;
      lTrackParametrization.getCircleParams(d_bz, lCircle, sna, csa);
      float lR = lCircle.rC;
      float lL = TMath::Sqrt(
        TMath::Power(lCircle.xC - collision.posX(), 2) +
        TMath::Power(lCircle.yC - collision.posY(), 2));
      float lCircleDCA = TMath::Sign(-1, d_bz) * lTrackParametrization.getSign() * (lL - lR); // signed dca
      //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*

      dcaInfo[0] = 999;
      dcaInfo[1] = 999;

      //*+-+*
      // Recalculate the propagation
      o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, lTrackParametrization, maxPropagStep, matCorr, &dcaInfo);
      float lRecalculatedDCA = dcaInfo[0];
      //*+-+*
      histos.fill(HIST("hDeltaTanLambdaVsPt"), track.tgl(), track.tgl() - lTrackParametrization.getTgl());
      histos.fill(HIST("hDeltaPtVsPt"), track.pt(), track.pt() - lTrackParametrization.getPt());

      histos.fill(HIST("hUpdateRadii"), lRadiusOfLastUpdate);
      histos.fill(HIST("hTrackX"), lTrackParametrization.getX());
      histos.fill(HIST("hdcaXYall"), lDCA);
      histos.fill(HIST("hCircleDCA"), lCircleDCA);
      histos.fill(HIST("hLastUpdateRadiusVsDCA"), lDCA, lRadiusOfLastUpdate);
      histos.fill(HIST("hTrackXVsDCA"), lDCA, lTrackParametrization.getX());
      histos.fill(HIST("hLastUpdateRadiusVsCircleDCA"), lCircleDCA, lRadiusOfLastUpdate);
      histos.fill(HIST("hTrackXVsCircleDCA"), lCircleDCA, lTrackParametrization.getX());
      histos.fill(HIST("hCircleDCAVsDCA"), lDCA, lCircleDCA);
      histos.fill(HIST("hDeltaDCAs"), lCircleDCA - lDCA);
      histos.fill(HIST("hDeltaDCAsVsPt"), track.pt(), lCircleDCA - lDCA);
      histos.fill(HIST("hRecalculatedDeltaDCAsVsPt"), track.pt(), lRecalculatedDCA - lDCA);

      // ITS cluster map
      float lMCCreation = 0.1; // dummy value, we don't know

      histos.fill(HIST("h2dITSCluMap"), static_cast<float>(track.itsClusterMap()), lMCCreation, track.pt());

      // A hack: use DCA as equiv to primary
      if (TMath::Abs(lDCA) < 0.05) { // 500 microns
        histos.fill(HIST("hPrimaryDeltaTanLambdaVsPt"), track.tgl(), track.tgl() - lTrackParametrization.getTgl());
        histos.fill(HIST("hPrimaryDeltaPtVsPt"), track.pt(), track.pt() - lTrackParametrization.getPt());
        histos.fill(HIST("hPrimaryUpdateRadii"), lRadiusOfLastUpdate);
        histos.fill(HIST("hPrimaryTrackX"), lTrackParametrization.getX());
        histos.fill(HIST("hPrimarydcaXYall"), lDCA);
        histos.fill(HIST("hPrimaryCircleDCA"), lCircleDCA);
        histos.fill(HIST("hPrimaryLastUpdateRadiusVsDCA"), lDCA, lRadiusOfLastUpdate);
        histos.fill(HIST("hPrimaryTrackXVsDCA"), lDCA, lTrackParametrization.getX());
        histos.fill(HIST("hPrimaryLastUpdateRadiusVsCircleDCA"), lCircleDCA, lRadiusOfLastUpdate);
        histos.fill(HIST("hPrimaryTrackXVsCircleDCA"), lCircleDCA, lTrackParametrization.getX());
        histos.fill(HIST("hPrimaryCircleDCAVsDCA"), lDCA, lCircleDCA);
        histos.fill(HIST("hPrimaryDeltaDCAs"), lCircleDCA - lDCA);
        histos.fill(HIST("hPrimaryDeltaDCAsVsPt"), track.pt(), lCircleDCA - lDCA);
        histos.fill(HIST("hPrimaryRecalculatedDeltaDCAsVsPt"), track.pt(), lRecalculatedDCA - lDCA);
        histos.fill(HIST("h2dITSCluMapPrimaries"), static_cast<float>(track.itsClusterMap()), lMCCreation, track.pt());
      }

      // determine if track was used in svertexer
      bool usedInSVertexer = false;
      bool lUsedByV0 = false, lUsedByCascade = false;
      for (const auto& V0 : V0s) {
        if (V0.posTrackId() == track.globalIndex()) {
          lUsedByV0 = true;
          break;
        }
        if (V0.negTrackId() == track.globalIndex()) {
          lUsedByV0 = true;
          break;
        }
      }
      for (const auto& cascade : cascades) {
        if (cascade.bachelorId() == track.globalIndex()) {
          lUsedByCascade = true;
          break;
        }
      }
      if (lUsedByV0 || lUsedByCascade)
        usedInSVertexer = true;

      if (usedInSVertexer) {
        histos.fill(HIST("hUpdateRadiiusedInSVertexer"), lRadiusOfLastUpdate);
        histos.fill(HIST("hdcaXYusedInSVertexer"), lDCA);
      }
    }
  }
  PROCESS_SWITCH(propagatorQa, processData, "process data", false);

  void processMatLUTTest(aod::Collision const& collision, soa::Join<aod::Tracks, aod::TracksCov, aod::TracksExtra, aod::TracksDCA> const& tracks, soa::Join<aod::TracksIU, aod::TracksCovIU, aod::TracksExtra> const& tracksIU, aod::BCsWithTimestamps const&)
  {
    /* check the previous run number */
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);
    std::array<float, 2> dcaInfo;

    for (const auto& trackIU : tracksIU) {
      if (trackIU.tpcNClsFound() < minTPCClustersRequired)
        continue; // skip if not enough TPC clusters

      if (trackIU.trackType() != aod::track::TrackIU && trackIU.x() > maxXtoConsider)
        continue; // skip if not track IU or if beyong the max X to be considered

      o2::track::TrackParCov trackParCov = getTrackParCov(trackIU);

      dcaInfo[0] = 999;
      dcaInfo[1] = 999;

      // Recalculate the propagation with this instance of the matLUT
      o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParCov, maxPropagStep, matCorr, &dcaInfo);

      auto track = tracks.iteratorAt(trackIU.globalIndex() - tracksIU.offset());

      histos.fill(HIST("hDeltaTanLambdaVsPt"), track.tgl(), track.tgl() - trackParCov.getTgl());
      histos.fill(HIST("hDeltaPtVsPt"), track.pt(), track.pt() - trackParCov.getPt());
      histos.fill(HIST("hDeltaDCAs"), track.dcaXY() - dcaInfo[0]);
      histos.fill(HIST("hDeltaDCAsVsPt"), track.pt(), track.dcaXY() - dcaInfo[0]);
    }
  }
  PROCESS_SWITCH(propagatorQa, processMatLUTTest, "process mat lut test", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<propagatorQa>(cfgc)};
}
