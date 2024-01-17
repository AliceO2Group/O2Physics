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

// Task to produce a table joinable to the jet tables for hf jet tagging
//
/// copy from Common/TableProducer/trackPropagation.cxx on 23.Nov.23
/// \author Hanseo Park <hanseo.park@cern.ch>

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"
#include "DetectorsBase/Propagator.h"
#include "ReconstructionDataFormats/DCA.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsCalibration/MeanVertexObject.h"
#include "CCDB/BasicCCDBManager.h"
#include "CCDB/CcdbApi.h"
#include "Common/Core/trackUtilities.h"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetTagging.h"
#include "PWGJE/Core/JetTaggingUtilities.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct JetTaggerHFTrackExtension {

  Produces<aod::JTrackTagDcas> jTracksTagDcaTable;
  Produces<aod::JTrackTagDcaCovs> jTracksTagDcaCovTable;

  // CCDM options
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;

  const o2::dataformats::MeanVertexObject* mVtx = nullptr;
  o2::parameters::GRPMagField* grpmag = nullptr;
  o2::base::MatLayerCylSet* lut = nullptr;

  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> mVtxPath{"mVtxPath", "GLO/Calib/MeanVertex", "Path of the mean vertex file"};
  Configurable<float> minPropagationRadius{"minPropagationDistance", o2::constants::geom::XTPCInnerRef + 0.1, "Only tracks which are at a smaller radius will be propagated, defaults to TPC inner wall"};

  // Operation and minimisation criteria
  Configurable<bool> doFillHistogram{"doFillHistogram", false, "fill histogram"};
  ConfigurableAxis binDcaX{"binDcaX", {1001, -1.5f, 1.5f}, ""};
  ConfigurableAxis binDcaY{"binDcaY", {1001, -1.5f, 1.5f}, ""};
  ConfigurableAxis binDcaXY{"binDcaXY", {1001, -1.5f, 1.5f}, ""};
  ConfigurableAxis binDcaZ{"binDcaZ", {1001, -1.5f, 1.5f}, ""};
  ConfigurableAxis binDcaXYZ{"binDcaXYZ", {1001, -1.5f, 1.5f}, ""};
  ConfigurableAxis binDcaXYSigma{"binDcaXYSigma", {1001, -0.5f, 1.5f}, ""};
  ConfigurableAxis binDcaXYZSigma{"binDcaXYZSigma", {1001, -0.5f, 1.5f}, ""};

  // Axis
  AxisSpec DcaXAxis = {binDcaX, "DCA_{X} [#mum]"};
  AxisSpec DcaYAxis = {binDcaY, "DCA_{Y} [#mum]"};
  AxisSpec DcaXYAxis = {binDcaXY, "DCA_{XY} [#mum]"};
  AxisSpec DcaZAxis = {binDcaZ, "DCA_{Z} [#mum]"};
  AxisSpec DcaXYZAxis = {binDcaXYZ, "DCA_{XYZ} [#mum]"};
  AxisSpec DcaXYSigmaAxis = {binDcaXY, "#sigma_{DCA_{XY}} [#mum]"};
  AxisSpec DcaXYZSigmaAxis = {binDcaXY, "#sigma_{DCA_{XYZ}} [#mum]"};

  int runNumber = -1;

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext& initContext)
  {
    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(lutPath));

    if (doFillHistogram) {
      registry.add("hDCAx", "DCA to PV; DCA_{x} (#mum)", {HistType::kTH1F, {DcaXAxis}});
      registry.add("hDCAy", "DCA to PV; DCA_{y} (#mum)", {HistType::kTH1F, {DcaYAxis}});
      registry.add("hDCAz", "DCA to PV; DCA_{z} (#mum)", {HistType::kTH1F, {DcaZAxis}});
      registry.add("hDCAxy", "DCA to PV; DCA_{xy} (#mum)", {HistType::kTH1F, {DcaXYAxis}});
      registry.add("hDCAxyz", "DCA to PV; DCA_{xyz} (#mum)", {HistType::kTH1F, {DcaXYZAxis}});
      registry.add("hDCAxyVsDCAxyFromExt", "DCA to PV; DCA_{xy} (#mum); DCA_{xy}^{ext} (#mum)", {HistType::kTH2F, {{DcaXYAxis}, {DcaXYAxis}}});
      registry.add("hDCAxy_uncertainty", "uncertainty of DCA to PV; #sigma_{DCA_{xy}} (#mum)", {HistType::kTH1F, {DcaXYSigmaAxis}});
      registry.add("hDCAxyz_uncertainty", "uncertainty of DCA to PV; #sigma_{DCA_{xyz}} (#mum)", {HistType::kTH1F, {DcaXYZSigmaAxis}});
    }
  }

  template <typename T>
  void initCCDB(T const& bc)
  {
    if (runNumber == bc.runNumber()) {
      return;
    }
    grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, bc.timestamp());
    LOG(info) << "Setting magnetic field to current " << grpmag->getL3Current() << " A for run " << bc.runNumber() << " from its GRPMagField CCDB object";
    o2::base::Propagator::initFieldFromGRP(grpmag);
    o2::base::Propagator::Instance()->setMatLUT(lut);
    mVtx = ccdb->getForTimeStamp<o2::dataformats::MeanVertexObject>(mVtxPath, bc.timestamp());
    runNumber = bc.runNumber();
  }

  void calculateDcaXYZ(float& dcaXYZ, float& sigmaDcaXYZ2, float dcaXY, float dcaZ, float cYY, float cZY, float cZZ, float sigmaDcaXY2, float sigmaDcaZ2)
  {
    dcaXYZ = std::sqrt(dcaXY * dcaXY + dcaZ * dcaZ);
    float dFdxy = 2 * dcaXY / dcaXYZ;
    float dFdz = 2 * dcaZ / dcaXYZ;
    sigmaDcaXYZ2 = std::abs(cYY * dFdxy * dFdxy + cZZ * dFdz * dFdz + 2 * cZY * dFdxy * dFdz);
  }

  void processDummy(aod::Collision const& collision)
  {
  }
  PROCESS_SWITCH(JetTaggerHFTrackExtension, processDummy, "Dummy process", true);

  void processTracks(soa::Join<aod::Tracks, aod::TracksCov, aod::TrackSelection, aod::TracksDCA, aod::TracksDCACov, aod::McTrackLabels>::iterator const& track, aod::BCsWithTimestamps&, aod::Collisions&)
  {
    if (track.has_collision()) {
      auto collisionbc = track.collision_as<aod::Collisions>();
      auto bc = collisionbc.bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      gpu::gpustd::array<float, 2> dcaInfo;
      dcaInfo[0] = 999;
      dcaInfo[1] = 999;
      auto trackPar = getTrackPar(track);
      if (track.trackType() == aod::track::TrackIU && track.x() < minPropagationRadius) {
        if (track.has_collision()) {
          auto const& collision = track.collision();
          o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackPar, 2.f, matCorr, &dcaInfo);
        } else {
          o2::base::Propagator::Instance()->propagateToDCABxByBz({mVtx->getX(), mVtx->getY(), mVtx->getZ()}, trackPar, 2.f, matCorr, &dcaInfo);
        }
      }
      auto xyz = trackPar.getXYZGlo();
      float dcaX = 999.;
      float dcaY = 999.;
      if (track.has_collision()) {
        auto const& collision = track.collision();
        dcaX = xyz.X() - collision.posX();
        dcaY = xyz.Y() - collision.posY();
      } else {
        dcaX = xyz.X() - mVtx->getX();
        dcaY = xyz.Y() - mVtx->getY();
      }

      float dcaXY = track.dcaXY();
      float dcaXYfromExt = std::sqrt(dcaX * dcaX + dcaY * dcaY);
      float absdcaXY = TMath::Abs(dcaXY);
      if (absdcaXY < dcaXYfromExt || absdcaXY > dcaXYfromExt) {
        LOGF(info, Form("DCA value is not same. abs(dcaXY): %f, dcaXYfromExt: %f", absdcaXY, dcaXYfromExt));
      }
      float dcaZ = track.dcaZ();
      float sigmaDcaXY2 = track.sigmaDcaXY2();
      float sigmaDcaZ2 = track.sigmaDcaZ2();
      float dcaXYZ, sigmaDcaXYZ2;
      calculateDcaXYZ(dcaXYZ, sigmaDcaXYZ2, dcaXY, dcaZ, track.cYY(), track.cZY(), track.cZZ(), sigmaDcaXY2, sigmaDcaZ2);

      jTracksTagDcaTable(dcaX, dcaY, dcaXY, dcaZ, dcaXYZ);
      jTracksTagDcaCovTable(sigmaDcaXY2, sigmaDcaZ2, sigmaDcaXYZ2);

      if (doFillHistogram) {
        registry.fill(HIST("hDCAx"), dcaX);
        registry.fill(HIST("hDCAy"), dcaY);
        registry.fill(HIST("hDCAz"), dcaY);
        registry.fill(HIST("hDCAxy"), dcaXY);
        registry.fill(HIST("hDCAxyz"), dcaXYZ);
        registry.fill(HIST("hDCAxyVsDCAxyFromExt"), absdcaXY, dcaXYfromExt);
        registry.fill(HIST("hDCAxy_uncertainty"), dcaXY / std::sqrt(sigmaDcaXY2));
        registry.fill(HIST("hDCAxyz_uncertainty"), dcaXYZ / std::sqrt(sigmaDcaXYZ2));
      }
    } else {
      // TODO: set size of JTracks that some track doesn't have collision
      jTracksTagDcaTable(0, 0, 0, 0, 0);
      jTracksTagDcaCovTable(0, 0, 0);
    }
  }
  PROCESS_SWITCH(JetTaggerHFTrackExtension, processTracks, "Fill track's information of extension for tagging jet", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<JetTaggerHFTrackExtension>(cfgc, TaskName{"jet-taggerhf-extension"})};
}
