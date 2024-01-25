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

// Task to produce a table extension JTracks for JTracksTagDca and JTracksTagDcaCov values
//
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

  Produces<aod::JTracksTagDca> jtracksTagDcaTable;
  Produces<aod::JTracksTagDcaCov> jtracksTagDcaCovTable;

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
  AxisSpec DcaXAxis = {binDcaX, "DCA_{X} [cm]"};
  AxisSpec DcaYAxis = {binDcaY, "DCA_{Y} [cm]"};
  AxisSpec DcaXYAxis = {binDcaXY, "DCA_{XY} [cm]"};
  AxisSpec DcaXYExtAxis = {binDcaXY, "DCA_{XY}^{Ext} [cm]"};
  AxisSpec DcaZAxis = {binDcaZ, "DCA_{Z} [cm]"};
  AxisSpec DcaXYZAxis = {binDcaXYZ, "DCA_{XYZ} [cm]"};
  AxisSpec DcaXYSigmaAxis = {binDcaXY, "#sigma_{DCA_{XY}} [cm]"};
  AxisSpec DcaXYZSigmaAxis = {binDcaXY, "#sigma_{DCA_{XYZ}} [cm]"};

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext& initContext)
  {
    if (doFillHistogram) {
      registry.add("hDCAx", "DCA to PV; DCA_{x} (cm)", {HistType::kTH1F, {DcaXAxis}});
      registry.add("hDCAy", "DCA to PV; DCA_{y} (cm)", {HistType::kTH1F, {DcaYAxis}});
      registry.add("hDCAz", "DCA to PV; DCA_{z} (cm)", {HistType::kTH1F, {DcaZAxis}});
      registry.add("hDCAxy", "DCA to PV; DCA_{xy} (cm)", {HistType::kTH1F, {DcaXYAxis}});
      registry.add("hDCAxyz", "DCA to PV; DCA_{xyz} (cm)", {HistType::kTH1F, {DcaXYZAxis}});
      registry.add("hDCAxyVsDCAxyFromExt", "DCA to PV; DCA_{xy} (cm); DCA_{xy}^{ext} (cm)", {HistType::kTH2F, {{DcaXYAxis}, {DcaXYExtAxis}}});
      registry.add("hDCAxy_uncertainty", "uncertainty of DCA to PV; #sigma_{DCA_{xy}} (cm)", {HistType::kTH1F, {DcaXYSigmaAxis}});
      registry.add("hDCAxyz_uncertainty", "uncertainty of DCA to PV; #sigma_{DCA_{xyz}} (cm)", {HistType::kTH1F, {DcaXYZSigmaAxis}});
    }
  }

  void processDummy(aod::Collision const& collision)
  {
  }
  PROCESS_SWITCH(JetTaggerHFTrackExtension, processDummy, "Dummy process", true);

  void processTracksDca(soa::Join<aod::Tracks, aod::TracksCov, aod::TrackSelection, aod::TracksDCA, aod::TracksDCACov>::iterator const& track, aod::Collisions&)
  {
    if (track.has_collision()) {
      auto trackPar = getTrackPar(track);
      auto xyz = trackPar.getXYZGlo();
      float dcaX = 999.;
      float dcaY = 999.;
      auto const& collision = track.collision();
      dcaX = xyz.X() - collision.posX();
      dcaY = xyz.Y() - collision.posY();

      float dcaXY = track.dcaXY();
      float dcaXYfromExt = std::sqrt(dcaX * dcaX + dcaY * dcaY);
      float absdcaXY = TMath::Abs(dcaXY);

      float dcaZ = track.dcaZ();
      float sigmaDcaXY2 = track.sigmaDcaXY2();
      float sigmaDcaZ2 = track.sigmaDcaZ2();
      float dcaXYZ, sigmaDcaXYZ2;
      jettaggingutilities::calculateDcaXYZ(dcaXYZ, sigmaDcaXYZ2, dcaXY, dcaZ, track.cYY(), track.cZY(), track.cZZ(), sigmaDcaXY2, sigmaDcaZ2);

      jtracksTagDcaTable(dcaX, dcaY, dcaXY, dcaZ, dcaXYZ);
      jtracksTagDcaCovTable(sigmaDcaXY2, sigmaDcaZ2, sigmaDcaXYZ2);

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
      jtracksTagDcaTable(0, 0, 0, 0, 0);
      jtracksTagDcaCovTable(0, 0, 0);
    }
  }
  PROCESS_SWITCH(JetTaggerHFTrackExtension, processTracksDca, "Fill track's information of extension for tagging jet", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<JetTaggerHFTrackExtension>(cfgc, TaskName{"jet-taggerhf-extension"})};
}
