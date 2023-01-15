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
//

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/trackUtilities.h"
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#include "Framework/runDataProcessing.h"

struct propagatorQa {
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  int mRunNumber;
  float d_bz;

  Configurable<float> windowDCA{"windowDCA", 50, "windowDCA"};
  Configurable<int> NbinsX{"NbinsX", 500, "NbinsX"};
  Configurable<int> NbinsDCA{"NbinsDCA", 2000, "NbinsDCA"};
  Configurable<int> NbinsPt{"NbinsPt", 10, "NbinsDCA"};
  Configurable<float> maxXtoConsider{"maxXtoConsider", 10000, "max X to consider"};
  // Operation and minimisation criteria
  Configurable<double> d_bz_input{"d_bz", -999, "bz field, -999 is automatic"};

  // CCDB options
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  o2::track::TrackPar lTrackParametrization;

  void init(InitContext& context)
  {
    mRunNumber = 0;
    d_bz = 0;

    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    // output objects
    const AxisSpec axisX{(int)NbinsX, 0.0f, +250.0f, "X value"};
    const AxisSpec axisDCAxy{(int)NbinsDCA, -windowDCA, windowDCA, "DCA_{xy} (cm)"};

    histos.add("hTrackX", "hTrackX", kTH1F, {axisX});
    histos.add("hUpdateRadii", "hUpdateRadii", kTH1F, {axisX});
    histos.add("hUpdateRadiiusedInSVertexer", "hUpdateRadiiusedInSVertexer", kTH1F, {axisX});

    histos.add("hdcaXYall", "hdcaXYall", kTH1F, {axisDCAxy});
    histos.add("hCircleDCA", "hCircleDCA", kTH1F, {axisDCAxy});
    histos.add("hdcaXYusedInSVertexer", "hdcaXYusedInSVertexer", kTH1F, {axisDCAxy});

    histos.add("hTrackXVsDCA", "hTrackXVsDCA", kTH2F, {axisDCAxy, axisX});
    histos.add("hLastUpdateRadiusVsDCA", "hLastUpdateRadiusVsDCA", kTH2F, {axisDCAxy, axisX});
    histos.add("hTrackXVsCircleDCA", "hTrackXVsCircleDCA", kTH2F, {axisDCAxy, axisX});
    histos.add("hLastUpdateRadiusVsCircleDCA", "hLastUpdateRadiusVsCircleDCA", kTH2F, {axisDCAxy, axisX});

    histos.add("hCircleDCAVsDCA", "hCircleDCAVsDCA", kTH2F, {axisDCAxy, axisDCAxy});
    histos.add("hDeltaDCAs", "hDeltaDCAs", kTH1F, {axisDCAxy});
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }
    auto run3grp_timestamp = bc.timestamp();

    o2::parameters::GRPObject* grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(grpPath, run3grp_timestamp);
    o2::parameters::GRPMagField* grpmag = 0x0;
    if (grpo) {
      if (d_bz_input < -990) {
        // Fetch magnetic field from ccdb for current collision
        d_bz = grpo->getNominalL3Field();
        LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
      } else {
        d_bz = d_bz_input;
      }
    } else {
      grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, run3grp_timestamp);
      if (!grpmag) {
        LOG(fatal) << "Got nullptr from CCDB for path " << grpmagPath << " of object GRPMagField and " << grpPath << " of object GRPObject for timestamp " << run3grp_timestamp;
      }
      if (d_bz_input < -990) {
        // Fetch magnetic field from ccdb for current collision
        d_bz = std::lround(5.f * grpmag->getL3Current() / 30000.f);
        LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
      } else {
        d_bz = d_bz_input;
      }
    }
    mRunNumber = bc.runNumber();
  }

  void process(aod::Collision const& collision, aod::V0s const& V0s, aod::Cascades const& cascades, soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksDCA> const& tracks, aod::BCsWithTimestamps const&)
  {
    /* check the previous run number */
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);

    for (auto& track : tracks) {
      if (track.trackType() != aod::track::TrackIU && track.x() > maxXtoConsider)
        continue;

      std::array<float, 3> pos;
      lTrackParametrization = getTrackPar(track);
      lTrackParametrization.getXYZGlo(pos);
      float lRadiusOfLastUpdate = TMath::Sqrt(pos[0] * pos[0] + pos[1] * pos[1]);
      histos.fill(HIST("hUpdateRadii"), lRadiusOfLastUpdate);

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

      histos.fill(HIST("hTrackX"), lTrackParametrization.getX());
      histos.fill(HIST("hdcaXYall"), lDCA);
      histos.fill(HIST("hCircleDCA"), lCircleDCA);
      histos.fill(HIST("hLastUpdateRadiusVsDCA"), lDCA, lRadiusOfLastUpdate);
      histos.fill(HIST("hTrackXVsDCA"), lDCA, lTrackParametrization.getX());
      histos.fill(HIST("hLastUpdateRadiusVsCircleDCA"), lCircleDCA, lRadiusOfLastUpdate);
      histos.fill(HIST("hTrackXVsCircleDCA"), lCircleDCA, lTrackParametrization.getX());
      histos.fill(HIST("hCircleDCAVsDCA"), lDCA, lCircleDCA);
      histos.fill(HIST("hDeltaDCAs"), lCircleDCA - lDCA);

      // determine if track was used in svertexer
      bool usedInSVertexer = false;
      bool lUsedByV0 = false, lUsedByCascade = false;
      for (auto& V0 : V0s) {
        if (V0.posTrackId() == track.globalIndex()) {
          lUsedByV0 = true;
          break;
        }
        if (V0.negTrackId() == track.globalIndex()) {
          lUsedByV0 = true;
          break;
        }
      }
      for (auto& cascade : cascades) {
        if (cascade.bachelorId() == track.globalIndex()) {
          lUsedByCascade = true;
          break;
        }
      }
      if (lUsedByV0 || lUsedByCascade)
        usedInSVertexer = true;

      if (usedInSVertexer)
        histos.fill(HIST("hUpdateRadiiusedInSVertexer"), lRadiusOfLastUpdate);
      if (usedInSVertexer)
        histos.fill(HIST("hdcaXYusedInSVertexer"), lDCA);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<propagatorQa>(cfgc)};
}
