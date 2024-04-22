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

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "CCDB/BasicCCDBManager.h"
#include "Framework/HistogramRegistry.h"
#include "TPCCalibration/TPCMShapeCorrection.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "DataFormatsParameters/GRPECSObject.h"

#include "TTree.h"

using namespace o2;
using namespace o2::framework;
using BCsRun3 = soa::Join<aod::BCs, aod::Timestamps>;
using BarrelTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA>;
const AxisSpec axisQoverPt{100, -5., 5., "q/p_{T}, 1/GeV"};
const AxisSpec axisDcaR{1000, -5., 5., "DCA_{r}, cm"};
const AxisSpec axisDcaZ{1000, -5., 5., "DCA_{z}, cm"};
const AxisSpec axisSparseQoverPt{20, -5., 5., "q/p_{T}, 1/GeV"};
const AxisSpec axisSparseDcaR{100, -5., 5., "DCA_{r}, cm"};
const AxisSpec axisSparseDcaZ{100, -5., 5., "DCA_{z}, cm"};

struct MshapeQaTask {
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  o2::tpc::TPCMShapeCorrection mshape; // object for simple access
  int lastRunNumber = -1;
  double maxSec = 1;
  double minSec = 0;
  void init(InitContext&)
  {
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    histos.add("hQoverPt", "", kTH1F, {axisQoverPt});
    histos.add("hDcaR", "", kTH1F, {axisDcaR});
    histos.add("hDcaZ", "", kTH1F, {axisDcaZ});
    histos.add("hQoverPtDcaR", "", kTH2F, {axisSparseQoverPt, axisSparseDcaR});
    histos.add("hQoverPtDcaZ", "", kTH2F, {axisSparseQoverPt, axisSparseDcaZ});
  }

  void process(aod::Collision const& col, BCsRun3 const& bcs, BarrelTracks const& tracks)
  {
    int runNumber = bcs.iteratorAt(0).runNumber();
    if (runNumber != lastRunNumber) {
      lastRunNumber = runNumber;
      std::map<std::string, std::string> metadata;
      metadata["runNumber"] = Form("%d", runNumber);
      auto grpecs = ccdb->getSpecific<o2::parameters::GRPECSObject>("GLO/Config/GRPECS", bcs.iteratorAt(0).timestamp(), metadata);
      minSec = floor(grpecs->getTimeStart() / 1000.);
      maxSec = ceil(grpecs->getTimeEnd() / 1000.);
      const AxisSpec axisSeconds{static_cast<int>(maxSec - minSec) * 10, 0, maxSec - minSec, "seconds"};
      histos.add("hSecondsAsideQoverPtSumDcaR", "", kTH2F, {axisSeconds, axisSparseQoverPt});
      histos.add("hSecondsAsideQoverPtSumDcaZ", "", kTH2F, {axisSeconds, axisSparseQoverPt});
      histos.add("hSecondsCsideQoverPtSumDcaR", "", kTH2F, {axisSeconds, axisSparseQoverPt});
      histos.add("hSecondsCsideQoverPtSumDcaZ", "", kTH2F, {axisSeconds, axisSparseQoverPt});
      histos.add("hSecondsQoverPtSumDcaR", "", kTH2F, {axisSeconds, axisSparseQoverPt});
      histos.add("hSecondsQoverPtSumDcaZ", "", kTH2F, {axisSeconds, axisSparseQoverPt});

      histos.add("hSecondsAsideSumDcaR", "", kTH1F, {axisSeconds});
      histos.add("hSecondsAsideSumDcaZ", "", kTH1F, {axisSeconds});
      histos.add("hSecondsCsideSumDcaR", "", kTH1F, {axisSeconds});
      histos.add("hSecondsCsideSumDcaZ", "", kTH1F, {axisSeconds});
      histos.add("hSecondsSumDcaR", "", kTH1F, {axisSeconds});
      histos.add("hSecondsSumDcaZ", "", kTH1F, {axisSeconds});
      histos.add("hSecondsTracks", "", kTH1F, {axisSeconds});
      histos.add("hSecondsTracksMshape", "", kTH1F, {axisSeconds});
      histos.add("hSecondsAsideITSTPCcontrib", "", kTH1F, {axisSeconds});
      histos.add("hSecondsCsideITSTPCcontrib", "", kTH1F, {axisSeconds});
      histos.add("hSecondsCollisions", "", kTH1F, {axisSeconds});

      const AxisSpec axisPhi{64, 0, TMath::TwoPi(), "#varphi"};
      histos.add("hSecondsITSlayer0vsPhi", "", kTH2F, {axisSeconds, axisPhi});
      histos.add("hSecondsITSlayer1vsPhi", "", kTH2F, {axisSeconds, axisPhi});
      histos.add("hSecondsITSlayer2vsPhi", "", kTH2F, {axisSeconds, axisPhi});
      histos.add("hSecondsITSlayer3vsPhi", "", kTH2F, {axisSeconds, axisPhi});
      histos.add("hSecondsITSlayer4vsPhi", "", kTH2F, {axisSeconds, axisPhi});
      histos.add("hSecondsITSlayer5vsPhi", "", kTH2F, {axisSeconds, axisPhi});
      histos.add("hSecondsITSlayer6vsPhi", "", kTH2F, {axisSeconds, axisPhi});
    }

    int64_t ts = col.bc_as<BCsRun3>().timestamp();
    auto mShapeTree = ccdb->getForTimeStamp<TTree>("TPC/Calib/MShapePotential", ts);
    mshape.setFromTree(*mShapeTree);
    bool isMshape = !mshape.getBoundaryPotential(ts).mPotential.empty();

    double secFromSOR = ts / 1000. - minSec;

    int nAsideITSTPCContrib = 0;
    int nCsideITSTPCContrib = 0;
    for (const auto& track : tracks) {
      if (!track.hasTPC() || !track.hasITS()) {
        continue;
      }
      float qpt = track.signed1Pt();
      float dcaR = track.dcaXY();
      float dcaZ = track.dcaZ();
      LOGP(debug, "dcaR = {} dcaZ = {}", dcaR, dcaZ);
      histos.fill(HIST("hQoverPt"), qpt);
      histos.fill(HIST("hDcaR"), dcaR);
      histos.fill(HIST("hDcaZ"), dcaZ);
      histos.fill(HIST("hQoverPtDcaR"), qpt, dcaR);
      histos.fill(HIST("hQoverPtDcaZ"), qpt, dcaZ);
      histos.fill(HIST("hSecondsSumDcaR"), secFromSOR, dcaR);
      histos.fill(HIST("hSecondsSumDcaZ"), secFromSOR, dcaZ);
      histos.fill(HIST("hSecondsQoverPtSumDcaR"), secFromSOR, qpt, dcaR);
      histos.fill(HIST("hSecondsQoverPtSumDcaZ"), secFromSOR, qpt, dcaZ);
      histos.fill(HIST("hSecondsTracks"), secFromSOR);
      if (track.tgl() > 0.) {
        histos.fill(HIST("hSecondsAsideQoverPtSumDcaR"), secFromSOR, qpt, dcaR);
        histos.fill(HIST("hSecondsAsideQoverPtSumDcaZ"), secFromSOR, qpt, dcaZ);
        histos.fill(HIST("hSecondsAsideSumDcaR"), secFromSOR, dcaR);
        histos.fill(HIST("hSecondsAsideSumDcaZ"), secFromSOR, dcaZ);

      } else {
        histos.fill(HIST("hSecondsCsideQoverPtSumDcaR"), secFromSOR, qpt, dcaR);
        histos.fill(HIST("hSecondsCsideQoverPtSumDcaZ"), secFromSOR, qpt, dcaZ);
        histos.fill(HIST("hSecondsCsideSumDcaR"), secFromSOR, dcaR);
        histos.fill(HIST("hSecondsCsideSumDcaZ"), secFromSOR, dcaZ);
      }
      if (isMshape) {
        histos.fill(HIST("hSecondsTracksMshape"), secFromSOR);
      }
      if (track.isPVContributor()) {
        if (track.tgl() > 0.) {
          nAsideITSTPCContrib++;
        } else {
          nCsideITSTPCContrib++;
        }

        // study ITS cluster pattern vs sec
        if (track.itsClusterMap() & (1 << 0))
          histos.fill(HIST("hSecondsITSlayer0vsPhi"), secFromSOR, track.phi());
        if (track.itsClusterMap() & (1 << 1))
          histos.fill(HIST("hSecondsITSlayer1vsPhi"), secFromSOR, track.phi());
        if (track.itsClusterMap() & (1 << 2))
          histos.fill(HIST("hSecondsITSlayer2vsPhi"), secFromSOR, track.phi());
        if (track.itsClusterMap() & (1 << 3))
          histos.fill(HIST("hSecondsITSlayer3vsPhi"), secFromSOR, track.phi());
        if (track.itsClusterMap() & (1 << 4))
          histos.fill(HIST("hSecondsITSlayer4vsPhi"), secFromSOR, track.phi());
        if (track.itsClusterMap() & (1 << 5))
          histos.fill(HIST("hSecondsITSlayer5vsPhi"), secFromSOR, track.phi());
        if (track.itsClusterMap() & (1 << 6))
          histos.fill(HIST("hSecondsITSlayer6vsPhi"), secFromSOR, track.phi());
      }
    }
    histos.fill(HIST("hSecondsCollisions"), secFromSOR);
    histos.fill(HIST("hSecondsAsideITSTPCcontrib"), secFromSOR, nAsideITSTPCContrib);
    histos.fill(HIST("hSecondsCsideITSTPCcontrib"), secFromSOR, nCsideITSTPCContrib);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<MshapeQaTask>(cfgc)};
}
