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
// ========================
//
// This code runs loop over ULS ee pars for virtual photon QC.
//    Please write to: daiki.sekihata@cern.ch

#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

struct pcmConverter1 {
  Produces<aod::V0PhotonsKF_001> v0photon_001;
  Produces<aod::V0Legs_001> v0leg_001;

  void process(aod::V0PhotonsKF_000 const& v0s, aod::V0Legs_000 const& v0legs)
  {
    for (auto& v0 : v0s) {
      v0photon_001(
        v0.collisionId(),
        v0.v0Id(),
        v0.posTrackId(),
        v0.negTrackId(),
        v0.vx(),
        v0.vy(),
        v0.vz(),
        v0.px(),
        v0.py(),
        v0.pz(),
        v0.mGamma(),
        v0.dcaXYtopv(),
        v0.dcaZtopv(),
        v0.cospa(),
        v0.cospaXY(),
        v0.cospaRZ(),
        v0.pca(),
        v0.alpha(),
        v0.qtarm(),
        v0.chiSquareNDF());
    } // end of v0 loop

    for (auto& v0leg : v0legs) {
      v0leg_001(
        v0leg.collisionId(),
        v0leg.trackId(),
        v0leg.sign(),
        v0leg.px(),
        v0leg.py(),
        v0leg.pz(),
        v0leg.dcaXY(),
        v0leg.dcaZ(),
        v0leg.tpcNClsFindable(),
        v0leg.tpcNClsFindableMinusFound(),
        v0leg.tpcNClsFindableMinusCrossedRows(),
        v0leg.tpcNClsShared(),
        v0leg.tpcChi2NCl(),
        v0leg.tpcInnerParam(),
        v0leg.tpcSignal(),
        v0leg.tpcNSigmaEl(),
        v0leg.tpcNSigmaPi(),
        v0leg.itsClusterSizes(),
        v0leg.itsChi2NCl(),
        v0leg.detectorMap());
    } // end of v0leg loop
  } // end of process
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<pcmConverter1>(cfgc, TaskName{"pcm-converter1"})};
}
