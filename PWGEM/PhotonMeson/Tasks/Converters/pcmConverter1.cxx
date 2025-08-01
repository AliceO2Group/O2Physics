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
        static_cast<uint16_t>(v0.mGamma() * 1e+5),
        static_cast<int16_t>(v0.dcaXYtopv() * 1e+2),
        static_cast<int16_t>(v0.dcaZtopv() * 1e+2),
        static_cast<uint16_t>(v0.cospa() * 1e+4),
        static_cast<uint16_t>(v0.cospaXY() * 1e+4),
        static_cast<uint16_t>(v0.cospaRZ() * 1e+4),
        static_cast<uint16_t>(v0.pca() * 1e+4),
        static_cast<int16_t>(v0.alpha() * 1e+4),
        static_cast<uint16_t>(v0.qtarm() * 1e+5),
        static_cast<uint16_t>(v0.chiSquareNDF() * 1e+2));
    } // end of v0 loop

    for (auto& v0leg : v0legs) {

      float itsChi2NCl = (v0leg.hasITS() && v0leg.itsChi2NCl() > 0.f) ? v0leg.itsChi2NCl() : -299.f;
      float tpcChi2NCl = (v0leg.hasTPC() && v0leg.tpcChi2NCl() > 0.f) ? v0leg.tpcChi2NCl() : -299.f;
      float tpcSignal = v0leg.hasTPC() ? v0leg.tpcSignal() : 0.f;
      float tpcNSigmaEl = v0leg.hasTPC() ? v0leg.tpcNSigmaEl() : -299.f;
      float tpcNSigmaPi = v0leg.hasTPC() ? v0leg.tpcNSigmaPi() : -299.f;

      v0leg_001(
        v0leg.collisionId(),
        v0leg.trackId(),
        v0leg.sign(),
        v0leg.px(),
        v0leg.py(),
        v0leg.pz(),
        static_cast<int16_t>(v0leg.dcaXY() * 1e+2),
        static_cast<int16_t>(v0leg.dcaZ() * 1e+2),
        v0leg.tpcNClsFindable(),
        v0leg.tpcNClsFindableMinusFound(),
        v0leg.tpcNClsFindableMinusCrossedRows(),
        v0leg.tpcNClsShared(),
        static_cast<int16_t>(tpcChi2NCl * 1e+2),
        v0leg.tpcInnerParam(),
        static_cast<uint16_t>(tpcSignal * 1e+2),
        static_cast<int16_t>(tpcNSigmaEl * 1e+2),
        static_cast<int16_t>(tpcNSigmaPi * 1e+2),
        v0leg.itsClusterSizes(),
        static_cast<int16_t>(itsChi2NCl * 1e+2),
        v0leg.detectorMap(),
        static_cast<uint16_t>(0),
        static_cast<uint16_t>(v0leg.x() * 1e+2),
        static_cast<int16_t>(v0leg.y() * 1e+2),
        static_cast<int16_t>(v0leg.z() * 1e+2),
        v0leg.tgl());
    } // end of v0leg loop
  } // end of process
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<pcmConverter1>(cfgc, TaskName{"pcm-converter1"})};
}
