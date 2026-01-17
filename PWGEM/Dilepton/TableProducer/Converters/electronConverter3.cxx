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

#include "PWGEM/Dilepton/DataModel/dileptonTables.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

struct electronConverter3 {
  Produces<aod::EMPrimaryElectrons_003> track_003;

  void process(aod::EMPrimaryElectrons_002 const& tracks)
  {
    for (const auto& track : tracks) {
      track_003(track.collisionId(),
                track.trackId(),
                track.sign(),
                track.pt(),
                track.eta(),
                track.phi(),
                track.dcaXY(),
                track.dcaZ(),
                track.tpcNClsFindable(),
                track.tpcNClsFindableMinusFound(),
                track.tpcNClsFindableMinusCrossedRows(),
                track.tpcNClsShared(),
                track.tpcChi2NCl(),
                track.tpcInnerParam(),
                track.tpcSignal(),
                track.tpcNSigmaEl(),
                // track.tpcNSigmaMu(),
                track.tpcNSigmaPi(),
                track.tpcNSigmaKa(),
                track.tpcNSigmaPr(),
                track.beta(),
                track.tofNSigmaEl(),
                // track.tofNSigmaMu(),
                track.tofNSigmaPi(),
                track.tofNSigmaKa(),
                track.tofNSigmaPr(),
                track.itsClusterSizes(),
                // track.itsNSigmaEl(),
                // track.itsNSigmaMu(),
                // track.itsNSigmaPi(),
                // track.itsNSigmaKa(),
                // track.itsNSigmaPr(),
                track.itsChi2NCl(),
                track.tofChi2(),
                track.detectorMap(),
                track.x(),
                track.alpha(),
                track.y(),
                track.z(),
                track.snp(),
                track.tgl(),
                track.isAssociatedToMPC(),
                -1);
    } // end of track loop
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<electronConverter3>(cfgc, TaskName{"electron-converter3"})};
}
