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

struct electronConverter4 {
  Produces<aod::EMPrimaryElectrons_004> track_004;

  using MyElectrons002 = soa::Join<aod::EMPrimaryElectrons_002, aod::EMPrimaryElectronsCov_000>;
  void process002to004(MyElectrons002 const& tracks)
  {
    for (const auto& track : tracks) {
      track_004(track.collisionId(),
                track.trackId(),
                track.sign(),
                track.pt(),
                track.eta(),
                track.phi(),
                track.dcaXY(),
                track.dcaZ(),
                track.cYY(),
                track.cZY(),
                track.cZZ(),
                track.tpcNClsFindable(),
                track.tpcNClsFindableMinusFound(),
                track.tpcNClsFindableMinusCrossedRows(),
                track.tpcNClsShared(),
                track.tpcChi2NCl(),
                track.tpcInnerParam(),
                track.tpcSignal(),
                track.tpcNSigmaEl(),
                track.tpcNSigmaPi(),
                track.tpcNSigmaKa(),
                track.tpcNSigmaPr(),
                track.beta(),
                track.tofNSigmaEl(),
                // track.tofNSigmaPi(),
                // track.tofNSigmaKa(),
                // track.tofNSigmaPr(),
                track.itsClusterSizes(),
                track.itsChi2NCl(),
                track.tofChi2(),
                track.detectorMap(),
                // track.tgl(),
                track.isAssociatedToMPC(),
                false,
                0.f,
                0.f);
    } // end of track loop
  }
  PROCESS_SWITCH(electronConverter4, process002to004, "convert from 002 into 004", false);

  using MyElectrons003 = soa::Join<aod::EMPrimaryElectrons_003, aod::EMPrimaryElectronsCov_000>;
  void process003to004(MyElectrons003 const& tracks)
  {
    for (const auto& track : tracks) {
      track_004(track.collisionId(),
                track.trackId(),
                track.sign(),
                track.pt(),
                track.eta(),
                track.phi(),
                track.dcaXY(),
                track.dcaZ(),
                track.cYY(),
                track.cZY(),
                track.cZZ(),
                track.tpcNClsFindable(),
                track.tpcNClsFindableMinusFound(),
                track.tpcNClsFindableMinusCrossedRows(),
                track.tpcNClsShared(),
                track.tpcChi2NCl(),
                track.tpcInnerParam(),
                track.tpcSignal(),
                track.tpcNSigmaEl(),
                track.tpcNSigmaPi(),
                track.tpcNSigmaKa(),
                track.tpcNSigmaPr(),
                track.beta(),
                track.tofNSigmaEl(),
                // track.tofNSigmaPi(),
                // track.tofNSigmaKa(),
                // track.tofNSigmaPr(),
                track.itsClusterSizes(),
                track.itsChi2NCl(),
                track.tofChi2(),
                track.detectorMap(),
                // track.tgl(),
                track.isAssociatedToMPC(),
                false,
                0.f,
                track.mcTunedTPCSignal());
    } // end of track loop
  }
  PROCESS_SWITCH(electronConverter4, process003to004, "convert from 003 into 004", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<electronConverter4>(cfgc, TaskName{"electron-converter4"})};
}
