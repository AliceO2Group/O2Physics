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
    for (auto& track : tracks) {
      float itsChi2NCl = (track.hasITS() && track.itsChi2NCl() > 0.f) ? track.itsChi2NCl() : -299.f;
      float tpcChi2NCl = (track.hasTPC() && track.tpcChi2NCl() > 0.f) ? track.tpcChi2NCl() : -299.f;
      float beta = track.hasTOF() ? track.beta() : -29.f;
      float tofNSigmaEl = track.hasTOF() ? track.tofNSigmaEl() : -299.f;
      float tofNSigmaPi = track.hasTOF() ? track.tofNSigmaPi() : -299.f;
      float tofNSigmaKa = track.hasTOF() ? track.tofNSigmaKa() : -299.f;
      float tofNSigmaPr = track.hasTOF() ? track.tofNSigmaPr() : -299.f;
      float tofChi2 = track.hasTOF() ? track.tofChi2() : -299.f;

      float tpcSignal = track.hasTPC() ? track.tpcSignal() : -299.f;
      float tpcNSigmaEl = track.hasTPC() ? track.tpcNSigmaEl() : -299.f;
      float tpcNSigmaPi = track.hasTPC() ? track.tpcNSigmaPi() : -299.f;
      float tpcNSigmaKa = track.hasTPC() ? track.tpcNSigmaKa() : -299.f;
      float tpcNSigmaPr = track.hasTPC() ? track.tpcNSigmaPr() : -299.f;

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

                static_cast<int16_t>(tpcChi2NCl * 1e+2),
                track.tpcInnerParam(),
                static_cast<uint16_t>(tpcSignal * 1e+2),
                static_cast<int16_t>(tpcNSigmaEl * 1e+2),
                static_cast<int16_t>(tpcNSigmaPi * 1e+2),
                static_cast<int16_t>(tpcNSigmaKa * 1e+2),
                static_cast<int16_t>(tpcNSigmaPr * 1e+2),
                static_cast<int16_t>(beta * 1e+3),
                static_cast<int16_t>(tofNSigmaEl * 1e+2),
                static_cast<int16_t>(tofNSigmaPi * 1e+2),
                static_cast<int16_t>(tofNSigmaKa * 1e+2),
                static_cast<int16_t>(tofNSigmaPr * 1e+2),
                track.itsClusterSizes(),
                static_cast<int16_t>(itsChi2NCl * 1e+2),
                static_cast<int16_t>(tofChi2 * 1e+2),
                track.detectorMap(),
                track.tgl(),
                track.isAssociatedToMPC(),
                false,
                0.f,
                static_cast<uint16_t>(0));
    } // end of track loop
  }
  PROCESS_SWITCH(electronConverter4, process002to004, "convert from 002 into 004", false);

  using MyElectrons003 = soa::Join<aod::EMPrimaryElectrons_003, aod::EMPrimaryElectronsCov_000>;
  void process003to004(MyElectrons003 const& tracks)
  {
    for (auto& track : tracks) {
      float itsChi2NCl = track.itsChi2NCl() > 0.f ? track.itsChi2NCl() : -299.f;
      float tpcChi2NCl = track.tpcChi2NCl() > 0.f ? track.tpcChi2NCl() : -299.f;
      float beta = track.hasTOF() ? track.beta() : -29.f;
      float tofNSigmaEl = track.hasTOF() ? track.tofNSigmaEl() : -299.f;
      float tofNSigmaPi = track.hasTOF() ? track.tofNSigmaPi() : -299.f;
      float tofNSigmaKa = track.hasTOF() ? track.tofNSigmaKa() : -299.f;
      float tofNSigmaPr = track.hasTOF() ? track.tofNSigmaPr() : -299.f;
      float tofChi2 = track.hasTOF() ? track.tofChi2() : -299.f;

      float tpcSignal = track.hasTPC() ? track.tpcSignal() : 0.f;
      float mcTunedTPCSignal = track.hasTPC() ? track.mcTunedTPCSignal() : 0.f;
      float tpcNSigmaEl = track.hasTPC() ? track.tpcNSigmaEl() : -299.f;
      float tpcNSigmaPi = track.hasTPC() ? track.tpcNSigmaPi() : -299.f;
      float tpcNSigmaKa = track.hasTPC() ? track.tpcNSigmaKa() : -299.f;
      float tpcNSigmaPr = track.hasTPC() ? track.tpcNSigmaPr() : -299.f;

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

                static_cast<int16_t>(tpcChi2NCl * 1e+2),
                track.tpcInnerParam(),
                static_cast<uint16_t>(tpcSignal * 1e+2),
                static_cast<int16_t>(tpcNSigmaEl * 1e+2),
                static_cast<int16_t>(tpcNSigmaPi * 1e+2),
                static_cast<int16_t>(tpcNSigmaKa * 1e+2),
                static_cast<int16_t>(tpcNSigmaPr * 1e+2),
                static_cast<int16_t>(beta * 1e+3),
                static_cast<int16_t>(tofNSigmaEl * 1e+2),
                static_cast<int16_t>(tofNSigmaPi * 1e+2),
                static_cast<int16_t>(tofNSigmaKa * 1e+2),
                static_cast<int16_t>(tofNSigmaPr * 1e+2),
                track.itsClusterSizes(),
                static_cast<int16_t>(itsChi2NCl * 1e+2),
                static_cast<int16_t>(tofChi2 * 1e+2),
                track.detectorMap(),
                track.tgl(),
                track.isAssociatedToMPC(),
                false,
                0.f,
                static_cast<uint16_t>(mcTunedTPCSignal));
    } // end of track loop
  }
  PROCESS_SWITCH(electronConverter4, process003to004, "convert from 003 into 004", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<electronConverter4>(cfgc, TaskName{"electron-converter4"})};
}
