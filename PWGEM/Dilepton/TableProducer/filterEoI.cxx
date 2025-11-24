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
// This code filters events that are interesting for dilepton analyses.
//    Please write to: daiki.sekihata@cern.ch

#include "PWGEM/Dilepton/DataModel/dileptonTables.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

struct filterEoI {
  enum SubSystem {
    kElectron = 0x1,
    kFwdMuon = 0x2,
    kPCM = 0x4,
    kElectronFromDalitz = 0x8,
  };
  Produces<o2::aod::EMEoIs> emeoi;
  Configurable<bool> inheritFromOtherTask{"inheritFromOtherTask", true, "Flag to iherit all common configurables from skimmerPrimaryElectron or skimmerPrimaryMuon"};
  Configurable<int> minNelectron{"minNelectron", -1, "min number of electron candidates per collision"};
  Configurable<int> minNmuon{"minNmuon", -1, "min number of muon candidates per collision"};

  HistogramRegistry fRegistry{"output"};
  void init(o2::framework::InitContext& initContext)
  {
    if (inheritFromOtherTask.value) { // Inheriting from other task
      getTaskOptionValue(initContext, "skimmer-primary-electron", "minNelectron", minNelectron.value, true);
      getTaskOptionValue(initContext, "skimmer-primary-muon", "minNmuon", minNmuon.value, true);
    }

    LOGF(info, "minNelectron = %d", minNelectron.value);
    LOGF(info, "minNmuon = %d", minNmuon.value);

    auto hEventCounter = fRegistry.add<TH1>("hEventCounter", "hEventCounter", kTH1D, {{8, 0.5f, 8.5f}});
    hEventCounter->GetXaxis()->SetBinLabel(1, "all");
    hEventCounter->GetXaxis()->SetBinLabel(2, "event with electron");
    hEventCounter->GetXaxis()->SetBinLabel(3, "event with forward muon");
    hEventCounter->GetXaxis()->SetBinLabel(4, "event with v0");
    hEventCounter->GetXaxis()->SetBinLabel(5, "event with electron or forward muon");
    hEventCounter->GetXaxis()->SetBinLabel(6, "event with electron and forward muon");
    hEventCounter->GetXaxis()->SetBinLabel(7, "event with electron or forward muon or v0");
    hEventCounter->GetXaxis()->SetBinLabel(8, "event with v0 or electrons from dalitz");
  }

  SliceCache cache;
  Preslice<aod::EMPrimaryElectrons> perCollision_el = aod::emprimaryelectron::collisionId;
  Preslice<aod::EMPrimaryMuons> perCollision_mu = aod::emprimarymuon::collisionId;
  Preslice<aod::V0PhotonsKF> perCollision_v0 = aod::v0photonkf::collisionId;
  Preslice<aod::EMPrimaryElectronsFromDalitz> perCollision_elda = aod::emprimaryelectron::collisionId;

  template <uint8_t system, typename TCollisions, typename TElectrons, typename TMuons, typename TV0s, typename TElectronsDA>
  void selectEoI(TCollisions const& collisions, TElectrons const& electrons, TMuons const& muons, TV0s const& v0s, TElectronsDA const& electronsda)
  {
    for (const auto& collision : collisions) {
      bool does_electron_exist = false;
      bool does_fwdmuon_exist = false;
      bool does_pcm_exist = false;
      bool does_electronda_exist = false;
      fRegistry.fill(HIST("hEventCounter"), 1);

      if constexpr (static_cast<bool>(system & kElectron)) {
        auto electrons_coll = electrons.sliceBy(perCollision_el, collision.globalIndex());
        if (electrons_coll.size() >= minNelectron) {
          does_electron_exist = true;
          fRegistry.fill(HIST("hEventCounter"), 2);
        }
      }
      if constexpr (static_cast<bool>(system & kFwdMuon)) {
        auto muons_coll = muons.sliceBy(perCollision_mu, collision.globalIndex());
        if (muons_coll.size() >= minNmuon) {
          does_fwdmuon_exist = true;
          fRegistry.fill(HIST("hEventCounter"), 3);
        }
      }
      if constexpr (static_cast<bool>(system & kPCM)) {
        auto v0s_coll = v0s.sliceBy(perCollision_v0, collision.globalIndex());
        if (v0s_coll.size() >= 1) {
          does_pcm_exist = true;
          fRegistry.fill(HIST("hEventCounter"), 4);
        }
      }
      if constexpr (static_cast<bool>(system & kElectronFromDalitz)) {
        auto electronsda_coll = electronsda.sliceBy(perCollision_elda, collision.globalIndex());
        if (electronsda_coll.size() >= 2) {
          does_electronda_exist = true;
        }
      }

      if (does_electron_exist || does_fwdmuon_exist) {
        fRegistry.fill(HIST("hEventCounter"), 5);
      }
      if (does_electron_exist && does_fwdmuon_exist) {
        fRegistry.fill(HIST("hEventCounter"), 6);
      }
      if (does_electron_exist || does_fwdmuon_exist || does_pcm_exist) {
        fRegistry.fill(HIST("hEventCounter"), 7);
      }
      if (does_pcm_exist || does_electronda_exist) {
        fRegistry.fill(HIST("hEventCounter"), 8);
      }

      emeoi(does_electron_exist || does_fwdmuon_exist || does_pcm_exist || does_electronda_exist);

    } // end of collision loop

  } // end of selectEoI

  void process_Electron(aod::Collisions const& collisions, aod::EMPrimaryElectrons const& electrons)
  {
    const uint8_t sysflag = kElectron;
    selectEoI<sysflag>(collisions, electrons, nullptr, nullptr, nullptr);
  }

  void process_FwdMuon(aod::Collisions const& collisions, aod::EMPrimaryMuons const& muons)
  {
    const uint8_t sysflag = kFwdMuon;
    selectEoI<sysflag>(collisions, nullptr, muons, nullptr, nullptr);
  }

  void process_Electron_FwdMuon(aod::Collisions const& collisions, aod::EMPrimaryElectrons const& electrons, aod::EMPrimaryMuons const& muons)
  {
    const uint8_t sysflag = kElectron | kFwdMuon;
    selectEoI<sysflag>(collisions, electrons, muons, nullptr, nullptr);
  }

  void process_PCM(aod::Collisions const& collisions, aod::V0PhotonsKF const& v0s)
  {
    const uint8_t sysflag = kPCM;
    selectEoI<sysflag>(collisions, nullptr, nullptr, v0s, nullptr);
  }

  void process_Electron_FwdMuon_PCM(aod::Collisions const& collisions, aod::EMPrimaryElectrons const& electrons, aod::EMPrimaryMuons const& muons, aod::V0PhotonsKF const& v0s)
  {
    const uint8_t sysflag = kElectron | kFwdMuon | kPCM;
    selectEoI<sysflag>(collisions, electrons, muons, v0s, nullptr);
  }

  void process_PCM_ElectronFromDalitz(aod::Collisions const& collisions, aod::V0PhotonsKF const& v0s, aod::EMPrimaryElectronsFromDalitz const& electronsda)
  {
    const uint8_t sysflag = kPCM | kElectronFromDalitz;
    selectEoI<sysflag>(collisions, nullptr, nullptr, v0s, electronsda);
  }

  void processDummy(aod::Collisions const& collisions)
  {
    for (int i = 0; i < collisions.size(); i++) {
      emeoi(true);
    }
  }

  PROCESS_SWITCH(filterEoI, process_Electron, "create filter bit for Electron", false);
  PROCESS_SWITCH(filterEoI, process_FwdMuon, "create filter bit for Forward Muon", false);
  PROCESS_SWITCH(filterEoI, process_PCM, "create filter bit for PCM", false);
  PROCESS_SWITCH(filterEoI, process_Electron_FwdMuon, "create filter bit for Electron, FwdMuon", false);
  PROCESS_SWITCH(filterEoI, process_Electron_FwdMuon_PCM, "create filter bit for Electron, FwdMuon, PCM", false);
  PROCESS_SWITCH(filterEoI, process_PCM_ElectronFromDalitz, "create filter bit for PCM, ElectronFromDalitz", false);
  PROCESS_SWITCH(filterEoI, processDummy, "processDummy", true);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<filterEoI>(cfgc, TaskName{"filter-eoi"})};
}
