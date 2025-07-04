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
  };
  Produces<o2::aod::EMEoIs> emeoi;
  Configurable<int> minNElectrons{"minNElectrons", 1, "min number of e+ and e- at midrapidity"};
  Configurable<int> minNMuons{"minNMuons", 1, "min number of mu+ and mu- at forward rapidity"};

  HistogramRegistry fRegistry{"output"};
  void init(o2::framework::InitContext&)
  {
    auto hEventCounter = fRegistry.add<TH1>("hEventCounter", "hEventCounter", kTH1D, {{5, 0.5f, 5.5f}});
    hEventCounter->GetXaxis()->SetBinLabel(1, "all");
    hEventCounter->GetXaxis()->SetBinLabel(2, "event with electron");
    hEventCounter->GetXaxis()->SetBinLabel(3, "event with forward muon");
    hEventCounter->GetXaxis()->SetBinLabel(4, "event with electron or forward muon");
    hEventCounter->GetXaxis()->SetBinLabel(5, "event with electron and forward muon");
  }

  SliceCache cache;
  Preslice<aod::EMPrimaryElectrons> perCollision_el = aod::emprimaryelectron::collisionId;
  Preslice<aod::EMPrimaryMuons> perCollision_mu = aod::emprimarymuon::collisionId;

  template <uint8_t system, typename TCollisions, typename TElectrons, typename TMuons>
  void selectEoI(TCollisions const& collisions, TElectrons const& electrons, TMuons const& muons)
  {
    for (const auto& collision : collisions) {
      bool does_electron_exist = false;
      bool does_fwdmuon_exist = false;
      fRegistry.fill(HIST("hEventCounter"), 1);

      if constexpr (static_cast<bool>(system & kElectron)) {
        auto electrons_coll = electrons.sliceBy(perCollision_el, collision.globalIndex());
        if (electrons_coll.size() >= minNElectrons) {
          does_electron_exist = true;
          fRegistry.fill(HIST("hEventCounter"), 2);
        }
      }
      if constexpr (static_cast<bool>(system & kFwdMuon)) {
        auto muons_coll = muons.sliceBy(perCollision_mu, collision.globalIndex());
        if (muons_coll.size() >= minNMuons) {
          does_fwdmuon_exist = true;
          fRegistry.fill(HIST("hEventCounter"), 3);
        }
      }

      if (does_electron_exist || does_fwdmuon_exist) {
        fRegistry.fill(HIST("hEventCounter"), 4);
      }
      if (does_electron_exist && does_fwdmuon_exist) {
        fRegistry.fill(HIST("hEventCounter"), 5);
      }

      emeoi(does_electron_exist || does_fwdmuon_exist);

    } // end of collision loop

  } // end of selectEoI

  void process_Electron(aod::Collisions const& collisions, aod::EMPrimaryElectrons const& electrons)
  {
    const uint8_t sysflag = kElectron;
    selectEoI<sysflag>(collisions, electrons, nullptr);
  }

  void process_FwdMuon(aod::Collisions const& collisions, aod::EMPrimaryMuons const& muons)
  {
    const uint8_t sysflag = kFwdMuon;
    selectEoI<sysflag>(collisions, nullptr, muons);
  }

  void process_Electron_FwdMuon(aod::Collisions const& collisions, aod::EMPrimaryElectrons const& electrons, aod::EMPrimaryMuons const& muons)
  {
    const uint8_t sysflag = kElectron | kFwdMuon;
    selectEoI<sysflag>(collisions, electrons, muons);
  }

  void processDummy(aod::Collisions const& collisions)
  {
    for (int i = 0; i < collisions.size(); i++) {
      emeoi(true);
    }
  }

  PROCESS_SWITCH(filterEoI, process_Electron, "create filter bit for Electron", false);
  PROCESS_SWITCH(filterEoI, process_FwdMuon, "create filter bit for Forward Muon", false);
  PROCESS_SWITCH(filterEoI, process_Electron_FwdMuon, "create filter bit for Electron, FwdMuon", false);
  PROCESS_SWITCH(filterEoI, processDummy, "processDummy", true);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<filterEoI>(cfgc, TaskName{"filter-eoi"})};
}
