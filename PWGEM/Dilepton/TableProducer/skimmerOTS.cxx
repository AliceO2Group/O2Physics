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
// This code produces trigger information. OTS = offline trigger selection.
//    Please write to: daiki.sekihata@cern.ch

#include "PWGEM/Dilepton/DataModel/dileptonTables.h"

#include "Common/Core/TableHelper.h"
#include "Common/Core/Zorro.h"

#include "CCDB/BasicCCDBManager.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

struct skimmerOTS {
  Produces<o2::aod::EMSWTriggerInfosTMP> swtinfo_tmp; // Join aod::Collision later.
  Produces<o2::aod::EMSWTriggerBitsTMP> swtbit_tmp;
  Produces<o2::aod::EMSWTriggerATCountersTMP> swtcounterAT_tmp;
  Produces<o2::aod::EMSWTriggerTOICountersTMP> swtcounterTOI_tmp;

  // CCDB options
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> cfg_swt_names{"cfg_swt_names", "fLMeeIMR,fLMeeHMR", "comma-separated software trigger names"}; // !trigger names have to be pre-registered in dileptonTable.h for bit operation!
  o2::framework::Configurable<std::string> ccdbPathSoftwareTrigger{"ccdbPathSoftwareTrigger", "EventFiltering/Zorro/", "ccdb path for ZORRO objects"};
  Configurable<uint64_t> bcMarginForSoftwareTrigger{"bcMarginForSoftwareTrigger", 100, "Number of BCs of margin for software triggers"};

  std::vector<std::string> swt_names;
  int mRunNumber;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Zorro zorro;
  std::vector<int> mTOIidx;
  uint64_t mNinspectedTVX{0};
  std::vector<uint64_t> mScalers;
  std::vector<uint64_t> mSelections;
  std::vector<int> mTOICounters;
  std::vector<int> mATCounters;

  HistogramRegistry registry{"registry"};

  void init(o2::framework::InitContext&)
  {
    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);
    mRunNumber = 0;

    LOGF(info, "enable software triggers : %s", cfg_swt_names.value.data());
    std::stringstream tokenizer(cfg_swt_names.value);
    std::string token;
    while (std::getline(tokenizer, token, ',')) {
      swt_names.emplace_back(token);
    }

    int nbin = swt_names.size();
    auto hCollisionCounter = registry.add<TH1>("hCollisionCounter", "hCollisionCounter;;Number of collisions", kTH1D, {{nbin + 1, 0.5f, nbin + 1 + 0.5f}});
    hCollisionCounter->GetXaxis()->SetBinLabel(1, "all");
    for (size_t idx = 0; idx < swt_names.size(); idx++) {
      hCollisionCounter->GetXaxis()->SetBinLabel(idx + 2, swt_names[idx].data());
    }

    const int ntrg = static_cast<int>(o2::aod::pwgem::dilepton::swt::swtAliases::kNaliases);
    mNinspectedTVX = 0;
    mScalers.resize(ntrg);
    mSelections.resize(ntrg);
    mTOICounters.resize(ntrg);
    mATCounters.resize(ntrg);
    for (int idx = 0; idx < ntrg; idx++) {
      mTOICounters[idx] = 0;
      mATCounters[idx] = 0;
      mScalers[idx] = 0;
      mSelections[idx] = 0;
    }
  }

  ~skimmerOTS()
  {
    swt_names.clear();
    swt_names.shrink_to_fit();
    mTOICounters.clear();
    mTOICounters.shrink_to_fit();
    mATCounters.clear();
    mATCounters.shrink_to_fit();
  }

  template <typename TBC>
  void initCCDB(TBC const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }

    zorro.setCCDBpath(ccdbPathSoftwareTrigger);
    zorro.setBCtolerance(bcMarginForSoftwareTrigger); // this does nothing.
    mTOIidx = zorro.initCCDB(ccdb.service, bc.runNumber(), bc.timestamp(), cfg_swt_names.value);
    zorro.populateHistRegistry(registry, bc.runNumber());

    mNinspectedTVX = zorro.getInspectedTVX()->GetBinContent(1);
    LOGF(info, "total inspected TVX events = %llu in run number %d", mNinspectedTVX, bc.runNumber());

    for (size_t idx = 0; idx < mTOIidx.size(); idx++) {
      auto swtname = swt_names[idx];
      int emswtId = o2::aod::pwgem::dilepton::swt::aliasLabels.at(swtname);
      mScalers[emswtId] = zorro.getScalers()->GetBinContent(mTOIidx[idx] + 2);
      mSelections[emswtId] = zorro.getSelections()->GetBinContent(mTOIidx[idx] + 2);
      LOGF(info, "Trigger of Interest : index = %d in Zorro, %d in EM, scaler = %llu, selection = %llu", mTOIidx[idx], emswtId, mScalers[emswtId], mSelections[emswtId]);
    }
    swtinfo_tmp(bc.runNumber(), mNinspectedTVX, mScalers, mSelections);
    mRunNumber = bc.runNumber();
  }

  void process(aod::Collisions const& collisions, aod::BCsWithTimestamps const&)
  {
    for (const auto& collision : collisions) {
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>(); // don't use foundBC.
      initCCDB(bc);

      uint16_t trigger_bitmap = 0;
      // uint16_t analyzed_bitmap = 0;
      // uint16_t analyzedToI_bitmap = 0;
      registry.fill(HIST("hCollisionCounter"), 1); // all

      if (zorro.isSelected(bc.globalBC(), bcMarginForSoftwareTrigger)) { // triggered event
        auto swt_bitset = zorro.getLastResult();                         // this has to be called after zorro::isSelected, or simply call zorro.fetch
        auto TOIcounters = zorro.getTOIcounters();                       // this has to be called after zorro::isSelected, or simply call zorro.fetch
        auto ATcounters = zorro.getATcounters();                         // this has to be called after zorro::isSelected, or simply call zorro.fetch

        // LOGF(info, "swt_bitset.to_string().c_str() = %s", swt_bitset.to_string().c_str());
        for (size_t idx = 0; idx < mTOIidx.size(); idx++) {
          if (swt_bitset.test(mTOIidx[idx])) {
            auto swtname = swt_names[idx];
            int emswtId = o2::aod::pwgem::dilepton::swt::aliasLabels.at(swtname);
            trigger_bitmap |= BIT(emswtId);
            // LOGF(info, "swtname = %s is fired. swt index in original swt table = %d, swt index for EM table = %d", swtname.data(), mTOIidx[idx], o2::aod::pwgem::dilepton::swt::aliasLabels.at(swtname));
            registry.fill(HIST("hCollisionCounter"), idx + 2); // fired trigger

            // LOGF(info, "ATcounters[mTOIidx[idx]] = %d, TOIcounters[idx] = %d", ATcounters[mTOIidx[idx]], TOIcounters[idx]);

            while (ATcounters[mTOIidx[idx]] > mATCounters[emswtId]) {
              mATCounters[emswtId]++;
              swtcounterAT_tmp(BIT(emswtId));
            }

            while (TOIcounters[idx] > mTOICounters[emswtId]) {
              mTOICounters[emswtId]++; // always incremented by 1 in zorro!!
              swtcounterTOI_tmp(BIT(emswtId));
            }

            // LOGF(info, "collision.globalIndex() = %d, bc.globalBC() = %llu, mTOICounters[%d] = %d, mATcounters[%d] = %d", collision.globalIndex(), bc.globalBC(), emswtId, mTOICounters[emswtId], emswtId, mATCounters[emswtId]);
          }
        } // end of TOI loop
      }
      swtbit_tmp(trigger_bitmap);
    } // end of collision loop
  } // end of process
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<skimmerOTS>(cfgc, TaskName{"skimmer-ots"})};
}
