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
#include "EventFiltering/Zorro.h"

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
  Produces<o2::aod::EMSWTriggerInfosTMP> swt_tmp;

  // CCDB options
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> cfg_swt_names{"cfg_swt_names", "fHighTrackMult,fHighFt0Mult", "comma-separated software trigger names"}; // !trigger names have to be pre-registered in dileptonTable.h for bit operation!

  std::vector<std::string> swt_names;
  int mRunNumber;
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  HistogramRegistry registry{"registry"};
  void init(o2::framework::InitContext&)
  {
    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    LOGF(info, "enable software triggers : %s", cfg_swt_names.value.data());
    std::stringstream tokenizer(cfg_swt_names.value);
    std::string token;
    while (std::getline(tokenizer, token, ',')) {
      swt_names.emplace_back(token);
    }

    const int nbin = swt_names.size();
    auto hEventCounter = registry.add<TH1>("hEventCounter", "hEventCounter;;Number of Events", kTH1D, {{nbin + 1, 0.5f, nbin + 1 + 0.5f}});
    hEventCounter->GetXaxis()->SetBinLabel(1, "all");
    for (int idx = 0; idx < nbin; idx++) {
      hEventCounter->GetXaxis()->SetBinLabel(idx + 2, swt_names[idx].data());
    }

    registry.add("hNInspectedTVX", "N inspected TVX;run number;N_{TVX}", kTProfile, {{80000, 520000.5, 600000.5}}, true);
  }

  ~skimmerOTS()
  {
    swt_names.clear();
    swt_names.shrink_to_fit();
  }

  Zorro zorro;
  std::vector<int> mTOIidx;
  uint64_t mNinspectedTVX{0};

  template <typename TBC>
  void initCCDB(TBC const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }

    mTOIidx = zorro.initCCDB(ccdb.service, bc.runNumber(), bc.timestamp(), cfg_swt_names.value);
    for (auto& idx : mTOIidx) {
      LOGF(info, "Trigger of Interest : index = %d", idx);
    }
    mNinspectedTVX = zorro.getInspectedTVX()->GetBinContent(1);
    LOGF(info, "total inspected TVX events = %d in run number %d", mNinspectedTVX, bc.runNumber());
    registry.fill(HIST("hNInspectedTVX"), bc.runNumber(), mNinspectedTVX);

    mRunNumber = bc.runNumber();
  }

  using MyBCs = soa::Join<aod::BCsWithTimestamps, aod::BcSels>;
  using MyCollisions = soa::Join<aod::Collisions, aod::EvSels>;

  void process(MyCollisions const& collisions, MyBCs const&)
  {
    for (auto& collision : collisions) {
      auto bc = collision.template bc_as<MyBCs>(); // don't use foundBC.
      initCCDB(bc);

      uint16_t trigger_bitmap = 0;
      registry.fill(HIST("hEventCounter"), 1);   // all
      if (zorro.isSelected(bc.globalBC())) {     // triggered event
        auto swt_bitset = zorro.getLastResult(); // this has to be called after zorro::isSelected, or simply call zorro.fetch
        // LOGF(info, "swt_bitset.to_string().c_str() = %s", swt_bitset.to_string().c_str());
        for (size_t idx = 0; idx < mTOIidx.size(); idx++) {
          if (swt_bitset.test(mTOIidx[idx])) {
            auto swtname = swt_names[idx];
            trigger_bitmap |= BIT(o2::aod::pwgem::dilepton::swt::aliasLabels.at(swtname));
            // LOGF(info, "swtname = %s is fired. swt index in original swt table = %d, swt index for EM table = %d", swtname.data(), mTOIidx[idx], o2::aod::pwgem::dilepton::swt::aliasLabels.at(swtname));
            registry.fill(HIST("hEventCounter"), idx + 2); // fired trigger
          }
        }
      }
      // LOGF(info, "trigger_bitmap = %d, mNinspectedTVX = %d", trigger_bitmap, mNinspectedTVX);
      swt_tmp(trigger_bitmap, mNinspectedTVX);
    } // end of collision loop
  }
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<skimmerOTS>(cfgc, TaskName{"skimmer-ots"})};
}
