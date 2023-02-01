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
// O2 includes

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Common/DataModel/CaloClusters.h"
#include "DataFormatsPHOS/TriggerRecord.h"

#include "../filterTables.h"

#include "Framework/HistogramRegistry.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct photonFilter {
  static constexpr int nTrigs{4};
  enum trigs { kPhot,
               kEl,
               kPair,
               kNbar };

  Produces<aod::PhotFilters> tags;

  Configurable<float> ePhot{"ePhot", 2., "Minimal photon energy (GeV)"};
  Configurable<float> eEl{"eEl", 2., "Minimal electron energy (GeV)"};
  Configurable<float> ePair{"ePair", 0.5, "Minimal photon pair mass (GeV)"};

  HistogramRegistry mHistManager{"events", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  void init(o2::framework::InitContext&)
  {
    auto scalers{std::get<std::shared_ptr<TH1>>(mHistManager.add("fProcessedEvents", "Number of filtered events", HistType::kTH1F, {{10, -0.5, 9.5}}))};
    scalers->GetXaxis()->SetBinLabel(1, "Processed events");
    scalers->GetXaxis()->SetBinLabel(2, "PHOS photon");
    scalers->GetXaxis()->SetBinLabel(3, "PHOS electron");
    scalers->GetXaxis()->SetBinLabel(4, "PHOS pair");
    scalers->GetXaxis()->SetBinLabel(5, "PHOS nbar");
    scalers->GetXaxis()->SetBinLabel(6, "PHOS photon & electron");
    scalers->GetXaxis()->SetBinLabel(7, "PHOS photon & pair");
    scalers->GetXaxis()->SetBinLabel(8, "events with PHOS");
    scalers->GetXaxis()->SetBinLabel(9, "several BC collision");
  }

  void process(o2::aod::CaloClusters const& clusters,
               o2::aod::CPVClusters const& cpvs)
  {
    // Process all clusters per TF and fill corresponding collisions
    bool keepEvent[nTrigs]{false};

    mHistManager.fill(HIST("fProcessedEvents"), 1);

    // to trace switching to new BC
    o2::InteractionRecord prevIr;
    int prevColId = -2;
    o2::InteractionRecord ir;
    // PHOS part
    int nPHOSclu = 0;
    for (const auto& clu : clusters) {
      if (clu.caloType() != 0)
        continue;
      nPHOSclu++;
      if (prevColId == -2) { // first cluster
        prevColId = clu.colId();
        prevIr.setFromLong(clu.bc().globalBC());
      }
      ir.setFromLong(clu.bc().globalBC());
      if (ir != prevIr) {                                // switch to next BC
        mHistManager.fill(HIST("fProcessedEvents"), 9.); // Double BCs???
        // reset everething
        for (int i = 0; i < nTrigs; i++) {
          keepEvent[i] = false;
        }
        prevIr = ir;
        prevColId = clu.colId();
      }

      // Scan current cluster
      //  photons
      keepEvent[kPhot] |= (clu.e() > ePhot);
      // charged clusters above threshold
      keepEvent[kEl] |= (clu.trackdist() < 2. && clu.e() > eEl); // 2: Distance to CPV cluster in sigmas
      // antineutrons
      keepEvent[kNbar] |= (clu.ncell() > 2 && clu.m02() > 0.2 && clu.e() > 0.7 && clu.trackdist() > 2.) &&
                          ((clu.e() < 2. && clu.m02() > 4.5 - clu.m20()) ||
                           (clu.e() > 2. && clu.m02() > 4. - clu.m20()));

      // inv mass
      auto clu2 = clu;
      ++clu2;
      for (; !keepEvent[kPair] && clu2 != clusters.end() && clu.bc().globalBC() == clu2.bc().globalBC(); clu2++) { // scan same event only
        // cluster selections
        if (clu.trackdist() < 1. || clu2.trackdist() < 1.) { // select neutral clusters. Disp, Ncell cuts?
          continue;
        }
        double m = pow(clu.e() + clu2.e(), 2) - pow(clu.px() + clu2.px(), 2) -
                   pow(clu.py() + clu2.py(), 2) - pow(clu.pz() + clu2.pz(), 2);
        if (m > ePair * ePair) {
          keepEvent[kPair] |= true;
          break;
        }
      }
    }

    // Normally, one BC should be processed, fill scalers here
    if (nPHOSclu) {
      mHistManager.fill(HIST("fProcessedEvents"), 8.);
    }
    // Can not fill with variable, have to fill manually
    if (keepEvent[kPhot]) {
      mHistManager.fill(HIST("fProcessedEvents"), 2.);
      if (keepEvent[kEl]) {
        mHistManager.fill(HIST("fProcessedEvents"), 6.);
      }
      if (keepEvent[kPair]) {
        mHistManager.fill(HIST("fProcessedEvents"), 7.);
      }
    }
    if (keepEvent[kEl]) {
      mHistManager.fill(HIST("fProcessedEvents"), 3.);
    }
    if (keepEvent[kPair]) {
      mHistManager.fill(HIST("fProcessedEvents"), 4.);
    }
    if (keepEvent[kNbar]) {
      mHistManager.fill(HIST("fProcessedEvents"), 5.);
    }
    // fill
    tags(keepEvent[kPhot], keepEvent[kEl], keepEvent[kPair], keepEvent[kNbar]);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{
  return WorkflowSpec{
    adaptAnalysisTask<photonFilter>(cfg)};
}
