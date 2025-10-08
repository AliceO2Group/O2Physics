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

/// \file EventHistograms.h
/// \brief header file for event histograms
/// \author daiki.sekihata@cern.ch

#ifndef PWGEM_PHOTONMESON_UTILS_EVENTHISTOGRAMS_H_
#define PWGEM_PHOTONMESON_UTILS_EVENTHISTOGRAMS_H_

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/TriggerAliases.h"

#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>

#include <TH1.h>

#include <string_view>

namespace o2::aod::pwgem::photonmeson::utils::eventhistogram
{
inline void addEventHistograms(o2::framework::HistogramRegistry* fRegistry)
{
  // event info
  auto hCollisionCounter = fRegistry->add<TH1>("Event/before/hCollisionCounter", "collision counter;;Number of events", o2::framework::kTH1D, {{12, 0.5, 12.5}}, false);
  hCollisionCounter->GetXaxis()->SetBinLabel(1, "all");
  hCollisionCounter->GetXaxis()->SetBinLabel(2, "No TF border");
  hCollisionCounter->GetXaxis()->SetBinLabel(3, "No ITS ROF border");
  hCollisionCounter->GetXaxis()->SetBinLabel(4, "No Same Bunch Pileup");
  hCollisionCounter->GetXaxis()->SetBinLabel(5, "Is Vertex ITSTPC");
  hCollisionCounter->GetXaxis()->SetBinLabel(6, "Is Good Zvtx FT0vsPV");
  hCollisionCounter->GetXaxis()->SetBinLabel(7, "FT0AND");
  hCollisionCounter->GetXaxis()->SetBinLabel(8, "sel8");
  hCollisionCounter->GetXaxis()->SetBinLabel(9, "|Z_{vtx}| < 10 cm");
  hCollisionCounter->GetXaxis()->SetBinLabel(10, "EMC MB Readout");
  hCollisionCounter->GetXaxis()->SetBinLabel(11, "EMC L0 Triggered");
  hCollisionCounter->GetXaxis()->SetBinLabel(12, "accepted");

  const o2::framework::AxisSpec axis_cent_ft0m{{0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
                                                20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39,
                                                40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59,
                                                60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79,
                                                80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110},
                                               "centrality FT0M (%)"};

  const o2::framework::AxisSpec axis_cent_ft0a{{0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
                                                20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39,
                                                40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59,
                                                60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79,
                                                80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110},
                                               "centrality FT0A (%)"};

  const o2::framework::AxisSpec axis_cent_ft0c{{0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
                                                20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39,
                                                40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59,
                                                60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79,
                                                80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110},
                                               "centrality FT0C (%)"};

  fRegistry->add("Event/before/hZvtx", "vertex z; Z_{vtx} (cm)", o2::framework::kTH1F, {{100, -50, +50}}, false);
  fRegistry->add("Event/before/hMultNTracksPV", "hMultNTracksPV; N_{track} to PV", o2::framework::kTH1F, {{6001, -0.5, 6000.5}}, false);
  fRegistry->add("Event/before/hMultNTracksPVeta1", "hMultNTracksPVeta1; N_{track} to PV", o2::framework::kTH1F, {{6001, -0.5, 6000.5}}, false);
  fRegistry->add("Event/before/hMultFT0", "hMultFT0;mult. FT0A;mult. FT0C", o2::framework::kTH2F, {{200, 0, 200000}, {60, 0, 60000}}, false);
  fRegistry->add("Event/before/hCentFT0A", "hCentFT0A;centrality FT0A (%)", o2::framework::kTH1F, {{axis_cent_ft0a}}, false);
  fRegistry->add("Event/before/hCentFT0C", "hCentFT0C;centrality FT0C (%)", o2::framework::kTH1F, {{axis_cent_ft0c}}, false);
  fRegistry->add("Event/before/hCentFT0M", "hCentFT0M;centrality FT0M (%)", o2::framework::kTH1F, {{axis_cent_ft0m}}, false);
  fRegistry->add("Event/before/hCentFT0CvsMultNTracksPV", "hCentFT0CvsMultNTracksPV;centrality FT0C (%);N_{track} to PV", o2::framework::kTH2F, {{110, 0, 110}, {500, 0, 5000}}, false);
  fRegistry->add("Event/before/hMultFT0CvsMultNTracksPV", "hMultFT0CvsMultNTracksPV;mult. FT0C;N_{track} to PV", o2::framework::kTH2F, {{60, 0, 60000}, {500, 0, 5000}}, false);
  fRegistry->add("Event/before/hMultFT0CvsOccupancy", "hMultFT0CvsOccupancy;mult. FT0C;N_{track} in time range", o2::framework::kTH2F, {{60, 0, 60000}, {2000, 0, 20000}}, false);
  fRegistry->addClone("Event/before/", "Event/after/");
}

template <const int ev_id, typename TCollision>
void fillEventInfo(o2::framework::HistogramRegistry* fRegistry, TCollision const& collision, float weight = 1.)
{
  const float maxZ = 10.f;
  static constexpr std::string_view kEventTypes[2] = {"before/", "after/"};
  fRegistry->fill(HIST("Event/") + HIST(kEventTypes[ev_id]) + HIST("hCollisionCounter"), 1.0, weight);
  if (collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
    fRegistry->fill(HIST("Event/") + HIST(kEventTypes[ev_id]) + HIST("hCollisionCounter"), 2.0, weight);
  }
  if (collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
    fRegistry->fill(HIST("Event/") + HIST(kEventTypes[ev_id]) + HIST("hCollisionCounter"), 3.0, weight);
  }
  if (collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
    fRegistry->fill(HIST("Event/") + HIST(kEventTypes[ev_id]) + HIST("hCollisionCounter"), 4.0, weight);
  }
  if (collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
    fRegistry->fill(HIST("Event/") + HIST(kEventTypes[ev_id]) + HIST("hCollisionCounter"), 5.0, weight);
  }
  if (collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
    fRegistry->fill(HIST("Event/") + HIST(kEventTypes[ev_id]) + HIST("hCollisionCounter"), 6.0, weight);
  }
  if (collision.selection_bit(o2::aod::evsel::kIsTriggerTVX)) {
    fRegistry->fill(HIST("Event/") + HIST(kEventTypes[ev_id]) + HIST("hCollisionCounter"), 7.0, weight);
  }
  if (collision.sel8()) {
    fRegistry->fill(HIST("Event/") + HIST(kEventTypes[ev_id]) + HIST("hCollisionCounter"), 8.0, weight);
  }
  if (std::abs(collision.posZ()) < maxZ) {
    fRegistry->fill(HIST("Event/") + HIST(kEventTypes[ev_id]) + HIST("hCollisionCounter"), 9.0, weight);
  }
  fRegistry->fill(HIST("Event/") + HIST(kEventTypes[ev_id]) + HIST("hZvtx"), collision.posZ(), weight);
  if (collision.alias_bit(kTVXinEMC)) {
    fRegistry->fill(HIST("Event/") + HIST(kEventTypes[ev_id]) + HIST("hCollisionCounter"), 10.0, weight);
  }
  if (collision.alias_bit(kEMC7) || collision.alias_bit(kDMC7)) {
    fRegistry->fill(HIST("Event/") + HIST(kEventTypes[ev_id]) + HIST("hCollisionCounter"), 11.0, weight);
  }

  fRegistry->fill(HIST("Event/") + HIST(kEventTypes[ev_id]) + HIST("hMultNTracksPV"), collision.multNTracksPV(), weight);
  fRegistry->fill(HIST("Event/") + HIST(kEventTypes[ev_id]) + HIST("hMultNTracksPVeta1"), collision.multNTracksPVeta1(), weight);
  fRegistry->fill(HIST("Event/") + HIST(kEventTypes[ev_id]) + HIST("hMultFT0"), collision.multFT0A(), collision.multFT0C(), weight);
  fRegistry->fill(HIST("Event/") + HIST(kEventTypes[ev_id]) + HIST("hCentFT0A"), collision.centFT0A(), weight);
  fRegistry->fill(HIST("Event/") + HIST(kEventTypes[ev_id]) + HIST("hCentFT0C"), collision.centFT0C(), weight);
  fRegistry->fill(HIST("Event/") + HIST(kEventTypes[ev_id]) + HIST("hCentFT0M"), collision.centFT0M(), weight);
  fRegistry->fill(HIST("Event/") + HIST(kEventTypes[ev_id]) + HIST("hCentFT0CvsMultNTracksPV"), collision.centFT0C(), collision.multNTracksPV(), weight);
  fRegistry->fill(HIST("Event/") + HIST(kEventTypes[ev_id]) + HIST("hMultFT0CvsMultNTracksPV"), collision.multFT0C(), collision.multNTracksPV(), weight);
  fRegistry->fill(HIST("Event/") + HIST(kEventTypes[ev_id]) + HIST("hMultFT0CvsOccupancy"), collision.multFT0C(), collision.trackOccupancyInTimeRange(), weight);
}

} // namespace o2::aod::pwgem::photonmeson::utils::eventhistogram

#endif // PWGEM_PHOTONMESON_UTILS_EVENTHISTOGRAMS_H_
