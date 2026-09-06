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

#include <Framework/Concepts.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>

#include <TH1.h>

#include <array>
#include <string_view>

namespace o2::aod::pwgem::photonmeson::utils::eventhistogram
{
inline void addEventHistograms(o2::framework::HistogramRegistry* fRegistry, bool useWeight = false)
{
  // event info
  auto hCollisionCounter = fRegistry->add<TH1>("Event/before/hCollisionCounter", "collision counter;;Number of events", o2::framework::HistType::kTH1D, {{12, 0.5, 12.5}}, useWeight);
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

  fRegistry->add("Event/before/hZvtx", "vertex z; Z_{vtx} (cm)", o2::framework::HistType::kTH1D, {{220, -11, +11}}, useWeight);
  fRegistry->add("Event/before/hMultNTracksPV", "hMultNTracksPV; N_{track} to PV", o2::framework::HistType::kTH1F, {{6001, -0.5, 6000.5}}, useWeight);
  fRegistry->add("Event/before/hMultNTracksPVeta1", "hMultNTracksPVeta1; N_{track} to PV", o2::framework::HistType::kTH1F, {{6001, -0.5, 6000.5}}, useWeight);
  fRegistry->add("Event/before/hMultFT0", "hMultFT0;mult. FT0A;mult. FT0C", o2::framework::HistType::kTH2F, {{200, 0, 200000}, {60, 0, 60000}}, useWeight);
  fRegistry->add("Event/before/hCentFT0A", "hCentFT0A;centrality FT0A (%)", o2::framework::HistType::kTH1F, {{100, 0., 100.}}, useWeight);
  fRegistry->add("Event/before/hCentFT0C", "hCentFT0C;centrality FT0C (%)", o2::framework::HistType::kTH1F, {{100, 0., 100.}}, useWeight);
  fRegistry->add("Event/before/hCentFT0M", "hCentFT0M;centrality FT0M (%)", o2::framework::HistType::kTH1F, {{100, 0., 100.}}, useWeight);
  fRegistry->add("Event/before/hCentFT0CvsMultNTracksPV", "hCentFT0CvsMultNTracksPV;centrality FT0C (%);N_{track} to PV", o2::framework::HistType::kTH2F, {{110, 0, 110}, {500, 0, 5000}}, useWeight);
  fRegistry->add("Event/before/hMultFT0CvsMultNTracksPV", "hMultFT0CvsMultNTracksPV;mult. FT0C;N_{track} to PV", o2::framework::HistType::kTH2F, {{60, 0, 60000}, {500, 0, 5000}}, useWeight);
  fRegistry->add("Event/before/hMultFT0CvsOccupancy", "hMultFT0CvsOccupancy;mult. FT0C;N_{track} in time range", o2::framework::HistType::kTH2F, {{60, 0, 60000}, {500, 0, 10000}}, useWeight);
  fRegistry->addClone("Event/before/", "Event/after/");
}

template <const int ev_id, o2::soa::is_iterator TCollision>
void fillEventInfo(o2::framework::HistogramRegistry* fRegistry, TCollision const& collision, float weight = 1.)
{
  const float maxZ = 10.f;
  static constexpr std::array<std::string_view, 2> EventTypes = {"before/", "after/"};
  fRegistry->fill(HIST("Event/") + HIST(EventTypes[ev_id]) + HIST("hCollisionCounter"), 1.0, weight);
  if (collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
    fRegistry->fill(HIST("Event/") + HIST(EventTypes[ev_id]) + HIST("hCollisionCounter"), 2.0, weight);
  }
  if (collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
    fRegistry->fill(HIST("Event/") + HIST(EventTypes[ev_id]) + HIST("hCollisionCounter"), 3.0, weight);
  }
  if (collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
    fRegistry->fill(HIST("Event/") + HIST(EventTypes[ev_id]) + HIST("hCollisionCounter"), 4.0, weight);
  }
  if (collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
    fRegistry->fill(HIST("Event/") + HIST(EventTypes[ev_id]) + HIST("hCollisionCounter"), 5.0, weight);
  }
  if (collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
    fRegistry->fill(HIST("Event/") + HIST(EventTypes[ev_id]) + HIST("hCollisionCounter"), 6.0, weight);
  }
  if (collision.selection_bit(o2::aod::evsel::kIsTriggerTVX)) {
    fRegistry->fill(HIST("Event/") + HIST(EventTypes[ev_id]) + HIST("hCollisionCounter"), 7.0, weight);
  }
  if (collision.sel8()) {
    fRegistry->fill(HIST("Event/") + HIST(EventTypes[ev_id]) + HIST("hCollisionCounter"), 8.0, weight);
  }
  if (std::abs(collision.posZ()) < maxZ) {
    fRegistry->fill(HIST("Event/") + HIST(EventTypes[ev_id]) + HIST("hCollisionCounter"), 9.0, weight);
  }
  fRegistry->fill(HIST("Event/") + HIST(EventTypes[ev_id]) + HIST("hZvtx"), collision.posZ(), weight);
  if (collision.alias_bit(kTVXinEMC)) {
    fRegistry->fill(HIST("Event/") + HIST(EventTypes[ev_id]) + HIST("hCollisionCounter"), 10.0, weight);
  }
  if (collision.alias_bit(kEMC7) || collision.alias_bit(kDMC7)) {
    fRegistry->fill(HIST("Event/") + HIST(EventTypes[ev_id]) + HIST("hCollisionCounter"), 11.0, weight);
  }

  fRegistry->fill(HIST("Event/") + HIST(EventTypes[ev_id]) + HIST("hMultNTracksPV"), collision.multNTracksPV(), weight);
  fRegistry->fill(HIST("Event/") + HIST(EventTypes[ev_id]) + HIST("hMultNTracksPVeta1"), collision.multNTracksPVeta1(), weight);
  fRegistry->fill(HIST("Event/") + HIST(EventTypes[ev_id]) + HIST("hMultFT0"), collision.multFT0A(), collision.multFT0C(), weight);
  fRegistry->fill(HIST("Event/") + HIST(EventTypes[ev_id]) + HIST("hCentFT0A"), collision.centFT0A(), weight);
  fRegistry->fill(HIST("Event/") + HIST(EventTypes[ev_id]) + HIST("hCentFT0C"), collision.centFT0C(), weight);
  fRegistry->fill(HIST("Event/") + HIST(EventTypes[ev_id]) + HIST("hCentFT0M"), collision.centFT0M(), weight);
  fRegistry->fill(HIST("Event/") + HIST(EventTypes[ev_id]) + HIST("hCentFT0CvsMultNTracksPV"), collision.centFT0C(), collision.multNTracksPV(), weight);
  fRegistry->fill(HIST("Event/") + HIST(EventTypes[ev_id]) + HIST("hMultFT0CvsMultNTracksPV"), collision.multFT0C(), collision.multNTracksPV(), weight);
  fRegistry->fill(HIST("Event/") + HIST(EventTypes[ev_id]) + HIST("hMultFT0CvsOccupancy"), collision.multFT0C(), collision.trackOccupancyInTimeRange(), weight);
}

} // namespace o2::aod::pwgem::photonmeson::utils::eventhistogram

#endif // PWGEM_PHOTONMESON_UTILS_EVENTHISTOGRAMS_H_
