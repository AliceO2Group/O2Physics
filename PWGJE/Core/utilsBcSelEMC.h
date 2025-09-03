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

/// \file utilsBcSelEMC.h
/// \brief BC selection utilities for EMCal QA
/// \author Marvin Hemmer <marvin.hemmer@cern.ch>, Goethe-University

#ifndef PWGJE_CORE_UTILSBCSELEMC_H_
#define PWGJE_CORE_UTILSBCSELEMC_H_

#include "Common/CCDB/EventSelectionParams.h"

#include <Framework/Configurable.h>
#include <Framework/DeviceSpec.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>

#include <TH1.h>

#include <Rtypes.h>

#include <cstddef>
#include <cstdint>
#include <memory> // std::shared_ptr
#include <string> // std::string

namespace o2::emc_evsel
{
// event rejection types
enum EventRejection {
  None = 0,
  Trigger,
  TvxTrigger,
  TimeFrameBorderCut,
  ItsRofBorderCut,
  NEventRejection
};

inline o2::framework::AxisSpec axisEvents = {EventRejection::NEventRejection, -0.5f, +EventRejection::NEventRejection - 0.5f, ""};

/// \brief Function to put labels on monitoring histogram
/// \param hRejection monitoring histogram
template <typename Histo>
void setEventRejectionLabels(Histo& hRejection)
{
  // Puts labels on the bc monitoring histogram.
  hRejection->GetXaxis()->SetBinLabel(EventRejection::None + 1, "All");
  hRejection->GetXaxis()->SetBinLabel(EventRejection::Trigger + 1, "Sel8");
  hRejection->GetXaxis()->SetBinLabel(EventRejection::TvxTrigger + 1, "TVX Trigger");
  hRejection->GetXaxis()->SetBinLabel(EventRejection::TimeFrameBorderCut + 1, "TF border");
  hRejection->GetXaxis()->SetBinLabel(EventRejection::ItsRofBorderCut + 1, "ITS ROF border");
}

struct EMCEventSelection : o2::framework::ConfigurableGroup {
  std::string prefix = "emcBcSel"; // JSON group name
  // event selection parameters (in chronological order of application)
  o2::framework::Configurable<bool> useSel8Trigger{"useSel8Trigger", true, "Apply the sel8 event selection"};
  o2::framework::Configurable<int> triggerClass{"triggerClass", -1, "Trigger class different from sel8 (e.g. kINT7 for Run2) used only if useSel8Trigger is false"};
  o2::framework::Configurable<bool> useTvxTrigger{"useTvxTrigger", true, "Apply TVX trigger sel"};
  o2::framework::Configurable<bool> useTimeFrameBorderCut{"useTimeFrameBorderCut", true, "Apply TF border cut"};
  o2::framework::Configurable<bool> useItsRofBorderCut{"useItsRofBorderCut", true, "Apply ITS ROF border cut"};
  // histogram names
  static constexpr char NameHistBCs[] = "hBCs";

  std::shared_ptr<TH1> hBCs;

  int currentRun{-1};

  /// \brief Adds bc monitoring histograms in the histogram registry.
  /// \param registry reference to the histogram registry
  void addHistograms(o2::framework::HistogramRegistry& registry)
  {
    hBCs = registry.add<TH1>(NameHistBCs, "EMC event counter;;# of accepted collisions", {o2::framework::HistType::kTH1D, {axisEvents}});
    setEventRejectionLabels(hBCs);
  }

  /// \brief Applies event selection.
  /// \tparam useEvSel use information from the EvSel table
  /// \param bc bc to test against the selection criteria
  /// \return bitmask with the event selection criteria not satisfied by the bc
  template <bool useEvSel, typename Bc>
  uint16_t getEMCCollisionRejectionMask(const Bc& bc)
  {
    uint16_t rejectionMask{0}; // 16 bits, in case new ev. selections will be added

    if constexpr (useEvSel) {
      /// trigger condition
      bool sel8 = bc.selection_bit(o2::aod::evsel::kIsTriggerTVX) && bc.selection_bit(o2::aod::evsel::kNoTimeFrameBorder) && bc.selection_bit(o2::aod::evsel::kNoITSROFrameBorder);
      if ((useSel8Trigger && !sel8) || (!useSel8Trigger && triggerClass > -1 && !bc.alias_bit(triggerClass))) {
        SETBIT(rejectionMask, EventRejection::Trigger);
      }
      /// TVX trigger selection
      if (useTvxTrigger && !bc.selection_bit(o2::aod::evsel::kIsTriggerTVX)) {
        SETBIT(rejectionMask, EventRejection::TvxTrigger);
      }
      /// time frame border cut
      if (useTimeFrameBorderCut && !bc.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
        SETBIT(rejectionMask, EventRejection::TimeFrameBorderCut);
      }
      /// ITS rof border cut
      if (useItsRofBorderCut && !bc.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
        SETBIT(rejectionMask, EventRejection::ItsRofBorderCut);
      }
    }
    return rejectionMask;
  }

  /// \brief Fills histograms for monitoring event selections satisfied by the bc.
  /// \param rejectionMask bitmask storing the info about which ev. selections are not satisfied by the bc
  void fillHistograms(const uint16_t rejectionMask)
  {
    hBCs->Fill(EventRejection::None);
    for (std::size_t reason = 1; reason < EventRejection::NEventRejection; reason++) {
      if (TESTBIT(rejectionMask, reason)) {
        return;
      }
      hBCs->Fill(reason);
    }
  }
};

struct EMCEventSelectionMc {
  // event selection parameters (in chronological order of application)
  bool useSel8Trigger{false};       // Apply the Sel8 selection
  bool useTvxTrigger{false};        // Apply the TVX trigger
  bool useTimeFrameBorderCut{true}; // Apply TF border cut
  bool useItsRofBorderCut{false};   // Apply the ITS RO frame border cut

  void configureFromDevice(const o2::framework::DeviceSpec& device)
  {
    for (const auto& option : device.options) {
      if (option.name.compare("emcEvSel.useSel8Trigger") == 0) {
        useSel8Trigger = option.defaultValue.get<bool>();
      } else if (option.name.compare("emcEvSel.useTvxTrigger") == 0) {
        useTvxTrigger = option.defaultValue.get<bool>();
      } else if (option.name.compare("emcEvSel.useTimeFrameBorderCut") == 0) {
        useTimeFrameBorderCut = option.defaultValue.get<bool>();
      } else if (option.name.compare("emcEvSel.useItsRofBorderCut") == 0) {
        useItsRofBorderCut = option.defaultValue.get<bool>();
      }
    }
  }
};
} // namespace o2::emc_evsel

#endif // PWGJE_CORE_UTILSBCSELEMC_H_
