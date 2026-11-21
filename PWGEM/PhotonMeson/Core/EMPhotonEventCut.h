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
// Class for em photon event selection
//

#ifndef PWGEM_PHOTONMESON_CORE_EMPHOTONEVENTCUT_H_
#define PWGEM_PHOTONMESON_CORE_EMPHOTONEVENTCUT_H_

#include "PWGEM/PhotonMeson/Utils/EventUtilities.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/TriggerAliases.h"

#include <Framework/Concepts.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/Logger.h>

#include <TH1.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <format>
#include <limits>
#include <string>
#include <string_view>
#include <vector>

class EMPhotonEventCut
{
 public:
  EMPhotonEventCut() = default;
  EMPhotonEventCut(const char* name, const char* title) : name(name), title(title) {}

  enum class EMPhotonEventCuts : int {
    kSel8 = 0,
    kFT0AND,
    kZvtx,
    kNoTFB,
    kNoITSROFB,
    kNoSameBunchPileup,
    kIsVertexITSTPC,
    kIsVertexTOFmatched,
    kIsGoodZvtxFT0vsPV,
    kNoCollInTimeRangeStandard,
    kNoCollInTimeRangeStrict,
    kNoCollInITSROFStandard,
    kNoCollInITSROFStrict,
    kNoHighMultCollInPrevRof,
    kIsGoodITSLayer3,
    kIsGoodITSLayer0123,
    kIsGoodITSLayersAll,
    kEMCReadoutInMB,
    kEMCHardwareTriggered,
    kTrackOccupancy,
    kFT0COccupancy,
    kGoodQvec,
    kNCuts
  };

  /// bins of the Event/{before,after}/hCollisionCounter histogram
  enum class EventCounterBin : uint8_t {
    kAll = 1,
    kNoTFBorder = 2,
    kNoITSROFBorder = 3,
    kNoSameBunchPileup = 4,
    kVertexITSTPC = 5,
    kGoodZVtxFT0vsPV = 6,
    kFT0AND = 7,
    kSel8 = 8,
    kZVtx = 9,
    kTrackOccupancy = 10,
    kFT0COccupancy = 11,
    kEMCMBReadout = 12,
    kEMCL0Trigger = 13,
    kGoodQvec = 14,
    kAccepted = 15
  };

  [[nodiscard]] const std::string& getName() const { return name; }
  [[nodiscard]] const std::string& getTitle() const { return title; }

  template <typename T>
  bool IsSelected(T const& collision) const
  {
    if (mRequireSel8 && !IsSelected(collision, EMPhotonEventCuts::kSel8)) {
      return false;
    }
    if (mRequireFT0AND && !IsSelected(collision, EMPhotonEventCuts::kFT0AND)) {
      return false;
    }
    if (!IsSelected(collision, EMPhotonEventCuts::kZvtx)) {
      return false;
    }
    if (mRequireNoTFB && !IsSelected(collision, EMPhotonEventCuts::kNoTFB)) {
      return false;
    }
    if (mRequireNoITSROFB && !IsSelected(collision, EMPhotonEventCuts::kNoITSROFB)) {
      return false;
    }
    if (mRequireNoSameBunchPileup && !IsSelected(collision, EMPhotonEventCuts::kNoSameBunchPileup)) {
      return false;
    }
    if (mRequireVertexITSTPC && !IsSelected(collision, EMPhotonEventCuts::kIsVertexITSTPC)) {
      return false;
    }
    if (mRequireVertexTOFmatched && !IsSelected(collision, EMPhotonEventCuts::kIsVertexTOFmatched)) {
      return false;
    }
    if (mRequireGoodZvtxFT0vsPV && !IsSelected(collision, EMPhotonEventCuts::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    if (mRequireNoCollInTimeRangeStandard && !IsSelected(collision, EMPhotonEventCuts::kNoCollInTimeRangeStandard)) {
      return false;
    }
    if (mRequireNoCollInTimeRangeStrict && !IsSelected(collision, EMPhotonEventCuts::kNoCollInTimeRangeStrict)) {
      return false;
    }
    if (mRequireNoCollInITSROFStandard && !IsSelected(collision, EMPhotonEventCuts::kNoCollInITSROFStandard)) {
      return false;
    }
    if (mRequireNoCollInITSROFStrict && !IsSelected(collision, EMPhotonEventCuts::kNoCollInITSROFStrict)) {
      return false;
    }
    if (mRequireNoHighMultCollInPrevRof && !IsSelected(collision, EMPhotonEventCuts::kNoHighMultCollInPrevRof)) {
      return false;
    }
    if (mRequireGoodITSLayer3 && !IsSelected(collision, EMPhotonEventCuts::kIsGoodITSLayer3)) {
      return false;
    }
    if (mRequireGoodITSLayer0123 && !IsSelected(collision, EMPhotonEventCuts::kIsGoodITSLayer0123)) {
      return false;
    }
    if (mRequireGoodITSLayersAll && !IsSelected(collision, EMPhotonEventCuts::kIsGoodITSLayersAll)) {
      return false;
    }
    if (mRequireEMCReadoutInMB && !IsSelected(collision, EMPhotonEventCuts::kEMCReadoutInMB)) {
      return false;
    }
    if (mRequireEMCHardwareTriggered && !IsSelected(collision, EMPhotonEventCuts::kEMCHardwareTriggered)) {
      return false;
    }
    if (!IsSelected(collision, EMPhotonEventCuts::kTrackOccupancy)) {
      return false;
    }
    if (!IsSelected(collision, EMPhotonEventCuts::kFT0COccupancy)) {
      return false;
    }
    if (mRequireGoodQvec && !IsSelected(collision, EMPhotonEventCuts::kGoodQvec)) {
      return false;
    }
    return true;
  }

  template <typename T>
  bool IsSelected(T const& collision, const EMPhotonEventCuts& cut) const
  {
    switch (cut) {
      case EMPhotonEventCuts::kSel8:
        return collision.sel8();

      case EMPhotonEventCuts::kFT0AND:
        return collision.selection_bit(o2::aod::evsel::kIsTriggerTVX);

      case EMPhotonEventCuts::kZvtx:
        return mMinZvtx < collision.posZ() && collision.posZ() < mMaxZvtx;

      case EMPhotonEventCuts::kNoTFB:
        return collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder);

      case EMPhotonEventCuts::kNoITSROFB:
        return collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder);

      case EMPhotonEventCuts::kNoSameBunchPileup:
        return collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup);

      case EMPhotonEventCuts::kIsVertexITSTPC:
        return collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC);

      case EMPhotonEventCuts::kIsVertexTOFmatched:
        return collision.selection_bit(o2::aod::evsel::kIsVertexTOFmatched);

      case EMPhotonEventCuts::kIsGoodZvtxFT0vsPV:
        return collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV);

      case EMPhotonEventCuts::kNoCollInTimeRangeStandard:
        return collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard);

      case EMPhotonEventCuts::kNoCollInTimeRangeStrict:
        return collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStrict);

      case EMPhotonEventCuts::kNoCollInITSROFStandard:
        return collision.selection_bit(o2::aod::evsel::kNoCollInRofStandard);

      case EMPhotonEventCuts::kNoCollInITSROFStrict:
        return collision.selection_bit(o2::aod::evsel::kNoCollInRofStrict);

      case EMPhotonEventCuts::kNoHighMultCollInPrevRof:
        return collision.selection_bit(o2::aod::evsel::kNoHighMultCollInPrevRof);

      case EMPhotonEventCuts::kIsGoodITSLayer3:
        return collision.selection_bit(o2::aod::evsel::kIsGoodITSLayer3);

      case EMPhotonEventCuts::kIsGoodITSLayer0123:
        return collision.selection_bit(o2::aod::evsel::kIsGoodITSLayer0123);

      case EMPhotonEventCuts::kIsGoodITSLayersAll:
        return collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll);

      case EMPhotonEventCuts::kEMCReadoutInMB:
        return (collision.alias_bit(kTVXinEMC));

      case EMPhotonEventCuts::kEMCHardwareTriggered:
        return (collision.alias_bit(kEMC7) || collision.alias_bit(kDMC7));

      case EMPhotonEventCuts::kTrackOccupancy:
        if constexpr (o2::aod::pwgem::pm::utils::event::HasTrackOccupancy<T>) {
          return mMinTrackOcc <= collision.trackOccupancyInTimeRange() && collision.trackOccupancyInTimeRange() < mMaxTrackOcc;
        } else {
          // the event table has no track occupancy, the cut can not be applied
          warnOccupancyNotAvailable(isTrackOccupancyCutSet(), "track");
          return true;
        }

      case EMPhotonEventCuts::kFT0COccupancy:
        if constexpr (o2::aod::pwgem::pm::utils::event::HasFT0COccupancy<T>) {
          return mMinFT0COcc <= collision.ft0cOccupancyInTimeRange() && collision.ft0cOccupancyInTimeRange() < mMaxFT0COcc;
        } else {
          // the event table has no FT0C occupancy, the cut can not be applied
          warnOccupancyNotAvailable(isFT0COccupancyCutSet(), "FT0C");
          return true;
        }

      case EMPhotonEventCuts::kGoodQvec:
        if constexpr (o2::aod::pwgem::pm::utils::event::HasQvec<T>) {
          return o2::aod::pwgem::pm::utils::event::isQvecGood(collision, mQvecDetectors, mQvecHarmonic, mMaxQvec);
        } else {
          // the event table has no Q vectors, the cut can not be applied
          if (mRequireGoodQvec && !mWarnedNoQvec) {
            LOG(warning) << "EM Photon Event Cut: Q vector cut requested, but the collision table has no Q vectors. Q vector cut is not applied!";
            mWarnedNoQvec = true;
          }
          return true;
        }

      default:
        return true;
    }
  }

  /// \brief add the event QA histograms to the given registry. Call this after all cut values are set!
  /// \param fRegistry pointer to the histogram registry of the task
  /// \param useWeight whether the histograms should store weights (Sumw2)
  void addEventHistograms(o2::framework::HistogramRegistry* fRegistry, bool useWeight = false) const
  {
    auto hCollisionCounter = fRegistry->add<TH1>("Event/before/hCollisionCounter", "collision counter;;Number of events", o2::framework::HistType::kTH1D, {{static_cast<int>(EventCounterBin::kAccepted), 0.5, static_cast<float>(EventCounterBin::kAccepted) + 0.5f}}, useWeight);
    hCollisionCounter->GetXaxis()->SetBinLabel(static_cast<int>(EventCounterBin::kAll), "all");
    hCollisionCounter->GetXaxis()->SetBinLabel(static_cast<int>(EventCounterBin::kNoTFBorder), "No TF border");
    hCollisionCounter->GetXaxis()->SetBinLabel(static_cast<int>(EventCounterBin::kNoITSROFBorder), "No ITS ROF border");
    hCollisionCounter->GetXaxis()->SetBinLabel(static_cast<int>(EventCounterBin::kNoSameBunchPileup), "No Same Bunch Pileup");
    hCollisionCounter->GetXaxis()->SetBinLabel(static_cast<int>(EventCounterBin::kVertexITSTPC), "Is Vertex ITSTPC");
    hCollisionCounter->GetXaxis()->SetBinLabel(static_cast<int>(EventCounterBin::kGoodZVtxFT0vsPV), "Is Good Zvtx FT0vsPV");
    hCollisionCounter->GetXaxis()->SetBinLabel(static_cast<int>(EventCounterBin::kFT0AND), "FT0AND");
    hCollisionCounter->GetXaxis()->SetBinLabel(static_cast<int>(EventCounterBin::kSel8), "sel8");
    hCollisionCounter->GetXaxis()->SetBinLabel(static_cast<int>(EventCounterBin::kZVtx), std::format("{:g} < Z_{{vtx}} < {:g} cm", mMinZvtx, mMaxZvtx).c_str());
    hCollisionCounter->GetXaxis()->SetBinLabel(static_cast<int>(EventCounterBin::kTrackOccupancy), "track occupancy");
    hCollisionCounter->GetXaxis()->SetBinLabel(static_cast<int>(EventCounterBin::kFT0COccupancy), "FT0C occupancy");
    hCollisionCounter->GetXaxis()->SetBinLabel(static_cast<int>(EventCounterBin::kEMCMBReadout), "EMC MB Readout");
    hCollisionCounter->GetXaxis()->SetBinLabel(static_cast<int>(EventCounterBin::kEMCL0Trigger), "EMC L0 Triggered");
    hCollisionCounter->GetXaxis()->SetBinLabel(static_cast<int>(EventCounterBin::kGoodQvec), "good Q vector");
    hCollisionCounter->GetXaxis()->SetBinLabel(static_cast<int>(EventCounterBin::kAccepted), "accepted");

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

  /// \brief fill the event QA histograms
  /// \tparam ev_id 0 for before the event selection, 1 for after the event selection
  /// \param fRegistry pointer to the histogram registry of the task
  /// \param collision collision to fill
  /// \param weight weight of the collision
  template <const int ev_id, o2::soa::is_iterator TCollision>
  void fillEventInfo(o2::framework::HistogramRegistry* fRegistry, TCollision const& collision, float weight = 1.f) const
  {
    static constexpr std::array<std::string_view, 2> EventTypes = {"before/", "after/"};
    auto fillCounter = [&](EventCounterBin bin) {
      fRegistry->fill(HIST("Event/") + HIST(EventTypes[ev_id]) + HIST("hCollisionCounter"), static_cast<float>(bin), weight);
    };

    fillCounter(EventCounterBin::kAll);
    if (IsSelected(collision, EMPhotonEventCuts::kNoTFB)) {
      fillCounter(EventCounterBin::kNoTFBorder);
    }
    if (IsSelected(collision, EMPhotonEventCuts::kNoITSROFB)) {
      fillCounter(EventCounterBin::kNoITSROFBorder);
    }
    if (IsSelected(collision, EMPhotonEventCuts::kNoSameBunchPileup)) {
      fillCounter(EventCounterBin::kNoSameBunchPileup);
    }
    if (IsSelected(collision, EMPhotonEventCuts::kIsVertexITSTPC)) {
      fillCounter(EventCounterBin::kVertexITSTPC);
    }
    if (IsSelected(collision, EMPhotonEventCuts::kIsGoodZvtxFT0vsPV)) {
      fillCounter(EventCounterBin::kGoodZVtxFT0vsPV);
    }
    if (IsSelected(collision, EMPhotonEventCuts::kFT0AND)) {
      fillCounter(EventCounterBin::kFT0AND);
    }
    if (IsSelected(collision, EMPhotonEventCuts::kSel8)) {
      fillCounter(EventCounterBin::kSel8);
    }
    if (IsSelected(collision, EMPhotonEventCuts::kZvtx)) {
      fillCounter(EventCounterBin::kZVtx);
    }
    if (IsSelected(collision, EMPhotonEventCuts::kTrackOccupancy)) {
      fillCounter(EventCounterBin::kTrackOccupancy);
    }
    if (IsSelected(collision, EMPhotonEventCuts::kFT0COccupancy)) {
      fillCounter(EventCounterBin::kFT0COccupancy);
    }
    if (IsSelected(collision, EMPhotonEventCuts::kEMCReadoutInMB)) {
      fillCounter(EventCounterBin::kEMCMBReadout);
    }
    if (IsSelected(collision, EMPhotonEventCuts::kEMCHardwareTriggered)) {
      fillCounter(EventCounterBin::kEMCL0Trigger);
    }
    if constexpr (o2::aod::pwgem::pm::utils::event::HasQvec<TCollision>) {
      if (IsSelected(collision, EMPhotonEventCuts::kGoodQvec)) {
        fillCounter(EventCounterBin::kGoodQvec);
      }
    }

    fRegistry->fill(HIST("Event/") + HIST(EventTypes[ev_id]) + HIST("hZvtx"), collision.posZ(), weight);
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

  /// \brief fill the "accepted" bin of the before and after collision counter
  /// \param fRegistry pointer to the histogram registry of the task
  /// \param weight weight of the collision
  void fillAccepted(o2::framework::HistogramRegistry* fRegistry, float weight = 1.f) const
  {
    fRegistry->fill(HIST("Event/before/hCollisionCounter"), static_cast<float>(EventCounterBin::kAccepted), weight);
    fRegistry->fill(HIST("Event/after/hCollisionCounter"), static_cast<float>(EventCounterBin::kAccepted), weight);
  }

  // Setters
  void SetRequireSel8(bool flag);
  void SetRequireFT0AND(bool flag);
  void SetZvtxRange(float min, float max);
  void SetRequireNoTFB(bool flag);
  void SetRequireNoITSROFB(bool flag);
  void SetRequireNoSameBunchPileup(bool flag);
  void SetRequireVertexITSTPC(bool flag);
  void SetRequireVertexTOFmatched(bool flag);
  void SetRequireGoodZvtxFT0vsPV(bool flag);
  void SetRequireNoCollInTimeRangeStandard(bool flag);
  void SetRequireNoCollInTimeRangeStrict(bool flag);
  void SetRequireNoCollInITSROFStandard(bool flag);
  void SetRequireNoCollInITSROFStrict(bool flag);
  void SetRequireNoHighMultCollInPrevRof(bool flag);
  void SetRequireGoodITSLayer3(bool flag);
  void SetRequireGoodITSLayer0123(bool flag);
  void SetRequireGoodITSLayersAll(bool flag);
  void SetRequireEMCReadoutInMB(bool flag);
  void SetRequireEMCHardwareTriggered(bool flag);
  void SetTrackOccupancyRange(float min, float max);
  void SetFT0COccupancyRange(float min, float max);
  void SetRequireGoodQvec(bool flag);
  void SetQvecHarmonic(int harmonic);
  void SetMaxQvec(float max);
  void SetQvecDetectors(std::vector<int> const& detectors);

 private:
  [[nodiscard]] bool isTrackOccupancyCutSet() const
  {
    return mMinTrackOcc != std::numeric_limits<float>::lowest() || mMaxTrackOcc != std::numeric_limits<float>::max();
  }

  [[nodiscard]] bool isFT0COccupancyCutSet() const
  {
    return mMinFT0COcc != std::numeric_limits<float>::lowest() || mMaxFT0COcc != std::numeric_limits<float>::max();
  }

  void warnOccupancyNotAvailable(bool cutSet, const char* estimator) const
  {
    if (cutSet && !mWarnedNoOccupancy) {
      LOG(warning) << "EM Photon Event Cut: " << estimator << " occupancy cut set, but the collision table has no " << estimator << " occupancy. Occupancy cut is not applied!";
      mWarnedNoOccupancy = true;
    }
  }

  std::string name;
  std::string title;
  bool mRequireSel8{false};
  bool mRequireFT0AND{true};
  float mMinZvtx{-10.f}, mMaxZvtx{+10.f};
  bool mRequireNoTFB{false};
  bool mRequireNoITSROFB{false};
  bool mRequireNoSameBunchPileup{false};
  bool mRequireVertexITSTPC{false};
  bool mRequireVertexTOFmatched{false};
  bool mRequireGoodZvtxFT0vsPV{false};
  bool mRequireNoCollInTimeRangeStandard{false};
  bool mRequireNoCollInTimeRangeStrict{false};
  bool mRequireNoCollInITSROFStandard{false};
  bool mRequireNoCollInITSROFStrict{false};
  bool mRequireNoHighMultCollInPrevRof{false};
  bool mRequireGoodITSLayer3{false};
  bool mRequireGoodITSLayer0123{false};
  bool mRequireGoodITSLayersAll{false};
  bool mRequireEMCReadoutInMB{false};
  bool mRequireEMCHardwareTriggered{false};

  // occupancy cuts, by default no cut is applied
  float mMinTrackOcc{std::numeric_limits<float>::lowest()}, mMaxTrackOcc{std::numeric_limits<float>::max()};
  float mMinFT0COcc{std::numeric_limits<float>::lowest()}, mMaxFT0COcc{std::numeric_limits<float>::max()};

  // Q vector cut, off by default and only applied if the collision table has Q vectors
  bool mRequireGoodQvec{false};
  int mQvecHarmonic{o2::aod::pwgem::pm::utils::event::kElliptic};
  float mMaxQvec{20.f};
  std::vector<int> mQvecDetectors;
  mutable bool mWarnedNoQvec{false};
  mutable bool mWarnedNoOccupancy{false};
};

#endif // PWGEM_PHOTONMESON_CORE_EMPHOTONEVENTCUT_H_
