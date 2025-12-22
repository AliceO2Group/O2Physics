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

/// \file utilsEvSelHf.h
/// \brief Event selection utilities for HF analyses
/// \author Mattia Faggin <mfaggin@cern.ch>, CERN
/// \author Vít Kučera <vit.kucera@cern.ch>, Inha University
/// \author Luigi Dello Stritto <luigi.dello.stritto@cern.ch>, CERN

#ifndef PWGHF_UTILS_UTILSEVSELHF_H_
#define PWGHF_UTILS_UTILSEVSELHF_H_

#include "PWGHF/Core/CentralityEstimation.h"
//
#include "PWGUD/Core/SGCutParHolder.h"
#include "PWGUD/Core/SGSelector.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/RCTSelectionFlags.h"
#include "Common/CCDB/ctpRateFetcher.h"
#include "Common/Core/CollisionTypeHelper.h"
#include "Common/Core/Zorro.h"
#include "Common/Core/ZorroSummary.h"

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPLHCIFData.h"
#include <Framework/AnalysisHelpers.h>
#include <Framework/Configurable.h>
#include <Framework/DeviceSpec.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/Logger.h>

#include <TH1.h>
#include <TH2.h>

#include <Rtypes.h>

#include <cstddef>
#include <cstdint>
#include <memory> // std::shared_ptr
#include <string> // std::string

namespace o2::hf_occupancy
{
// centrality selection estimators
enum OccupancyEstimator { None = 0,
                          Its,
                          Ft0c };

/// Get the occupancy
/// \param collision is the collision with the occupancy information
/// \return collision occupancy
template <typename TCollision>
float getOccupancyColl(TCollision const& collision, const int occEstimator)
{
  switch (occEstimator) {
    case OccupancyEstimator::Its:
      return static_cast<float>(collision.trackOccupancyInTimeRange());
    case OccupancyEstimator::Ft0c:
      return static_cast<float>(collision.ft0cOccupancyInTimeRange());
    default:
      LOG(fatal) << "Occupancy estimator not valid. See OccupancyEstimator for valid values.";
      break;
  }
  return -999.f;
};

/// \brief Function to get MC collision occupancy
/// \param collSlice collection of reconstructed collisions associated to a generated one
/// \return generated MC collision occupancy
template <typename TCollisions>
float getOccupancyGenColl(TCollisions const& collSlice, const int occEstimator)
{
  using TMult = uint16_t; // type of numContrib
  TMult multiplicity{};
  float occupancy{0.f};
  for (const auto& collision : collSlice) {
    const TMult collMult = collision.numContrib();
    if (collMult > multiplicity) {
      occupancy = getOccupancyColl(collision, occEstimator);
      multiplicity = collMult;
    }
  } // end loop over collisions
  return occupancy;
};
} // namespace o2::hf_occupancy

namespace o2::hf_evsel
{
// event rejection types
enum EventRejection {
  None = 0,
  Rct,
  SoftwareTrigger,
  Centrality,
  Trigger,
  TvxTrigger,
  TimeFrameBorderCut,
  ItsRofBorderCut,
  IsGoodZvtxFT0vsPV,
  NoSameBunchPileup,
  Occupancy,
  NContrib,
  NoCollInTimeRangeNarrow,
  NoCollInTimeRangeStandard,
  NoCollInRofStandard,
  UpcEventCut,
  Chi2,
  PositionZ,
  NEventRejection
};

using HfCollisionRejectionMask = uint32_t; // 32 bits, in case new ev. selections will be added

const o2::framework::AxisSpec axisEvents = {EventRejection::NEventRejection, -0.5f, +EventRejection::NEventRejection - 0.5f, ""};
const o2::framework::AxisSpec axisUpcEvents = {o2::aod::sgselector::DoubleGap + 1, -0.5f, +o2::aod::sgselector::DoubleGap + 0.5f, ""};

/// \brief Function to put labels on monitoring histogram
/// \param hRejection monitoring histogram
/// \param softwareTriggerLabel bin label for software trigger rejection
template <typename Histo>
void setEventRejectionLabels(Histo& hRejection, std::string const& softwareTriggerLabel = "")
{
  // Puts labels on the collision monitoring histogram.
  hRejection->GetXaxis()->SetBinLabel(EventRejection::None + 1, "All");
  hRejection->GetXaxis()->SetBinLabel(EventRejection::Rct + 1, "RCT");
  hRejection->GetXaxis()->SetBinLabel(EventRejection::SoftwareTrigger + 1, softwareTriggerLabel.data());
  hRejection->GetXaxis()->SetBinLabel(EventRejection::Centrality + 1, "Centrality");
  hRejection->GetXaxis()->SetBinLabel(EventRejection::Trigger + 1, "Trigger");
  hRejection->GetXaxis()->SetBinLabel(EventRejection::TvxTrigger + 1, "TVX Trigger");
  hRejection->GetXaxis()->SetBinLabel(EventRejection::TimeFrameBorderCut + 1, "TF border");
  hRejection->GetXaxis()->SetBinLabel(EventRejection::ItsRofBorderCut + 1, "ITS ROF border");
  hRejection->GetXaxis()->SetBinLabel(EventRejection::UpcEventCut + 1, "UPC event");
  hRejection->GetXaxis()->SetBinLabel(EventRejection::IsGoodZvtxFT0vsPV + 1, "PV #it{z} consistency FT0 timing");
  hRejection->GetXaxis()->SetBinLabel(EventRejection::NoSameBunchPileup + 1, "No same-bunch pile-up"); // POTENTIALLY BAD FOR BEAUTY ANALYSES
  hRejection->GetXaxis()->SetBinLabel(EventRejection::Occupancy + 1, "Occupancy");
  hRejection->GetXaxis()->SetBinLabel(EventRejection::NContrib + 1, "# of PV contributors");
  hRejection->GetXaxis()->SetBinLabel(EventRejection::Chi2 + 1, "PV #it{#chi}^{2}");
  hRejection->GetXaxis()->SetBinLabel(EventRejection::PositionZ + 1, "PV #it{z}");
  hRejection->GetXaxis()->SetBinLabel(EventRejection::NoCollInTimeRangeNarrow + 1, "No coll timerange narrow");
  hRejection->GetXaxis()->SetBinLabel(EventRejection::NoCollInTimeRangeStandard + 1, "No coll timerange strict");
  hRejection->GetXaxis()->SetBinLabel(EventRejection::NoCollInRofStandard + 1, "No coll in ROF std");
}

struct HfEventSelection : o2::framework::ConfigurableGroup {
  std::string prefix = "hfEvSel"; // JSON group name
  // event selection parameters (in chronological order of application)
  o2::framework::Configurable<float> centralityMin{"centralityMin", -10.f, "Minimum centrality (0 rejects gen. collisions with no reco. collision)"};
  o2::framework::Configurable<float> centralityMax{"centralityMax", 100.f, "Maximum centrality"};
  o2::framework::Configurable<bool> useSel8Trigger{"useSel8Trigger", true, "Apply the sel8 event selection"};
  o2::framework::Configurable<int> triggerClass{"triggerClass", -1, "Trigger class different from sel8 (e.g. kINT7 for Run2) used only if useSel8Trigger is false"};
  o2::framework::Configurable<bool> useTvxTrigger{"useTvxTrigger", true, "Apply TVX trigger sel"};
  o2::framework::Configurable<bool> useTimeFrameBorderCut{"useTimeFrameBorderCut", true, "Apply TF border cut"};
  o2::framework::Configurable<bool> useItsRofBorderCut{"useItsRofBorderCut", true, "Apply ITS ROF border cut"};
  o2::framework::Configurable<bool> useIsGoodZvtxFT0vsPV{"useIsGoodZvtxFT0vsPV", false, "Check consistency between PVz from central barrel with that from FT0 timing"};
  o2::framework::Configurable<bool> useNoSameBunchPileup{"useNoSameBunchPileup", false, "Exclude collisions in bunches with more than 1 reco. PV"}; // POTENTIALLY BAD FOR BEAUTY ANALYSES
  o2::framework::Configurable<bool> useOccupancyCut{"useOccupancyCut", false, "Apply occupancy selection (num. ITS tracks with at least 5 clusters or num. of signals in FT0c in +-100us from current collision)"};
  o2::framework::Configurable<int> occEstimator{"occEstimator", 1, "Occupancy estimation (1: ITS, 2: FT0C)"};
  o2::framework::Configurable<int> occupancyMin{"occupancyMin", 0, "Minimum occupancy"};
  o2::framework::Configurable<int> occupancyMax{"occupancyMax", 1000000, "Maximum occupancy"};
  o2::framework::Configurable<int> nPvContributorsMin{"nPvContributorsMin", 0, "Minimum number of PV contributors"};
  o2::framework::Configurable<float> chi2PvMax{"chi2PvMax", -1.f, "Maximum PV chi2"};
  o2::framework::Configurable<float> zPvPosMin{"zPvPosMin", -10.f, "Minimum PV posZ (cm)"};
  o2::framework::Configurable<float> zPvPosMax{"zPvPosMax", 10.f, "Maximum PV posZ (cm)"};
  o2::framework::Configurable<bool> useNoCollInTimeRangeNarrow{"useNoCollInTimeRangeNarrow", false, "Reject collisions in time range narrow"};
  o2::framework::Configurable<bool> useNoCollInTimeRangeStandard{"useNoCollInTimeRangeStandard", false, "Reject collisions in time range strict"};
  o2::framework::Configurable<bool> useNoCollInRofStandard{"useNoCollInRofStandard", false, "Reject collisions in ROF standard"};
  o2::framework::Configurable<std::string> softwareTrigger{"softwareTrigger", "", "Label of software trigger. Multiple triggers can be selected dividing them by a comma. Set None if you want bcs that are not selected by any trigger"};
  o2::framework::Configurable<uint64_t> bcMarginForSoftwareTrigger{"bcMarginForSoftwareTrigger", 100, "Number of BCs of margin for software triggers"};
  o2::framework::Configurable<std::string> ccdbPathSoftwareTrigger{"ccdbPathSoftwareTrigger", "EventFiltering/Zorro/", "ccdb path for ZORRO objects"};
  o2::framework::ConfigurableAxis th2ConfigAxisCent{"th2ConfigAxisCent", {100, 0., 100.}, ""};
  o2::framework::ConfigurableAxis th2ConfigAxisOccupancy{"th2ConfigAxisOccupancy", {100, 0, 100000}, ""};
  o2::framework::ConfigurableAxis th2ConfigAxisInteractionRate{"th2ConfigAxisInteractionRate", {500, 0, 50000}, ""};
  o2::framework::Configurable<bool> requireGoodRct{"requireGoodRct", false, "Flag to require good RCT"};
  o2::framework::Configurable<std::string> rctLabel{"rctLabel", "CBT_hadronPID", "RCT selection flag (CBT, CBT_hadronPID, CBT_electronPID, CBT_calo, CBT_muon, CBT_muon_glo)"};
  o2::framework::Configurable<bool> rctCheckZDC{"rctCheckZDC", false, "RCT flag to check whether the ZDC is present or not"};
  o2::framework::Configurable<bool> rctTreatLimitedAcceptanceAsBad{"rctTreatLimitedAcceptanceAsBad", false, "RCT flag to reject events with limited acceptance for selected detectors"};
  o2::framework::Configurable<std::string> irSource{"irSource", "", "Estimator of the interaction rate (Empty: automatically set. Otherwise recommended: pp --> T0VTX, Pb-Pb --> ZNC hadronic)"};

  //  SG selector
  SGSelector sgSelector;

  // histogram names
  static constexpr char NameHistCollisions[] = "hCollisions";
  static constexpr char NameHistSelCollisionsCent[] = "hSelCollisionsCent";
  static constexpr char NameHistPosZBeforeEvSel[] = "hPosZBeforeEvSel";
  static constexpr char NameHistPosZAfterEvSel[] = "hPosZAfterEvSel";
  static constexpr char NameHistPosXAfterEvSel[] = "hPosXAfterEvSel";
  static constexpr char NameHistPosYAfterEvSel[] = "hPosYAfterEvSel";
  static constexpr char NameHistNumPvContributorsAfterSel[] = "hNumPvContributorsAfterSel";
  static constexpr char NameHistCollisionsCentOcc[] = "hCollisionsCentOcc";
  static constexpr char NameHistCollisionsCentIR[] = "hCollisionsCentIR";
  static constexpr char NameHistUpCollisions[] = "hUpCollisions";

  std::shared_ptr<TH1> hCollisions, hSelCollisionsCent, hPosZBeforeEvSel, hPosZAfterEvSel, hPosXAfterEvSel, hPosYAfterEvSel, hNumPvContributorsAfterSel, hUpCollisions;
  std::shared_ptr<TH2> hCollisionsCentOcc;
  std::shared_ptr<TH2> hCollisionsCentIR;

  // util to retrieve the RCT info from CCDB
  o2::aod::rctsel::RCTFlagsChecker rctChecker;

  // util to retrieve trigger mask in case of software triggers
  Zorro zorro;
  int currentRun{-1};

  // util to retrieve IR
  ctpRateFetcher irFetcher;
  std::string irSourceForCptFetcher;

  // guard variable to guarantee full configuration
  // important for RCT, Zorro
  bool isInitCalled{false};

  // guard variable to guarantee that histograms are added to the registry before filling them
  bool areHistosInRegistry{false};

  /// Set standard preselection gap trigger (values taken from UD group)
  SGCutParHolder setSgPreselection()
  {
    SGCutParHolder sgCuts;
    sgCuts.SetNDtcoll(1);       // Minimum number of sigma around the collision
    sgCuts.SetMinNBCs(2);       // Minimum number of bunch crossings
    sgCuts.SetNTracks(2, 1000); // Minimum and maximum number of PV contributors
    sgCuts.SetMaxFITtime(34.f); // Maximum FIT time in ns

    // Set FIT amplitudes: FV0, FT0A, FT0C, FDDA, FDDC
    sgCuts.SetFITAmpLimits({-1.f, 1000.f, 1000.f, -1.f, -1.f});

    return sgCuts;
  }

  /// \brief Adds collision monitoring histograms in the histogram registry.
  /// \param registry reference to the histogram registry
  void addHistograms(o2::framework::HistogramRegistry& registry)
  {
    hCollisions = registry.add<TH1>(NameHistCollisions, "HF event counter;;# of accepted collisions", {o2::framework::HistType::kTH1D, {axisEvents}});
    hSelCollisionsCent = registry.add<TH1>(NameHistSelCollisionsCent, "HF event counter;T0M;# of accepted collisions", {o2::framework::HistType::kTH1D, {{100, 0., 100.}}});
    hPosZBeforeEvSel = registry.add<TH1>(NameHistPosZBeforeEvSel, "all events;#it{z}_{prim. vtx.} (cm);entries", {o2::framework::HistType::kTH1D, {{400, -20., 20.}}});
    hPosZAfterEvSel = registry.add<TH1>(NameHistPosZAfterEvSel, "selected events;#it{z}_{prim. vtx.} (cm);entries", {o2::framework::HistType::kTH1D, {{400, -20., 20.}}});
    hPosXAfterEvSel = registry.add<TH1>(NameHistPosXAfterEvSel, "selected events;#it{x}_{prim. vtx.} (cm);entries", {o2::framework::HistType::kTH1D, {{200, -0.5, 0.5}}});
    hPosYAfterEvSel = registry.add<TH1>(NameHistPosYAfterEvSel, "selected events;#it{y}_{prim. vtx.} (cm);entries", {o2::framework::HistType::kTH1D, {{200, -0.5, 0.5}}});
    hNumPvContributorsAfterSel = registry.add<TH1>(NameHistNumPvContributorsAfterSel, "selected events;number of prim. vtx. contributors;entries", {o2::framework::HistType::kTH1D, {{500, -0.5, 499.5}}});
    setEventRejectionLabels(hCollisions, softwareTrigger);
    hUpCollisions = registry.add<TH1>(NameHistUpCollisions, "HF UPC counter;;# of UPC events", {o2::framework::HistType::kTH1D, {axisUpcEvents}});
    const o2::framework::AxisSpec th2AxisCent{th2ConfigAxisCent, "Centrality"};
    const o2::framework::AxisSpec th2AxisOccupancy{th2ConfigAxisOccupancy, "Occupancy"};
    const o2::framework::AxisSpec th2AxisInteractionRate{th2ConfigAxisInteractionRate, "Interaction Rate [Hz]"};

    hCollisionsCentOcc = registry.add<TH2>(NameHistCollisionsCentOcc, "selected events;Centrality; Occupancy", {o2::framework::HistType::kTH2D, {th2AxisCent, th2AxisOccupancy}});
    hCollisionsCentIR = registry.add<TH2>(NameHistCollisionsCentIR, "selected events;Centrality; Interaction Rate [Hz]", {o2::framework::HistType::kTH2D, {th2AxisCent, th2AxisInteractionRate}});

    // histograms in registry
    // let's update the guard variable
    areHistosInRegistry = true;
  }

  /// \brief Inits the HF event selection object
  /// \param registry reference to the histogram registry
  void init(o2::framework::HistogramRegistry& registry, o2::framework::OutputObj<ZorroSummary>* zorroSummary = nullptr)
  {
    // we initialise the RCT checker
    if (requireGoodRct) {
      rctChecker.init(rctLabel.value, rctCheckZDC.value, rctTreatLimitedAcceptanceAsBad.value);
    }

    // we initialise the summary object
    if (!softwareTrigger.value.empty()) {
      if (zorroSummary == nullptr) {
        LOGP(fatal, "No OutputObj<ZorroSummary> provided to HF event selection object in your task. Add it if you want to get the normalisation from Zorro.");
        return;
      }
      zorroSummary->setObject(zorro.getZorroSummary());
    }

    // we initialise histograms
    addHistograms(registry);

    // we initialise IR fetcher
    if (!irSource.value.empty()) {
      irSourceForCptFetcher = irSource.value;
    }

    // full configuration complete: update the guard variable
    isInitCalled = true;
  }

  /// \brief Applies event selection.
  /// \tparam useEvSel use information from the EvSel table
  /// \tparam centEstimator centrality estimator
  /// \param collision collision to test against the selection criteria
  /// \param centrality collision centrality variable to be set in this function
  /// \param ccdb ccdb service needed to retrieve the needed info for zorro
  /// \param registry reference to the histogram registry needed for zorro
  /// \return bitmask with the event selection criteria not satisfied by the collision
  template <bool UseEvSel, o2::hf_centrality::CentralityEstimator CentEstimator, typename TBcs, typename TCollision>
  HfCollisionRejectionMask getHfCollisionRejectionMask(TCollision const& collision,
                                                       float& centrality,
                                                       o2::framework::Service<o2::ccdb::BasicCCDBManager> const& ccdb,
                                                       o2::framework::HistogramRegistry& registry)
  {
    HfCollisionRejectionMask rejectionMask{};

    if constexpr (CentEstimator != o2::hf_centrality::CentralityEstimator::None) {
      centrality = o2::hf_centrality::getCentralityColl(collision, CentEstimator);
      if (centrality < centralityMin || centrality > centralityMax) {
        SETBIT(rejectionMask, EventRejection::Centrality);
      }
    }

    if constexpr (UseEvSel) {
      /// RCT condition
      if (requireGoodRct) {
        if (!isInitCalled) {
          // protect against incomplete configuration
          LOG(fatal) << "Checking RCT flags w/o full HF event-selection configuration. Call the function HfEventSelection::init() to fix.";
        }
        if (!rctChecker.checkTable(collision)) {
          SETBIT(rejectionMask, EventRejection::Rct);
        }
      }
      /// trigger condition
      if ((useSel8Trigger && !collision.sel8()) || (!useSel8Trigger && triggerClass > -1 && !collision.alias_bit(triggerClass))) {
        SETBIT(rejectionMask, EventRejection::Trigger);
      }
      /// TVX trigger selection
      if (useTvxTrigger && !collision.selection_bit(o2::aod::evsel::kIsTriggerTVX)) {
        SETBIT(rejectionMask, EventRejection::TvxTrigger);
      }
      /// time frame border cut
      if (useTimeFrameBorderCut && !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
        SETBIT(rejectionMask, EventRejection::TimeFrameBorderCut);
      }
      /// ITS rof border cut
      if (useItsRofBorderCut && !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
        SETBIT(rejectionMask, EventRejection::ItsRofBorderCut);
      }
      /// PVz consistency tracking - FT0 timing
      if (useIsGoodZvtxFT0vsPV && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
        SETBIT(rejectionMask, EventRejection::IsGoodZvtxFT0vsPV);
      }
      /// remove collisions in bunches with more than 1 reco collision
      /// POTENTIALLY BAD FOR BEAUTY ANALYSES
      if (useNoSameBunchPileup && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
        SETBIT(rejectionMask, EventRejection::NoSameBunchPileup);
      }
      /// No collisions in time range narrow
      if (useNoCollInTimeRangeNarrow && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeNarrow)) {
        SETBIT(rejectionMask, EventRejection::NoCollInTimeRangeNarrow);
      }
      /// No collisions in time range strict
      if (useNoCollInTimeRangeStandard && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
        SETBIT(rejectionMask, EventRejection::NoCollInTimeRangeStandard);
      }
      /// No collisions in ROF standard
      if (useNoCollInRofStandard && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) {
        SETBIT(rejectionMask, EventRejection::NoCollInRofStandard);
      }
      if (useOccupancyCut) {
        const auto occupancy = o2::hf_occupancy::getOccupancyColl(collision, occEstimator);
        if (occupancy < occupancyMin || occupancy > occupancyMax) {
          SETBIT(rejectionMask, EventRejection::Occupancy);
        }
      }
    }

    /// number of PV contributors
    if (collision.numContrib() < nPvContributorsMin) {
      SETBIT(rejectionMask, EventRejection::NContrib);
    }

    /// max PV chi2
    if (chi2PvMax > 0. && collision.chi2() > chi2PvMax) {
      SETBIT(rejectionMask, EventRejection::Chi2);
    }

    /// primary vertex z
    if (collision.posZ() < zPvPosMin || collision.posZ() > zPvPosMax) {
      SETBIT(rejectionMask, EventRejection::PositionZ);
    }

    if (!softwareTrigger.value.empty()) {
      if (!isInitCalled) {
        // protect against incomplete configuration
        LOG(fatal) << "Using Zorro utility w/o full HF event-selection configuration. Call the function HfEventSelection::init() to fix.";
      }
      // we might have to update it from CCDB
      const auto bc = collision.template bc_as<TBcs>();
      const auto runNumber = bc.runNumber();
      if (runNumber != currentRun) { // We might need to update Zorro from CCDB if the run number changes
        zorro.setCCDBpath(ccdbPathSoftwareTrigger);
        zorro.setBCtolerance(bcMarginForSoftwareTrigger);
        zorro.initCCDB(ccdb.service, runNumber, bc.timestamp(), softwareTrigger.value);
        currentRun = runNumber;
      }
      zorro.populateHistRegistry(registry, runNumber);

      if (softwareTrigger.value != "None") {
        if (!zorro.isSelected(bc.globalBC(), bcMarginForSoftwareTrigger)) { /// Just let Zorro do the accounting
          SETBIT(rejectionMask, EventRejection::SoftwareTrigger);
        }
      } else {
        if (!zorro.isNotSelectedByAny(bc.globalBC(), bcMarginForSoftwareTrigger)) { /// Just let Zorro do the accounting of not selected BCs
          SETBIT(rejectionMask, EventRejection::SoftwareTrigger);
        }
      }
    }

    return rejectionMask;
  }

  template <bool UseEvSel, o2::hf_centrality::CentralityEstimator CentEstimator, typename TBcs, typename TCollision>
  HfCollisionRejectionMask getHfCollisionRejectionMaskWithUpc(TCollision const& collision,
                                                              float& centrality,
                                                              o2::framework::Service<o2::ccdb::BasicCCDBManager> const& ccdb,
                                                              o2::framework::HistogramRegistry& registry,
                                                              TBcs const& bcs)
  {
    auto rejectionMaskWithUpc = getHfCollisionRejectionMask<UseEvSel, CentEstimator, TBcs>(collision, centrality, ccdb, registry);

    if (UseEvSel) {
      const SGCutParHolder sgCuts = setSgPreselection();
      const auto bc = collision.template foundBC_as<TBcs>();
      const auto bcRange = udhelpers::compatibleBCs(collision, sgCuts.NDtcoll(), bcs, sgCuts.minNBCs());
      const auto sgSelectionResult = sgSelector.IsSelected(sgCuts, collision, bcRange, bc);
      const int upcEventType = sgSelectionResult.value;
      if (upcEventType > o2::aod::sgselector::DoubleGap) {
        SETBIT(rejectionMaskWithUpc, EventRejection::UpcEventCut);
      } else {
        hUpCollisions->Fill(upcEventType);
      }
    }

    return rejectionMaskWithUpc;
  }

  /// \brief Fills histograms for monitoring event selections satisfied by the collision.
  /// \param collision analysed collision
  /// \param rejectionMask bitmask storing the info about which ev. selections are not satisfied by the collision
  template <typename TCollision>
  void fillHistograms(TCollision const& collision,
                      const HfCollisionRejectionMask rejectionMask,
                      const float centrality,
                      const float occupancy = -1.f,
                      const float ir = -1.f)
  {
    if (!areHistosInRegistry) {
      // protect against missing histograms in registry
      LOG(fatal) << "You are trying to fill histograms, but they are not in the histogram registry. Call the function HfEventSelectionMc::addHistograms() to fix.";
    }

    hCollisions->Fill(EventRejection::None);
    const auto posZ = collision.posZ();
    hPosZBeforeEvSel->Fill(posZ);

    for (std::size_t reason = 1; reason < EventRejection::NEventRejection; reason++) {
      if (TESTBIT(rejectionMask, reason)) {
        return;
      }
      hCollisions->Fill(reason);
    }

    hPosXAfterEvSel->Fill(collision.posX());
    hPosYAfterEvSel->Fill(collision.posY());
    hPosZAfterEvSel->Fill(posZ);
    hNumPvContributorsAfterSel->Fill(collision.numContrib());
    hSelCollisionsCent->Fill(centrality);
    hCollisionsCentOcc->Fill(centrality, occupancy);
    hCollisionsCentIR->Fill(centrality, ir);
  }

  template <typename TBc>
  double getInteractionRate(TBc const& bc,
                            o2::framework::Service<o2::ccdb::BasicCCDBManager> const& ccdb)
  {
    if (irSourceForCptFetcher.empty()) {
      o2::parameters::GRPLHCIFData* grpo = ccdb.service->getSpecificForRun<o2::parameters::GRPLHCIFData>("GLO/Config/GRPLHCIF", bc.runNumber());
      auto collsys = o2::common::core::CollisionSystemType::getCollisionTypeFromGrp(grpo);
      if (collsys == o2::common::core::CollisionSystemType::kCollSyspp) {
        irSourceForCptFetcher = std::string("T0VTX");
      } else {
        irSourceForCptFetcher = std::string("ZNC hadronic");
      }
    }

    return irFetcher.fetch(ccdb.service, bc.timestamp(), bc.runNumber(), irSourceForCptFetcher, true);
  }
};

struct HfEventSelectionMc {
  // event selection parameters (in chronological order of application)
  bool useSel8Trigger{false};                 // Apply the Sel8 selection
  bool useTvxTrigger{false};                  // Apply the TVX trigger
  bool useTimeFrameBorderCut{true};           // Apply TF border cut
  bool useItsRofBorderCut{false};             // Apply the ITS RO frame border cut
  float zPvPosMin{-1000.f};                   // Minimum PV posZ (cm)
  float zPvPosMax{1000.f};                    // Maximum PV posZ (cm)
  float centralityMin{-10.f};                 // Minimum centrality
  float centralityMax{100.f};                 // Maximum centrality
  bool requireGoodRct{false};                 // Apply RCT selection
  std::string rctLabel;                       // RCT selection flag
  bool rctCheckZDC{false};                    // require ZDC from RCT
  bool rctTreatLimitedAcceptanceAsBad{false}; // RCT flag to reject events with limited acceptance for selected detectors
  bool isInitCalled{false};                   // guard variable to guarantee full configuration, important for RCT
  bool areHistosInRegistry{false};            // guard variable to guarantee that histograms are added to the registry before filling them

  // util to retrieve the RCT info from CCDB
  o2::aod::rctsel::RCTFlagsChecker rctChecker;

  // histogram names
  static constexpr char NameHistGenCollisionsCent[] = "hGenCollisionsCent";
  static constexpr char NameHistRecCollisionsCentMc[] = "hRecCollisionsCentMc";
  static constexpr char NameHistNSplitVertices[] = "hNSplitVertices";
  static constexpr char NameHistGenCollisions[] = "hGenCollisions";

  std::shared_ptr<TH1> hGenCollisionsCent;
  std::shared_ptr<TH1> hRecCollisionsCentMc;
  std::shared_ptr<TH1> hNSplitVertices;
  std::shared_ptr<TH1> hGenCollisions;

  /// \brief Adds collision monitoring histograms in the histogram registry.
  /// \param registry reference to the histogram registry
  void addHistograms(o2::framework::HistogramRegistry& registry)
  {
    hGenCollisionsCent = registry.add<TH1>(NameHistGenCollisionsCent, "HF event counter;T0M;# of generated collisions", {o2::framework::HistType::kTH1D, {{100, 0., 100.}}});
    hRecCollisionsCentMc = registry.add<TH1>(NameHistRecCollisionsCentMc, "HF event counter;T0M;# of reconstructed collisions", {o2::framework::HistType::kTH1D, {{100, 0., 100.}}});
    hNSplitVertices = registry.add<TH1>(NameHistNSplitVertices, "HF split vertices counter;;# of reconstructed collisions per mc collision", {o2::framework::HistType::kTH1D, {{4, 1., 5.}}});
    hGenCollisions = registry.add<TH1>(NameHistGenCollisions, "HF event counter;;# of accepted collisions", {o2::framework::HistType::kTH1D, {axisEvents}});
    // Puts labels on the collision monitoring histogram.
    setEventRejectionLabels(hGenCollisions);

    // histograms in registry
    // let's update the guard variable
    areHistosInRegistry = true;
  }

  /// \brief Configures the object from the reco workflow
  /// \param registry reference to the histogram registry
  /// \param device device spec to get the configs from the reco workflow
  void configureFromDevice(o2::framework::DeviceSpec const& device)
  {
    for (const auto& option : device.options) {
      if (option.name == "hfEvSel.useSel8Trigger") {
        useSel8Trigger = option.defaultValue.get<bool>();
      } else if (option.name == "hfEvSel.useTvxTrigger") {
        useTvxTrigger = option.defaultValue.get<bool>();
      } else if (option.name == "hfEvSel.useTimeFrameBorderCut") {
        useTimeFrameBorderCut = option.defaultValue.get<bool>();
      } else if (option.name == "hfEvSel.useItsRofBorderCut") {
        useItsRofBorderCut = option.defaultValue.get<bool>();
      } else if (option.name == "hfEvSel.zPvPosMin") {
        zPvPosMin = option.defaultValue.get<float>();
      } else if (option.name == "hfEvSel.zPvPosMax") {
        zPvPosMax = option.defaultValue.get<float>();
      } else if (option.name == "hfEvSel.centralityMin") {
        centralityMin = option.defaultValue.get<float>();
      } else if (option.name == "hfEvSel.centralityMax") {
        centralityMax = option.defaultValue.get<float>();
      } else if (option.name == "hfEvSel.requireGoodRct") {
        requireGoodRct = option.defaultValue.get<bool>();
      } else if (option.name == "hfEvSel.rctLabel") {
        rctLabel = option.defaultValue.get<std::string>();
      } else if (option.name == "hfEvSel.rctCheckZDC") {
        rctCheckZDC = option.defaultValue.get<bool>();
      } else if (option.name == "hfEvSel.rctTreatLimitedAcceptanceAsBad") {
        rctTreatLimitedAcceptanceAsBad = option.defaultValue.get<bool>();
      }
    }
  }

  /// \brief Inits the HF event selection object
  /// \param device device spec to get the configs from the reco workflow
  /// \param registry reference to the histogram registry
  void init(o2::framework::DeviceSpec const& device,
            o2::framework::HistogramRegistry& registry)
  {
    // we get the configuration from the reco workflow
    configureFromDevice(device);

    // we initialise the RCT checker
    if (requireGoodRct) {
      rctChecker.init(rctLabel, rctCheckZDC, rctTreatLimitedAcceptanceAsBad);
    }

    // we initialise histograms
    addHistograms(registry);

    // full configuration complete: update the guard variable
    isInitCalled = true;
  }

  /// \brief Function to apply event selections to generated MC collisions
  /// \param mcCollision MC collision to test against the selection criteria
  /// \param collSlice collection of reconstructed collisions
  /// \param centrality centrality variable to be set in this function
  /// \return a bitmask with the event selections not satisfied by the analysed collision
  template <typename TBcs, o2::hf_centrality::CentralityEstimator CentEstimator, typename TCollisions, typename TMcCollision>
  HfCollisionRejectionMask getHfMcCollisionRejectionMask(TMcCollision const& mcCollision,
                                                         TCollisions const& collSlice,
                                                         float& centrality)
  {
    HfCollisionRejectionMask rejectionMask{};
    const auto zPv = mcCollision.posZ();
    const auto bc = mcCollision.template bc_as<TBcs>();

    if constexpr (CentEstimator != o2::hf_centrality::CentralityEstimator::None) {
      centrality = o2::hf_centrality::getCentralityGenColl(collSlice, CentEstimator);
      /// centrality selection
      if (centrality < centralityMin || centrality > centralityMax) {
        SETBIT(rejectionMask, EventRejection::Centrality);
      }
    }

    /// RCT condition
    if (requireGoodRct) {
      if (!isInitCalled) {
        // protect against incomplete configuration
        LOG(fatal) << "Checking RCT flags w/o full HF event-selection configuration (MC). Call the function HfEventSelectionMc::init() to fix.";
      }
      if (!rctChecker.checkTable(bc)) {
        SETBIT(rejectionMask, EventRejection::Rct);
      }
    }
    /// Sel8 trigger selection
    if (useSel8Trigger && (!bc.selection_bit(o2::aod::evsel::kIsTriggerTVX) || !bc.selection_bit(o2::aod::evsel::kNoTimeFrameBorder) || !bc.selection_bit(o2::aod::evsel::kNoITSROFrameBorder))) {
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
    /// ITS RO frame border cut
    if (useItsRofBorderCut && !bc.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      SETBIT(rejectionMask, EventRejection::ItsRofBorderCut);
    }
    /// primary vertex z
    if (zPv < zPvPosMin || zPv > zPvPosMax) {
      SETBIT(rejectionMask, EventRejection::PositionZ);
    }

    return rejectionMask;
  }

  /// \brief Fills histogram for monitoring event selections satisfied by the collision.
  /// \param collision analysed collision
  /// \param rejectionMask bitmask storing the info about which ev. selections are not satisfied by the collision
  template <o2::hf_centrality::CentralityEstimator CentEstimator, typename TMcCollision>
  void fillHistograms(TMcCollision const& mcCollision,
                      const HfCollisionRejectionMask rejectionMask,
                      const int nSplitColl = 0)
  {
    if (!areHistosInRegistry) {
      // protect against missing histograms in registry
      LOG(fatal) << "You are trying to fill histograms, but they are not in the histogram registry. Call the function HfEventSelectionMc::addHistograms() to fix.";
    }

    hGenCollisions->Fill(EventRejection::None);

    if constexpr (CentEstimator == o2::hf_centrality::CentralityEstimator::FT0M) {
      if (!TESTBIT(rejectionMask, EventRejection::TimeFrameBorderCut) && !TESTBIT(rejectionMask, EventRejection::ItsRofBorderCut) && !TESTBIT(rejectionMask, EventRejection::PositionZ)) {
        hGenCollisionsCent->Fill(mcCollision.centFT0M());
      }
    }

    for (std::size_t reason = 1; reason < EventRejection::NEventRejection; reason++) {
      if (TESTBIT(rejectionMask, reason)) {
        return;
      }
      hGenCollisions->Fill(reason);
    }

    if constexpr (CentEstimator == o2::hf_centrality::CentralityEstimator::FT0M) {
      hNSplitVertices->Fill(nSplitColl);
      for (int nColl = 0; nColl < nSplitColl; nColl++) {
        hRecCollisionsCentMc->Fill(mcCollision.centFT0M());
      }
    }
  }
};
} // namespace o2::hf_evsel

#endif // PWGHF_UTILS_UTILSEVSELHF_H_
