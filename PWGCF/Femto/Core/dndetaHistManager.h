// Copyright 2019-2022 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file dndetaHistManager.h
/// \brief histogram manager for charged-particle pseudorapidity density measurements on femto derived data
/// \author anton.riedel@tum.de, TU München, anton.riedel@tum.de

#ifndef PWGCF_FEMTO_CORE_DNDETAHISTMANAGER_H_
#define PWGCF_FEMTO_CORE_DNDETAHISTMANAGER_H_

#include "PWGCF/Femto/Core/histManager.h"

#include <CommonConstants/MathConstants.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>

#include <TAxis.h>
#include <TH1.h>
#include <THnSparse.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <map>
#include <string>
#include <string_view>
#include <vector>

namespace o2::analysis::femto::dndetahistmanager
{

// ---------------------------------------------------------------------------
// categories used as histogram axes, with labels
// ---------------------------------------------------------------------------

/// event selection steps, used for the cutflow histogram
enum EventCut : int {
  kEventAll = 0,
  kEventVtxZ,
  kEventCollisionSelection, ///< multiplicity/centrality/magnetic field ranges and collision bitmask
  kEventInelGt0,
  kEventTrigger, ///< pair/triplet trigger (always passed for minimum bias)
  kEventCutLast
};
constexpr std::array<const char*, kEventCutLast> EventCutNames = {"All", "Vertex z", "Collision selection", "INEL>0", "Trigger"};

/// reconstructed track types in data
enum DataTrackType : int {
  kGlobalPlusIts = 0, ///< all selected tracks
  kGlobalOnly,        ///< tracks with TPC
  kItsOnly,           ///< tracks without TPC
  kDataTrackTypeLast
};
constexpr std::array<const char*, kDataTrackTypeLast> DataTrackTypeNames = {"Global+ITS", "Global", "ITS only"};

/// generated particle species
enum GenType : int {
  kGenAll = 0,
  kGenPion,
  kGenKaon,
  kGenProton,
  kGenOther,
  kGenTypeLast
};
constexpr std::array<const char*, kGenTypeLast> GenTypeNames = {"All", "Pion", "Kaon", "Proton", "Other"};

/// variations of the generated pT spectrum below the extrapolation threshold
enum GenPtVariation : int {
  kPtNominal = 0,
  kPtUp,
  kPtDown,
  kGenPtVariationLast
};
constexpr std::array<const char*, kGenPtVariationLast> GenPtVariationNames = {"Nominal", "p_{T} up", "p_{T} down"};

/// classification of reconstructed tracks in MC
enum RecoType : int {
  kRecoAll = 0,
  kRecoPion,
  kRecoKaon,
  kRecoProton,
  kRecoOther,
  kRecoSecondary,
  kRecoWeakDecay,
  kRecoFake,
  kRecoBkg,
  kRecoTypeLast
};
constexpr std::array<const char*, kRecoTypeLast> RecoTypeNames = {"All", "Pion", "Kaon", "Proton", "Other", "Secondary", "Weak decay", "Fake", "No MC particle"};

/// generated event classes for event/signal loss
enum LossEvent : int {
  kLossGenAll = 0,
  kLossGenWithReco,
  kLossEventLast
};
constexpr std::array<const char*, kLossEventLast> LossEventNames = {"Generated", "Generated with selected reco"};

// ---------------------------------------------------------------------------
// histograms
// ---------------------------------------------------------------------------

enum DndetaHist {
  // event
  kEventCutflow,
  kPosZ,
  kCent,
  kPosZVsCent,
  // data
  kDataDndeta,
  kDataDndetaMb,
  // mc: generated
  kGenPosZ,
  kGenPosZVsCent,
  kGenAssocRecoPosZ,
  kGenAssocRecoPosZVsCent,
  kGenDndeta,
  kGenAssocRecoDndeta,
  kGenAssocRecoDndetaMb,
  // mc: reconstructed
  kRecoPosZ,
  kRecoCent,
  kRecoPosZVsCent,
  kRecoDndeta,
  kRecoDndetaMb,
  // mc: event and signal loss
  kLossEvents,
  kLossCentGen,
  kLossCentGenAssocReco,
  kLossMultGen,
  kLossMultGenAssocReco,
  kLossEtaGen,
  kLossEtaGenAssocReco,
  kLossEtaVsCentGen,
  kLossEtaVsCentGenAssocReco,
  kLossEtaVsMultGen,
  kLossEtaVsMultGenAssocReco,
  // multiplicity correlations
  kCorrNchVsFt0a,
  kCorrNchVsFt0c,
  kCorrNchVsNpv,
  kCorrNpvVsFt0c,
  // strangeness
  kStrangePosZVsCent,
  kK0shortCentEtaMass,
  kLambdaCentEtaMass,
  kAntiLambdaCentEtaMass,
  kDndetaHistLast
};

constexpr std::string_view EventDir = "Dndeta/Event/";
constexpr std::string_view DataDir = "Dndeta/Data/";
constexpr std::string_view GenDir = "Dndeta/Gen/";
constexpr std::string_view RecoDir = "Dndeta/Reco/";
constexpr std::string_view LossDir = "Dndeta/Loss/";
constexpr std::string_view CorrelationDir = "Dndeta/Correlation/";
constexpr std::string_view StrangenessDir = "Dndeta/Strangeness/";

constexpr std::array<histmanager::HistInfo<DndetaHist>, kDndetaHistLast> HistTable = {
  {
    // event
    {kEventCutflow, o2::framework::HistType::kTH1F, "hEventCutflow", "Event selection; ; Entries"},
    {kPosZ, o2::framework::HistType::kTH1F, "hPosZ", "Vertex z; V_{Z} (cm); Entries"},
    {kCent, o2::framework::HistType::kTH1F, "hCent", "Centrality; Centrality (%); Entries"},
    {kPosZVsCent, o2::framework::HistType::kTH2F, "hPosZVsCent", "Vertex z vs centrality; V_{Z} (cm); Centrality (%)"},
    // data
    {kDataDndeta, o2::framework::HistType::kTHnSparseF, "hDndeta", "Tracks; V_{Z} (cm); Centrality (%); #eta; #varphi; Track type"},
    {kDataDndetaMb, o2::framework::HistType::kTHnSparseF, "hDndetaMb", "Tracks (MB); V_{Z} (cm); #eta; #varphi"},
    // mc: generated
    {kGenPosZ, o2::framework::HistType::kTH1F, "hPosZ", "Generated vertex z; V_{Z,gen} (cm); Entries"},
    {kGenPosZVsCent, o2::framework::HistType::kTH2F, "hPosZVsCent", "Generated vertex z vs centrality; V_{Z,gen} (cm); Centrality_{gen} (%)"},
    {kGenAssocRecoPosZ, o2::framework::HistType::kTH1F, "hAssocRecoPosZ", "Generated vertex z (with selected reco); V_{Z,gen} (cm); Entries"},
    {kGenAssocRecoPosZVsCent, o2::framework::HistType::kTH2F, "hAssocRecoPosZVsCent", "Generated vertex z vs centrality (with selected reco); V_{Z,gen} (cm); Centrality_{gen} (%)"},
    {kGenDndeta, o2::framework::HistType::kTHnSparseF, "hDndeta", "Generated primaries; V_{Z,gen} (cm); Centrality_{gen} (%); #eta; #varphi"},
    {kGenAssocRecoDndeta, o2::framework::HistType::kTHnSparseF, "hAssocRecoDndeta", "Generated primaries (with selected reco); V_{Z,gen} (cm); Centrality_{gen} (%); #eta; #varphi; Species; p_{T} variation"},
    {kGenAssocRecoDndetaMb, o2::framework::HistType::kTHnSparseF, "hAssocRecoDndetaMb", "Generated primaries (with selected reco, MB); V_{Z,gen} (cm); #eta; #varphi; Species"},
    // mc: reconstructed
    {kRecoPosZ, o2::framework::HistType::kTH1F, "hPosZ", "Reconstructed vertex z; V_{Z} (cm); Entries"},
    {kRecoCent, o2::framework::HistType::kTH1F, "hCent", "Reconstructed centrality; Centrality (%); Entries"},
    {kRecoPosZVsCent, o2::framework::HistType::kTH2F, "hPosZVsCent", "Reconstructed vertex z vs centrality; V_{Z} (cm); Centrality (%)"},
    {kRecoDndeta, o2::framework::HistType::kTHnSparseF, "hDndeta", "Reconstructed tracks; V_{Z} (cm); Centrality (%); #eta; #varphi; Track class"},
    {kRecoDndetaMb, o2::framework::HistType::kTHnSparseF, "hDndetaMb", "Reconstructed tracks (MB); V_{Z} (cm); #eta; #varphi; Track class"},
    // mc: event and signal loss
    {kLossEvents, o2::framework::HistType::kTH1F, "hEvents", "Generated events; ; Entries"},
    {kLossCentGen, o2::framework::HistType::kTH1F, "hCentGen", "Generated events; Centrality_{gen} (%); Entries"},
    {kLossCentGenAssocReco, o2::framework::HistType::kTH1F, "hCentGenAssocReco", "Generated events with selected reco; Centrality_{gen} (%); Entries"},
    {kLossMultGen, o2::framework::HistType::kTH1F, "hMultGen", "Generated events; N_{ch,gen}; Entries"},
    {kLossMultGenAssocReco, o2::framework::HistType::kTH1F, "hMultGenAssocReco", "Generated events with selected reco; N_{ch,gen}; Entries"},
    {kLossEtaGen, o2::framework::HistType::kTH1F, "hEtaGen", "Generated primaries; #eta; Entries"},
    {kLossEtaGenAssocReco, o2::framework::HistType::kTH1F, "hEtaGenAssocReco", "Generated primaries in events with selected reco; #eta; Entries"},
    {kLossEtaVsCentGen, o2::framework::HistType::kTH2F, "hEtaVsCentGen", "Generated primaries; #eta; Centrality_{gen} (%)"},
    {kLossEtaVsCentGenAssocReco, o2::framework::HistType::kTH2F, "hEtaVsCentGenAssocReco", "Generated primaries in events with selected reco; #eta; Centrality_{gen} (%)"},
    {kLossEtaVsMultGen, o2::framework::HistType::kTH2F, "hEtaVsMultGen", "Generated primaries; #eta; N_{ch,gen}"},
    {kLossEtaVsMultGenAssocReco, o2::framework::HistType::kTH2F, "hEtaVsMultGenAssocReco", "Generated primaries in events with selected reco; #eta; N_{ch,gen}"},
    // multiplicity correlations
    {kCorrNchVsFt0a, o2::framework::HistType::kTH2F, "hNchVsFt0a", "N_{ch} vs FT0A; N_{ch}; FT0A amplitude"},
    {kCorrNchVsFt0c, o2::framework::HistType::kTH2F, "hNchVsFt0c", "N_{ch} vs FT0C; N_{ch}; FT0C amplitude"},
    {kCorrNchVsNpv, o2::framework::HistType::kTH2F, "hNchVsNpv", "N_{ch} vs N_{PV}; N_{ch}; N_{PV contributors}"},
    {kCorrNpvVsFt0c, o2::framework::HistType::kTH2F, "hNpvVsFt0c", "N_{PV} vs FT0C; N_{PV contributors}; FT0C amplitude"},
    // strangeness
    {kStrangePosZVsCent, o2::framework::HistType::kTH2F, "hPosZVsCent", "Vertex z vs centrality; V_{Z} (cm); Centrality (%)"},
    {kK0shortCentEtaMass, o2::framework::HistType::kTH3F, "hK0shortCentEtaMass", "K^{0}_{S}; Centrality (%); #eta; m_{#pi#pi} (GeV/#it{c}^{2})"},
    {kLambdaCentEtaMass, o2::framework::HistType::kTH3F, "hLambdaCentEtaMass", "#Lambda; Centrality (%); #eta; m_{p#pi} (GeV/#it{c}^{2})"},
    {kAntiLambdaCentEtaMass, o2::framework::HistType::kTH3F, "hAntiLambdaCentEtaMass", "#bar{#Lambda}; Centrality (%); #eta; m_{p#pi} (GeV/#it{c}^{2})"},
  }};

// ---------------------------------------------------------------------------
// binning
// ---------------------------------------------------------------------------

struct ConfDndetaBinning : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("DndetaBinning");
  o2::framework::ConfigurableAxis vtxZ{"vtxZ", {40, -20, 20}, "Vertex z binning"};
  o2::framework::ConfigurableAxis cent{"cent", {o2::framework::VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100}, "Centrality binning for the dN/deta histograms"};
  o2::framework::ConfigurableAxis centFine{"centFine", {100, 0, 100}, "Centrality binning for event and loss histograms"};
  o2::framework::ConfigurableAxis eta{"eta", {40, -2, 2}, "Pseudorapidity binning"};
  o2::framework::ConfigurableAxis phi{"phi", {o2::framework::VARIABLE_WIDTH, 0, o2::constants::math::PIQuarter, o2::constants::math::PIHalf, o2::constants::math::PIQuarter * 3., o2::constants::math::PI, o2::constants::math::PIQuarter * 5., o2::constants::math::PIHalf * 3., o2::constants::math::PIQuarter * 7., o2::constants::math::TwoPI}, "Azimuth binning for the dN/deta histograms"};
  o2::framework::ConfigurableAxis multGen{"multGen", {500, 0, 500}, "Binning of the generated multiplicity estimator of the loss histograms (eta range: Dndeta.multEtaMax)"};
  o2::framework::ConfigurableAxis nch{"nch", {501, -0.5, 500.5}, "Number of selected tracks binning (correlations)"};
  o2::framework::ConfigurableAxis npv{"npv", {501, -0.5, 500.5}, "Number of PV contributors binning (correlations)"};
  o2::framework::ConfigurableAxis ft0a{"ft0a", {501, -0.5, 500.5}, "FT0A amplitude binning (correlations)"};
  o2::framework::ConfigurableAxis ft0c{"ft0c", {501, -0.5, 500.5}, "FT0C amplitude binning (correlations)"};
  o2::framework::ConfigurableAxis massK0short{"massK0short", {200, 0.4, 0.6}, "K0short mass binning"};
  o2::framework::ConfigurableAxis massLambda{"massLambda", {200, 1.07, 1.17}, "Lambda mass binning"};
};

namespace detail
{
/// one bin per category, centered on the integer value of the enum; labels are set after booking
template <std::size_t N>
o2::framework::AxisSpec makeCategoryAxis(std::array<const char*, N> const& /*names*/, const char* title)
{
  return {static_cast<int>(N), -0.5, static_cast<double>(N) - 0.5, title};
}
} // namespace detail

template <typename T>
auto makeDndetaHistSpecMap(T const& conf)
{
  const o2::framework::AxisSpec axisTrackType = detail::makeCategoryAxis(DataTrackTypeNames, "Track type");
  const o2::framework::AxisSpec axisGenType = detail::makeCategoryAxis(GenTypeNames, "Species");
  const o2::framework::AxisSpec axisPtVariation = detail::makeCategoryAxis(GenPtVariationNames, "p_{T} variation");
  const o2::framework::AxisSpec axisRecoType = detail::makeCategoryAxis(RecoTypeNames, "Track class");
  const o2::framework::AxisSpec axisEventCut = detail::makeCategoryAxis(EventCutNames, "");
  const o2::framework::AxisSpec axisLossEvent = detail::makeCategoryAxis(LossEventNames, "");

  return std::map<DndetaHist, std::vector<o2::framework::AxisSpec>>{
    // event
    {kEventCutflow, {axisEventCut}},
    {kPosZ, {conf.vtxZ}},
    {kCent, {conf.centFine}},
    {kPosZVsCent, {conf.vtxZ, conf.cent}},
    // data
    {kDataDndeta, {conf.vtxZ, conf.cent, conf.eta, conf.phi, axisTrackType}},
    {kDataDndetaMb, {conf.vtxZ, conf.eta, conf.phi}},
    // mc: generated
    {kGenPosZ, {conf.vtxZ}},
    {kGenPosZVsCent, {conf.vtxZ, conf.cent}},
    {kGenAssocRecoPosZ, {conf.vtxZ}},
    {kGenAssocRecoPosZVsCent, {conf.vtxZ, conf.cent}},
    {kGenDndeta, {conf.vtxZ, conf.cent, conf.eta, conf.phi}},
    {kGenAssocRecoDndeta, {conf.vtxZ, conf.cent, conf.eta, conf.phi, axisGenType, axisPtVariation}},
    {kGenAssocRecoDndetaMb, {conf.vtxZ, conf.eta, conf.phi, axisGenType}},
    // mc: reconstructed
    {kRecoPosZ, {conf.vtxZ}},
    {kRecoCent, {conf.centFine}},
    {kRecoPosZVsCent, {conf.vtxZ, conf.cent}},
    {kRecoDndeta, {conf.vtxZ, conf.cent, conf.eta, conf.phi, axisRecoType}},
    {kRecoDndetaMb, {conf.vtxZ, conf.eta, conf.phi, axisRecoType}},
    // mc: event and signal loss
    {kLossEvents, {axisLossEvent}},
    {kLossCentGen, {conf.centFine}},
    {kLossCentGenAssocReco, {conf.centFine}},
    {kLossMultGen, {conf.multGen}},
    {kLossMultGenAssocReco, {conf.multGen}},
    {kLossEtaGen, {conf.eta}},
    {kLossEtaGenAssocReco, {conf.eta}},
    {kLossEtaVsCentGen, {conf.eta, conf.centFine}},
    {kLossEtaVsCentGenAssocReco, {conf.eta, conf.centFine}},
    {kLossEtaVsMultGen, {conf.eta, conf.multGen}},
    {kLossEtaVsMultGenAssocReco, {conf.eta, conf.multGen}},
    // multiplicity correlations
    {kCorrNchVsFt0a, {conf.nch, conf.ft0a}},
    {kCorrNchVsFt0c, {conf.nch, conf.ft0c}},
    {kCorrNchVsNpv, {conf.nch, conf.npv}},
    {kCorrNpvVsFt0c, {conf.npv, conf.ft0c}},
    // strangeness
    {kStrangePosZVsCent, {conf.vtxZ, conf.cent}},
    {kK0shortCentEtaMass, {conf.centFine, conf.eta, conf.massK0short}},
    {kLambdaCentEtaMass, {conf.centFine, conf.eta, conf.massLambda}},
    {kAntiLambdaCentEtaMass, {conf.centFine, conf.eta, conf.massLambda}},
  };
}

/// blocks of histograms, enabled at runtime depending on the active process functions
enum Block : uint32_t {
  kBlockEvent = 1u << 0,
  kBlockData = 1u << 1,
  kBlockMcEfficiency = 1u << 2,
  kBlockLoss = 1u << 3,
  kBlockCorrelation = 1u << 4,
  kBlockStrangeness = 1u << 5,
};

class DndetaHistManager
{
 public:
  DndetaHistManager() = default;
  ~DndetaHistManager() = default;

  /// \param blocks bitwise OR of Block values
  void init(o2::framework::HistogramRegistry* registry,
            std::map<DndetaHist, std::vector<o2::framework::AxisSpec>> const& Specs,
            uint32_t blocks)
  {
    mHistogramRegistry = registry;
    mBlocks = blocks;
    if ((mBlocks & kBlockEvent) != 0u) {
      initEvent(Specs);
    }
    if ((mBlocks & kBlockData) != 0u) {
      initData(Specs);
    }
    if ((mBlocks & kBlockMcEfficiency) != 0u) {
      initMcEfficiency(Specs);
    }
    if ((mBlocks & kBlockLoss) != 0u) {
      initLoss(Specs);
    }
    if ((mBlocks & kBlockCorrelation) != 0u) {
      initCorrelation(Specs);
    }
    if ((mBlocks & kBlockStrangeness) != 0u) {
      initStrangeness(Specs);
    }
  }

  // ---- event --------------------------------------------------------------

  void fillEventCutflow(EventCut cut)
  {
    mHistogramRegistry->fill(HIST(EventDir) + HIST(getHistName(kEventCutflow, HistTable)), static_cast<float>(cut));
  }

  void fillEvent(float posZ, float cent)
  {
    mHistogramRegistry->fill(HIST(EventDir) + HIST(getHistName(kPosZ, HistTable)), posZ);
    mHistogramRegistry->fill(HIST(EventDir) + HIST(getHistName(kCent, HistTable)), cent);
    mHistogramRegistry->fill(HIST(EventDir) + HIST(getHistName(kPosZVsCent, HistTable)), posZ, cent);
  }

  // ---- data ---------------------------------------------------------------

  /// fill a selected data track (track QA is done with the track histogram managers)
  void fillDataTrack(float posZ, float cent, float eta, float phi, bool hasTpc)
  {
    mHistogramRegistry->fill(HIST(DataDir) + HIST(getHistName(kDataDndetaMb, HistTable)), posZ, eta, phi);
    mHistogramRegistry->fill(HIST(DataDir) + HIST(getHistName(kDataDndeta, HistTable)), posZ, cent, eta, phi, static_cast<double>(kGlobalPlusIts));
    const DataTrackType type = hasTpc ? kGlobalOnly : kItsOnly;
    mHistogramRegistry->fill(HIST(DataDir) + HIST(getHistName(kDataDndeta, HistTable)), posZ, cent, eta, phi, static_cast<double>(type));
  }

  // ---- mc: generated ------------------------------------------------------

  void fillGenCollision(float posZ, float cent, bool hasSelectedReco)
  {
    mHistogramRegistry->fill(HIST(GenDir) + HIST(getHistName(kGenPosZ, HistTable)), posZ);
    mHistogramRegistry->fill(HIST(GenDir) + HIST(getHistName(kGenPosZVsCent, HistTable)), posZ, cent);
    if (hasSelectedReco) {
      mHistogramRegistry->fill(HIST(GenDir) + HIST(getHistName(kGenAssocRecoPosZ, HistTable)), posZ);
      mHistogramRegistry->fill(HIST(GenDir) + HIST(getHistName(kGenAssocRecoPosZVsCent, HistTable)), posZ, cent);
    }
  }

  /// fill a selected generated primary
  /// weightUp/weightDown are the weights of the pT variations (uncertainty of the extrapolation to pT = 0),
  /// they are computed by the dndeta builder
  void fillGenParticle(float posZ, float cent, float eta, float phi, GenType species, bool hasSelectedReco, double weightUp, double weightDown)
  {
    mHistogramRegistry->fill(HIST(GenDir) + HIST(getHistName(kGenDndeta, HistTable)), posZ, cent, eta, phi);
    if (!hasSelectedReco) {
      return;
    }
    const auto all = static_cast<double>(kGenAll);
    mHistogramRegistry->fill(HIST(GenDir) + HIST(getHistName(kGenAssocRecoDndeta, HistTable)), posZ, cent, eta, phi, all, static_cast<double>(kPtNominal));
    mHistogramRegistry->fill(HIST(GenDir) + HIST(getHistName(kGenAssocRecoDndetaMb, HistTable)), posZ, eta, phi, all);
    mHistogramRegistry->fill(HIST(GenDir) + HIST(getHistName(kGenAssocRecoDndeta, HistTable)), posZ, cent, eta, phi, all, static_cast<double>(kPtUp), weightUp);
    mHistogramRegistry->fill(HIST(GenDir) + HIST(getHistName(kGenAssocRecoDndeta, HistTable)), posZ, cent, eta, phi, all, static_cast<double>(kPtDown), weightDown);

    const auto spec = static_cast<double>(species);
    mHistogramRegistry->fill(HIST(GenDir) + HIST(getHistName(kGenAssocRecoDndeta, HistTable)), posZ, cent, eta, phi, spec, static_cast<double>(kPtNominal));
    mHistogramRegistry->fill(HIST(GenDir) + HIST(getHistName(kGenAssocRecoDndetaMb, HistTable)), posZ, eta, phi, spec);
  }

  // ---- mc: reconstructed --------------------------------------------------

  void fillRecoCollision(float posZ, float cent)
  {
    mHistogramRegistry->fill(HIST(RecoDir) + HIST(getHistName(kRecoPosZ, HistTable)), posZ);
    mHistogramRegistry->fill(HIST(RecoDir) + HIST(getHistName(kRecoCent, HistTable)), cent);
    mHistogramRegistry->fill(HIST(RecoDir) + HIST(getHistName(kRecoPosZVsCent, HistTable)), posZ, cent);
  }

  /// eta/phi are the reconstructed ones for kRecoAll and kRecoBkg, the generated ones otherwise (as in the legacy task)
  void fillRecoTrack(float posZ, float cent, float eta, float phi, RecoType type)
  {
    mHistogramRegistry->fill(HIST(RecoDir) + HIST(getHistName(kRecoDndeta, HistTable)), posZ, cent, eta, phi, static_cast<double>(type));
    mHistogramRegistry->fill(HIST(RecoDir) + HIST(getHistName(kRecoDndetaMb, HistTable)), posZ, eta, phi, static_cast<double>(type));
  }

  // ---- mc: event and signal loss ------------------------------------------

  void fillLossEvent(float cent, float multGen, bool hasSelectedReco)
  {
    mHistogramRegistry->fill(HIST(LossDir) + HIST(getHistName(kLossEvents, HistTable)), static_cast<float>(kLossGenAll));
    mHistogramRegistry->fill(HIST(LossDir) + HIST(getHistName(kLossCentGen, HistTable)), cent);
    mHistogramRegistry->fill(HIST(LossDir) + HIST(getHistName(kLossMultGen, HistTable)), multGen);
    if (hasSelectedReco) {
      mHistogramRegistry->fill(HIST(LossDir) + HIST(getHistName(kLossEvents, HistTable)), static_cast<float>(kLossGenWithReco));
      mHistogramRegistry->fill(HIST(LossDir) + HIST(getHistName(kLossCentGenAssocReco, HistTable)), cent);
      mHistogramRegistry->fill(HIST(LossDir) + HIST(getHistName(kLossMultGenAssocReco, HistTable)), multGen);
    }
  }

  void fillLossParticle(float eta, float cent, float multGen, bool hasSelectedReco)
  {
    mHistogramRegistry->fill(HIST(LossDir) + HIST(getHistName(kLossEtaGen, HistTable)), eta);
    mHistogramRegistry->fill(HIST(LossDir) + HIST(getHistName(kLossEtaVsCentGen, HistTable)), eta, cent);
    mHistogramRegistry->fill(HIST(LossDir) + HIST(getHistName(kLossEtaVsMultGen, HistTable)), eta, multGen);
    if (hasSelectedReco) {
      mHistogramRegistry->fill(HIST(LossDir) + HIST(getHistName(kLossEtaGenAssocReco, HistTable)), eta);
      mHistogramRegistry->fill(HIST(LossDir) + HIST(getHistName(kLossEtaVsCentGenAssocReco, HistTable)), eta, cent);
      mHistogramRegistry->fill(HIST(LossDir) + HIST(getHistName(kLossEtaVsMultGenAssocReco, HistTable)), eta, multGen);
    }
  }

  // ---- multiplicity correlations ------------------------------------------

  void fillCorrelation(float nch, float npv, float ft0a, float ft0c)
  {
    mHistogramRegistry->fill(HIST(CorrelationDir) + HIST(getHistName(kCorrNchVsFt0a, HistTable)), nch, ft0a);
    mHistogramRegistry->fill(HIST(CorrelationDir) + HIST(getHistName(kCorrNchVsFt0c, HistTable)), nch, ft0c);
    mHistogramRegistry->fill(HIST(CorrelationDir) + HIST(getHistName(kCorrNchVsNpv, HistTable)), nch, npv);
    mHistogramRegistry->fill(HIST(CorrelationDir) + HIST(getHistName(kCorrNpvVsFt0c, HistTable)), npv, ft0c);
  }

  // ---- strangeness --------------------------------------------------------

  void fillStrangenessCollision(float posZ, float cent)
  {
    mHistogramRegistry->fill(HIST(StrangenessDir) + HIST(getHistName(kStrangePosZVsCent, HistTable)), posZ, cent);
  }

  void fillK0short(float cent, float eta, float mass)
  {
    mHistogramRegistry->fill(HIST(StrangenessDir) + HIST(getHistName(kK0shortCentEtaMass, HistTable)), cent, eta, mass);
  }

  void fillLambda(float cent, float eta, float mass, bool isAntiLambda)
  {
    if (isAntiLambda) {
      mHistogramRegistry->fill(HIST(StrangenessDir) + HIST(getHistName(kAntiLambdaCentEtaMass, HistTable)), cent, eta, mass);
    } else {
      mHistogramRegistry->fill(HIST(StrangenessDir) + HIST(getHistName(kLambdaCentEtaMass, HistTable)), cent, eta, mass);
    }
  }

 private:
  void add(std::string_view dir, DndetaHist hist, std::map<DndetaHist, std::vector<o2::framework::AxisSpec>> const& Specs)
  {
    mHistogramRegistry->add(std::string(dir) + getHistNameV2(hist, HistTable), getHistDesc(hist, HistTable), getHistType(hist, HistTable), {Specs.at(hist)});
  }

  template <std::size_t N>
  static void setLabels(TAxis* axis, std::array<const char*, N> const& names)
  {
    for (std::size_t i = 0; i < N; i++) {
      axis->SetBinLabel(static_cast<int>(i) + 1, names[i]);
    }
  }

  void initEvent(std::map<DndetaHist, std::vector<o2::framework::AxisSpec>> const& Specs)
  {
    add(EventDir, kEventCutflow, Specs);
    setLabels(mHistogramRegistry->get<TH1>(HIST(EventDir) + HIST(getHistName(kEventCutflow, HistTable)))->GetXaxis(), EventCutNames);
    add(EventDir, kPosZ, Specs);
    add(EventDir, kCent, Specs);
    add(EventDir, kPosZVsCent, Specs);
  }

  void initData(std::map<DndetaHist, std::vector<o2::framework::AxisSpec>> const& Specs)
  {
    add(DataDir, kDataDndeta, Specs);
    setLabels(mHistogramRegistry->get<THnSparse>(HIST(DataDir) + HIST(getHistName(kDataDndeta, HistTable)))->GetAxis(4), DataTrackTypeNames);
    add(DataDir, kDataDndetaMb, Specs);
  }

  void initMcEfficiency(std::map<DndetaHist, std::vector<o2::framework::AxisSpec>> const& Specs)
  {
    add(GenDir, kGenPosZ, Specs);
    add(GenDir, kGenPosZVsCent, Specs);
    add(GenDir, kGenAssocRecoPosZ, Specs);
    add(GenDir, kGenAssocRecoPosZVsCent, Specs);
    add(GenDir, kGenDndeta, Specs);
    add(GenDir, kGenAssocRecoDndeta, Specs);
    setLabels(mHistogramRegistry->get<THnSparse>(HIST(GenDir) + HIST(getHistName(kGenAssocRecoDndeta, HistTable)))->GetAxis(4), GenTypeNames);
    setLabels(mHistogramRegistry->get<THnSparse>(HIST(GenDir) + HIST(getHistName(kGenAssocRecoDndeta, HistTable)))->GetAxis(5), GenPtVariationNames);
    add(GenDir, kGenAssocRecoDndetaMb, Specs);
    setLabels(mHistogramRegistry->get<THnSparse>(HIST(GenDir) + HIST(getHistName(kGenAssocRecoDndetaMb, HistTable)))->GetAxis(3), GenTypeNames);

    add(RecoDir, kRecoPosZ, Specs);
    add(RecoDir, kRecoCent, Specs);
    add(RecoDir, kRecoPosZVsCent, Specs);
    add(RecoDir, kRecoDndeta, Specs);
    setLabels(mHistogramRegistry->get<THnSparse>(HIST(RecoDir) + HIST(getHistName(kRecoDndeta, HistTable)))->GetAxis(4), RecoTypeNames);
    add(RecoDir, kRecoDndetaMb, Specs);
    setLabels(mHistogramRegistry->get<THnSparse>(HIST(RecoDir) + HIST(getHistName(kRecoDndetaMb, HistTable)))->GetAxis(3), RecoTypeNames);
  }

  void initLoss(std::map<DndetaHist, std::vector<o2::framework::AxisSpec>> const& Specs)
  {
    add(LossDir, kLossEvents, Specs);
    setLabels(mHistogramRegistry->get<TH1>(HIST(LossDir) + HIST(getHistName(kLossEvents, HistTable)))->GetXaxis(), LossEventNames);
    add(LossDir, kLossCentGen, Specs);
    add(LossDir, kLossCentGenAssocReco, Specs);
    add(LossDir, kLossMultGen, Specs);
    add(LossDir, kLossMultGenAssocReco, Specs);
    add(LossDir, kLossEtaGen, Specs);
    add(LossDir, kLossEtaGenAssocReco, Specs);
    add(LossDir, kLossEtaVsCentGen, Specs);
    add(LossDir, kLossEtaVsCentGenAssocReco, Specs);
    add(LossDir, kLossEtaVsMultGen, Specs);
    add(LossDir, kLossEtaVsMultGenAssocReco, Specs);
  }

  void initCorrelation(std::map<DndetaHist, std::vector<o2::framework::AxisSpec>> const& Specs)
  {
    add(CorrelationDir, kCorrNchVsFt0a, Specs);
    add(CorrelationDir, kCorrNchVsFt0c, Specs);
    add(CorrelationDir, kCorrNchVsNpv, Specs);
    add(CorrelationDir, kCorrNpvVsFt0c, Specs);
  }

  void initStrangeness(std::map<DndetaHist, std::vector<o2::framework::AxisSpec>> const& Specs)
  {
    add(StrangenessDir, kStrangePosZVsCent, Specs);
    add(StrangenessDir, kK0shortCentEtaMass, Specs);
    add(StrangenessDir, kLambdaCentEtaMass, Specs);
    add(StrangenessDir, kAntiLambdaCentEtaMass, Specs);
  }

  o2::framework::HistogramRegistry* mHistogramRegistry = nullptr;
  uint32_t mBlocks = 0;
};

} // namespace o2::analysis::femto::dndetahistmanager

#endif // PWGCF_FEMTO_CORE_DNDETAHISTMANAGER_H_
