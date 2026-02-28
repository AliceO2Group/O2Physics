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

/// \file DataCreationCharmReso.h
/// \brief utility functions for charm-hadron resonance derived data creators
///
/// \author Luca Aglietta <luca.aglietta@cern.ch>, UniTO Turin
/// \author Fabrizio Grosa <fabrizio.grosa@cern.ch>, CERN

#ifndef PWGHF_D2H_CORE_DATACREATIONCHARMRESO_H_
#define PWGHF_D2H_CORE_DATACREATIONCHARMRESO_H_

#ifndef HomogeneousField
#define HomogeneousField // needed for KFParticle::SetField(magneticField);
#endif

#include "PWGEM/Dilepton/Utils/PairUtilities.h"
#include "PWGEM/PhotonMeson/Utils/PCMUtilities.h"
#include "PWGEM/PhotonMeson/Utils/TrackSelection.h"
#include "PWGHF/Core/CentralityEstimation.h"
#include "PWGHF/Core/DecayChannels.h"
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/Utils/utilsEvSelHf.h"
#include "PWGHF/Utils/utilsMcMatching.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TPCVDriftManager.h"
#include "Common/Core/trackUtilities.h"
#include "Tools/KFparticle/KFUtilities.h"

#include <CommonConstants/PhysicsConstants.h>
#include <DCAFitter/DCAFitterN.h>
#include <DetectorsBase/Propagator.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/O2DatabasePDGPlugin.h>

#include <TH1.h>
#include <TPDGCode.h>

#include <KFPTrack.h>
#include <KFPVertex.h>
#include <KFParticle.h>

#include <Rtypes.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <map>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <vector>

namespace o2::analysis
{
namespace hf_charm_reso
{

// event types
enum EventType : uint8_t {
  Processed = 0,
  NoDV0Selected,
  DV0Selected,
  NEventType
};

enum BachelorType : uint8_t {
  K0s = 0,
  Lambda,
  AntiLambda,
  Track,
  Gamma,
  NBachelorTypes
};

enum DMesonType : uint8_t {
  Dplus = 1,
  Dstar,
  D0,
  NDMesonType
};

enum PairingType : uint8_t {
  V0Only,
  TrackOnly,
  V0AndTrack,
  GammaOnly,
  NPairingType
};

enum D0Sel : uint8_t {
  SelectedD0 = 0,
  SelectedD0Bar,
  ND0Sel
};

enum TrackSel : uint8_t {
  GlobalTrack = 0,
  GlobalTrackWoDca,
  QualityTrackITS
};

// Helper structs to pass D and V0 informations

struct HfResoCandidateV0 {
  std::array<float, 3> pos = {0.f, 0.f, 0.f};
  std::array<float, 3> mom = {0.f, 0.f, 0.f};
  std::array<float, 3> momPos = {0.f, 0.f, 0.f};
  std::array<float, 3> momNeg = {0.f, 0.f, 0.f};
  float pT = -1.f;
  float cosPA = -2.f;
  float dcaV0ToPv = 1000.f;
  float dcaDau = 1000.f;
  float alpha = -1.f;
  float qt = -1.f;
  float eta = -999.f;
  float radius = 0.f;
  float mK0Short = 0.f;
  float mLambda = 0.f;
  uint8_t v0Type = 0u;
};

struct HfResoVarContainer {
  float invMassD = 0.f;
  float ptD = -1.f;
  float invMassD0 = 0.f;
  float invMassD0Bar = 0.f;
  float invMassReso = 0.f;
  float invMassResoBar = 0.f;
  float ptReso = -1.f;
  int8_t signD = 0;
  std::array<float, 3> pVectorProng0 = {0.f, 0.f, 0.f};
  std::array<float, 3> pVectorProng1 = {0.f, 0.f, 0.f};
  std::array<float, 3> pVectorProng2 = {0.f, 0.f, 0.f};
};

struct HfResoConfigV0Cuts : o2::framework::ConfigurableGroup {
  std::string prefix = "v0s"; // JSON group name
  o2::framework::Configurable<float> deltaMassK0s{"deltaMassK0s", 0.02, "delta mass cut for K0S"};
  o2::framework::Configurable<float> deltaMassLambda{"deltaMassLambda", 0.01, "delta mass cut for Lambda"};
  o2::framework::Configurable<float> etaMax{"etaMax", 0.8f, "maximum eta"};
  o2::framework::Configurable<float> etaMaxDau{"etaMaxDau", 5.f, "maximum eta V0 daughters"};
  o2::framework::Configurable<float> trackNclusItsCut{"trackNclusItsCut", 0, "Minimum number of ITS clusters for V0 daughter"};
  o2::framework::Configurable<int> trackNCrossedRowsTpc{"trackNCrossedRowsTpc", 50, "Minimum TPC crossed rows"};
  o2::framework::Configurable<float> trackNsharedClusTpc{"trackNsharedClusTpc", 1000, "Maximum number of shared TPC clusters for V0 daughter"};
  o2::framework::Configurable<float> trackFracMaxindableTpcCls{"trackFracMaxindableTpcCls", 0.8f, "Maximum fraction of findable TPC clusters for V0 daughter"};
  o2::framework::Configurable<float> dcaDau{"dcaDau", 1.f, "DCA V0 daughters"};
  o2::framework::Configurable<float> dcaMaxDauToPv{"dcaMaxDauToPv", 0.1f, "Maximum daughter's DCA to PV"};
  o2::framework::Configurable<float> dcaPv{"dcaPv", 1.f, "DCA V0 to PV"};
  o2::framework::Configurable<double> cosPa{"cosPa", 0.99f, "V0 CosPA"};
  o2::framework::Configurable<float> radiusMin{"radiusMin", 0.9f, "Minimum v0 radius accepted"};
  o2::framework::Configurable<float> nSigmaTpc{"nSigmaTpc", 4.f, "Nsigmatpc"};
  o2::framework::Configurable<float> nSigmaTofPr{"nSigmaTofPr", 4.f, "N sigma TOF for protons only"};
  o2::framework::Configurable<bool> propagateV0toPV{"propagateV0toPV", false, "Enable or disable V0 propagation to V0"};
};

struct HfResoConfigGammaCuts : o2::framework::ConfigurableGroup {
  std::string prefix = "gammas"; // JSON group name
  o2::framework::Configurable<float> etaMax{"etaMax", 0.8f, "maximum eta"};
  o2::framework::Configurable<float> ptMin{"ptMin", 0.1f, "minimum pT"};
  o2::framework::Configurable<float> ptMaxItsOnly{"ptMaxItsOnly", 0.3f, "maximum pT for ITS-only gammas"};
  o2::framework::Configurable<float> etaMaxDau{"etaMaxDau", 1.f, "maximum eta gamma daughters"};
  o2::framework::Configurable<float> trackNclusItsCut{"trackNclusItsCut", 0, "Minimum number of ITS clusters for gamma daughter"};
  o2::framework::Configurable<int> trackNCrossedRowsTpc{"trackNCrossedRowsTpc", 50, "Minimum TPC crossed rows"};
  o2::framework::Configurable<float> trackNsharedClusTpc{"trackNsharedClusTpc", 1000, "Maximum number of shared TPC clusters for gamma daughter"};
  o2::framework::Configurable<float> trackFracMaxindableTpcCls{"trackFracMaxindableTpcCls", 0.8f, "Maximum fraction of findable TPC clusters for gamma daughter"};
  o2::framework::Configurable<float> dcaDauIts{"dcaDauIts", 0.5f, "maximum DCA gamma daughters (ITS)"};
  o2::framework::Configurable<float> dcaDauItsIb{"dcaDauItsIb", 1.0f, "maximum DCA gamma daughters (ITS IB)"};
  o2::framework::Configurable<float> dcaDauTpc{"dcaDauTpc", 0.5f, "maximum DCA gamma daughters (TPC)"};
  o2::framework::Configurable<float> dcaDauTpcInner{"dcaDauTpcInner", 1.0f, "maximum DCA gamma daughters (TPC inner)"};
  o2::framework::Configurable<float> dcaMaxDauToPv{"dcaMaxDauToPv", 0.1f, "Maximum gamma daughter's DCA to PV"};
  o2::framework::Configurable<float> dcaPv{"dcaPv", 1.f, "DCA gamma to PV"};
  o2::framework::Configurable<double> cosPa{"cosPa", 0.99f, "gamma CosPA"};
  o2::framework::Configurable<float> radiusMin{"radiusMin", 1.0f, "Minimum gamma radius accepted"};
  o2::framework::Configurable<float> radiusMax{"radiusMax", 90.f, "Maximum gamma radius accepted"};
  o2::framework::Configurable<float> alphaApMax{"alphaApMax", 0.95f, "Maximum alpha AP"};
  o2::framework::Configurable<float> qtApMax{"qtApMax", 0.01f, "Maximum qt AP"};
  o2::framework::Configurable<float> nSigmaTpcEl{"nSigmaTpcEl", 4.f, "N sigma TPC for electrons"};
  o2::framework::Configurable<bool> propagateGammatoPV{"propagateGammatoPV", false, "Enable or disable V0 propagation to V0"};
};

struct HfResoConfigSingleTrackCuts : o2::framework::ConfigurableGroup {
  std::string prefix = "singleTracks"; // JSON group name
  o2::framework::Configurable<int> setTrackSelections{"setTrackSelections", 2, "flag to apply track selections: 0=none; 1=global track w/o DCA selection; 2=global track; 3=only ITS quality"};
  o2::framework::Configurable<float> maxEta{"maxEta", 0.8, "maximum pseudorapidity for single tracks to be paired with D mesons"};
  o2::framework::Configurable<float> minPt{"minPt", 0.1, "minimum pT for single tracks to be paired with D mesons"};
  o2::framework::Configurable<float> maxNsigmaTpcPi{"maxNsigmaTpcPi", -1., "maximum pion NSigma in TPC for single tracks to be paired with D mesons; set negative to reject"};
  o2::framework::Configurable<float> maxNsigmaTpcKa{"maxNsigmaTpcKa", -1., "maximum kaon NSigma in TPC for single tracks to be paired with D mesons; set negative to reject"};
  o2::framework::Configurable<float> maxNsigmaTpcPr{"maxNsigmaTpcPr", 3., "maximum proton NSigma in TPC for single tracks to be paired with D mesons; set negative to reject"};
};

struct HfResoConfigQaPlots : o2::framework::ConfigurableGroup {
  std::string prefix = "qaPlots"; // JSON group name
  o2::framework::Configurable<bool> applyCutsForQaHistograms{"applyCutsForQaHistograms", true, "flag to apply cuts to QA histograms"};
  o2::framework::Configurable<float> cutMassDMin{"cutMassDMin", 1.83, "minimum mass for D0 and Dplus candidates"};
  o2::framework::Configurable<float> cutMassDMax{"cutMassDMax", 1.92, "maximum mass for D0 and Dplus candidates"};
  o2::framework::Configurable<float> cutMassDstarMin{"cutMassDstarMin", 0.139, "minimum mass for Dstar candidates"}; // o2-linter: disable=pdg/explicit-mass (false positive)
  o2::framework::Configurable<float> cutMassDstarMax{"cutMassDstarMax", 0.175, "maximum mass for Dstar candidates"};
  o2::framework::Configurable<float> cutMassK0sMin{"cutMassK0sMin", 0.485, "minimum mass for K0s candidates"};
  o2::framework::Configurable<float> cutMassK0sMax{"cutMassK0sMax", 0.509, "maximum mass for K0s candidates"};
  o2::framework::Configurable<float> cutMassLambdaMin{"cutMassLambdaMin", 1.11, "minimum mass for Lambda candidates"};
  o2::framework::Configurable<float> cutMassLambdaMax{"cutMassLambdaMax", 1.12, "maximum mass for Lambda candidates"};
};

/// Helper method to add histograms to the registry
/// \param registry is the histogram registry
template <bool DoMc, DMesonType DType>
void addHistograms(o2::framework::HistogramRegistry& registry)
{
  constexpr uint8_t NumBinsEvents = EventType::NEventType;
  std::string labels[NumBinsEvents];
  labels[EventType::Processed] = "processed";
  labels[EventType::NoDV0Selected] = "without DV0 pairs";
  labels[EventType::DV0Selected] = "with DV0 pairs";
  const o2::framework::AxisSpec axisEvents = {NumBinsEvents, 0.5, NumBinsEvents + 0.5, ""};
  registry.add("hEvents", "Events;;entries", o2::framework::HistType::kTH1D, {axisEvents});
  for (auto iBin{0u}; iBin < NumBinsEvents; iBin++) {
    registry.get<TH1>(HIST("hEvents"))->GetXaxis()->SetBinLabel(iBin + 1, labels[iBin].data());
  }

  const o2::framework::AxisSpec axisPt{50, 0.f, 50.f, "#it{p}_{T} (GeV/#it{c})"};
  const o2::framework::AxisSpec axisP{100, 0.f, 10.f, "#it{p} (GeV/#it{c})"};
  const o2::framework::AxisSpec axisDeDx{500, 0.f, 1000.f, ""};
  const o2::framework::AxisSpec axisMassLambda{100, 1.05f, 1.35f, "inv. mass (GeV/#it{c}^{2})"};
  const o2::framework::AxisSpec axisMassKzero{100, 0.35f, 0.65f, "inv. mass (GeV/#it{c}^{2})"};
  const o2::framework::AxisSpec axisDeltaMassToK{500, 0.49, 1.49, "inv. mass (GeV/#it{c}^{2})"};
  const o2::framework::AxisSpec axisDeltaMassToPi{500, 0.13, 1.13, "inv. mass (GeV/#it{c}^{2})"};
  const o2::framework::AxisSpec axisDeltaMassToPr{500, 0.93, 1.93, "inv. mass (GeV/#it{c}^{2})"};
  const o2::framework::AxisSpec axisDeltaMassToLambda{500, 1.05, 2.05, "inv. mass (GeV/#it{c}^{2})"};
  const o2::framework::AxisSpec axisDeltaMassToGamma{500, 0., 0.25, "inv. mass (GeV/#it{c}^{2})"};
  const o2::framework::AxisSpec axisMassDsj{400, 0.49f, 0.89f, ""}; // Ds1 and Ds2Star legacy
  const o2::framework::AxisSpec axisAlpha{100, -1.f, 1.f};
  const o2::framework::AxisSpec axisQt{100, 0.f, 0.25f};
  const o2::framework::AxisSpec axisRadius{450, 0.f, 90.f};

  registry.add("hMassVsPtK0s", "K0^{s} candidates;#it{p}_{T} (GeV/#it{c});inv. mass (#pi^{#plus}#pi^{#minus}) (GeV/#it{c}^{2});entries", {o2::framework::HistType::kTH2D, {axisPt, axisMassKzero}});
  registry.add("hMassVsPtLambda", "Lambda candidates;#it{p}_{T} (GeV/#it{c});inv. mass (p #pi^{#minus}) (GeV/#it{c}^{2});entries", {o2::framework::HistType::kTH2D, {axisPt, axisMassLambda}});
  registry.add("hAP", "Aremnteros-Podolanski plot for V0 candidates;#it{#alpha};#it{q}_{T} (GeV/#it{c});entries", {o2::framework::HistType::kTH2D, {axisAlpha, axisQt}});
  registry.add("hRadius", "Radius of V0 candidates;#it{R} (cm);entries", {o2::framework::HistType::kTH1D, {axisRadius}});
  registry.add("hdEdxVsP", "Tracks;#it{p} (GeV/#it{c});d#it{E}/d#it{x};entries", {o2::framework::HistType::kTH2D, {axisP, axisDeDx}});

  if constexpr (DType == DMesonType::D0) {
    const o2::framework::AxisSpec axisMassD0{200, 1.7f, 2.1f, "inv. mass (GeV/#it{c}^{2})"};
    registry.add("hMassVsPtD0All", "D0 candidates (all, regardless the pairing with V0s);#it{p}_{T} (GeV/#it{c});inv. mass (GeV/#it{c}^{2});entries", {o2::framework::HistType::kTH2D, {axisPt, axisMassD0}});
    registry.add("hMassVsPtD0BarAll", "D0bar candidates (all, regardless the pairing with V0s);#it{p}_{T} (GeV/#it{c});inv. mass (GeV/#it{c}^{2});entries", {o2::framework::HistType::kTH2D, {axisPt, axisMassD0}});
    registry.add("hMassVsPtD0Paired", "D0 candidates (paired with V0s);#it{p}_{T} (GeV/#it{c});inv. mass (GeV/#it{c}^{2});entries", {o2::framework::HistType::kTH2D, {axisPt, axisMassD0}});
    registry.add("hMassVsPtD0BarPaired", "D0 candidates (paired with V0s);#it{p}_{T} (GeV/#it{c});inv. mass (GeV/#it{c}^{2});entries", {o2::framework::HistType::kTH2D, {axisPt, axisMassD0}});
    registry.add("hMassD0Pi", "D0Pi candidates; m_{D^{0}#pi^{+}} - m_{D^{0}} (GeV/#it{c}^{2});entries", {o2::framework::HistType::kTH2D, {axisPt, axisDeltaMassToPi}});
    registry.add("hMassD0K", "D0Kplus candidates; m_{D^{0}K^{+}} - m_{D^{0}} (GeV/#it{c}^{2});entries", {o2::framework::HistType::kTH2D, {axisPt, axisDeltaMassToK}});
    registry.add("hMassD0Proton", "D0Proton candidates; m_{D^{0}p} - m_{D^{0}} (GeV/#it{c}^{2});entries", {o2::framework::HistType::kTH2D, {axisPt, axisDeltaMassToPr}});
    registry.add("hMassD0Lambda", "D0Lambda candidates; m_{D^{0}#Lambda} - m_{D^{0}} (GeV/#it{c}^{2});entries", {o2::framework::HistType::kTH2D, {axisPt, axisDeltaMassToLambda}});
    registry.add("hMassD0Gamma", "D0Gamma candidates; m_{D^{0}#gamma} - m_{D^{0}} (GeV/#it{c}^{2});entries", {o2::framework::HistType::kTH2D, {axisPt, axisDeltaMassToGamma}});
  } else if constexpr (DType == DMesonType::Dplus) {
    const o2::framework::AxisSpec axisMassDplus{200, 1.7f, 2.1f, "inv. mass (GeV/#it{c}^{2})"};
    registry.add("hMassVsPtDplusAll", "Dplus candidates (all, regardless the pairing with V0s);#it{p}_{T} (GeV/#it{c});inv. mass (GeV/#it{c}^{2});entries", {o2::framework::HistType::kTH2D, {axisPt, axisMassDplus}});
    registry.add("hMassVsPtDplusPaired", "Dplus candidates (paired with V0s);#it{p}_{T} (GeV/#it{c});inv. mass (GeV/#it{c}^{2});entries", {o2::framework::HistType::kTH2D, {axisPt, axisMassDplus}});
    registry.add("hMassDplusK0s", "DplusK0s candidates; m_{D^{+}K^{0}_{S}} (GeV/#it{c}^{2});entries", {o2::framework::HistType::kTH2D, {axisPt, axisDeltaMassToK}});
    registry.add("hMassDplusPi", "DplusPi candidates; m_{D^{+}#pi^{-}} (GeV/#it{c}^{2});entries", {o2::framework::HistType::kTH2D, {axisPt, axisDeltaMassToPi}});
    registry.add("hMassDplusK", "DplusK candidates; m_{D^{+}#pi^{-}} (GeV/#it{c}^{2});entries", {o2::framework::HistType::kTH2D, {axisPt, axisDeltaMassToK}});
    registry.add("hMassDplusProton", "DplusProton candidates; m_{D^{+}p} (GeV/#it{c}^{2});entries", {o2::framework::HistType::kTH2D, {axisPt, axisDeltaMassToPr}});
    registry.add("hMassDplusLambda", "DplusLambda candidates; m_{D^{+}#Lambda} (GeV/#it{c}^{2});entries", {o2::framework::HistType::kTH2D, {axisPt, axisDeltaMassToLambda}});
  } else if constexpr (DType == DMesonType::Dstar) {
    const o2::framework::AxisSpec axisMassDstar{200, 0.139f, 0.179f, "delta inv. mass (GeV/#it{c}^{2})"}; // o2-linter: disable=pdg/explicit-mass (false positive)
    registry.add("hMassVsPtDstarAll", "Dstar candidates (all, regardless the pairing with V0s);#it{p}_{T} (GeV/#it{c});inv. mass (GeV/#it{c}^{2});entries", {o2::framework::HistType::kTH2D, {axisPt, axisMassDstar}});
    registry.add("hMassVsPtDstarPaired", "Dstar candidates (paired with V0s);#it{p}_{T} (GeV/#it{c});inv. mass (GeV/#it{c}^{2});entries", {o2::framework::HistType::kTH2D, {axisPt, axisMassDstar}});
    registry.add("hMassDstarPi", "DstarPi candidates; m_{D^{*+}#pi^{-}} (GeV/#it{c}^{2});entries", {o2::framework::HistType::kTH2D, {axisPt, axisDeltaMassToPi}});
    registry.add("hMassDstarK", "DstarK candidates; m_{D^{*+}#pi^{-}} (GeV/#it{c}^{2});entries", {o2::framework::HistType::kTH2D, {axisPt, axisDeltaMassToK}});
    registry.add("hMassDstarProton", "DstarProton candidates; m_{D^{*}p} (GeV/#it{c}^{2});entries", {o2::framework::HistType::kTH2D, {axisPt, axisDeltaMassToPr}});
    registry.add("hMassDstarK0s", "DstarK0s candidates; m_{D^{*}K^{0}_{S}} (GeV/#it{c}^{2});entries", {o2::framework::HistType::kTH2D, {axisPt, axisDeltaMassToK}});
    registry.add("hMassDstarLambda", "DstarLambda candidates; m_{D^{*}#Lambda} (GeV/#it{c}^{2});entries", {o2::framework::HistType::kTH2D, {axisPt, axisDeltaMassToLambda}});
  }

  if constexpr (DoMc) {
    // MC Rec
    int const nChannels = hf_decay::hf_cand_reso::DecayChannelMain::NChannelsMain;
    registry.add("hMCRecCounter", "Number of Reconstructed MC Matched candidates per channel", {o2::framework::HistType::kTH1D, {{2 * nChannels + 1, -(nChannels + 0.5), nChannels + 0.5}}});
    registry.add("hMCRecDebug", "Debug of MC Reco", {o2::framework::HistType::kTH1D, {{551, -0.5, 550.5}}});
    registry.add("hMCRecOrigin", "Origin of Matched particles", {o2::framework::HistType::kTH1D, {{3, -0.5, 2.5}}});
    registry.add("hMCRecMassGen", "Generated inv. mass of resoncances", {o2::framework::HistType::kTH1D, {{2000, 1.8, 3.8}}});
    registry.add("hMCRecCharmDau", "Charm daughter flag", {o2::framework::HistType::kTH1D, {{57, -28.5, 28.5}}});
    // MC Gen
    registry.add("hMCGenCounter", "Number of Generated particles; Decay Channel Flag; pT (GeV/#it{c})", {o2::framework::HistType::kTH2D, {{17, -8.5, 8.5}, {100, 0, 50}}});
    registry.add("hMCGenOrigin", "Origin of Generated particles", {o2::framework::HistType::kTH1D, {{3, -0.5, 2.5}}});
  }
}

/// Basic track quality selections for V0 daughters
/// \param Tr is a track
/// \param dDaughtersIds are the IDs of the D meson daughter tracks
/// \param cfgV0Cuts are the cuts to be applied to the V0
/// \param rejectPairsWithCommonDaughter is a flag to activate rejection of pairs sharing a daughter track
template <typename Tr, typename Cuts>
bool selectV0Daughter(Tr const& track, const std::array<int, 3>& dDaughtersIds, const Cuts& cfgV0Cuts, bool rejectPairsWithCommonDaughter)
{
  // acceptance selection
  if (std::abs(track.eta()) > cfgV0Cuts.etaMaxDau.value) {
    return false;
  }
  // Tpc Refit
  if (!(track.hasTPC())) {
    return false;
  }
  // track quality selection
  if (track.itsNCls() < cfgV0Cuts.trackNclusItsCut.value ||
      track.tpcNClsFound() < cfgV0Cuts.trackNCrossedRowsTpc.value ||
      track.tpcNClsCrossedRows() < cfgV0Cuts.trackNCrossedRowsTpc.value ||
      track.tpcNClsCrossedRows() < cfgV0Cuts.trackFracMaxindableTpcCls.value * track.tpcNClsFindable() ||
      track.tpcNClsShared() > cfgV0Cuts.trackNsharedClusTpc.value) {
    return false;
  }
  // rejection of tracks that share a daughter with the D meson
  if (rejectPairsWithCommonDaughter && std::find(dDaughtersIds.begin(), dDaughtersIds.end(), track.globalIndex()) != dDaughtersIds.end()) {
    return false;
  }
  return true;
}

/// Utility to find which v0 daughter carries the largest fraction of the mother longitudinal momentum
/// \param momV0 is the momentum of the V0
/// \param momDau0 is the momentum of first daughter
/// \param momDau1 is the momentum of second daughter
/// \return alphaAP
float alphaAP(std::array<float, 3> const& momV0, std::array<float, 3> const& momDau0, std::array<float, 3> const& momDau1)
{
  float const momTot = std::hypot(momV0[0], momV0[1], momV0[2]);
  float const lQlPos = (momDau0[0] * momV0[0] + momDau0[1] * momV0[1] + momDau0[2] * momV0[2]) / momTot;
  float const lQlNeg = (momDau1[0] * momV0[0] + momDau1[1] * momV0[1] + momDau1[2] * momV0[2]) / momTot;
  return (lQlPos - lQlNeg) / (lQlPos + lQlNeg);
}

/// Utility to compute qT
/// \param momDau0 is the momentum of first daughter
/// \param momDau1 is the momentum of second daughter
/// \return qtAP
//_______________________________________________________________________
inline float qtAP(std::array<float, 3> const& momDau0, std::array<float, 3> const& momDau1)
{
  float momTot = RecoDecay::p2(momDau0[0] + momDau1[0], momDau0[1] + momDau1[1], momDau0[2] + momDau1[2]);
  float dp = RecoDecay::dotProd(std::array{momDau1[0], momDau1[1], momDau1[2]}, std::array{momDau0[0] + momDau1[0], momDau0[1] + momDau1[1], momDau0[2] + momDau1[2]});
  return std::sqrt(RecoDecay::p2(momDau1[0], momDau1[1], momDau1[2]) - dp * dp / momTot); // qt of v0
}

/// Utility to find DCA of V0 to Primary vertex
/// \param x is the x-coordinate
/// \param y is the y-coordinate
/// \param z is the z-coordinate
/// \param px is the x-component of the momentum
/// \param py is the y-component of the momentum
/// \param pz is the z-component of the momentum
/// \param pvX is the x-coordinate of the PV
/// \param pvY is the y-coordinate of the PV
/// \param pvZ is the z-coordinate of the PV
/// \return the DCA
float calculateDCAStraightToPV(float x, float y, float z, float px, float py, float pz, float pvX, float pvY, float pvZ)
{
  return std::hypot((pvY - y) * pz - (pvZ - z) * py, (pvX - x) * pz - (pvZ - z) * px, (pvX - x) * py - (pvY - y) * px) / (px * px + py * py + pz * pz);
}

/// Basic selection of V0 candidates
/// \param collision is the current collision
/// \param dauTracks are the v0 daughter tracks
/// \param dDaughtersIds are the IDs of the D meson daughter tracks
/// \param fitter is the DCAFitter object
/// \param cfgV0Cuts are the cuts to be applied to the V0
/// \param v0 is the V0 candidate
/// \param rejectPairsWithCommonDaughter is a flag to activate rejection of pairs sharing a daughter track
/// \return a bitmap with mass hypotesis if passes all cuts
template <typename Coll, typename Tr, typename Cuts>
bool buildAndSelectV0(const Coll& collision, const std::array<int, 3>& dDaughtersIds, const std::array<Tr, 2>& dauTracks, const Cuts& cfgV0Cuts, o2::vertexing::DCAFitterN<2>& fitter, HfResoCandidateV0& v0, bool rejectPairsWithCommonDaughter)
{
  const auto& trackPos = dauTracks[0];
  const auto& trackNeg = dauTracks[1];
  // single-tracks selection
  if (!selectV0Daughter(trackPos, dDaughtersIds, cfgV0Cuts, rejectPairsWithCommonDaughter) || !selectV0Daughter(trackNeg, dDaughtersIds, cfgV0Cuts, rejectPairsWithCommonDaughter)) {
    return false;
  }
  // daughters DCA to V0's collision primary vertex
  std::array<float, 2> dcaInfo{};
  auto trackPosPar = getTrackPar(trackPos);
  o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackPosPar, 2.f, fitter.getMatCorrType(), &dcaInfo);
  auto trackPosDcaXY = dcaInfo[0];
  auto trackNegPar = getTrackPar(trackNeg);
  o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackNegPar, 2.f, fitter.getMatCorrType(), &dcaInfo);
  auto trackNegDcaXY = dcaInfo[0];
  if (std::fabs(trackPosDcaXY) < cfgV0Cuts.dcaMaxDauToPv.value || std::fabs(trackNegDcaXY) < cfgV0Cuts.dcaMaxDauToPv.value) {
    return false;
  }
  // vertex reconstruction
  auto trackPosCov = getTrackParCov(trackPos);
  auto trackNegCov = getTrackParCov(trackNeg);
  int nCand = 0;
  try {
    nCand = fitter.process(trackPosCov, trackNegCov);
  } catch (...) {
    return false;
  }
  if (nCand == 0) {
    return false;
  }
  // compute candidate momentum from tracks propagated to decay vertex
  auto& trackPosProp = fitter.getTrack(0);
  auto& trackNegProp = fitter.getTrack(1);
  trackPosProp.getPxPyPzGlo(v0.momPos);
  trackNegProp.getPxPyPzGlo(v0.momNeg);

  v0.mom = RecoDecay::pVec(v0.momPos, v0.momNeg);

  v0.pT = std::hypot(v0.mom[0], v0.mom[1]);
  // topological selections:
  // v0 eta
  v0.eta = RecoDecay::eta(v0.mom);
  if (std::abs(v0.eta) > cfgV0Cuts.etaMax.value) {
    return false;
  }
  // daughters DCA
  v0.dcaDau = std::sqrt(fitter.getChi2AtPCACandidate());
  if (v0.dcaDau > cfgV0Cuts.dcaDau.value) {
    return false;
  }
  // v0 radius
  const auto& vtx = fitter.getPCACandidate();
  v0.radius = std::hypot(vtx[0], vtx[1]);
  if (v0.radius < cfgV0Cuts.radiusMin.value) {
    return false;
  }
  std::copy(vtx.begin(), vtx.end(), v0.pos.begin());

  // v0 DCA to primary vertex
  v0.dcaV0ToPv = calculateDCAStraightToPV(
    vtx[0], vtx[1], vtx[2],
    v0.momPos[0] + v0.momNeg[0],
    v0.momPos[1] + v0.momNeg[1],
    v0.momPos[2] + v0.momNeg[2],
    collision.posX(), collision.posY(), collision.posZ());
  if (std::abs(v0.dcaV0ToPv) > cfgV0Cuts.dcaPv.value) {
    return false;
  }
  // v0 cosine of pointing angle
  std::array<float, 3> const primVtx = {collision.posX(), collision.posY(), collision.posZ()};
  v0.cosPA = RecoDecay::cpa(primVtx, vtx, v0.mom);
  if (v0.cosPA < cfgV0Cuts.cosPa.value) {
    return false;
  }
  // distinguish between K0s, and Lambda hypotesys
  v0.v0Type = {BIT(BachelorType::K0s) | BIT(BachelorType::Lambda) | BIT(BachelorType::AntiLambda)};
  // for lambda hypotesys define if its lambda or anti-lambda
  v0.alpha = alphaAP(v0.mom, v0.momPos, v0.momNeg);
  v0.qt = qtAP(v0.momPos, v0.momNeg);
  bool const matter = v0.alpha > 0;
  CLRBIT(v0.v0Type, matter ? BachelorType::AntiLambda : BachelorType::Lambda);
  auto massPos = matter ? o2::constants::physics::MassProton : o2::constants::physics::MassPionCharged;
  auto massNeg = matter ? o2::constants::physics::MassPionCharged : o2::constants::physics::MassProton;
  // mass hypotesis
  v0.mLambda = RecoDecay::m(std::array{v0.momPos, v0.momNeg}, std::array{massPos, massNeg});
  v0.mK0Short = RecoDecay::m(std::array{v0.momPos, v0.momNeg}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassPionCharged});
  if (std::fabs(v0.mK0Short - o2::constants::physics::MassK0) > cfgV0Cuts.deltaMassK0s.value) {
    CLRBIT(v0.v0Type, BachelorType::K0s);
  }
  if (std::fabs(v0.mLambda - o2::constants::physics::MassLambda0) > cfgV0Cuts.deltaMassLambda.value) {
    CLRBIT(v0.v0Type, BachelorType::Lambda);
    CLRBIT(v0.v0Type, BachelorType::AntiLambda);
  }
  // PID
  if (TESTBIT(v0.v0Type, BachelorType::K0s)) {
    if ((trackPos.hasTPC() && std::fabs(trackPos.tpcNSigmaPi()) > cfgV0Cuts.nSigmaTpc.value) ||
        (trackNeg.hasTPC() && std::fabs(trackNeg.tpcNSigmaPi()) > cfgV0Cuts.nSigmaTpc.value)) {
      CLRBIT(v0.v0Type, BachelorType::K0s);
    }
  }
  if (TESTBIT(v0.v0Type, BachelorType::Lambda)) {
    if ((trackPos.hasTPC() && std::fabs(trackPos.tpcNSigmaPr()) > cfgV0Cuts.nSigmaTpc.value) ||
        (trackPos.hasTOF() && std::fabs(trackPos.tofNSigmaPr()) > cfgV0Cuts.nSigmaTofPr.value) ||
        (trackNeg.hasTPC() && std::fabs(trackNeg.tpcNSigmaPi()) > cfgV0Cuts.nSigmaTpc.value)) {
      CLRBIT(v0.v0Type, BachelorType::Lambda);
    }
  }
  if (TESTBIT(v0.v0Type, BachelorType::AntiLambda)) {
    if ((trackPos.hasTPC() && std::fabs(trackPos.tpcNSigmaPi()) > cfgV0Cuts.nSigmaTpc.value) ||
        (trackNeg.hasTPC() && std::fabs(trackNeg.tpcNSigmaPr()) > cfgV0Cuts.nSigmaTpc.value) ||
        (trackNeg.hasTOF() && std::fabs(trackNeg.tofNSigmaPr()) > cfgV0Cuts.nSigmaTofPr.value)) {
      CLRBIT(v0.v0Type, BachelorType::AntiLambda);
    }
  }
  if (v0.v0Type == 0) {
    return false;
  }
  return true;
}

/// Basic selection of V0 candidates
/// \param collision is the current collision
/// \param dauTracks are the v0 daughter tracks
/// \param dDaughtersIds are the IDs of the D meson daughter tracks
/// \param cfgV0Cuts are the cuts to be applied to the V0
/// \param v0 is the V0 candidate
/// \param matCorr is the material correction type to be used in the track propagation
/// \param bz is the magnetic field
/// \param vDriftMgr is the TPC velocity drift manager
/// \param rejectPairsWithCommonDaughter is a flag to activate rejection of pairs sharing a daughter track
/// \return a bitmap with mass hypotesis if passes all cuts
template <class BCs, class Colls, typename Coll, typename Tr, typename Cuts>
bool buildAndSelectGamma(const Coll& collision, const std::array<int, 3>& dDaughtersIds, const std::array<Tr, 2>& dauTracks, const Cuts& cfgGammaCuts, HfResoCandidateV0& v0, o2::base::Propagator::MatCorrType const& matCorr, o2::aod::common::TPCVDriftManager* vDriftMgr, bool rejectPairsWithCommonDaughter)
{
  const auto& trackPos = dauTracks[0];
  const auto& trackNeg = dauTracks[1];
  if (trackPos.sign() * trackNeg.sign() > 0) { // reject same sign pair
    return false;
  }
  if (trackPos.globalIndex() == trackNeg.globalIndex()) {
    return false;
  }
  if (o2::pwgem::photonmeson::isITSonlyTrack(trackPos) && !trackNeg.hasITS()) {
    return false;
  }
  if (o2::pwgem::photonmeson::isITSonlyTrack(trackNeg) && !trackPos.hasITS()) {
    return false;
  }

  // single-tracks selection
  if (!selectV0Daughter(trackPos, dDaughtersIds, cfgGammaCuts, rejectPairsWithCommonDaughter) || !selectV0Daughter(trackNeg, dDaughtersIds, cfgGammaCuts, rejectPairsWithCommonDaughter)) {
    return false;
  }
  if ((trackPos.hasTPC() && std::abs(trackPos.tpcNSigmaEl()) > cfgGammaCuts.nSigmaTpcEl.value) || (trackNeg.hasTPC() && std::abs(trackNeg.tpcNSigmaEl()) > cfgGammaCuts.nSigmaTpcEl.value)) {
    return false;
  }

  std::array<float, 2> dcaInfo;
  auto trackParPos = getTrackParCov(trackPos);
  if (o2::pwgem::photonmeson::isTPConlyTrack(trackPos) && !vDriftMgr->moveTPCTrack<BCs, Colls>(collision, trackPos, trackParPos)) {
    LOGP(error, "failed correction for positive tpc track");
    return false;
  }
  auto trackParPropPos = trackParPos;
  trackParPropPos.setPID(o2::track::PID::Electron);
  o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParPropPos, 2.f, matCorr, &dcaInfo);
  auto trackPosDcaXY = dcaInfo[0];

  auto trackParNeg = getTrackParCov(trackNeg);
  if (o2::pwgem::photonmeson::isTPConlyTrack(trackNeg) && !vDriftMgr->moveTPCTrack<BCs, Colls>(collision, trackNeg, trackParNeg)) {
    LOGP(error, "failed correction for negative tpc track");
    return false;
  }
  auto trackParPropNeg = trackParNeg;
  trackParPropNeg.setPID(o2::track::PID::Electron);
  o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParPropNeg, 2.f, matCorr, &dcaInfo);
  auto trackNegDcaXY = dcaInfo[0];

  // daughters DCA to V0's collision primary vertex
  if (std::fabs(trackPosDcaXY) < cfgGammaCuts.dcaMaxDauToPv.value || std::fabs(trackNegDcaXY) < cfgGammaCuts.dcaMaxDauToPv.value) {
    return false;
  }

  float gammaVtx[3] = {0.f, 0.f, 0.f};
  Vtx_recalculationParCov(o2::base::Propagator::Instance(), trackParPropPos, trackParPropNeg, gammaVtx, matCorr);
  float radiusXy = std::hypot(gammaVtx[0], gammaVtx[1]);
  const float maxX{83.1f};    // max X for track IU
  const float marginTpc{7.f}; // margin for r cut in cm
  if (radiusXy > maxX + marginTpc) {
    return false;
  }
  if (radiusXy < std::fabs(gammaVtx[2]) * std::tan(2 * std::atan(std::exp(-cfgGammaCuts.etaMax.value))) - marginTpc) {
    return false; // RZ line cut
  }

  // vertex reconstruction
  KFPTrack kfpTrackPos = createKFPTrackFromTrackParCov(trackParPropPos, trackPos.sign(), trackPos.tpcNClsFound(), trackPos.tpcChi2NCl());
  KFPTrack kfpTrackNeg = createKFPTrackFromTrackParCov(trackParPropNeg, trackNeg.sign(), trackNeg.tpcNClsFound(), trackNeg.tpcChi2NCl());
  KFParticle kfPartPos(kfpTrackPos, kPositron);
  KFParticle kfPartNeg(kfpTrackNeg, kElectron);
  const KFParticle* gammaDaughters[2] = {&kfPartPos, &kfPartNeg};

  KFParticle gamma;
  gamma.SetConstructMethod(2);
  gamma.Construct(gammaDaughters, 2);
  KFPVertex kfpVertex = createKFPVertexFromCollision(collision);
  KFParticle KFPV(kfpVertex);

  // Transport the gamma to the recalculated decay vertex
  KFParticle gammaDecayVtx = gamma; // with respect to (0,0,0)
  gammaDecayVtx.TransportToPoint(gammaVtx);
  v0.cosPA = cpaFromKF(gammaDecayVtx, KFPV);
  if (v0.cosPA < cfgGammaCuts.cosPa.value) {
    return false;
  }

  v0.pos = {gammaDecayVtx.GetX(), gammaDecayVtx.GetY(), gammaDecayVtx.GetZ()};
  v0.radius = std::hypot(gammaDecayVtx.GetX(), gammaDecayVtx.GetY());
  if (v0.radius > maxX + marginTpc) {
    return false;
  }
  if (v0.radius < std::fabs(gammaDecayVtx.GetZ()) * std::tan(2 * std::atan(std::exp(-cfgGammaCuts.etaMax.value))) - marginTpc) {
    return false; // RZ line cut
  }
  if (v0.radius < cfgGammaCuts.radiusMin.value || cfgGammaCuts.radiusMax.value < v0.radius) {
    return false;
  }

  const float minRadTpcOnly{16.f};
  if ((!trackNeg.hasITS() && !trackNeg.hasITS()) && v0.radius < minRadTpcOnly) { // TPConly tracks can detect conversion points larger than minRadTpcOnly.
    return false;
  }

  // Apply a topological constraint of the gamma to the PV. Parameters will be given at the primary vertex.
  KFParticle gammaPvConstr = gamma;
  gammaPvConstr.SetProductionVertex(KFPV);
  v0.mom = RecoDecay::pVec(std::array{gammaPvConstr.GetPx(), gammaPvConstr.GetPy(), gammaPvConstr.GetPz()});
  v0.pT = std::hypot(v0.mom[0], v0.mom[1]);
  if (v0.pT < cfgGammaCuts.ptMin.value) {
    return false;
  }
  if (o2::pwgem::photonmeson::isITSonlyTrack(trackNeg) && o2::pwgem::photonmeson::isITSonlyTrack(trackPos) && v0.pT > cfgGammaCuts.ptMaxItsOnly.value) {
    return false;
  }
  v0.eta = RecoDecay::eta(v0.mom);
  if (std::abs(v0.eta) > cfgGammaCuts.etaMax.value) {
    return false;
  }

  KFParticle kfPartDecayVtxPos = kfPartPos;     // Don't set Primary Vertex
  KFParticle kfPartDecayVtxNeg = kfPartNeg;     // Don't set Primary Vertex
  kfPartDecayVtxPos.TransportToPoint(gammaVtx); // Don't set Primary Vertex
  kfPartDecayVtxNeg.TransportToPoint(gammaVtx); // Don't set Primary Vertex
  v0.dcaDau = kfPartDecayVtxPos.GetDistanceFromParticle(kfPartDecayVtxNeg);
  v0.momPos = RecoDecay::pVec(std::array{kfPartDecayVtxPos.GetPx(), kfPartDecayVtxPos.GetPy(), kfPartDecayVtxPos.GetPz()});
  v0.momNeg = RecoDecay::pVec(std::array{kfPartDecayVtxNeg.GetPx(), kfPartDecayVtxNeg.GetPy(), kfPartDecayVtxNeg.GetPz()});
  float ptItsOnlyMax{0.15f};
  if (o2::pwgem::photonmeson::isITSonlyTrack(trackPos) && std::hypot(v0.momPos[0], v0.momPos[1]) > ptItsOnlyMax) {
    return false;
  }
  if (o2::pwgem::photonmeson::isITSonlyTrack(trackNeg) && std::hypot(v0.momNeg[0], v0.momNeg[1]) > ptItsOnlyMax) {
    return false;
  }

  const float maxRItsMft{66.f};
  if (!trackNeg.hasITS() && !trackPos.hasITS()) { // V0s with TPConly-TPConly
    if (maxRItsMft < v0.radius && v0.radius < maxX + marginTpc) {
      if (v0.dcaDau > cfgGammaCuts.dcaDauTpcInner.value) {
        return false;
      }
    } else {
      if (v0.dcaDau > cfgGammaCuts.dcaDauTpc.value) {
        return false;
      }
    }
  } else { // V0s with ITS hits
    if (v0.radius < minRadTpcOnly) {
      if (v0.dcaDau > cfgGammaCuts.dcaDauItsIb.value) {
        return false;
      }
    } else {
      if (v0.dcaDau > cfgGammaCuts.dcaDauIts.value) {
        return false;
      }
    }
  }

  // v0 DCA to primary vertex
  v0.dcaV0ToPv = calculateDCAStraightToPV(
    v0.pos[0], v0.pos[1], v0.pos[2],
    v0.momPos[0] + v0.momNeg[0],
    v0.momPos[1] + v0.momNeg[1],
    v0.momPos[2] + v0.momNeg[2],
    collision.posX(), collision.posY(), collision.posZ());
  if (std::abs(v0.dcaV0ToPv) > cfgGammaCuts.dcaPv.value) {
    return false;
  }

  // distinguish V0 hypotheses
  v0.alpha = alphaAP(v0.mom, v0.momPos, v0.momNeg);
  v0.qt = qtAP(v0.momPos, v0.momNeg);
  ;
  if (!checkAP(v0.alpha, v0.qt, cfgGammaCuts.alphaApMax.value, cfgGammaCuts.qtApMax.value)) { // store only photon conversions
    return false;
  }
  v0.v0Type = BIT(BachelorType::Gamma);
  return true;
}

/// Basic selection of tracks
/// \param track is the track
/// \param dDaughtersIds are the IDs of the D meson daughter tracks
/// \param cfgV0Cuts are the cuts to be applied to the track
/// \param rejectPairsWithCommonDaughter is a flag to activate rejection of pairs sharing a daughter track
/// \return true if passes all cuts
template <typename Tr, typename Cuts>
bool isTrackSelected(const Tr& track, const std::array<int, 3>& dDaughtersIds, const Cuts& cfgSingleTrackCuts, bool rejectPairsWithCommonDaughter)
{
  if (rejectPairsWithCommonDaughter && std::find(dDaughtersIds.begin(), dDaughtersIds.end(), track.globalIndex()) != dDaughtersIds.end()) {
    return false;
  }
  switch (cfgSingleTrackCuts.setTrackSelections.value) {
    case TrackSel::GlobalTrack:
      if (!track.isGlobalTrack()) {
        return false;
      }
      break;
    case TrackSel::GlobalTrackWoDca:
      if (!track.isGlobalTrackWoDCA()) {
        return false;
      }
      break;
    case TrackSel::QualityTrackITS:
      if (!track.isQualityTrackITS()) {
        return false;
      }
      break;
  }
  if (track.pt() < cfgSingleTrackCuts.minPt.value) {
    return false;
  }
  if (std::abs(track.eta()) > cfgSingleTrackCuts.maxEta.value) {
    return false;
  }
  if (!track.hasTPC()) {
    return false;
  }
  bool const isPion = std::abs(track.tpcNSigmaPi()) < cfgSingleTrackCuts.maxNsigmaTpcPi.value;
  bool const isKaon = std::abs(track.tpcNSigmaKa()) < cfgSingleTrackCuts.maxNsigmaTpcKa.value;
  bool const isProton = std::abs(track.tpcNSigmaPr()) < cfgSingleTrackCuts.maxNsigmaTpcPr.value;
  return (isPion || isKaon || isProton); // we keep the track if is it compatible with at least one of the PID hypotheses selected
}

/// Matching of V0 candidates to MC truth
/// \param particlesMc is the table of MC particles
/// \param arrDaughtersV0 is the array of V0 daughter tracks
/// \return the MC matching flag for the V0
template <typename PParticles, typename TrIU>
int8_t getMatchingFlagV0(PParticles const& particlesMc, const std::array<TrIU, 2>& arrDaughtersV0)
{
  int8_t signV0{0};
  int indexRec{-1};
  int flagV0{0};
  indexRec = RecoDecay::getMatchedMCRec<false, true, false, true, true>(particlesMc, arrDaughtersV0, kK0, std::array{+kPiPlus, -kPiPlus}, true, &signV0, 2);
  if (indexRec > -1) {
    flagV0 = hf_decay::hf_cand_reso::PartialMatchMc::K0Matched;
  } else {
    indexRec = RecoDecay::getMatchedMCRec<false, true, false, true, true>(particlesMc, arrDaughtersV0, kLambda0, std::array{+kProton, -kPiPlus}, true, &signV0, 2);
    if (indexRec > -1) {
      flagV0 = signV0 * hf_decay::hf_cand_reso::PartialMatchMc::LambdaMatched;
    }
  }
  return flagV0; // Placeholder, should return the actual flag based on matching logic
}

/// Matching of V0 candidates to MC truth
/// \param bachTrack is the track
/// \return the MC matching flag for the track
template <typename Tr>
int8_t getMatchingFlagTrack(Tr const& bachTrack)
{
  auto particle = bachTrack.mcParticle();
  auto pdgCode = std::abs(particle.pdgCode());
  if (pdgCode == kPiPlus) {
    return hf_decay::hf_cand_reso::PartialMatchMc::PionMatched;
  }
  if (pdgCode == kKPlus) {
    return hf_decay::hf_cand_reso::PartialMatchMc::KaonMatched;
  }
  if (pdgCode == kProton) {
    return hf_decay::hf_cand_reso::PartialMatchMc::ProtonMatched;
  }
  return 0;
}

/// Matching of V0 candidates to MC truth
/// \param particlesMc is the table of MC particles
/// \param indexRec is the index of the MC particle associated to the reconstructed canditate
/// \param pdg is the O2DatabasePDG service
/// \return the generated invariant mass
template <typename PParticles>
float computeInvMassGen(PParticles const& particlesMc, int indexRec, o2::framework::Service<o2::framework::O2DatabasePDG> const& pdg)
{
  if (indexRec < 0) {
    return -1.f;
  }
  auto particleReso = particlesMc.iteratorAt(indexRec);
  auto dau1 = particlesMc.iteratorAt(particleReso.daughtersIds().front());
  auto dau2 = particlesMc.iteratorAt(particleReso.daughtersIds().back());
  std::array<std::array<double, 3>, 2> pArr = {{{dau1.px(), dau1.py(), dau1.pz()}, {dau2.px(), dau2.py(), dau2.pz()}}};
  std::array<double, 2> mArr = {pdg->Mass(dau1.pdgCode()), pdg->Mass(dau2.pdgCode())};
  return static_cast<float>(RecoDecay::m(pArr, mArr));
}

/// Function for filling MC reco information of DV0 candidates in the tables
/// \tparam dType is the D meson type (Dstar, Dplus or D0)
/// \param particlesMc is the table with MC particles
/// \param candCharmBach is the D meson candidate
/// \param bachelorV0 is the V0 candidate
/// \param tracks is the table with tracks
/// \param indexHfCandCharm is the index of the charm-hadron bachelor in the reduced table
/// \param indexCandV0TrBach is the index of the v0 bachelor in the reduced table
/// \param pdg is the O2DatabasePDG service
/// \param registry is the histogram registry
/// \param rowMcRecReduced is the table to be filled
template <DMesonType DType, typename PParticles, typename CCand, typename BBachV0, typename Tr, typename Table>
void fillMcRecoInfoDV0(PParticles const& particlesMc,
                       CCand const& candCharmBach,
                       BBachV0 const& bachelorV0,
                       Tr const& tracks,
                       int64_t& indexHfCandCharm,
                       int64_t& indexCandV0Bach,
                       o2::framework::Service<o2::framework::O2DatabasePDG> const& pdg,
                       o2::framework::HistogramRegistry& registry,
                       Table& rowMcRecReduced)
{
  std::vector<typename Tr::iterator> vecDaughtersReso{};
  int8_t sign{0}, nKinkedTracks{0}, origin{0}, flagCharmBach{0}, flagCharmBachInterm{0}, flagV0{0}, flagReso{0};
  int indexRec{-1}, debugMcRec{0};
  float ptGen{-1.f}, invMassGen{-1.f};
  if constexpr (DType == DMesonType::Dstar) {
    vecDaughtersReso.push_back(tracks.rawIteratorAt(candCharmBach.prong0Id()));
    vecDaughtersReso.push_back(tracks.rawIteratorAt(candCharmBach.prong1Id()));
    vecDaughtersReso.push_back(tracks.rawIteratorAt(candCharmBach.prongPiId()));
    // Check if D* is matched
    flagCharmBach = candCharmBach.flagMcMatchRec();
    if (flagCharmBach != 0) {
      SETBIT(debugMcRec, hf_decay::hf_cand_reso::PartialMatchMc::DstarMatched);
      origin = candCharmBach.originMcRec();
    }
    // Check if D0 is matched
    flagCharmBachInterm = candCharmBach.flagMcMatchRecD0();
    if (flagCharmBachInterm != 0) {
      SETBIT(debugMcRec, hf_decay::hf_cand_reso::PartialMatchMc::D0Matched);
    }
    // Check if V0 is matched
    vecDaughtersReso.push_back(tracks.rawIteratorAt(bachelorV0.posTrackId()));
    vecDaughtersReso.push_back(tracks.rawIteratorAt(bachelorV0.negTrackId()));
    flagV0 = getMatchingFlagV0(particlesMc, std::array{vecDaughtersReso[3], vecDaughtersReso[4]});
    if (flagV0 != 0) {
      SETBIT(debugMcRec, std::abs(flagV0));
    }
    // If both D* and K0s are matched, try to match resonance
    if (flagCharmBach != 0 && flagV0 == hf_decay::hf_cand_reso::PartialMatchMc::K0Matched) {
      std::array<int, 5> const pdgCodesDaughters = {+kPiPlus, -kKPlus, +kPiPlus, +kPiPlus, -kPiPlus};
      auto arrDaughtersReso = std::array{vecDaughtersReso[0], vecDaughtersReso[1], vecDaughtersReso[2], vecDaughtersReso[3], vecDaughtersReso[4]};
      for (const auto& [decayChannelFlag, pdgCodeReso] : hf_decay::hf_cand_reso::particlesToDstarK0s) {
        indexRec = RecoDecay::getMatchedMCRec<false, true, false, true, true>(particlesMc, arrDaughtersReso, pdgCodeReso, pdgCodesDaughters, true, &sign, 3, &nKinkedTracks);
        if (indexRec > -1) {
          flagReso = sign * decayChannelFlag;
          break;
        }
      }
    } else if (flagCharmBachInterm != 0 && flagV0 == hf_decay::hf_cand_reso::PartialMatchMc::K0Matched) {
      std::array<int, 4> const pdgCodesDaughters = {+kPiPlus, -kKPlus, +kPiPlus, -kPiPlus};
      auto arrDaughtersReso = std::array{vecDaughtersReso[0], vecDaughtersReso[1], vecDaughtersReso[3], vecDaughtersReso[4]};
      // Peaking background of D0K0s <- Ds* with spurious soft pion
      for (const auto& [decayChannelFlag, pdgCodeReso] : hf_decay::hf_cand_reso::particlesToDstarK0s) {
        indexRec = RecoDecay::getMatchedMCRec<false, true, true, true, true>(particlesMc, arrDaughtersReso, pdgCodeReso, pdgCodesDaughters, true, &sign, 3, &nKinkedTracks);
        if (indexRec > -1) {
          flagReso = sign * decayChannelFlag;
          SETBIT(debugMcRec, hf_decay::hf_cand_reso::PartialMatchMc::ResoPartlyMatched);
          break;
        }
      }
    }
    // No physical channel expected in D*Lambda
    if (indexRec > -1) {
      auto particleReso = particlesMc.iteratorAt(indexRec);
      ptGen = particleReso.pt();
      invMassGen = computeInvMassGen(particlesMc, indexRec, pdg);
    }
    rowMcRecReduced(indexHfCandCharm, indexCandV0Bach,
                    flagReso, flagCharmBach,
                    flagCharmBachInterm, debugMcRec,
                    origin, ptGen, invMassGen,
                    nKinkedTracks);
  } else if constexpr (DType == DMesonType::Dplus) {
    vecDaughtersReso.push_back(tracks.rawIteratorAt(candCharmBach.prong0Id()));
    vecDaughtersReso.push_back(tracks.rawIteratorAt(candCharmBach.prong1Id()));
    vecDaughtersReso.push_back(tracks.rawIteratorAt(candCharmBach.prong2Id()));
    // Check if D+ is matched
    flagCharmBach = candCharmBach.flagMcMatchRec();
    flagCharmBachInterm = candCharmBach.flagMcDecayChanRec();
    if (flagCharmBach != 0) {
      SETBIT(debugMcRec, hf_decay::hf_cand_reso::PartialMatchMc::DplusMatched);
      origin = candCharmBach.originMcRec();
    }
    // Check if V0 is matched
    vecDaughtersReso.push_back(tracks.rawIteratorAt(bachelorV0.posTrackId()));
    vecDaughtersReso.push_back(tracks.rawIteratorAt(bachelorV0.negTrackId()));
    flagV0 = getMatchingFlagV0(particlesMc, std::array{vecDaughtersReso[3], vecDaughtersReso[4]});
    if (flagV0 != 0) {
      SETBIT(debugMcRec, std::abs(flagV0));
    }
    // If both D+ and K0s are matched, try to match resonance
    if (hf_decay::hf_cand_3prong::daughtersDplusMain.contains(static_cast<hf_decay::hf_cand_3prong::DecayChannelMain>(std::abs(flagCharmBach))) && flagV0 == hf_decay::hf_cand_reso::PartialMatchMc::K0Matched) {
      auto arrDaughtersReso = std::array{vecDaughtersReso[0], vecDaughtersReso[1], vecDaughtersReso[2], vecDaughtersReso[3], vecDaughtersReso[4]};
      auto pdgCodesDplusDaughters = hf_decay::hf_cand_3prong::daughtersDplusMain.at(static_cast<hf_decay::hf_cand_3prong::DecayChannelMain>(std::abs(flagCharmBach)));
      auto pdgCodesDaughters = std::array{pdgCodesDplusDaughters[0], pdgCodesDplusDaughters[1], pdgCodesDplusDaughters[2], +kPiPlus, -kPiPlus};
      for (const auto& [decayChannelFlag, pdgCodeReso] : hf_decay::hf_cand_reso::particlesToDplusK0s) {
        indexRec = RecoDecay::getMatchedMCRec<false, true, false, true, true>(particlesMc, arrDaughtersReso, pdgCodeReso, pdgCodesDaughters, true, &sign, 3, &nKinkedTracks);
        if (indexRec > -1) {
          flagReso = sign * decayChannelFlag;
          break;
        }
      }
      // Partial matching of Dsj -> D*K0s -> (D+ pi0) (K0s) with missing neutral
      if (indexRec < 0) {
        for (const auto& [decayChannelFlag, pdgCodeReso] : hf_decay::hf_cand_reso::particlesToDstarK0s) {
          indexRec = RecoDecay::getMatchedMCRec<false, true, true, true, true>(particlesMc, arrDaughtersReso, pdgCodeReso, pdgCodesDaughters, true, &sign, 3, &nKinkedTracks);
          if (indexRec > -1) {
            flagReso = sign * decayChannelFlag;
            SETBIT(debugMcRec, hf_decay::hf_cand_reso::PartialMatchMc::ResoPartlyMatched);
            break;
          }
        }
      }

    } else if (hf_decay::hf_cand_3prong::daughtersDplusMain.contains(static_cast<hf_decay::hf_cand_3prong::DecayChannelMain>(std::abs(flagCharmBach))) && std::abs(flagV0) == hf_decay::hf_cand_reso::PartialMatchMc::LambdaMatched) {
      // Peaking background of D+Lambda <- Ds* with spurious soft pion
      auto arrDaughtersReso = std::array{vecDaughtersReso[0], vecDaughtersReso[1], vecDaughtersReso[2], vecDaughtersReso[3], vecDaughtersReso[4]};
      auto pdgCodesDplusDaughters = hf_decay::hf_cand_3prong::daughtersDplusMain.at(static_cast<hf_decay::hf_cand_3prong::DecayChannelMain>(std::abs(flagCharmBach)));
      auto pdgCodesDaughters = std::array{pdgCodesDplusDaughters[0], pdgCodesDplusDaughters[1], pdgCodesDplusDaughters[2], +kProton, -kPiPlus};
      for (const auto& [decayChannelFlag, pdgCodeReso] : hf_decay::hf_cand_reso::particlesToDplusLambda) {
        indexRec = RecoDecay::getMatchedMCRec<false, true, false, true, true>(particlesMc, arrDaughtersReso, pdgCodeReso, pdgCodesDaughters, true, &sign, 3, &nKinkedTracks);
        if (indexRec > -1) {
          flagReso = sign * decayChannelFlag;
          break;
        }
      }
    }
    if (indexRec > -1) {
      auto particleReso = particlesMc.iteratorAt(indexRec);
      ptGen = particleReso.pt();
      invMassGen = computeInvMassGen(particlesMc, indexRec, pdg);
    }
    rowMcRecReduced(indexHfCandCharm, indexCandV0Bach,
                    flagReso, flagCharmBach,
                    flagCharmBachInterm, debugMcRec,
                    origin, ptGen, invMassGen,
                    nKinkedTracks);
  } else if constexpr (DType == DMesonType::D0) {
    vecDaughtersReso.push_back(tracks.rawIteratorAt(candCharmBach.prong0Id()));
    vecDaughtersReso.push_back(tracks.rawIteratorAt(candCharmBach.prong1Id()));
    // Check if D0 is matched
    flagCharmBach = candCharmBach.flagMcMatchRec();
    flagCharmBachInterm = candCharmBach.flagMcDecayChanRec();
    if (flagCharmBach != 0) {
      SETBIT(debugMcRec, hf_decay::hf_cand_reso::PartialMatchMc::D0Matched);
      origin = candCharmBach.originMcRec();
    }
    // Check if V0 is matched
    vecDaughtersReso.push_back(tracks.rawIteratorAt(bachelorV0.posTrackId()));
    vecDaughtersReso.push_back(tracks.rawIteratorAt(bachelorV0.negTrackId()));
    flagV0 = getMatchingFlagV0(particlesMc, std::array{vecDaughtersReso[2], vecDaughtersReso[3]});
    if (flagV0 != 0) {
      SETBIT(debugMcRec, std::abs(flagV0));
    }
    // No physical channel expected in D0 K0s
    // If both D0 and Lambda are matched, try to match resonance
    if (hf_decay::hf_cand_2prong::daughtersD0Main.contains(static_cast<hf_decay::hf_cand_2prong::DecayChannelMain>(std::abs(flagCharmBach))) && std::abs(flagV0) == hf_decay::hf_cand_reso::PartialMatchMc::LambdaMatched) {
      auto arrDaughtersReso = std::array{vecDaughtersReso[0], vecDaughtersReso[1], vecDaughtersReso[2], vecDaughtersReso[3]};
      auto pdgCodesDzeroDaughters = hf_decay::hf_cand_2prong::daughtersD0Main.at(static_cast<hf_decay::hf_cand_2prong::DecayChannelMain>(std::abs(flagCharmBach)));
      auto pdgCodesDaughters = std::array{pdgCodesDzeroDaughters[0], pdgCodesDzeroDaughters[1], +kProton, -kPiPlus};
      for (const auto& [decayChannelFlag, pdgCodeReso] : hf_decay::hf_cand_reso::particlesToD0Lambda) {
        indexRec = RecoDecay::getMatchedMCRec<false, true, false, true, true>(particlesMc, arrDaughtersReso, pdgCodeReso, pdgCodesDaughters, true, &sign, 3, &nKinkedTracks);
        if (indexRec > -1) {
          flagReso = sign * decayChannelFlag;
          break;
        }
      }
    }
    if (indexRec > -1) {
      auto particleReso = particlesMc.iteratorAt(indexRec);
      ptGen = particleReso.pt();
      invMassGen = computeInvMassGen(particlesMc, indexRec, pdg);
    }
    rowMcRecReduced(indexHfCandCharm, indexCandV0Bach,
                    flagReso, flagCharmBach,
                    flagCharmBachInterm, debugMcRec,
                    origin, ptGen, invMassGen,
                    nKinkedTracks);
  }
  registry.fill(HIST("hMCRecDebug"), debugMcRec);
  if (indexRec > -1) {
    registry.fill(HIST("hMCRecCounter"), flagReso);
    registry.fill(HIST("hMCRecOrigin"), origin);
    registry.fill(HIST("hMCRecMassGen"), invMassGen);
  }
  if (flagCharmBach != 0) {
    registry.fill(HIST("hMCRecCharmDau"), flagCharmBach);
  }
}

// Function for filling MC reco information of D Track candidates in the tables
/// \tparam dType is the D meson type (Dstar, Dplus or D0)
/// \param particlesMc is the table with MC particles
/// \param candCharmBach is the D meson candidate
/// \param bachelorTrack is the bachelor track
/// \param tracks is the table with tracks
/// \param indexHfCandCharm is the index of the charm-hadron bachelor in the reduced table
/// \param indexCandTrBach is the index of the v0 bachelor in the reduced table
/// \param pdg is the O2DatabasePDG service
/// \param registry is the histogram registry
/// \param rowMcRecReduced is the table to be filled
template <DMesonType DType, typename PParticles, typename CCand, typename BBachTr, typename Tr, typename Table>
void fillMcRecoInfoDTrack(PParticles const& particlesMc,
                          CCand const& candCharmBach,
                          BBachTr const& bachelorTrack,
                          Tr const& tracks,
                          const int64_t indexHfCandCharm,
                          const int64_t indexCandTrBach,
                          o2::framework::Service<o2::framework::O2DatabasePDG> const& pdg,
                          o2::framework::HistogramRegistry& registry,
                          Table& rowMcRecReduced)
{
  std::vector<typename Tr::iterator> vecDaughtersReso{};
  int8_t sign{0}, nKinkedTracks{0}, origin{0}, flagCharmBach{0}, flagCharmBachInterm{0}, flagTrack{0}, flagReso{0};
  int indexRec{-1};
  uint16_t debugMcRec{0};
  float ptGen{-1.f}, invMassGen{-1.f};
  if constexpr (DType == DMesonType::Dstar) {
    vecDaughtersReso.push_back(tracks.rawIteratorAt(candCharmBach.prong0Id()));
    vecDaughtersReso.push_back(tracks.rawIteratorAt(candCharmBach.prong1Id()));
    vecDaughtersReso.push_back(tracks.rawIteratorAt(candCharmBach.prongPiId()));
    // Check if D* is matched
    flagCharmBach = candCharmBach.flagMcMatchRec();
    if (flagCharmBach != 0) {
      SETBIT(debugMcRec, hf_decay::hf_cand_reso::PartialMatchMc::DstarMatched);
      origin = candCharmBach.originMcRec();
    }
    // Check if D0 is matched
    flagCharmBachInterm = candCharmBach.flagMcMatchRecD0();
    if (flagCharmBachInterm != 0) {
      SETBIT(debugMcRec, hf_decay::hf_cand_reso::PartialMatchMc::D0Matched);
    }
    // Check if Track is matched
    flagTrack = getMatchingFlagTrack(bachelorTrack);
    if (flagTrack != 0) {
      SETBIT(debugMcRec, flagTrack);
    }
    // If both D* and Track are matched, try to match resonance
    if (flagCharmBach != 0 && flagTrack == hf_decay::hf_cand_reso::PartialMatchMc::PionMatched) {
      auto arrDaughtersReso = std::array{vecDaughtersReso[0], vecDaughtersReso[1], vecDaughtersReso[2], bachelorTrack};
      auto pdgCodesDaughters = std::array{+kPiPlus, -kKPlus, +kPiPlus, -kPiPlus};
      for (const auto& [decayChannelFlag, pdgCodeReso] : hf_decay::hf_cand_reso::particlesToDstarPi) {
        indexRec = RecoDecay::getMatchedMCRec<false, true, false, true, true>(particlesMc, arrDaughtersReso, pdgCodeReso, pdgCodesDaughters, true, &sign, 3, &nKinkedTracks);
        if (indexRec > -1) {
          flagReso = sign * decayChannelFlag;
          break;
        }
      }
    }
    // No channels in D*K+ or D*Pr
    if (indexRec > -1) {
      auto particleReso = particlesMc.iteratorAt(indexRec);
      ptGen = particleReso.pt();
      invMassGen = computeInvMassGen(particlesMc, indexRec, pdg);
    }
    rowMcRecReduced(indexHfCandCharm, indexCandTrBach,
                    flagReso, flagCharmBach,
                    flagCharmBachInterm, debugMcRec,
                    origin, ptGen, invMassGen,
                    nKinkedTracks);
  } else if constexpr (DType == DMesonType::Dplus) {
    vecDaughtersReso.push_back(tracks.rawIteratorAt(candCharmBach.prong0Id()));
    vecDaughtersReso.push_back(tracks.rawIteratorAt(candCharmBach.prong1Id()));
    vecDaughtersReso.push_back(tracks.rawIteratorAt(candCharmBach.prong2Id()));
    // Check if D+ is matched
    flagCharmBach = candCharmBach.flagMcMatchRec();
    flagCharmBachInterm = candCharmBach.flagMcDecayChanRec();
    if (flagCharmBach != 0) {
      SETBIT(debugMcRec, hf_decay::hf_cand_reso::PartialMatchMc::DplusMatched);
      origin = candCharmBach.originMcRec();
    }
    // Check if Track is matched
    flagTrack = getMatchingFlagTrack(bachelorTrack);
    if (flagTrack != 0) {
      SETBIT(debugMcRec, flagTrack);
    }
    // If both D+ and Track are matched, try to match resonance
    if (hf_decay::hf_cand_3prong::daughtersDplusMain.contains(static_cast<hf_decay::hf_cand_3prong::DecayChannelMain>(std::abs(flagCharmBach))) && flagTrack == hf_decay::hf_cand_reso::PartialMatchMc::PionMatched) {
      auto arrDaughtersReso = std::array{vecDaughtersReso[0], vecDaughtersReso[1], vecDaughtersReso[2], bachelorTrack};
      auto pdgCodesDplusDaughters = hf_decay::hf_cand_3prong::daughtersDplusMain.at(static_cast<hf_decay::hf_cand_3prong::DecayChannelMain>(std::abs(flagCharmBach)));
      auto pdgCodesDaughters = std::array{pdgCodesDplusDaughters[0], pdgCodesDplusDaughters[1], pdgCodesDplusDaughters[2], -kPiPlus};
      for (const auto& [decayChannelFlag, pdgCodeReso] : hf_decay::hf_cand_reso::particlesToDplusPi) {
        indexRec = RecoDecay::getMatchedMCRec<false, true, false, true, true>(particlesMc, arrDaughtersReso, pdgCodeReso, pdgCodesDaughters, true, &sign, 3, &nKinkedTracks);
        if (indexRec > -1) {
          flagReso = sign * decayChannelFlag;
          break;
        }
      }
      // Partial matching of Dj -> D*Pi -> (D+ pi0) (pi) with missing neutral
      if (indexRec < 0) {
        for (const auto& [decayChannelFlag, pdgCodeReso] : hf_decay::hf_cand_reso::particlesToDstarPi) {
          indexRec = RecoDecay::getMatchedMCRec<false, true, true, true, true>(particlesMc, arrDaughtersReso, pdgCodeReso, pdgCodesDaughters, true, &sign, 3, &nKinkedTracks);
          if (indexRec > -1) {
            flagReso = sign * decayChannelFlag;
            SETBIT(debugMcRec, hf_decay::hf_cand_reso::PartialMatchMc::ResoPartlyMatched);
            break;
          }
        }
      }
    }
    // No channels in D+K+ or D+Pr
    if (indexRec > -1) {
      auto particleReso = particlesMc.iteratorAt(indexRec);
      ptGen = particleReso.pt();
      invMassGen = computeInvMassGen(particlesMc, indexRec, pdg);
    }
    rowMcRecReduced(indexHfCandCharm, indexCandTrBach,
                    flagReso, flagCharmBach,
                    flagCharmBachInterm, debugMcRec,
                    origin, ptGen, invMassGen,
                    nKinkedTracks);
  } else if constexpr (DType == DMesonType::D0) {
    vecDaughtersReso.push_back(tracks.rawIteratorAt(candCharmBach.prong0Id()));
    vecDaughtersReso.push_back(tracks.rawIteratorAt(candCharmBach.prong1Id()));
    // Check if D0 is matched
    flagCharmBach = candCharmBach.flagMcMatchRec();
    flagCharmBachInterm = candCharmBach.flagMcDecayChanRec();
    if (flagCharmBach != 0) {
      SETBIT(debugMcRec, hf_decay::hf_cand_reso::PartialMatchMc::D0Matched);
      origin = candCharmBach.originMcRec();
    }
    flagTrack = getMatchingFlagTrack(bachelorTrack);
    if (flagTrack != 0) {
      SETBIT(debugMcRec, flagTrack);
    }
    if (hf_decay::hf_cand_2prong::daughtersD0Main.contains(static_cast<hf_decay::hf_cand_2prong::DecayChannelMain>(std::abs(flagCharmBach))) && flagTrack == hf_decay::hf_cand_reso::PartialMatchMc::PionMatched) {
      auto arrDaughtersReso = std::array{vecDaughtersReso[0], vecDaughtersReso[1], bachelorTrack};
      auto pdgCodesDzeroDaughters = hf_decay::hf_cand_2prong::daughtersD0Main.at(static_cast<hf_decay::hf_cand_2prong::DecayChannelMain>(std::abs(flagCharmBach)));
      auto pdgCodesDaughters = std::array{pdgCodesDzeroDaughters[0], pdgCodesDzeroDaughters[1], +kPiPlus};
      for (const auto& [decayChannelFlag, pdgCodeReso] : hf_decay::hf_cand_reso::particlesToD0Pi) {
        indexRec = RecoDecay::getMatchedMCRec<false, true, false, true, true>(particlesMc, arrDaughtersReso, pdgCodeReso, pdgCodesDaughters, true, &sign, 3, &nKinkedTracks);
        if (indexRec > -1) {
          flagReso = sign * decayChannelFlag;
          break;
        }
      }
      // Partial matching of Dj -> D*Pi -> (D0 pi) (pi) with missing pion
      if (indexRec < 0) {
        for (const auto& [decayChannelFlag, pdgCodeReso] : hf_decay::hf_cand_reso::particlesToDstarPi) {
          indexRec = RecoDecay::getMatchedMCRec<false, true, true, true, true>(particlesMc, arrDaughtersReso, pdgCodeReso, pdgCodesDaughters, true, &sign, 3, &nKinkedTracks);
          if (indexRec > -1) {
            flagReso = sign * decayChannelFlag;
            SETBIT(debugMcRec, hf_decay::hf_cand_reso::PartialMatchMc::ResoPartlyMatched);
            break;
          }
        }
      }
    } else if (hf_decay::hf_cand_2prong::daughtersD0Main.contains(static_cast<hf_decay::hf_cand_2prong::DecayChannelMain>(std::abs(flagCharmBach))) && flagTrack == hf_decay::hf_cand_reso::PartialMatchMc::KaonMatched) {
      auto arrDaughtersReso = std::array{vecDaughtersReso[0], vecDaughtersReso[1], bachelorTrack};
      auto pdgCodesDzeroDaughters = hf_decay::hf_cand_2prong::daughtersD0Main.at(static_cast<hf_decay::hf_cand_2prong::DecayChannelMain>(std::abs(flagCharmBach)));
      auto pdgCodesDaughters = std::array{pdgCodesDzeroDaughters[0], pdgCodesDzeroDaughters[1], +kKPlus};
      for (const auto& [decayChannelFlag, pdgCodeReso] : hf_decay::hf_cand_reso::particlesToD0Kplus) {
        indexRec = RecoDecay::getMatchedMCRec<false, true, false, true, true>(particlesMc, arrDaughtersReso, pdgCodeReso, pdgCodesDaughters, true, &sign, 3, &nKinkedTracks);
        if (indexRec > -1) {
          flagReso = sign * decayChannelFlag;
          break;
        }
      }
    }
    if (indexRec > -1) {
      auto particleReso = particlesMc.iteratorAt(indexRec);
      ptGen = particleReso.pt();
      invMassGen = computeInvMassGen(particlesMc, indexRec, pdg);
    }
    rowMcRecReduced(indexHfCandCharm, indexCandTrBach,
                    flagReso, flagCharmBach,
                    flagCharmBachInterm, debugMcRec,
                    origin, ptGen, invMassGen,
                    nKinkedTracks);
  }
  registry.fill(HIST("hMCRecDebug"), debugMcRec);
  if (indexRec > -1) {
    registry.fill(HIST("hMCRecCounter"), flagReso);
    registry.fill(HIST("hMCRecOrigin"), origin);
    registry.fill(HIST("hMCRecMassGen"), invMassGen);
  }
  if (flagCharmBach != 0) {
    registry.fill(HIST("hMCRecCharmDau"), flagCharmBach);
  }
} // fillMcRecoInfoDTrack

// Function for derived data creation
/// \tparam dType is the D meson type (Dstar, Dplus or D0)
/// \param collision is the current collision
/// \param candsD are the D meson candidates in the current collision
/// \param bachelorV0s are the V0 candidates in the current collision
/// \param bachelorTrks are the track ids in the current collision
/// \param tracks is the track table
/// \param tracksIU is the trackIU table
/// \param particlesMc is the MC particle table
/// \param hfRejMap is the event rejection map from the HF event selection util
/// \param bz is the magnetic field
/// \param pdg is the O2DatabasePDG service
/// \param registry is the histogram registry
/// \param matCorr is the material correction type to be used in the track propagation
/// \param fitter is the DCAFitter object
/// \param rejectPairsWithCommonDaughter is a flag to activate rejection of pairs sharing a daughter track
/// \param rowCollisionReduced is the collision reduced table to be filled
/// \param rowCandDmesReduced is the D-meson reduced table to be filled
/// \param rowCandV0Reduced is the V0 reduced table to be filled
/// \param rowTrkReduced is the track reduced table to be filled
/// \param rowMcRecV0Reduced is the MC reco D-V0 reduced table to be filled
/// \param rowMcRecTrkReduced is the MC reco D-track reduced table to be filled
/// \param rowCandDmesMlReduced is the ML reduced table to be filled
/// \param vDriftMgr is the TPC velocity drift manager object
template <bool WithMl, bool DoMc, DMesonType DType, PairingType PairType, class BCs, class Colls, typename Coll, typename CCands, typename Tr, typename TrIU, typename PParticles, typename BBachV0s, typename BBachTracks, typename DmesCuts, typename TrkCuts, typename V0Cuts, typename GammaCuts, typename QaConfig, typename TableCollRed, typename TableCandDRed, typename TableCandV0Red, typename TableTrkRed, typename TableMcRecV0Red, typename TableMcRecTrkRed, typename TableCandDMlRed>
void runDataCreation(Coll const& collision,
                     CCands const& candsD,
                     BBachV0s const& bachelorV0s,
                     BBachTracks const& bachelorTrks,
                     Tr const& tracks,
                     TrIU const& tracksIU,
                     PParticles const& particlesMc,
                     o2::hf_evsel::HfCollisionRejectionMask hfRejMap,
                     const float bz,
                     o2::framework::Service<o2::framework::O2DatabasePDG> const& pdg,
                     o2::framework::HistogramRegistry& registry,
                     o2::base::Propagator::MatCorrType const& matCorr,
                     o2::vertexing::DCAFitterN<2>& fitter,
                     DmesCuts const& cfgDmesCuts,
                     TrkCuts const& cfgSingleTrackCuts,
                     V0Cuts const& cfgV0Cuts,
                     GammaCuts const& cfgGammaCuts,
                     QaConfig const& cfgQaPlots,
                     bool rejectPairsWithCommonDaughter,
                     TableCollRed& rowCollisionReduced,
                     TableCandDRed& rowCandDmesReduced,
                     TableCandV0Red& rowCandV0Reduced,
                     TableTrkRed& rowTrkReduced,
                     TableMcRecV0Red& rowMcRecV0Reduced,
                     TableMcRecTrkRed& rowMcRecTrkReduced,
                     TableCandDMlRed& rowCandDmesMlReduced,
                     o2::aod::common::TPCVDriftManager* vDriftMgr = nullptr)
{
  int const indexHfReducedCollision = rowCollisionReduced.lastIndex() + 1;
  // std::map where the key is the V0.globalIndex() and
  // the value is the V0 index in the table of the selected v0s
  std::map<int64_t, int64_t> selectedV0s;
  std::map<int64_t, int64_t> selectedTracks;
  std::map<int64_t, int64_t> selectedGammas;
  bool fillHfReducedCollision = false;
  constexpr bool DoTracks = PairType == PairingType::TrackOnly || PairType == PairingType::V0AndTrack;
  constexpr bool DoV0s = PairType == PairingType::V0Only || PairType == PairingType::V0AndTrack;
  constexpr bool DoGammas = PairType == PairingType::GammaOnly;
  // loop on D candidates
  for (const auto& candD : candsD) {
    // initialize variables depending on D meson type
    bool fillHfCandD = false;
    std::array<float, 3> secondaryVertexD{};
    std::array<int, 3> prongIdsD{};
    std::array<float, 6> bdtScores = {-1.f, -1.f, -1.f, -1.f, -1.f, -1.f};
    std::vector<std::decay_t<typename TrIU::iterator>> charmHadDauTracks{};
    HfResoVarContainer varUtils;
    varUtils.ptD = candD.pt();
    if constexpr (DType == DMesonType::Dstar) {
      varUtils.signD = candD.signSoftPi();
      if (varUtils.signD > 0) {
        varUtils.invMassD = candD.invMassDstar();
        varUtils.invMassD0 = candD.invMassD0();
      } else {
        varUtils.invMassD = candD.invMassAntiDstar();
        varUtils.invMassD0 = candD.invMassD0Bar();
      }
      secondaryVertexD[0] = candD.xSecondaryVertexD0();
      secondaryVertexD[1] = candD.ySecondaryVertexD0();
      secondaryVertexD[2] = candD.zSecondaryVertexD0();
      prongIdsD[0] = candD.prong0Id();
      prongIdsD[1] = candD.prong1Id();
      prongIdsD[2] = candD.prongPiId();
      varUtils.pVectorProng0 = candD.pVectorProng0();
      varUtils.pVectorProng1 = candD.pVectorProng1();
      varUtils.pVectorProng2 = candD.pVecSoftPi();
      charmHadDauTracks.push_back(tracksIU.rawIteratorAt(candD.prong0Id()));
      charmHadDauTracks.push_back(tracksIU.rawIteratorAt(candD.prong1Id()));
      if constexpr (WithMl) {
        std::copy(candD.mlProbDstarToD0Pi().begin(), candD.mlProbDstarToD0Pi().end(), bdtScores.begin());
      }
      registry.fill(HIST("hMassVsPtDstarAll"), varUtils.ptD, varUtils.invMassD - varUtils.invMassD0);
    } else if constexpr (DType == DMesonType::Dplus) {
      auto prong0 = tracksIU.rawIteratorAt(candD.prong0Id());
      varUtils.invMassD = HfHelper::invMassDplusToPiKPi(candD);
      secondaryVertexD[0] = candD.xSecondaryVertex();
      secondaryVertexD[1] = candD.ySecondaryVertex();
      secondaryVertexD[2] = candD.zSecondaryVertex();
      prongIdsD[0] = candD.prong0Id();
      prongIdsD[1] = candD.prong1Id();
      prongIdsD[2] = candD.prong2Id();
      varUtils.signD = prong0.sign();
      varUtils.pVectorProng0 = candD.pVectorProng0();
      varUtils.pVectorProng1 = candD.pVectorProng1();
      varUtils.pVectorProng2 = candD.pVectorProng2();
      charmHadDauTracks.push_back(tracksIU.rawIteratorAt(candD.prong0Id()));
      charmHadDauTracks.push_back(tracksIU.rawIteratorAt(candD.prong1Id()));
      charmHadDauTracks.push_back(tracksIU.rawIteratorAt(candD.prong2Id()));
      if constexpr (WithMl) {
        std::copy(candD.mlProbDplusToPiKPi().begin(), candD.mlProbDplusToPiKPi().end(), bdtScores.begin());
      }
      registry.fill(HIST("hMassVsPtDplusAll"), varUtils.ptD, varUtils.invMassD);
    } else if constexpr (DType == DMesonType::D0) {
      varUtils.invMassD0 = HfHelper::invMassD0ToPiK(candD);
      varUtils.invMassD0Bar = HfHelper::invMassD0barToKPi(candD);
      secondaryVertexD[0] = candD.xSecondaryVertex();
      secondaryVertexD[1] = candD.ySecondaryVertex();
      secondaryVertexD[2] = candD.zSecondaryVertex();
      prongIdsD[0] = candD.prong0Id();
      prongIdsD[1] = candD.prong1Id();
      prongIdsD[2] = -1; // D0 does not have a third prong
      charmHadDauTracks.push_back(tracksIU.rawIteratorAt(candD.prong0Id()));
      charmHadDauTracks.push_back(tracksIU.rawIteratorAt(candD.prong1Id()));
      varUtils.pVectorProng0 = candD.pVectorProng0();
      varUtils.pVectorProng1 = candD.pVectorProng1();
      varUtils.pVectorProng2 = {0.f, 0.f, 0.f}; // D0 does not have a third prong
      if constexpr (WithMl) {
        std::copy(candD.mlProbD0().begin(), candD.mlProbD0().end(), bdtScores.begin());
        std::copy(candD.mlProbD0bar().begin(), candD.mlProbD0bar().end(), bdtScores.begin() + 3);
      }
      if (candD.isSelD0() >= cfgDmesCuts.selectionFlagD0.value) {
        registry.fill(HIST("hMassVsPtD0All"), varUtils.ptD, varUtils.invMassD0);
      }
      if (candD.isSelD0bar() >= cfgDmesCuts.selectionFlagD0Bar.value) {
        registry.fill(HIST("hMassVsPtD0BarAll"), varUtils.ptD, varUtils.invMassD0Bar);
      }
    } // end of dType switch

    // Get single track variables
    float chi2TpcDauMax = -1.f;
    int nItsClsDauMin = 8, nTpcCrossRowsDauMin = 200;
    float chi2TpcSoftPi = -1.f;
    int nItsClsSoftPi = 8, nTpcCrossRowsSoftPi = 200;
    for (const auto& charmHadTrack : charmHadDauTracks) {
      if (charmHadTrack.itsNCls() < nItsClsDauMin) {
        nItsClsDauMin = charmHadTrack.itsNCls();
      }
      if (charmHadTrack.tpcNClsCrossedRows() < nTpcCrossRowsDauMin) {
        nTpcCrossRowsDauMin = charmHadTrack.tpcNClsCrossedRows();
      }
      if (charmHadTrack.tpcChi2NCl() > chi2TpcDauMax) {
        chi2TpcDauMax = charmHadTrack.tpcChi2NCl();
      }
    }
    if constexpr (DType == DMesonType::Dstar) {
      auto softPi = tracksIU.rawIteratorAt(candD.prongPiId());
      nItsClsSoftPi = softPi.itsNCls();
      nTpcCrossRowsSoftPi = softPi.tpcNClsCrossedRows();
      chi2TpcSoftPi = softPi.tpcChi2NCl();
      charmHadDauTracks.push_back(softPi);
    }
    // Loop on the bachelor V0s
    if constexpr (DoV0s) {
      for (const auto& v0 : bachelorV0s) {
        auto trackPos = tracksIU.rawIteratorAt(v0.posTrackId());
        auto trackNeg = tracksIU.rawIteratorAt(v0.negTrackId());
        // Apply selsection
        auto v0DauTracks = std::array{trackPos, trackNeg};
        HfResoCandidateV0 candV0;
        if (!buildAndSelectV0(collision, prongIdsD, v0DauTracks, cfgV0Cuts, fitter, candV0, rejectPairsWithCommonDaughter)) {
          continue;
        }
        // Get single track variables
        float chi2TpcDauV0Max = -1.f;
        int nItsClsDauV0Min = 8, nTpcCrossRowsDauV0Min = 200;
        for (const auto& v0Track : v0DauTracks) {
          if (v0Track.itsNCls() < nItsClsDauV0Min) {
            nItsClsDauV0Min = v0Track.itsNCls();
          }
          if (v0Track.tpcNClsCrossedRows() < nTpcCrossRowsDauV0Min) {
            nTpcCrossRowsDauV0Min = v0Track.tpcNClsCrossedRows();
          }
          if (v0Track.tpcChi2NCl() > chi2TpcDauV0Max) {
            chi2TpcDauV0Max = v0Track.tpcChi2NCl();
          }
        }
        // propagate V0 to primary vertex (if enabled)
        if (cfgV0Cuts.propagateV0toPV.value) {
          std::array<float, 3> const pVecV0Orig = {candV0.mom[0], candV0.mom[1], candV0.mom[2]};
          std::array<float, 2> dcaInfo{};
          auto trackParK0 = o2::track::TrackPar(candV0.pos, pVecV0Orig, 0, true);
          trackParK0.setPID(o2::track::PID::K0);
          trackParK0.setAbsCharge(0);
          o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParK0, 2.f, matCorr, &dcaInfo);
          getPxPyPz(trackParK0, candV0.mom);
        }
        // compute resonance invariant mass and filling of QA histograms
        registry.fill(HIST("hAP"), candV0.alpha, candV0.qt);
        registry.fill(HIST("hV0Radius"), candV0.radius);
        if (TESTBIT(candV0.v0Type, BachelorType::K0s)) {
          registry.fill(HIST("hMassVsPtK0s"), candV0.pT, candV0.mK0Short);
          if constexpr (DType == DMesonType::Dstar) {
            varUtils.ptReso = RecoDecay::pt(RecoDecay::sumOfVec(varUtils.pVectorProng0, varUtils.pVectorProng1, varUtils.pVectorProng2, candV0.mom));
            if (varUtils.signD > 0) {
              varUtils.invMassReso = RecoDecay::m(std::array{varUtils.pVectorProng0, varUtils.pVectorProng1, varUtils.pVectorProng2, candV0.mom}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassKPlus, o2::constants::physics::MassPionCharged, o2::constants::physics::MassK0});
            } else {
              varUtils.invMassReso = RecoDecay::m(std::array{varUtils.pVectorProng1, varUtils.pVectorProng0, varUtils.pVectorProng2, candV0.mom}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassKPlus, o2::constants::physics::MassPionCharged, o2::constants::physics::MassK0});
            }
            if (!cfgQaPlots.applyCutsForQaHistograms.value ||
                (varUtils.invMassD - varUtils.invMassD0 > cfgQaPlots.cutMassDstarMin.value &&
                 varUtils.invMassD - varUtils.invMassD0 < cfgQaPlots.cutMassDstarMax.value &&
                 candV0.mK0Short > cfgQaPlots.cutMassK0sMin.value &&
                 candV0.mK0Short < cfgQaPlots.cutMassK0sMax.value)) {
              registry.fill(HIST("hMassDstarK0s"), varUtils.ptReso, varUtils.invMassReso - varUtils.invMassD);
            }
          } else if constexpr (DType == DMesonType::Dplus) {
            varUtils.invMassReso = RecoDecay::m(std::array{varUtils.pVectorProng0, varUtils.pVectorProng1, varUtils.pVectorProng2, candV0.mom}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassKPlus, o2::constants::physics::MassPionCharged, o2::constants::physics::MassK0});
            varUtils.ptReso = RecoDecay::pt(RecoDecay::sumOfVec(varUtils.pVectorProng0, varUtils.pVectorProng1, varUtils.pVectorProng2, candV0.mom));
            if (!cfgQaPlots.applyCutsForQaHistograms.value ||
                (varUtils.invMassD > cfgQaPlots.cutMassDMin.value &&
                 varUtils.invMassD < cfgQaPlots.cutMassDMax.value &&
                 candV0.mK0Short > cfgQaPlots.cutMassK0sMin.value &&
                 candV0.mK0Short < cfgQaPlots.cutMassK0sMax.value)) {
              registry.fill(HIST("hMassDplusK0s"), varUtils.ptReso, varUtils.invMassReso - varUtils.invMassD);
            }
          } // end of dType switch
        } // matched with K0s
        bool const isLambda = TESTBIT(candV0.v0Type, BachelorType::Lambda);
        bool const isAntiLambda = TESTBIT(candV0.v0Type, BachelorType::AntiLambda);
        if (isLambda || isAntiLambda) {
          registry.fill(HIST("hMassVsPtLambda"), candV0.pT, candV0.mLambda);
          if constexpr (DType == DMesonType::Dstar) {
            varUtils.ptReso = RecoDecay::pt(RecoDecay::sumOfVec(varUtils.pVectorProng0, varUtils.pVectorProng1, varUtils.pVectorProng2, candV0.mom));
            if (varUtils.signD > 0) {
              varUtils.invMassReso = RecoDecay::m(std::array{varUtils.pVectorProng0, varUtils.pVectorProng1, varUtils.pVectorProng2, candV0.mom}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassKPlus, o2::constants::physics::MassPionCharged, o2::constants::physics::MassLambda});
            } else {
              varUtils.invMassReso = RecoDecay::m(std::array{varUtils.pVectorProng1, varUtils.pVectorProng0, varUtils.pVectorProng2, candV0.mom}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassKPlus, o2::constants::physics::MassPionCharged, o2::constants::physics::MassLambda});
            }
            if (!cfgQaPlots.applyCutsForQaHistograms.value ||
                (varUtils.invMassD - varUtils.invMassD0 > cfgQaPlots.cutMassDstarMin.value &&
                 varUtils.invMassD - varUtils.invMassD0 < cfgQaPlots.cutMassDstarMax.value &&
                 candV0.mLambda > cfgQaPlots.cutMassLambdaMin.value &&
                 candV0.mLambda < cfgQaPlots.cutMassLambdaMax.value)) {
              registry.fill(HIST("hMassDstarLambda"), varUtils.ptReso, varUtils.invMassReso - varUtils.invMassD);
            }
          } else if constexpr (DType == DMesonType::Dplus) {
            varUtils.invMassReso = RecoDecay::m(std::array{varUtils.pVectorProng0, varUtils.pVectorProng1, varUtils.pVectorProng2, candV0.mom}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassKPlus, o2::constants::physics::MassPionCharged, o2::constants::physics::MassLambda});
            varUtils.ptReso = RecoDecay::pt(RecoDecay::sumOfVec(varUtils.pVectorProng0, varUtils.pVectorProng1, varUtils.pVectorProng2, candV0.mom));
            if (!cfgQaPlots.applyCutsForQaHistograms.value ||
                (varUtils.invMassD > cfgQaPlots.cutMassDMin.value &&
                 varUtils.invMassD < cfgQaPlots.cutMassDMax.value &&
                 candV0.mLambda > cfgQaPlots.cutMassLambdaMin.value &&
                 candV0.mLambda < cfgQaPlots.cutMassLambdaMax.value)) {
              registry.fill(HIST("hMassDplusLambda"), varUtils.ptReso, varUtils.invMassReso - varUtils.invMassD);
            }
          } else if constexpr (DType == DMesonType::D0) {
            if (isLambda) {
              varUtils.invMassReso = RecoDecay::m(std::array{varUtils.pVectorProng0, varUtils.pVectorProng1, candV0.mom}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassKPlus, o2::constants::physics::MassLambda});
            } else {
              varUtils.invMassReso = RecoDecay::m(std::array{varUtils.pVectorProng1, varUtils.pVectorProng0, candV0.mom}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassKPlus, o2::constants::physics::MassLambda});
            }
            varUtils.ptReso = RecoDecay::pt(RecoDecay::sumOfVec(varUtils.pVectorProng0, varUtils.pVectorProng1, candV0.mom));
            if (!cfgQaPlots.applyCutsForQaHistograms.value ||
                (((varUtils.invMassD0 > cfgQaPlots.cutMassDMin.value && varUtils.invMassD0 < cfgQaPlots.cutMassDMax.value) ||
                  (varUtils.invMassD0Bar > cfgQaPlots.cutMassDMin.value && varUtils.invMassD0Bar < cfgQaPlots.cutMassDMax.value)) &&
                 candV0.mLambda > cfgQaPlots.cutMassLambdaMin.value &&
                 candV0.mLambda < cfgQaPlots.cutMassLambdaMax.value)) {
              if (isLambda) {
                registry.fill(HIST("hMassD0Lambda"), varUtils.ptReso, varUtils.invMassReso - varUtils.invMassD0);
              } else {
                registry.fill(HIST("hMassD0Lambda"), varUtils.ptReso, varUtils.invMassReso - varUtils.invMassD0Bar);
              }
            }
          } // end of dType switch
        } // matched with Lambda or AntiLambda
        // fill V0 table
        // if information on V0 already stored, go to next V0
        if (!selectedV0s.count(v0.globalIndex())) {
          rowCandV0Reduced(trackPos.globalIndex(), trackNeg.globalIndex(),
                           indexHfReducedCollision,
                           candV0.pos[0], candV0.pos[1], candV0.pos[2],
                           candV0.momPos[0], candV0.momPos[1], candV0.momPos[2],
                           candV0.momNeg[0], candV0.momNeg[1], candV0.momNeg[2],
                           candV0.cosPA,
                           candV0.dcaV0ToPv,
                           nItsClsDauV0Min, nTpcCrossRowsDauV0Min, chi2TpcDauV0Max,
                           candV0.v0Type);
          selectedV0s[v0.globalIndex()] = rowCandV0Reduced.lastIndex();
        }
        fillHfCandD = true;
        // Optional filling of MC Rec table, for now only implemented for Ds1->D*K0s and Ds2*->D+K0s
        if constexpr (DoMc) {
          auto indexHfCandCharm = rowCandDmesReduced.lastIndex() + 1;
          fillMcRecoInfoDV0<DType>(particlesMc, candD, v0, tracksIU, indexHfCandCharm, selectedV0s[v0.globalIndex()], pdg, registry, rowMcRecV0Reduced);
        }
      } // end of loop on V0 candidates
    } // end of do V0s
    // Loop on the bachelor tracks
    if constexpr (DoTracks) {
      for (const auto& trackIndex : bachelorTrks) {
        auto track = tracks.rawIteratorAt(trackIndex.trackId());
        if (!isTrackSelected(track, prongIdsD, cfgSingleTrackCuts, rejectPairsWithCommonDaughter)) {
          continue;
        }
        // if the track has been reassociated, re-propagate it to PV (minor difference)
        auto trackParCovTrack = getTrackParCov(track);
        std::array<float, 2> dcaTrack{track.dcaXY(), track.dcaZ()};
        std::array<float, 3> pVecTrack = track.pVector();
        if (track.collisionId() != collision.globalIndex()) {
          o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParCovTrack, 2.f, matCorr, &dcaTrack);
          getPxPyPz(trackParCovTrack, pVecTrack);
        }
        registry.fill(HIST("hdEdxVsP"), track.p(), track.tpcSignal());
        // compute invariant mass and filling of QA histograms
        if constexpr (DType == DMesonType::Dstar) {
          // D* pi
          if (std::abs(track.tpcNSigmaPi()) < cfgSingleTrackCuts.maxNsigmaTpcPi.value) {
            if (varUtils.signD > 0 && track.sign() < 0) {
              varUtils.invMassReso = RecoDecay::m(std::array{varUtils.pVectorProng0, varUtils.pVectorProng1, varUtils.pVectorProng2, pVecTrack}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassKPlus, o2::constants::physics::MassPionCharged, o2::constants::physics::MassPionCharged});
            } else if (varUtils.signD < 0 && track.sign() > 0) {
              varUtils.invMassReso = RecoDecay::m(std::array{varUtils.pVectorProng1, varUtils.pVectorProng0, varUtils.pVectorProng2, pVecTrack}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassKPlus, o2::constants::physics::MassPionCharged, o2::constants::physics::MassPionCharged});
            } else {
              varUtils.invMassReso = -1.f; // invalid case
            }
            varUtils.ptReso = RecoDecay::pt(RecoDecay::sumOfVec(varUtils.pVectorProng0, varUtils.pVectorProng1, varUtils.pVectorProng2, pVecTrack));
            if (!cfgQaPlots.applyCutsForQaHistograms.value ||
                (varUtils.invMassD - varUtils.invMassD0 > cfgQaPlots.cutMassDstarMin.value &&
                 varUtils.invMassD - varUtils.invMassD0 < cfgQaPlots.cutMassDstarMax.value)) {
              registry.fill(HIST("hMassDstarPi"), varUtils.ptReso, varUtils.invMassReso - varUtils.invMassD);
            }
          }
          if (std::abs(track.tpcNSigmaKa()) < cfgSingleTrackCuts.maxNsigmaTpcKa.value) {
            if (varUtils.signD > 0 && track.sign() < 0) {
              varUtils.invMassReso = RecoDecay::m(std::array{varUtils.pVectorProng0, varUtils.pVectorProng1, varUtils.pVectorProng2, pVecTrack}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassKPlus, o2::constants::physics::MassPionCharged, o2::constants::physics::MassKPlus});
            } else if (varUtils.signD < 0 && track.sign() > 0) {
              varUtils.invMassReso = RecoDecay::m(std::array{varUtils.pVectorProng1, varUtils.pVectorProng0, varUtils.pVectorProng2, pVecTrack}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassKPlus, o2::constants::physics::MassPionCharged, o2::constants::physics::MassKPlus});
            } else {
              varUtils.invMassReso = -1.f; // invalid case
            }
            varUtils.ptReso = RecoDecay::pt(RecoDecay::sumOfVec(varUtils.pVectorProng0, varUtils.pVectorProng1, varUtils.pVectorProng2, pVecTrack));
            if (!cfgQaPlots.applyCutsForQaHistograms.value ||
                (varUtils.invMassD - varUtils.invMassD0 > cfgQaPlots.cutMassDstarMin.value &&
                 varUtils.invMassD - varUtils.invMassD0 < cfgQaPlots.cutMassDstarMax.value)) {
              registry.fill(HIST("hMassDstarK"), varUtils.ptReso, varUtils.invMassReso - varUtils.invMassD);
            }
          }
          // D* p
          if (std::abs(track.tpcNSigmaPr()) < cfgSingleTrackCuts.maxNsigmaTpcPr.value) {
            if (varUtils.signD > 0 && track.sign() > 0) {
              varUtils.invMassReso = RecoDecay::m(std::array{varUtils.pVectorProng0, varUtils.pVectorProng1, varUtils.pVectorProng2, pVecTrack}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassKPlus, o2::constants::physics::MassPionCharged, o2::constants::physics::MassProton});
            } else if (varUtils.signD < 0 && track.sign() < 0) {
              varUtils.invMassReso = RecoDecay::m(std::array{varUtils.pVectorProng1, varUtils.pVectorProng0, varUtils.pVectorProng2, pVecTrack}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassKPlus, o2::constants::physics::MassPionCharged, o2::constants::physics::MassProton});
            } else {
              varUtils.invMassReso = -1.f; // invalid case
            }
            varUtils.ptReso = RecoDecay::pt(RecoDecay::sumOfVec(varUtils.pVectorProng0, varUtils.pVectorProng1, varUtils.pVectorProng2, pVecTrack));
            if (!cfgQaPlots.applyCutsForQaHistograms.value ||
                (varUtils.invMassD - varUtils.invMassD0 > cfgQaPlots.cutMassDstarMin.value &&
                 varUtils.invMassD - varUtils.invMassD0 < cfgQaPlots.cutMassDstarMax.value)) {
              registry.fill(HIST("hMassDstarProton"), varUtils.ptReso, varUtils.invMassReso - varUtils.invMassD);
            }
          }
        } else if constexpr (DType == DMesonType::Dplus) {
          // D+ pi
          if (std::abs(track.tpcNSigmaPi()) < cfgSingleTrackCuts.maxNsigmaTpcPi.value) {
            if (varUtils.signD * track.sign() < 0) {
              varUtils.invMassReso = RecoDecay::m(std::array{varUtils.pVectorProng0, varUtils.pVectorProng1, varUtils.pVectorProng2, pVecTrack}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassKPlus, o2::constants::physics::MassPionCharged, o2::constants::physics::MassPionCharged});
            } else {
              varUtils.invMassReso = -1.f; // invalid case
            }
            varUtils.ptReso = RecoDecay::pt(RecoDecay::sumOfVec(varUtils.pVectorProng0, varUtils.pVectorProng1, varUtils.pVectorProng2, pVecTrack));
            if (!cfgQaPlots.applyCutsForQaHistograms.value ||
                (varUtils.invMassD > cfgQaPlots.cutMassDMin.value &&
                 varUtils.invMassD < cfgQaPlots.cutMassDMax.value)) {
              registry.fill(HIST("hMassDplusPi"), varUtils.ptReso, varUtils.invMassReso - varUtils.invMassD);
            }
          }
          // D+ K
          if (std::abs(track.tpcNSigmaKa()) < cfgSingleTrackCuts.maxNsigmaTpcKa.value) {
            if (varUtils.signD * track.sign() < 0) {
              varUtils.invMassReso = RecoDecay::m(std::array{varUtils.pVectorProng0, varUtils.pVectorProng1, varUtils.pVectorProng2, pVecTrack}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassKPlus, o2::constants::physics::MassPionCharged, o2::constants::physics::MassKPlus});
            } else {
              varUtils.invMassReso = -1.f; // invalid case
            }
            varUtils.ptReso = RecoDecay::pt(RecoDecay::sumOfVec(varUtils.pVectorProng0, varUtils.pVectorProng1, varUtils.pVectorProng2, pVecTrack));
            if (!cfgQaPlots.applyCutsForQaHistograms.value ||
                (varUtils.invMassD > cfgQaPlots.cutMassDMin.value &&
                 varUtils.invMassD < cfgQaPlots.cutMassDMax.value)) {
              registry.fill(HIST("hMassDplusK"), varUtils.ptReso, varUtils.invMassReso - varUtils.invMassD);
            }
          }
          // D+ pr
          if (std::abs(track.tpcNSigmaPr()) < cfgSingleTrackCuts.maxNsigmaTpcPr.value) {
            if (varUtils.signD * track.sign() < 0) {
              varUtils.invMassReso = RecoDecay::m(std::array{varUtils.pVectorProng0, varUtils.pVectorProng1, varUtils.pVectorProng2, pVecTrack}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassKPlus, o2::constants::physics::MassPionCharged, o2::constants::physics::MassProton});
            } else {
              varUtils.invMassReso = -1.f; // invalid case
            }
            varUtils.ptReso = RecoDecay::pt(RecoDecay::sumOfVec(varUtils.pVectorProng0, varUtils.pVectorProng1, varUtils.pVectorProng2, pVecTrack));
            if (!cfgQaPlots.applyCutsForQaHistograms.value ||
                (varUtils.invMassD > cfgQaPlots.cutMassDMin.value &&
                 varUtils.invMassD < cfgQaPlots.cutMassDMax.value)) {
              registry.fill(HIST("hMassDplusProton"), varUtils.ptReso, varUtils.invMassReso - varUtils.invMassD);
            }
          }
        } else if constexpr (DType == DMesonType::D0) {
          // D0 pi
          if (std::abs(track.tpcNSigmaPi()) < cfgSingleTrackCuts.maxNsigmaTpcPi.value) {
            if (track.sign() > 0) {
              varUtils.invMassReso = RecoDecay::m(std::array{varUtils.pVectorProng0, varUtils.pVectorProng1, pVecTrack}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassKPlus, o2::constants::physics::MassPionCharged});
            } else {
              varUtils.invMassReso = RecoDecay::m(std::array{varUtils.pVectorProng1, varUtils.pVectorProng0, pVecTrack}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassKPlus, o2::constants::physics::MassPionCharged});
            }
            varUtils.ptReso = RecoDecay::pt(RecoDecay::sumOfVec(varUtils.pVectorProng0, varUtils.pVectorProng1, pVecTrack));
            if (!cfgQaPlots.applyCutsForQaHistograms.value ||
                ((varUtils.invMassD0 > cfgQaPlots.cutMassDMin.value &&
                  varUtils.invMassD0 < cfgQaPlots.cutMassDMax.value) ||
                 (varUtils.invMassD0Bar > cfgQaPlots.cutMassDMin.value &&
                  varUtils.invMassD0Bar < cfgQaPlots.cutMassDMax.value))) {
              if (track.sign() > 0) {
                registry.fill(HIST("hMassD0Pi"), varUtils.ptReso, varUtils.invMassReso - varUtils.invMassD0);
              } else {
                registry.fill(HIST("hMassD0Pi"), varUtils.ptReso, varUtils.invMassReso - varUtils.invMassD0Bar);
              }
            }
          }
          // D0 K
          if (std::abs(track.tpcNSigmaKa()) < cfgSingleTrackCuts.maxNsigmaTpcKa.value) {
            if (track.sign() > 0) {
              varUtils.invMassReso = RecoDecay::m(std::array{varUtils.pVectorProng0, varUtils.pVectorProng1, pVecTrack}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassKPlus, o2::constants::physics::MassKPlus});
            } else {
              varUtils.invMassReso = RecoDecay::m(std::array{varUtils.pVectorProng1, varUtils.pVectorProng0, pVecTrack}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassKPlus, o2::constants::physics::MassKPlus});
            }
            varUtils.ptReso = RecoDecay::pt(RecoDecay::sumOfVec(varUtils.pVectorProng0, varUtils.pVectorProng1, pVecTrack));
            if (!cfgQaPlots.applyCutsForQaHistograms.value ||
                ((varUtils.invMassD0 > cfgQaPlots.cutMassDMin.value &&
                  varUtils.invMassD0 < cfgQaPlots.cutMassDMax.value) ||
                 (varUtils.invMassD0Bar > cfgQaPlots.cutMassDMin.value &&
                  varUtils.invMassD0Bar < cfgQaPlots.cutMassDMax.value))) {
              if (track.sign() > 0) {
                registry.fill(HIST("hMassD0K"), varUtils.ptReso, varUtils.invMassReso - varUtils.invMassD0);
              } else {
                registry.fill(HIST("hMassD0K"), varUtils.ptReso, varUtils.invMassReso - varUtils.invMassD0Bar);
              }
            }
          }
          // D0 p
          if (std::abs(track.tpcNSigmaPr()) < cfgSingleTrackCuts.maxNsigmaTpcPr.value) {
            if (track.sign() > 0) {
              varUtils.invMassReso = RecoDecay::m(std::array{varUtils.pVectorProng0, varUtils.pVectorProng1, pVecTrack}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassKPlus, o2::constants::physics::MassProton});
            } else {
              varUtils.invMassReso = RecoDecay::m(std::array{varUtils.pVectorProng1, varUtils.pVectorProng0, pVecTrack}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassKPlus, o2::constants::physics::MassProton});
            }
            varUtils.ptReso = RecoDecay::pt(RecoDecay::sumOfVec(varUtils.pVectorProng0, varUtils.pVectorProng1, pVecTrack));
            if (!cfgQaPlots.applyCutsForQaHistograms.value ||
                ((varUtils.invMassD0 > cfgQaPlots.cutMassDMin.value &&
                  varUtils.invMassD0 < cfgQaPlots.cutMassDMax.value) ||
                 (varUtils.invMassD0Bar > cfgQaPlots.cutMassDMin.value &&
                  varUtils.invMassD0Bar < cfgQaPlots.cutMassDMax.value))) {
              if (track.sign() > 0) {
                registry.fill(HIST("hMassD0Proton"), varUtils.ptReso, varUtils.invMassReso - varUtils.invMassD0);
              } else {
                registry.fill(HIST("hMassD0Proton"), varUtils.ptReso, varUtils.invMassReso - varUtils.invMassD0Bar);
              }
            }
          }
        } // end of DType switch
        // fill track table
        if (!selectedTracks.count(track.globalIndex())) {
          rowTrkReduced(track.globalIndex(),
                        indexHfReducedCollision,
                        track.px(), track.py(), track.pz(), track.sign(),
                        track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr(),
                        track.tofNSigmaPi(), track.tofNSigmaKa(), track.tofNSigmaPr(),
                        track.hasTOF(), track.hasTPC(), track.itsNCls(), track.tpcNClsCrossedRows(), track.tpcChi2NCl());
          selectedTracks[track.globalIndex()] = rowTrkReduced.lastIndex();
        }
        fillHfCandD = true;
        if constexpr (DoMc) {
          auto indexHfCandCharm = rowCandDmesReduced.lastIndex() + 1;
          fillMcRecoInfoDTrack<DType>(particlesMc, candD, track, tracks, indexHfCandCharm, selectedTracks[track.globalIndex()], pdg, registry, rowMcRecTrkReduced);
        }
      } // end of loop on bachelor tracks
    } // end of do tracks
    if constexpr (DoGammas) {
      for (const auto& gamma : bachelorV0s) {
        auto trackPos = tracksIU.rawIteratorAt(gamma.posTrackId());
        auto trackNeg = tracksIU.rawIteratorAt(gamma.negTrackId());
        // Apply selsection
        auto gammaDauTracks = std::array{trackPos, trackNeg};
        HfResoCandidateV0 candGamma;
        if (!buildAndSelectGamma<BCs, Colls>(collision, prongIdsD, gammaDauTracks, cfgGammaCuts, candGamma, matCorr, vDriftMgr, rejectPairsWithCommonDaughter)) {
          continue;
        }
        // Get single track variables
        float chi2TpcDauGammaMax = -1.f;
        int nItsClsDauGammaMin = 8, nTpcCrossRowsDauGammaMin = 200;
        for (const auto& gammaTrack : gammaDauTracks) {
          if (gammaTrack.itsNCls() < nItsClsDauGammaMin) {
            nItsClsDauGammaMin = gammaTrack.itsNCls();
          }
          if (gammaTrack.tpcNClsCrossedRows() < nTpcCrossRowsDauGammaMin) {
            nTpcCrossRowsDauGammaMin = gammaTrack.tpcNClsCrossedRows();
          }
          if (gammaTrack.tpcChi2NCl() > chi2TpcDauGammaMax) {
            chi2TpcDauGammaMax = gammaTrack.tpcChi2NCl();
          }
        }
        // propagate gamma to primary vertex (if enabled)
        if (cfgGammaCuts.propagateGammatoPV.value) {
          std::array<float, 3> const pVecGammaOrig = {candGamma.mom[0], candGamma.mom[1], candGamma.mom[2]};
          std::array<float, 2> dcaInfo{};
          auto trackParGamma = o2::track::TrackPar(candGamma.pos, pVecGammaOrig, 0, true);
          trackParGamma.setPID(o2::track::PID::Photon);
          trackParGamma.setAbsCharge(0);
          o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParGamma, 2.f, matCorr, &dcaInfo);
          getPxPyPz(trackParGamma, candGamma.mom);
        }
        registry.fill(HIST("hAP"), candGamma.alpha, candGamma.qt);
        registry.fill(HIST("hV0Radius"), candGamma.radius);
        if constexpr (DType == DMesonType::D0) {
          varUtils.invMassReso = RecoDecay::m(std::array{varUtils.pVectorProng0, varUtils.pVectorProng1, candGamma.mom}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassKPlus, o2::constants::physics::MassGamma});
          varUtils.invMassResoBar = RecoDecay::m(std::array{varUtils.pVectorProng0, varUtils.pVectorProng1, candGamma.mom}, std::array{o2::constants::physics::MassKPlus, o2::constants::physics::MassPionCharged, o2::constants::physics::MassGamma});
          varUtils.ptReso = RecoDecay::pt(RecoDecay::sumOfVec(varUtils.pVectorProng0, varUtils.pVectorProng1, candGamma.mom));
          if (!cfgQaPlots.applyCutsForQaHistograms.value ||
              ((varUtils.invMassD0 > cfgQaPlots.cutMassDMin.value && varUtils.invMassD0 < cfgQaPlots.cutMassDMax.value) ||
               (varUtils.invMassD0Bar > cfgQaPlots.cutMassDMin.value && varUtils.invMassD0Bar < cfgQaPlots.cutMassDMax.value))) {
            registry.fill(HIST("hMassD0Gamma"), varUtils.ptReso, varUtils.invMassReso - varUtils.invMassD0);
            registry.fill(HIST("hMassD0Gamma"), varUtils.ptReso, varUtils.invMassResoBar - varUtils.invMassD0Bar);
          }
        }
        // fill V0 table --> use same for V0s and gammas
        // if information on V0 already stored, go to next V0
        if (!selectedGammas.count(gamma.globalIndex())) {
          rowCandV0Reduced(trackPos.globalIndex(), trackNeg.globalIndex(),
                           indexHfReducedCollision,
                           candGamma.pos[0], candGamma.pos[1], candGamma.pos[2],
                           candGamma.momPos[0], candGamma.momPos[1], candGamma.momPos[2],
                           candGamma.momNeg[0], candGamma.momNeg[1], candGamma.momNeg[2],
                           candGamma.cosPA,
                           candGamma.dcaV0ToPv,
                           nItsClsDauGammaMin, nTpcCrossRowsDauGammaMin, chi2TpcDauGammaMax,
                           candGamma.v0Type);
          selectedGammas[gamma.globalIndex()] = rowCandV0Reduced.lastIndex();
        }
        fillHfCandD = true;
      } // end of loop on V0 candidates
    } // end of do gammas
    // fill D candidate table
    if (fillHfCandD) { // fill candDplus table only once per D candidate, only if at least one V0 is found
      if constexpr (DType == DMesonType::Dplus) {
        rowCandDmesReduced(prongIdsD[0], prongIdsD[1], prongIdsD[2],
                           indexHfReducedCollision,
                           secondaryVertexD[0], secondaryVertexD[1], secondaryVertexD[2],
                           candD.pxProng0(), candD.pyProng0(), candD.pzProng0(),
                           candD.pxProng1(), candD.pyProng1(), candD.pzProng1(),
                           varUtils.pVectorProng2[0], varUtils.pVectorProng2[1], varUtils.pVectorProng2[2],
                           nItsClsDauMin, nTpcCrossRowsDauMin, chi2TpcDauMax, varUtils.signD);
        if constexpr (WithMl) {
          rowCandDmesMlReduced(bdtScores[0], bdtScores[1], bdtScores[2], bdtScores[3], bdtScores[4], bdtScores[5]);
        }
      } else if constexpr (DType == DMesonType::D0) {
        uint8_t selFlagD0 = {BIT(D0Sel::SelectedD0) | BIT(D0Sel::SelectedD0Bar)};
        if (candD.isSelD0() < cfgDmesCuts.selectionFlagD0.value) {
          CLRBIT(selFlagD0, D0Sel::SelectedD0);
        }
        if (candD.isSelD0bar() < cfgDmesCuts.selectionFlagD0Bar.value) {
          CLRBIT(selFlagD0, D0Sel::SelectedD0Bar);
        }
        rowCandDmesReduced(prongIdsD[0], prongIdsD[1],
                           indexHfReducedCollision,
                           secondaryVertexD[0], secondaryVertexD[1], secondaryVertexD[2],
                           candD.pxProng0(), candD.pyProng0(), candD.pzProng0(),
                           candD.pxProng1(), candD.pyProng1(), candD.pzProng1(),
                           nItsClsDauMin, nTpcCrossRowsDauMin, chi2TpcDauMax,
                           selFlagD0);
        if constexpr (WithMl) {
          rowCandDmesMlReduced(bdtScores[0], bdtScores[1], bdtScores[2], bdtScores[3], bdtScores[4], bdtScores[5]);
        }
      } else if constexpr (DType == DMesonType::Dstar) {
        rowCandDmesReduced(prongIdsD[0], prongIdsD[1], prongIdsD[2],
                           indexHfReducedCollision,
                           secondaryVertexD[0], secondaryVertexD[1], secondaryVertexD[2],
                           candD.pxProng0(), candD.pyProng0(), candD.pzProng0(),
                           candD.pxProng1(), candD.pyProng1(), candD.pzProng1(),
                           varUtils.pVectorProng2[0], varUtils.pVectorProng2[1], varUtils.pVectorProng2[2],
                           nItsClsDauMin, nTpcCrossRowsDauMin, chi2TpcDauMax,
                           nItsClsSoftPi, nTpcCrossRowsSoftPi, chi2TpcSoftPi,
                           varUtils.signD);
        if constexpr (WithMl) {
          rowCandDmesMlReduced(bdtScores[0], bdtScores[1], bdtScores[2], bdtScores[3], bdtScores[4], bdtScores[5]);
        }
      }
      fillHfReducedCollision = true;
      if constexpr (DType == DMesonType::Dstar) {
        registry.fill(HIST("hMassVsPtDstarPaired"), candD.pt(), varUtils.invMassD - varUtils.invMassD0);
      } else if constexpr (DType == DMesonType::Dplus) {
        registry.fill(HIST("hMassVsPtDplusPaired"), candD.pt(), varUtils.invMassD);
      } else if constexpr (DType == DMesonType::D0) {
        if (candD.isSelD0() >= cfgDmesCuts.selectionFlagD0.value) {
          registry.fill(HIST("hMassVsPtD0Paired"), varUtils.ptD, varUtils.invMassD0);
        }
        if (candD.isSelD0bar() >= cfgDmesCuts.selectionFlagD0Bar.value) {
          registry.fill(HIST("hMassVsPtD0BarPaired"), varUtils.ptD, varUtils.invMassD0Bar);
        }
      }
    }
  } // candsD loop
  registry.fill(HIST("hEvents"), 1 + EventType::Processed);
  if (!fillHfReducedCollision) {
    registry.fill(HIST("hEvents"), 1 + EventType::NoDV0Selected);
    return;
  }
  registry.fill(HIST("hEvents"), 1 + EventType::DV0Selected);
  // fill collision table if it contains a DPi pair a minima
  rowCollisionReduced(collision.posX(), collision.posY(), collision.posZ(), collision.numContrib(), hfRejMap, bz);
} // end of runDataCreation function

// Function for derived data creation
/// \tparam dType is the D meson type (Dstar, Dplus or D0)
/// \param mcParticles is the MC particle table
/// \param collInfos is the reco collision table with MC info
/// \param mcCollisions is the MC collision table
/// \param hfEvSelMc is the HF event selection util object from MC
/// \param rejectCollisionsWithBadEvSel is the flag to not store collisions rejected by event selection
/// \param registry is the histogram registry
/// \param pdg is the O2DatabasePDG service
/// \param rowHfResoMcGenReduced is the MC gen reduced table
template <DMesonType DType, PairingType PairType, typename McParticles, typename McParticlesPerMcColl, typename CCs, typename CollPerMcColl, typename McCollisions, typename TableMcGenRed, typename BCsInfo>
void runMcGen(McParticles const& mcParticles,
              McParticlesPerMcColl const& mcParticlesPerMcCollision,
              CCs const& collInfos,
              CollPerMcColl const& colPerMcCollision,
              McCollisions const& mcCollisions,
              o2::hf_evsel::HfEventSelectionMc& hfEvSelMc,
              bool rejectCollisionsWithBadEvSel,
              o2::framework::HistogramRegistry& registry,
              o2::framework::Service<o2::framework::O2DatabasePDG> const& pdg,
              TableMcGenRed& rowHfResoMcGenReduced,
              BCsInfo const&)
{
  bool const doV0s = (PairType == PairingType::V0Only || PairType == PairingType::V0AndTrack);
  bool const doTracks = (PairType == PairingType::TrackOnly || PairType == PairingType::V0AndTrack);
  for (const auto& mcCollision : mcCollisions) {
    // Slice the particles table to get the particles for the current MC collision
    const auto mcParticlesPerMcColl = mcParticles.sliceBy(mcParticlesPerMcCollision, mcCollision.globalIndex());
    // Slice the collisions table to get the collision info for the current MC collision
    float centrality{-1.f};
    o2::hf_evsel::HfCollisionRejectionMask rejectionMask{};
    int const nSplitColl = 0;
    const auto collSlice = collInfos.sliceBy(colPerMcCollision, mcCollision.globalIndex());
    rejectionMask = hfEvSelMc.getHfMcCollisionRejectionMask<BCsInfo, o2::hf_centrality::CentralityEstimator::None>(mcCollision, collSlice, centrality);
    hfEvSelMc.fillHistograms<o2::hf_centrality::CentralityEstimator::None>(mcCollision, rejectionMask, nSplitColl);
    if (rejectCollisionsWithBadEvSel && rejectionMask != 0) {
      // at least one event selection not satisfied --> reject all gen particles from this collision
      continue;
    }
    for (const auto& particle : mcParticlesPerMcColl) {
      int8_t sign{0};
      int8_t flag{0};
      int8_t signD{0};
      int8_t signBach{0};
      int8_t origin{0};
      bool matchedReso{false}, matchedD{false}, matchedV0Tr{false};
      std::vector<int> idxBhadMothers{};
      if constexpr (DType == DMesonType::Dstar) {
        if (doV0s) {
          // D* K0s
          for (const auto& [decayChannelFlag, pdgCodeReso] : hf_decay::hf_cand_reso::particlesToDstarK0s) {
            matchedReso = RecoDecay::isMatchedMCGen<false, true>(mcParticlesPerMcColl, particle, pdgCodeReso, std::array{static_cast<int>(o2::constants::physics::Pdg::kDStar), +kK0}, true, &sign, 1);
            if (matchedReso) {
              flag = sign * decayChannelFlag;
              auto candV0MC = mcParticles.rawIteratorAt(particle.daughtersIds().back());
              matchedV0Tr = RecoDecay::isMatchedMCGen<false, true>(mcParticlesPerMcColl, candV0MC, kK0, std::array{+kPiPlus, -kPiPlus}, true, &signBach, 2);
              break;
            }
          }
        }
        if (doTracks && !matchedReso) {
          // D*+ pi-
          for (const auto& [decayChannelFlag, pdgCodeReso] : hf_decay::hf_cand_reso::particlesToDstarPi) {
            matchedReso = RecoDecay::isMatchedMCGen<false, true>(mcParticlesPerMcColl, particle, pdgCodeReso, std::array{static_cast<int>(o2::constants::physics::Pdg::kDStar), -static_cast<int>(kPiPlus)}, true, &sign, 1);
            if (matchedReso) {
              flag = sign * decayChannelFlag;
              matchedV0Tr = true;
              break;
            }
          }
        }
        if (matchedReso && matchedV0Tr) {
          auto candDstarMC = mcParticles.rawIteratorAt(particle.daughtersIds().front());
          matchedD = RecoDecay::isMatchedMCGen(mcParticlesPerMcColl, candDstarMC, o2::constants::physics::Pdg::kDStar, std::array{static_cast<int>(o2::constants::physics::Pdg::kD0), +static_cast<int>(kPiPlus)}, true, &signD, 1);
          if (matchedD) {
            auto candD0MC = mcParticles.rawIteratorAt(candDstarMC.daughtersIds().front());
            matchedD = RecoDecay::isMatchedMCGen(mcParticlesPerMcColl, candD0MC, o2::constants::physics::Pdg::kD0, std::array{-kKPlus, +kPiPlus}, true, &signD, 2);
          }
        }
      } else if constexpr (DType == DMesonType::Dplus) {
        if (doV0s) {
          // D+ K0s
          for (const auto& [decayChannelFlag, pdgCodeReso] : hf_decay::hf_cand_reso::particlesToDplusK0s) {
            matchedReso = RecoDecay::isMatchedMCGen<false, true>(mcParticlesPerMcColl, particle, pdgCodeReso, std::array{static_cast<int>(o2::constants::physics::Pdg::kDPlus), +kK0}, true, &sign, 1);
            if (matchedReso) {
              flag = sign * decayChannelFlag;
              auto candV0MC = mcParticles.rawIteratorAt(particle.daughtersIds().back());
              matchedV0Tr = RecoDecay::isMatchedMCGen<false, true>(mcParticlesPerMcColl, candV0MC, kK0, std::array{+kPiPlus, -kPiPlus}, true, &signBach, 2);
              break;
            }
          }
          if (!matchedReso) {
            // D+ lambda
            for (const auto& [decayChannelFlag, pdgCodeReso] : hf_decay::hf_cand_reso::particlesToDplusLambda) {
              matchedReso = RecoDecay::isMatchedMCGen<false, true>(mcParticlesPerMcColl, particle, pdgCodeReso, std::array{static_cast<int>(o2::constants::physics::Pdg::kDPlus), +kLambda0}, true, &sign, 1);
              if (matchedReso) {
                flag = sign * decayChannelFlag;
                auto candV0MC = mcParticles.rawIteratorAt(particle.daughtersIds().back());
                matchedV0Tr = RecoDecay::isMatchedMCGen<false, true>(mcParticlesPerMcColl, candV0MC, kLambda0, std::array{+kProton, -kPiPlus}, true, &signBach, 1);
                break;
              }
            }
          }
        }
        if (doTracks && !matchedReso) {
          // D+ pi-
          for (const auto& [decayChannelFlag, pdgCodeReso] : hf_decay::hf_cand_reso::particlesToDplusPi) {
            matchedReso = RecoDecay::isMatchedMCGen<false, true>(mcParticlesPerMcColl, particle, pdgCodeReso, std::array{static_cast<int>(o2::constants::physics::Pdg::kDPlus), -static_cast<int>(kPiPlus)}, true, &sign, 1);
            if (matchedReso) {
              flag = sign * decayChannelFlag;
              matchedV0Tr = true;
              break;
            }
          }
        }
        if (matchedReso && matchedV0Tr) {
          auto candDplusMC = mcParticles.rawIteratorAt(particle.daughtersIds().front());
          matchedD = RecoDecay::isMatchedMCGen(mcParticlesPerMcColl, candDplusMC, o2::constants::physics::Pdg::kDPlus, std::array{+kPiPlus, -kKPlus, +kPiPlus}, true, &signD, 2);
        }
      } else if constexpr (DType == DMesonType::D0) {
        if (doV0s) {
          // D0 Lambda
          for (const auto& [decayChannelFlag, pdgCodeReso] : hf_decay::hf_cand_reso::particlesToD0Lambda) {
            matchedReso = RecoDecay::isMatchedMCGen<false, true>(mcParticlesPerMcColl, particle, pdgCodeReso, std::array{static_cast<int>(o2::constants::physics::Pdg::kD0), +kLambda0}, true, &sign, 1);
            if (matchedReso) {
              flag = sign * decayChannelFlag;
              auto candV0MC = mcParticles.rawIteratorAt(particle.daughtersIds().back());
              matchedV0Tr = RecoDecay::isMatchedMCGen<false, true>(mcParticlesPerMcColl, candV0MC, kLambda0, std::array{+kProton, -kPiPlus}, true, &signBach, 1);
              break;
            }
          }
        }
        if (doTracks && !matchedReso) {
          // D0 pi+
          for (const auto& [decayChannelFlag, pdgCodeReso] : hf_decay::hf_cand_reso::particlesToD0Pi) {
            matchedReso = RecoDecay::isMatchedMCGen<false, true>(mcParticlesPerMcColl, particle, pdgCodeReso, std::array{static_cast<int>(o2::constants::physics::Pdg::kD0), +static_cast<int>(kPiPlus)}, true, &sign, 1);
            if (matchedReso) {
              flag = sign * decayChannelFlag;
              matchedV0Tr = true;
              break;
            }
          }
          // D0 K+
          if (!matchedReso) {
            for (const auto& [decayChannelFlag, pdgCodeReso] : hf_decay::hf_cand_reso::particlesToD0Kplus) {
              matchedReso = RecoDecay::isMatchedMCGen<false, true>(mcParticlesPerMcColl, particle, pdgCodeReso, std::array{static_cast<int>(o2::constants::physics::Pdg::kD0), +static_cast<int>(kKPlus)}, true, &sign, 1);
              if (matchedReso) {
                flag = sign * decayChannelFlag;
                matchedV0Tr = true;
                break;
              }
            }
          }
        }
        if (matchedReso && matchedV0Tr) {
          auto candD0MC = mcParticles.rawIteratorAt(particle.daughtersIds().front());
          matchedD = RecoDecay::isMatchedMCGen(mcParticlesPerMcColl, candD0MC, o2::constants::physics::Pdg::kD0, std::array{-kKPlus, +kPiPlus}, true, &signD, 2);
        }
      }
      if (matchedReso && matchedD && matchedV0Tr) {
        origin = RecoDecay::getCharmHadronOrigin(mcParticlesPerMcColl, particle, false, &idxBhadMothers);
        registry.fill(HIST("hMCGenOrigin"), origin);
        auto ptParticle = particle.pt();
        auto invMassGen = computeInvMassGen(mcParticles, particle.globalIndex(), pdg);
        auto yParticle = RecoDecay::y(particle.pVector(), invMassGen);
        auto etaParticle = particle.eta();

        std::array<float, 2> ptProngs{};
        std::array<float, 2> yProngs{};
        std::array<float, 2> etaProngs{};
        int counter = 0;
        for (const auto& daught : particle.template daughters_as<McParticles>()) {
          ptProngs[counter] = daught.pt();
          etaProngs[counter] = daught.eta();
          yProngs[counter] = RecoDecay::y(daught.pVector(), pdg->Mass(daught.pdgCode()));
          counter++;
        }
        registry.fill(HIST("hMCGenCounter"), flag, ptParticle);
        rowHfResoMcGenReduced(flag, origin, ptParticle, yParticle, etaParticle,
                              ptProngs[0], yProngs[0], etaProngs[0],
                              ptProngs[1], yProngs[1], etaProngs[1],
                              invMassGen, rejectionMask);
      }
    }
  }
}
} // namespace hf_charm_reso
} // namespace o2::analysis

#endif // PWGHF_D2H_CORE_DATACREATIONCHARMRESO_H_
