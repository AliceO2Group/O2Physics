// Copyright 2019-2025 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file pairHistManager.h
/// \brief histogram manager for pair tasks
/// \author anton.riedel@tum.de, TU München, anton.riedel@tum.de

#ifndef PWGCF_FEMTO_CORE_PAIRHISTMANAGER_H_
#define PWGCF_FEMTO_CORE_PAIRHISTMANAGER_H_

#include "PWGCF/Femto/Core/femtoUtils.h"
#include "PWGCF/Femto/Core/histManager.h"
#include "PWGCF/Femto/Core/modes.h"

#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/Logger.h>

#include <Math/GenVector/Boost.h>
#include <Math/Vector4D.h> // IWYU pragma: keep (do not replace with Math/Vector4Dfwd.h)
#include <Math/Vector4Dfwd.h>

#include <array>
#include <cmath>
#include <map>
#include <set>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace o2::analysis::femto
{
namespace pairhistmanager
{
// enum for pair histograms
enum PairHist {
  // standard 1D
  kKstar,
  kKt,
  kMt,
  kMinv,
  // standard 2D
  kPt1VsPt2,
  kPt1VsKstar,
  kPt2VsKstar,
  kPt1VsKt,
  kPt2VsKt,
  kPt1VsMt,
  kPt2VsMt,
  kKstarVsKt,
  kKstarVsMt,
  kKstarVsMult,
  kKstarVsCent,
  // 2D with mass
  kKstarVsMass1,
  kKstarVsMass2,
  kMass1VsMass2,
  // 2D with invariant mass
  kKstarVsMinv,
  kPt1VsMinv,
  kPt2VsMinv,
  // higher dimensions
  kKstarVsMtVsMult,
  kKstarVsMtVsMultVsCent,
  // higher dimensions with pt
  kKstarVsMtVsPt1VsPt2,
  kKstarVsMtVsPt1VsPt2VsMult,
  kKstarVsMtVsPt1VsPt2VsMultVsCent,
  // higher dimensions with mass
  kKstarVsMtVsMass1VsMass2,
  kKstarVsMtVsMass1VsMass2VsMult,
  kKstarVsMtVsMass1VsMass2VsMultVsCent,
  // higher dimension with pt and mass
  kKstarVsMtVsMass1VsMass2VsPt1VsPt2,
  kKstarVsMtVsMass1VsMass2VsPt1VsPt2VsMult,
  kKstarVsMtVsMass1VsMass2VsPt1VsPt2VsMultVsCent,
  // higher dimensions with pt and invariant mass
  kKstarVsMtVsMinvVsPt1VsPt2,
  kKstarVsMtVsMinvVsPt1VsPt2VsMult,
  kKstarVsMtVsMinvVsPt1VsPt2VsMultVsCent,
  // dalitz plots
  kDalitz, // between a track and pos/neg daughter of another particle
  // mc
  kTrueKstarVsKstar,
  kTrueKtVsKt,
  kTrueMtVsMt,
  kTrueMinvVsMinv,
  kTrueMultVsMult,
  kTrueCentVsCent,

  // mixing qa
  kSeNpart1VsNpart2,                         // number of particles 1 vs number of particles 2 in each same event
  kMeMixingDepth,                            // mixing depth
  kMeNpart1VsNpart2,                         // number of particles 1 vs number of particles 2 in each mixed event
  kMeNpart1,                                 // number of unique particles 1 in each mixing bin
  kMeNpart2,                                 // number of unique particles 2 in each mixing bin
  kMeVtz1VsMult1VsCent1VsVtz2VsMult2VsCent2, // correlation of event properties in each mixing bin

  kPairHistogramLast
};

enum MixingPolicy {
  kVtxMult,
  kVtxCent,
  kVtxMultCent,
  kMixingPolicyLast
};

// Mixing configurables
struct ConfMixing : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("Mixing");
  o2::framework::ConfigurableAxis multBins{"multBins", {o2::framework::VARIABLE_WIDTH, 0.0f, 4.0f, 8.0f, 12.0f, 16.0f, 20.0f, 24.0f, 28.0f, 32.0f, 36.0f, 40.0f, 44.0f, 48.0f, 52.0f, 56.0f, 60.0f, 64.0f, 68.0f, 72.0f, 76.0f, 80.0f, 84.0f, 88.0f, 92.0f, 96.0f, 100.0f, 200.0f}, "Mixing bins - multiplicity"};
  o2::framework::ConfigurableAxis centBins{"centBins", {o2::framework::VARIABLE_WIDTH, 0.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.0f}, "Mixing bins - centrality"};
  o2::framework::ConfigurableAxis vtxBins{"vtxBins", {o2::framework::VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  o2::framework::Configurable<int> depth{"depth", 5, "Number of events for mixing"};
  o2::framework::Configurable<int> policy{"policy", 0, "Binning policy for mixing (alywas in combination with z-vertex) -> 0: multiplicity, -> 1: centrality, -> 2: both"};
  o2::framework::Configurable<bool> sameSpecies{"sameSpecies", false, "Enable if particle 1 and particle 2 are the same"};
  o2::framework::Configurable<int> seed{"seed", -1, "Seed to randomize particle 1 and particle 2 (if they are identical). Set to negative value to deactivate. Set to 0 to generate unique seed in time."};
  o2::framework::Configurable<bool> enablePairCorrelationQa{"enablePairCorrelationQa", true, "Enable pair-level correlation QA (same-event + mixed-event)"};
  o2::framework::Configurable<bool> enableEventMixingQa{"enableEventMixingQa", false, "Enable QA of event properties used in event mixing (vtx, multiplicity, centrality)"};
  o2::framework::ConfigurableAxis particleBinning{"particleBinning", {50, -0.5f, 49.5f}, "Binning for particle number correlation in pairs"};
};

struct ConfPairBinning : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("PairBinning");
  o2::framework::Configurable<bool> usePdgMass{"usePdgMass", true, "Use PDF masses for 4-vectors. If false, use reconstructed mass (if available)"};
  o2::framework::Configurable<bool> plot1D{"plot1D", true, "Enable 1D histograms"};
  o2::framework::Configurable<bool> plot2D{"plot2D", true, "Enable 2D histograms"};
  o2::framework::Configurable<bool> plotKstarVsMtVsMult{"plotKstarVsMtVsMult", false, "Enable 3D histogram (Kstar Vs Mt Vs Mult)"};
  o2::framework::Configurable<bool> plotKstarVsMtVsMultVsCent{"plotKstarVsMtVsMultVsCent", false, "Enable 4D histogram (Kstar Vs Mt Vs Mult Vs Cent)"};
  o2::framework::Configurable<bool> plotKstarVsMtVsPt1VsPt2{"plotKstarVsMtVsPt1VsPt2", false, "Enable 4D histogram (Kstar Vs Mt Vs Pt1 Vs Pt2 Vs Mult)"};
  o2::framework::Configurable<bool> plotKstarVsMtVsPt1VsPt2VsMult{"plotKstarVsMtVsPt1VsPt2VsMult", false, "Enable 5D histogram (Kstar Vs Mt Vs Pt1 Vs Pt2 Vs Mult)"};
  o2::framework::Configurable<bool> plotKstarVsMtVsPt1VsPt2VsMultVsCent{"plotKstarVsMtVsPt1VsPt2VsMultVsCent", false, "Enable 6D histogram (Kstar Vs Mt Vs Pt1 Vs Pt2 Vs Mult Vs Cent)"};
  o2::framework::Configurable<bool> plotKstarVsMtVsMass1VsMass2{"plotKstarVsMtVsMass1VsMass2", false, "Enable 4D histogram (Kstar Vs Mt Vs Mass1 Vs Mass2)"};
  o2::framework::Configurable<bool> plotKstarVsMtVsMass1VsMass2VsMult{"plotKstarVsMtVsMass1VsMass2VsMult", false, "Enable 5D histogram (Kstar Vs Mt Vs Mass1 Vs Mass2 Vs Mult)"};
  o2::framework::Configurable<bool> plotKstarVsMtVsMass1VsMass2VsMultVsCent{"plotKstarVsMtVsMass1VsMass2VsMultVsCent", false, "Enable 6D histogram (Kstar Vs Mt Vs Mass1 Vs Mass2 Vs Mult Vs Cent)"};
  o2::framework::Configurable<bool> plotKstarVsMtVsMass1VsMass2VsPt1VsPt2{"plotKstarVsMtVsMass1VsMass2VsPt1VsPt2", false, "Enable 6D histogram (Kstar Vs Mt Vs Pt1 Vs Pt2 Vs Mass1 Vs Mass2)"};
  o2::framework::Configurable<bool> plotKstarVsMtVsMass1VsMass2VsPt1VsPt2VsMult{"plotKstarVsMtVsMass1VsMass2VsPt1VsPt2VsMult", false, "Enable 7D histogram (Kstar Vs Mt Vs Pt1 Vs Pt2 Vs Mass1 Vs Mass2 Vs Mult)"};
  o2::framework::Configurable<bool> plotKstarVsMtVsMass1VsMass2VsPt1VsPt2VsMultVsCent{"plotKstarVsMtVsMass1VsMass2VsPt1VsPt2VsMultVsCent", false, "Enable 8D histogram (Kstar Vs Mt Vs Pt1 Vs Pt2 Vs Mass1 Vs Mass2 Vs Mult Vs Cent)"};
  o2::framework::Configurable<bool> plotDalitz{"plotDalitz", false, "Enable dalitz plot"};
  o2::framework::ConfigurableAxis kstar{"kstar", {{600, 0, 6}}, "kstar"};
  o2::framework::ConfigurableAxis kt{"kt", {{600, 0, 6}}, "kt"};
  o2::framework::ConfigurableAxis mt{"mt", {{500, 0.8, 5.8}}, "mt"};
  o2::framework::ConfigurableAxis multiplicity{"multiplicity", {{50, 0, 200}}, "multiplicity"};
  o2::framework::ConfigurableAxis centrality{"centrality", {{10, 0, 100}}, "centrality (mult. percentile)"};
  o2::framework::ConfigurableAxis pt1{"pt1", {{100, 0, 6}}, "Pt binning for particle 1"};
  o2::framework::ConfigurableAxis pt2{"pt2", {{100, 0, 6}}, "Pt binning for particle 2"};
  o2::framework::ConfigurableAxis mass1{"mass1", {{100, 0, 2}}, "Mass binning for particle 1 (if particle has mass getter, otherwise PDG mass)"};
  o2::framework::ConfigurableAxis mass2{"mass2", {{100, 0, 2}}, "Mass binning for particle 2 (if particle has mass getter, otherwise PDG mass)"};
  o2::framework::ConfigurableAxis massInv{"massInv", {{100, 0, 2}}, "Invariant Mass binning"};
  o2::framework::ConfigurableAxis dalitzMtot{"dalitzMtot", {{100, 0, 10}}, "Total invariant mass squared binning in darlitz plot"};
  o2::framework::ConfigurableAxis dalitzM12{"dalitzM12", {{100, 0, 10}}, "Mass12 binning of darlitz plot"};
  o2::framework::ConfigurableAxis dalitzM13{"dalitzM13", {{100, 0, 10}}, "Mass13 binning of darlitz plot"};
  o2::framework::Configurable<int> transverseMassType{"transverseMassType", static_cast<int>(modes::TransverseMassType::kAveragePdgMass), "Type of transverse mass (0-> Average Pdg Mass, 1-> Reduced Pdg Mass, 2-> Mt from combined 4 vector)"};
};

struct ConfPairCuts : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("PairCuts");
  o2::framework::Configurable<float> kstarMax{"kstarMax", -1, "Maximal kstar (set to -1 to deactivate)"};
  o2::framework::Configurable<float> kstarMin{"kstarMin", -1, "Minimal kstar (set to -1 to deactivate)"};
  o2::framework::Configurable<float> ktMax{"ktMax", -1, "Maximal kt (set to -1 to deactivate)"};
  o2::framework::Configurable<float> ktMin{"ktMin", -1, "Minimal kt (set to -1 to deactivate)"};
  o2::framework::Configurable<float> mtMax{"mtMax", -1, "Maximal mt (set to -1 to deactivate)"};
  o2::framework::Configurable<float> mtMin{"mtMin", -1, "Minimal mt (set to -1 to deactivate)"};
  o2::framework::Configurable<float> massInvMin{"massInvMin", -1, "Minimal invariant mass (set to -1 to deactivate)"};
  o2::framework::Configurable<float> massInvMax{"massInvMax", -1, "Maximal invariant mass (set to -1 to deactivate)"};
  o2::framework::Configurable<bool> mixOnlyCommonAncestor{"mixOnlyCommonAncestor", false, "Require pair to have common anchestor (in the same event)"};
  o2::framework::Configurable<bool> mixOnlyNonCommonAncestor{"mixOnlyNonCommonAncestor", false, "Require pair to have non-common anchestor (in the same event)"};
};

// the enum gives the correct index in the array
constexpr std::array<histmanager::HistInfo<PairHist>, kPairHistogramLast>
  HistTable = {
    {
      // 1D
      {kKstar, o2::framework::HistType::kTH1F, "hKstar", "k*; k* (GeV/#it{c}); Entries"},
      {kKt, o2::framework::HistType::kTH1F, "hKt", "transverse momentum; k_{T} (GeV/#it{c}); Entries"},
      {kMt, o2::framework::HistType::kTH1F, "hMt", "transverse mass; m_{T} (GeV/#it{c}^{2}); Entries"},
      {kMinv, o2::framework::HistType::kTH1F, "hMinv", "invariant mass; m_{Inv} (GeV/#it{c}^{2}); Entries"},
      // 2D
      {kPt1VsPt2, o2::framework::HistType::kTH2F, "hPt1VsPt2", " p_{T,1} vs p_{T,2}; p_{T,1} (GeV/#it{c}); p_{T,2} (GeV/#it{c})"},
      {kPt1VsKstar, o2::framework::HistType::kTH2F, "hPt1VsKstar", "p_{T,1} vs k*; p_{T,2} (GeV/#it{c}); k* (GeV/#it{c})"},
      {kPt2VsKstar, o2::framework::HistType::kTH2F, "hPt2VsKstar", "p_{T,2} vs k*; p_{T,2} (GeV/#it{c}); k* (GeV/#it{c})"},
      {kPt1VsKt, o2::framework::HistType::kTH2F, "hPt1VsKt", "p_{T,1} vs k_{T}; p_{T,1} (GeV/#it{c}); k_{T} (GeV/#it{c})"},
      {kPt2VsKt, o2::framework::HistType::kTH2F, "hPt2VsKt", "p_{T,2} vs k_{T}; p_{T,2} (GeV/#it{c}); k_{T} (GeV/#it{c})"},
      {kPt1VsMt, o2::framework::HistType::kTH2F, "hPt1VsMt", "p_{T,1} vs m_{T}; p_{T,1} (GeV/#it{c}); m_{T} (GeV/#it{c}^{2})"},
      {kPt2VsMt, o2::framework::HistType::kTH2F, "hPt2VsMt", "p_{T,2} vs m_{T}; p_{T,2} (GeV/#it{c}); m_{T} (GeV/#it{c}^{2})"},
      {kKstarVsKt, o2::framework::HistType::kTH2F, "hKstarVsKt", "k* vs k_{T}; k* (GeV/#it{c}); k_{T} (GeV/#it{c})"},
      {kKstarVsMt, o2::framework::HistType::kTH2F, "hKstarVsMt", "k* vs m_{T}; k* (GeV/#it{c}); m_{T} (GeV/#it{c}^{2})"},
      {kKstarVsCent, o2::framework::HistType::kTH2F, "hKstarVsCent", "k* vs Centrality (Mult. Percentile); k* (GeV/#it{c}); Centrality (%)"},
      {kKstarVsMult, o2::framework::HistType::kTH2F, "hKstarVsMult", "k* vs Multiplicity; k* (GeV/#it{c}); Multiplicity"},
      // 2D with mass
      {kKstarVsMass1, o2::framework::HistType::kTH2F, "hKstarVsMass1", "k* vs m_{1}; k* (GeV/#it{c}); m_{1} (GeV/#it{c}^{2})"},
      {kKstarVsMass2, o2::framework::HistType::kTH2F, "hKstarVsMass2", "k* vs m_{2}; k* (GeV/#it{c}); m_{2} (GeV/#it{c}^{2})"},
      {kMass1VsMass2, o2::framework::HistType::kTH2F, "hMass1VsMass2", "m_{1} vs m_{2}; m_{1} (GeV/#it{c}^{2}); m_{2} (GeV/#it{c}^{2})"},
      {kKstarVsMinv, o2::framework::HistType::kTH2F, "hKstarVsMinv", "k* vs m_{Inv}; k* (GeV/#it{c}); m_{Inv} (GeV/#it{c}^{2})"},
      {kPt1VsMinv, o2::framework::HistType::kTH2F, "hPt1VsMinv", "p_{T,1} vs m_{Inv}; p_{T,1} (GeV/#it{c}); m_{Inv} (GeV/#it{c}^{2})"},
      {kPt2VsMinv, o2::framework::HistType::kTH2F, "hPt2VsMinv", "p_{T,2} vs m_{Inv}; p_{T,2} (GeV/#it{c}); m_{Inv} (GeV/#it{c}^{2})"},
      // n-D
      {kKstarVsMtVsMult, o2::framework::HistType::kTHnSparseF, "hKstarVsMtVsMult", "k* vs m_{T} vs multiplicity; k* (GeV/#it{c}); m_{T} (GeV/#it{c}^{2}); Multiplicity;"},
      {kKstarVsMtVsMultVsCent, o2::framework::HistType::kTHnSparseF, "hKstarVsMtVsMultVsCent", "k* vs m_{T} vs multiplicity vs centrality; k* (GeV/#it{c}); m_{T} (GeV/#it{c}^{2}); Multiplicity; Centrality (%);"},
      // n-D with pt
      {kKstarVsMtVsPt1VsPt2, o2::framework::HistType::kTHnSparseF, "hKstarVsMtVsPt1VsPt2", "k* vs m_{T} vs p_{T,1} vs p_{T,2}; k* (GeV/#it{c}); m_{T} (GeV/#it{c}^{2}); p_{T,1} (GeV/#it{c}); p_{T,2} (GeV/#it{c});"},
      {kKstarVsMtVsPt1VsPt2VsMult, o2::framework::HistType::kTHnSparseF, "hKstarVsMtVsPt1VsPt2VsMult", "k* vs m_{T} vs p_{T,1} vs p_{T,2} vs multiplicity; k* (GeV/#it{c}); m_{T} (GeV/#it{c}^{2}); p_{T,1} (GeV/#it{c}); p_{T,2} (GeV/#it{c}); Multiplicity;"},
      {kKstarVsMtVsPt1VsPt2VsMultVsCent, o2::framework::HistType::kTHnSparseF, "hKstarVsMtVsPt1VsPt2VsMultVsCent", "k* vs m_{T} vs p_{T,1} vs p_{T,2} vs multiplicity vs centrality; k* (GeV/#it{c}); m_{T} (GeV/#it{c}^{2}); p_{T,1} (GeV/#it{c}); p_{T,2} (GeV/#it{c}); Multiplicity; Centrality;"},
      // n-D with mass
      {kKstarVsMtVsMass1VsMass2, o2::framework::HistType::kTHnSparseF, "hKstarVsMtVsMass1VsMass2", "k* vs m_{T} vs m_{1} vs m_{2}; k* (GeV/#it{c}); m_{T} (GeV/#it{c}^{2}); m_{1} (GeV/#it{c}^{2}); m_{1} (GeV/#it{c}^{2});"},
      {kKstarVsMtVsMass1VsMass2VsMult, o2::framework::HistType::kTHnSparseF, "hKstarVsMtVsMass1VsMass2VsMult", "k* vs m_{T} vs m_{1} vs m_{2} vs multiplicity; k* (GeV/#it{c}); m_{T} (GeV/#it{c}^{2}); m_{1} (GeV/#it{c}^{2}); m_{1} (GeV/#it{c}^{2}); Multiplicity;"},
      {kKstarVsMtVsMass1VsMass2VsMultVsCent, o2::framework::HistType::kTHnSparseF, "hKstarVsMtVsMass1VsMass2VsMultVsCent", "k* vs m_{T} vs m_{1} vs m_{2} vs multiplicity vs centrality; k* (GeV/#it{c}); m_{T} (GeV/#it{c}^{2}); m_{1} (GeV/#it{c}^{2}); m_{1} (GeV/#it{c}^{2}); Multiplicity; Centrality (%);"},
      // n-D with pt and mass
      {kKstarVsMtVsMass1VsMass2VsPt1VsPt2, o2::framework::HistType::kTHnSparseF, "hKstarVsMtVsMass1VsMass2VsPt1VsPt2", "k* vs m_{T} vs m_{1} vs m_{2} vs p_{T,1} vs p_{T,2}; k* (GeV/#it{c}); m_{T} (GeV/#it{c}^{2}); m_{1} (GeV/#it{c}^{2}); m_{1} (GeV/#it{c}^{2}); p_{T,1} (GeV/#it{c}); p_{T,2} (GeV/#it{c});"},
      {kKstarVsMtVsMass1VsMass2VsPt1VsPt2VsMult, o2::framework::HistType::kTHnSparseF, "hKstarVsMtVsMass1VsMass2VsPt1VsPt2VsMult", "k* vs m_{T} vs m_{1} vs m_{2} vs p_{T,1} vs p_{T,2} vs multiplicity; k* (GeV/#it{c}); m_{T} (GeV/#it{c}^{2}); m_{1} (GeV/#it{c}^{2}); m_{1} (GeV/#it{c}^{2}); p_{T,1} (GeV/#it{c}); p_{T,2} (GeV/#it{c}); Multiplicity;"},
      {kKstarVsMtVsMass1VsMass2VsPt1VsPt2VsMultVsCent, o2::framework::HistType::kTHnSparseF, "hKstarVsMtVsMass1VsMass2VsPt1VsPt2VsMultVsCent", "k* vs m_{T} vs m_{1} vs m_{2} vs p_{T,1} vs p_{T,2} vs multiplicity vs centrality; k* (GeV/#it{c}); m_{T} (GeV/#it{c}^{2}); m_{1} (GeV/#it{c}^{2}); m_{1} (GeV/#it{c}^{2}); p_{T,1} (GeV/#it{c}); p_{T,2} (GeV/#it{c}); Multiplicity; Centrality (%);"},
      {kDalitz, o2::framework::HistType::kTHnSparseF, "hDalitz", "Dalitz plot; k* (GeV/#it{c}); m^{2}_{123} (GeV/#it{c}^{2})^{2}; m^{2}_{12} (GeV/#it{c}^{2})^{2}; m^{2}_{13} (GeV/#it{c}^{2})^{2};"},
      {kTrueKstarVsKstar, o2::framework::HistType::kTH2F, "hTrueKstarVsKstar", "k*_{True} vs k*; k*_{True} (GeV/#it{c});  k* (GeV/#it{c})"},
      {kTrueKtVsKt, o2::framework::HistType::kTH2F, "hTrueKtVsKt", "k_{T,True} vs k_{T}; k_{T,True} (GeV/#it{c});  k_{T} (GeV/#it{c})"},
      {kTrueMtVsMt, o2::framework::HistType::kTH2F, "hTrueMtVsMt", "m_{T,True} vs m_{T}; m_{T,True} (GeV/#it{c}^{2}); m_{T} (GeV/#it{c}^{2})"},
      {kTrueMinvVsMinv, o2::framework::HistType::kTH2F, "hTrueMinvVsMinv", "m_{Inv,True} vs m_{Inv}; m_{Inv,True} (GeV/#it{c}^{2}); m_{Inv} (GeV/#it{c}^{2})"},
      {kTrueMultVsMult, o2::framework::HistType::kTH2F, "hTrueMultVsMult", "Multiplicity_{True} vs Multiplicity; Multiplicity_{True} ;  Multiplicity"},
      {kTrueCentVsCent, o2::framework::HistType::kTH2F, "hTrueCentVsCent", "Centrality_{True} vs Centrality; Centrality_{True} (%); Centrality (%)"},
      {kSeNpart1VsNpart2, o2::framework::HistType::kTH2F, "hSeNpart1VsNpart2", "# particle 1 vs # particle 2 in each same event; # partilce 1; # particle 2"},
      {kMeMixingDepth, o2::framework::HistType::kTH1F, "hMeMixingDepth", "Mixing Depth; Mixing Depth ; Entries"},
      {kMeNpart1VsNpart2, o2::framework::HistType::kTH2F, "hMeNpart1VsNpart2", "# particle 1 vs # particle 2 in each mixed event pair; # partilce 1; # particle 2"},
      {kMeNpart1, o2::framework::HistType::kTH1F, "hMeNpart1", "# particle 1 in each mixing bin; # partilce 1; Entries"},
      {kMeNpart2, o2::framework::HistType::kTH1F, "hMeNpart2", "# particle 2 in each mixing bin; # particle 2; Entries"},
      {kMeVtz1VsMult1VsCent1VsVtz2VsMult2VsCent2, o2::framework::HistType::kTHnSparseF, "hVtz1VsMult1VsCent1VsVtz2VsMult2VsCent2", "Mixing bins; V_{z,1} (cm); multiplicity_{1}; centrality_{1} (%); V_{z,2} (cm); multiplicity_{2}; centrality_{2} (%)"},
    }};

#define PAIR_HIST_ANALYSIS_MAP(confAnalysis, confMixing)                                                                                                                                                                     \
  {kKstar, {confAnalysis.kstar}},                                                                                                                                                                                            \
    {kKt, {confAnalysis.kt}},                                                                                                                                                                                                \
    {kMt, {confAnalysis.mt}},                                                                                                                                                                                                \
    {kMinv, {confAnalysis.massInv}},                                                                                                                                                                                         \
    {kPt1VsPt2, {confAnalysis.pt1, confAnalysis.pt2}},                                                                                                                                                                       \
    {kPt1VsKstar, {confAnalysis.pt1, confAnalysis.kstar}},                                                                                                                                                                   \
    {kPt2VsKstar, {confAnalysis.pt2, confAnalysis.kstar}},                                                                                                                                                                   \
    {kPt1VsKt, {confAnalysis.pt1, confAnalysis.kt}},                                                                                                                                                                         \
    {kPt2VsKt, {confAnalysis.pt2, confAnalysis.kt}},                                                                                                                                                                         \
    {kPt1VsMt, {confAnalysis.pt1, confAnalysis.mt}},                                                                                                                                                                         \
    {kPt2VsMt, {confAnalysis.pt2, confAnalysis.mt}},                                                                                                                                                                         \
    {kKstarVsKt, {confAnalysis.kstar, confAnalysis.kt}},                                                                                                                                                                     \
    {kKstarVsMt, {confAnalysis.kstar, confAnalysis.mt}},                                                                                                                                                                     \
    {kKstarVsMult, {confAnalysis.kstar, confAnalysis.multiplicity}},                                                                                                                                                         \
    {kKstarVsCent, {confAnalysis.kstar, confAnalysis.centrality}},                                                                                                                                                           \
    {kKstarVsMass1, {confAnalysis.kstar, confAnalysis.mass1}},                                                                                                                                                               \
    {kKstarVsMass2, {confAnalysis.kstar, confAnalysis.mass2}},                                                                                                                                                               \
    {kMass1VsMass2, {confAnalysis.mass1, confAnalysis.mass2}},                                                                                                                                                               \
    {kKstarVsMinv, {confAnalysis.kstar, confAnalysis.massInv}},                                                                                                                                                              \
    {kPt1VsMinv, {confAnalysis.pt1, confAnalysis.massInv}},                                                                                                                                                                  \
    {kPt2VsMinv, {confAnalysis.pt2, confAnalysis.massInv}},                                                                                                                                                                  \
    {kKstarVsMtVsMult, {confAnalysis.kstar, confAnalysis.mt, confAnalysis.multiplicity}},                                                                                                                                    \
    {kKstarVsMtVsMultVsCent, {confAnalysis.kstar, confAnalysis.mt, confAnalysis.multiplicity, confAnalysis.centrality}},                                                                                                     \
    {kKstarVsMtVsPt1VsPt2, {confAnalysis.kstar, confAnalysis.mt, confAnalysis.pt1, confAnalysis.pt2}},                                                                                                                       \
    {kKstarVsMtVsPt1VsPt2VsMult, {confAnalysis.kstar, confAnalysis.mt, confAnalysis.pt1, confAnalysis.pt2, confAnalysis.multiplicity}},                                                                                      \
    {kKstarVsMtVsPt1VsPt2VsMultVsCent, {confAnalysis.kstar, confAnalysis.mt, confAnalysis.pt1, confAnalysis.pt2, confAnalysis.multiplicity, confAnalysis.centrality}},                                                       \
    {kKstarVsMtVsMass1VsMass2, {confAnalysis.kstar, confAnalysis.mt, confAnalysis.mass1, confAnalysis.mass2}},                                                                                                               \
    {kKstarVsMtVsMass1VsMass2VsMult, {confAnalysis.kstar, confAnalysis.mt, confAnalysis.mass1, confAnalysis.mass2, confAnalysis.multiplicity}},                                                                              \
    {kKstarVsMtVsMass1VsMass2VsMultVsCent, {confAnalysis.kstar, confAnalysis.mt, confAnalysis.mass1, confAnalysis.mass2, confAnalysis.multiplicity, confAnalysis.centrality}},                                               \
    {kKstarVsMtVsMass1VsMass2VsPt1VsPt2, {confAnalysis.kstar, confAnalysis.mt, confAnalysis.mass1, confAnalysis.mass2, confAnalysis.pt1, confAnalysis.pt2}},                                                                 \
    {kKstarVsMtVsMass1VsMass2VsPt1VsPt2VsMult, {confAnalysis.kstar, confAnalysis.mt, confAnalysis.mass1, confAnalysis.mass2, confAnalysis.pt1, confAnalysis.pt2, confAnalysis.multiplicity}},                                \
    {kKstarVsMtVsMass1VsMass2VsPt1VsPt2VsMultVsCent, {confAnalysis.kstar, confAnalysis.mt, confAnalysis.mass1, confAnalysis.mass2, confAnalysis.pt1, confAnalysis.pt2, confAnalysis.multiplicity, confAnalysis.centrality}}, \
    {kDalitz, {confAnalysis.kstar, confAnalysis.dalitzMtot, confAnalysis.dalitzM12, confAnalysis.dalitzM13}},                                                                                                                \
    {kSeNpart1VsNpart2, {confMixing.particleBinning, confMixing.particleBinning}},                                                                                                                                           \
    {kMeMixingDepth, {confMixing.particleBinning}},                                                                                                                                                                          \
    {kMeNpart1VsNpart2, {confMixing.particleBinning, confMixing.particleBinning}},                                                                                                                                           \
    {kMeNpart1, {confMixing.particleBinning}},                                                                                                                                                                               \
    {kMeNpart2, {confMixing.particleBinning}},                                                                                                                                                                               \
    {kMeVtz1VsMult1VsCent1VsVtz2VsMult2VsCent2, {confMixing.vtxBins, confMixing.multBins, confMixing.centBins, confMixing.vtxBins, confMixing.multBins, confMixing.centBins}},

#define PAIR_HIST_MC_MAP(conf)                                 \
  {kTrueKstarVsKstar, {conf.kstar, conf.kstar}},               \
    {kTrueKtVsKt, {conf.kt, conf.kt}},                         \
    {kTrueMtVsMt, {conf.mt, conf.mt}},                         \
    {kTrueMinvVsMinv, {conf.massInv, conf.massInv}},           \
    {kTrueMultVsMult, {conf.multiplicity, conf.multiplicity}}, \
    {kTrueCentVsCent, {conf.centrality, conf.centrality}},

template <typename T1, typename T2>
auto makePairHistSpecMap(T1 const& confPairBinning, T2 const& confMixing)
{
  return std::map<PairHist, std::vector<o2::framework::AxisSpec>>{
    PAIR_HIST_ANALYSIS_MAP(confPairBinning, confMixing)};
};

template <typename T1, typename T2>
auto makePairMcHistSpecMap(T1 const& confPairBinning, T2 const& confMixing)
{
  return std::map<PairHist, std::vector<o2::framework::AxisSpec>>{
    PAIR_HIST_ANALYSIS_MAP(confPairBinning, confMixing)
      PAIR_HIST_MC_MAP(confPairBinning)};
};

#undef PAIR_HIST_ANALYSIS_MAP
#undef PAIR_HIST_MC_MAP

constexpr char PrefixTrackTrackSe[] = "TrackTrack/SE/";
constexpr char PrefixTrackTrackMe[] = "TrackTrack/ME/";

constexpr char PrefixTrackV0Se[] = "TrackV0/SE/";
constexpr char PrefixTrackV0Me[] = "TrackV0/ME/";

constexpr char PrefixV0V0Se[] = "V0V0/SE/";
constexpr char PrefixV0V0Me[] = "V0V0/ME/";

constexpr char PrefixTrackResonanceSe[] = "TrackResonance/SE/";
constexpr char PrefixTrackResonanceMe[] = "TrackResonance/ME/";

constexpr char PrefixTrackCascadeSe[] = "TrackCascade/SE/";
constexpr char PrefixTrackCascadeMe[] = "TrackCascade/ME/";

constexpr char PrefixTrackKinkSe[] = "TrackKink/SE/";
constexpr char PrefixTrackKinkMe[] = "TrackKink/ME/";

constexpr std::string_view AnalysisDir = "Analysis/";
constexpr std::string_view MixingQaDir = "MixingQA/";
constexpr std::string_view McDir = "MC/";

template <const char* prefix,
          modes::Particle particleType1,
          modes::Particle particleType2>
class PairHistManager
{
 public:
  PairHistManager() = default;
  ~PairHistManager() = default;

  template <modes::Mode mode, typename T1, typename T2, typename T3>
  void init(o2::framework::HistogramRegistry* registry,
            std::map<PairHist, std::vector<o2::framework::AxisSpec>> const& Specs,
            T1 const& ConfPairBinning,
            T2 const& ConfPairCuts,
            T3 const& ConfMixing)
  {
    mHistogramRegistry = registry;

    mUsePdgMass = ConfPairBinning.usePdgMass.value;

    // flags for histograms
    mPlot1d = ConfPairBinning.plot1D.value;
    mPlot2d = ConfPairBinning.plot2D.value;
    mPlotKstarVsMtVsMult = ConfPairBinning.plotKstarVsMtVsMult.value;
    mPlotKstarVsMtVsMultVsCent = ConfPairBinning.plotKstarVsMtVsMultVsCent.value;

    mPlotKstarVsMtVsPt1VsPt2 = ConfPairBinning.plotKstarVsMtVsPt1VsPt2.value;
    mPlotKstarVsMtVsPt1VsPt2VsMult = ConfPairBinning.plotKstarVsMtVsPt1VsPt2VsMult.value;
    mPlotKstarVsMtVsPt1VsPt2VsMultVsCent = ConfPairBinning.plotKstarVsMtVsPt1VsPt2VsMultVsCent.value;

    mPlotKstarVsMtVsMass1VsMass2 = ConfPairBinning.plotKstarVsMtVsMass1VsMass2.value;
    mPlotKstarVsMtVsMass1VsMass2VsMult = ConfPairBinning.plotKstarVsMtVsMass1VsMass2VsMult.value;
    mPlotKstarVsMtVsMass1VsMass2VsMultVsCent = ConfPairBinning.plotKstarVsMtVsMass1VsMass2VsMultVsCent.value;

    mPlotKstarVsMtVsMass1VsMass2VsPt1VsPt2 = ConfPairBinning.plotKstarVsMtVsMass1VsMass2VsPt1VsPt2.value;
    mPlotKstarVsMtVsMass1VsMass2VsPt1VsPt2VsMult = ConfPairBinning.plotKstarVsMtVsMass1VsMass2VsPt1VsPt2VsMult.value;
    mPlotKstarVsMtVsMass1VsMass2VsPt1VsPt2VsMultVsCent = ConfPairBinning.plotKstarVsMtVsMass1VsMass2VsPt1VsPt2VsMultVsCent.value;

    mPlotDalitz = ConfPairBinning.plotDalitz.value;

    // transverse mass type
    mMtType = static_cast<modes::TransverseMassType>(ConfPairBinning.transverseMassType.value);

    // values for cuts
    mKstarMin = ConfPairCuts.kstarMin.value;
    mKstarMax = ConfPairCuts.kstarMax.value;
    mKtMin = ConfPairCuts.ktMin.value;
    mKtMax = ConfPairCuts.ktMax.value;
    mMtMin = ConfPairCuts.mtMin.value;
    mMtMax = ConfPairCuts.mtMax.value;
    mMassInvMin = ConfPairCuts.massInvMin.value;
    mMassInvMax = ConfPairCuts.massInvMax.value;

    mPairCorrelationQa = ConfMixing.enablePairCorrelationQa.value;
    mEventMixingQa = ConfMixing.enableEventMixingQa.value;

    if constexpr (isFlagSet(mode, modes::Mode::kAnalysis)) {
      initAnalysis(Specs);
    }

    if constexpr (isFlagSet(mode, modes::Mode::kMc)) {
      initMc(Specs);
    }

    if constexpr (isFlagSet(mode, modes::Mode::kSe)) {
      initSeMixingQa(Specs);
    }

    if constexpr (isFlagSet(mode, modes::Mode::kMe)) {
      initMeMixingQa(Specs);
    }
  }

  void setMass(int PdgParticle1, int PdgParticle2)
  {
    mPdgMass1 = utils::getPdgMass(PdgParticle1);
    mPdgMass2 = utils::getPdgMass(PdgParticle2);
  }

  void setMass(int PdgParticle1, int PdgPosDauParticle1, int PdgNegDauParticle1, int PdgParticle2, int PdgPosDauParticle2, int PdgNegDauParticle2)
  {
    mPdgMass1 = utils::getPdgMass(PdgParticle1);
    mPdgMassPosDau1 = utils::getPdgMass(PdgPosDauParticle1);
    mPdgMassNegDau1 = utils::getPdgMass(PdgNegDauParticle1);
    mPdgMass2 = utils::getPdgMass(PdgParticle2);
    mPdgMassPosDau2 = utils::getPdgMass(PdgPosDauParticle2);
    mPdgMassNegDau2 = utils::getPdgMass(PdgNegDauParticle2);
  }
  void setCharge(int chargeAbsParticle1, int chargeAbsParticle2)
  {
    // the pt stored is actually as pt/z for tracks, so in case of particles with z > 1, we have to rescale the pt (this is so far only for He3 the case)
    mAbsCharge1 = std::abs(chargeAbsParticle1);
    mAbsCharge2 = std::abs(chargeAbsParticle2);
  }

  template <typename T1, typename T2, typename T3>
  void setPair(T1 const& particle1, T2 const& particle2, T3 const& trackTable)
  {
    // if one of the particles has a mass getter (like lambda), we cache the value for the filling later
    // otherwise set it to the pdg mass
    if constexpr (utils::HasMass<T1>) {
      mRecoMass1 = particle1.mass();
    } else {
      mRecoMass1 = mPdgMass1;
    }
    if constexpr (utils::HasMass<T2>) {
      mRecoMass2 = particle2.mass();
    } else {
      mRecoMass2 = mPdgMass2;
    }

    // get mass for 4-vectors
    double mass1 = 0.f;
    double mass2 = 0.f;
    if (mUsePdgMass) {
      mass1 = mPdgMass1;
      mass2 = mPdgMass2;
    } else {
      mass1 = mRecoMass1;
      mass2 = mRecoMass2;
    }

    // pt in track table is calculated from 1/signedPt from the original track table
    // in case of He with Z=2, we have to rescale the pt with the absolute charge
    mParticle1 = ROOT::Math::PtEtaPhiMVector(mAbsCharge1 * particle1.pt(), particle1.eta(), particle1.phi(), mass1);
    mParticle2 = ROOT::Math::PtEtaPhiMVector(mAbsCharge2 * particle2.pt(), particle2.eta(), particle2.phi(), mass2);

    // set kT
    mKt = getKt(mParticle1, mParticle2);

    // set mT
    mMt = getMt(mParticle1, mParticle2);

    // set Minv
    mMassInv = getMinv(mParticle1, mParticle2);

    // set kstar
    mKstar = getKstar(mParticle1, mParticle2);

    if (mPlotDalitz) {
      if constexpr (modes::isEqual(particleType1, modes::Particle::kTrack) && modes::isEqual(particleType2, modes::Particle::kV0)) {
        auto posDaughter = trackTable.rawIteratorAt(particle2.posDauId() - trackTable.offset());
        auto negDaughter = trackTable.rawIteratorAt(particle2.negDauId() - trackTable.offset());
        ROOT::Math::PtEtaPhiMVector posDau4v = ROOT::Math::PtEtaPhiMVector(posDaughter.pt(), posDaughter.eta(), posDaughter.phi(), mPdgMassPosDau2);
        ROOT::Math::PtEtaPhiMVector negDau4v = ROOT::Math::PtEtaPhiMVector(negDaughter.pt(), negDaughter.eta(), negDaughter.phi(), mPdgMassNegDau2);
        mMassTot2 = (mParticle1 + posDau4v + negDau4v).M2();
        mMass12 = (mParticle1 + posDau4v).M2();
        mMass13 = (mParticle1 + negDau4v).M2();
      }
    }
  }

  template <typename T1, typename T2, typename T3, typename T4>
  void setPair(T1 const& particle1, T2 const& particle2, T3 const& trackTable, T4 const& col)
  {
    setPair(particle1, particle2, trackTable);
    mMult = col.mult();
    mCent = col.cent();
  }

  template <typename T1, typename T2, typename T3, typename T4, typename T5>
  void setPair(T1 const& particle1, T2 const& particle2, T3 const& trackTable, T4 const& col1, T5 const& col2)
  {
    setPair(particle1, particle2, trackTable);
    mMult = 0.5f * (col1.mult() + col2.mult()); // if mixing with multiplicity, should be in the same mixing bin
    mCent = 0.5f * (col1.cent() + col2.cent()); // if mixing with centrality, should be in the same mixing bin
  }

  template <typename T1, typename T2, typename T3>
  void setPairMc(T1 const& particle1, T2 const& particle2, const T3& /*mcParticles*/)
  {
    if (!particle1.has_fMcParticle() || !particle2.has_fMcParticle()) {
      mHasMcPair = false;
      return;
    }
    mHasMcPair = true;

    auto mcParticle1 = particle1.template fMcParticle_as<T3>();
    auto mcParticle2 = particle2.template fMcParticle_as<T3>();

    mTrueParticle1 = ROOT::Math::PtEtaPhiMVector(mAbsCharge1 * mcParticle1.pt(), mcParticle1.eta(), mcParticle1.phi(), mPdgMass1);
    mTrueParticle2 = ROOT::Math::PtEtaPhiMVector(mAbsCharge2 * mcParticle2.pt(), mcParticle2.eta(), mcParticle2.phi(), mPdgMass2);

    // compute true kinematics
    mTrueKt = getKt(mTrueParticle1, mTrueParticle2);
    mTrueMt = getMt(mTrueParticle1, mTrueParticle2);
    mTrueMt = getMinv(mTrueParticle1, mTrueParticle2);
    mTrueKstar = getKstar(mTrueParticle1, mTrueParticle2);
  }

  template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
  void setPairMc(T1 const& particle1, T2 const& particle2, T3 const& trackTable, T4 const& mcParticles, T5 const& col, T6 const& /*mcCols*/)
  {
    setPair(particle1, particle2, trackTable, col);
    setPairMc(particle1, particle2, mcParticles);
    if (!col.has_fMcCol()) {
      mHasMcCol = false;
      return;
    }
    mHasMcCol = true;
    auto mcCol = col.template fMcCol_as<T6>();
    mTrueMult = mcCol.mult();
    mTrueCent = mcCol.cent();
  }

  template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7>
  void setPairMc(T1 const& particle1, T2 const& particle2, T3 const& trackTable, T4 const& mcParticles, T5 const& col1, T6 const& col2, T7 const& /*mcCols*/)
  {
    setPair(particle1, particle2, trackTable, col1, col2);
    setPairMc(particle1, particle2, mcParticles);
    if (!col1.has_fMcCol() || !col2.has_fMcCol()) {
      mHasMcCol = false;
      return;
    }
    mHasMcCol = true;
    auto mcCol1 = col1.template fMcCol_as<T7>();
    auto mcCol2 = col2.template fMcCol_as<T7>();
    mTrueMult = 0.5f * (mcCol1.mult() + mcCol2.mult());
    mTrueCent = 0.5f * (mcCol1.cent() + mcCol2.cent());
  }

  bool checkPairCuts() const
  {
    return (!(mKstarMin > 0.f) || mKstar > mKstarMin) &&
           (!(mKstarMax > 0.f) || mKstar < mKstarMax) &&
           (!(mKtMin > 0.f) || mKt > mKtMin) &&
           (!(mKtMax > 0.f) || mKt < mKtMax) &&
           (!(mMtMin > 0.f) || mMt > mMtMin) &&
           (!(mMtMax > 0.f) || mMt < mMtMax) &&
           (!(mMassInvMin > 0.f) || mMassInv > mMassInvMin) &&
           (!(mMassInvMax > 0.f) || mMassInv < mMassInvMax)

      ;
  }

  template <modes::Mode mode>
  void fill()
  {
    if constexpr (isFlagSet(mode, modes::Mode::kAnalysis)) {
      fillAnalysis();
    }
    if constexpr (isFlagSet(mode, modes::Mode::kMc)) {
      fillMc();
    }
  }

  template <typename T1, typename T2>
  void trackParticlesPerMixingBin(T1 const& particle1, T2 const& particle2)
  {
    mUniqueParticles1PerMixingBin.emplace(particle1.globalIndex(), 0);
    mUniqueParticles1PerMixingBin.at(particle1.globalIndex())++;
    mUniqueParticles2PerMixingBin.emplace(particle2.globalIndex(), 0);
    mUniqueParticles2PerMixingBin.at(particle2.globalIndex())++;
  }

  void resetTrackedParticlesPerMixingBin()
  {
    mUniqueParticles1PerMixingBin.clear();
    mUniqueParticles2PerMixingBin.clear();
  }

  template <typename T1, typename T2>
  void trackParticlesPerEvent(T1 const& particle1, T2 const& particle2)
  {
    mUniqueParticles1PerEvent.insert(particle1.globalIndex());
    mUniqueParticles2PerEvent.insert(particle2.globalIndex());
  }

  template <typename T1, typename T2>
  void fillMixingQaMe(T1 const& col1, T2 const& col2)
  {
    if (mEventMixingQa) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(MixingQaDir) + HIST(getHistName(kMeVtz1VsMult1VsCent1VsVtz2VsMult2VsCent2, HistTable)), col1.posZ(), col1.mult(), col1.cent(), col2.posZ(), col2.mult(), col2.cent());
    }
  }

  void resetTrackedParticlesPerEvent()
  {
    mUniqueParticles1PerEvent.clear();
    mUniqueParticles2PerEvent.clear();
  }

  void fillMixingQaSe()
  {
    if (mPairCorrelationQa) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(MixingQaDir) + HIST(getHistName(kSeNpart1VsNpart2, HistTable)), mUniqueParticles1PerEvent.size(), mUniqueParticles2PerEvent.size());
    }
  }

  void fillMixingQaMePerEvent()
  {
    if (mPairCorrelationQa) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(MixingQaDir) + HIST(getHistName(kMeNpart1VsNpart2, HistTable)), mUniqueParticles1PerEvent.size(), mUniqueParticles2PerEvent.size());
    }
  }

  void fillMixingQaMePerMixingBin(int windowSize)
  {
    if (!mPairCorrelationQa || windowSize <= 0) {
      return;
    }
    mHistogramRegistry->fill(HIST(prefix) + HIST(MixingQaDir) + HIST(getHistName(kMeMixingDepth, HistTable)), windowSize);
    for (const auto& [_, nPart1] : mUniqueParticles1PerMixingBin) {
      mHistogramRegistry->fill(
        HIST(prefix) + HIST(MixingQaDir) + HIST(getHistName(kMeNpart1, HistTable)),
        nPart1);
    }
    for (const auto& [_, nPart2] : mUniqueParticles2PerMixingBin) {
      mHistogramRegistry->fill(
        HIST(prefix) + HIST(MixingQaDir) + HIST(getHistName(kMeNpart2, HistTable)),
        nPart2);
    }
  }

  float getKstar() const { return mKstar; }

 private:
  void initAnalysis(std::map<PairHist, std::vector<o2::framework::AxisSpec>> const& Specs)
  {
    std::string analysisDir = std::string(prefix) + std::string(AnalysisDir);
    if (mPlot1d) {
      mHistogramRegistry->add(analysisDir + getHistNameV2(kKstar, HistTable), getHistDesc(kKstar, HistTable), getHistType(kKstar, HistTable), {Specs.at(kKstar)});
      mHistogramRegistry->add(analysisDir + getHistNameV2(kKt, HistTable), getHistDesc(kKt, HistTable), getHistType(kKt, HistTable), {Specs.at(kKt)});
      mHistogramRegistry->add(analysisDir + getHistNameV2(kMt, HistTable), getHistDesc(kMt, HistTable), getHistType(kMt, HistTable), {Specs.at(kMt)});
      mHistogramRegistry->add(analysisDir + getHistNameV2(kMinv, HistTable), getHistDesc(kMinv, HistTable), getHistType(kMinv, HistTable), {Specs.at(kMinv)});
    }
    if (mPlot2d) {
      mHistogramRegistry->add(analysisDir + getHistNameV2(kPt1VsPt2, HistTable), getHistDesc(kPt1VsPt2, HistTable), getHistType(kPt1VsPt2, HistTable), {Specs.at(kPt1VsPt2)});
      mHistogramRegistry->add(analysisDir + getHistNameV2(kPt1VsKstar, HistTable), getHistDesc(kPt1VsKstar, HistTable), getHistType(kPt1VsKstar, HistTable), {Specs.at(kPt1VsKstar)});
      mHistogramRegistry->add(analysisDir + getHistNameV2(kPt2VsKstar, HistTable), getHistDesc(kPt2VsKstar, HistTable), getHistType(kPt2VsKstar, HistTable), {Specs.at(kPt2VsKstar)});
      mHistogramRegistry->add(analysisDir + getHistNameV2(kPt1VsKt, HistTable), getHistDesc(kPt1VsKt, HistTable), getHistType(kPt1VsKt, HistTable), {Specs.at(kPt1VsKt)});
      mHistogramRegistry->add(analysisDir + getHistNameV2(kPt2VsKt, HistTable), getHistDesc(kPt2VsKt, HistTable), getHistType(kPt2VsKt, HistTable), {Specs.at(kPt2VsKt)});
      mHistogramRegistry->add(analysisDir + getHistNameV2(kPt1VsMt, HistTable), getHistDesc(kPt1VsMt, HistTable), getHistType(kPt1VsMt, HistTable), {Specs.at(kPt1VsMt)});
      mHistogramRegistry->add(analysisDir + getHistNameV2(kPt2VsMt, HistTable), getHistDesc(kPt2VsMt, HistTable), getHistType(kPt2VsMt, HistTable), {Specs.at(kPt2VsMt)});

      mHistogramRegistry->add(analysisDir + getHistNameV2(kKstarVsKt, HistTable), getHistDesc(kKstarVsKt, HistTable), getHistType(kKstarVsKt, HistTable), {Specs.at(kKstarVsKt)});
      mHistogramRegistry->add(analysisDir + getHistNameV2(kKstarVsMt, HistTable), getHistDesc(kKstarVsMt, HistTable), getHistType(kKstarVsMt, HistTable), {Specs.at(kKstarVsMt)});
      mHistogramRegistry->add(analysisDir + getHistNameV2(kKstarVsMult, HistTable), getHistDesc(kKstarVsMult, HistTable), getHistType(kKstarVsMult, HistTable), {Specs.at(kKstarVsMult)});
      mHistogramRegistry->add(analysisDir + getHistNameV2(kKstarVsCent, HistTable), getHistDesc(kKstarVsCent, HistTable), getHistType(kKstarVsCent, HistTable), {Specs.at(kKstarVsCent)});
      mHistogramRegistry->add(analysisDir + getHistNameV2(kKstarVsMass1, HistTable), getHistDesc(kKstarVsMass1, HistTable), getHistType(kKstarVsMass1, HistTable), {Specs.at(kKstarVsMass1)});
      mHistogramRegistry->add(analysisDir + getHistNameV2(kKstarVsMass2, HistTable), getHistDesc(kKstarVsMass2, HistTable), getHistType(kKstarVsMass2, HistTable), {Specs.at(kKstarVsMass2)});
      mHistogramRegistry->add(analysisDir + getHistNameV2(kMass1VsMass2, HistTable), getHistDesc(kMass1VsMass2, HistTable), getHistType(kMass1VsMass2, HistTable), {Specs.at(kMass1VsMass2)});

      mHistogramRegistry->add(analysisDir + getHistNameV2(kKstarVsMinv, HistTable), getHistDesc(kKstarVsMinv, HistTable), getHistType(kKstarVsMinv, HistTable), {Specs.at(kKstarVsMinv)});
      mHistogramRegistry->add(analysisDir + getHistNameV2(kPt1VsMinv, HistTable), getHistDesc(kPt1VsMinv, HistTable), getHistType(kPt1VsMinv, HistTable), {Specs.at(kPt1VsMinv)});
      mHistogramRegistry->add(analysisDir + getHistNameV2(kPt2VsMinv, HistTable), getHistDesc(kPt2VsMinv, HistTable), getHistType(kPt2VsMinv, HistTable), {Specs.at(kPt2VsMinv)});
    }

    // higher dimensional histograms
    if (mPlotKstarVsMtVsMult) {
      mHistogramRegistry->add(analysisDir + getHistNameV2(kKstarVsMtVsMult, HistTable), getHistDesc(kKstarVsMtVsMult, HistTable), getHistType(kKstarVsMtVsMult, HistTable), {Specs.at(kKstarVsMtVsMult)});
    }
    if (mPlotKstarVsMtVsMultVsCent) {
      mHistogramRegistry->add(analysisDir + getHistNameV2(kKstarVsMtVsMultVsCent, HistTable), getHistDesc(kKstarVsMtVsMultVsCent, HistTable), getHistType(kKstarVsMtVsMultVsCent, HistTable), {Specs.at(kKstarVsMtVsMultVsCent)});
    }
    // add pt
    if (mPlotKstarVsMtVsPt1VsPt2) {
      mHistogramRegistry->add(analysisDir + getHistNameV2(kKstarVsMtVsPt1VsPt2, HistTable), getHistDesc(kKstarVsMtVsPt1VsPt2, HistTable), getHistType(kKstarVsMtVsPt1VsPt2, HistTable), {Specs.at(kKstarVsMtVsPt1VsPt2)});
    }
    if (mPlotKstarVsMtVsPt1VsPt2VsMult) {
      mHistogramRegistry->add(analysisDir + getHistNameV2(kKstarVsMtVsPt1VsPt2VsMult, HistTable), getHistDesc(kKstarVsMtVsPt1VsPt2VsMult, HistTable), getHistType(kKstarVsMtVsPt1VsPt2VsMult, HistTable), {Specs.at(kKstarVsMtVsPt1VsPt2VsMult)});
    }
    if (mPlotKstarVsMtVsPt1VsPt2VsMultVsCent) {
      mHistogramRegistry->add(analysisDir + getHistNameV2(kKstarVsMtVsPt1VsPt2VsMultVsCent, HistTable), getHistDesc(kKstarVsMtVsPt1VsPt2VsMultVsCent, HistTable), getHistType(kKstarVsMtVsPt1VsPt2VsMultVsCent, HistTable), {Specs.at(kKstarVsMtVsPt1VsPt2VsMultVsCent)});
    }
    // add mass
    if (mPlotKstarVsMtVsMass1VsMass2) {
      mHistogramRegistry->add(analysisDir + getHistNameV2(kKstarVsMtVsMass1VsMass2, HistTable), getHistDesc(kKstarVsMtVsMass1VsMass2, HistTable), getHistType(kKstarVsMtVsMass1VsMass2, HistTable), {Specs.at(kKstarVsMtVsMass1VsMass2)});
    }
    if (mPlotKstarVsMtVsMass1VsMass2VsMult) {
      mHistogramRegistry->add(analysisDir + getHistNameV2(kKstarVsMtVsMass1VsMass2VsMult, HistTable), getHistDesc(kKstarVsMtVsMass1VsMass2VsMult, HistTable), getHistType(kKstarVsMtVsMass1VsMass2VsMult, HistTable), {Specs.at(kKstarVsMtVsMass1VsMass2VsMult)});
    }
    if (mPlotKstarVsMtVsMass1VsMass2VsMultVsCent) {
      mHistogramRegistry->add(analysisDir + getHistNameV2(kKstarVsMtVsMass1VsMass2VsMultVsCent, HistTable), getHistDesc(kKstarVsMtVsMass1VsMass2VsMultVsCent, HistTable), getHistType(kKstarVsMtVsMass1VsMass2VsMultVsCent, HistTable), {Specs.at(kKstarVsMtVsMass1VsMass2VsMultVsCent)});
    }
    if (mPlotKstarVsMtVsMass1VsMass2VsPt1VsPt2) {
      mHistogramRegistry->add(analysisDir + getHistNameV2(kKstarVsMtVsMass1VsMass2VsPt1VsPt2, HistTable), getHistDesc(kKstarVsMtVsMass1VsMass2VsPt1VsPt2, HistTable), getHistType(kKstarVsMtVsMass1VsMass2VsPt1VsPt2, HistTable), {Specs.at(kKstarVsMtVsMass1VsMass2VsPt1VsPt2)});
    }
    if (mPlotKstarVsMtVsMass1VsMass2VsPt1VsPt2VsMult) {
      mHistogramRegistry->add(analysisDir + getHistNameV2(kKstarVsMtVsMass1VsMass2VsPt1VsPt2VsMult, HistTable), getHistDesc(kKstarVsMtVsMass1VsMass2VsPt1VsPt2VsMult, HistTable), getHistType(kKstarVsMtVsMass1VsMass2VsPt1VsPt2VsMult, HistTable), {Specs.at(kKstarVsMtVsMass1VsMass2VsPt1VsPt2VsMult)});
    }
    if (mPlotKstarVsMtVsMass1VsMass2VsPt1VsPt2VsMultVsCent) {
      mHistogramRegistry->add(analysisDir + getHistNameV2(kKstarVsMtVsMass1VsMass2VsPt1VsPt2VsMultVsCent, HistTable), getHistDesc(kKstarVsMtVsMass1VsMass2VsPt1VsPt2VsMultVsCent, HistTable), getHistType(kKstarVsMtVsMass1VsMass2VsPt1VsPt2VsMultVsCent, HistTable), {Specs.at(kKstarVsMtVsMass1VsMass2VsPt1VsPt2VsMultVsCent)});
    }
    if (mPlotDalitz) {
      mHistogramRegistry->add(analysisDir + getHistNameV2(kDalitz, HistTable), getHistDesc(kDalitz, HistTable), getHistType(kDalitz, HistTable), {Specs.at(kDalitz)});
    }
  }

  void initMc(std::map<PairHist, std::vector<o2::framework::AxisSpec>> const& Specs)
  {
    std::string mcDir = std::string(prefix) + std::string(McDir);
    mHistogramRegistry->add(mcDir + getHistNameV2(kTrueKstarVsKstar, HistTable), getHistDesc(kTrueKstarVsKstar, HistTable), getHistType(kTrueKstarVsKstar, HistTable), {Specs.at(kTrueKstarVsKstar)});
    mHistogramRegistry->add(mcDir + getHistNameV2(kTrueKtVsKt, HistTable), getHistDesc(kTrueKtVsKt, HistTable), getHistType(kTrueKtVsKt, HistTable), {Specs.at(kTrueKtVsKt)});
    mHistogramRegistry->add(mcDir + getHistNameV2(kTrueMtVsMt, HistTable), getHistDesc(kTrueMtVsMt, HistTable), getHistType(kTrueMtVsMt, HistTable), {Specs.at(kTrueMtVsMt)});
    mHistogramRegistry->add(mcDir + getHistNameV2(kTrueMinvVsMinv, HistTable), getHistDesc(kTrueMinvVsMinv, HistTable), getHistType(kTrueMinvVsMinv, HistTable), {Specs.at(kTrueMinvVsMinv)});
    mHistogramRegistry->add(mcDir + getHistNameV2(kTrueMultVsMult, HistTable), getHistDesc(kTrueMultVsMult, HistTable), getHistType(kTrueMultVsMult, HistTable), {Specs.at(kTrueMultVsMult)});
    mHistogramRegistry->add(mcDir + getHistNameV2(kTrueCentVsCent, HistTable), getHistDesc(kTrueCentVsCent, HistTable), getHistType(kTrueCentVsCent, HistTable), {Specs.at(kTrueCentVsCent)});
  }

  void initSeMixingQa(std::map<PairHist, std::vector<o2::framework::AxisSpec>> const& Specs)
  {
    std::string mcDir = std::string(prefix) + std::string(MixingQaDir);
    if (mPairCorrelationQa) {
      mHistogramRegistry->add(mcDir + getHistNameV2(kSeNpart1VsNpart2, HistTable), getHistDesc(kSeNpart1VsNpart2, HistTable), getHistType(kSeNpart1VsNpart2, HistTable), {Specs.at(kSeNpart1VsNpart2)});
    }
  }

  void initMeMixingQa(std::map<PairHist, std::vector<o2::framework::AxisSpec>> const& Specs)
  {
    std::string mcDir = std::string(prefix) + std::string(MixingQaDir);
    if (mPairCorrelationQa) {
      mHistogramRegistry->add(mcDir + getHistNameV2(kMeMixingDepth, HistTable), getHistDesc(kMeMixingDepth, HistTable), getHistType(kMeMixingDepth, HistTable), {Specs.at(kMeMixingDepth)});
      mHistogramRegistry->add(mcDir + getHistNameV2(kMeNpart1VsNpart2, HistTable), getHistDesc(kMeNpart1VsNpart2, HistTable), getHistType(kMeNpart1VsNpart2, HistTable), {Specs.at(kMeNpart1VsNpart2)});
      mHistogramRegistry->add(mcDir + getHistNameV2(kMeNpart1, HistTable), getHistDesc(kMeNpart1, HistTable), getHistType(kMeNpart1, HistTable), {Specs.at(kMeNpart1)});
      mHistogramRegistry->add(mcDir + getHistNameV2(kMeNpart2, HistTable), getHistDesc(kMeNpart2, HistTable), getHistType(kMeNpart2, HistTable), {Specs.at(kMeNpart2)});
    }
    if (mEventMixingQa) {
      mHistogramRegistry->add(mcDir + getHistNameV2(kMeVtz1VsMult1VsCent1VsVtz2VsMult2VsCent2, HistTable), getHistDesc(kMeVtz1VsMult1VsCent1VsVtz2VsMult2VsCent2, HistTable), getHistType(kMeVtz1VsMult1VsCent1VsVtz2VsMult2VsCent2, HistTable), {Specs.at(kMeVtz1VsMult1VsCent1VsVtz2VsMult2VsCent2)});
    }
  }

  void fillAnalysis()
  {
    if (mPlot1d) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kKstar, HistTable)), mKstar);
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kMt, HistTable)), mMt);
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kKt, HistTable)), mKt);
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kMinv, HistTable)), mMassInv);
    }
    if (mPlot2d) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kPt1VsPt2, HistTable)), mParticle1.Pt(), mParticle2.Pt());
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kPt1VsKstar, HistTable)), mParticle1.Pt(), mKstar);
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kPt1VsMt, HistTable)), mParticle1.Pt(), mMt);
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kPt1VsKt, HistTable)), mParticle1.Pt(), mKt);
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kPt2VsKstar, HistTable)), mParticle2.Pt(), mKstar);
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kPt2VsMt, HistTable)), mParticle2.Pt(), mMt);
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kPt2VsKt, HistTable)), mParticle2.Pt(), mKt);
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kKstarVsKt, HistTable)), mKstar, mKt);
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kKstarVsMt, HistTable)), mKstar, mMt);
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kKstarVsMult, HistTable)), mKstar, mMult);
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kKstarVsCent, HistTable)), mKstar, mCent);

      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kKstarVsMass1, HistTable)), mKstar, mRecoMass1);
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kKstarVsMass2, HistTable)), mKstar, mRecoMass2);
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kMass1VsMass2, HistTable)), mRecoMass1, mRecoMass2);

      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kKstarVsMinv, HistTable)), mKstar, mMassInv);
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kPt1VsMinv, HistTable)), mParticle1.Pt(), mMassInv);
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kPt2VsMinv, HistTable)), mParticle2.Pt(), mMassInv);
    }

    // n-D histograms are only filled if enabled
    // if "mass" getter does not exist for particle, it will be just set to 0
    // the user has to make sure that in this case the bin number of this dimension is set to 1
    if (mPlotKstarVsMtVsMult) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kKstarVsMtVsMult, HistTable)), mKstar, mMt, mMult);
    }
    if (mPlotKstarVsMtVsMultVsCent) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kKstarVsMtVsMultVsCent, HistTable)), mKstar, mMt, mMult, mCent);
    }
    if (mPlotKstarVsMtVsPt1VsPt2) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kKstarVsMtVsPt1VsPt2, HistTable)), mKstar, mMt, mParticle1.Pt(), mParticle2.Pt());
    }
    if (mPlotKstarVsMtVsPt1VsPt2VsMult) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kKstarVsMtVsPt1VsPt2VsMult, HistTable)), mKstar, mMt, mParticle1.Pt(), mParticle2.Pt(), mMult);
    }
    if (mPlotKstarVsMtVsPt1VsPt2VsMultVsCent) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kKstarVsMtVsPt1VsPt2VsMultVsCent, HistTable)), mKstar, mMt, mParticle1.Pt(), mParticle2.Pt(), mMult, mCent);
    }
    if (mPlotKstarVsMtVsMass1VsMass2) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kKstarVsMtVsMass1VsMass2, HistTable)), mKstar, mMt, mRecoMass1, mRecoMass2);
    }
    if (mPlotKstarVsMtVsMass1VsMass2VsMult) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kKstarVsMtVsMass1VsMass2VsMult, HistTable)), mKstar, mMt, mRecoMass1, mRecoMass2, mMult);
    }
    if (mPlotKstarVsMtVsMass1VsMass2VsMultVsCent) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kKstarVsMtVsMass1VsMass2VsMultVsCent, HistTable)), mKstar, mMt, mRecoMass1, mRecoMass2, mMult, mCent);
    }
    if (mPlotKstarVsMtVsMass1VsMass2VsPt1VsPt2) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kKstarVsMtVsMass1VsMass2VsPt1VsPt2, HistTable)), mKstar, mMt, mRecoMass1, mRecoMass2, mParticle1.Pt(), mParticle2.Pt());
    }
    if (mPlotKstarVsMtVsMass1VsMass2VsPt1VsPt2VsMult) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kKstarVsMtVsMass1VsMass2VsPt1VsPt2VsMult, HistTable)), mKstar, mMt, mRecoMass1, mRecoMass2, mParticle1.Pt(), mParticle2.Pt(), mMult);
    }
    if (mPlotKstarVsMtVsMass1VsMass2VsPt1VsPt2VsMultVsCent) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kKstarVsMtVsMass1VsMass2VsPt1VsPt2VsMultVsCent, HistTable)), mKstar, mMt, mRecoMass1, mRecoMass2, mParticle1.Pt(), mParticle2.Pt(), mMult, mCent);
    }
    if (mPlotDalitz) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kDalitz, HistTable)), mKstar, mMassTot2, mMass12, mMass13);
    }
  }

  void fillMc()
  {
    if (mHasMcPair) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kTrueKstarVsKstar, HistTable)), mTrueKstar, mKstar);
      mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kTrueKtVsKt, HistTable)), mTrueKt, mKt);
      mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kTrueMtVsMt, HistTable)), mTrueMt, mMt);
      mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kTrueMinvVsMinv, HistTable)), mTrueMinv, mMassInv);
    }
    if (mHasMcCol) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kTrueMultVsMult, HistTable)), mTrueMult, mMult);
      mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kTrueCentVsCent, HistTable)), mTrueCent, mCent);
    }
  }

  float getKt(ROOT::Math::PtEtaPhiMVector const& part1, ROOT::Math::PtEtaPhiMVector const& part2)
  {
    auto sum = (part1 + part2);
    double kt = 0.5 * sum.Pt();
    return static_cast<float>(kt);
  }

  float getMt(ROOT::Math::PtEtaPhiMVector const& part1, ROOT::Math::PtEtaPhiMVector const& part2)
  {
    auto sum = part1 + part2;
    double mt = 0;
    double averageMass = 0;
    double reducedMass = 0;
    switch (mMtType) {
      case modes::TransverseMassType::kAveragePdgMass:
        averageMass = 0.5 * (part1.M() + part2.M());
        mt = std::hypot(0.5 * sum.Pt(), averageMass);
        break;
      case modes::TransverseMassType::kReducedPdgMass:
        reducedMass = 2. * (part1.M() * part2.M()) / (part1.M() + part2.M());
        mt = std::hypot(0.5 * sum.Pt(), reducedMass);
        break;
      case modes::TransverseMassType::kMt4Vector:
        mt = 0.5 * sum.Mt();
        break;
      default:
        LOG(fatal) << "Invalid transverse mass type, breaking...";
    }
    return static_cast<float>(mt);
  }

  float getMinv(ROOT::Math::PtEtaPhiMVector const& part1, ROOT::Math::PtEtaPhiMVector const& part2)
  {
    return static_cast<float>((part1 + part2).M());
  }

  float getKstar(ROOT::Math::PtEtaPhiMVector const& part1, ROOT::Math::PtEtaPhiMVector const& part2)
  {
    // compute pair momentum
    auto sum = part1 + part2;
    // Boost particle 1 to the pair rest frame (Prf) and calculate k* (would be equivalent using particle 2)
    // make a copy of particle 1
    auto particle1Prf = ROOT::Math::PtEtaPhiMVector(part1);
    // get lorentz boost into pair rest frame
    ROOT::Math::Boost boostPrf(sum.BoostToCM());
    // boost particle 1 into pair rest frame and calculate its momentum, which has the same value as k*
    return static_cast<float>(boostPrf(particle1Prf).P());
  }

  o2::framework::HistogramRegistry* mHistogramRegistry = nullptr;
  bool mUsePdgMass = true;
  double mPdgMass1 = 0.;
  double mPdgMassPosDau1 = 0;
  double mPdgMassNegDau1 = 0;
  double mPdgMass2 = 0.;
  double mPdgMassPosDau2 = 0;
  double mPdgMassNegDau2 = 0;

  modes::TransverseMassType mMtType = modes::TransverseMassType::kAveragePdgMass;

  int mAbsCharge1 = 1;
  int mAbsCharge2 = 1;
  ROOT::Math::PtEtaPhiMVector mParticle1{};
  ROOT::Math::PtEtaPhiMVector mParticle2{};
  float mRecoMass1 = 0.f;
  float mRecoMass2 = 0.f;
  float mKstar = 0.f;
  float mKt = 0.f;
  float mMt = 0.f;
  float mMassInv = 0.f;
  float mMult = 0.f;
  float mCent = 0.f;
  double mMass12 = 0.;
  double mMass13 = 0.;
  double mMassTot2 = 0.;

  // mc
  ROOT::Math::PtEtaPhiMVector mTrueParticle1{};
  ROOT::Math::PtEtaPhiMVector mTrueParticle2{};
  float mTrueKstar = 0.f;
  float mTrueKt = 0.f;
  float mTrueMt = 0.f;
  float mTrueMinv = 0.f;
  float mTrueMult = 0.f;
  float mTrueCent = 0.f;

  // cuts
  bool mHasMcPair = false;
  bool mHasMcCol = false;
  float mKstarMin = -1.f;
  float mKstarMax = -1.f;
  float mKtMin = -1.f;
  float mKtMax = -1.f;
  float mMtMin = -1.f;
  float mMtMax = -1.f;
  float mMassInvMin = -1;
  float mMassInvMax = -1;

  // flags
  bool mPlot1d = true;
  bool mPlot2d = true;

  bool mPlotKstarVsMtVsMult = false;
  bool mPlotKstarVsMtVsMultVsCent = false;

  bool mPlotKstarVsMtVsPt1VsPt2 = false;
  bool mPlotKstarVsMtVsPt1VsPt2VsMult = false;
  bool mPlotKstarVsMtVsPt1VsPt2VsMultVsCent = false;

  bool mPlotKstarVsMtVsMass1VsMass2 = false;
  bool mPlotKstarVsMtVsMass1VsMass2VsMult = false;
  bool mPlotKstarVsMtVsMass1VsMass2VsMultVsCent = false;

  bool mPlotKstarVsMtVsMass1VsMass2VsPt1VsPt2 = false;
  bool mPlotKstarVsMtVsMass1VsMass2VsPt1VsPt2VsMult = false;
  bool mPlotKstarVsMtVsMass1VsMass2VsPt1VsPt2VsMultVsCent = false;

  bool mPlotDalitz = false;

  // qa
  bool mPairCorrelationQa = false;
  bool mEventMixingQa = false;

  std::set<int64_t> mUniqueParticles1PerEvent = {};
  std::set<int64_t> mUniqueParticles2PerEvent = {};

  std::unordered_map<int64_t, int> mUniqueParticles1PerMixingBin = {};
  std::unordered_map<int64_t, int> mUniqueParticles2PerMixingBin = {};
};

}; // namespace pairhistmanager
}; // namespace o2::analysis::femto
#endif // PWGCF_FEMTO_CORE_PAIRHISTMANAGER_H_
