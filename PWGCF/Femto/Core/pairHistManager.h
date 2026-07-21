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

#include "PWGCF/Femto/Core/femtoSpherHarMath.h"
#include "PWGCF/Femto/Core/femtoUtils.h"
#include "PWGCF/Femto/Core/histManager.h"
#include "PWGCF/Femto/Core/modes.h"

#include "Common/Core/RecoDecay.h"

#include <CommonConstants/MathConstants.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/Logger.h>

#include <Math/Vector4D.h> // IWYU pragma: keep (do not replace with Math/Vector4Dfwd.h)
#include <Math/Vector4Dfwd.h>
#include <TH1.h>
#include <TH3.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <map>
#include <memory>
#include <string>
#include <string_view>
#include <tuple>
#include <unordered_set>
#include <vector>

namespace o2::analysis::femto::pairhistmanager
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
  // higher dimensions with pt, invariant mass and pair pt
  kKstarVsMinvVsPtPairVsMult,
  kKstarVsMtVsMinvVsPtPairVsMult,
  kKstarVsMtVsMinvVsPtPairVsPt1VsPt2,
  kKstarVsMtVsMinvVsPtPairVsPt1VsPt2VsMult,
  kKstarVsMtVsMinvVsPtPairVsPt1VsPt2VsMultVsCent,
  // dalitz plots
  kDalitz, // between a track and pos/neg daughter of another particle
  // reco-vs-mc-truth correlation (requires BOTH a reco pair and matched mc info)
  kTrueKstarVsKstar,
  kTrueKtVsKt,
  kTrueMtVsMt,
  kTrueMinvVsMinv,
  kTrueMultVsMult,
  kTrueCentVsCent,
  // pure mc-truth pair (no reco counterpart, kMc without kReco)
  kTrueKstar,
  kTrueKt,
  kTrueMt,
  kTrueMinv1D,
  kTruePt1VsPt2,
  kTruePt1VsKstar,
  kTruePt2VsKstar,
  kTruePt1VsKt,
  kTruePt2VsKt,
  kTruePt1VsMt,
  kTruePt2VsMt,
  kTrueKstarVsKt,
  kTrueKstarVsMt,
  kTrueKstarVsMult,
  kTrueKstarVsCent,
  kTrueDeltaEtaDeltaPhi,
  kTrueQout,
  kTrueQside,
  kTrueQlong,
  kTrueQoutQsideQlong,
  kTrueKstarVsMtVsMult,
  kTrueKstarVsMtVsMultVsCent,
  kTrueKstarVsMtVsPt1VsPt2,
  kTrueKstarVsMtVsPt1VsPt2VsMult,
  kTrueKstarVsMtVsPt1VsPt2VsMultVsCent,
  kTrueQoutVsQout,
  kTrueQsideVsQside,
  kTrueQlongVsQlong,

  // mixing qa
  kSeNpart1VsNpart2,                         // number of unique particles 1 vs unique number of particles 2 in each same event
  kMeMixingWindowRaw,                        // mixing window size
  kMeMixingWindowEffective,                  // mixing window size, counting event pairs with particle pairs
  kMeNpart1VsNpart2,                         // number of unique particles 1 vs number of unique particles 2 in each mixed event
  kMeVtz1VsMult1VsCent1VsVtz2VsMult2VsCent2, // correlation of event properties in each mixing bin

  // angular
  kDeltaEtaDeltaPhi,
  // Bertsch-Pratt 3D decomposition in LCMS
  kQout,
  kQside,
  kQlong,
  kQoutQsideQlong,
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
  o2::framework::Configurable<bool> usePdgMass{"usePdgMass", true, "(Reco) Use PDF masses for 4-vectors. If false, use reconstructed mass (if available). Not consulted for pure mc-truth pairs, which always use PDG mass"};
  o2::framework::Configurable<bool> plot1D{"plot1D", true, "(Reco/Mc) Enable 1D histograms"};
  o2::framework::Configurable<bool> plot2D{"plot2D", true, "(Reco/Mc) Enable 2D histograms"};
  o2::framework::Configurable<bool> plotKstarVsMtVsMult{"plotKstarVsMtVsMult", false, "(Reco/Mc) Enable 3D histogram (Kstar Vs Mt Vs Mult)"};
  o2::framework::Configurable<bool> plotKstarVsMtVsMultVsCent{"plotKstarVsMtVsMultVsCent", false, "(Reco/Mc) Enable 4D histogram (Kstar Vs Mt Vs Mult Vs Cent)"};
  o2::framework::Configurable<bool> plotKstarVsMtVsPt1VsPt2{"plotKstarVsMtVsPt1VsPt2", false, "(Reco/Mc) Enable 4D histogram (Kstar Vs Mt Vs Pt1 Vs Pt2)"};
  o2::framework::Configurable<bool> plotKstarVsMtVsPt1VsPt2VsMult{"plotKstarVsMtVsPt1VsPt2VsMult", false, "(Reco/Mc) Enable 5D histogram (Kstar Vs Mt Vs Pt1 Vs Pt2 Vs Mult)"};
  o2::framework::Configurable<bool> plotKstarVsMtVsPt1VsPt2VsMultVsCent{"plotKstarVsMtVsPt1VsPt2VsMultVsCent", false, "(Reco/Mc) Enable 6D histogram (Kstar Vs Mt Vs Pt1 Vs Pt2 Vs Mult Vs Cent)"};
  o2::framework::Configurable<bool> plotKstarVsMtVsMass1VsMass2{"plotKstarVsMtVsMass1VsMass2", false, "(Reco) Enable 4D histogram (Kstar Vs Mt Vs Mass1 Vs Mass2)"};
  o2::framework::Configurable<bool> plotKstarVsMtVsMass1VsMass2VsMult{"plotKstarVsMtVsMass1VsMass2VsMult", false, "(Reco) Enable 5D histogram (Kstar Vs Mt Vs Mass1 Vs Mass2 Vs Mult)"};
  o2::framework::Configurable<bool> plotKstarVsMtVsMass1VsMass2VsMultVsCent{"plotKstarVsMtVsMass1VsMass2VsMultVsCent", false, "(Reco) Enable 6D histogram (Kstar Vs Mt Vs Mass1 Vs Mass2 Vs Mult Vs Cent)"};
  o2::framework::Configurable<bool> plotKstarVsMtVsMass1VsMass2VsPt1VsPt2{"plotKstarVsMtVsMass1VsMass2VsPt1VsPt2", false, "(Reco) Enable 6D histogram (Kstar Vs Mt Vs Pt1 Vs Pt2 Vs Mass1 Vs Mass2)"};
  o2::framework::Configurable<bool> plotKstarVsMtVsMass1VsMass2VsPt1VsPt2VsMult{"plotKstarVsMtVsMass1VsMass2VsPt1VsPt2VsMult", false, "(Reco) Enable 7D histogram (Kstar Vs Mt Vs Pt1 Vs Pt2 Vs Mass1 Vs Mass2 Vs Mult)"};
  o2::framework::Configurable<bool> plotKstarVsMtVsMass1VsMass2VsPt1VsPt2VsMultVsCent{"plotKstarVsMtVsMass1VsMass2VsPt1VsPt2VsMultVsCent", false, "(Reco) Enable 8D histogram (Kstar Vs Mt Vs Pt1 Vs Pt2 Vs Mass1 Vs Mass2 Vs Mult Vs Cent)"};
  o2::framework::Configurable<bool> plotKstarVsMtVsMinv1VsPt1VsPt2{"plotKstarVsMtVsMinv1VsPt1VsPt2", false, "(Reco) Enable 5D histogram (Kstar Vs Mt Vs Minv Vs Pt1 Vs Pt2)"};
  o2::framework::Configurable<bool> plotKstarVsMtVsMinv1VsPt1VsPt2VsMult{"plotKstarVsMtVsMinv1VsPt1VsPt2VsMult", false, "(Reco) Enable 6D histogram (Kstar Vs Mt Vs Minv Vs Pt1 Vs Pt2 Vs Mult)"};
  o2::framework::Configurable<bool> plotKstarVsMtVsMinv1VsPt1VsPt2VsMultVsCent{"plotKstarVsMtVsMinv1VsPt1VsPt2VsMultVsCent", false, "(Reco) Enable 7D histogram (Kstar Vs Mt Vs Minv Vs Pt1 Vs Pt2 Vs Mult Vs Cent)"};
  o2::framework::Configurable<bool> plotKstarVsMinvVsPtPairVsMult{"plotKstarVsMinvVsPtPairVsMult", false, "(Reco) Enable 4D histogram (Kstar Vs Minv Vs PtPair Vs Mult)"};
  o2::framework::Configurable<bool> plotKstarVsMtVsMinvVsPtPairVsMult{"plotKstarVsMtVsMinvVsPtPairVsMult", false, "(Reco) Enable 5D histogram (Kstar Vs Mt Vs Minv Vs PtPair Vs Mult)"};
  o2::framework::Configurable<bool> plotKstarVsMtVsMinvVsPtPairVsPt1VsPt2{"plotKstarVsMtVsMinvVsPtPairVsPt1VsPt2", false, "(Reco) Enable 6D histogram (Kstar Vs Mt Vs Minv Vs PtPair Vs Pt1 Vs Pt2)"};
  o2::framework::Configurable<bool> plotKstarVsMtVsMinvVsPtPairVsPt1VsPt2VsMult{"plotKstarVsMtVsMinvVsPtPairVsPt1VsPt2VsMult", false, "(Reco) Enable 7D histogram (Kstar Vs Mt Vs Minv Vs PtPair Vs Pt1 Vs Pt2 Vs Mult)"};
  o2::framework::Configurable<bool> plotKstarVsMtVsMinvVsPtPairVsPt1VsPt2VsMultVsCent{"plotKstarVsMtVsMinvVsPtPairVsPt1VsPt2VsMultVsCent", false, "(Reco) Enable 8D histogram (Kstar Vs Mt Vs Minv Vs PtPair Vs Pt1 Vs Pt2 Vs Mult Vs Cent)"};
  o2::framework::Configurable<bool> plotDalitz{"plotDalitz", false, "(Reco) Enable dalitz plot. Not supported for pure mc-truth pairs (no trackTable/daughter structure)"};
  o2::framework::Configurable<bool> plotDeltaEtaDeltaPhi{"plotDeltaEtaDeltaPhi", false, "(Reco/Mc) Plot #Delta#phi vs #Delta#eta"};
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
  o2::framework::ConfigurableAxis ptPair{"ptPair", {{120, 0, 12}}, "Pair transverse momentum binning (from summed px and py of both particles)"};
  o2::framework::ConfigurableAxis dalitzMtot{"dalitzMtot", {{100, 0, 10}}, "Total invariant mass squared binning in darlitz plot"};
  o2::framework::ConfigurableAxis dalitzM12{"dalitzM12", {{100, 0, 10}}, "Mass12 binning of darlitz plot"};
  o2::framework::ConfigurableAxis dalitzM13{"dalitzM13", {{100, 0, 10}}, "Mass13 binning of darlitz plot"};
  o2::framework::Configurable<int> transverseMassType{"transverseMassType", static_cast<int>(modes::TransverseMassType::kAveragePdgMass), "Type of transverse mass (0-> Average Pdg Mass, 1-> Reduced Pdg Mass, 2-> Mt from combined 4 vector)"};
  o2::framework::ConfigurableAxis binningDeltaEta{"binningDeltaEta", {{35, -1.6, 1.6}}, "Delta eta"};
  o2::framework::ConfigurableAxis binningDeltaPhi{"binningDeltaPhi", {{35, -o2::constants::math::PIHalf, 3 * o2::constants::math::PIHalf}}, "Delta phi"};
  o2::framework::Configurable<bool> plotBertschPratt{"plotBertschPratt", false, "(Reco/Mc) Enable 1D projections and 3D (q_out, q_side, q_long) Bertsch-Pratt histograms in LCMS"};
  o2::framework::ConfigurableAxis qout{"qout", {{300, -1.5f, 1.5f}}, "q_{out} (GeV/c) in LCMS"};
  o2::framework::ConfigurableAxis qside{"qside", {{300, -1.5f, 1.5f}}, "q_{side} (GeV/c) in LCMS"};
  o2::framework::ConfigurableAxis qlong{"qlong", {{300, -1.5f, 1.5f}}, "q_{long} (GeV/c) in LCMS"};
  o2::framework::Configurable<bool> plotSH{"plotSH", false, "(Reco) Enable spherical-harmonics decomposition of the pair momentum-difference vector"};
  o2::framework::Configurable<int> shLMax{"shLMax", 2, "Maximum l for SH decomposition (0..5). FemtoUniverse hard-codes 1."};
  o2::framework::Configurable<int> shFrame{"shFrame", 1, "SH reference frame/variable: 0=LCMS non-identical (k*), 1=LCMS identical (qinv, FemtoUniverse default), 2=PRF (q_PRF, matches FemtoUniverse isIdenPRF=true)"};
  o2::framework::ConfigurableAxis shKstar{"shKstar", {{60, 0.0f, 0.3f}}, "k*/qinv binning for SH histograms"};
  o2::framework::Configurable<bool> shUseCent{"shUseCent", false, "SH: bin by centrality instead of multiplicity"};
  o2::framework::ConfigurableAxis shCentBins{"shCentBins", {o2::framework::VARIABLE_WIDTH, 0.0f, 200.0f}, "SH: multiplicity/centrality bin edges (like FemtoUniverse confMultKstarBins)"};
  o2::framework::ConfigurableAxis shKtBins{"shKtBins", {o2::framework::VARIABLE_WIDTH, 0.1f, 0.2f, 0.3f, 0.4f}, "SH: kT bin edges (like FemtoUniverse confKtKstarBins)"};
  o2::framework::Configurable<bool> shPlot1D{"shPlot1D", false, "(SH) Also fill 1D qinv/k* numerator/denominator (h1D) per (mult,kT) bin"};
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
      {kPt1VsKstar, o2::framework::HistType::kTH2F, "hPt1VsKstar", "p_{T,1} vs k*; p_{T,1} (GeV/#it{c}); k* (GeV/#it{c})"},
      {kPt2VsKstar, o2::framework::HistType::kTH2F, "hPt2VsKstar", "p_{T,2} vs k*; p_{T,2} (GeV/#it{c}); k* (GeV/#it{c})"},
      {kPt1VsKt, o2::framework::HistType::kTH2F, "hPt1VsKt", "p_{T,1} vs k_{T}; p_{T,1} (GeV/#it{c}); k_{T} (GeV/#it{c})"},
      {kPt2VsKt, o2::framework::HistType::kTH2F, "hPt2VsKt", "p_{T,2} vs k_{T}; p_{T,2} (GeV/#it{c}); k_{T} (GeV/#it{c})"},
      {kPt1VsMt, o2::framework::HistType::kTH2F, "hPt1VsMt", "p_{T,1} vs m_{T}; p_{T,1} (GeV/#it{c}); m_{T} (GeV/#it{c}^{2})"},
      {kPt2VsMt, o2::framework::HistType::kTH2F, "hPt2VsMt", "p_{T,2} vs m_{T}; p_{T,2} (GeV/#it{c}); m_{T} (GeV/#it{c}^{2})"},
      {kKstarVsKt, o2::framework::HistType::kTH2F, "hKstarVsKt", "k* vs k_{T}; k* (GeV/#it{c}); k_{T} (GeV/#it{c})"},
      {kKstarVsMt, o2::framework::HistType::kTH2F, "hKstarVsMt", "k* vs m_{T}; k* (GeV/#it{c}); m_{T} (GeV/#it{c}^{2})"},
      {kKstarVsMult, o2::framework::HistType::kTH2F, "hKstarVsMult", "k* vs Multiplicity; k* (GeV/#it{c}); Multiplicity"},
      {kKstarVsCent, o2::framework::HistType::kTH2F, "hKstarVsCent", "k* vs Centrality (Mult. Percentile); k* (GeV/#it{c}); Centrality (%)"},
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
      {kKstarVsMtVsMass1VsMass2, o2::framework::HistType::kTHnSparseF, "hKstarVsMtVsMass1VsMass2", "k* vs m_{T} vs m_{1} vs m_{2}; k* (GeV/#it{c}); m_{T} (GeV/#it{c}^{2}); m_{1} (GeV/#it{c}^{2}); m_{2} (GeV/#it{c}^{2});"},
      {kKstarVsMtVsMass1VsMass2VsMult, o2::framework::HistType::kTHnSparseF, "hKstarVsMtVsMass1VsMass2VsMult", "k* vs m_{T} vs m_{1} vs m_{2} vs multiplicity; k* (GeV/#it{c}); m_{T} (GeV/#it{c}^{2}); m_{1} (GeV/#it{c}^{2}); m_{2} (GeV/#it{c}^{2}); Multiplicity;"},
      {kKstarVsMtVsMass1VsMass2VsMultVsCent, o2::framework::HistType::kTHnSparseF, "hKstarVsMtVsMass1VsMass2VsMultVsCent", "k* vs m_{T} vs m_{1} vs m_{2} vs multiplicity vs centrality; k* (GeV/#it{c}); m_{T} (GeV/#it{c}^{2}); m_{1} (GeV/#it{c}^{2}); m_{2} (GeV/#it{c}^{2}); Multiplicity; Centrality (%);"},
      // n-D with pt and mass
      {kKstarVsMtVsMass1VsMass2VsPt1VsPt2, o2::framework::HistType::kTHnSparseF, "hKstarVsMtVsMass1VsMass2VsPt1VsPt2", "k* vs m_{T} vs m_{1} vs m_{2} vs p_{T,1} vs p_{T,2}; k* (GeV/#it{c}); m_{T} (GeV/#it{c}^{2}); m_{1} (GeV/#it{c}^{2}); m_{2} (GeV/#it{c}^{2}); p_{T,1} (GeV/#it{c}); p_{T,2} (GeV/#it{c});"},
      {kKstarVsMtVsMass1VsMass2VsPt1VsPt2VsMult, o2::framework::HistType::kTHnSparseF, "hKstarVsMtVsMass1VsMass2VsPt1VsPt2VsMult", "k* vs m_{T} vs m_{1} vs m_{2} vs p_{T,1} vs p_{T,2} vs multiplicity; k* (GeV/#it{c}); m_{T} (GeV/#it{c}^{2}); m_{1} (GeV/#it{c}^{2}); m_{2} (GeV/#it{c}^{2}); p_{T,1} (GeV/#it{c}); p_{T,2} (GeV/#it{c}); Multiplicity;"},
      {kKstarVsMtVsMass1VsMass2VsPt1VsPt2VsMultVsCent, o2::framework::HistType::kTHnSparseF, "hKstarVsMtVsMass1VsMass2VsPt1VsPt2VsMultVsCent", "k* vs m_{T} vs m_{1} vs m_{2} vs p_{T,1} vs p_{T,2} vs multiplicity vs centrality; k* (GeV/#it{c}); m_{T} (GeV/#it{c}^{2}); m_{1} (GeV/#it{c}^{2}); m_{2} (GeV/#it{c}^{2}); p_{T,1} (GeV/#it{c}); p_{T,2} (GeV/#it{c}); Multiplicity; Centrality (%);"},
      {kKstarVsMtVsMinvVsPt1VsPt2, o2::framework::HistType::kTHnSparseF, "hKstarVsMtVsMinvVsPt1VsPt2", "k* vs m_{T} vs m_{Inv} vs p_{T,1} vs p_{T,2}; k* (GeV/#it{c}); m_{T} (GeV/#it{c}^{2}); m_{Inv} (GeV/#it{c}^{2}); p_{T,1} (GeV/#it{c}); p_{T,2} (GeV/#it{c})"},
      {kKstarVsMtVsMinvVsPt1VsPt2VsMult, o2::framework::HistType::kTHnSparseF, "hKstarVsMtVsMinvVsPt1VsPt2VsMult", "k* vs m_{T} vs m_{Inv} vs p_{T,1} vs p_{T,2} vs multiplicity; k* (GeV/#it{c}); m_{T} (GeV/#it{c}^{2}); m_{Inv} (GeV/#it{c}^{2}); p_{T,1} (GeV/#it{c}); p_{T,2} (GeV/#it{c}); Multiplicity"},
      {kKstarVsMtVsMinvVsPt1VsPt2VsMultVsCent, o2::framework::HistType::kTHnSparseF, "hKstarVsMtVsMinvVsPt1VsPt2VsMultVsCent", "k* vs m_{T} vs m_{Inv} vs p_{T,1} vs p_{T,2} vs multiplicity vs centrality; k* (GeV/#it{c}); m_{T} (GeV/#it{c}^{2}); m_{Inv} (GeV/#it{c}^{2}); p_{T,1} (GeV/#it{c}); p_{T,2} (GeV/#it{c}); Multiplicity; Centrality (%)"},
      {kKstarVsMinvVsPtPairVsMult, o2::framework::HistType::kTHnSparseF, "hKstarVsMinvVsPtPairVsMult", "k* vs m_{Inv} vs p_{T,pair} vs multiplicity; k* (GeV/#it{c}); m_{Inv} (GeV/#it{c}^{2}); p_{T,pair} (GeV/#it{c}); Multiplicity"},
      {kKstarVsMtVsMinvVsPtPairVsMult, o2::framework::HistType::kTHnSparseF, "hKstarVsMtVsMinvVsPtPairVsMult", "k* vs m_{T} vs m_{Inv} vs p_{T,pair} vs multiplicity; k* (GeV/#it{c}); m_{T} (GeV/#it{c}^{2}); m_{Inv} (GeV/#it{c}^{2}); p_{T,pair} (GeV/#it{c}); Multiplicity"},
      {kKstarVsMtVsMinvVsPtPairVsPt1VsPt2, o2::framework::HistType::kTHnSparseF, "hKstarVsMtVsMinvVsPtPairVsPt1VsPt2", "k* vs m_{T} vs m_{Inv} vs p_{T,pair} vs p_{T,1} vs p_{T,2}; k* (GeV/#it{c}); m_{T} (GeV/#it{c}^{2}); m_{Inv} (GeV/#it{c}^{2}); p_{T,pair} (GeV/#it{c}); p_{T,1} (GeV/#it{c}); p_{T,2} (GeV/#it{c})"},
      {kKstarVsMtVsMinvVsPtPairVsPt1VsPt2VsMult, o2::framework::HistType::kTHnSparseF, "hKstarVsMtVsMinvVsPtPairVsPt1VsPt2VsMult", "k* vs m_{T} vs m_{Inv} vs p_{T,pair} vs p_{T,1} vs p_{T,2} vs multiplicity; k* (GeV/#it{c}); m_{T} (GeV/#it{c}^{2}); m_{Inv} (GeV/#it{c}^{2}); p_{T,pair} (GeV/#it{c}); p_{T,1} (GeV/#it{c}); p_{T,2} (GeV/#it{c}); Multiplicity"},
      {kKstarVsMtVsMinvVsPtPairVsPt1VsPt2VsMultVsCent, o2::framework::HistType::kTHnSparseF, "hKstarVsMtVsMinvVsPtPairVsPt1VsPt2VsMultVsCent", "k* vs m_{T} vs m_{Inv} vs p_{T,pair} vs p_{T,1} vs p_{T,2} vs multiplicity vs centrality; k* (GeV/#it{c}); m_{T} (GeV/#it{c}^{2}); m_{Inv} (GeV/#it{c}^{2}); p_{T,pair} (GeV/#it{c}); p_{T,1} (GeV/#it{c}); p_{T,2} (GeV/#it{c}); Multiplicity; Centrality (%)"},
      {kDalitz, o2::framework::HistType::kTHnSparseF, "hDalitz", "Dalitz plot; k* (GeV/#it{c}); m^{2}_{123} (GeV/#it{c}^{2})^{2}; m^{2}_{12} (GeV/#it{c}^{2})^{2}; m^{2}_{13} (GeV/#it{c}^{2})^{2};"},
      // reco-vs-mc-truth correlation
      {kTrueKstarVsKstar, o2::framework::HistType::kTH2F, "hTrueKstarVsKstar", "k*_{True} vs k*; k*_{True} (GeV/#it{c});  k* (GeV/#it{c})"},
      {kTrueKtVsKt, o2::framework::HistType::kTH2F, "hTrueKtVsKt", "k_{T,True} vs k_{T}; k_{T,True} (GeV/#it{c});  k_{T} (GeV/#it{c})"},
      {kTrueMtVsMt, o2::framework::HistType::kTH2F, "hTrueMtVsMt", "m_{T,True} vs m_{T}; m_{T,True} (GeV/#it{c}^{2}); m_{T} (GeV/#it{c}^{2})"},
      {kTrueMinvVsMinv, o2::framework::HistType::kTH2F, "hTrueMinvVsMinv", "m_{Inv,True} vs m_{Inv}; m_{Inv,True} (GeV/#it{c}^{2}); m_{Inv} (GeV/#it{c}^{2})"},
      {kTrueMultVsMult, o2::framework::HistType::kTH2F, "hTrueMultVsMult", "Multiplicity_{True} vs Multiplicity; Multiplicity_{True} ;  Multiplicity"},
      {kTrueCentVsCent, o2::framework::HistType::kTH2F, "hTrueCentVsCent", "Centrality_{True} vs Centrality; Centrality_{True} (%); Centrality (%)"},
      // pure mc-truth pair (no reco counterpart)
      {kTrueKstar, o2::framework::HistType::kTH1F, "hTrueKstar", "k* (mc-truth pair); k*_{True} (GeV/#it{c}); Entries"},
      {kTrueKt, o2::framework::HistType::kTH1F, "hTrueKt", "transverse momentum (mc-truth pair); k_{T,True} (GeV/#it{c}); Entries"},
      {kTrueMt, o2::framework::HistType::kTH1F, "hTrueMt", "transverse mass (mc-truth pair); m_{T,True} (GeV/#it{c}^{2}); Entries"},
      {kTrueMinv1D, o2::framework::HistType::kTH1F, "hTrueMinv", "invariant mass (mc-truth pair); m_{Inv,True} (GeV/#it{c}^{2}); Entries"},
      {kTruePt1VsPt2, o2::framework::HistType::kTH2F, "hTruePt1VsPt2", "p_{T,1} vs p_{T,2} (mc-truth pair); p_{T,1,True} (GeV/#it{c}); p_{T,2,True} (GeV/#it{c})"},
      {kTruePt1VsKstar, o2::framework::HistType::kTH2F, "hTruePt1VsKstar", "p_{T,1} vs k* (mc-truth pair); p_{T,1,True} (GeV/#it{c}); k*_{True} (GeV/#it{c})"},
      {kTruePt2VsKstar, o2::framework::HistType::kTH2F, "hTruePt2VsKstar", "p_{T,2} vs k* (mc-truth pair); p_{T,2,True} (GeV/#it{c}); k*_{True} (GeV/#it{c})"},
      {kTruePt1VsKt, o2::framework::HistType::kTH2F, "hTruePt1VsKt", "p_{T,1} vs k_{T} (mc-truth pair); p_{T,1,True} (GeV/#it{c}); k_{T,True} (GeV/#it{c})"},
      {kTruePt2VsKt, o2::framework::HistType::kTH2F, "hTruePt2VsKt", "p_{T,2} vs k_{T} (mc-truth pair); p_{T,2,True} (GeV/#it{c}); k_{T,True} (GeV/#it{c})"},
      {kTruePt1VsMt, o2::framework::HistType::kTH2F, "hTruePt1VsMt", "p_{T,1} vs m_{T} (mc-truth pair); p_{T,1,True} (GeV/#it{c}); m_{T,True} (GeV/#it{c}^{2})"},
      {kTruePt2VsMt, o2::framework::HistType::kTH2F, "hTruePt2VsMt", "p_{T,2} vs m_{T} (mc-truth pair); p_{T,2,True} (GeV/#it{c}); m_{T,True} (GeV/#it{c}^{2})"},
      {kTrueKstarVsKt, o2::framework::HistType::kTH2F, "hTrueKstarVsKt", "k* vs k_{T} (mc-truth pair); k*_{True} (GeV/#it{c}); k_{T,True} (GeV/#it{c})"},
      {kTrueKstarVsMt, o2::framework::HistType::kTH2F, "hTrueKstarVsMt", "k* vs m_{T} (mc-truth pair); k*_{True} (GeV/#it{c}); m_{T,True} (GeV/#it{c}^{2})"},
      {kTrueKstarVsMult, o2::framework::HistType::kTH2F, "hTrueKstarVsMult", "k* vs Multiplicity (mc-truth pair); k*_{True} (GeV/#it{c}); Multiplicity_{True}"},
      {kTrueKstarVsCent, o2::framework::HistType::kTH2F, "hTrueKstarVsCent", "k* vs Centrality (mc-truth pair); k*_{True} (GeV/#it{c}); Centrality_{True} (%)"},
      {kTrueDeltaEtaDeltaPhi, o2::framework::HistType::kTH2F, "hTrueDeltaEtaDeltaPhi", "#Delta#phi vs #Delta#eta (mc-truth pair); #Delta#phi_{True}; #Delta#eta_{True}"},
      {kTrueQout, o2::framework::HistType::kTH1F, "hTrueQout", "q_{out} in LCMS (mc-truth pair); q_{out,True} (GeV/#it{c}); Entries"},
      {kTrueQside, o2::framework::HistType::kTH1F, "hTrueQside", "q_{side} in LCMS (mc-truth pair); q_{side,True} (GeV/#it{c}); Entries"},
      {kTrueQlong, o2::framework::HistType::kTH1F, "hTrueQlong", "q_{long} in LCMS (mc-truth pair); q_{long,True} (GeV/#it{c}); Entries"},
      {kTrueQoutQsideQlong, o2::framework::HistType::kTH3F, "hTrueQoutQsideQlong", "Bertsch-Pratt 3D (mc-truth pair); q_{out,True} (GeV/#it{c}); q_{side,True} (GeV/#it{c}); q_{long,True} (GeV/#it{c})"},
      {kTrueKstarVsMtVsMult, o2::framework::HistType::kTHnSparseF, "hTrueKstarVsMtVsMult", "k* vs m_{T} vs multiplicity (mc-truth pair); k*_{True} (GeV/#it{c}); m_{T,True} (GeV/#it{c}^{2}); Multiplicity_{True};"},
      {kTrueKstarVsMtVsMultVsCent, o2::framework::HistType::kTHnSparseF, "hTrueKstarVsMtVsMultVsCent", "k* vs m_{T} vs multiplicity vs centrality (mc-truth pair); k*_{True} (GeV/#it{c}); m_{T,True} (GeV/#it{c}^{2}); Multiplicity_{True}; Centrality_{True} (%);"},
      {kTrueKstarVsMtVsPt1VsPt2, o2::framework::HistType::kTHnSparseF, "hTrueKstarVsMtVsPt1VsPt2", "k* vs m_{T} vs p_{T,1} vs p_{T,2} (mc-truth pair); k*_{True} (GeV/#it{c}); m_{T,True} (GeV/#it{c}^{2}); p_{T,1,True} (GeV/#it{c}); p_{T,2,True} (GeV/#it{c});"},
      {kTrueKstarVsMtVsPt1VsPt2VsMult, o2::framework::HistType::kTHnSparseF, "hTrueKstarVsMtVsPt1VsPt2VsMult", "k* vs m_{T} vs p_{T,1} vs p_{T,2} vs multiplicity (mc-truth pair); k*_{True} (GeV/#it{c}); m_{T,True} (GeV/#it{c}^{2}); p_{T,1,True} (GeV/#it{c}); p_{T,2,True} (GeV/#it{c}); Multiplicity_{True};"},
      {kTrueKstarVsMtVsPt1VsPt2VsMultVsCent, o2::framework::HistType::kTHnSparseF, "hTrueKstarVsMtVsPt1VsPt2VsMultVsCent", "k* vs m_{T} vs p_{T,1} vs p_{T,2} vs multiplicity vs centrality (mc-truth pair); k*_{True} (GeV/#it{c}); m_{T,True} (GeV/#it{c}^{2}); p_{T,1,True} (GeV/#it{c}); p_{T,2,True} (GeV/#it{c}); Multiplicity_{True}; Centrality_{True} (%);"},
      // mixing qa
      {kSeNpart1VsNpart2, o2::framework::HistType::kTH2F, "hSeNpart1VsNpart2", "# unique particle 1 vs # unique particle 2 in each same event; # partilce 1; # particle 2"},
      {kMeMixingWindowRaw, o2::framework::HistType::kTH1F, "hMeMixingWindowRaw", "Raw Mixing Window; Raw Mixing Window  Entries"},
      {kMeMixingWindowEffective, o2::framework::HistType::kTH1F, "hMeMixingWindowEffective", "Effective Mixing Window; Effective Mixing Window; Entries"},
      {kMeNpart1VsNpart2, o2::framework::HistType::kTH2F, "hMeNpart1VsNpart2", "# unique particle 1 vs # unique partilce 2 in each mixing bin; # partilce 1; # particle 2"},
      {kMeVtz1VsMult1VsCent1VsVtz2VsMult2VsCent2, o2::framework::HistType::kTHnSparseF, "hVtz1VsMult1VsCent1VsVtz2VsMult2VsCent2", "Mixing bins; V_{z,1} (cm); multiplicity_{1}; centrality_{1} (%); V_{z,2} (cm); multiplicity_{2}; centrality_{2} (%)"},
      // angular
      {kDeltaEtaDeltaPhi, o2::framework::HistType::kTH2F, "hDeltaEtaDeltaPhi", "#Delta#phi vs #Delta#eta; #Delta#phi; #Delta#eta"},
      // Bertsch-Pratt 3D decomposition in LCMS
      {kQout, o2::framework::HistType::kTH1F, "hQout", "q_{out} in LCMS; q_{out} (GeV/#it{c}); Entries"},
      {kQside, o2::framework::HistType::kTH1F, "hQside", "q_{side} in LCMS; q_{side} (GeV/#it{c}); Entries"},
      {kQlong, o2::framework::HistType::kTH1F, "hQlong", "q_{long} in LCMS; q_{long} (GeV/#it{c}); Entries"},
      {kQoutQsideQlong, o2::framework::HistType::kTH3F, "hQoutQsideQlong", "Bertsch-Pratt 3D; q_{out} (GeV/#it{c}); q_{side} (GeV/#it{c}); q_{long} (GeV/#it{c})"},
      {kTrueQoutVsQout, o2::framework::HistType::kTH2F, "hTrueQoutVsQout", "q_{out,True} vs q_{out}; q_{out,True} (GeV/#it{c}); q_{out} (GeV/#it{c})"},
      {kTrueQsideVsQside, o2::framework::HistType::kTH2F, "hTrueQsideVsQside", "q_{side,True} vs q_{side}; q_{side,True} (GeV/#it{c}); q_{side} (GeV/#it{c})"},
      {kTrueQlongVsQlong, o2::framework::HistType::kTH2F, "hTrueQlongVsQlong", "q_{long,True} vs q_{long}; q_{long,True} (GeV/#it{c}); q_{long} (GeV/#it{c})"},
    }};

// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define PAIR_HIST_ANALYSIS_MAP(confAnalysis)                                                                                                                                                                                                    \
  {kKstar, {(confAnalysis).kstar}},                                                                                                                                                                                                             \
    {kKt, {(confAnalysis).kt}},                                                                                                                                                                                                                 \
    {kMt, {(confAnalysis).mt}},                                                                                                                                                                                                                 \
    {kMinv, {(confAnalysis).massInv}},                                                                                                                                                                                                          \
    {kPt1VsPt2, {(confAnalysis).pt1, (confAnalysis).pt2}},                                                                                                                                                                                      \
    {kPt1VsKstar, {(confAnalysis).pt1, (confAnalysis).kstar}},                                                                                                                                                                                  \
    {kPt2VsKstar, {(confAnalysis).pt2, (confAnalysis).kstar}},                                                                                                                                                                                  \
    {kPt1VsKt, {(confAnalysis).pt1, (confAnalysis).kt}},                                                                                                                                                                                        \
    {kPt2VsKt, {(confAnalysis).pt2, (confAnalysis).kt}},                                                                                                                                                                                        \
    {kPt1VsMt, {(confAnalysis).pt1, (confAnalysis).mt}},                                                                                                                                                                                        \
    {kPt2VsMt, {(confAnalysis).pt2, (confAnalysis).mt}},                                                                                                                                                                                        \
    {kKstarVsKt, {(confAnalysis).kstar, (confAnalysis).kt}},                                                                                                                                                                                    \
    {kKstarVsMt, {(confAnalysis).kstar, (confAnalysis).mt}},                                                                                                                                                                                    \
    {kKstarVsMult, {(confAnalysis).kstar, (confAnalysis).multiplicity}},                                                                                                                                                                        \
    {kKstarVsCent, {(confAnalysis).kstar, (confAnalysis).centrality}},                                                                                                                                                                          \
    {kKstarVsMass1, {(confAnalysis).kstar, (confAnalysis).mass1}},                                                                                                                                                                              \
    {kKstarVsMass2, {(confAnalysis).kstar, (confAnalysis).mass2}},                                                                                                                                                                              \
    {kMass1VsMass2, {(confAnalysis).mass1, (confAnalysis).mass2}},                                                                                                                                                                              \
    {kKstarVsMinv, {(confAnalysis).kstar, (confAnalysis).massInv}},                                                                                                                                                                             \
    {kPt1VsMinv, {(confAnalysis).pt1, (confAnalysis).massInv}},                                                                                                                                                                                 \
    {kPt2VsMinv, {(confAnalysis).pt2, (confAnalysis).massInv}},                                                                                                                                                                                 \
    {kKstarVsMtVsMult, {(confAnalysis).kstar, (confAnalysis).mt, (confAnalysis).multiplicity}},                                                                                                                                                 \
    {kKstarVsMtVsMultVsCent, {(confAnalysis).kstar, (confAnalysis).mt, (confAnalysis).multiplicity, (confAnalysis).centrality}},                                                                                                                \
    {kKstarVsMtVsPt1VsPt2, {(confAnalysis).kstar, (confAnalysis).mt, (confAnalysis).pt1, (confAnalysis).pt2}},                                                                                                                                  \
    {kKstarVsMtVsPt1VsPt2VsMult, {(confAnalysis).kstar, (confAnalysis).mt, (confAnalysis).pt1, (confAnalysis).pt2, (confAnalysis).multiplicity}},                                                                                               \
    {kKstarVsMtVsPt1VsPt2VsMultVsCent, {(confAnalysis).kstar, (confAnalysis).mt, (confAnalysis).pt1, (confAnalysis).pt2, (confAnalysis).multiplicity, (confAnalysis).centrality}},                                                              \
    {kKstarVsMtVsMass1VsMass2, {(confAnalysis).kstar, (confAnalysis).mt, (confAnalysis).mass1, (confAnalysis).mass2}},                                                                                                                          \
    {kKstarVsMtVsMass1VsMass2VsMult, {(confAnalysis).kstar, (confAnalysis).mt, (confAnalysis).mass1, (confAnalysis).mass2, (confAnalysis).multiplicity}},                                                                                       \
    {kKstarVsMtVsMass1VsMass2VsMultVsCent, {(confAnalysis).kstar, (confAnalysis).mt, (confAnalysis).mass1, (confAnalysis).mass2, (confAnalysis).multiplicity, (confAnalysis).centrality}},                                                      \
    {kKstarVsMtVsMass1VsMass2VsPt1VsPt2, {(confAnalysis).kstar, (confAnalysis).mt, (confAnalysis).mass1, (confAnalysis).mass2, (confAnalysis).pt1, (confAnalysis).pt2}},                                                                        \
    {kKstarVsMtVsMass1VsMass2VsPt1VsPt2VsMult, {(confAnalysis).kstar, (confAnalysis).mt, (confAnalysis).mass1, (confAnalysis).mass2, (confAnalysis).pt1, (confAnalysis).pt2, (confAnalysis).multiplicity}},                                     \
    {kKstarVsMtVsMass1VsMass2VsPt1VsPt2VsMultVsCent, {(confAnalysis).kstar, (confAnalysis).mt, (confAnalysis).mass1, (confAnalysis).mass2, (confAnalysis).pt1, (confAnalysis).pt2, (confAnalysis).multiplicity, (confAnalysis).centrality}},    \
    {kKstarVsMtVsMinvVsPt1VsPt2, {(confAnalysis).kstar, (confAnalysis).mt, (confAnalysis).massInv, (confAnalysis).pt1, (confAnalysis).pt2}},                                                                                                    \
    {kKstarVsMtVsMinvVsPt1VsPt2VsMult, {(confAnalysis).kstar, (confAnalysis).mt, (confAnalysis).massInv, (confAnalysis).pt1, (confAnalysis).pt2, (confAnalysis).multiplicity}},                                                                 \
    {kKstarVsMtVsMinvVsPt1VsPt2VsMultVsCent, {(confAnalysis).kstar, (confAnalysis).mt, (confAnalysis).massInv, (confAnalysis).pt1, (confAnalysis).pt2, (confAnalysis).multiplicity, (confAnalysis).centrality}},                                \
    {kKstarVsMinvVsPtPairVsMult, {(confAnalysis).kstar, (confAnalysis).massInv, (confAnalysis).ptPair, (confAnalysis).multiplicity}},                                                                                                           \
    {kKstarVsMtVsMinvVsPtPairVsMult, {(confAnalysis).kstar, (confAnalysis).mt, (confAnalysis).massInv, (confAnalysis).ptPair, (confAnalysis).multiplicity}},                                                                                    \
    {kKstarVsMtVsMinvVsPtPairVsPt1VsPt2, {(confAnalysis).kstar, (confAnalysis).mt, (confAnalysis).massInv, (confAnalysis).ptPair, (confAnalysis).pt1, (confAnalysis).pt2}},                                                                     \
    {kKstarVsMtVsMinvVsPtPairVsPt1VsPt2VsMult, {(confAnalysis).kstar, (confAnalysis).mt, (confAnalysis).massInv, (confAnalysis).ptPair, (confAnalysis).pt1, (confAnalysis).pt2, (confAnalysis).multiplicity}},                                  \
    {kKstarVsMtVsMinvVsPtPairVsPt1VsPt2VsMultVsCent, {(confAnalysis).kstar, (confAnalysis).mt, (confAnalysis).massInv, (confAnalysis).ptPair, (confAnalysis).pt1, (confAnalysis).pt2, (confAnalysis).multiplicity, (confAnalysis).centrality}}, \
    {kDalitz, {(confAnalysis).kstar, (confAnalysis).dalitzMtot, (confAnalysis).dalitzM12, (confAnalysis).dalitzM13}},                                                                                                                           \
    {kDeltaEtaDeltaPhi, {(confAnalysis).binningDeltaPhi, (confAnalysis).binningDeltaEta}},                                                                                                                                                      \
    {kQout, {(confAnalysis).qout}},                                                                                                                                                                                                             \
    {kQside, {(confAnalysis).qside}},                                                                                                                                                                                                           \
    {kQlong, {(confAnalysis).qlong}},                                                                                                                                                                                                           \
    {kQoutQsideQlong, {(confAnalysis).qout, (confAnalysis).qside, (confAnalysis).qlong}},

// mixing-qa entries are independent of reco vs mc-truth status — both the reco
// analysis path and the pure mc-truth path need them whenever kSe/kMe is set
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define PAIR_HIST_MIXING_QA_MAP(confMixing)                                            \
  {kSeNpart1VsNpart2, {(confMixing).particleBinning, (confMixing).particleBinning}},   \
    {kMeMixingWindowRaw, {(confMixing).particleBinning}},                              \
    {kMeMixingWindowEffective, {(confMixing).particleBinning}},                        \
    {kMeNpart1VsNpart2, {(confMixing).particleBinning, (confMixing).particleBinning}}, \
    {kMeVtz1VsMult1VsCent1VsVtz2VsMult2VsCent2, {(confMixing).vtxBins, (confMixing).multBins, (confMixing).centBins, (confMixing).vtxBins, (confMixing).multBins, (confMixing).centBins}},

// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define PAIR_HIST_MC_MAP(conf)                                     \
  {kTrueKstarVsKstar, {(conf).kstar, (conf).kstar}},               \
    {kTrueKtVsKt, {(conf).kt, (conf).kt}},                         \
    {kTrueMtVsMt, {(conf).mt, (conf).mt}},                         \
    {kTrueMinvVsMinv, {(conf).massInv, (conf).massInv}},           \
    {kTrueMultVsMult, {(conf).multiplicity, (conf).multiplicity}}, \
    {kTrueCentVsCent, {(conf).centrality, (conf).centrality}},     \
    {kTrueQoutVsQout, {(conf).qout, (conf).qout}},                 \
    {kTrueQsideVsQside, {(conf).qside, (conf).qside}},             \
    {kTrueQlongVsQlong, {(conf).qlong, (conf).qlong}},

// pure mc-truth pair (no reco counterpart) — reuses the same analysis binning,
// since there is no separate "true" axis (conf)iguration: the truth value IS the
// analysis-level value for this path.
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define PAIR_HIST_MC_TRUTH_MAP(conf)                                                                          \
  {kTrueKstar, {(conf).kstar}},                                                                               \
    {kTrueKt, {(conf).kt}},                                                                                   \
    {kTrueMt, {(conf).mt}},                                                                                   \
    {kTrueMinv1D, {(conf).massInv}},                                                                          \
    {kTruePt1VsPt2, {(conf).pt1, (conf).pt2}},                                                                \
    {kTruePt1VsKstar, {(conf).pt1, (conf).kstar}},                                                            \
    {kTruePt2VsKstar, {(conf).pt2, (conf).kstar}},                                                            \
    {kTruePt1VsKt, {(conf).pt1, (conf).kt}},                                                                  \
    {kTruePt2VsKt, {(conf).pt2, (conf).kt}},                                                                  \
    {kTruePt1VsMt, {(conf).pt1, (conf).mt}},                                                                  \
    {kTruePt2VsMt, {(conf).pt2, (conf).mt}},                                                                  \
    {kTrueKstarVsKt, {(conf).kstar, (conf).kt}},                                                              \
    {kTrueKstarVsMt, {(conf).kstar, (conf).mt}},                                                              \
    {kTrueKstarVsMult, {(conf).kstar, (conf).multiplicity}},                                                  \
    {kTrueKstarVsCent, {(conf).kstar, (conf).centrality}},                                                    \
    {kTrueDeltaEtaDeltaPhi, {(conf).binningDeltaPhi, (conf).binningDeltaEta}},                                \
    {kTrueQout, {(conf).qout}},                                                                               \
    {kTrueQside, {(conf).qside}},                                                                             \
    {kTrueQlong, {(conf).qlong}},                                                                             \
    {kTrueQoutQsideQlong, {(conf).qout, (conf).qside, (conf).qlong}},                                         \
    {kTrueKstarVsMtVsMult, {(conf).kstar, (conf).mt, (conf).multiplicity}},                                   \
    {kTrueKstarVsMtVsMultVsCent, {(conf).kstar, (conf).mt, (conf).multiplicity, (conf).centrality}},          \
    {kTrueKstarVsMtVsPt1VsPt2, {(conf).kstar, (conf).mt, (conf).pt1, (conf).pt2}},                            \
    {kTrueKstarVsMtVsPt1VsPt2VsMult, {(conf).kstar, (conf).mt, (conf).pt1, (conf).pt2, (conf).multiplicity}}, \
    {kTrueKstarVsMtVsPt1VsPt2VsMultVsCent, {(conf).kstar, (conf).mt, (conf).pt1, (conf).pt2, (conf).multiplicity, (conf).centrality}},

template <typename T1, typename T2>
auto makePairHistSpecMap(T1 const& confPairBinning, T2 const& confMixing)
{
  return std::map<PairHist, std::vector<o2::framework::AxisSpec>>{
    PAIR_HIST_ANALYSIS_MAP(confPairBinning)
      PAIR_HIST_MIXING_QA_MAP(confMixing)};
};

template <typename T1, typename T2>
auto makePairMcHistSpecMap(T1 const& confPairBinning, T2 const& confMixing)
{
  return std::map<PairHist, std::vector<o2::framework::AxisSpec>>{
    PAIR_HIST_ANALYSIS_MAP(confPairBinning)
      PAIR_HIST_MIXING_QA_MAP(confMixing)
        PAIR_HIST_MC_MAP(confPairBinning)};
};

// for the pure mc-truth-only pair builder (kMc without kReco): needs both the
// mc-truth binning AND the mixing-qa entries (kSe/kMe histograms are filled
// for this path too, via initSeMixingQa/initMeMixingQa being gated on kSe/kMe
// alone, independent of kReco/kMc).
template <typename T1, typename T2>
auto makePairMcTruthHistSpecMap(T1 const& confPairBinning, T2 const& confMixing)
{
  return std::map<PairHist, std::vector<o2::framework::AxisSpec>>{
    PAIR_HIST_MC_TRUTH_MAP(confPairBinning)
      PAIR_HIST_MIXING_QA_MAP(confMixing)};
};

#undef PAIR_HIST_ANALYSIS_MAP
#undef PAIR_HIST_MIXING_QA_MAP
#undef PAIR_HIST_MC_MAP
#undef PAIR_HIST_MC_TRUTH_MAP

constexpr char PrefixTrackTrackSe[] = "TrackTrack/SE/";
constexpr char PrefixTrackTrackMe[] = "TrackTrack/ME/";

constexpr char PrefixTrackV0Se[] = "TrackV0/SE/";
constexpr char PrefixTrackV0Me[] = "TrackV0/ME/";

constexpr char PrefixV0V0Se[] = "V0V0/SE/";
constexpr char PrefixV0V0Me[] = "V0V0/ME/";

constexpr char PrefixTrackResonanceSe[] = "TrackResonance/SE/";
constexpr char PrefixTrackResonanceMe[] = "TrackResonance/ME/";

constexpr char PrefixV0ResonanceSe[] = "V0Resonance/SE/";
constexpr char PrefixV0ResonanceMe[] = "V0Resonance/ME/";

constexpr char PrefixTrackCascadeSe[] = "TrackCascade/SE/";
constexpr char PrefixTrackCascadeMe[] = "TrackCascade/ME/";

constexpr char PrefixTrackKinkSe[] = "TrackKink/SE/";
constexpr char PrefixTrackKinkMe[] = "TrackKink/ME/";

constexpr char PrefixMcParticleMcParticleSe[] = "McParticleMcParticle/SE/";
constexpr char PrefixMcParticleMcParticleMe[] = "McParticleMcParticle/ME/";

constexpr std::string_view AnalysisDir = "Analysis/";
constexpr std::string_view MixingQaDir = "MixingQA/";
constexpr std::string_view McDir = "MC/";

template <auto& prefix,
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

    mPlotKstarVsMtVsMinvVsPt1VsPt2 = ConfPairBinning.plotKstarVsMtVsMinv1VsPt1VsPt2.value;
    mPlotKstarVsMtVsMinvVsPt1VsPt2VsMult = ConfPairBinning.plotKstarVsMtVsMinv1VsPt1VsPt2VsMult.value;
    mPlotKstarVsMtVsMinvVsPt1VsPt2VsMultVsCent = ConfPairBinning.plotKstarVsMtVsMinv1VsPt1VsPt2VsMultVsCent.value;

    mPlotKstarVsMinvVsPtPairVsMult = ConfPairBinning.plotKstarVsMinvVsPtPairVsMult.value;
    mPlotKstarVsMtVsMinvVsPtPairVsMult = ConfPairBinning.plotKstarVsMtVsMinvVsPtPairVsMult.value;
    mPlotKstarVsMtVsMinvVsPtPairVsPt1VsPt2 = ConfPairBinning.plotKstarVsMtVsMinvVsPtPairVsPt1VsPt2.value;
    mPlotKstarVsMtVsMinvVsPtPairVsPt1VsPt2VsMult = ConfPairBinning.plotKstarVsMtVsMinvVsPtPairVsPt1VsPt2VsMult.value;
    mPlotKstarVsMtVsMinvVsPtPairVsPt1VsPt2VsMultVsCent = ConfPairBinning.plotKstarVsMtVsMinvVsPtPairVsPt1VsPt2VsMultVsCent.value;

    mPlotDalitz = ConfPairBinning.plotDalitz.value;
    mPlotDeltaEtaDeltaPhi = ConfPairBinning.plotDeltaEtaDeltaPhi.value;
    mPlotBertschPratt = ConfPairBinning.plotBertschPratt.value;

    mPlotSH = ConfPairBinning.plotSH.value;
    mShUseCent = ConfPairBinning.shUseCent.value;
    mShLMax = ConfPairBinning.shLMax.value;
    mShFrame = ConfPairBinning.shFrame.value;
    mShPlot1D = ConfPairBinning.shPlot1D.value;
    if (mPlotSH) {
      mShKstarSpec = {ConfPairBinning.shKstar, "k* (GeV/#it{c})"};
      mYlm.initializeYlms();
      // copy bin edges, stripping the leading VARIABLE_WIDTH (0) marker
      mShCentEdges.assign(ConfPairBinning.shCentBins.value.begin() + 1, ConfPairBinning.shCentBins.value.end());
      mShKtEdges.assign(ConfPairBinning.shKtBins.value.begin() + 1, ConfPairBinning.shKtBins.value.end());
    }

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

    if constexpr (isFlagSet(mode, modes::Mode::kReco)) {
      initAnalysis(Specs);
    }

    // reco-vs-truth correlation: requires BOTH a reco pair and matched mc info
    if constexpr (isFlagSet(mode, modes::Mode::kReco) && isFlagSet(mode, modes::Mode::kMc)) {
      initMc(Specs);
    }

    // pure mc-truth pair: requires mc info WITHOUT a reco pair
    if constexpr (isFlagSet(mode, modes::Mode::kMc) && !isFlagSet(mode, modes::Mode::kReco)) {
      initMcTruth(Specs);
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

    // set pair pT, since kT is defined as half the transverse momentum of the summed 4-vector, this is just twice kT
    mPtPair = 2.f * mKt;

    // set mT
    mMt = getMt(mParticle1, mParticle2);

    // set Minv
    mMassInv = getMinv(mParticle1, mParticle2);

    // set kstar
    mKstar = getKstar(mParticle1, mParticle2);

    if (mPlotBertschPratt) {
      std::tie(mQout, mQside, mQlong) = computeBertschPrattLCMS(mParticle1, mParticle2);
    }

    if (mPlotSH) {
      std::tie(mShKv, mShOut, mShSide, mShLong) = computeShKinematics(mParticle1, mParticle2);
    }

    if (mPlotDeltaEtaDeltaPhi) {
      mDeltaEta = particle1.eta() - particle2.eta();
      mDeltaPhi = RecoDecay::constrainAngle(particle1.phi() - particle2.phi(), -o2::constants::math::PIHalf);
    }

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

    mTrueParticle1 = ROOT::Math::PtEtaPhiMVector(mcParticle1.pt(), mcParticle1.eta(), mcParticle1.phi(), mPdgMass1);
    mTrueParticle2 = ROOT::Math::PtEtaPhiMVector(mcParticle2.pt(), mcParticle2.eta(), mcParticle2.phi(), mPdgMass2);

    // compute true kinematics
    mTrueKt = getKt(mTrueParticle1, mTrueParticle2);
    mTrueMt = getMt(mTrueParticle1, mTrueParticle2);
    mTrueMinv = getMinv(mTrueParticle1, mTrueParticle2);
    mTrueKstar = getKstar(mTrueParticle1, mTrueParticle2);

    if (mPlotBertschPratt) {
      std::tie(mTrueQout, mTrueQside, mTrueQlong) =
        computeBertschPrattLCMS(mTrueParticle1, mTrueParticle2);
    }
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

  // pure mc-truth pair: particles ARE the truth, no reco track/trackTable exists.
  // Mass is always the PDG mass here — there is no "reconstructed mass" concept
  // for a truth particle, so mUsePdgMass is not consulted in this path.
  // NOTE: Dalitz plots are not supported here — there is no trackTable/daughter
  // structure for a pure mc-truth particle-particle pair.
  template <typename T1, typename T2>
  void setPairMcTruth(T1 const& particle1, T2 const& particle2)
  {
    mTrueParticle1 = ROOT::Math::PtEtaPhiMVector(particle1.pt(), particle1.eta(), particle1.phi(), mPdgMass1);
    mTrueParticle2 = ROOT::Math::PtEtaPhiMVector(particle2.pt(), particle2.eta(), particle2.phi(), mPdgMass2);

    mTrueKt = getKt(mTrueParticle1, mTrueParticle2);
    mTrueMt = getMt(mTrueParticle1, mTrueParticle2);
    mTrueMinv = getMinv(mTrueParticle1, mTrueParticle2);
    mTrueKstar = getKstar(mTrueParticle1, mTrueParticle2);

    if (mPlotBertschPratt) {
      std::tie(mTrueQout, mTrueQside, mTrueQlong) = computeBertschPrattLCMS(mTrueParticle1, mTrueParticle2);
    }
    if (mPlotDeltaEtaDeltaPhi) {
      mTrueDeltaEta = particle1.eta() - particle2.eta();
      mTrueDeltaPhi = RecoDecay::constrainAngle(particle1.phi() - particle2.phi(), -o2::constants::math::PIHalf);
    }
  }

  // same-event: single (truth) collision for true mult/cent
  template <typename T1, typename T2, typename T3>
  void setPairMcTruth(T1 const& particle1, T2 const& particle2, T3 const& col)
  {
    setPairMcTruth(particle1, particle2);
    mTrueMult = col.mult();
    mTrueCent = col.cent();
  }

  // mixed-event: two (truth) collisions, averaged mult/cent — same convention as setPair
  template <typename T1, typename T2, typename T3, typename T4>
  void setPairMcTruth(T1 const& particle1, T2 const& particle2, T3 const& col1, T4 const& col2)
  {
    setPairMcTruth(particle1, particle2);
    mTrueMult = 0.5f * (col1.mult() + col2.mult());
    mTrueCent = 0.5f * (col1.cent() + col2.cent());
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
    if constexpr (isFlagSet(mode, modes::Mode::kReco)) {
      fillAnalysis();
    }
    if constexpr (isFlagSet(mode, modes::Mode::kReco) && isFlagSet(mode, modes::Mode::kMc)) {
      fillMc();
    }
    if constexpr (isFlagSet(mode, modes::Mode::kMc) && !isFlagSet(mode, modes::Mode::kReco)) {
      fillMcTruth();
    }
  }

  template <typename T1, typename T2>
  void trackParticlesPerEvent(T1 const& particle1, T2 const& particle2)
  {
    if (!mPairCorrelationQa) {
      return;
    }
    mParticles1PerEvent.insert(particle1.globalIndex());
    mParticles2PerEvent.insert(particle2.globalIndex());
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
    mParticles1PerEvent.clear();
    mParticles1PerEvent.reserve(100);
    mParticles2PerEvent.clear();
    mParticles2PerEvent.reserve(100);
  }

  void fillMixingQaSe()
  {
    if (mPairCorrelationQa) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(MixingQaDir) + HIST(getHistName(kSeNpart1VsNpart2, HistTable)), mParticles1PerEvent.size(), mParticles2PerEvent.size());
    }
  }

  void fillMixingQaMePerEvent()
  {
    if (mPairCorrelationQa) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(MixingQaDir) + HIST(getHistName(kMeNpart1VsNpart2, HistTable)), mParticles1PerEvent.size(), mParticles2PerEvent.size());
    }
  }

  void fillMixingQaMePerMixingBin(int windowSizeRaw, int windowSizeEffective)
  {
    if (mPairCorrelationQa) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(MixingQaDir) + HIST(getHistName(kMeMixingWindowRaw, HistTable)), windowSizeRaw);
      mHistogramRegistry->fill(HIST(prefix) + HIST(MixingQaDir) + HIST(getHistName(kMeMixingWindowEffective, HistTable)), windowSizeEffective);
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

    if (mPlotKstarVsMtVsMinvVsPt1VsPt2) {
      mHistogramRegistry->add(analysisDir + getHistNameV2(kKstarVsMtVsMinvVsPt1VsPt2, HistTable), getHistDesc(kKstarVsMtVsMinvVsPt1VsPt2, HistTable), getHistType(kKstarVsMtVsMinvVsPt1VsPt2, HistTable), {Specs.at(kKstarVsMtVsMinvVsPt1VsPt2)});
    }
    if (mPlotKstarVsMtVsMinvVsPt1VsPt2VsMult) {
      mHistogramRegistry->add(analysisDir + getHistNameV2(kKstarVsMtVsMinvVsPt1VsPt2VsMult, HistTable), getHistDesc(kKstarVsMtVsMinvVsPt1VsPt2VsMult, HistTable), getHistType(kKstarVsMtVsMinvVsPt1VsPt2VsMult, HistTable), {Specs.at(kKstarVsMtVsMinvVsPt1VsPt2VsMult)});
    }
    if (mPlotKstarVsMtVsMinvVsPt1VsPt2VsMultVsCent) {
      mHistogramRegistry->add(analysisDir + getHistNameV2(kKstarVsMtVsMinvVsPt1VsPt2VsMultVsCent, HistTable), getHistDesc(kKstarVsMtVsMinvVsPt1VsPt2VsMultVsCent, HistTable), getHistType(kKstarVsMtVsMinvVsPt1VsPt2VsMultVsCent, HistTable), {Specs.at(kKstarVsMtVsMinvVsPt1VsPt2VsMultVsCent)});
    }

    if (mPlotKstarVsMinvVsPtPairVsMult) {
      mHistogramRegistry->add(analysisDir + getHistNameV2(kKstarVsMinvVsPtPairVsMult, HistTable), getHistDesc(kKstarVsMinvVsPtPairVsMult, HistTable), getHistType(kKstarVsMinvVsPtPairVsMult, HistTable), {Specs.at(kKstarVsMinvVsPtPairVsMult)});
    }
    if (mPlotKstarVsMtVsMinvVsPtPairVsMult) {
      mHistogramRegistry->add(analysisDir + getHistNameV2(kKstarVsMtVsMinvVsPtPairVsMult, HistTable), getHistDesc(kKstarVsMtVsMinvVsPtPairVsMult, HistTable), getHistType(kKstarVsMtVsMinvVsPtPairVsMult, HistTable), {Specs.at(kKstarVsMtVsMinvVsPtPairVsMult)});
    }
    if (mPlotKstarVsMtVsMinvVsPtPairVsPt1VsPt2) {
      mHistogramRegistry->add(analysisDir + getHistNameV2(kKstarVsMtVsMinvVsPtPairVsPt1VsPt2, HistTable), getHistDesc(kKstarVsMtVsMinvVsPtPairVsPt1VsPt2, HistTable), getHistType(kKstarVsMtVsMinvVsPtPairVsPt1VsPt2, HistTable), {Specs.at(kKstarVsMtVsMinvVsPtPairVsPt1VsPt2)});
    }
    if (mPlotKstarVsMtVsMinvVsPtPairVsPt1VsPt2VsMult) {
      mHistogramRegistry->add(analysisDir + getHistNameV2(kKstarVsMtVsMinvVsPtPairVsPt1VsPt2VsMult, HistTable), getHistDesc(kKstarVsMtVsMinvVsPtPairVsPt1VsPt2VsMult, HistTable), getHistType(kKstarVsMtVsMinvVsPtPairVsPt1VsPt2VsMult, HistTable), {Specs.at(kKstarVsMtVsMinvVsPtPairVsPt1VsPt2VsMult)});
    }
    if (mPlotKstarVsMtVsMinvVsPtPairVsPt1VsPt2VsMultVsCent) {
      mHistogramRegistry->add(analysisDir + getHistNameV2(kKstarVsMtVsMinvVsPtPairVsPt1VsPt2VsMultVsCent, HistTable), getHistDesc(kKstarVsMtVsMinvVsPtPairVsPt1VsPt2VsMultVsCent, HistTable), getHistType(kKstarVsMtVsMinvVsPtPairVsPt1VsPt2VsMultVsCent, HistTable), {Specs.at(kKstarVsMtVsMinvVsPtPairVsPt1VsPt2VsMultVsCent)});
    }

    if (mPlotDalitz) {
      mHistogramRegistry->add(analysisDir + getHistNameV2(kDalitz, HistTable), getHistDesc(kDalitz, HistTable), getHistType(kDalitz, HistTable), {Specs.at(kDalitz)});
    }
    if (mPlotDeltaEtaDeltaPhi) {
      mHistogramRegistry->add(analysisDir + getHistNameV2(kDeltaEtaDeltaPhi, HistTable), getHistDesc(kDeltaEtaDeltaPhi, HistTable), getHistType(kDeltaEtaDeltaPhi, HistTable), {Specs.at(kDeltaEtaDeltaPhi)});
    }
    if (mPlotBertschPratt) {
      mHistogramRegistry->add(analysisDir + getHistNameV2(kQout, HistTable), getHistDesc(kQout, HistTable), getHistType(kQout, HistTable), {Specs.at(kQout)});
      mHistogramRegistry->add(analysisDir + getHistNameV2(kQside, HistTable), getHistDesc(kQside, HistTable), getHistType(kQside, HistTable), {Specs.at(kQside)});
      mHistogramRegistry->add(analysisDir + getHistNameV2(kQlong, HistTable), getHistDesc(kQlong, HistTable), getHistType(kQlong, HistTable), {Specs.at(kQlong)});
      mHistogramRegistry->add(analysisDir + getHistNameV2(kQoutQsideQlong, HistTable), getHistDesc(kQoutQsideQlong, HistTable), getHistType(kQoutQsideQlong, HistTable), {Specs.at(kQoutQsideQlong)});
    }
    if (mPlotSH) {
      const int nJM = (mShLMax + 1) * (mShLMax + 1);
      const int nCent = static_cast<int>(mShCentEdges.size()) - 1; // n edges -> n - 1 bins
      const int nKt = static_cast<int>(mShKtEdges.size()) - 1;
      mShYlmBuffer.assign(nJM, {});

      mShReal.resize(nCent);
      mShImag.resize(nCent);
      mShCov.resize(nCent);
      mSh1D.resize(nCent);
      mShBinCount.resize(nCent);
      for (int iCent = 0; iCent < nCent; ++iCent) {
        mShReal[iCent].resize(nKt);
        mShImag[iCent].resize(nKt);
        mShCov[iCent].resize(nKt);
        mSh1D[iCent].resize(nKt);
        mShBinCount[iCent].resize(nKt);
        // folder name: mult_{low}_{high}
        const std::string centFolder = "mult_" + std::to_string(static_cast<int>(mShCentEdges[iCent])) +
                                       "_" + std::to_string(static_cast<int>(mShCentEdges[iCent + 1]));
        for (int iKt = 0; iKt < nKt; ++iKt) {
          mShReal[iCent][iKt].resize(nJM);
          mShImag[iCent][iKt].resize(nJM);
          // folder name: kT_{low*100}_{high*100}
          std::string ktFolder = "kT_";
          ktFolder += std::to_string(static_cast<int>(mShKtEdges[iKt] * 100.0));
          ktFolder += "_";
          ktFolder += std::to_string(static_cast<int>(mShKtEdges[iKt + 1] * 100.0));
          std::string dir = std::string(prefix) + std::string(AnalysisDir) + "SH/";
          dir += centFolder;
          dir += "/";
          dir += ktFolder;
          dir += "/";

          int ihist = 0;
          for (int l = 0; l <= mShLMax; ++l) {
            for (int m = -l; m <= l; ++m) {
              std::string lm = std::to_string(l);
              lm += (m < 0) ? std::to_string(l - m) : std::to_string(m);
              std::string nameRe = dir;
              nameRe += "ReYlm";
              nameRe += lm;
              std::string nameIm = dir;
              nameIm += "ImYlm";
              nameIm += lm;
              // shared "Y_{l}^{m}" suffix for both titles
              std::string ylmLabel = "Y_{";
              ylmLabel += std::to_string(l);
              ylmLabel += "}^{";
              ylmLabel += std::to_string(m);
              ylmLabel += "}";
              std::string titleRe = "Re ";
              titleRe += ylmLabel;
              titleRe += "; k* (GeV/#it{c}); Re[A_{l}^{m}]";
              std::string titleIm = "Im ";
              titleIm += ylmLabel;
              titleIm += "; k* (GeV/#it{c}); Im[A_{l}^{m}]";
              mShReal[iCent][iKt][ihist] = mHistogramRegistry->add<TH1>(nameRe.c_str(), titleRe.c_str(), o2::framework::kTH1D, {mShKstarSpec});
              mShImag[iCent][iKt][ihist] = mHistogramRegistry->add<TH1>(nameIm.c_str(), titleIm.c_str(), o2::framework::kTH1D, {mShKstarSpec});
              mShReal[iCent][iKt][ihist]->Sumw2();
              mShImag[iCent][iKt][ihist]->Sumw2();
              ++ihist;
            }
          }

          // SH covariance TH3D
          const int nAxisLM = 2 * nJM;
          const o2::framework::AxisSpec covLmAxis{nAxisLM, -0.5, static_cast<double>(nAxisLM) - 0.5, "l,m #times (re,im)"};
          std::string nameCov = dir;
          nameCov += "Cov";
          mShCov[iCent][iKt] = mHistogramRegistry->add<TH3>(nameCov.c_str(), "SH covariance; k* (GeV/#it{c}); l,m; l,m", o2::framework::kTH3D, {mShKstarSpec, covLmAxis, covLmAxis});
          mShCov[iCent][iKt]->Sumw2();

          std::string nameBinCount = dir;
          nameBinCount += "BinCount";
          mShBinCount[iCent][iKt] = mHistogramRegistry->add<TH1>(nameBinCount.c_str(), "SH bin occupancy; k* (GeV/#it{c}); Entries", o2::framework::kTH1D, {mShKstarSpec});

          if (mShPlot1D) {
            std::string name1D = dir;
            name1D += "h1D";
            mSh1D[iCent][iKt] = mHistogramRegistry->add<TH1>(name1D.c_str(), "1D distribution; k* (GeV/#it{c}); Entries", o2::framework::kTH1D, {mShKstarSpec});
            mSh1D[iCent][iKt]->Sumw2();
          }
        }
      }
    }
  }

  // reco-vs-truth correlation histograms (kReco and kMc both set)
  void initMc(std::map<PairHist, std::vector<o2::framework::AxisSpec>> const& Specs)
  {
    std::string mcDir = std::string(prefix) + std::string(McDir);
    mHistogramRegistry->add(mcDir + getHistNameV2(kTrueKstarVsKstar, HistTable), getHistDesc(kTrueKstarVsKstar, HistTable), getHistType(kTrueKstarVsKstar, HistTable), {Specs.at(kTrueKstarVsKstar)});
    mHistogramRegistry->add(mcDir + getHistNameV2(kTrueKtVsKt, HistTable), getHistDesc(kTrueKtVsKt, HistTable), getHistType(kTrueKtVsKt, HistTable), {Specs.at(kTrueKtVsKt)});
    mHistogramRegistry->add(mcDir + getHistNameV2(kTrueMtVsMt, HistTable), getHistDesc(kTrueMtVsMt, HistTable), getHistType(kTrueMtVsMt, HistTable), {Specs.at(kTrueMtVsMt)});
    mHistogramRegistry->add(mcDir + getHistNameV2(kTrueMinvVsMinv, HistTable), getHistDesc(kTrueMinvVsMinv, HistTable), getHistType(kTrueMinvVsMinv, HistTable), {Specs.at(kTrueMinvVsMinv)});
    mHistogramRegistry->add(mcDir + getHistNameV2(kTrueMultVsMult, HistTable), getHistDesc(kTrueMultVsMult, HistTable), getHistType(kTrueMultVsMult, HistTable), {Specs.at(kTrueMultVsMult)});
    mHistogramRegistry->add(mcDir + getHistNameV2(kTrueCentVsCent, HistTable), getHistDesc(kTrueCentVsCent, HistTable), getHistType(kTrueCentVsCent, HistTable), {Specs.at(kTrueCentVsCent)});

    if (mPlotBertschPratt) {
      mHistogramRegistry->add(mcDir + getHistNameV2(kTrueQoutVsQout, HistTable), getHistDesc(kTrueQoutVsQout, HistTable), getHistType(kTrueQoutVsQout, HistTable), {Specs.at(kTrueQoutVsQout)});
      mHistogramRegistry->add(mcDir + getHistNameV2(kTrueQsideVsQside, HistTable), getHistDesc(kTrueQsideVsQside, HistTable), getHistType(kTrueQsideVsQside, HistTable), {Specs.at(kTrueQsideVsQside)});
      mHistogramRegistry->add(mcDir + getHistNameV2(kTrueQlongVsQlong, HistTable), getHistDesc(kTrueQlongVsQlong, HistTable), getHistType(kTrueQlongVsQlong, HistTable), {Specs.at(kTrueQlongVsQlong)});
    }
  }

  // pure mc-truth pair (no reco counterpart, kMc without kReco) — reuses the
  // same mPlot1d/mPlot2d/mPlotDeltaEtaDeltaPhi/mPlotBertschPratt flags as the
  // reco analysis path, since these are histogram-content flags, not mode flags
  void initMcTruth(std::map<PairHist, std::vector<o2::framework::AxisSpec>> const& Specs)
  {
    std::string mcDir = std::string(prefix) + std::string(McDir);
    if (mPlot1d) {
      mHistogramRegistry->add(mcDir + getHistNameV2(kTrueKstar, HistTable), getHistDesc(kTrueKstar, HistTable), getHistType(kTrueKstar, HistTable), {Specs.at(kTrueKstar)});
      mHistogramRegistry->add(mcDir + getHistNameV2(kTrueKt, HistTable), getHistDesc(kTrueKt, HistTable), getHistType(kTrueKt, HistTable), {Specs.at(kTrueKt)});
      mHistogramRegistry->add(mcDir + getHistNameV2(kTrueMt, HistTable), getHistDesc(kTrueMt, HistTable), getHistType(kTrueMt, HistTable), {Specs.at(kTrueMt)});
      mHistogramRegistry->add(mcDir + getHistNameV2(kTrueMinv1D, HistTable), getHistDesc(kTrueMinv1D, HistTable), getHistType(kTrueMinv1D, HistTable), {Specs.at(kTrueMinv1D)});
    }
    if (mPlot2d) {
      mHistogramRegistry->add(mcDir + getHistNameV2(kTruePt1VsPt2, HistTable), getHistDesc(kTruePt1VsPt2, HistTable), getHistType(kTruePt1VsPt2, HistTable), {Specs.at(kTruePt1VsPt2)});
      mHistogramRegistry->add(mcDir + getHistNameV2(kTruePt1VsKstar, HistTable), getHistDesc(kTruePt1VsKstar, HistTable), getHistType(kTruePt1VsKstar, HistTable), {Specs.at(kTruePt1VsKstar)});
      mHistogramRegistry->add(mcDir + getHistNameV2(kTruePt2VsKstar, HistTable), getHistDesc(kTruePt2VsKstar, HistTable), getHistType(kTruePt2VsKstar, HistTable), {Specs.at(kTruePt2VsKstar)});
      mHistogramRegistry->add(mcDir + getHistNameV2(kTruePt1VsKt, HistTable), getHistDesc(kTruePt1VsKt, HistTable), getHistType(kTruePt1VsKt, HistTable), {Specs.at(kTruePt1VsKt)});
      mHistogramRegistry->add(mcDir + getHistNameV2(kTruePt2VsKt, HistTable), getHistDesc(kTruePt2VsKt, HistTable), getHistType(kTruePt2VsKt, HistTable), {Specs.at(kTruePt2VsKt)});
      mHistogramRegistry->add(mcDir + getHistNameV2(kTruePt1VsMt, HistTable), getHistDesc(kTruePt1VsMt, HistTable), getHistType(kTruePt1VsMt, HistTable), {Specs.at(kTruePt1VsMt)});
      mHistogramRegistry->add(mcDir + getHistNameV2(kTruePt2VsMt, HistTable), getHistDesc(kTruePt2VsMt, HistTable), getHistType(kTruePt2VsMt, HistTable), {Specs.at(kTruePt2VsMt)});
      mHistogramRegistry->add(mcDir + getHistNameV2(kTrueKstarVsKt, HistTable), getHistDesc(kTrueKstarVsKt, HistTable), getHistType(kTrueKstarVsKt, HistTable), {Specs.at(kTrueKstarVsKt)});
      mHistogramRegistry->add(mcDir + getHistNameV2(kTrueKstarVsMt, HistTable), getHistDesc(kTrueKstarVsMt, HistTable), getHistType(kTrueKstarVsMt, HistTable), {Specs.at(kTrueKstarVsMt)});
      mHistogramRegistry->add(mcDir + getHistNameV2(kTrueKstarVsMult, HistTable), getHistDesc(kTrueKstarVsMult, HistTable), getHistType(kTrueKstarVsMult, HistTable), {Specs.at(kTrueKstarVsMult)});
      mHistogramRegistry->add(mcDir + getHistNameV2(kTrueKstarVsCent, HistTable), getHistDesc(kTrueKstarVsCent, HistTable), getHistType(kTrueKstarVsCent, HistTable), {Specs.at(kTrueKstarVsCent)});
    }
    if (mPlotDeltaEtaDeltaPhi) {
      mHistogramRegistry->add(mcDir + getHistNameV2(kTrueDeltaEtaDeltaPhi, HistTable), getHistDesc(kTrueDeltaEtaDeltaPhi, HistTable), getHistType(kTrueDeltaEtaDeltaPhi, HistTable), {Specs.at(kTrueDeltaEtaDeltaPhi)});
    }
    if (mPlotBertschPratt) {
      mHistogramRegistry->add(mcDir + getHistNameV2(kTrueQout, HistTable), getHistDesc(kTrueQout, HistTable), getHistType(kTrueQout, HistTable), {Specs.at(kTrueQout)});
      mHistogramRegistry->add(mcDir + getHistNameV2(kTrueQside, HistTable), getHistDesc(kTrueQside, HistTable), getHistType(kTrueQside, HistTable), {Specs.at(kTrueQside)});
      mHistogramRegistry->add(mcDir + getHistNameV2(kTrueQlong, HistTable), getHistDesc(kTrueQlong, HistTable), getHistType(kTrueQlong, HistTable), {Specs.at(kTrueQlong)});
      mHistogramRegistry->add(mcDir + getHistNameV2(kTrueQoutQsideQlong, HistTable), getHistDesc(kTrueQoutQsideQlong, HistTable), getHistType(kTrueQoutQsideQlong, HistTable), {Specs.at(kTrueQoutQsideQlong)});
    }
    if (mPlotKstarVsMtVsMult) {
      mHistogramRegistry->add(mcDir + getHistNameV2(kTrueKstarVsMtVsMult, HistTable), getHistDesc(kTrueKstarVsMtVsMult, HistTable), getHistType(kTrueKstarVsMtVsMult, HistTable), {Specs.at(kTrueKstarVsMtVsMult)});
    }
    if (mPlotKstarVsMtVsMultVsCent) {
      mHistogramRegistry->add(mcDir + getHistNameV2(kTrueKstarVsMtVsMultVsCent, HistTable), getHistDesc(kTrueKstarVsMtVsMultVsCent, HistTable), getHistType(kTrueKstarVsMtVsMultVsCent, HistTable), {Specs.at(kTrueKstarVsMtVsMultVsCent)});
    }
    if (mPlotKstarVsMtVsPt1VsPt2) {
      mHistogramRegistry->add(mcDir + getHistNameV2(kTrueKstarVsMtVsPt1VsPt2, HistTable), getHistDesc(kTrueKstarVsMtVsPt1VsPt2, HistTable), getHistType(kTrueKstarVsMtVsPt1VsPt2, HistTable), {Specs.at(kTrueKstarVsMtVsPt1VsPt2)});
    }
    if (mPlotKstarVsMtVsPt1VsPt2VsMult) {
      mHistogramRegistry->add(mcDir + getHistNameV2(kTrueKstarVsMtVsPt1VsPt2VsMult, HistTable), getHistDesc(kTrueKstarVsMtVsPt1VsPt2VsMult, HistTable), getHistType(kTrueKstarVsMtVsPt1VsPt2VsMult, HistTable), {Specs.at(kTrueKstarVsMtVsPt1VsPt2VsMult)});
    }
    if (mPlotKstarVsMtVsPt1VsPt2VsMultVsCent) {
      mHistogramRegistry->add(mcDir + getHistNameV2(kTrueKstarVsMtVsPt1VsPt2VsMultVsCent, HistTable), getHistDesc(kTrueKstarVsMtVsPt1VsPt2VsMultVsCent, HistTable), getHistType(kTrueKstarVsMtVsPt1VsPt2VsMultVsCent, HistTable), {Specs.at(kTrueKstarVsMtVsPt1VsPt2VsMultVsCent)});
    }
  }

  void initSeMixingQa(std::map<PairHist, std::vector<o2::framework::AxisSpec>> const& Specs)
  {
    std::string dir = std::string(prefix) + std::string(MixingQaDir);
    if (mPairCorrelationQa) {
      mHistogramRegistry->add(dir + getHistNameV2(kSeNpart1VsNpart2, HistTable), getHistDesc(kSeNpart1VsNpart2, HistTable), getHistType(kSeNpart1VsNpart2, HistTable), {Specs.at(kSeNpart1VsNpart2)});
    }
  }

  void initMeMixingQa(std::map<PairHist, std::vector<o2::framework::AxisSpec>> const& Specs)
  {
    std::string dir = std::string(prefix) + std::string(MixingQaDir);
    if (mPairCorrelationQa) {
      mHistogramRegistry->add(dir + getHistNameV2(kMeMixingWindowRaw, HistTable), getHistDesc(kMeMixingWindowRaw, HistTable), getHistType(kMeMixingWindowRaw, HistTable), {Specs.at(kMeMixingWindowRaw)});
      mHistogramRegistry->add(dir + getHistNameV2(kMeMixingWindowEffective, HistTable), getHistDesc(kMeMixingWindowEffective, HistTable), getHistType(kMeMixingWindowEffective, HistTable), {Specs.at(kMeMixingWindowEffective)});
      mHistogramRegistry->add(dir + getHistNameV2(kMeNpart1VsNpart2, HistTable), getHistDesc(kMeNpart1VsNpart2, HistTable), getHistType(kMeNpart1VsNpart2, HistTable), {Specs.at(kMeNpart1VsNpart2)});
    }
    if (mEventMixingQa) {
      mHistogramRegistry->add(dir + getHistNameV2(kMeVtz1VsMult1VsCent1VsVtz2VsMult2VsCent2, HistTable), getHistDesc(kMeVtz1VsMult1VsCent1VsVtz2VsMult2VsCent2, HistTable), getHistType(kMeVtz1VsMult1VsCent1VsVtz2VsMult2VsCent2, HistTable), {Specs.at(kMeVtz1VsMult1VsCent1VsVtz2VsMult2VsCent2)});
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
    if (mPlotKstarVsMtVsMinvVsPt1VsPt2) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kKstarVsMtVsMinvVsPt1VsPt2, HistTable)), mKstar, mMt, mMassInv, mParticle1.Pt(), mParticle2.Pt());
    }
    if (mPlotKstarVsMtVsMinvVsPt1VsPt2VsMult) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kKstarVsMtVsMinvVsPt1VsPt2VsMult, HistTable)), mKstar, mMt, mMassInv, mParticle1.Pt(), mParticle2.Pt(), mMult);
    }
    if (mPlotKstarVsMtVsMinvVsPt1VsPt2VsMultVsCent) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kKstarVsMtVsMinvVsPt1VsPt2VsMultVsCent, HistTable)), mKstar, mMt, mMassInv, mParticle1.Pt(), mParticle2.Pt(), mMult, mCent);
    }
    if (mPlotKstarVsMinvVsPtPairVsMult) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kKstarVsMinvVsPtPairVsMult, HistTable)), mKstar, mMassInv, mPtPair, mMult);
    }
    if (mPlotKstarVsMtVsMinvVsPtPairVsMult) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kKstarVsMtVsMinvVsPtPairVsMult, HistTable)), mKstar, mMt, mMassInv, mPtPair, mMult);
    }
    if (mPlotKstarVsMtVsMinvVsPtPairVsPt1VsPt2) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kKstarVsMtVsMinvVsPtPairVsPt1VsPt2, HistTable)), mKstar, mMt, mMassInv, mPtPair, mParticle1.Pt(), mParticle2.Pt());
    }
    if (mPlotKstarVsMtVsMinvVsPtPairVsPt1VsPt2VsMult) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kKstarVsMtVsMinvVsPtPairVsPt1VsPt2VsMult, HistTable)), mKstar, mMt, mMassInv, mPtPair, mParticle1.Pt(), mParticle2.Pt(), mMult);
    }
    if (mPlotKstarVsMtVsMinvVsPtPairVsPt1VsPt2VsMultVsCent) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kKstarVsMtVsMinvVsPtPairVsPt1VsPt2VsMultVsCent, HistTable)), mKstar, mMt, mMassInv, mPtPair, mParticle1.Pt(), mParticle2.Pt(), mMult, mCent);
    }
    if (mPlotDalitz) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kDalitz, HistTable)), mKstar, mMassTot2, mMass12, mMass13);
    }
    if (mPlotDeltaEtaDeltaPhi) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kDeltaEtaDeltaPhi, HistTable)), mDeltaPhi, mDeltaEta);
    }
    if (mPlotBertschPratt) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kQout, HistTable)), mQout);
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kQside, HistTable)), mQside);
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kQlong, HistTable)), mQlong);
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kQoutQsideQlong, HistTable)), mQout, mQside, mQlong);
    }
    if (mPlotSH) {
      const int iCent = findShBin(mShUseCent ? mCent : mMult, mShCentEdges);
      const int iKt = findShBin(mKt, mShKtEdges);
      if (iCent >= 0 && iKt >= 0) {
        mYlm.doYlmUpToL(mShLMax, mShOut, mShSide, mShLong, mShYlmBuffer.data());
        for (std::size_t i = 0; i < mShYlmBuffer.size(); ++i) {
          mShReal[iCent][iKt][i]->Fill(mShKv, std::real(mShYlmBuffer[i]));
          mShImag[iCent][iKt][i]->Fill(mShKv, -std::imag(mShYlmBuffer[i]));
        }
        // covariance: outer product of the (re, -im) Ylm vector packed on 2*nJM axes
        // (each Ylm contributes two consecutive axis bins: even = real, odd = -imag)
        static constexpr int ComponentsPerLM = 2;
        const int nAxisLM = ComponentsPerLM * static_cast<int>(mShYlmBuffer.size());
        for (int iz = 0; iz < nAxisLM; ++iz) {
          const double vz = (iz % ComponentsPerLM == 0) ? std::real(mShYlmBuffer[iz / ComponentsPerLM]) : -std::imag(mShYlmBuffer[iz / ComponentsPerLM]);
          for (int ip = 0; ip < nAxisLM; ++ip) {
            const double vp = (ip % ComponentsPerLM == 0) ? std::real(mShYlmBuffer[ip / ComponentsPerLM]) : -std::imag(mShYlmBuffer[ip / ComponentsPerLM]);
            mShCov[iCent][iKt]->Fill(mShKv, static_cast<double>(iz), static_cast<double>(ip), vz * vp);
          }
        }

        mShBinCount[iCent][iKt]->Fill(mShKv, 1.0);
        if (mShPlot1D) {
          // FemtoUniverse h1D = f3d[0]: qinv (=2k*) for identical-LCMS, else k*.
          const float sh1DValue = (mShFrame == ShFrameLcmsIdentical) ? (2.0f * mKstar) : mKstar;
          mSh1D[iCent][iKt]->Fill(sh1DValue);
        }
      }
    }
  }

  // reco-vs-truth correlation fill (kReco and kMc both set)
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
    if (mHasMcPair && mPlotBertschPratt) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kTrueQoutVsQout, HistTable)), mTrueQout, mQout);
      mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kTrueQsideVsQside, HistTable)), mTrueQside, mQside);
      mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kTrueQlongVsQlong, HistTable)), mTrueQlong, mQlong);
    }
  }

  // pure mc-truth pair fill (kMc without kReco) — uses the same mTrue* members
  // that setPairMcTruth() populated; no mHasMcPair/mHasMcCol gating needed since
  // there is no missing-link case here (the particles ARE the truth particles)
  void fillMcTruth()
  {
    if (mPlot1d) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kTrueKstar, HistTable)), mTrueKstar);
      mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kTrueKt, HistTable)), mTrueKt);
      mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kTrueMt, HistTable)), mTrueMt);
      mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kTrueMinv1D, HistTable)), mTrueMinv);
    }
    if (mPlot2d) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kTruePt1VsPt2, HistTable)), mTrueParticle1.Pt(), mTrueParticle2.Pt());
      mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kTruePt1VsKstar, HistTable)), mTrueParticle1.Pt(), mTrueKstar);
      mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kTruePt2VsKstar, HistTable)), mTrueParticle2.Pt(), mTrueKstar);
      mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kTruePt1VsKt, HistTable)), mTrueParticle1.Pt(), mTrueKt);
      mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kTruePt2VsKt, HistTable)), mTrueParticle2.Pt(), mTrueKt);
      mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kTruePt1VsMt, HistTable)), mTrueParticle1.Pt(), mTrueMt);
      mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kTruePt2VsMt, HistTable)), mTrueParticle2.Pt(), mTrueMt);
      mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kTrueKstarVsKt, HistTable)), mTrueKstar, mTrueKt);
      mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kTrueKstarVsMt, HistTable)), mTrueKstar, mTrueMt);
      mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kTrueKstarVsMult, HistTable)), mTrueKstar, mTrueMult);
      mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kTrueKstarVsCent, HistTable)), mTrueKstar, mTrueCent);
    }
    if (mPlotDeltaEtaDeltaPhi) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kTrueDeltaEtaDeltaPhi, HistTable)), mTrueDeltaPhi, mTrueDeltaEta);
    }
    if (mPlotBertschPratt) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kTrueQout, HistTable)), mTrueQout);
      mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kTrueQside, HistTable)), mTrueQside);
      mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kTrueQlong, HistTable)), mTrueQlong);
      mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kTrueQoutQsideQlong, HistTable)), mTrueQout, mTrueQside, mTrueQlong);
    }
    if (mPlotKstarVsMtVsMult) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kTrueKstarVsMtVsMult, HistTable)), mTrueKstar, mTrueMt, mTrueMult);
    }
    if (mPlotKstarVsMtVsMultVsCent) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kTrueKstarVsMtVsMultVsCent, HistTable)), mTrueKstar, mTrueMt, mTrueMult, mTrueCent);
    }
    if (mPlotKstarVsMtVsPt1VsPt2) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kTrueKstarVsMtVsPt1VsPt2, HistTable)), mTrueKstar, mTrueMt, mTrueParticle1.Pt(), mTrueParticle2.Pt());
    }
    if (mPlotKstarVsMtVsPt1VsPt2VsMult) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kTrueKstarVsMtVsPt1VsPt2VsMult, HistTable)), mTrueKstar, mTrueMt, mTrueParticle1.Pt(), mTrueParticle2.Pt(), mTrueMult);
    }
    if (mPlotKstarVsMtVsPt1VsPt2VsMultVsCent) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kTrueKstarVsMtVsPt1VsPt2VsMultVsCent, HistTable)), mTrueKstar, mTrueMt, mTrueParticle1.Pt(), mTrueParticle2.Pt(), mTrueMult, mTrueCent);
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
    // Use Cartesian 4-vectors: addition/M2() become pure arithmetic
    const ROOT::Math::PxPyPzEVector p1(part1);
    const ROOT::Math::PxPyPzEVector p2(part2);

    // Mandelstam s = (p1 + p2)^2
    const double s = (p1 + p2).M2();
    const double m1sq = p1.M2();
    const double m2sq = p2.M2();

    // Källen function λ(s, m1^2, m2^2) = (s - m1^2 - m2^2)² - 4*m1^2*m2^2
    const double kallen = (s - m1sq - m2sq) * (s - m1sq - m2sq) - 4.0 * m1sq * m2sq;

    // k* = 0.5 * sqrt(λ/s)
    return static_cast<float>(0.5 * std::sqrt(std::max(0.0, kallen) / s));
  }

  std::tuple<float, float, float> computeBertschPrattLCMS(ROOT::Math::PtEtaPhiMVector const& part1, ROOT::Math::PtEtaPhiMVector const& part2)
  {
    const ROOT::Math::PxPyPzEVector p1(part1);
    const ROOT::Math::PxPyPzEVector p2(part2);
    const ROOT::Math::PxPyPzEVector pSum = p1 + p2;

    const double tPx = pSum.Px();
    const double tPy = pSum.Py();
    const double tPz = pSum.Pz();
    const double tE = pSum.E();

    const double tPt = std::sqrt(tPx * tPx + tPy * tPy);
    const double tMt = std::sqrt(tE * tE - tPz * tPz);

    static constexpr double MinTransverseMomentum = 1e-9;
    if (tPt < MinTransverseMomentum || tMt < MinTransverseMomentum) {
      return {0.0f, 0.0f, 0.0f};
    }

    const double betaL = tPz / tE;
    const double gammaL = tE / tMt;

    const double kout1 = (p1.Px() * tPx + p1.Py() * tPy) / tPt;
    const double kside1 = (-p1.Px() * tPy + p1.Py() * tPx) / tPt;
    const double klong1 = gammaL * (p1.Pz() - betaL * p1.E());

    const double kout2 = (p2.Px() * tPx + p2.Py() * tPy) / tPt;
    const double kside2 = (p2.Py() * tPx - p2.Px() * tPy) / tPt;
    const double klong2 = gammaL * (p2.Pz() - betaL * p2.E());

    auto qOut = static_cast<float>(kout1 - kout2);
    auto qSide = static_cast<float>(kside1 - kside2);
    auto qLong = static_cast<float>(klong1 - klong2);

    return {qOut, qSide, qLong};
  }

  // Return the bin index for value given ascending bin edges, or -1 if out of range.
  // edges = {e0, e1, ..., eN} defines N bins [e0,e1), [e1,e2), ..., [e_{N-1},eN).
  static int findShBin(double value, std::vector<double> const& edges)
  {
    static constexpr std::size_t MinEdgesForOneBin = 2;
    if (edges.size() < MinEdgesForOneBin || value < edges.front() || value >= edges.back()) {
      return -1;
    }
    for (std::size_t i = 0; i < edges.size() - 1; ++i) {
      if (value >= edges[i] && value < edges[i + 1]) {
        return static_cast<int>(i);
      }
    }
    return -1;
  }

  // Kinematics feeding the spherical-harmonics decomposition, ported 1:1 from
  // FemtoUniverseMath::newpairfunc so results are numerically comparable.
  // Returns {kv, out, side, long}, selected by mShFrame:
  //   ShFrameLcmsNonIdentical: {kstar,  fDKOut,  fDKSide,  fDKLong}   (LCMS components of particle 1)
  //   ShFrameLcmsIdentical:    {qinv,   outLCMS, sideLCMS, longLCMS}  (LCMS q-differences)
  //   ShFramePrf:              {|qPRF|, outPRF,  sidePRF,  longPRF}   (matches FemtoUniverse isIdenPRF=true)
  std::tuple<float, float, float, float> computeShKinematics(ROOT::Math::PtEtaPhiMVector const& part1,
                                                             ROOT::Math::PtEtaPhiMVector const& part2)
  {
    const ROOT::Math::PxPyPzEVector p1(part1);
    const ROOT::Math::PxPyPzEVector p2(part2);
    const ROOT::Math::PxPyPzEVector sum = p1 + p2;

    const double tPx = sum.Px();
    const double tPy = sum.Py();
    const double tPz = sum.Pz();
    const double tE = sum.E();
    const double tPtSq = tPx * tPx + tPy * tPy;
    const double tMtSq = tE * tE - tPz * tPz;

    static constexpr double MinTransverseMomentum = 1e-9;
    if (tPtSq < MinTransverseMomentum || tMtSq < MinTransverseMomentum) {
      return {0.f, 0.f, 0.f, 0.f};
    }
    const double tPt = std::sqrt(tPtSq);
    const double tMt = std::sqrt(tMtSq);
    const double tM = std::sqrt(std::max(0.0, tMtSq - tPtSq));

    const double beta = tPz / tE;
    const double gamma = tE / tMt;

    // LCMS components of particle 1 (fDKOut/Side/Long in newpairfunc)
    const double e1 = p1.E();
    const double fDKOut = (p1.Px() * tPx + p1.Py() * tPy) / tPt;
    const double fDKSide = (-p1.Px() * tPy + p1.Py() * tPx) / tPt;
    const double fDKLong = gamma * (p1.Pz() - beta * e1);
    const double fDE = gamma * (e1 - beta * p1.Pz());

    // LCMS components of particle 2
    const double e2 = p2.E();
    const double px2 = (p2.Px() * tPx + p2.Py() * tPy) / tPt;
    const double py2 = (p2.Py() * tPx - p2.Px() * tPy) / tPt;
    const double pz2 = gamma * (p2.Pz() - beta * e2);

    const double outLCMS = fDKOut - px2;
    const double sideLCMS = fDKSide - py2;
    const double longLCMS = fDKLong - pz2;

    if (mShFrame == ShFrameLcmsIdentical) {
      // LCMS identical: Ylm <- LCMS q-differences, axis <- qinv
      const double qinv = std::sqrt(outLCMS * outLCMS + sideLCMS * sideLCMS + longLCMS * longLCMS);
      return {static_cast<float>(qinv),
              static_cast<float>(outLCMS), static_cast<float>(sideLCMS), static_cast<float>(longLCMS)};
    }
    if (mShFrame == ShFramePrf) {
      // PRF: Ylm <- PRF q-differences, axis <- |q_PRF| (matches FemtoUniverse isIdenPRF=true)
      const double pE2LCMS = gamma * (e2 - beta * p2.Pz());
      const double betaOut = tPt / tMt;
      const double gammaOut = tMt / tM;
      const double outPRF = gammaOut * (outLCMS - betaOut * (fDE - pE2LCMS));
      const double sidePRF = sideLCMS;
      const double longPRF = longLCMS;
      const double qPRF = std::sqrt(outPRF * outPRF + sidePRF * sidePRF + longPRF * longPRF);
      return {static_cast<float>(qPRF),
              static_cast<float>(outPRF), static_cast<float>(sidePRF), static_cast<float>(longPRF)};
    }
    // LCMS non-identical (default): Ylm <- LCMS components of particle 1, axis <- k*
    const double betaOut = tPt / tMt;
    const double gammaOut = tMt / tM;
    const double fKOut = gammaOut * (fDKOut - betaOut * fDE);
    const double kstar = std::sqrt(fKOut * fKOut + fDKSide * fDKSide + fDKLong * fDKLong);
    return {static_cast<float>(kstar),
            static_cast<float>(fDKOut), static_cast<float>(fDKSide), static_cast<float>(fDKLong)};
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
  ROOT::Math::PtEtaPhiMVector mParticle1;
  ROOT::Math::PtEtaPhiMVector mParticle2;
  float mRecoMass1 = 0.f;
  float mRecoMass2 = 0.f;
  float mKstar = 0.f;
  float mKt = 0.f;
  float mMt = 0.f;
  float mMassInv = 0.f;
  float mPtPair = 0.f;
  float mMult = 0.f;
  float mCent = 0.f;
  double mMass12 = 0.;
  double mMass13 = 0.;
  double mMassTot2 = 0.;

  // mc (used for both reco-vs-truth correlation AND pure mc-truth-only pairs —
  // for the latter, these are simply the primary/only kinematic values, not a
  // "true" comparison against anything)
  ROOT::Math::PtEtaPhiMVector mTrueParticle1;
  ROOT::Math::PtEtaPhiMVector mTrueParticle2;
  float mTrueKstar = 0.f;
  float mTrueKt = 0.f;
  float mTrueMt = 0.f;
  float mTrueMinv = 0.f;
  float mTrueMult = 0.f;
  float mTrueCent = 0.f;
  float mTrueQout = 0.f;
  float mTrueQside = 0.f;
  float mTrueQlong = 0.f;
  float mTrueDeltaEta = 0.f;
  float mTrueDeltaPhi = 0.f;

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

  bool mPlotKstarVsMtVsMinvVsPt1VsPt2 = false;
  bool mPlotKstarVsMtVsMinvVsPt1VsPt2VsMult = false;
  bool mPlotKstarVsMtVsMinvVsPt1VsPt2VsMultVsCent = false;

  bool mPlotKstarVsMinvVsPtPairVsMult = false;
  bool mPlotKstarVsMtVsMinvVsPtPairVsMult = false;
  bool mPlotKstarVsMtVsMinvVsPtPairVsPt1VsPt2 = false;
  bool mPlotKstarVsMtVsMinvVsPtPairVsPt1VsPt2VsMult = false;
  bool mPlotKstarVsMtVsMinvVsPtPairVsPt1VsPt2VsMultVsCent = false;

  bool mPlotDalitz = false;
  bool mPlotDeltaEtaDeltaPhi = false;

  bool mPlotBertschPratt = false;

  float mQout = 0.f;
  float mQside = 0.f;
  float mQlong = 0.f;

  float mDeltaEta = 0.f;
  float mDeltaPhi = 0.f;

  // Spherical harmonics
  bool mPlotSH = false;
  bool mShUseCent = false;
  int mShLMax = 1;
  int mShFrame = 0;
  static constexpr int ShFrameLcmsNonIdentical = 0;
  static constexpr int ShFrameLcmsIdentical = 1;
  static constexpr int ShFramePrf = 2;

  o2::framework::AxisSpec mShKstarSpec{{60, 0.0f, 0.3f}, "k* (GeV/#it{c})"}; // set in init()

  // kinematics computed in setPair(): axis value + 3 components feeding Ylm
  float mShKv = 0.f; // kstar (non-identical) or qinv (identical)
  float mShOut = 0.f;
  float mShSide = 0.f;
  float mShLong = 0.f;

  // SH histograms binned in [iCent][iKt][ihist]; ihist = l*(l+1)+m
  std::vector<std::vector<std::vector<std::shared_ptr<TH1>>>> mShReal;
  std::vector<std::vector<std::vector<std::shared_ptr<TH1>>>> mShImag;
  // SH covariance matrix per [iCent][iKt]; TH3d: k* on X, 2*nJM (l,m x re/im)
  std::vector<std::vector<std::shared_ptr<TH3>>> mShCov;
  bool mShPlot1D = false;
  std::vector<std::vector<std::shared_ptr<TH1>>> mSh1D;
  std::vector<std::vector<std::shared_ptr<TH1>>> mShBinCount;
  std::vector<std::complex<double>> mShYlmBuffer; // reused, allocated once
  std::vector<double> mShCentEdges;
  std::vector<double> mShKtEdges;

  o2::analysis::femto::SpherHarMath mYlm{};

  // qa
  bool mPairCorrelationQa = false;
  bool mEventMixingQa = false;

  std::unordered_set<int64_t> mParticles1PerEvent;
  std::unordered_set<int64_t> mParticles2PerEvent;
};
}; // namespace o2::analysis::femto::pairhistmanager
#endif // PWGCF_FEMTO_CORE_PAIRHISTMANAGER_H_
