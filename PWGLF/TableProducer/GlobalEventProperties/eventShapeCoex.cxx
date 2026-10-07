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

/// \file eventShapeCoex.cxx
/// \brief Moment sums of the forward/backward sub-event mean pT of mid-rapidity tracks
/// \author Cristian Andrei <Cristian.Andrei@cern.ch>
///
/// First O2Physics task of the AliPhysics PWGLF/SPECTRA/MultEvShape analyses.

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CommonConstants/MathConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <THn.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <memory>
#include <string>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace
{
constexpr float PtMin = 0.15f; // GeV/c
constexpr float EtaMax = 0.8f;
constexpr double MeanPtShift = 0.6; // GeV/c

// Sub-events are unions of |eta| x pT rings; the last pT ring is open-ended
constexpr int NEtaRings = 5;
constexpr std::array<float, NEtaRings + 1> EtaRingEdges = {0.f, 0.2f, 0.3f, 0.4f, 0.5f, EtaMax};
constexpr int NPtRings = 6;
constexpr std::array<float, NPtRings> PtRingLowEdges = {PtMin, 1.f, 1.5f, 2.f, 3.f, 5.f};
constexpr int EtaFirstNominal = 1; // |eta| > 0.2
constexpr int PtEndNominal = 3;    // pT < 2 GeV/c

enum Side { kF = 0,
            kB,
            NSides };

// Monomial exponents of (xF, xB, vF, vB); the order is part of the output format
constexpr int NMono1 = 16;
constexpr int NMono = 81;
constexpr std::array<std::array<uint8_t, 4>, NMono> MonoExp = {{{0, 0, 0, 0}, {1, 0, 0, 0}, {0, 1, 0, 0}, {2, 0, 0, 0}, {0, 2, 0, 0}, {1, 1, 0, 0}, {2, 1, 0, 0}, {1, 2, 0, 0}, {2, 2, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}, {0, 1, 1, 0}, {1, 0, 0, 1}, {0, 2, 1, 0}, {2, 0, 0, 1}, {0, 0, 1, 1}, {0, 0, 0, 2}, {0, 0, 1, 2}, {0, 0, 2, 0}, {0, 0, 2, 1}, {0, 0, 2, 2}, {0, 1, 0, 1}, {0, 1, 1, 1}, {0, 1, 2, 0}, {0, 1, 2, 1}, {0, 2, 0, 1}, {0, 2, 1, 1}, {0, 2, 2, 0}, {0, 2, 2, 1}, {0, 3, 0, 0}, {0, 3, 1, 0}, {0, 3, 2, 0}, {0, 4, 0, 0}, {0, 4, 1, 0}, {0, 4, 2, 0}, {1, 0, 0, 2}, {1, 0, 1, 0}, {1, 0, 1, 1}, {1, 0, 1, 2}, {1, 1, 0, 1}, {1, 1, 1, 0}, {1, 1, 1, 1}, {1, 2, 0, 1}, {1, 2, 1, 0}, {1, 2, 1, 1}, {1, 3, 0, 0}, {1, 3, 1, 0}, {1, 4, 0, 0}, {1, 4, 1, 0}, {2, 0, 0, 2}, {2, 0, 1, 0}, {2, 0, 1, 1}, {2, 0, 1, 2}, {2, 1, 0, 1}, {2, 1, 1, 0}, {2, 1, 1, 1}, {2, 2, 0, 1}, {2, 2, 1, 0}, {2, 2, 1, 1}, {2, 3, 0, 0}, {2, 3, 1, 0}, {2, 4, 0, 0}, {2, 4, 1, 0}, {3, 0, 0, 0}, {3, 0, 0, 1}, {3, 0, 0, 2}, {3, 1, 0, 0}, {3, 1, 0, 1}, {3, 2, 0, 0}, {3, 2, 0, 1}, {3, 3, 0, 0}, {3, 4, 0, 0}, {4, 0, 0, 0}, {4, 0, 0, 1}, {4, 0, 0, 2}, {4, 1, 0, 0}, {4, 1, 0, 1}, {4, 2, 0, 0}, {4, 2, 0, 1}, {4, 3, 0, 0}, {4, 4, 0, 0}}};

// Event-selection bits and the pile-up selection variants built from them
constexpr int NQaBits = 12;
constexpr std::array<o2::aod::evsel::EventSelectionFlags, NQaBits> QaBitOrder = {
  o2::aod::evsel::kNoSameBunchPileup, o2::aod::evsel::kIsGoodZvtxFT0vsPV,
  o2::aod::evsel::kIsVertexITSTPC, o2::aod::evsel::kIsVertexTOFmatched,
  o2::aod::evsel::kNoCollInTimeRangeNarrow, o2::aod::evsel::kNoCollInTimeRangeStrict,
  o2::aod::evsel::kNoCollInTimeRangeStandard, o2::aod::evsel::kNoCollInRofStrict,
  o2::aod::evsel::kNoCollInRofStandard, o2::aod::evsel::kNoHighMultCollInPrevRof,
  o2::aod::evsel::kNoITSROFrameBorder, o2::aod::evsel::kNoTimeFrameBorder};
enum EvVariant { kV1 = 0,
                 kV0,
                 kV3,
                 kV6,
                 NEvVariants };
constexpr std::array<uint16_t, NEvVariants> EvVariantMask = {
  (1u << 0) | (1u << 1),
  0u,
  (1u << 0) | (1u << 1) | (1u << 6),
  (1u << 0) | (1u << 1) | (1u << 2) | (1u << 6) | (1u << 8) | (1u << 9)};
constexpr uint16_t AllQaBits = 0xFFFF;

// Track count defining the stratum
enum Conditioning { kCondOwn = 0,      // own sub-events
                    kCondAll,          // all tracks, gap included
                    kCondUnwindowed }; // nominal gap, no upper pT cut

struct Config {
  const char* name;
  int etaFirst;      // |eta| rings [etaFirst, NEtaRings)
  int ptEnd;         // pT rings [0, ptEnd)
  Conditioning cond; // stratum count
  int validEtaFirst; // sub-events that must hold MinTracksValid tracks
  int validPtEnd;
  EvVariant evVariant; // pile-up selection
};
constexpr int NConfigs = 13;
constexpr int INominal = 0;
constexpr std::array<Config, NConfigs> Configs = {{{"NOM", 1, 3, kCondOwn, 1, 3, kV1},
                                                   {"GAP02", 1, 3, kCondAll, 4, 3, kV1},
                                                   {"GAP03", 2, 3, kCondAll, 4, 3, kV1},
                                                   {"GAP04", 3, 3, kCondAll, 4, 3, kV1},
                                                   {"GAP05", 4, 3, kCondAll, 4, 3, kV1},
                                                   {"PT10", 1, 1, kCondUnwindowed, 1, 1, kV1},
                                                   {"PT15", 1, 2, kCondUnwindowed, 1, 1, kV1},
                                                   {"PT20", 1, 3, kCondUnwindowed, 1, 1, kV1},
                                                   {"PT30", 1, 4, kCondUnwindowed, 1, 1, kV1},
                                                   {"PT50", 1, 5, kCondUnwindowed, 1, 1, kV1},
                                                   {"EVV0", 1, 3, kCondOwn, 1, 3, kV0},
                                                   {"EVV3", 1, 3, kCondOwn, 1, 3, kV3},
                                                   {"EVV6", 1, 3, kCondOwn, 1, 3, kV6}}};

constexpr int NCountBins = 64;    // larger counts in the overflow bin
constexpr int NSubCountBins = 32; // per sub-event
constexpr int MinTracksValid = 2; // per sub-event

// Data halves and subsample groups, defined on blocks of bunch crossings
constexpr int NHalves = 2;
constexpr std::array<const char*, NHalves> HalfName = {"A", "B_SEALED"};
constexpr int HalfBlockShift = 26;  // 2^26 BC = 1.7 s
constexpr int GroupBlockShift = 20; // 2^20 BC = 26 ms
constexpr int NGroups = 20;

/// splitmix64 finaliser
constexpr uint64_t mix64(uint64_t z)
{
  z += 0x9e3779b97f4a7c15ULL;
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
  z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
  return z ^ (z >> 31);
}

/// Hash of the (run, block of 2^shift bunch crossings) of an event
uint64_t blockHash(int run, uint64_t globalBC, int shift)
{
  return mix64(mix64(static_cast<uint64_t>(run)) ^ (globalBC >> shift));
}

struct SubEvent {
  int n = 0;
  double sumPt = 0.;
  double sumPt2 = 0.;
};

/// Track sums of one event per (side, |eta| ring, pT ring)
struct Rings {
  std::array<SubEvent, NSides * NEtaRings * NPtRings> ring{};

  void add(float eta, double pt)
  {
    const float absEta = std::abs(eta);
    int iEta = 0;
    while (iEta < NEtaRings - 1 && absEta > EtaRingEdges[iEta + 1]) {
      ++iEta;
    }
    int iPt = 0;
    while (iPt < NPtRings - 1 && pt >= PtRingLowEdges[iPt + 1]) {
      ++iPt;
    }
    auto& r = ring[((eta > 0.f ? kF : kB) * NEtaRings + iEta) * NPtRings + iPt];
    ++r.n;
    r.sumPt += pt;
    r.sumPt2 += pt * pt;
  }

  /// Sub-event of one side: |eta| rings [etaFirst, NEtaRings), pT rings [0, ptEnd)
  [[nodiscard]] SubEvent sum(int side, int etaFirst, int ptEnd) const
  {
    SubEvent s;
    for (int iEta = etaFirst; iEta < NEtaRings; ++iEta) {
      for (int iPt = 0; iPt < ptEnd; ++iPt) {
        const auto& r = ring[(side * NEtaRings + iEta) * NPtRings + iPt];
        s.n += r.n;
        s.sumPt += r.sumPt;
        s.sumPt2 += r.sumPt2;
      }
    }
    return s;
  }

  [[nodiscard]] int count(int etaFirst, int ptEnd) const { return sum(kF, etaFirst, ptEnd).n + sum(kB, etaFirst, ptEnd).n; }
};

/// x = <pT> - MeanPtShift, v = pT variance / n (0 for n < 2)
void observe(const SubEvent& s, double& x, double& v)
{
  x = s.sumPt / s.n - MeanPtShift;
  v = 0.;
  if (s.n > 1) {
    const double var = (s.sumPt2 - s.sumPt * s.sumPt / s.n) / (s.n - 1);
    v = std::max(var, 0.) / s.n;
  }
}

/// The first nMono monomials of (xF, xB, vF, vB)
void monomials(double xF, double xB, double vF, double vB, int nMono, std::array<double, NMono>& m)
{
  const std::array<double, 5> pxF = {1., xF, xF * xF, xF * xF * xF, xF * xF * xF * xF};
  const std::array<double, 5> pxB = {1., xB, xB * xB, xB * xB * xB, xB * xB * xB * xB};
  const std::array<double, 3> pvF = {1., vF, vF * vF};
  const std::array<double, 3> pvB = {1., vB, vB * vB};
  for (int k = 0; k < nMono; ++k) {
    const auto& e = MonoExp[k];
    m[k] = pxF[e[0]] * pxB[e[1]] * pvF[e[2]] * pvB[e[3]];
  }
}

/// Moment sums of one reconstructed half or of the generated events
struct Sums {
  std::array<std::shared_ptr<TH3>, NConfigs> config{}; // (N_fwd, stratum count, monomial)
  std::shared_ptr<THn> strata3{};                      // (N_fwd, n_F, n_B, first-order monomial)
  std::shared_ptr<TH3> occupancy{};                    // (configuration, N_fwd, stratum count)
};

/// Nominal-configuration observables of an event
struct Nominal {
  bool valid = false;
  int nMid = 0;
  double xF = 0., xB = 0., vF = 0., vB = 0.;
};

enum EventFlowBin { kFlowAll = 0,
                    kFlowSel8,
                    kFlowVtxZ,
                    kFlowPileupV1,
                    kFlowNominalValid,
                    kFlowHalfA,
                    kFlowHalfB,
                    NEventFlowBins };
constexpr std::array<const char*, NEventFlowBins> EventFlowLabels = {"all", "sel8", "|z|", "pile-up V1", "nominal valid", "half A", "half B"};
} // namespace

struct EventShapeCoex {
  HistogramRegistry registry{"registry"};
  Service<o2::framework::O2DatabasePDG> pdg{};

  Configurable<float> vtxZCut{"vtxZCut", 10.0f, "Max |PV z| (cm)"};
  Configurable<float> etaFwdAMin{"etaFwdAMin", 3.5f, "FT0-A proxy: lower eta edge (generator level)"};
  Configurable<float> etaFwdAMax{"etaFwdAMax", 4.9f, "FT0-A proxy: upper eta edge (generator level)"};
  Configurable<float> etaFwdCMin{"etaFwdCMin", -3.3f, "FT0-C proxy: lower eta edge (generator level)"};
  Configurable<float> etaFwdCMax{"etaFwdCMax", -2.1f, "FT0-C proxy: upper eta edge (generator level)"};
  ConfigurableAxis axisNFwd{"axisNFwd", {50, 1., 30000.}, "N_fwd = FT0A + FT0C amplitude (a fixed-width setting is made logarithmic)"};
  ConfigurableAxis axisNFwdGen{"axisNFwdGen", {50, 1., 1000.}, "generator-level N_fwd = charged primaries in the FT0 acceptance (a fixed-width setting is made logarithmic)"};
  ConfigurableAxis axisRun{"axisRun", {60000, 519999.5, 579999.5}, "run number (QA)"};

  Filter trackFilter = (aod::track::pt > PtMin) && (nabs(aod::track::eta) < EtaMax) && requireGlobalTrackInFilter();

  using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults>;
  using MyCollisionsMc = soa::Join<MyCollisions, aod::McCollisionLabels>;
  using MyTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection>>;

  std::array<Sums, NHalves> sumsReco{};
  std::array<std::shared_ptr<TH3>, NHalves> designEffect{}; // (N_fwd, group, first-order monomial)
  Sums sumsGen{};

  /// Books the moment sums of one level under dir/<configuration>/<half>
  void bookSums(Sums& sums, const std::string& dir, const std::string& half, const AxisSpec& nFwdAxis, int nConfigs)
  {
    const AxisSpec countAxis{NCountBins, -0.5, NCountBins - 0.5, "stratum count"};
    const AxisSpec subCountAxis{NSubCountBins, -0.5, NSubCountBins - 0.5, "sub-event count"};
    const AxisSpec monoAxis{NMono, -0.5, NMono - 0.5, "monomial"};
    const AxisSpec mono1Axis{NMono1, -0.5, NMono1 - 0.5, "monomial"};
    const AxisSpec configAxis{NConfigs, -0.5, NConfigs - 0.5, "configuration"};
    const std::string label = "[half " + half + "] moment sums ";
    for (int c = 0; c < nConfigs; ++c) {
      sums.config[c] = registry.add<TH3>(dir + "/" + Configs[c].name + "/" + half, (label + Configs[c].name).c_str(), HistType::kTH3D, {nFwdAxis, countAxis, monoAxis});
      sums.config[c]->SetBit(TH1::kIsNotW); // no Sumw2
    }
    sums.strata3 = registry.add<THn>(dir + "/NOM3/" + half, (label + "NOM in (n_F, n_B) strata").c_str(), HistType::kTHnD, {nFwdAxis, subCountAxis, subCountAxis, mono1Axis});
    sums.occupancy = registry.add<TH3>(dir + "/occupancy/" + half, "event counts per configuration", HistType::kTH3D, {configAxis, nFwdAxis, countAxis});
  }

  void init(InitContext const&)
  {
    if (doprocessReco && doprocessRecoMC) {
      LOGF(fatal, "processReco and processRecoMC both fill the reconstructed-level sums; enable only one of them.");
    }

    if (doprocessReco || doprocessRecoMC) {
      AxisSpec nFwdAxis{axisNFwd, "#it{N}_{fwd} (FT0A+FT0C amplitude)"};
      if (nFwdAxis.nBins) {
        nFwdAxis.makeLogarithmic();
      }
      const AxisSpec groupAxis{NGroups, -0.5, NGroups - 0.5, "group"};
      const AxisSpec mono1Axis{NMono1, -0.5, NMono1 - 0.5, "monomial"};
      for (int h = 0; h < NHalves; ++h) {
        bookSums(sumsReco[h], "accum", HalfName[h], nFwdAxis, NConfigs);
        designEffect[h] = registry.add<TH3>(std::string("designEffect/") + HalfName[h], "nominal first-order sums per block group", HistType::kTH3D, {nFwdAxis, groupAxis, mono1Axis});
        designEffect[h]->SetBit(TH1::kIsNotW);
      }

      // QA
      constexpr int NFlowBins = NEventFlowBins;
      auto hFlow = registry.add<TH1>("qa/hEventFlow", "events", HistType::kTH1D, {{NFlowBins, -0.5, NFlowBins - 0.5}});
      for (int b = 0; b < NFlowBins; ++b) {
        hFlow->GetXaxis()->SetBinLabel(b + 1, EventFlowLabels[b]);
      }
      registry.add("qa/hQaBits", "events with qaBits bit i set (order: QaBitOrder)", HistType::kTH1D, {{NQaBits, -0.5, NQaBits - 0.5}});
      registry.add("qa/hNFwd", "selected events", HistType::kTH1D, {nFwdAxis});
      AxisSpec ptAxis{100, PtMin, 20., "#it{p}_{T} (GeV/#it{c})"};
      ptAxis.makeLogarithmic();
      registry.add("qa/hTrackPt", "selected tracks", HistType::kTH1D, {ptAxis});
      registry.add("qa/hTrackEta", "selected tracks", HistType::kTH1D, {{160, -EtaMax, EtaMax, "#eta"}});
      registry.add("qa/hTrackPhi", "selected tracks", HistType::kTH1D, {{72, 0., o2::constants::math::TwoPI, "#varphi"}});
      // quantities: count, a_F, a_B, v_F, v_B
      const AxisSpec quantityAxis{5, -0.5, 4.5, "quantity"};
      registry.add<TH3>("qa/hMarginalSums", "nominal single-side sums", HistType::kTH3D, {nFwdAxis, {NCountBins, -0.5, NCountBins - 0.5, "n_{mid}"}, quantityAxis})->SetBit(TH1::kIsNotW);
      registry.add<TH2>("qa/hRunSums", "nominal single-side sums per run", HistType::kTH2D, {{axisRun, "run"}, quantityAxis})->SetBit(TH1::kIsNotW);
    }

    if (doprocessMC) {
      AxisSpec nFwdGenAxis{axisNFwdGen, "#it{N}_{fwd} (charged primaries in FT0A+FT0C)"};
      if (nFwdGenAxis.nBins) {
        nFwdGenAxis.makeLogarithmic();
      }
      bookSums(sumsGen, "gen", "all", nFwdGenAxis, 1); // nominal configuration only
      registry.add("gen/hEvents", "generated events", HistType::kTH1D, {{1, -0.5, 0.5}});
    }
  }

  /// Fills the sums of every valid configuration; returns the nominal observables
  Nominal fillSums(Sums& sums, const Rings& rings, double nFwd, uint16_t qaBits, int nConfigs)
  {
    Nominal nominal;
    std::array<double, NMono> m{};
    double xF = 0., xB = 0., vF = 0., vB = 0.;
    for (int c = 0; c < nConfigs; ++c) {
      const auto& cfg = Configs[c];
      const uint16_t mask = EvVariantMask[cfg.evVariant];
      if ((qaBits & mask) != mask) {
        continue;
      }
      if (rings.sum(kF, cfg.validEtaFirst, cfg.validPtEnd).n < MinTracksValid || rings.sum(kB, cfg.validEtaFirst, cfg.validPtEnd).n < MinTracksValid) {
        continue;
      }
      const SubEvent fwd = rings.sum(kF, cfg.etaFirst, cfg.ptEnd);
      const SubEvent bwd = rings.sum(kB, cfg.etaFirst, cfg.ptEnd);
      int nCond = fwd.n + bwd.n;
      if (cfg.cond == kCondAll) {
        nCond = rings.count(0, PtEndNominal);
      } else if (cfg.cond == kCondUnwindowed) {
        nCond = rings.count(EtaFirstNominal, NPtRings);
      }
      observe(fwd, xF, vF);
      observe(bwd, xB, vB);
      monomials(xF, xB, vF, vB, NMono, m);
      for (int k = 0; k < NMono; ++k) {
        sums.config[c]->Fill(nFwd, nCond, k, m[k]);
      }
      sums.occupancy->Fill(c, nFwd, nCond);
      if (c == INominal) {
        nominal = {true, nCond, xF, xB, vF, vB};
      }
    }

    // (n_F, n_B) strata, n_F, n_B >= 1
    const uint16_t maskNominal = EvVariantMask[Configs[INominal].evVariant];
    const SubEvent fwd = rings.sum(kF, EtaFirstNominal, PtEndNominal);
    const SubEvent bwd = rings.sum(kB, EtaFirstNominal, PtEndNominal);
    if ((qaBits & maskNominal) == maskNominal && fwd.n > 0 && bwd.n > 0) {
      observe(fwd, xF, vF);
      observe(bwd, xB, vB);
      monomials(xF, xB, vF, vB, NMono1, m);
      std::array<double, 4> coord = {nFwd, static_cast<double>(fwd.n), static_cast<double>(bwd.n), 0.};
      for (int k = 0; k < NMono1; ++k) {
        coord[3] = k;
        sums.strata3->Fill(coord.data(), m[k]);
      }
    }
    return nominal;
  }

  template <typename TCollision>
  void fillReco(TCollision const& coll, MyTracks const& tracks)
  {
    registry.fill(HIST("qa/hEventFlow"), kFlowAll);
    if (!coll.sel8()) {
      return;
    }
    registry.fill(HIST("qa/hEventFlow"), kFlowSel8);
    if (std::abs(coll.posZ()) > vtxZCut) {
      return;
    }
    registry.fill(HIST("qa/hEventFlow"), kFlowVtxZ);

    uint16_t qaBits = 0;
    for (int i = 0; i < NQaBits; ++i) {
      if (coll.selection_bit(QaBitOrder[i])) {
        qaBits |= static_cast<uint16_t>(1u << i);
        registry.fill(HIST("qa/hQaBits"), i);
      }
    }
    if ((qaBits & EvVariantMask[kV1]) == EvVariantMask[kV1]) {
      registry.fill(HIST("qa/hEventFlow"), kFlowPileupV1);
    }

    Rings rings;
    for (const auto& track : tracks) {
      rings.add(track.eta(), track.pt());
      registry.fill(HIST("qa/hTrackPt"), track.pt());
      registry.fill(HIST("qa/hTrackEta"), track.eta());
      registry.fill(HIST("qa/hTrackPhi"), track.phi());
    }

    const auto& bc = coll.template bc_as<aod::BCs>();
    const double nFwd = coll.multFT0A() + coll.multFT0C();
    registry.fill(HIST("qa/hNFwd"), nFwd);
    const int half = static_cast<int>(blockHash(bc.runNumber(), bc.globalBC(), HalfBlockShift) & 1u);
    const Nominal nominal = fillSums(sumsReco[half], rings, nFwd, qaBits, NConfigs);
    if (!nominal.valid) {
      return;
    }
    registry.fill(HIST("qa/hEventFlow"), kFlowNominalValid);
    registry.fill(HIST("qa/hEventFlow"), kFlowHalfA + half);

    // subsample groups
    const int group = static_cast<int>(blockHash(bc.runNumber(), bc.globalBC(), GroupBlockShift) % NGroups);
    std::array<double, NMono> m{};
    monomials(nominal.xF, nominal.xB, nominal.vF, nominal.vB, NMono1, m);
    for (int k = 0; k < NMono1; ++k) {
      designEffect[half]->Fill(nFwd, group, k, m[k]);
    }

    const std::array<double, 5> quantities = {1., nominal.xF + MeanPtShift, nominal.xB + MeanPtShift, nominal.vF, nominal.vB};
    for (int q = 0; q < static_cast<int>(quantities.size()); ++q) {
      registry.fill(HIST("qa/hMarginalSums"), nFwd, nominal.nMid, q, quantities[q]);
      registry.fill(HIST("qa/hRunSums"), bc.runNumber(), q, quantities[q]);
    }
  }

  void processReco(MyCollisions::iterator const& coll, MyTracks const& tracks, aod::BCs const&)
  {
    fillReco(coll, tracks);
  }
  PROCESS_SWITCH(EventShapeCoex, processReco, "Reconstructed-level moment sums", true);

  void processRecoMC(MyCollisionsMc::iterator const& coll, MyTracks const& tracks, aod::BCs const&)
  {
    fillReco(coll, tracks);
  }
  PROCESS_SWITCH(EventShapeCoex, processRecoMC, "Reconstructed-level moment sums on MC (MC-labelled collisions)", false);

  void processMC(aod::McCollision const&, aod::McParticles const& particles)
  {
    static constexpr double ChargeTolerance = 1.e-3;

    registry.fill(HIST("gen/hEvents"), 0);
    uint32_t nFwdA = 0, nFwdC = 0;
    Rings rings;
    for (const auto& particle : particles) {
      if (!particle.isPhysicalPrimary()) {
        continue;
      }
      const auto* pdgParticle = pdg->GetParticle(particle.pdgCode());
      if (pdgParticle == nullptr || std::abs(pdgParticle->Charge()) < ChargeTolerance) {
        continue; // unknown PDG code or neutral
      }
      const float eta = particle.eta();
      // forward counts: no pT window
      if (eta > etaFwdAMin && eta < etaFwdAMax) {
        ++nFwdA;
      }
      if (eta > etaFwdCMin && eta < etaFwdCMax) {
        ++nFwdC;
      }
      if (particle.pt() <= PtMin || std::abs(eta) >= EtaMax) {
        continue;
      }
      rings.add(eta, particle.pt());
    }
    fillSums(sumsGen, rings, nFwdA + nFwdC, AllQaBits, 1);
  }
  PROCESS_SWITCH(EventShapeCoex, processMC, "Generator-level moment sums (nominal configuration)", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<EventShapeCoex>(cfgc)};
}
