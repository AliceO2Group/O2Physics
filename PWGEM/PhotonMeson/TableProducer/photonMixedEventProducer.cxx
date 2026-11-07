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

/// \file photonMixedEventProducer.cxx
/// \brief This task implements a new mixed event method and is the necessary step between builder run 1 and run 2.
//         At first it pools of the run-1 events (photons and free legs,copies with their legs) in z x occupancy bins;
//         per event one same-event artificial event (the complete event) and cfgOverlay mixed artificial events
//         (the complete event + one complete partner event of its pool, translated to the
///        current vertex), written as PCMMixCollisions / PCMMixLegs for photon-conversion-builder-mixed. Needs builder run 1 with
///        fillLegCov = true and fillFreeLegs = true.
/// \author Stefanie Mrozinski

#include "PWGEM/PhotonMeson/DataModel/EventTables.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"

#include <CommonConstants/MathConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/DataTypes.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <functional>
#include <map>
#include <memory>
#include <numeric>
#include <random>
#include <set>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;

namespace
{
constexpr int kMaxAncestorDepth = 100; // mother-chain steps searched for the physical primary

// copy of one leg
struct LegCopy {
  int8_t sign = 0;
  float x = 0.f, alpha = 0.f, y = 0.f, z = 0.f, snp = 0.f, tgl = 0.f, signed1Pt = 0.f, pt = 0.f, phi = 0.f;
  std::array<float, 15> cov{};
  uint8_t tpcNClsFindable = 0;
  int8_t tpcNClsFindableMinusFound = 0, tpcNClsFindableMinusCrossedRows = 0;
  uint8_t tpcNClsShared = 0;
  float tpcChi2NCl = 0.f, tpcInnerParam = 0.f, tpcSignal = 0.f, tpcNSigmaEl = 0.f, tpcNSigmaPi = 0.f;
  uint32_t itsClusterSizes = 0;
  float itsChi2NCl = 0.f;
  uint8_t detectorMap = 0;
  int origTrackId = -1, mcParticleId = -1, mcMotherId = -1, mcPdg = 0, mcMotherPdg = 0, mcGrandMotherId = -1, mcGrandMotherPdg = 0, mcAncestorId = -1;
  int mcGGMotherId = -1, mcGGMotherPdg = 0, mcRootId = -1, mcRootPdg = 0, mcRootDepth = -1;
  float mcMotherPx = 0.f, mcMotherPy = 0.f, mcMotherPz = 0.f;
  float mcVx = 0.f, mcVy = 0.f, mcVz = 0.f; // MC production vertex of the leg = true conversion point
};

struct PhotonCopy {
  int64_t eventUid = -1, photonUid = -1;
  int v0Id = -1;
  float px = 0.f, py = 0.f, pz = 0.f, vx = 0.f, vy = 0.f, vz = 0.f;
  float pvX = 0.f, pvY = 0.f, pvZ = 0.f, pvCovXX = 0.f, pvCovXY = 0.f, pvCovXZ = 0.f, pvCovYY = 0.f, pvCovYZ = 0.f, pvCovZZ = 0.f, pvChi2 = 0.f; // primary vertex of its event
  uint16_t pvNContrib = 0;
  std::vector<LegCopy> legs;
  int timesUsed = 0; // in mixed artificial events
};

struct EventCopy {
  int64_t eventUid = -1;
  float x = 0.f, y = 0.f, z = 0.f, occupancy = 0.f, cent = 0.f;
  std::vector<PhotonCopy> photons;
  std::vector<PhotonCopy> freeLegs;
};

// vertex and centralities of the event whose vertex an artificial event gets
struct VertexCopy {
  float x = 0.f, y = 0.f, z = 0.f, covXX = 0.f, covXY = 0.f, covXZ = 0.f, covYY = 0.f, covYZ = 0.f, covZZ = 0.f, chi2 = 0.f;
  uint16_t numContrib = 0;
  float centFT0M = 0.f, centFT0A = 0.f, centFT0C = 0.f;
  int runNumber = 0, multNTracksPVeta1 = 0; // run and multiplicities of the event whose vertex it is (photonhbt: CCDB, 2D CF axes)
  float multFT0M = 0.f;
  uint64_t timestamp = 0; // CCDB timestamp
};
} // namespace

struct PhotonMixedEventProducer {
  using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::PMEvSels, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::PVMults, aod::FT0Mults>;
  using MyTracksMC = soa::Join<aod::TracksIU, aod::McTrackLabels>;
  using MyEvents = soa::Join<aod::PMEvents, aod::EMEventsXY, aod::EMEventsCov, aod::EMEventsMult_000, aod::EMEventsCent_000>;
  using MyPhotonsDerived = soa::Join<aod::V0PhotonsKF, aod::V0KFEMEventIds>;
  using MyLegsMCDerived = soa::Join<aod::V0Legs, aod::V0LegMCLabels>;

  Produces<aod::PCMMixCollisions> mixCollisions;
  Produces<aod::PCMMixLegs> mixLegs;

  Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "centrality estimator: 0 FT0M, 1 FT0A, 2 FT0C"};
  Configurable<float> cfgCentMin{"cfgCentMin", 10.f, "min centrality (one class)"};
  Configurable<float> cfgCentMax{"cfgCentMax", 20.f, "max centrality"};
  Configurable<float> cfgZMax{"cfgZMax", 10.f, "max |z_vtx| (cm)"};
  Configurable<float> cfgZBinWidth{"cfgZBinWidth", 2.f, "pool z bin width (cm)"};
  Configurable<int> cfgOccupancyEstimator{"cfgOccupancyEstimator", 0, "occupancy: 0 trackOccupancyInTimeRange, 1 ft0cOccupancyInTimeRange"};
  Configurable<std::vector<float>> cfgOccupancyBins{"cfgOccupancyBins", {-1.f, 1.e9f}, "pool occupancy bin edges (events outside are not used)"};
  Configurable<int> cfgPoolDepth{"cfgPoolDepth", 50, "events per pool (has to exceed the photon multiplicity N)"};
  Configurable<int> cfgNMixPerEvent{"cfgNMixPerEvent", 1, "mixed artificial events per event"};
  Configurable<bool> cfgRequireFullPool{"cfgRequireFullPool", true, "skip an event completely (same and mixed) if its pool has fewer than N events; false: write the same-event one and no mixed one"};
  Configurable<int> cfgSeed{"cfgSeed", 12345, "random seed (choice of partner events and photons)"};
  Configurable<bool> cfgAlignToVertex{"cfgAlignToVertex", true, "translate the legs and the conversion point of each partner photon by (current vertex - partner vertex): every photon of the mixed artificial event originates from its vertex, as in the same event"};
  Configurable<bool> cfgPoolGroups{"cfgPoolGroups", false, "pool conflict groups instead of photons (run 1 with deduplicationMode 3 = keep all): candidates sharing a track form one group with all its tracks once; the decision is taken only in run 2, after mixing, identically for the same and the mixed event"};
  Configurable<bool> cfgPhotonWeighted{"cfgPhotonWeighted", true, "partner events drawn with probability proportional to their photon number (each pooled photon equally likely, still N distinct events); false: events equally likely"};
  Configurable<int> cfgInject{"cfgInject", 0, "injection mixed event: instead of N photons from N partner events, K mixed artificial events per event, each = the COMPLETE current event (its photons and, with cfgPoolFreeLegs, its free legs) + ONE photon of a distinct partner event (translated with cfgAlignToVertex). K = cfgInject, -1: K = N (mixed pairs = 2 x same-event pairs), -2: K = N/2 rounded up (mixed pairs = same-event pairs, half the CPU); the pool needs K partner events. Run 2 with mixed_group.skipSameSourceCrossed: the current event's pairings stay as run 1 decided them, only the injected legs form new candidates (also with the free legs, which the N-photon mixed event cannot offer)"};
  Configurable<int> cfgOverlay{"cfgOverlay", 1, "two-event leg overlay: K mixed artificial events per event, each = the COMPLETE current event + the COMPLETE partner event k (all its photons and, with cfgPoolFreeLegs, all its free legs, translated with cfgAlignToVertex). Every leg configuration of the same event exists across the two events, incl. true pairings that failed a run-1 cut. Run 2 with mixed_group.skipSameSourceCrossed: pairings inside one source event stay as run 1 decided them, only cross-event pairings are new. 0 = off (not with cfgInject)"};
  Configurable<bool> cfgPoolFreeLegs{"cfgPoolFreeLegs", true, "pool the free legs of run 1 (PCMFreeLegs, builder fillFreeLegs) with the photons: the same artificial event gets all M free legs of the event, the mixed one M free legs drawn from the N partner events, so that run 2 pairs free legs identically in both"};

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  std::map<std::pair<int, int>, std::deque<EventCopy>> mPools;
  std::mt19937 mRandom;
  int mRunNumber = -1;
  int64_t mEventUid = 0;
  int64_t mPhotonUid = 0;
  int64_t mNPooledPhotons = 0;

  enum CutFlow { kAll = 0,
                 kSelected,
                 kCent,
                 kZ,
                 kOccupancy,
                 kHasPhotons,
                 kPoolTooSmall,
                 kWritten };

  void init(InitContext&)
  {
    mRandom.seed(cfgSeed.value);
    if (cfgOccupancyBins->size() < 2) {
      LOG(fatal) << "cfgOccupancyBins needs at least two edges";
    }
    const AxisSpec axType{2, -0.5f, 1.5f, "artificial event (0 same, 1 mixed)"};
    const AxisSpec axN{101, -0.5f, 100.5f, "photons per event N"};
    const int nZ = static_cast<int>(std::lround(2.f * cfgZMax / cfgZBinWidth));
    registry.add("Event/hCutFlow", "events;0 all, 1 selected, 2 centrality, 3 z, 4 occupancy, 5 with photons, 6 pool too small (skipped), 7 written", kTH1D, {{8, -0.5f, 7.5f}});
    registry.add("Event/hNPhotons", "accepted events;photons per event N", kTH1D, {axN});
    registry.add("Event/hNFreeLegs", "accepted events (cfgPoolFreeLegs);free legs per event M", kTH1D, {{300, -0.5f, 299.5f}});
    registry.add("Mix/hFreeLegsMixedVsCurrent", "cfgPoolFreeLegs;M of the current event;free legs in the mixed artificial event", kTH2D, {{300, -0.5f, 299.5f}, {300, -0.5f, 299.5f}});
    registry.add("Event/hNPhotonsMixedVsCurrent", "mixed artificial events (must be the diagonal);N of the current event;photons in the mixed event", kTH2D, {axN, axN});
    registry.add("Event/hZ", "accepted events;z_{vtx} (cm)", kTH1D, {{nZ, -cfgZMax, cfgZMax}});
    registry.add("Event/hOccupancy", "accepted events;occupancy", kTH1D, {{200, 0.f, 20000.f}});
    registry.add("Pool/hDepthVsZ", "pool size when used;z_{vtx} (cm);events in the pool", kTH2D, {{nZ, -cfgZMax, cfgZMax}, {cfgPoolDepth + 1, -0.5f, cfgPoolDepth + 0.5f}});
    registry.add("Pool/hDepthVsOccBin", "pool size when used;occupancy bin;events in the pool", kTH2D, {{static_cast<int>(cfgOccupancyBins->size()) - 1, -0.5f, cfgOccupancyBins->size() - 1.5f}, {cfgPoolDepth + 1, -0.5f, cfgPoolDepth + 0.5f}});
    registry.add("Pool/hPartnersVsN", "partner events available vs needed;N of the current event;events in the pool", kTH2D, {axN, {cfgPoolDepth + 1, -0.5f, cfgPoolDepth + 0.5f}});
    registry.add("Pool/hNPooledPhotons", "photons in all pools (memory), per event;photons", kTH1D, {{500, 0.f, 5.e5f}});
    registry.add("Pool/hNPools", "filled pools, per event;pools", kTH1D, {{1000, 0.f, 10000.f}});
    registry.add("Pool/hResets", "pool resets (run change);resets", kTH1D, {{1, -0.5f, 0.5f}});
    registry.add("Pool/hTimesUsed", "pooled photons leaving the pool;times used in a mixed event", kTH1D, {{51, -0.5f, 50.5f}});
    registry.add("Mix/hDz", "mixed photons;z_{vtx}^{partner} - z_{vtx}^{current} (cm)", kTH1D, {{200, -0.5f, 0.5f}});
    registry.add("Mix/hDOccupancy", "mixed photons;occupancy^{partner} - occupancy^{current}", kTH1D, {{200, -5000.f, 5000.f}});
    registry.add("Mix/hDCent", "mixed photons;centrality^{partner} - centrality^{current} (%)", kTH1D, {{200, -10.f, 10.f}});
    registry.add("Mix/hPartnerAge", "mixed photons;position of the partner event in the pool (0 = newest)", kTH1D, {{cfgPoolDepth, -0.5f, cfgPoolDepth - 0.5f}});
    registry.add("Mix/hDistinctSourcesMinusN", "mixed artificial events (must be 0);distinct source events - N", kTH1D, {{21, -10.5f, 10.5f}});
    registry.add("Mix/hSourceIsCurrent", "mixed photons from the current event (must be 0);0 other event, 1 current event", kTH1D, {{2, -0.5f, 1.5f}});
    registry.add("Mix/hGroupLegs", "cfgPoolGroups: tracks per conflict group;tracks;groups", kTH1D, {{30, 0.5f, 30.5f}});
    registry.add("Event/hNCandidates", "cfgPoolGroups: run-1 candidates per event (before grouping);candidates;events", kTH1D, {{200, -0.5f, 199.5f}});
    registry.add("Mix/hLegsMinus2N", "artificial events (must be 0);artificial event;legs written - 2 N", kTH2D, {axType, {21, -10.5f, 10.5f}});
    registry.add("Mix/hPhotonPtEta", "photons of the artificial events;artificial event;p_{T} (GeV/c);#eta", kTH3D, {axType, {100, 0.f, 5.f}, {40, -1.f, 1.f}});
    registry.add("Mix/hLegTypes", "legs of the artificial events;artificial event;0 ITS-TPC, 1 TPC-only, 2 ITS-only, 3 other", kTH2D, {axType, {4, -0.5f, 3.5f}});
    // AO2D input: leg pT and phi of the track vs the values computed from V0LegsCov (what the derived-data mode uses), must be 0
    registry.add("Leg/hPtPhiFromCov", "legs (AO2D input);(p_{T}^{track} - p_{T}^{cov}) / p_{T}^{track};#varphi^{track} - #varphi^{cov} (rad)", kTH2D, {{100, -1.e-4f, 1.e-4f}, {100, -1.e-4f, 1.e-4f}});
    if (static_cast<int>(doprocessRec.value) + static_cast<int>(doprocessMC.value) + static_cast<int>(doprocessDerived.value) + static_cast<int>(doprocessDerivedMC.value) != 1) {
      LOG(fatal) << "enable exactly one of processRec, processMC (AO2D input) and processDerived, processDerivedMC (derived data)";
    }
    if (cfgOverlay != 0 && cfgInject != 0) {
      LOG(fatal) << "cfgOverlay and cfgInject are two different mixed events, enable one";
    }
    if (cfgOverlay != 0 && !cfgPoolFreeLegs) {
      LOG(warning) << "cfgOverlay without cfgPoolFreeLegs: the partner event brings its stored photons only, not its free legs";
    }
    if (cfgPoolFreeLegs && (doprocessDerived || doprocessDerivedMC)) {
      LOG(fatal) << "cfgPoolFreeLegs needs PCMFreeLegs, which are not in the derived data";
    }
  }

  template <typename TCollision>
  float centrality(TCollision const& collision) const
  {
    switch (cfgCentEstimator.value) {
      case 0:
        return collision.centFT0M();
      case 1:
        return collision.centFT0A();
      default:
        return collision.centFT0C();
    }
  }

  template <typename TCollision>
  float occupancy(TCollision const& collision) const
  {
    return cfgOccupancyEstimator.value == 1 ? collision.ft0cOccupancyInTimeRange() : static_cast<float>(collision.trackOccupancyInTimeRange());
  }

  int occupancyBin(const float occ) const
  {
    const auto& edges = cfgOccupancyBins.value;
    for (size_t i = 0; i + 1 < edges.size(); ++i) {
      if (edges[i] <= occ && occ < edges[i + 1]) {
        return static_cast<int>(i);
      }
    }
    return -1;
  }

  // leg pT and phi from the fit-input parameters
  template <typename TLegCov>
  static std::pair<float, float> ptPhiOf(TLegCov const& legCov)
  {
    const float pt = legCov.signed1Pt() != 0.f ? 1.f / std::fabs(legCov.signed1Pt()) : 0.f;
    float phi = std::asin(std::clamp(legCov.snp(), -1.f, 1.f)) + legCov.alpha();
    phi = std::fmod(phi, o2::constants::math::TwoPI);
    if (phi < 0.f) {
      phi += o2::constants::math::TwoPI;
    }
    return {pt, phi};
  }

  // MC history of a leg
  template <typename TMCParticle, typename TMCParticles>
  static void fillLegMC(LegCopy& c, TMCParticle const& mcLeg, TMCParticles const& mcParticles)
  {
    c.mcParticleId = mcLeg.globalIndex();
    c.mcPdg = mcLeg.pdgCode();
    c.mcVx = mcLeg.vx();
    c.mcVy = mcLeg.vy();
    c.mcVz = mcLeg.vz();
    if (!mcLeg.has_mothers()) {
      return;
    }
    const auto mother = mcParticles.rawIteratorAt(mcLeg.mothersIds()[0]);
    c.mcMotherId = mother.globalIndex();
    c.mcMotherPdg = mother.pdgCode();
    c.mcMotherPx = mother.px();
    c.mcMotherPy = mother.py();
    c.mcMotherPz = mother.pz();
    if (mother.has_mothers()) { // e.g. the e+- that radiated a bremsstrahlung photon
      const auto grandMother = mcParticles.rawIteratorAt(mother.mothersIds()[0]);
      c.mcGrandMotherId = grandMother.globalIndex();
      c.mcGrandMotherPdg = grandMother.pdgCode();
      if (grandMother.has_mothers()) { // e.g. the eta of eta -> 3 pi0
        const auto ggMother = mcParticles.rawIteratorAt(grandMother.mothersIds()[0]);
        c.mcGGMotherId = ggMother.globalIndex();
        c.mcGGMotherPdg = ggMother.pdgCode();
      }
    }
    auto root = mother;
    int rootDepth = 0;
    for (; root.has_mothers() && rootDepth < kMaxAncestorDepth; ++rootDepth) {
      root = mcParticles.rawIteratorAt(root.mothersIds()[0]);
    }
    c.mcRootId = static_cast<int>(root.globalIndex());
    c.mcRootPdg = root.pdgCode();
    c.mcRootDepth = rootDepth;
    auto anc = mother;
    for (int depth = 0; !anc.isPhysicalPrimary() && anc.has_mothers() && depth < kMaxAncestorDepth; ++depth) {
      anc = mcParticles.rawIteratorAt(anc.mothersIds()[0]);
    }
    c.mcAncestorId = anc.isPhysicalPrimary() ? static_cast<int>(anc.globalIndex()) : -1;
  }

  template <bool isMC, bool isDerived, typename TLeg, typename TLegCov, typename TTracks, typename TMCParticles>
  LegCopy copyLeg(TLeg const& leg, TLegCov const& legCov, [[maybe_unused]] TTracks const* tracks, [[maybe_unused]] TMCParticles const* mcParticles)
  {
    LegCopy c;
    c.sign = leg.sign();
    c.x = legCov.x();
    c.alpha = legCov.alpha();
    c.y = legCov.y();
    c.z = legCov.z();
    c.snp = legCov.snp();
    c.tgl = legCov.tgl();
    c.signed1Pt = legCov.signed1Pt();
    const auto& cov = legCov.cov();
    for (size_t i = 0; i < c.cov.size(); ++i) {
      c.cov[i] = cov[i];
    }
    c.tpcNClsFindable = leg.tpcNClsFindable();
    c.tpcNClsFindableMinusFound = leg.tpcNClsFindableMinusFound();
    c.tpcNClsFindableMinusCrossedRows = leg.tpcNClsFindableMinusCrossedRows();
    c.tpcNClsShared = leg.tpcNClsShared();
    c.tpcChi2NCl = leg.tpcChi2NCl();
    c.tpcInnerParam = leg.tpcInnerParam();
    c.tpcSignal = leg.tpcSignal();
    c.tpcNSigmaEl = leg.tpcNSigmaEl();
    c.tpcNSigmaPi = leg.tpcNSigmaPi();
    c.itsClusterSizes = leg.itsClusterSizes();
    c.itsChi2NCl = leg.itsChi2NCl();
    c.detectorMap = leg.detectorMap();
    c.origTrackId = leg.trackId();
    const auto [ptCov, phiCov] = ptPhiOf(legCov);
    if constexpr (isDerived) {
      c.pt = ptCov;
      c.phi = phiCov;
      if constexpr (isMC) {
        if (leg.has_emmcparticle()) {
          fillLegMC(c, leg.template emmcparticle_as<TMCParticles>(), *mcParticles);
        }
      }
    } else {
      const auto track = tracks->rawIteratorAt(leg.trackId());
      c.pt = track.pt(); // as the builder reads them from the AO2D track
      c.phi = track.phi();
      registry.fill(HIST("Leg/hPtPhiFromCov"), c.pt > 0.f ? (c.pt - ptCov) / c.pt : 0.f, RecoDecay::constrainAngle(c.phi - phiCov, -o2::constants::math::PI));
      if constexpr (isMC) {
        if (track.has_mcParticle()) {
          fillLegMC(c, track.template mcParticle_as<TMCParticles>(), *mcParticles);
        }
      }
    }
    return c;
  }

  void writeLeg(LegCopy const& l, PhotonCopy const& p)
  {
    mixLegs(mixCollisions.lastIndex(), l.sign,
            l.x, l.alpha, l.y, l.z, l.snp, l.tgl, l.signed1Pt, l.pt, l.phi,
            l.cov[0], l.cov[1], l.cov[2], l.cov[3], l.cov[4], l.cov[5], l.cov[6], l.cov[7], l.cov[8], l.cov[9], l.cov[10], l.cov[11], l.cov[12], l.cov[13], l.cov[14],
            l.tpcNClsFindable, l.tpcNClsFindableMinusFound, l.tpcNClsFindableMinusCrossedRows, l.tpcNClsShared,
            l.tpcChi2NCl, l.tpcInnerParam, l.tpcSignal, l.tpcNSigmaEl, l.tpcNSigmaPi,
            l.itsClusterSizes, l.itsChi2NCl, l.detectorMap,
            p.eventUid, p.photonUid, p.v0Id, p.px, p.py, p.pz, p.vx, p.vy, p.vz,
            p.pvX, p.pvY, p.pvZ, p.pvCovXX, p.pvCovXY, p.pvCovXZ, p.pvCovYY, p.pvCovYZ, p.pvCovZZ, p.pvChi2, p.pvNContrib,
            l.origTrackId, l.mcParticleId, l.mcMotherId, l.mcPdg, l.mcMotherPdg, l.mcMotherPx, l.mcMotherPy, l.mcMotherPz, l.mcGrandMotherId, l.mcGrandMotherPdg, l.mcAncestorId,
            l.mcGGMotherId, l.mcGGMotherPdg, l.mcRootId, l.mcRootPdg, l.mcRootDepth, l.mcVx, l.mcVy, l.mcVz);
  }

  static PhotonCopy translated(PhotonCopy const& p, const float dx, const float dy, const float dz)
  {
    PhotonCopy t = p;
    t.vx += dx;
    t.vy += dy;
    t.vz += dz;
    t.pvX += dx; // its own vertex moves with it
    t.pvY += dy;
    t.pvZ += dz;
    for (auto& l : t.legs) {
      const float ca = std::cos(l.alpha), sa = std::sin(l.alpha);
      l.x += dx * ca + dy * sa;
      l.y += -dx * sa + dy * ca;
      l.z += dz;
      l.mcVx += dx; // the true conversion point moves with the photon
      l.mcVy += dy;
      l.mcVz += dz;
    }
    return t;
  }

  void writeArtificialEvent(VertexCopy const& v, EventCopy const& current, const bool isMixed, const int zBin, const int occBin,
                            std::vector<PhotonCopy const*> const& photons, const int currentCollisionId)
  {
    const float type = isMixed ? 1.f : 0.f;
    const int firstLeg = mixLegs.lastIndex() + 1;
    mixCollisions(v.x, v.y, v.z, v.covXX, v.covXY, v.covXZ, v.covYY, v.covYZ, v.covZZ, v.chi2, v.numContrib,
                  v.centFT0M, v.centFT0A, v.centFT0C, isMixed, currentCollisionId, current.eventUid,
                  zBin, occBin, current.occupancy, firstLeg, nLegsOf(photons), v.runNumber, v.multNTracksPVeta1, v.multFT0M, v.timestamp);
    for (const auto* p : photons) {
      if (p->legs.size() > 1) { // not a free-leg stub
        registry.fill(HIST("Mix/hPhotonPtEta"), type, std::hypot(p->px, p->py), std::asinh(p->pz / std::hypot(p->px, p->py)));
      }
      for (const auto& l : p->legs) {
        const bool hasITS = (l.detectorMap & o2::aod::track::ITS) != 0;
        const bool hasTPC = (l.detectorMap & o2::aod::track::TPC) != 0;
        registry.fill(HIST("Mix/hLegTypes"), type, hasITS && hasTPC ? 0.f : (hasTPC ? 1.f : (hasITS ? 2.f : 3.f)));
        writeLeg(l, *p);
      }
    }
    registry.fill(HIST("Mix/hLegsMinus2N"), type, (mixLegs.lastIndex() + 1 - firstLeg) - nLegsOf(photons));
  }

  static int nLegsOf(std::vector<PhotonCopy const*> const& photons)
  {
    int n = 0;
    for (const auto* p : photons) {
      n += static_cast<int>(p->legs.size());
    }
    return n;
  }

  std::vector<PhotonCopy> toGroups(std::vector<PhotonCopy> const& cands)
  {
    const int n = cands.size();
    std::vector<int> parent(n);
    std::iota(parent.begin(), parent.end(), 0);
    std::function<int(int)> find = [&](int i) { return parent[i] == i ? i : parent[i] = find(parent[i]); };
    std::unordered_map<int, int> firstOfTrack;
    for (int i = 0; i < n; ++i) {
      for (const auto& l : cands[i].legs) {
        const auto [it, inserted] = firstOfTrack.emplace(l.origTrackId, i);
        if (!inserted) {
          parent[find(i)] = find(it->second);
        }
      }
    }
    std::map<int, PhotonCopy> groups; // root -> group
    std::map<int, std::set<int>> tracksOfGroup;
    for (int i = 0; i < n; ++i) {
      const int r = find(i);
      auto [it, inserted] = groups.try_emplace(r, cands[i]);
      if (inserted) {
        it->second.legs.clear();
        it->second.photonUid = mPhotonUid++;
      }
      for (const auto& l : cands[i].legs) {
        if (tracksOfGroup[r].insert(l.origTrackId).second) {
          it->second.legs.push_back(l);
        }
      }
    }
    std::vector<PhotonCopy> out;
    for (auto& [r, g] : groups) {
      registry.fill(HIST("Mix/hGroupLegs"), static_cast<float>(g.legs.size()));
      out.push_back(std::move(g));
    }
    return out;
  }

  template <bool isMC, bool isDerived, typename TCollisions, typename TPhotons, typename TLegs, typename TTracks, typename TMCParticles>
  void runMixer(TCollisions const& collisions, TPhotons const& photons, TLegs const& legs, aod::V0LegsCov const& legsCov,
                aod::PCMFreeLegs const* freeLegs, TTracks const* tracks, TMCParticles const* mcParticles)
  {
    if (legs.size() != 2 * photons.size() || legsCov.size() != legs.size()) {
      LOG(fatal) << "V0LegsCov (" << legsCov.size() << ") / V0Legs (" << legs.size() << ") / V0PhotonsKF (" << photons.size() << ") do not match: run builder run 1 with fillLegCov = true";
    }
    if (cfgPoolFreeLegs && (freeLegs == nullptr || (freeLegs->size() == 0 && photons.size() > 0))) {
      LOG(fatal) << "cfgPoolFreeLegs: PCMFreeLegs is empty, run builder run 1 with fillFreeLegs = true";
    }
    std::unordered_map<int, std::vector<int64_t>> photonsOfCollision, freeLegsOfCollision;
    for (const auto& photon : photons) {
      if constexpr (isDerived) {
        photonsOfCollision[photon.pmeventId()].push_back(photon.globalIndex());
      } else {
        photonsOfCollision[photon.collisionId()].push_back(photon.globalIndex());
      }
    }
    if (cfgPoolFreeLegs) {
      for (const auto& leg : *freeLegs) {
        freeLegsOfCollision[leg.collisionId()].push_back(leg.globalIndex());
      }
    }

    for (const auto& collision : collisions) {
      registry.fill(HIST("Event/hCutFlow"), kAll);
      int runNumber = 0;
      uint64_t timestamp = 0;
      if constexpr (isDerived) {
        runNumber = collision.runNumber();
        timestamp = collision.timestamp();
      } else {
        const auto bc = collision.template foundBC_as<aod::BCsWithTimestamps>();
        runNumber = bc.runNumber();
        timestamp = bc.timestamp();
      }
      if (runNumber != mRunNumber) {
        if (mRunNumber >= 0) {
          registry.fill(HIST("Pool/hResets"), 0.f);
        }
        mPools.clear();
        mNPooledPhotons = 0;
        mRunNumber = runNumber;
      }
      if constexpr (!isDerived) {
        if (!collision.isSelected()) {
          continue;
        }
      }
      registry.fill(HIST("Event/hCutFlow"), kSelected);
      const float cent = centrality(collision);
      if (cent < cfgCentMin || cfgCentMax <= cent) {
        continue;
      }
      registry.fill(HIST("Event/hCutFlow"), kCent);
      if (std::fabs(collision.posZ()) >= cfgZMax) {
        continue;
      }
      registry.fill(HIST("Event/hCutFlow"), kZ);
      const float occ = occupancy(collision);
      const int occBin = occupancyBin(occ);
      if (occBin < 0) {
        continue;
      }
      registry.fill(HIST("Event/hCutFlow"), kOccupancy);
      const auto it = photonsOfCollision.find(collision.globalIndex());
      if (it == photonsOfCollision.end()) {
        continue;
      }
      registry.fill(HIST("Event/hCutFlow"), kHasPhotons);
      const int zBin = static_cast<int>(std::floor((collision.posZ() + cfgZMax) / cfgZBinWidth));
      registry.fill(HIST("Event/hZ"), collision.posZ());
      registry.fill(HIST("Event/hOccupancy"), occ);

      // copies of this event's photons
      EventCopy current;
      current.eventUid = mEventUid++;
      current.x = collision.posX();
      current.y = collision.posY();
      current.z = collision.posZ();
      current.occupancy = occ;
      current.cent = cent;
      for (const auto& photonId : it->second) {
        const auto photon = photons.rawIteratorAt(photonId);
        PhotonCopy p;
        p.eventUid = current.eventUid;
        p.photonUid = mPhotonUid++;
        p.v0Id = photon.v0Id();
        p.px = photon.px();
        p.py = photon.py();
        p.pz = photon.pz();
        p.vx = photon.vx();
        p.vy = photon.vy();
        p.vz = photon.vz();
        p.pvX = collision.posX();
        p.pvY = collision.posY();
        p.pvZ = collision.posZ();
        p.pvCovXX = collision.covXX();
        p.pvCovXY = collision.covXY();
        p.pvCovXZ = collision.covXZ();
        p.pvCovYY = collision.covYY();
        p.pvCovYZ = collision.covYZ();
        p.pvCovZZ = collision.covZZ();
        p.pvChi2 = collision.chi2();
        p.pvNContrib = collision.numContrib();
        p.legs.push_back(copyLeg<isMC, isDerived>(legs.rawIteratorAt(photon.posTrackId()), legsCov.rawIteratorAt(photon.posTrackId()), tracks, mcParticles));
        p.legs.push_back(copyLeg<isMC, isDerived>(legs.rawIteratorAt(photon.negTrackId()), legsCov.rawIteratorAt(photon.negTrackId()), tracks, mcParticles));
        current.photons.push_back(p);
      }
      if (cfgPoolGroups) {
        registry.fill(HIST("Event/hNCandidates"), static_cast<float>(current.photons.size()));
        current.photons = toGroups(current.photons);
      }
      if constexpr (!isDerived) {
        if (cfgPoolFreeLegs) { // one-leg stubs: event and vertex of this event, no photon momentum, v0Id -2
          if (const auto itF = freeLegsOfCollision.find(collision.globalIndex()); itF != freeLegsOfCollision.end()) {
            for (const auto& legId : itF->second) {
              const auto leg = freeLegs->rawIteratorAt(legId);
              PhotonCopy f;
              f.eventUid = current.eventUid;
              f.photonUid = mPhotonUid++;
              f.v0Id = -2;
              f.pvX = collision.posX();
              f.pvY = collision.posY();
              f.pvZ = collision.posZ();
              f.pvCovXX = collision.covXX();
              f.pvCovXY = collision.covXY();
              f.pvCovXZ = collision.covXZ();
              f.pvCovYY = collision.covYY();
              f.pvCovYZ = collision.covYZ();
              f.pvCovZZ = collision.covZZ();
              f.pvChi2 = collision.chi2();
              f.pvNContrib = collision.numContrib();
              f.legs.push_back(copyLeg<isMC, false>(leg, leg, tracks, mcParticles));
              current.freeLegs.push_back(f);
            }
          }
          registry.fill(HIST("Event/hNFreeLegs"), static_cast<float>(current.freeLegs.size()));
        }
      }
      const int nPhotons = static_cast<int>(current.photons.size());
      registry.fill(HIST("Event/hNPhotons"), nPhotons);

      auto& pool = mPools[{zBin, occBin}];
      const int depth = static_cast<int>(pool.size());
      registry.fill(HIST("Pool/hPartnersVsN"), nPhotons, depth);
      registry.fill(HIST("Pool/hDepthVsZ"), collision.posZ(), depth);
      registry.fill(HIST("Pool/hDepthVsOccBin"), occBin, depth);

      VertexCopy v{collision.posX(), collision.posY(), collision.posZ(), collision.covXX(), collision.covXY(), collision.covXZ(),
                   collision.covYY(), collision.covYZ(), collision.covZZ(), collision.chi2(), collision.numContrib(),
                   collision.centFT0M(), collision.centFT0A(), collision.centFT0C(), runNumber, collision.multNTracksPVeta1(), collision.multFT0M(), timestamp};
      const int nPartnersNeeded = cfgOverlay > 0 ? cfgOverlay.value : (cfgInject > 0 ? cfgInject.value : (cfgInject == -2 ? std::max(1, (nPhotons + 1) / 2) : nPhotons)); // overlay / injection: K partner events
      const bool poolFull = depth >= nPartnersNeeded;
      if (!poolFull) {
        registry.fill(HIST("Event/hCutFlow"), kPoolTooSmall);
      }
      if (poolFull || !cfgRequireFullPool) {
        registry.fill(HIST("Event/hCutFlow"), kWritten);
      }
      std::vector<PhotonCopy const*> host; // the complete current event: its photons and (cfgPoolFreeLegs) its free legs
      host.reserve(current.photons.size() + current.freeLegs.size());
      for (const auto& p : current.photons) {
        host.push_back(&p);
      }
      for (const auto& f : current.freeLegs) {
        host.push_back(&f);
      }
      if (poolFull || !cfgRequireFullPool) {
        writeArtificialEvent(v, current, false, zBin, occBin, host, collision.globalIndex());
      }
      if (poolFull) {
        std::vector<int> order(depth);
        auto drawOrder = [&]() { // partner events in random order (cfgPhotonWeighted: weighted sampling without replacement, largest log(u) / w first)
          std::iota(order.begin(), order.end(), 0);
          if (cfgPhotonWeighted) {
            std::uniform_real_distribution<double> uni(0., 1.);
            std::vector<double> key(depth);
            for (int k = 0; k < depth; ++k) {
              key[k] = std::log(std::max(uni(mRandom), 1e-300)) / static_cast<double>(pool[k].photons.size()); // o2-linter: disable=magic-number (log floor)
            }
            std::sort(order.begin(), order.end(), [&](int x, int y) { return key[x] > key[y]; });
          } else {
            std::shuffle(order.begin(), order.end(), mRandom);
          }
        };
        auto fillPartnerQA = [&](EventCopy const& partner, int age) {
          registry.fill(HIST("Mix/hDz"), partner.z - current.z);
          registry.fill(HIST("Mix/hDOccupancy"), partner.occupancy - current.occupancy);
          registry.fill(HIST("Mix/hDCent"), partner.cent - current.cent);
          registry.fill(HIST("Mix/hPartnerAge"), age);
          registry.fill(HIST("Mix/hSourceIsCurrent"), partner.eventUid == current.eventUid ? 1.f : 0.f);
        };
        if (cfgOverlay > 0) { // overlay: K artificial events = the complete current event + the complete partner event k
          for (int iMix = 0; iMix < cfgNMixPerEvent; ++iMix) {
            drawOrder();
            for (int k = 0; k < std::min(nPartnersNeeded, depth); ++k) {
              auto& partner = pool[order[k]];
              fillPartnerQA(partner, order[k]);
              std::vector<PhotonCopy const*> overlay(host);
              std::vector<std::unique_ptr<PhotonCopy>> aligned; // translated copies; reserved, so the pointers stay valid
              aligned.reserve(partner.photons.size() + partner.freeLegs.size());
              for (auto* src : {&partner.photons, &partner.freeLegs}) {
                for (auto& p : *src) {
                  if (src == &partner.photons) {
                    p.timesUsed++;
                  }
                  if (cfgAlignToVertex) {
                    aligned.push_back(std::make_unique<PhotonCopy>(translated(p, current.x - partner.x, current.y - partner.y, current.z - partner.z)));
                    overlay.push_back(aligned.back().get());
                  } else {
                    overlay.push_back(&p);
                  }
                }
              }
              registry.fill(HIST("Mix/hDistinctSourcesMinusN"), 2 - nPhotons); // overlay: the current event + one partner
              registry.fill(HIST("Event/hNPhotonsMixedVsCurrent"), nPhotons, nPhotons + static_cast<int>(partner.photons.size()));
              writeArtificialEvent(v, current, true, zBin, occBin, overlay, collision.globalIndex());
            }
          }
        }
        if (cfgInject != 0) {
          for (int iMix = 0; iMix < cfgNMixPerEvent; ++iMix) {
            drawOrder();
            for (int k = 0; k < std::min(nPartnersNeeded, depth); ++k) {
              auto& partner = pool[order[k]];
              std::uniform_int_distribution<size_t> pick(0, partner.photons.size() - 1);
              auto& p = partner.photons[pick(mRandom)];
              p.timesUsed++;
              fillPartnerQA(partner, order[k]);
              std::vector<PhotonCopy const*> injected(host);
              PhotonCopy alignedCopy;
              if (cfgAlignToVertex) {
                alignedCopy = translated(p, current.x - partner.x, current.y - partner.y, current.z - partner.z);
                injected.push_back(&alignedCopy);
              } else {
                injected.push_back(&p);
              }
              registry.fill(HIST("Mix/hDistinctSourcesMinusN"), 2 - nPhotons);
              registry.fill(HIST("Event/hNPhotonsMixedVsCurrent"), nPhotons, nPhotons + 1);
              writeArtificialEvent(v, current, true, zBin, occBin, injected, collision.globalIndex());
            }
          }
        }
        for (int iMix = 0; iMix < ((cfgInject != 0 || cfgOverlay > 0) ? 0 : cfgNMixPerEvent.value); ++iMix) {
          drawOrder();
          std::vector<PhotonCopy const*> mixedPhotons;
          std::vector<std::unique_ptr<PhotonCopy>> aligned; // translated copies (cfgAlignToVertex); reserved, so the pointers stay valid
          aligned.reserve(nPhotons + (cfgPoolFreeLegs ? static_cast<int>(current.freeLegs.size()) : 0));
          std::set<int64_t> sources;
          for (int k = 0; k < nPhotons; ++k) { // N partner events, one random photon each
            auto& partner = pool[order[k]];
            std::uniform_int_distribution<size_t> pick(0, partner.photons.size() - 1);
            auto& p = partner.photons[pick(mRandom)];
            p.timesUsed++;
            if (cfgAlignToVertex) {
              aligned.push_back(std::make_unique<PhotonCopy>(translated(p, current.x - partner.x, current.y - partner.y, current.z - partner.z)));
              // cppcheck-suppress invalidContainer
              mixedPhotons.push_back(aligned.back().get());
            } else {
              mixedPhotons.push_back(&p);
            }
            sources.insert(partner.eventUid);
            fillPartnerQA(partner, order[k]);
          }
          registry.fill(HIST("Mix/hDistinctSourcesMinusN"), static_cast<int>(sources.size()) - nPhotons);
          registry.fill(HIST("Event/hNPhotonsMixedVsCurrent"), nPhotons, mixedPhotons.size());
          if (cfgPoolFreeLegs) { // M free legs of the current event -> M drawn without replacement from the free legs of the N partner events
            std::vector<PhotonCopy const*> pooledFree;
            for (int k = 0; k < nPhotons; ++k) {
              for (const auto& f : pool[order[k]].freeLegs) {
                pooledFree.push_back(&f);
              }
            }
            std::shuffle(pooledFree.begin(), pooledFree.end(), mRandom);
            const size_t nTake = std::min(pooledFree.size(), current.freeLegs.size());
            for (size_t k = 0; k < nTake; ++k) {
              if (cfgAlignToVertex) {
                const auto& partner = *pooledFree[k];
                aligned.push_back(std::make_unique<PhotonCopy>(translated(partner, current.x - partner.pvX, current.y - partner.pvY, current.z - partner.pvZ)));
                // cppcheck-suppress invalidContainer
                mixedPhotons.push_back(aligned.back().get());
              } else {
                mixedPhotons.push_back(pooledFree[k]);
              }
            }
            registry.fill(HIST("Mix/hFreeLegsMixedVsCurrent"), static_cast<float>(current.freeLegs.size()), static_cast<float>(nTake));
          }
          writeArtificialEvent(v, current, true, zBin, occBin, mixedPhotons, collision.globalIndex());
        }
      }

      // this event joins its pool after it was mixed (never mixed with itself)
      mNPooledPhotons += nPhotons;
      pool.push_front(std::move(current));
      if (static_cast<int>(pool.size()) > cfgPoolDepth) {
        for (const auto& p : pool.back().photons) {
          registry.fill(HIST("Pool/hTimesUsed"), p.timesUsed);
        }
        mNPooledPhotons -= static_cast<int64_t>(pool.back().photons.size());
        pool.pop_back();
      }
      registry.fill(HIST("Pool/hNPooledPhotons"), mNPooledPhotons);
      registry.fill(HIST("Pool/hNPools"), mPools.size());
    }
  }

  // AO2D input (builder run 1 in the same workflow)
  void processRec(MyCollisions const& collisions, aod::BCsWithTimestamps const&, aod::V0PhotonsKF const& photons, aod::V0Legs const& legs, aod::V0LegsCov const& legsCov, aod::PCMFreeLegs const& freeLegs, aod::TracksIU const& tracks)
  {
    runMixer<false, false>(collisions, photons, legs, legsCov, &freeLegs, &tracks, static_cast<aod::McParticles const*>(nullptr));
  }
  PROCESS_SWITCH(PhotonMixedEventProducer, processRec, "data, AO2D input", true);

  void processMC(MyCollisions const& collisions, aod::BCsWithTimestamps const&, aod::V0PhotonsKF const& photons, aod::V0Legs const& legs, aod::V0LegsCov const& legsCov, aod::PCMFreeLegs const& freeLegs, MyTracksMC const& tracks, aod::McParticles const& mcParticles)
  {
    runMixer<true, false>(collisions, photons, legs, legsCov, &freeLegs, &tracks, &mcParticles);
  }
  PROCESS_SWITCH(PhotonMixedEventProducer, processMC, "MC, AO2D input (leg truth into PCMMixLegs)", false);

  // derived data (EM skim with V0LegsCov and, from create-emevent-photon fillVertex, EMEventsXY + EMEventsCov)
  void processDerived(MyEvents const& collisions, MyPhotonsDerived const& photons, aod::V0Legs const& legs, aod::V0LegsCov const& legsCov)
  {
    runMixer<false, true>(collisions, photons, legs, legsCov, nullptr, static_cast<aod::TracksIU const*>(nullptr), static_cast<aod::EMMCParticles const*>(nullptr));
  }
  PROCESS_SWITCH(PhotonMixedEventProducer, processDerived, "data, derived-data input", false);

  void processDerivedMC(MyEvents const& collisions, MyPhotonsDerived const& photons, MyLegsMCDerived const& legs, aod::V0LegsCov const& legsCov, aod::EMMCParticles const& mcParticles)
  {
    runMixer<true, true>(collisions, photons, legs, legsCov, nullptr, static_cast<aod::TracksIU const*>(nullptr), &mcParticles);
  }
  PROCESS_SWITCH(PhotonMixedEventProducer, processDerivedMC, "MC, derived-data input (leg truth from V0LegMCLabels / EMMCParticles)", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<PhotonMixedEventProducer>(cfgc, TaskName{"photon-mixed-event-producer"})};
}
