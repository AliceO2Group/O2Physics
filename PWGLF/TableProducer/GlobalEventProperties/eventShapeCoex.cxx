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
/// \brief Per-collision forward/backward sub-event quantities of the mid-rapidity tracks.
/// \author Cristian Andrei <Cristian.Andrei@cern.ch>
///
/// First O2Physics task of the AliPhysics PWGLF/SPECTRA/MultEvShape analyses.

#include "PWGLF/DataModel/LFEventTopologyTables.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/runDataProcessing.h>

#include <array>
#include <cmath>
#include <cstdint>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace
{
/// Sums over one sub-event, shared by the reconstructed and generated levels.
struct SubEventSums {
  uint32_t n = 0;
  double sumPt = 0.;
  double sumPt2 = 0.;
  double qx2 = 0., qy2 = 0.;
  double qx4 = 0., qy4 = 0.;

  void add(double pt, float phi)
  {
    // phi stays float, as stored in the AOD
    const double c = std::cos(phi);
    const double s = std::sin(phi);
    const double cos2 = c * c - s * s;
    const double sin2 = 2. * s * c;
    ++n;
    sumPt += pt;
    sumPt2 += pt * pt;
    qx2 += cos2;
    qy2 += sin2;
    qx4 += cos2 * cos2 - sin2 * sin2;
    qy4 += 2. * sin2 * cos2;
  }

  [[nodiscard]] float meanPt() const { return (n > 0) ? static_cast<float>(sumPt / n) : std::nanf(""); }
};
} // namespace

struct EventShapeCoex {
  Produces<aod::EvShapeCoex> coex;
  Produces<aod::EvShapeCoexGen> coexGen;
  Produces<aod::EvShapeCoexMcLabels> coexMcLabels;

  Service<o2::framework::O2DatabasePDG> pdg{};

  Configurable<float> vtxZCut{"vtxZCut", 10.0f, "Max |PV z| (cm)"};
  Configurable<float> ptMin{"ptMin", 0.15f, "Minimum track pT (GeV/c)"};
  Configurable<float> ptMax{"ptMax", 2.0f, "Maximum track pT (GeV/c)"};
  Configurable<float> etaMax{"etaMax", 0.8f, "Max |eta| (mid-rapidity measurement region)"};
  Configurable<float> gapHalf{"gapHalf", 0.2f, "Central-gap half-width: |eta| <= gapHalf excluded from F/B"};
  Configurable<float> etaFwdAMin{"etaFwdAMin", 3.5f, "FT0-A proxy: lower eta edge (generator level)"};
  Configurable<float> etaFwdAMax{"etaFwdAMax", 4.9f, "FT0-A proxy: upper eta edge (generator level)"};
  Configurable<float> etaFwdCMin{"etaFwdCMin", -3.3f, "FT0-C proxy: lower eta edge (generator level)"};
  Configurable<float> etaFwdCMax{"etaFwdCMax", -2.1f, "FT0-C proxy: upper eta edge (generator level)"};

  Filter trackFilter = (aod::track::pt > ptMin) && (aod::track::pt < ptMax) && (nabs(aod::track::eta) < etaMax) && requireGlobalTrackInFilter();

  using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults>;
  using MyCollisionsMc = soa::Join<MyCollisions, aod::McCollisionLabels>;
  using MyTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection>>;

  void init(InitContext const&)
  {
    if (doprocessReco && doprocessRecoMC) {
      LOGF(fatal, "processReco and processRecoMC both fill EvShapeCoex; enable only one of them.");
    }
    if (doprocessRecoMC && !doprocessMC) {
      LOGF(fatal, "processRecoMC links to EvShapeCoexGen rows and requires processMC in the same job.");
    }
  }

  /// Fills one EvShapeCoex row; returns false if the collision is rejected.
  template <typename TCollision>
  bool fillReco(TCollision const& coll, MyTracks const& tracks)
  {
    if (!coll.sel8() || std::abs(coll.posZ()) > vtxZCut) {
      return false;
    }

    const float g = gapHalf;
    SubEventSums fwd, bwd;
    uint32_t nGap = 0, nPlus = 0, nMinus = 0, nPVC = 0;
    for (const auto& track : tracks) {
      if (track.isPVContributor()) {
        ++nPVC;
      }
      const float eta = track.eta();
      if (eta > g) {
        fwd.add(track.pt(), track.phi());
      } else if (eta < -g) {
        bwd.add(track.pt(), track.phi());
      } else {
        ++nGap;
        continue;
      }
      if (track.sign() > 0) {
        ++nPlus;
      } else if (track.sign() < 0) {
        ++nMinus;
      }
    }

    // bit i of qaBits = QaBitOrder[i]
    static constexpr int NQaBits = 12;
    static constexpr std::array<o2::aod::evsel::EventSelectionFlags, NQaBits> QaBitOrder = {
      o2::aod::evsel::kNoSameBunchPileup, o2::aod::evsel::kIsGoodZvtxFT0vsPV,
      o2::aod::evsel::kIsVertexITSTPC, o2::aod::evsel::kIsVertexTOFmatched,
      o2::aod::evsel::kNoCollInTimeRangeNarrow, o2::aod::evsel::kNoCollInTimeRangeStrict,
      o2::aod::evsel::kNoCollInTimeRangeStandard, o2::aod::evsel::kNoCollInRofStrict,
      o2::aod::evsel::kNoCollInRofStandard, o2::aod::evsel::kNoHighMultCollInPrevRof,
      o2::aod::evsel::kNoITSROFrameBorder, o2::aod::evsel::kNoTimeFrameBorder};
    uint16_t qaBits = 0;
    for (int i = 0; i < NQaBits; ++i) {
      if (coll.selection_bit(QaBitOrder[i])) {
        qaBits |= static_cast<uint16_t>(1u << i);
      }
    }

    const auto& bc = coll.template bc_as<aod::BCs>();

    coex(static_cast<int>(coll.globalIndex()), bc.runNumber(), bc.globalBC(), coll.posZ(),
         coll.multFT0A(), coll.multFT0C(),
         static_cast<uint16_t>(fwd.n), static_cast<uint16_t>(bwd.n), static_cast<uint16_t>(nGap),
         fwd.meanPt(), bwd.meanPt(), static_cast<float>(fwd.sumPt2), static_cast<float>(bwd.sumPt2),
         static_cast<float>(fwd.qx2), static_cast<float>(fwd.qy2),
         static_cast<float>(bwd.qx2), static_cast<float>(bwd.qy2),
         static_cast<float>(fwd.qx4), static_cast<float>(fwd.qy4),
         static_cast<float>(bwd.qx4), static_cast<float>(bwd.qy4),
         static_cast<uint16_t>(nPlus), static_cast<uint16_t>(nMinus),
         qaBits, coll.trackOccupancyInTimeRange(), coll.ft0cOccupancyInTimeRange(),
         coll.numContrib(), coll.collisionTimeRes(), static_cast<uint16_t>(nPVC));
    return true;
  }

  void processReco(MyCollisions::iterator const& coll, MyTracks const& tracks, aod::BCs const&)
  {
    fillReco(coll, tracks);
  }
  PROCESS_SWITCH(EventShapeCoex, processReco, "Reconstructed-level reduction (EvShapeCoex)", true);

  void processRecoMC(MyCollisionsMc::iterator const& coll, MyTracks const& tracks, aod::BCs const&)
  {
    if (fillReco(coll, tracks)) {
      // EvShapeCoexGen has one row per McCollision, in order: the McCollision index is its row index
      coexMcLabels(coll.mcCollisionId());
    }
  }
  PROCESS_SWITCH(EventShapeCoex, processRecoMC, "Reconstructed-level reduction with MC labels (EvShapeCoex + EvShapeCoexMcLabels)", false);

  void processMC(aod::McCollision const& mcCollision, aod::McParticles const& particles)
  {
    static constexpr double ChargeTolerance = 1.e-3;

    const float g = gapHalf;
    const float em = etaMax;
    SubEventSums fwd, bwd;
    uint32_t nFwdA = 0, nFwdC = 0, nGap = 0, nPlus = 0, nMinus = 0;
    for (const auto& particle : particles) {
      if (!particle.isPhysicalPrimary()) {
        continue;
      }
      const auto* pdgParticle = pdg->GetParticle(particle.pdgCode());
      if (pdgParticle == nullptr) {
        continue; // unknown PDG code, charge undetermined
      }
      const double charge = pdgParticle->Charge();
      if (std::abs(charge) < ChargeTolerance) {
        continue;
      }
      const float eta = particle.eta();
      // forward counts: no pT window
      if (eta > etaFwdAMin && eta < etaFwdAMax) {
        ++nFwdA;
      }
      if (eta > etaFwdCMin && eta < etaFwdCMax) {
        ++nFwdC;
      }
      const float pt = particle.pt();
      if (pt <= ptMin || pt >= ptMax || std::abs(eta) >= em) {
        continue;
      }
      if (eta > g) {
        fwd.add(pt, particle.phi());
      } else if (eta < -g) {
        bwd.add(pt, particle.phi());
      } else {
        ++nGap;
        continue;
      }
      if (charge > 0.) {
        ++nPlus;
      } else {
        ++nMinus;
      }
    }

    coexGen(static_cast<int>(mcCollision.globalIndex()), mcCollision.posZ(),
            static_cast<uint16_t>(nFwdA), static_cast<uint16_t>(nFwdC),
            static_cast<uint16_t>(fwd.n), static_cast<uint16_t>(bwd.n), static_cast<uint16_t>(nGap),
            fwd.meanPt(), bwd.meanPt(), static_cast<float>(fwd.sumPt2), static_cast<float>(bwd.sumPt2),
            static_cast<float>(fwd.qx2), static_cast<float>(fwd.qy2),
            static_cast<float>(bwd.qx2), static_cast<float>(bwd.qy2),
            static_cast<float>(fwd.qx4), static_cast<float>(fwd.qy4),
            static_cast<float>(bwd.qx4), static_cast<float>(bwd.qy4),
            static_cast<uint16_t>(nPlus), static_cast<uint16_t>(nMinus));
  }
  PROCESS_SWITCH(EventShapeCoex, processMC, "Generator-level reduction (EvShapeCoexGen)", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<EventShapeCoex>(cfgc)};
}
