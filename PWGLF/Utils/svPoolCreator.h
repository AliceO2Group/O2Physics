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

/// \file svPoolCreator.h
/// \brief Time-compatible secondary-vertex track-pair and triplet pools
/// \author ALICE Collaboration

#ifndef PWGLF_UTILS_SVPOOLCREATOR_H_
#define PWGLF_UTILS_SVPOOLCREATOR_H_

#include <CommonConstants/LHCConstants.h>
#include <DCAFitter/DCAFitterN.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/DataTypes.h>
#include <Framework/Logger.h>
#include <MathUtils/Primitive2D.h>

#include <Rtypes.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <unordered_map>
#include <utility>
#include <vector>

using CollBracket = o2::math_utils::Bracket<int>;

constexpr uint64_t bOffsetMax = 241; // track compatibility can never go beyond 6 mus (ITS)
constexpr int NChargeSigns = 2;

struct TrackCand {
  int Idxtr;
  CollBracket collBracket{};
};

struct SVCand {
  int tr0Idx;
  int tr1Idx;
  CollBracket collBracket{};
};

struct SVCand3 {
  int tr0Idx;
  int tr1Idx;
  int tr2Idx;
  CollBracket collBracket{};
};

class svPoolCreator
{
 public:
  svPoolCreator() = default;
  svPoolCreator(int track0Pdg, int track1Pdg)
    : track0Pdg(track0Pdg), track1Pdg(track1Pdg)
  {
  }

  void setPDGs(int track0Pdg, int track1Pdg)
  {
    this->track0Pdg = track0Pdg;
    this->track1Pdg = track1Pdg;
  }

  void clearPools()
  {
    for (auto& pool : trackCandPool) { // o2-linter: disable=const-ref-in-for-loop (The pool is cleared in place.)
      pool.clear();
    }
    tmap.clear();
    svCandPool.clear();
    bc2Coll.clear();
  }

  void setTimeMargin(float timeMargin) { timeMarginNS = timeMargin; }
  void setFitter(const o2::vertexing::DCAFitterN<2>& fitter) { this->fitter = fitter; }
  void setSkipAmbiTracks() { skipAmbiTracks = true; }
  o2::vertexing::DCAFitterN<2>* getFitter() { return &fitter; }
  std::array<std::vector<TrackCand>, 4> getTrackCandPool() { return trackCandPool; }

  template <typename C, typename BC>
  void fillBC2Coll(const C& collisions, BC const&)
  {
    for (unsigned i = 0; i < collisions.size(); i++) {
      auto collision = collisions.rawIteratorAt(i);
      if (!collision.has_bc()) {
        continue;
      }
      bc2Coll[collision.template bc_as<BC>().globalBC()] = i;
    }
  }

  template <typename T, typename C, typename AT, typename BC>
  void appendTrackCand(const T& trackCand, const C& collisions, int pdgHypo, const AT& ambiTracks, BC const&)
  {
    if (pdgHypo != track0Pdg && pdgHypo != track1Pdg) {
      LOGP(debug, "Wrong pdg hypothesis");
      return;
    }
    bool isDau0 = pdgHypo == track0Pdg;
    constexpr uint64_t BcInvalid = -1;
    uint64_t globalBC = BcInvalid;
    if (trackCand.has_collision()) {
      if (trackCand.template collision_as<C>().has_bc()) {
        globalBC = trackCand.template collision_as<C>().template bc_as<BC>().globalBC();
      }
    } else if (!skipAmbiTracks) {
      for (const auto& ambTrack : ambiTracks) {
        if (ambTrack.trackId() != trackCand.globalIndex()) {
          continue;
        }
        if (!ambTrack.has_bc() || ambTrack.template bc_as<BC>().size() == 0) {
          globalBC = BcInvalid;
          break;
        }
        globalBC = ambTrack.template bc_as<BC>().begin().globalBC();
        break;
      }
    } else {
      globalBC = BcInvalid;
    }

    if (globalBC == BcInvalid) {
      return;
    }

    uint64_t firstBC = globalBC < bOffsetMax ? 0 : globalBC - bOffsetMax;
    uint64_t lastBC = globalBC + bOffsetMax;
    int firstCollIdx = -1;
    while (firstBC < lastBC) {
      if (bc2Coll.find(firstBC) != bc2Coll.end()) {
        firstCollIdx = bc2Coll[firstBC];
        break;
      }
      firstBC++;
    }
    if (firstCollIdx == -1) {
      return;
    }

    // now loop over all the collisions to make the pool
    for (int collIdx = firstCollIdx; collIdx < collisions.size(); collIdx++) {
      const auto& collision = collisions.rawIteratorAt(collIdx);
      float collTime = collision.collisionTime();
      float collTimeRes2 = collision.collisionTimeRes() * collision.collisionTimeRes();
      uint64_t collBC = collision.template bc_as<BC>().globalBC();
      // int collIdx = collision.globalIndex();
      int64_t bcOffset = globalBC - static_cast<int64_t>(collBC);
      if (static_cast<uint64_t>(std::abs(bcOffset)) > bOffsetMax) {
        if (bcOffset < 0) {
          break;
        } else if (bcOffset > 0) {
          continue;
        }
      }

      float trackTime{0.};
      float trackTimeRes{0.};
      if (trackCand.isPVContributor()) {
        trackTime = trackCand.template collision_as<C>().collisionTime(); // if PV contributor, we assume the time to be the one of the collision
        trackTimeRes = o2::constants::lhc::LHCBunchSpacingNS;             // 1 BC
      } else {
        trackTime = trackCand.trackTime();
        trackTimeRes = trackCand.trackTimeRes();
      }

      const float deltaTime = trackTime - collTime + bcOffset * o2::constants::lhc::LHCBunchSpacingNS;
      float sigmaTimeRes2 = collTimeRes2 + trackTimeRes * trackTimeRes;

      float thresholdTime = 0.;
      if (trackCand.isPVContributor()) {
        thresholdTime = trackTimeRes;
      } else if (TESTBIT(trackCand.flags(), o2::aod::track::TrackTimeResIsRange)) {
        thresholdTime = std::sqrt(sigmaTimeRes2);
        thresholdTime += timeMarginNS;
      } else {
        thresholdTime = 4. * std::sqrt(sigmaTimeRes2);
        thresholdTime += timeMarginNS;
      }
      // LOG(info) << "Threshold time: " << thresholdTime << " isPVContributor: " << trackCand.isPVContributor() << " time margin: " << timeMarginNS;
      if (std::abs(deltaTime) > thresholdTime) {
        continue;
      }

      const auto& tref = tmap.find(trackCand.globalIndex());
      if (tref != tmap.end()) {
        LOG(debug) << "Track: " << trackCand.globalIndex() << " already processed with other vertex";
        trackCandPool[tref->second.second][tref->second.first].collBracket.setMax(static_cast<int>(collIdx)); // this track was already processed with other vertex, account the latter
        continue;
      }

      int poolIndex = (1 - isDau0) * 2 + (trackCand.sign() < 0);
      trForpool.Idxtr = trackCand.globalIndex();
      trForpool.collBracket = {static_cast<int>(collIdx), static_cast<int>(collIdx)};
      // LOG(info) << "Adding track to pool: " << trForpool.Idxtr << " with bracket: " << trForpool.collBracket.getMin() << " " << trForpool.collBracket.getMax() << " and pool index: " << poolIndex;
      trackCandPool[poolIndex].emplace_back(trForpool);
      tmap[trackCand.globalIndex()] = {trackCandPool[poolIndex].size() - 1, poolIndex};
    }

    // is Sorting Needed ? TBD
  }

  template <typename C>
  std::vector<SVCand>& getSVCandPool(const C& collisions, bool combineLikeSign = false)
  {
    gsl::span<std::vector<TrackCand>> track0Pool{trackCandPool.data(), NChargeSigns};
    gsl::span<std::vector<TrackCand>> track1Pool{trackCandPool.data() + NChargeSigns, NChargeSigns};
    std::array<std::vector<int>, NChargeSigns> mVtxTrack0{}; // 1st pos. and neg. track of the kink pool for each vertex

    for (int i = 0; i < NChargeSigns; i++) {
      mVtxTrack0[i].clear();
      mVtxTrack0[i].resize(collisions.size(), -1);
    }
    for (int iCharge = 0; iCharge < NChargeSigns; iCharge++) {
      auto& vtxFirstT = mVtxTrack0[iCharge];
      const auto& signTrack0Pool = track0Pool[iCharge];
      for (unsigned i = 0; i < signTrack0Pool.size(); i++) {
        const auto& track0Seed = signTrack0Pool[i];
        for (int j{track0Seed.collBracket.getMin()}; j <= track0Seed.collBracket.getMax(); ++j) {
          if (vtxFirstT[j] == -1) {
            vtxFirstT[j] = i;
          }
        }
      }
      int track1sign = combineLikeSign ? iCharge : 1 - iCharge;
      auto& signTrack1 = track1Pool[track1sign];
      for (unsigned iTrack1 = 0; iTrack1 < signTrack1.size(); iTrack1++) {
        auto& track1Seed = signTrack1[iTrack1];
        int firstOverlapIdx = -1;
        for (int j{track1Seed.collBracket.getMin()}; j <= track1Seed.collBracket.getMax(); ++j) {
          LOG(debug) << "Checking vtxFirstT at position " << j << " with value " << vtxFirstT[j];
          if (vtxFirstT[j] != -1) {
            firstOverlapIdx = vtxFirstT[j];
            break;
          }
        }
        if (firstOverlapIdx < 0) {
          continue;
        }
        for (unsigned iTrack0 = firstOverlapIdx; iTrack0 < signTrack0Pool.size(); iTrack0++) {
          auto& track0Seed = signTrack0Pool[iTrack0];

          if (track0Seed.collBracket.getMin() > track1Seed.collBracket.getMax()) {
            break;
          }

          if (track0Seed.collBracket.isOutside(track1Seed.collBracket)) {
            LOG(debug) << "Brackets do not match";
            continue;
          }
          auto overlapBracket = track0Seed.collBracket.getOverlap(track1Seed.collBracket);

          svCandPool.emplace_back(SVCand{track0Seed.Idxtr, track1Seed.Idxtr, overlapBracket});
        }
      }
    }
    return svCandPool;
  }

  template <typename T>
  bool fitSV(unsigned int idxDau0, unsigned int idxDau1, T& trackTable);

 private:
  o2::vertexing::DCAFitterN<2> fitter;
  int track0Pdg;
  int track1Pdg;
  float timeMarginNS = 600.;
  bool skipAmbiTracks = false;
  std::unordered_map<int, std::pair<int, int>> tmap;
  std::unordered_map<uint64_t, int> bc2Coll;

  std::array<std::vector<TrackCand>, 4> trackCandPool; // Sorting: dau0 pos, dau0 neg, dau1 pos, dau1 neg
  std::vector<SVCand> svCandPool;                      // index of the two tracks in the track table
  TrackCand trForpool;
};

/// Build time-compatible three-track combinations by joining two svPoolCreator
/// pools on a common first prong. The class deliberately only handles
/// combinatorics; the caller remains responsible for fitting and physics
/// selections.
class svPoolCreator3Body
{
 public:
  svPoolCreator3Body() = default;
  svPoolCreator3Body(int track0Pdg, int track1Pdg, int track2Pdg)
    : track0Pdg(track0Pdg),
      track1Pdg(track1Pdg),
      track2Pdg(track2Pdg),
      pool01(track0Pdg, track1Pdg),
      pool02(track0Pdg, track2Pdg)
  {
  }

  void setPDGs(int pdg0, int pdg1, int pdg2)
  {
    track0Pdg = pdg0;
    track1Pdg = pdg1;
    track2Pdg = pdg2;
    pool01.setPDGs(track0Pdg, track1Pdg);
    pool02.setPDGs(track0Pdg, track2Pdg);
  }

  void clearPools()
  {
    pool01.clearPools();
    pool02.clearPools();
    svCandPool.clear();
  }

  void setTimeMargin(float timeMargin)
  {
    pool01.setTimeMargin(timeMargin);
    pool02.setTimeMargin(timeMargin);
  }

  void setSkipAmbiTracks()
  {
    pool01.setSkipAmbiTracks();
    pool02.setSkipAmbiTracks();
  }

  template <typename C, typename BC>
  void fillBC2Coll(const C& collisions, BC const& bcs)
  {
    pool01.fillBC2Coll(collisions, bcs);
    pool02.fillBC2Coll(collisions, bcs);
  }

  template <typename T, typename C, typename AT, typename BC>
  void appendTrackCand(const T& trackCand, const C& collisions, int pdgHypo, const AT& ambiTracks, BC const& bcs)
  {
    if (pdgHypo == track0Pdg) {
      pool01.appendTrackCand(trackCand, collisions, pdgHypo, ambiTracks, bcs);
      pool02.appendTrackCand(trackCand, collisions, pdgHypo, ambiTracks, bcs);
    } else if (pdgHypo == track1Pdg) {
      pool01.appendTrackCand(trackCand, collisions, pdgHypo, ambiTracks, bcs);
    } else if (pdgHypo == track2Pdg) {
      pool02.appendTrackCand(trackCand, collisions, pdgHypo, ambiTracks, bcs);
    } else {
      LOGP(debug, "Wrong PDG hypothesis for three-body pool");
    }
  }

  template <typename C>
  std::vector<SVCand3>& getSVCandPool(const C& collisions, bool combineLikeSign01 = false, bool combineLikeSign02 = false)
  {
    auto& candidates01 = pool01.getSVCandPool(collisions, combineLikeSign01);
    auto& candidates02 = pool02.getSVCandPool(collisions, combineLikeSign02);
    std::unordered_multimap<int, const SVCand*> candidates02ByTrack0;
    candidates02ByTrack0.reserve(candidates02.size());
    for (const auto& candidate02 : candidates02) {
      candidates02ByTrack0.emplace(candidate02.tr0Idx, &candidate02);
    }

    for (const auto& candidate01 : candidates01) {
      const auto [begin, end] = candidates02ByTrack0.equal_range(candidate01.tr0Idx);
      for (auto candidate02It = begin; candidate02It != end; ++candidate02It) {
        const auto& candidate02 = *candidate02It->second;
        if (candidate01.tr1Idx == candidate02.tr1Idx ||
            candidate01.collBracket.isOutside(candidate02.collBracket)) {
          continue;
        }
        svCandPool.emplace_back(SVCand3{candidate01.tr0Idx,
                                        candidate01.tr1Idx,
                                        candidate02.tr1Idx,
                                        candidate01.collBracket.getOverlap(candidate02.collBracket)});
      }
    }
    return svCandPool;
  }

 private:
  int track0Pdg = 0;
  int track1Pdg = 0;
  int track2Pdg = 0;
  svPoolCreator pool01;
  svPoolCreator pool02;
  std::vector<SVCand3> svCandPool;
};

#endif // PWGLF_UTILS_SVPOOLCREATOR_H_
