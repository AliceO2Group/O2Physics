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

#ifndef PWGLF_UTILS_SVPOOLCREATOR_H_
#define PWGLF_UTILS_SVPOOLCREATOR_H_

#include <array>
#include <unordered_map>
#include <vector>
#include <utility>
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Common/Core/trackUtilities.h"
#include "DCAFitter/DCAFitterN.h"
#include "Framework/AnalysisDataModel.h"

using namespace o2;
using namespace o2::constants;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;
using CollBracket = o2::math_utils::Bracket<int>;

constexpr uint64_t bOffsetMax = 241; // track compatibility can never go beyond 6 mus (ITS)

struct TrackCand {
  int Idxtr;
  CollBracket collBracket{};
};

struct SVCand {
  int tr0Idx;
  int tr1Idx;
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
    for (auto& pool : trackCandPool) {
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

  template <typename T, typename C, typename BC>
  void appendTrackCand(const T& trackCand, const C& collisions, int pdgHypo, o2::aod::AmbiguousTracks const& ambiTracks, BC const&)
  {
    if (pdgHypo != track0Pdg && pdgHypo != track1Pdg) {
      LOG(debug) << "Wrong pdg hypothesis";
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
        if (!ambTrack.has_bc() || ambTrack.bc_as<BC>().size() == 0) {
          globalBC = BcInvalid;
          break;
        }
        globalBC = ambTrack.bc_as<BC>().begin().globalBC();
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
    for (int i = firstCollIdx; i < collisions.size(); i++) {
      const auto& collision = collisions.rawIteratorAt(i);
      float collTime = collision.collisionTime();
      float collTimeRes2 = collision.collisionTimeRes() * collision.collisionTimeRes();
      uint64_t collBC = collision.template bc_as<BC>().globalBC();
      // int collIdx = collision.globalIndex();
      int collIdx = i;
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
        trackTimeRes = constants::lhc::LHCBunchSpacingNS;                 // 1 BC
      } else {
        trackTime = trackCand.trackTime();
        trackTimeRes = trackCand.trackTimeRes();
      }

      const float deltaTime = trackTime - collTime + bcOffset * constants::lhc::LHCBunchSpacingNS;
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
    gsl::span<std::vector<TrackCand>> track0Pool{trackCandPool.data(), 2};
    gsl::span<std::vector<TrackCand>> track1Pool{trackCandPool.data() + 2, 2};
    std::array<std::vector<int>, 2> mVtxTrack0{}; // 1st pos. and neg. track of the kink pool for each vertex

    for (int i = 0; i < 2; i++) {
      mVtxTrack0[i].clear();
      mVtxTrack0[i].resize(collisions.size(), -1);
    }

    for (int pn = 0; pn < 2; pn++) {
      auto& vtxFirstT = mVtxTrack0[pn];
      const auto& signTrack0Pool = track0Pool[pn];
      for (unsigned i = 0; i < signTrack0Pool.size(); i++) {
        const auto& track0Seed = signTrack0Pool[i];
        LOG(debug) << "Processsing track0 with index: " << track0Seed.Idxtr << " min bracket: " << track0Seed.collBracket.getMin() << " max bracket: " << track0Seed.collBracket.getMax();
        for (int j{track0Seed.collBracket.getMin()}; j <= track0Seed.collBracket.getMax(); ++j) {
          if (vtxFirstT[j] == -1) {
            vtxFirstT[j] = i;
          }
        }
      }
      int track1sign = combineLikeSign ? pn : 1 - pn;
      auto& signTrack1 = track1Pool[track1sign];
      for (unsigned itp = 0; itp < signTrack1.size(); itp++) {
        auto& track1Seed = signTrack1[itp];
        LOG(debug) << "Processing track1 with index: " << track1Seed.Idxtr << " min bracket: " << track1Seed.collBracket.getMin() << " max bracket: " << track1Seed.collBracket.getMax();
        int firsOverlapIdx = -1;
        for (int j{track1Seed.collBracket.getMin()}; j <= track1Seed.collBracket.getMax(); ++j) {
          LOG(debug) << "Checking vtxFirstT at position " << j << " with value " << vtxFirstT[j];
          if (vtxFirstT[j] != -1) {
            firsOverlapIdx = vtxFirstT[j];
            break;
          }
        }
        if (firsOverlapIdx < 0) {
          continue;
        }
        for (unsigned itn = firsOverlapIdx; itn < signTrack0Pool.size(); itn++) {
          auto& track0Seed = signTrack0Pool[itn];

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

#endif // PWGLF_UTILS_SVPOOLCREATOR_H_
