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

/// \event mixing handler
/// \author daiki.sekihata@cern.ch

#ifndef PWGEM_DILEPTON_UTILS_EVENTMIXINGHANDLER_H_
#define PWGEM_DILEPTON_UTILS_EVENTMIXINGHANDLER_H_

#include <map>
#include <optional>
#include <stdexcept>
#include <vector>

namespace o2::aod::pwgem::dilepton::utils
{
template <typename T, typename U, typename V>
class EventMixingHandler
{
 public:
  EventMixingHandler()
  {
    fNdepth = 0;
    fMapMixBins.clear();
    fMap_Tracks_per_collision.clear();
  }

  explicit EventMixingHandler(int ndepth)
  {
    fNdepth = ndepth;
    fMapMixBins.clear();
    fMap_Tracks_per_collision.clear();
    if (fNdepth <= 0) {
      throw std::invalid_argument("mixing depth must be positive");
    }
  }

  ~EventMixingHandler()
  {
    fMapMixBins.clear();
    fMap_Tracks_per_collision.clear();
  }

  void SetNdepth(int ndepth)
  {
    if (ndepth <= 0) {
      throw std::invalid_argument("mixing depth must be positive");
    }
    fNdepth = ndepth;
  }

  void ReserveNTracksPerCollision(U key_df_collision, int ntrack)
  {
    fMap_Tracks_per_collision[key_df_collision].reserve(ntrack);
  }

  void AddTrackToEventPool(U key_df_collision, V obj)
  {
    fMap_Tracks_per_collision[key_df_collision].emplace_back(obj);
  }

  const std::vector<U>& GetCollisionIdsFromEventPool(const T& key_bin) const
  {
    const auto it = fMapMixBins.find(key_bin);
    static const std::vector<U> empty;
    return it == fMapMixBins.end() ? empty : it->second;
  }

  const std::vector<V>& GetTracksPerCollision(const U& key) const
  {
    const auto it = fMap_Tracks_per_collision.find(key);
    static const std::vector<V> empty;
    return it == fMap_Tracks_per_collision.end() ? empty : it->second;
  }

  // call this function at the end of collision loop
  std::optional<U> AddCollisionIdAtLast(T key_bin, U key_df_collision)
  {
    // LOGF(info, "fMapMixBins[key_bin].size() = %d", fMapMixBins[key_bin].size());
    auto& ids = fMapMixBins[key_bin];
    std::optional<U> evicted;
    if (static_cast<int>(ids.size()) >= fNdepth) {
      evicted = ids.front();
      fMap_Tracks_per_collision.erase(*evicted);
      ids.erase(ids.begin());
    }
    ids.emplace_back(key_df_collision);
    return evicted;
  }

 private:
  int fNdepth;                                           // depth of event mixing
  std::map<T, std::vector<U>> fMapMixBins;               // map : e.g. <zbin, centbin, epbin> -> pair<df index, global collision index>
  std::map<U, std::vector<V>> fMap_Tracks_per_collision; // map : e.g. pair<df index, global collision index> -> track array
};
} // namespace o2::aod::pwgem::dilepton::utils
#endif // PWGEM_DILEPTON_UTILS_EVENTMIXINGHANDLER_H_
