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

#ifndef PWGEM_PHOTONMESON_UTILS_EVENTMIXINGHANDLER_H_
#define PWGEM_PHOTONMESON_UTILS_EVENTMIXINGHANDLER_H_

#include <map>
#include <utility>
#include <vector>

namespace o2::aod::pwgem::utils
{
template <typename T, typename U, typename V>
class EventMixingHandler : public TNamed
{
  EventMixingHandler(const char* name, const char* title) : TNamed(name, title)
  {
    fNdepth = 0;
  }

  EventMixingHandler(const char* name, const char* title, int ndepth) : TNamed(name, title)
  {
    fNdepth = ndepth;
  }

  ~EventMixingHandler()
  {
    fMapMixBins.clear();
    fMap_Tracks_per_collision.clear();
  }

  void SetNdepth(int ndepth) { fNdepth = ndepth; }

  void AddLast(U key_df_collision, V obj)
  {
    fMap_Tracks_per_collision[key_df_collision].emplace_back(obj);
  }

  std::vector<U> GetCollisionIds(T key_bin) { return fMapMixBins[key_bin]; }
  std::vector<V> GetTracksPerCollision(T key_bin, int index) { return fMap_Tracks_per_collision[fMapMixBins[key_bin][index]]; }
  std::vector<V> GetTracksPerCollision(U key_df_collision) { return fMap_Tracks_per_collision[key_df_collision]; }

  void Update(T key_bin, U key_df_collision)
  {
    if (static_cast<int>(fMapMixBins[key_bin].size()) >= fNdepth) {
      fMap_Tracks_per_collision[fMapMixBins[key_bin][0]].clear();
      fMap_Tracks_per_collision[fMapMixBins[key_bin][0]].shrink_to_fit();
      fMapMixBins[key_bin].erase(fMapMixBins[key_bin].begin());
    }
    fMapMixBins[key_bin].emplace_back(key_df_collision);
  }

 private:
  int fNdepth;                                           // depth of event mixing
  std::map<T, std::vector<U>> fMapMixBins;               // map : e.g. <zbin, centbin, epbin> -> pair<df index, global collision index>
  std::map<U, std::vector<V>> fMap_Tracks_per_collision; // map : e.g. pair<df index, global collision index> -> track array
};
} // namespace o2::aod::pwgem::utils
#endif // PWGEM_PHOTONMESON_UTILS_EVENTMIXINGHANDLER_H_
