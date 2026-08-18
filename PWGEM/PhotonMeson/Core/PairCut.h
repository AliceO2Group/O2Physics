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

//
// Class for 2-photon pair selection
//

#ifndef PWGEM_PHOTONMESON_CORE_PAIRCUT_H_
#define PWGEM_PHOTONMESON_CORE_PAIRCUT_H_

#include <array>
#include <cmath>
#include <cstddef>
#include <string>

class PairCut
{
 public:
  PairCut() = default;
  PairCut(const char* name, const char* title) : name(name), title(title) {}

  enum class PairCuts : int {
    // v0 cut
    kAsym = 0,
    kNCuts
  };

  [[nodiscard]] const std::string& getName() const { return name; }
  [[nodiscard]] const std::string& getTitle() const { return title; }

  static const std::array<std::string, static_cast<std::size_t>(PairCuts::kNCuts)> mCutNames;

  template <typename G1, typename G2>
  bool IsSelected(G1 const& g1, G2 const& g2) const
  {
    return IsSelectedPair(g1, g2, PairCuts::kAsym);
  }

  // template <typename U1, typename U2, typename G1, typename G2>
  template <typename G1, typename G2>
  bool IsSelectedPair(G1 const& g1, G2 const& g2, const PairCuts& cut) const
  {
    switch (cut) {
      case PairCuts::kAsym: {
        float asym = std::abs(g1.e() - g2.e()) / (g1.e() + g2.e());
        // float asym = abs(g1.p() - g2.p()) / (g1.p() + g2.p());
        return mMinAsym < asym && asym < mMaxAsym;
      }
      default:
        return false;
    }
  }

  // Setters
  void SetAsymRange(float min = -1e+10f, float max = 1e10f);

  /// @brief Print the pair selection
  void print() const;

 private:
  std::string name;
  std::string title;
  float mMinAsym{-1e+10}, mMaxAsym{1e+10};
};

#endif // PWGEM_PHOTONMESON_CORE_PAIRCUT_H_
