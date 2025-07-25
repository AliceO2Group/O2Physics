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

#include <set>
#include <vector>
#include <utility>
#include <string>
#include "Framework/Logger.h"
#include "Framework/DataTypes.h"
#include "Rtypes.h"
#include "TNamed.h"
#include "TMath.h"

class PairCut : public TNamed
{
 public:
  PairCut() = default;
  PairCut(const char* name, const char* title) : TNamed(name, title) {}

  enum class PairCuts : int {
    // v0 cut
    kAsym = 0,
    kNCuts
  };

  static const char* mCutNames[static_cast<int>(PairCuts::kNCuts)];

  template <typename G1, typename G2>
  bool IsSelected(G1 const& g1, G2 const& g2) const
  {
    if (!IsSelectedPair(g1, g2, PairCuts::kAsym)) {
      return false;
    }

    return true;
  }

  // template <typename U1, typename U2, typename G1, typename G2>
  template <typename G1, typename G2>
  bool IsSelectedPair(G1 const& g1, G2 const& g2, const PairCuts& cut) const
  {
    switch (cut) {
      case PairCuts::kAsym: {
        float asym = abs(g1.e() - g2.e()) / (g1.e() + g2.e());
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
  float mMinAsym{-1e+10}, mMaxAsym{1e+10};

  ClassDef(PairCut, 1);
};

#endif // PWGEM_PHOTONMESON_CORE_PAIRCUT_H_
