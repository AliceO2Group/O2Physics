// Copyright 2019-2022 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file collisionSelection.h
/// \brief collision selection
/// \author Anton Riedel, TU München, anton.riedel@cern.ch

#ifndef PWGCF_FEMTOUNITED_CORE_COLLISIONSELECTION_H_
#define PWGCF_FEMTOUNITED_CORE_COLLISIONSELECTION_H_

#include "PWGCF/FemtoUnited/Core/femtoUtils.h"
#include "PWGCF/FemtoUnited/Core/modes.h"

#include "Framework/Configurable.h"

#include <cmath>
#include <string>

namespace o2::analysis::femtounited
{
namespace collisionselection
{

// configurables for collision selection
struct ConfCollisionSelection : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("CollisionSelection");
  o2::framework::Configurable<bool> useOfflineSelection{"useOfflineSelection", true, "Use offline selections (Run3 -> Sel8)"};
  o2::framework::Configurable<float> vtxZMin{"vtxZMin", -10.f, "Minimum vertex Z position (cm)"};
  o2::framework::Configurable<float> vtxZMax{"vtxZMax", 10.f, "Maximum vertex Z position (cm)"};
  o2::framework::Configurable<float> multMin{"multMin", 0.f, "Minimum multiplicity"};
  o2::framework::Configurable<float> multMax{"multMax", 999.f, "Maximum multiplicity"};
  o2::framework::Configurable<float> centMin{"centMin", 0.f, "Minimum centrality (multiplicity percentile)"};
  o2::framework::Configurable<float> centMax{"centMax", 999.f, "Maximum centrality (multiplicity percentile)"};
  o2::framework::Configurable<float> spherMin{"spherMin", 0.f, "Minimum centrality (multiplicity percentile)"};
  o2::framework::Configurable<float> spherMax{"spherMax", 2.f, "Maximum centrality (multiplicity percentile)"};
  o2::framework::Configurable<float> magFieldMin{"magFieldMin", -1.f, "Minimum magnetic field strength (T)"};
  o2::framework::Configurable<float> magFieldMax{"magFieldMax", 1.f, "Maximum magnetic field strength (T)"};
};

/// \class FemtoDreamTrackCuts
/// \brief Cut class to contain and execute all cuts applied to tracks
class CollisionSelection
{
 public:
  CollisionSelection() {}
  virtual ~CollisionSelection() = default;

  template <typename T>
  void configure(T const& config)
  {
    mOfflineSelection = config.useOfflineSelection.value;
    mMagFieldMin = config.magFieldMin.value;
    mMagFieldMax = config.magFieldMax.value;
    mMultMin = config.multMin.value;
    mMultMax = config.multMax.value;
    mCentMin = config.centMin.value;
    mCentMax = config.centMax.value;
    mSphericityMin = config.spherMin.value;
    mSphericityMax = config.spherMax.value;
  };

  void setMagneticField(float MagField)
  {
    mMagField = MagField;
  }

  float getMagneticField()
  {
    return mMagField;
  }

  template <typename T>
  void setSphericity(T tracks)
  {
    mSphericity = utils::sphericity(tracks);
  }

  float getSphericity()
  {
    return mSphericity;
  }

  template <modes::System system, typename T>
  bool checkCuts(T const& col)
  {
    if constexpr (modes::isFlagSet(system, modes::System::kPP) && modes::isFlagSet(system, modes::System::kRun3)) {
      if (mOfflineSelection && !col.sel8()) {
        return false;
      }
      if (col.multNTracksPV() < mMultMin || col.multNTracksPV() > mMultMax) {
        return false;
      }
      if (col.centFT0M() < mCentMin || col.centFT0M() > mCentMax) {
        return false;
      }
      if (mMagField < mMagFieldMin || mMagField > mMagFieldMax) {
        return false;
      }
      if (mSphericity < mSphericityMin || mSphericity > mSphericityMax) {
        return false;
      }
    }
    return true;
  }

 private:
  bool mOfflineSelection = false;
  float mSphericityMin = 0.f;
  float mSphericityMax = 2.f;
  float mMagFieldMin = -1.f;
  float mMagFieldMax = 1.f;
  float mMultMin = 0.f;
  float mMultMax = 999.f;
  float mCentMin = 0.f;
  float mCentMax = 999.f;

  float mMagField = 0.f;
  float mSphericity = 0.f;
};
}; // namespace collisionselection
}; // namespace o2::analysis::femtounited
#endif // PWGCF_FEMTOUNITED_CORE_COLLISIONSELECTION_H_
