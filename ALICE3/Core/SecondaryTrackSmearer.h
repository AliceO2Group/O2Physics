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

#ifndef ALICE3_CORE_SECONDARYTRACKSMEARER_H_
#define ALICE3_CORE_SECONDARYTRACKSMEARER_H_

#include "ALICE3/Core/FlatTrackSmearer.h"
#include "ALICE3/Core/SecondaryLut.h"

#include <memory>
#include <unordered_map>

namespace o2::fastsim
{

/// Track smearing using ROOT-persistent secondary-particle LUTs.
class SecondaryTrackSmearer
{
 public:
  using EfficiencyType = TrackSmearer::EfficiencyType;
  static constexpr EfficiencyType kWhatEfficiencyReco = TrackSmearer::kWhatEfficiencyReco;
  static constexpr EfficiencyType kWhatEfficiencyRecoAndTOF = TrackSmearer::kWhatEfficiencyRecoAndTOF;

  /// Load an object from a ROOT file (a local path or a ccdb: path).
  bool loadTable(int pdg, const char* filename, bool forceReload = false, const char* objectName = "SecondaryLut");
  /// Copy an in-memory table; subsequent changes to the original do not affect the smearer.
  bool adoptTable(int pdg, const SecondaryLut& table, bool forceReload = false);
  [[nodiscard]] bool hasTable(int pdg) const;
  [[nodiscard]] const SecondaryLut* getLUT(int pdg) const;

  void useEfficiency(bool value)
  {
    mUseEfficiency = value;
    mSmearer.useEfficiency(value);
  }
  void interpolateEfficiency(bool value)
  {
    mInterpolateEfficiency = value;
    mSmearer.interpolateEfficiency(value);
  }
  void skipUnreconstructed(bool value)
  {
    mSkipUnreconstructed = value;
    mSmearer.skipUnreconstructed(value);
  }
  void setWhatEfficiency(EfficiencyType value);
  void setCcdbManager(o2::ccdb::BasicCCDBManager* manager) { mCcdbManager = manager; }

  /// Lookup clamps each coordinate to the table range; nHits must be nonnegative.
  /// The returned pointer is invalidated when its table is reloaded or destroyed.
  [[nodiscard]] const SecondaryLutEntry* getLUTEntry(int pdg, float nch, float radius, float eta, float pt, int nHits, float& interpolatedEff) const;
  /// wasSmeared reports whether the track was modified, including rejected
  /// tracks when skipUnreconstructed is false.
  bool smearTrack(O2Track& track, const SecondaryLutEntry* entry, float interpolatedEff, bool* wasSmeared = nullptr);
  bool smearTrack(O2Track& track, int pdg, float nch, float radius, int nHits, bool* wasSmeared = nullptr);

  /// Queries return NaN if no valid entry is available.
  [[nodiscard]] double getEfficiency(int pdg, float nch, float radius, float eta, float pt, int nHits) const;
  [[nodiscard]] double getPtRes(int pdg, float nch, float radius, float eta, float pt, int nHits) const;
  [[nodiscard]] double getEtaRes(int pdg, float nch, float radius, float eta, float pt, int nHits) const;
  [[nodiscard]] double getAbsPtRes(int pdg, float nch, float radius, float eta, float pt, int nHits) const;
  [[nodiscard]] double getAbsEtaRes(int pdg, float nch, float radius, float eta, float pt, int nHits) const;

 private:
  bool installTable(int pdg, std::unique_ptr<SecondaryLut> table, bool forceReload);
  [[nodiscard]] float efficiency(const SecondaryLutEntry& entry) const;

  std::unordered_map<int, std::unique_ptr<SecondaryLut>> mTables;
  TrackSmearer mSmearer;
  o2::ccdb::BasicCCDBManager* mCcdbManager = nullptr;
  bool mUseEfficiency = true;
  bool mInterpolateEfficiency = false;
  bool mSkipUnreconstructed = true;
  EfficiencyType mWhatEfficiency = kWhatEfficiencyReco;
};

} // namespace o2::fastsim

#endif // ALICE3_CORE_SECONDARYTRACKSMEARER_H_
