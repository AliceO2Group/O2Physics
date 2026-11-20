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

#ifndef ALICE3_CORE_FLATTRACKSMEARER_H_
#define ALICE3_CORE_FLATTRACKSMEARER_H_

#include "FlatLutEntry.h"

#include <CCDB/BasicCCDBManager.h>
#include <ReconstructionDataFormats/Track.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <span>

namespace o2::fastsim
{
/**
 * @brief Track smearing with flat LUT backend
 */
using O2Track = o2::track::TrackParCov;
class TrackSmearer
{
 public:
  TrackSmearer() = default;

  /** LUT methods **/
  bool loadTable(int pdg, const char* filename, bool forceReload = false);
  bool adoptTable(int pdg, const uint8_t* buffer, size_t size, bool forceReload = false);
  bool viewTable(int pdg, const uint8_t* buffer, size_t size, bool forceReload = false);
  bool viewTable(int pdg, std::span<std::byte> const& span, bool forceReload = false);
  [[nodiscard]] bool hasTable(int pdg) const;

  void useEfficiency(bool val) { mUseEfficiency = val; }
  void interpolateEfficiency(bool val) { mInterpolateEfficiency = val; }
  void skipUnreconstructed(bool val) { mSkipUnreconstructed = val; }
  enum EfficiencyType : uint8_t {
    kWhatEfficiencyReco = 0,
    kWhatEfficiencyRecoAndTOF = 1
  };
  void setWhatEfficiency(EfficiencyType val) { mWhatEfficiency = val; }

  [[nodiscard]] const lutHeader_t* getLUTHeader(const int pdg) const;
  const lutEntry_t* getLUTEntry(const int pdg, const float nch, const float radius, const float eta, const float pt, float& interpolatedEff) const;

  bool smearTrack(O2Track& o2track, const lutEntry_t* lutEntry, const float interpolatedEff);
  bool smearTrack(O2Track& o2track, const int pdg, const float nch);

  [[nodiscard]] double getPtRes(const int pdg, const float nch, const float eta, const float pt) const;
  [[nodiscard]] double getEtaRes(const int pdg, const float nch, const float eta, const float pt) const;
  [[nodiscard]] double getAbsPtRes(const int pdg, const float nch, const float eta, const float pt) const;
  [[nodiscard]] double getAbsEtaRes(const int pdg, const float nch, const float eta, const float pt) const;
  [[nodiscard]] double getEfficiency(const int pdg, const float nch, const float eta, const float pt) const;

  static int getIndexPDG(const int pdg);
  static const char* getParticleName(const int pdg);

  void setCcdbManager(o2::ccdb::BasicCCDBManager* mgr) { mCcdbManager = mgr; }

 protected:
  static constexpr unsigned int nLUTs = 9; // Number of LUT available, une per particle species (electron, muon, pion, kaon, proton, deuteron, triton, helium3, alpha)
  std::array<FlatLutData, nLUTs> mLUTData; // Flat data storage

  bool mUseEfficiency = true;
  bool mInterpolateEfficiency = false;
  bool mSkipUnreconstructed = true; // don't smear tracks that are not reco'ed
  EfficiencyType mWhatEfficiency = kWhatEfficiencyReco;

 private:
  o2::ccdb::BasicCCDBManager* mCcdbManager = nullptr;

  static bool checkSpecialCase(const int pdg, lutHeader_t const& header);
};

} // namespace o2::fastsim

#endif // ALICE3_CORE_FLATTRACKSMEARER_H_
