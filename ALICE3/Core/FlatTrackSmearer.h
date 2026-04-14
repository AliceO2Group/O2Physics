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

namespace o2::delphes
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
  bool hasTable(int pdg) const;

  void useEfficiency(bool val) { mUseEfficiency = val; }
  void interpolateEfficiency(bool val) { mInterpolateEfficiency = val; }
  void skipUnreconstructed(bool val) { mSkipUnreconstructed = val; }
  void setWhatEfficiency(int val);

  const lutHeader_t* getLUTHeader(int pdg) const;
  const lutEntry_t* getLUTEntry(int pdg, float nch, float radius, float eta, float pt, float& interpolatedEff) const;

  bool smearTrack(O2Track& o2track, const lutEntry_t* lutEntry, float interpolatedEff);
  bool smearTrack(O2Track& o2track, int pdg, float nch);

  double getPtRes(int pdg, float nch, float eta, float pt) const;
  double getEtaRes(int pdg, float nch, float eta, float pt) const;
  double getAbsPtRes(int pdg, float nch, float eta, float pt) const;
  double getAbsEtaRes(int pdg, float nch, float eta, float pt) const;
  double getEfficiency(int pdg, float nch, float eta, float pt) const;

  static int getIndexPDG(int pdg);
  static const char* getParticleName(int pdg);

  void setdNdEta(float val) { mdNdEta = val; }
  void setCcdbManager(o2::ccdb::BasicCCDBManager* mgr) { mCcdbManager = mgr; }

 protected:
  static constexpr unsigned int nLUTs = 9; // Number of LUT available
  FlatLutData mLUTData[nLUTs];             // Flat data storage

  bool mUseEfficiency = true;
  bool mInterpolateEfficiency = false;
  bool mSkipUnreconstructed = true; // don't smear tracks that are not reco'ed
  int mWhatEfficiency = 1;
  float mdNdEta = 1600.f;

 private:
  o2::ccdb::BasicCCDBManager* mCcdbManager = nullptr;

  static bool checkSpecialCase(int pdg, lutHeader_t const& header);
};

} // namespace o2::delphes

#endif // ALICE3_CORE_FLATTRACKSMEARER_H_
