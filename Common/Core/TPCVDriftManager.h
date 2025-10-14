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

#ifndef COMMON_CORE_TPCVDRIFTMANAGER_H_
#define COMMON_CORE_TPCVDRIFTMANAGER_H_

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/LHCConstants.h>
#include <DataFormatsTPC/VDriftCorrFact.h>
#include <Framework/DataTypes.h>
#include <Framework/Logger.h>
#include <ReconstructionDataFormats/Track.h>
#include <ReconstructionDataFormats/TrackParametrization.h>

#include <cstdint>
#include <string>

namespace o2::aod::common
{

// Thin wrapper for vdrift ccdb queries should partially mirror VDriftHelper class.
// Allows to move TPC standalone tracks under the assumption of a different
// collision than the track is associated to.
class TPCVDriftManager
{
 public:
  void init(o2::ccdb::BasicCCDBManager* ccdb) noexcept
  {
    mCCDB = ccdb;
  }

  void update(uint64_t timestamp) noexcept
  {
    // Check validity of already present obj, otherwise update
    if (mVD != nullptr && (timestamp > static_cast<uint64_t>(mVD->firstTime) || timestamp < static_cast<uint64_t>(mVD->lastTime))) {
      return;
    }

    // Update Obj
    mVD = mCCDB->getForTimeStamp<o2::tpc::VDriftCorrFact>("TPC/Calib/VDriftTgl", timestamp);
    if (mVD == nullptr || mVD->firstTime < 0 || mVD->lastTime < 0) {
      LOGP(error, "Got invalid VDriftCorrFact for {}", timestamp);
      mValid = false;
      return;
    }

    // TODO account for laser calib

    // Update factors
    mTPCVDriftNS = mVD->refVDrift * mVD->corrFact * 1e-3;

    mValid = true;
    LOGP(info, "Updated VDrift for timestamp {} with vdrift={:.7f} (cm/ns)", mVD->creationTime, mTPCVDriftNS);
  }

  template <typename BCs, typename Collisions, typename Collision, typename TrackExtra, typename Track>
  [[nodiscard]] bool moveTPCTrack(const Collision& col, const TrackExtra& trackExtra, Track& track) noexcept
  {
    ++mCalls;

    // Check if there is a good object available otherwise pretend everything is fine
    if (!mValid) {
      if (mInvalid < mWarningLimit) {
        LOGP(warn, "No VDrift object available, pretending track to be correct");
        if (mInvalid == mWarningLimit - 1) {
          LOGP(warn, "Silencing further warnings!");
        }
      }
      ++mInvalid;
      return true;
    }

    // track is fine, or cannot be moved has information is not available
    if (!(trackExtra.flags() & o2::aod::track::TrackFlags::TrackTimeAsym)) {
      ++mNoFlag;
      return true;
    }

    // TPC time is given relative to the closest BC in ns
    float tTB, tTBErr;
    if (col.collisionTimeRes() < 0.f) { // use track data
      ++mColResNeg;
      tTB = trackExtra.trackTime();
      o2::aod::track::extensions::TPCTimeErrEncoding enc;
      enc.encoding.timeErr = trackExtra.trackTimeRes();
      tTBErr = 0.5f * (enc.getDeltaTFwd() + enc.getDeltaTBwd());
    } else {
      ++mColResPos;
      // The TPC track can be associated to a different BC than the one the collision under assumption is;
      // we need to calculate the difference and subtract this from the trackTime()
      const auto trackBC = trackExtra.template collision_as<Collisions>().template foundBC_as<BCs>().globalBC();
      const auto colBC = col.template foundBC_as<BCs>().globalBC();
      float sign{1.f};
      uint64_t diffBC{0};
      if (colBC < trackBC) {
        sign = -1.f;
        diffBC = (trackBC - colBC);
      } else {
        diffBC = (colBC - trackBC);
      }
      float diffBCNS = sign * static_cast<float>(diffBC) * static_cast<float>(o2::constants::lhc::LHCBunchSpacingNS);
      tTB = col.collisionTime() + diffBCNS;
      tTBErr = col.collisionTimeRes();
    }
    float dTime = tTB - trackExtra.trackTime();
    float dDrift = dTime * mTPCVDriftNS;
    float dDriftErr = tTBErr * mTPCVDriftNS;
    if (dDriftErr < 0.f || dDrift > 250.f) { // we cannot move a track outside the drift volume
      if (mOutside < mWarningLimit) {
        LOGP(warn, "Skipping correction outside of tpc volume with dDrift={} +- {}", dDrift, dDriftErr);
        const auto trackBC = trackExtra.template collision_as<Collisions>().template foundBC_as<BCs>().globalBC();
        const auto colBC = col.template foundBC_as<BCs>().globalBC();
        int diffBC = colBC - trackBC;
        LOGP(info, "ct={}; ctr={}; tTB={}; t0={}; dTime={}; dDrift={}; tgl={}:   colBC={}   trackBC={}  diffBC={}", col.collisionTime(), col.collisionTimeRes(), tTB, trackExtra.trackTime(), dTime, dDrift, track.getTgl(), colBC, trackBC, diffBC);
        if (mOutside == mWarningLimit - 1) {
          LOGP(warn, "Silencing further warnings!");
        }
      }
      ++mOutside;
      return false;
    }

    // impose new Z coordinate
    track.setZ(track.getZ() + ((track.getTgl() < 0.) ? -dDrift : dDrift));
    if constexpr (std::is_base_of_v<o2::track::TrackParCov, Track>) {
      track.setCov(track.getSigmaZ2() + dDriftErr * dDriftErr, o2::track::kSigZ2);
    }

    ++mMovedTrks;

    return true;
  }

  void print() noexcept
  {
    LOGP(info, "TPC corrections called: {}; Moved Tracks: {}; Constrained Tracks={}; No Flag: {}; NULL: {}; Outside: {}; ColResPos {}; ColResNeg {};", mCalls, mMovedTrks, mConstrained, mNoFlag, mInvalid, mOutside, mColResPos, mColResNeg);
  }

 private:
  bool mValid{false};
  // Factors
  float mTPCVDriftNS{0.f}; // drift velocity in cm/ns

  // CCDB
  const o2::tpc::VDriftCorrFact* mVD{}; // reference to drift correction
  o2::ccdb::BasicCCDBManager* mCCDB{};  // reference to initialized ccdb manager

  static constexpr unsigned int mWarningLimit{10};

  // Counters
  unsigned int mCalls{0};       // total number of calls
  unsigned int mMovedTrks{0};   // number of moved tracks
  unsigned int mInvalid{0};     // number of tracks where no drift object was available
  unsigned int mColResNeg{0};   // number of collisions with negative resolution
  unsigned int mColResPos{0};   // number of collisions with positive resolution
  unsigned int mNoFlag{0};      // number of tracks without flag set
  unsigned int mOutside{0};     // number of tracks moved but outside of sensible volume
  unsigned int mConstrained{0}; // number of constrained tracks
};

} // namespace o2::aod::common

#endif // COMMON_CORE_TPCVDRIFTMANAGER_H_
