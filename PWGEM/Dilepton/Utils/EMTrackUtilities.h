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

/// \commonly used to calculate track variables
/// \author daiki.sekihata@cern.ch

#ifndef PWGEM_DILEPTON_UTILS_EMTRACKUTILITIES_H_
#define PWGEM_DILEPTON_UTILS_EMTRACKUTILITIES_H_

#include "Framework/DataTypes.h"
#include "Framework/Logger.h"

#include <algorithm>
#include <map>
#include <string>
#include <vector>

//_______________________________________________________________________
namespace o2::aod::pwgem::dilepton::utils::emtrackutil
{

enum class RefTrackBit : uint16_t { // This is not for leptons, but charged tracks for reference flow.
  kNclsITS5 = 1,
  kNclsITS6 = 2,
  kNcrTPC70 = 4,
  kNcrTPC90 = 8,
  kNclsTPC50 = 16, // (not necessary, if ncr is used.)
  kNclsTPC70 = 32, // (not necessary, if ncr is used.)
  kNclsTPC90 = 64, // (not necessary, if ncr is used.)
  kChi2TPC4 = 128,
  kChi2TPC3 = 256,
  kFracSharedTPC07 = 512,
  kDCAxy05cm = 1024, // default is 1 cm
  kDCAxy03cm = 2048,
  kDCAz05cm = 4096, // default is 1cm
  kDCAz03cm = 8192,
};

enum class RefMFTTrackBit : uint16_t { // This is not for leptons, but charged tracks for reference flow.
  kNclsMFT7 = 1,                       // default is 6
  kNclsMFT8 = 2,
  kChi2MFT4 = 4, // default is 5
  kChi2MFT3 = 8,
  kDCAxy004cm = 16, // default is 0.05 cm
  kDCAxy003cm = 32,
  kDCAxy002cm = 64,
  kDCAxy001cm = 128,
};

//_______________________________________________________________________
template <typename T>
float dca3DinSigma(T const& track)
{
  float cYY = track.cYY();
  float cZZ = track.cZZ();
  float cZY = track.cZY();
  float dcaXY = track.dcaXY(); // in cm
  float dcaZ = track.dcaZ();   // in cm

  float det = cYY * cZZ - cZY * cZY; // determinant
  if (det < 0) {
    return 999.f;
  } else {
    return std::sqrt(std::fabs((dcaXY * dcaXY * cZZ + dcaZ * dcaZ * cYY - 2. * dcaXY * dcaZ * cZY) / det / 2.)); // dca 3d in sigma
  }
}
//_______________________________________________________________________
template <typename T>
float sigmaDca3D(T const& track)
{
  float dcaXY = track.dcaXY();                          // in cm
  float dcaZ = track.dcaZ();                            // in cm
  float dca3d = std::sqrt(dcaXY * dcaXY + dcaZ * dcaZ); // in cm
  return dca3d / dca3DinSigma(track);
}
//_______________________________________________________________________
template <typename T>
float dcaXYinSigma(T const& track)
{
  return track.dcaXY() / std::sqrt(track.cYY());
}
//_______________________________________________________________________
template <typename T>
float dcaZinSigma(T const& track)
{
  return track.dcaZ() / std::sqrt(track.cZZ());
}
//_______________________________________________________________________
template <typename T>
float fwdDcaXYinSigma(T const& track)
{
  float cXX = track.cXXatDCA();      // in cm^2
  float cYY = track.cYYatDCA();      // in cm^2
  float cXY = track.cXYatDCA();      // in cm^2
  float dcaX = track.fwdDcaX();      // in cm
  float dcaY = track.fwdDcaY();      // in cm
  float det = cXX * cYY - cXY * cXY; // determinant

  if (det < 0) {
    return 999.f;
  } else {
    return std::sqrt(std::fabs((dcaX * dcaX * cYY + dcaY * dcaY * cXX - 2. * dcaX * dcaY * cXY) / det / 2.)); // dca xy in sigma
  }
}
//_______________________________________________________________________
template <typename T>
float sigmaFwdDcaXY(T const& track)
{
  float dcaX = track.fwdDcaX();                       // in cm
  float dcaY = track.fwdDcaY();                       // in cm
  float dcaXY = std::sqrt(dcaX * dcaX + dcaY * dcaY); // in cm
  return dcaXY / fwdDcaXYinSigma(track);
}
//_______________________________________________________________________
template <int begin = 0, int end = 9, typename T>
bool checkMFTHitMap(T const& track)
{
  // logical-OR
  uint64_t mftClusterSizesAndTrackFlags = track.mftClusterSizesAndTrackFlags();
  uint16_t clmap = 0;
  for (unsigned int layer = begin; layer <= end; layer++) {
    if ((mftClusterSizesAndTrackFlags >> (layer * 6)) & 0x3f) {
      clmap |= (1 << layer);
    }
  }
  return (clmap > 0);
}
//_______________________________________________________________________
template <bool is_wo_acc = false, typename TTrack, typename TCut, typename TTracks>
bool isBestMatch(TTrack const& track, TCut const& cut, TTracks const& tracks)
{
  // this is only for global muons at forward rapidity
  // Be careful! tracks are fwdtracks per DF.
  if (track.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) {
    std::map<int64_t, float> map_chi2MCHMFT;
    map_chi2MCHMFT[track.globalIndex()] = track.chi2MatchMCHMFT(); // add myself
    for (const auto& glmuonId : track.globalMuonsWithSameMFTIds()) {
      const auto& candidate = tracks.rawIteratorAt(glmuonId);
      if (track.mchtrackId() != candidate.mchtrackId() && track.mfttrackId() == candidate.mfttrackId()) {
        if (candidate.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) {
          if (cut.template IsSelectedTrack<is_wo_acc>(candidate)) {
            map_chi2MCHMFT[candidate.globalIndex()] = candidate.chi2MatchMCHMFT();
          }
        }
      }
    } // end of glmuonId

    auto it = std::min_element(map_chi2MCHMFT.begin(), map_chi2MCHMFT.end(), [](decltype(map_chi2MCHMFT)::value_type& l, decltype(map_chi2MCHMFT)::value_type& r) -> bool { return l.second < r.second; }); // search for minimum matching-chi2

    if (it->first == track.globalIndex()) {
      map_chi2MCHMFT.clear();
      return true;
    } else {
      map_chi2MCHMFT.clear();
      return false;
    }
  } else {
    return true;
  }
}
//_______________________________________________________________________
// template <typename T>
// float sigmaPt(T const& track)
// {
//   return std::sqrt(track.c1Pt21Pt2()) / std::pow(track.signed1Pt(), 2); // pT resolution
// }
// //_______________________________________________________________________
// template <typename T>
// float sigmaPhi(T const& track)
// {
//   return std::sqrt(track.cSnpSnp()) / std::sqrt(1.f - std::pow(track.snp(), 2)); // phi resolution
// }
// //_______________________________________________________________________
// template <typename T>
// float sigmaTheta(T const& track)
// {
//   return std::sqrt(track.cTglTgl()) / (1.f + std::pow(track.tgl(), 2)); // theta resolution = lambda resolution. // lambda = pi/2 - theta. theta is polar angle.
// }
// //_______________________________________________________________________
// template <typename T>
// float sigmaEta(T const& track)
// {
//   return std::sqrt(track.cTglTgl()) / std::sqrt(1.f + std::pow(track.tgl(), 2));
// }
// //_______________________________________________________________________
// template <typename T>
// float sigmaP(T const& track)
// {
//   // p = 1/1/pT x 1/cos(lambda);
//   return std::sqrt(std::pow(1.f / track.signed1Pt(), 4) * ((1.f + std::pow(track.tgl(), 2)) * track.c1Pt21Pt2() + 1.f / (1.f + std::pow(track.tgl(), 2)) * std::pow(track.signed1Pt() * track.tgl(), 2) * track.cTglTgl() - 2.f * track.signed1Pt() * track.tgl() * track.c1PtTgl()));
// }
//_______________________________________________________________________
} // namespace o2::aod::pwgem::dilepton::utils::emtrackutil
#endif // PWGEM_DILEPTON_UTILS_EMTRACKUTILITIES_H_
