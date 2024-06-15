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

/// \file MlResponseDiletponPair.h
/// \brief Class to compute the ML response for dielectron analyses at the pair level
/// \author Daniel Samitz <daniel.samitz@cern.ch>, SMI Vienna
///         Elisa Meninno, <elisa.meninno@cern.ch>, SMI Vienna

#ifndef PWGEM_DILEPTON_UTILS_MLRESPONSEDIELECTRONPAIR_H_
#define PWGEM_DILEPTON_UTILS_MLRESPONSEDIELECTRONPAIR_H_

#include <map>
#include <string>
#include <vector>

#include "Math/Vector4D.h"
#include "Tools/ML/MlResponse.h"

// Fill the map of available input features
// the key is the feature's name (std::string)
// the value is the corresponding value in EnumInputFeatures
#define FILL_MAP_DIELECTRON_PAIR(FEATURE)                                  \
  {                                                                        \
#FEATURE, static_cast < uint8_t>(InputFeaturesDielectronPair::FEATURE) \
  }

// Check if the index of mCachedIndices (index associated to a FEATURE)
// matches the entry in EnumInputFeatures associated to this FEATURE
// if so, the inputFeatures vector is filled with the FEATURE's value
// by calling the corresponding GETTER from pair fourmomentum v12
#define CHECK_AND_FILL_VEC_DIELECTRON_PAIR(FEATURE, GETTER)          \
  case static_cast<uint8_t>(InputFeaturesDielectronPair::FEATURE): { \
    inputFeatures.emplace_back(v12.GETTER());                        \
    break;                                                           \
  }

// Check if the index of mCachedIndices (index associated to a FEATURE)
// matches the entry in EnumInputFeatures associated to this FEATURE
// if so, the inputFeatures vector is filled with the FEATURE's value
// by calling the corresponding FUNCTION on the tracks t1 and t2
#define CHECK_AND_FILL_VEC_DIELECTRON_PAIR_FUNC(FEATURE, FUNCTION)   \
  case static_cast<uint8_t>(InputFeaturesDielectronPair::FEATURE): { \
    inputFeatures.emplace_back(FUNCTION(t1, t2));                    \
    break;                                                           \
  }

namespace o2::analysis
{
// possible input features for ML
enum class InputFeaturesDielectronPair : uint8_t {
  m,
  pt,
  eta,
  phi,
  phiv,
  pairDcaXY,
  pairDcaZ
};

template <typename TypeOutputScore = float>
class MlResponseDielectronPair : public MlResponse<TypeOutputScore>
{
 public:
  /// Default constructor
  MlResponseDielectronPair() = default;
  /// Default destructor
  virtual ~MlResponseDielectronPair() = default;

  template <typename T>
  float pair_dca_xy(T const& t1, T const& t2)
  {
    return sqrt((pow(t1.dcaXY() / sqrt(t1.cYY()), 2) + pow(t2.dcaXY() / sqrt(t2.cYY()), 2)) / 2.);
  }

  template <typename T>
  float pair_dca_z(T const& t1, T const& t2)
  {
    return sqrt((pow(t1.dcaZ() / sqrt(t1.cZZ()), 2) + pow(t2.dcaZ() / sqrt(t2.cZZ()), 2)) / 2.);
  }

  template <typename T>
  float get_phiv(T const& t1, T const& t2)
  {
    // cos(phiv) = w*a /|w||a|
    // with w = u x v
    // and  a = u x z / |u x z|   , unit vector perpendicular to v12 and z-direction (magnetic field)
    // u = v12 / |v12|            , the unit vector of v12
    // v = v1 x v2 / |v1 x v2|    , unit vector perpendicular to v1 and v2

    ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), o2::constants::physics::MassElectron);
    ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), o2::constants::physics::MassElectron);
    ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;

    bool swapTracks = false;
    if (v1.Pt() < v2.Pt()) { // ordering of track, pt1 > pt2
      ROOT::Math::PtEtaPhiMVector v3 = v1;
      v1 = v2;
      v2 = v3;
      swapTracks = true;
    }

    // momentum of e+ and e- in (ax,ay,az) axis. Note that az=0 by definition.
    // vector product of pep X pem
    float vpx = 0, vpy = 0, vpz = 0;
    if (t1.sign() * t2.sign() > 0) { // Like Sign
      if (!swapTracks) {
        if (d_bz * t1.sign() < 0) {
          vpx = v1.Py() * v2.Pz() - v1.Pz() * v2.Py();
          vpy = v1.Pz() * v2.Px() - v1.Px() * v2.Pz();
          vpz = v1.Px() * v2.Py() - v1.Py() * v2.Px();
        } else {
          vpx = v2.Py() * v1.Pz() - v2.Pz() * v1.Py();
          vpy = v2.Pz() * v1.Px() - v2.Px() * v1.Pz();
          vpz = v2.Px() * v1.Py() - v2.Py() * v1.Px();
        }
      } else { // swaped tracks
        if (d_bz * t2.sign() < 0) {
          vpx = v1.Py() * v2.Pz() - v1.Pz() * v2.Py();
          vpy = v1.Pz() * v2.Px() - v1.Px() * v2.Pz();
          vpz = v1.Px() * v2.Py() - v1.Py() * v2.Px();
        } else {
          vpx = v2.Py() * v1.Pz() - v2.Pz() * v1.Py();
          vpy = v2.Pz() * v1.Px() - v2.Px() * v1.Pz();
          vpz = v2.Px() * v1.Py() - v2.Py() * v1.Px();
        }
      }
    } else { // Unlike Sign
      if (!swapTracks) {
        if (d_bz * t1.sign() > 0) {
          vpx = v1.Py() * v2.Pz() - v1.Pz() * v2.Py();
          vpy = v1.Pz() * v2.Px() - v1.Px() * v2.Pz();
          vpz = v1.Px() * v2.Py() - v1.Py() * v2.Px();
        } else {
          vpx = v2.Py() * v1.Pz() - v2.Pz() * v1.Py();
          vpy = v2.Pz() * v1.Px() - v2.Px() * v1.Pz();
          vpz = v2.Px() * v1.Py() - v2.Py() * v1.Px();
        }
      } else { // swaped tracks
        if (d_bz * t2.sign() > 0) {
          vpx = v1.Py() * v2.Pz() - v1.Pz() * v2.Py();
          vpy = v1.Pz() * v2.Px() - v1.Px() * v2.Pz();
          vpz = v1.Px() * v2.Py() - v1.Py() * v2.Px();
        } else {
          vpx = v2.Py() * v1.Pz() - v2.Pz() * v1.Py();
          vpy = v2.Pz() * v1.Px() - v2.Px() * v1.Pz();
          vpz = v2.Px() * v1.Py() - v2.Py() * v1.Px();
        }
      }
    }

    // unit vector of pep X pem
    float vx = vpx / TMath::Sqrt(vpx * vpx + vpy * vpy + vpz * vpz);
    float vy = vpy / TMath::Sqrt(vpx * vpx + vpy * vpy + vpz * vpz);
    float vz = vpz / TMath::Sqrt(vpx * vpx + vpy * vpy + vpz * vpz);

    float px = v12.Px();
    float py = v12.Py();
    float pz = v12.Pz();

    // unit vector of (pep+pem)
    float ux = px / TMath::Sqrt(px * px + py * py + pz * pz);
    float uy = py / TMath::Sqrt(px * px + py * py + pz * pz);
    float uz = pz / TMath::Sqrt(px * px + py * py + pz * pz);
    float ax = uy / TMath::Sqrt(ux * ux + uy * uy);
    float ay = -ux / TMath::Sqrt(ux * ux + uy * uy);

    // The third axis defined by vector product (ux,uy,uz)X(vx,vy,vz)
    float wx = uy * vz - uz * vy;
    float wy = uz * vx - ux * vz;
    // by construction, (wx,wy,wz) must be a unit vector. Measure angle between (wx,wy,wz) and (ax,ay,0).
    // The angle between them should be small if the pair is conversion. This function then returns values close to pi!
    return TMath::ACos(wx * ax + wy * ay); // phiv in [0,pi] //cosPhiV = wx * ax + wy * ay;
  }

  /// Method to get the input features vector needed for ML inference
  /// \param t1 is the first track
  /// \param t2 is the second track
  /// \return inputFeatures vector
  template <typename T>
  std::vector<float> getInputFeatures(T const& t1, T const& t2)
  {
    std::vector<float> inputFeatures;
    ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), o2::constants::physics::MassElectron);
    ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), o2::constants::physics::MassElectron);
    ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;

    for (const auto& idx : MlResponse<TypeOutputScore>::mCachedIndices) {
      switch (idx) {
        CHECK_AND_FILL_VEC_DIELECTRON_PAIR(m, M);
        CHECK_AND_FILL_VEC_DIELECTRON_PAIR(pt, Pt);
        CHECK_AND_FILL_VEC_DIELECTRON_PAIR(eta, Eta);
        CHECK_AND_FILL_VEC_DIELECTRON_PAIR(phi, Phi);
        CHECK_AND_FILL_VEC_DIELECTRON_PAIR_FUNC(phiv, get_phiv);
        CHECK_AND_FILL_VEC_DIELECTRON_PAIR_FUNC(pairDcaXY, pair_dca_xy);
        CHECK_AND_FILL_VEC_DIELECTRON_PAIR_FUNC(pairDcaZ, pair_dca_z);
      }
    }

    return inputFeatures;
  }

  void setBz(float bz)
  {
    d_bz = bz;
  }

 protected:
  /// Method to fill the map of available input features
  void setAvailableInputFeatures()
  {
    MlResponse<TypeOutputScore>::mAvailableInputFeatures = {
      FILL_MAP_DIELECTRON_PAIR(m),
      FILL_MAP_DIELECTRON_PAIR(pt),
      FILL_MAP_DIELECTRON_PAIR(eta),
      FILL_MAP_DIELECTRON_PAIR(phi),
      FILL_MAP_DIELECTRON_PAIR(phiv),
      FILL_MAP_DIELECTRON_PAIR(pairDcaXY),
      FILL_MAP_DIELECTRON_PAIR(pairDcaZ)};
  }

  float d_bz = 0.;
};

} // namespace o2::analysis

#undef FILL_MAP_DIELECTRON_PAIR
#undef CHECK_AND_FILL_VEC_DIELECTRON_PAIR
#undef CHECK_AND_FILL_VEC_DIELECTRON_PAIR_FUNC

#endif // PWGEM_DILEPTON_UTILS_MLRESPONSEDIELECTRONPAIR_H_
