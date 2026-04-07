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

/// \file EMNonLin.h
/// \brief Header file of NonLin class for photons.
/// \author M. Hemmer, marvin.hemmer@cern.ch

#ifndef PWGEM_PHOTONMESON_CORE_EMNONLIN_H_
#define PWGEM_PHOTONMESON_CORE_EMNONLIN_H_

#include <TMatrixD.h> // IWYU pragma: keep (do not replace with TMatrixDfwd.h)
#include <TMatrixDfwd.h>

#include <algorithm>
#include <cstdint>

namespace o2::pwgem::nonlin
{

constexpr int MaxCent = 100.f;

/// \class EMNonLin
/// \brief Class to obtain non linear correction factors for PbPb.
/// Parameters are loaded from CCDB (TMatrixD: rows = iterations, cols = 6 pol1 coefficients).
/// Falls back to hardcoded static table if CCDB object is not available.
/// The correction is of type a + b * refVal^(-c), where a, b and c are themselves described by linear functions over centrality.
class EMNonLin
{
 public:
  static constexpr int MaxIter = 2; // hard cap on iterations
  static constexpr int PhotonN = 3;
  static constexpr int NCols = 6; // a0, a1, b0, b1, c0, c1

  enum class PhotonType : uint8_t {
    kEMC = 0,
    kPCM = 1,
    kPHOS = 2
  };

  struct NonLinParams {
    float a0{1.f}, a1{0.f}; // pol1 params for a: asymptote
    float b0{0.f}, b1{0.f}; // pol1 params for b: magnitude
    float c0{0.f}, c1{0.f}; // pol1 params for c: exponent
  };

  struct Context {
    const NonLinParams* params = nullptr;
    int nIter = 0;
    float cent = 0.f;

    /// \brief Sets parameters for the NonLin. Used with EMNonLin::resolveParams()
    /// \param newParams pointer to new NonLinParams
    void setParams(const NonLinParams* newParams)
    {
      params = newParams;
    }

    /// \brief Sets iteration used for the NonLin.
    /// \param iter iteration
    void setIter(int iter)
    {
      nIter = (iter < 0 || iter >= MaxIter) ? MaxIter - 1 : iter;
    }

    /// \brief Sets current centrality.
    /// \param centrality centrality
    void setCent(float centrality)
    {
      cent = (centrality >= MaxCent) ? MaxCent : centrality;
    }
  };

  /// \brief Load parameters from a TMatrixD fetched from CCDB.
  /// Rows = iterations, cols = {a0, a1, b0, b1, c0, c1}.
  /// Overwrites the static fallback for this photon type.
  /// \param mat  pointer to TMatrixD from CCDB (may be nullptr)
  /// \param type photon type to fill
  void getFromCCDBObject(const TMatrixD* mat, PhotonType type)
  {
    int iType = static_cast<int>(type);
    if (!mat || mat->GetNcols() != NCols) {
      mCCDBLoaded[iType] = false;
      return;
    }
    int nIter = std::min(mat->GetNrows(), MaxIter);
    for (int i = 0; i < nIter; ++i) {
      mCCDBParams[iType][i] = {
        static_cast<float>((*mat)(i, 0)), static_cast<float>((*mat)(i, 1)),
        static_cast<float>((*mat)(i, 2)), static_cast<float>((*mat)(i, 3)),
        static_cast<float>((*mat)(i, 4)), static_cast<float>((*mat)(i, 5))};
    }
    mCCDBLoaded[iType] = true;
  }

  /// \brief Compute the multiplicative correction factor for energy/pT.
  /// \param par energy or pT of the photon
  /// \param ctx context holding centrality, iteration level, and params pointer
  static float getCorrectionFactor(float var, const Context& ctx);

  /// \brief Return pointer to the params array for a given photon type.
  /// \returns CCDB-loaded params if available, otherwise the static fallback.
  const NonLinParams* resolveParams(PhotonType type) const
  {
    int iType = static_cast<int>(type);
    return mCCDBLoaded[iType] ? &mCCDBParams[iType][0] : &FallbackTable[iType][0];
  }

 private:
  NonLinParams mCCDBParams[PhotonN][MaxIter] = {}; // Runtime params loaded from CCDB (per photon type, per iteration)
  bool mCCDBLoaded[PhotonN] = {false, false, false};

  // -------------------------------------------------------
  // Static fallback tables (used when CCDB object is absent)
  // -------------------------------------------------------
  static constexpr NonLinParams FallbackTable[PhotonN][MaxIter] = {
    // kEMC
    {{1.f, 0.f, 0.f, 0.f, 0.f, 0.f},
     {1.f, 0.f, 0.f, 0.f, 0.f, 0.f}},
    // kPCM
    {{1.f, 0.f, 0.010417f, -1.09508e-05f, 0.355795f, 0.00427618f},
     {1.f, 0.f, 0.f, 0.f, 0.f, 0.f}},
    // kPHOS
    {{1.f, 0.f, 0.f, 0.f, 0.f, 0.f},
     {1.f, 0.f, 0.f, 0.f, 0.f, 0.f}},
  };

  static constexpr NonLinParams FallbackTableMC[PhotonN][MaxIter] = {
    // kEMC
    {{1.f, 0.f, 0.f, 0.f, 0.f, 0.f},
     {1.f, 0.f, 0.f, 0.f, 0.f, 0.f}},
    // kPCM
    {{1.f, 0.f, 0.010417f, -1.09508e-05f, 0.355795f, 0.00427618f},
     {1.f, 0.f, 0.f, 0.f, 0.f, 0.f}},
    // kPHOS
    {{1.f, 0.f, 0.f, 0.f, 0.f, 0.f},
     {1.f, 0.f, 0.f, 0.f, 0.f, 0.f}},
  };
};

} // namespace o2::pwgem::nonlin

#endif // PWGEM_PHOTONMESON_CORE_EMNONLIN_H_
