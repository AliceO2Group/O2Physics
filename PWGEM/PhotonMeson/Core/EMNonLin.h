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

#include <Framework/Configurable.h>

#include <cstdint> // uint8_t

namespace o2::pwgem::nonlin
{

/// \class EMNonLin
/// \brief Class to obtain non linear correction factors for PbPb.
class EMNonLin
{
 public:
  static constexpr int MaxIter = 2;
  static constexpr int PhotonN = 3;
  static constexpr int CentBins = 10;

  enum class PhotonType : uint8_t {
    kEMC = 0,
    kPCM = 1,
    kPHOS = 2 // just in case
  };

  struct NonLinParams {
    float par0{0.f};
    float par1{0.f};
    float par2{0.f};
  };

  struct Context {
    const NonLinParams* params = nullptr;
    int nIter = 0;

    void setParams(const NonLinParams* newParams)
    {
      params = newParams;
    }

    void setIter(int iter)
    {
      if (iter < 0 || iter >= MaxIter) {
        nIter = MaxIter - 1;
        return;
      }

      nIter = iter;
    }
  };

  /// \brief gets the correction value for energy or pT for a specific
  /// \param inputCalibValue pT or energy of the photon that needs calibration
  /// \param ctx Context which has the centrality, photontype and number of iterations stored inside
  static float getCorrectionFactor(float inputCalibValue, const Context& ctx);

  /// \brief sets the parameters accordingly to the photon type, centrality and the wanted iteration level
  /// \param photonType type of the photon (e.g. 0 for EMC)
  /// \param cent centrality of the current collision in case the correction is centrality dependent
  static const NonLinParams* resolveParams(PhotonType type, float cent);

  /// \brief sets the parameters accordingly to the photon type, centrality and the wanted iteration level for MC
  /// \param photonType type of the photon (e.g. 0 for EMC)
  /// \param cent centrality of the current collision in case the correction is centrality dependent
  static const NonLinParams* resolveParamsMC(PhotonType type, float cent);

 private:
  static constexpr NonLinParams kNonLinTable
    [PhotonN][CentBins][MaxIter] =
      {
        // ============================
        // PhotonType::kEMC  (0)
        // ============================
        {
          // 00–10
          {
            {-5.33426e-01f, 1.40144e-02f, -5.24434e-01f}, // iter 0
            {0.f, 0.f, 0.f}                               // iter 1
          },
          // 10–20
          {
            {-5.33426e-01f, 1.40144e-02f, -5.24434e-01f},
            {0.f, 0.f, 0.f}},
          // 20–30
          {
            {-5.33426e-01f, 1.40144e-02f, -5.24434e-01f},
            {0.f, 0.f, 0.f}},
          // 30–40
          {
            {-5.33426e-01f, 1.40144e-02f, -5.24434e-01f},
            {0.f, 0.f, 0.f}},
          // 40–50
          {
            {-5.33426e-01f, 1.40144e-02f, -5.24434e-01f},
            {0.f, 0.f, 0.f}},
          // 50–60
          {
            {-5.33426e-01f, 1.40144e-02f, -5.24434e-01f},
            {0.f, 0.f, 0.f}},
          // 60–70
          {
            {-5.33426e-01f, 1.40144e-02f, -5.24434e-01f},
            {0.f, 0.f, 0.f}},
          // 70–80
          {
            {-5.33426e-01f, 1.40144e-02f, -5.24434e-01f},
            {0.f, 0.f, 0.f}},
          // 80–90
          {
            {-5.33426e-01f, 1.40144e-02f, -5.24434e-01f},
            {0.f, 0.f, 0.f}},
          // 90–100
          {
            {-5.33426e-01f, 1.40144e-02f, -5.24434e-01f},
            {0.f, 0.f, 0.f}}},

        // ============================
        // PhotonType::kPCM  (1)
        // ============================
        {
          // 00–10
          {
            {10.7203f, 0.0383968f, 10.6025f}, // iter 0
            {7.84549f, 0.0250021f, 7.86976f}  // iter 1
          },
          // 10–20
          {
            {10.7203f, 0.0383968f, 10.6025f},
            {7.84549f, 0.0250021f, 7.86976f}},
          // 20–30
          {
            {10.7203f, 0.0383968f, 10.6025f},
            {7.84549f, 0.0250021f, 7.86976f}},
          // 30–40
          {
            {10.7203f, 0.0383968f, 10.6025f},
            {7.84549f, 0.0250021f, 7.86976f}},
          // 40–50
          {
            {10.7203f, 0.0383968f, 10.6025f},
            {7.84549f, 0.0250021f, 7.86976f}},
          // 50–60
          {
            {10.7203f, 0.0383968f, 10.6025f},
            {7.84549f, 0.0250021f, 7.86976f}},
          // 60–70
          {
            {10.7203f, 0.0383968f, 10.6025f},
            {7.84549f, 0.0250021f, 7.86976f}},
          // 70–80
          {
            {10.7203f, 0.0383968f, 10.6025f},
            {7.84549f, 0.0250021f, 7.86976f}},
          // 80–90
          {
            {10.7203f, 0.0383968f, 10.6025f},
            {7.84549f, 0.0250021f, 7.86976f}},
          // 90–100
          {
            {10.7203f, 0.0383968f, 10.6025f},
            {7.84549f, 0.0250021f, 7.86976f}}},

        // ============================
        // PhotonType::kPHOS  (2)
        // ============================
        {
          // All centralities identical
          {{0.f, 0.f, 0.f}, {0.f, 0.f, 0.f}},
          {{0.f, 0.f, 0.f}, {0.f, 0.f, 0.f}},
          {{0.f, 0.f, 0.f}, {0.f, 0.f, 0.f}},
          {{0.f, 0.f, 0.f}, {0.f, 0.f, 0.f}},
          {{0.f, 0.f, 0.f}, {0.f, 0.f, 0.f}},
          {{0.f, 0.f, 0.f}, {0.f, 0.f, 0.f}},
          {{0.f, 0.f, 0.f}, {0.f, 0.f, 0.f}},
          {{0.f, 0.f, 0.f}, {0.f, 0.f, 0.f}},
          {{0.f, 0.f, 0.f}, {0.f, 0.f, 0.f}},
          {{0.f, 0.f, 0.f}, {0.f, 0.f, 0.f}}}};

  static constexpr NonLinParams kNonLinTableMC
    [PhotonN][CentBins][MaxIter] =
      {
        // ============================
        // PhotonType::kEMC  (0)
        // ============================
        {
          // 00–10
          {
            {-5.33426e-01f, 1.40144e-02f, -5.24434e-01f}, // iter 0
            {0.f, 0.f, 0.f}                               // iter 1
          },
          // 10–20
          {
            {-5.33426e-01f, 1.40144e-02f, -5.24434e-01f},
            {0.f, 0.f, 0.f}},
          // 20–30
          {
            {-5.33426e-01f, 1.40144e-02f, -5.24434e-01f},
            {0.f, 0.f, 0.f}},
          // 30–40
          {
            {-5.33426e-01f, 1.40144e-02f, -5.24434e-01f},
            {0.f, 0.f, 0.f}},
          // 40–50
          {
            {-5.33426e-01f, 1.40144e-02f, -5.24434e-01f},
            {0.f, 0.f, 0.f}},
          // 50–60
          {
            {-5.33426e-01f, 1.40144e-02f, -5.24434e-01f},
            {0.f, 0.f, 0.f}},
          // 60–70
          {
            {-5.33426e-01f, 1.40144e-02f, -5.24434e-01f},
            {0.f, 0.f, 0.f}},
          // 70–80
          {
            {-5.33426e-01f, 1.40144e-02f, -5.24434e-01f},
            {0.f, 0.f, 0.f}},
          // 80–90
          {
            {-5.33426e-01f, 1.40144e-02f, -5.24434e-01f},
            {0.f, 0.f, 0.f}},
          // 90–100
          {
            {-5.33426e-01f, 1.40144e-02f, -5.24434e-01f},
            {0.f, 0.f, 0.f}}},

        // ============================
        // PhotonType::kPCM  (1)
        // ============================
        {
          // 00–10
          {
            {10.7203f, 0.0383968f, 10.6025f}, // iter 0
            {7.84549f, 0.0250021f, 7.86976f}  // iter 1
          },
          // 10–20
          {
            {10.7203f, 0.0383968f, 10.6025f},
            {7.84549f, 0.0250021f, 7.86976f}},
          // 20–30
          {
            {10.7203f, 0.0383968f, 10.6025f},
            {7.84549f, 0.0250021f, 7.86976f}},
          // 30–40
          {
            {10.7203f, 0.0383968f, 10.6025f},
            {7.84549f, 0.0250021f, 7.86976f}},
          // 40–50
          {
            {10.7203f, 0.0383968f, 10.6025f},
            {7.84549f, 0.0250021f, 7.86976f}},
          // 50–60
          {
            {10.7203f, 0.0383968f, 10.6025f},
            {7.84549f, 0.0250021f, 7.86976f}},
          // 60–70
          {
            {10.7203f, 0.0383968f, 10.6025f},
            {7.84549f, 0.0250021f, 7.86976f}},
          // 70–80
          {
            {10.7203f, 0.0383968f, 10.6025f},
            {7.84549f, 0.0250021f, 7.86976f}},
          // 80–90
          {
            {10.7203f, 0.0383968f, 10.6025f},
            {7.84549f, 0.0250021f, 7.86976f}},
          // 90–100
          {
            {10.7203f, 0.0383968f, 10.6025f},
            {7.84549f, 0.0250021f, 7.86976f}}},

        // ============================
        // PhotonType::kPHOS  (2)
        // ============================
        {
          // All centralities identical
          {{0.f, 0.f, 0.f}, {0.f, 0.f, 0.f}},
          {{0.f, 0.f, 0.f}, {0.f, 0.f, 0.f}},
          {{0.f, 0.f, 0.f}, {0.f, 0.f, 0.f}},
          {{0.f, 0.f, 0.f}, {0.f, 0.f, 0.f}},
          {{0.f, 0.f, 0.f}, {0.f, 0.f, 0.f}},
          {{0.f, 0.f, 0.f}, {0.f, 0.f, 0.f}},
          {{0.f, 0.f, 0.f}, {0.f, 0.f, 0.f}},
          {{0.f, 0.f, 0.f}, {0.f, 0.f, 0.f}},
          {{0.f, 0.f, 0.f}, {0.f, 0.f, 0.f}},
          {{0.f, 0.f, 0.f}, {0.f, 0.f, 0.f}}}};
};
} // namespace o2::pwgem::nonlin

#endif // PWGEM_PHOTONMESON_CORE_EMNONLIN_H_
