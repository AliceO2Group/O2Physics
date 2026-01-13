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

/// \file EmMLResponsePCM.h
/// \brief Class to compute the ML response for PCM analysis selections
/// \author Isabel Kantak <isabel.kantak@cern.ch>, University of Heidelberg

#ifndef PWGEM_PHOTONMESON_CORE_EMMLRESPONSEPCM_H_
#define PWGEM_PHOTONMESON_CORE_EMMLRESPONSEPCM_H_

#include "PWGEM/PhotonMeson/Core/EmMlResponse.h"

#include "Tools/ML/MlResponse.h"

#include <cstdint>
#include <vector>

// Fill the map of available input features
// the key is the feature's name (std::string)
// the value is the corresponding value in EnumInputFeatures
#define FILL_MAP_PCM(FEATURE) \
  {                           \
    #FEATURE, static_cast<uint8_t>(InputFeaturesPCM::FEATURE)}

// Check if the index of mCachedIndices (index associated to a FEATURE)
// matches the entry in EnumInputFeatures associated to this FEATURE
// if so, the inputFeatures vector is filled with the FEATURE's value
// by calling the corresponding GETTER from OBJECT
#define CHECK_AND_FILL_VEC_PCM_FULL(OBJECT, FEATURE, GETTER) \
  case static_cast<uint8_t>(InputFeaturesPCM::FEATURE): {    \
    inputFeatures.emplace_back(OBJECT.GETTER());             \
    break;                                                   \
  }

// Specific case of CHECK_AND_FILL_VEC_PCM_FULL(OBJECT, FEATURE, GETTER)
// where OBJECT is named candidate and FEATURE = GETTER
#define CHECK_AND_FILL_VEC_PCM(GETTER)                   \
  case static_cast<uint8_t>(InputFeaturesPCM::GETTER): { \
    inputFeatures.emplace_back(candidate.GETTER());      \
    break;                                               \
  }

namespace o2::analysis
{

enum class InputFeaturesPCM : uint8_t {
  V0PhotonCandidate_fDCAxyToPV,
  V0PhotonCandidate_fDCAzToPV,
  V0PhotonCandidate_fPCA,
  V0PhotonCandidate_fAlpha,
  V0PhotonCandidate_fQtArm,
  V0PhotonCandidate_fChiSquareNDF,
  V0PhotonCandidate_fCosPA,
  Pos_V0Leg_fTPCNSigmaEl,
  Pos_V0Leg_fTPCNSigmaPi,
  Neg_V0Leg_fTPCNSigmaEl,
  Neg_V0Leg_fTPCNSigmaPi
};

template <typename TypeOutputScore = float>
class EmMlResponsePCM : public EmMlResponse<TypeOutputScore>
{
 public:
  /// Default constructor
  EmMlResponsePCM() = default;
  /// Default destructor
  virtual ~EmMlResponsePCM() = default;

  /// Method to get the input features vector needed for ML inference
  /// \param candidate is the V0photon candidate
  /// \return inputFeatures vector
  template <typename T1, typename T2>
  std::vector<float> getInputFeatures(T1 const& candidate, T2 const& posLeg, T2 const& negLeg)
  {
    std::vector<float> inputFeatures;

    for (const auto& idx : MlResponse<TypeOutputScore>::mCachedIndices) {
      switch (idx) {
        CHECK_AND_FILL_VEC_PCM_FULL(candidate, V0PhotonCandidate_fDCAxyToPV, GetDcaXYToPV);
        CHECK_AND_FILL_VEC_PCM_FULL(candidate, V0PhotonCandidate_fDCAzToPV, GetDcaZToPV);
        CHECK_AND_FILL_VEC_PCM_FULL(candidate, V0PhotonCandidate_fPCA, GetPCA);
        CHECK_AND_FILL_VEC_PCM_FULL(candidate, V0PhotonCandidate_fAlpha, GetAlpha);
        CHECK_AND_FILL_VEC_PCM_FULL(candidate, V0PhotonCandidate_fQtArm, GetQt);
        CHECK_AND_FILL_VEC_PCM_FULL(candidate, V0PhotonCandidate_fChiSquareNDF, GetChi2NDF);
        CHECK_AND_FILL_VEC_PCM_FULL(candidate, V0PhotonCandidate_fCosPA, GetCosPA);
        CHECK_AND_FILL_VEC_PCM_FULL(posLeg, Pos_V0Leg_fTPCNSigmaEl, tpcNSigmaEl);
        CHECK_AND_FILL_VEC_PCM_FULL(posLeg, Pos_V0Leg_fTPCNSigmaPi, tpcNSigmaPi);
        CHECK_AND_FILL_VEC_PCM_FULL(negLeg, Neg_V0Leg_fTPCNSigmaEl, tpcNSigmaEl);
        CHECK_AND_FILL_VEC_PCM_FULL(negLeg, Neg_V0Leg_fTPCNSigmaPi, tpcNSigmaPi);
      }
    }
    return inputFeatures;
  }

 protected:
  /// Method to fill the map of available input features
  void setAvailableInputFeatures()
  {
    MlResponse<TypeOutputScore>::mAvailableInputFeatures = {
      FILL_MAP_PCM(V0PhotonCandidate_fDCAxyToPV),
      FILL_MAP_PCM(V0PhotonCandidate_fDCAzToPV),
      FILL_MAP_PCM(V0PhotonCandidate_fPCA),
      FILL_MAP_PCM(V0PhotonCandidate_fAlpha),
      FILL_MAP_PCM(V0PhotonCandidate_fQtArm),
      FILL_MAP_PCM(V0PhotonCandidate_fChiSquareNDF),
      FILL_MAP_PCM(V0PhotonCandidate_fCosPA),
      FILL_MAP_PCM(Pos_V0Leg_fTPCNSigmaEl),
      FILL_MAP_PCM(Pos_V0Leg_fTPCNSigmaPi),
      FILL_MAP_PCM(Neg_V0Leg_fTPCNSigmaEl),
      FILL_MAP_PCM(Neg_V0Leg_fTPCNSigmaPi)};
  }
};

} // namespace o2::analysis

#undef FILL_MAP_PCM
#undef CHECK_AND_FILL_VEC_PCM_FULL
#undef CHECK_AND_FILL_VEC_PCM

#endif // PWGEM_PHOTONMESON_CORE_EMMLRESPONSEPCM_H_
