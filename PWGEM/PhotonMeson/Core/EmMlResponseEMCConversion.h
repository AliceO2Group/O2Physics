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

/// \file EmMlResponseEMCConversion.h
/// \brief Class to compute the ML response for EMC conversion selections
/// \author Marvin Hemmer <marvin.hemmer@cern.ch>

#ifndef PWGEM_PHOTONMESON_CORE_EMMLRESPONSEEMCCONVERSION_H_
#define PWGEM_PHOTONMESON_CORE_EMMLRESPONSEEMCCONVERSION_H_

#include "PWGEM/PhotonMeson/Core/EMCConversionCandidate.h"

#include "Tools/ML/MlResponse.h"

#include <cstdint>
#include <vector>

// Fill the map of available input features
// the key is the feature's name (std::string)
// the value is the corresponding value in EnumInputFeatures
#define FILL_MAP_EMC_CONV(FEATURE) \
  {                                \
    #FEATURE, static_cast<uint8_t>(InputFeaturesEMCConversion::FEATURE)}

// Check if the index of mCachedIndices (index associated to a FEATURE)
// matches the entry in EnumInputFeatures associated to this FEATURE
// if so, the inputFeatures vector is filled with the FEATURE's value
// by calling the corresponding GETTER from OBJECT
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define CHECK_AND_FILL_VEC_EMC_CONV(GETTER)                        \
  case static_cast<uint8_t>(InputFeaturesEMCConversion::GETTER): { \
    inputFeatures.emplace_back(candidate.GETTER());                \
    break;                                                         \
  }

namespace o2::analysis::em::emcconv
{

// input feature used in the ml model
enum class InputFeaturesEMCConversion : uint8_t {
  minv,
  deltaEta,
  deltaR,
  phiv,
  rConv,
  totE,
  e2,
  e1,
  deltaPhi
};

template <typename TypeOutputScore = float>
class EmMlResponseEMCConversion : public MlResponse<TypeOutputScore>
{
 public:
  EmMlResponseEMCConversion() = default;
  virtual ~EmMlResponseEMCConversion() = default;

  template <o2::analysis::em::IsEmcConversionCandidate TCandidate>
  std::vector<float> getInputFeatures(TCandidate const& candidate)
  {
    std::vector<float> inputFeatures;
    for (const auto& idx : MlResponse<TypeOutputScore>::mCachedIndices) {
      switch (idx) {
        CHECK_AND_FILL_VEC_EMC_CONV(minv)
        CHECK_AND_FILL_VEC_EMC_CONV(deltaEta)
        CHECK_AND_FILL_VEC_EMC_CONV(deltaR)
        CHECK_AND_FILL_VEC_EMC_CONV(phiv)
        CHECK_AND_FILL_VEC_EMC_CONV(rConv)
        CHECK_AND_FILL_VEC_EMC_CONV(totE)
        CHECK_AND_FILL_VEC_EMC_CONV(e2)
        CHECK_AND_FILL_VEC_EMC_CONV(e1)
        CHECK_AND_FILL_VEC_EMC_CONV(deltaPhi)
      }
    }
    return inputFeatures;
  }

 protected:
  void setAvailableInputFeatures()
  {
    MlResponse<TypeOutputScore>::mAvailableInputFeatures = {
      FILL_MAP_EMC_CONV(minv), FILL_MAP_EMC_CONV(deltaEta), FILL_MAP_EMC_CONV(deltaR),
      FILL_MAP_EMC_CONV(phiv), FILL_MAP_EMC_CONV(rConv), FILL_MAP_EMC_CONV(totE),
      FILL_MAP_EMC_CONV(e2), FILL_MAP_EMC_CONV(e1), FILL_MAP_EMC_CONV(deltaPhi)};
  }
};

} // namespace o2::analysis::em::emcconv

#undef FILL_MAP_EMC_CONV
#undef CHECK_AND_FILL_VEC_EMC_CONV

#endif // PWGEM_PHOTONMESON_CORE_EMMLRESPONSEEMCCONVERSION_H_
