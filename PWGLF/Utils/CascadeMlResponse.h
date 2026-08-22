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

/// \file CascadeMlResponse.h
/// \brief Class to compute the ML response for Strangeness-analysis selections
/// \author Gianni Shigeru Setoue Liveraro Catalano <gianni.shigeru.setoue.liveraro@cern.ch>, UNICAMP
/// \author Romain Schotter <romain.schotter@cern.ch>, Austrian Academy of Sciences
/// \author David Dobrigkeit Chinellato <david.dobrigkeit.chinellato@cern.ch>, Austrian Academy of Sciences

#ifndef PWGLF_UTILS_CASCADEMLRESPONSE_H_
#define PWGLF_UTILS_CASCADEMLRESPONSE_H_

#include "Tools/ML/MlResponse.h"

#include <vector>

namespace o2::analysis
{
// list of input features that can be requested via the mlConfigurations.namesInputFeatures configurable
enum class InputFeaturesCasc : uint8_t {
  cascradius = 0,
  v0radius,
  casccosPA,
  v0cosPA,
  dcapostopv,
  dcanegtopv,
  dcabachtopv,
  dcacascdaughters,
  dcaV0daughters,
  dcav0topv,
  bachBaryonCosPA,
  bachBaryonDCAxyToPV
};

template <typename TypeOutputScore = float>
class CascadeMlResponse : public MlResponse<TypeOutputScore>
{
 public:
  CascadeMlResponse() = default;
  ~CascadeMlResponse() override = default;

  template <typename TCascade, typename TCollision>
  std::vector<float> getInputFeatures(TCascade const& casc, TCollision const& coll)
  {
    std::vector<float> inputFeatures;
    inputFeatures.reserve(MlResponse<TypeOutputScore>::mCachedIndices.size());

    for (const auto& idx : MlResponse<TypeOutputScore>::mCachedIndices) {
      switch (static_cast<InputFeaturesCasc>(idx)) {
        case InputFeaturesCasc::cascradius:
          inputFeatures.emplace_back(casc.cascradius());
          break;
        case InputFeaturesCasc::v0radius:
          inputFeatures.emplace_back(casc.v0radius());
          break;
        case InputFeaturesCasc::casccosPA:
          inputFeatures.emplace_back(casc.casccosPA(coll.posX(), coll.posY(), coll.posZ()));
          break;
        case InputFeaturesCasc::v0cosPA:
          inputFeatures.emplace_back(casc.v0cosPA(coll.posX(), coll.posY(), coll.posZ()));
          break;
        case InputFeaturesCasc::dcapostopv:
          inputFeatures.emplace_back(casc.dcapostopv());
          break;
        case InputFeaturesCasc::dcanegtopv:
          inputFeatures.emplace_back(casc.dcanegtopv());
          break;
        case InputFeaturesCasc::dcabachtopv:
          inputFeatures.emplace_back(casc.dcabachtopv());
          break;
        case InputFeaturesCasc::dcacascdaughters:
          inputFeatures.emplace_back(casc.dcacascdaughters());
          break;
        case InputFeaturesCasc::dcaV0daughters:
          inputFeatures.emplace_back(casc.dcaV0daughters());
          break;
        case InputFeaturesCasc::dcav0topv:
          inputFeatures.emplace_back(casc.dcav0topv(coll.posX(), coll.posY(), coll.posZ()));
          break;
        case InputFeaturesCasc::bachBaryonCosPA:
          inputFeatures.emplace_back(casc.bachBaryonCosPA());
          break;
        case InputFeaturesCasc::bachBaryonDCAxyToPV:
          inputFeatures.emplace_back(casc.bachBaryonDCAxyToPV());
          break;
      }
    }
    return inputFeatures;
  }

 protected:
  void setAvailableInputFeatures() override
  {
    MlResponse<TypeOutputScore>::mAvailableInputFeatures = {
      {"cascradius", static_cast<uint8_t>(InputFeaturesCasc::cascradius)},
      {"v0radius", static_cast<uint8_t>(InputFeaturesCasc::v0radius)},
      {"casccosPA", static_cast<uint8_t>(InputFeaturesCasc::casccosPA)},
      {"v0cosPA", static_cast<uint8_t>(InputFeaturesCasc::v0cosPA)},
      {"dcapostopv", static_cast<uint8_t>(InputFeaturesCasc::dcapostopv)},
      {"dcanegtopv", static_cast<uint8_t>(InputFeaturesCasc::dcanegtopv)},
      {"dcabachtopv", static_cast<uint8_t>(InputFeaturesCasc::dcabachtopv)},
      {"dcacascdaughters", static_cast<uint8_t>(InputFeaturesCasc::dcacascdaughters)},
      {"dcaV0daughters", static_cast<uint8_t>(InputFeaturesCasc::dcaV0daughters)},
      {"dcav0topv", static_cast<uint8_t>(InputFeaturesCasc::dcav0topv)},
      {"bachBaryonCosPA", static_cast<uint8_t>(InputFeaturesCasc::bachBaryonCosPA)},
      {"bachBaryonDCAxyToPV", static_cast<uint8_t>(InputFeaturesCasc::bachBaryonDCAxyToPV)}};
  }
};

} // namespace o2::analysis

#endif // PWGLF_UTILS_CASCADEMLRESPONSE_H_
