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

/// \file MulticharmMlResponse.h
/// \brief Class to compute the ML response for multi-charm candidates
/// \author Jesper Karlsson Gumprecht

#ifndef ALICE3_ML_MULTICHARMMLRESPONSE_H_
#define ALICE3_ML_MULTICHARMMLRESPONSE_H_

#include "Tools/ML/MlResponse.h"

#include <array>
#include <cstdint>
#include <string>
#include <vector>

namespace multi_charm_ml
{

static constexpr int NBinsPt = 10;
static constexpr int NClasses = 2;
static constexpr std::array<std::array<double, NClasses>, NBinsPt> Cuts = {{{0.5, 0.5},
                                                                            {0.5, 0.5},
                                                                            {0.5, 0.5},
                                                                            {0.5, 0.5},
                                                                            {0.5, 0.5},
                                                                            {0.5, 0.5},
                                                                            {0.5, 0.5},
                                                                            {0.5, 0.5},
                                                                            {0.5, 0.5},
                                                                            {0.5, 0.5}}};

static const std::vector<std::string> labelsPt = {"pT bin 0",
                                                  "pT bin 1",
                                                  "pT bin 2",
                                                  "pT bin 3",
                                                  "pT bin 4",
                                                  "pT bin 5",
                                                  "pT bin 6",
                                                  "pT bin 7",
                                                  "pT bin 8",
                                                  "pT bin 9"};

static const std::vector<std::string> labelsCutScore = {"score class 1", "score class 2"};
static const std::vector<std::string> namesInputFeatures{"xicDauDCA",
                                                         "xiccDauDCA",
                                                         "xiDCAxy",
                                                         "xicDCAxy",
                                                         "xiccDCAxy",
                                                         "xiDCAz",
                                                         "xicDCAz",
                                                         "xiccDCAz",
                                                         "pi1cDCAxy",
                                                         "pi2cDCAxy",
                                                         "piccDCAxy",
                                                         "pi1cDCAz",
                                                         "pi2cDCAz",
                                                         "piccDCAz",
                                                         "xicDecayRadius2D",
                                                         "xiccDecayRadius2D",
                                                         "xicProperLength",
                                                         "xicDistanceFromPV",
                                                         "xiccProperLength"};

} // namespace multi_charm_ml

namespace o2::analysis
{
// list of input features that can be requested via the mlConfigurations.namesInputFeatures configurable
enum class InputFeaturesMulticharm : uint8_t {
  xicDauDCA = 0,
  xiccDauDCA,
  xiDCAxy,
  xicDCAxy,
  xiccDCAxy,
  xiDCAz,
  xicDCAz,
  xiccDCAz,
  pi1cDCAxy,
  pi2cDCAxy,
  piccDCAxy,
  pi1cDCAz,
  pi2cDCAz,
  piccDCAz,
  xicDecayRadius2D,
  xiccDecayRadius2D,
  xicProperLength,
  xicDistanceFromPV,
  xiccProperLength
};

template <typename TypeOutputScore = float>
class MulticharmMlResponse : public MlResponse<TypeOutputScore>
{
 public:
  MulticharmMlResponse() = default;
  ~MulticharmMlResponse() override = default;

  template <typename TMulticharm, typename TCollision>
  std::vector<float> getInputFeatures(TMulticharm const& multicharm)
  {
    std::vector<float> inputFeatures;
    inputFeatures.reserve(MlResponse<TypeOutputScore>::mCachedIndices.size());

    for (const auto& idx : MlResponse<TypeOutputScore>::mCachedIndices) {
      switch (idx) {
        case InputFeaturesMulticharm::xicDauDCA:
          inputFeatures.emplace_back(multicharm.xicDauDCA());
          break;
        case InputFeaturesMulticharm::xiccDauDCA:
          inputFeatures.emplace_back(multicharm.xiccDauDCA());
          break;
        case InputFeaturesMulticharm::xiDCAxy:
          inputFeatures.emplace_back(multicharm.xiDCAxy());
          break;
        case InputFeaturesMulticharm::xicDCAxy:
          inputFeatures.emplace_back(multicharm.xicDCAxy());
          break;
        case InputFeaturesMulticharm::xiccDCAxy:
          inputFeatures.emplace_back(multicharm.xiccDCAxy());
          break;
        case InputFeaturesMulticharm::xiDCAz:
          inputFeatures.emplace_back(multicharm.xiDCAz());
          break;
        case InputFeaturesMulticharm::xicDCAz:
          inputFeatures.emplace_back(multicharm.xicDCAz());
          break;
        case InputFeaturesMulticharm::xiccDCAz:
          inputFeatures.emplace_back(multicharm.xiccDCAz());
          break;
        case InputFeaturesMulticharm::pi1cDCAxy:
          inputFeatures.emplace_back(multicharm.pi1cDCAxy());
          break;
        case InputFeaturesMulticharm::pi2cDCAxy:
          inputFeatures.emplace_back(multicharm.pi2cDCAxy());
          break;
        case InputFeaturesMulticharm::piccDCAxy:
          inputFeatures.emplace_back(multicharm.piccDCAxy());
          break;
        case InputFeaturesMulticharm::pi1cDCAz:
          inputFeatures.emplace_back(multicharm.pi1cDCAz());
          break;
        case InputFeaturesMulticharm::pi2cDCAz:
          inputFeatures.emplace_back(multicharm.pi2cDCAz());
          break;
        case InputFeaturesMulticharm::piccDCAz:
          inputFeatures.emplace_back(multicharm.piccDCAz());
          break;
        case InputFeaturesMulticharm::xicDecayRadius2D:
          inputFeatures.emplace_back(multicharm.xicDecayRadius2D());
          break;
        case InputFeaturesMulticharm::xiccDecayRadius2D:
          inputFeatures.emplace_back(multicharm.xiccDecayRadius2D());
          break;
        case InputFeaturesMulticharm::xicProperLength:
          inputFeatures.emplace_back(multicharm.xicProperLength());
          break;
        case InputFeaturesMulticharm::xicDistanceFromPV:
          inputFeatures.emplace_back(multicharm.xicDistanceFromPV());
          break;
        case InputFeaturesMulticharm::xiccProperLength:
          inputFeatures.emplace_back(multicharm.xiccProperLength());
          break;
      }
    }
    return inputFeatures;
  }

 protected:
  void setAvailableInputFeatures() override
  {
    MlResponse<TypeOutputScore>::mAvailableInputFeatures = {
      {"xicDauDCA", static_cast<uint8_t>(InputFeaturesMulticharm::xicDauDCA)},
      {"xiccDauDCA", static_cast<uint8_t>(InputFeaturesMulticharm::xiccDauDCA)},
      {"xiDCAxy", static_cast<uint8_t>(InputFeaturesMulticharm::xiDCAxy)},
      {"xicDCAxy", static_cast<uint8_t>(InputFeaturesMulticharm::xicDCAxy)},
      {"xiccDCAxy", static_cast<uint8_t>(InputFeaturesMulticharm::xiccDCAxy)},
      {"xiDCAz", static_cast<uint8_t>(InputFeaturesMulticharm::xiDCAz)},
      {"xicDCAz", static_cast<uint8_t>(InputFeaturesMulticharm::xicDCAz)},
      {"xiccDCAz", static_cast<uint8_t>(InputFeaturesMulticharm::xiccDCAz)},
      {"pi1cDCAxy", static_cast<uint8_t>(InputFeaturesMulticharm::pi1cDCAxy)},
      {"pi2cDCAxy", static_cast<uint8_t>(InputFeaturesMulticharm::pi2cDCAxy)},
      {"piccDCAxy", static_cast<uint8_t>(InputFeaturesMulticharm::piccDCAxy)},
      {"pi1cDCAz", static_cast<uint8_t>(InputFeaturesMulticharm::pi1cDCAz)},
      {"pi2cDCAz", static_cast<uint8_t>(InputFeaturesMulticharm::pi2cDCAz)},
      {"piccDCAz", static_cast<uint8_t>(InputFeaturesMulticharm::piccDCAz)},
      {"xicDecayRadius2D", static_cast<uint8_t>(InputFeaturesMulticharm::xicDecayRadius2D)},
      {"xiccDecayRadius2D", static_cast<uint8_t>(InputFeaturesMulticharm::xiccDecayRadius2D)},
      {"xicProperLength", static_cast<uint8_t>(InputFeaturesMulticharm::xicProperLength)},
      {"xicDistanceFromPV", static_cast<uint8_t>(InputFeaturesMulticharm::xicDistanceFromPV)},
      {"xiccProperLength", static_cast<uint8_t>(InputFeaturesMulticharm::xiccProperLength)}};
  }
};

} // namespace o2::analysis

#endif // ALICE3_ML_MULTICHARMMLRESPONSE_H_
