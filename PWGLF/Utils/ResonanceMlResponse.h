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

/// \file ResonanceMlResponse.h
/// \brief Class to compute the ML response for Lambda + photon resonance selections (Sigma0, Lambda(1520)),
///        either from sigma0builder candidates or from a Lambda and a photon V0 of the strangeness derived data
/// \author Oussama Benchikhi

#ifndef PWGLF_UTILS_RESONANCEMLRESPONSE_H_
#define PWGLF_UTILS_RESONANCEMLRESPONSE_H_

#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Tools/ML/MlResponse.h"

#include <Framework/ASoA.h>

#include <cstdint>
#include <vector>

namespace o2::analysis
{
// list of input features that can be requested via the namesInputFeatures configurable (sigma0builder getter names)
enum class InputFeaturesResonance : uint8_t {
  // Lambda
  lambdaQt = 0,
  lambdaAlpha,
  lambdaRadius,
  lambdaCosPA,
  lambdaDCADau,
  lambdaDCANegPV,
  lambdaDCAPosPV,
  lambdaPosEta,
  lambdaNegEta,
  lambdaPosPrTPCNSigma,
  lambdaPosPiTPCNSigma,
  lambdaNegPrTPCNSigma,
  lambdaNegPiTPCNSigma,
  // Photon
  photonQt,
  photonAlpha,
  photonCosPA,
  photonDCADau,
  photonDCANegPV,
  photonDCAPosPV,
  photonRadius,
  photonZconv,
  photonPsiPair,
  photonPosEta,
  photonNegEta,
  photonPosTPCNSigmaEl,
  photonNegTPCNSigmaEl,
  // Photon-Lambda pair
  opAngle
};

template <typename TypeOutputScore = float>
class ResonanceMlResponse : public MlResponse<TypeOutputScore>
{
 public:
  ResonanceMlResponse() = default;
  ~ResonanceMlResponse() override = default;

  /// Input features of a photon-Lambda pair
  template <typename TDauTracks = o2::soa::Join<o2::aod::DauTrackExtras, o2::aod::DauTrackTPCPIDs>, typename TLambda, typename TPhoton>
  std::vector<float> getInputFeatures(TLambda const& lambda, TPhoton const& photon, float opAngle)
  {
    // sigma0builder candidates carry prefixed columns (lambdaQt(), photonQt(), ...), derived V0s the plain V0 getters
    constexpr bool LambdaFromSigma0 = requires(TLambda const& cand) { cand.lambdaQt(); };
    constexpr bool PhotonFromSigma0 = requires(TPhoton const& cand) { cand.photonQt(); };

    std::vector<float> inputFeatures;
    inputFeatures.reserve(MlResponse<TypeOutputScore>::mCachedIndices.size());

    for (const auto& idx : MlResponse<TypeOutputScore>::mCachedIndices) {
      switch (static_cast<InputFeaturesResonance>(idx)) {
        // Lambda
        case InputFeaturesResonance::lambdaQt:
          if constexpr (LambdaFromSigma0) {
            inputFeatures.emplace_back(lambda.lambdaQt());
          } else {
            inputFeatures.emplace_back(lambda.qtarm());
          }
          break;
        case InputFeaturesResonance::lambdaAlpha:
          if constexpr (LambdaFromSigma0) {
            inputFeatures.emplace_back(lambda.lambdaAlpha());
          } else {
            inputFeatures.emplace_back(lambda.alpha());
          }
          break;
        case InputFeaturesResonance::lambdaRadius:
          if constexpr (LambdaFromSigma0) {
            inputFeatures.emplace_back(lambda.lambdaRadius());
          } else {
            inputFeatures.emplace_back(lambda.v0radius());
          }
          break;
        case InputFeaturesResonance::lambdaCosPA:
          if constexpr (LambdaFromSigma0) {
            inputFeatures.emplace_back(lambda.lambdaCosPA());
          } else {
            inputFeatures.emplace_back(lambda.v0cosPA());
          }
          break;
        case InputFeaturesResonance::lambdaDCADau:
          if constexpr (LambdaFromSigma0) {
            inputFeatures.emplace_back(lambda.lambdaDCADau());
          } else {
            inputFeatures.emplace_back(lambda.dcaV0daughters());
          }
          break;
        case InputFeaturesResonance::lambdaDCANegPV:
          if constexpr (LambdaFromSigma0) {
            inputFeatures.emplace_back(lambda.lambdaDCANegPV());
          } else {
            inputFeatures.emplace_back(lambda.dcanegtopv());
          }
          break;
        case InputFeaturesResonance::lambdaDCAPosPV:
          if constexpr (LambdaFromSigma0) {
            inputFeatures.emplace_back(lambda.lambdaDCAPosPV());
          } else {
            inputFeatures.emplace_back(lambda.dcapostopv());
          }
          break;
        case InputFeaturesResonance::lambdaPosEta:
          if constexpr (LambdaFromSigma0) {
            inputFeatures.emplace_back(lambda.lambdaPosEta());
          } else {
            inputFeatures.emplace_back(lambda.positiveeta());
          }
          break;
        case InputFeaturesResonance::lambdaNegEta:
          if constexpr (LambdaFromSigma0) {
            inputFeatures.emplace_back(lambda.lambdaNegEta());
          } else {
            inputFeatures.emplace_back(lambda.negativeeta());
          }
          break;
        case InputFeaturesResonance::lambdaPosPrTPCNSigma:
          if constexpr (LambdaFromSigma0) {
            inputFeatures.emplace_back(lambda.lambdaPosPrTPCNSigma());
          } else {
            inputFeatures.emplace_back(lambda.template posTrackExtra_as<TDauTracks>().tpcNSigmaPr());
          }
          break;
        case InputFeaturesResonance::lambdaPosPiTPCNSigma:
          if constexpr (LambdaFromSigma0) {
            inputFeatures.emplace_back(lambda.lambdaPosPiTPCNSigma());
          } else {
            inputFeatures.emplace_back(lambda.template posTrackExtra_as<TDauTracks>().tpcNSigmaPi());
          }
          break;
        case InputFeaturesResonance::lambdaNegPrTPCNSigma:
          if constexpr (LambdaFromSigma0) {
            inputFeatures.emplace_back(lambda.lambdaNegPrTPCNSigma());
          } else {
            inputFeatures.emplace_back(lambda.template negTrackExtra_as<TDauTracks>().tpcNSigmaPr());
          }
          break;
        case InputFeaturesResonance::lambdaNegPiTPCNSigma:
          if constexpr (LambdaFromSigma0) {
            inputFeatures.emplace_back(lambda.lambdaNegPiTPCNSigma());
          } else {
            inputFeatures.emplace_back(lambda.template negTrackExtra_as<TDauTracks>().tpcNSigmaPi());
          }
          break;
        // Photon
        case InputFeaturesResonance::photonQt:
          if constexpr (PhotonFromSigma0) {
            inputFeatures.emplace_back(photon.photonQt());
          } else {
            inputFeatures.emplace_back(photon.qtarm());
          }
          break;
        case InputFeaturesResonance::photonAlpha:
          if constexpr (PhotonFromSigma0) {
            inputFeatures.emplace_back(photon.photonAlpha());
          } else {
            inputFeatures.emplace_back(photon.alpha());
          }
          break;
        case InputFeaturesResonance::photonCosPA:
          if constexpr (PhotonFromSigma0) {
            inputFeatures.emplace_back(photon.photonCosPA());
          } else {
            inputFeatures.emplace_back(photon.v0cosPA());
          }
          break;
        case InputFeaturesResonance::photonDCADau:
          if constexpr (PhotonFromSigma0) {
            inputFeatures.emplace_back(photon.photonDCADau());
          } else {
            inputFeatures.emplace_back(photon.dcaV0daughters());
          }
          break;
        case InputFeaturesResonance::photonDCANegPV:
          if constexpr (PhotonFromSigma0) {
            inputFeatures.emplace_back(photon.photonDCANegPV());
          } else {
            inputFeatures.emplace_back(photon.dcanegtopv());
          }
          break;
        case InputFeaturesResonance::photonDCAPosPV:
          if constexpr (PhotonFromSigma0) {
            inputFeatures.emplace_back(photon.photonDCAPosPV());
          } else {
            inputFeatures.emplace_back(photon.dcapostopv());
          }
          break;
        case InputFeaturesResonance::photonRadius:
          if constexpr (PhotonFromSigma0) {
            inputFeatures.emplace_back(photon.photonRadius());
          } else {
            inputFeatures.emplace_back(photon.v0radius());
          }
          break;
        case InputFeaturesResonance::photonZconv:
          if constexpr (PhotonFromSigma0) {
            inputFeatures.emplace_back(photon.photonZconv());
          } else {
            inputFeatures.emplace_back(photon.z());
          }
          break;
        case InputFeaturesResonance::photonPsiPair:
          if constexpr (PhotonFromSigma0) {
            inputFeatures.emplace_back(photon.photonPsiPair());
          } else {
            inputFeatures.emplace_back(photon.psipair());
          }
          break;
        case InputFeaturesResonance::photonPosEta:
          if constexpr (PhotonFromSigma0) {
            inputFeatures.emplace_back(photon.photonPosEta());
          } else {
            inputFeatures.emplace_back(photon.positiveeta());
          }
          break;
        case InputFeaturesResonance::photonNegEta:
          if constexpr (PhotonFromSigma0) {
            inputFeatures.emplace_back(photon.photonNegEta());
          } else {
            inputFeatures.emplace_back(photon.negativeeta());
          }
          break;
        case InputFeaturesResonance::photonPosTPCNSigmaEl:
          if constexpr (PhotonFromSigma0) {
            inputFeatures.emplace_back(photon.photonPosTPCNSigmaEl());
          } else {
            inputFeatures.emplace_back(photon.template posTrackExtra_as<TDauTracks>().tpcNSigmaEl());
          }
          break;
        case InputFeaturesResonance::photonNegTPCNSigmaEl:
          if constexpr (PhotonFromSigma0) {
            inputFeatures.emplace_back(photon.photonNegTPCNSigmaEl());
          } else {
            inputFeatures.emplace_back(photon.template negTrackExtra_as<TDauTracks>().tpcNSigmaEl());
          }
          break;
        // Photon-Lambda pair
        case InputFeaturesResonance::opAngle:
          inputFeatures.emplace_back(opAngle);
          break;
      }
    }
    return inputFeatures;
  }

 protected:
  /// Method to fill the map of available input features
  void setAvailableInputFeatures() override
  {
    MlResponse<TypeOutputScore>::mAvailableInputFeatures = {
      {"lambdaQt", static_cast<uint8_t>(InputFeaturesResonance::lambdaQt)},
      {"lambdaAlpha", static_cast<uint8_t>(InputFeaturesResonance::lambdaAlpha)},
      {"lambdaRadius", static_cast<uint8_t>(InputFeaturesResonance::lambdaRadius)},
      {"lambdaCosPA", static_cast<uint8_t>(InputFeaturesResonance::lambdaCosPA)},
      {"lambdaDCADau", static_cast<uint8_t>(InputFeaturesResonance::lambdaDCADau)},
      {"lambdaDCANegPV", static_cast<uint8_t>(InputFeaturesResonance::lambdaDCANegPV)},
      {"lambdaDCAPosPV", static_cast<uint8_t>(InputFeaturesResonance::lambdaDCAPosPV)},
      {"lambdaPosEta", static_cast<uint8_t>(InputFeaturesResonance::lambdaPosEta)},
      {"lambdaNegEta", static_cast<uint8_t>(InputFeaturesResonance::lambdaNegEta)},
      {"lambdaPosPrTPCNSigma", static_cast<uint8_t>(InputFeaturesResonance::lambdaPosPrTPCNSigma)},
      {"lambdaPosPiTPCNSigma", static_cast<uint8_t>(InputFeaturesResonance::lambdaPosPiTPCNSigma)},
      {"lambdaNegPrTPCNSigma", static_cast<uint8_t>(InputFeaturesResonance::lambdaNegPrTPCNSigma)},
      {"lambdaNegPiTPCNSigma", static_cast<uint8_t>(InputFeaturesResonance::lambdaNegPiTPCNSigma)},
      {"photonQt", static_cast<uint8_t>(InputFeaturesResonance::photonQt)},
      {"photonAlpha", static_cast<uint8_t>(InputFeaturesResonance::photonAlpha)},
      {"photonCosPA", static_cast<uint8_t>(InputFeaturesResonance::photonCosPA)},
      {"photonDCADau", static_cast<uint8_t>(InputFeaturesResonance::photonDCADau)},
      {"photonDCANegPV", static_cast<uint8_t>(InputFeaturesResonance::photonDCANegPV)},
      {"photonDCAPosPV", static_cast<uint8_t>(InputFeaturesResonance::photonDCAPosPV)},
      {"photonRadius", static_cast<uint8_t>(InputFeaturesResonance::photonRadius)},
      {"photonZconv", static_cast<uint8_t>(InputFeaturesResonance::photonZconv)},
      {"photonPsiPair", static_cast<uint8_t>(InputFeaturesResonance::photonPsiPair)},
      {"photonPosEta", static_cast<uint8_t>(InputFeaturesResonance::photonPosEta)},
      {"photonNegEta", static_cast<uint8_t>(InputFeaturesResonance::photonNegEta)},
      {"photonPosTPCNSigmaEl", static_cast<uint8_t>(InputFeaturesResonance::photonPosTPCNSigmaEl)},
      {"photonNegTPCNSigmaEl", static_cast<uint8_t>(InputFeaturesResonance::photonNegTPCNSigmaEl)},
      {"opAngle", static_cast<uint8_t>(InputFeaturesResonance::opAngle)}};
  }
};

} // namespace o2::analysis

#endif // PWGLF_UTILS_RESONANCEMLRESPONSE_H_
