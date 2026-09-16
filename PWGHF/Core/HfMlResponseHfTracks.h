// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file HfMlResponseHfTracks.h
/// \brief Class to compute the ML response for single-track selection of HF daughters
/// \author Fabrizio Chinu <fabrizio.chinu@cern.ch>, Universita and INFN Torino

#ifndef PWGHF_CORE_HFMLRESPONSEHFTRACKS_H_
#define PWGHF_CORE_HFMLRESPONSEHFTRACKS_H_

#include "PWGHF/Core/HfMlResponse.h"

#include "Tools/ML/MlResponse.h"

#include <algorithm>
#include <cstdint>
#include <iterator>
#include <vector>

// Fill the map of available input features
// the key is the feature's name (std::string)
// the value is the corresponding value in EnumInputFeatures
#define FILL_MAP_HF_TRACK(FEATURE) \
  {                                \
    #FEATURE, static_cast<uint8_t>(InputFeaturesTracks::FEATURE)}

// Check if the index of mCachedIndices (index associated to a FEATURE)
// matches the entry in EnumInputFeatures associated to this FEATURE
// if so, the inputFeatures vector is filled with the FEATURE's value
#define CHECK_AND_FILL_VEC_HF_TRACK(FEATURE)               \
  case static_cast<uint8_t>(InputFeaturesTracks::FEATURE): \
    inputFeatures.emplace_back(track.FEATURE);             \
    break;

namespace o2::analysis
{

/// Output classes of the track-level HF model.
/// The ordering must match the one used at training time.
enum class HfTrackMlClass : uint8_t {
  Background = 0, ///< track is not a daughter of a D+ → πKπ / Ds → φπ → KKπ decay
  Pion,           ///< track is the π± of a D+ → πKπ / Ds → φπ → KKπ decay
  Kaon,           ///< track is the K∓ of a D+ → πKπ / Ds → φπ → KKπ decay
  NClasses
};

enum class InputFeaturesTracks : uint8_t {
  // kinematics
  pt = 0,
  eta,
  // impact parameter and its resolution
  dcaXY,
  dcaZ,
  sigmaDcaXY,
  sigmaDcaZ,
  normDcaXY,
  normDcaZ,
  // track parameters
  signed1Pt,
  tgl,
  sign,
  isPvContributor,
  // ITS quality
  itsNCls,
  itsNClsInnerBarrel,
  itsChi2NCl,
  // TPC quality
  tpcNClsFound,
  tpcCrossedRowsOverFindableCls,
  tpcChi2NCl,
  tpcFractionSharedCls,
  // TPC PID
  tpcNSigmaPi,
  tpcNSigmaKa,
  // event
  centrality
};

/// Plain container with every quantity the model can be fed.
/// It is filled once per (track, collision) association by the task and then reduced to the
/// configured subset of features by getInputFeatures().
/// Everything is a float so that the training tree and the inference input are bit-identical.
struct HfTrackMlFeatures {
  // kinematics
  float pt{0.f};
  float eta{0.f};
  // impact parameter and its resolution
  float dcaXY{0.f};
  float dcaZ{0.f};
  float sigmaDcaXY{0.f};
  float sigmaDcaZ{0.f};
  float normDcaXY{0.f};
  float normDcaZ{0.f};
  // track parameters
  float signed1Pt{0.f};
  float tgl{0.f};
  float sign{0.f};
  float isPvContributor{0.f};
  // ITS quality
  float itsNCls{0.f};
  float itsNClsInnerBarrel{0.f};
  float itsChi2NCl{0.f};
  // TPC quality
  float tpcNClsFound{0.f};
  float tpcCrossedRowsOverFindableCls{0.f};
  float tpcChi2NCl{0.f};
  float tpcFractionSharedCls{0.f};
  // TPC PID
  float tpcNSigmaPi{0.f};
  float tpcNSigmaKa{0.f};
  // event
  float centrality{0.f};
};

template <typename TypeOutputScore = float>
class HfMlResponseHfTracks : public HfMlResponse<TypeOutputScore>
{
 public:
  /// Default constructor
  HfMlResponseHfTracks() = default;
  /// Default destructor
  virtual ~HfMlResponseHfTracks() = default;

  /// Index of the model to be used for a given value of the binning variable (the track pT).
  /// Needed because the task applies its own per-class thresholds instead of the single-decision
  /// logic of MlResponse::isSelectedMl, whose private findBin is not reachable from here.
  /// Follows the same convention as the base class: mBinsLimits stores the bin edges.
  /// \param value is the value of the binning variable
  /// \return index of the model to be used, -1 if the value is outside the configured range
  int getModelBin(float value) const
  {
    const auto& binsLimits = MlResponse<TypeOutputScore>::mBinsLimits;
    const auto valueDouble = static_cast<double>(value);
    if (binsLimits.empty() || valueDouble < binsLimits.front() || valueDouble >= binsLimits.back()) {
      return -1;
    }
    return std::distance(binsLimits.begin(), std::upper_bound(binsLimits.begin(), binsLimits.end(), valueDouble)) - 1;
  }

  /// Method to get the input features vector needed for ML inference
  /// \param track is the container with all the candidate track quantities
  /// \return inputFeatures vector, in the order configured via cacheInputFeaturesIndices
  std::vector<float> getInputFeatures(HfTrackMlFeatures const& track)
  {
    std::vector<float> inputFeatures;
    inputFeatures.reserve(MlResponse<TypeOutputScore>::mCachedIndices.size());

    for (const auto& idx : MlResponse<TypeOutputScore>::mCachedIndices) {
      switch (idx) {
        CHECK_AND_FILL_VEC_HF_TRACK(pt);
        CHECK_AND_FILL_VEC_HF_TRACK(eta);
        CHECK_AND_FILL_VEC_HF_TRACK(dcaXY);
        CHECK_AND_FILL_VEC_HF_TRACK(dcaZ);
        CHECK_AND_FILL_VEC_HF_TRACK(sigmaDcaXY);
        CHECK_AND_FILL_VEC_HF_TRACK(sigmaDcaZ);
        CHECK_AND_FILL_VEC_HF_TRACK(normDcaXY);
        CHECK_AND_FILL_VEC_HF_TRACK(normDcaZ);
        CHECK_AND_FILL_VEC_HF_TRACK(signed1Pt);
        CHECK_AND_FILL_VEC_HF_TRACK(tgl);
        CHECK_AND_FILL_VEC_HF_TRACK(sign);
        CHECK_AND_FILL_VEC_HF_TRACK(isPvContributor);
        CHECK_AND_FILL_VEC_HF_TRACK(itsNCls);
        CHECK_AND_FILL_VEC_HF_TRACK(itsNClsInnerBarrel);
        CHECK_AND_FILL_VEC_HF_TRACK(itsChi2NCl);
        CHECK_AND_FILL_VEC_HF_TRACK(tpcNClsFound);
        CHECK_AND_FILL_VEC_HF_TRACK(tpcCrossedRowsOverFindableCls);
        CHECK_AND_FILL_VEC_HF_TRACK(tpcChi2NCl);
        CHECK_AND_FILL_VEC_HF_TRACK(tpcFractionSharedCls);
        CHECK_AND_FILL_VEC_HF_TRACK(tpcNSigmaPi);
        CHECK_AND_FILL_VEC_HF_TRACK(tpcNSigmaKa);
        CHECK_AND_FILL_VEC_HF_TRACK(centrality);
      }
    }

    return inputFeatures;
  }

 protected:
  /// Method to fill the map of available input features
  void setAvailableInputFeatures()
  {
    MlResponse<TypeOutputScore>::mAvailableInputFeatures = {
      FILL_MAP_HF_TRACK(pt),
      FILL_MAP_HF_TRACK(eta),
      FILL_MAP_HF_TRACK(dcaXY),
      FILL_MAP_HF_TRACK(dcaZ),
      FILL_MAP_HF_TRACK(sigmaDcaXY),
      FILL_MAP_HF_TRACK(sigmaDcaZ),
      FILL_MAP_HF_TRACK(normDcaXY),
      FILL_MAP_HF_TRACK(normDcaZ),
      FILL_MAP_HF_TRACK(signed1Pt),
      FILL_MAP_HF_TRACK(tgl),
      FILL_MAP_HF_TRACK(sign),
      FILL_MAP_HF_TRACK(isPvContributor),
      FILL_MAP_HF_TRACK(itsNCls),
      FILL_MAP_HF_TRACK(itsNClsInnerBarrel),
      FILL_MAP_HF_TRACK(itsChi2NCl),
      FILL_MAP_HF_TRACK(tpcNClsFound),
      FILL_MAP_HF_TRACK(tpcCrossedRowsOverFindableCls),
      FILL_MAP_HF_TRACK(tpcChi2NCl),
      FILL_MAP_HF_TRACK(tpcFractionSharedCls),
      FILL_MAP_HF_TRACK(tpcNSigmaPi),
      FILL_MAP_HF_TRACK(tpcNSigmaKa),
      FILL_MAP_HF_TRACK(centrality)};
  }
};

} // namespace o2::analysis

#undef FILL_MAP_HF_TRACK
#undef CHECK_AND_FILL_VEC_HF_TRACK

#endif // PWGHF_CORE_HFMLRESPONSEHFTRACKS_H_
