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

/// \file MlResponseHfTagging.h
/// \brief Class to compute the ML response for b-jet analysis
/// \author Hadi Hassan <hadi.hassan@cern.ch>, University of Jyväskylä

#ifndef PWGJE_CORE_MLRESPONSEHFTAGGING_H_
#define PWGJE_CORE_MLRESPONSEHFTAGGING_H_

#include "Tools/ML/MlResponse.h"

#include <Framework/Logger.h>

#include <onnxruntime_c_api.h>
#include <onnxruntime_cxx_api.h>

#include <cstddef>
#include <cstdint>
#include <cmath>
#include <utility>
#include <vector>

// Fill the map of available input features
// the key is the feature's name (std::string)
// the value is the corresponding value in EnumInputFeatures
#define FILL_MAP_BJET(FEATURE)                                 \
  {                                                            \
    #FEATURE, static_cast<uint8_t>(InputFeaturesBTag::FEATURE) \
  }

// Check if the index of mCachedIndices (index associated to a FEATURE)
// matches the entry in EnumInputFeatures associated to this FEATURE
// if so, the VECTOR vector is filled with the FEATURE's value
// by calling the corresponding GETTER from OBJECT
#define CHECK_AND_FILL_VEC_BTAG_FULL(VECTOR, OBJECT, FEATURE, GETTER) \
  case static_cast<uint8_t>(InputFeaturesBTag::FEATURE): {            \
    VECTOR.emplace_back(OBJECT.GETTER);                               \
    break;                                                            \
  }

// Specific case of CHECK_AND_FILL_VEC_BTAG_FULL(VECTOR, OBJECT, FEATURE, GETTER)
// where FEATURE = GETTER
#define CHECK_AND_FILL_VEC_BTAG(VECTOR, OBJECT, GETTER)   \
  case static_cast<uint8_t>(InputFeaturesBTag::GETTER): { \
    VECTOR.emplace_back(OBJECT.GETTER);                   \
    break;                                                \
  }

namespace o2::analysis
{
enum class InputFeaturesBTag : uint8_t {
  jetpT = 0,
  jetEta,
  jetPhi,
  nTracks,
  nSV,
  jetMass,
  trackpT,
  trackEta,
  dotProdTrackJet,
  dotProdTrackJetOverJet,
  deltaRJetTrack,
  signedIP2D,
  signedIP2DSign,
  signedIPz,
  signedIPzSign,
  signedIP3D,
  signedIP3DSign,
  momFraction,
  deltaRTrackVertex,
  trackPhi,
  trackCharge,
  trackITSChi2NCl,
  trackTPCChi2NCl,
  trackITSNCls,
  trackTPCNCls,
  trackTPCNCrossedRows,
  trackOrigin,
  trackVtxIndex,
  svpT,
  deltaRSVJet,
  svMass,
  svfE,
  svIPxy,
  svCPA,
  svChi2PCA,
  dispersion,
  decayLength2D,
  decayLength2DError,
  decayLength3D,
  decayLength3DError,
};

template <typename TypeOutputScore = float>
class MlResponseHfTagging : public MlResponse<TypeOutputScore>
{
 public:
  /// Default constructor
  MlResponseHfTagging() = default;
  /// Default destructor
  virtual ~MlResponseHfTagging() = default;

  /// @brief Method to get the input shape of the model
  /// @return A vector of input shapes
  std::vector<std::vector<int64_t>> getInputShape() const { return this->mModels[0].getInputShapes(); }

  /// @brief Method to get the output shape of a model
  /// \param imod is the index of the model
  /// @return number of output nodes
  int getOutputNodes(int imod = 0) const { return this->mModels[imod].getNumOutputNodes(); }

  /// Method to fill the inputs of jet, tracks and secondary vertices
  /// \param jet is the b-jet candidate
  /// \param tracks is the vector of tracks associated to the jet
  /// \param svs is the vector of secondary vertices associated to the jet
  /// \param jetInput the jet input features vector to be filled
  /// \param trackInput the tracks input features vector to be filled
  /// \param svInput the SVs input features vector to be filled
  template <typename T1, typename T2, typename T3>
  void fillInputFeatures(T1 const& jet, std::vector<T2> const& tracks, std::vector<T3> const& svs, std::vector<float>& jetInput, std::vector<float>& trackInput, std::vector<float>& svInput)
  {

    // Jet features
    for (const auto& idx : MlResponse<TypeOutputScore>::mCachedIndices) {
      switch (idx) {
        CHECK_AND_FILL_VEC_BTAG(jetInput, jet, jetpT)
        CHECK_AND_FILL_VEC_BTAG(jetInput, jet, jetEta)
        CHECK_AND_FILL_VEC_BTAG(jetInput, jet, jetPhi)
        CHECK_AND_FILL_VEC_BTAG(jetInput, jet, nTracks)
        CHECK_AND_FILL_VEC_BTAG(jetInput, jet, nSV)
        CHECK_AND_FILL_VEC_BTAG(jetInput, jet, jetMass)

        default:
          break;
      }
    }

    // Track features
    for (const auto& track : tracks) {
      for (const auto& idx : MlResponse<TypeOutputScore>::mCachedIndices) {
        switch (idx) {
          CHECK_AND_FILL_VEC_BTAG(trackInput, track, trackpT)
          CHECK_AND_FILL_VEC_BTAG(trackInput, track, trackEta)
          CHECK_AND_FILL_VEC_BTAG(trackInput, track, dotProdTrackJet)
          CHECK_AND_FILL_VEC_BTAG(trackInput, track, dotProdTrackJetOverJet)
          CHECK_AND_FILL_VEC_BTAG(trackInput, track, deltaRJetTrack)
          CHECK_AND_FILL_VEC_BTAG(trackInput, track, signedIP2D)
          CHECK_AND_FILL_VEC_BTAG(trackInput, track, signedIP2DSign)
          CHECK_AND_FILL_VEC_BTAG(trackInput, track, signedIPz)
          CHECK_AND_FILL_VEC_BTAG(trackInput, track, signedIPzSign)
          CHECK_AND_FILL_VEC_BTAG(trackInput, track, signedIP3D)
          CHECK_AND_FILL_VEC_BTAG(trackInput, track, signedIP3DSign)
          CHECK_AND_FILL_VEC_BTAG(trackInput, track, momFraction)
          CHECK_AND_FILL_VEC_BTAG(trackInput, track, deltaRTrackVertex)
          CHECK_AND_FILL_VEC_BTAG(trackInput, track, trackPhi)
          CHECK_AND_FILL_VEC_BTAG(trackInput, track, trackCharge)
          CHECK_AND_FILL_VEC_BTAG(trackInput, track, trackITSChi2NCl)
          CHECK_AND_FILL_VEC_BTAG(trackInput, track, trackTPCChi2NCl)
          CHECK_AND_FILL_VEC_BTAG(trackInput, track, trackITSNCls)
          CHECK_AND_FILL_VEC_BTAG(trackInput, track, trackTPCNCls)
          CHECK_AND_FILL_VEC_BTAG(trackInput, track, trackTPCNCrossedRows)
          CHECK_AND_FILL_VEC_BTAG(trackInput, track, trackOrigin)
          CHECK_AND_FILL_VEC_BTAG(trackInput, track, trackVtxIndex)

          default:
            break;
        }
      }
    }

    // Secondary vertex features
    for (const auto& sv : svs) {
      for (const auto& idx : MlResponse<TypeOutputScore>::mCachedIndices) {
        switch (idx) {
          CHECK_AND_FILL_VEC_BTAG(svInput, sv, svpT)
          CHECK_AND_FILL_VEC_BTAG(svInput, sv, deltaRSVJet)
          CHECK_AND_FILL_VEC_BTAG(svInput, sv, svMass)
          CHECK_AND_FILL_VEC_BTAG(svInput, sv, svfE)
          CHECK_AND_FILL_VEC_BTAG(svInput, sv, svIPxy)
          CHECK_AND_FILL_VEC_BTAG(svInput, sv, svCPA)
          CHECK_AND_FILL_VEC_BTAG(svInput, sv, svChi2PCA)
          CHECK_AND_FILL_VEC_BTAG(svInput, sv, dispersion)
          CHECK_AND_FILL_VEC_BTAG(svInput, sv, decayLength2D)
          CHECK_AND_FILL_VEC_BTAG(svInput, sv, decayLength2DError)
          CHECK_AND_FILL_VEC_BTAG(svInput, sv, decayLength3D)
          CHECK_AND_FILL_VEC_BTAG(svInput, sv, decayLength3DError)

          default:
            break;
        }
      }
    }
  }

  /// @brief Method to replace NaN and infinity values in a vector with a specified value
  /// @param vec is the vector to be processed
  /// @param value is the value to replace NaN values with
  /// @return the number of NaN values replaced
  template <typename T>
  static int replaceNaN(std::vector<T>& vec, T value)
  {
    int numNaN = 0;
    for (auto& el : vec) {
      if (std::isnan(el) || std::isinf(el)) {
        el = value;
        ++numNaN;
      }
    }
    return numNaN;
  }

  /// Method to get the input features vector needed for ML inference in a 2D vector
  /// \param jet is the b-jet candidate
  /// \param tracks is the vector of tracks associated to the jet
  /// \param svs is the vector of secondary vertices associated to the jet
  /// \return inputFeatures vector
  template <typename T1, typename T2, typename T3>
  std::vector<std::vector<float>> getInputFeatures2D(T1 const& jet, std::vector<T2> const& tracks, std::vector<T3> const& svs)
  {

    std::vector<float> jetInput;
    std::vector<float> trackInput;
    std::vector<float> svInput;

    fillInputFeatures(jet, tracks, svs, jetInput, trackInput, svInput);

    std::vector<std::vector<float>> inputFeatures;

    replaceNaN(jetInput, 0.f);
    replaceNaN(trackInput, 0.f);
    replaceNaN(svInput, 0.f);

    inputFeatures.push_back(jetInput);
    inputFeatures.push_back(trackInput);
    inputFeatures.push_back(svInput);

    return inputFeatures;
  }

  /// Method to get the input features vector needed for ML inference in a 1D vector
  /// \param jet is the b-jet candidate
  /// \param tracks is the vector of tracks associated to the jet
  /// \param svs is the vector of secondary vertices associated to the jet
  /// \return inputFeatures vector
  template <typename T1, typename T2, typename T3>
  std::vector<float> getInputFeatures1D(T1 const& jet, std::vector<T2> const& tracks, std::vector<T3> const& svs)
  {

    std::vector<float> jetInput;
    std::vector<float> trackInput;
    std::vector<float> svInput;

    fillInputFeatures(jet, tracks, svs, jetInput, trackInput, svInput);

    std::vector<float> inputFeatures;

    inputFeatures.insert(inputFeatures.end(), jetInput.begin(), jetInput.end());
    inputFeatures.insert(inputFeatures.end(), trackInput.begin(), trackInput.end());
    inputFeatures.insert(inputFeatures.end(), svInput.begin(), svInput.end());

    replaceNaN(inputFeatures, 0.f);

    return inputFeatures;
  }

 protected:
  /// Method to fill the map of available input features
  void setAvailableInputFeatures()
  {
    MlResponse<TypeOutputScore>::mAvailableInputFeatures = {
      // Jet features
      FILL_MAP_BJET(jetpT),
      FILL_MAP_BJET(jetEta),
      FILL_MAP_BJET(jetPhi),
      FILL_MAP_BJET(nTracks),
      FILL_MAP_BJET(nSV),
      FILL_MAP_BJET(jetMass),

      // Track features
      FILL_MAP_BJET(trackpT),
      FILL_MAP_BJET(trackEta),
      FILL_MAP_BJET(dotProdTrackJet),
      FILL_MAP_BJET(dotProdTrackJetOverJet),
      FILL_MAP_BJET(deltaRJetTrack),
      FILL_MAP_BJET(signedIP2D),
      FILL_MAP_BJET(signedIP2DSign),
      FILL_MAP_BJET(signedIPz),
      FILL_MAP_BJET(signedIPzSign),
      FILL_MAP_BJET(signedIP3D),
      FILL_MAP_BJET(signedIP3DSign),
      FILL_MAP_BJET(momFraction),
      FILL_MAP_BJET(deltaRTrackVertex),
      FILL_MAP_BJET(trackPhi),
      FILL_MAP_BJET(trackCharge),
      FILL_MAP_BJET(trackITSChi2NCl),
      FILL_MAP_BJET(trackTPCChi2NCl),
      FILL_MAP_BJET(trackITSNCls),
      FILL_MAP_BJET(trackTPCNCls),
      FILL_MAP_BJET(trackTPCNCrossedRows),
      FILL_MAP_BJET(trackOrigin),
      FILL_MAP_BJET(trackVtxIndex),

      // Secondary vertex features
      FILL_MAP_BJET(svpT),
      FILL_MAP_BJET(deltaRSVJet),
      FILL_MAP_BJET(svMass),
      FILL_MAP_BJET(svfE),
      FILL_MAP_BJET(svIPxy),
      FILL_MAP_BJET(svCPA),
      FILL_MAP_BJET(svChi2PCA),
      FILL_MAP_BJET(dispersion),
      FILL_MAP_BJET(decayLength2D),
      FILL_MAP_BJET(decayLength2DError),
      FILL_MAP_BJET(decayLength3D),
      FILL_MAP_BJET(decayLength3DError)};
  }
};

// ONNX Runtime tensor (Ort::Value) allocator for using customized inputs of ML models.
class TensorAllocator
{
 protected:
  Ort::MemoryInfo memInfo;

 public:
  TensorAllocator()
    : memInfo(Ort::MemoryInfo::CreateCpu(OrtAllocatorType::OrtArenaAllocator, OrtMemType::OrtMemTypeDefault))
  {
  }
  ~TensorAllocator() = default;
  template <typename T>
  Ort::Value createTensor(std::vector<T>& input, std::vector<int64_t>& inputShape)
  {
    return Ort::Value::CreateTensor<T>(memInfo, input.data(), input.size(), inputShape.data(), inputShape.size());
  }
};

// TensorAllocator for GNN b-jet tagger
class GNNBjetAllocator : public TensorAllocator
{
 private:
  int64_t nJetFeat;
  int64_t nTrkFeat;
  int64_t nFlav;
  int64_t nTrkOrigin;
  int64_t maxNNodes;

  std::vector<float> tfJetMean;
  std::vector<float> tfJetStdev;
  std::vector<float> tfTrkMean;
  std::vector<float> tfTrkStdev;

  std::vector<std::vector<int64_t>> edgesList;

  // Jet feature normalization
  template <typename T>
  T jetFeatureTransform(T feat, int idx) const
  {
    return std::tanh((feat - tfJetMean[idx]) / tfJetStdev[idx]);
  }

  // Track feature normalization
  template <typename T>
  T trkFeatureTransform(T feat, int idx) const
  {
    return std::tanh((feat - tfTrkMean[idx]) / tfTrkStdev[idx]);
  }

  // Edge input of GNN (fully-connected graph)
  void setEdgesList(void)
  {
    for (int64_t nNodes = 0; nNodes <= maxNNodes; ++nNodes) {
      std::vector<std::pair<int64_t, int64_t>> edges;
      // Generate all permutations of (i, j) where i != j
      for (int64_t i = 0; i < nNodes; ++i) {
        for (int64_t j = 0; j < nNodes; ++j) {
          if (i != j) {
            edges.emplace_back(i, j);
          }
        }
      }
      // Add self-loops (i, i)
      for (int64_t i = 0; i < nNodes; ++i) {
        edges.emplace_back(i, i);
      }
      // Flatten
      std::vector<int64_t> flattenedEdges;
      for (const auto& edge : edges) {
        flattenedEdges.push_back(edge.first);
      }
      for (const auto& edge : edges) {
        flattenedEdges.push_back(edge.second);
      }
      edgesList.push_back(flattenedEdges);
    }
  }

  // Replace NaN in a vector into value
  template <typename T>
  static int replaceNaN(std::vector<T>& vec, T value)
  {
    int numNaN = 0;
    for (auto& el : vec) { // o2-linter: disable=const-ref-in-for-loop
      if (std::isnan(el)) {
        el = value;
        ++numNaN;
      }
    }
    return numNaN;
  }

 public:
  GNNBjetAllocator() : TensorAllocator(), nJetFeat(4), nTrkFeat(13), nFlav(3), nTrkOrigin(5), maxNNodes(40) {}
  GNNBjetAllocator(int64_t nJetFeat, int64_t nTrkFeat, int64_t nFlav, int64_t nTrkOrigin, std::vector<float>& tfJetMean, std::vector<float>& tfJetStdev, std::vector<float>& tfTrkMean, std::vector<float>& tfTrkStdev, int64_t maxNNodes = 40)
    : TensorAllocator(), nJetFeat(nJetFeat), nTrkFeat(nTrkFeat), nFlav(nFlav), nTrkOrigin(nTrkOrigin), maxNNodes(maxNNodes), tfJetMean(tfJetMean), tfJetStdev(tfJetStdev), tfTrkMean(tfTrkMean), tfTrkStdev(tfTrkStdev)
  {
    setEdgesList();
  }
  ~GNNBjetAllocator() = default;

  // Copy operator for initializing GNNBjetAllocator using o2::framework::Configurable values
  GNNBjetAllocator& operator=(const GNNBjetAllocator& other)
  {
    nJetFeat = other.nJetFeat;
    nTrkFeat = other.nTrkFeat;
    nFlav = other.nFlav;
    nTrkOrigin = other.nTrkOrigin;
    maxNNodes = other.maxNNodes;
    tfJetMean = other.tfJetMean;
    tfJetStdev = other.tfJetStdev;
    tfTrkMean = other.tfTrkMean;
    tfTrkStdev = other.tfTrkStdev;
    setEdgesList();
    return *this;
  }

  // Allocate & Return GNN input tensors (std::vector<Ort::Value>)
  template <typename T>
  void getGNNInput(std::vector<T>& jetFeat, std::vector<std::vector<T>>& trkFeat, std::vector<T>& feat, std::vector<Ort::Value>& gnnInput)
  {
    int64_t nNodes = trkFeat.size();

    std::vector<int64_t> edgesShape{2, nNodes * nNodes};
    gnnInput.emplace_back(createTensor(edgesList[nNodes], edgesShape));

    std::vector<int64_t> featShape{nNodes, nJetFeat + nTrkFeat};

    int numNaN = replaceNaN(jetFeat, 0.f);
    for (auto& aTrkFeat : trkFeat) { // o2-linter: disable=const-ref-in-for-loop
      for (size_t i = 0; i < jetFeat.size(); ++i)
        feat.push_back(jetFeatureTransform(jetFeat[i], i));
      numNaN += replaceNaN(aTrkFeat, 0.f);
      for (size_t i = 0; i < aTrkFeat.size(); ++i)
        feat.push_back(trkFeatureTransform(aTrkFeat[i], i));
    }

    gnnInput.emplace_back(createTensor(feat, featShape));

    if (numNaN > 0) {
      LOGF(info, "NaN found in GNN input feature, number of NaN: %d", numNaN);
    }
  }
};

} // namespace o2::analysis

#undef FILL_MAP_BJET
#undef CHECK_AND_FILL_VEC_BTAG_FULL
#undef CHECK_AND_FILL_VEC_BTAG

#endif // PWGJE_CORE_MLRESPONSEHFTAGGING_H_
