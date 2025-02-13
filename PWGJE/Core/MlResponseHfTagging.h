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

#include <map>
#include <string>
#include <vector>

#include "Tools/ML/MlResponse.h"
#include "PWGJE/Core/JetTaggingUtilities.h"

#if __has_include(<onnxruntime/core/session/onnxruntime_cxx_api.h>)
#include <onnxruntime/core/session/experimental_onnxruntime_cxx_api.h>
#else
#include <onnxruntime_cxx_api.h>
#endif

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
  JetpT = 0,
  JetEta,
  JetPhi,
  NTracks,
  NSV,
  JetMass,
  TrackpT,
  TrackEta,
  DotProdTrackJet,
  DotProdTrackJetOverJet,
  DeltaRJetTrack,
  SignedIP2D,
  SignedIP2DSign,
  SignedIP3D,
  SignedIP3DSign,
  MomFraction,
  DeltaRTrackVertex,
  TrackPhi,
  TrackCharge,
  TrackITSChi2NCl,
  TrackTPCChi2NCl,
  TrackITSNCls,
  TrackTPCNCls,
  TrackTPCNCrossedRows,
  TrackOrigin,
  TrackVtxIndex,
  SVpT,
  DeltaRSVJet,
  SVMass,
  SVfE,
  IPxy,
  CPA,
  Chi2PCA,
  Dispersion,
  DecayLength2D,
  DecayLength2DError,
  DecayLength3D,
  DecayLength3DError,
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
        CHECK_AND_FILL_VEC_BTAG(jetInput, jet, JetpT)
        CHECK_AND_FILL_VEC_BTAG(jetInput, jet, JetEta)
        CHECK_AND_FILL_VEC_BTAG(jetInput, jet, JetPhi)
        CHECK_AND_FILL_VEC_BTAG(jetInput, jet, NTracks)
        CHECK_AND_FILL_VEC_BTAG(jetInput, jet, NSV)
        CHECK_AND_FILL_VEC_BTAG(jetInput, jet, JetMass)

        default:
          break;
      }
    }

    // Track features
    for (const auto& track : tracks) {
      for (const auto& idx : MlResponse<TypeOutputScore>::mCachedIndices) {
        switch (idx) {
          CHECK_AND_FILL_VEC_BTAG(trackInput, track, TrackpT)
          CHECK_AND_FILL_VEC_BTAG(trackInput, track, TrackEta)
          CHECK_AND_FILL_VEC_BTAG(trackInput, track, DotProdTrackJet)
          CHECK_AND_FILL_VEC_BTAG(trackInput, track, DotProdTrackJetOverJet)
          CHECK_AND_FILL_VEC_BTAG(trackInput, track, DeltaRJetTrack)
          CHECK_AND_FILL_VEC_BTAG(trackInput, track, SignedIP2D)
          CHECK_AND_FILL_VEC_BTAG(trackInput, track, SignedIP2DSign)
          CHECK_AND_FILL_VEC_BTAG(trackInput, track, SignedIP3D)
          CHECK_AND_FILL_VEC_BTAG(trackInput, track, SignedIP3DSign)
          CHECK_AND_FILL_VEC_BTAG(trackInput, track, MomFraction)
          CHECK_AND_FILL_VEC_BTAG(trackInput, track, DeltaRTrackVertex)
          CHECK_AND_FILL_VEC_BTAG(trackInput, track, TrackPhi)
          CHECK_AND_FILL_VEC_BTAG(trackInput, track, TrackCharge)
          CHECK_AND_FILL_VEC_BTAG(trackInput, track, TrackITSChi2NCl)
          CHECK_AND_FILL_VEC_BTAG(trackInput, track, TrackTPCChi2NCl)
          CHECK_AND_FILL_VEC_BTAG(trackInput, track, TrackITSNCls)
          CHECK_AND_FILL_VEC_BTAG(trackInput, track, TrackTPCNCls)
          CHECK_AND_FILL_VEC_BTAG(trackInput, track, TrackTPCNCrossedRows)
          CHECK_AND_FILL_VEC_BTAG(trackInput, track, TrackOrigin)
          CHECK_AND_FILL_VEC_BTAG(trackInput, track, TrackVtxIndex)

          default:
            break;
        }
      }
    }

    // Secondary vertex features
    for (const auto& sv : svs) {
      for (const auto& idx : MlResponse<TypeOutputScore>::mCachedIndices) {
        switch (idx) {
          CHECK_AND_FILL_VEC_BTAG(svInput, sv, SVpT)
          CHECK_AND_FILL_VEC_BTAG(svInput, sv, DeltaRSVJet)
          CHECK_AND_FILL_VEC_BTAG(svInput, sv, SVMass)
          CHECK_AND_FILL_VEC_BTAG(svInput, sv, SVfE)
          CHECK_AND_FILL_VEC_BTAG(svInput, sv, IPxy)
          CHECK_AND_FILL_VEC_BTAG(svInput, sv, CPA)
          CHECK_AND_FILL_VEC_BTAG(svInput, sv, Chi2PCA)
          CHECK_AND_FILL_VEC_BTAG(svInput, sv, Dispersion)
          CHECK_AND_FILL_VEC_BTAG(svInput, sv, DecayLength2D)
          CHECK_AND_FILL_VEC_BTAG(svInput, sv, DecayLength2DError)
          CHECK_AND_FILL_VEC_BTAG(svInput, sv, DecayLength3D)
          CHECK_AND_FILL_VEC_BTAG(svInput, sv, DecayLength3DError)

          default:
            break;
        }
      }
    }
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

    return inputFeatures;
  }

 protected:
  /// Method to fill the map of available input features
  void setAvailableInputFeatures()
  {
    MlResponse<TypeOutputScore>::mAvailableInputFeatures = {
      // Jet features
      FILL_MAP_BJET(JetpT),
      FILL_MAP_BJET(JetEta),
      FILL_MAP_BJET(JetPhi),
      FILL_MAP_BJET(NTracks),
      FILL_MAP_BJET(NSV),
      FILL_MAP_BJET(JetMass),

      // Track features
      FILL_MAP_BJET(TrackpT),
      FILL_MAP_BJET(TrackEta),
      FILL_MAP_BJET(DotProdTrackJet),
      FILL_MAP_BJET(DotProdTrackJetOverJet),
      FILL_MAP_BJET(DeltaRJetTrack),
      FILL_MAP_BJET(SignedIP2D),
      FILL_MAP_BJET(SignedIP2DSign),
      FILL_MAP_BJET(SignedIP3D),
      FILL_MAP_BJET(SignedIP3DSign),
      FILL_MAP_BJET(MomFraction),
      FILL_MAP_BJET(DeltaRTrackVertex),
      FILL_MAP_BJET(TrackPhi),
      FILL_MAP_BJET(TrackCharge),
      FILL_MAP_BJET(TrackITSChi2NCl),
      FILL_MAP_BJET(TrackTPCChi2NCl),
      FILL_MAP_BJET(TrackITSNCls),
      FILL_MAP_BJET(TrackTPCNCls),
      FILL_MAP_BJET(TrackTPCNCrossedRows),
      FILL_MAP_BJET(TrackOrigin),
      FILL_MAP_BJET(TrackVtxIndex),

      // Secondary vertex features
      FILL_MAP_BJET(SVpT),
      FILL_MAP_BJET(DeltaRSVJet),
      FILL_MAP_BJET(SVMass),
      FILL_MAP_BJET(SVfE),
      FILL_MAP_BJET(IPxy),
      FILL_MAP_BJET(CPA),
      FILL_MAP_BJET(Chi2PCA),
      FILL_MAP_BJET(Dispersion),
      FILL_MAP_BJET(DecayLength2D),
      FILL_MAP_BJET(DecayLength2DError),
      FILL_MAP_BJET(DecayLength3D),
      FILL_MAP_BJET(DecayLength3DError)};
  }
};

// ONNX Runtime tensor (Ort::Value) allocator for using customized inputs of ML models.
class TensorAllocator
{
 protected:
#if !__has_include(<onnxruntime/core/session/onnxruntime_cxx_api.h>)
  Ort::MemoryInfo memInfo;
#endif
 public:
  TensorAllocator()
#if !__has_include(<onnxruntime/core/session/onnxruntime_cxx_api.h>)
    : memInfo(Ort::MemoryInfo::CreateCpu(OrtAllocatorType::OrtArenaAllocator, OrtMemType::OrtMemTypeDefault))
#endif
  {
  }
  ~TensorAllocator() = default;
  template <typename T>
  Ort::Value createTensor(std::vector<T>& input, std::vector<int64_t>& inputShape)
  {
#if __has_include(<onnxruntime/core/session/onnxruntime_cxx_api.h>)
    return Ort::Experimental::Value::CreateTensor<T>(input.data(), input.size(), inputShape);
#else
    return Ort::Value::CreateTensor<T>(memInfo, input.data(), input.size(), inputShape.data(), inputShape.size());
#endif
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
    return (feat - tfJetMean[idx]) / tfJetStdev[idx];
  }

  // Track feature normalization
  template <typename T>
  T trkFeatureTransform(T feat, int idx) const
  {
    return (feat - tfTrkMean[idx]) / tfTrkStdev[idx];
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

  // Copy operator for initializing GNNBjetAllocator using Configurable values
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
