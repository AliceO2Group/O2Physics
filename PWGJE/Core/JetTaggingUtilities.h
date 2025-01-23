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

/// \file JetTaggingUtilities.h
/// \brief Jet tagging related utilities
///
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>
/// \author Hanseo Park <hanseo.park@cern.ch>

#ifndef PWGJE_CORE_JETTAGGINGUTILITIES_H_
#define PWGJE_CORE_JETTAGGINGUTILITIES_H_

#include <cmath>
#include <limits>
#include <numeric>
#include <tuple>
#include <vector>
#include <algorithm>
#include <functional>
#include <memory>
#include <string>
#include <unordered_map>

#include <TPDGCode.h>
#include "CommonConstants/PhysicsConstants.h"

#include "TF1.h"
#include "Framework/Logger.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "PWGJE/Core/JetUtilities.h"

#if __has_include(<onnxruntime/core/session/onnxruntime_cxx_api.h>)
#include <onnxruntime/core/session/experimental_onnxruntime_cxx_api.h>
#else
#include <onnxruntime_cxx_api.h>
#endif

using namespace o2::constants::physics;

enum JetTaggingSpecies {
  none = 0,
  charm = 1,
  beauty = 2,
  lightflavour = 3,
  lightquark = 4,
  gluon = 5
};

enum BJetTaggingMethod {
  IPsN1 = 0,
  IPsN2 = 1,
  IPsN3 = 2,
  IPs3DN1 = 3,
  IPs3DN2 = 4,
  IPs3DN3 = 5,
  SV = 6,
  SV3D = 7
};

namespace jettaggingutilities
{
const int cmTomum = 10000; // using cm -> #mum for impact parameter (dca)

struct BJetParams {
  float mJetpT = 0.0;
  float mJetEta = 0.0;
  float mJetPhi = 0.0;
  int mNTracks = -1;
  int mNSV = -1;
  float mJetMass = 0.0;
};

struct BJetTrackParams {
  double mTrackpT = 0.0;
  double mTrackEta = 0.0;
  double mDotProdTrackJet = 0.0;
  double mDotProdTrackJetOverJet = 0.0;
  double mDeltaRJetTrack = 0.0;
  double mSignedIP2D = 0.0;
  double mSignedIP2DSign = 0.0;
  double mSignedIP3D = 0.0;
  double mSignedIP3DSign = 0.0;
  double mMomFraction = 0.0;
  double mDeltaRTrackVertex = 0.0;
};

struct BJetSVParams {
  double mSVpT = 0.0;
  double mDeltaRSVJet = 0.0;
  double mSVMass = 0.0;
  double mSVfE = 0.0;
  double mIPXY = 0.0;
  double mCPA = 0.0;
  double mChi2PCA = 0.0;
  double mDispersion = 0.0;
  double mDecayLength2D = 0.0;
  double mDecayLength2DError = 0.0;
  double mDecayLength3D = 0.0;
  double mDecayLength3DError = 0.0;
};

// ONNX Runtime tensor (Ort::Value) allocator for using customized inputs of ML models.
class TensorAllocator
{
 protected:
#if !__has_include(<onnxruntime/core/session/onnxruntime_cxx_api.h>)
  Ort::MemoryInfo mem_info;
#endif
 public:
  TensorAllocator()
#if !__has_include(<onnxruntime/core/session/onnxruntime_cxx_api.h>)
    : mem_info(Ort::MemoryInfo::CreateCpu(OrtAllocatorType::OrtArenaAllocator, OrtMemType::OrtMemTypeDefault))
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
    return Ort::Value::CreateTensor<T>(mem_info, input.data(), input.size(), inputShape.data(), inputShape.size());
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
    for (auto& el : vec) {
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
    for (auto& aTrkFeat : trkFeat) {
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

//________________________________________________________________________
bool isBHadron(int pc)
{
  std::vector<int> bPdG = {Pdg::kB0, Pdg::kBPlus, 10511, 10521, 513, 523, 10513, 10523, 20513, 20523, 20513, 20523, 515, 525, Pdg::kBS, 10531, 533, 10533,
                           20533, 535, 541, 10541, 543, 10543, 20543, 545, 551, 10551, 100551, 110551, 200551, 210551, 553, 10553, 20553,
                           30553, 100553, 110553, 120553, 130553, 200553, 210553, 220553, 300553, 9000533, 9010553, 555, 10555, 20555,
                           100555, 110555, 120555, 200555, 557, 100557, Pdg::kLambdaB0, 5112, 5212, 5222, 5114, 5214, 5224, 5132, kXiB0, 5312, 5322,
                           5314, 5324, 5332, 5334, 5142, 5242, 5412, 5422, 5414, 5424, 5342, 5432, 5434, 5442, 5444, 5512, 5522, 5514, 5524,
                           5532, 5534, 5542, 5544, 5554};

  return (std::find(bPdG.begin(), bPdG.end(), std::abs(pc)) != bPdG.end());
}
//________________________________________________________________________
bool isCHadron(int pc)
{
  std::vector<int> bPdG = {Pdg::kDPlus, Pdg::kD0, Pdg::kD0StarPlus, Pdg::kD0Star0, 413, 423, 10413, 10423, 20431, 20423, Pdg::kD2StarPlus, Pdg::kD2Star0, Pdg::kDS, 10431, Pdg::kDSStar, Pdg::kDS1, 20433, Pdg::kDS2Star, 441,
                           10441, 100441, Pdg::kJPsi, 10443, Pdg::kChiC1, 100443, 30443, 9000443, 9010443, 9020443, 445, 100445, Pdg::kLambdaCPlus, Pdg::kSigmaCPlusPlus, 4212, Pdg::kSigmaC0,
                           4224, 4214, 4114, Pdg::kXiCPlus, Pdg::kXiC0, 4322, 4312, 4324, 4314, Pdg::kOmegaC0, 4334, 4412, Pdg::kXiCCPlusPlus, 4414, 4424, 4432, 4434, 4444};

  return (std::find(bPdG.begin(), bPdG.end(), std::abs(pc)) != bPdG.end());
}

/**
 * returns the globalIndex of the earliest mother of a particle in the shower. returns -1 if a suitable mother is not found
 *
 * @param particle MCParticle whose mother is to be found
 */
template <typename T>
int getOriginalMotherIndex(const typename T::iterator& particle)
{

  // if (!particle) {
  //   return -1.0;
  // }
  auto mother = particle;

  while (mother.has_mothers()) {

    mother = mother.template mothers_first_as<T>();

    int motherStatusCode = std::abs(mother.getGenStatusCode());

    if (motherStatusCode == 23 || motherStatusCode == 33 || motherStatusCode == 43 || motherStatusCode == 63) {
      return mother.globalIndex();
    }
  }
  return -1.0;
}

/**
 * returns the globalIndex of the earliest HF mother of a particle in the shower. returns -1 if a suitable mother is not found. Should be used only on already identified HF particles
 *
 * @param hfparticle MCParticle whose mother is to be found
 */

template <typename T>
int getOriginalHFMotherIndex(const typename T::iterator& hfparticle)
{

  // if (!hfparticle) {
  //   return -1.0;
  // }
  auto mother = hfparticle;

  while (mother.has_mothers()) {

    mother = mother.template mothers_first_as<T>();

    int motherStatusCode = std::abs(mother.getGenStatusCode());

    if (motherStatusCode == 23 || motherStatusCode == 33 || motherStatusCode == 43 || motherStatusCode == 63 || (motherStatusCode == 51 && mother.template mothers_first_as<T>().pdgCode() == 21)) {
      return mother.globalIndex();
    }
  }
  return -1.0;
}

/**
 * checks if atrack in a reco level jet originates from a HF shower. 0:no HF shower, 1:charm shower, 2:beauty shower. The first track originating from an HF shower can be extracted by reference
 *
 * @param jet
 * @param particles table of generator level particles to be searched through
 * @param hftrack track passed as reference which is then replaced by the first track that originated from an HF shower
 */
template <typename T, typename U, typename V>
int jetTrackFromHFShower(T const& jet, U const& /*tracks*/, V const& particles, typename U::iterator& hftrack, bool searchUpToQuark)
{

  bool hasMcParticle = false;
  int origin = -1;
  for (auto const& track : jet.template tracks_as<U>()) {
    hftrack = track; // for init if origin is 1 or 2, the track is not hftrack
    if (!track.has_mcParticle()) {
      continue;
    }
    hasMcParticle = true;
    auto const& particle = track.template mcParticle_as<V>();
    origin = RecoDecay::getParticleOrigin(particles, particle, searchUpToQuark);
    if (origin == 1 || origin == 2) { // 1=charm , 2=beauty
      hftrack = track;
      if (origin == 1) {
        return JetTaggingSpecies::charm;
      }
      if (origin == 2) {
        return JetTaggingSpecies::beauty;
      }
    }
  }

  if (hasMcParticle) {
    return JetTaggingSpecies::lightflavour;
  } else {
    return JetTaggingSpecies::none;
  }
}

/**
 * checks if a particle in a generator level jet originates from a HF shower. 0:no HF shower, 1:charm shower, 2:beauty shower. The first particle originating from an HF shower can be extracted by reference
 *
 * @param jet
 * @param particles table of generator level particles to be searched through
 * @param hfparticle particle passed as reference which is then replaced by the first track that originated from an HF shower
 */
template <typename T, typename U>
int jetParticleFromHFShower(T const& jet, U const& particles, typename U::iterator& hfparticle, bool searchUpToQuark)
{

  int origin = -1;
  for (const auto& particle : jet.template tracks_as<U>()) {
    hfparticle = particle; // for init if origin is 1 or 2, the particle is not hfparticle
    origin = RecoDecay::getParticleOrigin(particles, particle, searchUpToQuark);
    if (origin == 1 || origin == 2) { // 1=charm , 2=beauty
      hfparticle = particle;
      if (origin == 1) {
        return JetTaggingSpecies::charm;
      }
      if (origin == 2) {
        return JetTaggingSpecies::beauty;
      }
    }
  }
  return JetTaggingSpecies::lightflavour;
}

/**
 * returns if a reco level jet originates from a HF shower. 0:no HF shower, 1:charm shower, 2:beauty shower. The requirement is that the jet contains a particle from an HF shower and that the original HF quark is within dRMax of the jet axis in eta-phi
 *
 * @param jet
 * @param particles table of generator level particles to be searched through
 * @param dRMax maximum distance in eta-phi of initiating heavy-flavour quark from the jet axis
 */

template <typename T, typename U, typename V>
int mcdJetFromHFShower(T const& jet, U const& tracks, V const& particles, float dRMax = 0.25, bool searchUpToQuark = false)
{

  typename U::iterator hftrack;
  int origin = jetTrackFromHFShower(jet, tracks, particles, hftrack, searchUpToQuark);
  if (origin == JetTaggingSpecies::charm || origin == JetTaggingSpecies::beauty) {
    if (!hftrack.has_mcParticle()) {
      return JetTaggingSpecies::none;
    }
    auto const& hfparticle = hftrack.template mcParticle_as<V>();

    int originalHFMotherIndex = getOriginalHFMotherIndex<V>(hfparticle);
    if (originalHFMotherIndex > -1.0) {

      if (jetutilities::deltaR(jet, particles.iteratorAt(originalHFMotherIndex)) < dRMax) {

        return origin;

      } else {
        return JetTaggingSpecies::none;
      }

    } else {
      return JetTaggingSpecies::none;
    }

  } else {

    return JetTaggingSpecies::lightflavour;
  }
}

/**
 * checks if a generator level jet originates from a HF shower. 0:no HF shower, 1:charm shower, 2:beauty shower. The requirement is that the jet contains a particle from an HF shower and that the original HF quark is within dRMax of the jet axis in eta-phi
 *
 * @param jet
 * @param particles table of generator level particles to be searched through
 * @param dRMax maximum distance in eta-phi of initiating heavy-flavour quark from the jet axis
 */

template <typename T, typename U>
int mcpJetFromHFShower(T const& jet, U const& particles, float dRMax = 0.25, bool searchUpToQuark = false)
{

  typename U::iterator hfparticle;
  int origin = jetParticleFromHFShower(jet, particles, hfparticle, searchUpToQuark);
  if (origin == JetTaggingSpecies::charm || origin == JetTaggingSpecies::beauty) {

    int originalHFMotherIndex = getOriginalHFMotherIndex<U>(hfparticle);
    if (originalHFMotherIndex > -1.0) {

      if (jetutilities::deltaR(jet, particles.iteratorAt(originalHFMotherIndex)) < dRMax) {

        return origin;

      } else {
        return JetTaggingSpecies::none;
      }

    } else {
      return JetTaggingSpecies::none;
    }

  } else {

    return JetTaggingSpecies::lightflavour;
  }
}

/**
 * returns the pdg code of the original scattered parton closest to the jet axis, with the restriction that the parton and jet axis must be within dRMax in eta-phi
 *
 * @param jet
 * @param particles table of generator level particles to be searched through
 * @param dRMax maximum distance in eta-phi of initiating heavy-flavour quark from the jet axis
 */

template <typename T, typename U>
int jetOrigin(T const& jet, U const& particles, float dRMax = 0.25)
{
  bool firstPartonFound = false;
  typename U::iterator parton1;
  typename U::iterator parton2;
  for (auto const& particle : particles) {
    if (std::abs(particle.getGenStatusCode() == 23)) {
      if (!firstPartonFound) {
        parton1 = particle;
        firstPartonFound = true;
      } else {
        parton2 = particle;
      }
    }
  }

  float dR1 = jetutilities::deltaR(jet, parton1);
  float dR2 = jetutilities::deltaR(jet, parton2);

  if (dR1 <= dR2 && dR1 < dRMax) {

    return parton1.pdgCode();
  }
  if (dR2 <= dR1 && dR2 < dRMax) {

    return parton2.pdgCode();
  }

  return 0;
}

/**
 * return the jet flavor: 0 for lf-jet, 1 for c-jet, 2 for b-jet
 *
 * @param AnyJet the jet that we need to study its flavor
 * @param AllMCParticles a vector of all the mc particles stack
 */
template <typename AnyJet, typename AllMCParticles>
int16_t getJetFlavor(AnyJet const& jet, AllMCParticles const& mcparticles)
{
  bool charmQuark = false;
  for (auto const& mcpart : mcparticles) {
    int pdgcode = mcpart.pdgCode();
    if (std::abs(pdgcode) == 21 || (std::abs(pdgcode) >= 1 && std::abs(pdgcode) <= 5)) {
      double dR = jetutilities::deltaR(jet, mcpart);

      if (dR < jet.r() / 100.f) {
        if (std::abs(pdgcode) == 5) {
          return JetTaggingSpecies::beauty; // Beauty jet
        } else if (std::abs(pdgcode) == 4) {
          charmQuark = true;
        }
      }
    }
  }

  if (charmQuark) {
    return JetTaggingSpecies::charm; // Charm jet
  }

  return JetTaggingSpecies::lightflavour; // Light flavor jet
}

/**
 * return the jet flavor if it finds a HF hadron inside the jet: 0 for lf-jet, 1 for c-jet, 2 for b-jet
 *
 * @param AnyJet the jet that we need to study its flavor
 * @param AllMCParticles a vector of all the mc particles stack
 */
template <typename AnyJet, typename AllMCParticles>
int16_t getJetFlavorHadron(AnyJet const& jet, AllMCParticles const& mcparticles)
{
  bool charmHadron = false;

  for (auto const& mcpart : mcparticles) {
    int pdgcode = mcpart.pdgCode();
    if (isBHadron(pdgcode) || isCHadron(pdgcode)) {
      double dR = jetutilities::deltaR(jet, mcpart);

      if (dR < jet.r() / 100.f) {
        if (isBHadron(pdgcode)) {
          return JetTaggingSpecies::beauty; // Beauty jet
        } else if (isCHadron(pdgcode)) {
          charmHadron = true;
        }
      }
    }
  }

  if (charmHadron) {
    return JetTaggingSpecies::charm; // Charm jet
  }

  return JetTaggingSpecies::lightflavour; // Light flavor jet
}

/**
 * return acceptance of track about DCA xy and z due to cut for QualityTracks
 */
template <typename T>
bool trackAcceptanceWithDca(T const& track, float trackDcaXYMax, float trackDcaZMax)
{
  if (std::abs(track.dcaXY()) > trackDcaXYMax)
    return false;
  if (std::abs(track.dcaZ()) > trackDcaZMax)
    return false;
  return true;
}

/**
 * retrun acceptance of prong due to cut for high quality secondary vertex
 */
template <typename T>
bool prongAcceptance(T const& prong, float prongChi2PCAMin, float prongChi2PCAMax, float prongsigmaLxyMax, float prongIPxyMin, float prongIPxyMax, bool doXYZ)
{
  if (prong.chi2PCA() < prongChi2PCAMin)
    return false;
  if (prong.chi2PCA() > prongChi2PCAMax)
    return false;
  if (std::abs(prong.impactParameterXY()) < prongIPxyMin)
    return false;
  if (std::abs(prong.impactParameterXY()) > prongIPxyMax)
    return false;

  if (!doXYZ) {
    if (prong.errorDecayLengthXY() > prongsigmaLxyMax)
      return false;
  } else {
    if (prong.errorDecayLength() > prongsigmaLxyMax)
      return false;
  }
  return true;
}

/**
 * retrun acceptance of secondary vertex due to cut for high quality secondary vertex
 */
template <typename T>
bool svAcceptance(T const& sv, float svDispersionMax)
{
  if (sv.dispersion() > svDispersionMax)
    return false;
  return true;
}

/**
 * return geometric sign which is calculated scalar product between jet axis with DCA (track propagated to PV )
 * positive and negative value are expected from primary vertex
 * positive value is expected from secondary vertex
 *
 * @param jet
 * @param track which is needed aod::JTrackExtras
 */
template <typename T, typename U>
int getGeoSign(T const& jet, U const& track)
{
  auto sign = TMath::Sign(1, track.dcaX() * jet.px() + track.dcaY() * jet.py() + track.dcaZ() * jet.pz());
  if (sign < -1 || sign > 1)
    LOGF(info, Form("Sign is %d", sign));
  return sign;
}

/**
 * Orders the tracks associated with a jet based on signed impact parameter significance and stores them
 * in a vector in descending order.
 */
template <typename T, typename U, typename Vec = std::vector<float>>
void orderForIPJetTracks(T const& jet, U const& /*tracks*/, float trackDcaXYMax, float trackDcaZMax, Vec& vecSignImpSig, bool useIPxyz)
{
  for (auto const& track : jet.template tracks_as<U>()) {
    if (!trackAcceptanceWithDca(track, trackDcaXYMax, trackDcaZMax))
      continue;
    auto geoSign = getGeoSign(jet, track);
    float varSignImpSig;
    if (!useIPxyz) {
      varSignImpSig = geoSign * std::abs(track.dcaXY()) / track.sigmadcaXY();
    } else {
      varSignImpSig = geoSign * std::abs(track.dcaXYZ()) / track.sigmadcaXYZ();
    }
    vecSignImpSig.push_back(varSignImpSig);
  }
  std::sort(vecSignImpSig.begin(), vecSignImpSig.end(), std::greater<float>());
}

/**
 * Checks if a jet is greater than the given tagging working point based on the signed impact parameter significances
 * return (true, true, true) if the jet is tagged by the 1st, 2nd and 3rd largest IPs
 */
template <typename T, typename U>
std::tuple<bool, bool, bool> isGreaterThanTaggingPoint(T const& jet, U const& tracks, float trackDcaXYMax, float trackDcaZMax, float taggingPoint = 1.0, bool useIPxyz = false)
{
  bool taggedIPsN1 = false;
  bool taggedIPsN2 = false;
  bool taggedIPsN3 = false;
  std::vector<float> vecSignImpSig;
  orderForIPJetTracks(jet, tracks, trackDcaXYMax, trackDcaZMax, vecSignImpSig, useIPxyz);
  if (vecSignImpSig.size() > 0) {
    if (vecSignImpSig[0] > taggingPoint) { // tagger point set
      taggedIPsN1 = true;
    }
  }
  if (vecSignImpSig.size() > 1) {
    if (vecSignImpSig[1] > taggingPoint) { // tagger point set
      taggedIPsN2 = true;
    }
  }
  if (vecSignImpSig.size() > 2) {
    if (vecSignImpSig[2] > taggingPoint) { // tagger point set
      taggedIPsN3 = true;
    }
  }
  return std::make_tuple(taggedIPsN1, taggedIPsN2, taggedIPsN3);
}

/**
 * Creates and sets the parameters of a resolution function (TF1) for constituents of jet based on signed impact parameter significnace in plane XY
 * This function is typically used to set up a resolution function for jet tagging purposes.
 *
 * @param vecParams A vector containing the parameters for the resolution function.
 * @return A unique pointer to the TF1 resolution function with the parameters set.
 */
template <typename T>
std::unique_ptr<TF1> setResolutionFunction(T const& vecParams)
{
  std::unique_ptr<TF1> fResoFunc(new TF1("fResoFunc", "gaus(0)+expo(3)+expo(5)+expo(7)", -40, 0));
  for (typename T::size_type i = 0; i < vecParams.size(); i++) {
    fResoFunc->SetParameter(i, vecParams[i]);
  }

  return fResoFunc;
}

/**
 * Calculates the probability of a given track being associated with a jet, based on the geometric
 * sign and the resolution function of the jet's impact parameter significance. This probability
 * helps in distinguishing between tracks likely originating from the primary vertex and those from
 * secondary vertices, aiding in jet flavor tagging.
 *
 * @param fResoFuncjet The resolution function for the jet, used to model the distribution of impact
 *                     parameter significances for tracks associated with the jet.
 * @param track The track for which the probability is being calculated.
 * @param minSignImpXYSig The minimum significance of the impact parameter in the XY plane, used as
 *                        the lower limit for integration of the resolution function. Defaults to -40.
 * @return The calculated probability of the track being associated with the jet, based on its
 *         impact parameter significance.
 */
template <typename T, typename U>
float getTrackProbability(T const& fResoFuncjet, U const& track, float minSignImpXYSig = -40)
{
  float probTrack = 0;
  auto varSignImpXYSig = std::abs(track.dcaXY()) / track.sigmadcaXY();
  if (-varSignImpXYSig < minSignImpXYSig)
    varSignImpXYSig = -minSignImpXYSig - 0.01; // To avoid overflow for integral
  probTrack = fResoFuncjet->Integral(minSignImpXYSig, -varSignImpXYSig) / fResoFuncjet->Integral(minSignImpXYSig, 0);

  return probTrack;
}

/**
 * Computes the jet probability (JP) for a given jet, considering only tracks with a positive geometric
 * sign. JP is calculated using the product of individual track probabilities and the sum of logarithmic
 * terms derived from these probabilities, providing a measure for the likelihood of the jet being
 * associated with a particular flavor based on its constituent tracks' impact parameters.
 *
 * @param fResoFuncjet: The resolution function for the jet, applied to each track within the jet to
 *                     assess its probability based on the impact parameter significance.
 * @param collision: The collision event data, necessary for geometric sign calculations.
 * @param jet: The jet for which the probability is being calculated.
 * @param tracks: Tracks in jets
 * @param cnt: ordering number of impact parameter cnt=0: untagged, cnt=1: first, cnt=2: seconde, cnt=3: third.
 * @param tagPoint: tagging working point which is selected by condiered efficiency and puriy
 * @param minSignImpXYSig: To avoid over fitting
 * @return The jet probability (JP), indicating the likelihood of the jet's association with a
 *         specific flavor. Returns -1 if the jet contains fewer than two tracks with a positive
 *         geometric sign.
 */
template <typename T, typename U, typename V>
float getJetProbability(T const& fResoFuncjet, U const& jet, V const& /*tracks*/, float trackDcaXYMax, float trackDcaZMax, float minSignImpXYSig = -10)
{
  std::vector<float> jetTracksPt;
  float trackjetProb = 1.;

  for (auto const& track : jet.template tracks_as<V>()) {
    if (!trackAcceptanceWithDca(track, trackDcaXYMax, trackDcaZMax))
      continue;

    float probTrack = getTrackProbability(fResoFuncjet, track, minSignImpXYSig);

    auto geoSign = getGeoSign(jet, track);
    if (geoSign > 0) { // only take positive sign track for JP calculation
      trackjetProb *= probTrack;
      jetTracksPt.push_back(track.pt());
    }
  }

  float jetProb = -1.;
  if (jetTracksPt.size() < 2)
    return -1;

  float sumjetProb = 0.;
  for (std::vector<float>::size_type i = 0; i < jetTracksPt.size(); i++) {
    sumjetProb += (std::pow(-1 * std::log(trackjetProb), static_cast<int>(i)) / TMath::Factorial(i));
  }

  jetProb = trackjetProb * sumjetProb;
  return jetProb;
}

// overloading for the case of using resolution function for each pt range
template <typename T, typename U, typename V>
float getJetProbability(std::vector<std::unique_ptr<T>> const& fResoFuncjets, U const& jet, V const& /*tracks*/, float trackDcaXYMax, float trackDcaZMax, float minSignImpXYSig = -10)
{
  std::vector<float> jetTracksPt;
  float trackjetProb = 1.;

  for (auto const& track : jet.template tracks_as<V>()) {
    if (!trackAcceptanceWithDca(track, trackDcaXYMax, trackDcaZMax))
      continue;

    float probTrack = -1;
    // choose the proper resolution function for the track based on its pt.
    if (track.pt() >= 0.0 && track.pt() < 0.5) {
      probTrack = getTrackProbability(fResoFuncjets.at(0), track, minSignImpXYSig);
    } else if (track.pt() >= 0.5 && track.pt() < 1.0) {
      probTrack = getTrackProbability(fResoFuncjets.at(1), track, minSignImpXYSig);
    } else if (track.pt() >= 1.0 && track.pt() < 2.0) {
      probTrack = getTrackProbability(fResoFuncjets.at(2), track, minSignImpXYSig);
    } else if (track.pt() >= 2.0 && track.pt() < 4.0) {
      probTrack = getTrackProbability(fResoFuncjets.at(3), track, minSignImpXYSig);
    } else if (track.pt() >= 4.0 && track.pt() < 6.0) {
      probTrack = getTrackProbability(fResoFuncjets.at(4), track, minSignImpXYSig);
    } else if (track.pt() >= 6.0 && track.pt() < 9.0) {
      probTrack = getTrackProbability(fResoFuncjets.at(5), track, minSignImpXYSig);
    } else if (track.pt() >= 9.0) {
      probTrack = getTrackProbability(fResoFuncjets.at(6), track, minSignImpXYSig);
    }

    auto geoSign = getGeoSign(jet, track);
    if (geoSign > 0) { // only take positive sign track for JP calculation
      trackjetProb *= probTrack;
      jetTracksPt.push_back(track.pt());
    }
  }

  float jetProb = -1.;
  if (jetTracksPt.size() < 2)
    return -1;

  float sumjetProb = 0.;
  for (std::vector<float>::size_type i = 0; i < jetTracksPt.size(); i++) {
    sumjetProb += (std::pow(-1 * std::log(trackjetProb), static_cast<int>(i)) / TMath::Factorial(i));
  }

  jetProb = trackjetProb * sumjetProb;
  return jetProb;
}

// For secaondy vertex method utilites
template <typename ProngType, typename JetType>
typename ProngType::iterator jetFromProngMaxDecayLength(const JetType& jet, float const& prongChi2PCAMin, float prongChi2PCAMax, float prongsigmaLxyMax, float prongIPxyMin, float prongIPxyMax, bool doXYZ = false, bool* checkSv = nullptr)
{
  if (checkSv)
    *checkSv = false;
  float maxSxy = -1.0f;
  typename ProngType::iterator bjetCand;
  for (const auto& prong : jet.template secondaryVertices_as<ProngType>()) {
    if (!prongAcceptance(prong, prongChi2PCAMin, prongChi2PCAMax, prongsigmaLxyMax, prongIPxyMin, prongIPxyMax, doXYZ))
      continue;
    *checkSv = true;
    float sxy = -1.0f;
    if (!doXYZ) {
      sxy = prong.decayLengthXY() / prong.errorDecayLengthXY();
    } else {
      sxy = prong.decayLength() / prong.errorDecayLength();
    }
    if (maxSxy < sxy) {
      maxSxy = sxy;
      bjetCand = prong;
    }
  }
  return bjetCand;
}

template <typename T, typename U>
bool isTaggedJetSV(T const jet, U const& /*prongs*/, float prongChi2PCAMin, float prongChi2PCAMax, float prongsigmaLxyMax, float prongIPxyMin, float prongIPxyMax, float svDispersionMax, float doXYZ = false, float tagPointForSV = 15.)
{
  bool checkSv = false;
  auto bjetCand = jetFromProngMaxDecayLength<U>(jet, prongChi2PCAMin, prongChi2PCAMax, prongsigmaLxyMax, prongIPxyMin, prongIPxyMax, doXYZ, &checkSv);
  if (!(checkSv && svAcceptance(bjetCand, svDispersionMax)))
    return false;
  if (!doXYZ) {
    auto maxSxy = bjetCand.decayLengthXY() / bjetCand.errorDecayLengthXY();
    if (maxSxy < tagPointForSV)
      return false;
  } else {
    auto maxSxyz = bjetCand.decayLength() / bjetCand.errorDecayLength();
    if (maxSxyz < tagPointForSV)
      return false;
  }
  return true;
}

template <typename T, typename U, typename V = float>
uint16_t setTaggingIPBit(T const& jet, U const& tracks, V trackDcaXYMax, V trackDcaZMax, V tagPointForIP)
{
  uint16_t bit = 0;
  auto [taggedIPsN1, taggedIPsN2, taggedIPsN3] = isGreaterThanTaggingPoint(jet, tracks, trackDcaXYMax, trackDcaZMax, tagPointForIP, false);
  if (taggedIPsN1) {
    SETBIT(bit, BJetTaggingMethod::IPsN1);
  }
  if (taggedIPsN2) {
    SETBIT(bit, BJetTaggingMethod::IPsN2);
  }
  if (taggedIPsN3) {
    SETBIT(bit, BJetTaggingMethod::IPsN3);
  }

  auto [taggedIPs3DN1, taggedIPs3DN2, taggedIPs3DN3] = isGreaterThanTaggingPoint(jet, tracks, trackDcaXYMax, trackDcaZMax, tagPointForIP, true);
  if (taggedIPs3DN1) {
    SETBIT(bit, BJetTaggingMethod::IPs3DN1);
  }
  if (taggedIPs3DN2) {
    SETBIT(bit, BJetTaggingMethod::IPs3DN2);
  }
  if (taggedIPs3DN3) {
    SETBIT(bit, BJetTaggingMethod::IPs3DN3);
  }
  return bit;
}

template <typename T, typename U, typename V = float>
uint16_t setTaggingSVBit(T const& jet, U const& prongs, V prongChi2PCAMin, V prongChi2PCAMax, V prongsigmaLxyMax, float prongIPxyMin, float prongIPxyMax, V svDispersionMax, V tagPointForSV)
{
  uint16_t bit = 0;
  if (isTaggedJetSV(jet, prongs, prongChi2PCAMin, prongChi2PCAMax, prongsigmaLxyMax, prongIPxyMin, prongIPxyMax, svDispersionMax, false, tagPointForSV)) {
    SETBIT(bit, BJetTaggingMethod::SV);
  }
  if (isTaggedJetSV(jet, prongs, prongChi2PCAMin, prongChi2PCAMax, prongsigmaLxyMax, prongIPxyMin, prongIPxyMax, svDispersionMax, true, tagPointForSV)) {
    SETBIT(bit, BJetTaggingMethod::SV3D);
  }
  return bit;
}

/**
 * Clusters jet constituent tracks into groups of tracks originating from same mcParticle position (trkVtxIndex), and finds each track origin (trkOrigin). (for GNN b-jet tagging)
 * @param trkLabels Track labels for GNN vertex and track origin predictions. trkVtxIndex: The index value of each vertex (cluster) which is determined by the function. trkOrigin: The category of the track origin (0: not physical primary, 1: charm, 2: beauty, 3: primary vertex, 4: other secondary vertex).
 * @param vtxResParam Vertex resolution parameter which determines the cluster size. (cm)
 * @param trackPtMin Minimum value of track pT.
 * @return The number of vertices (clusters) in the jet.
 */
template <typename AnyCollision, typename AnalysisJet, typename AnyTracks, typename AnyParticles, typename AnyOriginalParticles>
int vertexClustering(AnyCollision const& collision, AnalysisJet const& jet, AnyTracks const&, AnyParticles const& particles, AnyOriginalParticles const&, std::unordered_map<std::string, std::vector<int>>& trkLabels, bool searchUpToQuark, float vtxResParam = 0.01 /* 0.01cm = 100um */, float trackPtMin = 0.5)
{
  const auto& tracks = jet.template tracks_as<AnyTracks>();
  const int nTrks = tracks.size();

  // trkVtxIndex

  std::vector<int> tempTrkVtxIndex;

  int i = 0;
  for (const auto& constituent : tracks) {
    if (!constituent.has_mcParticle() || !constituent.template mcParticle_as<AnyParticles>().isPhysicalPrimary() || constituent.pt() < trackPtMin)
      tempTrkVtxIndex.push_back(-1);
    else
      tempTrkVtxIndex.push_back(i++);
  }
  tempTrkVtxIndex.push_back(i); // temporary index for PV
  if (nTrks < 1) {              // the process should be done for nTrks == 1 as well
    trkLabels["trkVtxIndex"] = tempTrkVtxIndex;
    return nTrks;
  }

  int nPos = nTrks + 1;
  std::vector<float> dists(nPos * (nPos - 1) / 2);
  auto trkPairIdx = [nPos](int ti, int tj) {
    if (ti == tj || ti >= nPos || tj >= nPos || ti < 0 || tj < 0) {
      LOGF(info, "Track pair index out of range");
      return -1;
    } else {
      return (ti < tj) ? (ti * nPos - (ti * (ti + 1)) / 2 + tj - ti - 1) : (tj * nPos - (tj * (tj + 1)) / 2 + ti - tj - 1);
    }
  }; // index nTrks is for PV

  for (int ti = 0; ti < nPos - 1; ti++)
    for (int tj = ti + 1; tj < nPos; tj++) {
      std::array<float, 3> posi, posj;

      if (tj < nTrks) {
        if (tracks[tj].has_mcParticle()) {
          const auto& pj = tracks[tj].template mcParticle_as<AnyParticles>().template mcParticle_as<AnyOriginalParticles>();
          posj = std::array<float, 3>{pj.vx(), pj.vy(), pj.vz()};
        } else {
          dists[trkPairIdx(ti, tj)] = std::numeric_limits<float>::max();
          continue;
        }
      } else {
        posj = std::array<float, 3>{collision.posX(), collision.posY(), collision.posZ()};
      }

      if (tracks[ti].has_mcParticle()) {
        const auto& pi = tracks[ti].template mcParticle_as<AnyParticles>().template mcParticle_as<AnyOriginalParticles>();
        posi = std::array<float, 3>{pi.vx(), pi.vy(), pi.vz()};
      } else {
        dists[trkPairIdx(ti, tj)] = std::numeric_limits<float>::max();
        continue;
      }

      dists[trkPairIdx(ti, tj)] = RecoDecay::distance(posi, posj);
    }

  int clusteri = -1, clusterj = -1;
  float minMinDist = -1.f; // If there is an not-merge-able minDist pair, check the 2nd-minDist pair.
  while (true) {

    float minDist = -1.f; // Get minDist pair
    for (int ti = 0; ti < nPos - 1; ti++)
      for (int tj = ti + 1; tj < nPos; tj++)
        if (tempTrkVtxIndex[ti] != tempTrkVtxIndex[tj] && tempTrkVtxIndex[ti] >= 0 && tempTrkVtxIndex[tj] >= 0) {
          float dist = dists[trkPairIdx(ti, tj)];
          if ((dist < minDist || minDist < 0.f) && dist > minMinDist) {
            minDist = dist;
            clusteri = ti;
            clusterj = tj;
          }
        }
    if (clusteri < 0 || clusterj < 0)
      break;

    bool mrg = true; // Merge-ability check
    for (int ti = 0; ti < nPos && mrg; ti++)
      if (tempTrkVtxIndex[ti] == tempTrkVtxIndex[clusteri] && tempTrkVtxIndex[ti] >= 0) {
        for (int tj = 0; tj < nPos && mrg; tj++)
          if (tj != ti && tempTrkVtxIndex[tj] == tempTrkVtxIndex[clusterj] && tempTrkVtxIndex[tj] >= 0) {
            if (dists[trkPairIdx(ti, tj)] > vtxResParam) { // If there is more distant pair compared to vtx_res between two clusters, they cannot be merged.
              mrg = false;
              minMinDist = minDist;
            }
          }
      }
    if (minDist > vtxResParam || minDist < 0.f)
      break;

    if (mrg) { // Merge two clusters
      int oldIndex = tempTrkVtxIndex[clusterj];
      for (int t = 0; t < nPos; t++)
        if (tempTrkVtxIndex[t] == oldIndex)
          tempTrkVtxIndex[t] = tempTrkVtxIndex[clusteri];
    }
  }

  int nVertices = 0;

  // Sort the indices from PV (as 0) to the most distant SV (as 1~).
  int idxPV = tempTrkVtxIndex[nTrks];
  for (int t = 0; t < nTrks; t++)
    if (tempTrkVtxIndex[t] == idxPV) {
      tempTrkVtxIndex[t] = -2;
      nVertices = 1; // There is a track originating from PV
    }

  std::unordered_map<int, float> avgDistances;
  std::unordered_map<int, int> count;
  for (int t = 0; t < nTrks; t++) {
    if (tempTrkVtxIndex[t] >= 0) {
      avgDistances[tempTrkVtxIndex[t]] += dists[trkPairIdx(t, nTrks)];
      count[tempTrkVtxIndex[t]]++;
    }
  }

  trkLabels["trkVtxIndex"] = std::vector<int>(nTrks, -1);
  if (count.size() != 0) {                        // If there is any SV cluster not only PV cluster
    for (auto& [idx, avgDistance] : avgDistances) // o2-linter: disable=const-ref-in-for-loop
      avgDistance /= count[idx];

    nVertices += avgDistances.size();

    std::vector<std::pair<int, float>> sortedIndices(avgDistances.begin(), avgDistances.end());
    std::sort(sortedIndices.begin(), sortedIndices.end(), [](const auto& a, const auto& b) { return a.second < b.second; });
    int rank = 1;
    for (auto const& [idx, avgDistance] : sortedIndices) {
      bool found = false;
      for (int t = 0; t < nTrks; t++)
        if (tempTrkVtxIndex[t] == idx) {
          trkLabels["trkVtxIndex"][t] = rank;
          found = true;
        }
      rank += found;
    }
  }

  for (int t = 0; t < nTrks; t++)
    if (tempTrkVtxIndex[t] == -2)
      trkLabels["trkVtxIndex"][t] = 0;

  // trkOrigin

  int trkIdx = 0;
  for (auto const& constituent : jet.template tracks_as<AnyTracks>()) {
    if (!constituent.has_mcParticle() || !constituent.template mcParticle_as<AnyParticles>().isPhysicalPrimary() || constituent.pt() < trackPtMin) {
      trkLabels["trkOrigin"].push_back(0);
    } else {
      const auto& particle = constituent.template mcParticle_as<AnyParticles>();
      int orig = RecoDecay::getParticleOrigin(particles, particle, searchUpToQuark);
      trkLabels["trkOrigin"].push_back((orig > 0) ? orig : (trkLabels["trkVtxIndex"][trkIdx] == 0) ? 3
                                                                                                   : 4);
    }

    trkIdx++;
  }

  return nVertices;
}

std::vector<std::vector<float>> getInputsForML(BJetParams jetparams, std::vector<BJetTrackParams>& tracksParams, std::vector<BJetSVParams>& svsParams, int maxJetConst = 10)
{
  std::vector<float> jetInput = {jetparams.mJetpT, jetparams.mJetEta, jetparams.mJetPhi, static_cast<float>(jetparams.mNTracks), static_cast<float>(jetparams.mNSV), jetparams.mJetMass};
  std::vector<float> tracksInputFlat;
  std::vector<float> svsInputFlat;

  for (int iconstit = 0; iconstit < maxJetConst; iconstit++) {

    tracksInputFlat.push_back(tracksParams[iconstit].mTrackpT);
    tracksInputFlat.push_back(tracksParams[iconstit].mTrackEta);
    tracksInputFlat.push_back(tracksParams[iconstit].mDotProdTrackJet);
    tracksInputFlat.push_back(tracksParams[iconstit].mDotProdTrackJetOverJet);
    tracksInputFlat.push_back(tracksParams[iconstit].mDeltaRJetTrack);
    tracksInputFlat.push_back(tracksParams[iconstit].mSignedIP2D);
    tracksInputFlat.push_back(tracksParams[iconstit].mSignedIP2DSign);
    tracksInputFlat.push_back(tracksParams[iconstit].mSignedIP3D);
    tracksInputFlat.push_back(tracksParams[iconstit].mSignedIP3DSign);
    tracksInputFlat.push_back(tracksParams[iconstit].mMomFraction);
    tracksInputFlat.push_back(tracksParams[iconstit].mDeltaRTrackVertex);

    svsInputFlat.push_back(svsParams[iconstit].mSVpT);
    svsInputFlat.push_back(svsParams[iconstit].mDeltaRSVJet);
    svsInputFlat.push_back(svsParams[iconstit].mSVMass);
    svsInputFlat.push_back(svsParams[iconstit].mSVfE);
    svsInputFlat.push_back(svsParams[iconstit].mIPXY);
    svsInputFlat.push_back(svsParams[iconstit].mCPA);
    svsInputFlat.push_back(svsParams[iconstit].mChi2PCA);
    svsInputFlat.push_back(svsParams[iconstit].mDispersion);
    svsInputFlat.push_back(svsParams[iconstit].mDecayLength2D);
    svsInputFlat.push_back(svsParams[iconstit].mDecayLength2DError);
    svsInputFlat.push_back(svsParams[iconstit].mDecayLength3D);
    svsInputFlat.push_back(svsParams[iconstit].mDecayLength3DError);
  }

  std::vector<std::vector<float>> totalInput;
  totalInput.push_back(jetInput);
  totalInput.push_back(tracksInputFlat);
  totalInput.push_back(svsInputFlat);

  return totalInput;
}

// Looping over the SV info and putting them in the input vector
template <typename AnalysisJet, typename AnyTracks, typename SecondaryVertices>
void analyzeJetSVInfo4ML(AnalysisJet const& myJet, AnyTracks const& /*allTracks*/, SecondaryVertices const& /*allSVs*/, std::vector<BJetSVParams>& svsParams, float svPtMin = 1.0, int svReductionFactor = 3)
{
  using SVType = typename SecondaryVertices::iterator;

  // Min-heap to store the top 30 SVs by decayLengthXY/errorDecayLengthXY
  auto compare = [](SVType& sv1, SVType& sv2) {
    return (sv1.decayLengthXY() / sv1.errorDecayLengthXY()) > (sv2.decayLengthXY() / sv2.errorDecayLengthXY());
  };

  auto svs = myJet.template secondaryVertices_as<SecondaryVertices>();

  // Sort the SVs based on their decay length significance in descending order
  // This is needed in order to select longest SVs since some jets could have thousands of SVs
  std::sort(svs.begin(), svs.end(), compare);

  for (const auto& candSV : svs) {

    if (candSV.pt() < svPtMin) {
      continue;
    }

    double deltaRJetSV = jetutilities::deltaR(myJet, candSV);
    double massSV = candSV.m();
    double energySV = candSV.e();

    if (svsParams.size() < (svReductionFactor * myJet.template tracks_as<AnyTracks>().size())) {
      svsParams.emplace_back(BJetSVParams{candSV.pt(), deltaRJetSV, massSV, energySV / myJet.energy(), candSV.impactParameterXY(), candSV.cpa(), candSV.chi2PCA(), candSV.dispersion(), candSV.decayLengthXY(), candSV.errorDecayLengthXY(), candSV.decayLength(), candSV.errorDecayLength()});
    }
  }
}

// Looping over the track info and putting them in the input vector
template <typename AnalysisJet, typename AnyTracks, typename SecondaryVertices>
void analyzeJetTrackInfo4ML(AnalysisJet const& analysisJet, AnyTracks const& /*allTracks*/, SecondaryVertices const& /*allSVs*/, std::vector<BJetTrackParams>& tracksParams, float trackPtMin = 0.5)
{
  for (const auto& constituent : analysisJet.template tracks_as<AnyTracks>()) {

    if (constituent.pt() < trackPtMin) {
      continue;
    }

    double deltaRJetTrack = jetutilities::deltaR(analysisJet, constituent);
    double dotProduct = RecoDecay::dotProd(std::array<float, 3>{analysisJet.px(), analysisJet.py(), analysisJet.pz()}, std::array<float, 3>{constituent.px(), constituent.py(), constituent.pz()});
    int sign = jettaggingutilities::getGeoSign(analysisJet, constituent);

    float rClosestSV = 10.;
    for (const auto& candSV : analysisJet.template secondaryVertices_as<SecondaryVertices>()) {
      double deltaRTrackSV = jetutilities::deltaR(constituent, candSV);
      if (deltaRTrackSV < rClosestSV) {
        rClosestSV = deltaRTrackSV;
      }
    }

    tracksParams.emplace_back(BJetTrackParams{constituent.pt(), constituent.eta(), dotProduct, dotProduct / analysisJet.p(), deltaRJetTrack, std::abs(constituent.dcaXY()) * sign, constituent.sigmadcaXY(), std::abs(constituent.dcaXYZ()) * sign, constituent.sigmadcaXYZ(), constituent.p() / analysisJet.p(), rClosestSV});
  }

  auto compare = [](BJetTrackParams& tr1, BJetTrackParams& tr2) {
    return (tr1.mSignedIP2D / tr1.mSignedIP2DSign) > (tr2.mSignedIP2D / tr2.mSignedIP2DSign);
  };

  // Sort the tracks based on their IP significance in descending order
  std::sort(tracksParams.begin(), tracksParams.end(), compare);
}

// Looping over the track info and putting them in the input vector (for GNN b-jet tagging)
template <typename AnalysisJet, typename AnyTracks, typename AnyOriginalTracks>
void analyzeJetTrackInfo4GNN(AnalysisJet const& analysisJet, AnyTracks const& /*allTracks*/, AnyOriginalTracks const& /*origTracks*/, std::vector<std::vector<float>>& tracksParams, float trackPtMin = 0.5, int64_t nMaxConstit = 40)
{
  for (const auto& constituent : analysisJet.template tracks_as<AnyTracks>()) {

    if (constituent.pt() < trackPtMin) {
      continue;
    }

    int sign = jettaggingutilities::getGeoSign(analysisJet, constituent);

    auto origConstit = constituent.template track_as<AnyOriginalTracks>();

    if (static_cast<int64_t>(tracksParams.size()) < nMaxConstit) {
      tracksParams.emplace_back(std::vector<float>{constituent.pt(), origConstit.phi(), constituent.eta(), static_cast<float>(constituent.sign()), std::abs(constituent.dcaXY()) * sign, constituent.sigmadcaXY(), std::abs(constituent.dcaXYZ()) * sign, constituent.sigmadcaXYZ(), static_cast<float>(origConstit.itsNCls()), static_cast<float>(origConstit.tpcNClsFound()), static_cast<float>(origConstit.tpcNClsCrossedRows()), origConstit.itsChi2NCl(), origConstit.tpcChi2NCl()});
    } else {
      // If there are more than nMaxConstit constituents in the jet, select only nMaxConstit constituents with the highest DCA_XY significance.
      size_t minIdx = 0;
      for (size_t i = 0; i < tracksParams.size(); ++i) {
        if (tracksParams[i][4] / tracksParams[i][5] < tracksParams[minIdx][4] / tracksParams[minIdx][5])
          minIdx = i;
      }
      if (std::abs(constituent.dcaXY()) * sign / constituent.sigmadcaXY() > tracksParams[minIdx][4] / tracksParams[minIdx][5])
        tracksParams[minIdx] = std::vector<float>{constituent.pt(), origConstit.phi(), constituent.eta(), static_cast<float>(constituent.sign()), std::abs(constituent.dcaXY()) * sign, constituent.sigmadcaXY(), std::abs(constituent.dcaXYZ()) * sign, constituent.sigmadcaXYZ(), static_cast<float>(origConstit.itsNCls()), static_cast<float>(origConstit.tpcNClsFound()), static_cast<float>(origConstit.tpcNClsCrossedRows()), origConstit.itsChi2NCl(), origConstit.tpcChi2NCl()};
    }
  }
}

// Discriminant value for GNN b-jet tagging
template <typename T>
T Db(const std::vector<T>& logits, double fC = 0.018)
{
  auto softmax = [](const std::vector<T>& logits) {
    std::vector<T> res;
    T maxLogit = *std::max_element(logits.begin(), logits.end());
    T sumLogit = 0.;
    for (size_t i = 0; i < logits.size(); ++i) {
      res.push_back(std::exp(logits[i] - maxLogit));
      sumLogit += res[i];
    }
    for (size_t i = 0; i < logits.size(); ++i) {
      res[i] /= sumLogit;
    }
    return res;
  };

  std::vector<T> softmaxLogits = softmax(logits);

  if (softmaxLogits[1] == 0. && softmaxLogits[2] == 0.) {
    LOG(debug) << "jettaggingutilities::Db, Divide by zero: softmaxLogits = (" << softmaxLogits[0] << ", " << softmaxLogits[1] << ", " << softmaxLogits[2] << ")";
  }

  return std::log(softmaxLogits[0] / (fC * softmaxLogits[1] + (1. - fC) * softmaxLogits[2]));
}

}; // namespace jettaggingutilities

#endif // PWGJE_CORE_JETTAGGINGUTILITIES_H_
