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
// O2 includes

/// \file HFFilterHelpers.h
/// \brief Header file with definition of variables, methods, and tables used in the HFFilter.cxx task
///
/// \author Fabrizio Grosa <fabrizio.grosa@cern.ch>, CERN
/// \author Marcel Lesch <marcel.lesch@tum.de>, TUM
/// \author Alexandre Bigot <alexandre.bigot@cern.ch>, Strasbourg University
/// \author Biao Zhang <biao.zhang@cern.ch>, CCNU

#ifndef O2_ANALYSIS_HF_FILTER_HELPERS_H_
#define O2_ANALYSIS_HF_FILTER_HELPERS_H_

#include "Framework/DataTypes.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/Core/RecoDecay.h"

#include <vector>
#include <array>
#include <string>
#include <cmath>

#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "Math/GenVector/Boost.h"

// CCDB
#include "CCDB/CcdbApi.h"

namespace o2::aod
{

namespace hffilters
{

enum HfTriggers {
  kHighPt2P = 0,
  kHighPt3P,
  kBeauty3P,
  kBeauty4P,
  kFemto2P,
  kFemto3P,
  kDoubleCharm2P,
  kDoubleCharm3P,
  kDoubleCharmMix,
  kNtriggersHF
};

enum charmParticles {
  kD0 = 0,
  kDplus,
  kDs,
  kLc,
  kXic,
  kNCharmParticles
};

enum beautyParticles {
  kBplus = 0,
  kB0toDStar,
  kB0,
  kBs,
  kLb,
  kXib,
  kNBeautyParticles
};

enum beautyTrackSelection {
  kRejected = 0,
  kSoftPion,
  kRegular
};

static const std::array<std::string, kNtriggersHF> HfTriggerNames{"highPt", "beauty", "femto", "doubleCharm"};
static const std::array<std::string, kNCharmParticles> charmParticleNames{"D0", "Dplus", "Ds", "Lc", "Xic"};
static const std::array<std::string, kNBeautyParticles> beautyParticleNames{"Bplus", "B0toDStar", "B0", "Bs", "Lb", "Xib"};
static const std::array<int, kNCharmParticles> pdgCodesCharm{421, 411, 431, 4122, 4232};

static const std::tuple pdgCharmDaughters{
  std::array{-321, 211},        // D0
  std::array{-321, 211, 211},   // Dplus
  std::array{321, -321, 211},   // Ds
  std::array{2212, -321, 211},  // Lc
  std::array{2212, -321, 211}}; // Xic

static const float massPi = RecoDecay::getMassPDG(211);
static const float massK = RecoDecay::getMassPDG(321);
static const float massProton = RecoDecay::getMassPDG(2212);
static const float massPhi = RecoDecay::getMassPDG(333);
static const float massD0 = RecoDecay::getMassPDG(421);
static const float massDPlus = RecoDecay::getMassPDG(411);
static const float massDs = RecoDecay::getMassPDG(431);
static const float massLc = RecoDecay::getMassPDG(4122);
static const float massXic = RecoDecay::getMassPDG(4232);
static const float massDStar = RecoDecay::getMassPDG(413);
static const float massBPlus = RecoDecay::getMassPDG(511);
static const float massB0 = RecoDecay::getMassPDG(521);
static const float massBs = RecoDecay::getMassPDG(531);
static const float massLb = RecoDecay::getMassPDG(5122);
static const float massXib = RecoDecay::getMassPDG(5232);

/// Computation of the relative momentum between particle pairs
/// \param track is a track
/// \param ProtonMass is the mass of a proton
/// \param CharmCandMomentum is the three momentum of a charm candidate
/// \param CharmMass is the mass of the charm hadron
/// \return relative momentum of pair
template <typename T>
float computeRelativeMomentum(const T& track, const std::array<float, 3>& CharmCandMomentum, const float& CharmMass)
{
  ROOT::Math::PxPyPzMVector part1(track.px(), track.py(), track.pz(), massProton);
  ROOT::Math::PxPyPzMVector part2(CharmCandMomentum[0], CharmCandMomentum[1], CharmCandMomentum[2], CharmMass);

  ROOT::Math::PxPyPzMVector trackSum = part1 + part2;
  ROOT::Math::Boost boostv12{trackSum.BoostToCM()};
  ROOT::Math::PxPyPzMVector part1CM = boostv12(part1);
  ROOT::Math::PxPyPzMVector part2CM = boostv12(part2);
  ROOT::Math::PxPyPzMVector trackRelK = part1CM - part2CM;

  float kStar = 0.5 * trackRelK.P();
  return kStar;
} // float computeRelativeMomentum(const T& track, const std::array<float, 3>& CharmCandMomentum, const float& CharmMass)

/// Computation of the number of candidates in an event that do not share daughter tracks
/// \return 0 or 1 in case of less than 2 independent candidates in a single event, 2 otherwise
template <typename T>
int computeNumberOfCandidates(std::vector<std::vector<T>> indices)
{
  if (indices.size() < 2) {
    return indices.size();
  }

  std::vector<int> numIndependentCand{};
  for (auto iCand{0u}; iCand < indices.size(); ++iCand) {
    int nIndependent = 0;
    for (auto iCandSecond{0u}; iCandSecond < indices.size(); ++iCandSecond) {
      if (iCand == iCandSecond) {
        continue;
      } else {
        bool hasOverlap = false;
        for (auto idxFirst{0u}; idxFirst < indices[iCand].size(); ++idxFirst) {
          for (auto idxSecond{0u}; idxSecond < indices[iCandSecond].size(); ++idxSecond) {
            if (indices[iCand][idxFirst] == indices[iCandSecond][idxSecond]) {
              hasOverlap = true;
              break;
            }
          }
        }
        if (!hasOverlap) {
          nIndependent++;
        }
      }
    }
    numIndependentCand.push_back(nIndependent);
  }
  std::sort(numIndependentCand.begin(), numIndependentCand.end());

  if (numIndependentCand.back() == 0) {
    return numIndependentCand.back();
  }

  return 2;
}

/// ML helper methods

/// Iinitialisation of ONNX session
/// \param onnxFile is the onnx file name
/// \param partName is the particle name
/// \param inputNames is a vector of input names
/// \param inputShapes is a vector of input shapes
/// \param outputNames is a vector of output names
/// \param dataType is the data type (1=float, 11=double)
/// \param loadModelsFromCCDB is the flag to decide whether the ONNX file is read from CCDB or not
/// \param ccdbApi is the CCDB API
/// \param mlModelPathCCDB is the model path in CCDB
/// \param timestampCCDB is the CCDB timestamp
/// \return the ONNX Ort::Experimental::Session
std::shared_ptr<Ort::Experimental::Session> InitONNXSession(std::string& onnxFile, std::string partName, std::vector<std::string>& inputNames, std::vector<std::vector<int64_t>>& inputShapes, std::vector<std::string>& outputNames, int& dataType, bool loadModelsFromCCDB, o2::ccdb::CcdbApi& ccdbApi, std::string mlModelPathCCDB, long timestampCCDB)
{
  // hard coded, we do not let the user change this
  Ort::Env env{ORT_LOGGING_LEVEL_ERROR, Form("ml-model-%s-triggers", partName.data())};
  Ort::SessionOptions sessionOpt{Ort::SessionOptions()};
  std::shared_ptr<Ort::Experimental::Session> session = nullptr;
  sessionOpt.SetIntraOpNumThreads(1);
  sessionOpt.SetInterOpNumThreads(1);

  std::map<std::string, std::string> metadata;
  bool retrieveSuccess = true;
  if (loadModelsFromCCDB && timestampCCDB > 0) {
    retrieveSuccess = ccdbApi.retrieveBlob(mlModelPathCCDB + partName, ".", metadata, timestampCCDB, false, onnxFile);
  }
  if (retrieveSuccess) {
    session.reset(new Ort::Experimental::Session{env, onnxFile, sessionOpt});
    inputNames = session->GetInputNames();
    inputShapes = session->GetInputShapes();
    if (inputShapes[0][0] < 0) {
      LOGF(warning, Form("Model for %s with negative input shape likely because converted with ummingbird, setting it to 1.", partName.data()));
      inputShapes[0][0] = 1;
    }
    outputNames = session->GetOutputNames();

    Ort::TypeInfo typeInfo = session->GetInputTypeInfo(0);
    auto tensorInfo = typeInfo.GetTensorTypeAndShapeInfo();
    dataType = tensorInfo.GetElementType();
  } else {
    LOG(fatal) << "Error encountered while fetching/loading the network from CCDB! Maybe the network doesn't exist yet for this runnumber/timestamp?";
  }

  return session;
};

/// Iinitialisation of ONNX session
/// \param inputFeatures is the vector with input features
/// \param session is the ONNX Ort::Experimental::Session
/// \param inputNames is a vector of input names
/// \param inputShapes is a vector of input shapes
/// \param outputNames is a vector of output names
/// \return the array with the three output scores
template <typename T>
std::array<T, 3> PredictONNX(std::vector<T>& inputFeatures, std::shared_ptr<Ort::Experimental::Session>& session, std::vector<std::string>& inputNames, std::vector<std::vector<int64_t>>& inputShapes, std::vector<std::string>& outputNames)
{
  std::array<T, 3> scores{-1., 2., 2.};
  std::vector<Ort::Value> inputTensor;
  inputTensor.push_back(Ort::Experimental::Value::CreateTensor<T>(inputFeatures.data(), inputFeatures.size(), inputShapes[0]));

  // double-check the dimensions of the input tensor
  if (inputTensor[0].GetTensorTypeAndShapeInfo().GetShape()[0] > 0) { // vectorial models can have negative shape if the shape is unknown
    assert(inputTensor[0].IsTensor() && inputTensor[0].GetTensorTypeAndShapeInfo().GetShape() == inputShapes[0]);
  }
  try {
    auto outputTensor = session->Run(inputNames, inputTensor, outputNames);
    assert(outputTensor.size() == outputNames.size() && outputTensor[1].IsTensor());
    auto typeInfo = outputTensor[1].GetTensorTypeAndShapeInfo();
    assert(typeInfo.GetElementCount() == 3); // we need multiclass
    scores[0] = outputTensor[1].GetTensorMutableData<T>()[0];
    scores[1] = outputTensor[1].GetTensorMutableData<T>()[1];
    scores[2] = outputTensor[1].GetTensorMutableData<T>()[2];
  } catch (const Ort::Exception& exception) {
    // LOG(error) << "Error running model inference: " << exception.what();
  }

  return scores;
};
} // namespace hffilters

/// definition of tables

namespace extra2Prong
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);
} // namespace extra2Prong
namespace extra3Prong
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);
} // namespace extra3Prong
DECLARE_SOA_TABLE(Colls2Prong, "AOD", "COLLSID2P", o2::aod::extra2Prong::CollisionId);
DECLARE_SOA_TABLE(Colls3Prong, "AOD", "COLLSID3P", o2::aod::extra3Prong::CollisionId);

namespace hftraining
{
DECLARE_SOA_COLUMN(InvMassD0, invMassD0, float);                 //!
DECLARE_SOA_COLUMN(InvMassD0bar, invMassD0bar, float);           //!
DECLARE_SOA_COLUMN(InvMassDplus, invMassDplus, float);           //!
DECLARE_SOA_COLUMN(InvMassDsToKKPi, invMassDsToKKPi, float);     //!
DECLARE_SOA_COLUMN(InvMassDsToPiKK, invMassDsToPiKK, float);     //!
DECLARE_SOA_COLUMN(InvMassLcToPKPi, invMassLcToPKPi, float);     //!
DECLARE_SOA_COLUMN(InvMassLcToPiKP, invMassLcToPiKP, float);     //!
DECLARE_SOA_COLUMN(InvMassXicToPKPi, invMassXicToPKPi, float);   //!
DECLARE_SOA_COLUMN(InvMassXicToPiKP, invMassXicToPiKP, float);   //!
DECLARE_SOA_COLUMN(PT2Prong, pT2Prong, float);                   //!
DECLARE_SOA_COLUMN(PT3Prong, pT3Prong, float);                   //!
DECLARE_SOA_COLUMN(DeltaMassKKFirst, deltaMassKKFirst, float);   //!
DECLARE_SOA_COLUMN(DeltaMassKKSecond, deltaMassKKSecond, float); //!
DECLARE_SOA_COLUMN(PT1, pT1, float);                             //!
DECLARE_SOA_COLUMN(DCAPrimXY1, dcaPrimXY1, float);               //!
DECLARE_SOA_COLUMN(DCAPrimZ1, dcaPrimZ1, float);                 //!
DECLARE_SOA_COLUMN(NsigmaPiTPC1, nsigmaPiTPC1, float);           //!
DECLARE_SOA_COLUMN(NsigmaKaTPC1, nsigmaKaTPC1, float);           //!
DECLARE_SOA_COLUMN(NsigmaPrTPC1, nsigmaPrTPC1, float);           //!
DECLARE_SOA_COLUMN(NsigmaPiTOF1, nsigmaPiTOF1, float);           //!
DECLARE_SOA_COLUMN(NsigmaKaTOF1, nsigmaKaTOF1, float);           //!
DECLARE_SOA_COLUMN(NsigmaPrTOF1, nsigmaPrTOF1, float);           //!
DECLARE_SOA_COLUMN(PT2, pT2, float);                             //!
DECLARE_SOA_COLUMN(DCAPrimXY2, dcaPrimXY2, float);               //!
DECLARE_SOA_COLUMN(DCAPrimZ2, dcaPrimZ2, float);                 //!
DECLARE_SOA_COLUMN(NsigmaPiTPC2, nsigmaPiTPC2, float);           //!
DECLARE_SOA_COLUMN(NsigmaKaTPC2, nsigmaKaTPC2, float);           //!
DECLARE_SOA_COLUMN(NsigmaPrTPC2, nsigmaPrTPC2, float);           //!
DECLARE_SOA_COLUMN(NsigmaPiTOF2, nsigmaPiTOF2, float);           //!
DECLARE_SOA_COLUMN(NsigmaKaTOF2, nsigmaKaTOF2, float);           //!
DECLARE_SOA_COLUMN(NsigmaPrTOF2, nsigmaPrTOF2, float);           //!
DECLARE_SOA_COLUMN(PT3, pT3, float);                             //!
DECLARE_SOA_COLUMN(DCAPrimXY3, dcaPrimXY3, float);               //!
DECLARE_SOA_COLUMN(DCAPrimZ3, dcaPrimZ3, float);                 //!
DECLARE_SOA_COLUMN(NsigmaPiTPC3, nsigmaPiTPC3, float);           //!
DECLARE_SOA_COLUMN(NsigmaKaTPC3, nsigmaKaTPC3, float);           //!
DECLARE_SOA_COLUMN(NsigmaPrTPC3, nsigmaPrTPC3, float);           //!
DECLARE_SOA_COLUMN(NsigmaPiTOF3, nsigmaPiTOF3, float);           //!
DECLARE_SOA_COLUMN(NsigmaKaTOF3, nsigmaKaTOF3, float);           //!
DECLARE_SOA_COLUMN(NsigmaPrTOF3, nsigmaPrTOF3, float);           //!
DECLARE_SOA_COLUMN(FlagOrigin, flagOrigin, int8_t);              //!
DECLARE_SOA_COLUMN(Channel, channel, int8_t);                    //!
DECLARE_SOA_COLUMN(HFSelBit, hfselbit, int8_t);                  //!
} // namespace hftraining

DECLARE_SOA_TABLE(HFTrigTrain2P, "AOD", "HFTRIGTRAIN2P", //!
                  hftraining::InvMassD0,
                  hftraining::InvMassD0bar,
                  hftraining::PT2Prong,
                  hftraining::PT1,
                  hftraining::DCAPrimXY1,
                  hftraining::DCAPrimZ1,
                  hftraining::NsigmaPiTPC1,
                  hftraining::NsigmaKaTPC1,
                  hftraining::NsigmaPiTOF1,
                  hftraining::NsigmaKaTOF1,
                  hftraining::PT2,
                  hftraining::DCAPrimXY2,
                  hftraining::DCAPrimZ2,
                  hftraining::NsigmaPiTPC2,
                  hftraining::NsigmaKaTPC2,
                  hftraining::NsigmaPiTOF2,
                  hftraining::NsigmaKaTOF2,
                  hftraining::FlagOrigin);
DECLARE_SOA_TABLE(HFTrigTrain3P, "AOD", "HFTRIGTRAIN3P", //!
                  hftraining::InvMassDplus,
                  hftraining::InvMassDsToKKPi,
                  hftraining::InvMassDsToPiKK,
                  hftraining::InvMassLcToPKPi,
                  hftraining::InvMassLcToPiKP,
                  hftraining::InvMassXicToPKPi,
                  hftraining::InvMassXicToPiKP,
                  hftraining::PT3Prong,
                  hftraining::DeltaMassKKFirst,
                  hftraining::DeltaMassKKSecond,
                  hftraining::PT1,
                  hftraining::DCAPrimXY1,
                  hftraining::DCAPrimZ1,
                  hftraining::NsigmaPiTPC1,
                  hftraining::NsigmaKaTPC1,
                  hftraining::NsigmaPrTPC1,
                  hftraining::NsigmaPiTOF1,
                  hftraining::NsigmaKaTOF1,
                  hftraining::NsigmaPrTOF1,
                  hftraining::PT2,
                  hftraining::DCAPrimXY2,
                  hftraining::DCAPrimZ2,
                  hftraining::NsigmaPiTPC2,
                  hftraining::NsigmaKaTPC2,
                  hftraining::NsigmaPrTPC2,
                  hftraining::NsigmaPiTOF2,
                  hftraining::NsigmaKaTOF2,
                  hftraining::NsigmaPrTOF2,
                  hftraining::PT3,
                  hftraining::DCAPrimXY3,
                  hftraining::DCAPrimZ3,
                  hftraining::NsigmaPiTPC3,
                  hftraining::NsigmaKaTPC3,
                  hftraining::NsigmaPrTPC3,
                  hftraining::NsigmaPiTOF3,
                  hftraining::NsigmaKaTOF3,
                  hftraining::NsigmaPrTOF3,
                  hftraining::FlagOrigin,
                  hftraining::Channel,
                  hftraining::HFSelBit);

namespace hfoptimisationTree
{
DECLARE_SOA_COLUMN(CollisionIndex, collisionIndex, int); //!
DECLARE_SOA_COLUMN(ParticleID, particleID, int);         //!
DECLARE_SOA_COLUMN(Pt, pt, float);                       //!
DECLARE_SOA_COLUMN(BkgBDT, bkgBDT, float);               //!
DECLARE_SOA_COLUMN(PromptBDT, promptBDT, float);         //!
DECLARE_SOA_COLUMN(NonpromptBDT, nonpromptBDT, float);   //!
DECLARE_SOA_COLUMN(DCAXY, dcaXY, float);                 //!
DECLARE_SOA_COLUMN(KStar, kStar, float);                 //!
DECLARE_SOA_COLUMN(NsigmaPrTPC, nsigmaPrTPC, float);     //!
DECLARE_SOA_COLUMN(NsigmaPrTOF, nsigmaPrTOF, float);     //!
} // namespace hfoptimisationTree

DECLARE_SOA_TABLE(HFOptimisationTreeBeauty, "AOD", "HFOPTIMTREEB", //!
                  hfoptimisationTree::CollisionIndex,
                  hfoptimisationTree::ParticleID,
                  hfoptimisationTree::Pt,
                  hfoptimisationTree::BkgBDT,
                  hfoptimisationTree::PromptBDT,
                  hfoptimisationTree::NonpromptBDT,
                  hfoptimisationTree::DCAXY);
DECLARE_SOA_TABLE(HFOptimisationTreeCharm, "AOD", "HFOPTIMTREEC", //!
                  hfoptimisationTree::CollisionIndex,
                  hfoptimisationTree::ParticleID,
                  hfoptimisationTree::Pt,
                  hfoptimisationTree::BkgBDT,
                  hfoptimisationTree::PromptBDT,
                  hfoptimisationTree::NonpromptBDT);
DECLARE_SOA_TABLE(HFOptimisationTreeFemto, "AOD", "HFOPTIMTREEF", //!
                  hfoptimisationTree::CollisionIndex,
                  hfoptimisationTree::ParticleID,
                  hfoptimisationTree::Pt,
                  hfoptimisationTree::BkgBDT,
                  hfoptimisationTree::PromptBDT,
                  hfoptimisationTree::NonpromptBDT,
                  hfoptimisationTree::KStar,
                  hfoptimisationTree::NsigmaPrTPC,
                  hfoptimisationTree::NsigmaPrTOF);
DECLARE_SOA_TABLE(HFOptimisationTreeCollisions, "AOD", "HFOPTIMTREECOLL", //!
                  hfoptimisationTree::CollisionIndex)
} // namespace o2::aod

#endif // O2_ANALYSIS_HF_FILTER_HELPERS_