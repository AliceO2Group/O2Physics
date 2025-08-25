

// Copyright 2019-2025 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file femtoDreamResoSelection.h
/// \brief Definition of the FemtoDreamResoSelection
/// \author Christopher Klumm, TU München, christopher.klumm@cern.ch
/// \author Nils Fabian Konert, TU München, nils.fabian.konert@cern.ch

#ifndef PWGCF_FEMTODREAM_CORE_FEMTODREAMRESOSELECTION_H_
#define PWGCF_FEMTODREAM_CORE_FEMTODREAMRESOSELECTION_H_

#include "PWGCF/DataModel/FemtoDerived.h"
#include "PWGCF/FemtoDream/Core/femtoDreamObjectSelection.h"
#include "PWGCF/FemtoDream/Core/femtoDreamSelection.h"
#include "PWGCF/FemtoDream/Core/femtoDreamTrackSelection.h"

#include "Common/Core/RecoDecay.h"

#include "Framework/HistogramRegistry.h"
#include "ReconstructionDataFormats/PID.h"

#include "Math/Vector4D.h"
#include "TMath.h"

#include <array>
#include <cstdint>
#include <fstream>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;

namespace o2::analysis::femtoDream
{
namespace femtoDreamResoSelection
{
enum ResoSel {
  kResoSign
};
/// If you add a new selection, adjust kNresoSelection

enum Daughtertype {
  kPosdaugh,
  kNegdaugh
};

enum ResoMothers {
  kPhi,
  kKStar
};
} // namespace femtoDreamResoSelection

class FemtoDreamResoSelection
  : public FemtoDreamObjectSelection<float, femtoDreamResoSelection::ResoSel>
{

 public:
  FemtoDreamResoSelection() /// initialization currently kind of random change this!!!
    : mDaughPTPCThr(99.f), mPIDoffsetTPC(0.f), mPIDoffsetTOF(0.f)
  {
  }

  virtual ~FemtoDreamResoSelection() = default;

  template <typename V>
  uint32_t getType(V const& track1, V const& track2);

  template <typename V>
  size_t numBitsUsed(V const& origvalue);

  template <aod::femtodreamparticle::ParticleType part,
            aod::femtodreamparticle::ParticleType PartDaugh>
  void init(HistogramRegistry* QAregistry, HistogramRegistry* Registry);

  template <aod::femtodreamparticle::ParticleType part,
            aod::femtodreamparticle::TrackType trackType1,
            aod::femtodreamparticle::TrackType trackType2, typename T>
  void fillQA(T const& track1, T const& track2);

  template <typename T, typename V>
  void setDaughterCuts(femtoDreamResoSelection::Daughtertype child, T selVal,
                       V selVar, femtoDreamSelection::SelectionType selType);

  template <typename T, typename V>
  void setDaughterPIDSpecies(T const& daugh, V& pids);

  template <typename V>
  bool daughterSelectionPos(V const& track1, bool useThreshold);

  template <typename V>
  bool daughterSelectionNeg(V const& track2, bool useThreshold);

  template <typename cutContainerType, typename V>
  std::array<cutContainerType, 5> getCutContainer(V const& track1, V const& track2, float sign);

  template <typename T>
  std::pair<bool, bool> checkCombination(T const& PosTrack, T const& NegTrack, femtoDreamResoSelection::ResoMothers mother);

  std::pair<o2::track::PID::ID, o2::track::PID::ID> getPIDPairFromMother(femtoDreamResoSelection::ResoMothers mother);

  template <typename T>
  bool checkPID(T const& Track, float nSigTPC, float nSigTOF, float nSig2TPC, float nSig2TOF);

  void updateThreshold()
  {
    mDaughPTPCThr = assignedValue;
  };

  void setDaughternSigmaPIDOffset(femtoDreamResoSelection::Daughtertype daugh, float offsetTPC, float offsetTOF)
  {
    if (daugh == femtoDreamResoSelection::kPosdaugh) {
      posDaughTrack.setnSigmaPIDOffset(offsetTPC, offsetTOF);
    } else if (daugh == femtoDreamResoSelection::kNegdaugh) {
      negDaughTrack.setnSigmaPIDOffset(offsetTPC, offsetTOF);
    }
    mPIDoffsetTPC = offsetTPC;
    mPIDoffsetTOF = offsetTOF;
  };

  float getMass(o2::track::PID::ID pid)
  {
    switch (pid) {
      case (o2::track::PID::Kaon):
        return o2::constants::physics::MassKPlus;
      case (o2::track::PID::Pion):
        return o2::constants::physics::MassPiPlus;
      default:
        LOG(warn) << "PID not implemented in femtoDreamResoSelection.getMass";
        return 0.;
    }
  }

  std::pair<float, float> getMassDaughters(femtoDreamResoSelection::ResoMothers mother);

  /// The following functions might not be needed, as right now there is only one ResoSel (sign).
  /// However all the other selections are implemented this way (also in the CutCulator).
  /// So for now this is implemented analogous (migth also be beneficial if further ResoSels want to be implemented).

  /// Helper function to obtain the name of a given selection criterion for consistent naming of the configurables
  /// \param iSel Reso selection variable to be examined
  /// \param prefix Additional prefix for the name of the configurable
  /// \param suffix Additional suffix for the name of the configurable
  static std::string getSelectionName(femtoDreamResoSelection::ResoSel iSel,
                                      std::string_view prefix = "",
                                      std::string_view suffix = "")
  {
    std::string outString = static_cast<std::string>(prefix);
    outString += static_cast<std::string>(kSelectionNames[iSel]);
    outString += suffix;
    return outString;
  }

  /// Helper function to obtain the index of a given selection variable for consistent naming of the configurables
  /// \param obs Reso selection variable (together with prefix) got from file
  /// \param prefix Additional prefix for the output of the configurable
  static int findSelectionIndex(const std::string_view& obs,
                                std::string_view prefix = "")
  {
    for (int index = 0; index < kNresoSelection; index++) {
      std::string comp = static_cast<std::string>(prefix) +
                         static_cast<std::string>(kSelectionNames[index]);
      std::string_view cmp{comp};
      if (obs.compare(cmp) == 0)
        return index;
    }
    LOGF(info, "Variable %s not found", obs);
    return -1;
  }

  /// Helper function to obtain the type of a given selection variable for consistent naming of the configurables
  /// \param iSel Reso selection variable whose type is returned
  static femtoDreamSelection::SelectionType
    getSelectionType(femtoDreamResoSelection::ResoSel iSel)
  {
    return kSelectionTypes[iSel];
  }

  /// for consistent description of the configurables
  /// \param iSel Track selection variable to be examined
  /// \param prefix Additional prefix for the output of the configurable
  static std::string getSelectionHelper(femtoDreamResoSelection::ResoSel iSel,
                                        std::string_view prefix = "")
  {
    std::string outString = static_cast<std::string>(prefix);
    outString += static_cast<std::string>(kSelectionHelper[iSel]);
    return outString;
  }

 private:
  float mDaughPTPCThr;
  float mPIDoffsetTPC;
  float mPIDoffsetTOF;

  FemtoDreamTrackSelection posDaughTrack;
  FemtoDreamTrackSelection negDaughTrack;

  static constexpr int kNresoSelection = 1;

  static constexpr std::string_view kSelectionNames[kNresoSelection] = {"Sign"};

  static constexpr femtoDreamSelection::SelectionType kSelectionTypes[kNresoSelection]{
    femtoDreamSelection::kEqual};

  static constexpr std::string_view kSelectionHelper[kNresoSelection] = {
    "+1 for Reso, -1 for AntiReso"};

}; // namespace femtoDream

template <typename V>
uint32_t FemtoDreamResoSelection::getType(V const& track1, V const& track2)
{
  if (track1.pt() <= mDaughPTPCThr && track2.pt() <= mDaughPTPCThr) {
    return aod::femtodreamparticle::kPhiPosdaughTPC_NegdaughTPC;
  }
  if (track1.pt() <= mDaughPTPCThr && track2.pt() > mDaughPTPCThr) {
    return aod::femtodreamparticle::kPhiPosdaughTPC_NegdaughTOF;
  }
  if (track1.pt() > mDaughPTPCThr && track2.pt() <= mDaughPTPCThr) {
    return aod::femtodreamparticle::kPhiPosdaughTOF_NegdaughTPC;
  }
  if (track1.pt() > mDaughPTPCThr && track2.pt() > mDaughPTPCThr) {
    return aod::femtodreamparticle::kPhiPosdaughTOF_NegdaughTOF;
  }
  return 255; // as error filler
}

template <typename V>
size_t FemtoDreamResoSelection::numBitsUsed(V const& origvalue)
{
  size_t bits = 0;
  auto value = origvalue;
  while (value != 0) {
    ++bits;
    value >>= 1;
  }
  return bits;
}

template <aod::femtodreamparticle::ParticleType part,
          aod::femtodreamparticle::ParticleType PartDaugh>
void FemtoDreamResoSelection::init(HistogramRegistry* QAregistry, HistogramRegistry* Registry)
{
  if (QAregistry && Registry) {
    this->mHistogramRegistry = Registry;
    this->mQAHistogramRegistry = QAregistry;

    posDaughTrack.init<PartDaugh,
                       aod::femtodreamparticle::kPosChild,
                       aod::femtodreamparticle::cutContainerType>(
      mQAHistogramRegistry, mHistogramRegistry);

    negDaughTrack.init<PartDaugh,
                       aod::femtodreamparticle::kNegChild,
                       aod::femtodreamparticle::cutContainerType>(
      mQAHistogramRegistry, mHistogramRegistry);
  }
}

template <aod::femtodreamparticle::ParticleType part,
          aod::femtodreamparticle::TrackType trackType1,
          aod::femtodreamparticle::TrackType trackType2, typename T>
void FemtoDreamResoSelection::fillQA(T const& track1, T const& track2)
{
  posDaughTrack.fillQA<part, trackType1>(track1);
  negDaughTrack.fillQA<part, trackType2>(track2);
}

template <typename T, typename V>
void FemtoDreamResoSelection::setDaughterCuts(femtoDreamResoSelection::Daughtertype daugh, T selVal,
                                              V selVar, femtoDreamSelection::SelectionType selType)
{
  if (daugh == femtoDreamResoSelection::kPosdaugh) {
    posDaughTrack.setSelection(selVal, selVar, selType);
  };
  if (daugh == femtoDreamResoSelection::kNegdaugh) {
    negDaughTrack.setSelection(selVal, selVar, selType);
  }
}

template <typename T, typename V>
void FemtoDreamResoSelection::setDaughterPIDSpecies(T const& daugh, V& pids)
{
  if (daugh == femtoDreamResoSelection::kPosdaugh) {
    posDaughTrack.setPIDSpecies(pids);
  }
  if (daugh == femtoDreamResoSelection::kNegdaugh) {
    negDaughTrack.setPIDSpecies(pids);
  }
}

template <typename V>
bool FemtoDreamResoSelection::daughterSelectionPos(V const& track1, bool useThreshold)
{
  return posDaughTrack.isSelectedMinimal(track1, useThreshold);
}

template <typename V>
bool FemtoDreamResoSelection::daughterSelectionNeg(V const& track2, bool useThreshold)
{
  return negDaughTrack.isSelectedMinimal(track2, useThreshold);
}

template <typename T>
std::pair<bool, bool> FemtoDreamResoSelection::checkCombination(T const& PosTrack, T const& NegTrack, femtoDreamResoSelection::ResoMothers mother)
{
  /// first bool: normal or anti
  /// second bool: is not a valid combination

  auto [part1, part2] = getPIDPairFromMother(mother);

  float nSigPosTPC1 = o2::aod::pidutils::tpcNSigma(part1, PosTrack) - mPIDoffsetTPC;
  float nSigPosTOF1 = posDaughTrack.getNsigmaTOF(PosTrack, part1) - mPIDoffsetTOF; /// for TOF use function in TrackSelection, because it also checks hasTOF()
  float nSigPosTPC2 = o2::aod::pidutils::tpcNSigma(part2, PosTrack) - mPIDoffsetTPC;
  float nSigPosTOF2 = posDaughTrack.getNsigmaTOF(PosTrack, part2) - mPIDoffsetTOF;
  float nSigNegTPC1 = o2::aod::pidutils::tpcNSigma(part1, NegTrack) - mPIDoffsetTPC;
  float nSigNegTOF1 = negDaughTrack.getNsigmaTOF(NegTrack, part1) - mPIDoffsetTOF;
  float nSigNegTPC2 = o2::aod::pidutils::tpcNSigma(part2, NegTrack) - mPIDoffsetTPC;
  float nSigNegTOF2 = negDaughTrack.getNsigmaTOF(NegTrack, part2) - mPIDoffsetTOF;

  if (checkPID(PosTrack, nSigPosTPC1, nSigPosTOF1, nSigPosTPC2, nSigPosTOF2) && checkPID(NegTrack, nSigNegTPC2, nSigNegTOF2, nSigNegTPC1, nSigNegTOF1)) {
    return {true, false};
  } else if (checkPID(PosTrack, nSigPosTPC2, nSigPosTOF2, nSigPosTPC1, nSigPosTOF1) && checkPID(NegTrack, nSigNegTPC1, nSigNegTOF1, nSigNegTPC2, nSigNegTOF2)) {
    return {false, false};
  } else {
    return {false, true};
  }
}

template <typename T>
bool FemtoDreamResoSelection::checkPID(T const& Track, float nSig1TPC, float nSig1TOF, float nSig2TPC, float nSig2TOF)
{
  if (Track.pt() < mDaughPTPCThr) {
    return (std::abs(nSig1TPC) <= std::abs(nSig2TPC));
  } else {
    return (std::sqrt(nSig1TPC * nSig1TPC + nSig1TOF * nSig1TOF) <= std::sqrt(nSig2TPC * nSig2TPC + nSig2TOF * nSig2TOF));
  }
}

std::pair<o2::track::PID::ID, o2::track::PID::ID> FemtoDreamResoSelection::getPIDPairFromMother(femtoDreamResoSelection::ResoMothers mother)
{
  /// return is structured this way:
  /// The mother particle is assumed to be normal (not antiparticle). Then:
  /// 1. return value is positive daughter
  /// 2. return value is negative daughter
  switch (mother) {
    case (femtoDreamResoSelection::kPhi):
      return {o2::track::PID::Kaon, o2::track::PID::Kaon};
    case (femtoDreamResoSelection::kKStar):
      return {o2::track::PID::Kaon, o2::track::PID::Pion};
    default:
      LOG(warn) << "MotherPID not implemented in femtoDreamResoSelection.getPIDPairFromMother";
      return {o2::track::PID::Kaon, o2::track::PID::Kaon};
  }
}

std::pair<float, float> FemtoDreamResoSelection::getMassDaughters(femtoDreamResoSelection::ResoMothers mother)
{
  switch (mother) {
    case (femtoDreamResoSelection::kPhi):
      return {o2::constants::physics::MassKPlus, o2::constants::physics::MassKPlus};
    case (femtoDreamResoSelection::kKStar):
      return {o2::constants::physics::MassKPlus, o2::constants::physics::MassPiPlus};
    default:
      LOG(warn) << "MotherPID not implemented in femtoDreamResoSelection.getMassDauhters";
      return {o2::constants::physics::MassKPlus, o2::constants::physics::MassKPlus};
  }
}

//// new getCutContainer
template <typename cutContainerType, typename V>
std::array<cutContainerType, 5> FemtoDreamResoSelection::getCutContainer(V const& track1, V const& track2, float sign)
{
  cutContainerType outputSign = 0;
  cutContainerType outputPID = 0;
  size_t counter = 0;
  for (auto& sel : mSelections) { /// it should just be a 1D vector with sign
    const auto selVariable = sel.getSelectionVariable();
    if (selVariable == femtoDreamResoSelection::kResoSign) {
      sel.checkSelectionSetBit(sign, outputSign, counter, nullptr);
    }
  }

  const auto dCA1 = std::sqrt(track1.dcaXY() * track1.dcaXY() + track1.dcaZ() * track1.dcaZ());
  const auto dCA2 = std::sqrt(track2.dcaXY() * track2.dcaXY() + track2.dcaZ() * track2.dcaZ());

  auto outputPosTrack = posDaughTrack.getCutContainer<false, cutContainerType>(track1, track1.pt(), track1.eta(), dCA1); // false for useItsPid
  auto outputNegTrack = negDaughTrack.getCutContainer<false, cutContainerType>(track2, track2.pt(), track2.eta(), dCA2);

  const auto shiftvalue = numBitsUsed(outputSign);
  outputPID = (outputNegTrack.at(femtoDreamTrackSelection::TrackContainerPosition::kPID) << shiftvalue) | outputSign;

  std::array<cutContainerType, 5> bitmask = {outputPID,
                                             outputPosTrack.at(femtoDreamTrackSelection::TrackContainerPosition::kCuts),
                                             outputPosTrack.at(femtoDreamTrackSelection::TrackContainerPosition::kPID),
                                             outputNegTrack.at(femtoDreamTrackSelection::TrackContainerPosition::kCuts),
                                             outputNegTrack.at(femtoDreamTrackSelection::TrackContainerPosition::kPID)};
  return bitmask;
}
} // namespace o2::analysis::femtoDream

#endif // PWGCF_FEMTODREAM_CORE_FEMTODREAMRESOSELECTION_H_
