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
#include <string>
#include <utility>
#include <vector>

namespace o2::analysis::femtoDream // o2-linter: disable=name/namespace (Previously defined namespace)
{
namespace femto_dream_reso_selection
{

enum ResoSel {
  kResoSign
};
/// If you add a new selection, adjust kNresoSelection

enum Daughtertype {
  kPosdaugh,
  kNegdaugh
};
} // namespace femto_dream_reso_selection

class FemtoDreamResoSelection
  : public FemtoDreamObjectSelection<float, femto_dream_reso_selection::ResoSel>
{

 public:
  FemtoDreamResoSelection() /// initialization currently kind of random change this!!!
    : mDaughPTPCThr{99.f, 99.f}, mPIDoffsetTPC(0.f), mPIDoffsetTOF(0.f), mSigmaPIDMax(99.f)
  {
  }

  virtual ~FemtoDreamResoSelection() = default;

  template <typename V>
  int getType(V const& track1, V const& track2, bool resoIsNotAnti);

  /// assigns value from configurbale to private class member
  template <typename V>
  void assign(V& selVals)
  {
    mDaughPTPCThr = selVals;
  };

  template <typename V>
  size_t numBitsUsed(V const& origvalue);

  template <aod::femtodreamparticle::ParticleType part, aod::femtodreamparticle::ParticleType PartDaugh>
  void init(HistogramRegistry* QAregistry, HistogramRegistry* Registry);

  template <aod::femtodreamparticle::ParticleType part,
            aod::femtodreamparticle::TrackType trackType1,
            aod::femtodreamparticle::TrackType trackType2, typename T>
  void fillQA(T const& track1, T const& track2);

  template <typename T, typename V>
  void setDaughterCuts(femto_dream_reso_selection::Daughtertype child, T selVal,
                       V selVar, femtoDreamSelection::SelectionType selType);

  template <typename T, typename V>
  void setDaughterPIDSpecies(T const& daugh, V& pids);

  template <typename V>
  bool daughterSelectionPos(V const& track1);

  template <typename V>
  bool daughterSelectionNeg(V const& track2);

  template <typename V, typename T>
  bool isSelectedMinimalPIDPos(V const& track1, T const& pidVector);

  template <typename V, typename T>
  bool isSelectedMinimalPIDNeg(V const& track2, T const& pidVector);

  template <typename cutContainerType, typename V>
  std::array<cutContainerType, 5> getCutContainer(V const& track1, V const& track2, float sign);

  template <typename T, typename V>
  std::pair<bool, bool> checkCombination(T const& PosTrack, T const& NegTrack, V const& pidVector, float massDiff, float massDiffAnti, bool useMassDiff);

  template <typename T, typename V>
  float getNSigTotal(T const& track, V const& pid, float const& threshold);

  void updateThreshold()
  {
    mSigmaPIDMax = posDaughTrack.getMinimalSelection(o2::analysis::femtoDream::femtoDreamTrackSelection::kPIDnSigmaMax, femtoDreamSelection::kAbsUpperLimit); // the same for pos and neg
  };

  void setDaughternSigmaPIDOffset(femto_dream_reso_selection::Daughtertype daugh, float offsetTPC, float offsetTOF)
  {
    if (daugh == femto_dream_reso_selection::kPosdaugh) {
      posDaughTrack.setnSigmaPIDOffset(offsetTPC, offsetTOF);
    } else if (daugh == femto_dream_reso_selection::kNegdaugh) {
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
        LOG(warn) << "PID not implemented in femto_dream_reso_selection.getMass";
        return 0.;
    }
  }

  /// The following functions might not be needed, as right now there is only one ResoSel (sign).
  /// However all the other selections are implemented this way (also in the CutCulator).
  /// So for now this is implemented analogous (migth also be beneficial if further ResoSels want to be implemented).

  /// Helper function to obtain the name of a given selection criterion for consistent naming of the configurables
  /// \param iSel Reso selection variable to be examined
  /// \param prefix Additional prefix for the name of the configurable
  /// \param suffix Additional suffix for the name of the configurable
  static std::string getSelectionName(femto_dream_reso_selection::ResoSel iSel,
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
  static int findSelectionIndex(std::string_view obs,
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
  static femtoDreamSelection::SelectionType // o2-linter: disable=name/function-variable (defined with UpperCamelCase in femtoDreamSelection)
    getSelectionType(femto_dream_reso_selection::ResoSel iSel)
  {
    return kSelectionTypes[iSel];
  }

  /// for consistent description of the configurables
  /// \param iSel Track selection variable to be examined
  /// \param prefix Additional prefix for the output of the configurable
  static std::string getSelectionHelper(femto_dream_reso_selection::ResoSel iSel,
                                        std::string_view prefix = "")
  {
    std::string outString = static_cast<std::string>(prefix);
    outString += static_cast<std::string>(kSelectionHelper[iSel]);
    return outString;
  }

 private:
  std::vector<float> mDaughPTPCThr;
  float mPIDoffsetTPC;
  float mPIDoffsetTOF;
  float mSigmaPIDMax;

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
int FemtoDreamResoSelection::getType(V const& track1, V const& track2, bool resoIsNotAnti)
{
  float posThresh = 0.;
  float negThresh = 0.;
  if (resoIsNotAnti) {
    posThresh = mDaughPTPCThr[0];
    negThresh = mDaughPTPCThr[1];
  } else {
    posThresh = mDaughPTPCThr[1];
    negThresh = mDaughPTPCThr[0];
  }

  if (track1.pt() <= posThresh && track2.pt() <= negThresh) {
    return aod::femtodreamparticle::kResoPosdaughTPC_NegdaughTPC;
  }
  if (track1.pt() <= posThresh && track2.pt() > negThresh) {
    return aod::femtodreamparticle::kResoPosdaughTPC_NegdaughTOF;
  }
  if (track1.pt() > posThresh && track2.pt() <= negThresh) {
    return aod::femtodreamparticle::kResoPosdaughTOF_NegdaughTPC;
  }
  if (track1.pt() > posThresh && track2.pt() > negThresh) {
    return aod::femtodreamparticle::kResoPosdaughTOF_NegdaughTOF;
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
void FemtoDreamResoSelection::setDaughterCuts(femto_dream_reso_selection::Daughtertype daugh, T selVal,
                                              V selVar, femtoDreamSelection::SelectionType selType)
{
  if (daugh == femto_dream_reso_selection::kPosdaugh) {
    posDaughTrack.setSelection(selVal, selVar, selType);
  }
  if (daugh == femto_dream_reso_selection::kNegdaugh) {
    negDaughTrack.setSelection(selVal, selVar, selType);
  }
}

template <typename T, typename V>
void FemtoDreamResoSelection::setDaughterPIDSpecies(T const& daugh, V& pids)
{
  if (daugh == femto_dream_reso_selection::kPosdaugh) {
    posDaughTrack.setPIDSpecies(pids);
  }
  if (daugh == femto_dream_reso_selection::kNegdaugh) {
    negDaughTrack.setPIDSpecies(pids);
  }
}

template <typename V>
bool FemtoDreamResoSelection::daughterSelectionPos(V const& track1)
{
  return posDaughTrack.isSelectedMinimal(track1);
}

template <typename V>
bool FemtoDreamResoSelection::daughterSelectionNeg(V const& track2)
{
  return negDaughTrack.isSelectedMinimal(track2);
}

template <typename V, typename T>
bool FemtoDreamResoSelection::isSelectedMinimalPIDPos(V const& track1, T const& pidVector)
{
  int pidVecSize = pidVector.size();
  for (int i = 0; i < pidVecSize; i++) {
    const float pidTPC = posDaughTrack.getNsigmaTPC(track1, pidVector[i]);
    const float pidTOF = posDaughTrack.getNsigmaTOF(track1, pidVector[i]);

    if (track1.pt() < mDaughPTPCThr[i]) {
      if (std::fabs(pidTPC) < mSigmaPIDMax) {
        return true;
      }
    } else if ((std::sqrt(pidTPC * pidTPC + pidTOF * pidTOF) < mSigmaPIDMax)) {
      return true;
    }
  }
  return false;
}

template <typename V, typename T>
bool FemtoDreamResoSelection::isSelectedMinimalPIDNeg(V const& track2, T const& pidVector)
{
  int pidVecSize = pidVector.size();
  for (int i = 0; i < pidVecSize; i++) {
    const float pidTPC = negDaughTrack.getNsigmaTPC(track2, pidVector[i]);
    const float pidTOF = negDaughTrack.getNsigmaTOF(track2, pidVector[i]);

    if (track2.pt() < mDaughPTPCThr[i]) {
      if (std::fabs(pidTPC) < mSigmaPIDMax) {
        return true;
      }
    } else if ((std::sqrt(pidTPC * pidTPC + pidTOF * pidTOF) < mSigmaPIDMax)) {
      return true;
    }
  }
  return false;
}

template <typename T, typename V>
std::pair<bool, bool> FemtoDreamResoSelection::checkCombination(T const& PosTrack, T const& NegTrack, V const& pidVector, float massDiff, float massDiffAnti, bool useMassDiff)
{
  /// first bool: true (normal resonance) / false (anti resonance)
  /// second bool: is not a valid combination

  const auto part1 = pidVector[0]; /// particle type 1
  const auto part2 = pidVector[1]; /// particle type 2

  float nSigPosPart1Total = getNSigTotal(PosTrack, part1, mDaughPTPCThr[0]); /// Total propability that PosTrack is of particle type 1
  float nSigPosPart2Total = getNSigTotal(PosTrack, part2, mDaughPTPCThr[1]); /// Total propability that PosTrack is of particle type 2
  float nSigNegPart1Total = getNSigTotal(NegTrack, part1, mDaughPTPCThr[0]);
  float nSigNegPart2Total = getNSigTotal(NegTrack, part2, mDaughPTPCThr[1]);

  // check if PosTrack is more likely to be part1 than part2 (and vice versa for NegTrack) -> normal resonance
  bool couldBeNormal = nSigPosPart1Total < nSigPosPart2Total && nSigNegPart2Total < nSigNegPart1Total;
  // check if PosTrack is more likely to be part2 than part1 (and vice versa for NegTrack) -> anti resonance
  bool couldBeAnti = nSigPosPart2Total < nSigPosPart1Total && nSigNegPart1Total < nSigNegPart2Total;

  if (useMassDiff) {
    couldBeNormal = couldBeNormal && massDiff < massDiffAnti;
    couldBeAnti = couldBeAnti && massDiffAnti < massDiff;
  }

  if (couldBeNormal && !couldBeAnti) {
    return {true, false};
  }
  if (!couldBeNormal && couldBeAnti) {
    return {false, false};
  }
  // if ambiguous (both true) or invalid (both false)
  return {false, true};
}

template <typename T, typename V>
float FemtoDreamResoSelection::getNSigTotal(T const& track, V const& pid, float const& threshold)
{
  float nSigTPC = o2::aod::pidutils::tpcNSigma(pid, track);
  float pTtrack = track.pt();

  if (pTtrack < threshold) {
    return std::abs(nSigTPC);
  }

  float nSigTOF = track.hasTOF() ? o2::aod::pidutils::tofNSigma(pid, track) : 999.f;

  return std::hypot(nSigTPC, nSigTOF);
}

//// new getCutContainer
template <typename cutContainerType, typename V>
std::array<cutContainerType, 5> FemtoDreamResoSelection::getCutContainer(V const& track1, V const& track2, float sign)
{
  cutContainerType outputSign = 0;
  cutContainerType outputPID = 0;
  size_t counter = 0;
  for (auto& sel : mSelections) { // o2-linter: disable=const-ref-in-for-loop (femtoDreamObjectSelection has no const getter)
    const auto selVariable = sel.getSelectionVariable();
    if (selVariable == femto_dream_reso_selection::kResoSign) {
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
