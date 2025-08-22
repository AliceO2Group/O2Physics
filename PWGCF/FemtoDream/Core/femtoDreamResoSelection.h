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

/// \file FemtoDreamResoSelection.h
/// \brief Definition of the FemtoDreamResoSelection
/// \author Valentina Mantovani Sarti, TU München valentina.mantovani-sarti@tum.de
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de
/// \author Luca Barioglio, TU München, luca.barioglio@cern.ch

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
#include <iostream>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;

/// std::ofstream wtftest("logResoHeader.txt", std::ios::app);

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

} // namespace femtoDreamResoSelection

class FemtoDreamResoSelection
  : public FemtoDreamObjectSelection<float, femtoDreamResoSelection::ResoSel>
{

 public:
  FemtoDreamResoSelection() /// initialization currently kind of random change this!!!
    : mDaughPTPCThr(99.f)
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
  bool DaughterSelection1(V const& track1, bool UseThreshold);

  template <typename V>
  bool DaughterSelection2(V const& track2, bool UseThreshold);

  template <typename cutContainerType, typename V>
  std::array<cutContainerType, 5> getCutContainer(V const& track1, V const& track2, float sign);

  void updateMembersMinimal()
  {
    mDaughPTPCThr = assignedValue;
  };

  void setDaughternSigmaPIDOffset(femtoDreamResoSelection::Daughtertype daugh, float offsetTPC, float offsetTOF)
  {
    if (daugh == femtoDreamResoSelection::kPosdaugh) {
      PosDaughTrack.setnSigmaPIDOffset(offsetTPC, offsetTOF);
    } else if (daugh == femtoDreamResoSelection::kNegdaugh) {
      NegDaughTrack.setnSigmaPIDOffset(offsetTPC, offsetTOF);
    }
  };

  float getMass(o2::track::PID::ID pid)
  {
    switch (pid) {
      case (o2::track::PID::Kaon):
        return o2::constants::physics::MassKPlus;
      case (o2::track::PID::Pion):
        return o2::constants::physics::MassPiPlus;
      default:
        LOG(fatal) << "PID not implemented in femtoDreamResoSelection.getMass";
        return 0.;
    }
  }

  /// @brief checks if the reso-particle is either particle or antiparticle (dependent on isParticle)
  /// @param PID vector of PIDs of daughter particles
  /// @param isParticle true: function checks if Reso is particle, false: function checks if Reso is Antiparticle
  template <typename T>
  bool isResoOfKind(T const& PosTrack, T const& NegTrack, std::vector<int> const& PID, bool isParticle);

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
    outString += static_cast<std::string>(mSelectionNames[iSel]);
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
                         static_cast<std::string>(mSelectionNames[index]);
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
    return mSelectionTypes[iSel];
  }

  /// for consistent description of the configurables
  /// \param iSel Track selection variable to be examined
  /// \param prefix Additional prefix for the output of the configurable
  static std::string getSelectionHelper(femtoDreamResoSelection::ResoSel iSel,
                                        std::string_view prefix = "")
  {
    std::string outString = static_cast<std::string>(prefix);
    outString += static_cast<std::string>(mSelectionHelper[iSel]);
    return outString;
  }

 private:
  float mDaughPTPCThr;

  FemtoDreamTrackSelection PosDaughTrack;
  FemtoDreamTrackSelection NegDaughTrack;

  static constexpr int kNresoSelection = 1;

  static constexpr std::string_view mSelectionNames[kNresoSelection] = {"Sign"};

  static constexpr femtoDreamSelection::SelectionType mSelectionTypes[kNresoSelection]{
    femtoDreamSelection::kEqual};

  static constexpr std::string_view mSelectionHelper[kNresoSelection] = {
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

    PosDaughTrack.init<PartDaugh,
                       aod::femtodreamparticle::kPosChild,
                       aod::femtodreamparticle::cutContainerType>(
      mQAHistogramRegistry, mHistogramRegistry);

    NegDaughTrack.init<PartDaugh,
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
  PosDaughTrack.fillQA<part, trackType1>(track1);
  NegDaughTrack.fillQA<part, trackType2>(track2);
}

template <typename T, typename V>
void FemtoDreamResoSelection::setDaughterCuts(femtoDreamResoSelection::Daughtertype daugh, T selVal,
                                              V selVar, femtoDreamSelection::SelectionType selType)
{
  if (daugh == femtoDreamResoSelection::kPosdaugh) {
    PosDaughTrack.setSelection(selVal, selVar, selType);
  }
  if (daugh == femtoDreamResoSelection::kNegdaugh) {
    NegDaughTrack.setSelection(selVal, selVar, selType);
  }
}

template <typename T, typename V>
void FemtoDreamResoSelection::setDaughterPIDSpecies(T const& daugh, V& pids)
{
  if (daugh == femtoDreamResoSelection::kPosdaugh) {
    PosDaughTrack.setPIDSpecies(pids);
  }
  if (daugh == femtoDreamResoSelection::kNegdaugh) {
    NegDaughTrack.setPIDSpecies(pids);
  }
}

template <typename V>
bool FemtoDreamResoSelection::DaughterSelection1(V const& track1, bool UseThreshold)
{
  if (!PosDaughTrack.isSelectedMinimal(track1, UseThreshold)) {
    return false;
  }
  return true;
}

template <typename V>
bool FemtoDreamResoSelection::DaughterSelection2(V const& track2, bool UseThreshold)
{
  if (!NegDaughTrack.isSelectedMinimal(track2, UseThreshold)) {
    return false;
  }
  return true;
}

template <typename T>
bool FemtoDreamResoSelection::isResoOfKind(T const& PosTrack, T const& NegTrack, std::vector<int> const& PID, bool isParticle)
{
  int part1 = PID[0];
  int part2 = PID[1];
  if (!isParticle) {
    part1 = PID[1];
    part2 = PID[0];
  }

  bool resultPos = false;
  float nSigPosTPC1 = o2::aod::pidutils::tpcNSigma(part1, PosTrack);
  float nSigPosTOF1 = PosDaughTrack.getNsigmaTOF(PosTrack, part1); /// for TOF use function in TrackSelection, because it also checks hasTOF()
  float nSigPosTPC2 = o2::aod::pidutils::tpcNSigma(part2, PosTrack);
  float nSigPosTOF2 = PosDaughTrack.getNsigmaTOF(PosTrack, part2);
  if (PosTrack.pt() < mDaughPTPCThr) {
    resultPos = (std::abs(nSigPosTPC1) <= std::abs(nSigPosTPC2));
  } else {
    resultPos = (std::sqrt(nSigPosTPC1 * nSigPosTPC1 + nSigPosTOF1 * nSigPosTOF1) <= std::sqrt(nSigPosTPC2 * nSigPosTPC2 + nSigPosTOF2 * nSigPosTOF2));
  }

  bool resultNeg = false;
  float nSigNegTPC1 = o2::aod::pidutils::tpcNSigma(part1, NegTrack);
  float nSigNegTOF1 = NegDaughTrack.getNsigmaTOF(NegTrack, part1);
  float nSigNegTPC2 = o2::aod::pidutils::tpcNSigma(part2, NegTrack);
  float nSigNegTOF2 = NegDaughTrack.getNsigmaTOF(NegTrack, part2);
  if (NegTrack.pt() < mDaughPTPCThr) {
    resultNeg = (std::abs(nSigNegTPC2) <= std::abs(nSigNegTPC1));
  } else {
    resultNeg = (std::sqrt(nSigNegTPC2 * nSigNegTPC2 + nSigNegTOF2 * nSigNegTOF2) <= std::sqrt(nSigNegTPC1 * nSigNegTPC1 + nSigNegTOF1 * nSigNegTOF1));
  }
  return (resultPos && resultNeg);
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

  const auto DCA1 = std::sqrt(track1.dcaXY() * track1.dcaXY() + track1.dcaZ() * track1.dcaZ());
  const auto DCA2 = std::sqrt(track2.dcaXY() * track2.dcaXY() + track2.dcaZ() * track2.dcaZ());

  auto outputPosTrack = PosDaughTrack.getCutContainer<false, cutContainerType>(track1, track1.pt(), track1.eta(), DCA1); // false for useItsPid
  auto outputNegTrack = NegDaughTrack.getCutContainer<false, cutContainerType>(track2, track2.pt(), track2.eta(), DCA2);

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
