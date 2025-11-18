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

/// \file femtoDreamV0SelectionK0Short.h
/// \brief Definition of the FemtoDreamV0Selection
/// \author Valentina Mantovani Sarti, TU München valentina.mantovani-sarti@tum.de
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de
/// \author Luca Barioglio, TU München, luca.barioglio@cern.ch
/// \author Christopher Klumm, TU München, christopher.klumm@cern.ch
/// \author Nils Fabian Konert, TU München, nils.fabian.konert@cern.ch

#ifndef PWGCF_FEMTODREAM_CORE_FEMTODREAMV0SELECTIONK0SHORT_H_
#define PWGCF_FEMTODREAM_CORE_FEMTODREAMV0SELECTIONK0SHORT_H_

#include "PWGCF/FemtoDream/Core/femtoDreamObjectSelection.h"
#include "PWGCF/FemtoDream/Core/femtoDreamSelection.h"
#include "PWGCF/FemtoDream/Core/femtoDreamTrackSelection.h"

#include "Common/Core/RecoDecay.h"

#include "Framework/HistogramRegistry.h"
#include "ReconstructionDataFormats/PID.h"

#include <string>
#include <vector>

namespace o2::analysis::femtoDream
{

namespace femto_dream_v0_selection
{
/// The different selections this task is capable of doing
enum V0Sel {
  kV0Sign, ///< +1 particle, -1 antiparticle
  kV0pTMin,
  kV0pTMax,
  kV0etaMax,
  kV0DCADaughMax,
  kV0CPAMin,
  kV0TranRadMin,
  kV0TranRadMax,
  kV0DecVtxMax,
};

enum ChildTrackType { kPosTrack,
                      kNegTrack };

enum V0ContainerPosition {
  kV0,
  kPosCuts,
  kPosPID,
  kNegCuts,
  kNegPID,
}; /// Position in the full VO cut container

} // namespace  femto_dream_v0_selection

/// \class FemtoDreamV0Selection
/// \brief Cut class to contain and execute all cuts applied to V0s
class FemtoDreamV0Selection
  : public FemtoDreamObjectSelection<float, femto_dream_v0_selection::V0Sel>
{
 public:
  FemtoDreamV0Selection()
    : nPtV0MinSel(0), nPtV0MaxSel(0), nEtaV0MaxSel(0), nDCAV0DaughMax(0), nCPAV0Min(0), nTranRadV0Min(0), nTranRadV0Max(0), nDecVtxMax(0), pTV0Min(9999999.), pTV0Max(-9999999.), etaV0Max(-9999999.), dCAV0DaughMax(-9999999.), cPAV0Min(9999999.), tranRadV0Min(9999999.), tranRadV0Max(-9999999.), decVtxMax(-9999999.), fInvMassLowLimit(1.05), fInvMassUpLimit(1.3), fRejectKaon(false), fRejectLambda(false), fInvMassKaonLowLimit(0.48), fInvMassKaonUpLimit(0.515), nSigmaPIDOffsetTPC(0.), fMotherIsLambda(true) {}

  /// Initializes histograms for the task
  template <o2::aod::femtodreamparticle::ParticleType part, o2::aod::femtodreamparticle::ParticleType daugh, typename cutContainerType>
  void init(HistogramRegistry* QAregistry, HistogramRegistry* Registry);

  template <typename C, typename V, typename T>
  bool isSelectedMinimal(C const& col, V const& v0, T const& posTrack,
                         T const& negTrack);

  template <o2::aod::femtodreamparticle::ParticleType part,
            typename C, typename V, typename T>
  void fillLambdaQA(C const& col, V const& v0, T const& posTrack,
                    T const& negTrack);

  /// \todo for the moment the PID of the tracks is factored out into a separate
  /// field, hence 5 values in total \\ASK: what does it mean?
  template <typename cutContainerType, typename C, typename V, typename T>
  std::array<cutContainerType, 5> getCutContainer(C const& col, V const& v0,
                                                  T const& posTrack,
                                                  T const& negTrack);

  template <o2::aod::femtodreamparticle::ParticleType part,
            o2::aod::femtodreamparticle::ParticleType daugh, typename C,
            typename V, typename T>
  void fillQA(C const& col, V const& v0, T const& posTrack, T const& negTrack);

  template <typename T1, typename T2>
  void setChildCuts(femto_dream_v0_selection::ChildTrackType child, T1 selVal, T2 selVar, femtoDreamSelection::SelectionType selType)
  {
    if (child == femto_dream_v0_selection::kPosTrack) {
      posDaughTrack.setSelection(selVal, selVar, selType);
    } else if (child == femto_dream_v0_selection::kNegTrack) {
      negDaughTrack.setSelection(selVal, selVar, selType);
    }
  }

  template <typename T>
  void setChildPIDSpecies(femto_dream_v0_selection::ChildTrackType child, T& pids)
  {
    if (child == femto_dream_v0_selection::kPosTrack) {
      posDaughTrack.setPIDSpecies(pids);
    } else if (child == femto_dream_v0_selection::kNegTrack) {
      negDaughTrack.setPIDSpecies(pids);
    }
  }

  /// Helper function to obtain the name of a given selection criterion for consistent naming of the configurables
  /// \param iSel Track selection variable to be examined
  /// \param prefix Additional prefix for the name of the configurable
  /// \param suffix Additional suffix for the name of the configurable
  static std::string getSelectionName(femto_dream_v0_selection::V0Sel iSel, std::string_view prefix = "", std::string_view suffix = "")
  {
    std::string outString = static_cast<std::string>(prefix);
    outString += static_cast<std::string>(SelectionNames[iSel]);
    outString += suffix;
    return outString;
  }

  /// Helper function to obtain the index of a given selection variable for consistent naming of the configurables
  /// \param obs V0 selection variable (together with prefix) got from file
  /// \param prefix Additional prefix for the output of the configurable
  static int findSelectionIndex(const std::string_view& obs,
                                std::string_view prefix = "")
  {
    for (int index = 0; index < kNv0Selection; index++) {
      std::string comp = static_cast<std::string>(prefix) +
                         static_cast<std::string>(SelectionNames[index]);
      std::string_view cmp{comp};
      if (obs.compare(cmp) == 0)
        return index;
    }
    LOGF(info, "Variable %s not found", obs);
    return -1;
  }

  /// Helper function to obtain the type of a given selection variable for consistent naming of the configurables
  /// \param iSel V0 selection variable whose type is returned
  static femtoDreamSelection::SelectionType getSelectionType(femto_dream_v0_selection::V0Sel iSel)
  { // o2-linter: disable=name/function-variable (defined with UpperCamelCase in femtoDreamSelection)
    return mSelectionTypes[iSel];
  }

  /// Helper function to obtain the helper string of a given selection criterion
  /// for consistent description of the configurables
  /// \param iSel Track selection variable to be examined
  /// \param prefix Additional prefix for the output of the configurable
  static std::string getSelectionHelper(femto_dream_v0_selection::V0Sel iSel,
                                        std::string_view prefix = "")
  {
    std::string outString = static_cast<std::string>(prefix);
    outString += static_cast<std::string>(SelectionHelper[iSel]);
    return outString;
  }

  /// Set limit for the selection on the invariant mass
  /// \param lowLimit Lower limit for the invariant mass distribution
  /// \param upLimit Upper limit for the invariant mass distribution
  void setInvMassLimits(float lowLimit, float upLimit)
  {
    fInvMassLowLimit = lowLimit;
    fInvMassUpLimit = upLimit;
  }

  /// Set limit for the kaon rejection on the invariant mass
  /// \param lowLimit Lower limit for the invariant mass distribution
  /// \param upLimit Upper limit for the invariant mass distribution
  void setKaonInvMassLimits(float lowLimit, float upLimit)
  {
    fRejectKaon = true;
    fInvMassKaonLowLimit = lowLimit;
    fInvMassKaonUpLimit = upLimit;
  }

  void setnSigmaPIDOffsetTPC(float offsetTPC)
  {
    nSigmaPIDOffsetTPC = offsetTPC;
  }

  void setChildRejectNotPropagatedTracks(femto_dream_v0_selection::ChildTrackType child, bool reject)
  {
    if (child == femto_dream_v0_selection::kPosTrack) {
      posDaughTrack.setRejectNotPropagatedTracks(reject);
    } else if (child == femto_dream_v0_selection::kNegTrack) {
      negDaughTrack.setRejectNotPropagatedTracks(reject);
    }
  }

  void setChildnSigmaPIDOffset(femto_dream_v0_selection::ChildTrackType child, float offsetTPC, float offsetTOF)
  {
    if (child == femto_dream_v0_selection::kPosTrack) {
      posDaughTrack.setnSigmaPIDOffset(offsetTPC, offsetTOF);
    } else if (child == femto_dream_v0_selection::kNegTrack) {
      negDaughTrack.setnSigmaPIDOffset(offsetTPC, offsetTOF);
    }
  }

  void setIsMother(bool IsMother)
  {
    fMotherIsLambda = IsMother;
  }

  void setRejectLambda(bool reject)
  {
    fRejectLambda = reject;
  }

 private:
  int nPtV0MinSel;
  int nPtV0MaxSel;
  int nEtaV0MaxSel;
  int nDCAV0DaughMax;
  int nCPAV0Min;
  int nTranRadV0Min;
  int nTranRadV0Max;
  int nDecVtxMax;
  float pTV0Min;
  float pTV0Max;
  float etaV0Max;
  float dCAV0DaughMax;
  float cPAV0Min;
  float tranRadV0Min;
  float tranRadV0Max;
  float decVtxMax;

  float fInvMassLowLimit;
  float fInvMassUpLimit;

  bool fRejectKaon;
  bool fRejectLambda;
  float fInvMassKaonLowLimit;
  float fInvMassKaonUpLimit;

  float nSigmaPIDOffsetTPC;

  bool fMotherIsLambda;

  FemtoDreamTrackSelection posDaughTrack;
  FemtoDreamTrackSelection negDaughTrack;

  static constexpr int kNv0Selection = 9;

  static constexpr std::string_view SelectionNames[kNv0Selection] = {
    "Sign", "PtMin", "PtMax", "EtaMax", "DCAdaughMax", "CPAMin",
    "TranRadMin", "TranRadMax", "DecVecMax"}; ///< Name of the different
                                              ///< selections

  static constexpr femtoDreamSelection::SelectionType
    mSelectionTypes[kNv0Selection]{
      femtoDreamSelection::kEqual,
      femtoDreamSelection::kLowerLimit,
      femtoDreamSelection::kUpperLimit,
      femtoDreamSelection::kUpperLimit,
      femtoDreamSelection::kUpperLimit,
      femtoDreamSelection::kLowerLimit,
      femtoDreamSelection::kLowerLimit,
      femtoDreamSelection::kUpperLimit,
      femtoDreamSelection::kUpperLimit}; ///< Map to match a variable with
                                         ///< its type

  static constexpr std::string_view SelectionHelper[kNv0Selection] = {
    "+1 for lambda, -1 for antilambda",
    "Minimum pT (GeV/c)",
    "Maximum pT (GeV/c)",
    "Maximum |Eta|",
    "Maximum DCA between daughters (cm)",
    "Minimum Cosine of Pointing Angle",
    "Minimum transverse radius (cm)",
    "Maximum transverse radius (cm)",
    "Maximum distance from primary vertex"}; ///< Helper information for the
                                             ///< different selections

}; // namespace femtoDream

template <o2::aod::femtodreamparticle::ParticleType part,
          o2::aod::femtodreamparticle::ParticleType daugh,
          typename cutContainerType>
void FemtoDreamV0Selection::init(HistogramRegistry* QAregistry, HistogramRegistry* Registry)
{
  if (QAregistry && Registry) {
    mHistogramRegistry = Registry;
    mQAHistogramRegistry = QAregistry;
    fillSelectionHistogram<part>();
    fillSelectionHistogram<daugh>();

    AxisSpec massAxisLambda = {600, 0.0f, 3.0f, "m_{#Lambda} (GeV/#it{c}^{2})"}; /// paramters for K0Short
    AxisSpec massAxisAntiLambda = {600, 0.0f, 3.0f,
                                   "m_{#bar{#Lambda}} (GeV/#it{c}^{2})"};

    /// \todo this should be an automatic check in the parent class, and the
    /// return type should be templated
    size_t nSelections = getNSelections();
    if (nSelections > 8 * sizeof(cutContainerType)) { // o2-linter: disable=magic-number (number represents bitsize)
      LOG(fatal) << "FemtoDreamV0Cuts: Number of selections to large for your "
                    "container - quitting!";
    }

    std::string folderName = static_cast<std::string>(
      o2::aod::femtodreamparticle::ParticleTypeName[part]);
    /// \todo initialize histograms for children tracks of v0s
    mQAHistogramRegistry->add((folderName + "/hPt").c_str(),
                              "; #it{p}_{T} (GeV/#it{c}); Entries", kTH1F,
                              {{1000, 0, 10}});
    mQAHistogramRegistry->add((folderName + "/hEta").c_str(), "; #eta; Entries",
                              kTH1F, {{1000, -1, 1}});
    mQAHistogramRegistry->add((folderName + "/hPhi").c_str(), "; #phi; Entries",
                              kTH1F, {{1000, 0, o2::constants::math::TwoPI}});
    mQAHistogramRegistry->add((folderName + "/hDaughDCA").c_str(),
                              "; DCA^{daugh} (cm); Entries", kTH1F,
                              {{1000, 0, 10}});
    mQAHistogramRegistry->add((folderName + "/hTransRadius").c_str(),
                              "; #it{r}_{xy} (cm); Entries", kTH1F,
                              {{1500, 0, 150}});
    mQAHistogramRegistry->add((folderName + "/hDecayVtxX").c_str(),
                              "; #it{Vtx}_{x} (cm); Entries", kTH1F,
                              {{2000, 0, 200}});
    mQAHistogramRegistry->add((folderName + "/hDecayVtxY").c_str(),
                              "; #it{Vtx}_{y} (cm)); Entries", kTH1F,
                              {{2000, 0, 200}});
    mQAHistogramRegistry->add((folderName + "/hDecayVtxZ").c_str(),
                              "; #it{Vtx}_{z} (cm); Entries", kTH1F,
                              {{2000, 0, 200}});
    mQAHistogramRegistry->add((folderName + "/hCPA").c_str(),
                              "; #it{cos #theta_{p}}; Entries", kTH1F,
                              {{1000, 0.9, 1.}});
    mQAHistogramRegistry->add((folderName + "/hCPAvsPt").c_str(),
                              "; #it{p}_{T} (GeV/#it{c}); #it{cos #theta_{p}}",
                              kTH2F, {{8, 0.3, 4.3}, {1000, 0.9, 1.}});
    mQAHistogramRegistry->add((folderName + "/hInvMassLambda").c_str(), "", kTH1F,
                              {massAxisLambda});
    mQAHistogramRegistry->add((folderName + "/hInvMassAntiLambda").c_str(), "",
                              kTH1F, {massAxisAntiLambda});
    mQAHistogramRegistry->add((folderName + "/hInvMassLambdaAntiLambda").c_str(),
                              "", kTH2F, {massAxisLambda, massAxisAntiLambda});

    posDaughTrack.init<daugh,
                       aod::femtodreamparticle::TrackType::kPosChild,
                       aod::femtodreamparticle::cutContainerType>(
      mQAHistogramRegistry, mHistogramRegistry);
    negDaughTrack.init<daugh,
                       aod::femtodreamparticle::TrackType::kNegChild,
                       aod::femtodreamparticle::cutContainerType>(
      mQAHistogramRegistry, mHistogramRegistry);

    mQAHistogramRegistry->add(folderName + "QA/hInvMassLambdaNoCuts", "No cuts", kTH1F,
                              {massAxisLambda});
    mQAHistogramRegistry->add(folderName + "QA/hInvMassLambdaInvMassCut",
                              "Invariant mass cut", kTH1F, {massAxisLambda});
    mQAHistogramRegistry->add(folderName + "QA/hInvMassLambdaPtMin", "Minimum Pt cut",
                              kTH1F, {massAxisLambda});
    mQAHistogramRegistry->add(folderName + "QA/hInvMassLambdaPtMax", "Maximum Pt cut",
                              kTH1F, {massAxisLambda});
    mQAHistogramRegistry->add(folderName + "QA/hInvMassLambdaEtaMax", "Maximum Eta cut",
                              kTH1F, {massAxisLambda});
    mQAHistogramRegistry->add(folderName + "QA/hInvMassLambdaDCAV0Daugh",
                              "V0-daughters DCA cut", kTH1F, {massAxisLambda});
    mQAHistogramRegistry->add(folderName + "QA/hInvMassLambdaCPA", "CPA cut", kTH1F,
                              {massAxisLambda});
    mQAHistogramRegistry->add(folderName + "QA/hInvMassLambdaTranRadMin",
                              "Minimum transverse radius cut", kTH1F,
                              {massAxisLambda});
    mQAHistogramRegistry->add(folderName + "QA/hInvMassLambdaTranRadMax",
                              "Maximum transverse radius cut", kTH1F,
                              {massAxisLambda});
    mQAHistogramRegistry->add(folderName + "QA/hInvMassLambdaDecVtxMax",
                              "Maximum distance on  decay vertex cut", kTH1F,
                              {massAxisLambda});
  }

  /// check whether the most open cuts are fulfilled - most of this should have
  /// already be done by the filters
  nPtV0MinSel = getNSelections(femto_dream_v0_selection::kV0pTMin);
  nPtV0MaxSel = getNSelections(femto_dream_v0_selection::kV0pTMax);
  nEtaV0MaxSel = getNSelections(femto_dream_v0_selection::kV0etaMax);
  nDCAV0DaughMax = getNSelections(femto_dream_v0_selection::kV0DCADaughMax);
  nCPAV0Min = getNSelections(femto_dream_v0_selection::kV0CPAMin);
  nTranRadV0Min = getNSelections(femto_dream_v0_selection::kV0TranRadMin);
  nTranRadV0Max = getNSelections(femto_dream_v0_selection::kV0TranRadMax);
  nDecVtxMax = getNSelections(femto_dream_v0_selection::kV0DecVtxMax);

  pTV0Min = getMinimalSelection(femto_dream_v0_selection::kV0pTMin,
                                femtoDreamSelection::kLowerLimit);
  pTV0Max = getMinimalSelection(femto_dream_v0_selection::kV0pTMax,
                                femtoDreamSelection::kUpperLimit);
  etaV0Max = getMinimalSelection(femto_dream_v0_selection::kV0etaMax,
                                 femtoDreamSelection::kAbsUpperLimit);
  dCAV0DaughMax = getMinimalSelection(femto_dream_v0_selection::kV0DCADaughMax,
                                      femtoDreamSelection::kUpperLimit);
  cPAV0Min = getMinimalSelection(femto_dream_v0_selection::kV0CPAMin,
                                 femtoDreamSelection::kLowerLimit);
  tranRadV0Min = getMinimalSelection(femto_dream_v0_selection::kV0TranRadMin,
                                     femtoDreamSelection::kLowerLimit);
  tranRadV0Max = getMinimalSelection(femto_dream_v0_selection::kV0TranRadMax,
                                     femtoDreamSelection::kUpperLimit);
  decVtxMax = getMinimalSelection(femto_dream_v0_selection::kV0DecVtxMax,
                                  femtoDreamSelection::kAbsUpperLimit);
}

template <typename C, typename V, typename T>
bool FemtoDreamV0Selection::isSelectedMinimal(C const& /*col*/, V const& v0,
                                              T const& posTrack,
                                              T const& negTrack)
{
  const auto signPos = posTrack.sign();
  const auto signNeg = negTrack.sign();
  if (signPos < 0 || signNeg > 0) {
    LOG(warn) << "Something wrong in isSelectedMinimal";
    LOG(warn) << "ERROR - Wrong sign for V0 daughters";
  }
  // asfaf
  const float pT = v0.pt();
  const float eta = v0.eta();
  const std::vector<float> decVtx = {v0.x(), v0.y(), v0.z()};
  const float tranRad = v0.v0radius();
  const float dcaDaughv0 = v0.dcaV0daughters();
  const float cpav0 = v0.v0cosPA();

  const float invMassLambda = v0.mLambda();
  const float invMassAntiLambda = v0.mAntiLambda();

  const float invMassKaon = v0.mK0Short();

  if (fMotherIsLambda) { /// Lambda
    if ((invMassLambda < fInvMassLowLimit || invMassLambda > fInvMassUpLimit) &&
        (invMassAntiLambda < fInvMassLowLimit ||
         invMassAntiLambda > fInvMassUpLimit)) {
      return false;
    }
    if (fRejectKaon) {
      if (invMassKaon > fInvMassKaonLowLimit &&
          invMassKaon < fInvMassKaonUpLimit) {
        return false;
      }
    }
  } else { /// K0Short
    if ((invMassKaon < fInvMassKaonLowLimit || invMassKaon > fInvMassKaonUpLimit)) {
      return false;
    }
    if (fRejectLambda) {
      if ((invMassLambda > fInvMassLowLimit &&
           invMassLambda < fInvMassUpLimit) ||
          (invMassAntiLambda > fInvMassLowLimit &&
           invMassAntiLambda < fInvMassUpLimit)) {
        return false;
      }
    }
  }

  if (nPtV0MinSel > 0 && pT < pTV0Min) {
    return false;
  }
  if (nPtV0MaxSel > 0 && pT > pTV0Max) {
    return false;
  }
  if (nEtaV0MaxSel > 0 && std::abs(eta) > etaV0Max) {
    return false;
  }
  if (nDCAV0DaughMax > 0 && dcaDaughv0 > dCAV0DaughMax) {
    return false;
  }
  if (nCPAV0Min > 0 && cpav0 < cPAV0Min) {
    return false;
  }
  if (nTranRadV0Min > 0 && tranRad < tranRadV0Min) {
    return false;
  }
  if (nTranRadV0Max > 0 && tranRad > tranRadV0Max) {
    return false;
  }
  for (size_t i = 0; i < decVtx.size(); i++) {
    if (nDecVtxMax > 0 && decVtx.at(i) > decVtxMax) {
      return false;
    }
  }
  if (!posDaughTrack.isSelectedMinimal(posTrack)) {
    return false;
  }
  if (!negDaughTrack.isSelectedMinimal(negTrack)) {
    return false;
  }

  // check that track combinations for V0 or antiV0 would be fulfilling PID
  int nSigmaPIDMax = posDaughTrack.getSigmaPIDMax();
  if (fMotherIsLambda) { /// Lambda
    // antiV0
    auto nSigmaPrNeg = negTrack.tpcNSigmaPr();
    auto nSigmaPiPos = posTrack.tpcNSigmaPi();
    // v0
    auto nSigmaPiNeg = negTrack.tpcNSigmaPi();
    auto nSigmaPrPos = posTrack.tpcNSigmaPr();
    if (!(std::abs(nSigmaPrNeg - nSigmaPIDOffsetTPC) < nSigmaPIDMax &&
          std::abs(nSigmaPiPos - nSigmaPIDOffsetTPC) < nSigmaPIDMax) &&
        !(std::abs(nSigmaPrPos - nSigmaPIDOffsetTPC) < nSigmaPIDMax &&
          std::abs(nSigmaPiNeg - nSigmaPIDOffsetTPC) < nSigmaPIDMax)) {
      return false;
    }
  } else { /// K0SHort
    auto nSigmaPiNeg = negTrack.tpcNSigmaPi();
    auto nSigmaPiPos = posTrack.tpcNSigmaPi();
    if (!(std::abs(nSigmaPiPos - nSigmaPIDOffsetTPC) < nSigmaPIDMax &&
          std::abs(nSigmaPiNeg - nSigmaPIDOffsetTPC) < nSigmaPIDMax)) {
      return false;
    }
  }

  return true;
}

template <o2::aod::femtodreamparticle::ParticleType part,
          typename C, typename V, typename T>
void FemtoDreamV0Selection::fillLambdaQA(C const& /*col*/, V const& v0,
                                         T const& posTrack, T const& negTrack)
{
  const auto signPos = posTrack.sign();
  const auto signNeg = negTrack.sign();
  if (signPos < 0 || signNeg > 0) {
    LOG(warn) << "Something wrong in isSelectedMinimal";
    LOG(warn) << "ERROR - Wrong sign for V0 daughters";
  }
  const float pT = v0.pt();
  const float eta = v0.eta();
  const std::vector<float> decVtx = {v0.x(), v0.y(), v0.z()};
  const float tranRad = v0.v0radius();
  const float dcaDaughv0 = v0.dcaV0daughters();
  const float cpav0 = v0.v0cosPA();

  const float invMassLambda = v0.mLambda();
  const float invMassKaon = v0.mK0Short();

  float fillMass = 0.;
  if (fMotherIsLambda) {
    fillMass = v0.mLambda();
  } else {
    fillMass = v0.mK0Short();
  }

  mQAHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) + HIST("QA/hInvMassLambdaNoCuts"), fillMass);

  if (fMotherIsLambda) { /// Lambda
    if (invMassLambda > fInvMassLowLimit && invMassLambda < fInvMassUpLimit) {
      mQAHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) + HIST("QA/hInvMassLambdaInvMassCut"),
                                 v0.mLambda());
    }
  } else { /// K0Short
    if (invMassKaon > fInvMassKaonLowLimit && invMassKaon < fInvMassKaonUpLimit) {
      mQAHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) + HIST("QA/hInvMassLambdaInvMassCut"),
                                 v0.mK0Short());
    }
  }

  if (pT > pTV0Min) {
    mQAHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) + HIST("QA/hInvMassLambdaPtMin"),
                               fillMass);
  }
  if (pT < pTV0Max) {
    mQAHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) + HIST("QA/hInvMassLambdaPtMax"),
                               fillMass);
  }
  if (std::abs(eta) < etaV0Max) {
    mQAHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) + HIST("QA/hInvMassLambdaEtaMax"),
                               fillMass);
  }
  if (dcaDaughv0 < dCAV0DaughMax) {
    mQAHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) + HIST("QA/hInvMassLambdaDCAV0Daugh"),
                               fillMass);
  }
  if (cpav0 > cPAV0Min) {
    mQAHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) + HIST("QA/hInvMassLambdaCPA"), fillMass);
  }
  if (tranRad > tranRadV0Min) {
    mQAHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) + HIST("QA/hInvMassLambdaTranRadMin"),
                               fillMass);
  }
  if (tranRad < tranRadV0Max) {
    mQAHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) + HIST("QA/hInvMassLambdaTranRadMax"),
                               fillMass);
  }
  bool write = true;
  for (size_t i = 0; i < decVtx.size(); i++) {
    write = write && (decVtx.at(i) < decVtxMax);
  }
  if (write) {
    mQAHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) + HIST("QA/hInvMassLambdaDecVtxMax"),
                               fillMass);
  }
}

/// the CosPA of V0 needs as argument the posXYZ of collisions vertex so we need
/// to pass the collsion as well
template <typename cutContainerType, typename C, typename V, typename T>
std::array<cutContainerType, 5>
  FemtoDreamV0Selection::getCutContainer(C const& /*col*/, V const& v0, T const& posTrack, T const& negTrack)
{
  auto outputPosTrack = posDaughTrack.getCutContainer<false, cutContainerType>(posTrack, v0.positivept(), v0.positiveeta(), v0.dcapostopv());
  auto outputNegTrack = negDaughTrack.getCutContainer<false, cutContainerType>(negTrack, v0.negativept(), v0.negativeeta(), v0.dcanegtopv());
  cutContainerType output = 0;
  size_t counter = 0;

  float sign = 0.;
  /// Lambda
  if (fMotherIsLambda) {
    auto lambdaMassNominal = o2::constants::physics::MassLambda;
    auto lambdaMassHypothesis = v0.mLambda();
    auto antiLambdaMassHypothesis = v0.mAntiLambda();
    auto diffLambda = std::abs(lambdaMassNominal - lambdaMassHypothesis);
    auto diffAntiLambda = std::abs(antiLambdaMassHypothesis - lambdaMassHypothesis);

    int nSigmaPIDMax = posDaughTrack.getSigmaPIDMax();
    auto nSigmaPrNeg = negTrack.tpcNSigmaPr();
    auto nSigmaPiPos = posTrack.tpcNSigmaPi();
    auto nSigmaPiNeg = negTrack.tpcNSigmaPi();
    auto nSigmaPrPos = posTrack.tpcNSigmaPr();
    // check the mass and the PID of daughters
    if (std::abs(nSigmaPrNeg - nSigmaPIDOffsetTPC) < nSigmaPIDMax && std::abs(nSigmaPiPos - nSigmaPIDOffsetTPC) < nSigmaPIDMax && diffAntiLambda > diffLambda) {
      sign = -1.;
    } else if (std::abs(nSigmaPrPos - nSigmaPIDOffsetTPC) < nSigmaPIDMax && std::abs(nSigmaPiNeg - nSigmaPIDOffsetTPC) < nSigmaPIDMax && diffAntiLambda < diffLambda) {
      sign = 1.;
    } else {
      // if it happens that none of these are true, ignore the invariant mass
      if (std::abs(nSigmaPrNeg - nSigmaPIDOffsetTPC) < nSigmaPIDMax && std::abs(nSigmaPiPos - nSigmaPIDOffsetTPC) < nSigmaPIDMax) {
        sign = -1.;
      } else if (std::abs(nSigmaPrPos - nSigmaPIDOffsetTPC) < nSigmaPIDMax && std::abs(nSigmaPiNeg - nSigmaPIDOffsetTPC) < nSigmaPIDMax) {
        sign = 1.;
      }
    }
  } else {
    sign = 1.; // for the K0Short arbitrarily set the sign to 1
  }

  const auto pT = v0.pt();
  const auto eta = v0.eta();
  const auto tranRad = v0.v0radius();
  const auto dcaDaughv0 = v0.dcaV0daughters();
  const auto cpav0 = v0.v0cosPA();
  const std::vector<float> decVtx = {v0.x(), v0.y(), v0.z()};

  float observable = 0.;
  for (auto& sel : mSelections) { // o2-linter: disable=const-ref-in-for-loop (femtoDreamObjectSelection has no const getter)
    const auto selVariable = sel.getSelectionVariable();
    if (selVariable == femto_dream_v0_selection::kV0DecVtxMax) {
      for (size_t i = 0; i < decVtx.size(); ++i) {
        auto decVtxValue = decVtx.at(i);
        sel.checkSelectionSetBit(decVtxValue, output, counter, nullptr);
      }
    } else {
      switch (selVariable) {
        case (femto_dream_v0_selection::kV0Sign):
          observable = sign;
          break;
        case (femto_dream_v0_selection::kV0pTMin):
          observable = pT;
          break;
        case (femto_dream_v0_selection::kV0pTMax):
          observable = pT;
          break;
        case (femto_dream_v0_selection::kV0etaMax):
          observable = eta;
          break;
        case (femto_dream_v0_selection::kV0DCADaughMax):
          observable = dcaDaughv0;
          break;
        case (femto_dream_v0_selection::kV0CPAMin):
          observable = cpav0;
          break;
        case (femto_dream_v0_selection::kV0TranRadMin):
          observable = tranRad;
          break;
        case (femto_dream_v0_selection::kV0TranRadMax):
          observable = tranRad;
          break;
        case (femto_dream_v0_selection::kV0DecVtxMax):
          break;
      }
      sel.checkSelectionSetBit(observable, output, counter, nullptr);
    }
  }
  return {
    output,
    outputPosTrack.at(femtoDreamTrackSelection::TrackContainerPosition::kCuts),
    outputPosTrack.at(femtoDreamTrackSelection::TrackContainerPosition::kPID),
    outputNegTrack.at(femtoDreamTrackSelection::TrackContainerPosition::kCuts),
    outputNegTrack.at(femtoDreamTrackSelection::TrackContainerPosition::kPID)};
}

template <o2::aod::femtodreamparticle::ParticleType part,
          o2::aod::femtodreamparticle::ParticleType daugh, typename C,
          typename V, typename T>
void FemtoDreamV0Selection::fillQA(C const& /*col*/, V const& v0, T const& posTrack,
                                   T const& negTrack)
{
  if (mQAHistogramRegistry) {
    mQAHistogramRegistry->fill(
      HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) +
        HIST("/hPt"),
      v0.pt());
    mQAHistogramRegistry->fill(
      HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) +
        HIST("/hEta"),
      v0.eta());
    mQAHistogramRegistry->fill(
      HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) +
        HIST("/hPhi"),
      v0.phi());
    mQAHistogramRegistry->fill(
      HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) +
        HIST("/hDaughDCA"),
      v0.dcaV0daughters());
    mQAHistogramRegistry->fill(
      HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) +
        HIST("/hTransRadius"),
      v0.v0radius());
    mQAHistogramRegistry->fill(
      HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) +
        HIST("/hDecayVtxX"),
      v0.x());
    mQAHistogramRegistry->fill(
      HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) +
        HIST("/hDecayVtxY"),
      v0.y());
    mQAHistogramRegistry->fill(
      HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) +
        HIST("/hDecayVtxZ"),
      v0.z());
    mQAHistogramRegistry->fill(
      HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) +
        HIST("/hCPA"),
      v0.v0cosPA());
    mQAHistogramRegistry->fill(
      HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) +
        HIST("/hCPAvsPt"),
      v0.pt(), v0.v0cosPA());

    if (fMotherIsLambda) { /// Lambda
      mQAHistogramRegistry->fill(
        HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) +
          HIST("/hInvMassLambda"),
        v0.mLambda());
      mQAHistogramRegistry->fill(
        HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) +
          HIST("/hInvMassAntiLambda"),
        v0.mAntiLambda());
      mQAHistogramRegistry->fill(
        HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) +
          HIST("/hInvMassLambdaAntiLambda"),
        v0.mLambda(), v0.mAntiLambda());
    } else { /// K0Short
      mQAHistogramRegistry->fill(
        HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) +
          HIST("/hInvMassLambda"),
        v0.mK0Short());
    }
  }

  posDaughTrack.fillQA<daugh,
                       aod::femtodreamparticle::TrackType::kPosChild>(posTrack);
  negDaughTrack.fillQA<daugh,
                       aod::femtodreamparticle::TrackType::kNegChild>(negTrack);
}

} // namespace o2::analysis::femtoDream

#endif // PWGCF_FEMTODREAM_CORE_FEMTODREAMV0SELECTIONK0SHORT_H_
