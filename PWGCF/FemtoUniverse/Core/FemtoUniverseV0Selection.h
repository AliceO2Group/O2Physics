// Copyright 2019-2022 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file FemtoUniverseV0Selection.h
/// \brief Definition of the FemtoUniverseV0Selection
/// \author Valentina Mantovani Sarti, TU München valentina.mantovani-sarti@tum.de
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de
/// \author Luca Barioglio, TU München, luca.barioglio@cern.ch
/// \author Zuzanna Chochulska, WUT Warsaw & CTU Prague, zchochul@cern.ch

#ifndef PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSEV0SELECTION_H_
#define PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSEV0SELECTION_H_

#include <iostream>
#include <string>
#include <vector>

#include <TDatabasePDG.h> // FIXME

#include "PWGCF/FemtoUniverse/Core/FemtoUniverseObjectSelection.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseSelection.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseTrackSelection.h"

#include "Common/Core/RecoDecay.h"
#include "Framework/HistogramRegistry.h"
#include "ReconstructionDataFormats/PID.h"

using namespace o2::framework;

namespace o2::analysis::femtoUniverse
{
namespace femtoUniverseV0Selection
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
  kV0DecVtxMax
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

} // namespace femtoUniverseV0Selection

/// \class FemtoUniverseV0Selection
/// \brief Cut class to contain and execute all cuts applied to V0s
class FemtoUniverseV0Selection
  : public FemtoUniverseObjectSelection<float, femtoUniverseV0Selection::V0Sel>
{
 public:
  FemtoUniverseV0Selection()
    : nPtV0MinSel(0), nPtV0MaxSel(0), nEtaV0MaxSel(0), nDCAV0DaughMax(0), nCPAV0Min(0), nTranRadV0Min(0), nTranRadV0Max(0), nDecVtxMax(0), pTV0Min(9999999.), pTV0Max(-9999999.), etaV0Max(-9999999.), DCAV0DaughMax(-9999999.), CPAV0Min(9999999.), TranRadV0Min(9999999.), TranRadV0Max(-9999999.), DecVtxMax(-9999999.), fInvMassLowLimit(1.05), fInvMassUpLimit(1.3), fRejectKaon(false), fInvMassKaonLowLimit(0.48), fInvMassKaonUpLimit(0.515), nSigmaPIDOffsetTPC(0.) {}
  /// Initializes histograms for the task
  template <o2::aod::femtouniverseparticle::ParticleType part,
            o2::aod::femtouniverseparticle::ParticleType daugh,
            typename cutContainerType>
  void init(HistogramRegistry* registry);

  template <typename C, typename V, typename T>
  bool isSelectedMinimal(C const& col, V const& v0, T const& posTrack,
                         T const& negTrack);

  template <typename C, typename V, typename T>
  void fillLambdaQA(C const& col, V const& v0, T const& posTrack,
                    T const& negTrack);

  /// \todo for the moment the PID of the tracks is factored out into a separate
  /// field, hence 5 values in total \\ASK: what does it mean?
  template <typename cutContainerType, typename C, typename V, typename T>
  std::array<cutContainerType, 5> getCutContainer(C const& col, V const& v0,
                                                  T const& posTrack,
                                                  T const& negTrack);

  template <o2::aod::femtouniverseparticle::ParticleType part,
            o2::aod::femtouniverseparticle::ParticleType daugh, typename C,
            typename V, typename T>
  void fillQA(C const& col, V const& v0, T const& posTrack, T const& negTrack);

  template <typename T1, typename T2>
  void setChildCuts(femtoUniverseV0Selection::ChildTrackType child, T1 selVal,
                    T2 selVar, femtoUniverseSelection::SelectionType selType)
  {
    if (child == femtoUniverseV0Selection::kPosTrack) {
      PosDaughTrack.setSelection(selVal, selVar, selType);
    } else if (child == femtoUniverseV0Selection::kNegTrack) {
      NegDaughTrack.setSelection(selVal, selVar, selType);
    }
  }
  template <typename T>
  void setChildPIDSpecies(femtoUniverseV0Selection::ChildTrackType child,
                          T& pids)
  {
    if (child == femtoUniverseV0Selection::kPosTrack) {
      PosDaughTrack.setPIDSpecies(pids);
    } else if (child == femtoUniverseV0Selection::kNegTrack) {
      NegDaughTrack.setPIDSpecies(pids);
    }
  }

  /// Helper function to obtain the name of a given selection criterion for consistent naming of the configurables
  /// \param iSel Track selection variable to be examined
  /// \param prefix Additional prefix for the name of the configurable
  /// \param suffix Additional suffix for the name of the configurable
  static std::string getSelectionName(femtoUniverseV0Selection::V0Sel iSel,
                                      std::string_view prefix = "",
                                      std::string_view suffix = "")
  {
    std::string outString = static_cast<std::string>(prefix);
    outString += static_cast<std::string>(mSelectionNames[iSel]);
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
                         static_cast<std::string>(mSelectionNames[index]);
      std::string_view cmp{comp};
      if (obs.compare(cmp) == 0)
        return index;
    }
    LOGF(info, "Variable %s not found", obs);
    return -1;
  }

  /// Helper function to obtain the type of a given selection variable for consistent naming of the configurables
  /// \param iSel V0 selection variable whose type is returned
  static femtoUniverseSelection::SelectionType
    getSelectionType(femtoUniverseV0Selection::V0Sel iSel)
  {
    return mSelectionTypes[iSel];
  }

  /// Helper function to obtain the helper string of a given selection criterion
  /// for consistent description of the configurables
  /// \param iSel Track selection variable to be examined
  /// \param prefix Additional prefix for the output of the configurable
  static std::string getSelectionHelper(femtoUniverseV0Selection::V0Sel iSel,
                                        std::string_view prefix = "")
  {
    std::string outString = static_cast<std::string>(prefix);
    outString += static_cast<std::string>(mSelectionHelper[iSel]);
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

  void setChildRejectNotPropagatedTracks(femtoUniverseV0Selection::ChildTrackType child, bool reject)
  {
    if (child == femtoUniverseV0Selection::kPosTrack) {
      PosDaughTrack.setRejectNotPropagatedTracks(reject);
    } else if (child == femtoUniverseV0Selection::kNegTrack) {
      NegDaughTrack.setRejectNotPropagatedTracks(reject);
    }
  }

  void setChildnSigmaPIDOffset(femtoUniverseV0Selection::ChildTrackType child, float offsetTPC, float offsetTOF)
  {
    if (child == femtoUniverseV0Selection::kPosTrack) {
      PosDaughTrack.setnSigmaPIDOffset(offsetTPC, offsetTOF);
    } else if (child == femtoUniverseV0Selection::kNegTrack) {
      NegDaughTrack.setnSigmaPIDOffset(offsetTPC, offsetTOF);
    }
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
  float DCAV0DaughMax;
  float CPAV0Min;
  float TranRadV0Min;
  float TranRadV0Max;
  float DecVtxMax;

  float fInvMassLowLimit;
  float fInvMassUpLimit;

  bool fRejectKaon;
  float fInvMassKaonLowLimit;
  float fInvMassKaonUpLimit;

  float nSigmaPIDOffsetTPC;

  FemtoUniverseTrackSelection PosDaughTrack;
  FemtoUniverseTrackSelection NegDaughTrack;

  static constexpr int kNv0Selection = 9;

  static constexpr std::string_view mSelectionNames[kNv0Selection] = {
    "Sign", "PtMin", "PtMax", "EtaMax", "DCAdaughMax", "CPAMin",
    "TranRadMin", "TranRadMax", "DecVecMax"}; ///< Name of the different
                                              ///< selections

  static constexpr femtoUniverseSelection::SelectionType
    mSelectionTypes[kNv0Selection]{
      femtoUniverseSelection::kEqual,
      femtoUniverseSelection::kLowerLimit,
      femtoUniverseSelection::kUpperLimit,
      femtoUniverseSelection::kUpperLimit,
      femtoUniverseSelection::kUpperLimit,
      femtoUniverseSelection::kLowerLimit,
      femtoUniverseSelection::kLowerLimit,
      femtoUniverseSelection::kUpperLimit,
      femtoUniverseSelection::kUpperLimit}; ///< Map to match a variable with
                                            ///< its type

  static constexpr std::string_view mSelectionHelper[kNv0Selection] = {
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

}; // namespace femtoUniverse

template <o2::aod::femtouniverseparticle::ParticleType part,
          o2::aod::femtouniverseparticle::ParticleType daugh,
          typename cutContainerType>
void FemtoUniverseV0Selection::init(HistogramRegistry* registry)
{
  if (registry) {
    mHistogramRegistry = registry;
    fillSelectionHistogram<part>();
    fillSelectionHistogram<daugh>();

    AxisSpec massAxisLambda = {600, 0.0f, 3.0f, "m_{#Lambda} (GeV/#it{c}^{2})"};
    AxisSpec massAxisAntiLambda = {600, 0.0f, 3.0f,
                                   "m_{#bar{#Lambda}} (GeV/#it{c}^{2})"};

    /// \todo this should be an automatic check in the parent class, and the
    /// return type should be templated
    size_t nSelections = getNSelections();
    if (nSelections > 8 * sizeof(cutContainerType)) {
      LOG(fatal) << "FemtoUniverseV0Cuts: Number of selections to large for your "
                    "container - quitting!";
    }
    std::string folderName = static_cast<std::string>(
      o2::aod::femtouniverseparticle::ParticleTypeName[part]);
    /// \todo initialize histograms for children tracks of v0s
    mHistogramRegistry->add((folderName + "/hPt").c_str(),
                            "; #it{p}_{T} (GeV/#it{c}); Entries", kTH1F,
                            {{1000, 0, 10}});
    mHistogramRegistry->add((folderName + "/hEta").c_str(), "; #eta; Entries",
                            kTH1F, {{1000, -1, 1}});
    mHistogramRegistry->add((folderName + "/hPhi").c_str(), "; #phi; Entries",
                            kTH1F, {{1000, 0, 2. * M_PI}});
    mHistogramRegistry->add((folderName + "/hDaughDCA").c_str(),
                            "; DCA^{daugh} (cm); Entries", kTH1F,
                            {{1000, 0, 10}});
    mHistogramRegistry->add((folderName + "/hTransRadius").c_str(),
                            "; #it{r}_{xy} (cm); Entries", kTH1F,
                            {{1500, 0, 150}});
    mHistogramRegistry->add((folderName + "/hDecayVtxX").c_str(),
                            "; #it{Vtx}_{x} (cm); Entries", kTH1F,
                            {{2000, 0, 200}});
    mHistogramRegistry->add((folderName + "/hDecayVtxY").c_str(),
                            "; #it{Vtx}_{y} (cm)); Entries", kTH1F,
                            {{2000, 0, 200}});
    mHistogramRegistry->add((folderName + "/hDecayVtxZ").c_str(),
                            "; #it{Vtx}_{z} (cm); Entries", kTH1F,
                            {{2000, 0, 200}});
    mHistogramRegistry->add((folderName + "/hCPA").c_str(),
                            "; #it{cos #theta_{p}}; Entries", kTH1F,
                            {{1000, 0.9, 1.}});
    mHistogramRegistry->add((folderName + "/hCPAvsPt").c_str(),
                            "; #it{p}_{T} (GeV/#it{c}); #it{cos #theta_{p}}",
                            kTH2F, {{8, 0.3, 4.3}, {1000, 0.9, 1.}});
    mHistogramRegistry->add((folderName + "/hInvMassLambda").c_str(), "", kTH1F,
                            {massAxisLambda});
    mHistogramRegistry->add((folderName + "/hInvMassAntiLambda").c_str(), "",
                            kTH1F, {massAxisAntiLambda});
    mHistogramRegistry->add((folderName + "/hInvMassLambdaAntiLambda").c_str(),
                            "", kTH2F, {massAxisLambda, massAxisAntiLambda});

    PosDaughTrack.init<aod::femtouniverseparticle::ParticleType::kV0Child,
                       aod::femtouniverseparticle::TrackType::kPosChild,
                       aod::femtouniverseparticle::cutContainerType>(
      mHistogramRegistry);
    NegDaughTrack.init<aod::femtouniverseparticle::ParticleType::kV0Child,
                       aod::femtouniverseparticle::TrackType::kNegChild,
                       aod::femtouniverseparticle::cutContainerType>(
      mHistogramRegistry);

    mHistogramRegistry->add("LambdaQA/hInvMassLambdaNoCuts", "No cuts", kTH1F,
                            {massAxisLambda});
    mHistogramRegistry->add("LambdaQA/hInvMassLambdaInvMassCut",
                            "Invariant mass cut", kTH1F, {massAxisLambda});
    mHistogramRegistry->add("LambdaQA/hInvMassLambdaPtMin", "Minimum Pt cut",
                            kTH1F, {massAxisLambda});
    mHistogramRegistry->add("LambdaQA/hInvMassLambdaPtMax", "Maximum Pt cut",
                            kTH1F, {massAxisLambda});
    mHistogramRegistry->add("LambdaQA/hInvMassLambdaEtaMax", "Maximum Eta cut",
                            kTH1F, {massAxisLambda});
    mHistogramRegistry->add("LambdaQA/hInvMassLambdaDCAV0Daugh",
                            "V0-daughters DCA cut", kTH1F, {massAxisLambda});
    mHistogramRegistry->add("LambdaQA/hInvMassLambdaCPA", "CPA cut", kTH1F,
                            {massAxisLambda});
    mHistogramRegistry->add("LambdaQA/hInvMassLambdaTranRadMin",
                            "Minimum transverse radius cut", kTH1F,
                            {massAxisLambda});
    mHistogramRegistry->add("LambdaQA/hInvMassLambdaTranRadMax",
                            "Maximum transverse radius cut", kTH1F,
                            {massAxisLambda});
    mHistogramRegistry->add("LambdaQA/hInvMassLambdaDecVtxMax",
                            "Maximum distance on  decay vertex cut", kTH1F,
                            {massAxisLambda});
  }
  /// check whether the most open cuts are fulfilled - most of this should have
  /// already be done by the filters
  nPtV0MinSel = getNSelections(femtoUniverseV0Selection::kV0pTMin);
  nPtV0MaxSel = getNSelections(femtoUniverseV0Selection::kV0pTMax);
  nEtaV0MaxSel = getNSelections(femtoUniverseV0Selection::kV0etaMax);
  nDCAV0DaughMax = getNSelections(femtoUniverseV0Selection::kV0DCADaughMax);
  nCPAV0Min = getNSelections(femtoUniverseV0Selection::kV0CPAMin);
  nTranRadV0Min = getNSelections(femtoUniverseV0Selection::kV0TranRadMin);
  nTranRadV0Max = getNSelections(femtoUniverseV0Selection::kV0TranRadMax);
  nDecVtxMax = getNSelections(femtoUniverseV0Selection::kV0DecVtxMax);

  pTV0Min = getMinimalSelection(femtoUniverseV0Selection::kV0pTMin,
                                femtoUniverseSelection::kLowerLimit);
  pTV0Max = getMinimalSelection(femtoUniverseV0Selection::kV0pTMax,
                                femtoUniverseSelection::kUpperLimit);
  etaV0Max = getMinimalSelection(femtoUniverseV0Selection::kV0etaMax,
                                 femtoUniverseSelection::kAbsUpperLimit);
  DCAV0DaughMax = getMinimalSelection(femtoUniverseV0Selection::kV0DCADaughMax,
                                      femtoUniverseSelection::kUpperLimit);
  CPAV0Min = getMinimalSelection(femtoUniverseV0Selection::kV0CPAMin,
                                 femtoUniverseSelection::kLowerLimit);
  TranRadV0Min = getMinimalSelection(femtoUniverseV0Selection::kV0TranRadMin,
                                     femtoUniverseSelection::kLowerLimit);
  TranRadV0Max = getMinimalSelection(femtoUniverseV0Selection::kV0TranRadMax,
                                     femtoUniverseSelection::kUpperLimit);
  DecVtxMax = getMinimalSelection(femtoUniverseV0Selection::kV0DecVtxMax,
                                  femtoUniverseSelection::kAbsUpperLimit);
}

template <typename C, typename V, typename T>
bool FemtoUniverseV0Selection::isSelectedMinimal(C const& col, V const& v0,
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

  if ((invMassLambda < fInvMassLowLimit || invMassLambda > fInvMassUpLimit) &&
      (invMassAntiLambda < fInvMassLowLimit ||
       invMassAntiLambda > fInvMassUpLimit)) {
    return false;
  }
  if (fRejectKaon) {
    const float invMassKaon = v0.mK0Short();
    if (invMassKaon > fInvMassKaonLowLimit &&
        invMassKaon < fInvMassKaonUpLimit) {
      return false;
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
  if (nDCAV0DaughMax > 0 && dcaDaughv0 > DCAV0DaughMax) {
    return false;
  }
  if (nCPAV0Min > 0 && cpav0 < CPAV0Min) {
    return false;
  }
  if (nTranRadV0Min > 0 && tranRad < TranRadV0Min) {
    return false;
  }
  if (nTranRadV0Max > 0 && tranRad > TranRadV0Max) {
    return false;
  }
  for (size_t i = 0; i < decVtx.size(); i++) {
    if (nDecVtxMax > 0 && decVtx.at(i) > DecVtxMax) {
      return false;
    }
  }
  if (!PosDaughTrack.isSelectedMinimal(posTrack)) {
    return false;
  }
  if (!NegDaughTrack.isSelectedMinimal(negTrack)) {
    return false;
  }

  // check that track combinations for V0 or antiV0 would be fulfilling PID
  float nSigmaPIDMax = PosDaughTrack.getSigmaPIDMax();
  // antiV0
  auto nSigmaPrNeg = negTrack.tpcNSigmaPr();
  auto nSigmaPiPos = posTrack.tpcNSigmaPi();
  // v0
  auto nSigmaPiNeg = negTrack.tpcNSigmaPi();
  auto nSigmaPrPos = posTrack.tpcNSigmaPr();
  if (!(abs(nSigmaPrNeg - nSigmaPIDOffsetTPC) < nSigmaPIDMax &&
        abs(nSigmaPiPos - nSigmaPIDOffsetTPC) < nSigmaPIDMax) &&
      !(abs(nSigmaPrPos - nSigmaPIDOffsetTPC) < nSigmaPIDMax &&
        abs(nSigmaPiNeg - nSigmaPIDOffsetTPC) < nSigmaPIDMax)) {
    return false;
  }

  return true;
}

template <typename C, typename V, typename T>
void FemtoUniverseV0Selection::fillLambdaQA(C const& col, V const& v0,
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

  mHistogramRegistry->fill(HIST("LambdaQA/hInvMassLambdaNoCuts"), v0.mLambda());

  if (invMassLambda > fInvMassLowLimit && invMassLambda < fInvMassUpLimit) {
    mHistogramRegistry->fill(HIST("LambdaQA/hInvMassLambdaInvMassCut"),
                             v0.mLambda());
  }

  if (pT > pTV0Min) {
    mHistogramRegistry->fill(HIST("LambdaQA/hInvMassLambdaPtMin"),
                             v0.mLambda());
  }
  if (pT < pTV0Max) {
    mHistogramRegistry->fill(HIST("LambdaQA/hInvMassLambdaPtMax"),
                             v0.mLambda());
  }
  if (std::abs(eta) < etaV0Max) {
    mHistogramRegistry->fill(HIST("LambdaQA/hInvMassLambdaEtaMax"),
                             v0.mLambda());
  }
  if (dcaDaughv0 < DCAV0DaughMax) {
    mHistogramRegistry->fill(HIST("LambdaQA/hInvMassLambdaDCAV0Daugh"),
                             v0.mLambda());
  }
  if (cpav0 > CPAV0Min) {
    mHistogramRegistry->fill(HIST("LambdaQA/hInvMassLambdaCPA"), v0.mLambda());
  }
  if (tranRad > TranRadV0Min) {
    mHistogramRegistry->fill(HIST("LambdaQA/hInvMassLambdaTranRadMin"),
                             v0.mLambda());
  }
  if (tranRad < TranRadV0Max) {
    mHistogramRegistry->fill(HIST("LambdaQA/hInvMassLambdaTranRadMax"),
                             v0.mLambda());
  }
  bool write = true;
  for (size_t i = 0; i < decVtx.size(); i++) {
    write = write && (decVtx.at(i) < DecVtxMax);
  }
  if (write) {
    mHistogramRegistry->fill(HIST("LambdaQA/hInvMassLambdaDecVtxMax"),
                             v0.mLambda());
  }
}

/// the CosPA of V0 needs as argument the posXYZ of collisions vertex so we need
/// to pass the collsion as well
template <typename cutContainerType, typename C, typename V, typename T>
std::array<cutContainerType, 5>
  FemtoUniverseV0Selection::getCutContainer(C const& col, V const& v0, T const& posTrack, T const& negTrack)
{
  auto outputPosTrack = PosDaughTrack.getCutContainer<cutContainerType>(posTrack);
  auto outputNegTrack = NegDaughTrack.getCutContainer<cutContainerType>(negTrack);
  cutContainerType output = 0;
  size_t counter = 0;

  auto lambdaMassNominal = TDatabasePDG::Instance()->GetParticle(3122)->Mass(); // FIXME: Get from the common header
  auto lambdaMassHypothesis = v0.mLambda();
  auto antiLambdaMassHypothesis = v0.mAntiLambda();
  auto diffLambda = abs(lambdaMassNominal - lambdaMassHypothesis);
  auto diffAntiLambda = abs(antiLambdaMassHypothesis - lambdaMassHypothesis);

  float sign = 0.;
  float nSigmaPIDMax = PosDaughTrack.getSigmaPIDMax();
  auto nSigmaPrNeg = negTrack.tpcNSigmaPr();
  auto nSigmaPiPos = posTrack.tpcNSigmaPi();
  auto nSigmaPiNeg = negTrack.tpcNSigmaPi();
  auto nSigmaPrPos = posTrack.tpcNSigmaPr();
  // check the mass and the PID of daughters
  if (abs(nSigmaPrNeg - nSigmaPIDOffsetTPC) < nSigmaPIDMax && abs(nSigmaPiPos - nSigmaPIDOffsetTPC) < nSigmaPIDMax && diffAntiLambda > diffLambda) {
    sign = -1.;
  } else if (abs(nSigmaPrPos - nSigmaPIDOffsetTPC) < nSigmaPIDMax && abs(nSigmaPiNeg - nSigmaPIDOffsetTPC) < nSigmaPIDMax && diffAntiLambda < diffLambda) {
    sign = 1.;
  } else {
    // if it happens that none of these are true, ignore the invariant mass
    if (abs(nSigmaPrNeg - nSigmaPIDOffsetTPC) < nSigmaPIDMax && abs(nSigmaPiPos - nSigmaPIDOffsetTPC) < nSigmaPIDMax) {
      sign = -1.;
    } else if (abs(nSigmaPrPos - nSigmaPIDOffsetTPC) < nSigmaPIDMax && abs(nSigmaPiNeg - nSigmaPIDOffsetTPC) < nSigmaPIDMax) {
      sign = 1.;
    }
  }

  const auto pT = v0.pt();
  const auto eta = v0.eta();
  const auto tranRad = v0.v0radius();
  const auto dcaDaughv0 = v0.dcaV0daughters();
  const auto cpav0 = v0.v0cosPA();
  const std::vector<float> decVtx = {v0.x(), v0.y(), v0.z()};

  float observable = 0.;
  for (auto& sel : mSelections) {
    const auto selVariable = sel.getSelectionVariable();
    if (selVariable == femtoUniverseV0Selection::kV0DecVtxMax) {
      for (size_t i = 0; i < decVtx.size(); ++i) {
        auto decVtxValue = decVtx.at(i);
        sel.checkSelectionSetBit(decVtxValue, output, counter);
      }
    } else {
      switch (selVariable) {
        case (femtoUniverseV0Selection::kV0Sign):
          observable = sign;
          break;
        case (femtoUniverseV0Selection::kV0pTMin):
          observable = pT;
          break;
        case (femtoUniverseV0Selection::kV0pTMax):
          observable = pT;
          break;
        case (femtoUniverseV0Selection::kV0etaMax):
          observable = eta;
          break;
        case (femtoUniverseV0Selection::kV0DCADaughMax):
          observable = dcaDaughv0;
          break;
        case (femtoUniverseV0Selection::kV0CPAMin):
          observable = cpav0;
          break;
        case (femtoUniverseV0Selection::kV0TranRadMin):
          observable = tranRad;
          break;
        case (femtoUniverseV0Selection::kV0TranRadMax):
          observable = tranRad;
          break;
        case (femtoUniverseV0Selection::kV0DecVtxMax):
          break;
      }
      sel.checkSelectionSetBit(observable, output, counter);
    }
  }
  return {
    output,
    outputPosTrack.at(femtoUniverseTrackSelection::TrackContainerPosition::kCuts),
    outputPosTrack.at(femtoUniverseTrackSelection::TrackContainerPosition::kPID),
    outputNegTrack.at(femtoUniverseTrackSelection::TrackContainerPosition::kCuts),
    outputNegTrack.at(femtoUniverseTrackSelection::TrackContainerPosition::kPID)};
}

template <o2::aod::femtouniverseparticle::ParticleType part,
          o2::aod::femtouniverseparticle::ParticleType daugh, typename C,
          typename V, typename T>
void FemtoUniverseV0Selection::fillQA(C const& col, V const& v0, T const& posTrack,
                                      T const& negTrack)
{
  if (mHistogramRegistry) {
    mHistogramRegistry->fill(
      HIST(o2::aod::femtouniverseparticle::ParticleTypeName[part]) +
        HIST("/hPt"),
      v0.pt());
    mHistogramRegistry->fill(
      HIST(o2::aod::femtouniverseparticle::ParticleTypeName[part]) +
        HIST("/hEta"),
      v0.eta());
    mHistogramRegistry->fill(
      HIST(o2::aod::femtouniverseparticle::ParticleTypeName[part]) +
        HIST("/hPhi"),
      v0.phi());
    mHistogramRegistry->fill(
      HIST(o2::aod::femtouniverseparticle::ParticleTypeName[part]) +
        HIST("/hDaughDCA"),
      v0.dcaV0daughters());
    mHistogramRegistry->fill(
      HIST(o2::aod::femtouniverseparticle::ParticleTypeName[part]) +
        HIST("/hTransRadius"),
      v0.v0radius());
    mHistogramRegistry->fill(
      HIST(o2::aod::femtouniverseparticle::ParticleTypeName[part]) +
        HIST("/hDecayVtxX"),
      v0.x());
    mHistogramRegistry->fill(
      HIST(o2::aod::femtouniverseparticle::ParticleTypeName[part]) +
        HIST("/hDecayVtxY"),
      v0.y());
    mHistogramRegistry->fill(
      HIST(o2::aod::femtouniverseparticle::ParticleTypeName[part]) +
        HIST("/hDecayVtxZ"),
      v0.z());
    mHistogramRegistry->fill(
      HIST(o2::aod::femtouniverseparticle::ParticleTypeName[part]) +
        HIST("/hCPA"),
      v0.v0cosPA());
    mHistogramRegistry->fill(
      HIST(o2::aod::femtouniverseparticle::ParticleTypeName[part]) +
        HIST("/hCPAvsPt"),
      v0.pt(), v0.v0cosPA());
    mHistogramRegistry->fill(
      HIST(o2::aod::femtouniverseparticle::ParticleTypeName[part]) +
        HIST("/hInvMassLambda"),
      v0.mLambda());
    mHistogramRegistry->fill(
      HIST(o2::aod::femtouniverseparticle::ParticleTypeName[part]) +
        HIST("/hInvMassAntiLambda"),
      v0.mAntiLambda());
    mHistogramRegistry->fill(
      HIST(o2::aod::femtouniverseparticle::ParticleTypeName[part]) +
        HIST("/hInvMassLambdaAntiLambda"),
      v0.mLambda(), v0.mAntiLambda());
  }

  PosDaughTrack.fillQA<aod::femtouniverseparticle::ParticleType::kV0Child,
                       aod::femtouniverseparticle::TrackType::kPosChild>(posTrack);
  NegDaughTrack.fillQA<aod::femtouniverseparticle::ParticleType::kV0Child,
                       aod::femtouniverseparticle::TrackType::kNegChild>(negTrack);
}

} // namespace o2::analysis::femtoUniverse

#endif // PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSEV0SELECTION_H_
