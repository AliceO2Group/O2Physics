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

/// \file FemtoUniversePhiSelection.h
/// \brief Definition of the FemtoUniversePhiSelection
/// \author Valentina Mantovani Sarti, TU München valentina.mantovani-sarti@tum.de
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de
/// \author Luca Barioglio, TU München, luca.barioglio@cern.ch
/// \author Zuzanna Chochulska, WUT Warsaw, zuzanna.chochulska.stud@pw.edu.pl

#ifndef PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSEPHISELECTION_H_
#define PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSEPHISELECTION_H_

#include <iostream>
#include <string>
#include <vector>

#include "PWGCF/FemtoUniverse/Core/FemtoUniverseObjectSelection.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseSelection.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseTrackSelection.h"

#include "Common/Core/RecoDecay.h"
#include "Framework/HistogramRegistry.h"
#include "ReconstructionDataFormats/PID.h"
#include "TLorentzVector.h"

using namespace o2::framework;

namespace o2::analysis::femtoUniverse
{
namespace femtoUniversePhiSelection
{
/// The different selections this task is capable of doing
enum PhiSel {
  kPhiSign, ///< +1 particle, -1 antiparticle
  kPhipTMin,
  kPhipTMax,
  kPhietaMax,
  kPhiDCADaughMax,
  kPhiCPAMin,
  kPhiTranRadMin,
  kPhiTranRadMax,
  kPhiDecVtxMax
};

enum ChildTrackType { kPosTrack,
                      kNegTrack };

enum PhiContainerPosition {
  kPhi,
  kPosCuts,
  kPosPID,
  kNegCuts,
  kNegPID,
}; /// Position in the full VO cut container

} // namespace femtoUniversePhiSelection

/// \class FemtoUniversePhiSelection
/// \brief Cut class to contain and execute all cuts applied to Phis
class FemtoUniversePhiSelection
  : public FemtoUniverseObjectSelection<float, femtoUniversePhiSelection::PhiSel>
{
 public:
  FemtoUniversePhiSelection()
    : nPtPhiMinSel(0), nPtPhiMaxSel(0), nEtaPhiMaxSel(0), nDCAPhiDaughMax(0), nCPAPhiMin(0), nTranRadPhiMin(0), nTranRadPhiMax(0), nDecVtxMax(0), pTPhiMin(9999999.), pTPhiMax(-9999999.), etaPhiMax(-9999999.), DCAPhiDaughMax(-9999999.), CPAPhiMin(9999999.), TranRadPhiMin(9999999.), TranRadPhiMax(-9999999.), DecVtxMax(-9999999.), fInvMassLowLimit(1.05), fInvMassUpLimit(1.3), fRejectKaon(false), fInvMassKaonLowLimit(0.48), fInvMassKaonUpLimit(0.515), nSigmaPIDOffsetTPC(0.) {}
  /// Initializes histograms for the task
  template <o2::aod::femtouniverseparticle::ParticleType part,
            o2::aod::femtouniverseparticle::ParticleType daugh,
            typename cutContainerType>
  void init(HistogramRegistry* registry);

  template <typename C, typename V, typename T>
  bool isSelectedMinimal(C const& col, V const& phi, T const& posTrack,
                         T const& negTrack);

  template <typename C, typename V, typename T>
  void fillLambdaQA(C const& col, V const& phi, T const& posTrack,
                    T const& negTrack);

  /// \todo for the moment the PID of the tracks is factored out into a separate
  /// field, hence 5 values in total \\ASK: what does it mean?
  template <typename cutContainerType, typename C, typename V, typename T>
  std::array<cutContainerType, 5> getCutContainer(C const& col, V const& phi,
                                                  T const& posTrack,
                                                  T const& negTrack);

  template <o2::aod::femtouniverseparticle::ParticleType part,
            o2::aod::femtouniverseparticle::ParticleType daugh, typename C,
            typename V, typename T, typename Q>
  void fillQA(C const& col, V const& phi, T const& posTrack, T const& negTrack, Q const& posPID, Q const& negPID);

  template <typename T1, typename T2>
  void setChildCuts(femtoUniversePhiSelection::ChildTrackType child, T1 selVal,
                    T2 selVar, femtoUniverseSelection::SelectionType selType)
  {
    if (child == femtoUniversePhiSelection::kPosTrack) {
      PosDaughTrack.setSelection(selVal, selVar, selType);
    } else if (child == femtoUniversePhiSelection::kNegTrack) {
      NegDaughTrack.setSelection(selVal, selVar, selType);
    }
  }
  template <typename T>
  void setChildPIDSpecies(femtoUniversePhiSelection::ChildTrackType child,
                          T& pids)
  {
    if (child == femtoUniversePhiSelection::kPosTrack) {
      PosDaughTrack.setPIDSpecies(pids);
    } else if (child == femtoUniversePhiSelection::kNegTrack) {
      NegDaughTrack.setPIDSpecies(pids);
    }
  }

  /// Helper function to obtain the name of a given selection criterion for consistent naming of the configurables
  /// \param iSel Track selection variable to be examined
  /// \param prefix Additional prefix for the name of the configurable
  /// \param suffix Additional suffix for the name of the configurable
  static std::string getSelectionName(femtoUniversePhiSelection::PhiSel iSel,
                                      std::string_view prefix = "",
                                      std::string_view suffix = "")
  {
    std::string outString = static_cast<std::string>(prefix);
    outString += static_cast<std::string>(mSelectionNames[iSel]);
    outString += suffix;
    return outString;
  }

  /// Helper function to obtain the index of a given selection variable for consistent naming of the configurables
  /// \param obs Phi selection variable (together with prefix) got from file
  /// \param prefix Additional prefix for the output of the configurable
  static int findSelectionIndex(const std::string_view& obs,
                                std::string_view prefix = "")
  {
    for (int index = 0; index < kNphiSelection; index++) {
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
  /// \param iSel Phi selection variable whose type is returned
  static femtoUniverseSelection::SelectionType
    getSelectionType(femtoUniversePhiSelection::PhiSel iSel)
  {
    return mSelectionTypes[iSel];
  }

  /// Helper function to obtain the helper string of a given selection criterion
  /// for consistent description of the configurables
  /// \param iSel Track selection variable to be examined
  /// \param prefix Additional prefix for the output of the configurable
  static std::string getSelectionHelper(femtoUniversePhiSelection::PhiSel iSel,
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

  void setChildRejectNotPropagatedTracks(femtoUniversePhiSelection::ChildTrackType child, bool reject)
  {
    if (child == femtoUniversePhiSelection::kPosTrack) {
      PosDaughTrack.setRejectNotPropagatedTracks(reject);
    } else if (child == femtoUniversePhiSelection::kNegTrack) {
      NegDaughTrack.setRejectNotPropagatedTracks(reject);
    }
  }

  void setChildnSigmaPIDOffset(femtoUniversePhiSelection::ChildTrackType child, float offsetTPC, float offsetTOF)
  {
    if (child == femtoUniversePhiSelection::kPosTrack) {
      PosDaughTrack.setnSigmaPIDOffset(offsetTPC, offsetTOF);
    } else if (child == femtoUniversePhiSelection::kNegTrack) {
      NegDaughTrack.setnSigmaPIDOffset(offsetTPC, offsetTOF);
    }
  }

 private:
  int nPtPhiMinSel;
  int nPtPhiMaxSel;
  int nEtaPhiMaxSel;
  int nDCAPhiDaughMax;
  int nCPAPhiMin;
  int nTranRadPhiMin;
  int nTranRadPhiMax;
  int nDecVtxMax;
  float pTPhiMin;
  float pTPhiMax;
  float etaPhiMax;
  float DCAPhiDaughMax;
  float CPAPhiMin;
  float TranRadPhiMin;
  float TranRadPhiMax;
  float DecVtxMax;

  float fInvMassLowLimit;
  float fInvMassUpLimit;

  bool fRejectKaon;
  float fInvMassKaonLowLimit;
  float fInvMassKaonUpLimit;

  float nSigmaPIDOffsetTPC;

  FemtoUniverseTrackSelection PosDaughTrack;
  FemtoUniverseTrackSelection NegDaughTrack;

  static constexpr int kNphiSelection = 9;

  static constexpr std::string_view mSelectionNames[kNphiSelection] = {
    "Sign", "PtMin", "PtMax", "EtaMax", "DCAdaughMax", "CPAMin",
    "TranRadMin", "TranRadMax", "DecVecMax"}; ///< Name of the different
                                              ///< selections

  static constexpr femtoUniverseSelection::SelectionType
    mSelectionTypes[kNphiSelection]{
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

  static constexpr std::string_view mSelectionHelper[kNphiSelection] = {
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
void FemtoUniversePhiSelection::init(HistogramRegistry* registry)
{

  if (registry) {
    mHistogramRegistry = registry;
    fillSelectionHistogram<part>();
    fillSelectionHistogram<daugh>();

    AxisSpec massAxisPhi = {6000, 0.9f, 3.0f, "m_{#Phi} (GeV/#it{c}^{2})"};
    AxisSpec massAxisLambda = {600, 0.0f, 3.0f, "m_{#Lambda} (GeV/#it{c}^{2})"};
    AxisSpec massAxisAntiLambda = {600, 0.0f, 3.0f,
                                   "m_{#bar{#Lambda}} (GeV/#it{c}^{2})"};

    /// \todo this should be an automatic check in the parent class, and the
    /// return type should be templated
    size_t nSelections = getNSelections();
    if (nSelections > 8 * sizeof(cutContainerType)) {
      LOG(fatal) << "FemtoUniversePhiCuts: Number of selections to large for your "
                    "container - quitting!";
    }
    std::string folderName = static_cast<std::string>(
      o2::aod::femtouniverseparticle::ParticleTypeName[part]);
    /// \todo initialize histograms for children tracks of phis
    mHistogramRegistry->add((folderName + "/hPt").c_str(),
                            "; #it{p}_{T} (GeV/#it{c}); Entries", kTH1F,
                            {{1000, 0, 10}});
    mHistogramRegistry->add((folderName + "/hEta").c_str(), "; #eta; Entries",
                            kTH1F, {{1000, -1, 1}});
    mHistogramRegistry->add((folderName + "/hPhi").c_str(), "; #phi; Entries",
                            kTH1F, {{1000, 0, 2. * M_PI}});
    mHistogramRegistry->add((folderName + "/hInvMassPhi").c_str(), "", kTH1F,
                            {massAxisPhi});

    PosDaughTrack.init<aod::femtouniverseparticle::ParticleType::kPhiChild,
                       aod::femtouniverseparticle::TrackType::kPosChild,
                       aod::femtouniverseparticle::cutContainerType>(
      mHistogramRegistry);
    NegDaughTrack.init<aod::femtouniverseparticle::ParticleType::kPhiChild,
                       aod::femtouniverseparticle::TrackType::kNegChild,
                       aod::femtouniverseparticle::cutContainerType>(
      mHistogramRegistry);

    // mHistogramRegistry->add("LambdaQA/hInvMassLambdaNoCuts", "No cuts", kTH1F,
    //                         {massAxisLambda});
    // mHistogramRegistry->add("LambdaQA/hInvMassLambdaInvMassCut",
    //                         "Invariant mass cut", kTH1F, {massAxisLambda});
    // mHistogramRegistry->add("LambdaQA/hInvMassLambdaPtMin", "Minimum Pt cut",
    //                         kTH1F, {massAxisLambda});
    // mHistogramRegistry->add("LambdaQA/hInvMassLambdaPtMax", "Maximum Pt cut",
    //                         kTH1F, {massAxisLambda});
    // mHistogramRegistry->add("LambdaQA/hInvMassLambdaEtaMax", "Maximum Eta cut",
    //                         kTH1F, {massAxisLambda});
    // mHistogramRegistry->add("LambdaQA/hInvMassLambdaDCAPhiDaugh",
    //                         "Phi-daughters DCA cut", kTH1F, {massAxisLambda});
    // mHistogramRegistry->add("LambdaQA/hInvMassLambdaCPA", "CPA cut", kTH1F,
    //                         {massAxisLambda});
    // mHistogramRegistry->add("LambdaQA/hInvMassLambdaTranRadMin",
    //                         "Minimum transverse radius cut", kTH1F,
    //                         {massAxisLambda});
    // mHistogramRegistry->add("LambdaQA/hInvMassLambdaTranRadMax",
    //                         "Maximum transverse radius cut", kTH1F,
    //                         {massAxisLambda});
    // mHistogramRegistry->add("LambdaQA/hInvMassLambdaDecVtxMax",
    //                         "Maximum distance on  decay vertex cut", kTH1F,
    //                         {massAxisLambda});
  }
  /// check whether the most open cuts are fulfilled - most of this should have
  /// already be done by the filters
  nPtPhiMinSel = getNSelections(femtoUniversePhiSelection::kPhipTMin);
  nPtPhiMaxSel = getNSelections(femtoUniversePhiSelection::kPhipTMax);
  nEtaPhiMaxSel = getNSelections(femtoUniversePhiSelection::kPhietaMax);
  nDCAPhiDaughMax = getNSelections(femtoUniversePhiSelection::kPhiDCADaughMax);
  nCPAPhiMin = getNSelections(femtoUniversePhiSelection::kPhiCPAMin);
  nTranRadPhiMin = getNSelections(femtoUniversePhiSelection::kPhiTranRadMin);
  nTranRadPhiMax = getNSelections(femtoUniversePhiSelection::kPhiTranRadMax);
  nDecVtxMax = getNSelections(femtoUniversePhiSelection::kPhiDecVtxMax);

  pTPhiMin = getMinimalSelection(femtoUniversePhiSelection::kPhipTMin,
                                 femtoUniverseSelection::kLowerLimit);
  pTPhiMax = getMinimalSelection(femtoUniversePhiSelection::kPhipTMax,
                                 femtoUniverseSelection::kUpperLimit);
  etaPhiMax = getMinimalSelection(femtoUniversePhiSelection::kPhietaMax,
                                  femtoUniverseSelection::kAbsUpperLimit);
  DCAPhiDaughMax = getMinimalSelection(femtoUniversePhiSelection::kPhiDCADaughMax,
                                       femtoUniverseSelection::kUpperLimit);
  CPAPhiMin = getMinimalSelection(femtoUniversePhiSelection::kPhiCPAMin,
                                  femtoUniverseSelection::kLowerLimit);
  TranRadPhiMin = getMinimalSelection(femtoUniversePhiSelection::kPhiTranRadMin,
                                      femtoUniverseSelection::kLowerLimit);
  TranRadPhiMax = getMinimalSelection(femtoUniversePhiSelection::kPhiTranRadMax,
                                      femtoUniverseSelection::kUpperLimit);
  DecVtxMax = getMinimalSelection(femtoUniversePhiSelection::kPhiDecVtxMax,
                                  femtoUniverseSelection::kAbsUpperLimit);
}

template <typename C, typename V, typename T>
bool FemtoUniversePhiSelection::isSelectedMinimal(C const& col, V const& phi,
                                                  T const& posTrack,
                                                  T const& negTrack)
{
  const auto signPos = posTrack.sign();
  const auto signNeg = negTrack.sign();
  if (signPos < 0 || signNeg > 0) {
    LOG(warn) << "Something wrong in isSelectedMinimal";
    LOG(warn) << "ERROR - Wrong sign for Phi daughters";
  }
  // asfaf
  const float pT = phi.pt();
  const float eta = phi.eta();
  const std::vector<float> decVtx = {phi.x(), phi.y(), phi.z()};
  const float tranRad = phi.phiradius();
  const float dcaDaughphi = phi.dcaPhidaughters();
  const float cpaphi = phi.phicosPA(col.posX(), col.posY(), col.posZ());

  const float invMassLambda = phi.mLambda();
  const float invMassAntiLambda = phi.mAntiLambda();

  if ((invMassLambda < fInvMassLowLimit || invMassLambda > fInvMassUpLimit) &&
      (invMassAntiLambda < fInvMassLowLimit ||
       invMassAntiLambda > fInvMassUpLimit)) {
    return false;
  }
  if (fRejectKaon) {
    const float invMassKaon = phi.mK0Short();
    if (invMassKaon > fInvMassKaonLowLimit &&
        invMassKaon < fInvMassKaonUpLimit) {
      return false;
    }
  }
  if (nPtPhiMinSel > 0 && pT < pTPhiMin) {
    return false;
  }
  if (nPtPhiMaxSel > 0 && pT > pTPhiMax) {
    return false;
  }
  if (nEtaPhiMaxSel > 0 && std::abs(eta) > etaPhiMax) {
    return false;
  }
  if (nDCAPhiDaughMax > 0 && dcaDaughphi > DCAPhiDaughMax) {
    return false;
  }
  if (nCPAPhiMin > 0 && cpaphi < CPAPhiMin) {
    return false;
  }
  if (nTranRadPhiMin > 0 && tranRad < TranRadPhiMin) {
    return false;
  }
  if (nTranRadPhiMax > 0 && tranRad > TranRadPhiMax) {
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

  // check that track combinations for Phi or antiPhi would be fulfilling PID
  int nSigmaPIDMax = PosDaughTrack.getSigmaPIDMax();
  // antiPhi
  auto nSigmaPrNeg = negTrack.tpcNSigmaPr();
  auto nSigmaPiPos = posTrack.tpcNSigmaPi();
  // phi
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
void FemtoUniversePhiSelection::fillLambdaQA(C const& col, V const& phi,
                                             T const& posTrack, T const& negTrack)
{
  const auto signPos = posTrack.sign();
  const auto signNeg = negTrack.sign();
  if (signPos < 0 || signNeg > 0) {
    LOG(warn) << "Something wrong in isSelectedMinimal";
    LOG(warn) << "ERROR - Wrong sign for Phi daughters";
  }
  const float pT = phi.pt();
  const float eta = phi.eta();
  const std::vector<float> decVtx = {phi.x(), phi.y(), phi.z()};
  const float tranRad = phi.phiradius();
  const float dcaDaughphi = phi.dcaPhidaughters();
  const float cpaphi = phi.phicosPA(col.posX(), col.posY(), col.posZ());

  const float invMassLambda = phi.mLambda();

  mHistogramRegistry->fill(HIST("LambdaQA/hInvMassLambdaNoCuts"), phi.mLambda());

  if (invMassLambda > fInvMassLowLimit && invMassLambda < fInvMassUpLimit) {
    mHistogramRegistry->fill(HIST("LambdaQA/hInvMassLambdaInvMassCut"),
                             phi.mLambda());
  }

  if (pT > pTPhiMin) {
    mHistogramRegistry->fill(HIST("LambdaQA/hInvMassLambdaPtMin"),
                             phi.mLambda());
  }
  if (pT < pTPhiMax) {
    mHistogramRegistry->fill(HIST("LambdaQA/hInvMassLambdaPtMax"),
                             phi.mLambda());
  }
  if (std::abs(eta) < etaPhiMax) {
    mHistogramRegistry->fill(HIST("LambdaQA/hInvMassLambdaEtaMax"),
                             phi.mLambda());
  }
  if (dcaDaughphi < DCAPhiDaughMax) {
    mHistogramRegistry->fill(HIST("LambdaQA/hInvMassLambdaDCAPhiDaugh"),
                             phi.mLambda());
  }
  if (cpaphi > CPAPhiMin) {
    mHistogramRegistry->fill(HIST("LambdaQA/hInvMassLambdaCPA"), phi.mLambda());
  }
  if (tranRad > TranRadPhiMin) {
    mHistogramRegistry->fill(HIST("LambdaQA/hInvMassLambdaTranRadMin"),
                             phi.mLambda());
  }
  if (tranRad < TranRadPhiMax) {
    mHistogramRegistry->fill(HIST("LambdaQA/hInvMassLambdaTranRadMax"),
                             phi.mLambda());
  }
  bool write = true;
  for (size_t i = 0; i < decVtx.size(); i++) {
    write = write && (decVtx.at(i) < DecVtxMax);
  }
  if (write) {
    mHistogramRegistry->fill(HIST("LambdaQA/hInvMassLambdaDecVtxMax"),
                             phi.mLambda());
  }
}

/// the CosPA of Phi needs as argument the posXYZ of collisions vertex so we need
/// to pass the collsion as well
template <typename cutContainerType, typename C, typename V, typename T>
std::array<cutContainerType, 5>
  FemtoUniversePhiSelection::getCutContainer(C const& col, V const& phi, T const& posTrack, T const& negTrack)
{
  auto outputPosTrack = PosDaughTrack.getCutContainer<cutContainerType>(posTrack);
  auto outputNegTrack = NegDaughTrack.getCutContainer<cutContainerType>(negTrack);
  cutContainerType output = 0;
  size_t counter = 0;

  auto lambdaMassNominal = TDatabasePDG::Instance()->GetParticle(3122)->Mass();
  auto lambdaMassHypothesis = phi.mLambda();
  auto antiLambdaMassHypothesis = phi.mAntiLambda();
  auto diffLambda = abs(lambdaMassNominal - lambdaMassHypothesis);
  auto diffAntiLambda = abs(antiLambdaMassHypothesis - lambdaMassHypothesis);

  float sign = 0.;
  int nSigmaPIDMax = PosDaughTrack.getSigmaPIDMax();
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

  const auto pT = phi.pt();
  const auto eta = phi.eta();
  const auto tranRad = phi.phiradius();
  const auto dcaDaughphi = phi.dcaPhidaughters();
  const auto cpaphi = phi.phicosPA(col.posX(), col.posY(), col.posZ());
  const std::vector<float> decVtx = {phi.x(), phi.y(), phi.z()};

  float observable = 0.;
  for (auto& sel : mSelections) {
    const auto selVariable = sel.getSelectionVariable();
    if (selVariable == femtoUniversePhiSelection::kPhiDecVtxMax) {
      for (size_t i = 0; i < decVtx.size(); ++i) {
        auto decVtxValue = decVtx.at(i);
        sel.checkSelectionSetBit(decVtxValue, output, counter);
      }
    } else {
      switch (selVariable) {
        case (femtoUniversePhiSelection::kPhiSign):
          observable = sign;
          break;
        case (femtoUniversePhiSelection::kPhipTMin):
          observable = pT;
          break;
        case (femtoUniversePhiSelection::kPhipTMax):
          observable = pT;
          break;
        case (femtoUniversePhiSelection::kPhietaMax):
          observable = eta;
          break;
        case (femtoUniversePhiSelection::kPhiDCADaughMax):
          observable = dcaDaughphi;
          break;
        case (femtoUniversePhiSelection::kPhiCPAMin):
          observable = cpaphi;
          break;
        case (femtoUniversePhiSelection::kPhiTranRadMin):
          observable = tranRad;
          break;
        case (femtoUniversePhiSelection::kPhiTranRadMax):
          observable = tranRad;
          break;
        case (femtoUniversePhiSelection::kPhiDecVtxMax):
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
          typename V, typename T, typename Q>
void FemtoUniversePhiSelection::fillQA(C const& col, V const& phi, T const& posTrack,
                                       T const& negTrack, Q const& posPID, Q const& negPID)
{
  if (mHistogramRegistry) {
    TLorentzVector part1Vec;
    TLorentzVector part2Vec;
    float mMassOne = TDatabasePDG::Instance()->GetParticle(321)->Mass();
    float mMassTwo = TDatabasePDG::Instance()->GetParticle(321)->Mass();

    part1Vec.SetPtEtaPhiM(posTrack.pt(), posTrack.eta(), posTrack.phi(), mMassOne);
    part2Vec.SetPtEtaPhiM(negTrack.pt(), negTrack.eta(), negTrack.phi(), mMassTwo);

    TLorentzVector sumVec(part1Vec);
    sumVec += part2Vec;
    float phiEta = sumVec.Eta();
    float phiPt = sumVec.Pt();
    // float phiP = sumVec.P();
    float phiPhi = sumVec.Phi();
    if (sumVec.Phi() < 0) {
      phiPhi = sumVec.Phi() + 2 * o2::constants::math::PI;
    } else if (sumVec.Phi() >= 0) {
      phiPhi = sumVec.Phi();
    }

    float phiM = sumVec.M();

    mHistogramRegistry->fill(
      HIST(o2::aod::femtouniverseparticle::ParticleTypeName[part]) +
        HIST("/hPt"),
      phiPt);
    mHistogramRegistry->fill(
      HIST(o2::aod::femtouniverseparticle::ParticleTypeName[part]) +
        HIST("/hEta"),
      phiEta);
    mHistogramRegistry->fill(
      HIST(o2::aod::femtouniverseparticle::ParticleTypeName[part]) +
        HIST("/hPhi"),
      phiPhi);

    mHistogramRegistry->fill(
      HIST(o2::aod::femtouniverseparticle::ParticleTypeName[part]) +
        HIST("/hInvMassPhi"),
      phiM);
  }

  PosDaughTrack.fillQA<aod::femtouniverseparticle::ParticleType::kPhiChild,
                       aod::femtouniverseparticle::TrackType::kPosChild>(posTrack);
  NegDaughTrack.fillQA<aod::femtouniverseparticle::ParticleType::kPhiChild,
                       aod::femtouniverseparticle::TrackType::kNegChild>(negTrack);
}

} // namespace o2::analysis::femtoUniverse

#endif // PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSEPHISELECTION_H_
