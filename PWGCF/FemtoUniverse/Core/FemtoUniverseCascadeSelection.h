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

/// \file femtoUniverseCascadeSelection.h
/// \brief Definition of the femtoUniverseCascadeSelection
/// \author Valentina Mantovani Sarti, TU München valentina.mantovani-sarti@tum.de
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de
/// \author Luca Barioglio, TU München, luca.barioglio@cern.ch
/// \author Zuzanna Chochulska, WUT Warsaw & CTU Prague, zchochul@cern.ch
/// \author Barbara Chytla, WUT Warsaw, barbara.chytla@cern.ch

#ifndef PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSECASCADESELECTION_H_
#define PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSECASCADESELECTION_H_

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
namespace femtoUniverseCascadeSelection
{
/// The different selections this task is capable of doing
enum CascadeSel {
  kCascadeSign, ///< +1 particle, -1 antiparticle
  kCascadepTMin,
  kCascadepTMax,
  kCascadeetaMax,
  kCascadeV0DCADaughMax,
  kCascadeV0CPAMin,
  kCascadeV0TranRadMin,
  kCascadeV0TranRadMax,
  kCascadeV0DecVtxMax,
  kCascadeDCADaughMax,
  kCascadeCPAMin,
  kCascadeTranRadMin,
  kCascadeTranRadMax,
  kCascadeDecVtxMax,
  kCascadeDCAPosToPV,
  kCascadeDCANegToPV,
  kCascadeDCABachToPV,
  kCascadeDCAV0ToPV,
  kCascadeV0MassMin,
  kCascadeV0MassMax
};

enum ChildTrackType { kPosTrack,
                      kNegTrack,
                      kBachTrack };

/*enum CascadeContainerPosition {
  kCascade,
  kPosCuts,
  kPosPID,
  kNegCuts,
  kNegPID,
}; /// Position in the full VO cut container (for cutculator)
*/
} // namespace femtoUniverseCascadeSelection

/// \class FemtoUniverseCascadeSelection
/// \brief Cut class to contain and execute all cuts applied to Cascades
class FemtoUniverseCascadeSelection
  : public FemtoUniverseObjectSelection<float, femtoUniverseCascadeSelection::CascadeSel>
{
  // do cutculatora

 public:
  FemtoUniverseCascadeSelection()
    : nPtCascadeMinSel(0), nPtCascadeMaxSel(0), nEtaCascadeMaxSel(0), nDCAV0DaughMax(0), nCPAV0Min(0), nTranRadV0Min(0), nTranRadV0Max(0), nV0DecVtxMax(0), nDCACascadeDaughMax(0), nCPACascadeMin(0), nTranRadCascadeMin(0), nTranRadCascadeMax(0), nDecVtxMax(0), nDCAPosToPV(0), nDCANegToPV(0), nDCABachToPV(0), nDCAV0ToPV(0), pTCascadeMin(9999999.), pTCascadeMax(-9999999.), etaCascadeMax(-9999999.), DCAV0DaughMax(-9999999.), CPAV0Min(9999999.), TranRadV0Min(9999999.), TranRadV0Max(-9999999.), V0DecVtxMax(-9999999.), DCACascadeDaughMax(-9999999.), CPACascadeMin(9999999.), TranRadCascadeMin(9999999.), TranRadCascadeMax(-9999999.), DecVtxMax(-9999999.), DCAPosToPV(9999999.), DCANegToPV(9999999.), DCABachToPV(9999999.), DCAV0ToPV(9999999.), fV0InvMassLowLimit(1.05), fV0InvMassUpLimit(1.3), fInvMassLowLimit(1.25), fInvMassUpLimit(1.4), fRejectOmega(false), fInvMassOmegaLowLimit(1.5), fInvMassOmegaUpLimit(2.0) /*, nSigmaPIDOffsetTPC(0.)*/
  {
  }

  /// Initializes histograms for the task
  template <o2::aod::femtouniverseparticle::ParticleType part, o2::aod::femtouniverseparticle::ParticleType daugh, o2::aod::femtouniverseparticle::ParticleType bach, typename cutContainerType>
  void init(HistogramRegistry* registry);

  template <typename Col, typename Casc, typename Track>
  bool isSelectedMinimal(Col const& col, Casc const& cascade, Track const& posTrack, Track const& negTrack, Track const& bachTrack);

  template <typename Col, typename Casc, typename Track>
  void fillCascadeQA(Col const& col, Casc const& cascade, Track const& posTrack, Track const& negTrack);

  template <typename Col, typename Casc, typename Track>
  void fillQA(Col const& col, Casc const& cascade, Track const& posTrack, Track const& negTrack, Track const& bachTrack);

  template <typename T1, typename T2>
  void setChildCuts(femtoUniverseCascadeSelection::ChildTrackType child, T1 selVal,
                    T2 selVar, femtoUniverseSelection::SelectionType selType)
  {
    if (child == femtoUniverseCascadeSelection::kPosTrack) {
      PosDaughTrack.setSelection(selVal, selVar, selType);
    } else if (child == femtoUniverseCascadeSelection::kNegTrack) {
      NegDaughTrack.setSelection(selVal, selVar, selType);
    } else if (child == femtoUniverseCascadeSelection::kBachTrack) {
      BachTrack.setSelection(selVal, selVar, selType);
    }
  }

  template <typename T>
  void setChildPIDSpecies(femtoUniverseCascadeSelection::ChildTrackType child,
                          T& pids)
  {
    if (child == femtoUniverseCascadeSelection::kPosTrack) {
      PosDaughTrack.setPIDSpecies(pids);
    } else if (child == femtoUniverseCascadeSelection::kNegTrack) {
      NegDaughTrack.setPIDSpecies(pids);
    } else if (child == femtoUniverseCascadeSelection::kBachTrack) {
      BachTrack.setPIDSpecies(pids);
    }
  }

  /// Helper function to obtain the name of a given selection criterion for consistent naming of the configurables
  /// \param iSel Track selection variable to be examined
  /// \param prefix Additional prefix for the name of the configurable
  /// \param suffix Additional suffix for the name of the configurable
  static std::string getSelectionName(femtoUniverseCascadeSelection::CascadeSel iSel,
                                      std::string_view prefix = "",
                                      std::string_view suffix = "")
  {
    std::string outString = static_cast<std::string>(prefix);
    outString += static_cast<std::string>(mSelectionNames[iSel]);
    outString += suffix;
    return outString;
  }

  /// Helper function to obtain the helper string of a given selection criterion
  /// for consistent description of the configurables
  /// \param iSel Track selection variable to be examined
  /// \param prefix Additional prefix for the output of the configurable
  static std::string getSelectionHelper(femtoUniverseCascadeSelection::CascadeSel iSel,
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

  /// Set limit for the omega rejection on the invariant mass
  /// \param lowLimit Lower limit for the invariant mass distribution
  /// \param upLimit Upper limit for the invariant mass distribution
  void setOmegaInvMassLimits(float lowLimit, float upLimit)
  {
    fRejectOmega = true;
    fInvMassOmegaLowLimit = lowLimit;
    fInvMassOmegaUpLimit = upLimit;
  }

 private:
  int nPtCascadeMinSel;
  int nPtCascadeMaxSel;
  int nEtaCascadeMaxSel;
  int nDCAV0DaughMax;
  int nCPAV0Min;
  int nTranRadV0Min;
  int nTranRadV0Max;
  int nV0DecVtxMax;
  int nDCACascadeDaughMax;
  int nCPACascadeMin;
  int nTranRadCascadeMin;
  int nTranRadCascadeMax;
  int nDecVtxMax;
  int nDCAPosToPV;
  int nDCANegToPV;
  int nDCABachToPV;
  int nDCAV0ToPV;
  float pTCascadeMin;
  float pTCascadeMax;
  float etaCascadeMax;
  float DCAV0DaughMax;
  float CPAV0Min;
  float TranRadV0Min;
  float TranRadV0Max;
  float V0DecVtxMax;
  float DCACascadeDaughMax;
  float CPACascadeMin;
  float TranRadCascadeMin;
  float TranRadCascadeMax;
  float DecVtxMax;
  float DCAPosToPV;
  float DCANegToPV;
  float DCABachToPV;
  float DCAV0ToPV;

  float fV0InvMassLowLimit;
  float fV0InvMassUpLimit;

  float fInvMassLowLimit;
  float fInvMassUpLimit;

  float fRejectOmega;
  float fInvMassOmegaLowLimit;
  float fInvMassOmegaUpLimit;

  // float nSigmaPIDOffsetTPC;

  FemtoUniverseTrackSelection PosDaughTrack;
  FemtoUniverseTrackSelection NegDaughTrack;
  FemtoUniverseTrackSelection BachTrack;

  static constexpr int kNcascadeSelection = 20; // can I do less ?

  static constexpr std::string_view mSelectionNames[kNcascadeSelection] = {
    "Sign", "PtMin", "PtMax", "EtaMax", "DCAv0daughMax", "v0CPAMin",
    "v0TranRadMin", "v0TranRadMax", "v0DecVecMax", "DCAcascDaugh",
    "CPAMin", "TranRadMin", "TranRadMax", "DecVtxMax",
    "DCAPosToPV", "DCANegToPV", "DCABachToPV", "DCAV0ToPV",
    "kV0MassMin", "V0MassMax"}; ///< Name of the different
                                ///< selections

  static constexpr femtoUniverseSelection::SelectionType
    mSelectionTypes[kNcascadeSelection]{
      femtoUniverseSelection::kEqual,      // sign
      femtoUniverseSelection::kLowerLimit, // pt min
      femtoUniverseSelection::kUpperLimit, // pt max
      femtoUniverseSelection::kUpperLimit, // eta max
      femtoUniverseSelection::kUpperLimit, // DCA v0 daughters max
      femtoUniverseSelection::kLowerLimit, // v0 cos PA min
      femtoUniverseSelection::kLowerLimit, // v0 tran rad min
      femtoUniverseSelection::kUpperLimit, // v0 tran rad max
      femtoUniverseSelection::kUpperLimit, // v0 maximum distance of decay vertex to PV
      femtoUniverseSelection::kUpperLimit, // DCA cascade daughters max
      femtoUniverseSelection::kLowerLimit, // cascade cos PA min
      femtoUniverseSelection::kLowerLimit, // cascade tran rad min
      femtoUniverseSelection::kUpperLimit, // cascade tran rad max
      femtoUniverseSelection::kUpperLimit, // cascade maximum distance of decay vertex to PV
      femtoUniverseSelection::kLowerLimit, // DCA pos to PV max
      femtoUniverseSelection::kLowerLimit, // DCA neg to PV max
      femtoUniverseSelection::kLowerLimit, // DCA bach to PV max
      femtoUniverseSelection::kLowerLimit, // DCA v0 to PV max
      femtoUniverseSelection::kLowerLimit, // v0 mass min
      femtoUniverseSelection::kUpperLimit, // v0 mass max
    };                                     ///< Map to match a variable with
                                           ///< its type

  static constexpr std::string_view mSelectionHelper[kNcascadeSelection] = {
    "Cascade particle sign (+1 or -1)",
    "Minimum pT (GeV/c)",
    "Maximum pT (GeV/c)",
    "Maximum |Eta|",
    "Maximum DCA between v0 daughters (cm)",
    "Minimum Cosine of Pointing Angle for v0",
    "Minimum v0 transverse radius (cm)",
    "Maximum v0 transverse radius (cm)",
    "Maximum distance of v0 from primary vertex",
    "Maximum DCA between cascade daughters (cm)",
    "Minimum Cosine of Pointing Angle for cascade",
    "Minimum cascade transverse radius (cm)",
    "Maximum cascade transverse radius (cm)",
    "Maximum distance of cascade from primary vertex",
    "Maximum DCA of positive track form primary vertex",
    "Maximum DCA of negative track form primary vertex",
    "Maximum DCA of bachelor track form primary vertex",
    "Maximum DCA of v0 form primary vertex",
    "Minimum V0 mass",
    "Maximum V0 mass"}; ///< Helper information for the
                        ///< different selections

}; // namespace femtoUniverse

template <o2::aod::femtouniverseparticle::ParticleType part, o2::aod::femtouniverseparticle::ParticleType daugh, o2::aod::femtouniverseparticle::ParticleType bach, typename cutContainerType>
void FemtoUniverseCascadeSelection::init(HistogramRegistry* registry)
{

  if (registry) {
    mHistogramRegistry = registry;
    fillSelectionHistogram<part>();  // cascade
    fillSelectionHistogram<daugh>(); // pos, neg
    fillSelectionHistogram<bach>();  // bach

    AxisSpec massAxisCascade = {600, 1.25f, 1.4f, "m_{#Cascade} (GeV/#it{c}^{2})"};
    AxisSpec massAxisV0 = {600, 0.0f, 3.0f, "m_{#V0} (GeV/#it{c}^{2})"};
    AxisSpec DCADaughAxis = {1000, 0.0f, 2.0f, "DCA (cm)"};
    AxisSpec DCAToPVAxis = {1000, -10.0f, 10.0f, "DCA to PV (cm)"};
    AxisSpec ptAxis = {100, 0.0f, 10.0f, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec etaAxis = {100, -2.0f, 2.0f, "#it{#eta}"};
    AxisSpec phiAxis = {100, 0.0f, 6.0f, "#it{#phi}"};
    AxisSpec CPAAxis = {1000, 0.95f, 1.0f, "#it{cos #theta_{p}}"};
    AxisSpec tranRadAxis = {1000, 0.0f, 100.0f, "#it{r}_{xy} (cm)"};

    /// \todo this should be an automatic check in the parent class, and the
    /// return type should be templated
    size_t nSelections = getNSelections();
    if (nSelections > 17 * sizeof(cutContainerType)) {
      LOG(fatal) << "FemtoUniverseCascadeCuts: Number of selections to large for your "
                    "container - quitting!";
    }

    PosDaughTrack.init<aod::femtouniverseparticle::ParticleType::kV0Child,
                       aod::femtouniverseparticle::TrackType::kPosChild,
                       aod::femtouniverseparticle::cutContainerType>(
      mHistogramRegistry);
    NegDaughTrack.init<aod::femtouniverseparticle::ParticleType::kV0Child,
                       aod::femtouniverseparticle::TrackType::kNegChild,
                       aod::femtouniverseparticle::cutContainerType>(
      mHistogramRegistry);
    BachTrack.init<aod::femtouniverseparticle::ParticleType::kCascadeBachelor,
                   aod::femtouniverseparticle::TrackType::kBachelor,
                   aod::femtouniverseparticle::cutContainerType>(
      mHistogramRegistry);

    // V0 (Lambda)
    // mHistogramRegistry->add("CascadeQA/hInvMassV0NoCuts", "No cuts", kTH1F, {massAxisV0});
    mHistogramRegistry->add("CascadeQA/hInvMassV0Cut", "Invariant mass cut", kTH1F, {massAxisV0});
    mHistogramRegistry->add("CascadeQA/hDCAV0Daugh", "V0-daughters DCA", kTH1F, {DCADaughAxis});
    mHistogramRegistry->add("CascadeQA/hV0CPA", "V0 cos PA", kTH1F, {CPAAxis});
    mHistogramRegistry->add("CascadeQA/hV0TranRad", "V0 transverse radius", kTH1F, {tranRadAxis});
    // mHistogramRegistry->add("CascadeQA/hV0DecVtxMax", "V0 maximum distance on decay vertex", kTH1F, {massAxisV0});

    // Cascade (Xi, Omega)
    // mHistogramRegistry->add("CascadeQA/hInvMassCascadeNoCuts", "No cuts", kTH1F, {massAxisCascade});
    mHistogramRegistry->add("CascadeQA/hInvMassCascadeCut", "Invariant mass with cut", kTH1F, {massAxisCascade});
    mHistogramRegistry->add("CascadeQA/hCascadePt", "pT distribution", kTH1F, {ptAxis});
    mHistogramRegistry->add("CascadeQA/hCascadeEta", "Eta distribution", kTH1F, {etaAxis});
    mHistogramRegistry->add("CascadeQA/hCascadePhi", "Phi distribution", kTH1F, {phiAxis});
    mHistogramRegistry->add("CascadeQA/hDCACascadeDaugh", "Cascade-daughters DCA", kTH1F, {DCADaughAxis});
    mHistogramRegistry->add("CascadeQA/hCascadeCPA", "Cos PA", kTH1F, {CPAAxis});
    mHistogramRegistry->add("CascadeQA/hCascadeTranRad", "Transverse radius", kTH1F, {tranRadAxis});
    mHistogramRegistry->add("CascadeQA/hDCAPosToPV", "Pos V0 daughter DCA to primary vertex", kTH1F, {DCAToPVAxis});
    mHistogramRegistry->add("CascadeQA/hDCANegToPV", "Neg V0 daughter DCA to primary vertex", kTH1F, {DCAToPVAxis});
    mHistogramRegistry->add("CascadeQA/hDCABachToPV", "Bachelor DCA to primary vertex", kTH1F, {DCAToPVAxis});
    mHistogramRegistry->add("CascadeQA/hDCAV0ToPV", "V0 DCA to primary vertex", kTH1F, {DCAToPVAxis});
  }

  /// check whether the most open cuts are fulfilled - most of this should have
  /// already be done by the filters
  nPtCascadeMinSel = getNSelections(femtoUniverseCascadeSelection::kCascadepTMin);
  nPtCascadeMaxSel = getNSelections(femtoUniverseCascadeSelection::kCascadepTMax);
  nEtaCascadeMaxSel = getNSelections(femtoUniverseCascadeSelection::kCascadeetaMax);
  nDCAV0DaughMax = getNSelections(femtoUniverseCascadeSelection::kCascadeV0DCADaughMax);
  nCPAV0Min = getNSelections(femtoUniverseCascadeSelection::kCascadeV0CPAMin);
  nTranRadV0Min = getNSelections(femtoUniverseCascadeSelection::kCascadeV0TranRadMin);
  nTranRadV0Max = getNSelections(femtoUniverseCascadeSelection::kCascadeV0TranRadMax);
  nV0DecVtxMax = getNSelections(femtoUniverseCascadeSelection::kCascadeV0DecVtxMax);
  nDCACascadeDaughMax = getNSelections(femtoUniverseCascadeSelection::kCascadeDCADaughMax);
  nCPACascadeMin = getNSelections(femtoUniverseCascadeSelection::kCascadeCPAMin);
  nTranRadCascadeMin = getNSelections(femtoUniverseCascadeSelection::kCascadeTranRadMin);
  nTranRadCascadeMax = getNSelections(femtoUniverseCascadeSelection::kCascadeTranRadMax);
  nDecVtxMax = getNSelections(femtoUniverseCascadeSelection::kCascadeDecVtxMax);
  nDCAPosToPV = getNSelections(femtoUniverseCascadeSelection::kCascadeDCAPosToPV);
  nDCANegToPV = getNSelections(femtoUniverseCascadeSelection::kCascadeDCANegToPV);
  nDCABachToPV = getNSelections(femtoUniverseCascadeSelection::kCascadeDCABachToPV);
  nDCAV0ToPV = getNSelections(femtoUniverseCascadeSelection::kCascadeDCAV0ToPV);
  // dodac V0 mass min i max

  pTCascadeMin = getMinimalSelection(femtoUniverseCascadeSelection::kCascadepTMin,
                                     femtoUniverseSelection::kLowerLimit);
  pTCascadeMax = getMinimalSelection(femtoUniverseCascadeSelection::kCascadepTMax,
                                     femtoUniverseSelection::kUpperLimit);
  etaCascadeMax = getMinimalSelection(femtoUniverseCascadeSelection::kCascadeetaMax,
                                      femtoUniverseSelection::kAbsUpperLimit);
  DCAV0DaughMax = getMinimalSelection(femtoUniverseCascadeSelection::kCascadeV0DCADaughMax,
                                      femtoUniverseSelection::kUpperLimit);
  CPAV0Min = getMinimalSelection(femtoUniverseCascadeSelection::kCascadeV0CPAMin,
                                 femtoUniverseSelection::kLowerLimit);
  TranRadV0Min = getMinimalSelection(femtoUniverseCascadeSelection::kCascadeV0TranRadMin,
                                     femtoUniverseSelection::kLowerLimit);
  TranRadV0Max = getMinimalSelection(femtoUniverseCascadeSelection::kCascadeV0TranRadMax,
                                     femtoUniverseSelection::kUpperLimit);
  V0DecVtxMax = getMinimalSelection(femtoUniverseCascadeSelection::kCascadeV0DecVtxMax,
                                    femtoUniverseSelection::kAbsUpperLimit);
  DCACascadeDaughMax = getMinimalSelection(femtoUniverseCascadeSelection::kCascadeDCADaughMax,
                                           femtoUniverseSelection::kUpperLimit);
  CPACascadeMin = getMinimalSelection(femtoUniverseCascadeSelection::kCascadeCPAMin,
                                      femtoUniverseSelection::kLowerLimit);
  TranRadCascadeMin = getMinimalSelection(femtoUniverseCascadeSelection::kCascadeTranRadMin,
                                          femtoUniverseSelection::kLowerLimit);
  TranRadCascadeMax = getMinimalSelection(femtoUniverseCascadeSelection::kCascadeTranRadMax,
                                          femtoUniverseSelection::kUpperLimit);
  DecVtxMax = getMinimalSelection(femtoUniverseCascadeSelection::kCascadeDecVtxMax,
                                  femtoUniverseSelection::kAbsUpperLimit);
  DCAPosToPV = getMinimalSelection(femtoUniverseCascadeSelection::kCascadeDCAPosToPV,
                                   femtoUniverseSelection::kLowerLimit);
  DCANegToPV = getMinimalSelection(femtoUniverseCascadeSelection::kCascadeDCANegToPV,
                                   femtoUniverseSelection::kLowerLimit);
  DCABachToPV = getMinimalSelection(femtoUniverseCascadeSelection::kCascadeDCABachToPV,
                                    femtoUniverseSelection::kLowerLimit);
  DCAV0ToPV = getMinimalSelection(femtoUniverseCascadeSelection::kCascadeDCAV0ToPV,
                                  femtoUniverseSelection::kLowerLimit);
  fV0InvMassLowLimit = getMinimalSelection(femtoUniverseCascadeSelection::kCascadeV0MassMin,
                                           femtoUniverseSelection::kLowerLimit);
  fV0InvMassUpLimit = getMinimalSelection(femtoUniverseCascadeSelection::kCascadeV0MassMax,
                                          femtoUniverseSelection::kUpperLimit);
}

template <typename Col, typename Casc, typename Track>
bool FemtoUniverseCascadeSelection::isSelectedMinimal(Col const& col, Casc const& cascade, Track const& posTrack, Track const& negTrack, Track const& bachTrack)
{
  const auto signPos = posTrack.sign();
  const auto signNeg = negTrack.sign();

  if (signPos < 0 || signNeg > 0) {
    LOG(warn) << "Something wrong in isSelectedMinimal";
    LOG(warn) << "ERROR - Wrong sign for V0 daughters";
  }

  const std::vector<float> decVtx = {cascade.x(), cascade.y(), cascade.z()};

  const float cpav0 = cascade.v0cosPA(col.posX(), col.posY(), col.posZ());
  const float cpaCasc = cascade.casccosPA(col.posX(), col.posY(), col.posZ());
  const float dcav0topv = cascade.dcav0topv(col.posX(), col.posY(), col.posZ());
  const float invMassLambda = cascade.mLambda();
  const float invMassXi = cascade.mXi();

  if (invMassLambda < fV0InvMassLowLimit || invMassLambda > fV0InvMassUpLimit) {
    return false;
  }
  if (invMassXi < fInvMassLowLimit || invMassXi > fInvMassUpLimit) {
    return false;
  }
  if (fRejectOmega) {
    const float invMassOmega = cascade.mOmega();
    if (invMassOmega > fInvMassOmegaLowLimit &&
        invMassOmega < fInvMassOmegaUpLimit) {
      return false;
    }
  }
  if (nPtCascadeMinSel > 0 && cascade.pt() < pTCascadeMin) {
    return false;
  }
  if (nPtCascadeMaxSel > 0 && cascade.pt() > pTCascadeMax) {
    return false;
  }
  if (nEtaCascadeMaxSel > 0 && std::abs(cascade.eta()) > etaCascadeMax) {
    return false;
  }
  if (nDCAV0DaughMax > 0 && cascade.dcaV0daughters() > DCAV0DaughMax) {
    return false;
  }
  if (nCPAV0Min > 0 && cpav0 < CPAV0Min) {
    return false;
  }
  if (nTranRadV0Min > 0 && cascade.v0radius() < TranRadV0Min) {
    return false;
  }
  if (nTranRadV0Max > 0 && cascade.v0radius() > TranRadV0Max) {
    return false;
  }
  if (nDCACascadeDaughMax > 0 && cascade.dcacascdaughters() > DCACascadeDaughMax) {
    return false;
  }
  if (nCPACascadeMin > 0 && cpaCasc < CPACascadeMin) {
    return false;
  }
  if (nTranRadCascadeMin > 0 && cascade.cascradius() < TranRadCascadeMin) {
    return false;
  }
  if (nTranRadCascadeMax > 0 && cascade.cascradius() > TranRadCascadeMax) {
    return false;
  }
  for (size_t i = 0; i < decVtx.size(); i++) {
    if (nDecVtxMax > 0 && decVtx.at(i) > DecVtxMax) {
      return false;
    }
  }
  if (nDCAPosToPV > 0 && abs(cascade.dcapostopv()) < DCAPosToPV) {
    return false;
  }
  if (nDCANegToPV > 0 && abs(cascade.dcanegtopv()) < DCANegToPV) {
    return false;
  }
  if (nDCABachToPV > 0 && abs(cascade.dcabachtopv()) < DCABachToPV) {
    return false;
  }
  if (nDCAV0ToPV > 0 && abs(dcav0topv) < DCAV0ToPV) {
    return false;
  }

  if (!PosDaughTrack.isSelectedMinimal(posTrack)) {
    return false;
  }
  if (!NegDaughTrack.isSelectedMinimal(negTrack)) {
    return false;
  }
  if (!BachTrack.isSelectedMinimal(bachTrack)) {
    return false;
  }
  /*
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
  */
  return true;
}

template <typename Col, typename Casc, typename Track>
void FemtoUniverseCascadeSelection::fillCascadeQA(Col const& col, Casc const& cascade, Track const& posTrack, Track const& negTrack)
{
  const auto signPos = posTrack.sign();
  const auto signNeg = negTrack.sign();
  if (signPos < 0 || signNeg > 0) {
    LOG(warn) << "Something wrong in isSelectedMinimal";
    LOG(warn) << "ERROR - Wrong sign for V0 daughters";
  }

  // const std::vector<float> decVtx = {cascade.x(), cascade.y(), cascade.z()};
  const float cpav0 = cascade.v0cosPA(col.posX(), col.posY(), col.posZ());
  const float cpaCasc = cascade.casccosPA(col.posX(), col.posY(), col.posZ());
  const float dcav0topv = cascade.dcav0topv(col.posX(), col.posY(), col.posZ());

  const float invMassLambda = cascade.mLambda();
  const float invMassXi = cascade.mXi();

  mHistogramRegistry->fill(HIST("CascadeQA/hInvMassV0Cut"), invMassLambda);
  mHistogramRegistry->fill(HIST("CascadeQA/hInvMassCascadeCut"), invMassXi);
  mHistogramRegistry->fill(HIST("CascadeQA/hCascadePt"), cascade.pt());
  mHistogramRegistry->fill(HIST("CascadeQA/hCascadeEta"), cascade.eta());
  mHistogramRegistry->fill(HIST("CascadeQA/hCascadePhi"), cascade.phi());
  mHistogramRegistry->fill(HIST("CascadeQA/hDCAV0Daugh"), cascade.dcaV0daughters());
  mHistogramRegistry->fill(HIST("CascadeQA/hV0CPA"), cpav0);
  mHistogramRegistry->fill(HIST("CascadeQA/hV0TranRad"), cascade.v0radius());
  mHistogramRegistry->fill(HIST("CascadeQA/hCascadeCPA"), cpaCasc);
  mHistogramRegistry->fill(HIST("CascadeQA/hDCACascadeDaugh"), cascade.dcacascdaughters());
  mHistogramRegistry->fill(HIST("CascadeQA/hCascadeTranRad"), cascade.cascradius());
  mHistogramRegistry->fill(HIST("CascadeQA/hDCAPosToPV"), cascade.dcapostopv());
  mHistogramRegistry->fill(HIST("CascadeQA/hDCANegToPV"), cascade.dcanegtopv());
  mHistogramRegistry->fill(HIST("CascadeQA/hDCABachToPV"), cascade.dcabachtopv());
  mHistogramRegistry->fill(HIST("CascadeQA/hDCAV0ToPV"), dcav0topv);

  // is this necessary
  /*
  bool write = true;
  for (size_t i = 0; i < decVtx.size(); i++) {
    write = write && (decVtx.at(i) < DecVtxMax);
  }
  if (write) {
    mHistogramRegistry->fill(HIST("CAscadeQA/hInvMassCascadeDecVtxMax"),
                             cascade.mXi());
  }
  */
}

template <typename Col, typename Casc, typename Track>
void FemtoUniverseCascadeSelection::fillQA(Col const& /*col*/, Casc const& /*cascade*/, Track const& posTrack, Track const& negTrack, Track const& bachTrack)
{
  PosDaughTrack.fillQA<aod::femtouniverseparticle::ParticleType::kV0Child,
                       aod::femtouniverseparticle::TrackType::kPosChild>(posTrack);
  NegDaughTrack.fillQA<aod::femtouniverseparticle::ParticleType::kV0Child,
                       aod::femtouniverseparticle::TrackType::kNegChild>(negTrack);
  BachTrack.fillQA<aod::femtouniverseparticle::ParticleType::kCascadeBachelor,
                   aod::femtouniverseparticle::TrackType::kBachelor>(bachTrack);
}

} // namespace o2::analysis::femtoUniverse

#endif // PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSECASCADESELECTION_H_
