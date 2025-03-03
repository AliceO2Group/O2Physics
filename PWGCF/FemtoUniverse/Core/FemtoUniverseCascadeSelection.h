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

/// \file FemtoUniverseCascadeSelection.h
/// \brief Definition of the femto_universe_cascade_selection
/// \author Valentina Mantovani Sarti, TU München valentina.mantovani-sarti@tum.de
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de
/// \author Luca Barioglio, TU München, luca.barioglio@cern.ch
/// \author Zuzanna Chochulska, WUT Warsaw & CTU Prague, zchochul@cern.ch
/// \author Barbara Chytla, WUT Warsaw, barbara.chytla@cern.ch
/// \author Shirajum Monira, WUT Warsaw, shirajum.monira@cern.ch

#ifndef PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSECASCADESELECTION_H_
#define PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSECASCADESELECTION_H_

#include <string>
#include <vector>
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseObjectSelection.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseSelection.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseTrackSelection.h"
#include "Common/Core/RecoDecay.h"
#include "Framework/HistogramRegistry.h"
#include "ReconstructionDataFormats/PID.h"

namespace o2::analysis::femto_universe
{
namespace femto_universe_cascade_selection
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
} // namespace femto_universe_cascade_selection

/// \class FemtoUniverseCascadeSelection
/// \brief Cut class to contain and execute all cuts applied to Cascades
class FemtoUniverseCascadeSelection
  : public FemtoUniverseObjectSelection<float, femto_universe_cascade_selection::CascadeSel>
{
  // do cutculatora

 public:
  FemtoUniverseCascadeSelection()
    : nPtCascadeMinSel(0), nPtCascadeMaxSel(0), nEtaCascadeMaxSel(0), nDCAV0DaughMax(0), nCPAV0Min(0), nTranRadV0Min(0), nTranRadV0Max(0), nV0DecVtxMax(0), nDCACascadeDaughMax(0), nCPACascadeMin(0), nTranRadCascadeMin(0), nTranRadCascadeMax(0), nDecVtxMax(0), nDCAPosToPV(0), nDCANegToPV(0), nDCABachToPV(0), nDCAV0ToPV(0), pTCascadeMin(9999999.), pTCascadeMax(-9999999.), etaCascadeMax(-9999999.), fDCAV0DaughMax(-9999999.), fCPAV0Min(9999999.), fTranRadV0Min(9999999.), fTranRadV0Max(-9999999.), fV0DecVtxMax(-9999999.), fDCACascadeDaughMax(-9999999.), fCPACascadeMin(9999999.), fTranRadCascadeMin(9999999.), fTranRadCascadeMax(-9999999.), fDecVtxMax(-9999999.), fDCAPosToPV(9999999.), fDCANegToPV(9999999.), fDCABachToPV(9999999.), fDCAV0ToPV(9999999.), fV0InvMassLowLimit(1.05), fV0InvMassUpLimit(1.3), fInvMassLowLimit(1.25), fInvMassUpLimit(1.4), fRejectCompetingMass(false), fInvMassCompetingLowLimit(1.5), fInvMassCompetingUpLimit(2.0), isCascOmega(false) /*, nSigmaPIDOffsetTPC(0.)*/
  {
  }

  /// Initializes histograms for the task
  template <o2::aod::femtouniverseparticle::ParticleType part, o2::aod::femtouniverseparticle::ParticleType daugh, o2::aod::femtouniverseparticle::ParticleType bach, typename CutContainerType>
  void init(HistogramRegistry* registry, bool isSelectCascOmega = false);

  template <typename Col, typename Casc, typename Track>
  bool isSelectedMinimal(Col const& col, Casc const& cascade, Track const& posTrack, Track const& negTrack, Track const& bachTrack);

  template <typename Col, typename Casc, typename Track>
  void fillCascadeQA(Col const& col, Casc const& cascade, Track const& posTrack, Track const& negTrack);

  template <typename Col, typename Casc, typename Track>
  void fillQA(Col const& col, Casc const& cascade, Track const& posTrack, Track const& negTrack, Track const& bachTrack);

  template <typename T1, typename T2>
  void setChildCuts(femto_universe_cascade_selection::ChildTrackType child, T1 selVal,
                    T2 selVar, femto_universe_selection::SelectionType selType)
  {
    if (child == femto_universe_cascade_selection::kPosTrack) {
      posDaughTrack.setSelection(selVal, selVar, selType);
    } else if (child == femto_universe_cascade_selection::kNegTrack) {
      negDaughTrack.setSelection(selVal, selVar, selType);
    } else if (child == femto_universe_cascade_selection::kBachTrack) {
      bachTrackSel.setSelection(selVal, selVar, selType);
    }
  }

  template <typename T>
  void setChildPIDSpecies(femto_universe_cascade_selection::ChildTrackType child,
                          T& pids)
  {
    if (child == femto_universe_cascade_selection::kPosTrack) {
      posDaughTrack.setPIDSpecies(pids);
    } else if (child == femto_universe_cascade_selection::kNegTrack) {
      negDaughTrack.setPIDSpecies(pids);
    } else if (child == femto_universe_cascade_selection::kBachTrack) {
      bachTrackSel.setPIDSpecies(pids);
    }
  }

  /// Helper function to obtain the name of a given selection criterion for consistent naming of the configurables
  /// \param iSel Track selection variable to be examined
  /// \param prefix Additional prefix for the name of the configurable
  /// \param suffix Additional suffix for the name of the configurable
  static std::string getSelectionName(femto_universe_cascade_selection::CascadeSel iSel,
                                      std::string_view prefix = "",
                                      std::string_view suffix = "")
  {
    std::string outString = static_cast<std::string>(prefix);
    outString += static_cast<std::string>(SelectionNames[iSel]);
    outString += suffix;
    return outString;
  }

  /// Helper function to obtain the helper string of a given selection criterion
  /// for consistent description of the configurables
  /// \param iSel Track selection variable to be examined
  /// \param prefix Additional prefix for the output of the configurable
  static std::string getSelectionHelper(femto_universe_cascade_selection::CascadeSel iSel,
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

  /// Set limit for the omega rejection on the invariant mass
  /// \param lowLimit Lower limit for the invariant mass distribution
  /// \param upLimit Upper limit for the invariant mass distribution
  void setCompetingInvMassLimits(float lowLimit, float upLimit)
  {
    fRejectCompetingMass = true;
    fInvMassCompetingLowLimit = lowLimit;
    fInvMassCompetingUpLimit = upLimit;
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
  float fDCAV0DaughMax;
  float fCPAV0Min;
  float fTranRadV0Min;
  float fTranRadV0Max;
  float fV0DecVtxMax;
  float fDCACascadeDaughMax;
  float fCPACascadeMin;
  float fTranRadCascadeMin;
  float fTranRadCascadeMax;
  float fDecVtxMax;
  float fDCAPosToPV;
  float fDCANegToPV;
  float fDCABachToPV;
  float fDCAV0ToPV;

  float fV0InvMassLowLimit;
  float fV0InvMassUpLimit;

  float fInvMassLowLimit;
  float fInvMassUpLimit;

  float fRejectCompetingMass;
  float fInvMassCompetingLowLimit;
  float fInvMassCompetingUpLimit;

  bool isCascOmega;

  // float nSigmaPIDOffsetTPC;

  FemtoUniverseTrackSelection posDaughTrack;
  FemtoUniverseTrackSelection negDaughTrack;
  FemtoUniverseTrackSelection bachTrackSel;

  static constexpr int kNcascadeSelection = 20; // can I do less ?

  static constexpr std::string_view SelectionNames[kNcascadeSelection] = {
    "Sign", "PtMin", "PtMax", "EtaMax", "DCAv0daughMax", "v0CPAMin",
    "v0TranRadMin", "v0TranRadMax", "v0DecVecMax", "DCAcascDaugh",
    "CPAMin", "TranRadMin", "TranRadMax", "DecVtxMax",
    "DCAPosToPV", "DCANegToPV", "DCABachToPV", "DCAV0ToPV",
    "kV0MassMin", "V0MassMax"}; ///< Name of the different
                                ///< selections

  static constexpr femto_universe_selection::SelectionType
    mSelectionTypes[kNcascadeSelection]{
      femto_universe_selection::kEqual,      // sign
      femto_universe_selection::kLowerLimit, // pt min
      femto_universe_selection::kUpperLimit, // pt max
      femto_universe_selection::kUpperLimit, // eta max
      femto_universe_selection::kUpperLimit, // DCA v0 daughters max
      femto_universe_selection::kLowerLimit, // v0 cos PA min
      femto_universe_selection::kLowerLimit, // v0 tran rad min
      femto_universe_selection::kUpperLimit, // v0 tran rad max
      femto_universe_selection::kUpperLimit, // v0 maximum distance of decay vertex to PV
      femto_universe_selection::kUpperLimit, // DCA cascade daughters max
      femto_universe_selection::kLowerLimit, // cascade cos PA min
      femto_universe_selection::kLowerLimit, // cascade tran rad min
      femto_universe_selection::kUpperLimit, // cascade tran rad max
      femto_universe_selection::kUpperLimit, // cascade maximum distance of decay vertex to PV
      femto_universe_selection::kLowerLimit, // DCA pos to PV max
      femto_universe_selection::kLowerLimit, // DCA neg to PV max
      femto_universe_selection::kLowerLimit, // DCA bach to PV max
      femto_universe_selection::kLowerLimit, // DCA v0 to PV max
      femto_universe_selection::kLowerLimit, // v0 mass min
      femto_universe_selection::kUpperLimit, // v0 mass max
    }; ///< Map to match a variable with
       ///< its type

  static constexpr std::string_view SelectionHelper[kNcascadeSelection] = {
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

}; // namespace femto_universe

template <o2::aod::femtouniverseparticle::ParticleType part, o2::aod::femtouniverseparticle::ParticleType daugh, o2::aod::femtouniverseparticle::ParticleType bach, typename CutContainerType>
void FemtoUniverseCascadeSelection::init(HistogramRegistry* registry, bool isSelectCascOmega)
{

  if (registry) {
    mHistogramRegistry = registry;
    fillSelectionHistogram<part>();  // cascade
    fillSelectionHistogram<daugh>(); // pos, neg
    fillSelectionHistogram<bach>();  // bach

    AxisSpec massAxisCascade = {2200, 1.25f, 1.8f, "m_{Cascade} (GeV/#it{c}^{2})"};
    AxisSpec massAxisV0 = {600, 0.0f, 3.0f, "m_{V0} (GeV/#it{c}^{2})"};
    AxisSpec aDCADaughAxis = {1000, 0.0f, 2.0f, "DCA (cm)"};
    AxisSpec aDCAToPVAxis = {1000, -10.0f, 10.0f, "DCA to PV (cm)"};
    AxisSpec ptAxis = {100, 0.0f, 10.0f, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec etaAxis = {100, -2.0f, 2.0f, "#it{#eta}"};
    AxisSpec phiAxis = {100, 0.0f, 6.0f, "#it{#phi}"};
    AxisSpec aCPAAxis = {1000, 0.95f, 1.0f, "#it{cos #theta_{p}}"};
    AxisSpec tranRadAxis = {1000, 0.0f, 100.0f, "#it{r}_{xy} (cm)"};

    /// \todo this should be an automatic check in the parent class, and the
    /// return type should be templated
    size_t nSelections = getNSelections();
    if (nSelections > 17 * sizeof(CutContainerType)) {
      LOG(fatal) << "FemtoUniverseCascadeCuts: Number of selections to large for your "
                    "container - quitting!";
    }

    posDaughTrack.init<aod::femtouniverseparticle::ParticleType::kV0Child,
                       aod::femtouniverseparticle::TrackType::kPosChild,
                       aod::femtouniverseparticle::CutContainerType>(
      mHistogramRegistry);
    negDaughTrack.init<aod::femtouniverseparticle::ParticleType::kV0Child,
                       aod::femtouniverseparticle::TrackType::kNegChild,
                       aod::femtouniverseparticle::CutContainerType>(
      mHistogramRegistry);
    bachTrackSel.init<aod::femtouniverseparticle::ParticleType::kCascadeBachelor,
                      aod::femtouniverseparticle::TrackType::kBachelor,
                      aod::femtouniverseparticle::CutContainerType>(
      mHistogramRegistry);

    // V0 (Lambda)
    // mHistogramRegistry->add("CascadeQA/hInvMassV0NoCuts", "No cuts", kTH1F, {massAxisV0});
    mHistogramRegistry->add("CascadeQA/hInvMassV0Cut", "Invariant mass cut", kTH1F, {massAxisV0});
    mHistogramRegistry->add("CascadeQA/hDCAV0Daugh", "V0-daughters DCA", kTH1F, {aDCADaughAxis});
    mHistogramRegistry->add("CascadeQA/hV0CPA", "V0 cos PA", kTH1F, {aCPAAxis});
    mHistogramRegistry->add("CascadeQA/hV0TranRad", "V0 transverse radius", kTH1F, {tranRadAxis});
    // mHistogramRegistry->add("CascadeQA/hV0DecVtxMax", "V0 maximum distance on decay vertex", kTH1F, {massAxisV0});

    // Cascade (Xi, Omega)
    // mHistogramRegistry->add("CascadeQA/hInvMassCascadeNoCuts", "No cuts", kTH1F, {massAxisCascade});
    mHistogramRegistry->add("CascadeQA/hInvMassCascadeCut", "Invariant mass with cut", kTH1F, {massAxisCascade});
    mHistogramRegistry->add("CascadeQA/hCascadePt", "pT distribution", kTH1F, {ptAxis});
    mHistogramRegistry->add("CascadeQA/hCascadeEta", "Eta distribution", kTH1F, {etaAxis});
    mHistogramRegistry->add("CascadeQA/hCascadePhi", "Phi distribution", kTH1F, {phiAxis});
    mHistogramRegistry->add("CascadeQA/hDCACascadeDaugh", "Cascade-daughters DCA", kTH1F, {aDCADaughAxis});
    mHistogramRegistry->add("CascadeQA/hCascadeCPA", "Cos PA", kTH1F, {aCPAAxis});
    mHistogramRegistry->add("CascadeQA/hCascadeTranRad", "Transverse radius", kTH1F, {tranRadAxis});
    mHistogramRegistry->add("CascadeQA/hDCAPosToPV", "Pos V0 daughter DCA to primary vertex", kTH1F, {aDCAToPVAxis});
    mHistogramRegistry->add("CascadeQA/hDCANegToPV", "Neg V0 daughter DCA to primary vertex", kTH1F, {aDCAToPVAxis});
    mHistogramRegistry->add("CascadeQA/hDCABachToPV", "Bachelor DCA to primary vertex", kTH1F, {aDCAToPVAxis});
    mHistogramRegistry->add("CascadeQA/hDCAV0ToPV", "V0 DCA to primary vertex", kTH1F, {aDCAToPVAxis});
  }

  /// check whether the most open cuts are fulfilled - most of this should have
  /// already be done by the filters
  nPtCascadeMinSel = getNSelections(femto_universe_cascade_selection::kCascadepTMin);
  nPtCascadeMaxSel = getNSelections(femto_universe_cascade_selection::kCascadepTMax);
  nEtaCascadeMaxSel = getNSelections(femto_universe_cascade_selection::kCascadeetaMax);
  nDCAV0DaughMax = getNSelections(femto_universe_cascade_selection::kCascadeV0DCADaughMax);
  nCPAV0Min = getNSelections(femto_universe_cascade_selection::kCascadeV0CPAMin);
  nTranRadV0Min = getNSelections(femto_universe_cascade_selection::kCascadeV0TranRadMin);
  nTranRadV0Max = getNSelections(femto_universe_cascade_selection::kCascadeV0TranRadMax);
  nV0DecVtxMax = getNSelections(femto_universe_cascade_selection::kCascadeV0DecVtxMax);
  nDCACascadeDaughMax = getNSelections(femto_universe_cascade_selection::kCascadeDCADaughMax);
  nCPACascadeMin = getNSelections(femto_universe_cascade_selection::kCascadeCPAMin);
  nTranRadCascadeMin = getNSelections(femto_universe_cascade_selection::kCascadeTranRadMin);
  nTranRadCascadeMax = getNSelections(femto_universe_cascade_selection::kCascadeTranRadMax);
  nDecVtxMax = getNSelections(femto_universe_cascade_selection::kCascadeDecVtxMax);
  nDCAPosToPV = getNSelections(femto_universe_cascade_selection::kCascadeDCAPosToPV);
  nDCANegToPV = getNSelections(femto_universe_cascade_selection::kCascadeDCANegToPV);
  nDCABachToPV = getNSelections(femto_universe_cascade_selection::kCascadeDCABachToPV);
  nDCAV0ToPV = getNSelections(femto_universe_cascade_selection::kCascadeDCAV0ToPV);
  // dodac V0 mass min i max

  pTCascadeMin = getMinimalSelection(femto_universe_cascade_selection::kCascadepTMin,
                                     femto_universe_selection::kLowerLimit);
  pTCascadeMax = getMinimalSelection(femto_universe_cascade_selection::kCascadepTMax,
                                     femto_universe_selection::kUpperLimit);
  etaCascadeMax = getMinimalSelection(femto_universe_cascade_selection::kCascadeetaMax,
                                      femto_universe_selection::kAbsUpperLimit);
  fDCAV0DaughMax = getMinimalSelection(femto_universe_cascade_selection::kCascadeV0DCADaughMax,
                                       femto_universe_selection::kUpperLimit);
  fCPAV0Min = getMinimalSelection(femto_universe_cascade_selection::kCascadeV0CPAMin,
                                  femto_universe_selection::kLowerLimit);
  fTranRadV0Min = getMinimalSelection(femto_universe_cascade_selection::kCascadeV0TranRadMin,
                                      femto_universe_selection::kLowerLimit);
  fTranRadV0Max = getMinimalSelection(femto_universe_cascade_selection::kCascadeV0TranRadMax,
                                      femto_universe_selection::kUpperLimit);
  fV0DecVtxMax = getMinimalSelection(femto_universe_cascade_selection::kCascadeV0DecVtxMax,
                                     femto_universe_selection::kAbsUpperLimit);
  fDCACascadeDaughMax = getMinimalSelection(femto_universe_cascade_selection::kCascadeDCADaughMax,
                                            femto_universe_selection::kUpperLimit);
  fCPACascadeMin = getMinimalSelection(femto_universe_cascade_selection::kCascadeCPAMin,
                                       femto_universe_selection::kLowerLimit);
  fTranRadCascadeMin = getMinimalSelection(femto_universe_cascade_selection::kCascadeTranRadMin,
                                           femto_universe_selection::kLowerLimit);
  fTranRadCascadeMax = getMinimalSelection(femto_universe_cascade_selection::kCascadeTranRadMax,
                                           femto_universe_selection::kUpperLimit);
  fDecVtxMax = getMinimalSelection(femto_universe_cascade_selection::kCascadeDecVtxMax,
                                   femto_universe_selection::kAbsUpperLimit);
  fDCAPosToPV = getMinimalSelection(femto_universe_cascade_selection::kCascadeDCAPosToPV,
                                    femto_universe_selection::kLowerLimit);
  fDCANegToPV = getMinimalSelection(femto_universe_cascade_selection::kCascadeDCANegToPV,
                                    femto_universe_selection::kLowerLimit);
  fDCABachToPV = getMinimalSelection(femto_universe_cascade_selection::kCascadeDCABachToPV,
                                     femto_universe_selection::kLowerLimit);
  fDCAV0ToPV = getMinimalSelection(femto_universe_cascade_selection::kCascadeDCAV0ToPV,
                                   femto_universe_selection::kLowerLimit);
  fV0InvMassLowLimit = getMinimalSelection(femto_universe_cascade_selection::kCascadeV0MassMin,
                                           femto_universe_selection::kLowerLimit);
  fV0InvMassUpLimit = getMinimalSelection(femto_universe_cascade_selection::kCascadeV0MassMax,
                                          femto_universe_selection::kUpperLimit);

  isCascOmega = isSelectCascOmega;
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
  const float invMass = isCascOmega ? cascade.mOmega() : cascade.mXi();
  const float nSigmaPIDMax = bachTrackSel.getSigmaPIDMax();

  if (invMassLambda < fV0InvMassLowLimit || invMassLambda > fV0InvMassUpLimit) {
    return false;
  }
  if (invMass < fInvMassLowLimit || invMass > fInvMassUpLimit) {
    return false;
  }
  if (fRejectCompetingMass) {
    const float invMassCompeting = isCascOmega ? cascade.mXi() : cascade.mOmega();
    if (invMassCompeting > fInvMassCompetingLowLimit &&
        invMassCompeting < fInvMassCompetingUpLimit) {
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
  if (nDCAV0DaughMax > 0 && cascade.dcaV0daughters() > fDCAV0DaughMax) {
    return false;
  }
  if (nCPAV0Min > 0 && cpav0 < fCPAV0Min) {
    return false;
  }
  if (nTranRadV0Min > 0 && cascade.v0radius() < fTranRadV0Min) {
    return false;
  }
  if (nTranRadV0Max > 0 && cascade.v0radius() > fTranRadV0Max) {
    return false;
  }
  if (nDCACascadeDaughMax > 0 && cascade.dcacascdaughters() > fDCACascadeDaughMax) {
    return false;
  }
  if (nCPACascadeMin > 0 && cpaCasc < fCPACascadeMin) {
    return false;
  }
  if (nTranRadCascadeMin > 0 && cascade.cascradius() < fTranRadCascadeMin) {
    return false;
  }
  if (nTranRadCascadeMax > 0 && cascade.cascradius() > fTranRadCascadeMax) {
    return false;
  }
  for (size_t i = 0; i < decVtx.size(); i++) {
    if (nDecVtxMax > 0 && decVtx.at(i) > fDecVtxMax) {
      return false;
    }
  }
  if (nDCAPosToPV > 0 && std::abs(cascade.dcapostopv()) < fDCAPosToPV) {
    return false;
  }
  if (nDCANegToPV > 0 && std::abs(cascade.dcanegtopv()) < fDCANegToPV) {
    return false;
  }
  if (nDCABachToPV > 0 && std::abs(cascade.dcabachtopv()) < fDCABachToPV) {
    return false;
  }
  if (nDCAV0ToPV > 0 && std::abs(dcav0topv) < fDCAV0ToPV) {
    return false;
  }

  if (!posDaughTrack.isSelectedMinimal(posTrack)) {
    return false;
  }
  if (!negDaughTrack.isSelectedMinimal(negTrack)) {
    return false;
  }
  if (bachTrack.hasTOF() && ((isCascOmega && bachTrack.tofNSigmaKa() > nSigmaPIDMax) || (!isCascOmega && bachTrack.tofNSigmaPi() > nSigmaPIDMax))) {
    return false;
  }
  if (!bachTrackSel.isSelectedMinimal(bachTrack)) {
    return false;
  }
  /*
    // check that track combinations for V0 or antiV0 would be fulfilling PID
    float nSigmaPIDMax = posDaughTrack.getSigmaPIDMax();
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
  const float invMass = isCascOmega ? cascade.mOmega() : cascade.mXi();

  mHistogramRegistry->fill(HIST("CascadeQA/hInvMassV0Cut"), invMassLambda);
  mHistogramRegistry->fill(HIST("CascadeQA/hInvMassCascadeCut"), invMass);
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
  posDaughTrack.fillQA<aod::femtouniverseparticle::ParticleType::kV0Child,
                       aod::femtouniverseparticle::TrackType::kPosChild>(posTrack);
  negDaughTrack.fillQA<aod::femtouniverseparticle::ParticleType::kV0Child,
                       aod::femtouniverseparticle::TrackType::kNegChild>(negTrack);
  bachTrackSel.fillQA<aod::femtouniverseparticle::ParticleType::kCascadeBachelor,
                      aod::femtouniverseparticle::TrackType::kBachelor>(bachTrack);
}

} // namespace o2::analysis::femto_universe

#endif // PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSECASCADESELECTION_H_
