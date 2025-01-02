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

/// \file femtoDreamCascadeSelection.h
/// \brief Definition of the femtoDreamCascadeSelection
/// \author Valentina Mantovani Sarti, TU München valentina.mantovani-sarti@tum.de
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de
/// \author Luca Barioglio, TU München, luca.barioglio@cern.ch
/// \author Zuzanna Chochulska, WUT Warsaw & CTU Prague, zchochul@cern.ch
/// \author Barbara Chytla, WUT Warsaw, barbara.chytla@cern.ch
/// \author Shirajum Monira, WUT Warsaw, shirajum.monira@cern.ch

#ifndef PWGCF_FEMTODREAM_CORE_FEMTODREAMCASCADESELECTION_H_
#define PWGCF_FEMTODREAM_CORE_FEMTODREAMCASCADESELECTION_H_

#include <iostream>
#include <string>
#include <vector>

#include <TDatabasePDG.h> // FIXME

#include "PWGCF/FemtoDream/Core/femtoDreamObjectSelection.h"
#include "PWGCF/FemtoDream/Core/femtoDreamSelection.h"
#include "PWGCF/FemtoDream/Core/femtoDreamTrackSelection.h"

#include "Common/Core/RecoDecay.h"
#include "Framework/HistogramRegistry.h"
#include "ReconstructionDataFormats/PID.h"

using namespace o2::framework;

namespace o2::analysis::femtoDream
{
namespace femtoDreamCascadeSelection
{
/// The different selections this task is capable of doing
enum CascadeSel {
  kCascadeSign, ///< +1 particle, -1 antiparticle
  kCascadepTMin,
  kCascadepTMax,
};

enum ChildTrackType { kPosTrack,
                      kNegTrack,
                      kBachTrack };

enum CascadeContainerPosition {
  kCascade,
  kPosCuts,
  kPosPID,
  kNegCuts,
  kNegPID,
  kBachCuts,
  kBachPID,
}; /// Position in the full VO cut container (for cutculator)

} // namespace femtoDreamCascadeSelection

/// \class FemtoDreamCascadeSelection
/// \brief Cut class to contain and execute all cuts applied to Cascades
class FemtoDreamCascadeSelection
  : public FemtoDreamObjectSelection<float, femtoDreamCascadeSelection::CascadeSel>
{
 public:
  FemtoDreamCascadeSelection()
    : nPtCascadeMinSel(0), nPtCascadeMaxSel(0), pTCascadeMin(9999999.), pTCascadeMax(-9999999.) , fInvMassLowLimit(1.25), fInvMassUpLimit(1.4), fRejectCompetingMass(false), fInvMassCompetingLowLimit(1.5), fInvMassCompetingUpLimit(2.0), isCascOmega(false)
  {
  }

  /// Initializes histograms for the task
  template <o2::aod::femtodreamparticle::ParticleType part,
            o2::aod::femtodreamparticle::ParticleType v0daugh,
            o2::aod::femtodreamparticle::ParticleType bach,
            typename cutContainerType>
  void init(HistogramRegistry* QAregistry, HistogramRegistry* Registry, bool isSelectCascOmega = false);

  template <typename Col, typename Casc, typename Track>
  bool isSelectedMinimal(Col const& col, Casc const& cascade, Track const& posTrack, Track const& negTrack, Track const& bachTrack);

  template <typename Col, typename Casc, typename Track>
  void fillCascadeQA(Col const& col, Casc const& cascade, Track const& posTrack, Track const& negTrack);

  template <o2::aod::femtodreamparticle::ParticleType part, typename Col, typename Casc, typename Track>
  void fillQA(Col const& col, Casc const& cascade, Track const& posTrack, Track const& negTrack, Track const& bachTrack);

  template <typename cutContainerType, typename Col, typename Casc, typename Track>
  std::array<cutContainerType, 8> getCutContainer(Col const& /*col*/, Casc const& casc, Track const& posTrack, Track const& negTrack, Track const& bachTrack);

  template <typename T1, typename T2>
  void setChildCuts(femtoDreamCascadeSelection::ChildTrackType child,
                    T1 selVal,
                    T2 selVar,
                    femtoDreamSelection::SelectionType selType)
  {
    if (child == femtoDreamCascadeSelection::kPosTrack) {
      PosDaughTrack.setSelection(selVal, selVar, selType);
    } else if (child == femtoDreamCascadeSelection::kNegTrack) {
      NegDaughTrack.setSelection(selVal, selVar, selType);
    } else if (child == femtoDreamCascadeSelection::kBachTrack) {
      BachTrack.setSelection(selVal, selVar, selType);
    }
  }

  template <typename T>
  void setChildPIDSpecies(femtoDreamCascadeSelection::ChildTrackType child,
                          T& pids)
  {
    if (child == femtoDreamCascadeSelection::kPosTrack) {
      PosDaughTrack.setPIDSpecies(pids);
    } else if (child == femtoDreamCascadeSelection::kNegTrack) {
      NegDaughTrack.setPIDSpecies(pids);
    } else if (child == femtoDreamCascadeSelection::kBachTrack) {
      BachTrack.setPIDSpecies(pids);
    }
  }

  /// Helper function to obtain the name of a given selection criterion for consistent naming of the configurables
  /// \param iSel Track selection variable to be examined
  /// \param prefix Additional prefix for the name of the configurable
  /// \param suffix Additional suffix for the name of the configurable
  static std::string getSelectionName(femtoDreamCascadeSelection::CascadeSel iSel,
                                      std::string_view prefix = "",
                                      std::string_view suffix = "")
  {
    std::string outString = static_cast<std::string>(prefix);
    outString += static_cast<std::string>(mSelectionNames[iSel]);
    outString += suffix;
    return outString;
  }
  
  /// Helper function to obtain the index of a given selection variable for consistent naming of the configurables
  /// \param obs Cascade selection variable (together with prefix) got from file
  /// \param prefix Additional prefix for the output of the configurable
  static int findSelectionIndex(const std::string_view& obs,
                                std::string_view prefix = "")
  {
    for (int index = 0; index < kNcascadeSelection; index++) {
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
  /// \param iSel Casc selection variable whose type is returned
  static femtoDreamSelection::SelectionType getSelectionType(femtoDreamCascadeSelection::CascadeSel iSel)
  {
    return mSelectionTypes[iSel];
  }

  /// Helper function to obtain the helper string of a given selection criterion
  /// for consistent description of the configurables
  /// \param iSel Track selection variable to be examined
  /// \param prefix Additional prefix for the output of the configurable
  static std::string getSelectionHelper(femtoDreamCascadeSelection::CascadeSel iSel,
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
  void setCompetingInvMassLimits(float lowLimit, float upLimit)
  {
    fRejectCompetingMass = true;
    fInvMassCompetingLowLimit = lowLimit;
    fInvMassCompetingUpLimit = upLimit;
  }

 private:
  int nPtCascadeMinSel;
  int nPtCascadeMaxSel;
  float pTCascadeMin;
  float pTCascadeMax;
  
  float fInvMassLowLimit;
  float fInvMassUpLimit;

  float fRejectCompetingMass;
  float fInvMassCompetingLowLimit;
  float fInvMassCompetingUpLimit;

  bool isCascOmega;

  // float nSigmaPIDOffsetTPC;

  FemtoDreamTrackSelection PosDaughTrack;
  FemtoDreamTrackSelection NegDaughTrack;
  FemtoDreamTrackSelection BachTrack;

  static constexpr int kNcascadeSelection = 3; // can I do less ?

  static constexpr std::string_view mSelectionNames[kNcascadeSelection] = {
    "Sign", "PtMin", "PtMax"};
    /*
    "Sign", "PtMin", "PtMax", "EtaMax", "DCAv0daughMax", "v0CPAMin",
    "v0TranRadMin", "v0TranRadMax", "v0DecVecMax", "DCAcascDaugh",
    "CPAMin", "TranRadMin", "TranRadMax", "DecVtxMax",
    "DCAPosToPV", "DCANegToPV", "DCABachToPV", "DCAV0ToPV",
    "kV0MassMin", "V0MassMax"}; ///< Name of the different
                                ///< selections
    */ 

  static constexpr femtoDreamSelection::SelectionType
    mSelectionTypes[kNcascadeSelection]{
      femtoDreamSelection::kEqual,      // sign
      femtoDreamSelection::kLowerLimit, // pt min
      femtoDreamSelection::kUpperLimit, // pt max
    }; ///< Map to match a variable with
       ///< its type

  static constexpr std::string_view mSelectionHelper[kNcascadeSelection] = {
    "Cascade particle sign (+1 or -1)",
    "Minimum pT (GeV/c)",
    "Maximum pT (GeV/c)"}; ///< Helper information for the
                        ///< different selections

}; // namespace femtoDream

template <o2::aod::femtodreamparticle::ParticleType part, o2::aod::femtodreamparticle::ParticleType daugh, o2::aod::femtodreamparticle::ParticleType bach, typename cutContainerType>
void FemtoDreamCascadeSelection::init(HistogramRegistry* QAregistry, HistogramRegistry* Registry, bool isSelectCascOmega)
{

  if (QAregistry && Registry) {
    mHistogramRegistry = Registry;
    mQAHistogramRegistry = QAregistry;
    fillSelectionHistogram<part>();  // cascade
    fillSelectionHistogram<daugh>(); // pos, neg
    fillSelectionHistogram<bach>();  // bach

    AxisSpec massAxisCascade = {2200, 1.25f, 1.8f, "m_{#Cascade} (GeV/#it{c}^{2})"};
    AxisSpec ptAxis = {100, 0.0f, 10.0f, "#it{p}_{T} (GeV/#it{c})"};
    /*
    AxisSpec massAxisV0 = {600, 0.0f, 3.0f, "m_{#V0} (GeV/#it{c}^{2})"};
    AxisSpec DCADaughAxis = {1000, 0.0f, 2.0f, "DCA (cm)"};
    AxisSpec DCAToPVAxis = {1000, -10.0f, 10.0f, "DCA to PV (cm)"};
    AxisSpec etaAxis = {100, -2.0f, 2.0f, "#it{#eta}"};
    AxisSpec phiAxis = {100, 0.0f, 6.0f, "#it{#phi}"};
    AxisSpec CPAAxis = {1000, 0.95f, 1.0f, "#it{cos #theta_{p}}"};
    AxisSpec tranRadAxis = {1000, 0.0f, 100.0f, "#it{r}_{xy} (cm)"};
    */

    /// \todo this should be an automatic check in the parent class, and the
    /// return type should be templated
    size_t nSelections = getNSelections();
    if (nSelections > 3 * sizeof(cutContainerType)) {
      LOG(fatal) << "FemtoDreamCascadeCuts: Number of selections to large for your "
                    "container - quitting!";
    }

    std::string folderName = static_cast<std::string>(
      o2::aod::femtodreamparticle::ParticleTypeName[part]);
    mQAHistogramRegistry->add((folderName + "/hPt").c_str(),
                              "; #it{p}_{T} (GeV/#it{c}); Entries", kTH1F,
                              {{1000, 0, 10}});
    mQAHistogramRegistry->add((folderName + "/hEta").c_str(), "; #eta; Entries",
                              kTH1F, {{1000, -1, 1}});
    mQAHistogramRegistry->add((folderName + "/hPhi").c_str(), "; #phi; Entries",
                              kTH1F, {{1000, 0, 2. * M_PI}});

    PosDaughTrack.init<aod::femtodreamparticle::ParticleType::kCascadeV0Child,
                       aod::femtodreamparticle::TrackType::kPosChild,
                       aod::femtodreamparticle::cutContainerType>(
      mQAHistogramRegistry, mQAHistogramRegistry);
    NegDaughTrack.init<aod::femtodreamparticle::ParticleType::kCascadeV0Child,
                       aod::femtodreamparticle::TrackType::kNegChild,
                       aod::femtodreamparticle::cutContainerType>(
      mQAHistogramRegistry, mQAHistogramRegistry);
    BachTrack.init<aod::femtodreamparticle::ParticleType::kCascadeBachelor,
                   aod::femtodreamparticle::TrackType::kBachelor,
                   aod::femtodreamparticle::cutContainerType>(
      mQAHistogramRegistry, mQAHistogramRegistry);

    // V0 (Lambda)
    //mQAHistogramRegistry->add("CascadeQA/hInvMassV0NoCuts", "No cuts", kTH1F, {massAxisV0});
    // mQAHistogramRegistry->add("CascadeQA/hInvMassV0Cut", "Invariant mass cut", kTH1F, {massAxisV0});
    // mQAHistogramRegistry->add("CascadeQA/hDCAV0Daugh", "V0-daughters DCA", kTH1F, {DCADaughAxis});
    // mQAHistogramRegistry->add("CascadeQA/hV0CPA", "V0 cos PA", kTH1F, {CPAAxis});
    // mQAHistogramRegistry->add("CascadeQA/hV0TranRad", "V0 transverse radius", kTH1F, {tranRadAxis});
    // mQAHistogramRegistry->add("CascadeQA/hV0DecVtxMax", "V0 maximum distance on decay vertex", kTH1F, {massAxisV0});

    // Cascade (Xi, Omega)
    // mQAHistogramRegistry->add("CascadeQA/hInvMassCascadeNoCuts", "No cuts", kTH1F, {massAxisCascade});
    mQAHistogramRegistry->add("CascadeQA/hInvMassCascadeCut", "Invariant mass with cut", kTH1F, {massAxisCascade});
    mQAHistogramRegistry->add("CascadeQA/hCascadePt", "pT distribution", kTH1F, {ptAxis});
    // mQAHistogramRegistry->add("CascadeQA/hCascadeEta", "Eta distribution", kTH1F, {etaAxis});
    // mQAHistogramRegistry->add("CascadeQA/hCascadePhi", "Phi distribution", kTH1F, {phiAxis});
    // mQAHistogramRegistry->add("CascadeQA/hDCACascadeDaugh", "Cascade-daughters DCA", kTH1F, {DCADaughAxis});
    // mQAHistogramRegistry->add("CascadeQA/hCascadeCPA", "Cos PA", kTH1F, {CPAAxis});
    // mQAHistogramRegistry->add("CascadeQA/hCascadeTranRad", "Transverse radius", kTH1F, {tranRadAxis});
    // mQAHistogramRegistry->add("CascadeQA/hDCAPosToPV", "Pos V0 daughter DCA to primary vertex", kTH1F, {DCAToPVAxis});
    // mQAHistogramRegistry->add("CascadeQA/hDCANegToPV", "Neg V0 daughter DCA to primary vertex", kTH1F, {DCAToPVAxis});
    // mQAHistogramRegistry->add("CascadeQA/hDCABachToPV", "Bachelor DCA to primary vertex", kTH1F, {DCAToPVAxis});
    // mQAHistogramRegistry->add("CascadeQA/hDCAV0ToPV", "V0 DCA to primary vertex", kTH1F, {DCAToPVAxis});
  }

  /// check whether the most open cuts are fulfilled - most of this should have
  /// already be done by the filters
  nPtCascadeMinSel = getNSelections(femtoDreamCascadeSelection::kCascadepTMin);
  nPtCascadeMaxSel = getNSelections(femtoDreamCascadeSelection::kCascadepTMax);
  // nEtaCascadeMaxSel = getNSelections(femtoDreamCascadeSelection::kCascadeetaMax);
  // nDCAV0DaughMax = getNSelections(femtoDreamCascadeSelection::kCascadeV0DCADaughMax);
  // nCPAV0Min = getNSelections(femtoDreamCascadeSelection::kCascadeV0CPAMin);
  // nTranRadV0Min = getNSelections(femtoDreamCascadeSelection::kCascadeV0TranRadMin);
  // nTranRadV0Max = getNSelections(femtoDreamCascadeSelection::kCascadeV0TranRadMax);
  // nV0DecVtxMax = getNSelections(femtoDreamCascadeSelection::kCascadeV0DecVtxMax);
  // nDCACascadeDaughMax = getNSelections(femtoDreamCascadeSelection::kCascadeDCADaughMax);
  // nCPACascadeMin = getNSelections(femtoDreamCascadeSelection::kCascadeCPAMin);
  // nTranRadCascadeMin = getNSelections(femtoDreamCascadeSelection::kCascadeTranRadMin);
  // nTranRadCascadeMax = getNSelections(femtoDreamCascadeSelection::kCascadeTranRadMax);
  // nDecVtxMax = getNSelections(femtoDreamCascadeSelection::kCascadeDecVtxMax);
  // nDCAPosToPV = getNSelections(femtoDreamCascadeSelection::kCascadeDCAPosToPV);
  // nDCANegToPV = getNSelections(femtoDreamCascadeSelection::kCascadeDCANegToPV);
  // nDCABachToPV = getNSelections(femtoDreamCascadeSelection::kCascadeDCABachToPV);
  // nDCAV0ToPV = getNSelections(femtoDreamCascadeSelection::kCascadeDCAV0ToPV);
  // dodac V0 mass min i max

  pTCascadeMin = getMinimalSelection(femtoDreamCascadeSelection::kCascadepTMin,
                                     femtoDreamSelection::kLowerLimit);
  pTCascadeMax = getMinimalSelection(femtoDreamCascadeSelection::kCascadepTMax,
                                     femtoDreamSelection::kUpperLimit);
  // etaCascadeMax = getMinimalSelection(femtoDreamCascadeSelection::kCascadeetaMax,
  //                                     femtoDreamSelection::kAbsUpperLimit);
  // DCAV0DaughMax = getMinimalSelection(femtoDreamCascadeSelection::kCascadeV0DCADaughMax,
  //                                     femtoDreamSelection::kUpperLimit);
  // CPAV0Min = getMinimalSelection(femtoDreamCascadeSelection::kCascadeV0CPAMin,
  //                                femtoDreamSelection::kLowerLimit);
  // TranRadV0Min = getMinimalSelection(femtoDreamCascadeSelection::kCascadeV0TranRadMin,
  //                                    femtoDreamSelection::kLowerLimit);
  // TranRadV0Max = getMinimalSelection(femtoDreamCascadeSelection::kCascadeV0TranRadMax,
  //                                    femtoDreamSelection::kUpperLimit);
  // V0DecVtxMax = getMinimalSelection(femtoDreamCascadeSelection::kCascadeV0DecVtxMax,
  //                                   femtoDreamSelection::kAbsUpperLimit);
  // DCACascadeDaughMax = getMinimalSelection(femtoDreamCascadeSelection::kCascadeDCADaughMax,
  //                                          femtoDreamSelection::kUpperLimit);
  // CPACascadeMin = getMinimalSelection(femtoDreamCascadeSelection::kCascadeCPAMin,
  //                                     femtoDreamSelection::kLowerLimit);
  // TranRadCascadeMin = getMinimalSelection(femtoDreamCascadeSelection::kCascadeTranRadMin,
  //                                         femtoDreamSelection::kLowerLimit);
  // TranRadCascadeMax = getMinimalSelection(femtoDreamCascadeSelection::kCascadeTranRadMax,
  //                                         femtoDreamSelection::kUpperLimit);
  // DecVtxMax = getMinimalSelection(femtoDreamCascadeSelection::kCascadeDecVtxMax,
  //                                 femtoDreamSelection::kAbsUpperLimit);
  // DCAPosToPV = getMinimalSelection(femtoDreamCascadeSelection::kCascadeDCAPosToPV,
  //                                  femtoDreamSelection::kLowerLimit);
  // DCANegToPV = getMinimalSelection(femtoDreamCascadeSelection::kCascadeDCANegToPV,
  //                                  femtoDreamSelection::kLowerLimit);
  // DCABachToPV = getMinimalSelection(femtoDreamCascadeSelection::kCascadeDCABachToPV,
  //                                   femtoDreamSelection::kLowerLimit);
  // DCAV0ToPV = getMinimalSelection(femtoDreamCascadeSelection::kCascadeDCAV0ToPV,
  //                                 femtoDreamSelection::kLowerLimit);
  // fV0InvMassLowLimit = getMinimalSelection(femtoDreamCascadeSelection::kCascadeV0MassMin,
  //                                          femtoDreamSelection::kLowerLimit);
  // fV0InvMassUpLimit = getMinimalSelection(femtoDreamCascadeSelection::kCascadeV0MassMax,
  //                                           femtoDreamSelection::kUpperLimit);
  isCascOmega = isSelectCascOmega;
}

template <typename Col, typename Casc, typename Track>
bool FemtoDreamCascadeSelection::isSelectedMinimal(Col const& col, Casc const& cascade, Track const& posTrack, Track const& negTrack, Track const& bachTrack)
{
  const auto signPos = posTrack.sign();
  const auto signNeg = negTrack.sign();

  if (signPos < 0 || signNeg > 0) {
    LOG(warn) << "Something wrong in isSelectedMinimal";
    LOG(warn) << "ERROR - Wrong sign for V0 daughters";
  }

  const std::vector<float> decVtx = {cascade.x(), cascade.y(), cascade.z()};

  //const float cpav0 = cascade.v0cosPA(col.posX(), col.posY(), col.posZ());
  //const float cpaCasc = cascade.casccosPA(col.posX(), col.posY(), col.posZ());
  //const float dcav0topv = cascade.dcav0topv(col.posX(), col.posY(), col.posZ());
  // const float invMassLambda = cascade.mLambda();
  //const float invMass = isCascOmega ? cascade.mOmega() : cascade.mXi();
  const float invMass = cascade.mXi();
  /*
  if (invMassLambda < fV0InvMassLowLimit || invMassLambda > fV0InvMassUpLimit) {
    return false;
  }
  */
  if (invMass < fInvMassLowLimit || invMass > fInvMassUpLimit) {
    return false;
  }
  /* 
  if (fRejectCompetingMass) {
    const float invMassCompeting = isCascOmega ? cascade.mXi() : cascade.mOmega();
    if (invMassCompeting > fInvMassCompetingLowLimit &&
        invMassCompeting < fInvMassCompetingUpLimit) {
      return false;
    }
  }
  */
  if (nPtCascadeMinSel > 0 && cascade.pt() < pTCascadeMin) {
    return false;
  }
  if (nPtCascadeMaxSel > 0 && cascade.pt() > pTCascadeMax) {
    return false;
  }
  /*
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
  */

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
  LOGF(info, "GG CascadeSelection: A Xi is selected!"); //REMOVE COMMENT
  return true;
}

template <typename Col, typename Casc, typename Track>
void FemtoDreamCascadeSelection::fillCascadeQA(Col const& col, Casc const& cascade, Track const& posTrack, Track const& negTrack)
{
  const auto signPos = posTrack.sign();
  const auto signNeg = negTrack.sign();
  if (signPos < 0 || signNeg > 0) {
    LOG(warn) << "Something wrong in isSelectedMinimal";
    LOG(warn) << "ERROR - Wrong sign for V0 daughters";
  }

  // const std::vector<float> decVtx = {cascade.x(), cascade.y(), cascade.z()};
  // const float cpav0 = cascade.v0cosPA(col.posX(), col.posY(), col.posZ());
  // const float cpaCasc = cascade.casccosPA(col.posX(), col.posY(), col.posZ());
  // const float dcav0topv = cascade.dcav0topv(col.posX(), col.posY(), col.posZ());

  const float invMassLambda = cascade.mLambda();
  const float invMass = isCascOmega ? cascade.mOmega() : cascade.mXi();

  // mQAHistogramRegistry->fill(HIST("CascadeQA/hInvMassV0Cut"), invMassLambda);
  mQAHistogramRegistry->fill(HIST("CascadeQA/hInvMassCascadeCut"), invMass);
  mQAHistogramRegistry->fill(HIST("CascadeQA/hCascadePt"), cascade.pt());
  // mQAHistogramRegistry->fill(HIST("CascadeQA/hCascadeEta"), cascade.eta());
  // mQAHistogramRegistry->fill(HIST("CascadeQA/hCascadePhi"), cascade.phi());
  // mQAHistogramRegistry->fill(HIST("CascadeQA/hDCAV0Daugh"), cascade.dcaV0daughters());
  // mQAHistogramRegistry->fill(HIST("CascadeQA/hV0CPA"), cpav0);
  // mQAHistogramRegistry->fill(HIST("CascadeQA/hV0TranRad"), cascade.v0radius());
  // mQAHistogramRegistry->fill(HIST("CascadeQA/hCascadeCPA"), cpaCasc);
  // mQAHistogramRegistry->fill(HIST("CascadeQA/hDCACascadeDaugh"), cascade.dcacascdaughters());
  // mQAHistogramRegistry->fill(HIST("CascadeQA/hCascadeTranRad"), cascade.cascradius());
  // mQAHistogramRegistry->fill(HIST("CascadeQA/hDCAPosToPV"), cascade.dcapostopv());
  // mQAHistogramRegistry->fill(HIST("CascadeQA/hDCANegToPV"), cascade.dcanegtopv());
  // mQAHistogramRegistry->fill(HIST("CascadeQA/hDCABachToPV"), cascade.dcabachtopv());
  // mQAHistogramRegistry->fill(HIST("CascadeQA/hDCAV0ToPV"), dcav0topv);

  // is this necessary
  /*
  bool write = true;
  for (size_t i = 0; i < decVtx.size(); i++) {
    write = write && (decVtx.at(i) < DecVtxMax);
  }
  if (write) {
    mQAHistogramRegistry->fill(HIST("CAscadeQA/hInvMassCascadeDecVtxMax"),
                             cascade.mXi());
  }
  */
}

template <typename cutContainerType, typename Col, typename Casc, typename Track>
std::array<cutContainerType, 8> FemtoDreamCascadeSelection::getCutContainer(Col const& /*col*/, Casc const& casc, Track const& posTrack, Track const& negTrack, Track const& bachTrack)
{
  //Cut bit shenanigan
  auto outputPosTrack = PosDaughTrack.getCutContainer<cutContainerType>(posTrack, casc.positivept(), casc.positiveeta(), casc.dcapostopv());
  auto outputNegTrack = NegDaughTrack.getCutContainer<cutContainerType>(negTrack, casc.negativept(), casc.negativeeta(), casc.dcanegtopv());
  auto outputBachTrack = BachTrack.getCutContainer<cutContainerType>(bachTrack, casc.bachelorpt(), casc.bacheloreta(), casc.dcabachtopv());
  cutContainerType output = 0;
  size_t counter = 0;

  //auto xiMassNominal = o2::constants::physics::MassXiMinus;
  //auto xiMassHypothesis = casc.mXi();
  //auto antiXiMassHypothesis = casc.mAntiXi();
  //auto diffXi = abs(xiMassNominal - xiMassHypothesis);
  //auto diffAntiXi = abs(xiMassNominal - antiXiMassHypothesis);

  float sign = 0.;
  int nSigmaPIDMax = PosDaughTrack.getSigmaPIDMax();
  
  auto nSigmaPrNeg = negTrack.tpcNSigmaPr();
  auto nSigmaPiPos = posTrack.tpcNSigmaPi();
  auto nSigmaPiNeg = negTrack.tpcNSigmaPi();
  auto nSigmaPrPos = posTrack.tpcNSigmaPr();
  float nSigmaPIDOffsetTPC = 0.;

  //TODO: improve the selection of the Xi candidates (now I select only Xi/ antiXi based on the daughter Lambda)
  if (abs(nSigmaPrNeg - nSigmaPIDOffsetTPC) < nSigmaPIDMax && abs(nSigmaPiPos - nSigmaPIDOffsetTPC) < nSigmaPIDMax) {
    sign = -1.;
  } else if (abs(nSigmaPrPos - nSigmaPIDOffsetTPC) < nSigmaPIDMax && abs(nSigmaPiNeg - nSigmaPIDOffsetTPC) < nSigmaPIDMax) {
    sign = 1.;
  }

  const auto pT = casc.pt();

  float observable = 0.;
  for (auto& sel : mSelections) {
    const auto selVariable = sel.getSelectionVariable();
    switch (selVariable) {
      case (femtoDreamCascadeSelection::kCascadeSign):
        observable = sign;
        break;
      case (femtoDreamCascadeSelection::kCascadepTMin):
        observable = pT;
        break;
      case (femtoDreamCascadeSelection::kCascadepTMax):
        observable = pT;
        break;
    } //switch
    sel.checkSelectionSetBit(observable, output, counter, nullptr);
  }
  return {
    output,
    outputPosTrack.at(femtoDreamTrackSelection::TrackContainerPosition::kCuts),
    outputPosTrack.at(femtoDreamTrackSelection::TrackContainerPosition::kPID),
    outputNegTrack.at(femtoDreamTrackSelection::TrackContainerPosition::kCuts),
    outputNegTrack.at(femtoDreamTrackSelection::TrackContainerPosition::kPID),
    outputBachTrack.at(femtoDreamTrackSelection::TrackContainerPosition::kCuts),
    outputBachTrack.at(femtoDreamTrackSelection::TrackContainerPosition::kPID)};
}

template <o2::aod::femtodreamparticle::ParticleType part, typename Col, typename Casc, typename Track>
void FemtoDreamCascadeSelection::fillQA(Col const& /*col*/, Casc const& casc, Track const& posTrack, Track const& negTrack, Track const& bachTrack)
{

  if (mQAHistogramRegistry) {
    mQAHistogramRegistry->fill(
      HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) +
        HIST("/hPt"),
      casc.pt());
    mQAHistogramRegistry->fill(
      HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) +
        HIST("/hEta"),
      casc.eta());
    mQAHistogramRegistry->fill(
      HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) +
        HIST("/hPhi"),
      casc.phi());
  }
  PosDaughTrack.fillQA<aod::femtodreamparticle::ParticleType::kV0Child,
                       aod::femtodreamparticle::TrackType::kPosChild>(posTrack);
  NegDaughTrack.fillQA<aod::femtodreamparticle::ParticleType::kV0Child,
                       aod::femtodreamparticle::TrackType::kNegChild>(negTrack);
  BachTrack.fillQA<aod::femtodreamparticle::ParticleType::kCascadeBachelor,
                   aod::femtodreamparticle::TrackType::kBachelor>(bachTrack);
}

} // namespace o2::analysis::femtoDream

#endif // PWGCF_FEMTODREAM_CORE_FEMTODREAMCASCADESELECTION_H_
