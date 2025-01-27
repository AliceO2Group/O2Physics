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
/// \author Georgios Mantzaridis, TU München, georgios.mantzaridis@tum.de

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
  kCascadePtMin,
  kCascadePtMax,
  kCascadeEtaMax,
  kCascadeDCADaughMax,
  kCascadeCPAMin,
  kCascadeTranRadMin,
  kCascadeTranRadMax,
  kCascadeDecVtxMax,
  kCascadeV0DCADaughMax,
  kCascadeV0CPAMin,
  kCascadeV0TranRadMin,
  kCascadeV0TranRadMax,
  kCascadeV0DCAtoPVMin,
  kCascadeV0DCAtoPVMax

  // kNcascadeSelection
  // kCascadeV0MassMin,
  // kCascadeV0MassMax
};
/*
kCascadeDCAPosToPV,
kCascadeDCANegToPV,
kCascadeDCABachToPV,
*/

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
    : nCascadePtMin(0),
      nCascadePtMax(0),
      nCascadeEtaMax(0),
      nCascadeDCADaughMax(0),
      nCascadeCPAMin(0),
      nCascadeTranRadMin(0),
      nCascadeTranRadMax(0),
      nCascadeDecVtxMax(0),
      /*
      nCascadeDCAPosToPV(0),
      nCascadeDCANegToPV(0),
      nCascadeDCABachToPV(0),
      */
      nCascadeV0DCADaughMax(0),
      nCascadeV0CPAMin(0),
      nCascadeV0TranRadMin(0),
      nCascadeV0TranRadMax(0),
      nCascadeV0DCAToPVMin(0),
      nCascadeV0DCAToPVMax(0),

      fCascadePtMin(9999999),
      fCascadePtMax(-9999999),
      fCascadeEtaMax(-9999999),
      fCascadeDCADaughMax(-9999999),
      fCascadeCPAMin(9999999),
      fCascadeTranRadMin(9999999),
      fCascadeTranRadMax(-9999999),
      fCascadeDecVtxMax(-9999999),
      /*
      fCascadeDCAPosToPV(9999999),
      fCascadeDCANegToPV(9999999),
      fCascadeDCABachToPV(9999999),
      */
      fCascadeV0DCADaughMax(-9999999),
      fCascadeV0CPAMin(9999999),
      fCascadeV0TranRadMin(9999999),
      fCascadeV0TranRadMax(-9999999),
      fCascadeV0DCAToPVMin(9999999),
      fCascadeV0DCAToPVMax(-9999999),

      fV0InvMassLowLimit(1.05),
      fV0InvMassUpLimit(1.3),
      fInvMassLowLimit(1.25),
      fInvMassUpLimit(1.4),
      fRejectCompetingMass(false),
      fInvMassCompetingLowLimit(1.5),
      fInvMassCompetingUpLimit(2.0),
      isCascOmega(false)
  /*,nSigmaPIDOffsetTPC(0.)*/
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

  template <int cutstage = 0, o2::aod::femtodreamparticle::ParticleType part, typename Col, typename Casc, typename Track>
  void fillQA(Col const& col, Casc const& cascade, Track const& posTrack, Track const& negTrack, Track const& bachTrack);

  // template <typename cutContainerType, typename Col, typename Casc, typename V0, typename Track>
  // std::array<cutContainerType, 8> getCutContainer(Col const& col, Casc const& casc, V0 const& v0Daugh, Track const& posTrack, Track const& negTrack, Track const& bachTrack);
  template <typename cutContainerType, typename Col, typename Casc, typename Track>
  std::array<cutContainerType, 8> getCutContainer(Col const& col, Casc const& casc, Track const& posTrack, Track const& negTrack, Track const& bachTrack);

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
      BachDaugTrack.setSelection(selVal, selVar, selType);
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
      BachDaugTrack.setPIDSpecies(pids);
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

  void setV0InvMassLimits(float lowLimit, float upLimit)
  {
    fV0InvMassLowLimit = lowLimit;
    fV0InvMassUpLimit = upLimit;
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
  int nCascadePtMin;
  int nCascadePtMax;
  int nCascadeEtaMax;
  int nCascadeDCADaughMax;
  int nCascadeCPAMin;
  int nCascadeTranRadMin;
  int nCascadeTranRadMax;
  int nCascadeDecVtxMax;

  int nCascadeV0DCADaughMax;
  int nCascadeV0CPAMin;
  int nCascadeV0TranRadMin;
  int nCascadeV0TranRadMax;
  int nCascadeV0DCAToPVMin;
  int nCascadeV0DCAToPVMax;

  float fCascadePtMin;
  float fCascadePtMax;
  float fCascadeEtaMax;
  float fCascadeDCADaughMax;
  float fCascadeCPAMin;
  float fCascadeTranRadMin;
  float fCascadeTranRadMax;
  float fCascadeDecVtxMax;

  float fCascadeV0DCADaughMax;
  float fCascadeV0CPAMin;
  float fCascadeV0TranRadMin;
  float fCascadeV0TranRadMax;
  float fCascadeV0DCAToPVMin;
  float fCascadeV0DCAToPVMax;

  float fV0InvMassLowLimit;
  float fV0InvMassUpLimit;

  float fInvMassLowLimit;
  float fInvMassUpLimit;

  float fRejectCompetingMass;
  float fInvMassCompetingLowLimit;
  float fInvMassCompetingUpLimit;

  bool isCascOmega;

  // float nSigmaPIDOffsetTPC;

  FemtoDreamTrackSelection PosDaughTrack;
  FemtoDreamTrackSelection NegDaughTrack;
  FemtoDreamTrackSelection BachDaugTrack;

  static constexpr int kNcascadeSelection = 16;

  static constexpr std::string_view mSelectionNames[kNcascadeSelection] = {
    "Sign", "PtMin", "PtMax", "EtaMax", "DCADaughMax", "CPAMin", "TranRadMin", "TranRadMax", "DecVtxMax", // Cascade Selections
    "DCAv0daughMax", "v0CPAMin", "v0TranRadMin", "v0TranRadMax", "V0DCAToPVMin", "V0DCAToPVMax"};         // CascadeV0 selections

  static constexpr femtoDreamSelection::SelectionType
    mSelectionTypes[kNcascadeSelection]{
      femtoDreamSelection::kEqual,      // sign
      femtoDreamSelection::kLowerLimit, // pt min
      femtoDreamSelection::kUpperLimit, // pt max
      femtoDreamSelection::kUpperLimit, // eta max
      femtoDreamSelection::kUpperLimit, // DCA cascade daughters max
      femtoDreamSelection::kLowerLimit, // cascade cos PA min
      femtoDreamSelection::kLowerLimit, // cascade tran rad min
      femtoDreamSelection::kUpperLimit, // cascade tran rad max
      femtoDreamSelection::kUpperLimit, // cascade maximum distance of decay vertex to PV
      femtoDreamSelection::kUpperLimit, // v0 daughters DCA max
      femtoDreamSelection::kLowerLimit, // v0 cos PA min
      femtoDreamSelection::kLowerLimit, // v0 tran rad min
      femtoDreamSelection::kUpperLimit, // v0 tran rad max
      femtoDreamSelection::kLowerLimit, // v0 minimum distance of decay vertex to PV
      femtoDreamSelection::kUpperLimit  // v0 maximum distance of decay vertex to PV

    }; ///< Map to match a variable with
       ///< its type

  static constexpr std::string_view mSelectionHelper[kNcascadeSelection] = {
    "Cascade particle sign (+1 or -1)",
    "Minimum pT (GeV/c)",
    "Maximum pT (GeV/c)",
    "Maximum |Eta|",
    "Maximum DCA between cascade daughters (cm)",
    "Minimum Cosine of Pointing Angle for cascade",
    "Minimum cascade transverse radius (cm)",
    "Maximum cascade transverse radius (cm)",
    "Maximum distance of cascade from primary vertex",
    "Maximum DCA between v0 daughters (cm)",
    "Minimum Cosine of Pointing Angle for v0",
    "Minimum v0 transverse radius (cm)",
    "Maximum v0 transverse radius (cm)",
    "Minimum distance of v0 from primary vertex",
    "Maximum distance of v0 from primary vertex"

    //"Minimum V0 mass",
    //"Maximum V0 mass"
  }; ///< Helper information for the
     ///< different selections

  static constexpr int kNcutStages = 2;
  static constexpr std::string_view mCutStage[kNcutStages] = {"BeforeSel", "AfterSel"};
}; // namespace femtoDream

template <o2::aod::femtodreamparticle::ParticleType part, o2::aod::femtodreamparticle::ParticleType daugh, o2::aod::femtodreamparticle::ParticleType bach, typename cutContainerType>
void FemtoDreamCascadeSelection::init(HistogramRegistry* QAregistry, HistogramRegistry* Registry, bool isSelectCascOmega)
{

  if (QAregistry && Registry) {
    mHistogramRegistry = Registry;
    mQAHistogramRegistry = QAregistry;
    // fillSelectionHistogram<part>();  // cascade
    // fillSelectionHistogram<daugh>(); // pos, neg
    // fillSelectionHistogram<bach>();  // bach

    AxisSpec ptAxis = {100, 0.0f, 10.0f, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec etaAxis = {100, -2.0f, 2.0f, "#it{#eta}"};
    AxisSpec phiAxis = {100, 0.0f, 6.0f, "#it{#phi}"};
    AxisSpec DCADaughAxis = {1000, 0.0f, 2.0f, "DCA (cm)"};
    AxisSpec CPAAxis = {1000, 0.95f, 1.0f, "#it{cos #theta_{p}}"};
    AxisSpec tranRadAxis = {1000, 0.0f, 100.0f, "#it{r}_{xy} (cm)"};
    AxisSpec decVtxAxis = {2000, 0, 200, "#it{Vtx}_{z} (cm)"};
    AxisSpec massAxisCascade = {2200, 1.25f, 1.8f, "m_{#Cascade} (GeV/#it{c}^{2})"};

    AxisSpec DCAToPVAxis = {1000, -10.0f, 10.0f, "DCA to PV (cm)"};

    AxisSpec massAxisV0 = {600, 0.0f, 3.0f, "m_{#V0} (GeV/#it{c}^{2})"};

    /// \todo this should be an automatic check in the parent class, and the
    /// return type should be templated
    size_t nSelections = getNSelections();
    if (nSelections > 8 * sizeof(cutContainerType)) {
      LOGF(info, "Number of selections %i", nSelections);
      LOG(fatal) << "FemtoDreamCascadeCuts: Number of selections to large for your "
                    "container - quitting!";
    }

    std::string folderName = static_cast<std::string>(o2::aod::femtodreamparticle::ParticleTypeName[part]);
    for (int istage = 0; istage < kNcutStages; istage++) {
      mQAHistogramRegistry->add((folderName + "/" + static_cast<std::string>(mCutStage[istage]) + "/hSign").c_str(), "; Sign of the Cascade ; Entries", kTH1I, {{3, -1, 2}});
      mQAHistogramRegistry->add((folderName + "/" + static_cast<std::string>(mCutStage[istage]) + "/hPt").c_str(), "; #it{p}_{T} (GeV/#it{c}); Entries", kTH1F, {{1000, 0, 10}});
      mQAHistogramRegistry->add((folderName + "/" + static_cast<std::string>(mCutStage[istage]) + "/hEta").c_str(), "; #eta; Entries", kTH1F, {{1000, -1, 1}});
      mQAHistogramRegistry->add((folderName + "/" + static_cast<std::string>(mCutStage[istage]) + "/hPhi").c_str(), "; #phi; Entries", kTH1F, {{1000, 0, 2. * M_PI}});
      mQAHistogramRegistry->add((folderName + "/" + static_cast<std::string>(mCutStage[istage]) + "/hDCADaugh").c_str(), "; daughters DCA; Entries", kTH1F, {DCADaughAxis});
      mQAHistogramRegistry->add((folderName + "/" + static_cast<std::string>(mCutStage[istage]) + "/hCPA").c_str(), "; Cos PA; Entries", kTH1F, {CPAAxis});
      mQAHistogramRegistry->add((folderName + "/" + static_cast<std::string>(mCutStage[istage]) + "/hTranRad").c_str(), "; Transverse Radius; Entries", kTH1F, {tranRadAxis});
      mQAHistogramRegistry->add((folderName + "/" + static_cast<std::string>(mCutStage[istage]) + "/hDecVtxX").c_str(), "; Decay vertex x position; Entries", kTH1F, {tranRadAxis});
      mQAHistogramRegistry->add((folderName + "/" + static_cast<std::string>(mCutStage[istage]) + "/hDecVtxY").c_str(), "; Decay vertex y position; Entries", kTH1F, {tranRadAxis});
      mQAHistogramRegistry->add((folderName + "/" + static_cast<std::string>(mCutStage[istage]) + "/hDecVtxZ").c_str(), "; Decay vertex z position; Entries", kTH1F, {tranRadAxis});
      mQAHistogramRegistry->add((folderName + "/" + static_cast<std::string>(mCutStage[istage]) + "/hInvMass").c_str(), "; Invariant mass; Entries", kTH1F, {massAxisCascade});
      mQAHistogramRegistry->add((folderName + "/" + static_cast<std::string>(mCutStage[istage]) + "/hV0DCADaugh").c_str(), "; V0-daughters DCA; Entries", kTH1F, {DCADaughAxis});
      mQAHistogramRegistry->add((folderName + "/" + static_cast<std::string>(mCutStage[istage]) + "/hV0CPA").c_str(), "; V0 cos PA; Entries", kTH1F, {CPAAxis});
      mQAHistogramRegistry->add((folderName + "/" + static_cast<std::string>(mCutStage[istage]) + "/hV0TranRad").c_str(), "; V0 transverse radius; Entries", kTH1F, {tranRadAxis});
      mQAHistogramRegistry->add((folderName + "/" + static_cast<std::string>(mCutStage[istage]) + "/hV0DCAToPV").c_str(), "; DCA of the V0 to the PV; Entries", kTH1F, {DCAToPVAxis});
      mQAHistogramRegistry->add((folderName + "/" + static_cast<std::string>(mCutStage[istage]) + "/hV0InvMass").c_str(), "; Invariant mass Cascade V0; Entries", kTH1F, {massAxisV0});
    }

    PosDaughTrack.init<aod::femtodreamparticle::ParticleType::kCascadeV0Child,
                       aod::femtodreamparticle::TrackType::kPosChild,
                       aod::femtodreamparticle::cutContainerType>(mQAHistogramRegistry, mHistogramRegistry);

    NegDaughTrack.init<aod::femtodreamparticle::ParticleType::kCascadeV0Child,
                       aod::femtodreamparticle::TrackType::kNegChild,
                       aod::femtodreamparticle::cutContainerType>(mQAHistogramRegistry, mHistogramRegistry);

    BachDaugTrack.init<aod::femtodreamparticle::ParticleType::kCascadeBachelor,
                       aod::femtodreamparticle::TrackType::kBachelor,
                       aod::femtodreamparticle::cutContainerType>(mQAHistogramRegistry, mHistogramRegistry);
  }

  /// check whether the most open cuts are fulfilled - most of this should have
  /// already be done by the filters
  nCascadePtMin = getNSelections(femtoDreamCascadeSelection::kCascadePtMin);
  nCascadePtMax = getNSelections(femtoDreamCascadeSelection::kCascadePtMax);
  nCascadeEtaMax = getNSelections(femtoDreamCascadeSelection::kCascadeEtaMax);
  nCascadeDCADaughMax = getNSelections(femtoDreamCascadeSelection::kCascadeDCADaughMax);
  nCascadeCPAMin = getNSelections(femtoDreamCascadeSelection::kCascadeCPAMin);
  nCascadeTranRadMin = getNSelections(femtoDreamCascadeSelection::kCascadeTranRadMin);
  nCascadeTranRadMax = getNSelections(femtoDreamCascadeSelection::kCascadeTranRadMax);
  nCascadeDecVtxMax = getNSelections(femtoDreamCascadeSelection::kCascadeDecVtxMax);

  nCascadeV0DCADaughMax = getNSelections(femtoDreamCascadeSelection::kCascadeV0DCADaughMax);
  nCascadeV0CPAMin = getNSelections(femtoDreamCascadeSelection::kCascadeV0CPAMin);
  nCascadeV0TranRadMin = getNSelections(femtoDreamCascadeSelection::kCascadeV0TranRadMin);
  nCascadeV0TranRadMax = getNSelections(femtoDreamCascadeSelection::kCascadeV0TranRadMax);

  nCascadeV0DCAToPVMin = getNSelections(femtoDreamCascadeSelection::kCascadeV0DCAtoPVMin);
  nCascadeV0DCAToPVMax = getNSelections(femtoDreamCascadeSelection::kCascadeV0DCAtoPVMax);

  fCascadePtMin = getMinimalSelection(femtoDreamCascadeSelection::kCascadePtMin,
                                      femtoDreamSelection::kLowerLimit);
  fCascadePtMax = getMinimalSelection(femtoDreamCascadeSelection::kCascadePtMax,
                                      femtoDreamSelection::kUpperLimit);
  fCascadeEtaMax = getMinimalSelection(femtoDreamCascadeSelection::kCascadeEtaMax,
                                       femtoDreamSelection::kAbsUpperLimit);
  fCascadeDCADaughMax = getMinimalSelection(femtoDreamCascadeSelection::kCascadeDCADaughMax,
                                            femtoDreamSelection::kUpperLimit);
  fCascadeCPAMin = getMinimalSelection(femtoDreamCascadeSelection::kCascadeCPAMin,
                                       femtoDreamSelection::kLowerLimit);
  fCascadeTranRadMin = getMinimalSelection(femtoDreamCascadeSelection::kCascadeTranRadMin,
                                           femtoDreamSelection::kLowerLimit);
  fCascadeTranRadMax = getMinimalSelection(femtoDreamCascadeSelection::kCascadeTranRadMax,
                                           femtoDreamSelection::kUpperLimit);
  fCascadeDecVtxMax = getMinimalSelection(femtoDreamCascadeSelection::kCascadeDecVtxMax,
                                          femtoDreamSelection::kAbsUpperLimit);
  fCascadeV0DCADaughMax = getMinimalSelection(femtoDreamCascadeSelection::kCascadeV0DCADaughMax,
                                              femtoDreamSelection::kUpperLimit);
  fCascadeV0CPAMin = getMinimalSelection(femtoDreamCascadeSelection::kCascadeV0CPAMin,
                                         femtoDreamSelection::kLowerLimit);
  fCascadeV0TranRadMin = getMinimalSelection(femtoDreamCascadeSelection::kCascadeV0TranRadMin,
                                             femtoDreamSelection::kLowerLimit);
  fCascadeV0TranRadMax = getMinimalSelection(femtoDreamCascadeSelection::kCascadeV0TranRadMax,
                                             femtoDreamSelection::kUpperLimit);
  fCascadeV0DCAToPVMin = getMinimalSelection(femtoDreamCascadeSelection::kCascadeV0DCAtoPVMin,
                                             femtoDreamSelection::kLowerLimit);
  fCascadeV0DCAToPVMax = getMinimalSelection(femtoDreamCascadeSelection::kCascadeV0DCAtoPVMax,
                                             femtoDreamSelection::kUpperLimit);
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

  const float cpav0 = cascade.v0cosPA(col.posX(), col.posY(), col.posZ());
  const float cpaCasc = cascade.casccosPA(col.posX(), col.posY(), col.posZ());
  const float v0dcatopv = cascade.dcav0topv(col.posX(), col.posY(), col.posZ());
  const float invMassLambda = cascade.mLambda();
  const float invMass = isCascOmega ? cascade.mOmega() : cascade.mXi();
  // const float invMass = cascade.mXi();

  // LOGF(info, "GG producer: Charge %i", cascade.sign());
  if (invMassLambda < fV0InvMassLowLimit || invMassLambda > fV0InvMassUpLimit) {
    return false;
  }

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

  if (nCascadePtMin > 0 && cascade.pt() < fCascadePtMin) {
    return false;
  }
  if (nCascadePtMax > 0 && cascade.pt() > fCascadePtMax) {
    return false;
  }
  if (nCascadeEtaMax > 0 && std::fabs(cascade.eta()) > fCascadeEtaMax) {
    return false;
  }
  if (nCascadeDCADaughMax > 0 && cascade.dcacascdaughters() > fCascadeDCADaughMax) {
    return false;
  }
  if (fCascadeCPAMin > 0 && cpaCasc < fCascadeCPAMin) {
    return false;
  }
  if (nCascadeTranRadMin > 0 && cascade.cascradius() < fCascadeTranRadMin) {
    return false;
  }
  if (nCascadeTranRadMax > 0 && cascade.cascradius() > fCascadeTranRadMax) {
    return false;
  }
  for (size_t i = 0; i < decVtx.size(); i++) {
    if (nCascadeDecVtxMax > 0 && decVtx.at(i) > fCascadeDecVtxMax) {
      return false;
    }
  }

  // v0 criteria
  if (nCascadeV0DCADaughMax > 0 && cascade.dcaV0daughters() > fCascadeV0DCADaughMax) {
    return false;
  }
  if (nCascadeV0CPAMin > 0 && cpav0 < fCascadeV0CPAMin) {
    return false;
  }
  if (nCascadeV0TranRadMin > 0 && cascade.v0radius() < fCascadeV0TranRadMin) {
    return false;
  }
  if (nCascadeV0TranRadMax > 0 && cascade.v0radius() > fCascadeV0TranRadMax) {
    return false;
  }
  if (nCascadeV0DCAToPVMin > 0 && std::fabs(v0dcatopv) < fCascadeV0DCAToPVMin) {
    return false;
  }
  if (nCascadeV0DCAToPVMax > 0 && std::fabs(v0dcatopv) > fCascadeV0DCAToPVMax) {
    return false;
  }

  // Chech the selection criteria for the tracks as well
  if (!PosDaughTrack.isSelectedMinimal(posTrack)) {
    return false;
  }
  if (!NegDaughTrack.isSelectedMinimal(negTrack)) {
    return false;
  }
  if (!BachDaugTrack.isSelectedMinimal(bachTrack)) {
    return false;
  }

  return true;
}

template <typename cutContainerType, typename Col, typename Casc, typename Track>
std::array<cutContainerType, 8> FemtoDreamCascadeSelection::getCutContainer(Col const& col, Casc const& casc, Track const& posTrack, Track const& negTrack, Track const& bachTrack)
{
  // Cut bit
  auto outputPosTrack = PosDaughTrack.getCutContainer<false, cutContainerType>(posTrack, casc.positivept(), casc.positiveeta(), casc.dcapostopv());
  auto outputNegTrack = NegDaughTrack.getCutContainer<false, cutContainerType>(negTrack, casc.negativept(), casc.negativeeta(), casc.dcanegtopv());
  auto outputBachTrack = BachDaugTrack.getCutContainer<false, cutContainerType>(bachTrack, casc.bachelorpt(), casc.bacheloreta(), casc.dcabachtopv());
  cutContainerType output = 0;
  size_t counter = 0;

  float sign = 0.;
  if (casc.sign() < 0) {
    sign = -1.;
  } else {
    sign = 1.;
  }

  const auto cpaCasc = casc.casccosPA(col.posX(), col.posY(), col.posZ());
  const std::vector<float> decVtx = {casc.x(), casc.y(), casc.z()};
  const auto cpav0 = casc.v0cosPA(col.posX(), col.posY(), col.posZ());
  const auto v0dcatopv = casc.dcav0topv(col.posX(), col.posY(), col.posZ());

  // LOGF(info, "GG producer: New dcatoPV: %f", dcav0topv);
  float observable = 0.;
  for (auto& sel : mSelections) {

    const auto selVariable = sel.getSelectionVariable();
    switch (selVariable) {
      case (femtoDreamCascadeSelection::kCascadeSign):
        observable = sign;
        break;
      case (femtoDreamCascadeSelection::kCascadePtMin):
        observable = casc.pt();
        break;
      case (femtoDreamCascadeSelection::kCascadePtMax):
        observable = casc.pt();
        break;
      case (femtoDreamCascadeSelection::kCascadeEtaMax):
        observable = casc.eta();
        break;
      case (femtoDreamCascadeSelection::kCascadeDCADaughMax):
        observable = casc.dcacascdaughters();
        break;
      case (femtoDreamCascadeSelection::kCascadeCPAMin):
        observable = cpaCasc;
        break;
      case (femtoDreamCascadeSelection::kCascadeTranRadMin):
        observable = casc.cascradius();
        break;
      case (femtoDreamCascadeSelection::kCascadeTranRadMax):
        observable = casc.cascradius();
        break;
      // kCascadeDecVtxMax is done above
      case (femtoDreamCascadeSelection::kCascadeDecVtxMax):
        for (size_t i = 0; i < decVtx.size(); ++i) {
          auto decVtxValue = decVtx.at(i);
          sel.checkSelectionSetBit(decVtxValue, output, counter, nullptr);
        }
        continue;
        break;

      case (femtoDreamCascadeSelection::kCascadeV0DCADaughMax):
        observable = casc.dcaV0daughters();
        break;
      case (femtoDreamCascadeSelection::kCascadeV0CPAMin):
        observable = cpav0;
        break;
      case (femtoDreamCascadeSelection::kCascadeV0TranRadMin):
        observable = casc.v0radius();
        break;
      case (femtoDreamCascadeSelection::kCascadeV0TranRadMax):
        observable = casc.v0radius();
        break;
      case (femtoDreamCascadeSelection::kCascadeV0DCAtoPVMin):
        observable = v0dcatopv;
        // LOGF(info, "==> Now it is: %f", dcav0topv);
        break;
      case (femtoDreamCascadeSelection::kCascadeV0DCAtoPVMax):
        observable = v0dcatopv;
        break;
    } // switch
    sel.checkSelectionSetBit(observable, output, counter, nullptr);
    //}
  } // for loop

  return {
    output,
    outputPosTrack.at(femtoDreamTrackSelection::TrackContainerPosition::kCuts),
    outputPosTrack.at(femtoDreamTrackSelection::TrackContainerPosition::kPID),
    outputNegTrack.at(femtoDreamTrackSelection::TrackContainerPosition::kCuts),
    outputNegTrack.at(femtoDreamTrackSelection::TrackContainerPosition::kPID),
    outputBachTrack.at(femtoDreamTrackSelection::TrackContainerPosition::kCuts),
    outputBachTrack.at(femtoDreamTrackSelection::TrackContainerPosition::kPID)};
}

template <int cutstage, o2::aod::femtodreamparticle::ParticleType part, typename Col, typename Casc, typename Track>
void FemtoDreamCascadeSelection::fillQA(Col const& col, Casc const& casc, Track const& posTrack, Track const& negTrack, Track const& bachTrack)
{

  const std::vector<float> decVtx = {casc.x(), casc.y(), casc.z()};
  const float cpaCasc = casc.casccosPA(col.posX(), col.posY(), col.posZ());
  const float cpav0 = casc.v0cosPA(col.posX(), col.posY(), col.posZ());
  const float v0dcatopv = casc.dcav0topv(col.posX(), col.posY(), col.posZ());
  if (mQAHistogramRegistry) {
    mQAHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) + HIST("/") + HIST(mCutStage[cutstage]) + HIST("/hSign"), casc.sign());
    mQAHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) + HIST("/") + HIST(mCutStage[cutstage]) + HIST("/hPt"), casc.pt());
    mQAHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) + HIST("/") + HIST(mCutStage[cutstage]) + HIST("/hEta"), casc.eta());
    mQAHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) + HIST("/") + HIST(mCutStage[cutstage]) + HIST("/hPhi"), casc.phi());
    mQAHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) + HIST("/") + HIST(mCutStage[cutstage]) + HIST("/hDCADaugh"), casc.dcacascdaughters());
    mQAHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) + HIST("/") + HIST(mCutStage[cutstage]) + HIST("/hCPA"), cpaCasc);
    mQAHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) + HIST("/") + HIST(mCutStage[cutstage]) + HIST("/hTranRad"), casc.cascradius());
    mQAHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) + HIST("/") + HIST(mCutStage[cutstage]) + HIST("/hDecVtxX"), decVtx.at(0));
    mQAHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) + HIST("/") + HIST(mCutStage[cutstage]) + HIST("/hDecVtxY"), decVtx.at(1));
    mQAHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) + HIST("/") + HIST(mCutStage[cutstage]) + HIST("/hDecVtxZ"), decVtx.at(2));
    mQAHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) + HIST("/") + HIST(mCutStage[cutstage]) + HIST("/hInvMass"), casc.mXi());
    mQAHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) + HIST("/") + HIST(mCutStage[cutstage]) + HIST("/hV0DCADaugh"), casc.dcaV0daughters());
    mQAHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) + HIST("/") + HIST(mCutStage[cutstage]) + HIST("/hV0CPA"), cpav0);
    mQAHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) + HIST("/") + HIST(mCutStage[cutstage]) + HIST("/hV0TranRad"), casc.v0radius());
    mQAHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) + HIST("/") + HIST(mCutStage[cutstage]) + HIST("/hV0DCAToPV"), v0dcatopv);
    mQAHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) + HIST("/") + HIST(mCutStage[cutstage]) + HIST("/hV0InvMass"), casc.mLambda());

    PosDaughTrack.fillQA<aod::femtodreamparticle::ParticleType::kCascadeV0Child,
                         aod::femtodreamparticle::TrackType::kPosChild, false, cutstage>(posTrack);
    NegDaughTrack.fillQA<aod::femtodreamparticle::ParticleType::kCascadeV0Child,
                         aod::femtodreamparticle::TrackType::kNegChild, false, cutstage>(negTrack);
    BachDaugTrack.fillQA<aod::femtodreamparticle::ParticleType::kCascadeBachelor,
                         aod::femtodreamparticle::TrackType::kBachelor, false, cutstage>(bachTrack);
  }
}

} // namespace o2::analysis::femtoDream

#endif // PWGCF_FEMTODREAM_CORE_FEMTODREAMCASCADESELECTION_H_
