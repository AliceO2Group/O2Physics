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
  kCascadeV0DCAtoPVMax,
  kCascadeV0MassMin,
  kCascadeV0MassMax
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

  template <typename Col, typename Casc, typename Track>
  void fillCascadeQA(Col const& col, Casc const& cascade, Track const& posTrack, Track const& negTrack);

  template <o2::aod::femtodreamparticle::ParticleType part, typename Col, typename Casc, typename Track>
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
  int nCascadePtMin;
  int nCascadePtMax;
  int nCascadeEtaMax;
  int nCascadeDCADaughMax;
  int nCascadeCPAMin;
  int nCascadeTranRadMin;
  int nCascadeTranRadMax;
  int nCascadeDecVtxMax;
  /* 
  int nCascadeDCAPosToPV;
  int nCascadeDCANegToPV;
  int nCascadeDCABachToPV;
  */
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
  /* 
  float fCascadeDCAPosToPV;
  float fCascadeDCANegToPV;
  float fCascadeDCABachToPV;
  */  
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
  FemtoDreamTrackSelection BachTrack;
  //FemtoDreamV0Selection V0DaughSel;

  static constexpr int kNcascadeSelection = 17; //TODO can I do less ?

  static constexpr std::string_view mSelectionNames[kNcascadeSelection] = {
    "Sign", "PtMin", "PtMax", "EtaMax", "DCAcascDaugh", "CPAMin", "TranRadMin", "TranRadMax", "DecVtxMax",                   //Cascade Selections
    "DCAv0daughMax", "v0CPAMin", "v0TranRadMin", "v0TranRadMax", "DCAV0ToPVMin", "DCAV0ToPVMax", "kV0MassMin", "V0MassMax"}; //CascadeV0 selections
    // "DCAPosToPV", "DCANegToPV", "DCABachToPV",                                                             //Cascade daughter track selections
    // }; //<< Name of the different selections

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
      femtoDreamSelection::kUpperLimit, // v0 maximum distance of decay vertex to PV
      femtoDreamSelection::kLowerLimit, // v0 mass min
      femtoDreamSelection::kUpperLimit // v0 mass max
    }; ///< Map to match a variable with
       ///< its type
      
      /*
      femtoDreamSelection::kLowerLimit, // DCA pos to PV max
      femtoDreamSelection::kLowerLimit, // DCA neg to PV max
      femtoDreamSelection::kLowerLimit, // DCA bach to PV max
      */
      

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
    "Maximum distance of v0 from primary vertex",
    "Minimum V0 mass",
    "Maximum V0 mass"
  }; ///< Helper information for the
     ///< different selections
    
    /*
    "Maximum DCA of positive track form primary vertex",
    "Maximum DCA of negative track form primary vertex",
    "Maximum DCA of bachelor track form primary vertex",

    }; ///< Helper information for the
                        ///< different selections             ///< different selections
    */ 
}; // namespace femtoDream

template <o2::aod::femtodreamparticle::ParticleType part, o2::aod::femtodreamparticle::ParticleType daugh, o2::aod::femtodreamparticle::ParticleType bach, typename cutContainerType>
void FemtoDreamCascadeSelection::init(HistogramRegistry* QAregistry, HistogramRegistry* Registry, bool isSelectCascOmega)
{

  if (QAregistry && Registry) {
    mHistogramRegistry = Registry;
    mQAHistogramRegistry = QAregistry;
    //fillSelectionHistogram<part>();  // cascade
    //fillSelectionHistogram<daugh>(); // pos, neg
    //fillSelectionHistogram<bach>();  // bach

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
                       aod::femtodreamparticle::cutContainerType>(mQAHistogramRegistry, mQAHistogramRegistry);
    
    NegDaughTrack.init<aod::femtodreamparticle::ParticleType::kCascadeV0Child,
                       aod::femtodreamparticle::TrackType::kNegChild,
                       aod::femtodreamparticle::cutContainerType>(mQAHistogramRegistry, mQAHistogramRegistry);
    
    BachTrack.init<aod::femtodreamparticle::ParticleType::kCascadeBachelor,
                   aod::femtodreamparticle::TrackType::kBachelor,
                   aod::femtodreamparticle::cutContainerType>(mQAHistogramRegistry, mQAHistogramRegistry);

    // Cascade (Xi, Omega)
    mQAHistogramRegistry->add("CascadeQA/hCascadePt", "pT distribution", kTH1F, {ptAxis});
    mQAHistogramRegistry->add("CascadeQA/hCascadeEta", "Eta distribution", kTH1F, {etaAxis});
    mQAHistogramRegistry->add("CascadeQA/hCascadePhi", "Phi distribution", kTH1F, {phiAxis});
    mQAHistogramRegistry->add("CascadeQA/hCascadeDCADaugh", "Cascade-daughters DCA", kTH1F, {DCADaughAxis});
    mQAHistogramRegistry->add("CascadeQA/hCascadeCPA", "Cos PA", kTH1F, {CPAAxis});
    mQAHistogramRegistry->add("CascadeQA/hCascadeTranRad", "Transverse radius", kTH1F, {tranRadAxis});
    mQAHistogramRegistry->add("CascadeQA/hCascadeDecVtxX", "Decay vertex x position", kTH1F, {decVtxAxis});
    mQAHistogramRegistry->add("CascadeQA/hCascadeDecVtxY", "Decay vertex y position", kTH1F, {decVtxAxis});
    mQAHistogramRegistry->add("CascadeQA/hCascadeDecVtxZ", "Decay vertex z position", kTH1F, {decVtxAxis});
    mQAHistogramRegistry->add("CascadeQA/hInvMassCascade", "Invariant mass Cascade", kTH1F, {massAxisCascade});
    // V0 (Lambda)
    mQAHistogramRegistry->add("CascadeQA/hCascadeV0DCADaugh", "V0-daughters DCA", kTH1F, {DCADaughAxis});
    mQAHistogramRegistry->add("CascadeQA/hCascadeV0CPA", "V0 cos PA", kTH1F, {CPAAxis});
    mQAHistogramRegistry->add("CascadeQA/hCascadeV0TranRad", "V0 transverse radius", kTH1F, {tranRadAxis});
    mQAHistogramRegistry->add("CascadeQA/hCascadeV0DCAToPV", "DCA of the V0 to the PV", kTH1F, {massAxisV0});
    mQAHistogramRegistry->add("CascadeQA/hInvMassV0", "Invariant mass Cascade V0", kTH1F, {massAxisV0});
    
    /* 
    // Dauchter Tracks
    mQAHistogramRegistry->add("CascadeQA/hCascadeDCAPosToPV", "Pos V0 daughter DCA to primary vertex", kTH1F, {DCAToPVAxis});
    mQAHistogramRegistry->add("CascadeQA/hCascadeDCANegToPV", "Neg V0 daughter DCA to primary vertex", kTH1F, {DCAToPVAxis});
    mQAHistogramRegistry->add("CascadeQA/hCascadeDCABachToPV", "Bachelor DCA to primary vertex", kTH1F, {DCAToPVAxis});
    */
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
  
  /*
  nCascadeDCAPosToPV = getNSelections(femtoDreamCascadeSelection::kCascadeDCAPosToPV);
  nCascadeDCANegToPV = getNSelections(femtoDreamCascadeSelection::kCascadeDCANegToPV);
  nCascadeDCABachToPV = getNSelections(femtoDreamCascadeSelection::kCascadeDCABachToPV);
  */ 
  //TODO v0mass??? 



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
  fV0InvMassLowLimit = getMinimalSelection(femtoDreamCascadeSelection::kCascadeV0MassMin,
                                           femtoDreamSelection::kLowerLimit);
  fV0InvMassUpLimit = getMinimalSelection(femtoDreamCascadeSelection::kCascadeV0MassMax,
                                          femtoDreamSelection::kUpperLimit);
   
  /*
  fCascadeDCAPosToPV = getMinimalSelection(femtoDreamCascadeSelection::kCascadeDCAPosToPV,
                                           femtoDreamSelection::kLowerLimit);
  fCascadeDCANegToPV = getMinimalSelection(femtoDreamCascadeSelection::kCascadeDCANegToPV,
                                           femtoDreamSelection::kLowerLimit);
  fCascadeDCABachToPV = getMinimalSelection(femtoDreamCascadeSelection::kCascadeDCABachToPV,
                                            femtoDreamSelection::kLowerLimit);
  */ 
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
  const float dcav0topv = cascade.dcav0topv(col.posX(), col.posY(), col.posZ());
  const float invMassLambda = cascade.mLambda();
  //const float invMass = isCascOmega ? cascade.mOmega() : cascade.mXi();
  const float invMass = cascade.mXi();
  
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
  if (nCascadeEtaMax > 0 && std::abs(cascade.eta()) > fCascadeEtaMax) {
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

  //v0 criteria
  if (nCascadeV0DCADaughMax > 0 && cascade.dcaV0daughters() > fCascadeV0DCADaughMax) {
    return false;
  }
  if (nCascadeV0CPAMin> 0 && cpav0 < fCascadeV0CPAMin) {
    return false;
  }
  if (nCascadeV0TranRadMin> 0 && cascade.v0radius() < fCascadeV0TranRadMin) {
    return false;
  }
  if (nCascadeV0TranRadMax> 0 && cascade.v0radius() < fCascadeV0TranRadMax) {
    return false;
  }
  if (nCascadeV0DCAToPVMin > 0 && abs(dcav0topv) < fCascadeV0DCAToPVMin) {
    return false;
  }
  if (nCascadeV0DCAToPVMax > 0 && abs(dcav0topv) < fCascadeV0DCAToPVMax) {
    return false;
  }
  
  
  /*
  if (nCascadeDCAPosToPV > 0 && abs(cascade.dcapostopv()) < fCascadeDCAPosToPV) {
    return false;
  }
  if (nCascadeDCANegToPV > 0 && abs(cascade.dcanegtopv()) < fCascadeDCANegToPV) {
    return false;
  }
  if (nCascadeDCABachToPV > 0 && abs(cascade.dcabachtopv()) < fCascadeDCABachToPV) {
    return false;
  }
  //Chech the selection criteria for the tracks as well (TODO) 
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
  //LOGF(info, "GG CascadeSelection: A Xi is selected!"); //REMOVE COMMENT
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

  const std::vector<float> decVtx = {cascade.x(), cascade.y(), cascade.z()};
  const float cpaCasc = cascade.casccosPA(col.posX(), col.posY(), col.posZ());
  const float invMass = isCascOmega ? cascade.mOmega() : cascade.mXi();
  
  const float cpav0 = cascade.v0cosPA(col.posX(), col.posY(), col.posZ());
  const float dcav0topv = cascade.dcav0topv(col.posX(), col.posY(), col.posZ());
  const float invMassLambda = cascade.mLambda();

  //Cascade
  mQAHistogramRegistry->fill(HIST("CascadeQA/hCascadePt"), cascade.pt());
  mQAHistogramRegistry->fill(HIST("CascadeQA/hCascadeEta"), cascade.eta());
  mQAHistogramRegistry->fill(HIST("CascadeQA/hCascadePhi"), cascade.phi());
  mQAHistogramRegistry->fill(HIST("CascadeQA/hCascadeDCADaugh"), cascade.dcacascdaughters());
  mQAHistogramRegistry->fill(HIST("CascadeQA/hCascadeCPA"), cpaCasc);
  mQAHistogramRegistry->fill(HIST("CascadeQA/hCascadeTranRad"), cascade.cascradius());
  mQAHistogramRegistry->fill(HIST("CascadeQA/hCascadeDecVtxX"), decVtx.at(0));
  mQAHistogramRegistry->fill(HIST("CascadeQA/hCascadeDecVtxY"), decVtx.at(1));
  mQAHistogramRegistry->fill(HIST("CascadeQA/hCascadeDecVtxZ"), decVtx.at(2));
  mQAHistogramRegistry->fill(HIST("CascadeQA/hInvMassCascade"), invMass);
  /*
  // Daughter Tracks
  mQAHistogramRegistry->fill(HIST("CascadeQA/hCascadeDCAPosToPV"), cascade.dcapostopv());
  mQAHistogramRegistry->fill(HIST("CascadeQA/hCascadeDCANegToPV"), cascade.dcanegtopv());
  mQAHistogramRegistry->fill(HIST("CascadeQA/hCascadeDCABachToPV"), cascade.dcabachtopv());
  //V0 (Lambda)
  mQAHistogramRegistry->fill(HIST("CascadeQA/hCascadeV0DCADaugh"), cascade.dcaV0daughters());
  mQAHistogramRegistry->fill(HIST("CascadeQA/hCascadeV0CPA"), cpav0);
  mQAHistogramRegistry->fill(HIST("CascadeQA/hCascadeV0TranRad"), cascade.v0radius());
  mQAHistogramRegistry->fill(HIST("CascadeQA/hCascadeV0DCAToPV"), dcav0topv);
  mQAHistogramRegistry->fill(HIST("CascadeQA/hInvMassV0"), invMassLambda);
  */ 

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

// template <typename cutContainerType, typename Col, typename Casc, typename V0, typename Track>
// std::array<cutContainerType, 8> FemtoDreamCascadeSelection::getCutContainer(Col const& col, Casc const& casc, V0 const& v0Daugh, Track const& posTrack, Track const& negTrack, Track const& bachTrack)
template <typename cutContainerType, typename Col, typename Casc, typename Track>
std::array<cutContainerType, 8> FemtoDreamCascadeSelection::getCutContainer(Col const& col, Casc const& casc, Track const& posTrack, Track const& negTrack, Track const& bachTrack)
{
  //Cut bit
  //auto outputV0Daugh = V0DaughSel.getCutContainer<cutContainerType>(v0Daugh, posTrack, negTrack);
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

  const auto cpaCasc = casc.casccosPA(col.posX(), col.posY(), col.posZ());
  const std::vector<float> decVtx = {casc.x(), casc.y(), casc.z()};
  const auto cpav0 = casc.v0cosPA(col.posX(), col.posY(), col.posZ());
  const auto dcav0topv = casc.dcav0topv(col.posX(), col.posY(), col.posZ());


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
        observable = dcav0topv; 
        break;
      case (femtoDreamCascadeSelection::kCascadeV0DCAtoPVMax):
        observable = dcav0topv; 
        break;
      case (femtoDreamCascadeSelection::kCascadeV0MassMin):
        observable = casc.mLambda(); 
        break;
      case (femtoDreamCascadeSelection::kCascadeV0MassMax):
        observable = casc.mLambda(); 
        break;
      
      /*
      case (femtoDreamCascadeSelection::kCascadeDCAPosToPV):
        observable = casc.dcapostopv();
        break;
      case (femtoDreamCascadeSelection::kCascadeDCANegToPV):
        observable = casc.dcanegtopv();
        break;
      case (femtoDreamCascadeSelection::kCascadeDCABachToPV):
        observable = casc.dcabachtopv();
        break;
      */ 

    } //switch
    sel.checkSelectionSetBit(observable, output, counter, nullptr);
    //}
  } //for loop

    /*
    outputV0Daugh.at(0), //daughter V0
    outputV0Daugh.at(1), //posDaughterCuts
    outputV0Daugh.at(2), //posDaughterPID
    outputV0Daugh.at(3), //negDaughterCuts
    outputV0Daugh.at(4), //negDaugherPID
    */
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
  //PosDaughTrack.fillQA<aod::femtodreamparticle::ParticleType::kCascadeV0Child,
  //                     aod::femtodreamparticle::TrackType::kPosChild>(posTrack);
  //NegDaughTrack.fillQA<aod::femtodreamparticle::ParticleType::kCascadeV0Child,
  //                     aod::femtodreamparticle::TrackType::kNegChild>(negTrack);
  BachTrack.fillQA<aod::femtodreamparticle::ParticleType::kCascadeBachelor,
                   aod::femtodreamparticle::TrackType::kBachelor>(bachTrack);
}

} // namespace o2::analysis::femtoDream

#endif // PWGCF_FEMTODREAM_CORE_FEMTODREAMCASCADESELECTION_H_
