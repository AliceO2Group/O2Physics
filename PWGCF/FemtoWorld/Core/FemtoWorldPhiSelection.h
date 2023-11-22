// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file FemtoWorldPhiSelection.h
/// \brief Definition of the FemtoWorldPhiSelection
/// \author Valentina Mantovani Sarti, TU München valentina.mantovani-sarti@tum.de
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de
/// \author Luca Barioglio, TU München, luca.barioglio@cern.ch
/// \author Zuzanna Chochulska, WUT Warsaw, zchochul@cern.ch

#ifndef PWGCF_FEMTOWORLD_CORE_FEMTOWORLDPHISELECTION_H_
#define PWGCF_FEMTOWORLD_CORE_FEMTOWORLDPHISELECTION_H_

#include <string>
#include <vector>

#include <TDatabasePDG.h> // FIXME

#include "PWGCF/FemtoWorld/Core/FemtoWorldObjectSelection.h"
#include "PWGCF/FemtoWorld/Core/FemtoWorldTrackSelection.h"
#include "PWGCF/FemtoWorld/Core/FemtoWorldSelection.h"

#include "ReconstructionDataFormats/PID.h"
#include "Common/Core/RecoDecay.h"
#include "Framework/HistogramRegistry.h"
#include "TLorentzVector.h"

using namespace o2::framework;

namespace o2::analysis::femtoWorld
{
namespace femtoWorldPhiSelection
{
/// The different selections this task is capable of doing
enum PhiSel { kPhiSign, ///< +1 particle, -1 antiparticle
              kpTPhiMin,
              kpTPhiMax,
              kDCAPhiDaughMax,
              kCPAPhiMin,
              kTranRadPhiMin,
              kTranRadPhiMax,
              kDecVtxMax };

enum ChildTrackType { kPosTrackPhi,
                      kNegTrackPhi };

enum PhiContainerPosition {
  kPhi,
  kPosCuts,
  kPosPID,
  kNegCuts,
  kNegPID,
}; /// Position in the full VO cut container

} // namespace femtoWorldPhiSelection

/// \class FemtoWorldPhiSelection
/// \brief Cut class to contain and execute all cuts applied to Phis
class FemtoWorldPhiSelection : public FemtoWorldObjectSelection<float, femtoWorldPhiSelection::PhiSel>
{
 public:
  FemtoWorldPhiSelection() : nPtPhiMinSel(0),
                             nPtPhiMaxSel(0),
                             nDCAPhiDaughMax(0),
                             nCPAPhiMin(0),
                             nTranRadPhiMin(0),
                             nTranRadPhiMax(0),
                             nDecVtxMax(0),
                             pTPhiMin(9999999.),
                             pTPhiMax(-9999999.),
                             DCAPhiDaughMax(-9999999.),
                             CPAPhiMin(9999999.),
                             TranRadPhiMin(9999999.),
                             TranRadPhiMax(-9999999.),
                             DecVtxMax(-9999999.),
                             fInvMassLowLimit(1.05),
                             fInvMassUpLimit(1.3),
                             fRejectKaon(false),
                             fInvMassKaonLowLimit(0.48),
                             fInvMassKaonUpLimit(0.515) {}
  /// Initializes histograms for the task
  template <o2::aod::femtoworldparticle::ParticleType part, o2::aod::femtoworldparticle::ParticleType daugh, typename cutContainerType>
  void init(HistogramRegistry* registry);

  template <o2::aod::femtoworldparticle::ParticleType part, o2::aod::femtoworldparticle::ParticleType daugh, typename cutContainerType>
  void initPhi(HistogramRegistry* registry);

  template <typename C, typename V, typename T>
  bool isSelectedMinimal(C const& col, V const& v0, T const& posTrack, T const& negTrack);

  template <typename C, typename V, typename T>
  void fillPhiQA(C const& col, V const& v0, T const& posTrack, T const& negTrack);

  template <typename C, typename V, typename T, typename P>
  void fillPhiQAMass(C const& col, V const& v0, T const& posTrack, T const& negTrack, P const& ConfInvMassLowLimit, P const& ConfInvMassUpLimit);

  /// \todo for the moment the PID of the tracks is factored out into a separate field, hence 5 values in total \\ASK: what does it mean?
  template <typename cutContainerType, typename C, typename T>
  std::array<cutContainerType, 5> getCutContainer(C const& col, T const& posTrack, T const& negTrack);

  template <o2::aod::femtoworldparticle::ParticleType part, o2::aod::femtoworldparticle::ParticleType daugh, typename C, typename V, typename T>
  void fillQA(C const& col, V const& v0, T const& posTrack, T const& negTrack);

  template <o2::aod::femtoworldparticle::ParticleType part, o2::aod::femtoworldparticle::ParticleType daugh, typename C, typename V, typename T>
  void fillQAPhi(C const& col, V const& v0, T const& posTrack, T const& negTrack);

  template <typename T1, typename T2>
  void setChildCuts(femtoWorldPhiSelection::ChildTrackType child, T1 selVal, T2 selVar, femtoWorldSelection::SelectionType selType)
  {
    if (child == femtoWorldPhiSelection::kPosTrackPhi) {
      PosDaughTrack.setSelection(selVal, selVar, selType);
    } else if (child == femtoWorldPhiSelection::kNegTrackPhi) {
      NegDaughTrack.setSelection(selVal, selVar, selType);
    }
  }
  template <typename T>
  void setChildPIDSpecies(femtoWorldPhiSelection::ChildTrackType child, T& pids)
  {
    if (child == femtoWorldPhiSelection::kPosTrackPhi) {
      PosDaughTrack.setPIDSpecies(pids);
    } else if (child == femtoWorldPhiSelection::kNegTrackPhi) {
      NegDaughTrack.setPIDSpecies(pids);
    }
  }

  /// Helper function to obtain the name of a given selection criterion for consistent naming of the configurables
  /// \param iSel Track selection variable to be examined
  /// \param prefix Additional prefix for the name of the configurable
  /// \param suffix Additional suffix for the name of the configurable
  static std::string getSelectionName(femtoWorldPhiSelection::PhiSel iSel, std::string_view prefix = "", std::string_view suffix = "")
  {
    std::string outString = static_cast<std::string>(prefix);
    outString += static_cast<std::string>(mSelectionNames[iSel]);
    outString += suffix;
    return outString;
  }

  /// Helper function to obtain the index of a given selection variable for consistent naming of the configurables
  /// \param obs Phi selection variable (together with prefix) got from file
  /// \param prefix Additional prefix for the output of the configurable
  static int findSelectionIndex(const std::string_view& obs, std::string_view prefix = "")
  {
    for (int index = 0; index < kNv0Selection; index++) {
      std::string comp = static_cast<std::string>(prefix) + static_cast<std::string>(mSelectionNames[index]);
      std::string_view cmp{comp};
      if (obs.compare(cmp) == 0)
        return index;
    }
    LOGF(info, "Variable %s not found", obs);
    return -1;
  }

  /// Helper function to obtain the type of a given selection variable for consistent naming of the configurables
  /// \param iSel Phi selection variable whose type is returned
  static femtoWorldSelection::SelectionType getSelectionType(femtoWorldPhiSelection::PhiSel iSel)
  {
    return mSelectionTypes[iSel];
  }

  /// Helper function to obtain the helper string of a given selection criterion for consistent description of the configurables
  /// \param iSel Track selection variable to be examined
  /// \param prefix Additional prefix for the output of the configurable
  static std::string getSelectionHelper(femtoWorldPhiSelection::PhiSel iSel, std::string_view prefix = "")
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

  void setChildRejectNotPropagatedTracks(femtoWorldPhiSelection::ChildTrackType child, bool reject)
  {
    if (child == femtoWorldPhiSelection::kPosTrackPhi) {
      PosDaughTrack.setRejectNotPropagatedTracks(reject);
    } else if (child == femtoWorldPhiSelection::kNegTrackPhi) {
      NegDaughTrack.setRejectNotPropagatedTracks(reject);
    }
  }

 private:
  int nPtPhiMinSel;
  int nPtPhiMaxSel;
  int nDCAPhiDaughMax;
  int nCPAPhiMin;
  int nTranRadPhiMin;
  int nTranRadPhiMax;
  int nDecVtxMax;
  float pTPhiMin;
  float pTPhiMax;
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

  FemtoWorldTrackSelection PosDaughTrack;
  FemtoWorldTrackSelection NegDaughTrack;

  static constexpr int kNv0Selection = 8;

  static constexpr std::string_view mSelectionNames[kNv0Selection] = {
    "Sign",
    "PtMin",
    "PtMax",
    "DCAdaughMax",
    "CPAMin",
    "TranRadMin",
    "TranRadMax",
    "DecVecMax"}; ///< Name of the different selections

  static constexpr femtoWorldSelection::SelectionType mSelectionTypes[kNv0Selection]{
    femtoWorldSelection::kEqual,
    femtoWorldSelection::kLowerLimit,
    femtoWorldSelection::kUpperLimit,
    femtoWorldSelection::kUpperLimit,
    femtoWorldSelection::kLowerLimit,
    femtoWorldSelection::kLowerLimit,
    femtoWorldSelection::kLowerLimit,
    femtoWorldSelection::kUpperLimit}; ///< Map to match a variable with its type

  static constexpr std::string_view mSelectionHelper[kNv0Selection] = {
    "+1 for lambda, -1 for antilambda",
    "Minimum pT (GeV/c)",
    "Maximum pT (GeV/c)",
    "Maximum DCA between daughters (cm)",
    "Minimum Cosine of Pointing Angle",
    "Minimum transverse radius (cm)",
    "Maximum transverse radius (cm)",
    "Maximum distance from primary vertex"}; ///< Helper information for the different selections

}; // namespace femtoWorld

template <o2::aod::femtoworldparticle::ParticleType part, o2::aod::femtoworldparticle::ParticleType daugh, typename cutContainerType>
void FemtoWorldPhiSelection::init(HistogramRegistry* registry)
{
  if (registry) {
    mHistogramRegistry = registry;
    fillSelectionHistogram<part>();
    fillSelectionHistogram<daugh>();

    AxisSpec massAxisPhi = {60000, 0.0f, 3.0f, "m_{#Phi} (GeV/#it{c}^{2})"};
    AxisSpec massAxisAntiPhi = {60000, 0.0f, 3.0f, "m_{#bar{#Phi}} (GeV/#it{c}^{2})"};

    /// \todo this should be an automatic check in the parent class, and the return type should be templated
    size_t nSelections = getNSelections();
    if (nSelections > 8 * sizeof(cutContainerType)) {
      LOG(fatal) << "FemtoWorldPhiCuts: Number of selections to large for your container - quitting!";
    }
    std::string folderName = static_cast<std::string>(o2::aod::femtoworldparticle::ParticleTypeName[part]);
    /// \todo initialize histograms for children tracks of v0s
    mHistogramRegistry->add((folderName + "/hPt").c_str(), "; #it{p}_{T} (GeV/#it{c}); Entries", kTH1F, {{1000, 0, 10}});
    mHistogramRegistry->add((folderName + "/hEta").c_str(), "; #eta; Entries", kTH1F, {{1000, -1, 1}});
    mHistogramRegistry->add((folderName + "/hPhi").c_str(), "; #phi; Entries", kTH1F, {{1000, 0, 2. * M_PI}});
    // mHistogramRegistry->add((folderName + "/hDaughDCA").c_str(), "; DCA^{daugh} (cm); Entries", kTH1F, {{1000, 0, 10}});
    // mHistogramRegistry->add((folderName + "/hTransRadius").c_str(), "; #it{r}_{xy} (cm); Entries", kTH1F, {{1500, 0, 150}});
    // mHistogramRegistry->add((folderName + "/hDecayVtxX").c_str(), "; #it{Vtx}_{x} (cm); Entries", kTH1F, {{2000, 0, 200}});
    // mHistogramRegistry->add((folderName + "/hDecayVtxY").c_str(), "; #it{Vtx}_{y} (cm)); Entries", kTH1F, {{2000, 0, 200}});
    // mHistogramRegistry->add((folderName + "/hDecayVtxZ").c_str(), "; #it{Vtx}_{z} (cm); Entries", kTH1F, {{2000, 0, 200}});
    // mHistogramRegistry->add((folderName + "/hCPA").c_str(), "; #it{cos #theta_{p}}; Entries", kTH1F, {{1000, 0.9, 1.}});
    // mHistogramRegistry->add((folderName + "/hCPAvsPt").c_str(), "; #it{p}_{T} (GeV/#it{c}); #it{cos #theta_{p}}", kTH2F, {{8, 0.3, 4.3}, {1000, 0.9, 1.}});
    mHistogramRegistry->add((folderName + "/hInvMassPhi").c_str(), "", kTH1F, {massAxisPhi});
    // mHistogramRegistry->add((folderName + "/hInvMassAntiPhi").c_str(), "", kTH1F, {massAxisAntiPhi});
    // mHistogramRegistry->add((folderName + "/hInvMassPhiPhi").c_str(), "", kTH2F, {massAxisPhi, massAxisPhi});

    PosDaughTrack.init<aod::femtoworldparticle::ParticleType::kPhiChild, aod::femtoworldparticle::TrackType::kPosChild, aod::femtoworldparticle::cutContainerType>(mHistogramRegistry);
    NegDaughTrack.init<aod::femtoworldparticle::ParticleType::kPhiChild, aod::femtoworldparticle::TrackType::kNegChild, aod::femtoworldparticle::cutContainerType>(mHistogramRegistry);

    mHistogramRegistry->add("PhiQA/hInvMassPhiNoCuts", "No cuts Phi", kTH1F, {massAxisPhi});
    mHistogramRegistry->add("PhiQA/hInvMassPhiInvMassCut", "Invariant mass cut Phi", kTH1F, {massAxisPhi});
    mHistogramRegistry->add("PhiQA/hInvMassPhiPtMin", "Minimum Pt cut", kTH1F, {massAxisPhi});
    mHistogramRegistry->add("PhiQA/hInvMassPhiPtMax", "Maximum Pt cut", kTH1F, {massAxisPhi});
    mHistogramRegistry->add("PhiQA/hInvMassPhiDCAPhiDaugh", "Phi-daughters DCA cut", kTH1F, {massAxisPhi});
    mHistogramRegistry->add("PhiQA/hInvMassPhiCPA", "CPA cut", kTH1F, {massAxisPhi});
    mHistogramRegistry->add("PhiQA/hInvMassPhiTranRadMin", "Minimum transverse radius cut", kTH1F, {massAxisPhi});
    mHistogramRegistry->add("PhiQA/hInvMassPhiTranRadMax", "Maximum transverse radius cut", kTH1F, {massAxisPhi});
    mHistogramRegistry->add("PhiQA/hInvMassPhiDecVtxMax", "Maximum distance on  decay vertex cut", kTH1F, {massAxisPhi});
  }
  /// check whether the most open cuts are fulfilled - most of this should have already be done by the filters
  nPtPhiMinSel = getNSelections(femtoWorldPhiSelection::kpTPhiMin);
  nPtPhiMaxSel = getNSelections(femtoWorldPhiSelection::kpTPhiMax);
  nDCAPhiDaughMax = getNSelections(femtoWorldPhiSelection::kDCAPhiDaughMax);
  nCPAPhiMin = getNSelections(femtoWorldPhiSelection::kCPAPhiMin);
  nTranRadPhiMin = getNSelections(femtoWorldPhiSelection::kTranRadPhiMin);
  nTranRadPhiMax = getNSelections(femtoWorldPhiSelection::kTranRadPhiMax);
  nDecVtxMax = getNSelections(femtoWorldPhiSelection::kDecVtxMax);

  pTPhiMin = getMinimalSelection(femtoWorldPhiSelection::kpTPhiMin, femtoWorldSelection::kLowerLimit);
  pTPhiMax = getMinimalSelection(femtoWorldPhiSelection::kpTPhiMax, femtoWorldSelection::kUpperLimit);
  DCAPhiDaughMax = getMinimalSelection(femtoWorldPhiSelection::kDCAPhiDaughMax, femtoWorldSelection::kUpperLimit);
  CPAPhiMin = getMinimalSelection(femtoWorldPhiSelection::kCPAPhiMin, femtoWorldSelection::kLowerLimit);
  TranRadPhiMin = getMinimalSelection(femtoWorldPhiSelection::kTranRadPhiMin, femtoWorldSelection::kLowerLimit);
  TranRadPhiMax = getMinimalSelection(femtoWorldPhiSelection::kTranRadPhiMax, femtoWorldSelection::kUpperLimit);
  DecVtxMax = getMinimalSelection(femtoWorldPhiSelection::kDecVtxMax, femtoWorldSelection::kAbsUpperLimit);
}

// Phi initialization
template <o2::aod::femtoworldparticle::ParticleType part, o2::aod::femtoworldparticle::ParticleType daugh, typename cutContainerType>
void FemtoWorldPhiSelection::initPhi(HistogramRegistry* registry)
{
  if (registry) {
    mHistogramRegistry = registry;
    fillSelectionHistogram<part>();
    fillSelectionHistogram<daugh>();

    AxisSpec massAxisPhi = {60000, 0.0f, 3.0f, "m_{#Phi} (GeV/#it{c}^{2})"};
    // AxisSpec massAxisAntiPhi = {600, 0.0f, 3.0f, "m_{#bar{#Phi}} (GeV/#it{c}^{2})"};

    /// \todo this should be an automatic check in the parent class, and the return type should be templated
    size_t nSelections = getNSelections();
    if (nSelections > 8 * sizeof(cutContainerType)) {
      LOG(fatal) << "FemtoWorldPhiCuts: Number of selections to large for your container - quitting!";
    }
    std::string folderName = static_cast<std::string>(o2::aod::femtoworldparticle::ParticleTypeName[part]);
    /// \todo initialize histograms for children tracks of v0s
    mHistogramRegistry->add((folderName + "/hPtPhi").c_str(), "; #it{p}_{T} (GeV/#it{c}); Entries", kTH1F, {{1000, 0, 10}});
    mHistogramRegistry->add((folderName + "/hEtaPhi").c_str(), "; #eta; Entries", kTH1F, {{1000, -1, 1}});
    mHistogramRegistry->add((folderName + "/hPhiPhi").c_str(), "; #phi; Entries", kTH1F, {{1000, 0, 2. * M_PI}});
    // mHistogramRegistry->add((folderName + "/hDaughDCA").c_str(), "; DCA^{daugh} (cm); Entries", kTH1F, {{1000, 0, 10}});
    // mHistogramRegistry->add((folderName + "/hTransRadius").c_str(), "; #it{r}_{xy} (cm); Entries", kTH1F, {{1500, 0, 150}});
    // mHistogramRegistry->add((folderName + "/hDecayVtxX").c_str(), "; #it{Vtx}_{x} (cm); Entries", kTH1F, {{2000, 0, 200}});
    // mHistogramRegistry->add((folderName + "/hDecayVtxY").c_str(), "; #it{Vtx}_{y} (cm)); Entries", kTH1F, {{2000, 0, 200}});
    // mHistogramRegistry->add((folderName + "/hDecayVtxZ").c_str(), "; #it{Vtx}_{z} (cm); Entries", kTH1F, {{2000, 0, 200}});
    // mHistogramRegistry->add((folderName + "/hCPA").c_str(), "; #it{cos #theta_{p}}; Entries", kTH1F, {{1000, 0.9, 1.}});
    // mHistogramRegistry->add((folderName + "/hCPAvsPt").c_str(), "; #it{p}_{T} (GeV/#it{c}); #it{cos #theta_{p}}", kTH2F, {{8, 0.3, 4.3}, {1000, 0.9, 1.}});
    mHistogramRegistry->add((folderName + "/hInvMassPhi").c_str(), "", kTH1F, {massAxisPhi});
    PosDaughTrack.init<aod::femtoworldparticle::ParticleType::kPhiChild, aod::femtoworldparticle::TrackType::kPosChild, aod::femtoworldparticle::cutContainerType>(mHistogramRegistry);
    NegDaughTrack.init<aod::femtoworldparticle::ParticleType::kPhiChild, aod::femtoworldparticle::TrackType::kNegChild, aod::femtoworldparticle::cutContainerType>(mHistogramRegistry);

    mHistogramRegistry->add("PhiQA/hInvMasPhiNoCuts", "No cuts", kTH1F, {massAxisPhi});
    mHistogramRegistry->add("PhiQA/hInvMassPhiInvMassCut", "Invariant mass cut", kTH1F, {massAxisPhi});
    mHistogramRegistry->add("PhiQA/hInvMassPhiPtMin", "Minimum Pt cut", kTH1F, {massAxisPhi});
    mHistogramRegistry->add("PhiQA/hInvMassPhiPtMax", "Maximum Pt cut", kTH1F, {massAxisPhi});
    // mHistogramRegistry->add("PhiQA/hInvMassPhiDCAPhiDaugh", "Phi-daughters DCA cut", kTH1F, {massAxisPhi});
    // mHistogramRegistry->add("PhiQA/hInvMassPhiCPA", "CPA cut", kTH1F, {massAxisPhi});
    // mHistogramRegistry->add("PhiQA/hInvMassPhiTranRadMin", "Minimum transverse radius cut", kTH1F, {massAxisPhi});
    // mHistogramRegistry->add("PhiQA/hInvMassPhiTranRadMax", "Maximum transverse radius cut", kTH1F, {massAxisPhi});
    // mHistogramRegistry->add("PhiQA/hInvMassPhiDecVtxMax", "Maximum distance on  decay vertex cut", kTH1F, {massAxisPhi});
  }
  /// check whether the most open cuts are fulfilled - most of this should have already be done by the filters
  nPtPhiMinSel = getNSelections(femtoWorldPhiSelection::kpTPhiMin);
  nPtPhiMaxSel = getNSelections(femtoWorldPhiSelection::kpTPhiMax);
  nDCAPhiDaughMax = getNSelections(femtoWorldPhiSelection::kDCAPhiDaughMax);
  nCPAPhiMin = getNSelections(femtoWorldPhiSelection::kCPAPhiMin);
  nTranRadPhiMin = getNSelections(femtoWorldPhiSelection::kTranRadPhiMin);
  nTranRadPhiMax = getNSelections(femtoWorldPhiSelection::kTranRadPhiMax);
  nDecVtxMax = getNSelections(femtoWorldPhiSelection::kDecVtxMax);

  pTPhiMin = getMinimalSelection(femtoWorldPhiSelection::kpTPhiMin, femtoWorldSelection::kLowerLimit);
  pTPhiMax = getMinimalSelection(femtoWorldPhiSelection::kpTPhiMax, femtoWorldSelection::kUpperLimit);
  DCAPhiDaughMax = getMinimalSelection(femtoWorldPhiSelection::kDCAPhiDaughMax, femtoWorldSelection::kUpperLimit);
  CPAPhiMin = getMinimalSelection(femtoWorldPhiSelection::kCPAPhiMin, femtoWorldSelection::kLowerLimit);
  TranRadPhiMin = getMinimalSelection(femtoWorldPhiSelection::kTranRadPhiMin, femtoWorldSelection::kLowerLimit);
  TranRadPhiMax = getMinimalSelection(femtoWorldPhiSelection::kTranRadPhiMax, femtoWorldSelection::kUpperLimit);
  DecVtxMax = getMinimalSelection(femtoWorldPhiSelection::kDecVtxMax, femtoWorldSelection::kAbsUpperLimit);
}

template <typename C, typename V, typename T>
bool FemtoWorldPhiSelection::isSelectedMinimal(C const& col, V const& v0, T const& posTrack, T const& negTrack)
{
  const auto signPos = posTrack.sign();
  const auto signNeg = negTrack.sign();
  if (signPos < 0 || signNeg > 0) {
    LOGF(error, "-Something wrong in isSelectedMinimal--\n");
    LOGF(error, "ERROR - Wrong sign for Phi daughters\n");
  }
  const float pT = v0.pt();
  const std::vector<float> decVtx = {v0.x(), v0.y(), v0.z()};
  const float tranRad = v0.v0radius();
  const float dcaDaughv0 = v0.dcaPhidaughters();
  const float cpav0 = v0.v0cosPA();

  const float invMassPhi = v0.mPhi();
  const float invMassAntiPhi = v0.mAntiPhi();

  if (((invMassPhi < fInvMassLowLimit) || (invMassPhi > fInvMassUpLimit)) && ((invMassAntiPhi < fInvMassLowLimit) || (invMassAntiPhi > fInvMassUpLimit))) {
    return false;
  }
  if (fRejectKaon) {
    const float invMassKaon = v0.mK0Short();
    if (invMassKaon > fInvMassKaonLowLimit && invMassKaon < fInvMassKaonUpLimit) {
      return false;
    }
  }
  if (nPtPhiMinSel > 0 && pT < pTPhiMin) {
    return false;
  }
  if (nPtPhiMaxSel > 0 && pT > pTPhiMax) {
    return false;
  }
  if (nDCAPhiDaughMax > 0 && dcaDaughv0 > DCAPhiDaughMax) {
    return false;
  }
  if (nCPAPhiMin > 0 && cpav0 < CPAPhiMin) {
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
  // v0
  auto nSigmaPiNeg = negTrack.tpcNSigmaPi();
  auto nSigmaPrPos = posTrack.tpcNSigmaPr();
  if (!((abs(nSigmaPrNeg) < nSigmaPIDMax) && (abs(nSigmaPiPos) < nSigmaPIDMax)) && !((abs(nSigmaPrPos) < nSigmaPIDMax) && (abs(nSigmaPiNeg) < nSigmaPIDMax))) {
    return false;
  }

  return true;
}

template <typename C, typename V, typename T>
void FemtoWorldPhiSelection::fillPhiQA(C const& col, V const& v0, T const& posTrack, T const& negTrack)
{
  const auto signPos = posTrack.sign();
  const auto signNeg = negTrack.sign();
  if (signPos < 0 || signNeg > 0) {
    LOGF(error, "-Something wrong in isSelectedMinimal--\n");
    LOGF(error, "ERROR - Wrong sign for Phi daughters\n");
  }
  const float pT = v0.pt();
  const std::vector<float> decVtx = {v0.x(), v0.y(), v0.z()};
  const float tranRad = v0.v0radius();
  const float dcaDaughv0 = v0.dcaPhidaughters();
  const float cpav0 = v0.v0cosPA();

  const float invMassPhi = v0.mPhi();

  mHistogramRegistry->fill(HIST("PhiQA/hInvMassPhiNoCuts"), v0.mPhi());

  if ((invMassPhi > fInvMassLowLimit) && (invMassPhi < fInvMassUpLimit)) {
    mHistogramRegistry->fill(HIST("PhiQA/hInvMassPhiInvMassCut"), v0.mPhi());
  }

  if (pT > pTPhiMin) {
    mHistogramRegistry->fill(HIST("PhiQA/hInvMassPhiPtMin"), v0.mPhi());
  }
  if (pT < pTPhiMax) {
    mHistogramRegistry->fill(HIST("PhiQA/hInvMassPhiPtMax"), v0.mPhi());
  }
  if (dcaDaughv0 < DCAPhiDaughMax) {
    mHistogramRegistry->fill(HIST("PhiQA/hInvMassPhiDCAPhiDaugh"), v0.mPhi());
  }
  if (cpav0 > CPAPhiMin) {
    mHistogramRegistry->fill(HIST("PhiQA/hInvMassPhiCPA"), v0.mPhi());
  }
  if (tranRad > TranRadPhiMin) {
    mHistogramRegistry->fill(HIST("PhiQA/hInvMassPhiTranRadMin"), v0.mPhi());
  }
  if (tranRad < TranRadPhiMax) {
    mHistogramRegistry->fill(HIST("PhiQA/hInvMassPhiTranRadMax"), v0.mPhi());
  }
  bool write = true;
  for (size_t i = 0; i < decVtx.size(); i++) {
    write = write && (decVtx.at(i) < DecVtxMax);
  }
  if (write) {
    mHistogramRegistry->fill(HIST("PhiQA/hInvMassPhiDecVtxMax"), v0.mPhi());
  }
}

template <typename C, typename V, typename T, typename P>
void FemtoWorldPhiSelection::fillPhiQAMass(C const& col, V const& MassPhi, T const& posTrack, T const& negTrack, P const& ConfInvMassLowLimit, P const& ConfInvMassUpLimit)
{
  const auto signPos = posTrack.sign();
  const auto signNeg = negTrack.sign();
  if (signPos < 0 || signNeg > 0) {
    LOGF(error, "-Something wrong in isSelectedMinimal--\n");
    LOGF(error, "ERROR - Wrong sign for Phi daughters\n");
  }
  mHistogramRegistry->fill(HIST("PhiQA/hInvMassPhiNoCuts"), MassPhi);
}

/// the CosPA of Phi needs as argument the posXYZ of collisions vertex so we need to pass the collsion as well
template <typename cutContainerType, typename C, typename T>
std::array<cutContainerType, 5> FemtoWorldPhiSelection::getCutContainer(C const& col, T const& posTrack, T const& negTrack)
{
  auto outputPosTrack = PosDaughTrack.getCutContainer<cutContainerType>(posTrack);
  auto outputNegTrack = NegDaughTrack.getCutContainer<cutContainerType>(negTrack);
  cutContainerType output = 0;
  size_t counter = 0;

  // auto lambdaMassNominal = TDatabasePDG::Instance()->GetParticle(321)->Mass(); // FIXME: Get from the common header
  // auto lambdaMassHypothesis = v0.mPhi();
  // auto antiPhiMassHypothesis = v0.mAntiPhi();
  // auto diffPhi = abs(lambdaMassNominal - lambdaMassHypothesis);
  // auto diffAntiPhi = abs(antiPhiMassHypothesis - lambdaMassHypothesis);

  float sign = 0.;
  int nSigmaPIDMax = PosDaughTrack.getSigmaPIDMax();
  auto nSigmaPrNeg = negTrack.tpcNSigmaPr();
  auto nSigmaPiPos = posTrack.tpcNSigmaPi();
  auto nSigmaPiNeg = negTrack.tpcNSigmaPi();
  auto nSigmaPrPos = posTrack.tpcNSigmaPr();
  // check the mass and the PID of daughters
  if (abs(nSigmaPrNeg) < nSigmaPIDMax && abs(nSigmaPiPos) < nSigmaPIDMax) {
    sign = -1.;
  } else if (abs(nSigmaPrPos) < nSigmaPIDMax && abs(nSigmaPiNeg) < nSigmaPIDMax) {
    sign = 1.;
  } else { // if it happens that none of these are true, ignore the invariant mass
    if (abs(nSigmaPrNeg) < nSigmaPIDMax && abs(nSigmaPiPos) < nSigmaPIDMax) {
      sign = -1.;
    } else if (abs(nSigmaPrPos) < nSigmaPIDMax && abs(nSigmaPiNeg) < nSigmaPIDMax) {
      sign = 1.;
    }
  }

  TLorentzVector part1Vec;
  TLorentzVector part2Vec;
  Configurable<int> ConfPDGCodePartOne{"ConfPDGCodePartOne", 321, "Particle 1 - PDG code"};
  Configurable<int> ConfPDGCodePartTwo{"ConfPDGCodePartTwo", 321, "Particle 2 - PDG code"};
  float mMassOne = TDatabasePDG::Instance()->GetParticle(ConfPDGCodePartOne)->Mass(); // FIXME: Get from the PDG service of the common header
  float mMassTwo = TDatabasePDG::Instance()->GetParticle(ConfPDGCodePartTwo)->Mass(); // FIXME: Get from the PDG service of the common header
  part1Vec.SetPtEtaPhiM(posTrack.pt(), posTrack.eta(), posTrack.phi(), mMassOne);
  part2Vec.SetPtEtaPhiM(negTrack.pt(), negTrack.eta(), negTrack.phi(), mMassTwo);

  TLorentzVector sumVec(part1Vec);
  sumVec += part2Vec;

  const auto pT = sumVec.Pt();
  // const auto tranRad = v0.v0radius();
  // const auto dcaDaughv0 = v0.dcaPhidaughters();
  // const auto cpav0 = v0.v0cosPA();
  const std::vector<float> decVtx = {posTrack.x(), posTrack.y(), posTrack.z()};

  float observable = 0.;
  for (auto& sel : mSelections) {
    const auto selVariable = sel.getSelectionVariable();

    if (selVariable == femtoWorldPhiSelection::kDecVtxMax) {
      for (size_t i = 0; i < decVtx.size(); ++i) {
        auto decVtxValue = decVtx.at(i);
        sel.checkSelectionSetBit(decVtxValue, output, counter);
      }
    } else {
      switch (selVariable) {
        case (femtoWorldPhiSelection::kPhiSign):
          observable = sign;
          break;
        case (femtoWorldPhiSelection::kpTPhiMin):
        case (femtoWorldPhiSelection::kpTPhiMax):
          observable = pT;
          break;
        case (femtoWorldPhiSelection::kDCAPhiDaughMax):
          // observable = dcaDaughv0;
          break;
        case (femtoWorldPhiSelection::kCPAPhiMin):
          // observable = cpav0;
          break;
        case (femtoWorldPhiSelection::kTranRadPhiMin):
        case (femtoWorldPhiSelection::kTranRadPhiMax):
          // observable = tranRad;
          break;
        case (femtoWorldPhiSelection::kDecVtxMax):
          break;
      }
      sel.checkSelectionSetBit(observable, output, counter);
    }
  }
  return {output, outputPosTrack.at(femtoWorldTrackSelection::TrackContainerPosition::kCuts),
          outputPosTrack.at(femtoWorldTrackSelection::TrackContainerPosition::kPID),
          outputNegTrack.at(femtoWorldTrackSelection::TrackContainerPosition::kCuts),
          outputNegTrack.at(femtoWorldTrackSelection::TrackContainerPosition::kPID)};
}

template <o2::aod::femtoworldparticle::ParticleType part, o2::aod::femtoworldparticle::ParticleType daugh, typename C, typename V, typename T>
void FemtoWorldPhiSelection::fillQA(C const& col, V const& v0, T const& posTrack, T const& negTrack)
{
  if (mHistogramRegistry) {
    TLorentzVector part1Vec;
    TLorentzVector part2Vec;
    Configurable<int> ConfPDGCodePartOne{"ConfPDGCodePartOne", 321, "Particle 1 - PDG code"};
    Configurable<int> ConfPDGCodePartTwo{"ConfPDGCodePartTwo", 321, "Particle 2 - PDG code"};
    float mMassOne = TDatabasePDG::Instance()->GetParticle(ConfPDGCodePartOne)->Mass(); // FIXME: Get from the PDG service of the common header
    float mMassTwo = TDatabasePDG::Instance()->GetParticle(ConfPDGCodePartTwo)->Mass(); // FIXME: Get from the PDG service of the common header
    part1Vec.SetPtEtaPhiM(posTrack.pt(), posTrack.eta(), posTrack.phi(), mMassOne);
    part2Vec.SetPtEtaPhiM(negTrack.pt(), negTrack.eta(), negTrack.phi(), mMassTwo);

    TLorentzVector sumVec(part1Vec);
    sumVec += part2Vec;

    float phiEta = sumVec.Eta();
    float phiPhi = sumVec.Phi();
    float phiPt = sumVec.Pt();
    float phiInvMass = sumVec.M();

    mHistogramRegistry->fill(HIST(o2::aod::femtoworldparticle::ParticleTypeName[part]) + HIST("/hPt"), phiPt);
    mHistogramRegistry->fill(HIST(o2::aod::femtoworldparticle::ParticleTypeName[part]) + HIST("/hEta"), phiEta);
    mHistogramRegistry->fill(HIST(o2::aod::femtoworldparticle::ParticleTypeName[part]) + HIST("/hPhi"), phiPhi);
    // mHistogramRegistry->fill(HIST(o2::aod::femtoworldparticle::ParticleTypeName[part]) + HIST("/hDaughDCA"), v0.dcaPhidaughters());
    // mHistogramRegistry->fill(HIST(o2::aod::femtoworldparticle::ParticleTypeName[part]) + HIST("/hTransRadius"), v0.v0radius());
    //  mHistogramRegistry->fill(HIST(o2::aod::femtoworldparticle::ParticleTypeName[part]) + HIST("/hDecayVtxX"), v0.x());
    // mHistogramRegistry->fill(HIST(o2::aod::femtoworldparticle::ParticleTypeName[part]) + HIST("/hDecayVtxY"), v0.y());
    // mHistogramRegistry->fill(HIST(o2::aod::femtoworldparticle::ParticleTypeName[part]) + HIST("/hDecayVtxZ"), v0.z());
    mHistogramRegistry->fill(HIST(o2::aod::femtoworldparticle::ParticleTypeName[part]) + HIST("/hInvMassPhi"), phiInvMass);
    // mHistogramRegistry->fill(HIST(o2::aod::femtoworldparticle::ParticleTypeName[part]) + HIST("/hInvMassAntiPhi"), v0.mAntiPhi());
    // mHistogramRegistry->fill(HIST(o2::aod::femtoworldparticle::ParticleTypeName[part]) + HIST("/hInvMassPhiPhi"), v0.mPhi(), v0.mPhi());
  }

  PosDaughTrack.fillQA<aod::femtoworldparticle::ParticleType::kPhiChild, aod::femtoworldparticle::TrackType::kPosChild>(posTrack);
  NegDaughTrack.fillQA<aod::femtoworldparticle::ParticleType::kPhiChild, aod::femtoworldparticle::TrackType::kNegChild>(negTrack);
}

template <o2::aod::femtoworldparticle::ParticleType part, o2::aod::femtoworldparticle::ParticleType daugh, typename C, typename V, typename T>
void FemtoWorldPhiSelection::fillQAPhi(C const& col, V const& v0, T const& posTrack, T const& negTrack)
{
  if (mHistogramRegistry) {
    mHistogramRegistry->fill(HIST(o2::aod::femtoworldparticle::ParticleTypeName[part]) + HIST("/hPtPhi"), v0.pt());
    mHistogramRegistry->fill(HIST(o2::aod::femtoworldparticle::ParticleTypeName[part]) + HIST("/hEtaPhi"), v0.eta());
    mHistogramRegistry->fill(HIST(o2::aod::femtoworldparticle::ParticleTypeName[part]) + HIST("/hPhiPhi"), v0.phi());
    // mHistogramRegistry->fill(HIST(o2::aod::femtoworldparticle::ParticleTypeName[part]) + HIST("/hDaughDCA"), v0.dcaPhidaughters());
    // mHistogramRegistry->fill(HIST(o2::aod::femtoworldparticle::ParticleTypeName[part]) + HIST("/hTransRadius"), v0.v0radius());
    // mHistogramRegistry->fill(HIST(o2::aod::femtoworldparticle::ParticleTypeName[part]) + HIST("/hDecayVtxX"), v0.x());
    // mHistogramRegistry->fill(HIST(o2::aod::femtoworldparticle::ParticleTypeName[part]) + HIST("/hDecayVtxY"), v0.y());
    // mHistogramRegistry->fill(HIST(o2::aod::femtoworldparticle::ParticleTypeName[part]) + HIST("/hDecayVtxZ"), v0.z());
    // mHistogramRegistry->fill(HIST(o2::aod::femtoworldparticle::ParticleTypeName[part]) + HIST("/hCPA"), v0.v0cosPA());
    // mHistogramRegistry->fill(HIST(o2::aod::femtoworldparticle::ParticleTypeName[part]) + HIST("/hCPAvsPt"), v0.pt(), v0.v0cosPA());
    // mHistogramRegistry->fill(HIST(o2::aod::femtoworldparticle::ParticleTypeName[part]) + HIST("/hInvMassPhi"), v0.mPhi());
    // mHistogramRegistry->fill(HIST(o2::aod::femtoworldparticle::ParticleTypeName[part]) + HIST("/hInvMassAntiPhi"), v0.mAntiPhi());
    // mHistogramRegistry->fill(HIST(o2::aod::femtoworldparticle::ParticleTypeName[part]) + HIST("/hInvMassPhiAntiPhi"), v0.mPhi(), v0.mAntiPhi());
  }

  PosDaughTrack.fillQA<aod::femtoworldparticle::ParticleType::kPhiChild, aod::femtoworldparticle::TrackType::kPosChild>(posTrack);
  NegDaughTrack.fillQA<aod::femtoworldparticle::ParticleType::kPhiChild, aod::femtoworldparticle::TrackType::kNegChild>(negTrack);
}

} // namespace o2::analysis::femtoWorld

#endif // PWGCF_FEMTOWORLD_CORE_FEMTOWORLDPHISELECTION_H_
