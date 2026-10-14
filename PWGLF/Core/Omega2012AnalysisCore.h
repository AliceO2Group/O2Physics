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
///
/// \file Omega2012AnalysisCore.h
/// \brief Omega(2012) selection, candidate enumeration and truth classification, separately for the two decay modes
/// \author Bong-Hwi Lim <bong-hwi.lim@cern.ch>
///
/// Mode A (DecayMode::XiK0s): Omega(2012)- -> Xi- K0S.
/// Mode B (DecayMode::Xi1530K): Omega(2012)- -> Xi(1530)0 K- -> Xi- pi+ K- (charged kaon track).
/// The two modes are never merged: each has its own enumeration, hooks, cut flow, pass bits and truth classification.
/// The Xi selection is common to both modes. The core owns the selection and its cut-flow instrumentation
/// (CutFlow/*, ML/*); the tasks own their output histograms and fill them from the hooks.

#ifndef PWGLF_CORE_OMEGA2012ANALYSISCORE_H_
#define PWGLF_CORE_OMEGA2012ANALYSISCORE_H_

#include "PWGLF/Core/Omega2012MlFeatures.h"
#include "PWGLF/Core/ResoAnalysisSelectionCore.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/Logger.h>
#include <Framework/StructToTuple.h>

#include <Math/GenVector/VectorUtil.h>
#include <Math/Vector4D.h> // IWYU pragma: keep (do not replace with Math/Vector4Dfwd.h)
#include <Math/Vector4Dfwd.h>
#include <TH1.h>
#include <TH2.h>
#include <TPDGCode.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace o2::analysis::omega2012
{

using o2::analysis::resonance::EventCuts;
using o2::analysis::resonance::PIDCutConfig;
using o2::analysis::resonance::ResoAnalysisSelectionCore;
using o2::analysis::resonance::TrackCuts;
using o2::analysis::resonance::TrackStage;

enum class DecayMode : uint8_t {
  XiK0s = 0,  // Omega(2012)- -> Xi- K0S
  Xi1530K = 1 // Omega(2012)- -> Xi(1530)0 K- -> Xi- pi+ K-
};

inline constexpr int PdgXi1530Zero = 3324;   // Xi(1530)0, not in o2::constants::physics::Pdg
inline constexpr float SmallNumber = 1e-10f; // avoids division by zero (value of the original task)
inline constexpr float MaxDCAV0ToPV = 1.0f;  // maximum K0S DCA to PV (value of the original task)
inline constexpr float CascDCAzPtExponent = -1.1f;
inline constexpr const char* TrackCutsPrefix = "trk."; // JSON prefix of the mode-B track selection
inline constexpr int NXiK0sCandidateStages = 4;
inline constexpr int NXi1530KCandidateStages = 4;
inline constexpr int NGeneratedChannels = 3;

// Last stage passed by a cascade (Xi selection, common to both modes).
enum XiStage : int {
  kXiInput = 0,        // failed |eta| or pT
  kXiKinematics = 1,   // passed |eta| and pT
  kXiDCA = 2,          // passed the cascade DCA to PV
  kXiV0Topology = 3,   // passed the Lambda topology and mass window
  kXiCascTopology = 4, // passed the cascade topology
  kXiMass = 5,         // passed the Xi mass window
  kXiSelected = kXiMass,
  kXiNStages = 6
};

// Last stage passed by a V0 (K0S selection of mode A).
enum K0sStage : int {
  kK0sInput = 0,      // failed |eta| or pT
  kK0sKinematics = 1, // passed |eta| and pT
  kK0sTopology = 2,   // passed cosPA, daughter DCAs, radius and DCA to PV
  kK0sLifetime = 3,   // passed the proper lifetime and the minimum qT
  kK0sMass = 4,       // passed the K0S mass window and the (anti)Lambda rejection
  kK0sDaughters = 5,  // passed the daughter pion TPC nSigma and crossed rows
  kK0sArmenteros = 6, // passed the Armenteros qT > coefficient * |alpha| cut
  kK0sSelected = kK0sArmenteros,
  kK0sNStages = 7
};

// The value is the species index of the PID configuration in ResoAnalysisSelectionCore.
enum class Species : int {
  Pion = 0,
  Kaon = 1
};

// Cumulative pass bits of mode A (see LFOmega2012MlTables.h)
enum XiK0sPassBit : uint16_t {
  kXiK0sPassLoose = 1,
  kXiK0sPassXi = 2,
  kXiK0sPassK0s = 4,
  kXiK0sPassKinematic = 8
};
inline constexpr uint16_t XiK0sPassBitsSelected = kXiK0sPassLoose | kXiK0sPassXi | kXiK0sPassK0s | kXiK0sPassKinematic;

// Cumulative pass bits of mode B (see LFOmega2012MlTables.h)
enum Xi1530KPassBit : uint16_t {
  kXi1530KPassLoose = 1,
  kXi1530KPassXi = 2,
  kXi1530KPassQuality = 4,
  kXi1530KPassPID = 8,
  kXi1530KPassWindow = 16
};
inline constexpr uint16_t Xi1530KPassBitsSelected = kXi1530KPassLoose | kXi1530KPassXi | kXi1530KPassQuality | kXi1530KPassPID | kXi1530KPassWindow;

enum class XiK0sTruth : uint8_t {
  None = 0,
  Matched = 1
};

enum class Xi1530KTruth : uint8_t {
  None = 0,
  Matched = 1
};

// Immediate decay channel of a generated Omega(2012)
enum class GeneratedChannel : uint8_t {
  Other = 0,
  XiK0s = 1,
  Xi1530K = 2
};

// Configurable groups without prefix: the JSON keys are the plain configurable names of the original task.

/// Xi (cascade) selection, common to both modes. cMinPtcut is also the K0S minimum pT (as in the original task).
struct XiCuts : o2::framework::ConfigurableGroup {
  o2::framework::Configurable<float> cMinPtcut{"cMinPtcut", 0.15, "Minimum pT for candidates"};
  o2::framework::Configurable<float> cMaxEtaCut{"cMaxEtaCut", 0.8, "Maximum |eta|"};
  o2::framework::Configurable<float> cDCAxyToPVByPtCascP0{"cDCAxyToPVByPtCascP0", 999., "Cascade DCAxy p0"};
  o2::framework::Configurable<float> cDCAxyToPVByPtCascExp{"cDCAxyToPVByPtCascExp", 1., "Cascade DCAxy exp"};
  o2::framework::Configurable<bool> cDCAxyToPVAsPtForCasc{"cDCAxyToPVAsPtForCasc", true, "Use pt-dep DCAxy cut (casc)"};
  o2::framework::Configurable<bool> cDCAzToPVAsPtForCasc{"cDCAzToPVAsPtForCasc", true, "Use pt-dep DCAz cut (casc)"};
  // V0 topology inside cascade (Lambda)
  o2::framework::Configurable<float> cDCALambdaDaugtherscut{"cDCALambdaDaugtherscut", 0.7, "Λ daughters DCA cut"};
  o2::framework::Configurable<float> cDCALambdaToPVcut{"cDCALambdaToPVcut", 0.02, "Λ DCA to PV min"};
  o2::framework::Configurable<float> cDCAPionToPVcut{"cDCAPionToPVcut", 0.06, "π DCA to PV min"};
  o2::framework::Configurable<float> cDCAProtonToPVcut{"cDCAProtonToPVcut", 0.07, "p DCA to PV min"};
  o2::framework::Configurable<float> cV0CosPACutPtDepP0{"cV0CosPACutPtDepP0", 0.25, "V0 CosPA p0"};
  o2::framework::Configurable<float> cV0CosPACutPtDepP1{"cV0CosPACutPtDepP1", 0.022, "V0 CosPA p1"};
  o2::framework::Configurable<float> cMaxV0radiuscut{"cMaxV0radiuscut", 200., "V0 radius max"};
  o2::framework::Configurable<float> cMinV0radiuscut{"cMinV0radiuscut", 2.5, "V0 radius min"};
  o2::framework::Configurable<float> cMasswindowV0cut{"cMasswindowV0cut", 0.005, "Λ mass window for cascade V0"};
  // Cascade topology
  o2::framework::Configurable<float> cDCABachlorToPVcut{"cDCABachlorToPVcut", 0.06, "Bachelor DCA to PV min"};
  o2::framework::Configurable<float> cDCAXiDaugthersCutPtRangeLower{"cDCAXiDaugthersCutPtRangeLower", 1., "Xi pt low boundary"};
  o2::framework::Configurable<float> cDCAXiDaugthersCutPtRangeUpper{"cDCAXiDaugthersCutPtRangeUpper", 4., "Xi pt high boundary"};
  o2::framework::Configurable<float> cDCAXiDaugthersCutPtDepLower{"cDCAXiDaugthersCutPtDepLower", 0.8, "Xi daugh DCA (pt<low)"};
  o2::framework::Configurable<float> cDCAXiDaugthersCutPtDepMiddle{"cDCAXiDaugthersCutPtDepMiddle", 0.5, "Xi daugh DCA (low<=pt<high)"};
  o2::framework::Configurable<float> cDCAXiDaugthersCutPtDepUpper{"cDCAXiDaugthersCutPtDepUpper", 0.2, "Xi daugh DCA (pt>=high)"};
  o2::framework::Configurable<float> cCosPACascCutPtDepP0{"cCosPACascCutPtDepP0", 0.2, "Cascade CosPA p0"};
  o2::framework::Configurable<float> cCosPACascCutPtDepP1{"cCosPACascCutPtDepP1", 0.022, "Cascade CosPA p1"};
  o2::framework::Configurable<float> cMaxCascradiuscut{"cMaxCascradiuscut", 200., "Cascade radius max"};
  o2::framework::Configurable<float> cMinCascradiuscut{"cMinCascradiuscut", 1.1, "Cascade radius min"};
  o2::framework::Configurable<float> cMasswindowCasccut{"cMasswindowCasccut", 0.008, "Xi mass window"};
  o2::framework::Configurable<float> cMassXiminus{"cMassXiminus", 1.32171, "Xi mass (GeV/c^2)"}; // PDG
};

/// K0S selection and Xi-K0S kinematic (opening-angle) cut of mode A
struct K0sCuts : o2::framework::ConfigurableGroup {
  o2::framework::Configurable<bool> cKinCuts{"cKinCuts", false, "Kinematic cuts for Xi-K0s opening angle"};
  o2::framework::Configurable<std::vector<float>> cKinCutsPt{"cKinCutsPt", {0.0, 0.4, 0.6, 0.8, 1.0, 1.4, 1.8, 2.2, 2.6, 3.0, 4.0, 5.0, 6.0, 1e10}, "Omega(2012) pT bins for kinematic cuts"};
  o2::framework::Configurable<std::vector<float>> cKinLowerCutsAlpha{"cKinLowerCutsAlpha", {1.5, 1.0, 0.5, 0.3, 0.2, 0.15, 0.1, 0.08, 0.07, 0.06, 0.04, 0.02, 0.02}, "Lower cut on Xi-K0s opening angle"};
  o2::framework::Configurable<std::vector<float>> cKinUpperCutsAlpha{"cKinUpperCutsAlpha", {3.0, 2.0, 1.5, 1.4, 1.0, 0.8, 0.6, 0.5, 0.45, 0.35, 0.3, 0.25, 0.2}, "Upper cut on Xi-K0s opening angle"};
  o2::framework::Configurable<double> cK0sMinCosPA{"cK0sMinCosPA", 0.98, "K0s minimum pointing angle cosine"};
  o2::framework::Configurable<double> cK0sMaxDaughDCA{"cK0sMaxDaughDCA", 0.5, "K0s daughter DCA Maximum"};
  o2::framework::Configurable<double> cK0sMassWindow{"cK0sMassWindow", 0.025, "Mass window for K0s selection (GeV/c^2)"};
  o2::framework::Configurable<double> cMaxV0Etacut{"cMaxV0Etacut", 0.8, "V0 maximum eta cut"};
  o2::framework::Configurable<float> cK0sProperLifetimeMax{"cK0sProperLifetimeMax", 20.0, "K0s proper lifetime max (cm/c)"};
  o2::framework::Configurable<float> cK0sArmenterosQtMin{"cK0sArmenterosQtMin", 0.0, "K0s Armenteros qt min"};
  o2::framework::Configurable<float> cK0sArmenterosAlphaCoeff{"cK0sArmenterosAlphaCoeff", 0.2, "K0s Armenteros alpha coefficient"};
  o2::framework::Configurable<float> cK0sDauPosDCAtoPVMin{"cK0sDauPosDCAtoPVMin", 0.05, "K0s positive daughter DCA to PV min"};
  o2::framework::Configurable<float> cK0sDauNegDCAtoPVMin{"cK0sDauNegDCAtoPVMin", 0.05, "K0s negative daughter DCA to PV min"};
  o2::framework::Configurable<float> cK0sRadiusMin{"cK0sRadiusMin", 0.5, "K0s decay radius min"};
  o2::framework::Configurable<float> cK0sRadiusMax{"cK0sRadiusMax", 200.0, "K0s decay radius max"};
  o2::framework::Configurable<bool> cK0sCrossMassRejection{"cK0sCrossMassRejection", true, "Enable Lambda mass rejection for K0s"};
  o2::framework::Configurable<float> cK0sCrossMassRejectionWindow{"cK0sCrossMassRejectionWindow", 0.01, "Lambda mass rejection window for K0s (GeV/c^2)"};
  o2::framework::Configurable<float> cK0sDaughterPiTPCNSigmaMax{"cK0sDaughterPiTPCNSigmaMax", 5.0, "Maximum TPC NSigma for K0s daughter pions"};
  o2::framework::Configurable<int> cK0sPosDaughterMinCrossedRows{"cK0sPosDaughterMinCrossedRows", 50, "Minimum TPC crossed rows for K0s positive daughter"};
  o2::framework::Configurable<int> cK0sNegDaughterMinCrossedRows{"cK0sNegDaughterMinCrossedRows", 50, "Minimum TPC crossed rows for K0s negative daughter"};
};

/// Xi(1530)0 K- selection of mode B
struct Xi1530KCuts : o2::framework::ConfigurableGroup {
  o2::framework::Configurable<float> cXi1530Mass{"cXi1530Mass", 1.53, "Xi(1530) mass (GeV/c^2)"};
  o2::framework::Configurable<float> cXi1530MassWindow{"cXi1530MassWindow", 0.01, "Xi(1530) mass window (GeV/c^2)"};
  o2::framework::Configurable<bool> cXi1530UseMassWindow{"cXi1530UseMassWindow", true, "Require m(Xi pi) inside the Xi(1530) mass window"};
  o2::framework::Configurable<bool> cXi1530KFillWrongSign{"cXi1530KFillWrongSign", true, "Fill the wrong-sign control histograms (xi1530K_wrongSign/)"};
  o2::framework::Configurable<bool> cByPassTOF{"cByPassTOF", false, "Bypass the TOF nSigma selection of the pion and the kaon"};
};

/// Pion PID of mode B (keys of the original three-body pion selection)
struct PionPidCuts : o2::framework::ConfigurableGroup {
  o2::framework::Configurable<float> cPionTPCNSigmaMax{"cPionTPCNSigmaMax", 3.0, "Maximum TPC NSigma for pion"};
  o2::framework::Configurable<float> cPionTOFNSigmaMax{"cPionTOFNSigmaMax", 3.0, "Maximum TOF NSigma for pion"};
  o2::framework::Configurable<bool> cPionUsePtDepPID{"cPionUsePtDepPID", false, "Use pT-dependent PID cuts for pion"};
  o2::framework::Configurable<std::vector<float>> cPionPIDPtBins{"cPionPIDPtBins", {0.0f, 0.5f, 0.8f, 2.0f, 999.0f}, "pT bin edges for pion PID cuts"};
  o2::framework::Configurable<std::vector<float>> cPionTPCNSigmaCuts{"cPionTPCNSigmaCuts", {3.0f, 3.0f, 2.0f, 2.0f}, "TPC NSigma cuts per pT bin (pion)"};
  o2::framework::Configurable<std::vector<float>> cPionTOFNSigmaCuts{"cPionTOFNSigmaCuts", {3.0f, 3.0f, 3.0f, 3.0f}, "TOF NSigma cuts per pT bin (pion)"};
  o2::framework::Configurable<std::vector<int>> cPionTOFRequired{"cPionTOFRequired", {0, 0, 1, 1}, "Require TOF per pT bin (pion)"};
};

/// Kaon PID of mode B (same shape as the pion keys)
struct KaonPidCuts : o2::framework::ConfigurableGroup {
  o2::framework::Configurable<float> cKaonTPCNSigmaMax{"cKaonTPCNSigmaMax", 3.0, "Maximum TPC NSigma for kaon"};
  o2::framework::Configurable<float> cKaonTOFNSigmaMax{"cKaonTOFNSigmaMax", 3.0, "Maximum TOF NSigma for kaon"};
  o2::framework::Configurable<bool> cKaonUsePtDepPID{"cKaonUsePtDepPID", false, "Use pT-dependent PID cuts for kaon"};
  o2::framework::Configurable<std::vector<float>> cKaonPIDPtBins{"cKaonPIDPtBins", {0.0f, 0.5f, 0.8f, 2.0f, 999.0f}, "pT bin edges for kaon PID cuts"};
  o2::framework::Configurable<std::vector<float>> cKaonTPCNSigmaCuts{"cKaonTPCNSigmaCuts", {3.0f, 3.0f, 2.0f, 2.0f}, "TPC NSigma cuts per pT bin (kaon)"};
  o2::framework::Configurable<std::vector<float>> cKaonTOFNSigmaCuts{"cKaonTOFNSigmaCuts", {3.0f, 3.0f, 3.0f, 3.0f}, "TOF NSigma cuts per pT bin (kaon)"};
  o2::framework::Configurable<std::vector<int>> cKaonTOFRequired{"cKaonTOFRequired", {0, 0, 1, 1}, "Require TOF per pT bin (kaon)"};
};

/// Common candidate selection
struct CandidateCuts : o2::framework::ConfigurableGroup {
  o2::framework::Configurable<float> cfgRapidityCut{"cfgRapidityCut", 0.5, "Rapidity cut"};
};

/// Track selection of the mode-B pion and kaon: the common resonance TrackCuts with the JSON prefix "trk."
/// (the Xi selection owns cMinPtcut). The framework prefixes a group only through a `prefix` data member, which
/// TrackCuts does not have and a derived group cannot add (members must belong to one class for the structured
/// binding of the option registration); the names are therefore prefixed here in the same way ("trk.cMinPtcut").
/// Defaults that differ from TrackCuts follow the original pion selection: |eta| < 0.8 (cPionEtaMax),
/// 70 found TPC clusters for ResoTracks (cPionTPCNClusMin) and DCAz < 0.15 cm (cPionDCAzMax was 0.2 cm, which is
/// beyond the 0.15 cm range of the micro001 DCA grid).
inline TrackCuts makeXi1530KTrackCuts()
{
  TrackCuts cuts;
  cuts.cMaxEtacut.value = 0.8f;
  cuts.cMaxDCAzToPVcut.value = 0.15;
  cuts.cfgTPCcluster.value = 70;
  o2::framework::homogeneous_apply_refs<true>(
    [](auto& option) {
      if constexpr (requires { option.name; }) {
        option.name.insert(0, TrackCutsPrefix);
      }
      return true;
    },
    cuts);
  return cuts;
}

/// Process functions enabled in a task; they decide the configuration checks and the registered histograms.
struct ProcessModes {
  bool xiK0s = false;       // any mode-A process function
  bool xi1530K = false;     // any mode-B process function
  bool microTracks = false; // any mode-B process function reading micro tracks (quantised DCA and nSigma)
  bool mcReco = false;      // reconstructed MC
  bool mcGen = false;       // generated Omega(2012) parents (ResoMCParents_001)
  bool mixing = false;      // event mixing
};

/// Loose-stage traversal for the export hook (as the K1 LooseStageOptions).
/// With the defaults and without an export hook, the enumeration applies only the conventional selection.
struct LooseStageOptions {
  bool audit = false;          // fill ML/<mode>/looseCutflow and ML/<mode>/looseMassPtActivity
  bool exportSelected = false; // hand candidates to the export hook at the selected stage instead of the loose stage
};

/// Mode-A candidate handed to the candidate and export hooks.
/// alpha is computed when cKinCuts is on or the values were requested; passesKinCut is true when cKinCuts is off.
struct XiK0sCandidateValues {
  ROOT::Math::PxPyPzEVector omega; // xi + k0s
  ROOT::Math::PxPyPzEVector xi;    // (pt, eta, phi, mXi)
  ROOT::Math::PxPyPzEVector k0s;   // (pt, eta, phi, MassK0Short)
  float alpha = 0.f;               // Xi-K0S opening angle
  bool passesKinCut = true;
  bool inRapidity = false; // |y| < cfgRapidityCut
};

/// Mode-B candidate handed to the candidate and export hooks.
/// massXiK, massPiK, openingAngleXi1530K and cosThetaStar are computed only when the values were requested.
struct Xi1530KCandidateValues {
  ROOT::Math::PxPyPzEVector omega;  // xi + pion + kaon
  ROOT::Math::PxPyPzEVector xi;     // (pt, eta, phi, mXi)
  ROOT::Math::PxPyPzEVector pion;   // (pt, eta, phi, MassPionCharged)
  ROOT::Math::PxPyPzEVector kaon;   // (pt, eta, phi, MassKaonCharged)
  ROOT::Math::PxPyPzEVector xi1530; // xi + pion
  float massXiPi = 0.f;
  float massXiK = 0.f;
  float massPiK = 0.f;
  float openingAngleXi1530K = 0.f;
  float cosThetaStar = 0.f;
  uint8_t chargePattern = o2::analysis::omega2012ml::kSignalPattern;
  bool inXi1530Window = false;
  bool inRapidity = false; // |y| < cfgRapidityCut
};

// Truth classification, separately per mode

/// Mode A: reconstructed Xi and K0S from the same Omega(2012) (criteria of the original task).
/// A K0S from an intermediate K0bar has motherPDG 311 and is not matched; the truth table keeps motherPDG for audit.
template <typename Xi, typename V0>
XiK0sTruth classifyXiK0sTruth(const Xi& xi, const V0& v0)
{
  if (std::abs(xi.pdgCode()) != kXiMinus || std::abs(v0.pdgCode()) != kK0Short) {
    return XiK0sTruth::None;
  }
  if (xi.motherId() < 0 || xi.motherId() != v0.motherId()) {
    return XiK0sTruth::None;
  }
  if (std::abs(xi.motherPDG()) != o2::constants::physics::Pdg::kOmega2012Minus || xi.motherPDG() != v0.motherPDG()) {
    return XiK0sTruth::None;
  }
  return XiK0sTruth::Matched;
}

template <typename Track>
bool hasSibling(const Track& track, int particleId)
{
  if (particleId < 0) {
    return false;
  }
  const auto siblings = track.siblingIds();
  return siblings[0] == particleId || siblings[1] == particleId;
}

/// Mode B: Xi and pion from the same Xi(1530)0, kaon from the Omega(2012) whose other daughter is that Xi(1530)0.
/// Cascades carry no sibling IDs, so the kaon carries the link (siblingIds contains the Xi(1530)0 ID).
template <typename Xi, typename Track>
Xi1530KTruth classifyXi1530KTruth(const Xi& xi, const Track& pion, const Track& kaon)
{
  if (std::abs(xi.pdgCode()) != kXiMinus) {
    return Xi1530KTruth::None;
  }
  const int xiCharge = xi.pdgCode() > 0 ? -1 : 1; // Xi- has PDG code +3312
  if (pion.pdgCode() != -xiCharge * kPiPlus || kaon.pdgCode() != xiCharge * kKPlus) {
    return Xi1530KTruth::None;
  }
  if (xi.motherId() < 0 || xi.motherId() != pion.motherId() || std::abs(xi.motherPDG()) != PdgXi1530Zero || pion.motherPDG() != xi.motherPDG()) {
    return Xi1530KTruth::None;
  }
  // Omega(2012)- and Xi- both have positive PDG codes
  if (kaon.motherId() < 0 || kaon.motherPDG() != -xiCharge * o2::constants::physics::Pdg::kOmega2012Minus) {
    return Xi1530KTruth::None;
  }
  return hasSibling(kaon, xi.motherId()) ? Xi1530KTruth::Matched : Xi1530KTruth::None;
}

/// Immediate channel of a generated Omega(2012) from the PDG codes of its two daughters (charge conjugates included).
inline GeneratedChannel classifyGeneratedOmega2012(int pdg, int daughter1, int daughter2)
{
  const int sign = pdg > 0 ? 1 : -1;
  auto isPair = [&](int a, int b) {
    if (a == sign * PdgXi1530Zero && b == -sign * kKPlus) {
      return GeneratedChannel::Xi1530K;
    }
    if (a == sign * kXiMinus && (b == kK0Short || b == -sign * kK0)) {
      return GeneratedChannel::XiK0s;
    }
    return GeneratedChannel::Other;
  };
  const auto channel = isPair(daughter1, daughter2);
  return channel != GeneratedChannel::Other ? channel : isPair(daughter2, daughter1);
}

/// Omega(2012) selection and candidate enumeration of the Omega(2012) tasks.
/// The task owns the configurable groups and the histogram registry and passes them in init().
class Omega2012AnalysisCore
{
 public:
  void init(o2::framework::HistogramRegistry& histos,
            EventCuts const& eventCuts, XiCuts const& xiCuts, K0sCuts const& k0sCuts,
            Xi1530KCuts const& xi1530KCuts, TrackCuts const& trackCuts,
            PionPidCuts const& pionPidCuts, KaonPidCuts const& kaonPidCuts,
            CandidateCuts const& candidateCuts, ProcessModes const& modes, LooseStageOptions const& looseOptions = {})
  {
    mXiCuts = xiCuts;
    mK0sCuts = k0sCuts;
    mXi1530KCuts = xi1530KCuts;
    mCandidateCuts = candidateCuts;
    mLooseOptions = looseOptions;

    // The order follows Species
    std::vector<PIDCutConfig> pid(2);
    auto& pion = pid[static_cast<int>(Species::Pion)];
    pion.species = "Pion";
    pion.maxTPCnSigma = pionPidCuts.cPionTPCNSigmaMax.value;
    pion.maxTOFnSigma = pionPidCuts.cPionTOFNSigmaMax.value;
    pion.usePtDependent = pionPidCuts.cPionUsePtDepPID.value;
    pion.ptBins = pionPidCuts.cPionPIDPtBins.value;
    pion.tpcNSigmaCuts = pionPidCuts.cPionTPCNSigmaCuts.value;
    pion.tofNSigmaCuts = pionPidCuts.cPionTOFNSigmaCuts.value;
    pion.tofRequired = pionPidCuts.cPionTOFRequired.value;
    pion.maxTPCName = "cPionTPCNSigmaMax";
    pion.maxTOFName = "cPionTOFNSigmaMax";
    pion.tpcCutsName = "cPionTPCNSigmaCuts";
    pion.tofCutsName = "cPionTOFNSigmaCuts";
    auto& kaon = pid[static_cast<int>(Species::Kaon)];
    kaon.species = "Kaon";
    kaon.maxTPCnSigma = kaonPidCuts.cKaonTPCNSigmaMax.value;
    kaon.maxTOFnSigma = kaonPidCuts.cKaonTOFNSigmaMax.value;
    kaon.usePtDependent = kaonPidCuts.cKaonUsePtDepPID.value;
    kaon.ptBins = kaonPidCuts.cKaonPIDPtBins.value;
    kaon.tpcNSigmaCuts = kaonPidCuts.cKaonTPCNSigmaCuts.value;
    kaon.tofNSigmaCuts = kaonPidCuts.cKaonTOFNSigmaCuts.value;
    kaon.tofRequired = kaonPidCuts.cKaonTOFRequired.value;
    kaon.maxTPCName = "cKaonTPCNSigmaMax";
    kaon.maxTOFName = "cKaonTOFNSigmaMax";
    kaon.tpcCutsName = "cKaonTPCNSigmaCuts";
    kaon.tofCutsName = "cKaonTOFNSigmaCuts";
    mSelection.init(eventCuts, trackCuts, std::move(pid), mXi1530KCuts.cByPassTOF.value, modes.microTracks);

    registerHistograms(histos, modes);
  }

  template <typename CollisionType>
  bool passesEventCuts(const CollisionType& collision)
  {
    return mSelection.passesEventCuts(collision);
  }

  template <typename CollisionType>
  bool passesMCEventCuts(const CollisionType& collision)
  {
    return mSelection.passesMCEventCuts(collision);
  }

  [[nodiscard]] bool kinCutsEnabled() const { return mK0sCuts.cKinCuts.value; }
  [[nodiscard]] bool fillWrongSign() const { return mXi1530KCuts.cXi1530KFillWrongSign.value; }

  // Xi selection (common to both modes); equivalent to cascprimaryTrackCut && casctopCut of the original task.
  template <typename CascadeType>
  [[nodiscard]] int xiSelectionStage(const CascadeType& c) const
  {
    const auto& cut = mXiCuts;
    if (std::abs(c.eta()) > cut.cMaxEtaCut.value || std::abs(c.pt()) < cut.cMinPtcut.value) {
      return kXiInput;
    }
    if (cut.cDCAxyToPVAsPtForCasc.value && std::abs(c.dcaXYCascToPV()) > (cut.cDCAxyToPVByPtCascP0.value + cut.cDCAxyToPVByPtCascExp.value * c.pt())) {
      return kXiKinematics;
    }
    // The DCAz cut uses the DCAxy parameters, as in the original task
    if (cut.cDCAzToPVAsPtForCasc.value && std::abs(c.dcaZCascToPV()) > (cut.cDCAxyToPVByPtCascP0.value + cut.cDCAxyToPVByPtCascExp.value * std::pow(c.pt(), CascDCAzPtExponent))) {
      return kXiKinematics;
    }

    // V0 (Lambda) topology inside the cascade
    if (std::abs(c.daughDCA()) > cut.cDCALambdaDaugtherscut.value || std::abs(c.dcav0topv()) < cut.cDCALambdaToPVcut.value) {
      return kXiDCA;
    }
    if (c.sign() < 0) { // Xi-
      if (std::abs(c.dcanegtopv()) < cut.cDCAPionToPVcut.value || std::abs(c.dcapostopv()) < cut.cDCAProtonToPVcut.value) {
        return kXiDCA;
      }
    } else { // Anti-Xi
      if (std::abs(c.dcanegtopv()) < cut.cDCAProtonToPVcut.value || std::abs(c.dcapostopv()) < cut.cDCAPionToPVcut.value) {
        return kXiDCA;
      }
    }
    if (c.v0CosPA() < std::cos(cut.cV0CosPACutPtDepP0.value - cut.cV0CosPACutPtDepP1.value * c.pt())) {
      return kXiDCA;
    }
    if (c.transRadius() > cut.cMaxV0radiuscut.value || c.transRadius() < cut.cMinV0radiuscut.value) {
      return kXiDCA;
    }
    if (std::abs(c.mLambda() - o2::constants::physics::MassLambda) > cut.cMasswindowV0cut.value) {
      return kXiDCA;
    }

    // Cascade topology
    if (std::abs(c.dcabachtopv()) < cut.cDCABachlorToPVcut.value) {
      return kXiV0Topology;
    }
    if (c.pt() < cut.cDCAXiDaugthersCutPtRangeLower.value) {
      if (c.cascDaughDCA() > cut.cDCAXiDaugthersCutPtDepLower.value) {
        return kXiV0Topology;
      }
    } else if (c.pt() < cut.cDCAXiDaugthersCutPtRangeUpper.value) {
      if (c.cascDaughDCA() > cut.cDCAXiDaugthersCutPtDepMiddle.value) {
        return kXiV0Topology;
      }
    } else {
      if (c.cascDaughDCA() > cut.cDCAXiDaugthersCutPtDepUpper.value) {
        return kXiV0Topology;
      }
    }
    if (c.cascCosPA() < std::cos(cut.cCosPACascCutPtDepP0.value - cut.cCosPACascCutPtDepP1.value * c.pt())) {
      return kXiV0Topology;
    }
    if (c.cascTransRadius() > cut.cMaxCascradiuscut.value || c.cascTransRadius() < cut.cMinCascradiuscut.value) {
      return kXiV0Topology;
    }
    if (std::abs(c.mXi() - cut.cMassXiminus.value) > cut.cMasswindowCasccut.value) {
      return kXiCascTopology;
    }
    return kXiMass;
  }

  // K0S proper lifetime (cm/c), the expression of the original task
  template <typename CollisionType, typename V0Type>
  static double k0sProperLifetime(const CollisionType& collision, const V0Type& v0)
  {
    float dx = v0.decayVtxX() - collision.posX();
    float dy = v0.decayVtxY() - collision.posY();
    float dz = v0.decayVtxZ() - collision.posZ();
    float l = std::sqrt(dx * dx + dy * dy + dz * dz);
    float p = std::sqrt(v0.px() * v0.px() + v0.py() * v0.py() + v0.pz() * v0.pz());
    return (l / (p + SmallNumber)) * o2::constants::physics::MassK0Short;
  }

  // K0S selection of mode A; equivalent to v0CutEnhanced of the original task.
  // The collision is the one of the V0 (its primary vertex defines the proper lifetime).
  template <typename CollisionType, typename V0Type>
  [[nodiscard]] int k0sSelectionStage(const CollisionType& collision, const V0Type& v0) const
  {
    const auto& cut = mK0sCuts;
    if (std::abs(v0.eta()) > cut.cMaxV0Etacut.value || v0.pt() < mXiCuts.cMinPtcut.value) {
      return kK0sInput;
    }
    if (v0.v0CosPA() < cut.cK0sMinCosPA.value || v0.daughDCA() > cut.cK0sMaxDaughDCA.value) {
      return kK0sKinematics;
    }
    if (std::abs(v0.dcapostopv()) < cut.cK0sDauPosDCAtoPVMin.value || std::abs(v0.dcanegtopv()) < cut.cK0sDauNegDCAtoPVMin.value) {
      return kK0sKinematics;
    }
    auto radius = v0.transRadius();
    if (radius < cut.cK0sRadiusMin.value || radius > cut.cK0sRadiusMax.value) {
      return kK0sKinematics;
    }
    if (std::abs(v0.dcav0topv()) > MaxDCAV0ToPV) {
      return kK0sKinematics;
    }
    if (k0sProperLifetime(collision, v0) > cut.cK0sProperLifetimeMax.value) {
      return kK0sTopology;
    }
    if (v0.qtarm() < cut.cK0sArmenterosQtMin.value) {
      return kK0sTopology;
    }
    if (std::abs(v0.mK0Short() - o2::constants::physics::MassK0Short) > cut.cK0sMassWindow.value) {
      return kK0sLifetime;
    }
    if (cut.cK0sCrossMassRejection.value) {
      if (std::abs(v0.mLambda() - o2::constants::physics::MassLambda) < cut.cK0sCrossMassRejectionWindow.value ||
          std::abs(v0.mAntiLambda() - o2::constants::physics::MassLambda) < cut.cK0sCrossMassRejectionWindow.value) {
        return kK0sLifetime;
      }
    }
    if (std::abs(v0.daughterTPCNSigmaPosPi()) >= cut.cK0sDaughterPiTPCNSigmaMax.value ||
        std::abs(v0.daughterTPCNSigmaNegPi()) >= cut.cK0sDaughterPiTPCNSigmaMax.value) {
      return kK0sMass;
    }
    if (v0.nCrossedRowsPos() <= cut.cK0sPosDaughterMinCrossedRows.value || v0.nCrossedRowsNeg() <= cut.cK0sNegDaughterMinCrossedRows.value) {
      return kK0sMass;
    }
    if (v0.qtarm() < cut.cK0sArmenterosAlphaCoeff.value * std::fabs(v0.alpha())) {
      return kK0sDaughters;
    }
    return kK0sArmenteros;
  }

  // Xi-K0S opening angle and the pT-dependent opening-angle window of the original task (kinCuts)
  template <typename FirstVecT, typename SecondVecT, typename MotherVecT>
  [[nodiscard]] bool passesKinematicCut(const FirstVecT& firstDaughter, const SecondVecT& secondDaughter, const MotherVecT& mother, float& alpha) const
  {
    auto firstP = std::sqrt(firstDaughter.Px() * firstDaughter.Px() + firstDaughter.Py() * firstDaughter.Py() + firstDaughter.Pz() * firstDaughter.Pz());
    auto secondP = std::sqrt(secondDaughter.Px() * secondDaughter.Px() + secondDaughter.Py() * secondDaughter.Py() + secondDaughter.Pz() * secondDaughter.Pz());
    if (firstP < SmallNumber || secondP < SmallNumber) {
      alpha = 0.f;
      return false;
    }

    auto cosAlpha = (firstDaughter.Px() * secondDaughter.Px() + firstDaughter.Py() * secondDaughter.Py() + firstDaughter.Pz() * secondDaughter.Pz()) / (firstP * secondP);
    if (cosAlpha > 1.) {
      cosAlpha = 1.;
    } else if (cosAlpha < -1.) {
      cosAlpha = -1.;
    }
    alpha = std::acos(cosAlpha);

    const auto& kinCutsPt = mK0sCuts.cKinCutsPt.value;
    const auto& kinLowerCutsAlpha = mK0sCuts.cKinLowerCutsAlpha.value;
    const auto& kinUpperCutsAlpha = mK0sCuts.cKinUpperCutsAlpha.value;

    int kinCutsSize = static_cast<int>(kinUpperCutsAlpha.size());
    if (kinCutsSize > static_cast<int>(kinLowerCutsAlpha.size())) {
      kinCutsSize = static_cast<int>(kinLowerCutsAlpha.size());
    }
    if (kinCutsSize > static_cast<int>(kinCutsPt.size()) - 1) {
      kinCutsSize = static_cast<int>(kinCutsPt.size()) - 1;
    }

    for (int i = 0; i < kinCutsSize; ++i) {
      if ((mother.Pt() > kinCutsPt[i] && mother.Pt() <= kinCutsPt[i + 1]) && (alpha < kinLowerCutsAlpha[i] || alpha > kinUpperCutsAlpha[i])) {
        return false;
      }
    }
    return true;
  }

  // Full selection stage of a mode-B track (quality, TOF requirement, PID)
  template <bool IsResoMicrotrack, Species S, typename TrackType>
  int trackSelectionStage(const TrackType& track)
  {
    return speciesStage<IsResoMicrotrack, S>(track, mSelection.trackQualityStage<IsResoMicrotrack>(track));
  }

  template <bool IsResoMicrotrack, typename TrackType>
  int pionSelectionStage(const TrackType& track)
  {
    return trackSelectionStage<IsResoMicrotrack, Species::Pion>(track);
  }

  template <bool IsResoMicrotrack, typename TrackType>
  int kaonSelectionStage(const TrackType& track)
  {
    return trackSelectionStage<IsResoMicrotrack, Species::Kaon>(track);
  }

  // Mode A: (Xi, K0S) enumeration of one collision, or of one mixed pair (Xi from collision, K0S from v0Collision).
  // The order of the original task is kept: K0S selection, Xi selection, then Xi (outer) x K0S (inner).
  // Hooks (nullptr to skip):
  //  - onK0s(v0, properLifetime, selected): every V0, in table order,
  //  - onXi(xi, selected): every cascade, in table order,
  //  - onCandidate(xi, v0, XiK0sCandidateValues): every pair of selected objects without shared daughters
  //    (the daughter-ID check is skipped in mixed events); kinematic and rapidity flags are in the values,
  //  - onExport(collision, xi, v0, values, passBits): same-event candidates at the loose or selected stage.
  template <bool IsMC, bool IsMix, typename CollisionType, typename V0CollisionType, typename CascadesType, typename V0sType,
            typename XiHook = std::nullptr_t, typename K0sHook = std::nullptr_t, typename CandidateHook = std::nullptr_t, typename ExportHook = std::nullptr_t>
  void forEachXiK0sCandidate(o2::framework::HistogramRegistry& histos, const CollisionType& collision, const V0CollisionType& v0Collision,
                             const CascadesType& cascades, const V0sType& v0s, bool computeValues,
                             XiHook onXi = nullptr, K0sHook onK0s = nullptr, CandidateHook onCandidate = nullptr, ExportHook onExport = nullptr)
  {
    constexpr bool HasXiHook = !std::is_same_v<XiHook, std::nullptr_t>;
    constexpr bool HasK0sHook = !std::is_same_v<K0sHook, std::nullptr_t>;
    constexpr bool HasCandidateHook = !std::is_same_v<CandidateHook, std::nullptr_t>;
    constexpr bool HasExportHook = !std::is_same_v<ExportHook, std::nullptr_t>;
    constexpr bool FillCutFlow = !IsMix;
    bool visitLoose = false;
    bool checkValidity = false;
    if constexpr (!IsMix) {
      visitLoose = mLooseOptions.audit || (!mLooseOptions.exportSelected && HasExportHook);
      checkValidity = visitLoose || (mLooseOptions.exportSelected && HasExportHook);
    }

    // K0S selection cache
    using V0Row = std::decay_t<decltype(v0s.begin())>;
    std::vector<CachedObject<V0Row, 2>> k0sList;
    k0sList.reserve(v0s.size());
    for (const auto& v0 : v0s) {
      const int stage = k0sSelectionStage(v0Collision, v0);
      const bool selected = stage == kK0sSelected;
      if constexpr (HasK0sHook) {
        onK0s(v0, k0sProperLifetime(v0Collision, v0), selected);
      }
      if constexpr (FillCutFlow) {
        for (int i = 0; i <= stage; ++i) {
          histos.fill(HIST("CutFlow/xiK0s/k0s"), i);
        }
      }
      if (selected || visitLoose) {
        const auto indices = v0.indices();
        k0sList.push_back({v0, {indices[0], indices[1]}, selected});
      }
    }

    // Xi selection cache
    using XiRow = std::decay_t<decltype(cascades.begin())>;
    std::vector<CachedObject<XiRow, 3>> xiList;
    xiList.reserve(cascades.size());
    for (const auto& xi : cascades) {
      const int stage = xiSelectionStage(xi);
      const bool selected = stage == kXiSelected;
      if constexpr (HasXiHook) {
        onXi(xi, selected);
      }
      if constexpr (FillCutFlow) {
        for (int i = 0; i <= stage; ++i) {
          histos.fill(HIST("CutFlow/xiK0s/xi"), i);
        }
      }
      if (selected || visitLoose) {
        const auto indices = xi.cascadeIndices();
        xiList.push_back({xi, {indices[0], indices[1], indices[2]}, selected});
      }
    }

    const bool kinCutsOn = mK0sCuts.cKinCuts.value;
    const float rapidityMax = mCandidateCuts.cfgRapidityCut.value;
    for (const auto& cachedXi : xiList) {
      const auto& xi = cachedXi.row;
      for (const auto& cachedK0s : k0sList) {
        const auto& v0 = cachedK0s.row;
        const bool objectsSelected = cachedXi.selected && cachedK0s.selected;
        bool truthMatched = false;
        if constexpr (IsMC && !IsMix) {
          truthMatched = classifyXiK0sTruth(xi, v0) == XiK0sTruth::Matched;
        }
        auto countCandidate = [&](int stage) {
          if constexpr (FillCutFlow) {
            if (objectsSelected) {
              histos.fill(HIST("CutFlow/xiK0s/candidates"), stage, 0);
              if (truthMatched) {
                histos.fill(HIST("CutFlow/xiK0s/candidates"), stage, 1);
              }
            }
          }
        };
        countCandidate(0);
        if constexpr (!IsMix) {
          if (sharesAnyDaughterId(cachedXi.daughterIds, cachedK0s.daughterIds)) {
            continue;
          }
        }
        countCandidate(1);

        // 4-vectors (construction of the original task)
        XiK0sCandidateValues values;
        values.xi = ROOT::Math::PxPyPzEVector(ROOT::Math::PtEtaPhiMVector(xi.pt(), xi.eta(), xi.phi(), xi.mXi()));
        values.k0s = ROOT::Math::PxPyPzEVector(ROOT::Math::PtEtaPhiMVector(v0.pt(), v0.eta(), v0.phi(), o2::constants::physics::MassK0Short));
        values.omega = values.xi + values.k0s;
        if (kinCutsOn || computeValues || checkValidity) {
          float alpha = 0.f;
          const bool kinCutFlag = passesKinematicCut(values.xi, values.k0s, values.omega, alpha);
          values.alpha = alpha;
          if (kinCutsOn) {
            values.passesKinCut = kinCutFlag;
          }
        }
        values.inRapidity = !(std::abs(values.omega.Rapidity()) >= rapidityMax);
        if (values.passesKinCut) {
          countCandidate(2);
          if (values.inRapidity) {
            countCandidate(3);
          }
        }

        if constexpr (!IsMix) {
          if (checkValidity) {
            exportXiK0s(histos, collision, xi, v0, values, cachedXi.selected, cachedK0s.selected, truthMatched, visitLoose, onExport);
          }
        }
        if constexpr (HasCandidateHook) {
          if (objectsSelected) {
            onCandidate(xi, v0, values);
          }
        }
      }
    }
  }

  // Mode B: (Xi, pion, kaon) enumeration of one collision, or of one mixed pair (Xi from collision, tracks from the other).
  // trackIds: ResoTrackTracks for full ResoTracks (positional source track IDs), unused for ResoMicroTracks_001 (trackId()).
  // Hooks (nullptr to skip):
  //  - onXi(xi, selected): every cascade, in table order,
  //  - onTrack(track, pionStage, kaonStage): every track, in table order (TrackStage values),
  //  - onCandidate(xi, pion, kaon, Xi1530KCandidateValues): every triplet of selected objects with five distinct
  //    source track IDs (no Xi-daughter check in mixed events) inside the Xi(1530) window; all charge patterns,
  //  - onExport(collision, xi, pion, kaon, values, passBits): same-event micro candidates at the loose or selected stage.
  template <bool IsMC, bool IsMix, bool IsResoMicrotrack, typename CollisionType, typename CascadesType, typename TracksType, typename TrackIdsType,
            typename XiHook = std::nullptr_t, typename TrackHook = std::nullptr_t, typename CandidateHook = std::nullptr_t, typename ExportHook = std::nullptr_t>
  void forEachXi1530KCandidate(o2::framework::HistogramRegistry& histos, const CollisionType& collision,
                               const CascadesType& cascades, const TracksType& tracks, const TrackIdsType& trackIds, bool computeValues,
                               XiHook onXi = nullptr, TrackHook onTrack = nullptr, CandidateHook onCandidate = nullptr, ExportHook onExport = nullptr)
  {
    constexpr bool HasXiHook = !std::is_same_v<XiHook, std::nullptr_t>;
    constexpr bool HasTrackHook = !std::is_same_v<TrackHook, std::nullptr_t>;
    constexpr bool HasCandidateHook = !std::is_same_v<CandidateHook, std::nullptr_t>;
    constexpr bool HasExportHook = !std::is_same_v<ExportHook, std::nullptr_t>;
    constexpr bool FillCutFlow = !IsMix;
    bool visitLoose = false;
    bool checkValidity = false;
    if constexpr (IsResoMicrotrack && !IsMix) {
      visitLoose = mLooseOptions.audit || (!mLooseOptions.exportSelected && HasExportHook);
      checkValidity = visitLoose || (mLooseOptions.exportSelected && HasExportHook);
    }

    // Xi selection cache
    using XiRow = std::decay_t<decltype(cascades.begin())>;
    std::vector<CachedObject<XiRow, 3>> xiList;
    xiList.reserve(cascades.size());
    for (const auto& xi : cascades) {
      const int stage = xiSelectionStage(xi);
      const bool selected = stage == kXiSelected;
      if constexpr (HasXiHook) {
        onXi(xi, selected);
      }
      if constexpr (FillCutFlow) {
        for (int i = 0; i <= stage; ++i) {
          histos.fill(HIST("CutFlow/xi1530K/xi"), i);
        }
      }
      if (selected || visitLoose) {
        const auto indices = xi.cascadeIndices();
        xiList.push_back({xi, {indices[0], indices[1], indices[2]}, selected});
      }
    }

    // Track selection cache: every track is selected once per collision, for both species
    const std::size_t nTracks = tracks.size();
    if (nTracks == 0) {
      return;
    }
    const int64_t firstIndex = tracks.begin().index();
    std::vector<TrackCacheEntry> trackCache(nTracks);
    for (const auto& track : tracks) {
      auto& entry = trackCache[getCacheIndex(track, firstIndex, nTracks)];
      const int qualityStage = mSelection.trackQualityStage<IsResoMicrotrack>(track);
      const int pionStage = speciesStage<IsResoMicrotrack, Species::Pion>(track, qualityStage);
      const int kaonStage = speciesStage<IsResoMicrotrack, Species::Kaon>(track, qualityStage);
      entry.sourceId = sourceTrackId(track, trackIds);
      entry.quality = qualityStage == TrackStage::kTrkClusters;
      entry.pionSelected = pionStage == TrackStage::kTrkPID;
      entry.kaonSelected = kaonStage == TrackStage::kTrkPID;
      if constexpr (HasTrackHook) {
        onTrack(track, pionStage, kaonStage);
      }
      if constexpr (FillCutFlow) {
        for (int i = 0; i <= pionStage; ++i) {
          histos.fill(HIST("CutFlow/xi1530K/tracks"), i, static_cast<int>(Species::Pion));
        }
        for (int i = 0; i <= kaonStage; ++i) {
          histos.fill(HIST("CutFlow/xi1530K/tracks"), i, static_cast<int>(Species::Kaon));
        }
      }
    }

    const float rapidityMax = mCandidateCuts.cfgRapidityCut.value;
    const bool useWindow = mXi1530KCuts.cXi1530UseMassWindow.value;
    const float windowCenter = mXi1530KCuts.cXi1530Mass.value;
    const float windowWidth = mXi1530KCuts.cXi1530MassWindow.value;
    for (const auto& cachedXi : xiList) {
      const auto& xi = cachedXi.row;
      const int xiSign = xi.sign() < 0 ? -1 : 1;
      const ROOT::Math::PxPyPzEVector pXi(ROOT::Math::PtEtaPhiMVector(xi.pt(), xi.eta(), xi.phi(), xi.mXi()));
      for (const auto& pion : tracks) {
        const auto& pionEntry = trackCache[getCacheIndex(pion, firstIndex, nTracks)];
        if (!pionEntry.pionSelected && !visitLoose) {
          continue;
        }
        if constexpr (!IsMix) {
          if (sharesDaughterId(cachedXi.daughterIds, pionEntry.sourceId)) {
            continue;
          }
        }
        const bool pairSelected = cachedXi.selected && pionEntry.pionSelected;
        const bool pionSignal = pion.sign() != xiSign;
        const ROOT::Math::PxPyPzEVector pPion(ROOT::Math::PtEtaPhiMVector(pion.pt(), pion.eta(), pion.phi(), o2::constants::physics::MassPionCharged));
        const ROOT::Math::PxPyPzEVector pXi1530 = pXi + pPion;
        const float massXiPi = pXi1530.M();
        const bool inWindow = !useWindow || std::abs(massXiPi - windowCenter) < windowWidth;
        if constexpr (FillCutFlow) {
          if (pairSelected) {
            countXi1530K(histos, 0, pionSignal, false);
            if (inWindow) {
              countXi1530K(histos, 1, pionSignal, false);
            }
          }
        }
        if (!(pairSelected && inWindow) && !visitLoose) {
          continue;
        }
        for (const auto& kaon : tracks) {
          if (kaon.index() == pion.index()) {
            continue;
          }
          const auto& kaonEntry = trackCache[getCacheIndex(kaon, firstIndex, nTracks)];
          if (!kaonEntry.kaonSelected && !visitLoose) {
            continue;
          }
          if (kaonEntry.sourceId == pionEntry.sourceId) {
            continue;
          }
          if constexpr (!IsMix) {
            if (sharesDaughterId(cachedXi.daughterIds, kaonEntry.sourceId)) {
              continue;
            }
          }
          const bool tripletSelected = pairSelected && inWindow && kaonEntry.kaonSelected;
          Xi1530KCandidateValues values;
          values.xi = pXi;
          values.pion = pPion;
          values.kaon = ROOT::Math::PxPyPzEVector(ROOT::Math::PtEtaPhiMVector(kaon.pt(), kaon.eta(), kaon.phi(), o2::constants::physics::MassKaonCharged));
          values.xi1530 = pXi1530;
          values.omega = pXi1530 + values.kaon;
          values.massXiPi = massXiPi;
          values.chargePattern = o2::analysis::omega2012ml::chargePattern(xiSign, pion.sign(), kaon.sign());
          values.inXi1530Window = inWindow;
          values.inRapidity = !(std::abs(values.omega.Rapidity()) >= rapidityMax);
          if (computeValues || checkValidity) {
            values.massXiK = (pXi + values.kaon).M();
            values.massPiK = (pPion + values.kaon).M();
            values.openingAngleXi1530K = ROOT::Math::VectorUtil::Angle(pXi1530, values.kaon);
            values.cosThetaStar = static_cast<float>(o2::analysis::omega2012ml::detail::cosThetaStar(
              o2::analysis::omega2012ml::detail::LorentzVector(pXi1530), o2::analysis::omega2012ml::detail::LorentzVector(values.omega)));
          }
          bool truthMatched = false;
          if constexpr (IsMC && !IsMix) {
            truthMatched = classifyXi1530KTruth(xi, pion, kaon) == Xi1530KTruth::Matched;
          }
          if constexpr (FillCutFlow) {
            if (tripletSelected) {
              const bool signal = values.chargePattern == o2::analysis::omega2012ml::kSignalPattern;
              countXi1530K(histos, 2, signal, truthMatched);
              if (values.inRapidity) {
                countXi1530K(histos, 3, signal, truthMatched);
              }
            }
          }
          if constexpr (IsResoMicrotrack && !IsMix) {
            if (checkValidity) {
              exportXi1530K(histos, collision, xi, pion, kaon, values, cachedXi.selected,
                            pionEntry.quality && kaonEntry.quality, pionEntry.pionSelected && kaonEntry.kaonSelected,
                            truthMatched, visitLoose, onExport);
            }
          }
          if constexpr (HasCandidateHook) {
            if (tripletSelected) {
              onCandidate(xi, pion, kaon, values);
            }
          }
        }
      }
    }
  }

  // Generated Omega(2012) parents of a selected reconstructed MC collision (ResoMCParents_001). The optional callback
  // receives (parent, immediate channel) for the parents inside the rapidity window. Split reconstructed collisions
  // repeat parent sets: this is not an unconditional generated denominator.
  template <typename ParentsType, typename Callback = std::nullptr_t>
  void forEachGeneratedOmega2012(o2::framework::HistogramRegistry& histos, const ParentsType& resoParents, Callback callback = nullptr)
  {
    for (const auto& part : resoParents) {
      if (std::abs(part.pdgCode()) != o2::constants::physics::Pdg::kOmega2012Minus) {
        continue;
      }
      const GeneratedChannel channel = classifyGeneratedOmega2012(part.pdgCode(), part.daughterPDG1(), part.daughterPDG2());
      histos.fill(HIST("CutFlow/generated"), 0, static_cast<int>(channel));
      if (!(std::abs(part.y()) < mCandidateCuts.cfgRapidityCut.value)) {
        continue;
      }
      histos.fill(HIST("CutFlow/generated"), 1, static_cast<int>(channel));
      if constexpr (!std::is_same_v<Callback, std::nullptr_t>) {
        callback(part, channel);
      }
    }
  }

  template <std::size_t N, std::size_t M>
  static bool sharesAnyDaughterId(const std::array<int, N>& first, const std::array<int, M>& second)
  {
    for (const auto& firstId : first) {
      for (const auto& secondId : second) {
        if (firstId == secondId) {
          return true;
        }
      }
    }
    return false;
  }

  template <std::size_t N>
  static bool sharesDaughterId(const std::array<int, N>& daughters, int64_t trackId)
  {
    for (const auto& daughterId : daughters) {
      if (static_cast<int64_t>(daughterId) == trackId) {
        return true;
      }
    }
    return false;
  }

  // Source track ID: trackId() of ResoMicroTracks_001, or the positional ResoTrackTracks row of a ResoTracks row.
  template <typename TrackType, typename TrackIdsType>
  static int64_t sourceTrackId(const TrackType& track, const TrackIdsType& trackIds)
  {
    if constexpr (requires { track.trackId(); }) {
      return static_cast<int64_t>(track.trackId());
    } else {
      const auto rowIndex = track.globalIndex();
      if (rowIndex >= 0 && rowIndex < static_cast<int64_t>(trackIds.size())) {
        return static_cast<int64_t>(trackIds.rawIteratorAt(rowIndex).trackId());
      }
      return static_cast<int64_t>(rowIndex);
    }
  }

 private:
  template <typename Row, std::size_t N>
  struct CachedObject {
    Row row;
    std::array<int, N> daughterIds{};
    bool selected = false;
  };

  struct TrackCacheEntry {
    int64_t sourceId = -1;
    bool quality = false;
    bool pionSelected = false;
    bool kaonSelected = false;
  };

  template <bool IsResoMicrotrack, Species S, typename TrackType>
  int speciesStage(const TrackType& track, int qualityStage)
  {
    if (qualityStage < TrackStage::kTrkClusters) {
      return qualityStage;
    }
    constexpr int SpeciesIndex = static_cast<int>(S);
    if (!mSelection.passesTOFRequired(SpeciesIndex, track)) {
      return TrackStage::kTrkClusters;
    }
    const bool hasTOF = track.hasTOF();
    const double tpcNSigma = (S == Species::Pion) ? track.tpcNSigmaPi() : track.tpcNSigmaKa();
    double tofNSigma = std::numeric_limits<double>::quiet_NaN(); // TOF value is only valid with hasTOF
    if (hasTOF) {
      tofNSigma = (S == Species::Pion) ? track.tofNSigmaPi() : track.tofNSigmaKa();
    }
    if (!mSelection.passesPID<IsResoMicrotrack>(SpeciesIndex, track.pt(), hasTOF, tpcNSigma, tofNSigma)) {
      return TrackStage::kTrkTOFRequired;
    }
    return TrackStage::kTrkPID;
  }

  template <typename TrackType>
  static std::size_t getCacheIndex(const TrackType& track, int64_t firstIndex, std::size_t size)
  {
    const int64_t index = static_cast<int64_t>(track.index()) - firstIndex;
    if (index < 0 || index >= static_cast<int64_t>(size)) {
      LOG(fatal) << "Track index " << track.index() << " is outside the selection cache [" << firstIndex << ", " << firstIndex + static_cast<int64_t>(size) << ")";
    }
    return static_cast<std::size_t>(index);
  }

  static void countXi1530K(o2::framework::HistogramRegistry& histos, int stage, bool signalPattern, bool truthMatched)
  {
    histos.fill(HIST("CutFlow/xi1530K/candidates"), stage, 0);
    if (signalPattern) {
      histos.fill(HIST("CutFlow/xi1530K/candidates"), stage, 1);
    }
    if (truthMatched) {
      histos.fill(HIST("CutFlow/xi1530K/candidates"), stage, 2);
    }
  }

  // Loose / selected stage of a mode-A same-event candidate
  template <typename CollisionType, typename XiType, typename V0Type, typename ExportHook>
  void exportXiK0s(o2::framework::HistogramRegistry& histos, const CollisionType& collision, const XiType& xi, const V0Type& v0,
                   XiK0sCandidateValues const& values, bool xiSelected, bool k0sSelected, bool truthMatched,
                   bool visitLoose, ExportHook onExport)
  {
    constexpr bool HasExportHook = !std::is_same_v<ExportHook, std::nullptr_t>;
    namespace ml = o2::analysis::omega2012ml;
    const auto canonical = ml::canonicalizeXiK0s(ml::makeCascadeSnapshot(collision, xi), ml::makeV0Snapshot(collision, v0));
    const bool valid = canonical.status == ml::BuildStatus::Ok &&
                       std::isfinite(values.omega.M()) && std::isfinite(values.omega.Pt()) && std::isfinite(values.omega.Rapidity()) &&
                       ml::buildXiK0sFeatures(canonical.candidate).status == ml::BuildStatus::Ok;
    if (!valid) {
      return;
    }
    auto countMl = [&](int stage) {
      if (mLooseOptions.audit) {
        histos.fill(HIST("ML/xiK0s/looseCutflow"), stage, 0);
        if (truthMatched) {
          histos.fill(HIST("ML/xiK0s/looseCutflow"), stage, 1);
        }
      }
    };
    countMl(0);
    if (!values.inRapidity) {
      return;
    }
    const bool selected = xiSelected && k0sSelected && values.passesKinCut;
    if (visitLoose) {
      countMl(1);
      if (mLooseOptions.audit) {
        histos.fill(HIST("ML/xiK0s/looseMassPtActivity"), values.omega.M(), values.omega.Pt(), collision.cent());
      }
      uint16_t passBits = kXiK0sPassLoose;
      if (xiSelected) {
        passBits |= kXiK0sPassXi;
        countMl(2);
        if (k0sSelected) {
          passBits |= kXiK0sPassK0s;
          countMl(3);
          if (values.passesKinCut) {
            passBits |= kXiK0sPassKinematic;
            countMl(4);
          }
        }
      }
      if constexpr (HasExportHook) {
        if (!mLooseOptions.exportSelected) {
          onExport(collision, xi, v0, values, passBits);
        }
      }
    }
    if constexpr (HasExportHook) {
      if (mLooseOptions.exportSelected && selected) {
        onExport(collision, xi, v0, values, XiK0sPassBitsSelected);
      }
    }
  }

  // Loose / selected stage of a mode-B same-event micro candidate
  template <typename CollisionType, typename XiType, typename TrackType, typename ExportHook>
  void exportXi1530K(o2::framework::HistogramRegistry& histos, const CollisionType& collision, const XiType& xi,
                     const TrackType& pion, const TrackType& kaon, Xi1530KCandidateValues const& values,
                     bool xiSelected, bool qualityBoth, bool pidBoth, bool truthMatched, bool visitLoose, ExportHook onExport)
  {
    constexpr bool HasExportHook = !std::is_same_v<ExportHook, std::nullptr_t>;
    namespace ml = o2::analysis::omega2012ml;
    const auto canonical = ml::canonicalizeXi1530K(ml::makeCascadeSnapshot(collision, xi), ml::makeTrackSnapshot(pion), ml::makeTrackSnapshot(kaon));
    const bool valid = canonical.status == ml::BuildStatus::Ok &&
                       std::isfinite(values.omega.M()) && std::isfinite(values.omega.Pt()) && std::isfinite(values.omega.Rapidity()) &&
                       ml::buildXi1530KFeatures(canonical.candidate).status == ml::BuildStatus::Ok;
    if (!valid) {
      return;
    }
    const int category = values.chargePattern == ml::kSignalPattern ? 0 : 1;
    auto countMl = [&](int stage) {
      if (mLooseOptions.audit) {
        histos.fill(HIST("ML/xi1530K/looseCutflow"), stage, category);
        if (truthMatched) {
          histos.fill(HIST("ML/xi1530K/looseCutflow"), stage, 2);
        }
      }
    };
    countMl(0);
    if (!values.inRapidity) {
      return;
    }
    const bool selected = xiSelected && qualityBoth && pidBoth && values.inXi1530Window;
    if (visitLoose) {
      countMl(1);
      if (mLooseOptions.audit && category == 0) {
        histos.fill(HIST("ML/xi1530K/looseMassPtActivity"), values.omega.M(), values.omega.Pt(), collision.cent());
      }
      uint16_t passBits = kXi1530KPassLoose;
      if (xiSelected) {
        passBits |= kXi1530KPassXi;
        countMl(2);
        if (qualityBoth) {
          passBits |= kXi1530KPassQuality;
          countMl(3);
          if (pidBoth) {
            passBits |= kXi1530KPassPID;
            countMl(4);
            if (values.inXi1530Window) {
              passBits |= kXi1530KPassWindow;
              countMl(5);
            }
          }
        }
      }
      if constexpr (HasExportHook) {
        if (!mLooseOptions.exportSelected) {
          onExport(collision, xi, pion, kaon, values, passBits);
        }
      }
    }
    if constexpr (HasExportHook) {
      if (mLooseOptions.exportSelected && selected) {
        onExport(collision, xi, pion, kaon, values, Xi1530KPassBitsSelected);
      }
    }
  }

  // Cut-flow instrumentation of the selection, per decay mode; the output histograms belong to the tasks.
  void registerHistograms(o2::framework::HistogramRegistry& histos, ProcessModes const& modes)
  {
    using o2::framework::HistType;
    const std::array<const char*, kXiNStages> xiLabels{"input", "|#eta|, p_{T}", "DCA to PV", "#Lambda topology", "cascade topology", "#Xi mass"};
    if (modes.xiK0s) {
      auto xiFlow = histos.add<TH1>("CutFlow/xiK0s/xi", "XiK0s mode: cascades, once per selected collision;stage;cascades", HistType::kTH1D, {{kXiNStages, -0.5, static_cast<double>(kXiNStages) - 0.5}});
      for (std::size_t i = 0; i < xiLabels.size(); ++i) {
        xiFlow->GetXaxis()->SetBinLabel(i + 1, xiLabels[i]);
      }
      auto k0sFlow = histos.add<TH1>("CutFlow/xiK0s/k0s", "XiK0s mode: V0s, once per selected collision;stage;V0s", HistType::kTH1D, {{kK0sNStages, -0.5, static_cast<double>(kK0sNStages) - 0.5}});
      const std::array<const char*, kK0sNStages> k0sLabels{"input", "|#eta|, p_{T}", "topology", "lifetime, min q_{T}", "mass, #Lambda rejection", "daughters", "Armenteros"};
      for (std::size_t i = 0; i < k0sLabels.size(); ++i) {
        k0sFlow->GetXaxis()->SetBinLabel(i + 1, k0sLabels[i]);
      }
      auto candidateFlow = histos.add<TH2>("CutFlow/xiK0s/candidates", "XiK0s mode: selected Xi x selected K0S, same event;stage;category", HistType::kTH2D,
                                           {{NXiK0sCandidateStages, -0.5, NXiK0sCandidateStages - 0.5}, {2, -0.5, 1.5}});
      const std::array<const char*, NXiK0sCandidateStages> candidateLabels{"selected pairs", "distinct daughters", "kinematic cut", "rapidity"};
      for (std::size_t i = 0; i < candidateLabels.size(); ++i) {
        candidateFlow->GetXaxis()->SetBinLabel(i + 1, candidateLabels[i]);
      }
      candidateFlow->GetYaxis()->SetBinLabel(1, "all");
      candidateFlow->GetYaxis()->SetBinLabel(2, "XiK0s truth");
      if (mLooseOptions.audit) {
        auto flow = histos.add<TH2>("ML/xiK0s/looseCutflow", "XiK0s mode: valid canonical candidates;stage;category", HistType::kTH2D, {{5, -0.5, 4.5}, {2, -0.5, 1.5}});
        const std::array<const char*, 5> labels{"valid canonical", "loose acceptance", "Xi selection", "K0S selection", "kinematic cut"};
        for (std::size_t i = 0; i < labels.size(); ++i) {
          flow->GetXaxis()->SetBinLabel(i + 1, labels[i]);
        }
        flow->GetYaxis()->SetBinLabel(1, "all");
        flow->GetYaxis()->SetBinLabel(2, "XiK0s truth");
        histos.add("ML/xiK0s/looseMassPtActivity", "XiK0s mode loose candidates;mass (GeV/c^{2});pT (GeV/c);centrality", HistType::kTH3D,
                   {{300, 1.6, 2.8}, {{0., 0.5, 1., 2., 3., 5., 8., 15., 30., 100.}, "pT"}, {{0., 10., 30., 50., 70., 100., 110.}, "centrality"}});
      }
    }
    if (modes.xi1530K) {
      auto xiFlow = histos.add<TH1>("CutFlow/xi1530K/xi", "Xi1530K mode: cascades, once per selected collision;stage;cascades", HistType::kTH1D, {{kXiNStages, -0.5, static_cast<double>(kXiNStages) - 0.5}});
      for (std::size_t i = 0; i < xiLabels.size(); ++i) {
        xiFlow->GetXaxis()->SetBinLabel(i + 1, xiLabels[i]);
      }
      constexpr int NTrackStages = TrackStage::kTrkNStages;
      auto trackFlow = histos.add<TH2>("CutFlow/xi1530K/tracks", "Xi1530K mode: tracks, once per selected collision;stage;species", HistType::kTH2D,
                                       {{NTrackStages, -0.5, NTrackStages - 0.5}, {2, -0.5, 1.5}});
      const std::array<const char*, NTrackStages> trackLabels{"input", "pT", "eta", "DCAxy", "DCAz", "track flags", "clusters / crossed rows", "TOF required", "PID"};
      for (std::size_t i = 0; i < trackLabels.size(); ++i) {
        trackFlow->GetXaxis()->SetBinLabel(i + 1, trackLabels[i]);
      }
      trackFlow->GetYaxis()->SetBinLabel(1, "pion");
      trackFlow->GetYaxis()->SetBinLabel(2, "kaon");
      auto candidateFlow = histos.add<TH2>("CutFlow/xi1530K/candidates", "Xi1530K mode: same event;stage;category", HistType::kTH2D,
                                           {{NXi1530KCandidateStages, -0.5, NXi1530KCandidateStages - 0.5}, {3, -0.5, 2.5}});
      const std::array<const char*, NXi1530KCandidateStages> candidateLabels{"selected Xi #pi pairs", "Xi(1530) window", "selected Xi #pi K triplets", "rapidity"};
      for (std::size_t i = 0; i < candidateLabels.size(); ++i) {
        candidateFlow->GetXaxis()->SetBinLabel(i + 1, candidateLabels[i]);
      }
      candidateFlow->GetYaxis()->SetBinLabel(1, "all");
      candidateFlow->GetYaxis()->SetBinLabel(2, "signal charge pattern");
      candidateFlow->GetYaxis()->SetBinLabel(3, "Xi1530K truth");
      if (mLooseOptions.audit) {
        auto flow = histos.add<TH2>("ML/xi1530K/looseCutflow", "Xi1530K mode: valid canonical candidates;stage;category", HistType::kTH2D, {{6, -0.5, 5.5}, {3, -0.5, 2.5}});
        const std::array<const char*, 6> labels{"valid canonical", "loose acceptance", "Xi selection", "track quality", "TOF + PID", "Xi(1530) window"};
        for (std::size_t i = 0; i < labels.size(); ++i) {
          flow->GetXaxis()->SetBinLabel(i + 1, labels[i]);
        }
        flow->GetYaxis()->SetBinLabel(1, "signal charge pattern");
        flow->GetYaxis()->SetBinLabel(2, "wrong-sign control");
        flow->GetYaxis()->SetBinLabel(3, "Xi1530K truth");
        histos.add("ML/xi1530K/looseMassPtActivity", "Xi1530K mode loose signal-pattern candidates;mass (GeV/c^{2});pT (GeV/c);centrality", HistType::kTH3D,
                   {{300, 1.8, 3.0}, {{0., 0.5, 1., 2., 3., 5., 8., 15., 30., 100.}, "pT"}, {{0., 10., 30., 50., 70., 100., 110.}, "centrality"}});
      }
    }
    if (modes.mcGen) {
      auto generated = histos.add<TH2>("CutFlow/generated", "Generated Omega(2012) parent rows of selected reconstructed events;stage;immediate channel", HistType::kTH2D,
                                       {{2, -0.5, 1.5}, {NGeneratedChannels, -0.5, NGeneratedChannels - 0.5}});
      generated->GetXaxis()->SetBinLabel(1, "all parent rows");
      generated->GetXaxis()->SetBinLabel(2, "rapidity window");
      generated->GetYaxis()->SetBinLabel(1, "other / unresolved");
      generated->GetYaxis()->SetBinLabel(2, "Xi K0S");
      generated->GetYaxis()->SetBinLabel(3, "Xi(1530) K");
    }
  }

  ResoAnalysisSelectionCore mSelection;
  XiCuts mXiCuts;
  K0sCuts mK0sCuts;
  Xi1530KCuts mXi1530KCuts;
  CandidateCuts mCandidateCuts;
  LooseStageOptions mLooseOptions;
};

} // namespace o2::analysis::omega2012

#endif // PWGLF_CORE_OMEGA2012ANALYSISCORE_H_
