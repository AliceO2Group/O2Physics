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
/// \file K1AnalysisMicroCore.h
/// \brief Shared selection, truth classification and candidate loop of the K1(1270) resonance tasks
/// \author Su-Jeong Ji <su-jeong.ji@cern.ch>, Bong-Hwi Lim <bong-hwi.lim@cern.ch>
///

#ifndef PWGLF_CORE_K1ANALYSISMICROCORE_H_
#define PWGLF_CORE_K1ANALYSISMICROCORE_H_

#include "PWGLF/Core/K1MlFeatures.h"
#include "PWGLF/DataModel/LFResonanceTables.h"

#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/Logger.h>

#include <Math/GenVector/VectorUtil.h>
#include <Math/Vector4D.h> // IWYU pragma: keep (do not replace with Math/Vector4Dfwd.h)
#include <TH1.h>
#include <TH2.h>
#include <TPDGCode.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <set>
#include <type_traits>
#include <vector>

namespace o2::analysis::k1micro
{

enum class K1TruthChannel {
  None = 0,
  RhoK = 1,
  KStarPi = 2
};

enum BinAnti : unsigned int {
  kNormal = 0,
  kAnti,
  kNAEnd
};

enum BinType : unsigned int {
  kK1P = 0,
  kK1N,
  kK1P_Mix,
  kK1N_Mix,
  kK1P_GenINEL10,
  kK1N_GenINEL10,
  kK1P_GenINELgt10,
  kK1N_GenINELgt10,
  kK1P_GenTrig10,
  kK1N_GenTrig10,
  kK1P_GenEvtSel,
  kK1N_GenEvtSel,
  kK1P_Rec,
  kK1N_Rec,
  kTYEnd
};

enum class Species : int {
  Pion = 0,
  Kaon = 1
};

// Last stage passed by a track; the cut-flow histogram is filled directly from this value.
enum TrackStage : int {
  kTrkInput = 0,
  kTrkPt,
  kTrkEta,
  kTrkDCAxy,
  kTrkDCAz,
  kTrkFlags,
  kTrkClusters,
  kTrkTOFRequired,
  kTrkPID,
  kTrkNStages
};

enum class QAFolder {
  Before, // QA/*: before the candidate cuts
  After,  // QAcut/*: after the candidate cuts
  MC      // QAMC/*: matched K1 truth candidates
};

// Cumulative selection bits of an unlike-sign candidate handed to the candidate callback.
enum CandidatePassBit : uint16_t {
  kPassLoose = 1,     // valid canonical candidate inside the K1 rapidity window
  kPassQuality = 2,   // track quality of all three tracks
  kPassPID = 4,       // TOF requirement and PID of all three tracks
  kPassPair = 8,      // pion-pair pT and secondary mass window
  kPassCandidate = 16 // candidate cuts
};
inline constexpr uint16_t PassBitsSelected = kPassLoose | kPassQuality | kPassPID | kPassPair | kPassCandidate;

// Resolved PID cut of one species at a given pT.
struct PIDCut {
  double tpcMax = 0.;
  double tofMax = 0.;
  double combined = 0.;
  bool tofRequired = false;
};

inline constexpr float DisabledCut = -999.f;  // an optional cut with this value is off and not evaluated
inline constexpr double MassRho770 = 0.77526; // PDG 2024, not available in o2::constants::physics
inline constexpr double DCAGridStep = 0.025;  // v001 micro DCA encoding, lower-inclusive bins up to DCAGridMax
inline constexpr double DCAGridMax = 0.15;
inline constexpr double PIDGridStart = 2.0; // v001 micro nSigma encoding: 0.25 bins in [2.0, 3.5]
inline constexpr double PIDGridStep = 0.25;
inline constexpr double PIDGridMax = 3.5;
inline constexpr double GridTolerance = 1e-4;
inline constexpr std::size_t MinPtBinEdges = 2;     // a pT dependent PID table needs at least one bin
inline constexpr float ProducerDCAPtP0 = 0.004f;    // resonanceModuleInitializer cfgTightDCAOffset default
inline constexpr float ProducerDCAPtCoeff = 0.013f; // resonanceModuleInitializer cfgTightDCAPtCoefficient default
inline constexpr float ProducerDCAPtPower = 1.f;    // resonanceModuleInitializer cfgTightDCAPtPower default
inline constexpr float ConfigTolerance = 1e-6f;
inline constexpr int NCandidateStages = 12;
inline constexpr int NTruthChannels = 3;

// Configurable groups without prefix: the JSON keys are the plain configurable names.

/// Event selection
struct EventCuts : o2::framework::ConfigurableGroup {
  o2::framework::Configurable<bool> cRecoINELgt0{"cRecoINELgt0", false, "Apply reconstructed INEL>0 selection"};
  o2::framework::Configurable<bool> cMCINELgt0{"cMCINELgt0", false, "Require generator INEL>0 in MC processes"};
  o2::framework::Configurable<bool> cMCVtxIn10{"cMCVtxIn10", false, "Require generator |vz| < 10 cm in MC processes"};
};

/// Track selections (common for pion and kaon, -999 switches an optional cut off)
struct TrackCuts : o2::framework::ConfigurableGroup {
  o2::framework::Configurable<double> cMinPtcut{"cMinPtcut", 0.15, "Track minium pt cut"};
  o2::framework::Configurable<float> cMaxEtacut{"cMaxEtacut", -999.f, "Track maximum |eta| cut (-999: off)"};
  // DCAr to PV
  o2::framework::Configurable<double> cMaxDCArToPVcut{"cMaxDCArToPVcut", 0.1, "Track DCAr cut to PV Maximum"};
  // DCAz to PV
  o2::framework::Configurable<double> cMaxDCAzToPVcut{"cMaxDCAzToPVcut", 0.1, "Track DCAz cut to PV Maximum"};
  o2::framework::Configurable<double> cMinDCAzToPVcut{"cMinDCAzToPVcut", 0.0, "Track DCAz cut to PV Minimum"};
  o2::framework::Configurable<bool> cfgUsePtDepDCA{"cfgUsePtDepDCA", false, "Use pT dependent DCA cut instead of the fixed maximum"};
  o2::framework::Configurable<float> cDCAToPVByPtP0{"cDCAToPVByPtP0", 0.004f, "pT dependent DCA cut = P0 + coefficient / pT^power (cm)"};
  o2::framework::Configurable<float> cDCAToPVByPtCoeff{"cDCAToPVByPtCoeff", 0.013f, "Coefficient in the pT dependent DCA cut"};
  o2::framework::Configurable<float> cDCAToPVByPtPower{"cDCAToPVByPtPower", 1.f, "Power in the pT dependent DCA cut"};
  o2::framework::Configurable<bool> cfgPrimaryTrack{"cfgPrimaryTrack", true, "Primary track selection"};                    // kGoldenChi2 | kDCAxy | kDCAz
  o2::framework::Configurable<bool> cfgGlobalWoDCATrack{"cfgGlobalWoDCATrack", true, "Global track selection without DCA"}; // kQualityTracks (kTrackType | kTPCNCls | kTPCCrossedRows | kTPCCrossedRowsOverNCls | kTPCChi2NDF | kTPCRefit | kITSNCls | kITSChi2NDF | kITSRefit | kITSHits) | kInAcceptanceTracks (kPtRange | kEtaRange)
  o2::framework::Configurable<bool> cfgGlobalTrack{"cfgGlobalTrack", false, "Global track selection"};                      // kGoldenChi2 | kDCAxy | kDCAz
  o2::framework::Configurable<bool> cfgPVContributor{"cfgPVContributor", false, "PV contributor track selection"};          // PV Contriuibutor
  o2::framework::Configurable<bool> cfgUseTPCRefit{"cfgUseTPCRefit", false, "Require TPC Refit"};
  o2::framework::Configurable<bool> cfgUseITSRefit{"cfgUseITSRefit", false, "Require ITS Refit"};
  o2::framework::Configurable<int> cfgTPCcluster{"cfgTPCcluster", 0, "Number of TPC cluster (found clusters, ResoTracks only)"};
  o2::framework::Configurable<int> cfgTPCCrossedRowsMin{"cfgTPCCrossedRowsMin", 0, "Minimum number of TPC crossed rows"};
  o2::framework::Configurable<int> cfgITSNClsMin{"cfgITSNClsMin", 0, "Minimum number of ITS clusters (ResoMicroTracks only)"};
  o2::framework::Configurable<bool> cfgHasTOF{"cfgHasTOF", false, "Require TOF"};
};

/// Pion PID selection
struct PionPidCuts : o2::framework::ConfigurableGroup {
  o2::framework::Configurable<double> cMaxTPCnSigmaPion{"cMaxTPCnSigmaPion", 3.0, "TPC nSigma cut for Pion (-999: off)"};    // TPC
  o2::framework::Configurable<double> cMaxTOFnSigmaPion{"cMaxTOFnSigmaPion", 3.0, "TOF nSigma cut for Pion (-999: off)"};    // TOF
  o2::framework::Configurable<double> nsigmaCutCombinedPion{"nsigmaCutCombinedPion", -999, "Combined nSigma cut for Pion"};  // Combined
  o2::framework::Configurable<bool> cUseOnlyTOFTrackPi{"cUseOnlyTOFTrackPi", false, "Use only TOF track for PID selection"}; // Use only TOF track for Pion PID selection
  o2::framework::Configurable<bool> cPionUsePtDepPID{"cPionUsePtDepPID", false, "Use pT-dependent PID cuts for pion"};
  o2::framework::Configurable<std::vector<float>> cPionPIDPtBins{"cPionPIDPtBins", {0.0f, 0.5f, 0.8f, 2.0f, 999.0f}, "pT bin edges for pion PID cuts"};
  o2::framework::Configurable<std::vector<float>> cPionTPCNSigmaCuts{"cPionTPCNSigmaCuts", {3.0f, 3.0f, 2.0f, 2.0f}, "TPC NSigma cuts per pT bin (pion)"};
  o2::framework::Configurable<std::vector<float>> cPionTOFNSigmaCuts{"cPionTOFNSigmaCuts", {3.0f, 3.0f, 3.0f, 3.0f}, "TOF NSigma cuts per pT bin (pion)"};
  o2::framework::Configurable<std::vector<int>> cPionTOFRequired{"cPionTOFRequired", {0, 0, 1, 1}, "Require TOF per pT bin (pion)"};
};

/// Kaon PID selection
struct KaonPidCuts : o2::framework::ConfigurableGroup {
  o2::framework::Configurable<double> cMaxTPCnSigmaKaon{"cMaxTPCnSigmaKaon", 3.0, "TPC nSigma cut for Kaon (-999: off)"};    // TPC
  o2::framework::Configurable<double> cMaxTOFnSigmaKaon{"cMaxTOFnSigmaKaon", 3.0, "TOF nSigma cut for Kaon (-999: off)"};    // TOF
  o2::framework::Configurable<double> nsigmaCutCombinedKaon{"nsigmaCutCombinedKaon", -999, "Combined nSigma cut for Kaon"};  // Combined
  o2::framework::Configurable<bool> cUseOnlyTOFTrackKa{"cUseOnlyTOFTrackKa", false, "Use only TOF track for PID selection"}; // Use only TOF track for Kaon PID selection
  o2::framework::Configurable<bool> cKaonUsePtDepPID{"cKaonUsePtDepPID", false, "Use pT-dependent PID cuts for kaon"};
  o2::framework::Configurable<std::vector<float>> cKaonPIDPtBins{"cKaonPIDPtBins", {0.0f, 0.5f, 0.8f, 2.0f, 999.0f}, "pT bin edges for kaon PID cuts"};
  o2::framework::Configurable<std::vector<float>> cKaonTPCNSigmaCuts{"cKaonTPCNSigmaCuts", {3.0f, 3.0f, 2.0f, 2.0f}, "TPC NSigma cuts per pT bin (kaon)"};
  o2::framework::Configurable<std::vector<float>> cKaonTOFNSigmaCuts{"cKaonTOFNSigmaCuts", {3.0f, 3.0f, 3.0f, 3.0f}, "TOF NSigma cuts per pT bin (kaon)"};
  o2::framework::Configurable<std::vector<int>> cKaonTOFRequired{"cKaonTOFRequired", {0, 0, 1, 1}, "Require TOF per pT bin (kaon)"};
};

/// Secondary selection (-999 switches a cut off; the values it needs are then not computed)
struct SecondaryCuts : o2::framework::ConfigurableGroup {
  o2::framework::Configurable<double> cMinSecondaryPtCut{"cMinSecondaryPtCut", 0.5, "Min pT cut for secondary selection"};
  o2::framework::Configurable<bool> cfgModeK892orRho{"cfgModeK892orRho", false, "Secondary scenario for K892 (true) or Rho (false)"};
  o2::framework::Configurable<double> cSecondaryMasswindow{"cSecondaryMasswindow", -999, "Secondary inv mass selection window"};
  o2::framework::Configurable<double> cMinAnotherSecondaryMassCut{"cMinAnotherSecondaryMassCut", -999, "Min inv. mass selection of another secondary scenario"};
  o2::framework::Configurable<double> cMaxAnotherSecondaryMassCut{"cMaxAnotherSecondaryMassCut", -999, "MAx inv. mass selection of another secondary scenario"};
  o2::framework::Configurable<double> cMinPiKaMassCut{"cMinPiKaMassCut", -999, "bPion-Kaon pair inv mass selection minimum"};
  o2::framework::Configurable<double> cMaxPiKaMassCut{"cMaxPiKaMassCut", -999, "bPion-Kaon pair inv mass selection maximum"};
  o2::framework::Configurable<double> cMinAngle{"cMinAngle", -999, "Minimum angle between the secondary resonance and the bachelor"};
  o2::framework::Configurable<double> cMaxAngle{"cMaxAngle", -999, "Maximum angle between the secondary resonance and the bachelor"};
  o2::framework::Configurable<double> cMinPairAsym{"cMinPairAsym", -999, "Minimum pair asymmetry"};
  o2::framework::Configurable<double> cMaxPairAsym{"cMaxPairAsym", -999, "Maximum pair asymmetry"};
};

/// Common TOF switch and K1 selection
struct CandidateCuts : o2::framework::ConfigurableGroup {
  o2::framework::Configurable<bool> cByPassTOF{"cByPassTOF", false, "Bypass the TOF nSigma selection"};
  o2::framework::Configurable<double> cK1MaxRap{"cK1MaxRap", 0.5, "K1 maximum rapidity"};
  o2::framework::Configurable<double> cK1MinRap{"cK1MinRap", -0.5, "K1 minimum rapidity"};
};

/// Histogram binning, QA and debug output
struct HistogramOptions : o2::framework::ConfigurableGroup {
  o2::framework::Configurable<int> cNbinsDiv{"cNbinsDiv", 1, "Integer to divide the number of bins"};
  o2::framework::Configurable<bool> additionalQAplots{"additionalQAplots", true, "Additional QA plots"};
  o2::framework::Configurable<int> cfgTruthDebug{"cfgTruthDebug", 0, "Maximum logged matched candidates per truth channel"};
};

/// Process functions enabled in the task; they decide the configuration checks and the registered histograms.
struct ProcessModes {
  bool microTracks = false; // any process function reading micro tracks (quantised DCA and nSigma)
  bool mcReco = false;      // reconstructed MC with full or micro tracks
  bool mcRecoMicro = false; // reconstructed MC with micro tracks
  bool mcGen = false;       // generated K1 parents in selected reconstructed events
};

/// Loose-stage traversal of unlike-sign micro candidates for the candidate callback.
/// With the defaults and without a callback, the candidate loop applies only the conventional selection.
struct LooseStageOptions {
  bool audit = false;          // fill ML/looseCutflow and ML/looseMassPtActivity
  bool exportSelected = false; // hand candidates to the callback at the selected stage instead of the loose stage
};

// A cut is on unless it carries the disabled value (tolerant to the float parsing of the JSON value).
inline bool isCutEnabled(float value)
{
  return value > DisabledCut + 1.f;
}

// v001 micro values are lower-inclusive bin edges: a maximum cut on the grid keeps bins below it.
inline bool passesBinnedMax(double decoded, double cut)
{
  return decoded < cut - o2::constants::math::Epsilon;
}

// Minimum cut on the grid keeps the bin starting at the cut.
inline bool passesBinnedMin(double decoded, double cut)
{
  return decoded >= cut - o2::constants::math::Epsilon;
}

template <bool IsResoMicrotrack>
bool passesMax(double value, double cut)
{
  if constexpr (IsResoMicrotrack) {
    return passesBinnedMax(value, cut);
  } else {
    return value < cut;
  }
}

inline bool isInRange(double value, double minimum, double maximum)
{
  if (isCutEnabled(minimum) && value < minimum) {
    return false;
  }
  if (isCutEnabled(maximum) && value > maximum) {
    return false;
  }
  return true;
}

inline bool isInWindow(double value, double center, double width)
{
  return std::abs(value - center) < width;
}

// Preserve pT-bin membership [low, high).
inline int getPtBinIndex(float pt, const std::vector<float>& ptBins)
{
  for (std::size_t i = 1; i < ptBins.size(); ++i) {
    if (pt >= ptBins[i - 1] && pt < ptBins[i]) {
      return static_cast<int>(i - 1);
    }
  }
  return -1;
}

// Truth classification from the immediate mothers and the sibling IDs of the reconstructed daughters.
template <typename Track>
bool hasSibling(const Track& directDaughter, int resonanceId)
{
  if (resonanceId < 0) {
    return false;
  }
  const auto siblings = directDaughter.siblingIds();
  return siblings[0] == resonanceId || siblings[1] == resonanceId;
}

template <typename Track, typename Kaon>
bool matchesKStarPi(const Track& resonancePion, const Track& directPion, const Kaon& kaon)
{
  const int charge = kaon.pdgCode() > 0 ? 1 : -1;
  if (resonancePion.motherId() != kaon.motherId() ||
      resonancePion.motherId() == directPion.motherId()) {
    return false;
  }
  if (resonancePion.motherPDG() != charge * o2::constants::physics::Pdg::kK0Star892 || kaon.motherPDG() != charge * o2::constants::physics::Pdg::kK0Star892) {
    return false;
  }
  if (resonancePion.pdgCode() != -charge * kPiPlus || directPion.pdgCode() != charge * kPiPlus ||
      directPion.motherPDG() != charge * o2::constants::physics::Pdg::kK1_1270Plus) {
    return false;
  }
  return hasSibling(directPion, kaon.motherId());
}

template <typename Track, typename Kaon>
K1TruthChannel classifyK1Truth(const Track& pion1, const Track& pion2, const Kaon& kaon)
{
  if (std::abs(pion1.pdgCode()) != kPiPlus || std::abs(pion2.pdgCode()) != kPiPlus ||
      std::abs(kaon.pdgCode()) != kKPlus) {
    return K1TruthChannel::None;
  }
  if (pion1.motherId() < 0 || pion2.motherId() < 0 || kaon.motherId() < 0) {
    return K1TruthChannel::None;
  }
  const int charge = kaon.pdgCode() > 0 ? 1 : -1;
  const bool rhoPions = pion1.motherId() == pion2.motherId() &&
                        pion1.motherPDG() == kRho770_0 && pion2.motherPDG() == kRho770_0 &&
                        pion1.pdgCode() == -pion2.pdgCode();
  if (rhoPions && kaon.motherPDG() == charge * o2::constants::physics::Pdg::kK1_1270Plus &&
      kaon.motherId() != pion1.motherId() && hasSibling(kaon, pion1.motherId())) {
    return K1TruthChannel::RhoK;
  }
  if (matchesKStarPi(pion1, pion2, kaon) || matchesKStarPi(pion2, pion1, kaon)) {
    return K1TruthChannel::KStarPi;
  }
  return K1TruthChannel::None;
}

// Immediate decay channel of a generated K1 from the PDG codes of its two daughters.
inline K1TruthChannel classifyGeneratedK1(int charge, int daughter1, int daughter2)
{
  if ((daughter1 == kRho770_0 && daughter2 == charge * kKPlus) ||
      (daughter2 == kRho770_0 && daughter1 == charge * kKPlus)) {
    return K1TruthChannel::RhoK;
  }
  if ((daughter1 == charge * o2::constants::physics::Pdg::kK0Star892 && daughter2 == charge * kPiPlus) ||
      (daughter2 == charge * o2::constants::physics::Pdg::kK0Star892 && daughter1 == charge * kPiPlus)) {
    return K1TruthChannel::KStarPi;
  }
  return K1TruthChannel::None;
}

// Track QA of a pion; isPrimary selects the trkppion (first) or trkspion (second) histograms
template <QAFolder Folder, typename TrackType>
void fillPionQA(o2::framework::HistogramRegistry& histos, const TrackType& track, bool isPrimary)
{
  const bool hasTOF = track.hasTOF();
  if (isPrimary) {
    if constexpr (Folder == QAFolder::Before) {
      histos.fill(HIST("QA/trkppionTPCPID"), track.pt(), track.tpcNSigmaPi());
      if (hasTOF) {
        histos.fill(HIST("QA/trkppionTOFPID"), track.pt(), track.tofNSigmaPi());
        histos.fill(HIST("QA/trkppionTPCTOFPID"), track.tpcNSigmaPi(), track.tofNSigmaPi());
      }
      histos.fill(HIST("QA/trkppionpT"), track.pt());
      histos.fill(HIST("QA/trkppionDCAxy"), track.dcaXY());
      histos.fill(HIST("QA/trkppionDCAz"), track.dcaZ());
    } else if constexpr (Folder == QAFolder::After) {
      histos.fill(HIST("QAcut/trkppionTPCPID"), track.pt(), track.tpcNSigmaPi());
      if (hasTOF) {
        histos.fill(HIST("QAcut/trkppionTOFPID"), track.pt(), track.tofNSigmaPi());
        histos.fill(HIST("QAcut/trkppionTPCTOFPID"), track.tpcNSigmaPi(), track.tofNSigmaPi());
      }
      histos.fill(HIST("QAcut/trkppionpT"), track.pt());
      histos.fill(HIST("QAcut/trkppionDCAxy"), track.dcaXY());
      histos.fill(HIST("QAcut/trkppionDCAz"), track.dcaZ());
    } else {
      histos.fill(HIST("QAMC/trkppionTPCPID"), track.pt(), track.tpcNSigmaPi());
      if (hasTOF) {
        histos.fill(HIST("QAMC/trkppionTOFPID"), track.pt(), track.tofNSigmaPi());
        histos.fill(HIST("QAMC/trkppionTPCTOFPID"), track.tpcNSigmaPi(), track.tofNSigmaPi());
      }
      histos.fill(HIST("QAMC/trkppionpT"), track.pt());
      histos.fill(HIST("QAMC/trkppionDCAxy"), track.dcaXY());
      histos.fill(HIST("QAMC/trkppionDCAz"), track.dcaZ());
    }
  } else {
    if constexpr (Folder == QAFolder::Before) {
      histos.fill(HIST("QA/trkspionTPCPID"), track.pt(), track.tpcNSigmaPi());
      if (hasTOF) {
        histos.fill(HIST("QA/trkspionTOFPID"), track.pt(), track.tofNSigmaPi());
        histos.fill(HIST("QA/trkspionTPCTOFPID"), track.tpcNSigmaPi(), track.tofNSigmaPi());
      }
      histos.fill(HIST("QA/trkspionpT"), track.pt());
      histos.fill(HIST("QA/trkspionDCAxy"), track.dcaXY());
      histos.fill(HIST("QA/trkspionDCAz"), track.dcaZ());
    } else if constexpr (Folder == QAFolder::After) {
      histos.fill(HIST("QAcut/trkspionTPCPID"), track.pt(), track.tpcNSigmaPi());
      if (hasTOF) {
        histos.fill(HIST("QAcut/trkspionTOFPID"), track.pt(), track.tofNSigmaPi());
        histos.fill(HIST("QAcut/trkspionTPCTOFPID"), track.tpcNSigmaPi(), track.tofNSigmaPi());
      }
      histos.fill(HIST("QAcut/trkspionpT"), track.pt());
      histos.fill(HIST("QAcut/trkspionDCAxy"), track.dcaXY());
      histos.fill(HIST("QAcut/trkspionDCAz"), track.dcaZ());
    } else {
      histos.fill(HIST("QAMC/trkspionTPCPID"), track.pt(), track.tpcNSigmaPi());
      if (hasTOF) {
        histos.fill(HIST("QAMC/trkspionTOFPID"), track.pt(), track.tofNSigmaPi());
        histos.fill(HIST("QAMC/trkspionTPCTOFPID"), track.tpcNSigmaPi(), track.tofNSigmaPi());
      }
      histos.fill(HIST("QAMC/trkspionpT"), track.pt());
      histos.fill(HIST("QAMC/trkspionDCAxy"), track.dcaXY());
      histos.fill(HIST("QAMC/trkspionDCAz"), track.dcaZ());
    }
  }
}

// Track QA of the bachelor kaon
template <QAFolder Folder, typename TrackType>
void fillKaonQA(o2::framework::HistogramRegistry& histos, const TrackType& track)
{
  const bool hasTOF = track.hasTOF();
  if constexpr (Folder == QAFolder::Before) {
    histos.fill(HIST("QA/trkkaonTPCPID"), track.pt(), track.tpcNSigmaKa());
    if (hasTOF) {
      histos.fill(HIST("QA/trkkaonTOFPID"), track.pt(), track.tofNSigmaKa());
      histos.fill(HIST("QA/trkkaonTPCTOFPID"), track.tpcNSigmaKa(), track.tofNSigmaKa());
    }
    histos.fill(HIST("QA/trkkaonpT"), track.pt());
    histos.fill(HIST("QA/trkkaonDCAxy"), track.dcaXY());
    histos.fill(HIST("QA/trkkaonDCAz"), track.dcaZ());
  } else if constexpr (Folder == QAFolder::After) {
    histos.fill(HIST("QAcut/trkkaonTPCPID"), track.pt(), track.tpcNSigmaKa());
    if (hasTOF) {
      histos.fill(HIST("QAcut/trkkaonTOFPID"), track.pt(), track.tofNSigmaKa());
      histos.fill(HIST("QAcut/trkkaonTPCTOFPID"), track.tpcNSigmaKa(), track.tofNSigmaKa());
    }
    histos.fill(HIST("QAcut/trkkaonpT"), track.pt());
    histos.fill(HIST("QAcut/trkkaonDCAxy"), track.dcaXY());
    histos.fill(HIST("QAcut/trkkaonDCAz"), track.dcaZ());
  } else {
    histos.fill(HIST("QAMC/trkkaonTPCPID"), track.pt(), track.tpcNSigmaKa());
    if (hasTOF) {
      histos.fill(HIST("QAMC/trkkaonTOFPID"), track.pt(), track.tofNSigmaKa());
      histos.fill(HIST("QAMC/trkkaonTPCTOFPID"), track.tpcNSigmaKa(), track.tofNSigmaKa());
    }
    histos.fill(HIST("QAMC/trkkaonpT"), track.pt());
    histos.fill(HIST("QAMC/trkkaonDCAxy"), track.dcaXY());
    histos.fill(HIST("QAMC/trkkaonDCAz"), track.dcaZ());
  }
}

/// Selection and candidate loop shared by the K1 tasks.
/// The task owns the configurable groups and the histogram registry and passes them in init().
class K1AnalysisMicroCore
{
 public:
  void init(o2::framework::HistogramRegistry& histos,
            EventCuts const& eventCuts, TrackCuts const& trackCuts,
            PionPidCuts const& pionPidCuts, KaonPidCuts const& kaonPidCuts,
            SecondaryCuts const& secondaryCuts, CandidateCuts const& candidateCuts,
            HistogramOptions const& histogramOptions, ProcessModes const& modes,
            LooseStageOptions const& looseOptions = {})
  {
    mEventCuts = eventCuts;
    mTrackCuts = trackCuts;
    mPionPid = pionPidCuts;
    mKaonPid = kaonPidCuts;
    mSecondaryCuts = secondaryCuts;
    mCandidateCuts = candidateCuts;
    mHistogramOptions = histogramOptions;
    mLooseOptions = looseOptions;
    mTruthDebugCounts = {};

    mSecondaryWindowOn = isCutEnabled(mSecondaryCuts.cSecondaryMasswindow);
    mAnotherMassCutOn = isCutEnabled(mSecondaryCuts.cMinAnotherSecondaryMassCut) || isCutEnabled(mSecondaryCuts.cMaxAnotherSecondaryMassCut);
    mPiKaMassCutOn = isCutEnabled(mSecondaryCuts.cMinPiKaMassCut) || isCutEnabled(mSecondaryCuts.cMaxPiKaMassCut);
    mAngleCutOn = isCutEnabled(mSecondaryCuts.cMinAngle) || isCutEnabled(mSecondaryCuts.cMaxAngle);
    mPairAsymCutOn = isCutEnabled(mSecondaryCuts.cMinPairAsym) || isCutEnabled(mSecondaryCuts.cMaxPairAsym);

    checkConfiguration(modes);
    registerHistograms(histos, modes);
  }

  template <typename CollisionType>
  bool passesEventCuts(const CollisionType& collision)
  {
    return !(mEventCuts.cRecoINELgt0 && !collision.isRecINELgt0());
  }

  template <typename CollisionType>
  bool passesMCEventCuts(const CollisionType& collision)
  {
    if (mEventCuts.cMCINELgt0 && !collision.isINELgt0()) {
      return false;
    }
    if (mEventCuts.cMCVtxIn10 && !collision.isVtxIn10()) {
      return false;
    }
    return true;
  }

  // Resolve the PID cut of one species at a given pT; false if the pT is outside all pT-dependent bins.
  template <Species S>
  bool getPIDCut(float pt, PIDCut& cut)
  {
    if constexpr (S == Species::Pion) {
      cut.tpcMax = mPionPid.cMaxTPCnSigmaPion;
      cut.tofMax = mPionPid.cMaxTOFnSigmaPion;
      cut.combined = mPionPid.nsigmaCutCombinedPion;
      cut.tofRequired = false;
      if (mPionPid.cPionUsePtDepPID) {
        const int ptBin = getPtBinIndex(pt, mPionPid.cPionPIDPtBins.value);
        if (ptBin < 0) {
          return false;
        }
        const auto bin = static_cast<std::size_t>(ptBin);
        cut.tpcMax = mPionPid.cPionTPCNSigmaCuts.value[bin];
        cut.tofMax = mPionPid.cPionTOFNSigmaCuts.value[bin];
        cut.tofRequired = mPionPid.cPionTOFRequired.value[bin] != 0;
      }
    } else {
      cut.tpcMax = mKaonPid.cMaxTPCnSigmaKaon;
      cut.tofMax = mKaonPid.cMaxTOFnSigmaKaon;
      cut.combined = mKaonPid.nsigmaCutCombinedKaon;
      cut.tofRequired = false;
      if (mKaonPid.cKaonUsePtDepPID) {
        const int ptBin = getPtBinIndex(pt, mKaonPid.cKaonPIDPtBins.value);
        if (ptBin < 0) {
          return false;
        }
        const auto bin = static_cast<std::size_t>(ptBin);
        cut.tpcMax = mKaonPid.cKaonTPCNSigmaCuts.value[bin];
        cut.tofMax = mKaonPid.cKaonTOFNSigmaCuts.value[bin];
        cut.tofRequired = mKaonPid.cKaonTOFRequired.value[bin] != 0;
      }
    }
    return true;
  }

  // Track quality selection shared by pion and kaon. Returns the last stage that was passed.
  // Full tracks store exact values; micro tracks store quantised DCA (see LFResonanceTables.h).
  template <bool IsResoMicrotrack, typename TrackType>
  int trackQualityStage(const TrackType& track)
  {
    const double pt = track.pt();
    const double dcaXY = track.dcaXY();
    const double dcaZ = track.dcaZ();
    // Invalid micro DCA codes decode to NaN
    if (!std::isfinite(pt) || !std::isfinite(track.eta()) || !std::isfinite(dcaXY) || !std::isfinite(dcaZ)) {
      return kTrkInput;
    }
    if (std::abs(pt) < mTrackCuts.cMinPtcut) {
      return kTrkInput;
    }
    if (isCutEnabled(mTrackCuts.cMaxEtacut) && !(std::abs(track.eta()) < mTrackCuts.cMaxEtacut)) {
      return kTrkPt;
    }

    if (mTrackCuts.cfgUsePtDepDCA) {
      if constexpr (IsResoMicrotrack) {
        if (!track.passedPtDependentDCAxy()) {
          return kTrkEta;
        }
        if (!track.passedPtDependentDCAz()) {
          return kTrkDCAxy;
        }
      } else {
        const double dcaPtCut = mTrackCuts.cDCAToPVByPtP0 + mTrackCuts.cDCAToPVByPtCoeff * std::pow(pt, -static_cast<double>(mTrackCuts.cDCAToPVByPtPower));
        if (!(std::abs(dcaXY) < dcaPtCut)) {
          return kTrkEta;
        }
        if (!(std::abs(dcaZ) < dcaPtCut)) {
          return kTrkDCAxy;
        }
      }
    } else {
      if (isCutEnabled(mTrackCuts.cMaxDCArToPVcut)) {
        if constexpr (IsResoMicrotrack) {
          if (!passesBinnedMax(dcaXY, mTrackCuts.cMaxDCArToPVcut)) {
            return kTrkEta;
          }
        } else {
          if (!(std::abs(dcaXY) <= mTrackCuts.cMaxDCArToPVcut)) {
            return kTrkEta;
          }
        }
      }
      if (isCutEnabled(mTrackCuts.cMaxDCAzToPVcut)) {
        if constexpr (IsResoMicrotrack) {
          if (!passesBinnedMax(dcaZ, mTrackCuts.cMaxDCAzToPVcut)) {
            return kTrkDCAxy;
          }
        } else {
          if (!(std::abs(dcaZ) <= mTrackCuts.cMaxDCAzToPVcut)) {
            return kTrkDCAxy;
          }
        }
      }
    }
    if (isCutEnabled(mTrackCuts.cMinDCAzToPVcut)) {
      if constexpr (IsResoMicrotrack) {
        if (!passesBinnedMin(dcaZ, mTrackCuts.cMinDCAzToPVcut)) {
          return kTrkDCAxy;
        }
      } else {
        if (!(std::abs(dcaZ) >= mTrackCuts.cMinDCAzToPVcut)) {
          return kTrkDCAxy;
        }
      }
    }

    // Track flags
    if ((mTrackCuts.cfgPrimaryTrack && !track.isPrimaryTrack()) ||
        (mTrackCuts.cfgGlobalWoDCATrack && !track.isGlobalTrackWoDCA()) ||
        (mTrackCuts.cfgGlobalTrack && !track.isGlobalTrack()) ||
        (mTrackCuts.cfgPVContributor && !track.isPVContributor()) ||
        (mTrackCuts.cfgUseITSRefit && !track.passedITSRefit()) ||
        (mTrackCuts.cfgUseTPCRefit && !track.passedTPCRefit())) {
      return kTrkDCAz;
    }

    // Clusters: found clusters exist only in ResoTracks, ITS clusters only in ResoMicroTracks
    if constexpr (!IsResoMicrotrack) {
      if constexpr (requires { track.tpcNClsFound(); }) {
        if (track.tpcNClsFound() < mTrackCuts.cfgTPCcluster) {
          return kTrkFlags;
        }
      }
    }
    if constexpr (requires { track.tpcNClsCrossedRows(); }) {
      if (track.tpcNClsCrossedRows() < mTrackCuts.cfgTPCCrossedRowsMin) {
        return kTrkFlags;
      }
    }
    if constexpr (IsResoMicrotrack) {
      if constexpr (requires { track.itsNCls(); }) {
        if (track.itsNCls() < mTrackCuts.cfgITSNClsMin) {
          return kTrkFlags;
        }
      }
    }
    return kTrkClusters;
  }

  // TOF signal requirement of the track (global, per species, or per pT bin)
  template <Species S, typename TrackType>
  bool passesTOFRequired(const TrackType& track)
  {
    bool required = mTrackCuts.cfgHasTOF;
    if constexpr (S == Species::Pion) {
      required = required || mPionPid.cUseOnlyTOFTrackPi;
    } else {
      required = required || mKaonPid.cUseOnlyTOFTrackKa;
    }
    PIDCut cut;
    // A pT outside all bins is rejected by passesPID
    if (!mCandidateCuts.cByPassTOF && getPIDCut<S>(track.pt(), cut) && cut.tofRequired) {
      required = true;
    }
    return !required || track.hasTOF();
  }

  // PID selection: the same code for full and micro tracks, only the comparison is quantisation aware
  template <bool IsResoMicrotrack, Species S, typename TrackType>
  bool passesPID(const TrackType& track)
  {
    PIDCut cut;
    if (!getPIDCut<S>(track.pt(), cut)) {
      return false;
    }
    const bool hasTOF = track.hasTOF();
    const double tpcNSigma = (S == Species::Pion) ? track.tpcNSigmaPi() : track.tpcNSigmaKa();
    double tofNSigma = std::numeric_limits<double>::quiet_NaN(); // TOF value is only valid with hasTOF
    if constexpr (S == Species::Pion) {
      if (hasTOF) {
        tofNSigma = track.tofNSigmaPi();
      }
    } else {
      if (hasTOF) {
        tofNSigma = track.tofNSigmaKa();
      }
    }
    if (isCutEnabled(cut.tpcMax) && !passesMax<IsResoMicrotrack>(std::abs(tpcNSigma), cut.tpcMax)) {
      return false;
    }
    // Missing TOF is handled by passesTOFRequired; here the TPC alone decides
    if (mCandidateCuts.cByPassTOF || !hasTOF) {
      return true;
    }
    bool tofPassed = !isCutEnabled(cut.tofMax) || passesMax<IsResoMicrotrack>(std::abs(tofNSigma), cut.tofMax);
    if (!tofPassed && cut.combined > 0 && tpcNSigma * tpcNSigma + tofNSigma * tofNSigma < cut.combined * cut.combined) {
      tofPassed = true;
    }
    return tofPassed;
  }

  // Full selection stage of a track (quality, TOF requirement, PID)
  template <bool IsResoMicrotrack, Species S, typename TrackType>
  int trackSelectionStage(const TrackType& track)
  {
    const int qualityStage = trackQualityStage<IsResoMicrotrack>(track);
    if (qualityStage < kTrkClusters) {
      return qualityStage;
    }
    if (!passesTOFRequired<S>(track)) {
      return kTrkClusters;
    }
    if (!passesPID<IsResoMicrotrack, S>(track)) {
      return kTrkTOFRequired;
    }
    return kTrkPID;
  }

  // Unordered (pion, pion, kaon) candidate loop of one collision (or one mixed pair of collisions).
  // dTracks1: bachelor kaons, dTracks2: pions. The optional callback receives the unlike-sign micro
  // same-event candidates in the canonical roles (collision, kaon, same-sign pion, opposite-sign pion,
  // truth channel, pass bits) at the loose or selected stage configured by LooseStageOptions.
  template <bool IsMC, bool IsMix, bool IsResoMicrotrack, typename CollisionType, typename TracksType, typename Callback = std::nullptr_t>
  void fillHistograms(o2::framework::HistogramRegistry& histos, const CollisionType& collision, const TracksType& dTracks1, const TracksType& dTracks2, Callback callback = nullptr)
  {
    if (dTracks1.size() == 0 || dTracks2.size() == 0) {
      return;
    }
    constexpr bool HasCallback = !std::is_same_v<Callback, std::nullptr_t>;
    // Sets are local to this reconstructed collision: IDs cannot leak across DFs.
    // Source-file/DF deduplication across split collisions belongs in the audit.
    std::array<std::set<int>, NTruthChannels> matchedMothers;

    // Selection cache: every track is selected once, not once per pair x bachelor.
    // dTracks1: bachelor kaons, dTracks2: pions (different collisions in mixed events).
    constexpr bool FillCutFlow = IsResoMicrotrack && !IsMix;
    const int64_t firstKaonIndex = dTracks1.begin().index();
    const int64_t firstPionIndex = dTracks2.begin().index();
    const auto kaonSelected = buildSelectionCache<IsResoMicrotrack, Species::Kaon, FillCutFlow>(histos, dTracks1, firstKaonIndex);
    const auto pionSelected = buildSelectionCache<IsResoMicrotrack, Species::Pion, FillCutFlow>(histos, dTracks2, firstPionIndex);

    // Only micro same-event ML work needs traversal before conventional cuts.
    // The canonical-candidate validity is also required before a selected-stage callback.
    bool visitLoose = false;
    bool checkValidity = false;
    if constexpr (IsResoMicrotrack && !IsMix) {
      visitLoose = mLooseOptions.audit || (!mLooseOptions.exportSelected && HasCallback);
      checkValidity = visitLoose || (mLooseOptions.exportSelected && HasCallback);
    }
    std::vector<uint8_t> kaonQuality(dTracks1.size(), 0), pionQuality(dTracks2.size(), 0);
    if (visitLoose) {
      for (const auto& track : dTracks1) {
        kaonQuality[getCacheIndex(track, firstKaonIndex, kaonQuality.size())] = trackQualityStage<IsResoMicrotrack>(track) == kTrkClusters;
      }
      for (const auto& track : dTracks2) {
        pionQuality[getCacheIndex(track, firstPionIndex, pionQuality.size())] = trackQualityStage<IsResoMicrotrack>(track) == kTrkClusters;
      }
    }

    // Values needed only by switched-on cuts or QA are computed only then
    const bool isK892Mode = mSecondaryCuts.cfgModeK892orRho;
    const bool fillQA = !IsMix && mHistogramOptions.additionalQAplots;
    const bool needAngle = IsMC || fillQA || mAngleCutOn;
    const bool needPairAsym = IsMC || fillQA || mPairAsymCutOn;
    // K892 mode: the K* candidate is (trk1, K), rho mode: the rho is (trk1, trk2)
    const bool needMass13 = IsMC || fillQA || (isK892Mode ? mSecondaryWindowOn : mAnotherMassCutOn) || (isK892Mode && (needAngle || needPairAsym));
    const bool needMass23 = IsMC || fillQA || mPiKaMassCutOn;
    const bool rhoWindowOn = mSecondaryWindowOn && !isK892Mode;

    auto multiplicity = collision.cent();
    ROOT::Math::PxPyPzMVector lDecayDaughter1, lDecayDaughter2, lResonanceSecondary, lDecayDaughter_bach, lResonanceK1, lPair13, lPair23;
    // Unordered pion pairs: each (pion, pion, kaon) triplet is filled once.
    // Here trk1 is the pion with the lower index; the roles are assigned once the bachelor is known.
    for (const auto& [trk1, trk2] : o2::soa::combinations(o2::soa::CombinationsStrictlyUpperIndexPolicy(dTracks2, dTracks2))) {
      // trk1: pion, trk2: pion, bTrack: kaon
      const bool pionsSelected = pionSelected[getCacheIndex(trk1, firstPionIndex, pionSelected.size())] && pionSelected[getCacheIndex(trk2, firstPionIndex, pionSelected.size())];
      bool pairPt = false;
      bool rhoWindow = true;
      if (pionsSelected || visitLoose) {
        // Resonance reconstruction
        lDecayDaughter1.SetCoordinates(trk1.px(), trk1.py(), trk1.pz(), o2::constants::physics::MassPionCharged);
        lDecayDaughter2.SetCoordinates(trk2.px(), trk2.py(), trk2.pz(), o2::constants::physics::MassPionCharged);
        lResonanceSecondary = lDecayDaughter1 + lDecayDaughter2;
        pairPt = !(lResonanceSecondary.Pt() < mSecondaryCuts.cMinSecondaryPtCut);
        rhoWindow = !rhoWindowOn || isInWindow(lResonanceSecondary.M(), MassRho770, mSecondaryCuts.cSecondaryMasswindow);
      }
      if constexpr (FillCutFlow) {
        // Early stages count potential triplets: each pair carries N bachelor trials.
        // This preserves the pair-first reconstruction and avoids a new cubic data loop.
        // Distinct pion IDs are guaranteed by the strictly upper index policy (stage 1 is always passed).
        const int lastStage = !pionsSelected ? 1 : !pairPt    ? 3
                                                 : !rhoWindow ? 4
                                                              : 5;
        for (int stage = 0; stage <= lastStage; ++stage) {
          histos.fill(HIST("CutFlow/candidates"), stage, 0, static_cast<double>(dTracks1.size()));
        }
        if constexpr (IsMC) {
          // Match before rejecting quality/pT so both channels have an upstream numerator.
          if (std::abs(trk1.pdgCode()) == kPiPlus && trk1.pdgCode() == -trk2.pdgCode()) {
            for (const auto& bachelor : dTracks1) {
              const auto channel = classifyK1Truth(trk1, trk2, bachelor);
              if (channel != K1TruthChannel::None) {
                for (int stage = 0; stage <= lastStage; ++stage) {
                  histos.fill(HIST("CutFlow/candidates"), stage, static_cast<int>(channel));
                }
              }
            }
          }
        }
      }
      if (!pionsSelected && !visitLoose) {
        continue;
      }

      if (fillQA && pionsSelected) {
        fillPionQA<QAFolder::Before>(histos, trk1, true);
        fillPionQA<QAFolder::Before>(histos, trk2, false);
      }

      if (!pairPt && !visitLoose) {
        continue;
      }

      if (fillQA && pionsSelected && pairPt) {
        histos.fill(HIST("QA/hInvmassSecon"), lResonanceSecondary.M());
      }
      if constexpr (IsMC) {
        if (pionsSelected && pairPt) {
          histos.fill(HIST("QAMC/hpT_Secondary"), lResonanceSecondary.Pt());
        }
      }
      // Secondary mass window (rho mode): the bachelor loop is skipped for rejected pairs
      if (!rhoWindow && !visitLoose) {
        continue;
      }

      for (const auto& bTrack : dTracks1) {
        if (bTrack.index() == trk1.index() || bTrack.index() == trk2.index()) {
          continue;
        }
        K1TruthChannel flowChannel = K1TruthChannel::None;
        if constexpr (IsMC && IsResoMicrotrack && !IsMix) {
          flowChannel = classifyK1Truth(trk1, trk2, bTrack);
        }
        auto countCandidate = [&](int stage) {
          if constexpr (FillCutFlow) {
            histos.fill(HIST("CutFlow/candidates"), stage, 0);
            if constexpr (IsMC) {
              if (flowChannel != K1TruthChannel::None) {
                histos.fill(HIST("CutFlow/candidates"), stage, static_cast<int>(flowChannel));
              }
            }
          }
        };
        const bool pairSelected = pionsSelected && pairPt && rhoWindow;
        const bool bachelorSelected = kaonSelected[getCacheIndex(bTrack, firstKaonIndex, kaonSelected.size())];
        if (pairSelected) {
          countCandidate(6);
        }
        if ((!pairSelected || !bachelorSelected) && !visitLoose) {
          continue;
        }
        if (pairSelected && bachelorSelected) {
          countCandidate(7);
        }

        if (fillQA && pairSelected && bachelorSelected) {
          fillKaonQA<QAFolder::Before>(histos, bTrack);
        }

        // Canonical assignment of the pion roles, once the bachelor is known.
        // Unlike-sign pair: the pion with the sign opposite to the kaon is pion 1 (K*0 partner), the other is pion 2.
        // Like-sign pair (the rule is ambiguous): the pion with the lower index is pion 1.
        const bool isUnlikeSign = trk1.sign() * trk2.sign() < 0;
        const bool swapPions = isUnlikeSign && trk1.sign() == bTrack.sign();
        const auto& pion1 = swapPions ? trk2 : trk1;
        const auto& pion2 = swapPions ? trk1 : trk2;
        const auto& lPion1 = swapPions ? lDecayDaughter2 : lDecayDaughter1;
        const auto& lPion2 = swapPions ? lDecayDaughter1 : lDecayDaughter2;

        // K1 reconstruction
        lDecayDaughter_bach.SetCoordinates(bTrack.px(), bTrack.py(), bTrack.pz(), o2::constants::physics::MassKaonCharged);
        lResonanceK1 = lResonanceSecondary + lDecayDaughter_bach;

        auto countMl = [&](int stage) {
          if (mLooseOptions.audit && isUnlikeSign) {
            histos.fill(HIST("ML/looseCutflow"), stage, 0);
            if (flowChannel != K1TruthChannel::None) {
              const int stratum = 2 * static_cast<int>(flowChannel) - (bTrack.sign() > 0 ? 1 : 0);
              histos.fill(HIST("ML/looseCutflow"), stage, stratum);
            }
          }
        };
        bool validLoose = false;
        if constexpr (IsResoMicrotrack && !IsMix) {
          if (checkValidity && isUnlikeSign) {
            const auto canonical = o2::analysis::k1ml::canonicalizeUS(o2::analysis::k1ml::makeTrackSnapshot(bTrack), o2::analysis::k1ml::makeTrackSnapshot(pion2), o2::analysis::k1ml::makeTrackSnapshot(pion1));
            validLoose = canonical.status == o2::analysis::k1ml::BuildStatus::Ok &&
                         std::isfinite(lResonanceK1.M()) && std::isfinite(lResonanceK1.Pt()) &&
                         std::isfinite(lResonanceK1.Rapidity()) && std::isfinite(lResonanceK1.Eta()) &&
                         std::isfinite(lResonanceK1.Phi()) &&
                         o2::analysis::k1ml::buildMasterFeatures(canonical.candidate).status == o2::analysis::k1ml::BuildStatus::Ok;
            if (validLoose) {
              countMl(0);
            }
          }
        }

        // Stage L common acceptance uses the existing inclusive rapidity window.
        if (lResonanceK1.Rapidity() > mCandidateCuts.cK1MaxRap || lResonanceK1.Rapidity() < mCandidateCuts.cK1MinRap) {
          continue;
        }
        if (pairSelected && bachelorSelected) {
          countCandidate(8);
        }

        double mass13 = 0.;
        double mass23 = 0.;
        double lK1Angle = 0.;
        double lPairAsym = 0.;
        if (needMass13) {
          lPair13 = lPion1 + lDecayDaughter_bach;
          mass13 = lPair13.M();
        }
        if (needMass23) {
          lPair23 = lPion2 + lDecayDaughter_bach;
          mass23 = lPair23.M();
        }
        // Rho mode: secondary = (trk1, trk2) against the bachelor. K892 mode: secondary = (trk1, K) against trk2.
        if (needAngle) {
          lK1Angle = isK892Mode ? ROOT::Math::VectorUtil::Angle(lPair13, lPion2) : ROOT::Math::VectorUtil::Angle(lResonanceSecondary, lDecayDaughter_bach);
        }
        if (needPairAsym) {
          lPairAsym = isK892Mode ? (lPair13.E() - lPion2.E()) / (lPair13.E() + lPion2.E())
                                 : (lResonanceSecondary.E() - lDecayDaughter_bach.E()) / (lResonanceSecondary.E() + lDecayDaughter_bach.E());
        }

        // Candidate cuts (each one is evaluated only if switched on)
        const bool candidateCutsPass =
          !(isK892Mode && mSecondaryWindowOn && (!isInWindow(mass13, o2::constants::physics::MassK0Star892, mSecondaryCuts.cSecondaryMasswindow) || pion1.sign() == bTrack.sign())) &&
          !(mAnotherMassCutOn && !isInRange(isK892Mode ? lResonanceSecondary.M() : mass13, mSecondaryCuts.cMinAnotherSecondaryMassCut, mSecondaryCuts.cMaxAnotherSecondaryMassCut)) &&
          !(mPiKaMassCutOn && !isInRange(mass23, mSecondaryCuts.cMinPiKaMassCut, mSecondaryCuts.cMaxPiKaMassCut)) &&
          !(mAngleCutOn && !isInRange(lK1Angle, mSecondaryCuts.cMinAngle, mSecondaryCuts.cMaxAngle)) &&
          !(mPairAsymCutOn && !isInRange(lPairAsym, mSecondaryCuts.cMinPairAsym, mSecondaryCuts.cMaxPairAsym));
        auto emitCandidate = [&](uint16_t passBits) {
          if constexpr (IsResoMicrotrack && !IsMix && HasCallback) {
            callback(collision, bTrack, pion2, pion1, flowChannel, passBits);
          } else {
            static_cast<void>(passBits);
          }
        };
        if (visitLoose && validLoose) {
          countMl(1);
          if (mLooseOptions.audit) {
            histos.fill(HIST("ML/looseMassPtActivity"), lResonanceK1.M(), lResonanceK1.Pt(), collision.cent());
          }
          const bool qualityPass = kaonQuality[getCacheIndex(bTrack, firstKaonIndex, kaonQuality.size())] &&
                                   pionQuality[getCacheIndex(trk1, firstPionIndex, pionQuality.size())] &&
                                   pionQuality[getCacheIndex(trk2, firstPionIndex, pionQuality.size())];
          uint16_t passBits = kPassLoose;
          if (qualityPass) {
            passBits |= kPassQuality;
            countMl(2);
            if (pionsSelected && bachelorSelected) {
              passBits |= kPassPID;
              countMl(3);
              if (pairPt && rhoWindow) {
                passBits |= kPassPair;
                countMl(4);
                if (candidateCutsPass) {
                  passBits |= kPassCandidate;
                  countMl(5);
                }
              }
            }
          }
          if (!mLooseOptions.exportSelected) {
            emitCandidate(passBits);
          }
        }
        // Stage C retains the frozen conventional selections and QA population.
        if (!pairSelected || !bachelorSelected) {
          continue;
        }

        // QA histogram before the candidate cuts
        if (fillQA) {
          histos.fill(HIST("QA/K1OA"), lK1Angle);
          histos.fill(HIST("QA/K1PairAsym"), lPairAsym);
          histos.fill(HIST("QA/hInvmassK892_Rho"), mass13, lResonanceSecondary.M());
          histos.fill(HIST("QA/hInvmassSecon_PiKa"), lResonanceSecondary.M(), mass23);
          histos.fill(HIST("QA/hpT_Secondary"), lResonanceSecondary.Pt());
        }

        if (!candidateCutsPass) {
          continue;
        }
        countCandidate(9);

        // QA histograms after the candidate cuts
        if (fillQA) {
          fillPionQA<QAFolder::After>(histos, pion1, true);
          fillPionQA<QAFolder::After>(histos, pion2, false);
          fillKaonQA<QAFolder::After>(histos, bTrack);
          histos.fill(HIST("QAcut/K1OA"), lK1Angle);
          histos.fill(HIST("QAcut/K1PairAsym"), lPairAsym);
          histos.fill(HIST("QAcut/hInvmassK892_Rho"), mass13, lResonanceSecondary.M());
          histos.fill(HIST("QAcut/hInvmassSecon_PiKa"), lResonanceSecondary.M(), mass23);
          histos.fill(HIST("QAcut/hInvmassSecon"), lResonanceSecondary.M());
          histos.fill(HIST("QAcut/hpT_Secondary"), lResonanceSecondary.Pt());
        }

        countCandidate(isUnlikeSign ? 10 : 11);
        if (isUnlikeSign && mLooseOptions.exportSelected && validLoose) {
          emitCandidate(PassBitsSelected);
        }
        if constexpr (IsMC && IsResoMicrotrack && !IsMix) {
          if (flowChannel != K1TruthChannel::None) {
            const int mother = flowChannel == K1TruthChannel::RhoK ? bTrack.motherId() : std::abs(pion1.motherPDG()) == o2::constants::physics::Pdg::kK1_1270Plus ? pion1.motherId()
                                                                                                                                                                  : pion2.motherId();
            if (matchedMothers[static_cast<int>(flowChannel)].insert(mother).second) {
              histos.fill(HIST("CutFlow/uniqueMothersPerCollision"), static_cast<int>(flowChannel));
            }
          }
        }

        if constexpr (!IsMix) {
          unsigned int typeK1 = bTrack.sign() > 0 ? BinType::kK1P : BinType::kK1N;
          unsigned int typeNormal = BinAnti::kNormal;
          if (isUnlikeSign) {
            histos.fill(HIST("k1invmass"), lResonanceK1.M());
            histos.fill(HIST("hInvmass_K1"), typeNormal, typeK1, multiplicity, lResonanceK1.Pt(), lResonanceK1.M());
          } else {
            histos.fill(HIST("k1invmass_LS"), lResonanceK1.M());
            histos.fill(HIST("hInvmass_K1_LS"), typeNormal, typeK1, multiplicity, lResonanceK1.Pt(), lResonanceK1.M());
          }

          if constexpr (IsMC) {
            const auto channel = classifyK1Truth(pion1, pion2, bTrack);
            const int channelBin = static_cast<int>(channel);
            histos.fill(HIST("MCReco/channel"), channelBin);
            histos.fill(HIST("MCReco/mass"), channelBin, lResonanceK1.M());
            histos.fill(HIST("MCReco/pt"), channelBin, lResonanceK1.Pt());
            histos.fill(HIST("MCReco/piPiMass"), channelBin, lResonanceSecondary.M());
            histos.fill(HIST("MCReco/pi1KMass"), channelBin, mass13);
            histos.fill(HIST("MCReco/pi2KMass"), channelBin, mass23);
            if (channel != K1TruthChannel::None) {
              if (mTruthDebugCounts[channelBin] < mHistogramOptions.cfgTruthDebug) {
                ++mTruthDebugCounts[channelBin];
                LOGP(info, "K1Truth channel={} collision={} tracks=({},{},{}) pdg=({},{},{}) mothers=({},{},{}) motherPDG=({},{},{}) siblings=(({},{}),({},{}),({},{}))",
                     channelBin, collision.globalIndex(),
                     pion1.globalIndex(), pion2.globalIndex(), bTrack.globalIndex(),
                     pion1.pdgCode(), pion2.pdgCode(), bTrack.pdgCode(), pion1.motherId(), pion2.motherId(), bTrack.motherId(),
                     pion1.motherPDG(), pion2.motherPDG(), bTrack.motherPDG(),
                     pion1.siblingIds()[0], pion1.siblingIds()[1], pion2.siblingIds()[0], pion2.siblingIds()[1], bTrack.siblingIds()[0], bTrack.siblingIds()[1]);
              }
              typeK1 = bTrack.sign() > 0 ? BinType::kK1P_Rec : BinType::kK1N_Rec;
              histos.fill(HIST("hInvmass_K1"), typeNormal, typeK1, multiplicity, lResonanceK1.Pt(), lResonanceK1.M());
              histos.fill(HIST("k1invmass_MC"), lResonanceK1.M());
              histos.fill(HIST("QAMC/K1OA"), lK1Angle);
              histos.fill(HIST("QAMC/K1PairAsym"), lPairAsym);
              histos.fill(HIST("QAMC/hInvmassK892_Rho"), mass13, lResonanceSecondary.M());
              histos.fill(HIST("QAMC/hInvmassSecon_PiKa"), lResonanceSecondary.M(), mass23);
              histos.fill(HIST("QAMC/hInvmassSecon"), lResonanceSecondary.M());
              histos.fill(HIST("QAMC/hpT_Secondary"), lResonanceSecondary.Pt());

              // PID QA primary and secondary pion
              fillPionQA<QAFolder::MC>(histos, pion1, true);
              fillPionQA<QAFolder::MC>(histos, pion2, false);
              fillKaonQA<QAFolder::MC>(histos, bTrack);
            } else {
              histos.fill(HIST("k1invmass_MC_noK1"), lResonanceK1.M());
            }
          } // IsMC
        } else {
          unsigned int typeK1 = bTrack.sign() > 0 ? BinType::kK1P_Mix : BinType::kK1N_Mix;
          unsigned int typeNormal = BinAnti::kNormal;
          histos.fill(HIST("hInvmass_K1_Mix"), typeNormal, typeK1, multiplicity, lResonanceK1.Pt(), lResonanceK1.M());
          histos.fill(HIST("k1invmass_Mix"), lResonanceK1.M());
        }
      } // bTrack
    }
  } // fillHistograms

  // Generated K1 parents of a selected reconstructed MC collision. The optional callback receives
  // (parent, immediate channel) for the parents inside the K1 rapidity window.
  // Parents belong to selected reconstructed events; split reco collisions
  // repeat parent sets. This is not an unconditional generated denominator.
  template <typename ParentsType, typename Callback = std::nullptr_t>
  void fillGenerated(o2::framework::HistogramRegistry& histos, const ParentsType& resoParents, Callback callback = nullptr)
  {
    for (const auto& part : resoParents) {
      if (std::abs(part.pdgCode()) != o2::constants::physics::Pdg::kK1_1270Plus) {
        continue;
      }
      const int charge = part.pdgCode() > 0 ? 1 : -1;
      const K1TruthChannel channel = classifyGeneratedK1(charge, part.daughterPDG1(), part.daughterPDG2());
      histos.fill(HIST("CutFlow/generated"), 0, static_cast<int>(channel));
      if (part.y() < mCandidateCuts.cK1MinRap || part.y() > mCandidateCuts.cK1MaxRap) {
        continue;
      }
      histos.fill(HIST("CutFlow/generated"), 1, static_cast<int>(channel));
      // Keep other/unresolved immediate decays too; never require both pairs.
      histos.fill(HIST("MCGen/chargeChannel"), charge, static_cast<int>(channel));
      histos.fill(HIST("MCGen/ptChannel"), static_cast<int>(channel), part.pt());
      if constexpr (!std::is_same_v<Callback, std::nullptr_t>) {
        callback(part, channel);
      }
    }
  }

 private:
  // Selection cache of one track slice. The row number of a grouped slice is global, hence the offset.
  template <typename TrackType>
  static std::size_t getCacheIndex(const TrackType& track, int64_t firstIndex, std::size_t size)
  {
    const int64_t index = static_cast<int64_t>(track.index()) - firstIndex;
    if (index < 0 || index >= static_cast<int64_t>(size)) {
      LOG(fatal) << "Track index " << track.index() << " is outside the selection cache [" << firstIndex << ", " << firstIndex + static_cast<int64_t>(size) << ")";
    }
    return static_cast<std::size_t>(index);
  }

  template <bool IsResoMicrotrack, Species S, bool FillCutFlow, typename TracksType>
  std::vector<uint8_t> buildSelectionCache(o2::framework::HistogramRegistry& histos, const TracksType& tracks, int64_t firstIndex)
  {
    std::vector<uint8_t> selected(tracks.size(), 0);
    for (const auto& track : tracks) {
      const int stage = trackSelectionStage<IsResoMicrotrack, S>(track);
      selected[getCacheIndex(track, firstIndex, selected.size())] = (stage == kTrkPID) ? 1 : 0;
      if constexpr (FillCutFlow) {
        for (int i = 0; i <= stage; ++i) {
          histos.fill(HIST("CutFlow/tracks"), i, static_cast<int>(S));
        }
      }
    }
    return selected;
  }

  void checkConfiguration(ProcessModes const& modes)
  {
    // Consistency of the pT dependent PID configuration
    if (mPionPid.cPionUsePtDepPID) {
      const auto& bins = mPionPid.cPionPIDPtBins.value;
      if (bins.size() < MinPtBinEdges || mPionPid.cPionTPCNSigmaCuts.value.size() != bins.size() - 1 ||
          mPionPid.cPionTOFNSigmaCuts.value.size() != bins.size() - 1 || mPionPid.cPionTOFRequired.value.size() != bins.size() - 1) {
        LOG(fatal) << "Pion pT dependent PID vectors must have (number of pT bin edges - 1) entries";
      }
    }
    if (mKaonPid.cKaonUsePtDepPID) {
      const auto& bins = mKaonPid.cKaonPIDPtBins.value;
      if (bins.size() < MinPtBinEdges || mKaonPid.cKaonTPCNSigmaCuts.value.size() != bins.size() - 1 ||
          mKaonPid.cKaonTOFNSigmaCuts.value.size() != bins.size() - 1 || mKaonPid.cKaonTOFRequired.value.size() != bins.size() - 1) {
        LOG(fatal) << "Kaon pT dependent PID vectors must have (number of pT bin edges - 1) entries";
      }
    }
    if (mCandidateCuts.cByPassTOF && (mPionPid.cUseOnlyTOFTrackPi || mKaonPid.cUseOnlyTOFTrackKa)) {
      LOG(warning) << "cByPassTOF skips the TOF nSigma cut, but cUseOnlyTOFTrack* still requires a TOF signal";
    }

    // Micro tracks store quantised DCA and nSigma: a cut off the grid would silently act as a different cut.
    if (!modes.microTracks) {
      return;
    }
    auto checkDCAGrid = [](const char* name, double cut) {
      const double nearest = std::min(std::max(std::round(cut / DCAGridStep) * DCAGridStep, 0.), DCAGridMax);
      if (std::abs(cut - nearest) > GridTolerance) {
        LOG(fatal) << name << " = " << cut << " is not on the quantised DCA grid (multiples of " << DCAGridStep << " up to " << DCAGridMax << "); nearest value: " << nearest;
      }
    };
    auto checkPIDGrid = [](const char* name, double cut) {
      const double nearest = std::min(std::max(PIDGridStart + std::round((cut - PIDGridStart) / PIDGridStep) * PIDGridStep, PIDGridStart), PIDGridMax);
      if (std::abs(cut - nearest) > GridTolerance) {
        LOG(fatal) << name << " = " << cut << " is not on the quantised nSigma grid ({2.0, 2.25, ..., 3.5}); nearest value: " << nearest;
      }
    };
    if (mTrackCuts.cfgUsePtDepDCA) {
      LOG(info) << "Micro tracks use the producer pT dependent DCA flags (0.004 + 0.013 / pT); cDCAToPVByPt* are ignored";
      if (std::abs(mTrackCuts.cDCAToPVByPtP0 - ProducerDCAPtP0) > ConfigTolerance || std::abs(mTrackCuts.cDCAToPVByPtCoeff - ProducerDCAPtCoeff) > ConfigTolerance || std::abs(mTrackCuts.cDCAToPVByPtPower - ProducerDCAPtPower) > ConfigTolerance) {
        LOG(warning) << "cDCAToPVByPt* differ from the producer defaults, but micro tracks always use the producer formula";
      }
    } else {
      if (isCutEnabled(mTrackCuts.cMaxDCArToPVcut)) {
        checkDCAGrid("cMaxDCArToPVcut", mTrackCuts.cMaxDCArToPVcut);
      }
      if (isCutEnabled(mTrackCuts.cMaxDCAzToPVcut)) {
        checkDCAGrid("cMaxDCAzToPVcut", mTrackCuts.cMaxDCAzToPVcut);
      }
    }
    if (isCutEnabled(mTrackCuts.cMinDCAzToPVcut)) {
      checkDCAGrid("cMinDCAzToPVcut", mTrackCuts.cMinDCAzToPVcut);
    }
    if (isCutEnabled(mPionPid.cMaxTPCnSigmaPion) && !mPionPid.cPionUsePtDepPID) {
      checkPIDGrid("cMaxTPCnSigmaPion", mPionPid.cMaxTPCnSigmaPion);
    }
    if (isCutEnabled(mPionPid.cMaxTOFnSigmaPion) && !mPionPid.cPionUsePtDepPID) {
      checkPIDGrid("cMaxTOFnSigmaPion", mPionPid.cMaxTOFnSigmaPion);
    }
    if (isCutEnabled(mKaonPid.cMaxTPCnSigmaKaon) && !mKaonPid.cKaonUsePtDepPID) {
      checkPIDGrid("cMaxTPCnSigmaKaon", mKaonPid.cMaxTPCnSigmaKaon);
    }
    if (isCutEnabled(mKaonPid.cMaxTOFnSigmaKaon) && !mKaonPid.cKaonUsePtDepPID) {
      checkPIDGrid("cMaxTOFnSigmaKaon", mKaonPid.cMaxTOFnSigmaKaon);
    }
    if (mPionPid.cPionUsePtDepPID) {
      for (const auto& cut : mPionPid.cPionTPCNSigmaCuts.value) {
        if (isCutEnabled(cut)) {
          checkPIDGrid("cPionTPCNSigmaCuts", cut);
        }
      }
      for (const auto& cut : mPionPid.cPionTOFNSigmaCuts.value) {
        if (isCutEnabled(cut)) {
          checkPIDGrid("cPionTOFNSigmaCuts", cut);
        }
      }
    }
    if (mKaonPid.cKaonUsePtDepPID) {
      for (const auto& cut : mKaonPid.cKaonTPCNSigmaCuts.value) {
        if (isCutEnabled(cut)) {
          checkPIDGrid("cKaonTPCNSigmaCuts", cut);
        }
      }
      for (const auto& cut : mKaonPid.cKaonTOFNSigmaCuts.value) {
        if (isCutEnabled(cut)) {
          checkPIDGrid("cKaonTOFNSigmaCuts", cut);
        }
      }
    }
    if (mPionPid.nsigmaCutCombinedPion > 0 || mKaonPid.nsigmaCutCombinedKaon > 0) {
      LOG(warning) << "nsigmaCutCombined* on micro tracks uses quantised nSigma values (approximate)";
    }
  }

  void registerHistograms(o2::framework::HistogramRegistry& histos, ProcessModes const& modes)
  {
    using o2::framework::AxisSpec;
    using o2::framework::HistType;
    const int nBinsDiv = mHistogramOptions.cNbinsDiv;
    std::vector<double> centBinning = {0., 1., 5., 10., 15., 20., 25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 80., 90., 100., 200.};
    AxisSpec centAxis = {centBinning, "T0M (%)"};
    AxisSpec ptAxis = {150, 0, 15, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec dcaxyAxis = {300, 0, 3, "DCA_{#it{xy}} (cm)"};
    AxisSpec dcazAxis = {500, 0, 5, "DCA_{#it{z}} (cm)"};
    AxisSpec invMassAxisK892 = {1400 / nBinsDiv, 0.6, 2.0, "Invariant Mass (GeV/#it{c}^2)"};   // K(892)0
    AxisSpec invMassAxisRho = {2000 / nBinsDiv, 0.0, 2.0, "Invariant Mass (GeV/#it{c}^2)"};    // rho
    AxisSpec invMassAxisReso = {1600 / nBinsDiv, 0.9f, 2.5f, "Invariant Mass (GeV/#it{c}^2)"}; // K1
    AxisSpec pidQAAxis = {130, -6.5, 6.5};

    // THnSparse
    AxisSpec axisAnti = {BinAnti::kNAEnd, 0, BinAnti::kNAEnd, "Type of bin: Normal or Anti"};
    AxisSpec axisType = {BinType::kTYEnd, 0, BinType::kTYEnd, "Type of bin with charge and mix"};

    if (mLooseOptions.audit) {
      auto flow = histos.add<TH2>("ML/looseCutflow", "US triplets;stage;signal stratum", HistType::kTH2D, {{6, -0.5, 5.5}, {5, -0.5, 4.5}});
      const std::array<const char*, 6> labels{"structural US", "loose acceptance", "track quality", "TOF + PID", "pair requirements", "selected US"};
      const std::array<const char*, 5> strata{"all US", "rhoK+", "rhoK-", "KstarPi+", "KstarPi-"};
      for (std::size_t i = 0; i < labels.size(); ++i) {
        flow->GetXaxis()->SetBinLabel(i + 1, labels[i]);
      }
      for (std::size_t i = 0; i < strata.size(); ++i) {
        flow->GetYaxis()->SetBinLabel(i + 1, strata[i]);
      }
      histos.add("ML/looseMassPtActivity", "Loose US;mass (GeV/c^{2});pT (GeV/c);FT0M percentile", HistType::kTH3D,
                 {{300, 0.7, 3.7}, {{0., 0.5, 1., 2., 3., 5., 8., 15., 30., 100.}, "pT"}, {{0., 10., 30., 50., 70., 100., 110.}, "FT0M percentile"}});
    }

    // Micro-only instrumentation: category 0 includes all combinations, not just unmatched.
    auto trackFlow = histos.add<TH2>("CutFlow/tracks", "Micro tracks, once per selected collision;stage;species", HistType::kTH2D, {{static_cast<int>(kTrkNStages), -0.5, static_cast<int>(kTrkNStages) - 0.5}, {2, -0.5, 1.5}});
    const std::array<const char*, TrackStage::kTrkNStages> trackLabels{"input", "pT", "eta", "DCAxy", "DCAz", "track flags", "clusters / crossed rows", "TOF required", "PID"};
    for (std::size_t i = 0; i < trackLabels.size(); ++i) {
      trackFlow->GetXaxis()->SetBinLabel(i + 1, trackLabels[i]);
    }
    trackFlow->GetYaxis()->SetBinLabel(1, "pion");
    trackFlow->GetYaxis()->SetBinLabel(2, "kaon");
    auto candidateFlow = histos.add<TH2>("CutFlow/candidates", "Unordered micro triplets;stage;category", HistType::kTH2D, {{NCandidateStages, -0.5, NCandidateStages - 0.5}, {3, -0.5, 2.5}});
    const std::array<const char*, NCandidateStages> candidateLabels{"input unordered triplets", "distinct pion IDs", "pion selection (quality+PID)", "pion pair constructed", "pair pT", "secondary mass window (rho mode)", "three distinct IDs", "kaon selection (quality+PID)", "K1 rapidity", "candidate cuts", "final US", "final LS"};
    for (std::size_t i = 0; i < candidateLabels.size(); ++i) {
      candidateFlow->GetXaxis()->SetBinLabel(i + 1, candidateLabels[i]);
    }
    candidateFlow->GetYaxis()->SetBinLabel(1, "all");
    candidateFlow->GetYaxis()->SetBinLabel(2, "rho K");
    candidateFlow->GetYaxis()->SetBinLabel(3, "K* pi");
    if (modes.mcRecoMicro) {
      auto mothers = histos.add<TH1>("CutFlow/uniqueMothersPerCollision", "Final unique K1 IDs summed over reconstructed collisions (not globally deduplicated)", HistType::kTH1D, {{2, 0.5, 2.5}});
      mothers->GetXaxis()->SetBinLabel(1, "rho K");
      mothers->GetXaxis()->SetBinLabel(2, "K* pi");
    }
    if (modes.mcGen) {
      auto generated = histos.add<TH2>("CutFlow/generated", "K1 parent rows conditional on selected reconstructed events;stage;immediate channel", HistType::kTH2D, {{2, -0.5, 1.5}, {3, -0.5, 2.5}});
      generated->GetXaxis()->SetBinLabel(1, "all K1 parent rows");
      generated->GetXaxis()->SetBinLabel(2, "K1 rapidity window");
      generated->GetYaxis()->SetBinLabel(1, "other / unresolved");
      generated->GetYaxis()->SetBinLabel(2, "rho K");
      generated->GetYaxis()->SetBinLabel(3, "K* pi");
    }

    // DCA QA
    // Primary pion
    histos.add("QA/trkppionDCAxy", "DCAxy disstribution of primary pion candidates", HistType::kTH1F, {dcaxyAxis});
    histos.add("QA/trkppionDCAz", "DCAz disstribution of primary pion candidates", HistType::kTH1F, {dcazAxis});
    histos.add("QA/trkppionpT", "pT distribution of primary pion candidates", HistType::kTH1F, {ptAxis});
    histos.add("QA/trkppionTPCPID", "TPC PID of primary pion candidates", HistType::kTH2F, {ptAxis, pidQAAxis});
    histos.add("QA/trkppionTOFPID", "TOF PID of primary pion candidates", HistType::kTH2F, {ptAxis, pidQAAxis});
    histos.add("QA/trkppionTPCTOFPID", "TPC-TOF PID map of primary pion candidates", HistType::kTH2F, {pidQAAxis, pidQAAxis});

    histos.add("QAcut/trkppionDCAxy", "DCAxy distribution of primary pion candidates", HistType::kTH1F, {dcaxyAxis});
    histos.add("QAcut/trkppionDCAz", "DCAz distribution of primary pion candidates", HistType::kTH1F, {dcazAxis});
    histos.add("QAcut/trkppionpT", "pT distribution of primary pion candidates", HistType::kTH1F, {ptAxis});
    histos.add("QAcut/trkppionTPCPID", "TPC PID of primary pion candidates", HistType::kTH2F, {ptAxis, pidQAAxis});
    histos.add("QAcut/trkppionTOFPID", "TOF PID of primary pion candidates", HistType::kTH2F, {ptAxis, pidQAAxis});
    histos.add("QAcut/trkppionTPCTOFPID", "TPC-TOF PID map of primary pion candidates", HistType::kTH2F, {pidQAAxis, pidQAAxis});

    // Secondary pion
    histos.add("QA/trkspionDCAxy", "DCAxy distribution of secondary pion candidates", HistType::kTH1F, {dcaxyAxis});
    histos.add("QA/trkspionDCAz", "DCAz distribution of secondary pion candidates", HistType::kTH1F, {dcazAxis});
    histos.add("QA/trkspionpT", "pT distribution of secondary pion candidates", HistType::kTH1F, {ptAxis});
    histos.add("QA/trkspionTPCPID", "TPC PID of secondary pion candidates", HistType::kTH2F, {ptAxis, pidQAAxis});
    histos.add("QA/trkspionTOFPID", "TOF PID of secondary pion candidates", HistType::kTH2F, {ptAxis, pidQAAxis});
    histos.add("QA/trkspionTPCTOFPID", "TPC-TOF PID map of secondary pion candidates", HistType::kTH2F, {pidQAAxis, pidQAAxis});

    histos.add("QAcut/trkspionDCAxy", "DCAxy distribution of secondary pion candidates", HistType::kTH1F, {dcaxyAxis});
    histos.add("QAcut/trkspionDCAz", "DCAz distribution of secondary pion candidates", HistType::kTH1F, {dcazAxis});
    histos.add("QAcut/trkspionpT", "pT distribution of secondary pion candidates", HistType::kTH1F, {ptAxis});
    histos.add("QAcut/trkspionTPCPID", "TPC PID of secondary pion candidates", HistType::kTH2F, {ptAxis, pidQAAxis});
    histos.add("QAcut/trkspionTOFPID", "TOF PID of secondary pion candidates", HistType::kTH2F, {ptAxis, pidQAAxis});
    histos.add("QAcut/trkspionTPCTOFPID", "TPC-TOF PID map of secondary pion candidates", HistType::kTH2F, {pidQAAxis, pidQAAxis});

    // Kaon
    histos.add("QA/trkkaonDCAxy", "DCAxy distribution of kaon candidates", HistType::kTH1F, {dcaxyAxis});
    histos.add("QA/trkkaonDCAz", "DCAz distribution of kaon candidates", HistType::kTH1F, {dcazAxis});
    histos.add("QA/trkkaonpT", "pT distribution of kaon candidates", HistType::kTH1F, {ptAxis});
    histos.add("QA/trkkaonTPCPID", "TPC PID of kaon candidates", HistType::kTH2F, {ptAxis, pidQAAxis});
    histos.add("QA/trkkaonTOFPID", "TOF PID of kaon candidates", HistType::kTH2F, {ptAxis, pidQAAxis});
    histos.add("QA/trkkaonTPCTOFPID", "TPC-TOF PID map of kaon candidates", HistType::kTH2F, {pidQAAxis, pidQAAxis});

    histos.add("QAcut/trkkaonDCAxy", "DCAxy distribution of kaon candidates", HistType::kTH1F, {dcaxyAxis});
    histos.add("QAcut/trkkaonDCAz", "DCAz distribution of kaon candidates", HistType::kTH1F, {dcazAxis});
    histos.add("QAcut/trkkaonpT", "pT distribution of kaon candidates", HistType::kTH1F, {ptAxis});
    histos.add("QAcut/trkkaonTPCPID", "TPC PID of kaon candidates", HistType::kTH2F, {ptAxis, pidQAAxis});
    histos.add("QAcut/trkkaonTOFPID", "TOF PID of kaon candidates", HistType::kTH2F, {ptAxis, pidQAAxis});
    histos.add("QAcut/trkkaonTPCTOFPID", "TPC-TOF PID map of kaon candidates", HistType::kTH2F, {pidQAAxis, pidQAAxis});

    // K1
    histos.add("QA/K1OA", "Opening angle of K1(1270)", HistType::kTH1F, {AxisSpec{100, 0, 3.14, "Opening angle"}});
    histos.add("QA/K1PairAsym", "Pair asymmetry of K1(1270)", HistType::kTH1F, {AxisSpec{100, -1, 1, "Pair asymmetry"}});
    histos.add("QA/hInvmassK892_Rho", "Invariant mass of K(892)0 vs Rho(770)", HistType::kTH2F, {invMassAxisK892, invMassAxisRho});
    histos.add("QA/hInvmassSecon_PiKa", "Invariant mass of secondary resonance vs pion-kaon", HistType::kTH2F, {invMassAxisRho, invMassAxisK892});
    histos.add("QA/hInvmassSecon", "Invariant mass of secondary resonance", HistType::kTH1F, {invMassAxisRho});
    histos.add("QA/hpT_Secondary", "pT distribution of secondary resonance", HistType::kTH1F, {ptAxis});

    histos.add("QAcut/K1OA", "Opening angle of K1(1270)", HistType::kTH1F, {AxisSpec{100, 0, 3.14, "Opening angle"}});
    histos.add("QAcut/K1PairAsym", "Pair asymmetry of K1(1270)", HistType::kTH1F, {AxisSpec{100, -1, 1, "Pair asymmetry"}});
    histos.add("QAcut/hInvmassK892_Rho", "Invariant mass of K(892)0 vs Rho(770)", HistType::kTH2F, {invMassAxisK892, invMassAxisRho});
    histos.add("QAcut/hInvmassSecon_PiKa", "Invariant mass of secondary resonance vs pion-kaon", HistType::kTH2F, {invMassAxisRho, invMassAxisK892});
    histos.add("QAcut/hInvmassSecon", "Invariant mass of secondary resonance", HistType::kTH1F, {invMassAxisRho});
    histos.add("QAcut/hpT_Secondary", "pT distribution of secondary resonance", HistType::kTH1F, {ptAxis});

    // Invariant mass
    histos.add("hInvmass_K1", "Invariant mass of K1(1270) (US)", HistType::kTHnSparseD, {axisAnti, axisType, centAxis, ptAxis, invMassAxisReso});
    histos.add("hInvmass_K1_LS", "Invariant mass of K1(1270) (LS)", HistType::kTHnSparseD, {axisAnti, axisType, centAxis, ptAxis, invMassAxisReso});
    histos.add("hInvmass_K1_Mix", "Invariant mass of K1(1270) (ME)", HistType::kTHnSparseD, {axisAnti, axisType, centAxis, ptAxis, invMassAxisReso});
    // Mass QA (quick check)
    histos.add("k1invmass", "Invariant mass of K1(1270) (US)", HistType::kTH1F, {invMassAxisReso});
    histos.add("k1invmass_LS", "Invariant mass of K1(1270) (LS)", HistType::kTH1F, {invMassAxisReso});
    histos.add("k1invmass_Mix", "Invariant mass of K1(1270) (ME)", HistType::kTH1F, {invMassAxisReso});

    // MC
    if (modes.mcReco) {
      AxisSpec channelAxis = {3, -0.5, 2.5, "0: non-K1, 1: rho K, 2: K* pi"};
      histos.add("MCReco/collisions", "Selected reconstructed MC collisions", HistType::kTH1D, {{1, 0, 1}});
      histos.add("MCReco/microTracks", "Input micro tracks in selected MC collisions", HistType::kTH1D, {{1, 0, 1}});
      histos.add("MCReco/channel", "All selected pi-pi-K combinations by truth channel", HistType::kTH1D, {channelAxis});
      histos.add("MCReco/mass", "Reconstructed mass by truth channel", HistType::kTH2D, {channelAxis, invMassAxisReso});
      histos.add("MCReco/pt", "Reconstructed pT by truth channel", HistType::kTH2D, {channelAxis, ptAxis});
      histos.add("MCReco/piPiMass", "pi-pi mass by truth channel", HistType::kTH2D, {channelAxis, invMassAxisRho});
      histos.add("MCReco/pi1KMass", "First pion-kaon mass by truth channel", HistType::kTH2D, {channelAxis, invMassAxisK892});
      histos.add("MCReco/pi2KMass", "Second pion-kaon mass by truth channel", HistType::kTH2D, {channelAxis, invMassAxisK892});
      histos.add("k1invmass_MC", "Invariant mass of K1(1270)", HistType::kTH1F, {invMassAxisReso});
      histos.add("k1invmass_MC_noK1", "Invariant mass of K1(1270)", HistType::kTH1F, {invMassAxisReso});

      histos.add("QAMC/trkppionDCAxy", "DCAxy distribution of primary pion candidates", HistType::kTH1F, {dcaxyAxis});
      histos.add("QAMC/trkppionDCAz", "DCAz distribution of primary pion candidates", HistType::kTH1F, {dcazAxis});
      histos.add("QAMC/trkppionpT", "pT distribution of primary pion candidates", HistType::kTH1F, {ptAxis});
      histos.add("QAMC/trkppionTPCPID", "TPC PID of primary pion candidates", HistType::kTH2F, {ptAxis, pidQAAxis});
      histos.add("QAMC/trkppionTOFPID", "TOF PID of primary pion candidates", HistType::kTH2F, {ptAxis, pidQAAxis});
      histos.add("QAMC/trkppionTPCTOFPID", "TPC-TOF PID map of primary pion candidates", HistType::kTH2F, {pidQAAxis, pidQAAxis});

      histos.add("QAMC/trkspionDCAxy", "DCAxy distribution of secondary pion candidates", HistType::kTH1F, {dcaxyAxis});
      histos.add("QAMC/trkspionDCAz", "DCAz distribution of secondary pion candidates", HistType::kTH1F, {dcazAxis});
      histos.add("QAMC/trkspionpT", "pT distribution of secondary pion candidates", HistType::kTH1F, {ptAxis});
      histos.add("QAMC/trkspionTPCPID", "TPC PID of secondary pion candidates", HistType::kTH2F, {ptAxis, pidQAAxis});
      histos.add("QAMC/trkspionTOFPID", "TOF PID of secondary pion candidates", HistType::kTH2F, {ptAxis, pidQAAxis});
      histos.add("QAMC/trkspionTPCTOFPID", "TPC-TOF PID map of secondary pion candidates", HistType::kTH2F, {pidQAAxis, pidQAAxis});

      histos.add("QAMC/trkkaonDCAxy", "DCAxy distribution of kaon candidates", HistType::kTH1F, {dcaxyAxis});
      histos.add("QAMC/trkkaonDCAz", "DCAz distribution of kaon candidates", HistType::kTH1F, {dcazAxis});
      histos.add("QAMC/trkkaonpT", "pT distribution of kaon candidates", HistType::kTH1F, {ptAxis});
      histos.add("QAMC/trkkaonTPCPID", "TPC PID of kaon candidates", HistType::kTH2F, {ptAxis, pidQAAxis});
      histos.add("QAMC/trkkaonTOFPID", "TOF PID of kaon candidates", HistType::kTH2F, {ptAxis, pidQAAxis});
      histos.add("QAMC/trkkaonTPCTOFPID", "TPC-TOF PID map of kaon candidates", HistType::kTH2F, {pidQAAxis, pidQAAxis});

      histos.add("QAMC/K1OA", "Opening angle of K1(1270)", HistType::kTH1F, {AxisSpec{100, 0, 3.14, "Opening angle"}});
      histos.add("QAMC/K1PairAsym", "Pair asymmetry of K1(1270)", HistType::kTH1F, {AxisSpec{100, -1, 1, "Pair asymmetry"}});
      histos.add("QAMC/hInvmassK892_Rho", "Invariant mass of K(892)0 vs Rho(770)", HistType::kTH2F, {invMassAxisK892, invMassAxisRho});
      histos.add("QAMC/hInvmassSecon_PiKa", "Invariant mass of secondary resonance vs pion-kaon", HistType::kTH2F, {invMassAxisRho, invMassAxisK892});
      histos.add("QAMC/hInvmassSecon", "Invariant mass of secondary resonance", HistType::kTH1F, {invMassAxisRho});
      histos.add("QAMC/hpT_Secondary", "pT distribution of secondary resonance", HistType::kTH1F, {ptAxis});
    } // mcReco
    if (modes.mcGen) {
      AxisSpec channelAxis = {3, -0.5, 2.5, "0: other/unresolved, 1: rho K, 2: K* pi"};
      histos.add("MCGen/chargeChannel", "K1 parents in selected reconstructed events, inside the K1 rapidity window", HistType::kTH2D, {{2, -1.5, 1.5, "K1 charge"}, channelAxis});
      histos.add("MCGen/ptChannel", "Generated K1 pT by immediate decay channel", HistType::kTH2D, {channelAxis, ptAxis});
    }
  }

  EventCuts mEventCuts;
  TrackCuts mTrackCuts;
  PionPidCuts mPionPid;
  KaonPidCuts mKaonPid;
  SecondaryCuts mSecondaryCuts;
  CandidateCuts mCandidateCuts;
  HistogramOptions mHistogramOptions;
  LooseStageOptions mLooseOptions;

  // Derived once in init(): which candidate cuts are switched on.
  bool mSecondaryWindowOn = false;
  bool mAnotherMassCutOn = false;
  bool mPiKaMassCutOn = false;
  bool mAngleCutOn = false;
  bool mPairAsymCutOn = false;

  std::array<int, NTruthChannels> mTruthDebugCounts{};
};

} // namespace o2::analysis::k1micro

#endif // PWGLF_CORE_K1ANALYSISMICROCORE_H_
