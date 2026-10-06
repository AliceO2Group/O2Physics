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
/// \brief Shared K1(1270) selection, truth classification and candidate enumeration of the K1 resonance tasks
/// \author Su-Jeong Ji <su-jeong.ji@cern.ch>, Bong-Hwi Lim <bong-hwi.lim@cern.ch>
///
/// The core owns the selection and its cut-flow instrumentation (CutFlow/*, ML/*). The tasks own their
/// output histograms and fill them from the pair and candidate hooks of forEachCandidate().

#ifndef PWGLF_CORE_K1ANALYSISMICROCORE_H_
#define PWGLF_CORE_K1ANALYSISMICROCORE_H_

#include "PWGLF/Core/K1MlFeatures.h"
#include "PWGLF/Core/ResoAnalysisSelectionCore.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/Logger.h>

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
#include <set>
#include <type_traits>
#include <utility>
#include <vector>

namespace o2::analysis::k1micro
{

using o2::analysis::resonance::EventCuts;
using o2::analysis::resonance::isCutEnabled;
using o2::analysis::resonance::isInRange;
using o2::analysis::resonance::isInWindow;
using o2::analysis::resonance::PIDCutConfig;
using o2::analysis::resonance::ResoAnalysisSelectionCore;
using o2::analysis::resonance::TrackCuts;
using o2::analysis::resonance::TrackStage;

enum class K1TruthChannel {
  None = 0,
  RhoK = 1,
  KStarPi = 2
};

// The value is the species index of the PID configuration in ResoAnalysisSelectionCore.
enum class Species : int {
  Pion = 0,
  Kaon = 1
};

// Cumulative selection bits of an unlike-sign candidate handed to the export hook.
enum CandidatePassBit : uint16_t {
  kPassLoose = 1,     // valid canonical candidate inside the K1 rapidity window
  kPassQuality = 2,   // track quality of all three tracks
  kPassPID = 4,       // TOF requirement and PID of all three tracks
  kPassPair = 8,      // pion-pair pT and secondary mass window
  kPassCandidate = 16 // candidate cuts
};
inline constexpr uint16_t PassBitsSelected = kPassLoose | kPassQuality | kPassPID | kPassPair | kPassCandidate;

inline constexpr double MassRho770 = 0.77526; // PDG 2024, not available in o2::constants::physics
inline constexpr int NCandidateStages = 12;
inline constexpr int NTruthChannels = 3;

// Configurable groups without prefix: the JSON keys are the plain configurable names.

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

/// Process functions enabled in the task; they decide the configuration checks and the registered histograms.
struct ProcessModes {
  bool microTracks = false; // any process function reading micro tracks (quantised DCA and nSigma)
  bool mcReco = false;      // reconstructed MC with full or micro tracks
  bool mcRecoMicro = false; // reconstructed MC with micro tracks
  bool mcGen = false;       // generated K1 parents in selected reconstructed events
};

/// Loose-stage traversal of unlike-sign micro candidates for the export hook.
/// With the defaults and without an export hook, the candidate loop applies only the conventional selection.
struct LooseStageOptions {
  bool audit = false;          // fill ML/looseCutflow and ML/looseMassPtActivity
  bool exportSelected = false; // hand candidates to the export hook at the selected stage instead of the loose stage
};

/// Selected (pion, pion, kaon) triplet handed to the candidate hook.
/// mass13, mass23, angle and pairAsym are computed only inside the K1 rapidity window, and there only
/// when a candidate cut needs them or forEachCandidate() was asked for them; otherwise they are 0.
struct K1CandidateValues {
  ROOT::Math::PxPyPzMVector k1;        // pion1 + pion2 + kaon
  ROOT::Math::PxPyPzMVector secondary; // pion1 + pion2
  double mass13 = 0.;                  // pion1 + kaon (K*0 candidate in the unlike-sign case)
  double mass23 = 0.;                  // pion2 + kaon
  double angle = 0.;                   // opening angle between the secondary resonance and the bachelor
  double pairAsym = 0.;                // energy asymmetry between the secondary resonance and the bachelor
  bool isUnlikeSign = false;           // opposite-sign pion pair
  bool inRapidity = false;             // K1 rapidity window
  bool passesCandidateCuts = false;    // secondary mass, pi-K mass, angle and asymmetry cuts (inside the rapidity window)
};

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

/// K1 selection and candidate enumeration shared by the K1 tasks.
/// The task owns the configurable groups and the histogram registry and passes them in init().
class K1AnalysisMicroCore
{
 public:
  void init(o2::framework::HistogramRegistry& histos,
            EventCuts const& eventCuts, TrackCuts const& trackCuts,
            PionPidCuts const& pionPidCuts, KaonPidCuts const& kaonPidCuts,
            SecondaryCuts const& secondaryCuts, CandidateCuts const& candidateCuts,
            ProcessModes const& modes, LooseStageOptions const& looseOptions = {})
  {
    mSecondaryCuts = secondaryCuts;
    mCandidateCuts = candidateCuts;
    mLooseOptions = looseOptions;

    // The order follows Species
    std::vector<PIDCutConfig> pid(2);
    auto& pion = pid[static_cast<int>(Species::Pion)];
    pion.species = "Pion";
    pion.maxTPCnSigma = pionPidCuts.cMaxTPCnSigmaPion.value;
    pion.maxTOFnSigma = pionPidCuts.cMaxTOFnSigmaPion.value;
    pion.combinedNSigma = pionPidCuts.nsigmaCutCombinedPion.value;
    pion.onlyTOFTracks = pionPidCuts.cUseOnlyTOFTrackPi.value;
    pion.usePtDependent = pionPidCuts.cPionUsePtDepPID.value;
    pion.ptBins = pionPidCuts.cPionPIDPtBins.value;
    pion.tpcNSigmaCuts = pionPidCuts.cPionTPCNSigmaCuts.value;
    pion.tofNSigmaCuts = pionPidCuts.cPionTOFNSigmaCuts.value;
    pion.tofRequired = pionPidCuts.cPionTOFRequired.value;
    pion.maxTPCName = "cMaxTPCnSigmaPion";
    pion.maxTOFName = "cMaxTOFnSigmaPion";
    pion.tpcCutsName = "cPionTPCNSigmaCuts";
    pion.tofCutsName = "cPionTOFNSigmaCuts";
    auto& kaon = pid[static_cast<int>(Species::Kaon)];
    kaon.species = "Kaon";
    kaon.maxTPCnSigma = kaonPidCuts.cMaxTPCnSigmaKaon.value;
    kaon.maxTOFnSigma = kaonPidCuts.cMaxTOFnSigmaKaon.value;
    kaon.combinedNSigma = kaonPidCuts.nsigmaCutCombinedKaon.value;
    kaon.onlyTOFTracks = kaonPidCuts.cUseOnlyTOFTrackKa.value;
    kaon.usePtDependent = kaonPidCuts.cKaonUsePtDepPID.value;
    kaon.ptBins = kaonPidCuts.cKaonPIDPtBins.value;
    kaon.tpcNSigmaCuts = kaonPidCuts.cKaonTPCNSigmaCuts.value;
    kaon.tofNSigmaCuts = kaonPidCuts.cKaonTOFNSigmaCuts.value;
    kaon.tofRequired = kaonPidCuts.cKaonTOFRequired.value;
    kaon.maxTPCName = "cMaxTPCnSigmaKaon";
    kaon.maxTOFName = "cMaxTOFnSigmaKaon";
    kaon.tpcCutsName = "cKaonTPCNSigmaCuts";
    kaon.tofCutsName = "cKaonTOFNSigmaCuts";
    mSelection.init(eventCuts, trackCuts, std::move(pid), mCandidateCuts.cByPassTOF, modes.microTracks);

    mSecondaryWindowOn = isCutEnabled(mSecondaryCuts.cSecondaryMasswindow);
    mAnotherMassCutOn = isCutEnabled(mSecondaryCuts.cMinAnotherSecondaryMassCut) || isCutEnabled(mSecondaryCuts.cMaxAnotherSecondaryMassCut);
    mPiKaMassCutOn = isCutEnabled(mSecondaryCuts.cMinPiKaMassCut) || isCutEnabled(mSecondaryCuts.cMaxPiKaMassCut);
    mAngleCutOn = isCutEnabled(mSecondaryCuts.cMinAngle) || isCutEnabled(mSecondaryCuts.cMaxAngle);
    mPairAsymCutOn = isCutEnabled(mSecondaryCuts.cMinPairAsym) || isCutEnabled(mSecondaryCuts.cMaxPairAsym);

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

  // Full selection stage of a track (quality, TOF requirement, PID)
  template <bool IsResoMicrotrack, Species S, typename TrackType>
  int trackSelectionStage(const TrackType& track)
  {
    const int qualityStage = mSelection.trackQualityStage<IsResoMicrotrack>(track);
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
    if constexpr (S == Species::Pion) {
      if (hasTOF) {
        tofNSigma = track.tofNSigmaPi();
      }
    } else {
      if (hasTOF) {
        tofNSigma = track.tofNSigmaKa();
      }
    }
    if (!mSelection.passesPID<IsResoMicrotrack>(SpeciesIndex, track.pt(), hasTOF, tpcNSigma, tofNSigma)) {
      return TrackStage::kTrkTOFRequired;
    }
    return TrackStage::kTrkPID;
  }

  // Unordered (pion, pion, kaon) candidate enumeration of one collision (or one mixed pair of collisions).
  // dTracks1: bachelor kaons, dTracks2: pions. The core fills the cut-flow instrumentation; the hooks
  // (nullptr to skip) receive
  //  - onPair(trk1, trk2, secondary, passesPairPt): every pion pair passing the pion selection,
  //    trk1 being the pion with the lower index,
  //  - onCandidate(kaon, pion1, pion2, K1CandidateValues): every triplet passing the pion, pair and kaon
  //    selection, with the canonical pion roles (pion1: opposite sign to the kaon in the unlike-sign case),
  //  - onExport(collision, kaon, same-sign pion, opposite-sign pion, truth channel, pass bits): the unlike-sign
  //    micro same-event candidates at the loose or selected stage configured by LooseStageOptions.
  // computeValues requests mass13, mass23, angle and pairAsym for every candidate in the rapidity window.
  template <bool IsMC, bool IsMix, bool IsResoMicrotrack, typename CollisionType, typename TracksType,
            typename PairHook = std::nullptr_t, typename CandidateHook = std::nullptr_t, typename ExportHook = std::nullptr_t>
  void forEachCandidate(o2::framework::HistogramRegistry& histos, const CollisionType& collision, const TracksType& dTracks1, const TracksType& dTracks2,
                        bool computeValues, PairHook onPair = nullptr, CandidateHook onCandidate = nullptr, ExportHook onExport = nullptr)
  {
    if (dTracks1.size() == 0 || dTracks2.size() == 0) {
      return;
    }
    constexpr bool HasPairHook = !std::is_same_v<PairHook, std::nullptr_t>;
    constexpr bool HasCandidateHook = !std::is_same_v<CandidateHook, std::nullptr_t>;
    constexpr bool HasExportHook = !std::is_same_v<ExportHook, std::nullptr_t>;
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
    // The canonical-candidate validity is also required before a selected-stage export.
    bool visitLoose = false;
    bool checkValidity = false;
    if constexpr (IsResoMicrotrack && !IsMix) {
      visitLoose = mLooseOptions.audit || (!mLooseOptions.exportSelected && HasExportHook);
      checkValidity = visitLoose || (mLooseOptions.exportSelected && HasExportHook);
    }
    std::vector<uint8_t> kaonQuality(dTracks1.size(), 0), pionQuality(dTracks2.size(), 0);
    if (visitLoose) {
      for (const auto& track : dTracks1) {
        kaonQuality[getCacheIndex(track, firstKaonIndex, kaonQuality.size())] = mSelection.trackQualityStage<IsResoMicrotrack>(track) == TrackStage::kTrkClusters;
      }
      for (const auto& track : dTracks2) {
        pionQuality[getCacheIndex(track, firstPionIndex, pionQuality.size())] = mSelection.trackQualityStage<IsResoMicrotrack>(track) == TrackStage::kTrkClusters;
      }
    }

    // Values needed only by switched-on cuts or by the task are computed only then
    const bool isK892Mode = mSecondaryCuts.cfgModeK892orRho;
    const bool needAngle = computeValues || mAngleCutOn;
    const bool needPairAsym = computeValues || mPairAsymCutOn;
    // K892 mode: the K* candidate is (trk1, K), rho mode: the rho is (trk1, trk2)
    const bool needMass13 = computeValues || (isK892Mode ? mSecondaryWindowOn : mAnotherMassCutOn) || (isK892Mode && (needAngle || needPairAsym));
    const bool needMass23 = computeValues || mPiKaMassCutOn;
    const bool rhoWindowOn = mSecondaryWindowOn && !isK892Mode;

    ROOT::Math::PxPyPzMVector lDecayDaughter1, lDecayDaughter2, lResonanceSecondary, lDecayDaughter_bach, lResonanceK1, lPair13, lPair23;
    // Unordered pion pairs: each (pion, pion, kaon) triplet is visited once.
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
      if constexpr (HasPairHook) {
        if (pionsSelected) {
          onPair(trk1, trk2, lResonanceSecondary, pairPt);
        }
      }
      if (!pairPt && !visitLoose) {
        continue;
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
        const bool tripletSelected = pairSelected && bachelorSelected;
        if (tripletSelected) {
          countCandidate(7);
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
        K1CandidateValues values;
        values.k1 = lResonanceK1;
        values.secondary = lResonanceSecondary;
        values.isUnlikeSign = isUnlikeSign;

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
        values.inRapidity = lResonanceK1.Rapidity() >= mCandidateCuts.cK1MinRap && lResonanceK1.Rapidity() <= mCandidateCuts.cK1MaxRap;
        if (!values.inRapidity) {
          if constexpr (HasCandidateHook) {
            if (tripletSelected) {
              onCandidate(bTrack, pion1, pion2, values);
            }
          }
          continue;
        }
        if (tripletSelected) {
          countCandidate(8);
        }

        if (needMass13) {
          lPair13 = lPion1 + lDecayDaughter_bach;
          values.mass13 = lPair13.M();
        }
        if (needMass23) {
          lPair23 = lPion2 + lDecayDaughter_bach;
          values.mass23 = lPair23.M();
        }
        // Rho mode: secondary = (trk1, trk2) against the bachelor. K892 mode: secondary = (trk1, K) against trk2.
        if (needAngle) {
          values.angle = isK892Mode ? ROOT::Math::VectorUtil::Angle(lPair13, lPion2) : ROOT::Math::VectorUtil::Angle(lResonanceSecondary, lDecayDaughter_bach);
        }
        if (needPairAsym) {
          values.pairAsym = isK892Mode ? (lPair13.E() - lPion2.E()) / (lPair13.E() + lPion2.E())
                                       : (lResonanceSecondary.E() - lDecayDaughter_bach.E()) / (lResonanceSecondary.E() + lDecayDaughter_bach.E());
        }

        // Candidate cuts (each one is evaluated only if switched on)
        values.passesCandidateCuts =
          (!isK892Mode || !mSecondaryWindowOn || (isInWindow(values.mass13, o2::constants::physics::MassK0Star892, mSecondaryCuts.cSecondaryMasswindow) && pion1.sign() != bTrack.sign())) &&
          (!mAnotherMassCutOn || isInRange(isK892Mode ? lResonanceSecondary.M() : values.mass13, mSecondaryCuts.cMinAnotherSecondaryMassCut, mSecondaryCuts.cMaxAnotherSecondaryMassCut)) &&
          (!mPiKaMassCutOn || isInRange(values.mass23, mSecondaryCuts.cMinPiKaMassCut, mSecondaryCuts.cMaxPiKaMassCut)) &&
          (!mAngleCutOn || isInRange(values.angle, mSecondaryCuts.cMinAngle, mSecondaryCuts.cMaxAngle)) &&
          (!mPairAsymCutOn || isInRange(values.pairAsym, mSecondaryCuts.cMinPairAsym, mSecondaryCuts.cMaxPairAsym));
        auto exportCandidate = [&](uint16_t passBits) {
          if constexpr (IsResoMicrotrack && !IsMix && HasExportHook) {
            onExport(collision, bTrack, pion2, pion1, flowChannel, passBits);
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
                if (values.passesCandidateCuts) {
                  passBits |= kPassCandidate;
                  countMl(5);
                }
              }
            }
          }
          if (!mLooseOptions.exportSelected) {
            exportCandidate(passBits);
          }
        }
        // Stage C retains the frozen conventional selections.
        if (!tripletSelected) {
          continue;
        }
        if (values.passesCandidateCuts) {
          countCandidate(9);
          countCandidate(isUnlikeSign ? 10 : 11);
          if (isUnlikeSign && mLooseOptions.exportSelected && validLoose) {
            exportCandidate(PassBitsSelected);
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
        }
        if constexpr (HasCandidateHook) {
          onCandidate(bTrack, pion1, pion2, values);
        }
      } // bTrack
    }
  } // forEachCandidate

  // Generated K1 parents of a selected reconstructed MC collision. The optional callback receives
  // (parent, immediate channel) for the parents inside the K1 rapidity window.
  // Parents belong to selected reconstructed events; split reco collisions
  // repeat parent sets. This is not an unconditional generated denominator.
  template <typename ParentsType, typename Callback = std::nullptr_t>
  void forEachGeneratedK1(o2::framework::HistogramRegistry& histos, const ParentsType& resoParents, Callback callback = nullptr)
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
      selected[getCacheIndex(track, firstIndex, selected.size())] = (stage == TrackStage::kTrkPID) ? 1 : 0;
      if constexpr (FillCutFlow) {
        for (int i = 0; i <= stage; ++i) {
          histos.fill(HIST("CutFlow/tracks"), i, static_cast<int>(S));
        }
      }
    }
    return selected;
  }

  // Cut-flow instrumentation of the selection; the output histograms belong to the tasks.
  void registerHistograms(o2::framework::HistogramRegistry& histos, ProcessModes const& modes)
  {
    using o2::framework::HistType;
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
    constexpr int NTrackStages = TrackStage::kTrkNStages;
    auto trackFlow = histos.add<TH2>("CutFlow/tracks", "Micro tracks, once per selected collision;stage;species", HistType::kTH2D, {{NTrackStages, -0.5, NTrackStages - 0.5}, {2, -0.5, 1.5}});
    const std::array<const char*, NTrackStages> trackLabels{"input", "pT", "eta", "DCAxy", "DCAz", "track flags", "clusters / crossed rows", "TOF required", "PID"};
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
  }

  ResoAnalysisSelectionCore mSelection;
  SecondaryCuts mSecondaryCuts;
  CandidateCuts mCandidateCuts;
  LooseStageOptions mLooseOptions;

  // Derived once in init(): which candidate cuts are switched on.
  bool mSecondaryWindowOn = false;
  bool mAnotherMassCutOn = false;
  bool mPiKaMassCutOn = false;
  bool mAngleCutOn = false;
  bool mPairAsymCutOn = false;
};

} // namespace o2::analysis::k1micro

#endif // PWGLF_CORE_K1ANALYSISMICROCORE_H_
