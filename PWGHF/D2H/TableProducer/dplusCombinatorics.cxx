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

/// \file dplusCombinatorics.cxx
/// \brief D± → π± K∓ π± and D±s → K± K∓ π± secondary-vertex finding driven by track-level ML
///
/// Standalone alternative to the 3-prong branch of trackIndexSkimCreator.cxx. The linear,
/// pT-binned single-track selections of the skim creator are replaced by *two* multi-class BDTs,
/// one per species, each scoring every track as {not-a-daughter, pion-from-X, kaon-from-X}.
/// Inference runs once per (track, collision) association, so a tighter track selection
/// translates directly into fewer triplets entering the vertex fitter.
///
/// Two tasks live here:
///  - HfDplusTrackSelMl   : builds the track features and either dumps them for training (MC)
///                          or evaluates both models and flags pion-/kaon-like tracks per species
///  - HfDplusCombinatorics: the triple loop over the flagged tracks, filling aod::Hf3Prongs
///
/// The combinatorics offers three strategies, so that the cost of each can be measured:
///  - processDplus        : D+ only, the original implementation, kept so that the earlier
///                          efficiency/rate comparison against the skim creator stays reproducible
///  - processSeparateLoops: one triple loop per species, each over its own track sets. A triplet
///                          that satisfies both hypotheses is fitted -- and written -- twice.
///  - processMergedLoop   : a single triple loop over the union of the selected tracks, with a
///                          per-triplet hypothesis mask. A triplet satisfying both hypotheses is
///                          fitted once and written once, with both flags set. This is what the
///                          standard skim creator does, and whether it pays off depends entirely
///                          on how much the two species' track sets overlap.
///
/// Both new strategies cache the propagated track parameters once per collision instead of
/// re-propagating inside the innermost loop, so that they differ only in the loop structure.
///
/// The efficiency/rate QA used to compare this workflow with the standard track-index skim
/// creator lives in a separate executable, dplusCombinatoricsQa.cxx, so that it can be piped
/// after either skimmer without two devices declaring aod::Hf3Prongs as output.
///
/// \author Fabrizio Chinu <fabrizio.chinu@cern.ch>, Universita and INFN Torino

#include "PWGHF/Core/HfMlResponseHfTracks.h"
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/AliasTables.h"
#include "PWGHF/DataModel/CharmTrackScoreScanTables.h"
#include "PWGHF/DataModel/TrackIndexSkimmingTables.h"
#include "PWGHF/Utils/utilsAnalysis.h"
#include "PWGHF/Utils/utilsBfieldCCDB.h"
#include "PWGHF/Utils/utilsEvSelHf.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/ZorroSummary.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Tools/ML/MlResponse.h"

#include <CCDB/BasicCCDBManager.h>
#include <CCDB/CcdbApi.h>
#include <CommonConstants/PhysicsConstants.h>
#include <DCAFitter/DCAFitterN.h>
#include <DetectorsBase/MatLayerCylSet.h>
#include <DetectorsBase/Propagator.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Array2D.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/DCA.h>
#include <ReconstructionDataFormats/Track.h>
#include <ReconstructionDataFormats/Vertex.h>

#include <TH1.h>
#include <TPDGCode.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::analysis;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::hf_centrality;
using namespace o2::hf_evsel;
using namespace o2::constants::physics;

namespace o2::aod
{
namespace hf_skim_sel_collision
{
DECLARE_SOA_COLUMN(WhyRejectColl, whyRejectColl, o2::hf_evsel::HfCollisionRejectionMask); //! reason for rejection (0 = selected)
DECLARE_SOA_COLUMN(Centrality, centrality, float);                                        //! centrality percentile of the collision
} // namespace hf_skim_sel_collision

DECLARE_SOA_TABLE(HfSkimSelColls, "AOD", "HFSKIMSELCOLL", //! Collision selection for the charm combinatorics
                  hf_skim_sel_collision::WhyRejectColl,
                  hf_skim_sel_collision::Centrality);

namespace hf_skim_track_sel
{
// D0 -> K pi model
DECLARE_SOA_COLUMN(ScoreBkgDzero, scoreBkgDzero, float);   //! ML score, background class
DECLARE_SOA_COLUMN(ScorePionDzero, scorePionDzero, float); //! ML score, pion-from-D0 class
DECLARE_SOA_COLUMN(ScoreKaonDzero, scoreKaonDzero, float); //! ML score, kaon-from-D0 class
DECLARE_SOA_COLUMN(IsSelPionDzero, isSelPionDzero, bool);  //! passes the D0 pion-score threshold
DECLARE_SOA_COLUMN(IsSelKaonDzero, isSelKaonDzero, bool);  //! passes the D0 kaon-score threshold
  // D+ -> K pi pi model
DECLARE_SOA_COLUMN(ScoreBkgDplus, scoreBkgDplus, float);   //! ML score, background class
DECLARE_SOA_COLUMN(ScorePionDplus, scorePionDplus, float); //! ML score, pion-from-D+ class
DECLARE_SOA_COLUMN(ScoreKaonDplus, scoreKaonDplus, float); //! ML score, kaon-from-D+ class
DECLARE_SOA_COLUMN(IsSelPionDplus, isSelPionDplus, bool);  //! passes the D+ pion-score threshold
DECLARE_SOA_COLUMN(IsSelKaonDplus, isSelKaonDplus, bool);  //! passes the D+ kaon-score threshold
// Ds -> K K pi model
DECLARE_SOA_COLUMN(ScoreBkgDs, scoreBkgDs, float);   //! ML score, background class
DECLARE_SOA_COLUMN(ScorePionDs, scorePionDs, float); //! ML score, pion-from-Ds class
DECLARE_SOA_COLUMN(ScoreKaonDs, scoreKaonDs, float); //! ML score, kaon-from-Ds class
DECLARE_SOA_COLUMN(IsSelPionDs, isSelPionDs, bool);  //! passes the Ds pion-score threshold
DECLARE_SOA_COLUMN(IsSelKaonDs, isSelKaonDs, bool);  //! passes the Ds kaon-score threshold
// D* -> D0 pi -> pi K pi model
DECLARE_SOA_COLUMN(ScoreBkgDstar, scoreBkgDstar, float);   //! ML score, background class
DECLARE_SOA_COLUMN(ScorePionDstar, scorePionDstar, float); //! ML score, pion-from-D* class
DECLARE_SOA_COLUMN(ScoreKaonDstar, scoreKaonDstar, float); //! ML score, kaon-from-D* class
DECLARE_SOA_COLUMN(IsSelPionDstar, isSelPionDstar, bool);  //! passes the D* pion-score threshold
DECLARE_SOA_COLUMN(IsSelKaonDstar, isSelKaonDstar, bool);  //! passes the D* kaon-score threshold
// common
DECLARE_SOA_COLUMN(IsPositive, isPositive, bool);  //! charge sign of the track
DECLARE_SOA_COLUMN(IsSelAny, isSelAny, bool);      //! passes at least one of the four selections
} // namespace hf_skim_track_sel

DECLARE_SOA_TABLE(HfSkimTrackSels, "AOD", "HFSKIMTRKSEL", //! Track-level ML selection, aligned with aod::TrackAssoc
                  hf_skim_track_sel::ScoreBkgDzero,
                  hf_skim_track_sel::ScorePionDzero,
                  hf_skim_track_sel::ScoreKaonDzero,
                  hf_skim_track_sel::IsSelPionDzero,
                  hf_skim_track_sel::IsSelKaonDzero,
                  hf_skim_track_sel::ScoreBkgDplus,
                  hf_skim_track_sel::ScorePionDplus,
                  hf_skim_track_sel::ScoreKaonDplus,
                  hf_skim_track_sel::IsSelPionDplus,
                  hf_skim_track_sel::IsSelKaonDplus,
                  hf_skim_track_sel::ScoreBkgDs,
                  hf_skim_track_sel::ScorePionDs,
                  hf_skim_track_sel::ScoreKaonDs,
                  hf_skim_track_sel::IsSelPionDs,
                  hf_skim_track_sel::IsSelKaonDs,
                  hf_skim_track_sel::ScoreBkgDstar,
                  hf_skim_track_sel::ScorePionDstar,
                  hf_skim_track_sel::ScoreKaonDstar,
                  hf_skim_track_sel::IsSelPionDstar,
                  hf_skim_track_sel::IsSelKaonDstar,
                  hf_skim_track_sel::IsPositive,
                  hf_skim_track_sel::IsSelAny);

namespace hf_skim_track_sel_mask
{
DECLARE_SOA_COLUMN(CollIdx, collIdx, int);      //! global index of the collision the association refers to
DECLARE_SOA_COLUMN(SelMask, selMask, uint16_t);  //! bit0 -> pi(D0)
                                                 //! bit1 -> K(D0)
                                                 //! bit2 -> pi(D+)
                                                 //! bit3 -> K(D+)
                                                 //! bit4 -> pi(Ds)
                                                 //! bit5 -> K(Ds)
                                                 //! bit6 -> pi(D*)
                                                 //! bit7 -> K(D*)
                                                 //! bit8 -> positive charge
} // namespace hf_skim_track_sel_mask

/// One row per *selected* (track, collision) association. This is the smallest object from which
/// the overlap between the four role selections, and the number of triplets each combinatorics
/// strategy would enumerate, can be reconstructed exactly offline: the 16 mask multiplicities per
/// (collision, charge) are a sufficient statistic for both.
DECLARE_SOA_TABLE(HfSkimTrackMasks, "AOD", "HFSKIMTRKMASK", //! Compact dump of the per-track role selections
                  hf_skim_track_sel_mask::CollIdx,
                  hf_skim_track_sel_mask::SelMask);

namespace hf_skim_track_train
{
// features
DECLARE_SOA_COLUMN(Pt, pt, float);                                                     //!
DECLARE_SOA_COLUMN(Eta, eta, float);                                                   //!
DECLARE_SOA_COLUMN(DcaXY, dcaXY, float);                                               //!
DECLARE_SOA_COLUMN(DcaZ, dcaZ, float);                                                 //!
DECLARE_SOA_COLUMN(SigmaDcaXY, sigmaDcaXY, float);                                     //!
DECLARE_SOA_COLUMN(SigmaDcaZ, sigmaDcaZ, float);                                       //!
DECLARE_SOA_COLUMN(NormDcaXY, normDcaXY, float);                                       //!
DECLARE_SOA_COLUMN(NormDcaZ, normDcaZ, float);                                         //!
DECLARE_SOA_COLUMN(Signed1Pt, signed1Pt, float);                                       //!
DECLARE_SOA_COLUMN(Tgl, tgl, float);                                                   //!
DECLARE_SOA_COLUMN(Sign, sign, float);                                                 //!
DECLARE_SOA_COLUMN(IsPvContributor, isPvContributor, float);                           //!
DECLARE_SOA_COLUMN(ItsNCls, itsNCls, float);                                           //!
DECLARE_SOA_COLUMN(ItsNClsInnerBarrel, itsNClsInnerBarrel, float);                     //!
DECLARE_SOA_COLUMN(ItsChi2NCl, itsChi2NCl, float);                                     //!
DECLARE_SOA_COLUMN(TpcNClsFound, tpcNClsFound, float);                                 //!
DECLARE_SOA_COLUMN(TpcNClsCrossedRows, tpcNClsCrossedRows, float);                     //!
DECLARE_SOA_COLUMN(TpcCrossedRowsOverFindableCls, tpcCrossedRowsOverFindableCls, float); //!
DECLARE_SOA_COLUMN(TpcChi2NCl, tpcChi2NCl, float);                                     //!
DECLARE_SOA_COLUMN(TpcFractionSharedCls, tpcFractionSharedCls, float);                 //!
DECLARE_SOA_COLUMN(TpcNSigmaPi, tpcNSigmaPi, float);                                   //!
DECLARE_SOA_COLUMN(TpcNSigmaKa, tpcNSigmaKa, float);                                   //!
DECLARE_SOA_COLUMN(Centrality, centrality, float);                                     //!
// labels and auxiliary information, not meant to be used as features
DECLARE_SOA_COLUMN(TrackIndex, trackIndex, int64_t);                     //! global index of the track: the same track can be associated to several collisions, and therefore appear more than once, so the train/test split has to group by it
DECLARE_SOA_COLUMN(McMask, mcMask, uint16_t);                            //! MC label quality mask (bit 15 = fake track)
DECLARE_SOA_COLUMN(HasCharmBeautyAncestor, hasCharmBeautyAncestor, bool); //! any charm/beauty hadron in the ancestry
// Selection outcomes
DECLARE_SOA_COLUMN(Label, label, int8_t);       //! D+ model: -1 = not usable, 0 = bkg, 1 = pion from D+, 2 = kaon from D+
DECLARE_SOA_COLUMN(PtCand, ptCand, float);      //! pT of the mother D+ (-1 if none)
DECLARE_SOA_COLUMN(IsPrompt, isPrompt, int8_t); //! 1 prompt, 0 non-prompt, -1 not applicable
DECLARE_SOA_COLUMN(Channel, channel, int8_t);   //! -1 none, 0 D+ -> K pi pi, 1 D+ -> anti-K*0 pi
DECLARE_SOA_COLUMN(PdgMother, pdgMother, int);  //! PDG code of the mother particle
} // namespace hf_skim_track_train

DECLARE_SOA_TABLE(HfSkimTrackTrains, "AOD", "HFSKIMTRKTRAIN", //! Track-level training sample, one row per usable association
                  hf_skim_track_train::Pt,
                  hf_skim_track_train::Eta,
                  hf_skim_track_train::DcaXY,
                  hf_skim_track_train::DcaZ,
                  hf_skim_track_train::SigmaDcaXY,
                  hf_skim_track_train::SigmaDcaZ,
                  hf_skim_track_train::NormDcaXY,
                  hf_skim_track_train::NormDcaZ,
                  hf_skim_track_train::Signed1Pt,
                  hf_skim_track_train::Tgl,
                  hf_skim_track_train::Sign,
                  hf_skim_track_train::IsPvContributor,
                  hf_skim_track_train::ItsNCls,
                  hf_skim_track_train::ItsNClsInnerBarrel,
                  hf_skim_track_train::ItsChi2NCl,
                  hf_skim_track_train::TpcNClsFound,
                  hf_skim_track_train::TpcNClsCrossedRows,
                  hf_skim_track_train::TpcCrossedRowsOverFindableCls,
                  hf_skim_track_train::TpcChi2NCl,
                  hf_skim_track_train::TpcFractionSharedCls,
                  hf_skim_track_train::TpcNSigmaPi,
                  hf_skim_track_train::TpcNSigmaKa,
                  hf_skim_track_train::Centrality,
                  hf_skim_track_train::TrackIndex,
                  hf_skim_track_train::Label,
                  hf_skim_track_train::PtCand,
                  hf_skim_track_train::IsPrompt,
                  hf_skim_track_train::Channel,
                  hf_skim_track_train::PdgMother,
                  hf_skim_track_train::McMask,
                  hf_skim_track_train::HasCharmBeautyAncestor);
} // namespace o2::aod

/// Bit set in aod::McTrackLabels::mcMask when the MC label of the track is not reliable.
constexpr uint16_t McMaskFakeTrackBit = 15;

// HF mothers to investigate
constexpr std::array<int, 4> PdgMotherList = {Pdg::kD0, Pdg::kDPlus, Pdg::kDS, Pdg::kDStar};

/// Decay topologies used to identify the channels in the MC kinematic tree.
constexpr int NDaughtersTwoBody = 2;
constexpr int NDaughtersDzeroToPiK = 2;
constexpr int NPionsDzeroToPiK = 1;
constexpr int NKaonsDzeroToPiK = 1;
constexpr int NDaughtersDplusToPiKPi = 3;
constexpr int NPionsDplusToPiKPi = 2;
constexpr int NKaonsDplusToPiKPi = 1;
constexpr int NDaughtersDsToKKPi = 3;
constexpr int NKaonsDsToKKPi = 2;
constexpr int NPionsDsToKKPi = 1;
constexpr int NDaughtersDstarToD0Pi = 2;
constexpr int ND0DstarToD0Pi = 1;
constexpr int NPionsDstarToD0Pi = 1;

/// Resonant sub-channels recorded next to the training label, so that the signal definition can be
/// narrowed offline without re-running the workflow.
enum Channels : int {
  kNone = -1,        ///< no mother found, or mother not in the list of interest
  kDzeroDirect,
  kDplusDirect,
  kDplusKstar,
  kDsPhiPi,
  kDsKstarK,
  kDsDirect,
  kDstarToD0Pi,
};
// constexpr int8_t ChannelNone = -1;
// constexpr int8_t ChannelDzeroDirect = 0;  ///< D0 -> K- pi+
// constexpr int8_t ChannelDplusDirect = 0;  ///< D+ -> K- pi+ pi+
// constexpr int8_t ChannelDplusKstar = 1;   ///< D+ -> anti-K*(892)0 pi+, anti-K*0 -> K- pi+
// constexpr int8_t ChannelDsPhiPi = 0;      ///< Ds+ -> phi(1020) pi+, phi -> K+ K-
// constexpr int8_t ChannelDsKstarK = 1;     ///< Ds+ -> anti-K*(892)0 K+, anti-K*0 -> K- pi+
// constexpr int8_t ChannelDsDirect = 2;     ///< Ds+ -> K+ K- pi+ non resonant

/// Above this value the chi2 cut on the two-track vertex is considered disabled.
constexpr double MaxTwoTrackChi2PcaDisabled = 1.e9;

/// Margin applied to the kinematic cuts of the cheap pre-filter that runs before the track is
/// re-propagated. Propagation without material corrections preserves q/pT and tgl, so pT and eta
/// are in practice unchanged; the margin only absorbs numerical drift. The exact cut is applied
/// afterwards, on the propagated values.
constexpr float PreselKinematicMargin = 0.1f;

/// Bits of the compact per-association role mask, both in aod::HfSkimTrackMasks and in the
/// per-collision prong cache of the combinatorics.
/// Output-class count of the origin-aware (5-class) track model.
constexpr int NClassesWithOrigin = 5;

enum TrackRoleBit : uint16_t {
  BitPionDzero = 0, ///< pion-like for D0 -> K pi
  BitKaonDzero,     ///< kaon-like for D0 -> K pi
  BitPionDplus,     ///< pion-like for D+ -> K pi pi
  BitKaonDplus,     ///< kaon-like for D+ -> K pi pi
  BitPionDs,        ///< pion-like for Ds -> K K pi
  BitKaonDs,        ///< kaon-like for Ds -> K K pi
  BitPionDstar,     ///< pion-like for D* -> D0 pi
  BitKaonDstar,     ///< kaon-like for D* -> D0 pi
  BitPositive       ///< charge sign, so that one row fully describes the association
};

enum TrackType : uint8_t {
  Pion,
  Kaon
};

/// CPU accounting. steady_clock::now() costs ~20 ns through the vDSO, so the timers are placed
/// only around blocks that themselves cost microseconds -- a model evaluation, a vertex fit, one
/// collision's worth of looping -- and never inside the per-triplet enumeration, whose own cost
/// is comparable to a clock read. Everything not attributed to a measured block is reported as a
/// remainder, so the bins always add up to the total.
using StopWatch = std::chrono::steady_clock;

inline double elapsedSeconds(StopWatch::time_point const& since)
{
  return std::chrono::duration<double>(StopWatch::now() - since).count();
}

/// Timing bins of hf-dplus-track-sel-ml.
enum TimingTrackSel : int {
  TimeTrackSelTotal = 0,
  TimeTrackSelEvSel,
  TimeTrackSelFeatures,
  TimeTrackSelInferenceDzero,
  TimeTrackSelInferenceDplus,
  TimeTrackSelInferenceDs,
  TimeTrackSelInferenceDstar,
  TimeTrackSelRemainder,
  NTimingTrackSel
};

/// Timing bins of hf-dplus-combinatorics.
enum TimingCombinatorics : int {
  TimeCombTotal = 0,
  TimeCombCache,
  TimeCombFit2Prong,
  TimeCombFit3Prong,
  TimeCombRemainder,
  NTimingCombinatorics
};

/// Counter bins of hf-dplus-track-sel-ml.
enum CounterTrackSel : int {
  CountTrackSelCollisions = 0,
  CountAssociations,
  CountCheapPresel,
  CountPresel,
  CountInferenceDzero,
  CountInferenceDplus,
  CountInferenceDs,
  CountInferenceDstar,
  CountSelPionDzero,
  CountSelKaonDzero,
  CountSelPionDplus,
  CountSelKaonDplus,
  CountSelPionDs,
  CountSelKaonDs,
  CountSelPionDstar,
  CountSelKaonDstar,
  CountSelAny,
  NCounterTrackSel
};

/// Outcomes of the heavy-flavour ancestry walk of the written training rows (hHfAncestryWalk).
enum HfAncestryWalkOutcome : int {
  HfWalkNoMcParticle = 0,
  HfWalkNotFound,
  HfWalkFound,
  HfWalkDepthLimit,
  NHfAncestryWalkOutcomes
};

/// Counter bins of hf-dplus-combinatorics.
enum CounterCombinatorics : int {
  CountCollisions = 0,
  CountTriplets,      ///< triplets whose hypothesis mask is non-empty, i.e. reaching the preselection
  CountTripletsPresel,///< triplets passing the mass/pT preselection, i.e. reaching the 3-prong fitter
  CountFit2Prong,
  CountFit3Prong,
  CountRows,          ///< rows written to aod::Hf3Prongs
  CountRowsDplus,
  CountRowsDs,
  CountRowsBoth,
  NCounterCombinatorics
};

/// Deterministic pseudo-random number in [0, 1) from a 64-bit index.
/// Used to downsample the background class reproducibly, without any RNG state.
inline float deterministicUniform(uint64_t index, uint64_t seed)
{
  // splitmix64
  uint64_t z = index + seed + 0x9e3779b97f4a7c15ULL;
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
  z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
  z = z ^ (z >> 31);
  return static_cast<float>(z >> 11) / static_cast<float>(1ULL << 53);
}

/// Track-level feature building, training-sample production and ML inference for both species.
struct HfDplusTrackSelMl {
  Produces<aod::HfSkimSelColls> rowSelectedCollision;
  Produces<aod::HfSkimTrackSels> rowSelectedTrack;
  Produces<aod::HfSkimTrackMasks> rowSelMask;
  Produces<aod::HfSkimTrackTrains> rowTrackTrain;

  struct : ConfigurableGroup {
    Configurable<bool> fillHistograms{"fillHistograms", true, "fill QA histograms"};
    Configurable<bool> enableTiming{"enableTiming", false, "account the CPU time of the individual stages into hTiming"};
    Configurable<bool> fillSelMaskTable{"fillSelMaskTable", false, "write one aod::HfSkimTrackMasks row per selected association (for the offline overlap study)"};
    Configurable<bool> rejBkgCollisions{"rejBkgCollisions", true, "Reject tracks from background events"};
    // loose pre-filter applied before inference (kept minimal on purpose)
    Configurable<float> ptMinTrack{"ptMinTrack", 0.3f, "min. track pT for the charm combinatorics"};
    Configurable<float> etaMaxTrack{"etaMaxTrack", 0.8f, "max. |eta| for the charm combinatorics"};
    // Run 3 equivalent of the old ITS-refit requirement. It is NOT spelled
    //   TESTBIT(track.flags(), o2::aod::track::ITSrefit)
    // the way the track-index skim creator does it: ITSrefit belongs to TrackFlagsRun2Enum and
    // equals 0x1, so TESTBIT(flags, 0x1) tests bit *index* 1, i.e. mask 0x2 -- which in Run 3 is
    // PVContributor, not ITSrefit. That silently required charm daughters to contribute to the
    // primary vertex: verified on this sample, every selected track had isPvContributor == 1,
    // and it removed ~40% of the genuine signal, preferentially the displaced tracks the
    // selection is supposed to keep.
    Configurable<bool> requireIts{"requireIts", true, "require the track to have ITS information"};
    Configurable<bool> requireItsIbHit{"requireItsIbHit", true, "require a hit in one of the two innermost ITS layers"};
    // D0 model
    Configurable<bool> applyMlDzero{"applyMlDzero", true, "evaluate the D+ track model"};
    Configurable<int> nClassesDzero{"nClassesDzero", 3, "output classes of the D+ model: 3 (bkg/pi/K) or 5 (bkg/prompt pi/prompt K/non-prompt pi/non-prompt K)"};
    Configurable<std::vector<double>> binsPtMlDzero{"binsPtMlDzero", std::vector<double>{0., 1.e10}, "pT bin limits of the D0 track models"};
    Configurable<std::vector<std::string>> onnxFileNamesDzero{"onnxFileNamesDzero", std::vector<std::string>{"ModelHandler_DzeroTracks.onnx"}, "ONNX file names of the D0 models, one per pT bin"};
    Configurable<std::vector<std::string>> inputFeaturesDzero{"inputFeaturesDzero", std::vector<std::string>{"pt", "eta", "dcaXY", "dcaZ", "sigmaDcaXY", "sigmaDcaZ", "normDcaXY", "normDcaZ", "signed1Pt", "tgl", "sign", "isPvContributor", "itsNCls", "itsNClsInnerBarrel", "itsChi2NCl", "tpcNClsFound", "tpcCrossedRowsOverFindableCls", "tpcChi2NCl", "tpcFractionSharedCls", "tpcNSigmaPi", "tpcNSigmaKa"}, "input features of the D0 models, in the order expected by them"};
    Configurable<std::vector<double>> thresholdScorePionDzero{"thresholdScorePionDzero", std::vector<double>{0.}, "min. D0 pion-class score, one per pT bin"};
    Configurable<std::vector<double>> thresholdScoreKaonDzero{"thresholdScoreKaonDzero", std::vector<double>{0.}, "min. D0 kaon-class score, one per pT bin"};
    // D+ model
    Configurable<bool> applyMlDplus{"applyMlDplus", true, "evaluate the D+ track model"};
    Configurable<int> nClassesDplus{"nClassesDplus", 3, "output classes of the D+ model: 3 (bkg/pi/K) or 5 (bkg/prompt pi/prompt K/non-prompt pi/non-prompt K)"};
    Configurable<std::vector<double>> binsPtMlDplus{"binsPtMlDplus", std::vector<double>{0., 1.e10}, "pT bin limits of the D+ track models"};
    Configurable<std::vector<std::string>> onnxFileNamesDplus{"onnxFileNamesDplus", std::vector<std::string>{"ModelHandler_DplusTracks.onnx"}, "ONNX file names of the D+ models, one per pT bin"};
    Configurable<std::vector<std::string>> inputFeaturesDplus{"inputFeaturesDplus", std::vector<std::string>{"pt", "eta", "dcaXY", "dcaZ", "sigmaDcaXY", "sigmaDcaZ", "normDcaXY", "normDcaZ", "signed1Pt", "tgl", "sign", "isPvContributor", "itsNCls", "itsNClsInnerBarrel", "itsChi2NCl", "tpcNClsFound", "tpcCrossedRowsOverFindableCls", "tpcChi2NCl", "tpcFractionSharedCls", "tpcNSigmaPi", "tpcNSigmaKa"}, "input features of the D+ models, in the order expected by them"};
    Configurable<std::vector<double>> thresholdScorePionDplus{"thresholdScorePionDplus", std::vector<double>{0.}, "min. D+ pion-class score, one per pT bin"};
    Configurable<std::vector<double>> thresholdScoreKaonDplus{"thresholdScoreKaonDplus", std::vector<double>{0.}, "min. D+ kaon-class score, one per pT bin"};
    // Ds model
    Configurable<bool> applyMlDs{"applyMlDs", true, "evaluate the Ds track model"};
    Configurable<int> nClassesDs{"nClassesDs", 3, "output classes of the Ds model: 3 or 5, as for the D+"};
    Configurable<std::vector<double>> binsPtMlDs{"binsPtMlDs", std::vector<double>{0., 1.e10}, "pT bin limits of the Ds track models"};
    Configurable<std::vector<std::string>> onnxFileNamesDs{"onnxFileNamesDs", std::vector<std::string>{"ModelHandler_DsTracks.onnx"}, "ONNX file names of the Ds models, one per pT bin"};
    Configurable<std::vector<std::string>> inputFeaturesDs{"inputFeaturesDs", std::vector<std::string>{"pt", "eta", "dcaXY", "dcaZ", "sigmaDcaXY", "sigmaDcaZ", "normDcaXY", "normDcaZ", "signed1Pt", "tgl", "sign", "isPvContributor", "itsNCls", "itsNClsInnerBarrel", "itsChi2NCl", "tpcNClsFound", "tpcCrossedRowsOverFindableCls", "tpcChi2NCl", "tpcFractionSharedCls", "tpcNSigmaPi", "tpcNSigmaKa"}, "input features of the Ds models, in the order expected by them"};
    Configurable<std::vector<double>> thresholdScorePionDs{"thresholdScorePionDs", std::vector<double>{0.}, "min. Ds pion-class score, one per pT bin"};
    Configurable<std::vector<double>> thresholdScoreKaonDs{"thresholdScoreKaonDs", std::vector<double>{0.}, "min. Ds kaon-class score, one per pT bin"};
    // Ds model
    Configurable<bool> applyMlDstar{"applyMlDstar", true, "evaluate the Dstar track model"};
    Configurable<int> nClassesDstar{"nClassesDstar", 2, "output classes of the Dstar model: 2, bkg or pion"};
    Configurable<std::vector<double>> binsPtMlDstar{"binsPtMlDstar", std::vector<double>{0., 1.e10}, "pT bin limits of the Dstar track models"};
    Configurable<std::vector<std::string>> onnxFileNamesDstar{"onnxFileNamesDstar", std::vector<std::string>{"ModelHandler_DstarTracks.onnx"}, "ONNX file names of the Dstar models, one per pT bin"};
    Configurable<std::vector<std::string>> inputFeaturesDstar{"inputFeaturesDstar", std::vector<std::string>{"pt", "eta", "dcaXY", "dcaZ", "sigmaDcaXY", "sigmaDcaZ", "normDcaXY", "normDcaZ", "signed1Pt", "tgl", "sign", "isPvContributor", "itsNCls", "itsNClsInnerBarrel", "itsChi2NCl", "tpcNClsFound", "tpcCrossedRowsOverFindableCls", "tpcChi2NCl", "tpcFractionSharedCls", "tpcNSigmaPi", "tpcNSigmaKa"}, "input features of the Dstar models, in the order expected by them"};
    Configurable<std::vector<double>> thresholdScorePionDstar{"thresholdScorePionDstar", std::vector<double>{0.}, "min. Dstar pion-class score, one per pT bin"};
    // ONNX runtime
    Configurable<bool> loadModelsFromCcdb{"loadModelsFromCcdb", false, "load the ONNX models from CCDB instead of a local path"};
    Configurable<std::string> mlModelPathCcdb{"mlModelPathCcdb", "Users/f/fchinu/DplusTracks", "CCDB path of the ML models"};
    Configurable<int64_t> timestampCcdbForMlModels{"timestampCcdbForMlModels", -1, "timestamp of the ONNX files to be queried in CCDB"};
    Configurable<bool> enableOnnxOptimizations{"enableOnnxOptimizations", true, "enable the ONNX graph optimisations"};
    Configurable<int> onnxThreads{"onnxThreads", 1, "number of threads used by the ONNX runtime (0 = let onnxruntime decide)"};
    // training-sample production
    Configurable<float> bkgSamplingFraction{"bkgSamplingFraction", 0.05f, "fraction of background tracks written to the training tree"};
    Configurable<uint64_t> bkgSamplingSeed{"bkgSamplingSeed", 12345, "seed of the deterministic background downsampling"};
    Configurable<bool> keepFakeTracksInBkg{"keepFakeTracksInBkg", true, "keep tracks with an unreliable MC label in the background class"};
    // RecoDecay::getMother walks the ancestry breadth-first, keeping every stage in a vector of
    // vectors and deduplicating each stage with a linear std::find. With an unbounded depth that
    // walk runs back to the beam particles: in Pb-Pb the fan-out makes it quadratic and enormous
    // (tens of GB of retained heap, minutes per timeframe). It is physically wrong too -- through
    // the fragmentation record almost any track in an event containing a D+ acquires a D+
    // "ancestor". A real charm decay product is a direct daughter, or sits behind one resonance
    // or one strange decay, so depth 3 covers every genuine chain.
    Configurable<int> maxAncestryDepth{"maxAncestryDepth", 3, "max. ancestry depth for the same-species veto and the prompt/non-prompt flag"};
    // The charm/beauty ancestor flag walks the decay chain only (it stops at partons, see
    // continueAncestryWalk), so it does not suffer from the fan-out above and can go much deeper.
    Configurable<int> maxAncestryDepthHfFlag{"maxAncestryDepthHfFlag", 30, "max. ancestry depth of the charm/beauty ancestor flag (reaching it sets the flag)"};
    // CCDB
    Configurable<bool> isRun2{"isRun2", false, "enable Run 2 or Run 3 GRP objects for the magnetic field"};
    Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
    Configurable<std::string> ccdbPathLut{"ccdbPathLut", "GLO/Param/MatLUT", "path for LUT parametrization"};
    Configurable<std::string> ccdbPathGrp{"ccdbPathGrp", "GLO/GRP/GRP", "path of the grp file (Run 2)"};
    Configurable<std::string> ccdbPathGrpMag{"ccdbPathGrpMag", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object (Run 3)"};
  } config;

  HfEventSelection hfEvSel;                   // event selection and monitoring, shared with the skim creator
  Service<o2::ccdb::BasicCCDBManager> ccdb{}; // needed for the event selection and the magnetic field
  o2::ccdb::CcdbApi ccdbApi;
  o2::base::MatLayerCylSet* lut{};
  o2::base::Propagator::MatCorrType noMatCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;
  int runNumber{};

  /// One track-level model: the response, its per-pT-bin score thresholds and whether it runs.
  struct TrackModel {
    HfMlResponseHfTracks<float> response{};
    std::vector<double> thresholdPion{};
    std::vector<double> thresholdKaon{};
    bool enabled{false};
    /// 3 = bkg / pion / kaon, 5 = bkg / prompt pi / prompt K / non-prompt pi / non-prompt K.
    /// With 5 the two origin scores of a role are summed into that role's score, so the working
    /// point stays one threshold per role and the two model kinds are directly comparable.
    int nClasses{3};
  };
  TrackModel modelDzero;
  TrackModel modelDplus;
  TrackModel modelDs;
  TrackModel modelDstar;

  // running variables for the re-propagation of tracks to a non-default collision
  o2::dataformats::DCA dcaInfoCov;
  o2::dataformats::VertexBase vtx;

  // per-collision timing accumulators, in seconds
  std::array<double, NTimingTrackSel> timing{};

  using CollisionsWithCent = soa::Join<aod::Collisions, aod::EvSels, aod::PVMults, aod::CentFT0Cs>; // PVMults: INEL>0 flag of the HF event selection
  using TracksWithSelAndPid = soa::Join<aod::TracksWCovDcaExtra, aod::TracksDCACov, aod::TrackSelection, aod::pidTPCFullPi, aod::pidTPCFullKa>;
  using TracksWithSelAndPidMc = soa::Join<aod::TracksWCovDcaExtra, aod::TracksDCACov, aod::TrackSelection, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::McTrackLabels>;

  Preslice<aod::TrackAssoc> trackIndicesPerCollision = aod::track_association::collisionId;

  HistogramRegistry registry{"registry"};
  // OutputObj<ZorroSummary> zorroSummary{"zorroSummary"};

  /// Configure one model from its configurables. Returns false if the model is disabled.
  /// \param model is the model to configure
  /// \param binsPt is the pT binning
  /// \param onnxFiles are the ONNX file names, one per pT bin
  /// \param features are the input feature names, in the order the model expects
  /// \param thrPion, thrKaon are the per-bin score thresholds
  /// \param name is used in the error messages
  bool configureModel(TrackModel& model,
                      std::vector<double> const& binsPt,
                      std::vector<std::string> const& onnxFiles,
                      std::vector<std::string> const& features,
                      std::vector<double> const& thrPion,
                      std::vector<double> const& thrKaon,
                      int nClasses,
                      const char* name)
  {
    if (nClasses != static_cast<int>(HfTrackMlClass::NClasses) && nClasses != NClassesWithOrigin) {
      LOGP(fatal, "{}: nClasses must be {} or {}, got {}!", name,
           static_cast<int>(HfTrackMlClass::NClasses), NClassesWithOrigin, nClasses);
    }
    model.nClasses = nClasses;
    const std::size_t nModels = binsPt.size() - 1;
    if (onnxFiles.size() != nModels) {
      LOGP(fatal, "{}: number of ONNX files ({}) does not match the number of pT bins ({})!", name, onnxFiles.size(), nModels);
    }
    if (thrPion.size() != nModels || thrKaon.size() != nModels) {
      LOGP(fatal, "{}: number of score thresholds (pion {}, kaon {}) does not match the number of pT bins ({})!", name, thrPion.size(), thrKaon.size(), nModels);
    }
    model.thresholdPion = thrPion;
    model.thresholdKaon = thrKaon;

    // the per-class thresholds are applied by the task, so no cut is delegated to MlResponse
    const auto nOut = static_cast<std::size_t>(nClasses);
    const std::vector<int> cutDir(nOut, o2::cuts_ml::CutDirection::CutNot);
    const std::vector<double> dummyCutValues(nModels * nOut, 0.);
    const o2::framework::LabeledArray<double> dummyCuts{dummyCutValues.data(), static_cast<uint32_t>(nModels), static_cast<uint32_t>(nOut), {}, {}};

    model.response.configure(binsPt, dummyCuts, cutDir, static_cast<uint8_t>(nOut));
    if (config.loadModelsFromCcdb) {
      ccdbApi.init(config.ccdbUrl);
      model.response.setModelPathsCCDB(onnxFiles, ccdbApi, std::vector<std::string>(nModels, config.mlModelPathCcdb.value), config.timestampCcdbForMlModels);
    } else {
      model.response.setModelPathsLocal(onnxFiles);
    }
    model.response.cacheInputFeaturesIndices(features);
    model.response.init(config.enableOnnxOptimizations, config.onnxThreads);
    model.enabled = true;
    LOGP(info, "{}: {} model(s) configured with {} input features and {} output classes", name, nModels, features.size(), nClasses);
    return true;
  }

  void init(InitContext const&)
  {
    const std::array<int, 3> doProcess = {doprocessTraining, doprocessTrainingData, doprocessApplication};
    if (std::accumulate(doProcess.begin(), doProcess.end(), 0) != 1) {
      LOGP(fatal, "One and only one process function of HfDplusTrackSelMl must be enabled at a time!");
    }

    // hfEvSel.init(registry, &zorroSummary);
    hfEvSel.init(registry);

    ccdb->setURL(config.ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(config.ccdbPathLut));

    if (config.fillHistograms) {
      const AxisSpec axisPt{200, 0., 20., "#it{p}_{T}^{track} (GeV/#it{c})"};
      const AxisSpec axisScore{100, 0., 1., "ML score"};
      const AxisSpec axisNTracks{500, 0., 5000., "n. tracks"};

      registry.add("hPtAllTracks", "all associated tracks;#it{p}_{T}^{track} (GeV/#it{c});entries", {HistType::kTH1D, {axisPt}});
      registry.add("hPtPreselected", "tracks passing the loose pre-filter;#it{p}_{T}^{track} (GeV/#it{c});entries", {HistType::kTH1D, {axisPt}});

      if (doprocessApplication) {
        registry.add("hScoreBkgDzero", "D^{0} background score;score;entries", {HistType::kTH1D, {axisScore}});
        registry.add("hScorePionDzero", "D^{0} pion score;score;entries", {HistType::kTH1D, {axisScore}});
        registry.add("hScoreKaonDzero", "D^{0} kaon score;score;entries", {HistType::kTH1D, {axisScore}});
        registry.add("hPtSelPionDzero", "D^{0} pion-like tracks;#it{p}_{T}^{track} (GeV/#it{c});entries", {HistType::kTH1D, {axisPt}});
        registry.add("hPtSelKaonDzero", "D^{0} kaon-like tracks;#it{p}_{T}^{track} (GeV/#it{c});entries", {HistType::kTH1D, {axisPt}});
        registry.add("hNTracksSelPionDzeroPerColl", "D^{0} pion-like tracks per collision;n. tracks;entries", {HistType::kTH1D, {axisNTracks}});
        registry.add("hNTracksSelKaonDzeroPerColl", "D^{0} kaon-like tracks per collision;n. tracks;entries", {HistType::kTH1D, {axisNTracks}});
        registry.add("hScoreBkgDplus", "D^{+} background score;score;entries", {HistType::kTH1D, {axisScore}});
        registry.add("hScorePionDplus", "D^{+} pion score;score;entries", {HistType::kTH1D, {axisScore}});
        registry.add("hScoreKaonDplus", "D^{+} kaon score;score;entries", {HistType::kTH1D, {axisScore}});
        registry.add("hPtSelPionDplus", "D^{+} pion-like tracks;#it{p}_{T}^{track} (GeV/#it{c});entries", {HistType::kTH1D, {axisPt}});
        registry.add("hPtSelKaonDplus", "D^{+} kaon-like tracks;#it{p}_{T}^{track} (GeV/#it{c});entries", {HistType::kTH1D, {axisPt}});
        registry.add("hNTracksSelPionDplusPerColl", "D^{+} pion-like tracks per collision;n. tracks;entries", {HistType::kTH1D, {axisNTracks}});
        registry.add("hNTracksSelKaonDplusPerColl", "D^{+} kaon-like tracks per collision;n. tracks;entries", {HistType::kTH1D, {axisNTracks}});
        registry.add("hScoreBkgDs", "D_{s} background score;score;entries", {HistType::kTH1D, {axisScore}});
        registry.add("hScorePionDs", "D_{s} pion score;score;entries", {HistType::kTH1D, {axisScore}});
        registry.add("hScoreKaonDs", "D_{s} kaon score;score;entries", {HistType::kTH1D, {axisScore}});
        registry.add("hPtSelPionDs", "D_{s} pion-like tracks;#it{p}_{T}^{track} (GeV/#it{c});entries", {HistType::kTH1D, {axisPt}});
        registry.add("hPtSelKaonDs", "D_{s} kaon-like tracks;#it{p}_{T}^{track} (GeV/#it{c});entries", {HistType::kTH1D, {axisPt}});
        registry.add("hNTracksSelPionDsPerColl", "D_{s} pion-like tracks per collision;n. tracks;entries", {HistType::kTH1D, {axisNTracks}});
        registry.add("hNTracksSelKaonDsPerColl", "D_{s} kaon-like tracks per collision;n. tracks;entries", {HistType::kTH1D, {axisNTracks}});
        registry.add("hScoreBkgDstar", "D^{*} background score;score;entries", {HistType::kTH1D, {axisScore}});
        registry.add("hScorePionDstar", "D^{*} pion score;score;entries", {HistType::kTH1D, {axisScore}});
        registry.add("hScoreKaonDstar", "D^{*} kaon score;score;entries", {HistType::kTH1D, {axisScore}});
        registry.add("hPtSelPionDstar", "D^{*} pion-like tracks;#it{p}_{T}^{track} (GeV/#it{c});entries", {HistType::kTH1D, {axisPt}});
        registry.add("hPtSelKaonDstar", "D^{*} kaon-like tracks;#it{p}_{T}^{track} (GeV/#it{c});entries", {HistType::kTH1D, {axisPt}});
        registry.add("hNTracksSelPionDstarPerColl", "D^{*} pion-like tracks per collision;n. tracks;entries", {HistType::kTH1D, {axisNTracks}});
        registry.add("hNTracksSelKaonDstarPerColl", "D^{*} kaon-like tracks per collision;n. tracks;entries", {HistType::kTH1D, {axisNTracks}});
        registry.add("hNTracksPreselPerColl", "tracks passing the pre-filter per collision;n. tracks;entries", {HistType::kTH1D, {axisNTracks}});
        registry.add("hNTracksSelAnyPerColl", "tracks passing at least one role selection;n. tracks;entries", {HistType::kTH1D, {axisNTracks}});
        // Role overlap, summed over the whole sample. The 16 bins are the multiplicities of the
        // four-bit role mask: this is what tells whether the two species should share one triple
        // loop or run two of their own.
        registry.add("hRoleMask", "role mask of the selected associations;mask (bit 0 #pi D^{+}, 1 K D^{+}, 2 #pi D_{s}, 3 K D_{s});entries", {HistType::kTH1D, {{16, -0.5, 15.5}}});
      }
      if (doprocessTraining) {
        registry.add("hLabelDzero", "D^{0} training sample;label;entries", {HistType::kTH1D, {{3, -0.5, 2.5}}});
        registry.get<TH1>(HIST("hLabelDzero"))->GetXaxis()->SetBinLabel(1, "background");
        registry.get<TH1>(HIST("hLabelDzero"))->GetXaxis()->SetBinLabel(2, "#pi from D^{0}");
        registry.get<TH1>(HIST("hLabelDzero"))->GetXaxis()->SetBinLabel(3, "K from D^{0}");
        registry.add("hLabelDplus", "D^{+} training sample;label;entries", {HistType::kTH1D, {{3, -0.5, 2.5}}});
        registry.get<TH1>(HIST("hLabelDplus"))->GetXaxis()->SetBinLabel(1, "background");
        registry.get<TH1>(HIST("hLabelDplus"))->GetXaxis()->SetBinLabel(2, "#pi from D^{+}");
        registry.get<TH1>(HIST("hLabelDplus"))->GetXaxis()->SetBinLabel(3, "K from D^{+}");
        registry.add("hLabelDs", "D_{s} training sample;label;entries", {HistType::kTH1D, {{3, -0.5, 2.5}}});
        registry.get<TH1>(HIST("hLabelDs"))->GetXaxis()->SetBinLabel(1, "background");
        registry.get<TH1>(HIST("hLabelDs"))->GetXaxis()->SetBinLabel(2, "#pi from D_{s}");
        registry.get<TH1>(HIST("hLabelDs"))->GetXaxis()->SetBinLabel(3, "K from D_{s}");
        registry.add("hLabelDstar", "D^{*} training sample;label;entries", {HistType::kTH1D, {{3, -0.5, 2.5}}});
        registry.get<TH1>(HIST("hLabelDstar"))->GetXaxis()->SetBinLabel(1, "background");
        registry.get<TH1>(HIST("hLabelDstar"))->GetXaxis()->SetBinLabel(2, "#pi from D^{*}");
        registry.get<TH1>(HIST("hLabelDstar"))->GetXaxis()->SetBinLabel(3, "K from D^{*}");
        registry.add("hChannelDzero", "resonant sub-channel of the labelled D^{0};channel;entries", {HistType::kTH1D, {{2, -0.5, 1.5}}});
        registry.add("hChannelDplus", "resonant sub-channel of the labelled D^{+};channel;entries", {HistType::kTH1D, {{2, -0.5, 1.5}}});
        registry.add("hChannelDs", "resonant sub-channel of the labelled D_{s};channel;entries", {HistType::kTH1D, {{3, -0.5, 2.5}}});
        registry.add("hChannelDstar", "resonant sub-channel of the labelled D^{*};channel;entries", {HistType::kTH1D, {{2, -0.5, 1.5}}});
        float hfAncestryWalkOutcomes = static_cast<float>(NHfAncestryWalkOutcomes);
        registry.add("hHfAncestryWalk", "charm/beauty ancestry of the written rows;;entries", {HistType::kTH1D, {{hfAncestryWalkOutcomes, -0.5, hfAncestryWalkOutcomes - 0.5}}});
        const char* walkLabels[NHfAncestryWalkOutcomes] = {"no MC particle", "no HF ancestor", "HF ancestor", "depth limit"};
        for (int i = 0; i < NHfAncestryWalkOutcomes; ++i) {
          registry.get<TH1>(HIST("hHfAncestryWalk"))->GetXaxis()->SetBinLabel(i + 1, walkLabels[i]);
        }
      }
    }

    // CPU accounting. Always booked, so that a run with the timers off is recognisable by an
    // empty histogram rather than by a missing one.
    registry.add("hTiming", "CPU time by stage;;CPU time (s)", {HistType::kTH1D, {{NTimingTrackSel, -0.5, NTimingTrackSel - 0.5}}});
    const char* timingLabels[NTimingTrackSel] = {"total", "event selection", "features", "inference D0", "inference D+", "inference Ds", "inference Dstar", "remainder"};
    for (int i = 0; i < NTimingTrackSel; ++i) {
      registry.get<TH1>(HIST("hTiming"))->GetXaxis()->SetBinLabel(i + 1, timingLabels[i]);
    }
    registry.add("hCounters", "counters;;entries", {HistType::kTH1D, {{NCounterTrackSel, -0.5, NCounterTrackSel - 0.5}}});
    const char* counterLabels[NCounterTrackSel] = {"collisions", "associations", "cheap presel", "presel", "inference D0", "inference D+", "inference Ds", "inference Dstar",
                                                   "#pi D^{0}", "K D^{0}", "#pi D^{+}", "K D^{+}", "#pi D_{s}", "K D_{s}", "#pi D^{*}", "K D^{*}", "any"};
    for (int i = 0; i < NCounterTrackSel; ++i) {
      registry.get<TH1>(HIST("hCounters"))->GetXaxis()->SetBinLabel(i + 1, counterLabels[i]);
    }

    if (doprocessApplication) {
      if (config.applyMlDzero) {
        configureModel(modelDzero, config.binsPtMlDzero, config.onnxFileNamesDzero, config.inputFeaturesDzero,
                       config.thresholdScorePionDzero, config.thresholdScoreKaonDzero,
                       config.nClassesDzero, "D0 track model");
      }
      if (config.applyMlDplus) {
        configureModel(modelDplus, config.binsPtMlDplus, config.onnxFileNamesDplus, config.inputFeaturesDplus,
                       config.thresholdScorePionDplus, config.thresholdScoreKaonDplus,
                       config.nClassesDplus, "D+ track model");
      }
      if (config.applyMlDs) {
        configureModel(modelDs, config.binsPtMlDs, config.onnxFileNamesDs, config.inputFeaturesDs,
                       config.thresholdScorePionDs, config.thresholdScoreKaonDs,
                       config.nClassesDs, "Ds track model");
      }
      if (config.applyMlDstar) {
        configureModel(modelDstar, config.binsPtMlDstar, config.onnxFileNamesDstar, config.inputFeaturesDstar,
                       config.thresholdScorePionDstar, std::vector<double>{0.},
                       config.nClassesDstar, "Dstar track model");
      }
      if (!modelDzero.enabled && !modelDplus.enabled && !modelDs.enabled && !modelDstar.enabled) {
        LOGP(fatal, "processApplication is enabled but neither model is: nothing would be selected.");
      }
    }
  }

  /// Propagation-independent part of the preselection, evaluated on the un-propagated track.
  /// Re-propagating a track to a collision that is not its own costs a full covariance
  /// propagation, so it must not be paid for tracks that no charm analysis would ever use. The
  /// ITS requirements here are exact (propagation cannot change them); the kinematic ones carry a
  /// margin and are re-applied exactly, on the propagated values, by isPreselected.
  /// \param track is the track
  /// \return true if the track is worth propagating
  template <typename TTrack>
  bool passesCheapPreselection(TTrack const& track)
  {
    if (config.requireIts && !track.hasITS()) {
      return false;
    }
    if (config.requireItsIbHit) {
      const auto clusterMap = track.itsClusterMap();
      if (!(TESTBIT(clusterMap, 0) || TESTBIT(clusterMap, 1))) {
        return false;
      }
    }
    if (track.pt() < config.ptMinTrack * (1.f - PreselKinematicMargin)) {
      return false;
    }
    if (std::abs(track.eta()) > config.etaMaxTrack + PreselKinematicMargin) {
      return false;
    }
    return true;
  }

  template <typename TTrack>
  bool isPreselected(TTrack const& track, float trackPt, float trackEta)
  {
    if (trackPt < config.ptMinTrack) {
      return false;
    }
    if (std::abs(trackEta) > config.etaMaxTrack) {
      return false;
    }
    if (config.requireIts && !track.hasITS()) {
      return false;
    }
    if (config.requireItsIbHit) {
      const auto clusterMap = track.itsClusterMap();
      if (!(TESTBIT(clusterMap, 0) || TESTBIT(clusterMap, 1))) {
        return false;
      }
    }
    return true;
  }

  /// Fill the feature container for one (track, collision) association.
  /// Tracks that do not belong to the collision under study are re-propagated to it, so that the
  /// impact parameter and its resolution refer to the correct primary vertex.
  /// \param collision is the collision the track is associated to
  /// \param track is the track
  /// \param centrality is the centrality percentile of the collision
  /// \param features is the container to be filled
  template <typename TCollision, typename TTrack>
  void fillFeatures(TCollision const& collision, TTrack const& track, float centrality, HfTrackMlFeatures& features)
  {
    float trackPt = track.pt();
    float trackEta = track.eta();
    float dcaXY = track.dcaXY();
    float dcaZ = track.dcaZ();
    float sigmaDcaXY = std::sqrt(std::abs(track.sigmaDcaXY2()));
    float sigmaDcaZ = std::sqrt(std::abs(track.sigmaDcaZ2()));

    if (track.collisionId() != collision.globalIndex()) {
      // not the default collision of this track: re-propagate to get DCA and its covariance
      auto trackParCov = getTrackParCov(track);
      dcaInfoCov.set(999.f, 999.f, 999.f, 999.f, 999.f);
      vtx.setPos({collision.posX(), collision.posY(), collision.posZ()});
      vtx.setCov(collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ());
      if (o2::base::Propagator::Instance()->propagateToDCABxByBz(vtx, trackParCov, 2.f, noMatCorr, &dcaInfoCov)) {
        trackPt = trackParCov.getPt();
        trackEta = trackParCov.getEta();
        dcaXY = dcaInfoCov.getY();
        dcaZ = dcaInfoCov.getZ();
        sigmaDcaXY = std::sqrt(std::abs(dcaInfoCov.getSigmaY2()));
        sigmaDcaZ = std::sqrt(std::abs(dcaInfoCov.getSigmaZ2()));
      }
    }

    constexpr float MinSigmaDca = 1.e-6f; // guard against a null resolution before dividing

    features.pt = trackPt;
    features.eta = trackEta;
    features.dcaXY = dcaXY;
    features.dcaZ = dcaZ;
    features.sigmaDcaXY = sigmaDcaXY;
    features.sigmaDcaZ = sigmaDcaZ;
    features.normDcaXY = dcaXY / std::max(sigmaDcaXY, MinSigmaDca);
    features.normDcaZ = dcaZ / std::max(sigmaDcaZ, MinSigmaDca);
    features.signed1Pt = track.signed1Pt();
    features.tgl = track.tgl();
    features.sign = static_cast<float>(track.sign());
    features.isPvContributor = static_cast<float>(track.isPVContributor());
    features.itsNCls = static_cast<float>(track.itsNCls());
    features.itsNClsInnerBarrel = static_cast<float>(track.itsNClsInnerBarrel());
    features.itsChi2NCl = track.itsChi2NCl();
    features.tpcNClsFound = static_cast<float>(track.tpcNClsFound());
    features.tpcNClsCrossedRows = static_cast<float>(track.tpcNClsCrossedRows());
    features.tpcCrossedRowsOverFindableCls = track.tpcCrossedRowsOverFindableCls();
    features.tpcChi2NCl = track.tpcChi2NCl();
    features.tpcFractionSharedCls = track.tpcFractionSharedCls();
    features.tpcNSigmaPi = track.tpcNSigmaPi();
    features.tpcNSigmaKa = track.tpcNSigmaKa();
    features.centrality = centrality;
  }

  /// Index of the closest ancestor with the requested |PDG| code, searched breadth-first up to
  /// maxDepth generations, or -1 if there is none.
  ///
  /// Deliberately NOT RecoDecay::getMother. That helper iterates the mother indices as a *range*
  ///     for (i = particleMother.mothersIds().front(); i <= particleMother.mothersIds().back(); ++i)
  /// but mothersIds() is an array of the actual mother indices, not a contiguous slice. A hadron
  /// produced in string fragmentation has two mothers -- the distant endpoints of the string -- so
  /// that range spans hundreds of unrelated particles (median 650, up to 1880 in this sample) and
  /// any D+ sitting inside the window is returned as an "ancestor". Measured on AO2D_002, this
  /// mislabelled 2869 fragmentation pions and kaons as D+ daughters against 1594 genuine ones,
  /// pushing the pion:kaon ratio to 4.5 where the K pi pi topology demands 2.0.
  /// The same pattern appears in getCharmHadronOrigin, so the prompt/non-prompt flag is derived
  /// here too rather than through that helper.
  /// \param particle is the MC particle to start from
  /// \param mcParticles is the MC particle table
  /// \param pdgWanted is the |PDG| code looked for (particles and antiparticles both match)
  /// \param maxDepth is the maximum number of generations to walk
  /// \return the global index of the matched ancestor, or -1
  template <typename TMcParticle, typename TMcParticles>
  int findAncestorWithPdg(TMcParticle const& particle, TMcParticles const& mcParticles, int pdgWanted, int maxDepth)
  {
    std::vector<int64_t> current{particle.globalIndex()};
    std::vector<int64_t> next{};
    for (int depth = 0; depth < maxDepth && !current.empty(); ++depth) {
      next.clear();
      for (const auto idx : current) { // o2-linter: disable=const-ref-in-for-loop (int elements)
        const auto thisParticle = mcParticles.rawIteratorAt(idx - mcParticles.offset());
        if (!thisParticle.has_mothers()) {
          continue;
        }
        for (const auto motherId : thisParticle.mothersIds()) { // o2-linter: disable=const-ref-in-for-loop (int elements)
          if (motherId < 0) {
            continue;
          }
          const auto mother = mcParticles.rawIteratorAt(motherId - mcParticles.offset());
          if (std::abs(mother.pdgCode()) == pdgWanted) {
            return motherId;
          }
          next.push_back(motherId);
        }
      }
      current.swap(next);
    }
    return -1;
  }

  /// Whether a particle decays into exactly the two requested species.
  /// \param particle is the mother
  /// \param pdgA, pdgB are the |PDG| codes of the two daughters (they may be equal)
  /// \return true if the decay is a two-body one into pdgA pdgB
  template <typename TMcParticle, typename TMcParticles>
  bool decaysToTwoBody(TMcParticle const& particle, TMcParticles const&, int pdgA, int pdgB)
  {
    int nDaughters{0}, nA{0}, nB{0};
    for (const auto& daughter : particle.template daughters_as<TMcParticles>()) {
      ++nDaughters;
      const int pdgDaughter = std::abs(daughter.pdgCode());
      if (pdgDaughter == pdgA) {
        ++nA;
      }
      if (pdgDaughter == pdgB) {
        ++nB;
      }
    }
    if (nDaughters != NDaughtersTwoBody) {
      return false;
    }
    return (pdgA == pdgB) ? (nA == NDaughtersTwoBody) : (nA == 1 && nB == 1);
  }

  /// Sub-channel of a D0 decay into the K pi final state, or ChannelNone if it is not one.
  /// \param mother is the D0 MC particle
  /// \param mcParticles is the MC particle table
  /// \return ChannelDzeroDirect or ChannelNone
  template <typename TMcParticle, typename TMcParticles>
  int8_t getDzeroToKPiChannel(TMcParticle const& mother, TMcParticles const&) // mcParticles)
  {
    LOG(info) << "Running getDzeroToKPiChannel";
    int nDaughters{0}, nPions{0}, nKaons{0}, nOther{0};
    for (const auto& daughter : mother.template daughters_as<TMcParticles>()) {
      ++nDaughters;
      LOG(info) << "D0 daughter PDG code: " << daughter.pdgCode();
      switch (std::abs(daughter.pdgCode())) {
        case kPiPlus:
          ++nPions;
          break;
        case kKPlus:
          ++nKaons;
          break;
        default:
          break;
      }
    }
    LOG(info) << "D0 daughters: " << nDaughters << ", pions: " << nPions << ", kaons: " << nKaons << ", other: " << nOther;
    if (nDaughters == NDaughtersDzeroToPiK && nPions == NPionsDzeroToPiK && nKaons == NKaonsDzeroToPiK) {
      LOG(info) << "D0 decay channel identified as Channels::kDzeroDirect";
      return Channels::kDzeroDirect;
    }
    LOG(info) << "D0 decay channel not identified, returning Channels::kNone";
    return Channels::kNone;
  }

  /// Sub-channel of a D0 decay into the K pi final state, or ChannelNone if it is not one.
  /// \param mother is the D* MC particle
  /// \param mcParticles is the MC particle table
  /// \return ChannelDstarToD0Pi or ChannelNone
  template <typename TMcParticle, typename TMcParticles>
  int8_t getDstarToD0PiChannel(TMcParticle const& mother, TMcParticles const&) // mcParticles)
  {
    LOG(info) << "Running getDstarToD0PiChannel";
    int nDaughters{0}, nPions{0}, nD0{0}, nOther{0};
    for (const auto& daughter : mother.template daughters_as<TMcParticles>()) {
      ++nDaughters;
      LOG(info) << "D* daughter PDG code: " << daughter.pdgCode();
      switch (std::abs(daughter.pdgCode())) {
        case kPiPlus:
          ++nPions;
          break;
        case kD0:
          ++nD0;
          break;
        default:
          ++nOther;
          break;
      }
    }
    LOG(info) << "D* daughters: " << nDaughters << ", pions: " << nPions << ", D0: " << nD0 << ", other: " << nOther;
    if (nDaughters == NDaughtersDstarToD0Pi && nPions == NPionsDstarToD0Pi && nD0 == ND0DstarToD0Pi) {
      LOG(info) << "D* decay channel identified as Channels::kDstarToD0Pi";
      return Channels::kDstarToD0Pi;
    }
    LOG(info) << "D* decay channel not identified, returning Channels::kNone";
    return Channels::kNone;
  }

  /// Sub-channel of a D+ decay into the K pi pi final state, or Channels::kNone if it is not one.
  /// \param mother is the D+ MC particle
  /// \param mcParticles is the MC particle table
  /// \return Channels::kDplusDirect, Channels::kDplusKstar or Channels::kNone
  template <typename TMcParticle, typename TMcParticles>
  int8_t getDplusToKPiPiChannel(TMcParticle const& mother, TMcParticles const& mcParticles)
  {
    int nDaughters{0}, nPions{0}, nKaons{0}, nKstars{0};
    int64_t indexKstar{-1};
    for (const auto& daughter : mother.template daughters_as<TMcParticles>()) {
      ++nDaughters;
      switch (std::abs(daughter.pdgCode())) {
        case kPiPlus:
          ++nPions;
          break;
        case kKPlus:
          ++nKaons;
          break;
        case Pdg::kK0Star892:
          ++nKstars;
          indexKstar = daughter.globalIndex();
          break;
        default:
          break;
      }
    }
    if (nDaughters == NDaughtersDplusToPiKPi && nPions == NPionsDplusToPiKPi && nKaons == NKaonsDplusToPiKPi) {
      return Channels::kDplusDirect;
    }
    if (nDaughters == NDaughtersTwoBody && nPions == 1 && nKstars == 1) {
      const auto kstar = mcParticles.rawIteratorAt(indexKstar - mcParticles.offset());
      if (decaysToTwoBody(kstar, mcParticles, kKPlus, kPiPlus)) {
        return Channels::kDplusKstar;
      }
    }
    return Channels::kNone;
  }

  /// Sub-channel of a Ds decay into the K K pi final state, or Channels::kNone if it is not one.
  /// The three sub-channels share the same final state -- two kaons and one pion -- so the role of
  /// a given track follows from its own PDG code, whatever the intermediate resonance.
  /// \param mother is the Ds MC particle
  /// \param mcParticles is the MC particle table
  /// \return Channels::kDsPhiPi, Channels::kDsKstarK, Channels::kDsDirect or Channels::kNone
  template <typename TMcParticle, typename TMcParticles>
  int8_t getDsToKKPiChannel(TMcParticle const& mother, TMcParticles const& mcParticles)
  {
    int nDaughters{0}, nPions{0}, nKaons{0}, nPhis{0}, nKstars{0};
    int64_t indexPhi{-1}, indexKstar{-1};
    for (const auto& daughter : mother.template daughters_as<TMcParticles>()) {
      ++nDaughters;
      switch (std::abs(daughter.pdgCode())) {
        case kPiPlus:
          ++nPions;
          break;
        case kKPlus:
          ++nKaons;
          break;
        case Pdg::kPhi:
          ++nPhis;
          indexPhi = daughter.globalIndex();
          break;
        case Pdg::kK0Star892:
          ++nKstars;
          indexKstar = daughter.globalIndex();
          break;
        default:
          break;
      }
    }
    if (nDaughters == NDaughtersTwoBody && nPions == 1 && nPhis == 1) {
      const auto phi = mcParticles.rawIteratorAt(indexPhi - mcParticles.offset());
      if (decaysToTwoBody(phi, mcParticles, kKPlus, kKPlus)) {
        return Channels::kDsPhiPi;
      }
    }
    if (nDaughters == NDaughtersTwoBody && nKaons == 1 && nKstars == 1) {
      const auto kstar = mcParticles.rawIteratorAt(indexKstar - mcParticles.offset());
      if (decaysToTwoBody(kstar, mcParticles, kKPlus, kPiPlus)) {
        return Channels::kDsKstarK;
      }
    }
    if (nDaughters == NDaughtersDsToKKPi && nKaons == NKaonsDsToKKPi && nPions == NPionsDsToKKPi) {
      return Channels::kDsDirect;
    }
    return Channels::kNone;
  }

  /// Whether a PDG code is a hadron containing a charm or beauty (anti)quark, quarkonia included.
  /// Decided from the quark digits of the PDG numbering scheme rather than from a list of species,
  /// so that the excited and rarer states a list misses (D*0, Ds*, Sigma_c, Xi_c0, Omega_c, B*,
  /// Sigma_b, charmonia, ...) are covered too. Diquarks (e.g. 4101) and nuclei are not hadrons here.
  static bool isHeavyFlavourHadron(int pdgCode)
  {
    const int pdgAbs = std::abs(pdgCode);
    if (pdgAbs < 100 || pdgAbs >= 1000000000) { // partons, leptons, bosons, generator codes; nuclei
      return false;
    }
    const int nq1 = (pdgAbs / 1000) % 10;
    const int nq2 = (pdgAbs / 100) % 10;
    const int nq3 = (pdgAbs / 10) % 10;
    if (pdgAbs < 10000 && nq1 != 0 && nq3 == 0) { // diquark
      return false;
    }
    const auto isHeavyQuark = [](int nq) { return nq == 4 || nq == 5; };
    return isHeavyQuark(nq1) || isHeavyQuark(nq2) || isHeavyQuark(nq3);
  }

  /// Whether the ancestry walk continues through a particle. It stops at quarks, gluons, diquarks
  /// and generator-internal objects (strings, clusters, codes 81-100): above a hadron those are the
  /// parton shower, whose history fans out to the beams (see maxAncestryDepth). A decay product of
  /// a heavy-flavour hadron reaches it through hadrons, leptons and photons only -- unless the
  /// generator decays it into partons first (e.g. PYTHIA B decays), which this walk does not follow.
  static bool continueAncestryWalk(int pdgCode)
  {
    const int pdgAbs = std::abs(pdgCode);
    if (pdgAbs <= 8 || pdgAbs == 21 || (pdgAbs >= 81 && pdgAbs <= 100)) {
      return false;
    }
    return !(pdgAbs >= 1000 && pdgAbs < 10000 && (pdgAbs / 10) % 10 == 0);
  }

  /// Whether the track's MC particle is a heavy-flavour hadron or descends from one, at any depth
  /// of the decay chain (feed-down through strange hadrons, conversions and interactions in the
  /// detector included), up to maxAncestryDepthHfFlag generations.
  /// Stored as an auxiliary column, so that the background veto is applied offline. Deliberately
  /// called only for the tracks that are actually written out.
  /// \param track is the track, joined with the MC labels
  /// \param mcParticles is the MC particle table
  /// \return true if a charm or beauty hadron is found among the ancestors; if the depth limit is
  ///         reached first, also true (and counted in hHfAncestryWalk), so that the veto errs on
  ///         the side of dropping a track
  template <typename TTrack, typename TMcParticles>
  bool hasHeavyFlavourAncestor(TTrack const& track, TMcParticles const& mcParticles)
  {
    if (!track.has_mcParticle()) {
      fillHfAncestryWalk(HfWalkNoMcParticle);
      return false;
    }
    const auto mcParticle = track.template mcParticle_as<TMcParticles>();
    if (isHeavyFlavourHadron(mcParticle.pdgCode())) {
      fillHfAncestryWalk(HfWalkFound);
      return true;
    }
    std::vector<int64_t> current{mcParticle.globalIndex()};
    std::vector<int64_t> next{};
    std::vector<int64_t> visited{mcParticle.globalIndex()};
    for (int depth = 0; depth < config.maxAncestryDepthHfFlag; ++depth) {
      next.clear();
      for (const auto idx : current) { // o2-linter: disable=const-ref-in-for-loop (int elements)
        const auto particle = mcParticles.rawIteratorAt(idx - mcParticles.offset());
        if (!particle.has_mothers()) {
          continue;
        }
        for (const auto motherId : particle.mothersIds()) { // o2-linter: disable=const-ref-in-for-loop (int elements)
          if (motherId < 0 || std::find(visited.begin(), visited.end(), motherId) != visited.end()) {
            continue;
          }
          visited.push_back(motherId);
          const auto mother = mcParticles.rawIteratorAt(motherId - mcParticles.offset());
          if (isHeavyFlavourHadron(mother.pdgCode())) {
            fillHfAncestryWalk(HfWalkFound);
            return true;
          }
          if (continueAncestryWalk(mother.pdgCode())) {
            next.push_back(motherId);
          }
        }
      }
      if (next.empty()) {
        fillHfAncestryWalk(HfWalkNotFound);
        return false;
      }
      current.swap(next);
    }
    fillHfAncestryWalk(HfWalkDepthLimit);
    return true;
  }

  void fillHfAncestryWalk(int outcome)
  {
    if (config.fillHistograms) {
      registry.fill(HIST("hHfAncestryWalk"), outcome);
    }
  }

  /// Whether the charm hadron descends from a beauty hadron.
  /// getCharmHadronOrigin walks the mother indices the same broken way as getMother, so the
  /// search is done here instead.
  template <typename TMcParticle, typename TMcParticles>
  bool isFromBeauty(TMcParticle const& mother, TMcParticles const& mcParticles)
  {
    for (const int pdgBeauty : {Pdg::kBPlus, Pdg::kB0, Pdg::kBS, Pdg::kLambdaB0}) { // o2-linter: disable=const-ref-in-for-loop (int elements)
      if (findAncestorWithPdg(mother, mcParticles, pdgBeauty, config.maxAncestryDepth) > -1) {
        return true;
      }
    }
    return false;
  }

  /// Charm hadron the track's MC particle is a *prong* of, or -1 if it is not one.
  ///
  /// A prong is a direct daughter of the charm hadron, or a direct daughter of a phi(1020) or a
  /// K*(892)0 that is itself a direct daughter of the charm hadron. Anything deeper is a decay or
  /// interaction product and must not be labelled as a prong, even though a plain depth-2
  /// ancestry search would happily find the charm hadron: a kaon from D+ -> K pi pi decaying in
  /// flight gives a pion whose grandmother is the D+, and calling that pion a D+ prong would teach
  /// the model to accept tracks that are not there.
  /// \param particle is the MC particle of the track
  /// \param mcParticles is the MC particle table
  /// \param pdgMother is the |PDG| code of the charm hadron
  /// \return the global index of the charm hadron, or -1
  template <typename TMcParticle, typename TMcParticles>
  int64_t findCharmMotherOfProng(TMcParticle const& particle, TMcParticles const& mcParticles, int& pdgMother)
  {
    if (!particle.has_mothers()) {
      return -1;
    }
    for (const auto motherId : particle.mothersIds()) { // o2-linter: disable=const-ref-in-for-loop (int elements)
      if (motherId < 0 || motherId >= particle.globalIndex()) {
        continue;
      }
      const auto mother = mcParticles.rawIteratorAt(motherId - mcParticles.offset());
      const int pdgMotherAbs = std::abs(mother.pdgCode());
      if (std::find(PdgMotherList.begin(), PdgMotherList.end(), pdgMotherAbs) != PdgMotherList.end()) {
        LOG(info) << "Found charm mother of prong: PDG " << pdgMotherAbs << ", global index " << mother.globalIndex();
        pdgMother = pdgMotherAbs;
        return motherId;
      }
      if (pdgMotherAbs != Pdg::kPhi && pdgMotherAbs != Pdg::kK0Star892) {
        continue;
      }
      if (!mother.has_mothers()) {
        continue;
      }
      for (const auto grandMotherId : mother.mothersIds()) { // o2-linter: disable=const-ref-in-for-loop (int elements)
        if (grandMotherId < 0 || grandMotherId >= mother.globalIndex()) {
          continue;
        }
        const auto grandMother = mcParticles.rawIteratorAt(grandMotherId - mcParticles.offset());
        if (std::find(PdgMotherList.begin(), PdgMotherList.end(), std::abs(grandMother.pdgCode())) != PdgMotherList.end()) {
          LOG(info) << "Found charm grandmother of prong: PDG " << std::abs(grandMother.pdgCode()) << ", global index " << grandMother.globalIndex();
          pdgMother = std::abs(grandMother.pdgCode());
          return grandMotherId;
        }
      }
    }
    return -1;
  }

  /// Sub-channel of a charm hadron decay, or Channels::kNone if it is not one of the wanted final
  /// states. The sub-channel is used to narrow the signal definition offline, but is not part of
  /// the training label itself.
  /// \param mother is the charm hadron MC particle
  /// \param mcParticles is the MC particle table
  /// \return the sub-channel, or Channels::kNone
  template <typename TMcParticle, typename TMcParticles>
  int8_t flagSubChannel(TMcParticle const& mother, TMcParticles const& mcParticles)
  {
    const int pdgMotherAbs = std::abs(mother.pdgCode());
    LOG(info) << "Flagging sub-channel for charm mother PDG " << pdgMotherAbs;
    if (pdgMotherAbs == Pdg::kD0) {
      return getDzeroToKPiChannel(mother, mcParticles);
    }
    if (pdgMotherAbs == Pdg::kDPlus) {
      return getDplusToKPiPiChannel(mother, mcParticles);
    }
    if (pdgMotherAbs == Pdg::kDS) {
      return getDsToKKPiChannel(mother, mcParticles);
    }
    if (pdgMotherAbs == Pdg::kDStar) {
      LOG(info) << "D* found, checking daughters ... ";
      return getDstarToD0PiChannel(mother, mcParticles);
    }
    return Channels::kNone;
  }

  /// Training label of a track for one charm species.
  ///
  /// Signal: the track is MC-matched to a pion or a kaon that is a prong (see
  ///   findCharmMotherOfProng) of a charm hadron of the requested species, produced in the
  ///   injected (non-background) event and decaying into the wanted final state. The resonant
  ///   sub-channels (D+ -> K*0 pi, Ds -> phi pi, Ds -> K*0 K) are included and the sub-channel is
  ///   recorded, so that the signal definition can still be narrowed offline.
  /// Background: the track comes from the background event and has no hadron of that species
  ///   anywhere in its ancestry -- background events do contain real charm, so the veto matters.
  /// Anything else -- in particular an underlying-event track of the injected event, and a track
  ///   from a charm hadron of the *other* species -- is not labelled: it is neither a daughter nor
  ///   a clean negative example.
  ///
  /// \param track is the track, joined with the MC labels
  /// \param mcParticles is the MC particle table
  /// \param pdgMother is the |PDG| code of the charm hadron
  /// \param keepIfBackground is the outcome of the background downsampling, decided by the caller
  ///        beforehand so that the ancestry veto -- the expensive part of the labelling -- is not
  ///        paid for the background tracks that are discarded anyway
  /// \param ptMother, isPrompt, channel are filled with the properties of the mother, if any
  /// \return the label, or -1 if the track must not enter this species' training sample
  template <typename TTrack, typename TMcParticles>
  int8_t getTrainingLabel(TTrack const& track, TMcParticles const& mcParticles,
                          int& pdgMother, bool keepIfBackground,
                          float& ptMother, int8_t& isPrompt, int8_t& channel)
  {
    ptMother = -1.f;
    isPrompt = -1;
    channel = Channels::kNone;

    if (!track.has_mcParticle()) {
      // a track without an MC particle cannot be a charm daughter; treat it as background
      return (config.keepFakeTracksInBkg && keepIfBackground) ? static_cast<int8_t>(HfTrackMlClass::Background) : -1;
    }
    const auto mcParticle = track.template mcParticle_as<TMcParticles>();

    const int pdgAbs = std::abs(mcParticle.pdgCode());
    const bool isPionOrKaon = (pdgAbs == kPiPlus || pdgAbs == kKPlus);

    if (!config.rejBkgCollisions || (config.rejBkgCollisions && !mcParticle.fromBackgroundEvent())) {
      if (!isPionOrKaon) {
        return -1;
      }
      LOG(info) << "\n\nFlagging track with PDG code " << mcParticle.pdgCode() << " ... ";
      for (const auto motherId : mcParticle.mothersIds()) {
        const auto m = mcParticles.rawIteratorAt(motherId - mcParticles.offset());
        LOG(info) << "kaon " << mcParticle.globalIndex() << " mother id " << motherId
                  << " PDG " << m.pdgCode() << " nDaughters " << m.daughtersIds().size();
      }
      const int64_t indexMother = findCharmMotherOfProng(mcParticle, mcParticles, pdgMother);
      if (indexMother < 0) {
        return -1;
      }
      const auto mother = mcParticles.rawIteratorAt(indexMother - mcParticles.offset());
      channel = flagSubChannel(mother, mcParticles);
      if (channel == Channels::kNone) {
        return -1;
      }
      ptMother = mother.pt();
      isPrompt = isFromBeauty(mother, mcParticles) ? 0 : 1;
      if (pdgAbs == kKPlus && channel == Channels::kDstarToD0Pi) {
        LOG(info) << "Track is a kaon from D* -> D0 pi, which is not a prong of the D* itself";
      }
      return (pdgAbs == kKPlus) ? static_cast<int8_t>(HfTrackMlClass::Kaon) : static_cast<int8_t>(HfTrackMlClass::Pion);
    }

    if (!keepIfBackground) { // downsampled away: skip the ancestry walk entirely
      return -1;
    }
    if (findAncestorWithPdg(mcParticle, mcParticles, pdgMother, config.maxAncestryDepth) > -1) {
      return -1;
    }
    return static_cast<int8_t>(HfTrackMlClass::Background);
  }

  /// Evaluate one model on a set of features.
  /// \param model is the model to evaluate
  /// \param features are the track features
  /// \param scores is filled with the three class scores (-1 if the track is outside the pT range)
  /// \param isSelPion, isSelKaon are filled with the threshold decisions
  /// \return true if the model was actually evaluated
  bool evaluateModel(TrackModel& model, HfTrackMlFeatures const& features,
                     std::array<float, 3>& scores, bool& isSelPion, bool& isSelKaon)
  {
    scores = {-1.f, -1.f, -1.f};
    isSelPion = false;
    isSelKaon = false;
    if (!model.enabled) {
      return false;
    }
    const int nModel = model.response.getModelBin(features.pt);
    if (nModel < 0) { // outside the pT range covered by the models
      return false;
    }
    auto inputFeatures = model.response.getInputFeatures(features);
    const auto output = model.response.getModelOutput(inputFeatures, nModel);

    // Collapse whatever the model emits onto (background, pion-role, kaon-role). For the 5-class
    // model the two origin scores of a role are summed: P(pion) = P(prompt pion) + P(non-prompt
    // pion), which is the probability the track is that role's daughter irrespective of where the
    // charm hadron came from. That keeps one threshold per role and makes the 3- and 5-class
    // working points mean the same thing.
    scores[0] = output[0];
    if (model.nClasses == NClassesWithOrigin) {
      scores[1] = output[1] + output[3];
      scores[2] = output[2] + output[4];
    } else {
      scores[1] = output[1];
      scores[2] = output[2];
    }

    isSelPion = scores[1] > model.thresholdPion[nModel];
    isSelKaon = scores[2] > model.thresholdKaon[nModel];
    return true;
  }

  /// One row of aod::HfSkimTrackSels for a track that was not scored.
  void fillEmptyTrackRow(bool isPositive)
  {
    rowSelectedTrack(-1.f, -1.f, -1.f, false, false,
                     -1.f, -1.f, -1.f, false, false,
                     -1.f, -1.f, -1.f, false, false,
                     -1.f, -1.f, -1.f, false, false,
                     isPositive, false);
  }

  /// Common loop over the collisions: event selection, then per-track work.
  /// With IsTraining and !IsMc the training table is filled on data: no labels (all -1), and every
  /// preselected association is downsampled with bkgSamplingFraction, to study the feature
  /// distributions the model sees on data.
  template <bool IsTraining, bool IsMc = true, typename TCollisions, typename TTracks, typename TMcParticles>
  void runTrackSelection(TCollisions const& collisions,
                         aod::TrackAssoc const& trackIndices,
                         TTracks const& /*tracks*/,
                         TMcParticles const& mcParticles,
                         aod::BcFullInfos const& /*bcs*/)
  {
    const bool doTiming = config.enableTiming;
    rowSelectedCollision.reserve(collisions.size());

    for (const auto& collision : collisions) {
      const auto tStartCollision = StopWatch::now();
      timing.fill(0.);
      registry.fill(HIST("hCounters"), CountTrackSelCollisions);

      const auto tStartEvSel = StopWatch::now();
      float centrality{-1.f};
      const auto rejectionMask = hfEvSel.getHfCollisionRejectionMask<true, CentralityEstimator::FT0C, aod::BcFullInfos>(collision, centrality, ccdb, registry);
      if (config.fillHistograms) {
        hfEvSel.fillHistograms(collision, rejectionMask, centrality);
      }
      rowSelectedCollision(rejectionMask, centrality);
      if (doTiming) {
        timing[TimeTrackSelEvSel] += elapsedSeconds(tStartEvSel);
      }

      const auto thisCollId = collision.globalIndex();
      const auto groupedTrackIndices = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);

      if (rejectionMask != 0) {
        // collision rejected: still fill one row per association to keep the table aligned
        if constexpr (!IsTraining) {
          for (auto iAssoc{0u}; iAssoc < groupedTrackIndices.size(); ++iAssoc) {
            fillEmptyTrackRow(false);
          }
        }
        continue;
      }

      // set the magnetic field for this run
      const auto bc = collision.template bc_as<aod::BcFullInfos>();
      initCCDB(bc, runNumber, ccdb, config.isRun2 ? config.ccdbPathGrp : config.ccdbPathGrpMag, lut, config.isRun2);

      [[maybe_unused]] int nPresel{0};
      [[maybe_unused]] int nSelPionDzero{0};
      [[maybe_unused]] int nSelKaonDzero{0};
      [[maybe_unused]] int nSelPionDplus{0};
      [[maybe_unused]] int nSelKaonDplus{0};
      [[maybe_unused]] int nSelPionDs{0};
      [[maybe_unused]] int nSelKaonDs{0};
      [[maybe_unused]] int nSelPionDstar{0};
      [[maybe_unused]] int nSelKaonDstar{0};
      [[maybe_unused]] int nSelAny{0};

      for (const auto& trackId : groupedTrackIndices) {
        const auto track = trackId.template track_as<TTracks>();
        registry.fill(HIST("hCounters"), CountAssociations);

        if (config.fillHistograms) {
          registry.fill(HIST("hPtAllTracks"), track.pt());
        }

        if (!passesCheapPreselection(track)) {
          if constexpr (!IsTraining) {
            // one row per association, so that the table stays aligned with aod::TrackAssoc
            fillEmptyTrackRow(track.sign() > 0);
          }
          continue;
        }
        registry.fill(HIST("hCounters"), CountCheapPresel);

        const auto tStartFeatures = StopWatch::now();
        HfTrackMlFeatures features{};
        fillFeatures(collision, track, centrality, features);
        if (doTiming) {
          timing[TimeTrackSelFeatures] += elapsedSeconds(tStartFeatures);
        }

        const bool preselected = isPreselected(track, features.pt, features.eta);
        if (preselected) {
          ++nPresel;
          registry.fill(HIST("hCounters"), CountPresel);
          if (config.fillHistograms) {
            registry.fill(HIST("hPtPreselected"), features.pt);
          }
        }

        if constexpr (IsTraining) {
          if (!preselected) {
            continue;
          }
          // The downsampling verdict is drawn first, from a deterministic hash of the track index,
          // so that the labelling can skip the charm ancestry veto for the background tracks that
          // would be thrown away right afterwards. Signal tracks are never downsampled.
          const bool keepIfBackground =
            deterministicUniform(static_cast<uint64_t>(track.globalIndex()), config.bkgSamplingSeed) < config.bkgSamplingFraction;

          if constexpr (!IsMc) {
            // data: no labels, a uniform random subsample of the preselected associations
            if (!keepIfBackground) {
              continue;
            }
            rowTrackTrain(features.pt, features.eta, features.dcaXY, features.dcaZ,
                          features.sigmaDcaXY, features.sigmaDcaZ, features.normDcaXY, features.normDcaZ,
                          features.signed1Pt, features.tgl, features.sign, features.isPvContributor,
                          features.itsNCls, features.itsNClsInnerBarrel, features.itsChi2NCl,
                          features.tpcNClsFound, features.tpcNClsCrossedRows, features.tpcCrossedRowsOverFindableCls,
                          features.tpcChi2NCl, features.tpcFractionSharedCls,
                          features.tpcNSigmaPi, features.tpcNSigmaKa, features.centrality,
                          track.globalIndex(),
                          -1, -1.f, -1, Channels::kNone, -1,
                          0, false);
            continue;
          } else {
            float pt{-1.f};
            int pdgMother{-1};
            int8_t isPrompt{-1}, channel{Channels::kNone};
            const int8_t label = getTrainingLabel(
                          track, mcParticles, pdgMother, keepIfBackground,
                          pt, isPrompt, channel);

            if (label < 0) { // usable for neither model
              continue;
            }
            // only now, on the tracks that are actually written out, pay for the ancestry scan
            const bool hasCharmBeautyAncestor = hasHeavyFlavourAncestor(track, mcParticles);
            if (config.fillHistograms) {
              if (channel == Channels::kDzeroDirect) {
                registry.fill(HIST("hLabelDzero"), label);
                if (channel >= 0) {
                  registry.fill(HIST("hChannelDzero"), channel);
                }
              }
              if (channel == Channels::kDplusDirect ||
                  channel == Channels::kDplusKstar) {
                registry.fill(HIST("hLabelDplus"), label);
                if (channel >= 0) {
                  registry.fill(HIST("hChannelDplus"), channel);
                }
              }
              if (channel == Channels::kDsDirect ||
                  channel == Channels::kDsKstarK ||
                  channel == Channels::kDsPhiPi) {
                registry.fill(HIST("hLabelDs"), label);
                if (channel >= 0) {
                  registry.fill(HIST("hChannelDs"), channel);
                }
              }
              if (channel == Channels::kDstarToD0Pi) {
                registry.fill(HIST("hLabelDstar"), label);
                if (channel >= 0) {
                  registry.fill(HIST("hChannelDstar"), channel);
                }
              }
            }
            rowTrackTrain(features.pt, features.eta, features.dcaXY, features.dcaZ,
                          features.sigmaDcaXY, features.sigmaDcaZ, features.normDcaXY, features.normDcaZ,
                          features.signed1Pt, features.tgl, features.sign, features.isPvContributor,
                          features.itsNCls, features.itsNClsInnerBarrel, features.itsChi2NCl,
                          features.tpcNClsFound, features.tpcNClsCrossedRows, features.tpcCrossedRowsOverFindableCls,
                          features.tpcChi2NCl, features.tpcFractionSharedCls,
                          features.tpcNSigmaPi, features.tpcNSigmaKa, features.centrality,
                          track.globalIndex(),
                          label, pt, isPrompt, channel, pdgMother,
                          track.mcMask(), hasCharmBeautyAncestor);
          }
        } else {
          if (!preselected) {
            fillEmptyTrackRow(track.sign() > 0);
            continue;
          }

          std::array<float, 3> scoresDzero{-1.f, -1.f, -1.f};
          std::array<float, 3> scoresDplus{-1.f, -1.f, -1.f};
          std::array<float, 3> scoresDs{-1.f, -1.f, -1.f};
          std::array<float, 3> scoresDstar{-1.f, -1.f, -1.f};
          bool isSelPionDzero{false}, isSelKaonDzero{false};
          bool isSelPionDplus{false}, isSelKaonDplus{false};
          bool isSelPionDs{false}, isSelKaonDs{false};
          bool isSelPionDstar{false}, isSelKaonDstar{false};

          const auto tStartDzero = StopWatch::now();
          if (evaluateModel(modelDzero, features, scoresDzero, isSelPionDzero, isSelKaonDzero)) {
            registry.fill(HIST("hCounters"), CountInferenceDzero);
          }
          if (doTiming) {
            timing[TimeTrackSelInferenceDzero] += elapsedSeconds(tStartDzero);
          }

          const auto tStartDplus = StopWatch::now();
          if (evaluateModel(modelDplus, features, scoresDplus, isSelPionDplus, isSelKaonDplus)) {
            registry.fill(HIST("hCounters"), CountInferenceDplus);
          }
          if (doTiming) {
            timing[TimeTrackSelInferenceDplus] += elapsedSeconds(tStartDplus);
          }

          const auto tStartDs = StopWatch::now();
          if (evaluateModel(modelDs, features, scoresDs, isSelPionDs, isSelKaonDs)) {
            registry.fill(HIST("hCounters"), CountInferenceDs);
          }
          if (doTiming) {
            timing[TimeTrackSelInferenceDs] += elapsedSeconds(tStartDs);
          }

          const auto tStartDstar = StopWatch::now();
          if (evaluateModel(modelDstar, features, scoresDstar, isSelPionDstar, isSelKaonDstar)) {
            registry.fill(HIST("hCounters"), CountInferenceDstar);
          }
          if (doTiming) {
            timing[TimeTrackSelInferenceDstar] += elapsedSeconds(tStartDstar);
          }

          uint16_t mask{0};
          if (isSelPionDzero) {
            SETBIT(mask, BitPionDzero);
            ++nSelPionDzero;
            registry.fill(HIST("hCounters"), CountSelPionDzero);
          }
          if (isSelKaonDzero) {
            SETBIT(mask, BitKaonDzero);
            ++nSelKaonDzero;
            registry.fill(HIST("hCounters"), CountSelKaonDzero);
          }
          if (isSelPionDplus) {
            SETBIT(mask, BitPionDplus);
            ++nSelPionDplus;
            registry.fill(HIST("hCounters"), CountSelPionDplus);
          }
          if (isSelKaonDplus) {
            SETBIT(mask, BitKaonDplus);
            ++nSelKaonDplus;
            registry.fill(HIST("hCounters"), CountSelKaonDplus);
          }
          if (isSelPionDs) {
            SETBIT(mask, BitPionDs);
            ++nSelPionDs;
            registry.fill(HIST("hCounters"), CountSelPionDs);
          }
          if (isSelKaonDs) {
            SETBIT(mask, BitKaonDs);
            ++nSelKaonDs;
            registry.fill(HIST("hCounters"), CountSelKaonDs);
          }
          if (isSelPionDstar) {
            SETBIT(mask, BitPionDstar);
            ++nSelPionDstar;
            registry.fill(HIST("hCounters"), CountSelPionDstar);
          }
          if (isSelKaonDstar) {
            SETBIT(mask, BitKaonDstar);
            ++nSelKaonDstar;
            registry.fill(HIST("hCounters"), CountSelKaonDstar);
          }
          const bool isSelAny = mask != 0;
          if (isSelAny) {
            ++nSelAny;
            registry.fill(HIST("hCounters"), CountSelAny);
          }

          if (config.fillHistograms) {
            registry.fill(HIST("hScoreBkgDzero"), scoresDzero[0]);
            registry.fill(HIST("hScorePionDzero"), scoresDzero[1]);
            registry.fill(HIST("hScoreKaonDzero"), scoresDzero[2]);
            registry.fill(HIST("hScoreBkgDplus"), scoresDplus[0]);
            registry.fill(HIST("hScorePionDplus"), scoresDplus[1]);
            registry.fill(HIST("hScoreKaonDplus"), scoresDplus[2]);
            registry.fill(HIST("hScoreBkgDs"), scoresDs[0]);
            registry.fill(HIST("hScorePionDs"), scoresDs[1]);
            registry.fill(HIST("hScoreKaonDs"), scoresDs[2]);
            registry.fill(HIST("hScoreBkgDstar"), scoresDstar[0]);
            registry.fill(HIST("hScorePionDstar"), scoresDstar[1]);
            registry.fill(HIST("hScoreKaonDstar"), scoresDstar[2]);
            if (isSelPionDzero) {
              registry.fill(HIST("hPtSelPionDzero"), features.pt);
            }
            if (isSelKaonDzero) {
              registry.fill(HIST("hPtSelKaonDzero"), features.pt);
            }
            if (isSelPionDplus) {
              registry.fill(HIST("hPtSelPionDplus"), features.pt);
            }
            if (isSelKaonDplus) {
              registry.fill(HIST("hPtSelKaonDplus"), features.pt);
            }
            if (isSelPionDs) {
              registry.fill(HIST("hPtSelPionDs"), features.pt);
            }
            if (isSelKaonDs) {
              registry.fill(HIST("hPtSelKaonDs"), features.pt);
            }
            if (isSelPionDstar) {
              registry.fill(HIST("hPtSelPionDstar"), features.pt);
            }
            if (isSelKaonDstar) {
              registry.fill(HIST("hPtSelKaonDstar"), features.pt);
            }
            if (isSelAny) {
              registry.fill(HIST("hRoleMask"), mask);
            }
          }

          const bool isPositive = track.sign() > 0;
          rowSelectedTrack(scoresDzero[0], scoresDzero[1], scoresDzero[2], isSelPionDzero, isSelKaonDzero,
                           scoresDplus[0], scoresDplus[1], scoresDplus[2], isSelPionDplus, isSelKaonDplus,
                           scoresDs[0], scoresDs[1], scoresDs[2], isSelPionDs, isSelKaonDs,
                           scoresDstar[0], scoresDstar[1], scoresDstar[2], isSelPionDstar, isSelKaonDstar,
                           isPositive, isSelAny);
          if (config.fillSelMaskTable && isSelAny) {
            uint16_t maskOut = mask;
            if (isPositive) {
              SETBIT(maskOut, BitPositive);
            }
            rowSelMask(static_cast<int>(thisCollId), maskOut);
          }
        }
      }

      if constexpr (!IsTraining) {
        if (config.fillHistograms) {
          registry.fill(HIST("hNTracksPreselPerColl"), nPresel);
          registry.fill(HIST("hNTracksSelPionDzeroPerColl"), nSelPionDzero);
          registry.fill(HIST("hNTracksSelKaonDzeroPerColl"), nSelKaonDzero);
          registry.fill(HIST("hNTracksSelPionDplusPerColl"), nSelPionDplus);
          registry.fill(HIST("hNTracksSelKaonDplusPerColl"), nSelKaonDplus);
          registry.fill(HIST("hNTracksSelPionDsPerColl"), nSelPionDs);
          registry.fill(HIST("hNTracksSelKaonDsPerColl"), nSelKaonDs);
          registry.fill(HIST("hNTracksSelPionDstarPerColl"), nSelPionDstar);
          registry.fill(HIST("hNTracksSelKaonDstarPerColl"), nSelKaonDstar);
          registry.fill(HIST("hNTracksSelAnyPerColl"), nSelAny);
        }
      }

      if (doTiming) {
        timing[TimeTrackSelTotal] = elapsedSeconds(tStartCollision);
        timing[TimeTrackSelRemainder] = timing[TimeTrackSelTotal] - timing[TimeTrackSelEvSel] - timing[TimeTrackSelFeatures] - timing[TimeTrackSelInferenceDplus] - timing[TimeTrackSelInferenceDs];
        for (int iStage = 0; iStage < NTimingTrackSel; ++iStage) {
          registry.fill(HIST("hTiming"), iStage, timing[iStage]);
        }
      }
    }
  }

  void processTraining(CollisionsWithCent const& collisions,
                       aod::TrackAssoc const& trackIndices,
                       TracksWithSelAndPidMc const& tracks,
                       aod::McParticles const& mcParticles,
                       aod::BcFullInfos const& bcs)
  {
    runTrackSelection<true>(collisions, trackIndices, tracks, mcParticles, bcs);
  }
  PROCESS_SWITCH(HfDplusTrackSelMl, processTraining, "Produce the track-level training sample (MC only)", false);

  void processTrainingData(CollisionsWithCent const& collisions,
                           aod::TrackAssoc const& trackIndices,
                           TracksWithSelAndPid const& tracks,
                           aod::BcFullInfos const& bcs)
  {
    runTrackSelection<true, false>(collisions, trackIndices, tracks, tracks, bcs);
  }
  PROCESS_SWITCH(HfDplusTrackSelMl, processTrainingData, "Write the unlabelled feature table on data, downsampled with bkgSamplingFraction", false);

  void processApplication(CollisionsWithCent const& collisions,
                          aod::TrackAssoc const& trackIndices,
                          TracksWithSelAndPid const& tracks,
                          aod::BcFullInfos const& bcs)
  {
    runTrackSelection<false>(collisions, trackIndices, tracks, tracks, bcs);
  }
  PROCESS_SWITCH(HfDplusTrackSelMl, processApplication, "Apply the track-level ML selection", true);
};

/// Charm 3-prong combinatorics over the ML-selected tracks.
struct HfDplusCombinatorics {
  Produces<aod::Hf3Prongs> rowTrackIndexProng3;
  Produces<aod::HfCandScoreScans> rowCandScoreScan;

  struct : ConfigurableGroup {
    Configurable<bool> fillHistograms{"fillHistograms", true, "fill QA histograms"};
    Configurable<bool> fillScoreScanTable{"fillScoreScanTable", false, "write one aod::HfCandScoreScans row per candidate, so the working-point scan can be redone offline"};
    Configurable<bool> enableTiming{"enableTiming", false, "account the CPU time of the individual stages into hTiming"};
    Configurable<double> ptTolerance{"ptTolerance", 0.1, "pT tolerance in GeV/c for the preselections before vertex reconstruction"};
    // D+ preselections, same layout and defaults as in the track-index skim creator
    Configurable<std::vector<double>> binsPtDplusToPiKPi{"binsPtDplusToPiKPi", std::vector<double>{hf_cuts_presel_3prong::vecBinsPt}, "pT bin limits for D+->piKpi pT-dependent cuts"};
    Configurable<LabeledArray<double>> cutsDplusToPiKPi{"cutsDplusToPiKPi", {hf_cuts_presel_3prong::Cuts[0], hf_cuts_presel_3prong::NBinsPt, hf_cuts_presel_3prong::NCutVars, hf_cuts_presel_3prong::labelsPt, hf_cuts_presel_3prong::labelsCutVar}, "D+->piKpi selections per pT bin"};
    // Ds preselections, same layout and defaults as in the track-index skim creator
    Configurable<std::vector<double>> binsPtDsToKKPi{"binsPtDsToKKPi", std::vector<double>{hf_cuts_presel_ds::vecBinsPt}, "pT bin limits for Ds->KKpi pT-dependent cuts"};
    Configurable<LabeledArray<double>> cutsDsToKKPi{"cutsDsToKKPi", {hf_cuts_presel_ds::Cuts[0], hf_cuts_presel_ds::NBinsPt, hf_cuts_presel_ds::NCutVars, hf_cuts_presel_ds::labelsPt, hf_cuts_presel_ds::labelsCutVar}, "Ds->KKpi selections per pT bin"};
    // two-track speedups, carried over from the skim creator so that the CPU comparison is fair
    Configurable<double> minTwoTrackDecayLengthFor3Prongs{"minTwoTrackDecayLengthFor3Prongs", 0., "minimum decay length computed with 2 tracks, to speed up the combinatorics"};
    Configurable<double> maxTwoTrackChi2PcaFor3Prongs{"maxTwoTrackChi2PcaFor3Prongs", 1.e10, "maximum chi2 PCA computed with 2 tracks, to speed up the combinatorics"};
    // vertexing
    Configurable<bool> propagateToPCA{"propagateToPCA", true, "create tracks version propagated to PCA"};
    Configurable<bool> useAbsDCA{"useAbsDCA", false, "minimise abs. distance rather than chi2"};
    Configurable<bool> useWeightedFinalPCA{"useWeightedFinalPCA", false, "recalculate vertex position using track covariances, effective only if useAbsDCA is true"};
    Configurable<double> maxR{"maxR", 200., "reject PCA's above this radius"};
    Configurable<double> maxDZIni{"maxDZIni", 4., "reject (if>0) PCA candidate if tracks DZ exceeds threshold"};
    Configurable<double> minParamChange{"minParamChange", 1.e-3, "stop iterations if largest change of any X is smaller than this"};
    Configurable<double> minRelChi2Change{"minRelChi2Change", 0.9, "stop iterations if chi2/chi2old > this"};
    // CCDB
    Configurable<bool> isRun2{"isRun2", false, "enable Run 2 or Run 3 GRP objects for the magnetic field"};
    Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
    Configurable<std::string> ccdbPathLut{"ccdbPathLut", "GLO/Param/MatLUT", "path for LUT parametrization"};
    Configurable<std::string> ccdbPathGrp{"ccdbPathGrp", "GLO/GRP/GRP", "path of the grp file (Run 2)"};
    Configurable<std::string> ccdbPathGrpMag{"ccdbPathGrpMag", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object (Run 3)"};
  } config;

  SliceCache cache;
  o2::vertexing::DCAFitterN<2> df2; // 2-prong fitter, only for the two-track speedups
  o2::vertexing::DCAFitterN<3> df3; // 3-prong vertex fitter
  Service<o2::ccdb::BasicCCDBManager> ccdb{};
  o2::base::MatLayerCylSet* lut{};
  o2::base::Propagator::MatCorrType noMatCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;
  int runNumber{};

  /// Prong ordering conventions of aod::Hf3Prongs, shared with the track-index skim creator:
  /// prong 1 is always the opposite-sign track, prongs 0 and 2 the same-sign pair.
  std::array<double, 3> massDplusToPiKPi{MassPiPlus, MassKPlus, MassPiPlus};
  std::array<std::array<double, 3>, 2> massDsToKKPi{std::array{MassKPlus, MassKPlus, MassPiPlus},  // hypothesis 0: K K pi
                                                    std::array{MassPiPlus, MassKPlus, MassKPlus}}; // hypothesis 1: pi K K

  /// One selected track, propagated to the collision under study once and reused by every triplet
  /// it takes part in. The original D+-only implementation re-propagated the innermost prong
  /// inside the innermost loop, i.e. once per triplet rather than once per track.
  struct ProngCandidate {
    o2::track::TrackParCov parCov;
    std::array<float, 3> pVec;
    int64_t trackIndex;
    uint8_t mask;
    // kept so that every built candidate can record the per-role minimum score it came from,
    // which is what makes the offline working-point scan possible from a single run
    float scorePionDplus;
    float scoreKaonDplus;
    float scorePionDs;
    float scoreKaonDs;
  };

  /// Per-collision working area, allocated once and cleared for each collision.
  struct ProngCache {
    std::vector<ProngCandidate> prongs;
    std::vector<int> pionDplus;
    std::vector<int> kaonDplus;
    std::vector<int> pionDs;
    std::vector<int> kaonDs;
    std::vector<int> dsSameSign;  ///< can sit in a same-sign slot of a Ds candidate: pi(Ds) or K(Ds)
    std::vector<int> sameSignAny; ///< can sit in a same-sign slot of any candidate: pi(D+), pi(Ds) or K(Ds)
    std::vector<int> oppSignAny;  ///< can sit in the opposite-sign slot: K(D+) or K(Ds)

    void clear()
    {
      prongs.clear();
      pionDplus.clear();
      kaonDplus.clear();
      pionDs.clear();
      kaonDs.clear();
      dsSameSign.clear();
      sameSignAny.clear();
      oppSignAny.clear();
    }
  };
  ProngCache cachePos;
  ProngCache cacheNeg;

  // per-collision accumulators
  std::array<double, NTimingCombinatorics> timing{};
  std::array<int64_t, NCounterCombinatorics> counters{};

  using SelectedCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::HfSkimSelColls>>;
  using TracksWithCov = aod::TracksWCovDcaExtra;
  using FilteredTrackAssocSel = soa::Filtered<soa::Join<aod::TrackAssoc, aod::HfSkimTrackSels>>;

  Filter filterSelectCollisions = (aod::hf_skim_sel_collision::whyRejectColl == static_cast<o2::hf_evsel::HfCollisionRejectionMask>(0));
  Filter filterSelectTrackIds = (aod::hf_skim_track_sel::isSelAny == true);

  // Registers the slice cache for the filtered association table. Without it the partitions below
  // cannot be sliced per collision and the device dies with
  //   (AOD/TRACKASSOC/0) /fIndexCollisions not found neither in sorted or unsorted cache
  Preslice<FilteredTrackAssocSel> trackIndicesPerCollision = aod::track_association::collisionId;

  // Partitions used by the legacy D+-only process function.
  Partition<FilteredTrackAssocSel> positivePions = aod::hf_skim_track_sel::isPositive == true && aod::hf_skim_track_sel::isSelPionDplus == true;
  Partition<FilteredTrackAssocSel> negativePions = aod::hf_skim_track_sel::isPositive == false && aod::hf_skim_track_sel::isSelPionDplus == true;
  Partition<FilteredTrackAssocSel> positiveKaons = aod::hf_skim_track_sel::isPositive == true && aod::hf_skim_track_sel::isSelKaonDplus == true;
  Partition<FilteredTrackAssocSel> negativeKaons = aod::hf_skim_track_sel::isPositive == false && aod::hf_skim_track_sel::isSelKaonDplus == true;

  HistogramRegistry registry{"registry"};

  void init(InitContext const&)
  {
    const std::array<int, 5> doProcess = {doprocessDplus, doprocessDplusCached, doprocessSeparateLoops, doprocessMergedLoop, doprocessDummy};
    if (std::accumulate(doProcess.begin(), doProcess.end(), 0) != 1) {
      LOGP(fatal, "One and only one process function of HfDplusCombinatorics must be enabled at a time!");
    }

    // CPU accounting and counters are booked unconditionally: an empty histogram is a clearer
    // statement than a missing one.
    registry.add("hTiming", "CPU time by stage;;CPU time (s)", {HistType::kTH1D, {{NTimingCombinatorics, -0.5, NTimingCombinatorics - 0.5}}});
    const char* timingLabels[NTimingCombinatorics] = {"total", "prong cache", "2-prong fits", "3-prong fits", "remainder"};
    for (int i = 0; i < NTimingCombinatorics; ++i) {
      registry.get<TH1>(HIST("hTiming"))->GetXaxis()->SetBinLabel(i + 1, timingLabels[i]);
    }
    registry.add("hCounters", "counters;;entries", {HistType::kTH1D, {{NCounterCombinatorics, -0.5, NCounterCombinatorics - 0.5}}});
    const char* counterLabels[NCounterCombinatorics] = {"collisions", "triplets", "triplets presel", "2-prong fits", "3-prong fits",
                                                        "rows", "rows D^{+}", "rows D_{s}", "rows both"};
    for (int i = 0; i < NCounterCombinatorics; ++i) {
      registry.get<TH1>(HIST("hCounters"))->GetXaxis()->SetBinLabel(i + 1, counterLabels[i]);
    }

    if (doprocessDummy) {
      return;
    }

    ccdb->setURL(config.ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(config.ccdbPathLut));

    df2.setPropagateToPCA(config.propagateToPCA);
    df2.setMaxR(config.maxR);
    df2.setMaxDZIni(config.maxDZIni);
    df2.setMinParamChange(config.minParamChange);
    df2.setMinRelChi2Change(config.minRelChi2Change);
    df2.setUseAbsDCA(config.useAbsDCA);
    df2.setWeightedFinalPCA(config.useWeightedFinalPCA);

    df3.setPropagateToPCA(config.propagateToPCA);
    df3.setMaxR(config.maxR);
    df3.setMaxDZIni(config.maxDZIni);
    df3.setMinParamChange(config.minParamChange);
    df3.setMinRelChi2Change(config.minRelChi2Change);
    df3.setUseAbsDCA(config.useAbsDCA);
    df3.setWeightedFinalPCA(config.useWeightedFinalPCA);

    if (config.fillHistograms) {
      registry.add("hMassDplusToPiKPi", "D^{+} candidates;inv. mass (#pi K #pi) (GeV/#it{c}^{2});entries", {HistType::kTH1D, {{500, 1.6, 2.1}}});
      registry.add("hMassDsToKKPi", "D_{s} candidates;inv. mass (K K #pi) (GeV/#it{c}^{2});entries", {HistType::kTH1D, {{500, 1.6, 2.3}}});
      registry.add("hMassKK", "D_{s} candidates;inv. mass (K K) (GeV/#it{c}^{2});entries", {HistType::kTH1D, {{200, 0.98, 1.18}}});
      registry.add("hPtCand", "3-prong candidates;#it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1D, {{200, 0., 20.}}});
      registry.add("hNCandPerColl", "candidates per collision;n. candidates;entries", {HistType::kTH1D, {{500, -0.5, 499.5}}});
      registry.add("hVtxX", "secondary vertex;x (cm);entries", {HistType::kTH1D, {{200, -2., 2.}}});
    }
  }

  /// Preselection of a D+ candidate before vertex reconstruction: pT bin, invariant mass and
  /// minimum prong pT. Mirrors applyPreselection3Prong of the track-index skim creator, reduced
  /// to the D+ hypothesis.
  /// \param pVecTrack0, pVecTrack1, pVecTrack2 are the prong momenta (π, K, π ordering)
  /// \return true if the candidate survives
  template <typename TVec>
  bool applyPreselectionDplus(TVec const& pVecTrack0, TVec const& pVecTrack1, TVec const& pVecTrack2)
  {
    const auto pt = RecoDecay::pt(pVecTrack0, pVecTrack1, pVecTrack2) + config.ptTolerance;
    const auto binPt = findBin(config.binsPtDplusToPiKPi, pt);
    if (binPt == -1) {
      return false;
    }

    const double minMass = config.cutsDplusToPiKPi->get(binPt, 0u);
    const double maxMass = config.cutsDplusToPiKPi->get(binPt, 1u);
    if (minMass >= 0. && maxMass > 0.) {
      const std::array arrMom{pVecTrack0, pVecTrack1, pVecTrack2};
      const double mass2 = RecoDecay::m2(arrMom, massDplusToPiKPi);
      if (mass2 < minMass * minMass || mass2 >= maxMass * maxMass) {
        return false;
      }
    }

    const auto ptProngMin = config.cutsDplusToPiKPi->get(binPt, 4u);
    const auto pt2ProngMin = ptProngMin * ptProngMin;
    if (RecoDecay::pt2(pVecTrack0) < pt2ProngMin || RecoDecay::pt2(pVecTrack1) < pt2ProngMin || RecoDecay::pt2(pVecTrack2) < pt2ProngMin) {
      return false;
    }

    return true;
  }

  /// Preselection of a Ds candidate before vertex reconstruction, for the mass hypotheses still
  /// allowed by the track selection. Mirrors applyPreselection3Prong plus applyPreselectionPhiDecay
  /// of the track-index skim creator: the two hypotheses are K K π (bit 0) and π K K (bit 1), the
  /// opposite-sign prong being the middle one in both, and each surviving hypothesis must also
  /// have its kaon pair compatible with the φ(1020).
  /// \param pVecTrack0, pVecTrack1, pVecTrack2 are the prong momenta
  /// \param hypoMask is the incoming hypothesis mask, narrowed in place
  /// \return true if at least one hypothesis survives
  template <typename TVec>
  bool applyPreselectionDs(TVec const& pVecTrack0, TVec const& pVecTrack1, TVec const& pVecTrack2, uint8_t& hypoMask)
  {
    const auto pt = RecoDecay::pt(pVecTrack0, pVecTrack1, pVecTrack2) + config.ptTolerance;
    const auto binPt = findBin(config.binsPtDsToKKPi, pt);
    if (binPt == -1) {
      hypoMask = 0;
      return false;
    }

    const auto ptProngMin = config.cutsDsToKKPi->get(binPt, 4u);
    const auto pt2ProngMin = ptProngMin * ptProngMin;
    if (RecoDecay::pt2(pVecTrack0) < pt2ProngMin || RecoDecay::pt2(pVecTrack1) < pt2ProngMin || RecoDecay::pt2(pVecTrack2) < pt2ProngMin) {
      hypoMask = 0;
      return false;
    }

    const double minMass = config.cutsDsToKKPi->get(binPt, 0u);
    const double maxMass = config.cutsDsToKKPi->get(binPt, 1u);
    const std::array arrMom{pVecTrack0, pVecTrack1, pVecTrack2};
    if (minMass >= 0. && maxMass > 0.) {
      for (uint8_t iHypo = 0; iHypo < 2; ++iHypo) {
        if (!TESTBIT(hypoMask, iHypo)) {
          continue;
        }
        const double mass2 = RecoDecay::m2(arrMom, massDsToKKPi[iHypo]);
        if (mass2 < minMass * minMass || mass2 >= maxMass * maxMass) {
          CLRBIT(hypoMask, iHypo);
        }
      }
      if (hypoMask == 0) {
        return false;
      }
    }

    // φ(1020) window on the kaon pair: prongs (0, 1) for the K K π hypothesis, (1, 2) for π K K
    const double deltaMassKKMax = config.cutsDsToKKPi->get(binPt, 5u);
    if (deltaMassKKMax > 0.) {
      if (TESTBIT(hypoMask, 0)) {
        const double massKK = RecoDecay::m(std::array{pVecTrack0, pVecTrack1}, std::array{MassKPlus, MassKPlus});
        if (std::abs(massKK - MassPhi) > deltaMassKKMax) {
          CLRBIT(hypoMask, 0);
        }
      }
      if (TESTBIT(hypoMask, 1)) {
        const double massKK = RecoDecay::m(std::array{pVecTrack1, pVecTrack2}, std::array{MassKPlus, MassKPlus});
        if (std::abs(massKK - MassPhi) > deltaMassKKMax) {
          CLRBIT(hypoMask, 1);
        }
      }
    }

    return hypoMask != 0;
  }

  /// Selection after vertex reconstruction: cos(pointing angle) and decay length.
  /// Mirrors applySelection3Prong of the track-index skim creator.
  template <typename TBins, typename TCuts, typename TVec, typename TVtx, typename TPv>
  bool applySelectionAfterVertexing(TBins const& binsPt, TCuts const& cuts, TVec const& pVecCand, TVtx const& secVtx, TPv const& primVtx)
  {
    const auto binPt = findBin(binsPt, RecoDecay::pt(pVecCand));
    if (binPt == -1) {
      return false;
    }
    if (RecoDecay::cpa(primVtx, secVtx, pVecCand) < cuts->get(binPt, 2u)) {
      return false;
    }
    if (RecoDecay::distance(primVtx, secVtx) < cuts->get(binPt, 3u)) {
      return false;
    }
    return true;
  }

  /// Check whether the vertex built with the first two prongs is compatible with a 3-prong decay.
  /// Used only if the corresponding speedup configurables are set.
  template <typename TVtx, typename TPv, typename TFitter>
  bool isTwoTrackVertexSelected(TVtx const& secVtx, TPv const& primVtx, TFitter const& fitter)
  {
    if (config.maxTwoTrackChi2PcaFor3Prongs < MaxTwoTrackChi2PcaDisabled && fitter.getChi2AtPCACandidate() > config.maxTwoTrackChi2PcaFor3Prongs) {
      return false;
    }
    if (config.minTwoTrackDecayLengthFor3Prongs > 0.f && RecoDecay::distance(primVtx, secVtx) < config.minTwoTrackDecayLengthFor3Prongs) {
      return false;
    }
    return true;
  }

  /// Whether the two-track speedups are configured at all.
  bool twoTrackSpeedupsEnabled()
  {
    return config.minTwoTrackDecayLengthFor3Prongs > 0.f || config.maxTwoTrackChi2PcaFor3Prongs < MaxTwoTrackChi2PcaDisabled;
  }

  /// Fill the per-collision prong cache: propagate every selected association to this collision
  /// once, and build the per-role index lists the loops iterate over.
  /// \param collision is the collision under study
  /// \param groupedAssocs are its selected (track, collision) associations
  /// \param tracks is the track table
  template <typename TCollision, typename TAssocs, typename TTracks>
  void buildProngCache(TCollision const& collision, TAssocs const& groupedAssocs, TTracks const& /*tracks*/)
  {
    cachePos.clear();
    cacheNeg.clear();
    const auto thisCollId = collision.globalIndex();

    for (const auto& assocId : groupedAssocs) {
      const auto track = assocId.template track_as<TTracks>();

      uint8_t mask{0};
      if (assocId.isSelPionDplus()) {
        SETBIT(mask, BitPionDplus);
      }
      if (assocId.isSelKaonDplus()) {
        SETBIT(mask, BitKaonDplus);
      }
      if (assocId.isSelPionDs()) {
        SETBIT(mask, BitPionDs);
      }
      if (assocId.isSelKaonDs()) {
        SETBIT(mask, BitKaonDs);
      }
      if (mask == 0) {
        continue;
      }

      ProngCandidate prong{getTrackParCov(track), track.pVector(), track.globalIndex(), mask,
                           assocId.scorePionDplus(), assocId.scoreKaonDplus(),
                           assocId.scorePionDs(), assocId.scoreKaonDs()};
      if (thisCollId != track.collisionId()) {
        // not the "default" collision of this track: re-propagate it, once
        std::array dcaInfo{-999.f, -999.f};
        o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, prong.parCov, 2.f, noMatCorr, &dcaInfo);
        getPxPyPz(prong.parCov, prong.pVec);
      }

      auto& cacheThis = assocId.isPositive() ? cachePos : cacheNeg;
      const int index = static_cast<int>(cacheThis.prongs.size());
      cacheThis.prongs.push_back(prong);
      if (TESTBIT(mask, BitPionDplus)) {
        cacheThis.pionDplus.push_back(index);
      }
      if (TESTBIT(mask, BitKaonDplus)) {
        cacheThis.kaonDplus.push_back(index);
      }
      if (TESTBIT(mask, BitPionDs)) {
        cacheThis.pionDs.push_back(index);
      }
      if (TESTBIT(mask, BitKaonDs)) {
        cacheThis.kaonDs.push_back(index);
      }
      if (TESTBIT(mask, BitPionDs) || TESTBIT(mask, BitKaonDs)) {
        cacheThis.dsSameSign.push_back(index);
      }
      if (TESTBIT(mask, BitPionDplus) || TESTBIT(mask, BitPionDs) || TESTBIT(mask, BitKaonDs)) {
        cacheThis.sameSignAny.push_back(index);
      }
      if (TESTBIT(mask, BitKaonDplus) || TESTBIT(mask, BitKaonDs)) {
        cacheThis.oppSignAny.push_back(index);
      }
    }
  }

  /// Run the 2-prong speedup fit on a (same-sign, opposite-sign) pair.
  /// \return true if the pair is worth continuing with
  template <typename TPv>
  bool passesTwoTrackSpeedup(ProngCandidate const& prongSameSign, ProngCandidate const& prongOppSign, TPv const& pvCoord, bool doTiming)
  {
    const auto tStart = StopWatch::now();
    int nVtx{0};
    ++counters[CountFit2Prong];
    try {
      nVtx = df2.process(prongSameSign.parCov, prongOppSign.parCov);
    } catch (...) {
      nVtx = 0;
    }
    const bool ok = (nVtx != 0) && isTwoTrackVertexSelected(df2.getPCACandidate(), pvCoord, df2);
    if (doTiming) {
      timing[TimeCombFit2Prong] += elapsedSeconds(tStart);
    }
    return ok;
  }

  /// Fit the 3-prong vertex and, if it converges, hand back the refitted prong momenta.
  /// \return true if the fit converged
  bool fitThreeProngVertex(ProngCandidate const& p0, ProngCandidate const& p1, ProngCandidate const& p2,
                           std::array<float, 3>& pvec0, std::array<float, 3>& pvec1, std::array<float, 3>& pvec2,
                           bool doTiming)
  {
    const auto tStart = StopWatch::now();
    int nVtx{0};
    ++counters[CountFit3Prong];
    try {
      nVtx = df3.process(p0.parCov, p1.parCov, p2.parCov);
    } catch (...) {
      nVtx = 0;
    }
    if (nVtx != 0) {
      df3.getTrack(0).getPxPyPzGlo(pvec0);
      df3.getTrack(1).getPxPyPzGlo(pvec1);
      df3.getTrack(2).getPxPyPzGlo(pvec2);
    }
    if (doTiming) {
      timing[TimeCombFit3Prong] += elapsedSeconds(tStart);
    }
    return nVtx != 0;
  }

  /// The triple loop itself, shared by every strategy so that they differ only in the track lists
  /// they iterate over and in the hypotheses they are allowed to test.
  ///
  /// The prong ordering handed to the fitter and written to aod::Hf3Prongs is
  /// (same-sign 1, opposite-sign, same-sign 2), the convention of the track-index skim creator.
  ///
  /// \param collision is the collision under study
  /// \param sameSign, oppSign are the per-charge prong caches, in the roles this charge
  ///        combination gives them
  /// \param sameSignList, oppSignList are the index lists to loop over
  /// \param pairIsUseful(maskOpp, maskSame1) rejects a (same-sign, opposite-sign) pair before the
  ///        2-prong speedup fit is paid for it
  /// \param hypothesisOf(maskOpp, mask1, mask2, dplusOk, dsHypoMask) fills the hypotheses this
  ///        strategy is allowed to test on one triplet
  /// \param nCand counts the rows written
  template <typename TCollision, typename TPairFilter, typename THypothesis>
  void runTripleLoop(TCollision const& collision,
                     ProngCache const& sameSign, ProngCache const& oppSign,
                     std::vector<int> const& sameSignList, std::vector<int> const& oppSignList,
                     TPairFilter pairIsUseful, THypothesis hypothesisOf,
                     bool doTiming, int& nCand)
  {
    if (sameSignList.size() < 2 || oppSignList.empty()) {
      return;
    }
    const auto thisCollId = collision.globalIndex();
    const std::array pvCoord{collision.posX(), collision.posY(), collision.posZ()};
    const bool useTwoTrackSpeedups = twoTrackSpeedupsEnabled();

    for (const int iOpp : oppSignList) { // o2-linter: disable=const-ref-in-for-loop (int elements)
      const auto& prongOpp = oppSign.prongs[iOpp];

      for (std::size_t i1 = 0; i1 + 1 < sameSignList.size(); ++i1) {
        const auto& prong1 = sameSign.prongs[sameSignList[i1]];
        if (!pairIsUseful(prongOpp.mask, prong1.mask)) {
          continue;
        }
        if (useTwoTrackSpeedups && !passesTwoTrackSpeedup(prong1, prongOpp, pvCoord, doTiming)) {
          continue;
        }

        for (std::size_t i2 = i1 + 1; i2 < sameSignList.size(); ++i2) {
          const auto& prong2 = sameSign.prongs[sameSignList[i2]];

          bool dplusOk{false};
          uint8_t dsHypo{0};
          hypothesisOf(prongOpp.mask, prong1.mask, prong2.mask, dplusOk, dsHypo);
          if (!dplusOk && dsHypo == 0) {
            continue;
          }
          ++counters[CountTriplets];

          if (dplusOk && !applyPreselectionDplus(prong1.pVec, prongOpp.pVec, prong2.pVec)) {
            dplusOk = false;
          }
          if (dsHypo != 0) {
            applyPreselectionDs(prong1.pVec, prongOpp.pVec, prong2.pVec, dsHypo);
          }
          if (!dplusOk && dsHypo == 0) {
            continue;
          }
          ++counters[CountTripletsPresel];

          std::array<float, 3> pvec0{};
          std::array<float, 3> pvec1{};
          std::array<float, 3> pvec2{};
          if (!fitThreeProngVertex(prong1, prongOpp, prong2, pvec0, pvec1, pvec2, doTiming)) {
            continue;
          }

          const auto& secondaryVertex = df3.getPCACandidate();
          const auto pVecCand = RecoDecay::pVec(pvec0, pvec1, pvec2);

          uint8_t hfFlag{0};
          if (dplusOk && applySelectionAfterVertexing(config.binsPtDplusToPiKPi, config.cutsDplusToPiKPi, pVecCand, secondaryVertex, pvCoord)) {
            SETBIT(hfFlag, hf_cand_3prong::DecayType::DplusToPiKPi);
          }
          if (dsHypo != 0 && applySelectionAfterVertexing(config.binsPtDsToKKPi, config.cutsDsToKKPi, pVecCand, secondaryVertex, pvCoord)) {
            SETBIT(hfFlag, hf_cand_3prong::DecayType::DsToKKPi);
          }
          if (hfFlag == 0) {
            continue;
          }

          rowTrackIndexProng3(thisCollId, prong1.trackIndex, prongOpp.trackIndex, prong2.trackIndex, hfFlag);
          ++nCand;
          ++counters[CountRows];
          const bool isDplus = TESTBIT(hfFlag, hf_cand_3prong::DecayType::DplusToPiKPi);
          const bool isDs = TESTBIT(hfFlag, hf_cand_3prong::DecayType::DsToKKPi);

          if (config.fillScoreScanTable) {
            // Prong order is (sameSign1, oppSign, sameSign2). For the D+ the opposite-sign prong
            // is the kaon and both same-sign prongs are pions. For the Ds the opposite-sign prong
            // is always a kaon, and the two hypotheses differ in which same-sign prong is the
            // pion -- dsHypo carries only the hypotheses that survived the mass and phi windows.
            const bool hypoKKPi = isDs && TESTBIT(dsHypo, 0); // prong1 = K, prong2 = pi
            const bool hypoPiKK = isDs && TESTBIT(dsHypo, 1); // prong1 = pi, prong2 = K
            rowCandScoreScan(RecoDecay::pt(pVecCand),
                             isDplus ? std::min(prong1.scorePionDplus, prong2.scorePionDplus) : kScoreScanNone,
                             isDplus ? prongOpp.scoreKaonDplus : kScoreScanNone,
                             hypoKKPi ? std::min(prong1.scoreKaonDs, prongOpp.scoreKaonDs) : kScoreScanNone,
                             hypoKKPi ? prong2.scorePionDs : kScoreScanNone,
                             hypoPiKK ? std::min(prong2.scoreKaonDs, prongOpp.scoreKaonDs) : kScoreScanNone,
                             hypoPiKK ? prong1.scorePionDs : kScoreScanNone);
          }
          if (isDplus) {
            ++counters[CountRowsDplus];
          }
          if (isDs) {
            ++counters[CountRowsDs];
          }
          if (isDplus && isDs) {
            ++counters[CountRowsBoth];
          }

          if (config.fillHistograms) {
            const std::array arrMom{pvec0, pvec1, pvec2};
            if (isDplus) {
              registry.fill(HIST("hMassDplusToPiKPi"), RecoDecay::m(arrMom, massDplusToPiKPi));
            }
            if (isDs) {
              for (uint8_t iHypo = 0; iHypo < 2; ++iHypo) {
                if (!TESTBIT(dsHypo, iHypo)) {
                  continue;
                }
                registry.fill(HIST("hMassDsToKKPi"), RecoDecay::m(arrMom, massDsToKKPi[iHypo]));
                const auto& pVecKaonA = (iHypo == 0) ? pvec0 : pvec1;
                const auto& pVecKaonB = (iHypo == 0) ? pvec1 : pvec2;
                registry.fill(HIST("hMassKK"), RecoDecay::m(std::array{pVecKaonA, pVecKaonB}, std::array{MassKPlus, MassKPlus}));
              }
            }
            registry.fill(HIST("hPtCand"), RecoDecay::pt(pVecCand));
            registry.fill(HIST("hVtxX"), secondaryVertex[0]);
          }
        }
      }
    }
  }

  /// Hypothesis evaluators. The prong ordering is (same-sign 1, opposite-sign, same-sign 2), so
  /// D+ -> π K π needs both same-sign prongs pion-like, and Ds -> K K π needs the opposite-sign
  /// prong kaon-like plus one kaon and one pion among the same-sign pair -- two orderings, which
  /// become the two mass hypotheses of aod::Hf3Prongs.
  static void hypothesisDplusOnly(uint8_t maskOpp, uint8_t mask1, uint8_t mask2, bool& dplusOk, uint8_t& dsHypo)
  {
    dsHypo = 0;
    dplusOk = TESTBIT(maskOpp, BitKaonDplus) && TESTBIT(mask1, BitPionDplus) && TESTBIT(mask2, BitPionDplus);
  }

  static void hypothesisDsOnly(uint8_t maskOpp, uint8_t mask1, uint8_t mask2, bool& dplusOk, uint8_t& dsHypo)
  {
    dplusOk = false;
    dsHypo = 0;
    if (!TESTBIT(maskOpp, BitKaonDs)) {
      return;
    }
    if (TESTBIT(mask1, BitKaonDs) && TESTBIT(mask2, BitPionDs)) {
      SETBIT(dsHypo, 0);
    }
    if (TESTBIT(mask1, BitPionDs) && TESTBIT(mask2, BitKaonDs)) {
      SETBIT(dsHypo, 1);
    }
  }

  static void hypothesisBoth(uint8_t maskOpp, uint8_t mask1, uint8_t mask2, bool& dplusOk, uint8_t& dsHypo)
  {
    hypothesisDplusOnly(maskOpp, mask1, mask2, dplusOk, dsHypo);
    bool dummy{false};
    hypothesisDsOnly(maskOpp, mask1, mask2, dummy, dsHypo);
  }

  static bool pairUsefulDplus(uint8_t maskOpp, uint8_t mask1)
  {
    return TESTBIT(maskOpp, BitKaonDplus) && TESTBIT(mask1, BitPionDplus);
  }

  static bool pairUsefulDs(uint8_t maskOpp, uint8_t mask1)
  {
    return TESTBIT(maskOpp, BitKaonDs) && (TESTBIT(mask1, BitKaonDs) || TESTBIT(mask1, BitPionDs));
  }

  static bool pairUsefulBoth(uint8_t maskOpp, uint8_t mask1)
  {
    return pairUsefulDplus(maskOpp, mask1) || pairUsefulDs(maskOpp, mask1);
  }

  /// Per-collision preamble shared by the two new strategies: magnetic field, prong cache, timers.
  template <typename TCollision, typename TAssocs, typename TTracks>
  void prepareCollision(TCollision const& collision, TAssocs const& trackIndices, TTracks const& tracks, bool doTiming)
  {
    const auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc, runNumber, ccdb, config.isRun2 ? config.ccdbPathGrp : config.ccdbPathGrpMag, lut, config.isRun2);
    const auto bz = o2::base::Propagator::Instance()->getNominalBz();
    df2.setBz(bz);
    df3.setBz(bz);

    const auto tStart = StopWatch::now();
    const auto groupedAssocs = trackIndices.sliceBy(trackIndicesPerCollision, collision.globalIndex());
    buildProngCache(collision, groupedAssocs, tracks);
    if (doTiming) {
      timing[TimeCombCache] += elapsedSeconds(tStart);
    }
  }

  /// Book-keeping at the end of a collision: fill the timing and counter histograms.
  void finaliseCollision(StopWatch::time_point const& tStartCollision, bool doTiming, int nCand)
  {
    if (doTiming) {
      timing[TimeCombTotal] = elapsedSeconds(tStartCollision);
      timing[TimeCombRemainder] = timing[TimeCombTotal] - timing[TimeCombCache] - timing[TimeCombFit2Prong] - timing[TimeCombFit3Prong];
      for (int iStage = 0; iStage < NTimingCombinatorics; ++iStage) {
        registry.fill(HIST("hTiming"), iStage, timing[iStage]);
      }
    }
    for (int iCounter = 0; iCounter < NCounterCombinatorics; ++iCounter) {
      if (counters[iCounter] != 0) {
        registry.fill(HIST("hCounters"), iCounter, static_cast<double>(counters[iCounter]));
      }
    }
    if (config.fillHistograms) {
      registry.fill(HIST("hNCandPerColl"), nCand);
    }
  }

  /// Legacy D+-only loop, kept verbatim except for the renamed selection columns so that the
  /// earlier efficiency/rate/CPU comparison against the track-index skim creator stays
  /// reproducible. It re-propagates the innermost prong inside the innermost loop, which is
  /// exactly what the prong cache of the newer strategies removes.
  template <typename TCollision, typename TTracks, typename TSameSign, typename TOppSign>
  int buildCandidatesLegacy(TCollision const& collision,
                            TTracks const& /*tracks*/,
                            TSameSign const& sameSignTracks,
                            TOppSign const& oppSignTracks)
  {
    const auto thisCollId = collision.globalIndex();
    const std::array pvCoord{collision.posX(), collision.posY(), collision.posZ()};
    int nCand{0};

    for (const auto& oppSignId : oppSignTracks) {
      const auto trackKa = oppSignId.template track_as<TTracks>();
      auto trackParVarKa = getTrackParCov(trackKa);
      std::array pVecTrackKa{trackKa.pVector()};
      if (thisCollId != trackKa.collisionId()) {
        std::array dcaInfo{-999.f, -999.f};
        o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParVarKa, 2.f, noMatCorr, &dcaInfo);
        getPxPyPz(trackParVarKa, pVecTrackKa);
      }

      for (auto sameSignId1 = sameSignTracks.begin(); sameSignId1 != sameSignTracks.end(); ++sameSignId1) {
        const auto trackPi1 = sameSignId1.template track_as<TTracks>();
        auto trackParVarPi1 = getTrackParCov(trackPi1);
        std::array pVecTrackPi1{trackPi1.pVector()};
        if (thisCollId != trackPi1.collisionId()) {
          std::array dcaInfo{-999.f, -999.f};
          o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParVarPi1, 2.f, noMatCorr, &dcaInfo);
          getPxPyPz(trackParVarPi1, pVecTrackPi1);
        }

        // optional two-track speedups, evaluated once per (π, K) pair
        if (twoTrackSpeedupsEnabled()) {
          int nVtxFrom2ProngFitter{0};
          ++counters[CountFit2Prong];
          try {
            nVtxFrom2ProngFitter = df2.process(trackParVarPi1, trackParVarKa);
          } catch (...) {
            continue;
          }
          if (nVtxFrom2ProngFitter == 0) {
            continue;
          }
          if (!isTwoTrackVertexSelected(df2.getPCACandidate(), pvCoord, df2)) {
            continue;
          }
        }

        for (auto sameSignId2 = sameSignId1 + 1; sameSignId2 != sameSignTracks.end(); ++sameSignId2) {
          const auto trackPi2 = sameSignId2.template track_as<TTracks>();
          auto trackParVarPi2 = getTrackParCov(trackPi2);
          std::array pVecTrackPi2{trackPi2.pVector()};
          if (thisCollId != trackPi2.collisionId()) {
            std::array dcaInfo{-999.f, -999.f};
            o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParVarPi2, 2.f, noMatCorr, &dcaInfo);
            getPxPyPz(trackParVarPi2, pVecTrackPi2);
          }

          ++counters[CountTriplets];
          if (!applyPreselectionDplus(pVecTrackPi1, pVecTrackKa, pVecTrackPi2)) {
            continue;
          }
          ++counters[CountTripletsPresel];

          int nVtxFrom3ProngFitter{0};
          ++counters[CountFit3Prong];
          try {
            nVtxFrom3ProngFitter = df3.process(trackParVarPi1, trackParVarKa, trackParVarPi2);
          } catch (...) {
            continue;
          }
          if (nVtxFrom3ProngFitter == 0) {
            continue;
          }

          const auto& secondaryVertex = df3.getPCACandidate();
          std::array<float, 3> pvec0{};
          std::array<float, 3> pvec1{};
          std::array<float, 3> pvec2{};
          df3.getTrack(0).getPxPyPzGlo(pvec0);
          df3.getTrack(1).getPxPyPzGlo(pvec1);
          df3.getTrack(2).getPxPyPzGlo(pvec2);
          const auto pVecCand = RecoDecay::pVec(pvec0, pvec1, pvec2);

          if (!applySelectionAfterVertexing(config.binsPtDplusToPiKPi, config.cutsDplusToPiKPi, pVecCand, secondaryVertex, pvCoord)) {
            continue;
          }

          uint8_t hfFlag{0};
          SETBIT(hfFlag, hf_cand_3prong::DecayType::DplusToPiKPi);
          rowTrackIndexProng3(thisCollId, trackPi1.globalIndex(), trackKa.globalIndex(), trackPi2.globalIndex(), hfFlag);
          if (config.fillScoreScanTable) {
            // the legacy path re-propagates per triplet and keeps no prong cache, so it has no
            // scores to report; the row still has to exist to keep the two tables in lockstep
            rowCandScoreScan(RecoDecay::pt(pVecCand), kScoreScanNone, kScoreScanNone,
                             kScoreScanNone, kScoreScanNone, kScoreScanNone, kScoreScanNone);
          }
          ++nCand;
          ++counters[CountRows];
          ++counters[CountRowsDplus];

          if (config.fillHistograms) {
            const std::array arrMom{pvec0, pvec1, pvec2};
            registry.fill(HIST("hMassDplusToPiKPi"), RecoDecay::m(arrMom, massDplusToPiKPi));
            registry.fill(HIST("hPtCand"), RecoDecay::pt(pVecCand));
            registry.fill(HIST("hVtxX"), secondaryVertex[0]);
          }
        }
      }
    }

    return nCand;
  }

  void processDplus(SelectedCollisions const& collisions,
                    FilteredTrackAssocSel const& /*trackIndices*/,
                    TracksWithCov const& tracks,
                    aod::BCsWithTimestamps const& /*bcWithTimeStamps*/)
  {
    const bool doTiming = config.enableTiming;
    for (const auto& collision : collisions) {
      const auto tStartCollision = StopWatch::now();
      timing.fill(0.);
      counters.fill(0);
      counters[CountCollisions] = 1;
      const auto thisCollId = collision.globalIndex();

      const auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc, runNumber, ccdb, config.isRun2 ? config.ccdbPathGrp : config.ccdbPathGrpMag, lut, config.isRun2);
      const auto bz = o2::base::Propagator::Instance()->getNominalBz();
      df2.setBz(bz);
      df3.setBz(bz);

      const auto groupedPosPions = positivePions->sliceByCached(aod::track::collisionId, thisCollId, cache);
      const auto groupedNegPions = negativePions->sliceByCached(aod::track::collisionId, thisCollId, cache);
      const auto groupedPosKaons = positiveKaons->sliceByCached(aod::track::collisionId, thisCollId, cache);
      const auto groupedNegKaons = negativeKaons->sliceByCached(aod::track::collisionId, thisCollId, cache);

      int nCand{0};
      // D+ -> pi+ K- pi+
      nCand += buildCandidatesLegacy(collision, tracks, groupedPosPions, groupedNegKaons);
      // D- -> pi- K+ pi-
      nCand += buildCandidatesLegacy(collision, tracks, groupedNegPions, groupedPosKaons);

      finaliseCollision(tStartCollision, doTiming, nCand);
    }
  }
  PROCESS_SWITCH(HfDplusCombinatorics, processDplus, "D+ only, original implementation (no prong cache)", true);

  void processDplusCached(SelectedCollisions const& collisions,
                          FilteredTrackAssocSel const& trackIndices,
                          TracksWithCov const& tracks,
                          aod::BCsWithTimestamps const& /*bcWithTimeStamps*/)
  {
    const bool doTiming = config.enableTiming;
    for (const auto& collision : collisions) {
      const auto tStartCollision = StopWatch::now();
      timing.fill(0.);
      counters.fill(0);
      counters[CountCollisions] = 1;

      prepareCollision(collision, trackIndices, tracks, doTiming);

      int nCand{0};
      runTripleLoop(collision, cachePos, cacheNeg, cachePos.pionDplus, cacheNeg.kaonDplus, pairUsefulDplus, hypothesisDplusOnly, doTiming, nCand);
      runTripleLoop(collision, cacheNeg, cachePos, cacheNeg.pionDplus, cachePos.kaonDplus, pairUsefulDplus, hypothesisDplusOnly, doTiming, nCand);

      finaliseCollision(tStartCollision, doTiming, nCand);
    }
  }
  PROCESS_SWITCH(HfDplusCombinatorics, processDplusCached, "D+ only, with the per-collision prong cache", false);

  void processSeparateLoops(SelectedCollisions const& collisions,
                            FilteredTrackAssocSel const& trackIndices,
                            TracksWithCov const& tracks,
                            aod::BCsWithTimestamps const& /*bcWithTimeStamps*/)
  {
    const bool doTiming = config.enableTiming;
    for (const auto& collision : collisions) {
      const auto tStartCollision = StopWatch::now();
      timing.fill(0.);
      counters.fill(0);
      counters[CountCollisions] = 1;

      prepareCollision(collision, trackIndices, tracks, doTiming);

      int nCand{0};
      // D+ -> pi+ K- pi+ and charge conjugate
      runTripleLoop(collision, cachePos, cacheNeg, cachePos.pionDplus, cacheNeg.kaonDplus, pairUsefulDplus, hypothesisDplusOnly, doTiming, nCand);
      runTripleLoop(collision, cacheNeg, cachePos, cacheNeg.pionDplus, cachePos.kaonDplus, pairUsefulDplus, hypothesisDplusOnly, doTiming, nCand);
      // Ds -> K+ K- pi+ and charge conjugate, as a second and independent loop: a triplet good for
      // both species is fitted twice and written twice
      runTripleLoop(collision, cachePos, cacheNeg, cachePos.dsSameSign, cacheNeg.kaonDs, pairUsefulDs, hypothesisDsOnly, doTiming, nCand);
      runTripleLoop(collision, cacheNeg, cachePos, cacheNeg.dsSameSign, cachePos.kaonDs, pairUsefulDs, hypothesisDsOnly, doTiming, nCand);

      finaliseCollision(tStartCollision, doTiming, nCand);
    }
  }
  PROCESS_SWITCH(HfDplusCombinatorics, processSeparateLoops, "One independent triple loop per species", false);

  void processMergedLoop(SelectedCollisions const& collisions,
                         FilteredTrackAssocSel const& trackIndices,
                         TracksWithCov const& tracks,
                         aod::BCsWithTimestamps const& /*bcWithTimeStamps*/)
  {
    const bool doTiming = config.enableTiming;
    for (const auto& collision : collisions) {
      const auto tStartCollision = StopWatch::now();
      timing.fill(0.);
      counters.fill(0);
      counters[CountCollisions] = 1;

      prepareCollision(collision, trackIndices, tracks, doTiming);

      int nCand{0};
      // a single loop over the union of the selected tracks, with a per-triplet hypothesis mask
      runTripleLoop(collision, cachePos, cacheNeg, cachePos.sameSignAny, cacheNeg.oppSignAny, pairUsefulBoth, hypothesisBoth, doTiming, nCand);
      runTripleLoop(collision, cacheNeg, cachePos, cacheNeg.sameSignAny, cachePos.oppSignAny, pairUsefulBoth, hypothesisBoth, doTiming, nCand);

      finaliseCollision(tStartCollision, doTiming, nCand);
    }
  }
  PROCESS_SWITCH(HfDplusCombinatorics, processMergedLoop, "One triple loop over the union of the selected tracks", false);

  void processDummy(aod::Collisions const&)
  {
    // no-op, to be enabled when only the training sample is produced
  }
  PROCESS_SWITCH(HfDplusCombinatorics, processDummy, "Do nothing", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfDplusTrackSelMl>(cfgc),
    adaptAnalysisTask<HfDplusCombinatorics>(cfgc)};
}
