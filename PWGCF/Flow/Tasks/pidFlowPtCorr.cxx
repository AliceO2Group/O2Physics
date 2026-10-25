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

/// \file   pidFlowPtCorr.cxx
/// \author Fuchun Cui(fcui@cern.ch) Qiuyu Xia(qiuyu.xia@cern.ch)
/// \since  Nov/24/2025
/// \brief  This task is to caculate vn-[pt] correlation of PID particles

#include "PWGCF/GenericFramework/Core/FlowContainer.h"
#include "PWGCF/GenericFramework/Core/GFW.h"
#include "PWGCF/GenericFramework/Core/GFWWeights.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/TriggerAliases.h"
#include "Common/CCDB/ctpRateFetcher.h"
#include "Common/Core/TableHelper.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/LHCConstants.h>
#include <CommonConstants/MathConstants.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DataFormatsParameters/GRPObject.h>
#include <DetectorsBase/Propagator.h>
#include <Field/MagneticField.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/DataTypes.h>
#include <Framework/Expressions.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/StringHelpers.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/PID.h>

#include <TAxis.h>
#include <TF1.h>
#include <TGeoGlobalMagField.h>
#include <TH1.h>
#include <TH2.h>
#include <THnSparse.h>
#include <TNamed.h>
#include <TObjArray.h>
#include <TPDGCode.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TRandom3.h>
#include <TString.h>

#include <Rtypes.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <memory>
#include <string>
#include <unordered_set>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{
namespace pidflowtimetrack
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);
DECLARE_SOA_INDEX_COLUMN(Track, track);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(DcaXY, dcaXY, float);
DECLARE_SOA_COLUMN(DcaZ, dcaZ, float);
DECLARE_SOA_COLUMN(HasITS, hasITS, bool);
DECLARE_SOA_COLUMN(HasTPC, hasTPC, bool);
DECLARE_SOA_COLUMN(ItsNCls, itsNCls, uint8_t);
DECLARE_SOA_COLUMN(ItsChi2NCl, itsChi2NCl, float);
DECLARE_SOA_COLUMN(TpcNClsCrossedRows, tpcNClsCrossedRows, int);
DECLARE_SOA_COLUMN(TpcNClsFound, tpcNClsFound, int);
DECLARE_SOA_COLUMN(TpcChi2NCl, tpcChi2NCl, float);
DECLARE_SOA_COLUMN(NumCompatibleCollisions, numCompatibleCollisions, unsigned int);
} // namespace pidflowtimetrack

// One row per accepted (target collision, source track), ordered by collision.
// Kinematics/DCA belong to the target vertex; detector quality belongs to the source.
// Baseline selection is applied here; analysis quality/pT cuts stay in the consumers.
DECLARE_SOA_TABLE(PidFlowTimeTracks, "AOD", "PFTIMETRACK",
                  o2::soa::Index<>, pidflowtimetrack::CollisionId, pidflowtimetrack::TrackId,
                  pidflowtimetrack::Pt, pidflowtimetrack::Eta, pidflowtimetrack::Phi,
                  pidflowtimetrack::DcaXY, pidflowtimetrack::DcaZ,
                  pidflowtimetrack::HasITS, pidflowtimetrack::HasTPC,
                  pidflowtimetrack::ItsNCls, pidflowtimetrack::ItsChi2NCl,
                  pidflowtimetrack::TpcNClsCrossedRows, pidflowtimetrack::TpcNClsFound,
                  pidflowtimetrack::TpcChi2NCl, pidflowtimetrack::NumCompatibleCollisions);
} // namespace o2::aod

namespace
{
// Keep the corrconfigs vector layout in one named index table. Namespace
// grouping avoids storing compile-time constants as analysis-task members.
namespace corrConfigIndex
{
namespace charged
{
enum : std::size_t {
  kC22Gap = 0,
  kC24 = 1,
  kC22Full = 2,
  kC32Gap = 3,
  kC34Gap = 4,
  kC24Gap = 52
};
} // namespace charged

namespace pion
{
enum : std::size_t {
  kC22GapA = 5,
  kC22GapB = 6,
  kC24A = 11,
  kC24B = 12,
  kC32GapA = 17,
  kC32GapB = 18,
  kC34A = 23,
  kC34B = 24,
  kC22Pure = 29,
  kC32Pure = 32,
  kC22FullA = 35,
  kC22FullB = 36,
  kC24PIDChargedGap = 48
};
} // namespace pion

namespace kaon
{
enum : std::size_t {
  kC22GapA = 7,
  kC22GapB = 8,
  kC24A = 13,
  kC24B = 14,
  kC32GapA = 19,
  kC32GapB = 20,
  kC34A = 25,
  kC34B = 26,
  kC22Pure = 30,
  kC32Pure = 33,
  kC22FullA = 37,
  kC22FullB = 38,
  kC24PIDChargedGap = 49
};
} // namespace kaon

namespace proton
{
enum : std::size_t {
  kC22GapA = 9,
  kC22GapB = 10,
  kC24A = 15,
  kC24B = 16,
  kC32GapA = 21,
  kC32GapB = 22,
  kC34A = 27,
  kC34B = 28,
  kC22Pure = 31,
  kC32Pure = 34,
  kC22FullA = 39,
  kC22FullB = 40,
  kC24PIDChargedGap = 50
};
} // namespace proton

// Unidentified correlations are appended after the original configurations
// so the charged and identified-particle vector layout remains unchanged.
namespace unidentified
{
enum : std::size_t {
  kC22GapA = 41,
  kC22GapB = 42,
  kC24A = 43,
  kC24B = 44,
  kC22Pure = 45,
  kC22FullA = 46,
  kC22FullB = 47,
  kC24PIDChargedGap = 51
};
} // namespace unidentified

enum : std::size_t {
  kCount = 53
};
} // namespace corrConfigIndex

// GFW region bitmasks are grouped by particle type and region role. POI and
// overlap use distinct bits because the overlap weight can differ from POI.
namespace gfwBitMask
{
namespace charged
{
enum : int {
  kReference = 1 << 0
};
} // namespace charged

namespace pion
{
enum : int {
  kPOI = 1 << 1,
  kOverlap = 1 << 4
};
} // namespace pion

namespace kaon
{
enum : int {
  kPOI = 1 << 2,
  kOverlap = 1 << 5
};
} // namespace kaon

namespace proton
{
enum : int {
  kPOI = 1 << 3,
  kOverlap = 1 << 6
};
} // namespace proton

namespace unidentified
{
enum : int {
  kPOI = 1 << 7,
  kOverlap = 1 << 8
};
} // namespace unidentified
} // namespace gfwBitMask

// Only target-dependent kinematics change; detector quality and MC labels stay
// on the source row. This view deliberately has no PID-table dependency.
template <typename Track>
struct TimeAssociatedChargedTrack : Track {
  explicit TimeAssociatedChargedTrack(Track const& source)
    : Track(source), targetPt(source.pt()), targetEta(source.eta()), targetPhi(source.phi()), targetDCAxy(source.dcaXY()), targetDCAz(source.dcaZ())
  {
  }
  float targetPt, targetEta, targetPhi, targetDCAxy, targetDCAz;
  [[nodiscard]] float pt() const { return targetPt; }
  [[nodiscard]] float eta() const { return targetEta; }
  [[nodiscard]] float phi() const { return targetPhi; }
  [[nodiscard]] float dcaXY() const { return targetDCAxy; }
  [[nodiscard]] float dcaZ() const { return targetDCAz; }
};
template <typename T>
T timeAssociationOption(InitContext& context, std::string const& task, std::string const& option)
{
  T value{};
  if (!o2::common::core::getTaskOptionValue(context, task, option, value, false)) {
    LOGF(fatal, "Charged time association needs option %s from producer %s", option.c_str(), task.c_str());
  }
  return value;
}

} // namespace

// Shared charged-only producer: public TrackAssoc supplies time compatibility;
// this task evaluates each pair at its target vertex once for flow, NUA and NUE.
struct PidFlowChargedTimeProducer {
  Produces<aod::PidFlowTimeTracks> timeTracks;
  // Original-vertex DCA/eta must not preselect candidates for another vertex.
  using SourceTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA>;
  Preslice<aod::TrackAssoc> timePairsPerCollision = aod::track_association::collisionId;
  Service<ccdb::BasicCCDBManager> ccdb{};
  HistogramRegistry registry{"registry"};
  OutputObj<TNamed> chargedTimeAssociationConfig{TNamed("ChargedTimeAssociationConfig", "disabled")};
  TrackSelection timeGlobalSelection;
  TrackSelection timeSDDSelection;
  bool timeIsRun3 = false;
  int timeFieldRun = -1;
  uint64_t timeFieldTimestamp = 0;
  float timeL3Current = 0.f;
  float timeDipoleCurrent = 0.f;
  bool timeUniformField = false;
  bool cfgRequireUniqueTimeAssociation = false;
  bool cfgTimeAssociationQA = true;
  bool timeUsePVAssociation = true;
  float timeNSigma = 4.f;
  float timeMargin = 0.f;
  float cutEta = 0.8f;
  float cutVertex = 10.f;

  void init(InitContext& context)
  {
    if (doprocessDummy.value == doprocessTimeTracks.value) {
      LOGF(fatal, "Enable exactly one of processDummy and processTimeTracks in pid-flow-charged-time-producer");
    }
    const bool needed = timeAssociationOption<bool>(context, "pid-flow-pt-corr", "processDataTimeAssociation") ||
                        timeAssociationOption<bool>(context, "pid-flow-pt-corr", "fillCorrectionGraphTimeAssociation") ||
                        timeAssociationOption<bool>(context, "pid-flow-pt-corr", "processMCTimeAssociation") ||
                        timeAssociationOption<bool>(context, "pid-flow-pt-corr", "processTimeAssociationComparison");
    if (needed != doprocessTimeTracks.value) {
      LOGF(fatal, "Enable processTimeTracks (and disable processDummy) exactly when a charged time consumer is enabled");
    }
    if (!doprocessTimeTracks) {
      return;
    }
    cfgRequireUniqueTimeAssociation = timeAssociationOption<bool>(context, "pid-flow-pt-corr", "cfgRequireUniqueTimeAssociation");
    cfgTimeAssociationQA = timeAssociationOption<bool>(context, "pid-flow-pt-corr", "cfgTimeAssociationQA");
    timeUsePVAssociation = timeAssociationOption<bool>(context, "track-to-collision-association", "usePVAssociation");
    timeNSigma = timeAssociationOption<float>(context, "track-to-collision-association", "nSigmaForTimeCompat");
    timeMargin = timeAssociationOption<float>(context, "track-to-collision-association", "timeMargin");
    cutEta = timeAssociationOption<float>(context, "pid-flow-pt-corr", "trkQualityOpts.cfgCutEta");
    cutVertex = timeAssociationOption<float>(context, "pid-flow-pt-corr", "cfgCutVertex");
    ccdb->setURL(timeAssociationOption<std::string>(context, "pid-flow-pt-corr", "cfgurl"));
    ccdb->setCaching(true);
    ccdb->setCreatedNotAfter(timeAssociationOption<int64_t>(context, "pid-flow-pt-corr", "cfgnolaterthan"));
    if (!timeAssociationOption<bool>(context, "track-to-collision-association", "processAssocWithTime") ||
        timeAssociationOption<bool>(context, "track-to-collision-association", "processStandardAssoc") ||
        timeAssociationOption<int>(context, "track-to-collision-association", "setTrackSelections") != 0) {
      LOGF(fatal, "Use track-to-collision-association: processAssocWithTime=true, processStandardAssoc=false, setTrackSelections=0");
    }

    // Inherit the effective upstream selection rather than introducing a
    // second set of global-track cuts with different hidden DCA defaults.
    timeIsRun3 = timeAssociationOption<bool>(context, "track-selection", "isRun3");
    const int matching = timeAssociationOption<int>(context, "track-selection", "itsMatching");
    const int dcaSetup = timeAssociationOption<int>(context, "track-selection", "dcaSetup");
    const bool compatibilityIU = timeAssociationOption<bool>(context, "track-selection", "compatibilityIU");
    if (timeIsRun3) {
      using Matching = TrackSelection::GlobalTrackRun3ITSMatching;
      const std::array<int, 5> modes{Matching::Run3ITSibAny, Matching::Run3ITSibAny, Matching::Run3ITSallAny, Matching::Run3ITSall7Layers, Matching::Run3ITSibFirst};
      if (matching < 0 || matching >= static_cast<int>(modes.size())) {
        LOGF(fatal, "Unsupported track-selection itsMatching=%d", matching);
      }
      timeGlobalSelection = getGlobalTrackSelectionRun3ITSMatch(modes[matching], dcaSetup);
      timeGlobalSelection.SetTrackType(compatibilityIU ? aod::track::TrackTypeEnum::TrackIU : aod::track::TrackTypeEnum::Track);
    } else {
      if (matching != 0) {
        LOGF(fatal, "Run 2 time association requires track-selection itsMatching=0");
      }
      timeGlobalSelection = getGlobalTrackSelection();
    }
    timeSDDSelection = getGlobalTrackSelectionSDD();
    const auto ptMin = timeAssociationOption<float>(context, "track-selection", "ptMin");
    const auto ptMax = timeAssociationOption<float>(context, "track-selection", "ptMax");
    const auto etaMin = timeAssociationOption<float>(context, "track-selection", "etaMin");
    const auto etaMax = timeAssociationOption<float>(context, "track-selection", "etaMax");
    for (auto const& selection : {&timeGlobalSelection, &timeSDDSelection}) {
      selection->SetPtRange(ptMin, ptMax);
      selection->SetEtaRange(etaMin, etaMax);
    }
    std::string metadata = "association=TrackAssoc;table=PidFlowTimeTracks;charged-only;PID=none;local-density=none;source=Tracks;propagation=target-DCA,USEMatCorrNONE;PV=not-refitted;unique=" + std::to_string(static_cast<int>(cfgRequireUniqueTimeAssociation));
    metadata += ";isRun3=" + std::to_string(static_cast<int>(timeIsRun3)) + ";itsMatching=" + std::to_string(matching) + ";dcaSetup=" + std::to_string(dcaSetup) + ";compatibilityIU=" + std::to_string(static_cast<int>(compatibilityIU));
    metadata += ";pt=" + std::to_string(ptMin) + "," + std::to_string(ptMax) + ";eta=" + std::to_string(etaMin) + "," + std::to_string(etaMax);
    metadata += ";producer.nSigmaForTimeCompat=" + std::to_string(timeAssociationOption<float>(context, "track-to-collision-association", "nSigmaForTimeCompat"));
    metadata += ";producer.timeMargin=" + std::to_string(timeAssociationOption<float>(context, "track-to-collision-association", "timeMargin"));
    metadata += ";producer.bcWindowForOneSigma=" + std::to_string(timeAssociationOption<int>(context, "track-to-collision-association", "bcWindowForOneSigma"));
    metadata += ";producer.usePVAssociation=" + std::to_string(static_cast<int>(timeAssociationOption<bool>(context, "track-to-collision-association", "usePVAssociation")));
    metadata += ";producer.includeUnassigned=" + std::to_string(static_cast<int>(timeAssociationOption<bool>(context, "track-to-collision-association", "includeUnassigned")));
    chargedTimeAssociationConfig->SetTitle(metadata.c_str());
    LOG(info) << "Charged time association: " << metadata;
    if (!cfgTimeAssociationQA) {
      return;
    }
    registry.add("timeAssociation/candidatesPerTrack", "Time-compatible collisions per source track before event/track cuts;N_{coll}^{comp};Track count", {HistType::kTH1D, {{101, -0.5, 100.5}}});
    registry.add("timeAssociation/candidatesPerTrackSourceBaseline", "Source global/SDD + finite kinematics/DCA + eta only;N_{coll}^{comp};Track count", {HistType::kTH1D, {{101, -0.5, 100.5}}});
    registry.add("timeAssociation/candidatesVsSourceType", "All source tracks;N_{coll}^{comp};Source type", {HistType::kTH2D, {{101, -0.5, 100.5}, {3, -0.5, 2.5}}});
    registry.add("timeAssociation/candidatesVsAmbiguousBC", "AmbiguousTracks BC status (also for assigned tracks);N_{coll}^{comp};Ambiguous BC status", {HistType::kTH2D, {{101, -0.5, 100.5}, {3, -0.5, 2.5}}});
    registry.add("timeAssociation/candidatesVsRawTimeStatus", "Stored track time/resolution before PV substitution;N_{coll}^{comp};Raw time status", {HistType::kTH2D, {{101, -0.5, 100.5}, {4, -0.5, 3.5}}});
    registry.add("timeAssociation/candidatesVsOriginalTimeTest", "Fine time test at original collision; excludes coarse BC/search checks;N_{coll}^{comp};Original time test", {HistType::kTH2D, {{101, -0.5, 100.5}, {6, -0.5, 5.5}}});
    registry.add("timeAssociation/candidatesVsAcceptedTargets", "After unique cut, target vertex cut, propagation and baseline (before consumer cuts);N_{coll}^{comp};N accepted target collisions", {HistType::kTH2D, {{101, -0.5, 100.5}, {101, -0.5, 100.5}}});
    registry.add("timeAssociation/zeroCandidatePtEta", "Zero-candidate tracks with finite source pT/eta;Source p_{T} (GeV/c);Source #eta", {HistType::kTH2D, {{200, 0., 20.}, {100, -5., 5.}}});
    const auto setLabels = [](TAxis* axis, auto const& names) {
      for (std::size_t bin = 0; bin < names.size(); ++bin) {
        axis->SetBinLabel(bin + 1, names[bin]);
      }
    };
    setLabels(registry.get<TH2>(HIST("timeAssociation/candidatesVsSourceType"))->GetYaxis(),
              std::array{"no original collision", "assigned non-PV", "assigned PV contributor"});
    setLabels(registry.get<TH2>(HIST("timeAssociation/candidatesVsAmbiguousBC"))->GetYaxis(),
              std::array{"no ambiguous row", "invalid/empty BC slice", "valid BC slice"});
    setLabels(registry.get<TH2>(HIST("timeAssociation/candidatesVsRawTimeStatus"))->GetYaxis(),
              std::array{"finite time, res > 0", "finite time, res = 0", "finite time, res < 0", "nonfinite time or res"});
    setLabels(registry.get<TH2>(HIST("timeAssociation/candidatesVsOriginalTimeTest"))->GetYaxis(),
              std::array{"no original collision", "invalid effective time/res", "fails fine time test", "passes, original pair absent", "passes, original pair present", "invalid original collision index"});
    registry.add("timeAssociation/pairStatus", "Target-track pairs (stages overlap);stage;pairs", {HistType::kTH1D, {{9, 0, 9}}});
    const std::array<const char*, 9> labels{"producer pair", "duplicate pair", "nonunique rejected", "same collision", "different collision", "unassigned source", "propagation failed", "baseline rejected", "baseline accepted"};
    for (std::size_t bin = 0; bin < labels.size(); ++bin) {
      registry.get<TH1>(HIST("timeAssociation/pairStatus"))->GetXaxis()->SetBinLabel(bin + 1, labels[bin]);
    }
  }

  template <typename GRP>
  void setChargedTimeField(GRP const* grp)
  {
    auto* field = dynamic_cast<o2::field::MagneticField*>(TGeoGlobalMagField::Instance()->GetField());
    if (!field) {
      o2::base::Propagator::initFieldFromGRP(grp);
    } else if (!field->TestBit(o2::field::MagneticField::kOverrideGRP) &&
               (timeFieldRun < 0 || timeL3Current != grp->getL3Current() || timeDipoleCurrent != grp->getDipoleCurrent() || timeUniformField != grp->getFieldUniformity())) {
      // Keep the field object alive: Propagator can already hold its pointer.
      field->rescaleField(grp->getL3Current(), grp->getDipoleCurrent(), grp->getFieldUniformity());
    }
    o2::base::Propagator::Instance()->updateField();
    timeL3Current = grp->getL3Current();
    timeDipoleCurrent = grp->getDipoleCurrent();
    timeUniformField = grp->getFieldUniformity();
  }

  void prepareChargedTimeField(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (timeFieldRun == bc.runNumber() && timeFieldTimestamp == bc.timestamp()) {
      return;
    }
    if (timeIsRun3) {
      auto* grp = ccdb->getForTimeStamp<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", bc.timestamp());
      if (!grp) {
        LOGF(fatal, "Missing GRPMagField for charged time association");
      }
      setChargedTimeField(grp);
    } else {
      auto* grp = ccdb->getForTimeStamp<o2::parameters::GRPObject>("GLO/GRP/GRP", bc.timestamp());
      if (!grp) {
        LOGF(fatal, "Missing GRPObject for charged time association");
      }
      setChargedTimeField(grp);
    }
    timeFieldRun = bc.runNumber();
    timeFieldTimestamp = bc.timestamp();
  }

  template <typename TTracks>
  std::vector<unsigned int> chargedTimeMultiplicity(aod::TrackAssoc const& pairs, TTracks const& tracks,
                                                    aod::Collisions const& collisions, aod::BCsWithTimestamps const& bcs, aod::AmbiguousTracks const& ambiguousTracks)
  {
    std::vector<unsigned int> counts(tracks.size(), 0);
    std::vector<int64_t> lastCollision(tracks.size(), -1);
    std::vector<bool> hasOriginalPair(cfgTimeAssociationQA ? tracks.size() : 0, false);
    for (auto const& pair : pairs) {
      const int64_t index = pair.trackId() - tracks.offset();
      if (index < 0 || index >= static_cast<int64_t>(tracks.size())) {
        LOGF(fatal, "TrackAssoc track index is outside the full source track table");
      }
      // TrackAssoc is grouped by collision. Count distinct target collisions,
      // including those that will later fail event or target-track cuts.
      if (lastCollision[index] != pair.collisionId()) {
        ++counts[index];
        lastCollision[index] = pair.collisionId();
      }
      if (cfgTimeAssociationQA && pair.collisionId() == tracks.rawIteratorAt(index).collisionId()) {
        hasOriginalPair[index] = true;
      }
    }
    if (!cfgTimeAssociationQA) {
      return counts;
    }
    // Cache the first ambiguous row, matching the upstream lookup. Do not
    // dereference its BC slice until its endpoints have been checked.
    std::vector<uint8_t> ambiguousBCStatus(tracks.size(), 0);
    for (auto const& ambiguous : ambiguousTracks) {
      const int64_t index = ambiguous.trackId() - tracks.offset();
      if (index < 0 || index >= static_cast<int64_t>(tracks.size()) || ambiguousBCStatus[index] != 0) {
        continue;
      }
      auto const& ids = ambiguous.bcIds();
      const bool valid = ambiguous.has_bc() && ids[0] >= 0 && ids[1] >= ids[0] && ids[1] < bcs.size();
      ambiguousBCStatus[index] = valid ? 2 : 1;
    }
    // One unweighted entry per source track, including zero candidates, before
    // the unique-association cut or any target event/track selection.
    for (auto const& track : tracks) {
      const int64_t index = track.globalIndex() - tracks.offset();
      const auto count = counts[index];
      registry.fill(HIST("timeAssociation/candidatesPerTrack"), count);
      const int sourceType = !track.has_collision() ? 0 : (track.isPVContributor() ? 2 : 1);
      registry.fill(HIST("timeAssociation/candidatesVsSourceType"), count, sourceType);
      registry.fill(HIST("timeAssociation/candidatesVsAmbiguousBC"), count, ambiguousBCStatus[index]);
      const bool finiteTime = std::isfinite(track.trackTime()) && std::isfinite(track.trackTimeRes());
      const int rawTimeStatus = !finiteTime ? 3 : (track.trackTimeRes() > 0.f ? 0 : (track.trackTimeRes() == 0.f ? 1 : 2));
      registry.fill(HIST("timeAssociation/candidatesVsRawTimeStatus"), count, rawTimeStatus);
      const bool finite = std::isfinite(track.pt()) && std::isfinite(track.eta()) && std::isfinite(track.phi()) && std::isfinite(track.dcaXY()) && std::isfinite(track.dcaZ());
      if (finite && (track.isGlobalTrack() || track.isGlobalTrackSDD()) && std::fabs(track.eta()) < cutEta) {
        // Diagnostic at the source vertex only; not the final flow selection.
        registry.fill(HIST("timeAssociation/candidatesPerTrackSourceBaseline"), count);
      }
      if (count == 0 && std::isfinite(track.pt()) && std::isfinite(track.eta())) {
        registry.fill(HIST("timeAssociation/zeroCandidatePtEta"), track.pt(), track.eta());
      }
      int originalTimeStatus = 0;
      if (track.has_collision()) {
        const int64_t collisionIndex = track.collisionId() - collisions.offset();
        originalTimeStatus = 5;
        if (collisionIndex >= 0 && collisionIndex < collisions.size()) {
          auto collision = collisions.rawIteratorAt(collisionIndex);
          const bool usePVTime = timeUsePVAssociation && track.isPVContributor();
          const float trackTime = usePVTime ? collision.collisionTime() : track.trackTime();
          const float trackRes = usePVTime ? o2::constants::lhc::LHCBunchSpacingNS : track.trackTimeRes();
          const float collisionRes = collision.collisionTimeRes();
          originalTimeStatus = 1;
          if (std::isfinite(trackTime) && std::isfinite(trackRes) && trackRes >= 0.f && std::isfinite(collision.collisionTime()) && std::isfinite(collisionRes) && collisionRes >= 0.f) {
            float threshold = trackRes;
            if (!usePVTime) {
              threshold = TESTBIT(track.flags(), aod::track::TrackTimeResIsRange)
                            ? trackRes + timeNSigma * std::sqrt(collisionRes * collisionRes) + timeMargin
                            : timeNSigma * std::sqrt(collisionRes * collisionRes + trackRes * trackRes) + timeMargin;
            }
            // Same BC as the original collision, hence bcOffset = 0. This
            // checks only fine time compatibility, not the upstream BC cache
            // or iterator window. A pass with no pair needs upstream study.
            originalTimeStatus = std::fabs(trackTime - collision.collisionTime()) < threshold ? (hasOriginalPair[index] ? 4 : 3) : 2;
          }
        }
      }
      registry.fill(HIST("timeAssociation/candidatesVsOriginalTimeTest"), count, originalTimeStatus);
    }
    return counts;
  }

  template <typename TCollision, typename TTracks>
  void fillChargedTimeTracks(TCollision const& collision, aod::TrackAssoc const& pairs, TTracks const& tracks, std::vector<unsigned int> const& counts, std::vector<unsigned int>& acceptedCounts)
  {
    auto collisionPairs = pairs.sliceBy(timePairsPerCollision, collision.globalIndex());
    std::unordered_set<int64_t> seen;
    bool fieldReady = false;
    for (auto const& pair : collisionPairs) {
      if (cfgTimeAssociationQA) {
        registry.fill(HIST("timeAssociation/pairStatus"), 0.5);
      }
      if (!seen.insert(pair.trackId()).second) {
        if (cfgTimeAssociationQA) {
          registry.fill(HIST("timeAssociation/pairStatus"), 1.5);
        }
        continue;
      }
      const int64_t index = pair.trackId() - tracks.offset();
      if (cfgRequireUniqueTimeAssociation && counts[index] != 1) {
        if (cfgTimeAssociationQA) {
          registry.fill(HIST("timeAssociation/pairStatus"), 2.5);
        }
        continue;
      }
      auto source = tracks.rawIteratorAt(index);
      TimeAssociatedChargedTrack<typename TTracks::iterator> candidate(source);
      const bool sameCollision = source.collisionId() == collision.globalIndex();
      if (cfgTimeAssociationQA) {
        registry.fill(HIST("timeAssociation/pairStatus"), sameCollision ? 3.5 : 4.5);
      }
      if (cfgTimeAssociationQA && !source.has_collision()) {
        registry.fill(HIST("timeAssociation/pairStatus"), 5.5);
      }
      if (!sameCollision) {
        if (!fieldReady) {
          prepareChargedTimeField(collision.template bc_as<aod::BCsWithTimestamps>());
          fieldReady = true;
        }
        auto parameters = getTrackPar(source);
        std::array<float, 2> dca{};
        if (!o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, parameters, 2.f, o2::base::Propagator::MatCorrType::USEMatCorrNONE, &dca)) {
          if (cfgTimeAssociationQA) {
            registry.fill(HIST("timeAssociation/pairStatus"), 6.5);
          }
          continue;
        }
        candidate.targetPt = parameters.getPt();
        candidate.targetEta = parameters.getEta();
        candidate.targetPhi = parameters.getPhi();
        candidate.targetDCAxy = dca[0];
        candidate.targetDCAz = dca[1];
      }
      const bool finite = std::isfinite(candidate.pt()) && std::isfinite(candidate.eta()) && std::isfinite(candidate.phi()) && std::isfinite(candidate.dcaXY()) && std::isfinite(candidate.dcaZ());
      const bool baseline = finite && (sameCollision ? (source.isGlobalTrack() || source.isGlobalTrackSDD()) : (timeGlobalSelection.IsSelected(candidate) || (!timeIsRun3 && timeSDDSelection.IsSelected(candidate))));
      if (!baseline || std::fabs(candidate.eta()) >= cutEta) {
        if (cfgTimeAssociationQA) {
          registry.fill(HIST("timeAssociation/pairStatus"), 7.5);
        }
        continue;
      }
      if (cfgTimeAssociationQA) {
        registry.fill(HIST("timeAssociation/pairStatus"), 8.5);
        ++acceptedCounts[index];
      }
      timeTracks(collision.globalIndex(), source.globalIndex(),
                 candidate.pt(), candidate.eta(), candidate.phi(), candidate.dcaXY(), candidate.dcaZ(),
                 source.hasITS(), source.hasTPC(), source.itsNCls(), source.itsChi2NCl(),
                 source.tpcNClsCrossedRows(), source.tpcNClsFound(), source.tpcChi2NCl(), counts[index]);
    }
  }

  void processTimeTracks(aod::Collisions const& collisions, aod::BCsWithTimestamps const& bcs,
                         aod::TrackAssoc const& pairs, SourceTracks const& sources, aod::AmbiguousTracks const& ambiguousTracks)
  {
    // Count all producer candidates before vertex, event or target-track cuts.
    auto counts = chargedTimeMultiplicity(pairs, sources, collisions, bcs, ambiguousTracks);
    std::vector<unsigned int> acceptedCounts(cfgTimeAssociationQA ? sources.size() : 0, 0);
    for (auto const& collision : collisions) {
      if (std::fabs(collision.posZ()) < cutVertex) {
        fillChargedTimeTracks(collision, pairs, sources, counts, acceptedCounts);
      }
    }
    if (cfgTimeAssociationQA) {
      for (std::size_t index = 0; index < counts.size(); ++index) {
        registry.fill(HIST("timeAssociation/candidatesVsAcceptedTargets"), counts[index], acceptedCounts[index]);
      }
    }
  }
  PROCESS_SWITCH(PidFlowChargedTimeProducer, processTimeTracks, "Build shared charged time-associated track table; no PID", false);

  void processDummy(aod::Collisions const&) {}
  PROCESS_SWITCH(PidFlowChargedTimeProducer, processDummy, "Idle entry for original collisionId analysis", true);
};

struct PidFlowPtCorr {
  // configurable
  double minVal4Float = 1e-3;

  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCutChi2prTPCcls{"cfgCutChi2prTPCcls", 2.5, "Chi2 per TPC clusters"};
  Configurable<int> cfgCentEstimator{"cfgCentEstimator", 0, "0:FT0C; 1:FT0CVariant1; 2:FT0M; 3:FV0A; 4:NTPV; 5:NGlobal; 6:MFT"};
  Configurable<bool> cfgUseMultiplicity{"cfgUseMultiplicity", false, "Use raw PV-contributor count in |eta| < 0.8 instead of centrality for output axes and NUE lookup; configure axisMultiplicity and disable cfgDoMultPVCut"};
  Configurable<bool> cfgRequireUniqueTimeAssociation{"cfgRequireUniqueTimeAssociation", false, "Time entries: require exactly one producer candidate collision before event/track cuts"};
  Configurable<bool> cfgTimeAssociationQA{"cfgTimeAssociationQA", true, "Create and fill time-association producer, comparison and MC QA; keep event counters and flow/NUA/NUE outputs"};
  ConfigurableAxis cfgaxisAssociationNch{"cfgaxisAssociationNch", {601, -5., 6005.}, "Association QA Nch (default bin width: 10 tracks; includes zero)"};
  ConfigurableAxis cfgaxisAssociationChangedNch{"cfgaxisAssociationChangedNch", {1001, -0.5, 1000.5}, "Association QA lost/gained track count"};
  ConfigurableAxis cfgaxisAssociationDeltaNch{"cfgaxisAssociationDeltaNch", {4001, -2000.5, 2000.5}, "Association QA Nch(time) minus Nch(original)"};
  ConfigurableAxis cfgaxisAssociationOccupancy{"cfgaxisAssociationOccupancy", {200, -0.5, 19999.5}, "Association QA track occupancy (-1 goes to underflow)"};

  struct : ConfigurableGroup {
    std::string prefix = "trkQualityOpts";
    // track selections
    Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "Eta range for tracks"};
    Configurable<float> cfgRangeEta{"cfgRangeEta", 0.4f, "Independent |eta| half-width for mean-pT particles"};
    Configurable<float> cfgCutPtMin{"cfgCutPtMin", 0.2f, "Minimal pT for ref tracks"};
    Configurable<float> cfgCutPtMax{"cfgCutPtMax", 10.0f, "Maximal pT for ref tracks"};
    Configurable<float> cfgCutPtMinPi{"cfgCutPtMinPi", 0.2f, "Minimal pT for pion POI"};
    Configurable<float> cfgCutPtMaxPi{"cfgCutPtMaxPi", 10.0f, "Maximal pT for pion POI"};
    Configurable<float> cfgCutPtMinKa{"cfgCutPtMinKa", 0.2f, "Minimal pT for kaon POI"};
    Configurable<float> cfgCutPtMaxKa{"cfgCutPtMaxKa", 10.0f, "Maximal pT for kaon POI"};
    Configurable<float> cfgCutPtMinPr{"cfgCutPtMinPr", 0.2f, "Minimal pT for proton POI"};
    Configurable<float> cfgCutPtMaxPr{"cfgCutPtMaxPr", 10.0f, "Maximal pT for proton POI"};
    // track quality selections for daughter track
    Configurable<int> cfgITSNCls{"cfgITSNCls", 5, "check minimum number of ITS clusters"};
    Configurable<int> cfgTPCNCls{"cfgTPCNCls", 50, "check minimum number of TPC hits"};
    Configurable<int> cfgTPCCrossedRows{"cfgTPCCrossedRows", 70, "check minimum number of TPC crossed rows"};
    Configurable<double> cfgITSChi2NDF{"cfgITSChi2NDF", 2.5, "check ITS Chi2NDF"};
    Configurable<bool> cfgCheckGlobalTrack{"cfgCheckGlobalTrack", false, "check global track"};
  } trkQualityOpts;

  struct : ConfigurableGroup {
    std::string prefix = "evtSelOpts";
    Configurable<bool> cfgDoTVXinTRD{"cfgDoTVXinTRD", false, "check kTVXinTRD"};
    Configurable<bool> cfgDoNoTimeFrameBorder{"cfgDoNoTimeFrameBorder", true, "check kNoTimeFrameBorder"};
    Configurable<bool> cfgDoNoITSROFrameBorder{"cfgDoNoITSROFrameBorder", true, "check kNoITSROFrameBorder"};
    Configurable<bool> cfgDoNoSameBunchPileup{"cfgDoNoSameBunchPileup", true, "check kNoITSROFrameBorder"};
    Configurable<bool> cfgDoIsGoodZvtxFT0vsPV{"cfgDoIsGoodZvtxFT0vsPV", true, "check kIsGoodZvtxFT0vsPV"};
    Configurable<bool> cfgDoNoCollInTimeRangeStandard{"cfgDoNoCollInTimeRangeStandard", true, "check kNoCollInTimeRangeStandard"};
    Configurable<bool> cfgDoIsGoodITSLayersAll{"cfgDoIsGoodITSLayersAll", true, "check kIsGoodITSLayersAll"};
    Configurable<int> cfgCutOccupancyHigh{"cfgCutOccupancyHigh", 3000, "High cut on TPC occupancy"};
    Configurable<bool> cfgDoMultPVCut{"cfgDoMultPVCut", true, "do multNTracksPV vs cent cut"};
    Configurable<std::vector<float>> cfgMultPVCut{"cfgMultPVCut", (std::vector<float>{3074.43, -106.192, 1.46176, -0.00968364, 2.61923e-05, 182.128, -7.43492, 0.193901, -0.00256715, 1.22594e-05}), "Used MultPVCut function parameter"};
    Configurable<bool> cfgDoV0AT0Acut{"cfgDoV0AT0Acut", true, "do V0A-T0A cut"};
    Configurable<float> cfgCutminIR{"cfgCutminIR", -1, "cut min IR"};
    Configurable<float> cfgCutmaxIR{"cfgCutmaxIR", 3000, "cut max IR"};
    Configurable<std::string> cfgIRSource{"cfgIRSource", "ZNC hadronic", "Interaction-rate source for CTP rate fetcher. Use T0VTX for pp/OO and ZNC hadronic for Pb-Pb; empty disables fetching"};
  } evtSeleOpts;

  struct : ConfigurableGroup {
    std::string prefix = "dcaCutOpts";
    TF1* fPtDepDCAxy = nullptr;
    Configurable<float> cfgCutDCAz{"cfgCutDCAz", 2.0f, "max DCA to vertex z"};
    Configurable<float> cfgDCAxyNSigma{"cfgDCAxyNSigma", 0, "0: disable; Cut on number of sigma deviations from expected DCA in the transverse direction, nsigma=7 is the same with global track"};
    Configurable<std::string> cfgDCAxy{"cfgDCAxy", "(0.0026+0.005/(x^1.01))", "Functional form of pt-dependent DCAxy cut"};
  } dcaCutOpts;

  struct : ConfigurableGroup {
    std::string prefix = "correctionPathOpts";
    Configurable<std::vector<std::string>> cfgAcceptancePath{"cfgAcceptancePath", (std::vector<std::string>{"Users/f/fcui/NUA/NUAREFPartical", "Users/f/fcui/NUA/NUAK0s", "Users/f/fcui/NUA/NUALambda", "Users/f/fcui/NUA/NUAXi", "Users/f/fcui/NUA/NUAOmega"}), "CCDB path to acceptance object"};

    Configurable<std::vector<std::string>> cfgAcceptancePathWithPt{"cfgAcceptancePathWithPt", (std::vector<std::string>{"Users/f/fcui/NUA/NUAREFParticalWithPt"}), "CCDB path to acceptance THnSparse(phi,eta,vz,pt) object, same content as cfgAcceptancePath but with an extra pt axis"};

    Configurable<std::vector<std::string>> cfgEfficiencyPath{"cfgEfficiencyPath", (std::vector<std::string>{"PathtoRef"}), "CCDB path to efficiency object"};
  } correctionPathOpts;

  Configurable<std::vector<int>> cfgRunNumbers{"cfgRunNumbers", (std::vector<int>{544095, 544098, 544116, 544121, 544122, 544123, 544124}), "Preconfigured run numbers"};
  Configurable<float> cfgEtaGap{"cfgEtaGap", 0.4, "Eta-gap half-width for flow particles; 0.4 selects (-etaMax,-0.4) and (0.4,etaMax)"};
  Configurable<int> cfgFlowNbootstrap{"cfgFlowNbootstrap", 30, "Number of subsamples for bootstrap"};

  // switch
  struct : ConfigurableGroup {
    std::string prefix = "switchsOpts";
    Configurable<bool> cfgOutputrunbyrun{"cfgOutputrunbyrun", false, "OPEN IF USE FUNCTION(fillcorrectiongraph) Fill and output NUA weights run by run"};
    Configurable<bool> cfgOutPutMC{"cfgOutPutMC", false, "Fill MC graphs, note that if the processMCgen is open,this MUST be open"};
    Configurable<bool> cfgOutputQA{"cfgOutputQA", false, "OPEN IF USE FUNCTION(detectorPidQa) do QA"};
    Configurable<bool> cfgDebugMyCode{"cfgDebugMyCode", false, "output some graph for code debug"};
    Configurable<bool> cfgOutPutPtSpectra{"cfgOutPutPtSpectra", false, "output pt spectra for data, MC and RECO"};
    Configurable<bool> cfgCheck2MethodDiff{"cfgCheck2MethodDiff", false, "check difference between v2' && v2''"};
    Configurable<int> cfgClosureTest{"cfgClosureTest", 0, "Nonzero enables extra pion mean-pT weight QA profiles"};
    Configurable<bool> cfgOutPutMC1D{"cfgOutPutMC1D", true, "Fill MC graphs, note that if the processMCgen is open,this MUST be open"};
    Configurable<bool> cfgAddPidResponseMatrixHistograms{"cfgAddPidResponseMatrixHistograms", false, "Add PID response matrix histograms; enable together with processPidResponseMatrix"};
    Configurable<bool> cfgAddC22DeltaPtHistograms{"cfgAddC22DeltaPtHistograms", false, "Add histograms for processDataC22DeltaPt; enable together with processDataC22DeltaPt"};
    Configurable<bool> cfgAddPidPtCorrelationProfiles{"cfgAddPidPtCorrelationProfiles", false, "Add FlowContainer profiles for processPidPtCorrelations; enable together with processPidPtCorrelations"};
    Configurable<bool> cfgC22DeltaPtUsePure{"cfgC22DeltaPtUsePure", false, "true: use PID POI-POI Pure profiles; false: use PID POI-ref and ref-ref profiles"};
    Configurable<bool> cfgAddMeanPtCentNbsHistograms{"cfgAddMeanPtCentNbsHistograms", false, "Add meanptCentNbs TProfile3D histograms"};
    // These switches suppress result allocation and filling only. Charged reference
    // Q-vectors are still calculated when a PID or unidentified output needs them.
    Configurable<bool> cfgOutputCharged{"cfgOutputCharged", true, "Output charged-particle flow, mean-pT, and c22DeltaPt observables"};
    Configurable<bool> cfgOutputPion{"cfgOutputPion", true, "Output pion flow, mean-pT, and c22DeltaPt observables"};
    Configurable<bool> cfgOutputKaon{"cfgOutputKaon", true, "Output kaon flow, mean-pT, and c22DeltaPt observables"};
    Configurable<bool> cfgOutputProton{"cfgOutputProton", true, "Output proton flow, mean-pT, and c22DeltaPt observables"};
    Configurable<bool> cfgOutputUnidentified{"cfgOutputUnidentified", true, "Output PID-failed or ambiguous flow, mean-pT, and c22DeltaPt observables"};

    Configurable<bool> cfgProcessQAOutput{"cfgProcessQAOutput", false, "QA plots for processQA"};
    Configurable<bool> cfgUseNUAWithPt{"cfgUseNUAWithPt", false, "false: use original GFWWeights (phi,eta,vz) NUA; true: use THnSparse (phi,eta,vz,pt) NUA loaded from cfgAcceptancePathWithPt"};
  } switchsOpts;

  // TOF momentum ranges used by circle PID.
  struct : ConfigurableGroup {
    std::string prefix = "pidPtRangeOpts";
    Configurable<float> cfgPtMin4TOFPiKa{"cfgPtMin4TOFPiKa", 0.5, "pt min for TOF to separate pi and k"};
    Configurable<float> cfgPtMax4TOFPiKa{"cfgPtMax4TOFPiKa", 3.0, "pt max for TOF to separate pi and k"};
    Configurable<float> cfgPtMin4TOFKaPr{"cfgPtMin4TOFKaPr", 0.6, "pt min for TOF to separate k and p"};
    Configurable<float> cfgPtMax4TOFKaPr{"cfgPtMax4TOFKaPr", 3.0, "pt max for TOF to separate k and p"};
  } pidPtRangeOpts;

  struct : ConfigurableGroup {
    std::string prefix = "circleCutOpts";
    // TOF pT threshold: above this pt, TOF signal is required
    Configurable<float> cfgCircleCutTofPtCut{"cfgCircleCutTofPtCut", 0.5, "TOF pT threshold for requiring TOF signal"};

    // Circular cut (TPC+TOF) values (radius squared, default = 2^2 = 4)
    Configurable<float> cfgCircleCutSigmaPi{"cfgCircleCutSigmaPi", 4.0, "Circular cut radius squared for Pion (TPC+TOF)"};
    Configurable<float> cfgCircleCutSigmaKa{"cfgCircleCutSigmaKa", 4.0, "Circular cut radius squared for Kaon (TPC+TOF)"};
    Configurable<float> cfgCircleCutSigmaPr{"cfgCircleCutSigmaPr", 4.0, "Circular cut radius squared for Proton (TPC+TOF)"};

    // TPC fallback cuts, combined with ITS below cfgCircleCutITSPtMax.
    Configurable<float> cfgCircleCutTPCPi{"cfgCircleCutTPCPi", 2.0, "TPC fallback nsigma cut for Pion"};
    Configurable<float> cfgCircleCutTPCKa{"cfgCircleCutTPCKa", 2.0, "TPC fallback nsigma cut for Kaon"};
    Configurable<float> cfgCircleCutTPCPr{"cfgCircleCutTPCPr", 2.0, "TPC fallback nsigma cut for Proton"};
    Configurable<float> cfgCircleCutITSPtMax{"cfgCircleCutITSPtMax", 0.5, "Require ITS PID below this pT (GeV/c); 0 disables"};
    Configurable<float> cfgCircleCutITSPi{"cfgCircleCutITSPi", 2.0, "Low-pT ITS absolute nsigma cut for Pion"};
    Configurable<float> cfgCircleCutITSKa{"cfgCircleCutITSKa", 2.0, "Low-pT ITS absolute nsigma cut for Kaon"};
    Configurable<float> cfgCircleCutITSPr{"cfgCircleCutITSPr", 2.0, "Low-pT ITS absolute nsigma cut for Proton"};
  } circleCutOpts;

  Configurable<std::string> cfgTPCPIDClustersMinFunction{"cfgTPCPIDClustersMinFunction", "0", "Minimum TPC PID clusters as a ROOT formula; x = tpcNClsFound(), numeric coefficients, 0 imposes no lower bound"};
  std::unique_ptr<TF1> fTPCPIDClustersMin;

  struct : ConfigurableGroup {
    std::string prefix = "particleAbundanceOpts";
    ConfigurableAxis cfgaxisAbundancePi{"cfgaxisAbundancePi", {100, 0, 1100}, "axis for Abundance Pi"};
    ConfigurableAxis cfgaxisAbundanceKa{"cfgaxisAbundanceKa", {100, 0, 200}, "axis for Abundance ka"};
    ConfigurableAxis cfgaxisAbundancePr{"cfgaxisAbundancePr", {100, 0, 50}, "axis for Abundance Pr"};

    Configurable<bool> cfgOutPutAbundanceDis{"cfgOutPutAbundanceDis", false, "out put hists for pid particle Abundance QA"};
    Configurable<bool> cfgOutputAbundanceYields{"cfgOutputAbundanceYields", false, "Output detector-PID abundance yields for pion, kaon, proton, and unidentified"};
  } particleAbundanceOpts;

  ConfigurableAxis cfgaxisVertex{"cfgaxisVertex", {20, -10, 10}, "vertex axis for histograms"};
  ConfigurableAxis cfgaxisPhi{"cfgaxisPhi", {60, 0.0, constants::math::TwoPI}, "phi axis for histograms"};
  ConfigurableAxis cfgaxisEta{"cfgaxisEta", {40, -1., 1.}, "eta axis for histograms"};
  ConfigurableAxis cfgaxisPt{"cfgaxisPt", {VARIABLE_WIDTH, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.20, 2.40, 2.60, 2.80, 3.00, 3.50, 4.00, 4.50, 5.00, 5.50, 6.00, 10.0}, "pt (GeV)"};
  ConfigurableAxis cfgaxisTPCDedx{"cfgaxisTPCDedx", {500, 0.0f, 1000.0f}, "d#it{E}/d#it{x} axis for detector-PID QA"};
  ConfigurableAxis cfgaxisTOFDelta{"cfgaxisTOFDelta", {200, -1000.0f, 1000.0f}, "TOF #it{t} - #it{t}_{ev} - #it{t}_{exp} axis for detector-PID QA"};
  ConfigurableAxis cfgaxisTOFTime{"cfgaxisTOFTime", {10000, 0.0f, 2.0e6f}, "TOF #it{t} - #it{t}_{ev} (ps) axis for detector-PID QA"};
  ConfigurableAxis cfgaxisTPCPIDClusters{"cfgaxisTPCPIDClusters", {201, -0.5, 200.5}, "TPC clusters used for PID"};
  ConfigurableAxis cfgaxisPIDNSigma{"cfgaxisPIDNSigma", {240, -12., 12.}, "Signed TPC n-sigma axis for PID-quality QA"};
  ConfigurableAxis cfgaxisPIDMomentum{"cfgaxisPIDMomentum", {1000, 0., 10.}, "Total momentum for ITS cluster-size QA (GeV/c)"};
  ConfigurableAxis cfgaxisITSNSigma{"cfgaxisITSNSigma", {161, -12.075, 12.075}, "Signed ITS n-sigma; zero at bin center"};
  ConfigurableAxis cfgaxisITSPt{"cfgaxisITSPt", {300, 0., 3.}, "Fine pT axis for ITS response QA (GeV/c)"};
  ConfigurableAxis cfgaxisITSClusterSize{"cfgaxisITSClusterSize", {101, -0.075, 15.075}, "Truncated mean ITS cluster size divided by cosh(eta)"};
  ConfigurableAxis cfgaxisMeanPt{"cfgaxisMeanPt", {300, 0, 3}, "pt (GeV)"};
  ConfigurableAxis cfgaxisNch{"cfgaxisNch", {3000, 0.5, 3000.5}, "Nch"};
  ConfigurableAxis cfgaxisRun{"cfgaxisRun", {7, 0, 7}, "axis of runs in the data"};

  struct : ConfigurableGroup {
    std::string prefix = "meanptC22GraphOpts";
    ConfigurableAxis cfgaxisBootstrap{"cfgaxisBootstrap", {30, 0, 30}, "cfgaxisBootstrap"};
  } meanptC22GraphOpts;

  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {VARIABLE_WIDTH, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 80, 90}, "Centrality (%) or raw PV multiplicity when cfgUseMultiplicity is enabled"};

  // configurable

  // filter
  // filter and using
  // data
  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter trackFilter = ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t)true)) && (nabs(aod::track::eta) < trkQualityOpts.cfgCutEta.value);

  using TracksPID = soa::Join<aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr, aod::pidTOFPi, aod::pidTOFKa, aod::pidTOFPr>;
  using TracksPIDForDetectorPidQA = soa::Join<aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>;
  // data tracks filter
  using AodTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, o2::aod::TrackSelectionExtension, aod::TracksExtra, TracksPID, aod::TracksIU, aod::TracksDCA>>;
  using AodTracksForDetectorPidQA = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, o2::aod::TrackSelectionExtension, aod::TracksExtra, aod::TOFSignal, aod::TOFEvTime, TracksPIDForDetectorPidQA, aod::TracksIU, aod::TracksDCA>>;
  // data collisions filter
  using AodCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::CentFT0CVariant1s, aod::CentFT0Ms, aod::CentFV0As, aod::CentNTPVs, aod::CentNGlobals, aod::CentMFTs, aod::MultsRun3>>;

  using ChargedTimeCollisionsMC = soa::Filtered<soa::Join<AodCollisions, aod::McCollisionLabels>>;

  // MC
  Filter mccollisionFilter = nabs(aod::mccollision::posZ) < cfgCutVertex;
  using FilteredMcCollisions = soa::Filtered<aod::McCollisions>;

  Filter particleFilter = (nabs(aod::mcparticle::eta) < trkQualityOpts.cfgCutEta.value);
  using FilteredMcParticles = soa::Filtered<aod::McParticles>;

  using FilteredTracksWithMCLabel = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA, aod::McTrackLabels, TracksPID>>;
  using FilteredCollisionsWithMCLabel = soa::Filtered<soa::Join<AodCollisions, aod::EvSels, aod::McCollisionLabels>>;
  // end using and filter

  // filter

  // others

  Preslice<aod::Tracks> perCollision = aod::track::collisionId;
  Preslice<aod::PidFlowTimeTracks> timeTracksPerCollision = aod::pidflowtimetrack::collisionId;
  Preslice<aod::McParticles> particlesPerMCCollision = aod::mcparticle::mcCollisionId;

  // Connect to ccdb
  Service<ccdb::BasicCCDBManager> ccdb{};
  ctpRateFetcher rateFetcher;
  Configurable<int64_t> cfgnolaterthan{"cfgnolaterthan", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};
  Configurable<std::string> cfgurl{"cfgurl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  // Define output
  HistogramRegistry registry{"registry"};

  // val used for bootstrap
  TRandom3* fRndm = new TRandom3(0);
  OutputObj<FlowContainer> fFCCh{FlowContainer("FlowContainerCharged")};
  OutputObj<FlowContainer> fFCPi{FlowContainer("FlowContainerPi")};
  OutputObj<FlowContainer> fFCKa{FlowContainer("FlowContainerKa")};
  OutputObj<FlowContainer> fFCPr{FlowContainer("FlowContainerPr")};
  OutputObj<FlowContainer> fFCUnidentified{FlowContainer("FlowContainerUnidentified")};
  OutputObj<FlowContainer> fFCPidPtCorr{FlowContainer("FlowContainerPidPtCorr")};
  bool timeAcceptanceLoaded = false;
  bool timeEfficiencyLoaded = false;
  // end val used for bootstrap

  // define global variables
  GFW* fGFW = new GFW(); // GFW class used from main src
  std::vector<GFW::CorrConfig> corrconfigs;

  std::vector<std::string> cfgAcceptance;
  std::vector<std::string> cfgAcceptanceWithPt;

  std::vector<float> cfgMultPVCutPara;
  std::vector<int> runNumbers;

  enum MyParticleType {
    kCharged = 0,
    kPion,
    kKaon,
    kProton,
    kUnidentified,
    kNumberOfParticles
  };

  bool isParticleOutputEnabled(MyParticleType type)
  {
    switch (type) {
      case MyParticleType::kCharged:
        return switchsOpts.cfgOutputCharged.value;
      case MyParticleType::kPion:
        return switchsOpts.cfgOutputPion.value;
      case MyParticleType::kKaon:
        return switchsOpts.cfgOutputKaon.value;
      case MyParticleType::kProton:
        return switchsOpts.cfgOutputProton.value;
      case MyParticleType::kUnidentified:
        return switchsOpts.cfgOutputUnidentified.value;
      default:
        return false;
    }
  }

  enum CentEstimators {
    kCentFT0C = 0,
    kCentFT0CVariant1,
    kCentFT0M,
    kCentFV0A,
    kCentNTPV,
    kCentNGlobal,
    kCentMFT,
    // Count the total number of enum
    kCount_CentEstimators
  };

  enum OutputTH1Names {
    // here are TProfiles for vn-pt correlations that are not implemented in GFW
    hPhi = 0,
    hPhicorr,
    kCount_TH1Names
  };

  enum OutputTH3Names {
    hPhiEtaVtxz = 0,
    kCount_TH3Names
  };

  enum MyFunctionName {
    funcProcessData = 0,
    funcDetectorPidQA,
    funcFillCorrectionGraph,
    funcProcessReco,
    funcProcessSim,
    funcProcessQA,
    funcProcessPidResponseMatrix,
    funcProcessDataC22DeltaPt,
    funcProcessPidPtCorrelations,
    funcProcessDataTimeAssociation,
    funcFillCorrectionGraphTimeAssociation,
    funcProcessMCTimeAssociation,
    funcProcessTimeAssociationComparison,
    funcNumber
  };

  enum PidResponseBin {
    kPidRespPion = 0,
    kPidRespKaon = 1,
    kPidRespProton = 2,
    kPidRespOther = 3,
    kPidRespNTrueBins = 4,
    kPidRespUnidentified = kPidRespOther,
    kPidRespNRecoBins = kPidRespNTrueBins
  };

  // graphs for NUE / NUA
  std::vector<GFWWeights*> mAcceptance;
  std::vector<TH2*> mEfficiency;
  // NUA with extra pt axis: THnSparse(phi, eta, vz, pt), same content as mAcceptance but pt-differential
  std::vector<THnSparse*> mAcceptanceWithPt;

  bool correctionsLoaded = false;

  TF1* fMultPVCutLow = nullptr;
  TF1* fMultPVCutHigh = nullptr;
  TF1* fT0AV0AMean = nullptr;
  TF1* fT0AV0ASigma = nullptr;

  // Index order is pion, kaon, proton; stage 0 is before PID, stage 1 is
  // after the current unique PID assignment. The registry owns the objects.
  struct PidQualityQAHists {
    std::array<TH2*, 2> tpcPIDClustersVsPt{};
    std::array<THnSparse*, 2> tpcNSigmaVsPIDClusters{};
    std::array<TH2*, 2> itsNSigmaVsPt{};
    std::array<THnSparse*, 2> tpcVsITSNSigma{};
    TH2* itsClusterSizeVsP = nullptr;
    TH2* itsClusterSizeVsPt = nullptr;
    TProfile* meanITSClusterSizeVsP = nullptr;
    TProfile* meanITSClusterSizeVsPt = nullptr;
  };
  std::array<PidQualityQAHists, 3> pidQualityQAHists{};

  const char* getActivityAxisTitle() const
  {
    return cfgUseMultiplicity.value ? "N_{tracks}^{PV} (|#eta| < 0.8)" : "Centrality (%)";
  }

  void initPidQualityQA()
  {
    const AxisSpec ptAxis{cfgaxisPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec itsPtAxis{cfgaxisITSPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec etaAxis{cfgaxisEta, "#eta"};
    const AxisSpec pAxis{cfgaxisPIDMomentum, "#it{p} (GeV/#it{c})"};
    const AxisSpec centAxis{axisMultiplicity, getActivityAxisTitle()};
    const AxisSpec clusterAxis{cfgaxisTPCPIDClusters, "TPC clusters used for PID"};
    const AxisSpec foundClusterAxis{cfgaxisTPCPIDClusters, "TPC found clusters"};
    const AxisSpec tpcNSigmaAxis{cfgaxisPIDNSigma, "n#sigma_{TPC}"};
    const AxisSpec itsNSigmaAxis{cfgaxisITSNSigma, "n#sigma_{ITS}"};
    const AxisSpec itsSizeAxis{cfgaxisITSClusterSize, "#LTITS cluster size#GT_{trunc}/cosh(#eta)"};
    registry.add("DetectorPidQA/PIDQuality/tpcPIDClustersAll", "Before PID;#it{p}_{T} (GeV/#it{c});TPC clusters used for PID", {HistType::kTH2F, {ptAxis, clusterAxis}});
    registry.add("DetectorPidQA/PIDQuality/tpcPIDClustersUnidentified", "PID failed or ambiguous;#it{p}_{T} (GeV/#it{c});TPC clusters used for PID", {HistType::kTH2F, {ptAxis, clusterAxis}});
    registry.add("DetectorPidQA/PIDQuality/tpcPIDClustersVsFound", "Before PID;TPC found clusters;TPC clusters used for PID", {HistType::kTH2F, {foundClusterAxis, clusterAxis}});
    registry.add("DetectorPidQA/PIDQuality/itsClusterSizeAll", "Before PID;#it{p} (GeV/#it{c});#LTITS cluster size#GT_{trunc}/cosh(#eta)", {HistType::kTH2F, {pAxis, itsSizeAxis}});
    registry.add("DetectorPidQA/PIDQuality/itsClusterSizeAllVsPt", "Before PID;#it{p}_{T} (GeV/#it{c});#LTITS cluster size#GT_{trunc}/cosh(#eta)", {HistType::kTH2F, {itsPtAxis, itsSizeAxis}});
    registry.add("DetectorPidQA/PIDQuality/itsClusterSizeVsPEta", "Before PID: ITS response versus momentum and eta", {HistType::kTHnSparseF, {pAxis, etaAxis, itsSizeAxis}});
    registry.add("DetectorPidQA/PIDQuality/itsResponseStatus", "ITS QA validity;status;tracks", {HistType::kTH1D, {{4, 0., 4.}}});
    auto* statusAxis = registry.get<TH1>(HIST("DetectorPidQA/PIDQuality/itsResponseStatus"))->GetXaxis();
    statusAxis->SetBinLabel(1, "Fewer than 2 cluster sizes");
    statusAxis->SetBinLabel(2, "Invalid momentum or eta");
    statusAxis->SetBinLabel(3, "Nonfinite ITS response");
    statusAxis->SetBinLabel(4, "All ITS hypotheses finite");

    constexpr std::array<const char*, 3> species{"Pi", "Ka", "Pr"};
    constexpr std::array<const char*, 2> stages{"BeforePID", "AfterPID"};
    std::array<TProfile*, 3> expectedITSResponse{};
    for (std::size_t i = 0; i < species.size(); ++i) {
      auto& h = pidQualityQAHists[i];
      for (std::size_t stage = 0; stage < stages.size(); ++stage) {
        const std::string path = std::string("DetectorPidQA/PIDQuality/") + stages[stage] + "/" + species[i] + "/";
        h.tpcPIDClustersVsPt[stage] = registry.add<TH2>(path + "tpcPIDClustersVsPt", "TPC PID cluster count", HistType::kTH2F, {ptAxis, clusterAxis}).get();
        h.tpcNSigmaVsPIDClusters[stage] = registry.add<THnSparse>(path + "tpcNSigmaVsPIDClusters", "TPC PID response versus PID cluster count", HistType::kTHnSparseF, {ptAxis, centAxis, clusterAxis, tpcNSigmaAxis}).get();
        h.itsNSigmaVsPt[stage] = registry.add<TH2>(path + "itsNSigmaVsPt", "ITS response for this mass hypothesis", HistType::kTH2F, {itsPtAxis, itsNSigmaAxis}).get();
        h.tpcVsITSNSigma[stage] = registry.add<THnSparse>(path + "tpcVsITSNSigma", "TPC versus ITS response for this mass hypothesis", HistType::kTHnSparseF, {ptAxis, centAxis, tpcNSigmaAxis, itsNSigmaAxis}).get();
      }
      const std::string path = std::string("DetectorPidQA/PIDQuality/AfterPID/") + species[i] + "/";
      h.itsClusterSizeVsP = registry.add<TH2>(path + "itsClusterSizeVsP", "ITS cluster size after current PID", HistType::kTH2F, {pAxis, itsSizeAxis}).get();
      h.itsClusterSizeVsPt = registry.add<TH2>(path + "itsClusterSizeVsPt", "ITS cluster size after current PID", HistType::kTH2F, {itsPtAxis, itsSizeAxis}).get();
      h.meanITSClusterSizeVsP = registry.add<TProfile>(path + "meanITSClusterSizeVsP", "Measured mean ITS cluster size after current PID;#it{p} (GeV/#it{c});#LTITS signal#GT", HistType::kTProfile, {pAxis}).get();
      h.meanITSClusterSizeVsPt = registry.add<TProfile>(path + "meanITSClusterSizeVsPt", "Measured mean ITS cluster size after current PID;#it{p}_{T} (GeV/#it{c});#LTITS signal#GT", HistType::kTProfile, {itsPtAxis}).get();
      const std::string expectedPath = std::string("DetectorPidQA/PIDQuality/Expected/") + species[i] + "/itsClusterSizeVsP";
      expectedITSResponse[i] = registry.add<TProfile>(expectedPath, "ITS expected mean (model, not measured data)", HistType::kTProfile, {pAxis}).get();
    }
    // Identical axes and response parameters preserve model values on merge.
    // Entries count model samples per worker, not tracks; do not normalize.
    // Sample bin centers; the response is defined only at positive momentum.
    for (int bin = 1; bin <= expectedITSResponse[0]->GetNbinsX(); ++bin) {
      const float p = expectedITSResponse[0]->GetXaxis()->GetBinCenter(bin);
      if (p <= 0.f) {
        continue;
      }
      const std::array<float, 3> expected{
        o2::aod::ITSResponse::expSignal<o2::track::PID::Pion>(p),
        o2::aod::ITSResponse::expSignal<o2::track::PID::Kaon>(p),
        o2::aod::ITSResponse::expSignal<o2::track::PID::Proton>(p)};
      for (std::size_t i = 0; i < expected.size(); ++i) {
        if (std::isfinite(expected[i])) {
          expectedITSResponse[i]->Fill(p, expected[i]);
        }
      }
    }
  }

  template <typename TrackObject>
  void fillPidQualityQA(TrackObject const& track, float cent, int pid)
  {
    // Called only by detectorPidQA, once per quality-selected track. All
    // counts have unit weight; PID decisions come from getPidConfigurable.
    const float pt = track.pt();
    const int nPIDClusters = track.tpcNClsPID();
    registry.fill(HIST("DetectorPidQA/PIDQuality/tpcPIDClustersAll"), pt, nPIDClusters);
    registry.fill(HIST("DetectorPidQA/PIDQuality/tpcPIDClustersVsFound"), track.tpcNClsFound(), nPIDClusters);
    if (pid < 0) {
      registry.fill(HIST("DetectorPidQA/PIDQuality/tpcPIDClustersUnidentified"), pt, nPIDClusters);
    }

    constexpr std::array<int, 3> species{MyParticleType::kPion, MyParticleType::kKaon, MyParticleType::kProton};
    const std::array<float, 3> tpcNSigma{track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr()};
    for (std::size_t i = 0; i < species.size(); ++i) {
      if (!isWithinPOIPtRange(species[i], pt)) {
        continue;
      }
      for (std::size_t stage = 0; stage < 2; ++stage) {
        if (stage == 1 && pid != species[i]) {
          continue;
        }
        auto& h = pidQualityQAHists[i];
        h.tpcPIDClustersVsPt[stage]->Fill(pt, nPIDClusters);
        if (std::isfinite(tpcNSigma[i])) {
          const std::array<double, 4> values{pt, cent, static_cast<double>(nPIDClusters), tpcNSigma[i]};
          h.tpcNSigmaVsPIDClusters[stage]->Fill(values.data());
        }
      }
    }

    // averageClusterSize removes the largest cluster and divides by N-1.
    const uint32_t clusterSizes = track.itsClusterSizes();
    int nSizeLayers = 0;
    for (int layer = 0; layer < 7; ++layer) {
      nSizeLayers += ((clusterSizes >> (4 * layer)) & 0xf) != 0;
    }
    if (!track.hasITS() || nSizeLayers < 2) {
      registry.fill(HIST("DetectorPidQA/PIDQuality/itsResponseStatus"), 0.5);
      return;
    }
    const float momentum = track.p();
    if (!std::isfinite(momentum) || momentum <= 0.f || !std::isfinite(track.eta())) {
      registry.fill(HIST("DetectorPidQA/PIDQuality/itsResponseStatus"), 1.5);
      return;
    }
    const float clusterSize = o2::aod::ITSResponse::averageClusterSize(clusterSizes) / std::cosh(track.eta());
    // Reuse the response used by the low-pT ITS PID; do not change
    // its global parameters or initialize a different data/MC response here.
    const std::array<float, 3> itsNSigma{
      o2::aod::ITSResponse::nSigmaITS<o2::track::PID::Pion>(track),
      o2::aod::ITSResponse::nSigmaITS<o2::track::PID::Kaon>(track),
      o2::aod::ITSResponse::nSigmaITS<o2::track::PID::Proton>(track)};
    const bool allFinite = std::isfinite(clusterSize) && std::all_of(itsNSigma.begin(), itsNSigma.end(), [](float value) { return std::isfinite(value); });
    registry.fill(HIST("DetectorPidQA/PIDQuality/itsResponseStatus"), allFinite ? 3.5 : 2.5);
    if (std::isfinite(clusterSize)) {
      registry.fill(HIST("DetectorPidQA/PIDQuality/itsClusterSizeAll"), momentum, clusterSize);
      registry.fill(HIST("DetectorPidQA/PIDQuality/itsClusterSizeAllVsPt"), pt, clusterSize);
      registry.fill(HIST("DetectorPidQA/PIDQuality/itsClusterSizeVsPEta"), momentum, track.eta(), clusterSize);
    }
    for (std::size_t i = 0; i < species.size(); ++i) {
      if (!isWithinPOIPtRange(species[i], pt)) {
        continue;
      }
      auto& h = pidQualityQAHists[i];
      if (pid == species[i] && std::isfinite(clusterSize)) {
        h.itsClusterSizeVsP->Fill(momentum, clusterSize);
        h.itsClusterSizeVsPt->Fill(pt, clusterSize);
        h.meanITSClusterSizeVsP->Fill(momentum, clusterSize);
        h.meanITSClusterSizeVsPt->Fill(pt, clusterSize);
      }
      if (!std::isfinite(itsNSigma[i])) {
        continue;
      }
      for (std::size_t stage = 0; stage < 2; ++stage) {
        if (stage == 1 && pid != species[i]) {
          continue;
        }
        h.itsNSigmaVsPt[stage]->Fill(pt, itsNSigma[i]);
        if (std::isfinite(tpcNSigma[i])) {
          const std::array<double, 4> values{pt, cent, tpcNSigma[i], itsNSigma[i]};
          h.tpcVsITSNSigma[stage]->Fill(values.data());
        }
      }
    }
  }

  // others

  void initChargedTimeAssociation(InitContext& context)
  {
    if (doprocessData || dofillCorrectionGraph || doprocessReco || doprocessSim || doprocessReco1D || doprocessSim1D ||
        doprocessDataC22DeltaPt || doprocessPidPtCorrelations || doprocessMCClosure || doprocessQA || dodetectorPidQA || doprocessPidResponseMatrix) {
      LOGF(fatal, "Charged time mode: explicitly disable all original process switches to avoid mixed association/PID samples");
    }
    if (!switchsOpts.cfgOutputCharged || switchsOpts.cfgOutputPion || switchsOpts.cfgOutputKaon || switchsOpts.cfgOutputProton || switchsOpts.cfgOutputUnidentified) {
      LOGF(fatal, "Charged time mode requires cfgOutputCharged=true and all four PID/unidentified output switches=false");
    }
    if (dofillCorrectionGraphTimeAssociation && !switchsOpts.cfgOutputrunbyrun) {
      LOGF(fatal, "fillCorrectionGraphTimeAssociation requires switchsOpts.cfgOutputrunbyrun=true");
    }
    if (doprocessMCTimeAssociation && !switchsOpts.cfgOutPutMC) {
      LOGF(fatal, "processMCTimeAssociation requires switchsOpts.cfgOutPutMC=true");
    }
    if (!timeAssociationOption<bool>(context, "pid-flow-charged-time-producer", "processTimeTracks") ||
        timeAssociationOption<bool>(context, "pid-flow-charged-time-producer", "processDummy")) {
      LOGF(fatal, "Charged time mode requires pid-flow-charged-time-producer: processTimeTracks=true, processDummy=false");
    }
  }

  void init(InitContext& context) // Initialization
  {
    const bool timeMode = doprocessDataTimeAssociation || dofillCorrectionGraphTimeAssociation || doprocessMCTimeAssociation || doprocessTimeAssociationComparison;
    if (timeMode) {
      initChargedTimeAssociation(context);
    }
    const AxisSpec activityAxis{axisMultiplicity, getActivityAxisTitle()};
    // FlowContainer's AxisSpec overload expects a full array of bin edges.
    // Use its uniform-axis overload for fixed-width axisMultiplicity settings.
    const auto initializeFlowContainer = [&](auto& container, TObjArray* profiles) {
      if (activityAxis.nBins.has_value()) {
        container->Initialize(profiles, activityAxis.nBins.value(), activityAxis.binEdges.front(), activityAxis.binEdges.back(), cfgFlowNbootstrap);
      } else {
        container->Initialize(profiles, activityAxis, cfgFlowNbootstrap);
      }
      container->GetProfile()->GetXaxis()->SetTitle(getActivityAxisTitle());
      if (auto* subprofiles = container->GetSubProfiles()) {
        for (auto* profile : *subprofiles) {
          static_cast<TProfile2D*>(profile)->GetXaxis()->SetTitle(getActivityAxisTitle());
        }
      }
    };

    // init and add lots of graphs according to switch
    ccdb->setURL(cfgurl.value);
    ccdb->setCaching(true);
    ccdb->setCreatedNotAfter(cfgnolaterthan.value);

    cfgAcceptance = correctionPathOpts.cfgAcceptancePath.value;
    cfgAcceptanceWithPt = correctionPathOpts.cfgAcceptancePathWithPt.value;

    cfgMultPVCutPara = evtSeleOpts.cfgMultPVCut;

    if (!timeMode) {
      fTPCPIDClustersMin = std::make_unique<TF1>("pidFlowPtCorrTPCPIDClustersMin", cfgTPCPIDClustersMinFunction.value.c_str(), 0., 256.);
      if (!fTPCPIDClustersMin->IsValid() || fTPCPIDClustersMin->GetNpar() != 0) {
        LOGF(fatal, "Invalid cfgTPCPIDClustersMinFunction: use a ROOT formula in x = tpcNClsFound() with numeric coefficients and no free parameters");
      }
    }

    // Set the pt, mult and phi Axis;

    if (dcaCutOpts.cfgDCAxyNSigma > 0) {
      dcaCutOpts.fPtDepDCAxy = new TF1("ptDepDCAxy", Form("[0]*%s", dcaCutOpts.cfgDCAxy->c_str()), 0.001, 1000);
      dcaCutOpts.fPtDepDCAxy->SetParameter(0, dcaCutOpts.cfgDCAxyNSigma);
      LOGF(info, "DCAxy pt-dependence function: %s", Form("[0]*%s", dcaCutOpts.cfgDCAxy->c_str()));
    }

    // Add some output objects to the histogram registry
    registry.add("hPhi", "", {HistType::kTH1D, {cfgaxisPhi}});
    registry.add("hPhicorr", "", {HistType::kTH1D, {cfgaxisPhi}});
    registry.add("hPhiPi", "", {HistType::kTH1D, {cfgaxisPhi}});
    registry.add("hPhicorrPi", "", {HistType::kTH1D, {cfgaxisPhi}});
    registry.add("hPhiKa", "", {HistType::kTH1D, {cfgaxisPhi}});
    registry.add("hPhicorrKa", "", {HistType::kTH1D, {cfgaxisPhi}});
    registry.add("hPhiPr", "", {HistType::kTH1D, {cfgaxisPhi}});
    registry.add("hPhicorrPr", "", {HistType::kTH1D, {cfgaxisPhi}});
    registry.add("hPhiCorrNUANUE", "", {HistType::kTH1D, {cfgaxisPhi}});
    registry.add("hPhiCorrNUANUEPi", "", {HistType::kTH1D, {cfgaxisPhi}});
    registry.add("hPhiCorrNUANUEKa", "", {HistType::kTH1D, {cfgaxisPhi}});
    registry.add("hPhiCorrNUANUEPr", "", {HistType::kTH1D, {cfgaxisPhi}});
    registry.add("hEta", "", {HistType::kTH1D, {cfgaxisEta}});
    registry.add("hVtxZ", "", {HistType::kTH1D, {cfgaxisVertex}});
    registry.add("hMult", "", {HistType::kTH1D, {cfgaxisNch}});
    registry.add("hMultTPC", "", {HistType::kTH1D, {cfgaxisNch}});
    std::string hCentTitle = cfgUseMultiplicity.value ? "PV multiplicity distribution" : "Centrality distribution, Estimator " + std::to_string(cfgCentEstimator);
    const AxisSpec hCentAxis = cfgUseMultiplicity.value ? activityAxis : AxisSpec{90, 0., 90., getActivityAxisTitle()};
    registry.add("hCent", hCentTitle.c_str(), {HistType::kTH1D, {hCentAxis}});
    registry.add("hPt", "", {HistType::kTH1D, {cfgaxisPt}});
    registry.add("hPtCorr", "", {HistType::kTH1D, {cfgaxisPt}});
    registry.add("hPtPi", "", {HistType::kTH1D, {cfgaxisPt}});
    registry.add("hPtCorrPi", "", {HistType::kTH1D, {cfgaxisPt}});
    registry.add("hPtKa", "", {HistType::kTH1D, {cfgaxisPt}});
    registry.add("hPtCorrKa", "", {HistType::kTH1D, {cfgaxisPt}});
    registry.add("hPtPr", "", {HistType::kTH1D, {cfgaxisPt}});
    registry.add("hPtCorrPr", "", {HistType::kTH1D, {cfgaxisPt}});
    registry.add("hNTracksPVvsCentrality", "", {HistType::kTH2D, {{5000, 0, 5000}, activityAxis}});

    registry.add("hNchUnCorrectedVSNchCorrected", "", {HistType::kTH2D, {cfgaxisNch, cfgaxisNch}});
    runNumbers = cfgRunNumbers;
    // Minimal detector-PID performance QA.
    if (switchsOpts.cfgOutputQA.value) {
      initPidQualityQA();
      registry.add("DetectorPidQA/TPC/dEdxPi", "TPC d#it{E}/d#it{x} (pion)", {HistType::kTH2F, {cfgaxisPt, cfgaxisTPCDedx}});
      registry.add("DetectorPidQA/TPC/dEdxKa", "TPC d#it{E}/d#it{x} (kaon)", {HistType::kTH2F, {cfgaxisPt, cfgaxisTPCDedx}});
      registry.add("DetectorPidQA/TPC/dEdxPr", "TPC d#it{E}/d#it{x} (proton)", {HistType::kTH2F, {cfgaxisPt, cfgaxisTPCDedx}});
      registry.add("DetectorPidQA/TPC/dEdxAll", "TPC d#it{E}/d#it{x} (all particles)", {HistType::kTH2F, {cfgaxisPt, cfgaxisTPCDedx}});

      registry.add("DetectorPidQA/TOF/deltaPi", "TOF #it{t} - #it{t}_{ev} - #it{t}_{exp} (pion)", {HistType::kTH2F, {cfgaxisPt, cfgaxisTOFDelta}});
      registry.add("DetectorPidQA/TOF/deltaKa", "TOF #it{t} - #it{t}_{ev} - #it{t}_{exp} (kaon)", {HistType::kTH2F, {cfgaxisPt, cfgaxisTOFDelta}});
      registry.add("DetectorPidQA/TOF/deltaPr", "TOF #it{t} - #it{t}_{ev} - #it{t}_{exp} (proton)", {HistType::kTH2F, {cfgaxisPt, cfgaxisTOFDelta}});
      registry.add("DetectorPidQA/TOF/tofTimeAll", "TOF #it{t} - #it{t}_{ev} (all particles)", {HistType::kTH2F, {cfgaxisPt, cfgaxisTOFTime}});
    } // cfgoutputqa

    if (switchsOpts.cfgOutPutMC.value) {
      // hist for NUE correction
      registry.add("correction/hPtCentMcRec", "", {HistType::kTH2D, {cfgaxisPt, activityAxis}});
      registry.add("correction/hPtCentMcGen", "", {HistType::kTH2D, {cfgaxisPt, activityAxis}});
      // ITS only

      // hist for Pi eff
      registry.add("correction/hPtCentMcRecPi", "", {HistType::kTH2D, {cfgaxisPt, activityAxis}});
      registry.add("correction/hPtCentMcGenPi", "", {HistType::kTH2D, {cfgaxisPt, activityAxis}});
      // hist for Ka eff
      registry.add("correction/hPtCentMcRecKa", "", {HistType::kTH2D, {cfgaxisPt, activityAxis}});
      registry.add("correction/hPtCentMcGenKa", "", {HistType::kTH2D, {cfgaxisPt, activityAxis}});
      // hist for Pr eff
      registry.add("correction/hPtCentMcRecPr", "", {HistType::kTH2D, {cfgaxisPt, activityAxis}});
      registry.add("correction/hPtCentMcGenPr", "", {HistType::kTH2D, {cfgaxisPt, activityAxis}});
      // Charged-track reconstruction efficiency split by MC true species.
      // Detector PID is deliberately not used for these numerators.
      registry.add("correction/hPtCentMcRecTruePi", "", {HistType::kTH2D, {cfgaxisPt, activityAxis}});
      registry.add("correction/hPtCentMcRecTrueKa", "", {HistType::kTH2D, {cfgaxisPt, activityAxis}});
      registry.add("correction/hPtCentMcRecTruePr", "", {HistType::kTH2D, {cfgaxisPt, activityAxis}});
      registry.add("correction/hPtCentMcRecTrueOther", "", {HistType::kTH2D, {cfgaxisPt, activityAxis}});
      registry.add("correction/hPtCentMcGenTruePi", "", {HistType::kTH2D, {cfgaxisPt, activityAxis}});
      registry.add("correction/hPtCentMcGenTrueKa", "", {HistType::kTH2D, {cfgaxisPt, activityAxis}});
      registry.add("correction/hPtCentMcGenTruePr", "", {HistType::kTH2D, {cfgaxisPt, activityAxis}});
      registry.add("correction/hPtCentMcGenTrueOther", "", {HistType::kTH2D, {cfgaxisPt, activityAxis}});
    } // cfgoutputMC

    if (switchsOpts.cfgOutPutMC1D.value) {
      registry.add("correction1D/hPtMCRec1D", "", {HistType::kTH1D, {cfgaxisPt}});
      registry.add("correction1D/hPtMCRec1DPi", "", {HistType::kTH1D, {cfgaxisPt}});
      registry.add("correction1D/hPtMCRec1DKa", "", {HistType::kTH1D, {cfgaxisPt}});
      registry.add("correction1D/hPtMCRec1DPr", "", {HistType::kTH1D, {cfgaxisPt}});

      registry.add("correction1D/hPtMCGen1D", "", {HistType::kTH1D, {cfgaxisPt}});
      registry.add("correction1D/hPtMCGen1DPi", "", {HistType::kTH1D, {cfgaxisPt}});
      registry.add("correction1D/hPtMCGen1DKa", "", {HistType::kTH1D, {cfgaxisPt}});
      registry.add("correction1D/hPtMCGen1DPr", "", {HistType::kTH1D, {cfgaxisPt}});
    }

    if (switchsOpts.cfgOutPutPtSpectra.value) {
      registry.add("ptSpectra/hPtCentData", "", {HistType::kTH2D, {cfgaxisPt, activityAxis}});
      registry.add("ptSpectra/hCentEventCountData", "", {HistType::kTH1D, {activityAxis}});
      registry.add("ptSpectra/hCentEventCountMcGen", "", {HistType::kTH1D, {activityAxis}});
      registry.add("ptSpectra/hCentEventCountMcRec", "", {HistType::kTH1D, {activityAxis}});

      registry.add("c22PrimeVsc22/Pi", "", {HistType::kTH2D, {{100, 0., 0.01}, {100, 0., 0.01}}});
      registry.add("c22PrimeVsc22/Ka", "", {HistType::kTH2D, {{100, 0., 0.01}, {100, 0., 0.01}}});
      registry.add("c22PrimeVsc22/Pr", "", {HistType::kTH2D, {{100, 0., 0.01}, {100, 0., 0.01}}});
    } // cfgoutputptspectra

    if (switchsOpts.cfgOutputrunbyrun.value) {
      // hist for NUA
      registry.add("correction/hRunNumberPhiEtaVertex", "", {HistType::kTHnSparseF, {cfgaxisRun, cfgaxisPhi, cfgaxisEta, cfgaxisVertex, cfgaxisPt}});
      registry.add("correction/hRunNumberPhiEtaVertexPion", "", {HistType::kTHnSparseF, {cfgaxisRun, cfgaxisPhi, cfgaxisEta, cfgaxisVertex, cfgaxisPt}});
      registry.add("correction/hRunNumberPhiEtaVertexKaon", "", {HistType::kTHnSparseF, {cfgaxisRun, cfgaxisPhi, cfgaxisEta, cfgaxisVertex, cfgaxisPt}});
      registry.add("correction/hRunNumberPhiEtaVertexProton", "", {HistType::kTHnSparseF, {cfgaxisRun, cfgaxisPhi, cfgaxisEta, cfgaxisVertex, cfgaxisPt}});

      // set "correction/hRunNumberPhiEtaVertex" axis0 label
      for (uint64_t idx = 1; idx <= runNumbers.size(); idx++) {
        registry.get<THnSparse>(HIST("correction/hRunNumberPhiEtaVertex"))->GetAxis(0)->SetBinLabel(idx, std::to_string(runNumbers[idx - 1]).c_str());
        registry.get<THnSparse>(HIST("correction/hRunNumberPhiEtaVertexPion"))->GetAxis(0)->SetBinLabel(idx, std::to_string(runNumbers[idx - 1]).c_str());
        registry.get<THnSparse>(HIST("correction/hRunNumberPhiEtaVertexKaon"))->GetAxis(0)->SetBinLabel(idx, std::to_string(runNumbers[idx - 1]).c_str());
        registry.get<THnSparse>(HIST("correction/hRunNumberPhiEtaVertexProton"))->GetAxis(0)->SetBinLabel(idx, std::to_string(runNumbers[idx - 1]).c_str());
      }
      // end set "correction/hRunNumberPhiEtaVertex" axis0 label
    } // cfgooutputrunbyrun

    if (particleAbundanceOpts.cfgOutPutAbundanceDis) {
      registry.add("abundance/hNumOfPiEventCount", "", {HistType::kTH1D, {particleAbundanceOpts.cfgaxisAbundancePi}});
      registry.add("abundance/hNumOfKaEventCount", "", {HistType::kTH1D, {particleAbundanceOpts.cfgaxisAbundanceKa}});
      registry.add("abundance/hNumOfPrEventCount", "", {HistType::kTH1D, {particleAbundanceOpts.cfgaxisAbundancePr}});
    }
    if (particleAbundanceOpts.cfgOutputAbundanceYields.value) {
      registry.add("abundance/hPtCentDataPi", "", {HistType::kTH2D, {cfgaxisPt, activityAxis}});
      registry.add("abundance/hPtCentDataKa", "", {HistType::kTH2D, {cfgaxisPt, activityAxis}});
      registry.add("abundance/hPtCentDataPr", "", {HistType::kTH2D, {cfgaxisPt, activityAxis}});
      registry.add("abundance/hPtCentDataUnidentified", "", {HistType::kTH2D, {cfgaxisPt, activityAxis}});
    }

    // set bin label for hEventCount
    // processdata
    registry.add("hEventCount/processData", "", {HistType::kTH1D, {{14, 0, 14}}});
    registry.get<TH1>(HIST("hEventCount/processData"))->GetXaxis()->SetBinLabel(1, "Filtered event");
    registry.get<TH1>(HIST("hEventCount/processData"))->GetXaxis()->SetBinLabel(2, "after sel8");
    registry.get<TH1>(HIST("hEventCount/processData"))->GetXaxis()->SetBinLabel(3, "after kTVXinTRD");
    registry.get<TH1>(HIST("hEventCount/processData"))->GetXaxis()->SetBinLabel(4, "after kNoTimeFrameBorder");
    registry.get<TH1>(HIST("hEventCount/processData"))->GetXaxis()->SetBinLabel(5, "after kNoITSROFrameBorder");
    registry.get<TH1>(HIST("hEventCount/processData"))->GetXaxis()->SetBinLabel(6, "after kDoNoSameBunchPileup");
    registry.get<TH1>(HIST("hEventCount/processData"))->GetXaxis()->SetBinLabel(7, "after kIsGoodZvtxFT0vsPV");
    registry.get<TH1>(HIST("hEventCount/processData"))->GetXaxis()->SetBinLabel(8, "after kNoCollInTimeRangeStandard");
    registry.get<TH1>(HIST("hEventCount/processData"))->GetXaxis()->SetBinLabel(9, "after kIsGoodITSLayersAll");
    registry.get<TH1>(HIST("hEventCount/processData"))->GetXaxis()->SetBinLabel(10, "after MultPVCut");
    registry.get<TH1>(HIST("hEventCount/processData"))->GetXaxis()->SetBinLabel(11, "after TPC occupancy cut");
    registry.get<TH1>(HIST("hEventCount/processData"))->GetXaxis()->SetBinLabel(12, "after V0AT0Acut");
    registry.get<TH1>(HIST("hEventCount/processData"))->GetXaxis()->SetBinLabel(13, "after IRmincut");
    registry.get<TH1>(HIST("hEventCount/processData"))->GetXaxis()->SetBinLabel(14, "after IRmaxcut");
    // detector pid QA
    registry.addClone("hEventCount/processData", "hEventCount/detectorPidQA");
    // fillcorrectiongraph
    registry.addClone("hEventCount/processData", "hEventCount/fillCorrectionGraph");
    // processReco
    registry.addClone("hEventCount/processData", "hEventCount/processReco");
    // processSim
    registry.addClone("hEventCount/processData", "hEventCount/processSim");
    // processPidResponseMatrix
    registry.addClone("hEventCount/processData", "hEventCount/processPidResponseMatrix");
    if (timeMode) {
      auto originalCounter = registry.get<TH1>(HIST("hEventCount/processData"));
      for (auto const& name : {"hEventCount/processDataTimeAssociation", "hEventCount/fillCorrectionGraphTimeAssociation", "hEventCount/processMCTimeAssociation"}) {
        auto counter = registry.add<TH1>(name, "Target collisions;selection;collisions", {HistType::kTH1D, {{15, 0, 15}}});
        for (int bin = 1; bin <= originalCounter->GetNbinsX(); ++bin) {
          counter->GetXaxis()->SetBinLabel(bin, originalCounter->GetXaxis()->GetBinLabel(bin));
        }
        counter->GetXaxis()->SetBinLabel(15, "has associated baseline track");
      }
      if (doprocessMCTimeAssociation && cfgTimeAssociationQA) {
        registry.add("timeAssociation/mcEventStatus", "MC event sample;stage;target collisions", {HistType::kTH1D, {{4, 0, 4}}});
        const std::array<const char*, 4> mcLabels{"target collision", "has MC label", "one reco per MC in vertex acceptance", "accepted for reco and gen"};
        for (std::size_t bin = 0; bin < mcLabels.size(); ++bin) {
          registry.get<TH1>(HIST("timeAssociation/mcEventStatus"))->GetXaxis()->SetBinLabel(bin + 1, mcLabels[bin]);
        }
        registry.add("timeAssociation/mcTrackStatus", "Selected charged tracks;truth category;tracks", {HistType::kTH1D, {{4, 0, 4}}});
        const std::array<const char*, 4> mcLabelsTrack{"no MC particle", "other MC collision", "same MC outside primary truth sample", "matched primary numerator"};
        for (std::size_t bin = 0; bin < mcLabelsTrack.size(); ++bin) {
          registry.get<TH1>(HIST("timeAssociation/mcTrackStatus"))->GetXaxis()->SetBinLabel(bin + 1, mcLabelsTrack[bin]);
        }
        registry.add("timeAssociation/hPtCentMcSelected", "All selected reassociated tracks;p_{T};centrality", {HistType::kTH2D, {cfgaxisPt, axisMultiplicity}});
        registry.add("timeAssociation/hPtCentMcOtherCollision", "Tracks from a different true event;p_{T};centrality", {HistType::kTH2D, {cfgaxisPt, axisMultiplicity}});
      }
    }
    if (doprocessTimeAssociationComparison && cfgTimeAssociationQA) {
      registry.addClone("hEventCount/processData", "hEventCount/processTimeAssociationComparison");
      initTimeAssociationComparisonQA();
    }
    if (switchsOpts.cfgAddPidPtCorrelationProfiles.value) {
      registry.addClone("hEventCount/processData", "hEventCount/processPidPtCorrelations");
    }

    if (switchsOpts.cfgAddC22DeltaPtHistograms.value) {
      registry.addClone("hEventCount/processData", "hEventCount/processDataC22DeltaPt");

      // Store the ensemble mean pT required for the offline delta[pT]/<[pT]>
      // transformation. The profile weight is the NUE-corrected particle count.
      if (switchsOpts.cfgOutputCharged.value) {
        registry.add("c22DeltaPt/c22dmeanptCharged", Form("Charged c_{2}{2} versus charged event [p_{T}];%s;[p_{T}]_{ch} (GeV/#it{c});Bootstrap subsample;c_{2}{2}", getActivityAxisTitle()), {HistType::kTProfile3D, {activityAxis, cfgaxisMeanPt, meanptC22GraphOpts.cfgaxisBootstrap}});
        registry.add("c22DeltaPt/hMeanPtCharged", Form("Charged mean p_{T};%s;Bootstrap subsample;[p_{T}]_{ch} (GeV/#it{c})", getActivityAxisTitle()), {HistType::kTProfile2D, {activityAxis, meanptC22GraphOpts.cfgaxisBootstrap}});
      }
      if (switchsOpts.cfgOutputPion.value) {
        registry.add("c22DeltaPt/hMeanPtPion", Form("Pion mean p_{T};%s;Bootstrap subsample;[p_{T}]_{#pi} (GeV/#it{c})", getActivityAxisTitle()), {HistType::kTProfile2D, {activityAxis, meanptC22GraphOpts.cfgaxisBootstrap}});
        registry.add("c22DeltaPt/c22dmeanptPionRefRef", Form("Ref-ref c_{2}{2} versus pion event [p_{T}];%s;[p_{T}]_{#pi} (GeV/#it{c});Bootstrap subsample;c_{2}{2}", getActivityAxisTitle()), {HistType::kTProfile3D, {activityAxis, cfgaxisMeanPt, meanptC22GraphOpts.cfgaxisBootstrap}});
      }
      if (switchsOpts.cfgOutputKaon.value) {
        registry.add("c22DeltaPt/hMeanPtKaon", Form("Kaon mean p_{T};%s;Bootstrap subsample;[p_{T}]_{K} (GeV/#it{c})", getActivityAxisTitle()), {HistType::kTProfile2D, {activityAxis, meanptC22GraphOpts.cfgaxisBootstrap}});
        registry.add("c22DeltaPt/c22dmeanptKaonRefRef", Form("Ref-ref c_{2}{2} versus kaon event [p_{T}];%s;[p_{T}]_{K} (GeV/#it{c});Bootstrap subsample;c_{2}{2}", getActivityAxisTitle()), {HistType::kTProfile3D, {activityAxis, cfgaxisMeanPt, meanptC22GraphOpts.cfgaxisBootstrap}});
      }
      if (switchsOpts.cfgOutputProton.value) {
        registry.add("c22DeltaPt/hMeanPtProton", Form("Proton mean p_{T};%s;Bootstrap subsample;[p_{T}]_{p} (GeV/#it{c})", getActivityAxisTitle()), {HistType::kTProfile2D, {activityAxis, meanptC22GraphOpts.cfgaxisBootstrap}});
        registry.add("c22DeltaPt/c22dmeanptProtonRefRef", Form("Ref-ref c_{2}{2} versus proton event [p_{T}];%s;[p_{T}]_{p} (GeV/#it{c});Bootstrap subsample;c_{2}{2}", getActivityAxisTitle()), {HistType::kTProfile3D, {activityAxis, cfgaxisMeanPt, meanptC22GraphOpts.cfgaxisBootstrap}});
      }
      if (switchsOpts.cfgOutputUnidentified.value) {
        registry.add("c22DeltaPt/hMeanPtUnidentified", Form("Unidentified mean p_{T};%s;Bootstrap subsample;[p_{T}]_{unidentified} (GeV/#it{c})", getActivityAxisTitle()), {HistType::kTProfile2D, {activityAxis, meanptC22GraphOpts.cfgaxisBootstrap}});
        registry.add("c22DeltaPt/c22dmeanptUnidentifiedRefRef", Form("Ref-ref c_{2}{2} versus unidentified event [p_{T}];%s;[p_{T}]_{unidentified} (GeV/#it{c});Bootstrap subsample;c_{2}{2}", getActivityAxisTitle()), {HistType::kTProfile3D, {activityAxis, cfgaxisMeanPt, meanptC22GraphOpts.cfgaxisBootstrap}});
      }

      // Ref-ref is required in either PID mode: Pure/ref-ref or POI-ref/ref-ref.
      // Separate copies are binned by the event mean pT of each PID species.
      // Select one PID correlation method at initialization so disabled-mode
      // profiles are neither allocated nor filled.
      if (switchsOpts.cfgC22DeltaPtUsePure.value) {
        if (switchsOpts.cfgOutputPion.value) {
          registry.add("c22DeltaPt/c22dmeanptPionPure", Form("Pion POI-POI c_{2}{2} versus pion event [p_{T}];%s;[p_{T}]_{#pi} (GeV/#it{c});Bootstrap subsample;c_{2}{2}", getActivityAxisTitle()), {HistType::kTProfile3D, {activityAxis, cfgaxisMeanPt, meanptC22GraphOpts.cfgaxisBootstrap}});
        }
        if (switchsOpts.cfgOutputKaon.value) {
          registry.add("c22DeltaPt/c22dmeanptKaonPure", Form("Kaon POI-POI c_{2}{2} versus kaon event [p_{T}];%s;[p_{T}]_{K} (GeV/#it{c});Bootstrap subsample;c_{2}{2}", getActivityAxisTitle()), {HistType::kTProfile3D, {activityAxis, cfgaxisMeanPt, meanptC22GraphOpts.cfgaxisBootstrap}});
        }
        if (switchsOpts.cfgOutputProton.value) {
          registry.add("c22DeltaPt/c22dmeanptProtonPure", Form("Proton POI-POI c_{2}{2} versus proton event [p_{T}];%s;[p_{T}]_{p} (GeV/#it{c});Bootstrap subsample;c_{2}{2}", getActivityAxisTitle()), {HistType::kTProfile3D, {activityAxis, cfgaxisMeanPt, meanptC22GraphOpts.cfgaxisBootstrap}});
        }
        if (switchsOpts.cfgOutputUnidentified.value) {
          registry.add("c22DeltaPt/c22dmeanptUnidentifiedPure", Form("Unidentified POI-POI c_{2}{2} versus unidentified event [p_{T}];%s;[p_{T}]_{unidentified} (GeV/#it{c});Bootstrap subsample;c_{2}{2}", getActivityAxisTitle()), {HistType::kTProfile3D, {activityAxis, cfgaxisMeanPt, meanptC22GraphOpts.cfgaxisBootstrap}});
        }
      } else {
        // For each PID species, compare POI-ref and ref-ref using the same PID mean-pT axis.
        if (switchsOpts.cfgOutputPion.value) {
          registry.add("c22DeltaPt/c22dmeanptPionPOIRef", Form("Pion POI-ref c_{2}{2} versus pion event [p_{T}];%s;[p_{T}]_{#pi} (GeV/#it{c});Bootstrap subsample;c_{2}{2}", getActivityAxisTitle()), {HistType::kTProfile3D, {activityAxis, cfgaxisMeanPt, meanptC22GraphOpts.cfgaxisBootstrap}});
        }
        if (switchsOpts.cfgOutputKaon.value) {
          registry.add("c22DeltaPt/c22dmeanptKaonPOIRef", Form("Kaon POI-ref c_{2}{2} versus kaon event [p_{T}];%s;[p_{T}]_{K} (GeV/#it{c});Bootstrap subsample;c_{2}{2}", getActivityAxisTitle()), {HistType::kTProfile3D, {activityAxis, cfgaxisMeanPt, meanptC22GraphOpts.cfgaxisBootstrap}});
        }
        if (switchsOpts.cfgOutputProton.value) {
          registry.add("c22DeltaPt/c22dmeanptProtonPOIRef", Form("Proton POI-ref c_{2}{2} versus proton event [p_{T}];%s;[p_{T}]_{p} (GeV/#it{c});Bootstrap subsample;c_{2}{2}", getActivityAxisTitle()), {HistType::kTProfile3D, {activityAxis, cfgaxisMeanPt, meanptC22GraphOpts.cfgaxisBootstrap}});
        }
        if (switchsOpts.cfgOutputUnidentified.value) {
          registry.add("c22DeltaPt/c22dmeanptUnidentifiedPOIRef", Form("Unidentified POI-ref c_{2}{2} versus unidentified event [p_{T}];%s;[p_{T}]_{unidentified} (GeV/#it{c});Bootstrap subsample;c_{2}{2}", getActivityAxisTitle()), {HistType::kTProfile3D, {activityAxis, cfgaxisMeanPt, meanptC22GraphOpts.cfgaxisBootstrap}});
        }
      }
    }

    registry.add("hInteractionRate", "", {HistType::kTH1D, {{1000, 0, 1000}}});
    // end set bin label for eventcount

    if (switchsOpts.cfgAddPidResponseMatrixHistograms.value) {
      registry.add("pidResponseMatrix/hPidResponseMatrix", "PID response matrix;true PID;reconstructed PID", {HistType::kTH2D, {{4, -0.5, 3.5}, {4, -0.5, 3.5}}});
      registry.add("pidResponseMatrix/hPidResponseMatrixRecPtCent", Form("PID response matrix vs reconstructed pT and event activity;true PID;reconstructed PID;p_{T}^{rec} (GeV/#it{c});%s", getActivityAxisTitle()), {HistType::kTHnSparseF, {{4, -0.5, 3.5}, {4, -0.5, 3.5}, cfgaxisPt, activityAxis}});

      auto setPidResponseLabels = [](TAxis* axis, bool includeUnidentified) {
        axis->SetBinLabel(kPidRespPion + 1, "Pion");
        axis->SetBinLabel(kPidRespKaon + 1, "Kaon");
        axis->SetBinLabel(kPidRespProton + 1, "Proton");
        if (includeUnidentified) {
          axis->SetBinLabel(kPidRespUnidentified + 1, "Unidentified");
        } else {
          axis->SetBinLabel(kPidRespOther + 1, "Other");
        }
      };
      setPidResponseLabels(registry.get<TH2>(HIST("pidResponseMatrix/hPidResponseMatrix"))->GetXaxis(), false);
      setPidResponseLabels(registry.get<TH2>(HIST("pidResponseMatrix/hPidResponseMatrix"))->GetYaxis(), true);
      setPidResponseLabels(registry.get<THnSparse>(HIST("pidResponseMatrix/hPidResponseMatrixRecPtCent"))->GetAxis(0), false);
      setPidResponseLabels(registry.get<THnSparse>(HIST("pidResponseMatrix/hPidResponseMatrixRecPtCent"))->GetAxis(1), true);
    }

    // flow container setup
    // cumulant of flow
    // fill TObjArray for charged
    auto* oba4Ch = new TObjArray();
    oba4Ch->Add(new TNamed("c22", "c22"));
    oba4Ch->Add(new TNamed("c32", "c32"));
    oba4Ch->Add(new TNamed("c24", "c24"));
    oba4Ch->Add(new TNamed("c34", "c34"));
    oba4Ch->Add(new TNamed("c22Full", "c22Full"));
    oba4Ch->Add(new TNamed("c22TrackWeight", "c22TrackWeight"));
    oba4Ch->Add(new TNamed("c32TrackWeight", "c32TrackWeight"));
    oba4Ch->Add(new TNamed("c24TrackWeight", "c24TrackWeight"));
    oba4Ch->Add(new TNamed("c34TrackWeight", "c34TrackWeight"));
    oba4Ch->Add(new TNamed("c22FullTrackWeight", "c22FullTrackWeight"));
    oba4Ch->Add(new TNamed("covV2Pt", "covV2Pt"));
    oba4Ch->Add(new TNamed("covV3Pt", "covV3Pt"));
    oba4Ch->Add(new TNamed("ptSquareAve", "ptSquareAve"));
    oba4Ch->Add(new TNamed("ptAve", "ptAve"));
    oba4Ch->Add(new TNamed("hMeanPt", "hMeanPt"));

    oba4Ch->Add(new TNamed("covV2PtWeightOne", "covV2PtWeightOne"));
    oba4Ch->Add(new TNamed("c22TrackWeightOne", "c22TrackWeightOne"));
    oba4Ch->Add(new TNamed("hMeanPtWeightOne", "hMeanPtWeightOne"));
    oba4Ch->Add(new TNamed("ptAveWeightOne", "ptAveWeightOne"));
    oba4Ch->Add(new TNamed("ptSquareAveWeightOne", "ptSquareAveWeightOne"));

    oba4Ch->Add(new TNamed("hMeanPtWeightFull", "hMeanPtWeightFull"));
    // end fill TObjArray for charged

    // Clone the shared profiles before appending the charged-only SC baseline.
    auto* oba4PID = dynamic_cast<TObjArray*>(oba4Ch->Clone());
    oba4Ch->Add(new TNamed("c24Gap", "c24Gap"));

    // init fFCCh
    if (switchsOpts.cfgOutputCharged.value) {
      fFCCh->SetName("FlowContainerCharged");
      initializeFlowContainer(fFCCh, oba4Ch);
    }
    // end init fFCCh

    // init fFCPID
    // note that need to add c22pure and c32pure
    oba4PID->Add(new TNamed("c22pure", "c22pure"));
    oba4PID->Add(new TNamed("c32pure", "c32pure"));
    oba4PID->Add(new TNamed("covV2PtPID", "covV2PtPID"));
    oba4PID->Add(new TNamed("c22TrackWeightPID", "c22TrackWeightPID"));

    oba4PID->Add(new TNamed("hMeanPtWeightCharged", "hMeanPtWeightCharged"));
    // Raw mixed <4> for SC/NSC; reuse PID c22pure and charged c22 as <2>.
    oba4PID->Add(new TNamed("c24PIDChargedGap", "c24PIDChargedGap"));

    if (switchsOpts.cfgOutputPion.value) {
      fFCPi->SetName("FlowContainerPi");
      initializeFlowContainer(fFCPi, oba4PID);
    }
    if (switchsOpts.cfgOutputKaon.value) {
      fFCKa->SetName("FlowContainerKa");
      initializeFlowContainer(fFCKa, oba4PID);
    }
    if (switchsOpts.cfgOutputProton.value) {
      fFCPr->SetName("FlowContainerPr");
      initializeFlowContainer(fFCPr, oba4PID);
    }
    if (switchsOpts.cfgOutputUnidentified.value) {
      fFCUnidentified->SetName("FlowContainerUnidentified");
      initializeFlowContainer(fFCUnidentified, oba4PID);
    }
    // end init fFCPID

    if (switchsOpts.cfgAddPidPtCorrelationProfiles.value) {
      auto* oba4PidPtCorr = new TObjArray();
      oba4PidPtCorr->Add(new TNamed("meanPtPi", "meanPtPi"));
      oba4PidPtCorr->Add(new TNamed("meanPtKa", "meanPtKa"));
      oba4PidPtCorr->Add(new TNamed("meanPtPr", "meanPtPr"));
      oba4PidPtCorr->Add(new TNamed("ptProductPiPi", "ptProductPiPi"));
      oba4PidPtCorr->Add(new TNamed("ptPiInPiPiPairs", "ptPiInPiPiPairs"));
      oba4PidPtCorr->Add(new TNamed("ptProductKaKa", "ptProductKaKa"));
      oba4PidPtCorr->Add(new TNamed("ptKaInKaKaPairs", "ptKaInKaKaPairs"));
      oba4PidPtCorr->Add(new TNamed("ptProductPrPr", "ptProductPrPr"));
      oba4PidPtCorr->Add(new TNamed("ptPrInPrPrPairs", "ptPrInPrPrPairs"));
      oba4PidPtCorr->Add(new TNamed("ptProductPiKa", "ptProductPiKa"));
      oba4PidPtCorr->Add(new TNamed("ptPiInPiKaPairs", "ptPiInPiKaPairs"));
      oba4PidPtCorr->Add(new TNamed("ptKaInPiKaPairs", "ptKaInPiKaPairs"));
      oba4PidPtCorr->Add(new TNamed("ptProductPiPr", "ptProductPiPr"));
      oba4PidPtCorr->Add(new TNamed("ptPiInPiPrPairs", "ptPiInPiPrPairs"));
      oba4PidPtCorr->Add(new TNamed("ptPrInPiPrPairs", "ptPrInPiPrPairs"));
      oba4PidPtCorr->Add(new TNamed("ptProductKaPr", "ptProductKaPr"));
      oba4PidPtCorr->Add(new TNamed("ptKaInKaPrPairs", "ptKaInKaPrPairs"));
      oba4PidPtCorr->Add(new TNamed("ptPrInKaPrPairs", "ptPrInKaPrPairs"));
      fFCPidPtCorr->SetName("FlowContainerPidPtCorr");
      initializeFlowContainer(fFCPidPtCorr, oba4PidPtCorr);
    }

    if (switchsOpts.cfgAddMeanPtCentNbsHistograms.value) {
      // TProfile3D axes: event mean pT, centrality, and bootstrap subsample.
      if (switchsOpts.cfgOutputCharged.value) {
        registry.add("meanptCentNbs/hCharged", "", {HistType::kTProfile3D, {cfgaxisMeanPt, activityAxis, meanptC22GraphOpts.cfgaxisBootstrap}});
        registry.add("meanptCentNbs/hChargedMeanpt", "", {HistType::kTProfile3D, {cfgaxisMeanPt, activityAxis, meanptC22GraphOpts.cfgaxisBootstrap}});
      }
      if (switchsOpts.cfgOutputPion.value) {
        registry.add("meanptCentNbs/hChargedPionWithNpair", "", {HistType::kTProfile3D, {cfgaxisMeanPt, activityAxis, meanptC22GraphOpts.cfgaxisBootstrap}});
        registry.add("meanptCentNbs/hChargedPionFull", "", {HistType::kTProfile3D, {cfgaxisMeanPt, activityAxis, meanptC22GraphOpts.cfgaxisBootstrap}});
        registry.add("meanptCentNbs/hPion", "", {HistType::kTProfile3D, {cfgaxisMeanPt, activityAxis, meanptC22GraphOpts.cfgaxisBootstrap}});
        registry.add("meanptCentNbs/hPionMeanpt", "", {HistType::kTProfile3D, {cfgaxisMeanPt, activityAxis, meanptC22GraphOpts.cfgaxisBootstrap}});
        if (switchsOpts.cfgClosureTest.value != 0) {
          registry.add("meanptCentNbs/hPionMeanptWeightC22pure", "", {HistType::kTProfile3D, {cfgaxisMeanPt, activityAxis, meanptC22GraphOpts.cfgaxisBootstrap}});
          registry.add("meanptCentNbs/hPionMeanptWeightMeanpt", "", {HistType::kTProfile3D, {cfgaxisMeanPt, activityAxis, meanptC22GraphOpts.cfgaxisBootstrap}});
          registry.add("meanptCentNbs/hPionMeanptWeightC22prime", "", {HistType::kTProfile3D, {cfgaxisMeanPt, activityAxis, meanptC22GraphOpts.cfgaxisBootstrap}});
        }
      }
      if (switchsOpts.cfgOutputKaon.value) {
        registry.add("meanptCentNbs/hChargedKaonWithNpair", "", {HistType::kTProfile3D, {cfgaxisMeanPt, activityAxis, meanptC22GraphOpts.cfgaxisBootstrap}});
        registry.add("meanptCentNbs/hChargedKaonFull", "", {HistType::kTProfile3D, {cfgaxisMeanPt, activityAxis, meanptC22GraphOpts.cfgaxisBootstrap}});
        registry.add("meanptCentNbs/hKaon", "", {HistType::kTProfile3D, {cfgaxisMeanPt, activityAxis, meanptC22GraphOpts.cfgaxisBootstrap}});
        registry.add("meanptCentNbs/hKaonMeanpt", "", {HistType::kTProfile3D, {cfgaxisMeanPt, activityAxis, meanptC22GraphOpts.cfgaxisBootstrap}});
      }
      if (switchsOpts.cfgOutputProton.value) {
        registry.add("meanptCentNbs/hChargedProtonWithNpair", "", {HistType::kTProfile3D, {cfgaxisMeanPt, activityAxis, meanptC22GraphOpts.cfgaxisBootstrap}});
        registry.add("meanptCentNbs/hChargedProtonFull", "", {HistType::kTProfile3D, {cfgaxisMeanPt, activityAxis, meanptC22GraphOpts.cfgaxisBootstrap}});
        registry.add("meanptCentNbs/hProton", "", {HistType::kTProfile3D, {cfgaxisMeanPt, activityAxis, meanptC22GraphOpts.cfgaxisBootstrap}});
        registry.add("meanptCentNbs/hProtonMeanpt", "", {HistType::kTProfile3D, {cfgaxisMeanPt, activityAxis, meanptC22GraphOpts.cfgaxisBootstrap}});
      }
      if (switchsOpts.cfgOutputUnidentified.value) {
        registry.add("meanptCentNbs/hChargedUnidentifiedWithNpair", "", {HistType::kTProfile3D, {cfgaxisMeanPt, activityAxis, meanptC22GraphOpts.cfgaxisBootstrap}});
        registry.add("meanptCentNbs/hChargedUnidentifiedFull", "", {HistType::kTProfile3D, {cfgaxisMeanPt, activityAxis, meanptC22GraphOpts.cfgaxisBootstrap}});
        registry.add("meanptCentNbs/hUnidentified", "", {HistType::kTProfile3D, {cfgaxisMeanPt, activityAxis, meanptC22GraphOpts.cfgaxisBootstrap}});
        registry.add("meanptCentNbs/hUnidentifiedMeanpt", "", {HistType::kTProfile3D, {cfgaxisMeanPt, activityAxis, meanptC22GraphOpts.cfgaxisBootstrap}});
      }
    }

    /// @note init QA plot for processQA
    if (switchsOpts.cfgProcessQAOutput.value) {
      /// @note track QA
      registry.add("hProcessQA/hDCAz", "DCAz after collision cuts; DCAz (cm); Pt", {HistType::kTH2D, {{200, -5, 5}, {200, 0, 3}}});
      registry.add("hProcessQA/hDCAxy", "DCAxy after cuts; DCAxy (cm); Pt", {HistType::kTH2D, {{200, -0.5, 0.5}, {200, 0, 5}}});

      registry.add("hProcessQA/hChi2prTPCcls", "#chi^{2}/cluster for the TPC track segment", {HistType::kTH1D, {{100, 0., 5.}}});
      registry.add("hProcessQA/hChi2prITScls", "#chi^{2}/cluster for the ITS track", {HistType::kTH1D, {{100, 0., 50.}}});

      registry.add("hProcessQA/hnTPCClu", "Number of found TPC clusters", {HistType::kTH1D, {{100, 40, 180}}});
      registry.add("hProcessQA/hnITSClu", "Number of found ITS clusters", {HistType::kTH1D, {{10, 0, 10}}});
      registry.add("hProcessQA/hnTPCCrossedRow", "Number of crossed TPC Rows", {HistType::kTH1D, {{100, 40, 180}}});
      // end track QA

      /// @note evetn QA
      registry.add("hProcessQA/centVsMult", Form("Event activity versus input track count;%s;Input track count", getActivityAxisTitle()), {HistType::kTH2D, {activityAxis, cfgaxisNch}});
      registry.add("hProcessQA/IR", "", {HistType::kTH1D, {{100, 0, 100}}});
      registry.add("hProcessQA/Occupacy", "", {HistType::kTH1D, {{1000, 0, 10000}}});
      registry.add("hProcessQA/occupancyVsCent", Form("Track occupancy versus event activity;%s;Track occupancy in time range", getActivityAxisTitle()), {HistType::kTH2D, {activityAxis, {1000, 0, 10000}}});
      // end evetn QA
    }
    // end init QA plot for processQA

    // fgfw set up and correlation config setup
    // Data stored in fGFW
    double etaMax = trkQualityOpts.cfgCutEta.value;
    double etaGap = cfgEtaGap;
    // Charged reference regions.
    fGFW->AddRegion("reffull", -etaMax, etaMax, 1, gfwBitMask::charged::kReference); // ("name", etamin, etamax, ptbinnum, bitmask)eta region -0.8 to 0.8
    fGFW->AddRegion("refN08", -etaMax, -etaGap, 1, gfwBitMask::charged::kReference);
    fGFW->AddRegion("refP08", etaGap, etaMax, 1, gfwBitMask::charged::kReference);
    fGFW->AddRegion("refN", -etaMax, 0, 1, gfwBitMask::charged::kReference);
    fGFW->AddRegion("refP", 0, etaMax, 1, gfwBitMask::charged::kReference);

    // Pion POI and overlap regions.
    fGFW->AddRegion("poiPiN08", -etaMax, -etaGap, 1, gfwBitMask::pion::kPOI);
    fGFW->AddRegion("poiPiP08", etaGap, etaMax, 1, gfwBitMask::pion::kPOI);
    fGFW->AddRegion("poiPiN", -etaMax, 0, 1, gfwBitMask::pion::kPOI);
    fGFW->AddRegion("poiPiP", 0, etaMax, 1, gfwBitMask::pion::kPOI);

    fGFW->AddRegion("olPiN", -etaMax, 0, 1, gfwBitMask::pion::kOverlap);
    fGFW->AddRegion("olPiP", 0, etaMax, 1, gfwBitMask::pion::kOverlap);
    fGFW->AddRegion("olPiN08", -etaMax, -etaGap, 1, gfwBitMask::pion::kOverlap);
    fGFW->AddRegion("olPiP08", etaGap, etaMax, 1, gfwBitMask::pion::kOverlap);

    // Kaon POI and overlap regions.
    fGFW->AddRegion("poiKaN08", -etaMax, -etaGap, 1, gfwBitMask::kaon::kPOI);
    fGFW->AddRegion("poiKaP08", etaGap, etaMax, 1, gfwBitMask::kaon::kPOI);
    fGFW->AddRegion("poiKaN", -etaMax, 0, 1, gfwBitMask::kaon::kPOI);
    fGFW->AddRegion("poiKaP", 0, etaMax, 1, gfwBitMask::kaon::kPOI);

    fGFW->AddRegion("olKaN", -etaMax, 0, 1, gfwBitMask::kaon::kOverlap);
    fGFW->AddRegion("olKaP", 0, etaMax, 1, gfwBitMask::kaon::kOverlap);
    fGFW->AddRegion("olKaN08", -etaMax, -etaGap, 1, gfwBitMask::kaon::kOverlap);
    fGFW->AddRegion("olKaP08", etaGap, etaMax, 1, gfwBitMask::kaon::kOverlap);

    // Proton POI and overlap regions.
    fGFW->AddRegion("poiPrN08", -etaMax, -etaGap, 1, gfwBitMask::proton::kPOI);
    fGFW->AddRegion("poiPrP08", etaGap, etaMax, 1, gfwBitMask::proton::kPOI);
    fGFW->AddRegion("poiPrN", -etaMax, 0, 1, gfwBitMask::proton::kPOI);
    fGFW->AddRegion("poiPrP", 0, etaMax, 1, gfwBitMask::proton::kPOI);

    fGFW->AddRegion("olPrN", -etaMax, 0, 1, gfwBitMask::proton::kOverlap);
    fGFW->AddRegion("olPrP", 0, etaMax, 1, gfwBitMask::proton::kOverlap);
    fGFW->AddRegion("olPrN08", -etaMax, -etaGap, 1, gfwBitMask::proton::kOverlap);
    fGFW->AddRegion("olPrP08", etaGap, etaMax, 1, gfwBitMask::proton::kOverlap);

    // Tracks for which the PID selection returns -1. They use the charged
    // pT range and corrections, but have independent POI and overlap masks.
    fGFW->AddRegion("poiUnidentifiedN08", -etaMax, -etaGap, 1, gfwBitMask::unidentified::kPOI);
    fGFW->AddRegion("poiUnidentifiedP08", etaGap, etaMax, 1, gfwBitMask::unidentified::kPOI);
    fGFW->AddRegion("poiUnidentifiedN", -etaMax, 0, 1, gfwBitMask::unidentified::kPOI);
    fGFW->AddRegion("poiUnidentifiedP", 0, etaMax, 1, gfwBitMask::unidentified::kPOI);
    fGFW->AddRegion("olUnidentifiedN", -etaMax, 0, 1, gfwBitMask::unidentified::kOverlap);
    fGFW->AddRegion("olUnidentifiedP", 0, etaMax, 1, gfwBitMask::unidentified::kOverlap);
    fGFW->AddRegion("olUnidentifiedN08", -etaMax, -etaGap, 1, gfwBitMask::unidentified::kOverlap);
    fGFW->AddRegion("olUnidentifiedP08", etaGap, etaMax, 1, gfwBitMask::unidentified::kOverlap);
    // end data region add

    // pushback
    // Data
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refP08 {2} refN08 {-2}", "Ref08Gap22", false)); // 0
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("reffull {2 2 -2 -2}", "RefFull24", false));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("reffull {2 -2}", "RefFull22", false));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refP08 {3} refN08 {-3}", "Ref08Gap32", false));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refP08 {3 3} refN08 {-3 -3}", "Ref08Gap34", false));

    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiPiN08 {2} refP08 {-2}", "Pion08gap22a", false)); // 5
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiPiP08 {2} refN08 {-2}", "Pion08gap22b", false));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiKaN08 {2} refP08 {-2}", "Kaon08gap22a", false));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiKaP08 {2} refN08 {-2}", "Kaon08gap22b", false));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiPrN08 {2} refP08 {-2}", "Prot08gap22a", false));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiPrP08 {2} refN08 {-2}", "Prot08gap22b", false)); // 10
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiPiN reffull | olPiN {2 2 -2 -2}", "PionFull24a", false));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiPiP reffull | olPiP {2 2 -2 -2}", "PionFull24b", false));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiKaN reffull | olKaN {2 2 -2 -2}", "KaonFull24a", false));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiKaP reffull | olKaP {2 2 -2 -2}", "KaonFull24b", false));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiPrN reffull | olPrN {2 2 -2 -2}", "ProtFull24a", false)); // 15
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiPrP reffull | olPrP {2 2 -2 -2}", "ProtFull24b", false));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiPiN08 {3} refP08 {-3}", "Pion08gap32a", false));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiPiP08 {3} refN08 {-3}", "Pion08gap32b", false));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiKaN08 {3} refP08 {-3}", "Kaon08gap32a", false));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiKaP08 {3} refN08 {-3}", "Kaon08gap32b", false)); // 20
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiPrN08 {3} refP08 {-3}", "Prot08gap32a", false));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiPrP08 {3} refN08 {-3}", "Prot08gap32b", false));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiPiN refN | olPiN {3 3} refP {-3 -3}", "Pion0gap34a", false));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiPiP refP | olPiP {3 3} refN {-3 -3}", "Pion0gap34b", false));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiKaN refN | olKaN {3 3} refP {-3 -3}", "Kaon0gap34a", false)); // 25
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiKaP refP | olKaP {3 3} refN {-3 -3}", "Kaon0gap34b", false));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiPrN refN | olPrN {3 3} refP {-3 -3}", "Prot0gap34a", false));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiPrP refP | olPrP {3 3} refN {-3 -3}", "Prot0gap34b", false));

    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiPiN08 {2} poiPiP08 {-2}", "PiPi08gap22", false));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiKaN08 {2} poiKaP08 {-2}", "KaKa08gap22", false)); // 30
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiPrN08 {2} poiPrP08 {-2}", "PrPr08gap22", false));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiPiN08 {3} poiPiP08 {-3}", "PiPi08gap33", false));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiKaN08 {3} poiKaP08 {-3}", "KaKa08gap33", false));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiPrN08 {3} poiPrP08 {-3}", "PrPr08gap33", false));

    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiPiN reffull | olPiN {2 -2}", "PionFull22a", false)); // 35
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiPiP reffull | olPiP {2 -2}", "PionFull22b", false));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiKaN reffull | olKaN {2 -2}", "KaonFull22a", false));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiKaP reffull | olKaP {2 -2}", "KaonFull22b", false));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiPrN reffull | olPrN {2 -2}", "ProtFull22a", false));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiPrP reffull | olPrP {2 -2}", "ProtFull22b", false)); // 40

    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiUnidentifiedN08 {2} refP08 {-2}", "Unidentified08gap22a", false)); // 41
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiUnidentifiedP08 {2} refN08 {-2}", "Unidentified08gap22b", false));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiUnidentifiedN reffull | olUnidentifiedN {2 2 -2 -2}", "UnidentifiedFull24a", false));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiUnidentifiedP reffull | olUnidentifiedP {2 2 -2 -2}", "UnidentifiedFull24b", false));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiUnidentifiedN08 {2} poiUnidentifiedP08 {-2}", "UnidentifiedPure08gap22", false));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiUnidentifiedN reffull | olUnidentifiedN {2 -2}", "UnidentifiedFull22a", false));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiUnidentifiedP reffull | olUnidentifiedP {2 -2}", "UnidentifiedFull22b", false)); // 47

    // One POI and one reference on each eta side, with distinct tracks.
    // Each overlap region has the same eta range as its POI/reference pair.
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiPiN08 refN08 | olPiN08 {2 2} poiPiP08 refP08 | olPiP08 {-2 -2}", "PionPIDChargedGap24", false)); // 48
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiKaN08 refN08 | olKaN08 {2 2} poiKaP08 refP08 | olKaP08 {-2 -2}", "KaonPIDChargedGap24", false));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiPrN08 refN08 | olPrN08 {2 2} poiPrP08 refP08 | olPrP08 {-2 -2}", "ProtonPIDChargedGap24", false));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiUnidentifiedN08 refN08 | olUnidentifiedN08 {2 2} poiUnidentifiedP08 refP08 | olUnidentifiedP08 {-2 -2}", "UnidentifiedPIDChargedGap24", false)); // 51

    // Charged-charged SC baseline: two distinct reference tracks per eta side.
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN08 {2 2} refP08 {-2 -2}", "Ref08Gap24", false)); // 52

    if (corrconfigs.size() != corrConfigIndex::kCount) {
      LOGF(fatal, "Correlation index table expects %zu configurations, but %zu were created", static_cast<std::size_t>(corrConfigIndex::kCount), corrconfigs.size());
    }
    fGFW->CreateRegions(); // finalize the initialization

    // params
    // used for event selection
    fMultPVCutLow = new TF1("fMultPVCutLow", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x - 3.5*([5]+[6]*x+[7]*x*x+[8]*x*x*x+[9]*x*x*x*x)", 0, 100);
    fMultPVCutLow->SetParameters(cfgMultPVCutPara[0], cfgMultPVCutPara[1], cfgMultPVCutPara[2], cfgMultPVCutPara[3], cfgMultPVCutPara[4], cfgMultPVCutPara[5], cfgMultPVCutPara[6], cfgMultPVCutPara[7], cfgMultPVCutPara[8], cfgMultPVCutPara[9]);
    fMultPVCutHigh = new TF1("fMultPVCutHigh", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x + 3.5*([5]+[6]*x+[7]*x*x+[8]*x*x*x+[9]*x*x*x*x)", 0, 100);
    fMultPVCutHigh->SetParameters(cfgMultPVCutPara[0], cfgMultPVCutPara[1], cfgMultPVCutPara[2], cfgMultPVCutPara[3], cfgMultPVCutPara[4], cfgMultPVCutPara[5], cfgMultPVCutPara[6], cfgMultPVCutPara[7], cfgMultPVCutPara[8], cfgMultPVCutPara[9]);

    fT0AV0AMean = new TF1("fT0AV0AMean", "[0]+[1]*x", 0, 200000);
    fT0AV0AMean->SetParameters(-1601.0581, 9.417652e-01);
    fT0AV0ASigma = new TF1("fT0AV0ASigma", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x", 0, 200000);
    fT0AV0ASigma->SetParameters(463.4144, 6.796509e-02, -9.097136e-07, 7.971088e-12, -2.600581e-17);
  }

  // pid utils function

  template <typename TrackObject>
  int getPidWithCircleCut(TrackObject const& track)
  {
    const float pt = track.pt();

    // If above TOF pt threshold but no TOF signal, reject
    if (pt > circleCutOpts.cfgCircleCutTofPtCut.value && !track.hasTOF()) {
      return -1;
    }

    bool isPionBool = false;
    bool isKaonBool = false;
    bool isProtonBool = false;

    // ------------------------------
    // Pion hypothesis
    // ------------------------------
    {
      const float tpcNsigma = track.tpcNSigmaPi();
      const float tofNsigma = track.tofNSigmaPi();
      const float ptMin = pidPtRangeOpts.cfgPtMin4TOFPiKa.value;
      const float ptMax = pidPtRangeOpts.cfgPtMax4TOFPiKa.value;

      if (pt > ptMin && pt < ptMax) {
        isPionBool = std::hypot(tpcNsigma, tofNsigma) < circleCutOpts.cfgCircleCutSigmaPi.value;
      } else {
        // TPC fallback; the low-pT ITS requirement is applied below.
        isPionBool = std::fabs(tpcNsigma) < circleCutOpts.cfgCircleCutTPCPi.value;
      }
    }

    // ------------------------------
    // Kaon hypothesis (overlap region)
    // ------------------------------
    {
      const float tpcNsigma = track.tpcNSigmaKa();
      const float tofNsigma = track.tofNSigmaKa();
      const float ptMinPiKa = pidPtRangeOpts.cfgPtMin4TOFPiKa.value;
      const float ptMaxPiKa = pidPtRangeOpts.cfgPtMax4TOFPiKa.value;
      const float ptMinKaPr = pidPtRangeOpts.cfgPtMin4TOFKaPr.value;
      const float ptMaxKaPr = pidPtRangeOpts.cfgPtMax4TOFKaPr.value;

      const float ptMin = std::max(ptMinPiKa, ptMinKaPr);
      const float ptMax = std::min(ptMaxPiKa, ptMaxKaPr);

      if (pt > ptMin && pt < ptMax) {
        isKaonBool = std::hypot(tpcNsigma, tofNsigma) < circleCutOpts.cfgCircleCutSigmaKa.value;
      } else {
        isKaonBool = std::fabs(tpcNsigma) < circleCutOpts.cfgCircleCutTPCKa.value;
      }
    }

    // ------------------------------
    // Proton hypothesis
    // ------------------------------
    {
      const float tpcNsigma = track.tpcNSigmaPr();
      const float tofNsigma = track.tofNSigmaPr();
      const float ptMin = pidPtRangeOpts.cfgPtMin4TOFKaPr.value;
      const float ptMax = pidPtRangeOpts.cfgPtMax4TOFKaPr.value;

      if (pt > ptMin && pt < ptMax) {
        isProtonBool = std::hypot(tpcNsigma, tofNsigma) < circleCutOpts.cfgCircleCutSigmaPr.value;
      } else {
        isProtonBool = std::fabs(tpcNsigma) < circleCutOpts.cfgCircleCutTPCPr.value;
      }
    }

    // Apply ITS to every low-pT candidate before resolving ambiguities.
    if (pt < circleCutOpts.cfgCircleCutITSPtMax.value) {
      const uint32_t clusterSizes = track.itsClusterSizes();
      int nSizeLayers = 0;
      for (int layer = 0; layer < 7; ++layer) {
        nSizeLayers += ((clusterSizes >> (4 * layer)) & 0xf) != 0;
      }
      if (!track.hasITS() || nSizeLayers < 2 || !std::isfinite(track.p()) || track.p() <= 0.f || !std::isfinite(track.eta())) {
        return -1;
      }
      const auto passesITS = [](float nsigma, float cut) {
        return std::isfinite(nsigma) && std::fabs(nsigma) < cut;
      };
      isPionBool = isPionBool && passesITS(o2::aod::ITSResponse::nSigmaITS<o2::track::PID::Pion>(track), circleCutOpts.cfgCircleCutITSPi.value);
      isKaonBool = isKaonBool && passesITS(o2::aod::ITSResponse::nSigmaITS<o2::track::PID::Kaon>(track), circleCutOpts.cfgCircleCutITSKa.value);
      isProtonBool = isProtonBool && passesITS(o2::aod::ITSResponse::nSigmaITS<o2::track::PID::Proton>(track), circleCutOpts.cfgCircleCutITSPr.value);
    }

    // ------------------------------
    // Ambiguity rejection
    // ------------------------------
    int nCandidates = isPionBool + isKaonBool + isProtonBool;
    if (nCandidates > 1) {
      return -1; // Reject if multiple hypotheses satisfied
    }

    // ------------------------------
    // Final PID assignment
    // ------------------------------
    if (isPionBool) {
      return MyParticleType::kPion;
    }
    if (isKaonBool) {
      return MyParticleType::kKaon;
    }
    if (isProtonBool) {
      return MyParticleType::kProton;
    }
    return -1;
  }

  /**
   * @brief Apply the TPC PID cluster cut, then circle PID with low-pT ITS
   *
   * @tparam TrackObject
   * @param track
   * @return MyParticleType
   */
  template <typename TrackObject>
  int getPidConfigurable(TrackObject const& track)
  {
    // PID quality changes the species assignment, not charged track selection.
    const double minPIDClusters = fTPCPIDClustersMin->Eval(track.tpcNClsFound());
    if (!std::isfinite(minPIDClusters) || track.tpcNClsPID() < minPIDClusters) {
      return -1;
    }

    return getPidWithCircleCut(track);
  }

  bool isWithinRefPtRange(float pt)
  {
    return (pt > trkQualityOpts.cfgCutPtMin.value) && (pt < trkQualityOpts.cfgCutPtMax.value);
  }

  bool isWithinPOIPtRange(int pid, float pt)
  {
    switch (pid) {
      case MyParticleType::kPion:
        return (pt > trkQualityOpts.cfgCutPtMinPi.value) && (pt < trkQualityOpts.cfgCutPtMaxPi.value);
      case MyParticleType::kKaon:
        return (pt > trkQualityOpts.cfgCutPtMinKa.value) && (pt < trkQualityOpts.cfgCutPtMaxKa.value);
      case MyParticleType::kProton:
        return (pt > trkQualityOpts.cfgCutPtMinPr.value) && (pt < trkQualityOpts.cfgCutPtMaxPr.value);
      default:
        return false;
    }
  }

  int getPidResponseTrueBin(int pdgCode)
  {
    switch (std::abs(pdgCode)) {
      case PDG_t::kPiPlus:
        return kPidRespPion;
      case PDG_t::kKPlus:
        return kPidRespKaon;
      case PDG_t::kProton:
        return kPidRespProton;
      default:
        return kPidRespOther;
    }
  }

  int getPidResponseRecoBin(int recoPid)
  {
    switch (recoPid) {
      case MyParticleType::kPion:
        return kPidRespPion;
      case MyParticleType::kKaon:
        return kPidRespKaon;
      case MyParticleType::kProton:
        return kPidRespProton;
      default:
        return kPidRespUnidentified;
    }
  }

  template <typename McParticle>
  int getTrueSpeciesBin(McParticle const& particle)
  {
    return getPidResponseTrueBin(particle.pdgCode());
  }

  // pid util function

  // other utils

  double getPidC22InOneEvent(const GFW::CorrConfig& corrconfA, const GFW::CorrConfig& corrconfB)
  {
    double NpairA = fGFW->Calculate(corrconfA, 0, true).real();
    double NpairB = fGFW->Calculate(corrconfB, 0, true).real();

    if (NpairA == 0 && NpairB == 0) {
      return 0;
    }

    double ChC22A = NpairA > 0 ? fGFW->Calculate(corrconfA, 0, false).real() / NpairA : 0.;
    double ChC22B = NpairB > 0 ? fGFW->Calculate(corrconfB, 0, false).real() / NpairB : 0.;

    double ChC22 = (ChC22A * NpairA + ChC22B * NpairB) / (NpairA + NpairB);

    return ChC22;
  }

  /**
   * @brief get stable particle
   * @note stable particle include
   *        1. pi
   *        2. Ka
   *        3. proton
   *        4. electron
   *        5. muon
   *
   * @param pdg
   * @return true stable
   * @return false unstable
   */
  bool isStable(int pdg)
  {
    switch (std::abs(pdg)) {
      case PDG_t::kPiPlus:
      case PDG_t::kKPlus:
      case PDG_t::kProton:
      case PDG_t::kElectron:
      case PDG_t::kMuonMinus:
        return true;
        break;
      default:
        return false;
        break;
    }

    return false;
  }

  // other utils

  // fgfw filling helpers

  /**
   * @brief this function is used to fill fFCCh4PtC22
   * @details weight is nch * npair
   *
   * @param cent
   * @param ptSum
   * @param nch
   * @param rndm
   */
  void fillFC4PtC22(const double& cent, const double& ptSum, const double& nch, const double& rndm)
  {
    if (!isParticleOutputEnabled(MyParticleType::kCharged)) {
      return;
    }
    double dnx = 0.0, val = 0.0;

    dnx = fGFW->Calculate(corrconfigs.at(corrConfigIndex::charged::kC22Gap), 0, true).real();
    if (dnx == 0) {
      return;
    }

    // <2>
    val = fGFW->Calculate(corrconfigs.at(corrConfigIndex::charged::kC22Gap), 0, false).real() / dnx;
    if (std::fabs(val) >= 1) {
      return;
    }

    // These large TProfile3D objects are optional; keep FlowContainer output independent of this switch.
    if (switchsOpts.cfgAddMeanPtCentNbsHistograms.value) {
      registry.fill(HIST("meanptCentNbs/hCharged"), ptSum / nch, cent, rndm * cfgFlowNbootstrap, val, nch * dnx);
      registry.fill(HIST("meanptCentNbs/hChargedMeanpt"), ptSum / nch, cent, rndm * cfgFlowNbootstrap, ptSum / nch, 1.);
    }

    fFCCh->FillProfile("hMeanPtWeightFull", cent, (ptSum / nch), nch * dnx, rndm);
  }

  /**
   * @brief note that the graph's x axis is pid meanpt, for <2'> weight is nPid * npairPID, for <2> weight is nPid * npair
   *
   * @param cent
   * @param rndm
   * @param type
   * @param pidPtSum
   * @param nPid
   */
  void fillFC4PtC22(const double& cent, const double& rndm, MyParticleType type, const double& pidPtSum, const double& nPid)
  {
    if (!isParticleOutputEnabled(type)) {
      return;
    }
    // <2>
    double dnx = 0.0, val = 0.0;

    dnx = fGFW->Calculate(corrconfigs.at(corrConfigIndex::charged::kC22Gap), 0, true).real();
    if (dnx == 0) {
      return;
    }

    // <2>
    val = fGFW->Calculate(corrconfigs.at(corrConfigIndex::charged::kC22Gap), 0, false).real() / dnx;
    if (std::fabs(val) >= 1) {
      return;
    }

    double pidc22 = 0;
    double npairPid = 0;
    switch (type) {
      case MyParticleType::kPion:
        pidc22 = getPidC22InOneEvent(corrconfigs.at(corrConfigIndex::pion::kC22GapA), corrconfigs.at(corrConfigIndex::pion::kC22GapB));
        if (pidc22 == 0) {
          return;
        }

        npairPid = fGFW->Calculate(corrconfigs.at(corrConfigIndex::pion::kC22GapA), 0, true).real() + fGFW->Calculate(corrconfigs.at(corrConfigIndex::pion::kC22GapB), 0, true).real();
        if (npairPid == 0) {
          return;
        }

        if (switchsOpts.cfgAddMeanPtCentNbsHistograms.value) {
          registry.fill(HIST("meanptCentNbs/hPion"), pidPtSum / nPid, cent, rndm * cfgFlowNbootstrap, pidc22, nPid * npairPid);
          registry.fill(HIST("meanptCentNbs/hChargedPionFull"), pidPtSum / nPid, cent, rndm * cfgFlowNbootstrap, val, dnx * nPid);
          registry.fill(HIST("meanptCentNbs/hChargedPionWithNpair"), pidPtSum / nPid, cent, rndm * cfgFlowNbootstrap, val, dnx);
          registry.fill(HIST("meanptCentNbs/hPionMeanpt"), pidPtSum / nPid, cent, rndm * cfgFlowNbootstrap, pidPtSum / nPid, 1.);
        }

        fFCPi->FillProfile("hMeanPtWeightFull", cent, (pidPtSum / nPid), nPid * npairPid, rndm);
        fFCPi->FillProfile("hMeanPtWeightCharged", cent, (pidPtSum / nPid), dnx * nPid, rndm);

        if (switchsOpts.cfgAddMeanPtCentNbsHistograms.value && switchsOpts.cfgClosureTest.value != 0) {
          double npair4c22pure = fGFW->Calculate(corrconfigs.at(corrConfigIndex::pion::kC22Pure), 0, true).real();
          if (npair4c22pure > minVal4Float) {
            registry.fill(HIST("meanptCentNbs/hPionMeanptWeightC22pure"),
                          pidPtSum / nPid, cent, rndm * cfgFlowNbootstrap,
                          pidPtSum / nPid,
                          nPid * npairPid * fGFW->Calculate(corrconfigs.at(corrConfigIndex::pion::kC22Pure), 0, false).real() / npair4c22pure);
          }

          registry.fill(HIST("meanptCentNbs/hPionMeanptWeightMeanpt"),
                        pidPtSum / nPid, cent, rndm * cfgFlowNbootstrap,
                        pidPtSum / nPid,
                        nPid * npairPid * pidPtSum / nPid);
          registry.fill(HIST("meanptCentNbs/hPionMeanptWeightC22prime"),
                        pidPtSum / nPid, cent, rndm * cfgFlowNbootstrap,
                        pidPtSum / nPid,
                        nPid * npairPid * pidc22);
        }

        break;
        // end pion

      case MyParticleType::kKaon:
        pidc22 = getPidC22InOneEvent(corrconfigs.at(corrConfigIndex::kaon::kC22GapA), corrconfigs.at(corrConfigIndex::kaon::kC22GapB));
        if (pidc22 == 0) {
          return;
        }

        npairPid = fGFW->Calculate(corrconfigs.at(corrConfigIndex::kaon::kC22GapA), 0, true).real() + fGFW->Calculate(corrconfigs.at(corrConfigIndex::kaon::kC22GapB), 0, true).real();
        if (npairPid == 0) {
          return;
        }

        if (switchsOpts.cfgAddMeanPtCentNbsHistograms.value) {
          registry.fill(HIST("meanptCentNbs/hKaon"), pidPtSum / nPid, cent, rndm * cfgFlowNbootstrap, pidc22, nPid * npairPid);
          registry.fill(HIST("meanptCentNbs/hChargedKaonFull"), pidPtSum / nPid, cent, rndm * cfgFlowNbootstrap, val, dnx * nPid);
          registry.fill(HIST("meanptCentNbs/hChargedKaonWithNpair"), pidPtSum / nPid, cent, rndm * cfgFlowNbootstrap, val, dnx);
          registry.fill(HIST("meanptCentNbs/hKaonMeanpt"), pidPtSum / nPid, cent, rndm * cfgFlowNbootstrap, pidPtSum / nPid, 1.);
        }

        fFCKa->FillProfile("hMeanPtWeightFull", cent, (pidPtSum / nPid), nPid * npairPid, rndm);
        fFCKa->FillProfile("hMeanPtWeightCharged", cent, (pidPtSum / nPid), dnx * nPid, rndm);

        break;
        // end kaon

      case MyParticleType::kProton:
        pidc22 = getPidC22InOneEvent(corrconfigs.at(corrConfigIndex::proton::kC22GapA), corrconfigs.at(corrConfigIndex::proton::kC22GapB));
        if (pidc22 == 0) {
          return;
        }

        npairPid = fGFW->Calculate(corrconfigs.at(corrConfigIndex::proton::kC22GapA), 0, true).real() + fGFW->Calculate(corrconfigs.at(corrConfigIndex::proton::kC22GapB), 0, true).real();
        if (npairPid == 0) {
          return;
        }

        if (switchsOpts.cfgAddMeanPtCentNbsHistograms.value) {
          registry.fill(HIST("meanptCentNbs/hProton"), pidPtSum / nPid, cent, rndm * cfgFlowNbootstrap, pidc22, nPid * npairPid);
          registry.fill(HIST("meanptCentNbs/hChargedProtonFull"), pidPtSum / nPid, cent, rndm * cfgFlowNbootstrap, val, dnx * nPid);
          registry.fill(HIST("meanptCentNbs/hChargedProtonWithNpair"), pidPtSum / nPid, cent, rndm * cfgFlowNbootstrap, val, dnx);
          registry.fill(HIST("meanptCentNbs/hProtonMeanpt"), pidPtSum / nPid, cent, rndm * cfgFlowNbootstrap, pidPtSum / nPid, 1.);
        }

        fFCPr->FillProfile("hMeanPtWeightFull", cent, (pidPtSum / nPid), nPid * npairPid, rndm);
        fFCPr->FillProfile("hMeanPtWeightCharged", cent, (pidPtSum / nPid), dnx * nPid, rndm);

        break;
        // end proton

      case MyParticleType::kUnidentified:
        pidc22 = getPidC22InOneEvent(corrconfigs.at(corrConfigIndex::unidentified::kC22GapA), corrconfigs.at(corrConfigIndex::unidentified::kC22GapB));
        if (pidc22 == 0) {
          return;
        }

        npairPid = fGFW->Calculate(corrconfigs.at(corrConfigIndex::unidentified::kC22GapA), 0, true).real() + fGFW->Calculate(corrconfigs.at(corrConfigIndex::unidentified::kC22GapB), 0, true).real();
        if (npairPid == 0) {
          return;
        }

        if (switchsOpts.cfgAddMeanPtCentNbsHistograms.value) {
          registry.fill(HIST("meanptCentNbs/hUnidentified"), pidPtSum / nPid, cent, rndm * cfgFlowNbootstrap, pidc22, nPid * npairPid);
          registry.fill(HIST("meanptCentNbs/hChargedUnidentifiedFull"), pidPtSum / nPid, cent, rndm * cfgFlowNbootstrap, val, dnx * nPid);
          registry.fill(HIST("meanptCentNbs/hChargedUnidentifiedWithNpair"), pidPtSum / nPid, cent, rndm * cfgFlowNbootstrap, val, dnx);
          registry.fill(HIST("meanptCentNbs/hUnidentifiedMeanpt"), pidPtSum / nPid, cent, rndm * cfgFlowNbootstrap, pidPtSum / nPid, 1.);
        }

        fFCUnidentified->FillProfile("hMeanPtWeightFull", cent, pidPtSum / nPid, nPid * npairPid, rndm);
        fFCUnidentified->FillProfile("hMeanPtWeightCharged", cent, pidPtSum / nPid, dnx * nPid, rndm);
        break;

      default:
        return;
        break;
    }
    // end switch particle type
  }

  bool fillFC(MyParticleType type, const GFW::CorrConfig& corrconf, const double& cent, const double& rndm, const char* tarName)
  {
    if (!isParticleOutputEnabled(type)) {
      return false;
    }
    double dnx = 0.0;
    // calculate #sum exp{i * 0 (#phi_{i} - #phi_{j})} == N_{pairs}
    // note that weight is ignored in the formula but not in the calculation, for c24 is similar
    dnx = fGFW->Calculate(corrconf, 0, true).real();
    if (dnx == 0) {
      return false;
    }
    if (!corrconf.pTDif) {
      // #sum exp{i * 2 * (#phi_{i} - #phi_{j})} / N_{pairs} == < 2 >
      double val = fGFW->Calculate(corrconf, 0, false).real() / dnx;
      if (std::fabs(val) < 1) {
        // NOTE that dnx is WEIGHT
        switch (type) {
          case MyParticleType::kCharged:
            fFCCh->FillProfile(tarName, cent, val, dnx, rndm);
            break;

          case MyParticleType::kPion:
            fFCPi->FillProfile(tarName, cent, val, dnx, rndm);
            break;
          case MyParticleType::kKaon:
            fFCKa->FillProfile(tarName, cent, val, dnx, rndm);
            break;
          case MyParticleType::kProton:
            fFCPr->FillProfile(tarName, cent, val, dnx, rndm);
            break;
          case MyParticleType::kUnidentified:
            fFCUnidentified->FillProfile(tarName, cent, val, dnx, rndm);
            break;

          default:
            LOGF(warning, "particle not found");
            break;
        }
        return true;
      }
    }
    return true;
  }

  /**
   * @brief fill graphs like c22, c24, etc.
   *
   * @tparam chars
   * @param corrconf
   * @param tarName graph name
   * @param cent
   */
  template <char... chars>
  void fillProfile(const GFW::CorrConfig& corrconf, const ConstStr<chars...>& tarName, const double& cent)
  {
    // calculate #sum exp{i * 0 (#phi_{i} - #phi_{j})} == N_{pairs}
    // note that weight is ignored in the formula but not in the calculation, for c24 is similar
    double dnx = fGFW->Calculate(corrconf, 0, true).real();
    if (dnx == 0) {
      return;
    }
    if (!corrconf.pTDif) {
      // #sum exp{i * 2 * (#phi_{i} - #phi_{j})} / N_{pairs} == < 2 >
      double val = fGFW->Calculate(corrconf, 0, false).real() / dnx;
      if (std::fabs(val) < 1) {
        // NOTE that dnx is WEIGHT
        registry.fill(tarName, cent, val, dnx);
        return;
      }
    }
  }

  void fillFCvnpt(MyParticleType type,
                  const GFW::CorrConfig& corrconf,
                  const double& cent,
                  const double& rndm,
                  const double& ptSum,
                  const double& nch,
                  const char* tarName,
                  bool useWeightOne = false)
  {
    if (!isParticleOutputEnabled(type)) {
      return;
    }
    double dnx = fGFW->Calculate(corrconf, 0, true).real();
    if (dnx == 0) {
      return;
    }
    double val = fGFW->Calculate(corrconf, 0, false).real() / dnx;
    if (std::fabs(val) < 1) {
      switch (type) {
        case MyParticleType::kCharged:
          fFCCh->FillProfile(tarName, cent, val * (ptSum / nch), useWeightOne ? 1.0 : dnx * nch, rndm);
          break;

        case MyParticleType::kPion:
          fFCPi->FillProfile(tarName, cent, val * (ptSum / nch), useWeightOne ? 1.0 : dnx * nch, rndm);
          break;
        case MyParticleType::kKaon:
          fFCKa->FillProfile(tarName, cent, val * (ptSum / nch), useWeightOne ? 1.0 : dnx * nch, rndm);
          break;
        case MyParticleType::kProton:
          fFCPr->FillProfile(tarName, cent, val * (ptSum / nch), useWeightOne ? 1.0 : dnx * nch, rndm);
          break;
        case MyParticleType::kUnidentified:
          fFCUnidentified->FillProfile(tarName, cent, val * (ptSum / nch), useWeightOne ? 1.0 : dnx * nch, rndm);
          break;

        default:
          LOGF(warning, "particle not found");
          break;
      }
    }
  }

  /**
   * @brief this function is used to fill weighted profiles
   * @note why we need weightedc22? when calculating cov(v2,pt), we need N_{pair} * N_{charged} to be weight
   *       NOTICE!!! when filling weighted, value gave to param ptSum is nch, when filling diffpt, its just ptsum
   *
   * @tparam chars
   * @param corrconf
   * @param tarName
   * @param cent
   * @param ptSum
   * @param nch
   * @param meanPt
   */
  template <char... chars>
  void fillProfilevnpt(const GFW::CorrConfig& corrconf, const ConstStr<chars...>& tarName, const double& cent, const double& ptSum, const double& nch, const double& meanPt = 0)
  {
    double dnx = 0.0, val = 0.0;
    dnx = fGFW->Calculate(corrconf, 0, true).real();
    if (dnx == 0) {
      return;
    }
    val = fGFW->Calculate(corrconf, 0, false).real() / dnx;
    if (std::fabs(val) < 1) {
      registry.fill(tarName, cent, val * (ptSum / nch - meanPt), dnx * nch);
    }
  }

  template <char... chars>
  void fillProfilePOIvnpt(const GFW::CorrConfig& corrconf, const ConstStr<chars...>& tarName, const double& cent, const double& ptSum, const double& nch, const double& rndm)
  {
    double dnx = 0.0, val = 0.0;
    dnx = fGFW->Calculate(corrconf, 0, true).real();
    if (dnx == 0) {
      return;
    }
    val = fGFW->Calculate(corrconf, 0, false).real() / dnx;

    if (std::fabs(val) < 1) {
      registry.fill(tarName, cent, ptSum / nch, rndm * cfgFlowNbootstrap, val, dnx);
    }
  }

  /**
   * @brief Fill one event-level POI-ref value formed from the two eta orientations.
   * @details The two correlators share one profile entry. Their numerators and pair
   *          counts are added before division. The PID mean pT is only the profile
   *          coordinate, so the profile is weighted by the combined pair count alone.
   */
  template <char... chars>
  void fillProfilePOIvnpt(const GFW::CorrConfig& corrconfA, const GFW::CorrConfig& corrconfB, const ConstStr<chars...>& tarName, const double& cent, const double& ptSum, const double& nPid, const double& rndm)
  {
    const double dnxA = fGFW->Calculate(corrconfA, 0, true).real();
    const double dnxB = fGFW->Calculate(corrconfB, 0, true).real();
    const double dnx = dnxA + dnxB;
    if (dnx == 0.) {
      return;
    }

    const double numerator = fGFW->Calculate(corrconfA, 0, false).real() + fGFW->Calculate(corrconfB, 0, false).real();
    const double val = numerator / dnx;
    if (std::fabs(val) < 1.) {
      registry.fill(tarName, cent, ptSum / nPid, rndm * cfgFlowNbootstrap, val, dnx);
    }
  }

  // fill fgfw helper

  // correction apply functions

  /**
   * @brief load NUE(1D) NUE(2D) NUA graphs
   * @note if u write more than one path in cfg, the graph would not load, that's the strange way to close NUE/NUA corr separatly
   *
   * @param timestamp
   */
  void loadCorrections(uint64_t timestamp)
  {
    if (correctionsLoaded) {
      return;
    }
    int nspecies = 1;
    int totalSpecies = 4;
    /// @note if put 4 path in cfg then load charged, pi, ka, pr
    // else only load charegd
    // if path size != 1 or 4, nothing will be loaded

    // load NUA
    if (cfgAcceptance.size() == static_cast<uint64_t>(nspecies)) {
      for (int i = 0; i < nspecies; i++) {
        mAcceptance.push_back(ccdb->getForTimeStamp<GFWWeights>(cfgAcceptance[i], timestamp));
      }
      if (mAcceptance.size() == static_cast<uint64_t>(nspecies)) {
        LOGF(info, "Loaded acceptance weights");
      } else {
        LOGF(warning, "Could not load acceptance weights");
      }
    }

    if (cfgAcceptance.size() == static_cast<uint64_t>(totalSpecies)) {
      for (int i = 0; i <= totalSpecies - 1; i++) {
        mAcceptance.push_back(ccdb->getForTimeStamp<GFWWeights>(cfgAcceptance[i], timestamp));
      }
      if (mAcceptance.size() == static_cast<uint64_t>(totalSpecies)) {
        LOGF(info, "Loaded acceptance weights * 4");
      } else {
        LOGF(warning, "Could not load acceptance weights * 4");
      }
    }

    // end load NUA

    // load NUA with pt axis (THnSparse: phi, eta, vz, pt), only if cfgUseNUAWithPt is enabled
    // same convention as the GFWWeights NUA above: 1 path -> charged only, 4 paths -> charged, pi, ka, pr
    if (switchsOpts.cfgUseNUAWithPt.value) {
      if (cfgAcceptanceWithPt.size() == static_cast<uint64_t>(nspecies)) {
        for (int i = 0; i < nspecies; i++) {
          mAcceptanceWithPt.push_back(ccdb->getForTimeStamp<THnSparse>(cfgAcceptanceWithPt[i], timestamp));
        }
        if (mAcceptanceWithPt.size() == static_cast<uint64_t>(nspecies)) {
          LOGF(info, "Loaded acceptance weights with pt axis");
        } else {
          LOGF(warning, "Could not load acceptance weights with pt axis");
        }
      }

      if (cfgAcceptanceWithPt.size() == static_cast<uint64_t>(totalSpecies)) {
        for (int i = 0; i <= totalSpecies - 1; i++) {
          mAcceptanceWithPt.push_back(ccdb->getForTimeStamp<THnSparse>(cfgAcceptanceWithPt[i], timestamp));
        }
        if (mAcceptanceWithPt.size() == static_cast<uint64_t>(totalSpecies)) {
          LOGF(info, "Loaded acceptance weights with pt axis * 4");
        } else {
          LOGF(warning, "Could not load acceptance weights with pt axis * 4");
        }
      }
    }
    // end load NUA with pt axis

    /// @note dont modify load NUA!, just modify NUE
    // here needs to load 2 NUE graph selected from 6 ccdb path
    /// @note when using eff 2D, convert TH1* to (TH2D *)
    // pid eff is not available now...
    // Load NUE (PID Efficiency)
    std::vector<std::string> effPath = correctionPathOpts.cfgEfficiencyPath.value;

    if (effPath.size() == static_cast<uint64_t>(nspecies)) {
      for (int i = 0; i < nspecies; i++) {
        mEfficiency.push_back(ccdb->getForTimeStamp<TH2>(effPath[i], timestamp));
      }
      if (mEfficiency.size() == static_cast<uint64_t>(nspecies)) {
        LOGF(info, "Loaded PID efficiency histogram");
      } else {
        LOGF(warning, "Could not load PID efficiency histogram");
      }
    }

    if (effPath.size() == static_cast<uint64_t>(totalSpecies)) {
      for (int i = 0; i <= totalSpecies - 1; i++) {
        mEfficiency.push_back(ccdb->getForTimeStamp<TH2>(effPath[i], timestamp));
      }
      if (mEfficiency.size() == static_cast<uint64_t>(totalSpecies)) {
        LOGF(info, "Loaded PID efficiency histogram * 4");
      } else {
        LOGF(warning, "Could not load PID efficiency histogram * 4");
      }
    }

    // End load PID efficiency for ITS + global track

    correctionsLoaded = true;
  }

  /**
   * @brief get NUA weight from a THnSparse(phi, eta, vz, pt), same convention as GFWWeights::getNUA
   *        (bin content is the normalized acceptance, 0-1, so the correction weight is 1/content)
   *
   * @param hist THnSparse with axes (phi, eta, vz, pt), may be nullptr
   */
  float getNUAFromTHn(THnSparse* hist, float phi, float eta, float vz, float pt)
  {
    if (!hist) {
      return 1;
    }
    std::array<int, 4> bins{};
    bins[0] = hist->GetAxis(0)->FindBin(phi);
    bins[1] = hist->GetAxis(1)->FindBin(eta);
    bins[2] = hist->GetAxis(2)->FindBin(vz);
    bins[3] = hist->GetAxis(3)->FindBin(pt);
    double weight = hist->GetBinContent(bins.data());
    if (weight != 0) {
      return 1. / weight;
    }
    return 1;
  }

  /**
   * @brief Set the Particle NUA, note that it works only when macceptance size >= 1,
   * macceptance size can only be 1 or 4
   * @note controlled by switchsOpts.cfgUseNUAWithPt: false uses the original GFWWeights
   *       (phi,eta,vz) NUA (mAcceptance); true uses the THnSparse (phi,eta,vz,pt) NUA (mAcceptanceWithPt)
   *
   * @tparam TrackObject
   * @param weight_nua
   * @param track
   * @param vtxz
   * @return true
   * @return false
   */
  template <typename TrackObject>
  bool setParticleNUAWeight(float& weight_nua, TrackObject& track, float vtxz)
  {
    if (switchsOpts.cfgUseNUAWithPt.value) {
      if (mAcceptanceWithPt.size() == static_cast<uint64_t>(1) || mAcceptanceWithPt.size() == static_cast<uint64_t>(4)) {
        weight_nua = getNUAFromTHn(mAcceptanceWithPt[0], track.phi(), track.eta(), vtxz, track.pt());
      } else {
        weight_nua = 1;
      }
      return true;
    }
    if (mAcceptance.size() == static_cast<uint64_t>(1) || mAcceptance.size() == static_cast<uint64_t>(4)) {
      weight_nua = mAcceptance[0]->getNUA(track.phi(), track.eta(), vtxz);
    } else {
      weight_nua = 1;
    }
    return true;
  }

  template <typename TrackObject>
  bool setParticleNUAWeight(float& weight_nua, TrackObject& track, float vtxz, int pid)
  {
    if (switchsOpts.cfgUseNUAWithPt.value) {
      if (mAcceptanceWithPt.size() == static_cast<uint64_t>(4)) {
        switch (pid) {
          case MyParticleType::kPion:
            weight_nua = getNUAFromTHn(mAcceptanceWithPt[1], track.phi(), track.eta(), vtxz, track.pt());
            break;
          case MyParticleType::kKaon:
            weight_nua = getNUAFromTHn(mAcceptanceWithPt[2], track.phi(), track.eta(), vtxz, track.pt());
            break;
          case MyParticleType::kProton:
            weight_nua = getNUAFromTHn(mAcceptanceWithPt[3], track.phi(), track.eta(), vtxz, track.pt());
            break;

          default:
            weight_nua = 1;
            break;
        }
      } else {
        weight_nua = 1;
      }
      return true;
    }
    if (mAcceptance.size() == static_cast<uint64_t>(4)) {

      switch (pid) {
        case MyParticleType::kPion:
          weight_nua = mAcceptance[1]->getNUA(track.phi(), track.eta(), vtxz);
          break;
        case MyParticleType::kKaon:
          weight_nua = mAcceptance[2]->getNUA(track.phi(), track.eta(), vtxz);
          break;
        case MyParticleType::kProton:
          weight_nua = mAcceptance[3]->getNUA(track.phi(), track.eta(), vtxz);
          break;

        default:
          weight_nua = 1;
          break;
      }

    } else {
      weight_nua = 1;
    }
    return true;
  }

  template <typename TrackObject>
  bool setParticleNUEWeight(float& weight_nue, TrackObject& track, double cent)
  {
    if (mEfficiency.size() == static_cast<uint64_t>(1) || mEfficiency.size() == static_cast<uint64_t>(4)) {
      TH2* eff2D = mEfficiency[0];
      int ptBin = eff2D->GetXaxis()->FindBin(track.pt());
      int centBin = eff2D->GetYaxis()->FindBin(cent);
      float eff = eff2D->GetBinContent(ptBin, centBin);

      if (eff == 0.) {
        weight_nue = 1.;
        return true;
      }
      weight_nue = 1. / eff;
    } else {
      weight_nue = 1.;
    }
    return true;
  }

  template <typename TrackObject>
  bool setParticleNUEWeight(float& weight_nue, TrackObject& track, double cent, int pid)
  {
    if (mEfficiency.size() == static_cast<uint64_t>(4)) {
      TH2* eff2D = nullptr;
      switch (pid) {
        case MyParticleType::kPion:
          eff2D = mEfficiency[1];
          break;
        case MyParticleType::kKaon:
          eff2D = mEfficiency[2];
          break;
        case MyParticleType::kProton:
          eff2D = mEfficiency[3];
          break;
        default:
          weight_nue = 1.;
          return true;
      }

      int ptBin = eff2D->GetXaxis()->FindBin(track.pt());
      int centBin = eff2D->GetYaxis()->FindBin(cent);
      float eff = eff2D->GetBinContent(ptBin, centBin);

      if (eff == 0.) {
        weight_nue = 1.;
        return true;
      }
      weight_nue = 1. / eff;
    } else {
      weight_nue = 1.;
    }
    return true;
  }

  // correction apply function

  // track cut functions

  /**
   * @brief cut MC particles
   * @note include
   *        1. eta cut
   *        2. pt cut
   *        3. stable partccle
   *
   * @tparam mcParticle
   * @param particle
   * @return true
   * @return false
   */
  template <typename mcParticle>
  bool particleSelected(mcParticle& particle)
  {
    // eta
    if (std::fabs(particle.eta()) > trkQualityOpts.cfgCutEta.value) {
      return false;
    }
    // pt
    if (particle.pt() < trkQualityOpts.cfgCutPtMin.value) {
      return false;
    }
    if (particle.pt() > trkQualityOpts.cfgCutPtMax.value) {
      return false;
    }
    // stable particle
    if (!isStable(particle.pdgCode())) {
      return false;
    }

    return true;
  }

  /**
   * @brief global track selection, cut for specified detector is not include
   * @note include: 1. dcaxy
   *                2. pt cut
   *                3. eta
   *                4. dcaz
   *
   *
   * @tparam TTrack
   * @param track
   * @param checkPt if false, skip the internal ref-pt cut (cfgCutPtMin/cfgCutPtMax); needed
   *                so callers can apply an independent POI pt range instead (e.g. processData)
   * @return true pass selection
   * @return false reject track
   */
  template <typename TTrack>
  bool trackSelectedGlobal(TTrack& track, bool checkPt = true)
  {
    // dca xy cut
    // note that z cut is in filter
    if (dcaCutOpts.cfgDCAxyNSigma && (std::fabs(track.dcaXY()) > dcaCutOpts.fPtDepDCAxy->Eval(track.pt()))) {
      return false;
    }

    // pt
    if (checkPt && !((track.pt() > trkQualityOpts.cfgCutPtMin.value) && (track.pt() < trkQualityOpts.cfgCutPtMax.value))) {
      return false;
    }

    // eta
    if (std::fabs(track.eta()) > trkQualityOpts.cfgCutEta.value) {
      return false;
    }

    // dca z
    if (std::fabs(track.dcaZ()) > dcaCutOpts.cfgCutDCAz.value) {
      return false;
    }

    return true;
  }

  /**
   * @brief cut for ITS, include:
   *        1. its n cluster
   *        2. its chi2
   *
   * @tparam TTrack
   * @param track
   * @return true
   * @return false
   */
  template <typename TTrack>
  bool trackSelected4ITS(TTrack& track)
  {
    // its ncls cut
    if (track.itsNCls() <= trkQualityOpts.cfgITSNCls.value) {
      return false;
    }

    // its chi2
    if (track.itsChi2NCl() > trkQualityOpts.cfgITSChi2NDF.value) {
      return false;
    }

    return true;
  }

  /**
   * @brief cut for TPC, include
   *        1. TPC  cross rows
   *        2. tpc n cluster
   *        3. tpc chi2
   *
   * @tparam TTrack
   * @param track
   * @return true
   * @return false
   */
  template <typename TTrack>
  bool trackSelected4TPC(TTrack& track)
  {
    // tpc crossedrows cut
    if (track.tpcNClsCrossedRows() <= trkQualityOpts.cfgTPCCrossedRows.value) {
      return false;
    }
    // tpc ncls cut
    if (track.tpcNClsFound() <= trkQualityOpts.cfgTPCNCls.value) {
      return false;
    }
    // tpc chi2
    if (track.tpcChi2NCl() > cfgCutChi2prTPCcls) {
      return false;
    }

    return true;
  }

  template <typename TTrack>
  bool trackSelectedForFlow(TTrack& track)
  {
    return trackSelectedGlobal(track, false) && track.hasITS() && track.hasTPC() && trackSelected4ITS(track) && trackSelected4TPC(track);
  }

  // track cut

  // event selection functions

  /**
   * @brief fill eventCount for different function
   *
   * @param funcName
   * @param position
   */
  void fillEventCountHelper(MyFunctionName funcName, double position)
  {
    switch (funcName) {
      case MyFunctionName::funcProcessData:
        registry.fill(HIST("hEventCount/processData"), position);
        break;
      case MyFunctionName::funcDetectorPidQA:
        registry.fill(HIST("hEventCount/detectorPidQA"), position);
        break;
      case MyFunctionName::funcFillCorrectionGraph:
        registry.fill(HIST("hEventCount/fillCorrectionGraph"), position);
        break;
      case MyFunctionName::funcProcessReco:
        registry.fill(HIST("hEventCount/processReco"), position);
        break;
      case MyFunctionName::funcProcessSim:
        registry.fill(HIST("hEventCount/processSim"), position);
        break;
      case MyFunctionName::funcProcessPidResponseMatrix:
        registry.fill(HIST("hEventCount/processPidResponseMatrix"), position);
        break;
      case MyFunctionName::funcProcessDataC22DeltaPt:
        registry.fill(HIST("hEventCount/processDataC22DeltaPt"), position);
        break;
      case MyFunctionName::funcProcessPidPtCorrelations:
        registry.fill(HIST("hEventCount/processPidPtCorrelations"), position);
        break;
      case MyFunctionName::funcProcessDataTimeAssociation:
        registry.fill(HIST("hEventCount/processDataTimeAssociation"), position);
        break;
      case MyFunctionName::funcFillCorrectionGraphTimeAssociation:
        registry.fill(HIST("hEventCount/fillCorrectionGraphTimeAssociation"), position);
        break;
      case MyFunctionName::funcProcessMCTimeAssociation:
        registry.fill(HIST("hEventCount/processMCTimeAssociation"), position);
        break;
      case MyFunctionName::funcProcessTimeAssociationComparison:
        registry.fill(HIST("hEventCount/processTimeAssociationComparison"), position);
        break;

      default:
        // LOGF(warning, "could not find event count graph");
        break;
    }
  }

  /**
   * @brief Get the event activity used for output binning and NUE lookup.
   * @note Multiplicity mode uses raw PV contributors in |eta| < 0.8.
   *       It requires matching multiplicity axes and NUE maps; the caller
   *       must disable the centrality-calibrated cfgDoMultPVCut.
   */
  template <typename TCollision>
  float getCentrality(TCollision const& collision)
  {
    if (cfgUseMultiplicity.value) {
      return collision.multNTracksPV();
    }
    return getPercentileCentrality(collision);
  }

  template <typename TCollision>
  float getPercentileCentrality(TCollision const& collision)
  {
    float cent = -1;
    switch (cfgCentEstimator) {
      case kCentFT0C:
        cent = collision.centFT0C();
        break;
      case kCentFT0CVariant1:
        cent = collision.centFT0CVariant1();
        break;
      case kCentFT0M:
        cent = collision.centFT0M();
        break;
      case kCentFV0A:
        cent = collision.centFV0A();
        break;
      case kCentNTPV:
        cent = collision.centNTPV();
        break;
      case kCentNGlobal:
        cent = collision.centNGlobal();
        break;
      case kCentMFT:
        cent = collision.centMFT();
        break;
      default:
        cent = collision.centFT0C();
    }
    return cent;
  }

  /**
   * @brief collision selection
   * @note include: 1. TRD triggered
   *                2. collisions close to Time Frame borders
   *                3. ITS ROF border
   *                4. pile up
   *                5. mathc z of PV by tracks and z of PV from FT0 A-C
   *                6. no collisions in specified time range
   *                7. good ITS layer
   *                8. vz cut && multiplicity cut
   *                9. occupacy cut
   *                10. V0A T0A 5 sigma cut
   *                11. interaction rate
   *
   * @tparam TCollision
   * @param collision
   * @param centrality
   * @param interactionRate
   * @param functionName need to fill graph in different function, dont fill the same graph!
   * @return true
   * @return false
   */
  template <typename TCollision>
  bool eventSelected(TCollision& collision, const float centrality, float interactionRate = -1, MyFunctionName funcName = MyFunctionName::funcProcessData)
  {
    if (evtSeleOpts.cfgDoTVXinTRD.value && collision.alias_bit(kTVXinTRD)) {
      // TRD triggered
      return false;
    }
    fillEventCountHelper(funcName, 2.5);
    if (evtSeleOpts.cfgDoNoTimeFrameBorder.value && !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
      // reject collisions close to Time Frame borders
      // https://its.cern.ch/jira/browse/O2-4623
      return false;
    }
    fillEventCountHelper(funcName, 3.5);
    if (evtSeleOpts.cfgDoNoITSROFrameBorder.value && !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      // reject events affected by the ITS ROF border
      // https://its.cern.ch/jira/browse/O2-4309
      return false;
    }
    fillEventCountHelper(funcName, 4.5);
    if (evtSeleOpts.cfgDoNoSameBunchPileup.value && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      // rejects collisions which are associated with the same "found-by-T0" bunch crossing
      // https://indico.cern.ch/event/1396220/#1-event-selection-with-its-rof
      return false;
    }
    fillEventCountHelper(funcName, 5.5);
    if (evtSeleOpts.cfgDoIsGoodZvtxFT0vsPV.value && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      // removes collisions with large differences between z of PV by tracks and z of PV from FT0 A-C time difference
      // use this cut at low multiplicities with caution
      return false;
    }
    fillEventCountHelper(funcName, 6.5);
    if (evtSeleOpts.cfgDoNoCollInTimeRangeStandard.value && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      // no collisions in specified time range
      return 0;
    }
    fillEventCountHelper(funcName, 7.5);
    if (evtSeleOpts.cfgDoIsGoodITSLayersAll.value && !collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
      // cut time intervals with dead ITS staves
      return 0;
    }
    fillEventCountHelper(funcName, 8.5);
    float vtxz = -999;
    if (collision.numContrib() > 1) {
      vtxz = collision.posZ();
      float zRes = std::sqrt(collision.covZZ());
      double zResMin = 0.25;
      int numContMax = 20;
      if (zRes > zResMin && collision.numContrib() < numContMax) {
        vtxz = -999;
      }
    }
    auto multNTracksPV = collision.multNTracksPV();
    auto occupancy = collision.trackOccupancyInTimeRange();

    if (std::fabs(vtxz) > cfgCutVertex) {
      return false;
    }

    // The comparison has its own counter; do not double-fill physics event QA.
    if (funcName != funcProcessTimeAssociationComparison) {
      registry.fill(HIST("hNTracksPVvsCentrality"), multNTracksPV, centrality);
    }
    if (evtSeleOpts.cfgDoMultPVCut.value) {
      if (multNTracksPV < fMultPVCutLow->Eval(centrality)) {
        return false;
      }
      if (multNTracksPV > fMultPVCutHigh->Eval(centrality)) {
        return false;
      }
    }
    fillEventCountHelper(funcName, 9.5);

    if (occupancy > evtSeleOpts.cfgCutOccupancyHigh.value) {
      return 0;
    }
    fillEventCountHelper(funcName, 10.5);

    // V0A T0A 5 sigma cut
    if (evtSeleOpts.cfgDoV0AT0Acut.value) {
      int nsigma = 5;
      if (std::fabs(collision.multFV0A() - fT0AV0AMean->Eval(collision.multFT0A())) > nsigma * fT0AV0ASigma->Eval(collision.multFT0A())) {
        return 0;
      }
    }
    fillEventCountHelper(funcName, 11.5);

    if (funcName != funcProcessTimeAssociationComparison) {
      registry.fill(HIST("hInteractionRate"), interactionRate);
    }
    if (evtSeleOpts.cfgCutminIR.value >= 0 && interactionRate > 0 && interactionRate < evtSeleOpts.cfgCutminIR.value) {
      return false;
    }
    fillEventCountHelper(funcName, 12.5);
    if (evtSeleOpts.cfgCutmaxIR.value >= 0 && interactionRate > evtSeleOpts.cfgCutmaxIR.value) {
      return false;
    }
    fillEventCountHelper(funcName, 13.5);

    return true;
  }

  double getInteractionRate(uint64_t timestamp, int runNumber)
  {
    const bool useMinIRCut = evtSeleOpts.cfgCutminIR.value >= 0;
    const bool useMaxIRCut = evtSeleOpts.cfgCutmaxIR.value >= 0;
    if ((!useMinIRCut && !useMaxIRCut) || evtSeleOpts.cfgIRSource.value.empty()) {
      return -1.;
    }
    return rateFetcher.fetch(ccdb.service, timestamp, runNumber, evtSeleOpts.cfgIRSource.value) * 1.e-3;
  }

  // event selection

  // main functions

  void loadChargedTimeCorrections(uint64_t timestamp, bool acceptance = true)
  {
    // Preserve the existing one/four-path convention, but load only charged.
    constexpr std::size_t FullCorrectionSpecies = 4;
    if (!timeEfficiencyLoaded) {
      auto const& paths = correctionPathOpts.cfgEfficiencyPath.value;
      if (paths.size() == 1 || paths.size() == FullCorrectionSpecies) {
        auto* efficiency = ccdb->getForTimeStamp<TH2>(paths.front(), timestamp);
        if (!efficiency) {
          LOGF(fatal, "Missing charged time-association efficiency: %s", paths.front().c_str());
        }
        mEfficiency.push_back(efficiency);
      }
      timeEfficiencyLoaded = true;
    }
    if (acceptance && !timeAcceptanceLoaded) {
      if (switchsOpts.cfgUseNUAWithPt) {
        if (cfgAcceptanceWithPt.size() == 1 || cfgAcceptanceWithPt.size() == FullCorrectionSpecies) {
          auto* weights = ccdb->getForTimeStamp<THnSparse>(cfgAcceptanceWithPt.front(), timestamp);
          if (!weights) {
            LOGF(fatal, "Missing charged time-association NUA: %s", cfgAcceptanceWithPt.front().c_str());
          }
          mAcceptanceWithPt.push_back(weights);
        }
      } else if (cfgAcceptance.size() == 1 || cfgAcceptance.size() == FullCorrectionSpecies) {
        auto* weights = ccdb->getForTimeStamp<GFWWeights>(cfgAcceptance.front(), timestamp);
        if (!weights) {
          LOGF(fatal, "Missing charged time-association NUA: %s", cfgAcceptance.front().c_str());
        }
        mAcceptance.push_back(weights);
      }
      timeAcceptanceLoaded = true;
    }
  }

  template <typename TCollision>
  bool chargedTimeEventSelected(TCollision const& collision, MyFunctionName name)
  {
    fillEventCountHelper(name, 0.5);
    if (!collision.sel8()) {
      return false;
    }
    fillEventCountHelper(name, 1.5);
    auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
    return eventSelected(collision, getCentrality(collision), getInteractionRate(bc.timestamp(), bc.runNumber()), name);
  }

  void processData(AodCollisions::iterator const& collision, aod::BCsWithTimestamps const&, AodTracks const& tracks)
  {
    // sub region, init
    // init
    float rndm = fRndm->Rndm();
    int nTot = tracks.size();
    float nMultTPC = collision.multTPC();
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    int runNumber = bc.runNumber();
    double interactionRate = getInteractionRate(bc.timestamp(), runNumber);
    // end init

    // collision cut
    // include : 1.track.size 2.collision.sel8 3. evenSelected
    registry.fill(HIST("hEventCount/processData"), 0.5);
    if (nTot < 1) {
      return;
    }
    fGFW->Clear();
    const auto cent = getCentrality(collision);
    if (!collision.sel8()) {
      return;
    }
    registry.fill(HIST("hEventCount/processData"), 1.5);
    if (!eventSelected(collision, cent, interactionRate)) {
      return;
    }
    // end collision cut

    // pt spectra
    if (switchsOpts.cfgOutPutPtSpectra.value) {
      registry.fill(HIST("ptSpectra/hCentEventCountData"), cent);
    } // switchsOpts.cfgOutPutPtSpectra.value

    // correction
    loadCorrections(bc.timestamp());
    float vtxz = collision.posZ();
    registry.fill(HIST("hVtxZ"), vtxz);
    registry.fill(HIST("hMult"), nTot);
    registry.fill(HIST("hMultTPC"), nMultTPC);
    registry.fill(HIST("hCent"), cent);
    // end correction

    double ptSum = 0;
    double ptSumw2 = 0;
    double nch = 0;
    double nchSquare = 0;
    double ptSquareSum = 0;
    double nchCorrectedPassed = 0;
    double nchPassedSelection = 0;

    if (switchsOpts.cfgDebugMyCode.value) {
      LOGF(info, "====================================");
    }

    // val for pid particles
    double pionPtSum = 0;
    double kaonPtSum = 0;
    double protonPtSum = 0;
    double unidentifiedPtSum = 0;

    double nPionWeighted = 0;
    double nKaonWeighted = 0;
    double nProtonWeighted = 0;
    double nUnidentifiedWeighted = 0;

    double pionPtSumw2 = 0;
    double kaonPtSumw2 = 0;
    double protonPtSumw2 = 0;
    double unidentifiedPtSumw2 = 0;

    double nPionSquare = 0;
    double nKaonSquare = 0;
    double nProtonSquare = 0;
    double nUnidentifiedSquare = 0;

    double pionPtSquareSum = 0;
    double kaonPtSquareSum = 0;
    double protonPtSquareSum = 0;
    double unidentifiedPtSquareSum = 0;
    // end val for pid particles

    // Calculate mean-pT moments and PID-particle sums.
    for (const auto& track : tracks) {
      float weff = 1;

      // The reference and POI pT ranges are applied independently below.
      if (!trackSelectedForFlow(track)) {
        continue;
      }

      // do nue
      setParticleNUEWeight(weff, track, cent);
      // end do nue

      // calculate ncharged(nch with weight) and pt
      bool withinPtRef = isWithinRefPtRange(track.pt());
      if (withinPtRef && std::fabs(track.eta()) < trkQualityOpts.cfgRangeEta.value) {
        nch += weff;
        nchSquare += weff * weff;
        ptSum += weff * track.pt();
        ptSumw2 += weff * weff * track.pt();
        ptSquareSum += weff * weff * track.pt() * track.pt();
      }
      // end calculate nch and pt

      // Calculate identified POIs with species corrections. Tracks returning
      // pid == -1 form the unidentified POI and use charged pT/NUE settings.
      if (std::fabs(track.eta()) < trkQualityOpts.cfgRangeEta.value) {
        // ------------------------------
        // Unified PID logic (configurable)
        // ------------------------------
        int pid = getPidConfigurable(track);
        if (pid == -1 && withinPtRef) {
          nUnidentifiedWeighted += weff;
          nUnidentifiedSquare += weff * weff;
          unidentifiedPtSum += weff * track.pt();
          unidentifiedPtSumw2 += weff * weff * track.pt();
          unidentifiedPtSquareSum += weff * weff * track.pt() * track.pt();
        } else if (isWithinPOIPtRange(pid, track.pt())) {
          float weffPid = 1.;
          setParticleNUEWeight(weffPid, track, cent, pid);
          if (pid == MyParticleType::kPion) {
            nPionWeighted += weffPid;
            nPionSquare += weffPid * weffPid;
            pionPtSum += weffPid * track.pt();
            pionPtSumw2 += weffPid * weffPid * track.pt();
            pionPtSquareSum += weffPid * weffPid * track.pt() * track.pt();
          } else if (pid == MyParticleType::kKaon) {
            nKaonWeighted += weffPid;
            nKaonSquare += weffPid * weffPid;
            kaonPtSum += weffPid * track.pt();
            kaonPtSumw2 += weffPid * weffPid * track.pt();
            kaonPtSquareSum += weffPid * weffPid * track.pt() * track.pt();
          } else if (pid == MyParticleType::kProton) {
            nProtonWeighted += weffPid;
            nProtonSquare += weffPid * weffPid;
            protonPtSum += weffPid * track.pt();
            protonPtSumw2 += weffPid * weffPid * track.pt();
            protonPtSquareSum += weffPid * weffPid * track.pt() * track.pt();
          }
        }
      }
      // end calculate POI sums

      if (withinPtRef) {
        nchCorrectedPassed += weff;
        nchPassedSelection += 1;
      }
    }
    // end pt calculation using ITS only

    int totalGlobalTrack = 0;
    // calculate number of pid particle
    int numOfPi = 0;
    int numOfKa = 0;
    int numOfPr = 0;

    /// @note fill v2 flow
    // use global track
    for (const auto& track : tracks) {
      float weff = 1;
      float wacc = 1;

      // do NUE && NUA
      setParticleNUAWeight(wacc, track, vtxz);
      setParticleNUEWeight(weff, track, cent);
      // end do NUE && NUA

      // The reference and POI pT ranges are applied independently below.
      if (!trackSelectedForFlow(track)) {
        continue;
      }
      totalGlobalTrack++;

      if (switchsOpts.cfgDebugMyCode.value && weff == 1.) {
        LOGF(info, "weff for global track is 1, if NUE is open and this appears a lot, check!");
      }

      if (switchsOpts.cfgDebugMyCode.value && wacc == 1.) {
        LOGF(info, "wacc for global track is 1, if NUA is open and this appears alot, check!");
      }

      bool withinPtRefGlobal = isWithinRefPtRange(track.pt());
      if (withinPtRefGlobal) {
        // fill QA hist
        registry.fill(HIST("hPhi"), track.phi());
        registry.fill(HIST("hPhicorr"), track.phi(), wacc);
        registry.fill(HIST("hPhiCorrNUANUE"), track.phi(), wacc * weff);
        registry.fill(HIST("hEta"), track.eta());
        registry.fill(HIST("hPt"), track.pt());
        registry.fill(HIST("hPtCorr"), track.pt(), weff);
        // end fill QA hist

        // pt spectra
        if (switchsOpts.cfgOutPutPtSpectra.value) {
          registry.fill(HIST("ptSpectra/hPtCentData"), track.pt(), cent);
        }

        // Fill the charged reference regions.
        fGFW->Fill(track.eta(), 0, track.phi(), wacc * weff, gfwBitMask::charged::kReference); //(eta, ptbin, phi, wacc*weff, bitmask)
      }

      // ------------------------------
      // Unified PID logic (configurable)
      // ------------------------------
      int pid = getPidConfigurable(track);
      if (particleAbundanceOpts.cfgOutputAbundanceYields.value && withinPtRefGlobal) {
        if (pid == MyParticleType::kPion) {
          registry.fill(HIST("abundance/hPtCentDataPi"), track.pt(), cent);
        } else if (pid == MyParticleType::kKaon) {
          registry.fill(HIST("abundance/hPtCentDataKa"), track.pt(), cent);
        } else if (pid == MyParticleType::kProton) {
          registry.fill(HIST("abundance/hPtCentDataPr"), track.pt(), cent);
        } else {
          registry.fill(HIST("abundance/hPtCentDataUnidentified"), track.pt(), cent);
        }
      }
      if (pid == -1 && withinPtRefGlobal) {
        // Unidentified tracks have no PID-specific correction. Reuse exactly
        // the charged NUA/NUE weights.
        fGFW->Fill(track.eta(), 0, track.phi(), wacc * weff, gfwBitMask::unidentified::kPOI);
        fGFW->Fill(track.eta(), 0, track.phi(), wacc * weff, gfwBitMask::unidentified::kOverlap, wacc * weff);
      } else if (isWithinPOIPtRange(pid, track.pt())) {
        float waccPid = 1;
        float weffPid = 1;
        this->setParticleNUAWeight(waccPid, track, vtxz, pid);
        this->setParticleNUEWeight(weffPid, track, cent, pid);

        // Fill GFW and counters using the same species-specific POI pT range.
        if (pid == MyParticleType::kPion) {
          fGFW->Fill(track.eta(), 0, track.phi(), waccPid * weffPid, gfwBitMask::pion::kPOI);
          // Only tracks in both pT selections contribute to POI-ref overlap.
          if (withinPtRefGlobal) {
            fGFW->Fill(track.eta(), 0, track.phi(), waccPid * weffPid, gfwBitMask::pion::kOverlap, wacc * weff);
          }
          registry.fill(HIST("hPhiPi"), track.phi());
          registry.fill(HIST("hPhicorrPi"), track.phi(), waccPid);
          registry.fill(HIST("hPhiCorrNUANUEPi"), track.phi(), waccPid * weffPid);
          registry.fill(HIST("hPtPi"), track.pt());
          registry.fill(HIST("hPtCorrPi"), track.pt(), weffPid);
          numOfPi++;
        } else if (pid == MyParticleType::kKaon) {
          fGFW->Fill(track.eta(), 0, track.phi(), waccPid * weffPid, gfwBitMask::kaon::kPOI);
          if (withinPtRefGlobal) {
            fGFW->Fill(track.eta(), 0, track.phi(), waccPid * weffPid, gfwBitMask::kaon::kOverlap, wacc * weff);
          }
          registry.fill(HIST("hPhiKa"), track.phi());
          registry.fill(HIST("hPhicorrKa"), track.phi(), waccPid);
          registry.fill(HIST("hPhiCorrNUANUEKa"), track.phi(), waccPid * weffPid);
          registry.fill(HIST("hPtKa"), track.pt());
          registry.fill(HIST("hPtCorrKa"), track.pt(), weffPid);
          numOfKa++;
        } else if (pid == MyParticleType::kProton) {
          fGFW->Fill(track.eta(), 0, track.phi(), waccPid * weffPid, gfwBitMask::proton::kPOI);
          if (withinPtRefGlobal) {
            fGFW->Fill(track.eta(), 0, track.phi(), waccPid * weffPid, gfwBitMask::proton::kOverlap, wacc * weff);
          }
          registry.fill(HIST("hPhiPr"), track.phi());
          registry.fill(HIST("hPhicorrPr"), track.phi(), waccPid);
          registry.fill(HIST("hPhiCorrNUANUEPr"), track.phi(), waccPid * weffPid);
          registry.fill(HIST("hPtPr"), track.pt());
          registry.fill(HIST("hPtCorrPr"), track.pt(), weffPid);
          numOfPr++;
        }
      }
    } // end track loop for v2 calculation

    // sub region, fill graphs after 2 loop on all tracks
    if (particleAbundanceOpts.cfgOutPutAbundanceDis) {
      registry.fill(HIST("abundance/hNumOfPiEventCount"), numOfPi);
      registry.fill(HIST("abundance/hNumOfKaEventCount"), numOfKa);
      registry.fill(HIST("abundance/hNumOfPrEventCount"), numOfPr);
    } // outputabundacedis

    if (switchsOpts.cfgDebugMyCode.value) {
      LOGF(info, Form("global track num %d", totalGlobalTrack));
    }

    registry.fill(HIST("hNchUnCorrectedVSNchCorrected"), nchPassedSelection, nchCorrectedPassed);

    // fill hist using fGFW
    if (nch > 0) {
      fillFC(MyParticleType::kCharged, corrconfigs.at(corrConfigIndex::charged::kC22Gap), cent, rndm, "c22");
      fillFC(MyParticleType::kCharged, corrconfigs.at(corrConfigIndex::charged::kC24), cent, rndm, "c24");
      // Raw gap <4> for SC(ch,ch); reuse c22 for both two-particle factors.
      fillFC(MyParticleType::kCharged, corrconfigs.at(corrConfigIndex::charged::kC24Gap), cent, rndm, "c24Gap");
      fillFC(MyParticleType::kCharged, corrconfigs.at(corrConfigIndex::charged::kC22Full), cent, rndm, "c22Full");
      fillFC(MyParticleType::kCharged, corrconfigs.at(corrConfigIndex::charged::kC32Gap), cent, rndm, "c32");
      fillFC(MyParticleType::kCharged, corrconfigs.at(corrConfigIndex::charged::kC34Gap), cent, rndm, "c34");

      fillFC(MyParticleType::kPion, corrconfigs.at(corrConfigIndex::pion::kC22FullA), cent, rndm, "c22Full");
      fillFC(MyParticleType::kPion, corrconfigs.at(corrConfigIndex::pion::kC22FullB), cent, rndm, "c22Full");
      fillFC(MyParticleType::kKaon, corrconfigs.at(corrConfigIndex::kaon::kC22FullA), cent, rndm, "c22Full");
      fillFC(MyParticleType::kKaon, corrconfigs.at(corrConfigIndex::kaon::kC22FullB), cent, rndm, "c22Full");
      fillFC(MyParticleType::kProton, corrconfigs.at(corrConfigIndex::proton::kC22FullA), cent, rndm, "c22Full");
      fillFC(MyParticleType::kProton, corrconfigs.at(corrConfigIndex::proton::kC22FullB), cent, rndm, "c22Full");
      fillFC(MyParticleType::kUnidentified, corrconfigs.at(corrConfigIndex::unidentified::kC22FullA), cent, rndm, "c22Full");
      fillFC(MyParticleType::kUnidentified, corrconfigs.at(corrConfigIndex::unidentified::kC22FullB), cent, rndm, "c22Full");

      fillFC(MyParticleType::kPion, corrconfigs.at(corrConfigIndex::pion::kC24A), cent, rndm, "c24");
      fillFC(MyParticleType::kPion, corrconfigs.at(corrConfigIndex::pion::kC24B), cent, rndm, "c24");
      fillFC(MyParticleType::kKaon, corrconfigs.at(corrConfigIndex::kaon::kC24A), cent, rndm, "c24");
      fillFC(MyParticleType::kKaon, corrconfigs.at(corrConfigIndex::kaon::kC24B), cent, rndm, "c24");
      fillFC(MyParticleType::kProton, corrconfigs.at(corrConfigIndex::proton::kC24A), cent, rndm, "c24");
      fillFC(MyParticleType::kProton, corrconfigs.at(corrConfigIndex::proton::kC24B), cent, rndm, "c24");
      fillFC(MyParticleType::kUnidentified, corrconfigs.at(corrConfigIndex::unidentified::kC24A), cent, rndm, "c24");
      fillFC(MyParticleType::kUnidentified, corrconfigs.at(corrConfigIndex::unidentified::kC24B), cent, rndm, "c24");

      fillFC(MyParticleType::kPion, corrconfigs.at(corrConfigIndex::pion::kC32GapA), cent, rndm, "c32");
      fillFC(MyParticleType::kPion, corrconfigs.at(corrConfigIndex::pion::kC32GapB), cent, rndm, "c32");
      fillFC(MyParticleType::kKaon, corrconfigs.at(corrConfigIndex::kaon::kC32GapA), cent, rndm, "c32");
      fillFC(MyParticleType::kKaon, corrconfigs.at(corrConfigIndex::kaon::kC32GapB), cent, rndm, "c32");
      fillFC(MyParticleType::kProton, corrconfigs.at(corrConfigIndex::proton::kC32GapA), cent, rndm, "c32");
      fillFC(MyParticleType::kProton, corrconfigs.at(corrConfigIndex::proton::kC32GapB), cent, rndm, "c32");

      fillFC(MyParticleType::kPion, corrconfigs.at(corrConfigIndex::pion::kC34A), cent, rndm, "c34");
      fillFC(MyParticleType::kPion, corrconfigs.at(corrConfigIndex::pion::kC34B), cent, rndm, "c34");
      fillFC(MyParticleType::kKaon, corrconfigs.at(corrConfigIndex::kaon::kC34A), cent, rndm, "c34");
      fillFC(MyParticleType::kKaon, corrconfigs.at(corrConfigIndex::kaon::kC34B), cent, rndm, "c34");
      fillFC(MyParticleType::kProton, corrconfigs.at(corrConfigIndex::proton::kC34A), cent, rndm, "c34");
      fillFC(MyParticleType::kProton, corrconfigs.at(corrConfigIndex::proton::kC34B), cent, rndm, "c34");

      bool filledPi = fillFC(MyParticleType::kPion, corrconfigs.at(corrConfigIndex::pion::kC22Pure), cent, rndm, "c22pure");
      bool filledKa = fillFC(MyParticleType::kKaon, corrconfigs.at(corrConfigIndex::kaon::kC22Pure), cent, rndm, "c22pure");
      bool filledPr = fillFC(MyParticleType::kProton, corrconfigs.at(corrConfigIndex::proton::kC22Pure), cent, rndm, "c22pure");
      bool filledUnidentified = fillFC(MyParticleType::kUnidentified, corrconfigs.at(corrConfigIndex::unidentified::kC22Pure), cent, rndm, "c22pure");

      // SC/NSC components only: the mixed raw <4> uses its GFW tuple weight
      // and the same bootstrap draw as the existing c22pure and charged c22.
      // Existing <2> profiles retain their own valid-event samples/weights.
      fillFC(MyParticleType::kPion, corrconfigs.at(corrConfigIndex::pion::kC24PIDChargedGap), cent, rndm, "c24PIDChargedGap");
      fillFC(MyParticleType::kKaon, corrconfigs.at(corrConfigIndex::kaon::kC24PIDChargedGap), cent, rndm, "c24PIDChargedGap");
      fillFC(MyParticleType::kProton, corrconfigs.at(corrConfigIndex::proton::kC24PIDChargedGap), cent, rndm, "c24PIDChargedGap");
      fillFC(MyParticleType::kUnidentified, corrconfigs.at(corrConfigIndex::unidentified::kC24PIDChargedGap), cent, rndm, "c24PIDChargedGap");

      fillFC(MyParticleType::kPion, corrconfigs.at(corrConfigIndex::pion::kC32Pure), cent, rndm, "c32pure");
      fillFC(MyParticleType::kKaon, corrconfigs.at(corrConfigIndex::kaon::kC32Pure), cent, rndm, "c32pure");
      fillFC(MyParticleType::kProton, corrconfigs.at(corrConfigIndex::proton::kC32Pure), cent, rndm, "c32pure");

      if (filledPi || !switchsOpts.cfgCheck2MethodDiff.value) {
        fillFC(MyParticleType::kPion, corrconfigs.at(corrConfigIndex::pion::kC22GapA), cent, rndm, "c22");
        fillFC(MyParticleType::kPion, corrconfigs.at(corrConfigIndex::pion::kC22GapB), cent, rndm, "c22");
      }
      if (filledKa || !switchsOpts.cfgCheck2MethodDiff.value) {
        fillFC(MyParticleType::kKaon, corrconfigs.at(corrConfigIndex::kaon::kC22GapA), cent, rndm, "c22");
        fillFC(MyParticleType::kKaon, corrconfigs.at(corrConfigIndex::kaon::kC22GapB), cent, rndm, "c22");
      }
      if (filledPr || !switchsOpts.cfgCheck2MethodDiff.value) {
        fillFC(MyParticleType::kProton, corrconfigs.at(corrConfigIndex::proton::kC22GapA), cent, rndm, "c22");
        fillFC(MyParticleType::kProton, corrconfigs.at(corrConfigIndex::proton::kC22GapB), cent, rndm, "c22");
      }
      if (filledUnidentified || !switchsOpts.cfgCheck2MethodDiff.value) {
        fillFC(MyParticleType::kUnidentified, corrconfigs.at(corrConfigIndex::unidentified::kC22GapA), cent, rndm, "c22");
        fillFC(MyParticleType::kUnidentified, corrconfigs.at(corrConfigIndex::unidentified::kC22GapB), cent, rndm, "c22");
      }
      fillFCvnpt(MyParticleType::kCharged, corrconfigs.at(corrConfigIndex::charged::kC22Gap), cent, rndm, nch, nch, "c22TrackWeight");
      fillFCvnpt(MyParticleType::kCharged, corrconfigs.at(corrConfigIndex::charged::kC22Gap), cent, rndm, nch, nch, "c22TrackWeightOne", true);

      fillFCvnpt(MyParticleType::kCharged, corrconfigs.at(corrConfigIndex::charged::kC24), cent, rndm, nch, nch, "c24TrackWeight");
      fillFCvnpt(MyParticleType::kCharged, corrconfigs.at(corrConfigIndex::charged::kC22Full), cent, rndm, nch, nch, "c22FullTrackWeight");
      fillFCvnpt(MyParticleType::kCharged, corrconfigs.at(corrConfigIndex::charged::kC32Gap), cent, rndm, nch, nch, "c32TrackWeight");
      fillFCvnpt(MyParticleType::kCharged, corrconfigs.at(corrConfigIndex::charged::kC34Gap), cent, rndm, nch, nch, "c34TrackWeight");

      fillFCvnpt(MyParticleType::kPion, corrconfigs.at(corrConfigIndex::pion::kC22Pure), cent, rndm, nch, nch, "c22TrackWeight");
      fillFCvnpt(MyParticleType::kKaon, corrconfigs.at(corrConfigIndex::kaon::kC22Pure), cent, rndm, nch, nch, "c22TrackWeight");
      fillFCvnpt(MyParticleType::kProton, corrconfigs.at(corrConfigIndex::proton::kC22Pure), cent, rndm, nch, nch, "c22TrackWeight");
      fillFCvnpt(MyParticleType::kUnidentified, corrconfigs.at(corrConfigIndex::unidentified::kC22Pure), cent, rndm, nch, nch, "c22TrackWeight");

      fillFCvnpt(MyParticleType::kPion, corrconfigs.at(corrConfigIndex::pion::kC32Pure), cent, rndm, nch, nch, "c32TrackWeight");
      fillFCvnpt(MyParticleType::kKaon, corrconfigs.at(corrConfigIndex::kaon::kC32Pure), cent, rndm, nch, nch, "c32TrackWeight");
      fillFCvnpt(MyParticleType::kProton, corrconfigs.at(corrConfigIndex::proton::kC32Pure), cent, rndm, nch, nch, "c32TrackWeight");

      fillFCvnpt(MyParticleType::kCharged, corrconfigs.at(corrConfigIndex::charged::kC22Gap), cent, rndm, ptSum, nch, "covV2Pt");
      fillFCvnpt(MyParticleType::kCharged, corrconfigs.at(corrConfigIndex::charged::kC22Gap), cent, rndm, ptSum, nch, "covV2PtWeightOne", true);

      fillFCvnpt(MyParticleType::kPion, corrconfigs.at(corrConfigIndex::pion::kC22Pure), cent, rndm, ptSum, nch, "covV2Pt");
      fillFCvnpt(MyParticleType::kKaon, corrconfigs.at(corrConfigIndex::kaon::kC22Pure), cent, rndm, ptSum, nch, "covV2Pt");
      fillFCvnpt(MyParticleType::kProton, corrconfigs.at(corrConfigIndex::proton::kC22Pure), cent, rndm, ptSum, nch, "covV2Pt");
      fillFCvnpt(MyParticleType::kUnidentified, corrconfigs.at(corrConfigIndex::unidentified::kC22Pure), cent, rndm, ptSum, nch, "covV2Pt");

      fillFCvnpt(MyParticleType::kCharged, corrconfigs.at(corrConfigIndex::charged::kC32Gap), cent, rndm, ptSum, nch, "covV3Pt");
      fillFCvnpt(MyParticleType::kPion, corrconfigs.at(corrConfigIndex::pion::kC32Pure), cent, rndm, ptSum, nch, "covV3Pt");
      fillFCvnpt(MyParticleType::kKaon, corrconfigs.at(corrConfigIndex::kaon::kC32Pure), cent, rndm, ptSum, nch, "covV3Pt");
      fillFCvnpt(MyParticleType::kProton, corrconfigs.at(corrConfigIndex::proton::kC32Pure), cent, rndm, ptSum, nch, "covV3Pt");

      fillFC4PtC22(cent, ptSum, nch, rndm);
      if (nPionWeighted > 0) {
        fillFC4PtC22(cent, rndm, MyParticleType::kPion, pionPtSum, nPionWeighted);
      }

      if (nKaonWeighted > 0) {
        fillFC4PtC22(cent, rndm, MyParticleType::kKaon, kaonPtSum, nKaonWeighted);
      }

      if (nProtonWeighted > 0) {
        fillFC4PtC22(cent, rndm, MyParticleType::kProton, protonPtSum, nProtonWeighted);
      }
      if (nUnidentifiedWeighted > 0) {
        fillFC4PtC22(cent, rndm, MyParticleType::kUnidentified, unidentifiedPtSum, nUnidentifiedWeighted);
      }

      if (switchsOpts.cfgOutPutPtSpectra.value) {
        // charged calculation
        double NpairCharged = fGFW->Calculate(corrconfigs.at(corrConfigIndex::charged::kC22Gap), 0, true).real();
        double chargedC22 = NpairCharged > 0 ? fGFW->Calculate(corrconfigs.at(corrConfigIndex::charged::kC22Gap), 0, false).real() / NpairCharged : 0;
        // end charged calculation

        // pi
        double pidChargedC22Pi = getPidC22InOneEvent(corrconfigs.at(corrConfigIndex::pion::kC22GapA), corrconfigs.at(corrConfigIndex::pion::kC22GapB));
        if (pidChargedC22Pi > 0 && chargedC22 > 0) {
          registry.fill(HIST("c22PrimeVsc22/Pi"), pidChargedC22Pi, chargedC22);
        }
        // end pi

        // Ka
        double pidKaonC22 = getPidC22InOneEvent(corrconfigs.at(corrConfigIndex::kaon::kC22GapA), corrconfigs.at(corrConfigIndex::kaon::kC22GapB));
        if (pidKaonC22 > 0 && chargedC22 > 0) {
          registry.fill(HIST("c22PrimeVsc22/Ka"), pidKaonC22, chargedC22);
        }
        // end Ka

        // Pr
        double pidProtonC22 = getPidC22InOneEvent(corrconfigs.at(corrConfigIndex::proton::kC22GapA), corrconfigs.at(corrConfigIndex::proton::kC22GapB));
        if (pidProtonC22 > 0 && chargedC22 > 0) {
          registry.fill(HIST("c22PrimeVsc22/Pr"), pidProtonC22, chargedC22);
        }
        // end Pr
      }

      if (switchsOpts.cfgOutputCharged.value) {
        fFCCh->FillProfile("hMeanPt", cent, (ptSum / nch), nch, rndm);
        fFCCh->FillProfile("hMeanPtWeightOne", cent, (ptSum / nch), 1., rndm);
      }

      if (nPionWeighted > 0 && switchsOpts.cfgOutputPion.value) {
        fFCPi->FillProfile("hMeanPt", cent, (pionPtSum / nPionWeighted), nPionWeighted, rndm);
        fFCPi->FillProfile("hMeanPtWeightOne", cent, (pionPtSum / nPionWeighted), 1., rndm);
        fillFCvnpt(MyParticleType::kPion, corrconfigs.at(corrConfigIndex::pion::kC22Pure), cent, rndm, pionPtSum, nPionWeighted, "covV2PtPID");
        fillFCvnpt(MyParticleType::kPion, corrconfigs.at(corrConfigIndex::pion::kC22Pure), cent, rndm, nPionWeighted, nPionWeighted, "c22TrackWeightPID");
      }

      if (nKaonWeighted > 0 && switchsOpts.cfgOutputKaon.value) {
        fFCKa->FillProfile("hMeanPt", cent, (kaonPtSum / nKaonWeighted), nKaonWeighted, rndm);
        fFCKa->FillProfile("hMeanPtWeightOne", cent, (kaonPtSum / nKaonWeighted), 1., rndm);
        fillFCvnpt(MyParticleType::kKaon, corrconfigs.at(corrConfigIndex::kaon::kC22Pure), cent, rndm, kaonPtSum, nKaonWeighted, "covV2PtPID");
        fillFCvnpt(MyParticleType::kKaon, corrconfigs.at(corrConfigIndex::kaon::kC22Pure), cent, rndm, nKaonWeighted, nKaonWeighted, "c22TrackWeightPID");
      }

      if (nProtonWeighted > 0 && switchsOpts.cfgOutputProton.value) {
        fFCPr->FillProfile("hMeanPt", cent, (protonPtSum / nProtonWeighted), nProtonWeighted, rndm);
        fFCPr->FillProfile("hMeanPtWeightOne", cent, (protonPtSum / nProtonWeighted), 1., rndm);
        fillFCvnpt(MyParticleType::kProton, corrconfigs.at(corrConfigIndex::proton::kC22Pure), cent, rndm, protonPtSum, nProtonWeighted, "covV2PtPID");
        fillFCvnpt(MyParticleType::kProton, corrconfigs.at(corrConfigIndex::proton::kC22Pure), cent, rndm, nProtonWeighted, nProtonWeighted, "c22TrackWeightPID");
      }

      if (nUnidentifiedWeighted > 0 && switchsOpts.cfgOutputUnidentified.value) {
        fFCUnidentified->FillProfile("hMeanPt", cent, unidentifiedPtSum / nUnidentifiedWeighted, nUnidentifiedWeighted, rndm);
        fFCUnidentified->FillProfile("hMeanPtWeightOne", cent, unidentifiedPtSum / nUnidentifiedWeighted, 1., rndm);
        fillFCvnpt(MyParticleType::kUnidentified, corrconfigs.at(corrConfigIndex::unidentified::kC22Pure), cent, rndm, unidentifiedPtSum, nUnidentifiedWeighted, "covV2PtPID");
        fillFCvnpt(MyParticleType::kUnidentified, corrconfigs.at(corrConfigIndex::unidentified::kC22Pure), cent, rndm, nUnidentifiedWeighted, nUnidentifiedWeighted, "c22TrackWeightPID");
      }

      double nchDiff = nch * nch - nchSquare;
      if (nchDiff > minVal4Float && switchsOpts.cfgOutputCharged.value) {
        fFCCh->FillProfile("ptSquareAve", cent,
                           (ptSum * ptSum - ptSquareSum) / nchDiff,
                           nchDiff, rndm);

        fFCCh->FillProfile("ptSquareAveWeightOne", cent,
                           (ptSum * ptSum - ptSquareSum) / nchDiff,
                           1., rndm);

        fFCCh->FillProfile("ptAve", cent,
                           (nch * ptSum - ptSumw2) / nchDiff,
                           nchDiff, rndm);

        fFCCh->FillProfile("ptAveWeightOne", cent,
                           (nch * ptSum - ptSumw2) / nchDiff,
                           1., rndm);
      }

      double pionDiff = nPionWeighted * nPionWeighted - nPionSquare;
      if (pionDiff > minVal4Float && switchsOpts.cfgOutputPion.value) {
        fFCPi->FillProfile("ptSquareAve", cent,
                           (pionPtSum * pionPtSum - pionPtSquareSum) / pionDiff,
                           pionDiff, rndm);

        fFCPi->FillProfile("ptSquareAveWeightOne", cent,
                           (pionPtSum * pionPtSum - pionPtSquareSum) / pionDiff,
                           1., rndm);

        fFCPi->FillProfile("ptAve", cent,
                           (nPionWeighted * pionPtSum - pionPtSumw2) / pionDiff,
                           pionDiff, rndm);

        fFCPi->FillProfile("ptAveWeightOne", cent,
                           (nPionWeighted * pionPtSum - pionPtSumw2) / pionDiff,
                           1., rndm);
      }

      double kaonDiff = nKaonWeighted * nKaonWeighted - nKaonSquare;
      if (kaonDiff > minVal4Float && switchsOpts.cfgOutputKaon.value) {
        fFCKa->FillProfile("ptSquareAve", cent,
                           (kaonPtSum * kaonPtSum - kaonPtSquareSum) / kaonDiff,
                           kaonDiff, rndm);

        fFCKa->FillProfile("ptSquareAveWeightOne", cent,
                           (kaonPtSum * kaonPtSum - kaonPtSquareSum) / kaonDiff,
                           1., rndm);

        fFCKa->FillProfile("ptAve", cent,
                           (nKaonWeighted * kaonPtSum - kaonPtSumw2) / kaonDiff,
                           kaonDiff, rndm);

        fFCKa->FillProfile("ptAveWeightOne", cent,
                           (nKaonWeighted * kaonPtSum - kaonPtSumw2) / kaonDiff,
                           1., rndm);
      }

      double protonDiff = nProtonWeighted * nProtonWeighted - nProtonSquare;
      if (protonDiff > minVal4Float && switchsOpts.cfgOutputProton.value) {
        fFCPr->FillProfile("ptSquareAve", cent,
                           (protonPtSum * protonPtSum - protonPtSquareSum) / protonDiff,
                           protonDiff, rndm);

        fFCPr->FillProfile("ptSquareAveWeightOne", cent,
                           (protonPtSum * protonPtSum - protonPtSquareSum) / protonDiff,
                           1., rndm);

        fFCPr->FillProfile("ptAve", cent,
                           (nProtonWeighted * protonPtSum - protonPtSumw2) / protonDiff,
                           protonDiff, rndm);

        fFCPr->FillProfile("ptAveWeightOne", cent,
                           (nProtonWeighted * protonPtSum - protonPtSumw2) / protonDiff,
                           1., rndm);
      }

      double unidentifiedDiff = nUnidentifiedWeighted * nUnidentifiedWeighted - nUnidentifiedSquare;
      if (unidentifiedDiff > minVal4Float && switchsOpts.cfgOutputUnidentified.value) {
        fFCUnidentified->FillProfile("ptSquareAve", cent,
                                     (unidentifiedPtSum * unidentifiedPtSum - unidentifiedPtSquareSum) / unidentifiedDiff,
                                     unidentifiedDiff, rndm);
        fFCUnidentified->FillProfile("ptSquareAveWeightOne", cent,
                                     (unidentifiedPtSum * unidentifiedPtSum - unidentifiedPtSquareSum) / unidentifiedDiff,
                                     1., rndm);
        fFCUnidentified->FillProfile("ptAve", cent,
                                     (nUnidentifiedWeighted * unidentifiedPtSum - unidentifiedPtSumw2) / unidentifiedDiff,
                                     unidentifiedDiff, rndm);
        fFCUnidentified->FillProfile("ptAveWeightOne", cent,
                                     (nUnidentifiedWeighted * unidentifiedPtSum - unidentifiedPtSumw2) / unidentifiedDiff,
                                     1., rndm);
      }

    } // end fill hist using fillProfile
  }
  PROCESS_SWITCH(PidFlowPtCorr, processData, "", true);

  void initTimeAssociationComparisonQA()
  {
    constexpr const char* selection = "Common selected events, final charged selection, unit weights";
    registry.add("associationComparison/hEventAvailability", "Common events including empty track sets;Centrality (%);Track availability", {HistType::kTH2D, {axisMultiplicity, {4, -0.5, 3.5}}});
    registry.add("associationComparison/hNchOriginalVsTime", "Final charged multiplicity;N_{ch} original;N_{ch} time;Event count", {HistType::kTH2D, {cfgaxisAssociationNch, cfgaxisAssociationNch}});
    registry.add("associationComparison/hLostVsGained", "Track membership changes;Lost tracks;Gained tracks;Event count", {HistType::kTH2D, {cfgaxisAssociationChangedNch, cfgaxisAssociationChangedNch}});
    registry.add("associationComparison/hDeltaNchVsCent", "Multiplicity change;Centrality (%);N_{ch}^{time}-N_{ch}^{original};Event count", {HistType::kTH2D, {axisMultiplicity, cfgaxisAssociationDeltaNch}});
    registry.add("associationComparison/hDeltaNchVsOccupancy", "Multiplicity change;Track occupancy;N_{ch}^{time}-N_{ch}^{original};Event count", {HistType::kTH2D, {cfgaxisAssociationOccupancy, cfgaxisAssociationDeltaNch}});
    registry.add("associationComparison/hSetChangeVsCent", "Nonempty union only;Centrality (%);(lost+gained)/(common+lost+gained);Event count", {HistType::kTH2D, {axisMultiplicity, {101, -0.005, 1.005}}});
    registry.add("associationComparison/hSetChangeVsOccupancy", "Nonempty union only;Track occupancy;(lost+gained)/(common+lost+gained);Event count", {HistType::kTH2D, {cfgaxisAssociationOccupancy, {101, -0.005, 1.005}}});
    registry.add("associationComparison/pSetChangeVsCent", "Event mean, nonempty union only;Centrality (%);Mean set-change fraction", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("associationComparison/pSetChangeVsOccupancy", "Event mean, nonempty union only;Track occupancy;Mean set-change fraction", {HistType::kTProfile, {cfgaxisAssociationOccupancy}});
    registry.add("associationComparison/pLostFractionVsCent", "Event mean, original Nch > 0;Centrality (%);Mean lost/original fraction", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("associationComparison/pLostFractionVsOccupancy", "Event mean, original Nch > 0;Track occupancy;Mean lost/original fraction", {HistType::kTProfile, {cfgaxisAssociationOccupancy}});
    registry.add("associationComparison/pGainedFractionVsCent", "Event mean, time Nch > 0;Centrality (%);Mean gained/time fraction", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("associationComparison/pGainedFractionVsOccupancy", "Event mean, time Nch > 0;Track occupancy;Mean gained/time fraction", {HistType::kTProfile, {cfgaxisAssociationOccupancy}});
    // Category projections provide original/time denominators for lost/gained
    // spectra ratios. Entries count track-target pairs, not distinct sources.
    registry.add("associationComparison/hPtByCategory", "Final tracks per target collision;p_{T} (GeV/c);Track category;Track-target count", {HistType::kTH2D, {cfgaxisPt, {5, -0.5, 4.5}}});
    registry.add("associationComparison/hEtaByCategory", "Final tracks per target collision;#eta;Track category;Track-target count", {HistType::kTH2D, {cfgaxisEta, {5, -0.5, 4.5}}});
    registry.add("associationComparison/hPhiByCategory", "Final tracks per target collision;#varphi;Track category;Track-target count", {HistType::kTH2D, {cfgaxisPhi, {5, -0.5, 4.5}}});
    registry.add("associationComparison/hFinalCollisionsOriginal", "All source tracks once per dataframe;N selected collisions (original);Source track count", {HistType::kTH1D, {{2, -0.5, 1.5}}});
    registry.add("associationComparison/hFinalCollisionsTime", "All source tracks once per dataframe;N selected collisions (time);Source track count", {HistType::kTH1D, {{101, -0.5, 100.5}}});
    registry.add("associationComparison/hFinalCollisionsOriginalVsTime", "All source tracks once per dataframe;N selected collisions (original);N selected collisions (time);Source track count", {HistType::kTH2D, {{2, -0.5, 1.5}, {101, -0.5, 100.5}}});
    const auto setLabels = [](TAxis* axis, auto const& labels) {
      for (std::size_t index = 0; index < labels.size(); ++index) {
        axis->SetBinLabel(index + 1, labels[index]);
      }
    };
    setLabels(registry.get<TH2>(HIST("associationComparison/hEventAvailability"))->GetYaxis(),
              std::array{"both empty", "original only", "time only", "both nonempty"});
    const std::array categories{"original", "time", "common", "lost", "gained"};
    setLabels(registry.get<TH2>(HIST("associationComparison/hPtByCategory"))->GetYaxis(), categories);
    setLabels(registry.get<TH2>(HIST("associationComparison/hEtaByCategory"))->GetYaxis(), categories);
    setLabels(registry.get<TH2>(HIST("associationComparison/hPhiByCategory"))->GetYaxis(), categories);
    LOG(info) << "Association comparison: " << selection << "; original=collisionId, time=PidFlowTimeTracks; unique=" << cfgRequireUniqueTimeAssociation.value;
  }

  // Evaluate both associations in the same event loop without PID or correction
  // weights. Keep the full source table to include zero-use tracks in usage QA.
  void processTimeAssociationComparison(AodCollisions const& collisions, aod::BCsWithTimestamps const&,
                                        PidFlowChargedTimeProducer::SourceTracks const& sources, aod::PidFlowTimeTracks const& associatedTracks)
  {
    if (!cfgTimeAssociationQA) {
      return;
    }
    std::vector<unsigned int> originalUses(sources.size(), 0), timeUses(sources.size(), 0);
    const auto sourceIndex = [&](int64_t trackId) {
      const int64_t index = trackId - sources.offset();
      if (index < 0 || index >= sources.size()) {
        LOGF(fatal, "Association comparison track index is outside the full source track table");
      }
      return index;
    };
    const auto fillKinematics = [&](auto const& track, int category) {
      registry.fill(HIST("associationComparison/hPtByCategory"), track.pt(), category);
      registry.fill(HIST("associationComparison/hEtaByCategory"), track.eta(), category);
      registry.fill(HIST("associationComparison/hPhiByCategory"), track.phi(), category);
    };
    for (auto const& collision : collisions) {
      // No track-presence requirement: both methods share this event sample.
      if (!chargedTimeEventSelected(collision, funcProcessTimeAssociationComparison)) {
        continue;
      }
      auto originalTracks = sources.sliceBy(perCollision, collision.globalIndex());
      auto timeTracks = associatedTracks.sliceBy(timeTracksPerCollision, collision.globalIndex());
      std::unordered_set<int64_t> originalIds, timeIds;
      for (auto const& track : originalTracks) {
        // Reproduce processData's framework filter and final charged cuts,
        // without requiring or evaluating any PID tables.
        if (!(track.isGlobalTrack() || track.isGlobalTrackSDD() != 0) || !(std::fabs(track.eta()) < trkQualityOpts.cfgCutEta.value) ||
            !trackSelectedForFlow(track) || !isWithinRefPtRange(track.pt())) {
          continue;
        }
        if (originalIds.insert(track.globalIndex()).second) {
          ++originalUses[sourceIndex(track.globalIndex())];
          fillKinematics(track, 0);
        }
      }
      unsigned int common = 0;
      for (auto const& track : timeTracks) {
        if (!trackSelectedForFlow(track) || !isWithinRefPtRange(track.pt()) || !timeIds.insert(track.trackId()).second) {
          continue;
        }
        ++timeUses[sourceIndex(track.trackId())];
        fillKinematics(track, 1);
        if (originalIds.contains(track.trackId())) {
          ++common;
          fillKinematics(track, 2);
        } else {
          fillKinematics(track, 4);
        }
      }
      for (auto const& track : originalTracks) {
        if (originalIds.contains(track.globalIndex()) && !timeIds.contains(track.globalIndex())) {
          fillKinematics(track, 3);
        }
      }
      const auto nOriginal = originalIds.size();
      const auto nTime = timeIds.size();
      const auto lost = nOriginal - common;
      const auto gained = nTime - common;
      const auto unionSize = common + lost + gained;
      const double cent = getCentrality(collision);
      const auto occupancy = collision.trackOccupancyInTimeRange();
      const double deltaNch = static_cast<double>(nTime) - static_cast<double>(nOriginal);
      const int availability = (nOriginal > 0 ? 1 : 0) + (nTime > 0 ? 2 : 0);
      registry.fill(HIST("associationComparison/hEventAvailability"), cent, availability);
      registry.fill(HIST("associationComparison/hNchOriginalVsTime"), nOriginal, nTime);
      registry.fill(HIST("associationComparison/hLostVsGained"), lost, gained);
      registry.fill(HIST("associationComparison/hDeltaNchVsCent"), cent, deltaNch);
      registry.fill(HIST("associationComparison/hDeltaNchVsOccupancy"), occupancy, deltaNch);
      if (unionSize > 0) {
        const double changeFraction = static_cast<double>(lost + gained) / unionSize;
        registry.fill(HIST("associationComparison/hSetChangeVsCent"), cent, changeFraction);
        registry.fill(HIST("associationComparison/hSetChangeVsOccupancy"), occupancy, changeFraction);
        registry.fill(HIST("associationComparison/pSetChangeVsCent"), cent, changeFraction);
        registry.fill(HIST("associationComparison/pSetChangeVsOccupancy"), occupancy, changeFraction);
      }
      if (nOriginal > 0) {
        const double lostFraction = static_cast<double>(lost) / nOriginal;
        registry.fill(HIST("associationComparison/pLostFractionVsCent"), cent, lostFraction);
        registry.fill(HIST("associationComparison/pLostFractionVsOccupancy"), occupancy, lostFraction);
      }
      if (nTime > 0) {
        const double gainedFraction = static_cast<double>(gained) / nTime;
        registry.fill(HIST("associationComparison/pGainedFractionVsCent"), cent, gainedFraction);
        registry.fill(HIST("associationComparison/pGainedFractionVsOccupancy"), occupancy, gainedFraction);
      }
    }
    for (std::size_t index = 0; index < originalUses.size(); ++index) {
      registry.fill(HIST("associationComparison/hFinalCollisionsOriginal"), originalUses[index]);
      registry.fill(HIST("associationComparison/hFinalCollisionsTime"), timeUses[index]);
      registry.fill(HIST("associationComparison/hFinalCollisionsOriginalVsTime"), originalUses[index], timeUses[index]);
    }
  }
  PROCESS_SWITCH(PidFlowPtCorr, processTimeAssociationComparison, "Compare original collisionId and time-associated final charged tracks; no PID/NUA/NUE", false);

  // Independent charged-only counterpart of processData. All moment/profile
  // definitions below retain the original charged weights and bootstrap draw.
  void processDataTimeAssociation(AodCollisions::iterator const& collision, aod::BCsWithTimestamps const&, aod::PidFlowTimeTracks const& tracks)
  {
    if (!chargedTimeEventSelected(collision, funcProcessDataTimeAssociation)) {
      return;
    }
    if (tracks.size() == 0) {
      return;
    }
    fillEventCountHelper(funcProcessDataTimeAssociation, 14.5);
    const double cent = getCentrality(collision);
    const float vtxz = collision.posZ();
    const float rndm = fRndm->Rndm();
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    loadChargedTimeCorrections(bc.timestamp());
    fGFW->Clear();
    registry.fill(HIST("hVtxZ"), vtxz);
    registry.fill(HIST("hMult"), tracks.size());
    registry.fill(HIST("hMultTPC"), collision.multTPC());
    registry.fill(HIST("hCent"), cent);
    if (switchsOpts.cfgOutPutPtSpectra) {
      registry.fill(HIST("ptSpectra/hCentEventCountData"), cent);
    }

    double ptSum = 0., ptSumw2 = 0., nch = 0., nchSquare = 0., ptSquareSum = 0.;
    double nchPassedSelection = 0., nchCorrectedPassed = 0.;
    for (auto const& track : tracks) {
      if (!trackSelectedForFlow(track) || !isWithinRefPtRange(track.pt())) {
        continue;
      }
      float weff = 1.f;
      setParticleNUEWeight(weff, track, cent);
      nchPassedSelection += 1.;
      nchCorrectedPassed += weff;
      if (std::fabs(track.eta()) < trkQualityOpts.cfgRangeEta.value) {
        nch += weff;
        nchSquare += weff * weff;
        ptSum += weff * track.pt();
        ptSumw2 += weff * weff * track.pt();
        ptSquareSum += weff * weff * track.pt() * track.pt();
      }
    }
    for (auto const& track : tracks) {
      if (!trackSelectedForFlow(track) || !isWithinRefPtRange(track.pt())) {
        continue;
      }
      float weff = 1.f, wacc = 1.f;
      setParticleNUEWeight(weff, track, cent);
      setParticleNUAWeight(wacc, track, vtxz);
      registry.fill(HIST("hPhi"), track.phi());
      registry.fill(HIST("hPhicorr"), track.phi(), wacc);
      registry.fill(HIST("hPhiCorrNUANUE"), track.phi(), wacc * weff);
      registry.fill(HIST("hEta"), track.eta());
      registry.fill(HIST("hPt"), track.pt());
      registry.fill(HIST("hPtCorr"), track.pt(), weff);
      if (switchsOpts.cfgOutPutPtSpectra) {
        registry.fill(HIST("ptSpectra/hPtCentData"), track.pt(), cent);
      }
      fGFW->Fill(track.eta(), 0, track.phi(), wacc * weff, gfwBitMask::charged::kReference);
    }
    registry.fill(HIST("hNchUnCorrectedVSNchCorrected"), nchPassedSelection, nchCorrectedPassed);
    if (nch <= 0.) {
      return;
    }
    fillFC(MyParticleType::kCharged, corrconfigs.at(corrConfigIndex::charged::kC22Gap), cent, rndm, "c22");
    fillFC(MyParticleType::kCharged, corrconfigs.at(corrConfigIndex::charged::kC24), cent, rndm, "c24");
    fillFC(MyParticleType::kCharged, corrconfigs.at(corrConfigIndex::charged::kC24Gap), cent, rndm, "c24Gap");
    fillFC(MyParticleType::kCharged, corrconfigs.at(corrConfigIndex::charged::kC22Full), cent, rndm, "c22Full");
    fillFC(MyParticleType::kCharged, corrconfigs.at(corrConfigIndex::charged::kC32Gap), cent, rndm, "c32");
    fillFC(MyParticleType::kCharged, corrconfigs.at(corrConfigIndex::charged::kC34Gap), cent, rndm, "c34");
    fillFCvnpt(MyParticleType::kCharged, corrconfigs.at(corrConfigIndex::charged::kC22Gap), cent, rndm, nch, nch, "c22TrackWeight");
    fillFCvnpt(MyParticleType::kCharged, corrconfigs.at(corrConfigIndex::charged::kC22Gap), cent, rndm, nch, nch, "c22TrackWeightOne", true);
    fillFCvnpt(MyParticleType::kCharged, corrconfigs.at(corrConfigIndex::charged::kC24), cent, rndm, nch, nch, "c24TrackWeight");
    fillFCvnpt(MyParticleType::kCharged, corrconfigs.at(corrConfigIndex::charged::kC22Full), cent, rndm, nch, nch, "c22FullTrackWeight");
    fillFCvnpt(MyParticleType::kCharged, corrconfigs.at(corrConfigIndex::charged::kC32Gap), cent, rndm, nch, nch, "c32TrackWeight");
    fillFCvnpt(MyParticleType::kCharged, corrconfigs.at(corrConfigIndex::charged::kC34Gap), cent, rndm, nch, nch, "c34TrackWeight");
    fillFCvnpt(MyParticleType::kCharged, corrconfigs.at(corrConfigIndex::charged::kC22Gap), cent, rndm, ptSum, nch, "covV2Pt");
    fillFCvnpt(MyParticleType::kCharged, corrconfigs.at(corrConfigIndex::charged::kC22Gap), cent, rndm, ptSum, nch, "covV2PtWeightOne", true);
    fillFCvnpt(MyParticleType::kCharged, corrconfigs.at(corrConfigIndex::charged::kC32Gap), cent, rndm, ptSum, nch, "covV3Pt");
    fillFC4PtC22(cent, ptSum, nch, rndm);
    fFCCh->FillProfile("hMeanPt", cent, ptSum / nch, nch, rndm);
    fFCCh->FillProfile("hMeanPtWeightOne", cent, ptSum / nch, 1., rndm);
    const double nchDiff = nch * nch - nchSquare;
    if (nchDiff > minVal4Float) {
      fFCCh->FillProfile("ptSquareAve", cent, (ptSum * ptSum - ptSquareSum) / nchDiff, nchDiff, rndm);
      fFCCh->FillProfile("ptSquareAveWeightOne", cent, (ptSum * ptSum - ptSquareSum) / nchDiff, 1., rndm);
      fFCCh->FillProfile("ptAve", cent, (nch * ptSum - ptSumw2) / nchDiff, nchDiff, rndm);
      fFCCh->FillProfile("ptAveWeightOne", cent, (nch * ptSum - ptSumw2) / nchDiff, 1., rndm);
    }
    if (switchsOpts.cfgAddC22DeltaPtHistograms) {
      registry.fill(HIST("c22DeltaPt/hMeanPtCharged"), cent, rndm * cfgFlowNbootstrap, ptSum / nch, nch);
      fillProfilePOIvnpt(corrconfigs.at(corrConfigIndex::charged::kC22Gap), HIST("c22DeltaPt/c22dmeanptCharged"), cent, ptSum, nch, rndm);
    }
  }
  PROCESS_SWITCH(PidFlowPtCorr, processDataTimeAssociation, "Charged flow and mean-pT using time-based TrackAssoc; no PID", false);

  void fillCorrectionGraphTimeAssociation(AodCollisions::iterator const& collision, aod::BCsWithTimestamps const&, aod::PidFlowTimeTracks const& tracks)
  {
    if (!chargedTimeEventSelected(collision, funcFillCorrectionGraphTimeAssociation)) {
      return;
    }
    if (tracks.size() == 0) {
      return;
    }
    fillEventCountHelper(funcFillCorrectionGraphTimeAssociation, 14.5);
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    auto run = std::find(runNumbers.begin(), runNumbers.end(), bc.runNumber());
    if (run == runNumbers.end()) {
      return;
    }
    const int runIndex = std::distance(runNumbers.begin(), run);
    const double cent = getCentrality(collision);
    if (switchsOpts.cfgUseNUAWithPt) {
      loadChargedTimeCorrections(bc.timestamp(), false);
    }
    for (auto const& track : tracks) {
      if (!trackSelectedForFlow(track) || !isWithinRefPtRange(track.pt())) {
        continue;
      }
      float weff = 1.f;
      if (switchsOpts.cfgUseNUAWithPt) {
        setParticleNUEWeight(weff, track, cent);
      }
      registry.fill(HIST("correction/hRunNumberPhiEtaVertex"), runIndex, track.phi(), track.eta(), collision.posZ(), track.pt(), weff);
    }
  }
  PROCESS_SWITCH(PidFlowPtCorr, fillCorrectionGraphTimeAssociation, "Charged NUA source using time-based TrackAssoc; no PID", false);

  // Fill reco numerator and generated denominator together, so they use the
  // very same accepted target event. A different true collision is contamination,
  // not successful reconstruction of a particle in this denominator.
  void processMCTimeAssociation(ChargedTimeCollisionsMC const& collisions, aod::BCsWithTimestamps const&, aod::PidFlowTimeTracks const& associatedTracks,
                                aod::McCollisions const& mcCollisions, aod::McParticles const& particles, aod::McTrackLabels const& labels)
  {
    std::vector<unsigned int> recoPerMC(mcCollisions.size(), 0);
    for (auto const& collision : collisions) {
      if (collision.has_mcCollision()) {
        ++recoPerMC[collision.mcCollisionId() - mcCollisions.offset()];
      }
    }
    for (auto const& collision : collisions) {
      if (cfgTimeAssociationQA) {
        registry.fill(HIST("timeAssociation/mcEventStatus"), 0.5);
      }
      if (!collision.has_mcCollision()) {
        continue;
      }
      if (cfgTimeAssociationQA) {
        registry.fill(HIST("timeAssociation/mcEventStatus"), 1.5);
      }
      const int64_t mcIndex = collision.mcCollisionId() - mcCollisions.offset();
      auto mcCollision = mcCollisions.rawIteratorAt(mcIndex);
      if (recoPerMC[mcIndex] != 1 || std::fabs(mcCollision.posZ()) >= cfgCutVertex.value) {
        continue;
      }
      if (cfgTimeAssociationQA) {
        registry.fill(HIST("timeAssociation/mcEventStatus"), 2.5);
      }
      if (!chargedTimeEventSelected(collision, funcProcessMCTimeAssociation)) {
        continue;
      }
      auto tracks = associatedTracks.sliceBy(timeTracksPerCollision, collision.globalIndex());
      if (tracks.size() == 0) {
        continue;
      }
      fillEventCountHelper(funcProcessMCTimeAssociation, 14.5);
      if (cfgTimeAssociationQA) {
        registry.fill(HIST("timeAssociation/mcEventStatus"), 3.5);
      }
      const double cent = getCentrality(collision);
      if (switchsOpts.cfgOutPutPtSpectra) {
        registry.fill(HIST("ptSpectra/hCentEventCountMcRec"), cent);
        registry.fill(HIST("ptSpectra/hCentEventCountMcGen"), cent);
      }
      for (auto const& particle : particles.sliceBy(particlesPerMCCollision, collision.mcCollisionId())) {
        if (particle.isPhysicalPrimary() && particleSelected(particle)) {
          registry.fill(HIST("correction/hPtCentMcGen"), particle.pt(), cent);
        }
      }
      for (auto const& track : tracks) {
        if (!trackSelectedForFlow(track) || !isWithinRefPtRange(track.pt())) {
          continue;
        }
        if (cfgTimeAssociationQA) {
          registry.fill(HIST("timeAssociation/hPtCentMcSelected"), track.pt(), cent);
        }
        const int64_t labelIndex = track.trackId() - labels.offset();
        if (labelIndex < 0 || labelIndex >= labels.size()) {
          LOGF(fatal, "Time-associated source track index is outside McTrackLabels");
        }
        auto source = labels.rawIteratorAt(labelIndex);
        if (!source.has_mcParticle()) {
          if (cfgTimeAssociationQA) {
            registry.fill(HIST("timeAssociation/mcTrackStatus"), 0.5);
          }
          continue;
        }
        auto particle = source.mcParticle_as<aod::McParticles>();
        if (particle.mcCollisionId() != collision.mcCollisionId()) {
          if (cfgTimeAssociationQA) {
            registry.fill(HIST("timeAssociation/mcTrackStatus"), 1.5);
            registry.fill(HIST("timeAssociation/hPtCentMcOtherCollision"), track.pt(), cent);
          }
          continue;
        }
        if (!particle.isPhysicalPrimary() || !particleSelected(particle)) {
          if (cfgTimeAssociationQA) {
            registry.fill(HIST("timeAssociation/mcTrackStatus"), 2.5);
          }
          continue;
        }
        if (cfgTimeAssociationQA) {
          registry.fill(HIST("timeAssociation/mcTrackStatus"), 3.5);
        }
        registry.fill(HIST("correction/hPtCentMcRec"), track.pt(), cent);
      }
    }
  }
  PROCESS_SWITCH(PidFlowPtCorr, processMCTimeAssociation, "Charged NUE reco/gen with common time-association event sample; no detector PID", false);

  /**
   * @brief Fill charged, identified-particle, and unidentified-particle c22 versus event mean pT.
   * @note Each profile uses the mean pT of the same particle sample on its pT axis:
   *       charged, pion, kaon, proton, or PID-failed/ambiguous unidentified tracks.
   *       The corresponding FlowContainer hMeanPt
   *       observable provides the ensemble mean pT needed to convert that axis to
   *       delta[pT]/<[pT]> offline.
   */
  void processDataC22DeltaPt(AodCollisions::iterator const& collision, aod::BCsWithTimestamps const&, AodTracks const& tracks)
  {
    if (!switchsOpts.cfgAddC22DeltaPtHistograms.value) {
      return;
    }
    const double rndm = fRndm->Rndm();

    /// @note event selection
    registry.fill(HIST("hEventCount/processDataC22DeltaPt"), 0.5);
    if (tracks.size() < 1) {
      return;
    }

    fGFW->Clear();
    const auto cent = getCentrality(collision);
    if (!collision.sel8()) {
      return;
    }
    registry.fill(HIST("hEventCount/processDataC22DeltaPt"), 1.5);

    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    const double interactionRate = getInteractionRate(bc.timestamp(), bc.runNumber());
    if (!eventSelected(collision, cent, interactionRate, MyFunctionName::funcProcessDataC22DeltaPt)) {
      return;
    }
    // event selection

    loadCorrections(bc.timestamp());

    const float vtxz = collision.posZ();

    // Event-wise mean-pT numerators and NUE-corrected particle counts. Charged
    // particles use the reference pT range; each PID species uses its own POI pT range.
    double ptSum = 0.;
    double nch = 0.;
    double pionPtSum = 0.;
    double kaonPtSum = 0.;
    double protonPtSum = 0.;
    double unidentifiedPtSum = 0.;
    double nPionWeighted = 0.;
    double nKaonWeighted = 0.;
    double nProtonWeighted = 0.;
    double nUnidentifiedWeighted = 0.;

    for (const auto& track : tracks) {
      // The mean-pT samples use the same track-quality and eta selections. The pT
      // selection and NUE correction are then applied independently for charged and PID.
      if (!trackSelectedForFlow(track) || std::fabs(track.eta()) >= trkQualityOpts.cfgRangeEta.value) {
        continue;
      }

      float weff = 1.f;
      setParticleNUEWeight(weff, track, cent);
      const bool withinPtRef = isWithinRefPtRange(track.pt());
      if (withinPtRef) {
        nch += weff;
        ptSum += weff * track.pt();
      }

      const int pid = getPidConfigurable(track);
      if (pid == -1 && withinPtRef) {
        nUnidentifiedWeighted += weff;
        unidentifiedPtSum += weff * track.pt();
        continue;
      }
      if (!isWithinPOIPtRange(pid, track.pt())) {
        continue;
      }

      float weffPid = 1.f;
      setParticleNUEWeight(weffPid, track, cent, pid);
      if (pid == MyParticleType::kPion) {
        nPionWeighted += weffPid;
        pionPtSum += weffPid * track.pt();
      } else if (pid == MyParticleType::kKaon) {
        nKaonWeighted += weffPid;
        kaonPtSum += weffPid * track.pt();
      } else if (pid == MyParticleType::kProton) {
        nProtonWeighted += weffPid;
        protonPtSum += weffPid * track.pt();
      }
    }
    if (nch <= 0.) {
      return;
    }

    // Fill the charged reference and particle-specific POI regions in a separate pass.
    for (const auto& track : tracks) {
      if (!trackSelectedForFlow(track)) {
        continue;
      }

      float weff = 1.f;
      float wacc = 1.f;
      setParticleNUEWeight(weff, track, cent);
      setParticleNUAWeight(wacc, track, vtxz);

      const bool withinPtRef = isWithinRefPtRange(track.pt());
      if (withinPtRef) {
        /// @note passed selection: global cut + pt cut
        fGFW->Fill(track.eta(), 0, track.phi(), wacc * weff, gfwBitMask::charged::kReference);
      }

      const int pid = getPidConfigurable(track);
      if (pid == -1 && withinPtRef) {
        fGFW->Fill(track.eta(), 0, track.phi(), wacc * weff, gfwBitMask::unidentified::kPOI);
        continue;
      }
      if (!isWithinPOIPtRange(pid, track.pt())) {
        continue;
      }

      float weffPid = 1.f;
      float waccPid = 1.f;
      setParticleNUEWeight(weffPid, track, cent, pid);
      setParticleNUAWeight(waccPid, track, vtxz, pid);

      if (pid == MyParticleType::kPion) {
        fGFW->Fill(track.eta(), 0, track.phi(), waccPid * weffPid, gfwBitMask::pion::kPOI);
      } else if (pid == MyParticleType::kKaon) {
        fGFW->Fill(track.eta(), 0, track.phi(), waccPid * weffPid, gfwBitMask::kaon::kPOI);
      } else if (pid == MyParticleType::kProton) {
        fGFW->Fill(track.eta(), 0, track.phi(), waccPid * weffPid, gfwBitMask::proton::kPOI);
      }
    }

    // All profiles from one event use the same bootstrap subsample. Mean-pT profiles
    // retain the NUE-corrected particle-count weight used by FlowContainer hMeanPt.
    if (switchsOpts.cfgOutputCharged.value) {
      registry.fill(HIST("c22DeltaPt/hMeanPtCharged"), cent, rndm * cfgFlowNbootstrap, ptSum / nch, nch);
    }
    if (nPionWeighted > 0. && switchsOpts.cfgOutputPion.value) {
      registry.fill(HIST("c22DeltaPt/hMeanPtPion"), cent, rndm * cfgFlowNbootstrap, pionPtSum / nPionWeighted, nPionWeighted);
      fillProfilePOIvnpt(corrconfigs.at(corrConfigIndex::charged::kC22Gap), HIST("c22DeltaPt/c22dmeanptPionRefRef"), cent, pionPtSum, nPionWeighted, rndm);
    }
    if (nKaonWeighted > 0. && switchsOpts.cfgOutputKaon.value) {
      registry.fill(HIST("c22DeltaPt/hMeanPtKaon"), cent, rndm * cfgFlowNbootstrap, kaonPtSum / nKaonWeighted, nKaonWeighted);
      fillProfilePOIvnpt(corrconfigs.at(corrConfigIndex::charged::kC22Gap), HIST("c22DeltaPt/c22dmeanptKaonRefRef"), cent, kaonPtSum, nKaonWeighted, rndm);
    }
    if (nProtonWeighted > 0. && switchsOpts.cfgOutputProton.value) {
      registry.fill(HIST("c22DeltaPt/hMeanPtProton"), cent, rndm * cfgFlowNbootstrap, protonPtSum / nProtonWeighted, nProtonWeighted);
      fillProfilePOIvnpt(corrconfigs.at(corrConfigIndex::charged::kC22Gap), HIST("c22DeltaPt/c22dmeanptProtonRefRef"), cent, protonPtSum, nProtonWeighted, rndm);
    }
    if (nUnidentifiedWeighted > 0. && switchsOpts.cfgOutputUnidentified.value) {
      registry.fill(HIST("c22DeltaPt/hMeanPtUnidentified"), cent, rndm * cfgFlowNbootstrap, unidentifiedPtSum / nUnidentifiedWeighted, nUnidentifiedWeighted);
      fillProfilePOIvnpt(corrconfigs.at(corrConfigIndex::charged::kC22Gap), HIST("c22DeltaPt/c22dmeanptUnidentifiedRefRef"), cent, unidentifiedPtSum, nUnidentifiedWeighted, rndm);
    }

    // Charged and PID-binned ref-ref are common to both modes. The numerator is
    // selected between PID Pure and POI-ref by cfgC22DeltaPtUsePure.
    if (switchsOpts.cfgOutputCharged.value) {
      fillProfilePOIvnpt(corrconfigs.at(corrConfigIndex::charged::kC22Gap), HIST("c22DeltaPt/c22dmeanptCharged"), cent, ptSum, nch, rndm);
    }
    if (switchsOpts.cfgC22DeltaPtUsePure.value) {
      if (nPionWeighted > 0. && switchsOpts.cfgOutputPion.value) {
        fillProfilePOIvnpt(corrconfigs.at(corrConfigIndex::pion::kC22Pure), HIST("c22DeltaPt/c22dmeanptPionPure"), cent, pionPtSum, nPionWeighted, rndm);
      }
      if (nKaonWeighted > 0. && switchsOpts.cfgOutputKaon.value) {
        fillProfilePOIvnpt(corrconfigs.at(corrConfigIndex::kaon::kC22Pure), HIST("c22DeltaPt/c22dmeanptKaonPure"), cent, kaonPtSum, nKaonWeighted, rndm);
      }
      if (nProtonWeighted > 0. && switchsOpts.cfgOutputProton.value) {
        fillProfilePOIvnpt(corrconfigs.at(corrConfigIndex::proton::kC22Pure), HIST("c22DeltaPt/c22dmeanptProtonPure"), cent, protonPtSum, nProtonWeighted, rndm);
      }
      if (nUnidentifiedWeighted > 0. && switchsOpts.cfgOutputUnidentified.value) {
        fillProfilePOIvnpt(corrconfigs.at(corrConfigIndex::unidentified::kC22Pure), HIST("c22DeltaPt/c22dmeanptUnidentifiedPure"), cent, unidentifiedPtSum, nUnidentifiedWeighted, rndm);
      }
    } else {
      if (nPionWeighted > 0. && switchsOpts.cfgOutputPion.value) {
        fillProfilePOIvnpt(corrconfigs.at(corrConfigIndex::pion::kC22GapA), corrconfigs.at(corrConfigIndex::pion::kC22GapB), HIST("c22DeltaPt/c22dmeanptPionPOIRef"), cent, pionPtSum, nPionWeighted, rndm);
      }
      if (nKaonWeighted > 0. && switchsOpts.cfgOutputKaon.value) {
        fillProfilePOIvnpt(corrconfigs.at(corrConfigIndex::kaon::kC22GapA), corrconfigs.at(corrConfigIndex::kaon::kC22GapB), HIST("c22DeltaPt/c22dmeanptKaonPOIRef"), cent, kaonPtSum, nKaonWeighted, rndm);
      }
      if (nProtonWeighted > 0. && switchsOpts.cfgOutputProton.value) {
        fillProfilePOIvnpt(corrconfigs.at(corrConfigIndex::proton::kC22GapA), corrconfigs.at(corrConfigIndex::proton::kC22GapB), HIST("c22DeltaPt/c22dmeanptProtonPOIRef"), cent, protonPtSum, nProtonWeighted, rndm);
      }
      if (nUnidentifiedWeighted > 0. && switchsOpts.cfgOutputUnidentified.value) {
        fillProfilePOIvnpt(corrconfigs.at(corrConfigIndex::unidentified::kC22GapA), corrconfigs.at(corrConfigIndex::unidentified::kC22GapB), HIST("c22DeltaPt/c22dmeanptUnidentifiedPOIRef"), cent, unidentifiedPtSum, nUnidentifiedWeighted, rndm);
      }
    }
  }
  PROCESS_SWITCH(PidFlowPtCorr, processDataC22DeltaPt, "Fill c22dmeanpt profiles", false);

  void fillPidPtSelfMoments(const char* productName,
                            const char* singlePtName,
                            const double cent,
                            const double rndm,
                            const double sumWeight,
                            const double sumPt,
                            const double sumWeight2,
                            const double sumPtWeight2,
                            const double sumPt2Weight2)
  {
    const double pairWeight = sumWeight * sumWeight - sumWeight2;
    if (pairWeight <= minVal4Float) {
      return;
    }
    fFCPidPtCorr->FillProfile(productName, cent,
                              (sumPt * sumPt - sumPt2Weight2) / pairWeight,
                              pairWeight, rndm);
    fFCPidPtCorr->FillProfile(singlePtName, cent,
                              (sumWeight * sumPt - sumPtWeight2) / pairWeight,
                              pairWeight, rndm);
  }

  void fillPidPtCrossMoments(const char* productName,
                             const char* alphaPtName,
                             const char* betaPtName,
                             const double cent,
                             const double rndm,
                             const double alphaWeight,
                             const double alphaPtSum,
                             const double betaWeight,
                             const double betaPtSum)
  {
    if (alphaWeight <= 0. || betaWeight <= 0.) {
      return;
    }
    const double pairWeight = alphaWeight * betaWeight;
    fFCPidPtCorr->FillProfile(productName, cent,
                              alphaPtSum * betaPtSum / pairWeight,
                              pairWeight, rndm);
    fFCPidPtCorr->FillProfile(alphaPtName, cent,
                              alphaPtSum / alphaWeight,
                              pairWeight, rndm);
    fFCPidPtCorr->FillProfile(betaPtName, cent,
                              betaPtSum / betaWeight,
                              pairWeight, rndm);
  }

  void processPidPtCorrelations(AodCollisions::iterator const& collision,
                                aod::BCsWithTimestamps const&,
                                AodTracks const& tracks)
  {
    if (!switchsOpts.cfgAddPidPtCorrelationProfiles.value) {
      return;
    }
    registry.fill(HIST("hEventCount/processPidPtCorrelations"), 0.5);
    if (tracks.size() < 1 || !collision.sel8()) {
      return;
    }
    registry.fill(HIST("hEventCount/processPidPtCorrelations"), 1.5);

    const auto cent = getCentrality(collision);
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    const double interactionRate = getInteractionRate(bc.timestamp(), bc.runNumber());
    if (!eventSelected(collision, cent, interactionRate, MyFunctionName::funcProcessPidPtCorrelations)) {
      return;
    }

    loadCorrections(bc.timestamp());
    const double rndm = fRndm->Rndm();
    std::array<double, 3> sumWeight{};
    std::array<double, 3> sumPt{};
    std::array<double, 3> sumWeight2{};
    std::array<double, 3> sumPtWeight2{};
    std::array<double, 3> sumPt2Weight2{};

    for (const auto& track : tracks) {
      if (!trackSelectedForFlow(track)) {
        continue;
      }
      const int pid = getPidConfigurable(track);
      if (pid < MyParticleType::kPion || pid > MyParticleType::kProton || !isWithinPOIPtRange(pid, track.pt())) {
        continue;
      }

      float weight = 1.f;
      setParticleNUEWeight(weight, track, cent, pid);
      const auto index = static_cast<std::size_t>(pid - MyParticleType::kPion);
      const double particleWeight = weight;
      const double weight2 = particleWeight * particleWeight;
      const double pt = track.pt();
      sumWeight[index] += particleWeight;
      sumPt[index] += particleWeight * pt;
      sumWeight2[index] += weight2;
      sumPtWeight2[index] += weight2 * pt;
      sumPt2Weight2[index] += weight2 * pt * pt;
    }

    constexpr std::size_t kPi = 0;
    constexpr std::size_t kKa = 1;
    constexpr std::size_t kPr = 2;
    if (sumWeight[kPi] > 0.) {
      fFCPidPtCorr->FillProfile("meanPtPi", cent, sumPt[kPi] / sumWeight[kPi], sumWeight[kPi], rndm);
    }
    if (sumWeight[kKa] > 0.) {
      fFCPidPtCorr->FillProfile("meanPtKa", cent, sumPt[kKa] / sumWeight[kKa], sumWeight[kKa], rndm);
    }
    if (sumWeight[kPr] > 0.) {
      fFCPidPtCorr->FillProfile("meanPtPr", cent, sumPt[kPr] / sumWeight[kPr], sumWeight[kPr], rndm);
    }

    fillPidPtSelfMoments("ptProductPiPi", "ptPiInPiPiPairs", cent, rndm,
                         sumWeight[kPi], sumPt[kPi], sumWeight2[kPi], sumPtWeight2[kPi], sumPt2Weight2[kPi]);
    fillPidPtSelfMoments("ptProductKaKa", "ptKaInKaKaPairs", cent, rndm,
                         sumWeight[kKa], sumPt[kKa], sumWeight2[kKa], sumPtWeight2[kKa], sumPt2Weight2[kKa]);
    fillPidPtSelfMoments("ptProductPrPr", "ptPrInPrPrPairs", cent, rndm,
                         sumWeight[kPr], sumPt[kPr], sumWeight2[kPr], sumPtWeight2[kPr], sumPt2Weight2[kPr]);

    fillPidPtCrossMoments("ptProductPiKa", "ptPiInPiKaPairs", "ptKaInPiKaPairs", cent, rndm,
                          sumWeight[kPi], sumPt[kPi], sumWeight[kKa], sumPt[kKa]);
    fillPidPtCrossMoments("ptProductPiPr", "ptPiInPiPrPairs", "ptPrInPiPrPairs", cent, rndm,
                          sumWeight[kPi], sumPt[kPi], sumWeight[kPr], sumPt[kPr]);
    fillPidPtCrossMoments("ptProductKaPr", "ptKaInKaPrPairs", "ptPrInKaPrPairs", cent, rndm,
                          sumWeight[kKa], sumPt[kKa], sumWeight[kPr], sumPt[kPr]);
  }
  PROCESS_SWITCH(PidFlowPtCorr, processPidPtCorrelations, "Calculate self-contained PID mean-pT correlations", false);

  /**
   * @brief Run the flow calculation on generated MC particles for a closure test.
   * @note Reconstructed collisions only provide centrality or raw PV multiplicity. The selection is
   *       deliberately truth-level: MC vertex filter, physical-primary particles, kinematics,
   *       and PDG PID. No data event, track-quality, detector-PID, NUA, or NUE cuts are applied.
   */
  void processMCClosure(FilteredMcCollisions::iterator const& mcCollision,
                        aod::BCsWithTimestamps const&,
                        soa::SmallGroups<soa::Join<aod::McCollisionLabels, AodCollisions>> const& collisions,
                        FilteredMcParticles const& mcParticles,
                        FilteredTracksWithMCLabel const&)
  {
    registry.fill(HIST("hEventCount/processData"), 0.5);
    // A generated event must be associated with exactly one reconstructed collision.
    // Otherwise the same truth particles would be filled more than once.
    if (collisions.size() != 1 || mcParticles.size() <= 0) {
      return;
    }

    for (const auto& collision : collisions) {
      const float cent = getCentrality(collision);
      const float vtxz = mcCollision.posZ();
      const float rndm = fRndm->Rndm();
      const int nTot = mcParticles.size();

      fGFW->Clear();
      registry.fill(HIST("hEventCount/processData"), 1.5);
      registry.fill(HIST("hVtxZ"), vtxz);
      registry.fill(HIST("hMult"), nTot);
      registry.fill(HIST("hCent"), cent);
      if (switchsOpts.cfgOutPutPtSpectra.value) {
        registry.fill(HIST("ptSpectra/hCentEventCountData"), cent);
      }

      double ptSum = 0., ptSumw2 = 0., nch = 0., nchSquare = 0., ptSquareSum = 0.;
      double pionPtSum = 0., kaonPtSum = 0., protonPtSum = 0.;
      double nPionWeighted = 0., nKaonWeighted = 0., nProtonWeighted = 0.;
      double pionPtSumw2 = 0., kaonPtSumw2 = 0., protonPtSumw2 = 0.;
      double nPionSquare = 0., nKaonSquare = 0., nProtonSquare = 0.;
      double pionPtSquareSum = 0., kaonPtSquareSum = 0., protonPtSquareSum = 0.;
      int numOfPi = 0, numOfKa = 0, numOfPr = 0;

      for (const auto& mcParticle : mcParticles) {
        if (!mcParticle.isPhysicalPrimary() || !isStable(mcParticle.pdgCode()) || std::fabs(mcParticle.eta()) > trkQualityOpts.cfgCutEta.value) {
          continue;
        }

        const float pt = mcParticle.pt();
        const bool withinPtRef = pt > trkQualityOpts.cfgCutPtMin.value && pt < trkQualityOpts.cfgCutPtMax.value;
        const bool withinMeanPtEta = std::fabs(mcParticle.eta()) < trkQualityOpts.cfgRangeEta.value;
        const int absPdg = std::abs(mcParticle.pdgCode());
        int pid = -1;
        if (absPdg == PDG_t::kPiPlus) {
          pid = MyParticleType::kPion;
        } else if (absPdg == PDG_t::kKPlus) {
          pid = MyParticleType::kKaon;
        } else if (absPdg == PDG_t::kProton) {
          pid = MyParticleType::kProton;
        }

        if (withinPtRef && withinMeanPtEta) {
          nch += 1.;
          nchSquare += 1.;
          ptSum += pt;
          ptSumw2 += pt;
          ptSquareSum += pt * pt;
        }

        if (withinMeanPtEta && isWithinPOIPtRange(pid, pt)) {
          if (pid == MyParticleType::kPion) {
            nPionWeighted += 1.;
            nPionSquare += 1.;
            pionPtSum += pt;
            pionPtSumw2 += pt;
            pionPtSquareSum += pt * pt;
          } else if (pid == MyParticleType::kKaon) {
            nKaonWeighted += 1.;
            nKaonSquare += 1.;
            kaonPtSum += pt;
            kaonPtSumw2 += pt;
            kaonPtSquareSum += pt * pt;
          } else if (pid == MyParticleType::kProton) {
            nProtonWeighted += 1.;
            nProtonSquare += 1.;
            protonPtSum += pt;
            protonPtSumw2 += pt;
            protonPtSquareSum += pt * pt;
          }
        }

        if (withinPtRef) {
          registry.fill(HIST("hPhi"), mcParticle.phi());
          registry.fill(HIST("hPhicorr"), mcParticle.phi());
          registry.fill(HIST("hPhiCorrNUANUE"), mcParticle.phi());
          registry.fill(HIST("hEta"), mcParticle.eta());
          registry.fill(HIST("hPt"), pt);
          registry.fill(HIST("hPtCorr"), pt);
          if (switchsOpts.cfgOutPutPtSpectra.value) {
            registry.fill(HIST("ptSpectra/hPtCentData"), pt, cent);
          }
          fGFW->Fill(mcParticle.eta(), 0, mcParticle.phi(), 1., gfwBitMask::charged::kReference);
        }

        if (!isWithinPOIPtRange(pid, pt)) {
          continue;
        }
        if (pid == MyParticleType::kPion) {
          fGFW->Fill(mcParticle.eta(), 0, mcParticle.phi(), 1., gfwBitMask::pion::kPOI);
          fGFW->Fill(mcParticle.eta(), 0, mcParticle.phi(), 1., gfwBitMask::pion::kOverlap, 1.);
          registry.fill(HIST("hPhiPi"), mcParticle.phi());
          registry.fill(HIST("hPhicorrPi"), mcParticle.phi());
          registry.fill(HIST("hPhiCorrNUANUEPi"), mcParticle.phi());
          registry.fill(HIST("hPtPi"), pt);
          registry.fill(HIST("hPtCorrPi"), pt);
          ++numOfPi;
        } else if (pid == MyParticleType::kKaon) {
          fGFW->Fill(mcParticle.eta(), 0, mcParticle.phi(), 1., gfwBitMask::kaon::kPOI);
          fGFW->Fill(mcParticle.eta(), 0, mcParticle.phi(), 1., gfwBitMask::kaon::kOverlap, 1.);
          registry.fill(HIST("hPhiKa"), mcParticle.phi());
          registry.fill(HIST("hPhicorrKa"), mcParticle.phi());
          registry.fill(HIST("hPhiCorrNUANUEKa"), mcParticle.phi());
          registry.fill(HIST("hPtKa"), pt);
          registry.fill(HIST("hPtCorrKa"), pt);
          ++numOfKa;
        } else if (pid == MyParticleType::kProton) {
          fGFW->Fill(mcParticle.eta(), 0, mcParticle.phi(), 1., gfwBitMask::proton::kPOI);
          fGFW->Fill(mcParticle.eta(), 0, mcParticle.phi(), 1., gfwBitMask::proton::kOverlap, 1.);
          registry.fill(HIST("hPhiPr"), mcParticle.phi());
          registry.fill(HIST("hPhicorrPr"), mcParticle.phi());
          registry.fill(HIST("hPhiCorrNUANUEPr"), mcParticle.phi());
          registry.fill(HIST("hPtPr"), pt);
          registry.fill(HIST("hPtCorrPr"), pt);
          ++numOfPr;
        }
      }

      if (particleAbundanceOpts.cfgOutPutAbundanceDis) {
        registry.fill(HIST("abundance/hNumOfPiEventCount"), numOfPi);
        registry.fill(HIST("abundance/hNumOfKaEventCount"), numOfKa);
        registry.fill(HIST("abundance/hNumOfPrEventCount"), numOfPr);
      }
      registry.fill(HIST("hNchUnCorrectedVSNchCorrected"), nch, nch);

      if (nch <= 0.) {
        continue;
      }

      fillFC(MyParticleType::kCharged, corrconfigs.at(corrConfigIndex::charged::kC22Gap), cent, rndm, "c22");
      fillFC(MyParticleType::kCharged, corrconfigs.at(corrConfigIndex::charged::kC24), cent, rndm, "c24");
      fillFC(MyParticleType::kCharged, corrconfigs.at(corrConfigIndex::charged::kC22Full), cent, rndm, "c22Full");
      fillFC(MyParticleType::kCharged, corrconfigs.at(corrConfigIndex::charged::kC32Gap), cent, rndm, "c32");
      fillFC(MyParticleType::kCharged, corrconfigs.at(corrConfigIndex::charged::kC34Gap), cent, rndm, "c34");

      fillFC(MyParticleType::kPion, corrconfigs.at(corrConfigIndex::pion::kC22FullA), cent, rndm, "c22Full");
      fillFC(MyParticleType::kPion, corrconfigs.at(corrConfigIndex::pion::kC22FullB), cent, rndm, "c22Full");
      fillFC(MyParticleType::kKaon, corrconfigs.at(corrConfigIndex::kaon::kC22FullA), cent, rndm, "c22Full");
      fillFC(MyParticleType::kKaon, corrconfigs.at(corrConfigIndex::kaon::kC22FullB), cent, rndm, "c22Full");
      fillFC(MyParticleType::kProton, corrconfigs.at(corrConfigIndex::proton::kC22FullA), cent, rndm, "c22Full");
      fillFC(MyParticleType::kProton, corrconfigs.at(corrConfigIndex::proton::kC22FullB), cent, rndm, "c22Full");

      fillFC(MyParticleType::kPion, corrconfigs.at(corrConfigIndex::pion::kC24A), cent, rndm, "c24");
      fillFC(MyParticleType::kPion, corrconfigs.at(corrConfigIndex::pion::kC24B), cent, rndm, "c24");
      fillFC(MyParticleType::kKaon, corrconfigs.at(corrConfigIndex::kaon::kC24A), cent, rndm, "c24");
      fillFC(MyParticleType::kKaon, corrconfigs.at(corrConfigIndex::kaon::kC24B), cent, rndm, "c24");
      fillFC(MyParticleType::kProton, corrconfigs.at(corrConfigIndex::proton::kC24A), cent, rndm, "c24");
      fillFC(MyParticleType::kProton, corrconfigs.at(corrConfigIndex::proton::kC24B), cent, rndm, "c24");

      fillFC(MyParticleType::kPion, corrconfigs.at(corrConfigIndex::pion::kC32GapA), cent, rndm, "c32");
      fillFC(MyParticleType::kPion, corrconfigs.at(corrConfigIndex::pion::kC32GapB), cent, rndm, "c32");
      fillFC(MyParticleType::kKaon, corrconfigs.at(corrConfigIndex::kaon::kC32GapA), cent, rndm, "c32");
      fillFC(MyParticleType::kKaon, corrconfigs.at(corrConfigIndex::kaon::kC32GapB), cent, rndm, "c32");
      fillFC(MyParticleType::kProton, corrconfigs.at(corrConfigIndex::proton::kC32GapA), cent, rndm, "c32");
      fillFC(MyParticleType::kProton, corrconfigs.at(corrConfigIndex::proton::kC32GapB), cent, rndm, "c32");

      fillFC(MyParticleType::kPion, corrconfigs.at(corrConfigIndex::pion::kC34A), cent, rndm, "c34");
      fillFC(MyParticleType::kPion, corrconfigs.at(corrConfigIndex::pion::kC34B), cent, rndm, "c34");
      fillFC(MyParticleType::kKaon, corrconfigs.at(corrConfigIndex::kaon::kC34A), cent, rndm, "c34");
      fillFC(MyParticleType::kKaon, corrconfigs.at(corrConfigIndex::kaon::kC34B), cent, rndm, "c34");
      fillFC(MyParticleType::kProton, corrconfigs.at(corrConfigIndex::proton::kC34A), cent, rndm, "c34");
      fillFC(MyParticleType::kProton, corrconfigs.at(corrConfigIndex::proton::kC34B), cent, rndm, "c34");

      const bool filledPi = fillFC(MyParticleType::kPion, corrconfigs.at(corrConfigIndex::pion::kC22Pure), cent, rndm, "c22pure");
      const bool filledKa = fillFC(MyParticleType::kKaon, corrconfigs.at(corrConfigIndex::kaon::kC22Pure), cent, rndm, "c22pure");
      const bool filledPr = fillFC(MyParticleType::kProton, corrconfigs.at(corrConfigIndex::proton::kC22Pure), cent, rndm, "c22pure");
      fillFC(MyParticleType::kPion, corrconfigs.at(corrConfigIndex::pion::kC32Pure), cent, rndm, "c32pure");
      fillFC(MyParticleType::kKaon, corrconfigs.at(corrConfigIndex::kaon::kC32Pure), cent, rndm, "c32pure");
      fillFC(MyParticleType::kProton, corrconfigs.at(corrConfigIndex::proton::kC32Pure), cent, rndm, "c32pure");

      if (filledPi || !switchsOpts.cfgCheck2MethodDiff.value) {
        fillFC(MyParticleType::kPion, corrconfigs.at(corrConfigIndex::pion::kC22GapA), cent, rndm, "c22");
        fillFC(MyParticleType::kPion, corrconfigs.at(corrConfigIndex::pion::kC22GapB), cent, rndm, "c22");
      }
      if (filledKa || !switchsOpts.cfgCheck2MethodDiff.value) {
        fillFC(MyParticleType::kKaon, corrconfigs.at(corrConfigIndex::kaon::kC22GapA), cent, rndm, "c22");
        fillFC(MyParticleType::kKaon, corrconfigs.at(corrConfigIndex::kaon::kC22GapB), cent, rndm, "c22");
      }
      if (filledPr || !switchsOpts.cfgCheck2MethodDiff.value) {
        fillFC(MyParticleType::kProton, corrconfigs.at(corrConfigIndex::proton::kC22GapA), cent, rndm, "c22");
        fillFC(MyParticleType::kProton, corrconfigs.at(corrConfigIndex::proton::kC22GapB), cent, rndm, "c22");
      }

      fillFCvnpt(MyParticleType::kCharged, corrconfigs.at(corrConfigIndex::charged::kC22Gap), cent, rndm, nch, nch, "c22TrackWeight");
      fillFCvnpt(MyParticleType::kCharged, corrconfigs.at(corrConfigIndex::charged::kC22Gap), cent, rndm, nch, nch, "c22TrackWeightOne", true);
      fillFCvnpt(MyParticleType::kCharged, corrconfigs.at(corrConfigIndex::charged::kC24), cent, rndm, nch, nch, "c24TrackWeight");
      fillFCvnpt(MyParticleType::kCharged, corrconfigs.at(corrConfigIndex::charged::kC22Full), cent, rndm, nch, nch, "c22FullTrackWeight");
      fillFCvnpt(MyParticleType::kCharged, corrconfigs.at(corrConfigIndex::charged::kC32Gap), cent, rndm, nch, nch, "c32TrackWeight");
      fillFCvnpt(MyParticleType::kCharged, corrconfigs.at(corrConfigIndex::charged::kC34Gap), cent, rndm, nch, nch, "c34TrackWeight");
      fillFCvnpt(MyParticleType::kPion, corrconfigs.at(corrConfigIndex::pion::kC22Pure), cent, rndm, nch, nch, "c22TrackWeight");
      fillFCvnpt(MyParticleType::kKaon, corrconfigs.at(corrConfigIndex::kaon::kC22Pure), cent, rndm, nch, nch, "c22TrackWeight");
      fillFCvnpt(MyParticleType::kProton, corrconfigs.at(corrConfigIndex::proton::kC22Pure), cent, rndm, nch, nch, "c22TrackWeight");
      fillFCvnpt(MyParticleType::kPion, corrconfigs.at(corrConfigIndex::pion::kC32Pure), cent, rndm, nch, nch, "c32TrackWeight");
      fillFCvnpt(MyParticleType::kKaon, corrconfigs.at(corrConfigIndex::kaon::kC32Pure), cent, rndm, nch, nch, "c32TrackWeight");
      fillFCvnpt(MyParticleType::kProton, corrconfigs.at(corrConfigIndex::proton::kC32Pure), cent, rndm, nch, nch, "c32TrackWeight");
      fillFCvnpt(MyParticleType::kCharged, corrconfigs.at(corrConfigIndex::charged::kC22Gap), cent, rndm, ptSum, nch, "covV2Pt");
      fillFCvnpt(MyParticleType::kCharged, corrconfigs.at(corrConfigIndex::charged::kC22Gap), cent, rndm, ptSum, nch, "covV2PtWeightOne", true);
      fillFCvnpt(MyParticleType::kPion, corrconfigs.at(corrConfigIndex::pion::kC22Pure), cent, rndm, ptSum, nch, "covV2Pt");
      fillFCvnpt(MyParticleType::kKaon, corrconfigs.at(corrConfigIndex::kaon::kC22Pure), cent, rndm, ptSum, nch, "covV2Pt");
      fillFCvnpt(MyParticleType::kProton, corrconfigs.at(corrConfigIndex::proton::kC22Pure), cent, rndm, ptSum, nch, "covV2Pt");
      fillFCvnpt(MyParticleType::kCharged, corrconfigs.at(corrConfigIndex::charged::kC32Gap), cent, rndm, ptSum, nch, "covV3Pt");
      fillFCvnpt(MyParticleType::kPion, corrconfigs.at(corrConfigIndex::pion::kC32Pure), cent, rndm, ptSum, nch, "covV3Pt");
      fillFCvnpt(MyParticleType::kKaon, corrconfigs.at(corrConfigIndex::kaon::kC32Pure), cent, rndm, ptSum, nch, "covV3Pt");
      fillFCvnpt(MyParticleType::kProton, corrconfigs.at(corrConfigIndex::proton::kC32Pure), cent, rndm, ptSum, nch, "covV3Pt");

      fillFC4PtC22(cent, ptSum, nch, rndm);
      if (nPionWeighted > 0.) {
        fillFC4PtC22(cent, rndm, MyParticleType::kPion, pionPtSum, nPionWeighted);
      }
      if (nKaonWeighted > 0.) {
        fillFC4PtC22(cent, rndm, MyParticleType::kKaon, kaonPtSum, nKaonWeighted);
      }
      if (nProtonWeighted > 0.) {
        fillFC4PtC22(cent, rndm, MyParticleType::kProton, protonPtSum, nProtonWeighted);
      }

      if (switchsOpts.cfgOutputCharged.value) {
        fFCCh->FillProfile("hMeanPt", cent, ptSum / nch, nch, rndm);
        fFCCh->FillProfile("hMeanPtWeightOne", cent, ptSum / nch, 1., rndm);
      }
      if (nPionWeighted > 0. && switchsOpts.cfgOutputPion.value) {
        fFCPi->FillProfile("hMeanPt", cent, pionPtSum / nPionWeighted, nPionWeighted, rndm);
        fFCPi->FillProfile("hMeanPtWeightOne", cent, pionPtSum / nPionWeighted, 1., rndm);
        fillFCvnpt(MyParticleType::kPion, corrconfigs.at(corrConfigIndex::pion::kC22Pure), cent, rndm, pionPtSum, nPionWeighted, "covV2PtPID");
        fillFCvnpt(MyParticleType::kPion, corrconfigs.at(corrConfigIndex::pion::kC22Pure), cent, rndm, nPionWeighted, nPionWeighted, "c22TrackWeightPID");
      }
      if (nKaonWeighted > 0. && switchsOpts.cfgOutputKaon.value) {
        fFCKa->FillProfile("hMeanPt", cent, kaonPtSum / nKaonWeighted, nKaonWeighted, rndm);
        fFCKa->FillProfile("hMeanPtWeightOne", cent, kaonPtSum / nKaonWeighted, 1., rndm);
        fillFCvnpt(MyParticleType::kKaon, corrconfigs.at(corrConfigIndex::kaon::kC22Pure), cent, rndm, kaonPtSum, nKaonWeighted, "covV2PtPID");
        fillFCvnpt(MyParticleType::kKaon, corrconfigs.at(corrConfigIndex::kaon::kC22Pure), cent, rndm, nKaonWeighted, nKaonWeighted, "c22TrackWeightPID");
      }
      if (nProtonWeighted > 0. && switchsOpts.cfgOutputProton.value) {
        fFCPr->FillProfile("hMeanPt", cent, protonPtSum / nProtonWeighted, nProtonWeighted, rndm);
        fFCPr->FillProfile("hMeanPtWeightOne", cent, protonPtSum / nProtonWeighted, 1., rndm);
        fillFCvnpt(MyParticleType::kProton, corrconfigs.at(corrConfigIndex::proton::kC22Pure), cent, rndm, protonPtSum, nProtonWeighted, "covV2PtPID");
        fillFCvnpt(MyParticleType::kProton, corrconfigs.at(corrConfigIndex::proton::kC22Pure), cent, rndm, nProtonWeighted, nProtonWeighted, "c22TrackWeightPID");
      }

      if (switchsOpts.cfgOutPutPtSpectra.value) {
        const double nPairCharged = fGFW->Calculate(corrconfigs.at(corrConfigIndex::charged::kC22Gap), 0, true).real();
        const double chargedC22 = nPairCharged > 0. ? fGFW->Calculate(corrconfigs.at(corrConfigIndex::charged::kC22Gap), 0, false).real() / nPairCharged : 0.;
        const double pidChargedC22Pi = getPidC22InOneEvent(corrconfigs.at(corrConfigIndex::pion::kC22GapA), corrconfigs.at(corrConfigIndex::pion::kC22GapB));
        const double pidKaonC22 = getPidC22InOneEvent(corrconfigs.at(corrConfigIndex::kaon::kC22GapA), corrconfigs.at(corrConfigIndex::kaon::kC22GapB));
        const double pidProtonC22 = getPidC22InOneEvent(corrconfigs.at(corrConfigIndex::proton::kC22GapA), corrconfigs.at(corrConfigIndex::proton::kC22GapB));
        if (pidChargedC22Pi > 0. && chargedC22 > 0.) {
          registry.fill(HIST("c22PrimeVsc22/Pi"), pidChargedC22Pi, chargedC22);
        }
        if (pidKaonC22 > 0. && chargedC22 > 0.) {
          registry.fill(HIST("c22PrimeVsc22/Ka"), pidKaonC22, chargedC22);
        }
        if (pidProtonC22 > 0. && chargedC22 > 0.) {
          registry.fill(HIST("c22PrimeVsc22/Pr"), pidProtonC22, chargedC22);
        }
      }

      const double nchDiff = nch * nch - nchSquare;
      if (nchDiff > minVal4Float && switchsOpts.cfgOutputCharged.value) {
        fFCCh->FillProfile("ptSquareAve", cent, (ptSum * ptSum - ptSquareSum) / nchDiff, nchDiff, rndm);
        fFCCh->FillProfile("ptSquareAveWeightOne", cent, (ptSum * ptSum - ptSquareSum) / nchDiff, 1., rndm);
        fFCCh->FillProfile("ptAve", cent, (nch * ptSum - ptSumw2) / nchDiff, nchDiff, rndm);
        fFCCh->FillProfile("ptAveWeightOne", cent, (nch * ptSum - ptSumw2) / nchDiff, 1., rndm);
      }
      const double pionDiff = nPionWeighted * nPionWeighted - nPionSquare;
      if (pionDiff > minVal4Float && switchsOpts.cfgOutputPion.value) {
        fFCPi->FillProfile("ptSquareAve", cent, (pionPtSum * pionPtSum - pionPtSquareSum) / pionDiff, pionDiff, rndm);
        fFCPi->FillProfile("ptSquareAveWeightOne", cent, (pionPtSum * pionPtSum - pionPtSquareSum) / pionDiff, 1., rndm);
        fFCPi->FillProfile("ptAve", cent, (nPionWeighted * pionPtSum - pionPtSumw2) / pionDiff, pionDiff, rndm);
        fFCPi->FillProfile("ptAveWeightOne", cent, (nPionWeighted * pionPtSum - pionPtSumw2) / pionDiff, 1., rndm);
      }
      const double kaonDiff = nKaonWeighted * nKaonWeighted - nKaonSquare;
      if (kaonDiff > minVal4Float && switchsOpts.cfgOutputKaon.value) {
        fFCKa->FillProfile("ptSquareAve", cent, (kaonPtSum * kaonPtSum - kaonPtSquareSum) / kaonDiff, kaonDiff, rndm);
        fFCKa->FillProfile("ptSquareAveWeightOne", cent, (kaonPtSum * kaonPtSum - kaonPtSquareSum) / kaonDiff, 1., rndm);
        fFCKa->FillProfile("ptAve", cent, (nKaonWeighted * kaonPtSum - kaonPtSumw2) / kaonDiff, kaonDiff, rndm);
        fFCKa->FillProfile("ptAveWeightOne", cent, (nKaonWeighted * kaonPtSum - kaonPtSumw2) / kaonDiff, 1., rndm);
      }
      const double protonDiff = nProtonWeighted * nProtonWeighted - nProtonSquare;
      if (protonDiff > minVal4Float && switchsOpts.cfgOutputProton.value) {
        fFCPr->FillProfile("ptSquareAve", cent, (protonPtSum * protonPtSum - protonPtSquareSum) / protonDiff, protonDiff, rndm);
        fFCPr->FillProfile("ptSquareAveWeightOne", cent, (protonPtSum * protonPtSum - protonPtSquareSum) / protonDiff, 1., rndm);
        fFCPr->FillProfile("ptAve", cent, (nProtonWeighted * protonPtSum - protonPtSumw2) / protonDiff, protonDiff, rndm);
        fFCPr->FillProfile("ptAveWeightOne", cent, (nProtonWeighted * protonPtSum - protonPtSumw2) / protonDiff, 1., rndm);
      }
    }
  }
  PROCESS_SWITCH(PidFlowPtCorr, processMCClosure, "Run truth-level MC flow closure", false);

  /**
   * @brief Fill THnSparse histograms used to derive NUA weights
   * @details Histograms use (runNumberIDX, phi, eta, Vz, pT); cfgUseNUAWithPt enables NUE weighting before filling
   *
   * @param collision
   * @param tracks
   */
  void fillCorrectionGraph(AodCollisions::iterator const& collision, aod::BCsWithTimestamps const&, AodTracks const& tracks)
  {
    // init
    int nTot = tracks.size();
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    int runNumber = bc.runNumber();
    double interactionRate = getInteractionRate(bc.timestamp(), runNumber);
    // end init

    // collision cut
    // include : 1.track.size 2.collision.sel8 3. evenSelected
    registry.fill(HIST("hEventCount/fillCorrectionGraph"), 0.5);
    if (nTot < 1) {
      return;
    }
    fGFW->Clear();
    const auto cent = getCentrality(collision);
    if (!collision.sel8()) {
      return;
    }
    registry.fill(HIST("hEventCount/fillCorrectionGraph"), 1.5);
    if (!eventSelected(collision, cent, interactionRate, MyFunctionName::funcFillCorrectionGraph)) {
      return;
    }
    // end collision cut

    if (switchsOpts.cfgUseNUAWithPt.value) {
      loadCorrections(bc.timestamp());
    }

    // loop the vector, find the place to put (phi eta Vz)
    // if the run number is new, create one
    int matchedPosition = -1;
    for (uint64_t idxPosition = 0; idxPosition < runNumbers.size(); idxPosition++) {
      if (runNumbers[idxPosition] == runNumber) {
        matchedPosition = idxPosition;
        break;
      }
    }
    if (matchedPosition == -1) {
      return;
    }
    // end find place to put run data

    // loop all the track
    for (const auto& track : tracks) {
      // track cut
      if (!trackSelectedGlobal(track, false)) {
        continue;
      }
      if (!track.hasITS()) {
        continue;
      }
      if (!track.hasTPC()) {
        continue;
      }
      if (!trackSelected4ITS(track)) {
        continue;
      }
      if (!trackSelected4TPC(track)) {
        continue;
      }
      // end track cut

      // fill the THn
      if (isWithinRefPtRange(track.pt())) {
        float weightNUE = 1.f;
        if (switchsOpts.cfgUseNUAWithPt.value) {
          setParticleNUEWeight(weightNUE, track, cent);
        }
        registry.fill(HIST("correction/hRunNumberPhiEtaVertex"), matchedPosition, track.phi(), track.eta(), collision.posZ(), track.pt(), weightNUE);
      }

      int pid = this->getPidConfigurable(track);
      if (!isWithinPOIPtRange(pid, track.pt())) {
        continue;
      }
      float weightNUEPid = 1.f;
      if (switchsOpts.cfgUseNUAWithPt.value) {
        setParticleNUEWeight(weightNUEPid, track, cent, pid);
      }
      switch (pid) {
        case MyParticleType::kPion:
          registry.fill(HIST("correction/hRunNumberPhiEtaVertexPion"), matchedPosition, track.phi(), track.eta(), collision.posZ(), track.pt(), weightNUEPid);
          break;

        case MyParticleType::kKaon:
          registry.fill(HIST("correction/hRunNumberPhiEtaVertexKaon"), matchedPosition, track.phi(), track.eta(), collision.posZ(), track.pt(), weightNUEPid);
          break;

        case MyParticleType::kProton:
          registry.fill(HIST("correction/hRunNumberPhiEtaVertexProton"), matchedPosition, track.phi(), track.eta(), collision.posZ(), track.pt(), weightNUEPid);
          break;

        default:
          break;
      }
      // end fill the THn

    } // end loop all the track
  }
  PROCESS_SWITCH(PidFlowPtCorr, fillCorrectionGraph, "", true);

  void processQA(AodCollisions::iterator const& collision, aod::BCsWithTimestamps const&, AodTracks const& tracks)
  {
    if (!switchsOpts.cfgProcessQAOutput.value) {
      LOGF(fatal, "processQAOutput no open!");
    }

    /// @note init && collision cut
    int nTot = tracks.size();
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    int runNumber = bc.runNumber();
    double interactionRate = getInteractionRate(bc.timestamp(), runNumber);
    // end init

    /// @note collision cut
    // include : 1.track.size 2.collision.sel8 3. evenSelected
    if (nTot < 1) {
      return;
    }

    const auto cent = getCentrality(collision);
    if (!collision.sel8()) {
      return;
    }

    /// @note event qa
    registry.fill(HIST("hProcessQA/centVsMult"), cent, tracks.size());
    registry.fill(HIST("hProcessQA/IR"), interactionRate);
    registry.fill(HIST("hProcessQA/Occupacy"), collision.trackOccupancyInTimeRange());
    registry.fill(HIST("hProcessQA/occupancyVsCent"), cent, collision.trackOccupancyInTimeRange());
    // end event qa

    // i dont want to fill event count again as it's filled in other 5 function
    if (!eventSelected(collision, cent, interactionRate, MyFunctionName::funcProcessQA)) {
      return;
    }
    // end init && collision cut

    /// @note track loop
    for (const auto& track : tracks) {
      /// @note select global track
      // track cut isnt applied, this is track QA function!
      if (!track.hasITS()) {
        continue;
      }
      if (!track.hasTPC()) {
        continue;
      }
      // end select global track

      /// @note fill QA graphs
      registry.fill(HIST("hProcessQA/hDCAz"), track.dcaZ(), track.pt());
      registry.fill(HIST("hProcessQA/hDCAxy"), track.dcaXY(), track.pt());

      registry.fill(HIST("hProcessQA/hChi2prTPCcls"), track.tpcChi2NCl());
      registry.fill(HIST("hProcessQA/hChi2prITScls"), track.itsChi2NCl());

      registry.fill(HIST("hProcessQA/hnTPCClu"), track.tpcNClsFound());
      registry.fill(HIST("hProcessQA/hnITSClu"), track.itsNCls());
      registry.fill(HIST("hProcessQA/hnTPCCrossedRow"), track.tpcNClsCrossedRows());
      // end fill QA grahps
    }
    // end track loop
  }
  PROCESS_SWITCH(PidFlowPtCorr, processQA, "", true);

  /**
   * @brief PID-performance QA for TPC, TOF, PID clusters and ITS response.
   * @note open switch outputQA if use it
   *
   * @param collision
   * @param tracks
   */
  void detectorPidQA(AodCollisions::iterator const& collision, aod::BCsWithTimestamps const&, AodTracksForDetectorPidQA const& tracks)
  {
    if (!switchsOpts.cfgOutputQA.value) {
      return;
    }

    // cut and correction
    int nTot = tracks.size();
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    int runNumber = bc.runNumber();
    double interactionRate = getInteractionRate(bc.timestamp(), runNumber);
    // end init

    // collision cut
    registry.fill(HIST("hEventCount/detectorPidQA"), 0.5);
    if (nTot < 1) {
      return;
    }
    fGFW->Clear();
    const auto cent = getCentrality(collision);
    if (!collision.sel8()) {
      return;
    }
    registry.fill(HIST("hEventCount/detectorPidQA"), 1.5);
    if (!eventSelected(collision, cent, interactionRate, MyFunctionName::funcDetectorPidQA)) {
      return;
    }
    // loadCorrections(bc.timestamp());
    // end cut and correction

    // start filling graphs
    for (const auto& track : tracks) {
      // Use the same track-quality selection as processData, but do not
      // require a PID hypothesis for the all-particle QA histograms.
      if (!trackSelectedForFlow(track)) {
        continue;
      }
      const float pt = track.pt();
      if (std::abs(track.sign()) == 1) {
        registry.fill(HIST("DetectorPidQA/TPC/dEdxAll"), pt, track.tpcSignal());
        if (track.hasTOF()) {
          registry.fill(HIST("DetectorPidQA/TOF/tofTimeAll"), pt, track.tofSignal() - track.tofEvTime());
        }
      }

      const int pid = getPidConfigurable(track);
      fillPidQualityQA(track, cent, pid);
      // The existing species-dependent TPC/TOF histograms require the current
      // unique PID assignment and its corresponding POI pT range.
      if (!isWithinPOIPtRange(pid, pt)) {
        continue;
      }

      switch (pid) {
        case MyParticleType::kPion:
          registry.fill(HIST("DetectorPidQA/TPC/dEdxPi"), pt, track.tpcSignal());
          if (track.hasTOF()) {
            registry.fill(HIST("DetectorPidQA/TOF/deltaPi"), pt, o2::aod::pidutils::tofExpSignalDiff<o2::track::PID::Pion>(track));
          }
          break;
        case MyParticleType::kKaon:
          registry.fill(HIST("DetectorPidQA/TPC/dEdxKa"), pt, track.tpcSignal());
          if (track.hasTOF()) {
            registry.fill(HIST("DetectorPidQA/TOF/deltaKa"), pt, o2::aod::pidutils::tofExpSignalDiff<o2::track::PID::Kaon>(track));
          }
          break;
        case MyParticleType::kProton:
          registry.fill(HIST("DetectorPidQA/TPC/dEdxPr"), pt, track.tpcSignal());
          if (track.hasTOF()) {
            registry.fill(HIST("DetectorPidQA/TOF/deltaPr"), pt, o2::aod::pidutils::tofExpSignalDiff<o2::track::PID::Proton>(track));
          }
          break;
        default:
          break;
      }
    } // end filling graphs
  }
  PROCESS_SWITCH(PidFlowPtCorr, detectorPidQA, "", true);

  /**
   * @brief Fill PID response matrix with MC truth PID versus reconstructed PID.
   * @note The reconstructed PID uses getPidConfigurable(), so it follows the same PID cuts
   *       as the flow and efficiency parts of this task. Reco bin Unidentified contains
   *       tracks that pass quality cuts but fail or ambiguously pass PID selection.
   */
  void processPidResponseMatrix(FilteredCollisionsWithMCLabel::iterator const& collision,
                                aod::BCsWithTimestamps const&,
                                FilteredTracksWithMCLabel const& tracks,
                                aod::McParticles const&,
                                aod::McCollisions const&)
  {
    if (!switchsOpts.cfgAddPidResponseMatrixHistograms.value) {
      return;
    }
    registry.fill(HIST("hEventCount/processPidResponseMatrix"), 0.5);
    if (tracks.size() <= 0) {
      return;
    }
    if (!collision.sel8()) {
      return;
    }
    registry.fill(HIST("hEventCount/processPidResponseMatrix"), 1.5);

    const auto cent = getCentrality(collision);
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    int runNumber = bc.runNumber();
    double interactionRate = getInteractionRate(bc.timestamp(), runNumber);

    if (!eventSelected(collision, cent, interactionRate, MyFunctionName::funcProcessPidResponseMatrix)) {
      return;
    }

    for (const auto& track : tracks) {
      if (!trackSelectedGlobal(track)) {
        continue;
      }
      if (!track.hasITS()) {
        continue;
      }
      if (!track.hasTPC()) {
        continue;
      }
      if (!trackSelected4ITS(track)) {
        continue;
      }
      if (!trackSelected4TPC(track)) {
        continue;
      }
      if (!track.has_mcParticle()) {
        continue;
      }

      auto mcParticle = track.mcParticle();
      if (!mcParticle.isPhysicalPrimary() || !particleSelected(mcParticle)) {
        continue;
      }

      const int truePidBin = getPidResponseTrueBin(mcParticle.pdgCode());
      const int recoPidBin = getPidResponseRecoBin(getPidConfigurable(track));

      registry.fill(HIST("pidResponseMatrix/hPidResponseMatrix"), truePidBin, recoPidBin);
      registry.fill(HIST("pidResponseMatrix/hPidResponseMatrixRecPtCent"), truePidBin, recoPidBin, track.pt(), cent);
    }
  }
  PROCESS_SWITCH(PidFlowPtCorr, processPidResponseMatrix, "Fill PID response matrix from MC truth and reconstructed PID", false);

  /**
   * @brief this function is used to fill Reco particles in order to calculate pt eff
   * @note Efficiency = Efficiency(pt, cent)
   *
   * @param collision
   * @param tracks
   */
  void processReco(FilteredCollisionsWithMCLabel::iterator const& collision,
                   aod::BCsWithTimestamps const&,
                   FilteredTracksWithMCLabel const& tracks,
                   aod::McParticles const&,
                   aod::McCollisions const&)
  {
    // init && cut
    registry.fill(HIST("hEventCount/processReco"), 0.5);
    if (tracks.size() <= 0) {
      return;
    }
    if (!collision.sel8()) {
      return;
    }
    registry.fill(HIST("hEventCount/processReco"), 1.5);

    const auto cent = getCentrality(collision);
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    int runNumber = bc.runNumber();
    double interactionRate = getInteractionRate(bc.timestamp(), runNumber);

    if (!eventSelected(collision, cent, interactionRate, MyFunctionName::funcProcessReco)) {
      return;
    }
    // end init && cut

    if (switchsOpts.cfgOutPutPtSpectra.value) {
      registry.fill(HIST("ptSpectra/hCentEventCountMcRec"), cent);
    } // cfgoutputptspectra

    // loop tracks
    for (const auto& track : tracks) {
      // track cut
      if (!trackSelectedGlobal(track)) {
        continue;
      }
      // end track cut

      // loop track's matched MC particle
      if (track.has_mcParticle()) {
        auto mcParticle = track.mcParticle();
        // fill graph
        if (particleSelected(mcParticle) && mcParticle.isPhysicalPrimary()) {
          /// @note global track, fill rec hist
          if (track.hasITS() && track.hasTPC() && trackSelected4ITS(track) && trackSelected4TPC(track)) {
            // graph for all particles
            registry.fill(HIST("correction/hPtCentMcRec"), track.pt(), cent);

            switch (getTrueSpeciesBin(mcParticle)) {
              case kPidRespPion:
                registry.fill(HIST("correction/hPtCentMcRecTruePi"), track.pt(), cent);
                break;
              case kPidRespKaon:
                registry.fill(HIST("correction/hPtCentMcRecTrueKa"), track.pt(), cent);
                break;
              case kPidRespProton:
                registry.fill(HIST("correction/hPtCentMcRecTruePr"), track.pt(), cent);
                break;
              default:
                registry.fill(HIST("correction/hPtCentMcRecTrueOther"), track.pt(), cent);
                break;
            }

            // ------------------------------
            // Unified PID logic (configurable)
            // ------------------------------
            int pid = getPidConfigurable(track);

            // Fill MC reco histograms based on unified result
            if (pid == MyParticleType::kPion) {
              registry.fill(HIST("correction/hPtCentMcRecPi"), track.pt(), cent);
            } else if (pid == MyParticleType::kKaon) {
              registry.fill(HIST("correction/hPtCentMcRecKa"), track.pt(), cent);
            } else if (pid == MyParticleType::kProton) {
              registry.fill(HIST("correction/hPtCentMcRecPr"), track.pt(), cent);
            }
            // else: do nothing (ambiguous or not identified)
            // end identify particle and fill graph
          }
          // end global track, fill rec hist
        }
        // end fill graph
      }
      // end loop track's matched MC particle
    }
    // end loop tracks
  }
  PROCESS_SWITCH(PidFlowPtCorr, processReco, "function used to do pt eff, NOTE (OutPutMc, processReco, processSim) should be open", true);

  /**
   * @brief this function is used to calculate pt eff
   *
   * @param mcCollision
   * @param collisions
   * @param mcParticles
   */
  void processSim(FilteredMcCollisions::iterator const&,
                  aod::BCsWithTimestamps const&,
                  soa::SmallGroups<soa::Join<aod::McCollisionLabels, AodCollisions>> const& collisions,
                  FilteredMcParticles const& mcParticles,
                  FilteredTracksWithMCLabel const& tracks)
  {
    registry.fill(HIST("hEventCount/processSim"), 0.5);
    // Keep the generated denominator aligned with one reconstructed collision.
    // Multiple associations would otherwise fill the same MC particles repeatedly.
    if (collisions.size() != 1) {
      return;
    }

    // cut && init
    // loop all the collisions reco matched to MC

    for (const auto& oneColl : collisions) {
      auto bc = oneColl.bc_as<aod::BCsWithTimestamps>();
      int runNumber = bc.runNumber();
      double interactionRate = getInteractionRate(bc.timestamp(), runNumber);
      double cent = getCentrality(oneColl);

      auto groupedTracks = tracks.sliceBy(perCollision, oneColl.globalIndex());
      if (groupedTracks.size() <= 0) {
        continue;
      }

      // collision cut
      if (!oneColl.sel8()) {
        continue;
      }
      registry.fill(HIST("hEventCount/processSim"), 1.5);

      if (!eventSelected(oneColl, cent, interactionRate, MyFunctionName::funcProcessSim)) {
        continue;
      }
      // end collision cut

      if (switchsOpts.cfgOutPutPtSpectra.value) {
        registry.fill(HIST("ptSpectra/hCentEventCountMcGen"), cent);
      } // cfgoutputptspectra

      // loop mc particles
      for (const auto& mcParticle : mcParticles) {
        // fill graph
        if (particleSelected(mcParticle) && mcParticle.isPhysicalPrimary()) {
          // graph for all particles
          registry.fill(HIST("correction/hPtCentMcGen"), mcParticle.pt(), cent);

          switch (getTrueSpeciesBin(mcParticle)) {
            case kPidRespPion:
              registry.fill(HIST("correction/hPtCentMcGenTruePi"), mcParticle.pt(), cent);
              break;
            case kPidRespKaon:
              registry.fill(HIST("correction/hPtCentMcGenTrueKa"), mcParticle.pt(), cent);
              break;
            case kPidRespProton:
              registry.fill(HIST("correction/hPtCentMcGenTruePr"), mcParticle.pt(), cent);
              break;
            default:
              registry.fill(HIST("correction/hPtCentMcGenTrueOther"), mcParticle.pt(), cent);
              break;
          }

          // identify particle and fill graph
          if (std::abs(mcParticle.pdgCode()) == PDG_t::kPiPlus) {
            registry.fill(HIST("correction/hPtCentMcGenPi"), mcParticle.pt(), cent);
          }
          if (std::abs(mcParticle.pdgCode()) == PDG_t::kKPlus) {
            registry.fill(HIST("correction/hPtCentMcGenKa"), mcParticle.pt(), cent);
          }
          if (std::abs(mcParticle.pdgCode()) == PDG_t::kProton) {
            registry.fill(HIST("correction/hPtCentMcGenPr"), mcParticle.pt(), cent);
          }
          // end identify particle and fill graph
        }
        // end fill graph
      }
      // end loop mc particles
    }
    // end loop all the collisions reco matched to MC
    // end cut && init
  }
  PROCESS_SWITCH(PidFlowPtCorr, processSim, "function used to do pt eff, NOTE (OutPutMc, processReco, processSim) should be open", true);

  void processReco1D(FilteredCollisionsWithMCLabel::iterator const& collision,
                     aod::BCsWithTimestamps const&,
                     FilteredTracksWithMCLabel const& tracks,
                     aod::McParticles const&,
                     aod::McCollisions const&)
  {
    if (!collision.sel8()) {
      return;
    }
    if (tracks.size() < 1) {
      return;
    }

    for (const auto& track : tracks) {

      if (!trackSelectedGlobal(track)) {
        continue;
      }
      if (!track.hasITS()) {
        continue;
      }
      if (!track.hasTPC()) {
        continue;
      }
      if (!trackSelected4ITS(track)) {
        continue;
      }
      if (!trackSelected4TPC(track)) {
        continue;
      }

      if (track.has_mcParticle()) {
        auto mcParticle = track.mcParticle();

        if (isStable(mcParticle.pdgCode())) {
          registry.fill(HIST("correction1D/hPtMCRec1D"), track.pt());

          switch (std::abs(mcParticle.pdgCode())) {
            case PDG_t::kPiPlus:
              registry.fill(HIST("correction1D/hPtMCRec1DPi"), track.pt());
              break;
            case PDG_t::kKPlus:
              registry.fill(HIST("correction1D/hPtMCRec1DKa"), track.pt());
              break;
            case PDG_t::kProton:
              registry.fill(HIST("correction1D/hPtMCRec1DPr"), track.pt());
              break;
            default:
              break;
          } // switch PDGCODE
        } // stable
      } // track.has_mcparticle
    } // end track loop
  }
  PROCESS_SWITCH(PidFlowPtCorr, processReco1D, "process reconstructed information", true);

  void processSim1D(FilteredMcCollisions::iterator const&,
                    aod::BCsWithTimestamps const&,
                    soa::SmallGroups<soa::Join<aod::McCollisionLabels, AodCollisions>> const& collisions,
                    FilteredMcParticles const& mcParticles,
                    FilteredTracksWithMCLabel const&)
  {
    if (collisions.size() > -1) {
      for (const auto& mcParticle : mcParticles) {
        if (mcParticle.isPhysicalPrimary() && isStable(mcParticle.pdgCode())) {
          registry.fill(HIST("correction1D/hPtMCGen1D"), mcParticle.pt());
          switch (std::abs(mcParticle.pdgCode())) {
            case PDG_t::kPiPlus:
              registry.fill(HIST("correction1D/hPtMCGen1DPi"), mcParticle.pt());
              break;
            case PDG_t::kKPlus:
              registry.fill(HIST("correction1D/hPtMCGen1DKa"), mcParticle.pt());
              break;
            case PDG_t::kProton:
              registry.fill(HIST("correction1D/hPtMCGen1DPr"), mcParticle.pt());
              break;
            default:
              break;
          }
        }
      }
    }
  }
  PROCESS_SWITCH(PidFlowPtCorr, processSim1D, "process pure simulation information", true);

  // main function
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<PidFlowChargedTimeProducer>(cfgc),
    adaptAnalysisTask<PidFlowPtCorr>(cfgc)};
}
