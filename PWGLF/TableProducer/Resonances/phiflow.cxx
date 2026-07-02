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

/// \file phiflow.cxx
/// \brief Table producer for Phi FLow
///
/// \author prottay.das@cern.ch

#include "PWGLF/DataModel/LFPhiFlowTables.h"
#include "PWGLF/DataModel/SPCalibrationTables.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/RCTSelectionFlags.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/PID.h>

#include <Math/Vector4D.h> // IWYU pragma: keep (do not replace with Math/Vector4Dfwd.h)

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::rctsel;

struct phiflow {

  Produces<aod::KaonEvents> kaonEvent;
  Produces<aod::KaonTracks> kaonTrack;

  struct : ConfigurableGroup {
    Configurable<bool> requireRCTFlagChecker{"requireRCTFlagChecker", true, "Check event quality in run condition table"};
    Configurable<std::string> cfgEvtRCTFlagCheckerLabel{"cfgEvtRCTFlagCheckerLabel", "CBT_hadronPID", "Evt sel: RCT flag checker label"};
    Configurable<bool> cfgEvtRCTFlagCheckerZDCCheck{"cfgEvtRCTFlagCheckerZDCCheck", false, "Evt sel: RCT flag checker ZDC check"};
    Configurable<bool> cfgEvtRCTFlagCheckerLimitAcceptAsBad{"cfgEvtRCTFlagCheckerLimitAcceptAsBad", true, "Evt sel: RCT flag checker treat Limited Acceptance As Bad"};
  } rctCut;

  Configurable<bool> useNoCollInTimeRangeStandard{"useNoCollInTimeRangeStandard", false, "Apply kNoCollInTimeRangeStandard selection bit"};
  Configurable<bool> useGoodITSLayersAll{"useGoodITSLayersAll", true, "Apply kIsGoodITSLayersAll selection bit"};
  Configurable<int> cfgCutOccupancy{"cfgCutOccupancy", 2000, "Occupancy cut"};

  // events
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCutCentralityMax{"cfgCutCentralityMax", 80.0f, "Accepted maximum Centrality"};
  Configurable<float> cfgCutCentralityMin{"cfgCutCentralityMin", 0.0f, "Accepted minimum Centrality"};

  // Configs for track filtering
  Configurable<float> cfgCutCharge{"cfgCutCharge", 0.0, "Charge cut on daughter"};
  Configurable<float> cfgCutPt{"cfgCutPt", 0.2, "Pt cut on daughter track"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8, "Eta cut on daughter track"};
  Configurable<float> cfgCutDCAxy{"cfgCutDCAxy", 0.2, "DCAxy cut on daughter track"};
  Configurable<float> cfgCutDCAz{"cfgCutDCAz", 0.2, "DCAz cut on daughter track"};
  Configurable<float> nsigmaCutTPC{"nsigmaCutTPC", 3.0, "Maximum nsigma cut TPC for filtered kaon track"};

  // Configs for kaon
  struct : ConfigurableGroup {
    Configurable<bool> itsPIDSelection{"itsPIDSelection", true, "PID ITS"};
    Configurable<float> lowITSPIDNsigma{"lowITSPIDNsigma", -3.0, "lower cut on PID nsigma for ITS"};
    Configurable<float> highITSPIDNsigma{"highITSPIDNsigma", 3.0, "higher cut on PID nsigma for ITS"};
    Configurable<int> itsclusterKaMeson{"itsclusterKaMeson", 5, "Minimum number of ITS cluster for kaon meson track"};
    Configurable<int> tpcCrossedRowsKaMeson{"tpcCrossedRowsKaMeson", 80, "Minimum number of TPC Crossed Rows for kaon meson track"};
    Configurable<float> cutDCAxyKaMeson{"cutDCAxyKaMeson", 0.1, "Maximum DCAxy for kaon meson track"};
    Configurable<float> cutDCAzKaMeson{"cutDCAzKaMeson", 0.1, "Maximum DCAz for kaon meson track"};
    Configurable<float> cutEtaKaMeson{"cutEtaKaMeson", 0.8, "Maximum eta for kaon meson track"};
    Configurable<float> cutPTKaMeson{"cutPTKaMeson", 0.8, "Minimum pt for kaon meson track"};
    Configurable<bool> usePID{"usePID", true, "Flag for using PID selection for kaon meson track"};
    Configurable<float> nsigmaCutTPCKaMeson{"nsigmaCutTPCKaMeson", 3.0, "Maximum nsigma cut TPC for kaon meson track"};
    Configurable<float> nsigmaCutTOFKaMeson{"nsigmaCutTOFKaMeson", 3.0, "Maximum nsigma cut TOF for kaon meson track"};
    Configurable<float> cutTOFBetaKaMeson{"cutTOFBetaKaMeson", 0.5f, "Maximum beta cut for kaon meson track"};
    Configurable<float> pSwitchPID{"pSwitchPID", 0.5f, "pT switch for pT-dependent kaon PID"};
    Configurable<float> pSwitchAsymPID{"pSwitchAsymPID", 0.7f, "Momentum switch for asymmetric kaon PID"};
    Configurable<bool> applyTOFAsymPID{"applyTOFAsymPID", true, "Apply TOF requirement in PID7"};
  } grpKaon;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  RCTFlagsChecker rctChecker;
  void init(o2::framework::InitContext&)
  {
    rctChecker.init(rctCut.cfgEvtRCTFlagCheckerLabel, rctCut.cfgEvtRCTFlagCheckerZDCCheck, rctCut.cfgEvtRCTFlagCheckerLimitAcceptAsBad);

    histos.add("hCent", "hCent", kTH1F, {{8, 0, 80.0}});
    histos.add("hEvtSelInfo", "hEvtSelInfo", kTH1F, {{5, 0, 5.0}});
    histos.add("hTrkSelInfo", "hTrkSelInfo", kTH1F, {{10, 0, 10.0}});
  }

  template <typename T>
  bool selectionTrack(const T& candidate)
  {
    return candidate.isGlobalTrack() &&
           candidate.isPVContributor() &&
           candidate.itsNCls() >= grpKaon.itsclusterKaMeson &&
           candidate.tpcNClsCrossedRows() > grpKaon.tpcCrossedRowsKaMeson &&
           std::abs(candidate.dcaXY()) <= grpKaon.cutDCAxyKaMeson &&
           std::abs(candidate.dcaZ()) <= grpKaon.cutDCAzKaMeson &&
           std::abs(candidate.eta()) <= grpKaon.cutEtaKaMeson &&
           candidate.pt() >= grpKaon.cutPTKaMeson;
  }

  enum KaonPidBits : uint8_t {
    kPID1 = 1u << 0, // TOF-availability-dependent rectangular PID
    kPID2 = 1u << 1, // TOF-availability-dependent circular PID
    kPID3 = 1u << 2, // pT-dependent circular PID
    kPID4 = 1u << 3, // pT-dependent rectangular PID, phi-v2-like
    kPID5 = 1u << 4, // TOF-only PID
    kPID6 = 1u << 5, // Momentum-dependent asymmetric TOF PID
    kPID7 = 1u << 6  // Optional-TOF asymmetric PID
  };

  template <typename T>
  bool selectionPID(const T& candidate)
  {
    const float nTPC = candidate.tpcNSigmaKa();

    // No TOF hit: TPC-only PID at all pT.
    if (!candidate.hasTOF()) {
      return std::abs(nTPC) < grpKaon.nsigmaCutTPCKaMeson;
    }

    // TOF hit: require beta, TPC PID, and TOF PID.
    if (candidate.beta() <= grpKaon.cutTOFBetaKaMeson) {
      return false;
    }

    return std::abs(nTPC) < grpKaon.nsigmaCutTPCKaMeson &&
           std::abs(candidate.tofNSigmaKa()) < grpKaon.nsigmaCutTOFKaMeson;
  }

  template <typename T>
  bool selectionPID2(const T& candidate)
  {
    const float nTPC = candidate.tpcNSigmaKa();

    // No TOF hit: TPC-only PID at all pT.
    if (!candidate.hasTOF()) {
      return std::abs(nTPC) < grpKaon.nsigmaCutTPCKaMeson;
    }

    if (candidate.beta() <= grpKaon.cutTOFBetaKaMeson) {
      return false;
    }

    const float nTOF = candidate.tofNSigmaKa();
    const float nCombined = std::sqrt(nTPC * nTPC + nTOF * nTOF);

    return nCombined < grpKaon.nsigmaCutTOFKaMeson;
  }

  template <typename T>
  bool selectionPID3(const T& candidate)
  {
    const float pt = candidate.pt();
    const float nTPC = candidate.tpcNSigmaKa();

    // Low pT: TPC-only. TOF is ignored even if present.
    if (pt < grpKaon.pSwitchPID) {
      return std::abs(nTPC) < grpKaon.nsigmaCutTPCKaMeson;
    }

    // High pT: TOF is mandatory.
    if (!candidate.hasTOF()) {
      return false;
    }

    if (candidate.beta() <= grpKaon.cutTOFBetaKaMeson) {
      return false;
    }

    const float nTOF = candidate.tofNSigmaKa();
    const float nCombined = std::sqrt(nTPC * nTPC + nTOF * nTOF);

    return nCombined < grpKaon.nsigmaCutTOFKaMeson;
  }

  template <typename T>
  bool selectionPID4(const T& candidate)
  {
    const float pt = candidate.pt();
    const float nTPC = candidate.tpcNSigmaKa();

    // Low pT: TPC-only. TOF is ignored even if present.
    if (pt < grpKaon.pSwitchPID) {
      return std::abs(nTPC) < grpKaon.nsigmaCutTPCKaMeson;
    }

    // High pT: TOF is mandatory.
    if (!candidate.hasTOF()) {
      return false;
    }

    if (candidate.beta() <= grpKaon.cutTOFBetaKaMeson) {
      return false;
    }

    const float nTOF = candidate.tofNSigmaKa();

    return std::abs(nTPC) < grpKaon.nsigmaCutTPCKaMeson &&
           std::abs(nTOF) < grpKaon.nsigmaCutTOFKaMeson;
  }

  template <typename T>
  bool selectionPID5(const T& candidate)
  {
    // This matches the TOF-only selection from the other task.
    if (!candidate.hasTOF()) {
      return false;
    }

    if (candidate.beta() <= grpKaon.cutTOFBetaKaMeson) {
      return false;
    }

    return std::abs(candidate.tofNSigmaKa()) <
           grpKaon.nsigmaCutTOFKaMeson;
  }

  template <typename T>
  bool selectionPID6(const T& candidate)
  {
    const float p = candidate.p();
    const float nTPC = candidate.tpcNSigmaKa();

    // This follows selectionPIDpTdependent2 from the other task.
    // Note: it uses total momentum p, not transverse momentum pT.
    if (p < grpKaon.pSwitchAsymPID) {
      return std::abs(nTPC) < grpKaon.nsigmaCutTPCKaMeson;
    }

    if (!candidate.hasTOF()) {
      return false;
    }

    if (std::abs(nTPC) >= grpKaon.nsigmaCutTPCKaMeson) {
      return false;
    }

    if (candidate.beta() <= grpKaon.cutTOFBetaKaMeson) {
      return false;
    }

    const float nTOF = candidate.tofNSigmaKa();

    if (p < 1.6f) {
      return nTOF > -5.0f && nTOF < 10.0f;
    }

    if (p < 2.0f) {
      return nTOF > -3.0f && nTOF < 10.0f;
    }

    if (p < 2.5f) {
      return nTOF > -3.0f && nTOF < 6.0f;
    }

    if (p < 4.0f) {
      return nTOF > -2.5f && nTOF < 4.0f;
    }

    if (p < 5.0f) {
      return nTOF > -4.0f && nTOF < 3.0f;
    }

    if (p < 6.0f) {
      return nTOF > -4.0f && nTOF < 2.5f;
    }

    return nTOF > -3.0f && nTOF < 3.0f;
  }

  template <typename T>
  bool selectionPID7(const T& candidate)
  {
    const float p = candidate.p();
    const float nTPC = candidate.tpcNSigmaKa();

    // This follows selectionPID22 from the other task.
    // If TOF is disabled, use only the TPC kaon selection.
    if (!grpKaon.applyTOFAsymPID) {
      return std::abs(nTPC) < grpKaon.nsigmaCutTPCKaMeson;
    }

    // No TOF hit: allow TPC-only PID at all momenta.
    if (!candidate.hasTOF()) {
      return std::abs(nTPC) < grpKaon.nsigmaCutTPCKaMeson;
    }

    // For tracks with TOF, asymmetric TOF selection starts above p = 0.5.
    if (p <= 0.5f) {
      return false;
    }

    if (std::abs(nTPC) >= grpKaon.nsigmaCutTPCKaMeson) {
      return false;
    }

    if (candidate.beta() <= grpKaon.cutTOFBetaKaMeson) {
      return false;
    }

    const float nTOF = candidate.tofNSigmaKa();

    if (p < 1.6f) {
      return nTOF > -5.0f && nTOF < 10.0f;
    }

    if (p < 2.0f) {
      return nTOF > -3.0f && nTOF < 10.0f;
    }

    if (p < 2.5f) {
      return nTOF > -3.0f && nTOF < 6.0f;
    }

    if (p < 4.0f) {
      return nTOF > -2.5f && nTOF < 4.0f;
    }

    if (p < 5.0f) {
      return nTOF > -4.0f && nTOF < 3.0f;
    }

    if (p < 6.0f) {
      return nTOF > -4.0f && nTOF < 2.5f;
    }

    return nTOF > -3.0f && nTOF < 3.0f;
  }

  template <typename T>
  uint8_t kaonPidMask(const T& track)
  {
    uint8_t mask = 0;

    if (selectionPID(track)) {
      mask |= kPID1;
    }

    if (selectionPID2(track)) {
      mask |= kPID2;
    }

    if (selectionPID3(track)) {
      mask |= kPID3;
    }

    if (selectionPID4(track)) {
      mask |= kPID4;
    }

    if (selectionPID5(track)) {
      mask |= kPID5;
    }

    if (selectionPID6(track)) {
      mask |= kPID6;
    }

    if (selectionPID7(track)) {
      mask |= kPID7;
    }

    return mask;
  }

  struct StoredKaon {
    float px;
    float py;
    float pz;
    int8_t charge;
    int64_t trackId;
    uint8_t pidMask;
  };

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter centralityFilter = (aod::cent::centFT0C < cfgCutCentralityMax && aod::cent::centFT0C > cfgCutCentralityMin);
  Filter acceptanceFilter = (nabs(aod::track::eta) < cfgCutEta && (aod::track::pt) > cfgCutPt);
  Filter DCAcutFilter = (nabs(aod::track::dcaXY) < cfgCutDCAxy) && (nabs(aod::track::dcaZ) < cfgCutDCAz);
  Filter PIDcutFilter = nabs(aod::pidtpc::tpcNSigmaKa) < nsigmaCutTPC;

  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::FT0Mults, aod::TPCMults, aod::SPCalibrationTables>>;
  using AllTrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTOFbeta>>;
  SliceCache cache;
  Partition<AllTrackCandidates> posTracks = aod::track::signed1Pt > cfgCutCharge;
  Partition<AllTrackCandidates> negTracks = aod::track::signed1Pt < cfgCutCharge;

  void processData(EventCandidates::iterator const& collision,
                   AllTrackCandidates const& /*tracks*/)
  {
    o2::aod::ITSResponse itsResponse;

    const auto centrality = collision.centFT0C();
    const auto vz = collision.posZ();
    const int occupancy = collision.trackOccupancyInTimeRange();

    histos.fill(HIST("hEvtSelInfo"), 0.5);

    if (!collision.triggereventsp()) {
      return;
    }

    histos.fill(HIST("hEvtSelInfo"), 1.5);

    if (!((!rctCut.requireRCTFlagChecker || rctChecker(collision)) &&
          collision.selection_bit(aod::evsel::kNoSameBunchPileup) &&
          collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV) &&
          (!useNoCollInTimeRangeStandard ||
           collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) &&
          collision.sel8() &&
          (!useGoodITSLayersAll ||
           collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) &&
          occupancy < cfgCutOccupancy)) {
      return;
    }

    histos.fill(HIST("hEvtSelInfo"), 2.5);

    histos.fill(HIST("hCent"), centrality);

    std::vector<StoredKaon> selectedKaons;

    // Same slicing pattern that worked in your old producer.
    auto posThisColl = posTracks->sliceByCached(
      aod::track::collisionId,
      collision.globalIndex(),
      cache);

    auto negThisColl = negTracks->sliceByCached(
      aod::track::collisionId,
      collision.globalIndex(),
      cache);

    // ---------------- Positive kaons ----------------
    for (const auto& track1 : posThisColl) {
      histos.fill(HIST("hTrkSelInfo"), 0.5);

      if (!selectionTrack(track1)) {
        continue;
      }

      histos.fill(HIST("hTrkSelInfo"), 1.5);

      const float nSigmaITS1 =
        itsResponse.nSigmaITS<o2::track::PID::Kaon>(track1);

      if (grpKaon.itsPIDSelection && track1.p() < 1.0 &&
          (nSigmaITS1 <= grpKaon.lowITSPIDNsigma ||
           nSigmaITS1 >= grpKaon.highITSPIDNsigma)) {
        continue;
      }

      histos.fill(HIST("hTrkSelInfo"), 2.5);

      const uint8_t mask1 = kaonPidMask(track1);

      if (grpKaon.usePID && mask1 == 0) {
        continue;
      }

      histos.fill(HIST("hTrkSelInfo"), 3.5);

      selectedKaons.push_back({track1.px(),
                               track1.py(),
                               track1.pz(),
                               +1,
                               static_cast<int64_t>(track1.globalIndex()),
                               mask1});
    }

    // ---------------- Negative kaons ----------------
    for (const auto& track2 : negThisColl) {
      histos.fill(HIST("hTrkSelInfo"), 4.5);

      if (!selectionTrack(track2)) {
        continue;
      }

      histos.fill(HIST("hTrkSelInfo"), 5.5);

      const float nSigmaITS2 =
        itsResponse.nSigmaITS<o2::track::PID::Kaon>(track2);

      if (grpKaon.itsPIDSelection && track2.p() < 1.0 &&
          (nSigmaITS2 <= grpKaon.lowITSPIDNsigma ||
           nSigmaITS2 >= grpKaon.highITSPIDNsigma)) {
        continue;
      }

      histos.fill(HIST("hTrkSelInfo"), 6.5);

      const uint8_t mask2 = kaonPidMask(track2);

      if (grpKaon.usePID && mask2 == 0) {
        continue;
      }

      histos.fill(HIST("hTrkSelInfo"), 7.5);

      selectedKaons.push_back({track2.px(),
                               track2.py(),
                               track2.pz(),
                               -1,
                               static_cast<int64_t>(track2.globalIndex()),
                               mask2});
    }

    // No selected K+ or K- in this collision:
    // not writing a reduced event row.
    histos.fill(HIST("hEvtSelInfo"), 3.5);
    /*
    if (selectedKaons.empty()) {
      return;
    }
    */
    histos.fill(HIST("hEvtSelInfo"), 4.5);

    kaonEvent(centrality,
              vz,
              collision.qxZDCA(),
              collision.qxZDCC(),
              collision.qyZDCA(),
              collision.qyZDCC());

    const int64_t indexEvent = kaonEvent.lastIndex();

    // One table row per selected kaon.
    for (const auto& kaon : selectedKaons) {
      kaonTrack(indexEvent,
                kaon.px,
                kaon.py,
                kaon.pz,
                kaon.charge,
                kaon.trackId,
                kaon.pidMask);
    }
  }
  PROCESS_SWITCH(phiflow, processData, "Process data", true);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<phiflow>(cfgc)};
}
