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
#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/TriggerAliases.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TProfile.h>

#include <cstdlib>

using namespace o2;
using namespace o2::framework;

struct CentralityQa {
  HistogramRegistry histos{"histos"};

  bool isRun2 = false;
  bool isMC = false;

  Configurable<int> nBins{"nBins", 1050, "number of bins"};
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {1000, 0, 1000}, "Multiplicity"};
  ConfigurableAxis axisMultiplicityPV{"axisMultiplicityPV", {1000, 0, 1000}, "Multiplicity PV"};

  struct : ConfigurableGroup {
    std::string prefix = "eventSelections"; // JSON group name
    Configurable<bool> requireSel8{"requireSel8", true, "require sel8 event selection"};
    Configurable<bool> requireTriggerTVX{"requireTriggerTVX", true, "require FT0 vertex (acceptable FT0C-FT0A time difference) at trigger level"};
    Configurable<bool> rejectITSROFBorder{"rejectITSROFBorder", true, "reject events at ITS ROF border (Run 3 only)"};
    Configurable<bool> rejectTFBorder{"rejectTFBorder", true, "reject events at TF border (Run 3 only)"};
    Configurable<bool> requireIsVertexITSTPC{"requireIsVertexITSTPC", false, "require events with at least one ITS-TPC track (Run 3 only)"};
    Configurable<bool> requireIsGoodZvtxFT0VsPV{"requireIsGoodZvtxFT0VsPV", true, "require events with PV position along z consistent (within 1 cm) between PV reconstructed using tracks and PV using FT0 A-C time difference (Run 3 only)"};
    Configurable<bool> requireIsVertexTOFmatched{"requireIsVertexTOFmatched", false, "require events with at least one of vertex contributors matched to TOF (Run 3 only)"};
    Configurable<bool> requireIsVertexTRDmatched{"requireIsVertexTRDmatched", false, "require events with at least one of vertex contributors matched to TRD (Run 3 only)"};
    Configurable<bool> rejectSameBunchPileup{"rejectSameBunchPileup", true, "reject collisions in case of pileup with another collision in the same foundBC (Run 3 only)"};
    Configurable<bool> requireNoCollInTimeRangeStd{"requireNoCollInTimeRangeStd", false, "reject collisions corrupted by the cannibalism, with other collisions within +/- 2 microseconds or mult above a certain threshold in -4 - -2 microseconds (Run 3 only)"};
    Configurable<bool> requireNoCollInTimeRangeStrict{"requireNoCollInTimeRangeStrict", false, "reject collisions corrupted by the cannibalism, with other collisions within +/- 10 microseconds (Run 3 only)"};
    Configurable<bool> requireNoCollInTimeRangeNarrow{"requireNoCollInTimeRangeNarrow", false, "reject collisions corrupted by the cannibalism, with other collisions within +/- 2 microseconds (Run 3 only)"};
    Configurable<bool> requireNoCollInROFStd{"requireNoCollInROFStd", false, "reject collisions corrupted by the cannibalism, with other collisions within the same ITS ROF with mult. above a certain threshold (Run 3 only)"};
    Configurable<bool> requireNoCollInROFStrict{"requireNoCollInROFStrict", false, "reject collisions corrupted by the cannibalism, with other collisions within the same ITS ROF (Run 3 only)"};
    Configurable<bool> requireINEL0{"requireINEL0", true, "require INEL>0 event selection"};
    Configurable<bool> requireINEL1{"requireINEL1", false, "require INEL>1 event selection"};

    Configurable<float> maxZVtxPosition{"maxZVtxPosition", 10., "max Z vtx position"};

    Configurable<bool> useFT0CbasedOccupancy{"useFT0CbasedOccupancy", false, "Use sum of FT0-C amplitudes for estimating occupancy? (if not, use track-based definition)"};
    // fast check on occupancy
    Configurable<float> minOccupancy{"minOccupancy", -1, "minimum occupancy from neighbouring collisions"};
    Configurable<float> maxOccupancy{"maxOccupancy", -1, "maximum occupancy from neighbouring collisions"};

    Configurable<bool> requireIsBBT0A{"requireIsBBT0A", false, "Require beam-beam collisions based on timing information in FT0A"};
    Configurable<bool> requireIsBBT0C{"requireIsBBT0C", false, "Require beam-beam collisions based on timing information in FT0C"};

    Configurable<bool> rejectMismatchedBCs{"rejectMismatchedBCs", false, "Reject collision with BC different from MC BC"};
    Configurable<bool> rejectMismatchedFoundBCs{"rejectMismatchedFoundBCs", false, "Reject collision with found BC different from MC BC"};

    // Run 2 specific event selections
    Configurable<bool> requireSel7{"requireSel7", true, "require sel7 event selection (Run 2 only: event selection decision based on V0A & V0C)"};
    Configurable<bool> requireINT7{"requireINT7", true, "require INT7 trigger selection (Run 2 only)"};
    Configurable<bool> rejectIncompleteDAQ{"rejectIncompleteDAQ", true, "reject events with incomplete DAQ (Run 2 only)"};
    Configurable<bool> requireConsistentSPDAndTrackVtx{"requireConsistentSPDAndTrackVtx", true, "reject events with inconsistent in SPD and Track vertices (Run 2 only)"};
    Configurable<bool> rejectPileupFromSPD{"rejectPileupFromSPD", true, "reject events with pileup according to SPD vertexer (Run 2 only)"};
    Configurable<bool> rejectV0PFPileup{"rejectV0PFPileup", false, "reject events tagged as OOB pileup according to V0 past-future info (Run 2 only)"};
    Configurable<bool> rejectPileupInMultBins{"rejectPileupInMultBins", true, "reject events tagged as pileup according to multiplicity-differential pileup checks (Run 2 only)"};
    Configurable<bool> rejectPileupMV{"rejectPileupMV", true, "reject events tagged as pileup according to according to multi-vertexer (Run 2 only)"};
    Configurable<bool> rejectTPCPileup{"rejectTPCPileup", false, "reject events tagged as pileup according to pileup in TPC (Run 2 only)"};
    Configurable<bool> requireNoV0MOnVsOffPileup{"requireNoV0MOnVsOffPileup", false, "reject events tagged as OOB pileup according to online-vs-offline VOM correlation (Run 2 only)"};
    Configurable<bool> requireNoSPDOnVsOffPileup{"requireNoSPDOnVsOffPileup", false, "reject events tagged as pileup according to online-vs-offline SPD correlation (Run 2 only)"};
    Configurable<bool> requireNoSPDClsVsTklBG{"requireNoSPDClsVsTklBG", true, "reject events tagged as beam-gas and pileup according to cluster-vs-tracklet correlation (Run 2 only)"};
  } eventSelections;

  void init(o2::framework::InitContext& /*initContext*/)
  {
    if (doprocessRun2PP ||
        doprocessRun2PPb ||
        doprocessRun2PbPb) {
      isRun2 = true;
    } else {
      isRun2 = false;
    }

    if (doprocessMonteCarloRun3_FV0A ||
        doprocessMonteCarloRun3_FT0M ||
        doprocessMonteCarloRun3_FT0A ||
        doprocessMonteCarloRun3_FT0C ||
        doprocessMonteCarloRun3_FT0CVar1 ||
        doprocessMonteCarloRun3_FT0CVar2 ||
        doprocessMonteCarloRun3_MFT ||
        doprocessMonteCarloRun3_NGlobal ||
        doprocessMonteCarloRun3_NTPV) {
      isMC = true;
    } else {
      isMC = false;
    }

    if (isRun2) {
      histos.add("hCentRun2V0M", "V0M centrality (%)", kTH1D, {{nBins, 0, 105.}});
      histos.add("hCentRun2V0A", "V0A centrality (%)", kTH1D, {{nBins, 0, 105.}});
      histos.add("hCentRun2SPDTks", "SPD tracklet centrality (%)", kTH1D, {{nBins, 0, 105.}});
      histos.add("hCentRun2SPDCls", "SPD cluster centrality (%)", kTH1D, {{nBins, 0, 105.}});
      histos.add("hCentRun2CL0", "CL0 centrality (%)", kTH1D, {{nBins, 0, 105.}});
      histos.add("hCentRun2CL1", "CL1 centrality (%)", kTH1D, {{nBins, 0, 105.}});
    } else {
      histos.add("hCentFV0A", "FV0A centrality (%)", kTH1D, {{nBins, 0, 105.}});
      histos.add("hCentFT0M", "FT0M centrality (%)", kTH1D, {{nBins, 0, 105.}});
      histos.add("hCentFT0A", "FT0A centrality (%)", kTH1D, {{nBins, 0, 105.}});
      histos.add("hCentFT0C", "FT0C centrality (%)", kTH1D, {{nBins, 0, 105.}});
      histos.add("hCentFT0CVar1", "FT0CVar1 centrality (%)", kTH1D, {{nBins, 0, 105.}});
      histos.add("hCentFT0CVar2", "FT0CVar2 centrality (%)", kTH1D, {{nBins, 0, 105.}});
      histos.add("hCentFDDM", "FDDM centrality (%)", kTH1D, {{nBins, 0, 105.}});
      histos.add("hCentNTPV", "NTPV centrality (%)", kTH1D, {{nBins, 0, 105.}});
      histos.add("hCentNGlobal", "NGlobal centrality (%)", kTH1D, {{nBins, 0, 105.}});
      histos.add("hCentMFT", "MFT centrality (%)", kTH1D, {{nBins, 0, 105.}});

      // profiles of midrapidity multiplicity density
      histos.add("hCentProfileFV0A", "FV0A centrality (%)", kTProfile, {{nBins, 0, 105.}});
      histos.add("hCentProfileFT0M", "FT0M centrality (%)", kTProfile, {{nBins, 0, 105.}});
      histos.add("hCentProfileFT0A", "FT0A centrality (%)", kTProfile, {{nBins, 0, 105.}});
      histos.add("hCentProfileFT0C", "FT0C centrality (%)", kTProfile, {{nBins, 0, 105.}});
      histos.add("hCentProfileFT0CVar1", "FT0CVar1 centrality (%)", kTProfile, {{nBins, 0, 105.}});
      histos.add("hCentProfileFT0CVar2", "FT0CVar2 centrality (%)", kTProfile, {{nBins, 0, 105.}});
      histos.add("hCentProfileFDDM", "FDDM centrality (%)", kTProfile, {{nBins, 0, 105.}});
      histos.add("hCentProfileNTPV", "NTPV centrality (%)", kTProfile, {{nBins, 0, 105.}});
      histos.add("hCentProfileNGlobal", "NGlobal centrality (%)", kTProfile, {{nBins, 0, 105.}});
      histos.add("hCentProfileMFT", "MFT centrality (%)", kTProfile, {{nBins, 0, 105.}});

      histos.add("hMultEta05VsCentFV0A", "FV0A centrality (%); Multiplicity PV contributors (|#it{#eta}| < 0.5)", kTH2D, {{nBins, 0, 105.}, axisMultiplicityPV});
      histos.add("hMultEta05VsCentFT0M", "FT0M centrality (%); Multiplicity PV contributors (|#it{#eta}| < 0.5)", kTH2D, {{nBins, 0, 105.}, axisMultiplicityPV});
      histos.add("hMultEta05VsCentFT0A", "FT0A centrality (%); Multiplicity PV contributors (|#it{#eta}| < 0.5)", kTH2D, {{nBins, 0, 105.}, axisMultiplicityPV});
      histos.add("hMultEta05VsCentFT0C", "FT0C centrality (%); Multiplicity PV contributors (|#it{#eta}| < 0.5)", kTH2D, {{nBins, 0, 105.}, axisMultiplicityPV});
      histos.add("hMultEta05VsCentFT0CVar1", "FT0CVar1 centrality (%); Multiplicity PV contributors (|#it{#eta}| < 0.5)", kTH2D, {{nBins, 0, 105.}, axisMultiplicityPV});
      histos.add("hMultEta05VsCentFT0CVar2", "FT0CVar2 centrality (%); Multiplicity PV contributors (|#it{#eta}| < 0.5)", kTH2D, {{nBins, 0, 105.}, axisMultiplicityPV});
      histos.add("hMultEta05VsCentFDDM", "FDDM centrality (%); Multiplicity PV contributors (|#it{#eta}| < 0.5)", kTH2D, {{nBins, 0, 105.}, axisMultiplicityPV});
      histos.add("hMultEta05VsCentNTPV", "NTPV centrality (%); Multiplicity PV contributors (|#it{#eta}| < 0.5)", kTH2D, {{nBins, 0, 105.}, axisMultiplicityPV});
      histos.add("hMultEta05VsCentNGlobal", "NGlobal centrality (%); Multiplicity PV contributors (|#it{#eta}| < 0.5)", kTH2D, {{nBins, 0, 105.}, axisMultiplicityPV});
      histos.add("hMultEta05VsCentMFT", "MFT centrality (%); Multiplicity PV contributors (|#it{#eta}| < 0.5)", kTH2D, {{nBins, 0, 105.}, axisMultiplicityPV});

      if (isMC) {
        histos.add("hMultEta05VsGenMultFV0A", "Multiplicity FV0A; Multiplicity PV contributors (|#it{#eta}| < 0.5)", kTH2D, {axisMultiplicity, axisMultiplicityPV});
        histos.add("hMultEta05VsGenMultFT0M", "Multiplicity FT0M; Multiplicity PV contributors (|#it{#eta}| < 0.5)", kTH2D, {axisMultiplicity, axisMultiplicityPV});
        histos.add("hMultEta05VsGenMultFT0A", "Multiplicity FT0A; Multiplicity PV contributors (|#it{#eta}| < 0.5)", kTH2D, {axisMultiplicity, axisMultiplicityPV});
        histos.add("hMultEta05VsGenMultFT0C", "Multiplicity FT0C; Multiplicity PV contributors (|#it{#eta}| < 0.5)", kTH2D, {axisMultiplicity, axisMultiplicityPV});
        histos.add("hMultEta05VsGenMultFT0CVar1", "Multiplicity FT0CVar1; Multiplicity PV contributors (|#it{#eta}| < 0.5)", kTH2D, {axisMultiplicity, axisMultiplicityPV});
        histos.add("hMultEta05VsGenMultFT0CVar2", "Multiplicity FT0CVar2; Multiplicity PV contributors (|#it{#eta}| < 0.5)", kTH2D, {axisMultiplicity, axisMultiplicityPV});
        histos.add("hMultEta05VsGenMultFDDM", "Multiplicity FDDM; Multiplicity PV contributors (|#it{#eta}| < 0.5)", kTH2D, {axisMultiplicity, axisMultiplicityPV});
        histos.add("hMultEta05VsGenMultNTPV", "Multiplicity NTPV; Multiplicity PV contributors (|#it{#eta}| < 0.5)", kTH2D, {axisMultiplicity, axisMultiplicityPV});
        histos.add("hMultEta05VsGenMultNGlobal", "Multiplicity NGlobal; Multiplicity PV contributors (|#it{#eta}| < 0.5)", kTH2D, {axisMultiplicity, axisMultiplicityPV});
        histos.add("hMultEta05VsGenMultMFT", "Multiplicity MFT; Multiplicity PV contributors (|#it{#eta}| < 0.5)", kTH2D, {axisMultiplicity, axisMultiplicityPV});
      }
    }

    histos.print();
  }

  template <typename TCollision>
  bool isEventAccepted(TCollision const& collision)
  // check whether the collision passes our collision selections
  {
    if constexpr (
      requires { collision.centFV0A(); } ||
      requires { collision.centFT0M(); } ||
      requires { collision.centFT0A(); } ||
      requires { collision.centFT0C(); } ||
      requires { collision.centFT0CVariant1(); } ||
      requires { collision.centFT0CVariant2(); } ||
      requires { collision.centFDDM(); } ||
      requires { collision.centNTPV(); } ||
      requires { collision.centNGlobal(); } ||
      requires { collision.centMFT(); }) { // check if we are in Run 3
      if (eventSelections.requireSel8 && !collision.sel8()) {
        return false;
      }

      if (eventSelections.requireTriggerTVX && !collision.selection_bit(aod::evsel::kIsTriggerTVX)) {
        return false;
      }

      if (eventSelections.rejectITSROFBorder && !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
        return false;
      }

      if (eventSelections.rejectTFBorder && !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
        return false;
      }

      if (std::abs(collision.posZ()) > eventSelections.maxZVtxPosition) {
        return false;
      }

      if (eventSelections.requireIsBBT0A && !collision.selection_bit(aod::evsel::kIsBBT0A)) {
        return false;
      }

      if (eventSelections.requireIsBBT0C && !collision.selection_bit(aod::evsel::kIsBBT0C)) {
        return false;
      }

      if (eventSelections.requireIsVertexITSTPC && !collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
        return false;
      }

      if (eventSelections.requireIsGoodZvtxFT0VsPV && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
        return false;
      }

      if (eventSelections.requireIsVertexTOFmatched && !collision.selection_bit(o2::aod::evsel::kIsVertexTOFmatched)) {
        return false;
      }

      if (eventSelections.requireIsVertexTRDmatched && !collision.selection_bit(o2::aod::evsel::kIsVertexTRDmatched)) {
        return false;
      }

      if (eventSelections.rejectSameBunchPileup && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
        return false;
      }

      if (eventSelections.requireNoCollInTimeRangeStd && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
        return false;
      }

      if (eventSelections.requireNoCollInTimeRangeStrict && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStrict)) {
        return false;
      }

      if (eventSelections.requireNoCollInTimeRangeNarrow && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeNarrow)) {
        return false;
      }

      if (eventSelections.requireNoCollInROFStd && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) {
        return false;
      }

      if (eventSelections.requireNoCollInROFStrict && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStrict)) {
        return false;
      }

      if (eventSelections.requireINEL0 && collision.multNTracksPVeta1() < 1) {
        return false;
      }

      if (eventSelections.requireINEL1 && collision.multNTracksPVeta1() < 2) {
        return false;
      }

      float collisionOccupancy = eventSelections.useFT0CbasedOccupancy ? collision.ft0cOccupancyInTimeRange() : collision.trackOccupancyInTimeRange();
      if (eventSelections.minOccupancy >= 0 && collisionOccupancy < eventSelections.minOccupancy) {
        return false;
      }

      if (eventSelections.maxOccupancy >= 0 && collisionOccupancy > eventSelections.maxOccupancy) {
        return false;
      }

      if constexpr (requires { collision.has_mcCollision(); }) { // check if we are in MC
        if (!collision.has_mcCollision()) {
          return false;
        }

        const auto& mcCollision = collision.template mcCollision_as<soa::Join<aod::McCollisions, aod::MultMCExtras>>();
        const auto& recoBC = collision.template bc_as<soa::Join<aod::BCs, aod::Run3MatchedToBCSparse>>();
        const auto& foundBC = collision.template foundBC_as<soa::Join<aod::BCs, aod::Run3MatchedToBCSparse>>();
        const auto& mcBC = mcCollision.template bc_as<soa::Join<aod::BCs, aod::Run3MatchedToBCSparse>>();

        // Check that the BC in data and MC is the same
        if (eventSelections.rejectMismatchedBCs && recoBC.globalBC() != mcBC.globalBC()) {
          return false;
        }
        if (eventSelections.rejectMismatchedFoundBCs && foundBC.globalBC() != mcBC.globalBC()) {
          return false;
        }
      }
    } else { // we are in Run 2
      if (eventSelections.requireSel8 && !collision.sel8()) {
        return false;
      }

      if (eventSelections.requireSel7 && !collision.sel7()) {
        return false;
      }

      if (eventSelections.requireINT7 && !collision.alias_bit(kINT7)) {
        return false;
      }

      if (eventSelections.requireTriggerTVX && !collision.selection_bit(o2::aod::evsel::kIsTriggerTVX)) {
        return false;
      }

      if (eventSelections.rejectIncompleteDAQ && !collision.selection_bit(o2::aod::evsel::kNoIncompleteDAQ)) {
        return false;
      }

      if (std::abs(collision.posZ()) > eventSelections.maxZVtxPosition) {
        return false;
      }

      if (eventSelections.requireConsistentSPDAndTrackVtx && !collision.selection_bit(o2::aod::evsel::kNoInconsistentVtx)) {
        return false;
      }

      if (eventSelections.rejectPileupFromSPD && !collision.selection_bit(o2::aod::evsel::kNoPileupFromSPD)) {
        return false;
      }

      if (eventSelections.rejectV0PFPileup && !collision.selection_bit(o2::aod::evsel::kNoV0PFPileup)) {
        return false;
      }

      if (eventSelections.rejectPileupInMultBins && !collision.selection_bit(o2::aod::evsel::kNoPileupInMultBins)) {
        return false;
      }

      if (eventSelections.rejectPileupMV && !collision.selection_bit(o2::aod::evsel::kNoPileupMV)) {
        return false;
      }

      if (eventSelections.rejectTPCPileup && !collision.selection_bit(o2::aod::evsel::kNoPileupTPC)) {
        return false;
      }

      if (eventSelections.requireNoV0MOnVsOffPileup && !collision.selection_bit(o2::aod::evsel::kNoV0MOnVsOfPileup)) {
        return false;
      }

      if (eventSelections.requireNoSPDOnVsOffPileup && !collision.selection_bit(o2::aod::evsel::kNoSPDOnVsOfPileup)) {
        return false;
      }

      if (eventSelections.requireNoSPDClsVsTklBG && !collision.selection_bit(o2::aod::evsel::kNoSPDClsVsTklBG)) {
        return false;
      }

      if (eventSelections.requireINEL0 && collision.multNTracksPVeta1() < 1) {
        return false;
      }

      if (eventSelections.requireINEL1 && collision.multNTracksPVeta1() < 2) {
        return false;
      }
    }

    return true;
  }

  void processRun2PP(soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0Ms, aod::CentRun2SPDTrks, aod::CentRun2SPDClss, aod::Mults>::iterator const& col)
  {
    if (!isEventAccepted(col)) {
      return;
    }
    LOGF(debug, "centV0M=%.0f", col.centRun2V0M());
    LOGF(debug, "centSPDTracklets=%.0f", col.centRun2SPDTracklets());
    LOGF(debug, "centSPDClusters=%.0f", col.centRun2SPDClusters());
    // fill centrality histos
    histos.fill(HIST("hCentRun2V0M"), col.centRun2V0M());
    histos.fill(HIST("hCentRun2SPDTks"), col.centRun2SPDTracklets());
    histos.fill(HIST("hCentRun2SPDCls"), col.centRun2SPDClusters());
  }
  PROCESS_SWITCH(CentralityQa, processRun2PP, "Process with Run2 SPD clusters centrality/multiplicity estimation", false);

  void processRun2PbPb(soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0Ms, aod::CentRun2SPDTrks, aod::CentRun2CL0s, aod::CentRun2CL1s, aod::Mults>::iterator const& col)
  {
    if (!isEventAccepted(col)) {
      return;
    }
    LOGF(debug, "centV0M=%.0f", col.centRun2V0M());
    LOGF(debug, "centSPDTracklets=%.0f", col.centRun2SPDTracklets());
    LOGF(debug, "centCL0=%.0f", col.centRun2CL0());
    LOGF(debug, "centCL1=%.0f", col.centRun2CL1());
    // fill centrality histos
    histos.fill(HIST("hCentRun2V0M"), col.centRun2V0M());
    histos.fill(HIST("hCentRun2SPDTks"), col.centRun2SPDTracklets());
    histos.fill(HIST("hCentRun2CL0"), col.centRun2CL0());
    histos.fill(HIST("hCentRun2CL1"), col.centRun2CL1());
  }
  PROCESS_SWITCH(CentralityQa, processRun2PbPb, "Process with Run2 CL0 and CL1 multiplicities centrality/multiplicity  estimation", false);

  void processRun2PPb(soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0As, aod::Mults>::iterator const& col)
  {
    if (!isEventAccepted(col)) {
      return;
    }
    LOGF(debug, "centV0A=%.0f", col.centRun2V0A());
    // fill centrality histos
    histos.fill(HIST("hCentRun2V0A"), col.centRun2V0A());
  }
  PROCESS_SWITCH(CentralityQa, processRun2PPb, "Process with Run2 V0A multiplicitY centrality/multiplicity  estimation", false);

  void processRun3_FV0A(soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFV0As>::iterator const& col)
  {
    if (!isEventAccepted(col)) {
      return;
    }
    LOGF(debug, "centFV0A=%.0f", col.centFV0A());
    histos.fill(HIST("hCentFV0A"), col.centFV0A());
    histos.fill(HIST("hCentProfileFV0A"), col.centFV0A(), col.multNTracksPVetaHalf());
    histos.fill(HIST("hMultEta05VsCentFV0A"), col.centFV0A(), col.multNTracksPVetaHalf());
  }
  PROCESS_SWITCH(CentralityQa, processRun3_FV0A, "Process with Run 3 FV0A estimator", false);

  void processRun3_FT0M(soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Ms>::iterator const& col)
  {
    if (!isEventAccepted(col)) {
      return;
    }
    LOGF(debug, "centFT0M=%.0f", col.centFT0M());
    histos.fill(HIST("hCentFT0M"), col.centFT0M());
    histos.fill(HIST("hCentProfileFT0M"), col.centFT0M(), col.multNTracksPVetaHalf());
    histos.fill(HIST("hMultEta05VsCentFT0M"), col.centFT0M(), col.multNTracksPVetaHalf());
  }
  PROCESS_SWITCH(CentralityQa, processRun3_FT0M, "Process with Run 3 FT0M estimator", false);

  void processRun3_FT0A(soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0As>::iterator const& col)
  {
    if (!isEventAccepted(col)) {
      return;
    }
    histos.fill(HIST("hCentFT0A"), col.centFT0A());
    histos.fill(HIST("hCentProfileFT0A"), col.centFT0A(), col.multNTracksPVetaHalf());
    histos.fill(HIST("hMultEta05VsCentFT0A"), col.centFT0A(), col.multNTracksPVetaHalf());
  }
  PROCESS_SWITCH(CentralityQa, processRun3_FT0A, "Process with Run 3 FT0A estimator", false);

  void processRun3_FT0C(soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Cs>::iterator const& col)
  {
    if (!isEventAccepted(col)) {
      return;
    }
    histos.fill(HIST("hCentFT0C"), col.centFT0C());
    histos.fill(HIST("hCentProfileFT0C"), col.centFT0C(), col.multNTracksPVetaHalf());
    histos.fill(HIST("hMultEta05VsCentFT0C"), col.centFT0C(), col.multNTracksPVetaHalf());
  }
  PROCESS_SWITCH(CentralityQa, processRun3_FT0C, "Process with Run 3 FT0C estimator", false);

  void processRun3_FT0CVar1(soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0CVariant1s>::iterator const& col)
  {
    if (!isEventAccepted(col)) {
      return;
    }
    histos.fill(HIST("hCentFT0CVar1"), col.centFT0CVariant1());
    histos.fill(HIST("hCentProfileFT0CVar1"), col.centFT0CVariant1(), col.multNTracksPVetaHalf());
    histos.fill(HIST("hMultEta05VsCentFT0CVar1"), col.centFT0CVariant1(), col.multNTracksPVetaHalf());
  }
  PROCESS_SWITCH(CentralityQa, processRun3_FT0CVar1, "Process with Run 3 FT0CVar1 estimator", false);

  void processRun3_FT0CVar2(soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0CVariant2s>::iterator const& col)
  {
    if (!isEventAccepted(col)) {
      return;
    }
    histos.fill(HIST("hCentFT0CVar2"), col.centFT0CVariant2());
    histos.fill(HIST("hCentProfileFT0CVar2"), col.centFT0CVariant2(), col.multNTracksPVetaHalf());
    histos.fill(HIST("hMultEta05VsCentFT0CVar2"), col.centFT0CVariant2(), col.multNTracksPVetaHalf());
  }
  PROCESS_SWITCH(CentralityQa, processRun3_FT0CVar2, "Process with Run 3 FT0CVar2 estimator", false);

  void processRun3_FDDM(soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFDDMs>::iterator const& col)
  {
    if (!isEventAccepted(col)) {
      return;
    }
    histos.fill(HIST("hCentFDDM"), col.centFDDM());
    histos.fill(HIST("hCentProfileFDDM"), col.centFDDM(), col.multNTracksPVetaHalf());
    histos.fill(HIST("hMultEta05VsCentFDDM"), col.centFDDM(), col.multNTracksPVetaHalf());
  }
  PROCESS_SWITCH(CentralityQa, processRun3_FDDM, "Process with Run 3 FDDM estimator", false);

  void processRun3_NTPV(soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentNTPVs>::iterator const& col)
  {
    if (!isEventAccepted(col)) {
      return;
    }
    histos.fill(HIST("hCentNTPV"), col.centNTPV());
    histos.fill(HIST("hCentProfileNTPV"), col.centNTPV(), col.multNTracksPVetaHalf());
    histos.fill(HIST("hMultEta05VsCentNTPV"), col.centNTPV(), col.multNTracksPVetaHalf());
  }
  PROCESS_SWITCH(CentralityQa, processRun3_NTPV, "Process with Run 3 NTPV estimator", false);

  void processRun3_NGlobal(soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentNGlobals>::iterator const& col)
  {
    if (!isEventAccepted(col)) {
      return;
    }
    histos.fill(HIST("hCentNGlobal"), col.centNGlobal());
    histos.fill(HIST("hCentProfileNGlobal"), col.centNGlobal(), col.multNTracksPVetaHalf());
    histos.fill(HIST("hMultEta05VsCentNGlobal"), col.centNGlobal(), col.multNTracksPVetaHalf());
  }
  PROCESS_SWITCH(CentralityQa, processRun3_NGlobal, "Process with Run 3 NGlobal estimator", false);

  void processRun3_MFT(soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentMFTs>::iterator const& col)
  {
    if (!isEventAccepted(col)) {
      return;
    }
    histos.fill(HIST("hCentMFT"), col.centMFT());
    histos.fill(HIST("hCentProfileMFT"), col.centMFT(), col.multNTracksPVetaHalf());
    histos.fill(HIST("hMultEta05VsCentMFT"), col.centMFT(), col.multNTracksPVetaHalf());
  }
  PROCESS_SWITCH(CentralityQa, processRun3_MFT, "Process with Run 3 MFT estimator", false);

  void processMonteCarloRun3_FV0A(soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::Mults, aod::CentFV0As>::iterator const& col,
                                  soa::Join<aod::McCollisions, aod::MultMCExtras> const& /*mcCollisions*/,
                                  soa::Join<aod::BCs, aod::Run3MatchedToBCSparse> const& /*bcs*/)
  {
    if (!isEventAccepted(col)) {
      return;
    }
    const auto& mcCol = col.mcCollision_as<soa::Join<aod::McCollisions, aod::MultMCExtras>>();

    LOGF(debug, "centFV0A=%.0f", col.centFV0A());
    histos.fill(HIST("hCentFV0A"), col.centFV0A());
    histos.fill(HIST("hCentProfileFV0A"), col.centFV0A(), col.multNTracksPVetaHalf());
    histos.fill(HIST("hMultEta05VsCentFV0A"), col.centFV0A(), col.multNTracksPVetaHalf());
    histos.fill(HIST("hMultEta05VsGenMultFV0A"), mcCol.multMCFV0A(), col.multNTracksPVetaHalf());
  }
  PROCESS_SWITCH(CentralityQa, processMonteCarloRun3_FV0A, "Process with Run 3 FV0A estimator", false);

  void processMonteCarloRun3_FT0M(soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::Mults, aod::CentFT0Ms>::iterator const& col,
                                  soa::Join<aod::McCollisions, aod::MultMCExtras> const& /*mcCollisions*/,
                                  soa::Join<aod::BCs, aod::Run3MatchedToBCSparse> const& /*bcs*/)
  {
    if (!isEventAccepted(col)) {
      return;
    }
    const auto& mcCol = col.mcCollision_as<soa::Join<aod::McCollisions, aod::MultMCExtras>>();

    LOGF(debug, "centFT0M=%.0f", col.centFT0M());
    histos.fill(HIST("hCentFT0M"), col.centFT0M());
    histos.fill(HIST("hCentProfileFT0M"), col.centFT0M(), col.multNTracksPVetaHalf());
    histos.fill(HIST("hMultEta05VsCentFT0M"), col.centFT0M(), col.multNTracksPVetaHalf());
    histos.fill(HIST("hMultEta05VsGenMultFT0M"), mcCol.multMCFT0A() + mcCol.multMCFT0C(), col.multNTracksPVetaHalf());
  }
  PROCESS_SWITCH(CentralityQa, processMonteCarloRun3_FT0M, "Process with Run 3 FT0M estimator", false);

  void processMonteCarloRun3_FT0A(soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::Mults, aod::CentFT0As>::iterator const& col,
                                  soa::Join<aod::McCollisions, aod::MultMCExtras> const& /*mcCollisions*/,
                                  soa::Join<aod::BCs, aod::Run3MatchedToBCSparse> const& /*bcs*/)
  {
    if (!isEventAccepted(col)) {
      return;
    }
    const auto& mcCol = col.mcCollision_as<soa::Join<aod::McCollisions, aod::MultMCExtras>>();

    histos.fill(HIST("hCentFT0A"), col.centFT0A());
    histos.fill(HIST("hCentProfileFT0A"), col.centFT0A(), col.multNTracksPVetaHalf());
    histos.fill(HIST("hMultEta05VsCentFT0A"), col.centFT0A(), col.multNTracksPVetaHalf());
    histos.fill(HIST("hMultEta05VsGenMultFT0A"), mcCol.multMCFT0A(), col.multNTracksPVetaHalf());
  }
  PROCESS_SWITCH(CentralityQa, processMonteCarloRun3_FT0A, "Process with Run 3 FT0A estimator", false);

  void processMonteCarloRun3_FT0C(soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::Mults, aod::CentFT0Cs>::iterator const& col,
                                  soa::Join<aod::McCollisions, aod::MultMCExtras> const& /*mcCollisions*/,
                                  soa::Join<aod::BCs, aod::Run3MatchedToBCSparse> const& /*bcs*/)
  {
    if (!isEventAccepted(col)) {
      return;
    }
    const auto& mcCol = col.mcCollision_as<soa::Join<aod::McCollisions, aod::MultMCExtras>>();

    histos.fill(HIST("hCentFT0C"), col.centFT0C());
    histos.fill(HIST("hCentProfileFT0C"), col.centFT0C(), col.multNTracksPVetaHalf());
    histos.fill(HIST("hMultEta05VsCentFT0C"), col.centFT0C(), col.multNTracksPVetaHalf());
    histos.fill(HIST("hMultEta05VsGenMultFT0C"), mcCol.multMCFT0C(), col.multNTracksPVetaHalf());
  }
  PROCESS_SWITCH(CentralityQa, processMonteCarloRun3_FT0C, "Process with Run 3 FT0C estimator", false);

  void processMonteCarloRun3_FT0CVar1(soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::Mults, aod::CentFT0CVariant1s>::iterator const& col,
                                      soa::Join<aod::McCollisions, aod::MultMCExtras> const& /*mcCollisions*/,
                                      soa::Join<aod::BCs, aod::Run3MatchedToBCSparse> const& /*bcs*/)
  {
    if (!isEventAccepted(col)) {
      return;
    }
    const auto& mcCol = col.mcCollision_as<soa::Join<aod::McCollisions, aod::MultMCExtras>>();

    histos.fill(HIST("hCentFT0CVar1"), col.centFT0CVariant1());
    histos.fill(HIST("hCentProfileFT0CVar1"), col.centFT0CVariant1(), col.multNTracksPVetaHalf());
    histos.fill(HIST("hMultEta05VsCentFT0CVar1"), col.centFT0CVariant1(), col.multNTracksPVetaHalf());
    histos.fill(HIST("hMultEta05VsGenMultFT0CVar1"), mcCol.multMCFT0C(), col.multNTracksPVetaHalf());
  }
  PROCESS_SWITCH(CentralityQa, processMonteCarloRun3_FT0CVar1, "Process with Run 3 FT0CVar1 estimator", false);

  void processMonteCarloRun3_FT0CVar2(soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::Mults, aod::CentFT0CVariant2s>::iterator const& col,
                                      soa::Join<aod::McCollisions, aod::MultMCExtras> const& /*mcCollisions*/,
                                      soa::Join<aod::BCs, aod::Run3MatchedToBCSparse> const& /*bcs*/)
  {
    if (!isEventAccepted(col)) {
      return;
    }
    const auto& mcCol = col.mcCollision_as<soa::Join<aod::McCollisions, aod::MultMCExtras>>();

    histos.fill(HIST("hCentFT0CVar2"), col.centFT0CVariant2());
    histos.fill(HIST("hCentProfileFT0CVar2"), col.centFT0CVariant2(), col.multNTracksPVetaHalf());
    histos.fill(HIST("hMultEta05VsCentFT0CVar2"), col.centFT0CVariant2(), col.multNTracksPVetaHalf());
    histos.fill(HIST("hMultEta05VsGenMultFT0CVar2"), mcCol.multMCFT0C(), col.multNTracksPVetaHalf());
  }
  PROCESS_SWITCH(CentralityQa, processMonteCarloRun3_FT0CVar2, "Process with Run 3 FT0CVar2 estimator", false);

  void processMonteCarloRun3_FDDM(soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::Mults, aod::CentFDDMs>::iterator const& col,
                                  soa::Join<aod::McCollisions, aod::MultMCExtras> const& /*mcCollisions*/,
                                  soa::Join<aod::BCs, aod::Run3MatchedToBCSparse> const& /*bcs*/)
  {
    if (!isEventAccepted(col)) {
      return;
    }
    const auto& mcCol = col.mcCollision_as<soa::Join<aod::McCollisions, aod::MultMCExtras>>();

    histos.fill(HIST("hCentFDDM"), col.centFDDM());
    histos.fill(HIST("hCentProfileFDDM"), col.centFDDM(), col.multNTracksPVetaHalf());
    histos.fill(HIST("hMultEta05VsCentFDDM"), col.centFDDM(), col.multNTracksPVetaHalf());
    histos.fill(HIST("hMultEta05VsGenMultFDDM"), mcCol.multMCFDDA() + mcCol.multMCFDDC(), col.multNTracksPVetaHalf());
  }
  PROCESS_SWITCH(CentralityQa, processMonteCarloRun3_FDDM, "Process with Run 3 FDDM estimator", false);

  void processMonteCarloRun3_NTPV(soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::Mults, aod::CentNTPVs>::iterator const& col,
                                  soa::Join<aod::McCollisions, aod::MultMCExtras> const& /*mcCollisions*/,
                                  soa::Join<aod::BCs, aod::Run3MatchedToBCSparse> const& /*bcs*/)
  {
    if (!isEventAccepted(col)) {
      return;
    }
    const auto& mcCol = col.mcCollision_as<soa::Join<aod::McCollisions, aod::MultMCExtras>>();

    histos.fill(HIST("hCentFDDM"), col.centNTPV());
    histos.fill(HIST("hCentProfileNTPV"), col.centNTPV(), col.multNTracksPVetaHalf());
    histos.fill(HIST("hMultEta05VsCentNTPV"), col.centNTPV(), col.multNTracksPVetaHalf());
    histos.fill(HIST("hMultEta05VsGenMultNTPV"), mcCol.multMCNParticlesEta08(), col.multNTracksPVetaHalf());
  }
  PROCESS_SWITCH(CentralityQa, processMonteCarloRun3_NTPV, "Process with Run 3 NTPV estimator", false);

  void processMonteCarloRun3_NGlobal(soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::Mults, aod::CentNGlobals>::iterator const& col,
                                     soa::Join<aod::McCollisions, aod::MultMCExtras> const& /*mcCollisions*/,
                                     soa::Join<aod::BCs, aod::Run3MatchedToBCSparse> const& /*bcs*/)
  {
    if (!isEventAccepted(col)) {
      return;
    }
    const auto& mcCol = col.mcCollision_as<soa::Join<aod::McCollisions, aod::MultMCExtras>>();

    histos.fill(HIST("hCentNGlobal"), col.centNGlobal());
    histos.fill(HIST("hCentProfileNGlobal"), col.centNGlobal(), col.multNTracksPVetaHalf());
    histos.fill(HIST("hMultEta05VsCentNGlobal"), col.centNGlobal(), col.multNTracksPVetaHalf());
    histos.fill(HIST("hMultEta05VsGenMultNGlobal"), mcCol.multMCNParticlesEta08(), col.multNTracksPVetaHalf());
  }
  PROCESS_SWITCH(CentralityQa, processMonteCarloRun3_NGlobal, "Process with Run 3 NGlobal estimator", false);

  void processMonteCarloRun3_MFT(soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::Mults, aod::CentMFTs>::iterator const& col,
                                 soa::Join<aod::McCollisions, aod::MultMCExtras> const& /*mcCollisions*/,
                                 soa::Join<aod::BCs, aod::Run3MatchedToBCSparse> const& /*bcs*/)
  {
    if (!isEventAccepted(col)) {
      return;
    }
    // const auto& mcCol = col.mcCollision_as<soa::Join<aod::McCollisions, aod::MultMCExtras>>(); // FIXME: uncomment when MC MFT mult is added in aod::MultMCExtras

    histos.fill(HIST("hCentMFT"), col.centMFT());
    histos.fill(HIST("hCentProfileMFT"), col.centMFT(), col.multNTracksPVetaHalf());
    histos.fill(HIST("hMultEta05VsCentMFT"), col.centMFT(), col.multNTracksPVetaHalf());
    // histos.fill(HIST("hMultEta05VsGenMultMFT"), mcCol.multMCMFT(), col.multNTracksPVetaHalf()); // FIXME: uncomment when MC MFT mult is added in aod::MultMCExtras
  }
  PROCESS_SWITCH(CentralityQa, processMonteCarloRun3_MFT, "Process with Run 3 MFT estimator", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<CentralityQa>(cfgc)};
}
