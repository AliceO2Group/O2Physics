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
/// \brief This task pre-filters tracks to do h-resonance (phi, K*0)
///        correlations with an analysis task.
///
/// \file hResonanceCorrelationFilter.cxx
/// \author Hirak Kumar Koley (hirak.koley@cern.ch)

#include "PWGLF/DataModel/LFHResonanceCorrelationTables.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/Expressions.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <Math/Vector4D.h> // IWYU pragma: keep (do not replace with Math/Vector4Dfwd.h)
#include <Math/Vector4Dfwd.h>
#include <TH1.h>

#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::soa;
using namespace o2::constants::math;
using namespace o2::framework;
using namespace o2::framework::expressions;

using LorentzVectorPtEtaPhiMass = ROOT::Math::PtEtaPhiMVector;

enum PIDCutType {
  SquareType = 1,
  CircularType,
};

#define BIT_SET(var, nbit) ((var) |= (1 << (nbit)))
#define BIT_CHECK(var, nbit) ((var) & (1 << (nbit)))

struct HResonanceCorrelationFilter {

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Single, shared event-selection cut list used by every process function
  // (processTriggers, processAssocPions/Kaons/Hadrons, processPhis,
  // processKstars, and their MC counterparts) via isSelectedEvents(), so the
  // trigger table and every associated-track/resonance table are always
  // built over the exact same set of selected collisions.
  struct : ConfigurableGroup {
    std::string prefix = "configEvents";
    Configurable<float> cfgEvtZvtx{"cfgEvtZvtx", 10.0f, "Evt sel: Max. z-Vertex (cm)"};
    Configurable<bool> cfgEvtTriggerTVXSel{"cfgEvtTriggerTVXSel", true, "Evt sel: triggerTVX selection (MB)"};
    Configurable<bool> cfgEvtNoTFBorderCut{"cfgEvtNoTFBorderCut", true, "Evt sel: apply TF border cut"};
    Configurable<bool> cfgEvtNoITSROFrameBorderCut{"cfgEvtNoITSROFrameBorderCut", true, "Evt sel: apply NoITSRO border cut"};
    Configurable<bool> cfgEvtSel8{"cfgEvtSel8", true, "Evt Sel 8 check for offline selection"};
    Configurable<bool> cfgEvtNoSameBunchPileupCut{"cfgEvtNoSameBunchPileupCut", true, "Evt sel: reject collisions associated with the same found-by-T0 bunch crossing"};
    Configurable<bool> cfgEvtGoodZvtxFT0vsPVCut{"cfgEvtGoodZvtxFT0vsPVCut", true, "Evt sel: require small difference between z-vertex from PV and from FT0"};
    // Centrality source switch for hCentralitySelected: O-O collisions have a
    // well-defined FT0C centrality estimator; p-O (asymmetric system) does
    // not, so FV0A is used there instead. Only affects which estimator is
    // read for that one QA histogram -- doesn't change any selection cut.
    Configurable<bool> cfgUseFV0ACentrality{"cfgUseFV0ACentrality", false, "Centrality estimator for hCentralitySelected: false = O-O (FT0C), true = p-O (FV0A)"};
  } configEvents;

  // Merged from the formerly separate configTracks/generalSelections/
  // trackSelections groups -- all three were "cuts applied to a track or
  // track pair", just split by who last touched them. No member names
  // collided, so every Configurable below keeps its original key; only the
  // "<group>." prefix at each call site changed, to configTracks.
  // General, species-agnostic track-quality cuts: applied identically to
  // every track that reaches trackCut() -- trigger, associated hadron, AND
  // Phi/K*0 daughter candidates alike (see trackCut() below, and its call
  // sites at isValidTrigger(), isValidAssocTrack<Species>(), and both
  // processPhis()/processKstars()). cfgCutEta/cMinPtcut additionally gate
  // the framework-level `acceptanceFilter` Filter, applied to every track
  // in every process function before any of the above ever sees it.
  struct : ConfigurableGroup {
    std::string prefix = "configTracks";
    // Pre-selection Track cuts
    Configurable<float> cMinPtcut{"cMinPtcut", 0.15f, "Minimal pT for tracks"};
    Configurable<float> cMinTPCNClsFound{"cMinTPCNClsFound", 120, "minimum TPCNClsFound value for good track"};
    Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "Eta range for tracks"};
    Configurable<int> cfgMinCrossedRows{"cfgMinCrossedRows", 70, "min crossed rows for good track"};
    Configurable<float> cfgMaxTPCChi2NCl{"cfgMaxTPCChi2NCl", 4.0f, "max TPC chi2/clusters for good track"};
    Configurable<float> cfgMaxITSChi2NCl{"cfgMaxITSChi2NCl", 36.0f, "max ITS chi2/clusters for good track"};
    Configurable<float> cfgMinTPCCrossedRowsOverFindableCls{"cfgMinTPCCrossedRowsOverFindableCls", 0.8f, "min ratio of TPC crossed rows over findable clusters for good track"};
    Configurable<int> cfgMinITSNCls{"cfgMinITSNCls", 5, "min number of ITS hits (clusters) for good track"};

    // DCA Selections
    // DCAz to PV
    Configurable<float> cMaxDCAzToPVcut{"cMaxDCAzToPVcut", 0.1f, "Track DCAz cut to PV Maximum"};

    // Track selections
    Configurable<bool> cfgPrimaryTrack{"cfgPrimaryTrack", true, "Primary track selection"};                    // kGoldenChi2 | kDCAxy | kDCAz
    Configurable<bool> cfgGlobalWoDCATrack{"cfgGlobalWoDCATrack", true, "Global track selection without DCA"}; // kQualityTracks (kTrackType | kTPCNCls | kTPCCrossedRows | kTPCCrossedRowsOverNCls | kTPCChi2NDF | kTPCRefit | kITSNCls | kITSChi2NDF | kITSRefit | kITSHits) | kInAcceptanceTracks (kPtRange | kEtaRange)
    Configurable<bool> cfgGlobalTrack{"cfgGlobalTrack", false, "Global track selection"};                      // kGoldenChi2 | kDCAxy | kDCAz
    Configurable<bool> cfgPVContributor{"cfgPVContributor", false, "PV contributor track selection"};          // PV Contriuibutor
    Configurable<bool> cfgHasTOF{"cfgHasTOF", false, "Require TOF"};
    Configurable<bool> cTPCNClsFound{"cTPCNClsFound", false, "Switch to turn on/off TPCNClsFound cut"};

    // Track quality shared by trigger AND assoc selection (isValidTrigger()
    // and isValidAssocTrack<Species>() both read this same configurable --
    // it is NOT split into trigger-/assoc-only copies below).
    Configurable<int> minTPCNCrossedRows{"minTPCNCrossedRows", 70, "Minimum TPC crossed rows"};

    // primary particle DCAxy selections (formerly trackSelections)
    // formula: |DCAxy| <  0.004f + (0.013f / pt)
    Configurable<float> dcaXYconstant{"dcaXYconstant", 0.004, "[0] in |DCAxy| < [0]+[1]/pT"};
    Configurable<float> dcaXYpTdep{"dcaXYpTdep", 0.013, "[1] in |DCAxy| < [0]+[1]/pT"};
  } configTracks;

  // Phi/K*0 resonance-candidate cuts: applied to the RECONSTRUCTED candidate
  // (the K+K- or K-pi 4-vector), not to the daughter tracks themselves --
  // the daughters go through the shared configTracks/trackCut() cuts above
  // like every other track. Currently holds just the one configurable; kept
  // as its own group rather than folded into configTracks so it reads as
  // "resonance-level", not "track-level".
  struct : ConfigurableGroup {
    std::string prefix = "configResoDauTracks";
    Configurable<float> cfgCutRapidity{"cfgCutRapidity", 0.5f, "rapidity range for the reconstructed Phi/K*0 candidate"};
  } configResoDauTracks;

  // Trigger-hadron-only phase-space and quality cuts.
  struct : ConfigurableGroup {
    std::string prefix = "configTriggerTracks";
    Configurable<float> triggerEtaMin{"triggerEtaMin", -0.8, "triggeretamin"};
    Configurable<float> triggerEtaMax{"triggerEtaMax", 0.8, "triggeretamax"};
    Configurable<float> triggerPtCutMin{"triggerPtCutMin", 2, "triggerptmin"};
    Configurable<float> triggerPtCutMax{"triggerPtCutMax", 20, "triggerptmax"};
    Configurable<bool> triggerRequireITS{"triggerRequireITS", true, "require ITS signal in trigger tracks"};
    Configurable<int> triggerMaxTPCSharedClusters{"triggerMaxTPCSharedClusters", 200, "maximum number of shared TPC clusters (inclusive)"};
    Configurable<bool> triggerRequireL0{"triggerRequireL0", false, "require ITS L0 cluster for trigger"};
  } configTriggerTracks;

  // Associated-hadron-only phase-space and quality cuts.
  struct : ConfigurableGroup {
    std::string prefix = "configAssocTracks";
    Configurable<float> assocEtaMin{"assocEtaMin", -0.8, "triggeretamin"};
    Configurable<float> assocEtaMax{"assocEtaMax", 0.8, "triggeretamax"};
    Configurable<float> assocPtCutMin{"assocPtCutMin", 0.2, "assocptmin"};
    Configurable<float> assocPtCutMax{"assocPtCutMax", 10, "assocptmax"};
    Configurable<bool> assocRequireITS{"assocRequireITS", true, "require ITS signal in assoc tracks"};
  } configAssocTracks;

  struct : ConfigurableGroup {
    std::string prefix = "configPID";
    /// PID Selections
    Configurable<bool> cByPassTOF{"cByPassTOF", false, "By pass TOF PID selection"};                       // By pass TOF PID selection
    Configurable<int> cPIDcutType{"cPIDcutType", 2, "cPIDcutType = 1 for square cut, 2 for circular cut"}; // By pass TOF PID selection
    Configurable<bool> ispTdepPID{"ispTdepPID", false, "enable pT dependent PID"};
    // Were loose top-level Configurables; folded in here as unused (no call
    // sites reference them) and logically the same "Kaon PID" kind as the
    // kaonTPCPIDcuts/kaonTOFPIDcuts group below.
    Configurable<float> cfgTPCNsigmaKaon{"cfgTPCNsigmaKaon", 3.0f, "TPC Kaon PID"};
    Configurable<float> cfgTOFNsigmaKaon{"cfgTOFNsigmaKaon", 3.0f, "TOF Kaon PID"};

    // Kaon
    Configurable<std::vector<float>> kaonTPCPIDpTintv{"kaonTPCPIDpTintv", {0.5f}, "pT intervals for Kaon TPC PID cuts"};
    Configurable<std::vector<float>> kaonTPCPIDcuts{"kaonTPCPIDcuts", {2}, "nSigma list for Kaon TPC PID cuts"};
    Configurable<std::vector<float>> kaonTOFPIDpTintv{"kaonTOFPIDpTintv", {999.0f}, "pT intervals for Kaon TOF PID cuts"};
    Configurable<std::vector<float>> kaonTOFPIDcuts{"kaonTOFPIDcuts", {2}, "nSigma list for Kaon TOF PID cuts"};
    Configurable<std::vector<float>> kaonTPCTOFCombinedpTintv{"kaonTPCTOFCombinedpTintv", {999.0f}, "pT intervals for Kaon TPC-TOF PID cuts"};
    Configurable<std::vector<float>> kaonTPCTOFCombinedPIDcuts{"kaonTPCTOFCombinedPIDcuts", {2}, "nSigma list for Kaon TPC-TOF PID cuts"};

    // Pion (K*0 daughter)
    Configurable<std::vector<float>> pionTPCPIDpTintv{"pionTPCPIDpTintv", {0.5f}, "pT intervals for Pion TPC PID cuts"};
    Configurable<std::vector<float>> pionTPCPIDcuts{"pionTPCPIDcuts", {2}, "nSigma list for Pion TPC PID cuts"};
    Configurable<std::vector<float>> pionTOFPIDpTintv{"pionTOFPIDpTintv", {999.0f}, "pT intervals for Pion TOF PID cuts"};
    Configurable<std::vector<float>> pionTOFPIDcuts{"pionTOFPIDcuts", {2}, "nSigma list for Pion TOF PID cuts"};
    Configurable<std::vector<float>> pionTPCTOFCombinedpTintv{"pionTPCTOFCombinedpTintv", {999.0f}, "pT intervals for Pion TPC-TOF PID cuts"};
    Configurable<std::vector<float>> pionTPCTOFCombinedPIDcuts{"pionTPCTOFCombinedPIDcuts", {2}, "nSigma list for Pion TPC-TOF PID cuts"};

    // Associated kaon identification (formerly trackSelections; mirrors the
    // pion selection above, with the accept/reject roles swapped: accept
    // kaon-consistent tracks, reject pion-/proton-consistent ones). Selected
    // via the Species template tag on isValidAssocTrack<Species>() below --
    // production path (which table rows get made) is still separate per
    // process function, only the PID gate logic itself is shared.
    Configurable<float> assocKaonNSigmaTPCFOF{"assocKaonNSigmaTPCFOF", 3, "minimal n sigma in TOF and TPC for Kaon ID"};

    // Associated pion identification (formerly trackSelections)
    Configurable<float> assocPionNSigmaTPCFOF{"assocPionNSigmaTPCFOF", 3, "minimal n sigma in TOF and TPC for Pion ID"};
    Configurable<float> rejectSigma{"rejectSigma", 1, "n sigma for rejecting pion candidates"};
  } configPID;

  // must include windows for background and peak
  //
  // This filter-level cut is an acceptance pre-cut only -- it decides what
  // ever reaches AssocPhis/AssocKstars, not what counts as signal. It must
  // stay wider than whatever signal+background region hResonanceCorrelation.cxx
  // wants downstream (its massWindowConfigurationsPhi/Kstar go out to
  // maxBgNSigma=6 by default), or that background region gets silently
  // truncated before it ever reaches the analysis task. peakMass/sigma here
  // mirror the analysis task's own massWindowConfigurationsPhi/Kstar
  // (sigma = PDG Gamma/2.355 placeholder -- refit from your own peak and
  // update both files together, they must agree on what "sigma" means).
  struct : ConfigurableGroup {
    std::string prefix = "massWindowConfigurations";
    Configurable<float> maxMassNSigma{"maxMassNSigma", 12.0f, "max mass region to be considered for further analysis, in units of sigma"};
  } massWindowConfigurations;

  struct : ConfigurableGroup {
    std::string prefix = "massWindowConfigurationsPhi";
    Configurable<float> peakMass{"peakMass", 1.019455f, "Phi(1020) PDG mass (GeV)"};
    Configurable<float> sigma{"sigma", 0.0018f, "effective width (GeV) for the maxMassNSigma pre-cut -- keep in sync with hResonanceCorrelation.cxx's massWindowConfigurationsPhi.sigma"};
  } massWindowConfigurationsPhi;

  struct : ConfigurableGroup {
    std::string prefix = "massWindowConfigurationsKstar";
    Configurable<float> peakMass{"peakMass", 0.89555f, "K*0(892) PDG mass (GeV)"};
    Configurable<float> sigma{"sigma", 0.0201f, "effective width (GeV) for the maxMassNSigma pre-cut -- keep in sync with hResonanceCorrelation.cxx's massWindowConfigurationsKstar.sigma"};
  } massWindowConfigurationsKstar;

  // Derived acceptance-window bounds, computed once in init() from the
  // Configurables above (peakMass +/- maxMassNSigma * sigma) rather than
  // hand-entered GeV numbers, so there is exactly one place to widen the
  // pre-cut instead of four independent literals that can drift apart.
  float mMinMassPhi = 0.f;
  float mMaxMassPhi = 0.f;
  float mMinMassKstar = 0.f;
  float mMaxMassKstar = 0.f;

  // For extracting strangeness mass QA plots
  struct : ConfigurableGroup {
    std::string prefix = "axesConfigurations";
    ConfigurableAxis axisPtQA{"axisPtQA", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "pt axis for QA histograms"};
    // Daughter/trigger-track QA axes (Phi->KK, K*0->Kpi, trigger hadrons):
    // each is a 1D QA histogram's own axis, booked in init().
    ConfigurableAxis axisEtaQA{"axisEtaQA", {100, -1.0f, 1.0f}, "#eta"};
    ConfigurableAxis axisNSigmaQA{"axisNSigmaQA", {100, -10.0f, 10.0f}, "n#sigma"};
    ConfigurableAxis axisDCAxyQA{"axisDCAxyQA", {200, -0.5f, 0.5f}, "DCA_{xy} (cm)"};
    ConfigurableAxis axisDCAzQA{"axisDCAzQA", {200, -0.5f, 0.5f}, "DCA_{z} (cm)"};
    ConfigurableAxis axisTPCCrossedRowsQA{"axisTPCCrossedRowsQA", {160, 0, 160}, "TPC crossed rows"};
    ConfigurableAxis axisVertexZQA{"axisVertexZQA", {60, -15.0f, 15.0f}, "Vertex Z (cm)"};
    ConfigurableAxis axisCentralityQA{"axisCentralityQA", {100, 0.0f, 100.0f}, "Centrality (%)"};
  } axesConfigurations;

  // QA
  struct : ConfigurableGroup {
    std::string prefix = "qaConfigurations";
    Configurable<bool> doTrueSelectionInMass{"doTrueSelectionInMass", false, "Fill mass histograms only with true primary Particles for MC"};
  } qaConfigurations;

  // Do declarative selections for DCAs, if possible
  Filter acceptanceFilter = (nabs(aod::track::eta) < configTracks.cfgCutEta && aod::track::pt > configTracks.cMinPtcut) &&
                            (nabs(aod::track::dcaXY) < configTracks.dcaXYconstant + configTracks.dcaXYpTdep * nabs(aod::track::signed1Pt)) &&
                            (nabs(aod::track::dcaZ) < configTracks.cMaxDCAzToPVcut);

  // All four aliases below join TrackSelection + TrackSelectionExtension so
  // that trackCut() -- the shared "is this a good, well-reconstructed track"
  // baseline -- is callable (and applied) on trigger and associated-particle
  // tracks exactly as it already is on Phi/K*0 daughters (TrackCandidates),
  // not just on the latter.
  using FullTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::TrackSelectionExtension>;
  using FullTracksMC = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::TrackSelectionExtension, aod::McTrackLabels>;
  using DauTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::TracksDCA>;
  using DauTracksMC = soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::TracksDCA, aod::McTrackLabels>;
  // using IDTracks= soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidBayesPi, aod::pidBayesKa, aod::pidBayesPr, aod::TOFSignal>; // prepared for Bayesian PID
  using IDTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPr, aod::pidTOFFullPr, aod::pidTPCFullEl, aod::pidTOFFullEl, aod::TOFSignal, aod::TracksDCA, aod::TrackSelection, aod::TrackSelectionExtension>;
  using IDTracksMC = soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPr, aod::pidTOFFullPr, aod::pidTPCFullEl, aod::pidTOFFullEl, aod::TOFSignal, aod::TracksDCA, aod::TrackSelection, aod::TrackSelectionExtension, aod::McTrackLabels>;

  Produces<aod::TriggerTracks> triggerTrack;
  Produces<aod::TriggerTrackExtras> triggerTrackExtra;
  Produces<aod::AssocPhis> assocPhis;
  Produces<aod::AssocKstars> assocKstars;
  Produces<aod::AssocHadrons> assocHadrons;
  Produces<aod::AssocPID> assocPID;

  using EventCandidates = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0As, aod::Mults>;
  using TrackCandidates = soa::Filtered<soa::Join<aod::FullTracks, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::TracksDCA, aod::TrackSelection, aod::TrackSelectionExtension>>;

  // for MC reco
  using MCEventCandidates = soa::Join<EventCandidates, aod::McCollisionLabels>;
  using MCTrackCandidates = soa::Filtered<soa::Join<TrackCandidates, aod::McTrackLabels>>;

  struct TriggCandidate {
    float pt = 0.f;
    int collisionId = -1;
    int trackId = -1;
    bool isPhysicalPrimary = false;
    float origPt = 0.f;
  };
  TriggCandidate thisTrigg;

  std::vector<TriggCandidate> triggerCandidates;

  void init(InitContext const&)
  {
    histos.add("CollCutCounts", "No. of event after cuts", kTH1I, {{10, 0, 10}});
    // Z-vertex and centrality of selected events (filled in isSelectedEvents()
    // once a collision passes every cut, i.e. the same "All Passed Events"
    // point as CollCutCounts bin 8). Centrality source for the latter is
    // configEvents.cfgUseFV0ACentrality: FT0C for O-O, FV0A for p-O.
    histos.add("hVertexZSelected", "Vertex Z of selected events", kTH1F, {axesConfigurations.axisVertexZQA});
    histos.add("hCentralitySelected", "Centrality of selected events", kTH1F, {axesConfigurations.axisCentralityQA});
    histos.get<TH1>(HIST("CollCutCounts"))->GetXaxis()->SetBinLabel(1, "All Events");
    histos.get<TH1>(HIST("CollCutCounts"))->GetXaxis()->SetBinLabel(2, "|Vz| < cut");
    histos.get<TH1>(HIST("CollCutCounts"))->GetXaxis()->SetBinLabel(3, "kIsTriggerTVX");
    histos.get<TH1>(HIST("CollCutCounts"))->GetXaxis()->SetBinLabel(4, "kNoTimeFrameBorder");
    histos.get<TH1>(HIST("CollCutCounts"))->GetXaxis()->SetBinLabel(5, "kNoITSROFrameBorder");
    histos.get<TH1>(HIST("CollCutCounts"))->GetXaxis()->SetBinLabel(6, "sel8");
    histos.get<TH1>(HIST("CollCutCounts"))->GetXaxis()->SetBinLabel(7, "kNoSameBunchPileup");
    histos.get<TH1>(HIST("CollCutCounts"))->GetXaxis()->SetBinLabel(8, "kIsGoodZvtxFT0vsPV");
    histos.get<TH1>(HIST("CollCutCounts"))->GetXaxis()->SetBinLabel(9, "All Passed Events");

    mMinMassPhi = massWindowConfigurationsPhi.peakMass - massWindowConfigurations.maxMassNSigma * massWindowConfigurationsPhi.sigma;
    mMaxMassPhi = massWindowConfigurationsPhi.peakMass + massWindowConfigurations.maxMassNSigma * massWindowConfigurationsPhi.sigma;
    mMinMassKstar = massWindowConfigurationsKstar.peakMass - massWindowConfigurations.maxMassNSigma * massWindowConfigurationsKstar.sigma;
    mMaxMassKstar = massWindowConfigurationsKstar.peakMass + massWindowConfigurations.maxMassNSigma * massWindowConfigurationsKstar.sigma;
    LOGF(info, "Assoc mass pre-cut windows: Phi [%.4f, %.4f] GeV, K*0 [%.4f, %.4f] GeV",
         mMinMassPhi, mMaxMassPhi, mMinMassKstar, mMaxMassKstar);

    // QA histograms for Phi (K+K-) and K*0 (K-pi) daughter tracks: booked once
    // per daughter species (Phi's kaon, Kstar's kaon, Kstar's pion), filled at
    // the point each candidate passes all cuts (including the mass window) in
    // processPhis/processKstars (and their MC counterparts).
    auto bookDaughterQA = [&](const std::string& dir) {
      histos.add((dir + "/hPt").c_str(), "p_{T}", kTH1F, {axesConfigurations.axisPtQA});
      histos.add((dir + "/hEta").c_str(), "#eta", kTH1F, {axesConfigurations.axisEtaQA});
      histos.add((dir + "/hTPCNSigma").c_str(), "TPC n#sigma", kTH1F, {axesConfigurations.axisNSigmaQA});
      histos.add((dir + "/hTOFNSigma").c_str(), "TOF n#sigma", kTH1F, {axesConfigurations.axisNSigmaQA});
      histos.add((dir + "/hDCAxy").c_str(), "DCA_{xy}", kTH1F, {axesConfigurations.axisDCAxyQA});
      histos.add((dir + "/hDCAz").c_str(), "DCA_{z}", kTH1F, {axesConfigurations.axisDCAzQA});
      histos.add((dir + "/hTPCCrossedRows").c_str(), "TPC crossed rows", kTH1F, {axesConfigurations.axisTPCCrossedRowsQA});
      // 2D pT-binned QA, same variables as above, to see pT dependence
      histos.add((dir + "/hPtVsDCAxy").c_str(), "p_{T} vs DCA_{xy}", kTH2F, {axesConfigurations.axisPtQA, axesConfigurations.axisDCAxyQA});
      histos.add((dir + "/hPtVsDCAz").c_str(), "p_{T} vs DCA_{z}", kTH2F, {axesConfigurations.axisPtQA, axesConfigurations.axisDCAzQA});
      histos.add((dir + "/hPtVsTPCNSigma").c_str(), "p_{T} vs TPC n#sigma", kTH2F, {axesConfigurations.axisPtQA, axesConfigurations.axisNSigmaQA});
      histos.add((dir + "/hPtVsTOFNSigma").c_str(), "p_{T} vs TOF n#sigma", kTH2F, {axesConfigurations.axisPtQA, axesConfigurations.axisNSigmaQA});
      // TPC vs TOF nSigma correlation, for the usual PID "banana plot" check
      histos.add((dir + "/hTPCNSigmaVsTOFNSigma").c_str(), "TPC n#sigma vs TOF n#sigma", kTH2F, {axesConfigurations.axisNSigmaQA, axesConfigurations.axisNSigmaQA});
    };
    bookDaughterQA("QA/Phi/Kaon");
    bookDaughterQA("QA/Kstar/Kaon");
    bookDaughterQA("QA/Kstar/Pion");

    // Trigger hadron QA: no PID applied to trigger tracks, so only the 5
    // non-PID variables are booked (pt, eta, dcaXY, dcaZ, TPC crossed rows),
    // plus the pT-binned 2D versions of DCAxy/DCAz (no nSigma -- no PID here).
    histos.add("QA/TriggerHadron/hPt", "p_{T}", kTH1F, {axesConfigurations.axisPtQA});
    histos.add("QA/TriggerHadron/hEta", "#eta", kTH1F, {axesConfigurations.axisEtaQA});
    histos.add("QA/TriggerHadron/hDCAxy", "DCA_{xy}", kTH1F, {axesConfigurations.axisDCAxyQA});
    histos.add("QA/TriggerHadron/hDCAz", "DCA_{z}", kTH1F, {axesConfigurations.axisDCAzQA});
    histos.add("QA/TriggerHadron/hTPCCrossedRows", "TPC crossed rows", kTH1F, {axesConfigurations.axisTPCCrossedRowsQA});
    histos.add("QA/TriggerHadron/hPtVsDCAxy", "p_{T} vs DCA_{xy}", kTH2F, {axesConfigurations.axisPtQA, axesConfigurations.axisDCAxyQA});
    histos.add("QA/TriggerHadron/hPtVsDCAz", "p_{T} vs DCA_{z}", kTH2F, {axesConfigurations.axisPtQA, axesConfigurations.axisDCAzQA});
  }

  // HIST() needs a compile-time literal path, so each daughter species gets
  // its own explicit fill function rather than one generic/parameterized
  // helper. Kaon QA uses the Ka PID accessors, pion QA uses the Pi ones.
  template <typename Track>
  void fillPhiKaonQA(const Track& track)
  {
    histos.fill(HIST("QA/Phi/Kaon/hPt"), track.pt());
    histos.fill(HIST("QA/Phi/Kaon/hEta"), track.eta());
    histos.fill(HIST("QA/Phi/Kaon/hTPCNSigma"), track.tpcNSigmaKa());
    histos.fill(HIST("QA/Phi/Kaon/hTOFNSigma"), track.tofNSigmaKa());
    histos.fill(HIST("QA/Phi/Kaon/hDCAxy"), track.dcaXY());
    histos.fill(HIST("QA/Phi/Kaon/hDCAz"), track.dcaZ());
    histos.fill(HIST("QA/Phi/Kaon/hTPCCrossedRows"), track.tpcNClsCrossedRows());
    histos.fill(HIST("QA/Phi/Kaon/hPtVsDCAxy"), track.pt(), track.dcaXY());
    histos.fill(HIST("QA/Phi/Kaon/hPtVsDCAz"), track.pt(), track.dcaZ());
    histos.fill(HIST("QA/Phi/Kaon/hPtVsTPCNSigma"), track.pt(), track.tpcNSigmaKa());
    histos.fill(HIST("QA/Phi/Kaon/hPtVsTOFNSigma"), track.pt(), track.tofNSigmaKa());
    if (track.hasTOF()) {
      histos.fill(HIST("QA/Phi/Kaon/hTPCNSigmaVsTOFNSigma"), track.tpcNSigmaKa(), track.tofNSigmaKa());
    }
  }

  template <typename Track>
  void fillKstarKaonQA(const Track& track)
  {
    histos.fill(HIST("QA/Kstar/Kaon/hPt"), track.pt());
    histos.fill(HIST("QA/Kstar/Kaon/hEta"), track.eta());
    histos.fill(HIST("QA/Kstar/Kaon/hTPCNSigma"), track.tpcNSigmaKa());
    histos.fill(HIST("QA/Kstar/Kaon/hTOFNSigma"), track.tofNSigmaKa());
    histos.fill(HIST("QA/Kstar/Kaon/hDCAxy"), track.dcaXY());
    histos.fill(HIST("QA/Kstar/Kaon/hDCAz"), track.dcaZ());
    histos.fill(HIST("QA/Kstar/Kaon/hTPCCrossedRows"), track.tpcNClsCrossedRows());
    histos.fill(HIST("QA/Kstar/Kaon/hPtVsDCAxy"), track.pt(), track.dcaXY());
    histos.fill(HIST("QA/Kstar/Kaon/hPtVsDCAz"), track.pt(), track.dcaZ());
    histos.fill(HIST("QA/Kstar/Kaon/hPtVsTPCNSigma"), track.pt(), track.tpcNSigmaKa());
    histos.fill(HIST("QA/Kstar/Kaon/hPtVsTOFNSigma"), track.pt(), track.tofNSigmaKa());
    if (track.hasTOF()) {
      histos.fill(HIST("QA/Kstar/Kaon/hTPCNSigmaVsTOFNSigma"), track.tpcNSigmaKa(), track.tofNSigmaKa());
    }
  }

  template <typename Track>
  void fillKstarPionQA(const Track& track)
  {
    histos.fill(HIST("QA/Kstar/Pion/hPt"), track.pt());
    histos.fill(HIST("QA/Kstar/Pion/hEta"), track.eta());
    histos.fill(HIST("QA/Kstar/Pion/hTPCNSigma"), track.tpcNSigmaPi());
    histos.fill(HIST("QA/Kstar/Pion/hTOFNSigma"), track.tofNSigmaPi());
    histos.fill(HIST("QA/Kstar/Pion/hDCAxy"), track.dcaXY());
    histos.fill(HIST("QA/Kstar/Pion/hDCAz"), track.dcaZ());
    histos.fill(HIST("QA/Kstar/Pion/hTPCCrossedRows"), track.tpcNClsCrossedRows());
    histos.fill(HIST("QA/Kstar/Pion/hPtVsDCAxy"), track.pt(), track.dcaXY());
    histos.fill(HIST("QA/Kstar/Pion/hPtVsDCAz"), track.pt(), track.dcaZ());
    histos.fill(HIST("QA/Kstar/Pion/hPtVsTPCNSigma"), track.pt(), track.tpcNSigmaPi());
    histos.fill(HIST("QA/Kstar/Pion/hPtVsTOFNSigma"), track.pt(), track.tofNSigmaPi());
    if (track.hasTOF()) {
      histos.fill(HIST("QA/Kstar/Pion/hTPCNSigmaVsTOFNSigma"), track.tpcNSigmaPi(), track.tofNSigmaPi());
    }
  }

  // Trigger hadrons carry no PID selection, so only the 5 non-PID variables
  // are filled here.
  template <typename Track>
  void fillTriggerHadronQA(const Track& track)
  {
    histos.fill(HIST("QA/TriggerHadron/hPt"), track.pt());
    histos.fill(HIST("QA/TriggerHadron/hEta"), track.eta());
    histos.fill(HIST("QA/TriggerHadron/hDCAxy"), track.dcaXY());
    histos.fill(HIST("QA/TriggerHadron/hDCAz"), track.dcaZ());
    histos.fill(HIST("QA/TriggerHadron/hTPCCrossedRows"), track.tpcNClsCrossedRows());
    histos.fill(HIST("QA/TriggerHadron/hPtVsDCAxy"), track.pt(), track.dcaXY());
    histos.fill(HIST("QA/TriggerHadron/hPtVsDCAz"), track.pt(), track.dcaZ());
  }

  template <typename Coll>
  bool isSelectedEvents(const Coll& collision, bool fillHist = true)
  {
    auto applyCut = [&](bool enabled, bool condition, int bin) {
      if (!enabled) {
        return true;
      }
      if (!condition) {
        return false;
      }
      if (fillHist) {
        histos.fill(HIST("CollCutCounts"), bin);
      }
      return true;
    };

    if (fillHist) {
      histos.fill(HIST("CollCutCounts"), 0);
    }

    if (!applyCut(true, std::abs(collision.posZ()) <= configEvents.cfgEvtZvtx, 1)) {
      return false;
    }

    if (!applyCut(configEvents.cfgEvtTriggerTVXSel,
                  collision.selection_bit(aod::evsel::kIsTriggerTVX), 2)) {
      return false;
    }

    if (!applyCut(configEvents.cfgEvtNoTFBorderCut,
                  collision.selection_bit(aod::evsel::kNoTimeFrameBorder), 3)) {
      return false;
    }

    if (!applyCut(configEvents.cfgEvtNoITSROFrameBorderCut,
                  collision.selection_bit(aod::evsel::kNoITSROFrameBorder), 4)) {
      return false;
    }

    if (!applyCut(configEvents.cfgEvtSel8, collision.sel8(), 5)) {
      return false;
    }

    if (!applyCut(configEvents.cfgEvtNoSameBunchPileupCut,
                  collision.selection_bit(aod::evsel::kNoSameBunchPileup), 6)) {
      return false;
    }

    if (!applyCut(configEvents.cfgEvtGoodZvtxFT0vsPVCut,
                  collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV), 7)) {
      return false;
    }

    if (fillHist) {
      histos.fill(HIST("CollCutCounts"), 8);
      histos.fill(HIST("hVertexZSelected"), collision.posZ());
      // Guarded with if constexpr: isSelectedEvents() is called with several
      // different collision join types (some callers pass fillHist=false and
      // don't join centrality tables at all), so this only compiles/fills for
      // callers whose collision type actually has both estimators joined.
      if constexpr (requires { collision.centFT0C(); collision.centFV0A(); }) {
        histos.fill(HIST("hCentralitySelected"), configEvents.cfgUseFV0ACentrality ? collision.centFV0A() : collision.centFT0C());
      }
    }

    return true;
  }

  template <typename TrackType>
  bool trackCut(const TrackType& track)
  {
    // basic track cuts
    // NOTE: pT-dependent |DCAxy| cut is NOT re-applied here -- it's already
    // enforced unconditionally, for every track in every process function,
    // by acceptanceFilter (Filter on aod::track::dcaXY using the live
    // configTracks.dcaXYconstant/dcaXYpTdep). A second check here would be
    // either a no-op (same cut) or a stale one (if those configurables are
    // ever retuned and this duplicate doesn't track them).
    if (configTracks.cfgGlobalWoDCATrack && !track.isGlobalTrackWoDCA()) {
      return false;
    }
    if (configTracks.cfgPVContributor && !track.isPVContributor()) {
      return false;
    }
    if (track.tpcNClsCrossedRows() < configTracks.cfgMinCrossedRows) {
      return false;
    }
    if (track.tpcChi2NCl() > configTracks.cfgMaxTPCChi2NCl) {
      return false;
    }
    if (track.itsChi2NCl() > configTracks.cfgMaxITSChi2NCl) {
      return false;
    }
    if (track.tpcCrossedRowsOverFindableCls() < configTracks.cfgMinTPCCrossedRowsOverFindableCls) {
      return false;
    }
    if (track.itsNCls() < configTracks.cfgMinITSNCls) {
      return false;
    }

    return true;
  }

  template <typename T>
  bool ptDependentPidKaon(const T& candidate)
  {
    auto vKaonTPCPIDpTintv = configPID.kaonTPCPIDpTintv.value;
    vKaonTPCPIDpTintv.insert(vKaonTPCPIDpTintv.begin(), configTracks.cMinPtcut);
    auto vKaonTPCPIDcuts = configPID.kaonTPCPIDcuts.value;
    auto vKaonTOFPIDpTintv = configPID.kaonTOFPIDpTintv.value;
    auto vKaonTPCTOFCombinedpTintv = configPID.kaonTPCTOFCombinedpTintv.value;
    auto vKaonTPCTOFCombinedPIDcuts = configPID.kaonTPCTOFCombinedPIDcuts.value;
    auto vKaonTOFPIDcuts = configPID.kaonTOFPIDcuts.value;

    float pt = candidate.pt();
    float ptSwitchToTOF = vKaonTPCPIDpTintv.back();
    float tpcNsigmaKa = candidate.tpcNSigmaKa();
    float tofNsigmaKa = candidate.tofNSigmaKa();

    bool tpcPIDPassed = false;

    // TPC PID interval-based check
    for (size_t i = 0; i < vKaonTPCPIDpTintv.size() - 1; ++i) {
      if (pt > vKaonTPCPIDpTintv[i] && pt < vKaonTPCPIDpTintv[i + 1]) {
        if (std::abs(tpcNsigmaKa) < vKaonTPCPIDcuts[i]) {
          tpcPIDPassed = true;
          break;
        }
      }
    }

    // TOF bypass option
    if (configPID.cByPassTOF) {
      return std::abs(tpcNsigmaKa) < vKaonTPCPIDcuts.back();
    }

    // Case 1: No TOF and pt ≤ ptSwitch → use TPC-only
    if (!candidate.hasTOF() && pt <= ptSwitchToTOF) {
      return tpcPIDPassed;
    }

    // Case 2: No TOF but pt > ptSwitch → reject
    if (!candidate.hasTOF() && pt > ptSwitchToTOF) {
      return false;
    }

    // Case 3: TOF is available → apply TPC+TOF PID logic
    if (candidate.hasTOF()) {
      if (configPID.cPIDcutType == SquareType) {
        // Rectangular cut
        for (size_t i = 0; i < vKaonTOFPIDpTintv.size(); ++i) {
          if (pt < vKaonTOFPIDpTintv[i]) {
            if (std::abs(tofNsigmaKa) < vKaonTOFPIDcuts[i] &&
                std::abs(tpcNsigmaKa) < vKaonTPCPIDcuts.back()) {
              return true;
            }
          }
        }
      } else if (configPID.cPIDcutType == CircularType) {
        // Circular cut
        for (size_t i = 0; i < vKaonTPCTOFCombinedpTintv.size(); ++i) {
          if (pt < vKaonTPCTOFCombinedpTintv[i]) {
            float combinedSigma2 = tpcNsigmaKa * tpcNsigmaKa +
                                   tofNsigmaKa * tofNsigmaKa;
            if (combinedSigma2 < vKaonTPCTOFCombinedPIDcuts[i] * vKaonTPCTOFCombinedPIDcuts[i]) {
              return true;
            }
          }
        }
      }
    }

    return false;
  }

  template <typename T>
  bool selectionPID(const T& candidate)
  {
    auto vKaonTPCPIDcuts = configPID.kaonTPCPIDcuts.value;
    auto vKaonTPCTOFCombinedPIDcuts = configPID.kaonTPCTOFCombinedPIDcuts.value;

    if (!configPID.cByPassTOF && candidate.hasTOF() && (candidate.tofNSigmaKa() * candidate.tofNSigmaKa() + candidate.tpcNSigmaKa() * candidate.tpcNSigmaKa()) < (vKaonTPCTOFCombinedPIDcuts[0] * vKaonTPCTOFCombinedPIDcuts[0])) {
      return true;
    }
    if (!configPID.cByPassTOF && !candidate.hasTOF() && std::abs(candidate.tpcNSigmaKa()) < vKaonTPCPIDcuts[0]) {
      return true;
    }
    if (configPID.cByPassTOF && std::abs(candidate.tpcNSigmaKa()) < vKaonTPCPIDcuts[0]) {
      return true;
    }
    return false;
  }

  template <typename T>
  bool ptDependentPidPion(const T& candidate)
  {
    auto vPionTPCPIDpTintv = configPID.pionTPCPIDpTintv.value;
    vPionTPCPIDpTintv.insert(vPionTPCPIDpTintv.begin(), configTracks.cMinPtcut);
    auto vPionTPCPIDcuts = configPID.pionTPCPIDcuts.value;
    auto vPionTOFPIDpTintv = configPID.pionTOFPIDpTintv.value;
    auto vPionTPCTOFCombinedpTintv = configPID.pionTPCTOFCombinedpTintv.value;
    auto vPionTPCTOFCombinedPIDcuts = configPID.pionTPCTOFCombinedPIDcuts.value;
    auto vPionTOFPIDcuts = configPID.pionTOFPIDcuts.value;

    float pt = candidate.pt();
    float ptSwitchToTOF = vPionTPCPIDpTintv.back();
    float tpcNsigmaPi = candidate.tpcNSigmaPi();
    float tofNsigmaPi = candidate.tofNSigmaPi();

    bool tpcPIDPassed = false;

    // TPC PID interval-based check
    for (size_t i = 0; i < vPionTPCPIDpTintv.size() - 1; ++i) {
      if (pt > vPionTPCPIDpTintv[i] && pt < vPionTPCPIDpTintv[i + 1]) {
        if (std::abs(tpcNsigmaPi) < vPionTPCPIDcuts[i]) {
          tpcPIDPassed = true;
          break;
        }
      }
    }

    // TOF bypass option
    if (configPID.cByPassTOF) {
      return std::abs(tpcNsigmaPi) < vPionTPCPIDcuts.back();
    }

    // Case 1: No TOF and pt ≤ ptSwitch → use TPC-only
    if (!candidate.hasTOF() && pt <= ptSwitchToTOF) {
      return tpcPIDPassed;
    }

    // Case 2: No TOF but pt > ptSwitch → reject
    if (!candidate.hasTOF() && pt > ptSwitchToTOF) {
      return false;
    }

    // Case 3: TOF is available → apply TPC+TOF PID logic
    if (candidate.hasTOF()) {
      if (configPID.cPIDcutType == SquareType) {
        // Rectangular cut
        for (size_t i = 0; i < vPionTOFPIDpTintv.size(); ++i) {
          if (pt < vPionTOFPIDpTintv[i]) {
            if (std::abs(tofNsigmaPi) < vPionTOFPIDcuts[i] &&
                std::abs(tpcNsigmaPi) < vPionTPCPIDcuts.back()) {
              return true;
            }
          }
        }
      } else if (configPID.cPIDcutType == CircularType) {
        // Circular cut
        for (size_t i = 0; i < vPionTPCTOFCombinedpTintv.size(); ++i) {
          if (pt < vPionTPCTOFCombinedpTintv[i]) {
            float combinedSigma2 = tpcNsigmaPi * tpcNsigmaPi +
                                   tofNsigmaPi * tofNsigmaPi;
            if (combinedSigma2 < vPionTPCTOFCombinedPIDcuts[i] * vPionTPCTOFCombinedPIDcuts[i]) {
              return true;
            }
          }
        }
      }
    }

    return false;
  }

  template <typename T>
  bool selectionPIDPion(const T& candidate)
  {
    auto vPionTPCPIDcuts = configPID.pionTPCPIDcuts.value;
    auto vPionTPCTOFCombinedPIDcuts = configPID.pionTPCTOFCombinedPIDcuts.value;

    if (!configPID.cByPassTOF && candidate.hasTOF() && (candidate.tofNSigmaPi() * candidate.tofNSigmaPi() + candidate.tpcNSigmaPi() * candidate.tpcNSigmaPi()) < (vPionTPCTOFCombinedPIDcuts[0] * vPionTPCTOFCombinedPIDcuts[0])) {
      return true;
    }
    if (!configPID.cByPassTOF && !candidate.hasTOF() && std::abs(candidate.tpcNSigmaPi()) < vPionTPCPIDcuts[0]) {
      return true;
    }
    if (configPID.cByPassTOF && std::abs(candidate.tpcNSigmaPi()) < vPionTPCPIDcuts[0]) {
      return true;
    }
    return false;
  }

  // reco-level trigger quality checks (N.B.: DCA is filtered, not selected)
  template <class TTrack>
  bool isValidTrigger(TTrack const& track)
  {
    // Shared "good track" baseline, applied identically across trigger,
    // associated-particle, and resonance-daughter tracks (see trackCut()).
    if (!trackCut(track)) {
      return false;
    }
    if (track.eta() > configTriggerTracks.triggerEtaMax || track.eta() < configTriggerTracks.triggerEtaMin) {
      return false;
    }
    // if (track.sign()= 1 ) {continue;}
    if (track.pt() > configTriggerTracks.triggerPtCutMax || track.pt() < configTriggerTracks.triggerPtCutMin) {
      return false;
    }
    if (track.tpcNClsCrossedRows() < configTracks.minTPCNCrossedRows) {
      return false; // crossed rows
    }
    if (!track.hasITS() && configTriggerTracks.triggerRequireITS) {
      return false; // skip, doesn't have ITS signal (skips lots of TPC-only!)
    }
    if (track.tpcNClsShared() > configTriggerTracks.triggerMaxTPCSharedClusters) {
      return false; // skip, has shared clusters
    }
    if (!(BIT_CHECK(track.itsClusterMap(), 0)) && configTriggerTracks.triggerRequireL0) {
      return false; // skip, doesn't have cluster in ITS L0
    }
    return true;
  }

  // Species tag for isValidAssocTrack<Species>()'s PID gate below. AssocPion
  // and AssocKaon select which nSigma cut direction applies (accept that
  // species, reject the other two); AssocHadron carries no PID logic of its
  // own -- passed for clarity at Hadron call sites, but never actually reaches
  // the `if constexpr (Species == AssocKaon) ... else ...` branch, since the
  // outer `requires { assoc.tofSignal(); }` is already false for FullTracks
  // (Hadron's track type, which has no PID columns at all).
  enum AssocSpecies { AssocPion = 0,
                      AssocKaon = 1,
                      AssocHadron = 2 };

  // Merged predicate for the Pion/Kaon/Hadron associated-track pools. Pion and
  // Kaon both run over IDTracks/IDTracksMC (identical C++ type), so the two
  // PID directions cannot be told apart by a `requires{}` type check the way
  // Hadron (FullTracks, no PID columns) is told apart from the other two --
  // hence the explicit non-type template parameter instead. Everything except
  // the PID nSigma cut direction (phase-space cuts, MC-truth lookup, table
  // fills) is identical across all three species, so it stays unduplicated
  // here; only the accept/reject block branches on Species.
  template <int Species, class TTrack>
  bool isValidAssocTrack(TTrack const& assoc)
  {
    static_assert(Species == AssocPion || Species == AssocKaon || Species == AssocHadron,
                  "isValidAssocTrack: unknown species tag");
    // Shared "good track" baseline, applied identically across trigger,
    // associated-particle, and resonance-daughter tracks (see trackCut()).
    if (!trackCut(assoc)) {
      return false;
    }
    if (assoc.eta() > configAssocTracks.assocEtaMax || assoc.eta() < configAssocTracks.assocEtaMin) {
      return false;
    }
    if (assoc.pt() > configAssocTracks.assocPtCutMax || assoc.pt() < configAssocTracks.assocPtCutMin) {
      return false;
    }
    if (assoc.tpcNClsCrossedRows() < configTracks.minTPCNCrossedRows) {
      return false; // crossed rows
    }
    if (!assoc.hasITS() && configAssocTracks.assocRequireITS) {
      return false; // skip, doesn't have ITS signal (skips lots of TPC-only!)
    }

    // do this only if information is available
    float nSigmaTPCTOF[8] = {-10, -10, -10, -10, -10, -10, -10, -10};
    if constexpr (requires { assoc.tofSignal(); } && !requires { assoc.mcParticle(); }) {
      if (assoc.tofSignal() > 0) {
        if constexpr (Species == AssocKaon) {
          if (std::sqrt(assoc.tofNSigmaKa() * assoc.tofNSigmaKa() + assoc.tpcNSigmaKa() * assoc.tpcNSigmaKa()) > configPID.assocKaonNSigmaTPCFOF)
            return false;
          if (assoc.tofNSigmaPr() < configPID.rejectSigma)
            return false;
          if (assoc.tpcNSigmaPr() < configPID.rejectSigma)
            return false;
          if (assoc.tofNSigmaPi() < configPID.rejectSigma)
            return false;
          if (assoc.tpcNSigmaPi() < configPID.rejectSigma)
            return false;
        } else {
          if (std::sqrt(assoc.tofNSigmaPi() * assoc.tofNSigmaPi() + assoc.tpcNSigmaPi() * assoc.tpcNSigmaPi()) > configPID.assocPionNSigmaTPCFOF)
            return false;
          if (assoc.tofNSigmaPr() < configPID.rejectSigma)
            return false;
          if (assoc.tpcNSigmaPr() < configPID.rejectSigma)
            return false;
          if (assoc.tofNSigmaKa() < configPID.rejectSigma)
            return false;
          if (assoc.tpcNSigmaKa() < configPID.rejectSigma)
            return false;
        }
        nSigmaTPCTOF[4] = assoc.tofNSigmaPi();
        nSigmaTPCTOF[5] = assoc.tofNSigmaKa();
        nSigmaTPCTOF[6] = assoc.tofNSigmaPr();
        nSigmaTPCTOF[7] = assoc.tofNSigmaEl();
      } else {
        if constexpr (Species == AssocKaon) {
          if (assoc.tpcNSigmaKa() > configPID.assocKaonNSigmaTPCFOF)
            return false;
          if (assoc.tpcNSigmaPr() < configPID.rejectSigma)
            return false;
          if (assoc.tpcNSigmaPi() < configPID.rejectSigma)
            return false;
        } else {
          if (assoc.tpcNSigmaPi() > configPID.assocPionNSigmaTPCFOF)
            return false;
          if (assoc.tpcNSigmaPr() < configPID.rejectSigma)
            return false;
          if (assoc.tpcNSigmaKa() < configPID.rejectSigma)
            return false;
        }
      }
      nSigmaTPCTOF[0] = assoc.tpcNSigmaPi();
      nSigmaTPCTOF[1] = assoc.tpcNSigmaKa();
      nSigmaTPCTOF[2] = assoc.tpcNSigmaPr();
      nSigmaTPCTOF[3] = assoc.tpcNSigmaEl();
    }

    bool physicalPrimary = false;
    float origPt = -1;
    float code = -9999;
    if constexpr (requires { assoc.mcParticle(); }) {
      if (assoc.has_mcParticle()) {
        auto mcParticle = assoc.mcParticle();
        physicalPrimary = mcParticle.isPhysicalPrimary();
        origPt = mcParticle.pt();
        code = mcParticle.pdgCode();
      }
    }

    // Rapidity under the species mass hypothesis; AssocHadron has no mass
    // hypothesis to assume, so it's filled with a dummy value instead.
    float rapidity = -999.f;
    if constexpr (Species == AssocKaon) {
      rapidity = LorentzVectorPtEtaPhiMass(assoc.pt(), assoc.eta(), assoc.phi(),
                                           o2::constants::physics::MassKPlus)
                   .Rapidity();
    } else if constexpr (Species == AssocPion) {
      rapidity = LorentzVectorPtEtaPhiMass(assoc.pt(), assoc.eta(), assoc.phi(),
                                           o2::constants::physics::MassPiPlus)
                   .Rapidity();
    }

    assocHadrons(
      assoc.collisionId(),
      physicalPrimary,
      assoc.globalIndex(),
      origPt,
      code,
      rapidity);
    assocPID(
      nSigmaTPCTOF[0],
      nSigmaTPCTOF[1],
      nSigmaTPCTOF[2],
      nSigmaTPCTOF[3],
      nSigmaTPCTOF[4],
      nSigmaTPCTOF[5],
      nSigmaTPCTOF[6],
      nSigmaTPCTOF[7]);
    return true;
  }

  // for real data processing
  void processTriggers(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFV0As, aod::PVMults>::iterator const& collision, soa::Filtered<FullTracks> const& tracks)
  {
    triggerCandidates.clear();
    if (!isSelectedEvents(collision)) {
      return;
    }

    /// _________________________________________________
    /// Step 1: Populate table with trigger tracks
    double leadingPt = -1.;
    int leadingId = -1;
    for (auto const& track : tracks) {
      if (!isValidTrigger(track))
        continue;
      fillTriggerHadronQA(track);
      thisTrigg.pt = track.pt();
      thisTrigg.trackId = track.globalIndex();
      thisTrigg.collisionId = track.collisionId();
      thisTrigg.isPhysicalPrimary = false; // if you decide to check real data for primaries, you'll have a hard time
      thisTrigg.origPt = 0;
      triggerCandidates.push_back(thisTrigg);
      if (track.pt() > leadingPt) {
        leadingPt = track.pt();
        leadingId = track.globalIndex();
      }
    }
    for (auto const& TriggCandidate : triggerCandidates) {
      bool isLeading = (leadingId == TriggCandidate.trackId);
      triggerTrack(
        TriggCandidate.collisionId,
        TriggCandidate.isPhysicalPrimary,
        TriggCandidate.trackId,
        TriggCandidate.origPt,
        isLeading);
      triggerTrackExtra(1);
    }
  }

  // for MC processing
  void processTriggersMC(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFV0As, aod::PVMults>::iterator const& collision, soa::Filtered<FullTracksMC> const& tracks, aod::McParticles const&)
  {
    triggerCandidates.clear();
    if (!isSelectedEvents(collision)) {
      return;
    }

    /// _________________________________________________
    /// Step 1: Populate table with trigger tracks
    double leadingPt = -1.;
    int leadingId = -1;
    for (auto const& track : tracks) {
      if (!isValidTrigger(track))
        continue;
      fillTriggerHadronQA(track);
      thisTrigg.pt = track.pt();
      thisTrigg.trackId = track.globalIndex();
      thisTrigg.collisionId = track.collisionId();
      if (track.has_mcParticle()) {
        auto mcParticle = track.mcParticle();
        thisTrigg.isPhysicalPrimary = mcParticle.isPhysicalPrimary();
        thisTrigg.origPt = mcParticle.pt();
      }
      triggerCandidates.push_back(thisTrigg);
      if (track.pt() > leadingPt) {
        leadingPt = track.pt();
        leadingId = track.globalIndex();
      }
    }

    for (auto const& TriggCandidate : triggerCandidates) {
      bool isLeading = (leadingId == TriggCandidate.trackId);
      triggerTrack(
        TriggCandidate.collisionId,
        TriggCandidate.isPhysicalPrimary,
        TriggCandidate.trackId,
        TriggCandidate.origPt,
        isLeading);
      triggerTrackExtra(1);
    }
  }

  void processAssocPions(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, soa::Filtered<IDTracks> const& tracks)
  {
    // Perform basic event selection (fillHist=false: processTriggers/processTriggersMC own the shared CollCutCounts fill)
    if (!isSelectedEvents(collision, false)) {
      return;
    }

    /// _________________________________________________
    /// Step 1: Populate table with trigger tracks
    for (auto const& track : tracks) {
      if (!isValidAssocTrack<AssocPion>(track))
        continue;
    }
  }

  void processAssocPionsMC(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, soa::Filtered<IDTracksMC> const& tracks, aod::McParticles const&)
  {
    // Perform basic event selection (fillHist=false: processTriggers/processTriggersMC own the shared CollCutCounts fill)
    if (!isSelectedEvents(collision, false)) {
      return;
    }

    /// _________________________________________________
    /// Step 1: Populate table with trigger tracks
    for (auto const& track : tracks) {
      if (!isValidAssocTrack<AssocPion>(track))
        continue;
    }
  }

  void processAssocKaons(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, soa::Filtered<IDTracks> const& tracks)
  {
    // Perform basic event selection (fillHist=false: processTriggers/processTriggersMC own the shared CollCutCounts fill)
    if (!isSelectedEvents(collision, false)) {
      return;
    }

    /// _________________________________________________
    /// Step 1: Populate table with trigger tracks
    for (auto const& track : tracks) {
      if (!isValidAssocTrack<AssocKaon>(track))
        continue;
    }
  }

  void processAssocKaonsMC(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, soa::Filtered<IDTracksMC> const& tracks, aod::McParticles const&)
  {
    // Perform basic event selection (fillHist=false: processTriggers/processTriggersMC own the shared CollCutCounts fill)
    if (!isSelectedEvents(collision, false)) {
      return;
    }

    /// _________________________________________________
    /// Step 1: Populate table with trigger tracks
    for (auto const& track : tracks) {
      if (!isValidAssocTrack<AssocKaon>(track))
        continue;
    }
  }

  void processAssocHadrons(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, soa::Filtered<FullTracks> const& tracks)
  {
    // Perform basic event selection (fillHist=false: processTriggers/processTriggersMC own the shared CollCutCounts fill)
    if (!isSelectedEvents(collision, false)) {
      return;
    }

    /// _________________________________________________
    /// Step 1: Populate table with trigger tracks
    for (auto const& track : tracks) {
      if (!isValidAssocTrack<AssocHadron>(track))
        continue;
    }
  }
  void processAssocHadronsMC(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, soa::Filtered<FullTracksMC> const& tracks, aod::McParticles const&)
  {
    // Perform basic event selection (fillHist=false: processTriggers/processTriggersMC own the shared CollCutCounts fill)
    if (!isSelectedEvents(collision, false)) {
      return;
    }

    /// _________________________________________________
    /// Step 1: Populate table with trigger tracks
    for (auto const& track : tracks) {
      if (!isValidAssocTrack<AssocHadron>(track))
        continue;
    }
  }

  void processPhis(EventCandidates::iterator const& collision,
                   TrackCandidates const& tracks)
  {
    if (!isSelectedEvents(collision, false)) {
      return;
    }

    LorentzVectorPtEtaPhiMass kplus;
    LorentzVectorPtEtaPhiMass kminus;
    LorentzVectorPtEtaPhiMass phi;

    for (auto const& [trk1, trk2] :
         combinations(CombinationsFullIndexPolicy(tracks, tracks))) {

      if (trk2.index() <= trk1.index()) {
        continue;
      }

      // same collision
      if (trk1.collisionId() != collision.globalIndex() ||
          trk2.collisionId() != collision.globalIndex()) {
        continue;
      }

      // quality cuts
      if (!trackCut(trk1) || !trackCut(trk2)) {
        continue;
      }

      // unlike-sign only
      if (trk1.sign() * trk2.sign() > 0) {
        continue;
      }

      // kaon PID
      if (configPID.ispTdepPID) {
        if (!ptDependentPidKaon(trk1) ||
            !ptDependentPidKaon(trk2)) {
          continue;
        }
      } else {
        if (!selectionPID(trk1) ||
            !selectionPID(trk2)) {
          continue;
        }
      }

      // assign K+ and K-
      if (trk1.sign() > 0) {
        kplus = LorentzVectorPtEtaPhiMass(trk1.pt(), trk1.eta(), trk1.phi(),
                                          o2::constants::physics::MassKPlus);
        kminus = LorentzVectorPtEtaPhiMass(trk2.pt(), trk2.eta(), trk2.phi(),
                                           o2::constants::physics::MassKPlus);
      } else {
        kplus = LorentzVectorPtEtaPhiMass(trk2.pt(), trk2.eta(), trk2.phi(),
                                          o2::constants::physics::MassKPlus);
        kminus = LorentzVectorPtEtaPhiMass(trk1.pt(), trk1.eta(), trk1.phi(),
                                           o2::constants::physics::MassKPlus);
      }

      phi = kplus + kminus;

      float invMass = phi.M();

      // rapidity cut
      if (std::abs(phi.Rapidity()) > configResoDauTracks.cfgCutRapidity) {
        continue;
      }

      // optional mass window
      if (invMass < mMinMassPhi || invMass > mMaxMassPhi) {
        continue;
      }

      // Daughter QA (both legs are kaons for Phi)
      fillPhiKaonQA(trk1);
      fillPhiKaonQA(trk2);

      assocPhis(
        collision.globalIndex(),
        false,
        false,
        phi.Pt(),
        phi.Eta(),
        phi.Phi(),
        invMass,
        phi.Rapidity(),
        trk1.globalIndex(),
        trk2.globalIndex());
    }
  }

  void processPhisMC(
    MCEventCandidates::iterator const& collision,
    aod::McCollisions const&,
    MCTrackCandidates const& tracks,
    aod::McParticles const&)
  {
    if (!collision.has_mcCollision()) {
      return;
    }

    if (!isSelectedEvents(collision, false)) {
      return;
    }

    for (auto const& [trk1, trk2] :
         combinations(CombinationsFullIndexPolicy(tracks, tracks))) {

      // -----------------------------
      // Same cuts as phi1020analysis
      // -----------------------------

      if (trk2.index() <= trk1.index()) {
        continue;
      }

      if (trk1.collisionId() != collision.globalIndex() ||
          trk2.collisionId() != collision.globalIndex()) {
        continue;
      }

      if (!trackCut(trk1) || !trackCut(trk2)) {
        continue;
      }

      if (trk1.sign() * trk2.sign() > 0) {
        continue;
      }

      if (!selectionPID(trk1) || !selectionPID(trk2)) {
        continue;
      }

      LorentzVectorPtEtaPhiMass kplus;
      LorentzVectorPtEtaPhiMass kminus;
      LorentzVectorPtEtaPhiMass phi;

      // assign K+ and K-
      if (trk1.sign() > 0) {
        kplus = LorentzVectorPtEtaPhiMass(trk1.pt(), trk1.eta(), trk1.phi(),
                                          o2::constants::physics::MassKPlus);
        kminus = LorentzVectorPtEtaPhiMass(trk2.pt(), trk2.eta(), trk2.phi(),
                                           o2::constants::physics::MassKPlus);
      } else {
        kplus = LorentzVectorPtEtaPhiMass(trk2.pt(), trk2.eta(), trk2.phi(),
                                          o2::constants::physics::MassKPlus);
        kminus = LorentzVectorPtEtaPhiMass(trk1.pt(), trk1.eta(), trk1.phi(),
                                           o2::constants::physics::MassKPlus);
      }

      phi = kplus + kminus;

      float invMass = phi.M();

      // rapidity cut
      if (std::abs(phi.Rapidity()) > configResoDauTracks.cfgCutRapidity) {
        continue;
      }

      // optional mass window
      if (invMass < mMinMassPhi || invMass > mMaxMassPhi) {
        continue;
      }

      // Daughter QA (both legs are kaons for Phi)
      fillPhiKaonQA(trk1);
      fillPhiKaonQA(trk2);

      bool mcTruePhi = false;
      bool mcPhysicalPrimary = false;

      // -----------------------------
      // MC matching
      // -----------------------------

      if (trk1.has_mcParticle() && trk2.has_mcParticle()) {

        auto mc1 = trk1.mcParticle();
        auto mc2 = trk2.mcParticle();

        if (mc1.has_mothers() && mc2.has_mothers()) {

          for (auto const& mother1 : mc1.mothers_as<aod::McParticles>()) {
            for (auto const& mother2 : mc2.mothers_as<aod::McParticles>()) {

              if (mother1.globalIndex() != mother2.globalIndex()) {
                continue;
              }

              if (std::abs(mother1.pdgCode()) != 333) {
                continue;
              }

              mcTruePhi = true;
              mcPhysicalPrimary = mother1.isPhysicalPrimary();
            }
          }
        }
      }

      assocPhis(
        collision.globalIndex(),
        mcTruePhi,
        mcPhysicalPrimary,
        phi.Pt(),
        phi.Eta(),
        phi.Phi(),
        invMass,
        phi.Rapidity(),
        trk1.globalIndex(),
        trk2.globalIndex());
    }
  }

  // K*0 -> K+ pi- (and c.c.) reconstruction. Unlike phi (K+K-, symmetric daughter
  // mass), the two daughter tracks are different species, so both mass hypotheses
  // (trk1=kaon,trk2=pion) and (trk1=pion,trk2=kaon) are tried independently; either,
  // both, or neither may pass PID for a given unlike-sign pair.
  void processKstars(EventCandidates::iterator const& collision,
                     TrackCandidates const& tracks)
  {
    if (!isSelectedEvents(collision, false)) {
      return;
    }

    LorentzVectorPtEtaPhiMass kaon;
    LorentzVectorPtEtaPhiMass pion;
    LorentzVectorPtEtaPhiMass kstar;

    for (auto const& [trk1, trk2] :
         combinations(CombinationsFullIndexPolicy(tracks, tracks))) {

      if (trk2.index() <= trk1.index()) {
        continue;
      }

      // same collision
      if (trk1.collisionId() != collision.globalIndex() ||
          trk2.collisionId() != collision.globalIndex()) {
        continue;
      }

      // quality cuts
      if (!trackCut(trk1) || !trackCut(trk2)) {
        continue;
      }

      // unlike-sign only
      if (trk1.sign() * trk2.sign() > 0) {
        continue;
      }

      for (int hypothesis = 0; hypothesis < 2; hypothesis++) {
        auto const& kaonTrack = (hypothesis == 0) ? trk1 : trk2;
        auto const& pionTrack = (hypothesis == 0) ? trk2 : trk1;

        // kaon + pion PID
        bool passKaonPID = configPID.ispTdepPID ? ptDependentPidKaon(kaonTrack) : selectionPID(kaonTrack);
        bool passPionPID = configPID.ispTdepPID ? ptDependentPidPion(pionTrack) : selectionPIDPion(pionTrack);
        if (!passKaonPID || !passPionPID) {
          continue;
        }

        kaon = LorentzVectorPtEtaPhiMass(kaonTrack.pt(), kaonTrack.eta(), kaonTrack.phi(),
                                         o2::constants::physics::MassKPlus);
        pion = LorentzVectorPtEtaPhiMass(pionTrack.pt(), pionTrack.eta(), pionTrack.phi(),
                                         o2::constants::physics::MassPiPlus);

        kstar = kaon + pion;

        float invMass = kstar.M();

        // rapidity cut
        if (std::abs(kstar.Rapidity()) > configResoDauTracks.cfgCutRapidity) {
          continue;
        }

        // mass window
        if (invMass < mMinMassKstar || invMass > mMaxMassKstar) {
          continue;
        }

        // Daughter QA (kaon leg and pion leg filled separately)
        fillKstarKaonQA(kaonTrack);
        fillKstarPionQA(pionTrack);

        auto const& posTrack = (kaonTrack.sign() > 0) ? kaonTrack : pionTrack;
        auto const& negTrack = (kaonTrack.sign() > 0) ? pionTrack : kaonTrack;
        // Rapidity under whichever mass hypothesis (kaon or pion) was
        // actually assigned to that daughter for this hypothesis -- matches
        // posTrack/negTrack above, not a fixed species per pos/neg slot.

        assocKstars(
          collision.globalIndex(),
          false,
          false,
          kstar.Pt(),
          kstar.Eta(),
          kstar.Phi(),
          invMass,
          kstar.Rapidity(),
          posTrack.globalIndex(),
          negTrack.globalIndex());
      }
    }
  }

  void processKstarsMC(
    MCEventCandidates::iterator const& collision,
    aod::McCollisions const&,
    MCTrackCandidates const& tracks,
    aod::McParticles const&)
  {
    if (!collision.has_mcCollision()) {
      return;
    }

    if (!isSelectedEvents(collision, false)) {
      return;
    }

    for (auto const& [trk1, trk2] :
         combinations(CombinationsFullIndexPolicy(tracks, tracks))) {

      // -----------------------------
      // Same cuts as processKstars
      // -----------------------------

      if (trk2.index() <= trk1.index()) {
        continue;
      }

      if (trk1.collisionId() != collision.globalIndex() ||
          trk2.collisionId() != collision.globalIndex()) {
        continue;
      }

      if (!trackCut(trk1) || !trackCut(trk2)) {
        continue;
      }

      if (trk1.sign() * trk2.sign() > 0) {
        continue;
      }

      for (int hypothesis = 0; hypothesis < 2; hypothesis++) {
        auto const& kaonTrack = (hypothesis == 0) ? trk1 : trk2;
        auto const& pionTrack = (hypothesis == 0) ? trk2 : trk1;

        if (!selectionPID(kaonTrack) || !selectionPIDPion(pionTrack)) {
          continue;
        }

        LorentzVectorPtEtaPhiMass kaon = LorentzVectorPtEtaPhiMass(kaonTrack.pt(), kaonTrack.eta(), kaonTrack.phi(),
                                                                   o2::constants::physics::MassKPlus);
        LorentzVectorPtEtaPhiMass pion = LorentzVectorPtEtaPhiMass(pionTrack.pt(), pionTrack.eta(), pionTrack.phi(),
                                                                   o2::constants::physics::MassPiPlus);
        LorentzVectorPtEtaPhiMass kstar = kaon + pion;

        float invMass = kstar.M();

        // rapidity cut
        if (std::abs(kstar.Rapidity()) > configResoDauTracks.cfgCutRapidity) {
          continue;
        }

        // mass window
        if (invMass < mMinMassKstar || invMass > mMaxMassKstar) {
          continue;
        }

        // Daughter QA (kaon leg and pion leg filled separately)
        fillKstarKaonQA(kaonTrack);
        fillKstarPionQA(pionTrack);

        bool mcTrueKstar = false;
        bool mcPhysicalPrimary = false;

        // -----------------------------
        // MC matching: common mother with |pdgCode| == 313 (K*0/anti-K*0)
        // -----------------------------

        if (kaonTrack.has_mcParticle() && pionTrack.has_mcParticle()) {

          auto mcKaon = kaonTrack.mcParticle();
          auto mcPion = pionTrack.mcParticle();

          if (mcKaon.has_mothers() && mcPion.has_mothers()) {

            for (auto const& motherKaon : mcKaon.mothers_as<aod::McParticles>()) {
              for (auto const& motherPion : mcPion.mothers_as<aod::McParticles>()) {

                if (motherKaon.globalIndex() != motherPion.globalIndex()) {
                  continue;
                }

                if (std::abs(motherKaon.pdgCode()) != 313) {
                  continue;
                }

                mcTrueKstar = true;
                mcPhysicalPrimary = motherKaon.isPhysicalPrimary();
              }
            }
          }
        }

        auto const& posTrack = (kaonTrack.sign() > 0) ? kaonTrack : pionTrack;
        auto const& negTrack = (kaonTrack.sign() > 0) ? pionTrack : kaonTrack;

        assocKstars(
          collision.globalIndex(),
          mcTrueKstar,
          mcPhysicalPrimary,
          kstar.Pt(),
          kstar.Eta(),
          kstar.Phi(),
          invMass,
          kstar.Rapidity(),
          posTrack.globalIndex(),
          negTrack.globalIndex());
      }
    }
  }

  PROCESS_SWITCH(HResonanceCorrelationFilter, processTriggers, "Produce trigger tables", true);
  PROCESS_SWITCH(HResonanceCorrelationFilter, processTriggersMC, "Produce trigger tables for MC", false);
  PROCESS_SWITCH(HResonanceCorrelationFilter, processAssocPions, "Produce associated Pion tables", false);
  PROCESS_SWITCH(HResonanceCorrelationFilter, processAssocPionsMC, "Produce associated Pion tables for MC", false);
  PROCESS_SWITCH(HResonanceCorrelationFilter, processAssocKaons, "Produce associated Kaon tables", false);
  PROCESS_SWITCH(HResonanceCorrelationFilter, processAssocKaonsMC, "Produce associated Kaon tables for MC", false);
  PROCESS_SWITCH(HResonanceCorrelationFilter, processAssocHadrons, "Produce associated Hadron tables", true);
  PROCESS_SWITCH(HResonanceCorrelationFilter, processAssocHadronsMC, "Produce associated Hadron tables for MC", false);
  PROCESS_SWITCH(HResonanceCorrelationFilter, processPhis, "Produce associated phi tables", true);
  PROCESS_SWITCH(HResonanceCorrelationFilter, processPhisMC, "Produce associated phi tables for MC", false);
  PROCESS_SWITCH(HResonanceCorrelationFilter, processKstars, "Produce associated K*0 tables", true);
  PROCESS_SWITCH(HResonanceCorrelationFilter, processKstarsMC, "Produce associated K*0 tables for MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HResonanceCorrelationFilter>(cfgc)};
}
