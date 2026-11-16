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

/// \file chargedParticles.cxx
/// \brief Task for analysis of charged particle pt spectra vs multiplicity with 2d unfolding.
/// \author Mario Krüger <mario.kruger@cern.ch>

#include "PWGLF/DataModel/particleCompositionCorrectionTable.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/TriggerAliases.h"
#include "Common/Core/MetadataHelper.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CommonConstants/MathConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/runDataProcessing.h>

#include <TAxis.h>
#include <TH1.h>
#include <TH2.h>

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <string>
#include <unordered_set>
#include <vector>

using namespace o2;
using namespace o2::framework;
using aod::track::TrackSelectionFlags;

o2::common::core::MetadataHelper metadataInfo; // ao2d metadata (production, pass, anchor) for the qa

//--------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------
// Task declaration
//--------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------
struct ChargedParticles {
  // pcc repetitions: floor(weight) + 1 with probability frac(weight)
  // random number = hash of the particle -> same for all loops, reproducible, independent of processing order
  static uint64_t splitMix64(uint64_t x)
  {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
  }
  template <typename P>
  double getParticleUniform(const P& particle) const
  {
    uint64_t hash = seedOffset.value;
    for (const float val : {particle.px(), particle.py(), particle.pz(), particle.vx(), particle.vy(), particle.vz()}) {
      hash = splitMix64(hash ^ std::bit_cast<uint32_t>(val));
    }
    hash = splitMix64(hash ^ static_cast<uint32_t>(particle.pdgCode()));
    return std::ldexp(static_cast<double>(hash >> 11), -53); // uniform in [0, 1)
  }
  template <typename P>
  uint32_t getNRepetitons(float scalingFactor, const P& particle) const
  {
    if (!(scalingFactor > 0.f)) { // also catches nan
      return 0u;
    }
    auto nRepetitions = static_cast<uint32_t>(scalingFactor);
    const double rest = scalingFactor - nRepetitions;
    // hash only needed for fractional weights
    return nRepetitions + ((rest > 0. && getParticleUniform(particle) < rest) ? 1u : 0u);
  }

  HistogramRegistry histos;
  Service<o2::framework::O2DatabasePDG> pdg;

  Configurable<bool> isRun3{"isRun3", true, "is Run3 dataset"};
  Configurable<std::string> collisionSystem{"collisionSystem", "pp", "collision system (sets multiplicity axes and event selection default): pp, pPb, pO, OO, NeNe, XeXe, PbPb"};
  int maxMultMeas{100}; // from collisionSystem
  int maxMultTrue{100};
  Configurable<float> etaCut{"etaCut", 0.8f, "eta cut"};
  Configurable<float> ptMinCut{"ptMinCut", 0.15f, "pt min cut"};
  Configurable<float> ptMaxCut{"ptMaxCut", 10.f, "pt max cut"};
  Configurable<std::string> eventClass{"eventClass", "fiducial", "mc event class: fiducial (N_ch > 0), inelgt0, triggered (tvx), inel"};
  enum EventClass { kFiducial,
                    kINELGT0,
                    kTriggered,
                    kINEL };
  EventClass eventClassType{kFiducial};
  static constexpr int kNStagesMC = 8;
  static constexpr int kNLayersITS = 7;
  std::string metaLabel; // qa: ao2d metadata label
  int lastRunNumber{-1}; // qa: cached run label
  std::string lastRunLabel;
  Configurable<std::string> evSelBits{"evSelBits", "auto", "event selection bits on top of sel8 (names as in EventSelectionParams), comma separated; auto = default of the collision system (currently kNoSameBunchPileup for all)"};
  Configurable<float> maxVtxZ{"maxVtxZ", 10.f, "max |z| of the reco and true vertex"};
  uint64_t evSelMask{0};
  std::vector<int> evSelBitList; // configured bits in the given order (cut-flow)
  std::vector<int> trackCutBits; // single cuts of the track selection mask (track qa)
  bool isNominal{true};          // qa only for the nominal setting
  Configurable<uint32_t> seedOffset{"seedOffset", 0, "offset for the (deterministic) random numbers of the particle composition correction"};

  enum : uint32_t {
    kSystNominal = 100,
    kSystDownChi2PerClusterITS,
    kSystUpChi2PerClusterITS,
    kSystDownChi2PerClusterTPC,
    kSystUpChi2PerClusterTPC,
    kSystDownTPCCrossedRowsOverNCls,
    kSystUpTPCCrossedRowsOverNCls,
    kSystDownDCAxy = 111,
    kSystUpDCAxy,
    kSystDownDCAz,
    kSystUpDCAz,
    kSystITSHits, // only relevant for converted data
    kSystDownTPCCrossedRows,
    kSystUpTPCCrossedRows,
    kSystDownPCC = 120,
    kSystUpPCC
  };
  Configurable<uint32_t> systMode{"systMode", kSystNominal, "variation for systematic uncertainties"};
  uint16_t trackSelMask{TrackSelectionFlags::kGlobalTrackWoPtEta}; // track selection bitmask (without cut that is being varied)
  uint16_t cutVarFlag{0};
  TrackSelection trackSel;
  TrackSelection::TrackCuts trackSelFlag;

  // helper struct to store transient properties
  struct VarContainer {
    uint32_t multMeas{0u};
    uint32_t multTrue{0u};
    bool isAcceptedEvent{false};
    bool isAcceptedEventMC{false};
    bool isGoodZPosMC{false};  // true vertex within |z| cut
    bool isTriggeredMC{false}; // tvx in the bc of the mc collision
    bool isINELGT0EventMC{false};
    bool isChargedPrimary{false};
    uint32_t nRepetitions{1u};
    int64_t mcCollisionIndex{-1}; // current mc collision
  };
  VarContainer vars;

  void init(InitContext const&);

  template <typename P>
  bool initParticle(const P& particle);

  template <typename T>
  bool initTrack(const T& track);

  template <bool IS_MC, typename P, typename T>
  uint32_t getTrackRepetitions(const T& track, bool& isValidParticle);

  template <bool IS_MC, typename P = aod::McParticles, typename C, typename T>
  void initEvent(const C& collision, const T& tracks);

  template <typename C, typename P>
  void initEventMC(const C& collision, const P& particles);

  template <bool IS_MC, typename P = aod::McParticles, typename T>
  void processMeas(const T& tracks, bool isRepresentative = true);

  template <typename P>
  void processTrue(const P& particles);

  template <typename T>
  void fillTrackQA(const T& track);

  template <typename T>
  void fillTrackCutQA(const T& track);

  template <bool IS_MC, typename C>
  void fillEvSelQA(const C& collision, aod::FT0s const& ft0s, bool isRepresentative = false);

  using CollisionTableData = soa::Join<aod::Collisions, aod::EvSels>;
  using TrackTableData = soa::Join<aod::FullTracks, aod::TracksCov, aod::TracksDCA, aod::TrackSelection>;
  void processData(CollisionTableData::iterator const& collision, TrackTableData const& tracks, aod::FT0s const& ft0s, aod::BCs const&);
  PROCESS_SWITCH(ChargedParticles, processData, "process data", false);

  using CollisionTableMCTrue = aod::McCollisions;
  using CollisionTableMC = soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions, aod::EvSels>>;
  using TrackTableMC = soa::Join<aod::FullTracks, aod::TracksCov, aod::TracksDCA, aod::TrackSelection, aod::McTrackLabels>;
  using ParticleTableMC = soa::Join<aod::McParticles, aod::ParticleCompositionCorrection>;
  using BCsWithSel = soa::Join<aod::BCs, aod::BcSels>; // trigger info of the mc collision's bc
  Preslice<TrackTableMC> perCollision = aod::track::collisionId;
  template <typename P>
  void processMCImpl(CollisionTableMCTrue::iterator const& mcCollision, TrackTableMC const& tracks, CollisionTableMC const& collisions, P const& particles, aod::FT0s const& ft0s);

  void processMC(CollisionTableMCTrue::iterator const& mcCollision, TrackTableMC const& tracks, CollisionTableMC const& collisions, ParticleTableMC const& particles, BCsWithSel const&, aod::FT0s const& ft0s);
  PROCESS_SWITCH(ChargedParticles, processMC, "process mc", true);

  // mc without pcc: no pcc table subscribed, pcc producer not needed
  void processMCnoPCC(CollisionTableMCTrue::iterator const& mcCollision, TrackTableMC const& tracks, CollisionTableMC const& collisions, aod::McParticles const& particles, BCsWithSel const&, aod::FT0s const& ft0s);
  PROCESS_SWITCH(ChargedParticles, processMCnoPCC, "process mc without particle composition correction", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  metadataInfo.initMetadata(cfgc);
  return WorkflowSpec{adaptAnalysisTask<ChargedParticles>(cfgc)};
}

//--------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------
// Task implementation
//--------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------

//**************************************************************************************************
/**
 * Initialise the task and add histograms.
 */
//**************************************************************************************************
void ChargedParticles::init(InitContext const&)
{
  // collision system: multiplicity axes (run 2 values; light ions estimated) and event selection default
  // todo: measure the event selection default per system (data vs mc rejection, qa/events/multMeasBits)
  std::string evSelBitsDefault = "kNoSameBunchPileup";
  if (collisionSystem.value == "pp") {
    maxMultMeas = 100;
    maxMultTrue = 100;
  } else if (collisionSystem.value == "pPb" || collisionSystem.value == "pO") {
    maxMultMeas = 200;
    maxMultTrue = 200;
  } else if (collisionSystem.value == "OO" || collisionSystem.value == "NeNe") {
    maxMultMeas = 400;
    maxMultTrue = 600;
  } else if (collisionSystem.value == "XeXe" || collisionSystem.value == "PbPb") {
    maxMultMeas = 2500;
    maxMultTrue = 3800;
  } else {
    LOGP(fatal, "unknown collisionSystem '{}' (pp, pPb, pO, OO, NeNe, XeXe, PbPb)", collisionSystem.value);
  }

  if (eventClass.value == "fiducial") {
    eventClassType = kFiducial;
  } else if (eventClass.value == "inelgt0") {
    eventClassType = kINELGT0;
  } else if (eventClass.value == "triggered") {
    eventClassType = kTriggered;
  } else if (eventClass.value == "inel") {
    eventClassType = kINEL;
  } else {
    LOGP(fatal, "unknown eventClass '{}' (fiducial, inelgt0, triggered, inel)", eventClass.value);
  }

  // event selection bits required on top of sel8
  const std::string bits = (evSelBits.value == "auto") ? evSelBitsDefault : evSelBits.value;
  for (size_t start = 0; start < bits.size();) {
    size_t end = bits.find(',', start);
    if (end == std::string::npos) {
      end = bits.size();
    }
    std::string name = bits.substr(start, end - start);
    name.erase(0, name.find_first_not_of(' '));
    name.erase(name.find_last_not_of(' ') + 1);
    start = end + 1;
    if (name.empty()) {
      continue;
    }
    int bit = -1;
    for (int i = 0; i < aod::evsel::kNsel; ++i) {
      if (name == aod::evsel::selectionLabels[i]) {
        bit = i;
      }
    }
    if (bit < 0) {
      LOGP(fatal, "unknown event selection bit '{}' in evSelBits", name);
    }
    evSelMask |= (1ULL << bit);
    evSelBitList.push_back(bit);
  }

  std::vector<double> ptBinEdges = {0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75,
                                    0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9,
                                    2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0, 5.5,
                                    6.0, 6.5, 7.0, 8.0, 9.0, 10.0};

  const AxisSpec ptMeasAxis{ptBinEdges, "#it{p}^{ meas}_{T} (GeV/#it{c})", "pt_meas"};

  const int nBinsMultMeas = maxMultMeas + 1;
  const AxisSpec multMeasAxis = {nBinsMultMeas, -0.5, nBinsMultMeas - 0.5, "#it{N}^{ meas}_{ch}", "mult_meas"};

  histos.add("multDist_evt_meas", "", kTH1D, {multMeasAxis});               // measured event distribution (contains contamination from events not in specified class or with wrong vertex position)
  histos.add("multPtSpec_trk_meas", "", kTH2D, {multMeasAxis, ptMeasAxis}); // measured tracks (contains contamination from secondary particles, particles smeared into acceptance and tracks originating from background events as defined above )

  // event selection qa: each evsel bit separately on top of sel8 (data vs mc rejection vs multiplicity)
  isNominal = (systMode == kSystNominal);
  const AxisSpec bitAxis{aod::evsel::kNsel + 2, -0.5, static_cast<double>(aod::evsel::kNsel) + 1.5, "", "evsel_bit"};
  const int nStages = static_cast<int>(evSelBitList.size()) + 3; // all, sel8, bits (cumulative), |z|
  if (isNominal) {
    histos.add("QA/events/multMeasBits", "", kTH2D, {multMeasAxis, bitAxis});
    histos.add("QA/events/dzFT0PV", "", kTH2D, {multMeasAxis, {100, -5., 5., "#it{z}_{FT0} #minus #it{z}_{PV} (cm)"}});
    histos.add("QA/events/zVtxMeas", "", kTH1D, {{61, -30.5, 30.5, "#it{V}^{ meas}_{z} (cm)"}});
    histos.add("QA/events/eventSelection", "", kTH1D, {{nStages, -0.5, nStages - 0.5}});
    histos.add("QA/events/runStatistics", "", kTH1D, {{512, 0., 512.}});
    histos.add("QA/events/metaInfo", "", kTH1D, {{4, 0., 4.}});
    // track qa (selected tracks in selected events), run 2 names and binning (7 its layers, run 3 dca cut ranges)
    histos.add("QA/tracks/signed1Pt", "", kTH1D, {{200, -8., 8., "q/#it{p}_{T}"}});
    histos.add("QA/tracks/eta", "", kTH1D, {{180, -0.9, 0.9, "#eta"}});
    histos.add("QA/tracks/phi", "", kTH1D, {{180, 0., o2::constants::math::TwoPI, "#phi (rad)"}});
    histos.add("QA/tracks/itsNCls", "", kTH1D, {{8, -0.5, 7.5, "# clusters ITS"}});
    histos.add("QA/tracks/itsHits", "", kTH1D, {{7, -0.5, 6.5, "layer ITS"}});
    histos.add("QA/tracks/itsChi2NCl", "", kTH1D, {{100, 0., 40., "chi2 per cluster ITS"}});
    histos.add("QA/tracks/tpcNClsFindable", "", kTH1D, {{165, -0.5, 164.5, "# findable clusters TPC"}});
    histos.add("QA/tracks/tpcNClsFound", "", kTH1D, {{165, -0.5, 164.5, "# clusters TPC"}});
    histos.add("QA/tracks/tpcCrossedRows", "", kTH1D, {{165, -0.5, 164.5, "# crossed rows TPC"}});
    histos.add("QA/tracks/tpcCrossedRowsOverFindableCls", "", kTH1D, {{90, 0.7, 1.6, "crossed rows / findable clusters TPC"}});
    histos.add("QA/tracks/tpcChi2NCl", "", kTH1D, {{100, 0., 10., "chi2 / cluster TPC"}});
    histos.add("QA/tracks/sigmaPt", "", kTH2D, {ptMeasAxis, {100, 0., 0.3, "#sigma(#it{p}^{ meas}_{T}) / #it{p}^{ meas}_{T}"}});
    histos.add("QA/tracks/dcaXY", "", kTH2D, {ptMeasAxis, {500, -0.5, 0.5, "dca XY (cm)"}});
    histos.add("QA/tracks/dcaZ", "", kTH2D, {ptMeasAxis, {500, -2.5, 2.5, "dca Z (cm)"}});
    // track cut qa: n-1 and cumulative per single cut of the selection mask, n-1 dca distributions (no dca cut)
    // (nominal only, so the mask does not contain a varied cut yet)
    for (int i = 0; i < static_cast<int>(TrackSelection::TrackCuts::kNCuts); ++i) {
      if ((trackSelMask & (1 << i)) != 0) {
        trackCutBits.push_back(i);
      }
    }
    const int nCutBins = 2 * static_cast<int>(trackCutBits.size()) + 3;
    histos.add("QA/tracks/trackCuts", "", kTH2D, {ptMeasAxis, {nCutBins, -0.5, nCutBins - 0.5}});
    histos.add("QA/tracks/dcaXYnm1", "", kTH2D, {ptMeasAxis, {800, -1., 1., "dca XY (cm)"}});
    histos.add("QA/tracks/dcaZnm1", "", kTH2D, {ptMeasAxis, {600, -3., 3., "dca Z (cm)"}});
  }

  if (doprocessMC || doprocessMCnoPCC) {

    const AxisSpec ptTrueAxis{ptBinEdges, "#it{p}_{T} (GeV/c)", "pt_true"};

    const int nBinsMultTrue = maxMultTrue + 1;
    const AxisSpec multTrueAxis = {nBinsMultTrue, -0.5, nBinsMultTrue - 0.5, "#it{N}_{ch}", "mult_true"};

    histos.add("collision_ambiguity", "", kTH1D, {{6, -0.5, 5.5, "reco collisions per true collision"}}); // log the number of collisions that were reconstructed for a MC collision
    histos.add("track_ambiguity", "", kTH1D, {{6, 0.5, 6.5, "reco tracks per true particle"}});           // log the number of tracks that were reconstructed for a MC particle

    histos.add("multDist_evt_gen", "", kTH1D, {multTrueAxis}); // generated event distribution  (from events within specified class and with proper vertex position)
    if (isNominal) {
      histos.add("QA/events/multTrueBits", "", kTH2D, {multTrueAxis, bitAxis});
      histos.add("QA/events/zVtxGen", "", kTH1D, {{61, -30.5, 30.5, "#it{V}_{z} (cm)"}});
      histos.add("QA/events/eventSelectionMC", "", kTH1D, {{kNStagesMC, -0.5, kNStagesMC - 0.5}});
      histos.add("QA/tracks/deltaPt", "", kTH2D, {ptMeasAxis, {100, 0., 0.3, "#Delta(#it{p}_{T}) / #it{p}^{ meas}_{T}"}});
    }
    histos.add("multDist_evt_gen_trig", "", kTH1D, {multTrueAxis}); // generated event distribution (from events within specified class and with proper vertex position) that in addition fulfils the trigger condition [to disentangle trigger eff from reco eff ]

    histos.add("multCorrel_evt", "", kTHnSparseD, {multMeasAxis, multTrueAxis});  // multiplicity correlation of measured events (excluding background events)
    histos.add("multCorrel_prim", "", kTHnSparseD, {multMeasAxis, multTrueAxis}); // multiplicity correlation of measured primary charged particles (excluding particles from background events)
    histos.add("ptCorrel_prim", "", kTH2D, {ptMeasAxis, ptTrueAxis});             // pT correlation of measured primary charged particles  (excluding particles from background events)

    histos.add("multPtSpec_prim_gen", "", kTH2D, {multTrueAxis, ptTrueAxis});         // generated primary charged particles as function of true properties (from events within specified class and with proper vertex position)
    histos.add("multPtSpec_prim_gen_evtloss", "", kTH2D, {multTrueAxis, ptTrueAxis}); // generated primary charged particles of events that did not pass the event selection as function of multiplicity and pt
    histos.add("multPtSpec_prim_gen_notrig", "", kTH2D, {multTrueAxis, ptTrueAxis});  // generated primary charged particles of untriggered events (tvx) with good z, independent of the event class
    histos.add("multPtSpec_prim_meas", "", kTH2D, {multTrueAxis, ptTrueAxis});        // measured primary charged particles as function of true properties (no contamination from background events)

    histos.add("multPtSpec_trk_prim_meas", "", kTH2D, {multMeasAxis, ptMeasAxis});    // tracks from measured primaries (no contamination from secondaries, particles smeared into acceptance or background events)
    histos.add("multPtSpec_trk_sec_meas", "", kTH2D, {multMeasAxis, ptMeasAxis});     // tracks from measured secondaries (no contamination from particles smeared into acceptance or background events)  [for QA to disentangle secondaries from other contamination]
    histos.add("multPtSpec_trk_meas_evtcont", "", kTH2D, {multMeasAxis, ptMeasAxis}); // tracks from events that are measured, but do not belong to the desired class of events
    histos.add("multPtSpec_trk_inter", "", kTH2D, {multTrueAxis, ptMeasAxis});
  }

  auto setBitLabels = [](TAxis* axis) {
    for (int i = 0; i < aod::evsel::kNsel; ++i) {
      axis->SetBinLabel(i + 1, aod::evsel::selectionLabels[i]);
    }
    axis->SetBinLabel(aod::evsel::kNsel + 1, "sel8");
    axis->SetBinLabel(aod::evsel::kNsel + 2, "selected");
  };
  if (isNominal) {
    setBitLabels(histos.get<TH2>(HIST("QA/events/multMeasBits"))->GetYaxis());
    auto* axis = histos.get<TH1>(HIST("QA/events/eventSelection"))->GetXaxis();
    axis->SetBinLabel(1, "all");
    axis->SetBinLabel(2, "sel8");
    for (size_t i = 0; i < evSelBitList.size(); ++i) {
      axis->SetBinLabel(i + 3, aod::evsel::selectionLabels[evSelBitList[i]]);
    }
    axis->SetBinLabel(nStages, "|z| < maxVtxZ");
    auto* axisCuts = histos.get<TH2>(HIST("QA/tracks/trackCuts"))->GetYaxis();
    const int nCuts = trackCutBits.size();
    axisCuts->SetBinLabel(1, "candidates");
    for (int k = 0; k < nCuts; ++k) {
      axisCuts->SetBinLabel(k + 2, ("w/o " + TrackSelection::mCutNames[trackCutBits[k]]).c_str());
      axisCuts->SetBinLabel(nCuts + k + 4, ("after " + TrackSelection::mCutNames[trackCutBits[k]]).c_str());
    }
    axisCuts->SetBinLabel(nCuts + 2, "w/o dca");
    axisCuts->SetBinLabel(nCuts + 3, "all cuts");
    auto* axisHits = histos.get<TH1>(HIST("QA/tracks/itsHits"))->GetXaxis();
    for (int i = 0; i < kNLayersITS; ++i) {
      axisHits->SetBinLabel(i + 1, ("L" + std::to_string(i)).c_str());
    }
    metaLabel = metadataInfo.isInitialized() ? metadataInfo.makeMetadataLabel() + ", O2 " + metadataInfo.getO2Version() : std::string("no ao2d metadata");
    if (doprocessMC || doprocessMCnoPCC) {
      metaLabel += ", " + eventClass.value;
      setBitLabels(histos.get<TH2>(HIST("QA/events/multTrueBits"))->GetYaxis());
      auto* axisMC = histos.get<TH1>(HIST("QA/events/eventSelectionMC"))->GetXaxis();
      int bin = 1;
      for (const auto* label : {"all", "|z_{true}| < maxVtxZ", "event class", "tvx", "reco collision", "sel8", "evSelBits", "|z| < maxVtxZ"}) {
        axisMC->SetBinLabel(bin++, label);
      }
    }
  }

  trackSel = getGlobalTrackSelection();
  if (systMode == kSystDownChi2PerClusterITS) {
    trackSel.SetMaxChi2PerClusterITS(25.);
    cutVarFlag = TrackSelectionFlags::kITSChi2NDF;
    trackSelFlag = TrackSelection::TrackCuts::kITSChi2NDF;
  } else if (systMode == kSystUpChi2PerClusterITS) {
    trackSel.SetMaxChi2PerClusterITS(49.);
    cutVarFlag = TrackSelectionFlags::kITSChi2NDF;
    trackSelFlag = TrackSelection::TrackCuts::kITSChi2NDF;
  } else if (systMode == kSystDownChi2PerClusterTPC) {
    trackSel.SetMaxChi2PerClusterTPC(3.0);
    cutVarFlag = TrackSelectionFlags::kTPCChi2NDF;
    trackSelFlag = TrackSelection::TrackCuts::kTPCChi2NDF;
  } else if (systMode == kSystUpChi2PerClusterTPC) {
    trackSel.SetMaxChi2PerClusterTPC(5.0);
    cutVarFlag = TrackSelectionFlags::kTPCChi2NDF;
    trackSelFlag = TrackSelection::TrackCuts::kTPCChi2NDF;
  } else if (systMode == kSystDownTPCCrossedRowsOverNCls) {
    trackSel.SetMinNCrossedRowsOverFindableClustersTPC(0.7);
    cutVarFlag = TrackSelectionFlags::kTPCCrossedRowsOverNCls;
    trackSelFlag = TrackSelection::TrackCuts::kTPCCrossedRowsOverNCls;
  } else if (systMode == kSystUpTPCCrossedRowsOverNCls) {
    trackSel.SetMinNCrossedRowsOverFindableClustersTPC(0.9);
    cutVarFlag = TrackSelectionFlags::kTPCCrossedRowsOverNCls;
    trackSelFlag = TrackSelection::TrackCuts::kTPCCrossedRowsOverNCls;
  } else if (systMode == kSystDownDCAxy) {
    trackSel.SetMaxDcaXYPtDep([](float pt) { return 4. / 7. * (0.0105f + 0.0350f / std::pow(pt, 1.1f)); });
    cutVarFlag = TrackSelectionFlags::kDCAxy;
    trackSelFlag = TrackSelection::TrackCuts::kDCAxy;
  } else if (systMode == kSystUpDCAxy) {
    trackSel.SetMaxDcaXYPtDep([](float pt) { return 10. / 7. * (0.0105f + 0.0350f / std::pow(pt, 1.1f)); });
    cutVarFlag = TrackSelectionFlags::kDCAxy;
    trackSelFlag = TrackSelection::TrackCuts::kDCAxy;
  } else if (systMode == kSystDownDCAz) {
    trackSel.SetMaxDcaZ(1.0);
    cutVarFlag = TrackSelectionFlags::kDCAz;
    trackSelFlag = TrackSelection::TrackCuts::kDCAz;
  } else if (systMode == kSystUpDCAz) {
    trackSel.SetMaxDcaZ(5.0);
    cutVarFlag = TrackSelectionFlags::kDCAz;
    trackSelFlag = TrackSelection::TrackCuts::kDCAz;
  } else if (systMode == kSystITSHits) {
    trackSel.ResetITSRequirements();
    cutVarFlag = TrackSelectionFlags::kITSHits;
    trackSelFlag = TrackSelection::TrackCuts::kITSHits;
  } else if (systMode == kSystDownTPCCrossedRows) {
    trackSel.SetMinNCrossedRowsTPC(60);
    cutVarFlag = TrackSelectionFlags::kTPCCrossedRows;
    trackSelFlag = TrackSelection::TrackCuts::kTPCCrossedRows;
  } else if (systMode == kSystUpTPCCrossedRows) {
    trackSel.SetMinNCrossedRowsTPC(80);
    cutVarFlag = TrackSelectionFlags::kTPCCrossedRows;
    trackSelFlag = TrackSelection::TrackCuts::kTPCCrossedRows;
  }
  trackSelMask &= (~cutVarFlag);
}

//**************************************************************************************************
/**
 * Entrypoint to processes data.
 */
//**************************************************************************************************
void ChargedParticles::processData(CollisionTableData::iterator const& collision, TrackTableData const& tracks, aod::FT0s const& ft0s, aod::BCs const&)
{
  initEvent<false>(collision, tracks);
  fillEvSelQA<false>(collision, ft0s);
  processMeas<false>(tracks);
}

//**************************************************************************************************
/**
 * Entrypoint to processes mc.
 */
//**************************************************************************************************
void ChargedParticles::processMC(CollisionTableMCTrue::iterator const& mcCollision, TrackTableMC const& tracks, CollisionTableMC const& collisions, ParticleTableMC const& particles, BCsWithSel const&, aod::FT0s const& ft0s)
{
  processMCImpl(mcCollision, tracks, collisions, particles, ft0s);
}

void ChargedParticles::processMCnoPCC(CollisionTableMCTrue::iterator const& mcCollision, TrackTableMC const& tracks, CollisionTableMC const& collisions, aod::McParticles const& particles, BCsWithSel const&, aod::FT0s const& ft0s)
{
  processMCImpl(mcCollision, tracks, collisions, particles, ft0s);
}

template <typename P>
void ChargedParticles::processMCImpl(CollisionTableMCTrue::iterator const& mcCollision, TrackTableMC const& tracks, CollisionTableMC const& collisions, P const& particles, aod::FT0s const& ft0s)
{
  histos.fill(HIST("collision_ambiguity"), collisions.size());

  initEventMC(mcCollision, particles);

  // several reco collisions per mc collision (rare, mostly split-off vertices): the one with most contributors
  // represents the mc collision, the others count as event contamination (as in data)
  int64_t representativeIndex = -1;
  uint16_t maxContrib = 0;
  for (const auto& collision : collisions) {
    if (representativeIndex < 0 || collision.numContrib() > maxContrib) {
      representativeIndex = collision.globalIndex();
      maxContrib = collision.numContrib();
    }
  }
  bool isAcceptedRepresentative = false;
  bool isSel8Representative = false;
  bool hasBitsRepresentative = false;
  for (const auto& collision : collisions) {
    auto curTracks = tracks.sliceBy(perCollision, collision.globalIndex());
    initEvent<true, P>(collision, curTracks);
    const bool isRepresentative = (collision.globalIndex() == representativeIndex);
    fillEvSelQA<true>(collision, ft0s, isRepresentative);
    processMeas<true, P>(curTracks, isRepresentative);
    if (isRepresentative) {
      isAcceptedRepresentative = vars.isAcceptedEvent;
      isSel8Representative = isRun3 && collision.sel8();
      hasBitsRepresentative = (collision.selection_raw() & evSelMask) == evSelMask;
    }
  }
  vars.isAcceptedEvent = isAcceptedRepresentative; // event loss refers to the representative; processTrue must come after this

  if (isNominal) {
    histos.fill(HIST("QA/events/zVtxGen"), mcCollision.posZ());
    // cumulative mc event flow
    const std::array<bool, kNStagesMC> stages = {true, vars.isGoodZPosMC, vars.isAcceptedEventMC, vars.isTriggeredMC, collisions.size() > 0,
                                                 isSel8Representative, hasBitsRepresentative, isAcceptedRepresentative};
    for (int i = 0; i < kNStagesMC && stages[i]; ++i) {
      histos.fill(HIST("QA/events/eventSelectionMC"), i);
    }
  }
  processTrue(particles);
}

//**************************************************************************************************
/**
 * Check if particle is good.
 */
//**************************************************************************************************
template <typename P>
bool ChargedParticles::initParticle(const P& particle)
{
  vars.isChargedPrimary = false;
  vars.nRepetitions = 1u;
  auto pdgParticle = pdg->GetParticle(particle.pdgCode());
  if (!pdgParticle || pdgParticle->Charge() == 0.) {
    return false;
  }
  vars.isChargedPrimary = particle.isPhysicalPrimary();

  // repetitions before the kinematic cuts (needed for inel>0)
  float pccWeight = 1.f;
  if constexpr (requires { particle.pccWeight(); }) {
    pccWeight = particle.pccWeight();
    if (systMode == kSystDownPCC) {
      pccWeight = particle.pccWeightSysDown();
    } else if (systMode == kSystUpPCC) {
      pccWeight = particle.pccWeightSysUp();
    }
  }
  vars.nRepetitions = getNRepetitons(pccWeight, particle);

  if (std::abs(particle.eta()) >= etaCut) {
    return false;
  }
  if (particle.pt() <= ptMinCut || particle.pt() >= ptMaxCut) {
    return false;
  }
  return true;
}

//**************************************************************************************************
/**
 * how often a track counts: pcc repetitions of its particle (1 for neutral/unknown particles), 1 in data and for tracks without mc particle.
 * valid particle: passes initParticle and belongs to the current mc collision.
 */
//**************************************************************************************************
template <bool IS_MC, typename P, typename T>
uint32_t ChargedParticles::getTrackRepetitions(const T& track, bool& isValidParticle)
{
  isValidParticle = false;
  if constexpr (IS_MC) {
    if (track.has_mcParticle()) {
      const auto& particle = track.template mcParticle_as<P>();
      const bool isInAcceptance = initParticle(particle);
      // particle outside the true acceptance or from another mc collision (pileup) -> contamination, but keeps its repetitions
      isValidParticle = isInAcceptance && (particle.mcCollisionId() == vars.mcCollisionIndex);
      return vars.nRepetitions;
    }
  }
  return 1u;
}

//**************************************************************************************************
/**
 * Check if track is good.
 */
//**************************************************************************************************
template <typename T>
bool ChargedParticles::initTrack(const T& track)
{
  if (std::abs(track.eta()) >= etaCut) {
    return false;
  }
  if (track.pt() <= ptMinCut || track.pt() >= ptMaxCut) {
    return false;
  }
  if (!TrackSelectionFlags::checkFlag(track.trackCutFlag(), trackSelMask)) {
    return false;
  }
  // for systematic variation of standard selections, check if the varied cut is passed
  if (cutVarFlag && !trackSel.IsSelected(track, trackSelFlag)) {
    return false;
  }
  return true;
}

//**************************************************************************************************
/**
 * Check if event is good.
 */
//**************************************************************************************************
template <bool IS_MC, typename P, typename C, typename T>
void ChargedParticles::initEvent(const C& collision, const T& tracks)
{
  vars.multMeas = 0;
  bool isValidParticle{false};
  for (const auto& track : tracks) {
    if (initTrack(track)) {
      // same definition in data and mc: every selected track counts
      vars.multMeas += getTrackRepetitions<IS_MC, P>(track, isValidParticle);
    }
  }

  vars.isAcceptedEvent = false;
  if (std::abs(collision.posZ()) < maxVtxZ) {
    if (isRun3) {
      if (collision.sel8() && (collision.selection_raw() & evSelMask) == evSelMask) {
        vars.isAcceptedEvent = true;
      }
    } else {
      if (collision.sel7() && ((doprocessMC || doprocessMCnoPCC) ? true : collision.alias_bit(kINT7))) {
        vars.isAcceptedEvent = true;
      }
    }
  }
}

//**************************************************************************************************
/**
 * Check if MC event is good.
 */
//**************************************************************************************************
template <typename C, typename P>
void ChargedParticles::initEventMC(const C& collision, const P& particles)
{
  vars.isINELGT0EventMC = false;
  vars.mcCollisionIndex = collision.globalIndex();
  vars.multTrue = 0;
  for (const auto& particle : particles) {
    const bool isInAcceptance = initParticle(particle);
    if (!vars.isChargedPrimary || vars.nRepetitions == 0) {
      continue;
    }
    if (std::abs(particle.eta()) < 1.f) {
      vars.isINELGT0EventMC = true; // inel>0: charged primary in |eta| < 1
    }
    if (isInAcceptance) {
      vars.multTrue += vars.nRepetitions;
    }
  }
  vars.isGoodZPosMC = std::abs(collision.posZ()) < maxVtxZ;
  vars.isTriggeredMC = collision.template bc_as<BCsWithSel>().selection_bit(aod::evsel::kIsTriggerTVX);
  bool isGoodEventClass = vars.multTrue > 0;
  if (eventClassType == kINELGT0) {
    isGoodEventClass = vars.isINELGT0EventMC;
  } else if (eventClassType == kTriggered) {
    isGoodEventClass = vars.isTriggeredMC;
  } else if (eventClassType == kINEL) {
    isGoodEventClass = true; // all generated events
  }
  vars.isAcceptedEventMC = isGoodEventClass && vars.isGoodZPosMC;
}

//**************************************************************************************************
/**
 * Track cut qa: candidates = tracks in selected events within the eta/pt window (nominal setting only).
 * n-1 per single cut, cumulative cut-flow in bit order, dca distributions without the dca cuts.
 */
//**************************************************************************************************
template <typename T>
void ChargedParticles::fillTrackCutQA(const T& track)
{
  if (std::abs(track.eta()) >= etaCut || track.pt() <= ptMinCut || track.pt() >= ptMaxCut) {
    return;
  }
  const auto flags = track.trackCutFlag();
  const int nCuts = trackCutBits.size();
  const auto passes = [&](int bit) { return (flags & (1 << bit)) != 0; };
  int nFailed = 0;
  int failedBit = -1;
  for (const int bit : trackCutBits) {
    if (!passes(bit)) {
      ++nFailed;
      failedBit = bit;
    }
  }
  const float pt = track.pt();
  histos.fill(HIST("QA/tracks/trackCuts"), pt, 0);
  for (int k = 0; k < nCuts; ++k) {
    if (nFailed == 0 || (nFailed == 1 && failedBit == trackCutBits[k])) {
      histos.fill(HIST("QA/tracks/trackCuts"), pt, k + 1); // passes all other cuts
    }
  }
  const auto dcaMask = (1 << static_cast<int>(TrackSelection::TrackCuts::kDCAxy)) | (1 << static_cast<int>(TrackSelection::TrackCuts::kDCAz));
  bool passesWoDca = true;
  for (const int bit : trackCutBits) {
    if (!((1 << bit) & dcaMask) && !passes(bit)) {
      passesWoDca = false;
    }
  }
  if (passesWoDca) {
    histos.fill(HIST("QA/tracks/trackCuts"), pt, nCuts + 1);
    histos.fill(HIST("QA/tracks/dcaXYnm1"), pt, track.dcaXY());
    histos.fill(HIST("QA/tracks/dcaZnm1"), pt, track.dcaZ());
  }
  if (nFailed == 0) {
    histos.fill(HIST("QA/tracks/trackCuts"), pt, nCuts + 2);
  }
  for (int k = 0; k < nCuts; ++k) { // cumulative
    if (!passes(trackCutBits[k])) {
      break;
    }
    histos.fill(HIST("QA/tracks/trackCuts"), pt, nCuts + 3 + k);
  }
}

//**************************************************************************************************
/**
 * Track qa (selected tracks in selected events, nominal setting only).
 */
//**************************************************************************************************
template <typename T>
void ChargedParticles::fillTrackQA(const T& track)
{
  histos.fill(HIST("QA/tracks/signed1Pt"), track.signed1Pt());
  histos.fill(HIST("QA/tracks/eta"), track.eta());
  histos.fill(HIST("QA/tracks/phi"), track.phi());
  histos.fill(HIST("QA/tracks/itsNCls"), track.itsNCls());
  for (int i = 0; i < kNLayersITS; ++i) {
    if (track.itsClusterMap() & (1 << i)) {
      histos.fill(HIST("QA/tracks/itsHits"), i);
    }
  }
  histos.fill(HIST("QA/tracks/itsChi2NCl"), track.itsChi2NCl());
  histos.fill(HIST("QA/tracks/tpcNClsFindable"), track.tpcNClsFindable());
  histos.fill(HIST("QA/tracks/tpcNClsFound"), track.tpcNClsFound());
  histos.fill(HIST("QA/tracks/tpcCrossedRows"), track.tpcNClsCrossedRows());
  histos.fill(HIST("QA/tracks/tpcCrossedRowsOverFindableCls"), track.tpcCrossedRowsOverFindableCls());
  histos.fill(HIST("QA/tracks/tpcChi2NCl"), track.tpcChi2NCl());
  histos.fill(HIST("QA/tracks/sigmaPt"), track.pt(), track.sigma1Pt() * track.pt());
  histos.fill(HIST("QA/tracks/dcaXY"), track.pt(), track.dcaXY());
  histos.fill(HIST("QA/tracks/dcaZ"), track.pt(), track.dcaZ());
}

//**************************************************************************************************
/**
 * Event qa: cut-flow, z vertex, run statistics, metadata; for sel8 events in the z range every evsel bit separately and z_ft0 - z_pv.
 */
//**************************************************************************************************
template <bool IS_MC, typename C>
void ChargedParticles::fillEvSelQA(const C& collision, aod::FT0s const& ft0s, bool isRepresentative)
{
  if (!isNominal || !isRun3) {
    return;
  }
  histos.fill(HIST("QA/events/eventSelection"), 0);
  if (!collision.sel8()) {
    return;
  }
  // cumulative cut-flow of the applied selection
  int stage = 1;
  histos.fill(HIST("QA/events/eventSelection"), stage);
  bool passes = true;
  for (const int bit : evSelBitList) {
    passes = passes && ((collision.selection_raw() >> bit) & 1ULL);
    ++stage;
    if (passes) {
      histos.fill(HIST("QA/events/eventSelection"), stage);
    }
  }
  if (passes && std::abs(collision.posZ()) < maxVtxZ) {
    histos.fill(HIST("QA/events/eventSelection"), stage + 1);
  }
  histos.fill(HIST("QA/events/zVtxMeas"), collision.posZ());
  if (vars.isAcceptedEvent) {
    int runNumber = -1;
    if constexpr (IS_MC) {
      runNumber = collision.template bc_as<BCsWithSel>().runNumber();
    } else {
      runNumber = collision.template bc_as<aod::BCs>().runNumber();
    }
    if (runNumber != lastRunNumber) {
      lastRunNumber = runNumber;
      lastRunLabel = std::to_string(runNumber);
    }
    histos.get<TH1>(HIST("QA/events/runStatistics"))->Fill(lastRunLabel.c_str(), 1.);
    histos.get<TH1>(HIST("QA/events/metaInfo"))->Fill(metaLabel.c_str(), 1.);
  }
  if (std::abs(collision.posZ()) >= maxVtxZ) {
    return;
  }
  const uint64_t selection = collision.selection_raw();
  const bool fillTrue = IS_MC && isRepresentative && vars.isAcceptedEventMC;
  auto fillBit = [&](int bin) {
    histos.fill(HIST("QA/events/multMeasBits"), vars.multMeas, bin);
    if (fillTrue) {
      histos.fill(HIST("QA/events/multTrueBits"), vars.multTrue, bin);
    }
  };
  fillBit(aod::evsel::kNsel); // sel8
  for (int i = 0; i < aod::evsel::kNsel; ++i) {
    if ((selection >> i) & 1ULL) {
      fillBit(i);
    }
  }
  if (vars.isAcceptedEvent) {
    fillBit(aod::evsel::kNsel + 1); // selected
  }
  float dz = -999.f; // no ft0 in the found bc -> underflow
  if (collision.foundFT0Id() >= 0) {
    dz = ft0s.rawIteratorAt(collision.foundFT0Id()).posZ() - collision.posZ();
  }
  histos.fill(HIST("QA/events/dzFT0PV"), vars.multMeas, dz);
}

//**************************************************************************************************
/**
 * Function to processes MC truth info. Assumes initEvent and initEventMC have been called previously.
 */
//**************************************************************************************************
template <typename P>
void ChargedParticles::processTrue(const P& particles)
{
  // untriggered primaries independent of the event class (trigger bias)
  const bool isUntriggeredGoodZPos = vars.isGoodZPosMC && !vars.isTriggeredMC;
  if (!vars.isAcceptedEventMC && !isUntriggeredGoodZPos) {
    return;
  }

  if (vars.isAcceptedEventMC) {
    histos.fill(HIST("multDist_evt_gen"), vars.multTrue);
    if (vars.isTriggeredMC) {
      histos.fill(HIST("multDist_evt_gen_trig"), vars.multTrue);
    }
  }

  for (const auto& particle : particles) {
    if (initParticle(particle) && vars.isChargedPrimary) {
      for (auto i = 0u; i < vars.nRepetitions; ++i) {
        if (vars.isAcceptedEventMC) {
          histos.fill(HIST("multPtSpec_prim_gen"), vars.multTrue, particle.pt());
          if (!vars.isAcceptedEvent) {
            histos.fill(HIST("multPtSpec_prim_gen_evtloss"), vars.multTrue, particle.pt());
          }
        }
        if (isUntriggeredGoodZPos) {
          histos.fill(HIST("multPtSpec_prim_gen_notrig"), vars.multTrue, particle.pt());
        }
      }
    }
  }
}

//**************************************************************************************************
/**
 * Function to process reconstructed data and MC. Assumes initEvent (and initEventMC in case of MC) have been called previously.
 */
//**************************************************************************************************
template <bool IS_MC, typename P, typename T>
void ChargedParticles::processMeas(const T& tracks, bool isRepresentative)
{
  if (!vars.isAcceptedEvent) {
    return;
  }
  // mc: only the representative collision corresponds to the true event
  [[maybe_unused]] const bool isTrueEvent = vars.isAcceptedEventMC && isRepresentative;

  histos.fill(HIST("multDist_evt_meas"), vars.multMeas);

  if constexpr (IS_MC) {
    if (isTrueEvent) {
      histos.fill(HIST("multCorrel_evt"), vars.multMeas, vars.multTrue);
    }
  }

  std::vector<int> foundParticles;
  for (const auto& track : tracks) {
    if (isNominal) {
      fillTrackCutQA(track);
    }
    if (!initTrack(track)) {
      continue;
    }

    // all fills repeated as often as the particle (pcc)
    bool isValidParticle{false};
    const uint32_t nRepetitions = getTrackRepetitions<IS_MC, P>(track, isValidParticle);
    if (isNominal) {
      fillTrackQA(track);
      if constexpr (IS_MC) {
        if (isTrueEvent && isValidParticle && vars.isChargedPrimary) {
          histos.fill(HIST("QA/tracks/deltaPt"), track.pt(), std::abs(track.pt() - track.template mcParticle_as<P>().pt()) / track.pt());
        }
      }
    }
    if constexpr (IS_MC) {
      if (track.has_mcParticle()) {
        foundParticles.push_back(track.mcParticleId());
      }
    }

    for (auto i = 0u; i < nRepetitions; ++i) {
      histos.fill(HIST("multPtSpec_trk_meas"), vars.multMeas, track.pt());

      if constexpr (IS_MC) {
        if (!isTrueEvent) {
          histos.fill(HIST("multPtSpec_trk_meas_evtcont"), vars.multMeas, track.pt());
          continue;
        }

        histos.fill(HIST("multPtSpec_trk_inter"), vars.multTrue, track.pt());
        if (isValidParticle) {
          if (!vars.isChargedPrimary) {
            histos.fill(HIST("multPtSpec_trk_sec_meas"), vars.multMeas, track.pt());
          } else {
            const auto& particle = track.template mcParticle_as<P>();
            histos.fill(HIST("multCorrel_prim"), vars.multMeas, vars.multTrue);
            histos.fill(HIST("ptCorrel_prim"), track.pt(), particle.pt());
            histos.fill(HIST("multPtSpec_prim_meas"), vars.multTrue, particle.pt());
            histos.fill(HIST("multPtSpec_trk_prim_meas"), vars.multMeas, track.pt());
          }
        }
      }
    }
  }

  std::unordered_set<int32_t> uniqueIndices(foundParticles.begin(), foundParticles.end());
  for (const auto& mcParticleID : uniqueIndices) {
    histos.fill(HIST("track_ambiguity"), std::count(foundParticles.begin(), foundParticles.end(), mcParticleID));
  }
}
