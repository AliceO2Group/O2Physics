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

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <Framework/AnalysisTask.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/Track.h>

#include <random>
#include <unordered_set>
#include <vector>

using namespace o2;
using namespace o2::framework;
using aod::track::TrackSelectionFlags;

//--------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------
// Task declaration
//--------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------
struct ChargedParticles {
  std::random_device rnd;
  std::mt19937 gen{rnd()};
  std::uniform_real_distribution<float> dist{0.f, 1.f};
  uint32_t getNRepetitons(float scalingFactor)
  {
    uint32_t nRepetitions = static_cast<uint32_t>(scalingFactor);
    float rest = scalingFactor - nRepetitions;
    if (rest) {
      nRepetitions += (dist(gen) <= rest) ? 1u : 0u;
      // LOGP(info, "scalingFactor: {} -> {}", scalingFactor, nRepetitions);
    }
    return nRepetitions;
  };

  HistogramRegistry histos;
  Service<o2::framework::O2DatabasePDG> pdg;

  Configurable<bool> isRun3{"isRun3", true, "is Run3 dataset"};
  Configurable<uint32_t> maxMultMeas{"maxMultMeas", 100, "max measured multiplicity"};
  Configurable<uint32_t> maxMultTrue{"maxMultTrue", 100, "max true multiplicity"};
  Configurable<float> etaCut{"etaCut", 0.8f, "eta cut"};
  Configurable<float> ptMinCut{"ptMinCut", 0.15f, "pt min cut"};
  Configurable<float> ptMaxCut{"ptMaxCut", 10.f, "pt max cut"};
  Configurable<bool> normINELGT0{"normINELGT0", false, "normalize INEL>0 according to MC"};

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
    bool isINELGT0EventMC{false};
    bool isChargedPrimary{false};
    uint32_t nRepetitions{1u};
  };
  VarContainer vars;
  static constexpr float kMaxVtxZ = 10.f;

  void init(InitContext const&);

  template <typename P>
  bool initParticle(const P& particle);

  template <typename T>
  bool initTrack(const T& track);

  template <bool IS_MC, typename C, typename T>
  void initEvent(const C& collision, const T& tracks);

  template <typename C, typename P>
  void initEventMC(const C& collision, const P& particles);

  template <bool IS_MC, typename T>
  void processMeas(const T& tracks);

  template <typename P>
  void processTrue(const P& particles);

  using CollisionTableData = soa::Join<aod::Collisions, aod::EvSels>;
  using TrackTableData = soa::Join<aod::FullTracks, aod::TracksDCA, aod::TrackSelection>;
  void processData(CollisionTableData::iterator const& collision, TrackTableData const& tracks);
  PROCESS_SWITCH(ChargedParticles, processData, "process data", false);

  using CollisionTableMCTrue = aod::McCollisions;
  using CollisionTableMC = soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions, aod::EvSels>>;
  using TrackTableMC = soa::Join<aod::FullTracks, aod::TracksDCA, aod::TrackSelection, aod::McTrackLabels>;
  using ParticleTableMC = soa::Join<aod::McParticles, aod::ParticleCompositionCorrection>;
  Preslice<TrackTableMC> perCollision = aod::track::collisionId;
  void processMC(CollisionTableMCTrue::iterator const& mcCollision, TrackTableMC const& tracks, CollisionTableMC const& collisions, ParticleTableMC const& particles);
  PROCESS_SWITCH(ChargedParticles, processMC, "process mc", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
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
  std::vector<double> ptBinEdges = {0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75,
                                    0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9,
                                    2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0, 5.5,
                                    6.0, 6.5, 7.0, 8.0, 9.0, 10.0};

  const AxisSpec ptMeasAxis{ptBinEdges, "#it{p}^{ meas}_{T} (GeV/#it{c})", "pt_meas"};

  const int nBinsMultMeas = maxMultMeas + 1;
  const AxisSpec multMeasAxis = {nBinsMultMeas, -0.5, nBinsMultMeas - 0.5, "#it{N}^{ meas}_{ch}", "mult_meas"};

  histos.add("multDist_evt_meas", "", kTH1D, {multMeasAxis});               // measured event distribution (contains contamination from events not in specified class or with wrong vertex position)
  histos.add("multPtSpec_trk_meas", "", kTH2D, {multMeasAxis, ptMeasAxis}); // measured tracks (contains contamination from secondary particles, particles smeared into acceptance and tracks originating from background events as defined above )

  if (doprocessMC) {

    const AxisSpec ptTrueAxis{ptBinEdges, "#it{p}_{T} (GeV/c)", "pt_true"};

    const int nBinsMultTrue = maxMultTrue + 1;
    const AxisSpec multTrueAxis = {nBinsMultTrue, -0.5, nBinsMultTrue - 0.5, "#it{N}_{ch}", "mult_true"};

    histos.add("collision_ambiguity", "", kTH1D, {{6, -0.5, 5.5, "reco collisions per true collision"}}); // log the number of collisions that were reconstructed for a MC collision
    histos.add("track_ambiguity", "", kTH1D, {{6, 0.5, 6.5, "reco tracks per true particle"}});           // log the number of tracks that were reconstructed for a MC particle

    histos.add("multDist_evt_gen", "", kTH1D, {multTrueAxis});      // generated event distribution  (from events within specified class and with proper vertex position)
    histos.add("multDist_evt_gen_trig", "", kTH1D, {multTrueAxis}); // generated event distribution (from events within specified class and with proper vertex position) that in addition fulfils the trigger condition [to disentangle trigger eff from reco eff ]

    histos.add("multCorrel_evt", "", kTH2D, {multMeasAxis, multTrueAxis});  // multiplicity correlation of measured events (excluding background events)
    histos.add("multCorrel_prim", "", kTH2D, {multMeasAxis, multTrueAxis}); // multiplicity correlation of measured primary charged particles (excluding particles from background events)
    histos.add("ptCorrel_prim", "", kTH2D, {ptMeasAxis, ptTrueAxis});       // pT correlation of measured primary charged particles  (excluding particles from background events)

    histos.add("multPtSpec_prim_gen", "", kTH2D, {multTrueAxis, ptTrueAxis});         // generated primary charged particles as function of true properties (from events within specified class and with proper vertex position)
    histos.add("multPtSpec_prim_gen_evtloss", "", kTH2D, {multTrueAxis, ptTrueAxis}); // generated primary charged particles of events that did not pass the event selection as function of multiplicity and pt
    histos.add("multPtSpec_prim_gen_notrig", "", kTH2D, {multTrueAxis, ptTrueAxis});  // generated primary charged particles of events that did not fulfil physics selection and trigger condition as function of multiplicity and pt
    histos.add("multPtSpec_prim_meas", "", kTH2D, {multTrueAxis, ptTrueAxis});        // measured primary charged particles as function of true properties (no contamination from background events)

    histos.add("multPtSpec_trk_prim_meas", "", kTH2D, {multMeasAxis, ptMeasAxis});    // tracks from measured primaries (no contamination from secondaries, particles smeared into acceptance or background events)
    histos.add("multPtSpec_trk_sec_meas", "", kTH2D, {multMeasAxis, ptMeasAxis});     // tracks from measured secondaries (no contamination from particles smeared into acceptance or background events)  [for QA to disentangle secondaries from other contamination]
    histos.add("multPtSpec_trk_meas_evtcont", "", kTH2D, {multMeasAxis, ptMeasAxis}); // tracks from events that are measured, but do not belong to the desired class of events
    histos.add("multPtSpec_trk_inter", "", kTH2D, {multTrueAxis, ptMeasAxis});
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
void ChargedParticles::processData(CollisionTableData::iterator const& collision, TrackTableData const& tracks)
{
  initEvent<false>(collision, tracks);
  processMeas<false>(tracks);
}

//**************************************************************************************************
/**
 * Entrypoint to processes mc.
 */
//**************************************************************************************************
void ChargedParticles::processMC(CollisionTableMCTrue::iterator const& mcCollision, TrackTableMC const& tracks, CollisionTableMC const& collisions, ParticleTableMC const& particles)
{
  histos.fill(HIST("collision_ambiguity"), collisions.size());

  // TODO: process only most probable collision (run3)
  if (collisions.size() > 1) {
    // FIXME: for now skip all ambiguously reconstructed collisions as we do not know what to do with them
    return;
  }
  // MEMO: this ambiguity of the reconstructed collisions raises several questions:
  // - how to select most probable collision?
  // - how to avoid double or triple counting of an actual collision in data (or how to treat this as additional contamination of the measurement based on MC info)
  // - how does this pollute the event reconstruction efficiency

  initEventMC(mcCollision, particles);
  if (collisions.size() == 0) {
    vars.isAcceptedEvent = false;
  } else {
    for (const auto& collision : collisions) {
      auto curTracks = tracks.sliceBy(perCollision, collision.globalIndex());
      initEvent<true>(collision, curTracks);
      processMeas<true>(curTracks);
      break; // for now look only at first collision...
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
  auto pdgParticle = pdg->GetParticle(particle.pdgCode());
  if (!pdgParticle || pdgParticle->Charge() == 0.) {
    return false;
  }
  vars.isChargedPrimary = particle.isPhysicalPrimary();

  // event class is INEL>0 in case it has a charged particle in abs(eta) < 1
  vars.isINELGT0EventMC = vars.isINELGT0EventMC || (vars.isChargedPrimary && (std::abs(particle.eta()) < 1.));

  if (std::abs(particle.eta()) >= etaCut) {
    return false;
  }
  if (particle.pt() <= ptMinCut || particle.pt() >= ptMaxCut) {
    return false;
  }

  float pccWeight = particle.pccWeight();
  if (systMode == kSystDownPCC) {
    pccWeight = particle.pccWeightSysDown();
  } else if (systMode == kSystUpPCC) {
    pccWeight = particle.pccWeightSysUp();
  }
  vars.nRepetitions = getNRepetitons(pccWeight);
  // FIXME: in case of nRepetitions = 0 INELGT0 can be wrong
  return true;
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
template <bool IS_MC, typename C, typename T>
void ChargedParticles::initEvent(const C& collision, const T& tracks)
{
  vars.multMeas = 0;
  for (const auto& track : tracks) {
    if (initTrack(track)) {
      if constexpr (IS_MC) {
        if (!track.has_mcParticle()) {
          continue;
        }
        const auto& particle = track.template mcParticle_as<ParticleTableMC>();
        if (!initParticle(particle)) {
          continue;
        }
        vars.multMeas += vars.nRepetitions;
      } else {
        ++vars.multMeas;
      }
    }
  }

  vars.isAcceptedEvent = false;
  if (std::abs(collision.posZ()) < kMaxVtxZ) {
    if (isRun3) {
      if (collision.sel8() &&
          collision.selection_bit(aod::evsel::kNoSameBunchPileup) &&
          collision.selection_bit(aod::evsel::kIsVertexITSTPC) &&
          collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV))
      // !collision.selection_bit(o2::aod::evsel::kIsVertexTRDmatched) &&
      // collision.selection_bit(aod::evsel::kNoCollInTimeRangeStandard) &&
      // collision.selection_bit(aod::evsel::kNoCollInRofStandard) &&
      // collision.selection_bit(aod::evsel::kIsVertexTOFmatched)
      {
        vars.isAcceptedEvent = true;
      }
    } else {
      if (collision.sel7() && ((doprocessMC) ? true : collision.alias_bit(kINT7))) {
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
  vars.isINELGT0EventMC = false; // will be set to true in case a charged particle within eta +-1 is found
  vars.multTrue = 0;
  for (const auto& particle : particles) {
    if (!initParticle(particle) || !vars.isChargedPrimary) {
      continue;
    }
    vars.multTrue += vars.nRepetitions;
  }
  bool isGoodEventClass = (normINELGT0) ? vars.isINELGT0EventMC : (vars.multTrue > 0);
  vars.isAcceptedEventMC = isGoodEventClass && (std::abs(collision.posZ()) < kMaxVtxZ);
}

//**************************************************************************************************
/**
 * Function to processes MC truth info. Assumes initEvent and initEventMC have been called previously.
 */
//**************************************************************************************************
template <typename P>
void ChargedParticles::processTrue(const P& particles)
{
  if (!vars.isAcceptedEventMC) {
    return;
  }

  histos.fill(HIST("multDist_evt_gen"), vars.multTrue);

  for (const auto& particle : particles) {
    if (initParticle(particle) && vars.isChargedPrimary) {
      for (auto i = 0u; i < vars.nRepetitions; ++i) {
        histos.fill(HIST("multPtSpec_prim_gen"), vars.multTrue, particle.pt());
        if (!vars.isAcceptedEvent) {
          histos.fill(HIST("multPtSpec_prim_gen_evtloss"), vars.multTrue, particle.pt());
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
template <bool IS_MC, typename T>
void ChargedParticles::processMeas(const T& tracks)
{
  if (!vars.isAcceptedEvent) {
    return;
  }

  histos.fill(HIST("multDist_evt_meas"), vars.multMeas);

  if constexpr (IS_MC) {
    if (vars.isAcceptedEventMC) {
      histos.fill(HIST("multCorrel_evt"), vars.multMeas, vars.multTrue);
    }
  }

  std::vector<int> foundParticles;
  for (const auto& track : tracks) {
    if (!initTrack(track)) {
      continue;
    }

    histos.fill(HIST("multPtSpec_trk_meas"), vars.multMeas, track.pt());

    if constexpr (IS_MC) {
      if (!track.has_mcParticle()) {
        continue;
      }

      /*
      if(count(foundParticles.begin(), foundParticles.end(), track.mcParticleId()){
        //LOGP(info, "Multiple tracks reconstructed for particle with label {}", track.mcParticleId());
       // for now only consider first particle that is found...
        continue;
      }
      */

      foundParticles.push_back(track.mcParticleId());

      const auto& particle = track.template mcParticle_as<ParticleTableMC>();

      if (!vars.isAcceptedEventMC) {
        histos.fill(HIST("multPtSpec_trk_meas_evtcont"), vars.multMeas, track.pt());
        continue;
      }

      histos.fill(HIST("multPtSpec_trk_inter"), vars.multTrue, track.pt());
      if (initParticle(particle)) {
        for (auto i = 0u; i < vars.nRepetitions; ++i) {
          if (!vars.isChargedPrimary) {
            histos.fill(HIST("multPtSpec_trk_sec_meas"), vars.multMeas, track.pt());
          } else {
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
