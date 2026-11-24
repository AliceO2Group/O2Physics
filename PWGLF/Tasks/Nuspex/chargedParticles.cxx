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
#include "Common/CCDB/ctpRateFetcher.h"
#include "Common/Core/MetadataHelper.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <DataFormatsCTP/Configuration.h>
#include <DataFormatsCTP/Scalers.h>
#include <DataFormatsParameters/GRPLHCIFData.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/RawDeviceService.h>
#include <Framework/StaticFor.h>
#include <Framework/runDataProcessing.h>

#include <TAxis.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>

#include <fairmq/Device.h>

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <initializer_list>
#include <map>
#include <string>
#include <string_view>
#include <type_traits>
#include <unordered_map>
#include <utility>
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
struct ChargedParticles { // NOLINT(clang-analyzer-optin.performance.Padding): one task instance, member order kept readable
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
    for (const auto& val : {particle.px(), particle.py(), particle.pz(), particle.vx(), particle.vy(), particle.vz()}) {
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
  Service<o2::ccdb::BasicCCDBManager> ccdb{}; // ctp scalers for the interaction rate
  o2::ctpRateFetcher rateFetcher;

  enum DataRun : int { kRunAuto = 0,
                       kRun2 = 2,
                       kRun3 = 3 };
  Configurable<int> dataRun{"dataRun", kRunAuto, "lhc run of the data: 0 = auto (from ao2d metadata), 2 = run 2 (converted), 3 = run 3"};
  bool isRun3{true}; // from dataRun
  Configurable<std::string> collisionSystem{"collisionSystem", "pp", "collision system (sets multiplicity axes, event selection default and interaction rate source): pp, pPb, pO, OO, NeNe, XeXe, PbPb"};
  int maxMultMeas{100}; // from collisionSystem
  int maxMultTrue{100};
  std::string irSource{"T0VTX"}; // interaction rate estimator (ctp rate fetcher), from collisionSystem
  Configurable<float> etaCut{"etaCut", 0.8f, "eta cut"};
  Configurable<std::string> ptRange{"ptRange", "nominal", "pt window and binning of tracks and particles: nominal (0.15-10), lowPt (0.075-10, low field), highPt (0.15-50)"};
  Configurable<bool> fiducialMult{"fiducialMult", true, "count n_ch (measured and true, and so the event class fiducial) only in the fiducial window 0.15 < pt < 10, whatever ptRange"};
  Configurable<int> multEstimator{"multEstimator", 0, "measured multiplicity: 0 = selected tracks (ITS-TPC), 1 = also ITS-only tracks that pass all other selections (spectra stay ITS-TPC)"};
  static constexpr float FiducialPtMin = 0.15f;
  static constexpr float FiducialPtMax = 10.f;
  float ptMin{FiducialPtMin}; // from ptRange
  float ptMax{FiducialPtMax};
  bool isInMultWindow(float pt) const { return !fiducialMult.value || (pt > FiducialPtMin && pt < FiducialPtMax); }
  Configurable<std::string> eventClasses{"eventClasses", "fiducial", "mc event classes, comma separated: fiducial (N_ch > 0, N_ch as set by ptRange and fiducialMult), inelgt0, triggered (tvx; run 2: v0and, 2010/2011: int1), inel; one class: flat output, several: one subfolder per class"};
  enum EventClass { kFiducial,
                    kINELGT0,
                    kTriggered,
                    kINEL,
                    kNClasses };
  static constexpr std::array<std::string_view, kNClasses + 1> ClassDir = {"fiducial/", "inelgt0/", "triggered/", "inel/", ""}; // last: flat output
  static constexpr std::array<std::string_view, kNClasses> ClassName = {"fiducial", "inelgt0", "triggered", "inel"};
  uint8_t enabledDirs{0}; // output directories in use (bit kNClasses: flat output for a single class)
  int flatClass{kFiducial};
  // call f(std::integral_constant<int, dir>) for every output directory in use
  template <typename F>
  void forEachDir(F const& f)
  {
    static_for<0, kNClasses>([&](auto i) {
      constexpr int Dir = i.value;
      if (((enabledDirs >> Dir) & 1) != 0) {
        f(std::integral_constant<int, Dir>{});
      }
    });
  }
  // is the current mc event in the class of output directory dir
  bool isInDir(int dir) const { return ((vars.classMask >> (dir == kNClasses ? flatClass : dir)) & 1) != 0; }
  static constexpr int NStagesMC = 8;
  static constexpr int NLayersITS = 7;
  std::string metaLabel;      // qa: ao2d metadata label
  int lastRunNumber{-1};      // qa: cached run label
  TProfile* rateRun{nullptr}; // qa: interaction rate per run (run 3, nominal)
  bool hasRateInfo{false};    // qa: ccdb objects of the rate fetcher found for the current run
  struct ParticleTracks {
    float pt{0.f};
    int nTracks{0};
    int nCollisions{0};
    int64_t lastCollision{-1};
  };
  std::unordered_map<int64_t, ParticleTracks> tracksPerParticle; // qa: selected tracks per mc particle in the current mc event
  std::string lastRunLabel;
  int phiRunBin{0};           // qa: run bin of qa/tracks/phiRun (0 = no run label yet)
  int lastEvSelRunNumber{-1}; // qa: run label of the per-run cut flow (all candidates, not only accepted events)
  std::string evSelRunLabel;
  std::string mcRunLabel;                           // qa: run of the current mc collision (per-run efficiencies)
  std::array<double, kNClasses + 1> primGenRunX{};  // qa: x of the run bin of qa/tracks/primGenRun per class, set once per mc collision
  std::array<double, kNClasses + 1> primRecoRunX{}; // qa: same for qa/tracks/primRecoRun, set once per reco collision
  Configurable<std::string> evSelBits{"evSelBits", "auto", "event selection bits on top of sel8 (names as in EventSelectionParams), comma separated; auto = default of the collision system (currently kNoSameBunchPileup for all)"};
  Configurable<float> maxVtxZ{"maxVtxZ", 10.f, "max |z| of the reco and true vertex"};
  // own dca cuts |dca| < k * sigma(pt), sigma = p0 + p1 / pt^p2 (cm), k fixed for the analysis; run 3: measured
  // resolution; run 2 (converted data): the dca cuts of the central global track selection (dca_xy 0.0105 + 0.035 /
  // pt^1.1, dca_z 2 cm)
  struct DcaPars {
    float p0{0.f};
    float p1{0.f};
    float p2{1.f};
  };
  static constexpr DcaPars DefaultDcaXY{.p0 = 0.0015f, .p1 = 0.005f, .p2 = 1.1f};
  static constexpr DcaPars DefaultDcaZ{.p0 = 0.2f, .p1 = 0.f, .p2 = 1.f};
  static constexpr float NSigmaDcaXY = 7.f; // as in run 2; variations 4, 10
  static constexpr float NSigmaDcaZ = 10.f; // variations 6, 14
  DcaPars dcaCutXY{};                       // cut = p0 + p1 / pt^p2, i.e. k * sigma (from init)
  DcaPars dcaCutZ{};
  // run 3 dca resolution: core width of pp 13.6 tev data (2026, lhc26ac), used for all systems (within 20 % in pp at
  // low field, o-o and pb-pb mc)
  static constexpr DcaPars DcaResRun3XY{.p0 = 0.0006405f, .p1 = 0.002207f, .p2 = 0.950f};
  static constexpr DcaPars DcaResRun3Z{.p0 = 0.0005770f, .p1 = 0.002423f, .p2 = 0.9414f};
  uint64_t evSelMask{0};
  std::vector<int> evSelBitList; // configured bits in the given order (cut-flow)
  std::vector<int> trackCutBits; // single cuts of the track selection mask (track qa)
  uint32_t tpcCutMask{0};        // qa: tpc cuts of the selection mask (its-tpc matching)
  uint32_t itsCutMask{0};        // qa: its cuts of the selection mask (its-tpc matching)
  bool isNominal{true};          // qa only for the nominal setting
  Configurable<uint32_t> seedOffset{"seedOffset", 0, "offset for the (deterministic) random numbers of the particle composition correction"};

  // systematic variations, numbering as in run 2 where the cut exists (unfolding framework expects these ids);
  // 107: max shared tpc cluster fraction 0.4 (nominal: no cut; one-sided, 108 unused);
  // 109-110 (golden chi2) and 118-119 (geometric cut) do not exist in run 3 (fatal);
  // 116-117 are tpc crossed rows (run 2: geometric length cut); 130+ are new in run 3
  enum : uint32_t {
    kSystNominal = 100,
    kSystDownChi2PerClusterITS,      // 101
    kSystUpChi2PerClusterITS,        // 102
    kSystDownChi2PerClusterTPC,      // 103
    kSystUpChi2PerClusterTPC,        // 104
    kSystDownTPCCrossedRowsOverNCls, // 105
    kSystUpTPCCrossedRowsOverNCls,   // 106
    kSystDownTPCFracSharedCls,       // 107
    kSystDownDCAxy = 111,            // 4 sigma (nominal 7)
    kSystUpDCAxy,                    // 112, 10 sigma
    kSystDownDCAz,                   // 113, 6 sigma (nominal 10)
    kSystUpDCAz,                     // 114, 14 sigma
    kSystITSHits,                    // 115, only relevant for converted data
    kSystDownTPCCrossedRows,         // 116
    kSystUpTPCCrossedRows,           // 117
    kSystDownPCC = 120,
    kSystUpPCC,                 // 121
    kSystDownVtxZ = 130,        // |z| < 7 cm (reco and true vertex)
    kSystEvSelNoCollInRof,      // 131, + kNoCollInRofStandard (run 3 only)
    kSystEvSelNoCollInTimeRange // 132, + kNoCollInTimeRangeStandard (run 3 only)
  };
  Configurable<uint32_t> systMode{"systMode", kSystNominal, "variation for systematic uncertainties"};
  float vtxZCut() const { return (systMode.value == kSystDownVtxZ) ? 7.f : maxVtxZ.value; } // applied |z| cut (reco and true vertex)
  uint16_t trackSelMask{TrackSelectionFlags::kGlobalTrackWoPtEta};                          // track selection bitmask (without cut that is being varied)
  uint16_t cutVarFlag{0};
  bool isCutVar{false}; // a single cut is varied (checked with trackSel)
  TrackSelection trackSel;
  TrackSelection::TrackCuts trackSelFlag{};
  static bool isINT1Period(int runNumber)
  {
    // int1 trigger (spd or v0a or v0c) in 2010 and part of 2011: run ranges as in the event selection
    return runNumber <= 136377 || (runNumber >= 144871 && runNumber <= 159582); // o2-linter: disable=magic-number (int1 run ranges)
  }

  // helper struct to store transient properties
  struct VarContainer {
    uint32_t multMeas{0u};
    uint32_t multTrue{0u};
    bool isAcceptedEvent{false};
    uint8_t classMask{0};      // mc event classes the event belongs to (each with good true z)
    bool isGoodZPosMC{false};  // true vertex within |z| cut
    bool isTriggeredMC{false}; // tvx in the bc of the mc collision
    bool isINELGT0EventMC{false};
    bool isChargedPrimary{false};
    uint32_t nRepetitions{1u};
    int64_t mcCollisionIndex{-1}; // current mc collision
    float zVtx{0.f};              // qa: reco vertex z (its-tpc matching)
  };
  VarContainer vars;

  void init(InitContext& ic);

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
  template <typename T>
  void fillMatchQA(const T& tracks);
  template <typename T>
  bool isItsOnlyCandidate(const T& track);
  template <typename T>
  bool passesDcaXY(const T& track) const;
  template <typename T>
  bool passesDcaZ(const T& track) const;
  static bool isTpcCut(TrackSelection::TrackCuts cut)
  {
    return cut == TrackSelection::TrackCuts::kTPCNCls || cut == TrackSelection::TrackCuts::kTPCCrossedRows || cut == TrackSelection::TrackCuts::kTPCCrossedRowsOverNCls || cut == TrackSelection::TrackCuts::kTPCChi2NDF || cut == TrackSelection::TrackCuts::kTPCRefit || cut == TrackSelection::TrackCuts::kTPCFracSharedCls;
  }

  template <typename P>
  void processTrue(const P& particles);

  template <typename T>
  void fillTrackQA(const T& track);

  template <typename T>
  void fillTrackCutQA(const T& track);

  template <bool IS_MC, typename C>
  void fillEvSelQA(const C& collision, aod::FT0s const& ft0s, bool isRepresentative = false);
  bool hasCtpRateInfo(uint64_t timestamp, int runNumber);

  using CollisionTableData = soa::Join<aod::Collisions, aod::EvSels>;
  using TrackTableData = soa::Join<aod::FullTracks, aod::TracksCov, aod::TracksDCA, aod::TrackSelection>;
  void processData(CollisionTableData::iterator const& collision, TrackTableData const& tracks, aod::FT0s const& ft0s, aod::BCsWithTimestamps const&);
  PROCESS_SWITCH(ChargedParticles, processData, "process data", false);

  using CollisionTableMCTrue = aod::McCollisions;
  using CollisionTableMC = soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions, aod::EvSels>>;
  using TrackTableMC = soa::Join<aod::FullTracks, aod::TracksCov, aod::TracksDCA, aod::TrackSelection, aod::McTrackLabels>;
  using ParticleTableMC = soa::Join<aod::McParticles, aod::ParticleCompositionCorrection>;
  using BCsWithSel = soa::Join<aod::BCs, aod::BcSels, aod::Timestamps>; // trigger info of the mc collision's bc; timestamp for the interaction rate
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
void ChargedParticles::init(InitContext& ic)
{
  if (static_cast<int>(doprocessData.value) + static_cast<int>(doprocessMC.value) + static_cast<int>(doprocessMCnoPCC.value) != 1) {
    LOGP(fatal, "enable exactly one of processData, processMC, processMCnoPCC");
  }

  // lhc run: from the ao2d metadata unless set explicitly; a contradiction stops the task
  int runFromMetadata = kRunAuto;
  if (metadataInfo.isInitialized()) {
    const std::string run = metadataInfo.get("Run");
    if (run == "2" || run == "3") {
      runFromMetadata = std::stoi(run);
    }
  }
  if (dataRun != kRunAuto && dataRun != kRun2 && dataRun != kRun3) {
    LOGP(fatal, "invalid dataRun {} (0 = auto, 2, 3)", dataRun.value);
  }
  if (dataRun == kRunAuto && runFromMetadata == kRunAuto) {
    LOGP(fatal, "dataRun = 0 (auto) but the ao2d metadata has no run: set dataRun to 2 or 3");
  }
  if (dataRun != kRunAuto && runFromMetadata != kRunAuto && dataRun != runFromMetadata) {
    LOGP(fatal, "dataRun {} contradicts the ao2d metadata (run {})", dataRun.value, runFromMetadata);
  }
  isRun3 = ((dataRun != kRunAuto) ? dataRun.value : runFromMetadata) == kRun3;
  LOGP(info, "data from run {} ({})", isRun3 ? 3 : 2, (dataRun != kRunAuto) ? "dataRun" : "ao2d metadata");
  if (!isRun3 && evSelBits.value != "auto") {
    LOGP(fatal, "evSelBits '{}' is not applied to run 2 data (run 2 uses sel7 and kINT7): keep evSelBits = auto", evSelBits.value);
  }

  // collision system: multiplicity axis ranges (measured, true) and default event selection
  // pp, p-pb: as in run 2; p-o, o-o, ne-ne: cover the 2025 data; pb-pb, xe-xe: extended to cover central pb-pb collisions at 5.36 tev
  // event selection default: kNoSameBunchPileup for all systems (pb-pb under study)
  std::string evSelBitsDefault = "kNoSameBunchPileup";
  if (collisionSystem.value == "pp") {
    maxMultMeas = 100;
    maxMultTrue = 100;
  } else if (collisionSystem.value == "pPb") {
    maxMultMeas = 200;
    maxMultTrue = 200;
  } else if (collisionSystem.value == "pO") {
    maxMultMeas = 300;
    maxMultTrue = 350;
  } else if (collisionSystem.value == "OO" || collisionSystem.value == "NeNe") {
    maxMultMeas = 400;
    maxMultTrue = 600;
  } else if (collisionSystem.value == "XeXe" || collisionSystem.value == "PbPb") {
    maxMultMeas = 4000;
    maxMultTrue = 5000;
  } else {
    LOGP(fatal, "unknown collisionSystem '{}' (pp, pPb, pO, OO, NeNe, XeXe, PbPb)", collisionSystem.value);
  }
  irSource = (collisionSystem.value == "PbPb") ? "ZNC hadronic" : "T0VTX";

  // event classes: one class -> flat output, several -> one subfolder per class
  int nClasses = 0;
  uint8_t classes = 0;
  for (size_t start = 0; start < eventClasses.value.size();) {
    size_t end = eventClasses.value.find(',', start);
    if (end == std::string::npos) {
      end = eventClasses.value.size();
    }
    std::string name = eventClasses.value.substr(start, end - start);
    name.erase(0, name.find_first_not_of(' '));
    name.erase(name.find_last_not_of(' ') + 1);
    start = end + 1;
    int cls = -1;
    for (int c = 0; c < kNClasses; ++c) {
      if (name == ClassName[c]) {
        cls = c;
      }
    }
    if (cls < 0 || ((classes >> cls) & 1) != 0) {
      LOGP(fatal, "unknown or duplicate event class '{}' in eventClasses (fiducial, inelgt0, triggered, inel)", name);
    }
    classes |= (1 << cls);
    flatClass = cls;
    ++nClasses;
  }
  if (nClasses == 0) {
    LOGP(fatal, "no event class given in eventClasses");
  }
  enabledDirs = (nClasses == 1) ? (1 << kNClasses) : classes;

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
  // event selection variations: one more bit on top of the configured ones
  if (systMode == kSystEvSelNoCollInRof || systMode == kSystEvSelNoCollInTimeRange) {
    const int bit = (systMode == kSystEvSelNoCollInRof) ? aod::evsel::kNoCollInRofStandard : aod::evsel::kNoCollInTimeRangeStandard;
    if (!isRun3) {
      LOGP(fatal, "systMode {} (event selection variation) exists only for run 3", systMode.value);
    }
    if (((evSelMask >> bit) & 1ULL) != 0) {
      LOGP(fatal, "systMode {}: {} is already in evSelBits (variation would equal nominal)", systMode.value, aod::evsel::selectionLabels[bit]);
    }
    evSelMask |= (1ULL << bit);
    evSelBitList.push_back(bit);
  }

  std::vector<double> ptBinEdges = {0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75,
                                    0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9,
                                    2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0, 5.5,
                                    6.0, 6.5, 7.0, 8.0, 9.0, 10.0};
  // pt window = first and last edge
  if (ptRange.value == "lowPt") {
    ptBinEdges.insert(ptBinEdges.begin(), {0.075, 0.1, 0.125});
  } else if (ptRange.value == "highPt") {
    ptBinEdges.insert(ptBinEdges.end(), {11., 12., 13., 14., 15., 16., 18., 20., 22., 24., 26., 30., 34., 40., 50.});
  } else if (ptRange.value != "nominal") {
    LOGP(fatal, "unknown ptRange '{}' (nominal, lowPt, highPt)", ptRange.value);
  }
  ptMin = ptBinEdges.front();
  ptMax = ptBinEdges.back();

  // relative pt resolution (qa): logarithmic bins, 6 % wide, from 1e-4 to 1
  constexpr int NRelResBins = 160;
  constexpr double RelResLog10Min = -4.;
  constexpr double RelResLog10Max = 0.;
  std::vector<double> relResEdges;
  for (int i = 0; i <= NRelResBins; ++i) {
    relResEdges.push_back(std::pow(10., RelResLog10Min + (RelResLog10Max - RelResLog10Min) * i / NRelResBins));
  }

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
    histos.add("QA/events/eventSelectionRun", "", kTH2D, {{512, 0., 512.}, {nStages, -0.5, nStages - 0.5}}); // run axis labelled as runStatistics
    histos.add("QA/events/runStatistics", "", kTH1D, {{512, 0., 512.}});
    if (isRun3) {
      rateRun = histos.add<TProfile>("QA/events/rateRun", "", kTProfile, {{512, 0., 512.}}).get(); // mean interaction rate (khz) of accepted events per run; run axis labelled as runStatistics
    }
    histos.add("QA/events/metaInfo", "", kTH1D, {{4, 0., 4.}});
    // track qa (selected tracks in selected events), run 2 names and binning (7 its layers, run 3 dca cut ranges)
    histos.add("QA/tracks/signed1Pt", "", kTH1D, {{200, -8., 8., "q/#it{p}_{T}"}});
    histos.add("QA/tracks/eta", "", kTH1D, {{180, -0.9, 0.9, "#eta"}});
    histos.add("QA/tracks/phi", "", kTH1D, {{180, 0., o2::constants::math::TwoPI, "#phi (rad)"}});
    histos.add("QA/tracks/etaPhi", "", kTH2D, {{16, -0.8, 0.8, "#eta"}, {180, 0., o2::constants::math::TwoPI, "#phi (rad)"}}); // acceptance holes: tpc (phi stripes) vs its (eta-phi patches)
    const AxisSpec matchPtAxis{std::vector<double>{0.15, 0.5, 1., 2., 5., 10.}, "#it{p}_{T} of a selected track (GeV/#it{c})"};
    const AxisSpec matchCandAxis{3, 1.5, 4.5, "candidates (selected + ITS-only)"};
    const AxisSpec matchOtherAxis{std::vector<double>{0.15, 0.5, 10.}, "#it{p}_{T} of another candidate (GeV/#it{c})"};
    histos.add("QA/tracks/matchOthersSelected", "", kTH3D, {matchPtAxis, matchCandAxis, matchOtherAxis});              // its-tpc matching of the other tracks vs a selected track, collisions with 2-4 candidates: other candidate selected
    histos.add("QA/tracks/matchOthersItsOnly", "", kTH3D, {matchPtAxis, matchCandAxis, matchOtherAxis});               // same, other candidate its-only
    histos.add("QA/tracks/phiRun", "", kTH2D, {{180, 0., o2::constants::math::TwoPI, "#phi (rad)"}, {512, 0., 512.}}); // run axis labelled as runStatistics
    histos.add("QA/tracks/itsNCls", "", kTH1D, {{8, -0.5, 7.5, "# clusters ITS"}});
    histos.add("QA/tracks/itsHits", "", kTH1D, {{7, -0.5, 6.5, "layer ITS"}});
    histos.add("QA/tracks/itsChi2NCl", "", kTH1D, {{100, 0., 40., "chi2 per cluster ITS"}});
    histos.add("QA/tracks/tpcNClsFindable", "", kTH1D, {{165, -0.5, 164.5, "# findable clusters TPC"}});
    histos.add("QA/tracks/tpcNClsFound", "", kTH1D, {{165, -0.5, 164.5, "# clusters TPC"}});
    histos.add("QA/tracks/tpcCrossedRows", "", kTH1D, {{165, -0.5, 164.5, "# crossed rows TPC"}});
    histos.add("QA/tracks/tpcCrossedRowsOverFindableCls", "", kTH1D, {{90, 0.7, 1.6, "crossed rows / findable clusters TPC"}});
    histos.add("QA/tracks/tpcFractionSharedCls", "", kTH2D, {ptMeasAxis, {51, 0., 1.02, "fraction of shared clusters TPC"}});
    histos.add("QA/tracks/tpcChi2NCl", "", kTH1D, {{100, 0., 10., "chi2 / cluster TPC"}});
    histos.add("QA/tracks/sigmaPt", "", kTH2D, {ptMeasAxis, {relResEdges, "#sigma(#it{p}^{ meas}_{T}) / #it{p}^{ meas}_{T}"}});
    // one dca axis for all dca qa: a bin centred at 0, fine in the core (resolution down to ~10 um at high pt), coarse
    // in the tails; segments of {bin width (cm), bins per side}, each width 5x the previous and every boundary on the
    // coarser grids, so RebinY(5) and RebinY(25) give the next levels (pure, centred at 0)
    const std::vector<double> dcaEdges = [](double centralWidth, std::initializer_list<std::pair<double, int>> segments) {
      std::vector<double> edges{centralWidth / 2.};
      for (const auto& [width, nBins] : segments) {
        const double inner = edges.back();
        for (int j = 1; j <= nBins; ++j) {
          edges.push_back(inner + j * width);
        }
      }
      std::vector<double> all;
      for (auto it = edges.rbegin(); it != edges.rend(); ++it) {
        all.push_back(-*it);
      }
      all.insert(all.end(), edges.begin(), edges.end());
      return all;
    }(0.00025, {{0.00025, 62}, {0.00125, 75}, {0.00625, 75}, {0.03125, 100}, {0.15625, 50}}); // 2.5 um to 1.5625 mm, 725 bins in +-11.5 cm
    const AxisSpec dcaXYAxis{dcaEdges, "dca XY (cm)"};
    const AxisSpec dcaZAxis{dcaEdges, "dca Z (cm)"};
    histos.add("QA/tracks/dcaXY", "", kTH2D, {ptMeasAxis, dcaXYAxis});
    histos.add("QA/tracks/dcaZ", "", kTH2D, {ptMeasAxis, dcaZAxis});
    // track cut qa: n-1 and cumulative per single cut of the selection mask, n-1 dca distributions (no dca cut)
    // (nominal only, so the mask does not contain a varied cut yet)
    for (int i = 0; i < static_cast<int>(TrackSelection::TrackCuts::kNCuts); ++i) {
      if ((trackSelMask & (1 << i)) != 0) {
        trackCutBits.push_back(i);
      }
    }
    const int nCutBins = 2 * static_cast<int>(trackCutBits.size()) + 3;
    histos.add("QA/tracks/trackCuts", "", kTH2D, {ptMeasAxis, {nCutBins, -0.5, nCutBins - 0.5}});
    histos.add("QA/tracks/dcaXYnm1", "", kTH2D, {ptMeasAxis, dcaXYAxis});
    histos.add("QA/tracks/dcaZnm1", "", kTH2D, {ptMeasAxis, dcaZAxis});
    // its-tpc matching: tracks passing the tpc cuts of the mask (no track type, its or dca cut), with an its track, and passing the its cuts
    for (const auto& bit : trackCutBits) {
      const auto cut = static_cast<TrackSelection::TrackCuts>(bit);
      if (cut == TrackSelection::TrackCuts::kTPCNCls || cut == TrackSelection::TrackCuts::kTPCCrossedRows || cut == TrackSelection::TrackCuts::kTPCCrossedRowsOverNCls || cut == TrackSelection::TrackCuts::kTPCChi2NDF || cut == TrackSelection::TrackCuts::kTPCRefit) {
        tpcCutMask |= (1u << bit);
      } else if (cut == TrackSelection::TrackCuts::kITSNCls || cut == TrackSelection::TrackCuts::kITSChi2NDF || cut == TrackSelection::TrackCuts::kITSRefit || cut == TrackSelection::TrackCuts::kITSHits) {
        itsCutMask |= (1u << bit);
      }
    }
    histos.add("QA/tracks/itsTpcMatch", "", kTH2D, {ptMeasAxis, {5, -0.5, 4.5}});
  }

  if (doprocessMC || doprocessMCnoPCC) {

    const AxisSpec ptTrueAxis{ptBinEdges, "#it{p}_{T} (GeV/c)", "pt_true"};

    const int nBinsMultTrue = maxMultTrue + 1;
    const AxisSpec multTrueAxis = {nBinsMultTrue, -0.5, nBinsMultTrue - 0.5, "#it{N}_{ch}", "mult_true"};
    // per-run qa: unit bins up to n_ch = 100, then about 50 equal bins (keeps the 512-run histograms small in pb-pb)
    constexpr int MaxUnitBinsRun = 100;
    constexpr int NCoarseBinsRun = 50;
    std::vector<double> runMultEdges;
    for (int i = 0; i <= std::min(maxMultTrue, MaxUnitBinsRun) + 1; ++i) {
      runMultEdges.push_back(i - 0.5);
    }
    if (maxMultTrue > MaxUnitBinsRun) {
      const int step = (maxMultTrue - MaxUnitBinsRun + NCoarseBinsRun - 1) / NCoarseBinsRun;
      while (runMultEdges.back() < maxMultTrue + 0.5) {
        runMultEdges.push_back(runMultEdges.back() + step);
      }
    }
    const AxisSpec runMultTrueAxis{runMultEdges, "#it{N}_{ch}", "mult_true_run"};

    histos.add("collision_ambiguity", "", kTH1D, {{6, -0.5, 5.5, "reco collisions per true collision"}}); // log the number of collisions that were reconstructed for a MC collision
    histos.add("track_ambiguity", "", kTH1D, {{6, 0.5, 6.5, "reco tracks per true particle"}});           // log the number of tracks that were reconstructed for a MC particle

    if (isNominal) {
      histos.add("QA/events/zVtxGen", "", kTH1D, {{61, -30.5, 30.5, "#it{V}_{z} (cm)"}});
      histos.add("QA/tracks/tracksPerParticle", "", kTH2D, {ptTrueAxis, {6, 0.5, 6.5, "selected tracks per particle"}});
      histos.add("QA/tracks/collisionsPerParticle", "", kTH2D, {ptTrueAxis, {4, 0.5, 4.5, "reco collisions per particle"}});
    }
    // class-dependent histograms: flat (one class) or one subfolder per class
    for (int d = 0; d <= kNClasses; ++d) {
      if (((enabledDirs >> d) & 1) == 0) {
        continue;
      }
      const std::string dir(ClassDir[d]);
      histos.add((dir + "multDist_evt_gen").c_str(), "", kTH1D, {multTrueAxis});                        // generated event distribution  (from events within specified class and with proper vertex position)
      histos.add((dir + "multDist_evt_gen_trig").c_str(), "", kTH1D, {multTrueAxis});                   // generated event distribution (from events within specified class and with proper vertex position) that in addition fulfils the trigger condition [to disentangle trigger eff from reco eff ]
      histos.add((dir + "multCorrel_evt").c_str(), "", kTHnSparseD, {multMeasAxis, multTrueAxis});      // multiplicity correlation of measured events (excluding background events)
      histos.add((dir + "multCorrel_prim").c_str(), "", kTHnSparseD, {multMeasAxis, multTrueAxis});     // multiplicity correlation of measured primary charged particles (excluding particles from background events)
      histos.add((dir + "ptCorrel_prim").c_str(), "", kTH2D, {ptMeasAxis, ptTrueAxis});                 // pT correlation of measured primary charged particles  (excluding particles from background events)
      histos.add((dir + "multPtSpec_prim_gen").c_str(), "", kTH2D, {multTrueAxis, ptTrueAxis});         // generated primary charged particles as function of true properties (from events within specified class and with proper vertex position)
      histos.add((dir + "multPtSpec_prim_gen_evtloss").c_str(), "", kTH2D, {multTrueAxis, ptTrueAxis}); // generated primary charged particles of events that did not pass the event selection as function of multiplicity and pt
      histos.add((dir + "multPtSpec_prim_meas").c_str(), "", kTH2D, {multTrueAxis, ptTrueAxis});        // measured primary charged particles as function of true properties (no contamination from background events)
      histos.add((dir + "multPtSpec_trk_prim_meas").c_str(), "", kTH2D, {multMeasAxis, ptMeasAxis});    // tracks from measured primaries (no contamination from secondaries, particles smeared into acceptance or background events)
      histos.add((dir + "multPtSpec_trk_sec_meas").c_str(), "", kTH2D, {multMeasAxis, ptMeasAxis});     // tracks from measured secondaries (no contamination from particles smeared into acceptance or background events)  [for QA to disentangle secondaries from other contamination]
      histos.add((dir + "multPtSpec_trk_pileup_meas").c_str(), "", kTH2D, {multMeasAxis, ptMeasAxis});  // tracks from other mc collisions (pileup) in measured events of the class  [for QA to disentangle pileup from other contamination]
      histos.add((dir + "multPtSpec_trk_meas_evtcont").c_str(), "", kTH2D, {multMeasAxis, ptMeasAxis}); // tracks from events that are measured, but do not belong to the desired class of events
      histos.add((dir + "multPtSpec_trk_inter").c_str(), "", kTH2D, {multTrueAxis, ptMeasAxis});
      if (isNominal) {
        histos.add((dir + "QA/events/multTrueBits").c_str(), "", kTH2D, {multTrueAxis, bitAxis});
        histos.add((dir + "QA/events/eventSelectionMC").c_str(), "", kTH1D, {{NStagesMC, -0.5, NStagesMC - 0.5}});
        histos.add((dir + "QA/tracks/deltaPt").c_str(), "", kTH2D, {ptMeasAxis, {relResEdges, "#Delta(#it{p}_{T}) / #it{p}^{ meas}_{T}"}});
        // event and tracking efficiency per run, as in the unfolding (run axis labelled as runStatistics)
        histos.add((dir + "QA/events/evtGenRun").c_str(), "", kTH2D, {{512, 0., 512.}, runMultTrueAxis});
        histos.add((dir + "QA/events/evtRecoRun").c_str(), "", kTH2D, {{512, 0., 512.}, runMultTrueAxis});
        histos.add((dir + "QA/tracks/primGenRun").c_str(), "", kTH2D, {{512, 0., 512.}, ptTrueAxis});
        histos.add((dir + "QA/tracks/primRecoRun").c_str(), "", kTH2D, {{512, 0., 512.}, ptTrueAxis});
      }
    }

    histos.add("multPtSpec_prim_gen_notrig", "", kTH2D, {multTrueAxis, ptTrueAxis}); // generated primary charged particles of events with proper vertex position that did not fulfil the trigger condition (tvx; run 2: v0and, 2010/2011: int1), independent of the event class, as function of multiplicity and pt
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
    for (int i = 0; i < NLayersITS; ++i) {
      axisHits->SetBinLabel(i + 1, ("L" + std::to_string(i)).c_str());
    }
    auto* axisMatch = histos.get<TH2>(HIST("QA/tracks/itsTpcMatch"))->GetYaxis();
    axisMatch->SetBinLabel(1, "tpc cuts");
    axisMatch->SetBinLabel(2, "matchable (its or x < 90 cm)");
    axisMatch->SetBinLabel(3, "pointing in z");
    axisMatch->SetBinLabel(4, "with its track");
    axisMatch->SetBinLabel(5, "its cuts");
    metaLabel = metadataInfo.isInitialized() ? metadataInfo.makeMetadataLabel() : std::string("no ao2d metadata");
    const char* o2PhysicsRoot = std::getenv("O2PHYSICS_ROOT"); // grid: .../Packages/O2Physics/<tag>
    metaLabel += ", O2Physics " + (o2PhysicsRoot ? std::filesystem::path(o2PhysicsRoot).filename().string() : std::string("unknown"));
    const char* train = std::getenv("ALIEN_JDL_HYTRAIN"); // hyperloop train jobs
    const auto session = ic.services().get<RawDeviceService>().device()->GetConfig()->GetPropertyAsString("session", "");
    if (train) {
      metaLabel += ", train " + std::string(train);
    } else if (!session.empty() && !session.starts_with("dpl_")) {
      metaLabel += ", session " + session; // hyperloop tests: test number (dpl_<n>: random local default)
    }
    if (ptRange.value != "nominal") {
      metaLabel += ", ptRange " + ptRange.value + ", fiducialMult " + (fiducialMult.value ? "true" : "false");
    }
    if (multEstimator.value != 0) {
      metaLabel += ", multEstimator " + std::to_string(multEstimator.value);
    }
    if (doprocessMC || doprocessMCnoPCC) {
      metaLabel += ", " + eventClasses.value;
      forEachDir([&](auto dir) {
        constexpr int Dir = decltype(dir)::value;
        setBitLabels(histos.get<TH2>(HIST(ClassDir[Dir]) + HIST("QA/events/multTrueBits"))->GetYaxis());
        auto* axisMC = histos.get<TH1>(HIST(ClassDir[Dir]) + HIST("QA/events/eventSelectionMC"))->GetXaxis();
        int bin = 1;
        for (const auto& label : {"all", "|z_{true}| < maxVtxZ", "event class", "tvx", "reco collision", "sel8", "evSelBits", "|z| < maxVtxZ"}) {
          axisMC->SetBinLabel(bin++, label);
        }
      });
    }
    histos.get<TH1>(HIST("QA/events/metaInfo"))->Fill(metaLabel.c_str(), 1.); // once per task instance: number of jobs after merging
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
  } else if (systMode == kSystDownTPCFracSharedCls) {
    trackSel.SetMaxTPCFractionSharedCls(0.4f); // not in the nominal mask: nothing to remove from it
    trackSelFlag = TrackSelection::TrackCuts::kTPCFracSharedCls;
  } else if (systMode == kSystDownDCAxy || systMode == kSystUpDCAxy || systMode == kSystDownDCAz || systMode == kSystUpDCAz) {
    // own dca cuts with another number of sigmas (below)
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
  } else if (systMode == kSystDownPCC || systMode == kSystUpPCC) {
    if (!doprocessMC) {
      LOGP(fatal, "systMode {} (pcc variation) needs processMC (with pcc)", systMode.value);
    }
  } else if (systMode != kSystNominal && systMode != kSystDownVtxZ && systMode != kSystEvSelNoCollInRof && systMode != kSystEvSelNoCollInTimeRange) {
    LOGP(fatal, "undefined systMode {} (valid: 100-107, 111-117, 120-121, 130-132)", systMode.value);
  }
  isCutVar = (cutVarFlag != 0) || (systMode == kSystDownTPCFracSharedCls);
  trackSelMask &= (~cutVarFlag);

  // own dca cuts replace the dca bits of the central track selection
  trackSelMask &= ~(TrackSelectionFlags::kDCAxy | TrackSelectionFlags::kDCAz);
  float nSigmaXY = NSigmaDcaXY;
  float nSigmaZ = NSigmaDcaZ;
  if (systMode == kSystDownDCAxy) {
    nSigmaXY = 4.f;
  } else if (systMode == kSystUpDCAxy) {
    nSigmaXY = 10.f;
  } else if (systMode == kSystDownDCAz) {
    nSigmaZ = 6.f;
  } else if (systMode == kSystUpDCAz) {
    nSigmaZ = 14.f;
  }
  const DcaPars resXY = isRun3 ? DcaResRun3XY : DefaultDcaXY;
  const DcaPars resZ = isRun3 ? DcaResRun3Z : DefaultDcaZ;
  dcaCutXY = {.p0 = nSigmaXY * resXY.p0, .p1 = nSigmaXY * resXY.p1, .p2 = resXY.p2};
  dcaCutZ = {.p0 = nSigmaZ * resZ.p0, .p1 = nSigmaZ * resZ.p1, .p2 = resZ.p2};

  if (multEstimator.value != 0 && multEstimator.value != 1) {
    LOGP(fatal, "invalid multEstimator {} (0 = selected tracks, 1 = also its-only tracks)", multEstimator.value);
  }
  if (multEstimator.value == 1 && !isRun3) {
    LOGP(fatal, "multEstimator 1 (its-only tracks) is defined for run 3 only");
  }
  if (rateRun != nullptr) { // ctp scalers: run 3 only
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false); // missing objects: run without rate logging (hasCtpRateInfo)
  }
}

//**************************************************************************************************
/**
 * Entrypoint to processes data.
 */
//**************************************************************************************************
void ChargedParticles::processData(CollisionTableData::iterator const& collision, TrackTableData const& tracks, aod::FT0s const& ft0s, aod::BCsWithTimestamps const&)
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
  tracksPerParticle.clear();

  initEventMC(mcCollision, particles);
  if (isNominal && isRun3) {
    mcRunLabel = std::to_string(mcCollision.template bc_as<BCsWithSel>().runNumber()); // also for events without reco collision
  }

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
    for (const auto& [id, entry] : tracksPerParticle) {
      histos.fill(HIST("QA/tracks/tracksPerParticle"), entry.pt, entry.nTracks);
      histos.fill(HIST("QA/tracks/collisionsPerParticle"), entry.pt, entry.nCollisions);
    }
  }

  if (isNominal) {
    histos.fill(HIST("QA/events/zVtxGen"), mcCollision.posZ());
    // cumulative mc event flow
    forEachDir([&](auto dir) {
      constexpr int Dir = decltype(dir)::value;
      const std::array<bool, NStagesMC> stages = {true, vars.isGoodZPosMC, isInDir(Dir), vars.isTriggeredMC, collisions.size() > 0,
                                                  isSel8Representative, hasBitsRepresentative, isAcceptedRepresentative};
      for (int i = 0; i < NStagesMC && stages[i]; ++i) {
        histos.fill(HIST(ClassDir[Dir]) + HIST("QA/events/eventSelectionMC"), i);
      }
    });
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

  // run 2 mc with pileup: particles of the embedded pileup events are not part of this event
  // (in run 3 the flag marks the main event itself, so it is not used there)
  if (!isRun3 && particle.fromBackgroundEvent()) {
    vars.isChargedPrimary = false;
    return false;
  }

  if (std::abs(particle.eta()) >= etaCut) {
    return false;
  }
  if (particle.pt() <= ptMin || particle.pt() >= ptMax) {
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
 * ITS-only track (run 3): no TPC, passes every selection that does not concern the TPC (with the varied cut of a
 * systematic variant, if not a TPC cut) and the optional dca_z cut, in the eta/pt and multiplicity window.
 */
//**************************************************************************************************
template <typename T>
bool ChargedParticles::isItsOnlyCandidate(const T& track)
{
  if (track.hasTPC() || std::abs(track.eta()) >= etaCut || track.pt() <= ptMin || track.pt() >= ptMax || !isInMultWindow(track.pt())) {
    return false;
  }
  constexpr uint32_t TpcBits = TrackSelectionFlags::kTPCNCls | TrackSelectionFlags::kTPCCrossedRows | TrackSelectionFlags::kTPCCrossedRowsOverNCls | TrackSelectionFlags::kTPCChi2NDF | TrackSelectionFlags::kTPCRefit;
  if (((~static_cast<uint32_t>(track.trackCutFlag())) & trackSelMask & ~TpcBits) != 0) {
    return false;
  }
  if (!passesDcaXY(track) || !passesDcaZ(track)) {
    return false;
  }
  return !(isCutVar && !isTpcCut(trackSelFlag) && !trackSel.IsSelected(track, trackSelFlag));
}

//**************************************************************************************************
/**
 * Own DCA cuts |DCA| < k * sigma(pT) (dcaCutXY, dcaCutZ = k * parameters; same arithmetic as the central selection).
 */
//**************************************************************************************************
template <typename T>
bool ChargedParticles::passesDcaXY(const T& track) const
{
  return std::fabs(track.dcaXY()) <= dcaCutXY.p0 + dcaCutXY.p1 / std::pow(track.pt(), dcaCutXY.p2);
}

template <typename T>
bool ChargedParticles::passesDcaZ(const T& track) const
{
  return std::fabs(track.dcaZ()) <= dcaCutZ.p0 + dcaCutZ.p1 / std::pow(track.pt(), dcaCutZ.p2);
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
  if (track.pt() <= ptMin || track.pt() >= ptMax) {
    return false;
  }
  if (!TrackSelectionFlags::checkFlag(track.trackCutFlag(), trackSelMask)) {
    return false;
  }
  if (!passesDcaXY(track) || !passesDcaZ(track)) {
    return false;
  }
  // for systematic variation of standard selections, check if the varied cut is passed
  if (isCutVar && !trackSel.IsSelected(track, trackSelFlag)) {
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
    // same definition in data and mc: every selected track counts
    // multEstimator 1: also its-only candidates
    if ((initTrack(track) && isInMultWindow(track.pt())) || (multEstimator.value == 1 && isItsOnlyCandidate(track))) {
      vars.multMeas += getTrackRepetitions<IS_MC, P>(track, isValidParticle);
    }
  }

  vars.isAcceptedEvent = false;
  vars.zVtx = collision.posZ();
  if (std::abs(collision.posZ()) < vtxZCut()) {
    if (isRun3) {
      if (collision.sel8() && (collision.selection_raw() & evSelMask) == evSelMask) {
        vars.isAcceptedEvent = true;
      }
    } else {
      int runNumber = -1;
      if constexpr (IS_MC) {
        runNumber = collision.template bc_as<BCsWithSel>().runNumber();
      } else {
        runNumber = collision.template bc_as<aod::BCsWithTimestamps>().runNumber();
      }
      const bool isMC = doprocessMC || doprocessMCnoPCC;
      if (isINT1Period(runNumber)) {
        // int1 selection of the event selection (sel1, not stored there)
        const uint64_t sel = collision.selection_raw();
        const bool sel1 = std::ranges::all_of(std::array{aod::evsel::kIsINT1, aod::evsel::kNoBGV0A, aod::evsel::kNoBGV0C, aod::evsel::kNoTPCLaserWarmUp, aod::evsel::kNoTPCHVdip},
                                              [sel](int bit) { return ((sel >> bit) & 1ULL) != 0; });
        if (sel1 && (isMC || collision.alias_bit(kINT1))) {
          vars.isAcceptedEvent = true;
        }
      } else if (collision.sel7() && (isMC || collision.alias_bit(kINT7))) {
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
    if (isInAcceptance && isInMultWindow(particle.pt())) {
      vars.multTrue += vars.nRepetitions;
    }
  }
  vars.isGoodZPosMC = std::abs(collision.posZ()) < vtxZCut();
  // trigger of the mc collision's bunch crossing: tvx (run 3), v0and as for kINT7 (run 2), int1 (2010/2011)
  const auto& bc = collision.template bc_as<BCsWithSel>();
  if (isRun3) {
    vars.isTriggeredMC = bc.selection_bit(aod::evsel::kIsTriggerTVX);
  } else if (isINT1Period(bc.runNumber())) {
    vars.isTriggeredMC = bc.selection_bit(aod::evsel::kIsINT1);
  } else {
    vars.isTriggeredMC = bc.selection_bit(aod::evsel::kIsBBV0A) && bc.selection_bit(aod::evsel::kIsBBV0C);
  }
  vars.classMask = 0;
  if (vars.isGoodZPosMC) {
    vars.classMask |= (vars.multTrue > 0) ? (1 << kFiducial) : 0;
    vars.classMask |= vars.isINELGT0EventMC ? (1 << kINELGT0) : 0;
    vars.classMask |= vars.isTriggeredMC ? (1 << kTriggered) : 0;
    vars.classMask |= (1 << kINEL); // all generated events
  }
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
  if (std::abs(track.eta()) >= etaCut || track.pt() <= ptMin || track.pt() >= ptMax) {
    return;
  }
  const auto flags = track.trackCutFlag();
  const int nCuts = trackCutBits.size();
  const int dcaXYBit = static_cast<int>(TrackSelection::TrackCuts::kDCAxy);
  const int dcaZBit = static_cast<int>(TrackSelection::TrackCuts::kDCAz);
  const bool ownDcaXY = passesDcaXY(track);
  const bool ownDcaZ = passesDcaZ(track);
  const auto passes = [&](int bit) { // dca bits: own cuts
    return bit == dcaXYBit ? ownDcaXY : (bit == dcaZBit ? ownDcaZ : (flags & (1 << bit)) != 0);
  };
  int nFailed = 0;
  int failedBit = -1;
  for (const auto& bit : trackCutBits) {
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
  for (const auto& bit : trackCutBits) {
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
  // its-tpc matching; tpc-only tracks starting deep in the tpc (second arcs, loopers) are not matchable;
  // pointing in z: straight-line extrapolation in r-z to the beam line (= dca_z for propagated tracks) vs the vertex
  constexpr float MaxTpcStartX = 90.f; // cm, tpc inner wall at about 83 cm
  constexpr float MaxDzPointing = 3.f; // cm
  if ((flags & tpcCutMask) == tpcCutMask) {
    histos.fill(HIST("QA/tracks/itsTpcMatch"), pt, 0);
    if (track.hasITS() || track.x() < MaxTpcStartX) {
      histos.fill(HIST("QA/tracks/itsTpcMatch"), pt, 1);
      if (std::abs(track.z() - track.tgl() * track.x() - vars.zVtx) < MaxDzPointing) {
        histos.fill(HIST("QA/tracks/itsTpcMatch"), pt, 2);
        if (track.hasITS()) {
          histos.fill(HIST("QA/tracks/itsTpcMatch"), pt, 3);
          if ((flags & itsCutMask) == itsCutMask) {
            histos.fill(HIST("QA/tracks/itsTpcMatch"), pt, 4);
          }
        }
      }
    }
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
  histos.fill(HIST("QA/tracks/etaPhi"), track.eta(), track.phi());
  if (phiRunBin > 0) {
    histos.fill(HIST("QA/tracks/phiRun"), track.phi(), phiRunBin - 0.5);
  }
  histos.fill(HIST("QA/tracks/itsNCls"), track.itsNCls());
  for (int i = 0; i < NLayersITS; ++i) {
    if (track.itsClusterMap() & (1 << i)) {
      histos.fill(HIST("QA/tracks/itsHits"), i);
    }
  }
  histos.fill(HIST("QA/tracks/itsChi2NCl"), track.itsChi2NCl());
  histos.fill(HIST("QA/tracks/tpcNClsFindable"), track.tpcNClsFindable());
  histos.fill(HIST("QA/tracks/tpcNClsFound"), track.tpcNClsFound());
  histos.fill(HIST("QA/tracks/tpcCrossedRows"), track.tpcNClsCrossedRows());
  histos.fill(HIST("QA/tracks/tpcCrossedRowsOverFindableCls"), track.tpcCrossedRowsOverFindableCls());
  histos.fill(HIST("QA/tracks/tpcFractionSharedCls"), track.pt(), track.tpcFractionSharedCls());
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
/**
 * Whether the ccdb holds the objects the ctp rate fetcher needs for this run (it stops the task when one is missing).
 */
//**************************************************************************************************
bool ChargedParticles::hasCtpRateInfo(uint64_t timestamp, int runNumber)
{
  std::map<std::string, std::string> metadata;
  const bool hasLhcIf = ccdb->getSpecific<o2::parameters::GRPLHCIFData>("GLO/Config/GRPLHCIF", timestamp, metadata) != nullptr;
  metadata["runNumber"] = std::to_string(runNumber);
  const bool hasConfig = ccdb->getSpecific<o2::ctp::CTPConfiguration>("CTP/Config/Config", timestamp, metadata) != nullptr;
  const bool hasScalers = ccdb->getSpecific<o2::ctp::CTPRunScalers>("CTP/Calib/Scalers", timestamp, metadata) != nullptr;
  if (!(hasLhcIf && hasConfig && hasScalers)) {
    LOGP(warning, "run {}: no ctp rate information in the ccdb (lhcif {}, config {}, scalers {}), no interaction rate logged", runNumber, hasLhcIf, hasConfig, hasScalers);
    return false;
  }
  return true;
}

//**************************************************************************************************
template <bool IS_MC, typename C>
void ChargedParticles::fillEvSelQA(const C& collision, aod::FT0s const& ft0s, bool isRepresentative)
{
  if (!isNominal || !isRun3) {
    return;
  }
  int evSelRunNumber = -1;
  if constexpr (IS_MC) {
    evSelRunNumber = collision.template bc_as<BCsWithSel>().runNumber();
  } else {
    evSelRunNumber = collision.template bc_as<aod::BCsWithTimestamps>().runNumber();
  }
  if (evSelRunNumber != lastEvSelRunNumber) {
    lastEvSelRunNumber = evSelRunNumber;
    evSelRunLabel = std::to_string(evSelRunNumber);
  }
  auto* evSelRun = histos.get<TH2>(HIST("QA/events/eventSelectionRun")).get(); // filled by run label, as runStatistics
  histos.fill(HIST("QA/events/eventSelection"), 0);
  evSelRun->Fill(evSelRunLabel.c_str(), 0., 1.);
  if (!collision.sel8()) {
    return;
  }
  // cumulative cut-flow of the applied selection
  int stage = 1;
  histos.fill(HIST("QA/events/eventSelection"), stage);
  evSelRun->Fill(evSelRunLabel.c_str(), stage, 1.);
  bool passes = true;
  for (const auto& bit : evSelBitList) {
    passes = passes && ((collision.selection_raw() >> bit) & 1ULL);
    ++stage;
    if (passes) {
      histos.fill(HIST("QA/events/eventSelection"), stage);
      evSelRun->Fill(evSelRunLabel.c_str(), stage, 1.);
    }
  }
  if (passes && std::abs(collision.posZ()) < vtxZCut()) {
    histos.fill(HIST("QA/events/eventSelection"), stage + 1);
    evSelRun->Fill(evSelRunLabel.c_str(), stage + 1, 1.);
  }
  histos.fill(HIST("QA/events/zVtxMeas"), collision.posZ());
  if (vars.isAcceptedEvent) {
    int runNumber = -1;
    uint64_t timestamp = 0;
    if constexpr (IS_MC) {
      const auto& bc = collision.template bc_as<BCsWithSel>();
      runNumber = bc.runNumber();
      timestamp = bc.timestamp();
    } else {
      const auto& bc = collision.template bc_as<aod::BCsWithTimestamps>();
      runNumber = bc.runNumber();
      timestamp = bc.timestamp();
    }
    if (runNumber != lastRunNumber) {
      lastRunNumber = runNumber;
      lastRunLabel = std::to_string(runNumber);
      phiRunBin = histos.get<TH2>(HIST("QA/tracks/phiRun"))->GetYaxis()->FindBin(lastRunLabel.c_str()); // adds the label once per run
      hasRateInfo = (rateRun != nullptr) && hasCtpRateInfo(timestamp, runNumber);
    }
    histos.get<TH1>(HIST("QA/events/runStatistics"))->Fill(lastRunLabel.c_str(), 1.);
    if (hasRateInfo) {
      const double rate = rateFetcher.fetch(ccdb.service, timestamp, runNumber, irSource, false); // hz; negative: not available
      if (rate >= 0.) {
        rateRun->Fill(lastRunLabel.c_str(), rate * 1.e-3); // khz
      }
    }
  }
  if (std::abs(collision.posZ()) >= vtxZCut()) {
    return;
  }
  const uint64_t selection = collision.selection_raw();
  auto fillBit = [&](int bin) {
    histos.fill(HIST("QA/events/multMeasBits"), vars.multMeas, bin);
    if constexpr (IS_MC) {
      if (isRepresentative) {
        forEachDir([&](auto dir) {
          constexpr int Dir = decltype(dir)::value;
          if (isInDir(Dir)) {
            histos.fill(HIST(ClassDir[Dir]) + HIST("QA/events/multTrueBits"), vars.multTrue, bin);
          }
        });
      }
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
  if (vars.classMask == 0 && !isUntriggeredGoodZPos) {
    return;
  }
  const bool fillRunQA = isNominal && isRun3; // per-run efficiencies

  forEachDir([&](auto dir) {
    constexpr int Dir = decltype(dir)::value;
    if (isInDir(Dir)) {
      histos.fill(HIST(ClassDir[Dir]) + HIST("multDist_evt_gen"), vars.multTrue);
      if (fillRunQA) {
        histos.get<TH2>(HIST(ClassDir[Dir]) + HIST("QA/events/evtGenRun"))->Fill(mcRunLabel.c_str(), vars.multTrue, 1.);
        auto* axis = histos.get<TH2>(HIST(ClassDir[Dir]) + HIST("QA/tracks/primGenRun"))->GetXaxis();
        primGenRunX[Dir] = axis->GetBinCenter(axis->FindBin(mcRunLabel.c_str())); // adds the label once per run
      }
      if (vars.isTriggeredMC) {
        histos.fill(HIST(ClassDir[Dir]) + HIST("multDist_evt_gen_trig"), vars.multTrue);
      }
    }
  });

  for (const auto& particle : particles) {
    if (initParticle(particle) && vars.isChargedPrimary) {
      for (auto i = 0u; i < vars.nRepetitions; ++i) {
        forEachDir([&](auto dir) {
          constexpr int Dir = decltype(dir)::value;
          if (isInDir(Dir)) {
            histos.fill(HIST(ClassDir[Dir]) + HIST("multPtSpec_prim_gen"), vars.multTrue, particle.pt());
            if (!vars.isAcceptedEvent) {
              histos.fill(HIST(ClassDir[Dir]) + HIST("multPtSpec_prim_gen_evtloss"), vars.multTrue, particle.pt());
            } else if (fillRunQA) {
              histos.fill(HIST(ClassDir[Dir]) + HIST("QA/tracks/primGenRun"), primGenRunX[Dir], particle.pt());
            }
          }
        });
        if (isUntriggeredGoodZPos) {
          histos.fill(HIST("multPtSpec_prim_gen_notrig"), vars.multTrue, particle.pt());
        }
      }
    }
  }
}

//**************************************************************************************************
/**
 * ITS-TPC matching QA: in collisions with 2-4 candidates (selected tracks and ITS-only tracks), for every selected
 * track the other candidates, matched or not.
 */
//**************************************************************************************************
template <typename T>
void ChargedParticles::fillMatchQA(const T& tracks)
{
  constexpr size_t MinCandidates = 2;
  constexpr size_t MaxCandidates = 4;
  std::array<float, MaxCandidates> selectedPt{}, itsOnlyPt{};
  size_t nSelected = 0, nItsOnly = 0;
  for (const auto& track : tracks) {
    const bool isSelected = initTrack(track) && isInMultWindow(track.pt());
    if (!isSelected && !isItsOnlyCandidate(track)) {
      continue;
    }
    if (nSelected + nItsOnly == MaxCandidates) {
      return; // more than 4 candidates
    }
    if (isSelected) {
      selectedPt[nSelected++] = track.pt();
    } else {
      itsOnlyPt[nItsOnly++] = track.pt();
    }
  }
  const size_t nCandidates = nSelected + nItsOnly;
  if (nCandidates < MinCandidates) {
    return;
  }
  for (size_t i = 0; i < nSelected; ++i) {
    for (size_t j = 0; j < nSelected; ++j) {
      if (j != i) {
        histos.fill(HIST("QA/tracks/matchOthersSelected"), selectedPt[i], nCandidates, selectedPt[j]);
      }
    }
    for (size_t j = 0; j < nItsOnly; ++j) {
      histos.fill(HIST("QA/tracks/matchOthersItsOnly"), selectedPt[i], nCandidates, itsOnlyPt[j]);
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
  if (isNominal && isRun3) { // its-only = no tpc: run 3 tracking
    fillMatchQA(tracks);
  }
  // mc: only the representative collision corresponds to the true event (per event class)
  [[maybe_unused]] auto isTrueEvent = [&](int dir) { return isRepresentative && isInDir(dir); };

  histos.fill(HIST("multDist_evt_meas"), vars.multMeas);

  if constexpr (IS_MC) {
    forEachDir([&](auto dir) {
      constexpr int Dir = decltype(dir)::value;
      if (isTrueEvent(Dir)) {
        histos.fill(HIST(ClassDir[Dir]) + HIST("multCorrel_evt"), vars.multMeas, vars.multTrue);
        if (isNominal && isRun3) {
          histos.get<TH2>(HIST(ClassDir[Dir]) + HIST("QA/events/evtRecoRun"))->Fill(mcRunLabel.c_str(), vars.multTrue, 1.);
          auto* axis = histos.get<TH2>(HIST(ClassDir[Dir]) + HIST("QA/tracks/primRecoRun"))->GetXaxis();
          primRecoRunX[Dir] = axis->GetBinCenter(axis->FindBin(mcRunLabel.c_str())); // adds the label once per run
        }
      }
    });
  }

  std::vector<int> foundParticles;
  for (const auto& track : tracks) {
    if (isNominal) { // qa with the pcc repetitions (= 1 without pcc)
      bool isValidCandidate{false};
      const uint32_t nRepCandidate = getTrackRepetitions<IS_MC, P>(track, isValidCandidate);
      for (auto i = 0u; i < nRepCandidate; ++i) {
        fillTrackCutQA(track);
      }
    }
    if (!initTrack(track)) {
      continue;
    }

    // all fills repeated as often as the particle (pcc)
    bool isValidParticle{false};
    const uint32_t nRepetitions = getTrackRepetitions<IS_MC, P>(track, isValidParticle);
    [[maybe_unused]] bool isPileup{false}; // particle from another mc collision
    if constexpr (IS_MC) {
      isPileup = track.has_mcParticle() && track.template mcParticle_as<P>().mcCollisionId() != vars.mcCollisionIndex;
    }
    if (isNominal) {
      for (auto i = 0u; i < nRepetitions; ++i) {
        fillTrackQA(track);
      }
      if constexpr (IS_MC) {
        if (isValidParticle && vars.isChargedPrimary) {
          const float relDeltaPt = std::abs(track.pt() - track.template mcParticle_as<P>().pt()) / track.pt();
          forEachDir([&](auto dir) {
            constexpr int Dir = decltype(dir)::value;
            if (isTrueEvent(Dir)) {
              for (auto i = 0u; i < nRepetitions; ++i) {
                histos.fill(HIST(ClassDir[Dir]) + HIST("QA/tracks/deltaPt"), track.pt(), relDeltaPt);
              }
            }
          });
        }
      }
    }
    if constexpr (IS_MC) {
      if (track.has_mcParticle()) {
        foundParticles.push_back(track.mcParticleId());
        if (isNominal) { // tracks and collisions per particle over all reco collisions of the mc event
          auto& entry = tracksPerParticle[track.mcParticleId()];
          if (entry.nTracks == 0) {
            entry.pt = track.template mcParticle_as<P>().pt();
          }
          ++entry.nTracks;
          if (entry.lastCollision != track.collisionId()) {
            ++entry.nCollisions;
            entry.lastCollision = track.collisionId();
          }
        }
      }
    }

    for (auto i = 0u; i < nRepetitions; ++i) {
      histos.fill(HIST("multPtSpec_trk_meas"), vars.multMeas, track.pt());

      if constexpr (IS_MC) {
        forEachDir([&](auto dir) {
          constexpr int Dir = decltype(dir)::value;
          if (!isTrueEvent(Dir)) {
            histos.fill(HIST(ClassDir[Dir]) + HIST("multPtSpec_trk_meas_evtcont"), vars.multMeas, track.pt());
            return;
          }
          histos.fill(HIST(ClassDir[Dir]) + HIST("multPtSpec_trk_inter"), vars.multTrue, track.pt());
          if (isValidParticle) {
            if (!vars.isChargedPrimary) {
              histos.fill(HIST(ClassDir[Dir]) + HIST("multPtSpec_trk_sec_meas"), vars.multMeas, track.pt());
            } else {
              const auto& particle = track.template mcParticle_as<P>();
              histos.fill(HIST(ClassDir[Dir]) + HIST("multCorrel_prim"), vars.multMeas, vars.multTrue);
              histos.fill(HIST(ClassDir[Dir]) + HIST("ptCorrel_prim"), track.pt(), particle.pt());
              histos.fill(HIST(ClassDir[Dir]) + HIST("multPtSpec_prim_meas"), vars.multTrue, particle.pt());
              if (isNominal && isRun3) {
                histos.fill(HIST(ClassDir[Dir]) + HIST("QA/tracks/primRecoRun"), primRecoRunX[Dir], particle.pt());
              }
              histos.fill(HIST(ClassDir[Dir]) + HIST("multPtSpec_trk_prim_meas"), vars.multMeas, track.pt());
            }
          } else if (isPileup) {
            histos.fill(HIST(ClassDir[Dir]) + HIST("multPtSpec_trk_pileup_meas"), vars.multMeas, track.pt());
          }
        });
      }
    }
  }

  std::unordered_map<int32_t, int> nTracksPerParticle;
  for (const auto& mcParticleID : foundParticles) {
    ++nTracksPerParticle[mcParticleID];
  }
  for (const auto& entry : nTracksPerParticle) {
    histos.fill(HIST("track_ambiguity"), entry.second);
  }
}
