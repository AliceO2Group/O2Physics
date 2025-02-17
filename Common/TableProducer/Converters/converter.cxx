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
/// \file   converter.cxx
/// \since  2024-11-12
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \brief  Combined converter task
///

#include <vector>
#include <string>

// O2 includes
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"

// O2Physics includes
#include "TableHelper.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::zdc;

template <typename FlagType>
void autoSetProcessFunction(o2::framework::InitContext& initContext, const std::string& table, FlagType& flag)
{
  flag.value = false;
  enableFlagIfTableRequired(initContext, table, flag);
}

// Converts bc_000 into bc_001
struct converterBCs_000 {
  Produces<aod::BCs_001> BCs_001;
  void process(aod::BCs_000 const& bcTable) // BC converter is always needed
  {
    for (const auto& bc : bcTable) {
      constexpr uint64_t lEmptyTriggerInputs = 0;
      BCs_001(bc.runNumber(), bc.globalBC(), bc.triggerMask(), lEmptyTriggerInputs);
    }
  }
};

// Swaps covariance matrix elements if the data is known to be bogus (collision_000 is bogus)
struct converterCollisions_000 {
  Produces<aod::Collisions_001> Collisions_001;
  Configurable<bool> doNotSwap{"doNotSwap", false, "simple pass-through"};
  Configurable<bool> debug{"debug", false, "flag to save debug histo"};
  Configurable<int> nbins{"nbins", 1, "number of bins in debug histo"};
  Configurable<float> tolerance{"tolerance", 1e-3, "Tolerance for CYY check"};
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  void init(o2::framework::InitContext& initContext)
  {
    autoSetProcessFunction(initContext, "Collisions_001", doprocessConverter);
    const AxisSpec axisCYYdebug{nbins, -1.0f, +1.0f, ""};
    histos.add("hCYY", "hCYY", kTH1F, {axisCYYdebug});
  }
  void process(aod::BCs const&) {} // Dummy processor in case the other is disabled
  void processConverter(aod::Collisions_000 const& collisionTable)
  {
    float negtolerance = -1.0f * tolerance;
    for (auto& collision : collisionTable) {
      float lYY = collision.covXZ();
      float lXZ = collision.covYY();
      if (doNotSwap) {
        lYY = collision.covYY();
        lXZ = collision.covXZ();
      }
      if (debug)
        histos.fill(HIST("hCYY"), lYY);
      if (lYY < negtolerance) {
        // This happened by accident!
        if (!doNotSwap && !debug) {
          LOGF(info, "Collision converter task found negative YY element!");
          LOGF(info, "CYY = %.10f, exceeds tolerance of %.10f", lYY, negtolerance);
          LOGF(info, "This is an indication that you're looping over data");
          LOGF(info, "produced with an O2 version of late December 2022.");
          LOGF(info, "Unfortunately, O2 versions of late December 2022");
          LOGF(info, "have a mistake in them for which a special mode");
          LOGF(info, "of this task exists. ");
          LOGF(info, "For this data, please operate the collision converter");
          LOGF(info, "with the configurable 'doNotSwap' set to true.");
          LOGF(info, "This program will now crash. Please adjust your settings!");
          LOGF(fatal, "FATAL: please set doNotSwap to true!");
        }
        if (!debug) {
          LOGF(info, "Collision converter task found negative YY element!");
          LOGF(info, "CYY = %.10f, exceeds tolerance of %.10f", lYY, negtolerance);
          LOGF(info, "You're running with 'doNotSwap' enabled, but the ");
          LOGF(info, "data your're analysing requires it to be disabled. ");
          LOGF(info, "This program will now crash. Please adjust your settings!");
          LOGF(fatal, "FATAL: please set doNotSwap to false!");
        }
      }
      // Repopulate new table
      Collisions_001(
        collision.bcId(),
        collision.posX(), collision.posY(), collision.posZ(),
        collision.covXX(),
        collision.covXY(),
        lYY, // <- this is the fixed part
        lXZ, // <- this is the fixed part
        collision.covYZ(),
        collision.covZZ(),
        collision.flags(), collision.chi2(), collision.numContrib(),
        collision.collisionTime(), collision.collisionTimeRes());
    }
  }
  PROCESS_SWITCH(converterCollisions_000, processConverter, "Process converter (autoset)", true);
};

// Converts FDD table from version 000 to 001
struct converterFDDs_000 {
  Produces<aod::FDDs_001> FDDs_001;
  void init(o2::framework::InitContext& initContext) { autoSetProcessFunction(initContext, "FDDs_001", doprocessConverter); }
  void process(aod::BCs const&) {} // Dummy processor in case the other is disabled
  void processConverter(aod::FDDs_000 const& fdd_000)
  {
    for (auto& p : fdd_000) {
      int16_t chargeA[8] = {0u};
      int16_t chargeC[8] = {0u};

      for (int i = 0; i < 4; i++) {
        chargeA[i] = p.amplitudeA()[i];
        chargeA[i + 4] = p.amplitudeA()[i];

        chargeC[i] = p.amplitudeC()[i];
        chargeC[i + 4] = p.amplitudeC()[i];
      }

      FDDs_001(p.bcId(), chargeA, chargeC,
               p.timeA(), p.timeC(), p.triggerMask());
    }
  }
  PROCESS_SWITCH(converterFDDs_000, processConverter, "Process converter (autoset)", true);
};

struct converterHMPID_000 {
  Produces<aod::HMPID_001> HMPID_001;
  void init(o2::framework::InitContext& initContext) { autoSetProcessFunction(initContext, "HMPID_001", doprocessConverter); }
  void process(aod::BCs const&) {} // Dummy processor in case the other is disabled
  void processConverter(aod::HMPID_000 const& hmpLegacy, aod::Tracks const&)
  {
    for (auto& hmpData : hmpLegacy) {

      float phots[] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
      auto trackid = hmpData.trackId();
      auto hmpidSignal = hmpData.hmpidSignal();
      auto hmpidXTrack = -999.; // dummy
      auto hmpidYTrack = -999.; // dummy
      auto hmpidXMip = -999.;   // dummy
      auto hmpidYMip = -999.;   // dummy
      auto hmpidNPhotons = hmpData.hmpidNPhotons();
      auto hmpidQMip = hmpData.hmpidQMip();
      auto hmpidClusSize = -999; // dummy
      auto hmpidMom = -999;      // dummy
      auto hmpidPhotsCharge = phots;

      HMPID_001(trackid,
                hmpidSignal,
                hmpidXTrack,
                hmpidYTrack,
                hmpidXMip,
                hmpidYMip,
                hmpidNPhotons,
                hmpidQMip,
                hmpidClusSize,
                hmpidMom,
                hmpidPhotsCharge);
    }
  }
  PROCESS_SWITCH(converterHMPID_000, processConverter, "Process converter (autoset)", true);
};

// Converts the old McCaloLabels_000 table to the new McCaloLabels_001 table where we have a variable size array for associated MCParticles for each calo cell
struct converterMcCaloLabels_000 {
  Produces<aod::McCaloLabels_001> McCaloLabels_001;
  void init(o2::framework::InitContext& initContext) { autoSetProcessFunction(initContext, "McCaloLabels_001", doprocessConverter); }
  void process(aod::BCs const&) {} // Dummy processor in case the other is disabled
  void processConverter(aod::McCaloLabels_000 const& mccalolabelTable)
  {
    std::vector<float> amplitude = {0};
    std::vector<int32_t> particleId = {0};
    for (auto& mccalolabel : mccalolabelTable) {
      particleId[0] = mccalolabel.mcParticleId();
      // Repopulate new table
      McCaloLabels_001(particleId,
                       amplitude);
    }
  }
  PROCESS_SWITCH(converterMcCaloLabels_000, processConverter, "Process converter (autoset)", true);
};

// Converts MCParticle table from version 000 to 001
struct converterStoredMcParticles_000 {
  Produces<aod::StoredMcParticles_001> StoredMcParticles_001;
  void init(o2::framework::InitContext& initContext) { autoSetProcessFunction(initContext, "StoredMcParticles_001", doprocessConverter); }
  void process(aod::BCs const&) {} // Dummy processor in case the other is disabled
  void processConverter(aod::StoredMcParticles_000 const& mcParticles_000)
  {
    for (auto& p : mcParticles_000) {

      std::vector<int> mothers;
      if (p.mother0Id() >= 0) {
        mothers.push_back(p.mother0Id());
      }
      if (p.mother1Id() >= 0) {
        mothers.push_back(p.mother1Id());
      }

      int daughters[2] = {-1, -1};
      if (p.daughter0Id() >= 0 && p.daughter1Id() >= 0) {
        daughters[0] = p.daughter0Id();
        daughters[1] = p.daughter1Id();
      } else if (p.daughter0Id() >= 0) {
        daughters[0] = p.daughter0Id();
        daughters[1] = p.daughter0Id();
      }

      StoredMcParticles_001(p.mcCollisionId(), p.pdgCode(), p.statusCode(), p.flags(),
                            mothers, daughters, p.weight(), p.px(), p.py(), p.pz(), p.e(),
                            p.vx(), p.vy(), p.vz(), p.vt());
    }
  }
  PROCESS_SWITCH(converterStoredMcParticles_000, processConverter, "Process converter (autoset)", true);
};

// Tracks extra table 000 to 002
struct converterStoredTracksExtra_000 {
  Produces<aod::StoredTracksExtra_002> StoredTracksExtra_002;
  void init(o2::framework::InitContext& initContext) { autoSetProcessFunction(initContext, "StoredTracksExtra_002", doprocessConverter); }
  void process(aod::BCs const&) {} // Dummy processor in case the other is disabled
  void processConverter(aod::TracksExtra_000 const& inputTable)
  {
    const int8_t TPCNClsFindableMinusPID = 0;
    for (const auto& track0 : inputTable) {
      uint32_t itsClusterSizes = 0;
      for (int layer = 0; layer < 7; layer++) {
        if (track0.itsClusterMap() & (1 << layer)) {
          itsClusterSizes |= (0xf << (layer * 4));
        }
      }
      StoredTracksExtra_002(track0.tpcInnerParam(),
                            track0.flags(),
                            itsClusterSizes,
                            track0.tpcNClsFindable(),
                            track0.tpcNClsFindableMinusFound(),
                            TPCNClsFindableMinusPID,
                            track0.tpcNClsFindableMinusCrossedRows(),
                            track0.tpcNClsShared(),
                            track0.trdPattern(),
                            track0.itsChi2NCl(),
                            track0.tpcChi2NCl(),
                            track0.trdChi2(),
                            track0.tofChi2(),
                            track0.tpcSignal(),
                            track0.trdSignal(),
                            track0.length(),
                            track0.tofExpMom(),
                            track0.trackEtaEmcal(),
                            track0.trackPhiEmcal(),
                            track0.trackTime(),
                            track0.trackTimeRes());
    }
  }
  PROCESS_SWITCH(converterStoredTracksExtra_000, processConverter, "Process converter (autoset)", true);
};

// Tracks extra table 001 to 002
struct converterStoredTracksExtra_001 {
  Produces<aod::StoredTracksExtra_002> StoredTracksExtra_002;
  void init(o2::framework::InitContext& initContext) { autoSetProcessFunction(initContext, "StoredTracksExtra_002", doprocessConverter); }
  void process(aod::BCs const&) {} // Dummy processor in case the other is disabled
  void processConverter(aod::TracksExtra_001 const& inputTable)
  {
    const int8_t TPCNClsFindableMinusPID = 0;
    for (const auto& track1 : inputTable) {
      StoredTracksExtra_002(track1.tpcInnerParam(),
                            track1.flags(),
                            track1.itsClusterSizes(),
                            track1.tpcNClsFindable(),
                            track1.tpcNClsFindableMinusFound(),
                            TPCNClsFindableMinusPID,
                            track1.tpcNClsFindableMinusCrossedRows(),
                            track1.tpcNClsShared(),
                            track1.trdPattern(),
                            track1.itsChi2NCl(),
                            track1.tpcChi2NCl(),
                            track1.trdChi2(),
                            track1.tofChi2(),
                            track1.tpcSignal(),
                            track1.trdSignal(),
                            track1.length(),
                            track1.tofExpMom(),
                            track1.trackEtaEmcal(),
                            track1.trackPhiEmcal(),
                            track1.trackTime(),
                            track1.trackTimeRes());
    }
  }
  PROCESS_SWITCH(converterStoredTracksExtra_001, processConverter, "Process converter (autoset)", true);
};

/// Spawn the extended table for TracksExtra002 to avoid the call to the internal spawner and a consequent circular dependency
struct spawnerTracksExtra_002 {
  Spawns<aod::TracksExtra_002> TracksExtra_002;
};

struct converterZdcs_000 {
  Produces<aod::Zdcs_001> Zdcs_001;
  void init(o2::framework::InitContext& initContext) { autoSetProcessFunction(initContext, "Zdcs_001", doprocessConverter); }
  void process(aod::BCs const&) {} // Dummy processor in case the other is disabled
  void processConverter(aod::Zdcs_000 const& zdcLegacy, aod::BCs const&)
  {
    for (auto& zdcData : zdcLegacy) {
      // Get legacy information, please
      auto bc = zdcData.bc();
      auto energyZEM1 = zdcData.energyZEM1();
      auto energyZEM2 = zdcData.energyZEM2();
      auto energyCommonZNA = zdcData.energyCommonZNA();
      auto energyCommonZNC = zdcData.energyCommonZNC();
      auto energyCommonZPA = zdcData.energyCommonZPA();
      auto energyCommonZPC = zdcData.energyCommonZPC();
      auto energySectorZNA = zdcData.energySectorZNA();
      auto energySectorZNC = zdcData.energySectorZNC();
      auto energySectorZPA = zdcData.energySectorZPA();
      auto energySectorZPC = zdcData.energySectorZPC();
      auto timeZEM1 = zdcData.timeZEM1();
      auto timeZEM2 = zdcData.timeZEM2();
      auto timeZNA = zdcData.timeZNA();
      auto timeZNC = zdcData.timeZNC();
      auto timeZPA = zdcData.timeZPA();
      auto timeZPC = zdcData.timeZPC();

      // Create variables to initialize Zdcs_001 table
      std::vector<float> zdcEnergy, zdcAmplitudes, zdcTime;
      std::vector<uint8_t> zdcChannelsE, zdcChannelsT;

      // Tie variables in such that they get read correctly later
      zdcEnergy.emplace_back(energyZEM1);
      zdcChannelsE.emplace_back(IdZEM1);
      zdcAmplitudes.emplace_back(energyZEM1); // WARNING: DUMMY VALUE
      zdcTime.emplace_back(timeZEM1);
      zdcChannelsT.emplace_back(IdZEM1);

      zdcEnergy.emplace_back(energyZEM2);
      zdcChannelsE.emplace_back(IdZEM2);
      zdcAmplitudes.emplace_back(energyZEM2); // WARNING: DUMMY VALUE
      zdcTime.emplace_back(timeZEM2);
      zdcChannelsT.emplace_back(IdZEM2);

      zdcEnergy.emplace_back(energyCommonZNA);
      zdcChannelsE.emplace_back(IdZNAC);
      zdcAmplitudes.emplace_back(energyCommonZNA); // WARNING: DUMMY VALUE
      zdcTime.emplace_back(timeZNA);
      zdcChannelsT.emplace_back(IdZNAC);

      zdcEnergy.emplace_back(energyCommonZNC);
      zdcChannelsE.emplace_back(IdZNCC);
      zdcAmplitudes.emplace_back(energyCommonZNC); // WARNING: DUMMY VALUE
      zdcTime.emplace_back(timeZNC);
      zdcChannelsT.emplace_back(IdZNCC);

      zdcEnergy.emplace_back(energyCommonZPA);
      zdcChannelsE.emplace_back(IdZPAC);
      zdcAmplitudes.emplace_back(energyCommonZPA); // WARNING: DUMMY VALUE
      zdcTime.emplace_back(timeZPA);
      zdcChannelsT.emplace_back(IdZPAC);

      zdcEnergy.emplace_back(energyCommonZPC);
      zdcChannelsE.emplace_back(IdZPCC);
      zdcAmplitudes.emplace_back(energyCommonZPC); // WARNING: DUMMY VALUE
      zdcTime.emplace_back(timeZPC);
      zdcChannelsT.emplace_back(IdZPCC);

      for (uint64_t ic = 0; ic < 4; ic++) {
        zdcEnergy.emplace_back(energySectorZNA[ic]);
        zdcChannelsE.emplace_back(IdZNA1 + ic);
        zdcEnergy.emplace_back(energySectorZNC[ic]);
        zdcChannelsE.emplace_back(IdZNC1 + ic);
        zdcEnergy.emplace_back(energySectorZPA[ic]);
        zdcChannelsE.emplace_back(IdZPA1 + ic);
        zdcEnergy.emplace_back(energySectorZPC[ic]);
        zdcChannelsE.emplace_back(IdZPC1 + ic);
      }

      Zdcs_001(bc,
               zdcEnergy,
               zdcChannelsE,
               zdcAmplitudes,
               zdcTime,
               zdcChannelsT);
    }
  }
  PROCESS_SWITCH(converterZdcs_000, processConverter, "Process converter (autoset)", true);
};

struct converterMcCollisions_000 {
  Produces<aod::McCollisions_001> McCollisions_001;
  void init(o2::framework::InitContext& initContext) { autoSetProcessFunction(initContext, "McCollisions_001", doprocessConverter); }
  void process(aod::BCs const&) {} // Dummy processor in case the other is disabled
  void processConverter(aod::McCollisions_000 const& mcCollisionTable)
  {
    for (auto& mcCollision : mcCollisionTable) {

      // Repopulate new table
      McCollisions_001(
        mcCollision.bcId(),
        mcCollision.generatorsID(),
        mcCollision.posX(), mcCollision.posY(), mcCollision.posZ(),
        mcCollision.t(), mcCollision.weight(),
        mcCollision.impactParameter(),
        0.0f); // dummy event plane, not available in _000
    }
  }
  PROCESS_SWITCH(converterMcCollisions_000, processConverter, "Process converter (autoset)", true);
};

struct converterStoredMFTTracks_000 {
  Produces<aod::StoredMFTTracks_001> StoredMFTTracks_001;
  void init(o2::framework::InitContext& initContext) { autoSetProcessFunction(initContext, "StoredMFTTracks_001", doprocessConverter); }
  void process(aod::BCs const&) {} // Dummy processor in case the other is disabled
  void processConverter(aod::MFTTracks_000 const& mftTracks_000)
  {
    for (const auto& track0 : mftTracks_000) {
      uint64_t mftClusterSizesAndTrackFlags = 0;
      int8_t nClusters = track0.nClusters();

      for (int layer = 0; layer < 10; ++layer) {
        mftClusterSizesAndTrackFlags &= ~(0x3fULL << (layer * 6));
        mftClusterSizesAndTrackFlags |= (layer < nClusters) ? (1ULL << (layer * 6)) : 0;
      }
      StoredMFTTracks_001(track0.collisionId(),
                          track0.x(),
                          track0.y(),
                          track0.z(),
                          track0.phi(),
                          track0.tgl(),
                          track0.signed1Pt(),
                          mftClusterSizesAndTrackFlags,
                          track0.chi2(),
                          track0.trackTime(),
                          track0.trackTimeRes());
    }
  }
  PROCESS_SWITCH(converterStoredMFTTracks_000, processConverter, "Process converter (autoset)", true);
};

/// Spawn the extended table for MFTTracks001 to avoid the call to the internal spawner and a consequent circular dependency
struct spawnerMFTTracks {
  Spawns<aod::MFTTracks> MFTTracks;
};

struct converterV0s_000 {
  Produces<aod::V0s_001> V0s_001;
  void init(o2::framework::InitContext& initContext) { autoSetProcessFunction(initContext, "V0s_001", doprocessConverter); }
  void process(aod::BCs const&) {} // Dummy processor in case the other is disabled
  void processConverter(aod::V0s_000 const& v0s, aod::Tracks const&)
  {
    for (auto& v0 : v0s) {
      if (v0.posTrack().collisionId() != v0.negTrack().collisionId()) {
        LOGF(fatal, "V0 %d has inconsistent collision information (%d, %d)", v0.globalIndex(), v0.posTrack().collisionId(), v0.negTrack().collisionId());
      }
      V0s_001(v0.posTrack().collisionId(), v0.posTrackId(), v0.negTrackId());
    }
  }
  PROCESS_SWITCH(converterV0s_000, processConverter, "Process converter (autoset)", true);
};

struct converterV0s_001 {
  Produces<aod::V0s_002> V0s_002;
  void init(o2::framework::InitContext& initContext) { autoSetProcessFunction(initContext, "V0s_002", doprocessConverter); }
  void process(aod::BCs const&) {} // Dummy processor in case the other is disabled
  void processConverter(aod::V0s_001 const& v0s)
  {
    for (auto& v0 : v0s) {
      uint8_t bitMask = static_cast<uint8_t>(1); // first bit on
      V0s_002(v0.collisionId(), v0.posTrackId(), v0.negTrackId(), bitMask);
    }
  }
  PROCESS_SWITCH(converterV0s_001, processConverter, "Process converter (autoset)", true);
};

// NOTE These tasks have to be split because for the cascades, V0s and not V0s_000 are needed
struct converterCascades_000 {
  Produces<aod::Cascades_001> Cascades_001;
  void init(o2::framework::InitContext& initContext) { autoSetProcessFunction(initContext, "Cascades_001", doprocessConverter); }
  void process(aod::BCs const&) {} // Dummy processor in case the other is disabled
  void processConverter(aod::V0s const&, aod::Cascades_000 const& cascades, aod::Tracks const&)
  {
    for (auto& cascade : cascades) {
      if (cascade.bachelor().collisionId() != cascade.v0().posTrack().collisionId() || cascade.v0().posTrack().collisionId() != cascade.v0().negTrack().collisionId()) {
        LOGF(fatal, "Cascade %d has inconsistent collision information (%d, %d, %d) track ids %d %d %d", cascade.globalIndex(), cascade.bachelor().collisionId(),
             cascade.v0().posTrack().collisionId(), cascade.v0().negTrack().collisionId(), cascade.bachelorId(), cascade.v0().posTrackId(), cascade.v0().negTrackId());
      }
      Cascades_001(cascade.bachelor().collisionId(), cascade.v0Id(), cascade.bachelorId());
    }
  }
  PROCESS_SWITCH(converterCascades_000, processConverter, "Process converter (autoset)", true);
};

// Run 2 converters
struct converterRun2BCInfos_000 {
  Produces<aod::Run2BCInfos_001> Run2BCInfos_001;
  void init(o2::framework::InitContext& initContext) { autoSetProcessFunction(initContext, "Run2BCInfos_001", doprocessConverter); }
  void process(aod::BCs const&) {} // Dummy processor in case the other is disabled
  void processConverter(aod::Run2BCInfos_000 const& entries)
  {
    for (const auto& entry : entries) {
      Run2BCInfos_001(entry.eventCuts(),
                      entry.triggerMaskNext50(), entry.l0TriggerInputMask(),
                      entry.spdClustersL0(), entry.spdClustersL1(),
                      entry.spdFiredChipsL0(), entry.spdFiredChipsL1(),
                      entry.spdFiredFastOrL0(), entry.spdFiredFastOrL1(),
                      entry.v0TriggerChargeA(), entry.v0TriggerChargeC(),
                      0, 0);
    }
  }
  PROCESS_SWITCH(converterRun2BCInfos_000, processConverter, "Process converter (autoset)", true);
};

struct converterRun2TrackExtras_000 {
  Produces<aod::Run2TrackExtras_001> Run2TrackExtras_001;
  void init(o2::framework::InitContext& initContext) { autoSetProcessFunction(initContext, "Run2TrackExtras_001", doprocessConverter); }
  void process(aod::BCs const&) {} // Dummy processor in case the other is disabled
  void processConverter(aod::Run2TrackExtras_000 const& tracks)
  {
    for (const auto& track0 : tracks) {
      Run2TrackExtras_001(track0.itsSignal(), 0);
    }
  }
  PROCESS_SWITCH(converterRun2TrackExtras_000, processConverter, "Process converter (autoset)", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  auto workflow = WorkflowSpec{};
  std::unordered_set<std::string> addedConverters; // Track added converters

  // Check if 'aod-metadata-tables' option is available in the config context
  if (cfgc.options().hasOption("aod-metadata-tables")) {
    const std::vector<std::string> tables = cfgc.options().get<std::vector<std::string>>("aod-metadata-tables");

    // Map of table names to their corresponding converter task functions
    std::unordered_map<std::string, std::vector<std::function<void()>>> tableToTasks = {
      {"O2bc", {[&]() { workflow.push_back(adaptAnalysisTask<converterBCs_000>(cfgc)); }}},
      {"O2collision", {[&]() { workflow.push_back(adaptAnalysisTask<converterCollisions_000>(cfgc)); }}},
      {"O2fdd", {[&]() { workflow.push_back(adaptAnalysisTask<converterFDDs_000>(cfgc)); }}},
      {"O2hmpid", {[&]() { workflow.push_back(adaptAnalysisTask<converterHMPID_000>(cfgc)); }}},
      {"O2mccalolabel", {[&]() { workflow.push_back(adaptAnalysisTask<converterMcCaloLabels_000>(cfgc)); }}},
      {"O2mfttrack", {[&]() { workflow.push_back(adaptAnalysisTask<converterStoredMFTTracks_000>(cfgc)); }, [&]() { workflow.push_back(adaptAnalysisTask<spawnerMFTTracks>(cfgc)); }}},
      {"O2v0", {[&]() { workflow.push_back(adaptAnalysisTask<converterV0s_000>(cfgc)); }, [&]() { workflow.push_back(adaptAnalysisTask<converterV0s_000>(cfgc)); }}},
      {"O2v0_001", {[&]() { workflow.push_back(adaptAnalysisTask<converterV0s_001>(cfgc)); }}},
      {"O2cascades", {[&]() { workflow.push_back(adaptAnalysisTask<converterCascades_000>(cfgc)); }}},
      {"O2mccollision", {[&]() { workflow.push_back(adaptAnalysisTask<converterMcCollisions_000>(cfgc)); }}},
      {"O2mcparticle", {[&]() { workflow.push_back(adaptAnalysisTask<converterStoredMcParticles_000>(cfgc)); }}},
      {"O2run2bcinfos", {[&]() { workflow.push_back(adaptAnalysisTask<converterRun2BCInfos_000>(cfgc)); }}},
      {"O2run2trackextras", {[&]() { workflow.push_back(adaptAnalysisTask<converterRun2TrackExtras_000>(cfgc)); }}},
      {"O2trackextra", {[&]() { workflow.push_back(adaptAnalysisTask<converterStoredTracksExtra_000>(cfgc)); }, [&]() { workflow.push_back(adaptAnalysisTask<spawnerTracksExtra_002>(cfgc)); }}},
      {"O2trackextra_001", {[&]() { workflow.push_back(adaptAnalysisTask<converterStoredTracksExtra_001>(cfgc)); }, [&]() { workflow.push_back(adaptAnalysisTask<spawnerTracksExtra_002>(cfgc)); }}},
      {"O2zdc", {[&]() { workflow.push_back(adaptAnalysisTask<converterZdcs_000>(cfgc)); }}},
    };

    // Iterate through the tables and process based on the mapping
    for (auto const& table : tables) {
      LOG(info) << "AOD converter: Table " << table << " checking for converters";

      if (tableToTasks.find(table) != tableToTasks.end()) {
        for (auto const& task : tableToTasks[table]) {
          task();
        }
        LOG(info) << "  + AOD converter: for table " << table << " adding converter";
      } else {
        LOG(info) << "  - AOD converter: for table " << table << " needs no converter";
      }
    }
  } else {
    LOG(warning) << "AOD converter: No tables found in the meta data";
  }
  return workflow;
}
