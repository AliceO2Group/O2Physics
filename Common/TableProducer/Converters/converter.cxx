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
#include "Common/DataModel/Multiplicity.h"

// O2Physics includes
#include "TableHelper.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::zdc;

// Converts bc_000 into bc_001
struct bcConverter {
  Produces<aod::BCs_001> bc_001;
  static constexpr std::string tableTo = "BCs_001";
  void process(aod::BCs_000 const& bcTable) // BC converter is always needed
  {
    for (auto& bc : bcTable) {
      constexpr uint64_t lEmptyTriggerInputs = 0;
      bc_001(bc.runNumber(), bc.globalBC(), bc.triggerMask(), lEmptyTriggerInputs);
    }
  }
};

// Swaps covariance matrix elements if the data is known to be bogus (collision_000 is bogus)
struct collisionConverter {
  Produces<aod::Collisions_001> Collisions_001;
  static constexpr std::string tableTo = "Collisions_001";
  Configurable<bool> doNotSwap{"doNotSwap", false, "simple pass-through"};
  Configurable<bool> debug{"debug", false, "flag to save debug histo"};
  Configurable<int> nbins{"nbins", 1, "number of bins in debug histo"};
  Configurable<float> tolerance{"tolerance", 1e-3, "Tolerance for CYY check"};
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  void init(o2::framework::InitContext& initContext)
  {
    enableFlagIfTableRequired(initContext, tableTo, doprocessConverter);
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
  PROCESS_SWITCH(collisionConverter, processConverter, "Process converter (autoset)", false);
};

// Converts FDD table from version 000 to 001
struct FddConverter {
  Produces<aod::FDDs_001> fdd_001;
  static constexpr std::string tableTo = "FDDs_001";
  void init(o2::framework::InitContext& initContext) { enableFlagIfTableRequired(initContext, tableTo, doprocessConverter); }
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

      fdd_001(p.bcId(), chargeA, chargeC,
              p.timeA(), p.timeC(), p.triggerMask());
    }
  }
  PROCESS_SWITCH(FddConverter, processConverter, "Process converter (autoset if false)", false);
};

struct hmpConverter {
  Produces<aod::HMPID_001> HMPID_001;
  static constexpr std::string tableTo = "HMPID_001";
  void init(o2::framework::InitContext& initContext) { enableFlagIfTableRequired(initContext, tableTo, doprocessConverter); }
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
  PROCESS_SWITCH(hmpConverter, processConverter, "Process converter (autoset if false)", false);
};

// Converts the old McCaloLabels_000 table to the new McCaloLabels_001 table where we have a variable size array for associated MCParticles for each calo cell
struct caloLabelConverter {
  Produces<aod::McCaloLabels_001> McCaloLabels_001;
  static constexpr std::string tableTo = "McCaloLabels_001";
  void init(o2::framework::InitContext& initContext) { enableFlagIfTableRequired(initContext, tableTo, doprocessConverter); }
  void process(aod::BCs const&) {} // Dummy processor in case the other is disabled
  void processConverter(aod::McCaloLabels_000 const& mccalolabelTable)
  {
    std::vector<float> amplitude = {0};
    std::vector<int32_t> particleId = {0};
    for (auto& mccalolabel : mccalolabelTable) {
      particleId[0] = mccalolabel.mcParticleId();
      // Repopulate new table
      McCaloLabels_001(
        particleId,
        amplitude);
    }
  }
  PROCESS_SWITCH(caloLabelConverter, processConverter, "Process converter (autoset if false)", false);
};

// Converts MCParticle table from version 000 to 001
struct McConverter {
  Produces<aod::StoredMcParticles_001> mcParticles_001;
  static constexpr std::string tableTo = "StoredMcParticles_001";
  void init(o2::framework::InitContext& initContext) { enableFlagIfTableRequired(initContext, tableTo, doprocessConverter); }
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

      mcParticles_001(p.mcCollisionId(), p.pdgCode(), p.statusCode(), p.flags(),
                      mothers, daughters, p.weight(), p.px(), p.py(), p.pz(), p.e(),
                      p.vx(), p.vy(), p.vz(), p.vt());
    }
  }
  PROCESS_SWITCH(McConverter, processConverter, "Process converter (autoset if false)", false);
};

struct TracksExtraConverter {
  Produces<aod::StoredTracksExtra_001> tracksExtra_001;
  static constexpr std::string tableTo = "StoredTracksExtra_001";
  void init(o2::framework::InitContext& initContext) { enableFlagIfTableRequired(initContext, tableTo, doprocessConverter); }
  void process(aod::BCs const&) {} // Dummy processor in case the other is disabled
  void processConverter(aod::TracksExtra_000 const& tracksExtra_000)
  {

    for (const auto& track0 : tracksExtra_000) {

      uint32_t itsClusterSizes = 0;
      for (int layer = 0; layer < 7; layer++) {
        if (track0.itsClusterMap() & (1 << layer)) {
          itsClusterSizes |= (0xf << (layer * 4));
        }
      }

      tracksExtra_001(track0.tpcInnerParam(),
                      track0.flags(),
                      itsClusterSizes,
                      track0.tpcNClsFindable(),
                      track0.tpcNClsFindableMinusFound(),
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
  PROCESS_SWITCH(TracksExtraConverter, processConverter, "Process converter (autoset if false)", true);
};

struct TracksExtraSpawner {
  Spawns<aod::TracksExtra> tracksExtra;
};

struct zdcConverter {
  Produces<aod::Zdcs_001> Zdcs_001;
  static constexpr std::string tableTo = "Zdcs_001";
  void init(o2::framework::InitContext& initContext) { enableFlagIfTableRequired(initContext, tableTo, doprocessConverter); }
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
  PROCESS_SWITCH(zdcConverter, processConverter, "Process converter (autoset if false)", false);
};

struct mcCollisionConverter {
  Produces<aod::McCollisions_001> mcCollisions_001;
  static constexpr std::string tableTo = "McCollisions_001";
  void init(o2::framework::InitContext& initContext) { enableFlagIfTableRequired(initContext, tableTo, doprocessConverter); }
  void process(aod::BCs const&) {} // Dummy processor in case the other is disabled
  void processConverter(aod::McCollisions_000 const& mcCollisionTable)
  {
    for (auto& mcCollision : mcCollisionTable) {

      // Repopulate new table
      mcCollisions_001(
        mcCollision.bcId(),
        mcCollision.generatorsID(),
        mcCollision.posX(), mcCollision.posY(), mcCollision.posZ(),
        mcCollision.t(), mcCollision.weight(),
        mcCollision.impactParameter(),
        0.0f); // dummy event plane, not available in _000
    }
  }
  PROCESS_SWITCH(mcCollisionConverter, processConverter, "Process converter (autoset if false)", false);
};

struct MftTracksConverter {
  Produces<aod::StoredMFTTracks_001> mftTracks_001;
  static constexpr std::string tableTo = "StoredMFTTracks_001";
  void init(o2::framework::InitContext& initContext) { enableFlagIfTableRequired(initContext, tableTo, doprocessConverter); }
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

      mftTracks_001(track0.collisionId(),
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
  PROCESS_SWITCH(MftTracksConverter, processConverter, "Process converter (autoset if false)", false);
};

/// Spawn the extended table for MFTTracks001 to avoid the call to the internal spawner and a consequent circular dependency
struct MFTTracksSpawner {
  Spawns<aod::MFTTracks_001> mftTracks_001;
};

struct V0Converter {
  Produces<aod::V0s_002> v0s_002;
  static constexpr std::string tableTo = "V0s_002";
  void init(o2::framework::InitContext& initContext) { enableFlagIfTableRequired(initContext, tableTo, doprocessConverter); }
  void process(aod::BCs const&) {} // Dummy processor in case the other is disabled
  void processConverter(aod::V0s_001 const& v0s)
  {
    for (auto& v0 : v0s) {
      uint8_t bitMask = static_cast<uint8_t>(1); // first bit on
      v0s_002(v0.collisionId(), v0.posTrackId(), v0.negTrackId(), bitMask);
    }
  }
  PROCESS_SWITCH(V0Converter, processConverter, "Process converter (autoset if false)", false);
};

struct WeakDecayIndicesV0 {
  Produces<aod::V0s_001> v0s_001;
  static constexpr std::string tableTo = "V0s_001";
  void init(o2::framework::InitContext& initContext) { enableFlagIfTableRequired(initContext, tableTo, doprocessConverter); }
  void process(aod::BCs const&) {} // Dummy processor in case the other is disabled
  void processConverter(aod::V0s_000 const& v0s, aod::Tracks const&)
  {
    for (auto& v0 : v0s) {
      if (v0.posTrack().collisionId() != v0.negTrack().collisionId()) {
        LOGF(fatal, "V0 %d has inconsistent collision information (%d, %d)", v0.globalIndex(), v0.posTrack().collisionId(), v0.negTrack().collisionId());
      }
      v0s_001(v0.posTrack().collisionId(), v0.posTrackId(), v0.negTrackId());
    }
  }
  PROCESS_SWITCH(WeakDecayIndicesV0, processConverter, "Process converter (autoset if false)", false);
};

// NOTE These tasks have to be split because for the cascades, V0s and not V0s_000 are needed
struct WeakDecayIndicesCascades {
  Produces<aod::Cascades_001> cascades_001;
  static constexpr std::string tableTo = "Cascades_001";
  void init(o2::framework::InitContext& initContext) { enableFlagIfTableRequired(initContext, tableTo, doprocessConverter); }
  void process(aod::BCs const&) {} // Dummy processor in case the other is disabled
  void processConverter(aod::V0s const&, aod::Cascades_000 const& cascades, aod::Tracks const&)
  {
    for (auto& cascade : cascades) {
      if (cascade.bachelor().collisionId() != cascade.v0().posTrack().collisionId() || cascade.v0().posTrack().collisionId() != cascade.v0().negTrack().collisionId()) {
        LOGF(fatal, "Cascade %d has inconsistent collision information (%d, %d, %d) track ids %d %d %d", cascade.globalIndex(), cascade.bachelor().collisionId(),
             cascade.v0().posTrack().collisionId(), cascade.v0().negTrack().collisionId(), cascade.bachelorId(), cascade.v0().posTrackId(), cascade.v0().negTrackId());
      }
      cascades_001(cascade.bachelor().collisionId(), cascade.v0Id(), cascade.bachelorId());
    }
  }
  PROCESS_SWITCH(WeakDecayIndicesCascades, processConverter, "Process converter (autoset if false)", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  auto workflow = WorkflowSpec{};
  if (cfgc.options().hasOption("aod-metadata-tables")) {
    const std::vector<std::string> tables = cfgc.options().get<std::vector<std::string>>("aod-metadata-tables");
    for (auto t : tables) {
      LOG(info) << "AOD converter: Table " << t << " checking for converters";
      if (t == "O2bc") {
        LOG(info) << "  + AOD converter: Table " << t << " found, adding converter task";
        workflow.push_back(adaptAnalysisTask<bcConverter>(cfgc));
        continue;
      }
      if (t == "O2collision") {
        LOG(info) << "  + AOD converter: Table " << t << " found, adding converter task";
        workflow.push_back(adaptAnalysisTask<collisionConverter>(cfgc));
        continue;
      }
      if (t == "O2fdd") {
        LOG(info) << "  + AOD converter: Table " << t << " found, adding converter task";
        workflow.push_back(adaptAnalysisTask<FddConverter>(cfgc));
        continue;
      }
      if (t == "O2hmpid") {
        LOG(info) << "  + AOD converter: Table " << t << " found, adding converter task";
        workflow.push_back(adaptAnalysisTask<hmpConverter>(cfgc));
        continue;
      }
      if (t == "O2mccalolabel") {
        LOG(info) << "  + AOD converter: Table " << t << " found, adding converter task";
        workflow.push_back(adaptAnalysisTask<caloLabelConverter>(cfgc));
        continue;
      }
      if (t == "O2mfttrack") {
        LOG(info) << "  + AOD converter: Table " << t << " found, adding converter task";
        workflow.push_back(adaptAnalysisTask<MftTracksConverter>(cfgc));
        workflow.push_back(adaptAnalysisTask<MFTTracksSpawner>(cfgc));
        continue;
      }
      if (t == "O2v0") {
        LOG(info) << "  + AOD converter: Table " << t << " found, adding converter task";
        workflow.push_back(adaptAnalysisTask<WeakDecayIndicesV0>(cfgc));
        workflow.push_back(adaptAnalysisTask<V0Converter>(cfgc));
        continue;
      }
      if (t == "O2v0_001") {
        LOG(info) << "  + AOD converter: Table " << t << " found, adding converter task";
        workflow.push_back(adaptAnalysisTask<V0Converter>(cfgc));
        continue;
      }
      if (t == "O2cascades") {
        LOG(info) << "  + AOD converter: Table " << t << " found, adding converter task";
        workflow.push_back(adaptAnalysisTask<WeakDecayIndicesCascades>(cfgc));
        continue;
      }
      if (t == "O2mccollision") {
        LOG(info) << "  + AOD converter: Table " << t << " found, adding converter task mcCollisionConverter";
        workflow.push_back(adaptAnalysisTask<mcCollisionConverter>(cfgc));
        continue;
      }
      if (t == "O2mcparticle") {
        LOG(info) << "  + AOD converter: Table " << t << " found, adding converter task";
        workflow.push_back(adaptAnalysisTask<McConverter>(cfgc));
        continue;
      }
      if (t == "O2trackextra") {
        LOG(info) << "  + AOD converter: Table " << t << " found, adding converter task";
        workflow.push_back(adaptAnalysisTask<TracksExtraConverter>(cfgc));
        workflow.push_back(adaptAnalysisTask<TracksExtraSpawner>(cfgc));
        continue;
      }
      if (t == "O2zdc") {
        LOG(info) << "  + AOD converter: Table " << t << " found, adding converter task";
        workflow.push_back(adaptAnalysisTask<zdcConverter>(cfgc));
        continue;
      }
      LOG(info) << "  - AOD converter: Table " << t << " needs no converter";
    }
  } else {
    LOG(warning) << "AOD converter: No tables found in the meta data";
  }
  return workflow;
}
