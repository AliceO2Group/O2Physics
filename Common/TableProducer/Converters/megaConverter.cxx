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
/// \file   megaConverter.cxx
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

  void init(o2::framework::InitContext& initContext)
  {
    if (!doprocessConverter) {
      doprocessConverter.value = isTableRequiredInWorkflow(initContext, "BCs_001");
    }
  }
  void process(aod::Timestamps const&) {}
  void processConverter(aod::BCs_000 const& bcTable)
  {
    for (auto& bc : bcTable) {
      constexpr uint64_t lEmptyTriggerInputs = 0;
      bc_001(bc.runNumber(), bc.globalBC(), bc.triggerMask(), lEmptyTriggerInputs);
    }
  }
  PROCESS_SWITCH(bcConverter, processConverter, "Process converter (autoset if false)", false);
};

// Swaps covariance matrix elements if the data is known to be bogus (collision_000 is bogus)
struct collisionConverter {
  Produces<aod::Collisions_001> Collisions_001;

  Configurable<bool> doNotSwap{"doNotSwap", false, "simple pass-through"};
  Configurable<bool> debug{"debug", false, "flag to save debug histo"};
  Configurable<int> nbins{"nbins", 1, "number of bins in debug histo"};
  Configurable<float> tolerance{"tolerance", 1e-3, "Tolerance for CYY check"};

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&)
  {
    const AxisSpec axisCYYdebug{nbins, -1.0f, +1.0f, ""};
    histos.add("hCYY", "hCYY", kTH1F, {axisCYYdebug});
  }

  void process(aod::Collisions_000 const& collisionTable)
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
};

// Converts FDD table from version 000 to 001
struct FddConverter {
  Produces<aod::FDDs_001> fdd_001;

  void process(aod::FDDs_000 const& fdd_000)
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
};

struct hmpConverter {
  Produces<aod::HMPID_001> HMPID_001;

  void process(aod::HMPID_000 const& hmpLegacy, aod::Tracks const&)
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
};

// Converts the old McCaloLabels_000 table to the new McCaloLabels_001 table where we have a variable size array for associated MCParticles for each calo cell
struct caloLabelConverter {
  Produces<aod::McCaloLabels_001> McCaloLabels_001;

  void process(aod::McCaloLabels_000 const& mccalolabelTable)
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
};

// Converts MCParticle table from version 000 to 001
struct McConverter {
  Produces<aod::StoredMcParticles_001> mcParticles_001;

  void process(aod::StoredMcParticles_000 const& mcParticles_000)
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
};

struct TracksExtraConverter {
  Produces<aod::StoredTracksExtra_001> tracksExtra_001;
  void process(aod::TracksExtra_000 const& tracksExtra_000)
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
};

struct zdcConverter {
  Produces<aod::Zdcs_001> Zdcs_001;

  void process(aod::Zdcs_000 const& zdcLegacy, aod::BCs const&)
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
};

struct mcCollisionConverter {
  Produces<aod::McCollisions_001> mcCollisions_001;

  void process(aod::McCollisions_000 const& mcCollisionTable)
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
};

struct MftTracksConverter {
  Produces<aod::StoredMFTTracks_001> mftTracks_001;
  void process(aod::MFTTracks_000 const& mftTracks_000)
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
};

/// Spawn the extended table for MFTTracks001 to avoid the call to the internal spawner and a consequent circular dependency
struct MFTTracksSpawner {
  Spawns<aod::MFTTracks_001> mftTracks_001;
};

struct MultsExtraConverter {
  Produces<aod::MultsExtra_001> multsExtra_001;
  void process(aod::MultsExtra_000 const& multsExtra_000)
  {
    for (const auto& r : multsExtra_000) {
      multsExtra_001(r.multPVTotalContributors(), r.multPVChi2(),
                     r.multCollisionTimeRes(), r.multRunNumber(), r.multPVz(), r.multSel8(),
                     r.multNTracksHasITS(), r.multNTracksHasTPC(), r.multNTracksHasTOF(),
                     r.multNTracksHasTRD(), r.multNTracksITSOnly(),
                     r.multNTracksTPCOnly(), r.multNTracksITSTPC(),
                     r.multAllTracksTPCOnly(), r.multAllTracksITSTPC(),
                     r.trackOccupancyInTimeRange(),
                     0.0f,
                     r.flags());
    }
  }
};

struct V0Converter {
  Produces<aod::V0s_002> v0s_002;

  void process(aod::V0s_001 const& v0s)
  {
    for (auto& v0 : v0s) {
      uint8_t bitMask = static_cast<uint8_t>(1); // first bit on
      v0s_002(v0.collisionId(), v0.posTrackId(), v0.negTrackId(), bitMask);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  auto workflow = WorkflowSpec{};
  if (cfgc.options().hasOption("aod-metadata-tables")) {
    const std::vector<std::string> tables = cfgc.options().get<std::vector<std::string>>("aod-metadata-tables");
    for (auto t : tables) {
      if (t == "O2cascade_001") {
      } else if (t == "O2cascade_000") {
      } else if (t == "O2bc_001") {
      } else if (t == "O2bc_000") {
        workflow.push_back(adaptAnalysisTask<bcConverter>(cfgc));
      } else if (t == "O2collision_001") {
      } else if (t == "O2collision_000") {
        workflow.push_back(adaptAnalysisTask<collisionConverter>(cfgc));
      } else if (t == "O2fdd_001") {
      } else if (t == "O2fdd_000") {
        workflow.push_back(adaptAnalysisTask<FddConverter>(cfgc));
      } else if (t == "O2hmpid_001") {
      } else if (t == "O2hmpid_000") {
        workflow.push_back(adaptAnalysisTask<hmpConverter>(cfgc));
      } else if (t == "O2mccalolabel_001") {
      } else if (t == "O2mccalolabel_000") {
        workflow.push_back(adaptAnalysisTask<caloLabelConverter>(cfgc));
      } else if (t == "O2mfttrack_001") {
      } else if (t == "O2mfttrack_000") {
        workflow.push_back(adaptAnalysisTask<MftTracksConverter>(cfgc));
        workflow.push_back(adaptAnalysisTask<MFTTracksSpawner>(cfgc));
      } else if (t == "O2multextra_001") {
      } else if (t == "O2multextra_000") {
        workflow.push_back(adaptAnalysisTask<MultsExtraConverter>(cfgc));
      } else if (t == "O2v0_002") {
      } else if (t == "O2v0_001") {
        workflow.push_back(adaptAnalysisTask<V0Converter>(cfgc));
      } else if (t == "O2v0_000") {
        workflow.push_back(adaptAnalysisTask<V0Converter>(cfgc));
      } else if (t == "O2mccollision_001") {
      } else if (t == "O2mccollision_000") {
        workflow.push_back(adaptAnalysisTask<mcCollisionConverter>(cfgc));
      } else if (t == "O2mcparticle_001") {
      } else if (t == "O2mcparticle_000") {
        workflow.push_back(adaptAnalysisTask<McConverter>(cfgc));
      } else if (t == "O2trackextra_001") {
      } else if (t == "O2trackextra_000") {
        workflow.push_back(adaptAnalysisTask<TracksExtraConverter>(cfgc));
      } else if (t == "O2zdc_001") {
      } else if (t == "O2zdc_000") {
        workflow.push_back(adaptAnalysisTask<zdcConverter>(cfgc));
      } else if (t.find("_0") != std::string::npos) {
        LOG(warning) << "Versioned table not covered " << t;
      }
    }
  }
  return workflow;
}
