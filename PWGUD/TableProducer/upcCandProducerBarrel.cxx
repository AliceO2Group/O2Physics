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

/// \file upcCandProducerBarrel.cxx
/// \brief Compact candidate table producer for diffraction and UPC studies.
/// Requires: event selection, propagation, TPC PID and TOF PID services
/// \author Nazar Burmasov (JINR), Evgeny Kryshen (JINR)

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/RCTSelectionFlags.h"
#include "Common/Core/TableHelper.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/LHCConstants.h>
#include <DataFormatsFIT/Triggers.h>
#include <DataFormatsParameters/AggregatedRunInfo.h>
#include <DataFormatsParameters/GRPLHCIFData.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TH2.h>

#include <algorithm>
#include <bitset>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::aod::rctsel;

namespace
{
constexpr int NBCsPerOrbit = o2::constants::lhc::LHCMaxBunches;
constexpr int MaxStoredDistance = 15;
constexpr float MaxFITTime = 30.f;

using DistanceMap = std::vector<std::vector<uint32_t>>;

DistanceMap buildMinimumDistanceMap(std::vector<uint32_t> const& tfIDs,
                                    int64_t bcSOR, int64_t nBCsPerTF,
                                    std::vector<int64_t>& activeBCs, uint32_t maxDistance)
{
  const uint32_t overflowDistance = maxDistance + 1;
  DistanceMap distances(tfIDs.size(), std::vector<uint32_t>(nBCsPerTF, overflowDistance));
  if (activeBCs.empty()) {
    return distances;
  }
  std::sort(activeBCs.begin(), activeBCs.end());
  activeBCs.erase(std::unique(activeBCs.begin(), activeBCs.end()), activeBCs.end());
  size_t nextIndex = 0;
  for (size_t iTF = 0; iTF < tfIDs.size(); ++iTF) {
    const int64_t tfStartBC = bcSOR + tfIDs[iTF] * nBCsPerTF;
    for (int64_t bcInTF = 0; bcInTF < nBCsPerTF; ++bcInTF) {
      const int64_t currentBC = tfStartBC + bcInTF;
      while (nextIndex < activeBCs.size() && activeBCs[nextIndex] < currentBC) {
        ++nextIndex;
      }
      int64_t distance = overflowDistance;
      if (nextIndex < activeBCs.size()) {
        distance = std::min(distance, activeBCs[nextIndex] - currentBC);
      }
      if (nextIndex > 0) {
        distance = std::min(distance, currentBC - activeBCs[nextIndex - 1]);
      }
      distances[iTF][bcInTF] = static_cast<uint32_t>(distance);
    }
  }
  return distances;
}

void fillVetoHistograms(std::shared_ptr<TH2> const& hVetoT00, std::shared_ptr<TH2> const& hVetoTV0,
                        std::shared_ptr<TH2> const& hVetoTVD, int bcInOrbit,
                        uint32_t distanceFT0, uint32_t distanceFV0, uint32_t distanceFDD)
{
  hVetoT00->Fill(bcInOrbit, -1);
  hVetoTV0->Fill(bcInOrbit, -1);
  hVetoTVD->Fill(bcInOrbit, -1);
  for (uint32_t threshold = 0; threshold <= MaxStoredDistance; ++threshold) {
    if (distanceFT0 <= threshold) {
      continue;
    }
    hVetoT00->Fill(bcInOrbit, threshold);
    if (distanceFV0 <= threshold) {
      continue;
    }
    hVetoTV0->Fill(bcInOrbit, threshold);
    if (distanceFDD > threshold) {
      hVetoTVD->Fill(bcInOrbit, threshold);
    }
  }
}
} // namespace

namespace o2::aod::upc_cand_prod_bar
{
DECLARE_SOA_COLUMN(RunNumber, runNumber, int32_t);
DECLARE_SOA_COLUMN(GlobalBC, globalBC, uint64_t);
DECLARE_SOA_COLUMN(TfId, tfId, uint32_t);
DECLARE_SOA_COLUMN(Timestamp, timestamp, uint64_t);
DECLARE_SOA_COLUMN(PosX, posX, float);
DECLARE_SOA_COLUMN(PosY, posY, float);
DECLARE_SOA_COLUMN(PosZ, posZ, float);
DECLARE_SOA_COLUMN(Px, px, std::vector<float>);
DECLARE_SOA_COLUMN(Py, py, std::vector<float>);
DECLARE_SOA_COLUMN(Pz, pz, std::vector<float>);
DECLARE_SOA_COLUMN(TpcSignal, tpcSignal, std::vector<float>);
DECLARE_SOA_COLUMN(TpcNSigmaEl, tpcNSigmaEl, std::vector<float>);
DECLARE_SOA_COLUMN(TpcNSigmaPi, tpcNSigmaPi, std::vector<float>);
DECLARE_SOA_COLUMN(TpcNSigmaKa, tpcNSigmaKa, std::vector<float>);
DECLARE_SOA_COLUMN(TpcNSigmaPr, tpcNSigmaPr, std::vector<float>);
DECLARE_SOA_COLUMN(TofNSigmaEl, tofNSigmaEl, std::vector<float>);
DECLARE_SOA_COLUMN(TofNSigmaPi, tofNSigmaPi, std::vector<float>);
DECLARE_SOA_COLUMN(TofNSigmaKa, tofNSigmaKa, std::vector<float>);
DECLARE_SOA_COLUMN(TofNSigmaPr, tofNSigmaPr, std::vector<float>);
DECLARE_SOA_COLUMN(ItsClusterMap, itsClusterMap, std::vector<uint8_t>);
DECLARE_SOA_COLUMN(NClusters, nClusters, std::vector<uint8_t>);
DECLARE_SOA_COLUMN(Sign, sign, std::vector<int8_t>);
DECLARE_SOA_COLUMN(MinimumDistanceFT0, minimumDistanceFT0, int8_t);
DECLARE_SOA_COLUMN(MinimumDistanceFV0, minimumDistanceFV0, int8_t);
DECLARE_SOA_COLUMN(MinimumDistanceFDD, minimumDistanceFDD, int8_t);
} // namespace o2::aod::upc_cand_prod_bar

namespace o2::aod
{
DECLARE_SOA_TABLE(UPCBarrelCands, "AOD", "UPCBARRELCANDS",
                  upc_cand_prod_bar::RunNumber,
                  upc_cand_prod_bar::GlobalBC,
                  upc_cand_prod_bar::TfId,
                  upc_cand_prod_bar::Timestamp,
                  upc_cand_prod_bar::PosX,
                  upc_cand_prod_bar::PosY,
                  upc_cand_prod_bar::PosZ,
                  upc_cand_prod_bar::Px,
                  upc_cand_prod_bar::Py,
                  upc_cand_prod_bar::Pz,
                  upc_cand_prod_bar::TpcSignal,
                  upc_cand_prod_bar::TpcNSigmaEl,
                  upc_cand_prod_bar::TpcNSigmaPi,
                  upc_cand_prod_bar::TpcNSigmaKa,
                  upc_cand_prod_bar::TpcNSigmaPr,
                  upc_cand_prod_bar::TofNSigmaEl,
                  upc_cand_prod_bar::TofNSigmaPi,
                  upc_cand_prod_bar::TofNSigmaKa,
                  upc_cand_prod_bar::TofNSigmaPr,
                  upc_cand_prod_bar::ItsClusterMap,
                  upc_cand_prod_bar::NClusters,
                  upc_cand_prod_bar::Sign,
                  upc_cand_prod_bar::MinimumDistanceFT0,
                  upc_cand_prod_bar::MinimumDistanceFV0,
                  upc_cand_prod_bar::MinimumDistanceFDD);
} // namespace o2::aod

struct UpcCandProducerBarrel {
  Produces<aod::UPCBarrelCands> selectedCandidates;
  HistogramRegistry registry{"registry", {}};
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  int cachedRunNumber = -1;
  int64_t bcSOR = 0;
  int64_t nBCsPerTF = 0;
  int configuredOrbitsPerTF = -1;
  int configuredTFStartBorder = -1;
  int configuredTFEndBorder = -1;
  int tfStartBorder = 0;
  int tfEndBorder = 0;
  std::bitset<NBCsPerOrbit> collidingBCs;

  Configurable<int> nTracks{"nTracks", 2, "Required N tracks in selected collisions"};
  Configurable<float> maxAbsEta{"maxAbsEta", 0.8f, "Maximum |eta| of selected tracks"};
  Configurable<float> minPt{"minPt", 0.2f, "Minimum pT of selected tracks (GeV/c)"};
  Configurable<int> vetoBCWindow{"vetoBCWindow", 0, "FIT veto half-window in BC; <0 disables"};
  Configurable<int> vetoFT0{"vetoFT0", 1, "FT0 veto: 0=off, 1=on"};
  Configurable<int> vetoFV0{"vetoFV0", 0, "FV0 veto: 0=off, 1=on"};
  Configurable<int> vetoFDD{"vetoFDD", 0, "FDD veto: 0=off, 1=on"};

  using CollisionsWithSels = soa::Join<aod::Collisions, aod::EvSels>;
  using BCsWithSels = soa::Join<aod::BCsWithTimestamps, aod::BcSels>;
  using TracksWithPID = soa::Join<aod::Tracks, aod::TracksExtra,
                                  aod::pidTPCEl, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr,
                                  aod::pidTOFEl, aod::pidTOFPi, aod::pidTOFKa, aod::pidTOFPr>;
  Preslice<TracksWithPID> tracksPerCollision = aod::track::collisionId;

  RCTFlagsChecker rctChecker{kFDDBad, kFT0Bad, kFV0Bad, kITSBad, kITSLimAccMCRepr, kTPCBadTracking, kTPCLimAccMCRepr, kTPCBadPID, kTOFBad, kTOFLimAccMCRepr, kCcdbObjectLoaded};

  void init(InitContext& initContext)
  {
    ccdb->setLocalObjectValidityChecking();
    auto inheritEventSelectionOption = [&](char const* name, int& value) {
      if (!o2::common::core::getTaskOptionValue(initContext, "eventselection-run3", name, value, false)) {
        LOGF(fatal, "Could not inherit option %s from eventselection-run3", name);
      }
    };
    inheritEventSelectionOption("bcselOpts.NumberOfOrbitsPerTF", configuredOrbitsPerTF);
    inheritEventSelectionOption("bcselOpts.TimeFrameStartBorderMargin", configuredTFStartBorder);
    inheritEventSelectionOption("bcselOpts.TimeFrameEndBorderMargin", configuredTFEndBorder);
    registry.add("hProcessedTFs", "Processed TFs;;TFs", HistType::kTH1D, {{1, 0., 1.}});
    registry.add("hTVX", "TVX counts per run;run;TVX BCs", HistType::kTH1D, {{NBCsPerOrbit, 0., NBCsPerOrbit}});
    registry.add("hTVXRCT", "TVX counts per run;run;TVX BCs", HistType::kTH1D, {{NBCsPerOrbit, 0., NBCsPerOrbit}});
    registry.add("hVetoT00", ";BC in orbit;veto range (BC);colliding BC count, no RCT mask", HistType::kTH2D, {{NBCsPerOrbit, 0., NBCsPerOrbit}, {17, -1.5, 15.5}});
    registry.add("hVetoTV0", ";BC in orbit;veto range (BC);colliding BC count, no RCT mask", HistType::kTH2D, {{NBCsPerOrbit, 0., NBCsPerOrbit}, {17, -1.5, 15.5}});
    registry.add("hVetoTVD", ";BC in orbit;veto range (BC);colliding BC count, no RCT mask", HistType::kTH2D, {{NBCsPerOrbit, 0., NBCsPerOrbit}, {17, -1.5, 15.5}});
    registry.add("hVetoT00RCT", ";BC in orbit;veto range (BC);colliding BC count, RCT mask", HistType::kTH2D, {{NBCsPerOrbit, 0., NBCsPerOrbit}, {17, -1.5, 15.5}});
    registry.add("hVetoTV0RCT", ";BC in orbit;veto range (BC);colliding BC count, RCT mask", HistType::kTH2D, {{NBCsPerOrbit, 0., NBCsPerOrbit}, {17, -1.5, 15.5}});
    registry.add("hVetoTVDRCT", ";BC in orbit;veto range (BC);colliding BC count, RCT mask", HistType::kTH2D, {{NBCsPerOrbit, 0., NBCsPerOrbit}, {17, -1.5, 15.5}});
    LOGF(info, "Veto config: window=%d FT0=%d FV0=%d FDD=%d", vetoBCWindow.value, vetoFT0.value, vetoFV0.value, vetoFDD.value);
  }

  void updateRunInfo(int runNumber)
  {
    if (runNumber == cachedRunNumber) {
      return;
    }
    const auto runInfo = o2::parameters::AggregatedRunInfo::buildAggregatedRunInfo(ccdb->instance(), runNumber);
    const int64_t orbitsPerTF = configuredOrbitsPerTF < 0 ? runInfo.orbitsPerTF : configuredOrbitsPerTF;
    const bool needsEventSelectionParams = configuredTFStartBorder < 0 || configuredTFEndBorder < 0;
    const auto* eventSelectionParams = needsEventSelectionParams ? ccdb->getForTimeStamp<EventSelectionParams>("EventSelection/EventSelectionParams", runInfo.sor / 2 + runInfo.eor / 2) : nullptr;
    if (orbitsPerTF <= 0 || !runInfo.grpLHC || (needsEventSelectionParams && !eventSelectionParams)) {
      LOGF(fatal, "Incomplete run information for run %d", runNumber);
    }
    bcSOR = runInfo.orbitSOR * NBCsPerOrbit;
    nBCsPerTF = orbitsPerTF * NBCsPerOrbit;
    tfStartBorder = configuredTFStartBorder < 0 ? eventSelectionParams->fTimeFrameStartBorderMargin : configuredTFStartBorder;
    tfEndBorder = configuredTFEndBorder < 0 ? eventSelectionParams->fTimeFrameEndBorderMargin : configuredTFEndBorder;
    collidingBCs = runInfo.grpLHC->getBunchFilling().getBCPattern();
    cachedRunNumber = runNumber;
  }

  void process(CollisionsWithSels const& collisions, TracksWithPID const& tracks, BCsWithSels const& bcs, aod::FT0s const& ft0s, aod::FV0As const& fv0s, aod::FDDs const& fdds)
  {
    updateRunInfo(bcs.begin().runNumber());
    std::vector<uint32_t> tfIDs;
    std::vector<uint8_t> tfPassesRCT;
    std::vector<size_t> localTFIndex;
    localTFIndex.reserve(bcs.size());
    for (const auto& bc : bcs) {
      const auto tfID = static_cast<uint32_t>((static_cast<int64_t>(bc.globalBC()) - bcSOR) / nBCsPerTF);
      if (tfIDs.empty() || tfID != tfIDs.back()) {
        tfIDs.push_back(tfID);
        tfPassesRCT.push_back(static_cast<uint8_t>(rctChecker(bc)));
      }
      localTFIndex.push_back(tfIDs.size() - 1);
    }
    registry.fill(HIST("hProcessedTFs"), 0.5, static_cast<double>(tfIDs.size()));

    auto hTVX = registry.get<TH1>(HIST("hTVX"));
    auto hTVXRCT = registry.get<TH1>(HIST("hTVXRCT"));

    std::vector<int64_t> bcsWithFT0;
    std::vector<int64_t> bcsWithFV0;
    std::vector<int64_t> bcsWithFDD;
    bcsWithFT0.reserve(ft0s.size());
    bcsWithFV0.reserve(fv0s.size());
    bcsWithFDD.reserve(fdds.size());

    for (const auto& ft0 : ft0s) {
      const auto bc = ft0.bc_as<BCsWithSels>();
      const auto gbc = bc.globalBC();
      const auto bcInOrbit = gbc % NBCsPerOrbit;
      if (ft0.timeA() < MaxFITTime || ft0.timeC() < MaxFITTime) {
        bcsWithFT0.push_back(gbc);
      }
      if (!collidingBCs[bcInOrbit] || !TESTBIT(ft0.triggerMask(), o2::fit::Triggers::bitVertex)) {
        continue;
      }
      hTVX->Fill(bcInOrbit);
      if (!bc.selection_bit(aod::evsel::kNoTimeFrameBorder) || !rctChecker(bc)) {
        continue;
      }
      hTVXRCT->Fill(bcInOrbit);
    }
    for (const auto& fv0 : fv0s) {
      if (fv0.time() < MaxFITTime) {
        bcsWithFV0.push_back(fv0.bc_as<BCsWithSels>().globalBC());
      }
    }
    for (const auto& fdd : fdds) {
      if (fdd.timeA() < MaxFITTime || fdd.timeC() < MaxFITTime) {
        bcsWithFDD.push_back(fdd.bc_as<BCsWithSels>().globalBC());
      }
    }

    const uint32_t scanDistance = std::max(MaxStoredDistance, vetoBCWindow.value);
    const DistanceMap minDistanceFT0 = buildMinimumDistanceMap(tfIDs, bcSOR, nBCsPerTF, bcsWithFT0, scanDistance);
    const DistanceMap minDistanceFV0 = buildMinimumDistanceMap(tfIDs, bcSOR, nBCsPerTF, bcsWithFV0, scanDistance);
    const DistanceMap minDistanceFDD = buildMinimumDistanceMap(tfIDs, bcSOR, nBCsPerTF, bcsWithFDD, scanDistance);

    auto hVetoT00 = registry.get<TH2>(HIST("hVetoT00"));
    auto hVetoTV0 = registry.get<TH2>(HIST("hVetoTV0"));
    auto hVetoTVD = registry.get<TH2>(HIST("hVetoTVD"));
    auto hVetoT00RCT = registry.get<TH2>(HIST("hVetoT00RCT"));
    auto hVetoTV0RCT = registry.get<TH2>(HIST("hVetoTV0RCT"));
    auto hVetoTVDRCT = registry.get<TH2>(HIST("hVetoTVDRCT"));
    for (size_t iTF = 0; iTF < tfIDs.size(); ++iTF) {
      const int64_t tfStartBC = bcSOR + static_cast<int64_t>(tfIDs[iTF]) * nBCsPerTF;
      for (int64_t bcInTF = tfStartBorder + 1; bcInTF < nBCsPerTF - tfEndBorder; ++bcInTF) {
        const int64_t globalBC = tfStartBC + bcInTF;
        const int bcInOrbit = globalBC % NBCsPerOrbit;
        if (!collidingBCs[bcInOrbit]) {
          continue;
        }
        const auto distanceFT0 = minDistanceFT0[iTF][bcInTF];
        const auto distanceFV0 = minDistanceFV0[iTF][bcInTF];
        const auto distanceFDD = minDistanceFDD[iTF][bcInTF];
        fillVetoHistograms(hVetoT00, hVetoTV0, hVetoTVD, bcInOrbit, distanceFT0, distanceFV0, distanceFDD);
        if (tfPassesRCT[iTF] != 0u) {
          fillVetoHistograms(hVetoT00RCT, hVetoTV0RCT, hVetoTVDRCT, bcInOrbit, distanceFT0, distanceFV0, distanceFDD);
        }
      }
    }

    for (const auto& collision : collisions) {
      if (collision.numContrib() != nTracks || !collision.selection_bit(aod::evsel::kNoTimeFrameBorder) || !rctChecker(collision)) {
        continue;
      }
      auto bc = collision.bc_as<BCsWithSels>();
      auto gbc = bc.globalBC();
      const auto iTF = localTFIndex[bc.globalIndex()];
      const int64_t bcInTF = (static_cast<int64_t>(gbc) - bcSOR) % nBCsPerTF;
      const auto distFT0 = minDistanceFT0[iTF][bcInTF];
      const auto distFV0 = minDistanceFV0[iTF][bcInTF];
      const auto distFDD = minDistanceFDD[iTF][bcInTF];
      const auto vetoWindow = static_cast<uint32_t>(vetoBCWindow.value);
      if (vetoBCWindow.value >= 0 && (((vetoFT0.value != 0) && distFT0 <= vetoWindow) || ((vetoFV0.value != 0) && distFV0 <= vetoWindow) || ((vetoFDD.value != 0) && distFDD <= vetoWindow))) {
        continue;
      }
      std::vector<float> px, py, pz;
      std::vector<float> tpcSignal, tpcNSigmaEl, tpcNSigmaPi, tpcNSigmaKa, tpcNSigmaPr;
      std::vector<float> tofNSigmaEl, tofNSigmaPi, tofNSigmaKa, tofNSigmaPr;
      std::vector<uint8_t> itsClusterMap;
      std::vector<uint8_t> nClusters;
      std::vector<int8_t> sign;
      for (const auto& track : tracks.sliceBy(tracksPerCollision, collision.globalIndex())) {
        if (!track.isPVContributor() || !track.hasITS() || !track.hasTPC() || std::abs(track.eta()) > maxAbsEta.value || track.pt() < minPt.value) {
          continue;
        }
        px.push_back(track.px());
        py.push_back(track.py());
        pz.push_back(track.pz());
        tpcSignal.push_back(track.tpcSignal());
        tpcNSigmaEl.push_back(track.tpcNSigmaEl());
        tpcNSigmaPi.push_back(track.tpcNSigmaPi());
        tpcNSigmaKa.push_back(track.tpcNSigmaKa());
        tpcNSigmaPr.push_back(track.tpcNSigmaPr());
        tofNSigmaEl.push_back(track.tofNSigmaEl());
        tofNSigmaPi.push_back(track.tofNSigmaPi());
        tofNSigmaKa.push_back(track.tofNSigmaKa());
        tofNSigmaPr.push_back(track.tofNSigmaPr());
        itsClusterMap.push_back(track.itsClusterMap());
        nClusters.push_back(track.tpcNClsFound());
        sign.push_back(track.sign());
      }
      if (px.size() != static_cast<size_t>(nTracks)) {
        continue;
      }
      selectedCandidates(bc.runNumber(), gbc, tfIDs[iTF], bc.timestamp(), collision.posX(), collision.posY(), collision.posZ(),
                         px, py, pz, tpcSignal, tpcNSigmaEl, tpcNSigmaPi, tpcNSigmaKa, tpcNSigmaPr,
                         tofNSigmaEl, tofNSigmaPi, tofNSigmaKa, tofNSigmaPr, itsClusterMap, nClusters, sign,
                         std::min<uint32_t>(distFT0, MaxStoredDistance + 1),
                         std::min<uint32_t>(distFV0, MaxStoredDistance + 1),
                         std::min<uint32_t>(distFDD, MaxStoredDistance + 1));
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& context) { return WorkflowSpec{adaptAnalysisTask<UpcCandProducerBarrel>(context)}; }
