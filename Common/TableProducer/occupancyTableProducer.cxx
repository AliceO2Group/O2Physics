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
/// \brief  Occupancy Table Producer : TPC PID - Calibration
///         Occupancy calculater using tracks which have entry for collision and trackQA tables
///         Ambg tracks were not used
/// \author Rahul Verma (rahul.verma@iitb.ac.in) :: Marian I Ivanov (marian.ivanov@cern.ch)

#include <vector>
#include <unordered_map>
#include <algorithm>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"

#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"

#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/DataModel/OccupancyTables.h"

#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsFT0/Digit.h"
#include "DataFormatsParameters/GRPLHCIFData.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

int32_t nBCsPerOrbit = o2::constants::lhc::LHCMaxBunches;
// const int nBCinTFgrp80 = 1425;
// nOrbitsPerTF = run < 534133 ? 128 : 32;
// for 128 => nBCsPerTF = 456192 , for 32 => nBCsPerTF = 114048
const int nBCinTF = 114048;         /// CCDB value // to be obtained from CCDB in future
const int nBCinDrift = 114048 / 32; /// to get from ccdb in future
const int arraySize = 3;            // Max no timeframes that can be present in a dataframe

struct OccupancyTableProducer {

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  // declare production of tables
  Produces<aod::OccIndexTable> genOccIndexTable;
  Produces<aod::BCTFinfoTable> genBCTFinfoTable;
  Produces<aod::OccsBCsList> genOccsBCsList;

  Produces<aod::OccsDet> genOccsDet;
  Produces<aod::OccsTrackMult> genOccsTrackMult;
  Produces<aod::OccsMultExtra> genOccsMultExtra;
  Produces<aod::OccsRobust> genOccsRobust;

  Produces<aod::OccsMeanDet> genOccsMeanDet;
  Produces<aod::OccsMeanTrkMult> genOccsMeanTrkMult;
  Produces<aod::OccsMnMultExtra> genOccsMnMultExtra;
  Produces<aod::OccsMeanRobust> genOccsMeanRobust;

  Configurable<int> customOrbitOffset{"customOrbitOffset", 0, "customOrbitOffset for MC"};
  Configurable<int> grouping{"grouping", 80, "grouping of BCs"};

  Configurable<bool> buildOccIndexTable{"buildOccIndexTable", true, "builder of table OccIndexTable"};
  Configurable<bool> buildBCTFinfoTable{"buildBCTFinfoTable", true, "builder of table BCTFinfoTable"};
  Configurable<bool> buildOccsBCsList{"buildOccsBCsList", true, "builder of table OccsBCsList"};

  Configurable<bool> buildOccsDet{"buildOccsDet", true, "builder of table OccsDet"};
  Configurable<bool> buildOccsTrackMult{"buildOccsTrackMult", true, "builder of table OccsTrackMult"};
  Configurable<bool> buildOccsMultExtra{"buildOccsMultExtra", true, "builder of table OccsMultExtra"};
  Configurable<bool> buildOccsRobust{"buildOccsRobust", true, "builder of table OccsRobust"};

  Configurable<bool> buildOccsMeanDet{"buildOccsMeanDet", true, "builder of table OccsMeanDet"};
  Configurable<bool> buildOccsMeanTrkMult{"buildOccsMeanTrkMult", true, "builder of table OccsMeanTrkMult"};
  Configurable<bool> buildOccsMnMultExtra{"buildOccsMnMultExtra", true, "builder of table OccsMnMultExtra"};
  Configurable<bool> buildOccsMeanRobust{"buildOccsMeanRobust", true, "builder of table OccsMeanRobust"};

  // Histogram registry;
  HistogramRegistry recoEvent{"recoEvent", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // Data Structures for Occupancy estimation
  std::array<int, arraySize> tfList;
  std::array<std::vector<int64_t>, arraySize> bcTFMap;
  std::array<std::vector<float>, arraySize> occPrimUnfm80;
  std::array<std::vector<float>, arraySize> occFV0AUnfm80;
  std::array<std::vector<float>, arraySize> occFV0CUnfm80;
  std::array<std::vector<float>, arraySize> occFT0AUnfm80;
  std::array<std::vector<float>, arraySize> occFT0CUnfm80;
  std::array<std::vector<float>, arraySize> occFDDAUnfm80;
  std::array<std::vector<float>, arraySize> occFDDCUnfm80;

  std::array<std::vector<float>, arraySize> occNTrackITSUnfm80;
  std::array<std::vector<float>, arraySize> occNTrackTPCUnfm80;
  std::array<std::vector<float>, arraySize> occNTrackTRDUnfm80;
  std::array<std::vector<float>, arraySize> occNTrackTOFUnfm80;
  std::array<std::vector<float>, arraySize> occNTrackSizeUnfm80;
  std::array<std::vector<float>, arraySize> occNTrackTPCAUnfm80;
  std::array<std::vector<float>, arraySize> occNTrackTPCCUnfm80;
  std::array<std::vector<float>, arraySize> occNTrackITSTPCUnfm80;
  std::array<std::vector<float>, arraySize> occNTrackITSTPCAUnfm80;
  std::array<std::vector<float>, arraySize> occNTrackITSTPCCUnfm80;

  std::array<std::vector<float>, arraySize> occMultNTracksHasITSUnfm80;
  std::array<std::vector<float>, arraySize> occMultNTracksHasTPCUnfm80;
  std::array<std::vector<float>, arraySize> occMultNTracksHasTOFUnfm80;
  std::array<std::vector<float>, arraySize> occMultNTracksHasTRDUnfm80;
  std::array<std::vector<float>, arraySize> occMultNTracksITSOnlyUnfm80;
  std::array<std::vector<float>, arraySize> occMultNTracksTPCOnlyUnfm80;
  std::array<std::vector<float>, arraySize> occMultNTracksITSTPCUnfm80;
  std::array<std::vector<float>, arraySize> occMultAllTracksTPCOnlyUnfm80;

  std::vector<float> vecRobustOccT0V0PrimUnfm80;
  std::vector<float> vecRobustOccFDDT0V0PrimUnfm80;
  std::vector<float> vecRobustOccNtrackDetUnfm80;
  std::vector<float> vecRobustOccmultTableUnfm80;
  std::vector<std::array<int, 2>> vecRobustOccT0V0PrimUnfm80medianPosVec; // Median => one for odd and two for even entries
  std::vector<std::array<int, 2>> vecRobustOccFDDT0V0PrimUnfm80medianPosVec;
  std::vector<std::array<int, 2>> vecRobustOccNtrackDetUnfm80medianPosVec;
  std::vector<std::array<int, 2>> vecRobustOccmultTableUnfm80medianPosVec;

  void init(InitContext const&)
  {
    // Set size of the vectors
    for (int i = 0; i < arraySize; i++) {
      bcTFMap[i].resize(nBCinTF / 80);
      occPrimUnfm80[i].resize(nBCinTF / 80);
      occFV0AUnfm80[i].resize(nBCinTF / 80);
      occFV0CUnfm80[i].resize(nBCinTF / 80);
      occFT0AUnfm80[i].resize(nBCinTF / 80);
      occFT0CUnfm80[i].resize(nBCinTF / 80);
      occFDDAUnfm80[i].resize(nBCinTF / 80);
      occFDDCUnfm80[i].resize(nBCinTF / 80);

      occNTrackITSUnfm80[i].resize(nBCinTF / 80);
      occNTrackTPCUnfm80[i].resize(nBCinTF / 80);
      occNTrackTRDUnfm80[i].resize(nBCinTF / 80);
      occNTrackTOFUnfm80[i].resize(nBCinTF / 80);
      occNTrackSizeUnfm80[i].resize(nBCinTF / 80);
      occNTrackTPCAUnfm80[i].resize(nBCinTF / 80);
      occNTrackTPCCUnfm80[i].resize(nBCinTF / 80);
      occNTrackITSTPCUnfm80[i].resize(nBCinTF / 80);
      occNTrackITSTPCAUnfm80[i].resize(nBCinTF / 80);
      occNTrackITSTPCCUnfm80[i].resize(nBCinTF / 80);

      occMultNTracksHasITSUnfm80[i].resize(nBCinTF / 80);
      occMultNTracksHasTPCUnfm80[i].resize(nBCinTF / 80);
      occMultNTracksHasTOFUnfm80[i].resize(nBCinTF / 80);
      occMultNTracksHasTRDUnfm80[i].resize(nBCinTF / 80);
      occMultNTracksITSOnlyUnfm80[i].resize(nBCinTF / 80);
      occMultNTracksTPCOnlyUnfm80[i].resize(nBCinTF / 80);
      occMultNTracksITSTPCUnfm80[i].resize(nBCinTF / 80);
      occMultAllTracksTPCOnlyUnfm80[i].resize(nBCinTF / 80);
    }

    vecRobustOccT0V0PrimUnfm80.resize(nBCinTF / 80);
    vecRobustOccFDDT0V0PrimUnfm80.resize(nBCinTF / 80);
    vecRobustOccNtrackDetUnfm80.resize(nBCinTF / 80);
    vecRobustOccmultTableUnfm80.resize(nBCinTF / 80);

    vecRobustOccT0V0PrimUnfm80medianPosVec.resize(nBCinTF / 80); // Median => one for odd and two for even entries
    vecRobustOccFDDT0V0PrimUnfm80medianPosVec.resize(nBCinTF / 80);
    vecRobustOccNtrackDetUnfm80medianPosVec.resize(nBCinTF / 80);
    vecRobustOccmultTableUnfm80medianPosVec.resize(nBCinTF / 80);

    // Getting Info from CCDB, to be implemented Later
    recoEvent.add("h_nBCinTF", "h_nBCinTF(to check nBCinTF)", {HistType::kTH1F, {{100, 114040, 114060}}}); // 114048
    recoEvent.add("h_bcInTF", "h_bcInTF", {HistType::kTH1F, {{2000, 0, 200000}}});
    recoEvent.add("h_RO_T0V0PrimUnfm80", "h_RO_T0V0PrimUnfm80:median contributors", {HistType::kTH1F, {{12 * 2, -1, 11}}});
    recoEvent.add("h_RO_FDDT0V0PrimUnfm80", "h_RO_FDDT0V0PrimUnfm80:median contributors", {HistType::kTH1F, {{12 * 2, -1, 11}}});
    recoEvent.add("h_RO_NtrackDetUnfm80", "h_RO_NtrackDetITS/TPC/TRD/TOF_80:median contributors", {HistType::kTH1F, {{12 * 2, -1, 11}}});
    recoEvent.add("h_RO_multTableUnfm80", "h_RO_multTableExtra_80:median contributors", {HistType::kTH1F, {{12 * 2, -1, 11}}});

    recoEvent.print();
  }

  void normalizeVector(std::vector<float>& OriginalVec, const float& scaleFactor)
  {
    std::transform(OriginalVec.begin(), OriginalVec.end(), OriginalVec.begin(), [scaleFactor](float x) { return x * scaleFactor; });
  }

  template <typename... Vecs>
  void getMedianOccVect(
    std::vector<float>& medianVector,
    std::vector<std::array<int, 2>>& medianPosVec,
    const Vecs&... vectors)
  {
    const int n = sizeof...(Vecs);                             // Number of vectors
    const int size = std::get<0>(std::tie(vectors...)).size(); // Size of the first vector

    for (int i = 0; i < size; i++) {
      std::vector<std::array<double, 2>> data; // first element is entry, second is index
      int iEntry = 0;

      // Lambda to iterate over all vectors
      auto collect = [&](const auto& vec) {
        data.push_back({vec[i], static_cast<double>(iEntry)});
        iEntry++;
      };
      (collect(vectors), ...); // Unpack variadic arguments and apply lambda

      // Sort the data
      std::sort(data.begin(), data.end(), [](const std::array<double, 2>& a, const std::array<double, 2>& b) {
        return a[0] < b[0];
      });

      double median;
      // Find the median
      if (n % 2 == 0) {
        median = (data[(n - 1) / 2][0] + data[(n - 1) / 2 + 1][0]) / 2;
        medianPosVec[i][0] = static_cast<int>(data[(n - 1) / 2][1] + 0.001);
        medianPosVec[i][1] = static_cast<int>(data[(n - 1) / 2 + 1][1] + 0.001);
      } else {
        median = data[n / 2][0];
        medianPosVec[i][0] = static_cast<int>(data[n / 2][1] + 0.001);
        medianPosVec[i][1] = -10; // For odd entries, only one value can be the median
      }
      medianVector[i] = median;
    }
  }

  void getRunInfo(const int& run, int& nBCsPerTF, int64_t& bcSOR)
  {
    auto runDuration = ccdb->getRunDuration(run, true);
    int64_t tsSOR = runDuration.first;
    auto ctpx = ccdb->getForTimeStamp<std::vector<int64_t>>("CTP/Calib/OrbitReset", tsSOR);
    int64_t tsOrbitReset = (*ctpx)[0];
    uint32_t nOrbitsPerTF = run < 534133 ? 128 : 32;
    int64_t orbitSOR = (tsSOR * 1000 - tsOrbitReset) / o2::constants::lhc::LHCOrbitMUS;
    orbitSOR = orbitSOR / nOrbitsPerTF * nOrbitsPerTF;
    bcSOR = orbitSOR * nBCsPerOrbit + customOrbitOffset * nBCsPerOrbit; // customOrbitOffset is a configurable
    nBCsPerTF = nOrbitsPerTF * nBCsPerOrbit;
  }

  template <typename T>
  void getTimingInfo(const T& bc, int& lastRun, int32_t& nBCsPerTF, int64_t& bcSOR, uint64_t& time, int64_t& tfIdThis, int& bcInTF)
  {
    int run = bc.runNumber();
    if (run != lastRun) { // update run info
      lastRun = run;
      getRunInfo(run, nBCsPerTF, bcSOR); // update nBCsPerTF && bcSOR
    }
    // update the information
    time = bc.timestamp();
    tfIdThis = (bc.globalBC() - bcSOR) / nBCsPerTF;
    bcInTF = (bc.globalBC() - bcSOR) % nBCsPerTF;
  }

  using MyCollisions = soa::Join<aod::Collisions, aod::Mults, aod::MultsExtra>;
  using MyTracks = soa::Join<aod::Tracks, o2::aod::TracksCov, aod::TracksExtra>;

  Preslice<MyTracks> tracksPerCollisionPreslice = o2::aod::track::collisionId;

  int dfCount = 0;
  int32_t nBCsPerTF = -999;
  int64_t bcSOR = -999;
  uint64_t time = -1;
  int64_t tfIdThis = -1;
  int bcInTF = -1;
  uint tfCounted = 0;
  int tfIDX = 0;
  int lastRun = -999;

  // Process the Data
  void process(o2::aod::BCsWithTimestamps const& BCs, MyCollisions const& collisions, MyTracks const& tracks) // aod::TracksQA const& tracksQA, o2::aod::Origins const& Origins //tables only used during debugging
  {
    // dfCount++;LOG(info) << "DEBUG 1 :: df_" << dfCount ;//<< " :: DF_" << Origins.begin().dataframeID() << " :: collisions.size() = " << collisions.size() << " :: tracks.size() = " << tracks.size() << " :: tracksQA.size() = " << tracksQA.size() << " :: BCs.size() = " << BCs.size();

    if (collisions.size() == 0) {
      for (const auto& BC : BCs) { // For BCs and OccIndexTable to have same size for joining
        getTimingInfo(BC, lastRun, nBCsPerTF, bcSOR, time, tfIdThis, bcInTF);
        genOccIndexTable(BC.globalIndex(), -999); // BCId, OccId
        genBCTFinfoTable(tfIdThis, bcInTF);
      }
      return;
    }

    // Initialisze the vectors components to zero
    tfIDX = 0;
    tfCounted = 0;
    for (int i = 0; i < arraySize; i++) {
      tfList[i] = -1;
      bcTFMap[i].clear(); // list of BCs used in one time frame;
      std::fill(occPrimUnfm80[i].begin(), occPrimUnfm80[i].end(), 0.);
      std::fill(occFV0AUnfm80[i].begin(), occFV0AUnfm80[i].end(), 0.);
      std::fill(occFV0CUnfm80[i].begin(), occFV0CUnfm80[i].end(), 0.);
      std::fill(occFT0AUnfm80[i].begin(), occFT0AUnfm80[i].end(), 0.);
      std::fill(occFT0CUnfm80[i].begin(), occFT0CUnfm80[i].end(), 0.);
      std::fill(occFDDAUnfm80[i].begin(), occFDDAUnfm80[i].end(), 0.);
      std::fill(occFDDCUnfm80[i].begin(), occFDDCUnfm80[i].end(), 0.);

      std::fill(occNTrackITSUnfm80[i].begin(), occNTrackITSUnfm80[i].end(), 0.);
      std::fill(occNTrackTPCUnfm80[i].begin(), occNTrackTPCUnfm80[i].end(), 0.);
      std::fill(occNTrackTRDUnfm80[i].begin(), occNTrackTRDUnfm80[i].end(), 0.);
      std::fill(occNTrackTOFUnfm80[i].begin(), occNTrackTOFUnfm80[i].end(), 0.);
      std::fill(occNTrackSizeUnfm80[i].begin(), occNTrackSizeUnfm80[i].end(), 0.);
      std::fill(occNTrackTPCAUnfm80[i].begin(), occNTrackTPCAUnfm80[i].end(), 0.);
      std::fill(occNTrackTPCCUnfm80[i].begin(), occNTrackTPCCUnfm80[i].end(), 0.);
      std::fill(occNTrackITSTPCUnfm80[i].begin(), occNTrackITSTPCUnfm80[i].end(), 0.);
      std::fill(occNTrackITSTPCAUnfm80[i].begin(), occNTrackITSTPCAUnfm80[i].end(), 0.);
      std::fill(occNTrackITSTPCCUnfm80[i].begin(), occNTrackITSTPCCUnfm80[i].end(), 0.);

      std::fill(occMultNTracksHasITSUnfm80[i].begin(), occMultNTracksHasITSUnfm80[i].end(), 0.);
      std::fill(occMultNTracksHasTPCUnfm80[i].begin(), occMultNTracksHasTPCUnfm80[i].end(), 0.);
      std::fill(occMultNTracksHasTOFUnfm80[i].begin(), occMultNTracksHasTOFUnfm80[i].end(), 0.);
      std::fill(occMultNTracksHasTRDUnfm80[i].begin(), occMultNTracksHasTRDUnfm80[i].end(), 0.);
      std::fill(occMultNTracksITSOnlyUnfm80[i].begin(), occMultNTracksITSOnlyUnfm80[i].end(), 0.);
      std::fill(occMultNTracksTPCOnlyUnfm80[i].begin(), occMultNTracksTPCOnlyUnfm80[i].end(), 0.);
      std::fill(occMultNTracksITSTPCUnfm80[i].begin(), occMultNTracksITSTPCUnfm80[i].end(), 0.);
      std::fill(occMultAllTracksTPCOnlyUnfm80[i].begin(), occMultAllTracksTPCOnlyUnfm80[i].end(), 0.);
    }

    std::vector<int64_t> tfIDList;
    for (const auto& collision : collisions) {
      const auto& bc = collision.bc_as<aod::BCsWithTimestamps>();
      getTimingInfo(bc, lastRun, nBCsPerTF, bcSOR, time, tfIdThis, bcInTF);

      recoEvent.fill(HIST("h_nBCinTF"), nBCsPerTF);
      recoEvent.fill(HIST("h_bcInTF"), bcInTF);

      if (nBCsPerTF > nBCinTF) {
        LOG(error) << "DEBUG :: FATAL ERROR :: nBCsPerTF > nBCinTF i.e " << nBCsPerTF << " > " << nBCinTF << " will cause crash in further process";
        return;
      }

      const uint64_t collIdx = collision.globalIndex();
      const auto tracksTablePerColl = tracks.sliceBy(tracksPerCollisionPreslice, collIdx);

      int nTrackITS = 0;
      int nTrackTPC = 0;
      int nTrackTRD = 0;
      int nTrackTOF = 0;
      int nTrackTPCA = 0;
      int nTrackTPCC = 0;
      int nTrackITSTPCA = 0;
      int nTrackITSTPCC = 0;

      for (const auto& track : tracksTablePerColl) {
        if (track.hasITS()) {
          nTrackITS++;
        } // Flag to check if track has a ITS match
        if (track.hasTPC()) {
          nTrackTPC++;
          if (track.eta() <= 0.0) {
            nTrackTPCA++; // includes tracks at eta zero as well.
          } else {
            nTrackTPCC++;
          }
        } // Flag to check if track has a TPC match
        if (track.hasTRD()) {
          nTrackTRD++;
        } // Flag to check if track has a TRD match
        if (track.hasTOF()) {
          nTrackTOF++;
        } // Flag to check if track has a TOF measurement
        if (track.hasITS() && track.hasTPC()) {
          if (track.eta() <= 0.0) {
            nTrackITSTPCA++; // includes tracks at eta zero as well.
          } else {
            nTrackITSTPCC++;
          }
        }
      } // track loop

      // if (collision.multNTracksTPCOnly() != 0) {
      //   LOG(error) << "DEBUG :: ERROR = multNTracksTPCOnly != 0" << collision.multNTracksTPCOnly();
      //   return;
      // }
      // if (collision.multAllTracksITSTPC() != nTrackITSTPC) {
      //   LOG(error) << "DEBUG :: ERROR :: 10 multAllTracksITSTPC :: " << collision.multAllTracksITSTPC() << " != " << nTrackITSTPC;
      //   return;
      // }

      tfIDList.push_back(tfIdThis);

      if (tfList[tfIDX] != tfIdThis) {
        if (tfCounted != 0) {
          tfIDX++;
        } //
        tfList[tfIDX] = tfIdThis;
        tfCounted++;
      }

      bcTFMap[tfIDX].push_back(bc.globalIndex());
      auto& tfOccPrimUnfm80 = occPrimUnfm80[tfIDX];
      auto& tfOccFV0AUnfm80 = occFV0AUnfm80[tfIDX];
      auto& tfOccFV0CUnfm80 = occFV0CUnfm80[tfIDX];
      auto& tfOccFT0AUnfm80 = occFT0AUnfm80[tfIDX];
      auto& tfOccFT0CUnfm80 = occFT0CUnfm80[tfIDX];
      auto& tfOccFDDAUnfm80 = occFDDAUnfm80[tfIDX];
      auto& tfOccFDDCUnfm80 = occFDDCUnfm80[tfIDX];

      auto& tfOccNTrackITSUnfm80 = occNTrackITSUnfm80[tfIDX];
      auto& tfOccNTrackTPCUnfm80 = occNTrackTPCUnfm80[tfIDX];
      auto& tfOccNTrackTRDUnfm80 = occNTrackTRDUnfm80[tfIDX];
      auto& tfOccNTrackTOFUnfm80 = occNTrackTOFUnfm80[tfIDX];
      auto& tfOccNTrackSizeUnfm80 = occNTrackSizeUnfm80[tfIDX];
      auto& tfOccNTrackTPCAUnfm80 = occNTrackTPCAUnfm80[tfIDX];
      auto& tfOccNTrackTPCCUnfm80 = occNTrackTPCCUnfm80[tfIDX];
      auto& tfOccNTrackITSTPCUnfm80 = occNTrackITSTPCUnfm80[tfIDX];
      auto& tfOccNTrackITSTPCAUnfm80 = occNTrackITSTPCAUnfm80[tfIDX];
      auto& tfOccNTrackITSTPCCUnfm80 = occNTrackITSTPCCUnfm80[tfIDX];

      auto& tfOccMultNTracksHasITSUnfm80 = occMultNTracksHasITSUnfm80[tfIDX];
      auto& tfOccMultNTracksHasTPCUnfm80 = occMultNTracksHasTPCUnfm80[tfIDX];
      auto& tfOccMultNTracksHasTOFUnfm80 = occMultNTracksHasTOFUnfm80[tfIDX];
      auto& tfOccMultNTracksHasTRDUnfm80 = occMultNTracksHasTRDUnfm80[tfIDX];
      auto& tfOccMultNTracksITSOnlyUnfm80 = occMultNTracksITSOnlyUnfm80[tfIDX];
      auto& tfOccMultNTracksTPCOnlyUnfm80 = occMultNTracksTPCOnlyUnfm80[tfIDX];
      auto& tfOccMultNTracksITSTPCUnfm80 = occMultNTracksITSTPCUnfm80[tfIDX];
      auto& tfOccMultAllTracksTPCOnlyUnfm80 = occMultAllTracksTPCOnlyUnfm80[tfIDX];

      // current collision bin in 80/160 grouping.
      int bin80Zero = bcInTF / 80;
      // int bin160_0=bcInTF/160;

      // float fbin80Zero =float(bcInTF)/80;
      // float fbin160_0=float(bcInTF)/160;

      ushort fNumContrib = collision.numContrib();
      float fMultFV0A = collision.multFV0A(), fMultFV0C = collision.multFV0C();
      float fMultFT0A = collision.multFT0A(), fMultFT0C = collision.multFT0C();
      float fMultFDDA = collision.multFDDA(), fMultFDDC = collision.multFDDC();
      int fNTrackITS = nTrackITS;
      int fNTrackTPC = nTrackTPC, fNTrackTRD = nTrackTRD;
      int fNTrackTOF = nTrackTOF;
      int fNTrackTPCA = nTrackTPCA, fNTrackTPCC = nTrackTPCC;
      int fNTrackITSTPC = collision.multAllTracksITSTPC(), fNTrackSize = tracksTablePerColl.size();
      int fNTrackITSTPCA = nTrackITSTPCA, fNTrackITSTPCC = nTrackITSTPCC;

      // Processing for grouping of 80 BCs
      for (int deltaBin = 0; deltaBin < nBCinDrift / 80; deltaBin++) {
        tfOccPrimUnfm80[(bin80Zero + deltaBin) % (nBCinTF / 80)] += fNumContrib * 1;
        tfOccFV0AUnfm80[(bin80Zero + deltaBin) % (nBCinTF / 80)] += fMultFV0A * 1;
        tfOccFV0CUnfm80[(bin80Zero + deltaBin) % (nBCinTF / 80)] += fMultFV0C * 1;
        tfOccFT0AUnfm80[(bin80Zero + deltaBin) % (nBCinTF / 80)] += fMultFT0A * 1;
        tfOccFT0CUnfm80[(bin80Zero + deltaBin) % (nBCinTF / 80)] += fMultFT0C * 1;
        tfOccFDDAUnfm80[(bin80Zero + deltaBin) % (nBCinTF / 80)] += fMultFDDA * 1;
        tfOccFDDCUnfm80[(bin80Zero + deltaBin) % (nBCinTF / 80)] += fMultFDDC * 1;

        tfOccNTrackITSUnfm80[(bin80Zero + deltaBin) % (nBCinTF / 80)] += fNTrackITS * 1;
        tfOccNTrackTPCUnfm80[(bin80Zero + deltaBin) % (nBCinTF / 80)] += fNTrackTPC * 1;
        tfOccNTrackTRDUnfm80[(bin80Zero + deltaBin) % (nBCinTF / 80)] += fNTrackTRD * 1;
        tfOccNTrackTOFUnfm80[(bin80Zero + deltaBin) % (nBCinTF / 80)] += fNTrackTOF * 1;
        tfOccNTrackSizeUnfm80[(bin80Zero + deltaBin) % (nBCinTF / 80)] += fNTrackSize * 1;
        tfOccNTrackTPCAUnfm80[(bin80Zero + deltaBin) % (nBCinTF / 80)] += fNTrackTPCA * 1;
        tfOccNTrackTPCCUnfm80[(bin80Zero + deltaBin) % (nBCinTF / 80)] += fNTrackTPCC * 1;
        tfOccNTrackITSTPCUnfm80[(bin80Zero + deltaBin) % (nBCinTF / 80)] += fNTrackITSTPC * 1;
        tfOccNTrackITSTPCAUnfm80[(bin80Zero + deltaBin) % (nBCinTF / 80)] += fNTrackITSTPCA * 1;
        tfOccNTrackITSTPCCUnfm80[(bin80Zero + deltaBin) % (nBCinTF / 80)] += fNTrackITSTPCC * 1;

        tfOccMultNTracksHasITSUnfm80[(bin80Zero + deltaBin) % (nBCinTF / 80)] += collision.multNTracksHasITS() * 1;
        tfOccMultNTracksHasTPCUnfm80[(bin80Zero + deltaBin) % (nBCinTF / 80)] += collision.multNTracksHasTPC() * 1;
        tfOccMultNTracksHasTOFUnfm80[(bin80Zero + deltaBin) % (nBCinTF / 80)] += collision.multNTracksHasTOF() * 1;
        tfOccMultNTracksHasTRDUnfm80[(bin80Zero + deltaBin) % (nBCinTF / 80)] += collision.multNTracksHasTRD() * 1;
        tfOccMultNTracksITSOnlyUnfm80[(bin80Zero + deltaBin) % (nBCinTF / 80)] += collision.multNTracksITSOnly() * 1;
        tfOccMultNTracksTPCOnlyUnfm80[(bin80Zero + deltaBin) % (nBCinTF / 80)] += collision.multNTracksTPCOnly() * 1;
        tfOccMultNTracksITSTPCUnfm80[(bin80Zero + deltaBin) % (nBCinTF / 80)] += collision.multNTracksITSTPC() * 1;
        tfOccMultAllTracksTPCOnlyUnfm80[(bin80Zero + deltaBin) % (nBCinTF / 80)] += collision.multAllTracksTPCOnly() * 1;
      }
    }
    // collision Loop is over

    std::vector<int64_t> sortedTfIDList = tfIDList;
    std::sort(sortedTfIDList.begin(), sortedTfIDList.end());
    auto last = std::unique(sortedTfIDList.begin(), sortedTfIDList.end());
    sortedTfIDList.erase(last, sortedTfIDList.end());

    if (tfCounted != sortedTfIDList.size()) {
      LOG(error) << "DEBUG :: Number mismatch for tf counted and filled :: " << tfCounted << " != " << sortedTfIDList.size();
    }

    int totalBCcountSize = 0;
    for (int i = 0; i < arraySize; i++) {
      totalBCcountSize += bcTFMap[i].size();
      // check if the BCs are already sorted or not
      if (!std::is_sorted(bcTFMap[i].begin(), bcTFMap[i].end())) {
        LOG(debug) << "DEBUG :: ERROR :: BCs are not sorted";
      }
    }
    //
    if (totalBCcountSize != collisions.size()) {
      LOG(debug) << "DEBUG :: ERROR :: filled TF list and collision size mismatch ::  filledTF_Size = " << totalBCcountSize << " != " << collisions.size() << " = collisions.size()";
    }

    // Fill the Producers
    for (uint i = 0; i < tfCounted; i++) {

      auto& vecOccPrimUnfm80 = occPrimUnfm80[i];
      auto& vecOccFV0AUnfm80 = occFV0AUnfm80[i];
      auto& vecOccFV0CUnfm80 = occFV0CUnfm80[i];
      auto& vecOccFT0AUnfm80 = occFT0AUnfm80[i];
      auto& vecOccFT0CUnfm80 = occFT0CUnfm80[i];
      auto& vecOccFDDAUnfm80 = occFDDAUnfm80[i];
      auto& vecOccFDDCUnfm80 = occFDDCUnfm80[i];

      auto& vecOccNTrackITSUnfm80 = occNTrackITSUnfm80[i];
      auto& vecOccNTrackTPCUnfm80 = occNTrackTPCUnfm80[i];
      auto& vecOccNTrackTRDUnfm80 = occNTrackTRDUnfm80[i];
      auto& vecOccNTrackTOFUnfm80 = occNTrackTOFUnfm80[i];
      auto& vecOccNTrackSizeUnfm80 = occNTrackSizeUnfm80[i];
      auto& vecOccNTrackTPCAUnfm80 = occNTrackTPCAUnfm80[i];
      auto& vecOccNTrackTPCCUnfm80 = occNTrackTPCCUnfm80[i];
      auto& vecOccNTrackITSTPCUnfm80 = occNTrackITSTPCUnfm80[i];
      auto& vecOccNTrackITSTPCAUnfm80 = occNTrackITSTPCAUnfm80[i];
      auto& vecOccNTrackITSTPCCUnfm80 = occNTrackITSTPCCUnfm80[i];

      auto& vecOccMultNTracksHasITSUnfm80 = occMultNTracksHasITSUnfm80[i];
      auto& vecOccMultNTracksHasTPCUnfm80 = occMultNTracksHasTPCUnfm80[i];
      auto& vecOccMultNTracksHasTOFUnfm80 = occMultNTracksHasTOFUnfm80[i];
      auto& vecOccMultNTracksHasTRDUnfm80 = occMultNTracksHasTRDUnfm80[i];
      auto& vecOccMultNTracksITSOnlyUnfm80 = occMultNTracksITSOnlyUnfm80[i];
      auto& vecOccMultNTracksTPCOnlyUnfm80 = occMultNTracksTPCOnlyUnfm80[i];
      auto& vecOccMultNTracksITSTPCUnfm80 = occMultNTracksITSTPCUnfm80[i];
      auto& vecOccMultAllTracksTPCOnlyUnfm80 = occMultAllTracksTPCOnlyUnfm80[i];

      float meanOccPrimUnfm80 = TMath::Mean(vecOccPrimUnfm80.size(), vecOccPrimUnfm80.data());
      float meanOccFV0AUnfm80 = TMath::Mean(vecOccFV0AUnfm80.size(), vecOccFV0AUnfm80.data());
      float meanOccFV0CUnfm80 = TMath::Mean(vecOccFV0CUnfm80.size(), vecOccFV0CUnfm80.data());
      float meanOccFT0AUnfm80 = TMath::Mean(vecOccFT0AUnfm80.size(), vecOccFT0AUnfm80.data());
      float meanOccFT0CUnfm80 = TMath::Mean(vecOccFT0CUnfm80.size(), vecOccFT0CUnfm80.data());
      float meanOccFDDAUnfm80 = TMath::Mean(vecOccFDDAUnfm80.size(), vecOccFDDAUnfm80.data());
      float meanOccFDDCUnfm80 = TMath::Mean(vecOccFDDCUnfm80.size(), vecOccFDDCUnfm80.data());

      float meanOccNTrackITSUnfm80 = TMath::Mean(vecOccNTrackITSUnfm80.size(), vecOccNTrackITSUnfm80.data());
      float meanOccNTrackTPCUnfm80 = TMath::Mean(vecOccNTrackTPCUnfm80.size(), vecOccNTrackTPCUnfm80.data());
      float meanOccNTrackTRDUnfm80 = TMath::Mean(vecOccNTrackTRDUnfm80.size(), vecOccNTrackTRDUnfm80.data());
      float meanOccNTrackTOFUnfm80 = TMath::Mean(vecOccNTrackTOFUnfm80.size(), vecOccNTrackTOFUnfm80.data());
      float meanOccNTrackSizeUnfm80 = TMath::Mean(vecOccNTrackSizeUnfm80.size(), vecOccNTrackSizeUnfm80.data());
      float meanOccNTrackTPCAUnfm80 = TMath::Mean(vecOccNTrackTPCAUnfm80.size(), vecOccNTrackTPCAUnfm80.data());
      float meanOccNTrackTPCCUnfm80 = TMath::Mean(vecOccNTrackTPCCUnfm80.size(), vecOccNTrackTPCCUnfm80.data());
      float meanOccNTrackITSTPCUnfm80 = TMath::Mean(vecOccNTrackITSTPCUnfm80.size(), vecOccNTrackITSTPCUnfm80.data());
      float meanOccNTrackITSTPCAUnfm80 = TMath::Mean(vecOccNTrackITSTPCAUnfm80.size(), vecOccNTrackITSTPCAUnfm80.data());
      float meanOccNTrackITSTPCCUnfm80 = TMath::Mean(vecOccNTrackITSTPCCUnfm80.size(), vecOccNTrackITSTPCCUnfm80.data());

      float meanOccMultNTracksHasITSUnfm80 = TMath::Mean(vecOccMultNTracksHasITSUnfm80.size(), vecOccMultNTracksHasITSUnfm80.data());
      float meanOccMultNTracksHasTPCUnfm80 = TMath::Mean(vecOccMultNTracksHasTPCUnfm80.size(), vecOccMultNTracksHasTPCUnfm80.data());
      float meanOccMultNTracksHasTOFUnfm80 = TMath::Mean(vecOccMultNTracksHasTOFUnfm80.size(), vecOccMultNTracksHasTOFUnfm80.data());
      float meanOccMultNTracksHasTRDUnfm80 = TMath::Mean(vecOccMultNTracksHasTRDUnfm80.size(), vecOccMultNTracksHasTRDUnfm80.data());
      float meanOccMultNTracksITSOnlyUnfm80 = TMath::Mean(vecOccMultNTracksITSOnlyUnfm80.size(), vecOccMultNTracksITSOnlyUnfm80.data());
      float meanOccMultNTracksTPCOnlyUnfm80 = TMath::Mean(vecOccMultNTracksTPCOnlyUnfm80.size(), vecOccMultNTracksTPCOnlyUnfm80.data());
      float meanOccMultNTracksITSTPCUnfm80 = TMath::Mean(vecOccMultNTracksITSTPCUnfm80.size(), vecOccMultNTracksITSTPCUnfm80.data());
      float meanOccMultAllTracksTPCOnlyUnfm80 = TMath::Mean(vecOccMultAllTracksTPCOnlyUnfm80.size(), vecOccMultAllTracksTPCOnlyUnfm80.data());

      // Normalise the original vectors
      normalizeVector(vecOccPrimUnfm80, meanOccPrimUnfm80 / meanOccPrimUnfm80);
      normalizeVector(vecOccFV0AUnfm80, meanOccPrimUnfm80 / meanOccFV0AUnfm80);
      normalizeVector(vecOccFV0CUnfm80, meanOccPrimUnfm80 / meanOccFV0CUnfm80);
      normalizeVector(vecOccFT0AUnfm80, meanOccPrimUnfm80 / meanOccFT0AUnfm80);
      normalizeVector(vecOccFT0CUnfm80, meanOccPrimUnfm80 / meanOccFT0CUnfm80);
      normalizeVector(vecOccFDDAUnfm80, meanOccPrimUnfm80 / meanOccFDDAUnfm80);
      normalizeVector(vecOccFDDCUnfm80, meanOccPrimUnfm80 / meanOccFDDCUnfm80);

      normalizeVector(vecOccNTrackITSUnfm80, meanOccPrimUnfm80 / meanOccNTrackITSUnfm80);
      normalizeVector(vecOccNTrackTPCUnfm80, meanOccPrimUnfm80 / meanOccNTrackTPCUnfm80);
      normalizeVector(vecOccNTrackTRDUnfm80, meanOccPrimUnfm80 / meanOccNTrackTRDUnfm80);
      normalizeVector(vecOccNTrackTOFUnfm80, meanOccPrimUnfm80 / meanOccNTrackTOFUnfm80);
      normalizeVector(vecOccNTrackSizeUnfm80, meanOccPrimUnfm80 / meanOccNTrackSizeUnfm80);
      normalizeVector(vecOccNTrackTPCAUnfm80, meanOccPrimUnfm80 / meanOccNTrackTPCAUnfm80);
      normalizeVector(vecOccNTrackTPCCUnfm80, meanOccPrimUnfm80 / meanOccNTrackTPCCUnfm80);
      normalizeVector(vecOccNTrackITSTPCUnfm80, meanOccPrimUnfm80 / meanOccNTrackITSTPCUnfm80);
      normalizeVector(vecOccNTrackITSTPCAUnfm80, meanOccPrimUnfm80 / meanOccNTrackITSTPCAUnfm80);
      normalizeVector(vecOccNTrackITSTPCCUnfm80, meanOccPrimUnfm80 / meanOccNTrackITSTPCCUnfm80);

      normalizeVector(vecOccMultNTracksHasITSUnfm80, meanOccPrimUnfm80 / meanOccMultNTracksHasITSUnfm80);
      normalizeVector(vecOccMultNTracksHasTPCUnfm80, meanOccPrimUnfm80 / meanOccMultNTracksHasTPCUnfm80);
      normalizeVector(vecOccMultNTracksHasTOFUnfm80, meanOccPrimUnfm80 / meanOccMultNTracksHasTOFUnfm80);
      normalizeVector(vecOccMultNTracksHasTRDUnfm80, meanOccPrimUnfm80 / meanOccMultNTracksHasTRDUnfm80);
      normalizeVector(vecOccMultNTracksITSOnlyUnfm80, meanOccPrimUnfm80 / meanOccMultNTracksITSOnlyUnfm80);
      normalizeVector(vecOccMultNTracksTPCOnlyUnfm80, meanOccPrimUnfm80 / meanOccMultNTracksTPCOnlyUnfm80);
      normalizeVector(vecOccMultNTracksITSTPCUnfm80, meanOccPrimUnfm80 / meanOccMultNTracksITSTPCUnfm80);
      normalizeVector(vecOccMultAllTracksTPCOnlyUnfm80, meanOccPrimUnfm80 / meanOccMultAllTracksTPCOnlyUnfm80);

      // Find Robust estimators
      // T0A, T0C, V0A, Prim
      getMedianOccVect(vecRobustOccT0V0PrimUnfm80, vecRobustOccT0V0PrimUnfm80medianPosVec,
                       vecOccPrimUnfm80, vecOccFV0AUnfm80, vecOccFT0AUnfm80, vecOccFT0CUnfm80);

      // T0A, T0C, V0A, FDD, Prim
      getMedianOccVect(vecRobustOccFDDT0V0PrimUnfm80, vecRobustOccFDDT0V0PrimUnfm80medianPosVec,
                       vecOccPrimUnfm80, vecOccFV0AUnfm80, vecOccFT0AUnfm80, vecOccFT0CUnfm80, vecOccFDDAUnfm80, vecOccFDDCUnfm80);

      // NTrackDet
      getMedianOccVect(vecRobustOccNtrackDetUnfm80, vecRobustOccNtrackDetUnfm80medianPosVec,
                       vecOccNTrackITSUnfm80, vecOccNTrackTPCUnfm80, vecOccNTrackTRDUnfm80, vecOccNTrackTOFUnfm80);

      // multExtraTable
      getMedianOccVect(vecRobustOccmultTableUnfm80, vecRobustOccmultTableUnfm80medianPosVec,
                       vecOccPrimUnfm80, vecOccMultNTracksHasITSUnfm80, vecOccMultNTracksHasTPCUnfm80,
                       vecOccMultNTracksHasTOFUnfm80, vecOccMultNTracksHasTRDUnfm80, vecOccMultNTracksITSOnlyUnfm80,
                       vecOccMultNTracksTPCOnlyUnfm80, vecOccMultNTracksITSTPCUnfm80, vecOccMultAllTracksTPCOnlyUnfm80,
                       vecOccNTrackITSTPCUnfm80);

      for (const auto& vec : vecRobustOccT0V0PrimUnfm80medianPosVec) {
        recoEvent.fill(HIST("h_RO_T0V0PrimUnfm80"), vec[0]);
        recoEvent.fill(HIST("h_RO_T0V0PrimUnfm80"), vec[1]);
      }

      for (const auto& vec : vecRobustOccFDDT0V0PrimUnfm80medianPosVec) {
        recoEvent.fill(HIST("h_RO_FDDT0V0PrimUnfm80"), vec[0]);
        recoEvent.fill(HIST("h_RO_FDDT0V0PrimUnfm80"), vec[1]);
      }

      for (const auto& vec : vecRobustOccNtrackDetUnfm80medianPosVec) {
        recoEvent.fill(HIST("h_RO_NtrackDetUnfm80"), vec[0]);
        recoEvent.fill(HIST("h_RO_NtrackDetUnfm80"), vec[1]);
      }

      for (const auto& vec : vecRobustOccmultTableUnfm80medianPosVec) {
        recoEvent.fill(HIST("h_RO_multTableUnfm80"), vec[0]);
        recoEvent.fill(HIST("h_RO_multTableUnfm80"), vec[1]);
      }

      genOccsBCsList(tfList[i], bcTFMap[i]);

      if (buildOccsDet) {
        genOccsDet(vecOccPrimUnfm80,
                   vecOccFV0AUnfm80, vecOccFV0CUnfm80,
                   vecOccFT0AUnfm80, vecOccFT0CUnfm80,
                   vecOccFDDAUnfm80, vecOccFDDCUnfm80);
      }
      if (buildOccsTrackMult) {
        genOccsTrackMult(vecOccNTrackITSUnfm80, vecOccNTrackTPCUnfm80,
                         vecOccNTrackTRDUnfm80, vecOccNTrackTOFUnfm80,
                         vecOccNTrackSizeUnfm80, vecOccNTrackTPCAUnfm80,
                         vecOccNTrackTPCCUnfm80, vecOccNTrackITSTPCUnfm80,
                         vecOccNTrackITSTPCAUnfm80, vecOccNTrackITSTPCCUnfm80);
      }
      if (buildOccsMultExtra) {
        genOccsMultExtra(vecOccMultNTracksHasITSUnfm80, vecOccMultNTracksHasTPCUnfm80,
                         vecOccMultNTracksHasTOFUnfm80, vecOccMultNTracksHasTRDUnfm80,
                         vecOccMultNTracksITSOnlyUnfm80, vecOccMultNTracksTPCOnlyUnfm80,
                         vecOccMultNTracksITSTPCUnfm80, vecOccMultAllTracksTPCOnlyUnfm80);
      }

      if (buildOccsRobust) {
        genOccsRobust(vecRobustOccT0V0PrimUnfm80, vecRobustOccFDDT0V0PrimUnfm80, vecRobustOccNtrackDetUnfm80, vecRobustOccmultTableUnfm80);
      }
      if (buildOccsMeanDet) {
        genOccsMeanDet(meanOccPrimUnfm80,
                       meanOccFV0AUnfm80, meanOccFV0CUnfm80,
                       meanOccFT0AUnfm80, meanOccFT0CUnfm80,
                       meanOccFDDAUnfm80, meanOccFDDCUnfm80);
      }

      if (buildOccsMeanTrkMult) {
        genOccsMeanTrkMult(meanOccNTrackITSUnfm80,
                           meanOccNTrackTPCUnfm80,
                           meanOccNTrackTRDUnfm80,
                           meanOccNTrackTOFUnfm80,
                           meanOccNTrackSizeUnfm80,
                           meanOccNTrackTPCAUnfm80,
                           meanOccNTrackTPCCUnfm80,
                           meanOccNTrackITSTPCUnfm80,
                           meanOccNTrackITSTPCAUnfm80,
                           meanOccNTrackITSTPCCUnfm80);
      }

      if (buildOccsMnMultExtra) {
        genOccsMnMultExtra(meanOccMultNTracksHasITSUnfm80, meanOccMultNTracksHasTPCUnfm80,
                           meanOccMultNTracksHasTOFUnfm80, meanOccMultNTracksHasTRDUnfm80,
                           meanOccMultNTracksITSOnlyUnfm80, meanOccMultNTracksTPCOnlyUnfm80,
                           meanOccMultNTracksITSTPCUnfm80, meanOccMultAllTracksTPCOnlyUnfm80);
      }

      if (buildOccsMeanRobust) {
        genOccsMeanRobust(
          TMath::Mean(vecRobustOccT0V0PrimUnfm80.size(), vecRobustOccT0V0PrimUnfm80.data()),
          TMath::Mean(vecRobustOccFDDT0V0PrimUnfm80.size(), vecRobustOccFDDT0V0PrimUnfm80.data()),
          TMath::Mean(vecRobustOccNtrackDetUnfm80.size(), vecRobustOccNtrackDetUnfm80.data()),
          TMath::Mean(vecRobustOccmultTableUnfm80.size(), vecRobustOccmultTableUnfm80.data()));
      }
    }

    // Create a BC index table.
    int64_t occIDX = -1;
    int idx = -1;
    for (auto const& bc : BCs) {
      idx = -1;
      getTimingInfo(bc, lastRun, nBCsPerTF, bcSOR, time, tfIdThis, bcInTF);

      auto idxIt = std::find(tfList.begin(), tfList.end(), tfIdThis);
      if (idxIt != tfList.end()) {
        idx = std::distance(tfList.begin(), idxIt);
      } else {
        LOG(error) << "DEBUG :: SEVERE :: BC  Timeframe not in the list";
      }

      auto it = std::find(bcTFMap[idx].begin(), bcTFMap[idx].end(), bc.globalIndex()); // will find the iterator where object is placed.
      if (it != bcTFMap[idx].end()) {
        occIDX = idx; // Element is in the vector
      } else {
        occIDX = -1; // Element is not in the vector
      }

      genOccIndexTable(bc.globalIndex(), occIDX); // BCId, OccId
      genBCTFinfoTable(tfIdThis, bcInTF);
    }
  } // Process function ends
};

struct TrackMeanOccTableProducer {
  Produces<aod::TrackMeanOccs0> genTrackMeanOccs0;
  Produces<aod::TrackMeanOccs1> genTrackMeanOccs1;
  Produces<aod::TrackMeanOccs2> genTrackMeanOccs2;
  Produces<aod::TrackMeanOccs3> genTrackMeanOccs3;
  Produces<aod::TrackMeanOccs4> genTrackMeanOccs4;
  Produces<aod::TrackMeanOccs5> genTrackMeanOccs5;
  Produces<aod::TrackMeanOccs6> genTrackMeanOccs6;
  Produces<aod::TrackMeanOccs7> genTrackMeanOccs7;
  Produces<aod::TrackMeanOccs8> genTrackMeanOccs8;

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  // Configurables
  Configurable<int> customOrbitOffset{"customOrbitOffset", 0, "customOrbitOffset for MC"};

  Configurable<bool> buildTrackMeanOccs1{"buildTrackMeanOccs1", true, "builder TrackMeanOccs1"};
  Configurable<bool> buildTrackMeanOccs2{"buildTrackMeanOccs2", true, "builder TrackMeanOccs2"};
  Configurable<bool> buildTrackMeanOccs3{"buildTrackMeanOccs3", true, "builder TrackMeanOccs3"};
  Configurable<bool> buildTrackMeanOccs4{"buildTrackMeanOccs4", true, "builder TrackMeanOccs4"};
  Configurable<bool> buildTrackMeanOccs5{"buildTrackMeanOccs5", true, "builder TrackMeanOccs5"};
  Configurable<bool> buildTrackMeanOccs6{"buildTrackMeanOccs6", true, "builder TrackMeanOccs6"};
  Configurable<bool> buildTrackMeanOccs7{"buildTrackMeanOccs7", true, "builder TrackMeanOccs7"};
  Configurable<bool> buildTrackMeanOccs8{"buildTrackMeanOccs8", true, "builder TrackMeanOccs8"};

  // vectors to be used for occupancy estimation
  std::vector<float> occPrimUnfm80;
  std::vector<float> occFV0AUnfm80;
  std::vector<float> occFV0CUnfm80;
  std::vector<float> occFT0AUnfm80;
  std::vector<float> occFT0CUnfm80;
  std::vector<float> occFDDAUnfm80;
  std::vector<float> occFDDCUnfm80;

  std::vector<float> occNTrackITSUnfm80;
  std::vector<float> occNTrackTPCUnfm80;
  std::vector<float> occNTrackTRDUnfm80;
  std::vector<float> occNTrackTOFUnfm80;
  std::vector<float> occNTrackSizeUnfm80;
  std::vector<float> occNTrackTPCAUnfm80;
  std::vector<float> occNTrackTPCCUnfm80;
  std::vector<float> occNTrackITSTPCUnfm80;
  std::vector<float> occNTrackITSTPCAUnfm80;
  std::vector<float> occNTrackITSTPCCUnfm80;

  std::vector<float> occMultNTracksHasITSUnfm80;
  std::vector<float> occMultNTracksHasTPCUnfm80;
  std::vector<float> occMultNTracksHasTOFUnfm80;
  std::vector<float> occMultNTracksHasTRDUnfm80;
  std::vector<float> occMultNTracksITSOnlyUnfm80;
  std::vector<float> occMultNTracksTPCOnlyUnfm80;
  std::vector<float> occMultNTracksITSTPCUnfm80;
  std::vector<float> occMultAllTracksTPCOnlyUnfm80;

  std::vector<float> occRobustT0V0PrimUnfm80;
  std::vector<float> occRobustFDDT0V0PrimUnfm80;
  std::vector<float> occRobustNtrackDetUnfm80;
  std::vector<float> occRobustMultTableUnfm80;

  void init(InitContext const&)
  {
    // CCDB related part to be added later
    occPrimUnfm80.resize(nBCinTF / 80);
    occFV0AUnfm80.resize(nBCinTF / 80);
    occFV0CUnfm80.resize(nBCinTF / 80);
    occFT0AUnfm80.resize(nBCinTF / 80);
    occFT0CUnfm80.resize(nBCinTF / 80);
    occFDDAUnfm80.resize(nBCinTF / 80);
    occFDDCUnfm80.resize(nBCinTF / 80);

    occNTrackITSUnfm80.resize(nBCinTF / 80);
    occNTrackTPCUnfm80.resize(nBCinTF / 80);
    occNTrackTRDUnfm80.resize(nBCinTF / 80);
    occNTrackTOFUnfm80.resize(nBCinTF / 80);
    occNTrackSizeUnfm80.resize(nBCinTF / 80);
    occNTrackTPCAUnfm80.resize(nBCinTF / 80);
    occNTrackTPCCUnfm80.resize(nBCinTF / 80);
    occNTrackITSTPCUnfm80.resize(nBCinTF / 80);
    occNTrackITSTPCAUnfm80.resize(nBCinTF / 80);
    occNTrackITSTPCCUnfm80.resize(nBCinTF / 80);

    occMultNTracksHasITSUnfm80.resize(nBCinTF / 80);
    occMultNTracksHasTPCUnfm80.resize(nBCinTF / 80);
    occMultNTracksHasTOFUnfm80.resize(nBCinTF / 80);
    occMultNTracksHasTRDUnfm80.resize(nBCinTF / 80);
    occMultNTracksITSOnlyUnfm80.resize(nBCinTF / 80);
    occMultNTracksTPCOnlyUnfm80.resize(nBCinTF / 80);
    occMultNTracksITSTPCUnfm80.resize(nBCinTF / 80);
    occMultAllTracksTPCOnlyUnfm80.resize(nBCinTF / 80);

    occRobustT0V0PrimUnfm80.resize(nBCinTF / 80);
    occRobustFDDT0V0PrimUnfm80.resize(nBCinTF / 80);
    occRobustNtrackDetUnfm80.resize(nBCinTF / 80);
    occRobustMultTableUnfm80.resize(nBCinTF / 80);
  }

  void getRunInfo(const int& run, int& nBCsPerTF, int64_t& bcSOR)
  {
    auto runDuration = ccdb->getRunDuration(run, true);
    int64_t tsSOR = runDuration.first;
    auto ctpx = ccdb->getForTimeStamp<std::vector<int64_t>>("CTP/Calib/OrbitReset", tsSOR);
    int64_t tsOrbitReset = (*ctpx)[0];
    uint32_t nOrbitsPerTF = run < 534133 ? 128 : 32;
    int64_t orbitSOR = (tsSOR * 1000 - tsOrbitReset) / o2::constants::lhc::LHCOrbitMUS;
    orbitSOR = orbitSOR / nOrbitsPerTF * nOrbitsPerTF;
    bcSOR = orbitSOR * nBCsPerOrbit + customOrbitOffset * nBCsPerOrbit; // customOrbitOffset is a configurable
    nBCsPerTF = nOrbitsPerTF * nBCsPerOrbit;
  }

  template <typename T>
  void getTimingInfo(const T& bc, int& lastRun, int32_t& nBCsPerTF, int64_t& bcSOR, uint64_t& time, int64_t& tfIdThis, int& bcInTF)
  {
    int run = bc.runNumber();
    if (run != lastRun) { // update run info
      lastRun = run;
      getRunInfo(run, nBCsPerTF, bcSOR); // update nBCsPerTF && bcSOR
    }
    // update the information
    time = bc.timestamp();
    tfIdThis = (bc.globalBC() - bcSOR) / nBCsPerTF;
    bcInTF = (bc.globalBC() - bcSOR) % nBCsPerTF;
  }

  float getMeanOccupancy(int bcBegin, int bcEnd, const std::vector<float>& OccVector)
  {
    float sumOfBins = 0;
    int binStart, binEnd;
    if (bcBegin <= bcEnd) {
      binStart = bcBegin;
      binEnd = bcEnd;
    } else {
      binStart = bcEnd;
      binEnd = bcBegin;
    }
    for (int i = binStart; i <= binEnd; i++) {
      sumOfBins += OccVector[i];
    }
    float meanOccupancy = sumOfBins / static_cast<double>(binEnd - binStart + 1);
    return meanOccupancy;
  }

  float getWeightedMeanOccupancy(int bcBegin, int bcEnd, const std::vector<float>& OccVector)
  {
    float sumOfBins = 0;
    int binStart, binEnd;
    // Assuming linear dependence of R on bins
    float m;      // slope of the equation
    float c;      // some constant in linear
    float x1, x2; //, y1 = 90., y2 = 245.;

    if (bcBegin <= bcEnd) {
      binStart = bcBegin;
      binEnd = bcEnd;
      x1 = static_cast<float>(binStart);
      x2 = static_cast<float>(binEnd);
    } else {
      binStart = bcEnd;
      binEnd = bcBegin;
      x1 = static_cast<float>(binEnd);
      x2 = static_cast<float>(binStart);
    } //

    if (x2 == x1) {
      m = 0;
    } else {
      m = (245. - 90.) / (x2 - x1);
    }
    c = 245. - m * x2;
    float weightSum = 0;
    float wr = 0;
    float r = 0;
    for (int i = binStart; i <= binEnd; i++) {
      r = m * i + c;
      wr = 125. / r;
      if (x2 == x1) {
        wr = 1.0;
      }
      sumOfBins += OccVector[i] * wr;
      weightSum += wr;
    }
    float meanOccupancy = sumOfBins / weightSum;
    return meanOccupancy;
  }

  using MyCollisions = soa::Join<aod::Collisions, aod::Mults>;
  using MyTracks = soa::Join<aod::Tracks, o2::aod::TracksCov, aod::TracksExtra>;
  using MyTracksQA = aod::TracksQA_000; // aod::TracksQAVersion; //aod::TracksQA
  using MyBCTable = soa::Join<aod::BCsWithTimestamps, aod::OccIndexTable, aod::BCTFinfoTable>;

  using MyOccsDet = soa::Join<aod::OccsBCsList, aod::OccsDet>;
  using MyOccsTrackMult = soa::Join<aod::OccsBCsList, aod::OccsTrackMult>;
  using MyOccsMultExtra = soa::Join<aod::OccsBCsList, aod::OccsMultExtra>;
  using MyOccsRobust = soa::Join<aod::OccsBCsList, aod::OccsRobust>;

  // Process the Data
  int dfCount = 0;
  int32_t nBCsPerTF = -999;
  int64_t bcSOR = -999;
  int lastRun = -999;

  uint64_t time = -1;
  int64_t tfIdThis = -1;
  int bcInTF = -1;

  void process(MyBCTable const& BCs, MyCollisions const& collisions, MyTracks const& tracks, MyTracksQA const& tracksQA,
               //  o2::aod::AmbiguousTracks const& ambgTracks, o2::aod::Origins const& Origins, aod::OccsBCsList const& occsBCsList, //tables only used during debugging
               MyOccsDet const& occsDet, MyOccsTrackMult const& occsTrackMult,
               MyOccsMultExtra const& occsMultExtra, MyOccsRobust const& occsRobust)
  {
    // dfCount++;LOG(info) << "DEBUG 2 :: df_" << dfCount ;//<< " :: DF_" << Origins.begin().dataframeID() << " :: collisions.size() = " << collisions.size() << " :: tracks.size() = " << tracks.size() << " :: tracksQA.size() = " << tracksQA.size()
    //           << " :: MyBCTable.size() = " << BCs.size()
    //           << " :: occsBCsList.size() = "  <<occsBCsList.size()
    //           << " :: occsDet.size() = "      <<occsDet.size()
    //           << " :: occsTrackMult.size() = "<<occsTrackMult.size()
    //           <<"  :: occsMultExtra.size() = "<<occsMultExtra.size()
    //           << " :: occsRobust.size() = "   <<occsRobust.size()
    //           << " :: ambgTracks.size() = " <<ambgTracks.size()
    //           ;
    if (collisions.size() == 0) {
      return;
    }
    //
    // BCs.bindExternalIndices(&occsDet);
    // BCs.bindExternalIndices(&occsTrackMult);
    // BCs.bindExternalIndices(&occsRobust);

    auto bc = BCs.begin();

    int64_t oldTFid = -1;

    int64_t oldCollisionIndex = -100;
    bool hasCollision = false;
    bool isAmbgTrack = false;
    bool lastTrackHadCollision = false;
    bool doCollisionUpdate = false;
    bool doAmbgUpdate = false;

    double rBegin = 90., rEnd = 245.;
    double zBegin; // collision.posZ() + track.tgl()*rBegin;
    double zEnd;   // collision.posZ() + track.tgl()*rEnd;
    double vdrift = 2.64;

    double dTbegin; // ((250.- TMath::Abs(zBegin))/vdrift)/0.025;//bin
    double dTend;   // ((250.- TMath::Abs(zEnd))/vdrift)/0.025;  //bin

    double bcBegin; // tGlobalBC + dTbegin;
    double bcEnd;   // tGlobalBC + dTend  ;

    int binBCbegin;
    int binBCend;
    for (const auto& trackQA : tracksQA) {
      auto const& track = trackQA.track_as<MyTracks>();
      auto collision = collisions.begin();

      hasCollision = false;
      isAmbgTrack = false;

      if (track.collisionId() >= 0) {                   // track has collision
        collision = track.collision_as<MyCollisions>(); // It will build but crash while running for tracks with track.collisionId()= -1;//ambg tracks/orphan tracks
        if (track.collisionId() != collision.globalIndex()) {
          LOG(error) << "DEBUG :: ERROR :: track collId and collID Mismatch";
        }
        hasCollision = true;
      } else { // track is ambiguous/orphan
        isAmbgTrack = true;
      }

      // Checking out of the range errors
      if (trackQA.trackId() < 0 || tracks.size() <= trackQA.trackId()) {
        LOG(error) << "DEBUG :: ERROR :: trackQA has index out of scope :: trackQA.trackId() = " << trackQA.trackId() << " :: track.collisionId() = " << track.collisionId() << " :: track.signed1Pt() = " << track.signed1Pt();
      }
      if (!hasCollision && !isAmbgTrack) {
        LOG(error) << "DEBUG :: ERROR :: A track with no collsiion and is not Ambiguous";
      }
      if (hasCollision && isAmbgTrack) {
        LOG(error) << "DEBUG :: ERROR :: A track has collision and is also ambiguous";
      }

      if (hasCollision) {
        lastTrackHadCollision = true;
      }
      doCollisionUpdate = false; // default is false;
      doAmbgUpdate = false;
      if (hasCollision) {
        if (lastTrackHadCollision) {
          if (collision.globalIndex() == oldCollisionIndex) { // if collisions are same
            doCollisionUpdate = false;
          } else { // if collisions are different
            doCollisionUpdate = true;
          }
        } else { // LastTrackWasAmbiguous
          doCollisionUpdate = true;
        }
      } else if (isAmbgTrack) {
        doAmbgUpdate = true;
        // To be updated later
        //  if(LastTrackIsAmbg){
        //    if( haveSameInfo ) { doAmbgUpdate = false;}
        //    else              { doAmbgUpdate = true; }
        //  }
        //  else { doAmbgUpdate = true;} //Last track had Collisions
      }

      if (doAmbgUpdate) { // sKipping ambiguous tracks for now, will be updated in future
        continue;
      }
      if (doCollisionUpdate || doAmbgUpdate) { // collision.globalIndex() != oldCollisionIndex){ //don't update if info is same as old collision
        if (doCollisionUpdate) {
          oldCollisionIndex = collision.globalIndex();
          bc = collision.bc_as<MyBCTable>();
        }
        if (doAmbgUpdate) {
          // to be updated later
          //  bc = collisions.iteratorAt(2).bc_as<aod::BCsWithTimestamps>();
          //  bc = ambgTracks.iteratorAt(0).bc_as<aod::BCsWithTimestamps>();
        }
        // LOG(info)<<" What happens in the case when the collision id is = -1 and it tries to obtain bc"
        getTimingInfo(bc, lastRun, nBCsPerTF, bcSOR, time, tfIdThis, bcInTF);
      }

      if (tfIdThis != oldTFid) {
        oldTFid = tfIdThis;
        // auto OccList = Occs.iteratorAt(bc.occId());
        if (buildTrackMeanOccs1 || buildTrackMeanOccs5) {
          auto listOccsDet = occsDet.iteratorAt(bc.occId());
          std::copy(listOccsDet.occPrimUnfm80().begin(), listOccsDet.occPrimUnfm80().end(), occPrimUnfm80.begin());
          std::copy(listOccsDet.occFV0AUnfm80().begin(), listOccsDet.occFV0AUnfm80().end(), occFV0AUnfm80.begin());
          std::copy(listOccsDet.occFV0CUnfm80().begin(), listOccsDet.occFV0CUnfm80().end(), occFV0CUnfm80.begin());
          std::copy(listOccsDet.occFT0AUnfm80().begin(), listOccsDet.occFT0AUnfm80().end(), occFT0AUnfm80.begin());
          std::copy(listOccsDet.occFT0CUnfm80().begin(), listOccsDet.occFT0CUnfm80().end(), occFT0CUnfm80.begin());
          std::copy(listOccsDet.occFDDAUnfm80().begin(), listOccsDet.occFDDAUnfm80().end(), occFDDAUnfm80.begin());
          std::copy(listOccsDet.occFDDCUnfm80().begin(), listOccsDet.occFDDCUnfm80().end(), occFDDCUnfm80.begin());
        }

        if (buildTrackMeanOccs2 || buildTrackMeanOccs6) {
          auto listOccsTrackMult = occsTrackMult.iteratorAt(bc.occId());
          ;
          std::copy(listOccsTrackMult.occNTrackITSUnfm80().begin(), listOccsTrackMult.occNTrackITSUnfm80().end(), occNTrackITSUnfm80.begin());
          std::copy(listOccsTrackMult.occNTrackTPCUnfm80().begin(), listOccsTrackMult.occNTrackTPCUnfm80().end(), occNTrackTPCUnfm80.begin());
          std::copy(listOccsTrackMult.occNTrackTRDUnfm80().begin(), listOccsTrackMult.occNTrackTRDUnfm80().end(), occNTrackTRDUnfm80.begin());
          std::copy(listOccsTrackMult.occNTrackTOFUnfm80().begin(), listOccsTrackMult.occNTrackTOFUnfm80().end(), occNTrackTOFUnfm80.begin());
          std::copy(listOccsTrackMult.occNTrackSizeUnfm80().begin(), listOccsTrackMult.occNTrackSizeUnfm80().end(), occNTrackSizeUnfm80.begin());
          std::copy(listOccsTrackMult.occNTrackTPCAUnfm80().begin(), listOccsTrackMult.occNTrackTPCAUnfm80().end(), occNTrackTPCAUnfm80.begin());
          std::copy(listOccsTrackMult.occNTrackTPCCUnfm80().begin(), listOccsTrackMult.occNTrackTPCCUnfm80().end(), occNTrackTPCCUnfm80.begin());
          std::copy(listOccsTrackMult.occNTrackITSTPCUnfm80().begin(), listOccsTrackMult.occNTrackITSTPCUnfm80().end(), occNTrackITSTPCUnfm80.begin());
          std::copy(listOccsTrackMult.occNTrackITSTPCAUnfm80().begin(), listOccsTrackMult.occNTrackITSTPCAUnfm80().end(), occNTrackITSTPCAUnfm80.begin());
          std::copy(listOccsTrackMult.occNTrackITSTPCCUnfm80().begin(), listOccsTrackMult.occNTrackITSTPCCUnfm80().end(), occNTrackITSTPCCUnfm80.begin());
        }

        if (buildTrackMeanOccs3 || buildTrackMeanOccs7) {
          auto listOccsMultExtra = occsMultExtra.iteratorAt(bc.occId());
          ;
          std::copy(listOccsMultExtra.occMultNTracksHasITSUnfm80().begin(), listOccsMultExtra.occMultNTracksHasITSUnfm80().end(), occMultNTracksHasITSUnfm80.begin());
          std::copy(listOccsMultExtra.occMultNTracksHasTPCUnfm80().begin(), listOccsMultExtra.occMultNTracksHasTPCUnfm80().end(), occMultNTracksHasTPCUnfm80.begin());
          std::copy(listOccsMultExtra.occMultNTracksHasTOFUnfm80().begin(), listOccsMultExtra.occMultNTracksHasTOFUnfm80().end(), occMultNTracksHasTOFUnfm80.begin());
          std::copy(listOccsMultExtra.occMultNTracksHasTRDUnfm80().begin(), listOccsMultExtra.occMultNTracksHasTRDUnfm80().end(), occMultNTracksHasTRDUnfm80.begin());
          std::copy(listOccsMultExtra.occMultNTracksITSOnlyUnfm80().begin(), listOccsMultExtra.occMultNTracksITSOnlyUnfm80().end(), occMultNTracksITSOnlyUnfm80.begin());
          std::copy(listOccsMultExtra.occMultNTracksTPCOnlyUnfm80().begin(), listOccsMultExtra.occMultNTracksTPCOnlyUnfm80().end(), occMultNTracksTPCOnlyUnfm80.begin());
          std::copy(listOccsMultExtra.occMultNTracksITSTPCUnfm80().begin(), listOccsMultExtra.occMultNTracksITSTPCUnfm80().end(), occMultNTracksITSTPCUnfm80.begin());
          std::copy(listOccsMultExtra.occMultAllTracksTPCOnlyUnfm80().begin(), listOccsMultExtra.occMultAllTracksTPCOnlyUnfm80().end(), occMultAllTracksTPCOnlyUnfm80.begin());
        }

        if (buildTrackMeanOccs4 || buildTrackMeanOccs8) {
          auto listOccsRobust = occsRobust.iteratorAt(bc.occId());
          std::copy(listOccsRobust.occRobustT0V0PrimUnfm80().begin(), listOccsRobust.occRobustT0V0PrimUnfm80().end(), occRobustT0V0PrimUnfm80.begin());
          std::copy(listOccsRobust.occRobustFDDT0V0PrimUnfm80().begin(), listOccsRobust.occRobustFDDT0V0PrimUnfm80().end(), occRobustFDDT0V0PrimUnfm80.begin());
          std::copy(listOccsRobust.occRobustNtrackDetUnfm80().begin(), listOccsRobust.occRobustNtrackDetUnfm80().end(), occRobustNtrackDetUnfm80.begin());
          std::copy(listOccsRobust.occRobustMultExtraTableUnfm80().begin(), listOccsRobust.occRobustMultExtraTableUnfm80().end(), occRobustMultTableUnfm80.begin());
        }
      }

      // Timebc = TGlobalBC+ΔTdrift
      // ΔTdrift=((250(cm)-abs(z))/vdrift)
      // vdrift=2.64 cm/μs
      // z=zv+tgl*Radius

      rBegin = 90., rEnd = 245.;                        // in cm
      zBegin = collision.posZ() + track.tgl() * rBegin; // in cm
      zEnd = collision.posZ() + track.tgl() * rEnd;     // in cm
      vdrift = 2.64;                                    // cm/μs

      // clip the result at 250
      if (zBegin > 250) {
        zBegin = 250;
      } else if (zBegin < -250) {
        zBegin = -250;
      }

      if (zEnd > 250) {
        zEnd = 250;
      } else if (zEnd < -250) {
        zEnd = -250;
      }

      dTbegin = ((250. - std::abs(zBegin)) / vdrift) / 0.025;
      dTend = ((250. - std::abs(zEnd)) / vdrift) / 0.025;

      bcBegin = bcInTF + dTbegin;
      bcEnd = bcInTF + dTend;

      binBCbegin = bcBegin / 80;
      binBCend = bcEnd / 80;

      genTrackMeanOccs0(track.globalIndex());

      if (buildTrackMeanOccs1) {
        genTrackMeanOccs1(getMeanOccupancy(binBCbegin, binBCend, occPrimUnfm80),
                          getMeanOccupancy(binBCbegin, binBCend, occFV0AUnfm80),
                          getMeanOccupancy(binBCbegin, binBCend, occFV0CUnfm80),
                          getMeanOccupancy(binBCbegin, binBCend, occFT0AUnfm80),
                          getMeanOccupancy(binBCbegin, binBCend, occFT0CUnfm80),
                          getMeanOccupancy(binBCbegin, binBCend, occFDDAUnfm80),
                          getMeanOccupancy(binBCbegin, binBCend, occFDDCUnfm80));
      }
      if (buildTrackMeanOccs2) {
        genTrackMeanOccs2(
          getMeanOccupancy(binBCbegin, binBCend, occNTrackITSUnfm80),
          getMeanOccupancy(binBCbegin, binBCend, occNTrackTPCUnfm80),
          getMeanOccupancy(binBCbegin, binBCend, occNTrackTRDUnfm80),
          getMeanOccupancy(binBCbegin, binBCend, occNTrackTOFUnfm80),
          getMeanOccupancy(binBCbegin, binBCend, occNTrackSizeUnfm80),
          getMeanOccupancy(binBCbegin, binBCend, occNTrackTPCAUnfm80),
          getMeanOccupancy(binBCbegin, binBCend, occNTrackTPCCUnfm80),
          getMeanOccupancy(binBCbegin, binBCend, occNTrackITSTPCUnfm80),
          getMeanOccupancy(binBCbegin, binBCend, occNTrackITSTPCAUnfm80),
          getMeanOccupancy(binBCbegin, binBCend, occNTrackITSTPCCUnfm80));
      }
      if (buildTrackMeanOccs3) {
        genTrackMeanOccs3(getMeanOccupancy(binBCbegin, binBCend, occMultNTracksHasITSUnfm80),
                          getMeanOccupancy(binBCbegin, binBCend, occMultNTracksHasTPCUnfm80),
                          getMeanOccupancy(binBCbegin, binBCend, occMultNTracksHasTOFUnfm80),
                          getMeanOccupancy(binBCbegin, binBCend, occMultNTracksHasTRDUnfm80),
                          getMeanOccupancy(binBCbegin, binBCend, occMultNTracksITSOnlyUnfm80),
                          getMeanOccupancy(binBCbegin, binBCend, occMultNTracksTPCOnlyUnfm80),
                          getMeanOccupancy(binBCbegin, binBCend, occMultNTracksITSTPCUnfm80),
                          getMeanOccupancy(binBCbegin, binBCend, occMultAllTracksTPCOnlyUnfm80));
      }
      if (buildTrackMeanOccs4) {
        genTrackMeanOccs4(getMeanOccupancy(binBCbegin, binBCend, occRobustT0V0PrimUnfm80),
                          getMeanOccupancy(binBCbegin, binBCend, occRobustFDDT0V0PrimUnfm80),
                          getMeanOccupancy(binBCbegin, binBCend, occRobustNtrackDetUnfm80),
                          getMeanOccupancy(binBCbegin, binBCend, occRobustMultTableUnfm80));
      }

      if (buildTrackMeanOccs5) {
        genTrackMeanOccs5(getWeightedMeanOccupancy(binBCbegin, binBCend, occPrimUnfm80),
                          getWeightedMeanOccupancy(binBCbegin, binBCend, occFV0AUnfm80),
                          getWeightedMeanOccupancy(binBCbegin, binBCend, occFV0CUnfm80),
                          getWeightedMeanOccupancy(binBCbegin, binBCend, occFT0AUnfm80),
                          getWeightedMeanOccupancy(binBCbegin, binBCend, occFT0CUnfm80),
                          getWeightedMeanOccupancy(binBCbegin, binBCend, occFDDAUnfm80),
                          getWeightedMeanOccupancy(binBCbegin, binBCend, occFDDCUnfm80));
      }

      if (buildTrackMeanOccs6) {
        genTrackMeanOccs6(
          getWeightedMeanOccupancy(binBCbegin, binBCend, occNTrackITSUnfm80),
          getWeightedMeanOccupancy(binBCbegin, binBCend, occNTrackTPCUnfm80),
          getWeightedMeanOccupancy(binBCbegin, binBCend, occNTrackTRDUnfm80),
          getWeightedMeanOccupancy(binBCbegin, binBCend, occNTrackTOFUnfm80),
          getWeightedMeanOccupancy(binBCbegin, binBCend, occNTrackSizeUnfm80),
          getWeightedMeanOccupancy(binBCbegin, binBCend, occNTrackTPCAUnfm80),
          getWeightedMeanOccupancy(binBCbegin, binBCend, occNTrackTPCCUnfm80),
          getWeightedMeanOccupancy(binBCbegin, binBCend, occNTrackITSTPCUnfm80),
          getWeightedMeanOccupancy(binBCbegin, binBCend, occNTrackITSTPCAUnfm80),
          getWeightedMeanOccupancy(binBCbegin, binBCend, occNTrackITSTPCCUnfm80));
      }

      if (buildTrackMeanOccs7) {
        genTrackMeanOccs7(
          getWeightedMeanOccupancy(binBCbegin, binBCend, occMultNTracksHasITSUnfm80),
          getWeightedMeanOccupancy(binBCbegin, binBCend, occMultNTracksHasTPCUnfm80),
          getWeightedMeanOccupancy(binBCbegin, binBCend, occMultNTracksHasTOFUnfm80),
          getWeightedMeanOccupancy(binBCbegin, binBCend, occMultNTracksHasTRDUnfm80),
          getWeightedMeanOccupancy(binBCbegin, binBCend, occMultNTracksITSOnlyUnfm80),
          getWeightedMeanOccupancy(binBCbegin, binBCend, occMultNTracksTPCOnlyUnfm80),
          getWeightedMeanOccupancy(binBCbegin, binBCend, occMultNTracksITSTPCUnfm80),
          getWeightedMeanOccupancy(binBCbegin, binBCend, occMultAllTracksTPCOnlyUnfm80));
      }

      if (buildTrackMeanOccs8) {
        genTrackMeanOccs8(getWeightedMeanOccupancy(binBCbegin, binBCend, occRobustT0V0PrimUnfm80),
                          getWeightedMeanOccupancy(binBCbegin, binBCend, occRobustFDDT0V0PrimUnfm80),
                          getWeightedMeanOccupancy(binBCbegin, binBCend, occRobustNtrackDetUnfm80),
                          getWeightedMeanOccupancy(binBCbegin, binBCend, occRobustMultTableUnfm80));
      }
    } // end of trackQA loop
  } // Process function ends
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<OccupancyTableProducer>(cfgc),
    adaptAnalysisTask<TrackMeanOccTableProducer>(cfgc)};
}
