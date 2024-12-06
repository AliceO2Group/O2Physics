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

#include "OccupancyTables.h"

#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsFT0/Digit.h"
#include "DataFormatsParameters/GRPLHCIFData.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

int32_t nBCsPerOrbit = o2::constants::lhc::LHCMaxBunches;
const int nBCinTF_grp80 = 1425;
// for 128 => nBCsPerTF = 456192 , for 32 => nBCsPerTF = 114048 // nOrbitsPerTF = run < 534133 ? 128 : 32;
const int nBCinTF = 114048;         /// CCDB value // to be obtained from CCDB in future
const int nBCinDrift = 114048 / 32; /// to get from ccdb in future
const int arraySize = 3;            // Max no timeframes that can be present in a dataframe

struct occTableProducer {

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  // declare production of tables
  Produces<aod::OccIndexTable> GenOccIndexTable;
  Produces<aod::BCTFinfoTable> GenBCTFinfoTable;
  Produces<aod::OccsBCsList> GenOccsBCsList;

  Produces<aod::OccsDet> GenOccsDet;
  Produces<aod::OccsTrackMult> GenOccsTrackMult;
  Produces<aod::OccsMultExtra> GenOccsMultExtra;
  Produces<aod::OccsRobust> GenOccsRobust;

  Produces<aod::OccsMeanDet> GenOccsMeanDet;
  Produces<aod::OccsMeanTrkMult> GenOccsMeanTrkMult;
  Produces<aod::OccsMnMultExtra> GenOccsMnMultExtra;
  Produces<aod::OccsMeanRobust> GenOccsMeanRobust;

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
  std::array<int, arraySize> TFlist;
  std::array<std::vector<int64_t>, arraySize> BC_TF_Map;
  std::array<std::vector<float>, arraySize> occ_Prim_Unfm_80;
  std::array<std::vector<float>, arraySize> occ_FV0A_Unfm_80;
  std::array<std::vector<float>, arraySize> occ_FV0C_Unfm_80;
  std::array<std::vector<float>, arraySize> occ_FT0A_Unfm_80;
  std::array<std::vector<float>, arraySize> occ_FT0C_Unfm_80;
  std::array<std::vector<float>, arraySize> occ_FDDA_Unfm_80;
  std::array<std::vector<float>, arraySize> occ_FDDC_Unfm_80;

  std::array<std::vector<float>, arraySize> occ_NTrack_ITS_Unfm_80;
  std::array<std::vector<float>, arraySize> occ_NTrack_TPC_Unfm_80;
  std::array<std::vector<float>, arraySize> occ_NTrack_TRD_Unfm_80;
  std::array<std::vector<float>, arraySize> occ_NTrack_TOF_Unfm_80;
  std::array<std::vector<float>, arraySize> occ_NTrackSize_Unfm_80;
  std::array<std::vector<float>, arraySize> occ_NTrackTPC_A_Unfm_80;
  std::array<std::vector<float>, arraySize> occ_NTrackTPC_C_Unfm_80;
  std::array<std::vector<float>, arraySize> occ_NTrackITS_TPC_Unfm_80;
  std::array<std::vector<float>, arraySize> occ_NTrackITS_TPC_A_Unfm_80;
  std::array<std::vector<float>, arraySize> occ_NTrackITS_TPC_C_Unfm_80;

  std::array<std::vector<float>, arraySize> occ_multNTracksHasITS_Unfm_80;
  std::array<std::vector<float>, arraySize> occ_multNTracksHasTPC_Unfm_80;
  std::array<std::vector<float>, arraySize> occ_multNTracksHasTOF_Unfm_80;
  std::array<std::vector<float>, arraySize> occ_multNTracksHasTRD_Unfm_80;
  std::array<std::vector<float>, arraySize> occ_multNTracksITSOnly_Unfm_80;
  std::array<std::vector<float>, arraySize> occ_multNTracksTPCOnly_Unfm_80;
  std::array<std::vector<float>, arraySize> occ_multNTracksITSTPC_Unfm_80;
  std::array<std::vector<float>, arraySize> occ_multAllTracksTPCOnly_Unfm_80;

  std::vector<float> vecRobustOcc_T0V0Prim_Unfm_80;
  std::vector<float> vecRobustOcc_FDDT0V0Prim_Unfm_80;
  std::vector<float> vecRobustOcc_NtrackDet_Unfm_80;
  std::vector<float> vecRobustOcc_multTable_Unfm_80;
  std::vector<std::array<int, 2>> vecRobustOcc_T0V0Prim_Unfm_80_medianPosVec; // Median => one for odd and two for even entries
  std::vector<std::array<int, 2>> vecRobustOcc_FDDT0V0Prim_Unfm_80_medianPosVec;
  std::vector<std::array<int, 2>> vecRobustOcc_NtrackDet_Unfm_80_medianPosVec;
  std::vector<std::array<int, 2>> vecRobustOcc_multTable_Unfm_80_medianPosVec;

  void init(InitContext const&)
  {
    // Set size of the vectors
    for (int i = 0; i < arraySize; i++) {
      BC_TF_Map[i].resize(nBCinTF / 80);
      occ_Prim_Unfm_80[i].resize(nBCinTF / 80);
      occ_FV0A_Unfm_80[i].resize(nBCinTF / 80);
      occ_FV0C_Unfm_80[i].resize(nBCinTF / 80);
      occ_FT0A_Unfm_80[i].resize(nBCinTF / 80);
      occ_FT0C_Unfm_80[i].resize(nBCinTF / 80);
      occ_FDDA_Unfm_80[i].resize(nBCinTF / 80);
      occ_FDDC_Unfm_80[i].resize(nBCinTF / 80);

      occ_NTrack_ITS_Unfm_80[i].resize(nBCinTF / 80);
      occ_NTrack_TPC_Unfm_80[i].resize(nBCinTF / 80);
      occ_NTrack_TRD_Unfm_80[i].resize(nBCinTF / 80);
      occ_NTrack_TOF_Unfm_80[i].resize(nBCinTF / 80);
      occ_NTrackSize_Unfm_80[i].resize(nBCinTF / 80);
      occ_NTrackTPC_A_Unfm_80[i].resize(nBCinTF / 80);
      occ_NTrackTPC_C_Unfm_80[i].resize(nBCinTF / 80);
      occ_NTrackITS_TPC_Unfm_80[i].resize(nBCinTF / 80);
      occ_NTrackITS_TPC_A_Unfm_80[i].resize(nBCinTF / 80);
      occ_NTrackITS_TPC_C_Unfm_80[i].resize(nBCinTF / 80);

      occ_multNTracksHasITS_Unfm_80[i].resize(nBCinTF / 80);
      occ_multNTracksHasTPC_Unfm_80[i].resize(nBCinTF / 80);
      occ_multNTracksHasTOF_Unfm_80[i].resize(nBCinTF / 80);
      occ_multNTracksHasTRD_Unfm_80[i].resize(nBCinTF / 80);
      occ_multNTracksITSOnly_Unfm_80[i].resize(nBCinTF / 80);
      occ_multNTracksTPCOnly_Unfm_80[i].resize(nBCinTF / 80);
      occ_multNTracksITSTPC_Unfm_80[i].resize(nBCinTF / 80);
      occ_multAllTracksTPCOnly_Unfm_80[i].resize(nBCinTF / 80);
    }

    vecRobustOcc_T0V0Prim_Unfm_80.resize(nBCinTF / 80);
    vecRobustOcc_FDDT0V0Prim_Unfm_80.resize(nBCinTF / 80);
    vecRobustOcc_NtrackDet_Unfm_80.resize(nBCinTF / 80);
    vecRobustOcc_multTable_Unfm_80.resize(nBCinTF / 80);

    vecRobustOcc_T0V0Prim_Unfm_80_medianPosVec.resize(nBCinTF / 80); // Median => one for odd and two for even entries
    vecRobustOcc_FDDT0V0Prim_Unfm_80_medianPosVec.resize(nBCinTF / 80);
    vecRobustOcc_NtrackDet_Unfm_80_medianPosVec.resize(nBCinTF / 80);
    vecRobustOcc_multTable_Unfm_80_medianPosVec.resize(nBCinTF / 80);

    // Getting Info from CCDB, to be implemented Later
    recoEvent.add("h_nBCinTF", "h_nBCinTF(to check nBCinTF)", {HistType::kTH1F, {{100, 114040, 114060}}}); // 114048
    recoEvent.add("h_bcInTF", "h_bcInTF", {HistType::kTH1F, {{2000, 0, 200000}}});
    recoEvent.add("h_RO_T0V0Prim_Unfm_80", "h_RO_T0V0Prim_Unfm_80:median contributors", {HistType::kTH1F, {{12 * 2, -1, 11}}});
    recoEvent.add("h_RO_FDDT0V0Prim_Unfm_80", "h_RO_FDDT0V0Prim_Unfm_80:median contributors", {HistType::kTH1F, {{12 * 2, -1, 11}}});
    recoEvent.add("h_RO_NtrackDet_Unfm_80", "h_RO_NtrackDetITS/TPC/TRD/TOF_80:median contributors", {HistType::kTH1F, {{12 * 2, -1, 11}}});
    recoEvent.add("h_RO_multTable_Unfm_80", "h_RO_multTableExtra_80:median contributors", {HistType::kTH1F, {{12 * 2, -1, 11}}});

    recoEvent.print();
  }

  void NormalizeVector(std::vector<float>& OriginalVec, const float& scaleFactor)
  {
    std::transform(OriginalVec.begin(), OriginalVec.end(), OriginalVec.begin(), [scaleFactor](float x) { return x * scaleFactor; });
  }

  template <typename... Vecs>
  void GetMedianOccVect(
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

  void GetRunInfo(const int& run, int& nBCsPerTF, int64_t& bcSOR)
  {
    auto runDuration = ccdb->getRunDuration(run, true);
    int64_t tsSOR = runDuration.first;
    auto ctpx = ccdb->getForTimeStamp<std::vector<Long64_t>>("CTP/Calib/OrbitReset", tsSOR);
    int64_t tsOrbitReset = (*ctpx)[0];
    uint32_t nOrbitsPerTF = run < 534133 ? 128 : 32;
    int64_t orbitSOR = (tsSOR * 1000 - tsOrbitReset) / o2::constants::lhc::LHCOrbitMUS;
    orbitSOR = orbitSOR / nOrbitsPerTF * nOrbitsPerTF;
    bcSOR = orbitSOR * nBCsPerOrbit + customOrbitOffset * nBCsPerOrbit; // customOrbitOffset is a configurable
    nBCsPerTF = nOrbitsPerTF * nBCsPerOrbit;
  }

  template <typename T>
  void GetTimingInfo(const T& bc, int& lastRun, int32_t& nBCsPerTF, int64_t& bcSOR, uint64_t& time, int64_t& TFidThis, int& bcInTF)
  {
    int run = bc.runNumber();
    if (run != lastRun) { // update run info
      lastRun = run;
      GetRunInfo(run, nBCsPerTF, bcSOR); // update nBCsPerTF && bcSOR
    }
    // update the information
    time = bc.timestamp();
    TFidThis = (bc.globalBC() - bcSOR) / nBCsPerTF;
    bcInTF = (bc.globalBC() - bcSOR) % nBCsPerTF;
  }

  // using myCollisions = soa::Join<aod::Collisions, aod::Mults>;//, aod::MultsExtra>;
  using myCollisions = soa::Join<aod::Collisions, aod::Mults, aod::MultsExtra>;
  using myTracks = soa::Join<aod::Tracks, o2::aod::TracksCov, aod::TracksExtra>;

  Preslice<myTracks> TracksPerCollisionPreslice = o2::aod::track::collisionId;

  // int      dfCount = 0;
  int32_t nBCsPerTF = -1;
  int64_t bcSOR = -1;
  uint64_t time = -1;
  int64_t TFidThis = -1;
  int bcInTF = -1;
  uint tfCounted = 0;
  int tfIDX = 0;
  int lastRun = -999;

  // Process the Data
  void process(o2::aod::BCsWithTimestamps const& BCs, myCollisions const& collisions, myTracks const& tracks) // aod::TracksQA const& tracksQA, o2::aod::Origins const& Origins //tables only used during debugging
  {
    // dfCount++;
    // LOG(info) << "DEBUG :: df_" << dfCount << " :: DF_" << Origins.begin().dataframeID() << " :: collisions.size() = " << collisions.size() << " :: tracks.size() = " << tracks.size() << " :: tracksQA.size() = " << tracksQA.size() << " :: BCs.size() = " << BCs.size();

    // Defining fix value for {lastRun, nBCsPerTF, bcSOR} at beginning of each dataframe to avoid wierd behaviour//They will be calculated once per DF
    lastRun = -999;
    nBCsPerTF = -999;
    bcSOR = -999;

    if (collisions.size() == 0) {
      for (const auto& BC : BCs) { // For BCs and OccIndexTable to have same size for joining
        GetTimingInfo(BC, lastRun, nBCsPerTF, bcSOR, time, TFidThis, bcInTF);
        GenOccIndexTable(BC.globalIndex(), -999); // BCId, OccId
        GenBCTFinfoTable(TFidThis, bcInTF);
      }
      return;
    }

    // Initialisze the vectors components to zero
    tfIDX = 0;
    tfCounted = 0;
    for (int i = 0; i < arraySize; i++) {
      TFlist[i] = -1;
      BC_TF_Map[i].clear(); // list of BCs used in one time frame;
      std::fill(occ_Prim_Unfm_80[i].begin(), occ_Prim_Unfm_80[i].end(), 0.);
      std::fill(occ_FV0A_Unfm_80[i].begin(), occ_FV0A_Unfm_80[i].end(), 0.);
      std::fill(occ_FV0C_Unfm_80[i].begin(), occ_FV0C_Unfm_80[i].end(), 0.);
      std::fill(occ_FT0A_Unfm_80[i].begin(), occ_FT0A_Unfm_80[i].end(), 0.);
      std::fill(occ_FT0C_Unfm_80[i].begin(), occ_FT0C_Unfm_80[i].end(), 0.);
      std::fill(occ_FDDA_Unfm_80[i].begin(), occ_FDDA_Unfm_80[i].end(), 0.);
      std::fill(occ_FDDC_Unfm_80[i].begin(), occ_FDDC_Unfm_80[i].end(), 0.);

      std::fill(occ_NTrack_ITS_Unfm_80[i].begin(), occ_NTrack_ITS_Unfm_80[i].end(), 0.);
      std::fill(occ_NTrack_TPC_Unfm_80[i].begin(), occ_NTrack_TPC_Unfm_80[i].end(), 0.);
      std::fill(occ_NTrack_TRD_Unfm_80[i].begin(), occ_NTrack_TRD_Unfm_80[i].end(), 0.);
      std::fill(occ_NTrack_TOF_Unfm_80[i].begin(), occ_NTrack_TOF_Unfm_80[i].end(), 0.);
      std::fill(occ_NTrackSize_Unfm_80[i].begin(), occ_NTrackSize_Unfm_80[i].end(), 0.);
      std::fill(occ_NTrackTPC_A_Unfm_80[i].begin(), occ_NTrackTPC_A_Unfm_80[i].end(), 0.);
      std::fill(occ_NTrackTPC_C_Unfm_80[i].begin(), occ_NTrackTPC_C_Unfm_80[i].end(), 0.);
      std::fill(occ_NTrackITS_TPC_Unfm_80[i].begin(), occ_NTrackITS_TPC_Unfm_80[i].end(), 0.);
      std::fill(occ_NTrackITS_TPC_A_Unfm_80[i].begin(), occ_NTrackITS_TPC_A_Unfm_80[i].end(), 0.);
      std::fill(occ_NTrackITS_TPC_C_Unfm_80[i].begin(), occ_NTrackITS_TPC_C_Unfm_80[i].end(), 0.);

      std::fill(occ_multNTracksHasITS_Unfm_80[i].begin(), occ_multNTracksHasITS_Unfm_80[i].end(), 0.);
      std::fill(occ_multNTracksHasTPC_Unfm_80[i].begin(), occ_multNTracksHasTPC_Unfm_80[i].end(), 0.);
      std::fill(occ_multNTracksHasTOF_Unfm_80[i].begin(), occ_multNTracksHasTOF_Unfm_80[i].end(), 0.);
      std::fill(occ_multNTracksHasTRD_Unfm_80[i].begin(), occ_multNTracksHasTRD_Unfm_80[i].end(), 0.);
      std::fill(occ_multNTracksITSOnly_Unfm_80[i].begin(), occ_multNTracksITSOnly_Unfm_80[i].end(), 0.);
      std::fill(occ_multNTracksTPCOnly_Unfm_80[i].begin(), occ_multNTracksTPCOnly_Unfm_80[i].end(), 0.);
      std::fill(occ_multNTracksITSTPC_Unfm_80[i].begin(), occ_multNTracksITSTPC_Unfm_80[i].end(), 0.);
      std::fill(occ_multAllTracksTPCOnly_Unfm_80[i].begin(), occ_multAllTracksTPCOnly_Unfm_80[i].end(), 0.);
    }

    int iColl = -1;
    std::vector<int64_t> TFIDList;
    for (const auto& collision : collisions) {
      iColl++;
      const auto& bc = collision.bc_as<aod::BCsWithTimestamps>();
      GetTimingInfo(bc, lastRun, nBCsPerTF, bcSOR, time, TFidThis, bcInTF);

      recoEvent.fill(HIST("h_nBCinTF"), nBCsPerTF);
      recoEvent.fill(HIST("h_bcInTF"), bcInTF);

      if (nBCsPerTF > nBCinTF) {
        LOG(error) << "DEBUG :: FATAL ERROR :: nBCsPerTF > nBCinTF i.e " << nBCsPerTF << " > " << nBCinTF << " will cause crash in further process";
        return;
      }

      const uint64_t collIdx = collision.globalIndex();
      const auto TracksTable_perColl = tracks.sliceBy(TracksPerCollisionPreslice, collIdx);

      int nTrack_ITS = 0;
      int nTrack_TPC = 0;
      int nTrack_TRD = 0;
      int nTrack_TOF = 0;
      int nTrackTPC_A = 0;
      int nTrackTPC_C = 0;
      int nTrackITS_TPC = 0;
      int nTrackITS_TPC_A = 0;
      int nTrackITS_TPC_C = 0;

      for (const auto& track : TracksTable_perColl) {
        if (track.hasITS()) {
          nTrack_ITS++;
        } // Flag to check if track has a ITS match
        if (track.hasTPC()) {
          nTrack_TPC++;
          if (track.eta() <= 0.0) {
            nTrackTPC_A++; // includes tracks at eta zero as well.
          } else {
            nTrackTPC_C++;
          }
        } // Flag to check if track has a TPC match
        if (track.hasTRD()) {
          nTrack_TRD++;
        } // Flag to check if track has a TRD match
        if (track.hasTOF()) {
          nTrack_TOF++;
        } // Flag to check if track has a TOF measurement
        if (track.hasITS() && track.hasTPC()) {
          nTrackITS_TPC++;
          if (track.eta() <= 0.0) {
            nTrackITS_TPC_A++; // includes tracks at eta zero as well.
          } else {
            nTrackITS_TPC_C++;
          }
        }
      } // track loop

      // if (collision.multNTracksTPCOnly() != 0) {
      //   LOG(error) << "DEBUG :: ERROR = multNTracksTPCOnly != 0" << collision.multNTracksTPCOnly();
      //   return;
      // }
      // if (collision.multAllTracksITSTPC() != nTrackITS_TPC) {
      //   LOG(error) << "DEBUG :: ERROR :: 10 multAllTracksITSTPC :: " << collision.multAllTracksITSTPC() << " != " << nTrackITS_TPC;
      //   return;
      // }

      TFIDList.push_back(TFidThis);

      if (TFlist[tfIDX] != TFidThis) {
        if (tfCounted != 0) {
          tfIDX++;
        } //
        TFlist[tfIDX] = TFidThis;
        tfCounted++;
      }

      BC_TF_Map[tfIDX].push_back(bc.globalIndex());
      auto& TFocc_Prim_Unfm_80 = occ_Prim_Unfm_80[tfIDX];
      auto& TFocc_FV0A_Unfm_80 = occ_FV0A_Unfm_80[tfIDX];
      auto& TFocc_FV0C_Unfm_80 = occ_FV0C_Unfm_80[tfIDX];
      auto& TFocc_FT0A_Unfm_80 = occ_FT0A_Unfm_80[tfIDX];
      auto& TFocc_FT0C_Unfm_80 = occ_FT0C_Unfm_80[tfIDX];
      auto& TFocc_FDDA_Unfm_80 = occ_FDDA_Unfm_80[tfIDX];
      auto& TFocc_FDDC_Unfm_80 = occ_FDDC_Unfm_80[tfIDX];

      auto& TFocc_NTrack_ITS_Unfm_80 = occ_NTrack_ITS_Unfm_80[tfIDX];
      auto& TFocc_NTrack_TPC_Unfm_80 = occ_NTrack_TPC_Unfm_80[tfIDX];
      auto& TFocc_NTrack_TRD_Unfm_80 = occ_NTrack_TRD_Unfm_80[tfIDX];
      auto& TFocc_NTrack_TOF_Unfm_80 = occ_NTrack_TOF_Unfm_80[tfIDX];
      auto& TFocc_NTrackSize_Unfm_80 = occ_NTrackSize_Unfm_80[tfIDX];
      auto& TFocc_NTrackTPC_A_Unfm_80 = occ_NTrackTPC_A_Unfm_80[tfIDX];
      auto& TFocc_NTrackTPC_C_Unfm_80 = occ_NTrackTPC_C_Unfm_80[tfIDX];
      auto& TFocc_NTrackITS_TPC_Unfm_80 = occ_NTrackITS_TPC_Unfm_80[tfIDX];
      auto& TFocc_NTrackITS_TPC_A_Unfm_80 = occ_NTrackITS_TPC_A_Unfm_80[tfIDX];
      auto& TFocc_NTrackITS_TPC_C_Unfm_80 = occ_NTrackITS_TPC_C_Unfm_80[tfIDX];

      auto& TFocc_multNTracksHasITS_Unfm_80 = occ_multNTracksHasITS_Unfm_80[tfIDX];
      auto& TFocc_multNTracksHasTPC_Unfm_80 = occ_multNTracksHasTPC_Unfm_80[tfIDX];
      auto& TFocc_multNTracksHasTOF_Unfm_80 = occ_multNTracksHasTOF_Unfm_80[tfIDX];
      auto& TFocc_multNTracksHasTRD_Unfm_80 = occ_multNTracksHasTRD_Unfm_80[tfIDX];
      auto& TFocc_multNTracksITSOnly_Unfm_80 = occ_multNTracksITSOnly_Unfm_80[tfIDX];
      auto& TFocc_multNTracksTPCOnly_Unfm_80 = occ_multNTracksTPCOnly_Unfm_80[tfIDX];
      auto& TFocc_multNTracksITSTPC_Unfm_80 = occ_multNTracksITSTPC_Unfm_80[tfIDX];
      auto& TFocc_multAllTracksTPCOnly_Unfm_80 = occ_multAllTracksTPCOnly_Unfm_80[tfIDX];

      // current collision bin in 80/160 grouping.
      int bin80_0 = bcInTF / 80;
      // int bin160_0=bcInTF/160;

      // float fbin80_0 =float(bcInTF)/80;
      // float fbin160_0=float(bcInTF)/160;

      ushort fNumContrib = collision.numContrib();
      float fMultFV0A = collision.multFV0A(), fMultFV0C = collision.multFV0C();
      float fMultFT0A = collision.multFT0A(), fMultFT0C = collision.multFT0C();
      float fMultFDDA = collision.multFDDA(), fMultFDDC = collision.multFDDC();
      int fNTrack_ITS = nTrack_ITS;
      int fNTrack_TPC = nTrack_TPC, fNTrack_TRD = nTrack_TRD;
      int fNTrack_TOF = nTrack_TOF;
      int fNTrackTPC_A = nTrackTPC_A, fNTrackTPC_C = nTrackTPC_C;
      int fNTrackITS_TPC = collision.multAllTracksITSTPC(), fNTrackSize = TracksTable_perColl.size();
      int fNTrackITS_TPC_A = nTrackITS_TPC_A, fNTrackITS_TPC_C = nTrackITS_TPC_C;

      // Processing for grouping of 80 BCs
      for (int deltaBin = 0; deltaBin < nBCinDrift / 80; deltaBin++) {
        TFocc_Prim_Unfm_80[(bin80_0 + deltaBin) % (nBCinTF / 80)] += fNumContrib * 1;
        TFocc_FV0A_Unfm_80[(bin80_0 + deltaBin) % (nBCinTF / 80)] += fMultFV0A * 1;
        TFocc_FV0C_Unfm_80[(bin80_0 + deltaBin) % (nBCinTF / 80)] += fMultFV0C * 1;
        TFocc_FT0A_Unfm_80[(bin80_0 + deltaBin) % (nBCinTF / 80)] += fMultFT0A * 1;
        TFocc_FT0C_Unfm_80[(bin80_0 + deltaBin) % (nBCinTF / 80)] += fMultFT0C * 1;
        TFocc_FDDA_Unfm_80[(bin80_0 + deltaBin) % (nBCinTF / 80)] += fMultFDDA * 1;
        TFocc_FDDC_Unfm_80[(bin80_0 + deltaBin) % (nBCinTF / 80)] += fMultFDDC * 1;

        TFocc_NTrack_ITS_Unfm_80[(bin80_0 + deltaBin) % (nBCinTF / 80)] += fNTrack_ITS * 1;
        TFocc_NTrack_TPC_Unfm_80[(bin80_0 + deltaBin) % (nBCinTF / 80)] += fNTrack_TPC * 1;
        TFocc_NTrack_TRD_Unfm_80[(bin80_0 + deltaBin) % (nBCinTF / 80)] += fNTrack_TRD * 1;
        TFocc_NTrack_TOF_Unfm_80[(bin80_0 + deltaBin) % (nBCinTF / 80)] += fNTrack_TOF * 1;
        TFocc_NTrackSize_Unfm_80[(bin80_0 + deltaBin) % (nBCinTF / 80)] += fNTrackSize * 1;
        TFocc_NTrackTPC_A_Unfm_80[(bin80_0 + deltaBin) % (nBCinTF / 80)] += fNTrackTPC_A * 1;
        TFocc_NTrackTPC_C_Unfm_80[(bin80_0 + deltaBin) % (nBCinTF / 80)] += fNTrackTPC_C * 1;
        TFocc_NTrackITS_TPC_Unfm_80[(bin80_0 + deltaBin) % (nBCinTF / 80)] += fNTrackITS_TPC * 1;
        TFocc_NTrackITS_TPC_A_Unfm_80[(bin80_0 + deltaBin) % (nBCinTF / 80)] += fNTrackITS_TPC_A * 1;
        TFocc_NTrackITS_TPC_C_Unfm_80[(bin80_0 + deltaBin) % (nBCinTF / 80)] += fNTrackITS_TPC_C * 1;

        TFocc_multNTracksHasITS_Unfm_80[(bin80_0 + deltaBin) % (nBCinTF / 80)] += collision.multNTracksHasITS() * 1;
        TFocc_multNTracksHasTPC_Unfm_80[(bin80_0 + deltaBin) % (nBCinTF / 80)] += collision.multNTracksHasTPC() * 1;
        TFocc_multNTracksHasTOF_Unfm_80[(bin80_0 + deltaBin) % (nBCinTF / 80)] += collision.multNTracksHasTOF() * 1;
        TFocc_multNTracksHasTRD_Unfm_80[(bin80_0 + deltaBin) % (nBCinTF / 80)] += collision.multNTracksHasTRD() * 1;
        TFocc_multNTracksITSOnly_Unfm_80[(bin80_0 + deltaBin) % (nBCinTF / 80)] += collision.multNTracksITSOnly() * 1;
        TFocc_multNTracksTPCOnly_Unfm_80[(bin80_0 + deltaBin) % (nBCinTF / 80)] += collision.multNTracksTPCOnly() * 1;
        TFocc_multNTracksITSTPC_Unfm_80[(bin80_0 + deltaBin) % (nBCinTF / 80)] += collision.multNTracksITSTPC() * 1;
        TFocc_multAllTracksTPCOnly_Unfm_80[(bin80_0 + deltaBin) % (nBCinTF / 80)] += collision.multAllTracksTPCOnly() * 1;
      }
    }
    // collision Loop is over

    std::vector<int64_t> SortedTFIDList = TFIDList;
    std::sort(SortedTFIDList.begin(), SortedTFIDList.end());
    auto last = std::unique(SortedTFIDList.begin(), SortedTFIDList.end());
    SortedTFIDList.erase(last, SortedTFIDList.end());

    if (tfCounted != SortedTFIDList.size()) {
      LOG(error) << "DEBUG :: Number mismatch for tf counted and filled :: " << tfCounted << " != " << SortedTFIDList.size();
    }

    int totalBCcountSize = 0;
    for (int i = 0; i < arraySize; i++) {
      totalBCcountSize += BC_TF_Map[i].size();
      // check if the BCs are already sorted or not
      if (!std::is_sorted(BC_TF_Map[i].begin(), BC_TF_Map[i].end())) {
        LOG(debug) << "DEBUG :: ERROR :: BCs are not sorted";
      }
    }
    //
    if (totalBCcountSize != collisions.size()) {
      LOG(debug) << "DEBUG :: ERROR :: filled TF list and collision size mismatch ::  filledTF_Size = " << totalBCcountSize << " != " << collisions.size() << " = collisions.size()";
    }

    // Fill the Producers
    for (uint i = 0; i < tfCounted; i++) {

      auto& vecOcc_Prim_Unfm_80 = occ_Prim_Unfm_80[i];
      auto& vecOcc_FV0A_Unfm_80 = occ_FV0A_Unfm_80[i];
      auto& vecOcc_FV0C_Unfm_80 = occ_FV0C_Unfm_80[i];
      auto& vecOcc_FT0A_Unfm_80 = occ_FT0A_Unfm_80[i];
      auto& vecOcc_FT0C_Unfm_80 = occ_FT0C_Unfm_80[i];
      auto& vecOcc_FDDA_Unfm_80 = occ_FDDA_Unfm_80[i];
      auto& vecOcc_FDDC_Unfm_80 = occ_FDDC_Unfm_80[i];

      auto& vecOcc_NTrack_ITS_Unfm_80 = occ_NTrack_ITS_Unfm_80[i];
      auto& vecOcc_NTrack_TPC_Unfm_80 = occ_NTrack_TPC_Unfm_80[i];
      auto& vecOcc_NTrack_TRD_Unfm_80 = occ_NTrack_TRD_Unfm_80[i];
      auto& vecOcc_NTrack_TOF_Unfm_80 = occ_NTrack_TOF_Unfm_80[i];
      auto& vecOcc_NTrackSize_Unfm_80 = occ_NTrackSize_Unfm_80[i];
      auto& vecOcc_NTrackTPC_A_Unfm_80 = occ_NTrackTPC_A_Unfm_80[i];
      auto& vecOcc_NTrackTPC_C_Unfm_80 = occ_NTrackTPC_C_Unfm_80[i];
      auto& vecOcc_NTrackITS_TPC_Unfm_80 = occ_NTrackITS_TPC_Unfm_80[i];
      auto& vecOcc_NTrackITS_TPC_A_Unfm_80 = occ_NTrackITS_TPC_A_Unfm_80[i];
      auto& vecOcc_NTrackITS_TPC_C_Unfm_80 = occ_NTrackITS_TPC_C_Unfm_80[i];

      auto& vecOcc_multNTracksHasITS_Unfm_80 = occ_multNTracksHasITS_Unfm_80[i];
      auto& vecOcc_multNTracksHasTPC_Unfm_80 = occ_multNTracksHasTPC_Unfm_80[i];
      auto& vecOcc_multNTracksHasTOF_Unfm_80 = occ_multNTracksHasTOF_Unfm_80[i];
      auto& vecOcc_multNTracksHasTRD_Unfm_80 = occ_multNTracksHasTRD_Unfm_80[i];
      auto& vecOcc_multNTracksITSOnly_Unfm_80 = occ_multNTracksITSOnly_Unfm_80[i];
      auto& vecOcc_multNTracksTPCOnly_Unfm_80 = occ_multNTracksTPCOnly_Unfm_80[i];
      auto& vecOcc_multNTracksITSTPC_Unfm_80 = occ_multNTracksITSTPC_Unfm_80[i];
      auto& vecOcc_multAllTracksTPCOnly_Unfm_80 = occ_multAllTracksTPCOnly_Unfm_80[i];

      float meanOcc_Prim_Unfm_80 = TMath::Mean(vecOcc_Prim_Unfm_80.size(), vecOcc_Prim_Unfm_80.data());
      float meanOcc_FV0A_Unfm_80 = TMath::Mean(vecOcc_FV0A_Unfm_80.size(), vecOcc_FV0A_Unfm_80.data());
      float meanOcc_FV0C_Unfm_80 = TMath::Mean(vecOcc_FV0C_Unfm_80.size(), vecOcc_FV0C_Unfm_80.data());
      float meanOcc_FT0A_Unfm_80 = TMath::Mean(vecOcc_FT0A_Unfm_80.size(), vecOcc_FT0A_Unfm_80.data());
      float meanOcc_FT0C_Unfm_80 = TMath::Mean(vecOcc_FT0C_Unfm_80.size(), vecOcc_FT0C_Unfm_80.data());
      float meanOcc_FDDA_Unfm_80 = TMath::Mean(vecOcc_FDDA_Unfm_80.size(), vecOcc_FDDA_Unfm_80.data());
      float meanOcc_FDDC_Unfm_80 = TMath::Mean(vecOcc_FDDC_Unfm_80.size(), vecOcc_FDDC_Unfm_80.data());

      float meanOcc_NTrack_ITS_Unfm_80 = TMath::Mean(vecOcc_NTrack_ITS_Unfm_80.size(), vecOcc_NTrack_ITS_Unfm_80.data());
      float meanOcc_NTrack_TPC_Unfm_80 = TMath::Mean(vecOcc_NTrack_TPC_Unfm_80.size(), vecOcc_NTrack_TPC_Unfm_80.data());
      float meanOcc_NTrack_TRD_Unfm_80 = TMath::Mean(vecOcc_NTrack_TRD_Unfm_80.size(), vecOcc_NTrack_TRD_Unfm_80.data());
      float meanOcc_NTrack_TOF_Unfm_80 = TMath::Mean(vecOcc_NTrack_TOF_Unfm_80.size(), vecOcc_NTrack_TOF_Unfm_80.data());
      float meanOcc_NTrackSize_Unfm_80 = TMath::Mean(vecOcc_NTrackSize_Unfm_80.size(), vecOcc_NTrackSize_Unfm_80.data());
      float meanOcc_NTrackTPC_A_Unfm_80 = TMath::Mean(vecOcc_NTrackTPC_A_Unfm_80.size(), vecOcc_NTrackTPC_A_Unfm_80.data());
      float meanOcc_NTrackTPC_C_Unfm_80 = TMath::Mean(vecOcc_NTrackTPC_C_Unfm_80.size(), vecOcc_NTrackTPC_C_Unfm_80.data());
      float meanOcc_NTrackITS_TPC_Unfm_80 = TMath::Mean(vecOcc_NTrackITS_TPC_Unfm_80.size(), vecOcc_NTrackITS_TPC_Unfm_80.data());
      float meanOcc_NTrackITS_TPC_A_Unfm_80 = TMath::Mean(vecOcc_NTrackITS_TPC_A_Unfm_80.size(), vecOcc_NTrackITS_TPC_A_Unfm_80.data());
      float meanOcc_NTrackITS_TPC_C_Unfm_80 = TMath::Mean(vecOcc_NTrackITS_TPC_C_Unfm_80.size(), vecOcc_NTrackITS_TPC_C_Unfm_80.data());

      float meanOcc_multNTracksHasITS_Unfm_80 = TMath::Mean(vecOcc_multNTracksHasITS_Unfm_80.size(), vecOcc_multNTracksHasITS_Unfm_80.data());
      float meanOcc_multNTracksHasTPC_Unfm_80 = TMath::Mean(vecOcc_multNTracksHasTPC_Unfm_80.size(), vecOcc_multNTracksHasTPC_Unfm_80.data());
      float meanOcc_multNTracksHasTOF_Unfm_80 = TMath::Mean(vecOcc_multNTracksHasTOF_Unfm_80.size(), vecOcc_multNTracksHasTOF_Unfm_80.data());
      float meanOcc_multNTracksHasTRD_Unfm_80 = TMath::Mean(vecOcc_multNTracksHasTRD_Unfm_80.size(), vecOcc_multNTracksHasTRD_Unfm_80.data());
      float meanOcc_multNTracksITSOnly_Unfm_80 = TMath::Mean(vecOcc_multNTracksITSOnly_Unfm_80.size(), vecOcc_multNTracksITSOnly_Unfm_80.data());
      float meanOcc_multNTracksTPCOnly_Unfm_80 = TMath::Mean(vecOcc_multNTracksTPCOnly_Unfm_80.size(), vecOcc_multNTracksTPCOnly_Unfm_80.data());
      float meanOcc_multNTracksITSTPC_Unfm_80 = TMath::Mean(vecOcc_multNTracksITSTPC_Unfm_80.size(), vecOcc_multNTracksITSTPC_Unfm_80.data());
      float meanOcc_multAllTracksTPCOnly_Unfm_80 = TMath::Mean(vecOcc_multAllTracksTPCOnly_Unfm_80.size(), vecOcc_multAllTracksTPCOnly_Unfm_80.data());

      // Normalise the original vectors
      NormalizeVector(vecOcc_Prim_Unfm_80, meanOcc_Prim_Unfm_80 / meanOcc_Prim_Unfm_80);
      NormalizeVector(vecOcc_FV0A_Unfm_80, meanOcc_Prim_Unfm_80 / meanOcc_FV0A_Unfm_80);
      NormalizeVector(vecOcc_FV0C_Unfm_80, meanOcc_Prim_Unfm_80 / meanOcc_FV0C_Unfm_80);
      NormalizeVector(vecOcc_FT0A_Unfm_80, meanOcc_Prim_Unfm_80 / meanOcc_FT0A_Unfm_80);
      NormalizeVector(vecOcc_FT0C_Unfm_80, meanOcc_Prim_Unfm_80 / meanOcc_FT0C_Unfm_80);
      NormalizeVector(vecOcc_FDDA_Unfm_80, meanOcc_Prim_Unfm_80 / meanOcc_FDDA_Unfm_80);
      NormalizeVector(vecOcc_FDDC_Unfm_80, meanOcc_Prim_Unfm_80 / meanOcc_FDDC_Unfm_80);

      NormalizeVector(vecOcc_NTrack_ITS_Unfm_80, meanOcc_Prim_Unfm_80 / meanOcc_NTrack_ITS_Unfm_80);
      NormalizeVector(vecOcc_NTrack_TPC_Unfm_80, meanOcc_Prim_Unfm_80 / meanOcc_NTrack_TPC_Unfm_80);
      NormalizeVector(vecOcc_NTrack_TRD_Unfm_80, meanOcc_Prim_Unfm_80 / meanOcc_NTrack_TRD_Unfm_80);
      NormalizeVector(vecOcc_NTrack_TOF_Unfm_80, meanOcc_Prim_Unfm_80 / meanOcc_NTrack_TOF_Unfm_80);
      NormalizeVector(vecOcc_NTrackSize_Unfm_80, meanOcc_Prim_Unfm_80 / meanOcc_NTrackSize_Unfm_80);
      NormalizeVector(vecOcc_NTrackTPC_A_Unfm_80, meanOcc_Prim_Unfm_80 / meanOcc_NTrackTPC_A_Unfm_80);
      NormalizeVector(vecOcc_NTrackTPC_C_Unfm_80, meanOcc_Prim_Unfm_80 / meanOcc_NTrackTPC_C_Unfm_80);
      NormalizeVector(vecOcc_NTrackITS_TPC_Unfm_80, meanOcc_Prim_Unfm_80 / meanOcc_NTrackITS_TPC_Unfm_80);
      NormalizeVector(vecOcc_NTrackITS_TPC_A_Unfm_80, meanOcc_Prim_Unfm_80 / meanOcc_NTrackITS_TPC_A_Unfm_80);
      NormalizeVector(vecOcc_NTrackITS_TPC_C_Unfm_80, meanOcc_Prim_Unfm_80 / meanOcc_NTrackITS_TPC_C_Unfm_80);

      NormalizeVector(vecOcc_multNTracksHasITS_Unfm_80, meanOcc_Prim_Unfm_80 / meanOcc_multNTracksHasITS_Unfm_80);
      NormalizeVector(vecOcc_multNTracksHasTPC_Unfm_80, meanOcc_Prim_Unfm_80 / meanOcc_multNTracksHasTPC_Unfm_80);
      NormalizeVector(vecOcc_multNTracksHasTOF_Unfm_80, meanOcc_Prim_Unfm_80 / meanOcc_multNTracksHasTOF_Unfm_80);
      NormalizeVector(vecOcc_multNTracksHasTRD_Unfm_80, meanOcc_Prim_Unfm_80 / meanOcc_multNTracksHasTRD_Unfm_80);
      NormalizeVector(vecOcc_multNTracksITSOnly_Unfm_80, meanOcc_Prim_Unfm_80 / meanOcc_multNTracksITSOnly_Unfm_80);
      NormalizeVector(vecOcc_multNTracksTPCOnly_Unfm_80, meanOcc_Prim_Unfm_80 / meanOcc_multNTracksTPCOnly_Unfm_80);
      NormalizeVector(vecOcc_multNTracksITSTPC_Unfm_80, meanOcc_Prim_Unfm_80 / meanOcc_multNTracksITSTPC_Unfm_80);
      NormalizeVector(vecOcc_multAllTracksTPCOnly_Unfm_80, meanOcc_Prim_Unfm_80 / meanOcc_multAllTracksTPCOnly_Unfm_80);

      // Find Robust estimators
      // T0A, T0C, V0A, Prim
      GetMedianOccVect(vecRobustOcc_T0V0Prim_Unfm_80, vecRobustOcc_T0V0Prim_Unfm_80_medianPosVec,
                       vecOcc_Prim_Unfm_80, vecOcc_FV0A_Unfm_80, vecOcc_FT0A_Unfm_80, vecOcc_FT0C_Unfm_80);

      // T0A, T0C, V0A, FDD, Prim
      GetMedianOccVect(vecRobustOcc_FDDT0V0Prim_Unfm_80, vecRobustOcc_FDDT0V0Prim_Unfm_80_medianPosVec,
                       vecOcc_Prim_Unfm_80, vecOcc_FV0A_Unfm_80, vecOcc_FT0A_Unfm_80, vecOcc_FT0C_Unfm_80, vecOcc_FDDA_Unfm_80, vecOcc_FDDC_Unfm_80);

      // NTrackDet
      GetMedianOccVect(vecRobustOcc_NtrackDet_Unfm_80, vecRobustOcc_NtrackDet_Unfm_80_medianPosVec,
                       vecOcc_NTrack_ITS_Unfm_80, vecOcc_NTrack_TPC_Unfm_80, vecOcc_NTrack_TRD_Unfm_80, vecOcc_NTrack_TOF_Unfm_80);

      // multExtraTable
      GetMedianOccVect(vecRobustOcc_multTable_Unfm_80, vecRobustOcc_multTable_Unfm_80_medianPosVec,
                       vecOcc_Prim_Unfm_80, vecOcc_multNTracksHasITS_Unfm_80, vecOcc_multNTracksHasTPC_Unfm_80,
                       vecOcc_multNTracksHasTOF_Unfm_80, vecOcc_multNTracksHasTRD_Unfm_80, vecOcc_multNTracksITSOnly_Unfm_80,
                       vecOcc_multNTracksTPCOnly_Unfm_80, vecOcc_multNTracksITSTPC_Unfm_80, vecOcc_multAllTracksTPCOnly_Unfm_80,
                       vecOcc_NTrackITS_TPC_Unfm_80);

      for (const auto& vec : vecRobustOcc_T0V0Prim_Unfm_80_medianPosVec) {
        recoEvent.fill(HIST("h_RO_T0V0Prim_Unfm_80"), vec[0]);
        recoEvent.fill(HIST("h_RO_T0V0Prim_Unfm_80"), vec[1]);
      }

      for (const auto& vec : vecRobustOcc_FDDT0V0Prim_Unfm_80_medianPosVec) {
        recoEvent.fill(HIST("h_RO_FDDT0V0Prim_Unfm_80"), vec[0]);
        recoEvent.fill(HIST("h_RO_FDDT0V0Prim_Unfm_80"), vec[1]);
      }

      for (const auto& vec : vecRobustOcc_NtrackDet_Unfm_80_medianPosVec) {
        recoEvent.fill(HIST("h_RO_NtrackDet_Unfm_80"), vec[0]);
        recoEvent.fill(HIST("h_RO_NtrackDet_Unfm_80"), vec[1]);
      }

      for (const auto& vec : vecRobustOcc_multTable_Unfm_80_medianPosVec) {
        recoEvent.fill(HIST("h_RO_multTable_Unfm_80"), vec[0]);
        recoEvent.fill(HIST("h_RO_multTable_Unfm_80"), vec[1]);
      }

      GenOccsBCsList(TFlist[i], BC_TF_Map[i]);

      if (buildOccsDet) {
        GenOccsDet(vecOcc_Prim_Unfm_80,
                   vecOcc_FV0A_Unfm_80, vecOcc_FV0C_Unfm_80,
                   vecOcc_FT0A_Unfm_80, vecOcc_FT0C_Unfm_80,
                   vecOcc_FDDA_Unfm_80, vecOcc_FDDC_Unfm_80);
      }
      if (buildOccsTrackMult) {
        GenOccsTrackMult(vecOcc_NTrack_ITS_Unfm_80, vecOcc_NTrack_TPC_Unfm_80,
                         vecOcc_NTrack_TRD_Unfm_80, vecOcc_NTrack_TOF_Unfm_80,
                         vecOcc_NTrackSize_Unfm_80, vecOcc_NTrackTPC_A_Unfm_80,
                         vecOcc_NTrackTPC_C_Unfm_80, vecOcc_NTrackITS_TPC_Unfm_80,
                         vecOcc_NTrackITS_TPC_A_Unfm_80, vecOcc_NTrackITS_TPC_C_Unfm_80);
      }
      if (buildOccsMultExtra) {
        GenOccsMultExtra(vecOcc_multNTracksHasITS_Unfm_80, vecOcc_multNTracksHasTPC_Unfm_80,
                         vecOcc_multNTracksHasTOF_Unfm_80, vecOcc_multNTracksHasTRD_Unfm_80,
                         vecOcc_multNTracksITSOnly_Unfm_80, vecOcc_multNTracksTPCOnly_Unfm_80,
                         vecOcc_multNTracksITSTPC_Unfm_80, vecOcc_multAllTracksTPCOnly_Unfm_80);
      }

      if (buildOccsRobust) {
        GenOccsRobust(vecRobustOcc_T0V0Prim_Unfm_80, vecRobustOcc_FDDT0V0Prim_Unfm_80, vecRobustOcc_NtrackDet_Unfm_80, vecRobustOcc_multTable_Unfm_80);
      }
      if (buildOccsMeanDet) {
        GenOccsMeanDet(meanOcc_Prim_Unfm_80,
                       meanOcc_FV0A_Unfm_80, meanOcc_FV0C_Unfm_80,
                       meanOcc_FT0A_Unfm_80, meanOcc_FT0C_Unfm_80,
                       meanOcc_FDDA_Unfm_80, meanOcc_FDDC_Unfm_80);
      }

      if (buildOccsMeanTrkMult) {
        GenOccsMeanTrkMult(meanOcc_NTrack_ITS_Unfm_80,
                           meanOcc_NTrack_TPC_Unfm_80,
                           meanOcc_NTrack_TRD_Unfm_80,
                           meanOcc_NTrack_TOF_Unfm_80,
                           meanOcc_NTrackSize_Unfm_80,
                           meanOcc_NTrackTPC_A_Unfm_80,
                           meanOcc_NTrackTPC_C_Unfm_80,
                           meanOcc_NTrackITS_TPC_Unfm_80,
                           meanOcc_NTrackITS_TPC_A_Unfm_80,
                           meanOcc_NTrackITS_TPC_C_Unfm_80);
      }

      if (buildOccsMnMultExtra) {
        GenOccsMnMultExtra(meanOcc_multNTracksHasITS_Unfm_80, meanOcc_multNTracksHasTPC_Unfm_80,
                           meanOcc_multNTracksHasTOF_Unfm_80, meanOcc_multNTracksHasTRD_Unfm_80,
                           meanOcc_multNTracksITSOnly_Unfm_80, meanOcc_multNTracksTPCOnly_Unfm_80,
                           meanOcc_multNTracksITSTPC_Unfm_80, meanOcc_multAllTracksTPCOnly_Unfm_80);
      }

      if (buildOccsMeanRobust) {
        GenOccsMeanRobust(
          TMath::Mean(vecRobustOcc_T0V0Prim_Unfm_80.size(), vecRobustOcc_T0V0Prim_Unfm_80.data()),
          TMath::Mean(vecRobustOcc_FDDT0V0Prim_Unfm_80.size(), vecRobustOcc_FDDT0V0Prim_Unfm_80.data()),
          TMath::Mean(vecRobustOcc_NtrackDet_Unfm_80.size(), vecRobustOcc_NtrackDet_Unfm_80.data()),
          TMath::Mean(vecRobustOcc_multTable_Unfm_80.size(), vecRobustOcc_multTable_Unfm_80.data()));
      }
    }

    // Create a BC index table.
    int64_t occIDX = -1;
    int idx = -1;
    for (auto const& bc : BCs) {
      idx = -1;
      GetTimingInfo(bc, lastRun, nBCsPerTF, bcSOR, time, TFidThis, bcInTF);

      auto idxIt = std::find(TFlist.begin(), TFlist.end(), TFidThis);
      if (idxIt != TFlist.end()) {
        idx = std::distance(TFlist.begin(), idxIt);
      } else {
        LOG(error) << "DEBUG :: SEVERE :: BC  Timeframe not in the list";
      }

      auto it = std::find(BC_TF_Map[idx].begin(), BC_TF_Map[idx].end(), bc.globalIndex()); // will find the iterator where object is placed.
      if (it != BC_TF_Map[idx].end()) {                                                    // Element is in the vector
        occIDX = idx;
      } else { // Element is not in the vector
        occIDX = -1;
      }

      GenOccIndexTable(bc.globalIndex(), occIDX); // BCId, OccId
      GenBCTFinfoTable(TFidThis, bcInTF);
    }
  } // Process function ends
};

struct trackMeanOccTableProducer {
  Produces<aod::TrackMeanOccs0> GenTrackMeanOccs0;
  Produces<aod::TrackMeanOccs1> GenTrackMeanOccs1;
  Produces<aod::TrackMeanOccs2> GenTrackMeanOccs2;
  Produces<aod::TrackMeanOccs3> GenTrackMeanOccs3;
  Produces<aod::TrackMeanOccs4> GenTrackMeanOccs4;
  Produces<aod::TrackMeanOccs5> GenTrackMeanOccs5;
  Produces<aod::TrackMeanOccs6> GenTrackMeanOccs6;
  Produces<aod::TrackMeanOccs7> GenTrackMeanOccs7;
  Produces<aod::TrackMeanOccs8> GenTrackMeanOccs8;

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
  std::vector<float> Occ_Prim_Unfm_80;
  std::vector<float> Occ_FV0A_Unfm_80;
  std::vector<float> Occ_FV0C_Unfm_80;
  std::vector<float> Occ_FT0A_Unfm_80;
  std::vector<float> Occ_FT0C_Unfm_80;
  std::vector<float> Occ_FDDA_Unfm_80;
  std::vector<float> Occ_FDDC_Unfm_80;

  std::vector<float> Occ_NTrack_ITS_Unfm_80;
  std::vector<float> Occ_NTrack_TPC_Unfm_80;
  std::vector<float> Occ_NTrack_TRD_Unfm_80;
  std::vector<float> Occ_NTrack_TOF_Unfm_80;
  std::vector<float> Occ_NTrackSize_Unfm_80;
  std::vector<float> Occ_NTrackTPC_A_Unfm_80;
  std::vector<float> Occ_NTrackTPC_C_Unfm_80;
  std::vector<float> Occ_NTrackITS_TPC_Unfm_80;
  std::vector<float> Occ_NTrackITS_TPC_A_Unfm_80;
  std::vector<float> Occ_NTrackITS_TPC_C_Unfm_80;

  std::vector<float> Occ_multNTracksHasITS_Unfm_80;
  std::vector<float> Occ_multNTracksHasTPC_Unfm_80;
  std::vector<float> Occ_multNTracksHasTOF_Unfm_80;
  std::vector<float> Occ_multNTracksHasTRD_Unfm_80;
  std::vector<float> Occ_multNTracksITSOnly_Unfm_80;
  std::vector<float> Occ_multNTracksTPCOnly_Unfm_80;
  std::vector<float> Occ_multNTracksITSTPC_Unfm_80;
  std::vector<float> Occ_multAllTracksTPCOnly_Unfm_80;

  std::vector<float> OccRobust_T0V0Prim_Unfm_80;
  std::vector<float> OccRobust_FDDT0V0Prim_Unfm_80;
  std::vector<float> OccRobust_NtrackDet_Unfm_80;
  std::vector<float> OccRobust_multTable_Unfm_80;

  void init(InitContext const&)
  {
    // CCDB related part to be added later
    Occ_Prim_Unfm_80.resize(nBCinTF / 80);
    Occ_FV0A_Unfm_80.resize(nBCinTF / 80);
    Occ_FV0C_Unfm_80.resize(nBCinTF / 80);
    Occ_FT0A_Unfm_80.resize(nBCinTF / 80);
    Occ_FT0C_Unfm_80.resize(nBCinTF / 80);
    Occ_FDDA_Unfm_80.resize(nBCinTF / 80);
    Occ_FDDC_Unfm_80.resize(nBCinTF / 80);

    Occ_NTrack_ITS_Unfm_80.resize(nBCinTF / 80);
    Occ_NTrack_TPC_Unfm_80.resize(nBCinTF / 80);
    Occ_NTrack_TRD_Unfm_80.resize(nBCinTF / 80);
    Occ_NTrack_TOF_Unfm_80.resize(nBCinTF / 80);
    Occ_NTrackSize_Unfm_80.resize(nBCinTF / 80);
    Occ_NTrackTPC_A_Unfm_80.resize(nBCinTF / 80);
    Occ_NTrackTPC_C_Unfm_80.resize(nBCinTF / 80);
    Occ_NTrackITS_TPC_Unfm_80.resize(nBCinTF / 80);
    Occ_NTrackITS_TPC_A_Unfm_80.resize(nBCinTF / 80);
    Occ_NTrackITS_TPC_C_Unfm_80.resize(nBCinTF / 80);

    Occ_multNTracksHasITS_Unfm_80.resize(nBCinTF / 80);
    Occ_multNTracksHasTPC_Unfm_80.resize(nBCinTF / 80);
    Occ_multNTracksHasTOF_Unfm_80.resize(nBCinTF / 80);
    Occ_multNTracksHasTRD_Unfm_80.resize(nBCinTF / 80);
    Occ_multNTracksITSOnly_Unfm_80.resize(nBCinTF / 80);
    Occ_multNTracksTPCOnly_Unfm_80.resize(nBCinTF / 80);
    Occ_multNTracksITSTPC_Unfm_80.resize(nBCinTF / 80);
    Occ_multAllTracksTPCOnly_Unfm_80.resize(nBCinTF / 80);

    OccRobust_T0V0Prim_Unfm_80.resize(nBCinTF / 80);
    OccRobust_FDDT0V0Prim_Unfm_80.resize(nBCinTF / 80);
    OccRobust_NtrackDet_Unfm_80.resize(nBCinTF / 80);
    OccRobust_multTable_Unfm_80.resize(nBCinTF / 80);
  }

  void GetRunInfo(const int& run, int& nBCsPerTF, int64_t& bcSOR)
  {
    auto runDuration = ccdb->getRunDuration(run, true);
    int64_t tsSOR = runDuration.first;
    auto ctpx = ccdb->getForTimeStamp<std::vector<Long64_t>>("CTP/Calib/OrbitReset", tsSOR);
    int64_t tsOrbitReset = (*ctpx)[0];
    uint32_t nOrbitsPerTF = run < 534133 ? 128 : 32;
    int64_t orbitSOR = (tsSOR * 1000 - tsOrbitReset) / o2::constants::lhc::LHCOrbitMUS;
    orbitSOR = orbitSOR / nOrbitsPerTF * nOrbitsPerTF;
    bcSOR = orbitSOR * nBCsPerOrbit + customOrbitOffset * nBCsPerOrbit; // customOrbitOffset is a configurable
    nBCsPerTF = nOrbitsPerTF * nBCsPerOrbit;
  }

  template <typename T>
  void GetTimingInfo(const T& bc, int& lastRun, int32_t& nBCsPerTF, int64_t& bcSOR, uint64_t& time, int64_t& TFidThis, int& bcInTF)
  {
    int run = bc.runNumber();
    if (run != lastRun) { // update run info
      lastRun = run;
      GetRunInfo(run, nBCsPerTF, bcSOR); // update nBCsPerTF && bcSOR
    }
    // update the information
    time = bc.timestamp();
    TFidThis = (bc.globalBC() - bcSOR) / nBCsPerTF;
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
    float wR = 0;
    float R = 0;
    for (int i = binStart; i <= binEnd; i++) {
      R = m * i + c;
      wR = 125. / R;
      if (x2 == x1) {
        wR = 1.0;
      }
      sumOfBins += OccVector[i] * wR;
      weightSum += wR;
    }
    float meanOccupancy = sumOfBins / weightSum;
    return meanOccupancy;
  }

  using myCollisions = soa::Join<aod::Collisions, aod::Mults>;
  using myTracks = soa::Join<aod::Tracks, o2::aod::TracksCov, aod::TracksExtra>;
  using myBCTable = soa::Join<aod::BCsWithTimestamps, aod::OccIndexTable, aod::BCTFinfoTable>;

  using myOccsDet = soa::Join<aod::OccsBCsList, aod::OccsDet>;
  using myOccsTrackMult = soa::Join<aod::OccsBCsList, aod::OccsTrackMult>;
  using myOccsMultExtra = soa::Join<aod::OccsBCsList, aod::OccsMultExtra>;
  using myOccsRobust = soa::Join<aod::OccsBCsList, aod::OccsRobust>;

  // Process the Data
  // int dfCount = 0;
  int32_t nBCsPerTF = -1;
  int64_t bcSOR = -1;
  int lastRun = -999;

  uint64_t time = -1;
  int64_t TFidThis = -1;
  int bcInTF = -1;

  void process(myBCTable& BCs, myCollisions const& collisions, myTracks const& tracks, aod::TracksQA const& tracksQA,
               //  o2::aod::AmbiguousTracks const& ambgTracks, o2::aod::Origins const& Origins, aod::OccsBCsList const& occsBCsList, //tables only used during debugging
               myOccsDet const& occsDet, myOccsTrackMult const& occsTrackMult,
               myOccsMultExtra const& occsMultExtra, myOccsRobust const& occsRobust)
  {
    // dfCount++;
    // LOG(info) << "DEBUG :: df_" << dfCount << " :: DF_" << Origins.begin().dataframeID() << " :: collisions.size() = " << collisions.size() << " :: tracks.size() = " << tracks.size() << " :: tracksQA.size() = " << tracksQA.size()
    //           << " :: myBCTable.size() = " << BCs.size()
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

    // Defining fix value for lastRun, nBCsPerTF, bcSOR, at beginning of each dataframe to avoid wierd behaviour//They will be calculated once per DF
    lastRun = -999;
    nBCsPerTF = -999;
    bcSOR = -999;

    auto bc = BCs.begin();

    int64_t oldTFid = -1;

    int64_t oldCollisionIndex = -100;
    bool hasCollision = false;
    bool isAmbgTrack = false;
    bool LastTrackHadCollision = false;
    bool doCollisionUpdate = false;
    bool doAmbgUpdate = false;

    double Rbegin = 90., Rend = 245.;
    double Zbegin; // collision.posZ() + track.tgl()*Rbegin;
    double Zend;   // collision.posZ() + track.tgl()*Rend;
    double vdrift = 2.64;

    double dTbegin; // ((250.- TMath::Abs(Zbegin))/vdrift)/0.025;//bin
    double dTend;   // ((250.- TMath::Abs(Zend))/vdrift)/0.025;  //bin

    double bcBegin; // tGlobalBC + dTbegin;
    double bcEnd;   // tGlobalBC + dTend  ;

    int BinBCbegin;
    int BinBCend;

    for (const auto& trackQA : tracksQA) {
      auto const& track = trackQA.track_as<myTracks>();
      auto collision = collisions.begin();

      hasCollision = false;
      isAmbgTrack = false;

      if (track.collisionId() >= 0) {                   // track has collision
        collision = track.collision_as<myCollisions>(); // It will build but crash while running for tracks with track.collisionId()= -1;//ambg tracks/orphan tracks
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
        LastTrackHadCollision = true;
      }
      doCollisionUpdate = false; // default is false;
      doAmbgUpdate = false;
      if (hasCollision) {
        if (LastTrackHadCollision) {
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
          bc = collision.bc_as<myBCTable>();
        }
        if (doAmbgUpdate) {
          // to be updated later
          //  bc = collisions.iteratorAt(2).bc_as<aod::BCsWithTimestamps>();
          //  bc = ambgTracks.iteratorAt(0).bc_as<aod::BCsWithTimestamps>();
        }
        // LOG(info)<<" What happens in the case when the collision id is = -1 and it tries to obtain bc"
        GetTimingInfo(bc, lastRun, nBCsPerTF, bcSOR, time, TFidThis, bcInTF);
      }

      if (TFidThis != oldTFid) {
        oldTFid = TFidThis;
        // auto OccList = Occs.iteratorAt(bc.occId());
        if (buildTrackMeanOccs1 || buildTrackMeanOccs5) {
          auto ListOccsDet = occsDet.iteratorAt(bc.occId());
          std::copy(ListOccsDet.occ_Prim_Unfm_80().begin(), ListOccsDet.occ_Prim_Unfm_80().end(), Occ_Prim_Unfm_80.begin());
          std::copy(ListOccsDet.occ_FV0A_Unfm_80().begin(), ListOccsDet.occ_FV0A_Unfm_80().end(), Occ_FV0A_Unfm_80.begin());
          std::copy(ListOccsDet.occ_FV0C_Unfm_80().begin(), ListOccsDet.occ_FV0C_Unfm_80().end(), Occ_FV0C_Unfm_80.begin());
          std::copy(ListOccsDet.occ_FT0A_Unfm_80().begin(), ListOccsDet.occ_FT0A_Unfm_80().end(), Occ_FT0A_Unfm_80.begin());
          std::copy(ListOccsDet.occ_FT0C_Unfm_80().begin(), ListOccsDet.occ_FT0C_Unfm_80().end(), Occ_FT0C_Unfm_80.begin());
          std::copy(ListOccsDet.occ_FDDA_Unfm_80().begin(), ListOccsDet.occ_FDDA_Unfm_80().end(), Occ_FDDA_Unfm_80.begin());
          std::copy(ListOccsDet.occ_FDDC_Unfm_80().begin(), ListOccsDet.occ_FDDC_Unfm_80().end(), Occ_FDDC_Unfm_80.begin());
        }

        if (buildTrackMeanOccs2 || buildTrackMeanOccs6) {
          auto ListOccsTrackMult = occsTrackMult.iteratorAt(bc.occId());
          ;
          std::copy(ListOccsTrackMult.occ_NTrack_ITS_Unfm_80().begin(), ListOccsTrackMult.occ_NTrack_ITS_Unfm_80().end(), Occ_NTrack_ITS_Unfm_80.begin());
          std::copy(ListOccsTrackMult.occ_NTrack_TPC_Unfm_80().begin(), ListOccsTrackMult.occ_NTrack_TPC_Unfm_80().end(), Occ_NTrack_TPC_Unfm_80.begin());
          std::copy(ListOccsTrackMult.occ_NTrack_TRD_Unfm_80().begin(), ListOccsTrackMult.occ_NTrack_TRD_Unfm_80().end(), Occ_NTrack_TRD_Unfm_80.begin());
          std::copy(ListOccsTrackMult.occ_NTrack_TOF_Unfm_80().begin(), ListOccsTrackMult.occ_NTrack_TOF_Unfm_80().end(), Occ_NTrack_TOF_Unfm_80.begin());
          std::copy(ListOccsTrackMult.occ_NTrackSize_Unfm_80().begin(), ListOccsTrackMult.occ_NTrackSize_Unfm_80().end(), Occ_NTrackSize_Unfm_80.begin());
          std::copy(ListOccsTrackMult.occ_NTrackTPC_A_Unfm_80().begin(), ListOccsTrackMult.occ_NTrackTPC_A_Unfm_80().end(), Occ_NTrackTPC_A_Unfm_80.begin());
          std::copy(ListOccsTrackMult.occ_NTrackTPC_C_Unfm_80().begin(), ListOccsTrackMult.occ_NTrackTPC_C_Unfm_80().end(), Occ_NTrackTPC_C_Unfm_80.begin());
          std::copy(ListOccsTrackMult.occ_NTrackITS_TPC_Unfm_80().begin(), ListOccsTrackMult.occ_NTrackITS_TPC_Unfm_80().end(), Occ_NTrackITS_TPC_Unfm_80.begin());
          std::copy(ListOccsTrackMult.occ_NTrackITS_TPC_A_Unfm_80().begin(), ListOccsTrackMult.occ_NTrackITS_TPC_A_Unfm_80().end(), Occ_NTrackITS_TPC_A_Unfm_80.begin());
          std::copy(ListOccsTrackMult.occ_NTrackITS_TPC_C_Unfm_80().begin(), ListOccsTrackMult.occ_NTrackITS_TPC_C_Unfm_80().end(), Occ_NTrackITS_TPC_C_Unfm_80.begin());
        }

        if (buildTrackMeanOccs3 || buildTrackMeanOccs7) {
          auto ListOccsMultExtra = occsMultExtra.iteratorAt(bc.occId());
          ;
          std::copy(ListOccsMultExtra.occ_multNTracksHasITS_Unfm_80().begin(), ListOccsMultExtra.occ_multNTracksHasITS_Unfm_80().end(), Occ_multNTracksHasITS_Unfm_80.begin());
          std::copy(ListOccsMultExtra.occ_multNTracksHasTPC_Unfm_80().begin(), ListOccsMultExtra.occ_multNTracksHasTPC_Unfm_80().end(), Occ_multNTracksHasTPC_Unfm_80.begin());
          std::copy(ListOccsMultExtra.occ_multNTracksHasTOF_Unfm_80().begin(), ListOccsMultExtra.occ_multNTracksHasTOF_Unfm_80().end(), Occ_multNTracksHasTOF_Unfm_80.begin());
          std::copy(ListOccsMultExtra.occ_multNTracksHasTRD_Unfm_80().begin(), ListOccsMultExtra.occ_multNTracksHasTRD_Unfm_80().end(), Occ_multNTracksHasTRD_Unfm_80.begin());
          std::copy(ListOccsMultExtra.occ_multNTracksITSOnly_Unfm_80().begin(), ListOccsMultExtra.occ_multNTracksITSOnly_Unfm_80().end(), Occ_multNTracksITSOnly_Unfm_80.begin());
          std::copy(ListOccsMultExtra.occ_multNTracksTPCOnly_Unfm_80().begin(), ListOccsMultExtra.occ_multNTracksTPCOnly_Unfm_80().end(), Occ_multNTracksTPCOnly_Unfm_80.begin());
          std::copy(ListOccsMultExtra.occ_multNTracksITSTPC_Unfm_80().begin(), ListOccsMultExtra.occ_multNTracksITSTPC_Unfm_80().end(), Occ_multNTracksITSTPC_Unfm_80.begin());
          std::copy(ListOccsMultExtra.occ_multAllTracksTPCOnly_Unfm_80().begin(), ListOccsMultExtra.occ_multAllTracksTPCOnly_Unfm_80().end(), Occ_multAllTracksTPCOnly_Unfm_80.begin());
        }

        if (buildTrackMeanOccs4 || buildTrackMeanOccs8) {
          auto ListOccsRobust = occsRobust.iteratorAt(bc.occId());
          std::copy(ListOccsRobust.occRobust_T0V0Prim_Unfm_80().begin(), ListOccsRobust.occRobust_T0V0Prim_Unfm_80().end(), OccRobust_T0V0Prim_Unfm_80.begin());
          std::copy(ListOccsRobust.occRobust_FDDT0V0Prim_Unfm_80().begin(), ListOccsRobust.occRobust_FDDT0V0Prim_Unfm_80().end(), OccRobust_FDDT0V0Prim_Unfm_80.begin());
          std::copy(ListOccsRobust.occRobust_NtrackDet_Unfm_80().begin(), ListOccsRobust.occRobust_NtrackDet_Unfm_80().end(), OccRobust_NtrackDet_Unfm_80.begin());
          std::copy(ListOccsRobust.occRobust_multExtraTable_Unfm_80().begin(), ListOccsRobust.occRobust_multExtraTable_Unfm_80().end(), OccRobust_multTable_Unfm_80.begin());
        }
      }

      // Timebc = TGlobalBC+ΔTdrift
      // ΔTdrift=((250(cm)-abs(z))/vdrift)
      // vdrift=2.64 cm/μs
      // z=zv+tgl*Radius

      Rbegin = 90., Rend = 245.;                        // in cm
      Zbegin = collision.posZ() + track.tgl() * Rbegin; // in cm
      Zend = collision.posZ() + track.tgl() * Rend;     // in cm
      vdrift = 2.64;                                    // cm/μs

      // clip the result at 250
      if (Zbegin > 250) {
        Zbegin = 250;
      } else if (Zbegin < -250) {
        Zbegin = -250;
      }

      if (Zend > 250) {
        Zend = 250;
      } else if (Zend < -250) {
        Zend = -250;
      }

      dTbegin = ((250. - TMath::Abs(Zbegin)) / vdrift) / 0.025;
      dTend = ((250. - TMath::Abs(Zend)) / vdrift) / 0.025;

      bcBegin = bcInTF + dTbegin;
      bcEnd = bcInTF + dTend;

      BinBCbegin = bcBegin / 80;
      BinBCend = bcEnd / 80;

      GenTrackMeanOccs0(track.globalIndex());

      if (buildTrackMeanOccs1) {
        GenTrackMeanOccs1(getMeanOccupancy(BinBCbegin, BinBCend, Occ_Prim_Unfm_80),
                          getMeanOccupancy(BinBCbegin, BinBCend, Occ_FV0A_Unfm_80),
                          getMeanOccupancy(BinBCbegin, BinBCend, Occ_FV0C_Unfm_80),
                          getMeanOccupancy(BinBCbegin, BinBCend, Occ_FT0A_Unfm_80),
                          getMeanOccupancy(BinBCbegin, BinBCend, Occ_FT0C_Unfm_80),
                          getMeanOccupancy(BinBCbegin, BinBCend, Occ_FDDA_Unfm_80),
                          getMeanOccupancy(BinBCbegin, BinBCend, Occ_FDDC_Unfm_80));
      }
      if (buildTrackMeanOccs2) {
        GenTrackMeanOccs2(
          getMeanOccupancy(BinBCbegin, BinBCend, Occ_NTrack_ITS_Unfm_80),
          getMeanOccupancy(BinBCbegin, BinBCend, Occ_NTrack_TPC_Unfm_80),
          getMeanOccupancy(BinBCbegin, BinBCend, Occ_NTrack_TRD_Unfm_80),
          getMeanOccupancy(BinBCbegin, BinBCend, Occ_NTrack_TOF_Unfm_80),
          getMeanOccupancy(BinBCbegin, BinBCend, Occ_NTrackSize_Unfm_80),
          getMeanOccupancy(BinBCbegin, BinBCend, Occ_NTrackTPC_A_Unfm_80),
          getMeanOccupancy(BinBCbegin, BinBCend, Occ_NTrackTPC_C_Unfm_80),
          getMeanOccupancy(BinBCbegin, BinBCend, Occ_NTrackITS_TPC_Unfm_80),
          getMeanOccupancy(BinBCbegin, BinBCend, Occ_NTrackITS_TPC_A_Unfm_80),
          getMeanOccupancy(BinBCbegin, BinBCend, Occ_NTrackITS_TPC_C_Unfm_80));
      }
      if (buildTrackMeanOccs3) {
        GenTrackMeanOccs3(getMeanOccupancy(BinBCbegin, BinBCend, Occ_multNTracksHasITS_Unfm_80),
                          getMeanOccupancy(BinBCbegin, BinBCend, Occ_multNTracksHasTPC_Unfm_80),
                          getMeanOccupancy(BinBCbegin, BinBCend, Occ_multNTracksHasTOF_Unfm_80),
                          getMeanOccupancy(BinBCbegin, BinBCend, Occ_multNTracksHasTRD_Unfm_80),
                          getMeanOccupancy(BinBCbegin, BinBCend, Occ_multNTracksITSOnly_Unfm_80),
                          getMeanOccupancy(BinBCbegin, BinBCend, Occ_multNTracksTPCOnly_Unfm_80),
                          getMeanOccupancy(BinBCbegin, BinBCend, Occ_multNTracksITSTPC_Unfm_80),
                          getMeanOccupancy(BinBCbegin, BinBCend, Occ_multAllTracksTPCOnly_Unfm_80));
      }
      if (buildTrackMeanOccs4) {
        GenTrackMeanOccs4(getMeanOccupancy(BinBCbegin, BinBCend, OccRobust_T0V0Prim_Unfm_80),
                          getMeanOccupancy(BinBCbegin, BinBCend, OccRobust_FDDT0V0Prim_Unfm_80),
                          getMeanOccupancy(BinBCbegin, BinBCend, OccRobust_NtrackDet_Unfm_80),
                          getMeanOccupancy(BinBCbegin, BinBCend, OccRobust_multTable_Unfm_80));
      }

      if (buildTrackMeanOccs5) {
        GenTrackMeanOccs5(getWeightedMeanOccupancy(BinBCbegin, BinBCend, Occ_Prim_Unfm_80),
                          getWeightedMeanOccupancy(BinBCbegin, BinBCend, Occ_FV0A_Unfm_80),
                          getWeightedMeanOccupancy(BinBCbegin, BinBCend, Occ_FV0C_Unfm_80),
                          getWeightedMeanOccupancy(BinBCbegin, BinBCend, Occ_FT0A_Unfm_80),
                          getWeightedMeanOccupancy(BinBCbegin, BinBCend, Occ_FT0C_Unfm_80),
                          getWeightedMeanOccupancy(BinBCbegin, BinBCend, Occ_FDDA_Unfm_80),
                          getWeightedMeanOccupancy(BinBCbegin, BinBCend, Occ_FDDC_Unfm_80));
      }

      if (buildTrackMeanOccs6) {
        GenTrackMeanOccs6(
          getWeightedMeanOccupancy(BinBCbegin, BinBCend, Occ_NTrack_ITS_Unfm_80),
          getWeightedMeanOccupancy(BinBCbegin, BinBCend, Occ_NTrack_TPC_Unfm_80),
          getWeightedMeanOccupancy(BinBCbegin, BinBCend, Occ_NTrack_TRD_Unfm_80),
          getWeightedMeanOccupancy(BinBCbegin, BinBCend, Occ_NTrack_TOF_Unfm_80),
          getWeightedMeanOccupancy(BinBCbegin, BinBCend, Occ_NTrackSize_Unfm_80),
          getWeightedMeanOccupancy(BinBCbegin, BinBCend, Occ_NTrackTPC_A_Unfm_80),
          getWeightedMeanOccupancy(BinBCbegin, BinBCend, Occ_NTrackTPC_C_Unfm_80),
          getWeightedMeanOccupancy(BinBCbegin, BinBCend, Occ_NTrackITS_TPC_Unfm_80),
          getWeightedMeanOccupancy(BinBCbegin, BinBCend, Occ_NTrackITS_TPC_A_Unfm_80),
          getWeightedMeanOccupancy(BinBCbegin, BinBCend, Occ_NTrackITS_TPC_C_Unfm_80));
      }

      if (buildTrackMeanOccs7) {
        GenTrackMeanOccs7(
          getWeightedMeanOccupancy(BinBCbegin, BinBCend, Occ_multNTracksHasITS_Unfm_80),
          getWeightedMeanOccupancy(BinBCbegin, BinBCend, Occ_multNTracksHasTPC_Unfm_80),
          getWeightedMeanOccupancy(BinBCbegin, BinBCend, Occ_multNTracksHasTOF_Unfm_80),
          getWeightedMeanOccupancy(BinBCbegin, BinBCend, Occ_multNTracksHasTRD_Unfm_80),
          getWeightedMeanOccupancy(BinBCbegin, BinBCend, Occ_multNTracksITSOnly_Unfm_80),
          getWeightedMeanOccupancy(BinBCbegin, BinBCend, Occ_multNTracksTPCOnly_Unfm_80),
          getWeightedMeanOccupancy(BinBCbegin, BinBCend, Occ_multNTracksITSTPC_Unfm_80),
          getWeightedMeanOccupancy(BinBCbegin, BinBCend, Occ_multAllTracksTPCOnly_Unfm_80));
      }

      if (buildTrackMeanOccs8) {
        GenTrackMeanOccs8(getWeightedMeanOccupancy(BinBCbegin, BinBCend, OccRobust_T0V0Prim_Unfm_80),
                          getWeightedMeanOccupancy(BinBCbegin, BinBCend, OccRobust_FDDT0V0Prim_Unfm_80),
                          getWeightedMeanOccupancy(BinBCbegin, BinBCend, OccRobust_NtrackDet_Unfm_80),
                          getWeightedMeanOccupancy(BinBCbegin, BinBCend, OccRobust_multTable_Unfm_80));
      }
    } // end of trackQA loop
  } // Process function ends
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<occTableProducer>(cfgc),
    adaptAnalysisTask<trackMeanOccTableProducer>(cfgc)};
}
