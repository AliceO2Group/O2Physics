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
const int nBCinTF_grp80 = 1425;

struct occTableProducer {

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  // declare production of tables
  Produces<aod::OccIndexTable> GenOccIndexTable;
  Produces<aod::Occs> GenOccTable;

  Configurable<int> customOrbitOffset{"customOrbitOffset", 0, "customOrbitOffset for MC"};

  // Histogram registry;
  HistogramRegistry recoEvent{"recoEvent", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  void init(InitContext const&)
  {
    // Getting Info from CCDB, to be implemented Later
    recoEvent.add("h_nBCinTF", "h_nBCinTF(to check nBCinTF)", {HistType::kTH1F, {{100, 114040, 114060}}}); // 114048
    recoEvent.add("h_bcInTF", "h_bcInTF", {HistType::kTH1F, {{2000, 0, 200000}}});
    recoEvent.add("h_RO_T0V0Prim_Unfm_80", "h_RO_T0V0Prim_Unfm_80:median contributors", {HistType::kTH1F, {{12 * 2, -1, 11}}});
    recoEvent.add("h_RO_FDDT0V0Prim_Unfm_80", "h_RO_FDDT0V0Prim_Unfm_80:median contributors", {HistType::kTH1F, {{12 * 2, -1, 11}}});
    recoEvent.add("h_RO_NtrackDet_Unfm_80", "h_RO_NtrackDetITS/TPC/TRD/TOF_80:median contributors", {HistType::kTH1F, {{12 * 2, -1, 11}}});

    recoEvent.print();
  }

  template <typename T>
  bool vectorAscendingSortCheck(const T& vec)
  {
    for (uint i = 1; i < vec.size(); i++) {
      if (vec[i] < vec[i - 1]) {
        LOG(debug) << "DEBUG :: Vector unsorted at Position = " << i;
        return false;
      }
    }
    return true;
  }

  void GetNormalisedVector(const std::vector<float>& OriginalVec, std::vector<float>& NormVec, const float& scaleFactor)
  {
    NormVec.resize(OriginalVec.size());
    std::transform(OriginalVec.begin(), OriginalVec.end(), NormVec.begin(), [scaleFactor](float x) { return x * scaleFactor; });
  }

  void GetMedianOccVect(std::vector<float>& medianVector, std::vector<std::vector<int>>& medianPosVec, const std::vector<std::vector<float>>& vectors)
  {
    const int n = vectors.size();   // no of vectors/arrays
    const int size = nBCinTF_grp80; // size of each vector

    for (int i = 0; i < size; i++) {
      std::vector<std::vector<double>> data;
      double iEntry = 0;
      for (const auto& vec : vectors) {
        data.push_back({vec[i], iEntry});
        iEntry++;
      }
      // Now sort the vector;
      std::sort(data.begin(), data.end(), [](const std::vector<double>& a, const std::vector<double>& b) {
        return a[0] < b[0]; // Compare the first elements
      });

      double median;
      std::vector<int> medianPos;
      // Find the median
      if (n % 2 == 0) {
        // cout<<"Even Case ::m=";
        median = (data[(n - 1) / 2][0] + data[(n - 1) / 2 + 1][0]) / 2;
        medianPos = {static_cast<int>(data[(n - 1) / 2][1] + 0.001), static_cast<int>(data[(n - 1) / 2 + 1][1] + 0.001)};
      } else {
        // cout<<"Odd  Case ::m=";
        median = data[(n) / 2][0];
        medianPos = {static_cast<int>(data[(n) / 2][1] + 0.001)};
      }
      medianVector.push_back(median);
      medianPosVec.push_back(medianPos);
    }
  }

  int BinarySearchVector(int64_t Key, std::vector<int64_t> List, int low, int high)
  {
    while (low <= high) {
      int mid = low + (high - low) / 2;
      if (Key == List[mid]) {
        return mid;
      }

      if (Key > List[mid]) {
        low = mid + 1; // If Key is greater, ignore left  half, update the low
      } else {
        high = mid - 1;
      } // If Key is smaller, ignore right half, update the high
    }
    return -1; // Element is not present
  }

  template <typename T>
  void GetTimingInfo(const T& bc, uint64_t& time, int64_t& TFidThis, int& bcInTF)
  {
    int run = bc.runNumber();
    time = bc.timestamp();
    auto runDuration = ccdb->getRunDuration(run, true);
    int64_t tsSOR = runDuration.first;

    auto ctpx = ccdb->getForTimeStamp<std::vector<Long64_t>>("CTP/Calib/OrbitReset", tsSOR);
    int64_t tsOrbitReset = (*ctpx)[0];
    uint32_t nOrbitsPerTF = run < 534133 ? 128 : 32;
    int64_t orbitSOR = (tsSOR * 1000 - tsOrbitReset) / o2::constants::lhc::LHCOrbitMUS;
    orbitSOR = orbitSOR / nOrbitsPerTF * nOrbitsPerTF;
    int64_t bcSOR = orbitSOR * nBCsPerOrbit + customOrbitOffset * nBCsPerOrbit;
    int32_t nBCsPerTF = nOrbitsPerTF * nBCsPerOrbit;

    TFidThis = (bc.globalBC() - bcSOR) / nBCsPerTF;
    bcInTF = (bc.globalBC() - bcSOR) % nBCsPerTF;
  }

  using myCollisions = soa::Join<aod::Collisions, aod::Mults>;
  using myTracks = soa::Join<aod::Tracks, o2::aod::TracksCov, aod::TracksExtra>;

  Preslice<myTracks> TracksPerCollisionPreslice = o2::aod::track::collisionId;

  // Process the Data
  // int dfCount = 0;
  void process(o2::aod::BCsWithTimestamps const& BCs, myCollisions const& collisions, myTracks const& tracks, aod::TracksQA const& tracksQA, o2::aod::Origins const& Origins)
  {
    // dfCount++;
    // LOG(info) << "DEBUG :: df_" << dfCount << " :: DF_" << Origins.iteratorAt(0).dataframeID() << " :: collisions.size() = " << collisions.size() << " :: tracks.size() = " << tracks.size() << " :: tracksQA.size() = " << tracksQA.size() << " :: BCs.size() = " << BCs.size();

    if (collisions.size() == 0) {
      for (const auto& BC : BCs) {                            // uint i = 0; i < BCs.size() ; i++){
        GenOccIndexTable(BC.globalIndex(), -999, -999, -999); // For BCs and OccIndexTable to have same size();
      }
      return;
    }

    const int nBCinTF = 114048;         /// CCDB
    const int nBCinDrift = 114048 / 32; /// to get from ccdb in future

    // Occupancy Maps Per DF // a DF contains one or more timesFrames (TFs)
    std::unordered_map<int64_t, std::vector<int64_t>> BC_TF_Map;

    std::unordered_map<int64_t, std::vector<float>> occ_Prim_Unfm_80;
    std::unordered_map<int64_t, std::vector<float>> occ_FV0A_Unfm_80;
    std::unordered_map<int64_t, std::vector<float>> occ_FV0C_Unfm_80;
    std::unordered_map<int64_t, std::vector<float>> occ_FT0A_Unfm_80;
    std::unordered_map<int64_t, std::vector<float>> occ_FT0C_Unfm_80;
    std::unordered_map<int64_t, std::vector<float>> occ_FDDA_Unfm_80;
    std::unordered_map<int64_t, std::vector<float>> occ_FDDC_Unfm_80;

    std::unordered_map<int64_t, std::vector<float>> occ_NTrack_PVC_Unfm_80;
    std::unordered_map<int64_t, std::vector<float>> occ_NTrack_ITS_Unfm_80;
    std::unordered_map<int64_t, std::vector<float>> occ_NTrack_TPC_Unfm_80;
    std::unordered_map<int64_t, std::vector<float>> occ_NTrack_TRD_Unfm_80;
    std::unordered_map<int64_t, std::vector<float>> occ_NTrack_TOF_Unfm_80;
    std::unordered_map<int64_t, std::vector<float>> occ_NTrackSize_Unfm_80;
    std::unordered_map<int64_t, std::vector<float>> occ_NTrackTPC_A_Unfm_80;
    std::unordered_map<int64_t, std::vector<float>> occ_NTrackTPC_C_Unfm_80;
    std::unordered_map<int64_t, std::vector<float>> occ_NTrackITS_TPC_Unfm_80;
    std::unordered_map<int64_t, std::vector<float>> occ_NTrackITS_TPC_A_Unfm_80;
    std::unordered_map<int64_t, std::vector<float>> occ_NTrackITS_TPC_C_Unfm_80;
    //
    // Calculation of TFid and BCinTF
    auto bc = collisions.iteratorAt(0).bc_as<aod::BCsWithTimestamps>();
    int run = bc.runNumber();
    auto time = bc.timestamp();

    auto runDuration = ccdb->getRunDuration(run, true);
    int64_t tsSOR = runDuration.first;
    auto ctpx = ccdb->getForTimeStamp<std::vector<Long64_t>>("CTP/Calib/OrbitReset", tsSOR);
    int64_t tsOrbitReset = (*ctpx)[0];
    uint32_t nOrbitsPerTF = run < 534133 ? 128 : 32;
    int64_t orbitSOR = (tsSOR * 1000 - tsOrbitReset) / o2::constants::lhc::LHCOrbitMUS;
    orbitSOR = orbitSOR / nOrbitsPerTF * nOrbitsPerTF;
    int64_t bcSOR = orbitSOR * nBCsPerOrbit + customOrbitOffset * nBCsPerOrbit;
    int32_t nBCsPerTF = nOrbitsPerTF * nBCsPerOrbit;
    // nBCsPerITSROF = (run >= 543437 && run <= 545367) ? 594 : 198;

    int64_t TFidThis = (bc.globalBC() - bcSOR) / nBCsPerTF;
    int bcInTF = (bc.globalBC() - bcSOR) % nBCsPerTF;
    //
    int iColl = -1;
    int lastRun = -1;
    std::vector<int64_t> TFIDList;

    for (const auto& collision : collisions) {
      iColl++;

      bc = collision.bc_as<aod::BCsWithTimestamps>();
      run = bc.runNumber();
      if (run != lastRun) {
        // lastRun = run;
        time = bc.timestamp();
        runDuration = ccdb->getRunDuration(run, true);
        tsSOR = runDuration.first;
        ctpx = ccdb->getForTimeStamp<std::vector<Long64_t>>("CTP/Calib/OrbitReset", tsSOR);
        tsOrbitReset = (*ctpx)[0];
        nOrbitsPerTF = run < 534133 ? 128 : 32;
        orbitSOR = (tsSOR * 1000 - tsOrbitReset) / o2::constants::lhc::LHCOrbitMUS;
        orbitSOR = orbitSOR / nOrbitsPerTF * nOrbitsPerTF;
        bcSOR = orbitSOR * nBCsPerOrbit + customOrbitOffset * nBCsPerOrbit;
        nBCsPerTF = nOrbitsPerTF * nBCsPerOrbit;
        // nBCsPerITSROF = (run >= 543437 && run <= 545367) ? 594 : 198;

        TFidThis = (bc.globalBC() - bcSOR) / nBCsPerTF;
        bcInTF = (bc.globalBC() - bcSOR) % nBCsPerTF;

        recoEvent.fill(HIST("h_nBCinTF"), nBCsPerTF);
        recoEvent.fill(HIST("h_bcInTF"), bcInTF);

        if (nBCsPerTF > nBCinTF) {
          LOG(error) << "DEBUG :: FATAL ERROR :: nBCsPerTF > nBCinTF i.e " << nBCsPerTF << " > " << nBCinTF << " will cause crash in further process";
          return;
        }
      }

      TFIDList.push_back(TFidThis);
      BC_TF_Map[TFidThis].push_back(bc.globalIndex());

      const uint64_t collIdx = collision.globalIndex();
      const auto TracksTable_perColl = tracks.sliceBy(TracksPerCollisionPreslice, collIdx);

      int nTrack_PVC = 0;
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
        if (track.isPVContributor()) {
          nTrack_PVC++;
        } // Run 3: Has this track contributed to the collision vertex fit
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

      std::vector<float>& TFocc_Prim_Unfm_80 = occ_Prim_Unfm_80[TFidThis];
      std::vector<float>& TFocc_FV0A_Unfm_80 = occ_FV0A_Unfm_80[TFidThis];
      std::vector<float>& TFocc_FV0C_Unfm_80 = occ_FV0C_Unfm_80[TFidThis];
      std::vector<float>& TFocc_FT0A_Unfm_80 = occ_FT0A_Unfm_80[TFidThis];
      std::vector<float>& TFocc_FT0C_Unfm_80 = occ_FT0C_Unfm_80[TFidThis];
      std::vector<float>& TFocc_FDDA_Unfm_80 = occ_FDDA_Unfm_80[TFidThis];
      std::vector<float>& TFocc_FDDC_Unfm_80 = occ_FDDC_Unfm_80[TFidThis];

      std::vector<float>& TFocc_NTrack_PVC_Unfm_80 = occ_NTrack_PVC_Unfm_80[TFidThis];
      std::vector<float>& TFocc_NTrack_ITS_Unfm_80 = occ_NTrack_ITS_Unfm_80[TFidThis];
      std::vector<float>& TFocc_NTrack_TPC_Unfm_80 = occ_NTrack_TPC_Unfm_80[TFidThis];
      std::vector<float>& TFocc_NTrack_TRD_Unfm_80 = occ_NTrack_TRD_Unfm_80[TFidThis];
      std::vector<float>& TFocc_NTrack_TOF_Unfm_80 = occ_NTrack_TOF_Unfm_80[TFidThis];
      std::vector<float>& TFocc_NTrackSize_Unfm_80 = occ_NTrackSize_Unfm_80[TFidThis];
      std::vector<float>& TFocc_NTrackTPC_A_Unfm_80 = occ_NTrackTPC_A_Unfm_80[TFidThis];
      std::vector<float>& TFocc_NTrackTPC_C_Unfm_80 = occ_NTrackTPC_C_Unfm_80[TFidThis];
      std::vector<float>& TFocc_NTrackITS_TPC_Unfm_80 = occ_NTrackITS_TPC_Unfm_80[TFidThis];
      std::vector<float>& TFocc_NTrackITS_TPC_A_Unfm_80 = occ_NTrackITS_TPC_A_Unfm_80[TFidThis];
      std::vector<float>& TFocc_NTrackITS_TPC_C_Unfm_80 = occ_NTrackITS_TPC_C_Unfm_80[TFidThis];

      if (TFocc_Prim_Unfm_80.size() <= 0) {
        TFocc_Prim_Unfm_80.resize(nBCinTF / 80);
      }
      if (TFocc_FV0A_Unfm_80.size() <= 0) {
        TFocc_FV0A_Unfm_80.resize(nBCinTF / 80);
      }
      if (TFocc_FV0C_Unfm_80.size() <= 0) {
        TFocc_FV0C_Unfm_80.resize(nBCinTF / 80);
      }
      if (TFocc_FT0A_Unfm_80.size() <= 0) {
        TFocc_FT0A_Unfm_80.resize(nBCinTF / 80);
      }
      if (TFocc_FT0C_Unfm_80.size() <= 0) {
        TFocc_FT0C_Unfm_80.resize(nBCinTF / 80);
      }
      if (TFocc_FDDA_Unfm_80.size() <= 0) {
        TFocc_FDDA_Unfm_80.resize(nBCinTF / 80);
      }
      if (TFocc_FDDC_Unfm_80.size() <= 0) {
        TFocc_FDDC_Unfm_80.resize(nBCinTF / 80);
      }

      if (TFocc_NTrack_PVC_Unfm_80.size() <= 0) {
        TFocc_NTrack_PVC_Unfm_80.resize(nBCinTF / 80);
      }
      if (TFocc_NTrack_ITS_Unfm_80.size() <= 0) {
        TFocc_NTrack_ITS_Unfm_80.resize(nBCinTF / 80);
      }
      if (TFocc_NTrack_TPC_Unfm_80.size() <= 0) {
        TFocc_NTrack_TPC_Unfm_80.resize(nBCinTF / 80);
      }
      if (TFocc_NTrack_TRD_Unfm_80.size() <= 0) {
        TFocc_NTrack_TRD_Unfm_80.resize(nBCinTF / 80);
      }
      if (TFocc_NTrack_TOF_Unfm_80.size() <= 0) {
        TFocc_NTrack_TOF_Unfm_80.resize(nBCinTF / 80);
      }
      if (TFocc_NTrackSize_Unfm_80.size() <= 0) {
        TFocc_NTrackSize_Unfm_80.resize(nBCinTF / 80);
      }
      if (TFocc_NTrackTPC_A_Unfm_80.size() <= 0) {
        TFocc_NTrackTPC_A_Unfm_80.resize(nBCinTF / 80);
      }
      if (TFocc_NTrackTPC_C_Unfm_80.size() <= 0) {
        TFocc_NTrackTPC_C_Unfm_80.resize(nBCinTF / 80);
      }
      if (TFocc_NTrackITS_TPC_Unfm_80.size() <= 0) {
        TFocc_NTrackITS_TPC_Unfm_80.resize(nBCinTF / 80);
      }
      if (TFocc_NTrackITS_TPC_A_Unfm_80.size() <= 0) {
        TFocc_NTrackITS_TPC_A_Unfm_80.resize(nBCinTF / 80);
      }
      if (TFocc_NTrackITS_TPC_C_Unfm_80.size() <= 0) {
        TFocc_NTrackITS_TPC_C_Unfm_80.resize(nBCinTF / 80);
      }

      // current collision bin in 80/160 grouping.
      int bin80_0 = bcInTF / 80;
      // int bin160_0=bcInTF/160;

      // float fbin80_0 =float(bcInTF)/80;
      // float fbin160_0=float(bcInTF)/160;

      ushort fNumContrib = collision.numContrib();
      float fMultFV0A = collision.multFV0A(), fMultFV0C = collision.multFV0C();
      float fMultFT0A = collision.multFT0A(), fMultFT0C = collision.multFT0C();
      float fMultFDDA = collision.multFDDA(), fMultFDDC = collision.multFDDC();
      int fNTrack_PVC = nTrack_PVC, fNTrack_ITS = nTrack_ITS;
      int fNTrack_TPC = nTrack_TPC, fNTrack_TRD = nTrack_TRD;
      int fNTrack_TOF = nTrack_TOF;
      int fNTrackTPC_A = nTrackTPC_A, fNTrackTPC_C = nTrackTPC_C;
      int fNTrackITS_TPC = nTrackITS_TPC, fNTrackSize = TracksTable_perColl.size();
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

        TFocc_NTrack_PVC_Unfm_80[(bin80_0 + deltaBin) % (nBCinTF / 80)] += fNTrack_PVC * 1;
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
      }
    }
    // collision Loop is over

    std::vector<int64_t> SortedTFIDList = TFIDList;
    std::sort(SortedTFIDList.begin(), SortedTFIDList.end());
    auto last = std::unique(SortedTFIDList.begin(), SortedTFIDList.end());
    SortedTFIDList.erase(last, SortedTFIDList.end());

    int totalSize = 0;
    for (auto x : SortedTFIDList) {
      totalSize += BC_TF_Map[x].size();
      // check if the BCs are already sorted or not
      if (!vectorAscendingSortCheck(BC_TF_Map[x])) {
        LOG(debug) << "DEBUG :: ERROR :: BCs are not sorted";
      }
    }

    if (totalSize != collisions.size()) {
      LOG(debug) << "DEBUG :: ERROR :: filled TF list and collision size mismatch ::  filledTF_Size = " << totalSize << " :: " << collisions.size() << " = collisions.size()";
    }

    // Fill the Producer
    int ikey = -1;
    std::vector<int64_t> keyList;
    for (const auto& pair : occ_Prim_Unfm_80) {
      ikey++;
      int key = pair.first;
      keyList.push_back(key);
    }
    // Sort the keyList; //keys are not in ascending order while filling, but TFs are in ascending order
    std::sort(keyList.begin(), keyList.end());
    int keyCounter = -1;
    // Write the table
    for (const auto& key : keyList) {
      const std::vector<float>& vecOcc_Prim_Unfm_80 = occ_Prim_Unfm_80[key];
      const std::vector<float>& vecOcc_FV0A_Unfm_80 = occ_FV0A_Unfm_80[key];
      const std::vector<float>& vecOcc_FV0C_Unfm_80 = occ_FV0C_Unfm_80[key];
      const std::vector<float>& vecOcc_FT0A_Unfm_80 = occ_FT0A_Unfm_80[key];
      const std::vector<float>& vecOcc_FT0C_Unfm_80 = occ_FT0C_Unfm_80[key];
      const std::vector<float>& vecOcc_FDDA_Unfm_80 = occ_FDDA_Unfm_80[key];
      const std::vector<float>& vecOcc_FDDC_Unfm_80 = occ_FDDC_Unfm_80[key];

      const std::vector<float>& vecOcc_NTrack_PVC_Unfm_80 = occ_NTrack_PVC_Unfm_80[key];
      const std::vector<float>& vecOcc_NTrack_ITS_Unfm_80 = occ_NTrack_ITS_Unfm_80[key];
      const std::vector<float>& vecOcc_NTrack_TPC_Unfm_80 = occ_NTrack_TPC_Unfm_80[key];
      const std::vector<float>& vecOcc_NTrack_TRD_Unfm_80 = occ_NTrack_TRD_Unfm_80[key];
      const std::vector<float>& vecOcc_NTrack_TOF_Unfm_80 = occ_NTrack_TOF_Unfm_80[key];
      const std::vector<float>& vecOcc_NTrackSize_Unfm_80 = occ_NTrackSize_Unfm_80[key];
      const std::vector<float>& vecOcc_NTrackTPC_A_Unfm_80 = occ_NTrackTPC_A_Unfm_80[key];
      const std::vector<float>& vecOcc_NTrackTPC_C_Unfm_80 = occ_NTrackTPC_C_Unfm_80[key];
      const std::vector<float>& vecOcc_NTrackITS_TPC_Unfm_80 = occ_NTrackITS_TPC_Unfm_80[key];
      const std::vector<float>& vecOcc_NTrackITS_TPC_A_Unfm_80 = occ_NTrackITS_TPC_A_Unfm_80[key];
      const std::vector<float>& vecOcc_NTrackITS_TPC_C_Unfm_80 = occ_NTrackITS_TPC_C_Unfm_80[key];

      float meanOcc_Prim_Unfm_80 = TMath::Mean(vecOcc_Prim_Unfm_80.size(), vecOcc_Prim_Unfm_80.data());
      float meanOcc_FV0A_Unfm_80 = TMath::Mean(vecOcc_FV0A_Unfm_80.size(), vecOcc_FV0A_Unfm_80.data());
      float meanOcc_FV0C_Unfm_80 = TMath::Mean(vecOcc_FV0C_Unfm_80.size(), vecOcc_FV0C_Unfm_80.data());
      float meanOcc_FT0A_Unfm_80 = TMath::Mean(vecOcc_FT0A_Unfm_80.size(), vecOcc_FT0A_Unfm_80.data());
      float meanOcc_FT0C_Unfm_80 = TMath::Mean(vecOcc_FT0C_Unfm_80.size(), vecOcc_FT0C_Unfm_80.data());
      float meanOcc_FDDA_Unfm_80 = TMath::Mean(vecOcc_FDDA_Unfm_80.size(), vecOcc_FDDA_Unfm_80.data());
      float meanOcc_FDDC_Unfm_80 = TMath::Mean(vecOcc_FDDC_Unfm_80.size(), vecOcc_FDDC_Unfm_80.data());
      float meanOcc_NTrack_PVC_Unfm_80 = TMath::Mean(vecOcc_NTrack_PVC_Unfm_80.size(), vecOcc_NTrack_PVC_Unfm_80.data());
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

      // Get and Store the normalised vectors
      std::vector<float> vecOccNorm_Prim_Unfm_80;
      GetNormalisedVector(vecOcc_Prim_Unfm_80, vecOccNorm_Prim_Unfm_80, meanOcc_Prim_Unfm_80 / meanOcc_Prim_Unfm_80);
      std::vector<float> vecOccNorm_FV0A_Unfm_80;
      GetNormalisedVector(vecOcc_FV0A_Unfm_80, vecOccNorm_FV0A_Unfm_80, meanOcc_Prim_Unfm_80 / meanOcc_FV0A_Unfm_80);
      std::vector<float> vecOccNorm_FV0C_Unfm_80;
      GetNormalisedVector(vecOcc_FV0C_Unfm_80, vecOccNorm_FV0C_Unfm_80, meanOcc_Prim_Unfm_80 / meanOcc_FV0C_Unfm_80);
      std::vector<float> vecOccNorm_FT0A_Unfm_80;
      GetNormalisedVector(vecOcc_FT0A_Unfm_80, vecOccNorm_FT0A_Unfm_80, meanOcc_Prim_Unfm_80 / meanOcc_FT0A_Unfm_80);
      std::vector<float> vecOccNorm_FT0C_Unfm_80;
      GetNormalisedVector(vecOcc_FT0C_Unfm_80, vecOccNorm_FT0C_Unfm_80, meanOcc_Prim_Unfm_80 / meanOcc_FT0C_Unfm_80);
      std::vector<float> vecOccNorm_FDDA_Unfm_80;
      GetNormalisedVector(vecOcc_FDDA_Unfm_80, vecOccNorm_FDDA_Unfm_80, meanOcc_Prim_Unfm_80 / meanOcc_FDDA_Unfm_80);
      std::vector<float> vecOccNorm_FDDC_Unfm_80;
      GetNormalisedVector(vecOcc_FDDC_Unfm_80, vecOccNorm_FDDC_Unfm_80, meanOcc_Prim_Unfm_80 / meanOcc_FDDC_Unfm_80);
      std::vector<float> vecOccNorm_NTrack_PVC_Unfm_80;
      GetNormalisedVector(vecOcc_NTrack_PVC_Unfm_80, vecOccNorm_NTrack_PVC_Unfm_80, meanOcc_Prim_Unfm_80 / meanOcc_NTrack_PVC_Unfm_80);
      std::vector<float> vecOccNorm_NTrack_ITS_Unfm_80;
      GetNormalisedVector(vecOcc_NTrack_ITS_Unfm_80, vecOccNorm_NTrack_ITS_Unfm_80, meanOcc_Prim_Unfm_80 / meanOcc_NTrack_ITS_Unfm_80);
      std::vector<float> vecOccNorm_NTrack_TPC_Unfm_80;
      GetNormalisedVector(vecOcc_NTrack_TPC_Unfm_80, vecOccNorm_NTrack_TPC_Unfm_80, meanOcc_Prim_Unfm_80 / meanOcc_NTrack_TPC_Unfm_80);
      std::vector<float> vecOccNorm_NTrack_TRD_Unfm_80;
      GetNormalisedVector(vecOcc_NTrack_TRD_Unfm_80, vecOccNorm_NTrack_TRD_Unfm_80, meanOcc_Prim_Unfm_80 / meanOcc_NTrack_TRD_Unfm_80);
      std::vector<float> vecOccNorm_NTrack_TOF_Unfm_80;
      GetNormalisedVector(vecOcc_NTrack_TOF_Unfm_80, vecOccNorm_NTrack_TOF_Unfm_80, meanOcc_Prim_Unfm_80 / meanOcc_NTrack_TOF_Unfm_80);
      std::vector<float> vecOccNorm_NTrackSize_Unfm_80;
      GetNormalisedVector(vecOcc_NTrackSize_Unfm_80, vecOccNorm_NTrackSize_Unfm_80, meanOcc_Prim_Unfm_80 / meanOcc_NTrackSize_Unfm_80);
      std::vector<float> vecOccNorm_NTrackTPC_A_Unfm_80;
      GetNormalisedVector(vecOcc_NTrackTPC_A_Unfm_80, vecOccNorm_NTrackTPC_A_Unfm_80, meanOcc_Prim_Unfm_80 / meanOcc_NTrackTPC_A_Unfm_80);
      std::vector<float> vecOccNorm_NTrackTPC_C_Unfm_80;
      GetNormalisedVector(vecOcc_NTrackTPC_C_Unfm_80, vecOccNorm_NTrackTPC_C_Unfm_80, meanOcc_Prim_Unfm_80 / meanOcc_NTrackTPC_C_Unfm_80);
      std::vector<float> vecOccNorm_NTrackITS_TPC_Unfm_80;
      GetNormalisedVector(vecOcc_NTrackITS_TPC_Unfm_80, vecOccNorm_NTrackITS_TPC_Unfm_80, meanOcc_Prim_Unfm_80 / meanOcc_NTrackITS_TPC_Unfm_80);
      std::vector<float> vecOccNorm_NTrackITS_TPC_A_Unfm_80;
      GetNormalisedVector(vecOcc_NTrackITS_TPC_A_Unfm_80, vecOccNorm_NTrackITS_TPC_A_Unfm_80, meanOcc_Prim_Unfm_80 / meanOcc_NTrackITS_TPC_A_Unfm_80);
      std::vector<float> vecOccNorm_NTrackITS_TPC_C_Unfm_80;
      GetNormalisedVector(vecOcc_NTrackITS_TPC_C_Unfm_80, vecOccNorm_NTrackITS_TPC_C_Unfm_80, meanOcc_Prim_Unfm_80 / meanOcc_NTrackITS_TPC_C_Unfm_80);

      // Find Robust estimators
      // T0A, T0C, V0A, Prim
      std::vector<float> vecRobustOcc_T0V0Prim_Unfm_80;
      std::vector<std::vector<int>> vecRobustOcc_T0V0Prim_Unfm_80_medianPosVec;
      GetMedianOccVect(vecRobustOcc_T0V0Prim_Unfm_80, vecRobustOcc_T0V0Prim_Unfm_80_medianPosVec,
                       {vecOccNorm_Prim_Unfm_80, vecOccNorm_FV0A_Unfm_80, vecOccNorm_FT0A_Unfm_80, vecOccNorm_FT0C_Unfm_80});

      // T0A, T0C, V0A, FDD, Prim
      std::vector<float> vecRobustOcc_FDDT0V0Prim_Unfm_80;
      std::vector<std::vector<int>> vecRobustOcc_FDDT0V0Prim_Unfm_80_medianPosVec;
      GetMedianOccVect(vecRobustOcc_FDDT0V0Prim_Unfm_80, vecRobustOcc_FDDT0V0Prim_Unfm_80_medianPosVec,
                       {vecOccNorm_Prim_Unfm_80, vecOccNorm_FV0A_Unfm_80, vecOccNorm_FT0A_Unfm_80, vecOccNorm_FT0C_Unfm_80, vecOccNorm_FDDA_Unfm_80, vecOccNorm_FDDC_Unfm_80});

      // NTrackDet
      std::vector<float> vecRobustOcc_NtrackDet_Unfm_80;
      std::vector<std::vector<int>> vecRobustOcc_NtrackDet_Unfm_80_medianPosVec;
      GetMedianOccVect(vecRobustOcc_NtrackDet_Unfm_80, vecRobustOcc_NtrackDet_Unfm_80_medianPosVec,
                       {vecOccNorm_NTrack_ITS_Unfm_80, vecOccNorm_NTrack_TPC_Unfm_80, vecOccNorm_NTrack_TRD_Unfm_80, vecOccNorm_NTrack_TOF_Unfm_80});

      for (const auto& vec : vecRobustOcc_T0V0Prim_Unfm_80_medianPosVec) {
        for (const auto& val : vec) {
          recoEvent.fill(HIST("h_RO_T0V0Prim_Unfm_80"), val);
        }
      }

      for (const auto& vec : vecRobustOcc_FDDT0V0Prim_Unfm_80_medianPosVec) {
        for (const auto& val : vec) {
          recoEvent.fill(HIST("h_RO_FDDT0V0Prim_Unfm_80"), val);
        }
      }

      for (const auto& vec : vecRobustOcc_NtrackDet_Unfm_80_medianPosVec) {
        for (const auto& val : vec) {
          recoEvent.fill(HIST("h_RO_NtrackDet_Unfm_80"), val);
        }
      }

      keyCounter++;
      GenOccTable(
        key, BC_TF_Map[key],
        vecOccNorm_Prim_Unfm_80,
        vecOccNorm_FV0A_Unfm_80,
        vecOccNorm_FV0C_Unfm_80,
        vecOccNorm_FT0A_Unfm_80,
        vecOccNorm_FT0C_Unfm_80,
        vecOccNorm_FDDA_Unfm_80,
        vecOccNorm_FDDC_Unfm_80,
        vecOccNorm_NTrack_PVC_Unfm_80,
        vecOccNorm_NTrack_ITS_Unfm_80,
        vecOccNorm_NTrack_TPC_Unfm_80,
        vecOccNorm_NTrack_TRD_Unfm_80,
        vecOccNorm_NTrack_TOF_Unfm_80,
        vecOccNorm_NTrackSize_Unfm_80,
        vecOccNorm_NTrackTPC_A_Unfm_80,
        vecOccNorm_NTrackTPC_C_Unfm_80,
        vecOccNorm_NTrackITS_TPC_Unfm_80,
        vecOccNorm_NTrackITS_TPC_A_Unfm_80,
        vecOccNorm_NTrackITS_TPC_C_Unfm_80,
        vecRobustOcc_T0V0Prim_Unfm_80,
        vecRobustOcc_FDDT0V0Prim_Unfm_80,
        vecRobustOcc_NtrackDet_Unfm_80,
        TMath::Mean(vecOcc_Prim_Unfm_80.size(), vecOcc_Prim_Unfm_80.data()),
        TMath::Mean(vecOcc_FV0A_Unfm_80.size(), vecOcc_FV0A_Unfm_80.data()),
        TMath::Mean(vecOcc_FV0C_Unfm_80.size(), vecOcc_FV0C_Unfm_80.data()),
        TMath::Mean(vecOcc_FT0A_Unfm_80.size(), vecOcc_FT0A_Unfm_80.data()),
        TMath::Mean(vecOcc_FT0C_Unfm_80.size(), vecOcc_FT0C_Unfm_80.data()),
        TMath::Mean(vecOcc_FDDA_Unfm_80.size(), vecOcc_FDDA_Unfm_80.data()),
        TMath::Mean(vecOcc_FDDC_Unfm_80.size(), vecOcc_FDDC_Unfm_80.data()),
        TMath::Mean(vecOcc_NTrack_PVC_Unfm_80.size(), vecOcc_NTrack_PVC_Unfm_80.data()),
        TMath::Mean(vecOcc_NTrack_ITS_Unfm_80.size(), vecOcc_NTrack_ITS_Unfm_80.data()),
        TMath::Mean(vecOcc_NTrack_TPC_Unfm_80.size(), vecOcc_NTrack_TPC_Unfm_80.data()),
        TMath::Mean(vecOcc_NTrack_TRD_Unfm_80.size(), vecOcc_NTrack_TRD_Unfm_80.data()),
        TMath::Mean(vecOcc_NTrack_TOF_Unfm_80.size(), vecOcc_NTrack_TOF_Unfm_80.data()),
        TMath::Mean(vecOcc_NTrackSize_Unfm_80.size(), vecOcc_NTrackSize_Unfm_80.data()),
        TMath::Mean(vecOcc_NTrackTPC_A_Unfm_80.size(), vecOcc_NTrackTPC_A_Unfm_80.data()),
        TMath::Mean(vecOcc_NTrackTPC_C_Unfm_80.size(), vecOcc_NTrackTPC_C_Unfm_80.data()),
        TMath::Mean(vecOcc_NTrackITS_TPC_Unfm_80.size(), vecOcc_NTrackITS_TPC_Unfm_80.data()),
        TMath::Mean(vecOcc_NTrackITS_TPC_A_Unfm_80.size(), vecOcc_NTrackITS_TPC_A_Unfm_80.data()),
        TMath::Mean(vecOcc_NTrackITS_TPC_C_Unfm_80.size(), vecOcc_NTrackITS_TPC_C_Unfm_80.data()),
        TMath::Mean(vecRobustOcc_T0V0Prim_Unfm_80.size(), vecRobustOcc_T0V0Prim_Unfm_80.data()),
        TMath::Mean(vecRobustOcc_FDDT0V0Prim_Unfm_80.size(), vecRobustOcc_FDDT0V0Prim_Unfm_80.data()),
        TMath::Mean(vecRobustOcc_NtrackDet_Unfm_80.size(), vecRobustOcc_NtrackDet_Unfm_80.data()));
    }

    if (static_cast<int>(ikey + 1) != static_cast<int>(SortedTFIDList.size())) {
      LOG(error) << "DEBUG :: ERROR :: #keys and SortedTFIdList have different sizes " << static_cast<int>(ikey + 1) << " :: " << SortedTFIDList.size();
      return;
    }

    // Create a BC index table.
    for (auto const& bc : BCs) {
      GetTimingInfo(bc, time, TFidThis, bcInTF);
      // check if this BC was used
      int pos = BinarySearchVector(bc.globalIndex(), BC_TF_Map[TFidThis], 0, BC_TF_Map[TFidThis].size() - 1);
      int64_t occIDX = -1;
      if (pos == -1) {
        occIDX = -1;
      } else if (pos >= 0) {
        occIDX = BinarySearchVector(TFidThis, keyList, 0, keyList.size() - 1);
      } else {
        LOG(error) << "DEBUG :: ERROR :: For BC, occIDX = -2 (BC pos = -2 when searched) :";
        occIDX = -2;
      }
      GenOccIndexTable(bc.globalIndex(), occIDX, TFidThis, bcInTF);
    }
  } // Process function ends
};

struct trackMeanOccTableProducer {
  Produces<aod::TrackMeanOccs> meanOccTable;

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  // Configurables
  Configurable<int> customOrbitOffset{"customOrbitOffset", 0, "customOrbitOffset for MC"};

  void init(InitContext const&)
  {
    // CCDB related part to be added later
  }

  template <typename T>
  void GetTimingInfo(const T& bc, uint64_t& time, int64_t& TFidThis, int& bcInTF)
  {
    int run = bc.runNumber();
    time = bc.timestamp();
    auto runDuration = ccdb->getRunDuration(run, true);
    int64_t tsSOR = runDuration.first;

    auto ctpx = ccdb->getForTimeStamp<std::vector<Long64_t>>("CTP/Calib/OrbitReset", tsSOR);
    int64_t tsOrbitReset = (*ctpx)[0];
    uint32_t nOrbitsPerTF = run < 534133 ? 128 : 32;
    int64_t orbitSOR = (tsSOR * 1000 - tsOrbitReset) / o2::constants::lhc::LHCOrbitMUS;
    orbitSOR = orbitSOR / nOrbitsPerTF * nOrbitsPerTF;
    int64_t bcSOR = orbitSOR * nBCsPerOrbit + customOrbitOffset * nBCsPerOrbit;
    int32_t nBCsPerTF = nOrbitsPerTF * nBCsPerOrbit;

    TFidThis = (bc.globalBC() - bcSOR) / nBCsPerTF;
    bcInTF = (bc.globalBC() - bcSOR) % nBCsPerTF;
  }

  template <typename T, typename U>
  int BinarySearchTable(T Key, U Table, int low, int high)
  {
    while (low <= high) {
      int mid = low + (high - low) / 2;
      auto midElement = Table.iteratorAt(mid).trackId();
      if (Key == midElement) {
        return mid;
      }
      if (Key > midElement) {
        low = mid + 1; // If Key is greater, ignore left  half, update the low
      } else {
        high = mid - 1; // If Key is smaller, ignore right half, update the high
      }
    }
    return -1; // Element is not present
  }

  template <typename T, typename U>
  int FindInTable(T key, U Table)
  {
    return BinarySearchTable(key, Table, 0, Table.size() - 1);
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
  using myBCTable = soa::Join<aod::BCsWithTimestamps, aod::OccIndexTable>;

  // Process the Data
  // int dfCount = 0;
  void process(myBCTable const& BCs, myCollisions const& collisions, myTracks const& tracks,
               aod::TracksQA const& tracksQA, o2::aod::AmbiguousTracks const& ambgTracks, o2::aod::Origins const& Origins,
               aod::Occs const& Occs, aod::OccIndexTable const& OccIndices)
  {
    // dfCount++;
    // if(dfCount > 10) {return;}
    // LOG(info) << "DEBUG :: df_" << dfCount << " :: DF_" << Origins.iteratorAt(0).dataframeID() << " :: collisions.size() = " << collisions.size() << " :: tracks.size() = " << tracks.size() << " :: tracksQA.size() = " << tracksQA.size()
    //           << " :: myBCTable.size() = " << BCs.size()
    //           << " :: Occs.size() = " << Occs.size()
    //           << " :: OccIndices() = " << OccIndices.size();
    if (collisions.size() == 0) {
      return;
    }
    //
    auto bc = collisions.iteratorAt(0).bc_as<myBCTable>();
    int run = bc.runNumber();
    auto time = bc.timestamp();

    auto runDuration = ccdb->getRunDuration(run, true);
    int64_t tsSOR = runDuration.first;
    auto ctpx = ccdb->getForTimeStamp<std::vector<Long64_t>>("CTP/Calib/OrbitReset", tsSOR);
    int64_t tsOrbitReset = (*ctpx)[0];
    uint32_t nOrbitsPerTF = run < 534133 ? 128 : 32;
    int64_t orbitSOR = (tsSOR * 1000 - tsOrbitReset) / o2::constants::lhc::LHCOrbitMUS;
    orbitSOR = orbitSOR / nOrbitsPerTF * nOrbitsPerTF;
    int64_t bcSOR = orbitSOR * nBCsPerOrbit + customOrbitOffset * nBCsPerOrbit;
    int32_t nBCsPerTF = nOrbitsPerTF * nBCsPerOrbit;
    // nBCsPerITSROF = (run >= 543437 && run <= 545367) ? 594 : 198;

    int64_t TFidThis = (bc.globalBC() - bcSOR) / nBCsPerTF;
    int bcInTF = (bc.globalBC() - bcSOR) % nBCsPerTF;
    //

    std::vector<float> Occ_Prim_Unfm_80;
    std::vector<float> Occ_FV0A_Unfm_80;
    std::vector<float> Occ_FV0C_Unfm_80;
    std::vector<float> Occ_FT0A_Unfm_80;
    std::vector<float> Occ_FT0C_Unfm_80;
    std::vector<float> Occ_FDDA_Unfm_80;
    std::vector<float> Occ_FDDC_Unfm_80;

    std::vector<float> Occ_NTrack_PVC_Unfm_80;
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

    std::vector<float> OccRobust_T0V0Prim_Unfm_80;
    std::vector<float> OccRobust_FDDT0V0Prim_Unfm_80;
    std::vector<float> OccRobust_NtrackDet_Unfm_80;

    int64_t oldTFid = -1;

    int64_t oldCollisionIndex = -100;
    int CollisionId_negError = 0;
    int CollisionId_posError = 0;
    int CollisionId_NoError = 0;
    int nAmbgTracks = 0;
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
      // auto track = trackQA.track_as<myTracks>;  // dereferncing not working // only option is either use iterator way or Table slicing

      auto const& track = tracks.iteratorAt(trackQA.trackId()); // Find the track in track table
      if (track.globalIndex() != trackQA.trackId()) {
        LOG(error) << "DEBUG :: ERROR :: Track and TrackQA Mismatch";
      }

      auto const& collision = collisions.iteratorAt(track.collisionId()); // Find the collision in collision table
      if (track.collisionId() != collision.globalIndex()) {
        LOG(error) << "DEBUG :: ERROR :: track collId and collID Mismatch";
      }

      // Checking out of the range errors
      if (trackQA.trackId() < 0 || tracks.size() <= trackQA.trackId()) {
        LOG(error) << "DEBUG :: ERROR :: trackQA has index out of scope :: trackQA.trackId() = " << trackQA.trackId() << " :: track.collisionId() = " << track.collisionId() << " :: track.signed1Pt() = " << track.signed1Pt();
      }
      hasCollision = false;
      isAmbgTrack = false;
      if (track.collisionId() < 0 || collisions.size() <= track.collisionId()) {
        if (track.collisionId() < 0) {
          CollisionId_negError++;
          // check if it is an ambiguous track
          int ambgPos = FindInTable(track.globalIndex(), ambgTracks);
          if (ambgPos >= 0 && (ambgTracks.iteratorAt(ambgPos).trackId() == track.globalIndex())) {
            nAmbgTracks++;
            isAmbgTrack = true;
          } else {
            LOG(error) << "DEBUG :: ERROR :: Not an ambiguous track either ::";
          }
        } else {
          CollisionId_posError++;
        }
      } else {
        CollisionId_NoError++;
        hasCollision = true;
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
          } else {
            doCollisionUpdate = true;
          } // if collisions are different
        } else {
          doCollisionUpdate = true;
        } // LastTrackWasAmbiguous
      } else if (isAmbgTrack) {
        doAmbgUpdate = true;
        // To be updated later
        //  if(LastTrackIsAmbg){
        //    if( haveSameInfo ) { doAmbgUpdate = false;}
        //    else              { doAmbgUpdate = true; }
        //  }
        //  else { doAmbgUpdate = true;} //Last track had Collisions
      }

      if (doAmbgUpdate) {
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
        GetTimingInfo(bc, time, TFidThis, bcInTF);
      }

      if (TFidThis != oldTFid) {
        oldTFid = TFidThis;
        auto OccList = Occs.iteratorAt(bc.occId());
        Occ_Prim_Unfm_80 = std::vector<float>(OccList.occ_Prim_Unfm_80().begin(), OccList.occ_Prim_Unfm_80().end());
        Occ_FV0A_Unfm_80 = std::vector<float>(OccList.occ_FV0A_Unfm_80().begin(), OccList.occ_FV0A_Unfm_80().end());
        Occ_FV0C_Unfm_80 = std::vector<float>(OccList.occ_FV0C_Unfm_80().begin(), OccList.occ_FV0C_Unfm_80().end());
        Occ_FT0A_Unfm_80 = std::vector<float>(OccList.occ_FT0A_Unfm_80().begin(), OccList.occ_FT0A_Unfm_80().end());
        Occ_FT0C_Unfm_80 = std::vector<float>(OccList.occ_FT0C_Unfm_80().begin(), OccList.occ_FT0C_Unfm_80().end());
        Occ_FDDA_Unfm_80 = std::vector<float>(OccList.occ_FDDA_Unfm_80().begin(), OccList.occ_FDDA_Unfm_80().end());
        Occ_FDDC_Unfm_80 = std::vector<float>(OccList.occ_FDDC_Unfm_80().begin(), OccList.occ_FDDC_Unfm_80().end());

        Occ_NTrack_PVC_Unfm_80 = std::vector<float>(OccList.occ_NTrack_PVC_Unfm_80().begin(), OccList.occ_NTrack_PVC_Unfm_80().end());
        Occ_NTrack_ITS_Unfm_80 = std::vector<float>(OccList.occ_NTrack_ITS_Unfm_80().begin(), OccList.occ_NTrack_ITS_Unfm_80().end());
        Occ_NTrack_TPC_Unfm_80 = std::vector<float>(OccList.occ_NTrack_TPC_Unfm_80().begin(), OccList.occ_NTrack_TPC_Unfm_80().end());
        Occ_NTrack_TRD_Unfm_80 = std::vector<float>(OccList.occ_NTrack_TRD_Unfm_80().begin(), OccList.occ_NTrack_TRD_Unfm_80().end());
        Occ_NTrack_TOF_Unfm_80 = std::vector<float>(OccList.occ_NTrack_TOF_Unfm_80().begin(), OccList.occ_NTrack_TOF_Unfm_80().end());
        Occ_NTrackSize_Unfm_80 = std::vector<float>(OccList.occ_NTrackSize_Unfm_80().begin(), OccList.occ_NTrackSize_Unfm_80().end());
        Occ_NTrackTPC_A_Unfm_80 = std::vector<float>(OccList.occ_NTrackTPC_A_Unfm_80().begin(), OccList.occ_NTrackTPC_A_Unfm_80().end());
        Occ_NTrackTPC_C_Unfm_80 = std::vector<float>(OccList.occ_NTrackTPC_C_Unfm_80().begin(), OccList.occ_NTrackTPC_C_Unfm_80().end());
        Occ_NTrackITS_TPC_Unfm_80 = std::vector<float>(OccList.occ_NTrackITS_TPC_Unfm_80().begin(), OccList.occ_NTrackITS_TPC_Unfm_80().end());
        Occ_NTrackITS_TPC_A_Unfm_80 = std::vector<float>(OccList.occ_NTrackITS_TPC_A_Unfm_80().begin(), OccList.occ_NTrackITS_TPC_A_Unfm_80().end());
        Occ_NTrackITS_TPC_C_Unfm_80 = std::vector<float>(OccList.occ_NTrackITS_TPC_C_Unfm_80().begin(), OccList.occ_NTrackITS_TPC_C_Unfm_80().end());

        OccRobust_T0V0Prim_Unfm_80 = std::vector<float>(OccList.occRobust_T0V0Prim_Unfm_80().begin(), OccList.occRobust_T0V0Prim_Unfm_80().end());
        OccRobust_FDDT0V0Prim_Unfm_80 = std::vector<float>(OccList.occRobust_FDDT0V0Prim_Unfm_80().begin(), OccList.occRobust_FDDT0V0Prim_Unfm_80().end());
        OccRobust_NtrackDet_Unfm_80 = std::vector<float>(OccList.occRobust_NtrackDet_Unfm_80().begin(), OccList.occRobust_NtrackDet_Unfm_80().end());
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

      meanOccTable(
        track.globalIndex(),
        getMeanOccupancy(BinBCbegin, BinBCend, Occ_Prim_Unfm_80),
        getMeanOccupancy(BinBCbegin, BinBCend, Occ_FV0A_Unfm_80),
        getMeanOccupancy(BinBCbegin, BinBCend, Occ_FV0C_Unfm_80),
        getMeanOccupancy(BinBCbegin, BinBCend, Occ_FT0A_Unfm_80),
        getMeanOccupancy(BinBCbegin, BinBCend, Occ_FT0C_Unfm_80),
        getMeanOccupancy(BinBCbegin, BinBCend, Occ_FDDA_Unfm_80),
        getMeanOccupancy(BinBCbegin, BinBCend, Occ_FDDC_Unfm_80),
        getMeanOccupancy(BinBCbegin, BinBCend, Occ_NTrack_PVC_Unfm_80),
        getMeanOccupancy(BinBCbegin, BinBCend, Occ_NTrack_ITS_Unfm_80),
        getMeanOccupancy(BinBCbegin, BinBCend, Occ_NTrack_TPC_Unfm_80),
        getMeanOccupancy(BinBCbegin, BinBCend, Occ_NTrack_TRD_Unfm_80),
        getMeanOccupancy(BinBCbegin, BinBCend, Occ_NTrack_TOF_Unfm_80),
        getMeanOccupancy(BinBCbegin, BinBCend, Occ_NTrackSize_Unfm_80),
        getMeanOccupancy(BinBCbegin, BinBCend, Occ_NTrackTPC_A_Unfm_80),
        getMeanOccupancy(BinBCbegin, BinBCend, Occ_NTrackTPC_C_Unfm_80),
        getMeanOccupancy(BinBCbegin, BinBCend, Occ_NTrackITS_TPC_Unfm_80),
        getMeanOccupancy(BinBCbegin, BinBCend, Occ_NTrackITS_TPC_A_Unfm_80),
        getMeanOccupancy(BinBCbegin, BinBCend, Occ_NTrackITS_TPC_C_Unfm_80),
        getMeanOccupancy(BinBCbegin, BinBCend, OccRobust_T0V0Prim_Unfm_80),
        getMeanOccupancy(BinBCbegin, BinBCend, OccRobust_FDDT0V0Prim_Unfm_80),
        getMeanOccupancy(BinBCbegin, BinBCend, OccRobust_NtrackDet_Unfm_80),
        getWeightedMeanOccupancy(BinBCbegin, BinBCend, Occ_Prim_Unfm_80),
        getWeightedMeanOccupancy(BinBCbegin, BinBCend, Occ_FV0A_Unfm_80),
        getWeightedMeanOccupancy(BinBCbegin, BinBCend, Occ_FV0C_Unfm_80),
        getWeightedMeanOccupancy(BinBCbegin, BinBCend, Occ_FT0A_Unfm_80),
        getWeightedMeanOccupancy(BinBCbegin, BinBCend, Occ_FT0C_Unfm_80),
        getWeightedMeanOccupancy(BinBCbegin, BinBCend, Occ_FDDA_Unfm_80),
        getWeightedMeanOccupancy(BinBCbegin, BinBCend, Occ_FDDC_Unfm_80),
        getWeightedMeanOccupancy(BinBCbegin, BinBCend, Occ_NTrack_PVC_Unfm_80),
        getWeightedMeanOccupancy(BinBCbegin, BinBCend, Occ_NTrack_ITS_Unfm_80),
        getWeightedMeanOccupancy(BinBCbegin, BinBCend, Occ_NTrack_TPC_Unfm_80),
        getWeightedMeanOccupancy(BinBCbegin, BinBCend, Occ_NTrack_TRD_Unfm_80),
        getWeightedMeanOccupancy(BinBCbegin, BinBCend, Occ_NTrack_TOF_Unfm_80),
        getWeightedMeanOccupancy(BinBCbegin, BinBCend, Occ_NTrackSize_Unfm_80),
        getWeightedMeanOccupancy(BinBCbegin, BinBCend, Occ_NTrackTPC_A_Unfm_80),
        getWeightedMeanOccupancy(BinBCbegin, BinBCend, Occ_NTrackTPC_C_Unfm_80),
        getWeightedMeanOccupancy(BinBCbegin, BinBCend, Occ_NTrackITS_TPC_Unfm_80),
        getWeightedMeanOccupancy(BinBCbegin, BinBCend, Occ_NTrackITS_TPC_A_Unfm_80),
        getWeightedMeanOccupancy(BinBCbegin, BinBCend, Occ_NTrackITS_TPC_C_Unfm_80),
        getWeightedMeanOccupancy(BinBCbegin, BinBCend, OccRobust_T0V0Prim_Unfm_80),
        getWeightedMeanOccupancy(BinBCbegin, BinBCend, OccRobust_FDDT0V0Prim_Unfm_80),
        getWeightedMeanOccupancy(BinBCbegin, BinBCend, OccRobust_NtrackDet_Unfm_80));
    } // end of trackQA loop
  } // Process function ends
  // PROCESS_SWITCH(trackMeanOccTableProducer, processTrackOccTable, "Calculate Mean Occupancy for tracks having entry for trackQA & collisions", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<occTableProducer>(cfgc),
    adaptAnalysisTask<trackMeanOccTableProducer>(cfgc)};
}
