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
/// \file   occupancyTableProducer.cxx
/// \brief  Occupancy Table Producer : TPC PID - Calibration
///         Occupancy calculater using tracks which have entry for collision and trackQA tables
///         Ambg tracks were not used
/// \author Rahul Verma (rahul.verma@iitb.ac.in) :: Marian I Ivanov (marian.ivanov@cern.ch)

#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/OccupancyTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/LHCConstants.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TMath.h>

#include <sys/types.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

int32_t nBCsPerOrbit = o2::constants::lhc::LHCMaxBunches;
// const int nBCinTFgrp80 = 1425;
// nOrbitsPerTF = run < 534133 ? 128 : 32;
// for 128 => nBCsPerTF = 456192 , for 32 => nBCsPerTF = 114048
// const int nBCinTF = 114048;         /// CCDB value // to be obtained from CCDB in future
const int nBCinDrift = 114048 / 32; /// to get from ccdb in future

template <typename T, std::size_t N>
void sortVectorOfArray(std::vector<std::array<T, N>>& myVector, const int& myIDX)
{
  std::sort(myVector.begin(), myVector.end(), [myIDX](const std::array<T, N>& a, const std::array<T, N>& b) {
    return a[myIDX] < b[myIDX]; // sort at the required index
  });
}

template <typename T, std::size_t N>
void checkUniqueness(const std::vector<std::array<T, N>>& myVector, const int& myIDX)
{
  for (size_t i = 1; i < myVector.size(); i++) {
    if (myVector[i][myIDX] <= myVector[i - 1][myIDX]) {
      LOG(error) << "Duplicate Entries while creating Index tables :: (vec[" << i << "][" << myIDX << "]) " << myVector[i][myIDX] << " >= " << myVector[i - 1][myIDX] << " (vec[" << i - 1 << "][" << myIDX << "])";
    }
  }
}

struct OccupancyTableProducer {

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  // declare production of tables
  Produces<aod::BCTFinfoTable> genBCTFinfoTable;
  // 0
  Produces<aod::OccIndexTable> genOccIndexTable;
  Produces<aod::OccsBCsList> genOccsBCsList;
  // 1
  Produces<aod::OccsPrim> genOccsPrim;
  Produces<aod::OccsMeanPrim> genOccsMeanPrim;
  // 2
  Produces<aod::OccsT0V0> genOccsT0V0;
  Produces<aod::OccsMeanT0V0> genOccsMeanT0V0;
  ;
  Produces<aod::ORT0V0Prim> genORT0V0Prim;
  Produces<aod::OMRT0V0Prim> genOccsMeanRobustT0V0Prim;
  // 3
  Produces<aod::OccsFDD> genOccsFDD;
  Produces<aod::OccsMeanFDD> genOccsMeanFDD;
  Produces<aod::ORFDDT0V0Prim> genORFDDT0V0Prim;
  Produces<aod::OMRFDDT0V0Prim> genOccsMeanRobustFDDT0V0Prim;
  // 4
  Produces<aod::OccsNTrackDet> genOccsNTrackDet;
  Produces<aod::OccsMeanNTrkDet> genOccsMeanNTrkDet;
  Produces<aod::ORNtrackDet> genORNtrackDet;
  Produces<aod::OMRNtrackDet> genOccsMeanRobustNtrackDet;
  // 5
  Produces<aod::OccsMultExtra> genOccsMultExtra;
  Produces<aod::OccsMnMultExtra> genOccsMnMultExtra;
  Produces<aod::ORMultExtra> genORMultExtra;
  Produces<aod::OMRMultExtra> genOccsMeanRobustMultExtraTable;

  Configurable<int> customOrbitOffset{"customOrbitOffset", 0, "customOrbitOffset for MC"};
  Configurable<int> bcGrouping{"bcGrouping", 80, "bcGrouping of BCs"};
  Configurable<int> nBCinTF{"nBCinTF", 114048, "nBCinTF"};
  Configurable<int> occVecArraySize{"occVecArraySize", 10, "occVecArraySize"};

  Configurable<int> cfgNOrbitsPerTF0RunValue{"cfgNOrbitsPerTF0RunValue", 534133, "cfgNOrbitsPerTF0RunValue"};
  Configurable<int> cfgNOrbitsPerTF1TrueValue{"cfgNOrbitsPerTF1TrueValue", 128, "cfgNOrbitsPerTF1TrueValue"};
  Configurable<int> cfgNOrbitsPerTF2FalseValue{"cfgNOrbitsPerTF2FalseValue", 32, "ccfgNOrbitsPerTF2FalseValue"};

  // declare production of tables
  Configurable<bool> buildOnlyOccsPrim{"buildOnlyOccsPrim", true, "builder of table OccsPrim"};
  Configurable<bool> buildOnlyOccsT0V0Prim{"buildOnlyOccsT0V0Prim", true, "builder of table OccsT0V0Prim"};
  Configurable<bool> buildOnlyOccsFDDT0V0Prim{"buildOnlyOccsFDDT0V0Prim", true, "builder of table OccsFDDT0V0Prim"};
  Configurable<bool> buildOnlyOccsNtrackDet{"buildOnlyOccsNtrackDet", true, "builder of table OccsNtrackDet"};
  Configurable<bool> buildOnlyOccsMultExtra{"buildOnlyOccsMultExtra", true, "builder of table OccsMultExtra"};
  Configurable<bool> buildFullOccTableProducer{"buildFullOccTableProducer", true, "builder of all Occupancy Tables"};

  Configurable<bool> buildFlag00OccTable{"buildFlag00OccTable", true, "switch of table Occ Table"};
  Configurable<bool> buildFlag01OccMeanTable{"buildFlag01OccMeanTable", true, "switch of table Occ MeanTable"};
  Configurable<bool> buildFlag02OccRobustTable{"buildFlag02OccRobustTable", true, "switch of table Occ RobustTable"};
  Configurable<bool> buildFlag03OccMeanRobustTable{"buildFlag03OccMeanRobustTable", true, "switch of table Occ MeanRobustTable"};

  // Histogram registry;
  HistogramRegistry recoEvent{"recoEvent", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry occupancyQA{"occupancyQA", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  std::vector<int> tfList;
  std::vector<std::vector<int64_t>> bcTFMap;

  std::vector<std::vector<float>> occPrimUnfm80;
  std::vector<std::vector<float>> occFV0AUnfm80;
  std::vector<std::vector<float>> occFV0CUnfm80;
  std::vector<std::vector<float>> occFT0AUnfm80;
  std::vector<std::vector<float>> occFT0CUnfm80;

  std::vector<std::vector<float>> occFDDAUnfm80;
  std::vector<std::vector<float>> occFDDCUnfm80;

  std::vector<std::vector<float>> occNTrackITSUnfm80;
  std::vector<std::vector<float>> occNTrackTPCUnfm80;
  std::vector<std::vector<float>> occNTrackTRDUnfm80;
  std::vector<std::vector<float>> occNTrackTOFUnfm80;
  std::vector<std::vector<float>> occNTrackSizeUnfm80;
  std::vector<std::vector<float>> occNTrackTPCAUnfm80;
  std::vector<std::vector<float>> occNTrackTPCCUnfm80;
  std::vector<std::vector<float>> occNTrackITSTPCUnfm80;
  std::vector<std::vector<float>> occNTrackITSTPCAUnfm80;
  std::vector<std::vector<float>> occNTrackITSTPCCUnfm80;

  std::vector<std::vector<float>> occMultNTracksHasITSUnfm80;
  std::vector<std::vector<float>> occMultNTracksHasTPCUnfm80;
  std::vector<std::vector<float>> occMultNTracksHasTOFUnfm80;
  std::vector<std::vector<float>> occMultNTracksHasTRDUnfm80;
  std::vector<std::vector<float>> occMultNTracksITSOnlyUnfm80;
  std::vector<std::vector<float>> occMultNTracksTPCOnlyUnfm80;
  std::vector<std::vector<float>> occMultNTracksITSTPCUnfm80;
  std::vector<std::vector<float>> occMultAllTracksTPCOnlyUnfm80;

  std::vector<float> vecRobustOccT0V0PrimUnfm80;
  std::vector<float> vecRobustOccFDDT0V0PrimUnfm80;
  std::vector<float> vecRobustOccNtrackDetUnfm80;
  std::vector<float> vecRobustOccmultTableUnfm80;
  std::vector<std::array<int, 2>> vecRobustOccT0V0PrimUnfm80medianPosVec; // Median => one for odd and two for even entries
  std::vector<std::array<int, 2>> vecRobustOccFDDT0V0PrimUnfm80medianPosVec;
  std::vector<std::array<int, 2>> vecRobustOccNtrackDetUnfm80medianPosVec;
  std::vector<std::array<int, 2>> vecRobustOccmultTableUnfm80medianPosVec;

  std::vector<bool> processStatus;
  Configurable<uint> processStatusSize{"processStatusSize", 10, "processStatusSize"};
  void init(InitContext const&)
  {
    processStatus.resize(processStatusSize);
    for (uint i = 0; i < processStatusSize; i++) {
      processStatus[i] = false;
    }

    // outer vector resized at runtime
    tfList.resize(occVecArraySize);
    bcTFMap.resize(occVecArraySize);

    if (buildFullOccTableProducer || buildOnlyOccsPrim || buildOnlyOccsT0V0Prim || buildOnlyOccsFDDT0V0Prim || buildOnlyOccsNtrackDet || buildOnlyOccsMultExtra) {
      occPrimUnfm80.resize(occVecArraySize);
    }
    if (buildFullOccTableProducer || buildOnlyOccsT0V0Prim || buildOnlyOccsFDDT0V0Prim) {
      occFV0AUnfm80.resize(occVecArraySize);
      occFV0CUnfm80.resize(occVecArraySize);
      occFT0AUnfm80.resize(occVecArraySize);
      occFT0CUnfm80.resize(occVecArraySize);
    }
    if (buildFullOccTableProducer || buildOnlyOccsFDDT0V0Prim) {
      occFDDAUnfm80.resize(occVecArraySize);
      occFDDCUnfm80.resize(occVecArraySize);
    }
    if (buildFullOccTableProducer || buildOnlyOccsNtrackDet) {
      occNTrackITSUnfm80.resize(occVecArraySize);
      occNTrackTPCUnfm80.resize(occVecArraySize);
      occNTrackTRDUnfm80.resize(occVecArraySize);
      occNTrackTOFUnfm80.resize(occVecArraySize);
      occNTrackSizeUnfm80.resize(occVecArraySize);
      occNTrackTPCAUnfm80.resize(occVecArraySize);
      occNTrackTPCCUnfm80.resize(occVecArraySize);
      occNTrackITSTPCAUnfm80.resize(occVecArraySize);
      occNTrackITSTPCCUnfm80.resize(occVecArraySize);
    }
    if (buildFullOccTableProducer || buildOnlyOccsNtrackDet || buildOnlyOccsMultExtra) {
      occNTrackITSTPCUnfm80.resize(occVecArraySize);
    }
    if (buildFullOccTableProducer || buildOnlyOccsMultExtra) {
      occMultNTracksHasITSUnfm80.resize(occVecArraySize);
      occMultNTracksHasTPCUnfm80.resize(occVecArraySize);
      occMultNTracksHasTOFUnfm80.resize(occVecArraySize);
      occMultNTracksHasTRDUnfm80.resize(occVecArraySize);
      occMultNTracksITSOnlyUnfm80.resize(occVecArraySize);
      occMultNTracksTPCOnlyUnfm80.resize(occVecArraySize);
      occMultNTracksITSTPCUnfm80.resize(occVecArraySize);
      occMultAllTracksTPCOnlyUnfm80.resize(occVecArraySize);
    }

    for (int i = 0; i < occVecArraySize; i++) {
      bcTFMap[i].resize(nBCinTF / bcGrouping);
      if (buildFullOccTableProducer || buildOnlyOccsPrim || buildOnlyOccsT0V0Prim || buildOnlyOccsFDDT0V0Prim || buildOnlyOccsNtrackDet || buildOnlyOccsMultExtra) {
        occPrimUnfm80[i].resize(nBCinTF / bcGrouping);
      }
      if (buildFullOccTableProducer || buildOnlyOccsT0V0Prim || buildOnlyOccsFDDT0V0Prim) {
        occFV0AUnfm80[i].resize(nBCinTF / bcGrouping);
        occFV0CUnfm80[i].resize(nBCinTF / bcGrouping);
        occFT0AUnfm80[i].resize(nBCinTF / bcGrouping);
        occFT0CUnfm80[i].resize(nBCinTF / bcGrouping);
      }
      if (buildFullOccTableProducer || buildOnlyOccsFDDT0V0Prim) {
        occFDDAUnfm80[i].resize(nBCinTF / bcGrouping);
        occFDDCUnfm80[i].resize(nBCinTF / bcGrouping);
      }
      if (buildFullOccTableProducer || buildOnlyOccsNtrackDet) {
        occNTrackITSUnfm80[i].resize(nBCinTF / bcGrouping);
        occNTrackTPCUnfm80[i].resize(nBCinTF / bcGrouping);
        occNTrackTRDUnfm80[i].resize(nBCinTF / bcGrouping);
        occNTrackTOFUnfm80[i].resize(nBCinTF / bcGrouping);
        occNTrackSizeUnfm80[i].resize(nBCinTF / bcGrouping);
        occNTrackTPCAUnfm80[i].resize(nBCinTF / bcGrouping);
        occNTrackTPCCUnfm80[i].resize(nBCinTF / bcGrouping);
        occNTrackITSTPCAUnfm80[i].resize(nBCinTF / bcGrouping);
        occNTrackITSTPCCUnfm80[i].resize(nBCinTF / bcGrouping);
      }
      if (buildFullOccTableProducer || buildOnlyOccsNtrackDet || buildOnlyOccsMultExtra) {
        occNTrackITSTPCUnfm80[i].resize(nBCinTF / bcGrouping);
      }
      if (buildFullOccTableProducer || buildOnlyOccsMultExtra) {
        occMultNTracksHasITSUnfm80[i].resize(nBCinTF / bcGrouping);
        occMultNTracksHasTPCUnfm80[i].resize(nBCinTF / bcGrouping);
        occMultNTracksHasTOFUnfm80[i].resize(nBCinTF / bcGrouping);
        occMultNTracksHasTRDUnfm80[i].resize(nBCinTF / bcGrouping);
        occMultNTracksITSOnlyUnfm80[i].resize(nBCinTF / bcGrouping);
        occMultNTracksTPCOnlyUnfm80[i].resize(nBCinTF / bcGrouping);
        occMultNTracksITSTPCUnfm80[i].resize(nBCinTF / bcGrouping);
        occMultAllTracksTPCOnlyUnfm80[i].resize(nBCinTF / bcGrouping);
      }
    }

    if (buildFullOccTableProducer || buildOnlyOccsT0V0Prim || buildFlag02OccRobustTable || buildFlag03OccMeanRobustTable) {
      vecRobustOccT0V0PrimUnfm80.resize(nBCinTF / bcGrouping);
      vecRobustOccT0V0PrimUnfm80medianPosVec.resize(nBCinTF / bcGrouping); // Median => one for odd and two for even entries
    }
    if (buildFullOccTableProducer || buildOnlyOccsFDDT0V0Prim || buildFlag02OccRobustTable || buildFlag03OccMeanRobustTable) {
      vecRobustOccFDDT0V0PrimUnfm80.resize(nBCinTF / bcGrouping);
      vecRobustOccFDDT0V0PrimUnfm80medianPosVec.resize(nBCinTF / bcGrouping);
    }
    if (buildFullOccTableProducer || buildOnlyOccsNtrackDet || buildFlag02OccRobustTable || buildFlag03OccMeanRobustTable) {
      vecRobustOccNtrackDetUnfm80.resize(nBCinTF / bcGrouping);
      vecRobustOccNtrackDetUnfm80medianPosVec.resize(nBCinTF / bcGrouping);
    }
    if (buildFullOccTableProducer || buildOnlyOccsMultExtra || buildFlag02OccRobustTable || buildFlag03OccMeanRobustTable) {
      vecRobustOccmultTableUnfm80.resize(nBCinTF / bcGrouping);
      vecRobustOccmultTableUnfm80medianPosVec.resize(nBCinTF / bcGrouping);
    }

    // Getting Info from CCDB, to be implemented Later
    recoEvent.add("h_nBCinTF", "h_nBCinTF(to check nBCinTF)", {HistType::kTH1F, {{100, 114040, 114060}}}); // 114048
    recoEvent.add("h_bcInTF", "h_bcInTF", {HistType::kTH1F, {{2000, 0, 200000}}});
    recoEvent.add("h_RO_T0V0PrimUnfm80", "h_RO_T0V0PrimUnfm80:median contributors", {HistType::kTH1F, {{12 * 2, -1, 11}}});
    recoEvent.add("h_RO_FDDT0V0PrimUnfm80", "h_RO_FDDT0V0PrimUnfm80:median contributors", {HistType::kTH1F, {{12 * 2, -1, 11}}});
    recoEvent.add("h_RO_NtrackDetUnfm80", "h_RO_NtrackDetITS/TPC/TRD/TOF_80:median contributors", {HistType::kTH1F, {{12 * 2, -1, 11}}});
    recoEvent.add("h_RO_multTableUnfm80", "h_RO_multTableExtra_80:median contributors", {HistType::kTH1F, {{12 * 2, -1, 11}}});
    occupancyQA.add("h_TF_in_DataFrame", "h_TF_in_DataFrame", kTH1F, {{50, -1, 49}});
    occupancyQA.add("h_DFcount_Lvl0", "h_DFcount_Lvl0", kTH1F, {{1, 0, 1}});
    occupancyQA.add("h_DFcount_Lvl1", "h_DFcount_Lvl1", kTH1F, {{1, 0, 1}});
    occupancyQA.add("h_DFcount_Lvl2", "h_DFcount_Lvl2", kTH1F, {{1, 0, 1}});

    recoEvent.print();
    occupancyQA.print();
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
      int two = 2;
      // Find the median
      if (n % two == 0) {
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
    uint32_t nOrbitsPerTF = run < cfgNOrbitsPerTF0RunValue ? cfgNOrbitsPerTF1TrueValue : cfgNOrbitsPerTF2FalseValue;
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

  using MyTracks = aod::Tracks;

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
  //_________________________Full Exection Function_______________________________________________________________________________________

  enum ProcessTags {
    kProcessOnlyBCTFinfoTable = 0,
    kProcessOnlyOccPrim,
    kProcessOnlyOccT0V0Prim,
    kProcessOnlyOccFDDT0V0Prim,
    kProcessOnlyOccNtrackDet,
    kProcessOnlyOccMultExtra,
    kProcessFullOccTableProducer
  };

  static constexpr std::string_view ProcessNames[]{
    "processOnlyBCTFinfoTable",
    "processOnlyOccPrimUnfm",
    "processOnlyOccT0V0PrimUnfm",
    "processOnlyOccFDDT0V0PrimUnfm",
    "processOnlyOccNtrackDet",
    "processOnlyOccMultExtra",
    "processFullOccTableProduer"};

  enum FillMode {
    checkTableMode = 0,
    doNotFill,
    fillOccTable,
    fillMeanOccTable,
    fillOccRobustTable,
    fillOccMeanRobustTable
  };

  template <typename B, typename C>
  void executeCollisionCheckAndBCprocessing(B const& BCs, C const& collisions, bool& collisionsSizeIsZero)
  {
    if (collisions.size() == 0) {
      for (const auto& BC : BCs) { // For BCs and OccIndexTable to have same size for joining
        getTimingInfo(BC, lastRun, nBCsPerTF, bcSOR, time, tfIdThis, bcInTF);
        genBCTFinfoTable(tfIdThis, bcInTF);
        genOccIndexTable(BC.globalIndex(), -999); // BCId, OccId
      }
      collisionsSizeIsZero = true;
    }
  }

  int processTimeCounter = 0;
  template <int processMode, int tableMode, int meanTableMode, int robustTableMode, int meanRobustTableMode, typename B, typename C, typename T>
  void executeOccProducerProcessing(B const& BCs, C const& collisions, T const& tracks)
  {
    if (tableMode == checkTableMode) {
      if (buildFlag00OccTable) {
        executeOccProducerProcessing<processMode, fillOccTable, meanTableMode, robustTableMode, meanRobustTableMode>(BCs, collisions, tracks);
      } else {
        executeOccProducerProcessing<processMode, doNotFill, meanTableMode, robustTableMode, meanRobustTableMode>(BCs, collisions, tracks);
      }
    }
    if constexpr (tableMode == checkTableMode) {
      return;
    }

    if (meanTableMode == checkTableMode) {
      if (buildFlag01OccMeanTable) {
        executeOccProducerProcessing<processMode, tableMode, fillMeanOccTable, robustTableMode, meanRobustTableMode>(BCs, collisions, tracks);
      } else {
        executeOccProducerProcessing<processMode, tableMode, doNotFill, robustTableMode, meanRobustTableMode>(BCs, collisions, tracks);
      }
    }
    if constexpr (meanTableMode == checkTableMode) {
      return;
    }

    if (robustTableMode == checkTableMode) {
      if (buildFlag02OccRobustTable) {
        executeOccProducerProcessing<processMode, tableMode, meanTableMode, fillOccRobustTable, meanRobustTableMode>(BCs, collisions, tracks);
      } else {
        executeOccProducerProcessing<processMode, tableMode, meanTableMode, doNotFill, meanRobustTableMode>(BCs, collisions, tracks);
      }
    }
    if constexpr (robustTableMode == checkTableMode) {
      return;
    }

    if (meanRobustTableMode == checkTableMode) {
      if (buildFlag03OccMeanRobustTable) {
        executeOccProducerProcessing<processMode, tableMode, meanTableMode, robustTableMode, fillOccMeanRobustTable>(BCs, collisions, tracks);
      } else {
        executeOccProducerProcessing<processMode, tableMode, meanTableMode, robustTableMode, doNotFill>(BCs, collisions, tracks);
      }
    }
    if constexpr (meanRobustTableMode == checkTableMode) {
      return;
    }

    if constexpr (tableMode == checkTableMode || meanTableMode == checkTableMode || robustTableMode == checkTableMode || meanRobustTableMode == checkTableMode) {
      return;
    } else {

      occupancyQA.fill(HIST("h_DFcount_Lvl2"), 0.5);

      // Initialisze the vectors components to zero
      tfIDX = 0;
      tfCounted = 0;
      for (int i = 0; i < occVecArraySize; i++) {
        tfList[i] = -1;
        bcTFMap[i].clear(); // list of BCs used in one time frame;
        if constexpr (processMode == kProcessFullOccTableProducer || processMode == kProcessOnlyOccPrim || processMode == kProcessOnlyOccT0V0Prim || processMode == kProcessOnlyOccFDDT0V0Prim || processMode == kProcessOnlyOccNtrackDet || processMode == kProcessOnlyOccMultExtra) {
          std::fill(occPrimUnfm80[i].begin(), occPrimUnfm80[i].end(), 0.);
        }
        if constexpr (processMode == kProcessFullOccTableProducer || processMode == kProcessOnlyOccT0V0Prim || processMode == kProcessOnlyOccFDDT0V0Prim) {
          std::fill(occFV0AUnfm80[i].begin(), occFV0AUnfm80[i].end(), 0.);
          std::fill(occFV0CUnfm80[i].begin(), occFV0CUnfm80[i].end(), 0.);
          std::fill(occFT0AUnfm80[i].begin(), occFT0AUnfm80[i].end(), 0.);
          std::fill(occFT0CUnfm80[i].begin(), occFT0CUnfm80[i].end(), 0.);
        }
        if constexpr (processMode == kProcessFullOccTableProducer || processMode == kProcessOnlyOccFDDT0V0Prim) {
          std::fill(occFDDAUnfm80[i].begin(), occFDDAUnfm80[i].end(), 0.);
          std::fill(occFDDCUnfm80[i].begin(), occFDDCUnfm80[i].end(), 0.);
        }
        if constexpr (processMode == kProcessFullOccTableProducer || processMode == kProcessOnlyOccNtrackDet) {
          std::fill(occNTrackITSUnfm80[i].begin(), occNTrackITSUnfm80[i].end(), 0.);
          std::fill(occNTrackTPCUnfm80[i].begin(), occNTrackTPCUnfm80[i].end(), 0.);
          std::fill(occNTrackTRDUnfm80[i].begin(), occNTrackTRDUnfm80[i].end(), 0.);
          std::fill(occNTrackTOFUnfm80[i].begin(), occNTrackTOFUnfm80[i].end(), 0.);
          std::fill(occNTrackSizeUnfm80[i].begin(), occNTrackSizeUnfm80[i].end(), 0.);
          std::fill(occNTrackTPCAUnfm80[i].begin(), occNTrackTPCAUnfm80[i].end(), 0.);
          std::fill(occNTrackTPCCUnfm80[i].begin(), occNTrackTPCCUnfm80[i].end(), 0.);
          std::fill(occNTrackITSTPCAUnfm80[i].begin(), occNTrackITSTPCAUnfm80[i].end(), 0.);
          std::fill(occNTrackITSTPCCUnfm80[i].begin(), occNTrackITSTPCCUnfm80[i].end(), 0.);
        }
        if constexpr (processMode == kProcessFullOccTableProducer || processMode == kProcessOnlyOccNtrackDet || processMode == kProcessOnlyOccMultExtra) {
          std::fill(occNTrackITSTPCUnfm80[i].begin(), occNTrackITSTPCUnfm80[i].end(), 0.);
        }
        if constexpr (processMode == kProcessFullOccTableProducer || processMode == kProcessOnlyOccMultExtra) {
          std::fill(occMultNTracksHasITSUnfm80[i].begin(), occMultNTracksHasITSUnfm80[i].end(), 0.);
          std::fill(occMultNTracksHasTPCUnfm80[i].begin(), occMultNTracksHasTPCUnfm80[i].end(), 0.);
          std::fill(occMultNTracksHasTOFUnfm80[i].begin(), occMultNTracksHasTOFUnfm80[i].end(), 0.);
          std::fill(occMultNTracksHasTRDUnfm80[i].begin(), occMultNTracksHasTRDUnfm80[i].end(), 0.);
          std::fill(occMultNTracksITSOnlyUnfm80[i].begin(), occMultNTracksITSOnlyUnfm80[i].end(), 0.);
          std::fill(occMultNTracksTPCOnlyUnfm80[i].begin(), occMultNTracksTPCOnlyUnfm80[i].end(), 0.);
          std::fill(occMultNTracksITSTPCUnfm80[i].begin(), occMultNTracksITSTPCUnfm80[i].end(), 0.);
          std::fill(occMultAllTracksTPCOnlyUnfm80[i].begin(), occMultAllTracksTPCOnlyUnfm80[i].end(), 0.);
        }
      }

      std::vector<int64_t> tfIDList;
      int nTrackITS = 0;
      int nTrackTPC = 0;
      int nTrackTRD = 0;
      int nTrackTOF = 0;
      int nTrackTPCA = 0;
      int nTrackTPCC = 0;
      int nTrackITSTPCA = 0;
      int nTrackITSTPCC = 0;

      ushort fNumContrib = 0;

      float fMultFV0A = -99999, fMultFV0C = -99999;
      float fMultFT0A = -99999, fMultFT0C = -99999;
      float fMultFDDA = -99999, fMultFDDC = -99999;

      int fNTrackITS = -9999;
      int fNTrackTPC = -9999;
      int fNTrackTRD = -9999;
      int fNTrackTOF = -9999;
      int fNTrackTPCA = -9999;
      int fNTrackTPCC = -9999;
      int fNTrackSize = -9999;
      int fNTrackITSTPC = -9999;
      int fNTrackITSTPCA = -9999;
      int fNTrackITSTPCC = -9999;

      decltype(&occPrimUnfm80[0]) tfOccPrimUnfm80 = nullptr;
      decltype(&occFV0AUnfm80[0]) tfOccFV0AUnfm80 = nullptr;
      decltype(&occFV0CUnfm80[0]) tfOccFV0CUnfm80 = nullptr;
      decltype(&occFT0AUnfm80[0]) tfOccFT0AUnfm80 = nullptr;
      decltype(&occFT0CUnfm80[0]) tfOccFT0CUnfm80 = nullptr;

      decltype(&occFDDAUnfm80[0]) tfOccFDDAUnfm80 = nullptr;
      decltype(&occFDDCUnfm80[0]) tfOccFDDCUnfm80 = nullptr;

      decltype(&occNTrackITSTPCUnfm80[0]) tfOccNTrackITSTPCUnfm80 = nullptr;

      decltype(&occNTrackITSUnfm80[0]) tfOccNTrackITSUnfm80 = nullptr;
      decltype(&occNTrackTPCUnfm80[0]) tfOccNTrackTPCUnfm80 = nullptr;
      decltype(&occNTrackTRDUnfm80[0]) tfOccNTrackTRDUnfm80 = nullptr;
      decltype(&occNTrackTOFUnfm80[0]) tfOccNTrackTOFUnfm80 = nullptr;
      decltype(&occNTrackSizeUnfm80[0]) tfOccNTrackSizeUnfm80 = nullptr;
      decltype(&occNTrackTPCAUnfm80[0]) tfOccNTrackTPCAUnfm80 = nullptr;
      decltype(&occNTrackTPCCUnfm80[0]) tfOccNTrackTPCCUnfm80 = nullptr;
      decltype(&occNTrackITSTPCAUnfm80[0]) tfOccNTrackITSTPCAUnfm80 = nullptr;
      decltype(&occNTrackITSTPCCUnfm80[0]) tfOccNTrackITSTPCCUnfm80 = nullptr;

      decltype(&occMultNTracksHasITSUnfm80[0]) tfOccMultNTracksHasITSUnfm80 = nullptr;
      decltype(&occMultNTracksHasTPCUnfm80[0]) tfOccMultNTracksHasTPCUnfm80 = nullptr;
      decltype(&occMultNTracksHasTOFUnfm80[0]) tfOccMultNTracksHasTOFUnfm80 = nullptr;
      decltype(&occMultNTracksHasTRDUnfm80[0]) tfOccMultNTracksHasTRDUnfm80 = nullptr;
      decltype(&occMultNTracksITSOnlyUnfm80[0]) tfOccMultNTracksITSOnlyUnfm80 = nullptr;
      decltype(&occMultNTracksTPCOnlyUnfm80[0]) tfOccMultNTracksTPCOnlyUnfm80 = nullptr;
      decltype(&occMultNTracksITSTPCUnfm80[0]) tfOccMultNTracksITSTPCUnfm80 = nullptr;
      decltype(&occMultAllTracksTPCOnlyUnfm80[0]) tfOccMultAllTracksTPCOnlyUnfm80 = nullptr;

      for (const auto& collision : collisions) {
        const auto& bc = collision.template bc_as<B>();
        getTimingInfo(bc, lastRun, nBCsPerTF, bcSOR, time, tfIdThis, bcInTF);

        recoEvent.fill(HIST("h_nBCinTF"), nBCsPerTF);
        recoEvent.fill(HIST("h_bcInTF"), bcInTF);

        if (nBCsPerTF > nBCinTF) {
          LOG(error) << "DEBUG :: FATAL ERROR :: nBCsPerTF > nBCinTF i.e " << nBCsPerTF << " > " << nBCinTF << " will cause crash in further process";
        }

        if constexpr (processMode == kProcessFullOccTableProducer || processMode == kProcessOnlyOccNtrackDet) {
          const uint64_t collIdx = collision.globalIndex();
          const auto tracksTablePerColl = tracks.sliceBy(tracksPerCollisionPreslice, collIdx);

          fNTrackSize = tracksTablePerColl.size();

          nTrackITS = 0;
          nTrackTPC = 0;
          nTrackTRD = 0;
          nTrackTOF = 0;
          nTrackTPCA = 0;
          nTrackTPCC = 0;
          nTrackITSTPCA = 0;
          nTrackITSTPCC = 0;
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
        }
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
        if constexpr (processMode == kProcessFullOccTableProducer || processMode == kProcessOnlyOccPrim || processMode == kProcessOnlyOccT0V0Prim || processMode == kProcessOnlyOccFDDT0V0Prim || processMode == kProcessOnlyOccNtrackDet || processMode == kProcessOnlyOccMultExtra) {
          tfOccPrimUnfm80 = &occPrimUnfm80[tfIDX];
        }
        if constexpr (processMode == kProcessFullOccTableProducer || processMode == kProcessOnlyOccT0V0Prim || processMode == kProcessOnlyOccFDDT0V0Prim) {
          tfOccFV0AUnfm80 = &occFV0AUnfm80[tfIDX];
          tfOccFV0CUnfm80 = &occFV0CUnfm80[tfIDX];
          tfOccFT0AUnfm80 = &occFT0AUnfm80[tfIDX];
          tfOccFT0CUnfm80 = &occFT0CUnfm80[tfIDX];
        }
        if constexpr (processMode == kProcessFullOccTableProducer || processMode == kProcessOnlyOccFDDT0V0Prim) {
          tfOccFDDAUnfm80 = &occFDDAUnfm80[tfIDX];
          tfOccFDDCUnfm80 = &occFDDCUnfm80[tfIDX];
        }
        if constexpr (processMode == kProcessFullOccTableProducer || processMode == kProcessOnlyOccNtrackDet) {
          tfOccNTrackITSUnfm80 = &occNTrackITSUnfm80[tfIDX];
          tfOccNTrackTPCUnfm80 = &occNTrackTPCUnfm80[tfIDX];
          tfOccNTrackTRDUnfm80 = &occNTrackTRDUnfm80[tfIDX];
          tfOccNTrackTOFUnfm80 = &occNTrackTOFUnfm80[tfIDX];
          tfOccNTrackSizeUnfm80 = &occNTrackSizeUnfm80[tfIDX];
          tfOccNTrackTPCAUnfm80 = &occNTrackTPCAUnfm80[tfIDX];
          tfOccNTrackTPCCUnfm80 = &occNTrackTPCCUnfm80[tfIDX];
          tfOccNTrackITSTPCAUnfm80 = &occNTrackITSTPCAUnfm80[tfIDX];
          tfOccNTrackITSTPCCUnfm80 = &occNTrackITSTPCCUnfm80[tfIDX];
        }
        if constexpr (processMode == kProcessFullOccTableProducer || processMode == kProcessOnlyOccNtrackDet || processMode == kProcessOnlyOccMultExtra) {
          tfOccNTrackITSTPCUnfm80 = &occNTrackITSTPCUnfm80[tfIDX];
        }
        if constexpr (processMode == kProcessFullOccTableProducer || processMode == kProcessOnlyOccMultExtra) {
          tfOccMultNTracksHasITSUnfm80 = &occMultNTracksHasITSUnfm80[tfIDX];
          tfOccMultNTracksHasTPCUnfm80 = &occMultNTracksHasTPCUnfm80[tfIDX];
          tfOccMultNTracksHasTOFUnfm80 = &occMultNTracksHasTOFUnfm80[tfIDX];
          tfOccMultNTracksHasTRDUnfm80 = &occMultNTracksHasTRDUnfm80[tfIDX];
          tfOccMultNTracksITSOnlyUnfm80 = &occMultNTracksITSOnlyUnfm80[tfIDX];
          tfOccMultNTracksTPCOnlyUnfm80 = &occMultNTracksTPCOnlyUnfm80[tfIDX];
          tfOccMultNTracksITSTPCUnfm80 = &occMultNTracksITSTPCUnfm80[tfIDX];
          tfOccMultAllTracksTPCOnlyUnfm80 = &occMultAllTracksTPCOnlyUnfm80[tfIDX];
        }

        // current collision bin in 80/160 bcGrouping.
        int bin80Zero = bcInTF / bcGrouping;
        // int bin160_0=bcInTF/160;

        // float fbin80Zero =float(bcInTF)/80;
        // float fbin160_0=float(bcInTF)/160;

        if constexpr (processMode == kProcessFullOccTableProducer || processMode == kProcessOnlyOccPrim || processMode == kProcessOnlyOccT0V0Prim || processMode == kProcessOnlyOccFDDT0V0Prim || processMode == kProcessOnlyOccNtrackDet || processMode == kProcessOnlyOccMultExtra) {
          fNumContrib = collision.numContrib(); // only aod::Collisions will be needed
        }
        if constexpr (processMode == kProcessFullOccTableProducer || processMode == kProcessOnlyOccT0V0Prim || processMode == kProcessOnlyOccFDDT0V0Prim) {
          fMultFV0A = collision.multFV0A();
          fMultFV0C = collision.multFV0C(); // o2::aod::Mults will be needed
          fMultFT0A = collision.multFT0A();
          fMultFT0C = collision.multFT0C();
        }
        if constexpr (processMode == kProcessFullOccTableProducer || processMode == kProcessOnlyOccFDDT0V0Prim) {
          fMultFDDA = collision.multFDDA();
          fMultFDDC = collision.multFDDC();
        }
        if constexpr (processMode == kProcessFullOccTableProducer || processMode == kProcessOnlyOccNtrackDet) {
          fNTrackITS = nTrackITS;
          fNTrackTPC = nTrackTPC;
          fNTrackTRD = nTrackTRD;
          fNTrackTOF = nTrackTOF;
          fNTrackTPCA = nTrackTPCA;
          fNTrackTPCC = nTrackTPCC;
          fNTrackITSTPC = collision.multAllTracksITSTPC();
          fNTrackITSTPCA = nTrackITSTPCA;
          fNTrackITSTPCC = nTrackITSTPCC;
        }
        // Processing for bcGrouping of 80 BCs
        for (int deltaBin = 0; deltaBin < nBCinDrift / bcGrouping; deltaBin++) {

          if constexpr (processMode == kProcessFullOccTableProducer || processMode == kProcessOnlyOccPrim || processMode == kProcessOnlyOccT0V0Prim || processMode == kProcessOnlyOccFDDT0V0Prim || processMode == kProcessOnlyOccNtrackDet || processMode == kProcessOnlyOccMultExtra) {
            (*tfOccPrimUnfm80)[(bin80Zero + deltaBin) % (nBCinTF / bcGrouping)] += fNumContrib * 1;
          }
          if constexpr (processMode == kProcessFullOccTableProducer || processMode == kProcessOnlyOccT0V0Prim || processMode == kProcessOnlyOccFDDT0V0Prim) {
            (*tfOccFV0AUnfm80)[(bin80Zero + deltaBin) % (nBCinTF / bcGrouping)] += fMultFV0A * 1;
            (*tfOccFV0CUnfm80)[(bin80Zero + deltaBin) % (nBCinTF / bcGrouping)] += fMultFV0C * 1;
            (*tfOccFT0AUnfm80)[(bin80Zero + deltaBin) % (nBCinTF / bcGrouping)] += fMultFT0A * 1;
            (*tfOccFT0CUnfm80)[(bin80Zero + deltaBin) % (nBCinTF / bcGrouping)] += fMultFT0C * 1;
          }
          if constexpr (processMode == kProcessFullOccTableProducer || processMode == kProcessOnlyOccFDDT0V0Prim) {
            (*tfOccFDDAUnfm80)[(bin80Zero + deltaBin) % (nBCinTF / bcGrouping)] += fMultFDDA * 1;
            (*tfOccFDDCUnfm80)[(bin80Zero + deltaBin) % (nBCinTF / bcGrouping)] += fMultFDDC * 1;
          }
          if constexpr (processMode == kProcessFullOccTableProducer || processMode == kProcessOnlyOccNtrackDet) {
            (*tfOccNTrackITSUnfm80)[(bin80Zero + deltaBin) % (nBCinTF / bcGrouping)] += fNTrackITS * 1;
            (*tfOccNTrackTPCUnfm80)[(bin80Zero + deltaBin) % (nBCinTF / bcGrouping)] += fNTrackTPC * 1;
            (*tfOccNTrackTRDUnfm80)[(bin80Zero + deltaBin) % (nBCinTF / bcGrouping)] += fNTrackTRD * 1;
            (*tfOccNTrackTOFUnfm80)[(bin80Zero + deltaBin) % (nBCinTF / bcGrouping)] += fNTrackTOF * 1;
            (*tfOccNTrackSizeUnfm80)[(bin80Zero + deltaBin) % (nBCinTF / bcGrouping)] += fNTrackSize * 1;
            (*tfOccNTrackTPCAUnfm80)[(bin80Zero + deltaBin) % (nBCinTF / bcGrouping)] += fNTrackTPCA * 1;
            (*tfOccNTrackTPCCUnfm80)[(bin80Zero + deltaBin) % (nBCinTF / bcGrouping)] += fNTrackTPCC * 1;
            (*tfOccNTrackITSTPCAUnfm80)[(bin80Zero + deltaBin) % (nBCinTF / bcGrouping)] += fNTrackITSTPCA * 1;
            (*tfOccNTrackITSTPCCUnfm80)[(bin80Zero + deltaBin) % (nBCinTF / bcGrouping)] += fNTrackITSTPCC * 1;
          }
          if constexpr (processMode == kProcessFullOccTableProducer || processMode == kProcessOnlyOccNtrackDet || processMode == kProcessOnlyOccMultExtra) {
            (*tfOccNTrackITSTPCUnfm80)[(bin80Zero + deltaBin) % (nBCinTF / bcGrouping)] += fNTrackITSTPC * 1;
          }

          if constexpr (processMode == kProcessFullOccTableProducer || processMode == kProcessOnlyOccMultExtra) {
            (*tfOccMultNTracksHasITSUnfm80)[(bin80Zero + deltaBin) % (nBCinTF / bcGrouping)] += collision.multNTracksHasITS() * 1;
            (*tfOccMultNTracksHasTPCUnfm80)[(bin80Zero + deltaBin) % (nBCinTF / bcGrouping)] += collision.multNTracksHasTPC() * 1;
            (*tfOccMultNTracksHasTOFUnfm80)[(bin80Zero + deltaBin) % (nBCinTF / bcGrouping)] += collision.multNTracksHasTOF() * 1;
            (*tfOccMultNTracksHasTRDUnfm80)[(bin80Zero + deltaBin) % (nBCinTF / bcGrouping)] += collision.multNTracksHasTRD() * 1;
            (*tfOccMultNTracksITSOnlyUnfm80)[(bin80Zero + deltaBin) % (nBCinTF / bcGrouping)] += collision.multNTracksITSOnly() * 1;
            (*tfOccMultNTracksTPCOnlyUnfm80)[(bin80Zero + deltaBin) % (nBCinTF / bcGrouping)] += collision.multNTracksTPCOnly() * 1;
            (*tfOccMultNTracksITSTPCUnfm80)[(bin80Zero + deltaBin) % (nBCinTF / bcGrouping)] += collision.multNTracksITSTPC() * 1;
            (*tfOccMultAllTracksTPCOnlyUnfm80)[(bin80Zero + deltaBin) % (nBCinTF / bcGrouping)] += collision.multAllTracksTPCOnly() * 1;
          }
        }
      }
      // collision Loop is over

      occupancyQA.fill(HIST("h_TF_in_DataFrame"), tfCounted);

      std::vector<int64_t> sortedTfIDList = tfIDList;
      std::sort(sortedTfIDList.begin(), sortedTfIDList.end());
      auto last = std::unique(sortedTfIDList.begin(), sortedTfIDList.end());
      sortedTfIDList.erase(last, sortedTfIDList.end());

      if (tfCounted != sortedTfIDList.size()) {
        LOG(error) << "DEBUG :: Number mismatch for tf counted and filled :: " << tfCounted << " != " << sortedTfIDList.size();
      }

      int totalBCcountSize = 0;
      for (int i = 0; i < occVecArraySize; i++) {
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

        genOccsBCsList(tfList[i], bcTFMap[i]);

        auto& vecOccPrimUnfm80 = occPrimUnfm80[i];
        float meanOccPrimUnfm80 = TMath::Mean(vecOccPrimUnfm80.size(), vecOccPrimUnfm80.data());
        normalizeVector(vecOccPrimUnfm80, meanOccPrimUnfm80 / meanOccPrimUnfm80);

        if constexpr (processMode == kProcessFullOccTableProducer || processMode == kProcessOnlyOccPrim || processMode == kProcessOnlyOccT0V0Prim || processMode == kProcessOnlyOccFDDT0V0Prim || processMode == kProcessOnlyOccNtrackDet || processMode == kProcessOnlyOccMultExtra) {
          if constexpr (tableMode == fillOccTable) {
            genOccsPrim(vecOccPrimUnfm80);
          }
          if constexpr (meanTableMode == fillMeanOccTable) {
            genOccsMeanPrim(meanOccPrimUnfm80);
          }
        }

        if constexpr (processMode == kProcessFullOccTableProducer || processMode == kProcessOnlyOccT0V0Prim || processMode == kProcessOnlyOccFDDT0V0Prim) {
          auto& vecOccFV0AUnfm80 = occFV0AUnfm80[i];
          auto& vecOccFV0CUnfm80 = occFV0CUnfm80[i];
          auto& vecOccFT0AUnfm80 = occFT0AUnfm80[i];
          auto& vecOccFT0CUnfm80 = occFT0CUnfm80[i];

          float meanOccFV0AUnfm80 = TMath::Mean(vecOccFV0AUnfm80.size(), vecOccFV0AUnfm80.data());
          float meanOccFV0CUnfm80 = TMath::Mean(vecOccFV0CUnfm80.size(), vecOccFV0CUnfm80.data());
          float meanOccFT0AUnfm80 = TMath::Mean(vecOccFT0AUnfm80.size(), vecOccFT0AUnfm80.data());
          float meanOccFT0CUnfm80 = TMath::Mean(vecOccFT0CUnfm80.size(), vecOccFT0CUnfm80.data());

          // Normalise the original vectors
          normalizeVector(vecOccFV0AUnfm80, meanOccPrimUnfm80 / meanOccFV0AUnfm80);
          normalizeVector(vecOccFV0CUnfm80, meanOccPrimUnfm80 / meanOccFV0CUnfm80);
          normalizeVector(vecOccFT0AUnfm80, meanOccPrimUnfm80 / meanOccFT0AUnfm80);
          normalizeVector(vecOccFT0CUnfm80, meanOccPrimUnfm80 / meanOccFT0CUnfm80);

          // Find Robust estimators
          // T0A, T0C, V0A, Prim
          if constexpr (robustTableMode == fillOccRobustTable || meanRobustTableMode == fillOccMeanRobustTable) {
            getMedianOccVect(vecRobustOccT0V0PrimUnfm80, vecRobustOccT0V0PrimUnfm80medianPosVec,
                             vecOccPrimUnfm80, vecOccFV0AUnfm80, vecOccFT0AUnfm80, vecOccFT0CUnfm80);
            for (const auto& vec : vecRobustOccT0V0PrimUnfm80medianPosVec) {
              recoEvent.fill(HIST("h_RO_T0V0PrimUnfm80"), vec[0]);
              recoEvent.fill(HIST("h_RO_T0V0PrimUnfm80"), vec[1]);
            }
          }
          if constexpr (processMode == kProcessFullOccTableProducer || processMode == kProcessOnlyOccT0V0Prim || processMode == kProcessOnlyOccFDDT0V0Prim) {
            if constexpr (tableMode == fillOccTable) {
              genOccsT0V0(vecOccFV0AUnfm80, vecOccFV0CUnfm80, vecOccFT0AUnfm80, vecOccFT0CUnfm80);
            }
            if constexpr (meanTableMode == fillMeanOccTable) {
              genOccsMeanT0V0(meanOccFV0AUnfm80, meanOccFV0CUnfm80, meanOccFT0AUnfm80, meanOccFT0CUnfm80);
            }
            if constexpr (robustTableMode == fillOccRobustTable) {
              genORT0V0Prim(vecRobustOccT0V0PrimUnfm80);
            }
            if constexpr (meanRobustTableMode == fillOccMeanRobustTable) {
              genOccsMeanRobustT0V0Prim(TMath::Mean(vecRobustOccT0V0PrimUnfm80.size(), vecRobustOccT0V0PrimUnfm80.data()));
            }

            if constexpr (processMode == kProcessFullOccTableProducer || processMode == kProcessOnlyOccFDDT0V0Prim) {
              auto& vecOccFDDAUnfm80 = occFDDAUnfm80[i];
              auto& vecOccFDDCUnfm80 = occFDDCUnfm80[i];
              float meanOccFDDAUnfm80 = TMath::Mean(vecOccFDDAUnfm80.size(), vecOccFDDAUnfm80.data());
              float meanOccFDDCUnfm80 = TMath::Mean(vecOccFDDCUnfm80.size(), vecOccFDDCUnfm80.data());
              normalizeVector(vecOccFDDAUnfm80, meanOccPrimUnfm80 / meanOccFDDAUnfm80);
              normalizeVector(vecOccFDDCUnfm80, meanOccPrimUnfm80 / meanOccFDDCUnfm80);

              // T0A, T0C, V0A, FDD, Prim
              if constexpr (robustTableMode == fillOccRobustTable || meanRobustTableMode == fillOccMeanRobustTable) {
                getMedianOccVect(vecRobustOccFDDT0V0PrimUnfm80, vecRobustOccFDDT0V0PrimUnfm80medianPosVec,
                                 vecOccPrimUnfm80, vecOccFV0AUnfm80, vecOccFT0AUnfm80, vecOccFT0CUnfm80, vecOccFDDAUnfm80, vecOccFDDCUnfm80);
                for (const auto& vec : vecRobustOccFDDT0V0PrimUnfm80medianPosVec) {
                  recoEvent.fill(HIST("h_RO_FDDT0V0PrimUnfm80"), vec[0]);
                  recoEvent.fill(HIST("h_RO_FDDT0V0PrimUnfm80"), vec[1]);
                }

                if constexpr (tableMode == fillOccTable) {
                  genOccsFDD(vecOccFDDAUnfm80, vecOccFDDCUnfm80);
                }
                if constexpr (meanTableMode == fillMeanOccTable) {
                  genOccsMeanFDD(meanOccFDDAUnfm80, meanOccFDDCUnfm80);
                }
                if constexpr (robustTableMode == fillOccRobustTable) {
                  genORFDDT0V0Prim(vecRobustOccFDDT0V0PrimUnfm80);
                }
                if constexpr (meanRobustTableMode == fillOccMeanRobustTable) {
                  genOccsMeanRobustFDDT0V0Prim(TMath::Mean(vecRobustOccFDDT0V0PrimUnfm80.size(), vecRobustOccFDDT0V0PrimUnfm80.data()));
                }
              }
            } // Block for FDDT0V0Prim
          } // For T0V0Prim only and FDDT0V0Prim
        } // Detector Occupancy block

        if constexpr (processMode == kProcessFullOccTableProducer || processMode == kProcessOnlyOccNtrackDet) {
          // NTrackDet
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

          if constexpr (robustTableMode == fillOccRobustTable || meanRobustTableMode == fillOccMeanRobustTable) {
            getMedianOccVect(vecRobustOccNtrackDetUnfm80, vecRobustOccNtrackDetUnfm80medianPosVec,
                             vecOccNTrackITSUnfm80, vecOccNTrackTPCUnfm80, vecOccNTrackTRDUnfm80, vecOccNTrackTOFUnfm80);
            for (const auto& vec : vecRobustOccNtrackDetUnfm80medianPosVec) {
              recoEvent.fill(HIST("h_RO_NtrackDetUnfm80"), vec[0]);
              recoEvent.fill(HIST("h_RO_NtrackDetUnfm80"), vec[1]);
            }
          }

          if constexpr (tableMode == fillOccTable) {
            genOccsNTrackDet(vecOccNTrackITSUnfm80, vecOccNTrackTPCUnfm80,
                             vecOccNTrackTRDUnfm80, vecOccNTrackTOFUnfm80,
                             vecOccNTrackSizeUnfm80, vecOccNTrackTPCAUnfm80,
                             vecOccNTrackTPCCUnfm80, vecOccNTrackITSTPCUnfm80,
                             vecOccNTrackITSTPCAUnfm80, vecOccNTrackITSTPCCUnfm80);
          }
          if constexpr (meanTableMode == fillMeanOccTable) {
            genOccsMeanNTrkDet(meanOccNTrackITSUnfm80, meanOccNTrackTPCUnfm80,
                               meanOccNTrackTRDUnfm80, meanOccNTrackTOFUnfm80,
                               meanOccNTrackSizeUnfm80, meanOccNTrackTPCAUnfm80,
                               meanOccNTrackTPCCUnfm80, meanOccNTrackITSTPCUnfm80,
                               meanOccNTrackITSTPCAUnfm80, meanOccNTrackITSTPCCUnfm80);
          }
          if constexpr (robustTableMode == fillOccRobustTable) {
            genORNtrackDet(vecRobustOccNtrackDetUnfm80);
          }
          if constexpr (meanRobustTableMode == fillOccMeanRobustTable) {
            genOccsMeanRobustNtrackDet(TMath::Mean(vecRobustOccNtrackDetUnfm80.size(), vecRobustOccNtrackDetUnfm80.data()));
          }
        }

        if constexpr (processMode == kProcessFullOccTableProducer || processMode == kProcessOnlyOccMultExtra) {
          auto& vecOccNTrackITSTPCUnfm80 = occNTrackITSTPCUnfm80[i];
          float meanOccNTrackITSTPCUnfm80 = TMath::Mean(vecOccNTrackITSTPCUnfm80.size(), vecOccNTrackITSTPCUnfm80.data());
          normalizeVector(vecOccNTrackITSTPCUnfm80, meanOccPrimUnfm80 / meanOccNTrackITSTPCUnfm80);

          auto& vecOccMultNTracksHasITSUnfm80 = occMultNTracksHasITSUnfm80[i];
          auto& vecOccMultNTracksHasTPCUnfm80 = occMultNTracksHasTPCUnfm80[i];
          auto& vecOccMultNTracksHasTOFUnfm80 = occMultNTracksHasTOFUnfm80[i];
          auto& vecOccMultNTracksHasTRDUnfm80 = occMultNTracksHasTRDUnfm80[i];
          auto& vecOccMultNTracksITSOnlyUnfm80 = occMultNTracksITSOnlyUnfm80[i];
          auto& vecOccMultNTracksTPCOnlyUnfm80 = occMultNTracksTPCOnlyUnfm80[i];
          auto& vecOccMultNTracksITSTPCUnfm80 = occMultNTracksITSTPCUnfm80[i];
          auto& vecOccMultAllTracksTPCOnlyUnfm80 = occMultAllTracksTPCOnlyUnfm80[i];
          float meanOccMultNTracksHasITSUnfm80 = TMath::Mean(vecOccMultNTracksHasITSUnfm80.size(), vecOccMultNTracksHasITSUnfm80.data());
          float meanOccMultNTracksHasTPCUnfm80 = TMath::Mean(vecOccMultNTracksHasTPCUnfm80.size(), vecOccMultNTracksHasTPCUnfm80.data());
          float meanOccMultNTracksHasTOFUnfm80 = TMath::Mean(vecOccMultNTracksHasTOFUnfm80.size(), vecOccMultNTracksHasTOFUnfm80.data());
          float meanOccMultNTracksHasTRDUnfm80 = TMath::Mean(vecOccMultNTracksHasTRDUnfm80.size(), vecOccMultNTracksHasTRDUnfm80.data());
          float meanOccMultNTracksITSOnlyUnfm80 = TMath::Mean(vecOccMultNTracksITSOnlyUnfm80.size(), vecOccMultNTracksITSOnlyUnfm80.data());
          float meanOccMultNTracksTPCOnlyUnfm80 = TMath::Mean(vecOccMultNTracksTPCOnlyUnfm80.size(), vecOccMultNTracksTPCOnlyUnfm80.data());
          float meanOccMultNTracksITSTPCUnfm80 = TMath::Mean(vecOccMultNTracksITSTPCUnfm80.size(), vecOccMultNTracksITSTPCUnfm80.data());
          float meanOccMultAllTracksTPCOnlyUnfm80 = TMath::Mean(vecOccMultAllTracksTPCOnlyUnfm80.size(), vecOccMultAllTracksTPCOnlyUnfm80.data());
          // multExtraTable
          normalizeVector(vecOccMultNTracksHasITSUnfm80, meanOccPrimUnfm80 / meanOccMultNTracksHasITSUnfm80);
          normalizeVector(vecOccMultNTracksHasTPCUnfm80, meanOccPrimUnfm80 / meanOccMultNTracksHasTPCUnfm80);
          normalizeVector(vecOccMultNTracksHasTOFUnfm80, meanOccPrimUnfm80 / meanOccMultNTracksHasTOFUnfm80);
          normalizeVector(vecOccMultNTracksHasTRDUnfm80, meanOccPrimUnfm80 / meanOccMultNTracksHasTRDUnfm80);
          normalizeVector(vecOccMultNTracksITSOnlyUnfm80, meanOccPrimUnfm80 / meanOccMultNTracksITSOnlyUnfm80);
          normalizeVector(vecOccMultNTracksTPCOnlyUnfm80, meanOccPrimUnfm80 / meanOccMultNTracksTPCOnlyUnfm80);
          normalizeVector(vecOccMultNTracksITSTPCUnfm80, meanOccPrimUnfm80 / meanOccMultNTracksITSTPCUnfm80);
          normalizeVector(vecOccMultAllTracksTPCOnlyUnfm80, meanOccPrimUnfm80 / meanOccMultAllTracksTPCOnlyUnfm80);

          if constexpr (robustTableMode == fillOccRobustTable || meanRobustTableMode == fillOccMeanRobustTable) {
            getMedianOccVect(vecRobustOccmultTableUnfm80, vecRobustOccmultTableUnfm80medianPosVec,
                             vecOccPrimUnfm80, vecOccMultNTracksHasITSUnfm80, vecOccMultNTracksHasTPCUnfm80,
                             vecOccMultNTracksHasTOFUnfm80, vecOccMultNTracksHasTRDUnfm80, vecOccMultNTracksITSOnlyUnfm80,
                             vecOccMultNTracksTPCOnlyUnfm80, vecOccMultNTracksITSTPCUnfm80, vecOccMultAllTracksTPCOnlyUnfm80,
                             vecOccNTrackITSTPCUnfm80);
            for (const auto& vec : vecRobustOccmultTableUnfm80medianPosVec) {
              recoEvent.fill(HIST("h_RO_multTableUnfm80"), vec[0]);
              recoEvent.fill(HIST("h_RO_multTableUnfm80"), vec[1]);
            }
          }

          if constexpr (tableMode == fillOccTable) {
            genOccsMultExtra(vecOccMultNTracksHasITSUnfm80, vecOccMultNTracksHasTPCUnfm80,
                             vecOccMultNTracksHasTOFUnfm80, vecOccMultNTracksHasTRDUnfm80,
                             vecOccMultNTracksITSOnlyUnfm80, vecOccMultNTracksTPCOnlyUnfm80,
                             vecOccMultNTracksITSTPCUnfm80, vecOccMultAllTracksTPCOnlyUnfm80);
          }
          if constexpr (meanTableMode == fillMeanOccTable) {
            genOccsMnMultExtra(meanOccMultNTracksHasITSUnfm80, meanOccMultNTracksHasTPCUnfm80,
                               meanOccMultNTracksHasTOFUnfm80, meanOccMultNTracksHasTRDUnfm80,
                               meanOccMultNTracksITSOnlyUnfm80, meanOccMultNTracksTPCOnlyUnfm80,
                               meanOccMultNTracksITSTPCUnfm80, meanOccMultAllTracksTPCOnlyUnfm80);
          }
          if constexpr (robustTableMode == fillOccRobustTable) {
            genORMultExtra(vecRobustOccmultTableUnfm80);
          }
          if constexpr (meanRobustTableMode == fillOccMeanRobustTable) {
            genOccsMeanRobustMultExtraTable(TMath::Mean(vecRobustOccmultTableUnfm80.size(), vecRobustOccmultTableUnfm80.data()));
          }
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

        genBCTFinfoTable(tfIdThis, bcInTF);
        genOccIndexTable(bc.globalIndex(), occIDX); // BCId, OccId
      }
    } // else block for constexpr
  }

  void checkAllProcessFunctionStatus(std::vector<bool> const& processStatusVector, bool& singleProcessOn)
  {
    int nProcessOn = 0;
    const uint size = processStatusVector.size();
    for (uint i = 0; i < size; i++) {
      if (processStatusVector[i]) {
        nProcessOn++;
      }
    }

    if (nProcessOn > 1) {
      singleProcessOn = false;
      std::ostringstream warningLine;
      warningLine << "DEBUG :: More than one track-mean-occ-table-producer process function is on :: ";
      for (uint i = 0; i < size; i++) {
        if (processStatusVector[i]) {
          warningLine << std::string(ProcessNames[processStatusVector[i]]) << " == true :: ";
        }
      }
      LOG(error) << warningLine.str();
    } // check nProcess
  }

  //________________________________________End of Exection Function_________________________________________________________________________

  void processOnlyBCTFinfoTable(o2::aod::BCsWithTimestamps const& BCs)
  {
    occupancyQA.fill(HIST("h_DFcount_Lvl0"), 0.5);
    processStatus[kProcessOnlyBCTFinfoTable] = true;
    bool singleProcessOn = true;
    checkAllProcessFunctionStatus(processStatus, singleProcessOn);
    if (!singleProcessOn) {
      return;
    }

    for (const auto& BC : BCs) {
      getTimingInfo(BC, lastRun, nBCsPerTF, bcSOR, time, tfIdThis, bcInTF);
      genBCTFinfoTable(tfIdThis, bcInTF);
    }
    occupancyQA.fill(HIST("h_DFcount_Lvl1"), 0.5);
  }
  PROCESS_SWITCH(OccupancyTableProducer, processOnlyBCTFinfoTable, "processOnlyBCTFinfoTable", true);

  // // Process the Data
  void processOnlyOccPrimUnfm(o2::aod::BCsWithTimestamps const& BCs, aod::Collisions const& collisions)
  {
    occupancyQA.fill(HIST("h_DFcount_Lvl0"), 0.5);
    if (!buildOnlyOccsPrim) {
      LOG(error) << " DEBUG :: ERROR ERROR ERROR :: buildOnlyOccsPrim == false";
    }
    processStatus[kProcessOnlyOccPrim] = true;
    bool singleProcessOn = true;
    checkAllProcessFunctionStatus(processStatus, singleProcessOn);
    if (!singleProcessOn) {
      return;
    }

    bool collisionsSizeIsZero = false;
    executeCollisionCheckAndBCprocessing(BCs, collisions, collisionsSizeIsZero);
    if (collisionsSizeIsZero) {
      return;
    }
    executeOccProducerProcessing<kProcessOnlyOccPrim, checkTableMode, checkTableMode, checkTableMode, checkTableMode>(BCs, collisions, nullptr);
    occupancyQA.fill(HIST("h_DFcount_Lvl1"), 0.5);
  }
  PROCESS_SWITCH(OccupancyTableProducer, processOnlyOccPrimUnfm, "processOnlyOccPrimUnfm", false);

  void processOnlyOccT0V0PrimUnfm(o2::aod::BCsWithTimestamps const& BCs, soa::Join<aod::Collisions, aod::Mults> const& collisions)
  {
    occupancyQA.fill(HIST("h_DFcount_Lvl0"), 0.5);
    if (!buildOnlyOccsT0V0Prim) {
      LOG(error) << " DEBUG :: ERROR ERROR ERROR :: buildOnlyOccsT0V0Prim == false";
    }
    processStatus[kProcessOnlyOccT0V0Prim] = true;
    bool singleProcessOn = true;
    checkAllProcessFunctionStatus(processStatus, singleProcessOn);
    if (!singleProcessOn) {
      return;
    }

    bool collisionsSizeIsZero = false;
    executeCollisionCheckAndBCprocessing(BCs, collisions, collisionsSizeIsZero);
    if (collisionsSizeIsZero) {
      return;
    }
    executeOccProducerProcessing<kProcessOnlyOccT0V0Prim, checkTableMode, checkTableMode, checkTableMode, checkTableMode>(BCs, collisions, nullptr);
    occupancyQA.fill(HIST("h_DFcount_Lvl1"), 0.5);
  }
  PROCESS_SWITCH(OccupancyTableProducer, processOnlyOccT0V0PrimUnfm, "processOnlyOccT0V0PrimUnfm", false);

  void processOnlyOccFDDT0V0PrimUnfm(o2::aod::BCsWithTimestamps const& BCs, soa::Join<aod::Collisions, aod::Mults> const& collisions)
  {
    occupancyQA.fill(HIST("h_DFcount_Lvl0"), 0.5);
    if (!buildOnlyOccsFDDT0V0Prim) {
      LOG(error) << " DEBUG :: ERROR ERROR ERROR :: buildOnlyOccsFDDT0V0Prim == false";
    }
    processStatus[kProcessOnlyOccFDDT0V0Prim] = true;
    bool singleProcessOn = true;
    checkAllProcessFunctionStatus(processStatus, singleProcessOn);
    if (!singleProcessOn) {
      return;
    }

    bool collisionsSizeIsZero = false;
    executeCollisionCheckAndBCprocessing(BCs, collisions, collisionsSizeIsZero);
    if (collisionsSizeIsZero) {
      return;
    }
    executeOccProducerProcessing<kProcessOnlyOccFDDT0V0Prim, checkTableMode, checkTableMode, checkTableMode, checkTableMode>(BCs, collisions, nullptr);
    occupancyQA.fill(HIST("h_DFcount_Lvl1"), 0.5);
  }
  PROCESS_SWITCH(OccupancyTableProducer, processOnlyOccFDDT0V0PrimUnfm, "processOnlyOccFDDT0V0PrimUnfm", false);

  void processOnlyOccNtrackDet(o2::aod::BCsWithTimestamps const& BCs, soa::Join<aod::Collisions, aod::MultsExtra> const& collisions, soa::Join<aod::Tracks, o2::aod::TracksExtra> const& tracks)
  {
    occupancyQA.fill(HIST("h_DFcount_Lvl0"), 0.5);
    if (!buildOnlyOccsNtrackDet) {
      LOG(error) << " DEBUG :: ERROR ERROR ERROR :: buildOnlyOccsNtrackDet == false";
    }
    processStatus[kProcessOnlyOccNtrackDet] = true;
    bool singleProcessOn = true;
    checkAllProcessFunctionStatus(processStatus, singleProcessOn);
    if (!singleProcessOn) {
      return;
    }

    bool collisionsSizeIsZero = false;
    executeCollisionCheckAndBCprocessing(BCs, collisions, collisionsSizeIsZero);
    if (collisionsSizeIsZero) {
      return;
    }
    executeOccProducerProcessing<kProcessOnlyOccNtrackDet, checkTableMode, checkTableMode, checkTableMode, checkTableMode>(BCs, collisions, tracks);
    occupancyQA.fill(HIST("h_DFcount_Lvl1"), 0.5);
  }
  PROCESS_SWITCH(OccupancyTableProducer, processOnlyOccNtrackDet, "processOnlyOccNtrackDet", false);

  void processOnlyOccMultExtra(o2::aod::BCsWithTimestamps const& BCs, soa::Join<aod::Collisions, aod::MultsExtra> const& collisions)
  {
    occupancyQA.fill(HIST("h_DFcount_Lvl0"), 0.5);
    if (!buildOnlyOccsMultExtra) {
      LOG(error) << " DEBUG :: ERROR ERROR ERROR :: buildOnlyOccsMultExtra == false";
    }
    processStatus[kProcessOnlyOccMultExtra] = true;
    bool singleProcessOn = true;
    checkAllProcessFunctionStatus(processStatus, singleProcessOn);
    if (!singleProcessOn) {
      return;
    }

    bool collisionsSizeIsZero = false;
    executeCollisionCheckAndBCprocessing(BCs, collisions, collisionsSizeIsZero);
    if (collisionsSizeIsZero) {
      return;
    }
    executeOccProducerProcessing<kProcessOnlyOccMultExtra, checkTableMode, checkTableMode, checkTableMode, checkTableMode>(BCs, collisions, nullptr);
    occupancyQA.fill(HIST("h_DFcount_Lvl1"), 0.5);
  }
  PROCESS_SWITCH(OccupancyTableProducer, processOnlyOccMultExtra, "processOnlyOccMultExtra", false);

  void processFullOccTableProduer(o2::aod::BCsWithTimestamps const& BCs, soa::Join<aod::Collisions, aod::Mults, aod::MultsExtra> const& collisions, soa::Join<aod::Tracks, o2::aod::TracksExtra> const& tracks)
  {
    occupancyQA.fill(HIST("h_DFcount_Lvl0"), 0.5);
    if (!buildFullOccTableProducer) {
      LOG(error) << " DEBUG :: ERROR ERROR ERROR :: buildFullOccTableProducer == false";
    }
    processStatus[kProcessFullOccTableProducer] = true;
    bool singleProcessOn = true;
    checkAllProcessFunctionStatus(processStatus, singleProcessOn);
    if (!singleProcessOn) {
      return;
    }

    bool collisionsSizeIsZero = false;
    executeCollisionCheckAndBCprocessing(BCs, collisions, collisionsSizeIsZero);
    if (collisionsSizeIsZero) {
      return;
    }
    executeOccProducerProcessing<kProcessFullOccTableProducer, checkTableMode, checkTableMode, checkTableMode, checkTableMode>(BCs, collisions, tracks);
    occupancyQA.fill(HIST("h_DFcount_Lvl1"), 0.5);
  } // Process function ends
  PROCESS_SWITCH(OccupancyTableProducer, processFullOccTableProduer, "processFullOccTableProduer", false);
};

struct TrackMeanOccTableProducer {

  // //declare production of tables
  Produces<aod::TmoTrackIds> genTmoTrackId;
  Produces<aod::TmoToTrackQA> genTmoToTrackQA;
  Produces<aod::TrackQAToTmo> genTrackQAToTmo;

  Produces<aod::TmoPrim> genTmoPrim;
  Produces<aod::TmoT0V0> genTmoT0V0;
  Produces<aod::TmoFDD> genTmoFDD;
  Produces<aod::TmoNTrackDet> genTmoNTrackDet;
  Produces<aod::TmoMultExtra> genTmoMultExtra;
  Produces<aod::TmoRT0V0Prim> genTmoRT0V0Prim;
  Produces<aod::TmoRFDDT0V0Prim> genTmoRFDDT0V0Prim;
  Produces<aod::TmoRNtrackDet> genTmoRNtrackDet;
  Produces<aod::TmoRMultExtra> genTmoRMultExtra;

  Produces<aod::TwmoPrim> genTwmoPrim;
  Produces<aod::TwmoT0V0> genTwmoT0V0;
  Produces<aod::TwmoFDD> genTwmoFDD;
  Produces<aod::TwmoNTrackDet> genTwmoNTrackDet;
  Produces<aod::TwmoMultExtra> genTwmoMultExtra;
  Produces<aod::TwmoRT0V0Prim> genTwmoRT0V0Prim;
  Produces<aod::TwmoRFDDT0V0Pri> genTwmoRFDDT0V0Pri;
  Produces<aod::TwmoRNtrackDet> genTwmoRNtrackDet;
  Produces<aod::TwmoRMultExtra> genTwmoRMultExtra;

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  HistogramRegistry occupancyQA{"occupancyQA", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // Configurables
  Configurable<int> customOrbitOffset{"customOrbitOffset", 0, "customOrbitOffset for MC"};
  Configurable<int> bcGrouping{"bcGrouping", 80, "bcGrouping of BCs"};
  Configurable<int> nBCinTF{"nBCinTF", 114048, "nBCinTF"};

  Configurable<int> cfgNOrbitsPerTF0RunValue{"cfgNOrbitsPerTF0RunValue", 534133, "cfgNOrbitsPerTF0RunValue"};
  Configurable<int> cfgNOrbitsPerTF1TrueValue{"cfgNOrbitsPerTF1TrueValue", 128, "cfgNOrbitsPerTF1TrueValue"};
  Configurable<int> cfgNOrbitsPerTF2FalseValue{"cfgNOrbitsPerTF2FalseValue", 32, "ccfgNOrbitsPerTF2FalseValue"};

  Configurable<bool> buildOnlyOccsPrim{"buildOnlyOccsPrim", true, "builder of table OccsPrim"};
  Configurable<bool> buildOnlyOccsT0V0{"buildOnlyOccsT0V0", true, "builder of table OccsT0V0Prim"};
  Configurable<bool> buildOnlyOccsFDD{"buildOnlyOccsFDD", true, "builder of table OccsFDDT0V0Prim"};
  Configurable<bool> buildOnlyOccsNtrackDet{"buildOnlyOccsNtrackDet", true, "builder of table OccsNtrackDet"};
  Configurable<bool> buildOnlyOccsMultExtra{"buildOnlyOccsMultExtra", true, "builder of table OccsMultExtra"};

  Configurable<bool> buildOnlyOccsRobustT0V0Prim{"buildOnlyOccsRobustT0V0Prim", true, "build  buildOnlyOccsRobustT0V0Prim"};
  Configurable<bool> buildOnlyOccsRobustFDDT0V0Prim{"buildOnlyOccsRobustFDDT0V0Prim", true, "build  buildOnlyOccsRobustFDDT0V0Prim"};
  Configurable<bool> buildOnlyOccsRobustNtrackDet{"buildOnlyOccsRobustNtrackDet", true, "build  buildOnlyOccsRobustNtrackDet"};
  Configurable<bool> buildOnlyOccsRobustMultExtra{"buildOnlyOccsRobustMultExtra", true, "build  buildOnlyOccsRobustMultExtra"};

  Configurable<bool> buildFullOccTableProducer{"buildFullOccTableProducer", true, "builder of all Occupancy Tables"};

  Configurable<bool> buildFlag00MeanTable{"buildFlag00MeanTable", true, "build Flag00MeanTable"};
  Configurable<bool> buildFlag01WeightMeanTable{"buildFlag01WeightMeanTable", true, "build Flag01WeightMeanTable"};

  Configurable<bool> fillQA1{"fillQA1", true, "fill QA LOG Ratios"};
  Configurable<bool> fillQA2{"fillQA2", true, "fill QA condition dependent QAs"};

  Configurable<bool> buildPointerTrackQAToTMOTable{"buildPointerTrackQAToTMOTable", true, "buildPointerTrackQAToTMOTable"};
  Configurable<bool> buildPointerTMOToTrackQATable{"buildPointerTMOToTrackQATable", true, "buildPointerTMOToTrackQATable"};

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

  std::vector<bool> processStatus;
  std::vector<bool> processInThisBlock;
  void init(InitContext const&)
  {
    // CCDB related part to be added later
    processStatus.resize(11);
    processInThisBlock.resize(11);

    for (uint i = 0; i < processStatus.size(); i++) {
      processStatus[i] = false;
      processInThisBlock[i] = false;
    }

    if (buildFullOccTableProducer || buildOnlyOccsPrim) {
      occPrimUnfm80.resize(nBCinTF / bcGrouping);
    }
    if (buildFullOccTableProducer || buildOnlyOccsT0V0) {
      occFV0AUnfm80.resize(nBCinTF / bcGrouping);
      occFV0CUnfm80.resize(nBCinTF / bcGrouping);
      occFT0AUnfm80.resize(nBCinTF / bcGrouping);
      occFT0CUnfm80.resize(nBCinTF / bcGrouping);
    }
    if (buildFullOccTableProducer || buildOnlyOccsFDD) {
      occFDDAUnfm80.resize(nBCinTF / bcGrouping);
      occFDDCUnfm80.resize(nBCinTF / bcGrouping);
    }
    if (buildFullOccTableProducer || buildOnlyOccsNtrackDet) {
      occNTrackITSUnfm80.resize(nBCinTF / bcGrouping);
      occNTrackTPCUnfm80.resize(nBCinTF / bcGrouping);
      occNTrackTRDUnfm80.resize(nBCinTF / bcGrouping);
      occNTrackTOFUnfm80.resize(nBCinTF / bcGrouping);
      occNTrackSizeUnfm80.resize(nBCinTF / bcGrouping);
      occNTrackTPCAUnfm80.resize(nBCinTF / bcGrouping);
      occNTrackTPCCUnfm80.resize(nBCinTF / bcGrouping);
      occNTrackITSTPCUnfm80.resize(nBCinTF / bcGrouping);
      occNTrackITSTPCAUnfm80.resize(nBCinTF / bcGrouping);
      occNTrackITSTPCCUnfm80.resize(nBCinTF / bcGrouping);
    }

    if (buildFullOccTableProducer || buildOnlyOccsMultExtra) {
      occMultNTracksHasITSUnfm80.resize(nBCinTF / bcGrouping);
      occMultNTracksHasTPCUnfm80.resize(nBCinTF / bcGrouping);
      occMultNTracksHasTOFUnfm80.resize(nBCinTF / bcGrouping);
      occMultNTracksHasTRDUnfm80.resize(nBCinTF / bcGrouping);
      occMultNTracksITSOnlyUnfm80.resize(nBCinTF / bcGrouping);
      occMultNTracksTPCOnlyUnfm80.resize(nBCinTF / bcGrouping);
      occMultNTracksITSTPCUnfm80.resize(nBCinTF / bcGrouping);
      occMultAllTracksTPCOnlyUnfm80.resize(nBCinTF / bcGrouping);
    }

    if (buildFullOccTableProducer || buildOnlyOccsRobustT0V0Prim || fillQA1 || fillQA2) {
      occRobustT0V0PrimUnfm80.resize(nBCinTF / bcGrouping);
    }
    if (buildFullOccTableProducer || buildOnlyOccsRobustFDDT0V0Prim) {
      occRobustFDDT0V0PrimUnfm80.resize(nBCinTF / bcGrouping);
    }
    if (buildFullOccTableProducer || buildOnlyOccsRobustNtrackDet) {
      occRobustNtrackDetUnfm80.resize(nBCinTF / bcGrouping);
    }
    if (buildFullOccTableProducer || buildOnlyOccsRobustMultExtra) {
      occRobustMultTableUnfm80.resize(nBCinTF / bcGrouping);
    }

    const AxisSpec axisQA1 = {500, 0, 50000};
    const AxisSpec axisQA2 = {200, -2, 2};
    const AxisSpec axisQA3 = {200, -20, 20};

    occupancyQA.add("occTrackQA/Mean/OccPrimUnfm80", "OccPrimUnfm80", kTH1F, {axisQA1});
    occupancyQA.add("occTrackQA/Mean/OccFV0AUnfm80", "OccFV0AUnfm80", kTH1F, {axisQA1});
    occupancyQA.add("occTrackQA/Mean/OccFV0CUnfm80", "OccFV0CUnfm80", kTH1F, {axisQA1});
    occupancyQA.add("occTrackQA/Mean/OccFT0AUnfm80", "OccFT0AUnfm80", kTH1F, {axisQA1});
    occupancyQA.add("occTrackQA/Mean/OccFT0CUnfm80", "OccFT0CUnfm80", kTH1F, {axisQA1});
    occupancyQA.add("occTrackQA/Mean/OccFDDAUnfm80", "OccFDDAUnfm80", kTH1F, {axisQA1});
    occupancyQA.add("occTrackQA/Mean/OccFDDCUnfm80", "OccFDDCUnfm80", kTH1F, {axisQA1});

    occupancyQA.add("occTrackQA/Mean/OccNTrackITSUnfm80", "OccNTrackITSUnfm80", kTH1F, {axisQA1});
    occupancyQA.add("occTrackQA/Mean/OccNTrackTPCUnfm80", "OccNTrackTPCUnfm80", kTH1F, {axisQA1});
    occupancyQA.add("occTrackQA/Mean/OccNTrackTRDUnfm80", "OccNTrackTRDUnfm80", kTH1F, {axisQA1});
    occupancyQA.add("occTrackQA/Mean/OccNTrackTOFUnfm80", "OccNTrackTOFUnfm80", kTH1F, {axisQA1});
    occupancyQA.add("occTrackQA/Mean/OccNTrackSizeUnfm80", "OccNTrackSizeUnfm80", kTH1F, {axisQA1});
    occupancyQA.add("occTrackQA/Mean/OccNTrackTPCAUnfm80", "OccNTrackTPCAUnfm80", kTH1F, {axisQA1});
    occupancyQA.add("occTrackQA/Mean/OccNTrackTPCCUnfm80", "OccNTrackTPCCUnfm80", kTH1F, {axisQA1});
    occupancyQA.add("occTrackQA/Mean/OccNTrackITSTPCUnfm80", "OccNTrackITSTPCUnfm80", kTH1F, {axisQA1});
    occupancyQA.add("occTrackQA/Mean/OccNTrackITSTPCAUnfm80", "OccNTrackITSTPCAUnfm80", kTH1F, {axisQA1});
    occupancyQA.add("occTrackQA/Mean/OccNTrackITSTPCCUnfm80", "OccNTrackITSTPCCUnfm80", kTH1F, {axisQA1});

    occupancyQA.add("occTrackQA/Mean/OccMultNTracksHasITSUnfm80", "OccMultNTracksHasITSUnfm80", kTH1F, {axisQA1});
    occupancyQA.add("occTrackQA/Mean/OccMultNTracksHasTPCUnfm80", "OccMultNTracksHasTPCUnfm80", kTH1F, {axisQA1});
    occupancyQA.add("occTrackQA/Mean/OccMultNTracksHasTOFUnfm80", "OccMultNTracksHasTOFUnfm80", kTH1F, {axisQA1});
    occupancyQA.add("occTrackQA/Mean/OccMultNTracksHasTRDUnfm80", "OccMultNTracksHasTRDUnfm80", kTH1F, {axisQA1});
    occupancyQA.add("occTrackQA/Mean/OccMultNTracksITSOnlyUnfm80", "OccMultNTracksITSOnlyUnfm80", kTH1F, {axisQA1});
    occupancyQA.add("occTrackQA/Mean/OccMultNTracksTPCOnlyUnfm80", "OccMultNTracksTPCOnlyUnfm80", kTH1F, {axisQA1});
    occupancyQA.add("occTrackQA/Mean/OccMultNTracksITSTPCUnfm80", "OccMultNTracksITSTPCUnfm80", kTH1F, {axisQA1});
    occupancyQA.add("occTrackQA/Mean/OccMultAllTracksTPCOnlyUnfm80", "OccMultAllTracksTPCOnlyUnfm80", kTH1F, {axisQA1});

    occupancyQA.add("occTrackQA/Mean/OccRobustT0V0PrimUnfm80", "OccRobustT0V0PrimUnfm80", kTH1F, {axisQA1});
    occupancyQA.add("occTrackQA/Mean/OccRobustFDDT0V0PrimUnfm80", "OccRobustFDDT0V0PrimUnfm80", kTH1F, {axisQA1});
    occupancyQA.add("occTrackQA/Mean/OccRobustNtrackDetUnfm80", "OccRobustNtrackDetUnfm80", kTH1F, {axisQA1});
    occupancyQA.add("occTrackQA/Mean/OccRobustMultExtraTableUnfm80", "OccRobustMultExtraTableUnfm80", kTH1F, {axisQA1});

    occupancyQA.addClone("occTrackQA/Mean/", "occTrackQA/WeightMean/");

    if (fillQA1) {
      occupancyQA.add("occTrackQA/LogRatio/RobustT0V0Prim/Mean/OccPrimUnfm80", "OccPrimUnfm80", kTH1F, {axisQA2});
      occupancyQA.add("occTrackQA/LogRatio/RobustT0V0Prim/Mean/OccFV0AUnfm80", "OccFV0AUnfm80", kTH1F, {axisQA2});
      occupancyQA.add("occTrackQA/LogRatio/RobustT0V0Prim/Mean/OccFV0CUnfm80", "OccFV0CUnfm80", kTH1F, {axisQA2});
      occupancyQA.add("occTrackQA/LogRatio/RobustT0V0Prim/Mean/OccFT0AUnfm80", "OccFT0AUnfm80", kTH1F, {axisQA2});
      occupancyQA.add("occTrackQA/LogRatio/RobustT0V0Prim/Mean/OccFT0CUnfm80", "OccFT0CUnfm80", kTH1F, {axisQA2});
      occupancyQA.add("occTrackQA/LogRatio/RobustT0V0Prim/Mean/OccFDDAUnfm80", "OccFDDAUnfm80", kTH1F, {axisQA2});
      occupancyQA.add("occTrackQA/LogRatio/RobustT0V0Prim/Mean/OccFDDCUnfm80", "OccFDDCUnfm80", kTH1F, {axisQA2});

      occupancyQA.add("occTrackQA/LogRatio/RobustT0V0Prim/Mean/OccNTrackITSUnfm80", "OccNTrackITSUnfm80", kTH1F, {axisQA2});
      occupancyQA.add("occTrackQA/LogRatio/RobustT0V0Prim/Mean/OccNTrackTPCUnfm80", "OccNTrackTPCUnfm80", kTH1F, {axisQA2});
      occupancyQA.add("occTrackQA/LogRatio/RobustT0V0Prim/Mean/OccNTrackTRDUnfm80", "OccNTrackTRDUnfm80", kTH1F, {axisQA2});
      occupancyQA.add("occTrackQA/LogRatio/RobustT0V0Prim/Mean/OccNTrackTOFUnfm80", "OccNTrackTOFUnfm80", kTH1F, {axisQA2});
      occupancyQA.add("occTrackQA/LogRatio/RobustT0V0Prim/Mean/OccNTrackSizeUnfm80", "OccNTrackSizeUnfm80", kTH1F, {axisQA2});
      occupancyQA.add("occTrackQA/LogRatio/RobustT0V0Prim/Mean/OccNTrackTPCAUnfm80", "OccNTrackTPCAUnfm80", kTH1F, {axisQA2});
      occupancyQA.add("occTrackQA/LogRatio/RobustT0V0Prim/Mean/OccNTrackTPCCUnfm80", "OccNTrackTPCCUnfm80", kTH1F, {axisQA2});
      occupancyQA.add("occTrackQA/LogRatio/RobustT0V0Prim/Mean/OccNTrackITSTPCUnfm80", "OccNTrackITSTPCUnfm80", kTH1F, {axisQA2});
      occupancyQA.add("occTrackQA/LogRatio/RobustT0V0Prim/Mean/OccNTrackITSTPCAUnfm80", "OccNTrackITSTPCAUnfm80", kTH1F, {axisQA2});
      occupancyQA.add("occTrackQA/LogRatio/RobustT0V0Prim/Mean/OccNTrackITSTPCCUnfm80", "OccNTrackITSTPCCUnfm80", kTH1F, {axisQA2});

      occupancyQA.add("occTrackQA/LogRatio/RobustT0V0Prim/Mean/OccMultNTracksHasITSUnfm80", "OccMultNTracksHasITSUnfm80", kTH1F, {axisQA2});
      occupancyQA.add("occTrackQA/LogRatio/RobustT0V0Prim/Mean/OccMultNTracksHasTPCUnfm80", "OccMultNTracksHasTPCUnfm80", kTH1F, {axisQA2});
      occupancyQA.add("occTrackQA/LogRatio/RobustT0V0Prim/Mean/OccMultNTracksHasTOFUnfm80", "OccMultNTracksHasTOFUnfm80", kTH1F, {axisQA2});
      occupancyQA.add("occTrackQA/LogRatio/RobustT0V0Prim/Mean/OccMultNTracksHasTRDUnfm80", "OccMultNTracksHasTRDUnfm80", kTH1F, {axisQA2});
      occupancyQA.add("occTrackQA/LogRatio/RobustT0V0Prim/Mean/OccMultNTracksITSOnlyUnfm80", "OccMultNTracksITSOnlyUnfm80", kTH1F, {axisQA2});
      occupancyQA.add("occTrackQA/LogRatio/RobustT0V0Prim/Mean/OccMultNTracksTPCOnlyUnfm80", "OccMultNTracksTPCOnlyUnfm80", kTH1F, {axisQA2});
      occupancyQA.add("occTrackQA/LogRatio/RobustT0V0Prim/Mean/OccMultNTracksITSTPCUnfm80", "OccMultNTracksITSTPCUnfm80", kTH1F, {axisQA2});
      occupancyQA.add("occTrackQA/LogRatio/RobustT0V0Prim/Mean/OccMultAllTracksTPCOnlyUnfm80", "OccMultAllTracksTPCOnlyUnfm80", kTH1F, {axisQA2});

      occupancyQA.add("occTrackQA/LogRatio/RobustT0V0Prim/Mean/OccRobustT0V0PrimUnfm80", "OccRobustT0V0PrimUnfm80", kTH1F, {axisQA2});
      occupancyQA.add("occTrackQA/LogRatio/RobustT0V0Prim/Mean/OccRobustFDDT0V0PrimUnfm80", "OccRobustFDDT0V0PrimUnfm80", kTH1F, {axisQA2});
      occupancyQA.add("occTrackQA/LogRatio/RobustT0V0Prim/Mean/OccRobustNtrackDetUnfm80", "OccRobustNtrackDetUnfm80", kTH1F, {axisQA2});
      occupancyQA.add("occTrackQA/LogRatio/RobustT0V0Prim/Mean/OccRobustMultExtraTableUnfm80", "OccRobustMultExtraTableUnfm80", kTH1F, {axisQA2});

      occupancyQA.addClone("occTrackQA/LogRatio/RobustT0V0Prim/Mean/", "occTrackQA/LogRatio/RobustT0V0Prim/WeightMean/");
      occupancyQA.addClone("occTrackQA/LogRatio/RobustT0V0Prim/WeightMean/", "occTrackQA/LogRatio/weightRobustT0V0Prim/WeightMean/");
    }

    if (fillQA2) {
      occupancyQA.add("occTrackQA/Condition1/RobustT0V0Prim/Mean/OccPrimUnfm80", "OccPrimUnfm80", kTH1F, {axisQA3});
      occupancyQA.add("occTrackQA/Condition1/RobustT0V0Prim/Mean/OccFV0AUnfm80", "OccFV0AUnfm80", kTH1F, {axisQA3});
      occupancyQA.add("occTrackQA/Condition1/RobustT0V0Prim/Mean/OccFV0CUnfm80", "OccFV0CUnfm80", kTH1F, {axisQA3});
      occupancyQA.add("occTrackQA/Condition1/RobustT0V0Prim/Mean/OccFT0AUnfm80", "OccFT0AUnfm80", kTH1F, {axisQA3});
      occupancyQA.add("occTrackQA/Condition1/RobustT0V0Prim/Mean/OccFT0CUnfm80", "OccFT0CUnfm80", kTH1F, {axisQA3});
      occupancyQA.add("occTrackQA/Condition1/RobustT0V0Prim/Mean/OccFDDAUnfm80", "OccFDDAUnfm80", kTH1F, {axisQA3});
      occupancyQA.add("occTrackQA/Condition1/RobustT0V0Prim/Mean/OccFDDCUnfm80", "OccFDDCUnfm80", kTH1F, {axisQA3});

      occupancyQA.add("occTrackQA/Condition1/RobustT0V0Prim/Mean/OccNTrackITSUnfm80", "OccNTrackITSUnfm80", kTH1F, {axisQA3});
      occupancyQA.add("occTrackQA/Condition1/RobustT0V0Prim/Mean/OccNTrackTPCUnfm80", "OccNTrackTPCUnfm80", kTH1F, {axisQA3});
      occupancyQA.add("occTrackQA/Condition1/RobustT0V0Prim/Mean/OccNTrackTRDUnfm80", "OccNTrackTRDUnfm80", kTH1F, {axisQA3});
      occupancyQA.add("occTrackQA/Condition1/RobustT0V0Prim/Mean/OccNTrackTOFUnfm80", "OccNTrackTOFUnfm80", kTH1F, {axisQA3});
      occupancyQA.add("occTrackQA/Condition1/RobustT0V0Prim/Mean/OccNTrackSizeUnfm80", "OccNTrackSizeUnfm80", kTH1F, {axisQA3});
      occupancyQA.add("occTrackQA/Condition1/RobustT0V0Prim/Mean/OccNTrackTPCAUnfm80", "OccNTrackTPCAUnfm80", kTH1F, {axisQA3});
      occupancyQA.add("occTrackQA/Condition1/RobustT0V0Prim/Mean/OccNTrackTPCCUnfm80", "OccNTrackTPCCUnfm80", kTH1F, {axisQA3});
      occupancyQA.add("occTrackQA/Condition1/RobustT0V0Prim/Mean/OccNTrackITSTPCUnfm80", "OccNTrackITSTPCUnfm80", kTH1F, {axisQA3});
      occupancyQA.add("occTrackQA/Condition1/RobustT0V0Prim/Mean/OccNTrackITSTPCAUnfm80", "OccNTrackITSTPCAUnfm80", kTH1F, {axisQA3});
      occupancyQA.add("occTrackQA/Condition1/RobustT0V0Prim/Mean/OccNTrackITSTPCCUnfm80", "OccNTrackITSTPCCUnfm80", kTH1F, {axisQA3});

      occupancyQA.add("occTrackQA/Condition1/RobustT0V0Prim/Mean/OccMultNTracksHasITSUnfm80", "OccMultNTracksHasITSUnfm80", kTH1F, {axisQA3});
      occupancyQA.add("occTrackQA/Condition1/RobustT0V0Prim/Mean/OccMultNTracksHasTPCUnfm80", "OccMultNTracksHasTPCUnfm80", kTH1F, {axisQA3});
      occupancyQA.add("occTrackQA/Condition1/RobustT0V0Prim/Mean/OccMultNTracksHasTOFUnfm80", "OccMultNTracksHasTOFUnfm80", kTH1F, {axisQA3});
      occupancyQA.add("occTrackQA/Condition1/RobustT0V0Prim/Mean/OccMultNTracksHasTRDUnfm80", "OccMultNTracksHasTRDUnfm80", kTH1F, {axisQA3});
      occupancyQA.add("occTrackQA/Condition1/RobustT0V0Prim/Mean/OccMultNTracksITSOnlyUnfm80", "OccMultNTracksITSOnlyUnfm80", kTH1F, {axisQA3});
      occupancyQA.add("occTrackQA/Condition1/RobustT0V0Prim/Mean/OccMultNTracksTPCOnlyUnfm80", "OccMultNTracksTPCOnlyUnfm80", kTH1F, {axisQA3});
      occupancyQA.add("occTrackQA/Condition1/RobustT0V0Prim/Mean/OccMultNTracksITSTPCUnfm80", "OccMultNTracksITSTPCUnfm80", kTH1F, {axisQA3});
      occupancyQA.add("occTrackQA/Condition1/RobustT0V0Prim/Mean/OccMultAllTracksTPCOnlyUnfm80", "OccMultAllTracksTPCOnlyUnfm80", kTH1F, {axisQA3});

      occupancyQA.add("occTrackQA/Condition1/RobustT0V0Prim/Mean/OccRobustT0V0PrimUnfm80", "OccRobustT0V0PrimUnfm80", kTH1F, {axisQA3});
      occupancyQA.add("occTrackQA/Condition1/RobustT0V0Prim/Mean/OccRobustFDDT0V0PrimUnfm80", "OccRobustFDDT0V0PrimUnfm80", kTH1F, {axisQA3});
      occupancyQA.add("occTrackQA/Condition1/RobustT0V0Prim/Mean/OccRobustNtrackDetUnfm80", "OccRobustNtrackDetUnfm80", kTH1F, {axisQA3});
      occupancyQA.add("occTrackQA/Condition1/RobustT0V0Prim/Mean/OccRobustMultExtraTableUnfm80", "OccRobustMultExtraTableUnfm80", kTH1F, {axisQA3});

      occupancyQA.addClone("occTrackQA/Condition1/RobustT0V0Prim/Mean/", "occTrackQA/Condition1/RobustT0V0Prim/WeightMean/");
      occupancyQA.addClone("occTrackQA/Condition1/RobustT0V0Prim/WeightMean/", "occTrackQA/Condition1/weightRobustT0V0Prim/WeightMean/");

      occupancyQA.addClone("occTrackQA/Condition1/", "occTrackQA/Condition2/");
      occupancyQA.addClone("occTrackQA/Condition1/", "occTrackQA/Condition3/");
      occupancyQA.addClone("occTrackQA/Condition1/", "occTrackQA/Condition4/");
    }

    occupancyQA.add("h_DFcount_Lvl0", "h_DFcount_Lvl0", kTH1F, {{13, -1, 12}});
    occupancyQA.add("h_DFcount_Lvl1", "h_DFcount_Lvl1", kTH1F, {{13, -1, 12}});
    occupancyQA.add("h_DFcount_Lvl2", "h_DFcount_Lvl2", kTH1F, {{13, -1, 12}});

    occupancyQA.print();
  }

  enum OccNamesEnum {
    kOccPrimUnfm80 = 0,
    kOccFV0AUnfm80,
    kOccFV0CUnfm80,
    kOccFT0AUnfm80,
    kOccFT0CUnfm80,
    kOccFDDAUnfm80,
    kOccFDDCUnfm80,

    kOccNTrackITSUnfm80,
    kOccNTrackTPCUnfm80,
    kOccNTrackTRDUnfm80,
    kOccNTrackTOFUnfm80,
    kOccNTrackSizeUnfm80,
    kOccNTrackTPCAUnfm80,
    kOccNTrackTPCCUnfm80,
    kOccNTrackITSTPCUnfm80,
    kOccNTrackITSTPCAUnfm80,
    kOccNTrackITSTPCCUnfm80,

    kOccMultNTracksHasITSUnfm80,
    kOccMultNTracksHasTPCUnfm80,
    kOccMultNTracksHasTOFUnfm80,
    kOccMultNTracksHasTRDUnfm80,
    kOccMultNTracksITSOnlyUnfm80,
    kOccMultNTracksTPCOnlyUnfm80,
    kOccMultNTracksITSTPCUnfm80,
    kOccMultAllTracksTPCOnlyUnfm80,

    kOccRobustT0V0PrimUnfm80,
    kOccRobustFDDT0V0PrimUnfm80,
    kOccRobustNtrackDetUnfm80,
    kOccRobustMultTableUnfm80
  };

  static constexpr std::string_view OccNames[]{
    "OccPrimUnfm80",
    "OccFV0AUnfm80",
    "OccFV0CUnfm80",
    "OccFT0AUnfm80",
    "OccFT0CUnfm80",
    "OccFDDAUnfm80",
    "OccFDDCUnfm80",

    "OccNTrackITSUnfm80",
    "OccNTrackTPCUnfm80",
    "OccNTrackTRDUnfm80",
    "OccNTrackTOFUnfm80",
    "OccNTrackSizeUnfm80",
    "OccNTrackTPCAUnfm80",
    "OccNTrackTPCCUnfm80",
    "OccNTrackITSTPCUnfm80",
    "OccNTrackITSTPCAUnfm80",
    "OccNTrackITSTPCCUnfm80",

    "OccMultNTracksHasITSUnfm80",
    "OccMultNTracksHasTPCUnfm80",
    "OccMultNTracksHasTOFUnfm80",
    "OccMultNTracksHasTRDUnfm80",
    "OccMultNTracksITSOnlyUnfm80",
    "OccMultNTracksTPCOnlyUnfm80",
    "OccMultNTracksITSTPCUnfm80",
    "OccMultAllTracksTPCOnlyUnfm80",

    "OccRobustT0V0PrimUnfm80",
    "OccRobustFDDT0V0PrimUnfm80",
    "OccRobustNtrackDetUnfm80",
    "OccRobustMultExtraTableUnfm80"};

  enum OccDirEnum {
    kMean = 0,
    kWeightMean,
    kLogRatio,
    kRobustT0V0Prim,
    kWeightRobustT0V0Prim,
    kCondition1,
    kCondition2,
    kCondition3,
    kCondition4
  };

  static constexpr std::string_view OccDire[] = {
    "Mean/",
    "WeightMean/",
    "LogRatio/",
    "RobustT0V0Prim/",
    "weightRobustT0V0Prim/",
    "Condition1/",
    "Condition2/",
    "Condition3/",
    "Condition4/"};

  void getRunInfo(const int& run, int& nBCsPerTF, int64_t& bcSOR)
  {
    auto runDuration = ccdb->getRunDuration(run, true);
    int64_t tsSOR = runDuration.first;
    auto ctpx = ccdb->getForTimeStamp<std::vector<int64_t>>("CTP/Calib/OrbitReset", tsSOR);
    int64_t tsOrbitReset = (*ctpx)[0];
    uint32_t nOrbitsPerTF = run < cfgNOrbitsPerTF0RunValue ? cfgNOrbitsPerTF1TrueValue : cfgNOrbitsPerTF2FalseValue;
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

  template <int occMode, int occRobustMode, int occName>
  void fillQAInfo(const float& occValue, const float& occRobustValue)
  {
    occupancyQA.fill(HIST("occTrackQA/") + HIST(OccDire[occMode]) + HIST(OccNames[occName]), occValue);
    if (fillQA1) {
      occupancyQA.fill(HIST("occTrackQA/LogRatio/") + HIST(OccDire[occRobustMode]) + HIST(OccDire[occMode]) + HIST(OccNames[occName]), std::log(std::abs(occValue / occRobustValue)));
      if (fillQA2) {
        int two = 2, twenty = 20, fifty = 50, twoHundred = 200;
        if (std::abs(std::log(occValue / occRobustValue)) < two) { // conditional filling start
          occupancyQA.fill(HIST("occTrackQA/Condition1/") + HIST(OccDire[occRobustMode]) + HIST(OccDire[occMode]) + HIST(OccNames[occName]), (std::log(occValue / occRobustValue)) * std::sqrt(occValue + occRobustValue));
          if (std::abs(occRobustValue + occValue) > twoHundred) {
            occupancyQA.fill(HIST("occTrackQA/Condition4/") + HIST(OccDire[occRobustMode]) + HIST(OccDire[occMode]) + HIST(OccNames[occName]), (std::log(occValue / occRobustValue)) * std::sqrt(occValue + occRobustValue));
            occupancyQA.fill(HIST("occTrackQA/Condition3/") + HIST(OccDire[occRobustMode]) + HIST(OccDire[occMode]) + HIST(OccNames[occName]), (std::log(occValue / occRobustValue)) * std::sqrt(occValue + occRobustValue));
            occupancyQA.fill(HIST("occTrackQA/Condition2/") + HIST(OccDire[occRobustMode]) + HIST(OccDire[occMode]) + HIST(OccNames[occName]), (std::log(occValue / occRobustValue)) * std::sqrt(occValue + occRobustValue));
          } else if (std::abs(occRobustValue + occValue) > fifty) {
            occupancyQA.fill(HIST("occTrackQA/Condition3/") + HIST(OccDire[occRobustMode]) + HIST(OccDire[occMode]) + HIST(OccNames[occName]), (std::log(occValue / occRobustValue)) * std::sqrt(occValue + occRobustValue));
            occupancyQA.fill(HIST("occTrackQA/Condition2/") + HIST(OccDire[occRobustMode]) + HIST(OccDire[occMode]) + HIST(OccNames[occName]), (std::log(occValue / occRobustValue)) * std::sqrt(occValue + occRobustValue));
          } else if (std::abs(occRobustValue + occValue) > twenty) {
            occupancyQA.fill(HIST("occTrackQA/Condition2/") + HIST(OccDire[occRobustMode]) + HIST(OccDire[occMode]) + HIST(OccNames[occName]), (std::log(occValue / occRobustValue)) * std::sqrt(occValue + occRobustValue));
          }
        } // conditional filling end
      }
    }
  }

  using MyTracksQA = aod::TracksQAVersion; // using MyTracksQA = aod::TracksQA_002;

  // Process the Data
  int dfCount = 0;
  int32_t nBCsPerTF = -999;
  int64_t bcSOR = -999;
  int lastRun = -999;

  uint64_t time = -1;
  int64_t tfIdThis = -1;
  int bcInTF = -1;

  enum ProcessTags {
    kProcessNothing = 0,
    kProcessOnlyOccPrim,
    kProcessOnlyOccT0V0,
    kProcessOnlyOccFDD,
    kProcessOnlyOccNtrackDet,
    kProcessOnlyOccMultExtra,
    kProcessOnlyRobustT0V0Prim,
    kProcessOnlyRobustFDDT0V0Prim,
    kProcessOnlyRobustNtrackDet,
    kProcessOnlyRobustMultExtra,
    kProcessFullOccTableProducer
  };

  enum FillMode {
    checkTableMode = 0,
    checkQAMode,
    doNotFill,
    fillOccRobustT0V0dependentQA,
    fillMeanOccTable,
    fillWeightMeanOccTable
  };

  std::vector<std::array<int64_t, 2>> trackQAGIListforTMOList;
  template <int processMode, int meanTableMode, int weightMeanTableMode, int qaMode, typename B, typename C, typename T, typename U, typename O, typename V>
  void executeTrackOccProducerProcessing(B const& BCs, C const& collisions, T const& tracks, U const& tracksQA, O const& occsRobustT0V0Prim, V const& occs, bool const& executeInThisBlock)
  {
    if (collisions.size() == 0) {
      return;
    }

    if (meanTableMode == checkTableMode) {
      if (buildFlag00MeanTable) {
        executeTrackOccProducerProcessing<processMode, fillMeanOccTable, weightMeanTableMode, qaMode>(BCs, collisions, tracks, tracksQA, occsRobustT0V0Prim, occs, executeInThisBlock);
      } else {
        executeTrackOccProducerProcessing<processMode, doNotFill, weightMeanTableMode, qaMode>(BCs, collisions, tracks, tracksQA, occsRobustT0V0Prim, occs, executeInThisBlock);
      }
    }
    if constexpr (meanTableMode == checkTableMode) {
      return;
    }

    if (weightMeanTableMode == checkTableMode) {
      if (buildFlag01WeightMeanTable) {
        executeTrackOccProducerProcessing<processMode, meanTableMode, fillWeightMeanOccTable, qaMode>(BCs, collisions, tracks, tracksQA, occsRobustT0V0Prim, occs, executeInThisBlock);
      } else {
        executeTrackOccProducerProcessing<processMode, meanTableMode, doNotFill, qaMode>(BCs, collisions, tracks, tracksQA, occsRobustT0V0Prim, occs, executeInThisBlock);
      }
    }
    if constexpr (weightMeanTableMode == checkTableMode) {
      return;
    }

    if (qaMode == checkQAMode) {
      if (fillQA1 || fillQA2) {
        if (occsRobustT0V0Prim.size() == 0) {
          LOG(error) << "DEBUG :: ERROR ERROR ERROR :: OccsRobustT0V0Prim.size() == 0 :: Check \"occupancy-table-producer\" for \"buildOnlyOccsT0V0Prim == true\" & \"processOnlyOccT0V0PrimUnfm == true\"";
          return;
        }
        executeTrackOccProducerProcessing<processMode, meanTableMode, weightMeanTableMode, fillOccRobustT0V0dependentQA>(BCs, collisions, tracks, tracksQA, occsRobustT0V0Prim, occs, executeInThisBlock);
      } else {
        executeTrackOccProducerProcessing<processMode, meanTableMode, weightMeanTableMode, doNotFill>(BCs, collisions, tracks, tracksQA, occsRobustT0V0Prim, occs, executeInThisBlock);
      }
    }
    if constexpr (qaMode == checkQAMode) {
      return;
    }

    // BCs.bindExternalIndices(&occsDet);
    // BCs.bindExternalIndices(&occsNTrackDet);
    // BCs.bindExternalIndices(&occsRobust);

    if constexpr (meanTableMode == checkTableMode || weightMeanTableMode == checkTableMode || qaMode == checkQAMode) {
      return;
    } else {
      occupancyQA.fill(HIST("h_DFcount_Lvl2"), processMode);

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

      float meanOccPrimUnfm80 = 0;
      float meanOccFV0AUnfm80 = 0;
      float meanOccFV0CUnfm80 = 0;
      float meanOccFT0AUnfm80 = 0;
      float meanOccFT0CUnfm80 = 0;
      float meanOccFDDAUnfm80 = 0;
      float meanOccFDDCUnfm80 = 0;

      float meanOccNTrackITSUnfm80 = 0;
      float meanOccNTrackTPCUnfm80 = 0;
      float meanOccNTrackTRDUnfm80 = 0;
      float meanOccNTrackTOFUnfm80 = 0;
      float meanOccNTrackSizeUnfm80 = 0;
      float meanOccNTrackTPCAUnfm80 = 0;
      float meanOccNTrackTPCCUnfm80 = 0;
      float meanOccNTrackITSTPCUnfm80 = 0;
      float meanOccNTrackITSTPCAUnfm80 = 0;
      float meanOccNTrackITSTPCCUnfm80 = 0;

      float meanOccMultNTracksHasITSUnfm80 = 0;
      float meanOccMultNTracksHasTPCUnfm80 = 0;
      float meanOccMultNTracksHasTOFUnfm80 = 0;
      float meanOccMultNTracksHasTRDUnfm80 = 0;
      float meanOccMultNTracksITSOnlyUnfm80 = 0;
      float meanOccMultNTracksTPCOnlyUnfm80 = 0;
      float meanOccMultNTracksITSTPCUnfm80 = 0;
      float meanOccMultAllTracksTPCOnlyUnfm80 = 0;

      float meanOccRobustT0V0PrimUnfm80 = 0;
      float meanOccRobustFDDT0V0PrimUnfm80 = 0;
      float meanOccRobustNtrackDetUnfm80 = 0;
      float meanOccRobustMultTableUnfm80 = 0;

      float weightMeanOccPrimUnfm80 = 0;
      float weightMeanOccFV0AUnfm80 = 0;
      float weightMeanOccFV0CUnfm80 = 0;
      float weightMeanOccFT0AUnfm80 = 0;
      float weightMeanOccFT0CUnfm80 = 0;
      float weightMeanOccFDDAUnfm80 = 0;
      float weightMeanOccFDDCUnfm80 = 0;

      float weightMeanOccNTrackITSUnfm80 = 0;
      float weightMeanOccNTrackTPCUnfm80 = 0;
      float weightMeanOccNTrackTRDUnfm80 = 0;
      float weightMeanOccNTrackTOFUnfm80 = 0;
      float weightMeanOccNTrackSizeUnfm80 = 0;
      float weightMeanOccNTrackTPCAUnfm80 = 0;
      float weightMeanOccNTrackTPCCUnfm80 = 0;
      float weightMeanOccNTrackITSTPCUnfm80 = 0;
      float weightMeanOccNTrackITSTPCAUnfm80 = 0;
      float weightMeanOccNTrackITSTPCCUnfm80 = 0;

      float weightMeanOccMultNTracksHasITSUnfm80 = 0;
      float weightMeanOccMultNTracksHasTPCUnfm80 = 0;
      float weightMeanOccMultNTracksHasTOFUnfm80 = 0;
      float weightMeanOccMultNTracksHasTRDUnfm80 = 0;
      float weightMeanOccMultNTracksITSOnlyUnfm80 = 0;
      float weightMeanOccMultNTracksTPCOnlyUnfm80 = 0;
      float weightMeanOccMultNTracksITSTPCUnfm80 = 0;
      float weightMeanOccMultAllTracksTPCOnlyUnfm80 = 0;

      float weightMeanOccRobustT0V0PrimUnfm80 = 0;
      float weightMeanOccRobustFDDT0V0PrimUnfm80 = 0;
      float weightMeanOccRobustNtrackDetUnfm80 = 0;
      float weightMeanOccRobustMultTableUnfm80 = 0;

      int trackTMOcounter = -1;
      trackQAGIListforTMOList.clear();

      for (const auto& trackQA : tracksQA) {
        auto const& track = trackQA.template track_as<T>();
        auto collision = collisions.begin();

        hasCollision = false;
        isAmbgTrack = false;

        if (track.collisionId() >= 0) {                 // track has collision
          collision = track.template collision_as<C>(); // It will build but crash while running for tracks with track.collisionId()= -1;//ambg tracks/orphan tracks
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
            bc = collision.template bc_as<B>();
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
          auto occsList = occs.iteratorAt(bc.occId());

          if constexpr (qaMode == fillOccRobustT0V0dependentQA) {
            std::copy(occsRobustT0V0Prim.iteratorAt(bc.occId()).occRobustT0V0PrimUnfm80().begin(), occsRobustT0V0Prim.iteratorAt(bc.occId()).occRobustT0V0PrimUnfm80().end(), occRobustT0V0PrimUnfm80.begin());
          }

          if constexpr (processMode == kProcessFullOccTableProducer || processMode == kProcessOnlyOccPrim) {
            std::copy(occsList.occPrimUnfm80().begin(), occsList.occPrimUnfm80().end(), occPrimUnfm80.begin());
          }
          if constexpr (processMode == kProcessFullOccTableProducer || processMode == kProcessOnlyOccT0V0) {
            std::copy(occsList.occFV0AUnfm80().begin(), occsList.occFV0AUnfm80().end(), occFV0AUnfm80.begin());
            std::copy(occsList.occFV0CUnfm80().begin(), occsList.occFV0CUnfm80().end(), occFV0CUnfm80.begin());
            std::copy(occsList.occFT0AUnfm80().begin(), occsList.occFT0AUnfm80().end(), occFT0AUnfm80.begin());
            std::copy(occsList.occFT0CUnfm80().begin(), occsList.occFT0CUnfm80().end(), occFT0CUnfm80.begin());
          }
          if constexpr (processMode == kProcessFullOccTableProducer || processMode == kProcessOnlyOccFDD) {
            std::copy(occsList.occFDDAUnfm80().begin(), occsList.occFDDAUnfm80().end(), occFDDAUnfm80.begin());
            std::copy(occsList.occFDDCUnfm80().begin(), occsList.occFDDCUnfm80().end(), occFDDCUnfm80.begin());
          }

          if constexpr (processMode == kProcessFullOccTableProducer || processMode == kProcessOnlyOccNtrackDet) {
            std::copy(occsList.occNTrackITSUnfm80().begin(), occsList.occNTrackITSUnfm80().end(), occNTrackITSUnfm80.begin());
            std::copy(occsList.occNTrackTPCUnfm80().begin(), occsList.occNTrackTPCUnfm80().end(), occNTrackTPCUnfm80.begin());
            std::copy(occsList.occNTrackTRDUnfm80().begin(), occsList.occNTrackTRDUnfm80().end(), occNTrackTRDUnfm80.begin());
            std::copy(occsList.occNTrackTOFUnfm80().begin(), occsList.occNTrackTOFUnfm80().end(), occNTrackTOFUnfm80.begin());
            std::copy(occsList.occNTrackSizeUnfm80().begin(), occsList.occNTrackSizeUnfm80().end(), occNTrackSizeUnfm80.begin());
            std::copy(occsList.occNTrackTPCAUnfm80().begin(), occsList.occNTrackTPCAUnfm80().end(), occNTrackTPCAUnfm80.begin());
            std::copy(occsList.occNTrackTPCCUnfm80().begin(), occsList.occNTrackTPCCUnfm80().end(), occNTrackTPCCUnfm80.begin());
            std::copy(occsList.occNTrackITSTPCUnfm80().begin(), occsList.occNTrackITSTPCUnfm80().end(), occNTrackITSTPCUnfm80.begin());
            std::copy(occsList.occNTrackITSTPCAUnfm80().begin(), occsList.occNTrackITSTPCAUnfm80().end(), occNTrackITSTPCAUnfm80.begin());
            std::copy(occsList.occNTrackITSTPCCUnfm80().begin(), occsList.occNTrackITSTPCCUnfm80().end(), occNTrackITSTPCCUnfm80.begin());
          }

          if constexpr (processMode == kProcessFullOccTableProducer || processMode == kProcessOnlyOccMultExtra) {
            std::copy(occsList.occMultNTracksHasITSUnfm80().begin(), occsList.occMultNTracksHasITSUnfm80().end(), occMultNTracksHasITSUnfm80.begin());
            std::copy(occsList.occMultNTracksHasTPCUnfm80().begin(), occsList.occMultNTracksHasTPCUnfm80().end(), occMultNTracksHasTPCUnfm80.begin());
            std::copy(occsList.occMultNTracksHasTOFUnfm80().begin(), occsList.occMultNTracksHasTOFUnfm80().end(), occMultNTracksHasTOFUnfm80.begin());
            std::copy(occsList.occMultNTracksHasTRDUnfm80().begin(), occsList.occMultNTracksHasTRDUnfm80().end(), occMultNTracksHasTRDUnfm80.begin());
            std::copy(occsList.occMultNTracksITSOnlyUnfm80().begin(), occsList.occMultNTracksITSOnlyUnfm80().end(), occMultNTracksITSOnlyUnfm80.begin());
            std::copy(occsList.occMultNTracksTPCOnlyUnfm80().begin(), occsList.occMultNTracksTPCOnlyUnfm80().end(), occMultNTracksTPCOnlyUnfm80.begin());
            std::copy(occsList.occMultNTracksITSTPCUnfm80().begin(), occsList.occMultNTracksITSTPCUnfm80().end(), occMultNTracksITSTPCUnfm80.begin());
            std::copy(occsList.occMultAllTracksTPCOnlyUnfm80().begin(), occsList.occMultAllTracksTPCOnlyUnfm80().end(), occMultAllTracksTPCOnlyUnfm80.begin());
          }
          if constexpr (processMode == kProcessFullOccTableProducer || processMode == kProcessOnlyRobustT0V0Prim) {
            std::copy(occsList.occRobustT0V0PrimUnfm80().begin(), occsList.occRobustT0V0PrimUnfm80().end(), occRobustT0V0PrimUnfm80.begin());
          }
          if constexpr (processMode == kProcessFullOccTableProducer || processMode == kProcessOnlyRobustFDDT0V0Prim) {
            std::copy(occsList.occRobustFDDT0V0PrimUnfm80().begin(), occsList.occRobustFDDT0V0PrimUnfm80().end(), occRobustFDDT0V0PrimUnfm80.begin());
          }
          if constexpr (processMode == kProcessFullOccTableProducer || processMode == kProcessOnlyRobustNtrackDet) {
            std::copy(occsList.occRobustNtrackDetUnfm80().begin(), occsList.occRobustNtrackDetUnfm80().end(), occRobustNtrackDetUnfm80.begin());
          }
          if constexpr (processMode == kProcessFullOccTableProducer || processMode == kProcessOnlyRobustMultExtra) {
            std::copy(occsList.occRobustMultExtraTableUnfm80().begin(), occsList.occRobustMultExtraTableUnfm80().end(), occRobustMultTableUnfm80.begin());
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
        float length = 250.0;
        // clip the result at 250
        if (zBegin > length) {
          zBegin = 250;
        } else if (zBegin < -length) {
          zBegin = -250;
        }

        if (zEnd > length) {
          zEnd = 250;
        } else if (zEnd < -length) {
          zEnd = -250;
        }

        dTbegin = ((length - std::abs(zBegin)) / vdrift) / 0.025;
        dTend = ((length - std::abs(zEnd)) / vdrift) / 0.025;

        bcBegin = bcInTF + dTbegin;
        bcEnd = bcInTF + dTend;

        binBCbegin = bcBegin / 80;
        binBCend = bcEnd / 80;

        // If multiple process are on, fill this table only once
        if (executeInThisBlock) {
          trackTMOcounter++;
          genTmoTrackId(track.globalIndex());
          trackQAGIListforTMOList.push_back({trackQA.globalIndex(), trackTMOcounter});
        }

        if constexpr (qaMode == fillOccRobustT0V0dependentQA) {
          if constexpr (meanTableMode == fillMeanOccTable) {
            meanOccRobustT0V0PrimUnfm80 = getMeanOccupancy(binBCbegin, binBCend, occRobustT0V0PrimUnfm80);
          }
          if constexpr (weightMeanTableMode == fillWeightMeanOccTable) {
            weightMeanOccRobustT0V0PrimUnfm80 = getWeightedMeanOccupancy(binBCbegin, binBCend, occRobustT0V0PrimUnfm80);
          }
        }

        if constexpr (processMode == kProcessFullOccTableProducer || processMode == kProcessOnlyOccPrim) {
          if constexpr (meanTableMode == fillMeanOccTable) {
            meanOccPrimUnfm80 = getMeanOccupancy(binBCbegin, binBCend, occPrimUnfm80);
            genTmoPrim(meanOccPrimUnfm80);
            fillQAInfo<kMean, kRobustT0V0Prim, kOccPrimUnfm80>(meanOccPrimUnfm80, meanOccRobustT0V0PrimUnfm80);
          }
          if constexpr (weightMeanTableMode == fillWeightMeanOccTable) {
            weightMeanOccPrimUnfm80 = getWeightedMeanOccupancy(binBCbegin, binBCend, occPrimUnfm80);
            genTwmoPrim(weightMeanOccPrimUnfm80);
            fillQAInfo<kWeightMean, kRobustT0V0Prim, kOccPrimUnfm80>(weightMeanOccPrimUnfm80, meanOccRobustT0V0PrimUnfm80);
            fillQAInfo<kWeightMean, kWeightRobustT0V0Prim, kOccPrimUnfm80>(weightMeanOccPrimUnfm80, weightMeanOccRobustT0V0PrimUnfm80);
          }
        }

        if constexpr (processMode == kProcessFullOccTableProducer || processMode == kProcessOnlyOccT0V0) {
          if constexpr (meanTableMode == fillMeanOccTable) {
            meanOccFV0AUnfm80 = getMeanOccupancy(binBCbegin, binBCend, occFV0AUnfm80);
            meanOccFV0CUnfm80 = getMeanOccupancy(binBCbegin, binBCend, occFV0CUnfm80);
            meanOccFT0AUnfm80 = getMeanOccupancy(binBCbegin, binBCend, occFT0AUnfm80);
            meanOccFT0CUnfm80 = getMeanOccupancy(binBCbegin, binBCend, occFT0CUnfm80);
            genTmoT0V0(meanOccFV0AUnfm80,
                       meanOccFV0CUnfm80,
                       meanOccFT0AUnfm80,
                       meanOccFT0CUnfm80);
            fillQAInfo<kMean, kRobustT0V0Prim, kOccFV0AUnfm80>(meanOccFV0AUnfm80, meanOccRobustT0V0PrimUnfm80);
            fillQAInfo<kMean, kRobustT0V0Prim, kOccFV0CUnfm80>(meanOccFV0CUnfm80, meanOccRobustT0V0PrimUnfm80);
            fillQAInfo<kMean, kRobustT0V0Prim, kOccFT0AUnfm80>(meanOccFT0AUnfm80, meanOccRobustT0V0PrimUnfm80);
            fillQAInfo<kMean, kRobustT0V0Prim, kOccFT0CUnfm80>(meanOccFT0CUnfm80, meanOccRobustT0V0PrimUnfm80);
          }
          if constexpr (weightMeanTableMode == fillWeightMeanOccTable) {
            weightMeanOccFV0AUnfm80 = getWeightedMeanOccupancy(binBCbegin, binBCend, occFV0AUnfm80);
            weightMeanOccFV0CUnfm80 = getWeightedMeanOccupancy(binBCbegin, binBCend, occFV0CUnfm80);
            weightMeanOccFT0AUnfm80 = getWeightedMeanOccupancy(binBCbegin, binBCend, occFT0AUnfm80);
            weightMeanOccFT0CUnfm80 = getWeightedMeanOccupancy(binBCbegin, binBCend, occFT0CUnfm80);
            genTwmoT0V0(weightMeanOccFV0AUnfm80,
                        weightMeanOccFV0CUnfm80,
                        weightMeanOccFT0AUnfm80,
                        weightMeanOccFT0CUnfm80);
            fillQAInfo<kWeightMean, kRobustT0V0Prim, kOccFV0AUnfm80>(weightMeanOccFV0AUnfm80, meanOccRobustT0V0PrimUnfm80);
            fillQAInfo<kWeightMean, kRobustT0V0Prim, kOccFV0CUnfm80>(weightMeanOccFV0CUnfm80, meanOccRobustT0V0PrimUnfm80);
            fillQAInfo<kWeightMean, kRobustT0V0Prim, kOccFT0AUnfm80>(weightMeanOccFT0AUnfm80, meanOccRobustT0V0PrimUnfm80);
            fillQAInfo<kWeightMean, kRobustT0V0Prim, kOccFT0CUnfm80>(weightMeanOccFT0CUnfm80, meanOccRobustT0V0PrimUnfm80);

            fillQAInfo<kWeightMean, kWeightRobustT0V0Prim, kOccFV0AUnfm80>(weightMeanOccFV0AUnfm80, weightMeanOccRobustT0V0PrimUnfm80);
            fillQAInfo<kWeightMean, kWeightRobustT0V0Prim, kOccFV0CUnfm80>(weightMeanOccFV0CUnfm80, weightMeanOccRobustT0V0PrimUnfm80);
            fillQAInfo<kWeightMean, kWeightRobustT0V0Prim, kOccFT0AUnfm80>(weightMeanOccFT0AUnfm80, weightMeanOccRobustT0V0PrimUnfm80);
            fillQAInfo<kWeightMean, kWeightRobustT0V0Prim, kOccFT0CUnfm80>(weightMeanOccFT0CUnfm80, weightMeanOccRobustT0V0PrimUnfm80);
          }
        }

        if constexpr (processMode == kProcessFullOccTableProducer || processMode == kProcessOnlyOccFDD) {
          if constexpr (meanTableMode == fillMeanOccTable) {
            meanOccFDDAUnfm80 = getMeanOccupancy(binBCbegin, binBCend, occFDDAUnfm80);
            meanOccFDDCUnfm80 = getMeanOccupancy(binBCbegin, binBCend, occFDDCUnfm80);
            genTmoFDD(meanOccFDDAUnfm80,
                      meanOccFDDCUnfm80);
            fillQAInfo<kMean, kRobustT0V0Prim, kOccFDDAUnfm80>(meanOccFDDAUnfm80, meanOccRobustT0V0PrimUnfm80);
            fillQAInfo<kMean, kRobustT0V0Prim, kOccFDDCUnfm80>(meanOccFDDCUnfm80, meanOccRobustT0V0PrimUnfm80);
          }
          if constexpr (weightMeanTableMode == fillWeightMeanOccTable) {
            weightMeanOccFDDAUnfm80 = getWeightedMeanOccupancy(binBCbegin, binBCend, occFDDAUnfm80);
            weightMeanOccFDDCUnfm80 = getWeightedMeanOccupancy(binBCbegin, binBCend, occFDDCUnfm80);
            genTwmoFDD(weightMeanOccFDDAUnfm80,
                       weightMeanOccFDDCUnfm80);
            fillQAInfo<kWeightMean, kRobustT0V0Prim, kOccFDDAUnfm80>(weightMeanOccFDDAUnfm80, meanOccRobustT0V0PrimUnfm80);
            fillQAInfo<kWeightMean, kRobustT0V0Prim, kOccFDDCUnfm80>(weightMeanOccFDDCUnfm80, meanOccRobustT0V0PrimUnfm80);

            fillQAInfo<kWeightMean, kWeightRobustT0V0Prim, kOccFDDAUnfm80>(weightMeanOccFDDAUnfm80, weightMeanOccRobustT0V0PrimUnfm80);
            fillQAInfo<kWeightMean, kWeightRobustT0V0Prim, kOccFDDCUnfm80>(weightMeanOccFDDCUnfm80, weightMeanOccRobustT0V0PrimUnfm80);
          }
        }

        if constexpr (processMode == kProcessFullOccTableProducer || processMode == kProcessOnlyOccNtrackDet) {
          if constexpr (meanTableMode == fillMeanOccTable) {
            meanOccNTrackITSUnfm80 = getMeanOccupancy(binBCbegin, binBCend, occNTrackITSUnfm80);
            meanOccNTrackTPCUnfm80 = getMeanOccupancy(binBCbegin, binBCend, occNTrackTPCUnfm80);
            meanOccNTrackTRDUnfm80 = getMeanOccupancy(binBCbegin, binBCend, occNTrackTRDUnfm80);
            meanOccNTrackTOFUnfm80 = getMeanOccupancy(binBCbegin, binBCend, occNTrackTOFUnfm80);
            meanOccNTrackSizeUnfm80 = getMeanOccupancy(binBCbegin, binBCend, occNTrackSizeUnfm80);
            meanOccNTrackTPCAUnfm80 = getMeanOccupancy(binBCbegin, binBCend, occNTrackTPCAUnfm80);
            meanOccNTrackTPCCUnfm80 = getMeanOccupancy(binBCbegin, binBCend, occNTrackTPCCUnfm80);
            meanOccNTrackITSTPCUnfm80 = getMeanOccupancy(binBCbegin, binBCend, occNTrackITSTPCUnfm80);
            meanOccNTrackITSTPCAUnfm80 = getMeanOccupancy(binBCbegin, binBCend, occNTrackITSTPCAUnfm80);
            meanOccNTrackITSTPCCUnfm80 = getMeanOccupancy(binBCbegin, binBCend, occNTrackITSTPCCUnfm80);
            genTmoNTrackDet(meanOccNTrackITSUnfm80,
                            meanOccNTrackTPCUnfm80,
                            meanOccNTrackTRDUnfm80,
                            meanOccNTrackTOFUnfm80,
                            meanOccNTrackSizeUnfm80,
                            meanOccNTrackTPCAUnfm80,
                            meanOccNTrackTPCCUnfm80,
                            meanOccNTrackITSTPCUnfm80,
                            meanOccNTrackITSTPCAUnfm80,
                            meanOccNTrackITSTPCCUnfm80);
            fillQAInfo<kMean, kRobustT0V0Prim, kOccNTrackITSUnfm80>(meanOccNTrackITSUnfm80, meanOccRobustT0V0PrimUnfm80);
            fillQAInfo<kMean, kRobustT0V0Prim, kOccNTrackTPCUnfm80>(meanOccNTrackTPCUnfm80, meanOccRobustT0V0PrimUnfm80);
            fillQAInfo<kMean, kRobustT0V0Prim, kOccNTrackTRDUnfm80>(meanOccNTrackTRDUnfm80, meanOccRobustT0V0PrimUnfm80);
            fillQAInfo<kMean, kRobustT0V0Prim, kOccNTrackTOFUnfm80>(meanOccNTrackTOFUnfm80, meanOccRobustT0V0PrimUnfm80);
            fillQAInfo<kMean, kRobustT0V0Prim, kOccNTrackSizeUnfm80>(meanOccNTrackSizeUnfm80, meanOccRobustT0V0PrimUnfm80);
            fillQAInfo<kMean, kRobustT0V0Prim, kOccNTrackTPCAUnfm80>(meanOccNTrackTPCAUnfm80, meanOccRobustT0V0PrimUnfm80);
            fillQAInfo<kMean, kRobustT0V0Prim, kOccNTrackTPCCUnfm80>(meanOccNTrackTPCCUnfm80, meanOccRobustT0V0PrimUnfm80);
            fillQAInfo<kMean, kRobustT0V0Prim, kOccNTrackITSTPCUnfm80>(meanOccNTrackITSTPCUnfm80, meanOccRobustT0V0PrimUnfm80);
            fillQAInfo<kMean, kRobustT0V0Prim, kOccNTrackITSTPCAUnfm80>(meanOccNTrackITSTPCAUnfm80, meanOccRobustT0V0PrimUnfm80);
            fillQAInfo<kMean, kRobustT0V0Prim, kOccNTrackITSTPCCUnfm80>(meanOccNTrackITSTPCCUnfm80, meanOccRobustT0V0PrimUnfm80);
          }
          if constexpr (weightMeanTableMode == fillWeightMeanOccTable) {
            weightMeanOccNTrackITSUnfm80 = getWeightedMeanOccupancy(binBCbegin, binBCend, occNTrackITSUnfm80);
            weightMeanOccNTrackTPCUnfm80 = getWeightedMeanOccupancy(binBCbegin, binBCend, occNTrackTPCUnfm80);
            weightMeanOccNTrackTRDUnfm80 = getWeightedMeanOccupancy(binBCbegin, binBCend, occNTrackTRDUnfm80);
            weightMeanOccNTrackTOFUnfm80 = getWeightedMeanOccupancy(binBCbegin, binBCend, occNTrackTOFUnfm80);
            weightMeanOccNTrackSizeUnfm80 = getWeightedMeanOccupancy(binBCbegin, binBCend, occNTrackSizeUnfm80);
            weightMeanOccNTrackTPCAUnfm80 = getWeightedMeanOccupancy(binBCbegin, binBCend, occNTrackTPCAUnfm80);
            weightMeanOccNTrackTPCCUnfm80 = getWeightedMeanOccupancy(binBCbegin, binBCend, occNTrackTPCCUnfm80);
            weightMeanOccNTrackITSTPCUnfm80 = getWeightedMeanOccupancy(binBCbegin, binBCend, occNTrackITSTPCUnfm80);
            weightMeanOccNTrackITSTPCAUnfm80 = getWeightedMeanOccupancy(binBCbegin, binBCend, occNTrackITSTPCAUnfm80);
            weightMeanOccNTrackITSTPCCUnfm80 = getWeightedMeanOccupancy(binBCbegin, binBCend, occNTrackITSTPCCUnfm80);

            genTwmoNTrackDet(weightMeanOccNTrackITSUnfm80,
                             weightMeanOccNTrackTPCUnfm80,
                             weightMeanOccNTrackTRDUnfm80,
                             weightMeanOccNTrackTOFUnfm80,
                             weightMeanOccNTrackSizeUnfm80,
                             weightMeanOccNTrackTPCAUnfm80,
                             weightMeanOccNTrackTPCCUnfm80,
                             weightMeanOccNTrackITSTPCUnfm80,
                             weightMeanOccNTrackITSTPCAUnfm80,
                             weightMeanOccNTrackITSTPCCUnfm80);

            fillQAInfo<kWeightMean, kRobustT0V0Prim, kOccNTrackITSUnfm80>(weightMeanOccNTrackITSUnfm80, meanOccRobustT0V0PrimUnfm80);
            fillQAInfo<kWeightMean, kRobustT0V0Prim, kOccNTrackTPCUnfm80>(weightMeanOccNTrackTPCUnfm80, meanOccRobustT0V0PrimUnfm80);
            fillQAInfo<kWeightMean, kRobustT0V0Prim, kOccNTrackTRDUnfm80>(weightMeanOccNTrackTRDUnfm80, meanOccRobustT0V0PrimUnfm80);
            fillQAInfo<kWeightMean, kRobustT0V0Prim, kOccNTrackTOFUnfm80>(weightMeanOccNTrackTOFUnfm80, meanOccRobustT0V0PrimUnfm80);
            fillQAInfo<kWeightMean, kRobustT0V0Prim, kOccNTrackSizeUnfm80>(weightMeanOccNTrackSizeUnfm80, meanOccRobustT0V0PrimUnfm80);
            fillQAInfo<kWeightMean, kRobustT0V0Prim, kOccNTrackTPCAUnfm80>(weightMeanOccNTrackTPCAUnfm80, meanOccRobustT0V0PrimUnfm80);
            fillQAInfo<kWeightMean, kRobustT0V0Prim, kOccNTrackTPCCUnfm80>(weightMeanOccNTrackTPCCUnfm80, meanOccRobustT0V0PrimUnfm80);
            fillQAInfo<kWeightMean, kRobustT0V0Prim, kOccNTrackITSTPCUnfm80>(weightMeanOccNTrackITSTPCUnfm80, meanOccRobustT0V0PrimUnfm80);
            fillQAInfo<kWeightMean, kRobustT0V0Prim, kOccNTrackITSTPCAUnfm80>(weightMeanOccNTrackITSTPCAUnfm80, meanOccRobustT0V0PrimUnfm80);
            fillQAInfo<kWeightMean, kRobustT0V0Prim, kOccNTrackITSTPCCUnfm80>(weightMeanOccNTrackITSTPCCUnfm80, meanOccRobustT0V0PrimUnfm80);

            fillQAInfo<kWeightMean, kWeightRobustT0V0Prim, kOccNTrackITSUnfm80>(weightMeanOccNTrackITSUnfm80, weightMeanOccRobustT0V0PrimUnfm80);
            fillQAInfo<kWeightMean, kWeightRobustT0V0Prim, kOccNTrackTPCUnfm80>(weightMeanOccNTrackTPCUnfm80, weightMeanOccRobustT0V0PrimUnfm80);
            fillQAInfo<kWeightMean, kWeightRobustT0V0Prim, kOccNTrackTRDUnfm80>(weightMeanOccNTrackTRDUnfm80, weightMeanOccRobustT0V0PrimUnfm80);
            fillQAInfo<kWeightMean, kWeightRobustT0V0Prim, kOccNTrackTOFUnfm80>(weightMeanOccNTrackTOFUnfm80, weightMeanOccRobustT0V0PrimUnfm80);
            fillQAInfo<kWeightMean, kWeightRobustT0V0Prim, kOccNTrackSizeUnfm80>(weightMeanOccNTrackSizeUnfm80, weightMeanOccRobustT0V0PrimUnfm80);
            fillQAInfo<kWeightMean, kWeightRobustT0V0Prim, kOccNTrackTPCAUnfm80>(weightMeanOccNTrackTPCAUnfm80, weightMeanOccRobustT0V0PrimUnfm80);
            fillQAInfo<kWeightMean, kWeightRobustT0V0Prim, kOccNTrackTPCCUnfm80>(weightMeanOccNTrackTPCCUnfm80, weightMeanOccRobustT0V0PrimUnfm80);
            fillQAInfo<kWeightMean, kWeightRobustT0V0Prim, kOccNTrackITSTPCUnfm80>(weightMeanOccNTrackITSTPCUnfm80, weightMeanOccRobustT0V0PrimUnfm80);
            fillQAInfo<kWeightMean, kWeightRobustT0V0Prim, kOccNTrackITSTPCAUnfm80>(weightMeanOccNTrackITSTPCAUnfm80, weightMeanOccRobustT0V0PrimUnfm80);
            fillQAInfo<kWeightMean, kWeightRobustT0V0Prim, kOccNTrackITSTPCCUnfm80>(weightMeanOccNTrackITSTPCCUnfm80, weightMeanOccRobustT0V0PrimUnfm80);
          }
        }

        if constexpr (processMode == kProcessFullOccTableProducer || processMode == kProcessOnlyOccMultExtra) {
          if constexpr (meanTableMode == fillMeanOccTable) {
            meanOccMultNTracksHasITSUnfm80 = getMeanOccupancy(binBCbegin, binBCend, occMultNTracksHasITSUnfm80);
            meanOccMultNTracksHasTPCUnfm80 = getMeanOccupancy(binBCbegin, binBCend, occMultNTracksHasTPCUnfm80);
            meanOccMultNTracksHasTOFUnfm80 = getMeanOccupancy(binBCbegin, binBCend, occMultNTracksHasTOFUnfm80);
            meanOccMultNTracksHasTRDUnfm80 = getMeanOccupancy(binBCbegin, binBCend, occMultNTracksHasTRDUnfm80);
            meanOccMultNTracksITSOnlyUnfm80 = getMeanOccupancy(binBCbegin, binBCend, occMultNTracksITSOnlyUnfm80);
            meanOccMultNTracksTPCOnlyUnfm80 = getMeanOccupancy(binBCbegin, binBCend, occMultNTracksTPCOnlyUnfm80);
            meanOccMultNTracksITSTPCUnfm80 = getMeanOccupancy(binBCbegin, binBCend, occMultNTracksITSTPCUnfm80);
            meanOccMultAllTracksTPCOnlyUnfm80 = getMeanOccupancy(binBCbegin, binBCend, occMultAllTracksTPCOnlyUnfm80);
            genTmoMultExtra(meanOccMultNTracksHasITSUnfm80,
                            meanOccMultNTracksHasTPCUnfm80,
                            meanOccMultNTracksHasTOFUnfm80,
                            meanOccMultNTracksHasTRDUnfm80,
                            meanOccMultNTracksITSOnlyUnfm80,
                            meanOccMultNTracksTPCOnlyUnfm80,
                            meanOccMultNTracksITSTPCUnfm80,
                            meanOccMultAllTracksTPCOnlyUnfm80);
            fillQAInfo<kMean, kRobustT0V0Prim, kOccMultNTracksHasITSUnfm80>(meanOccMultNTracksHasITSUnfm80, meanOccRobustT0V0PrimUnfm80);
            fillQAInfo<kMean, kRobustT0V0Prim, kOccMultNTracksHasTPCUnfm80>(meanOccMultNTracksHasTPCUnfm80, meanOccRobustT0V0PrimUnfm80);
            fillQAInfo<kMean, kRobustT0V0Prim, kOccMultNTracksHasTOFUnfm80>(meanOccMultNTracksHasTOFUnfm80, meanOccRobustT0V0PrimUnfm80);
            fillQAInfo<kMean, kRobustT0V0Prim, kOccMultNTracksHasTRDUnfm80>(meanOccMultNTracksHasTRDUnfm80, meanOccRobustT0V0PrimUnfm80);
            fillQAInfo<kMean, kRobustT0V0Prim, kOccMultNTracksITSOnlyUnfm80>(meanOccMultNTracksITSOnlyUnfm80, meanOccRobustT0V0PrimUnfm80);
            fillQAInfo<kMean, kRobustT0V0Prim, kOccMultNTracksTPCOnlyUnfm80>(meanOccMultNTracksTPCOnlyUnfm80, meanOccRobustT0V0PrimUnfm80);
            fillQAInfo<kMean, kRobustT0V0Prim, kOccMultNTracksITSTPCUnfm80>(meanOccMultNTracksITSTPCUnfm80, meanOccRobustT0V0PrimUnfm80);
            fillQAInfo<kMean, kRobustT0V0Prim, kOccMultAllTracksTPCOnlyUnfm80>(meanOccMultAllTracksTPCOnlyUnfm80, meanOccRobustT0V0PrimUnfm80);
          }
          if constexpr (weightMeanTableMode == fillWeightMeanOccTable) {
            weightMeanOccMultNTracksHasITSUnfm80 = getWeightedMeanOccupancy(binBCbegin, binBCend, occMultNTracksHasITSUnfm80);
            weightMeanOccMultNTracksHasTPCUnfm80 = getWeightedMeanOccupancy(binBCbegin, binBCend, occMultNTracksHasTPCUnfm80);
            weightMeanOccMultNTracksHasTOFUnfm80 = getWeightedMeanOccupancy(binBCbegin, binBCend, occMultNTracksHasTOFUnfm80);
            weightMeanOccMultNTracksHasTRDUnfm80 = getWeightedMeanOccupancy(binBCbegin, binBCend, occMultNTracksHasTRDUnfm80);
            weightMeanOccMultNTracksITSOnlyUnfm80 = getWeightedMeanOccupancy(binBCbegin, binBCend, occMultNTracksITSOnlyUnfm80);
            weightMeanOccMultNTracksTPCOnlyUnfm80 = getWeightedMeanOccupancy(binBCbegin, binBCend, occMultNTracksTPCOnlyUnfm80);
            weightMeanOccMultNTracksITSTPCUnfm80 = getWeightedMeanOccupancy(binBCbegin, binBCend, occMultNTracksITSTPCUnfm80);
            weightMeanOccMultAllTracksTPCOnlyUnfm80 = getWeightedMeanOccupancy(binBCbegin, binBCend, occMultAllTracksTPCOnlyUnfm80);

            genTwmoMultExtra(weightMeanOccMultNTracksHasITSUnfm80,
                             weightMeanOccMultNTracksHasTPCUnfm80,
                             weightMeanOccMultNTracksHasTOFUnfm80,
                             weightMeanOccMultNTracksHasTRDUnfm80,
                             weightMeanOccMultNTracksITSOnlyUnfm80,
                             weightMeanOccMultNTracksTPCOnlyUnfm80,
                             weightMeanOccMultNTracksITSTPCUnfm80,
                             weightMeanOccMultAllTracksTPCOnlyUnfm80);

            fillQAInfo<kWeightMean, kRobustT0V0Prim, kOccMultNTracksHasITSUnfm80>(weightMeanOccMultNTracksHasITSUnfm80, meanOccRobustT0V0PrimUnfm80);
            fillQAInfo<kWeightMean, kRobustT0V0Prim, kOccMultNTracksHasTPCUnfm80>(weightMeanOccMultNTracksHasTPCUnfm80, meanOccRobustT0V0PrimUnfm80);
            fillQAInfo<kWeightMean, kRobustT0V0Prim, kOccMultNTracksHasTOFUnfm80>(weightMeanOccMultNTracksHasTOFUnfm80, meanOccRobustT0V0PrimUnfm80);
            fillQAInfo<kWeightMean, kRobustT0V0Prim, kOccMultNTracksHasTRDUnfm80>(weightMeanOccMultNTracksHasTRDUnfm80, meanOccRobustT0V0PrimUnfm80);
            fillQAInfo<kWeightMean, kRobustT0V0Prim, kOccMultNTracksITSOnlyUnfm80>(weightMeanOccMultNTracksITSOnlyUnfm80, meanOccRobustT0V0PrimUnfm80);
            fillQAInfo<kWeightMean, kRobustT0V0Prim, kOccMultNTracksTPCOnlyUnfm80>(weightMeanOccMultNTracksTPCOnlyUnfm80, meanOccRobustT0V0PrimUnfm80);
            fillQAInfo<kWeightMean, kRobustT0V0Prim, kOccMultNTracksITSTPCUnfm80>(weightMeanOccMultNTracksITSTPCUnfm80, meanOccRobustT0V0PrimUnfm80);
            fillQAInfo<kWeightMean, kRobustT0V0Prim, kOccMultAllTracksTPCOnlyUnfm80>(weightMeanOccMultAllTracksTPCOnlyUnfm80, meanOccRobustT0V0PrimUnfm80);

            fillQAInfo<kWeightMean, kWeightRobustT0V0Prim, kOccMultNTracksHasITSUnfm80>(weightMeanOccMultNTracksHasITSUnfm80, weightMeanOccRobustT0V0PrimUnfm80);
            fillQAInfo<kWeightMean, kWeightRobustT0V0Prim, kOccMultNTracksHasTPCUnfm80>(weightMeanOccMultNTracksHasTPCUnfm80, weightMeanOccRobustT0V0PrimUnfm80);
            fillQAInfo<kWeightMean, kWeightRobustT0V0Prim, kOccMultNTracksHasTOFUnfm80>(weightMeanOccMultNTracksHasTOFUnfm80, weightMeanOccRobustT0V0PrimUnfm80);
            fillQAInfo<kWeightMean, kWeightRobustT0V0Prim, kOccMultNTracksHasTRDUnfm80>(weightMeanOccMultNTracksHasTRDUnfm80, weightMeanOccRobustT0V0PrimUnfm80);
            fillQAInfo<kWeightMean, kWeightRobustT0V0Prim, kOccMultNTracksITSOnlyUnfm80>(weightMeanOccMultNTracksITSOnlyUnfm80, weightMeanOccRobustT0V0PrimUnfm80);
            fillQAInfo<kWeightMean, kWeightRobustT0V0Prim, kOccMultNTracksTPCOnlyUnfm80>(weightMeanOccMultNTracksTPCOnlyUnfm80, weightMeanOccRobustT0V0PrimUnfm80);
            fillQAInfo<kWeightMean, kWeightRobustT0V0Prim, kOccMultNTracksITSTPCUnfm80>(weightMeanOccMultNTracksITSTPCUnfm80, weightMeanOccRobustT0V0PrimUnfm80);
            fillQAInfo<kWeightMean, kWeightRobustT0V0Prim, kOccMultAllTracksTPCOnlyUnfm80>(weightMeanOccMultAllTracksTPCOnlyUnfm80, weightMeanOccRobustT0V0PrimUnfm80);
          }
        }

        if constexpr (processMode == kProcessFullOccTableProducer || processMode == kProcessOnlyRobustT0V0Prim) {
          if constexpr (meanTableMode == fillMeanOccTable) {
            meanOccRobustT0V0PrimUnfm80 = getMeanOccupancy(binBCbegin, binBCend, occRobustT0V0PrimUnfm80);
            genTmoRT0V0Prim(meanOccRobustT0V0PrimUnfm80);
            fillQAInfo<kMean, kRobustT0V0Prim, kOccRobustT0V0PrimUnfm80>(meanOccRobustT0V0PrimUnfm80, meanOccRobustT0V0PrimUnfm80);
          }
          if constexpr (weightMeanTableMode == fillWeightMeanOccTable) {
            weightMeanOccRobustT0V0PrimUnfm80 = getWeightedMeanOccupancy(binBCbegin, binBCend, occRobustT0V0PrimUnfm80);
            genTwmoRT0V0Prim(weightMeanOccRobustT0V0PrimUnfm80);
            fillQAInfo<kWeightMean, kRobustT0V0Prim, kOccRobustT0V0PrimUnfm80>(weightMeanOccRobustT0V0PrimUnfm80, meanOccRobustT0V0PrimUnfm80);
            fillQAInfo<kWeightMean, kWeightRobustT0V0Prim, kOccRobustT0V0PrimUnfm80>(weightMeanOccRobustT0V0PrimUnfm80, weightMeanOccRobustT0V0PrimUnfm80);
          }
        }

        if constexpr (processMode == kProcessFullOccTableProducer || processMode == kProcessOnlyRobustFDDT0V0Prim) {
          if constexpr (meanTableMode == fillMeanOccTable) {
            meanOccRobustFDDT0V0PrimUnfm80 = getMeanOccupancy(binBCbegin, binBCend, occRobustFDDT0V0PrimUnfm80);
            genTmoRFDDT0V0Prim(meanOccRobustFDDT0V0PrimUnfm80);
            fillQAInfo<kMean, kRobustT0V0Prim, kOccRobustFDDT0V0PrimUnfm80>(meanOccRobustFDDT0V0PrimUnfm80, meanOccRobustT0V0PrimUnfm80);
          }
          if constexpr (weightMeanTableMode == fillWeightMeanOccTable) {
            weightMeanOccRobustFDDT0V0PrimUnfm80 = getWeightedMeanOccupancy(binBCbegin, binBCend, occRobustFDDT0V0PrimUnfm80);
            genTwmoRFDDT0V0Pri(weightMeanOccRobustFDDT0V0PrimUnfm80);
            fillQAInfo<kWeightMean, kRobustT0V0Prim, kOccRobustFDDT0V0PrimUnfm80>(weightMeanOccRobustFDDT0V0PrimUnfm80, meanOccRobustT0V0PrimUnfm80);
            fillQAInfo<kWeightMean, kWeightRobustT0V0Prim, kOccRobustFDDT0V0PrimUnfm80>(weightMeanOccRobustFDDT0V0PrimUnfm80, weightMeanOccRobustT0V0PrimUnfm80);
          }
        }

        if constexpr (processMode == kProcessFullOccTableProducer || processMode == kProcessOnlyRobustNtrackDet) {
          if constexpr (meanTableMode == fillMeanOccTable) {
            meanOccRobustNtrackDetUnfm80 = getMeanOccupancy(binBCbegin, binBCend, occRobustNtrackDetUnfm80);
            genTmoRNtrackDet(meanOccRobustNtrackDetUnfm80);
            fillQAInfo<kMean, kRobustT0V0Prim, kOccRobustNtrackDetUnfm80>(meanOccRobustNtrackDetUnfm80, meanOccRobustT0V0PrimUnfm80);
          }
          if constexpr (weightMeanTableMode == fillWeightMeanOccTable) {
            weightMeanOccRobustNtrackDetUnfm80 = getWeightedMeanOccupancy(binBCbegin, binBCend, occRobustNtrackDetUnfm80);
            genTwmoRNtrackDet(weightMeanOccRobustNtrackDetUnfm80);
            fillQAInfo<kWeightMean, kRobustT0V0Prim, kOccRobustNtrackDetUnfm80>(weightMeanOccRobustNtrackDetUnfm80, meanOccRobustT0V0PrimUnfm80);
            fillQAInfo<kWeightMean, kWeightRobustT0V0Prim, kOccRobustNtrackDetUnfm80>(weightMeanOccRobustNtrackDetUnfm80, weightMeanOccRobustT0V0PrimUnfm80);
          }
        }

        if constexpr (processMode == kProcessFullOccTableProducer || processMode == kProcessOnlyRobustMultExtra) {
          if constexpr (meanTableMode == fillMeanOccTable) {
            meanOccRobustMultTableUnfm80 = getMeanOccupancy(binBCbegin, binBCend, occRobustMultTableUnfm80);
            genTmoRMultExtra(meanOccRobustMultTableUnfm80);
            fillQAInfo<kMean, kRobustT0V0Prim, kOccRobustMultTableUnfm80>(meanOccRobustMultTableUnfm80, meanOccRobustT0V0PrimUnfm80);
          }
          if constexpr (weightMeanTableMode == fillWeightMeanOccTable) {
            weightMeanOccRobustMultTableUnfm80 = getWeightedMeanOccupancy(binBCbegin, binBCend, occRobustMultTableUnfm80);
            genTwmoRMultExtra(weightMeanOccRobustMultTableUnfm80);
            fillQAInfo<kWeightMean, kRobustT0V0Prim, kOccRobustMultTableUnfm80>(weightMeanOccRobustMultTableUnfm80, meanOccRobustT0V0PrimUnfm80);
            fillQAInfo<kWeightMean, kWeightRobustT0V0Prim, kOccRobustMultTableUnfm80>(weightMeanOccRobustMultTableUnfm80, weightMeanOccRobustT0V0PrimUnfm80);
          }
        }
      } // end of trackQA loop

      // build the IndexTables here
      if (executeInThisBlock) {
        if (buildPointerTrackQAToTMOTable) {
          // create pointer table from trackQA to TrackMeanOcc
          sortVectorOfArray(trackQAGIListforTMOList, 0); // sort the list //Its easy to search in a sorted list
          checkUniqueness(trackQAGIListforTMOList, 0);   // check the uniqueness of track.globalIndex()

          int currentIDXforCheck = 0;
          int listSize = trackQAGIListforTMOList.size();
          for (const auto& trackQA : tracksQA) {
            while (trackQA.globalIndex() > trackQAGIListforTMOList[currentIDXforCheck][0]) {
              currentIDXforCheck++; // increment the currentIDXforCheck for missing or invalid cases e.g. value = -1;
              if (currentIDXforCheck >= listSize) {
                break;
              }
            }
            if (trackQA.globalIndex() == trackQAGIListforTMOList[currentIDXforCheck][0]) {
              genTrackQAToTmo(trackQAGIListforTMOList[currentIDXforCheck][1]);
            } else {
              genTrackQAToTmo(-1); // put a dummy index when track is not found in trackQA
            }
          }
        }
        if (buildPointerTMOToTrackQATable) {
          // create pointer table from TrackMeanOcc to trackQA
          sortVectorOfArray(trackQAGIListforTMOList, 1); // sort the list //Its easy to search in a sorted list
          checkUniqueness(trackQAGIListforTMOList, 1);   // check the uniqueness of track.globalIndex()

          int currentIDXforCheck = 0;
          int listSize = trackQAGIListforTMOList.size();
          for (int iCounter = 0; iCounter <= trackTMOcounter; iCounter++) {
            while (iCounter > trackQAGIListforTMOList[currentIDXforCheck][1]) {
              currentIDXforCheck++; // increment the currentIDXforCheck for missing or invalid cases e.g. value = -1;
              if (currentIDXforCheck >= listSize) {
                break;
              }
            }
            if (iCounter == trackQAGIListforTMOList[currentIDXforCheck][1]) {
              genTmoToTrackQA(trackQAGIListforTMOList[currentIDXforCheck][0]);
            } else {
              genTmoToTrackQA(-1); // put a dummy index when track is not found in trackQA
            }
          }
        }
      } // end of executeInThisBlock
    } // end of else block of constexpr
  }

  void checkAllProcessFunctionStatus(std::vector<bool> const& processStatusVector, bool& singleProcessOn)
  {
    int nProcessOn = 0;
    const uint size = processStatusVector.size();
    for (uint i = 0; i < size; i++) {
      if (processStatusVector[i]) {
        nProcessOn++;
      }
    }
    if (nProcessOn > 1) {
      singleProcessOn = false;
    } // check nProcess
  }

  //_________________________________Process Functions start from here______________________________________________________________________________________
  void processNothing(aod::Collisions const&)
  {
    occupancyQA.fill(HIST("h_DFcount_Lvl0"), kProcessNothing);
    return;
    occupancyQA.fill(HIST("h_DFcount_Lvl1"), kProcessNothing);
  }
  PROCESS_SWITCH(TrackMeanOccTableProducer, processNothing, "process Nothing From Track Mean Occ Table Producer", true);

  void processOnlyOccPrim(soa::Join<aod::BCsWithTimestamps, aod::OccIndexTable> const& BCs, aod::Collisions const& collisions, aod::Tracks const& tracks, MyTracksQA const& tracksQA, aod::ORT0V0Prim const& occsRobustT0V0Prim, aod::OccsPrim const& occs)
  {
    processStatus[kProcessOnlyOccPrim] = true;
    bool singleProcessOn = true;
    checkAllProcessFunctionStatus(processStatus, singleProcessOn);
    if (singleProcessOn) {
      processInThisBlock[kProcessOnlyOccPrim] = true;
    }
    occupancyQA.fill(HIST("h_DFcount_Lvl0"), kProcessOnlyOccPrim);
    if (collisions.size() == 0) {
      return;
    }
    if (!buildOnlyOccsPrim && !buildFullOccTableProducer) {
      LOG(error) << "DEBUG :: ERROR ERROR ERROR :: buildOnlyOccsPrim == false";
      return;
    }
    if (occs.size() == 0) {
      LOG(error) << "DEBUG :: ERROR ERROR ERROR :: OccsPrim.size() == 0 :: Check \"occupancy-table-producer\" for \"buildOnlyOccsPrim == true\" & \"processOnlyOccPrimUnfm == true\"";
      return;
    }
    executeTrackOccProducerProcessing<kProcessOnlyOccPrim, checkTableMode, checkTableMode, checkQAMode>(BCs, collisions, tracks, tracksQA, occsRobustT0V0Prim, occs, processInThisBlock[kProcessOnlyOccPrim]);
    occupancyQA.fill(HIST("h_DFcount_Lvl1"), kProcessOnlyOccPrim);
  }
  PROCESS_SWITCH(TrackMeanOccTableProducer, processOnlyOccPrim, "processOnlyOccPrim", false);

  void processOnlyOccT0V0(soa::Join<aod::BCsWithTimestamps, aod::OccIndexTable> const& BCs, aod::Collisions const& collisions, aod::Tracks const& tracks, MyTracksQA const& tracksQA, aod::ORT0V0Prim const& occsRobustT0V0Prim, aod::OccsT0V0 const& occs)
  {
    processStatus[kProcessOnlyOccT0V0] = true;
    bool singleProcessOn = true;
    checkAllProcessFunctionStatus(processStatus, singleProcessOn);
    if (singleProcessOn) {
      processInThisBlock[kProcessOnlyOccT0V0] = true;
    }
    occupancyQA.fill(HIST("h_DFcount_Lvl0"), kProcessOnlyOccT0V0);
    if (collisions.size() == 0) {
      return;
    }
    if (!buildOnlyOccsT0V0 && !buildFullOccTableProducer) {
      LOG(error) << "DEBUG :: ERROR ERROR ERROR :: buildOnlyOccsT0V0 == false";
      return;
    }
    if (occs.size() == 0) {
      LOG(error) << "DEBUG :: ERROR ERROR ERROR :: OccsT0V0.size() == 0 :: Check \"occupancy-table-producer\" for \"buildOnlyOccsT0V0Prim == true\" & \"processOnlyOccT0V0PrimUnfm == true\"";
      return;
    }
    executeTrackOccProducerProcessing<kProcessOnlyOccT0V0, checkTableMode, checkTableMode, checkQAMode>(BCs, collisions, tracks, tracksQA, occsRobustT0V0Prim, occs, processInThisBlock[kProcessOnlyOccT0V0]);
    occupancyQA.fill(HIST("h_DFcount_Lvl1"), kProcessOnlyOccT0V0);
  }
  PROCESS_SWITCH(TrackMeanOccTableProducer, processOnlyOccT0V0, "processOnlyOccT0V0", false);

  void processOnlyOccFDD(soa::Join<aod::BCsWithTimestamps, aod::OccIndexTable> const& BCs, aod::Collisions const& collisions, aod::Tracks const& tracks, MyTracksQA const& tracksQA, aod::ORT0V0Prim const& occsRobustT0V0Prim, aod::OccsFDD const& occs)
  {
    processStatus[kProcessOnlyOccFDD] = true;
    bool singleProcessOn = true;
    checkAllProcessFunctionStatus(processStatus, singleProcessOn);
    if (singleProcessOn) {
      processInThisBlock[kProcessOnlyOccFDD] = true;
    }
    occupancyQA.fill(HIST("h_DFcount_Lvl0"), kProcessOnlyOccFDD);
    if (collisions.size() == 0) {
      return;
    }
    if (!buildOnlyOccsFDD && !buildFullOccTableProducer) {
      LOG(error) << "DEBUG :: ERROR ERROR ERROR :: buildOnlyOccsFDD == false";
      return;
    }
    if (occs.size() == 0) {
      LOG(error) << "DEBUG :: ERROR ERROR ERROR :: OccsFDD.size() == 0 :: Check \"occupancy-table-producer\" for \"buildOnlyOccsFDDT0V0Prim == true\" & \"processOnlyOccFDDT0V0PrimUnfm == true\"";
      return;
    }
    executeTrackOccProducerProcessing<kProcessOnlyOccFDD, checkTableMode, checkTableMode, checkQAMode>(BCs, collisions, tracks, tracksQA, occsRobustT0V0Prim, occs, processInThisBlock[kProcessOnlyOccFDD]);
    occupancyQA.fill(HIST("h_DFcount_Lvl1"), kProcessOnlyOccFDD);
  }
  PROCESS_SWITCH(TrackMeanOccTableProducer, processOnlyOccFDD, "processOnlyOccFDD", false);

  void processOnlyOccNtrackDet(soa::Join<aod::BCsWithTimestamps, aod::OccIndexTable> const& BCs, aod::Collisions const& collisions, aod::Tracks const& tracks, MyTracksQA const& tracksQA, aod::ORT0V0Prim const& occsRobustT0V0Prim, aod::OccsNTrackDet const& occs)
  {
    processStatus[kProcessOnlyOccNtrackDet] = true;
    bool singleProcessOn = true;
    checkAllProcessFunctionStatus(processStatus, singleProcessOn);
    if (singleProcessOn) {
      processInThisBlock[kProcessOnlyOccNtrackDet] = true;
    }
    occupancyQA.fill(HIST("h_DFcount_Lvl0"), kProcessOnlyOccNtrackDet);
    if (collisions.size() == 0) {
      return;
    }
    if (!buildOnlyOccsNtrackDet && !buildFullOccTableProducer) {
      LOG(error) << "DEBUG :: ERROR ERROR ERROR :: buildOnlyOccsNtrackDet == false";
      return;
    }
    if (occs.size() == 0) {
      LOG(error) << "DEBUG :: ERROR ERROR ERROR :: OccsNtrackDet.size() == 0 :: Check \"occupancy-table-producer\" for \"buildOnlyOccsNtrackDet == true\" & \"processOnlyOccNtrackDet == true\"";
      return;
    }
    executeTrackOccProducerProcessing<kProcessOnlyOccNtrackDet, checkTableMode, checkTableMode, checkQAMode>(BCs, collisions, tracks, tracksQA, occsRobustT0V0Prim, occs, processInThisBlock[kProcessOnlyOccNtrackDet]);
    occupancyQA.fill(HIST("h_DFcount_Lvl1"), kProcessOnlyOccNtrackDet);
  }
  PROCESS_SWITCH(TrackMeanOccTableProducer, processOnlyOccNtrackDet, "processOnlyOccNtrackDet", false);

  void processOnlyOccMultExtra(soa::Join<aod::BCsWithTimestamps, aod::OccIndexTable> const& BCs, aod::Collisions const& collisions, aod::Tracks const& tracks, MyTracksQA const& tracksQA, aod::ORT0V0Prim const& occsRobustT0V0Prim, aod::OccsMultExtra const& occs)
  {
    processStatus[kProcessOnlyOccMultExtra] = true;
    bool singleProcessOn = true;
    checkAllProcessFunctionStatus(processStatus, singleProcessOn);
    if (singleProcessOn) {
      processInThisBlock[kProcessOnlyOccMultExtra] = true;
    }
    occupancyQA.fill(HIST("h_DFcount_Lvl0"), kProcessOnlyOccMultExtra);
    if (collisions.size() == 0) {
      return;
    }
    if (!buildOnlyOccsMultExtra && !buildFullOccTableProducer) {
      LOG(error) << "DEBUG :: ERROR ERROR ERROR :: buildOnlyOccsMultExtra == false";
      return;
    }
    if (occs.size() == 0) {
      LOG(error) << "DEBUG :: ERROR ERROR ERROR :: OccsMultExtra.size() == 0 :: Check \"occupancy-table-producer\" for \"buildOnlyOccsMultExtra == true\" & \"processOnlyOccMultExtra == true\"";
      return;
    }
    executeTrackOccProducerProcessing<kProcessOnlyOccMultExtra, checkTableMode, checkTableMode, checkQAMode>(BCs, collisions, tracks, tracksQA, occsRobustT0V0Prim, occs, processInThisBlock[kProcessOnlyOccMultExtra]);
    occupancyQA.fill(HIST("h_DFcount_Lvl1"), kProcessOnlyOccMultExtra);
  }
  PROCESS_SWITCH(TrackMeanOccTableProducer, processOnlyOccMultExtra, "processOnlyOccMultExtra", false);

  void processOnlyRobustT0V0Prim(soa::Join<aod::BCsWithTimestamps, aod::OccIndexTable> const& BCs, aod::Collisions const& collisions, aod::Tracks const& tracks, MyTracksQA const& tracksQA, aod::ORT0V0Prim const& occsRobustT0V0Prim)
  {
    processStatus[kProcessOnlyRobustT0V0Prim] = true;
    bool singleProcessOn = true;
    checkAllProcessFunctionStatus(processStatus, singleProcessOn);
    if (singleProcessOn) {
      processInThisBlock[kProcessOnlyRobustT0V0Prim] = true;
    }
    occupancyQA.fill(HIST("h_DFcount_Lvl0"), kProcessOnlyRobustT0V0Prim);
    if (collisions.size() == 0) {
      return;
    }
    if (!buildOnlyOccsRobustT0V0Prim && !buildFullOccTableProducer) {
      LOG(error) << "DEBUG :: ERROR ERROR ERROR :: buildOnlyOccsRobustT0V0Prim == false";
      return;
    }
    if (occsRobustT0V0Prim.size() == 0) {
      LOG(error) << "DEBUG :: ERROR ERROR ERROR :: OccsRobustT0V0Prim.size() == 0 :: Check \"occupancy-table-producer\" for \"buildOnlyOccsT0V0Prim == true\" & \"processOnlyOccT0V0PrimUnfm == true\"";
      return;
    }
    executeTrackOccProducerProcessing<kProcessOnlyRobustT0V0Prim, checkTableMode, checkTableMode, checkQAMode>(BCs, collisions, tracks, tracksQA, occsRobustT0V0Prim, occsRobustT0V0Prim, processInThisBlock[kProcessOnlyRobustT0V0Prim]);
    occupancyQA.fill(HIST("h_DFcount_Lvl1"), kProcessOnlyRobustT0V0Prim);
  }
  PROCESS_SWITCH(TrackMeanOccTableProducer, processOnlyRobustT0V0Prim, "processOnlyRobustT0V0Prim", false);

  void processOnlyRobustFDDT0V0Prim(soa::Join<aod::BCsWithTimestamps, aod::OccIndexTable> const& BCs, aod::Collisions const& collisions, aod::Tracks const& tracks, MyTracksQA const& tracksQA, aod::ORT0V0Prim const& occsRobustT0V0Prim, aod::ORFDDT0V0Prim const& occs)
  {
    processStatus[kProcessOnlyRobustFDDT0V0Prim] = true;
    bool singleProcessOn = true;
    checkAllProcessFunctionStatus(processStatus, singleProcessOn);
    if (singleProcessOn) {
      processInThisBlock[kProcessOnlyRobustFDDT0V0Prim] = true;
    }
    occupancyQA.fill(HIST("h_DFcount_Lvl0"), kProcessOnlyRobustFDDT0V0Prim);
    if (collisions.size() == 0) {
      return;
    }
    if (!buildOnlyOccsRobustFDDT0V0Prim && !buildFullOccTableProducer) {
      LOG(error) << "DEBUG :: ERROR ERROR ERROR :: buildOnlyOccsRobustFDDT0V0Prim == false";
      return;
    }
    if (occs.size() == 0) {
      LOG(error) << "DEBUG :: ERROR ERROR ERROR :: OccsRobustFDDT0V0Prim.size() == 0 :: Check \"occupancy-table-producer\" for \"buildOnlyOccsFDDT0V0Prim == true\" & \"processOnlyOccFDDT0V0PrimUnfm == true\"";
      return;
    }
    executeTrackOccProducerProcessing<kProcessOnlyRobustFDDT0V0Prim, checkTableMode, checkTableMode, checkQAMode>(BCs, collisions, tracks, tracksQA, occsRobustT0V0Prim, occs, processInThisBlock[kProcessOnlyRobustFDDT0V0Prim]);
    occupancyQA.fill(HIST("h_DFcount_Lvl1"), kProcessOnlyRobustFDDT0V0Prim);
  }
  PROCESS_SWITCH(TrackMeanOccTableProducer, processOnlyRobustFDDT0V0Prim, "processOnlyRobustFDDT0V0Prim", false);

  void processOnlyRobustNtrackDet(soa::Join<aod::BCsWithTimestamps, aod::OccIndexTable> const& BCs, aod::Collisions const& collisions, aod::Tracks const& tracks, MyTracksQA const& tracksQA, aod::ORT0V0Prim const& occsRobustT0V0Prim, aod::ORNtrackDet const& occs)
  {
    processStatus[kProcessOnlyRobustNtrackDet] = true;
    bool singleProcessOn = true;
    checkAllProcessFunctionStatus(processStatus, singleProcessOn);
    if (singleProcessOn) {
      processInThisBlock[kProcessOnlyRobustNtrackDet] = true;
    }
    occupancyQA.fill(HIST("h_DFcount_Lvl0"), kProcessOnlyRobustNtrackDet);
    if (collisions.size() == 0) {
      return;
    }
    if (!buildOnlyOccsRobustNtrackDet && !buildFullOccTableProducer) {
      LOG(error) << "DEBUG :: ERROR ERROR ERROR :: buildOnlyOccsRobustNtrackDet == false";
      return;
    }
    if (occs.size() == 0) {
      LOG(error) << "DEBUG :: ERROR ERROR ERROR :: OccsRobustNtrackDet.size() == 0 :: Check \"occupancy-table-producer\" for \"buildOnlyOccsNtrackDet == true\" & \"processOnlyOccNtrackDet == true\"";
      return;
    }
    executeTrackOccProducerProcessing<kProcessOnlyRobustNtrackDet, checkTableMode, checkTableMode, checkQAMode>(BCs, collisions, tracks, tracksQA, occsRobustT0V0Prim, occs, processInThisBlock[kProcessOnlyRobustNtrackDet]);
    occupancyQA.fill(HIST("h_DFcount_Lvl1"), kProcessOnlyRobustNtrackDet);
  }
  PROCESS_SWITCH(TrackMeanOccTableProducer, processOnlyRobustNtrackDet, "processOnlyRobustNtrackDet", false);

  void processOnlyRobustMultExtra(soa::Join<aod::BCsWithTimestamps, aod::OccIndexTable> const& BCs, aod::Collisions const& collisions, aod::Tracks const& tracks, MyTracksQA const& tracksQA, aod::ORT0V0Prim const& occsRobustT0V0Prim, aod::ORMultExtra const& occs)
  {
    processStatus[kProcessOnlyRobustMultExtra] = true;
    bool singleProcessOn = true;
    checkAllProcessFunctionStatus(processStatus, singleProcessOn);
    if (singleProcessOn) {
      processInThisBlock[kProcessOnlyRobustMultExtra] = true;
    }
    occupancyQA.fill(HIST("h_DFcount_Lvl0"), kProcessOnlyRobustMultExtra);
    if (collisions.size() == 0) {
      return;
    }
    if (!buildOnlyOccsRobustMultExtra && !buildFullOccTableProducer) {
      LOG(error) << "DEBUG :: ERROR ERROR ERROR :: buildOnlyOccsRobustMultExtra == false";
      return;
    }
    if (occs.size() == 0) {
      LOG(error) << "DEBUG :: ERROR ERROR ERROR :: OccsRobustMultExtra.size() == 0 :: Check \"occupancy-table-producer\" for \"buildOnlyOccsMultExtra == true\" & \"processOnlyOccMultExtra == true\"";
      return;
    }
    executeTrackOccProducerProcessing<kProcessOnlyRobustMultExtra, checkTableMode, checkTableMode, checkQAMode>(BCs, collisions, tracks, tracksQA, occsRobustT0V0Prim, occs, processInThisBlock[kProcessOnlyRobustMultExtra]);
    occupancyQA.fill(HIST("h_DFcount_Lvl1"), kProcessOnlyRobustMultExtra);
  }
  PROCESS_SWITCH(TrackMeanOccTableProducer, processOnlyRobustMultExtra, "processOnlyRobustMultExtra", false);

  using JoinedOccTables = soa::Join<aod::OccsPrim, aod::OccsT0V0, aod::OccsFDD, aod::OccsNTrackDet, aod::OccsMultExtra, aod::ORT0V0Prim, aod::ORFDDT0V0Prim, aod::ORNtrackDet, aod::ORMultExtra>;
  void processFullOccTableProduer(soa::Join<aod::BCsWithTimestamps, aod::OccIndexTable> const& BCs, aod::Collisions const& collisions, aod::Tracks const& tracks, MyTracksQA const& tracksQA, aod::ORT0V0Prim const& occsRobustT0V0Prim, JoinedOccTables const& occs)
  {
    processStatus[kProcessFullOccTableProducer] = true;
    bool singleProcessOn = true;
    checkAllProcessFunctionStatus(processStatus, singleProcessOn);
    if (singleProcessOn) {
      processInThisBlock[kProcessFullOccTableProducer] = true;
    }
    if (!singleProcessOn) {
      LOG(error) << "More than one process functions are on in track-mean-occ-table-producer";
      return;
    }
    occupancyQA.fill(HIST("h_DFcount_Lvl0"), kProcessFullOccTableProducer);
    if (collisions.size() == 0) {
      return;
    }
    if (!buildFullOccTableProducer) {
      LOG(error) << "DEBUG :: ERROR ERROR ERROR :: buildFullOccTableProducer == false";
      return;
    }
    if (occs.size() == 0) {
      LOG(error) << "DEBUG :: ERROR ERROR ERROR :: Full Occ Table Join size is 0 :: Check \"occupancy-table-producer\" for \"buildFullOccTableProducer == true\" & \"processFullOccTableProduer == true\"";
      return;
    }
    executeTrackOccProducerProcessing<kProcessFullOccTableProducer, checkTableMode, checkTableMode, checkQAMode>(BCs, collisions, tracks, tracksQA, occsRobustT0V0Prim, occs, processInThisBlock[kProcessFullOccTableProducer]);
    occupancyQA.fill(HIST("h_DFcount_Lvl1"), kProcessFullOccTableProducer);
  } // Process function ends
  PROCESS_SWITCH(TrackMeanOccTableProducer, processFullOccTableProduer, "processFullOccTableProduer", false);
};

struct CreatePointerTables {

  Produces<aod::TrackToTracksQA> genTrackToTracksQA;
  Produces<aod::TrackToTmo> genTrackToTmo;

  void processNothing(aod::Collisions const&)
  {
    return;
  }
  PROCESS_SWITCH(CreatePointerTables, processNothing, "process Nothing", true);

  std::vector<std::array<int64_t, 2>> trackGIForTrackQAIndexList;
  using MyTracksQA = aod::TracksQAVersion;
  void processTrackToTrackQAPointer(aod::Tracks const& tracks, MyTracksQA const& tracksQA)
  {
    trackGIForTrackQAIndexList.clear();
    for (const auto& trackQA : tracksQA) {
      auto const& track = trackQA.template track_as<aod::Tracks>();
      trackGIForTrackQAIndexList.push_back({track.globalIndex(), trackQA.globalIndex()});
    }

    sortVectorOfArray(trackGIForTrackQAIndexList, 0); // sort the list //Its easy to search in a sorted list
    checkUniqueness(trackGIForTrackQAIndexList, 0);   // check the uniqueness of track.globalIndex()

    // create pointer table
    int currentIDXforCheck = 0;
    int listSize = trackGIForTrackQAIndexList.size();
    bool breakOnOverflow = false;

    for (const auto& track : tracks) {
      while (!breakOnOverflow && track.globalIndex() > trackGIForTrackQAIndexList[currentIDXforCheck][0]) {
        currentIDXforCheck++; // increment the currentIDXforCheck for missing or invalid cases e.g. value = -1;
        if (currentIDXforCheck >= listSize) {
          breakOnOverflow = true;
          break;
        }
      }
      if (!breakOnOverflow && track.globalIndex() == trackGIForTrackQAIndexList[currentIDXforCheck][0]) {
        genTrackToTracksQA(trackGIForTrackQAIndexList[currentIDXforCheck][1]);
      } else {
        genTrackToTracksQA(-1); // put a dummy index when track is not found in trackQA
      }
    }
  }
  PROCESS_SWITCH(CreatePointerTables, processTrackToTrackQAPointer, "processTrackToTrackQAPointer", false);

  std::vector<std::array<int64_t, 2>> trackGIForTMOIndexList;
  void processTrackToTrackMeanOccsPointer(aod::Tracks const& tracks, aod::TmoTrackIds const& tmoTrackIds)
  {
    trackGIForTMOIndexList.clear();
    int tmoCounter = -1;
    for (const auto& tmoTrackId : tmoTrackIds) {
      tmoCounter++;
      auto const& track = tmoTrackId.template track_as<aod::Tracks>();
      trackGIForTMOIndexList.push_back({track.globalIndex(), tmoCounter}); // tmoTrackId Global Index is not working :: tmoTrackId.globalIndex()});
    }
    sortVectorOfArray(trackGIForTMOIndexList, 0); // sort the list //Its easy to search in a sorted list
    checkUniqueness(trackGIForTMOIndexList, 0);   // check the uniqueness of track.globalIndex()

    // create pointer table
    int currentIDXforCheck = 0;
    int listSize = trackGIForTMOIndexList.size();
    bool breakOnOverflow = false;

    for (const auto& track : tracks) {
      while (!breakOnOverflow && track.globalIndex() > trackGIForTMOIndexList[currentIDXforCheck][0]) {
        currentIDXforCheck++; // increment the currentIDXforCheck for missing or invalid cases e.g. value = -1;
        if (currentIDXforCheck >= listSize) {
          breakOnOverflow = true;
          break;
        }
      }
      if (!breakOnOverflow && track.globalIndex() == trackGIForTMOIndexList[currentIDXforCheck][0]) {
        genTrackToTmo(trackGIForTMOIndexList[currentIDXforCheck][1]);
      } else {
        genTrackToTmo(-1); // put a dummy index when track is not found in trackQA
      }
    }
  }
  PROCESS_SWITCH(CreatePointerTables, processTrackToTrackMeanOccsPointer, "processTrackToTrackMeanOccsPointer", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<OccupancyTableProducer>(cfgc),
    adaptAnalysisTask<TrackMeanOccTableProducer>(cfgc),
    adaptAnalysisTask<CreatePointerTables>(cfgc)};
}
