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

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/CCDB/EventSelectionParams.h"
#include "PWGUD/DataModel/UDTables.h"

using namespace o2::framework;
using namespace o2::framework::expressions;

struct UpcCandProducer {
  Produces<o2::aod::UDCollisions> eventCandidates;
  Produces<o2::aod::UDCollisionsSels> eventCandidatesSels;
  Produces<o2::aod::UDTrackCollisionIDs> barrelCandIds;
  Produces<o2::aod::UDFwdTrackCollisionIDs> muonCandIds;

  bool fDoSemiFwd{false};

  Configurable<int> fCheckTPCPID{"checkTPCPID", 0, "Check TPC PID. Useful for central selection -- see `tpcPIDSwitch` option"};
  Configurable<int> fTPCPIDSwitch{"tpcPIDSwitch", 0, "PID switch: 0 -- two muons/pions, 1 -- two electrons, 2 -- electron + muon/pion"};

  Configurable<int> fNFwdProngs{"nFwdProngs", 2, "Matched forward tracks per candidate"};
  Configurable<int> fNBarProngs{"nBarProngs", 0, "Matched barrel tracks per candidate"};

  // todo: get parameters from CCDB
  EventSelectionParams par;

  // helper struct
  struct FITInfo {
    float ampFT0A = -1;
    float ampFT0C = -1;
    float timeFT0A = -999.;
    float timeFT0C = -999.;
    uint8_t triggerMaskFT0 = 0;
    float ampFDDA = -1;
    float ampFDDC = -1;
    float timeFDDA = -999.;
    float timeFDDC = -999.;
    uint8_t triggerMaskFDD = 0;
    float ampFV0A = -1;
    float timeFV0A = -999.;
    uint8_t triggerMaskFV0A = 0;
    // selection flags
    bool isBBFT0A = false;
    bool isBGFT0A = false;
    bool isBBFT0C = false;
    bool isBGFT0C = false;
    bool isBBFV0A = false;
    bool isBGFV0A = false;
    bool isBBFDDA = false;
    bool isBGFDDA = false;
    bool isBBFDDC = false;
    bool isBGFDDC = false;
  };

  void init(InitContext&)
  {
    // use "default" parameters
    par.fV0ABBlower = -3.0;  // ns
    par.fV0ABBupper = +2.0;  // ns
    par.fV0ABGlower = 2.0;   // ns
    par.fV0ABGupper = 5.0;   // ns
    par.fFDABBlower = -3.0;  // ns
    par.fFDABBupper = +3.0;  // ns
    par.fFDABGlower = 10.0;  // ns
    par.fFDABGupper = 13.0;  // ns
    par.fFDCBBlower = -3.0;  // ns
    par.fFDCBBupper = +3.0;  // ns
    par.fFDCBGlower = -10.0; // ns
    par.fFDCBGupper = -3.0;  // ns
    par.fT0ABBlower = -1.0;  // ns
    par.fT0ABBupper = +1.0;  // ns
    par.fT0CBBlower = -1.0;  // ns
    par.fT0CBBupper = +1.0;  // ns
  }

  // filter tracks in the central barrel using TPC PID
  // return `true` if a candidate passes "number of tracks" requirement from `nBarProngs`, `false` otherwise
  // if candidate passes, `filteredTrackIDs` contain track IDs of passed tracks
  template <typename TBarrelTracks>
  bool checkTPCPID(TBarrelTracks* tracks, std::vector<int32_t> const& trackIDs, std::vector<int32_t>& filteredTrackIDs)
  {
    int32_t nIDs = trackIDs.size();
    std::vector<bool> pidFlagsEl(nIDs, false);
    std::vector<bool> pidFlagsPi(nIDs, false);
    int32_t nEl = 0;
    int32_t nPi = 0;

    for (int32_t itr = 0; itr < nIDs; itr++) {
      const auto& tr = tracks->iteratorAt(trackIDs[itr]);
      float tpcNSigmaEl = tr.tpcNSigmaEl();
      float tpcNSigmaPi = tr.tpcNSigmaPi();

      bool isEl = std::abs(tpcNSigmaEl) < 3.0f && std::abs(tpcNSigmaEl) < std::abs(tpcNSigmaPi);
      bool isPi = std::abs(tpcNSigmaPi) < 3.0f && std::abs(tpcNSigmaEl) > std::abs(tpcNSigmaPi);

      pidFlagsEl[itr] = isEl;
      if (isEl) {
        nEl++;
      }

      pidFlagsPi[itr] = isPi;
      if (isPi) {
        nPi++;
      }
    }

    bool pass = false;

    // two muons/pions
    if (fTPCPIDSwitch == 0 && nPi == 2 && nEl == 0) {
      for (int32_t itr = 0; itr < nIDs; itr++) {
        if (pidFlagsPi[itr]) {
          filteredTrackIDs.push_back(trackIDs[itr]);
        }
      }
      pass = true;
    }

    // two electrons
    if (fTPCPIDSwitch == 1 && nEl == 2) {
      for (int32_t itr = 0; itr < nIDs; itr++) {
        if (pidFlagsEl[itr]) {
          filteredTrackIDs.push_back(trackIDs[itr]);
        }
      }
      pass = true;
    }

    // electron + muon/pion
    if (fTPCPIDSwitch == 2 && nEl == 1 && nPi == 1) {
      for (int32_t itr = 0; itr < nIDs; itr++) {
        if (pidFlagsEl[itr]) {
          filteredTrackIDs.push_back(trackIDs[itr]);
        }
        if (pidFlagsPi[itr]) {
          filteredTrackIDs.push_back(trackIDs[itr]);
        }
      }
      pass = true;
    }

    return pass;
  }

  void processFITInfo(uint64_t bc,
                      std::map<uint64_t, int32_t>& bcsWithFT0,
                      std::map<uint64_t, int32_t>& bcsWithFDD,
                      std::map<uint64_t, int32_t>& bcsWithFV0A,
                      o2::aod::FT0s const& ft0s,
                      o2::aod::FDDs const& fdds,
                      o2::aod::FV0As const& fv0as,
                      FITInfo& fitInfo)
  {
    float timeV0ABG = -999.f;
    float timeT0ABG = -999.f;
    float timeT0CBG = -999.f;
    float timeFDABG = -999.f;
    float timeFDCBG = -999.f;

    // check FIT info in the same BC
    auto it = bcsWithFT0.find(bc);
    if (it != bcsWithFT0.end()) {
      const auto& ft0 = ft0s.iteratorAt(it->second);
      fitInfo.timeFT0A = ft0.timeA();
      fitInfo.timeFT0C = ft0.timeC();
      const auto& ampsA = ft0.amplitudeA();
      const auto& ampsC = ft0.amplitudeC();
      fitInfo.ampFT0A = 0.;
      for (auto amp : ampsA) {
        fitInfo.ampFT0A += amp;
      }
      fitInfo.ampFT0C = 0.;
      for (auto amp : ampsC) {
        fitInfo.ampFT0C += amp;
      }
      fitInfo.triggerMaskFT0 = ft0.triggerMask();
    }

    it = bcsWithFDD.find(bc);
    if (it != bcsWithFDD.end()) {
      const auto& fdd = fdds.iteratorAt(it->second);
      fitInfo.timeFDDA = fdd.timeA();
      fitInfo.timeFDDC = fdd.timeC();
      const auto& ampsA = fdd.chargeA();
      const auto& ampsC = fdd.chargeC();
      fitInfo.ampFDDA = 0.;
      for (auto amp : ampsA) {
        fitInfo.ampFDDA += amp;
      }
      fitInfo.ampFDDC = 0.;
      for (auto amp : ampsC) {
        fitInfo.ampFDDC += amp;
      }
      fitInfo.triggerMaskFDD = fdd.triggerMask();
    }

    it = bcsWithFV0A.find(bc);
    if (it != bcsWithFV0A.end()) {
      const auto& fv0a = fv0as.iteratorAt(it->second);
      fitInfo.timeFV0A = fv0a.time();
      const auto& amps = fv0a.amplitude();
      fitInfo.ampFV0A = 0.;
      for (auto amp : amps) {
        fitInfo.ampFV0A += amp;
      }
      fitInfo.triggerMaskFV0A = fv0a.triggerMask();
    }

    // check beam-gas
    it = bcsWithFT0.find(bc - 1);
    if (it != bcsWithFT0.end()) {
      const auto& ft0 = ft0s.iteratorAt(it->second);
      timeT0ABG = ft0.timeA();
      timeT0CBG = ft0.timeC();
    }

    it = bcsWithFDD.find(bc - 5);
    if (it != bcsWithFDD.end()) {
      const auto& ft0 = fdds.iteratorAt(it->second);
      timeFDABG = ft0.timeA();
      timeFDCBG = ft0.timeC();
    }

    it = bcsWithFV0A.find(bc - 1);
    if (it != bcsWithFV0A.end()) {
      const auto& fv0a = fv0as.iteratorAt(it->second);
      timeV0ABG = fv0a.time();
    }

    // beam-gas flags
    fitInfo.isBGFV0A = timeV0ABG > par.fV0ABGlower && timeV0ABG < par.fV0ABGupper;
    fitInfo.isBGFDDA = timeFDABG > par.fFDABGlower && timeFDABG < par.fFDABGupper;
    fitInfo.isBGFDDC = timeFDCBG > par.fFDCBGlower && timeFDCBG < par.fFDCBGupper;
    fitInfo.isBGFT0A = timeT0ABG > par.fT0ABGlower && timeT0ABG < par.fT0ABGupper;
    fitInfo.isBGFT0C = timeT0CBG > par.fT0CBGlower && timeT0CBG < par.fT0CBGupper;

    // beam-beam flags
    fitInfo.isBBFT0A = fitInfo.timeFT0A > par.fT0ABBlower && fitInfo.timeFT0A < par.fT0ABBupper;
    fitInfo.isBBFT0C = fitInfo.timeFT0C > par.fT0CBBlower && fitInfo.timeFT0C < par.fT0CBBupper;
    fitInfo.isBBFV0A = fitInfo.timeFV0A > par.fV0ABBlower && fitInfo.timeFV0A < par.fV0ABBupper;
    fitInfo.isBBFDDA = fitInfo.timeFDDA > par.fFDABBlower && fitInfo.timeFDDA < par.fFDABBupper;
    fitInfo.isBBFDDC = fitInfo.timeFDDC > par.fFDCBBlower && fitInfo.timeFDDC < par.fFDCBBupper;
  }

  template <typename TFwdTracks, typename TBarrelTracks, typename TBCs>
  void createCandidates(TFwdTracks* fwdTracks,
                        TBarrelTracks* barTracks,
                        TBCs const& bcs,
                        o2::aod::FT0s const& ft0s,
                        o2::aod::FDDs const& fdds,
                        o2::aod::FV0As const& fv0as)
  {
    // map track IDs to the respective event candidate IDs
    std::vector<int32_t> barTrackCandIds;
    std::vector<int32_t> fwdTrackCandIds;

    std::map<uint64_t, int32_t> BCsWithFT0;
    // collect BCs with FT0 signals
    for (const auto& ft0 : ft0s) {
      uint64_t bc = ft0.bc().globalBC();
      BCsWithFT0[bc] = ft0.globalIndex();
    }

    std::map<uint64_t, int32_t> BCsWithFDD;
    // collect BCs with FDD signals
    for (const auto& fdd : fdds) {
      uint64_t bc = fdd.bc().globalBC();
      BCsWithFDD[bc] = fdd.globalIndex();
    }

    std::map<uint64_t, int32_t> BCsWithFV0A;
    // collect BCs with FV0A signals
    for (const auto& fv0a : fv0as) {
      uint64_t bc = fv0a.bc().globalBC();
      BCsWithFV0A[bc] = fv0a.globalIndex();
    }

    // pairs of global BCs and vectors of matched track IDs:
    // global BC <-> <vector of fwd. trackIDs, vector of barrel trackIDs>
    std::map<uint64_t, std::pair<std::vector<int32_t>, std::vector<int32_t>>> bcsMatchedTrIds;

    // forward matching
    if (fwdTracks != nullptr) {
      fwdTrackCandIds.resize(fwdTracks->size(), -1);
      for (const auto& fwdTr : *fwdTracks) {
        uint64_t bc = fwdTr.globalBC();
        // search for BC:
        //  if found -> store track ID to vector of matched tracks
        //  else make a new vector of matched tracks and store track ID
        auto it = bcsMatchedTrIds.find(bc);
        if (it != bcsMatchedTrIds.end()) {
          it->second.first.emplace_back(fwdTr.globalIndex());
        } else {
          bcsMatchedTrIds[bc] = std::make_pair(std::vector<int32_t>(1, fwdTr.globalIndex()), std::vector<int32_t>());
        }
      }
    }

    // central barrel tracks
    if (barTracks != nullptr) {
      barTrackCandIds.resize(barTracks->size(), -1);
      for (const auto& barTr : *barTracks) {
        uint64_t bc = barTr.globalBC();
        // search for BC:
        //  if found -> store track ID to vector of matched tracks
        //  else make a new vector of matched tracks and store track ID
        auto it = bcsMatchedTrIds.find(bc);
        if (it != bcsMatchedTrIds.end()) {
          it->second.second.emplace_back(barTr.globalIndex());
        } else if (!fDoSemiFwd) { // tag central-barrel tracks to forward tracks in semiforward case
          bcsMatchedTrIds[bc] = std::make_pair(std::vector<int32_t>(), std::vector<int32_t>(1, barTr.globalIndex()));
        }
      }
    }

    // todo: calculate position of UD collision?
    float dummyX = 0.;
    float dummyY = 0.;
    float dummyZ = 0.;

    // storing n-prong matches
    int32_t candID = 0;
    for (const auto& item : bcsMatchedTrIds) {
      uint64_t bc = item.first;
      std::vector<int32_t> fwdTrackIDs = item.second.first;
      std::vector<int32_t> barTrackIDs = item.second.second;
      int32_t nFwdTracks = fwdTrackIDs.size();
      int32_t nBarTracks = barTrackIDs.size();
      // check number of tracks in a candidate
      bool checkForward = nFwdTracks == fNFwdProngs;
      bool checkCentral = nBarTracks == fNBarProngs;
      // check TPC PID if needed
      std::vector<int32_t> filteredTrackIDs;
      if (fCheckTPCPID) {
        checkCentral = checkTPCPID(barTracks, barTrackIDs, filteredTrackIDs);
        if (checkCentral) {
          barTrackIDs.swap(filteredTrackIDs);
          nBarTracks = barTrackIDs.size();
        }
      }
      if (!checkForward || !checkCentral) {
        continue;
      }
      float RgtrwTOF = 0.;
      for (auto id : barTrackIDs) {
        const auto& tr = barTracks->iteratorAt(id);
        if (tr.hasTOF()) {
          RgtrwTOF++;
        }
      }
      RgtrwTOF = nBarTracks != 0 ? RgtrwTOF / (float)nBarTracks : 0.;
      if (RgtrwTOF == 0 && fNBarProngs != 0) { // require at least 1 TOF track in central and semiforward cases
        continue;
      }
      int8_t netCharge = 0;
      uint16_t numContrib = nFwdTracks + nBarTracks;
      for (auto id : fwdTrackIDs) {
        fwdTrackCandIds[id] = candID;
        const auto& tr = fwdTracks->iteratorAt(id);
        netCharge += tr.sign();
      }
      for (auto id : barTrackIDs) {
        barTrackCandIds[id] = candID;
        const auto& tr = barTracks->iteratorAt(id);
        netCharge += tr.sign();
      }
      // fetching FT0, FDD, FV0 information
      // if there is no relevant signal, dummy info will be used
      FITInfo fitInfo;
      processFITInfo(bc, BCsWithFT0, BCsWithFDD, BCsWithFV0A, ft0s, fdds, fv0as, fitInfo);
      int32_t runNumber = bcs.iteratorAt(0).runNumber();
      eventCandidates(bc, runNumber, dummyX, dummyY, dummyZ, numContrib, netCharge, RgtrwTOF);
      eventCandidatesSels(fitInfo.ampFT0A, fitInfo.ampFT0C, fitInfo.timeFT0A, fitInfo.timeFT0C, fitInfo.triggerMaskFT0,
                          fitInfo.ampFDDA, fitInfo.ampFDDC, fitInfo.timeFDDA, fitInfo.timeFDDC, fitInfo.triggerMaskFDD,
                          fitInfo.ampFV0A, fitInfo.timeFV0A, fitInfo.triggerMaskFV0A,
                          fitInfo.isBBFT0A, fitInfo.isBBFT0C, fitInfo.isBGFT0A, fitInfo.isBGFT0C,
                          fitInfo.isBBFV0A, fitInfo.isBGFV0A,
                          fitInfo.isBBFDDA, fitInfo.isBBFDDC, fitInfo.isBGFDDA, fitInfo.isBGFDDC);
      candID++;
    }

    bcsMatchedTrIds.clear();
    BCsWithFT0.clear();
    BCsWithFDD.clear();
    BCsWithFV0A.clear();

    if (fwdTracks != nullptr) {
      for (const auto& fwdTr : *fwdTracks) {
        muonCandIds(fwdTrackCandIds[fwdTr.globalIndex()]);
      }
    }

    fwdTrackCandIds.clear();

    if (barTracks != nullptr) {
      for (const auto& barTr : *barTracks) {
        barrelCandIds(barTrackCandIds[barTr.globalIndex()]);
      }
    }

    barTrackCandIds.clear();
  }

  using BarrelTracks = o2::soa::Join<o2::aod::UDTracks, o2::aod::UDTracksExtra, o2::aod::UDTracksPID>;

  // create candidates for forward region
  void processFwd(o2::aod::UDFwdTracks const& muonTracks,
                  o2::aod::BCs const& bcs,
                  o2::aod::FT0s const& ft0s,
                  o2::aod::FDDs const& fdds,
                  o2::aod::FV0As const& fv0as)
  {
    fDoSemiFwd = false;
    createCandidates(&muonTracks, (BarrelTracks*)nullptr, bcs, ft0s, fdds, fv0as);
  }

  // create candidates for semiforward region
  void processSemiFwd(o2::aod::UDFwdTracks const& muonTracks,
                      BarrelTracks const& barTracks,
                      o2::aod::BCs const& bcs,
                      o2::aod::FT0s const& ft0s,
                      o2::aod::FDDs const& fdds,
                      o2::aod::FV0As const& fv0as)
  {
    fDoSemiFwd = true;
    createCandidates(&muonTracks, &barTracks, bcs, ft0s, fdds, fv0as);
  }

  // create candidates for central region
  void processCentral(o2::aod::UDFwdTracks const& muonTracks,
                      BarrelTracks const& barTracks,
                      o2::aod::BCs const& bcs,
                      o2::aod::FT0s const& ft0s,
                      o2::aod::FDDs const& fdds,
                      o2::aod::FV0As const& fv0as)
  {
    fDoSemiFwd = false;
    createCandidates((o2::aod::UDFwdTracks*)nullptr, &barTracks, bcs, ft0s, fdds, fv0as);
  }

  PROCESS_SWITCH(UpcCandProducer, processFwd, "Produce candidates for forward rapidities", false);
  PROCESS_SWITCH(UpcCandProducer, processSemiFwd, "Produce candidates in semiforward region", false);
  PROCESS_SWITCH(UpcCandProducer, processCentral, "Produce candidates in central region", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<UpcCandProducer>(cfgc)};
}
