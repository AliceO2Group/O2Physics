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

/// \file caloClusterProducer.cxx
/// \brief Produces PHOS clusters from PHOS cells
///
/// \author Dmitri Peresunko <Dmitri.Peresunko@cern.ch>

#include <algorithm>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "Framework/ConfigParamSpec.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/CaloClusters.h"
#include "Common/Core/trackUtilities.h"
#include "ReconstructionDataFormats/TrackParametrization.h"
#include "DetectorsBase/Propagator.h"

#include "CommonUtils/NameConf.h"
#include "CCDB/BasicCCDBManager.h"
#include "SimulationDataFormat/MCTruthContainer.h"

#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsPHOS/Cell.h"
#include "DataFormatsPHOS/Cluster.h"
#include "DataFormatsPHOS/TriggerRecord.h"
#include "DataFormatsPHOS/MCLabel.h"
#include "DataFormatsPHOS/BadChannelsMap.h"
#include "DataFormatsPHOS/CalibParams.h"
#include "PHOSBase/Geometry.h"
#include "PHOSReconstruction/Clusterer.h"

using namespace o2::framework;
using namespace o2;

using McCells = o2::soa::Join<aod::Calos, aod::McCaloLabels_001>;

struct CaloClusterProducer {
  Produces<aod::CaloClusters> clucursor;
  Produces<aod::CaloAmbiguousClusters> cluambcursor;
  Produces<aod::PHOSMatchedTracks> matchedTracks;
  Produces<aod::PHOSCluLabels> clumccursor;
  Produces<aod::PHOSAmbCluLabels> cluambmccursor;

  Configurable<bool> isMC{"isMC", 0, "0 - data, 1 - MC"};
  Configurable<bool> useCoreE{"useCoreE", 0, "0 - full energy, 1 - core energy"};
  Configurable<bool> skipL1phase{"skipL1phase", false, "skip or apply L1phase time correction"};
  Configurable<int> nonlinType{"nonlinType", 1, "0:no corr, 1: data, 2: MC"};
  Configurable<std::vector<double>> cpvMinE{"cpvMinE", {20., 50., 50.}, "minimal CPV cluster amplitude per module"};
  Configurable<std::string> badMapPath{"badMapPath", "PHS/Calib/BadMap", "path to BadMap snapshot"};
  Configurable<std::string> calibPath{"calibPath", "PHS/Calib/CalibParams", "path to Calibration snapshot"};
  Configurable<std::string> l1phasePath{"l1phasePath", "PHS/Calib/L1phase", "path to L1phase snapshot"};

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  std::unique_ptr<o2::phos::Geometry> geomPHOS;
  std::unique_ptr<o2::phos::Clusterer> clusterizerPHOS;
  std::vector<o2::phos::Cell> phosCells;
  std::vector<o2::phos::TriggerRecord> phosCellTRs;
  std::vector<o2::phos::CluElement> outputCluElements;
  std::vector<o2::phos::Cluster> outputPHOSClusters;
  std::vector<o2::phos::TriggerRecord> outputPHOSClusterTrigRecs;
  std::vector<int> mclabels;
  std::vector<float> mcamplitudes;

  int mRunNumber{0};

  static constexpr int16_t kCpvX = 7; // grid 6 steps along z and 7 along phi as largest match ellips 20x20 cm
  static constexpr int16_t kCpvZ = 6;
  static constexpr int16_t kCpvCells = 4 * kCpvX * kCpvZ; // 4 modules
  static constexpr float kCpvMaxX = 73;                   // max CPV coordinate phi
  static constexpr float kCpvMaxZ = 63;                   // max CPV coordinate z

  class TrackMatch
  {
   public:
    TrackMatch() = default;
    TrackMatch(float x, float z, int i) : pX(x), pZ(z), indx(i) {}
    ~TrackMatch() = default;

   public:
    float pX = 9999.; // X (phi) track coordinate in PHOS plane
    float pZ = 9999.; // Z (theta) track coordinate in PHOS plane
    int indx = -1;    // track global index
  };

  class TrackTrigRec
  {
   public:
    TrackTrigRec() = default;
    ~TrackTrigRec() = default;

   public:
    int64_t mTR;           // BC ref
    int mStart[kCpvCells]; // X (phi) track coordinate in PHOS plane
    int mEnd[kCpvCells];   // Z (theta) track coordinate in PHOS plane
  };

  void init(o2::framework::InitContext&)
  {
    ccdb->setURL(o2::base::NameConf::getCCDBServer());
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    geomPHOS = std::make_unique<o2::phos::Geometry>("PHOS");
    clusterizerPHOS = std::make_unique<o2::phos::Clusterer>();
  }

  void processStandalone(o2::aod::BCsWithTimestamps const& bcs,
                         o2::aod::Collisions const& colls,
                         o2::aod::Calos const& cells,
                         o2::aod::CaloTriggers const&,
                         o2::aod::CPVClusters const& cpvs)
  {

    int64_t timestamp = 0;
    if (bcs.begin() != bcs.end()) {
      timestamp = bcs.begin().timestamp(); // timestamp for CCDB object retrieval
    } else {
      return;
    }
    std::map<int64_t, int> bcMap;
    int bcId = 0;
    for (const auto& bc : bcs) {
      bcMap[bc.globalBC()] = bcId;
      bcId++;
    }

    // If several collisions appear in BC, choose one with largers number of contributors
    std::map<int64_t, int> colMap;
    int colId = 0;
    for (const auto& cl : colls) {
      auto colbc = colMap.find(cl.bc_as<aod::BCsWithTimestamps>().globalBC());
      if (colbc == colMap.end()) { // single collision per BC
        colMap[cl.bc_as<aod::BCsWithTimestamps>().globalBC()] = colId;
      } else { // not unique collision per BC
        auto coll2 = colls.begin() + colbc->second;
        if (cl.numContrib() > coll2.numContrib()) {
          colMap[cl.bc_as<aod::BCsWithTimestamps>().globalBC()] = colId;
        }
      }
      colId++;
    }

    // Fill list of cells and cell TrigRecs per TF as an input for clusterizer
    // clusterize
    // Fill output table

    // calibration may be updated by CCDB fetcher
    const o2::phos::BadChannelsMap* badMap = ccdb->getForTimeStamp<o2::phos::BadChannelsMap>(badMapPath, timestamp);
    const o2::phos::CalibParams* calibParams = ccdb->getForTimeStamp<o2::phos::CalibParams>(calibPath, timestamp);

    if (!isMC && !skipL1phase) {
      const std::vector<int>* vec = ccdb->getForTimeStamp<std::vector<int>>(l1phasePath, timestamp);
      if (vec) {
        clusterizerPHOS->setL1phase((*vec)[0]);
      } else {
        LOG(fatal) << "Can not get PHOS L1phase calibration";
      }
    }

    if (badMap) {
      clusterizerPHOS->setBadMap(badMap);
    } else {
      LOG(fatal) << "Can not get PHOS Bad Map";
    }
    if (calibParams) {
      clusterizerPHOS->setCalibration(calibParams);
    } else {
      LOG(fatal) << "Can not get PHOS calibration";
    }

    phosCells.clear();
    phosCells.reserve(cells.size());
    phosCellTRs.clear();
    phosCellTRs.reserve(bcs.size());
    outputCluElements.clear();
    outputPHOSClusters.clear();
    outputPHOSClusterTrigRecs.clear();

    o2::InteractionRecord ir;
    const int kPHOS = 0;
    for (const auto& c : cells) {
      if (c.caloType() != kPHOS) // PHOS
        continue;
      // Fix for bug in trigger digits
      if ((c.cellType() == phos::TRU2x2 || c.cellType() == phos::TRU4x4) && c.cellNumber() == 0) {
        continue;
      }
      if (phosCellTRs.size() == 0) { // first cell, first TrigRec
        ir.setFromLong(c.bc_as<aod::BCsWithTimestamps>().globalBC());
        phosCellTRs.emplace_back(ir, 0, 0); // BC,first cell, ncells
      }
      if (static_cast<uint64_t>(phosCellTRs.back().getBCData().toLong()) != c.bc_as<aod::BCsWithTimestamps>().globalBC()) { // switch to new BC
        // switch to another BC: set size and create next TriRec
        phosCellTRs.back().setNumberOfObjects(phosCells.size() - phosCellTRs.back().getFirstEntry());
        // Next event/trig rec.
        ir.setFromLong(c.bc_as<aod::BCsWithTimestamps>().globalBC());
        phosCellTRs.emplace_back(ir, phosCells.size(), 0);
      }
      phosCells.emplace_back(c.cellNumber(), c.amplitude(), c.time(),
                             static_cast<o2::phos::ChannelType_t>(c.cellType()));
    }
    // Set number of cells in last TrigRec
    if (phosCellTRs.size() > 0) {
      phosCellTRs.back().setNumberOfObjects(phosCells.size() - phosCellTRs.back().getFirstEntry());
    }

    // clusterize
    o2::dataformats::MCTruthContainer<o2::phos::MCLabel> dummyMC;
    clusterizerPHOS->processCells(phosCells, phosCellTRs, nullptr,
                                  outputPHOSClusters, outputCluElements, outputPHOSClusterTrigRecs, dummyMC);

    // Find  CPV clusters corresponding to PHOS trigger records
    std::vector<std::pair<float, float>> cpvMatchPoints[kCpvCells];
    // Number of entries in each cell per TrigRecord
    std::vector<TrackTrigRec> cpvNMatchPoints;
    cpvNMatchPoints.reserve(outputPHOSClusterTrigRecs.size());
    int64_t curBC = -1;
    if (cpvs.begin() != cpvs.end()) {
      curBC = cpvs.begin().bc_as<aod::BCsWithTimestamps>().globalBC();
      cpvNMatchPoints.emplace_back();
      cpvNMatchPoints.back().mTR = curBC;
      for (int i = kCpvCells; i--;) {
        cpvNMatchPoints.back().mStart[i] = 0;
      }
    }

    for (const auto& cpvclu : cpvs) {
      if (static_cast<int64_t>(cpvclu.bc_as<aod::BCsWithTimestamps>().globalBC()) != curBC) { // new BC
        // mark last entry in previous range
        for (int i = kCpvCells; i--;) {
          cpvNMatchPoints.back().mEnd[i] = cpvMatchPoints[i].size();
        }
        curBC = cpvclu.bc_as<aod::BCsWithTimestamps>().globalBC();
        cpvNMatchPoints.back() = cpvNMatchPoints.emplace_back();
        cpvNMatchPoints.back().mTR = curBC;
        for (int i = kCpvCells; i--;) {
          cpvNMatchPoints.back().mStart[i] = cpvMatchPoints[i].size();
        }
      }
      // anyway add coordinates
      if (cpvclu.amplitude() < static_cast<std::vector<double>>(cpvMinE)[static_cast<int>(cpvclu.moduleNumber()) - 2]) {
        continue;
      }
      int index = cpvMatchIndex(cpvclu.moduleNumber(), cpvclu.posX(), cpvclu.posZ());
      cpvMatchPoints[index].emplace_back(cpvclu.posX(), cpvclu.posZ());
    }
    if (cpvNMatchPoints.size() > 0) {
      for (int i = kCpvCells; i--;) {
        cpvNMatchPoints.back().mEnd[i] = cpvMatchPoints[i].size();
      }
    }

    // Fill output
    for (const auto& cluTR : outputPHOSClusterTrigRecs) {
      int firstClusterInEvent = cluTR.getFirstEntry();
      int lastClusterInEvent = firstClusterInEvent + cluTR.getNumberOfObjects();

      // Extract primary vertex
      TVector3 vtx = {0., 0., 0.}; // default, if not collision will be found
      int colId = -1;
      auto coliter = colMap.find(cluTR.getBCData().toLong());
      if (coliter != colMap.end()) { // get vertex from collision
        // find collision corresponding to current BC
        auto clvtx = colls.begin() + coliter->second;
        vtx.SetXYZ(clvtx.posX(), clvtx.posY(), clvtx.posZ());
        colId = coliter->second;
      }

      bool cpvExist = false;
      // find cpvTR for this BC
      auto cpvPoints = cpvNMatchPoints.begin();
      while (cpvPoints != cpvNMatchPoints.end()) {
        if (cpvPoints->mTR == cluTR.getBCData().toLong()) {
          cpvExist = true;
          break;
        }
        cpvPoints++;
      }

      for (int i = firstClusterInEvent; i < lastClusterInEvent; i++) {
        o2::phos::Cluster& clu = outputPHOSClusters[i];
        float e = (useCoreE) ? clu.getCoreEnergy() : clu.getEnergy();
        if (e == 0) {
          continue;
        }
        float posX, posZ;
        clu.getLocalPosition(posX, posZ);

        // Correction for the depth of the shower starting point (TDR p 127)
        const float para = 0.925;
        const float parb = 6.52;
        float depth = para * std::log(e) + parb;
        posX -= posX * depth / 460.;
        posZ -= (posZ - vtx.Z()) * depth / 460.;

        int mod = clu.module();
        TVector3 globaPos;
        geomPHOS->local2Global(mod, posX, posZ, globaPos);

        TVector3 mom = globaPos - vtx;
        if (mom.Mag() == 0) { // should not happpen
          continue;
        }

        e = nonlinearity(e);

        mom.SetMag(e);

        float cpvdist = 99.;
        const float cellSizeX = 2 * kCpvMaxX / kCpvX;
        const float cellSizeZ = 2 * kCpvMaxZ / kCpvZ;
        // look 9 CPV regions around PHOS cluster

        if (mod >= 2 && cpvExist) { // CPV exist in mods 2,3,4
          int phosIndex = cpvMatchIndex(mod, posX, posZ);
          std::vector<int> regions;
          regions.push_back(phosIndex);
          if (posX > -kCpvMaxX + cellSizeX) {
            if (posZ > -kCpvMaxZ + cellSizeZ) { // bottom left
              regions.push_back(phosIndex - kCpvZ - 1);
            }
            regions.push_back(phosIndex - kCpvZ);
            if (posZ < kCpvMaxZ - cellSizeZ) { // top left
              regions.push_back(phosIndex - kCpvZ + 1);
            }
          }
          if (posZ > -kCpvMaxZ + cellSizeZ) { // bottom
            regions.push_back(phosIndex - 1);
          }
          if (posZ < kCpvMaxZ - cellSizeZ) { // top
            regions.push_back(phosIndex + 1);
          }
          if (posX < kCpvMaxX - cellSizeX) {
            if (posZ > -kCpvMaxZ + cellSizeZ) { // bottom right
              regions.push_back(phosIndex + kCpvZ - 1);
            }
            regions.push_back(phosIndex + kCpvZ);
            if (posZ < kCpvMaxZ - cellSizeZ) { // top right
              regions.push_back(phosIndex + kCpvZ + 1);
            }
          }
          float sigmaX = 1. / std::min(5.2, 1.111 + 0.56 * std::exp(-0.031 * e * e) + 4.8 / std::pow(e + 0.61, 3)); // inverse sigma X
          float sigmaZ = 1. / std::min(3.3, 1.12 + 0.35 * std::exp(-0.032 * e * e) + 0.75 / std::pow(e + 0.24, 3)); // inverse sigma Z

          for (const int& indx : regions) {
            if (indx >= 0 && indx < kCpvCells) {
              for (int ii = cpvPoints->mStart[indx]; ii < cpvPoints->mEnd[indx]; ii++) {
                auto p = cpvMatchPoints[indx][ii];
                float d = std::pow((p.first - posX) * sigmaX, 2) + std::pow((p.second - posZ) * sigmaZ, 2);
                if (d < cpvdist) {
                  cpvdist = d;
                }
              }
            }
          }
        }
        if (cpvdist != 99.) {           // was evaluated
          cpvdist = std::sqrt(cpvdist); // was squared
        }
        int cpvindex = -2; // -2 no CPV in event
        if (cpvExist) {
          cpvindex = -1; // there were CPV clusters
        }

        float lambdaShort = 0., lambdaLong = 0.;
        clu.getElipsAxis(lambdaShort, lambdaLong);

        // Clear Collision assignment
        if (colId == -1) {
          // Ambiguos Collision assignment
          cluambcursor(
            bcMap[cluTR.getBCData().toLong()],
            mom.X(), mom.Y(), mom.Z(), e,
            mod, clu.getMultiplicity(), posX, posZ,
            globaPos.X(), globaPos.Y(), globaPos.Z(),
            clu.getTime(), clu.getNExMax(),
            lambdaLong, lambdaShort,
            cpvdist, cpvindex,
            clu.firedTrigger(),
            clu.getDistanceToBadChannel());

        } else { // Normal collision
          auto col = colls.begin() + colId;
          clucursor(
            col,
            mom.X(), mom.Y(), mom.Z(), e,
            mod, clu.getMultiplicity(), posX, posZ,
            globaPos.X(), globaPos.Y(), globaPos.Z(),
            clu.getTime(), clu.getNExMax(),
            lambdaLong, lambdaShort,
            cpvdist, cpvindex,
            clu.firedTrigger(),
            clu.getDistanceToBadChannel());
        }
      }
    }
  }

  PROCESS_SWITCH(CaloClusterProducer, processStandalone, "Process PHOS and CPV only", true);

  void processStandaloneMC(o2::aod::BCsWithTimestamps const& bcs,
                           o2::aod::Collisions const& colls,
                           McCells const& cells,
                           o2::aod::CaloTriggers const&,
                           o2::aod::CPVClusters const& cpvs)
  {

    int64_t timestamp = 0;
    if (bcs.begin() != bcs.end()) {
      timestamp = bcs.begin().timestamp(); // timestamp for CCDB object retrieval
    } else {
      return;
    }
    std::map<int64_t, int> bcMap;
    int bcId = 0;
    for (auto const& bc : bcs) {
      bcMap[bc.globalBC()] = bcId;
      bcId++;
    }

    // If several collisions appear in BC, choose one with largers number of contributors
    std::map<int64_t, int> colMap;
    int colId = 0;
    for (auto const& cl : colls) {
      auto colbc = colMap.find(cl.bc_as<aod::BCsWithTimestamps>().globalBC());
      if (colbc == colMap.end()) { // single collision per BC
        colMap[cl.bc_as<aod::BCsWithTimestamps>().globalBC()] = colId;
      } else { // not unique collision per BC
        auto coll2 = colls.begin() + colbc->second;
        if (cl.numContrib() > coll2.numContrib()) {
          colMap[cl.bc_as<aod::BCsWithTimestamps>().globalBC()] = colId;
        }
      }
      colId++;
    }

    // Fill list of cells and cell TrigRecs per TF as an input for clusterizer
    // clusterize
    // Fill output table

    // calibration may be updated by CCDB fetcher
    const o2::phos::BadChannelsMap* badMap = ccdb->getForTimeStamp<o2::phos::BadChannelsMap>(badMapPath, timestamp);
    const o2::phos::CalibParams* calibParams = ccdb->getForTimeStamp<o2::phos::CalibParams>(calibPath, timestamp);

    if (!isMC && !skipL1phase) {
      const std::vector<int>* vec = ccdb->getForTimeStamp<std::vector<int>>(l1phasePath, timestamp);
      if (vec) {
        clusterizerPHOS->setL1phase((*vec)[0]);
      } else {
        LOG(fatal) << "Can not get PHOS L1phase calibration";
      }
    }

    if (badMap) {
      clusterizerPHOS->setBadMap(badMap);
    } else {
      LOG(fatal) << "Can not get PHOS Bad Map";
    }
    if (calibParams) {
      clusterizerPHOS->setCalibration(calibParams);
    } else {
      LOG(fatal) << "Can not get PHOS calibration";
    }

    phosCells.clear();
    phosCells.reserve(cells.size());
    phosCellTRs.clear();
    phosCellTRs.reserve(bcs.size());
    outputCluElements.clear();
    outputPHOSClusters.clear();
    outputPHOSClusterTrigRecs.clear();

    o2::InteractionRecord ir;
    const int kPHOS = 0;
    o2::dataformats::MCTruthContainer<o2::phos::MCLabel> cellTruth;
    for (const auto& c : cells) {
      if (c.caloType() != kPHOS) // PHOS
        continue;
      // Fix for bug in trigger digits
      if ((c.cellType() == phos::TRU2x2 || c.cellType() == phos::TRU4x4) && c.cellNumber() == 0) {
        continue;
      }
      if (phosCellTRs.size() == 0) { // first cell, first TrigRec
        ir.setFromLong(c.bc_as<aod::BCsWithTimestamps>().globalBC());
        phosCellTRs.emplace_back(ir, 0, 0); // BC,first cell, ncells
      }
      if (static_cast<uint64_t>(phosCellTRs.back().getBCData().toLong()) != c.bc_as<aod::BCsWithTimestamps>().globalBC()) { // switch to new BC
        // switch to another BC: set size and create next TriRec
        phosCellTRs.back().setNumberOfObjects(phosCells.size() - phosCellTRs.back().getFirstEntry());
        // Next event/trig rec.
        ir.setFromLong(c.bc_as<aod::BCsWithTimestamps>().globalBC());
        phosCellTRs.emplace_back(ir, phosCells.size(), 0);
      }
      phosCells.emplace_back(c.cellNumber(), c.amplitude(), c.time(),
                             static_cast<o2::phos::ChannelType_t>(c.cellType()));
      // process MC info
      auto edep = c.amplitudeA();
      // find particle with non-zero E deposited fraction; if several, use one with largest label (may be daughter - to keep full history)
      float ed = 0;
      int indx = c.mcParticleIds()[0]; // first alwais exist
      for (uint32_t iii = 0; iii < edep.size(); iii++) {
        if (edep[iii] > 0) {
          if (ed == 0) { // first nontrivial parent
            ed = edep[iii];
            indx = c.mcParticleIds()[iii];
          } else {
            if (indx < c.mcParticleIds()[iii]) { // this might be parent? then take daughter
              ed = edep[iii];
              indx = c.mcParticleIds()[iii];
            }
          }
        }
      }
      int labelIndex = cellTruth.getIndexedSize();
      o2::phos::MCLabel label(indx, 0, 0, (ed == 0.), ed); // MCLabel(Int_t trackID, Int_t eventID, Int_t srcID, bool fake, float edep): fake if deposited energy zero
      cellTruth.addElement(labelIndex, label);
    }
    // Set number of cells in last TrigRec
    if (phosCellTRs.size() > 0) {
      phosCellTRs.back().setNumberOfObjects(phosCells.size() - phosCellTRs.back().getFirstEntry());
    }

    // clusterize
    o2::dataformats::MCTruthContainer<o2::phos::MCLabel> outputTruthCont;
    clusterizerPHOS->processCells(phosCells, phosCellTRs, &cellTruth,
                                  outputPHOSClusters, outputCluElements, outputPHOSClusterTrigRecs, outputTruthCont);

    // Find  CPV clusters corresponding to PHOS trigger records
    std::vector<std::pair<float, float>> cpvMatchPoints[kCpvCells];
    // Number of entries in each cell per TrigRecord
    std::vector<TrackTrigRec> cpvNMatchPoints;
    cpvNMatchPoints.reserve(outputPHOSClusterTrigRecs.size());
    int64_t curBC = -1;
    if (cpvs.begin() != cpvs.end()) {
      curBC = cpvs.begin().bc_as<aod::BCsWithTimestamps>().globalBC();
      cpvNMatchPoints.emplace_back();
      cpvNMatchPoints.back().mTR = curBC;
      for (int i = kCpvCells; i--;) {
        cpvNMatchPoints.back().mStart[i] = 0;
      }
    }

    for (const auto& cpvclu : cpvs) {
      if (static_cast<int64_t>(cpvclu.bc_as<aod::BCsWithTimestamps>().globalBC()) != curBC) { // new BC
        // mark last entry in previous range
        for (int i = kCpvCells; i--;) {
          cpvNMatchPoints.back().mEnd[i] = cpvMatchPoints[i].size();
        }
        curBC = cpvclu.bc_as<aod::BCsWithTimestamps>().globalBC();
        cpvNMatchPoints.back() = cpvNMatchPoints.emplace_back();
        cpvNMatchPoints.back().mTR = curBC;
        for (int i = kCpvCells; i--;) {
          cpvNMatchPoints.back().mStart[i] = cpvMatchPoints[i].size();
        }
      }
      // anyway add coordinates
      if (cpvclu.amplitude() < static_cast<std::vector<double>>(cpvMinE)[static_cast<int>(cpvclu.moduleNumber()) - 2]) {
        continue;
      }
      int index = cpvMatchIndex(cpvclu.moduleNumber(), cpvclu.posX(), cpvclu.posZ());
      cpvMatchPoints[index].emplace_back(cpvclu.posX(), cpvclu.posZ());
    }
    if (cpvNMatchPoints.size() > 0) {
      for (int i = kCpvCells; i--;) {
        cpvNMatchPoints.back().mEnd[i] = cpvMatchPoints[i].size();
      }
    }

    // Fill output
    for (const auto& cluTR : outputPHOSClusterTrigRecs) {
      int firstClusterInEvent = cluTR.getFirstEntry();
      int lastClusterInEvent = firstClusterInEvent + cluTR.getNumberOfObjects();

      // Extract primary vertex
      TVector3 vtx = {0., 0., 0.}; // default, if not collision will be found
      int colId = -1;
      auto coliter = colMap.find(cluTR.getBCData().toLong());
      if (coliter != colMap.end()) { // get vertex from collision
        // find collision corresponding to current BC
        auto clvtx = colls.begin() + coliter->second;
        vtx.SetXYZ(clvtx.posX(), clvtx.posY(), clvtx.posZ());
        colId = coliter->second;
      }

      bool cpvExist = false;
      // find cpvTR for this BC
      auto cpvPoints = cpvNMatchPoints.begin();
      while (cpvPoints != cpvNMatchPoints.end()) {
        if (cpvPoints->mTR == cluTR.getBCData().toLong()) {
          cpvExist = true;
          break;
        }
        cpvPoints++;
      }

      for (int i = firstClusterInEvent; i < lastClusterInEvent; i++) {
        o2::phos::Cluster& clu = outputPHOSClusters[i];
        float e = (useCoreE) ? clu.getCoreEnergy() : clu.getEnergy();
        if (e == 0) {
          continue;
        }
        float posX, posZ;
        clu.getLocalPosition(posX, posZ);

        // Correction for the depth of the shower starting point (TDR p 127)
        const float para = 0.925;
        const float parb = 6.52;
        float depth = para * std::log(e) + parb;
        posX -= posX * depth / 460.;
        posZ -= (posZ - vtx.Z()) * depth / 460.;

        int mod = clu.module();
        TVector3 globaPos;
        geomPHOS->local2Global(mod, posX, posZ, globaPos);

        TVector3 mom = globaPos - vtx;
        if (mom.Mag() == 0) { // should not happpen
          continue;
        }

        e = nonlinearity(e);

        mom.SetMag(e);

        float cpvdist = 99.;
        const float cellSizeX = 2 * kCpvMaxX / kCpvX;
        const float cellSizeZ = 2 * kCpvMaxZ / kCpvZ;
        // look 9 CPV regions around PHOS cluster

        if (mod >= 2 && cpvExist) { // CPV exist in mods 2,3,4
          int phosIndex = cpvMatchIndex(mod, posX, posZ);
          std::vector<int> regions;
          regions.push_back(phosIndex);
          if (posX > -kCpvMaxX + cellSizeX) {
            if (posZ > -kCpvMaxZ + cellSizeZ) { // bottom left
              regions.push_back(phosIndex - kCpvZ - 1);
            }
            regions.push_back(phosIndex - kCpvZ);
            if (posZ < kCpvMaxZ - cellSizeZ) { // top left
              regions.push_back(phosIndex - kCpvZ + 1);
            }
          }
          if (posZ > -kCpvMaxZ + cellSizeZ) { // bottom
            regions.push_back(phosIndex - 1);
          }
          if (posZ < kCpvMaxZ - cellSizeZ) { // top
            regions.push_back(phosIndex + 1);
          }
          if (posX < kCpvMaxX - cellSizeX) {
            if (posZ > -kCpvMaxZ + cellSizeZ) { // bottom right
              regions.push_back(phosIndex + kCpvZ - 1);
            }
            regions.push_back(phosIndex + kCpvZ);
            if (posZ < kCpvMaxZ - cellSizeZ) { // top right
              regions.push_back(phosIndex + kCpvZ + 1);
            }
          }
          float sigmaX = 1. / std::min(5.2, 1.111 + 0.56 * std::exp(-0.031 * e * e) + 4.8 / std::pow(e + 0.61, 3)); // inverse sigma X
          float sigmaZ = 1. / std::min(3.3, 1.12 + 0.35 * std::exp(-0.032 * e * e) + 0.75 / std::pow(e + 0.24, 3)); // inverse sigma Z

          for (const int& indx : regions) {
            if (indx >= 0 && indx < kCpvCells) {
              for (int ii = cpvPoints->mStart[indx]; ii < cpvPoints->mEnd[indx]; ii++) {
                auto p = cpvMatchPoints[indx][ii];
                float d = std::pow((p.first - posX) * sigmaX, 2) + std::pow((p.second - posZ) * sigmaZ, 2);
                if (d < cpvdist) {
                  cpvdist = d;
                }
              }
            }
          }
        }
        if (cpvdist != 99.) {           // was evaluated
          cpvdist = std::sqrt(cpvdist); // was squared
        }
        int cpvindex = -2; // -2 no CPV in event
        if (cpvExist) {
          cpvindex = -1; // there were CPV clusters
        }

        float lambdaShort = 0., lambdaLong = 0.;
        clu.getElipsAxis(lambdaShort, lambdaLong);

        // MC info
        mclabels.clear();
        mcamplitudes.clear();
        gsl::span<const o2::phos::MCLabel> spDigList = outputTruthCont.getLabels(i);
        for (const auto& cellLab : spDigList) {
          mclabels.push_back(cellLab.getTrackID()); // Track ID in current event?
          mcamplitudes.push_back(cellLab.getEdep());
        }
        // Clear Collision assignment
        if (colId == -1) {
          // Ambiguos Collision assignment
          cluambcursor(
            bcMap[cluTR.getBCData().toLong()],
            mom.X(), mom.Y(), mom.Z(), e,
            mod, clu.getMultiplicity(), posX, posZ,
            globaPos.X(), globaPos.Y(), globaPos.Z(),
            clu.getTime(), clu.getNExMax(),
            lambdaLong, lambdaShort,
            cpvdist, cpvindex,
            clu.firedTrigger(),
            clu.getDistanceToBadChannel());
          cluambmccursor(
            mclabels,
            mcamplitudes);

        } else { // Normal collision
          auto col = colls.begin() + colId;
          clucursor(
            col,
            mom.X(), mom.Y(), mom.Z(), e,
            mod, clu.getMultiplicity(), posX, posZ,
            globaPos.X(), globaPos.Y(), globaPos.Z(),
            clu.getTime(), clu.getNExMax(),
            lambdaLong, lambdaShort,
            cpvdist, cpvindex,
            clu.firedTrigger(),
            clu.getDistanceToBadChannel());
          clumccursor(
            mclabels,
            mcamplitudes);
        }
      }
    }
  }

  PROCESS_SWITCH(CaloClusterProducer, processStandaloneMC, "Process MC, PHOS and CPV only", false);

  //------------------------------------------------------------
  void processFull(o2::aod::BCsWithTimestamps const& bcs,
                   o2::aod::Collisions const& colls,
                   o2::aod::Calos const& cells,
                   o2::aod::CaloTriggers const&,
                   o2::aod::CPVClusters const& cpvs,
                   o2::aod::FullTracks const& tracks)
  {

    int64_t timestamp = 0;
    if (bcs.begin() != bcs.end()) {
      timestamp = bcs.begin().timestamp(); // timestamp for CCDB object retrieval
    } else {
      return;
    }

    if (mRunNumber != bcs.begin().runNumber()) {
      mRunNumber = bcs.begin().runNumber();
      LOG(info) << ">>>>>>>>>>>> Current run number: " << mRunNumber;
      o2::parameters::GRPMagField* grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", timestamp);
      if (grpo == nullptr) {
        LOGF(fatal, "Run 3 GRP object (type o2::parameters::GRPMagField) is not available in CCDB for run=%d at timestamp=%llu", mRunNumber, timestamp);
      }
      o2::base::Propagator::initFieldFromGRP(grpo);
    }

    std::map<int64_t, int> bcMap;
    int bcId = 0;
    for (const auto& bc : bcs) {
      bcMap[bc.globalBC()] = bcId;
      bcId++;
    }

    // If several collisions appear in BC, choose one with largers number of contributors
    std::map<int64_t, int> colMap;
    int colId = 0;
    for (const auto& cl : colls) {
      auto colbc = colMap.find(cl.bc_as<aod::BCsWithTimestamps>().globalBC());
      if (colbc == colMap.end()) { // single collision per BC
        colMap[cl.bc_as<aod::BCsWithTimestamps>().globalBC()] = colId;
      } else { // not unique collision per BC
        auto coll2 = colls.begin() + colbc->second;
        if (cl.numContrib() > coll2.numContrib()) {
          colMap[cl.bc_as<aod::BCsWithTimestamps>().globalBC()] = colId;
        }
      }
      colId++;
    }
    // Fill list of cells and cell TrigRecs per TF as an input for clusterizer
    // clusterize
    // Fill output table

    // calibration may be updated by CCDB fetcher
    const o2::phos::BadChannelsMap* badMap = ccdb->getForTimeStamp<o2::phos::BadChannelsMap>(badMapPath, timestamp);
    const o2::phos::CalibParams* calibParams = ccdb->getForTimeStamp<o2::phos::CalibParams>(calibPath, timestamp);

    if (!isMC && !skipL1phase) {
      const std::vector<int>* vec = ccdb->getForTimeStamp<std::vector<int>>(l1phasePath, timestamp);
      if (vec) {
        clusterizerPHOS->setL1phase((*vec)[0]);
      } else {
        LOG(fatal) << "Can not get PHOS L1phase calibration";
      }
    }

    if (badMap) {
      clusterizerPHOS->setBadMap(badMap);
    } else {
      LOG(fatal) << "Can not get PHOS Bad Map";
    }
    if (calibParams) {
      clusterizerPHOS->setCalibration(calibParams);
    } else {
      LOG(fatal) << "Can not get PHOS calibration";
    }

    phosCells.clear();
    phosCells.reserve(cells.size());
    phosCellTRs.clear();
    phosCellTRs.reserve(bcs.size());
    outputCluElements.clear();
    outputPHOSClusters.clear();
    outputPHOSClusterTrigRecs.clear();

    o2::InteractionRecord ir;
    const int kPHOS = 0;
    for (const auto& c : cells) {
      if (c.caloType() != kPHOS) // PHOS
        continue;
      // Fix for bug in trigger digits
      if ((c.cellType() == phos::TRU2x2 || c.cellType() == phos::TRU4x4) && c.cellNumber() == 0) {
        continue;
      }
      if (phosCellTRs.size() == 0) { // first cell, first TrigRec
        ir.setFromLong(c.bc_as<aod::BCsWithTimestamps>().globalBC());
        phosCellTRs.emplace_back(ir, 0, 0); // BC,first cell, ncells
      }
      if (static_cast<uint64_t>(phosCellTRs.back().getBCData().toLong()) != c.bc_as<aod::BCsWithTimestamps>().globalBC()) { // switch to new BC
        // switch to another BC: set size and create next TriRec
        phosCellTRs.back().setNumberOfObjects(phosCells.size() - phosCellTRs.back().getFirstEntry());
        // Next event/trig rec.
        ir.setFromLong(c.bc_as<aod::BCsWithTimestamps>().globalBC());
        phosCellTRs.emplace_back(ir, phosCells.size(), 0);
      }
      phosCells.emplace_back(c.cellNumber(), c.amplitude(), c.time(),
                             static_cast<o2::phos::ChannelType_t>(c.cellType()));
    }
    // Set number of cells in last TrigRec
    if (phosCellTRs.size() > 0) {
      phosCellTRs.back().setNumberOfObjects(phosCells.size() - phosCellTRs.back().getFirstEntry());
    }

    // clusterize
    o2::dataformats::MCTruthContainer<o2::phos::MCLabel> dummyMC;
    clusterizerPHOS->processCells(phosCells, phosCellTRs, nullptr,
                                  outputPHOSClusters, outputCluElements, outputPHOSClusterTrigRecs, dummyMC);

    // Find  CPV clusters corresponding to PHOS trigger records
    std::vector<std::pair<float, float>> cpvMatchPoints[kCpvCells];
    // Number of entries in each cell per TrigRecord
    std::vector<TrackTrigRec> cpvNMatchPoints;
    cpvNMatchPoints.reserve(outputPHOSClusterTrigRecs.size());

    int64_t curBC = -1;
    if (cpvs.begin() != cpvs.end()) {
      curBC = cpvs.begin().bc_as<aod::BCsWithTimestamps>().globalBC();
      cpvNMatchPoints.emplace_back();
      cpvNMatchPoints.back().mTR = curBC;
      for (int i = kCpvCells; i--;) {
        cpvNMatchPoints.back().mStart[i] = 0;
      }
    }

    for (const auto& cpvclu : cpvs) {
      if (static_cast<int64_t>(cpvclu.bc_as<aod::BCsWithTimestamps>().globalBC()) != curBC) { // new BC
        // mark last entry in previous range
        for (int i = kCpvCells; i--;) {
          cpvNMatchPoints.back().mEnd[i] = cpvMatchPoints[i].size();
        }
        curBC = cpvclu.bc_as<aod::BCsWithTimestamps>().globalBC();
        cpvNMatchPoints.back() = cpvNMatchPoints.emplace_back();
        cpvNMatchPoints.back().mTR = curBC;
        for (int i = kCpvCells; i--;) {
          cpvNMatchPoints.back().mStart[i] = cpvMatchPoints[i].size();
        }
      }
      if (cpvclu.amplitude() < static_cast<std::vector<double>>(cpvMinE)[static_cast<int>(cpvclu.moduleNumber()) - 2]) {
        continue;
      }
      int index = cpvMatchIndex(cpvclu.moduleNumber(), cpvclu.posX(), cpvclu.posZ());
      cpvMatchPoints[index].emplace_back(cpvclu.posX(), cpvclu.posZ());
    }
    if (cpvNMatchPoints.size()) {
      for (int i = kCpvCells; i--;) {
        cpvNMatchPoints.back().mEnd[i] = cpvMatchPoints[i].size();
      }
    }
    // same for tracks
    std::vector<TrackMatch> trackMatchPoints[kCpvCells]; // tracks hit in grid/cell in PHOS
    // Number of entries in each cell per TrigRecord
    std::vector<TrackTrigRec> trackNMatchPoints;
    trackNMatchPoints.reserve(outputPHOSClusterTrigRecs.size());

    curBC = 0;
    for (const auto& track : tracks) {
      if (track.has_collision()) { // ignore orphan tracks without collision
        curBC = track.collision().bc_as<aod::BCsWithTimestamps>().globalBC();
        break;
      }
    }
    bool keepBC = false;
    for (const auto& cluTR : outputPHOSClusterTrigRecs) {
      if (cluTR.getBCData().toLong() == curBC) {
        keepBC = true;
        break;
      }
    }
    if (keepBC) {
      trackNMatchPoints.emplace_back();
      trackNMatchPoints.back().mTR = curBC;
      for (int i = kCpvCells; i--;) {
        trackNMatchPoints.back().mStart[i] = 0;
      }
    }
    for (const auto& track : tracks) {
      if (!track.has_collision()) {
        continue;
      }
      if (static_cast<int64_t>(track.collision().bc_as<aod::BCsWithTimestamps>().globalBC()) != curBC) { // new BC
        // close previous BC if exist
        if (keepBC) {
          // mark last entry in previous range
          for (int i = kCpvCells; i--;) {
            trackNMatchPoints.back().mEnd[i] = trackMatchPoints[i].size();
          }
          curBC = track.collision().bc_as<aod::BCsWithTimestamps>().globalBC();
        }
        keepBC = false;
        for (const auto& cluTR : outputPHOSClusterTrigRecs) {
          if (cluTR.getBCData().toLong() == curBC) {
            keepBC = true;
            break;
          }
        }
        if (!keepBC) {
          continue;
        }
        trackNMatchPoints.emplace_back();
        trackNMatchPoints.back().mTR = curBC;
        for (int i = kCpvCells; i--;) {
          trackNMatchPoints.back().mStart[i] = trackMatchPoints[i].size();
        }
      }
      // if (!keepBC || !track.isGlobalTrack()) {  // only global tracks
      if (!keepBC) {
        continue;
      }
      // calculate coordinate in PHOS plane
      int16_t module;
      float trackX, trackZ;
      auto trackPar = getTrackPar(track);
      if (impactOnPHOS(trackPar, track.trackEtaEmcal(), track.trackPhiEmcal(), track.collision().posZ(), module, trackX, trackZ)) {
        int index = cpvMatchIndex(module, trackX, trackZ);
        trackMatchPoints[index].emplace_back(trackX, trackZ, track.globalIndex());
      }
    }
    if (keepBC) {
      for (int i = kCpvCells; i--;) {
        trackNMatchPoints.back().mEnd[i] = trackMatchPoints[i].size();
      }
    }

    // Fill output tables
    for (const auto& cluTR : outputPHOSClusterTrigRecs) {
      int firstClusterInEvent = cluTR.getFirstEntry();
      int lastClusterInEvent = firstClusterInEvent + cluTR.getNumberOfObjects();

      // Extract primary vertex
      TVector3 vtx = {0., 0., 0.}; // default, if not collision will be found
      int colId = -1;
      auto coliter = colMap.find(cluTR.getBCData().toLong());
      if (coliter != colMap.end()) { // get vertex from collision
        // find collision corresponding to current BC
        auto clvtx = colls.begin() + coliter->second;
        vtx.SetXYZ(clvtx.posX(), clvtx.posY(), clvtx.posZ());
        colId = coliter->second;
      }

      bool cpvExist = false;
      // find cpvTR for this BC
      auto cpvPoints = cpvNMatchPoints.begin();
      while (cpvPoints != cpvNMatchPoints.end()) {
        if (cpvPoints->mTR == cluTR.getBCData().toLong()) {
          cpvExist = true;
          break;
        }
        cpvPoints++;
      }

      // find cpvTR for this BC
      auto trackPoints = trackNMatchPoints.begin();
      while (trackPoints != trackNMatchPoints.end()) {
        if (trackPoints->mTR == cluTR.getBCData().toLong()) {
          break;
        }
        trackPoints++;
      }

      for (int i = firstClusterInEvent; i < lastClusterInEvent; i++) {
        o2::phos::Cluster& clu = outputPHOSClusters[i];
        float e = (useCoreE) ? clu.getCoreEnergy() : clu.getEnergy();
        if (e == 0) {
          continue;
        }
        float posX, posZ;
        clu.getLocalPosition(posX, posZ);

        // Correction for the depth of the shower starting point (TDR p 127)
        const float para = 0.925;
        const float parb = 6.52;
        float depth = para * std::log(e) + parb;
        posX -= posX * depth / 460.;
        posZ -= (posZ - vtx.Z()) * depth / 460.;

        int mod = clu.module();
        TVector3 globaPos;
        geomPHOS->local2Global(mod, posX, posZ, globaPos);

        TVector3 mom = globaPos - vtx;
        if (mom.Mag() == 0) { // should not happpen
          continue;
        }

        e = nonlinearity(e);

        mom.SetMag(e);

        // CPV and track match
        const float cellSizeX = 2 * kCpvMaxX / kCpvX;
        const float cellSizeZ = 2 * kCpvMaxZ / kCpvZ;
        // look 9 CPV regions around PHOS cluster
        int phosIndex = cpvMatchIndex(mod, posX, posZ);
        std::vector<int> regions;
        regions.push_back(phosIndex);
        if (posX > -kCpvMaxX + cellSizeX) {
          if (posZ > -kCpvMaxZ + cellSizeZ) { // bottom left
            regions.push_back(phosIndex - kCpvZ - 1);
          }
          regions.push_back(phosIndex - kCpvZ);
          if (posZ < kCpvMaxZ - cellSizeZ) { // top left
            regions.push_back(phosIndex - kCpvZ + 1);
          }
        }
        if (posZ > -kCpvMaxZ + cellSizeZ) { // bottom
          regions.push_back(phosIndex - 1);
        }
        if (posZ < kCpvMaxZ - cellSizeZ) { // top
          regions.push_back(phosIndex + 1);
        }
        if (posX < kCpvMaxX - cellSizeX) {
          if (posZ > -kCpvMaxZ + cellSizeZ) { // bottom right
            regions.push_back(phosIndex + kCpvZ - 1);
          }
          regions.push_back(phosIndex + kCpvZ);
          if (posZ < kCpvMaxZ - cellSizeZ) { // top right
            regions.push_back(phosIndex + kCpvZ + 1);
          }
        }
        float sigmaX = 1. / std::min(5.2, 1.111 + 0.56 * std::exp(-0.031 * e * e) + 4.8 / std::pow(e + 0.61, 3)); // inverse sigma X
        float sigmaZ = 1. / std::min(3.3, 1.12 + 0.35 * std::exp(-0.032 * e * e) + 0.75 / std::pow(e + 0.24, 3)); // inverse sigma Z
        float cpvdist = 99., trackdist = 99.;
        // float cpvDx = 0., cpvDz = 0.;
        float trackDx = 9999., trackDz = 9999.;
        int trackindex = -1;
        for (const int& indx : regions) {
          if (cpvPoints != cpvNMatchPoints.end()) {
            if (indx >= 0 && indx < kCpvCells) {
              for (int ii = cpvPoints->mStart[indx]; ii < cpvPoints->mEnd[indx]; ii++) {
                auto p = cpvMatchPoints[indx][ii];
                float d = std::pow((p.first - posX) * sigmaX, 2) + std::pow((p.second - posZ) * sigmaZ, 2);
                if (d < cpvdist) {
                  cpvdist = d;
                }
              }
            }
          }

          // same for tracks
          if (trackPoints != trackNMatchPoints.end()) {
            for (int ii = trackPoints->mStart[indx]; ii < trackPoints->mEnd[indx]; ii++) {
              auto pp = trackMatchPoints[indx][ii];
              float d = std::pow((pp.pX - posX) * sigmaX, 2) + std::pow((pp.pZ - posZ) * sigmaZ, 2); // TODO different sigma for tracks
              if (d < trackdist) {
                trackdist = d;
                trackDx = pp.pX - posX;
                trackDz = pp.pZ - posZ;
                trackindex = pp.indx;
              }
            }
          }
        }

        if (cpvdist != 99.) {           // was evaluated
          cpvdist = std::sqrt(cpvdist); // was squared
        }
        if (trackdist != 99.) {             // was evaluated
          trackdist = std::sqrt(trackdist); // was squared
        }

        float lambdaShort = 0., lambdaLong = 0.;
        clu.getElipsAxis(lambdaShort, lambdaLong);

        int cpvindex = -2; // -2 no CPV in event
        if (cpvExist) {
          cpvindex = -1; // there were CPV clusters
        }
        if (colId == -1) {
          // Ambiguos Collision assignment
          cluambcursor(
            bcMap[cluTR.getBCData().toLong()],
            mom.X(), mom.Y(), mom.Z(), e,
            mod, clu.getMultiplicity(), posX, posZ,
            globaPos.X(), globaPos.Y(), globaPos.Z(),
            clu.getTime(), clu.getNExMax(),
            lambdaLong, lambdaShort,
            cpvdist, cpvindex,
            clu.firedTrigger(),
            clu.getDistanceToBadChannel());

        } else { // Normal collision
          auto col = colls.begin() + colId;
          clucursor(
            col,
            mom.X(), mom.Y(), mom.Z(), e,
            mod, clu.getMultiplicity(), posX, posZ,
            globaPos.X(), globaPos.Y(), globaPos.Z(),
            clu.getTime(), clu.getNExMax(),
            lambdaLong, lambdaShort,
            cpvdist, cpvindex,
            clu.firedTrigger(),
            clu.getDistanceToBadChannel());

          matchedTracks(clucursor.lastIndex(), trackindex, trackDx, trackDz);
        }
      }
    }
  }

  PROCESS_SWITCH(CaloClusterProducer, processFull, "Process with track matching", false);

  //------------------------------------------------------------
  void processFullMC(o2::aod::BCsWithTimestamps const& bcs,
                     o2::aod::Collisions const& colls,
                     McCells const& cells,
                     o2::aod::CaloTriggers const&,
                     o2::aod::CPVClusters const& cpvs,
                     o2::aod::FullTracks const& tracks)
  {
    int64_t timestamp = 0;
    if (bcs.begin() != bcs.end()) {
      timestamp = bcs.begin().timestamp(); // timestamp for CCDB object retrieval
    } else {
      return;
    }

    if (mRunNumber != bcs.begin().runNumber()) {
      mRunNumber = bcs.begin().runNumber();
      LOG(info) << ">>>>>>>>>>>> Current run number: " << mRunNumber;
      o2::parameters::GRPMagField* grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", timestamp);
      if (grpo == nullptr) {
        LOGF(fatal, "Run 3 GRP object (type o2::parameters::GRPMagField) is not available in CCDB for run=%d at timestamp=%llu", mRunNumber, timestamp);
      }
      o2::base::Propagator::initFieldFromGRP(grpo);
    }

    std::map<int64_t, int> bcMap;
    int bcId = 0;
    for (const auto& bc : bcs) {
      bcMap[bc.globalBC()] = bcId;
      bcId++;
    }

    // If several collisions appear in BC, choose one with largers number of contributors
    std::map<int64_t, int> colMap;
    int colId = 0;
    for (const auto& cl : colls) {
      auto colbc = colMap.find(cl.bc_as<aod::BCsWithTimestamps>().globalBC());
      if (colbc == colMap.end()) { // single collision per BC
        colMap[cl.bc_as<aod::BCsWithTimestamps>().globalBC()] = colId;
      } else { // not unique collision per BC
        auto coll2 = colls.begin() + colbc->second;
        if (cl.numContrib() > coll2.numContrib()) {
          colMap[cl.bc_as<aod::BCsWithTimestamps>().globalBC()] = colId;
        }
      }
      colId++;
    }
    // Fill list of cells and cell TrigRecs per TF as an input for clusterizer
    // clusterize
    // Fill output table

    // calibration may be updated by CCDB fetcher
    const o2::phos::BadChannelsMap* badMap = ccdb->getForTimeStamp<o2::phos::BadChannelsMap>(badMapPath, timestamp);
    const o2::phos::CalibParams* calibParams = ccdb->getForTimeStamp<o2::phos::CalibParams>(calibPath, timestamp);

    if (!isMC && !skipL1phase) {
      const std::vector<int>* vec = ccdb->getForTimeStamp<std::vector<int>>(l1phasePath, timestamp);
      if (vec) {
        clusterizerPHOS->setL1phase((*vec)[0]);
      } else {
        LOG(fatal) << "Can not get PHOS L1phase calibration";
      }
    }

    if (badMap) {
      clusterizerPHOS->setBadMap(badMap);
    } else {
      LOG(fatal) << "Can not get PHOS Bad Map";
    }
    if (calibParams) {
      clusterizerPHOS->setCalibration(calibParams);
    } else {
      LOG(fatal) << "Can not get PHOS calibration";
    }

    phosCells.clear();
    phosCells.reserve(cells.size());
    phosCellTRs.clear();
    phosCellTRs.reserve(bcs.size());
    outputCluElements.clear();
    outputPHOSClusters.clear();
    outputPHOSClusterTrigRecs.clear();
    o2::dataformats::MCTruthContainer<o2::phos::MCLabel> cellTruth;

    o2::InteractionRecord ir;
    const int kPHOS = 0;
    for (const auto& c : cells) {
      if (c.caloType() != kPHOS) // PHOS
        continue;
      // Fix for bug in trigger digits
      if ((c.cellType() == phos::TRU2x2 || c.cellType() == phos::TRU4x4) && c.cellNumber() == 0) {
        continue;
      }
      if (phosCellTRs.size() == 0) { // first cell, first TrigRec
        ir.setFromLong(c.bc_as<aod::BCsWithTimestamps>().globalBC());
        phosCellTRs.emplace_back(ir, 0, 0); // BC,first cell, ncells
      }
      if (static_cast<uint64_t>(phosCellTRs.back().getBCData().toLong()) != c.bc_as<aod::BCsWithTimestamps>().globalBC()) { // switch to new BC
        // switch to another BC: set size and create next TriRec
        phosCellTRs.back().setNumberOfObjects(phosCells.size() - phosCellTRs.back().getFirstEntry());
        // Next event/trig rec.
        ir.setFromLong(c.bc_as<aod::BCsWithTimestamps>().globalBC());
        phosCellTRs.emplace_back(ir, phosCells.size(), 0);
      }
      phosCells.emplace_back(c.cellNumber(), c.amplitude(), c.time(),
                             static_cast<o2::phos::ChannelType_t>(c.cellType()));

      // process MC info
      auto edep = c.amplitudeA();
      // find particle with non-zero E deposited fraction; if several, use one with largest label (may be daughter - to keep full history)
      float ed = 0;
      int indx = c.mcParticleIds()[0]; // first alwais exist
      for (uint32_t iii = 0; iii < edep.size(); iii++) {
        if (edep[iii] > 0) {
          if (ed == 0) { // first nontrivial parent
            ed = edep[iii];
            indx = c.mcParticleIds()[iii];
          } else {
            if (indx < c.mcParticleIds()[iii]) { // this might be parent? then take daughter
              ed = edep[iii];
              indx = c.mcParticleIds()[iii];
            }
          }
        }
      }
      int labelIndex = cellTruth.getIndexedSize();
      o2::phos::MCLabel label(indx, 0, 0, (ed == 0.), ed); // MCLabel(Int_t trackID, Int_t eventID, Int_t srcID, bool fake, float edep): fake if deposited energy zero
      cellTruth.addElement(labelIndex, label);
    }
    // Set number of cells in last TrigRec
    if (phosCellTRs.size() > 0) {
      phosCellTRs.back().setNumberOfObjects(phosCells.size() - phosCellTRs.back().getFirstEntry());
    }
    // clusterize
    o2::dataformats::MCTruthContainer<o2::phos::MCLabel> outputTruthCont;
    clusterizerPHOS->processCells(phosCells, phosCellTRs, &cellTruth,
                                  outputPHOSClusters, outputCluElements, outputPHOSClusterTrigRecs, outputTruthCont);

    // Find  CPV clusters corresponding to PHOS trigger records
    std::vector<std::pair<float, float>> cpvMatchPoints[kCpvCells];
    // Number of entries in each cell per TrigRecord
    std::vector<TrackTrigRec> cpvNMatchPoints;
    cpvNMatchPoints.reserve(outputPHOSClusterTrigRecs.size());

    int64_t curBC = -1;
    if (cpvs.begin() != cpvs.end()) {
      curBC = cpvs.begin().bc_as<aod::BCsWithTimestamps>().globalBC();
      cpvNMatchPoints.emplace_back();
      cpvNMatchPoints.back().mTR = curBC;
      for (int i = kCpvCells; i--;) {
        cpvNMatchPoints.back().mStart[i] = 0;
      }
    }

    for (const auto& cpvclu : cpvs) {
      if (static_cast<int64_t>(cpvclu.bc_as<aod::BCsWithTimestamps>().globalBC()) != curBC) { // new BC
        // mark last entry in previous range
        for (int i = kCpvCells; i--;) {
          cpvNMatchPoints.back().mEnd[i] = cpvMatchPoints[i].size();
        }
        curBC = cpvclu.bc_as<aod::BCsWithTimestamps>().globalBC();
        cpvNMatchPoints.emplace_back();
        cpvNMatchPoints.back().mTR = curBC;
        for (int i = kCpvCells; i--;) {
          cpvNMatchPoints.back().mStart[i] = cpvMatchPoints[i].size();
        }
      }
      if (cpvclu.amplitude() < static_cast<std::vector<double>>(cpvMinE)[static_cast<int>(cpvclu.moduleNumber()) - 2]) {
        continue;
      }
      int index = cpvMatchIndex(cpvclu.moduleNumber(), cpvclu.posX(), cpvclu.posZ());
      cpvMatchPoints[index].emplace_back(cpvclu.posX(), cpvclu.posZ());
    }
    if (cpvNMatchPoints.size()) {
      for (int i = kCpvCells; i--;) {
        cpvNMatchPoints.back().mEnd[i] = cpvMatchPoints[i].size();
      }
    }
    // same for tracks
    std::vector<TrackMatch> trackMatchPoints[kCpvCells]; // tracks hit in grid/cell in PHOS
    // Number of entries in each cell per TrigRecord
    std::vector<TrackTrigRec> trackNMatchPoints;
    trackNMatchPoints.reserve(outputPHOSClusterTrigRecs.size());

    curBC = -1;
    for (const auto& track : tracks) {
      if (track.has_collision()) { // ignore orphan tracks without collision
        curBC = track.collision().bc_as<aod::BCsWithTimestamps>().globalBC();
        break;
      }
    }
    bool keepBC = false;
    for (const auto& cluTR : outputPHOSClusterTrigRecs) {
      if (cluTR.getBCData().toLong() == curBC) {
        keepBC = true;
        break;
      }
    }
    if (keepBC) {
      trackNMatchPoints.emplace_back();
      trackNMatchPoints.back().mTR = curBC;
      for (int i = kCpvCells; i--;) {
        trackNMatchPoints.back().mStart[i] = 0;
      }
    }
    for (const auto& track : tracks) {
      if (!track.has_collision()) {
        continue;
      }
      if (static_cast<int64_t>(track.collision().bc_as<aod::BCsWithTimestamps>().globalBC()) != curBC) { // new BC
        // close previous BC if exist
        if (keepBC) {
          // mark last entry in previous range
          for (int i = kCpvCells; i--;) {
            trackNMatchPoints.back().mEnd[i] = trackMatchPoints[i].size();
          }
          curBC = track.collision().bc_as<aod::BCsWithTimestamps>().globalBC();
        }
        keepBC = false;
        for (const auto& cluTR : outputPHOSClusterTrigRecs) {
          if (cluTR.getBCData().toLong() == curBC) {
            keepBC = true;
            break;
          }
        }
        if (!keepBC) {
          continue;
        }
        trackNMatchPoints.emplace_back();
        trackNMatchPoints.back().mTR = curBC;
        for (int i = kCpvCells; i--;) {
          trackNMatchPoints.back().mStart[i] = trackMatchPoints[i].size();
        }
      }
      // if (!keepBC || !track.isGlobalTrack()) {  // only global tracks
      if (!keepBC) {
        continue;
      }
      // calculate coordinate in PHOS plane
      int16_t module;
      float trackX, trackZ;
      auto trackPar = getTrackPar(track);
      if (impactOnPHOS(trackPar, track.trackEtaEmcal(), track.trackPhiEmcal(), track.collision().posZ(), module, trackX, trackZ)) {
        int index = cpvMatchIndex(module, trackX, trackZ);
        trackMatchPoints[index].emplace_back(trackX, trackZ, track.globalIndex());
      }
    }
    if (keepBC) {
      for (int i = kCpvCells; i--;) {
        trackNMatchPoints.back().mEnd[i] = trackMatchPoints[i].size();
      }
    }

    // Fill output tables
    for (const auto& cluTR : outputPHOSClusterTrigRecs) {
      int firstClusterInEvent = cluTR.getFirstEntry();
      int lastClusterInEvent = firstClusterInEvent + cluTR.getNumberOfObjects();

      // Extract primary vertex
      TVector3 vtx = {0., 0., 0.}; // default, if not collision will be found
      int colId = -1;
      auto coliter = colMap.find(cluTR.getBCData().toLong());
      if (coliter != colMap.end()) { // get vertex from collision
        // find collision corresponding to current BC
        auto clvtx = colls.begin() + coliter->second;
        vtx.SetXYZ(clvtx.posX(), clvtx.posY(), clvtx.posZ());
        colId = coliter->second;
      }

      bool cpvExist = false;
      // find cpvTR for this BC
      auto cpvPoints = cpvNMatchPoints.begin();
      while (cpvPoints != cpvNMatchPoints.end()) {
        if (cpvPoints->mTR == cluTR.getBCData().toLong()) {
          cpvExist = true;
          break;
        }
        cpvPoints++;
      }
      // find cpvTR for this BC
      auto trackPoints = trackNMatchPoints.begin();
      while (trackPoints != trackNMatchPoints.end()) {
        if (trackPoints->mTR == cluTR.getBCData().toLong()) {
          break;
        }
        trackPoints++;
      }

      for (int i = firstClusterInEvent; i < lastClusterInEvent; i++) {
        o2::phos::Cluster& clu = outputPHOSClusters[i];
        float e = (useCoreE) ? clu.getCoreEnergy() : clu.getEnergy();
        if (e == 0) {
          continue;
        }
        float posX, posZ;
        clu.getLocalPosition(posX, posZ);

        // Correction for the depth of the shower starting point (TDR p 127)
        const float para = 0.925;
        const float parb = 6.52;
        float depth = para * std::log(e) + parb;
        posX -= posX * depth / 460.;
        posZ -= (posZ - vtx.Z()) * depth / 460.;

        int mod = clu.module();
        TVector3 globaPos;
        geomPHOS->local2Global(mod, posX, posZ, globaPos);

        TVector3 mom = globaPos - vtx;
        if (mom.Mag() == 0) { // should not happpen
          continue;
        }

        e = nonlinearity(e);

        mom.SetMag(e);
        // CPV and track match
        const float cellSizeX = 2 * kCpvMaxX / kCpvX;
        const float cellSizeZ = 2 * kCpvMaxZ / kCpvZ;
        // look 9 CPV regions around PHOS cluster
        int phosIndex = cpvMatchIndex(mod, posX, posZ);
        std::vector<int> regions;
        regions.push_back(phosIndex);
        if (posX > -kCpvMaxX + cellSizeX) {
          if (posZ > -kCpvMaxZ + cellSizeZ) { // bottom left
            regions.push_back(phosIndex - kCpvZ - 1);
          }
          regions.push_back(phosIndex - kCpvZ);
          if (posZ < kCpvMaxZ - cellSizeZ) { // top left
            regions.push_back(phosIndex - kCpvZ + 1);
          }
        }
        if (posZ > -kCpvMaxZ + cellSizeZ) { // bottom
          regions.push_back(phosIndex - 1);
        }
        if (posZ < kCpvMaxZ - cellSizeZ) { // top
          regions.push_back(phosIndex + 1);
        }
        if (posX < kCpvMaxX - cellSizeX) {
          if (posZ > -kCpvMaxZ + cellSizeZ) { // bottom right
            regions.push_back(phosIndex + kCpvZ - 1);
          }
          regions.push_back(phosIndex + kCpvZ);
          if (posZ < kCpvMaxZ - cellSizeZ) { // top right
            regions.push_back(phosIndex + kCpvZ + 1);
          }
        }
        float sigmaX = 1. / std::min(5.2, 1.111 + 0.56 * std::exp(-0.031 * e * e) + 4.8 / std::pow(e + 0.61, 3)); // inverse sigma X
        float sigmaZ = 1. / std::min(3.3, 1.12 + 0.35 * std::exp(-0.032 * e * e) + 0.75 / std::pow(e + 0.24, 3)); // inverse sigma Z
        float cpvdist = 99., trackdist = 99.;
        // float cpvDx = 0., cpvDz = 0.;
        float trackDx = 9999., trackDz = 9999.;
        int trackindex = -1;
        for (const int& indx : regions) {
          if (cpvPoints != cpvNMatchPoints.end()) {
            if (indx >= 0 && indx < kCpvCells) {
              for (int ii = cpvPoints->mStart[indx]; ii < cpvPoints->mEnd[indx]; ii++) {
                auto p = cpvMatchPoints[indx][ii];
                float d = std::pow((p.first - posX) * sigmaX, 2) + std::pow((p.second - posZ) * sigmaZ, 2);
                if (d < cpvdist) {
                  cpvdist = d;
                }
              }
            }
          }
          // same for tracks
          if (trackPoints != trackNMatchPoints.end()) {
            for (int ii = trackPoints->mStart[indx]; ii < trackPoints->mEnd[indx]; ii++) {
              auto pp = trackMatchPoints[indx][ii];
              float d = std::pow((pp.pX - posX) * sigmaX, 2) + std::pow((pp.pZ - posZ) * sigmaZ, 2); // TODO different sigma for tracks
              if (d < trackdist) {
                trackdist = d;
                trackDx = pp.pX - posX;
                trackDz = pp.pZ - posZ;
                trackindex = pp.indx;
              }
            }
          }
        }

        if (cpvdist != 99.) {           // was evaluated
          cpvdist = std::sqrt(cpvdist); // was squared
        }
        if (trackdist != 99.) {             // was evaluated
          trackdist = std::sqrt(trackdist); // was squared
        }

        float lambdaShort = 0., lambdaLong = 0.;
        clu.getElipsAxis(lambdaShort, lambdaLong);

        int cpvindex = -2; // -2 no CPV in event
        if (cpvExist) {
          cpvindex = -1; // there were CPV clusters
        }
        // MC info
        mclabels.clear();
        mcamplitudes.clear();
        gsl::span<const o2::phos::MCLabel> spDigList = outputTruthCont.getLabels(i);
        for (const auto& cellLab : spDigList) {
          mclabels.push_back(cellLab.getTrackID()); // Track ID in current event?
          mcamplitudes.push_back(cellLab.getEdep());
        }
        if (colId == -1) {
          // Ambiguos Collision assignment
          cluambcursor(
            bcMap[cluTR.getBCData().toLong()],
            mom.X(), mom.Y(), mom.Z(), e,
            mod, clu.getMultiplicity(), posX, posZ,
            globaPos.X(), globaPos.Y(), globaPos.Z(),
            clu.getTime(), clu.getNExMax(),
            lambdaLong, lambdaShort,
            cpvdist, cpvindex,
            clu.firedTrigger(),
            clu.getDistanceToBadChannel());
          cluambmccursor(
            mclabels,
            mcamplitudes);
        } else { // Normal collision
          auto col = colls.begin() + colId;
          clucursor(
            col,
            mom.X(), mom.Y(), mom.Z(), e,
            mod, clu.getMultiplicity(), posX, posZ,
            globaPos.X(), globaPos.Y(), globaPos.Z(),
            clu.getTime(), clu.getNExMax(),
            lambdaLong, lambdaShort,
            cpvdist, cpvindex,
            clu.firedTrigger(),
            clu.getDistanceToBadChannel());

          matchedTracks(clucursor.lastIndex(), trackindex, trackDx, trackDz);
          clumccursor(
            mclabels,
            mcamplitudes);
        }
      }
    }
  }

  PROCESS_SWITCH(CaloClusterProducer, processFullMC, "Process MC with track matching", false);

  int cpvMatchIndex(int16_t module, float x, float z)
  {
    // calculate cell index in grid over PHOS detector
    const float cellSizeX = 2 * kCpvMaxX / kCpvX;
    const float cellSizeZ = 2 * kCpvMaxZ / kCpvZ;
    // in track matching tracks can be beyond CPV surface
    // assign these tracks to the closest cell
    int ix = std::max(0, static_cast<int>((x + kCpvMaxX) / cellSizeX));
    int iz = std::max(0, static_cast<int>((z + kCpvMaxZ) / cellSizeZ));
    if (ix >= kCpvX) {
      ix = kCpvX - 1;
    }
    if (iz >= kCpvZ) {
      iz = kCpvZ - 1;
    }
    return (module - 1) * kCpvX * kCpvZ + ix * kCpvZ + iz; // modules: 1,2,3,4
  }

  bool impactOnPHOS(o2::track::TrackParametrization<float>& trackPar, float trackEta, float trackPhi, float /*zvtx*/, int16_t& module, float& trackX, float& trackZ)
  {
    // eta,phi was calculated at EMCAL radius.
    // Extrapolate to PHOS assuming zeroB and current vertex
    // Check if direction in PHOS acceptance+20cm and return phos module number and coordinates in PHOS module plane
    const float phiMin = 240. * 0.017453293; // PHOS phi coverage-10degree *degToRad
    const float phiMax = 323. * 0.017453293; // PHOS+20 cm * degToRad
    const float xmax = 85.;                  // Maximal x distance per module
    const float etaMax = 0.178266;
    double bz = o2::base::Propagator::Instance()->getNominalBz(); // magnetic field

    trackPhi = RecoDecay::constrainAngle(trackPhi, 0., 1);                       // constrain angle to range 0,twoPi
    if (trackPhi < phiMin || trackPhi > phiMax || std::abs(trackEta) > etaMax) { // do not match even approximately
      return false;
    }

    const float dphi = 20. * 0.017453293;
    module = 1 + static_cast<int16_t>((trackPhi - phiMin) / dphi);
    if (module < 1) {
      module = 1;
    }
    if (module > 4) {
      module = 4;
    }

    // get PHOS radius
    const double shiftY = -1.26;       // Depth-optimized
    double posL[3] = {0., 0., shiftY}; // local position at the center of module
    double posG[3] = {0};
    geomPHOS->getAlignmentMatrix(module)->LocalToMaster(posL, posG);
    double rPHOS = std::sqrt(posG[0] * posG[0] + posG[1] * posG[1]);
    double alpha = (230. + 20. * module) * 0.017453293;

    // During main reconstruction track was propagated to radius 460 cm with accounting material
    // now material is not available. Therefore, start from main rec. position and extrapoate to actual radius without material
    // Get track parameters at point where main reconstruction stop
    float xPHOS = 460.f, xtrg = 0.f;
    if (!trackPar.getXatLabR(xPHOS, xtrg, bz, o2::track::DirType::DirOutward)) {
      return false;
    }
    auto prop = o2::base::Propagator::Instance();
    if (!trackPar.rotate(alpha) ||
        !prop->PropagateToXBxByBz(trackPar, xtrg, 0.95, 10, o2::base::Propagator::MatCorrType::USEMatCorrNONE)) {
      return false;
    }
    // calculate xyz from old (Phi, eta) and new r
    float r = std::sqrt(trackPar.getX() * trackPar.getX() + trackPar.getY() * trackPar.getY());
    trackPar.setX(r * std::cos(trackPhi - alpha));
    trackPar.setY(r * std::sin(trackPhi - alpha));
    trackPar.setZ(r / std::tan(2. * std::atan(std::exp(-trackEta))));

    if (!prop->PropagateToXBxByBz(trackPar, rPHOS, 0.95, 10, o2::base::Propagator::MatCorrType::USEMatCorrNONE)) {
      return false;
    }
    alpha = trackPar.getAlpha();
    double ca = std::cos(alpha), sa = std::sin(alpha);
    posG[0] = trackPar.getX() * ca - trackPar.getY() * sa;
    posG[1] = trackPar.getY() * ca + trackPar.getX() * sa;
    posG[2] = trackPar.getZ();

    geomPHOS->getAlignmentMatrix(module)->MasterToLocal(posG, posL);
    trackX = posL[0];
    trackZ = posL[1];
    // If trackX beyond the module, switch to the next one
    if (std::abs(trackX) < xmax || (trackX < -xmax && module == 1) || (trackX > xmax && module == 4)) {
      return true;
    }
    // re-do extrapolation to correct module
    if (trackX < xmax && module > 1) {
      --module;
    } else {
      if (trackX > xmax && module < 4) {
        ++module;
      }
    }
    // repeat extrapolation for correct module
    alpha = (230. + 20. * module) * 0.017453293;
    posL[0] = 0.;
    posL[1] = 0.;
    posL[2] = shiftY; // local position at the center of module
    geomPHOS->getAlignmentMatrix(module)->LocalToMaster(posL, posG);
    rPHOS = std::sqrt(posG[0] * posG[0] + posG[1] * posG[1]);

    if (!trackPar.rotate(alpha) ||
        !prop->PropagateToXBxByBz(trackPar, xtrg, 0.95, 10, o2::base::Propagator::MatCorrType::USEMatCorrNONE)) {
      return false;
    }
    // calculate xyz from old (Phi, eta) and new r
    r = std::sqrt(trackPar.getX() * trackPar.getX() + trackPar.getY() * trackPar.getY());
    trackPar.setX(r * std::cos(trackPhi - alpha));
    trackPar.setY(r * std::sin(trackPhi - alpha));
    trackPar.setZ(r / std::tan(2. * std::atan(std::exp(-trackEta))));

    if (!prop->PropagateToXBxByBz(trackPar, rPHOS, 0.95, 10, o2::base::Propagator::MatCorrType::USEMatCorrNONE)) {
      return false;
    }
    alpha = trackPar.getAlpha();
    ca = std::cos(alpha);
    sa = std::sin(alpha);
    posG[0] = trackPar.getX() * ca - trackPar.getY() * sa;
    posG[1] = trackPar.getY() * ca + trackPar.getX() * sa;
    posG[2] = trackPar.getZ();

    geomPHOS->getAlignmentMatrix(module)->MasterToLocal(posG, posL);
    trackX = posL[0];
    trackZ = posL[1];
    return true;
  }

  float nonlinearity(float en)
  {
    // Correct for non-linearity
    switch (nonlinType) {
      case 0:
        return en;
      case 1: { // Data Run3
        const double a = 9.3494e-01;
        const double b = 1.00526e-02;
        const double c = 8.45164e-02;
        const double d = -1.03364e-02;
        const double f = 5.4803e-03;
        const double g = 0.779983;
        const double h = 0.622282;
        const double k = 8.0182e-05;
        double eMin = std::max(static_cast<float>(0.1), en); // Parameterization valid down to 100 MeV
        return en * (a + b * eMin + c / eMin + d / (eMin * eMin) + f / ((eMin - g) * (eMin - g) + h * h) + k / std::pow(eMin, 4));
      }
      case 2: { // MC
        const double a = 1.14875;
        const double b = -1.24286e-04;
        const double c = -0.0498217;
        const double d = -0.00215362;
        const double f = 0.886539;
        const double g = -1.98282;
        const double h = 0.0178562;
        const double k = 5.03164e-04;
        double eMin = std::max(static_cast<float>(0.1), en); // Parameterization valid down to 100 MeV
        return en * (a + b * eMin + c / eMin + d / (eMin * eMin) + f / ((eMin - g) * (eMin - g) + h * h) + k / std::pow(eMin, 4));
      }
      case 3: { // Obsolete data Run3
        const double a = 9.34913e-01;
        const double b = 2.33e-03;
        const double c = -8.10e-05;
        const double d = 3.2e-02;
        const double f = -8.0e-03;
        const double g = 1.e-01;
        const double h = 2.e-01;
        const double k = -1.48e-04;
        const double l = 0.194;
        const double m = 0.0025;
        return en * (a + b * en + c * en * en + d / en + f / ((en - g) * (en - g) + h) + k / ((en - l) * (en - l) + m));
      }
      default:
        return en;
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<CaloClusterProducer>(cfgc)};
}
