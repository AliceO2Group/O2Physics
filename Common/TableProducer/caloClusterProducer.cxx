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

#include "Framework/ConfigParamSpec.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/DataModel/CaloClusters.h"

#include "CommonUtils/NameConf.h"
#include "CCDB/BasicCCDBManager.h"
#include "SimulationDataFormat/MCTruthContainer.h"

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

struct caloClusterProducerTask {
  Produces<aod::CaloClusters> clucursor;
  Produces<aod::CaloAmbiguousClusters> cluambcursor;
  Configurable<bool> isMC{"isMC", 0, "0 - data, 1 - MC"};
  Configurable<bool> useCoreE{"coreE", 0, "0 - full energy, 1 - core energy"};
  Configurable<std::vector<double>> cpvMinE{"cpvCluMinAmp", {20., 50., 50.}, "minimal CPV cluster amplitude per module"};

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  std::unique_ptr<o2::phos::Geometry> geomPHOS;
  std::unique_ptr<o2::phos::Clusterer> clusterizerPHOS;
  std::vector<o2::phos::Cell> phosCells;
  std::vector<o2::phos::TriggerRecord> phosCellTRs;
  std::vector<o2::phos::CluElement> outputCluElements;
  std::vector<o2::phos::Cluster> outputPHOSClusters;
  std::vector<o2::phos::TriggerRecord> outputPHOSClusterTrigRecs;
  static constexpr int16_t nCpvX = 7; // grid 13 steps along z and 7 along phi as largest match ellips 20x10 cm
  static constexpr int16_t nCpvZ = 13;
  static constexpr int16_t nCpvCells = 3 * nCpvX * nCpvZ; // 3 modules
  static constexpr float cpvMaxX = 73;                    // max CPV coordinate phi
  static constexpr float cpvMaxZ = 63;                    // max CPV coordinate z

  void init(o2::framework::InitContext&)
  {
    ccdb->setURL(o2::base::NameConf::getCCDBServer());
    geomPHOS = std::make_unique<o2::phos::Geometry>("PHOS");
    clusterizerPHOS = std::make_unique<o2::phos::Clusterer>();
  }

  void process(o2::aod::BCs const& bcs,
               o2::aod::Collisions const& colls,
               o2::aod::Calos const& cells,
               o2::aod::CaloTriggers const& ctrs,
               o2::aod::CPVClusters const& cpvs)
  {

    std::map<int64_t, int> bcMap;
    int bcId = 0;
    for (auto bc : bcs) {
      bcMap[bc.globalBC()] = bcId;
      bcId++;
    }

    // If several collisions appear in BC, choose one with largers number of contributors
    std::map<int64_t, int> colMap;
    int colId = 0;
    for (auto cl : colls) {
      auto colbc = colMap.find(cl.bc().globalBC());
      if (colbc == colMap.end()) { // single collision per BC
        colMap[cl.bc().globalBC()] = colId;
      } else { // not unique collision per BC
        auto coll2 = colls.begin() + colbc->second;
        if (cl.numContrib() > coll2.numContrib()) {
          colMap[cl.bc().globalBC()] = colId;
        }
      }
      colId++;
    }

    // Fill list of cells and cell TrigRecs per TF as an input for clusterizer
    // clusterize
    // Fill output table

    // calibration may be updated by CCDB fetcher
    const o2::phos::BadChannelsMap* badMap = ccdb->get<o2::phos::BadChannelsMap>("PHS/Calib/BadMap");
    const o2::phos::CalibParams* calibParams = ccdb->get<o2::phos::CalibParams>("PHS/Calib/CalibParams");
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
    for (auto& c : cells) {
      if (c.caloType() != kPHOS) // PHOS
        continue;
      if (phosCellTRs.size() == 0) { // first cell, first TrigRec
        ir.setFromLong(c.bc().globalBC());
        phosCellTRs.emplace_back(ir, 0, 0); // BC,first cell, ncells
      }
      if (static_cast<uint64_t>(phosCellTRs.back().getBCData().toLong()) != c.bc().globalBC()) { // switch to new BC
        // switch to another BC: set size and create next TriRec
        phosCellTRs.back().setNumberOfObjects(phosCells.size() - phosCellTRs.back().getFirstEntry());
        // Next event/trig rec.
        ir.setFromLong(c.bc().globalBC());
        phosCellTRs.emplace_back(ir, phosCells.size(), 0);
      }
      phosCells.emplace_back(c.cellNumber(), c.amplitude(), c.time(),
                             static_cast<o2::phos::ChannelType_t>(c.cellType()));
      // TODO process MC info
    }
    // Set number of cells in last TrigRec
    if (phosCellTRs.size() > 0) {
      phosCellTRs.back().setNumberOfObjects(phosCells.size() - phosCellTRs.back().getFirstEntry());
    }

    // clusterize
    if (isMC) {
      o2::dataformats::MCTruthContainer<o2::phos::MCLabel> cellTruth;
      o2::dataformats::MCTruthContainer<o2::phos::MCLabel> outputTruthCont;
      clusterizerPHOS->processCells(phosCells, phosCellTRs, &cellTruth,
                                    outputPHOSClusters, outputCluElements, outputPHOSClusterTrigRecs, outputTruthCont);
    } else {
      o2::dataformats::MCTruthContainer<o2::phos::MCLabel> dummyMC;
      clusterizerPHOS->processCells(phosCells, phosCellTRs, nullptr,
                                    outputPHOSClusters, outputCluElements, outputPHOSClusterTrigRecs, dummyMC);
    }
    // Fill output
    for (auto& cluTR : outputPHOSClusterTrigRecs) {
      int firstClusterInEvent = cluTR.getFirstEntry();
      int lastClusterInEvent = firstClusterInEvent + cluTR.getNumberOfObjects();

      // Prepare arrays with CPV clusters for match for this BC
      // No garantie that order of BCs is the same as in PHOS
      // CPV cluster positions stored within arrays: module:cellx:cellz
      // where cellx=20 cm, cellZ=10 cm
      std::vector<std::pair<float, float>> cpvMatchPoints[nCpvCells];
      // for(int icpv=nCpvCells; icpv--; ){
      //   cpvMatchPoints[icpv].clear();
      // }
      bool scanned = false;
      for (const auto& cpvclu : cpvs) {
        if (cpvclu.bc().globalBC() != static_cast<uint64_t>(cluTR.getBCData().toLong())) {
          if (scanned) {
            break;
          } else {
            continue;
          }
        } else {
          scanned = true;
        }
        if (cpvclu.amplitude() < static_cast<std::vector<double>>(cpvMinE)[static_cast<int>(cpvclu.moduleNumber()) - 2]) {
          continue;
        }
        int index = CpvMatchIndex(cpvclu.moduleNumber(), cpvclu.posX(), cpvclu.posZ());
        cpvMatchPoints[index].emplace_back(cpvclu.posX(), cpvclu.posZ());
      }

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
        float depth = para * TMath::Log(e) + parb;
        posX -= posX * depth / 460.;
        posZ -= (posZ - vtx.Z()) * depth / 460.;

        int mod = clu.module();
        TVector3 globaPos;
        geomPHOS->local2Global(mod, posX, posZ, globaPos);

        TVector3 mom = globaPos - vtx;
        if (mom.Mag() == 0) { // should not happpen
          continue;
        }

        e = Nonlinearity(e);

        mom.SetMag(e);

        // CPV match will be done in independent task
        float trackdist = 99.;
        const float cellSizeX = 2 * cpvMaxX / nCpvX;
        const float cellSizeZ = 2 * cpvMaxZ / nCpvZ;
        int trackindex = -1; // -1:no CPV, -2 CPV in event
        // look 9 CPV regions around PHOS cluster
        if (mod >= 2) { // CPV exist in mods 2,3,4
          int phosIndex = CpvMatchIndex(mod, posX, posZ);
          std::vector<int> regions;
          regions.push_back(phosIndex);
          if (posX > -cpvMaxX + cellSizeX) {
            if (posZ > -cpvMaxZ + cellSizeZ) { // bottom left
              regions.push_back(phosIndex - nCpvZ - 1);
            }
            regions.push_back(phosIndex - nCpvZ);
            if (posZ < cpvMaxZ - cellSizeZ) { // top left
              regions.push_back(phosIndex - nCpvZ + 1);
            }
          }
          if (posZ > -cpvMaxZ + cellSizeZ) { // bottom
            regions.push_back(phosIndex - 1);
          }
          if (posZ < cpvMaxZ - cellSizeZ) { // top
            regions.push_back(phosIndex + 1);
          }
          if (posX < cpvMaxX - cellSizeX) {
            if (posZ > -cpvMaxZ + cellSizeZ) { // bottom right
              regions.push_back(phosIndex + nCpvZ - 1);
            }
            regions.push_back(phosIndex + nCpvZ);
            if (posZ < cpvMaxZ - cellSizeZ) { // top right
              regions.push_back(phosIndex + nCpvZ + 1);
            }
          }
          float sigmaX = 1. / TMath::Min(5.2, 1.111 + 0.56 * TMath::Exp(-0.031 * e * e) + 4.8 / TMath::Power(e + 0.61, 3)); // inverse sigma X
          float sigmaZ = 1. / TMath::Min(3.3, 1.12 + 0.35 * TMath::Exp(-0.032 * e * e) + 0.75 / TMath::Power(e + 0.24, 3)); // inverse sigma Z

          for (int indx : regions) {
            for (auto p : cpvMatchPoints[indx]) {
              float d = pow((p.first - posX) * sigmaX, 2) + pow((p.second - posZ) * sigmaZ, 2);
              if (d < trackdist) {
                trackdist = d;
              }
            }
          }
        }
        if (trackdist != 99.) {        // was evaluated
          trackdist = sqrt(trackdist); // was squared
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
            trackdist, trackindex,
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
            trackdist, trackindex,
            clu.firedTrigger(),
            clu.getDistanceToBadChannel());
        }
      }
    }
  }

  int CpvMatchIndex(int16_t module, float x, float z)
  {
    // calculatex cell index in CPV grid
    const float cellSizeX = 2 * cpvMaxX / nCpvX;
    const float cellSizeZ = 2 * cpvMaxZ / nCpvZ;
    return (module - 2) * nCpvX * nCpvZ +
           static_cast<int>((x + cpvMaxX) / cellSizeX) * nCpvZ +
           static_cast<int>((z + cpvMaxZ) / cellSizeZ);
  }

  float Nonlinearity(float en)
  {
    // Correct for non-linearity
    // Parameters to be read from ccdb
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
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<caloClusterProducerTask>(cfgc)};
}
