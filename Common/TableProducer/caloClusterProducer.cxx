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
  Produces<aod::CaloClusters> clusters;
  Configurable<std::string> calorimeter{"caloType", "BOTH", "PHOS, EMCAL, BOTH"};
  Configurable<bool> isMC{"isMC", 0, "0 - data, 1 - MC"};
  Configurable<bool> useCoreE{"coreE", 0, "0 - full energy, 1 - core energy"};

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  std::unique_ptr<o2::phos::Geometry> geomPHOS;
  std::unique_ptr<o2::phos::Clusterer> clusterizerPHOS;
  std::vector<o2::phos::Cell> phosCells;
  std::vector<o2::phos::TriggerRecord> phosCellTRs;
  std::vector<o2::phos::CluElement> outputCluElements;
  std::vector<o2::phos::Cluster> outputPHOSClusters;
  std::vector<o2::phos::TriggerRecord> outputPHOSClusterTrigRecs;

  void init(o2::framework::InitContext&)
  {
    ccdb->setURL(o2::base::NameConf::getCCDBServer());
    if (calorimeter->compare("PHOS") == 0 || calorimeter->compare("BOTH") == 0) {
      geomPHOS = std::make_unique<o2::phos::Geometry>("PHOS");
      clusterizerPHOS = std::make_unique<o2::phos::Clusterer>();
    }
  }

  void process(o2::aod::BCs const& bcs,
               o2::aod::Collisions const& colls,
               o2::aod::Calos const& cells,
               o2::aod::CaloTriggers const& ctrs)
  {
    // Make map between collision and BC tables
    //  map: (bcId_long,collision index)
    std::map<int64_t, int> bcMap;
    int collId = 0;
    // TODO! handle several collisions assigned to same BC
    for (auto cl : colls) {
      bcMap[cl.bc().globalBC()] = collId;
      collId++;
    }

    if (calorimeter->compare("PHOS") == 0 || calorimeter->compare("BOTH") == 0) {
      const int kPHOS = 0;
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
        // find collision corresponding to current BC
        auto clvtx = colls.begin() + (bcMap[cluTR.getBCData().toLong()]);

        // Extract primary vertex
        TVector3 vtx = {clvtx.posX(), clvtx.posY(), clvtx.posZ()};

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

          // Track/CPV match will be done in independent task
          float trackdist = 999.;
          int trackindex = -1;

          float lambdaShort = 0., lambdaLong = 0.;
          clu.getElipsAxis(lambdaShort, lambdaLong);

          clusters(bcMap[cluTR.getBCData().toLong()], kPHOS, mom.X(), mom.Y(), mom.Z(), e,
                   mod, clu.getMultiplicity(), globaPos.X(), globaPos.Y(), globaPos.Z(),
                   clu.getTime(), clu.getNExMax(), lambdaShort, lambdaLong, trackdist, trackindex,
                   clu.firedTrigger(), clu.getDistanceToBadChannel());
        }
      }
    } // end isPHOS
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
