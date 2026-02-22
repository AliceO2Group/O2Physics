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

/// \file nonLinProducer.cxx
/// \brief produces nonLin tables for PCM and EMCal
/// \author marvin.hemmer@cern.ch
/// dependencies: skimmer-gamma-calo

#include "PWGEM/PhotonMeson/Core/EMNonLin.h"
#include "PWGEM/PhotonMeson/DataModel/GammaTablesRedux.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/emcalHistoDefinitions.h"

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <cstdint>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::emcdownscaling;
using namespace o2::pwgem::nonlin;

struct NonLinProducer {

  enum CentralityEstimator : uint8_t {
    None = 0,
    CFT0A = 1,
    CFT0C,
    CFT0M,
    NCentralityEstimators
  };

  Produces<aod::NonLinV0s> tableNonLinV0s;
  Produces<aod::NonLinEmcClusters> tableNonLinClusters;

  Configurable<int> centEstimator{"centEstimator", 2, "Centrality estimation (FT0A: 1, FT0C: 2, FT0M: 3)"};

  HistogramRegistry historeg{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  using EMCalPhotons = soa::Join<aod::EMCEMEventIds, aod::MinClusters>;
  using PcmPhotons = soa::Join<aod::V0PhotonsKF, aod::V0KFEMEventIds>;

  using Colls = soa::Join<aod::EMEvents_004, aod::EMEventsCent_000>;

  EMNonLin emNonLinEMC;
  EMNonLin emNonLinPCM;

  void init(o2::framework::InitContext&)
  {
    historeg.add("QA/EMC/EIn", "Energy of clusters before cuts", gHistoSpecClusterE);
    historeg.add("QA/EMC/EOut", "Energy of clusters after cuts", gHistoSpecClusterE);

    historeg.add("QA/EMC/PtIn", "Energy of clusters before cuts", gHistoSpecPt);
    historeg.add("QA/EMC/PtOut", "Energy of clusters after cuts", gHistoSpecPt);

    historeg.add("QA/PCM/PtIn", "Energy of clusters before cuts", gHistoSpecPt);
    historeg.add("QA/PCM/PtOut", "Energy of clusters after cuts", gHistoSpecPt);
  }

  /// Get the centrality
  /// \param collision is the collision with the centrality information
  template <o2::soa::is_iterator TCollision>
  float getCentrality(TCollision const& collision)
  {
    float cent = -999.;
    switch (centEstimator) {
      case CentralityEstimator::CFT0M:
        cent = collision.centFT0M();
        break;
      case CentralityEstimator::CFT0A:
        cent = collision.centFT0A();
        break;
      case CentralityEstimator::CFT0C:
        cent = collision.centFT0C();
        break;
      default:
        LOG(warning) << "Centrality estimator not valid. Possible values are T0M, T0A, T0C. Fallback to T0C";
        cent = collision.centFT0C();
        break;
    }
    return cent;
  }

  template <o2::soa::is_table TClusters, o2::soa::is_iterator TCollisio>
  void runEMC(TClusters const& clusters, TCollisio& collision)
  {
    float nonLinE = 0.f;
    float nonLinPt = 0.f;

    float nonLinFactor = 1.f;

    int32_t collIndex = collision.globalIndex();
    for (const auto& cluster : clusters) {

      // check that we are at the correct collision
      if (cluster.emphotoneventId() != collIndex) {
        collIndex = cluster.emphotoneventId();
        collision.setCursor(collIndex);
      }

      // fill before non lin histograms
      historeg.fill(HIST("QA/EMC/EIn"), cluster.e());
      historeg.fill(HIST("QA/EMC/PtIn"), cluster.pt());

      // get NonLin factor from class dependent on the centrality
      nonLinFactor = emNonLinEMC.getCorrectionFactor(cluster.e(), o2::pwgem::nonlin::EMNonLin::PhotonType::kEMC, getCentrality(collision));

      nonLinE = nonLinFactor * cluster.e();
      nonLinPt = nonLinFactor * cluster.pt();

      // fill after non lin histograms
      historeg.fill(HIST("QA/EMC/EIn"), nonLinE);
      historeg.fill(HIST("QA/EMC/PtIn"), nonLinPt);

      tableNonLinClusters(nonLinE, nonLinPt);
    }
  }

  template <o2::soa::is_table TV0, o2::soa::is_iterator TCollisio>
  void runPCM(TV0 const& v0s, TCollisio& collision)
  {
    float nonLinPt = 0.f;

    float nonLinFactor = 1.f;

    int32_t collIndex = collision.globalIndex();
    for (const auto& v0 : v0s) {

      // check that we are at the correct collision
      if (v0.emphotoneventId() != collIndex) {
        collIndex = v0.emphotoneventId();
        collision.setCursor(collIndex);
      }

      // fill before non lin histograms
      historeg.fill(HIST("QA/PCM/PtIn"), v0.pt());

      // get NonLin factor from class dependent on the centrality
      nonLinFactor = emNonLinEMC.getCorrectionFactor(v0.pt(), o2::pwgem::nonlin::EMNonLin::PhotonType::kPCM, getCentrality(collision));

      nonLinPt = nonLinFactor * v0.pt();

      // fill after non lin histograms
      historeg.fill(HIST("QA/PCM/PtIn"), nonLinPt);

      tableNonLinV0s(nonLinPt);
    }
  }

  void processEMC(Colls const& collisions, EMCalPhotons const& emcclusters)
  {

    if (emcclusters.size() == 0) {
      return;
    }

    auto collision = collisions.begin();
    runEMC(emcclusters, collision);
  }
  PROCESS_SWITCH(NonLinProducer, processEMC, "Create Non Lin table for EMC.", false);

  void processPCM(Colls const& collisions, PcmPhotons const& pcmPhotons)
  {
    if (pcmPhotons.size() == 0) {
      return;
    }
    auto collision = collisions.begin();
    runPCM(pcmPhotons, collision);
  }
  PROCESS_SWITCH(NonLinProducer, processPCM, "Create Non Lin table for PCM.", false);

  void processEMCDummy(EMCalPhotons::iterator const& emccluster)
  {
    tableNonLinClusters(emccluster.e(), emccluster.pt());
  }
  PROCESS_SWITCH(NonLinProducer, processEMCDummy, "Create dummy Non Lin table for EMC.", true);

  void processPCMDummy(PcmPhotons::iterator const& pcmPhoton)
  {
    tableNonLinV0s(pcmPhoton.pt());
  }
  PROCESS_SWITCH(NonLinProducer, processPCMDummy, "Createdumy  Non Lin table for PCM.", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<NonLinProducer>(cfgc)};
  return workflow;
}
