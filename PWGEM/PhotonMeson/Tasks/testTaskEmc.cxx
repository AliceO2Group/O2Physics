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

/// \file testTaskEmc.cxx
/// \brief Task to test values for EMCal
/// \author M. Hemmer, marvin.hemmer@cern.ch

#include "PWGEM/PhotonMeson/DataModel/GammaTablesRedux.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/emcalHistoDefinitions.h"
//

#include <CommonConstants/MathConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <string>
#include <vector>

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

struct TestTaskEmc {

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  void init(InitContext&)
  {
    const AxisSpec energyAxis{10000, 0., 10., "#it{E}_{clus} (GeV)"};
    const AxisSpec energyAxis2{10000, 0.699, 0.701, "#it{E}_{clus} (GeV)"};
    const AxisSpec m02axis{200, 0., +2, "#it{M}_{02} (a.u.)"};
    const AxisSpec phiAxis{200, 0., o2::constants::math::TwoPI, "#it{#varphi} (rad)"};
    const AxisSpec etaAxis{140, -0.7, +0.7, "#it{#eta}"};
    const AxisSpec nCellsAxis{140, -0.7, +0.7, "#it{N}_{cell}"};
    const AxisSpec isExoticAxis{2, -0.5, +1.5, "isExotic"};
    const AxisSpec timeAxis{2000, -1000., +1000., "#it{N}_{cell}"};

    const AxisSpec deltaEnergyAxis{1000, -0.01, 0.01, "#Delta#it{E}_{clus} (GeV)"};
    const AxisSpec deltaM02axis{1000, -0.01, 0.01, "#Delta#it{M}_{02} (a.u.)"};
    const AxisSpec deltaNCellsAxis{1000, -0.01, 0.01, "#Delta#it{N}_{cell}"};
    const AxisSpec deltaTimeAxis{1000, -0.01, 0.01, "#Delta#it{N}_{cell}"};

    const AxisSpec deltaPhiAxis{400, -0.05, +0.05, "#Delta#it{#varphi} (rad)"};
    const AxisSpec deltaEtaAxis{400, -0.05, +0.05, "#Delta#it{#eta}"};

    registry.add("Old/hEnergy", "hEnergy", HistType::kTH1D, {energyAxis});
    registry.add("Old/hEnergy2", "hEnergy2", HistType::kTH1D, {energyAxis2});
    registry.add("Old/hM02", "hM02", HistType::kTH1D, {m02axis});
    registry.add("Old/hEta", "hEta", HistType::kTH1D, {etaAxis});
    registry.add("Old/hPhi", "hPhi", HistType::kTH1D, {phiAxis});
    registry.add("Old/hNCells", "hNCells", HistType::kTH1D, {nCellsAxis});
    registry.add("Old/hExotic", "hExotic", HistType::kTH1D, {isExoticAxis});
    registry.add("Old/hTime", "hTime", HistType::kTH1D, {timeAxis});

    registry.add("Old/Prim/hDeltaEta", "hDeltaEta", HistType::kTH1D, {deltaEtaAxis});
    registry.add("Old/Prim/hDeltaPhi", "hDeltaPhi", HistType::kTH1D, {deltaPhiAxis});
    registry.add("Old/Prim/hP", "hP", HistType::kTH1D, {energyAxis});
    registry.add("Old/Prim/hPt", "hPt", HistType::kTH1D, {energyAxis});

    registry.add("Old/Sec/hDeltaEta", "hDeltaEta", HistType::kTH1D, {deltaEtaAxis});
    registry.add("Old/Sec/hDeltaPhi", "hDeltaPhi", HistType::kTH1D, {deltaPhiAxis});
    registry.add("Old/Sec/hP", "hP", HistType::kTH1D, {energyAxis});
    registry.add("Old/Sec/hPt", "hPt", HistType::kTH1D, {energyAxis});

    registry.addClone("Old/", "Redux/");

    registry.add("hDeltaEnergy", "hDeltaEnergy", HistType::kTH1D, {deltaEnergyAxis});
    registry.add("hDeltaM02", "hDeltaM02", HistType::kTH1D, {deltaM02axis});
    registry.add("hDeltaEta", "hDeltaEta", HistType::kTH1D, {deltaEtaAxis});
    registry.add("hDeltaPhi", "hDeltaPhi", HistType::kTH1D, {deltaPhiAxis});
    registry.add("hDeltaNCells", "hDeltaNCells", HistType::kTH1D, {deltaNCellsAxis});
    registry.add("hDeltaTime", "hDeltaTime", HistType::kTH1D, {deltaTimeAxis});

  }; // end init

  // EMCal calibration same event
  void processEMCalCalib(aod::SkimEMCClusters_001 const& oldClusters, aod::MinClusters const& reduxClusters, aod::MinMTracks const& reduxMatchedTracks, aod::MinMSTracks const& reduxMatchedSecondaries)
  {
    auto oldCluster = oldClusters.begin();
    auto reduxCluster = reduxClusters.begin();

    while (oldCluster != oldClusters.end() && reduxCluster != reduxClusters.end()) {
      registry.fill(HIST("hDeltaEnergy"), oldCluster.e() - reduxCluster.e());
      registry.fill(HIST("hDeltaM02"), oldCluster.m02() - reduxCluster.m02());
      registry.fill(HIST("hDeltaEta"), oldCluster.eta() - reduxCluster.eta());
      registry.fill(HIST("hDeltaPhi"), oldCluster.phi() - reduxCluster.phi());
      registry.fill(HIST("hDeltaNCells"), oldCluster.nCells() - reduxCluster.nCells());
      registry.fill(HIST("hDeltaTime"), oldCluster.time() - reduxCluster.time());
      ++oldCluster;
      ++reduxCluster;
    }
    for (const auto& cluster : oldClusters) {
      registry.fill(HIST("Old/hEnergy"), cluster.e());
      registry.fill(HIST("Old/hEnergy2"), cluster.e());
      registry.fill(HIST("Old/hM02"), cluster.m02());
      registry.fill(HIST("Old/hEta"), cluster.eta());
      registry.fill(HIST("Old/hPhi"), cluster.phi());
      registry.fill(HIST("Old/hNCells"), cluster.nCells());
      registry.fill(HIST("Old/hExotic"), cluster.isExotic());
      registry.fill(HIST("Old/hTime"), cluster.time());

      for (const auto& var : cluster.deltaEta()) {
        registry.fill(HIST("Old/Prim/hDeltaEta"), var);
      }
      for (const auto& var : cluster.deltaPhi()) {
        registry.fill(HIST("Old/Prim/hDeltaPhi"), var);
      }
      for (const auto& dPhi : cluster.trackp()) {
        registry.fill(HIST("Old/Prim/hP"), dPhi);
      }
      for (const auto& var : cluster.trackpt()) {
        registry.fill(HIST("Old/Prim/hPt"), var);
      }
      for (const auto& var : cluster.deltaEtaSec()) {
        registry.fill(HIST("Old/Sec/hDeltaEta"), var);
      }
      for (const auto& var : cluster.deltaPhiSec()) {
        registry.fill(HIST("Old/Sec/hDeltaPhi"), var);
      }
      for (const auto& var : cluster.trackpSec()) {
        registry.fill(HIST("Old/Sec/hP"), var);
      }
      for (const auto& var : cluster.trackptSec()) {
        registry.fill(HIST("Old/Sec/hPt"), var);
      }
    }

    for (const auto& cluster : reduxClusters) {
      registry.fill(HIST("Redux/hEnergy"), cluster.e());
      registry.fill(HIST("Redux/hEnergy2"), cluster.e());
      registry.fill(HIST("Redux/hM02"), cluster.m02());
      registry.fill(HIST("Redux/hEta"), cluster.eta());
      registry.fill(HIST("Redux/hPhi"), cluster.phi());
      registry.fill(HIST("Redux/hNCells"), cluster.nCells());
      registry.fill(HIST("Redux/hExotic"), cluster.isExotic());
      registry.fill(HIST("Redux/hTime"), cluster.time());
    }

    for (const auto& prim : reduxMatchedTracks) {
      registry.fill(HIST("Redux/Prim/hDeltaEta"), prim.deltaEta());
      registry.fill(HIST("Redux/Prim/hDeltaPhi"), prim.deltaPhi());
      registry.fill(HIST("Redux/Prim/hP"), prim.trackP());
      registry.fill(HIST("Redux/Prim/hPt"), prim.trackPt());
    }

    for (const auto& sec : reduxMatchedSecondaries) {
      registry.fill(HIST("Redux/Sec/hDeltaEta"), sec.deltaEta());
      registry.fill(HIST("Redux/Sec/hDeltaPhi"), sec.deltaPhi());
      registry.fill(HIST("Redux/Sec/hP"), sec.trackP());
      registry.fill(HIST("Redux/Sec/hPt"), sec.trackPt());
    }
  }

  PROCESS_SWITCH(TestTaskEmc, processEMCalCalib, "Process EMCal calibration same event", true);
}; // End struct TestTaskEmc

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<TestTaskEmc>(cfgc)};
}
