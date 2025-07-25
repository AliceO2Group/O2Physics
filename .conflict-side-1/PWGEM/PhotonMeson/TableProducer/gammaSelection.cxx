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

/// \brief skim cluster information to write photon cluster table in AO2D.root
/// dependencies: skimmergammacalo, skimmergammaconversions, skimmer-phos
/// \author marvin.hemmer@cern.ch

// TODO: add PCM table
#include "PWGEM/PhotonMeson/Utils/gammaSelectionCuts.h"

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include <typeinfo>

// includes for the R recalculation
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/emcalHistoDefinitions.h"

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/GeometryManager.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct gammaSelection {

  uint64_t EMC_CutModeBit;

  Preslice<o2::aod::SkimEMCMTs> perEMCClusterMT = o2::aod::caloextra::clusterId;

  Produces<aod::SkimGammas> tableGammaReco;
  Produces<aod::SkimEMCCuts> tableEMCCuts;

  // Configurable for filter/cuts
  Configurable<float> EMC_minTime{"EMC_minTime", -20., "Minimum cluster time for EMCal time cut"};
  Configurable<float> EMC_maxTime{"EMC_maxTime", +25., "Maximum cluster time for EMCal time cut"};
  Configurable<float> EMC_minM02{"EMC_minM02", 0.1, "Minimum M02 for EMCal M02 cut"};
  Configurable<float> EMC_maxM02{"EMC_maxM02", 0.7, "Maximum M02 for EMCal M02 cut"};
  Configurable<float> EMC_minE{"EMC_minE", 0.7, "Minimum cluster energy for EMCal energy cut"};
  Configurable<int> EMC_minNCell{"EMC_minNCell", 1, "Minimum number of cells per cluster for EMCal NCell cut"};
  Configurable<std::vector<float>> EMC_TM_Eta{"EMC_TM_Eta", {0.01f, 4.07f, -2.5f}, "|eta| <= [0]+(pT+[1])^[2] for EMCal track matching"};
  Configurable<std::vector<float>> EMC_TM_Phi{"EMC_TM_Phi", {0.015f, 3.65f, -2.f}, "|phi| <= [0]+(pT+[1])^[2] for EMCal track matching"};
  Configurable<float> EMC_Eoverp{"EMC_Eoverp", 1.75, "Minimum cluster energy over track momentum for EMCal track matching"};
  Configurable<std::string> EMC_CutMode{"EMC_CutMode", "0", "Cut Mode that is run. Each bit is a different setting. The first bit will use the cuts from the configurables."};

  Configurable<float> PHOS_minTime{"PHOS_minTime", -30., "Minimum cluster time for PHOS time cut"};
  Configurable<float> PHOS_maxTime{"PHOS_maxTime", +30., "Maximum cluster time for PHOS time cut"};
  Configurable<float> PHOS_minM02{"PHOS_minM02", 0.1, "Minimum M02 for PHOS M02 cut"};
  Configurable<float> PHOS_minE{"PHOS_minE", 0.3, "Minimum cluster energy for PHOS energy cut"};
  Configurable<float> PHOS_minENCell{"PHOS_minENCell", 0.1, "Threshold cluster energy for switch for PHOS NCell and M02 cut"};
  Configurable<int> PHOS_minNCell{"PHOS_minNCell", 1, "Minimum number of cells per cluster for PHOS NCell cut"};
  Configurable<float> PHOS_TM_Eta{"PHOS_TM_Eta", 0.02f, "|eta| <= value for PHOS track matching"};
  Configurable<float> PHOS_TM_Phi{"PHOS_TM_Phi", 0.08f, "|phi| <= value for PHOS track matching"};

  Configurable<int> PHOS_QA{"PHOS_QA", 0b0, "Flag to enable PHOS related QA plots. 1st bit for TM QA."};

  HistogramRegistry EMCHistos{
    "EMCHistos",
    {},
    OutputObjHandlingPolicy::QAObject,
    true,
    true};

  HistogramRegistry PHOSHistos{
    "PHOSHistos",
    {},
    OutputObjHandlingPolicy::QAObject,
    true,
    true};

  void init(o2::framework::InitContext&)
  {
    EMC_CutModeBit = stoi(EMC_CutMode, 0, 2);
    std::bitset<64> EMC_CutModeBitSet(EMC_CutModeBit);
    // EMCal
    EMCHistos.add("hClusterEIn", "hClusterEIn", gHistoSpec_clusterECuts);
    EMCHistos.add("hClusterEOut", "hClusterEOut", gHistoSpec_clusterECuts);
    auto hCaloCuts_EMC = EMCHistos.add<TH2>("hCaloCuts_EMC", "hCaloCuts_EMC", kTH2I, {{7, -0.5, 6.5}, {64, -0.5, 63.5}});
    hCaloCuts_EMC->GetXaxis()->SetBinLabel(1, "in");
    hCaloCuts_EMC->GetXaxis()->SetBinLabel(2, "#it{t}_{cluster} cut");
    hCaloCuts_EMC->GetXaxis()->SetBinLabel(3, "#it{M}_{02} cut");
    hCaloCuts_EMC->GetXaxis()->SetBinLabel(4, "#it{E} cut");
    hCaloCuts_EMC->GetXaxis()->SetBinLabel(5, "#it{N}_{cell} cut");
    hCaloCuts_EMC->GetXaxis()->SetBinLabel(6, "TM");
    hCaloCuts_EMC->GetXaxis()->SetBinLabel(7, "out");

    LOG(info) << "| ECMal cluster cut settings:";
    LOG(info) << "|\t Timing cut: " << EMC_minTime << " < t < " << EMC_maxTime;
    LOG(info) << "|\t M02 cut: " << EMC_minM02 << " < M02 < " << EMC_maxM02;
    LOG(info) << "|\t E_min cut: E_cluster > " << EMC_minE;
    LOG(info) << "|\t N_cell cut: N_cell > " << EMC_minNCell;
    LOG(info) << "|\t TM |eta|: |eta| <= " << EMC_TM_Eta->at(0) << " + (pT + " << EMC_TM_Eta->at(1) << ")^" << EMC_TM_Eta->at(2);
    LOG(info) << "|\t TM |phi|: |phi| <= " << EMC_TM_Phi->at(0) << " + (pT + " << EMC_TM_Phi->at(1) << ")^" << EMC_TM_Phi->at(2);
    LOG(info) << "|\t TM E/p: E/p < " << EMC_Eoverp;
    LOG(info) << "|\t Cut bit is set to: " << EMC_CutModeBitSet << std::endl;

    gatherCutsEMC(EMC_minTime, EMC_maxTime, EMC_minM02, EMC_maxM02, EMC_minE, EMC_minNCell, EMC_TM_Eta, EMC_TM_Phi, EMC_Eoverp);

    // PHOS
    PHOSHistos.add("hClusterEIn", "hClusterEIn", gHistoSpec_clusterECuts);
    PHOSHistos.add("hClusterEOut", "hClusterEOut", gHistoSpec_clusterECuts);
    auto hCaloCuts_PHOS = PHOSHistos.add<TH1>("hCaloCuts_PHOS", "hCaloCuts_PHOS", kTH1I, {{7, -0.5, 6.5}});
    hCaloCuts_PHOS->GetXaxis()->SetBinLabel(1, "in");
    hCaloCuts_PHOS->GetXaxis()->SetBinLabel(2, "#it{t}_{cluster} cut");
    hCaloCuts_PHOS->GetXaxis()->SetBinLabel(3, "#it{M}_{02} cut");
    hCaloCuts_PHOS->GetXaxis()->SetBinLabel(4, "#it{E} cut");
    hCaloCuts_PHOS->GetXaxis()->SetBinLabel(5, "#it{N}_{cell} cut");
    hCaloCuts_PHOS->GetXaxis()->SetBinLabel(6, "TM");
    hCaloCuts_PHOS->GetXaxis()->SetBinLabel(7, "out");
    if (PHOS_QA & 0b1) {
      PHOSHistos.add("clusterTM_dEtadPhi", "cluster trackmatching dEta/dPhi;d#it{#eta};d#it{#varphi} (rad)", kTH2F, {{100, -0.2, 0.2}, {100, -0.2, 0.2}}); // dEta dPhi map of matched tracks
    }

    LOG(info) << "| PHOS cluster cut settings:";
    LOG(info) << "|\t Timing cut: " << PHOS_minTime << " < t < " << PHOS_maxTime;
    LOG(info) << "|\t NCell cut: " << PHOS_minNCell << " <= NCell for E >= " << PHOS_minENCell;
    LOG(info) << "|\t M02 cut: " << PHOS_minM02 << " < M02 for E >= " << PHOS_minENCell;
    LOG(info) << "|\t E_min cut: E_cluster > " << PHOS_minE;
    LOG(info) << "|\t TM |eta|: |eta| <= " << PHOS_TM_Eta;
    LOG(info) << "|\t TM |phi|: |phi| <= " << PHOS_TM_Phi << std::endl;
  }

  void processRec(aod::EMEvents const&, aod::SkimEMCClusters const& emcclusters, aod::SkimEMCMTs const& matchedtracks, aod::PHOSClusters const& phosclusters)
  {
    for (const auto& emccluster : emcclusters) { // loop of EMC clusters
      uint64_t EMC_CutBit = doPhotonCutsEMC(EMC_CutModeBit, emccluster, matchedtracks, perEMCClusterMT, EMCHistos);
      tableEMCCuts(emccluster.globalIndex(), EMC_CutBit);
    } // end loop of EMC clusters

    for (const auto& phoscluster : phosclusters) { // loop over PHOS clusters
      PHOSHistos.fill(HIST("hClusterEIn"), phoscluster.e(), 0);
      PHOSHistos.fill(HIST("hCaloCuts_PHOS"), 0);

      if (phoscluster.time() > PHOS_maxTime || phoscluster.time() < PHOS_minTime) {
        PHOSHistos.fill(HIST("hCaloCuts_PHOS"), 1);
        continue;
      }
      if (!(phoscluster.e() >= PHOS_minENCell && phoscluster.m02() > PHOS_minM02)) {
        PHOSHistos.fill(HIST("hCaloCuts_PHOS"), 2);
        continue;
      }
      if (phoscluster.e() <= PHOS_minE) {
        PHOSHistos.fill(HIST("hCaloCuts_PHOS"), 3);
        continue;
      }
      if (!(phoscluster.e() >= PHOS_minENCell && phoscluster.nCells() > PHOS_minNCell)) {
        PHOSHistos.fill(HIST("hCaloCuts_PHOS"), 4);
        continue;
      }

      // TODO: add track matching for PHOS when available!
      // track matching
      bool hasMatchedTrack_PHOS = false;
      // double dEta_PHOS, dPhi_PHOS;
      // // only consider closest match
      // dEta_PHOS = phoscluster.tracketa() - phoscluster.eta();
      // dPhi_PHOS = phoscluster.trackphi() - phoscluster.phi();
      // if ((fabs(dEta_PHOS) < PHOS_TM_Eta) && (fabs(dPhi_PHOS) < PHOS_TM_Phi)) {
      //   hasMatchedTrack_PHOS = true;
      //   if (PHOS_QA & 0b1) {
      //     EMCHistos.fill(HIST("clusterTM_dEtadPhi"), dEta_PHOS, dPhi_PHOS);
      //   }
      // }
      if (hasMatchedTrack_PHOS) {
        PHOSHistos.fill(HIST("hCaloCuts_PHOS"), 5);
      } else {
        PHOSHistos.fill(HIST("hClusterEOut"), phoscluster.e(), 0);
        PHOSHistos.fill(HIST("hCaloCuts_PHOS"), 6);
        tableGammaReco(phoscluster.collisionId(), 2,
                       phoscluster.e(), phoscluster.eta(), phoscluster.phi(), 0, phoscluster.globalIndex());
      }
    } // end loop of PHOS clusters
  }
  PROCESS_SWITCH(gammaSelection, processRec, "process only reconstructed info", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<gammaSelection>(cfgc)};
  return workflow;
}
