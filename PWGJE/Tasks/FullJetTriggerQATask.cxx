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

// jet trigger QA task
//
// Author: Gijs van Weelden
//
#include <bitset>

#include "TH1F.h"
#include "TTree.h"

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/RunningWorkflowInfo.h"

#include "Common/DataModel/EventSelection.h"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/Core/JetFinder.h"

#include "EventFiltering/filterTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using selectedClusters = o2::soa::Filtered<o2::aod::EMCALClusters>;

struct JetTriggerQA {

  enum TriggerType_t {
    kMinBias,
    kEmcal,
    kEmcalJetFull,
    kEmcalJetNeutral,
    kEmcalGammaHigh,
    kDcalGammaHigh,
    kEmcalGammaLow,
    kDcalGammaLow,
    kNTriggers
  };

  // Preslice<aod::JetTrackConstituents> perJetTrackConstituents = o2::aod::jetconstituents::jetId;
  // Preslice<aod::JetClusterConstituents> perJetClusterConstituents = o2::aod::jetconstituents::jetId;

  HistogramRegistry registry{"registry"};

  Configurable<float> f_jetPtMin{"f_jetPtMin", 0.0, "minimum jet pT cut"};
  Configurable<float> f_SD_zCut{"f_SD_zCut", 0.1, "soft drop z cut"};
  Configurable<float> f_SD_beta{"f_SD_beta", 0.0, "soft drop beta"};
  Configurable<float> f_jetR{"f_jetR", 0.4, "jet resolution parameter that you have triggered on"};
  Configurable<float> f_PhiEmcalOrDcal{"f_PhiEmcalOrDcal", 4, "if cluster phi is less than this value, count it to be EMCAL"};
  Configurable<float> f_minClusterTime{"f_minClusterTime", -999, "Min. cluster time for gamma trigger (ns)"};
  Configurable<float> f_maxClusterTime{"f_maxClusterTime", 999, "Max. cluster time for gamma trigger (ns)"};
  Configurable<std::vector<float>> f_ang_kappa{"f_ang_kappa", {1.0, 1.0, 2.0}, "angularity momentum exponent"};
  Configurable<std::vector<float>> f_ang_alpha{"f_ang_alpha", {1.0, 2.0, 1.0}, "angularity angle exponent"};
  Configurable<bool> b_JetsInEmcalOnly{"b_JetsInEmcalOnly", true, "fill histograms only for jets inside the EMCAL"};
  Configurable<bool> b_IgnoreEmcalFlag{"b_IgnoreEmcalFlag", false, "ignore the EMCAL live flag check"};
  Configurable<bool> b_DoFiducialCut{"b_DoFiducialCut", false, "do a fiducial cut on jets to check if they are in the emcal"};
  Configurable<bool> b_RejectExoticClusters{"b_RejectExoticClusters", true, "Reject exotic clusters"};
  Configurable<int> f_GammaObservable{"f_gammaObservable", 0, "Observable for gamma trigger (0 - energy, 1 - pt)"};

  Configurable<float> cfgVertexCut{"cfgVertexCut", 10.0f, "Accepted z-vertex range"};
  Configurable<std::string> mClusterDefinition{"clusterDefinition", "kV3Default", "cluster definition to be selected, e.g. V3Default"};

  std::vector<fastjet::PseudoJet> jetConstituents;
  std::vector<fastjet::PseudoJet> jetClusterConstituents;
  std::vector<fastjet::PseudoJet> jetReclustered;
  JetFinder jetReclusterer;

  void init(InitContext& initContext)
  {
    jetReclusterer.isReclustering = true;
    jetReclusterer.algorithm = fastjet::JetAlgorithm::cambridge_algorithm;
    jetReclusterer.jetR = f_jetR * 2.5; // Use larger R for CA reclustering to prevent losing particles

    Int_t nPtBins = 200;
    Float_t kMinPt = 0.;
    Float_t kMaxPt = 201.;
    Int_t nPhiBins = 18 * 8;
    Float_t kMinPhi = 0.;
    Float_t kMaxPhi = 2. * TMath::Pi();
    Int_t nEtaBins = 100;
    Float_t kMinEta = -1.;
    Float_t kMaxEta = 1.;
    Int_t nRBins = 6;
    Float_t kMinR = 0.05;
    Float_t kMaxR = 0.65;

    // EMCAL gamma observables (for histogram and axis title)
    std::string observableName = (f_GammaObservable == 0) ? "Energy" : "#it{p}_{T}",
                observableAxisTitle = (f_GammaObservable == 0) ? "E (GeV)" : "#it{p}_{T} (GeV/c)",
                observableAxisTitleMax = (f_GammaObservable == 0) ? "E^{max} (GeV)" : "#it{p}_{T}^{max} (GeV/c)";

    AxisSpec rAxis = {nRBins, kMinR, kMaxR, "#it{R}"};
    AxisSpec jetPtAxis = {nPtBins, kMinPt, kMaxPt, "#it{p}_{T}^{jet} (GeV/c)"};
    AxisSpec jetMaxPtAxis = {nPtBins, kMinPt, kMaxPt, "#it{p}_{T,max}^{jet} (GeV/c)"};
    AxisSpec etaAxis = {nEtaBins, kMinEta, kMaxEta, "#eta"};
    AxisSpec phiAxis = {nPhiBins, kMinPhi, kMaxPhi, "#phi"};
    AxisSpec observableAxisCluster{nPtBins, kMinPt, kMaxPt / 2, observableAxisTitle},
      observableAxisMaxCluster{nPtBins, kMinPt, kMaxPt / 2, observableAxisTitleMax},
      ptAxisTrackInJet(nPtBins, kMinPt, kMaxPt / 2, "#it{p}_{T}^{track} (GeV/c)"),
      ptAxisClusterInJet(nPtBins, kMinPt, kMaxPt / 2, "#it{p}_{T}^{clus} (GeV/c)");

    HistogramConfigSpec hJetRPtEtaPhi({HistType::kTHnF, {rAxis, jetPtAxis, etaAxis, phiAxis}});
    HistogramConfigSpec hJetRMaxPtEtaPhi({HistType::kTHnF, {rAxis, jetPtAxis, etaAxis, phiAxis}});
    HistogramConfigSpec hJetRPtEtaPhiNoFiducial({HistType::kTHnF, {rAxis, jetPtAxis, etaAxis, phiAxis}});
    HistogramConfigSpec hJetRMaxPtEtaPhiNoFiducial({HistType::kTHnF, {rAxis, jetPtAxis, etaAxis, phiAxis}});

    registry.add("jetRPtEtaPhi", "JetRPtEtaPhi", hJetRPtEtaPhi);
    registry.add("jetRMaxPtEtaPhi", "JetRMaxPtEtaPhi", hJetRMaxPtEtaPhi);
    registry.add("jetRPtEtaPhiNoFiducial", "JetRPtEtaPhiNoFiducial", hJetRPtEtaPhiNoFiducial);
    registry.add("jetRMaxPtEtaPhiNoFiducial", "JetRMaxPtEtaPhiNoFiducial", hJetRMaxPtEtaPhiNoFiducial);

    registry.add("hProcessedEvents", "Processed events", HistType::kTH1I, {{14, -0.5, 13.5, "Trigger type"}});
    auto histProcessed = registry.get<TH1>(HIST("hProcessedEvents"));
    histProcessed->GetXaxis()->SetBinLabel(1, "MB");
    histProcessed->GetXaxis()->SetBinLabel(2, "EMC");
    histProcessed->GetXaxis()->SetBinLabel(3, "Selected Full Jet");
    histProcessed->GetXaxis()->SetBinLabel(4, "Selected Neutral Jet");
    histProcessed->GetXaxis()->SetBinLabel(5, "Selected Gamma high EMCAL");
    histProcessed->GetXaxis()->SetBinLabel(6, "Selected Gamma high DCAL");
    histProcessed->GetXaxis()->SetBinLabel(7, "Selected Gamma low EMCAL");
    histProcessed->GetXaxis()->SetBinLabel(8, "Selected Gamma low DCAL");
    histProcessed->GetXaxis()->SetBinLabel(9, "Selected full jet and Gamma high EMCAL");
    histProcessed->GetXaxis()->SetBinLabel(10, "Selected full jet and Gamma high DCAL");
    histProcessed->GetXaxis()->SetBinLabel(11, "Selected neutral jet and Gamma high EMCAL");
    histProcessed->GetXaxis()->SetBinLabel(12, "Selected neutral jet and Gamma high DCAL");
    histProcessed->GetXaxis()->SetBinLabel(13, "Selected Gamma high EMCAL and high Gamma DCAL");
    histProcessed->GetXaxis()->SetBinLabel(14, "Selected full jet and neutral jet");

    std::array<std::string, TriggerType_t::kNTriggers> triggerlabels = {{"MB", "EMC", "EMC jet full", "EMC jet neutral", "EMC gamma high", "DCL gamma high", "EMC gamma low", "DCL gamma low"}};
    registry.add("hTriggerCorrelation", "Correlation between EMCAL triggers", HistType::kTH2D, {{TriggerType_t::kNTriggers, -0.5, TriggerType_t::kNTriggers - 0.5, "Main trigger"}, {TriggerType_t::kNTriggers, -0.5, TriggerType_t::kNTriggers - 0.5, "Associated trigger"}});
    auto triggerCorrelation = registry.get<TH2>(HIST("hTriggerCorrelation"));
    for (std::size_t triggertype = 0; triggertype < TriggerType_t::kNTriggers; triggertype++) {
      triggerCorrelation->GetXaxis()->SetBinLabel(triggertype + 1, triggerlabels[triggertype].data());
      triggerCorrelation->GetYaxis()->SetBinLabel(triggertype + 1, triggerlabels[triggertype].data());
    }

    // Histograms for events where the EMCAL is live
    registry.add("hJetRPtEta", "Jets #it{p}_{T} and #eta", HistType::kTH3F, {rAxis, jetPtAxis, etaAxis});
    registry.add("hJetRPtPhi", "Jets #it{p}_{T} and #phi", HistType::kTH3F, {rAxis, jetPtAxis, phiAxis});
    registry.add("hJetRMaxPtEta", "Leading jets #it{p}_{T} and #eta", HistType::kTH3F, {rAxis, jetPtAxis, etaAxis});
    registry.add("hJetRMaxPtPhi", "Leading jets #it{p}_{T} and #phi", HistType::kTH3F, {rAxis, jetPtAxis, phiAxis});
    registry.add("hClusterPtEta", Form("Cluster %s and #eta", observableName.data()), HistType::kTH2F, {observableAxisCluster, etaAxis});
    registry.add("hClusterPtPhi", Form("Cluster %s and #phi", observableName.data()), HistType::kTH2F, {observableAxisCluster, phiAxis});
    registry.add("hClusterPtEtaPhi", Form("Cluster %s, #eta and #phi", observableName.data()), HistType::kTH3F, {observableAxisCluster, etaAxis, phiAxis});
    registry.add("hClusterEMCALMaxPtEtaPhi", Form("Leading cluster %s, #eta and #phi (EMCAL)", observableName.data()), HistType::kTH3F, {observableAxisCluster, etaAxis, phiAxis});
    registry.add("hClusterDCALMaxPtEtaPhi", Form("Leading cluster %s, #eta and #phi (DCAL)", observableName.data()), HistType::kTH3F, {observableAxisCluster, etaAxis, phiAxis});
    registry.add("hClusterMaxPtEta", Form("Leading clusters %s and #eta", observableName.data()), HistType::kTH2F, {observableAxisCluster, etaAxis});
    registry.add("hClusterMaxPtPhi", Form("Leading clusters %s and #phi", observableName.data()), HistType::kTH2F, {observableAxisCluster, phiAxis});
    registry.add("hJetRPtTrackPt", "Jets", HistType::kTH3F, {rAxis, jetPtAxis, ptAxisTrackInJet});
    registry.add("hJetRPtClusterPt", "Jets", HistType::kTH3F, {rAxis, jetPtAxis, ptAxisClusterInJet});
    registry.add("hJetRPtPtd", "Jets", HistType::kTH3F, {rAxis, jetPtAxis, {nPtBins / 2, 0., 1., "p_{t,D}"}});
    registry.add("hJetRPtChargeFrag", "Jets", HistType::kTH3F, {rAxis, jetPtAxis, {nPtBins / 2, 0., 1., "z"}});
    registry.add("hJetRPtNEF", "Jets", HistType::kTH3F, {rAxis, jetPtAxis, {nPtBins / 2, 0., 1., "NEF"}});
    registry.add("hJetRPtZTheta", "Jets", HistType::kTH3F, {rAxis, jetPtAxis, {nPtBins / 2, 0., 1., "z#theta"}});
    registry.add("hJetRPtZSqTheta", "Jets", HistType::kTH3F, {rAxis, jetPtAxis, {nPtBins / 2, 0., 1., "z^{2} #theta"}});
    registry.add("hJetRPtZThetaSq", "Jets", HistType::kTH3F, {rAxis, jetPtAxis, {nPtBins / 2, 0., 1., "z #theta^{2}"}});
    registry.add("hJetRMaxPtClusterMaxPt", "Leading jets and clusters", HistType::kTH3F, {rAxis, jetPtAxis, observableAxisCluster});

    // Histograms for triggered events
    registry.add("hSelectedJetRPtEta", "Selected jets #it{p}_{T} and #eta", HistType::kTH3F, {rAxis, jetPtAxis, etaAxis});
    registry.add("hSelectedJetRPtPhi", "Selected jets #it{p}_{T} and #phi", HistType::kTH3F, {rAxis, jetPtAxis, phiAxis});
    registry.add("hSelectedJetRMaxPtEta", "Leading selected jets #it{p}_{T} and #eta", HistType::kTH3F, {rAxis, jetPtAxis, etaAxis});
    registry.add("hSelectedJetRMaxPtPhi", "Leading selected jets #it{p}_{T} and #phi", HistType::kTH3F, {rAxis, jetPtAxis, phiAxis});
    registry.add("hSelectedClusterPtEta", Form("Selected Cluster %s and #eta", observableName.data()), HistType::kTH2F, {observableAxisCluster, etaAxis});
    registry.add("hSelectedClusterPtPhi", Form("Selected Cluster %s and #phi", observableName.data()), HistType::kTH2F, {observableAxisCluster, phiAxis});
    registry.add("hSelectedClusterPtEtaPhi", Form("Selected cluster %s, #eta and #phi", observableName.data()), HistType::kTH3F, {observableAxisCluster, etaAxis, phiAxis});
    registry.add("hSelectedClusterMaxPtEtaPhi", Form("Leading Selected cluster %s, #eta and #phi", observableName.data()), HistType::kTH3F, {observableAxisCluster, etaAxis, phiAxis});
    registry.add("hSelectedClusterMaxPtEta", Form("Leading selected clusters %s and #eta", observableName.data()), HistType::kTH2F, {observableAxisCluster, etaAxis});
    registry.add("hSelectedClusterMaxPtPhi", Form("Leading selected clusters %s and #phi", observableName.data()), HistType::kTH2F, {observableAxisCluster, phiAxis});
    registry.add("hSelectedJetRPtTrackPt", "Selected jets", HistType::kTH3F, {rAxis, jetPtAxis, ptAxisTrackInJet});
    registry.add("hSelectedJetRPtClusterPt", "Selected jets", HistType::kTH3F, {rAxis, jetPtAxis, ptAxisClusterInJet});
    registry.add("hSelectedJetRPtPtd", "Selected jets", HistType::kTH3F, {rAxis, jetPtAxis, {nPtBins / 2, 0., 1., "p_{t,D}"}});
    registry.add("hSelectedJetRPtChargeFrag", "Selected jets", HistType::kTH3F, {rAxis, jetPtAxis, {nPtBins / 2, 0., 1., "z"}});
    registry.add("hSelectedJetRPtNEF", "Selected jets", HistType::kTH3F, {rAxis, jetPtAxis, {nPtBins / 2, 0., 1., "NEF"}});
    registry.add("hSelectedJetRPtZTheta", "Selected jets", HistType::kTH3F, {rAxis, jetPtAxis, {nPtBins / 2, 0., 1., "z#theta"}});
    registry.add("hSelectedJetRPtZSqTheta", "Selected jets", HistType::kTH3F, {rAxis, jetPtAxis, {nPtBins / 2, 0., 1., "z^{2} #theta"}});
    registry.add("hSelectedJetRPtZThetaSq", "Selected jets", HistType::kTH3F, {rAxis, jetPtAxis, {nPtBins / 2, 0., 1., "z #theta^{2}"}});
    registry.add("hSelectedGammaEMCALPtEta", Form("Selected Gamma %s and #eta (EMCAL, high threshold)", observableName.data()), HistType::kTH2F, {observableAxisCluster, etaAxis});
    registry.add("hSelectedGammaEMCALPtPhi", Form("Selected Gamma %s and #phi (EMCAL, high threshold)", observableName.data()), HistType::kTH2F, {observableAxisCluster, phiAxis});
    registry.add("hSelectedGammaDCALPtEta", Form("Selected Gamma %s and #eta (DCAL, high threshold)", observableName.data()), HistType::kTH2F, {observableAxisCluster, etaAxis});
    registry.add("hSelectedGammaDCALPtPhi", Form("Selected Gamma %s and #phi (DCAL, high threshold)", observableName.data()), HistType::kTH2F, {observableAxisCluster, phiAxis});
    registry.add("hSelectedGammaEMCALPtEtaPhi", Form("Selected Gamma %s, #eta and #phi (EMCAL, high threshold)", observableName.data()), HistType::kTH3F, {observableAxisCluster, etaAxis, phiAxis});
    registry.add("hSelectedGammaEMCALMaxPtEtaPhi", Form("Leading selected gamma %s, #eta and #phi (EMCAL, high treshold)", observableName.data()), HistType::kTH3F, {observableAxisCluster, etaAxis, phiAxis});
    registry.add("hSelectedGammaEMCALMaxPtEta", Form("Leading selected gamma %s and #eta (EMCAL, high threshold)", observableName.data()), HistType::kTH2F, {observableAxisCluster, etaAxis});
    registry.add("hSelectedGammaEMCALMaxPtPhi", Form("Leading selected gamma %s and #phi (EMCAL, high threshold)", observableName.data()), HistType::kTH2F, {observableAxisCluster, phiAxis});
    registry.add("hSelectedGammaDCALMaxPtEta", Form("Leading selected gamma %s and #eta (DCAL, high threshold)", observableName.data()), HistType::kTH2F, {observableAxisCluster, etaAxis});
    registry.add("hSelectedGammaDCALMaxPtPhi", Form("Leading selected gamma %s and #phi (DCAL, high threshold)", observableName.data()), HistType::kTH2F, {observableAxisCluster, phiAxis});
    registry.add("hSelectedGammaDCALPtEtaPhi", Form("Selected gamma %s, #eta and #phi (DCAL, high treshold)", observableName.data()), HistType::kTH3F, {observableAxisCluster, etaAxis, phiAxis});
    registry.add("hSelectedGammaDCALMaxPtEtaPhi", Form("Leading selected gamma %s, #eta and #phi (DCAL, high treshold)", observableName.data()), HistType::kTH3F, {observableAxisCluster, etaAxis, phiAxis});
    registry.add("hSelectedGammaEMCALPtEtaLow", Form("Selected Gamma %s and #eta (EMCAL, low threshold)", observableName.data()), HistType::kTH2F, {observableAxisCluster, etaAxis});
    registry.add("hSelectedGammaEMCALPtPhiLow", Form("Selected Gamma %s and #phi (EMCAL, low threshold)", observableName.data()), HistType::kTH2F, {observableAxisCluster, phiAxis});
    registry.add("hSelectedGammaDCALPtEtaLow", Form("Selected Gamma %s and #eta (DCAL, low threshold)", observableName.data()), HistType::kTH2F, {observableAxisCluster, etaAxis});
    registry.add("hSelectedGammaDCALPtPhiLow", Form("Selected Gamma %s and #phi (DCAL, low threshold)", observableName.data()), HistType::kTH2F, {observableAxisCluster, phiAxis});
    registry.add("hSelectedGammaEMCALPtEtaPhiLow", Form("Selected gamma %s, #eta and #phi (EMCAL, low threshold)", observableName.data()), HistType::kTH3F, {observableAxisCluster, etaAxis, phiAxis});
    registry.add("hSelectedGammaEMCALMaxPtEtaPhiLow", Form("Leading selected gamma %s, #eta and #phi (EMCAL, low threshold)", observableName.data()), HistType::kTH3F, {observableAxisCluster, etaAxis, phiAxis});
    registry.add("hSelectedGammaEMCALMaxPtEtaLow", Form("Leading selected gamma %s and #eta (EMCAL, low threshold)", observableName.data()), HistType::kTH2F, {observableAxisCluster, etaAxis});
    registry.add("hSelectedGammaEMCALMaxPtPhiLow", Form("Leading selected gamma %s and #phi (EMCAL, low threshold)", observableName.data()), HistType::kTH2F, {observableAxisCluster, phiAxis});
    registry.add("hSelectedGammaDCALMaxPtEtaLow", Form("Leading selected gamma %s and #eta (DCAL, low threshold)", observableName.data()), HistType::kTH2F, {observableAxisCluster, etaAxis});
    registry.add("hSelectedGammaDCALMaxPtPhiLow", Form("Leading selected gamma %s and #phi (DCAL, low threshold)", observableName.data()), HistType::kTH2F, {observableAxisCluster, phiAxis});
    registry.add("hSelectedGammaDCALPtEtaPhiLow", Form("Selected gamma %s, #eta and #phi (DCAL, low threshold)", observableName.data()), HistType::kTH3F, {observableAxisCluster, etaAxis, phiAxis});
    registry.add("hSelectedGammaDCALMaxPtEtaPhiLow", Form("Leading selected gamma %s, #eta and #phi (DCAL, low threshold)", observableName.data()), HistType::kTH3F, {observableAxisCluster, etaAxis, phiAxis});
    registry.add("hSelectedJetRMaxPtClusterMaxPt", "Leading selected jets and clusters", HistType::kTH3F, {rAxis, jetPtAxis, observableAxisCluster});

    registry.add("hJetRMaxPtJetPt", "Leading jet #it{p}_{T} vs jet #it{p}_{T}", HistType::kTH3F, {rAxis, jetMaxPtAxis, jetPtAxis});
    registry.add("hJetRMaxPtJetPtNoFiducial", "Leading jet #it{p}_{T} vs jet #it{p}_{T} (no fiducial cut)", HistType::kTH3F, {rAxis, jetMaxPtAxis, jetPtAxis});
    registry.add("hClusterEMCALMaxPtClusterEMCALPt", Form("Leading clusterEMCAL %s vs clusterEMCAL %s", observableName.data(), observableName.data()), HistType::kTH2F, {observableAxisMaxCluster, observableAxisCluster});
    registry.add("hClusterDCALMaxPtClusterDCALPt", Form("Leading clusterDCAL %s vs clusterDCAL %s", observableName.data(), observableName.data()), HistType::kTH2F, {observableAxisMaxCluster, observableAxisCluster});

    if (b_JetsInEmcalOnly) {
      registry.get<TH3>(HIST("hJetRPtEta"))->SetTitle("Jets (in emcal only) #it{p}_{T} and #eta");
      registry.get<TH3>(HIST("hJetRPtPhi"))->SetTitle("Jets (in emcal only) #it{p}_{T} and #phi");
      registry.get<TH3>(HIST("hJetRPtTrackPt"))->SetTitle("Jets (in emcal only)");
      registry.get<TH3>(HIST("hJetRPtClusterPt"))->SetTitle("Jets (in emcal only)");
      registry.get<TH3>(HIST("hJetRPtPtd"))->SetTitle("Jets (in emcal only)");
      registry.get<TH3>(HIST("hJetRPtChargeFrag"))->SetTitle("Jets (in emcal only)");
      registry.get<TH3>(HIST("hJetRPtNEF"))->SetTitle("Jets (in emcal only)");
      registry.get<TH3>(HIST("hJetRPtZTheta"))->SetTitle("Jets (in emcal only)");
      registry.get<TH3>(HIST("hJetRPtZSqTheta"))->SetTitle("Jets (in emcal only)");
      registry.get<TH3>(HIST("hJetRPtZThetaSq"))->SetTitle("Jets (in emcal only)");
      registry.get<TH3>(HIST("hJetRMaxPtEta"))->SetTitle("Leading jets (in emcal only) #it{p}_{T} and #eta");
      registry.get<TH3>(HIST("hJetRMaxPtPhi"))->SetTitle("Leading jets (in emcal only) #it{p}_{T} and #phi");
      registry.get<TH3>(HIST("hJetRMaxPtClusterMaxPt"))->SetTitle("Leading jets (in emcal only) and clusters");

      registry.get<TH3>(HIST("hSelectedJetRPtEta"))->SetTitle("Selected jets (in emcal only) #it{p}_{T} and #eta");
      registry.get<TH3>(HIST("hSelectedJetRPtPhi"))->SetTitle("Selected jets (in emcal only) #it{p}_{T} and #phi");
      registry.get<TH3>(HIST("hSelectedJetRPtTrackPt"))->SetTitle("Selected jets (in emcal only)");
      registry.get<TH3>(HIST("hSelectedJetRPtClusterPt"))->SetTitle("Selected jets (in emcal only)");
      registry.get<TH3>(HIST("hSelectedJetRPtPtd"))->SetTitle("Selected jets (in emcal only)");
      registry.get<TH3>(HIST("hSelectedJetRPtChargeFrag"))->SetTitle("Selected jets (in emcal only)");
      registry.get<TH3>(HIST("hSelectedJetRPtNEF"))->SetTitle("Selected jets (in emcal only)");
      registry.get<TH3>(HIST("hSelectedJetRPtZTheta"))->SetTitle("Selected jets (in emcal only)");
      registry.get<TH3>(HIST("hSelectedJetRPtZSqTheta"))->SetTitle("Selected jets (in emcal only)");
      registry.get<TH3>(HIST("hSelectedJetRPtZThetaSq"))->SetTitle("Selected jets (in emcal only)");
      registry.get<TH3>(HIST("hSelectedJetRMaxPtEta"))->SetTitle("Leading selected jets (in emcal only) #it{p}_{T} and #eta");
      registry.get<TH3>(HIST("hSelectedJetRMaxPtPhi"))->SetTitle("Leading selected jets (in emcal only) #it{p}_{T} and #phi");
      registry.get<TH3>(HIST("hSelectedJetRMaxPtClusterMaxPt"))->SetTitle("Leading selected jets (in emcal only) and clusters");

      registry.get<TH3>(HIST("hJetRMaxPtJetPt"))->SetTitle("Leading jet #it{p}_{T} vs jet #it{p}_{T} (in emcal only)");
    }
  } // init

  o2::aod::EMCALClusterDefinition clusDef = o2::aod::emcalcluster::getClusterDefinitionFromString(mClusterDefinition.value);
  Filter clusterDefinitionSelection = o2::aod::emcalcluster::definition == static_cast<int>(clusDef);

  double check_dphi(double dphi)
  {
    if (dphi > TMath::Pi()) {
      dphi -= TMath::Pi();
    } else if (dphi <= -1 * TMath::Pi()) {
      dphi += TMath::Pi();
    }
    return dphi;
  }

  template <typename T>
  Bool_t isJetInEmcal(T const& jet)
  {
    double emcalEtaMin = -0.7, emcalEtaMax = 0.7, emcalPhiMin = 1.40, emcalPhiMax = 3.26; // Phi: 80 - 187 deg
    double R = jet.r() * 1e-2;                                                            // Jet R is saved as round(100*R)
    if ((jet.eta() >= emcalEtaMin + R) && (jet.eta() <= emcalEtaMax - R) && (jet.phi() >= emcalPhiMin + R) && (jet.phi() <= emcalPhiMax - R)) {
      return true;
    }
    return false;
  }

  Bool_t isClusterInEmcal(selectedClusters::iterator const& cluster)
  {
    if (cluster.phi() < f_PhiEmcalOrDcal) {
      return true;
    }
    return false;
  }

  template <typename T, typename U>
  void check_maxJetPt(T const jet, U& vecMaxJet)
  {
    for (unsigned int i = 0; i < vecMaxJet.size(); i++) {
      auto maxJet = vecMaxJet[i];
      if (jet.r() != maxJet.r()) {
        continue;
      }
      if (jet.pt() > maxJet.pt()) {
        vecMaxJet[i] = jet;
      }
      return;
    }
    // If there is no jet of this radius yet, this is the leading one and we add it to the vector
    vecMaxJet.push_back(jet);
  }

  void fillEventSelectionCounter(int triggerID)
  {
    registry.fill(HIST("hProcessedEvents"), triggerID);
  }

  void fillCorrelationMatrix(TriggerType_t mainTrigger, TriggerType_t assocTrigger)
  {
    registry.fill(HIST("hTriggerCorrelation"), mainTrigger, assocTrigger);
  }

  void process(soa::Join<aod::Collisions, aod::EvSels, aod::FullJetFilters>::iterator const& collision,
               soa::Join<aod::FullJets, aod::FullJetConstituents> const& jets,
               aod::Tracks const& tracks,
               selectedClusters const& clusters)
  {
    std::bitset<TriggerType_t::kNTriggers> triggerstatus;
    auto isTrigger = [&triggerstatus](TriggerType_t triggertype) -> bool {
      return triggerstatus.test(triggertype);
    };
    auto setTrigger = [&triggerstatus](TriggerType_t triggertype) {
      triggerstatus.set(triggertype, true);
    };

    fillEventSelectionCounter(0);
    setTrigger(TriggerType_t::kMinBias);

    // fill event counters and correlation matrix without constraint on the EMCAL trigger flag
    if (collision.alias()[kTVXinEMC]) {
      fillEventSelectionCounter(1);
      setTrigger(TriggerType_t::kEmcal);
    }

    bool isEvtSelected = false, isEventSelectedNeutralJet = false, isEvtSelectedGammaEMCAL = false, isEvtSelectedGammaDCAL = false, isEvtSelectedGammaLowEMCAL = false, isEvtSelectedGammaLowDCAL = false;
    if (collision.hasJetFullHighPt()) {
      isEvtSelected = true;
      fillEventSelectionCounter(2);
      setTrigger(TriggerType_t::kEmcalJetFull);
    }
    if (collision.hasJetNeutralHighPt()) {
      isEventSelectedNeutralJet = true;
      fillEventSelectionCounter(3);
      setTrigger(TriggerType_t::kEmcalJetNeutral);
    }
    if (collision.hasGammaHighPtEMCAL()) {
      isEvtSelectedGammaEMCAL = true;
      fillEventSelectionCounter(4);
      setTrigger(TriggerType_t::kEmcalGammaHigh);
    }
    if (collision.hasGammaHighPtDCAL()) {
      isEvtSelectedGammaDCAL = true;
      fillEventSelectionCounter(5);
      setTrigger(TriggerType_t::kDcalGammaHigh);
    }
    if (collision.hasGammaLowPtEMCAL()) {
      isEvtSelectedGammaLowEMCAL = true;
      fillEventSelectionCounter(6);
      setTrigger(TriggerType_t::kEmcalGammaLow);
    }
    if (collision.hasGammaLowPtDCAL()) {
      isEvtSelectedGammaLowDCAL = true;
      fillEventSelectionCounter(7);
      setTrigger(TriggerType_t::kDcalGammaLow);
    }
    if (collision.hasJetFullHighPt() && collision.hasGammaHighPtEMCAL()) {
      fillEventSelectionCounter(8);
    }
    if (collision.hasJetFullHighPt() && collision.hasGammaHighPtDCAL()) {
      fillEventSelectionCounter(9);
    }
    if (collision.hasJetNeutralHighPt() && collision.hasGammaHighPtEMCAL()) {
      fillEventSelectionCounter(10);
    }
    if (collision.hasJetNeutralHighPt() && collision.hasGammaHighPtDCAL()) {
      fillEventSelectionCounter(11);
    }
    if (collision.hasJetFullHighPt() && collision.hasJetNeutralHighPt()) {
      fillEventSelectionCounter(12);
    }
    if (collision.hasGammaHighPtEMCAL() && collision.hasGammaHighPtDCAL()) {
      fillEventSelectionCounter(13);
    }

    // Create correlationMatrix
    for (std::size_t maintriggertype = 0; maintriggertype < TriggerType_t::kNTriggers; maintriggertype++) {
      if (isTrigger(static_cast<TriggerType_t>(maintriggertype))) {
        for (std::size_t assoctriggertype = 0; assoctriggertype < TriggerType_t::kNTriggers; assoctriggertype++) {
          if (isTrigger(static_cast<TriggerType_t>(assoctriggertype))) {
            fillCorrelationMatrix(static_cast<TriggerType_t>(maintriggertype), static_cast<TriggerType_t>(assoctriggertype));
          }
        }
      }
    }

    // Discard collisions without EMCAL if it is not respecifically required to discard the the EMCAL flag
    // If ignore is true, we don't check for the flag and accept all events
    if (!b_IgnoreEmcalFlag) {
      if (!collision.alias()[kTVXinEMC]) {
        return; // Only consider events where EMCAL is live
      }
    }

    double maxClusterObservableEMCAL = -1., maxClusterObservableDCAL = -1.;
    selectedClusters::iterator maxClusterEMCAL, maxClusterDCAL;
    std::vector<soa::Join<aod::FullJets, aod::FullJetConstituents>::iterator> vecMaxJet;
    std::vector<soa::Join<aod::FullJets, aod::FullJetConstituents>::iterator> vecMaxJetNoFiducial;

    for (const auto& jet : jets) {
      double jetPt = jet.pt(), jetR = jet.r() * 1e-2;
      registry.get<THn>(HIST("jetRPtEtaPhiNoFiducial"))->Fill(jetR, jetPt, jet.eta(), jet.phi());
      check_maxJetPt(jet, vecMaxJetNoFiducial);
      if (b_JetsInEmcalOnly && !isJetInEmcal(jet)) {
        continue;
      }
      float neutralEnergyFraction = 0., ptD = 0.;
      float zTheta = 0., zSqTheta = 0., zThetaSq = 0; // Jet angularities (1,1), (1,2), (2,1)
      check_maxJetPt(jet, vecMaxJet);

      // Slice JetTrackConstituents table to obtain a list of trackId belonging to this jet only (and the same for clusters)
      // This gives us access to all jet substructure information
      // auto tracksInJet = jetTrackConstituents.sliceBy(perJetTrackConstituents, jet.globalIndex());
      // for (const auto& trackList : tracksInJet) {
      for (auto& track : jet.tracks_as<aod::Tracks>()) {
        auto trackPt = track.pt();
        auto chargeFrag = track.px() * jet.px() + track.py() * jet.py() + track.pz() * jet.pz();
        chargeFrag /= (jet.p() * jet.p());
        auto z = trackPt / jetPt;
        auto dphi = jet.phi() - track.phi();
        dphi = check_dphi(dphi);
        auto dR = TMath::Sqrt(dphi * dphi + TMath::Power(jet.eta() - track.eta(), 2));
        zTheta += z * dR / jetR;
        zSqTheta += z * z * dR / jetR;
        zThetaSq += z * dR * dR / (jetR * jetR);
        ptD += z * z;
        registry.fill(HIST("hJetRPtTrackPt"), jetR, jetPt, trackPt);
        registry.fill(HIST("hJetRPtChargeFrag"), jetR, jetPt, chargeFrag);
        if (isEvtSelected || isEventSelectedNeutralJet) {
          registry.fill(HIST("hSelectedJetRPtTrackPt"), jetR, jetPt, trackPt);
          registry.fill(HIST("hSelectedJetRPtChargeFrag"), jetR, jetPt, chargeFrag);
        }
      } // for tracks in jet

      // auto clustersInJet = jetClusterConstituents.sliceBy(perJetClusterConstituents, jet.globalIndex());
      // for (const auto& clusterList : clustersInJet) {
      for (auto& cluster : jet.tracks_as<selectedClusters>()) {
        auto clusterPt = cluster.energy() / std::cosh(cluster.eta());
        neutralEnergyFraction += cluster.energy();
        auto z = clusterPt / jetPt;
        auto dphi = jet.phi() - cluster.phi();
        dphi = check_dphi(dphi);
        auto dR = TMath::Sqrt(dphi * dphi + TMath::Power(jet.eta() - cluster.eta(), 2));
        zTheta += z * dR / jetR;
        zSqTheta += z * z * dR / jetR;
        zThetaSq += z * dR * dR / (jetR * jetR);
        ptD += z * z;
        registry.fill(HIST("hJetRPtClusterPt"), jetR, jetPt, clusterPt);
        if (isEvtSelected || isEventSelectedNeutralJet) {
          registry.fill(HIST("hSelectedJetRPtClusterPt"), jetR, jetPt, clusterPt);
        }
      } // for clusters in jet
      neutralEnergyFraction /= jet.energy();
      ptD = TMath::Sqrt(ptD);

      // Fillng histograms
      registry.fill(HIST("hJetRPtEta"), jetR, jetPt, jet.eta());
      registry.fill(HIST("hJetRPtPhi"), jetR, jetPt, jet.phi());
      registry.fill(HIST("hJetRPtPtd"), jetR, jetPt, ptD);
      registry.fill(HIST("hJetRPtNEF"), jetR, jetPt, neutralEnergyFraction);
      registry.fill(HIST("hJetRPtZTheta"), jetR, jetPt, zTheta);
      registry.fill(HIST("hJetRPtZSqTheta"), jetR, jetPt, zSqTheta);
      registry.fill(HIST("hJetRPtZThetaSq"), jetR, jetPt, zThetaSq);
      registry.get<THn>(HIST("jetRPtEtaPhi"))->Fill(jetR, jetPt, jet.eta(), jet.phi());
      if (isEvtSelected || isEventSelectedNeutralJet) {
        registry.fill(HIST("hSelectedJetRPtEta"), jetR, jetPt, jet.eta());
        registry.fill(HIST("hSelectedJetRPtPhi"), jetR, jetPt, jet.phi());
        registry.fill(HIST("hSelectedJetRPtPtd"), jetR, jetPt, ptD);
        registry.fill(HIST("hSelectedJetRPtNEF"), jetR, jetPt, neutralEnergyFraction);
        registry.fill(HIST("hSelectedJetRPtZTheta"), jetR, jetPt, zTheta);
        registry.fill(HIST("hSelectedJetRPtZSqTheta"), jetR, jetPt, zSqTheta);
        registry.fill(HIST("hSelectedJetRPtZThetaSq"), jetR, jetPt, zThetaSq);
      }
    } // for jets

    struct ClusterData {
      float mTriggerObservable;
      float mEta;
      float mPhi;
      bool mEMCALcluster;
    };
    std::vector<ClusterData> analysedClusters;

    for (const auto& cluster : clusters) {
      if (b_RejectExoticClusters && cluster.isExotic()) {
        continue;
      }
      if (cluster.time() < f_minClusterTime || cluster.time() > f_maxClusterTime) {
        continue;
      }
      bool emcalCluster = isClusterInEmcal(cluster);
      double clusterObservable = (f_GammaObservable == 0) ? cluster.energy() : cluster.energy() / std::cosh(cluster.eta());
      analysedClusters.push_back({static_cast<float>(clusterObservable), cluster.eta(), cluster.phi(), emcalCluster});

      if (emcalCluster && (clusterObservable > maxClusterObservableEMCAL)) {
        maxClusterObservableEMCAL = clusterObservable;
        maxClusterEMCAL = cluster;
      }
      if (!emcalCluster && (clusterObservable > maxClusterObservableDCAL)) {
        maxClusterObservableDCAL = clusterObservable;
        maxClusterDCAL = cluster;
      }
      registry.fill(HIST("hClusterPtEta"), clusterObservable, cluster.eta());
      registry.fill(HIST("hClusterPtPhi"), clusterObservable, cluster.phi());
      registry.fill(HIST("hClusterPtEtaPhi"), clusterObservable, cluster.eta(), cluster.phi());
      if (isEvtSelected || isEventSelectedNeutralJet) {
        registry.fill(HIST("hSelectedClusterPtEta"), clusterObservable, cluster.eta());
        registry.fill(HIST("hSelectedClusterPtPhi"), clusterObservable, cluster.phi());
        registry.fill(HIST("hSelectedClusterPtEtaPhi"), clusterObservable, cluster.eta(), cluster.phi());
      }
      if (isEvtSelectedGammaEMCAL && isClusterInEmcal(cluster)) { // Only fill EMCAL clusters
        registry.fill(HIST("hSelectedGammaEMCALPtEta"), clusterObservable, cluster.eta());
        registry.fill(HIST("hSelectedGammaEMCALPtPhi"), clusterObservable, cluster.phi());
        registry.fill(HIST("hSelectedGammaEMCALPtEtaPhi"), clusterObservable, cluster.eta(), cluster.phi());
      }
      if (isEvtSelectedGammaDCAL && !isClusterInEmcal(cluster)) { // Only fill DCAL clusters
        registry.fill(HIST("hSelectedGammaDCALPtEta"), clusterObservable, cluster.eta());
        registry.fill(HIST("hSelectedGammaDCALPtPhi"), clusterObservable, cluster.phi());
        registry.fill(HIST("hSelectedGammaDCALPtEtaPhi"), clusterObservable, cluster.eta(), cluster.phi());
      }
      if (isEvtSelectedGammaLowEMCAL && isClusterInEmcal(cluster)) { // Only fill EMCAL clusters
        registry.fill(HIST("hSelectedGammaEMCALPtEtaLow"), clusterObservable, cluster.eta());
        registry.fill(HIST("hSelectedGammaEMCALPtPhiLow"), clusterObservable, cluster.phi());
        registry.fill(HIST("hSelectedGammaEMCALPtEtaPhiLow"), clusterObservable, cluster.eta(), cluster.phi());
      }
      if (isEvtSelectedGammaLowDCAL && !isClusterInEmcal(cluster)) { // Only fill DCAL clusters
        registry.fill(HIST("hSelectedGammaDCALPtEtaLow"), clusterObservable, cluster.eta());
        registry.fill(HIST("hSelectedGammaDCALPtPhiLow"), clusterObservable, cluster.phi());
        registry.fill(HIST("hSelectedGammaDCALPtEtaPhiLow"), clusterObservable, cluster.eta(), cluster.phi());
      }
    } // for clusters

    if (maxClusterObservableEMCAL > 0) {
      // hClusterMaxPtEta->Fill(maxClusterPt, maxCluster.eta());
      // hClusterMaxPtPhi->Fill(maxClusterPt, maxCluster.phi());
      registry.fill(HIST("hClusterEMCALMaxPtEtaPhi"), maxClusterObservableEMCAL, maxClusterEMCAL.eta(), maxClusterEMCAL.phi());
      for (const auto& cluster : analysedClusters) {
        if (cluster.mEMCALcluster) {
          registry.fill(HIST("hClusterEMCALMaxPtClusterEMCALPt"), maxClusterObservableEMCAL, cluster.mTriggerObservable);
        }
      }
      if (isEvtSelected) {
        registry.fill(HIST("hSelectedClusterMaxPtEta"), maxClusterObservableEMCAL, maxClusterEMCAL.eta());
        registry.fill(HIST("hSelectedClusterMaxPtPhi"), maxClusterObservableEMCAL, maxClusterEMCAL.phi());
        registry.fill(HIST("hSelectedClusterMaxPtEtaPhi"), maxClusterObservableEMCAL, maxClusterEMCAL.eta(), maxClusterEMCAL.phi());
      }
      if (isEvtSelectedGammaEMCAL) {
        registry.fill(HIST("hSelectedGammaEMCALMaxPtEta"), maxClusterObservableEMCAL, maxClusterEMCAL.eta());
        registry.fill(HIST("hSelectedGammaEMCALMaxPtPhi"), maxClusterObservableEMCAL, maxClusterEMCAL.phi());
        registry.fill(HIST("hSelectedGammaEMCALMaxPtEtaPhi"), maxClusterObservableEMCAL, maxClusterEMCAL.eta(), maxClusterEMCAL.phi());
      }
      if (isEvtSelectedGammaLowEMCAL) {
        registry.fill(HIST("hSelectedGammaEMCALMaxPtEtaLow"), maxClusterObservableEMCAL, maxClusterEMCAL.eta());
        registry.fill(HIST("hSelectedGammaEMCALMaxPtPhiLow"), maxClusterObservableEMCAL, maxClusterEMCAL.phi());
        registry.fill(HIST("hSelectedGammaEMCALMaxPtEtaPhiLow"), maxClusterObservableEMCAL, maxClusterEMCAL.eta(), maxClusterEMCAL.phi());
      }
    }
    if (maxClusterObservableDCAL > 0) {
      registry.fill(HIST("hClusterDCALMaxPtEtaPhi"), maxClusterObservableDCAL, maxClusterDCAL.eta(), maxClusterDCAL.phi());
      for (const auto& cluster : analysedClusters) {
        if (!cluster.mEMCALcluster) {
          registry.fill(HIST("hClusterDCALMaxPtClusterDCALPt"), maxClusterObservableDCAL, cluster.mTriggerObservable);
        }
      }
      if (isEvtSelectedGammaDCAL) {
        registry.fill(HIST("hSelectedGammaDCALMaxPtEta"), maxClusterObservableDCAL, maxClusterDCAL.eta());
        registry.fill(HIST("hSelectedGammaDCALMaxPtPhi"), maxClusterObservableDCAL, maxClusterDCAL.phi());
        registry.fill(HIST("hSelectedGammaDCALMaxPtEtaPhi"), maxClusterObservableDCAL, maxClusterDCAL.eta(), maxClusterDCAL.phi());
      }
      if (isEvtSelectedGammaLowDCAL) {
        registry.fill(HIST("hSelectedGammaDCALMaxPtEtaLow"), maxClusterObservableDCAL, maxClusterDCAL.eta());
        registry.fill(HIST("hSelectedGammaDCALMaxPtPhiLow"), maxClusterObservableDCAL, maxClusterDCAL.phi());
        registry.fill(HIST("hSelectedGammaDCALMaxPtEtaPhiLow"), maxClusterObservableDCAL, maxClusterDCAL.eta(), maxClusterDCAL.phi());
      }
    }

    for (auto maxJet : vecMaxJet) {
      double jetR = maxJet.r() * 1e-2, jetPt = maxJet.pt(), jetEta = maxJet.eta(), jetPhi = maxJet.phi();
      registry.fill(HIST("hJetRMaxPtEta"), jetR, jetPt, jetEta);
      registry.fill(HIST("hJetRMaxPtPhi"), jetR, jetPt, jetPhi);
      // hJetRMaxPtEtaPhi->Fill(jetR, jetPt, jetEta, jetPhi);
      registry.get<THn>(HIST("jetRMaxPtEtaPhi"))->Fill(jetR, jetPt, jetEta, jetPhi);
      if (isEvtSelected) {
        registry.fill(HIST("hSelectedJetRMaxPtEta"), jetR, jetPt, jetEta);
        registry.fill(HIST("hSelectedJetRMaxPtPhi"), jetR, jetPt, jetPhi);
      }
      if (maxClusterObservableEMCAL > 0) {
        registry.fill(HIST("hJetRMaxPtClusterMaxPt"), jetR, jetPt, maxClusterObservableEMCAL);
        if (isEvtSelected) {
          registry.fill(HIST("hSelectedJetRMaxPtClusterMaxPt"), jetR, jetPt, maxClusterObservableEMCAL);
        }
      } // if maxClusterPt
      if (maxJet.r() == std::round(f_jetR * 100)) {
        for (const auto& jet : jets) {
          if (isJetInEmcal(jet)) {
            registry.fill(HIST("hJetRMaxPtJetPt"), jet.r() * 1e-2, jetPt, jet.pt());
          }
        } // for jets
      }   // if maxJet.r() == std::round(f_jetR * 100)
    }     // for maxJet
    for (auto maxJet : vecMaxJetNoFiducial) {
      double jetR = maxJet.r() * 1e-2, jetPt = maxJet.pt(), jetEta = maxJet.eta(), jetPhi = maxJet.phi();
      // hJetRMaxPtEtaPhiNoFiducial->Fill(jetR, jetPt, jetEta, jetPhi);
      registry.get<THn>(HIST("jetRMaxPtEtaPhiNoFiducial"))->Fill(jetR, jetPt, jetEta, jetPhi);
      if (maxJet.r() == std::round(f_jetR * 100)) {
        for (const auto& jet : jets) {
          registry.fill(HIST("hJetRMaxPtJetPtNoFiducial"), jet.r() * 1e-2, jetPt, jet.pt());
        } // for jets
      }   // if maxJet.r() == std::round(f_jetR * 100)
    }     // for maxjet no fiducial
  }       // process
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<JetTriggerQA>(cfgc, TaskName{"jet-full-trigger-qa"})};
}
