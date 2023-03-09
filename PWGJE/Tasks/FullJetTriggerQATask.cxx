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
  Preslice<aod::JetTrackConstituents> perJetTrackConstituents = o2::aod::jetconstituents::jetId;
  Preslice<aod::JetClusterConstituents> perJetClusterConstituents = o2::aod::jetconstituents::jetId;

  HistogramRegistry registry{"registry"};

  OutputObj<TH1I> hProcessedEvents{"hProcessedEvents"};

  OutputObj<TH3F> hJetRPtEta{"hJetRPtEta"};
  OutputObj<TH3F> hJetRPtPhi{"hJetRPtPhi"};
  OutputObj<TH3F> hJetRMaxPtEta{"hJetRMaxPtEta"};
  OutputObj<TH3F> hJetRMaxPtPhi{"hJetRMaxPtPhi"};
  OutputObj<TH2F> hClusterPtEta{"hClusterPtEta"};
  OutputObj<TH2F> hClusterPtPhi{"hClusterPtPhi"};
  OutputObj<TH3F> hClusterPtEtaPhi{"hClusterPtEtaPhi"};
  OutputObj<TH3F> hClusterEMCALMaxPtEtaPhi{"hClusterEMCALMaxPtEtaPhi"};
  OutputObj<TH3F> hClusterDCALMaxPtEtaPhi{"hClusterDCALMaxPtEtaPhi"};
  OutputObj<TH2F> hClusterMaxPtEta{"hClusterMaxPtEta"};
  OutputObj<TH2F> hClusterMaxPtPhi{"hClusterMaxPtPhi"};
  OutputObj<TH3F> hJetRPtTrackPt{"hJetRPtTrackPt"};
  OutputObj<TH3F> hJetRPtClusterPt{"hJetRPtClusterPt"};
  OutputObj<TH3F> hJetRPtPtd{"hJetRPtPtd"};
  OutputObj<TH3F> hJetRPtChargeFrag{"hJetRPtChargeFrag"};
  OutputObj<TH3F> hJetRPtNEF{"hJetRPtNEF"};
  OutputObj<TH3F> hJetRPtZTheta{"hJetRPtZTheta"};
  OutputObj<TH3F> hJetRPtZSqTheta{"hJetRPtZSqTheta"};
  OutputObj<TH3F> hJetRPtZThetaSq{"hJetRPtZThetaSq"};
  OutputObj<TH3F> hJetRMaxPtClusterMaxPt{"hJetRMaxPtClusterMaxPt"};

  OutputObj<TH3F> hSelectedJetRPtEta{"hSelectedJetRPtEta"};
  OutputObj<TH3F> hSelectedJetRPtPhi{"hSelectedJetRPtPhi"};
  OutputObj<TH3F> hSelectedJetRMaxPtEta{"hSelectedJetRMaxPtEta"};
  OutputObj<TH3F> hSelectedJetRMaxPtPhi{"hSelectedJetRMaxPtPhi"};
  OutputObj<TH2F> hSelectedClusterPtEta{"hSelectedClusterPtEta"};
  OutputObj<TH2F> hSelectedClusterPtPhi{"hSelectedClusterPtPhi"};
  OutputObj<TH3F> hSelectedClusterPtEtaPhi{"hSelectedClusterPtEtaPhi"};
  OutputObj<TH3F> hSelectedClusterMaxPtEtaPhi{"hSelectedClusterMaxPtEtaPhi"};
  OutputObj<TH2F> hSelectedClusterMaxPtEta{"hSelectedClusterMaxPtEta"};
  OutputObj<TH2F> hSelectedClusterMaxPtPhi{"hSelectedClusterMaxPtPhi"};
  OutputObj<TH2F> hSelectedClusterPtEtaLow{"hSelectedClusterPtEtaLow"};
  OutputObj<TH2F> hSelectedClusterPtPhiLow{"hSelectedClusterPtPhiLow"};
  OutputObj<TH3F> hSelectedClusterPtEtaPhiLow{"hSelectedClusterPtEtaPhiLow"};
  OutputObj<TH3F> hSelectedClusterMaxPtEtaPhiLow{"hSelectedClusterMaxPtEtaPhiLow"};
  OutputObj<TH2F> hSelectedClusterMaxPtEtaLow{"hSelectedClusterMaxPtEtaLow"};
  OutputObj<TH2F> hSelectedClusterMaxPtPhiLow{"hSelectedClusterMaxPtPhiLow"};
  OutputObj<TH3F> hSelectedJetRPtTrackPt{"hSelectedJetRPtTrackPt"};
  OutputObj<TH3F> hSelectedJetRPtClusterPt{"hSelectedJetRPtClusterPt"};
  OutputObj<TH3F> hSelectedJetRPtPtd{"hSelectedJetRPtPtd"};
  OutputObj<TH3F> hSelectedJetRPtChargeFrag{"hSelectedJetRPtChargeFrag"};
  OutputObj<TH3F> hSelectedJetRPtNEF{"hSelectedJetRPtNEF"};
  OutputObj<TH3F> hSelectedJetRPtZTheta{"hSelectedJetRPtZTheta"};
  OutputObj<TH3F> hSelectedJetRPtZSqTheta{"hSelectedJetRPtZSqTheta"};
  OutputObj<TH3F> hSelectedJetRPtZThetaSq{"hSelectedJetRPtZThetaSq"};
  OutputObj<TH3F> hSelectedJetRMaxPtClusterMaxPt{"hSelectedJetRMaxPtClusterMaxPt"};

  OutputObj<TH2F> hSelectedGammaEMCALPtEta{"hSelectedGammaEMCALPtEta"};
  OutputObj<TH2F> hSelectedGammaEMCALPtPhi{"hSelectedGammaEMCALPtPhi"};
  OutputObj<TH2F> hSelectedGammaEMCALMaxPtEta{"hSelectedGammaEMCALMaxPtEta"};
  OutputObj<TH2F> hSelectedGammaEMCALMaxPtPhi{"hSelectedGammaEMCALMaxPtPhi"};
  OutputObj<TH2F> hSelectedGammaDCALPtEta{"hSelectedGammaDCALPtEta"};
  OutputObj<TH2F> hSelectedGammaDCALPtPhi{"hSelectedGammaDCALPtPhi"};
  OutputObj<TH2F> hSelectedGammaDCALMaxPtEta{"hSelectedGammaDCALMaxPtEta"};
  OutputObj<TH2F> hSelectedGammaDCALMaxPtPhi{"hSelectedGammaDCALMaxPtPhi"};

  OutputObj<TH3F> hSelectedGammaEMCALPtEtaPhi{"hSelectedGammaEMCALPtEtaPhi"};
  OutputObj<TH3F> hSelectedGammaEMCALMaxPtEtaPhi{"hSelectedGammaEMCALMaxPtEtaPhi"};
  OutputObj<TH3F> hSelectedGammaDCALPtEtaPhi{"hSelectedGammaDCALPtEtaPhi"};
  OutputObj<TH3F> hSelectedGammaDCALMaxPtEtaPhi{"hSelectedGammaDCALMaxPtEtaPhi"};

  OutputObj<TH2F> hSelectedGammaEMCALPtEtaLow{"hSelectedGammaEMCALPtEtaLow"};
  OutputObj<TH2F> hSelectedGammaEMCALPtPhiLow{"hSelectedGammaEMCALPtPhiLow"};
  OutputObj<TH2F> hSelectedGammaEMCALMaxPtEtaLow{"hSelectedGammaEMCALMaxPtEtaLow"};
  OutputObj<TH2F> hSelectedGammaEMCALMaxPtPhiLow{"hSelectedGammaEMCALMaxPtPhiLow"};
  OutputObj<TH2F> hSelectedGammaDCALPtEtaLow{"hSelectedGammaDCALPtEtaLow"};
  OutputObj<TH2F> hSelectedGammaDCALPtPhiLow{"hSelectedGammaDCALPtPhiLow"};
  OutputObj<TH2F> hSelectedGammaDCALMaxPtEtaLow{"hSelectedGammaDCALMaxPtEtaLow"};
  OutputObj<TH2F> hSelectedGammaDCALMaxPtPhiLow{"hSelectedGammaDCALMaxPtPhiLow"};

  OutputObj<TH3F> hSelectedGammaEMCALPtEtaPhiLow{"hSelectedGammaEMCALPtEtaPhiLow"};
  OutputObj<TH3F> hSelectedGammaEMCALMaxPtEtaPhiLow{"hSelectedGammaEMCALMaxPtEtaPhiLow"};
  OutputObj<TH3F> hSelectedGammaDCALPtEtaPhiLow{"hSelectedGammaDCALPtEtaPhiLow"};
  OutputObj<TH3F> hSelectedGammaDCALMaxPtEtaPhiLow{"hSelectedGammaDCALMaxPtEtaPhiLow"};

  OutputObj<TH3F> hJetRMaxPtJetPt{"hJetRMaxPtJetPt"};
  OutputObj<TH3F> hJetRMaxPtJetPtNoFiducial{"hJetRMaxPtJetPtNoFiducial"};
  OutputObj<TH2F> hClusterEMCALMaxPtClusterEMCALPt{"hClusterEMCALMaxPtClusterEMCALPt"};
  OutputObj<TH2F> hClusterDCALMaxPtClusterDCALPt{"hClusterDCALMaxPtClusterDCALPt"};

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

    AxisSpec rAxis = {nRBins, kMinR, kMaxR, "#it{R}"};
    AxisSpec jetPtAxis = {nPtBins, kMinPt, kMaxPt, "#it{p}_{T}^{jet}"};
    AxisSpec etaAxis = {nEtaBins, kMinEta, kMaxEta, "#eta"};
    AxisSpec phiAxis = {nPhiBins, kMinPhi, kMaxPhi, "#phi"};

    HistogramConfigSpec hJetRPtEtaPhi({HistType::kTHnF, {rAxis, jetPtAxis, etaAxis, phiAxis}});
    HistogramConfigSpec hJetRMaxPtEtaPhi({HistType::kTHnF, {rAxis, jetPtAxis, etaAxis, phiAxis}});
    HistogramConfigSpec hJetRPtEtaPhiNoFiducial({HistType::kTHnF, {rAxis, jetPtAxis, etaAxis, phiAxis}});
    HistogramConfigSpec hJetRMaxPtEtaPhiNoFiducial({HistType::kTHnF, {rAxis, jetPtAxis, etaAxis, phiAxis}});

    registry.add("jetRPtEtaPhi", "JetRPtEtaPhi", hJetRPtEtaPhi);
    registry.add("jetRMaxPtEtaPhi", "JetRMaxPtEtaPhi", hJetRMaxPtEtaPhi);
    registry.add("jetRPtEtaPhiNoFiducial", "JetRPtEtaPhiNoFiducial", hJetRPtEtaPhiNoFiducial);
    registry.add("jetRMaxPtEtaPhiNoFiducial", "JetRMaxPtEtaPhiNoFiducial", hJetRMaxPtEtaPhiNoFiducial);

    hProcessedEvents.setObject(new TH1I("hProcessedEvents", "Processed events", 14, -0.5, 13.5));
    hProcessedEvents->GetXaxis()->SetBinLabel(1, "MB");
    hProcessedEvents->GetXaxis()->SetBinLabel(2, "EMC");
    hProcessedEvents->GetXaxis()->SetBinLabel(3, "Selected Full Jet");
    hProcessedEvents->GetXaxis()->SetBinLabel(4, "Selected Neutral Jet");
    hProcessedEvents->GetXaxis()->SetBinLabel(5, "Selected Gamma high EMCAL");
    hProcessedEvents->GetXaxis()->SetBinLabel(6, "Selected Gamma high DCAL");
    hProcessedEvents->GetXaxis()->SetBinLabel(7, "Selected Gamma low EMCAL");
    hProcessedEvents->GetXaxis()->SetBinLabel(8, "Selected Gamma low DCAL");
    hProcessedEvents->GetXaxis()->SetBinLabel(9, "Selected full jet and Gamma high EMCAL");
    hProcessedEvents->GetXaxis()->SetBinLabel(10, "Selected full jet and Gamma high DCAL");
    hProcessedEvents->GetXaxis()->SetBinLabel(11, "Selected neutral jet and Gamma high EMCAL");
    hProcessedEvents->GetXaxis()->SetBinLabel(12, "Selected neutral jet and Gamma high DCAL");
    hProcessedEvents->GetXaxis()->SetBinLabel(13, "Selected Gamma high EMCAL and high Gamma DCAL");
    hProcessedEvents->GetXaxis()->SetBinLabel(14, "Selected full jet and neutral jet");

    // Histograms for events where the EMCAL is live
    std::string observableName = (f_GammaObservable == 0) ? "Energy" : "#it{p}_{T}",
                observableAxisTitle = (f_GammaObservable == 0) ? "E (GeV)" : "#it{p}_{T} (GeV/c)",
                observableAxisTitleMax = (f_GammaObservable == 0) ? "E^{max} (GeV)" : "#it{p}_{T}^{max} (GeV/c)";
    hJetRPtEta.setObject(new TH3F("hJetRPtEta", "Jets #it{p}_{T} and #eta;#it{R};#it{p}_{T};#eta", nRBins, kMinR, kMaxR, nPtBins, kMinPt, kMaxPt, nEtaBins, kMinEta, kMaxEta));
    hJetRPtPhi.setObject(new TH3F("hJetRPtPhi", "Jets #it{p}_{T} and #phi;#it{R};#it{p}_{T};#phi", nRBins, kMinR, kMaxR, nPtBins, kMinPt, kMaxPt, nPhiBins, kMinPhi, kMaxPhi));
    hJetRMaxPtEta.setObject(new TH3F("hJetRMaxPtEta", "Leading jets #it{p}_{T} and #eta;#it{R};#it{p}_{T};#eta", nRBins, kMinR, kMaxR, nPtBins, kMinPt, kMaxPt, nEtaBins, kMinEta, kMaxEta));
    hJetRMaxPtPhi.setObject(new TH3F("hJetRMaxPtPhi", "Leading jets #it{p}_{T} and #phi;#it{R};#it{p}_{T};#phi", nRBins, kMinR, kMaxR, nPtBins, kMinPt, kMaxPt, nPhiBins, kMinPhi, kMaxPhi));
    hClusterPtEta.setObject(new TH2F("hClusterPtEta", Form("Cluster %s and #eta;%s;#eta", observableName.data(), observableAxisTitle.data()), nPtBins, kMinPt, kMaxPt / 2, nEtaBins, kMinEta, kMaxEta));
    hClusterPtPhi.setObject(new TH2F("hClusterPtPhi", Form("Cluster %s and #phi;%s;#phi", observableName.data(), observableAxisTitle.data()), nPtBins, kMinPt, kMaxPt / 2, nPhiBins, kMinPhi, kMaxPhi));
    hClusterPtEtaPhi.setObject(new TH3F("hClusterPtEtaPhi", Form("Cluster %s, #eta and #phi;%s;#eta;#phi", observableName.data(), observableAxisTitle.data()), nPtBins, kMinPt, kMaxPt / 2, nEtaBins, kMinEta, kMaxEta, nPhiBins, kMinPhi, kMaxPhi));
    hClusterEMCALMaxPtEtaPhi.setObject(new TH3F("hClusterEMCALMaxPtEtaPhi", Form("Leading clusterEMCAL %s, #eta and #phi;%s;#eta;#phi", observableName.data(), observableAxisTitle.data()), nPtBins, kMinPt, kMaxPt / 2, nEtaBins, kMinEta, kMaxEta, nPhiBins, kMinPhi, kMaxPhi));
    hClusterDCALMaxPtEtaPhi.setObject(new TH3F("hClusterDCALMaxPtEtaPhi", Form("Leading clusterDCAL %s, #eta and #phi;%s;#eta;#phi", observableName.data(), observableAxisTitle.data()), nPtBins, kMinPt, kMaxPt / 2, nEtaBins, kMinEta, kMaxEta, nPhiBins, kMinPhi, kMaxPhi));
    hClusterMaxPtEta.setObject(new TH2F("hClusterMaxPtEta", Form("Leading clusters %s and #eta;%s;#eta", observableName.data(), observableAxisTitle.data()), nPtBins, kMinPt, kMaxPt / 2, nEtaBins, kMinEta, kMaxEta));
    hClusterMaxPtPhi.setObject(new TH2F("hClusterMaxPtPhi", Form("Leading clusters %s and #phi;%s;#phi", observableName.data(), observableAxisTitle.data()), nPtBins, kMinPt, kMaxPt / 2, nPhiBins, kMinPhi, kMaxPhi));
    hJetRPtTrackPt.setObject(new TH3F("hJetRPtTrackPt", "Jets;#it{p}_{T};#it{R};#it{p}_{T}^{track}", nRBins, kMinR, kMaxR, nPtBins, kMinPt, kMaxPt, nPtBins, kMinPt, kMaxPt / 2));
    hJetRPtClusterPt.setObject(new TH3F("hJetRPtClusterPt", "Jets;#it{R};#it{p}_{T};#it{p}_{T}^{clus}", nRBins, kMinR, kMaxR, nPtBins, kMinPt, kMaxPt, nPtBins, kMinPt, kMaxPt / 2));
    hJetRPtPtd.setObject(new TH3F("hJetRPtPtd", "Jets;#it{R};#it{p}_{T};ptD", nRBins, kMinR, kMaxR, nPtBins, kMinPt, kMaxPt, nPtBins / 2, 0., 1.));
    hJetRPtChargeFrag.setObject(new TH3F("hJetRPtChargeFrag", "Jets;#it{R};#it{p}_{T};z", nRBins, kMinR, kMaxR, nPtBins, kMinPt, kMaxPt, nPtBins / 2, 0., 1.));
    hJetRPtNEF.setObject(new TH3F("hJetRPtNEF", "Jets;#it{R};#it{p}_{T};NEF", nRBins, kMinR, kMaxR, nPtBins, kMinPt, kMaxPt, nPtBins / 2, 0., 1.));
    hJetRPtZTheta.setObject(new TH3F("hJetRPtZTheta", "Jets;#it{R};#it{p}_{T};z#theta", nRBins, kMinR, kMaxR, nPtBins, kMinPt, kMaxPt, nPtBins / 2, 0., 1.));
    hJetRPtZSqTheta.setObject(new TH3F("hJetRPtZSqTheta", "Jets;#it{R};#it{p}_{T};z^{2} #theta", nRBins, kMinR, kMaxR, nPtBins, kMinPt, kMaxPt, nPtBins / 2, 0., 1.));
    hJetRPtZThetaSq.setObject(new TH3F("hJetRPtZThetaSq", "Jets;#it{R};#it{p}_{T};z #theta^{2}", nRBins, kMinR, kMaxR, nPtBins, kMinPt, kMaxPt, nPtBins / 2, 0., 1.));
    hJetRMaxPtClusterMaxPt.setObject(new TH3F("hJetRMaxPtClusterMaxPt", "Leading jets and clusters;#it{R};#it{p}_{T};#it{p}_{T}^{clus}", nRBins, kMinR, kMaxR, nPtBins, kMinPt, kMaxPt, nPtBins, kMinPt, kMaxPt / 2));

    // Histograms for triggered events
    hSelectedJetRPtEta.setObject(new TH3F("hSelectedJetRPtEta", "Selected jets #it{p}_{T} and #eta;#it{R};#it{p}_{T};#eta", nRBins, kMinR, kMaxR, nPtBins, kMinPt, kMaxPt, nEtaBins, kMinEta, kMaxEta));
    hSelectedJetRPtPhi.setObject(new TH3F("hSelectedJetRPtPhi", "Selected jets #it{p}_{T} and #phi;#it{R};#it{p}_{T};#phi", nRBins, kMinR, kMaxR, nPtBins, kMinPt, kMaxPt, nPhiBins, kMinPhi, kMaxPhi));
    hSelectedJetRMaxPtEta.setObject(new TH3F("hSelectedJetRMaxPtEta", "Leading selected jets #it{p}_{T} and #eta;#it{R};#it{p}_{T};#eta", nRBins, kMinR, kMaxR, nPtBins, kMinPt, kMaxPt, nEtaBins, kMinEta, kMaxEta));
    hSelectedJetRMaxPtPhi.setObject(new TH3F("hSelectedJetRMaxPtPhi", "Leading selected jets #it{p}_{T} and #phi;#it{R};#it{p}_{T};#phi", nRBins, kMinR, kMaxR, nPtBins, kMinPt, kMaxPt, nPhiBins, kMinPhi, kMaxPhi));
    hSelectedClusterPtEta.setObject(new TH2F("hSelectedClusterPtEta", Form("Selected Cluster %s and #eta;%s;#eta", observableName.data(), observableAxisTitle.data()), nPtBins, kMinPt, kMaxPt / 2, nEtaBins, kMinEta, kMaxEta));
    hSelectedClusterPtPhi.setObject(new TH2F("hSelectedClusterPtPhi", Form("Selected Cluster %s and #phi;%s;#phi", observableName.data(), observableAxisTitle.data()), nPtBins, kMinPt, kMaxPt / 2, nPhiBins, kMinPhi, kMaxPhi));
    hSelectedClusterPtEtaPhi.setObject(new TH3F("hSelectedClusterPtEtaPhi", Form("Selected cluster %s, #eta and #phi;%s;#eta;#phi", observableName.data(), observableAxisTitle.data()), nPtBins, kMinPt, kMaxPt / 2, nEtaBins, kMinEta, kMaxEta, nPhiBins, kMinPhi, kMaxPhi));
    hSelectedClusterMaxPtEtaPhi.setObject(new TH3F("hSelectedClusterMaxPtEtaPhi", Form("Leading Selected cluster %s, #eta and #phi;%s;#eta;#phi", observableName.data(), observableAxisTitle.data()), nPtBins, kMinPt, kMaxPt / 2, nEtaBins, kMinEta, kMaxEta, nPhiBins, kMinPhi, kMaxPhi));
    hSelectedClusterMaxPtEta.setObject(new TH2F("hSelectedClusterMaxPtEta", Form("Leading selected clusters %s and #eta;%s;#eta", observableName.data(), observableAxisTitle.data()), nPtBins, kMinPt, kMaxPt / 2, nEtaBins, kMinEta, kMaxEta));
    hSelectedClusterMaxPtPhi.setObject(new TH2F("hSelectedClusterMaxPtPhi", Form("Leading selected clusters %s and #phi;%s;#phi", observableName.data(), observableAxisTitle.data()), nPtBins, kMinPt, kMaxPt / 2, nPhiBins, kMinPhi, kMaxPhi));
    hSelectedJetRPtTrackPt.setObject(new TH3F("hSelectedJetRPtTrackPt", "Selected jets;#it{R};#it{p}_{T};#it{p}_{T}^{track}", nRBins, kMinR, kMaxR, nPtBins, kMinPt, kMaxPt, nPtBins, kMinPt, kMaxPt / 2));
    hSelectedJetRPtClusterPt.setObject(new TH3F("hSelectedJetRPtClusterPt", "Selected jets;#it{R};#it{p}_{T};#it{p}_{T}^{clus}", nRBins, kMinR, kMaxR, nPtBins, kMinPt, kMaxPt, nPtBins, kMinPt, kMaxPt / 2));
    hSelectedJetRPtPtd.setObject(new TH3F("hSelectedJetRPtPtd", "Selected jets;#it{R};#it{p}_{T};ptD", nRBins, kMinR, kMaxR, nPtBins, kMinPt, kMaxPt, nPtBins / 2, 0., 1.));
    hSelectedJetRPtChargeFrag.setObject(new TH3F("hSelectedJetRPtChargeFrag", "Selected jets;#it{R};#it{p}_{T};z", nRBins, kMinR, kMaxR, nPtBins, kMinPt, kMaxPt, nPtBins / 2, 0., 1.));
    hSelectedJetRPtNEF.setObject(new TH3F("hSelectedJetRPtNEF", "Selected jets;#it{R};#it{p}_{T};NEF", nRBins, kMinR, kMaxR, nPtBins, kMinPt, kMaxPt, nPtBins / 2, 0., 1.));
    hSelectedJetRPtZTheta.setObject(new TH3F("hSelectedJetRPtZTheta", "Selected jets;#it{R};#it{p}_{T};z#theta", nRBins, kMinR, kMaxR, nPtBins, kMinPt, kMaxPt, nPtBins / 2, 0., 1.));
    hSelectedJetRPtZSqTheta.setObject(new TH3F("hSelectedJetRPtZSqTheta", "Selected jets;#it{R};#it{p}_{T};z^{2} #theta", nRBins, kMinR, kMaxR, nPtBins, kMinPt, kMaxPt, nPtBins / 2, 0., 1.));
    hSelectedJetRPtZThetaSq.setObject(new TH3F("hSelectedJetRPtZThetaSq", "Selected jets;#it{R};#it{p}_{T};z #theta^{2}", nRBins, kMinR, kMaxR, nPtBins, kMinPt, kMaxPt, nPtBins / 2, 0., 1.));
    hSelectedGammaEMCALPtEta.setObject(new TH2F("hSelectedGammaEMCALPtEta", Form("Selected GammaEMCAL %s and #eta;%s;#eta", observableName.data(), observableAxisTitle.data()), nPtBins, kMinPt, kMaxPt / 2, nEtaBins, kMinEta, kMaxEta));
    hSelectedGammaEMCALPtPhi.setObject(new TH2F("hSelectedGammaEMCALPtPhi", Form("Selected GammaEMCAL %s and #phi;%s;#phi", observableName.data(), observableAxisTitle.data()), nPtBins, kMinPt, kMaxPt / 2, nPhiBins, kMinPhi, kMaxPhi));
    hSelectedGammaDCALPtEta.setObject(new TH2F("hSelectedGammaDCALPtEta", Form("Selected GammaDCAL %s and #eta;%s;#eta", observableName.data(), observableAxisTitle.data()), nPtBins, kMinPt, kMaxPt / 2, nEtaBins, kMinEta, kMaxEta));
    hSelectedGammaDCALPtPhi.setObject(new TH2F("hSelectedGammaDCALPtPhi", Form("Selected GammaDCAL %s and #phi;%s;#phi", observableName.data(), observableAxisTitle.data()), nPtBins, kMinPt, kMaxPt / 2, nPhiBins, kMinPhi, kMaxPhi));
    hSelectedGammaEMCALPtEtaPhi.setObject(new TH3F("hSelectedGammaEMCALPtEtaPhi", Form("Selected GammaEMCAL %s, #eta and #phi;%s;#eta;#phi", observableName.data(), observableAxisTitle.data()), nPtBins, kMinPt, kMaxPt / 2, nEtaBins, kMinEta, kMaxEta, nPhiBins, kMinPhi, kMaxPhi));
    hSelectedGammaEMCALMaxPtEtaPhi.setObject(new TH3F("hSelectedGammaEMCALMaxPtEtaPhi", Form("Leading selected gammaEMCALs %s, #eta and #phi;%s;#eta;#phi", observableName.data(), observableAxisTitle.data()), nPtBins, kMinPt, kMaxPt / 2, nEtaBins, kMinEta, kMaxEta, nPhiBins, kMinPhi, kMaxPhi));
    hSelectedGammaEMCALMaxPtEta.setObject(new TH2F("hSelectedGammaEMCALMaxPtEta", Form("Leading selected gammaEMCALs %s and #eta;%s;#eta", observableName.data(), observableAxisTitle.data()), nPtBins, kMinPt, kMaxPt / 2, nEtaBins, kMinEta, kMaxEta));
    hSelectedGammaEMCALMaxPtPhi.setObject(new TH2F("hSelectedGammaEMCALMaxPtPhi", Form("Leading selected gammaEMCALs %s and #phi;%s;#phi", observableName.data(), observableAxisTitle.data()), nPtBins, kMinPt, kMaxPt / 2, nPhiBins, kMinPhi, kMaxPhi));
    hSelectedGammaDCALMaxPtEta.setObject(new TH2F("hSelectedGammaDCALMaxPtEta", Form("Leading selected gammaDCALs %s and #eta;%s;#eta", observableName.data(), observableAxisTitle.data()), nPtBins, kMinPt, kMaxPt / 2, nEtaBins, kMinEta, kMaxEta));
    hSelectedGammaDCALMaxPtPhi.setObject(new TH2F("hSelectedGammaDCALMaxPtPhi", Form("Leading selected gammaDCALs %s and #phi;%s;#phi", observableName.data(), observableAxisTitle.data()), nPtBins, kMinPt, kMaxPt / 2, nPhiBins, kMinPhi, kMaxPhi));
    hSelectedGammaDCALPtEtaPhi.setObject(new TH3F("hSelectedGammaDCALPtEtaPhi", Form("Selected GammaDCAL %s, #eta and #phi;%s;#eta;#phi", observableName.data(), observableAxisTitle.data()), nPtBins, kMinPt, kMaxPt / 2, nEtaBins, kMinEta, kMaxEta, nPhiBins, kMinPhi, kMaxPhi));
    hSelectedGammaDCALMaxPtEtaPhi.setObject(new TH3F("hSelectedGammaDCALMaxPtEtaPhi", Form("Leading selected GammaDCALs %s, #eta and #phi;%s;#eta;#phi", observableName.data(), observableAxisTitle.data()), nPtBins, kMinPt, kMaxPt / 2, nEtaBins, kMinEta, kMaxEta, nPhiBins, kMinPhi, kMaxPhi));
    hSelectedGammaEMCALPtEtaLow.setObject(new TH2F("hSelectedGammaEMCALPtEtaLow", Form("Selected GammaEMCAL %s and #eta (low threshold);%s;#eta", observableName.data(), observableAxisTitle.data()), nPtBins, kMinPt, kMaxPt / 2, nEtaBins, kMinEta, kMaxEta));
    hSelectedGammaEMCALPtPhiLow.setObject(new TH2F("hSelectedGammaEMCALPtPhiLow", Form("Selected GammaEMCAL %s and #phi (low threshold);%s;#phi", observableName.data(), observableAxisTitle.data()), nPtBins, kMinPt, kMaxPt / 2, nPhiBins, kMinPhi, kMaxPhi));
    hSelectedGammaDCALPtEtaLow.setObject(new TH2F("hSelectedGammaDCALPtEtaLow", Form("Selected GammaDCAL %s and #eta (low threshold);%s;#eta", observableName.data(), observableAxisTitle.data()), nPtBins, kMinPt, kMaxPt / 2, nEtaBins, kMinEta, kMaxEta));
    hSelectedGammaDCALPtPhiLow.setObject(new TH2F("hSelectedGammaDCALPtPhiLow", Form("Selected GammaDCAL %s and #phi (low threshold);%s;#phi", observableName.data(), observableAxisTitle.data()), nPtBins, kMinPt, kMaxPt / 2, nPhiBins, kMinPhi, kMaxPhi));
    hSelectedGammaEMCALPtEtaPhiLow.setObject(new TH3F("hSelectedGammaEMCALPtEtaPhiLow", Form("Selected GammaEMCAL %s, #eta and #phi (low threshold);%s;#eta;#phi", observableName.data(), observableAxisTitle.data()), nPtBins, kMinPt, kMaxPt / 2, nEtaBins, kMinEta, kMaxEta, nPhiBins, kMinPhi, kMaxPhi));
    hSelectedGammaEMCALMaxPtEtaPhiLow.setObject(new TH3F("hSelectedGammaEMCALMaxPtEtaPhiLow", Form("Leading selected gammaEMCALs %s, #eta and #phi (low threshold);%s;#eta;#phi", observableName.data(), observableAxisTitle.data()), nPtBins, kMinPt, kMaxPt / 2, nEtaBins, kMinEta, kMaxEta, nPhiBins, kMinPhi, kMaxPhi));
    hSelectedGammaEMCALMaxPtEtaLow.setObject(new TH2F("hSelectedGammaEMCALMaxPtEtaLow", Form("Leading selected gammaEMCALs %s and #eta (low threshold);%s;#eta", observableName.data(), observableAxisTitle.data()), nPtBins, kMinPt, kMaxPt / 2, nEtaBins, kMinEta, kMaxEta));
    hSelectedGammaEMCALMaxPtPhiLow.setObject(new TH2F("hSelectedGammaEMCALMaxPtPhiLow", Form("Leading selected gammaEMCALs %s and #phi (low threshold);%s;#phi", observableName.data(), observableAxisTitle.data()), nPtBins, kMinPt, kMaxPt / 2, nPhiBins, kMinPhi, kMaxPhi));
    hSelectedGammaDCALMaxPtEtaLow.setObject(new TH2F("hSelectedGammaDCALMaxPtEtaLow", Form("Leading selected gammaDCALs %s and #eta (low threshold);%s;#eta", observableName.data(), observableAxisTitle.data()), nPtBins, kMinPt, kMaxPt / 2, nEtaBins, kMinEta, kMaxEta));
    hSelectedGammaDCALMaxPtPhiLow.setObject(new TH2F("hSelectedGammaDCALMaxPtPhiLow", Form("Leading selected gammaDCALs %s and #phi (low threshold);%s;#phi", observableName.data(), observableAxisTitle.data()), nPtBins, kMinPt, kMaxPt / 2, nPhiBins, kMinPhi, kMaxPhi));
    hSelectedGammaDCALPtEtaPhiLow.setObject(new TH3F("hSelectedGammaDCALPtEtaPhiLow", Form("Selected GammaDCAL %s, #eta and #phi (low threshold);%s;#eta;#phi", observableName.data(), observableAxisTitle.data()), nPtBins, kMinPt, kMaxPt / 2, nEtaBins, kMinEta, kMaxEta, nPhiBins, kMinPhi, kMaxPhi));
    hSelectedGammaDCALMaxPtEtaPhiLow.setObject(new TH3F("hSelectedGammaDCALMaxPtEtaPhiLow", Form("Leading selected GammaDCALs %s, #eta and #phi (low threshold);%s;#eta;#phi", observableName.data(), observableAxisTitle.data()), nPtBins, kMinPt, kMaxPt / 2, nEtaBins, kMinEta, kMaxEta, nPhiBins, kMinPhi, kMaxPhi));
    hSelectedJetRMaxPtClusterMaxPt.setObject(new TH3F("hSelectedJetRMaxPtClusterMaxPt", "Leading selected jets and clusters;#it{R};#it{p}_{T};#it{p}_{T}^{clus}", nRBins, kMinR, kMaxR, nPtBins, kMinPt, kMaxPt, nPtBins, kMinPt, kMaxPt / 2));

    hJetRMaxPtJetPt.setObject(new TH3F("hJetRMaxPtJetPt", "Leading jet #it{p}_{T} vs jet #it{p}_{T};#it{R};#it{p}_{T}^{max};#it{p}_{T}", nRBins, kMinR, kMaxR, nPtBins, kMinPt, kMaxPt, nPtBins, kMinPt, kMaxPt));
    hJetRMaxPtJetPtNoFiducial.setObject(new TH3F("hJetRMaxPtJetPtNoFiducial", "Leading jet #it{p}_{T} vs jet #it{p}_{T} (no fiducial cut);#it{R};#it{p}_{T}^{max};#it{p}_{T}", nRBins, kMinR, kMaxR, nPtBins, kMinPt, kMaxPt, nPtBins, kMinPt, kMaxPt));
    hClusterEMCALMaxPtClusterEMCALPt.setObject(new TH2F("hClusterEMCALMaxPtClusterEMCALPt", Form("Leading clusterEMCAL %s vs clusterEMCAL %s;%s;%s", observableName.data(), observableName.data(), observableAxisTitleMax.data(), observableAxisTitle.data()), nPtBins, kMinPt, kMaxPt / 2, nPtBins, kMinPt, kMaxPt / 2));
    hClusterDCALMaxPtClusterDCALPt.setObject(new TH2F("hClusterDCALMaxPtClusterDCALPt", Form("Leading clusterDCAL %s vs clusterDCAL %s;%s;%s", observableName.data(), observableName.data(), observableAxisTitleMax.data(), observableAxisTitle.data()), nPtBins, kMinPt, kMaxPt / 2, nPtBins, kMinPt, kMaxPt / 2));

    if (b_JetsInEmcalOnly) {
      hJetRPtEta->SetTitle("Jets (in emcal only) #it{p}_{T} and #eta;it{R};#it{p}_{T};#eta");
      hJetRPtPhi->SetTitle("Jets (in emcal only) #it{p}_{T} and #phi;it{R};#it{p}_{T};#phi");
      hJetRPtTrackPt->SetTitle("Jets (in emcal only);it{R};#it{p}_{T};#it{p}_{T}^{track}");
      hJetRPtClusterPt->SetTitle("Jets (in emcal only);it{R};#it{p}_{T};#it{p}_{T}^{clus}");
      hJetRPtPtd->SetTitle("Jets (in emcal only);it{R};#it{p}_{T};ptD");
      hJetRPtChargeFrag->SetTitle("Jets (in emcal only);it{R};#it{p}_{T};z");
      hJetRPtNEF->SetTitle("Jets (in emcal only);it{R};#it{p}_{T};NEF");
      hJetRPtZTheta->SetTitle("Jets (in emcal only);it{R};#it{p}_{T};z#theta");
      hJetRPtZSqTheta->SetTitle("Jets (in emcal only);it{R};#it{p}_{T};z^{2} #theta");
      hJetRPtZThetaSq->SetTitle("Jets (in emcal only);it{R};#it{p}_{T};z #theta^{2}");
      hJetRMaxPtEta->SetTitle("Leading jets (in emcal only) #it{p}_{T} and #eta;#it{R};#it{p}_{T};#eta");
      hJetRMaxPtPhi->SetTitle("Leading jets (in emcal only) #it{p}_{T} and #phi;#it{R};#it{p}_{T};#phi");
      hJetRMaxPtClusterMaxPt->SetTitle("Leading jets (in emcal only) and clusters;#it{R};#it{p}_{T};#it{p}_{T}^{clus}");

      hSelectedJetRPtEta->SetTitle("Selected jets (in emcal only) #it{p}_{T} and #eta;it{R};#it{p}_{T};#eta");
      hSelectedJetRPtPhi->SetTitle("Selected jets (in emcal only) #it{p}_{T} and #phi;it{R};#it{p}_{T};#phi");
      hSelectedJetRPtTrackPt->SetTitle("Selected jets (in emcal only);it{R};#it{p}_{T};#it{p}_{T}^{track}");
      hSelectedJetRPtClusterPt->SetTitle("Selected jets (in emcal only);it{R};#it{p}_{T};#it{p}_{T}^{clus}");
      hSelectedJetRPtPtd->SetTitle("Selected jets (in emcal only);it{R};#it{p}_{T};ptD");
      hSelectedJetRPtChargeFrag->SetTitle("Selected jets (in emcal only);it{R};#it{p}_{T};z");
      hSelectedJetRPtNEF->SetTitle("Selected jets (in emcal only);it{R};#it{p}_{T};NEF");
      hSelectedJetRPtZTheta->SetTitle("Selected jets (in emcal only);it{R};#it{p}_{T};z#theta");
      hSelectedJetRPtZSqTheta->SetTitle("Selected jets (in emcal only);it{R};#it{p}_{T};z^{2} #theta");
      hSelectedJetRPtZThetaSq->SetTitle("Selected jets (in emcal only);it{R};#it{p}_{T};z #theta^{2}");
      hSelectedJetRMaxPtEta->SetTitle("Leading selected jets (in emcal only) #it{p}_{T} and #eta;#it{R};#it{p}_{T};#eta");
      hSelectedJetRMaxPtPhi->SetTitle("Leading selected jets (in emcal only) #it{p}_{T} and #phi;#it{R};#it{p}_{T};#phi");
      hSelectedJetRMaxPtClusterMaxPt->SetTitle("Leading selected jets (in emcal only) and clusters;#it{R};#it{p}_{T};#it{p}_{T}^{clus}");

      hJetRMaxPtJetPt.setObject(new TH3F("hJetRMaxPtJetPt", "Leading jet #it{p}_{T} vs jet #it{p}_{T};#it{R};#it{p}_{T}^{max};#it{p}_{T}", nRBins, kMinR, kMaxR, nPtBins, kMinPt, kMaxPt, nPtBins, kMinPt, kMaxPt));
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

  Bool_t isJetInEmcal(aod::Jet const& jet)
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

  void check_maxJetPt(aod::Jet jet, std::vector<aod::Jet>& vecMaxJet)
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

  void process(soa::Join<aod::Collisions, aod::EvSels, aod::FullJetFilters>::iterator const& collision,
               aod::Jets const& jets,
               aod::JetTrackConstituents const& jetTrackConstituents,
               aod::JetClusterConstituents const& jetClusterConstituents,
               aod::Tracks const& tracks,
               selectedClusters const& clusters)
  {
    hProcessedEvents->Fill(0);

    // If ignore is true, we don't check for the flag
    // If ignore is true, we don't fill hProcessedEvents for EMCAL
    if (!b_IgnoreEmcalFlag) {
      if (!collision.alias()[kTVXinEMC]) {
        return; // Only consider events where EMCAL is live
      } else {
        hProcessedEvents->Fill(1);
      }
    }

    bool isEvtSelected = false, isEventSelectedNeutralJet = false, isEvtSelectedGammaEMCAL = false, isEvtSelectedGammaDCAL = false, isEvtSelectedGammaLowEMCAL = false, isEvtSelectedGammaLowDCAL = false;
    if (collision.hasJetFullHighPt()) {
      isEvtSelected = true;
      hProcessedEvents->Fill(2);
    }
    if (collision.hasJetNeutralHighPt()) {
      isEventSelectedNeutralJet = true;
      hProcessedEvents->Fill(3);
    }
    if (collision.hasGammaHighPtEMCAL()) {
      isEvtSelectedGammaEMCAL = true;
      hProcessedEvents->Fill(4);
    }
    if (collision.hasGammaHighPtDCAL()) {
      isEvtSelectedGammaDCAL = true;
      hProcessedEvents->Fill(5);
    }
    if (collision.hasGammaLowPtEMCAL()) {
      isEvtSelectedGammaLowEMCAL = true;
      hProcessedEvents->Fill(6);
    }
    if (collision.hasGammaLowPtDCAL()) {
      isEvtSelectedGammaLowDCAL = true;
      hProcessedEvents->Fill(7);
    }
    if (collision.hasJetFullHighPt() && collision.hasGammaHighPtEMCAL()) {
      hProcessedEvents->Fill(8);
    }
    if (collision.hasJetFullHighPt() && collision.hasGammaHighPtDCAL()) {
      hProcessedEvents->Fill(9);
    }
    if (collision.hasJetNeutralHighPt() && collision.hasGammaHighPtEMCAL()) {
      hProcessedEvents->Fill(10);
    }
    if (collision.hasJetNeutralHighPt() && collision.hasGammaHighPtDCAL()) {
      hProcessedEvents->Fill(11);
    }
    if (collision.hasJetFullHighPt() && collision.hasJetNeutralHighPt()) {
      hProcessedEvents->Fill(12);
    }
    if (collision.hasGammaHighPtEMCAL() && collision.hasGammaHighPtDCAL()) {
      hProcessedEvents->Fill(13);
    }

    double maxClusterObservableEMCAL = -1., maxClusterObservableDCAL = -1.;
    selectedClusters::iterator maxClusterEMCAL, maxClusterDCAL;
    std::vector<aod::Jet> vecMaxJet;
    std::vector<aod::Jet> vecMaxJetNoFiducial;

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
      auto tracksInJet = jetTrackConstituents.sliceBy(perJetTrackConstituents, jet.globalIndex());
      for (const auto& trackList : tracksInJet) {
        auto trackPt = trackList.track().pt();
        auto chargeFrag = trackList.track().px() * jet.px() + trackList.track().py() * jet.py() + trackList.track().pz() * jet.pz();
        chargeFrag /= (jet.p() * jet.p());
        auto z = trackPt / jetPt;
        auto dphi = jet.phi() - trackList.track().phi();
        dphi = check_dphi(dphi);
        auto dR = TMath::Sqrt(dphi * dphi + TMath::Power(jet.eta() - trackList.track().eta(), 2));
        zTheta += z * dR / jetR;
        zSqTheta += z * z * dR / jetR;
        zThetaSq += z * dR * dR / (jetR * jetR);
        ptD += z * z;
        hJetRPtTrackPt->Fill(jetR, jetPt, trackPt);
        hJetRPtChargeFrag->Fill(jetR, jetPt, chargeFrag);
        if (isEvtSelected) {
          hSelectedJetRPtTrackPt->Fill(jetR, jetPt, trackPt);
          hSelectedJetRPtChargeFrag->Fill(jetR, jetPt, chargeFrag);
        }
      } // for tracks in jet

      auto clustersInJet = jetClusterConstituents.sliceBy(perJetClusterConstituents, jet.globalIndex());
      for (const auto& clusterList : clustersInJet) {
        auto clusterPt = clusterList.cluster().energy() / std::cosh(clusterList.cluster().eta());
        neutralEnergyFraction += clusterList.cluster().energy();
        auto z = clusterPt / jetPt;
        auto dphi = jet.phi() - clusterList.cluster().phi();
        dphi = check_dphi(dphi);
        auto dR = TMath::Sqrt(dphi * dphi + TMath::Power(jet.eta() - clusterList.cluster().eta(), 2));
        zTheta += z * dR / jetR;
        zSqTheta += z * z * dR / jetR;
        zThetaSq += z * dR * dR / (jetR * jetR);
        ptD += z * z;
        hJetRPtClusterPt->Fill(jetR, jetPt, clusterPt);
        if (isEvtSelected) {
          hSelectedJetRPtClusterPt->Fill(jetR, jetPt, clusterPt);
        }
      } // for clusters in jet
      neutralEnergyFraction /= jet.energy();
      ptD = TMath::Sqrt(ptD);

      // Fillng histograms
      hJetRPtEta->Fill(jetR, jetPt, jet.eta());
      hJetRPtPhi->Fill(jetR, jetPt, jet.phi());
      hJetRPtPtd->Fill(jetR, jetPt, ptD);
      hJetRPtNEF->Fill(jetR, jetPt, neutralEnergyFraction);
      hJetRPtZTheta->Fill(jetR, jetPt, zTheta);
      hJetRPtZSqTheta->Fill(jetR, jetPt, zSqTheta);
      hJetRPtZThetaSq->Fill(jetR, jetPt, zThetaSq);
      registry.get<THn>(HIST("jetRPtEtaPhi"))->Fill(jetR, jetPt, jet.eta(), jet.phi());
      if (isEvtSelected) {
        hSelectedJetRPtEta->Fill(jetR, jetPt, jet.eta());
        hSelectedJetRPtPhi->Fill(jetR, jetPt, jet.phi());
        hSelectedJetRPtPtd->Fill(jetR, jetPt, ptD);
        hSelectedJetRPtNEF->Fill(jetR, jetPt, neutralEnergyFraction);
        hSelectedJetRPtZTheta->Fill(jetR, jetPt, zTheta);
        hSelectedJetRPtZSqTheta->Fill(jetR, jetPt, zSqTheta);
        hSelectedJetRPtZThetaSq->Fill(jetR, jetPt, zThetaSq);
      }
    } // for jets

    struct ClusterData {
      float mTriggerObservable;
      float mEta;
      float mPhi;
    };
    std::vector<ClusterData> analysedClusters;

    for (const auto& cluster : clusters) {
      if (b_RejectExoticClusters && cluster.isExotic()) {
        continue;
      }
      if (cluster.time() < f_minClusterTime || cluster.time() > f_maxClusterTime) {
        continue;
      }
      double clusterObservable = (f_GammaObservable == 0) ? cluster.energy() : cluster.energy() / std::cosh(cluster.eta());
      analysedClusters.push_back({static_cast<float>(clusterObservable), cluster.eta(), cluster.phi()});

      if (isClusterInEmcal(cluster) && clusterObservable > maxClusterObservableEMCAL) {
        maxClusterObservableEMCAL = clusterObservable;
        maxClusterEMCAL = cluster;
      }
      if (!isClusterInEmcal(cluster) && clusterObservable > maxClusterObservableDCAL) {
        maxClusterObservableDCAL = clusterObservable;
        maxClusterDCAL = cluster;
      }
      hClusterPtEta->Fill(clusterObservable, cluster.eta());
      hClusterPtPhi->Fill(clusterObservable, cluster.phi());
      hClusterPtEtaPhi->Fill(clusterObservable, cluster.eta(), cluster.phi());
      if (isEvtSelected || isEventSelectedNeutralJet) {
        hSelectedClusterPtEta->Fill(clusterObservable, cluster.eta());
        hSelectedClusterPtPhi->Fill(clusterObservable, cluster.phi());
        hSelectedClusterPtEtaPhi->Fill(clusterObservable, cluster.eta(), cluster.phi());
      }
      if (isEvtSelectedGammaEMCAL && isClusterInEmcal(cluster)) { // Only fill EMCAL clusters
        hSelectedGammaEMCALPtEta->Fill(clusterObservable, cluster.eta());
        hSelectedGammaEMCALPtPhi->Fill(clusterObservable, cluster.phi());
        hSelectedGammaEMCALPtEtaPhi->Fill(clusterObservable, cluster.eta(), cluster.phi());
      }
      if (isEvtSelectedGammaDCAL && !isClusterInEmcal(cluster)) { // Only fill DCAL clusters
        hSelectedGammaDCALPtEta->Fill(clusterObservable, cluster.eta());
        hSelectedGammaDCALPtPhi->Fill(clusterObservable, cluster.phi());
        hSelectedGammaDCALPtEtaPhi->Fill(clusterObservable, cluster.eta(), cluster.phi());
      }
      if (isEvtSelectedGammaLowEMCAL && isClusterInEmcal(cluster)) { // Only fill EMCAL clusters
        hSelectedGammaEMCALPtEtaLow->Fill(clusterObservable, cluster.eta());
        hSelectedGammaEMCALPtPhiLow->Fill(clusterObservable, cluster.phi());
        hSelectedGammaEMCALPtEtaPhiLow->Fill(clusterObservable, cluster.eta(), cluster.phi());
      }
      if (isEvtSelectedGammaLowDCAL && !isClusterInEmcal(cluster)) { // Only fill DCAL clusters
        hSelectedGammaDCALPtEtaLow->Fill(clusterObservable, cluster.eta());
        hSelectedGammaDCALPtPhiLow->Fill(clusterObservable, cluster.phi());
        hSelectedGammaDCALPtEtaPhiLow->Fill(clusterObservable, cluster.eta(), cluster.phi());
      }
    } // for clusters

    if (maxClusterObservableEMCAL > 0) {
      // hClusterMaxPtEta->Fill(maxClusterPt, maxCluster.eta());
      // hClusterMaxPtPhi->Fill(maxClusterPt, maxCluster.phi());
      hClusterEMCALMaxPtEtaPhi->Fill(maxClusterObservableEMCAL, maxClusterEMCAL.eta(), maxClusterEMCAL.phi());
      for (const auto& cluster : analysedClusters) {
        hClusterEMCALMaxPtClusterEMCALPt->Fill(maxClusterObservableEMCAL, cluster.mTriggerObservable);
      }
      if (isEvtSelected) {
        hSelectedClusterMaxPtEta->Fill(maxClusterObservableEMCAL, maxClusterEMCAL.eta());
        hSelectedClusterMaxPtPhi->Fill(maxClusterObservableEMCAL, maxClusterEMCAL.phi());
        hSelectedClusterMaxPtEtaPhi->Fill(maxClusterObservableEMCAL, maxClusterEMCAL.eta(), maxClusterEMCAL.phi());
      }
      if (isEvtSelectedGammaEMCAL) {
        hSelectedGammaEMCALMaxPtEta->Fill(maxClusterObservableEMCAL, maxClusterEMCAL.eta());
        hSelectedGammaEMCALMaxPtPhi->Fill(maxClusterObservableEMCAL, maxClusterEMCAL.phi());
        hSelectedGammaEMCALMaxPtEtaPhi->Fill(maxClusterObservableEMCAL, maxClusterEMCAL.eta(), maxClusterEMCAL.phi());
      }
      if (isEvtSelectedGammaLowEMCAL) {
        hSelectedGammaEMCALMaxPtEtaLow->Fill(maxClusterObservableEMCAL, maxClusterEMCAL.eta());
        hSelectedGammaEMCALMaxPtPhiLow->Fill(maxClusterObservableEMCAL, maxClusterEMCAL.phi());
        hSelectedGammaEMCALMaxPtEtaPhiLow->Fill(maxClusterObservableEMCAL, maxClusterEMCAL.eta(), maxClusterEMCAL.phi());
      }
    }
    if (maxClusterObservableDCAL > 0) {
      hClusterDCALMaxPtEtaPhi->Fill(maxClusterObservableDCAL, maxClusterDCAL.eta(), maxClusterDCAL.phi());
      for (const auto& cluster : analysedClusters) {
        hClusterDCALMaxPtClusterDCALPt->Fill(maxClusterObservableDCAL, cluster.mTriggerObservable);
      }
      if (isEvtSelectedGammaDCAL) {
        hSelectedGammaDCALMaxPtEta->Fill(maxClusterObservableDCAL, maxClusterDCAL.eta());
        hSelectedGammaDCALMaxPtPhi->Fill(maxClusterObservableDCAL, maxClusterDCAL.phi());
        hSelectedGammaDCALMaxPtEtaPhi->Fill(maxClusterObservableDCAL, maxClusterDCAL.eta(), maxClusterDCAL.phi());
      }
      if (isEvtSelectedGammaLowDCAL) {
        hSelectedGammaDCALMaxPtEtaLow->Fill(maxClusterObservableDCAL, maxClusterDCAL.eta());
        hSelectedGammaDCALMaxPtPhiLow->Fill(maxClusterObservableDCAL, maxClusterDCAL.phi());
        hSelectedGammaDCALMaxPtEtaPhiLow->Fill(maxClusterObservableDCAL, maxClusterDCAL.eta(), maxClusterDCAL.phi());
      }
    }

    for (auto maxJet : vecMaxJet) {
      double jetR = maxJet.r() * 1e-2, jetPt = maxJet.pt(), jetEta = maxJet.eta(), jetPhi = maxJet.phi();
      hJetRMaxPtEta->Fill(jetR, jetPt, jetEta);
      hJetRMaxPtPhi->Fill(jetR, jetPt, jetPhi);
      // hJetRMaxPtEtaPhi->Fill(jetR, jetPt, jetEta, jetPhi);
      registry.get<THn>(HIST("jetRMaxPtEtaPhi"))->Fill(jetR, jetPt, jetEta, jetPhi);
      if (isEvtSelected) {
        hSelectedJetRMaxPtEta->Fill(jetR, jetPt, jetEta);
        hSelectedJetRMaxPtPhi->Fill(jetR, jetPt, jetPhi);
      }
      if (maxClusterObservableEMCAL > 0) {
        hJetRMaxPtClusterMaxPt->Fill(jetR, jetPt, maxClusterObservableEMCAL);
        if (isEvtSelected) {
          hSelectedJetRMaxPtClusterMaxPt->Fill(jetR, jetPt, maxClusterObservableEMCAL);
        }
      } // if maxClusterPt
      if (maxJet.r() == std::round(f_jetR * 100)) {
        for (const auto& jet : jets) {
          if (isJetInEmcal(jet)) {
            hJetRMaxPtJetPt->Fill(jet.r() * 1e-2, jetPt, jet.pt());
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
          hJetRMaxPtJetPtNoFiducial->Fill(jet.r() * 1e-2, jetPt, jet.pt());
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
