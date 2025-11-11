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
/// \author Gijs van Weelden <g.van.weelden@cern.ch>
//
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/DataModel/EMCALClusterDefinition.h"
#include "PWGJE/DataModel/EMCALClusters.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"

#include "Common/CCDB/TriggerAliases.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisTask.h"
#include <CommonConstants/MathConstants.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include "TTree.h"
#include <TH2.h>
#include <TH3.h>
#include <THn.h>
#include <TString.h>

#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>

#include <algorithm>
#include <array>
#include <bitset>
#include <cmath>
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct JetTriggerQA {
  using selectedClusters = o2::soa::Filtered<aod::JetClusters>;
  using fullJetInfos = soa::Join<aod::FullJets, aod::FullJetConstituents>;
  using neutralJetInfos = soa::Join<aod::NeutralJets, aod::NeutralJetConstituents>;
  using collisionWithTrigger = soa::Join<aod::JCollisions, aod::JFullTrigSels>::iterator;

  enum TriggerType_t {
    kMinBias,
    kEmcalAny,
    kEmcalMB,
    kEmcalJetFull,
    kEmcalJetFullLow,
    kEmcalJetNeutral,
    kEmcalJetNeutralLow,
    kEmcalGammaVeryHigh,
    kDcalGammaVeryHigh,
    kEmcalGammaHigh,
    kDcalGammaHigh,
    kEmcalGammaLow,
    kDcalGammaLow,
    kEmcalGammaVeryLow,
    kDcalGammaVeryLow,
    kNTriggers
  };

  enum EMCALHardwareTrigger {
    TRG_MB,
    TRG_EMC7,
    TRG_DMC7,
    TRG_NTriggers
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
  Configurable<bool> b_doLightOutput{"b_doLightOutput", true, "do light output (disables most histograms not needed for QA)"};

  Configurable<float> cfgVertexCut{"cfgVertexCut", 10.0f, "Accepted z-vertex range"};
  Configurable<std::string> mClusterDefinition{"clusterDefinition", "kV3Default", "cluster definition to be selected, e.g. V3Default"};

  JetFinder jetReclusterer;

  void init(InitContext&)
  {
    jetReclusterer.isReclustering = true;
    jetReclusterer.algorithm = fastjet::JetAlgorithm::cambridge_algorithm;
    jetReclusterer.jetR = f_jetR * 2.5; // Use larger R for CA reclustering to prevent losing particles

    int nPtBins = 200;
    float kMinPt = 0.;
    float kMaxPt = 201.;
    int nPhiBins = 18 * 8;
    float kMinPhi = 0.;
    float kMaxPhi = o2::constants::math::TwoPI;
    int nEtaBins = 100;
    float kMinEta = -1.;
    float kMaxEta = 1.;
    int nRBins = 6;
    float kMinR = 0.05;
    float kMaxR = 0.65;

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
    if (!b_doLightOutput)
      registry.add("jetRMaxPtEtaPhi", "JetRMaxPtEtaPhi", hJetRMaxPtEtaPhi);
    registry.add("jetRPtEtaPhiNoFiducial", "JetRPtEtaPhiNoFiducial", hJetRPtEtaPhiNoFiducial);
    if (!b_doLightOutput)
      registry.add("jetRMaxPtEtaPhiNoFiducial", "JetRMaxPtEtaPhiNoFiducial", hJetRMaxPtEtaPhiNoFiducial);

    registry.add("hProcessedEvents", "Processed events", HistType::kTH1D, {{15, -0.5, 14.5, "Trigger type"}});
    auto histProcessed = registry.get<TH1>(HIST("hProcessedEvents"));
    histProcessed->GetXaxis()->SetBinLabel(1, "MB");
    histProcessed->GetXaxis()->SetBinLabel(2, "Any EMC trigger");
    histProcessed->GetXaxis()->SetBinLabel(3, "EMC Min. bias");
    histProcessed->GetXaxis()->SetBinLabel(4, "Selected Full Jet high");
    histProcessed->GetXaxis()->SetBinLabel(5, "Selected Full Jet low");
    histProcessed->GetXaxis()->SetBinLabel(6, "Selected Neutral Jet high");
    histProcessed->GetXaxis()->SetBinLabel(7, "Selected Neutral Jet low");
    histProcessed->GetXaxis()->SetBinLabel(8, "Selected Gamma very high EMCAL");
    histProcessed->GetXaxis()->SetBinLabel(9, "Selected Gamma very high DCAL");
    histProcessed->GetXaxis()->SetBinLabel(10, "Selected Gamma high EMCAL");
    histProcessed->GetXaxis()->SetBinLabel(11, "Selected Gamma high DCAL");
    histProcessed->GetXaxis()->SetBinLabel(12, "Selected Gamma low EMCAL");
    histProcessed->GetXaxis()->SetBinLabel(13, "Selected Gamma low DCAL");
    histProcessed->GetXaxis()->SetBinLabel(14, "Selected Gamma very low EMCAL");
    histProcessed->GetXaxis()->SetBinLabel(15, "Selected Gamma very low DCAL");

    std::array<std::string, TriggerType_t::kNTriggers> triggerlabels = {{"MB", "EMC Any", "EMC MB", "EMC jet full high", "EMC jet full low", "EMC jet neutral high", "EMC jet neutral low", "EMC gamma very high", "DCL gamma very high", "EMC gamma high", "DCL gamma high", "EMC gamma low", "DCL gamma low", "EMC gamma very low", "DCL gamma very low"}};
    registry.add("hTriggerCorrelation", "Correlation between EMCAL triggers", HistType::kTH2D, {{TriggerType_t::kNTriggers, -0.5, static_cast<double>(TriggerType_t::kNTriggers) - 0.5, "Main trigger"}, {TriggerType_t::kNTriggers, -0.5, static_cast<double>(TriggerType_t::kNTriggers) - 0.5, "Associated trigger"}});
    auto triggerCorrelation = registry.get<TH2>(HIST("hTriggerCorrelation"));
    for (std::size_t triggertype = 0; triggertype < TriggerType_t::kNTriggers; triggertype++) {
      triggerCorrelation->GetXaxis()->SetBinLabel(triggertype + 1, triggerlabels[triggertype].data());
      triggerCorrelation->GetYaxis()->SetBinLabel(triggertype + 1, triggerlabels[triggertype].data());
    }

    // Histograms for events where the EMCAL is live
    registry.add("hJetRPtEta", "Jets #it{p}_{T} and #eta", HistType::kTH3F, {rAxis, jetPtAxis, etaAxis});
    registry.add("hJetRPtPhi", "Jets #it{p}_{T} and #phi", HistType::kTH3F, {rAxis, jetPtAxis, phiAxis});
    if (!b_doLightOutput) {
      registry.add("hJetRMaxPtEta", "Leading jets #it{p}_{T} and #eta", HistType::kTH3F, {rAxis, jetPtAxis, etaAxis});
      registry.add("hJetRMaxPtPhi", "Leading jets #it{p}_{T} and #phi", HistType::kTH3F, {rAxis, jetPtAxis, phiAxis});
      registry.add("hJetRMaxPtEtaMinBias", "Leading jets #it{p}_{T} and #eta (min. bias)", HistType::kTH3F, {rAxis, jetPtAxis, etaAxis});
      registry.add("hJetRMaxPtPhiMinBias", "Leading jets #it{p}_{T} and #phi (min. bias)", HistType::kTH3F, {rAxis, jetPtAxis, phiAxis});
      registry.add("hJetRMaxPtEtaLevel0", "Leading jets #it{p}_{T} and #eta (level0)", HistType::kTH3F, {rAxis, jetPtAxis, etaAxis});
      registry.add("hJetRMaxPtPhiLevel0", "Leading jets #it{p}_{T} and #phi (level0)", HistType::kTH3F, {rAxis, jetPtAxis, phiAxis});
    }
    registry.add("hClusterPtEtaPhi", Form("Cluster %s, #eta and #phi", observableName.data()), HistType::kTH3F, {observableAxisCluster, etaAxis, phiAxis});
    registry.add("hClusterPtEtaPhiMinBias", Form("Cluster %s (Min. bias trigger), #eta and #phi", observableName.data()), HistType::kTH3F, {observableAxisCluster, etaAxis, phiAxis});
    registry.add("hClusterPtEtaPhiLevel0", Form("Cluster %s (Level-0 trigger), #eta and #phi", observableName.data()), HistType::kTH3F, {observableAxisCluster, etaAxis, phiAxis});

    if (!b_doLightOutput) {
      registry.add("hClusterEMCALMaxPtEtaPhi", Form("Leading cluster %s, #eta and #phi (EMCAL)", observableName.data()), HistType::kTH3F, {observableAxisCluster, etaAxis, phiAxis});
      registry.add("hClusterEMCALMaxPtEtaPhiMinBias", Form("Leading cluster %s, #eta and #phi (EMCAL, min. bias trigger)", observableName.data()), HistType::kTH3F, {observableAxisCluster, etaAxis, phiAxis});
      registry.add("hClusterEMCALMaxPtEtaPhiLevel0", Form("Leading cluster %s, #eta and #phi (EMCAL, Level-0 trigger)", observableName.data()), HistType::kTH3F, {observableAxisCluster, etaAxis, phiAxis});
      registry.add("hClusterDCALMaxPtEtaPhi", Form("Leading cluster %s, #eta and #phi (DCAL)", observableName.data()), HistType::kTH3F, {observableAxisCluster, etaAxis, phiAxis});
      registry.add("hClusterDCALMaxPtEtaPhiMinBias", Form("Leading cluster %s, #eta and #phi (DCAL, min, bias trigger)", observableName.data()), HistType::kTH3F, {observableAxisCluster, etaAxis, phiAxis});
      registry.add("hClusterDCALMaxPtEtaPhiLevel0", Form("Leading cluster %s, #eta and #phi (DCAL, Level-0 trigger)", observableName.data()), HistType::kTH3F, {observableAxisCluster, etaAxis, phiAxis});
      registry.add("hEventsNoMaxClusterEMCAL", "Events with no max. cluster in EMCAL", HistType::kTH1D, {{1, 0.5, 1.5}});
      registry.add("hEventsNoMaxClusterEMCALMinBias", "Events with no max. cluster in EMCAL (min. bias trigger)", HistType::kTH1D, {{1, 0.5, 1.5}});
      registry.add("hEventsNoMaxClusterEMCALLevel0", "Events with no max. cluster in EMCAL (level0 trigger)", HistType::kTH1D, {{1, 0.5, 1.5}});
      registry.add("hEventsNoMaxClusterDCAL", ("Events with no max. cluster in DCAL"), HistType::kTH1D, {{1, 0.5, 1.5}});
      registry.add("hEventsNoMaxClusterDCALMinBias", "Events with no max. cluster in DCAL (min. bias trigger)", HistType::kTH1D, {{1, 0.5, 1.5}});
      registry.add("hEventsNoMaxClusterDCALLevel0", "Events with no max. cluster in DCAL (level0 trigger)", HistType::kTH1D, {{1, 0.5, 1.5}});
      registry.add("hEventsNoMaxJet", "Events with no max. jet", HistType::kTH1D, {{rAxis}});
      registry.add("hEventsNoMaxJetMinBias", "Events with no max. jet (min. bias trigger)", HistType::kTH1D, {{rAxis}});
      registry.add("hEventsNoMaxJetLevel0", "Events with no max. jet (level0 trigger)", HistType::kTH1D, {{rAxis}});

      registry.add("hJetRPtTrackPt", "Jets", HistType::kTH3F, {rAxis, jetPtAxis, ptAxisTrackInJet});
      registry.add("hJetRPtClusterPt", "Jets", HistType::kTH3F, {rAxis, jetPtAxis, ptAxisClusterInJet});
      registry.add("hJetRPtPtd", "Jets", HistType::kTH3F, {rAxis, jetPtAxis, {nPtBins / 2, 0., 1., "p_{t,D}"}});
      registry.add("hJetRPtChargeFrag", "Jets", HistType::kTH3F, {rAxis, jetPtAxis, {nPtBins / 2, 0., 1., "z"}});
      registry.add("hJetRPtNEF", "Jets", HistType::kTH3F, {rAxis, jetPtAxis, {nPtBins / 2, 0., 1., "NEF"}});
      registry.add("hJetRPtZTheta", "Jets", HistType::kTH3F, {rAxis, jetPtAxis, {nPtBins / 2, 0., 1., "z#theta"}});
      registry.add("hJetRPtZSqTheta", "Jets", HistType::kTH3F, {rAxis, jetPtAxis, {nPtBins / 2, 0., 1., "z^{2} #theta"}});
      registry.add("hJetRPtZThetaSq", "Jets", HistType::kTH3F, {rAxis, jetPtAxis, {nPtBins / 2, 0., 1., "z #theta^{2}"}});

      registry.add("hJetRMaxPtClusterMaxPt", "Leading jets and clusters", HistType::kTH3F, {rAxis, jetPtAxis, observableAxisCluster});
    }

    // Histograms for triggered events
    registry.add("hSelectedClusterPtEtaPhi", Form("Selected cluster %s, #eta and #phi", observableName.data()), HistType::kTH3F, {observableAxisCluster, etaAxis, phiAxis});

    if (!b_doLightOutput) {
      registry.add("hSelectedClusterMaxPtEtaPhi", Form("Leading Selected cluster %s, #eta and #phi", observableName.data()), HistType::kTH3F, {observableAxisCluster, etaAxis, phiAxis});
    }

    // Jet high trigger
    if (!b_doLightOutput) {
      registry.add("hSelectedJetRMaxPtEta", "Leading selected jets #it{p}_{T} and #eta", HistType::kTH3F, {rAxis, jetPtAxis, etaAxis});
      registry.add("hSelectedJetRMaxPtPhi", "Leading selected jets #it{p}_{T} and #phi", HistType::kTH3F, {rAxis, jetPtAxis, phiAxis});
    }
    registry.add("hSelectedJetRPtEta", "Selected jets #it{p}_{T} and #eta", HistType::kTH3F, {rAxis, jetPtAxis, etaAxis});
    registry.add("hSelectedJetRPtPhi", "Selected jets #it{p}_{T} and #phi", HistType::kTH3F, {rAxis, jetPtAxis, phiAxis});
    if (!b_doLightOutput) {
      registry.add("hSelectedJetRPtTrackPt", "Selected jets", HistType::kTH3F, {rAxis, jetPtAxis, ptAxisTrackInJet});
      registry.add("hSelectedJetRPtClusterPt", "Selected jets", HistType::kTH3F, {rAxis, jetPtAxis, ptAxisClusterInJet});
      registry.add("hSelectedJetRPtPtd", "Selected jets", HistType::kTH3F, {rAxis, jetPtAxis, {nPtBins / 2, 0., 1., "p_{t,D}"}});
      registry.add("hSelectedJetRPtChargeFrag", "Selected jets", HistType::kTH3F, {rAxis, jetPtAxis, {nPtBins / 2, 0., 1., "z"}});
      registry.add("hSelectedJetRPtNEF", "Selected jets", HistType::kTH3F, {rAxis, jetPtAxis, {nPtBins / 2, 0., 1., "NEF"}});
      registry.add("hSelectedJetRPtZTheta", "Selected jets", HistType::kTH3F, {rAxis, jetPtAxis, {nPtBins / 2, 0., 1., "z#theta"}});
      registry.add("hSelectedJetRPtZSqTheta", "Selected jets", HistType::kTH3F, {rAxis, jetPtAxis, {nPtBins / 2, 0., 1., "z^{2} #theta"}});
      registry.add("hSelectedJetRPtZThetaSq", "Selected jets", HistType::kTH3F, {rAxis, jetPtAxis, {nPtBins / 2, 0., 1., "z #theta^{2}"}});
    }

    // Jet low trigger
    if (!b_doLightOutput) {
      registry.add("hSelectedJetLowRMaxPtEta", "Leading selected jets (low threshold) #it{p}_{T} and #eta", HistType::kTH3F, {rAxis, jetPtAxis, etaAxis});
      registry.add("hSelectedJetLowRMaxPtPhi", "Leading selected jets (low threshold) #it{p}_{T} and #phi", HistType::kTH3F, {rAxis, jetPtAxis, phiAxis});
    }
    registry.add("hSelectedJetLowRPtEta", "Selected jets (low threshold) #it{p}_{T} and #eta", HistType::kTH3F, {rAxis, jetPtAxis, etaAxis});
    registry.add("hSelectedJetLowRPtPhi", "Selected jets (low threshold) #it{p}_{T} and #phi", HistType::kTH3F, {rAxis, jetPtAxis, phiAxis});

    if (!b_doLightOutput) {
      registry.add("hSelectedJetLowRPtTrackPt", "Selected jets (low threshold)", HistType::kTH3F, {rAxis, jetPtAxis, ptAxisTrackInJet});
      registry.add("hSelectedJetLowRPtClusterPt", "Selected jets (low threshold)", HistType::kTH3F, {rAxis, jetPtAxis, ptAxisClusterInJet});
      registry.add("hSelectedJetLowRPtPtd", "Selected jets (low threshold)", HistType::kTH3F, {rAxis, jetPtAxis, {nPtBins / 2, 0., 1., "p_{t,D}"}});
      registry.add("hSelectedJetLowRPtChargeFrag", "Selected jets (low threshold)", HistType::kTH3F, {rAxis, jetPtAxis, {nPtBins / 2, 0., 1., "z"}});
      registry.add("hSelectedJetLowRPtNEF", "Selected jets (low threshold)", HistType::kTH3F, {rAxis, jetPtAxis, {nPtBins / 2, 0., 1., "NEF"}});
      registry.add("hSelectedJetLowRPtZTheta", "Selected jets (low threshold)", HistType::kTH3F, {rAxis, jetPtAxis, {nPtBins / 2, 0., 1., "z#theta"}});
      registry.add("hSelectedJetLowRPtZSqTheta", "Selected jets (low threshold)", HistType::kTH3F, {rAxis, jetPtAxis, {nPtBins / 2, 0., 1., "z^{2} #theta"}});
      registry.add("hSelectedJetLowRPtZThetaSq", "Selected jets (low threshold)", HistType::kTH3F, {rAxis, jetPtAxis, {nPtBins / 2, 0., 1., "z #theta^{2}"}});
    }

    // EMCAL gamma very-high trigger
    registry.add("hSelectedGammaEMCALPtEtaPhiVeryHigh", Form("Selected Gamma %s, #eta and #phi (EMCAL, very high threshold)", observableName.data()), HistType::kTH3F, {observableAxisCluster, etaAxis, phiAxis});
    if (!b_doLightOutput)
      registry.add("hSelectedGammaEMCALMaxPtEtaPhiVeryHigh", Form("Leading selected gamma %s, #eta and #phi (EMCAL, very high treshold)", observableName.data()), HistType::kTH3F, {observableAxisCluster, etaAxis, phiAxis});

    // DCAL gamma very-high trigger
    registry.add("hSelectedGammaDCALPtEtaPhiVeryHigh", Form("Selected gamma %s, #eta and #phi (DCAL, very high treshold)", observableName.data()), HistType::kTH3F, {observableAxisCluster, etaAxis, phiAxis});
    if (!b_doLightOutput)
      registry.add("hSelectedGammaDCALMaxPtEtaPhiVeryHigh", Form("Leading selected gamma %s, #eta and #phi (DCAL, very high treshold)", observableName.data()), HistType::kTH3F, {observableAxisCluster, etaAxis, phiAxis});

    // EG1 trigger
    registry.add("hSelectedGammaEMCALPtEtaPhi", Form("Selected Gamma %s, #eta and #phi (EMCAL, high threshold)", observableName.data()), HistType::kTH3F, {observableAxisCluster, etaAxis, phiAxis});
    if (!b_doLightOutput)
      registry.add("hSelectedGammaEMCALMaxPtEtaPhi", Form("Leading selected gamma %s, #eta and #phi (EMCAL, high treshold)", observableName.data()), HistType::kTH3F, {observableAxisCluster, etaAxis, phiAxis});

    // DG1 trigger
    registry.add("hSelectedGammaDCALPtEtaPhi", Form("Selected gamma %s, #eta and #phi (DCAL, high treshold)", observableName.data()), HistType::kTH3F, {observableAxisCluster, etaAxis, phiAxis});
    if (!b_doLightOutput)
      registry.add("hSelectedGammaDCALMaxPtEtaPhi", Form("Leading selected gamma %s, #eta and #phi (DCAL, high treshold)", observableName.data()), HistType::kTH3F, {observableAxisCluster, etaAxis, phiAxis});

    // EG2 trigger
    registry.add("hSelectedGammaEMCALPtEtaPhiLow", Form("Selected gamma %s, #eta and #phi (EMCAL, low threshold)", observableName.data()), HistType::kTH3F, {observableAxisCluster, etaAxis, phiAxis});
    if (!b_doLightOutput)
      registry.add("hSelectedGammaEMCALMaxPtEtaPhiLow", Form("Leading selected gamma %s, #eta and #phi (EMCAL, low threshold)", observableName.data()), HistType::kTH3F, {observableAxisCluster, etaAxis, phiAxis});

    // DG2 trigger
    registry.add("hSelectedGammaDCALPtEtaPhiLow", Form("Selected gamma %s, #eta and #phi (DCAL, low threshold)", observableName.data()), HistType::kTH3F, {observableAxisCluster, etaAxis, phiAxis});
    if (!b_doLightOutput)
      registry.add("hSelectedGammaDCALMaxPtEtaPhiLow", Form("Leading selected gamma %s, #eta and #phi (DCAL, low threshold)", observableName.data()), HistType::kTH3F, {observableAxisCluster, etaAxis, phiAxis});

    // EMCAL gamma very-low trigger
    registry.add("hSelectedGammaEMCALPtEtaPhiVeryLow", Form("Selected gamma %s, #eta and #phi (EMCAL, very low threshold)", observableName.data()), HistType::kTH3F, {observableAxisCluster, etaAxis, phiAxis});
    if (!b_doLightOutput)
      registry.add("hSelectedGammaEMCALMaxPtEtaPhiVeryLow", Form("Leading selected gamma %s, #eta and #phi (EMCAL, very low threshold)", observableName.data()), HistType::kTH3F, {observableAxisCluster, etaAxis, phiAxis});

    //  DCAL gamma very-low trigger
    registry.add("hSelectedGammaDCALPtEtaPhiVeryLow", Form("Selected gamma %s, #eta and #phi (DCAL, low threshold)", observableName.data()), HistType::kTH3F, {observableAxisCluster, etaAxis, phiAxis});
    if (!b_doLightOutput)
      registry.add("hSelectedGammaDCALMaxPtEtaPhiVeryLow", Form("Leading selected gamma %s, #eta and #phi (DCAL, low threshold)", observableName.data()), HistType::kTH3F, {observableAxisCluster, etaAxis, phiAxis});

    if (!b_doLightOutput) {
      registry.add("hSelectedJetRMaxPtClusterMaxPt", "Leading selected jets and clusters", HistType::kTH3F, {rAxis, jetPtAxis, observableAxisCluster});
      registry.add("hJetRMaxPtJetPt", "Leading jet #it{p}_{T} vs jet #it{p}_{T}", HistType::kTH3F, {rAxis, jetMaxPtAxis, jetPtAxis});
      registry.add("hJetRMaxPtJetPtNoFiducial", "Leading jet #it{p}_{T} vs jet #it{p}_{T} (no fiducial cut)", HistType::kTH3F, {rAxis, jetMaxPtAxis, jetPtAxis});
      registry.add("hClusterEMCALMaxPtClusterEMCALPt", Form("Leading clusterEMCAL %s vs clusterEMCAL %s", observableName.data(), observableName.data()), HistType::kTH2F, {observableAxisMaxCluster, observableAxisCluster});
      registry.add("hClusterDCALMaxPtClusterDCALPt", Form("Leading clusterDCAL %s vs clusterDCAL %s", observableName.data(), observableName.data()), HistType::kTH2F, {observableAxisMaxCluster, observableAxisCluster});
    }

    if (b_JetsInEmcalOnly) {
      registry.get<TH3>(HIST("hJetRPtEta"))->SetTitle("Jets (in emcal only) #it{p}_{T} and #eta");
      registry.get<TH3>(HIST("hJetRPtPhi"))->SetTitle("Jets (in emcal only) #it{p}_{T} and #phi");
      if (!b_doLightOutput) {
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
        registry.get<TH3>(HIST("hJetRMaxPtEtaMinBias"))->SetTitle("Leading jets (in emcal only, min. bias) #it{p}_{T} and #eta");
        registry.get<TH3>(HIST("hJetRMaxPtPhiMinBias"))->SetTitle("Leading jets (in emcal only, min. bias) #it{p}_{T} and #phi");
        registry.get<TH3>(HIST("hJetRMaxPtEtaLevel0"))->SetTitle("Leading jets (in emcal only, level-0) #it{p}_{T} and #eta");
        registry.get<TH3>(HIST("hJetRMaxPtPhiLevel0"))->SetTitle("Leading jets (in emcal only, level-0) #it{p}_{T} and #phi");
        registry.get<TH3>(HIST("hJetRMaxPtClusterMaxPt"))->SetTitle("Leading jets (in emcal only) and clusters");
      }

      registry.get<TH3>(HIST("hSelectedJetRPtEta"))->SetTitle("Selected jets (in emcal only) #it{p}_{T} and #eta");
      registry.get<TH3>(HIST("hSelectedJetRPtPhi"))->SetTitle("Selected jets (in emcal only) #it{p}_{T} and #phi");
      if (!b_doLightOutput) {
        registry.get<TH3>(HIST("hSelectedJetRPtTrackPt"))->SetTitle("Selected jets (in emcal only)");
        registry.get<TH3>(HIST("hSelectedJetRPtClusterPt"))->SetTitle("Selected jets (in emcal only)");
        registry.get<TH3>(HIST("hSelectedJetRPtPtd"))->SetTitle("Selected jets (in emcal only)");
        registry.get<TH3>(HIST("hSelectedJetRPtChargeFrag"))->SetTitle("Selected jets (in emcal only)");
        registry.get<TH3>(HIST("hSelectedJetRPtNEF"))->SetTitle("Selected jets (in emcal only)");
        registry.get<TH3>(HIST("hSelectedJetRPtZTheta"))->SetTitle("Selected jets (in emcal only)");
        registry.get<TH3>(HIST("hSelectedJetRPtZSqTheta"))->SetTitle("Selected jets (in emcal only)");
        registry.get<TH3>(HIST("hSelectedJetRPtZThetaSq"))->SetTitle("Selected jets (in emcal only)");
      }

      registry.get<TH3>(HIST("hSelectedJetLowRPtEta"))->SetTitle("Selected jets (low threshold, in emcal only) #it{p}_{T} and #eta");
      registry.get<TH3>(HIST("hSelectedJetLowRPtPhi"))->SetTitle("Selected jets (low threshold, in emcal only) #it{p}_{T} and #phi");
      if (!b_doLightOutput) {
        registry.get<TH3>(HIST("hSelectedJetLowRPtTrackPt"))->SetTitle("Selected jets (low threshold, in emcal only)");
        registry.get<TH3>(HIST("hSelectedJetLowRPtClusterPt"))->SetTitle("Selected jets (low threshold, in emcal only)");
        registry.get<TH3>(HIST("hSelectedJetLowRPtPtd"))->SetTitle("Selected jets (low threshold, in emcal only)");
        registry.get<TH3>(HIST("hSelectedJetLowRPtChargeFrag"))->SetTitle("Selected jets (low threshold, in emcal only)");
        registry.get<TH3>(HIST("hSelectedJetLowRPtNEF"))->SetTitle("Selected jets (low threshold, in emcal only)");
        registry.get<TH3>(HIST("hSelectedJetLowRPtZTheta"))->SetTitle("Selected jets (low threshold, in emcal only)");
        registry.get<TH3>(HIST("hSelectedJetLowRPtZSqTheta"))->SetTitle("Selected jets (low threshold, in emcal only)");
        registry.get<TH3>(HIST("hSelectedJetLowRPtZThetaSq"))->SetTitle("Selected jets (low threshold, in emcal only)");

        registry.get<TH3>(HIST("hSelectedJetRMaxPtEta"))->SetTitle("Leading selected jets (in emcal only) #it{p}_{T} and #eta");
        registry.get<TH3>(HIST("hSelectedJetRMaxPtPhi"))->SetTitle("Leading selected jets (in emcal only) #it{p}_{T} and #phi");
        registry.get<TH3>(HIST("hSelectedJetRMaxPtClusterMaxPt"))->SetTitle("Leading selected jets (in emcal only) and clusters");

        registry.get<TH3>(HIST("hJetRMaxPtJetPt"))->SetTitle("Leading jet #it{p}_{T} vs jet #it{p}_{T} (in emcal only)");
      }
    }
  } // init

  o2::aod::EMCALClusterDefinition clusDef = o2::aod::emcalcluster::getClusterDefinitionFromString(mClusterDefinition.value);
  Filter clusterDefinitionSelection = o2::aod::jcluster::definition == static_cast<int>(clusDef);

  double check_dphi(double dphi) const
  {
    if (dphi > o2::constants::math::PI) {
      dphi -= o2::constants::math::PI;
    } else if (dphi <= -1 * o2::constants::math::PI) {
      dphi += o2::constants::math::PI;
    }
    return dphi;
  }

  template <typename T>
  bool isJetInEmcal(T const& jet)
  {
    double emcalEtaMin = -0.7, emcalEtaMax = 0.7, emcalPhiMin = 1.40, emcalPhiMax = 3.26; // Phi: 80 - 187 deg
    double R = jet.r() * 1e-2;                                                            // Jet R is saved as round(100*R)
    if ((jet.eta() >= emcalEtaMin + R) && (jet.eta() <= emcalEtaMax - R) && (jet.phi() >= emcalPhiMin + R) && (jet.phi() <= emcalPhiMax - R)) {
      return true;
    }
    return false;
  }

  bool hasEMCAL(collisionWithTrigger const& collision) const
  {
    std::array<triggerAliases, 11> selectAliases = {{triggerAliases::kTVXinEMC, triggerAliases::kEMC7, triggerAliases::kDMC7, triggerAliases::kEG1, triggerAliases::kEG2, triggerAliases::kDG1, triggerAliases::kDG2, triggerAliases::kEJ1, triggerAliases::kEJ2, triggerAliases::kDJ1, triggerAliases::kDJ2}};
    bool found = false;
    for (const auto& alias : selectAliases) {
      if (collision.alias_bit(alias)) {
        found = true;
        break;
      }
    }
    return found;
  }

  bool isClusterInEmcal(selectedClusters::iterator const& cluster)
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

  std::pair<double, double> fillGammaQA(const selectedClusters& clusters, const std::bitset<EMCALHardwareTrigger::TRG_NTriggers>& hwtrg, const std::bitset<TriggerType_t::kNTriggers>& triggerstatus)
  {
    auto isTrigger = [&triggerstatus](TriggerType_t triggertype) -> bool {
      return triggerstatus.test(triggertype);
    };

    double maxClusterObservableEMCAL = -1., maxClusterObservableDCAL = -1.;
    selectedClusters::iterator maxClusterEMCAL, maxClusterDCAL;

    struct ClusterData {
      float mTriggerObservable;
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
      analysedClusters.push_back({static_cast<float>(clusterObservable), emcalCluster});

      if (emcalCluster && (clusterObservable > maxClusterObservableEMCAL)) {
        maxClusterObservableEMCAL = clusterObservable;
        maxClusterEMCAL = cluster;
      }
      if (!emcalCluster && (clusterObservable > maxClusterObservableDCAL)) {
        maxClusterObservableDCAL = clusterObservable;
        maxClusterDCAL = cluster;
      }
      registry.fill(HIST("hClusterPtEtaPhi"), clusterObservable, cluster.eta(), cluster.phi());
      if (hwtrg.test(EMCALHardwareTrigger::TRG_MB)) {
        registry.fill(HIST("hClusterPtEtaPhiMinBias"), clusterObservable, cluster.eta(), cluster.phi());
      }
      if (hwtrg.test(EMCALHardwareTrigger::TRG_EMC7) || hwtrg.test(EMCALHardwareTrigger::TRG_DMC7)) {
        registry.fill(HIST("hClusterPtEtaPhiLevel0"), clusterObservable, cluster.eta(), cluster.phi());
      }
      if (isTrigger(TriggerType_t::kEmcalJetFull) || isTrigger(TriggerType_t::kEmcalJetNeutral)) {
        registry.fill(HIST("hSelectedClusterPtEtaPhi"), clusterObservable, cluster.eta(), cluster.phi());
      }
      if (isTrigger(TriggerType_t::kEmcalGammaVeryHigh) && isClusterInEmcal(cluster)) { // Only fill EMCAL clusters
        registry.fill(HIST("hSelectedGammaEMCALPtEtaPhiVeryHigh"), clusterObservable, cluster.eta(), cluster.phi());
      }
      if (isTrigger(TriggerType_t::kDcalGammaVeryHigh) && !isClusterInEmcal(cluster)) { // Only fill DCAL clusters
        registry.fill(HIST("hSelectedGammaDCALPtEtaPhiVeryHigh"), clusterObservable, cluster.eta(), cluster.phi());
      }
      if (isTrigger(TriggerType_t::kEmcalGammaHigh) && isClusterInEmcal(cluster)) { // Only fill EMCAL clusters
        registry.fill(HIST("hSelectedGammaEMCALPtEtaPhi"), clusterObservable, cluster.eta(), cluster.phi());
      }
      if (isTrigger(TriggerType_t::kDcalGammaHigh) && !isClusterInEmcal(cluster)) { // Only fill DCAL clusters
        registry.fill(HIST("hSelectedGammaDCALPtEtaPhi"), clusterObservable, cluster.eta(), cluster.phi());
      }
      if (isTrigger(TriggerType_t::kEmcalGammaLow) && isClusterInEmcal(cluster)) { // Only fill EMCAL clusters
        registry.fill(HIST("hSelectedGammaEMCALPtEtaPhiLow"), clusterObservable, cluster.eta(), cluster.phi());
      }
      if (isTrigger(TriggerType_t::kDcalGammaLow) && !isClusterInEmcal(cluster)) { // Only fill DCAL clusters
        registry.fill(HIST("hSelectedGammaDCALPtEtaPhiLow"), clusterObservable, cluster.eta(), cluster.phi());
      }
      if (isTrigger(TriggerType_t::kEmcalGammaVeryLow) && isClusterInEmcal(cluster)) { // Only fill EMCAL clusters
        registry.fill(HIST("hSelectedGammaEMCALPtEtaPhiVeryLow"), clusterObservable, cluster.eta(), cluster.phi());
      }
      if (isTrigger(TriggerType_t::kDcalGammaVeryLow) && !isClusterInEmcal(cluster)) { // Only fill DCAL clusters
        registry.fill(HIST("hSelectedGammaDCALPtEtaPhiVeryLow"), clusterObservable, cluster.eta(), cluster.phi());
      }
    } // for clusters
    if (!b_doLightOutput) {
      if (maxClusterObservableEMCAL > 0) {
        registry.fill(HIST("hClusterEMCALMaxPtEtaPhi"), maxClusterObservableEMCAL, maxClusterEMCAL.eta(), maxClusterEMCAL.phi());
        if (hwtrg.test(EMCALHardwareTrigger::TRG_MB)) {
          registry.fill(HIST("hClusterEMCALMaxPtEtaPhiMinBias"), maxClusterObservableEMCAL, maxClusterEMCAL.eta(), maxClusterEMCAL.phi());
        }
        if (hwtrg.test(EMCALHardwareTrigger::TRG_EMC7)) {
          registry.fill(HIST("hClusterEMCALMaxPtEtaPhiLevel0"), maxClusterObservableEMCAL, maxClusterEMCAL.eta(), maxClusterEMCAL.phi());
        }
        for (const auto& cluster : analysedClusters) {
          if (cluster.mEMCALcluster) {
            registry.fill(HIST("hClusterEMCALMaxPtClusterEMCALPt"), maxClusterObservableEMCAL, cluster.mTriggerObservable);
          }
        }
        if (isTrigger(TriggerType_t::kEmcalJetFull) || isTrigger(TriggerType_t::kEmcalJetNeutral)) {
          registry.fill(HIST("hSelectedClusterMaxPtEtaPhi"), maxClusterObservableEMCAL, maxClusterEMCAL.eta(), maxClusterEMCAL.phi());
        }
        if (isTrigger(TriggerType_t::kEmcalGammaVeryHigh)) {
          registry.fill(HIST("hSelectedGammaEMCALMaxPtEtaPhiVeryHigh"), maxClusterObservableEMCAL, maxClusterEMCAL.eta(), maxClusterEMCAL.phi());
        }
        if (isTrigger(TriggerType_t::kEmcalGammaHigh)) {
          registry.fill(HIST("hSelectedGammaEMCALMaxPtEtaPhi"), maxClusterObservableEMCAL, maxClusterEMCAL.eta(), maxClusterEMCAL.phi());
        }
        if (isTrigger(TriggerType_t::kEmcalGammaLow)) {
          registry.fill(HIST("hSelectedGammaEMCALMaxPtEtaPhiLow"), maxClusterObservableEMCAL, maxClusterEMCAL.eta(), maxClusterEMCAL.phi());
        }
        if (isTrigger(TriggerType_t::kEmcalGammaVeryLow)) {
          registry.fill(HIST("hSelectedGammaEMCALMaxPtEtaPhiVeryLow"), maxClusterObservableEMCAL, maxClusterEMCAL.eta(), maxClusterEMCAL.phi());
        }
      } else {
        // count events where max. cluster has not been found
        registry.fill(HIST("hEventsNoMaxClusterEMCAL"), 1.);
        if (hwtrg.test(EMCALHardwareTrigger::TRG_MB)) {
          registry.fill(HIST("hEventsNoMaxClusterEMCALMinBias"), 1.);
        }
        if (hwtrg.test(EMCALHardwareTrigger::TRG_EMC7)) {
          registry.fill(HIST("hEventsNoMaxClusterEMCALLevel0"), 1.);
        }
      }
      if (maxClusterObservableDCAL > 0) {
        registry.fill(HIST("hClusterDCALMaxPtEtaPhi"), maxClusterObservableDCAL, maxClusterDCAL.eta(), maxClusterDCAL.phi());
        if (hwtrg.test(EMCALHardwareTrigger::TRG_MB)) {
          registry.fill(HIST("hClusterDCALMaxPtEtaPhiMinBias"), maxClusterObservableDCAL, maxClusterDCAL.eta(), maxClusterDCAL.phi());
        }
        if (hwtrg.test(EMCALHardwareTrigger::TRG_DMC7)) {
          registry.fill(HIST("hClusterDCALMaxPtEtaPhiLevel0"), maxClusterObservableDCAL, maxClusterDCAL.eta(), maxClusterDCAL.phi());
        }
        for (const auto& cluster : analysedClusters) {
          if (!cluster.mEMCALcluster) {
            registry.fill(HIST("hClusterDCALMaxPtClusterDCALPt"), maxClusterObservableDCAL, cluster.mTriggerObservable);
          }
        }
        if (isTrigger(TriggerType_t::kDcalGammaVeryHigh)) {
          registry.fill(HIST("hSelectedGammaDCALMaxPtEtaPhiVeryHigh"), maxClusterObservableDCAL, maxClusterDCAL.eta(), maxClusterDCAL.phi());
        }
        if (isTrigger(TriggerType_t::kDcalGammaHigh)) {
          registry.fill(HIST("hSelectedGammaDCALMaxPtEtaPhi"), maxClusterObservableDCAL, maxClusterDCAL.eta(), maxClusterDCAL.phi());
        }
        if (isTrigger(TriggerType_t::kDcalGammaLow)) {
          registry.fill(HIST("hSelectedGammaDCALMaxPtEtaPhiLow"), maxClusterObservableDCAL, maxClusterDCAL.eta(), maxClusterDCAL.phi());
        }
        if (isTrigger(TriggerType_t::kDcalGammaVeryLow)) {
          registry.fill(HIST("hSelectedGammaDCALMaxPtEtaPhiVeryLow"), maxClusterObservableDCAL, maxClusterDCAL.eta(), maxClusterDCAL.phi());
        }
      } else {
        registry.fill(HIST("hEventsNoMaxClusterDCAL"), 1.);
        // count events where max. cluster has not been found
        if (hwtrg.test(EMCALHardwareTrigger::TRG_MB)) {
          registry.fill(HIST("hEventsNoMaxClusterDCALMinBias"), 1.);
        }
        if (hwtrg.test(EMCALHardwareTrigger::TRG_DMC7)) {
          registry.fill(HIST("hEventsNoMaxClusterDCALLevel0"), 1.);
        }
      }
    }
    return std::make_pair(maxClusterObservableEMCAL, maxClusterObservableDCAL);
  }

  template <typename JetCollection>
  std::pair<std::vector<typename JetCollection::iterator>, std::vector<typename JetCollection::iterator>> fillJetQA(const JetCollection& jets, aod::JetTracks const& /*tracks*/, selectedClusters const& /*clusters*/, std::bitset<EMCALHardwareTrigger::TRG_NTriggers> /*hwtrg*/, const std::bitset<TriggerType_t::kNTriggers>& triggerstatus)
  {
    auto isTrigger = [&triggerstatus](TriggerType_t triggertype) -> bool {
      return triggerstatus.test(triggertype);
    };

    std::vector<typename JetCollection::iterator> vecMaxJet;
    std::vector<typename JetCollection::iterator> vecMaxJetNoFiducial;

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
      if (!b_doLightOutput) {
        for (const auto& track : jet.template tracks_as<aod::JetTracks>()) {
          auto trackPt = track.pt();
          auto chargeFrag = track.px() * jet.px() + track.py() * jet.py() + track.pz() * jet.pz();
          chargeFrag /= (jet.p() * jet.p());
          auto z = trackPt / jetPt;
          auto dphi = jet.phi() - track.phi();
          dphi = check_dphi(dphi);
          auto dR = std::sqrt(dphi * dphi + std::pow(jet.eta() - track.eta(), 2));
          zTheta += z * dR / jetR;
          zSqTheta += z * z * dR / jetR;
          zThetaSq += z * dR * dR / (jetR * jetR);
          ptD += z * z;
          registry.fill(HIST("hJetRPtTrackPt"), jetR, jetPt, trackPt);
          registry.fill(HIST("hJetRPtChargeFrag"), jetR, jetPt, chargeFrag);
          if (isTrigger(TriggerType_t::kEmcalJetFull) || isTrigger(TriggerType_t::kEmcalJetNeutral)) {
            registry.fill(HIST("hSelectedJetRPtTrackPt"), jetR, jetPt, trackPt);
            registry.fill(HIST("hSelectedJetRPtChargeFrag"), jetR, jetPt, chargeFrag);
          }
          if (isTrigger(TriggerType_t::kEmcalJetFullLow) || isTrigger(TriggerType_t::kEmcalJetNeutralLow)) {
            registry.fill(HIST("hSelectedJetLowRPtTrackPt"), jetR, jetPt, trackPt);
            registry.fill(HIST("hSelectedJetLowRPtChargeFrag"), jetR, jetPt, chargeFrag);
          }
        } // for tracks in jet
      }

      // auto clustersInJet = jetClusterConstituents.sliceBy(perJetClusterConstituents, jet.globalIndex());
      // for (const auto& clusterList : clustersInJet) {
      if (!b_doLightOutput) {
        for (const auto& cluster : jet.template clusters_as<selectedClusters>()) {
          auto clusterPt = cluster.energy() / std::cosh(cluster.eta());
          neutralEnergyFraction += cluster.energy();
          auto z = clusterPt / jetPt;
          auto dphi = jet.phi() - cluster.phi();
          dphi = check_dphi(dphi);
          auto dR = std::sqrt(dphi * dphi + std::pow(jet.eta() - cluster.eta(), 2));
          zTheta += z * dR / jetR;
          zSqTheta += z * z * dR / jetR;
          zThetaSq += z * dR * dR / (jetR * jetR);
          ptD += z * z;
          registry.fill(HIST("hJetRPtClusterPt"), jetR, jetPt, clusterPt);
          if (isTrigger(TriggerType_t::kEmcalJetFull) || isTrigger(TriggerType_t::kEmcalJetNeutral)) {
            registry.fill(HIST("hSelectedJetRPtClusterPt"), jetR, jetPt, clusterPt);
          }
          if (isTrigger(TriggerType_t::kEmcalJetFullLow) || isTrigger(TriggerType_t::kEmcalJetNeutralLow)) {
            registry.fill(HIST("hSelectedJetLowRPtClusterPt"), jetR, jetPt, clusterPt);
          }
        } // for clusters in jet
      }
      neutralEnergyFraction /= jet.energy();
      ptD = std::sqrt(ptD);

      // Fillng histograms
      registry.fill(HIST("hJetRPtEta"), jetR, jetPt, jet.eta());
      registry.fill(HIST("hJetRPtPhi"), jetR, jetPt, jet.phi());
      if (!b_doLightOutput) {
        registry.fill(HIST("hJetRPtPtd"), jetR, jetPt, ptD);
        registry.fill(HIST("hJetRPtNEF"), jetR, jetPt, neutralEnergyFraction);
        registry.fill(HIST("hJetRPtZTheta"), jetR, jetPt, zTheta);
        registry.fill(HIST("hJetRPtZSqTheta"), jetR, jetPt, zSqTheta);
        registry.fill(HIST("hJetRPtZThetaSq"), jetR, jetPt, zThetaSq);
      }
      registry.get<THn>(HIST("jetRPtEtaPhi"))->Fill(jetR, jetPt, jet.eta(), jet.phi());
      if (isTrigger(TriggerType_t::kEmcalJetFull) || isTrigger(TriggerType_t::kEmcalJetNeutral)) {
        registry.fill(HIST("hSelectedJetRPtEta"), jetR, jetPt, jet.eta());
        registry.fill(HIST("hSelectedJetRPtPhi"), jetR, jetPt, jet.phi());
        if (!b_doLightOutput) {
          registry.fill(HIST("hSelectedJetRPtPtd"), jetR, jetPt, ptD);
          registry.fill(HIST("hSelectedJetRPtNEF"), jetR, jetPt, neutralEnergyFraction);
          registry.fill(HIST("hSelectedJetRPtZTheta"), jetR, jetPt, zTheta);
          registry.fill(HIST("hSelectedJetRPtZSqTheta"), jetR, jetPt, zSqTheta);
          registry.fill(HIST("hSelectedJetRPtZThetaSq"), jetR, jetPt, zThetaSq);
        }
      }
      if (isTrigger(TriggerType_t::kEmcalJetFullLow) || isTrigger(TriggerType_t::kEmcalJetNeutralLow)) {
        registry.fill(HIST("hSelectedJetLowRPtEta"), jetR, jetPt, jet.eta());
        registry.fill(HIST("hSelectedJetLowRPtPhi"), jetR, jetPt, jet.phi());
        if (!b_doLightOutput) {
          registry.fill(HIST("hSelectedJetLowRPtPtd"), jetR, jetPt, ptD);
          registry.fill(HIST("hSelectedJetLowRPtNEF"), jetR, jetPt, neutralEnergyFraction);
          registry.fill(HIST("hSelectedJetLowRPtZTheta"), jetR, jetPt, zTheta);
          registry.fill(HIST("hSelectedJetLowRPtZSqTheta"), jetR, jetPt, zSqTheta);
          registry.fill(HIST("hSelectedJetLowRPtZThetaSq"), jetR, jetPt, zThetaSq);
        }
      }
    } // for jets
    return std::make_pair(vecMaxJet, vecMaxJetNoFiducial);
  }

  using JetCollisionsTable = soa::Join<aod::JetCollisions, aod::JFullTrigSels>;

  template <typename JetCollection>
  void runQA(collisionWithTrigger const& collision,
             JetCollection const& jets,
             aod::JetTracks const& tracks,
             selectedClusters const& clusters)
  {
    std::bitset<TriggerType_t::kNTriggers> triggerstatus;
    auto isTrigger = [&triggerstatus](TriggerType_t triggertype) -> bool {
      return triggerstatus.test(triggertype);
    };
    auto setTrigger = [&triggerstatus](TriggerType_t triggertype) {
      triggerstatus.set(triggertype, true);
    };

    std::bitset<EMCALHardwareTrigger::TRG_NTriggers> hardwaretriggers;
    if (collision.alias_bit(triggerAliases::kTVXinEMC)) {
      hardwaretriggers.set(EMCALHardwareTrigger::TRG_MB);
    }
    if (collision.alias_bit(triggerAliases::kEMC7)) {
      hardwaretriggers.set(EMCALHardwareTrigger::TRG_EMC7);
    }
    if (collision.alias_bit(triggerAliases::kDMC7)) {
      hardwaretriggers.set(EMCALHardwareTrigger::TRG_DMC7);
    }

    fillEventSelectionCounter(0);
    setTrigger(TriggerType_t::kMinBias);

    if (hasEMCAL(collision)) {
      fillEventSelectionCounter(1);
      setTrigger(TriggerType_t::kEmcalAny);
    }

    // fill event counters and correlation matrix without constraint on the EMCAL trigger flag
    if (collision.alias_bit(triggerAliases::kTVXinEMC)) {
      fillEventSelectionCounter(2);
      setTrigger(TriggerType_t::kEmcalMB);
    }

    if (jetderiveddatautilities::selectFullTrigger(collision, jetderiveddatautilities::JTrigSelFull::jetFullHighPt)) {
      fillEventSelectionCounter(3);
      setTrigger(TriggerType_t::kEmcalJetFull);
    }
    if (jetderiveddatautilities::selectFullTrigger(collision, jetderiveddatautilities::JTrigSelFull::jetFullLowPt)) {
      fillEventSelectionCounter(4);
      setTrigger(TriggerType_t::kEmcalJetFullLow);
    }
    if (jetderiveddatautilities::selectFullTrigger(collision, jetderiveddatautilities::JTrigSelFull::jetNeutralHighPt)) {
      fillEventSelectionCounter(5);
      setTrigger(TriggerType_t::kEmcalJetNeutral);
    }
    if (jetderiveddatautilities::selectFullTrigger(collision, jetderiveddatautilities::JTrigSelFull::jetNeutralLowPt)) {
      fillEventSelectionCounter(6);
      setTrigger(TriggerType_t::kEmcalJetNeutralLow);
    }
    if (jetderiveddatautilities::selectFullTrigger(collision, jetderiveddatautilities::JTrigSelFull::gammaVeryHighPtEMCAL)) {
      fillEventSelectionCounter(7);
      setTrigger(TriggerType_t::kEmcalGammaVeryHigh);
    }
    if (jetderiveddatautilities::selectFullTrigger(collision, jetderiveddatautilities::JTrigSelFull::gammaVeryHighPtDCAL)) {
      fillEventSelectionCounter(8);
      setTrigger(TriggerType_t::kDcalGammaVeryHigh);
    }
    if (jetderiveddatautilities::selectFullTrigger(collision, jetderiveddatautilities::JTrigSelFull::gammaHighPtEMCAL)) {
      fillEventSelectionCounter(9);
      setTrigger(TriggerType_t::kEmcalGammaHigh);
    }
    if (jetderiveddatautilities::selectFullTrigger(collision, jetderiveddatautilities::JTrigSelFull::gammaHighPtDCAL)) {
      fillEventSelectionCounter(10);
      setTrigger(TriggerType_t::kDcalGammaHigh);
    }
    if (jetderiveddatautilities::selectFullTrigger(collision, jetderiveddatautilities::JTrigSelFull::gammaLowPtEMCAL)) {
      fillEventSelectionCounter(11);
      setTrigger(TriggerType_t::kEmcalGammaLow);
    }
    if (jetderiveddatautilities::selectFullTrigger(collision, jetderiveddatautilities::JTrigSelFull::gammaLowPtDCAL)) {
      fillEventSelectionCounter(12);
      setTrigger(TriggerType_t::kDcalGammaLow);
    }
    if (jetderiveddatautilities::selectFullTrigger(collision, jetderiveddatautilities::JTrigSelFull::gammaVeryLowPtEMCAL)) {
      fillEventSelectionCounter(13);
      setTrigger(TriggerType_t::kEmcalGammaVeryLow);
    }
    if (jetderiveddatautilities::selectFullTrigger(collision, jetderiveddatautilities::JTrigSelFull::gammaVeryLowPtDCAL)) {
      fillEventSelectionCounter(14);
      setTrigger(TriggerType_t::kDcalGammaVeryLow);
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
      if (!isTrigger(TriggerType_t::kEmcalAny)) {
        return; // Only consider events where EMCAL is live
      }
    }

    auto [maxClusterObservableEMCAL, maxClusterObservableDCAL] = fillGammaQA(clusters, hardwaretriggers, triggerstatus);
    auto [vecMaxJet, vecMaxJetNoFiducial] = fillJetQA(jets, tracks, clusters, hardwaretriggers, triggerstatus);

    std::array<bool, 5> foundMaxJet;
    std::fill(foundMaxJet.begin(), foundMaxJet.end(), false);
    if (!b_doLightOutput) {
      for (const auto& maxJet : vecMaxJet) {
        double jetR = maxJet.r() * 1e-2, jetPt = maxJet.pt(), jetEta = maxJet.eta(), jetPhi = maxJet.phi();
        foundMaxJet[static_cast<int>(maxJet.r() * 1e-1) - 2] = true;
        registry.fill(HIST("hJetRMaxPtEta"), jetR, jetPt, jetEta);
        registry.fill(HIST("hJetRMaxPtPhi"), jetR, jetPt, jetPhi);
        if (hardwaretriggers.test(EMCALHardwareTrigger::TRG_MB)) {
          registry.fill(HIST("hJetRMaxPtEtaMinBias"), jetR, jetPt, jetEta);
          registry.fill(HIST("hJetRMaxPtPhiMinBias"), jetR, jetPt, jetPhi);
        }
        if (hardwaretriggers.test(EMCALHardwareTrigger::TRG_EMC7)) {
          registry.fill(HIST("hJetRMaxPtEtaLevel0"), jetR, jetPt, jetEta);
          registry.fill(HIST("hJetRMaxPtPhiLevel0"), jetR, jetPt, jetPhi);
        }
        // hJetRMaxPtEtaPhi->Fill(jetR, jetPt, jetEta, jetPhi);
        registry.get<THn>(HIST("jetRMaxPtEtaPhi"))->Fill(jetR, jetPt, jetEta, jetPhi);
        if (isTrigger(TriggerType_t::kEmcalJetFull) || isTrigger(TriggerType_t::kEmcalJetNeutral)) {
          registry.fill(HIST("hSelectedJetRMaxPtEta"), jetR, jetPt, jetEta);
          registry.fill(HIST("hSelectedJetRMaxPtPhi"), jetR, jetPt, jetPhi);
        }
        if (isTrigger(TriggerType_t::kEmcalJetFullLow) || isTrigger(TriggerType_t::kEmcalJetNeutralLow)) {
          registry.fill(HIST("hSelectedJetLowRMaxPtEta"), jetR, jetPt, jetEta);
          registry.fill(HIST("hSelectedJetLowRMaxPtPhi"), jetR, jetPt, jetPhi);
        }
        if (maxClusterObservableEMCAL > 0) {
          registry.fill(HIST("hJetRMaxPtClusterMaxPt"), jetR, jetPt, maxClusterObservableEMCAL);
          if (isTrigger(TriggerType_t::kEmcalJetFull) || isTrigger(TriggerType_t::kEmcalJetNeutral)) {
            registry.fill(HIST("hSelectedJetRMaxPtClusterMaxPt"), jetR, jetPt, maxClusterObservableEMCAL);
          }
        } // if maxClusterPt
        if (maxJet.r() == std::round(f_jetR * 100)) {
          for (const auto& jet : jets) {
            if (isJetInEmcal(jet)) {
              registry.fill(HIST("hJetRMaxPtJetPt"), jet.r() * 1e-2, jetPt, jet.pt());
            }
          } // for jets
        } // if maxJet.r() == std::round(f_jetR * 100)
      } // for maxJet
    }
    // Fill counters for events without max jets
    if (!b_doLightOutput) {
      for (std::size_t ir = 0; ir < foundMaxJet.size(); ir++) {
        if (!foundMaxJet[ir]) {
          double rval = static_cast<double>(ir) / 10.;
          registry.fill(HIST("hEventsNoMaxJet"), rval);
          if (hardwaretriggers.test(EMCALHardwareTrigger::TRG_MB)) {
            registry.fill(HIST("hEventsNoMaxJetMinBias"), rval);
          }
          if (hardwaretriggers.test(EMCALHardwareTrigger::TRG_EMC7)) {
            registry.fill(HIST("hEventsNoMaxJetLevel0"), rval);
          }
        }
      }
    }

    if (!b_doLightOutput) {
      for (const auto& maxJet : vecMaxJetNoFiducial) {
        double jetR = maxJet.r() * 1e-2, jetPt = maxJet.pt(), jetEta = maxJet.eta(), jetPhi = maxJet.phi();
        // hJetRMaxPtEtaPhiNoFiducial->Fill(jetR, jetPt, jetEta, jetPhi);
        registry.get<THn>(HIST("jetRMaxPtEtaPhiNoFiducial"))->Fill(jetR, jetPt, jetEta, jetPhi);
        if (maxJet.r() == std::round(f_jetR * 100)) {
          for (const auto& jet : jets) {
            registry.fill(HIST("hJetRMaxPtJetPtNoFiducial"), jet.r() * 1e-2, jetPt, jet.pt());
          } // for jets
        } // if maxJet.r() == std::round(f_jetR * 100)
      } // for maxjet no fiducial
    }
  } // process

  void processFullJets(collisionWithTrigger const& collision,
                       fullJetInfos const& jets,
                       aod::JetTracks const& tracks,
                       selectedClusters const& clusters)
  {
    runQA(collision, jets, tracks, clusters);
  }
  PROCESS_SWITCH(JetTriggerQA, processFullJets, "Run QA for full jets", true);

  void processNeutralJets(collisionWithTrigger const& collision,
                          neutralJetInfos const& jets,
                          aod::JetTracks const& tracks,
                          selectedClusters const& clusters)
  {
    runQA(collision, jets, tracks, clusters);
  }
  PROCESS_SWITCH(JetTriggerQA, processNeutralJets, "Run QA for neutral jets", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<JetTriggerQA>(cfgc, TaskName{"jet-full-trigger-qa"})};
}
