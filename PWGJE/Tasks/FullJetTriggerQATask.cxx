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

  OutputObj<TH1I> hProcessedEvents{"hProcessedEvents"};

  OutputObj<TH3F> hJetRPtEta{"hJetRPtEta"};
  OutputObj<TH3F> hJetRPtPhi{"hJetRPtPhi"};
  OutputObj<TH2F> hClusterPtEta{"hClusterPtEta"};
  OutputObj<TH2F> hClusterPtPhi{"hClusterPtPhi"};
  OutputObj<TH3F> hJetRPtTrackPt{"hJetRPtTrackPt"};
  OutputObj<TH3F> hJetRPtClusterPt{"hJetRPtClusterPt"};
  OutputObj<TH3F> hJetRPtPtd{"hJetRPtPtd"};
  OutputObj<TH3F> hJetRPtChargeFrag{"hJetRPtChargeFrag"};
  OutputObj<TH3F> hJetRPtNEF{"hJetRPtNEF"};
  OutputObj<TH3F> hJetRPtZTheta{"hJetRPtZTheta"};
  OutputObj<TH3F> hJetRPtZSqTheta{"hJetRPtZSqTheta"};
  OutputObj<TH3F> hJetRPtZThetaSq{"hJetRPtZThetaSq"};

  OutputObj<TH3F> hSelectedJetRPtEta{"hSelectedJetRPtEta"};
  OutputObj<TH3F> hSelectedJetRPtPhi{"hSelectedJetRPtPhi"};
  OutputObj<TH2F> hSelectedClusterPtEta{"hSelectedClusterPtEta"};
  OutputObj<TH2F> hSelectedClusterPtPhi{"hSelectedClusterPtPhi"};
  OutputObj<TH3F> hSelectedJetRPtTrackPt{"hSelectedJetRPtTrackPt"};
  OutputObj<TH3F> hSelectedJetRPtClusterPt{"hSelectedJetRPtClusterPt"};
  OutputObj<TH3F> hSelectedJetRPtPtd{"hSelectedJetRPtPtd"};
  OutputObj<TH3F> hSelectedJetRPtChargeFrag{"hSelectedJetRPtChargeFrag"};
  OutputObj<TH3F> hSelectedJetRPtNEF{"hSelectedJetRPtNEF"};
  OutputObj<TH3F> hSelectedJetRPtZTheta{"hSelectedJetRPtZTheta"};
  OutputObj<TH3F> hSelectedJetRPtZSqTheta{"hSelectedJetRPtZSqTheta"};
  OutputObj<TH3F> hSelectedJetRPtZThetaSq{"hSelectedJetRPtZThetaSq"};

  Configurable<float> f_jetPtMin{"f_jetPtMin", 0.0, "minimum jet pT cut"};
  Configurable<float> f_SD_zCut{"f_SD_zCut", 0.1, "soft drop z cut"};
  Configurable<float> f_SD_beta{"f_SD_beta", 0.0, "soft drop beta"};
  Configurable<float> f_jetR{"f_jetR", 0.4, "jet resolution parameter"};
  Configurable<std::vector<float>> f_ang_kappa{"f_ang_kappa", {1.0, 1.0, 2.0}, "angularity momentum exponent"};
  Configurable<std::vector<float>> f_ang_alpha{"f_ang_alpha", {1.0, 2.0, 1.0}, "angularity angle exponent"};
  Configurable<bool> b_DoConstSub{"b_DoConstSub", false, "do constituent subtraction"};
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
    Float_t kMaxPt = 200.;
    Int_t nPhiBins = 18 * 8;
    Float_t kMinPhi = 0.;
    Float_t kMaxPhi = 2. * TMath::Pi();
    Int_t nEtaBins = 100;
    Float_t kMinEta = -1.;
    Float_t kMaxEta = 1.;
    Int_t nRBins = 6;
    Float_t kMinR = 0.05;
    Float_t kMaxR = 0.65;

    hProcessedEvents.setObject(new TH1I("hProcessedEvents", "Processed events", 3, -0.5, 2.5));
    hProcessedEvents->GetXaxis()->SetBinLabel(1, "MB");
    hProcessedEvents->GetXaxis()->SetBinLabel(2, "EMC");
    hProcessedEvents->GetXaxis()->SetBinLabel(3, "Selected");

    // Histograms for events where the EMCAL is live
    hJetRPtEta.setObject(new TH3F("hJetRPtEta", "Jets #it{p}_{T} and #eta;#it{R};#it{p}_{T};#eta", nRBins, kMinR, kMaxR, nPtBins, kMinPt, kMaxPt, nEtaBins, kMinEta, kMaxEta));
    hJetRPtPhi.setObject(new TH3F("hJetRPtPhi", "Jets #it{p}_{T} and #phi;#it{R};#it{p}_{T};#phi", nRBins, kMinR, kMaxR, nPtBins, kMinPt, kMaxPt, nPhiBins, kMinPhi, kMaxPhi));
    hClusterPtEta.setObject(new TH2F("hClusterPtEta", "Cluster #it{p}_{T} and #eta;#it{p}_{T};#eta", nPtBins, kMinPt, kMaxPt / 2, nEtaBins, kMinEta, kMaxEta));
    hClusterPtPhi.setObject(new TH2F("hClusterPtPhi", "Cluster #it{p}_{T} and #phi;#it{p}_{T};#phi", nPtBins, kMinPt, kMaxPt / 2, nPhiBins, kMinPhi, kMaxPhi));
    hJetRPtTrackPt.setObject(new TH3F("hJetRPtTrackPt", "Jets;#it{p}_{T};#it{R};#it{p}_{T}^{track}", nRBins, kMinR, kMaxR, nPtBins, kMinPt, kMaxPt, nPtBins, kMinPt, kMaxPt / 2));
    hJetRPtClusterPt.setObject(new TH3F("hJetRPtClusterPt", "Jets;#it{R};#it{p}_{T};#it{p}_{T}^{clus}", nRBins, kMinR, kMaxR, nPtBins, kMinPt, kMaxPt, nPtBins, kMinPt, kMaxPt / 2));
    hJetRPtPtd.setObject(new TH3F("hJetRPtPtd", "Jets;#it{R};#it{p}_{T};ptD", nRBins, kMinR, kMaxR, nPtBins, kMinPt, kMaxPt, nPtBins / 2, 0., 1.));
    hJetRPtChargeFrag.setObject(new TH3F("hJetRPtChargeFrag", "Jets;#it{R};#it{p}_{T};z", nRBins, kMinR, kMaxR, nPtBins, kMinPt, kMaxPt, nPtBins / 2, 0., 1.));
    hJetRPtNEF.setObject(new TH3F("hJetRPtNEF", "Jets;#it{R};#it{p}_{T};NEF", nRBins, kMinR, kMaxR, nPtBins, kMinPt, kMaxPt, nPtBins / 2, 0., 1.));
    hJetRPtZTheta.setObject(new TH3F("hJetRPtZTheta", "Jets;#it{R};#it{p}_{T};z#theta", nRBins, kMinR, kMaxR, nPtBins, kMinPt, kMaxPt, nPtBins / 2, 0., 1.));
    hJetRPtZSqTheta.setObject(new TH3F("hJetRPtZSqTheta", "Jets;#it{R};#it{p}_{T};z^{2} #theta", nRBins, kMinR, kMaxR, nPtBins, kMinPt, kMaxPt, nPtBins / 2, 0., 1.));
    hJetRPtZThetaSq.setObject(new TH3F("hJetRPtZThetaSq", "Jets;#it{R};#it{p}_{T};z #theta^{2}", nRBins, kMinR, kMaxR, nPtBins, kMinPt, kMaxPt, nPtBins / 2, 0., 1.));

    // Histograms for triggered events
    hSelectedJetRPtEta.setObject(new TH3F("hSelectedJetRPtEta", "Selected jets #it{p}_{T} and #eta;#it{R};#it{p}_{T};#eta", nRBins, kMinR, kMaxR, nPtBins, kMinPt, kMaxPt, nEtaBins, kMinEta, kMaxEta));
    hSelectedJetRPtPhi.setObject(new TH3F("hSelectedJetRPtPhi", "Selected jets #it{p}_{T} and #phi;#it{R};#it{p}_{T};#phi", nRBins, kMinR, kMaxR, nPtBins, kMinPt, kMaxPt, nPhiBins, kMinPhi, kMaxPhi));
    hSelectedClusterPtEta.setObject(new TH2F("hSelectedClusterPtEta", "Selected Cluster #it{p}_{T} and #eta;#it{p}_{T};#eta", nPtBins, kMinPt, kMaxPt / 2, nEtaBins, kMinEta, kMaxEta));
    hSelectedClusterPtPhi.setObject(new TH2F("hSelectedClusterPtPhi", "Selected Cluster #it{p}_{T} and #phi;#it{p}_{T};#phi", nPtBins, kMinPt, kMaxPt / 2, nPhiBins, kMinPhi, kMaxPhi));
    hSelectedJetRPtTrackPt.setObject(new TH3F("hSelectedJetRPtTrackPt", "Selected jets;#it{R};#it{p}_{T};#it{p}_{T}^{track}", nRBins, kMinR, kMaxR, nPtBins, kMinPt, kMaxPt, nPtBins, kMinPt, kMaxPt / 2));
    hSelectedJetRPtClusterPt.setObject(new TH3F("hSelectedJetRPtClusterPt", "Selected jets;#it{R};#it{p}_{T};#it{p}_{T}^{clus}", nRBins, kMinR, kMaxR, nPtBins, kMinPt, kMaxPt, nPtBins, kMinPt, kMaxPt / 2));
    hSelectedJetRPtPtd.setObject(new TH3F("hSelectedJetRPtPtd", "Selected jets;#it{R};#it{p}_{T};ptD", nRBins, kMinR, kMaxR, nPtBins, kMinPt, kMaxPt, nPtBins / 2, 0., 1.));
    hSelectedJetRPtChargeFrag.setObject(new TH3F("hSelectedJetRPtChargeFrag", "Selected jets;#it{R};#it{p}_{T};z", nRBins, kMinR, kMaxR, nPtBins, kMinPt, kMaxPt, nPtBins / 2, 0., 1.));
    hSelectedJetRPtNEF.setObject(new TH3F("hSelectedJetRPtNEF", "Selected jets;#it{R};#it{p}_{T};NEF", nRBins, kMinR, kMaxR, nPtBins, kMinPt, kMaxPt, nPtBins / 2, 0., 1.));
    hSelectedJetRPtZTheta.setObject(new TH3F("hSelectedJetRPtZTheta", "Selected jets;#it{R};#it{p}_{T};z#theta", nRBins, kMinR, kMaxR, nPtBins, kMinPt, kMaxPt, nPtBins / 2, 0., 1.));
    hSelectedJetRPtZSqTheta.setObject(new TH3F("hSelectedJetRPtZSqTheta", "Selected jets;#it{R};#it{p}_{T};z^{2} #theta", nRBins, kMinR, kMaxR, nPtBins, kMinPt, kMaxPt, nPtBins / 2, 0., 1.));
    hSelectedJetRPtZThetaSq.setObject(new TH3F("hSelectedJetRPtZThetaSq", "Selected jets;#it{R};#it{p}_{T};z #theta^{2}", nRBins, kMinR, kMaxR, nPtBins, kMinPt, kMaxPt, nPtBins / 2, 0., 1.));
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

  void process(soa::Join<aod::Collisions, aod::EvSels, aod::FullJetFilters>::iterator const& collision,
               aod::Jets const& jets,
               aod::JetTrackConstituents const& jetTrackConstituents,
               aod::JetClusterConstituents const& jetClusterConstituents,
               aod::Tracks const& tracks,
               selectedClusters const& clusters)
  {
    hProcessedEvents->Fill(0);
    if (!collision.alias()[kTVXinEMC]) {
      return; // Only consider events where EMCAL is live
    }
    hProcessedEvents->Fill(1);

    bool isEvtSelected = false;
    if (collision.hasJetFullHighPt()) {
      isEvtSelected = true;
    }
    if (isEvtSelected) {
      hProcessedEvents->Fill(2);
    }

    for (const auto& jet : jets) {
      float neutralEnergyFraction = 0., ptD = 0.;
      float zTheta = 0., zSqTheta = 0., zThetaSq = 0; // Jet angularities (1,1), (1,2), (2,1)
      double jetPt = jet.pt(), jetR = jet.r() * 1e-2;

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

    for (const auto& cluster : clusters) {
      double clusterPt = cluster.energy() / std::cosh(cluster.eta());
      hClusterPtEta->Fill(clusterPt, cluster.eta());
      hClusterPtPhi->Fill(clusterPt, cluster.phi());
      if (isEvtSelected) {
        hSelectedClusterPtEta->Fill(clusterPt, cluster.eta());
        hSelectedClusterPtPhi->Fill(clusterPt, cluster.phi());
      }
    } // for clusters
  }   // process
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<JetTriggerQA>(cfgc, TaskName{"jet-full-trigger-qa"})};
}
