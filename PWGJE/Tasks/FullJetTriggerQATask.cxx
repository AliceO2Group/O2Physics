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

  OutputObj<TH2F> hJetPtEta{"hJetPtEta"};
  OutputObj<TH2F> hJetPtPhi{"hJetPtPhi"};
  OutputObj<TH2F> hJetSelectedPtEta{"hJetSelectedPtEta"};
  OutputObj<TH2F> hJetSelectedPtPhi{"hJetSelectedPtPhi"};

  OutputObj<TH2F> hClusterPtEta{"hClusterPtEta"};
  OutputObj<TH2F> hClusterPtPhi{"hClusterPtPhi"};
  OutputObj<TH2F> hClusterSelectedPtEta{"hClusterSelectedPtEta"};
  OutputObj<TH2F> hClusterSelectedPtPhi{"hClusterSelectedPtPhi"};

  OutputObj<TH3F> hJetRPtEta{"hJetRPtEta"};

  OutputObj<TH2F> hEmcJetPtTrackPt{"hEmcJetPtTrackPt"};
  OutputObj<TH2F> hEmcJetPtClusterPt{"hEmcJetPtClusterPt"};
  OutputObj<TH2F> hSelectedJetPtTrackPt{"hSelectedJetPtTrackPt"};
  OutputObj<TH2F> hSelectedJetPtClusterPt{"hSelectedJetPtClusterPt"};
  OutputObj<TH2F> hEmcJetPtPtd{"hEmcJetPtPtd"};

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
    Float_t kMinPt   = 0.;
    Float_t kMaxPt   = 200.;
    Int_t nPhiBins = 18*8;
    Float_t kMinPhi   = 0.;
    Float_t kMaxPhi   = 2.*TMath::Pi();
    Int_t nEtaBins = 100;
    Float_t kMinEta = -1.;
    Float_t kMaxEta =  1.;

    hProcessedEvents.setObject( new TH1I("hProcessedEvents", "Processed events", 3, -0.5, 2.5));
    hProcessedEvents->GetXaxis()->SetBinLabel(1, "MB");
    hProcessedEvents->GetXaxis()->SetBinLabel(2, "EMC");
    hProcessedEvents->GetXaxis()->SetBinLabel(3, "Selected");

    hJetPtEta.setObject( new TH2F("hJetPtEta", "Jet #it{p}_{T} and #eta;#it{p}_{T};#eta",
                                  nPtBins, kMinPt, kMaxPt, nEtaBins, kMinEta, kMaxEta));
    hJetPtPhi.setObject( new TH2F("hJetPtPhi", "Jet #it{p}_{T} and #phi;#it{p}_{T};#phi",
                                  nPtBins, kMinPt, kMaxPt, nPhiBins, kMinPhi, kMaxPhi));
    hJetSelectedPtEta.setObject( new TH2F("hJetSelectedPtEta", "Selected jet #it{p}_{T} and #eta;#it{p}_{T};#eta",
                                  nPtBins, kMinPt, kMaxPt, nEtaBins, kMinEta, kMaxEta));
    hJetSelectedPtPhi.setObject( new TH2F("hJetSelectedPtPhi", "Selected jet #it{p}_{T} and #phi;#it{p}_{T};#phi",
                                  nPtBins, kMinPt, kMaxPt, nPhiBins, kMinPhi, kMaxPhi));

    hClusterPtEta.setObject( new TH2F("hClusterPtEta", "Cluster #it{p}_{T} and #eta;#it{p}_{T};#eta",
                                  nPtBins, kMinPt, kMaxPt/2, nEtaBins, kMinEta, kMaxEta));
    hClusterPtPhi.setObject( new TH2F("hClusterPtPhi", "Cluster #it{p}_{T} and #phi;#it{p}_{T};#phi",
                                  nPtBins, kMinPt, kMaxPt/2, nPhiBins, kMinPhi, kMaxPhi));
    hClusterSelectedPtEta.setObject( new TH2F("hClusterSelectedPtEta", "Selected Cluster #it{p}_{T} and #eta;#it{p}_{T};#eta",
                                  nPtBins, kMinPt, kMaxPt/2, nEtaBins, kMinEta, kMaxEta));
    hClusterSelectedPtPhi.setObject( new TH2F("hClusterSelectedPtPhi", "Selected Cluster #it{p}_{T} and #phi;#it{p}_{T};#phi",
                                  nPtBins, kMinPt, kMaxPt/2, nPhiBins, kMinPhi, kMaxPhi));

    hJetRPtEta.setObject( new TH3F("hJetRPtEta", "Jet #it{R}, #it{p}_{T} and #eta;#it{R};#it{p}_{T};#eta",
                                   6, 0.05, 0.65, nPtBins, kMinPt, kMaxPt, nEtaBins, kMinEta, kMaxEta));

    hEmcJetPtTrackPt.setObject( new TH2F("hEmcJetPtTrackPt", "EMC Jets;#it{p}_{T};#it{p}_{T}^{track}",
                                         nPtBins, kMinPt, kMaxPt, nPtBins, kMinPt, kMaxPt));
    hEmcJetPtClusterPt.setObject( new TH2F("hEmcJetPtClusterPt", "EMC Jets;#it{p}_{T};#it{p}_{T}^{clus}",
                                           nPtBins, kMinPt, kMaxPt, nPtBins, kMinPt, kMaxPt));
    hSelectedJetPtTrackPt.setObject( new TH2F("hSelectedJetPtTrackPt", "Selected Jets;#it{p}_{T};#it{p}_{T}^{track}",
                                         nPtBins, kMinPt, kMaxPt, nPtBins, kMinPt, kMaxPt));
    hSelectedJetPtClusterPt.setObject( new TH2F("hSelectedJetPtClusterPt", "Selected Jets;#it{p}_{T};#it{p}_{T}^{clus}",
                                           nPtBins, kMinPt, kMaxPt, nPtBins, kMinPt, kMaxPt));
    hEmcJetPtPtd.setObject( new TH2F("hEmcJetPtPtd", "EMC Jets;#it{p}_{T};ptD",
                                     nPtBins, kMinPt, kMaxPt, nPtBins, kMinPt, kMaxPt));
  }

  o2::aod::EMCALClusterDefinition clusDef = o2::aod::emcalcluster::getClusterDefinitionFromString(mClusterDefinition.value);
  Filter clusterDefinitionSelection = o2::aod::emcalcluster::definition == static_cast<int>(clusDef);

  void process(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision,
               aod::Jets const& jets,
               aod::JetTrackConstituents const& jetTrackConstituents,
               aod::JetClusterConstituents const& jetClusterConstituents,
               aod::Tracks const& tracks,
               selectedClusters const& clusters,
               aod::FullJetFilters const& fullJetFilters)
  {
    hProcessedEvents->Fill(0);
    if (!collision.alias()[kTVXinEMC]){
      return; // Only consider events where EMCAL is live
    }
    hProcessedEvents->Fill(1);

    bool isEvtSelected = false;
    for (const auto& obj : fullJetFilters){
      if (obj.hasJetFullHighPt()){
        isEvtSelected = true;
      }
    }
    if (isEvtSelected) {
      hProcessedEvents->Fill(2);
    }

    for (const auto& jet : jets){
      float neutralEnergyFraction = 0., ptD = 0., chargeFrag = 0.;
      float zTheta = 0., zSqTheta = 0., zThetaSq = 0; // Jet angularities (1,1), (1,2), (2,1)

      // Slice JetTrackConstituents table to obtain a list of trackId belonging to this jet only (and the same for clusters)
      // This gives us access to all jet substructure information
      auto tracksInJet = jetTrackConstituents.sliceBy(perJetTrackConstituents, jet.globalIndex());
      for (const auto& trackList : tracksInJet) {
        auto trackPt = trackList.track().pt();
        hEmcJetPtTrackPt->Fill(jet.pt(), trackPt);
        if (isEvtSelected){
          hSelectedJetPtTrackPt->Fill(jet.pt(), trackPt);
        }
      } // for tracks in jet
      auto clustersInJet = jetClusterConstituents.sliceBy(perJetClusterConstituents, jet.globalIndex());
      for (const auto& clusterList : clustersInJet) {
        auto clusterPt = clusterList.cluster().energy() / std::cosh(clusterList.cluster().eta());
        neutralEnergyFraction += clusterList.cluster().energy();
        hEmcJetPtClusterPt->Fill(jet.pt(), clusterPt);
        if (isEvtSelected){
          hSelectedJetPtClusterPt->Fill(jet.pt(), clusterPt);
        }
      } // for clusters in jet

      neutralEnergyFraction /= jet.energy();
      // Filling histograms
      hJetPtEta->Fill(jet.pt(), jet.eta());
      hJetPtPhi->Fill(jet.pt(), jet.phi());
      hJetRPtEta->Fill(jet.r() * 1e-2, jet.pt(), jet.eta());
      if (isEvtSelected){
        hJetSelectedPtEta->Fill(jet.pt(), jet.eta());
        hJetSelectedPtPhi->Fill(jet.pt(), jet.phi());
      }
    } // for jets

    for (const auto& cluster : clusters){
      double clusterPt = cluster.energy() / std::cosh(cluster.eta());
      hClusterPtEta->Fill(clusterPt, cluster.eta());
      hClusterPtPhi->Fill(clusterPt, cluster.phi());
      if (isEvtSelected){
        hClusterSelectedPtEta->Fill(clusterPt, cluster.eta());
        hClusterSelectedPtPhi->Fill(clusterPt, cluster.phi());
      }
    } // for clusters
  }
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<JetTriggerQA>(cfgc, TaskName{"jet-full-trigger-qa"})};
}
