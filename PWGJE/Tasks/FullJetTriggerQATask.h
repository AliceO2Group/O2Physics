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

// jet trigger QA task (subscribing to jet finder task)
//
// Author: Gijs van Weelden
//

#ifndef O2_ANALYSIS_JETTRIGGER_QA_H
#define O2_ANALYSIS_JETTRIGGER_QA_H

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
// #include "EventFiltering/PWGJE/jetFilter.cxx"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/Core/JetFinder.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using selectedClusters = o2::soa::Filtered<o2::aod::EMCALClusters>;

// namespace o2::aod
// {
// namespace jettriggerqa
// {
// DECLARE_SOA_COLUMN(Zg, zg, float);
// DECLARE_SOA_COLUMN(Rg, rg, float);
// DECLARE_SOA_COLUMN(Nsd, nsd, float);
// } // namespace jettriggerqa
// DECLARE_SOA_TABLE(JetTriggerQA, "AOD", "JETTRIGGERQA", jettriggerqa::Zg, jettriggerqa::Rg, jettriggerqa::Nsd);
// } // namespace o2::aod

struct JetTriggerQA {
  // Produces<aod::JetTriggerQA> jetTriggerQA;
  // OutputObj<TH1I> hNEvents{"hNEvents"};
  // OutputObj<TH1I> hNjets{"hNjets"};
  // OutputObj<TH1I> hNclusters{"hNclusters"};
  // OutputObj<TH1I> hNclusterConstituents{"hNclusterConstituents"};

  // // OutputObj<TH2F> hMBSpectrum{"hMBSpectrum"};
  // // OutputObj<TH2F> hClusterComparison{"hClusterComparison"};

  // OutputObj<TH2F> hJetEtaPhi{"hJetEtaPhi"};
  // OutputObj<TH2F> hJetSelectedEtaPhi{"hJetSelectedEtaPhi"};
  // OutputObj<TH1F> hJetPt{"hJetPt"};
  // OutputObj<TH2F> hClusterEtaPhi{"hClusterEtaPhi"};
  // OutputObj<TH2F> hClusterSelectedEtaPhi{"hClusterSelectedEtaPhi"};
  // OutputObj<TH1F> hClusterPt{"hClusterPt"};

  OutputObj<TH2F> hJetPtEta{"hJetPtEta"};
  OutputObj<TH2F> hJetPtPhi{"hJetPtPhi"};
  OutputObj<TH2F> hJetSelectedPtEta{"hJetSelectedPtEta"};
  OutputObj<TH2F> hJetSelectedPtPhi{"hJetSelectedPtPhi"};

  // HistogramRegistry spectra{"spectra", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  Configurable<float> f_jetPtMin{"f_jetPtMin", 0.0, "minimum jet pT cut"};
  Configurable<float> f_SD_zCut{"f_SD_zCut", 0.1, "soft drop z cut"};
  Configurable<float> f_SD_beta{"f_SD_beta", 0.0, "soft drop beta"};
  Configurable<float> f_jetR{"f_jetR", 0.4, "jet resolution parameter"}; //possible to get configurable from another task? jetR. Or use jet.r()?

  Configurable<std::vector<float>> f_ang_kappa{"f_ang_kappa", {1.0, 1.0, 2.0}, "angularity momentum exponent"}; // This is usually 1.0
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

    // Histograms
    // hNEvents.setObject( new TH1I("hNEvents", "Number of events", 1, 0, 1));
    // hNjets.setObject( new TH1I("hNjets", "Number of jets per event", 100, 0, 100));
    // hNclusters.setObject( new TH1I("hNclusters", "Number of clusters per event", 100, 0, 100));
    // hNclusterConstituents.setObject( new TH1I("hNclusterConstituents", "Number of cluster constituents per event", 100, 0, 100));

    // hMBSpectrum.setObject(
    //   new TH2F("hMBSpectrum",
    //            "MB jet/cluster spectra;#it{p}_{T}^{jet} (GeV);#it{p}_{T}^{cl} (GeV)",
    //            nPtBins, kMinPt, kMaxPt, nPtBins, kMinPt, kMaxPt/2
    //            ));
    // hClusterComparison.setObject(
    //   new TH2F("hClusterComparison",
    //            "MB jet cluster/event cluster spectra;#it{p}_{T}^{jet cl} (GeV);#it{p}_{T}^{evt cl} (GeV)",
    //            nPtBins, kMinPt, kMaxPt/2, nPtBins, kMinPt, kMaxPt/2
    //            ));

    // hJetEtaPhi.setObject(
    //   new TH2F("hJetEtaPhi",
    //            "MB jet eta phi;#eta;#phi",
    //            nEtaBins, kMinEta, kMaxEta, nPhiBins, kMinPhi, kMaxPhi
    //            ));
    // hJetSelectedEtaPhi.setObject(
    //   new TH2F("hJetSelectedEtaPhi",
    //            "Selected jet eta phi;#eta;#phi",
    //            nEtaBins, kMinEta, kMaxEta, nPhiBins, kMinPhi, kMaxPhi
    //            ));
    // hJetPt.setObject(
    //   new TH1F("hJetPt",
    //            "All jet #it{p}_{T};#it{p}_{T}",
    //            nPtBins, kMinPt, kMaxPt
    //            ));

    // hClusterEtaPhi.setObject(
    //   new TH2F("hClusterEtaPhi",
    //            "MB cluster eta phi;#eta;#phi",
    //            nEtaBins, kMinEta, kMaxEta, nPhiBins, kMinPhi, kMaxPhi
    //            ));
    // hClusterSelectedEtaPhi.setObject(
    //   new TH2F("hClusterSelectedEtaPhi",
    //            "MB cluster eta phi;#eta;#phi",
    //            nEtaBins, kMinEta, kMaxEta, nPhiBins, kMinPhi, kMaxPhi
    //            ));
    // hClusterPt.setObject(
    //   new TH1F("hClusterPt",
    //            "All cluster #it{E}_{T};#it{E}_{T}",
    //            nPtBins, kMinPt, kMaxPt/2
    //            ));

    hJetPtEta.setObject( new TH2F("hJetPtEta", "Jet #it{p}_{T} and #eta;#it{p}_{T};#eta",
                                  nPtBins, kMinPt, kMaxPt, nEtaBins, kMinEta, kMaxEta));
    hJetPtPhi.setObject( new TH2F("hJetPtPhi", "Jet #it{p}_{T} and #phi;#it{p}_{T};#phi",
                                  nPtBins, kMinPt, kMaxPt, nPhiBins, kMinPhi, kMaxPhi));
    hJetSelectedPtEta.setObject( new TH2F("hJetSelectedPtEta", "Selected jet #it{p}_{T} and #eta;#it{p}_{T};#eta",
                                  nPtBins, kMinPt, kMaxPt, nEtaBins, kMinEta, kMaxEta));
    hJetSelectedPtPhi.setObject( new TH2F("hJetSelectedPtPhi", "Selected jet #it{p}_{T} and #phi;#it{p}_{T};#phi",
                                  nPtBins, kMinPt, kMaxPt, nPhiBins, kMinPhi, kMaxPhi));
    /*
    spectra.add("hJetZg", "Jet zg", HistType::kTH2F, {
      {fgkNPtBins, kMinPt, kMaxPt, "#it{p}_{T}^{jet} (GeV)"},
      {10, 0., .5, "#it{z}_{g}"}
      });
    spectra.add("hJetRg", "Jet Rg", HistType::kTH2F, {
      {fgkNPtBins, kMinPt, kMaxPt, "#it{p}_{T}^{jet} (GeV)"},
      {20, 0., 1., "#it{R}_{g}"}
      });
    spectra.add("hJetnsd", "Jet nsd", HistType::kTH2F, {
      {fgkNPtBins, kMinPt, kMaxPt, "#it{p}_{T}^{jet} (GeV)"},
      {10, 0., 10., "#it{n}_{SD}"}
      });
    spectra.add("hJetMass", "Jet mass", HistType::kTH2F, {
      {fgkNPtBins, kMinPt, kMaxPt, "#it{p}_{T}^{jet} (GeV)"},
      {50, 0., 50., "#it{m} (GeV)"}
      });
    spectra.add("hJetNEF", "Jet NEF", HistType::kTH2F, {
      {fgkNPtBins, kMinPt, kMaxPt, "#it{p}_{T}^{jet} (GeV)"},
      {50, 0., 1., "NEF"}
      });
    spectra.add("hJetnconst", "Jet nconst", HistType::kTH2F, {
      {fgkNPtBins, kMinPt, kMaxPt, "#it{p}_{T}^{jet} (GeV)"},
      {50, 0., 50., "#it{n}_{const}"}
      });
    spectra.add("hJetPtD", "Jet ptD", HistType::kTH2F, {
      {fgkNPtBins, kMinPt, kMaxPt, "#it{p}_{T}^{jet} (GeV)"},
      {50, 0., 1., "#it{p}_{T}^{D}"}
      });
    spectra.add("hJetScaledMass", "Jet scaled mass", HistType::kTH2F, {
      {fgkNPtBins, kMinPt, kMaxPt, "#it{p}_{T}^{jet} (GeV)"},
      {200, 0., 2., "#it{m}/#it{p}_{T}"}
      });
    spectra.add("hJetChargeFrag", "Jet Charge Fragmentation", HistType::kTH2F, {
      {fgkNPtBins, kMinPt, kMaxPt, "#it{p}_{T}^{jet} (GeV)"},
      {100, 0., 1., "#it{F}^{jet}"}
      });
    spectra.add("hJetAngularities_11", "Jet Angularities 1,1: girth", HistType::kTH2F, {
      {fgkNPtBins, kMinPt, kMaxPt, "#it{p}_{T}^{jet} (GeV)"},
      {20, 0., 1., "#it{z}^{1}#theta^{1}"}
      });
    spectra.add("hJetAngularities_21", "Jet Angularities 2,1", HistType::kTH2F, {
      {fgkNPtBins, kMinPt, kMaxPt, "#it{p}_{T}^{jet} (GeV)"},
      {20, 0., 1., "#it{z}^{2}#theta{1}"}
      });
    spectra.add("hJetAngularities_12", "Jet Angularities 1,2", HistType::kTH2F, {
      {fgkNPtBins, kMinPt, kMaxPt, "#it{p}_{T}^{jet} (GeV)"},
      {20, 0., 1., "#it{z}^{1}#theta{2}"}
      });
    spectra.add("hClusterE", "Cluster E", HistType::kTH2F, {
      {fgkNPtBins, kMinPt, kMaxPt, "#it{p}_{T}^{jet} (GeV)"},
      {fgkNPtBins, kMinPt, kMaxPt, "#it{E}"}
      });
    spectra.add("hTrackPt", "Track pt", HistType::kTH2F, {
      {fgkNPtBins, kMinPt, kMaxPt, "#it{p}_{T}^{jet} (GeV)"},
      {fgkNPtBins, kMinPt, kMaxPt, "#it{p}_{T}^{track} (GeV)"}
      });
    spectra.add("hMBSpectrum", "MB jet/cluster spectra", HistType::kTH2F, {
      {fgkNPtBins, kMinPt, kMaxPt, "#it{p}_{T}^{jet} (GeV)"},
      {fgkNPtBins, kMinPt, kMaxPt/2, "#it{E} (GeV)"}
      });
    */
  }

  void SoftDrop(aod::Jet const& jet,
                std::vector<fastjet::PseudoJet> const& jetReclustered,
                float zg = 0,
                float rg = 0,
                int nsd = 0);
  void NEF(aod::Jet const& jet, std::vector<fastjet::PseudoJet> const& jetClusterConstituents);
  void ChargeFragmentation(aod::Jet const& jet,
                           std::vector<fastjet::PseudoJet> const& jetConstituents);
  float Angularity(aod::Jet const& jet,
                   std::vector<fastjet::PseudoJet> const& jetConstituents,
                   std::vector<fastjet::PseudoJet> const& jetClusterConstituents,
                   float alpha = 1.0,
                   float kappa = 1.0);
  void ptD(aod::Jet const& jet,
           std::vector<fastjet::PseudoJet> const& jetConstituents,
           std::vector<fastjet::PseudoJet> const& jetClusterConstituents);

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgVertexCut;
  // Filters for EMCAL clusters
  o2::aod::EMCALClusterDefinition clusDef = o2::aod::emcalcluster::getClusterDefinitionFromString(mClusterDefinition.value);
  Filter clusterDefinitionSelection = o2::aod::emcalcluster::definition == static_cast<int>(clusDef);

  void process(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& collision,
               aod::Jets const& jets,
               selectedClusters const& clusters,
               aod::JetClusterConstituents const& clusterConstituents,
               aod::FullJetFilters const& fullJetFilters)
  {

    // double L0emcal = 2.5, L1emcal = 19; // EMCAL hardware JE trigger thresholds
    // double maxJetPt = -1., maxClusterPt = -1., maxJetClusterPt = -1.;
    int nJets = 0, nClusters = 0, nClusterConstituents = 0;

    // hNEvents->Fill(0.5); // All events, not just ones with jets and clusters

    // for (const auto& jet : jets){
    //   hJetPt->Fill(jet.pt());
    //   hJetEtaPhi->Fill(jet.eta(), jet.phi());
    //   nJets++;
    //   if (jet.pt() > maxJetPt){
    //     maxJetPt = jet.pt();
    //     for (const auto& clusterConstituent : clusterConstituents){
    //       const auto& jetCluster = clusterConstituent.cluster();
    //       double jetClusterPt = jetCluster.energy() / std::cosh(jetCluster.eta());
    //       if (jetClusterPt > maxJetClusterPt) maxJetClusterPt = jetClusterPt;
    //     }
    //   }
    // }
    for (const auto& obj : fullJetFilters){
      if (obj.hasJetFullHighPt()){
        // nJets++;
        LOG(info) << "Found accepted event";
      }
    }
    // for (const auto& cluster : clusters){
    //   double clusterPt = cluster.energy() / std::cosh(cluster.eta());
    //   hClusterPt->Fill(clusterPt);
    //   hClusterEtaPhi->Fill(cluster.eta(), cluster.phi());
    //   nClusters++;
    //   if (clusterPt > maxClusterPt) maxClusterPt = clusterPt;
    // }
    // hNjets->Fill(nJets);
    // hNclusters->Fill(nClusters);
    // if (maxJetPt >= 0 && maxClusterPt >= 0){
    //   hMBSpectrum->Fill(maxJetPt, maxClusterPt);
    //   for (const auto& jet : jets){
    //     hJetSelectedEtaPhi->Fill(jet.eta(), jet.phi());
    //   }
    //   for (const auto& cluster : clusters){
    //     hClusterSelectedEtaPhi->Fill(cluster.eta(), cluster.phi());
    //   }
    // }
    // if (maxClusterPt >= 0 && maxJetClusterPt >= 0) hClusterComparison->Fill(maxJetClusterPt, maxClusterPt);

    /*
    spectra.fill(HIST("hJetMass"), jet.pt(), jet.mass());
    spectra.fill(HIST("hJetScaledMass"), jet.pt(), jet.mass() / jet.pt());

    jetConstituents.clear();
    jetClusterConstituents.clear();
    jetReclustered.clear();
    if (b_DoConstSub) {
      for (const auto& constituent : constituentsSub) {
        fillConstituents(constituent, jetConstituents);
      }
    } else {
      for (const auto& constituentIndex : constituents) {
        auto constituent = constituentIndex.track();
        fillConstituents(constituent, jetConstituents);
      }
      int count = 0;
      for (const auto& clusterConstituent : clusterConstituents) {
        // This is based on https://github.com/AliceO2Group/O2Physics/blob/90b01104988b5697ac108e51ea4d60429d2e9e40/PWGJE/TableProducer/jetfinder.cxx#L231
        count++;
        auto cluster = clusterConstituent.cluster();
        double pt = cluster.energy() / std::cosh(cluster.eta());
        jetClusterConstituents.emplace_back(
          pt * std::cos(cluster.phi()),
          pt * std::sin(cluster.phi()),
          pt * std::sinh(cluster.eta()),
          cluster.energy()
        );
      }
    }

    spectra.fill(HIST("hJetnconst"), jet.pt(), jetConstituents.size() + jetClusterConstituents.size());
    ChargeFragmentation(jet, jetConstituents);
    NEF(jet, jetClusterConstituents);
    ptD(jet, jetConstituents, jetClusterConstituents);

    if (f_ang_alpha->size() != f_ang_kappa->size()){
      std::cout << "Warning: different amount of alpha (" << f_ang_alpha->size() << "), kappa (" << f_ang_kappa->size() << ") values. Truncating to shortest array." << std::endl;
    }
    float ang = -999.;
    ang = Angularity(jet, jetConstituents,
                      jetClusterConstituents,
                      f_ang_alpha->at(0), f_ang_kappa->at(0));
    spectra.fill(HIST("hJetAngularities_11"), jet.pt(), ang);
    ang = Angularity(jet, jetConstituents,
                      jetClusterConstituents,
                      f_ang_alpha->at(1), f_ang_kappa->at(1));
    spectra.fill(HIST("hJetAngularities_12"), jet.pt(), ang);
    ang = Angularity(jet, jetConstituents,
                      jetClusterConstituents,
                      f_ang_alpha->at(2), f_ang_kappa->at(2));
    spectra.fill(HIST("hJetAngularities_21"), jet.pt(), ang);

    std::vector<fastjet::PseudoJet> allJetConstituents;
    allJetConstituents.reserve(jetConstituents.size() + jetClusterConstituents.size());
    allJetConstituents.insert(allJetConstituents.end(), jetConstituents.begin(), jetConstituents.end());
    allJetConstituents.insert(allJetConstituents.end(), jetClusterConstituents.begin(), jetClusterConstituents.end());
    fastjet::ClusterSequenceArea clusterSeq(jetReclusterer.findJets(allJetConstituents, jetReclustered)); // TODO: Include cluster constituents here
    jetReclustered = sorted_by_pt(jetReclustered);
    */
    // int nsd = 0.0;
    // float zg = -1.0;
    // float rg = -1.0;
    // SoftDrop(jet, jetReclustered, zg, rg, nsd);

    // jetTriggerQA(zg, rg, nsd);
  }
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<JetTriggerQA>(cfgc, TaskName{"jet-full-trigger-qa"})};
}

#endif
