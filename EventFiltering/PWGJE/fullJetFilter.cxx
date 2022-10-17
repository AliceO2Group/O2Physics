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
// O2 includes

#include "ReconstructionDataFormats/Track.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/ASoAHelpers.h"
//FK #include "Common/Core/PID/PIDResponse.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/CCDB/TriggerAliases.h"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/EMCALClusters.h"
#include "PWGJE/Core/JetFinder.h"

// #include "EMCALCalib/BadChannelMap.h"

// #include "DataFormatsEMCAL/Cell.h"

// #include "Commmon/CCDB/TriggerAliases.h"
#include "/Users/gijsvanweelden/alice/O2/DataFormats/Detectors/EMCAL/include/DataFormatsEMCAL/AnalysisCluster.h"
#include "/Users/gijsvanweelden/alice/O2/DataFormats/Detectors/EMCAL/include/DataFormatsEMCAL/Cluster.h"

#include "../filterTables.h"

#include "Framework/HistogramRegistry.h"

#include <cmath>
#include <string>
#include <TMath.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using selectedClusters = o2::soa::Filtered<o2::aod::EMCALClusters>;

static const std::vector<std::string> matchInfoNames {"Patch", "MatchedJet"}; //{"MB", "Patch", "MatchedJet"};

struct fullJetFilter {
  enum { //kMinimumBias = 0,
    kPatch = 0,
    kMatchedJet,
    kCategories };

  // Produces<aod::FJetFilters> tags;
  OutputObj<TH1I> hProcessedEvents{"hProcessedEvents"};
  OutputObj<TH1I> hNEvents{"hNEvents"};
  OutputObj<TH1I> hNEmcEvents{"hNEmcEvents"};
  OutputObj<TH1I> hNjets{"hNjets"};
  OutputObj<TH1I> hNclusters{"hNclusters"};
  OutputObj<TH1I> hNclusterConstituents{"hNclusterConstituents"};

  OutputObj<TH2F> hMBSpectrum{"hMBSpectrum"};
  OutputObj<TH2F> hClusterComparison{"hClusterComparison"};

  OutputObj<TH1F> hCellAmp{"hCellAmp"};
  OutputObj<TH1F> hCellMaxAmp{"hCellMaxAmp"};
  OutputObj<TH2F> hJetEtaPhi{"hJetEtaPhi"};
  OutputObj<TH2F> hJetSelectedEtaPhi{"hJetSelectedEtaPhi"};
  OutputObj<TH1F> hJetPt{"hJetPt"};
  OutputObj<TH2F> hClusterEtaPhi{"hClusterEtaPhi"};
  OutputObj<TH2F> hClusterSelectedEtaPhi{"hClusterSelectedEtaPhi"};
  OutputObj<TH1F> hClusterPt{"hClusterPt"};

  Configurable<float> f_jetPtMin{"f_jetPtMin", 0.0, "minimum jet pT cut"};
  Configurable<float> cfgVertexCut{"cfgVertexCut", 10.0f, "Accepted z-vertex range"};
  Configurable<double> minCellAmplitude{"minCellAmplitude", 0.1, "Minimum cell amplitude for histograms."};
  Configurable<std::string> mClusterDefinition{"clusterDefinition", "kV3Default", "cluster definition to be selected, e.g. V3Default"};

  // std::shared_ptr<o2::emcal::BadChannelMap> mBadChannels;

  void init(o2::framework::InitContext&)
  {
    // Bins here
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
    hNEvents.setObject( new TH1I("hNEvents", "Number of events", 1, 0, 1));
    hNEmcEvents.setObject( new TH1I("hNEmcEvents", "Number of events where EMCAL is live", 1, 0, 1));
    hNjets.setObject( new TH1I("hNjets", "Number of jets per event", 100, 0, 100));
    hNclusters.setObject( new TH1I("hNclusters", "Number of clusters per event", 100, 0, 100));
    hNclusterConstituents.setObject( new TH1I("hNclusterConstituents", "Number of cluster constituents per event", 100, 0, 100));

    hCellAmp.setObject(
      new TH1F("hCellAmp",
               "Emc cell amp;amp",
               10*nPtBins, kMinPt, kMaxPt
               ));
    hCellMaxAmp.setObject(
      new TH1F("hCellMaxAmp",
               "Highest amp emc cell;amp",
               10*nPtBins, kMinPt, kMaxPt
               ));

    hMBSpectrum.setObject(
      new TH2F("hMBSpectrum",
               "MB jet/cluster spectra;#it{p}_{T}^{jet} (GeV);#it{p}_{T}^{cl} (GeV)",
               nPtBins, kMinPt, kMaxPt, nPtBins, kMinPt, kMaxPt/2
               ));
    hClusterComparison.setObject(
      new TH2F("hClusterComparison",
               "MB jet cluster/event cluster spectra;#it{p}_{T}^{jet cl} (GeV);#it{p}_{T}^{evt cl} (GeV)",
               nPtBins, kMinPt, kMaxPt/2, nPtBins, kMinPt, kMaxPt/2
               ));

    hJetEtaPhi.setObject(
      new TH2F("hJetEtaPhi",
               "MB jet eta phi;#eta;#phi",
               nEtaBins, kMinEta, kMaxEta, nPhiBins, kMinPhi, kMaxPhi
               ));
    hJetSelectedEtaPhi.setObject(
      new TH2F("hJetSelectedEtaPhi",
               "Selected jet eta phi;#eta;#phi",
               nEtaBins, kMinEta, kMaxEta, nPhiBins, kMinPhi, kMaxPhi
               ));
    hJetPt.setObject(
      new TH1F("hJetPt",
               "All jet #it{p}_{T};#it{p}_{T}",
               nPtBins, kMinPt, kMaxPt
               ));

    hClusterEtaPhi.setObject(
      new TH2F("hClusterEtaPhi",
               "MB cluster eta phi;#eta;#phi",
               nEtaBins, kMinEta, kMaxEta, nPhiBins, kMinPhi, kMaxPhi
               ));
    hClusterSelectedEtaPhi.setObject(
      new TH2F("hClusterSelectedEtaPhi",
               "MB cluster eta phi;#eta;#phi",
               nEtaBins, kMinEta, kMaxEta, nPhiBins, kMinPhi, kMaxPhi
               ));
    hClusterPt.setObject(
      new TH1F("hClusterPt",
               "All cluster #it{E}_{T};#it{E}_{T}",
               nPtBins, kMinPt, kMaxPt/2
               ));
  } // init()

  // Declare filters
  o2::aod::EMCALClusterDefinition clusDef = o2::aod::emcalcluster::getClusterDefinitionFromString(mClusterDefinition.value);
  Filter clusterDefinitionSelection = o2::aod::emcalcluster::definition == static_cast<int>(clusDef);

  // Filter collisionFilter = nabs(aod::collision::posZ) < cfgVertexCut;

  // In the backend, the process function matches the tables via their shared identifiers. Here: the bcID
  // This means we only get events that also have cells
  // We can circumvent this by using two process functions: one that knows about the cells and one that doesn't
  void process(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision,
               aod::Jets const& jets,
               selectedClusters const& clusters
               )
  {
    bool keepEvent[kCategories]{false};

    double L0emcal = 2.5, L1emcal = 19; // EMCAL hardware JE trigger thresholds
    double maxJetPt = -1., maxClusterPt = -1., maxJetClusterPt = -1., maxCellAmp = -1.;
    int nJets = 0, nClusters = 0, nClustersMethod2 = 0, nClusterConstituents = 0;

    // for (const auto& coll : collisions){
    hNEvents->Fill(0.5); // All events, not just ones with jets and clusters
    // }

    // /*
    for (const auto& jet : jets){
      hJetPt->Fill(jet.pt());
      hJetEtaPhi->Fill(jet.eta(), jet.phi());
      nJets++;
      if (jet.pt() > maxJetPt){
        maxJetPt = jet.pt();
        // for (const auto& clusterConstituent : clusterConstituents){
        //   const auto& jetCluster = clusterConstituent.cluster();
        //   double jetClusterPt = jetCluster.energy() / std::cosh(jetCluster.eta());
        //   if (jetClusterPt > maxJetClusterPt) maxJetClusterPt = jetClusterPt;
        // }
      }
    }
    // */
    // /*
    for (const auto& cluster : clusters){
      double clusterPt = cluster.energy() / std::cosh(cluster.eta());
      hClusterPt->Fill(clusterPt);
      hClusterEtaPhi->Fill(cluster.eta(), cluster.phi());
      nClusters++;
      if (clusterPt > maxClusterPt) maxClusterPt = clusterPt;
    }
    // */

    // /*
    hNjets->Fill(nJets);
    hNclusters->Fill(nClusters);

    // For now, use events with a cluster as nEmcEvents
    if (nClusters > 0) hNEmcEvents->Fill(0.5);
    nClustersMethod2 = clusters.size();
    if (nClustersMethod2 > 0) hNEmcEvents->Fill(1.5);

    if (maxJetPt >= 0 && maxClusterPt >= 0){
      hMBSpectrum->Fill(maxJetPt, maxClusterPt);
      for (const auto& jet : jets){
        hJetSelectedEtaPhi->Fill(jet.eta(), jet.phi());
      }
      for (const auto& cluster : clusters){
        hClusterSelectedEtaPhi->Fill(cluster.eta(), cluster.phi());
      }
    }
    if (maxClusterPt >= 0 && maxJetClusterPt >= 0) hClusterComparison->Fill(maxJetClusterPt, maxClusterPt);
    // */

    // /*
    for (int iDecision{0}; iDecision < kCategories; iDecision++){
      if (keepEvent[iDecision]) {
        hProcessedEvents->Fill(iDecision);
      }
    }
    // */
    // tags(keepEvent[0], keepEvent[1]);
  } // process()
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{
  return WorkflowSpec{
    adaptAnalysisTask<fullJetFilter>(cfg)};
}
