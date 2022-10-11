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

  // bool isCellMasked(int towerID)
  // {
  //   bool masked = false;
  //   if (mBadChannels) {
  //     auto maskStatus = mBadChannels->getChannelStatus(towerID);
  //     masked = (maskStatus != o2::emcal::BadChannelMap::MaskType_t::GOOD_CELL);
  //   }
  //   return masked;
  // }

  // Declare filters
  o2::aod::EMCALClusterDefinition clusDef = o2::aod::emcalcluster::getClusterDefinitionFromString(mClusterDefinition.value);
  Filter clusterDefinitionSelection = o2::aod::emcalcluster::definition == static_cast<int>(clusDef);

  // std::vector<o2::emcal::Cell> mEmcalCells;

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgVertexCut;
  // TODO: Should be fiducial jet cut. I.e.: eta_jet < eta_max - R
  // Filter jetEtaFilter = (nabs(aod::jet::eta) < cfgJetEtaCut);
  // Filter jetPhiFilter = ((aod::jet::phi < cfgJetPhiMax) && (aod::jet::phi > cfgJetPhiMin));
  // Filter trackFilter = (nabs(aod::track::eta) < cfgTrackEtaCut) && (aod::track::isGlobalTrack == (uint8_t) true) && (aod::track::pt > cfgTrackLowPtCut);

  //using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection>>;

  // void process(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& collision,
  //              aod::Jets const& jets,
  //              aod::JetClusterConstituents const& clusterConstituents,
  //              selectedClusters const& clusters
  //             )
  // {

    // In the backend, the process function matches the tables via their shared identifiers. Here: the bcID
    // This means we only get events that also have cells
    // We can circumvent this by using two process functions: one that knows about the cells and one that doesn't
  void process(//aod::BC const& bc,
               soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& collision,
               aod::BCs const& bcs,
               aod::Jets const& jets,
               selectedClusters const& clusters,
               aod::JetClusterConstituents const& clusterConstituents,
               o2::aod::Calos const& allCells)
                // aod::Jet const& jet, // Should not loop over jets but over events
              //  aod::Tracks const& tracks,
              //  aod::JetTrackConstituents const& constituents,
              //  aod::JetClusterConstituents const& clusterConstituents,
              //  aod::EMCALClusters const& emcalClusters,
              //  aod::JetConstituentsSub const& constituentsSub,
              //  aod::JetFilters const& jetFilter)
  {
    bool keepEvent[kCategories]{false};

    double L0emcal = 2.5, L1emcal = 19; // EMCAL hardware JE trigger thresholds

    // for colliosons in vc
    hNEvents->Fill(0.5); // All events, not just ones with jets and clusters
    // if (collisions.size() > 1){
    //   LOG(error) << "More than one collision in the bc. This is not supported.";
    // }
    // else{
    //   for (const auto& coll : collisions){
        double maxJetPt = -1., maxClusterPt = -1., maxJetClusterPt = -1.;//, maxCellAmp = -1.;
        int nJets = 0, nClusters = 0, nClusterConstituents = 0;
        double maxCellAmp = -1.;
        // if (cell.bc() != coll.bc()) continue;
        for (const auto& cell : allCells){
          if (cell.bc() != collision.bc()) continue;
          if (cell.caloType() != 1) continue; // Check if EMCAL cell
          // if (isCellMasked(cell.cellNumber())) continue;
          // TODO: isCellMasked?
          hCellAmp->Fill(cell.amplitude());
          if (cell.amplitude() < minCellAmplitude) continue;
          if (cell.amplitude() > maxCellAmp){
            maxCellAmp = cell.amplitude();
          }
        }
        if (maxCellAmp > minCellAmplitude){
          hCellMaxAmp->Fill(maxCellAmp);
          hNEmcEvents->Fill(0.5);
        }
        else return;

        for (const auto& jet : jets){
          hJetPt->Fill(jet.pt());
          hJetEtaPhi->Fill(jet.eta(), jet.phi());
          nJets++;
          if (jet.pt() > maxJetPt){
            maxJetPt = jet.pt();
            for (const auto& clusterConstituent : clusterConstituents){
              const auto& jetCluster = clusterConstituent.cluster();
              double jetClusterPt = jetCluster.energy() / std::cosh(jetCluster.eta());
              if (jetClusterPt > maxJetClusterPt) maxJetClusterPt = jetClusterPt;
            }
          }
        }
        for (const auto& cluster : clusters){
          double clusterPt = cluster.energy() / std::cosh(cluster.eta());
          hClusterPt->Fill(clusterPt);
          hClusterEtaPhi->Fill(cluster.eta(), cluster.phi());
          nClusters++;
          if (clusterPt > maxClusterPt) maxClusterPt = clusterPt;
        }
        hNjets->Fill(nJets);
        hNclusters->Fill(nClusters);
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
        // TODO:
        // Check if EMCAL trigger is fired. Skip otherwise.
        // -> Filter on collision table?
        // Construct patch proxy from cells
        // -> Maybe do this in a different task?
        // Match cluster with cells
        // Set flags
        // collision process loop



        // if (collision.alias()[kINT7]){
        //   LOG(debug) << "Event is not Minimum Bias. Skipping";
        //   return;
        // }
        // if ( !(collision.alias()[kEJ1] || collision.alias()[kEJ2]) ){
        // if ( !(collision.alias()[kEG1] || collision.alias()[kEG2]) ){
        //   LOG(debug) << "Event did not trigger EMCAL gamma trigger. Skipping";
        //   tags(keepEvent[0], keepEvent[1]); // Minimum Bias event
        //   return;
        // }

        // for (const auto& jet : jets){
        //   spectra.fill(HIST("fJet"), 0., jet.phi(), jet.eta());
        //   spectra.fill(HIST("fJetMB"), jet.pt());
        // }
        // for (const auto& cluster : clusters){
        //   spectra.fill(HIST("fCluster"), cluster.energy(), cluster.phi(), cluster.eta());
        //   spectra.fill(HIST("fClusterMB"), cluster.energy());
        // }

        // TODO: Use jet cluster constituents
        // for (const auto& jet : jets){
        //   for (const auto& clusterConstituent : clusterConstituents){
        //     const auto& cluster = clusterConstituent.cluster();
        //     double e = cluster.energy();
        //     std::cout << "E: " << e << std::endl;
        //   }
        // }

        // TODO: check if cluster is above jet patch threshold. Then it should always be found
        //  if clus.E > 19 GeV, it must fire the jet patch
        // for (const auto& cluster : clusters){
        //   if (cluster.energy() > selectionHighECluster){
        //     keepEvent[kPatch] = true;
        //     spectra.fill(HIST("fClusterESelected"), cluster.energy());
        //     break;
        //   }
        // }

        // for (const auto& jet : jets){
        //   for (const auto& cluster : clusters){
        //     double dist = TMath::Sqrt(TMath::Power(cluster.phi() - jet.phi(), 2) + TMath::Power(cluster.eta() - jet.eta(), 2));
        //     if (dist > jet.r()) continue;
        //     if (cluster.energy() > selectionHighECluster){
        //       keepEvent[kPatch] = true;
        //       spectra.fill(HIST("fClusterDistSelected"), cluster.energy(), dist);
        //       break;
        //     }
        //   }
        // }

        // for (const auto& jet : jets){
        //   if (jet.pt() < selectionJetPt) continue;
        //   if (TMath::Abs(jet.eta() > cfgJetEtaCut || jet.phi() < cfgJetPhiMin || jet.phi() > cfgJetPhiMax)) continue;
        //   keepEvent[kMatchedJet] = true;
        //   spectra.fill(HIST("fJetPtSelected"), jet.pt());
        // }

        // for (const auto& jet : jets){
        //   for (const auto& clusterConstituent : clusterConstituents){
        //     const auto& cluster = clusterConstituent.cluster();
        //     double e = cluster.energy();
        //     std::cout << "E: " << e << std::endl;
        //   }
        // }

        for (int iDecision{0}; iDecision < kCategories; iDecision++){
          if (keepEvent[iDecision]) {
            hProcessedEvents->Fill(iDecision);
            // spectra.fill(HIST("fProcessedEvents"), iDecision);
          }
        }
      // }
    // }

    // tags(keepEvent[0], keepEvent[1]);
  } // process()
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{
  return WorkflowSpec{
    adaptAnalysisTask<fullJetFilter>(cfg)};
}
