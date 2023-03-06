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

// Full Jet Filter
// Author: Gijs van Weelden

#include <cmath>
#include <string>
#include <TMath.h>

#include "ReconstructionDataFormats/Track.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/CCDB/TriggerAliases.h"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/EMCALClusters.h"
#include "PWGJE/Core/JetFinder.h"

#include "../filterTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using selectedClusters = o2::soa::Filtered<o2::aod::EMCALClusters>;
using filteredJets = o2::soa::Filtered<o2::aod::Jets>;

struct fullJetFilter {
  enum {
    kJetFullHighPt = 0,
    kGammaHighPtEMCAL,
    kGammaHighPtDCAL,
    kCategories
  };

  Produces<aod::FullJetFilters> tags;

  OutputObj<TH1I> hProcessedEvents{"hProcessedEvents"};
  OutputObj<TH2F> hEmcClusterPtEta{"hEmcClusterPtEta"};
  OutputObj<TH2F> hEmcClusterPtPhi{"hEmcClusterPtPhi"};
  OutputObj<TH2F> hSelectedClusterPtEta{"hSelectedClusterEta"};
  OutputObj<TH2F> hSelectedClusterPtPhi{"hSelectedClusterPhi"};
  OutputObj<TH2F> hEmcJetPtEta{"hEmcJetEta"};
  OutputObj<TH2F> hEmcJetPtPhi{"hEmcJetPhi"};
  OutputObj<TH2F> hSelectedJetPtEta{"hSelectedJetEta"};
  OutputObj<TH2F> hSelectedJetPtPhi{"hSelectedJetPhi"};
  OutputObj<TH2F> hSelectedGammaEMCALPtEta{"hSelectedGammaEMCALEta"};
  OutputObj<TH2F> hSelectedGammaEMCALPtPhi{"hSelectedGammaEMCALPhi"};
  OutputObj<TH2F> hSelectedGammaDCALPtEta{"hSelectedGammaDCALEta"};
  OutputObj<TH2F> hSelectedGammaDCALPtPhi{"hSelectedGammaDCALPhi"};

  // Configurables
  Configurable<float> f_jetPtMin{"f_jetPtMin", 0.0, "minimum jet pT cut"};
  Configurable<float> f_clusterPtMin{"f_clusterPtMin", 0.0, "minimum cluster pT cut"};
  Configurable<double> f_jetR{"f_jetR", 0.2, "jet R to trigger on"};
  Configurable<float> f_gammaPtMinEMCAL{"f_gammaPtMinEMCAL", 4.0, "minimum gamma pT cut in EMCAL"};
  Configurable<float> f_gammaPtMinDCAL{"f_gammaPtMinDCAL", 4.0, "minimum gamma pT cut in DCAL"};
  Configurable<float> f_PhiEmcalOrDcal{"f_PhiEmcalOrDcal", 4, "if cluster phi is less than this value, count it to be EMCAL"};

  Configurable<std::string> mClusterDefinition{"clusterDefinition", "kV3Default", "cluster definition to be selected, e.g. V3Default"};
  Configurable<bool> b_doJetTrigger{"b_doJetTrigger", true, "run the full jet trigger"};
  Configurable<bool> b_doGammaTrigger{"b_doGammaTrigger", true, "run the gamma trigger"};
  Configurable<bool> b_IgnoreEmcalFlag{"b_IgnoreEmcalFlag", false, "ignore the EMCAL live flag check"};
  Configurable<bool> b_DoFiducialCut{"b_DoFiducialCut", false, "do a fiducial cut on jets to check if they are in the emcal"};

  void init(o2::framework::InitContext&)
  {
    // Bins here
    Int_t nPtBins = 200;
    Float_t kMinPt = 0.;
    Float_t kMaxPt = 200.;
    Int_t nPhiBins = 18 * 8;
    Float_t kMinPhi = 0.;
    Float_t kMaxPhi = 2. * TMath::Pi();
    Int_t nEtaBins = 100;
    Float_t kMinEta = -1.;
    Float_t kMaxEta = 1.;

    hProcessedEvents.setObject(new TH1I("hProcessedEvents", ";;Number of filtered events", kCategories, -0.5, kCategories - 0.5));

    hEmcClusterPtEta.setObject(new TH2F("hEmcClusterPtEta", "Emc Clusters;#it{p}_{T};#eta", nPtBins, kMinPt, kMaxPt / 2, nEtaBins, kMinEta, kMaxEta));
    hEmcClusterPtPhi.setObject(new TH2F("hEmcClusterPtPhi", "Emc Clusters;#it{p}_{T};#phi", nPtBins, kMinPt, kMaxPt / 2, nPhiBins, kMinPhi, kMaxPhi));
    hSelectedClusterPtEta.setObject(new TH2F("hSelectedClusterPtEta", "Selected Clusters;#it{p}_{T};#eta", nPtBins, kMinPt, kMaxPt / 2, nEtaBins, kMinEta, kMaxEta));
    hSelectedClusterPtPhi.setObject(new TH2F("hSelectedClusterPtPhi", "Selected Clusters;#it{p}_{T};#phi", nPtBins, kMinPt, kMaxPt / 2, nPhiBins, kMinPhi, kMaxPhi));
    hEmcJetPtEta.setObject(new TH2F("hEmcJetPtEta", "Emc Jets;#it{p}_{T};#eta", nPtBins, kMinPt, kMaxPt, nEtaBins, kMinEta, kMaxEta));
    hEmcJetPtPhi.setObject(new TH2F("hEmcJetPtPhi", "Emc Jets;#it{p}_{T};#phi", nPtBins, kMinPt, kMaxPt, nPhiBins, kMinPhi, kMaxPhi));
    hSelectedJetPtEta.setObject(new TH2F("hSelectedJetPtEta", "Selected Jets;#it{p}_{T};#eta", nPtBins, kMinPt, kMaxPt, nEtaBins, kMinEta, kMaxEta));
    hSelectedJetPtPhi.setObject(new TH2F("hSelectedJetPtPhi", "Selected Jets;#it{p}_{T};#phi", nPtBins, kMinPt, kMaxPt, nPhiBins, kMinPhi, kMaxPhi));
    hSelectedGammaEMCALPtEta.setObject(new TH2F("hSelectedGammaEMCALPtEta", "Selected Gammas EMCAL;#it{p}_{T};#eta", nPtBins, kMinPt, kMaxPt / 2, nEtaBins, kMinEta, kMaxEta));
    hSelectedGammaEMCALPtPhi.setObject(new TH2F("hSelectedGammaEMCALPtPhi", "Selected Gammas EMCAL;#it{p}_{T};#phi", nPtBins, kMinPt, kMaxPt / 2, nPhiBins, kMinPhi, kMaxPhi));
    hSelectedGammaDCALPtEta.setObject(new TH2F("hSelectedGammaDCALPtEta", "Selected Gammas DCAL;#it{p}_{T};#eta", nPtBins, kMinPt, kMaxPt / 2, nEtaBins, kMinEta, kMaxEta));
    hSelectedGammaDCALPtPhi.setObject(new TH2F("hSelectedGammaDCALPtPhi", "Selected Gammas DCAL;#it{p}_{T};#phi", nPtBins, kMinPt, kMaxPt / 2, nPhiBins, kMinPhi, kMaxPhi));
  } // init()

  // Declare filters
  o2::aod::EMCALClusterDefinition clusDef = o2::aod::emcalcluster::getClusterDefinitionFromString(mClusterDefinition.value);
  Filter clusterDefinitionSelection = o2::aod::emcalcluster::definition == static_cast<int>(clusDef);
  Filter jetRadiusSelection = o2::aod::jet::r == std::round(f_jetR * 100);

  Bool_t isJetInEmcal(filteredJets::iterator const& jet)
  {
    double emcalEtaMin = -0.7, emcalEtaMax = 0.7, emcalPhiMin = 1.40, emcalPhiMax = 3.26; // Phi: 80 - 187 deg
    double R = jet.r() * 1e-2;                                                            // Jet R is saved as round(100*R)
    if ((jet.eta() >= emcalEtaMin + R) && (jet.eta() <= emcalEtaMax - R) && (jet.phi() >= emcalPhiMin + R) && (jet.phi() <= emcalPhiMax - R)) {
      return true;
    }
    return false;
  }

  Bool_t isEvtSelected(double const& jetpt)
  {
    if (jetpt > f_jetPtMin) {
      return true;
    }
    return false;
  }

  Bool_t isEvtSelectedGammaEMCAL(double const& gammapt)
  {
    if (gammapt > f_gammaPtMinEMCAL) {
      return true;
    }
    return false;
  }

  Bool_t isEvtSelectedGammaDCAL(double const& gammapt)
  {
    if (gammapt > f_gammaPtMinDCAL) {
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

  void process(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, filteredJets const& jets, selectedClusters const& clusters)
  {
    bool keepEvent[kCategories]{false};
    double maxClusterPtEMCAL = -1., maxClusterPtDCAL = -1., maxSelectedJetPt = -1.;

    if (!collision.alias()[kTVXinEMC] && !b_IgnoreEmcalFlag) {
      tags(keepEvent[0], keepEvent[1], keepEvent[2]);
      return; // Skip events where EMCAL is not live
    }

    if (b_doJetTrigger) {
      // Store pt of leading jet inside of the emcal to trigger on
      for (const auto& jet : jets) {
        hEmcJetPtEta->Fill(jet.pt(), jet.eta());
        hEmcJetPtPhi->Fill(jet.pt(), jet.phi());
        if (jet.pt() > maxSelectedJetPt) {
          if (b_DoFiducialCut) {
            if (isJetInEmcal(jet)) {
              maxSelectedJetPt = jet.pt();
            }
          } else {
            maxSelectedJetPt = jet.pt();
          }
        }
      }
      if (isEvtSelected(maxSelectedJetPt)) {
        keepEvent[kJetFullHighPt] = true;
        for (const auto& jet : jets) {
          if (isJetInEmcal(jet)) {
            hSelectedJetPtEta->Fill(jet.pt(), jet.eta());
            hSelectedJetPtPhi->Fill(jet.pt(), jet.phi());
          }
        }
        for (const auto& cluster : clusters) {
          double clusterPt = cluster.energy() / std::cosh(cluster.eta());
          hSelectedClusterPtEta->Fill(clusterPt, cluster.eta());
          hSelectedClusterPtPhi->Fill(clusterPt, cluster.phi());
        }
      }
    }

    if (b_doGammaTrigger) {
      for (const auto& cluster : clusters) {
        double clusterPt = cluster.energy() / std::cosh(cluster.eta());
        if (cluster.phi() < 4 && clusterPt > maxClusterPtEMCAL) {
          maxClusterPtEMCAL = clusterPt;
        } else if (cluster.phi() > 4 && clusterPt > maxClusterPtDCAL) {
          maxClusterPtDCAL = clusterPt;
        }
        hEmcClusterPtEta->Fill(clusterPt, cluster.eta());
        hEmcClusterPtPhi->Fill(clusterPt, cluster.phi());
      }
      if (isEvtSelectedGammaEMCAL(maxClusterPtEMCAL)) {
        keepEvent[kGammaHighPtEMCAL] = true;
        for (const auto& cluster : clusters) {
          double clusterPt = cluster.energy() / std::cosh(cluster.eta());
          hSelectedGammaEMCALPtEta->Fill(clusterPt, cluster.eta());
          hSelectedGammaEMCALPtPhi->Fill(clusterPt, cluster.phi());
        }
      }
      if (isEvtSelectedGammaDCAL(maxClusterPtDCAL)) {
        keepEvent[kGammaHighPtDCAL] = true;
        for (const auto& cluster : clusters) {
          double clusterPt = cluster.energy() / std::cosh(cluster.eta());
          hSelectedGammaDCALPtEta->Fill(clusterPt, cluster.eta());
          hSelectedGammaDCALPtPhi->Fill(clusterPt, cluster.phi());
        }
      }
    }

    for (int iDecision{0}; iDecision < kCategories; iDecision++) {
      if (keepEvent[iDecision]) {
        hProcessedEvents->Fill(iDecision);
      }
    }
    tags(keepEvent[0], keepEvent[1], keepEvent[2]);
  } // process()
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{
  return WorkflowSpec{
    adaptAnalysisTask<fullJetFilter>(cfg)};
}
