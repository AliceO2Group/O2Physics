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

#include "EMCALBase/Geometry.h"
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
#include "PWGJE/Core/FastJetUtilities.h"

#include "../filterTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using selectedClusters = o2::soa::Filtered<o2::aod::EMCALClusters>;
using filteredJets = o2::soa::Filtered<o2::aod::FullJets>;

struct fullJetFilter {
  enum {
    kEMCALReadout = 0,
    kJetFullHighPt,
    kJetNeutralHighPt,
    kGammaHighPtEMCAL,
    kGammaHighPtDCAL,
    kGammaLowPtEMCAL,
    kGammaLowPtDCAL,
    kCategories
  };

  enum class ThresholdType_t {
    HIGH_THRESHOLD,
    LOW_THRESHOLD
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
  OutputObj<TH2F> hSelectedGammaEMCALPtEtaLow{"hSelectedGammaEMCALEtaLow"};
  OutputObj<TH2F> hSelectedGammaEMCALPtPhiLow{"hSelectedGammaEMCALPhiLow"};
  OutputObj<TH2F> hSelectedGammaDCALPtEtaLow{"hSelectedGammaDCALEtaLow"};
  OutputObj<TH2F> hSelectedGammaDCALPtPhiLow{"hSelectedGammaDCALPhiLow"};

  // Configurables
  Configurable<float> f_jetPtMin{"f_jetPtMin", 0.0, "minimum jet pT cut"};
  Configurable<float> f_clusterPtMin{"f_clusterPtMin", 0.0, "minimum cluster pT cut"};
  Configurable<int> f_jetR{"f_jetR", 20, "jet R * 100 to trigger on"};
  Configurable<int> f_JetType{"f_JetType", 0, "Jet type used for the selection (0 - full jets, 1 - neutral jets)"};
  Configurable<int> f_ObservalbeGammaTrigger{"fObservableGammaTrigger", 0, "Observable for the gamma trigger (0 - Energy, 1 - pt)"};
  Configurable<bool> b_PublishReadoutTrigger{"b_publishReadoutTrigger", false, "Publish EMCAL readout status as trigger flag"};
  Configurable<bool> b_PublishNeutralJetTrigger{"b_publishNeutralJetTrigger", false, "Publish trigger on neutral jets"};
  Configurable<float> f_gammaPtMinEMCALHigh{"f_gammaPtMinEMCALHigh", 4.0, "minimum gamma pT cut in EMCAL high threshold"};
  Configurable<float> f_gammaPtMinEMCALLow{"f_gammaPtMinEMCALLow", 1.5, "minimum gamma pT cut in EMCAL low threshold"};
  Configurable<float> f_gammaPtMinDCALHigh{"f_gammaPtMinDCALHigh", 4.0, "minimum gamma pT cut in DCAL high threshold"};
  Configurable<float> f_gammaPtMinDCALLow{"f_gammaPtMinDCALLow", 1.5, "minimum gamma pT cut in DCAL low threshold"};
  Configurable<float> f_minClusterTime{"f_minClusterTime", -999, "Min. cluster time for gamma trigger (ns)"};
  Configurable<float> f_maxClusterTime{"f_maxClusterTime", 999, "Max. cluster time for gamma trigger (ns)"};
  Configurable<float> f_PhiEmcalOrDcal{"f_PhiEmcalOrDcal", 4, "if cluster phi is less than this value, count it to be EMCAL"};

  Configurable<std::string> mClusterDefinition{"clusterDefinition", "kV3Default", "cluster definition to be selected, e.g. V3Default"};
  Configurable<bool> b_doJetTrigger{"b_doJetTrigger", true, "run the full jet trigger"};
  Configurable<bool> b_doGammaTrigger{"b_doGammaTrigger", true, "run the gamma trigger"};
  Configurable<bool> b_IgnoreEmcalFlag{"b_IgnoreEmcalFlag", false, "ignore the EMCAL live flag check"};
  Configurable<bool> b_DoFiducialCut{"b_DoFiducialCut", false, "do a fiducial cut on jets to check if they are in the emcal"};
  Configurable<bool> b_RejectExoticClusters{"b_RejectExoticClusters", true, "Reject exotic clusters"};

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

    hEmcClusterPtEta.setObject(new TH2F("hEmcClusterPtEta", Form("Emc Clusters;%s;#eta", (f_ObservalbeGammaTrigger == 0) ? "E (GeV)" : "#it{p}_{T}"), nPtBins, kMinPt, kMaxPt / 2, nEtaBins, kMinEta, kMaxEta));
    hEmcClusterPtPhi.setObject(new TH2F("hEmcClusterPtPhi", Form("Emc Clusters;%s;#phi", (f_ObservalbeGammaTrigger == 0) ? "E (GeV)" : "#it{p}_{T}"), nPtBins, kMinPt, kMaxPt / 2, nPhiBins, kMinPhi, kMaxPhi));
    hSelectedClusterPtEta.setObject(new TH2F("hSelectedClusterPtEta", Form("Selected Clusters;%s;#eta", (f_ObservalbeGammaTrigger == 0) ? "E (GeV)" : "#it{p}_{T}"), nPtBins, kMinPt, kMaxPt / 2, nEtaBins, kMinEta, kMaxEta));
    hSelectedClusterPtPhi.setObject(new TH2F("hSelectedClusterPtPhi", Form("Selected Clusters;%s;#phi", (f_ObservalbeGammaTrigger == 0) ? "E (GeV)" : "#it{p}_{T}"), nPtBins, kMinPt, kMaxPt / 2, nPhiBins, kMinPhi, kMaxPhi));
    hEmcJetPtEta.setObject(new TH2F("hEmcJetPtEta", "Emc Jets;#it{p}_{T};#eta", nPtBins, kMinPt, kMaxPt, nEtaBins, kMinEta, kMaxEta));
    hEmcJetPtPhi.setObject(new TH2F("hEmcJetPtPhi", "Emc Jets;#it{p}_{T};#phi", nPtBins, kMinPt, kMaxPt, nPhiBins, kMinPhi, kMaxPhi));
    hSelectedJetPtEta.setObject(new TH2F("hSelectedJetPtEta", "Selected Jets;#it{p}_{T};#eta", nPtBins, kMinPt, kMaxPt, nEtaBins, kMinEta, kMaxEta));
    hSelectedJetPtPhi.setObject(new TH2F("hSelectedJetPtPhi", "Selected Jets;#it{p}_{T};#phi", nPtBins, kMinPt, kMaxPt, nPhiBins, kMinPhi, kMaxPhi));
    hSelectedGammaEMCALPtEta.setObject(new TH2F("hSelectedGammaEMCALPtEta", Form("Selected Gammas EMCAL;%s;#eta", (f_ObservalbeGammaTrigger == 0) ? "E (GeV)" : "#it{p}_{T}"), nPtBins, kMinPt, kMaxPt / 2, nEtaBins, kMinEta, kMaxEta));
    hSelectedGammaEMCALPtPhi.setObject(new TH2F("hSelectedGammaEMCALPtPhi", Form("Selected Gammas EMCAL;%s;#phi", (f_ObservalbeGammaTrigger == 0) ? "E (GeV)" : "#it{p}_{T}"), nPtBins, kMinPt, kMaxPt / 2, nPhiBins, kMinPhi, kMaxPhi));
    hSelectedGammaDCALPtEta.setObject(new TH2F("hSelectedGammaDCALPtEta", Form("Selected Gammas DCAL;%s;#eta", (f_ObservalbeGammaTrigger == 0) ? "E (GeV)" : "#it{p}_{T}"), nPtBins, kMinPt, kMaxPt / 2, nEtaBins, kMinEta, kMaxEta));
    hSelectedGammaDCALPtPhi.setObject(new TH2F("hSelectedGammaDCALPtPhi", Form("Selected Gammas DCAL;%s;#phi", (f_ObservalbeGammaTrigger == 0) ? "E (GeV)" : "#it{p}_{T}"), nPtBins, kMinPt, kMaxPt / 2, nPhiBins, kMinPhi, kMaxPhi));
    hSelectedGammaEMCALPtEtaLow.setObject(new TH2F("hSelectedGammaEMCALPtEtaLow", Form("Selected Gammas EMCAL (low threshold);%s;#eta", (f_ObservalbeGammaTrigger == 0) ? "E (GeV)" : "#it{p}_{T}"), nPtBins, kMinPt, kMaxPt / 2, nEtaBins, kMinEta, kMaxEta));
    hSelectedGammaEMCALPtPhiLow.setObject(new TH2F("hSelectedGammaEMCALPtPhiLow", Form("Selected Gammas EMCAL (low threshold);%s;#phi", (f_ObservalbeGammaTrigger == 0) ? "E (GeV)" : "#it{p}_{T}"), nPtBins, kMinPt, kMaxPt / 2, nPhiBins, kMinPhi, kMaxPhi));
    hSelectedGammaDCALPtEtaLow.setObject(new TH2F("hSelectedGammaDCALPtEtaLow", Form("Selected Gammas DCAL (low threshold);%s;#eta", (f_ObservalbeGammaTrigger == 0) ? "E (GeV)" : "#it{p}_{T}"), nPtBins, kMinPt, kMaxPt / 2, nEtaBins, kMinEta, kMaxEta));
    hSelectedGammaDCALPtPhiLow.setObject(new TH2F("hSelectedGammaDCALPtPhiLow", Form("Selected Gammas DCAL (low threshold);%s;#phi", (f_ObservalbeGammaTrigger == 0) ? "E (GeV)" : "#it{p}_{T}"), nPtBins, kMinPt, kMaxPt / 2, nPhiBins, kMinPhi, kMaxPhi));

    LOG(info) << "Jet trigger: " << (b_doJetTrigger ? "on" : "off");
    LOG(info) << "Gamma trigger: " << (b_doJetTrigger ? "on" : "off");
    LOG(info) << "Thresholds gamma trigger: EG1 " << f_gammaPtMinEMCALHigh << " GeV, DG1 " << f_gammaPtMinDCALHigh << " GeV, EG2 " << f_gammaPtMinEMCALLow << " GeV, DG2 " << f_gammaPtMinDCALLow << " GeV";
    LOG(info) << "Gamma trigger observable: " << (f_ObservalbeGammaTrigger == 0 ? "Energy" : "pt");
    LOG(info) << "Jet trigger: Type: " << (f_JetType == 0 ? "Full jets" : "Neutral jets") << ", R: " << (static_cast<double>(f_jetR) / 100.) << ", pt > " << f_jetPtMin << " cluster(" << f_clusterPtMin << ")";
    LOG(info) << "Ignore EMCAL flag: " << (b_IgnoreEmcalFlag ? "yes" : "no");
    LOG(info) << "Publishing neutral jet trigger: " << (b_PublishNeutralJetTrigger ? "yes" : "no");
    LOG(info) << "Publishing EMCAL trigger: " << (b_PublishReadoutTrigger ? "yes" : "no");
    LOG(info) << "Cluster time: " << f_minClusterTime << " ns < t < " << f_maxClusterTime << "ns";
    LOG(info) << "Exotics rejection: " << (b_RejectExoticClusters ? "yes" : "no");
  } // init()

  // Declare filters
  o2::aod::EMCALClusterDefinition clusDef = o2::aod::emcalcluster::getClusterDefinitionFromString(mClusterDefinition.value);
  Filter clusterDefinitionSelection = o2::aod::emcalcluster::definition == static_cast<int>(clusDef);
  Filter jetRadiusSelection = o2::aod::jet::r == f_jetR;

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

  Bool_t isEvtSelectedGamma(double const& gammapt, o2::emcal::AcceptanceType_t subdet, ThresholdType_t thresholdt)
  {
    float threshold = 0;
    switch (subdet) {
      case o2::emcal::AcceptanceType_t::EMCAL_ACCEPTANCE: {
        switch (thresholdt) {
          case ThresholdType_t::HIGH_THRESHOLD: {
            threshold = f_gammaPtMinEMCALHigh;
            break;
          }
          case ThresholdType_t::LOW_THRESHOLD: {
            threshold = f_gammaPtMinEMCALLow;
            break;
          }
        }
        break;
      }
      case o2::emcal::AcceptanceType_t::DCAL_ACCEPTANCE: {
        switch (thresholdt) {
          case ThresholdType_t::HIGH_THRESHOLD: {
            threshold = f_gammaPtMinDCALHigh;
            break;
          }
          case ThresholdType_t::LOW_THRESHOLD: {
            threshold = f_gammaPtMinDCALLow;
            break;
          }
        }
        break;
      }
      case o2::emcal::AcceptanceType_t::NON_ACCEPTANCE:
        return false;
    }
    if (gammapt > threshold) {
      return true;
    }
    return false;
  }

  Bool_t isClusterInEmcal(selectedClusters::iterator const& cluster)
  {
    if (TVector2::Phi_0_2pi(cluster.phi()) < f_PhiEmcalOrDcal) {
      return true;
    }
    return false;
  }

  void process(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, filteredJets const& jets, selectedClusters const& clusters)
  {
    bool keepEvent[kCategories]{false};
    double maxClusterObservableEMCAL = -1., maxClusterObservableDCAL = -1., maxSelectedJetPt = -1.;

    if (!b_IgnoreEmcalFlag && !collision.alias()[kTVXinEMC]) {
      tags(keepEvent[0], keepEvent[1], keepEvent[2], keepEvent[3], keepEvent[4], keepEvent[5], keepEvent[6]);
      return; // Skip events where EMCAL is not live
    }

    if (b_PublishReadoutTrigger && collision.alias()[kTVXinEMC]) {
      keepEvent[kEMCALReadout] = true;
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
        keepEvent[f_JetType == 0 ? kJetFullHighPt : kJetNeutralHighPt] = true;
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
      static constexpr std::array<ThresholdType_t, 2> thresholds = {{ThresholdType_t::HIGH_THRESHOLD, ThresholdType_t::LOW_THRESHOLD}};
      static constexpr std::array<o2::emcal::AcceptanceType_t, 2> subdets = {{o2::emcal::AcceptanceType_t::EMCAL_ACCEPTANCE, o2::emcal::AcceptanceType_t::DCAL_ACCEPTANCE}};
      std::array<TH2*, 4> acceptanceHistsPtEta{{hSelectedGammaEMCALPtEta.object.get(), hSelectedGammaDCALPtEta.object.get(), hSelectedGammaEMCALPtEtaLow.object.get(), hSelectedGammaDCALPtEtaLow.object.get()}},
        acceptanceHistsPtPhi{{hSelectedGammaEMCALPtPhi.object.get(), hSelectedGammaDCALPtPhi.object.get(), hSelectedGammaEMCALPtPhiLow.object.get(), hSelectedGammaDCALPtPhiLow.object.get()}};
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
        double observableGamma = (f_ObservalbeGammaTrigger == 0) ? cluster.energy() : cluster.energy() / std::cosh(cluster.eta());
        if (TVector2::Phi_0_2pi(cluster.phi()) < 4 && observableGamma > maxClusterObservableEMCAL) {
          maxClusterObservableEMCAL = observableGamma;
        } else if (TVector2::Phi_0_2pi(cluster.phi()) > 4 && observableGamma > maxClusterObservableDCAL) {
          maxClusterObservableDCAL = observableGamma;
        }
        hEmcClusterPtEta->Fill(observableGamma, cluster.eta());
        hEmcClusterPtPhi->Fill(observableGamma, cluster.phi());
        analysedClusters.push_back({static_cast<float>(observableGamma), cluster.eta(), cluster.phi()});
      }
      for (int ithreshold = 0; ithreshold < thresholds.size(); ithreshold++) {
        for (int isubdet = 0; isubdet < subdets.size(); isubdet++) {
          if (isEvtSelectedGamma(subdets[isubdet] == o2::emcal::AcceptanceType_t::EMCAL_ACCEPTANCE ? maxClusterObservableEMCAL : maxClusterObservableDCAL, subdets[isubdet], thresholds[ithreshold])) {
            keepEvent[kGammaHighPtEMCAL + ithreshold * subdets.size() + isubdet] = true;
            for (auto& cluster : analysedClusters) {
              acceptanceHistsPtEta[ithreshold * subdets.size() + isubdet]->Fill(cluster.mTriggerObservable, cluster.mEta);
              acceptanceHistsPtPhi[ithreshold * subdets.size() + isubdet]->Fill(cluster.mTriggerObservable, cluster.mPhi);
            }
          }
        }
      }
    }

    for (int iDecision{0}; iDecision < kCategories; iDecision++) {
      if (keepEvent[iDecision]) {
        hProcessedEvents->Fill(iDecision);
      }
    }
    tags(keepEvent[0], keepEvent[1], keepEvent[2], keepEvent[3], keepEvent[4], keepEvent[5], keepEvent[6]);
  } // process()
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{
  return WorkflowSpec{
    adaptAnalysisTask<fullJetFilter>(cfg)};
}
