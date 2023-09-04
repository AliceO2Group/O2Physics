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

#include <array>
#include <cmath>
#include <string>
#include <string_view>
#include <TMath.h>

#include <boost/algorithm/string/case_conv.hpp>

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsCTP/Configuration.h"
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

struct fullJetFilter {
  using collisionInfo = soa::Join<aod::Collisions, aod::EvSels>::iterator;
  using BCsWithBcSelsRun3 = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels>;
  using selectedClusters = o2::soa::Filtered<o2::aod::EMCALClusters>;
  using filteredFullJets = o2::soa::Filtered<o2::aod::FullJets>;
  using filteredNeutralJets = o2::soa::Filtered<o2::aod::NeutralJets>;

  enum {
    kEMCALReadout = 0,
    kJetFullHighPt,
    kJetFullLowPt,
    kJetNeutralHighPt,
    kJetNeutralLowPt,
    kGammaVeryHighPtEMCAL,
    kGammaVeryHighPtDCAL,
    kGammaHighPtEMCAL,
    kGammaHighPtDCAL,
    kGammaLowPtEMCAL,
    kGammaLowPtDCAL,
    kGammaVeryLowPtEMCAL,
    kGammaVeryLowPtDCAL,
    kCategories
  };

  enum class ThresholdType_t {
    VERY_HIGH_THRESHOLD,
    HIGH_THRESHOLD,
    LOW_THRESHOLD,
    VERY_LOW_THRESHOLD
  };

  enum class EMCALHWTriggerConfiguration {
    MB_ONLY,
    EMC_TRIGGERD,
    UNKNOWN
  };

  enum EMCALHardwareTrigger {
    TRG_MB,
    TRG_EMC7,
    TRG_DMC7,
    TRG_NTriggers
  };

  Produces<aod::FullJetFilters> tags;

  OutputObj<TH1D> hProcessedEvents{"hProcessedEvents"};
  OutputObj<TH2F> hEmcClusterPtEta{"hEmcClusterPtEta"};
  OutputObj<TH2F> hEmcClusterPtPhi{"hEmcClusterPtPhi"};
  OutputObj<TH2F> hEmcClusterPtEtaMinBias{"hEmcClusterPtEtaMinBias"};
  OutputObj<TH2F> hEmcClusterPtPhiMinBias{"hEmcClusterPtPhiMinBias"};
  OutputObj<TH2F> hEmcClusterPtEtaLevel0{"hEmcClusterPtEtaLevel0"};
  OutputObj<TH2F> hEmcClusterPtPhiLevel0{"hEmcClusterPtPhiLevel0"};
  OutputObj<TH2F> hSelectedClusterPtEta{"hSelectedClusterEta"};
  OutputObj<TH2F> hSelectedClusterPtPhi{"hSelectedClusterPhi"};
  OutputObj<TH2F> hSelectedClusterPtEtaLow{"hSelectedClusterEtaLow"};
  OutputObj<TH2F> hSelectedClusterPtPhiLow{"hSelectedClusterPhiLow"};
  OutputObj<TH2F> hEmcJetPtEta{"hEmcJetEta"};
  OutputObj<TH2F> hEmcJetPtPhi{"hEmcJetPhi"};
  OutputObj<TH2F> hSelectedJetPtEta{"hSelectedJetEta"};
  OutputObj<TH2F> hSelectedJetPtPhi{"hSelectedJetPhi"};
  OutputObj<TH2F> hSelectedJetPtEtaLow{"hSelectedJetEtaLow"};
  OutputObj<TH2F> hSelectedJetPtPhiLow{"hSelectedJetPhiLow"};
  OutputObj<TH2F> hSelectedGammaEMCALPtEtaVeryHigh{"hSelectedGammaEMCALEtaVeryHigh"};
  OutputObj<TH2F> hSelectedGammaEMCALPtPhiVeryHigh{"hSelectedGammaEMCALPhiVeryHigh"};
  OutputObj<TH2F> hSelectedGammaDCALPtEtaVeryHigh{"hSelectedGammaDCALEtaVeryHigh"};
  OutputObj<TH2F> hSelectedGammaDCALPtPhiVeryHigh{"hSelectedGammaDCALPhiVeryHigh"};
  OutputObj<TH2F> hSelectedGammaEMCALPtEta{"hSelectedGammaEMCALEta"};
  OutputObj<TH2F> hSelectedGammaEMCALPtPhi{"hSelectedGammaEMCALPhi"};
  OutputObj<TH2F> hSelectedGammaDCALPtEta{"hSelectedGammaDCALEta"};
  OutputObj<TH2F> hSelectedGammaDCALPtPhi{"hSelectedGammaDCALPhi"};
  OutputObj<TH2F> hSelectedGammaEMCALPtEtaLow{"hSelectedGammaEMCALEtaLow"};
  OutputObj<TH2F> hSelectedGammaEMCALPtPhiLow{"hSelectedGammaEMCALPhiLow"};
  OutputObj<TH2F> hSelectedGammaDCALPtEtaLow{"hSelectedGammaDCALEtaLow"};
  OutputObj<TH2F> hSelectedGammaDCALPtPhiLow{"hSelectedGammaDCALPhiLow"};
  OutputObj<TH2F> hSelectedGammaEMCALPtEtaVeryLow{"hSelectedGammaEMCALEtaVeryLow"};
  OutputObj<TH2F> hSelectedGammaEMCALPtPhiVeryLow{"hSelectedGammaEMCALPhiVeryLow"};
  OutputObj<TH2F> hSelectedGammaDCALPtEtaVeryLow{"hSelectedGammaDCALEtaVeryLow"};
  OutputObj<TH2F> hSelectedGammaDCALPtPhiVeryLow{"hSelectedGammaDCALPhiVeryLow"};
  OutputObj<TH1D> hMaxJetPt{"hMaxJetPt"};
  OutputObj<TH1D> hSelectMaxJetPt{"hSelectMaxJetPt"};
  OutputObj<TH1D> hSelectMaxJetPtLow{"hSelectMaxJetPtLow"};
  OutputObj<TH1D> hMaxClusterEMCAL{"hMaxClusterEMCAL"};
  OutputObj<TH1D> hMaxClusterDCAL{"hMaxClusterDCAL"};
  OutputObj<TH1D> hMaxClusterEMCALMinBias{"hMaxClusterEMCALMinBias"};
  OutputObj<TH1D> hMaxClusterDCALMinBias{"hMaxClusterDCALMinBias"};
  OutputObj<TH1D> hMaxClusterEMCALLevel0{"hMaxClusterEMCALLevel0"};
  OutputObj<TH1D> hMaxClusterDCALLevel0{"hMaxClusterDCALLevel0"};
  OutputObj<TH1D> hSelectGammaVeryHighMaxClusterEMCAL{"hSelectGammaVeryHighMaxClusterEMCAL"};
  OutputObj<TH1D> hSelectGammaMaxClusterEMCAL{"hSelectGammaMaxClusterEMCAL"};
  OutputObj<TH1D> hSelectGammaLowMaxClusterEMCAL{"hSelectGammaLowMaxClusterEMCAL"};
  OutputObj<TH1D> hSelectGammaVeryLowMaxClusterEMCAL{"hSelectGammaVeryLowMaxClusterEMCAL"};
  OutputObj<TH1D> hSelectGammaVeryHighMaxClusterDCAL{"hSelectGammaVeryHighMaxClusterDCAL"};
  OutputObj<TH1D> hSelectGammaMaxClusterDCAL{"hSelectGammaMaxClusterDCAL"};
  OutputObj<TH1D> hSelectGammaLowMaxClusterDCAL{"hSelectGammaLowMaxClusterDCAL"};
  OutputObj<TH1D> hSelectGammaVeryLowMaxClusterDCAL{"hSelectGammaVertLowMaxClusterDCAL"};

  // Configurables
  Configurable<float> f_jetPtMin{"f_jetPtMin", 0.0, "minimum jet pT cut"};
  Configurable<float> f_jetPtMinLow{"f_jetPtMinLow", 0.0, "minimum jet pT cut (low threshold)"};
  Configurable<float> f_jetPtMinMB{"f_jetPtMinMB", 0.0, "minimum jet pT cut (MB-only runs)"};
  Configurable<float> f_jetPtMinLowMB{"f_jetPtMinLowMB", 0.0, "minimum jet pT cut (low threshold, MB-only runs)"};
  Configurable<float> f_clusterPtMin{"f_clusterPtMin", 0.0, "minimum cluster pT cut"};
  Configurable<int> f_jetR{"f_jetR", 20, "jet R * 100 to trigger on"};
  Configurable<int> f_ObservalbeGammaTrigger{"fObservableGammaTrigger", 0, "Observable for the gamma trigger (0 - Energy, 1 - pt)"};
  Configurable<bool> b_PublishReadoutTrigger{"b_publishReadoutTrigger", false, "Publish EMCAL readout status as trigger flag"};
  Configurable<bool> b_PublishNeutralJetTrigger{"b_publishNeutralJetTrigger", false, "Publish trigger on neutral jets"};
  Configurable<float> f_gammaPtMinEMCALVeryHigh{"f_gammaPtMinEMCALVeryHigh", 9.0, "minimum gamma pT cut in EMCAL very high threshold"};
  Configurable<float> f_gammaPtMinEMCALHigh{"f_gammaPtMinEMCALHigh", 4.0, "minimum gamma pT cut in EMCAL high threshold"};
  Configurable<float> f_gammaPtMinEMCALLow{"f_gammaPtMinEMCALLow", 1.5, "minimum gamma pT cut in EMCAL low threshold"};
  Configurable<float> f_gammaPtMinEMCALVeryLow{"f_gammaPtMinEMCALVeryLow", 0.5, "minimum gamma pT cut in EMCAL very low threshold"};
  Configurable<float> f_gammaPtMinDCALVeryHigh{"f_gammaPtMinDCALVeryHigh", 9.0, "minimum gamma pT cut in DCAL very high threshold"};
  Configurable<float> f_gammaPtMinDCALHigh{"f_gammaPtMinDCALHigh", 4.0, "minimum gamma pT cut in DCAL high threshold"};
  Configurable<float> f_gammaPtMinDCALLow{"f_gammaPtMinDCALLow", 1.5, "minimum gamma pT cut in DCAL low threshold"};
  Configurable<float> f_gammaPtMinDCALVeryLow{"f_gammaPtMinDCALVeryLow", 0.5, "minimum gamma pT cut in DCAL very low threshold"};
  Configurable<float> f_gammaPtMinEMCALHighMB{"f_gammaPtMinEMCALHighMB", 4.0, "minimum gamma pT cut in EMCAL high threshold (MB-only runs)"};
  Configurable<float> f_gammaPtMinEMCALLowMB{"f_gammaPtMinEMCALLowMB", 1.5, "minimum gamma pT cut in EMCAL low threshold (MB-only runs)"};
  Configurable<float> f_gammaPtMinDCALHighMB{"f_gammaPtMinDCALHighMB", 4.0, "minimum gamma pT cut in DCAL high threshold (MB-only runs)"};
  Configurable<float> f_gammaPtMinDCALLowMB{"f_gammaPtMinDCALLowMB", 1.5, "minimum gamma pT cut in DCAL low threshold (MB-only runs)"};
  Configurable<float> f_minClusterTime{"f_minClusterTime", -999, "Min. cluster time for gamma trigger (ns)"};
  Configurable<float> f_maxClusterTime{"f_maxClusterTime", 999, "Max. cluster time for gamma trigger (ns)"};
  Configurable<float> f_PhiEmcalOrDcal{"f_PhiEmcalOrDcal", 4, "if cluster phi is less than this value, count it to be EMCAL"};

  Configurable<std::string> mClusterDefinition{"clusterDefinition", "kV3Default", "cluster definition to be selected, e.g. V3Default"};
  Configurable<bool> b_doJetTrigger{"b_doJetTrigger", true, "run the full jet trigger"};
  Configurable<bool> b_doGammaTrigger{"b_doGammaTrigger", true, "run the gamma trigger"};
  Configurable<bool> b_IgnoreEmcalFlag{"b_IgnoreEmcalFlag", false, "ignore the EMCAL live flag check"};
  Configurable<bool> b_DoFiducialCut{"b_DoFiducialCut", false, "do a fiducial cut on jets to check if they are in the emcal"};
  Configurable<bool> b_RejectExoticClusters{"b_RejectExoticClusters", true, "Reject exotic clusters"};

  int lastRun = -1;
  EMCALHWTriggerConfiguration mHardwareTriggerConfig = EMCALHWTriggerConfiguration::UNKNOWN;

  Service<o2::ccdb::BasicCCDBManager> ccdb;

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

    hProcessedEvents.setObject(new TH1D("hProcessedEvents", ";;Number of filtered events", kCategories, -0.5, kCategories - 0.5));

    hEmcClusterPtEta.setObject(new TH2F("hEmcClusterPtEta", Form("Emc Clusters;%s;#eta", (f_ObservalbeGammaTrigger == 0) ? "E (GeV)" : "#it{p}_{T}"), nPtBins, kMinPt, kMaxPt / 2, nEtaBins, kMinEta, kMaxEta));
    hEmcClusterPtPhi.setObject(new TH2F("hEmcClusterPtPhi", Form("Emc Clusters;%s;#phi", (f_ObservalbeGammaTrigger == 0) ? "E (GeV)" : "#it{p}_{T}"), nPtBins, kMinPt, kMaxPt / 2, nPhiBins, kMinPhi, kMaxPhi));
    hEmcClusterPtEtaMinBias.setObject(new TH2F("hEmcClusterPtEtaMinBias", Form("Emc Clusters (MB);%s;#eta", (f_ObservalbeGammaTrigger == 0) ? "E (GeV)" : "#it{p}_{T}"), nPtBins, kMinPt, kMaxPt / 2, nEtaBins, kMinEta, kMaxEta));
    hEmcClusterPtPhiMinBias.setObject(new TH2F("hEmcClusterPtPhiMinBias", Form("Emc Clusters (MB);%s;#phi", (f_ObservalbeGammaTrigger == 0) ? "E (GeV)" : "#it{p}_{T}"), nPtBins, kMinPt, kMaxPt / 2, nPhiBins, kMinPhi, kMaxPhi));
    hEmcClusterPtEtaLevel0.setObject(new TH2F("hEmcClusterPtEtaLevel0", Form("Emc Clusters (L0);%s;#eta", (f_ObservalbeGammaTrigger == 0) ? "E (GeV)" : "#it{p}_{T}"), nPtBins, kMinPt, kMaxPt / 2, nEtaBins, kMinEta, kMaxEta));
    hEmcClusterPtPhiLevel0.setObject(new TH2F("hEmcClusterPtPhiLevel0", Form("Emc Clusters (L0);%s;#phi", (f_ObservalbeGammaTrigger == 0) ? "E (GeV)" : "#it{p}_{T}"), nPtBins, kMinPt, kMaxPt / 2, nPhiBins, kMinPhi, kMaxPhi));
    hSelectedClusterPtEta.setObject(new TH2F("hSelectedClusterPtEta", Form("Selected Clusters;%s;#eta", (f_ObservalbeGammaTrigger == 0) ? "E (GeV)" : "#it{p}_{T}"), nPtBins, kMinPt, kMaxPt / 2, nEtaBins, kMinEta, kMaxEta));
    hSelectedClusterPtPhi.setObject(new TH2F("hSelectedClusterPtPhi", Form("Selected Clusters;%s;#phi", (f_ObservalbeGammaTrigger == 0) ? "E (GeV)" : "#it{p}_{T}"), nPtBins, kMinPt, kMaxPt / 2, nPhiBins, kMinPhi, kMaxPhi));
    hSelectedClusterPtEtaLow.setObject(new TH2F("hSelectedClusterPtEtaLow", Form("Selected Clusters (low threshold);%s;#eta", (f_ObservalbeGammaTrigger == 0) ? "E (GeV)" : "#it{p}_{T}"), nPtBins, kMinPt, kMaxPt / 2, nEtaBins, kMinEta, kMaxEta));
    hSelectedClusterPtPhiLow.setObject(new TH2F("hSelectedClusterPtPhiLow", Form("Selected Clusters (low threshold);%s;#phi", (f_ObservalbeGammaTrigger == 0) ? "E (GeV)" : "#it{p}_{T}"), nPtBins, kMinPt, kMaxPt / 2, nPhiBins, kMinPhi, kMaxPhi));
    hEmcJetPtEta.setObject(new TH2F("hEmcJetPtEta", "Emc Jets;#it{p}_{T};#eta", nPtBins, kMinPt, kMaxPt, nEtaBins, kMinEta, kMaxEta));
    hEmcJetPtPhi.setObject(new TH2F("hEmcJetPtPhi", "Emc Jets;#it{p}_{T};#phi", nPtBins, kMinPt, kMaxPt, nPhiBins, kMinPhi, kMaxPhi));
    hSelectedJetPtEta.setObject(new TH2F("hSelectedJetPtEta", "Selected Jets;#it{p}_{T};#eta", nPtBins, kMinPt, kMaxPt, nEtaBins, kMinEta, kMaxEta));
    hSelectedJetPtPhi.setObject(new TH2F("hSelectedJetPtPhi", "Selected Jets;#it{p}_{T};#phi", nPtBins, kMinPt, kMaxPt, nPhiBins, kMinPhi, kMaxPhi));
    hSelectedJetPtEtaLow.setObject(new TH2F("hSelectedJetPtEtaLow", "Selected Jets (low threshold);#it{p}_{T};#eta", nPtBins, kMinPt, kMaxPt, nEtaBins, kMinEta, kMaxEta));
    hSelectedJetPtPhiLow.setObject(new TH2F("hSelectedJetPtPhiLow", "Selected Jets (low threshold);#it{p}_{T};#phi", nPtBins, kMinPt, kMaxPt, nPhiBins, kMinPhi, kMaxPhi));
    hSelectedGammaEMCALPtEtaVeryHigh.setObject(new TH2F("hSelectedGammaEMCALEtaVeryHigh", Form("Selected Gammas EMCAL (very high threshold);%s;#eta", (f_ObservalbeGammaTrigger == 0) ? "E (GeV)" : "#it{p}_{T}"), nPtBins, kMinPt, kMaxPt / 2, nEtaBins, kMinEta, kMaxEta));
    hSelectedGammaEMCALPtPhiVeryHigh.setObject(new TH2F("hSelectedGammaEMCALPhiVeryHigh", Form("Selected Gammas EMCAL (very high threshold);%s;#phi", (f_ObservalbeGammaTrigger == 0) ? "E (GeV)" : "#it{p}_{T}"), nPtBins, kMinPt, kMaxPt / 2, nPhiBins, kMinPhi, kMaxPhi));
    hSelectedGammaDCALPtEtaVeryHigh.setObject(new TH2F("hSelectedGammaDCALEtaVeryHigh", Form("Selected Gammas DCAL (very high threshold);%s;#eta", (f_ObservalbeGammaTrigger == 0) ? "E (GeV)" : "#it{p}_{T}"), nPtBins, kMinPt, kMaxPt / 2, nEtaBins, kMinEta, kMaxEta));
    hSelectedGammaDCALPtPhiVeryHigh.setObject(new TH2F("hSelectedGammaDCALPhiVeryHigh", Form("Selected Gammas DCAL (very high threshold);%s;#phi", (f_ObservalbeGammaTrigger == 0) ? "E (GeV)" : "#it{p}_{T}"), nPtBins, kMinPt, kMaxPt / 2, nPhiBins, kMinPhi, kMaxPhi));
    hSelectedGammaEMCALPtEta.setObject(new TH2F("hSelectedGammaEMCALPtEta", Form("Selected Gammas EMCAL;%s;#eta", (f_ObservalbeGammaTrigger == 0) ? "E (GeV)" : "#it{p}_{T}"), nPtBins, kMinPt, kMaxPt / 2, nEtaBins, kMinEta, kMaxEta));
    hSelectedGammaEMCALPtPhi.setObject(new TH2F("hSelectedGammaEMCALPtPhi", Form("Selected Gammas EMCAL;%s;#phi", (f_ObservalbeGammaTrigger == 0) ? "E (GeV)" : "#it{p}_{T}"), nPtBins, kMinPt, kMaxPt / 2, nPhiBins, kMinPhi, kMaxPhi));
    hSelectedGammaDCALPtEta.setObject(new TH2F("hSelectedGammaDCALPtEta", Form("Selected Gammas DCAL;%s;#eta", (f_ObservalbeGammaTrigger == 0) ? "E (GeV)" : "#it{p}_{T}"), nPtBins, kMinPt, kMaxPt / 2, nEtaBins, kMinEta, kMaxEta));
    hSelectedGammaDCALPtPhi.setObject(new TH2F("hSelectedGammaDCALPtPhi", Form("Selected Gammas DCAL;%s;#phi", (f_ObservalbeGammaTrigger == 0) ? "E (GeV)" : "#it{p}_{T}"), nPtBins, kMinPt, kMaxPt / 2, nPhiBins, kMinPhi, kMaxPhi));
    hSelectedGammaEMCALPtEtaLow.setObject(new TH2F("hSelectedGammaEMCALPtEtaLow", Form("Selected Gammas EMCAL (low threshold);%s;#eta", (f_ObservalbeGammaTrigger == 0) ? "E (GeV)" : "#it{p}_{T}"), nPtBins, kMinPt, kMaxPt / 2, nEtaBins, kMinEta, kMaxEta));
    hSelectedGammaEMCALPtPhiLow.setObject(new TH2F("hSelectedGammaEMCALPtPhiLow", Form("Selected Gammas EMCAL (low threshold);%s;#phi", (f_ObservalbeGammaTrigger == 0) ? "E (GeV)" : "#it{p}_{T}"), nPtBins, kMinPt, kMaxPt / 2, nPhiBins, kMinPhi, kMaxPhi));
    hSelectedGammaDCALPtEtaLow.setObject(new TH2F("hSelectedGammaDCALPtEtaLow", Form("Selected Gammas DCAL (low threshold);%s;#eta", (f_ObservalbeGammaTrigger == 0) ? "E (GeV)" : "#it{p}_{T}"), nPtBins, kMinPt, kMaxPt / 2, nEtaBins, kMinEta, kMaxEta));
    hSelectedGammaDCALPtPhiLow.setObject(new TH2F("hSelectedGammaDCALPtPhiLow", Form("Selected Gammas DCAL (low threshold);%s;#phi", (f_ObservalbeGammaTrigger == 0) ? "E (GeV)" : "#it{p}_{T}"), nPtBins, kMinPt, kMaxPt / 2, nPhiBins, kMinPhi, kMaxPhi));
    hSelectedGammaEMCALPtEtaVeryLow.setObject(new TH2F("hSelectedGammaEMCALEtaVeryLow", Form("Selected Gammas EMCAL (very low threshold);%s;#eta", (f_ObservalbeGammaTrigger == 0) ? "E (GeV)" : "#it{p}_{T}"), nPtBins, kMinPt, kMaxPt / 2, nEtaBins, kMinEta, kMaxEta));
    hSelectedGammaEMCALPtPhiVeryLow.setObject(new TH2F("hSelectedGammaEMCALPhiVeryLow", Form("Selected Gammas EMCAL (very low threshold);%s;#phi", (f_ObservalbeGammaTrigger == 0) ? "E (GeV)" : "#it{p}_{T}"), nPtBins, kMinPt, kMaxPt / 2, nPhiBins, kMinPhi, kMaxPhi));
    hSelectedGammaDCALPtEtaVeryLow.setObject(new TH2F("hSelectedGammaDCALEtaVeryLow", Form("Selected Gammas DCAL (very low threshold);%s;#eta", (f_ObservalbeGammaTrigger == 0) ? "E (GeV)" : "#it{p}_{T}"), nPtBins, kMinPt, kMaxPt / 2, nEtaBins, kMinEta, kMaxEta));
    hSelectedGammaDCALPtPhiVeryLow.setObject(new TH2F("hSelectedGammaDCALPhiVeryLow", Form("Selected Gammas DCAL (very low threshold);%s;#phi", (f_ObservalbeGammaTrigger == 0) ? "E (GeV)" : "#it{p}_{T}"), nPtBins, kMinPt, kMaxPt / 2, nPhiBins, kMinPhi, kMaxPhi));

    hMaxJetPt.setObject(new TH1D("hMaxJetPt", "Max. jet pt", nPtBins, kMinPt, kMaxPt));
    hSelectMaxJetPt.setObject(new TH1D("hSelectMaxJetPt", "Max. jet pt selected collisions", nPtBins, kMinPt, kMaxPt));
    hSelectMaxJetPtLow.setObject(new TH1D("hSelectMaxJetPtLow", "Max. jet pt selected collisions (low threshold)", nPtBins, kMinPt, kMaxPt));
    hMaxClusterEMCAL.setObject(new TH1D("hMaxClusterEMCAL", "Max.cluster pt EMCAL", nPtBins, kMinPt, kMaxPt));
    hMaxClusterDCAL.setObject(new TH1D("hMaxClusterDCAL", "Max. cluster pt DCAL", nPtBins, kMinPt, kMaxPt));
    hMaxClusterEMCALMinBias.setObject(new TH1D("hMaxClusterEMCALMinBias", "Max.cluster pt EMCAL (min. bias)", nPtBins, kMinPt, kMaxPt));
    hMaxClusterDCALMinBias.setObject(new TH1D("hMaxClusterDCALMinBias", "Max. cluster pt DCAL (mib. bias)", nPtBins, kMinPt, kMaxPt));
    hMaxClusterEMCALLevel0.setObject(new TH1D("hMaxClusterEMCALLevel0", "Max.cluster pt EMCAL (Level-0)", nPtBins, kMinPt, kMaxPt));
    hMaxClusterDCALLevel0.setObject(new TH1D("hMaxClusterDCALLevel0", "Max. cluster pt DCAL (Level-0)", nPtBins, kMinPt, kMaxPt));
    hSelectGammaVeryHighMaxClusterEMCAL.setObject(new TH1D("hSelectGammaVeryHighMaxClusterEMCAL", "Max. cluster pt selected Gamms EMCAL (very high threshold)", nPtBins, kMinPt, kMaxPt));
    hSelectGammaMaxClusterEMCAL.setObject(new TH1D("hSelectGammaMaxClusterEMCAL", "Max. cluster pt selected Gamms EMCAL", nPtBins, kMinPt, kMaxPt));
    hSelectGammaLowMaxClusterEMCAL.setObject(new TH1D("hSelectGammaLowMaxClusterEMCAL", "Max. cluster pt selected Gamms EMCAL (low threshold)", nPtBins, kMinPt, kMaxPt));
    hSelectGammaVeryLowMaxClusterEMCAL.setObject(new TH1D("hSelectGammaVeryLowMaxClusterEMCAL", "Max. cluster pt selected Gamms EMCAL (very low threshold)", nPtBins, kMinPt, kMaxPt));
    hSelectGammaVeryHighMaxClusterDCAL.setObject(new TH1D("hSelectGammaVeryHighMaxClusterDCAL", "Max. cluster pt selected Gamms DCAL (very high threshold)", nPtBins, kMinPt, kMaxPt));
    hSelectGammaMaxClusterDCAL.setObject(new TH1D("hSelectGammaMaxClusterDCAL", "Max. cluster pt selected Gamms DCAL", nPtBins, kMinPt, kMaxPt));
    hSelectGammaLowMaxClusterDCAL.setObject(new TH1D("hSelectGammaLowMaxClusterDCAL", "Max. cluster pt selected Gamms DCAL (low threshold)", nPtBins, kMinPt, kMaxPt));
    hSelectGammaVeryLowMaxClusterDCAL.setObject(new TH1D("hSelectGammaVertLowMaxClusterDCAL", "Max. cluster pt selected Gamms DCAL (very low threshold)", nPtBins, kMinPt, kMaxPt));

    LOG(info) << "Jet trigger: " << (b_doJetTrigger ? "on" : "off");
    LOG(info) << "Gamma trigger: " << (b_doJetTrigger ? "on" : "off");
    LOG(info) << "Thresholds gamma trigger (L0-triggered runs): EG1 " << f_gammaPtMinEMCALHigh << " GeV, DG1 " << f_gammaPtMinDCALHigh << " GeV, EG2 " << f_gammaPtMinEMCALLow << " GeV, DG2 " << f_gammaPtMinDCALLow << " GeV";
    LOG(info) << "Thresholds gamma trigger (MB-triggered runs): EG1 " << f_gammaPtMinEMCALHighMB << " GeV, DG1 " << f_gammaPtMinDCALHighMB << " GeV, EG2 " << f_gammaPtMinEMCALLowMB << " GeV, DG2 " << f_gammaPtMinDCALLowMB << " GeV";
    LOG(info) << "Gamma trigger observable: " << (f_ObservalbeGammaTrigger == 0 ? "Energy" : "pt");
    LOG(info) << "Jet trigger: R: " << (static_cast<double>(f_jetR) / 100.) << ", pt (high threshold) > " << f_jetPtMin << " GeV/c, (low threshold) > " << f_jetPtMinLow << " GeV/c, cluster(" << f_clusterPtMin << ")";
    LOG(info) << "Jet thresholds min-bias triggered runs: High " << f_jetPtMinMB << " GeV/c, low " << f_jetPtMinLowMB << " GeV/c";
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

  template <typename JetIterator>
  Bool_t isJetInEmcal(JetIterator const& jet)
  {
    double emcalEtaMin = -0.7, emcalEtaMax = 0.7, emcalPhiMin = 1.40, emcalPhiMax = 3.26; // Phi: 80 - 187 deg
    double R = jet.r() * 1e-2;                                                            // Jet R is saved as round(100*R)
    if ((jet.eta() >= emcalEtaMin + R) && (jet.eta() <= emcalEtaMax - R) && (jet.phi() >= emcalPhiMin + R) && (jet.phi() <= emcalPhiMax - R)) {
      return true;
    }
    return false;
  }

  EMCALHWTriggerConfiguration getHWTriggerConfiguration(o2::ctp::CTPConfiguration& ctpconfig)
  {
    EMCALHWTriggerConfiguration result = EMCALHWTriggerConfiguration::UNKNOWN;
    bool hasMinBias = false, hasL0 = false;
    for (auto& cls : ctpconfig.getCTPClasses()) {
      auto trgclsname = boost::algorithm::to_upper_copy(cls.name);
      if (trgclsname.find("-EMC") == std::string::npos) {
        // Not an EMCAL trigger class
        continue;
      }
      std::vector<std::string> tokens;
      std::stringstream tokenizer(trgclsname.data());
      std::string buf;
      while (std::getline(tokenizer, buf, '-')) {
        tokens.emplace_back(buf);
      }
      if (tokens[1] != "B") {
        // Not bucket in both beams
        continue;
      }
      LOG(info) << "Found trigger class: " << trgclsname;
      if (tokens[0] == "C0TVX") {
        hasMinBias = true;
      }
      if (tokens[0] == "CTVXEMC" || tokens[0] == "CTVXDMC") {
        hasL0 = true;
      }
    }
    if (hasL0) {
      result = EMCALHWTriggerConfiguration::EMC_TRIGGERD;
    } else if (hasMinBias) {
      result = EMCALHWTriggerConfiguration::MB_ONLY;
    }
    return result;
  }

  float getGammaThreshold(o2::emcal::AcceptanceType_t subdet, ThresholdType_t thresholdt, bool mblike)
  {
    float threshold = 0;
    // Very-high and very-low threshold only needed and supported for trigger-like data
    switch (subdet) {
      case o2::emcal::AcceptanceType_t::EMCAL_ACCEPTANCE: {
        switch (thresholdt) {
          case ThresholdType_t::VERY_HIGH_THRESHOLD: {
            threshold = mblike ? FLT_MAX : f_gammaPtMinEMCALVeryHigh;
            break;
          }
          case ThresholdType_t::HIGH_THRESHOLD: {
            threshold = mblike ? f_gammaPtMinEMCALHighMB : f_gammaPtMinEMCALHigh;
            break;
          }
          case ThresholdType_t::LOW_THRESHOLD: {
            threshold = mblike ? f_gammaPtMinEMCALLowMB : f_gammaPtMinEMCALLow;
            break;
          }
          case ThresholdType_t::VERY_LOW_THRESHOLD: {
            threshold = mblike ? FLT_MAX : f_gammaPtMinEMCALVeryLow;
            break;
          }
        }
        break;
      }
      case o2::emcal::AcceptanceType_t::DCAL_ACCEPTANCE: {
        switch (thresholdt) {
          case ThresholdType_t::VERY_HIGH_THRESHOLD: {
            threshold = mblike ? FLT_MAX : f_gammaPtMinDCALVeryHigh;
            break;
          }
          case ThresholdType_t::HIGH_THRESHOLD: {
            threshold = mblike ? f_gammaPtMinDCALHighMB : f_gammaPtMinDCALHigh;
            break;
          }
          case ThresholdType_t::LOW_THRESHOLD: {
            threshold = mblike ? f_gammaPtMinDCALLowMB : f_gammaPtMinDCALLow;
            break;
          }
          case ThresholdType_t::VERY_LOW_THRESHOLD: {
            threshold = mblike ? FLT_MAX : f_gammaPtMinDCALVeryLow;
            break;
          }
        }
        break;
      }
      case o2::emcal::AcceptanceType_t::NON_ACCEPTANCE: {
        threshold = FLT_MAX;
        break;
      }
    }
    return threshold;
  }

  double getJetThreshold(ThresholdType_t thresholdt, bool mblike)
  {
    float threshold = 0.;
    switch (thresholdt) {
      case ThresholdType_t::HIGH_THRESHOLD:
        threshold = mblike ? f_jetPtMinMB : f_jetPtMin;
        break;
      case ThresholdType_t::LOW_THRESHOLD:
        threshold = mblike ? f_jetPtMinLowMB : f_jetPtMinLow;
        break;
      default:
        // Very low and very high threshold not supported for jet trigger
        threshold = FLT_MAX;
    }
    return threshold;
  }

  Bool_t isEvtSelectedJet(double const& jetpt, ThresholdType_t thresholdt)
  {
    float threshold = 0;
    switch (mHardwareTriggerConfig) {
      case EMCALHWTriggerConfiguration::EMC_TRIGGERD: {
        threshold = getJetThreshold(thresholdt, false);
        break;
      }
      case EMCALHWTriggerConfiguration::MB_ONLY: {
        threshold = getJetThreshold(thresholdt, true);
        break;
      }
      case EMCALHWTriggerConfiguration::UNKNOWN: {
        // No EMCAL trigger present in run - cannot select thresholds
        return false;
      }
    }
    if (jetpt > threshold) {
      return true;
    }
    return false;
  }

  Bool_t isEvtSelectedGamma(double const& gammapt, o2::emcal::AcceptanceType_t subdet, ThresholdType_t thresholdt)
  {
    float threshold = 0;
    switch (mHardwareTriggerConfig) {
      case EMCALHWTriggerConfiguration::EMC_TRIGGERD: {
        threshold = getGammaThreshold(subdet, thresholdt, false);
        break;
      }
      case EMCALHWTriggerConfiguration::MB_ONLY: {
        threshold = getGammaThreshold(subdet, thresholdt, true);
        break;
      }
      case EMCALHWTriggerConfiguration::UNKNOWN: {
        // No EMCAL trigger present in run - cannot select thresholds
        return false;
      }
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

  bool isEMCALMinBias(collisionInfo const& collision) const
  {
    return collision.alias_bit(kTVXinEMC);
  }

  bool isEMCALLevel0(collisionInfo const& collision) const
  {
    return collision.alias_bit(kEMC7) || collision.alias_bit(kDMC7);
  }

  bool isEMCALLevel1(collisionInfo const& collision) const
  {
    return collision.alias_bit(kEG1) || collision.alias_bit(kEG2) || collision.alias_bit(kDG1) || collision.alias_bit(kDG2) || collision.alias_bit(kEJ1) || collision.alias_bit(kEJ2) || collision.alias_bit(kDJ1) || collision.alias_bit(kDJ2);
  }

  bool hasEMCALData(collisionInfo const& collision) const
  {
    return isEMCALMinBias(collision) || isEMCALLevel0(collision) || isEMCALLevel1(collision);
  }

  void runGammaTrigger(const selectedClusters& clusters, std::bitset<EMCALHardwareTrigger::TRG_NTriggers> hardwaretriggers, std::array<bool, kCategories>& keepEvent)
  {
    double maxClusterObservableEMCAL = -1., maxClusterObservableDCAL = -1.;
    static constexpr std::array<ThresholdType_t, 4> thresholds = {{ThresholdType_t::VERY_HIGH_THRESHOLD, ThresholdType_t::HIGH_THRESHOLD, ThresholdType_t::LOW_THRESHOLD, ThresholdType_t::VERY_LOW_THRESHOLD}};
    static constexpr std::array<o2::emcal::AcceptanceType_t, 2> subdets = {{o2::emcal::AcceptanceType_t::EMCAL_ACCEPTANCE, o2::emcal::AcceptanceType_t::DCAL_ACCEPTANCE}};
    std::array<TH2*, 8> acceptanceHistsPtEta{{hSelectedGammaEMCALPtEtaVeryHigh.object.get(), hSelectedGammaDCALPtEtaVeryHigh.object.get(), hSelectedGammaEMCALPtEta.object.get(), hSelectedGammaDCALPtEta.object.get(), hSelectedGammaEMCALPtEtaLow.object.get(), hSelectedGammaDCALPtEtaLow.object.get(), hSelectedGammaEMCALPtEtaVeryLow.object.get(), hSelectedGammaDCALPtEtaVeryLow.object.get()}},
      acceptanceHistsPtPhi{{hSelectedGammaEMCALPtPhiVeryHigh.object.get(), hSelectedGammaDCALPtPhiVeryHigh.object.get(), hSelectedGammaEMCALPtPhi.object.get(), hSelectedGammaDCALPtPhi.object.get(), hSelectedGammaEMCALPtPhiLow.object.get(), hSelectedGammaDCALPtPhiLow.object.get(), hSelectedGammaEMCALPtPhiVeryLow.object.get(), hSelectedGammaDCALPtPhiVeryLow.object.get()}};
    std::array<TH1*, 8> maxClusterPtHists = {{hSelectGammaVeryHighMaxClusterEMCAL.object.get(), hSelectGammaVeryHighMaxClusterDCAL.object.get(), hSelectGammaMaxClusterEMCAL.object.get(), hSelectGammaMaxClusterDCAL.object.get(), hSelectGammaLowMaxClusterEMCAL.object.get(), hSelectGammaLowMaxClusterDCAL.object.get(), hSelectGammaVeryLowMaxClusterEMCAL.object.get(), hSelectGammaVeryLowMaxClusterDCAL.object.get()}};
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
      if (hardwaretriggers.test(EMCALHardwareTrigger::TRG_EMC7) || hardwaretriggers.test(EMCALHardwareTrigger::TRG_DMC7)) {
        hEmcClusterPtEtaLevel0->Fill(observableGamma, cluster.eta());
        hEmcClusterPtPhiLevel0->Fill(observableGamma, cluster.phi());
      }
      if (hardwaretriggers.test(EMCALHardwareTrigger::TRG_MB)) {
        hEmcClusterPtEtaMinBias->Fill(observableGamma, cluster.eta());
        hEmcClusterPtPhiMinBias->Fill(observableGamma, cluster.phi());
      }
      analysedClusters.push_back({static_cast<float>(observableGamma), cluster.eta(), cluster.phi()});
    }
    hMaxClusterEMCAL->Fill(maxClusterObservableEMCAL);
    hMaxClusterDCAL->Fill(maxClusterObservableDCAL);
    if (hardwaretriggers.test(EMCALHardwareTrigger::TRG_MB)) {
      hMaxClusterEMCALMinBias->Fill(maxClusterObservableEMCAL);
      hMaxClusterDCALMinBias->Fill(maxClusterObservableDCAL);
    }
    if (hardwaretriggers.test(EMCALHardwareTrigger::TRG_EMC7)) {
      hMaxClusterEMCALLevel0->Fill(maxClusterObservableEMCAL);
    }
    if (hardwaretriggers.test(EMCALHardwareTrigger::TRG_DMC7)) {
      hMaxClusterDCALLevel0->Fill(maxClusterObservableDCAL);
    }
    for (decltype(thresholds.size()) ithreshold = 0; ithreshold < thresholds.size(); ithreshold++) {
      if (thresholds[ithreshold] == ThresholdType_t::VERY_LOW_THRESHOLD) {
        // Accept very-low threshold only in case of min. bias input to EMCAL/DCAL
        if (!hardwaretriggers.test(EMCALHardwareTrigger::TRG_MB)) {
          continue;
        }
      }
      for (decltype(thresholds.size()) isubdet = 0; isubdet < subdets.size(); isubdet++) {
        if (isEvtSelectedGamma(subdets[isubdet] == o2::emcal::AcceptanceType_t::EMCAL_ACCEPTANCE ? maxClusterObservableEMCAL : maxClusterObservableDCAL, subdets[isubdet], thresholds[ithreshold])) {
          maxClusterPtHists[ithreshold * subdets.size() + isubdet]->Fill(subdets[isubdet] == o2::emcal::AcceptanceType_t::EMCAL_ACCEPTANCE ? maxClusterObservableEMCAL : maxClusterObservableDCAL);
          keepEvent[kGammaVeryHighPtEMCAL + ithreshold * subdets.size() + isubdet] = true;
          for (auto& cluster : analysedClusters) {
            acceptanceHistsPtEta[ithreshold * subdets.size() + isubdet]->Fill(cluster.mTriggerObservable, cluster.mEta);
            acceptanceHistsPtPhi[ithreshold * subdets.size() + isubdet]->Fill(cluster.mTriggerObservable, cluster.mPhi);
          }
        }
      }
    }
  }

  template <typename JetCollection>
  void runJetTrigger(JetCollection const& jets, selectedClusters const& clusters, std::array<bool, kCategories>& keepEvent)
  {
    int jettype = 0;
    if constexpr (std::is_same<JetCollection, filteredFullJets>::value) {
      LOG(detail) << "Running full jet filter";
      jettype = 0;
    } else if constexpr (std::is_same<JetCollection, filteredNeutralJets>::value) {
      LOG(detail) << "Running neutral jet filter";
      jettype = 1;
    } else {
      LOG(error) << "Unknown jet type";
      return;
    }
    double maxSelectedJetPt = 0.;
    // Store pt of leading jet inside of the emcal to trigger on
    static constexpr std::array<ThresholdType_t, 2>
      thresholds = {{ThresholdType_t::HIGH_THRESHOLD, ThresholdType_t::LOW_THRESHOLD}};
    std::array<TH2*, 2> acceptanceHistsPtEtaSel = {{hSelectedJetPtEta.object.get(), hSelectedJetPtEtaLow.object.get()}},
                        acceptanceHistsPtPhiSel = {{hSelectedJetPtPhi.object.get(), hSelectedJetPtPhiLow.object.get()}},
                        acceptanceHistsClusterPtEtaSel = {{hSelectedClusterPtEta.object.get(), hSelectedClusterPtEtaLow.object.get()}},
                        acceptanceHistsClusterPtPhiSel = {{hSelectedClusterPtPhi.object.get(), hSelectedClusterPtPhiLow.object.get()}};
    std::array<TH1*, 2> maxJetPtHists = {{hSelectMaxJetPt.object.get(), hSelectMaxJetPtLow.object.get()}};
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
    hMaxJetPt->Fill(maxSelectedJetPt);
    for (decltype(thresholds.size()) ithreshold = 0; ithreshold < thresholds.size(); ithreshold++) {
      if (isEvtSelectedJet(maxSelectedJetPt, thresholds[ithreshold])) {
        maxJetPtHists[ithreshold]->Fill(maxSelectedJetPt);
        std::size_t triggerIndex = kCategories;
        if (jettype == 0) {
          switch (thresholds[ithreshold]) {
            case ThresholdType_t::HIGH_THRESHOLD:
              triggerIndex = kJetFullHighPt;
              break;
            case ThresholdType_t::LOW_THRESHOLD:
              triggerIndex = kJetFullLowPt;
              break;
            default:
              // Very high / very low thresholds not supported for jet trigger
              break;
          }
        } else {
          switch (thresholds[ithreshold]) {
            case ThresholdType_t::HIGH_THRESHOLD:
              triggerIndex = kJetNeutralHighPt;
              break;
            case ThresholdType_t::LOW_THRESHOLD:
              triggerIndex = kJetNeutralLowPt;
              break;
            default:
              // Very high / very low thresholds not supported for jet trigger
              break;
          }
        }
        keepEvent[triggerIndex] = true;
        for (const auto& jet : jets) {
          if (isJetInEmcal(jet)) {
            acceptanceHistsPtEtaSel[ithreshold]->Fill(jet.pt(), jet.eta());
            acceptanceHistsPtPhiSel[ithreshold]->Fill(jet.pt(), jet.phi());
          }
        }
        for (const auto& cluster : clusters) {
          double clusterPt = cluster.energy() / std::cosh(cluster.eta());
          acceptanceHistsClusterPtEtaSel[ithreshold]->Fill(clusterPt, cluster.eta());
          acceptanceHistsClusterPtPhiSel[ithreshold]->Fill(clusterPt, cluster.phi());
        }
      }
    }
  }

  template <typename JetCollection>
  void runTrigger(collisionInfo const& collision, JetCollection const& jets, selectedClusters const& clusters, BCsWithBcSelsRun3 const& bcs)
  {
    int run = bcs.iteratorAt(0).runNumber();
    // extract bc pattern from CCDB for data or anchored MC only
    if (run != lastRun && run >= 500000) {
      lastRun = run;
      constexpr bool isFullJets = std::is_same<JetCollection, filteredFullJets>::value;
      constexpr bool isNeutralJets = std::is_same<JetCollection, filteredNeutralJets>::value;
      int64_t ts = bcs.iteratorAt(0).timestamp();
      // Request run number in addition to timestamp in order to handle concurrent runs
      std::map<std::string, std::string> metadata;
      metadata["runNumber"] = std::to_string(run);
      mHardwareTriggerConfig = getHWTriggerConfiguration(*(ccdb->getSpecific<o2::ctp::CTPConfiguration>("CTP/Config/Config", ts, metadata)));
      switch (mHardwareTriggerConfig) {
        case EMCALHWTriggerConfiguration::MB_ONLY: {
          LOG(info) << "Found hardware trigger configuration Min. bias only";
          LOG(info) << "Gamma thresholds: High " << f_gammaPtMinEMCALHighMB << ", Low " << f_gammaPtMinEMCALLowMB << " (EMCAL), High " << f_gammaPtMinDCALHighMB << ", Low " << f_gammaPtMinDCALLowMB << " (DCAL)";
          LOG(info) << "No extreme gamma thresholds used";
          LOG(info) << "Jet thresholds: High " << f_jetPtMinMB << ", Low " << f_jetPtMinLowMB;
          LOG(info) << "Jet type: full jets " << (isFullJets ? "yes" : "no") << ", neutral jets " << (isNeutralJets ? "yes" : "no");
          break;
        }
        case EMCALHWTriggerConfiguration::EMC_TRIGGERD: {
          LOG(info) << "Found hardware trigger configuration L0-triggered";
          LOG(info) << "Gamma thresholds: High " << f_gammaPtMinEMCALHigh << ", Low " << f_gammaPtMinEMCALLow << " (EMCAL), High " << f_gammaPtMinDCALHigh << ", Low " << f_gammaPtMinDCALLow << " (DCAL)";
          LOG(info) << "Extreme gamma thresholds: Very high " << f_gammaPtMinEMCALVeryHigh << ", very low " << f_gammaPtMinEMCALVeryLow << " (EMCAL), very high " << f_gammaPtMinDCALVeryHigh << ", very low " << f_gammaPtMinDCALVeryLow << " (DCAL)";
          LOG(info) << "Jet thresholds: High " << f_jetPtMin << ", Low " << f_jetPtMinLow;
          LOG(info) << "Jet type: full jets " << (isFullJets ? "yes" : "no") << ", neutral jets " << (isNeutralJets ? "yes" : "no");
          break;
        }
        case EMCALHWTriggerConfiguration::UNKNOWN: {
          std::cout << "No EMCAL trigger class found for run, event selection will not be possible" << std::endl;
          break;
        }
      }
    }

    std::array<bool, kCategories> keepEvent;
    std::fill(keepEvent.begin(), keepEvent.end(), false);

    if (!b_IgnoreEmcalFlag && !hasEMCALData(collision)) {
      tags(keepEvent[0], keepEvent[1], keepEvent[2], keepEvent[3], keepEvent[4], keepEvent[5], keepEvent[6], keepEvent[7], keepEvent[8], keepEvent[9], keepEvent[10], keepEvent[11], keepEvent[12]);
      return; // Skip events where EMCAL is not live
    }

    if (b_PublishReadoutTrigger && isEMCALMinBias(collision)) {
      keepEvent[kEMCALReadout] = true;
    }

    std::bitset<EMCALHardwareTrigger::TRG_NTriggers> hardwaretriggers;
    if (collision.alias_bit(triggerAliases::kTVXinEMC)) {
      hardwaretriggers.set(EMCALHardwareTrigger::TRG_MB, true);
    }
    if (collision.alias_bit(triggerAliases::kEMC7)) {
      hardwaretriggers.set(EMCALHardwareTrigger::TRG_EMC7, true);
    }
    if (collision.alias_bit(triggerAliases::kDMC7)) {
      hardwaretriggers.set(EMCALHardwareTrigger::TRG_DMC7, true);
    }

    if (b_doJetTrigger) {
      runJetTrigger(jets, clusters, keepEvent);
    }

    if (b_doGammaTrigger) {
      runGammaTrigger(clusters, hardwaretriggers, keepEvent);
    }

    for (int iDecision{0}; iDecision < kCategories; iDecision++) {
      if (keepEvent[iDecision]) {
        hProcessedEvents->Fill(iDecision);
      }
    }
    tags(keepEvent[0], keepEvent[1], keepEvent[2], keepEvent[3], keepEvent[4], keepEvent[5], keepEvent[6], keepEvent[7], keepEvent[8], keepEvent[9], keepEvent[10], keepEvent[11], keepEvent[12]);
  }

  void processFullJetTrigger(collisionInfo const& collision, filteredFullJets const& jets, selectedClusters const& clusters, BCsWithBcSelsRun3 const& bcs)
  {
    // Trigger selection (full jet case)
    runTrigger(collision, jets, clusters, bcs);
  }
  PROCESS_SWITCH(fullJetFilter, processFullJetTrigger, "run full jet triggere code", true);

  void processNeutralJetTrigger(collisionInfo const& collision, filteredNeutralJets const& jets, selectedClusters const& clusters, BCsWithBcSelsRun3 const& bcs)
  {
    // Trigger selection (neutral jet case)
    runTrigger(collision, jets, clusters, bcs);
  }
  PROCESS_SWITCH(fullJetFilter, processNeutralJetTrigger, "run neutral jet triggere code", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{
  return WorkflowSpec{
    adaptAnalysisTask<fullJetFilter>(cfg)};
}
