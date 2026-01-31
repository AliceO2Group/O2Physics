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

/// \file statPromptPhoton.cxx
/// \brief Reconstruction of Phi yield through track-track Minv correlations for resonance hadrochemistry analysis.
///
///
///  \author Adrian Fereydon Nassirpour <adrian.fereydon.nassirpour@cern.ch>

#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/DataModel/EMCALClusters.h"
#include "PWGJE/DataModel/Jet.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CommonConstants/PhysicsConstants.h"
#include "CommonDataFormat/InteractionRecord.h"
#include "DataFormatsEMCAL/AnalysisCluster.h"
#include "DataFormatsEMCAL/Cell.h"
#include "DataFormatsEMCAL/Constants.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/Propagator.h"
#include "EMCALBase/Geometry.h"
#include "EMCALCalib/BadChannelMap.h"
#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"
#include <CCDB/BasicCCDBManager.h>

#include <TLorentzVector.h>
#include <TVector2.h>

#include <iostream>
#include <optional>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct statPromptPhoton {

  SliceCache cache;
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Configurable<double> cfgMaxDCArToPVcut{"cfgMaxDCArToPVcut", 0.5, "Track DCAr cut to PV Maximum"};
  Configurable<double> cfgMaxDCAzToPVcut{"cfgMaxDCAzToPVcut", 2.0, "Track DCAz cut to PV Maximum"};
  Configurable<bool> cfgPrimaryTrack{"cfgPrimaryTrack", true, "Primary track selection"};
  Configurable<bool> cfgConnectedToPV{"cfgConnectedToPV", true, "PV contributor track selection"};
  Configurable<bool> cfgGlobalWoDCATrack{"cfgGlobalWoDCATrack", true, "Global track selection without DCA"};
  Configurable<double> cfgnFindableTPCClusters{"cfgnFindableTPCClusters", 50, "nFindable TPC Clusters"};
  Configurable<double> cfgnTPCCrossedRows{"cfgnTPCCrossedRows", 70, "nCrossed TPC Rows"};
  Configurable<double> cfgnRowsOverFindable{"cfgnRowsOverFindable", 1.2, "nRowsOverFindable TPC CLusters"};
  Configurable<double> cfgnTPCChi2{"cfgnTPChi2", 4.0, "nTPC Chi2 per Cluster"};
  Configurable<double> cfgnITSChi2{"cfgnITShi2", 36.0, "nITS Chi2 per Cluster"};
  Configurable<int> cfgClusterDefinition{"clusterDefinition", 10, "cluster definition to be selected, e.g. 10=kV3Default, 0 = kV1Default"};
  Configurable<float> cfgMinTime{"MinTime", -30., "Minimum cluster time for time cut"};
  Configurable<float> cfgMaxTime{"MaxTime", +35., "Maximum cluster time for time cut"};
  Configurable<float> cfgMinClusterEnergy{"MinClusterEnergy", 0.7f, "Minimal cluster energy"};
  Configurable<int> cfgMinNCells{"MinNCelss", 2, "Minimal amount of cells per cluster"};
  Configurable<int> cfgMaxNLM{"MaxNLM", 2, "Maximal amount of local Maxima per cluster"};
  Configurable<bool> cfgExoticContribution{"ExoticContribution", false, "Exotic cluster in the data"};
  Configurable<double> cfgtrkMinPt{"cfgtrkMinPt", 0.15, "set track min pT"};
  Configurable<double> cfgtrkMaxEta{"cfgtrkMaxEta", 0.6, "set track max Eta"};
  Configurable<float> cfgMinR{"MinR", 0.1, "Min. Radii of Delta R cone around photon trigger"};
  Configurable<float> cfgMaxR{"MaxR", 0.4, "Max. Radii of Delta R cone around photon trigger"};
  Configurable<float> cfgMinTrig{"MinTrig", 1, "Min. Trigger energy/momentum"};
  Configurable<float> cfgMaxTrig{"MaxTrig", 5, "Max. Trigger energy/momentum"};
  Configurable<float> cfgVtxCut{"cfgVtxCut", 10.0, "V_z cut selection"};
  Configurable<float> cfgLowM02{"cfgLowM02", 0.1, "Lower-bound M02 cut"};
  Configurable<float> cfgHighM02{"cfgHighM02", 0.3, "Higher-bound M02 cut"};
  Configurable<float> cfgLowClusterE{"cfgLowClusterE", 0.5, "Higher-bound Cluster E cut"};
  Configurable<float> cfgHighClusterE{"cfgHighClusterE", 500, "Lower-bound Cluster E cut"};
  Configurable<bool> cfgEmcTrigger{"cfgEmcTrigger", true, "Require EMC readout for event"};
  Configurable<bool> cfgGeoCut{"cfgGeoCut", true, "Performs Geometric TPC cut"};
  Configurable<bool> cfgPtClusterCut{"cfgPtClusterCut", true, "Performs Pt-dependent cluster-track matching"};
  Configurable<std::string> cfgTrackFilter{"cfgTrackFilter", "globalTracks", "set track selections"};
  Configurable<bool> cfgJETracks{"cfgJETracks", false, "Enables running on derived JE data"};
  Configurable<bool> cfgGenHistograms{"cfgGenHistograms", false, "Enables Generated histograms"};
  Configurable<bool> cfgRecHistograms{"cfgRecHistograms", false, "Enables Reconstructed histograms"};
  Configurable<bool> cfgDataHistograms{"cfgDataHistograms", false, "Enables Data histograms"};
  Configurable<std::string> cfgTriggerMasks{"cfgTriggerMasks", "", "possible JE Trigger masks: fJetChLowPt,fJetChHighPt,fTrackLowPt,fTrackHighPt,fJetD0ChLowPt,fJetD0ChHighPt,fJetLcChLowPt,fJetLcChHighPt,fEMCALReadout,fJetFullHighPt,fJetFullLowPt,fJetNeutralHighPt,fJetNeutralLowPt,fGammaVeryHighPtEMCAL,fGammaVeryHighPtDCAL,fGammaHighPtEMCAL,fGammaHighPtDCAL,fGammaLowPtEMCAL,fGammaLowPtDCAL,fGammaVeryLowPtEMCAL,fGammaVeryLowPtDCAL"};
  Configurable<bool> cfgDebug{"cfgDebug", false, "Enables debug information for local running"};

  int trackFilter = -1;
  std::vector<int> triggerMaskBits;

  // INIT
  void init(InitContext const&)
  {
    std::vector<double> ptBinning = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 12.0, 16.0, 20.0, 25.0, 30.0, 40.0, 50.0, 75.0, 100.0, 150.0, 200.0, 300.0, 500.0};
    AxisSpec pthadAxis = {ptBinning, "#it{p}_{T}^{had sum} [GeV/c]"};

    triggerMaskBits = jetderiveddatautilities::initialiseTriggerMaskBits(cfgTriggerMasks);
    if (cfgJETracks) {
      trackFilter = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(cfgTrackFilter));
    }

    if (cfgRecHistograms) {
      histos.add("REC_nEvents", "REC_nEvents", kTH1F, {{4, 0.0, 4.0}});
      histos.add("REC_Cluster_QA", "REC_Cluster_QA", kTH1F, {{10, -0.5, 9.5}});
      histos.add("REC_PtHadSum_Photon", "REC_PtHadSum_Photon", kTH1F, {pthadAxis});
      histos.add("REC_TrackPhi_photontrigger", "REC_TrackPhi_photontrigger", kTH1F, {{64, 0, 2 * TMath::Pi()}});
      histos.add("REC_TrackEta_photontrigger", "REC_TrackEta_photontrigger", kTH1F, {{100, -1, 1}});
      histos.add("REC_ClusterPhi", "REC_ClusterPhi", kTH1F, {{640 * 2, 0, 2 * TMath::Pi()}});
      histos.add("REC_ClusterEta", "REC_ClusterEta", kTH1F, {{100, -1, 1}});
      histos.add("REC_Track_Pt", "REC_Track_Pt", kTH1F, {{82, -1.0, 40.0}});
      histos.add("REC_Track_Phi", "REC_Track_Phi", kTH1F, {{640 * 2, 0, 2 * TMath::Pi()}});
      histos.add("REC_Track_PhiPrime_Pt", "REC_Track_PhiPrime_Pt", kTH2F, {{640, 0, 2 * TMath::Pi()}, {82, -1.0, 40.0}});
      histos.add("REC_Cluster_PhiPrime_Pt", "REC_Cluster_PhiPrime_Pt", kTH2F, {{640, 0, 2 * TMath::Pi()}, {82, -1.0, 40.0}});
      histos.add("REC_Cluster_PhiPrime_Pt_AC", "REC_Cluster_PhiPrime_Pt_AC", kTH2F, {{640, 0, 2 * TMath::Pi()}, {82, -1.0, 40.0}});
      histos.add("REC_Cluster_PhiPrime_Pt_C", "REC_Cluster_PhiPrime_Pt_C", kTH2F, {{640, 0, 2 * TMath::Pi()}, {82, -1.0, 40.0}});
      histos.add("REC_Cluster_Particle_Pt", "REC_Cluster_Particle_Pt", kTH1F, {{82, -1.0, 40.0}});
      histos.add("REC_Cluster_ParticleWITHtrack_Pt", "REC_Cluster_ParticleWITHtrack_Pt", kTH1F, {{82, -1.0, 40.0}});
      histos.add("REC_Cluster_ParticleWITHtrack_Phi", "REC_Cluster_ParticleWITHtrack_Phi", kTH1F, {{628, -2 * TMath::Pi(), 2 * TMath::Pi()}});
      histos.add("REC_Cluster_ParticleWITHtrack_Eta", "REC_Cluster_ParticleWITHtrack_Eta", kTH1F, {{628, -2 * TMath::Pi(), 2 * TMath::Pi()}});
      histos.add("REC_Cluster_ParticleWITHtrack_TrackPt", "REC_Cluster_ParticleWITHtrack_TrackPt", kTH1F, {{82, -1.0, 40.0}});
      histos.add("REC_Cluster_ParticleWITHtrack_Pt_Phi", "REC_Cluster_ParticleWITHtrack_Pt_Phi", kTH2F, {{82, -1.0, 40.0}, {628, -2 * TMath::Pi(), 2 * TMath::Pi()}});
      histos.add("REC_Cluster_ParticleWITHtrack_Pt_PhiPrime", "REC_Cluster_ParticleWITHtrack_Pt_PhiPrime", kTH2F, {{82, -1.0, 40.0}, {628, -2 * TMath::Pi(), 2 * TMath::Pi()}});
      histos.add("REC_Impurity_ParticleWITHtrack_Pt_PhiPrime", "REC_Impurity_ParticleWITHtrack_Pt_PhiPrime", kTH2F, {{82, -1.0, 40.0}, {628, -2 * TMath::Pi(), 2 * TMath::Pi()}});
      histos.add("REC_Cluster_ParticleWITHtrack_Pt_Eta", "REC_Cluster_ParticleWITHtrack_Pt_Eta", kTH2F, {{82, -1.0, 40.0}, {628, -2 * TMath::Pi(), 2 * TMath::Pi()}});
      histos.add("REC_Cluster_ParticleWITHOUTtrack_Pt", "REC_Cluster_ParticleWITHOUTtrack_Pt", kTH1F, {{82, -1.0, 40.0}});
      histos.add("REC_Cluster_ParticleWITHOUTtrack_Phi", "REC_Cluster_ParticleWITHOUTtrack_Phi", kTH1F, {{628, -2 * TMath::Pi(), 2 * TMath::Pi()}});
      histos.add("REC_Cluster_ParticleWITHOUTtrack_Eta", "REC_Cluster_ParticleWITHOUTtrack_Eta", kTH1F, {{628, -2 * TMath::Pi(), 2 * TMath::Pi()}});
      histos.add("REC_Cluster_ParticleWITHOUTtrack_Pt_Phi", "REC_Cluster_ParticleWITHOUTtrack_Pt_Phi", kTH2F, {{82, -1.0, 40.0}, {628, -2 * TMath::Pi(), 2 * TMath::Pi()}});
      histos.add("REC_Cluster_ParticleWITHOUTtrack_Pt_PhiPrime", "REC_Cluster_ParticleWITHOUTtrack_Pt_PhiPrime", kTH2F, {{82, -1.0, 40.0}, {628, -2 * TMath::Pi(), 2 * TMath::Pi()}});
      histos.add("REC_Cluster_ParticleWITHOUTtrack_Pt_Eta", "REC_Cluster_ParticleWITHOUTtrack_Pt_Eta", kTH2F, {{82, -1.0, 40.0}, {628, -2 * TMath::Pi(), 2 * TMath::Pi()}});
      histos.add("REC_Impurity_ParticleWITHOUTtrack_Pt_PhiPrime", "REC_Impurity_ParticleWITHOUTtrack_Pt_PhiPrime", kTH2F, {{82, -1.0, 40.0}, {628, -2 * TMath::Pi(), 2 * TMath::Pi()}});
      histos.add("REC_TrackPt_ClusterE", "REC_TrackPt_ClusterE", kTH2F, {{82, -1.0, 40.0}, {82, -1.0, 40.0}});
      histos.add("REC_ParticlePt_ClusterE", "REC_ParticlePt_ClusterE", kTH2F, {{82, -1.0, 40.0}, {82, -1.0, 40.0}});
      histos.add("REC_TrackPt_Phi_Eta", "REC_TrackPt_Phi_Eta", kTH2F, {{628, -2 * TMath::Pi(), 2 * TMath::Pi()}, {628, -2 * TMath::Pi(), 2 * TMath::Pi()}});
      histos.add("REC_ParticlePt_Phi_Eta", "REC_ParticlePt_Phi_Eta", kTH2F, {{628, -2 * TMath::Pi(), 2 * TMath::Pi()}, {628, -2 * TMath::Pi(), 2 * TMath::Pi()}});
      histos.add("REC_True_v_Cluster_Phi", "REC_True_v_Cluster_Phi", kTH1F, {{628, -2 * TMath::Pi(), 2 * TMath::Pi()}});
      histos.add("REC_True_v_Cluster_Eta", "REC_True_v_Cluster_Eta", kTH1F, {{628, -2 * TMath::Pi(), 2 * TMath::Pi()}});
      histos.add("REC_TrueImpurity_v_Cluster_Phi", "REC_TrueImpurity_v_Cluster_Phi", kTH1F, {{628, -2 * TMath::Pi(), 2 * TMath::Pi()}});
      histos.add("REC_TrueImpurity_v_Cluster_Eta", "REC_TrueImpurity_v_Cluster_Eta", kTH1F, {{628, -2 * TMath::Pi(), 2 * TMath::Pi()}});
      histos.add("REC_TrueImpurity_v_Cluster_PhiAbs", "REC_TrueImpurity_v_Cluster_PhiAbs", kTH1F, {{628, -2 * TMath::Pi(), 2 * TMath::Pi()}});
      histos.add("REC_TrueImpurity_v_Cluster_EtaAbs", "REC_TrueImpurity_v_Cluster_EtaAbs", kTH1F, {{628, -2 * TMath::Pi(), 2 * TMath::Pi()}});
      histos.add("REC_TrueImpurity_v_Cluster_Phi_Eta", "REC_TrueImpurity_v_Cluster_Phi_Eta", kTH2F, {{628, -2 * TMath::Pi(), 2 * TMath::Pi()}, {628, -2 * TMath::Pi(), 2 * TMath::Pi()}});
      histos.add("REC_Track_v_Cluster_Phi", "REC_Track_v_Cluster_Phi", kTH1F, {{628, -2 * TMath::Pi(), 2 * TMath::Pi()}});
      histos.add("REC_Track_v_Cluster_Eta", "REC_Track_v_Cluster_Eta", kTH1F, {{628, -2 * TMath::Pi(), 2 * TMath::Pi()}});
      histos.add("REC_Track_v_Cluster_Phi_Pt", "REC_Track_v_Cluster_Phi_Pt", kTH2F, {{628, -2 * TMath::Pi(), 2 * TMath::Pi()}, {82, -1.0, 40.0}});
      histos.add("REC_Track_v_Cluster_Eta_Pt", "REC_Track_v_Cluster_Eta_Pt", kTH2F, {{628, -2 * TMath::Pi(), 2 * TMath::Pi()}, {82, -1.0, 40.0}});
      histos.add("REC_Track_v_Cluster_Phi_Eta", "REC_Track_v_Cluster_Phi_Eta", kTH2F, {{628, -2 * TMath::Pi(), 2 * TMath::Pi()}, {628, -2 * TMath::Pi(), 2 * TMath::Pi()}});
      histos.add("REC_Track_v_Cluster_Phi_AC", "REC_Track_v_Cluster_Phi_AC", kTH1F, {{628, -2 * TMath::Pi(), 2 * TMath::Pi()}});
      histos.add("REC_Track_v_Cluster_Eta_AC", "REC_Track_v_Cluster_Eta_AC", kTH1F, {{628, -2 * TMath::Pi(), 2 * TMath::Pi()}});
      histos.add("REC_Track_v_Cluster_Phi_Eta_AC", "REC_Track_v_Cluster_Phi_Eta_AC", kTH2F, {{628, -2 * TMath::Pi(), 2 * TMath::Pi()}, {628, -2 * TMath::Pi(), 2 * TMath::Pi()}});
      histos.add("REC_Track_v_Cluster_Phi_C", "REC_Track_v_Cluster_Phi_C", kTH1F, {{628, -2 * TMath::Pi(), 2 * TMath::Pi()}});
      histos.add("REC_Track_v_Cluster_Eta_C", "REC_Track_v_Cluster_Eta_C", kTH1F, {{628, -2 * TMath::Pi(), 2 * TMath::Pi()}});
      histos.add("REC_Track_v_Cluster_Phi_Eta_C", "REC_Track_v_Cluster_Phi_Eta_C", kTH2F, {{628, -2 * TMath::Pi(), 2 * TMath::Pi()}, {628, -2 * TMath::Pi(), 2 * TMath::Pi()}});
      histos.add("REC_SumPt_BC", "REC_SumPt_BC", kTH1F, {{628, -2 * TMath::Pi(), 2 * TMath::Pi()}});
      histos.add("REC_SumPt_AC", "REC_SumPt_AC", kTH1F, {{628, -2 * TMath::Pi(), 2 * TMath::Pi()}});
      histos.add("REC_M02_BC", "REC_M02_BC", kTH1F, {{628, -2 * TMath::Pi(), 2 * TMath::Pi()}});
      histos.add("REC_M02_AC", "REC_M02_AC", kTH1F, {{628, -2 * TMath::Pi(), 2 * TMath::Pi()}});
      histos.add("REC_Trigger_Purity", "REC_Trigger_Purity", kTH1F, {{4, 0.0, 4.0}});
      histos.add("REC_Trigger_Energy", "REC_Trigger_Energy", kTH1F, {{82, -1.0, 40.0}});
      histos.add("REC_Trigger_Purity_v_Energy", "REC_Trigger_Purity_v_Energy", kTH2F, {{4, 0.0, 4.0}, {82, -1.0, 40.0}});
      histos.add("REC_Trigger_Energy_GOOD", "REC_Trigger_Energy_GOOD", kTH1F, {{82, -1.0, 40.0}});
      histos.add("REC_Trigger_Energy_MISS", "REC_Trigger_Energy_MISS", kTH1F, {{82, -1.0, 40.0}});
      histos.add("REC_Trigger_Energy_FAKE", "REC_Trigger_Energy_FAKE", kTH1F, {{82, -1.0, 40.0}});
      histos.add("REC_All_Energy", "REC_All_Energy", kTH1F, {{82, -1.0, 40.0}});
      histos.add("REC_Impurity_Energy", "REC_Impurity_Energy", kTH1F, {{82, -1.0, 40.0}});
      histos.add("REC_Impurity_Energy_v_Cluster_Phi", "REC_Impurity_Energy_v_Cluster_Phi", kTH2F, {{82, -1.0, 40.0}, {628, -2 * TMath::Pi(), 2 * TMath::Pi()}});
      histos.add("REC_Impurity_Energy_v_ClusterE_Phi", "REC_Impurity_Energy_v_ClusterE_Phi", kTH2F, {{82, -1.0, 40.0}, {628, -2 * TMath::Pi(), 2 * TMath::Pi()}});
      histos.add("REC_Impurity_Energy_v_ClusterEoP_Phi", "REC_Impurity_Energy_v_ClusterEoP_Phi", kTH2F, {{82, -1.0, 40.0}, {628, -2 * TMath::Pi(), 2 * TMath::Pi()}});
      histos.add("REC_Impurity_Energy_v_Cluster_Energy", "REC_Impurity_Energy_v_Cluster_Energy", kTH2F, {{82, -1.0, 40.0}, {82, -1.0, 40.0}});
      histos.add("REC_True_Trigger_Energy", "REC_True_Trigger_Energy", kTH1F, {{82, -1.0, 40.0}});
      histos.add("REC_True_Prompt_Trigger_Energy", "REC_True_Prompt_Trigger_Energy", kTH1F, {{82, -1.0, 40.0}});
      histos.add("REC_Trigger_V_PtHadSum_Stern", "REC_Trigger_V_PtHadSum_Stern", kTH2F, {{100, 0, 100}, pthadAxis});
      histos.add("REC_Trigger_V_PtHadSum_Nch", "REC_Trigger_V_PtHadSum_Nch", kTH2F, {{100, 0, 100}, pthadAxis});
      histos.add("REC_Trigger_V_PtHadSum_Photon", "REC_Trigger_V_PtHadSum_Photon", kTH2F, {{100, 0, 100}, pthadAxis});
      histos.add("REC_TrueTrigger_V_PtHadSum_Photon", "REC_Trigger_V_PtHadSum_Photon", kTH2F, {{100, 0, 100}, pthadAxis});
      histos.add("REC_dR_Photon", "REC_dR_Photon", kTH1F, {{628, 0.0, 2 * TMath::Pi()}});
      histos.add("REC_dR_Stern", "REC_dR_Stern", kTH1F, {{628, 0.0, 2 * TMath::Pi()}});
    }
    if (cfgGenHistograms) {
      histos.add("GEN_nEvents", "GEN_nEvents", kTH1F, {{4, 0.0, 4.0}});
      histos.add("GEN_True_Trigger_Energy", "GEN_True_Trigger_Energy", kTH1F, {{82, -1.0, 40.0}});
      histos.add("GEN_Particle_Pt", "GEN_Particle_Pt", kTH1F, {{82, -1.0, 40.0}});
      histos.add("GEN_True_Photon_Energy", "GEN_True_Photon_Energy", kTH1F, {{8200, -1.0, 40.0}});
      histos.add("GEN_True_Prompt_Photon_Energy", "GEN_True_Prompt_Photon_Energy", kTH1F, {{8200, -1.0, 40.0}});
      histos.add("GEN_Trigger_V_PtHadSum_Stern", "GEN_Trigger_V_PtHadSum_Stern", kTH2F, {{100, 0, 100}, pthadAxis});
      histos.add("GEN_Trigger_V_PtHadSum_Photon", "GEN_Trigger_V_PtHadSum_Photon", kTH2F, {{100, 0, 100}, pthadAxis});
      histos.add("GEN_TrueTrigger_V_PtHadSum_Photon", "GEN_Trigger_V_PtHadSum_Photon", kTH2F, {{100, 0, 100}, pthadAxis});
      histos.add("GEN_dR_Photon", "GEN_dR_Photon", kTH1F, {{628, 0.0, 2 * TMath::Pi()}});
      histos.add("GEN_dR_Stern", "GEN_dR_Stern", kTH1F, {{628, 0.0, 2 * TMath::Pi()}});
    }
    if (cfgDataHistograms) {
      histos.add("DATA_nEvents", "DATA_nEvents", kTH1F, {{4, 0.0, 4.0}});
      histos.add("DATA_M02_BC", "DATA_M02_BC", kTH1F, {{628, -2 * TMath::Pi(), 2 * TMath::Pi()}});
      histos.add("DATA_M02_AC", "DATA_M02_AC", kTH1F, {{628, -2 * TMath::Pi(), 2 * TMath::Pi()}});
      histos.add("DATA_Cluster_QA", "DATA_Cluster_QA", kTH1F, {{10, -0.5, 9.5}});
      histos.add("DATA_ClusterPhi", "DATA_ClusterPhi", kTH1F, {{640 * 2, 0, 2 * TMath::Pi()}});
      histos.add("DATA_ClusterEta", "DATA_ClusterEta", kTH1F, {{100, -1, 1}});
      histos.add("DATA_All_Energy", "DATA_All_Energy", kTH1F, {{82, -1.0, 40.0}});
      histos.add("DATA_Cluster_PhiPrime_Pt", "DATA_Cluster_PhiPrime_Pt", kTH2F, {{640, 0, 2 * TMath::Pi()}, {82, -1.0, 40.0}});
      histos.add("DATA_Cluster_PhiPrime_Pt_AC", "DATA_Cluster_PhiPrime_Pt_AC", kTH2F, {{640, 0, 2 * TMath::Pi()}, {82, -1.0, 40.0}});
      histos.add("DATA_Cluster_PhiPrime_Pt_C", "DATA_Cluster_PhiPrime_Pt_C", kTH2F, {{640, 0, 2 * TMath::Pi()}, {82, -1.0, 40.0}});
      histos.add("DATA_Track_v_Cluster_Phi", "DATA_Track_v_Cluster_Phi", kTH1F, {{628, -2 * TMath::Pi(), 2 * TMath::Pi()}});
      histos.add("DATA_Track_v_Cluster_Eta", "DATA_Track_v_Cluster_Eta", kTH1F, {{628, -2 * TMath::Pi(), 2 * TMath::Pi()}});
      histos.add("DATA_Track_v_Cluster_Phi_Pt", "DATA_Track_v_Cluster_Phi_Pt", kTH2F, {{628, -2 * TMath::Pi(), 2 * TMath::Pi()}, {82, -1.0, 40.0}});
      histos.add("DATA_Track_v_Cluster_Phi_Eta", "DATA_Track_v_Cluster_Phi_Eta", kTH2F, {{628, -2 * TMath::Pi(), 2 * TMath::Pi()}, {628, -2 * TMath::Pi(), 2 * TMath::Pi()}});
      histos.add("DATA_SumPt_BC", "DATA_SumPt_BC", kTH1F, {{628, -2 * TMath::Pi(), 2 * TMath::Pi()}});
      histos.add("DATA_SumPt_AC", "DATA_SumPt_AC", kTH1F, {{628, -2 * TMath::Pi(), 2 * TMath::Pi()}});
      histos.add("DATA_Trigger_V_PtHadSum_Photon", "DATA_Trigger_V_PtHadSum_Photon", kTH2F, {{100, 0, 100}, pthadAxis});
      histos.add("DATA_PtHadSum_Photon", "DATA_PtHadSum_Photon", kTH1F, {pthadAxis});
      histos.add("DATA_Trigger_Energy", "DATA_Trigger_Energy", kTH1F, {{82, -1.0, 40.0}});
      histos.add("DATA_Track_PhiPrime_Pt", "DATA_Track_PhiPrime_Pt", kTH2F, {{640, 0, 2 * TMath::Pi()}, {82, -1.0, 40.0}});
      histos.add("DATA_Track_Pt", "DATA_Track_Pt", kTH1F, {{82, -1.0, 40.0}});
      histos.add("DATA_Track_Phi", "DATA_Track_Phi", kTH1F, {{640 * 2, 0, 2 * TMath::Pi()}});
      histos.add("DATA_TrackPhi_photontrigger", "DATA_TrackPhi_photontrigger", kTH1F, {{64, 0, 2 * TMath::Pi()}});
      histos.add("DATA_TrackEta_photontrigger", "DATA_TrackEta_photontrigger", kTH1F, {{100, -1, 1}});
      histos.add("DATA_Trigger_V_PtHadSum_Stern", "DATA_Trigger_V_PtHadSum_Stern", kTH2F, {{100, 0, 100}, pthadAxis});
      histos.add("DATA_Trigger_V_PtHadSum_Nch", "DATA_Trigger_V_PtHadSum_Nch", kTH2F, {{100, 0, 100}, pthadAxis});
      histos.add("DATA_dR_Photon", "DATA_dR_Photon", kTH1F, {{628, 0.0, 2 * TMath::Pi()}});
      histos.add("DATA_dR_Stern", "DATA_dR_Stern", kTH1F, {{628, 0.0, 2 * TMath::Pi()}});
    }
  } // end of init

  Filter PosZFilter_JE = nabs(aod::jcollision::posZ) < cfgVtxCut;
  Filter mcPosZFilter = nabs(aod::jmccollision::posZ) < cfgVtxCut;
  Filter clusterDefinitionSelection_JE = (o2::aod::jcluster::definition == cfgClusterDefinition) && (o2::aod::jcluster::time >= cfgMinTime) && (o2::aod::jcluster::time <= cfgMaxTime) && (o2::aod::jcluster::energy > cfgMinClusterEnergy) && (o2::aod::jcluster::nCells >= cfgMinNCells) && (o2::aod::jcluster::nlm <= cfgMaxNLM) && (o2::aod::jcluster::isExotic == cfgExoticContribution);

  using selectedMCCollisions = aod::JMcCollisions;
  using filteredMCCollisions = soa::Filtered<selectedMCCollisions>;
  using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>;
  using BcCandidates = soa::Join<aod::JBCs, aod::JBCPIs>;

  using jTrackCandidates = soa::Join<aod::JTracks, aod::JTrackPIs, aod::JMcTrackLbs>;
  using jEMCtracks = aod::JEMCTracks;

  using jDataTrackCandidates = soa::Join<aod::JTracks, aod::JTrackPIs>;

  using jMCClusters = o2::soa::Join<o2::aod::JMcClusterLbs, o2::aod::JClusters, o2::aod::JClusterTracks>;
  using jClusters = o2::soa::Join<o2::aod::JClusters, o2::aod::JClusterTracks>;
  using jselectedCollisions = soa::Join<aod::JCollisions, aod::JCollisionBCs, aod::JCollisionPIs, aod::EvSels, aod::JEMCCollisionLbs, aod::JMcCollisionLbs>;
  using jselectedDataCollisions = soa::Join<aod::JCollisions, aod::JCollisionBCs, aod::JCollisionPIs, aod::EvSels, aod::JEMCCollisionLbs>;
  //  using jselectedDataCollisions = soa::Join<aod::JCollisions, aod::JCollisionBCs, aod::JCollisionPIs, aod::JCollisionMcInfos, aod::EvSels, aod::JEMCCollisionLbs>;
  using jfilteredCollisions = soa::Filtered<jselectedCollisions>;
  using jfilteredDataCollisions = soa::Filtered<jselectedDataCollisions>;
  using jfilteredMCClusters = soa::Filtered<jMCClusters>;
  using jfilteredClusters = soa::Filtered<jClusters>;

  Preslice<o2::aod::EMCALMatchedTracks> perClusterMatchedTracks = o2::aod::emcalclustercell::emcalclusterId;

  // Helper functions
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <typename Tracks, typename Trigger>

  double GetPtHadSum(const Tracks& tracks, const Trigger& trigger, double MinR, double MaxR, bool IsStern, bool IsParticle, bool DodR)
  {
    double eta_trigger, phi_trigger;
    if constexpr (requires { trigger.eta(); }) {
      eta_trigger = trigger.eta();
      phi_trigger = trigger.phi();
    } else if constexpr (requires { trigger.Eta(); }) {
      eta_trigger = trigger.Eta();
      phi_trigger = trigger.Phi();
    }
    double pthadsum = 0;

    for (auto& track : tracks) {
      double phi_track = track.phi();
      double eta_track = track.eta();
      double pt_track = track.pt();

      if (!IsParticle) {
        if constexpr (requires { track.trackId(); }) {
          if (cfgJETracks) {
            if (!jetderiveddatautilities::selectTrack(track, trackFilter)) {
              continue;
            }
          } else {
            auto originaltrack = track.template track_as<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>>();
            if (!trackSelection(originaltrack)) {
              continue;
            }
          } // reject track
        } else if constexpr (requires { track.sign(); }) { // checking for JTrack
          if (cfgJETracks) {
            if (!jetderiveddatautilities::selectTrack(track, trackFilter)) {
              continue;
            }
          } else {
            if (!trackSelection(track)) {
              continue;
            }
          } // reject track
        } // done checking for JTrack
      } else {
        if constexpr (requires { track.isPhysicalPrimary(); }) {
          if (track.pt() < 0.15) {
            continue;
          }
          if (std::abs(track.eta()) > cfgtrkMaxEta) {
            continue;
          }
          if (track.getGenStatusCode() < 20) {
            continue;
          }
          if (!track.isPhysicalPrimary()) {
            continue;
          }
          int pdg = std::abs(track.pdgCode());
          if (pdg != 211 && pdg != 321 && pdg != 2212) {
            continue;
          }
        }
      }
      if (IsStern || IsParticle) {
        if constexpr (requires { trigger.globalIndex(); }) {
          if (trigger.globalIndex() == track.globalIndex())
            continue;
        }
      }
      double phidiff = TVector2::Phi_mpi_pi(phi_track - phi_trigger);
      double etadiff = std::abs(eta_track - eta_trigger);
      double dR = TMath::Sqrt((etadiff * etadiff) + (phidiff * phidiff));

      if (DodR) {
        if (dR > MinR && dR < MaxR) {
          if (!IsParticle) {
            if (cfgRecHistograms) {
              if (IsStern) {
                histos.fill(HIST("REC_dR_Stern"), dR);
              }
              if (!IsStern) {
                histos.fill(HIST("REC_dR_Photon"), dR);
              }
            } else if (cfgDataHistograms) {
              if (IsStern) {
                histos.fill(HIST("DATA_dR_Stern"), dR);
              }
              if (!IsStern) {
                histos.fill(HIST("DATA_dR_Photon"), dR);
              }
            }
          } else {
            if (IsStern) {
              histos.fill(HIST("GEN_dR_Stern"), dR);
            }
            if (!IsStern) {
              histos.fill(HIST("GEN_dR_Photon"), dR);
            }
          }
        }
      }
      if (dR > MinR && dR < MaxR) {
        pthadsum += pt_track;
      }

    } // end of track loop
    return pthadsum;
  } // end of GetPtHadSum
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <typename TrackType>
  bool trackSelection(const TrackType track)
  {
    // basic track cuts
    if (track.pt() < cfgtrkMinPt)
      return false;

    if (std::abs(track.eta()) > cfgtrkMaxEta)
      return false;

    if (std::abs(track.dcaXY()) > cfgMaxDCArToPVcut)
      return false;

    if (std::abs(track.dcaZ()) > cfgMaxDCAzToPVcut)
      return false;

    if (cfgPrimaryTrack && !track.isPrimaryTrack())
      return false;

    if (track.tpcNClsFindable() < cfgnFindableTPCClusters)
      return false;

    if (track.tpcNClsCrossedRows() < cfgnTPCCrossedRows)
      return false;

    if (track.tpcCrossedRowsOverFindableCls() > cfgnRowsOverFindable)
      return false;

    if (track.tpcChi2NCl() > cfgnTPCChi2)
      return false;

    if (track.itsChi2NCl() > cfgnITSChi2)
      return false;

    if (cfgConnectedToPV && !track.isPVContributor())
      return false;

    return true;
  }; // end of track selection
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////
  // PROCESS
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  int nEventsGenMC = 0;

  // void processMCGen(filteredMCCollisions::iterator const& collision, soa::SmallGroups<soa::Join<aod::McCollisionLabels, filteredCollisions>> const& recocolls, aod::McParticles const& mcParticles, filteredMCClusters const&)
  void processMCGen(filteredMCCollisions::iterator const& collision, soa::SmallGroups<soa::Join<aod::JMcCollisionLbs, jfilteredCollisions>> const& recocolls, aod::JMcParticles const& mcParticles, jfilteredMCClusters const&)
  {
    nEventsGenMC++;
    if (cfgDebug) {
      if ((nEventsGenMC + 1) % 10000 == 0) {
        std::cout << "Processed Gen MC Events: " << nEventsGenMC << std::endl;
      }
    }
    histos.fill(HIST("GEN_nEvents"), 0.5);
    if (fabs(collision.posZ()) > cfgVtxCut)
      return;

    if (recocolls.size() <= 0) // not reconstructed
      return;
    for (auto& recocoll : recocolls) { // poorly reconstructed
      if (!recocoll.sel8())
        return;
      if (fabs(recocoll.posZ()) > cfgVtxCut)

        return;
      histos.fill(HIST("GEN_nEvents"), 1.5);

      if (cfgEmcTrigger) {
        if (!recocoll.isEmcalReadout())
          return;
      }
      histos.fill(HIST("GEN_nEvents"), 2.5);
    }

    for (auto& mcPhoton : mcParticles) {
      bool photontrigger = false;
      if (mcPhoton.pt() < 0.15)
        continue;
      if (std::abs(mcPhoton.eta()) > cfgtrkMaxEta)
        continue;
      double pdgcode = fabs(mcPhoton.pdgCode());
      if (mcPhoton.isPhysicalPrimary()) {
        if (pdgcode == 211 || pdgcode == 321 || pdgcode == 2212 || pdgcode == 11) {
          histos.fill(HIST("GEN_Particle_Pt"), mcPhoton.pt());
        }
      }
      if (mcPhoton.getGenStatusCode() < 20)
        continue;

      // first we check for pthadsums for all charged particles a la sternheimer
      if (mcPhoton.isPhysicalPrimary()) {
        int pdg = std::abs(mcPhoton.pdgCode());
        if (pdg == 211 || pdg == 321 || pdg == 2212) {
          bool sterntrigger = false;
          double sternPt = 0.0;
          if (mcPhoton.pt() > cfgMinTrig && mcPhoton.pt() < cfgMaxTrig) {
            if (fabs(mcPhoton.eta()) <= cfgtrkMaxEta) {
              sterntrigger = true;
              sternPt = mcPhoton.pt();
            }
          }
          // stern trigger
          if (sterntrigger) {
            bool doStern = true;
            double sterncount = 1.0;
            while (doStern) {
              TLorentzVector lParticleTrigger;
              lParticleTrigger.SetPxPyPzE(mcPhoton.px(), mcPhoton.py(), mcPhoton.pz(), mcPhoton.e());
              double pthadsum = GetPtHadSum(mcParticles, lParticleTrigger, cfgMinR, cfgMaxR, true, true, true);
              histos.fill(HIST("GEN_Trigger_V_PtHadSum_Stern"), sterncount, pthadsum, 2.0 / sternPt);
              if (sterncount < sternPt) {
                sterncount++;
              } else {
                doStern = false;
              }
            } // While sternin'
          } // stern trigger loop
        } // check if charged pikp
      } // check for primary particles

      // now we do all photons
      if (mcPhoton.pdgCode() == 22) {
        histos.fill(HIST("GEN_True_Photon_Energy"), mcPhoton.e());
        if (mcPhoton.pt() > cfgMinTrig && mcPhoton.pt() < cfgMaxTrig) {
          if (fabs(mcPhoton.eta()) <= cfgtrkMaxEta) {
            photontrigger = true;
          }
        } // check for photon trigger
        if (photontrigger) {
          TLorentzVector lRealPhoton;
          lRealPhoton.SetPxPyPzE(mcPhoton.px(), mcPhoton.py(), mcPhoton.pz(), mcPhoton.e());
          double truepthadsum = GetPtHadSum(mcParticles, lRealPhoton, cfgMinR, cfgMaxR, false, true, false);
          histos.fill(HIST("GEN_Trigger_V_PtHadSum_Photon"), mcPhoton.e(), truepthadsum);
        }

        // now we do all PROMPT photons
        histos.fill(HIST("GEN_True_Trigger_Energy"), mcPhoton.e());

        int mompdg1 = 0;
        int momindex1 = 0;
        int momstatus1 = 0;
        for (auto& photon_mom : mcPhoton.mothers_as<aod::JMcParticles>()) {
          if (mompdg1 == 0) {
            mompdg1 = photon_mom.pdgCode();
            momindex1 = photon_mom.globalIndex();
            momstatus1 = photon_mom.getGenStatusCode();
          }
        } // first photon loop

        if (std::fabs(mompdg1) < 40 && std::fabs(mompdg1) > 0) {
          int mompdg2 = 0;
          int momindex2 = 0;
          int momstatus2 = 0;
          int mompdg3 = 0;
          int momindex3 = 0;
          int momstatus3 = 0;
          for (auto& mcPhoton_mom : mcParticles) {
            if (mcPhoton_mom.globalIndex() == momindex1) {
              for (auto& photon_momom : mcPhoton_mom.mothers_as<aod::JMcParticles>()) {
                if (mompdg2 == 0) {
                  mompdg2 = photon_momom.pdgCode();
                  momindex2 = photon_momom.globalIndex();
                  momstatus2 = photon_momom.getGenStatusCode();
                }
              }
              break;
            }
          } // 2nd photon loop
          if (std::fabs(mompdg2) < 40 && std::fabs(mompdg2) > 0) {
            for (auto& mcPhoton_mom : mcParticles) {
              if (mcPhoton_mom.globalIndex() == momindex2) {
                for (auto& photon_momom : mcPhoton_mom.mothers_as<aod::JMcParticles>()) {
                  if (mompdg3 == 0) {
                    mompdg3 = photon_momom.pdgCode();
                    momindex3 = photon_momom.globalIndex();
                    momstatus3 = photon_momom.getGenStatusCode();
                  }
                }
                break;
              }
            } // 3rd photon loop
          } // 2nd photon check

          if (cfgDebug) {
            std::cout << "We have a GEN prompt photon" << std::endl;
            std::cout << "Photon gen status code chain: " << std::endl;
            std::cout << "Photon stat: " << mcPhoton.getGenStatusCode() << std::endl;
            std::cout << "Photon index: " << mcPhoton.globalIndex() << std::endl;
            std::cout << "Photon mompdg 1: " << mompdg1 << std::endl;
            std::cout << "Photon momstatus 1: " << momstatus1 << std::endl;
            std::cout << "Photon momindex 1: " << momindex1 << std::endl;
            std::cout << "Photon mompdg 2: " << mompdg2 << std::endl;
            std::cout << "Photon momstatus 2: " << momstatus2 << std::endl;
            std::cout << "Photon momindex 2: " << momindex2 << std::endl;
            std::cout << "Photon mompdg 3: " << mompdg3 << std::endl;
            std::cout << "Photon momstatus 3: " << momstatus3 << std::endl;
            std::cout << "Photon momindex 3: " << momindex3 << std::endl;
          }
        } else {
          continue;
        }

        if (std::abs(mcPhoton.getGenStatusCode()) > 19 && std::abs(mcPhoton.getGenStatusCode()) < 90) {
          if (mcPhoton.isPhysicalPrimary()) {
            histos.fill(HIST("GEN_True_Prompt_Photon_Energy"), mcPhoton.e());
            if (photontrigger) {
              TLorentzVector lRealPromptPhoton;
              lRealPromptPhoton.SetPxPyPzE(mcPhoton.px(), mcPhoton.py(), mcPhoton.pz(), mcPhoton.e());
              double truepthadsum = GetPtHadSum(mcParticles, lRealPromptPhoton, cfgMinR, cfgMaxR, false, true, true);
              histos.fill(HIST("GEN_TrueTrigger_V_PtHadSum_Photon"), mcPhoton.e(), truepthadsum);
            } // photontrigger
          } // check for primary photons
        } // prompt photon check

      } // photon check

    } // loop over mc particles

  } // end of process

  PROCESS_SWITCH(statPromptPhoton, processMCGen, "process MC Gen", true);

  PresliceUnsorted<jEMCtracks> EMCTrackPerTrack = aod::jemctrack::trackId;
  int nEventsRecMC_JE = 0;
  void processMCRec_JE(jfilteredCollisions::iterator const& collision, jfilteredMCClusters const& mcclusters, jTrackCandidates const& tracks, soa::Join<aod::JTracks, aod::JTrackExtras, aod::JTrackPIs> const&, TrackCandidates const&, aod::JMcParticles const&, BcCandidates const&, jEMCtracks const& emctracks, aod::JetMcCollisions const&)
  {

    nEventsRecMC_JE++;
    if (cfgDebug) {
      if ((nEventsRecMC_JE + 1) % 10000 == 0) {
        std::cout << "Processed JE Rec MC Events: " << nEventsRecMC_JE << std::endl;
      }
    }
    histos.fill(HIST("REC_nEvents"), 0.5);

    // required cuts
    if (fabs(collision.posZ()) > cfgVtxCut)
      return;
    if (!collision.sel8())
      return;

    histos.fill(HIST("REC_nEvents"), 1.5);

    if (cfgEmcTrigger) {
      if (!collision.isEmcalReadout())
        return;
    }
    histos.fill(HIST("REC_nEvents"), 2.5);

    if (!jetderiveddatautilities::selectTrigger(collision, triggerMaskBits)) {
      return;
    }

    histos.fill(HIST("REC_nEvents"), 3.5);

    double weight = 1;
    if (collision.has_mcCollision()) {
      weight = collision.mcCollision().weight();
    }

    bool noTrk = true;
    for (auto& track : tracks) {
      if (cfgJETracks) {
        if (!jetderiveddatautilities::selectTrack(track, trackFilter)) {
          continue;
        }
      } else {
        auto ogtrack = track.track_as<TrackCandidates>();
        if (!trackSelection(ogtrack)) {
          continue;
        }
        if (!ogtrack.isGlobalTrack()) {
          continue;
        }
      }
      noTrk = false;
      break;
    }

    if (noTrk)
      return;

    // now we do clusters
    bool clustertrigger = false;
    for (auto& mccluster : mcclusters) {
      histos.fill(HIST("REC_M02_BC"), mccluster.m02());
      if (mccluster.m02() < cfgLowM02)
        continue;
      if (mccluster.m02() > cfgHighM02)
        continue;
      if (mccluster.energy() < cfgLowClusterE)
        continue;
      if (mccluster.energy() > cfgHighClusterE)
        continue;
      if (fabs(mccluster.eta()) > cfgtrkMaxEta)
        continue;

      histos.fill(HIST("REC_Cluster_QA"), 0.5);
      histos.fill(HIST("REC_M02_AC"), mccluster.m02());
      histos.fill(HIST("REC_All_Energy"), mccluster.energy());
      histos.fill(HIST("REC_ClusterPhi"), mccluster.phi());
      histos.fill(HIST("REC_ClusterEta"), mccluster.eta());

      bool photontrigger = false; // is a neutral cluster
      bool chargetrigger = false; // is definitely not a neutral cluster
      // double photonPt = 0.0;
      double truephotonPt = 0.0;
      auto tracksofcluster = mccluster.matchedTracks_as<soa::Join<aod::JTracks, aod::JTrackExtras, aod::JTrackPIs>>();
      ///////////////
      ///////////////
      ///////////////

      double phiPrimeC = mccluster.phi();
      phiPrimeC = phiPrimeC + TMath::Pi() / 18.;
      phiPrimeC = fmod(phiPrimeC, 2 * TMath::Pi() / 18.);
      double ptC = mccluster.energy();
      bool geocut = false;
      histos.fill(HIST("REC_Cluster_PhiPrime_Pt"), phiPrimeC, mccluster.energy());
      if (phiPrimeC > (0.12 / ptC + TMath::Pi() / 18. + 0.035) ||
          phiPrimeC < (0.1 / ptC / ptC + TMath::Pi() / 18. - 0.025)) {
        geocut = false;
      } else {
        geocut = true;
      }

      if (cfgGeoCut) {
        if (geocut) {
          histos.fill(HIST("REC_Cluster_PhiPrime_Pt_C"), phiPrimeC, mccluster.energy());
          continue;
        }
      }

      histos.fill(HIST("REC_Cluster_PhiPrime_Pt_AC"), phiPrimeC, mccluster.energy());

      // first, we check if veto is required
      double sumptT = 0;
      bool clusterqa = false;
      for (auto& ctrack : tracksofcluster) {
        double etaT, phiT;

        if (cfgJETracks) {
          if (!jetderiveddatautilities::selectTrack(ctrack, trackFilter)) {
            continue;
          }
          auto emctracksPerTrack = emctracks.sliceBy(EMCTrackPerTrack, ctrack.globalIndex());
          auto emctrack = emctracksPerTrack.iteratorAt(0);
          etaT = emctrack.etaEmcal();
          phiT = emctrack.phiEmcal();
        } else {
          auto ogtrack = ctrack.track_as<TrackCandidates>();
          if (!trackSelection(ogtrack)) {
            continue;
          }
          if (!ogtrack.isGlobalTrack()) {
            continue;
          }
          etaT = ogtrack.trackEtaEmcal();
          phiT = ogtrack.trackPhiEmcal();
        }

        double etaC = mccluster.eta();
        double phiC = mccluster.phi();
        double ptT = ctrack.pt();
        bool etatrigger = false;
        bool phitrigger = false;
        double phidiff = TVector2::Phi_mpi_pi(mccluster.phi() - ctrack.phi());
        double etadiff = mccluster.eta() - ctrack.eta();

        if (cfgPtClusterCut) {
          if (fabs(etaT - etaC) < (0.010 + pow(ptT + 4.07, -2.5))) {
            etatrigger = true;
          }

          if (fabs(TVector2::Phi_mpi_pi(phiT - phiC)) < (0.015 + pow(ptT + 3.65, -2.0))) {
            phitrigger = true;
          }
        } else {
          if (fabs(etadiff) < 0.05) {
            etatrigger = true;
          }

          if (fabs(phidiff) < 0.05) {
            phitrigger = true;
          }
        }

        if (etatrigger && phitrigger) {
          chargetrigger = true;
          sumptT += ptT;
        }
        if (chargetrigger) {
          if (!clusterqa) {
            histos.fill(HIST("REC_Cluster_QA"), 1.5);
            clusterqa = true;
          }
        }
        histos.fill(HIST("REC_Track_v_Cluster_Phi"), phidiff);
        histos.fill(HIST("REC_Track_v_Cluster_Eta"), etadiff);

        histos.fill(HIST("REC_Track_v_Cluster_Phi_Eta"), phidiff, etadiff);

      } // track of cluster loop

      if (chargetrigger && sumptT > 0) {
        double mccluster_over_sumptT = mccluster.energy() / sumptT;
        histos.fill(HIST("REC_SumPt_BC"), mccluster_over_sumptT);
        if (mccluster_over_sumptT < 1.7) {
          histos.fill(HIST("REC_Cluster_QA"), 2.5); // veto fails, cluster is charged
        } else {
          histos.fill(HIST("REC_Cluster_QA"), 3.5); // veto is good, cluster is converted to neutral cluster
          // chargetrigger = false;
          histos.fill(HIST("REC_SumPt_AC"), mccluster_over_sumptT);
        }
      } // sumptT check

      if (!chargetrigger) {
        photontrigger = true;
      }

      ///////////////
      ///////////////
      ///////////////

      // check if cluster is good
      for (auto& ctrack : tracksofcluster) {
        if (cfgJETracks) {
          if (!jetderiveddatautilities::selectTrack(ctrack, trackFilter)) {
            continue;
          }
        } else {
          auto ogtrack = ctrack.track_as<TrackCandidates>();
          if (!trackSelection(ogtrack)) {
            continue;
          }
          if (!ogtrack.isGlobalTrack()) {
            continue;
          }
        }
        bool etatrigger = false;
        bool phitrigger = false;
        // double ptT = ctrack.pt();
        double phidiff = TVector2::Phi_mpi_pi(mccluster.phi() - ctrack.phi());
        double etadiff = mccluster.eta() - ctrack.eta();
        if (fabs(etadiff) < 0.05) {
          etatrigger = true;
        }

        if (fabs(phidiff) < 0.05) {
          phitrigger = true;
        }

        if (chargetrigger) {
          histos.fill(HIST("REC_Track_v_Cluster_Phi_C"), phidiff);
          histos.fill(HIST("REC_Track_v_Cluster_Eta_C"), etadiff);
          histos.fill(HIST("REC_Track_v_Cluster_Phi_Eta_C"), phidiff, etadiff);
        } else {
          if ((etatrigger || phitrigger) && chargetrigger) {
            if (cfgDebug) {
              std::cout << "????????????????????" << std::endl;
            }
          }
          histos.fill(HIST("REC_Track_v_Cluster_Phi_AC"), phidiff);
          histos.fill(HIST("REC_Track_v_Cluster_Eta_AC"), etadiff);
          histos.fill(HIST("REC_Track_v_Cluster_Phi_Eta_AC"), phidiff, etadiff);
        }
      } // tracks

      ///////////////
      ///////////////
      ///////////////

      if (photontrigger) { // if no charge trigger, cluster is good!
        histos.fill(HIST("REC_Cluster_QA"), 4.5);
        clustertrigger = true;
        double pthadsum = GetPtHadSum(tracks, mccluster, cfgMinR, cfgMaxR, false, false, true);
        histos.fill(HIST("REC_Trigger_V_PtHadSum_Photon"), mccluster.energy(), pthadsum, weight);
        histos.fill(HIST("REC_PtHadSum_Photon"), pthadsum, weight);
        histos.fill(HIST("REC_Trigger_Energy"), mccluster.energy(), weight);
      }

      auto ClusterParticles = mccluster.mcParticles_as<aod::JMcParticles>();

      // now we check the realness of our prompt photons
      bool goodgentrigger = true;
      double chPe = 0;
      for (auto& clusterparticle : ClusterParticles) {
        int cindex = clusterparticle.globalIndex();
        double pdgcode = fabs(clusterparticle.pdgCode());
        if (!clusterparticle.isPhysicalPrimary()) {
          continue;
        }
        if (pdgcode == 211 || pdgcode == 321 || pdgcode == 2212 || pdgcode == 11) {
          bool notrack = true;
          histos.fill(HIST("REC_Cluster_Particle_Pt"), clusterparticle.pt());

          double phiPrimeP = clusterparticle.phi();
          if (clusterparticle.pdgCode() < 0) {
            phiPrimeP = 2 * TMath::Pi() - phiPrimeP;
          }
          phiPrimeP = phiPrimeP + TMath::Pi() / 18.;
          phiPrimeP = fmod(phiPrimeP, 2 * TMath::Pi() / 18.);
          double ptP = clusterparticle.pt();
          for (auto& track : tracks) {
            if (!track.has_mcParticle())
              continue;

            if (cfgJETracks) {
              if (!jetderiveddatautilities::selectTrack(track, trackFilter)) {
                continue;
              }
            } else {
              auto ogtrack = track.track_as<TrackCandidates>();
              if (!trackSelection(ogtrack)) {
                continue;
              }
              if (!ogtrack.isGlobalTrack()) {
                continue;
              }
            }

            int tindex = track.mcParticleId();
            if (tindex == cindex) {
              histos.fill(HIST("REC_Cluster_ParticleWITHtrack_Pt"), clusterparticle.pt());
              histos.fill(HIST("REC_Cluster_ParticleWITHtrack_Phi"), clusterparticle.phi());
              histos.fill(HIST("REC_Cluster_ParticleWITHtrack_Eta"), clusterparticle.eta());
              histos.fill(HIST("REC_Cluster_ParticleWITHtrack_Pt_Phi"), clusterparticle.pt(), clusterparticle.phi());
              histos.fill(HIST("REC_Cluster_ParticleWITHtrack_Pt_PhiPrime"), ptP, phiPrimeP);
              if (photontrigger) {
                histos.fill(HIST("REC_Impurity_ParticleWITHtrack_Pt_PhiPrime"), ptP, phiPrimeP);
              }
              // }//geo cut
              histos.fill(HIST("REC_Cluster_ParticleWITHtrack_Pt_Eta"), clusterparticle.pt(), clusterparticle.eta());

              histos.fill(HIST("REC_Cluster_ParticleWITHtrack_TrackPt"), track.pt());
              notrack = false;
              break;
            }
          } // track loop

          if (notrack) {
            histos.fill(HIST("REC_Cluster_ParticleWITHOUTtrack_Pt"), clusterparticle.pt());
            histos.fill(HIST("REC_Cluster_ParticleWITHOUTtrack_Phi"), clusterparticle.phi());
            histos.fill(HIST("REC_Cluster_ParticleWITHOUTtrack_Eta"), clusterparticle.eta());
            histos.fill(HIST("REC_Cluster_ParticleWITHOUTtrack_Pt_Phi"), clusterparticle.pt(), clusterparticle.phi());
            histos.fill(HIST("REC_Cluster_ParticleWITHOUTtrack_Pt_PhiPrime"), ptP, phiPrimeP);
            if (photontrigger) {
              histos.fill(HIST("REC_Impurity_ParticleWITHOUTtrack_Pt_PhiPrime"), ptP, phiPrimeP);
            }
            //  }//geo cut
            histos.fill(HIST("REC_Cluster_ParticleWITHOUTtrack_Pt_Eta"), clusterparticle.pt(), clusterparticle.eta());
          }
        } // pdg code check

        double phidiff = TVector2::Phi_mpi_pi(mccluster.phi() - clusterparticle.phi());
        double etadiff = mccluster.eta() - clusterparticle.eta();

        if (pdgcode == 211 || pdgcode == 321 || pdgcode == 2212 || pdgcode == 11) {
          if (clusterparticle.e() > 0.01) {
            chPe += clusterparticle.e();
            goodgentrigger = false;
            if (photontrigger) {
              histos.fill(HIST("REC_Impurity_Energy_v_Cluster_Phi"), clusterparticle.e(), phidiff);
              histos.fill(HIST("REC_Impurity_Energy_v_ClusterE_Phi"), mccluster.energy(), phidiff);
              histos.fill(HIST("REC_Impurity_Energy_v_ClusterEoP_Phi"), mccluster.energy() / clusterparticle.e(), phidiff);
              histos.fill(HIST("REC_Impurity_Energy_v_Cluster_Energy"), mccluster.energy(), clusterparticle.e());

              histos.fill(HIST("REC_TrueImpurity_v_Cluster_Phi"), phidiff);
              histos.fill(HIST("REC_TrueImpurity_v_Cluster_PhiAbs"), clusterparticle.phi());
              histos.fill(HIST("REC_TrueImpurity_v_Cluster_Eta"), etadiff);
              histos.fill(HIST("REC_TrueImpurity_v_Cluster_EtaAbs"), clusterparticle.eta());
              histos.fill(HIST("REC_TrueImpurity_v_Cluster_Phi_Eta"), phidiff, etadiff);
            }
          }
        }
        histos.fill(HIST("REC_True_v_Cluster_Phi"), phidiff);
        histos.fill(HIST("REC_True_v_Cluster_Eta"), etadiff);

        if (!photontrigger) {
          continue;
        }
        if (clusterparticle.pdgCode() == 22) {
          histos.fill(HIST("REC_True_Trigger_Energy"), clusterparticle.e());
          int mom1 = 0;
          int mom2 = 0;
          for (auto& photon_mom : clusterparticle.mothers_as<aod::JMcParticles>()) {
            if (mom1 == 0) {
              mom1 = photon_mom.pdgCode();
            }
            if (mom1 != 0) {
              mom2 = photon_mom.pdgCode();
            }
          }
          if (std::fabs(mom1) > 40 && std::fabs(mom1) > 0)
            continue;

          if (cfgDebug) {
            std::cout << "We have a REC prompt photon" << std::endl;
            std::cout << "Photon gen status code: " << clusterparticle.getGenStatusCode() << std::endl;
            std::cout << "Photon mom 1: " << mom1 << std::endl;
            std::cout << "Photon mom 2: " << mom2 << std::endl;
          }
          if (std::abs(clusterparticle.getGenStatusCode()) > 19 && std::abs(clusterparticle.getGenStatusCode()) < 90) {
            histos.fill(HIST("REC_True_Prompt_Trigger_Energy"), clusterparticle.e(), weight);
            TLorentzVector lRealPhoton;
            lRealPhoton.SetPxPyPzE(clusterparticle.px(), clusterparticle.py(), clusterparticle.pz(), clusterparticle.e());
            double truepthadsum = GetPtHadSum(tracks, lRealPhoton, cfgMinR, cfgMaxR, false, false, false);
            truephotonPt = clusterparticle.e();
            histos.fill(HIST("REC_TrueTrigger_V_PtHadSum_Photon"), truephotonPt, truepthadsum, weight);
          }
        } // photon check
      } // clusterparticle loop

      if (cfgDebug) {
        if (chPe > 0) {
          if (photontrigger) {
            if (chPe / mccluster.energy() < 0.50) {
              goodgentrigger = true;
            }
          }
        }
      }
      if (goodgentrigger && photontrigger) {
        histos.fill(HIST("REC_Trigger_Purity"), 0.5);
        histos.fill(HIST("REC_Trigger_Energy_GOOD"), mccluster.energy());
        histos.fill(HIST("REC_Trigger_Purity_v_Energy"), 0.5, mccluster.energy());
      }
      if (goodgentrigger && !photontrigger) {
        histos.fill(HIST("REC_Trigger_Purity"), 1.5);
        histos.fill(HIST("REC_Trigger_Energy_MISS"), mccluster.energy());
        histos.fill(HIST("REC_Trigger_Purity_v_Energy"), 1.5, mccluster.energy());
      }
      if (!goodgentrigger && photontrigger) {
        histos.fill(HIST("REC_Trigger_Purity"), 2.5);
        histos.fill(HIST("REC_Trigger_Energy_FAKE"), mccluster.energy());
        histos.fill(HIST("REC_Trigger_Purity_v_Energy"), 2.5, mccluster.energy());
      }
    } // cluster loop

    // clusters done, now we do the sternheimer tracks
    for (auto& track : tracks) {
      bool sterntrigger = false;
      double sternPt = 0.0;
      if (cfgJETracks) {
        if (!jetderiveddatautilities::selectTrack(track, trackFilter)) {
          continue;
        }
      } else {
        auto ogtrack = track.track_as<TrackCandidates>();
        if (!trackSelection(ogtrack)) {
          continue;
        }
        if (!ogtrack.isGlobalTrack()) {
          continue;
        }
      }

      // Do stuff with geometric cuts
      double phiPrime = track.phi();
      if (track.sign() < 0) {
        phiPrime = 2 * TMath::Pi() - phiPrime;
      }

      phiPrime = phiPrime + TMath::Pi() / 18.;
      phiPrime = fmod(phiPrime, 2 * TMath::Pi() / 18.);
      histos.fill(HIST("REC_Track_PhiPrime_Pt"), phiPrime, track.pt());
      histos.fill(HIST("REC_Track_Pt"), track.pt());
      histos.fill(HIST("REC_Track_Phi"), track.phi());
      if (clustertrigger) {
        histos.fill(HIST("REC_TrackPhi_photontrigger"), track.phi());
        histos.fill(HIST("REC_TrackEta_photontrigger"), track.eta());
      }
      if (track.pt() > cfgMinTrig && track.pt() < cfgMaxTrig) {
        if (fabs(track.eta()) <= cfgtrkMaxEta) {
          sterntrigger = true;
          sternPt = track.pt();
        }
      }
      double pthadsum = GetPtHadSum(tracks, track, cfgMinR, cfgMaxR, true, false, true);
      histos.fill(HIST("REC_Trigger_V_PtHadSum_Nch"), sternPt, pthadsum, weight);
      if (sterntrigger) {
        bool doStern = true;
        double sterncount = 1.0;
        while (doStern) {
          histos.fill(HIST("REC_Trigger_V_PtHadSum_Stern"), sterncount, pthadsum, (2.0 / sternPt) * weight);
          if (sterncount < sternPt) {
            sterncount++;
          } else {
            doStern = false;
          }
        } // While sternin'
      } // stern trigger loop
    } // track loop

  } // end of process

  PROCESS_SWITCH(statPromptPhoton, processMCRec_JE, "processJE  MC data", false);

  int nEventsData = 0;
  void processData(jfilteredDataCollisions::iterator const& collision, jfilteredClusters const& clusters, jDataTrackCandidates const& tracks, soa::Join<aod::JTracks, aod::JTrackExtras, aod::JTrackPIs> const&, TrackCandidates const&, BcCandidates const&, jEMCtracks const& emctracks)
  {
    nEventsData++;
    if (cfgDebug) {
      if (nEventsData == 1) {
        std::cout << "Starting Data Processing: " << nEventsData << std::endl;
      }
      if ((nEventsData + 1) % 10000 == 0) {
        std::cout << "Processed Data Events: " << nEventsData << std::endl;
        std::cout << "Events Trigger Bit: " << collision.triggerSel() << std::endl;
        std::cout << "Trigger Mask Bit: " << triggerMaskBits[0] << std::endl;
        std::cout << "Trigger Mask Cfg Line: " << cfgTriggerMasks << std::endl;
      }
    }

    histos.fill(HIST("DATA_nEvents"), 0.5);

    // required cuts
    if (fabs(collision.posZ()) > cfgVtxCut)
      return;
    if (!collision.sel8())
      return;

    histos.fill(HIST("DATA_nEvents"), 1.5);
    if (cfgEmcTrigger) {
      if (!collision.isEmcalReadout())
        return;
    }

    histos.fill(HIST("DATA_nEvents"), 2.5);

    if (!jetderiveddatautilities::selectTrigger(collision, triggerMaskBits)) {
      return;
    }

    histos.fill(HIST("DATA_nEvents"), 3.5);

    bool noTrk = true;
    for (auto& track : tracks) {

      if (cfgJETracks) {
        if (!jetderiveddatautilities::selectTrack(track, trackFilter)) {
          continue;
        }
      } else {
        auto ogtrack = track.track_as<TrackCandidates>();
        if (!trackSelection(ogtrack)) {
          continue;
        }
        if (!ogtrack.isGlobalTrack()) {
          continue;
        }
      }
      noTrk = false;
      break;
    }
    if (noTrk)
      return;

    // now we do clusters
    bool clustertrigger = false;
    for (auto& cluster : clusters) {
      histos.fill(HIST("DATA_M02_BC"), cluster.m02());
      if (cluster.m02() < cfgLowM02)
        continue;
      if (cluster.m02() > cfgHighM02)
        continue;
      if (cluster.energy() < cfgLowClusterE)
        continue;
      if (cluster.energy() > cfgHighClusterE)
        continue;
      if (fabs(cluster.eta()) > cfgtrkMaxEta)
        continue;

      histos.fill(HIST("DATA_Cluster_QA"), 0.5);
      histos.fill(HIST("DATA_M02_AC"), cluster.m02());
      histos.fill(HIST("DATA_All_Energy"), cluster.energy());
      histos.fill(HIST("DATA_ClusterPhi"), cluster.phi());
      histos.fill(HIST("DATA_ClusterEta"), cluster.eta());

      bool photontrigger = false; // is a neutral cluster
      bool chargetrigger = false; // is definitely not a neutral cluster
      auto tracksofcluster = cluster.matchedTracks_as<soa::Join<aod::JTracks, aod::JTrackExtras, aod::JTrackPIs>>();

      ///*GEOMETRICAL CUT*///

      double phiPrimeC = cluster.phi();
      phiPrimeC = phiPrimeC + TMath::Pi() / 18.;
      phiPrimeC = fmod(phiPrimeC, 2 * TMath::Pi() / 18.);
      double ptC = cluster.energy();
      bool geocut = false;
      histos.fill(HIST("DATA_Cluster_PhiPrime_Pt"), phiPrimeC, cluster.energy());
      if (phiPrimeC > (0.12 / ptC + TMath::Pi() / 18. + 0.035) ||
          phiPrimeC < (0.1 / ptC / ptC + TMath::Pi() / 18. - 0.025)) {
        geocut = false;
      } else {
        geocut = true;
      }

      if (cfgGeoCut) {
        if (geocut) {
          histos.fill(HIST("DATA_Cluster_PhiPrime_Pt_C"), phiPrimeC, cluster.energy());
          continue;
        }
      }

      histos.fill(HIST("DATA_Cluster_PhiPrime_Pt_AC"), phiPrimeC, cluster.energy());

      ///*GEOMETRICAL CUT*///

      ///*CHECK FOR PHOTON CANDIDATE*///

      // first, we check if veto is required
      double sumptT = 0;
      bool clusterqa = false;
      for (auto& ctrack : tracksofcluster) {
        double etaT, phiT;
        if (cfgJETracks) {
          if (!jetderiveddatautilities::selectTrack(ctrack, trackFilter)) {
            continue;
          }
          auto emctracksPerTrack = emctracks.sliceBy(EMCTrackPerTrack, ctrack.globalIndex());
          auto emctrack = emctracksPerTrack.iteratorAt(0);
          etaT = emctrack.etaEmcal();
          phiT = emctrack.phiEmcal();
        } else {
          auto ogtrack = ctrack.track_as<TrackCandidates>();
          if (!trackSelection(ogtrack)) {
            continue;
          }
          if (!ogtrack.isGlobalTrack()) {
            continue;
          }
          etaT = ogtrack.trackEtaEmcal();
          phiT = ogtrack.trackPhiEmcal();
        }

        double etaC = cluster.eta();
        double phiC = cluster.phi();
        double ptT = ctrack.pt();
        bool etatrigger = false;
        bool phitrigger = false;
        double phidiff = TVector2::Phi_mpi_pi(cluster.phi() - ctrack.phi());
        double etadiff = cluster.eta() - ctrack.eta();

        if (cfgPtClusterCut) {
          if (fabs(etaT - etaC) < (0.010 + pow(ptT + 4.07, -2.5))) {
            etatrigger = true;
          }

          if (fabs(TVector2::Phi_mpi_pi(phiT - phiC)) < (0.015 + pow(ptT + 3.65, -2.0))) {
            phitrigger = true;
          }
        } else {
          if (fabs(etadiff) < 0.05) {
            etatrigger = true;
          }

          if (fabs(phidiff) < 0.05) {
            phitrigger = true;
          }
        }

        if (etatrigger && phitrigger) {
          chargetrigger = true;
          sumptT += ptT;
        }
        if (chargetrigger) {
          if (!clusterqa) {
            histos.fill(HIST("DATA_Cluster_QA"), 1.5);
            clusterqa = true;
          }
        }
        histos.fill(HIST("DATA_Track_v_Cluster_Phi"), phidiff);
        histos.fill(HIST("DATA_Track_v_Cluster_Eta"), etadiff);
        histos.fill(HIST("DATA_Track_v_Cluster_Phi_Eta"), phidiff, etadiff);

      } // track of cluster loop

      if (chargetrigger && sumptT > 0) {
        double cluster_over_sumptT = cluster.energy() / sumptT;
        histos.fill(HIST("DATA_SumPt_BC"), cluster_over_sumptT);
        if (cluster_over_sumptT < 1.7) {
          histos.fill(HIST("DATA_Cluster_QA"), 2.5); // veto fails, cluster is charged
        } else {
          histos.fill(HIST("DATA_Cluster_QA"), 3.5); // veto is good, cluster is converted to neutral cluster
          chargetrigger = false;
          histos.fill(HIST("DATA_SumPt_AC"), cluster_over_sumptT);
        }
      } // sumptT check

      if (!chargetrigger) {
        photontrigger = true;
      }

      ///*CHECK FOR PHOTON CANDIDATE*///

      ///*CALCULATE PTHAD SUM*///

      if (photontrigger) { // if no charge trigger, cluster is good!
        histos.fill(HIST("DATA_Cluster_QA"), 4.5);
        clustertrigger = true;
        double pthadsum = GetPtHadSum(tracks, cluster, cfgMinR, cfgMaxR, false, false, true);
        histos.fill(HIST("DATA_Trigger_V_PtHadSum_Photon"), cluster.energy(), pthadsum);
        histos.fill(HIST("DATA_PtHadSum_Photon"), pthadsum);
        histos.fill(HIST("DATA_Trigger_Energy"), cluster.energy());
      }

      ///*CALCULATE PTHAD SUM*///

    } // cluster loop

    ///*CALCULATE STERNHEIMER*///

    for (auto& track : tracks) {
      bool sterntrigger = false;
      double sternPt = 0.0;
      if (cfgJETracks) {
        if (!jetderiveddatautilities::selectTrack(track, trackFilter)) {
          continue;
        }
      } else {
        auto ogtrack = track.track_as<TrackCandidates>();
        if (!trackSelection(ogtrack)) {
          continue;
        }
        if (!ogtrack.isGlobalTrack()) {
          continue;
        }
      }

      // Do stuff with geometric cuts
      double phiPrime = track.phi();
      if (track.sign() < 0) {
        phiPrime = 2 * TMath::Pi() - phiPrime;
      }

      phiPrime = phiPrime + TMath::Pi() / 18.;
      phiPrime = fmod(phiPrime, 2 * TMath::Pi() / 18.);
      double pt = track.pt();
      if (phiPrime > (0.12 / pt + TMath::Pi() / 18. + 0.035) ||
          phiPrime < (0.1 / pt / pt + TMath::Pi() / 18. - 0.025)) {
        histos.fill(HIST("DATA_Track_PhiPrime_Pt"), phiPrime, track.pt());
      } // geo cut
      // Done with geometric cuts

      histos.fill(HIST("DATA_Track_Pt"), track.pt());
      histos.fill(HIST("DATA_Track_Phi"), track.phi());
      if (clustertrigger) {
        histos.fill(HIST("DATA_TrackPhi_photontrigger"), track.phi());
        histos.fill(HIST("DATA_TrackEta_photontrigger"), track.eta());
      }
      if (track.pt() > cfgMinTrig && track.pt() < cfgMaxTrig) {
        if (fabs(track.eta()) <= cfgtrkMaxEta) {
          sterntrigger = true;
          sternPt = track.pt();
        }
      }

      if (sterntrigger) {
        bool doStern = true;
        double sterncount = 1.0;
        double pthadsum = GetPtHadSum(tracks, track, cfgMinR, cfgMaxR, true, false, true);
        histos.fill(HIST("DATA_Trigger_V_PtHadSum_Nch"), sternPt, pthadsum);
        while (doStern) {
          double pthadsum = GetPtHadSum(tracks, track, cfgMinR, cfgMaxR, true, false, true);
          histos.fill(HIST("DATA_Trigger_V_PtHadSum_Stern"), sterncount, pthadsum, 2.0 / sternPt);
          if (sterncount < sternPt) {
            sterncount++;
          } else {
            doStern = false;
          }
        } // While sternin'
      } // stern trigger loop
    } // track loop

    ///*CALCULATE STERNHEIMER*///

  } // end of process

  PROCESS_SWITCH(statPromptPhoton, processData, "processJE data", false);

}; // end of main struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<statPromptPhoton>(cfgc)};
};
