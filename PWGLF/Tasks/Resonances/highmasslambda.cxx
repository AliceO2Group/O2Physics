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
// Phi meson spin alignment task
// sourav.kundu@cern.ch

#include <TH1F.h>
#include <TDirectory.h>
#include <THn.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TFile.h>
#include <TH2F.h>
#include <TLorentzVector.h>
#include <TPDGCode.h>
#include <TDatabasePDG.h>
#include <cmath>
#include <array>
#include <cstdlib>

#include "TRandom3.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "Math/GenVector/Boost.h"
#include "TF1.h"

#include "PWGLF/DataModel/EPCalibrationTables.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/StepTHn.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/Core/trackUtilities.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Common/Core/TrackSelection.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;
struct highmasslambda {

  int mRunNumber;
  int multEstimator;
  float d_bz;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Service<o2::framework::O2DatabasePDG> pdg;

  // CCDB options
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};

  // fill output
  Configurable<bool> additionalEvSel{"additionalEvSel", true, "additionalEvSel"};
  Configurable<bool> fillDefault{"fillDefault", true, "fill Occupancy"};
  Configurable<bool> fillOccupancy{"fillOccupancy", true, "fill Occupancy"};
  Configurable<bool> fillDecayLength{"fillDecayLength", false, "fill decay length"};
  Configurable<bool> fillPolarization{"fillPolarization", false, "fill polarization"};
  Configurable<bool> fillRotation{"fillRotation", false, "fill rotation"};
  // events
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCutCentralityMax{"cfgCutCentralityMax", 50.0f, "Accepted maximum Centrality"};
  Configurable<float> cfgCutCentralityMin{"cfgCutCentralityMin", 30.0f, "Accepted minimum Centrality"};
  // proton track cut
  Configurable<bool> isPVContributor{"isPVContributor", true, "is PV contributor"};
  Configurable<bool> rejectPID{"rejectPID", true, "Reject PID"};
  Configurable<float> confRapidity{"confRapidity", 0.8, "cut on Rapidity"};
  Configurable<float> cfgCutPT{"cfgCutPT", 0.3, "PT cut on daughter track"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8, "Eta cut on daughter track"};
  Configurable<float> cfgCutDCAxymin1{"cfgCutDCAxymin1", 0.005f, "Minimum DCAxy range for tracks pt 0 to 0.5"};
  Configurable<float> cfgCutDCAxymin2{"cfgCutDCAxymin2", 0.003f, "Minimum DCAxy range for tracks pt 0.5 to 1"};
  Configurable<float> cfgCutDCAxymin3{"cfgCutDCAxymin3", 0.003f, "Minimum DCAxy range for tracks pt 1.0 to 1.5"};
  Configurable<float> cfgCutDCAxymin4{"cfgCutDCAxymin4", 0.002f, "Minimum DCAxy range for tracks pt 1.5 to 2.0"};
  Configurable<float> cfgCutDCAxymin5{"cfgCutDCAxymin5", 0.001f, "Minimum DCAxy range for tracks pt 2.0 to 2.5"};
  Configurable<float> cfgCutDCAxymin6{"cfgCutDCAxymin6", 0.0003f, "Minimum DCAxy range for tracks pt 2.5 to 3.0"};
  Configurable<float> cfgCutDCAxymin7{"cfgCutDCAxymin7", 0.0003f, "Minimum DCAxy range for tracks pt 3.0 to 4.0"};
  Configurable<float> cfgCutDCAxymin8{"cfgCutDCAxymin8", 0.0003f, "Minimum DCAxy range for tracks pt 4.0 to 10.0"};
  Configurable<float> cfgCutDCAxy{"cfgCutDCAxy", 0.1f, "DCAxy range for tracks"};
  Configurable<float> cfgCutDCAz{"cfgCutDCAz", 1.0f, "DCAz range for tracks"};
  Configurable<int> cfgITScluster{"cfgITScluster", 5, "Number of ITS cluster"};
  Configurable<int> cfgTPCcluster{"cfgTPCcluster", 70, "Number of TPC cluster"};
  Configurable<bool> ispTdepPID{"ispTdepPID", true, "pT dependent PID"};
  Configurable<float> nsigmaCutTPC{"nsigmacutTPC", 3.0, "Value of the TPC Nsigma cut"};
  Configurable<float> nsigmaCutTOF{"nsigmacutTOF", 3.0, "Value of the TOF Nsigma cut"};
  Configurable<float> nsigmaCutTPCPre{"nsigmacutTPCPre", 5.0, "Value of the TPC Nsigma cut Pre filter"};
  // Configs for V0
  Configurable<float> ConfV0PtMin{"ConfV0PtMin", 0.f, "Minimum transverse momentum of V0"};
  Configurable<double> ConfV0DCADaughMax{"ConfV0DCADaughMax", 0.2f, "Maximum DCA between the V0 daughters"};
  Configurable<double> ConfV0CPAMin{"ConfV0CPAMin", 0.9998f, "Minimum CPA of V0"};
  Configurable<float> ConfV0TranRadV0Min{"ConfV0TranRadV0Min", 1.5f, "Minimum transverse radius"};
  Configurable<float> ConfV0TranRadV0Max{"ConfV0TranRadV0Max", 100.f, "Maximum transverse radius"};
  Configurable<double> cMaxV0DCA{"cMaxV0DCA", 0.2, "Maximum V0 DCA to PV"};
  Configurable<float> cMaxV0LifeTime{"cMaxV0LifeTime", 20, "Maximum V0 life time"};
  Configurable<float> cSigmaMassKs0{"cSigmaMassKs0", 0.006, "Sigma cut on KS0 mass"};
  // config for V0 daughters
  Configurable<float> ConfDaughEta{"ConfDaughEta", 0.8f, "V0 Daugh sel: max eta"};
  Configurable<float> ConfDaughPt{"ConfDaughPt", 0.1f, "V0 Daugh sel: min pt"};
  Configurable<float> ConfDaughTPCnclsMin{"ConfDaughTPCnclsMin", 50.f, "V0 Daugh sel: Min. nCls TPC"};
  Configurable<double> ConfDaughDCAMin{"ConfDaughDCAMin", 0.08f, "V0 Daugh sel:  Max. DCA Daugh to PV (cm)"};
  Configurable<float> ConfDaughPIDCuts{"ConfDaughPIDCuts", 3, "PID selections for KS0 daughters"};
  // Fill strategy
  Configurable<int> cfgSelectDaughterTopology{"cfgSelectDaughterTopology", 2, "Select daughter for topology"};
  // Mixed event
  Configurable<int> cfgNoMixedEvents{"cfgNoMixedEvents", 1, "Number of mixed events per event"};
  /// activate rotational background
  Configurable<int> nBkgRotations{"nBkgRotations", 9, "Number of rotated copies (background) per each original candidate"};

  // THnsparse bining
  ConfigurableAxis configThnAxisInvMass{"configThnAxisInvMass", {60, 2.15, 2.45}, "#it{M} (GeV/#it{c}^{2})"};
  ConfigurableAxis configThnAxisPt{"configThnAxisPt", {24, 1.0, 25.}, "#it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis configThnAxisCosThetaStar{"configThnAxisCosThetaStar", {10, -1.0, 1.}, "cos(#vartheta)"};
  ConfigurableAxis configThnAxisCentrality{"configThnAxisCentrality", {1, 30., 50}, "Centrality"};
  ConfigurableAxis configThnAxisPhiminusPsi{"configThnAxisPhiminusPsi", {6, 0.0, TMath::Pi()}, "#phi - #psi"};
  ConfigurableAxis configThnAxisV2{"configThnAxisV2", {80, -1, 1}, "V2"};
  ConfigurableAxis configThnAxisSA{"configThnAxisSA", {100, -1, 1}, "SA"};
  ConfigurableAxis cnfigThnAxisDecayLength{"cnfigThnAxisDecayLength", {150, 0.0, 0.3}, "SA"};
  ConfigurableAxis cnfigThnAxisPtProton{"cnfigThnAxisPtProton", {16, 0.0, 8.0}, "pT"};
  ConfigurableAxis configThnAxisSP{"configThnAxisSP", {400, -12, 12}, "SP"};

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter centralityFilter = (nabs(aod::cent::centFT0C) < cfgCutCentralityMax && nabs(aod::cent::centFT0C) > cfgCutCentralityMin);
  Filter acceptanceFilter = (nabs(aod::track::eta) < cfgCutEta && nabs(aod::track::pt) > cfgCutPT);
  Filter dcaCutFilter = (nabs(aod::track::dcaXY) < cfgCutDCAxy) && (nabs(aod::track::dcaZ) < cfgCutDCAz);
  Filter pidFilter = nabs(aod::pidtpc::tpcNSigmaPr) < nsigmaCutTPCPre;

  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::FV0Mults, aod::TPCMults, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0As, aod::EPCalibrationTables, aod::Mults>>;
  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPr, aod::pidTOFFullPr, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa>>;
  using AllTrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPi>;
  using ResoV0s = aod::V0Datas;

  SliceCache cache;
  // Partition<TrackCandidates> posTracks = aod::track::signed1Pt > cfgCutCharge;
  // Partition<TrackCandidates> negTracks = aod::track::signed1Pt < cfgCutCharge;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext&)
  {
    std::vector<double> occupancyBinning = {0.0, 500.0, 1000.0, 3000.0, 6000.0, 50000.0};
    // std::vector<double> dcaBinning = {0.0, 0.0005, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01, 0.012, 0.014, 0.016, 0.02, 0.03, 0.05, 0.1, 0.5, 1.0};
    std::vector<double> dcaBinning = {0.0, 0.0005, 0.001, 0.0015, 0.002, 0.003, 0.004, 0.006, 0.008, 0.01, 0.015, 0.02, 0.04, 1.0};
    std::vector<double> ptProtonBinning = {0.2, 0.3, 0.5, 1.0, 1.5, 100.0};

    const AxisSpec thnAxisInvMass{configThnAxisInvMass, "#it{M} (GeV/#it{c}^{2})"};
    const AxisSpec thnAxisPt{configThnAxisPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec thnAxisCosThetaStar{configThnAxisCosThetaStar, "cos(#vartheta)"};
    const AxisSpec thnAxisPhiminusPsi{configThnAxisPhiminusPsi, "#phi - #psi"};
    const AxisSpec thnAxisCentrality{configThnAxisCentrality, "Centrality (%)"};
    const AxisSpec thnAxisV2{configThnAxisV2, "V2"};
    const AxisSpec thnAxisSA{configThnAxisSA, "SA"};
    const AxisSpec thnAxisDecayLength{cnfigThnAxisDecayLength, "Decay Length"};
    const AxisSpec thnAxisPtProton{cnfigThnAxisPtProton, "Proton Pt"};
    const AxisSpec thnAxisSP{configThnAxisSP, "SP"};

    AxisSpec phiAxis = {500, -6.28, 6.28, "phi"};
    AxisSpec resAxis = {1000, -10, 10, "Res"};
    AxisSpec centAxis = {8, 0, 80, "V0M (%)"};
    // AxisSpec ptProtonAxis = {16, 0.0, 8, "V0M (%)"};
    AxisSpec dcaAxis = {dcaBinning, "DCAxy"};
    AxisSpec dcatoPVAxis = {50, 0.0, 0.4, "V0M (%)"};
    AxisSpec occupancyAxis = {occupancyBinning, "occupancy"};
    AxisSpec ptProtonAxis = {ptProtonBinning, "pT proton"};

    histos.add("hRejectPID", "hRejectPID", kTH1F, {{2, 0.0f, 2.0f}});
    histos.add("hMomCorr", "hMomCorr", kTH3F, {{200, -10.0f, 10.0f}, {200, -10.0f, 10.0f}, {8, 0.0f, 80.0f}});
    histos.add("hInvMassKs0", "hInvMassKs0", kTH1F, {{200, 0.4f, 0.6f}});
    histos.add("hV0Dca", "hV0Dca", kTH1F, {{2000, -1.0f, 1.0f}});
    histos.add("hpTvsRapidity", "pT vs Rapidity", kTH2F, {{100, 0.0f, 10.0f}, {300, -1.5f, 1.5f}});
    histos.add("hFTOCvsTPCNoCut", "Mult correlation FT0C vs. TPC without any cut", kTH2F, {{80, 0.0f, 80.0f}, {100, -0.5f, 5999.5f}});
    histos.add("hFTOCvsTPC", "Mult correlation FT0C vs. TPC", kTH2F, {{80, 0.0f, 80.0f}, {100, -0.5f, 5999.5f}});
    histos.add("hFTOCvsTPCSelected", "Mult correlation FT0C vs. TPC after selection", kTH2F, {{80, 0.0f, 80.0f}, {100, -0.5f, 5999.5f}});
    histos.add("hCentrality", "Centrality distribution", kTH1F, {{200, 0.0, 200.0}});
    histos.add("hVtxZ", "Vertex distribution in Z;Z (cm)", kTH1F, {{400, -20.0, 20.0}});
    histos.add("hOccupancy", "Occupancy", kTH1F, {{5000, 0.0, 50000.0}});
    histos.add("hEta", "Eta distribution", kTH1F, {{200, -1.0f, 1.0f}});
    histos.add("hDcaxy", "Dcaxy distribution", kTH1F, {{1000, -0.5f, 0.5f}});
    histos.add("hDcaz", "Dcaz distribution", kTH1F, {{1000, -0.5f, 0.5f}});
    histos.add("hNsigmaProtonTPCDiff", "Difference NsigmaProton NsigmaKaon TPC distribution", kTH3F, {{100, -5.0f, 5.0f}, {100, -5.0f, 5.0f}, {80, 0.0f, 8.0f}});
    histos.add("hNsigmaProtonTPC", "NsigmaProton TPC distribution", kTH2F, {{100, -5.0f, 5.0f}, {80, 0.0f, 8.0f}});
    histos.add("hNsigmaProtonTOF", "NsigmaProton TOF distribution", kTH2F, {{100, -5.0f, 5.0f}, {80, 0.0f, 8.0f}});
    histos.add("hPsiFT0C", "PsiFT0C", kTH2F, {centAxis, phiAxis});
    histos.add("hPsiFT0A", "PsiFT0A", kTH2F, {centAxis, phiAxis});
    histos.add("hPsiTPC", "PsiTPC", kTH2F, {centAxis, phiAxis});
    histos.add("hPsiTPCR", "PsiTPCR", kTH2F, {centAxis, phiAxis});
    histos.add("hPsiTPCL", "PsiTPCL", kTH2F, {centAxis, phiAxis});
    if (fillDefault) {
      if (fillDecayLength) {
        histos.add("hSparseV2SASameEvent_V2", "hSparseV2SASameEvent_V2", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisSP, dcaAxis, ptProtonAxis});
        histos.add("hSparseV2SAMixedEvent_V2", "hSparseV2SAMixedEvent_V2", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisSP, dcaAxis, ptProtonAxis});
        histos.add("hSparseV2SASameEventRotational_V2", "hSparseV2SASameEventRotational_V2", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisSP, dcaAxis, ptProtonAxis});
      }
      histos.add("hSparseV2SASameEvent_V2_new", "hSparseV2SASameEvent_V2_new", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisSP, dcaAxis, ptProtonAxis});
      histos.add("hSparseV2SAMixedEvent_V2_new", "hSparseV2SAMixedEvent_V2_new", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisSP, dcaAxis, ptProtonAxis});
      histos.add("hSparseV2SASameEventRotational_V2_new", "hSparseV2SASameEventRotational_V2_new", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisSP, dcaAxis, ptProtonAxis});
    }
    if (fillOccupancy) {
      if (fillDecayLength) {
        histos.add("hSparseV2SASameEvent_V2_occupancy", "hSparseV2SASameEvent_V2_occupancy", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisSP, dcaAxis, ptProtonAxis, occupancyAxis});
        histos.add("hSparseV2SASameEventRotational_V2_occupancy", "hSparseV2SASameEventRotational_V2_occupancy", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisSP, dcaAxis, ptProtonAxis, occupancyAxis});
        histos.add("hSparseV2SAMixedEvent_V2_occupancy", "hSparseV2SAMixedEvent_V2", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisSP, dcaAxis, ptProtonAxis, occupancyAxis});
      }
      histos.add("hSparseV2SASameEvent_V2_new_occupancy", "hSparseV2SASameEvent_V2_new_occupancy", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisSP, dcaAxis, ptProtonAxis, occupancyAxis});
      histos.add("hSparseV2SASameEventRotational_V2_new_occupancy", "hSparseV2SASameEventRotational_V2_new_occupancy", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisSP, dcaAxis, ptProtonAxis, occupancyAxis});
      histos.add("hSparseV2SAMixedEvent_V2_new_occupancy", "hSparseV2SAMixedEvent_V2_new", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisSP, dcaAxis, ptProtonAxis, occupancyAxis});
    }
    if (fillPolarization) {
      histos.add("hSparseV2SASameEventplus_SA", "hSparseV2SASameEventplus_SA", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisSA, thnAxisPhiminusPsi, thnAxisCentrality});
      histos.add("hSparseV2SASameEventplus_SA_A0", "hSparseV2SASameEventplus_SA_A0", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisSA, thnAxisPhiminusPsi, thnAxisCentrality});
      histos.add("hSparseV2SASameEventplus_SA_azimuth", "hSparseV2SASameEventplus_SA_azimuth", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisSA, thnAxisCentrality});
      histos.add("hSparseV2SASameEventminus_SA", "hSparseV2SASameEventminus_SA", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisSA, thnAxisPhiminusPsi, thnAxisCentrality});
      histos.add("hSparseV2SASameEventminus_SA_A0", "hSparseV2SASameEventminus_SA_A0", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisSA, thnAxisPhiminusPsi, thnAxisCentrality});
      histos.add("hSparseV2SASameEventminus_SA_azimuth", "hSparseV2SASameEventminus_SA_azimuth", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisSA, thnAxisCentrality});

      histos.add("hSparseV2SAMixedEventplus_SA", "hSparseV2SAMixedEventplus_SA", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisSA, thnAxisPhiminusPsi, thnAxisCentrality});
      histos.add("hSparseV2SAMixedEventplus_SA_A0", "hSparseV2SAMixedEventplus_SA_A0", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisSA, thnAxisPhiminusPsi, thnAxisCentrality});
      histos.add("hSparseV2SAMixedEventplus_SA_azimuth", "hSparseV2SAMixedEventplus_SA_azimuth", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisSA, thnAxisCentrality});
      histos.add("hSparseV2SAMixedEventminus_SA", "hSparseV2SAMixedEventminus_SA", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisSA, thnAxisPhiminusPsi, thnAxisCentrality});
      histos.add("hSparseV2SAMixedEventminus_SA_A0", "hSparseV2SAMixedEventminus_SA_A0", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisSA, thnAxisPhiminusPsi, thnAxisCentrality});
      histos.add("hSparseV2SAMixedEventminus_SA_azimuth", "hSparseV2SAMixedEventminus_SA_azimuth", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisSA, thnAxisCentrality});
    }

    // histogram for resolution
    histos.add("ResFT0CTPC", "ResFT0CTPC", kTH2F, {centAxis, resAxis});
    histos.add("ResFT0CTPCR", "ResFT0CTPCR", kTH2F, {centAxis, resAxis});
    histos.add("ResFT0CTPCL", "ResFT0CTPCL", kTH2F, {centAxis, resAxis});
    histos.add("ResTPCRTPCL", "ResTPCRTPCL", kTH2F, {centAxis, resAxis});
    histos.add("ResFT0CFT0A", "ResFT0CFT0A", kTH2F, {centAxis, resAxis});
    histos.add("ResFT0ATPC", "ResFT0ATPC", kTH2F, {centAxis, resAxis});
  }
  template <typename T>
  bool selectionTrack(const T& candidate)
  {
    if (isPVContributor && !(candidate.isGlobalTrack() && candidate.isPVContributor() && candidate.itsNCls() > cfgITScluster && candidate.tpcNClsFound() > cfgTPCcluster)) {
      return false;
    }
    if (!isPVContributor && !(candidate.isGlobalTrackWoDCA() && candidate.itsNCls() > cfgITScluster && candidate.tpcNClsFound() > cfgTPCcluster)) {
      return false;
    }
    if (candidate.pt() > 0.0 && candidate.pt() < 0.5 && TMath::Abs(candidate.dcaXY()) < cfgCutDCAxymin1) {
      return false;
    }
    if (candidate.pt() >= 0.5 && candidate.pt() < 1.0 && TMath::Abs(candidate.dcaXY()) < cfgCutDCAxymin2) {
      return false;
    }
    if (candidate.pt() >= 1.0 && candidate.pt() < 1.5 && TMath::Abs(candidate.dcaXY()) < cfgCutDCAxymin3) {
      return false;
    }
    if (candidate.pt() >= 1.5 && candidate.pt() < 2.0 && TMath::Abs(candidate.dcaXY()) < cfgCutDCAxymin4) {
      return false;
    }
    if (candidate.pt() >= 2.0 && candidate.pt() < 2.5 && TMath::Abs(candidate.dcaXY()) < cfgCutDCAxymin5) {
      return false;
    }
    if (candidate.pt() >= 2.5 && candidate.pt() < 3.0 && TMath::Abs(candidate.dcaXY()) < cfgCutDCAxymin6) {
      return false;
    }
    if (candidate.pt() >= 3.0 && candidate.pt() < 4.0 && TMath::Abs(candidate.dcaXY()) < cfgCutDCAxymin7) {
      return false;
    }
    if (candidate.pt() >= 4.0 && candidate.pt() < 10.0 && TMath::Abs(candidate.dcaXY()) < cfgCutDCAxymin8) {
      return false;
    }
    return true;
  }

  template <typename T>
  bool selectionPIDpTdependent(const T& candidate)
  {
    if (candidate.p() <= 0.5 && TMath::Abs(candidate.tpcNSigmaPr()) < 5.0) {
      return true;
    }
    if (candidate.p() > 0.5 && candidate.p() <= 0.8 && TMath::Abs(candidate.tpcNSigmaPr()) < 3.0) {
      return true;
    }
    if (candidate.p() > 0.8 && candidate.hasTOF() && TMath::Sqrt(candidate.tofNSigmaPr() * candidate.tofNSigmaPr() + candidate.tpcNSigmaPr() * candidate.tpcNSigmaPr()) < nsigmaCutTOF) {
      return true;
    }
    return false;
  }

  template <typename T>
  bool selectionPID(const T& candidate)
  {
    if (candidate.hasTOF()) {
      if (candidate.pt() < 0.8 && candidate.hasTOF() && TMath::Abs(candidate.tpcNSigmaPr()) < 3.0 && TMath::Abs(candidate.tofNSigmaPr()) < 10.0) {
        return true;
      }
      if (candidate.pt() >= 0.8 && candidate.pt() < 3.0 && candidate.hasTOF() && TMath::Abs(candidate.tpcNSigmaPr()) < 3.0 && TMath::Abs(candidate.tofNSigmaPr()) < 3.0) {
        return true;
      }
      if (candidate.pt() >= 3.0 && candidate.pt() < 4.0 && candidate.hasTOF() && TMath::Abs(candidate.tpcNSigmaPr()) < 3.0 && candidate.tofNSigmaPr() > -2.0 && candidate.tofNSigmaPr() < 3.0) {
        return true;
      }
      if (candidate.pt() >= 4.0 && candidate.pt() < 5.0 && candidate.hasTOF() && TMath::Abs(candidate.tpcNSigmaPr()) < 3.0 && candidate.tofNSigmaPr() > -1.5 && candidate.tofNSigmaPr() < 3.0) {
        return true;
      }
      if (candidate.pt() >= 5.0 && candidate.hasTOF() && TMath::Abs(candidate.tpcNSigmaPr()) < 3.0 && candidate.tofNSigmaPr() > -1.0 && candidate.tofNSigmaPr() < 3.0) {
        return true;
      }
    }
    if (!candidate.hasTOF()) {
      if (candidate.pt() < 0.7 && TMath::Abs(candidate.tpcNSigmaPr()) < 3.0) {
        return true;
      }
      if (candidate.pt() >= 0.7 && candidate.pt() < 0.8 && candidate.tpcNSigmaPr() > -2.0 && candidate.tpcNSigmaPr() < 3.0) {
        return true;
      }
      if (candidate.pt() >= 0.8 && candidate.tpcNSigmaPr() > -1.0 && candidate.tpcNSigmaPr() < 3.0) {
        return true;
      }
    }
    return false;
  }

  template <typename T>
  bool selectionPIDNew(const T& candidate)
  {
    if (candidate.pt() < 0.7 && !candidate.hasTOF() && TMath::Abs(candidate.tpcNSigmaPr()) < 3.0) {
      return true;
    }
    if (candidate.pt() >= 0.7 && !candidate.hasTOF() && candidate.pt() < 0.8 && candidate.tpcNSigmaPr() > -2.0 && candidate.tpcNSigmaPr() < 3.0) {
      return true;
    }
    if (candidate.pt() < 0.8 && candidate.hasTOF() && TMath::Abs(candidate.tpcNSigmaPr()) < 3.0 && TMath::Abs(candidate.tofNSigmaPr()) < 10.0) {
      return true;
    }

    if (candidate.pt() >= 0.8 && candidate.pt() < 1.0 && !candidate.hasTOF() && candidate.tpcNSigmaPr() > -1.0 && candidate.tpcNSigmaPr() < 3.0) {
      return true;
    }
    if (candidate.pt() >= 0.8 && candidate.pt() < 3.0 && candidate.hasTOF() && TMath::Abs(candidate.tpcNSigmaPr()) < 3.0 && TMath::Abs(candidate.tofNSigmaPr()) < 3.0) {
      return true;
    }

    if (candidate.pt() >= 3.0 && candidate.pt() < 4.0 && candidate.hasTOF() && TMath::Abs(candidate.tpcNSigmaPr()) < 3.0 && candidate.tofNSigmaPr() > -2.0 && candidate.tofNSigmaPr() < 3.0) {
      return true;
    }
    if (candidate.pt() >= 4.0 && candidate.pt() < 5.0 && candidate.hasTOF() && TMath::Abs(candidate.tpcNSigmaPr()) < 3.0 && candidate.tofNSigmaPr() > -1.5 && candidate.tofNSigmaPr() < 3.0) {
      return true;
    }

    if (candidate.p() > 5.0 && !candidate.hasTOF() && candidate.tpcNSigmaPr() > -1.0 && candidate.tpcNSigmaPr() < 3.0) {
      return true;
    }
    if (candidate.pt() >= 5.0 && candidate.hasTOF() && TMath::Abs(candidate.tpcNSigmaPr()) < 3.0 && candidate.tofNSigmaPr() > -1.0 && candidate.tofNSigmaPr() < 3.0) {
      return true;
    }
    return false;
  }

  template <typename T>
  bool RejectPion(const T& candidate)
  {
    if (candidate.p() > 1.0 && candidate.p() < 2.0 && !candidate.hasTOF() && candidate.tpcNSigmaPi() < 2) {
      return false;
    }
    return true;
  }

  template <typename T>
  bool RejectKaon(const T& candidate)
  {
    if (candidate.p() > 1.0 && candidate.p() < 2.0 && !candidate.hasTOF() && candidate.tpcNSigmaKa() < 2) {
      return false;
    }
    return true;
  }

  template <typename Collision, typename V0>
  bool SelectionV0(Collision const& collision, V0 const& candidate)
  {
    if (fabs(candidate.dcav0topv()) > cMaxV0DCA) {
      return false;
    }
    const float pT = candidate.pt();
    const std::vector<float> decVtx = {candidate.x(), candidate.y(), candidate.z()};
    const float tranRad = candidate.v0radius();
    const double dcaDaughv0 = TMath::Abs(candidate.dcaV0daughters());
    const double cpav0 = candidate.v0cosPA();

    float CtauK0s = candidate.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * TDatabasePDG::Instance()->GetParticle(kK0Short)->Mass(); // FIXME: Get from the common header
    float lowmasscutks0 = 0.497 - 2.0 * cSigmaMassKs0;
    float highmasscutks0 = 0.497 + 2.0 * cSigmaMassKs0;

    if (pT < ConfV0PtMin) {
      return false;
    }
    if (dcaDaughv0 > ConfV0DCADaughMax) {
      return false;
    }
    if (cpav0 < ConfV0CPAMin) {
      return false;
    }
    if (tranRad < ConfV0TranRadV0Min) {
      return false;
    }
    if (tranRad > ConfV0TranRadV0Max) {
      return false;
    }
    if (fabs(CtauK0s) > cMaxV0LifeTime || candidate.mK0Short() < lowmasscutks0 || candidate.mK0Short() > highmasscutks0) {
      return false;
    }
    return true;
  }
  template <typename T>
  bool isSelectedV0Daughter(T const& track, float charge)
  {
    const auto eta = track.eta();
    const auto pt = track.pt();
    const auto tpcNClsF = track.tpcNClsFound();
    const auto dcaXY = track.dcaXY();
    const auto sign = track.sign();
    if (charge < 0 && sign > 0) {
      return false;
    }
    if (charge > 0 && sign < 0) {
      return false;
    }
    if (std::abs(eta) > ConfDaughEta) {
      return false;
    }
    if (std::abs(pt) < ConfDaughPt) {
      return false;
    }
    if (tpcNClsF < ConfDaughTPCnclsMin) {
      return false;
    }
    if (std::abs(dcaXY) < ConfDaughDCAMin) {
      return false;
    }
    if (std::abs(track.tpcNSigmaPi()) > ConfDaughPIDCuts) {
      return false;
    }
    return true;
  }

  double GetPhiInRange(double phi)
  {
    double result = phi;
    while (result < 0) {
      result = result + 2. * TMath::Pi() / 2;
    }
    while (result > 2. * TMath::Pi() / 2) {
      result = result - 2. * TMath::Pi() / 2;
    }
    return result;
  }

  ConfigurableAxis axisVertex{"axisVertex", {10, -10, 10}, "vertex axis for bin"};
  ConfigurableAxis axisMultiplicityClass{"axisMultiplicityClass", {10, 0, 100}, "multiplicity percentile for bin"};
  ConfigurableAxis axisEPAngle{"axisEPAngle", {1, -TMath::Pi() / 2, TMath::Pi() / 2}, "event plane angle"};

  using BinningTypeVertexContributor = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C, aod::epcalibrationtable::PsiFT0C>;
  ROOT::Math::PxPyPzMVector Lambdac, Proton, Kshort, LambdacRot, ProtonRot, fourVecDauCM;
  ROOT::Math::XYZVector threeVecDauCM, threeVecDauCMXY, eventplaneVec, eventplaneVecNorm, beamvector;
  double massPi = TDatabasePDG::Instance()->GetParticle(kPiPlus)->Mass();   // FIXME: Get from the common header
  double massPr = TDatabasePDG::Instance()->GetParticle(kProton)->Mass();   // FIXME: Get from the common header
  double massK0s = TDatabasePDG::Instance()->GetParticle(kK0Short)->Mass(); // FIXME: Get from the common header
  void processSameEvent(EventCandidates::iterator const& collision, TrackCandidates const& tracks, AllTrackCandidates const&, ResoV0s const& V0s, aod::BCs const&)
  {
    if (!collision.sel8()) {
      return;
    }
    auto centrality = collision.centFT0C();
    auto multTPC = collision.multNTracksPV();
    histos.fill(HIST("hFTOCvsTPCNoCut"), centrality, multTPC);
    if (!collision.triggereventep()) {
      return;
    }
    if (additionalEvSel && (!collision.selection_bit(aod::evsel::kNoSameBunchPileup) || !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV))) {
      return;
    }
    int occupancy = collision.trackOccupancyInTimeRange();
    auto psiFT0C = collision.psiFT0C();
    auto psiFT0A = collision.psiFT0A();
    auto psiTPC = collision.psiTPC();
    auto psiTPCR = collision.psiTPCR();
    auto psiTPCL = collision.psiTPCL();

    auto QFT0C = collision.qFT0C();
    auto QFT0A = collision.qFT0A();
    auto QTPC = collision.qTPC();
    auto QTPCR = collision.qTPCR();
    auto QTPCL = collision.qTPCL();

    histos.fill(HIST("hFTOCvsTPC"), centrality, multTPC);
    histos.fill(HIST("hFTOCvsTPCSelected"), centrality, multTPC);
    histos.fill(HIST("hPsiFT0C"), centrality, psiFT0C);
    histos.fill(HIST("hPsiFT0A"), centrality, psiFT0A);
    histos.fill(HIST("hPsiTPC"), centrality, psiTPC);
    histos.fill(HIST("hPsiTPCR"), centrality, psiTPCR);
    histos.fill(HIST("hPsiTPCL"), centrality, psiTPCL);
    histos.fill(HIST("ResFT0CTPC"), centrality, QFT0C * QTPC * TMath::Cos(2.0 * (psiFT0C - psiTPC)));
    histos.fill(HIST("ResFT0CTPCR"), centrality, QFT0C * QTPCR * TMath::Cos(2.0 * (psiFT0C - psiTPCR)));
    histos.fill(HIST("ResFT0CTPCL"), centrality, QFT0C * QTPCL * TMath::Cos(2.0 * (psiFT0C - psiTPCL)));
    histos.fill(HIST("ResTPCRTPCL"), centrality, QTPCR * QTPCL * TMath::Cos(2.0 * (psiTPCR - psiTPCL)));
    histos.fill(HIST("ResFT0CFT0A"), centrality, QFT0C * QFT0A * TMath::Cos(2.0 * (psiFT0C - psiFT0A)));
    histos.fill(HIST("ResFT0ATPC"), centrality, QTPC * QFT0A * TMath::Cos(2.0 * (psiTPC - psiFT0A)));
    histos.fill(HIST("hCentrality"), centrality);
    histos.fill(HIST("hVtxZ"), collision.posZ());
    histos.fill(HIST("hOccupancy"), occupancy);
    auto firstprimarytrack = 0;
    for (auto track1 : tracks) {
      if (!selectionTrack(track1)) {
        continue;
      }
      // PID check
      if (ispTdepPID && !selectionPIDNew(track1)) {
        continue;
      }
      if (!ispTdepPID && !selectionPID(track1)) {
        continue;
      }
      if (!ispTdepPID && rejectPID && !RejectPion(track1)) {
        histos.fill(HIST("hRejectPID"), 0.5);
        continue;
      }
      if (!ispTdepPID && rejectPID && !RejectKaon(track1)) {
        histos.fill(HIST("hRejectPID"), 1.5);
        continue;
      }
      histos.fill(HIST("hMomCorr"), track1.p() / track1.sign(), track1.p() - track1.tpcInnerParam(), centrality);
      histos.fill(HIST("hEta"), track1.eta());
      histos.fill(HIST("hDcaxy"), track1.dcaXY());
      histos.fill(HIST("hDcaz"), track1.dcaZ());
      histos.fill(HIST("hNsigmaProtonTPCDiff"), track1.tpcNSigmaPr(), track1.tpcNSigmaKa(), track1.pt());
      histos.fill(HIST("hNsigmaProtonTPC"), track1.tpcNSigmaPr(), track1.pt());
      histos.fill(HIST("hNsigmaProtonTOF"), track1.tofNSigmaPr(), track1.pt());
      auto track1ID = track1.globalIndex();
      for (auto v0 : V0s) {
        if (firstprimarytrack == 0) {
          histos.fill(HIST("hV0Dca"), v0.dcav0topv());
        }
        if (!SelectionV0(collision, v0)) {
          continue;
        }
        auto postrack = v0.template posTrack_as<AllTrackCandidates>();
        auto negtrack = v0.template negTrack_as<AllTrackCandidates>();
        if (!isSelectedV0Daughter(postrack, 1)) {
          continue;
        }
        if (!isSelectedV0Daughter(negtrack, -1)) {
          continue;
        }
        if (track1ID == postrack.globalIndex()) {
          continue;
        }
        if (track1ID == negtrack.globalIndex()) {
          continue;
        }
        if (firstprimarytrack == 0) {
          histos.fill(HIST("hInvMassKs0"), v0.mK0Short());
        }
        firstprimarytrack = firstprimarytrack + 1;
        Proton = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massPr);
        Kshort = ROOT::Math::PxPyPzMVector(v0.px(), v0.py(), v0.pz(), massK0s);
        Lambdac = Proton + Kshort;
        if (TMath::Abs(Lambdac.Rapidity()) > confRapidity) {
          continue;
        }
        auto phiminuspsi = GetPhiInRange(Lambdac.Phi() - psiFT0C);
        auto v2 = TMath::Cos(2.0 * phiminuspsi) * QFT0C;

        auto diffangle = Proton.Phi() - Lambdac.Phi();
        auto decaylength = std::abs((track1.dcaXY() / TMath::Sin(diffangle)) / (Lambdac.P() / 2.286));
        if (fillDefault && Lambdac.M() > 2.15 && Lambdac.M() <= 2.45) {
          if (fillDecayLength) {
            histos.fill(HIST("hSparseV2SASameEvent_V2"), Lambdac.M(), Lambdac.Pt(), v2, decaylength, Proton.Pt(), centrality);
          }
          histos.fill(HIST("hSparseV2SASameEvent_V2_new"), Lambdac.M(), Lambdac.Pt(), v2, std::abs(track1.dcaXY()), Proton.Pt(), centrality);
        }
        if (fillOccupancy && Lambdac.M() > 2.15 && Lambdac.M() <= 2.45) {
          if (fillDecayLength) {
            histos.fill(HIST("hSparseV2SASameEvent_V2_occupancy"), Lambdac.M(), Lambdac.Pt(), v2, decaylength, Proton.Pt(), occupancy);
          }
          histos.fill(HIST("hSparseV2SASameEvent_V2_new_occupancy"), Lambdac.M(), Lambdac.Pt(), v2, std::abs(track1.dcaXY()), Proton.Pt(), occupancy);
        }
        if (fillRotation) {
          for (int nrotbkg = 1; nrotbkg <= nBkgRotations; nrotbkg++) {
            auto anglestep = nrotbkg * (TMath::Pi() / nBkgRotations);
            auto rotProtonPx = track1.px() * std::cos(anglestep) - track1.py() * std::sin(anglestep);
            auto rotProtonPy = track1.px() * std::sin(anglestep) + track1.py() * std::cos(anglestep);
            ProtonRot = ROOT::Math::PxPyPzMVector(rotProtonPx, rotProtonPy, track1.pz(), massPr);
            LambdacRot = ProtonRot + Kshort;
            auto phiminuspsiRot = GetPhiInRange(LambdacRot.Phi() - psiFT0C);
            auto v2Rot = TMath::Cos(2.0 * phiminuspsiRot) * QFT0C;
            auto diffangleRot = ProtonRot.Phi() - LambdacRot.Phi();
            auto decaylengthRot = std::abs((track1.dcaXY() / TMath::Sin(diffangleRot)) / (Lambdac.P() / 2.286));
            if (fillDefault && LambdacRot.M() > 2.15 && LambdacRot.M() <= 2.45) {
              if (fillDecayLength) {
                histos.fill(HIST("hSparseV2SASameEventRotational_V2"), LambdacRot.M(), LambdacRot.Pt(), v2Rot, decaylengthRot, Proton.Pt());
              }
              histos.fill(HIST("hSparseV2SASameEventRotational_V2_new"), LambdacRot.M(), LambdacRot.Pt(), v2Rot, std::abs(track1.dcaXY()), Proton.Pt());
            }
            if (fillOccupancy && LambdacRot.M() > 2.15 && LambdacRot.M() <= 2.45) {
              if (fillDecayLength) {
                histos.fill(HIST("hSparseV2SASameEventRotational_V2_occupancy"), LambdacRot.M(), LambdacRot.Pt(), v2Rot, decaylengthRot, Proton.Pt(), occupancy);
              }
              histos.fill(HIST("hSparseV2SASameEventRotational_V2_new_occupancy"), LambdacRot.M(), LambdacRot.Pt(), v2Rot, std::abs(track1.dcaXY()), Proton.Pt(), occupancy);
            }
          }
        }
        ROOT::Math::Boost boost{Lambdac.BoostToCM()};
        fourVecDauCM = boost(Kshort);
        threeVecDauCM = fourVecDauCM.Vect();
        threeVecDauCMXY = ROOT::Math::XYZVector(threeVecDauCM.X(), threeVecDauCM.Y(), 0.);
        beamvector = ROOT::Math::XYZVector(0, 0, 1);
        auto cosThetaStar = beamvector.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(beamvector.Mag2());
        auto SA = cosThetaStar * TMath::Sin(2.0 * phiminuspsi);
        if (fillPolarization) {
          if (track1.sign() > 0) {
            histos.fill(HIST("hSparseV2SASameEventplus_SA"), Lambdac.M(), Lambdac.Pt(), cosThetaStar, phiminuspsi, centrality);
            histos.fill(HIST("hSparseV2SASameEventplus_SA_A0"), Lambdac.M(), Lambdac.Pt(), cosThetaStar * cosThetaStar, phiminuspsi, centrality);
            histos.fill(HIST("hSparseV2SASameEventplus_SA_azimuth"), Lambdac.M(), Lambdac.Pt(), SA, centrality);
          }
          if (track1.sign() < 0) {
            histos.fill(HIST("hSparseV2SASameEventminus_SA"), Lambdac.M(), Lambdac.Pt(), cosThetaStar, phiminuspsi, centrality);
            histos.fill(HIST("hSparseV2SASameEventminus_SA_A0"), Lambdac.M(), Lambdac.Pt(), cosThetaStar * cosThetaStar, phiminuspsi, centrality);
            histos.fill(HIST("hSparseV2SASameEventminus_SA_azimuth"), Lambdac.M(), Lambdac.Pt(), SA, centrality);
          }
        }
      }
    }
  }
  PROCESS_SWITCH(highmasslambda, processSameEvent, "Process Same event", true);
  void processMixedEventOpti(EventCandidates const& collisions, TrackCandidates const& tracks, AllTrackCandidates const&, ResoV0s const& V0s)
  {
    auto tracksV0sTuple = std::make_tuple(tracks, V0s);
    BinningTypeVertexContributor binningOnPositions{{axisVertex, axisMultiplicityClass, axisEPAngle}, true};
    Pair<EventCandidates, TrackCandidates, ResoV0s, BinningTypeVertexContributor> pairs{binningOnPositions, cfgNoMixedEvents, -1, collisions, tracksV0sTuple, &cache}; // -1 is the number of the bin to skip
    for (auto& [collision1, tracks1, collision2, tracks2] : pairs) {
      if (!collision1.sel8() || !collision2.sel8()) {
        continue;
      }
      if (!collision1.triggereventep() || !collision2.triggereventep()) {
        continue;
      }
      if (additionalEvSel && (!collision1.selection_bit(aod::evsel::kNoSameBunchPileup) || !collision1.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV))) {
        continue;
      }
      if (additionalEvSel && (!collision2.selection_bit(aod::evsel::kNoSameBunchPileup) || !collision2.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV))) {
        continue;
      }
      auto centrality = collision1.centFT0C();
      auto psiFT0C = collision1.psiFT0C();
      auto QFT0C = collision1.qFT0C();
      int occupancy = collision1.trackOccupancyInTimeRange();
      for (auto& [track1, v0] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
        if (!selectionTrack(track1)) {
          continue;
        }
        // PID check
        if (ispTdepPID && !selectionPIDNew(track1)) {
          continue;
        }
        if (!ispTdepPID && !selectionPID(track1)) {
          continue;
        }
        if (!SelectionV0(collision2, v0)) {
          continue;
        }
        auto postrack = v0.template posTrack_as<AllTrackCandidates>();
        auto negtrack = v0.template negTrack_as<AllTrackCandidates>();
        if (!isSelectedV0Daughter(postrack, 1)) {
          continue;
        }
        if (!isSelectedV0Daughter(negtrack, -1)) {
          continue;
        }

        Proton = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massPr);
        Kshort = ROOT::Math::PxPyPzMVector(v0.px(), v0.py(), v0.pz(), massK0s);
        Lambdac = Proton + Kshort;
        if (TMath::Abs(Lambdac.Rapidity()) > confRapidity) {
          continue;
        }
        auto phiminuspsi = GetPhiInRange(Lambdac.Phi() - psiFT0C);
        auto v2 = TMath::Cos(2.0 * phiminuspsi) * QFT0C;
        auto diffangle = Proton.Phi() - Lambdac.Phi();
        auto decaylength = std::abs((track1.dcaXY() / TMath::Sin(diffangle)) / (Lambdac.P() / 2.286));
        if (fillDefault && Lambdac.M() > 2.15 && Lambdac.M() <= 2.45) {
          if (fillDecayLength) {
            histos.fill(HIST("hSparseV2SAMixedEvent_V2"), Lambdac.M(), Lambdac.Pt(), v2, decaylength, Proton.Pt());
          }
          histos.fill(HIST("hSparseV2SAMixedEvent_V2_new"), Lambdac.M(), Lambdac.Pt(), v2, std::abs(track1.dcaXY()), Proton.Pt());
        }
        if (fillOccupancy && Lambdac.M() > 2.15 && Lambdac.M() <= 2.45) {
          if (fillDecayLength) {
            histos.fill(HIST("hSparseV2SAMixedEvent_V2_occupancy"), Lambdac.M(), Lambdac.Pt(), v2, decaylength, Proton.Pt(), occupancy);
          }
          histos.fill(HIST("hSparseV2SAMixedEvent_V2_new_occupancy"), Lambdac.M(), Lambdac.Pt(), v2, std::abs(track1.dcaXY()), Proton.Pt(), occupancy);
        }
        ROOT::Math::Boost boost{Lambdac.BoostToCM()};
        fourVecDauCM = boost(Kshort);
        threeVecDauCM = fourVecDauCM.Vect();
        threeVecDauCMXY = ROOT::Math::XYZVector(threeVecDauCM.X(), threeVecDauCM.Y(), 0.);
        beamvector = ROOT::Math::XYZVector(0, 0, 1);
        auto cosThetaStar = beamvector.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(beamvector.Mag2());
        auto SA = cosThetaStar * TMath::Sin(2.0 * phiminuspsi);
        if (fillPolarization) {
          if (track1.sign() > 0) {
            histos.fill(HIST("hSparseV2SAMixedEventplus_SA"), Lambdac.M(), Lambdac.Pt(), cosThetaStar, phiminuspsi, centrality);
            histos.fill(HIST("hSparseV2SAMixedEventplus_SA_A0"), Lambdac.M(), Lambdac.Pt(), cosThetaStar * cosThetaStar, phiminuspsi, centrality);
            histos.fill(HIST("hSparseV2SAMixedEventplus_SA_azimuth"), Lambdac.M(), Lambdac.Pt(), SA, centrality);
          }
          if (track1.sign() < 0) {
            histos.fill(HIST("hSparseV2SAMixedEventminus_SA"), Lambdac.M(), Lambdac.Pt(), cosThetaStar, phiminuspsi, centrality);
            histos.fill(HIST("hSparseV2SAMixedEventminus_SA_A0"), Lambdac.M(), Lambdac.Pt(), cosThetaStar * cosThetaStar, phiminuspsi, centrality);
            histos.fill(HIST("hSparseV2SAMixedEventminus_SA_azimuth"), Lambdac.M(), Lambdac.Pt(), SA, centrality);
          }
        }
      }
    }
  }
  PROCESS_SWITCH(highmasslambda, processMixedEventOpti, "Process Mixed event new", true);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<highmasslambda>(cfgc, TaskName{"highmasslambda"})};
}
