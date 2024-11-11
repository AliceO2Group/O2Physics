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
#include <vector>

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
  // Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  // Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  // Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  // Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  // Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};

  // fill output
  Configurable<bool> useSP{"useSP", false, "useSP"};
  Configurable<bool> useSignDCAV0{"useSignDCAV0", true, "useSignDCAV0"};
  Configurable<bool> fillDefault{"fillDefault", false, "fill Occupancy"};
  Configurable<int> cfgOccupancyCut{"cfgOccupancyCut", 2500, "Occupancy cut"};
  Configurable<bool> fillDecayLength{"fillDecayLength", true, "fill decay length"};
  Configurable<bool> fillPolarization{"fillPolarization", false, "fill polarization"};
  Configurable<bool> fillRotation{"fillRotation", false, "fill rotation"};
  // events
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCutCentralityMax{"cfgCutCentralityMax", 50.0f, "Accepted maximum Centrality"};
  Configurable<float> cfgCutCentralityMin{"cfgCutCentralityMin", 30.0f, "Accepted minimum Centrality"};
  // proton track cut
  Configurable<bool> rejectPID{"reject PID", true, "pion, kaon, electron rejection"};
  Configurable<float> cfgCutTOFBeta{"cfgCutTOFBeta", 0.0, "cut TOF beta"};
  Configurable<bool> ispTdifferentialDCA{"ispTdifferentialDCA", true, "is pT differential DCA"};
  Configurable<bool> isPVContributor{"isPVContributor", true, "is PV contributor"};
  Configurable<float> confMinRot{"confMinRot", 5.0 * TMath::Pi() / 6.0, "Minimum of rotation"};
  Configurable<float> confMaxRot{"confMaxRot", 7.0 * TMath::Pi() / 6.0, "Maximum of rotation"};
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
  Configurable<int> PIDstrategy{"PIDstrategy", 0, "0: TOF Veto, 1: TOF Veto opti, 2: TOF, 3: TOF loose 1, 4: TOF loose 2, 5: old pt dep"};
  Configurable<float> nsigmaCutTPC{"nsigmacutTPC", 3.0, "Value of the TPC Nsigma cut"};
  Configurable<float> nsigmaCutTOF{"nsigmaCutTOF", 3.0, "TOF PID"};
  Configurable<float> nsigmaCutTPCPre{"nsigmacutTPCPre", 3.0, "Value of the TPC Nsigma cut Pre filter"};
  Configurable<float> kaonrejpar{"kaonrejpar", 1.0, "Kaon rej. par"};
  // Configs for V0
  Configurable<float> ConfV0PtMin{"ConfV0PtMin", 0.f, "Minimum transverse momentum of V0"};
  Configurable<double> ConfV0DCADaughMax{"ConfV0DCADaughMax", 0.2f, "Maximum DCA between the V0 daughters"};
  Configurable<double> ConfV0CPAMin{"ConfV0CPAMin", 0.9998f, "Minimum CPA of V0"};
  Configurable<float> ConfV0TranRadV0Min{"ConfV0TranRadV0Min", 1.5f, "Minimum transverse radius"};
  Configurable<float> ConfV0TranRadV0Max{"ConfV0TranRadV0Max", 100.f, "Maximum transverse radius"};
  Configurable<double> cMaxV0DCA{"cMaxV0DCA", 0.2, "Maximum V0 DCA to PV"};
  Configurable<float> cMaxV0LifeTime{"cMaxV0LifeTime", 20, "Maximum V0 life time"};
  Configurable<float> cSigmaMassKs0{"cSigmaMassKs0", 0.006, "Sigma cut on KS0 mass"};
  Configurable<float> cMinLambdaMass{"cMinLambdaMass", 2.18, "Minimum lambda mass"};
  Configurable<float> cMaxLambdaMass{"cMaxLambdaMass", 2.42, "Maximum lambda mass"};
  // config for V0 daughters
  Configurable<float> ConfDaughEta{"ConfDaughEta", 0.8f, "V0 Daugh sel: max eta"};
  Configurable<float> ConfDaughPt{"ConfDaughPt", 0.1f, "V0 Daugh sel: min pt"};
  Configurable<float> ConfDaughTPCnclsMin{"ConfDaughTPCnclsMin", 50.f, "V0 Daugh sel: Min. nCls TPC"};
  Configurable<double> ConfDaughDCAMin{"ConfDaughDCAMin", 0.08f, "V0 Daugh sel:  Max. DCA Daugh to PV (cm)"};
  Configurable<float> ConfDaughPIDCuts{"ConfDaughPIDCuts", 3, "PID selections for KS0 daughters"};
  // Fill strategy
  // Configurable<int> cfgSelectDaughterTopology{"cfgSelectDaughterTopology", 2, "Select daughter for topology"};
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
  ConfigurableAxis cnfigThnAxisDCASum{"cnfigThnAxisDCASum", {150, -0.3, 0.3}, "SA"};
  ConfigurableAxis cnfigThnAxisPtProton{"cnfigThnAxisPtProton", {16, 0.0, 8.0}, "pT"};
  ConfigurableAxis configThnAxisSP{"configThnAxisSP", {400, -12, 12}, "SP"};

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter centralityFilter = (nabs(aod::cent::centFT0C) < cfgCutCentralityMax && nabs(aod::cent::centFT0C) > cfgCutCentralityMin);
  Filter acceptanceFilter = (nabs(aod::track::eta) < cfgCutEta && nabs(aod::track::pt) > cfgCutPT);
  Filter dcaCutFilter = (nabs(aod::track::dcaXY) < cfgCutDCAxy) && (nabs(aod::track::dcaZ) < cfgCutDCAz);
  Filter pidFilter = nabs(aod::pidtpc::tpcNSigmaPr) < nsigmaCutTPCPre;

  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::FV0Mults, aod::TPCMults, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0As, aod::EPCalibrationTables, aod::Mults>>;
  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPr, aod::pidTOFbeta, aod::pidTOFFullPr, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullEl, aod::pidTOFFullEl>>;
  using AllTrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPi>;
  using ResoV0s = aod::V0Datas;

  SliceCache cache;
  // Partition<TrackCandidates> posTracks = aod::track::signed1Pt > cfgCutCharge;
  // Partition<TrackCandidates> negTracks = aod::track::signed1Pt < cfgCutCharge;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext&)
  {
    std::vector<double> occupancyBinning = {0.0, 500.0, 1000.0, 1500.0, 2000.0, 3000.0, 4000.0, 5000.0, 50000.0};
    // std::vector<double> dcaBinning = {0.0, 0.0005, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01, 0.012, 0.014, 0.016, 0.02, 0.03, 0.05, 0.1, 0.5, 1.0};
    std::vector<double> dcaBinning = {0.0, 0.0005, 0.001, 0.0012, 0.0014, 0.0016, 0.002, 0.0025, 0.003, 0.004, 0.005, 0.006, 0.008, 0.01, 0.015, 0.02, 0.04, 0.05, 0.06, 0.08, 0.1, 0.3, 1.0};
    std::vector<double> ptProtonBinning = {0.2, 0.3, 0.5, 0.6, 0.8, 1.2, 1.4, 1.6, 2.0, 3.0, 4.0, 6.0};
    std::vector<double> ptLambdaBinning = {2.0, 3.0, 4.0, 5.0, 6.0};

    const AxisSpec thnAxisInvMass{configThnAxisInvMass, "#it{M} (GeV/#it{c}^{2})"};
    const AxisSpec thnAxisPt{configThnAxisPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec thnAxisCosThetaStar{configThnAxisCosThetaStar, "cos(#vartheta)"};
    const AxisSpec thnAxisPhiminusPsi{configThnAxisPhiminusPsi, "#phi - #psi"};
    const AxisSpec thnAxisCentrality{configThnAxisCentrality, "Centrality (%)"};
    const AxisSpec thnAxisV2{configThnAxisV2, "V2"};
    const AxisSpec thnAxisSA{configThnAxisSA, "SA"};
    const AxisSpec thnAxisDecayLength{cnfigThnAxisDecayLength, "Decay Length"};
    const AxisSpec thnAxisDCASum{cnfigThnAxisDCASum, "DCA Sum"};
    const AxisSpec thnAxisPtProton{cnfigThnAxisPtProton, "Proton Pt"};
    const AxisSpec thnAxisSP{configThnAxisSP, "SP"};

    AxisSpec phiAxis = {500, -6.28, 6.28, "phi"};
    AxisSpec resAxis = {1000, -10, 10, "Res"};
    AxisSpec centAxis = {8, 0, 80, "V0M (%)"};
    // AxisSpec ptProtonAxis = {16, 0.0, 8, "V0M (%)"};
    AxisSpec dcaAxis = {dcaBinning, "DCAxy"};
    AxisSpec dcatoPVAxis = {50, 0.0, 0.5, "V0M (%)"};
    AxisSpec occupancyAxis = {occupancyBinning, "occupancy"};
    AxisSpec ptProtonAxis = {ptProtonBinning, "pT proton"};

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
    histos.add("hRotation", "hRotation", kTH1F, {{360, 0.0, 2.0 * TMath::Pi()}});
    histos.add("hEta", "Eta distribution", kTH1F, {{200, -1.0f, 1.0f}});
    histos.add("hDcaxy", "Dcaxy distribution", kTH1F, {{1000, -0.5f, 0.5f}});
    histos.add("hDcaz", "Dcaz distribution", kTH1F, {{1000, -0.5f, 0.5f}});
    histos.add("hNsigmaProtonTPCDiff", "Difference NsigmaProton NsigmaKaon TPC distribution", kTH3F, {{100, -5.0f, 5.0f}, {100, -5.0f, 5.0f}, {80, 0.0f, 8.0f}});

    histos.add("hNsigmaProtonElectronTPC", "NsigmaProton-Electron TPC distribution", kTH3F, {{60, -3.0f, 3.0f}, {200, -10.0f, 10.0f}, {60, 0.0f, 6.0f}});
    histos.add("hNsigmaProtonPionTPC", "NsigmaProton-Pion TPC distribution", kTH3F, {{60, -3.0f, 3.0f}, {200, -10.0f, 10.0f}, {60, 0.0f, 6.0f}});
    histos.add("hNsigmaProtonKaonTPC", "NsigmaProton-Kaon TPC distribution", kTH3F, {{60, -3.0f, 3.0f}, {200, -10.0f, 10.0f}, {60, 0.0f, 6.0f}});

    histos.add("hNsigmaProtonElectronTPC_afterPi", "NsigmaProton-Electron TPC distribution", kTH3F, {{60, -3.0f, 3.0f}, {200, -10.0f, 10.0f}, {60, 0.0f, 6.0f}});
    histos.add("hNsigmaProtonPionTPC_afterPi", "NsigmaProton-Pion TPC distribution", kTH3F, {{60, -3.0f, 3.0f}, {200, -10.0f, 10.0f}, {60, 0.0f, 6.0f}});
    histos.add("hNsigmaProtonKaonTPC_afterPi", "NsigmaProton-Kaon TPC distribution", kTH3F, {{60, -3.0f, 3.0f}, {200, -10.0f, 10.0f}, {60, 0.0f, 6.0f}});

    histos.add("hNsigmaProtonElectronTPC_afterEl", "NsigmaProton-Electron TPC distribution", kTH3F, {{60, -3.0f, 3.0f}, {200, -10.0f, 10.0f}, {60, 0.0f, 6.0f}});
    histos.add("hNsigmaProtonPionTPC_afterEl", "NsigmaProton-Pion TPC distribution", kTH3F, {{60, -3.0f, 3.0f}, {200, -10.0f, 10.0f}, {60, 0.0f, 6.0f}});
    histos.add("hNsigmaProtonKaonTPC_afterEl", "NsigmaProton-Kaon TPC distribution", kTH3F, {{60, -3.0f, 3.0f}, {200, -10.0f, 10.0f}, {60, 0.0f, 6.0f}});

    histos.add("hNsigmaProtonElectronTPC_afterKa", "NsigmaProton-Electron TPC distribution", kTH3F, {{60, -3.0f, 3.0f}, {200, -10.0f, 10.0f}, {60, 0.0f, 6.0f}});
    histos.add("hNsigmaProtonPionTPC_afterKa", "NsigmaProton-Pion TPC distribution", kTH3F, {{60, -3.0f, 3.0f}, {200, -10.0f, 10.0f}, {60, 0.0f, 6.0f}});
    histos.add("hNsigmaProtonKaonTPC_afterKa", "NsigmaProton-Kaon TPC distribution", kTH3F, {{60, -3.0f, 3.0f}, {200, -10.0f, 10.0f}, {60, 0.0f, 6.0f}});

    histos.add("hNsigmaProtonTPC", "NsigmaProton TPC distribution", kTH3F, {{100, -5.0f, 5.0f}, {60, 0.0f, 6.0f}, occupancyBinning});
    histos.add("hNsigmaProtonTOF", "NsigmaProton TOF distribution", kTH2F, {{1000, -50.0f, 50.0f}, {60, 0.0f, 6.0f}});
    histos.add("hNsigmaProtonTOFPre", "NsigmaProton TOF distribution Pre sel", kTH2F, {{1000, -50.0f, 50.0f}, {60, 0.0f, 6.0f}});
    histos.add("hMassvsDecaySum", "hMassvsDecaySum", kTH2F, {thnAxisInvMass, thnAxisDCASum});
    histos.add("hPsiFT0C", "PsiFT0C", kTH3F, {centAxis, phiAxis, occupancyAxis});
    histos.add("hPsiFT0A", "PsiFT0A", kTH3F, {centAxis, phiAxis, occupancyAxis});
    histos.add("hPsiTPC", "PsiTPC", kTH3F, {centAxis, phiAxis, occupancyAxis});
    histos.add("hPsiTPCR", "PsiTPCR", kTH3F, {centAxis, phiAxis, occupancyAxis});
    histos.add("hPsiTPCL", "PsiTPCL", kTH3F, {centAxis, phiAxis, occupancyAxis});

    if (fillDefault) {
      histos.add("hSparseV2SASameEvent_V2", "hSparseV2SASameEvent_V2", HistType::kTHnSparseF, {thnAxisInvMass, ptLambdaBinning, thnAxisSP, dcaAxis, ptProtonAxis});
      histos.add("hSparseV2SASameEventRotational_V2", "hSparseV2SASameEventRotational", HistType::kTHnSparseF, {thnAxisInvMass, ptLambdaBinning, thnAxisSP, dcaAxis, ptProtonAxis});
      histos.add("hSparseV2SAMixedEvent_V2", "hSparseV2SAMixedEvent_V2", HistType::kTHnSparseF, {thnAxisInvMass, ptLambdaBinning, thnAxisSP, dcaAxis, ptProtonAxis});
    }
    if (fillDecayLength) {
      histos.add("hSparseV2SASameEvent_V2_dcatopv", "hSparseV2SASameEvent_V2_dcatopv", HistType::kTHnSparseF, {thnAxisInvMass, ptLambdaBinning, thnAxisSP, dcatoPVAxis, ptProtonAxis});
      histos.add("hSparseV2SASameEventRotational_V2_dcatopv", "hSparseV2SASameEventRotational_V2_dcatopv", HistType::kTHnSparseF, {thnAxisInvMass, ptLambdaBinning, thnAxisSP, dcatoPVAxis, ptProtonAxis});
      histos.add("hSparseV2SAMixedEvent_V2_dcatopv", "hSparseV2SAMixedEvent_V2_dcatopv", HistType::kTHnSparseF, {thnAxisInvMass, ptLambdaBinning, thnAxisSP, dcatoPVAxis, ptProtonAxis});
    }
    if (fillPolarization) {
      histos.add("hSparseV2SASameEventplus_SA", "hSparseV2SASameEventplus_SA", HistType::kTHnSparseF, {thnAxisInvMass, ptLambdaBinning, thnAxisSA, thnAxisPhiminusPsi, thnAxisCentrality});
      histos.add("hSparseV2SASameEventplus_SA_A0", "hSparseV2SASameEventplus_SA_A0", HistType::kTHnSparseF, {thnAxisInvMass, ptLambdaBinning, thnAxisSA, thnAxisPhiminusPsi, thnAxisCentrality});
      histos.add("hSparseV2SASameEventplus_SA_azimuth", "hSparseV2SASameEventplus_SA_azimuth", HistType::kTHnSparseF, {thnAxisInvMass, ptLambdaBinning, thnAxisSA, thnAxisCentrality});
      histos.add("hSparseV2SASameEventminus_SA", "hSparseV2SASameEventminus_SA", HistType::kTHnSparseF, {thnAxisInvMass, ptLambdaBinning, thnAxisSA, thnAxisPhiminusPsi, thnAxisCentrality});
      histos.add("hSparseV2SASameEventminus_SA_A0", "hSparseV2SASameEventminus_SA_A0", HistType::kTHnSparseF, {thnAxisInvMass, ptLambdaBinning, thnAxisSA, thnAxisPhiminusPsi, thnAxisCentrality});
      histos.add("hSparseV2SASameEventminus_SA_azimuth", "hSparseV2SASameEventminus_SA_azimuth", HistType::kTHnSparseF, {thnAxisInvMass, ptLambdaBinning, thnAxisSA, thnAxisCentrality});

      histos.add("hSparseV2SAMixedEventplus_SA", "hSparseV2SAMixedEventplus_SA", HistType::kTHnSparseF, {thnAxisInvMass, ptLambdaBinning, thnAxisSA, thnAxisPhiminusPsi, thnAxisCentrality});
      histos.add("hSparseV2SAMixedEventplus_SA_A0", "hSparseV2SAMixedEventplus_SA_A0", HistType::kTHnSparseF, {thnAxisInvMass, ptLambdaBinning, thnAxisSA, thnAxisPhiminusPsi, thnAxisCentrality});
      histos.add("hSparseV2SAMixedEventplus_SA_azimuth", "hSparseV2SAMixedEventplus_SA_azimuth", HistType::kTHnSparseF, {thnAxisInvMass, ptLambdaBinning, thnAxisSA, thnAxisCentrality});
      histos.add("hSparseV2SAMixedEventminus_SA", "hSparseV2SAMixedEventminus_SA", HistType::kTHnSparseF, {thnAxisInvMass, ptLambdaBinning, thnAxisSA, thnAxisPhiminusPsi, thnAxisCentrality});
      histos.add("hSparseV2SAMixedEventminus_SA_A0", "hSparseV2SAMixedEventminus_SA_A0", HistType::kTHnSparseF, {thnAxisInvMass, ptLambdaBinning, thnAxisSA, thnAxisPhiminusPsi, thnAxisCentrality});
      histos.add("hSparseV2SAMixedEventminus_SA_azimuth", "hSparseV2SAMixedEventminus_SA_azimuth", HistType::kTHnSparseF, {thnAxisInvMass, ptLambdaBinning, thnAxisSA, thnAxisCentrality});
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
    if (isPVContributor && !(candidate.isGlobalTrack() && candidate.isPVContributor() && candidate.itsNCls() > cfgITScluster && candidate.tpcNClsFound() > cfgTPCcluster && candidate.itsNClsInnerBarrel() >= 1)) {
      return false;
    }
    if (!isPVContributor && !(candidate.isGlobalTrackWoDCA() && candidate.itsNCls() > cfgITScluster && candidate.tpcNClsFound() > cfgTPCcluster && candidate.itsNClsInnerBarrel() >= 1)) {
      return false;
    }

    if (ispTdifferentialDCA) {
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
    }

    if (!ispTdifferentialDCA) {
      if (candidate.pt() > 0.0 && candidate.pt() < 0.5 && TMath::Abs(candidate.dcaXY()) < 0.004) {
        return false;
      }
      if (candidate.pt() >= 0.5 && candidate.pt() < 1.0 && TMath::Abs(candidate.dcaXY()) < 0.003) {
        return false;
      }
      if (candidate.pt() >= 1.0 && candidate.pt() < 2.0 && TMath::Abs(candidate.dcaXY()) < 0.002) {
        return false;
      }
    }
    return true;
  }

  template <typename T>
  bool rejectPi(const T& candidate)
  {
    if (candidate.tpcInnerParam() > 0.9 && candidate.tpcInnerParam() < 1.0 && candidate.tpcNSigmaPi() < 6.0) {
      return false;
    }
    if (candidate.tpcInnerParam() > 1.0 && candidate.tpcInnerParam() < 1.1 && candidate.tpcNSigmaPi() < 4.0) {
      return false;
    }
    if (candidate.tpcInnerParam() > 1.1 && candidate.tpcInnerParam() < 1.2 && candidate.tpcNSigmaPi() < 3.0) {
      return false;
    }
    if (candidate.tpcInnerParam() > 1.2 && candidate.tpcInnerParam() < 1.4 && candidate.tpcNSigmaPi() < 1.0) {
      return false;
    }
    if (candidate.tpcInnerParam() > 1.4 && candidate.tpcInnerParam() < 1.5 && candidate.tpcNSigmaPi() < 0.5) {
      return false;
    }
    return true;
  }

  template <typename T>
  bool rejectEl(const T& candidate)
  {

    if (candidate.tpcInnerParam() > 0.7 && candidate.tpcInnerParam() < 0.8 && candidate.tpcNSigmaEl() < 2.0) {
      return false;
    }
    if (candidate.tpcInnerParam() > 0.8 && candidate.tpcInnerParam() < 0.9 && candidate.tpcNSigmaEl() < 0.0) {
      return false;
    }
    if (candidate.tpcInnerParam() > 0.9 && candidate.tpcInnerParam() < 1.0 && candidate.tpcNSigmaEl() < -1.0) {
      return false;
    }
    if (candidate.tpcInnerParam() > 1.0 && candidate.tpcInnerParam() < 1.1 && candidate.tpcNSigmaEl() < -2.0) {
      return false;
    }
    if (candidate.tpcInnerParam() > 1.1 && candidate.tpcInnerParam() < 1.2 && candidate.tpcNSigmaEl() < -3.0) {
      return false;
    }

    return true;
  }

  template <typename T>
  bool rejectKa(const T& candidate)
  {
    if (candidate.tpcInnerParam() > 0.7 && candidate.tpcInnerParam() < 0.8 && candidate.tpcNSigmaKa() < 7.5) {
      return false;
    }
    if (candidate.tpcInnerParam() > 0.8 && candidate.tpcInnerParam() < 0.9 && candidate.tpcNSigmaKa() < 6.0) {
      return false;
    }
    if (candidate.tpcInnerParam() > 0.9 && candidate.tpcInnerParam() < 1.1 && candidate.tpcNSigmaKa() < 5.0) {
      return false;
    }
    if (candidate.tpcInnerParam() > 1.1 && candidate.tpcInnerParam() < 1.2 && candidate.tpcNSigmaKa() < 3.5) {
      return false;
    }
    if (candidate.tpcInnerParam() > 1.2 && candidate.tpcInnerParam() < 1.4 && candidate.tpcNSigmaKa() < 3.0) {
      return false;
    }
    if (candidate.tpcInnerParam() > 1.4 && candidate.tpcInnerParam() < 1.5 && candidate.tpcNSigmaKa() < 2.5) {
      return false;
    }
    return true;
  }

  // TPC TOF
  template <typename T>
  bool selectionPID1(const T& candidate)
  {
    if (candidate.tpcInnerParam() < 0.9 && candidate.tpcNSigmaPr() > -2.0 && candidate.tpcNSigmaPr() < nsigmaCutTPC) {
      return true;
    }
    if (candidate.tpcInnerParam() >= 0.9) {
      if (candidate.hasTOF()) {
        if (candidate.beta() > cfgCutTOFBeta && candidate.tpcNSigmaPr() > -2.0 && candidate.tpcNSigmaPr() < nsigmaCutTPC && TMath::Abs(candidate.tofNSigmaPr()) < nsigmaCutTOF) {
          return true;
        }
      }
      if (!candidate.hasTOF()) {
        if (candidate.tpcNSigmaPr() > -2.0 && candidate.tpcNSigmaPr() < nsigmaCutTPC) {
          return true;
        }
      }
    }
    return false;
  }

  // TOF Veto
  template <typename T>
  bool selectionPID2(const T& candidate)
  {
    if (candidate.tpcInnerParam() < 0.9 && candidate.tpcNSigmaPr() > -2.0 && candidate.tpcNSigmaPr() < nsigmaCutTPC) {
      return true;
    }
    if (candidate.tpcInnerParam() >= 0.9 && candidate.beta() > cfgCutTOFBeta && candidate.tpcNSigmaPr() > -2.0 && candidate.tpcNSigmaPr() < nsigmaCutTPC && TMath::Abs(candidate.tofNSigmaPr()) < nsigmaCutTOF) {
      return true;
    }
    return false;
  }

  // TOF veto loose
  template <typename T>
  bool selectionPID3(const T& candidate)
  {
    if (candidate.tpcInnerParam() < 0.9 && candidate.tpcNSigmaPr() > -2.0 && candidate.tpcNSigmaPr() < nsigmaCutTPC) {
      return true;
    }
    if (candidate.tpcInnerParam() >= 0.9) {
      if (candidate.hasTOF()) {
        if (candidate.beta() > cfgCutTOFBeta && candidate.tpcNSigmaPr() > -2.0 && candidate.tpcNSigmaPr() < nsigmaCutTPC && TMath::Abs(candidate.tofNSigmaPr()) < nsigmaCutTOF) {
          return true;
        }
      }
      if (!candidate.hasTOF()) {
        if (candidate.tpcInnerParam() < 1.5 && candidate.tpcNSigmaPr() > -2.0 && candidate.tpcNSigmaPr() < nsigmaCutTPC) {
          return true;
        }
      }
    }
    return false;
  }

  // TOF veto very loose
  template <typename T>
  bool selectionPID4(const T& candidate)
  {
    if (candidate.tpcInnerParam() < 0.9 && candidate.tpcNSigmaPr() > -2.0 && candidate.tpcNSigmaPr() < nsigmaCutTPC) {
      return true;
    }
    if (candidate.tpcInnerParam() >= 0.9) {
      if (candidate.hasTOF()) {
        if (candidate.beta() > cfgCutTOFBeta && candidate.tpcNSigmaPr() > -2.0 && candidate.tpcNSigmaPr() < nsigmaCutTPC && TMath::Abs(candidate.tofNSigmaPr()) < nsigmaCutTOF) {
          return true;
        }
      }
      if (!candidate.hasTOF()) {
        if (candidate.tpcInnerParam() < 1.5 && candidate.tpcInnerParam() < 1.8 && candidate.tpcNSigmaPr() > -2.0 && candidate.tpcNSigmaPr() < nsigmaCutTPC) {
          return true;
        }
      }
    }
    return false;
  }

  template <typename Collision, typename V0>
  bool SelectionV0(Collision const& collision, V0 const& candidate)
  {
    if (TMath::Abs(candidate.dcav0topv()) > cMaxV0DCA) {
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
    if (TMath::Abs(CtauK0s) > cMaxV0LifeTime || candidate.mK0Short() < lowmasscutks0 || candidate.mK0Short() > highmasscutks0) {
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
    if (TMath::Abs(eta) > ConfDaughEta) {
      return false;
    }
    if (TMath::Abs(pt) < ConfDaughPt) {
      return false;
    }
    if (tpcNClsF < ConfDaughTPCnclsMin) {
      return false;
    }
    if (TMath::Abs(dcaXY) < ConfDaughDCAMin) {
      return false;
    }
    if (TMath::Abs(track.tpcNSigmaPi()) > ConfDaughPIDCuts) {
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
  ROOT::Math::PxPyPzMVector Lambdac, Proton, Kshort, LambdacRot, KshortRot, fourVecDauCM;
  ROOT::Math::XYZVector threeVecDauCM, threeVecDauCMXY, eventplaneVec, eventplaneVecNorm, beamvector;
  double massPi = TDatabasePDG::Instance()->GetParticle(kPiPlus)->Mass();   // FIXME: Get from the common header
  double massPr = TDatabasePDG::Instance()->GetParticle(kProton)->Mass();   // FIXME: Get from the common header
  double massK0s = TDatabasePDG::Instance()->GetParticle(kK0Short)->Mass(); // FIXME: Get from the common header
  double v2, v2Rot;
  void processSameEvent(EventCandidates::iterator const& collision, TrackCandidates const& tracks, AllTrackCandidates const&, ResoV0s const& V0s, aod::BCs const&)
  {
    if (!collision.sel8() || !collision.triggereventep() || !collision.selection_bit(aod::evsel::kNoTimeFrameBorder) || !collision.selection_bit(aod::evsel::kNoITSROFrameBorder) || !collision.selection_bit(aod::evsel::kNoSameBunchPileup) || !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV) || !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      return;
    }

    auto centrality = collision.centFT0C();
    auto multTPC = collision.multNTracksPV();
    histos.fill(HIST("hFTOCvsTPCNoCut"), centrality, multTPC);

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
    histos.fill(HIST("hPsiFT0C"), centrality, psiFT0C, occupancy);
    histos.fill(HIST("hPsiFT0A"), centrality, psiFT0A, occupancy);
    histos.fill(HIST("hPsiTPC"), centrality, psiTPC, occupancy);
    histos.fill(HIST("hPsiTPCR"), centrality, psiTPCR, occupancy);
    histos.fill(HIST("hPsiTPCL"), centrality, psiTPCL, occupancy);
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
      if (track1.hasTOF()) {
        histos.fill(HIST("hNsigmaProtonTOFPre"), track1.tofNSigmaPr(), track1.pt());
      }
      if (!track1.hasTOF()) {
        histos.fill(HIST("hNsigmaProtonElectronTPC"), track1.tpcNSigmaPr(), track1.tpcNSigmaEl(), track1.tpcInnerParam());
        histos.fill(HIST("hNsigmaProtonPionTPC"), track1.tpcNSigmaPr(), track1.tpcNSigmaPi(), track1.tpcInnerParam());
        histos.fill(HIST("hNsigmaProtonKaonTPC"), track1.tpcNSigmaPr(), track1.tpcNSigmaKa(), track1.tpcInnerParam());
        if (rejectPID && !rejectPi(track1)) {
          continue;
        }
        histos.fill(HIST("hNsigmaProtonElectronTPC_afterPi"), track1.tpcNSigmaPr(), track1.tpcNSigmaEl(), track1.tpcInnerParam());
        histos.fill(HIST("hNsigmaProtonPionTPC_afterPi"), track1.tpcNSigmaPr(), track1.tpcNSigmaPi(), track1.tpcInnerParam());
        histos.fill(HIST("hNsigmaProtonKaonTPC_afterPi"), track1.tpcNSigmaPr(), track1.tpcNSigmaKa(), track1.tpcInnerParam());
        if (rejectPID && !rejectEl(track1)) {
          continue;
        }
        histos.fill(HIST("hNsigmaProtonElectronTPC_afterEl"), track1.tpcNSigmaPr(), track1.tpcNSigmaEl(), track1.tpcInnerParam());
        histos.fill(HIST("hNsigmaProtonPionTPC_afterEl"), track1.tpcNSigmaPr(), track1.tpcNSigmaPi(), track1.tpcInnerParam());
        histos.fill(HIST("hNsigmaProtonKaonTPC_afterEl"), track1.tpcNSigmaPr(), track1.tpcNSigmaKa(), track1.tpcInnerParam());
        if (rejectPID && !rejectKa(track1)) {
          continue;
        }
        histos.fill(HIST("hNsigmaProtonElectronTPC_afterKa"), track1.tpcNSigmaPr(), track1.tpcNSigmaEl(), track1.tpcInnerParam());
        histos.fill(HIST("hNsigmaProtonPionTPC_afterKa"), track1.tpcNSigmaPr(), track1.tpcNSigmaPi(), track1.tpcInnerParam());
        histos.fill(HIST("hNsigmaProtonKaonTPC_afterKa"), track1.tpcNSigmaPr(), track1.tpcNSigmaKa(), track1.tpcInnerParam());
      }

      // PID check
      if (PIDstrategy == 0 && !selectionPID1(track1)) {
        continue;
      }
      if (PIDstrategy == 1 && !selectionPID2(track1)) {
        continue;
      }
      if (PIDstrategy == 2 && !selectionPID3(track1)) {
        continue;
      }
      if (PIDstrategy == 3 && !selectionPID4(track1)) {
        continue;
      }

      histos.fill(HIST("hMomCorr"), track1.p() / track1.sign(), track1.p() - track1.tpcInnerParam(), centrality);
      histos.fill(HIST("hEta"), track1.eta());
      histos.fill(HIST("hDcaxy"), track1.dcaXY());
      histos.fill(HIST("hDcaz"), track1.dcaZ());
      histos.fill(HIST("hNsigmaProtonTPCDiff"), track1.tpcNSigmaPr(), track1.tpcNSigmaKa(), track1.pt());
      histos.fill(HIST("hNsigmaProtonTPC"), track1.tpcNSigmaPr(), track1.pt(), occupancy);
      if (track1.hasTOF()) {
        histos.fill(HIST("hNsigmaProtonTOF"), track1.tofNSigmaPr(), track1.pt());
      }
      auto track1ID = track1.globalIndex();
      for (auto v0 : V0s) {
        if (!SelectionV0(collision, v0)) {
          continue;
        }
        auto anglesign = (v0.x() - collision.posX()) * v0.px() + (v0.y() - collision.posY()) * v0.py() + (v0.z() - collision.posZ()) * v0.pz();
        anglesign = anglesign / TMath::Abs(anglesign);
        if (firstprimarytrack == 0) {
          histos.fill(HIST("hV0Dca"), anglesign * v0.dcav0topv());
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
        Kshort = ROOT::Math::PxPyPzMVector(v0.px(), v0.py(), v0.pz(), v0.mK0Short());
        Lambdac = Proton + Kshort;
        auto phiminuspsi = GetPhiInRange(Lambdac.Phi() - psiFT0C);

        if (useSP) {
          v2 = TMath::Cos(2.0 * phiminuspsi) * QFT0C;
        }
        if (!useSP) {
          v2 = TMath::Cos(2.0 * phiminuspsi);
        }
        auto dcasum = 0.0;
        if (useSignDCAV0) {
          dcasum = anglesign * (v0.dcav0topv()) - track1.dcaXY();
        }
        if (!useSignDCAV0) {
          dcasum = v0.dcav0topv() - track1.dcaXY();
        }
        histos.fill(HIST("hMassvsDecaySum"), Lambdac.M(), dcasum);
        if (occupancy < cfgOccupancyCut && Lambdac.M() > cMinLambdaMass && Lambdac.M() <= cMaxLambdaMass && TMath::Abs(Lambdac.Rapidity()) < confRapidity && Lambdac.Pt() > 2.0 && Lambdac.Pt() <= 6.0) {
          if (fillDefault) {
            histos.fill(HIST("hSparseV2SASameEvent_V2"), Lambdac.M(), Lambdac.Pt(), v2, TMath::Abs(track1.dcaXY()), Proton.Pt());
          }
          if (fillDecayLength) {
            histos.fill(HIST("hSparseV2SASameEvent_V2_dcatopv"), Lambdac.M(), Lambdac.Pt(), v2, v0.dcav0topv(), Proton.Pt());
          }
        }

        if (fillRotation) {
          for (int nrotbkg = 0; nrotbkg < nBkgRotations; nrotbkg++) {
            auto anglestart = confMinRot;
            auto angleend = confMaxRot;
            auto anglestep = (angleend - anglestart) / (1.0 * (nBkgRotations - 1));
            auto rotangle = anglestart + nrotbkg * anglestep;
            histos.fill(HIST("hRotation"), rotangle);
            auto rotKaonPx = Kshort.px() * std::cos(rotangle) - Kshort.py() * std::sin(rotangle);
            auto rotKaonPy = Kshort.px() * std::sin(rotangle) + Kshort.py() * std::cos(rotangle);
            KshortRot = ROOT::Math::PxPyPzMVector(rotKaonPx, rotKaonPy, Kshort.pz(), massK0s);
            LambdacRot = Proton + KshortRot;
            auto phiminuspsiRot = GetPhiInRange(LambdacRot.Phi() - psiFT0C);
            if (useSP) {
              v2Rot = TMath::Cos(2.0 * phiminuspsiRot) * QFT0C;
            }
            if (!useSP) {
              v2Rot = TMath::Cos(2.0 * phiminuspsiRot);
            }

            if (occupancy < cfgOccupancyCut && LambdacRot.M() > cMinLambdaMass && LambdacRot.M() <= cMaxLambdaMass && TMath::Abs(LambdacRot.Rapidity()) < confRapidity && LambdacRot.Pt() > 2.0 && LambdacRot.Pt() <= 6.0) {
              if (fillDefault) {
                histos.fill(HIST("hSparseV2SASameEventRotational_V2"), LambdacRot.M(), LambdacRot.Pt(), v2Rot, TMath::Abs(track1.dcaXY()), Proton.Pt());
              }
              if (fillDecayLength) {
                histos.fill(HIST("hSparseV2SASameEventRotational_V2_dcatopv"), LambdacRot.M(), LambdacRot.Pt(), v2Rot, v0.dcav0topv(), Proton.Pt());
              }
            }
          }
        }
        if (fillPolarization) {
          ROOT::Math::Boost boost{Lambdac.BoostToCM()};
          fourVecDauCM = boost(Kshort);
          threeVecDauCM = fourVecDauCM.Vect();
          threeVecDauCMXY = ROOT::Math::XYZVector(threeVecDauCM.X(), threeVecDauCM.Y(), 0.);
          beamvector = ROOT::Math::XYZVector(0, 0, 1);
          auto cosThetaStar = beamvector.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(beamvector.Mag2());
          auto SA = cosThetaStar * TMath::Sin(2.0 * phiminuspsi);

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
      if (!collision1.sel8() || !collision1.triggereventep() || !collision1.selection_bit(aod::evsel::kNoTimeFrameBorder) || !collision1.selection_bit(aod::evsel::kNoITSROFrameBorder) || !collision1.selection_bit(aod::evsel::kNoSameBunchPileup) || !collision1.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV) || !collision1.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
        continue;
      }
      if (!collision2.sel8() || !collision2.triggereventep() || !collision2.selection_bit(aod::evsel::kNoTimeFrameBorder) || !collision2.selection_bit(aod::evsel::kNoITSROFrameBorder) || !collision2.selection_bit(aod::evsel::kNoSameBunchPileup) || !collision2.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV) || !collision2.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
        continue;
      }
      if (collision1.bcId() == collision2.bcId()) {
        continue;
      }
      auto centrality = collision1.centFT0C();
      auto psiFT0C = collision1.psiFT0C();
      auto QFT0C = collision1.qFT0C();
      int occupancy1 = collision1.trackOccupancyInTimeRange();
      int occupancy2 = collision1.trackOccupancyInTimeRange();
      for (auto& [track1, v0] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
        if (!selectionTrack(track1)) {
          continue;
        }
        if (!track1.hasTOF()) {
          if (rejectPID && !rejectPi(track1)) {
            continue;
          }
          if (rejectPID && !rejectEl(track1)) {
            continue;
          }
          if (rejectPID && !rejectKa(track1)) {
            continue;
          }
        }

        // PID check
        if (PIDstrategy == 0 && !selectionPID1(track1)) {
          continue;
        }
        if (PIDstrategy == 1 && !selectionPID2(track1)) {
          continue;
        }
        if (PIDstrategy == 2 && !selectionPID3(track1)) {
          continue;
        }
        if (PIDstrategy == 3 && !selectionPID4(track1)) {
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
        Kshort = ROOT::Math::PxPyPzMVector(v0.px(), v0.py(), v0.pz(), v0.mK0Short());
        Lambdac = Proton + Kshort;
        if (Lambdac.Pt() > 6.0 || Lambdac.Pt() < 2.0) {
          continue;
        }
        if (TMath::Abs(Lambdac.Rapidity()) > confRapidity) {
          continue;
        }
        auto phiminuspsi = GetPhiInRange(Lambdac.Phi() - psiFT0C);
        if (useSP) {
          v2 = TMath::Cos(2.0 * phiminuspsi) * QFT0C;
        }
        if (!useSP) {
          v2 = TMath::Cos(2.0 * phiminuspsi);
        }
        auto anglesign = (v0.x() - collision1.posX()) * v0.px() + (v0.y() - collision1.posY()) * v0.py() + (v0.z() - collision1.posZ()) * v0.pz();
        anglesign = anglesign / TMath::Abs(anglesign);
        auto dcasum = 0.0;
        if (useSignDCAV0) {
          dcasum = anglesign * (v0.dcav0topv()) - track1.dcaXY();
        }
        if (!useSignDCAV0) {
          dcasum = v0.dcav0topv() - track1.dcaXY();
        }

        if (occupancy1 < cfgOccupancyCut && occupancy2 < cfgOccupancyCut && Lambdac.M() > cMinLambdaMass && Lambdac.M() <= cMaxLambdaMass && TMath::Abs(Lambdac.Rapidity()) < confRapidity && Lambdac.Pt() > 2.0 && Lambdac.Pt() <= 6.0) {
          if (fillDefault) {
            histos.fill(HIST("hSparseV2SAMixedEvent_V2"), Lambdac.M(), Lambdac.Pt(), v2, TMath::Abs(track1.dcaXY()), Proton.Pt());
          }
          if (fillDecayLength) {
            histos.fill(HIST("hSparseV2SAMixedEvent_V2_dcatopv"), Lambdac.M(), Lambdac.Pt(), v2, v0.dcav0topv(), Proton.Pt());
          }
        }
        if (fillPolarization) {
          ROOT::Math::Boost boost{Lambdac.BoostToCM()};
          fourVecDauCM = boost(Kshort);
          threeVecDauCM = fourVecDauCM.Vect();
          threeVecDauCMXY = ROOT::Math::XYZVector(threeVecDauCM.X(), threeVecDauCM.Y(), 0.);
          beamvector = ROOT::Math::XYZVector(0, 0, 1);
          auto cosThetaStar = beamvector.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(beamvector.Mag2());
          auto SA = cosThetaStar * TMath::Sin(2.0 * phiminuspsi);

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
