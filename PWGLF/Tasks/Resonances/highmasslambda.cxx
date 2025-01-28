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
#include <string>

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
#include "Common/DataModel/PIDResponseITS.h"
// #include "PWGHF/Utils/utilsBfieldCCDB.h"
#include "PWGHF/Utils/utilsEvSelHf.h"
#include "PWGHF/Utils/utilsTrkCandHf.h"
#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/V0.h"
#include "DCAFitter/DCAFitterN.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;
struct highmasslambda {
  int multEstimator;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Service<o2::framework::O2DatabasePDG> pdg;
  o2::base::Propagator::MatCorrType noMatCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;
  o2::vertexing::DCAFitterN<2> df;
  int runNumber{0};
  double bz{0.};

  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> ccdbPathLut{"ccdbPathLut", "GLO/Param/MatLUT", "Path for LUT parametrization"};
  Configurable<std::string> ccdbPathGrp{"ccdbPathGrp", "GLO/GRP/GRP", "Path of the grp file (Run 2)"};
  Configurable<std::string> ccdbPathGrpMag{"ccdbPathGrpMag", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object (Run 3)"};
  Configurable<bool> isRun2{"isRun2", false, "enable Run 2 or Run 3 GRP objects for magnetic field"};
  Configurable<bool> cnfabsdca{"cnfabsdca", false, "Use Abs DCA for secondary vertex fitting"};

  // fill output
  Configurable<int> cfgOccupancyCut{"cfgOccupancyCut", 2500, "Occupancy cut"};
  Configurable<bool> fillRotation{"fillRotation", false, "fill rotation"};
  // events
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCutCentralityMax{"cfgCutCentralityMax", 50.0f, "Accepted maximum Centrality"};
  Configurable<float> cfgCutCentralityMin{"cfgCutCentralityMin", 30.0f, "Accepted minimum Centrality"};
  Configurable<bool> additionalEvSel{"additionalEvSel", true, "additionalEvSel"};
  // proton track cut
  Configurable<float> confMinRot{"confMinRot", 5.0 * TMath::Pi() / 6.0, "Minimum of rotation"};
  Configurable<float> confMaxRot{"confMaxRot", 7.0 * TMath::Pi() / 6.0, "Maximum of rotation"};
  Configurable<float> confRapidity{"confRapidity", 0.8, "cut on Rapidity"};
  Configurable<float> cfgCutPT{"cfgCutPT", 0.4, "PT cut on daughter track"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8, "Eta cut on daughter track"};
  Configurable<float> cfgCutDCAxymin1{"cfgCutDCAxymin1", 0.005f, "Minimum DCAxy range for tracks pt 0 to 0.5"};
  Configurable<float> cfgCutDCAxymin2{"cfgCutDCAxymin2", 0.003f, "Minimum DCAxy range for tracks pt 0.5 to 1"};
  Configurable<float> cfgCutDCAxymin3{"cfgCutDCAxymin3", 0.003f, "Minimum DCAxy range for tracks pt 1.0 to 1.5"};
  Configurable<float> cfgCutDCAxymin4{"cfgCutDCAxymin4", 0.002f, "Minimum DCAxy range for tracks pt 1.5 to 2.0"};
  Configurable<float> cfgCutDCAxymin5{"cfgCutDCAxymin5", 0.001f, "Minimum DCAxy range for tracks pt 2.0 to 1000.5"};
  Configurable<float> cfgCutDCAxy{"cfgCutDCAxy", 0.1f, "DCAxy range for tracks"};
  Configurable<float> cfgCutDCAz{"cfgCutDCAz", 1.0f, "DCAz range for tracks"};
  Configurable<int> cfgITScluster{"cfgITScluster", 5, "Number of ITS cluster"};
  Configurable<int> cfgTPCcluster{"cfgTPCcluster", 70, "Number of TPC cluster"};
  Configurable<int> PIDstrategy{"PIDstrategy", 0, "0: TOF Veto, 1: TOF Veto opti, 2: TOF, 3: TOF loose 1, 4: TOF loose 2, 5: old pt dep"};
  Configurable<float> nsigmaCutTPC{"nsigmacutTPC", 3.0, "Value of the TPC Nsigma cut"};
  Configurable<float> nsigmaCutTOF{"nsigmaCutTOF", 3.0, "TOF PID"};
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
  Configurable<double> ConfDaughDCAMin{"ConfDaughDCAMin", 0.08f, "V0 Daugh sel:  Max. DCA Daugh to PV (cm)"};
  Configurable<float> ConfDaughPIDCuts{"ConfDaughPIDCuts", 3, "PID selections for KS0 daughters"};
  // config SVx
  Configurable<double> ConfMaxDecayLength{"ConfMaxDecayLength", 0.1f, "Maximum decay length (cm)"};
  Configurable<double> ConfMinCPA{"ConfMinCPA", 0.9f, "Minimum CPA"};
  // Mixed event
  Configurable<int> cfgNoMixedEvents{"cfgNoMixedEvents", 1, "Number of mixed events per event"};
  /// activate rotational background
  Configurable<int> nBkgRotations{"nBkgRotations", 9, "Number of rotated copies (background) per each original candidate"};
  // THnsparse bining
  ConfigurableAxis configThnAxisInvMass{"configThnAxisInvMass", {60, 2.15, 2.45}, "#it{M} (GeV/#it{c}^{2})"};
  ConfigurableAxis configThnAxisPt{"configThnAxisPt", {5, 1.0, 6.}, "#it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis configThnAxisV2{"configThnAxisV2", {80, -1, 1}, "V2"};
  ConfigurableAxis cnfigThnAxisDCA{"cnfigThnAxisDCA", {100, 0.0, 0.1}, "DCA"};
  ConfigurableAxis cnfigThnAxisDecayLength{"cnfigThnAxisDecayLength", {150, 0.0, 0.3}, "decay length"};
  ConfigurableAxis cnfigThnAxisPtProton{"cnfigThnAxisPtProton", {16, 0.0, 8.0}, "pT"};
  ConfigurableAxis cnfigThnAxisCPA{"cnfigThnAxisCPA", {300, 0.8, 1.1}, "CPA"};
  // ConfigurableAxis configThnAxisCosThetaStar{"configThnAxisCosThetaStar", {10, -1.0, 1.}, "cos(#vartheta)"};
  // ConfigurableAxis configThnAxisSA{"configThnAxisSA", {100, -1, 1}, "SA"};

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter centralityFilter = (nabs(aod::cent::centFT0C) < cfgCutCentralityMax && nabs(aod::cent::centFT0C) > cfgCutCentralityMin);
  Filter acceptanceFilter = (nabs(aod::track::eta) < cfgCutEta && nabs(aod::track::pt) > cfgCutPT);
  Filter dcaCutFilter = (nabs(aod::track::dcaXY) < cfgCutDCAxy) && (nabs(aod::track::dcaZ) < cfgCutDCAz);
  Filter pidFilter = nabs(aod::pidtpc::tpcNSigmaPr) < nsigmaCutTPC;

  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::FV0Mults, aod::TPCMults, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0As, aod::EPCalibrationTables, aod::Mults>>;

  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPr, aod::pidTOFFullPr>>;
  using AllTrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPi>;
  using ResoV0s = aod::V0Datas;

  using TrackCandidatesSvx = soa::Filtered<soa::Join<aod::TracksIU, aod::TracksCovIU, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPr, aod::pidTOFFullPr>>;
  using AllTrackCandidatesSvx = soa::Join<aod::TracksIU, aod::TracksCovIU, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPi>;
  using ResoV0sSvx = soa::Join<aod::V0Datas, aod::V0Covs, aod::V0DauCovs>;

  SliceCache cache;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext&)
  {

    // std::vector<double> dcaBinning = {0.0, 0.0005, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01, 0.012, 0.014, 0.016, 0.02, 0.03, 0.05, 0.1, 0.5, 1.0};
    // std::vector<double> dcaBinning = {0.0, 0.0005, 0.001, 0.0012, 0.0014, 0.0016, 0.002, 0.0025, 0.003, 0.004, 0.005, 0.006, 0.008, 0.01, 0.015, 0.02, 0.04, 0.05, 0.06, 0.08, 0.1, 0.3, 1.0};
    // std::vector<double> ptProtonBinning = {0.2, 0.3, 0.5, 0.6, 0.8, 1.2, 1.4, 1.6, 2.0, 3.0, 4.0, 6.0};
    // std::vector<double> ptLambdaBinning = {2.0, 3.0, 4.0, 5.0, 6.0};

    std::vector<double> occupancyBinning = {-0.5, 500.0, 1000.0, 1500.0, 2000.0, 3000.0, 4000.0, 5000.0, 50000.0};
    AxisSpec resAxis = {1600, -16, 16, "Res"};
    AxisSpec phiAxis = {500, -6.28, 6.28, "phi"};
    AxisSpec centAxis = {8, 0, 80, "V0M (%)"};
    const AxisSpec thnAxisInvMass{configThnAxisInvMass, "#it{M} (GeV/#it{c}^{2})"};
    const AxisSpec thnAxisPt{configThnAxisPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec thnAxisV2{configThnAxisV2, "V2"};
    const AxisSpec thnAxisDCA{cnfigThnAxisDCA, "DCAxy"};
    const AxisSpec thnAxisDecayLength{cnfigThnAxisDecayLength, "Decay Length"};
    const AxisSpec thnAxisPtProton{cnfigThnAxisPtProton, "Proton Pt"};
    const AxisSpec thnAxisCPA{cnfigThnAxisCPA, "CPA"};
    AxisSpec occupancyAxis = {occupancyBinning, "occupancy"};

    // const AxisSpec thnAxisCosThetaStar{configThnAxisCosThetaStar, "cos(#vartheta)"};
    // const AxisSpec thnAxisPhiminusPsi{configThnAxisPhiminusPsi, "#phi - #psi"};
    // const AxisSpec thnAxisCentrality{configThnAxisCentrality, "Centrality (%)"};
    // const AxisSpec thnAxisSA{configThnAxisSA, "SA"};

    histos.add("hMomCorr", "hMomCorr", kTH2F, {{200, -10.0f, 10.0f}, {200, -10.0f, 10.0f}});
    histos.add("hInvMassKs0", "hInvMassKs0", kTH1F, {{200, 0.4f, 0.6f}});
    histos.add("hV0Dca", "hV0Dca", kTH1F, {{2000, -1.0f, 1.0f}});
    histos.add("hpTvsRapidity", "pT vs Rapidity", kTH2F, {{100, 0.0f, 10.0f}, {300, -1.5f, 1.5f}});
    histos.add("hCentrality", "Centrality distribution", kTH1F, {{200, 0.0, 200.0}});
    histos.add("hVtxZ", "Vertex distribution in Z;Z (cm)", kTH1F, {{400, -20.0, 20.0}});
    histos.add("hOccupancy", "Occupancy", kTH1F, {{5000, -0.5, 50000.5}});
    histos.add("hRotation", "hRotation", kTH1F, {{360, 0.0, 2.0 * TMath::Pi()}});
    histos.add("hEta", "Eta distribution", kTH1F, {{200, -1.0f, 1.0f}});
    histos.add("hDcaxy", "Dcaxy distribution", kTH1F, {{1000, -0.5f, 0.5f}});
    histos.add("hDcaz", "Dcaz distribution", kTH1F, {{1000, -0.5f, 0.5f}});
    histos.add("hNsigmaProtonTPC", "NsigmaProton TPC distribution", kTH2F, {{100, -5.0f, 5.0f}, {60, 0.0f, 6.0f}});
    histos.add("hNsigmaProtonTOF", "NsigmaProton TOF distribution", kTH2F, {{1000, -50.0f, 50.0f}, {60, 0.0f, 6.0f}});
    histos.add("hNsigmaProtonTPCPre", "NsigmaProton TPC distribution Pre sel", kTH2F, {{1000, -50.0f, 50.0f}, {60, 0.0f, 6.0f}});
    histos.add("hPsiFT0C", "PsiFT0C", kTH3F, {centAxis, phiAxis, occupancyAxis});
    histos.add("hPsiFT0A", "PsiFT0A", kTH3F, {centAxis, phiAxis, occupancyAxis});
    histos.add("hPsiTPC", "PsiTPC", kTH3F, {centAxis, phiAxis, occupancyAxis});

    // SVX histo
    histos.add("hDecayLengthxy", "Decay length xy", kTH1F, {{500, 0.0f, 0.1f}});
    histos.add("hDecayLength", "Decay length", kTH1F, {{500, 0.0f, 0.1f}});
    histos.add("hImpactPar0", "hImpactPar0", kTH1F, {{500, 0.0f, 0.1f}});
    histos.add("hImpactPar1", "hImpactPar1", kTH1F, {{500, 0.0f, 0.1f}});
    histos.add("hCPA", "hCPA", kTH1F, {{220, -1.1f, 1.1f}});

    histos.add("hSparseV2SASameEvent_V2", "hSparseV2SASameEvent_V2", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisV2, thnAxisDCA});
    histos.add("hSparseV2SASameEventRotational_V2", "hSparseV2SASameEventRotational", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisV2, thnAxisDCA});
    histos.add("hSparseV2SAMixedEvent_V2", "hSparseV2SAMixedEvent_V2", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisV2, thnAxisDCA});

    histos.add("hSparseV2SASameEvent_V2_SVX", "hSparseV2SASameEvent_V2_SVX", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisV2, thnAxisDecayLength, thnAxisCPA});
    histos.add("hSparseV2SASameEventRotational_V2_SVX", "hSparseV2SASameEventRotational_SVX", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisV2, thnAxisDecayLength, thnAxisCPA});

    // histogram for resolution
    histos.add("ResFT0CTPC", "ResFT0CTPC", kTH2F, {centAxis, resAxis});
    histos.add("ResFT0CTPCR", "ResFT0CTPCR", kTH2F, {centAxis, resAxis});
    histos.add("ResFT0CTPCL", "ResFT0CTPCL", kTH2F, {centAxis, resAxis});
    histos.add("ResTPCRTPCL", "ResTPCRTPCL", kTH2F, {centAxis, resAxis});
    histos.add("ResFT0CFT0A", "ResFT0CFT0A", kTH2F, {centAxis, resAxis});
    histos.add("ResFT0ATPC", "ResFT0ATPC", kTH2F, {centAxis, resAxis});

    df.setPropagateToPCA(true);
    df.setMaxR(200);
    df.setMaxDZIni(4);
    df.setMinParamChange(1.e-3);
    df.setMinRelChi2Change(0.9);
    df.setUseAbsDCA(cnfabsdca);
    df.setWeightedFinalPCA(cnfabsdca);
    df.setMatCorrType(noMatCorr);

    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    runNumber = 0;
    bz = 0;
  }
  template <typename T>
  bool selectionTrack(const T& candidate)
  {
    if (!(candidate.isGlobalTrackWoDCA() && candidate.itsNCls() > cfgITScluster && candidate.tpcNClsFound() > cfgTPCcluster)) {
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
    if (candidate.pt() >= 2.0 && candidate.pt() < 10000000.0 && TMath::Abs(candidate.dcaXY()) < cfgCutDCAxymin5) {
      return false;
    }
    return true;
  }

  // TOF Veto
  template <typename T>
  bool selectionPID1(const T& candidate)
  {
    if (candidate.hasTOF()) {
      if (candidate.p() < 0.5 && TMath::Abs(candidate.tpcNSigmaPr()) < nsigmaCutTPC) {
        return true;
      }
      if (candidate.p() >= 0.5 && TMath::Abs(candidate.tpcNSigmaPr()) < nsigmaCutTPC && TMath::Abs(candidate.tofNSigmaPr()) < nsigmaCutTOF) {
        return true;
      }
    }
    if (!candidate.hasTOF()) {
      if (TMath::Abs(candidate.tpcNSigmaPr()) < nsigmaCutTPC) {
        return true;
      }
    }
    return false;
  }

  // TPC TOF
  template <typename T>
  bool selectionPID2(const T& candidate)
  {
    if (candidate.tpcInnerParam() < 1.0 && TMath::Abs(candidate.tpcNSigmaPr()) < nsigmaCutTPC) {
      return true;
    }
    if (candidate.tpcInnerParam() >= 1.0 && candidate.tpcInnerParam() < 2.0 && TMath::Abs(candidate.tpcNSigmaPr()) < nsigmaCutTPC && candidate.hasTOF() && TMath::Abs(candidate.tofNSigmaPr()) < nsigmaCutTOF) {
      return true;
    }
    if (candidate.tpcInnerParam() >= 2.0) {
      if (candidate.hasTOF() && TMath::Abs(candidate.tpcNSigmaPr()) < nsigmaCutTPC && TMath::Abs(candidate.tofNSigmaPr()) < nsigmaCutTOF) {
        return true;
      }
      if (!candidate.hasTOF() && TMath::Abs(candidate.tpcNSigmaPr()) < nsigmaCutTPC) {
        return true;
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
    if (TMath::Abs(eta) > 0.8) {
      return false;
    }
    if (TMath::Abs(pt) < 0.15) {
      return false;
    }
    if (tpcNClsF < 50) {
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
  double massPr = TDatabasePDG::Instance()->GetParticle(kProton)->Mass();   // FIXME: Get from the common header
  double massK0s = TDatabasePDG::Instance()->GetParticle(kK0Short)->Mass(); // FIXME: Get from the common header
  double v2, v2Rot;
  void processSameEvent(EventCandidates::iterator const& collision, TrackCandidates const& tracks, AllTrackCandidates const&, ResoV0s const& V0s, aod::BCs const&)
  {
    if (!collision.sel8() || !collision.triggereventep() || !collision.selection_bit(aod::evsel::kNoTimeFrameBorder) || !collision.selection_bit(aod::evsel::kNoITSROFrameBorder) || !collision.selection_bit(aod::evsel::kNoSameBunchPileup) || !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV) || !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      return;
    }
    if (additionalEvSel && !collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
      return;
    }
    o2::aod::ITSResponse itsResponse;
    auto centrality = collision.centFT0C();
    // auto multTPC = collision.multNTracksPV();
    int occupancy = collision.trackOccupancyInTimeRange();
    if (occupancy > cfgOccupancyCut) {
      return;
    }
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
    histos.fill(HIST("hPsiFT0C"), centrality, psiFT0C, occupancy);
    histos.fill(HIST("hPsiFT0A"), centrality, psiFT0A, occupancy);
    histos.fill(HIST("hPsiTPC"), centrality, psiTPC, occupancy);
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
      histos.fill(HIST("hNsigmaProtonTPCPre"), track1.tpcNSigmaPr(), track1.pt());
      // PID check
      if (PIDstrategy == 0 && !selectionPID1(track1)) {
        continue;
      }
      if (PIDstrategy == 1 && !selectionPID2(track1)) {
        continue;
      }
      if (track1.p() < 1.0 && !(itsResponse.nSigmaITS<o2::track::PID::Proton>(track1) > -2.5 && itsResponse.nSigmaITS<o2::track::PID::Proton>(track1) < 2.5)) {
        continue;
      }
      histos.fill(HIST("hEta"), track1.eta());
      histos.fill(HIST("hDcaxy"), track1.dcaXY());
      histos.fill(HIST("hDcaz"), track1.dcaZ());
      histos.fill(HIST("hNsigmaProtonTPC"), track1.tpcNSigmaPr(), track1.pt());
      if (track1.hasTOF()) {
        histos.fill(HIST("hNsigmaProtonTOF"), track1.tofNSigmaPr(), track1.pt());
      }
      auto track1ID = track1.globalIndex();
      for (auto v0 : V0s) {
        if (!SelectionV0(collision, v0)) {
          continue;
        }
        if (firstprimarytrack == 0) {
          histos.fill(HIST("hV0Dca"), v0.dcav0topv());
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
        v2 = TMath::Cos(2.0 * phiminuspsi) * QFT0C;
        if (Lambdac.M() > cMinLambdaMass && Lambdac.M() <= cMaxLambdaMass && TMath::Abs(Lambdac.Rapidity()) < confRapidity && Lambdac.Pt() > 2.0 && Lambdac.Pt() <= 6.0) {
          histos.fill(HIST("hSparseV2SASameEvent_V2"), Lambdac.M(), Lambdac.Pt(), v2, TMath::Abs(track1.dcaXY()));
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
            v2Rot = TMath::Cos(2.0 * phiminuspsiRot) * QFT0C;
            if (LambdacRot.M() > cMinLambdaMass && LambdacRot.M() <= cMaxLambdaMass && TMath::Abs(LambdacRot.Rapidity()) < confRapidity && LambdacRot.Pt() > 2.0 && LambdacRot.Pt() <= 6.0) {
              histos.fill(HIST("hSparseV2SASameEventRotational_V2"), LambdacRot.M(), LambdacRot.Pt(), v2Rot, TMath::Abs(track1.dcaXY()));
            }
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
      if (additionalEvSel && !collision1.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
        continue;
      }
      if (additionalEvSel && !collision2.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
        continue;
      }
      o2::aod::ITSResponse itsResponse;

      auto psiFT0C = collision1.psiFT0C();
      auto QFT0C = collision1.qFT0C();
      int occupancy1 = collision1.trackOccupancyInTimeRange();
      int occupancy2 = collision1.trackOccupancyInTimeRange();
      if (occupancy1 > cfgOccupancyCut) {
        continue;
      }
      if (occupancy2 > cfgOccupancyCut) {
        continue;
      }
      for (auto& [track1, v0] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {

        if (!selectionTrack(track1)) {
          continue;
        }
        // PID check
        if (PIDstrategy == 0 && !selectionPID1(track1)) {
          continue;
        }
        if (PIDstrategy == 1 && !selectionPID2(track1)) {
          continue;
        }
        if (track1.p() < 1.0 && !(itsResponse.nSigmaITS<o2::track::PID::Proton>(track1) > -2.5 && itsResponse.nSigmaITS<o2::track::PID::Proton>(track1) < 2.5)) {
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
        v2 = TMath::Cos(2.0 * phiminuspsi) * QFT0C;
        if (occupancy1 < cfgOccupancyCut && occupancy2 < cfgOccupancyCut && Lambdac.M() > cMinLambdaMass && Lambdac.M() <= cMaxLambdaMass && TMath::Abs(Lambdac.Rapidity()) < confRapidity && Lambdac.Pt() > 2.0 && Lambdac.Pt() <= 6.0) {
          histos.fill(HIST("hSparseV2SAMixedEvent_V2"), Lambdac.M(), Lambdac.Pt(), v2, TMath::Abs(track1.dcaXY()));
        }
      }
    }
  }
  PROCESS_SWITCH(highmasslambda, processMixedEventOpti, "Process Mixed event new", false);
  void processSameEventSvx(EventCandidates::iterator const& collision, TrackCandidatesSvx const& tracks, AllTrackCandidatesSvx const&, ResoV0sSvx const& V0s, aod::BCsWithTimestamps const&)
  {
    if (!collision.sel8() || !collision.triggereventep() || !collision.selection_bit(aod::evsel::kNoTimeFrameBorder) || !collision.selection_bit(aod::evsel::kNoITSROFrameBorder) || !collision.selection_bit(aod::evsel::kNoSameBunchPileup) || !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV) || !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      return;
    }
    if (additionalEvSel && !collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
      return;
    }

    /// Set the magnetic field from ccdb.
    auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
    if (runNumber != bc.runNumber()) {
      o2::parameters::GRPMagField* grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(ccdbPathGrpMag, bc.timestamp());
      if (grpo == nullptr) {
        LOGF(fatal, "Run 3 GRP object (type o2::parameters::GRPMagField) is not available in CCDB for run=%d at timestamp=%llu", bc.runNumber(), bc.timestamp());
      }
      o2::base::Propagator::initFieldFromGRP(grpo);
      bz = o2::base::Propagator::Instance()->getNominalBz();
      runNumber = bc.runNumber();
    }
    df.setBz(bz);
    int occupancy = collision.trackOccupancyInTimeRange();
    if (occupancy > cfgOccupancyCut) {
      return;
    }
    o2::aod::ITSResponse itsResponse;
    auto centrality = collision.centFT0C();
    // auto multTPC = collision.multNTracksPV();
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
    histos.fill(HIST("hPsiFT0C"), centrality, psiFT0C, occupancy);
    histos.fill(HIST("hPsiFT0A"), centrality, psiFT0A, occupancy);
    histos.fill(HIST("hPsiTPC"), centrality, psiTPC, occupancy);
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
      histos.fill(HIST("hNsigmaProtonTPCPre"), track1.tpcNSigmaPr(), track1.pt());
      // PID check
      if (PIDstrategy == 0 && !selectionPID1(track1)) {
        continue;
      }
      if (PIDstrategy == 1 && !selectionPID2(track1)) {
        continue;
      }
      if (track1.p() < 1.0 && !(itsResponse.nSigmaITS<o2::track::PID::Proton>(track1) > -2.5 && itsResponse.nSigmaITS<o2::track::PID::Proton>(track1) < 2.5)) {
        continue;
      }
      histos.fill(HIST("hEta"), track1.eta());
      histos.fill(HIST("hDcaxy"), track1.dcaXY());
      histos.fill(HIST("hDcaz"), track1.dcaZ());
      histos.fill(HIST("hNsigmaProtonTPC"), track1.tpcNSigmaPr(), track1.pt());
      if (track1.hasTOF()) {
        histos.fill(HIST("hNsigmaProtonTOF"), track1.tofNSigmaPr(), track1.pt());
      }
      auto track1ID = track1.globalIndex();
      auto trackParCovBach = getTrackParCov(track1);

      for (auto v0 : V0s) {
        if (!SelectionV0(collision, v0)) {
          continue;
        }
        if (firstprimarytrack == 0) {
          histos.fill(HIST("hV0Dca"), v0.dcav0topv());
        }
        auto postrack = v0.template posTrack_as<AllTrackCandidatesSvx>();
        auto negtrack = v0.template negTrack_as<AllTrackCandidatesSvx>();
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
        float v0x, v0y, v0z, v0px, v0py, v0pz;
        float posTrackX, negTrackX;
        o2::track::TrackParCov trackParCovV0DaughPos;
        o2::track::TrackParCov trackParCovV0DaughNeg;
        trackParCovV0DaughPos = getTrackParCov(postrack); // check that aod::TracksWCov does not need TracksDCA!
        trackParCovV0DaughNeg = getTrackParCov(negtrack); // check that aod::TracksWCov does not need TracksDCA!
        posTrackX = v0.posX();
        negTrackX = v0.negX();
        v0x = v0.x();
        v0y = v0.y();
        v0z = v0.z();
        const std::array<float, 3> vertexV0 = {v0x, v0y, v0z};

        v0px = v0.px();
        v0py = v0.py();
        v0pz = v0.pz();
        const std::array<float, 3> momentumV0 = {v0px, v0py, v0pz};

        std::array<float, 6> covV0Pos = {0.};
        for (int i = 0; i < 6; i++) {
          covV0Pos[i] = v0.positionCovMat()[i];
        }
        trackParCovV0DaughPos.propagateTo(posTrackX, bz); // propagate the track to the X closest to the V0 vertex
        trackParCovV0DaughNeg.propagateTo(negTrackX, bz); // propagate the track to the X closest to the V0 vertex

        // we build the neutral track to then build the cascade
        // auto trackV0 = o2::dataformats::V0(vertexV0, momentumV0, {0, 0, 0, 0, 0, 0}, trackParCovV0DaughPos, trackParCovV0DaughNeg); // build the V0 track (indices for v0 daughters set to 0 for now)
        auto trackV0 = o2::dataformats::V0(vertexV0, momentumV0, covV0Pos, trackParCovV0DaughPos, trackParCovV0DaughNeg); // build the V0 track (indices for v0 daughters set to 0 for now)
        std::array<float, 3> pVecV0 = {0., 0., 0.};
        std::array<float, 3> pVecBach = {0., 0., 0.};
        std::array<float, 3> pVecCand = {0., 0., 0.};
        try {
          if (df.process(trackV0, trackParCovBach) == 0) {
            continue;
          } else {
          }
        } catch (const std::runtime_error& error) {
          continue;
        }
        df.propagateTracksToVertex();           // propagate the bachelor and V0 to the Lambdac vertex
        trackV0.getPxPyPzGlo(pVecV0);           // momentum of D0 at the Lambdac vertex
        trackParCovBach.getPxPyPzGlo(pVecBach); // momentum of proton at the Lambdac vertex
        pVecCand = RecoDecay::pVec(pVecV0, pVecBach);
        const auto& secondaryVertex = df.getPCACandidate();

        // get track impact parameters
        // This modifies track momenta!
        auto primaryVertex = getPrimaryVertex(collision);
        o2::dataformats::DCA impactParameter0;
        o2::dataformats::DCA impactParameter1;
        trackV0.propagateToDCA(primaryVertex, bz, &impactParameter0);
        trackParCovBach.propagateToDCA(primaryVertex, bz, &impactParameter1);
        double phi, theta;
        getPointDirection(std::array{collision.posX(), collision.posY(), collision.posZ()}, secondaryVertex, phi, theta);
        // auto errorDecayLength = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, theta) + getRotatedCovMatrixXX(covMatrixPCA, phi, theta));
        // auto errorDecayLengthXY = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, 0.) + getRotatedCovMatrixXX(covMatrixPCA, phi, 0.));

        Kshort = ROOT::Math::PxPyPzMVector(pVecV0[0], pVecV0[1], pVecV0[2], massK0s);
        Proton = ROOT::Math::PxPyPzMVector(pVecBach[0], pVecBach[1], pVecBach[2], massPr);
        Lambdac = Proton + Kshort;
        auto phiminuspsi = GetPhiInRange(Lambdac.Phi() - psiFT0C);
        v2 = TMath::Cos(2.0 * phiminuspsi) * QFT0C;

        double protonimpactparameter = impactParameter1.getY();
        double kshortimpactparameter = impactParameter0.getY();

        double decaylengthx = secondaryVertex[0] - collision.posX();
        double decaylengthy = secondaryVertex[1] - collision.posY();
        double decaylengthz = secondaryVertex[2] - collision.posZ();
        double decaylength = TMath::Sqrt(decaylengthx * decaylengthx + decaylengthy * decaylengthy + decaylengthz * decaylengthz);
        double decaylengthxy = TMath::Sqrt(decaylengthx * decaylengthx + decaylengthy * decaylengthy);
        double anglesign = decaylengthx * Lambdac.Px() + decaylengthy * Lambdac.Py() + decaylengthz * Lambdac.Pz();
        double CPAlambdac = anglesign / (decaylength * Lambdac.P());

        histos.fill(HIST("hDecayLengthxy"), decaylengthxy);
        histos.fill(HIST("hDecayLength"), decaylength);
        histos.fill(HIST("hImpactPar0"), protonimpactparameter);
        histos.fill(HIST("hImpactPar1"), kshortimpactparameter);
        histos.fill(HIST("hCPA"), CPAlambdac);
        histos.fill(HIST("hMomCorr"), Proton.P() - track1.p(), v0.p() - Kshort.P());
        if (decaylength > ConfMaxDecayLength) {
          continue;
        }
        if (CPAlambdac < ConfMinCPA) {
          continue;
        }

        if (Lambdac.M() > cMinLambdaMass && Lambdac.M() <= cMaxLambdaMass && TMath::Abs(Lambdac.Rapidity()) < confRapidity && Lambdac.Pt() > 2.0 && Lambdac.Pt() <= 6.0) {
          histos.fill(HIST("hSparseV2SASameEvent_V2_SVX"), Lambdac.M(), Lambdac.Pt(), v2, decaylengthxy, CPAlambdac);
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
            v2Rot = TMath::Cos(2.0 * phiminuspsiRot) * QFT0C;
            if (LambdacRot.M() > cMinLambdaMass && LambdacRot.M() <= cMaxLambdaMass && TMath::Abs(LambdacRot.Rapidity()) < confRapidity && LambdacRot.Pt() > 2.0 && LambdacRot.Pt() <= 6.0) {
              histos.fill(HIST("hSparseV2SASameEventRotational_V2_SVX"), LambdacRot.M(), LambdacRot.Pt(), v2Rot, decaylengthxy, CPAlambdac);
            }
          }
        }
      }
    }
  }
  PROCESS_SWITCH(highmasslambda, processSameEventSvx, "Process Same event SVX", true);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<highmasslambda>(cfgc, TaskName{"highmasslambda"})};
}
