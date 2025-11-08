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

#include "PWGHF/Utils/utilsBfieldCCDB.h"
#include "PWGHF/Utils/utilsEvSelHf.h"
#include "PWGHF/Utils/utilsTrkCandHf.h"
#include "PWGLF/DataModel/EPCalibrationTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/PhysicsConstants.h"
#include "DCAFitter/DCAFitterN.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/StepTHn.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/Track.h"
#include "ReconstructionDataFormats/V0.h"

#include "Math/GenVector/Boost.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "TF1.h"
#include "TRandom3.h"
#include <TDatabasePDG.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <THn.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TPDGCode.h>

#include <array>
#include <cmath>
#include <cstdlib>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;
struct highmasslambdasvx {

  int mRunNumber;
  int multEstimator;
  float d_bz;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Service<o2::framework::O2DatabasePDG> pdg;

  o2::vertexing::DCAFitterN<2> df; // 2-prong vertex fitter
  o2::base::MatLayerCylSet* lut;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
  int runNumber{0};
  double bz = 0.;

  // CCDB options
  // Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  // Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  // Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  // Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  // Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};

  // fill output
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> ccdbPathLut{"ccdbPathLut", "GLO/Param/MatLUT", "Path for LUT parametrization"};
  Configurable<std::string> ccdbPathGrp{"ccdbPathGrp", "GLO/GRP/GRP", "Path of the grp file (Run 2)"};
  Configurable<std::string> ccdbPathGrpMag{"ccdbPathGrpMag", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object (Run 3)"};
  Configurable<bool> isRun2{"isRun2", false, "enable Run 2 or Run 3 GRP objects for magnetic field"};
  Configurable<bool> useSP{"useSP", false, "useSP"};
  Configurable<bool> additionalEvSel{"additionalEvSel", true, "additionalEvSel"};
  Configurable<bool> additionalEvSel2{"additionalEvSel2", false, "additionalEvSel2"};
  // events
  Configurable<bool> cnfabsdca{"cnfabsdca", false, "Use Abs DCA for secondary vertex fitting"};
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCutCentralityMax{"cfgCutCentralityMax", 50.0f, "Accepted maximum Centrality"};
  Configurable<float> cfgCutCentralityMin{"cfgCutCentralityMin", 30.0f, "Accepted minimum Centrality"};
  // proton track cut
  Configurable<bool> useDecayLengthxy{"useDecayLengthxy", true, "use decay length xy"};
  Configurable<bool> ispTdifferentialDCA{"ispTdifferentialDCA", true, "is pT differential DCA"};
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
  Configurable<int> cfgITSclusterInnerlayer{"cfgITSclusterInnerlayer", 1, "Minimum Number of ITS cluster in inner barrel"};
  Configurable<int> cfgTPCcluster{"cfgTPCcluster", 70, "Number of TPC cluster"};
  Configurable<int> PIDstrategy{"PIDstrategy", 0, "0: TOF Veto, 1: TOF Veto opti, 2: TOF, 3: TOF loose 1, 4: TOF loose 2, 5: old pt dep"};
  Configurable<float> nsigmaCutTPC{"nsigmacutTPC", 3.0, "Value of the TPC Nsigma cut"};
  Configurable<float> nsigmaCutCombined{"nsigmaCutCombined", 3.0, "TPC TOF combined PID"};
  Configurable<float> nsigmaCutTPCPre{"nsigmacutTPCPre", 3.0, "Value of the TPC Nsigma cut Pre filter"};
  // Configs for V0
  Configurable<double> ConfV0DCADaughMax{"ConfV0DCADaughMax", 0.2f, "Maximum DCA between the V0 daughters"};
  Configurable<double> ConfV0CPAMin{"ConfV0CPAMin", 0.9998f, "Minimum CPA of V0"};
  Configurable<double> cMaxV0DCA{"cMaxV0DCA", 0.1, "Maximum V0 DCA to PV"};
  Configurable<float> cMaxV0LifeTime{"cMaxV0LifeTime", 20, "Maximum V0 life time"};
  Configurable<float> cSigmaMassKs0{"cSigmaMassKs0", 0.006, "Sigma cut on KS0 mass"};
  Configurable<float> cMinLambdaMass{"cMinLambdaMass", 2.18, "Minimum lambda mass"};
  Configurable<float> cMaxLambdaMass{"cMaxLambdaMass", 2.42, "Maximum lambda mass"};

  // config for V0 daughters
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
  Configurable<double> cutchi2PCA{"cutchi2PCA", 0.1f, "cut on chi2PCA"};
  // THnsparse bining
  ConfigurableAxis configThnAxisInvMass{"configThnAxisInvMass", {60, 2.15, 2.45}, "#it{M} (GeV/#it{c}^{2})"};
  ConfigurableAxis configThnAxisV2{"configThnAxisV2", {80, -1, 1}, "V2"};
  ConfigurableAxis configThnAxisSA{"configThnAxisSA", {100, -1, 1}, "SA"};
  ConfigurableAxis configThnAxisPhiminusPsi{"configThnAxisPhiminusPsi", {6, 0.0, TMath::Pi()}, "#phi - #psi"};
  ConfigurableAxis configThnAxisCentrality{"configThnAxisCentrality", {1, 30., 50}, "Centrality"};
  ConfigurableAxis configThnAxisDecayLength{"configThnAxisDecayLength", {60, 0.0, 0.06}, "Decay length"};

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter centralityFilter = (nabs(aod::cent::centFT0C) < cfgCutCentralityMax && nabs(aod::cent::centFT0C) > cfgCutCentralityMin);
  Filter acceptanceFilter = (nabs(aod::track::eta) < cfgCutEta && nabs(aod::track::pt) > cfgCutPT);
  Filter dcaCutFilter = (nabs(aod::track::dcaXY) < cfgCutDCAxy) && (nabs(aod::track::dcaZ) < cfgCutDCAz);
  Filter pidFilter = nabs(aod::pidtpc::tpcNSigmaPr) < nsigmaCutTPCPre;

  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::FV0Mults, aod::TPCMults, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0As, aod::EPCalibrationTables, aod::Mults>>;
  // using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksIU, aod::TracksCovIU, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPr, aod::pidTOFFullPr, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa>>;
  // using AllTrackCandidates = soa::Join<aod::Tracks, aod::TracksIU, aod::TracksCovIU, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPi>;

  using TrackCandidates = soa::Filtered<soa::Join<aod::TracksIU, aod::TracksCovIU, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPr, aod::pidTOFFullPr, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullEl, aod::pidTOFFullEl>>;
  using AllTrackCandidates = soa::Join<aod::TracksIU, aod::TracksCovIU, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPi>;
  using ResoV0s = soa::Join<aod::V0Datas, aod::V0Covs, aod::V0DauCovs>;

  SliceCache cache;
  // Partition<TrackCandidates> posTracks = aod::track::signed1Pt > cfgCutCharge;
  // Partition<TrackCandidates> negTracks = aod::track::signed1Pt < cfgCutCharge;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext&)
  {
    std::vector<double> occupancyBinning = {0.0, 500.0, 1000.0, 1500.0, 2000.0, 3000.0, 4000.0, 5000.0, 50000.0};
    std::vector<double> dcaBinning = {0.0, 0.0005, 0.001, 0.002, 0.003, 0.004, 0.006, 0.3};
    std::vector<double> ptProtonBinning = {0.0, 0.3, 0.5, 0.8, 1.2, 6.0};
    std::vector<double> ptLambdaBinning = {2.0, 3.0, 4.0, 5.0, 6.0};

    AxisSpec phiAxis = {500, -6.28, 6.28, "phi"};
    AxisSpec decaylengthAxis = {configThnAxisDecayLength, "decaylength"};
    AxisSpec resAxis = {1000, -10, 10, "Res"};
    AxisSpec centAxis = {8, 0, 80, "V0M (%)"};
    const AxisSpec thnAxisInvMass{configThnAxisInvMass, "#it{M} (GeV/#it{c}^{2})"};
    const AxisSpec thnAxisV2{configThnAxisV2, "V2"};
    const AxisSpec thnAxisCentrality{configThnAxisCentrality, "Centrality (%)"};
    const AxisSpec thnAxisPhiminusPsi{configThnAxisPhiminusPsi, "#phi - #psi"};
    AxisSpec occupancyAxis = {occupancyBinning, "occupancy"};
    AxisSpec ptAxis = {ptLambdaBinning, "pt"};
    AxisSpec dcaAxis = {dcaBinning, "dca"};
    AxisSpec ptProtonAxis = {ptProtonBinning, "daughter pt"};
    histos.add("hSparseV2SASameEvent_V2_EP", "hSparseV2SASameEvent_V2_EP", HistType::kTHnSparseF, {thnAxisInvMass, ptAxis, thnAxisV2, ptProtonAxis, decaylengthAxis, dcaAxis});
    histos.add("hSparseV2SASameEventRotational_V2_EP", "hSparseV2SASameEventRotational_V2_EP", HistType::kTHnSparseF, {thnAxisInvMass, ptAxis, thnAxisV2, ptProtonAxis, decaylengthAxis, dcaAxis});
    histos.add("hV0decaylength", "hV0decaylength", kTH1F, {{1000, 0.0f, 1000.0f}});
    histos.add("hMomCorr", "hMomCorr", kTH3F, {{200, -10.0f, 10.0f}, {200, -10.0f, 10.0f}, {8, 0.0f, 80.0f}});
    histos.add("hInvMassKs0", "hInvMassKs0", kTH1F, {{200, 0.4f, 0.6f}});
    histos.add("hchi2PCA", "hchi2PCA", kTH1F, {{1000, 0.0f, 1.f}});
    histos.add("hFTOCvsTPCNoCut", "Mult correlation FT0C vs. TPC without any cut", kTH2F, {{80, 0.0f, 80.0f}, {100, -0.5f, 5999.5f}});
    histos.add("hFTOCvsTPC", "Mult correlation FT0C vs. TPC", kTH2F, {{80, 0.0f, 80.0f}, {100, -0.5f, 5999.5f}});
    histos.add("hFTOCvsTPCSelected", "Mult correlation FT0C vs. TPC after selection", kTH2F, {{80, 0.0f, 80.0f}, {100, -0.5f, 5999.5f}});
    histos.add("hCentrality", "Centrality distribution", kTH1F, {{200, 0.0, 200.0}});
    histos.add("hVtxZ", "Vertex distribution in Z;Z (cm)", kTH1F, {{400, -20.0, 20.0}});
    histos.add("hOccupancy", "Occupancy", kTH1F, {{5000, 0.0, 50000.0}});
    histos.add("hRotation", "hRotation", kTH1F, {{360, 0.0, 2.0 * TMath::Pi()}});

    histos.add("hNsigmaProtonTPCDiff", "Difference NsigmaProton NsigmaKaon TPC distribution", kTH3F, {{100, -5.0f, 5.0f}, {100, -5.0f, 5.0f}, {80, 0.0f, 8.0f}});
    histos.add("hNsigmaProtonTPC", "NsigmaProton TPC distribution", kTH2F, {{100, -5.0f, 5.0f}, {80, 0.0f, 8.0f}});
    histos.add("hNsigmaProtonTOF", "NsigmaProton TOF distribution", kTH2F, {{100, -5.0f, 5.0f}, {80, 0.0f, 8.0f}});

    histos.add("hEta", "Eta distribution", kTH1F, {{200, -1.0f, 1.0f}});
    histos.add("hDcaz", "Dcaz distribution", kTH1F, {{1000, -0.5f, 0.5f}});
    histos.add("hDcaxy", "Dcaxy distribution", kTH1F, {{1000, -0.5f, 0.5f}});
    histos.add("hsignprotonDca", "hsignprotonDca", kTH1F, {{2000, -1.0f, 1.0f}});
    histos.add("hunsignprotonDca", "hunsignprotonDca", kTH1F, {{2000, -1.0f, 1.0f}});
    histos.add("hsignCPALambdac", "hsignCPALambdac", kTH1F, {{2000, -1.0f, 1.0f}});

    histos.add("hPsiFT0C", "PsiFT0C", kTH3F, {centAxis, phiAxis, occupancyAxis});
    histos.add("hPsiFT0A", "PsiFT0A", kTH3F, {centAxis, phiAxis, occupancyAxis});
    histos.add("hPsiTPC", "PsiTPC", kTH3F, {centAxis, phiAxis, occupancyAxis});
    histos.add("hPsiTPCR", "PsiTPCR", kTH3F, {centAxis, phiAxis, occupancyAxis});
    histos.add("hPsiTPCL", "PsiTPCL", kTH3F, {centAxis, phiAxis, occupancyAxis});

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
    df.setWeightedFinalPCA(true);
    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(ccdbPathLut));
    runNumber = 0;
  }
  template <typename T>
  bool selectionTrack(const T& candidate)
  {
    if (!(candidate.isGlobalTrackWoDCA() && candidate.itsNCls() > cfgITScluster && candidate.tpcNClsFound() > cfgTPCcluster && candidate.itsNClsInnerBarrel() >= cfgITSclusterInnerlayer)) {
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
    if (candidate.pt() >= 2.0 && candidate.pt() < 3.0 && TMath::Abs(candidate.dcaXY()) < cfgCutDCAxymin5) {
      return false;
    }
    if (candidate.pt() >= 3.0 && TMath::Abs(candidate.dcaXY()) < cfgCutDCAxymin6) {
      return false;
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
    if (candidate.tpcInnerParam() < 0.7 && TMath::Abs(candidate.tpcNSigmaPr()) < 3.0) {
      return true;
    }
    if (candidate.tpcInnerParam() >= 0.7) {
      // printf("I am here: %.3f\n", candidate.tpcInnerParam());
      if (candidate.hasTOF()) {
        auto combinedPID = TMath::Sqrt(candidate.tpcNSigmaPr() * candidate.tpcNSigmaPr() + candidate.tofNSigmaPr() * candidate.tofNSigmaPr()) / TMath::Sqrt(2.0);
        // printf("combine PIDA: %.3f\n", combinedPID);
        if (combinedPID < nsigmaCutCombined) {
          return true;
        }
      }
      if (!candidate.hasTOF()) {
        if (candidate.tpcInnerParam() < 1.5 && TMath::Abs(candidate.tpcNSigmaPr()) < 3.0) {
          return true;
        }
        if (candidate.tpcInnerParam() >= 1.5 && candidate.tpcNSigmaPr() > -2.0 && candidate.tpcNSigmaPr() < 2.0) {
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
    if (candidate.tpcInnerParam() < 0.7 && TMath::Abs(candidate.tpcNSigmaPr()) < 3.0) {
      return true;
    }
    if (candidate.tpcInnerParam() >= 0.7) {
      if (candidate.hasTOF()) {
        auto combinedPID = TMath::Sqrt(candidate.tpcNSigmaPr() * candidate.tpcNSigmaPr() + candidate.tofNSigmaPr() * candidate.tofNSigmaPr()) / TMath::Sqrt(2.0);
        if (combinedPID < nsigmaCutCombined) {
          return true;
        }
      }
    }
    return false;
  }

  // TOF veto loose
  template <typename T>
  bool selectionPID3(const T& candidate)
  {
    if (candidate.tpcInnerParam() < 0.7 && TMath::Abs(candidate.tpcNSigmaPr()) < 3.0) {
      return true;
    }
    if (candidate.tpcInnerParam() >= 0.7) {
      if (candidate.hasTOF()) {
        auto combinedPID = TMath::Sqrt(candidate.tpcNSigmaPr() * candidate.tpcNSigmaPr() + candidate.tofNSigmaPr() * candidate.tofNSigmaPr()) / TMath::Sqrt(2.0);
        if (combinedPID < nsigmaCutCombined) {
          return true;
        }
      }
      if (!candidate.hasTOF()) {
        if (candidate.tpcInnerParam() < 1.5 && TMath::Abs(candidate.tpcNSigmaPr()) < 3.0) {
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
    if (candidate.tpcInnerParam() < 0.7 && TMath::Abs(candidate.tpcNSigmaPr()) < 3.0) {
      return true;
    }
    if (candidate.tpcInnerParam() >= 0.7) {
      if (candidate.hasTOF()) {
        auto combinedPID = TMath::Sqrt(candidate.tpcNSigmaPr() * candidate.tpcNSigmaPr() + candidate.tofNSigmaPr() * candidate.tofNSigmaPr()) / TMath::Sqrt(2.0);
        if (combinedPID < nsigmaCutCombined) {
          return true;
        }
      }
      if (!candidate.hasTOF()) {
        if (candidate.tpcInnerParam() < 1.5 && TMath::Abs(candidate.tpcNSigmaPr()) < 3.0) {
          return true;
        }
        if (candidate.tpcInnerParam() >= 1.5 && candidate.tpcInnerParam() < 1.8 && candidate.tpcNSigmaPr() > -1.5 && candidate.tpcNSigmaPr() < 2.0) {
          return true;
        }
      }
    }
    return false;
  }

  template <typename Collision, typename V0>
  bool SelectionV0(Collision const& collision, V0 const& candidate)
  {
    // const float pT = candidate.pt();
    const std::vector<float> decVtx = {candidate.x(), candidate.y(), candidate.z()};
    const float tranRad = candidate.v0radius();
    const double dcaDaughv0 = TMath::Abs(candidate.dcaV0daughters());
    const double cpav0 = candidate.v0cosPA();
    float decaylength = TMath::Sqrt(TMath::Power(collision.posX() - candidate.x(), 2.0) + TMath::Power(collision.posY() - candidate.y(), 2.0) + TMath::Power(collision.posZ() - candidate.z(), 2.0));
    float CtauK0s = candidate.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * TDatabasePDG::Instance()->GetParticle(kK0Short)->Mass(); // FIXME: Get from the common header
    float lowmasscutks0 = 0.497 - 2.0 * cSigmaMassKs0;
    float highmasscutks0 = 0.497 + 2.0 * cSigmaMassKs0;

    if (TMath::Abs(CtauK0s) < 2.0 || TMath::Abs(CtauK0s) > cMaxV0LifeTime || candidate.mK0Short() < lowmasscutks0 || candidate.mK0Short() > highmasscutks0) {
      return false;
    }
    if (dcaDaughv0 > ConfV0DCADaughMax) {
      return false;
    }
    if (TMath::Abs(candidate.dcav0topv()) > cMaxV0DCA) {
      return false;
    }
    if (cpav0 < ConfV0CPAMin) {
      return false;
    }
    if (decaylength > 100) {
      return false;
    }
    if (tranRad > 100) {
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
  ROOT::Math::PxPyPzMVector Lambdac, Proton, Kshort, LambdacRot, KshortRot;
  // ROOT::Math::PxPyPzMVector fourVecDauCM;
  // ROOT::Math::XYZVector threeVecDauCM, threeVecDauCMXY, eventplaneVec, eventplaneVecNorm, beamvector;
  double massPr = TDatabasePDG::Instance()->GetParticle(kProton)->Mass();   // FIXME: Get from the common header
  double massK0s = TDatabasePDG::Instance()->GetParticle(kK0Short)->Mass(); // FIXME: Get from the common header
  double v2, v2Rot;
  void processSameEvent(EventCandidates::iterator const& collision, TrackCandidates const& tracks, AllTrackCandidates const&, ResoV0s const& V0s, aod::BCsWithTimestamps const&)
  {
    if (!collision.sel8()) {
      return;
    }
    auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
    if (runNumber != bc.runNumber()) {
      initCCDB(bc, runNumber, ccdb, isRun2 ? ccdbPathGrp : ccdbPathGrpMag, lut, isRun2);
      bz = o2::base::Propagator::Instance()->getNominalBz();
      // df.setBz(bz); /// put it outside the 'if'! Otherwise we have a difference wrt bz Configurable (< 1 permille) in Run2 conv. data
    }
    df.setBz(bz);
    auto centrality = collision.centFT0C();
    auto multTPC = collision.multNTracksPV();
    histos.fill(HIST("hFTOCvsTPCNoCut"), centrality, multTPC);
    if (!collision.triggereventep()) {
      return;
    }
    if (additionalEvSel && (!collision.selection_bit(aod::evsel::kNoSameBunchPileup) || !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV))) {
      return;
    }
    if (additionalEvSel2 && (!collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard))) {
      return;
    }
    int occupancy = collision.trackOccupancyInTimeRange();
    auto psiFT0C = collision.psiFT0C();
    auto psiFT0A = collision.psiFT0A();
    auto psiTPC = collision.psiTPC();
    auto psiTPCR = collision.psiTPCR();
    auto psiTPCL = collision.psiTPCL();

    // auto QFT0C = collision.qFT0C();
    // auto QFT0A = collision.qFT0A();
    // auto QTPC = collision.qTPC();
    // auto QTPCR = collision.qTPCR();
    // auto QTPCL = collision.qTPCL();

    histos.fill(HIST("hFTOCvsTPC"), centrality, multTPC);
    histos.fill(HIST("hFTOCvsTPCSelected"), centrality, multTPC);
    histos.fill(HIST("hPsiFT0C"), centrality, psiFT0C, occupancy);
    histos.fill(HIST("hPsiFT0A"), centrality, psiFT0A, occupancy);
    histos.fill(HIST("hPsiTPC"), centrality, psiTPC, occupancy);
    histos.fill(HIST("hPsiTPCR"), centrality, psiTPCR, occupancy);
    histos.fill(HIST("hPsiTPCL"), centrality, psiTPCL, occupancy);
    histos.fill(HIST("ResFT0CTPC"), centrality, TMath::Cos(2.0 * (psiFT0C - psiTPC)));
    histos.fill(HIST("ResFT0CTPCR"), centrality, TMath::Cos(2.0 * (psiFT0C - psiTPCR)));
    histos.fill(HIST("ResFT0CTPCL"), centrality, TMath::Cos(2.0 * (psiFT0C - psiTPCL)));
    histos.fill(HIST("ResTPCRTPCL"), centrality, TMath::Cos(2.0 * (psiTPCR - psiTPCL)));
    histos.fill(HIST("ResFT0CFT0A"), centrality, TMath::Cos(2.0 * (psiFT0C - psiFT0A)));
    histos.fill(HIST("ResFT0ATPC"), centrality, TMath::Cos(2.0 * (psiTPC - psiFT0A)));
    histos.fill(HIST("hCentrality"), centrality);
    histos.fill(HIST("hVtxZ"), collision.posZ());
    histos.fill(HIST("hOccupancy"), occupancy);
    auto firstprimarytrack = 0;
    for (auto track1 : tracks) {
      if (!selectionTrack(track1)) {
        continue;
      }

      if (!track1.hasTOF()) {
        if (!rejectPi(track1)) {
          continue;
        }
        if (!rejectEl(track1)) {
          continue;
        }
        if (!rejectKa(track1)) {
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
      histos.fill(HIST("hMomCorr"), track1.p() / track1.sign(), track1.p() - track1.tpcInnerParam(), centrality);
      histos.fill(HIST("hEta"), track1.eta());
      histos.fill(HIST("hDcaz"), track1.dcaZ());
      histos.fill(HIST("hNsigmaProtonTPCDiff"), track1.tpcNSigmaPr(), track1.tpcNSigmaKa(), track1.pt());
      histos.fill(HIST("hNsigmaProtonTPC"), track1.tpcNSigmaPr(), track1.pt());
      histos.fill(HIST("hNsigmaProtonTOF"), track1.tofNSigmaPr(), track1.pt());
      auto track1ID = track1.globalIndex();
      auto trackParCovBach = getTrackParCov(track1);
      for (auto v0 : V0s) {
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
        auto v0decaylength = TMath::Sqrt(TMath::Power(collision.posX() - v0.x(), 2.0) + TMath::Power(collision.posY() - v0.y(), 2.0) + TMath::Power(collision.posZ() - v0.z(), 2.0));
        if (firstprimarytrack == 0) {
          histos.fill(HIST("hV0decaylength"), v0decaylength);
        }
        float v0x, v0y, v0z, v0px, v0py, v0pz;
        // float v0PosPx, v0PosPy, v0PosPz, v0NegPx, v0NegPy, v0NegPz;
        // float dcaV0dau, dcaPosToPV, dcaNegToPV, v0cosPA;
        float posTrackX, negTrackX;
        o2::track::TrackParCov trackParCovV0DaughPos;
        o2::track::TrackParCov trackParCovV0DaughNeg;

        // const auto& trackV0DaughPos = v0.posTrack_as<aod::TracksWCov>();
        // const auto& trackV0DaughNeg = v0.negTrack_as<aod::TracksWCov>();

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
        // dcaV0dau = v0.dcaV0daughters();
        // dcaPosToPV = v0.dcapostopv();
        // dcaNegToPV = v0.dcanegtopv();
        // v0cosPA = v0.v0cosPA();

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
            // LOG(info) << "Vertexing succeeded for Lc candidate";
          }
        } catch (const std::runtime_error& error) {
          // LOG(info) << "Run time error found: " << error.what() << ". DCAFitterN cannot work, skipping the candidate.";
          continue;
        }
        if (firstprimarytrack == 0) {
          histos.fill(HIST("hInvMassKs0"), v0.mK0Short());
        }

        df.propagateTracksToVertex();           // propagate the bachelor and V0 to the Lambdac vertex
        trackV0.getPxPyPzGlo(pVecV0);           // momentum of D0 at the Lambdac vertex
        trackParCovBach.getPxPyPzGlo(pVecBach); // momentum of proton at the Lambdac vertex
        pVecCand = RecoDecay::pVec(pVecV0, pVecBach);

        const auto& secondaryVertex = df.getPCACandidate();
        auto chi2PCA = df.getChi2AtPCACandidate();
        if (chi2PCA > cutchi2PCA) {
          continue;
        }
        // auto covMatrixPCA = df.calcPCACovMatrixFlat();
        if (firstprimarytrack == 0) {
          histos.fill(HIST("hchi2PCA"), chi2PCA);
        }
        // get track impact parameters
        // This modifies track momenta!
        auto primaryVertex = getPrimaryVertex(collision);
        // auto covMatrixPV = primaryVertex.getCov();
        o2::dataformats::DCA impactParameter0;
        o2::dataformats::DCA impactParameter1;
        trackV0.propagateToDCA(primaryVertex, bz, &impactParameter0);
        trackParCovBach.propagateToDCA(primaryVertex, bz, &impactParameter1);

        // get uncertainty of the decay length
        double phi, theta;
        getPointDirection(std::array{collision.posX(), collision.posY(), collision.posZ()}, secondaryVertex, phi, theta);
        // auto errorDecayLength = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, theta) + getRotatedCovMatrixXX(covMatrixPCA, phi, theta));
        // auto errorDecayLengthXY = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, 0.) + getRotatedCovMatrixXX(covMatrixPCA, phi, 0.));

        Kshort = ROOT::Math::PxPyPzMVector(pVecV0[0], pVecV0[1], pVecV0[2], massK0s);
        Proton = ROOT::Math::PxPyPzMVector(pVecBach[0], pVecBach[1], pVecBach[2], massPr);
        Lambdac = Proton + Kshort;

        double protonimpactparameter = impactParameter1.getY();
        // double kshortimpactparameter=impactParameter0.getY();

        double decaylengthx = secondaryVertex[0] - collision.posX();
        double decaylengthy = secondaryVertex[1] - collision.posY();
        double decaylengthz = secondaryVertex[2] - collision.posZ();
        double decaylength = TMath::Sqrt(decaylengthx * decaylengthx + decaylengthy * decaylengthy + decaylengthz * decaylengthz);
        double decaylengthxy = TMath::Sqrt(decaylengthx * decaylengthx + decaylengthy * decaylengthy);
        // double lambdaclifetime = (decaylength * Lambdac.M()) / Lambdac.P();
        // if (lambdaclifetime > cutmaxctaulambdac) {
        //   continue;
        // }
        double anglesign = decaylengthx * Lambdac.Px() + decaylengthy * Lambdac.Py() + decaylengthz * Lambdac.Pz();
        double CPAlambdac = anglesign / (decaylength * Lambdac.P());
        anglesign = anglesign / TMath::Abs(anglesign);
        auto signprotonimpactparameter = protonimpactparameter * anglesign;
        if (firstprimarytrack == 0) {
          histos.fill(HIST("hDcaxy"), track1.dcaXY());
          histos.fill(HIST("hunsignprotonDca"), protonimpactparameter);
          histos.fill(HIST("hsignCPALambdac"), TMath::Abs(CPAlambdac));
          histos.fill(HIST("hsignprotonDca"), signprotonimpactparameter);
        }
        firstprimarytrack = firstprimarytrack + 1;
        auto phiminuspsi = GetPhiInRange(Lambdac.Phi() - psiFT0C);
        v2 = TMath::Cos(2.0 * phiminuspsi);
        // if (TMath::Abs(CPAlambdac) > cutCPAlambdac && Lambdac.M() > 2.18 && Lambdac.M() <= 2.42) {
        if (Lambdac.M() > 2.18 && Lambdac.M() <= 2.42 && TMath::Abs(Lambdac.Rapidity()) < confRapidity && Lambdac.Pt() > 2 && Lambdac.Pt() <= 6.0) {
          if (!useDecayLengthxy) {
            histos.fill(HIST("hSparseV2SASameEvent_V2_EP"), Lambdac.M(), Lambdac.Pt(), v2, Proton.Pt(), decaylength, TMath::Abs(track1.dcaXY()));
          }
          if (useDecayLengthxy) {
            histos.fill(HIST("hSparseV2SASameEvent_V2_EP"), Lambdac.M(), Lambdac.Pt(), v2, Proton.Pt(), decaylengthxy, TMath::Abs(track1.dcaXY()));
          }
        }
        for (int nrotbkg = 0; nrotbkg < nBkgRotations; nrotbkg++) {
          auto anglestart = confMinRot;
          auto angleend = confMaxRot;
          auto anglestep = (angleend - anglestart) / (1.0 * (nBkgRotations - 1));
          auto rotangle = anglestart + nrotbkg * anglestep;
          histos.fill(HIST("hRotation"), rotangle);
          auto rotKshortPx = Kshort.Px() * std::cos(rotangle) - Kshort.Py() * std::sin(rotangle);
          auto rotKshortPy = Kshort.Px() * std::sin(rotangle) + Kshort.Py() * std::cos(rotangle);
          KshortRot = ROOT::Math::PxPyPzMVector(rotKshortPx, rotKshortPy, Kshort.pz(), massK0s);
          LambdacRot = Proton + KshortRot;
          auto phiminuspsiRot = GetPhiInRange(LambdacRot.Phi() - psiFT0C);
          v2Rot = TMath::Cos(2.0 * phiminuspsiRot);
          // double CPAlambdacRot = (decaylengthx * LambdacRot.Px() + decaylengthy * LambdacRot.Py() + decaylengthz * LambdacRot.Pz()) / (decaylength * LambdacRot.P());
          // if (TMath::Abs(CPAlambdacRot) > cutCPAlambdac && LambdacRot.M() > 2.18 && LambdacRot.M() <= 2.42) {
          if (LambdacRot.M() > 2.18 && LambdacRot.M() <= 2.42 && TMath::Abs(LambdacRot.Rapidity()) < confRapidity && LambdacRot.Pt() > 2 && LambdacRot.Pt() <= 6.0) {
            if (!useDecayLengthxy) {
              histos.fill(HIST("hSparseV2SASameEventRotational_V2_EP"), LambdacRot.M(), LambdacRot.Pt(), v2Rot, Proton.Pt(), decaylength, TMath::Abs(track1.dcaXY()));
            }
            if (useDecayLengthxy) {
              histos.fill(HIST("hSparseV2SASameEventRotational_V2_EP"), LambdacRot.M(), LambdacRot.Pt(), v2Rot, Proton.Pt(), decaylengthxy, TMath::Abs(track1.dcaXY()));
            }
          }
        }

        // ROOT::Math::Boost boost{Lambdac.BoostToCM()};
        // fourVecDauCM = boost(Kshort);
        // threeVecDauCM = fourVecDauCM.Vect();
        // threeVecDauCMXY = ROOT::Math::XYZVector(threeVecDauCM.X(), threeVecDauCM.Y(), 0.);
        // beamvector = ROOT::Math::XYZVector(0, 0, 1);
        // auto cosThetaStar = beamvector.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(beamvector.Mag2());
        // auto SA = cosThetaStar * TMath::Sin(2.0 * phiminuspsi);
        // if (fillPolarization) {
        // if (track1.sign() > 0) {
        // histos.fill(HIST("hSparseV2SASameEventplus_SA"), Lambdac.M(), Lambdac.Pt(), cosThetaStar, phiminuspsi, centrality);
        // histos.fill(HIST("hSparseV2SASameEventplus_SA_A0"), Lambdac.M(), Lambdac.Pt(), cosThetaStar * cosThetaStar, phiminuspsi, centrality);
        // histos.fill(HIST("hSparseV2SASameEventplus_SA_azimuth"), Lambdac.M(), Lambdac.Pt(), SA, centrality);
        // }
        // if (track1.sign() < 0) {
        // histos.fill(HIST("hSparseV2SASameEventminus_SA"), Lambdac.M(), Lambdac.Pt(), cosThetaStar, phiminuspsi, centrality);
        // histos.fill(HIST("hSparseV2SASameEventminus_SA_A0"), Lambdac.M(), Lambdac.Pt(), cosThetaStar * cosThetaStar, phiminuspsi, centrality);
        // histos.fill(HIST("hSparseV2SASameEventminus_SA_azimuth"), Lambdac.M(), Lambdac.Pt(), SA, centrality);
        // }
        // }
      }
    }
  }
  PROCESS_SWITCH(highmasslambdasvx, processSameEvent, "Process Same event", true);
  /* void processMixedEventOpti(EventCandidates const& collisions, TrackCandidates const& tracks, AllTrackCandidates const&, ResoV0s const& V0s)
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
      if (additionalEvSel2 && (!collision1.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard))) {
        return;
      }
      if (additionalEvSel2 && (!collision2.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard))) {
        return;
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
          dcasum = TMath::Sqrt((track1.dcaXY() + (v0.dcav0topv())) * (track1.dcaXY() + (v0.dcav0topv())));
        }
        // auto diffangle = Proton.Phi() - Lambdac.Phi();
        // auto decaylength = TMath::Abs((track1.dcaXY() / TMath::Sin(diffangle)) / (Lambdac.P() / 2.286));
        // auto dcasum = TMath::Sqrt(track1.dcaXY() * track1.dcaXY() + v0.dcav0topv() * v0.dcav0topv());
        if (fillDefault && Lambdac.M() > 2.18 && Lambdac.M() <= 2.42) {
          if (fillDecayLength) {
            histos.fill(HIST("hSparseV2SAMixedEvent_V2"), Lambdac.M(), Lambdac.Pt(), v2, dcasum);
          }
          histos.fill(HIST("hSparseV2SAMixedEvent_V2_new"), Lambdac.M(), Lambdac.Pt(), v2, TMath::Abs(track1.dcaXY()), Proton.Pt());
        }
        if (fillOccupancy && Lambdac.M() > 2.18 && Lambdac.M() <= 2.42) {
          if (fillDecayLength) {
            histos.fill(HIST("hSparseV2SAMixedEvent_V2_occupancy"), Lambdac.M(), Lambdac.Pt(), v2, dcasum, TMath::Abs(track1.dcaXY()), occupancy);
          }
          histos.fill(HIST("hSparseV2SAMixedEvent_V2_new_occupancy"), Lambdac.M(), Lambdac.Pt(), v2, TMath::Abs(track1.dcaXY()), Proton.Pt(), occupancy);
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
  PROCESS_SWITCH(highmasslambdasvx, processMixedEventOpti, "Process Mixed event new", true);*/
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<highmasslambdasvx>(cfgc, TaskName{"highmasslambdasvx"})};
}
