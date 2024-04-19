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
  Configurable<bool> fillPolarization{"fillPolarization", false, "fill polarization"};

  // events
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCutCentralityMax{"cfgCutCentralityMax", 80.0f, "Accepted maximum Centrality"};
  Configurable<float> cfgCutCentralityMin{"cfgCutCentralityMin", 20.0f, "Accepted minimum Centrality"};

  // proton track cut
  Configurable<float> confRapidity{"confRapidity", 0.5, "cut on Rapidity"};
  Configurable<float> cfgCutPT{"cfgCutPT", 0.3, "PT cut on daughter track"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8, "Eta cut on daughter track"};
  Configurable<float> cfgCutMinDCAxy{"cfgCutMinDCAxy", 0.0f, "Minimum DCAxy range for proton"};
  Configurable<float> cfgCutDCAxy{"cfgCutDCAxy", 0.1f, "DCAxy range for tracks"};
  Configurable<float> cfgCutDCAz{"cfgCutDCAz", 1.0f, "DCAz range for tracks"};
  Configurable<int> cfgITScluster{"cfgITScluster", 0, "Number of ITS cluster"};
  Configurable<int> cfgTPCcluster{"cfgTPCcluster", 70, "Number of TPC cluster"};
  Configurable<bool> ispTdepPID{"ispTdepPID", true, "pT dependent PID"};
  Configurable<float> nsigmaCutTPC{"nsigmacutTPC", 3.0, "Value of the TPC Nsigma cut"};
  Configurable<float> nsigmaCutTOF{"nsigmacutTOF", 3.0, "Value of the TOF Nsigma cut"};
  Configurable<float> nsigmaCutTPCPre{"nsigmacutTPCPre", 5.0, "Value of the TPC Nsigma cut Pre filter"};
  // Configs for V0
  Configurable<float> ConfV0PtMin{"ConfV0PtMin", 0.f, "Minimum transverse momentum of V0"};
  Configurable<double> ConfV0DCADaughMax{"ConfV0DCADaughMax", 0.4f, "Maximum DCA between the V0 daughters"};
  Configurable<double> ConfV0CPAMin{"ConfV0CPAMin", 0.996f, "Minimum CPA of V0"};
  Configurable<float> ConfV0TranRadV0Min{"ConfV0TranRadV0Min", 1.5f, "Minimum transverse radius"};
  Configurable<float> ConfV0TranRadV0Max{"ConfV0TranRadV0Max", 100.f, "Maximum transverse radius"};
  Configurable<double> cMaxV0DCA{"cMaxV0DCA", 0.5, "Maximum V0 DCA to PV"};
  Configurable<float> cMaxV0LifeTime{"cMaxV0LifeTime", 20, "Maximum V0 life time"};
  Configurable<float> cSigmaMassKs0{"cSigmaMassKs0", 0.006, "Sigma cut on KS0 mass"};

  // config for V0 daughters
  Configurable<float> ConfDaughEta{"ConfDaughEta", 0.8f, "V0 Daugh sel: max eta"};
  Configurable<float> ConfDaughPt{"ConfDaughPt", 0.1f, "V0 Daugh sel: min pt"};
  Configurable<float> ConfDaughTPCnclsMin{"ConfDaughTPCnclsMin", 70.f, "V0 Daugh sel: Min. nCls TPC"};
  Configurable<double> ConfDaughDCAMin{"ConfDaughDCAMin", 0.08f, "V0 Daugh sel:  Max. DCA Daugh to PV (cm)"};
  Configurable<float> ConfDaughPIDCuts{"ConfDaughPIDCuts", 3, "PID selections for KS0 daughters"};

  // Mixed event
  Configurable<int> cfgNoMixedEvents{"cfgNoMixedEvents", 1, "Number of mixed events per event"};

  // THnsparse bining
  ConfigurableAxis configThnAxisInvMass{"configThnAxisInvMass", {60, 2.15, 2.45}, "#it{M} (GeV/#it{c}^{2})"};
  ConfigurableAxis configThnAxisPt{"configThnAxisPt", {30, 1.0, 31.}, "#it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis configThnAxisCosThetaStar{"configThnAxisCosThetaStar", {10, -1.0, 1.}, "cos(#vartheta)"};
  ConfigurableAxis configThnAxisCentrality{"configThnAxisCentrality", {8, 0., 80}, "Centrality"};
  ConfigurableAxis configThnAxisPhiminusPsi{"configThnAxisPhiminusPsi", {6, 0.0, TMath::Pi()}, "#phi - #psi"};
  ConfigurableAxis configThnAxisV2{"configThnAxisV2", {100, -1, 1}, "V2"};
  ConfigurableAxis configThnAxisSA{"configThnAxisSA", {100, -1, 1}, "SA"};

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter centralityFilter = (nabs(aod::cent::centFT0C) < cfgCutCentralityMax && nabs(aod::cent::centFT0C) > cfgCutCentralityMin);
  Filter acceptanceFilter = (nabs(aod::track::eta) < cfgCutEta && nabs(aod::track::pt) > cfgCutPT);
  Filter dcaCutFilter = (nabs(aod::track::dcaXY) < cfgCutDCAxy) && (nabs(aod::track::dcaZ) < cfgCutDCAz);
  Filter pidFilter = nabs(aod::pidtpc::tpcNSigmaPr) < nsigmaCutTPCPre;

  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::FV0Mults, aod::TPCMults, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0As, aod::EPCalibrationTables, aod::Mults>>;
  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPr, aod::pidTOFFullPr>>;
  using AllTrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPi>;
  using ResoV0s = aod::V0Datas;

  SliceCache cache;
  // Partition<TrackCandidates> posTracks = aod::track::signed1Pt > cfgCutCharge;
  // Partition<TrackCandidates> negTracks = aod::track::signed1Pt < cfgCutCharge;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext&)
  {
    const AxisSpec thnAxisInvMass{configThnAxisInvMass, "#it{M} (GeV/#it{c}^{2})"};
    const AxisSpec thnAxisPt{configThnAxisPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec thnAxisCosThetaStar{configThnAxisCosThetaStar, "cos(#vartheta)"};
    const AxisSpec thnAxisPhiminusPsi{configThnAxisPhiminusPsi, "#phi - #psi"};
    const AxisSpec thnAxisCentrality{configThnAxisCentrality, "Centrality (%)"};
    const AxisSpec thnAxisV2{configThnAxisV2, "V2"};
    const AxisSpec thnAxisSA{configThnAxisSA, "SA"};

    AxisSpec phiAxis = {500, -6.28, 6.28, "phi"};
    AxisSpec resAxis = {400, -2, 2, "Res"};
    AxisSpec centAxis = {8, 0, 80, "V0M (%)"};
    AxisSpec dcaAxis = {100, 0.0, 0.1, "V0M (%)"};
    AxisSpec dcatoPVAxis = {10, 0.0, 0.5, "V0M (%)"};

    histos.add("hInvMassKs0", "hInvMassKs0", kTH1F, {{200, 0.4f, 0.6f}});
    histos.add("hV0Dca", "hV0Dca", kTH1F, {{2000, -1.0f, 1.0f}});
    histos.add("hpTvsRapidity", "pT vs Rapidity", kTH2F, {{100, 0.0f, 10.0f}, {300, -1.5f, 1.5f}});
    histos.add("hFTOCvsTPCNoCut", "Mult correlation FT0C vs. TPC without any cut", kTH2F, {{80, 0.0f, 80.0f}, {100, -0.5f, 5999.5f}});
    histos.add("hFTOCvsTPC", "Mult correlation FT0C vs. TPC", kTH2F, {{80, 0.0f, 80.0f}, {100, -0.5f, 5999.5f}});
    histos.add("hFTOCvsTPCSelected", "Mult correlation FT0C vs. TPC after selection", kTH2F, {{80, 0.0f, 80.0f}, {100, -0.5f, 5999.5f}});
    histos.add("hCentrality", "Centrality distribution", kTH1F, {{200, 0.0, 200.0}});
    histos.add("hVtxZ", "Vertex distribution in Z;Z (cm)", kTH1F, {{400, -20.0, 20.0}});
    histos.add("hEta", "Eta distribution", kTH1F, {{200, -1.0f, 1.0f}});
    histos.add("hDcaxy", "Dcaxy distribution", kTH1F, {{1000, -0.5f, 0.5f}});
    histos.add("hDcaz", "Dcaz distribution", kTH1F, {{1000, -0.5f, 0.5f}});
    histos.add("hNsigmaProtonTPC", "NsigmaProton TPC distribution", kTH1F, {{200, -10.0f, 10.0f}});
    histos.add("hNsigmaProtonTOF", "NsigmaProton TOF distribution", kTH1F, {{200, -10.0f, 10.0f}});
    histos.add("hPsiFT0C", "PsiFT0C", kTH2F, {centAxis, phiAxis});
    histos.add("hPsiFT0A", "PsiFT0A", kTH2F, {centAxis, phiAxis});
    histos.add("hPsiTPC", "PsiTPC", kTH2F, {centAxis, phiAxis});
    histos.add("hPsiTPCR", "PsiTPCR", kTH2F, {centAxis, phiAxis});
    histos.add("hPsiTPCL", "PsiTPCL", kTH2F, {centAxis, phiAxis});

    histos.add("hSparseV2SASameEvent_V2", "hSparseV2SASameEvent_V2", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisV2, dcaAxis, dcatoPVAxis, thnAxisCentrality});
    histos.add("hSparseV2SAMixedEvent_V2", "hSparseV2SAMixedEvent_V2", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisV2, dcaAxis, dcatoPVAxis, thnAxisCentrality});

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
    if (!(candidate.isGlobalTrack() && candidate.isPVContributor() && candidate.itsNCls() > cfgITScluster && candidate.tpcNClsFound() > cfgTPCcluster && TMath::Abs(candidate.dcaXY()) > cfgCutMinDCAxy)) {
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
    if (candidate.p() > 0.8 && candidate.hasTOF() && TMath::Sqrt(candidate.tofNSigmaPr() * candidate.tofNSigmaPr() + candidate.tpcNSigmaPr() * candidate.tpcNSigmaPr()) < 2.0 * nsigmaCutTOF * nsigmaCutTOF) {
      return true;
    }
    return false;
  }

  template <typename T>
  bool selectionPID(const T& candidate)
  {
    if (candidate.p() <= 0.5 && TMath::Abs(candidate.tpcNSigmaPr()) < 5.0) {
      return true;
    }
    if (candidate.p() > 0.5 && candidate.p() <= 0.8 && TMath::Abs(candidate.tpcNSigmaPr()) < 3.0) {
      return true;
    }
    if (candidate.p() > 0.8 && !candidate.hasTOF() && TMath::Abs(candidate.tpcNSigmaPr()) < nsigmaCutTPC) {
      return true;
    }
    if (candidate.p() > 0.8 && candidate.hasTOF() && TMath::Sqrt(candidate.tofNSigmaPr() * candidate.tofNSigmaPr() + candidate.tpcNSigmaPr() * candidate.tpcNSigmaPr()) < 2.0 * nsigmaCutTOF * nsigmaCutTOF) {
      return true;
    }
    return false;
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
    const double dcaDaughv0 = candidate.dcaV0daughters();
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
  ROOT::Math::PxPyPzMVector Lambdac, Proton, Kshort, fourVecDauCM;
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
    auto psiFT0C = collision.psiFT0C();
    auto psiFT0A = collision.psiFT0A();
    auto psiTPC = collision.psiTPC();
    auto psiTPCR = collision.psiTPCR();
    auto psiTPCL = collision.psiTPCL();
    histos.fill(HIST("hFTOCvsTPC"), centrality, multTPC);
    histos.fill(HIST("hFTOCvsTPCSelected"), centrality, multTPC);
    histos.fill(HIST("hPsiFT0C"), centrality, psiFT0C);
    histos.fill(HIST("hPsiFT0A"), centrality, psiFT0A);
    histos.fill(HIST("hPsiTPC"), centrality, psiTPC);
    histos.fill(HIST("hPsiTPCR"), centrality, psiTPCR);
    histos.fill(HIST("hPsiTPCL"), centrality, psiTPCL);
    histos.fill(HIST("ResFT0CTPC"), centrality, TMath::Cos(2.0 * (psiFT0C - psiTPC)));
    histos.fill(HIST("ResFT0CTPCR"), centrality, TMath::Cos(2.0 * (psiFT0C - psiTPCR)));
    histos.fill(HIST("ResFT0CTPCL"), centrality, TMath::Cos(2.0 * (psiFT0C - psiTPCL)));
    histos.fill(HIST("ResTPCRTPCL"), centrality, TMath::Cos(2.0 * (psiTPCR - psiTPCL)));
    histos.fill(HIST("ResFT0CFT0A"), centrality, TMath::Cos(2.0 * (psiFT0C - psiFT0A)));
    histos.fill(HIST("ResFT0ATPC"), centrality, TMath::Cos(2.0 * (psiTPC - psiFT0A)));
    histos.fill(HIST("hCentrality"), centrality);
    histos.fill(HIST("hVtxZ"), collision.posZ());
    auto firstprimarytrack = 0;
    for (auto track1 : tracks) {
      if (!selectionTrack(track1)) {
        continue;
      }
      // PID check
      if (ispTdepPID && !selectionPIDpTdependent(track1)) {
        continue;
      }
      if (!ispTdepPID && !selectionPID(track1)) {
        continue;
      }
      histos.fill(HIST("hEta"), track1.eta());
      histos.fill(HIST("hDcaxy"), track1.dcaXY());
      histos.fill(HIST("hDcaz"), track1.dcaZ());
      histos.fill(HIST("hNsigmaProtonTPC"), track1.tpcNSigmaPr());
      histos.fill(HIST("hNsigmaProtonTOF"), track1.tofNSigmaPr());
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
        auto v2 = TMath::Cos(2.0 * phiminuspsi);
        histos.fill(HIST("hSparseV2SASameEvent_V2"), Lambdac.M(), Lambdac.Pt(), v2, std::abs(track1.dcaXY()), std::abs(v0.dcav0topv()), centrality);

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
      auto centrality = collision1.centFT0C();
      auto psiFT0C = collision1.psiFT0C();

      for (auto& [track1, v0] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
        if (!selectionTrack(track1)) {
          continue;
        }
        // PID check
        if (ispTdepPID && !selectionPIDpTdependent(track1)) {
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
        auto v2 = TMath::Cos(2.0 * phiminuspsi);
        histos.fill(HIST("hSparseV2SAMixedEvent_V2"), Lambdac.M(), Lambdac.Pt(), v2, std::abs(track1.dcaXY()), std::abs(v0.dcav0topv()), centrality);

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
