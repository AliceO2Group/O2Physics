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

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;
struct phipbpb {

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

  // events
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  // track
  Configurable<float> cfgCutCharge{"cfgCutCharge", 0.0, "cut on Charge"};
  Configurable<float> cfgCutPT{"cfgCutPT", 0.2, "PT cut on daughter track"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8, "Eta cut on daughter track"};
  Configurable<float> cfgCutDCAxy{"cfgCutDCAxy", 2.0f, "DCAxy range for tracks"};
  Configurable<float> cfgCutDCAz{"cfgCutDCAz", 2.0f, "DCAz range for tracks"};
  Configurable<float> nsigmaCutTPC{"nsigmacutTPC", 3.0, "Value of the TPC Nsigma cut"};
  Configurable<float> nsigmaCutCombined{"nsigmaCutCombined", 3.0, "Value of the TOF Nsigma cut"};
  Configurable<int> cfgNoMixedEvents{"cfgNoMixedEvents", 1, "Number of mixed events per event"};
  Configurable<int> cfgITScluster{"cfgITScluster", 0, "Number of ITS cluster"};
  Configurable<bool> isDeepAngle{"isDeepAngle", false, "Deep Angle cut"};
  Configurable<bool> ispTdepPID{"ispTdepPID", true, "pT dependent PID"};
  Configurable<double> cfgDeepAngle{"cfgDeepAngle", 0.04, "Deep Angle cut value"};
  Configurable<double> confRapidity{"confRapidity", 0.5, "Rapidity cut"};
  ConfigurableAxis configThnAxisInvMass{"configThnAxisInvMass", {120, 0.98, 1.1}, "#it{M} (GeV/#it{c}^{2})"};
  ConfigurableAxis configThnAxisPt{"configThnAxisPt", {100, 0.0, 10.}, "#it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis configThnAxisCosThetaStar{"configThnAxisCosThetaStar", {10, 0.0, 1.}, "cos(#vartheta)"};
  ConfigurableAxis configThnAxisCentrality{"configThnAxisCentrality", {8, 0., 80}, "Centrality"};
  ConfigurableAxis configThnAxisPhiminusPsi{"configThnAxisPhiminusPsi", {6, 0.0, TMath::Pi()}, "#phi - #psi"};
  ConfigurableAxis configThnAxisV2{"configThnAxisV2", {200, -1, 1}, "V2"};
  ConfigurableAxis configThnAxisSA{"configThnAxisSA", {200, -1, 1}, "SA"};
  Configurable<bool> additionalEvsel{"additionalEvsel", false, "Additional event selcection"};
  Configurable<bool> timFrameEvsel{"timFrameEvsel", false, "TPC Time frame boundary cut"};

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter acceptanceFilter = (nabs(aod::track::eta) < cfgCutEta && nabs(aod::track::pt) > cfgCutPT);
  Filter DCAcutFilter = (nabs(aod::track::dcaXY) < cfgCutDCAxy) && (nabs(aod::track::dcaZ) < cfgCutDCAz);
  Filter PIDcutFilter = nabs(aod::pidtpc::tpcNSigmaKa) < nsigmaCutTPC;

  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::FV0Mults, aod::TPCMults, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0As, aod::EPCalibrationTables, aod::Mults>>;
  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullKa, aod::pidTOFFullKa>>;

  using EventCandidatesMC = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>;
  using TrackCandidatesMC = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::McTrackLabels>>;

  SliceCache cache;
  Partition<TrackCandidates> posTracks = aod::track::signed1Pt > cfgCutCharge;
  Partition<TrackCandidates> negTracks = aod::track::signed1Pt < cfgCutCharge;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Event selection cuts - Alex
  TF1* fMultPVCutLow = nullptr;
  TF1* fMultPVCutHigh = nullptr;
  TF1* fMultCutLow = nullptr;
  TF1* fMultCutHigh = nullptr;
  TF1* fMultMultPVCut = nullptr;

  void init(o2::framework::InitContext&)
  {
    const AxisSpec thnAxisInvMass{configThnAxisInvMass, "#it{M} (GeV/#it{c}^{2})"};
    const AxisSpec thnAxisPt{configThnAxisPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec thnAxisCosThetaStarOP{configThnAxisCosThetaStar, "cos(#vartheta_{OP})"};
    const AxisSpec thnAxisCosThetaStarIP{configThnAxisCosThetaStar, "cos(#vartheta_{IP})"};
    const AxisSpec thnAxisPhiminusPsi{configThnAxisPhiminusPsi, "#phi - #psi"};
    const AxisSpec thnAxisCentrality{configThnAxisCentrality, "Centrality (%)"};
    const AxisSpec thnAxisV2{configThnAxisV2, "V2"};
    const AxisSpec thnAxisSA{configThnAxisSA, "SA"};
    AxisSpec phiAxis = {500, -6.28, 6.28, "phi"};
    AxisSpec resAxis = {400, -2, 2, "Res"};
    AxisSpec centAxis = {8, 0, 80, "V0M (%)"};

    histos.add("hpTvsRapidity", "pT vs Rapidity", kTH2F, {{100, 0.0f, 10.0f}, {300, -1.5f, 1.5f}});
    histos.add("hFTOCvsTPC", "Mult correlation FT0C vs. TPC", kTH2F, {{80, 0.0f, 80.0f}, {100, -0.5f, 5999.5f}});
    histos.add("hFTOCvsTPCSelected", "Mult correlation FT0C vs. TPC after selection", kTH2F, {{80, 0.0f, 80.0f}, {100, -0.5f, 5999.5f}});
    histos.add("hCentrality", "Centrality distribution", kTH1F, {{200, 0.0, 200.0}});
    histos.add("hVtxZ", "Vertex distribution in Z;Z (cm)", kTH1F, {{400, -20.0, 20.0}});
    histos.add("hEta", "Eta distribution", kTH1F, {{200, -1.0f, 1.0f}});
    histos.add("hDcaxy", "Dcaxy distribution", kTH1F, {{200, -1.0f, 1.0f}});
    histos.add("hDcaz", "Dcaz distribution", kTH1F, {{200, -1.0f, 1.0f}});
    histos.add("hNsigmaKaonTPC", "NsigmaKaon TPC distribution", kTH1F, {{200, -10.0f, 10.0f}});
    histos.add("hNsigmaKaonTOF", "NsigmaKaon TOF distribution", kTH1F, {{200, -10.0f, 10.0f}});
    histos.add("hPsiFT0C", "PsiFT0C", kTH2F, {centAxis, phiAxis});
    histos.add("hPsiFT0A", "PsiFT0A", kTH2F, {centAxis, phiAxis});
    histos.add("hPsiTPC", "PsiTPC", kTH2F, {centAxis, phiAxis});
    histos.add("hSparseV2SASameEvent_costhetastarOP", "hSparseV2SASameEvent_costhetastarOP", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisCosThetaStarOP, thnAxisPhiminusPsi, thnAxisCentrality});
    histos.add("hSparseV2SASameEvent_costhetastarIP", "hSparseV2SASameEvent_costhetastarIP", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisCosThetaStarOP, thnAxisPhiminusPsi, thnAxisCentrality});
    histos.add("hSparseV2SASameEvent_SA", "hSparseV2SASameEvent_SA", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisSA, thnAxisPhiminusPsi, thnAxisCentrality});
    histos.add("hSparseV2SASameEvent_V2", "hSparseV2SASameEvent_V2", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisV2, thnAxisCentrality});

    histos.add("hSparseV2SAMixedEvent_costhetastarOP", "hSparseV2SAMixedEvent_costhetastarOP", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisCosThetaStarOP, thnAxisPhiminusPsi, thnAxisCentrality});
    histos.add("hSparseV2SAMixedEvent_costhetastarIP", "hSparseV2SAMixedEvent_costhetastarIP", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisCosThetaStarOP, thnAxisPhiminusPsi, thnAxisCentrality});
    histos.add("hSparseV2SAMixedEvent_SA", "hSparseV2SAMixedEvent_SA", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisSA, thnAxisPhiminusPsi, thnAxisCentrality});
    histos.add("hSparseV2SAMixedEvent_V2", "hSparseV2SAMixedEvent_V2", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisV2, thnAxisCentrality});

    // histogram for resolution
    histos.add("ResFT0CTPC", "ResFT0CTPC", kTH2F, {centAxis, resAxis});
    histos.add("ResFT0CFT0A", "ResFT0CFT0A", kTH2F, {centAxis, resAxis});
    histos.add("ResFT0ATPC", "ResFT0ATPC", kTH2F, {centAxis, resAxis});

    // Event selection cut additional - Alex
    if (additionalEvsel) {
      fMultPVCutLow = new TF1("fMultPVCutLow", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 2.5*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x)", 0, 100);
      fMultPVCutLow->SetParameters(2834.66, -87.0127, 0.915126, -0.00330136, 332.513, -12.3476, 0.251663, -0.00272819, 1.12242e-05);
      fMultPVCutHigh = new TF1("fMultPVCutHigh", "[0]+[1]*x+[2]*x*x+[3]*x*x*x + 2.5*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x)", 0, 100);
      fMultPVCutHigh->SetParameters(2834.66, -87.0127, 0.915126, -0.00330136, 332.513, -12.3476, 0.251663, -0.00272819, 1.12242e-05);
      fMultCutLow = new TF1("fMultCutLow", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 2.5*([4]+[5]*x)", 0, 100);
      fMultCutLow->SetParameters(1893.94, -53.86, 0.502913, -0.0015122, 109.625, -1.19253);
      fMultCutHigh = new TF1("fMultCutHigh", "[0]+[1]*x+[2]*x*x+[3]*x*x*x + 3.*([4]+[5]*x)", 0, 100);
      fMultCutHigh->SetParameters(1893.94, -53.86, 0.502913, -0.0015122, 109.625, -1.19253);
      fMultMultPVCut = new TF1("fMultMultPVCut", "[0]+[1]*x+[2]*x*x", 0, 5000);
      fMultMultPVCut->SetParameters(-0.1, 0.785, -4.7e-05);
    }
  }

  double massKa = o2::constants::physics::MassKPlus;

  template <typename TCollision>
  bool eventSelected(TCollision collision, const int& multTrk, const float& centrality)
  {
    if (collision.alias_bit(kTVXinTRD)) {
      // TRD triggered
      // return 0;
    }
    auto multNTracksPV = collision.multNTracksPV();
    if (multNTracksPV < fMultPVCutLow->Eval(centrality))
      return 0;
    if (multNTracksPV > fMultPVCutHigh->Eval(centrality))
      return 0;
    // if (multTrk < fMultCutLow->Eval(centrality))
    //  return 0;
    // if (multTrk > fMultCutHigh->Eval(centrality))
    //  return 0;
    // if (multTrk > fMultMultPVCut->Eval(multNTracksPV))
    //  return 0;

    return 1;
  }

  template <typename T>
  bool selectionTrack(const T& candidate)
  {
    if (!(candidate.isGlobalTrack() || candidate.isPVContributor() || candidate.itsNCls() > cfgITScluster)) {
      return false;
    }
    return true;
  }
  template <typename T>
  bool selectionPIDpTdependent(const T& candidate)
  {
    if (candidate.pt() < 0.5 && TMath::Abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC) {
      return true;
    }
    if (candidate.pt() >= 0.5 && candidate.hasTOF() && ((candidate.tofNSigmaKa() * candidate.tofNSigmaKa()) + (candidate.tpcNSigmaKa() * candidate.tpcNSigmaKa())) < (nsigmaCutCombined * nsigmaCutCombined)) {
      return true;
    }
    return false;
  }
  template <typename T>
  bool selectionPID(const T& candidate)
  {
    if (!candidate.hasTOF() && TMath::Abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC) {
      return true;
    }
    if (candidate.hasTOF() && ((candidate.tofNSigmaKa() * candidate.tofNSigmaKa()) + (candidate.tpcNSigmaKa() * candidate.tpcNSigmaKa())) < (nsigmaCutCombined * nsigmaCutCombined)) {
      return true;
    }
    return false;
  }
  // deep angle cut on pair to remove photon conversion
  template <typename T1, typename T2>
  bool selectionPair(const T1& candidate1, const T2& candidate2)
  {
    double pt1, pt2, pz1, pz2, p1, p2, angle;
    pt1 = candidate1.pt();
    pt2 = candidate2.pt();
    pz1 = candidate1.pz();
    pz2 = candidate2.pz();
    p1 = candidate1.p();
    p2 = candidate2.p();
    angle = TMath::ACos((pt1 * pt2 + pz1 * pz2) / (p1 * p2));
    if (isDeepAngle && angle < cfgDeepAngle) {
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

  double GetDeltaPsiSubInRange(double psi1, double psi2)
  {
    double delta = psi1 - psi2;
    if (TMath::Abs(delta) > TMath::Pi() / 2) {
      if (delta > 0.)
        delta -= 2. * TMath::Pi() / 2;
      else
        delta += 2. * TMath::Pi() / 2;
    }
    return delta;
  }

  ConfigurableAxis axisVertex{"axisVertex", {20, -10, 10}, "vertex axis for bin"};
  ConfigurableAxis axisMultiplicityClass{"axisMultiplicityClass", {20, 0, 100}, "multiplicity percentile for bin"};
  ConfigurableAxis axisEPAngle{"axisEPAngle", {6, -TMath::Pi() / 2, TMath::Pi() / 2}, "event plane angle"};

  using BinningTypeVertexContributor = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C, aod::epcalibrationtable::PsiFT0C>;
  ROOT::Math::PxPyPzMVector PhiMesonMother, KaonPlus, KaonMinus, fourVecDauCM;
  ROOT::Math::XYZVector threeVecDauCM, threeVecDauCMXY, eventplaneVec, eventplaneVecNorm;

  void processSameEvent(EventCandidates::iterator const& collision, TrackCandidates const& tracks, aod::BCs const&)
  {
    if (!collision.sel8()) {
      return;
    }

    if (!collision.triggereventep()) {
      return;
    }
    if (timFrameEvsel && !collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
      return;
    }
    auto posThisColl = posTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto negThisColl = negTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto centrality = collision.centFT0C();
    auto multTPC = collision.multNTracksPV();
    auto psiFT0C = collision.psiFT0C();
    auto psiFT0A = collision.psiFT0A();
    auto psiTPC = collision.psiTPC();
    histos.fill(HIST("hFTOCvsTPC"), centrality, multTPC);
    if (additionalEvsel && !eventSelected(collision, tracks.size(), centrality)) {
      return;
    }
    histos.fill(HIST("hFTOCvsTPCSelected"), centrality, multTPC);
    histos.fill(HIST("hPsiFT0C"), centrality, psiFT0C);
    histos.fill(HIST("hPsiFT0A"), centrality, psiFT0A);
    histos.fill(HIST("hPsiTPC"), centrality, psiTPC);
    histos.fill(HIST("ResFT0CTPC"), centrality, TMath::Cos(2.0 * (psiFT0C - psiTPC)));
    histos.fill(HIST("ResFT0CFT0A"), centrality, TMath::Cos(2.0 * (psiFT0C - psiFT0A)));
    histos.fill(HIST("ResFT0ATPC"), centrality, TMath::Cos(2.0 * (psiTPC - psiFT0A)));
    histos.fill(HIST("hCentrality"), centrality);
    histos.fill(HIST("hVtxZ"), collision.posZ());
    for (auto track1 : posThisColl) {
      // track selection
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
      histos.fill(HIST("hNsigmaKaonTPC"), track1.tpcNSigmaKa());
      histos.fill(HIST("hNsigmaKaonTOF"), track1.tofNSigmaKa());
      auto track1ID = track1.globalIndex();
      for (auto track2 : negThisColl) {
        // track selection
        if (!selectionTrack(track2)) {
          continue;
        }
        // PID check
        if (ispTdepPID && !selectionPIDpTdependent(track2)) {
          continue;
        }
        if (!ispTdepPID && !selectionPID(track2)) {
          continue;
        }
        auto track2ID = track2.globalIndex();
        if (track2ID == track1ID) {
          continue;
        }
        if (!selectionPair(track1, track2)) {
          continue;
        }
        KaonPlus = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massKa);
        KaonMinus = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massKa);
        PhiMesonMother = KaonPlus + KaonMinus;
        histos.fill(HIST("hpTvsRapidity"), PhiMesonMother.Pt(), PhiMesonMother.Rapidity());
        if (TMath::Abs(PhiMesonMother.Rapidity()) > confRapidity) {
          continue;
        }
        ROOT::Math::Boost boost{PhiMesonMother.BoostToCM()};
        fourVecDauCM = boost(KaonMinus);
        threeVecDauCM = fourVecDauCM.Vect();
        threeVecDauCMXY = ROOT::Math::XYZVector(threeVecDauCM.X(), threeVecDauCM.Y(), 0.);
        eventplaneVec = ROOT::Math::XYZVector(std::cos(2.0 * psiFT0C), std::sin(2.0 * psiFT0C), 0);
        eventplaneVecNorm = ROOT::Math::XYZVector(std::sin(2.0 * psiFT0C), -std::cos(2.0 * psiFT0C), 0);

        // auto cosinephidaughterstarminuspsi = eventplaneVec.Dot(threeVecDauCMXY) / std::sqrt(threeVecDauCMXY.Mag2()) / std::sqrt(eventplaneVec.Mag2());
        // auto SA = (2.0 * cosinephidaughterstarminuspsi * cosinephidaughterstarminuspsi) - 1.0;
        auto cosPhistarminuspsi = GetPhiInRange(fourVecDauCM.Phi() - psiFT0C);
        auto SA = TMath::Cos(2.0 * cosPhistarminuspsi);
        auto cosThetaStarOP = TMath::Abs(eventplaneVecNorm.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(eventplaneVecNorm.Mag2()));
        auto cosThetaStarIP = TMath::Abs(eventplaneVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(eventplaneVec.Mag2()));
        auto phiminuspsi = GetPhiInRange(PhiMesonMother.Phi() - psiFT0C);
        auto v2 = TMath::Cos(2.0 * phiminuspsi);
        histos.fill(HIST("hSparseV2SASameEvent_costhetastarOP"), PhiMesonMother.M(), PhiMesonMother.Pt(), cosThetaStarOP, phiminuspsi, centrality);
        histos.fill(HIST("hSparseV2SASameEvent_costhetastarIP"), PhiMesonMother.M(), PhiMesonMother.Pt(), cosThetaStarIP, phiminuspsi, centrality);
        histos.fill(HIST("hSparseV2SASameEvent_SA"), PhiMesonMother.M(), PhiMesonMother.Pt(), SA, phiminuspsi, centrality);
        histos.fill(HIST("hSparseV2SASameEvent_V2"), PhiMesonMother.M(), PhiMesonMother.Pt(), v2, centrality);
      }
    }
  }
  PROCESS_SWITCH(phipbpb, processSameEvent, "Process Same event", true);
  void processMixedEvent(EventCandidates const& collisions, TrackCandidates const& tracks)
  {
    BinningTypeVertexContributor binningOnPositions{{axisVertex, axisMultiplicityClass, axisEPAngle}, true};
    for (auto const& [collision1, collision2] : o2::soa::selfCombinations(binningOnPositions, cfgNoMixedEvents, -1, collisions, collisions)) {
      if (!collision1.sel8() || !collision2.sel8()) {
        continue;
      }
      if (!collision1.triggereventep() || !collision2.triggereventep()) {
        return;
      }
      if (timFrameEvsel && (!collision1.selection_bit(aod::evsel::kNoTimeFrameBorder) || !collision2.selection_bit(aod::evsel::kNoTimeFrameBorder))) {
        return;
      }
      auto posThisColl = posTracks->sliceByCached(aod::track::collisionId, collision1.globalIndex(), cache);
      auto negThisColl = negTracks->sliceByCached(aod::track::collisionId, collision2.globalIndex(), cache);
      auto centrality = collision1.centFT0C();
      auto centrality2 = collision2.centFT0C();
      auto psiFT0C = collision1.psiFT0C();

      if (additionalEvsel && !eventSelected(collision1, tracks.size(), centrality)) {
        return;
      }
      if (additionalEvsel && !eventSelected(collision2, tracks.size(), centrality2)) {
        return;
      }
      for (auto& [track1, track2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(posThisColl, negThisColl))) {
        // track selection
        if (!selectionTrack(track1) || !selectionTrack(track2)) {
          continue;
        }
        // PID check
        if ((ispTdepPID && !selectionPIDpTdependent(track1)) || (ispTdepPID && !selectionPIDpTdependent(track2))) {
          continue;
        }
        if ((!ispTdepPID && !selectionPID(track1)) || (!ispTdepPID && !selectionPID(track2))) {
          continue;
        }
        if (!selectionPair(track1, track2)) {
          continue;
        }
        KaonPlus = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massKa);
        KaonMinus = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massKa);
        PhiMesonMother = KaonPlus + KaonMinus;
        if (TMath::Abs(PhiMesonMother.Rapidity()) > confRapidity) {
          continue;
        }
        ROOT::Math::Boost boost{PhiMesonMother.BoostToCM()};
        fourVecDauCM = boost(KaonMinus);
        threeVecDauCM = fourVecDauCM.Vect();
        threeVecDauCMXY = ROOT::Math::XYZVector(threeVecDauCM.X(), threeVecDauCM.Y(), 0.);
        eventplaneVec = ROOT::Math::XYZVector(std::cos(2.0 * psiFT0C), std::sin(2.0 * psiFT0C), 0);
        eventplaneVecNorm = ROOT::Math::XYZVector(std::sin(2.0 * psiFT0C), -std::cos(2.0 * psiFT0C), 0);

        // auto cosinephidaughterstarminuspsi = eventplaneVec.Dot(threeVecDauCMXY) / std::sqrt(threeVecDauCMXY.Mag2()) / std::sqrt(eventplaneVec.Mag2());
        // auto SA = (2.0 * cosinephidaughterstarminuspsi * cosinephidaughterstarminuspsi) - 1.0;
        auto cosPhistarminuspsi = GetPhiInRange(fourVecDauCM.Phi() - psiFT0C);
        auto SA = TMath::Cos(2.0 * cosPhistarminuspsi);
        auto cosThetaStarOP = TMath::Abs(eventplaneVecNorm.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(eventplaneVecNorm.Mag2()));
        auto cosThetaStarIP = TMath::Abs(eventplaneVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(eventplaneVec.Mag2()));
        auto phiminuspsi = GetPhiInRange(PhiMesonMother.Phi() - psiFT0C);
        auto v2 = TMath::Cos(2.0 * phiminuspsi);
        histos.fill(HIST("hSparseV2SAMixedEvent_costhetastarOP"), PhiMesonMother.M(), PhiMesonMother.Pt(), cosThetaStarOP, phiminuspsi, centrality);
        histos.fill(HIST("hSparseV2SAMixedEvent_costhetastarIP"), PhiMesonMother.M(), PhiMesonMother.Pt(), cosThetaStarIP, phiminuspsi, centrality);
        histos.fill(HIST("hSparseV2SAMixedEvent_SA"), PhiMesonMother.M(), PhiMesonMother.Pt(), SA, phiminuspsi, centrality);
        histos.fill(HIST("hSparseV2SAMixedEvent_V2"), PhiMesonMother.M(), PhiMesonMother.Pt(), v2, centrality);
      }
    }
  }
  PROCESS_SWITCH(phipbpb, processMixedEvent, "Process Mixed event", true);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<phipbpb>(cfgc, TaskName{"phipbpb"})};
}
