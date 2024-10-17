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
// Lambda polarisation task
// prottay.das@cern.ch

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

#include "PWGLF/DataModel/SPCalibrationTables.h"
// #include "SPCalibrationTableswrite.h"
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
#include "Common/DataModel/FT0Corrected.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

struct lambdapolsp {

  int mRunNumber;
  int multEstimator;
  float d_bz;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Service<o2::framework::O2DatabasePDG> pdg;

  // fill output
  Configurable<bool> additionalEvSel{"additionalEvSel", false, "additionalEvSel"};
  Configurable<bool> additionalEvSel2{"additionalEvSel2", false, "additionalEvSel2"};
  Configurable<bool> additionalEvSel3{"additionalEvSel3", false, "additionalEvSel3"};
  Configurable<bool> correction{"correction", false, "fill histograms including corrections"};
  // events
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCutCentralityMax{"cfgCutCentralityMax", 50.0f, "Accepted maximum Centrality"};
  Configurable<float> cfgCutCentralityMin{"cfgCutCentralityMin", 30.0f, "Accepted minimum Centrality"};
  // proton track cut
  Configurable<float> confRapidity{"confRapidity", 0.8, "cut on Rapidity"};
  Configurable<float> cfgCutPT{"cfgCutPT", 0.15, "PT cut on daughter track"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8, "Eta cut on daughter track"};
  Configurable<float> cfgCutDCAxy{"cfgCutDCAxy", 0.1f, "DCAxy range for tracks"};
  Configurable<float> cfgCutDCAz{"cfgCutDCAz", 0.1f, "DCAz range for tracks"};
  Configurable<int> cfgITScluster{"cfgITScluster", 5, "Number of ITS cluster"};
  Configurable<int> cfgTPCcluster{"cfgTPCcluster", 70, "Number of TPC cluster"};
  Configurable<bool> isPVContributor{"isPVContributor", true, "is PV contributor"};
  Configurable<bool> checkwithpub{"checkwithpub", true, "checking results with published"};
  // Configs for V0
  Configurable<float> ConfV0PtMin{"ConfV0PtMin", 0.f, "Minimum transverse momentum of V0"};
  Configurable<float> ConfV0Rap{"ConfV0Rap", 0.8f, "Rapidity range of V0"};
  Configurable<double> ConfV0DCADaughMax{"ConfV0DCADaughMax", 0.2f, "Maximum DCA between the V0 daughters"};
  Configurable<double> ConfV0CPAMin{"ConfV0CPAMin", 0.9998f, "Minimum CPA of V0"};
  Configurable<float> ConfV0TranRadV0Min{"ConfV0TranRadV0Min", 1.5f, "Minimum transverse radius"};
  Configurable<float> ConfV0TranRadV0Max{"ConfV0TranRadV0Max", 100.f, "Maximum transverse radius"};
  Configurable<double> cMinV0DCA{"cMinV0DCA", 0.05, "Minimum V0 daughters DCA to PV"};
  Configurable<float> cMaxV0LifeTime{"cMaxV0LifeTime", 20, "Maximum V0 life time"};
  Configurable<float> cSigmaMassKs0{"cSigmaMassKs0", 0.006, "Sigma cut on KS0 mass"};
  Configurable<float> cMinLambdaMass{"cMinLambdaMass", 1.0, "Minimum lambda mass"};
  Configurable<float> cMaxLambdaMass{"cMaxLambdaMass", 1.2, "Maximum lambda mass"};
  // config for V0 daughters
  Configurable<float> ConfDaughEta{"ConfDaughEta", 0.8f, "V0 Daugh sel: max eta"};
  Configurable<float> ConfDaughPt{"ConfDaughPt", 0.1f, "V0 Daugh sel: min pt"};
  Configurable<float> ConfDaughTPCnclsMin{"ConfDaughTPCnclsMin", 50.f, "V0 Daugh sel: Min. nCls TPC"};
  Configurable<double> ConfDaughDCAMin{"ConfDaughDCAMin", 0.08f, "V0 Daugh sel:  Max. DCA Daugh to PV (cm)"};
  Configurable<float> ConfDaughPIDCuts{"ConfDaughPIDCuts", 3, "PID selections for KS0 daughters"};

  Configurable<int> CentNbins{"CentNbins", 16, "Number of bins in cent histograms"};
  Configurable<float> lbinCent{"lbinCent", 0.0, "lower bin value in cent histograms"};
  Configurable<float> hbinCent{"hbinCent", 80.0, "higher bin value in cent histograms"};
  Configurable<int> SANbins{"SANbins", 20, "Number of bins in costhetastar"};
  Configurable<float> lbinSA{"lbinSA", -1.0, "lower bin value in costhetastar histograms"};
  Configurable<float> hbinSA{"hbinSA", 1.0, "higher bin value in costhetastar histograms"};
  Configurable<int> PolNbins{"PolNbins", 20, "Number of bins in polarisation"};
  Configurable<float> lbinPol{"lbinPol", -1.0, "lower bin value in #phi-#psi histograms"};
  Configurable<float> hbinPol{"hbinPol", 1.0, "higher bin value in #phi-#psi histograms"};
  Configurable<int> IMNbins{"IMNbins", 100, "Number of bins in invariant mass"};
  Configurable<float> lbinIM{"lbinIM", 1.0, "lower bin value in IM histograms"};
  Configurable<float> hbinIM{"hbinIM", 1.2, "higher bin value in IM histograms"};
  Configurable<int> ptNbins{"ptNbins", 50, "Number of bins in pt"};
  Configurable<float> lbinpt{"lbinpt", 0.0, "lower bin value in pt histograms"};
  Configurable<float> hbinpt{"hbinpt", 10.0, "higher bin value in pt histograms"};

  SliceCache cache;
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext&)
  {
    AxisSpec thnAxispT{ptNbins, lbinpt, hbinpt, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec thnAxisInvMass{IMNbins, lbinIM, hbinIM, "#it{M} (GeV/#it{c}^{2})"};
    AxisSpec thnAxisPol{PolNbins, lbinPol, hbinPol, "Sin(#phi - #psi)"};
    AxisSpec thnAxisCosThetaStar{SANbins, lbinSA, hbinSA, "SA"};
    AxisSpec centAxis = {CentNbins, lbinCent, hbinCent, "V0M (%)"};
    AxisSpec etaAxis = {8, -0.8, 0.8, "Eta"};

    if (checkwithpub) {
      histos.add("hpuxQxpvscent", "hpuxQxpvscent", kTProfile, {centAxis});
      histos.add("hpuyQypvscent", "hpuyQypvscent", kTProfile, {centAxis});
      histos.add("hpuxQxtvscent", "hpuxQxtvscent", kTProfile, {centAxis});
      histos.add("hpuyQytvscent", "hpuyQytvscent", kTProfile, {centAxis});
      histos.add("hpQxtQxpvscent", "hpQxtQxpvscent", kTProfile, {centAxis});
      histos.add("hpQytQypvscent", "hpQytQypvscent", kTProfile, {centAxis});

      histos.add("hpposuxyQxytvseta", "hpposuxyQxytvseta", kTProfile, {etaAxis});
      histos.add("hpposuxyQxypvseta", "hpposuxyQxypvseta", kTProfile, {etaAxis});
      histos.add("hpposQxytpvseta", "hpposQxytpvseta", kTProfile, {etaAxis});
      histos.add("hpneguxyQxytvseta", "hpneguxyQxytvseta", kTProfile, {etaAxis});
      histos.add("hpneguxyQxypvseta", "hpneguxyQxypvseta", kTProfile, {etaAxis});
      histos.add("hpnegQxytpvseta", "hpnegQxytpvseta", kTProfile, {etaAxis});
    }

    histos.add("hCentrality", "Centrality distribution", kTH1F, {{centAxis}});
    histos.add("hVtxZ", "Vertex distribution in Z;Z (cm)", kTH1F, {{20, -10.0, 10.0}});
    histos.add("hpRes", "hpRes", kTProfile, {centAxis});
    histos.add("hpResSin", "hpResSin", kTProfile, {centAxis});
    histos.add("hpCosPsiA", "hpCosPsiA", kTProfile, {centAxis});
    histos.add("hpCosPsiC", "hpCosPsiC", kTProfile, {centAxis});
    histos.add("hpSinPsiA", "hpSinPsiA", kTProfile, {centAxis});
    histos.add("hpSinPsiC", "hpSinPsiC", kTProfile, {centAxis});

    histos.add("hSparseLambdaPolA", "hSparseLambdaPolA", HistType::kTHnSparseF, {thnAxisInvMass, thnAxispT, thnAxisPol, centAxis}, true);
    histos.add("hSparseLambdaPolC", "hSparseLambdaPolC", HistType::kTHnSparseF, {thnAxisInvMass, thnAxispT, thnAxisPol, centAxis}, true);
    histos.add("hSparseAntiLambdaPolA", "hSparseAntiLambdaPolA", HistType::kTHnSparseF, {thnAxisInvMass, thnAxispT, thnAxisPol, centAxis}, true);
    histos.add("hSparseAntiLambdaPolC", "hSparseAntiLambdaPolC", HistType::kTHnSparseF, {thnAxisInvMass, thnAxispT, thnAxisPol, centAxis}, true);
    if (correction) {
      histos.add("hSparseLambdaPolA_corr", "hSparseLambdaPolA_corr", HistType::kTHnSparseF, {thnAxisInvMass, thnAxispT, thnAxisCosThetaStar, thnAxisPol, thnAxisPol, thnAxisPol, centAxis}, true);
      histos.add("hSparseLambdaPolC_corr", "hSparseLambdaPolC_corr", HistType::kTHnSparseF, {thnAxisInvMass, thnAxispT, thnAxisCosThetaStar, thnAxisPol, thnAxisPol, thnAxisPol, centAxis}, true);
      histos.add("hSparseAntiLambdaPolA_corr", "hSparseAntiLambdaPolA_corr", HistType::kTHnSparseF, {thnAxisInvMass, thnAxispT, thnAxisCosThetaStar, thnAxisPol, thnAxisPol, thnAxisPol, centAxis}, true);
      histos.add("hSparseAntiLambdaPolC_corr", "hSparseAntiLambdaPolC_corr", HistType::kTHnSparseF, {thnAxisInvMass, thnAxispT, thnAxisCosThetaStar, thnAxisPol, thnAxisPol, thnAxisPol, centAxis}, true);
    }
  }

  template <typename T>
  bool selectionTrack(const T& candidate)
  {
    if (!isPVContributor || !candidate.isGlobalTrackWoDCA() || !(candidate.itsNCls() > cfgITScluster) || !(candidate.tpcNClsFound() > cfgTPCcluster)) {
      return false;
    }
    return true;
  }

  template <typename Collision, typename V0>
  bool SelectionV0(Collision const& collision, V0 const& candidate)
  {

    const float pT = candidate.pt();
    // const std::vector<float> decVtx = {candidate.x(), candidate.y(), candidate.z()};
    const float tranRad = candidate.v0radius();
    const float dcaDaughv0 = TMath::Abs(candidate.dcaV0daughters());
    const float cpav0 = candidate.v0cosPA();

    float CtauLambda = candidate.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * massLambda;
    // float lowmasscutlambda = cMinLambdaMass;
    // float highmasscutlambda = cMaxLambdaMass;

    if (pT < ConfV0PtMin) {
      return false;
    }
    if (TMath::Abs(candidate.dcapostopv()) < cMinV0DCA) {
      return false;
    }
    if (TMath::Abs(candidate.dcanegtopv()) < cMinV0DCA) {
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
    if (TMath::Abs(CtauLambda) > cMaxV0LifeTime) {
      return false;
    }
    if (TMath::Abs(candidate.yLambda()) > ConfV0Rap) {
      return false;
    }
    return true;
  }
  template <typename T>
  bool isSelectedV0Daughter(T const& track, int pid)
  {
    const auto eta = track.eta();
    const auto pt = track.pt();
    const auto tpcNClsF = track.tpcNClsFound();
    // const auto dcaXY = track.dcaXY();
    // const auto sign = track.sign();
    /*
    if (charge < 0 && sign > 0) {
      return false;
    }
    if (charge > 0 && sign < 0) {
      return false;
      }*/
    if (TMath::Abs(eta) > ConfDaughEta) {
      return false;
    }
    if (pt < ConfDaughPt) {
      return false;
    }
    if (tpcNClsF < ConfDaughTPCnclsMin) {
      return false;
    }
    /*
    if (track.tpcCrossedRowsOverFindableCls() < 0.8) {
      return false;
      }
    if (TMath::Abs(dcaXY) < ConfDaughDCAMin) {
      return false;
      }*/
    if (pid == 0 && TMath::Abs(track.tpcNSigmaPr()) > ConfDaughPIDCuts) {
      return false;
    }
    if (pid == 1 && TMath::Abs(track.tpcNSigmaPi()) > ConfDaughPIDCuts) {
      return false;
    }

    return true;
  }

  double GetPhiInRange(double phi)
  {
    double result = phi;
    while (result < 0) {
      result = result + 2. * TMath::Pi();
    }
    while (result > 2. * TMath::Pi()) {
      result = result - 2. * TMath::Pi();
    }
    return result;
  }

  ROOT::Math::PxPyPzMVector Lambda, Proton, Pion, fourVecDauCM;
  // ROOT::Math::XYZVector threeVecDauCM, threeVecDauCMXY, eventplaneVec, eventplaneVecNorm, beamvector;
  ROOT::Math::XYZVector threeVecDauCM, threeVecDauCMXY;
  float phiangle = 0.0;
  // double massPi = TDatabasePDG::Instance()->GetParticle(kPiPlus)->Mass();
  // double massPr = TDatabasePDG::Instance()->GetParticle(kProton)->Mass();
  // double massLambda = TDatabasePDG::Instance()->GetParticle(kLambda0)->Mass();
  double massLambda = o2::constants::physics::MassLambda;
  double massPr = o2::constants::physics::MassProton;
  double massPi = o2::constants::physics::MassPionCharged;

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter centralityFilter = (nabs(aod::cent::centFT0C) < cfgCutCentralityMax && nabs(aod::cent::centFT0C) > cfgCutCentralityMin);
  Filter acceptanceFilter = (nabs(aod::track::eta) < cfgCutEta && nabs(aod::track::pt) > cfgCutPT);
  Filter dcaCutFilter = (nabs(aod::track::dcaXY) < cfgCutDCAxy) && (nabs(aod::track::dcaZ) < cfgCutDCAz);

  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::FV0Mults, aod::TPCMults, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0As, aod::SPCalibrationTables, aod::Mults>>;
  // using AllTrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPi, aod::pidTPCFullPr>;
  using AllTrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPi, aod::pidTPCFullPr>>;
  using ResoV0s = aod::V0Datas;

  // void processData(EventCandidates::iterator const& collision, AllTrackCandidates const&, ResoV0s const& V0s, aod::BCs const&)
  void processData(EventCandidates::iterator const& collision, AllTrackCandidates const& tracks, ResoV0s const& V0s, aod::BCs const&)
  {

    if (!collision.sel8()) {
      return;
    }
    auto centrality = collision.centFT0C();

    if (!collision.triggerevent()) {
      return;
    }

    if (additionalEvSel && (!collision.selection_bit(aod::evsel::kNoSameBunchPileup) || !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV))) {
      return;
    }
    if (additionalEvSel2 && (!collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard))) {
      return;
    }

    if (additionalEvSel3 && (!collision.selection_bit(aod::evsel::kNoTimeFrameBorder) || !collision.selection_bit(aod::evsel::kNoITSROFrameBorder))) {
      return;
    }

    auto psiZDCC = collision.psiZDCC();
    auto psiZDCA = collision.psiZDCA();

    histos.fill(HIST("hCentrality"), centrality);
    histos.fill(HIST("hVtxZ"), collision.posZ());
    histos.fill(HIST("hpRes"), centrality, (TMath::Cos(psiZDCA - psiZDCC)));
    histos.fill(HIST("hpResSin"), centrality, (TMath::Sin(psiZDCA - psiZDCC)));
    histos.fill(HIST("hpCosPsiA"), centrality, (TMath::Cos(psiZDCA)));
    histos.fill(HIST("hpCosPsiC"), centrality, (TMath::Cos(psiZDCC)));
    histos.fill(HIST("hpSinPsiA"), centrality, (TMath::Sin(psiZDCA)));
    histos.fill(HIST("hpSinPsiC"), centrality, (TMath::Sin(psiZDCC)));

    ///////////checking v1 and v2////////////////////////////////

    if (checkwithpub) {
      auto qxZDCA = collision.qxZDCA();
      auto qxZDCC = collision.qxZDCC();
      auto qyZDCA = collision.qyZDCA();
      auto qyZDCC = collision.qyZDCC();

      for (auto track : tracks) {
        if (!selectionTrack(track)) {
          continue;
        }

        float sign = track.sign();
        if (sign == 0.0) // removing neutral particles
          continue;

        auto ux = GetPhiInRange(TMath::Cos(track.phi()));
        auto uy = GetPhiInRange(TMath::Sin(track.phi()));

        auto uxQxp = ux * qxZDCA;
        auto uyQyp = uy * qyZDCA;
        auto uxyQxyp = uxQxp + uyQyp;
        auto uxQxt = ux * qxZDCC;
        auto uyQyt = uy * qyZDCC;
        auto uxyQxyt = uxQxt + uyQyt;

        auto QxtQxp = qxZDCA * qxZDCC;
        auto QytQyp = qyZDCA * qyZDCC;
        auto Qxytp = QxtQxp + QytQyp;
        histos.fill(HIST("hpuxQxpvscent"), centrality, uxQxp);
        histos.fill(HIST("hpuyQypvscent"), centrality, uyQyp);
        histos.fill(HIST("hpuxQxtvscent"), centrality, uxQxt);
        histos.fill(HIST("hpuyQytvscent"), centrality, uyQyt);
        histos.fill(HIST("hpQxtQxpvscent"), centrality, QxtQxp);
        histos.fill(HIST("hpQytQypvscent"), centrality, QytQyp);

        if (centrality > 5.0 && centrality < 40.0) {
          if (track.pt() > 0.2) {
            if (sign > 0.0) {
              histos.fill(HIST("hpposuxyQxytvseta"), track.eta(), uxyQxyt);
              histos.fill(HIST("hpposuxyQxypvseta"), track.eta(), uxyQxyp);
              histos.fill(HIST("hpposQxytpvseta"), track.eta(), Qxytp);
            } else if (sign < 0.0) {
              histos.fill(HIST("hpneguxyQxytvseta"), track.eta(), uxyQxyt);
              histos.fill(HIST("hpneguxyQxypvseta"), track.eta(), uxyQxyp);
              histos.fill(HIST("hpnegQxytpvseta"), track.eta(), Qxytp);
            }
          }
        }
      }
    }

    //////////////////////////////////////////////////////////////////////////

    for (auto v0 : V0s) {

      auto postrack = v0.template posTrack_as<AllTrackCandidates>();
      auto negtrack = v0.template negTrack_as<AllTrackCandidates>();

      int LambdaTag = 0;
      int aLambdaTag = 0;

      if (isSelectedV0Daughter(postrack, 0) && isSelectedV0Daughter(negtrack, 1)) {
        LambdaTag = 1;
      }
      if (isSelectedV0Daughter(negtrack, 0) && isSelectedV0Daughter(postrack, 1)) {
        aLambdaTag = 1;
      }

      if (LambdaTag == aLambdaTag)
        continue;

      if (!SelectionV0(collision, v0)) {
        continue;
      }

      if (LambdaTag) {
        Proton = ROOT::Math::PxPyPzMVector(postrack.px(), postrack.py(), postrack.pz(), massPr);
        Pion = ROOT::Math::PxPyPzMVector(negtrack.px(), negtrack.py(), negtrack.pz(), massPi);
      }
      if (aLambdaTag) {
        Proton = ROOT::Math::PxPyPzMVector(negtrack.px(), negtrack.py(), negtrack.pz(), massPr);
        Pion = ROOT::Math::PxPyPzMVector(postrack.px(), postrack.py(), postrack.pz(), massPi);
      }
      Lambda = Proton + Pion;
      Lambda.SetM(massLambda);

      ROOT::Math::Boost boost{Lambda.BoostToCM()};
      fourVecDauCM = boost(Proton);
      threeVecDauCM = fourVecDauCM.Vect();
      // beamvector = ROOT::Math::XYZVector(0, 0, 1);
      // eventplaneVec = ROOT::Math::XYZVector(collision.qFT0C(), collision.qFT0A(), 0); //this needs to be changed
      // eventplaneVecNorm = eventplaneVec.Cross(beamvector); //z'
      phiangle = TMath::ATan2(fourVecDauCM.Py(), fourVecDauCM.Px());

      auto phiminuspsiC = GetPhiInRange(phiangle - psiZDCC);
      auto phiminuspsiA = GetPhiInRange(phiangle - psiZDCA);
      // auto cosThetaStar = eventplaneVecNorm.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(eventplaneVecNorm.Mag2());
      auto cosThetaStar = fourVecDauCM.Pz() / fourVecDauCM.P(); // A0 correction
      auto PolC = TMath::Sin(phiminuspsiC);
      auto PolA = TMath::Sin(phiminuspsiA);

      // needed for corrections
      auto sinPhiStar = TMath::Sin(GetPhiInRange(phiangle));
      auto cosPhiStar = TMath::Cos(GetPhiInRange(phiangle));
      // auto sinThetaStarcosphiphiStar=sinThetaStar*TMath::Cos(2* GetPhiInRange((Lambda.Phi()-phiangle))); //A2 correction

      if (LambdaTag) {
        if (correction) {
          histos.fill(HIST("hSparseLambdaPolA_corr"), v0.mLambda(), v0.pt(), cosThetaStar, sinPhiStar, cosPhiStar, PolA, centrality);
          histos.fill(HIST("hSparseLambdaPolC_corr"), v0.mLambda(), v0.pt(), cosThetaStar, sinPhiStar, cosPhiStar, PolC, centrality);
        } else {
          histos.fill(HIST("hSparseLambdaPolA"), v0.mLambda(), v0.pt(), PolA, centrality);
          histos.fill(HIST("hSparseLambdaPolC"), v0.mLambda(), v0.pt(), PolC, centrality);
        }
      }
      if (aLambdaTag) {
        if (correction) {
          histos.fill(HIST("hSparseAntiLambdaPolA_corr"), v0.mAntiLambda(), v0.pt(), cosThetaStar, sinPhiStar, cosPhiStar, PolA, centrality);
          histos.fill(HIST("hSparseAntiLambdaPolC_corr"), v0.mAntiLambda(), v0.pt(), cosThetaStar, sinPhiStar, cosPhiStar, PolC, centrality);
        } else {
          histos.fill(HIST("hSparseAntiLambdaPolA"), v0.mAntiLambda(), v0.pt(), PolA, centrality);
          histos.fill(HIST("hSparseAntiLambdaPolC"), v0.mAntiLambda(), v0.pt(), PolC, centrality);
        }
      }
    }
  }
  PROCESS_SWITCH(lambdapolsp, processData, "Process data", true);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<lambdapolsp>(cfgc, TaskName{"lambdapolsp"})};
}
