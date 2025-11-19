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
// Cascade polarisation task
// prottay.das@cern.ch

#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/SPCalibrationTables.h"
#include "PWGLF/DataModel/cascqaanalysis.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/PhysicsConstants.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/StepTHn.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

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

using dauTracks = soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>;
using v0Candidates = soa::Join<aod::V0CollRefs, aod::V0Cores, aod::V0Extras>;

struct cascpolsp {

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Service<o2::framework::O2DatabasePDG> pdg;

  // fill output
  Configurable<bool> additionalEvSel{"additionalEvSel", false, "additionalEvSel"};
  Configurable<bool> additionalEvSel3{"additionalEvSel3", false, "additionalEvSel3"};
  Configurable<int> QxyNbins{"QxyNbins", 100, "Number of bins in QxQy histograms"};
  Configurable<float> lbinQxy{"lbinQxy", -5.0, "lower bin value in QxQy histograms"};
  Configurable<float> hbinQxy{"hbinQxy", 5.0, "higher bin value in QxQy histograms"};

  // events
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCutCentralityMax{"cfgCutCentralityMax", 50.0f, "Accepted maximum Centrality"};
  Configurable<float> cfgCutCentralityMin{"cfgCutCentralityMin", 30.0f, "Accepted minimum Centrality"};
  // track cut
  Configurable<float> cfgCutPT{"cfgCutPT", 0.15, "PT cut on daughter track"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8, "Eta cut on daughter track"};

  // Configs for V0
  Configurable<float> ConfV0PtMin{"ConfV0PtMin", 0.f, "Minimum transverse momentum of V0"};
  Configurable<float> ConfV0Rap{"ConfV0Rap", 0.8f, "Rapidity range of V0"};
  Configurable<double> ConfV0DCADaughMax{"ConfV0DCADaughMax", 0.2f, "Maximum DCA between the V0 daughters"};
  Configurable<double> ConfV0CPAMin{"ConfV0CPAMin", 0.9998f, "Minimum CPA of V0"};
  Configurable<float> ConfV0TranRadV0Min{"ConfV0TranRadV0Min", 1.5f, "Minimum transverse radius"};
  Configurable<float> ConfV0TranRadV0Max{"ConfV0TranRadV0Max", 100.f, "Maximum transverse radius"};
  Configurable<double> cMaxV0DCA{"cMaxV0DCA", 1.2, "Maximum V0 DCA to PV"};
  Configurable<double> cMinV0DCA{"cMinV0DCA", 0.05, "Minimum V0 daughters DCA to PV"};
  Configurable<float> cMaxV0LifeTime{"cMaxV0LifeTime", 20, "Maximum V0 life time"};

  // config for V0 daughters
  Configurable<float> ConfDaughEta{"ConfDaughEta", 0.8f, "V0 Daugh sel: max eta"};
  Configurable<float> cfgDaughPrPt{"cfgDaughPrPt", 0.4, "minimum daughter proton pt"};
  Configurable<float> cfgDaughPiPt{"cfgDaughPiPt", 0.2, "minimum daughter pion pt"};
  Configurable<float> ConfDaughTPCnclsMin{"ConfDaughTPCnclsMin", 50.f, "V0 Daugh sel: Min. nCls TPC"};
  Configurable<double> ConfDaughDCAMin{"ConfDaughDCAMin", 0.08f, "V0 Daugh sel:  Max. DCA Daugh to PV (cm)"};
  Configurable<float> ConfDaughPIDCuts{"ConfDaughPIDCuts", 3, "PID selections for Lambda daughters"};

  // config for cascades
  Configurable<float> cfgcasc_radius{"cfgcasc_radius", 0.8f, "Cascade radius"};
  Configurable<float> cfgcasc_casccospa{"cfgcasc_casccospa", 0.99f, "Cascade cosPA"};
  Configurable<float> cfgcasc_v0cospa{"cfgcasc_v0cospa", 0.99f, "Cascade v0cosPA"};
  Configurable<float> cfgcasc_dcav0topv{"cfgcasc_dcav0topv", 0.3f, "Cascade dcav0"};
  Configurable<float> cfgcasc_dcabachtopv{"cfgcasc_dcabachtopv", 0.3f, "Cascade dcabach"};
  Configurable<float> cfgcasc_dcacascdau{"cfgcasc_dcacascdau", 0.3f, "Cascade dcadau"};
  Configurable<float> cfgcasc_dcav0dau{"cfgcasc_dcav0dau", 0.3f, "Cascade dcav0dau"};
  Configurable<float> cfgcasc_mlambdawindow{"cfgcasc_mlambdawindow", 0.43f, "Cascade masswindow"};
  Configurable<float> cfgcascv0_radius{"cfgcascv0_radius", 0.8f, "Cascade v0 radius"};
  Configurable<float> cfgcasc_dcapostopv{"cfgcasc_dcapostopv", 0.3f, "Cascade dcapostoPV"};
  Configurable<float> cfgcasc_dcanegtopv{"cfgcasc_dcanegtopv", 0.3f, "Cascade dcanegtoPV"};
  Configurable<float> cfgcasc_lowmass{"cfgcasc_lowmass", 1.311, "Cascade lowmass cut"};
  Configurable<float> cfgcasc_highmass{"cfgcasc_highmass", 1.331, "Cascade highmass cut"};

  // config for histograms
  Configurable<int> IMNbins{"IMNbins", 100, "Number of bins in invariant mass"};
  Configurable<float> lbinIM{"lbinIM", 1.0, "lower bin value in IM histograms"};
  Configurable<float> hbinIM{"hbinIM", 1.2, "higher bin value in IM histograms"};
  Configurable<int> resNbins{"resNbins", 50, "Number of bins in reso"};
  Configurable<float> lbinres{"lbinres", 0.0, "lower bin value in reso histograms"};
  Configurable<float> hbinres{"hbinres", 10.0, "higher bin value in reso histograms"};

  ConfigurableAxis configcentAxis{"configcentAxis", {VARIABLE_WIDTH, 0.0, 10.0, 40.0, 80.0}, "Cent V0M"};
  ConfigurableAxis configthnAxispT{"configthnAxisPt", {VARIABLE_WIDTH, 0.2, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.5, 8.0, 10.0, 100.0}, "#it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis configetaAxis{"configetaAxis", {VARIABLE_WIDTH, -0.8, -0.4, -0.2, 0, 0.2, 0.4, 0.8}, "Eta"};
  ConfigurableAxis configthnAxisPol{"configthnAxisPol", {VARIABLE_WIDTH, -1.0, -0.6, -0.2, 0, 0.2, 0.4, 0.8}, "Pol"};
  ConfigurableAxis configphiAxis{"configphiAxis", {VARIABLE_WIDTH, 0.0, 0.2, 0.4, 0.8, 1.0, 2.0, 2.5, 3.0, 4.0, 5.0, 5.5, 6.28}, "PhiAxis"};

  SliceCache cache;
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext&)
  {
    AxisSpec thnAxisres{resNbins, lbinres, hbinres, "Reso"};
    AxisSpec thnAxisInvMass{IMNbins, lbinIM, hbinIM, "#it{M} (GeV/#it{c}^{2})"};

    histos.add("hCentrality", "Centrality distribution", kTH1F, {{configcentAxis}});
    histos.add("hpRes", "hpRes", HistType::kTHnSparseF, {configcentAxis, thnAxisres});

    histos.add("hSparseLambdaPolA", "hSparseLambdaPolA", HistType::kTHnSparseF, {thnAxisInvMass, configthnAxispT, configetaAxis, configthnAxisPol, configcentAxis}, true);
    histos.add("hSparseLambdaPolC", "hSparseLambdaPolC", HistType::kTHnSparseF, {thnAxisInvMass, configthnAxispT, configetaAxis, configthnAxisPol, configcentAxis}, true);
    histos.add("hSparseAntiLambdaPolA", "hSparseAntiLambdaPolA", HistType::kTHnSparseF, {thnAxisInvMass, configthnAxispT, configetaAxis, configthnAxisPol, configcentAxis}, true);
    histos.add("hSparseAntiLambdaPolC", "hSparseAntiLambdaPolC", HistType::kTHnSparseF, {thnAxisInvMass, configthnAxispT, configetaAxis, configthnAxisPol, configcentAxis}, true);

    histos.add("hSparseLambda_corr1a", "hSparseLambda_corr1a", HistType::kTHnSparseF, {thnAxisInvMass, configthnAxispT, configetaAxis, configthnAxisPol, configcentAxis}, true);
    histos.add("hSparseLambda_corr1b", "hSparseLambda_corr1b", HistType::kTHnSparseF, {thnAxisInvMass, configthnAxispT, configetaAxis, configthnAxisPol, configcentAxis}, true);
    histos.add("hSparseLambda_corr1c", "hSparseLambda_corr1c", HistType::kTHnSparseF, {thnAxisInvMass, configthnAxispT, configetaAxis, configphiAxis, configcentAxis}, true);
    histos.add("hSparseAntiLambda_corr1a", "hSparseAntiLambda_corr1a", HistType::kTHnSparseF, {thnAxisInvMass, configthnAxispT, configetaAxis, configthnAxisPol, configcentAxis}, true);
    histos.add("hSparseAntiLambda_corr1b", "hSparseAntiLambda_corr1b", HistType::kTHnSparseF, {thnAxisInvMass, configthnAxispT, configetaAxis, configthnAxisPol, configcentAxis}, true);
    histos.add("hSparseAntiLambda_corr1c", "hSparseAntiLambda_corr1c", HistType::kTHnSparseF, {thnAxisInvMass, configthnAxispT, configetaAxis, configphiAxis, configcentAxis}, true);

    histos.add("hSparseLambda_corr2a", "hSparseLambda_corr2a", HistType::kTHnSparseF, {thnAxisInvMass, configthnAxispT, configetaAxis, configthnAxisPol, configcentAxis}, true);
    histos.add("hSparseLambda_corr2b", "hSparseLambda_corr2b", HistType::kTHnSparseF, {thnAxisInvMass, configthnAxispT, configetaAxis, configthnAxisPol, configcentAxis}, true);
    histos.add("hSparseAntiLambda_corr2a", "hSparseAntiLambda_corr2a", HistType::kTHnSparseF, {thnAxisInvMass, configthnAxispT, configetaAxis, configthnAxisPol, configcentAxis}, true);
    histos.add("hSparseAntiLambda_corr2b", "hSparseAntiLambda_corr2b", HistType::kTHnSparseF, {thnAxisInvMass, configthnAxispT, configetaAxis, configthnAxisPol, configcentAxis}, true);
  }

  template <typename Collision, typename V0>
  bool SelectionV0(Collision const& collision, V0 const& candidate)
  {
    if (TMath::Abs(candidate.dcav0topv()) > cMaxV0DCA) {
      return false;
    }
    const float pT = candidate.pt();
    const float tranRad = candidate.v0radius();
    const float dcaDaughv0 = TMath::Abs(candidate.dcaV0daughters());
    const float cpav0 = candidate.v0cosPA();

    float CtauLambda = candidate.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * massLambda;

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
    if (track.tpcNClsCrossedRows() < 70) {
      return false;
    }
    if (TMath::Abs(eta) > ConfDaughEta) {
      return false;
    }
    if (tpcNClsF < ConfDaughTPCnclsMin) {
      return false;
    }
    if (track.tpcCrossedRowsOverFindableCls() < 0.8) {
      return false;
    }
    /*
      if (TMath::Abs(dcaXY) < ConfDaughDCAMin) {
      return false;
      }*/

    if (pid == 0 && TMath::Abs(track.tpcNSigmaPr()) > ConfDaughPIDCuts) {
      return false;
    }
    if (pid == 1 && TMath::Abs(track.tpcNSigmaPi()) > ConfDaughPIDCuts) {
      return false;
    }

    if (pid == 0 && pt < cfgDaughPrPt) {
      return false;
    }
    if (pid == 1 && pt < cfgDaughPiPt) {
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

  template <typename TCascade, typename collision_t>
  bool IsCascAccepted(TCascade casc, collision_t collision)
  {

    if (casc.cascradius() < cfgcasc_radius)
      return false;
    if (casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()) < cfgcasc_casccospa)
      return false;
    if (casc.dcabachtopv() < cfgcasc_dcabachtopv)
      return false;
    if (casc.dcacascdaughters() > cfgcasc_dcacascdau)
      return false;
    if (std::fabs(casc.mLambda() - o2::constants::physics::MassLambda0) > cfgcasc_mlambdawindow)
      return false;

    if (casc.v0radius() < cfgcascv0_radius)
      return false;
    if (casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) < cfgcasc_v0cospa)
      return false;
    if (TMath::Abs(casc.dcav0topv(collision.posX(), collision.posY(), collision.posZ())) < cfgcasc_dcav0topv)
      return false;
    if (TMath::Abs(casc.dcaV0daughters()) > cfgcasc_dcav0dau)
      return false;
    if (TMath::Abs(casc.dcapostopv()) < cfgcasc_dcapostopv)
      return false;
    if (TMath::Abs(casc.dcanegtopv()) < cfgcasc_dcanegtopv)
      return false;
    if (casc.mXi() < cfgcasc_lowmass || casc.mXi() > cfgcasc_highmass)
      return false;

    return true;
  }

  bool shouldReject(bool LambdaTag, bool aLambdaTag,
                    const ROOT::Math::PxPyPzMVector& Lambdadummy,
                    const ROOT::Math::PxPyPzMVector& AntiLambdadummy)
  {
    const double minMass = 1.105;
    const double maxMass = 1.125;
    return (LambdaTag && aLambdaTag &&
            (Lambdadummy.M() > minMass && Lambdadummy.M() < maxMass) &&
            (AntiLambdadummy.M() > minMass && AntiLambdadummy.M() < maxMass));
  }

  void fillHistograms(bool tag1, bool tag2, const ROOT::Math::PxPyPzMVector& particle,
                      const ROOT::Math::PxPyPzMVector& daughter,
                      double psiZDCC, double psiZDCA, double centrality,
                      double candmass, double candpt, double candeta)
  {

    ROOT::Math::Boost boost{particle.BoostToCM()};
    auto fourVecDauCM = boost(daughter);
    auto phiangle = TMath::ATan2(fourVecDauCM.Py(), fourVecDauCM.Px());

    auto phiminuspsiC = GetPhiInRange(phiangle - psiZDCC);
    auto phiminuspsiA = GetPhiInRange(phiangle - psiZDCA);
    auto cosThetaStar = fourVecDauCM.Pz() / fourVecDauCM.P();
    auto sinThetaStar = TMath::Sqrt(1 - (cosThetaStar * cosThetaStar));
    auto PolC = TMath::Sin(phiminuspsiC);
    auto PolA = TMath::Sin(phiminuspsiA);

    auto sinPhiStar = TMath::Sin(GetPhiInRange(phiangle));
    auto cosPhiStar = TMath::Cos(GetPhiInRange(phiangle));
    auto sinThetaStarcosphiphiStar = sinThetaStar * TMath::Cos(2 * GetPhiInRange(particle.Phi() - phiangle));
    auto phiphiStar = GetPhiInRange(particle.Phi() - phiangle);

    // Fill histograms using constructed names
    if (tag2) {
      histos.fill(HIST("hSparseAntiLambdaPolA"), candmass, candpt, candeta, PolA, centrality);
      histos.fill(HIST("hSparseAntiLambdaPolC"), candmass, candpt, candeta, PolC, centrality);
      histos.fill(HIST("hSparseAntiLambda_corr1a"), candmass, candpt, candeta, sinPhiStar, centrality);
      histos.fill(HIST("hSparseAntiLambda_corr1b"), candmass, candpt, candeta, cosPhiStar, centrality);
      histos.fill(HIST("hSparseAntiLambda_corr1c"), candmass, candpt, candeta, phiphiStar, centrality);
      histos.fill(HIST("hSparseAntiLambda_corr2a"), candmass, candpt, candeta, sinThetaStar, centrality);
      histos.fill(HIST("hSparseAntiLambda_corr2b"), candmass, candpt, candeta, sinThetaStarcosphiphiStar, centrality);
    }
    if (tag1) {
      histos.fill(HIST("hSparseLambdaPolA"), candmass, candpt, candeta, PolA, centrality);
      histos.fill(HIST("hSparseLambdaPolC"), candmass, candpt, candeta, PolC, centrality);
      histos.fill(HIST("hSparseLambda_corr1a"), candmass, candpt, candeta, sinPhiStar, centrality);
      histos.fill(HIST("hSparseLambda_corr1b"), candmass, candpt, candeta, cosPhiStar, centrality);
      histos.fill(HIST("hSparseLambda_corr1c"), candmass, candpt, candeta, phiphiStar, centrality);
      histos.fill(HIST("hSparseLambda_corr2a"), candmass, candpt, candeta, sinThetaStar, centrality);
      histos.fill(HIST("hSparseLambda_corr2b"), candmass, candpt, candeta, sinThetaStarcosphiphiStar, centrality);
    }
  }

  ROOT::Math::PxPyPzMVector Lambda, AntiLambda, Lambdadummy, AntiLambdadummy, Proton, Pion, AntiProton, AntiPion, fourVecDauCM, LC;
  ROOT::Math::XYZVector threeVecDauCM, threeVecDauCMXY;
  double phiangle = 0.0;
  double massLambda = o2::constants::physics::MassLambda;
  double massPr = o2::constants::physics::MassProton;
  double massPi = o2::constants::physics::MassPionCharged;

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter centralityFilter = (nabs(aod::cent::centFT0C) < cfgCutCentralityMax && nabs(aod::cent::centFT0C) > cfgCutCentralityMin);
  Filter acceptanceFilter = (nabs(aod::track::eta) < cfgCutEta && nabs(aod::track::pt) > cfgCutPT);

  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::FV0Mults, aod::TPCMults, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0As, aod::SPCalibrationTables, aod::Mults>>;
  using AllTrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPi, aod::pidTPCFullPr, aod::pidTPCFullKa>>;
  using CascCandidates = soa::Join<aod::CascCollRefs, aod::CascCores, aod::CascExtras, aod::CascBBs>;

  void processcascData(EventCandidates::iterator const& collision, aod::CascDataExt const& Cascades, AllTrackCandidates const&, aod::BCs const&)
  {

    if (!collision.sel8()) {
      return;
    }
    auto centrality = collision.centFT0C();
    if (!collision.triggereventsp()) {
      return;
    }

    if (additionalEvSel && (!collision.selection_bit(aod::evsel::kNoSameBunchPileup) || !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV))) {
      return;
    }

    if (additionalEvSel3 && (!collision.selection_bit(aod::evsel::kNoTimeFrameBorder) || !collision.selection_bit(aod::evsel::kNoITSROFrameBorder))) {
      return;
    }

    auto psiZDCC = collision.psiZDCC();
    auto psiZDCA = collision.psiZDCA();

    histos.fill(HIST("hCentrality"), centrality);
    histos.fill(HIST("hpRes"), centrality, (TMath::Cos(GetPhiInRange(psiZDCA - psiZDCC))));

    for (auto& casc : Cascades) {

      auto negtrack = casc.negTrack_as<AllTrackCandidates>();
      auto postrack = casc.posTrack_as<AllTrackCandidates>();
      auto bachtrack = casc.bachelor_as<AllTrackCandidates>();

      bool CascAcc = IsCascAccepted(casc, collision);

      if (!CascAcc)
        continue;

      int LambdaTag = 0;
      int aLambdaTag = 0;

      const auto signpos = postrack.sign();
      const auto signneg = negtrack.sign();
      const auto signbach = bachtrack.sign();

      if (signpos < 0 || signneg > 0 || signbach == 0) {
        continue;
      }

      if (isSelectedV0Daughter(postrack, 0) && isSelectedV0Daughter(negtrack, 1)) {
        LambdaTag = 1;
      }
      if (isSelectedV0Daughter(negtrack, 0) && isSelectedV0Daughter(postrack, 1)) {
        aLambdaTag = 1;
      }

      if (!isSelectedV0Daughter(bachtrack, 1)) // quality track selection for bachelor track
        continue;

      if (!LambdaTag && !aLambdaTag)
        continue;

      if (casc.sign() != signbach) // cascade sign to be equal to bachelor sign
        continue;

      if (casc.sign() > 0) {
        AntiProton = ROOT::Math::PxPyPzMVector(casc.pxneg(), casc.pyneg(), casc.pzneg(), massPr);
        Pion = ROOT::Math::PxPyPzMVector(casc.pxpos(), casc.pypos(), casc.pzpos(), massPi);
        AntiLambdadummy = AntiProton + Pion;

      } else {
        Proton = ROOT::Math::PxPyPzMVector(casc.pxpos(), casc.pypos(), casc.pzpos(), massPr);
        AntiPion = ROOT::Math::PxPyPzMVector(casc.pxneg(), casc.pyneg(), casc.pzneg(), massPi);
        Lambdadummy = Proton + AntiPion;
      }

      if (shouldReject(LambdaTag, aLambdaTag, Lambdadummy, AntiLambdadummy)) {
        continue;
      }

      int taga = LambdaTag;
      int tagb = aLambdaTag;

      if (LambdaTag) {
        Lambda = Proton + AntiPion;
        tagb = 0;
        fillHistograms(taga, tagb, Lambda, Proton, psiZDCC, psiZDCA, centrality, Lambda.M(), Lambda.Pt(), Lambda.Eta());
      }

      tagb = aLambdaTag;
      if (aLambdaTag) {
        AntiLambda = AntiProton + Pion;
        taga = 0;
        fillHistograms(taga, tagb, AntiLambda, AntiProton, psiZDCC, psiZDCA, centrality, AntiLambda.M(), AntiLambda.Pt(), AntiLambda.Eta());
      }
    }
  }
  PROCESS_SWITCH(cascpolsp, processcascData, "Process cascade data", false);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<cascpolsp>(cfgc, TaskName{"cascpolsp"})};
}
