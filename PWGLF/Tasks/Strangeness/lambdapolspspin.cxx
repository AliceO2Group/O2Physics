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
// Lambda spin spin correlation task
// prottay.das@cern.ch, sourav.kundu@cern.ch

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
#include <cmath>
#include <array>
#include <cstdlib>

#include "TRandom3.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "Math/GenVector/Boost.h"
#include "TF1.h"

// #include "Common/DataModel/Qvectors.h"
#include "PWGLF/DataModel/SPCalibrationTables.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/StepTHn.h"
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
#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "Common/DataModel/FT0Corrected.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

using dauTracks = soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>;
using v0Candidates = soa::Join<aod::V0CollRefs, aod::V0Cores, aod::V0Extras>;

struct lambdapolspspin {

  int mRunNumber;
  int multEstimator;
  float d_bz;
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  // fill output
  Configurable<bool> additionalEvSel{"additionalEvSel", false, "additionalEvSel"};
  Configurable<bool> additionalEvSel3{"additionalEvSel3", false, "additionalEvSel3"};

  // events
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCutCentralityMax{"cfgCutCentralityMax", 50.0f, "Accepted maximum Centrality"};
  Configurable<float> cfgCutCentralityMin{"cfgCutCentralityMin", 30.0f, "Accepted minimum Centrality"};

  // Configs for V0
  Configurable<float> ConfV0PtMin{"ConfV0PtMin", 0.f, "Minimum transverse momentum of V0"};
  Configurable<float> ConfV0Rap{"ConfV0Rap", 0.8f, "Rapidity range of V0"};
  Configurable<double> ConfV0DCADaughMax{"ConfV0DCADaughMax", 0.2f, "Maximum DCA between the V0 daughters"};
  Configurable<double> ConfV0CPAMin{"ConfV0CPAMin", 0.9998f, "Minimum CPA of V0"};
  Configurable<float> ConfV0TranRadV0Min{"ConfV0TranRadV0Min", 1.5f, "Minimum transverse radius"};
  Configurable<float> ConfV0TranRadV0Max{"ConfV0TranRadV0Max", 100.f, "Maximum transverse radius"};
  Configurable<double> cMaxV0DCA{"cMaxV0DCA", 1.2, "Maximum V0 DCA to PV"};
  Configurable<double> cMinV0DCAPr{"cMinV0DCAPr", 0.05, "Minimum V0 daughters DCA to PV for Pr"};
  Configurable<double> cMinV0DCAPi{"cMinV0DCAPi", 0.05, "Minimum V0 daughters DCA to PV for Pi"};
  Configurable<float> cMaxV0LifeTime{"cMaxV0LifeTime", 20, "Maximum V0 life time"};

  // config for V0 daughters
  Configurable<float> ConfDaughEta{"ConfDaughEta", 0.8f, "V0 Daugh sel: max eta"};
  Configurable<float> cfgDaughPrPt{"cfgDaughPrPt", 0.4, "minimum daughter proton pt"};
  Configurable<float> cfgDaughPiPt{"cfgDaughPiPt", 0.2, "minimum daughter pion pt"};
  Configurable<float> ConfDaughTPCnclsMin{"ConfDaughTPCnclsMin", 50.f, "V0 Daugh sel: Min. nCls TPC"};
  Configurable<double> ConfDaughDCAMin{"ConfDaughDCAMin", 0.08f, "V0 Daugh sel:  Max. DCA Daugh to PV (cm)"};
  Configurable<float> ConfDaughPIDCuts{"ConfDaughPIDCuts", 3, "PID selections for Lambda daughters"};

  Configurable<int> CentNbins{"CentNbins", 16, "Number of bins in cent histograms"};
  Configurable<float> lbinCent{"lbinCent", 0.0, "lower bin value in cent histograms"};
  Configurable<float> hbinCent{"hbinCent", 80.0, "higher bin value in cent histograms"};
  Configurable<int> PolNbins{"PolNbins", 20, "Number of bins in polarisation"};
  Configurable<float> lbinPol{"lbinPol", -1.0, "lower bin value in #phi-#psi histograms"};
  Configurable<float> hbinPol{"hbinPol", 1.0, "higher bin value in #phi-#psi histograms"};
  Configurable<int> IMNbins{"IMNbins", 100, "Number of bins in invariant mass"};
  Configurable<float> lbinIM{"lbinIM", 1.0, "lower bin value in IM histograms"};
  Configurable<float> hbinIM{"hbinIM", 1.2, "higher bin value in IM histograms"};
  Configurable<int> ptNbins{"ptNbins", 50, "Number of bins in pt"};
  Configurable<float> lbinpt{"lbinpt", 0.0, "lower bin value in pt histograms"};
  Configurable<float> hbinpt{"hbinpt", 10.0, "higher bin value in pt histograms"};
  Configurable<int> resNbins{"resNbins", 50, "Number of bins in reso"};
  Configurable<float> lbinres{"lbinres", 0.0, "lower bin value in reso histograms"};
  Configurable<float> hbinres{"hbinres", 10.0, "higher bin value in reso histograms"};
  Configurable<int> etaNbins{"etaNbins", 20, "Number of bins in eta"};
  Configurable<float> lbineta{"lbineta", -1.0, "lower bin value in eta histograms"};
  Configurable<float> hbineta{"hbineta", 1.0, "higher bin value in eta histograms"};
  Configurable<int> spNbins{"spNbins", 2000, "Number of bins in sp"};
  Configurable<float> lbinsp{"lbinsp", -1.0, "lower bin value in sp histograms"};
  Configurable<float> hbinsp{"hbinsp", 1.0, "higher bin value in sp histograms"};

  ConfigurableAxis configcentAxis{"configcentAxis", {VARIABLE_WIDTH, 0.0, 10.0, 40.0, 80.0}, "Cent V0M"};
  ConfigurableAxis configthnAxispT{"configthnAxisPt", {VARIABLE_WIDTH, 0.2, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.5, 8.0, 10.0, 100.0}, "#it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis configthnAxisPol{"configthnAxisPol", {VARIABLE_WIDTH, -1.0, -0.6, -0.2, 0, 0.2, 0.4, 0.8}, "Pol"};

  SliceCache cache;
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext&)
  {
    AxisSpec thnAxisres{resNbins, lbinres, hbinres, "Reso"};
    AxisSpec thnAxisInvMass{IMNbins, lbinIM, hbinIM, "#it{M} (GeV/#it{c}^{2})"};

    histos.add("hCentrality", "Centrality distribution", kTH1F, {{configcentAxis}});

    histos.add("hSparseLambdaLambda", "hSparseLambdaLambda", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisInvMass, configthnAxisPol, configcentAxis, configthnAxispT}, true);
    histos.add("hSparseLambdaAntiLambda", "hSparseLambdaAntiLambda", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisInvMass, configthnAxisPol, configcentAxis, configthnAxispT}, true);
    histos.add("hSparseAntiLambdaAntiLambda", "hSparseAntiLambdaAntiLambda", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisInvMass, configthnAxisPol, configcentAxis, configthnAxispT}, true);
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
  template <typename V0, typename T>
  bool isSelectedV0Daughter(V0 const& candidate, T const& track, int pid)
  {
    const auto tpcNClsF = track.tpcNClsFound();
    if (track.tpcNClsCrossedRows() < 70) {
      return false;
    }
    if (tpcNClsF < ConfDaughTPCnclsMin) {
      return false;
    }
    if (track.tpcCrossedRowsOverFindableCls() < 0.8) {
      return false;
    }

    if (pid == 0 && TMath::Abs(track.tpcNSigmaPr()) > ConfDaughPIDCuts) {
      return false;
    }
    if (pid == 1 && TMath::Abs(track.tpcNSigmaPi()) > ConfDaughPIDCuts) {
      return false;
    }
    if (pid == 0 && (candidate.positivept() < cfgDaughPrPt || candidate.negativept() < cfgDaughPiPt)) {
      return false;
    }
    if (pid == 1 && (candidate.positivept() < cfgDaughPiPt || candidate.negativept() < cfgDaughPrPt)) {
      return false;
    }
    if (std::abs(candidate.positiveeta()) > ConfDaughEta || std::abs(candidate.negativeeta()) > ConfDaughEta) {
      return false;
    }

    if (pid == 0 && (TMath::Abs(candidate.dcapostopv()) < cMinV0DCAPr || TMath::Abs(candidate.dcanegtopv()) < cMinV0DCAPi)) {
      return false;
    }
    if (pid == 1 && (TMath::Abs(candidate.dcapostopv()) < cMinV0DCAPi || TMath::Abs(candidate.dcanegtopv()) < cMinV0DCAPr)) {
      return false;
    }

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

  void fillHistograms(bool tag1, bool tag2, bool tag3, bool tag4, const ROOT::Math::PxPyPzMVector& particlepair,
                      const ROOT::Math::PxPyPzMVector& particle1, const ROOT::Math::PxPyPzMVector& particle2,
                      const ROOT::Math::PxPyPzMVector& daughpart1, const ROOT::Math::PxPyPzMVector& daughpart2,
                      double centrality, double candmass, double candpt)
  {

    ROOT::Math::Boost boostPairToCM{particlepair.BoostToCM()}; // boosting vector for pair CM
    // Boosting both Lambdas to Lambda-Lambda pair rest frame
    auto lambda1_CM = boostPairToCM(particle1);
    auto lambda2_CM = boostPairToCM(particle2);

    // Step 2: Boost Each Lambda to its Own Rest Frame
    ROOT::Math::Boost boostLambda1ToCM{lambda1_CM.BoostToCM()};
    ROOT::Math::Boost boostLambda2ToCM{lambda2_CM.BoostToCM()};

    // Also boost the daughter protons to the same frame
    auto proton1_pairCM = boostPairToCM(daughpart1); // proton1 to pair CM
    auto proton2_pairCM = boostPairToCM(daughpart2); // proton2 to pair CM

    // Boost protons into their respective Lambda rest frames
    auto proton1_LambdaRF = boostLambda1ToCM(proton1_pairCM);
    auto proton2_LambdaRF = boostLambda2ToCM(proton2_pairCM);

    // Method2
    ROOT::Math::XYZVector quantizationAxis = lambda1_CM.Vect().Unit(); // Unit vector along Lambda1's direction in pair rest frame
    double cosTheta1 = proton1_LambdaRF.Vect().Unit().Dot(quantizationAxis);
    double cosTheta2 = proton2_LambdaRF.Vect().Unit().Dot(-quantizationAxis); // Opposite for Lambda2

    double theta1 = acos(cosTheta1); // angle in radians
    double theta2 = acos(cosTheta2); // angle in radians
    // Step 2: Compute sin(theta1) and sin(theta2)
    double sinTheta1 = std::sqrt(1 - cosTheta1 * cosTheta1);
    double sinTheta2 = std::sqrt(1 - cosTheta2 * cosTheta2);

    // Step 3: Calculate cos(theta1 - theta2) using the trigonometric identity
    // double cosThetaDiff = cosTheta1 * cosTheta2 + sinTheta1 * sinTheta2;
    double cosThetaDiff = TMath::Cos(theta1 - theta2);

    if (tag1 && tag3)
      histos.fill(HIST("hSparseLambdaLambda"), particle1.M(), particle2.M(), cosThetaDiff, centrality, candpt);
    if (tag1 && tag4)
      histos.fill(HIST("hSparseLambdaAntiLambda"), particle1.M(), particle2.M(), cosThetaDiff, centrality, candpt);
    if (tag2 && tag4)
      histos.fill(HIST("hSparseAntiLambdaAntiLambda"), particle1.M(), particle2.M(), cosThetaDiff, centrality, candpt);
  }

  std::tuple<int, int, bool> getLambdaTags(const auto& v0, const auto& collision)
  {
    auto postrack = v0.template posTrack_as<AllTrackCandidates>();
    auto negtrack = v0.template negTrack_as<AllTrackCandidates>();

    int LambdaTag = 0;
    int aLambdaTag = 0;

    const auto signpos = postrack.sign();
    const auto signneg = negtrack.sign();

    if (signpos < 0 || signneg > 0) {
      return {0, 0, false}; // Invalid candidate
    }

    if (isSelectedV0Daughter(v0, postrack, 0) && isSelectedV0Daughter(v0, negtrack, 1)) {
      LambdaTag = 1;
    }
    if (isSelectedV0Daughter(v0, negtrack, 0) && isSelectedV0Daughter(v0, postrack, 1)) {
      aLambdaTag = 1;
    }

    if (!LambdaTag && !aLambdaTag) {
      return {0, 0, false}; // No valid tags
    }

    if (!SelectionV0(collision, v0)) {
      return {0, 0, false}; // Fails selection
    }

    if (TMath::Abs(v0.eta()) > 0.8) {
      return {0, 0, false}; // Fails selection
    }

    return {LambdaTag, aLambdaTag, true}; // Valid candidate
  }

  ROOT::Math::PxPyPzMVector Lambda, AntiLambda, Lambdadummy, AntiLambdadummy, Proton, Pion, AntiProton, AntiPion, fourVecDauCM;
  ROOT::Math::PxPyPzMVector Lambda2, AntiLambda2, Lambdadummy2, AntiLambdadummy2, Proton2, Pion2, AntiProton2, AntiPion2;
  ROOT::Math::PxPyPzMVector LambdaLambdapair, LambdaAntiLambdapair, AntiLambdaAntiLambdapair;
  ROOT::Math::XYZVector threeVecDauCM, threeVecDauCMXY;
  double massLambda = o2::constants::physics::MassLambda;
  double massPr = o2::constants::physics::MassProton;
  double massPi = o2::constants::physics::MassPionCharged;

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter centralityFilter = (nabs(aod::cent::centFT0C) < cfgCutCentralityMax && nabs(aod::cent::centFT0C) > cfgCutCentralityMin);

  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::FV0Mults, aod::TPCMults, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0As, aod::SPCalibrationTables, aod::Mults>>;
  using AllTrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPi, aod::pidTPCFullPr>;
  using ResoV0s = aod::V0Datas;

  void processData(EventCandidates::iterator const& collision, AllTrackCandidates const& tracks, ResoV0s const& V0s, aod::BCs const&)
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

    histos.fill(HIST("hCentrality"), centrality);

    for (auto v0 : V0s) {

      auto [LambdaTag, aLambdaTag, isValid] = getLambdaTags(v0, collision);
      if (!isValid)
        continue;

      if (LambdaTag) {
        Proton = ROOT::Math::PxPyPzMVector(v0.pxpos(), v0.pypos(), v0.pzpos(), massPr);
        AntiPion = ROOT::Math::PxPyPzMVector(v0.pxneg(), v0.pyneg(), v0.pzneg(), massPi);
        Lambdadummy = Proton + AntiPion;
      }
      if (aLambdaTag) {
        AntiProton = ROOT::Math::PxPyPzMVector(v0.pxneg(), v0.pyneg(), v0.pzneg(), massPr);
        Pion = ROOT::Math::PxPyPzMVector(v0.pxpos(), v0.pypos(), v0.pzpos(), massPi);
        AntiLambdadummy = AntiProton + Pion;
      }

      if (shouldReject(LambdaTag, aLambdaTag, Lambdadummy, AntiLambdadummy)) {
        continue;
      }

      int taga = LambdaTag;
      int tagb = aLambdaTag;

      // 2nd loop for combination of lambda lambda
      for (auto v02 : V0s) {

        if (v0.v0Id() >= v02.v0Id())
          continue;

        auto [LambdaTag2, aLambdaTag2, isValid2] = getLambdaTags(v02, collision);
        if (!isValid2)
          continue;

        if (LambdaTag2) {
          Proton2 = ROOT::Math::PxPyPzMVector(v02.pxpos(), v02.pypos(), v02.pzpos(), massPr);
          AntiPion2 = ROOT::Math::PxPyPzMVector(v02.pxneg(), v02.pyneg(), v02.pzneg(), massPi);
          Lambdadummy2 = Proton2 + AntiPion2;
        }
        if (aLambdaTag2) {
          AntiProton2 = ROOT::Math::PxPyPzMVector(v02.pxneg(), v02.pyneg(), v02.pzneg(), massPr);
          Pion2 = ROOT::Math::PxPyPzMVector(v02.pxpos(), v02.pypos(), v02.pzpos(), massPi);
          AntiLambdadummy2 = AntiProton2 + Pion2;
        }

        if (shouldReject(LambdaTag2, aLambdaTag2, Lambdadummy2, AntiLambdadummy2)) {
          continue;
        }

        int taga2 = LambdaTag2;
        int tagb2 = aLambdaTag2;

        if (LambdaTag && LambdaTag2) {
          LambdaLambdapair = Lambdadummy + Lambdadummy2;
          tagb = 0;
          tagb2 = 0;
          fillHistograms(taga, tagb, taga2, tagb2, LambdaLambdapair, Lambdadummy, Lambdadummy2, Proton, Proton2, centrality, LambdaLambdapair.M(), LambdaLambdapair.Pt());
        }

        tagb2 = aLambdaTag2;

        if (LambdaTag && aLambdaTag2) {
          LambdaAntiLambdapair = Lambdadummy + AntiLambdadummy2;
          tagb = 0;
          taga2 = 0;
          fillHistograms(taga, tagb, taga2, tagb2, LambdaAntiLambdapair, Lambdadummy, AntiLambdadummy2, Proton, AntiProton2, centrality, LambdaAntiLambdapair.M(), LambdaAntiLambdapair.Pt());
        }

        tagb = aLambdaTag;
        taga2 = LambdaTag2;

        if (aLambdaTag && aLambdaTag2) {
          AntiLambdaAntiLambdapair = AntiLambdadummy + AntiLambdadummy2;
          taga = 0;
          taga2 = 0;
          fillHistograms(taga, tagb, taga2, tagb2, AntiLambdaAntiLambdapair, AntiLambdadummy, AntiLambdadummy2, AntiProton, AntiProton2, centrality, AntiLambdaAntiLambdapair.M(), AntiLambdaAntiLambdapair.Pt());
        }
      }
    }
  }
  PROCESS_SWITCH(lambdapolspspin, processData, "Process data", true);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<lambdapolspspin>(cfgc, TaskName{"lambdapolspspin"})};
}
