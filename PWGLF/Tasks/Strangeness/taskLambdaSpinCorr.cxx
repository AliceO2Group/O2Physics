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

/// \file taskLambdaSpinCorr.cxx
/// \brief Analysis task for Lambda spin spin correlation
///
/// \author prottay.das@cern.ch

#include <tuple>
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "Math/GenVector/Boost.h"

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
#include "CCDB/BasicCCDBManager.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "Common/DataModel/FT0Corrected.h"
#include "PWGMM/Mult/DataModel/Index.h" // for Particles2Tracks table

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

struct LfTaskLambdaSpinCorr {

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  // fill output
  Configurable<bool> additionalEvSel{"additionalEvSel", false, "additionalEvSel"};
  Configurable<bool> additionalEvSel3{"additionalEvSel3", false, "additionalEvSel3"};
  Configurable<bool> fillGEN{"fillGEN", true, "filling generated histograms"};

  // events
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCutCentralityMax{"cfgCutCentralityMax", 50.0f, "Accepted maximum Centrality"};
  Configurable<float> cfgCutCentralityMin{"cfgCutCentralityMin", 30.0f, "Accepted minimum Centrality"};

  // Configs for track
  Configurable<float> cfgCutPt{"cfgCutPt", 0.2, "Pt cut on daughter track"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8, "Eta cut on daughter track"};
  // Configs for V0
  Configurable<float> confV0PtMin{"confV0PtMin", 0.f, "Minimum transverse momentum of V0"};
  Configurable<float> confV0Rap{"confV0Rap", 0.8f, "Rapidity range of V0"};
  Configurable<double> confV0DCADaughMax{"confV0DCADaughMax", 0.2f, "Maximum DCA between the V0 daughters"};
  Configurable<double> confV0CPAMin{"confV0CPAMin", 0.9998f, "Minimum CPA of V0"};
  Configurable<float> confV0TranRadV0Min{"confV0TranRadV0Min", 1.5f, "Minimum transverse radius"};
  Configurable<float> confV0TranRadV0Max{"confV0TranRadV0Max", 100.f, "Maximum transverse radius"};
  Configurable<double> cMaxV0DCA{"cMaxV0DCA", 1.2, "Maximum V0 DCA to PV"};
  Configurable<double> cMinV0DCAPr{"cMinV0DCAPr", 0.05, "Minimum V0 daughters DCA to PV for Pr"};
  Configurable<double> cMinV0DCAPi{"cMinV0DCAPi", 0.05, "Minimum V0 daughters DCA to PV for Pi"};
  Configurable<float> cMaxV0LifeTime{"cMaxV0LifeTime", 20, "Maximum V0 life time"};

  // config for V0 daughters
  Configurable<float> confDaughEta{"confDaughEta", 0.8f, "V0 Daugh sel: max eta"};
  Configurable<float> cfgDaughPrPt{"cfgDaughPrPt", 0.4, "minimum daughter proton pt"};
  Configurable<float> cfgDaughPiPt{"cfgDaughPiPt", 0.2, "minimum daughter pion pt"};
  Configurable<float> confDaughTPCnclsMin{"confDaughTPCnclsMin", 50.f, "V0 Daugh sel: Min. nCls TPC"};
  Configurable<double> confDaughDCAMin{"confDaughDCAMin", 0.08f, "V0 Daugh sel:  Max. DCA Daugh to PV (cm)"};
  Configurable<float> confDaughPIDCuts{"confDaughPIDCuts", 3, "PID selections for Lambda daughters"};

  Configurable<int> iMNbins{"iMNbins", 100, "Number of bins in invariant mass"};
  Configurable<float> lbinIM{"lbinIM", 1.0, "lower bin value in IM histograms"};
  Configurable<float> hbinIM{"hbinIM", 1.2, "higher bin value in IM histograms"};
  Configurable<int> iMNbinspair{"iMNbinspair", 100, "Number of bins in invariant mass pair"};
  Configurable<float> lbinIMpair{"lbinIMpair", 1.0, "lower bin value in IMpair histograms"};
  Configurable<float> hbinIMpair{"hbinIMpair", 1.2, "higher bin value in IMpair histograms"};

  ConfigurableAxis configcentAxis{"configcentAxis", {VARIABLE_WIDTH, 0.0, 10.0, 40.0, 80.0}, "Cent V0M"};
  ConfigurableAxis configthnAxisPt{"configthnAxisPt", {VARIABLE_WIDTH, 0.2, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.5, 8.0, 10.0, 100.0}, "#it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis configthnAxisPol{"configthnAxisPol", {VARIABLE_WIDTH, -1.0, -0.6, -0.2, 0, 0.2, 0.4, 0.8}, "Pol"};

  SliceCache cache;
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext&)
  {
    AxisSpec thnAxisInvMass{iMNbins, lbinIM, hbinIM, "#it{M} (GeV/#it{c}^{2})"};
    AxisSpec thnAxisInvMasspair{iMNbinspair, lbinIMpair, hbinIMpair, "#it{M} (GeV/#it{c}^{2})"};

    histos.add("hCentrality", "Centrality distribution", kTH1F, {{configcentAxis}});

    histos.add("hSparseLambdaLambda", "hSparseLambdaLambda", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisInvMass, configthnAxisPol, configcentAxis, thnAxisInvMasspair}, true);
    histos.add("hSparseLambdaAntiLambda", "hSparseLambdaAntiLambda", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisInvMass, configthnAxisPol, configcentAxis, thnAxisInvMasspair}, true);
    histos.add("hSparseAntiLambdaAntiLambda", "hSparseAntiLambdaAntiLambda", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisInvMass, configthnAxisPol, configcentAxis, thnAxisInvMasspair}, true);
    if (fillGEN) {
      histos.add("hSparseLambdaLambdaMC", "hSparseLambdaLambdaMC", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisInvMass, configthnAxisPol, configcentAxis, thnAxisInvMasspair}, true);
      histos.add("hSparseLambdaAntiLambdaMC", "hSparseLambdaAntiLambdaMC", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisInvMass, configthnAxisPol, configcentAxis, thnAxisInvMasspair}, true);
      histos.add("hSparseAntiLambdaAntiLambdaMC", "hSparseAntiLambdaAntiLambdaMC", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisInvMass, configthnAxisPol, configcentAxis, thnAxisInvMasspair}, true);
    }
  }

  template <typename Collision, typename V0>
  bool selectionV0(Collision const& collision, V0 const& candidate)
  {
    if (std::abs(candidate.dcav0topv()) > cMaxV0DCA) {
      return false;
    }
    const float pT = candidate.pt();
    const float tranRad = candidate.v0radius();
    const float dcaDaughv0 = std::abs(candidate.dcaV0daughters());
    const float cpav0 = candidate.v0cosPA();

    float ctauLambda = candidate.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * (o2::constants::physics::MassLambda);

    if (pT < confV0PtMin) {
      return false;
    }
    if (dcaDaughv0 > confV0DCADaughMax) {
      return false;
    }
    if (cpav0 < confV0CPAMin) {
      return false;
    }
    if (tranRad < confV0TranRadV0Min) {
      return false;
    }
    if (tranRad > confV0TranRadV0Max) {
      return false;
    }
    if (std::abs(ctauLambda) > cMaxV0LifeTime) {
      return false;
    }
    if (std::abs(candidate.yLambda()) > confV0Rap) {
      return false;
    }
    return true;
  }

  template <typename V0, typename T>
  bool isSelectedV0Daughter(V0 const& candidate, T const& track, int pid)
  {
    const auto tpcNClsF = track.tpcNClsFound();
    const auto ncr = 70;
    const auto ncrfc = 0.8;

    if (track.tpcNClsCrossedRows() < ncr) {
      return false;
    }
    if (tpcNClsF < confDaughTPCnclsMin) {
      return false;
    }
    if (track.tpcCrossedRowsOverFindableCls() < ncrfc) {
      return false;
    }

    if (pid == 0 && std::abs(track.tpcNSigmaPr()) > confDaughPIDCuts) {
      return false;
    }
    if (pid == 1 && std::abs(track.tpcNSigmaPi()) > confDaughPIDCuts) {
      return false;
    }
    if (pid == 0 && (candidate.positivept() < cfgDaughPrPt || candidate.negativept() < cfgDaughPiPt)) {
      return false;
    }
    if (pid == 1 && (candidate.positivept() < cfgDaughPiPt || candidate.negativept() < cfgDaughPrPt)) {
      return false;
    }
    if (std::abs(candidate.positiveeta()) > confDaughEta || std::abs(candidate.negativeeta()) > confDaughEta) {
      return false;
    }

    if (pid == 0 && (std::abs(candidate.dcapostopv()) < cMinV0DCAPr || std::abs(candidate.dcanegtopv()) < cMinV0DCAPi)) {
      return false;
    }
    if (pid == 1 && (std::abs(candidate.dcapostopv()) < cMinV0DCAPi || std::abs(candidate.dcanegtopv()) < cMinV0DCAPr)) {
      return false;
    }

    return true;
  }

  bool shouldReject(bool lambdaTag, bool aLambdaTag,
                    const ROOT::Math::PxPyPzMVector& Lambdadummy,
                    const ROOT::Math::PxPyPzMVector& AntiLambdadummy)
  {
    const double minMass = 1.105;
    const double maxMass = 1.125;
    return (lambdaTag && aLambdaTag &&
            (Lambdadummy.M() > minMass && Lambdadummy.M() < maxMass) &&
            (AntiLambdadummy.M() > minMass && AntiLambdadummy.M() < maxMass));
  }

  void fillHistograms(bool tag1, bool tag2, bool tag3, bool tag4, const ROOT::Math::PxPyPzMVector& particlepair,
                      const ROOT::Math::PxPyPzMVector& particle1, const ROOT::Math::PxPyPzMVector& particle2,
                      const ROOT::Math::PxPyPzMVector& daughpart1, const ROOT::Math::PxPyPzMVector& daughpart2,
                      double centrality, bool datatype)
  {

    ROOT::Math::Boost boostPairToCM{particlepair.BoostToCM()}; // boosting vector for pair CM
    // Boosting both Lambdas to Lambda-Lambda pair rest frame
    auto lambda1CM = boostPairToCM(particle1);
    auto lambda2CM = boostPairToCM(particle2);

    // Step 2: Boost Each Lambda to its Own Rest Frame
    ROOT::Math::Boost boostLambda1ToCM{lambda1CM.BoostToCM()};
    ROOT::Math::Boost boostLambda2ToCM{lambda2CM.BoostToCM()};

    // Also boost the daughter protons to the same frame
    auto proton1pairCM = boostPairToCM(daughpart1); // proton1 to pair CM
    auto proton2pairCM = boostPairToCM(daughpart2); // proton2 to pair CM

    // Boost protons into their respective Lambda rest frames
    auto proton1LambdaRF = boostLambda1ToCM(proton1pairCM);
    auto proton2LambdaRF = boostLambda2ToCM(proton2pairCM);

    // Method2
    /*
    ROOT::Math::XYZVector quantizationAxis = lambda1CM.Vect().Unit(); // Unit vector along Lambda1's direction in pair rest frame
    double cosTheta1 = proton1LambdaRF.Vect().Unit().Dot(quantizationAxis);
    double cosTheta2 = proton2LambdaRF.Vect().Unit().Dot(-quantizationAxis); // Opposite for Lambda2

    double theta1 = std::acos(cosTheta1); // angle in radians
    double theta2 = std::acos(cosTheta2); // angle in radians
    double cosThetaDiff = std::cos(theta1 - theta2);
    */
    // STAR method
    double cosThetaDiff = proton1LambdaRF.Vect().Unit().Dot(proton2LambdaRF.Vect().Unit());

    auto lowptcut = 0.5;
    auto highptcut = 10.0;

    if (datatype == 1) {
      if (tag1 && tag3)
        histos.fill(HIST("hSparseLambdaLambdaMC"), particle1.M(), particle2.M(), cosThetaDiff, centrality, particlepair.M());
      if (tag1 && tag4)
        histos.fill(HIST("hSparseLambdaAntiLambdaMC"), particle1.M(), particle2.M(), cosThetaDiff, centrality, particlepair.M());
      if (tag2 && tag4)
        histos.fill(HIST("hSparseAntiLambdaAntiLambdaMC"), particle1.M(), particle2.M(), cosThetaDiff, centrality, particlepair.M());
    } else {
      if (particle1.Pt() > lowptcut && particle1.Pt() < highptcut && particle2.Pt() > lowptcut && particle2.Pt() < highptcut) {
        if (tag1 && tag3)
          histos.fill(HIST("hSparseLambdaLambda"), particle1.M(), particle2.M(), cosThetaDiff, centrality, particlepair.M());
        if (tag1 && tag4)
          histos.fill(HIST("hSparseLambdaAntiLambda"), particle1.M(), particle2.M(), cosThetaDiff, centrality, particlepair.M());
        if (tag2 && tag4)
          histos.fill(HIST("hSparseAntiLambdaAntiLambda"), particle1.M(), particle2.M(), cosThetaDiff, centrality, particlepair.M());
      }
    }
  }

  std::tuple<int, int, bool> getLambdaTags(const auto& v0, const auto& collision)
  {
    auto postrack = v0.template posTrack_as<AllTrackCandidates>();
    auto negtrack = v0.template negTrack_as<AllTrackCandidates>();

    int lambdaTag = 0;
    int aLambdaTag = 0;

    const auto signpos = postrack.sign();
    const auto signneg = negtrack.sign();

    if (signpos < 0 || signneg > 0) {
      return {0, 0, false}; // Invalid candidate
    }

    if (isSelectedV0Daughter(v0, postrack, 0) && isSelectedV0Daughter(v0, negtrack, 1)) {
      lambdaTag = 1;
    }
    if (isSelectedV0Daughter(v0, negtrack, 0) && isSelectedV0Daughter(v0, postrack, 1)) {
      aLambdaTag = 1;
    }

    if (!lambdaTag && !aLambdaTag) {
      return {0, 0, false}; // No valid tags
    }

    if (!selectionV0(collision, v0)) {
      return {0, 0, false}; // Fails selection
    }

    const auto netav = 0.8;
    if (std::abs(v0.eta()) > netav) {
      return {0, 0, false}; // Fails selection
    }

    return {lambdaTag, aLambdaTag, true}; // Valid candidate
  }

  std::tuple<int, int, bool> getLambdaTagsMC(const auto& v0, const auto& collision)
  {
    auto postrack = v0.template posTrack_as<TrackMCRecTable>();
    auto negtrack = v0.template negTrack_as<TrackMCRecTable>();

    int lambdaTag = 0;
    int aLambdaTag = 0;

    const auto signpos = postrack.sign();
    const auto signneg = negtrack.sign();

    if (signpos < 0 || signneg > 0) {
      return {0, 0, false}; // Invalid candidate
    }

    if (isSelectedV0Daughter(v0, postrack, 0) && isSelectedV0Daughter(v0, negtrack, 1)) {
      lambdaTag = 1;
    }
    if (isSelectedV0Daughter(v0, negtrack, 0) && isSelectedV0Daughter(v0, postrack, 1)) {
      aLambdaTag = 1;
    }

    if (!lambdaTag && !aLambdaTag) {
      return {0, 0, false}; // No valid tags
    }

    if (!selectionV0(collision, v0)) {
      return {0, 0, false}; // Fails selection
    }

    const auto netav = 0.8;
    if (std::abs(v0.eta()) > netav) {
      return {0, 0, false}; // Fails selection
    }

    return {lambdaTag, aLambdaTag, true}; // Valid candidate
  }

  ROOT::Math::PxPyPzMVector lambda, antiLambda, lambdadummy, antiLambdadummy, proton, pion, antiProton, antiPion, fourVecDauCM;
  ROOT::Math::PxPyPzMVector lambda2, antiLambda2, lambdadummy2, antiLambdadummy2, proton2, pion2, antiProton2, antiPion2;
  ROOT::Math::PxPyPzMVector lambdaLambdapair, lambdaAntiLambdapair, antiLambdaAntiLambdapair;
  ROOT::Math::PxPyPzMVector lambdamc, antiLambdamc, lambdadummymc, antiLambdadummymc, protonmc, pionmc, antiProtonmc, antiPionmc;
  ROOT::Math::PxPyPzMVector lambda2mcmc, antiLambda2mc, lambdadummy2mc, antiLambdadummy2mc, proton2mc, pion2mc, antiProton2mc, antiPion2mc;
  ROOT::Math::PxPyPzMVector lambdaLambdapairmc, lambdaAntiLambdapairmc, antiLambdaAntiLambdapairmc;

  ROOT::Math::XYZVector threeVecDauCM, threeVecDauCMXY;

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter centralityFilter = (nabs(aod::cent::centFT0C) < cfgCutCentralityMax && nabs(aod::cent::centFT0C) > cfgCutCentralityMin);
  // Filter acceptanceFilter = (nabs(aod::track::eta) < cfgCutEta && nabs(aod::track::pt) > cfgCutPt);

  // using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::FV0Mults, aod::TPCMults, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0As, aod::Mults>>;
  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs>>;
  using AllTrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPi, aod::pidTPCFullPr>;
  using ResoV0s = aod::V0Datas;

  void processData(EventCandidates::iterator const& collision, AllTrackCandidates const& /*tracks*/, ResoV0s const& V0s, aod::BCs const&)
  {

    if (!collision.sel8()) {
      return;
    }
    auto centrality = collision.centFT0C();
    /*if (!collision.triggereventsp()) {
      return;
      }*/

    if (additionalEvSel && (!collision.selection_bit(aod::evsel::kNoSameBunchPileup) || !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV))) {
      return;
    }

    if (additionalEvSel3 && (!collision.selection_bit(aod::evsel::kNoTimeFrameBorder) || !collision.selection_bit(aod::evsel::kNoITSROFrameBorder))) {
      return;
    }

    histos.fill(HIST("hCentrality"), centrality);

    for (const auto& v0 : V0s) {

      auto [lambdaTag, aLambdaTag, isValid] = getLambdaTags(v0, collision);
      if (!isValid)
        continue;

      if (lambdaTag) {
        proton = ROOT::Math::PxPyPzMVector(v0.pxpos(), v0.pypos(), v0.pzpos(), o2::constants::physics::MassProton);
        antiPion = ROOT::Math::PxPyPzMVector(v0.pxneg(), v0.pyneg(), v0.pzneg(), o2::constants::physics::MassPionCharged);
        lambdadummy = proton + antiPion;
      }
      if (aLambdaTag) {
        antiProton = ROOT::Math::PxPyPzMVector(v0.pxneg(), v0.pyneg(), v0.pzneg(), o2::constants::physics::MassProton);
        pion = ROOT::Math::PxPyPzMVector(v0.pxpos(), v0.pypos(), v0.pzpos(), o2::constants::physics::MassPionCharged);
        antiLambdadummy = antiProton + pion;
      }

      if (shouldReject(lambdaTag, aLambdaTag, lambdadummy, antiLambdadummy)) {
        continue;
      }

      int taga = lambdaTag;
      int tagb = aLambdaTag;

      auto postrack1 = v0.template posTrack_as<AllTrackCandidates>();
      auto negtrack1 = v0.template negTrack_as<AllTrackCandidates>();

      // 2nd loop for combination of lambda lambda
      for (const auto& v02 : V0s) {

        if (v0.v0Id() >= v02.v0Id())
          continue;

        auto [lambdaTag2, aLambdaTag2, isValid2] = getLambdaTags(v02, collision);
        if (!isValid2)
          continue;

        if (lambdaTag2) {
          proton2 = ROOT::Math::PxPyPzMVector(v02.pxpos(), v02.pypos(), v02.pzpos(), o2::constants::physics::MassProton);
          antiPion2 = ROOT::Math::PxPyPzMVector(v02.pxneg(), v02.pyneg(), v02.pzneg(), o2::constants::physics::MassPionCharged);
          lambdadummy2 = proton2 + antiPion2;
        }
        if (aLambdaTag2) {
          antiProton2 = ROOT::Math::PxPyPzMVector(v02.pxneg(), v02.pyneg(), v02.pzneg(), o2::constants::physics::MassProton);
          pion2 = ROOT::Math::PxPyPzMVector(v02.pxpos(), v02.pypos(), v02.pzpos(), o2::constants::physics::MassPionCharged);
          antiLambdadummy2 = antiProton2 + pion2;
        }

        if (shouldReject(lambdaTag2, aLambdaTag2, lambdadummy2, antiLambdadummy2)) {
          continue;
        }

        int taga2 = lambdaTag2;
        int tagb2 = aLambdaTag2;

        auto postrack2 = v02.template posTrack_as<AllTrackCandidates>();
        auto negtrack2 = v02.template negTrack_as<AllTrackCandidates>();

        if (postrack1.globalIndex() == postrack2.globalIndex() || negtrack1.globalIndex() == negtrack2.globalIndex()) {
          continue; // no shared decay products
        }

        if (lambdaTag && lambdaTag2) {
          lambdaLambdapair = lambdadummy + lambdadummy2;
          tagb = 0;
          tagb2 = 0;
          fillHistograms(taga, tagb, taga2, tagb2, lambdaLambdapair, lambdadummy, lambdadummy2, proton, proton2, centrality, 0);
        }

        tagb2 = aLambdaTag2;

        if (lambdaTag && aLambdaTag2) {
          lambdaAntiLambdapair = lambdadummy + antiLambdadummy2;
          tagb = 0;
          taga2 = 0;
          fillHistograms(taga, tagb, taga2, tagb2, lambdaAntiLambdapair, lambdadummy, antiLambdadummy2, proton, antiProton2, centrality, 0);
        }

        tagb = aLambdaTag;
        taga2 = lambdaTag2;

        if (aLambdaTag && aLambdaTag2) {
          antiLambdaAntiLambdapair = antiLambdadummy + antiLambdadummy2;
          taga = 0;
          taga2 = 0;
          fillHistograms(taga, tagb, taga2, tagb2, antiLambdaAntiLambdapair, antiLambdadummy, antiLambdadummy2, antiProton, antiProton2, centrality, 0);
        }
      }
    }
  }
  PROCESS_SWITCH(LfTaskLambdaSpinCorr, processData, "Process data", true);

  using CollisionMCTrueTable = aod::McCollisions;
  using TrackMCTrueTable = aod::McParticles;

  using CollisionMCRecTableCentFT0C = soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions, aod::CentFT0Cs, aod::EvSels, aod::PVMults>>;
  using TrackMCRecTable = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::McTrackLabels, aod::TrackSelection, aod::pidTPCFullPi, aod::pidTPCFullPr>;
  // using FilTrackMCRecTable = soa::Filtered<TrackMCRecTable>;
  using FilTrackMCRecTable = TrackMCRecTable;
  Preslice<TrackMCRecTable> perCollision = aod::track::collisionId;
  using V0TrackCandidatesMC = soa::Join<aod::V0Datas, aod::McV0Labels>;

  void processMC(CollisionMCTrueTable::iterator const& /*TrueCollision*/, CollisionMCRecTableCentFT0C const& RecCollisions, TrackMCTrueTable const& GenParticles, FilTrackMCRecTable const& /*RecTracks*/, V0TrackCandidatesMC const& V0s)
  {

    for (const auto& RecCollision : RecCollisions) {
      if (!RecCollision.sel8()) {
        continue;
      }

      if (!RecCollision.selection_bit(aod::evsel::kNoTimeFrameBorder) || !RecCollision.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
        continue;
      }

      if (std::abs(RecCollision.posZ()) > cfgCutVertex) {
        continue;
      }

      auto centrality = RecCollision.centFT0C();
      histos.fill(HIST("hCentrality"), centrality);

      for (const auto& v0 : V0s) {

        auto [lambdaTag, aLambdaTag, isValid] = getLambdaTagsMC(v0, RecCollision);
        if (!isValid)
          continue;

        if (lambdaTag) {
          proton = ROOT::Math::PxPyPzMVector(v0.pxpos(), v0.pypos(), v0.pzpos(), o2::constants::physics::MassProton);
          antiPion = ROOT::Math::PxPyPzMVector(v0.pxneg(), v0.pyneg(), v0.pzneg(), o2::constants::physics::MassPionCharged);
          lambdadummy = proton + antiPion;
        }
        if (aLambdaTag) {
          antiProton = ROOT::Math::PxPyPzMVector(v0.pxneg(), v0.pyneg(), v0.pzneg(), o2::constants::physics::MassProton);
          pion = ROOT::Math::PxPyPzMVector(v0.pxpos(), v0.pypos(), v0.pzpos(), o2::constants::physics::MassPionCharged);
          antiLambdadummy = antiProton + pion;
        }

        if (shouldReject(lambdaTag, aLambdaTag, lambdadummy, antiLambdadummy)) {
          continue;
        }

        int taga = lambdaTag;
        int tagb = aLambdaTag;

        auto postrack1 = v0.template posTrack_as<TrackMCRecTable>();
        auto negtrack1 = v0.template negTrack_as<TrackMCRecTable>();

        // 2nd loop for combination of lambda lambda
        for (const auto& v02 : V0s) {

          if (v0.v0Id() >= v02.v0Id())
            continue;

          auto [lambdaTag2, aLambdaTag2, isValid2] = getLambdaTagsMC(v02, RecCollision);
          if (!isValid2)
            continue;

          if (lambdaTag2) {
            proton2 = ROOT::Math::PxPyPzMVector(v02.pxpos(), v02.pypos(), v02.pzpos(), o2::constants::physics::MassProton);
            antiPion2 = ROOT::Math::PxPyPzMVector(v02.pxneg(), v02.pyneg(), v02.pzneg(), o2::constants::physics::MassPionCharged);
            lambdadummy2 = proton2 + antiPion2;
          }
          if (aLambdaTag2) {
            antiProton2 = ROOT::Math::PxPyPzMVector(v02.pxneg(), v02.pyneg(), v02.pzneg(), o2::constants::physics::MassProton);
            pion2 = ROOT::Math::PxPyPzMVector(v02.pxpos(), v02.pypos(), v02.pzpos(), o2::constants::physics::MassPionCharged);
            antiLambdadummy2 = antiProton2 + pion2;
          }

          if (shouldReject(lambdaTag2, aLambdaTag2, lambdadummy2, antiLambdadummy2)) {
            continue;
          }

          int taga2 = lambdaTag2;
          int tagb2 = aLambdaTag2;

          auto postrack2 = v02.template posTrack_as<TrackMCRecTable>();
          auto negtrack2 = v02.template negTrack_as<TrackMCRecTable>();

          if (postrack1.globalIndex() == postrack2.globalIndex() || negtrack1.globalIndex() == negtrack2.globalIndex()) {
            continue; // no shared decay products
          }

          if (lambdaTag && lambdaTag2) {
            lambdaLambdapair = lambdadummy + lambdadummy2;
            tagb = 0;
            tagb2 = 0;
            fillHistograms(taga, tagb, taga2, tagb2, lambdaLambdapair, lambdadummy, lambdadummy2, proton, proton2, centrality, 0);
          }

          tagb2 = aLambdaTag2;

          if (lambdaTag && aLambdaTag2) {
            lambdaAntiLambdapair = lambdadummy + antiLambdadummy2;
            tagb = 0;
            taga2 = 0;
            fillHistograms(taga, tagb, taga2, tagb2, lambdaAntiLambdapair, lambdadummy, antiLambdadummy2, proton, antiProton2, centrality, 0);
          }

          tagb = aLambdaTag;
          taga2 = lambdaTag2;

          if (aLambdaTag && aLambdaTag2) {
            antiLambdaAntiLambdapair = antiLambdadummy + antiLambdadummy2;
            taga = 0;
            taga2 = 0;
            fillHistograms(taga, tagb, taga2, tagb2, antiLambdaAntiLambdapair, antiLambdadummy, antiLambdadummy2, antiProton, antiProton2, centrality, 0);
          }
        }
      }

      //*******generated****************
      for (const auto& mcParticle : GenParticles) {
        if (std::abs(mcParticle.y()) > confV0Rap) {
          continue;
        }
        if (std::abs(mcParticle.pdgCode()) != PDG_t::kLambda0) {
          continue;
        }

        int tagamc = 0;
        int tagbmc = 0;
        int taga2mc = 0;
        int tagb2mc = 0;

        auto pdg1 = mcParticle.pdgCode();
        auto kDaughters = mcParticle.daughters_as<aod::McParticles>();
        int daughsize = 2;
        if (kDaughters.size() != daughsize) {
          continue;
        }

        for (const auto& kCurrentDaughter : kDaughters) {

          if (std::abs(kCurrentDaughter.pdgCode()) != PDG_t::kProton && std::abs(kCurrentDaughter.pdgCode()) != PDG_t::kPiPlus) {
            continue;
          }

          if (kCurrentDaughter.pdgCode() == PDG_t::kProton) {
            protonmc = ROOT::Math::PxPyPzMVector(kCurrentDaughter.px(), kCurrentDaughter.py(), kCurrentDaughter.pz(), o2::constants::physics::MassProton);
          }
          if (kCurrentDaughter.pdgCode() == PDG_t::kPiMinus) {
            antiPionmc = ROOT::Math::PxPyPzMVector(kCurrentDaughter.px(), kCurrentDaughter.py(), kCurrentDaughter.pz(), o2::constants::physics::MassPionCharged);
          }

          if (kCurrentDaughter.pdgCode() == PDG_t::kProtonBar) {
            antiProtonmc = ROOT::Math::PxPyPzMVector(kCurrentDaughter.px(), kCurrentDaughter.py(), kCurrentDaughter.pz(), o2::constants::physics::MassProton);
          }
          if (kCurrentDaughter.pdgCode() == PDG_t::kPiPlus) {
            pionmc = ROOT::Math::PxPyPzMVector(kCurrentDaughter.px(), kCurrentDaughter.py(), kCurrentDaughter.pz(), o2::constants::physics::MassPionCharged);
          }
        }
        if (pdg1 == PDG_t::kLambda0) {
          tagamc = 1;
          lambdadummymc = protonmc + antiPionmc;
        }

        if (pdg1 == PDG_t::kLambda0Bar) {
          tagbmc = 1;
          antiLambdadummymc = antiProtonmc + pionmc;
        }

        for (const auto& mcParticle2 : GenParticles) {
          if (std::abs(mcParticle2.y()) > confV0Rap) {
            continue;
          }
          if (std::abs(mcParticle2.pdgCode()) != PDG_t::kLambda0) {
            continue;
          }
          if (mcParticle.globalIndex() >= mcParticle2.globalIndex()) {
            continue;
          }

          auto pdg2 = mcParticle2.pdgCode();
          auto kDaughters2 = mcParticle2.daughters_as<aod::McParticles>();

          if (kDaughters2.size() != daughsize) {
            continue;
          }

          for (const auto& kCurrentDaughter2 : kDaughters2) {
            if (std::abs(kCurrentDaughter2.pdgCode()) != PDG_t::kProton && std::abs(kCurrentDaughter2.pdgCode()) != PDG_t::kPiPlus) {
              continue;
            }

            if (kCurrentDaughter2.pdgCode() == PDG_t::kProton) {
              proton2mc = ROOT::Math::PxPyPzMVector(kCurrentDaughter2.px(), kCurrentDaughter2.py(), kCurrentDaughter2.pz(), o2::constants::physics::MassProton);
            }
            if (kCurrentDaughter2.pdgCode() == PDG_t::kPiMinus) {
              antiPion2mc = ROOT::Math::PxPyPzMVector(kCurrentDaughter2.px(), kCurrentDaughter2.py(), kCurrentDaughter2.pz(), o2::constants::physics::MassPionCharged);
            }

            if (kCurrentDaughter2.pdgCode() == PDG_t::kProtonBar) {
              antiProton2mc = ROOT::Math::PxPyPzMVector(kCurrentDaughter2.px(), kCurrentDaughter2.py(), kCurrentDaughter2.pz(), o2::constants::physics::MassProton);
            }
            if (kCurrentDaughter2.pdgCode() == PDG_t::kPiPlus) {
              pion2mc = ROOT::Math::PxPyPzMVector(kCurrentDaughter2.px(), kCurrentDaughter2.py(), kCurrentDaughter2.pz(), o2::constants::physics::MassPionCharged);
            }
          }

          if (pdg2 == PDG_t::kLambda0) {
            taga2mc = 1;
            lambdadummy2mc = proton2mc + antiPion2mc;
          }

          if (pdg2 == PDG_t::kLambda0Bar) {
            tagb2mc = 1;
            antiLambdadummy2mc = antiProton2mc + pion2mc;
          }

          if (tagamc && taga2mc) {
            lambdaLambdapairmc = lambdadummymc + lambdadummy2mc;
            fillHistograms(tagamc, tagbmc, taga2mc, tagb2mc, lambdaLambdapairmc, lambdadummymc, lambdadummy2mc, protonmc, proton2mc, centrality, 1);
          }

          if (tagamc && tagb2mc) {
            lambdaAntiLambdapairmc = lambdadummymc + antiLambdadummy2mc;
            fillHistograms(tagamc, tagbmc, taga2mc, tagb2mc, lambdaAntiLambdapairmc, lambdadummymc, antiLambdadummy2mc, protonmc, antiProton2mc, centrality, 1);
          }

          if (tagbmc && tagb2mc) {
            antiLambdaAntiLambdapairmc = antiLambdadummymc + antiLambdadummy2mc;
            fillHistograms(tagamc, tagbmc, taga2mc, tagb2mc, antiLambdaAntiLambdapairmc, antiLambdadummymc, antiLambdadummy2mc, antiProtonmc, antiProton2mc, centrality, 1);
          }
        }
      }
    }
  }

  PROCESS_SWITCH(LfTaskLambdaSpinCorr, processMC, "Process montecarlo", false);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<LfTaskLambdaSpinCorr>(cfgc)};
}
