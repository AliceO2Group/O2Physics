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
/// \author sourav.kundu@cern.ch

#include "PWGLF/DataModel/LFSpincorrelationTables.h"

#include "Common/Core/trackUtilities.h"

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/BinningPolicy.h"
#include "Framework/StepTHn.h"
#include "Framework/runDataProcessing.h"
#include <Framework/Configurable.h>

#include <Math/GenVector/Boost.h>
#include <Math/Vector3D.h>
#include <Math/Vector4D.h>
#include <TLorentzVector.h>
#include <TMath.h>

#include <fairlogger/Logger.h>

#include <algorithm>
#include <cmath> // for std::fabs
#include <deque>
#include <iostream>
#include <iterator>
#include <set> // <<< CHANGED: for dedup sets
#include <string>
#include <type_traits>
#include <unordered_map> // <<< CHANGED: for seenMap
#include <utility>
#include <vector>

// o2 includes.
#include "CCDB/BasicCCDBManager.h"
#include "CCDB/CcdbApi.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

struct lambdaspincorrderived {
  // BinningType colBinning;
  struct : ConfigurableGroup {
    Configurable<std::string> cfgURL{"cfgURL", "http://alice-ccdb.cern.ch", "Address of the CCDB to browse"};
    Configurable<int64_t> nolaterthan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "Latest acceptable timestamp of creation for the object"};
  } cfgCcdbParam;

  // Enable access to the CCDB for the offset and correction constants and save them in dedicated variables.
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::ccdb::CcdbApi ccdbApi;
  TH3D* hweight1;
  TH3D* hweight2;
  TH3D* hweight3;
  TH3D* hweight4;

  Configurable<std::string> ConfWeightPathLL{"ConfWeightPathLL", "Users/s/skundu/My/Object/spincorr/cent010LL", "Weight path"};
  Configurable<std::string> ConfWeightPathALAL{"ConfWeightPathALAL", "Users/s/skundu/My/Object/spincorr/cent010LL", "Weight path"};
  Configurable<std::string> ConfWeightPathLAL{"ConfWeightPathLAL", "Users/s/skundu/My/Object/spincorr/cent010LL", "Weight path"};
  Configurable<std::string> ConfWeightPathALL{"ConfWeightPathALL", "Users/s/skundu/My/Object/spincorr/cent010LL", "Weight path"};

  // event sel/////////
  Configurable<float> centMin{"centMin", 0, "Minimum Centrality"};
  Configurable<float> centMax{"centMax", 80, "Maximum Centrality"};

  // Lambda selection ////////////
  Configurable<unsigned> harmonic{"harmonic", 1, "Harmonic delta phi"};
  Configurable<bool> useweight{"useweight", 1, "Use weight"};
  Configurable<bool> usePDGM{"usePDGM", 1, "Use PDG mass"};
  Configurable<bool> checkDoubleStatus{"checkDoubleStatus", 0, "Check Double status"};
  Configurable<float> cosPA{"cosPA", 0.995, "Cosine Pointing Angle"};
  Configurable<float> radiusMin{"radiusMin", 3, "Minimum V0 radius"};
  Configurable<float> radiusMax{"radiusMax", 30, "Maximum V0 radius"};
  Configurable<float> dcaProton{"dcaProton", 0.1, "DCA Proton"};
  Configurable<float> dcaPion{"dcaPion", 0.2, "DCA Pion"};
  Configurable<float> dcaDaughters{"dcaDaughters", 1.0, "DCA between daughters"};
  Configurable<float> ptMin{"ptMin", 0.5, "V0 Pt minimum"};
  Configurable<float> ptMax{"ptMax", 3.0, "V0 Pt maximum"};
  Configurable<float> rapidity{"rapidity", 0.5, "Rapidity cut on lambda"};
  Configurable<float> v0eta{"v0eta", 0.8, "Eta cut on lambda"};

  // Event Mixing
  Configurable<int> cosDef{"cosDef", 1, "Defination of cos"};
  Configurable<int> nEvtMixing{"nEvtMixing", 10, "Number of events to mix"};
  ConfigurableAxis CfgVtxBins{"CfgVtxBins", {10, -10, 10}, "Mixing bins - z-vertex"};
  ConfigurableAxis CfgMultBins{"CfgMultBins", {8, 0.0, 80}, "Mixing bins - centrality"};
  Configurable<float> etaMix{"etaMix", 0.1, "Eta cut on event mixing"};
  Configurable<float> ptMix{"ptMix", 0.1, "Pt cut on event mixing"};
  Configurable<float> phiMix{"phiMix", 0.1, "Phi cut on event mixing"};
  Configurable<float> massMix{"massMix", 0.0028, "Masscut on event mixing"};

  // THnsparse bining
  ConfigurableAxis configThnAxisInvMass{"configThnAxisInvMass", {50, 1.09, 1.14}, "#it{M} (GeV/#it{c}^{2})"};
  ConfigurableAxis configThnAxisR{"configThnAxisR", {80, 0.0, 8.0}, "#it{R}"};
  ConfigurableAxis configThnAxisPol{"configThnAxisPol", {80, 0.0, 8.0}, "cos#it{#theta *}"};
  ConfigurableAxis configThnAxisCentrality{"configThnAxisCentrality", {8, 0.0, 80.0}, "Centrality"};
  ConfigurableAxis configThnAxisRapidity{"configThnAxisRapidity", {5, 0.0, 1.0}, "Rapidity"};
  ConfigurableAxis configThnAxisPhi{"configThnAxisPhi", {18, 0.0, 2.0 * TMath::Pi()}, "Phi"};
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext&)
  {
    histos.add("hPtYSame", "hPtYSame", kTH2F, {{100, 0.0, 10.0}, {200, -1.0, 1.0}});
    histos.add("hPtYMix", "hPtYMix", kTH2F, {{100, 0.0, 10.0}, {200, -1.0, 1.0}});
    histos.add("hCentrality", "Centrality distribution", kTH1F, {{configThnAxisCentrality}});
    histos.add("deltaPhiSame", "deltaPhiSame", HistType::kTH1D, {{72, 0.0, 2.0 * TMath::Pi()}}, true);
    histos.add("deltaPhiMix", "deltaPhiMix", HistType::kTH1D, {{72, 0.0, 2.0 * TMath::Pi()}}, true);
    histos.add("ptCent", "ptCent", HistType::kTH2D, {{100, 0.0, 10.0}, {8, 0.0, 80.0}}, true);
    histos.add("etaCent", "etaCent", HistType::kTH2D, {{32, -0.8, 0.8}, {8, 0.0, 80.0}}, true);

    histos.add("hLambdaSameForLL", "hLambdaSameForLL", HistType::kTH3D, {{50, 0.0, 5.0}, {32, -0.8, 0.8}, {72, 0.0, 2.0 * TMath::Pi()}}, true);
    histos.add("hLambdaSameForLAL", "hLambdaSameForLAL", HistType::kTH3D, {{50, 0.0, 5.0}, {32, -0.8, 0.8}, {72, 0.0, 2.0 * TMath::Pi()}}, true);
    histos.add("hLambdaSameForALL", "hLambdaSameForALL", HistType::kTH3D, {{50, 0.0, 5.0}, {32, -0.8, 0.8}, {72, 0.0, 2.0 * TMath::Pi()}}, true);
    histos.add("hAntiLambdaSameForALAL", "hAntiLambdaSameForALAL", HistType::kTH3D, {{50, 0.0, 5.0}, {32, -0.8, 0.8}, {72, 0.0, 2.0 * TMath::Pi()}}, true);

    histos.add("hLambdaMixForLL", "hLambdaMixForLL", HistType::kTH3D, {{50, 0.0, 5.0}, {32, -0.8, 0.8}, {72, 0.0, 2.0 * TMath::Pi()}}, true);
    histos.add("hLambdaMixForLAL", "hLambdaMixForLAL", HistType::kTH3D, {{50, 0.0, 5.0}, {32, -0.8, 0.8}, {72, 0.0, 2.0 * TMath::Pi()}}, true);
    histos.add("hLambdaMixForALL", "hLambdaMixForALL", HistType::kTH3D, {{50, 0.0, 5.0}, {32, -0.8, 0.8}, {72, 0.0, 2.0 * TMath::Pi()}}, true);
    histos.add("hAntiLambdaMixForALAL", "hAntiLambdaMixForALAL", HistType::kTH3D, {{50, 0.0, 5.0}, {32, -0.8, 0.8}, {72, 0.0, 2.0 * TMath::Pi()}}, true);

    histos.add("hSparseLambdaLambda", "hSparseLambdaLambda", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisInvMass, configThnAxisPol, configThnAxisR}, true);
    histos.add("hSparseLambdaAntiLambda", "hSparseLambdaAntiLambda", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisInvMass, configThnAxisPol, configThnAxisR}, true);
    histos.add("hSparseAntiLambdaLambda", "hSparseAntiLambdLambda", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisInvMass, configThnAxisPol, configThnAxisR}, true);
    histos.add("hSparseAntiLambdaAntiLambda", "hSparseAntiLambdaAntiLambda", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisInvMass, configThnAxisPol, configThnAxisR}, true);

    histos.add("hSparseLambdaLambdaMixed", "hSparseLambdaLambdaMixed", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisInvMass, configThnAxisPol, configThnAxisR}, true);
    histos.add("hSparseLambdaAntiLambdaMixed", "hSparseLambdaAntiLambdaMixed", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisInvMass, configThnAxisPol, configThnAxisR}, true);
    histos.add("hSparseAntiLambdaLambdaMixed", "hSparseAntiLambdaLambdaMixed", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisInvMass, configThnAxisPol, configThnAxisR}, true);
    histos.add("hSparseAntiLambdaAntiLambdaMixed", "hSparseAntiLambdaAntiLambdaMixed", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisInvMass, configThnAxisPol, configThnAxisR}, true);

    histos.add("hSparseRapLambdaLambda", "hSparseRapLambdaLambda", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisInvMass, configThnAxisPol, configThnAxisRapidity}, true);
    histos.add("hSparseRapLambdaAntiLambda", "hSparseRapLambdaAntiLambda", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisInvMass, configThnAxisPol, configThnAxisRapidity}, true);
    histos.add("hSparseRapAntiLambdaLambda", "hSparseRapAntiLambdLambda", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisInvMass, configThnAxisPol, configThnAxisRapidity}, true);
    histos.add("hSparseRapAntiLambdaAntiLambda", "hSparseRapAntiLambdaAntiLambda", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisInvMass, configThnAxisPol, configThnAxisRapidity}, true);

    histos.add("hSparseRapLambdaLambdaMixed", "hSparseRapLambdaLambdaMixed", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisInvMass, configThnAxisPol, configThnAxisRapidity}, true);
    histos.add("hSparseRapLambdaAntiLambdaMixed", "hSparseRapLambdaAntiLambdaMixed", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisInvMass, configThnAxisPol, configThnAxisRapidity}, true);
    histos.add("hSparseRapAntiLambdaLambdaMixed", "hSparseRapAntiLambdaLambdaMixed", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisInvMass, configThnAxisPol, configThnAxisRapidity}, true);
    histos.add("hSparseRapAntiLambdaAntiLambdaMixed", "hSparseRapAntiLambdaAntiLambdaMixed", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisInvMass, configThnAxisPol, configThnAxisRapidity}, true);

    histos.add("hSparsePhiLambdaLambda", "hSparsePhiLambdaLambda", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisInvMass, configThnAxisPol, configThnAxisPhi}, true);
    histos.add("hSparsePhiLambdaAntiLambda", "hSparsePhiLambdaAntiLambda", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisInvMass, configThnAxisPol, configThnAxisPhi}, true);
    histos.add("hSparsePhiAntiLambdaLambda", "hSparsePhiAntiLambdLambda", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisInvMass, configThnAxisPol, configThnAxisPhi}, true);
    histos.add("hSparsePhiAntiLambdaAntiLambda", "hSparsePhiAntiLambdaAntiLambda", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisInvMass, configThnAxisPol, configThnAxisPhi}, true);

    histos.add("hSparsePhiLambdaLambdaMixed", "hSparsePhiLambdaLambdaMixed", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisInvMass, configThnAxisPol, configThnAxisPhi}, true);
    histos.add("hSparsePhiLambdaAntiLambdaMixed", "hSparsePhiLambdaAntiLambdaMixed", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisInvMass, configThnAxisPol, configThnAxisPhi}, true);
    histos.add("hSparsePhiAntiLambdaLambdaMixed", "hSparsePhiAntiLambdaLambdaMixed", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisInvMass, configThnAxisPol, configThnAxisPhi}, true);
    histos.add("hSparsePhiAntiLambdaAntiLambdaMixed", "hSparsePhiAntiLambdaAntiLambdaMixed", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisInvMass, configThnAxisPol, configThnAxisPhi}, true);

    ccdb->setURL(cfgCcdbParam.cfgURL);
    ccdbApi.init("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    LOGF(info, "Getting alignment offsets from the CCDB...");
    if (useweight) {
      hweight1 = ccdb->getForTimeStamp<TH3D>(ConfWeightPathLL.value, cfgCcdbParam.nolaterthan.value);
      hweight2 = ccdb->getForTimeStamp<TH3D>(ConfWeightPathLAL.value, cfgCcdbParam.nolaterthan.value);
      hweight3 = ccdb->getForTimeStamp<TH3D>(ConfWeightPathALL.value, cfgCcdbParam.nolaterthan.value);
      hweight4 = ccdb->getForTimeStamp<TH3D>(ConfWeightPathALAL.value, cfgCcdbParam.nolaterthan.value);
    }
  }

  template <typename T>
  bool selectionV0(T const& candidate)
  {
    auto particle = ROOT::Math::PtEtaPhiMVector(candidate.lambdaPt(), candidate.lambdaEta(), candidate.lambdaPhi(), candidate.lambdaMass());
    if (std::abs(particle.Rapidity()) > rapidity || std::abs(particle.Eta()) > v0eta) {
      return false;
    }
    if (candidate.v0Cospa() < cosPA) {
      return false;
    }
    if (checkDoubleStatus && candidate.doubleStatus()) {
      return false;
    }
    if (candidate.v0Radius() > radiusMax) {
      return false;
    }
    if (candidate.v0Radius() < radiusMin) {
      return false;
    }
    if (candidate.dcaBetweenDaughter() > dcaDaughters) {
      return false;
    }
    if (candidate.v0Status() == 0 && std::abs(candidate.dcaPositive()) < dcaProton && std::abs(candidate.dcaNegative()) < dcaPion) {
      return false;
    }
    if (candidate.v0Status() == 1 && std::abs(candidate.dcaPositive()) < dcaPion && std::abs(candidate.dcaNegative()) < dcaProton) {
      return false;
    }
    if (candidate.lambdaPt() < ptMin) {
      return false;
    }
    if (candidate.lambdaPt() > ptMax) {
      return false;
    }
    return true;
  }

  template <typename T1, typename T2>
  bool checkKinematics(T1 const& candidate1, T2 const& candidate2)
  {
    if (candidate1.v0Status() != candidate2.v0Status()) {
      return false;
    }
    if (std::abs(candidate1.lambdaPt() - candidate2.lambdaPt()) > ptMix) {
      return false;
    }
    if (std::abs(candidate1.lambdaEta() - candidate2.lambdaEta()) > etaMix) {
      return false;
    }
    if (std::abs(RecoDecay::constrainAngle(candidate1.lambdaPhi(), 0.0F, harmonic) - RecoDecay::constrainAngle(candidate2.lambdaPhi(), 0.0F, harmonic)) > phiMix) {
      return false;
    }
    if (std::abs(candidate1.lambdaMass() - candidate2.lambdaMass()) > massMix) {
      return false;
    }
    return true;
  }

  void fillHistograms(int tag1, int tag2,
                      const ROOT::Math::PtEtaPhiMVector& particle1, const ROOT::Math::PtEtaPhiMVector& particle2,
                      const ROOT::Math::PtEtaPhiMVector& daughpart1, const ROOT::Math::PtEtaPhiMVector& daughpart2,
                      int datatype, float mixpairweight)
  {

    auto lambda1Mass = 0.0;
    auto lambda2Mass = 0.0;
    if (!usePDGM) {
      lambda1Mass = particle1.M();
      lambda2Mass = particle2.M();
    } else {
      lambda1Mass = o2::constants::physics::MassLambda;
      lambda2Mass = o2::constants::physics::MassLambda;
    }
    auto particle1Dummy = ROOT::Math::PtEtaPhiMVector(particle1.Pt(), particle1.Eta(), particle1.Phi(), lambda1Mass);
    auto particle2Dummy = ROOT::Math::PtEtaPhiMVector(particle2.Pt(), particle2.Eta(), particle2.Phi(), lambda2Mass);
    auto pairDummy = particle1Dummy + particle2Dummy;
    ROOT::Math::Boost boostPairToCM{pairDummy.BoostToCM()}; // boosting vector for pair CM

    // Step1: Boosting both Lambdas to Lambda-Lambda pair rest frame
    auto lambda1CM = boostPairToCM(particle1Dummy);
    auto lambda2CM = boostPairToCM(particle2Dummy);

    // Step 2: Boost Each Lambda to its Own Rest Frame
    ROOT::Math::Boost boostLambda1ToCM{lambda1CM.BoostToCM()};
    ROOT::Math::Boost boostLambda2ToCM{lambda2CM.BoostToCM()};

    // Also boost the daughter protons to the same frame
    auto proton1pairCM = boostPairToCM(daughpart1); // proton1 to pair CM
    auto proton2pairCM = boostPairToCM(daughpart2); // proton2 to pair CM

    // Boost protons into their respective Lambda rest frames
    auto proton1LambdaRF = boostLambda1ToCM(proton1pairCM);
    auto proton2LambdaRF = boostLambda2ToCM(proton2pairCM);

    // =================== Opening-angle correlator: cos(Δθ) for helicity-z and beam-z ===================

    // Proton unit directions in Λ rest frames
    TVector3 k1(proton1LambdaRF.Px(), proton1LambdaRF.Py(), proton1LambdaRF.Pz());
    k1 = k1.Unit();
    TVector3 k2(proton2LambdaRF.Px(), proton2LambdaRF.Py(), proton2LambdaRF.Pz());
    k2 = k2.Unit();

    // Helper: boost a spacelike axis (t=0) from PRF into a Λ rest frame
    auto transport = [](const TVector3& v, const ROOT::Math::Boost& B) -> TVector3 {
      ROOT::Math::PxPyPzEVector a(v.X(), v.Y(), v.Z(), 0.0);
      auto ar = B(a);
      TVector3 out(ar.Px(), ar.Py(), ar.Pz());
      return (out.Mag2() > 0) ? out.Unit() : out;
    };

    // ----------------------------- (1) Helicity-z construction -----------------------------
    // z along Λ1 in PRF
    TVector3 zPRF(lambda1CM.Px(), lambda1CM.Py(), lambda1CM.Pz());
    if (zPRF.Mag2() == 0)
      zPRF = TVector3(0, 0, 1);
    zPRF = zPRF.Unit();

    // transverse axes in PRF
    TVector3 ref(0, 0, 1);
    if (std::abs(zPRF.Dot(ref)) > 0.999)
      ref = TVector3(1, 0, 0);
    TVector3 xPRF = (ref - (ref.Dot(zPRF)) * zPRF).Unit();
    TVector3 yPRF = (zPRF.Cross(xPRF)).Unit();

    // carry PRF triad to Λ rest frames (flip triad for Λ2 to keep same PRF-handedness)
    TVector3 z1_h = transport(zPRF, boostLambda1ToCM);
    TVector3 x1_h = transport(xPRF, boostLambda1ToCM);
    TVector3 y1_h = transport(yPRF, boostLambda1ToCM);

    TVector3 z2_h = transport(-zPRF, boostLambda2ToCM);
    TVector3 x2_h = transport(-xPRF, boostLambda2ToCM);
    TVector3 y2_h = transport(-yPRF, boostLambda2ToCM);

    // angles and cosΔθ (helicity)
    double c1_h = k1.Dot(z1_h);
    double s1_h = std::sqrt(std::max(0.0, 1.0 - c1_h * c1_h));
    double phi1_h = std::atan2(k1.Dot(y1_h), k1.Dot(x1_h));

    double c2_h = k2.Dot(z2_h);
    double s2_h = std::sqrt(std::max(0.0, 1.0 - c2_h * c2_h));
    double phi2_h = std::atan2(k2.Dot(y2_h), k2.Dot(x2_h));

    double cosDeltaTheta_hel = c1_h * c2_h + s1_h * s2_h * std::cos(phi1_h - phi2_h);
    if (cosDeltaTheta_hel > 1.0)
      cosDeltaTheta_hel = 1.0;
    if (cosDeltaTheta_hel < -1.0)
      cosDeltaTheta_hel = -1.0;

    // ------------------------------- (2) Beam-z construction -------------------------------
    // z along beam in PRF; choose x by projecting Λ1 onto the ⟂ plane to fix azimuth zero
    TVector3 zB(0, 0, 1);
    TVector3 L1dir(lambda1CM.Px(), lambda1CM.Py(), lambda1CM.Pz());
    L1dir = L1dir.Unit();
    TVector3 xB = L1dir - (L1dir.Dot(zB)) * zB;
    if (xB.Mag2() < 1e-12)
      xB = TVector3(1, 0, 0);
    xB = xB.Unit();
    TVector3 yB = (zB.Cross(xB)).Unit();

    // carry beam triad to Λ rest frames (no flip for a common external axis)
    TVector3 z1_b = transport(zB, boostLambda1ToCM);
    TVector3 x1_b = transport(xB, boostLambda1ToCM);
    TVector3 y1_b = transport(yB, boostLambda1ToCM);

    TVector3 z2_b = transport(zB, boostLambda2ToCM);
    TVector3 x2_b = transport(xB, boostLambda2ToCM);
    TVector3 y2_b = transport(yB, boostLambda2ToCM);

    // angles and cosΔθ (beam)
    double c1_b = k1.Dot(z1_b);
    double s1_b = std::sqrt(std::max(0.0, 1.0 - c1_b * c1_b));
    double phi1_b = std::atan2(k1.Dot(y1_b), k1.Dot(x1_b));

    double c2_b = k2.Dot(z2_b);
    double s2_b = std::sqrt(std::max(0.0, 1.0 - c2_b * c2_b));
    double phi2_b = std::atan2(k2.Dot(y2_b), k2.Dot(x2_b));

    double cosDeltaTheta_beam = c1_b * c2_b + s1_b * s2_b * std::cos(phi1_b - phi2_b);
    if (cosDeltaTheta_beam > 1.0)
      cosDeltaTheta_beam = 1.0;
    if (cosDeltaTheta_beam < -1.0)
      cosDeltaTheta_beam = -1.0;

    // --- STAR-style Δθ (as written: dot product of proton directions in their own Λ RFs) ---

    // Boost each proton into its parent's rest frame
    ROOT::Math::Boost boostL1_LabToRF{particle1Dummy.BoostToCM()}; // Λ1 velocity in lab
    ROOT::Math::Boost boostL2_LabToRF{particle2Dummy.BoostToCM()}; // Λ2 velocity in lab

    auto p1_LRF = boostL1_LabToRF(daughpart1);
    auto p2_LRF = boostL2_LabToRF(daughpart2);

    // Unit 3-vectors (in different rest frames!)
    TVector3 u1 = TVector3(p1_LRF.Px(), p1_LRF.Py(), p1_LRF.Pz()).Unit();
    TVector3 u2 = TVector3(p2_LRF.Px(), p2_LRF.Py(), p2_LRF.Pz()).Unit();

    // STAR-style cosΔθ definition
    double cosDeltaTheta_STAR_naive = u1.Dot(u2);
    if (cosDeltaTheta_STAR_naive > 1.0)
      cosDeltaTheta_STAR_naive = 1.0;
    if (cosDeltaTheta_STAR_naive < -1.0)
      cosDeltaTheta_STAR_naive = -1.0;

    auto cosThetaDiff = -999.0;
    auto costhetaz1costhetaz2 = -999.0;
    if (cosDef == 0) {
      cosThetaDiff = cosDeltaTheta_STAR_naive;
      costhetaz1costhetaz2 = (proton1LambdaRF.Pz() * proton2LambdaRF.Pz()) / (proton1LambdaRF.P() * proton2LambdaRF.P());
    } else {
      cosThetaDiff = cosDeltaTheta_hel;
      costhetaz1costhetaz2 = cosDeltaTheta_beam;
    }

    double deltaPhi = std::abs(RecoDecay::constrainAngle(particle1Dummy.Phi(), 0.0F, harmonic) - RecoDecay::constrainAngle(particle2Dummy.Phi(), 0.0F, harmonic));
    double deltaEta = particle1Dummy.Eta() - particle2Dummy.Eta();
    double deltaRap = std::abs(particle1Dummy.Rapidity() - particle2Dummy.Rapidity());
    double deltaR = TMath::Sqrt(deltaEta * deltaEta + deltaPhi * deltaPhi);

    if (datatype == 0) {
      mixpairweight = 1.0;
      histos.fill(HIST("hPtYSame"), particle1.Pt(), particle1.Rapidity(), mixpairweight);
      if (tag1 == 0 && tag2 == 0) {
        histos.fill(HIST("hSparseLambdaLambda"), particle1.M(), particle2.M(), cosThetaDiff, deltaR, mixpairweight);
        histos.fill(HIST("hSparseRapLambdaLambda"), particle1.M(), particle2.M(), cosThetaDiff, deltaRap, mixpairweight);
        histos.fill(HIST("hSparsePhiLambdaLambda"), particle1.M(), particle2.M(), costhetaz1costhetaz2, deltaPhi, mixpairweight);
        histos.fill(HIST("hLambdaSameForLL"), particle1.Pt(), particle1.Eta(), RecoDecay::constrainAngle(particle1.Phi(), 0.0F, harmonic), mixpairweight);
      } else if (tag1 == 0 && tag2 == 1) {
        histos.fill(HIST("hSparseLambdaAntiLambda"), particle1.M(), particle2.M(), cosThetaDiff, deltaR, mixpairweight);
        histos.fill(HIST("hSparseRapLambdaAntiLambda"), particle1.M(), particle2.M(), cosThetaDiff, deltaRap, mixpairweight);
        histos.fill(HIST("hSparsePhiLambdaAntiLambda"), particle1.M(), particle2.M(), costhetaz1costhetaz2, deltaPhi, mixpairweight);
        histos.fill(HIST("hLambdaSameForLAL"), particle1.Pt(), particle1.Eta(), RecoDecay::constrainAngle(particle1.Phi(), 0.0F, harmonic), mixpairweight);
      } else if (tag1 == 1 && tag2 == 0) {
        histos.fill(HIST("hSparseAntiLambdaLambda"), particle1.M(), particle2.M(), cosThetaDiff, deltaR, mixpairweight);
        histos.fill(HIST("hSparseRapAntiLambdaLambda"), particle1.M(), particle2.M(), cosThetaDiff, deltaRap, mixpairweight);
        histos.fill(HIST("hSparsePhiAntiLambdaLambda"), particle1.M(), particle2.M(), costhetaz1costhetaz2, deltaPhi, mixpairweight);
        histos.fill(HIST("hLambdaSameForALL"), particle1.Pt(), particle1.Eta(), RecoDecay::constrainAngle(particle1.Phi(), 0.0F, harmonic), mixpairweight);
      } else if (tag1 == 1 && tag2 == 1) {
        histos.fill(HIST("hSparseAntiLambdaAntiLambda"), particle1.M(), particle2.M(), cosThetaDiff, deltaR, mixpairweight);
        histos.fill(HIST("hSparseRapAntiLambdaAntiLambda"), particle1.M(), particle2.M(), cosThetaDiff, deltaRap, mixpairweight);
        histos.fill(HIST("hSparsePhiAntiLambdaAntiLambda"), particle1.M(), particle2.M(), costhetaz1costhetaz2, deltaPhi, mixpairweight);
        histos.fill(HIST("hAntiLambdaSameForALAL"), particle1.Pt(), particle1.Eta(), RecoDecay::constrainAngle(particle1.Phi(), 0.0F, harmonic), mixpairweight);
      }
    } else if (datatype == 1) {
      double weight1 = mixpairweight;
      double weight2 = mixpairweight;
      double weight3 = mixpairweight;
      double weight4 = mixpairweight;
      if (useweight) {
        weight1 = mixpairweight * hweight1->GetBinContent(hweight1->FindBin(particle1.Pt(), particle1.Eta(), RecoDecay::constrainAngle(particle1.Phi(), 0.0F, harmonic)));
        weight2 = mixpairweight * hweight2->GetBinContent(hweight2->FindBin(particle1.Pt(), particle1.Eta(), RecoDecay::constrainAngle(particle1.Phi(), 0.0F, harmonic)));
        weight3 = mixpairweight * hweight3->GetBinContent(hweight3->FindBin(particle1.Pt(), particle1.Eta(), RecoDecay::constrainAngle(particle1.Phi(), 0.0F, harmonic)));
        weight4 = mixpairweight * hweight4->GetBinContent(hweight4->FindBin(particle1.Pt(), particle1.Eta(), RecoDecay::constrainAngle(particle1.Phi(), 0.0F, harmonic)));
      }

      if (tag1 == 0 && tag2 == 0) {
        histos.fill(HIST("hPtYMix"), particle1.Pt(), particle1.Rapidity(), weight1);
        histos.fill(HIST("hSparseLambdaLambdaMixed"), particle1.M(), particle2.M(), cosThetaDiff, deltaR, weight1);
        histos.fill(HIST("hSparseRapLambdaLambdaMixed"), particle1.M(), particle2.M(), cosThetaDiff, deltaRap, weight1);
        histos.fill(HIST("hSparsePhiLambdaLambdaMixed"), particle1.M(), particle2.M(), costhetaz1costhetaz2, deltaPhi, weight1);
        histos.fill(HIST("hLambdaMixForLL"), particle1.Pt(), particle1.Eta(), RecoDecay::constrainAngle(particle1.Phi(), 0.0F, harmonic), weight1);
      } else if (tag1 == 0 && tag2 == 1) {
        histos.fill(HIST("hPtYMix"), particle1.Pt(), particle1.Rapidity(), weight2);
        histos.fill(HIST("hSparseLambdaAntiLambdaMixed"), particle1.M(), particle2.M(), cosThetaDiff, deltaR, weight2);
        histos.fill(HIST("hSparseRapLambdaAntiLambdaMixed"), particle1.M(), particle2.M(), cosThetaDiff, deltaRap, weight2);
        histos.fill(HIST("hSparsePhiLambdaAntiLambdaMixed"), particle1.M(), particle2.M(), costhetaz1costhetaz2, deltaPhi, weight2);
        histos.fill(HIST("hLambdaMixForLAL"), particle1.Pt(), particle1.Eta(), RecoDecay::constrainAngle(particle1.Phi(), 0.0F, harmonic), weight2);
      } else if (tag1 == 1 && tag2 == 0) {
        histos.fill(HIST("hPtYMix"), particle1.Pt(), particle1.Rapidity(), weight3);
        histos.fill(HIST("hSparseAntiLambdaLambdaMixed"), particle1.M(), particle2.M(), cosThetaDiff, deltaR, weight3);
        histos.fill(HIST("hSparseRapAntiLambdaLambdaMixed"), particle1.M(), particle2.M(), cosThetaDiff, deltaRap, weight3);
        histos.fill(HIST("hSparsePhiAntiLambdaLambdaMixed"), particle1.M(), particle2.M(), costhetaz1costhetaz2, deltaPhi, weight3);
        histos.fill(HIST("hLambdaMixForALL"), particle1.Pt(), particle1.Eta(), RecoDecay::constrainAngle(particle1.Phi(), 0.0F, harmonic), weight3);
      } else if (tag1 == 1 && tag2 == 1) {
        histos.fill(HIST("hPtYMix"), particle1.Pt(), particle1.Rapidity(), weight4);
        histos.fill(HIST("hSparseAntiLambdaAntiLambdaMixed"), particle1.M(), particle2.M(), cosThetaDiff, deltaR, weight4);
        histos.fill(HIST("hSparseRapAntiLambdaAntiLambdaMixed"), particle1.M(), particle2.M(), cosThetaDiff, deltaRap, weight4);
        histos.fill(HIST("hSparsePhiAntiLambdaAntiLambdaMixed"), particle1.M(), particle2.M(), costhetaz1costhetaz2, deltaPhi, weight4);
        histos.fill(HIST("hAntiLambdaMixForALAL"), particle1.Pt(), particle1.Eta(), RecoDecay::constrainAngle(particle1.Phi(), 0.0F, harmonic), weight4);
      }
    }
  }

  ROOT::Math::PtEtaPhiMVector lambda0, proton0;
  ROOT::Math::PtEtaPhiMVector lambda, proton;
  ROOT::Math::PtEtaPhiMVector lambda2, proton2;

  Filter centralityFilter = (nabs(aod::lambdaevent::cent) < centMax && nabs(aod::lambdaevent::cent) > centMin);

  using EventCandidates = soa::Filtered<aod::LambdaEvents>;
  using AllTrackCandidates = aod::LambdaPairs;

  void processData(EventCandidates::iterator const& collision, AllTrackCandidates const& V0s)
  {
    auto centrality = collision.cent();
    for (const auto& v0 : V0s) {
      if (!selectionV0(v0)) {
        continue;
      }
      histos.fill(HIST("ptCent"), v0.lambdaPt(), centrality);
      histos.fill(HIST("etaCent"), v0.lambdaEta(), centrality);
      proton = ROOT::Math::PtEtaPhiMVector(v0.protonPt(), v0.protonEta(), v0.protonPhi(), o2::constants::physics::MassProton);
      lambda = ROOT::Math::PtEtaPhiMVector(v0.lambdaPt(), v0.lambdaEta(), v0.lambdaPhi(), v0.lambdaMass());
      for (const auto& v02 : V0s) {
        if (v02.index() <= v0.index()) {
          continue;
        }
        if (!selectionV0(v02)) {
          continue;
        }
        if (v0.protonIndex() == v02.protonIndex()) {
          continue;
        }
        if (v0.pionIndex() == v02.pionIndex()) {
          continue;
        }
        proton2 = ROOT::Math::PtEtaPhiMVector(v02.protonPt(), v02.protonEta(), v02.protonPhi(), o2::constants::physics::MassProton);
        lambda2 = ROOT::Math::PtEtaPhiMVector(v02.lambdaPt(), v02.lambdaEta(), v02.lambdaPhi(), v02.lambdaMass());
        histos.fill(HIST("deltaPhiSame"), std::abs(RecoDecay::constrainAngle(v0.lambdaPhi(), 0.0F, harmonic) - RecoDecay::constrainAngle(v02.lambdaPhi(), 0.0F, harmonic)));
        if (v0.v0Status() == 0 && v02.v0Status() == 0) {
          fillHistograms(0, 0, lambda, lambda2, proton, proton2, 0, 1.0);
        }
        if (v0.v0Status() == 0 && v02.v0Status() == 1) {
          fillHistograms(0, 1, lambda, lambda2, proton, proton2, 0, 1.0);
        }
        if (v0.v0Status() == 1 && v02.v0Status() == 0) {
          fillHistograms(1, 0, lambda, lambda2, proton, proton2, 0, 1.0);
        }
        if (v0.v0Status() == 1 && v02.v0Status() == 1) {
          fillHistograms(1, 1, lambda, lambda2, proton, proton2, 0, 1.0);
        }
      }
    }
  }
  PROCESS_SWITCH(lambdaspincorrderived, processData, "Process data", true);

  // Processing Event Mixing
  SliceCache cache;
  using BinningType = ColumnBinningPolicy<aod::lambdaevent::Posz, aod::lambdaevent::Cent>;
  BinningType colBinning{{CfgVtxBins, CfgMultBins}, true};
  Preslice<aod::LambdaPairs> tracksPerCollisionV0 = aod::lambdapair::lambdaeventId;
  void processME(EventCandidates const& collisions, AllTrackCandidates const& V0s)
  {
    auto collOldIndex = -999;
    std::vector<bool> t1Used;
    for (auto& [collision1, collision2] : selfCombinations(colBinning, nEvtMixing, -1, collisions, collisions)) {
      // LOGF(info, "Mixed event collisions: (%d, %d)", collision1.index(), collision2.index());
      // auto centrality = collision1.cent();
      auto groupV01 = V0s.sliceBy(tracksPerCollisionV0, collision1.index());
      auto groupV02 = V0s.sliceBy(tracksPerCollisionV0, collision1.index());
      auto groupV03 = V0s.sliceBy(tracksPerCollisionV0, collision2.index());
      auto collNewIndex = collision1.index();
      // LOGF(info, "Mixed event collisions: (%d, %d)", collNewIndex, collOldIndex);
      if (collOldIndex != collNewIndex) {
        t1Used.resize(groupV01.size(), false);
        // std::fill(t1Used.begin(), t1Used.end(), false);
        // std::vector<bool> t1Used(groupV01.size(), false); // <-- reset here
        collOldIndex = collNewIndex;
      }
      for (auto& [t1, t3] : soa::combinations(o2::soa::CombinationsFullIndexPolicy(groupV01, groupV03))) {
        if (t1Used[t1.index()]) {
          continue;
        }
        if (!checkKinematics(t1, t3)) {
          continue;
        }
        if (!selectionV0(t1)) {
          continue;
        }
        if (!selectionV0(t3)) {
          continue;
        }
        t1Used[t1.index()] = true;
        for (const auto& t2 : groupV02) {
          if (t2.index() <= t1.index()) {
            continue;
          }
          if (!selectionV0(t2)) {
            continue;
          }
          if (t1.protonIndex() == t2.protonIndex()) {
            continue;
          }
          if (t1.pionIndex() == t2.pionIndex()) {
            continue;
          }
          proton = ROOT::Math::PtEtaPhiMVector(t3.protonPt(), t3.protonEta(), t3.protonPhi(), o2::constants::physics::MassProton);
          lambda = ROOT::Math::PtEtaPhiMVector(t3.lambdaPt(), t3.lambdaEta(), t3.lambdaPhi(), t3.lambdaMass());
          proton2 = ROOT::Math::PtEtaPhiMVector(t2.protonPt(), t2.protonEta(), t2.protonPhi(), o2::constants::physics::MassProton);
          lambda2 = ROOT::Math::PtEtaPhiMVector(t2.lambdaPt(), t2.lambdaEta(), t2.lambdaPhi(), t2.lambdaMass());
          histos.fill(HIST("deltaPhiMix"), std::abs(RecoDecay::constrainAngle(t3.lambdaPhi(), 0.0F, harmonic) - RecoDecay::constrainAngle(t2.lambdaPhi(), 0.0F, harmonic)));
          if (t3.v0Status() == 0 && t2.v0Status() == 0) {
            fillHistograms(0, 0, lambda, lambda2, proton, proton2, 1, 1.0);
          }
          if (t3.v0Status() == 0 && t2.v0Status() == 1) {
            fillHistograms(0, 1, lambda, lambda2, proton, proton2, 1, 1.0);
          }
          if (t3.v0Status() == 1 && t2.v0Status() == 0) {
            fillHistograms(1, 0, lambda, lambda2, proton, proton2, 1, 1.0);
          }
          if (t3.v0Status() == 1 && t2.v0Status() == 1) {
            fillHistograms(1, 1, lambda, lambda2, proton, proton2, 1, 1.0);
          }
        }
      } // replacement track pair
    } // collision pair
  }
  PROCESS_SWITCH(lambdaspincorrderived, processME, "Process data ME", false);

  void processMEV2(EventCandidates const& collisions, AllTrackCandidates const& V0s)
  {
    auto nBins = colBinning.getAllBinsCount();
    std::vector<std::deque<std::pair<int, AllTrackCandidates>>> eventPools(nBins);

    for (auto& collision1 : collisions) {
      int bin = colBinning.getBin(std::make_tuple(collision1.posz(), collision1.cent()));
      auto poolA = V0s.sliceBy(tracksPerCollisionV0, collision1.index());
      // float centrality = collision1.cent();

      // <<< CHANGED: map old collision index → set of (t2.idx, t3.idx) we've already filled
      std::unordered_map<int, std::set<std::pair<int, int>>> seenMap;

      for (auto& [t1, t2] : soa::combinations(o2::soa::CombinationsFullIndexPolicy(poolA, poolA))) {
        if (!selectionV0(t1) || !selectionV0(t2))
          continue;
        if (t2.index() <= t1.index())
          continue;
        if (t1.protonIndex() == t2.protonIndex())
          continue;
        if (t1.pionIndex() == t2.pionIndex())
          continue;

        int mixes = 0;
        for (auto it = eventPools[bin].rbegin(); it != eventPools[bin].rend() && mixes < nEvtMixing; ++it, ++mixes) {
          int collision2idx = it->first;
          AllTrackCandidates& poolB = it->second;

          int nRepl = 0;
          for (auto& t3 : poolB) {
            if (selectionV0(t3) && checkKinematics(t1, t3)) {
              ++nRepl;
            }
          }
          if (nRepl == 0)
            continue;
          float invN = 1.0f / static_cast<float>(nRepl);

          for (auto& t3 : poolB) {
            if (!(selectionV0(t3) && checkKinematics(t1, t3))) {
              continue;
            }
            if (collision1.index() == collision2idx) {
              continue;
            }

            // <<< CHANGED: dedupe (t2, t3) pairs per prior collision
            auto key = std::make_pair(t2.index(), t3.index());
            auto& seen = seenMap[collision2idx];
            if (!seen.insert(key).second) {
              continue;
            }

            // reconstruct 4-vectors
            proton = ROOT::Math::PtEtaPhiMVector(t3.protonPt(), t3.protonEta(), t3.protonPhi(), o2::constants::physics::MassProton);
            lambda = ROOT::Math::PtEtaPhiMVector(t3.lambdaPt(), t3.lambdaEta(), t3.lambdaPhi(), t3.lambdaMass());
            proton2 = ROOT::Math::PtEtaPhiMVector(t2.protonPt(), t2.protonEta(), t2.protonPhi(), o2::constants::physics::MassProton);
            lambda2 = ROOT::Math::PtEtaPhiMVector(t2.lambdaPt(), t2.lambdaEta(), t2.lambdaPhi(), t2.lambdaMass());

            float dPhi = std::fabs(RecoDecay::constrainAngle(lambda.Phi(), 0.0F, harmonic) - RecoDecay::constrainAngle(lambda2.Phi(), 0.0F, harmonic));
            histos.fill(HIST("deltaPhiMix"), dPhi, invN);

            if (t3.v0Status() == 0 && t2.v0Status() == 0) {
              fillHistograms(0, 0, lambda, lambda2, proton, proton2, 1, invN);
            }
            if (t3.v0Status() == 0 && t2.v0Status() == 1) {
              fillHistograms(0, 1, lambda, lambda2, proton, proton2, 1, invN);
            }
            if (t3.v0Status() == 1 && t2.v0Status() == 0) {
              fillHistograms(1, 0, lambda, lambda2, proton, proton2, 1, invN);
            }
            if (t3.v0Status() == 1 && t2.v0Status() == 1) {
              fillHistograms(1, 1, lambda, lambda2, proton, proton2, 1, invN);
            }
          }
        } // end mixing-event loop
      } // end same-event pair loop

      auto sliced = V0s.sliceBy(tracksPerCollisionV0, collision1.index());
      eventPools[bin].emplace_back(collision1.index(), std::move(sliced));
      if (static_cast<int>(eventPools[bin].size()) > nEvtMixing) {
        eventPools[bin].pop_front();
      }
    } // end primary-event loop
  }
  PROCESS_SWITCH(lambdaspincorrderived, processMEV2, "Process data ME", false);

  void processMEV3(EventCandidates const& collisions, AllTrackCandidates const& V0s)
  {
    // one pool (deque) per mixing bin; each entry holds (collision index, slice of its V0s)
    auto nBins = colBinning.getAllBinsCount();
    std::vector<std::deque<std::pair<int, AllTrackCandidates>>> eventPools(nBins);

    for (auto& collision1 : collisions) {
      // select mixing bin for this event
      const int bin = colBinning.getBin(std::make_tuple(collision1.posz(), collision1.cent()));

      // all V0s from the current event
      auto poolA = V0s.sliceBy(tracksPerCollisionV0, collision1.index());

      // loop over same-event candidate pairs (t1,t2)
      for (auto& [t1, t2] : soa::combinations(o2::soa::CombinationsFullIndexPolicy(poolA, poolA))) {
        if (!selectionV0(t1) || !selectionV0(t2))
          continue;
        if (t2.index() <= t1.index())
          continue; // unique unordered pairs
        if (t1.protonIndex() == t2.protonIndex())
          continue; // no shared daughter
        if (t1.pionIndex() == t2.pionIndex())
          continue;

        // --- First pass over previous events in this bin: count replacements and build a list of usable pools
        int mixes = 0;
        struct PoolView {
          AllTrackCandidates* pool;
          int nRepl;
          int collIdx;
        };
        std::vector<PoolView> usable;
        int totalRepl = 0;

        for (auto it = eventPools[bin].rbegin();
             it != eventPools[bin].rend() && mixes < nEvtMixing; ++it, ++mixes) {
          const int collision2idx = it->first;
          auto& poolB = it->second;

          // (defensive; shouldn't happen because we push the current event after mixing)
          if (collision2idx == collision1.index())
            continue;

          int nRepl = 0;
          for (auto& t3 : poolB) {
            if (selectionV0(t3) && checkKinematics(t1, t3))
              ++nRepl;
          }
          if (nRepl > 0) {
            usable.push_back(PoolView{&poolB, nRepl, collision2idx});
            totalRepl += nRepl;
          }
        }

        if (totalRepl == 0)
          continue;
        const float w = 1.0f / static_cast<float>(totalRepl); // global normalization: sum of weights over all replacements = 1

        // --- Second pass: fill with normalized weight w
        for (auto& pv : usable) {
          auto& poolB = *pv.pool;
          for (auto& t3 : poolB) {
            if (!(selectionV0(t3) && checkKinematics(t1, t3)))
              continue;

            // build 4-vectors for the mixed pair (t3 from prior event replaces t1; t2 stays from current event)
            proton = ROOT::Math::PtEtaPhiMVector(t3.protonPt(), t3.protonEta(), t3.protonPhi(), o2::constants::physics::MassProton);
            lambda = ROOT::Math::PtEtaPhiMVector(t3.lambdaPt(), t3.lambdaEta(), t3.lambdaPhi(), t3.lambdaMass());
            proton2 = ROOT::Math::PtEtaPhiMVector(t2.protonPt(), t2.protonEta(), t2.protonPhi(), o2::constants::physics::MassProton);
            lambda2 = ROOT::Math::PtEtaPhiMVector(t2.lambdaPt(), t2.lambdaEta(), t2.lambdaPhi(), t2.lambdaMass());

            float dPhi = std::fabs(
              RecoDecay::constrainAngle(lambda.Phi(), 0.0F, harmonic) -
              RecoDecay::constrainAngle(lambda2.Phi(), 0.0F, harmonic));
            histos.fill(HIST("deltaPhiMix"), dPhi, w);

            if (t3.v0Status() == 0 && t2.v0Status() == 0) {
              fillHistograms(0, 0, lambda, lambda2, proton, proton2, 1, w);
            } else if (t3.v0Status() == 0 && t2.v0Status() == 1) {
              fillHistograms(0, 1, lambda, lambda2, proton, proton2, 1, w);
            } else if (t3.v0Status() == 1 && t2.v0Status() == 0) {
              fillHistograms(1, 0, lambda, lambda2, proton, proton2, 1, w);
            } else if (t3.v0Status() == 1 && t2.v0Status() == 1) {
              fillHistograms(1, 1, lambda, lambda2, proton, proton2, 1, w);
            }
          }
        }
      } // end same-event pair loop

      // after mixing with prior events, push current event into the pool
      auto sliced = V0s.sliceBy(tracksPerCollisionV0, collision1.index());
      eventPools[bin].emplace_back(collision1.index(), std::move(sliced));
      if (static_cast<int>(eventPools[bin].size()) > nEvtMixing) {
        eventPools[bin].pop_front();
      }
    } // end primary-event loop
  }
  PROCESS_SWITCH(lambdaspincorrderived, processMEV3, "Process data ME", false);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<lambdaspincorrderived>(cfgc)};
}
