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
#include <TRandom3.h>

#include <fairlogger/Logger.h>

#include <algorithm>
#include <cmath> // for std::fabs
#include <cstdint>
#include <deque>
#include <iostream>
#include <iterator>
#include <limits>
#include <random>
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

  TH3D* hweight12;
  TH3D* hweight22;
  TH3D* hweight32;
  TH3D* hweight42;

  Configurable<std::string> ConfWeightPathLL{"ConfWeightPathLL", "Users/s/skundu/My/Object/spincorr/cent010LL", "Weight path"};
  Configurable<std::string> ConfWeightPathALAL{"ConfWeightPathALAL", "Users/s/skundu/My/Object/spincorr/cent010LL", "Weight path"};
  Configurable<std::string> ConfWeightPathLAL{"ConfWeightPathLAL", "Users/s/skundu/My/Object/spincorr/cent010LL", "Weight path"};
  Configurable<std::string> ConfWeightPathALL{"ConfWeightPathALL", "Users/s/skundu/My/Object/spincorr/cent010LL", "Weight path"};

  Configurable<std::string> ConfWeightPathLL2{"ConfWeightPathLL2", "Users/s/skundu/My/Object/spincorr/cent010LL", "Weight path 2"};
  Configurable<std::string> ConfWeightPathALAL2{"ConfWeightPathALAL2", "Users/s/skundu/My/Object/spincorr/cent010LL", "Weight path 2"};
  Configurable<std::string> ConfWeightPathLAL2{"ConfWeightPathLAL2", "Users/s/skundu/My/Object/spincorr/cent010LL", "Weight path 2"};
  Configurable<std::string> ConfWeightPathALL2{"ConfWeightPathALL2", "Users/s/skundu/My/Object/spincorr/cent010LL", "Weight path 2"};

  // event sel/////////
  Configurable<float> centMin{"centMin", 0, "Minimum Centrality"};
  Configurable<float> centMax{"centMax", 80, "Maximum Centrality"};
  Configurable<int> rngSeed{"rngSeed", 12345, "Seed for random mixing (reproducible)"};
  std::mt19937 rng{12345};
  // Lambda selection ////////////
  Configurable<unsigned> harmonic{"harmonic", 1, "Harmonic phi"};
  Configurable<unsigned> harmonicDphi{"harmonicDphi", 2, "Harmonic delta phi"};
  Configurable<bool> useweight{"useweight", 0, "Use weight"};
  Configurable<bool> usebothweight{"usebothweight", 1, "Use both weight"};
  // Configurable<bool> useNUA{"useNUA", 0, "Use NUA weight"};
  Configurable<bool> usePDGM{"usePDGM", 1, "Use PDG mass"};
  Configurable<bool> useAdditionalHisto{"useAdditionalHisto", 0, "Use additional histogram"};
  Configurable<bool> checkDoubleStatus{"checkDoubleStatus", 0, "Check Double status"};
  Configurable<float> cosPA{"cosPA", 0.995, "Cosine Pointing Angle"};
  Configurable<float> radiusMin{"radiusMin", 3, "Minimum V0 radius"};
  Configurable<float> radiusMax{"radiusMax", 30, "Maximum V0 radius"};
  Configurable<float> dcaProton{"dcaProton", 0.1, "DCA Proton"};
  Configurable<float> dcaPion{"dcaPion", 0.2, "DCA Pion"};
  Configurable<float> dcaDaughters{"dcaDaughters", 1.0, "DCA between daughters"};
  Configurable<float> ptMin{"ptMin", 0.5, "V0 Pt minimum"};
  Configurable<float> ptMax{"ptMax", 3.0, "V0 Pt maximum"};
  Configurable<float> MassMin{"MassMin", 1.09, "V0 Mass minimum"};
  Configurable<float> MassMax{"MassMax", 1.14, "V0 Mass maximum"};
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

  ConfigurableAxis ax_dphi_h{"ax_dphi_h", {72, 0.0, 2.0 * TMath::Pi()}, "Δφ_h"};
  ConfigurableAxis ax_deta{"ax_deta", {40, -1.0, 1.0}, "Δη"};
  ConfigurableAxis ax_ptpair{"ax_ptpair", {100, 0.0, 10.0}, "p_{T,pair} (GeV/c)"};

  // THnsparse bining
  ConfigurableAxis configThnAxisInvMass{"configThnAxisInvMass", {50, 1.09, 1.14}, "#it{M} (GeV/#it{c}^{2})"};
  ConfigurableAxis configThnAxisR{"configThnAxisR", {80, 0.0, 8.0}, "#it{R}"};
  ConfigurableAxis configThnAxisPol{"configThnAxisPol", {80, 0.0, 8.0}, "cos#it{#theta *}"};
  ConfigurableAxis configThnAxisCentrality{"configThnAxisCentrality", {8, 0.0, 80.0}, "Centrality"};
  ConfigurableAxis configThnAxisRapidity{"configThnAxisRapidity", {5, 0.0, 1.0}, "Rapidity"};
  ConfigurableAxis configThnAxisPairMass{"configThnAxisPairMass", {100, 2.0, 3.0}, "PairMass"};
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

    // --- 3D SE/ME pair-space maps per category (LL, LAL, ALL, ALAL)
    histos.add("SE_LL", "SE pairs", HistType::kTH3D, {ax_dphi_h, ax_deta, ax_ptpair}, true);
    histos.add("SE_LAL", "SE pairs", HistType::kTH3D, {ax_dphi_h, ax_deta, ax_ptpair}, true);
    histos.add("SE_ALL", "SE pairs", HistType::kTH3D, {ax_dphi_h, ax_deta, ax_ptpair}, true);
    histos.add("SE_ALAL", "SE pairs", HistType::kTH3D, {ax_dphi_h, ax_deta, ax_ptpair}, true);

    histos.add("ME_LL", "ME pairs", HistType::kTH3D, {ax_dphi_h, ax_deta, ax_ptpair}, true);
    histos.add("ME_LAL", "ME pairs", HistType::kTH3D, {ax_dphi_h, ax_deta, ax_ptpair}, true);
    histos.add("ME_ALL", "ME pairs", HistType::kTH3D, {ax_dphi_h, ax_deta, ax_ptpair}, true);
    histos.add("ME_ALAL", "ME pairs", HistType::kTH3D, {ax_dphi_h, ax_deta, ax_ptpair}, true);

    histos.add("SE_LL2", "SE pairs 2", HistType::kTH3D, {ax_dphi_h, ax_deta, ax_ptpair}, true);
    histos.add("SE_LAL2", "SE pairs 2", HistType::kTH3D, {ax_dphi_h, ax_deta, ax_ptpair}, true);
    histos.add("SE_ALL2", "SE pairs 2", HistType::kTH3D, {ax_dphi_h, ax_deta, ax_ptpair}, true);
    histos.add("SE_ALAL2", "SE pairs 2", HistType::kTH3D, {ax_dphi_h, ax_deta, ax_ptpair}, true);

    histos.add("ME_LL2", "ME pairs 2", HistType::kTH3D, {ax_dphi_h, ax_deta, ax_ptpair}, true);
    histos.add("ME_LAL2", "ME pairs 2", HistType::kTH3D, {ax_dphi_h, ax_deta, ax_ptpair}, true);
    histos.add("ME_ALL2", "ME pairs 2", HistType::kTH3D, {ax_dphi_h, ax_deta, ax_ptpair}, true);
    histos.add("ME_ALAL2", "ME pairs 2", HistType::kTH3D, {ax_dphi_h, ax_deta, ax_ptpair}, true);

    histos.add("hSparseLambdaLambda", "hSparseLambdaLambda", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisInvMass, configThnAxisPol, configThnAxisR}, true);
    histos.add("hSparseLambdaAntiLambda", "hSparseLambdaAntiLambda", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisInvMass, configThnAxisPol, configThnAxisR}, true);
    histos.add("hSparseAntiLambdaLambda", "hSparseAntiLambdLambda", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisInvMass, configThnAxisPol, configThnAxisR}, true);
    histos.add("hSparseAntiLambdaAntiLambda", "hSparseAntiLambdaAntiLambda", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisInvMass, configThnAxisPol, configThnAxisR}, true);

    histos.add("hSparseLambdaLambdaMixed", "hSparseLambdaLambdaMixed", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisInvMass, configThnAxisPol, configThnAxisR}, true);
    histos.add("hSparseLambdaAntiLambdaMixed", "hSparseLambdaAntiLambdaMixed", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisInvMass, configThnAxisPol, configThnAxisR}, true);
    histos.add("hSparseAntiLambdaLambdaMixed", "hSparseAntiLambdaLambdaMixed", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisInvMass, configThnAxisPol, configThnAxisR}, true);
    histos.add("hSparseAntiLambdaAntiLambdaMixed", "hSparseAntiLambdaAntiLambdaMixed", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisInvMass, configThnAxisPol, configThnAxisR}, true);
    if (useAdditionalHisto) {
      histos.add("hSparseRapLambdaLambda", "hSparseRapLambdaLambda", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisInvMass, configThnAxisPol, configThnAxisRapidity}, true);
      histos.add("hSparseRapLambdaAntiLambda", "hSparseRapLambdaAntiLambda", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisInvMass, configThnAxisPol, configThnAxisRapidity}, true);
      histos.add("hSparseRapAntiLambdaLambda", "hSparseRapAntiLambdLambda", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisInvMass, configThnAxisPol, configThnAxisRapidity}, true);
      histos.add("hSparseRapAntiLambdaAntiLambda", "hSparseRapAntiLambdaAntiLambda", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisInvMass, configThnAxisPol, configThnAxisRapidity}, true);

      histos.add("hSparseRapLambdaLambdaMixed", "hSparseRapLambdaLambdaMixed", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisInvMass, configThnAxisPol, configThnAxisRapidity}, true);
      histos.add("hSparseRapLambdaAntiLambdaMixed", "hSparseRapLambdaAntiLambdaMixed", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisInvMass, configThnAxisPol, configThnAxisRapidity}, true);
      histos.add("hSparseRapAntiLambdaLambdaMixed", "hSparseRapAntiLambdaLambdaMixed", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisInvMass, configThnAxisPol, configThnAxisRapidity}, true);
      histos.add("hSparseRapAntiLambdaAntiLambdaMixed", "hSparseRapAntiLambdaAntiLambdaMixed", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisInvMass, configThnAxisPol, configThnAxisRapidity}, true);

      histos.add("hSparsePhiLambdaLambda", "hSparsePhiLambdaLambda", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisInvMass, configThnAxisPol, ax_dphi_h}, true);
      histos.add("hSparsePhiLambdaAntiLambda", "hSparsePhiLambdaAntiLambda", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisInvMass, configThnAxisPol, ax_dphi_h}, true);
      histos.add("hSparsePhiAntiLambdaLambda", "hSparsePhiAntiLambdLambda", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisInvMass, configThnAxisPol, ax_dphi_h}, true);
      histos.add("hSparsePhiAntiLambdaAntiLambda", "hSparsePhiAntiLambdaAntiLambda", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisInvMass, configThnAxisPol, ax_dphi_h}, true);

      histos.add("hSparsePhiLambdaLambdaMixed", "hSparsePhiLambdaLambdaMixed", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisInvMass, configThnAxisPol, ax_dphi_h}, true);
      histos.add("hSparsePhiLambdaAntiLambdaMixed", "hSparsePhiLambdaAntiLambdaMixed", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisInvMass, configThnAxisPol, ax_dphi_h}, true);
      histos.add("hSparsePhiAntiLambdaLambdaMixed", "hSparsePhiAntiLambdaLambdaMixed", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisInvMass, configThnAxisPol, ax_dphi_h}, true);
      histos.add("hSparsePhiAntiLambdaAntiLambdaMixed", "hSparsePhiAntiLambdaAntiLambdaMixed", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisInvMass, configThnAxisPol, ax_dphi_h}, true);

      histos.add("hSparsePairMassLambdaLambda", "hSparsePairMassLambdaLambda", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisInvMass, configThnAxisPol, configThnAxisPairMass}, true);
      histos.add("hSparsePairMassLambdaAntiLambda", "hSparsePairMassLambdaAntiLambda", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisInvMass, configThnAxisPol, configThnAxisPairMass}, true);
      histos.add("hSparsePairMassAntiLambdaLambda", "hSparsePairMassAntiLambdLambda", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisInvMass, configThnAxisPol, configThnAxisPairMass}, true);
      histos.add("hSparsePairMassAntiLambdaAntiLambda", "hSparsePairMassAntiLambdaAntiLambda", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisInvMass, configThnAxisPol, configThnAxisPairMass}, true);

      histos.add("hSparsePairMassLambdaLambdaMixed", "hSparsePairMassLambdaLambdaMixed", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisInvMass, configThnAxisPol, configThnAxisPairMass}, true);
      histos.add("hSparsePairMassLambdaAntiLambdaMixed", "hSparsePairMassLambdaAntiLambdaMixed", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisInvMass, configThnAxisPol, configThnAxisPairMass}, true);
      histos.add("hSparsePairMassAntiLambdaLambdaMixed", "hSparsePairMassAntiLambdaLambdaMixed", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisInvMass, configThnAxisPol, configThnAxisPairMass}, true);
      histos.add("hSparsePairMassAntiLambdaAntiLambdaMixed", "hSparsePairMassAntiLambdaAntiLambdaMixed", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisInvMass, configThnAxisPol, configThnAxisPairMass}, true);
    }
    rng.seed(static_cast<uint32_t>(rngSeed.value));
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

      hweight12 = ccdb->getForTimeStamp<TH3D>(ConfWeightPathLL2.value, cfgCcdbParam.nolaterthan.value);
      hweight22 = ccdb->getForTimeStamp<TH3D>(ConfWeightPathLAL2.value, cfgCcdbParam.nolaterthan.value);
      hweight32 = ccdb->getForTimeStamp<TH3D>(ConfWeightPathALL2.value, cfgCcdbParam.nolaterthan.value);
      hweight42 = ccdb->getForTimeStamp<TH3D>(ConfWeightPathALAL2.value, cfgCcdbParam.nolaterthan.value);
    }
  }

  template <typename T>
  bool selectionV0(T const& candidate)
  {
    auto particle = ROOT::Math::PtEtaPhiMVector(candidate.lambdaPt(), candidate.lambdaEta(), candidate.lambdaPhi(), candidate.lambdaMass());
    if (std::abs(particle.Rapidity()) > rapidity || std::abs(particle.Eta()) > v0eta) {
      return false;
    }
    if (candidate.lambdaMass() < MassMin || candidate.lambdaMass() > MassMax) {
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

    // --- STAR-style Δθ (as written: dot product of proton directions in their own Λ RFs) ---

    // Boost each proton into its parent's rest frame
    ROOT::Math::Boost boostL1_LabToRF{particle1Dummy.BoostToCM()}; // Λ1 velocity in lab
    ROOT::Math::Boost boostL2_LabToRF{particle2Dummy.BoostToCM()}; // Λ2 velocity in lab

    auto p1_LRF = boostL1_LabToRF(daughpart1);
    auto p2_LRF = boostL2_LabToRF(daughpart2);

    // Unit 3-vectors (in different rest frames!)
    TVector3 u1 = TVector3(p1_LRF.Px(), p1_LRF.Py(), p1_LRF.Pz()).Unit();
    TVector3 u2 = TVector3(p2_LRF.Px(), p2_LRF.Py(), p2_LRF.Pz()).Unit();

    // Proton unit directions in Λ rest frames
    TVector3 k1(proton1LambdaRF.Px(), proton1LambdaRF.Py(), proton1LambdaRF.Pz());
    k1 = k1.Unit();
    TVector3 k2(proton2LambdaRF.Px(), proton2LambdaRF.Py(), proton2LambdaRF.Pz());
    k2 = k2.Unit();

    // STAR-style cosΔθ definition
    double cosDeltaTheta_STAR_naive = u1.Dot(u2);
    if (cosDeltaTheta_STAR_naive > 1.0)
      cosDeltaTheta_STAR_naive = 111.0;
    if (cosDeltaTheta_STAR_naive < -1.0)
      cosDeltaTheta_STAR_naive = -111.0;

    double cosDeltaTheta_hel = k1.Dot(k2);
    if (cosDeltaTheta_hel > 1.0)
      cosDeltaTheta_hel = 111.0;
    if (cosDeltaTheta_hel < -1.0)
      cosDeltaTheta_hel = -111.0;

    auto cosThetaDiff = -999.0;
    if (cosDef == 0) {
      cosThetaDiff = cosDeltaTheta_STAR_naive;
    } else {
      cosThetaDiff = cosDeltaTheta_hel;
    }

    double pt1 = particle1.Pt();
    double dphi1 = RecoDecay::constrainAngle(particle1.Phi(), 0.0F, harmonic);
    double deta1 = particle1.Eta();

    double pt2 = particle2.Pt();
    double dphi2 = RecoDecay::constrainAngle(particle2.Phi(), 0.0F, harmonic);
    double deta2 = particle2.Eta();

    double deta_pair = std::abs(deta1 - deta2);
    double dphi_pair = RecoDecay::constrainAngle(particle1.Phi() - particle2.Phi(), 0.0F, harmonicDphi);

    double deltaR = TMath::Sqrt(deta_pair * deta_pair + dphi_pair * dphi_pair);
    double deltaRap = std::abs(particle1.Rapidity() - particle2.Rapidity());

    double epsWeight1 = 1.0;
    double epsWeight2 = 1.0;

    if (useweight && datatype == 1) {
      if (tag1 == 0 && tag2 == 0) {
        epsWeight1 = hweight1->GetBinContent(hweight1->FindBin(dphi1, deta1, pt1));
        epsWeight2 = hweight12->GetBinContent(hweight12->FindBin(dphi2, deta2, pt2));
      } else if (tag1 == 0 && tag2 == 1) {
        epsWeight1 = hweight2->GetBinContent(hweight2->FindBin(dphi1, deta1, pt1));
        epsWeight2 = hweight22->GetBinContent(hweight22->FindBin(dphi2, deta2, pt2));
      } else if (tag1 == 1 && tag2 == 0) {
        epsWeight1 = hweight3->GetBinContent(hweight3->FindBin(dphi1, deta1, pt1));
        epsWeight2 = hweight32->GetBinContent(hweight32->FindBin(dphi2, deta2, pt2));
      } else if (tag1 == 1 && tag2 == 1) {
        epsWeight1 = hweight4->GetBinContent(hweight4->FindBin(dphi1, deta1, pt1));
        epsWeight2 = hweight42->GetBinContent(hweight42->FindBin(dphi2, deta2, pt2));
      }
    }

    if (datatype == 0) {
      mixpairweight = 1.0;
      histos.fill(HIST("hPtYSame"), particle1.Pt(), particle1.Rapidity(), mixpairweight);
      if (tag1 == 0 && tag2 == 0) {
        histos.fill(HIST("SE_LL"), dphi1, deta1, pt1, mixpairweight);
        histos.fill(HIST("SE_LL2"), dphi2, deta2, pt2, mixpairweight);
        histos.fill(HIST("hSparseLambdaLambda"), particle1.M(), particle2.M(), cosThetaDiff, deltaR, mixpairweight);
        if (useAdditionalHisto) {
          histos.fill(HIST("hSparseRapLambdaLambda"), particle1.M(), particle2.M(), cosThetaDiff, deltaRap, mixpairweight);
          histos.fill(HIST("hSparsePhiLambdaLambda"), particle1.M(), particle2.M(), cosThetaDiff, dphi_pair, mixpairweight);
          histos.fill(HIST("hSparsePairMassLambdaLambda"), particle1.M(), particle2.M(), cosThetaDiff, pairDummy.M(), mixpairweight);
        }
      } else if (tag1 == 0 && tag2 == 1) {
        histos.fill(HIST("SE_LAL"), dphi1, deta1, pt1, mixpairweight);
        histos.fill(HIST("SE_LAL2"), dphi2, deta2, pt2, mixpairweight);
        histos.fill(HIST("hSparseLambdaAntiLambda"), particle1.M(), particle2.M(), cosThetaDiff, deltaR, mixpairweight);
        if (useAdditionalHisto) {
          histos.fill(HIST("hSparseRapLambdaAntiLambda"), particle1.M(), particle2.M(), cosThetaDiff, deltaRap, mixpairweight);
          histos.fill(HIST("hSparsePhiLambdaAntiLambda"), particle1.M(), particle2.M(), cosThetaDiff, dphi_pair, mixpairweight);
          histos.fill(HIST("hSparsePairMassLambdaAntiLambda"), particle1.M(), particle2.M(), cosThetaDiff, pairDummy.M(), mixpairweight);
        }
      } else if (tag1 == 1 && tag2 == 0) {
        histos.fill(HIST("hSparseAntiLambdaLambda"), particle1.M(), particle2.M(), cosThetaDiff, deltaR, mixpairweight);
        histos.fill(HIST("SE_ALL"), dphi1, deta1, pt1, mixpairweight);
        histos.fill(HIST("SE_ALL2"), dphi2, deta2, pt2, mixpairweight);
        if (useAdditionalHisto) {
          histos.fill(HIST("hSparseRapAntiLambdaLambda"), particle1.M(), particle2.M(), cosThetaDiff, deltaRap, mixpairweight);
          histos.fill(HIST("hSparsePhiAntiLambdaLambda"), particle1.M(), particle2.M(), cosThetaDiff, dphi_pair, mixpairweight);
          histos.fill(HIST("hSparsePairMassAntiLambdaLambda"), particle1.M(), particle2.M(), cosThetaDiff, pairDummy.M(), mixpairweight);
        }
      } else if (tag1 == 1 && tag2 == 1) {
        histos.fill(HIST("hSparseAntiLambdaAntiLambda"), particle1.M(), particle2.M(), cosThetaDiff, deltaR, mixpairweight);
        histos.fill(HIST("SE_ALAL"), dphi1, deta1, pt1, mixpairweight);
        histos.fill(HIST("SE_ALAL2"), dphi2, deta2, pt2, mixpairweight);
        if (useAdditionalHisto) {
          histos.fill(HIST("hSparseRapAntiLambdaAntiLambda"), particle1.M(), particle2.M(), cosThetaDiff, deltaRap, mixpairweight);
          histos.fill(HIST("hSparsePhiAntiLambdaAntiLambda"), particle1.M(), particle2.M(), cosThetaDiff, dphi_pair, mixpairweight);
          histos.fill(HIST("hSparsePairMassAntiLambdaAntiLambda"), particle1.M(), particle2.M(), cosThetaDiff, pairDummy.M(), mixpairweight);
        }
      }
    } else if (datatype == 1) {
      double weight = mixpairweight;
      if (useweight) {
        if (usebothweight) {
          weight = mixpairweight / (epsWeight1 * epsWeight2);
        } else {
          weight = mixpairweight / (epsWeight1);
        }
      }
      if (weight <= 0.0) {
        weight = 1.0;
      }
      histos.fill(HIST("hPtYMix"), particle1.Pt(), particle1.Rapidity(), weight);
      if (tag1 == 0 && tag2 == 0) {
        histos.fill(HIST("ME_LL"), dphi1, deta1, pt1, mixpairweight);
        histos.fill(HIST("ME_LL2"), dphi2, deta2, pt2, mixpairweight);
        histos.fill(HIST("hSparseLambdaLambdaMixed"), particle1.M(), particle2.M(), cosThetaDiff, deltaR, weight);
        if (useAdditionalHisto) {
          histos.fill(HIST("hSparseRapLambdaLambdaMixed"), particle1.M(), particle2.M(), cosThetaDiff, deltaRap, weight);
          histos.fill(HIST("hSparsePhiLambdaLambdaMixed"), particle1.M(), particle2.M(), cosThetaDiff, dphi_pair, weight);
          histos.fill(HIST("hSparsePairMassLambdaLambdaMixed"), particle1.M(), particle2.M(), cosThetaDiff, pairDummy.M(), weight);
        }
      } else if (tag1 == 0 && tag2 == 1) {
        histos.fill(HIST("ME_LAL"), dphi1, deta1, pt1, mixpairweight);
        histos.fill(HIST("ME_LAL2"), dphi2, deta2, pt2, mixpairweight);
        histos.fill(HIST("hSparseLambdaAntiLambdaMixed"), particle1.M(), particle2.M(), cosThetaDiff, deltaR, weight);
        if (useAdditionalHisto) {
          histos.fill(HIST("hSparseRapLambdaAntiLambdaMixed"), particle1.M(), particle2.M(), cosThetaDiff, deltaRap, weight);
          histos.fill(HIST("hSparsePhiLambdaAntiLambdaMixed"), particle1.M(), particle2.M(), cosThetaDiff, dphi_pair, weight);
          histos.fill(HIST("hSparsePairMassLambdaAntiLambdaMixed"), particle1.M(), particle2.M(), cosThetaDiff, pairDummy.M(), weight);
        }
      } else if (tag1 == 1 && tag2 == 0) {
        histos.fill(HIST("ME_ALL"), dphi1, deta1, pt1, mixpairweight);
        histos.fill(HIST("ME_ALL2"), dphi2, deta2, pt2, mixpairweight);
        histos.fill(HIST("hSparseAntiLambdaLambdaMixed"), particle1.M(), particle2.M(), cosThetaDiff, deltaR, weight);
        if (useAdditionalHisto) {
          histos.fill(HIST("hSparseRapAntiLambdaLambdaMixed"), particle1.M(), particle2.M(), cosThetaDiff, deltaRap, weight);
          histos.fill(HIST("hSparsePhiAntiLambdaLambdaMixed"), particle1.M(), particle2.M(), cosThetaDiff, dphi_pair, weight);
          histos.fill(HIST("hSparsePairMassAntiLambdaLambdaMixed"), particle1.M(), particle2.M(), cosThetaDiff, pairDummy.M(), weight);
        }
      } else if (tag1 == 1 && tag2 == 1) {
        histos.fill(HIST("ME_ALAL"), dphi1, deta1, pt1, mixpairweight);
        histos.fill(HIST("ME_ALAL2"), dphi2, deta2, pt2, mixpairweight);
        histos.fill(HIST("hSparseAntiLambdaAntiLambdaMixed"), particle1.M(), particle2.M(), cosThetaDiff, deltaR, weight);
        if (useAdditionalHisto) {
          histos.fill(HIST("hSparseRapAntiLambdaAntiLambdaMixed"), particle1.M(), particle2.M(), cosThetaDiff, deltaRap, weight);
          histos.fill(HIST("hSparsePhiAntiLambdaAntiLambdaMixed"), particle1.M(), particle2.M(), cosThetaDiff, dphi_pair, weight);
          histos.fill(HIST("hSparsePairMassAntiLambdaAntiLambdaMixed"), particle1.M(), particle2.M(), cosThetaDiff, pairDummy.M(), weight);
        }
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
        if (v0.protonIndex() == v02.pionIndex()) {
          continue;
        }
        if (v0.pionIndex() == v02.protonIndex()) {
          continue;
        }
        proton2 = ROOT::Math::PtEtaPhiMVector(v02.protonPt(), v02.protonEta(), v02.protonPhi(), o2::constants::physics::MassProton);
        lambda2 = ROOT::Math::PtEtaPhiMVector(v02.lambdaPt(), v02.lambdaEta(), v02.lambdaPhi(), v02.lambdaMass());
        histos.fill(HIST("deltaPhiSame"), RecoDecay::constrainAngle(v0.lambdaPhi() - v02.lambdaPhi(), 0.0F, harmonicDphi));
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
          histos.fill(HIST("deltaPhiMix"), RecoDecay::constrainAngle(t3.lambdaPhi() - t2.lambdaPhi(), 0.0F, harmonicDphi));
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
    auto nBins = colBinning.getAllBinsCount();
    std::vector<std::deque<std::pair<int, AllTrackCandidates>>> eventPools(nBins);

    for (auto& collision1 : collisions) {
      const int bin = colBinning.getBin(std::make_tuple(collision1.posz(), collision1.cent()));

      // if pool empty, push and continue
      if (eventPools[bin].empty()) {
        auto sliced = V0s.sliceBy(tracksPerCollisionV0, collision1.index());
        eventPools[bin].emplace_back(collision1.index(), std::move(sliced));
        if ((int)eventPools[bin].size() > nEvtMixing)
          eventPools[bin].pop_front();
        continue;
      }

      // current event slice
      auto poolA = V0s.sliceBy(tracksPerCollisionV0, collision1.index());

      // loop over SE unordered pairs (t1,t2)
      for (auto& [t1, t2] : soa::combinations(o2::soa::CombinationsFullIndexPolicy(poolA, poolA))) {
        if (!selectionV0(t1) || !selectionV0(t2))
          continue;
        if (t2.index() <= t1.index())
          continue;
        if (t1.protonIndex() == t2.protonIndex())
          continue;
        if (t1.pionIndex() == t2.pionIndex())
          continue;
        if (t1.protonIndex() == t2.pionIndex())
          continue;
        if (t1.pionIndex() == t2.protonIndex())
          continue;

        // scan prior events for replacements for t1
        struct PV {
          AllTrackCandidates* pool;
          int nRepl;
        };
        std::vector<PV> usable;
        int totalRepl = 0;

        int mixes = 0;
        for (auto it = eventPools[bin].rbegin();
             it != eventPools[bin].rend() && mixes < nEvtMixing; ++it, ++mixes) {
          const int collision2idx = it->first;
          auto& poolB = it->second;
          if (collision2idx == collision1.index())
            continue;

          int nRepl = 0;
          for (auto& tX : poolB) {
            if (!selectionV0(tX))
              continue;
            if (checkKinematics(t1, tX))
              ++nRepl;
          }
          if (nRepl > 0) {
            usable.push_back(PV{&poolB, nRepl});
            totalRepl += nRepl;
          }
        }

        if (totalRepl == 0)
          continue;
        const float wBase = 1.0f / static_cast<float>(totalRepl);

        // emit mixed pairs: tX replaces t1; t2 stays
        for (auto& pv : usable) {
          auto& poolB = *pv.pool;
          for (auto& tX : poolB) {
            if (!selectionV0(tX))
              continue;
            if (!checkKinematics(t1, tX))
              continue;

            auto proton = ROOT::Math::PtEtaPhiMVector(tX.protonPt(), tX.protonEta(), tX.protonPhi(), o2::constants::physics::MassProton);
            auto lambda = ROOT::Math::PtEtaPhiMVector(tX.lambdaPt(), tX.lambdaEta(), tX.lambdaPhi(), tX.lambdaMass());
            auto proton2 = ROOT::Math::PtEtaPhiMVector(t2.protonPt(), t2.protonEta(), t2.protonPhi(), o2::constants::physics::MassProton);
            auto lambda2 = ROOT::Math::PtEtaPhiMVector(t2.lambdaPt(), t2.lambdaEta(), t2.lambdaPhi(), t2.lambdaMass());

            const float dPhi = std::fabs(RecoDecay::constrainAngle(lambda.Phi(), 0.0F, harmonic) - RecoDecay::constrainAngle(lambda2.Phi(), 0.0F, harmonic));
            histos.fill(HIST("deltaPhiMix"), dPhi, wBase);
            fillHistograms(tX.v0Status(), t2.v0Status(), lambda, lambda2, proton, proton2, 1, wBase);
          }
        }
      }
      // push current event into pool
      auto sliced = V0s.sliceBy(tracksPerCollisionV0, collision1.index());
      eventPools[bin].emplace_back(collision1.index(), std::move(sliced));
      if ((int)eventPools[bin].size() > nEvtMixing)
        eventPools[bin].pop_front();
    }
  }
  PROCESS_SWITCH(lambdaspincorrderived, processMEV3, "Process data ME (first-leg, pair-3D maps)", false);

  static constexpr int N_STATUS = 2; // v0Status ∈ {0,1}

  struct MixBinner {
    // constructed from the task's configurables; φ is assumed already constrained into [0, 2π)
    float ptMin, ptMax, ptStep;
    float etaMin, etaMax, etaStep;
    float phiMin, phiMax, phiStep;

    // Mass binning: [1.09, 1.14) with 50 bins (1e-3 GeV/c^2)
    static constexpr float mMin = 1.09f;
    static constexpr float mMax = 1.14f; // exclusive
    static constexpr int nM_ = 1;
    static constexpr float mStep = (mMax - mMin) / nM_;

    int nPt_, nEta_, nPhi_;

    MixBinner(float ptMin_, float ptMax_, float ptStep_,
              float etaAbsMax, float etaStep_,
              float phiStep_)
      : ptMin(ptMin_), ptMax(ptMax_), ptStep(ptStep_), etaMin(-etaAbsMax), etaMax(+etaAbsMax), etaStep(etaStep_), phiMin(0.f), phiMax(static_cast<float>(2.0 * TMath::Pi())), phiStep(phiStep_)
    {
      ptStep = (ptStep > 0.f ? ptStep : 0.1f);
      etaStep = (etaStep > 0.f ? etaStep : 0.1f);
      phiStep = (phiStep > 0.f ? phiStep : 0.1f);

      nPt_ = std::max(1, static_cast<int>(std::floor((ptMax - ptMin) / ptStep + 0.5f)));
      nEta_ = std::max(1, static_cast<int>(std::floor((etaMax - etaMin) / etaStep + 0.5f)));
      nPhi_ = std::max(1, static_cast<int>(std::ceil((phiMax - phiMin) / phiStep)));
    }

    inline int nPt() const { return nPt_; }
    inline int nEta() const { return nEta_; }
    inline int nPhi() const { return nPhi_; }
    inline int nM() const { return nM_; }

    inline int binFromValue(float v, float vmin, float step, int nBins) const
    {
      if (!std::isfinite(v))
        return -1;
      const float x = (v - vmin) / step;
      int b = static_cast<int>(std::floor(x + 1e-6f));
      if (b < 0)
        return -1;
      if (b >= nBins)
        b = nBins - 1; // clamp exact-top edge
      return b;
    }

    inline int ptBin(float pt) const { return binFromValue(pt, ptMin, ptStep, nPt_); }
    inline int etaBin(float eta) const { return binFromValue(eta, etaMin, etaStep, nEta_); }
    inline int phiBin(float phi) const { return binFromValue(phi, phiMin, phiStep, nPhi_); } // φ already constrained upstream
    inline int massBin(float m) const { return binFromValue(m, mMin, mStep, nM_); }
  };

  struct BufferCand {
    int64_t collisionIdx; // from col.index()
    int64_t rowIndex;     // global row id in V0s
    uint8_t v0Status;
    uint16_t ptBin, etaBin, phiBin, mBin;
  };

  struct MatchRef {
    int64_t collisionIdx;
    int64_t rowIndex;
  };

  // 6D key: (colBin, status, pt, eta, phi, mass)
  static inline size_t linearKey(int colBin, int statBin,
                                 int ptBin, int etaBin, int phiBin, int mBin,
                                 int nStatus, int nPt, int nEta, int nPhi, int nM)
  {
    return ((((((static_cast<size_t>(colBin) * nStatus + statBin) * nPt + ptBin) * nEta + etaBin) * nPhi + phiBin) * nM + mBin));
  }

  // ===================== Main mixing (with mass-bin + random unique sampling) =====================
  void processMEV4(EventCandidates const& collisions, AllTrackCandidates const& V0s)
  {
    // Build binner from your existing configurables
    MixBinner mb{
      ptMin.value, ptMax.value, ptMix.value, // pT range & step
      v0eta.value, etaMix.value,             // |eta| max & step
      phiMix.value                           // φ step; φ range fixed to [0, 2π)
    };

    const int nCol = colBinning.getAllBinsCount(); // event-class bins (vz, centrality)
    const int nStat = N_STATUS;                    // 2
    const int nPt = mb.nPt();
    const int nEta = mb.nEta();
    const int nPhi = mb.nPhi();
    const int nM = mb.nM();

    const size_t nKeys = static_cast<size_t>(nCol) * nStat * nPt * nEta * nPhi * nM;
    std::vector<std::vector<BufferCand>> buffer(nKeys);

    // ---- PASS 1: fill 6D buffer ----
    for (auto const& col : collisions) {
      const int colBin = colBinning.getBin(std::make_tuple(col.posz(), col.cent()));
      auto slice = V0s.sliceBy(tracksPerCollisionV0, col.index());

      for (auto const& t : slice) {
        if (!selectionV0(t))
          continue;

        const int status = static_cast<int>(t.v0Status());
        if (status < 0 || status >= nStat)
          continue;

        // Bin kinematics (φ already constrained via your call-site)
        const int ptB = mb.ptBin(t.lambdaPt());
        const int etaB = mb.etaBin(t.lambdaEta());
        const int phiB = mb.phiBin(RecoDecay::constrainAngle(t.lambdaPhi(), 0.0F, harmonic));
        const int mB = mb.massBin(t.lambdaMass());
        if (ptB < 0 || etaB < 0 || phiB < 0 || mB < 0)
          continue;

        const size_t key = linearKey(colBin, status, ptB, etaB, phiB, mB,
                                     nStat, nPt, nEta, nPhi, nM);

        buffer[key].push_back(BufferCand{
          .collisionIdx = static_cast<int64_t>(col.index()),
          .rowIndex = static_cast<int64_t>(t.globalIndex()), // adapt accessor if needed
          .v0Status = static_cast<uint8_t>(status),
          .ptBin = static_cast<uint16_t>(ptB),
          .etaBin = static_cast<uint16_t>(etaB),
          .phiBin = static_cast<uint16_t>(phiB),
          .mBin = static_cast<uint16_t>(mB)});
      }
    }

    // ---- PASS 2: mixing over same-event pairs ----
    for (auto const& collision1 : collisions) {
      const int colBin = colBinning.getBin(std::make_tuple(collision1.posz(), collision1.cent()));
      auto poolA = V0s.sliceBy(tracksPerCollisionV0, collision1.index());

      for (auto const& [t1, t2] :
           soa::combinations(o2::soa::CombinationsFullIndexPolicy(poolA, poolA))) {

        if (!selectionV0(t1) || !selectionV0(t2))
          continue;
        if (t2.index() <= t1.index())
          continue;

        // no shared daughters
        if (t1.protonIndex() == t2.protonIndex())
          continue;
        if (t1.pionIndex() == t2.pionIndex())
          continue;
        if (t1.protonIndex() == t2.pionIndex())
          continue;
        if (t1.pionIndex() == t2.protonIndex())
          continue;

        const int status = static_cast<int>(t1.v0Status());
        if (status < 0 || status >= nStat)
          continue;

        // Bin of t1 defines where to search (exact 6D bin)
        const int ptB = mb.ptBin(t1.lambdaPt());
        const int etaB = mb.etaBin(t1.lambdaEta());
        const int phiB = mb.phiBin(RecoDecay::constrainAngle(t1.lambdaPhi(), 0.0F, harmonic)); // φ already constrained upstream
        const int mB = mb.massBin(t1.lambdaMass());
        if (ptB < 0 || etaB < 0 || phiB < 0 || mB < 0)
          continue;

        const size_t key = linearKey(colBin, status, ptB, etaB, phiB, mB,
                                     nStat, nPt, nEta, nPhi, nM);
        auto const& binVec = buffer[key];
        if (binVec.empty())
          continue;

        // Collect all partners from this 6D bin but different collision
        std::vector<MatchRef> matches;
        matches.reserve(binVec.size());
        const int64_t curColIdx = static_cast<int64_t>(collision1.index());

        for (const auto& bc : binVec) {
          if (bc.collisionIdx == curColIdx)
            continue; // ensure different event
          matches.push_back(MatchRef{bc.collisionIdx, bc.rowIndex});
        }
        if (matches.empty())
          continue;

        // ---------- YOUR PREFERRED RANDOM UNIQUE SAMPLING BLOCK ----------
        const int cap = nEvtMixing.value;
        const int n = static_cast<int>(matches.size());
        if (cap > 0 && cap < n) {
          std::uniform_int_distribution<int> dist(0, n - 1);
          // pick cap unique indices
          std::unordered_set<int> chosen;
          chosen.reserve(cap * 2);
          while ((int)chosen.size() < cap) {
            chosen.insert(dist(rng));
          }
          std::vector<MatchRef> subset;
          subset.reserve(cap);
          for (int idx : chosen)
            subset.push_back(matches[idx]);
          matches.swap(subset);
        } else {
          std::shuffle(matches.begin(), matches.end(), rng);
        }
        // ----------------------------------------------------------------

        const float wBase = 1.0f / static_cast<float>(matches.size());

        // Emit mixed pairs: tX replaces t1; keep t2
        for (const auto& m : matches) {
          auto tX = V0s.iteratorAt(m.rowIndex); // replace accessor if different
          if (!selectionV0(tX))
            continue; // optional extra guard
          if (!checkKinematics(t1, tX))
            continue;

          auto proton = ROOT::Math::PtEtaPhiMVector(tX.protonPt(), tX.protonEta(), tX.protonPhi(), o2::constants::physics::MassProton);
          auto lambda = ROOT::Math::PtEtaPhiMVector(tX.lambdaPt(), tX.lambdaEta(), tX.lambdaPhi(), tX.lambdaMass());
          auto proton2 = ROOT::Math::PtEtaPhiMVector(t2.protonPt(), t2.protonEta(), t2.protonPhi(), o2::constants::physics::MassProton);
          auto lambda2 = ROOT::Math::PtEtaPhiMVector(t2.lambdaPt(), t2.lambdaEta(), t2.lambdaPhi(), t2.lambdaMass());

          const float dPhi = std::fabs(RecoDecay::constrainAngle(lambda.Phi() - lambda2.Phi(), 0.0F, harmonicDphi));
          histos.fill(HIST("deltaPhiMix"), dPhi, wBase);
          fillHistograms(tX.v0Status(), t2.v0Status(), lambda, lambda2, proton, proton2, 1, wBase);
        }
      }
    }
  }
  PROCESS_SWITCH(lambdaspincorrderived, processMEV4, "Process data ME (5d buffer)", false);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<lambdaspincorrderived>(cfgc)};
}
