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
namespace mcacc
{
// event
template <typename Coll>
static inline float cent(const Coll& c)
{
  return c.centmc();
}

template <typename Coll>
static inline float posz(const Coll& c)
{
  return c.poszmc();
}

// pair / v0 candidate
template <typename T>
static inline int v0Status(const T& t)
{
  return t.v0Statusmc();
}

template <typename T>
static inline bool doubleStatus(const T& t)
{
  return t.doubleStatusmc();
}

template <typename T>
static inline float v0CosPA(const T& t)
{
  return t.v0Cospamc();
}

template <typename T>
static inline float v0Radius(const T& t)
{
  return t.v0Radiusmc();
}

template <typename T>
static inline float dcaPos(const T& t)
{
  return t.dcaPositivemc();
}

template <typename T>
static inline float dcaNeg(const T& t)
{
  return t.dcaNegativemc();
}

template <typename T>
static inline float dcaDau(const T& t)
{
  return t.dcaBetweenDaughtermc();
}

template <typename T>
static inline float lamPt(const T& t)
{
  return t.lambdaPtmc();
}

template <typename T>
static inline float lamEta(const T& t)
{
  return t.lambdaEtamc();
}

template <typename T>
static inline float lamPhi(const T& t)
{
  return t.lambdaPhimc();
}

template <typename T>
static inline float lamMass(const T& t)
{
  return t.lambdaMassmc();
}

template <typename T>
static inline float prPt(const T& t)
{
  return t.protonPtmc();
}

template <typename T>
static inline float prEta(const T& t)
{
  return t.protonEtamc();
}

template <typename T>
static inline float prPhi(const T& t)
{
  return t.protonPhimc();
}

template <typename T>
static inline int prIdx(const T& t)
{
  return t.protonIndexmc();
}

template <typename T>
static inline int piIdx(const T& t)
{
  return t.pionIndexmc();
}
} // namespace mcacc

struct lambdaspincorrderived {
  // BinningType colBinning;
  struct : ConfigurableGroup {
    Configurable<std::string> cfgURL{"cfgURL", "http://alice-ccdb.cern.ch", "Address of the CCDB to browse"};
    Configurable<int64_t> nolaterthan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "Latest acceptable timestamp of creation for the object"};
  } cfgCcdbParam;

  struct : ConfigurableGroup {
    ConfigurableAxis cfgMixRadiusBins{"cfgMixRadiusBins", {VARIABLE_WIDTH, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0, 15.0, 20.0, 25.0, 30.0}, "Radius bins for V6 radius buffer"};
  } cfgMixRadiusParam;

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
  TH2D* hweightCentPair = nullptr;

  Configurable<std::string> ConfWeightPathCentPair{"ConfWeightPathCentPair", "", "Centrality x pair-type weight path"};
  Configurable<std::string> ConfWeightPathLL{"ConfWeightPathLL", "Users/s/skundu/My/Object/spincorr/cent010LL", "Weight path"};
  Configurable<std::string> ConfWeightPathALAL{"ConfWeightPathALAL", "Users/s/skundu/My/Object/spincorr/cent010LL", "Weight path"};
  Configurable<std::string> ConfWeightPathLAL{"ConfWeightPathLAL", "Users/s/skundu/My/Object/spincorr/cent010LL", "Weight path"};
  Configurable<std::string> ConfWeightPathALL{"ConfWeightPathALL", "Users/s/skundu/My/Object/spincorr/cent010LL", "Weight path"};

  Configurable<std::string> ConfWeightPathLL2{"ConfWeightPathLL2", "Users/s/skundu/My/Object/spincorr/cent010LL", "Weight path 2"};
  Configurable<std::string> ConfWeightPathALAL2{"ConfWeightPathALAL2", "Users/s/skundu/My/Object/spincorr/cent010LL", "Weight path 2"};
  Configurable<std::string> ConfWeightPathLAL2{"ConfWeightPathLAL2", "Users/s/skundu/My/Object/spincorr/cent010LL", "Weight path 2"};
  Configurable<std::string> ConfWeightPathALL2{"ConfWeightPathALL2", "Users/s/skundu/My/Object/spincorr/cent010LL", "Weight path 2"};

  // Mixing /////////
  struct : ConfigurableGroup {
    Configurable<double> nKinematicPt{"nKinematicPt", 1.0, "Number of pT buffer bins"};
    Configurable<double> nKinematicEta{"nKinematicEta", 1.0, "Number of eta buffer bins"};
    Configurable<double> nKinematicPhi{"nKinematicPhi", 1.0, "Number of phi buffer bins"};
  } cfgKinematicBins;

  Configurable<float> ptMinMixBuffer{"ptMinMixBuffer", 0.7, "Minimum V0 pT for mix buffer"};
  Configurable<float> ptMaxMixBuffer{"ptMaxMixBuffer", 4.1, "Maximum V0 pT for mix buffer"};
  Configurable<int> cfgMixLegMode{"cfgMixLegMode", 0, "0=replace leg-1 only, 1=replace leg-2 only, 2=do both one-leg replacements"};
  Configurable<int> cfgV5MassBins{"cfgV5MassBins", 5, "Number of fixed mass bins for V5 mixing"};
  Configurable<int> cfgV5NeighborPt{"cfgV5NeighborPt", 0, "v5: neighbor bins in pT (use symmetric ±N, edge-safe)"};
  Configurable<int> cfgV5NeighborEta{"cfgV5NeighborEta", 0, "v5: neighbor bins in eta (use symmetric ±N, edge-safe)"};
  Configurable<int> cfgV5NeighborPhi{"cfgV5NeighborPhi", 0, "v5: neighbor bins in phi (use symmetric ±N, periodic wrap)"};
  Configurable<bool> usePairKineMatch{"usePairKineMatch", true, "Require pair-level matching between (A,B) and (C,B)"};
  Configurable<int> cfgV5MaxMatches{"cfgV5MaxMatches", 50, "v5: max ME replacements per SE pair (after all cuts)"};
  Configurable<uint64_t> cfgMixSeed{"cfgMixSeed", 0xdecafbadULL, "RNG seed for downsampling matches (deterministic)"};
  Configurable<float> centMin{"centMin", 0, "Minimum Centrality"};
  Configurable<float> centMax{"centMax", 80, "Maximum Centrality"};
  Configurable<int> rngSeed{"rngSeed", 12345, "Seed for random mixing (reproducible)"};
  std::mt19937 rng{12345};
  Configurable<int> nEvtMixing{"nEvtMixing", 10, "Number of events to mix"};
  ConfigurableAxis CfgVtxBins{"CfgVtxBins", {VARIABLE_WIDTH, -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10}, "Mixing bins - z-vertex"};
  ConfigurableAxis CfgMultBins{"CfgMultBins", {VARIABLE_WIDTH, 0, 110}, "Mixing bins - centrality"};
  Configurable<float> etaMix{"etaMix", 0.1, "Eta cut on event mixing"};
  Configurable<float> ptMix{"ptMix", 0.1, "Pt cut on event mixing"};
  Configurable<float> phiMix{"phiMix", 0.1, "Phi cut on event mixing"};
  Configurable<float> massMix{"massMix", 0.0028, "Masscut on event mixing"};
  Configurable<bool> userapidity{"userapidity", 1, "Use Rapidity for mixing"};

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
  Configurable<float> v0etaMixBuffer{"v0etaMixBuffer", 0.8, "Eta cut on mix event buffer"};
  Configurable<float> v0eta{"v0eta", 0.8, "Eta cut on lambda"};

  // Event Mixing
  Configurable<int> cosDef{"cosDef", 1, "Defination of cos"};

  ConfigurableAxis ax_dphi_h{"ax_dphi_h", {VARIABLE_WIDTH, 0.0, 2.0 * TMath::Pi()}, "Δφ_h"};
  ConfigurableAxis ax_deta{"ax_deta", {VARIABLE_WIDTH, -1.0, 1.0}, "Δη"};
  ConfigurableAxis ax_ptpair{"ax_ptpair", {VARIABLE_WIDTH, 0.0, 10.0}, "p_{T,pair} (GeV/c)"};

  // THnsparse bining
  ConfigurableAxis configThnAxisInvMass{"configThnAxisInvMass", {50, 1.09, 1.14}, "#it{M} (GeV/#it{c}^{2})"};
  ConfigurableAxis configThnAxisR{"configThnAxisR", {VARIABLE_WIDTH, 0.0, 8.0}, "#it{R}"};
  ConfigurableAxis configThnAxisPol{"configThnAxisPol", {VARIABLE_WIDTH, 0.0, 8.0}, "cos#it{#theta *}"};
  ConfigurableAxis configThnAxisCentrality{"configThnAxisCentrality", {VARIABLE_WIDTH, 0.0, 80.0}, "Centrality"};
  ConfigurableAxis configThnAxisRapidity{"configThnAxisRapidity", {VARIABLE_WIDTH, 0.0, 1.0}, "Rapidity"};
  ConfigurableAxis configThnAxisPairMass{"configThnAxisPairMass", {VARIABLE_WIDTH, 2.0, 3.0}, "PairMass"};
  ConfigurableAxis configThnAxisPhi{"configThnAxisPhi", {VARIABLE_WIDTH, 0.0, 2.0 * TMath::Pi()}, "Phi"};

  ConfigurableAxis configThnAxisDeltaPhi{"configThnAxisDeltaPhi", {VARIABLE_WIDTH, 0.0, TMath::Pi() / 6, 2.0 * TMath::Pi() / 6, 3.0 * TMath::Pi() / 6, 4.0 * TMath::Pi() / 6, 5.0 * TMath::Pi() / 6, TMath::Pi()}, "Delta Phi"};
  ConfigurableAxis configThnAxisDeltaR{"configThnAxisDeltaR", {VARIABLE_WIDTH, 0.0, 0.5, 1.2, 2.0, 3.1, 4.0}, "Delta R"};
  ConfigurableAxis configThnAxisDeltaRap{"configThnAxisDeltaRap", {VARIABLE_WIDTH, 0.0, 0.2, 0.5, 1.0, 1.6}, "Delta Rap"};

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext&)
  {
    histos.add("hPtRadiusV0", "V0 QA;#it{p}_{T}^{V0} (GeV/#it{c});V0 decay radius (cm)", kTH2F, {{100, 0.0, 10.0}, {120, 0.0, 45.0}});
    histos.add("hPtYSame", "hPtYSame", kTH2F, {{100, 0.0, 10.0}, {200, -1.0, 1.0}});
    histos.add("hPtYMix", "hPtYMix", kTH2F, {{100, 0.0, 10.0}, {200, -1.0, 1.0}});
    histos.add("hCentrality", "Centrality distribution", kTH1F, {{configThnAxisCentrality}});
    histos.add("deltaPhiSame", "deltaPhiSame", HistType::kTH1D, {{72, -TMath::Pi(), TMath::Pi()}}, true);
    histos.add("deltaPhiMix", "deltaPhiMix", HistType::kTH1D, {{72, -TMath::Pi(), TMath::Pi()}}, true);
    histos.add("ptCent", "ptCent", HistType::kTH2D, {{100, 0.0, 10.0}, {8, 0.0, 80.0}}, true);
    histos.add("etaCent", "etaCent", HistType::kTH2D, {{32, -0.8, 0.8}, {8, 0.0, 80.0}}, true);

    histos.add("hCentPairTypeSE", "SE pair-weighted centrality;Centrality;PairType", kTH2D, {{110, 0.0, 110.0}, {4, -0.5, 3.5}});
    histos.add("hCentPairTypeME", "ME pair-weighted centrality;Centrality;PairType", kTH2D, {{110, 0.0, 110.0}, {4, -0.5, 3.5}});

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

    histos.add("hSparseLambdaLambdaAnalysis", "hSparseLambdaLambdaAnalysis", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisInvMass, configThnAxisPol, configThnAxisDeltaR, configThnAxisDeltaRap, configThnAxisDeltaPhi}, true);
    histos.add("hSparseLambdaAntiLambdaAnalysis", "hSparseLambdaAntiLambdaAnalysis", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisInvMass, configThnAxisPol, configThnAxisDeltaR, configThnAxisDeltaRap, configThnAxisDeltaPhi}, true);
    histos.add("hSparseAntiLambdaLambdaAnalysis", "hSparseAntiLambdLambdaAnalysis", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisInvMass, configThnAxisPol, configThnAxisDeltaR, configThnAxisDeltaRap, configThnAxisDeltaPhi}, true);
    histos.add("hSparseAntiLambdaAntiLambdaAnalysis", "hSparseAntiLambdaAntiLambdaAnalysis", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisInvMass, configThnAxisPol, configThnAxisDeltaR, configThnAxisDeltaRap, configThnAxisDeltaPhi}, true);

    histos.add("hSparseLambdaLambdaMixedAnalysis", "hSparseLambdaLambdaMixedAnalysis", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisInvMass, configThnAxisPol, configThnAxisDeltaR, configThnAxisDeltaRap, configThnAxisDeltaPhi}, true);
    histos.add("hSparseLambdaAntiLambdaMixedAnalysis", "hSparseLambdaAntiLambdaMixedAnalysis", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisInvMass, configThnAxisPol, configThnAxisDeltaR, configThnAxisDeltaRap, configThnAxisDeltaPhi}, true);
    histos.add("hSparseAntiLambdaLambdaMixedAnalysis", "hSparseAntiLambdaLambdaMixedAnalysis", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisInvMass, configThnAxisPol, configThnAxisDeltaR, configThnAxisDeltaRap, configThnAxisDeltaPhi}, true);
    histos.add("hSparseAntiLambdaAntiLambdaMixedAnalysis", "hSparseAntiLambdaAntiLambdaMixedAnalysis", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisInvMass, configThnAxisPol, configThnAxisDeltaR, configThnAxisDeltaRap, configThnAxisDeltaPhi}, true);

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
    if (!ConfWeightPathCentPair.value.empty()) {
      hweightCentPair = ccdb->getForTimeStamp<TH2D>(ConfWeightPathCentPair.value, cfgCcdbParam.nolaterthan.value);
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
    if (candidate.v0Status() == 0 && (std::abs(candidate.dcaPositive()) < dcaProton || std::abs(candidate.dcaNegative()) < dcaPion)) {
      return false;
    }
    if (candidate.v0Status() == 1 && (std::abs(candidate.dcaPositive()) < dcaPion || std::abs(candidate.dcaNegative()) < dcaProton)) {
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

  template <typename T>
  bool selectionV0Buffer(T const& candidate)
  {
    auto particle = ROOT::Math::PtEtaPhiMVector(candidate.lambdaPt(), candidate.lambdaEta(), candidate.lambdaPhi(), candidate.lambdaMass());

    if (std::abs(particle.Rapidity()) > rapidity || std::abs(particle.Eta()) > v0etaMixBuffer) {
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
    if (candidate.v0Status() == 0 && (std::abs(candidate.dcaPositive()) < dcaProton || std::abs(candidate.dcaNegative()) < dcaPion)) {
      return false;
    }
    if (candidate.v0Status() == 1 && (std::abs(candidate.dcaPositive()) < dcaPion || std::abs(candidate.dcaNegative()) < dcaProton)) {
      return false;
    }
    if (candidate.lambdaPt() < ptMinMixBuffer) {
      return false;
    }
    if (candidate.lambdaPt() > ptMaxMixBuffer) {
      return false;
    }
    return true;
  }

  template <typename T>
  bool selectionV0BufferMC(T const& candidate)
  {
    auto particle = ROOT::Math::PtEtaPhiMVector(mcacc::lamPt(candidate),
                                                mcacc::lamEta(candidate),
                                                mcacc::lamPhi(candidate),
                                                mcacc::lamMass(candidate));

    if (std::abs(particle.Rapidity()) > rapidity || std::abs(particle.Eta()) > v0etaMixBuffer) {
      return false;
    }
    if (mcacc::lamMass(candidate) < MassMin || mcacc::lamMass(candidate) > MassMax) {
      return false;
    }
    if (mcacc::v0CosPA(candidate) < cosPA) {
      return false;
    }
    if (checkDoubleStatus && mcacc::doubleStatus(candidate)) {
      return false;
    }
    if (mcacc::v0Radius(candidate) > radiusMax) {
      return false;
    }
    if (mcacc::v0Radius(candidate) < radiusMin) {
      return false;
    }
    if (mcacc::dcaDau(candidate) > dcaDaughters) {
      return false;
    }
    if (mcacc::v0Status(candidate) == 0 && (std::abs(mcacc::dcaPos(candidate)) < dcaProton || std::abs(mcacc::dcaNeg(candidate)) < dcaPion)) {
      return false;
    }
    if (mcacc::v0Status(candidate) == 1 && (std::abs(mcacc::dcaPos(candidate)) < dcaPion || std::abs(mcacc::dcaNeg(candidate)) < dcaProton)) {
      return false;
    }
    if (mcacc::lamPt(candidate) < ptMinMixBuffer) {
      return false;
    }
    if (mcacc::lamPt(candidate) > ptMaxMixBuffer) {
      return false;
    }
    return true;
  }

  template <typename T1, typename T2>
  bool checkKinematics(T1 const& c1, T2 const& c2)
  {
    if (c1.v0Status() != c2.v0Status()) {
      return false;
    }

    if (std::abs(c1.lambdaPt() - c2.lambdaPt()) > ptMix) {
      return false;
    }

    if (!userapidity) {
      if (std::abs(c1.lambdaEta() - c2.lambdaEta()) > etaMix) {
        return false;
      }
    } else {
      const auto l1 = ROOT::Math::PtEtaPhiMVector(c1.lambdaPt(), c1.lambdaEta(), c1.lambdaPhi(), c1.lambdaMass());
      const auto l2 = ROOT::Math::PtEtaPhiMVector(c2.lambdaPt(), c2.lambdaEta(), c2.lambdaPhi(), c2.lambdaMass());
      if (std::abs(l1.Rapidity() - l2.Rapidity()) > etaMix) { // etaMix used as Δy
        return false;
      }
    }

    const float dphi = deltaPhiMinusPiToPi((float)c1.lambdaPhi(), (float)c2.lambdaPhi());
    if (std::abs(dphi) > phiMix) {
      return false;
    }

    if (std::abs(c1.lambdaMass() - c2.lambdaMass()) > massMix) {
      return false;
    }

    return true;
  }

  template <typename TA, typename TB, typename TC>
  bool checkPairKinematics(TA const& A, TB const& B, TC const& C)
  {
    if (!usePairKineMatch) {
      return true;
    }

    const auto lA = ROOT::Math::PtEtaPhiMVector(A.lambdaPt(), A.lambdaEta(), A.lambdaPhi(), A.lambdaMass());
    const auto lB = ROOT::Math::PtEtaPhiMVector(B.lambdaPt(), B.lambdaEta(), B.lambdaPhi(), B.lambdaMass());
    const auto lC = ROOT::Math::PtEtaPhiMVector(C.lambdaPt(), C.lambdaEta(), C.lambdaPhi(), C.lambdaMass());

    // relative pT inside the pair: |pT1 - pT2|
    const float dPtAB = std::abs(A.lambdaPt() - B.lambdaPt());
    const float dPtCB = std::abs(C.lambdaPt() - B.lambdaPt());
    if (std::abs(dPtAB - dPtCB) > ptMix) {
      return false;
    }

    // relative longitudinal kinematics: |Δy| or |Δη|
    if (userapidity) {
      const float dYAB = std::abs(lA.Rapidity() - lB.Rapidity());
      const float dYCB = std::abs(lC.Rapidity() - lB.Rapidity());
      if (std::abs(dYAB - dYCB) > etaMix) {
        return false;
      }
    } else {
      const float dEtaAB = std::abs(A.lambdaEta() - B.lambdaEta());
      const float dEtaCB = std::abs(C.lambdaEta() - B.lambdaEta());
      if (std::abs(dEtaAB - dEtaCB) > etaMix) {
        return false;
      }
    }

    // relative azimuth inside the pair: |Δφ|
    const float dPhiAB = std::abs(deltaPhiMinusPiToPi((float)A.lambdaPhi(), (float)B.lambdaPhi()));
    const float dPhiCB = std::abs(deltaPhiMinusPiToPi((float)C.lambdaPhi(), (float)B.lambdaPhi()));
    if (std::abs(dPhiAB - dPhiCB) > phiMix) {
      return false;
    }

    return true;
  }
  template <typename TA, typename TB, typename TC>
  bool checkPairKinematicsMC(TA const& A, TB const& B, TC const& C)
  {
    if (!usePairKineMatch) {
      return true;
    }

    const auto lA = ROOT::Math::PtEtaPhiMVector(mcacc::lamPt(A), mcacc::lamEta(A), mcacc::lamPhi(A), mcacc::lamMass(A));
    const auto lB = ROOT::Math::PtEtaPhiMVector(mcacc::lamPt(B), mcacc::lamEta(B), mcacc::lamPhi(B), mcacc::lamMass(B));
    const auto lC = ROOT::Math::PtEtaPhiMVector(mcacc::lamPt(C), mcacc::lamEta(C), mcacc::lamPhi(C), mcacc::lamMass(C));

    const float dPtAB = std::abs(mcacc::lamPt(A) - mcacc::lamPt(B));
    const float dPtCB = std::abs(mcacc::lamPt(C) - mcacc::lamPt(B));
    if (std::abs(dPtAB - dPtCB) > ptMix) {
      return false;
    }

    if (userapidity) {
      const float dYAB = std::abs(lA.Rapidity() - lB.Rapidity());
      const float dYCB = std::abs(lC.Rapidity() - lB.Rapidity());
      if (std::abs(dYAB - dYCB) > etaMix) {
        return false;
      }
    } else {
      const float dEtaAB = std::abs(mcacc::lamEta(A) - mcacc::lamEta(B));
      const float dEtaCB = std::abs(mcacc::lamEta(C) - mcacc::lamEta(B));
      if (std::abs(dEtaAB - dEtaCB) > etaMix) {
        return false;
      }
    }

    const float dPhiAB = std::abs(deltaPhiMinusPiToPi((float)mcacc::lamPhi(A), (float)mcacc::lamPhi(B)));
    const float dPhiCB = std::abs(deltaPhiMinusPiToPi((float)mcacc::lamPhi(C), (float)mcacc::lamPhi(B)));
    if (std::abs(dPhiAB - dPhiCB) > phiMix) {
      return false;
    }

    return true;
  }

  void fillHistograms(int tag1, int tag2,
                      const ROOT::Math::PtEtaPhiMVector& particle1, const ROOT::Math::PtEtaPhiMVector& particle2,
                      const ROOT::Math::PtEtaPhiMVector& daughpart1, const ROOT::Math::PtEtaPhiMVector& daughpart2,
                      int datatype, float mixpairweight, int replacedLeg = 1)
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
    ROOT::Math::Boost boostPairToCM{pairDummy.BoostToCM()};

    // Step1: Boost both Lambdas to pair rest frame
    auto lambda1CM = boostPairToCM(particle1Dummy);
    auto lambda2CM = boostPairToCM(particle2Dummy);

    // Step2: Boost each Lambda to its own rest frame
    ROOT::Math::Boost boostLambda1ToCM{lambda1CM.BoostToCM()};
    ROOT::Math::Boost boostLambda2ToCM{lambda2CM.BoostToCM()};

    // Also boost daughter protons to pair CM
    auto proton1pairCM = boostPairToCM(daughpart1);
    auto proton2pairCM = boostPairToCM(daughpart2);

    // Then into each Lambda rest frame
    auto proton1LambdaRF = boostLambda1ToCM(proton1pairCM);
    auto proton2LambdaRF = boostLambda2ToCM(proton2pairCM);

    // STAR-style alternative
    ROOT::Math::Boost boostL1_LabToRF{particle1Dummy.BoostToCM()};
    ROOT::Math::Boost boostL2_LabToRF{particle2Dummy.BoostToCM()};

    auto p1_LRF = boostL1_LabToRF(daughpart1);
    auto p2_LRF = boostL2_LabToRF(daughpart2);

    TVector3 u1 = TVector3(p1_LRF.Px(), p1_LRF.Py(), p1_LRF.Pz()).Unit();
    TVector3 u2 = TVector3(p2_LRF.Px(), p2_LRF.Py(), p2_LRF.Pz()).Unit();

    TVector3 k1(proton1LambdaRF.Px(), proton1LambdaRF.Py(), proton1LambdaRF.Pz());
    k1 = k1.Unit();
    TVector3 k2(proton2LambdaRF.Px(), proton2LambdaRF.Py(), proton2LambdaRF.Pz());
    k2 = k2.Unit();

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

    double cosThetaDiff = (cosDef == 0) ? cosDeltaTheta_STAR_naive : cosDeltaTheta_hel;

    double pt1 = particle1.Pt();
    double dphi1 = RecoDecay::constrainAngle(particle1.Phi(), 0.0F, harmonic);
    double deta1 = particle1.Eta();

    double pt2 = particle2.Pt();
    double dphi2 = RecoDecay::constrainAngle(particle2.Phi(), 0.0F, harmonic);
    double deta2 = particle2.Eta();

    double dphi_pair = RecoDecay::constrainAngle(dphi1 - dphi2, -TMath::Pi(), harmonicDphi);
    double deltaRap = std::abs(particle1.Rapidity() - particle2.Rapidity());
    double deltaR = TMath::Sqrt(deltaRap * deltaRap + dphi_pair * dphi_pair);

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
      const double weight = 1.0;

      if (tag1 == 0 && tag2 == 0) {
        if (!userapidity) {
          histos.fill(HIST("hPtYSame"), particle1.Pt(), particle1.Rapidity(), weight);
          histos.fill(HIST("SE_LL"), dphi1, deta1, pt1, weight);
          histos.fill(HIST("SE_LL2"), dphi2, deta2, pt2, weight);
        } else {
          histos.fill(HIST("hPtYSame"), particle1.Pt(), particle1.Rapidity(), weight);
          histos.fill(HIST("SE_LL"), dphi1, particle1.Rapidity(), pt1, weight);
          histos.fill(HIST("SE_LL2"), dphi2, particle2.Rapidity(), pt2, weight);
        }
        histos.fill(HIST("hSparseLambdaLambda"), particle1.M(), particle2.M(), cosThetaDiff, deltaR, weight);
        histos.fill(HIST("hSparseLambdaLambdaAnalysis"), particle1.M(), particle2.M(), cosThetaDiff, deltaR, deltaRap, std::abs(dphi_pair), weight);
        if (useAdditionalHisto) {
          histos.fill(HIST("hSparseRapLambdaLambda"), particle1.M(), particle2.M(), cosThetaDiff, deltaRap, weight);
          histos.fill(HIST("hSparsePhiLambdaLambda"), particle1.M(), particle2.M(), cosThetaDiff, dphi_pair, weight);
          histos.fill(HIST("hSparsePairMassLambdaLambda"), particle1.M(), particle2.M(), cosThetaDiff, pairDummy.M(), weight);
        }
      } else if (tag1 == 0 && tag2 == 1) {
        if (!userapidity) {
          histos.fill(HIST("SE_LAL"), dphi1, deta1, pt1, weight);
          histos.fill(HIST("SE_LAL2"), dphi2, deta2, pt2, weight);
        } else {
          histos.fill(HIST("SE_LAL"), dphi1, particle1.Rapidity(), pt1, weight);
          histos.fill(HIST("SE_LAL2"), dphi2, particle2.Rapidity(), pt2, weight);
        }
        histos.fill(HIST("hSparseLambdaAntiLambda"), particle1.M(), particle2.M(), cosThetaDiff, deltaR, weight);
        histos.fill(HIST("hSparseLambdaAntiLambdaAnalysis"), particle1.M(), particle2.M(), cosThetaDiff, deltaR, deltaRap, std::abs(dphi_pair), weight);
        if (useAdditionalHisto) {
          histos.fill(HIST("hSparseRapLambdaAntiLambda"), particle1.M(), particle2.M(), cosThetaDiff, deltaRap, weight);
          histos.fill(HIST("hSparsePhiLambdaAntiLambda"), particle1.M(), particle2.M(), cosThetaDiff, dphi_pair, weight);
          histos.fill(HIST("hSparsePairMassLambdaAntiLambda"), particle1.M(), particle2.M(), cosThetaDiff, pairDummy.M(), weight);
        }
      } else if (tag1 == 1 && tag2 == 0) {
        if (!userapidity) {
          histos.fill(HIST("SE_ALL"), dphi1, deta1, pt1, weight);
          histos.fill(HIST("SE_ALL2"), dphi2, deta2, pt2, weight);
        } else {
          histos.fill(HIST("SE_ALL"), dphi1, particle1.Rapidity(), pt1, weight);
          histos.fill(HIST("SE_ALL2"), dphi2, particle2.Rapidity(), pt2, weight);
        }
        histos.fill(HIST("hSparseAntiLambdaLambda"), particle1.M(), particle2.M(), cosThetaDiff, deltaR, weight);
        histos.fill(HIST("hSparseAntiLambdaLambdaAnalysis"), particle1.M(), particle2.M(), cosThetaDiff, deltaR, deltaRap, std::abs(dphi_pair), weight);
        if (useAdditionalHisto) {
          histos.fill(HIST("hSparseRapAntiLambdaLambda"), particle1.M(), particle2.M(), cosThetaDiff, deltaRap, weight);
          histos.fill(HIST("hSparsePhiAntiLambdaLambda"), particle1.M(), particle2.M(), cosThetaDiff, dphi_pair, weight);
          histos.fill(HIST("hSparsePairMassAntiLambdaLambda"), particle1.M(), particle2.M(), cosThetaDiff, pairDummy.M(), weight);
        }
      } else if (tag1 == 1 && tag2 == 1) {
        if (!userapidity) {
          histos.fill(HIST("SE_ALAL"), dphi1, deta1, pt1, weight);
          histos.fill(HIST("SE_ALAL2"), dphi2, deta2, pt2, weight);
        } else {
          histos.fill(HIST("SE_ALAL"), dphi1, particle1.Rapidity(), pt1, weight);
          histos.fill(HIST("SE_ALAL2"), dphi2, particle2.Rapidity(), pt2, weight);
        }
        histos.fill(HIST("hSparseAntiLambdaAntiLambda"), particle1.M(), particle2.M(), cosThetaDiff, deltaR, weight);
        histos.fill(HIST("hSparseAntiLambdaAntiLambdaAnalysis"), particle1.M(), particle2.M(), cosThetaDiff, deltaR, deltaRap, std::abs(dphi_pair), weight);
        if (useAdditionalHisto) {
          histos.fill(HIST("hSparseRapAntiLambdaAntiLambda"), particle1.M(), particle2.M(), cosThetaDiff, deltaRap, weight);
          histos.fill(HIST("hSparsePhiAntiLambdaAntiLambda"), particle1.M(), particle2.M(), cosThetaDiff, dphi_pair, weight);
          histos.fill(HIST("hSparsePairMassAntiLambdaAntiLambda"), particle1.M(), particle2.M(), cosThetaDiff, pairDummy.M(), weight);
        }
      }

    } else if (datatype == 1) {
      double weight = mixpairweight;

      if (useweight) {
        const double epsWeightReplaced = (replacedLeg == 2) ? epsWeight2 : epsWeight1;
        if (!std::isfinite(epsWeightReplaced) || epsWeightReplaced <= 0.0) {
          return;
        }
        weight = mixpairweight / epsWeightReplaced;
      }

      if (!std::isfinite(weight) || weight <= 0.0) {
        return;
      }

      if (tag1 == 0 && tag2 == 0) {
        if (replacedLeg == 1) {
          if (!userapidity) {
            histos.fill(HIST("hPtYMix"), particle1.Pt(), particle1.Rapidity(), weight);
            histos.fill(HIST("ME_LL"), dphi1, deta1, pt1, weight);
          } else {
            histos.fill(HIST("hPtYMix"), particle1.Pt(), particle1.Rapidity(), weight);
            histos.fill(HIST("ME_LL"), dphi1, particle1.Rapidity(), pt1, weight);
          }
        } else {
          if (!userapidity) {
            histos.fill(HIST("ME_LL2"), dphi2, deta2, pt2, weight);
          } else {
            histos.fill(HIST("ME_LL2"), dphi2, particle2.Rapidity(), pt2, weight);
          }
        }
        histos.fill(HIST("hSparseLambdaLambdaMixed"), particle1.M(), particle2.M(), cosThetaDiff, deltaR, weight);
        histos.fill(HIST("hSparseLambdaLambdaMixedAnalysis"), particle1.M(), particle2.M(), cosThetaDiff, deltaR, deltaRap, std::abs(dphi_pair), weight);
        if (useAdditionalHisto) {
          histos.fill(HIST("hSparseRapLambdaLambdaMixed"), particle1.M(), particle2.M(), cosThetaDiff, deltaRap, weight);
          histos.fill(HIST("hSparsePhiLambdaLambdaMixed"), particle1.M(), particle2.M(), cosThetaDiff, dphi_pair, weight);
          histos.fill(HIST("hSparsePairMassLambdaLambdaMixed"), particle1.M(), particle2.M(), cosThetaDiff, pairDummy.M(), weight);
        }

      } else if (tag1 == 0 && tag2 == 1) {
        if (replacedLeg == 1) {
          if (!userapidity) {
            histos.fill(HIST("ME_LAL"), dphi1, deta1, pt1, weight);
          } else {
            histos.fill(HIST("ME_LAL"), dphi1, particle1.Rapidity(), pt1, weight);
          }
        } else {
          if (!userapidity) {
            histos.fill(HIST("ME_LAL2"), dphi2, deta2, pt2, weight);
          } else {
            histos.fill(HIST("ME_LAL2"), dphi2, particle2.Rapidity(), pt2, weight);
          }
        }
        histos.fill(HIST("hSparseLambdaAntiLambdaMixed"), particle1.M(), particle2.M(), cosThetaDiff, deltaR, weight);
        histos.fill(HIST("hSparseLambdaAntiLambdaMixedAnalysis"), particle1.M(), particle2.M(), cosThetaDiff, deltaR, deltaRap, std::abs(dphi_pair), weight);
        if (useAdditionalHisto) {
          histos.fill(HIST("hSparseRapLambdaAntiLambdaMixed"), particle1.M(), particle2.M(), cosThetaDiff, deltaRap, weight);
          histos.fill(HIST("hSparsePhiLambdaAntiLambdaMixed"), particle1.M(), particle2.M(), cosThetaDiff, dphi_pair, weight);
          histos.fill(HIST("hSparsePairMassLambdaAntiLambdaMixed"), particle1.M(), particle2.M(), cosThetaDiff, pairDummy.M(), weight);
        }

      } else if (tag1 == 1 && tag2 == 0) {
        if (replacedLeg == 1) {
          if (!userapidity) {
            histos.fill(HIST("ME_ALL"), dphi1, deta1, pt1, weight);
          } else {
            histos.fill(HIST("ME_ALL"), dphi1, particle1.Rapidity(), pt1, weight);
          }
        } else {
          if (!userapidity) {
            histos.fill(HIST("ME_ALL2"), dphi2, deta2, pt2, weight);
          } else {
            histos.fill(HIST("ME_ALL2"), dphi2, particle2.Rapidity(), pt2, weight);
          }
        }
        histos.fill(HIST("hSparseAntiLambdaLambdaMixed"), particle1.M(), particle2.M(), cosThetaDiff, deltaR, weight);
        histos.fill(HIST("hSparseAntiLambdaLambdaMixedAnalysis"), particle1.M(), particle2.M(), cosThetaDiff, deltaR, deltaRap, std::abs(dphi_pair), weight);
        if (useAdditionalHisto) {
          histos.fill(HIST("hSparseRapAntiLambdaLambdaMixed"), particle1.M(), particle2.M(), cosThetaDiff, deltaRap, weight);
          histos.fill(HIST("hSparsePhiAntiLambdaLambdaMixed"), particle1.M(), particle2.M(), cosThetaDiff, dphi_pair, weight);
          histos.fill(HIST("hSparsePairMassAntiLambdaLambdaMixed"), particle1.M(), particle2.M(), cosThetaDiff, pairDummy.M(), weight);
        }

      } else if (tag1 == 1 && tag2 == 1) {
        if (replacedLeg == 1) {
          if (!userapidity) {
            histos.fill(HIST("ME_ALAL"), dphi1, deta1, pt1, weight);
          } else {
            histos.fill(HIST("ME_ALAL"), dphi1, particle1.Rapidity(), pt1, weight);
          }
        } else {
          if (!userapidity) {
            histos.fill(HIST("ME_ALAL2"), dphi2, deta2, pt2, weight);
          } else {
            histos.fill(HIST("ME_ALAL2"), dphi2, particle2.Rapidity(), pt2, weight);
          }
        }
        histos.fill(HIST("hSparseAntiLambdaAntiLambdaMixed"), particle1.M(), particle2.M(), cosThetaDiff, deltaR, weight);
        histos.fill(HIST("hSparseAntiLambdaAntiLambdaMixedAnalysis"), particle1.M(), particle2.M(), cosThetaDiff, deltaR, deltaRap, std::abs(dphi_pair), weight);
        if (useAdditionalHisto) {
          histos.fill(HIST("hSparseRapAntiLambdaAntiLambdaMixed"), particle1.M(), particle2.M(), cosThetaDiff, deltaRap, weight);
          histos.fill(HIST("hSparsePhiAntiLambdaAntiLambdaMixed"), particle1.M(), particle2.M(), cosThetaDiff, dphi_pair, weight);
          histos.fill(HIST("hSparsePairMassAntiLambdaAntiLambdaMixed"), particle1.M(), particle2.M(), cosThetaDiff, pairDummy.M(), weight);
        }
      }
    }
  }
  static inline int pairTypeCode(int tag1, int tag2)
  {
    if (tag1 == 0 && tag2 == 0) {
      return 0; // LL
    } else if (tag1 == 0 && tag2 == 1) {
      return 1; // LAL
    } else if (tag1 == 1 && tag2 == 0) {
      return 2; // ALL
    } else {
      return 3; // ALAL
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
      histos.fill(HIST("hPtRadiusV0"), v0.lambdaPt(), v0.v0Radius());
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
        histos.fill(HIST("deltaPhiSame"), RecoDecay::constrainAngle(v0.lambdaPhi() - v02.lambdaPhi(), -TMath::Pi(), harmonicDphi));
        const int ptype = pairTypeCode(v0.v0Status(), v02.v0Status());
        histos.fill(HIST("hCentPairTypeSE"), collision.cent(), ptype, 1.0);
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

  void processMEV3(EventCandidates const& collisions, AllTrackCandidates const& V0s)
  {
    auto nBins = colBinning.getAllBinsCount();
    std::vector<std::deque<std::pair<int, AllTrackCandidates>>> eventPools(nBins);

    for (auto& collision1 : collisions) {
      const int bin = colBinning.getBin(std::make_tuple(collision1.posz(), collision1.cent()));
      if (bin < 0) {
        continue;
      }

      auto poolA = V0s.sliceBy(tracksPerCollisionV0, collision1.index());

      // if pool empty, push and continue
      if (eventPools[bin].empty()) {
        eventPools[bin].emplace_back(collision1.index(), std::move(poolA));
        if ((int)eventPools[bin].size() > nEvtMixing) {
          eventPools[bin].pop_front();
        }
        continue;
      }

      for (auto& [t1, t2] : soa::combinations(o2::soa::CombinationsFullIndexPolicy(poolA, poolA))) {
        if (!selectionV0(t1) || !selectionV0(t2)) {
          continue;
        }
        if (t2.index() <= t1.index()) {
          continue;
        }

        if (t1.protonIndex() == t2.protonIndex()) {
          continue;
        }
        if (t1.pionIndex() == t2.pionIndex()) {
          continue;
        }
        if (t1.protonIndex() == t2.pionIndex()) {
          continue;
        }
        if (t1.pionIndex() == t2.protonIndex()) {
          continue;
        }

        const bool doMixLeg1 = (cfgMixLegMode.value == 0 || cfgMixLegMode.value == 2);
        const bool doMixLeg2 = (cfgMixLegMode.value == 1 || cfgMixLegMode.value == 2);

        struct PV {
          AllTrackCandidates* pool;
          int nRepl1 = 0;
          int nRepl2 = 0;
        };

        std::vector<PV> usable;
        int totalRepl = 0;

        int mixes = 0;
        for (auto it = eventPools[bin].rbegin(); it != eventPools[bin].rend() && mixes < nEvtMixing; ++it, ++mixes) {
          const int collision2idx = it->first;
          auto& poolB = it->second;

          if (collision2idx == collision1.index()) {
            continue;
          }

          int nRepl1 = 0;
          int nRepl2 = 0;

          for (auto& tX : poolB) {
            if (!selectionV0(tX)) {
              continue;
            }

            if (doMixLeg1) {
              if (checkKinematics(t1, tX) && checkPairKinematics(t1, t2, tX)) {
                ++nRepl1;
              }
            }

            if (doMixLeg2) {
              if (checkKinematics(t2, tX) && checkPairKinematics(t2, t1, tX)) {
                ++nRepl2;
              }
            }
          }

          if (nRepl1 > 0 || nRepl2 > 0) {
            usable.push_back(PV{&poolB, nRepl1, nRepl2});
            totalRepl += nRepl1 + nRepl2;
          }
        }

        if (totalRepl <= 0) {
          continue;
        }

        const float wBase = 1.0f / static_cast<float>(totalRepl);

        for (auto& pv : usable) {
          auto& poolB = *pv.pool;

          for (auto& tX : poolB) {
            if (!selectionV0(tX)) {
              continue;
            }

            // -------- leg-1 replacement: (tX, t2)
            if (doMixLeg1) {
              if (checkKinematics(t1, tX) && checkPairKinematics(t1, t2, tX)) {
                auto proton = ROOT::Math::PtEtaPhiMVector(tX.protonPt(), tX.protonEta(), tX.protonPhi(),
                                                          o2::constants::physics::MassProton);
                auto lambda = ROOT::Math::PtEtaPhiMVector(tX.lambdaPt(), tX.lambdaEta(), tX.lambdaPhi(),
                                                          tX.lambdaMass());
                auto proton2 = ROOT::Math::PtEtaPhiMVector(t2.protonPt(), t2.protonEta(), t2.protonPhi(),
                                                           o2::constants::physics::MassProton);
                auto lambda2 = ROOT::Math::PtEtaPhiMVector(t2.lambdaPt(), t2.lambdaEta(), t2.lambdaPhi(),
                                                           t2.lambdaMass());

                const float dPhi = RecoDecay::constrainAngle(
                  RecoDecay::constrainAngle(lambda.Phi(), 0.0F, harmonic) -
                    RecoDecay::constrainAngle(lambda2.Phi(), 0.0F, harmonic),
                  -TMath::Pi(), harmonicDphi);

                histos.fill(HIST("deltaPhiMix"), dPhi, wBase);
                fillHistograms(tX.v0Status(), t2.v0Status(),
                               lambda, lambda2, proton, proton2,
                               1, wBase, 1);
              }
            }

            // -------- leg-2 replacement: (t1, tX)
            if (doMixLeg2) {
              if (checkKinematics(t2, tX) && checkPairKinematics(t2, t1, tX)) {
                auto proton = ROOT::Math::PtEtaPhiMVector(t1.protonPt(), t1.protonEta(), t1.protonPhi(),
                                                          o2::constants::physics::MassProton);
                auto lambda = ROOT::Math::PtEtaPhiMVector(t1.lambdaPt(), t1.lambdaEta(), t1.lambdaPhi(),
                                                          t1.lambdaMass());
                auto proton2 = ROOT::Math::PtEtaPhiMVector(tX.protonPt(), tX.protonEta(), tX.protonPhi(),
                                                           o2::constants::physics::MassProton);
                auto lambda2 = ROOT::Math::PtEtaPhiMVector(tX.lambdaPt(), tX.lambdaEta(), tX.lambdaPhi(),
                                                           tX.lambdaMass());

                const float dPhi = RecoDecay::constrainAngle(
                  RecoDecay::constrainAngle(lambda.Phi(), 0.0F, harmonic) -
                    RecoDecay::constrainAngle(lambda2.Phi(), 0.0F, harmonic),
                  -TMath::Pi(), harmonicDphi);

                histos.fill(HIST("deltaPhiMix"), dPhi, wBase);
                fillHistograms(t1.v0Status(), tX.v0Status(),
                               lambda, lambda2, proton, proton2,
                               1, wBase, 2);
              }
            }
          }
        }
      }

      // push current event into pool
      auto sliced = V0s.sliceBy(tracksPerCollisionV0, collision1.index());
      eventPools[bin].emplace_back(collision1.index(), std::move(sliced));
      if ((int)eventPools[bin].size() > nEvtMixing) {
        eventPools[bin].pop_front();
      }
    }
  }
  PROCESS_SWITCH(lambdaspincorrderived, processMEV3, "Process data ME (first-leg, pair-3D maps)", false);

  static constexpr int N_STATUS = 2; // v0Status ∈ {0,1}
  struct MatchRef {
    int64_t collisionIdx;
    int64_t rowIndex;
  };

  static inline void limitMatchesToNEvents(std::vector<MatchRef>& matches, int nMixEvents)
  {
    if (nMixEvents <= 0 || matches.empty()) {
      return;
    }

    std::vector<MatchRef> kept;
    kept.reserve(matches.size());

    std::unordered_set<int64_t> usedEvents;
    usedEvents.reserve(nMixEvents * 2);

    for (const auto& m : matches) {
      if (usedEvents.count(m.collisionIdx) || (int)usedEvents.size() < nMixEvents) {
        kept.push_back(m);
        usedEvents.insert(m.collisionIdx);
      }
    }

    matches.swap(kept);
  }

  struct MixBinnerR {
    float ptMin, ptMax, ptStep;
    float etaMin, etaMax, etaStep;
    float phiMin, phiMax, phiStep;

    float mMin, mMax, mStep;
    int nM_;

    std::vector<double> rEdges;
    int nR_;

    int nPt_, nEta_, nPhi_;

    MixBinnerR(float ptMin_, float ptMax_, float ptStep_,
               float etaAbsMax, float etaStep_,
               float phiStep_,
               float mMin_, float mMax_, int nMassBins_,
               std::vector<double> rEdges_)
      : ptMin(ptMin_),
        ptMax(ptMax_),
        ptStep(ptStep_),
        etaMin(-etaAbsMax),
        etaMax(+etaAbsMax),
        etaStep(etaStep_),
        phiMin(-static_cast<float>(TMath::Pi())),
        phiMax(+static_cast<float>(TMath::Pi())),
        phiStep(phiStep_),
        mMin(mMin_),
        mMax(mMax_),
        mStep(0.f),
        nM_(std::max(1, nMassBins_)),
        rEdges(std::move(rEdges_)),
        nR_(0),
        nPt_(0),
        nEta_(0),
        nPhi_(0)
    {
      ptStep = (ptStep > 0.f ? ptStep : 0.1f);
      etaStep = (etaStep > 0.f ? etaStep : 0.1f);
      phiStep = (phiStep > 0.f ? phiStep : 0.1f);

      if (!(mMax > mMin)) {
        mMin = 1.09f;
        mMax = 1.14f;
      }
      mStep = (mMax - mMin) / static_cast<float>(nM_);
      if (!(mStep > 0.f)) {
        nM_ = 5;
        mMin = 1.09f;
        mMax = 1.14f;
        mStep = (mMax - mMin) / static_cast<float>(nM_);
      }

      if (rEdges.size() < 2) {
        rEdges = {3.0, 5.0, 7.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0};
      }
      nR_ = static_cast<int>(rEdges.size()) - 1;

      nPt_ = std::max(1, static_cast<int>(std::floor((ptMax - ptMin) / ptStep + 0.5f)));
      nEta_ = std::max(1, static_cast<int>(std::floor((etaMax - etaMin) / etaStep + 0.5f)));
      nPhi_ = std::max(1, static_cast<int>(std::ceil((phiMax - phiMin) / phiStep)));
    }

    inline int nPt() const { return nPt_; }
    inline int nEta() const { return nEta_; }
    inline int nPhi() const { return nPhi_; }
    inline int nM() const { return nM_; }
    inline int nR() const { return nR_; }

    inline int binFromValue(float v, float vmin, float step, int nBins) const
    {
      if (!std::isfinite(v) || !std::isfinite(vmin) || !std::isfinite(step) || step <= 0.f || nBins <= 0) {
        return -1;
      }
      const float x = (v - vmin) / step;
      int b = static_cast<int>(std::floor(x + 1e-6f));
      if (b < 0) {
        return -1;
      }
      if (b >= nBins) {
        b = nBins - 1;
      }
      return b;
    }

    inline int ptBin(float pt) const { return binFromValue(pt, ptMin, ptStep, nPt_); }
    inline int etaBin(float eta) const { return binFromValue(eta, etaMin, etaStep, nEta_); }
    inline int phiBin(float phi) const { return binFromValue(phi, phiMin, phiStep, nPhi_); }
    inline int massBin(float m) const { return binFromValue(m, mMin, mStep, nM_); }

    inline int radiusBin(float r) const
    {
      if (!std::isfinite(r) || nR_ <= 0) {
        return -1;
      }
      if (r < rEdges.front() || r >= rEdges.back()) {
        return -1;
      }
      auto it = std::upper_bound(rEdges.begin(), rEdges.end(), static_cast<double>(r));
      return static_cast<int>(it - rEdges.begin()) - 1;
    }
  };

  struct BufferCandR {
    int64_t collisionIdx;
    int64_t rowIndex;
    uint8_t v0Status;
    uint16_t ptBin, etaBin, phiBin, mBin, rBin;
  };

  static inline size_t linearKeyR(int colBin, int statBin,
                                  int ptBin, int etaBin, int phiBin, int mBin, int rBin,
                                  int nStatus, int nPt, int nEta, int nPhi, int nM, int nR)
  {
    return (((((((static_cast<size_t>(colBin) * nStatus + statBin) * nPt + ptBin) * nEta + etaBin) * nPhi + phiBin) * nM + mBin) * nR + rBin));
  }

  // -------------------------------------
  // 2) MC-only selection + kinematics cuts
  // -------------------------------------
  template <typename T>
  bool selectionV0MC(T const& candidate)
  {
    auto particle = ROOT::Math::PtEtaPhiMVector(mcacc::lamPt(candidate),
                                                mcacc::lamEta(candidate),
                                                mcacc::lamPhi(candidate),
                                                mcacc::lamMass(candidate));
    if (std::abs(particle.Rapidity()) > rapidity || std::abs(particle.Eta()) > v0eta) {
      return false;
    }
    if (mcacc::lamMass(candidate) < MassMin || mcacc::lamMass(candidate) > MassMax) {
      return false;
    }
    if (mcacc::v0CosPA(candidate) < cosPA) {
      return false;
    }
    if (checkDoubleStatus && mcacc::doubleStatus(candidate)) {
      return false;
    }
    if (mcacc::v0Radius(candidate) > radiusMax) {
      return false;
    }
    if (mcacc::v0Radius(candidate) < radiusMin) {
      return false;
    }
    if (mcacc::dcaDau(candidate) > dcaDaughters) {
      return false;
    }
    if (mcacc::v0Status(candidate) == 0 && (std::abs(mcacc::dcaPos(candidate)) < dcaProton || std::abs(mcacc::dcaNeg(candidate)) < dcaPion)) {
      return false;
    }
    if (mcacc::v0Status(candidate) == 1 && (std::abs(mcacc::dcaPos(candidate)) < dcaPion || std::abs(mcacc::dcaNeg(candidate)) < dcaProton)) {
      return false;
    }
    if (mcacc::lamPt(candidate) < ptMin) {
      return false;
    }
    if (mcacc::lamPt(candidate) > ptMax) {
      return false;
    }
    return true;
  }

  template <typename T1, typename T2>
  bool checkKinematicsMC(T1 const& candidate1, T2 const& candidate2)
  {
    // keep same species/status
    if (mcacc::v0Status(candidate1) != mcacc::v0Status(candidate2)) {
      return false;
    }

    // pT window
    if (std::abs(mcacc::lamPt(candidate1) - mcacc::lamPt(candidate2)) > ptMix) {
      return false;
    }

    // eta or rapidity window (etaMix used as Δη or Δy)
    if (!userapidity) {
      if (std::abs(mcacc::lamEta(candidate1) - mcacc::lamEta(candidate2)) > etaMix) {
        return false;
      }
    } else {
      const auto l1 = ROOT::Math::PtEtaPhiMVector(mcacc::lamPt(candidate1), mcacc::lamEta(candidate1),
                                                  mcacc::lamPhi(candidate1), mcacc::lamMass(candidate1));
      const auto l2 = ROOT::Math::PtEtaPhiMVector(mcacc::lamPt(candidate2), mcacc::lamEta(candidate2),
                                                  mcacc::lamPhi(candidate2), mcacc::lamMass(candidate2));
      if (std::abs(l1.Rapidity() - l2.Rapidity()) > etaMix) {
        return false;
      }
    }

    // delta-phi window (wrapped)
    const float dphi = deltaPhiMinusPiToPi((float)mcacc::lamPhi(candidate1),
                                           (float)mcacc::lamPhi(candidate2));
    if (std::abs(dphi) > phiMix) {
      return false;
    }

    // mass window (optional but consistent with data)
    if (std::abs(mcacc::lamMass(candidate1) - mcacc::lamMass(candidate2)) > massMix) {
      return false;
    }

    return true;
  }

  // -----------------------------------------
  // 3) MC filter + aliases (distinct from data)
  // -----------------------------------------
  Filter centralityFilterMC = (nabs(aod::lambdaeventmc::centmc) < centMax && nabs(aod::lambdaeventmc::centmc) > centMin);

  using EventCandidatesMC = soa::Filtered<aod::LambdaEventmcs>;
  using AllTrackCandidatesMC = aod::LambdaPairmcs;

  // IMPORTANT: MC preslice uses the MC event index column
  Preslice<aod::LambdaPairmcs> tracksPerCollisionV0mc = aod::lambdapairmc::lambdaeventmcId;

  // -----------------------------------------
  // 4) MC Same-event processing (like processData)
  // -----------------------------------------
  void processMC(EventCandidatesMC::iterator const& collision, AllTrackCandidatesMC const& V0sMC)
  {
    const float centrality = mcacc::cent(collision);

    for (const auto& v0 : V0sMC) {
      if (!selectionV0MC(v0)) {
        continue;
      }
      histos.fill(HIST("hPtRadiusV0"), mcacc::lamPt(v0), mcacc::v0Radius(v0));
      histos.fill(HIST("ptCent"), mcacc::lamPt(v0), centrality);
      histos.fill(HIST("etaCent"), mcacc::lamEta(v0), centrality);

      proton = ROOT::Math::PtEtaPhiMVector(mcacc::prPt(v0), mcacc::prEta(v0), mcacc::prPhi(v0),
                                           o2::constants::physics::MassProton);
      lambda = ROOT::Math::PtEtaPhiMVector(mcacc::lamPt(v0), mcacc::lamEta(v0), mcacc::lamPhi(v0),
                                           mcacc::lamMass(v0));

      for (const auto& v02 : V0sMC) {
        if (v02.index() <= v0.index()) {
          continue;
        }
        if (!selectionV0MC(v02)) {
          continue;
        }

        // no shared daughters
        if (mcacc::prIdx(v0) == mcacc::prIdx(v02)) {
          continue;
        }
        if (mcacc::piIdx(v0) == mcacc::piIdx(v02)) {
          continue;
        }
        if (mcacc::prIdx(v0) == mcacc::piIdx(v02)) {
          continue;
        }
        if (mcacc::piIdx(v0) == mcacc::prIdx(v02)) {
          continue;
        }

        proton2 = ROOT::Math::PtEtaPhiMVector(mcacc::prPt(v02), mcacc::prEta(v02), mcacc::prPhi(v02),
                                              o2::constants::physics::MassProton);
        lambda2 = ROOT::Math::PtEtaPhiMVector(mcacc::lamPt(v02), mcacc::lamEta(v02), mcacc::lamPhi(v02),
                                              mcacc::lamMass(v02));

        histos.fill(HIST("deltaPhiSame"),
                    RecoDecay::constrainAngle(mcacc::lamPhi(v0) - mcacc::lamPhi(v02),
                                              -TMath::Pi(), harmonicDphi));

        const int ptype = pairTypeCode(mcacc::v0Status(v0), mcacc::v0Status(v02));
        histos.fill(HIST("hCentPairTypeSE"), mcacc::cent(collision), ptype, 1.0);
        // datatype=0 (same event)
        fillHistograms(mcacc::v0Status(v0), mcacc::v0Status(v02),
                       lambda, lambda2, proton, proton2,
                       /*datatype=*/0, /*mixpairweight=*/1.0f);
      }
    }
  }
  PROCESS_SWITCH(lambdaspincorrderived, processMC, "Process MC (SE)", false);

  void processMCMEV3(EventCandidatesMC const& collisions, AllTrackCandidatesMC const& V0sMC)
  {
    auto nBins = colBinning.getAllBinsCount();
    std::vector<std::deque<std::pair<int, AllTrackCandidatesMC>>> eventPools(nBins);

    for (auto& collision1 : collisions) {
      const int bin = colBinning.getBin(std::make_tuple(mcacc::posz(collision1), mcacc::cent(collision1)));
      if (bin < 0) {
        continue;
      }

      auto poolA = V0sMC.sliceBy(tracksPerCollisionV0mc, collision1.index());

      if (eventPools[bin].empty()) {
        eventPools[bin].emplace_back(collision1.index(), std::move(poolA));
        if ((int)eventPools[bin].size() > nEvtMixing) {
          eventPools[bin].pop_front();
        }
        continue;
      }

      for (auto& [t1, t2] : soa::combinations(o2::soa::CombinationsFullIndexPolicy(poolA, poolA))) {
        if (!selectionV0MC(t1) || !selectionV0MC(t2)) {
          continue;
        }
        if (t2.index() <= t1.index()) {
          continue;
        }

        if (mcacc::prIdx(t1) == mcacc::prIdx(t2)) {
          continue;
        }
        if (mcacc::piIdx(t1) == mcacc::piIdx(t2)) {
          continue;
        }
        if (mcacc::prIdx(t1) == mcacc::piIdx(t2)) {
          continue;
        }
        if (mcacc::piIdx(t1) == mcacc::prIdx(t2)) {
          continue;
        }

        const bool doMixLeg1 = (cfgMixLegMode.value == 0 || cfgMixLegMode.value == 2);
        const bool doMixLeg2 = (cfgMixLegMode.value == 1 || cfgMixLegMode.value == 2);

        struct PV {
          AllTrackCandidatesMC* pool;
          int nRepl1 = 0;
          int nRepl2 = 0;
        };

        std::vector<PV> usable;
        int totalRepl = 0;

        int mixes = 0;
        for (auto it = eventPools[bin].rbegin(); it != eventPools[bin].rend() && mixes < nEvtMixing; ++it, ++mixes) {
          const int collision2idx = it->first;
          auto& poolB = it->second;

          if (collision2idx == collision1.index()) {
            continue;
          }

          int nRepl1 = 0;
          int nRepl2 = 0;

          for (auto& tX : poolB) {
            if (!selectionV0MC(tX)) {
              continue;
            }

            if (doMixLeg1) {
              if (checkKinematicsMC(t1, tX) && checkPairKinematicsMC(t1, t2, tX)) {
                ++nRepl1;
              }
            }

            if (doMixLeg2) {
              if (checkKinematicsMC(t2, tX) && checkPairKinematicsMC(t2, t1, tX)) {
                ++nRepl2;
              }
            }
          }

          if (nRepl1 > 0 || nRepl2 > 0) {
            usable.push_back(PV{&poolB, nRepl1, nRepl2});
            totalRepl += nRepl1 + nRepl2;
          }
        }

        if (totalRepl <= 0) {
          continue;
        }

        const float wBase = 1.0f / static_cast<float>(totalRepl);

        for (auto& pv : usable) {
          auto& poolB = *pv.pool;

          for (auto& tX : poolB) {
            if (!selectionV0MC(tX)) {
              continue;
            }

            // -------- leg-1 replacement: (tX, t2)
            if (doMixLeg1) {
              if (checkKinematicsMC(t1, tX) && checkPairKinematicsMC(t1, t2, tX)) {
                auto pX = ROOT::Math::PtEtaPhiMVector(mcacc::prPt(tX), mcacc::prEta(tX), mcacc::prPhi(tX),
                                                      o2::constants::physics::MassProton);
                auto lX = ROOT::Math::PtEtaPhiMVector(mcacc::lamPt(tX), mcacc::lamEta(tX), mcacc::lamPhi(tX),
                                                      mcacc::lamMass(tX));
                auto p2 = ROOT::Math::PtEtaPhiMVector(mcacc::prPt(t2), mcacc::prEta(t2), mcacc::prPhi(t2),
                                                      o2::constants::physics::MassProton);
                auto l2 = ROOT::Math::PtEtaPhiMVector(mcacc::lamPt(t2), mcacc::lamEta(t2), mcacc::lamPhi(t2),
                                                      mcacc::lamMass(t2));

                const float dPhi = RecoDecay::constrainAngle(
                  RecoDecay::constrainAngle(lX.Phi(), 0.0F, harmonic) -
                    RecoDecay::constrainAngle(l2.Phi(), 0.0F, harmonic),
                  -TMath::Pi(), harmonicDphi);

                histos.fill(HIST("deltaPhiMix"), dPhi, wBase);
                fillHistograms(mcacc::v0Status(tX), mcacc::v0Status(t2),
                               lX, l2, pX, p2,
                               1, wBase, 1);
              }
            }

            // -------- leg-2 replacement: (t1, tX)
            if (doMixLeg2) {
              if (checkKinematicsMC(t2, tX) && checkPairKinematicsMC(t2, t1, tX)) {
                auto p1 = ROOT::Math::PtEtaPhiMVector(mcacc::prPt(t1), mcacc::prEta(t1), mcacc::prPhi(t1),
                                                      o2::constants::physics::MassProton);
                auto l1 = ROOT::Math::PtEtaPhiMVector(mcacc::lamPt(t1), mcacc::lamEta(t1), mcacc::lamPhi(t1),
                                                      mcacc::lamMass(t1));
                auto pX = ROOT::Math::PtEtaPhiMVector(mcacc::prPt(tX), mcacc::prEta(tX), mcacc::prPhi(tX),
                                                      o2::constants::physics::MassProton);
                auto lX = ROOT::Math::PtEtaPhiMVector(mcacc::lamPt(tX), mcacc::lamEta(tX), mcacc::lamPhi(tX),
                                                      mcacc::lamMass(tX));

                const float dPhi = RecoDecay::constrainAngle(
                  RecoDecay::constrainAngle(l1.Phi(), 0.0F, harmonic) -
                    RecoDecay::constrainAngle(lX.Phi(), 0.0F, harmonic),
                  -TMath::Pi(), harmonicDphi);

                histos.fill(HIST("deltaPhiMix"), dPhi, wBase);
                fillHistograms(mcacc::v0Status(t1), mcacc::v0Status(tX),
                               l1, lX, p1, pX,
                               1, wBase, 2);
              }
            }
          }
        }
      }

      auto sliced = V0sMC.sliceBy(tracksPerCollisionV0mc, collision1.index());
      eventPools[bin].emplace_back(collision1.index(), std::move(sliced));
      if ((int)eventPools[bin].size() > nEvtMixing) {
        eventPools[bin].pop_front();
      }
    }
  }
  PROCESS_SWITCH(lambdaspincorrderived, processMCMEV3, "Process MC ME v3 FIFO", false);
  static inline float phi0To2Pi(float phi)
  {
    // harmonic=1, min=0 => [0, 2pi)
    return RecoDecay::constrainAngle(phi, 0.0f, 1);
  }

  static inline float deltaPhiMinusPiToPi(float phiA, float phiB)
  {
    // returns in [-pi, pi)
    const float d = phi0To2Pi(phiA) - phi0To2Pi(phiB);
    return RecoDecay::constrainAngle(d, -TMath::Pi(), 1);
  }

  static inline float absDeltaPhi(float phiA, float phiB)
  {
    return std::abs(deltaPhiMinusPiToPi(phiA, phiB));
  }

  // symmetric neighbors for phi: periodic wrap
  static inline void collectNeighborBinsPhi(int b, int nPhi, int nNeighbor, std::vector<int>& out)
  {
    out.clear();
    out.reserve(2 * nNeighbor + 1);
    for (int d = -nNeighbor; d <= nNeighbor; ++d) {
      int bb = b + d;
      bb %= nPhi;
      if (bb < 0) {
        bb += nPhi;
      }
      out.push_back(bb);
    }
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
  }

  static inline void collectNeighborBinsClamp(int b, int nBins, int nNeighbor, std::vector<int>& out)
  {
    out.clear();
    out.reserve(2 * nNeighbor + 1);
    for (int d = -nNeighbor; d <= nNeighbor; ++d) {
      const int bb = b + d;
      if (bb >= 0 && bb < nBins) {
        out.push_back(bb);
      }
    }
  }

  static inline uint64_t splitmix64(uint64_t x)
  {
    // simple deterministic hash for reproducible shuffling
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
  }

  void processMEV6(EventCandidates const& collisions, AllTrackCandidates const& V0s)
  {
    MixBinnerR mb{
      ptMinMixBuffer.value, ptMaxMixBuffer.value, (ptMaxMixBuffer.value - ptMinMixBuffer.value) / cfgKinematicBins.nKinematicPt.value,
      v0etaMixBuffer.value, (2.0 * v0etaMixBuffer.value) / cfgKinematicBins.nKinematicEta.value,
      2.0 * TMath::Pi() / cfgKinematicBins.nKinematicPhi.value,
      MassMin.value, MassMax.value, cfgV5MassBins.value,
      cfgMixRadiusParam.cfgMixRadiusBins.value};

    const int nCol = colBinning.getAllBinsCount();
    const int nStat = N_STATUS;
    const int nPt = mb.nPt();
    const int nEta = mb.nEta();
    const int nPhi = mb.nPhi();
    const int nM = mb.nM();
    const int nR = mb.nR();

    const size_t nKeys = static_cast<size_t>(nCol) * nStat * nPt * nEta * nPhi * nM * nR;
    std::vector<std::vector<BufferCandR>> buffer(nKeys);

    // -------- PASS 1: fill buffer --------
    for (auto const& col : collisions) {
      const int colBin = colBinning.getBin(std::make_tuple(col.posz(), col.cent()));
      if (colBin < 0) {
        continue;
      }

      auto slice = V0s.sliceBy(tracksPerCollisionV0, col.index());

      for (auto const& t : slice) {
        if (!selectionV0Buffer(t)) {
          continue;
        }

        const int status = static_cast<int>(t.v0Status());
        if (status < 0 || status >= nStat) {
          continue;
        }

        const int ptB = mb.ptBin(t.lambdaPt());

        int etaB = mb.etaBin(t.lambdaEta());
        if (userapidity) {
          const auto lv = ROOT::Math::PtEtaPhiMVector(t.lambdaPt(), t.lambdaEta(), t.lambdaPhi(), t.lambdaMass());
          etaB = mb.etaBin(lv.Rapidity());
        }

        const int phiB = mb.phiBin(RecoDecay::constrainAngle(t.lambdaPhi(), -TMath::Pi(), harmonic));
        const int mB = mb.massBin(t.lambdaMass());
        const int rB = mb.radiusBin(t.v0Radius());

        if (ptB < 0 || etaB < 0 || phiB < 0 || mB < 0 || rB < 0) {
          continue;
        }

        const size_t key = linearKeyR(colBin, status, ptB, etaB, phiB, mB, rB,
                                      nStat, nPt, nEta, nPhi, nM, nR);

        buffer[key].push_back(BufferCandR{
          .collisionIdx = static_cast<int64_t>(col.index()),
          .rowIndex = static_cast<int64_t>(t.globalIndex()),
          .v0Status = static_cast<uint8_t>(status),
          .ptBin = static_cast<uint16_t>(ptB),
          .etaBin = static_cast<uint16_t>(etaB),
          .phiBin = static_cast<uint16_t>(phiB),
          .mBin = static_cast<uint16_t>(mB),
          .rBin = static_cast<uint16_t>(rB)});
      }
    }

    const int nN_pt = std::max(0, cfgV5NeighborPt.value);
    const int nN_eta = std::max(0, cfgV5NeighborEta.value);
    const int nN_phi = std::max(0, cfgV5NeighborPhi.value);

    std::vector<int> ptBins, etaBins, phiBins;
    std::vector<MatchRef> matches1, matches2;
    matches1.reserve(256);
    matches2.reserve(256);

    auto collectMatchesForReplacedLeg = [&](auto const& tRep, auto const& tKeep, int colBin, int64_t curColIdx, std::vector<MatchRef>& matches) {
      matches.clear();

      const int status = static_cast<int>(tRep.v0Status());
      if (status < 0 || status >= nStat) {
        return;
      }

      const int ptB = mb.ptBin(tRep.lambdaPt());

      int etaB = mb.etaBin(tRep.lambdaEta());
      if (userapidity) {
        const auto lv = ROOT::Math::PtEtaPhiMVector(tRep.lambdaPt(), tRep.lambdaEta(), tRep.lambdaPhi(), tRep.lambdaMass());
        etaB = mb.etaBin(lv.Rapidity());
      }

      const int phiB = mb.phiBin(RecoDecay::constrainAngle(tRep.lambdaPhi(), -TMath::Pi(), harmonic));
      const int mB = mb.massBin(tRep.lambdaMass());
      const int rB = mb.radiusBin(tRep.v0Radius());

      if (ptB < 0 || etaB < 0 || phiB < 0 || mB < 0 || rB < 0) {
        return;
      }

      collectNeighborBinsClamp(ptB, nPt, nN_pt, ptBins);
      collectNeighborBinsClamp(etaB, nEta, nN_eta, etaBins);
      collectNeighborBinsPhi(phiB, nPhi, nN_phi, phiBins);

      for (int ptUse : ptBins) {
        for (int etaUse : etaBins) {
          for (int phiUse : phiBins) {
            const auto& vec = buffer[linearKeyR(colBin, status, ptUse, etaUse, phiUse, mB, rB,
                                                nStat, nPt, nEta, nPhi, nM, nR)];

            for (auto const& bc : vec) {
              if (bc.collisionIdx == curColIdx) {
                continue;
              }

              auto tX = V0s.iteratorAt(static_cast<uint64_t>(bc.rowIndex));

              if (!selectionV0(tX)) {
                continue;
              }
              if (!checkKinematics(tRep, tX)) {
                continue;
              }
              if (!checkPairKinematics(tRep, tKeep, tX)) {
                continue;
              }

              if (tX.globalIndex() == tRep.globalIndex()) {
                continue;
              }
              if (tX.globalIndex() == tKeep.globalIndex()) {
                continue;
              }

              matches.push_back(MatchRef{bc.collisionIdx, bc.rowIndex});
            }
          }
        }
      }

      std::sort(matches.begin(), matches.end(),
                [](auto const& a, auto const& b) {
                  return std::tie(a.collisionIdx, a.rowIndex) < std::tie(b.collisionIdx, b.rowIndex);
                });
      matches.erase(std::unique(matches.begin(), matches.end(),
                                [](auto const& a, auto const& b) {
                                  return a.collisionIdx == b.collisionIdx && a.rowIndex == b.rowIndex;
                                }),
                    matches.end());
    };

    auto downsampleMatches = [&](std::vector<MatchRef>& matches, uint64_t seedBase) {
      if (cfgV5MaxMatches.value > 0 && (int)matches.size() > cfgV5MaxMatches.value) {
        uint64_t seed = cfgMixSeed.value ^ splitmix64(seedBase);
        const int K = cfgV5MaxMatches.value;
        for (int i = 0; i < K; ++i) {
          seed = splitmix64(seed);
          const int j = i + (int)(seed % (uint64_t)(matches.size() - i));
          std::swap(matches[i], matches[j]);
        }
        matches.resize(K);
      }
    };

    // -------- PASS 2: configurable one-leg / two-leg mixing --------
    for (auto const& col1 : collisions) {
      const int colBin = colBinning.getBin(std::make_tuple(col1.posz(), col1.cent()));
      if (colBin < 0) {
        continue;
      }

      const int64_t curColIdx = static_cast<int64_t>(col1.index());
      auto poolA = V0s.sliceBy(tracksPerCollisionV0, col1.index());

      for (auto const& [t1, t2] : soa::combinations(o2::soa::CombinationsFullIndexPolicy(poolA, poolA))) {
        if (!selectionV0(t1) || !selectionV0(t2)) {
          continue;
        }
        if (t2.index() <= t1.index()) {
          continue;
        }

        if (t1.protonIndex() == t2.protonIndex()) {
          continue;
        }
        if (t1.pionIndex() == t2.pionIndex()) {
          continue;
        }
        if (t1.protonIndex() == t2.pionIndex()) {
          continue;
        }
        if (t1.pionIndex() == t2.protonIndex()) {
          continue;
        }

        const bool doMixLeg1 = (cfgMixLegMode.value == 0 || cfgMixLegMode.value == 2);
        const bool doMixLeg2 = (cfgMixLegMode.value == 1 || cfgMixLegMode.value == 2);

        if (doMixLeg1) {
          collectMatchesForReplacedLeg(t1, t2, colBin, curColIdx, matches1);
          limitMatchesToNEvents(matches1, nEvtMixing.value);
          downsampleMatches(matches1, (uint64_t)t1.globalIndex() ^ (splitmix64((uint64_t)t2.globalIndex()) + 0x111ULL) ^ splitmix64((uint64_t)curColIdx));
        } else {
          matches1.clear();
        }

        if (doMixLeg2) {
          collectMatchesForReplacedLeg(t2, t1, colBin, curColIdx, matches2);
          limitMatchesToNEvents(matches2, nEvtMixing.value);
          downsampleMatches(matches2, (uint64_t)t2.globalIndex() ^ (splitmix64((uint64_t)t1.globalIndex()) + 0x222ULL) ^ splitmix64((uint64_t)curColIdx));
        } else {
          matches2.clear();
        }

        const int nReuse = static_cast<int>(matches1.size() + matches2.size());
        if (nReuse <= 0) {
          continue;
        }

        const float wSE = 1.0f / static_cast<float>(nReuse);

        if (doMixLeg1) {
          for (auto const& m : matches1) {
            auto tX = V0s.iteratorAt(static_cast<uint64_t>(m.rowIndex));

            auto proton = ROOT::Math::PtEtaPhiMVector(tX.protonPt(), tX.protonEta(), tX.protonPhi(), o2::constants::physics::MassProton);
            auto lambda = ROOT::Math::PtEtaPhiMVector(tX.lambdaPt(), tX.lambdaEta(), tX.lambdaPhi(), tX.lambdaMass());
            auto proton2 = ROOT::Math::PtEtaPhiMVector(t2.protonPt(), t2.protonEta(), t2.protonPhi(), o2::constants::physics::MassProton);
            auto lambda2 = ROOT::Math::PtEtaPhiMVector(t2.lambdaPt(), t2.lambdaEta(), t2.lambdaPhi(), t2.lambdaMass());

            const int ptype = pairTypeCode(tX.v0Status(), t2.v0Status());
            double centPairWeight = 1.0;
            if (hweightCentPair) {
              const int bin = hweightCentPair->FindBin(col1.cent(), ptype);
              centPairWeight = hweightCentPair->GetBinContent(bin);
              if (centPairWeight <= 0.0) {
                centPairWeight = 1.0;
              }
            }

            const float meWeight = wSE * centPairWeight;
            const float dPhi = deltaPhiMinusPiToPi((float)lambda.Phi(), (float)lambda2.Phi());
            histos.fill(HIST("deltaPhiMix"), dPhi, wSE);
            histos.fill(HIST("hCentPairTypeME"), col1.cent(), ptype, wSE);

            fillHistograms(tX.v0Status(), t2.v0Status(), lambda, lambda2, proton, proton2, 1, meWeight, 1);
          }
        }

        if (doMixLeg2) {
          for (auto const& m : matches2) {
            auto tY = V0s.iteratorAt(static_cast<uint64_t>(m.rowIndex));

            auto proton = ROOT::Math::PtEtaPhiMVector(t1.protonPt(), t1.protonEta(), t1.protonPhi(), o2::constants::physics::MassProton);
            auto lambda = ROOT::Math::PtEtaPhiMVector(t1.lambdaPt(), t1.lambdaEta(), t1.lambdaPhi(), t1.lambdaMass());
            auto proton2 = ROOT::Math::PtEtaPhiMVector(tY.protonPt(), tY.protonEta(), tY.protonPhi(), o2::constants::physics::MassProton);
            auto lambda2 = ROOT::Math::PtEtaPhiMVector(tY.lambdaPt(), tY.lambdaEta(), tY.lambdaPhi(), tY.lambdaMass());

            const int ptype = pairTypeCode(t1.v0Status(), tY.v0Status());
            double centPairWeight = 1.0;
            if (hweightCentPair) {
              const int bin = hweightCentPair->FindBin(col1.cent(), ptype);
              centPairWeight = hweightCentPair->GetBinContent(bin);
              if (centPairWeight <= 0.0) {
                centPairWeight = 1.0;
              }
            }

            const float meWeight = wSE * centPairWeight;
            const float dPhi = deltaPhiMinusPiToPi((float)lambda.Phi(), (float)lambda2.Phi());
            histos.fill(HIST("deltaPhiMix"), dPhi, wSE);
            histos.fill(HIST("hCentPairTypeME"), col1.cent(), ptype, wSE);

            fillHistograms(t1.v0Status(), tY.v0Status(), lambda, lambda2, proton, proton2, 1, meWeight, 2);
          }
        }
      }
    }
  }
  PROCESS_SWITCH(lambdaspincorrderived, processMEV6, "Process data ME v6 with radius buffer", false);
  void processMCMEV6(EventCandidatesMC const& collisions, AllTrackCandidatesMC const& V0sMC)
  {
    MixBinnerR mb{
      ptMinMixBuffer.value, ptMaxMixBuffer.value, (ptMaxMixBuffer.value - ptMinMixBuffer.value) / cfgKinematicBins.nKinematicPt.value,
      v0etaMixBuffer.value, (2.0 * v0etaMixBuffer.value) / cfgKinematicBins.nKinematicEta.value,
      2.0 * TMath::Pi() / cfgKinematicBins.nKinematicPhi.value,
      MassMin.value, MassMax.value, cfgV5MassBins.value,
      cfgMixRadiusParam.cfgMixRadiusBins.value};

    const int nCol = colBinning.getAllBinsCount();
    const int nStat = N_STATUS;
    const int nPt = mb.nPt();
    const int nEta = mb.nEta();
    const int nPhi = mb.nPhi();
    const int nM = mb.nM();
    const int nR = mb.nR();

    const size_t nKeys = static_cast<size_t>(nCol) * nStat * nPt * nEta * nPhi * nM * nR;
    std::vector<std::vector<BufferCandR>> buffer(nKeys);

    // -------- PASS 1: fill buffer --------
    for (auto const& col : collisions) {
      const int colBin = colBinning.getBin(std::make_tuple(mcacc::posz(col), mcacc::cent(col)));
      if (colBin < 0) {
        continue;
      }

      auto slice = V0sMC.sliceBy(tracksPerCollisionV0mc, col.index());

      for (auto const& t : slice) {
        if (!selectionV0BufferMC(t)) {
          continue;
        }

        const int status = mcacc::v0Status(t);
        if (status < 0 || status >= nStat) {
          continue;
        }

        const int ptB = mb.ptBin(mcacc::lamPt(t));

        int etaB = mb.etaBin(mcacc::lamEta(t));
        if (userapidity) {
          const auto lv = ROOT::Math::PtEtaPhiMVector(mcacc::lamPt(t), mcacc::lamEta(t), mcacc::lamPhi(t), mcacc::lamMass(t));
          etaB = mb.etaBin(lv.Rapidity());
        }

        const int phiB = mb.phiBin(RecoDecay::constrainAngle(mcacc::lamPhi(t), -TMath::Pi(), harmonic));
        const int mB = mb.massBin(mcacc::lamMass(t));
        const int rB = mb.radiusBin(mcacc::v0Radius(t));

        if (ptB < 0 || etaB < 0 || phiB < 0 || mB < 0 || rB < 0) {
          continue;
        }

        const size_t key = linearKeyR(colBin, status, ptB, etaB, phiB, mB, rB,
                                      nStat, nPt, nEta, nPhi, nM, nR);

        buffer[key].push_back(BufferCandR{
          .collisionIdx = static_cast<int64_t>(col.index()),
          .rowIndex = static_cast<int64_t>(t.globalIndex()),
          .v0Status = static_cast<uint8_t>(status),
          .ptBin = static_cast<uint16_t>(ptB),
          .etaBin = static_cast<uint16_t>(etaB),
          .phiBin = static_cast<uint16_t>(phiB),
          .mBin = static_cast<uint16_t>(mB),
          .rBin = static_cast<uint16_t>(rB)});
      }
    }

    const int nN_pt = std::max(0, cfgV5NeighborPt.value);
    const int nN_eta = std::max(0, cfgV5NeighborEta.value);
    const int nN_phi = std::max(0, cfgV5NeighborPhi.value);

    std::vector<int> ptBins, etaBins, phiBins;
    std::vector<MatchRef> matches1, matches2;
    matches1.reserve(256);
    matches2.reserve(256);

    auto collectMatchesForReplacedLeg = [&](auto const& tRep, auto const& tKeep, int colBin, int64_t curColIdx, std::vector<MatchRef>& matches) {
      matches.clear();

      const int status = mcacc::v0Status(tRep);
      if (status < 0 || status >= nStat) {
        return;
      }

      const int ptB = mb.ptBin(mcacc::lamPt(tRep));

      int etaB = mb.etaBin(mcacc::lamEta(tRep));
      if (userapidity) {
        const auto lv = ROOT::Math::PtEtaPhiMVector(mcacc::lamPt(tRep), mcacc::lamEta(tRep), mcacc::lamPhi(tRep), mcacc::lamMass(tRep));
        etaB = mb.etaBin(lv.Rapidity());
      }

      const int phiB = mb.phiBin(RecoDecay::constrainAngle(mcacc::lamPhi(tRep), -TMath::Pi(), harmonic));
      const int mB = mb.massBin(mcacc::lamMass(tRep));
      const int rB = mb.radiusBin(mcacc::v0Radius(tRep));

      if (ptB < 0 || etaB < 0 || phiB < 0 || mB < 0 || rB < 0) {
        return;
      }

      collectNeighborBinsClamp(ptB, nPt, nN_pt, ptBins);
      collectNeighborBinsClamp(etaB, nEta, nN_eta, etaBins);
      collectNeighborBinsPhi(phiB, nPhi, nN_phi, phiBins);

      for (int ptUse : ptBins) {
        for (int etaUse : etaBins) {
          for (int phiUse : phiBins) {
            const auto& vec = buffer[linearKeyR(colBin, status, ptUse, etaUse, phiUse, mB, rB,
                                                nStat, nPt, nEta, nPhi, nM, nR)];

            for (auto const& bc : vec) {
              if (bc.collisionIdx == curColIdx) {
                continue;
              }

              auto tX = V0sMC.iteratorAt(static_cast<uint64_t>(bc.rowIndex));

              if (!selectionV0MC(tX)) {
                continue;
              }
              if (!checkKinematicsMC(tRep, tX)) {
                continue;
              }
              if (!checkPairKinematicsMC(tRep, tKeep, tX)) {
                continue;
              }

              if (tX.globalIndex() == tRep.globalIndex()) {
                continue;
              }
              if (tX.globalIndex() == tKeep.globalIndex()) {
                continue;
              }

              matches.push_back(MatchRef{bc.collisionIdx, bc.rowIndex});
            }
          }
        }
      }

      std::sort(matches.begin(), matches.end(),
                [](auto const& a, auto const& b) {
                  return std::tie(a.collisionIdx, a.rowIndex) < std::tie(b.collisionIdx, b.rowIndex);
                });
      matches.erase(std::unique(matches.begin(), matches.end(),
                                [](auto const& a, auto const& b) {
                                  return a.collisionIdx == b.collisionIdx && a.rowIndex == b.rowIndex;
                                }),
                    matches.end());
    };

    auto downsampleMatches = [&](std::vector<MatchRef>& matches, uint64_t seedBase) {
      if (cfgV5MaxMatches.value > 0 && (int)matches.size() > cfgV5MaxMatches.value) {
        uint64_t seed = cfgMixSeed.value ^ splitmix64(seedBase);
        const int K = cfgV5MaxMatches.value;
        for (int i = 0; i < K; ++i) {
          seed = splitmix64(seed);
          const int j = i + (int)(seed % (uint64_t)(matches.size() - i));
          std::swap(matches[i], matches[j]);
        }
        matches.resize(K);
      }
    };

    // -------- PASS 2: configurable one-leg / two-leg mixing --------
    for (auto const& col1 : collisions) {
      const int colBin = colBinning.getBin(std::make_tuple(mcacc::posz(col1), mcacc::cent(col1)));
      if (colBin < 0) {
        continue;
      }

      const int64_t curColIdx = static_cast<int64_t>(col1.index());
      auto poolA = V0sMC.sliceBy(tracksPerCollisionV0mc, col1.index());

      for (auto const& [t1, t2] : soa::combinations(o2::soa::CombinationsFullIndexPolicy(poolA, poolA))) {
        if (!selectionV0MC(t1) || !selectionV0MC(t2)) {
          continue;
        }
        if (t2.index() <= t1.index()) {
          continue;
        }

        if (mcacc::prIdx(t1) == mcacc::prIdx(t2)) {
          continue;
        }
        if (mcacc::piIdx(t1) == mcacc::piIdx(t2)) {
          continue;
        }
        if (mcacc::prIdx(t1) == mcacc::piIdx(t2)) {
          continue;
        }
        if (mcacc::piIdx(t1) == mcacc::prIdx(t2)) {
          continue;
        }

        const bool doMixLeg1 = (cfgMixLegMode.value == 0 || cfgMixLegMode.value == 2);
        const bool doMixLeg2 = (cfgMixLegMode.value == 1 || cfgMixLegMode.value == 2);

        if (doMixLeg1) {
          collectMatchesForReplacedLeg(t1, t2, colBin, curColIdx, matches1);
          limitMatchesToNEvents(matches1, nEvtMixing.value);
          downsampleMatches(matches1, (uint64_t)t1.globalIndex() ^ (splitmix64((uint64_t)t2.globalIndex()) + 0x111ULL) ^ splitmix64((uint64_t)curColIdx));
        } else {
          matches1.clear();
        }
        if (doMixLeg2) {
          collectMatchesForReplacedLeg(t2, t1, colBin, curColIdx, matches2);
          limitMatchesToNEvents(matches2, nEvtMixing.value);
          downsampleMatches(matches2, (uint64_t)t2.globalIndex() ^ (splitmix64((uint64_t)t1.globalIndex()) + 0x222ULL) ^ splitmix64((uint64_t)curColIdx));
        } else {
          matches2.clear();
        }
        const int nReuse = static_cast<int>(matches1.size() + matches2.size());
        if (nReuse <= 0) {
          continue;
        }

        const float wSE = 1.0f / static_cast<float>(nReuse);

        if (doMixLeg1) {
          for (auto const& m : matches1) {
            auto tX = V0sMC.iteratorAt(static_cast<uint64_t>(m.rowIndex));

            auto pX = ROOT::Math::PtEtaPhiMVector(mcacc::prPt(tX), mcacc::prEta(tX), mcacc::prPhi(tX),
                                                  o2::constants::physics::MassProton);
            auto lX = ROOT::Math::PtEtaPhiMVector(mcacc::lamPt(tX), mcacc::lamEta(tX), mcacc::lamPhi(tX),
                                                  mcacc::lamMass(tX));
            auto p2 = ROOT::Math::PtEtaPhiMVector(mcacc::prPt(t2), mcacc::prEta(t2), mcacc::prPhi(t2),
                                                  o2::constants::physics::MassProton);
            auto l2 = ROOT::Math::PtEtaPhiMVector(mcacc::lamPt(t2), mcacc::lamEta(t2), mcacc::lamPhi(t2),
                                                  mcacc::lamMass(t2));

            const int ptype = pairTypeCode(mcacc::v0Status(tX), mcacc::v0Status(t2));
            double centPairWeight = 1.0;
            if (hweightCentPair) {
              const int bin = hweightCentPair->FindBin(mcacc::cent(col1), ptype);
              centPairWeight = hweightCentPair->GetBinContent(bin);
              if (centPairWeight <= 0.0) {
                centPairWeight = 1.0;
              }
            }

            const float meWeight = wSE * centPairWeight;
            const float dPhi = deltaPhiMinusPiToPi((float)lX.Phi(), (float)l2.Phi());
            histos.fill(HIST("deltaPhiMix"), dPhi, wSE);
            histos.fill(HIST("hCentPairTypeME"), mcacc::cent(col1), ptype, wSE);

            fillHistograms(mcacc::v0Status(tX), mcacc::v0Status(t2),
                           lX, l2, pX, p2,
                           1, meWeight, 1);
          }
        }

        if (doMixLeg2) {
          for (auto const& m : matches2) {
            auto tY = V0sMC.iteratorAt(static_cast<uint64_t>(m.rowIndex));

            auto p1 = ROOT::Math::PtEtaPhiMVector(mcacc::prPt(t1), mcacc::prEta(t1), mcacc::prPhi(t1),
                                                  o2::constants::physics::MassProton);
            auto l1 = ROOT::Math::PtEtaPhiMVector(mcacc::lamPt(t1), mcacc::lamEta(t1), mcacc::lamPhi(t1),
                                                  mcacc::lamMass(t1));
            auto pY = ROOT::Math::PtEtaPhiMVector(mcacc::prPt(tY), mcacc::prEta(tY), mcacc::prPhi(tY),
                                                  o2::constants::physics::MassProton);
            auto lY = ROOT::Math::PtEtaPhiMVector(mcacc::lamPt(tY), mcacc::lamEta(tY), mcacc::lamPhi(tY),
                                                  mcacc::lamMass(tY));

            const int ptype = pairTypeCode(mcacc::v0Status(t1), mcacc::v0Status(tY));
            double centPairWeight = 1.0;
            if (hweightCentPair) {
              const int bin = hweightCentPair->FindBin(mcacc::cent(col1), ptype);
              centPairWeight = hweightCentPair->GetBinContent(bin);
              if (centPairWeight <= 0.0) {
                centPairWeight = 1.0;
              }
            }

            const float meWeight = wSE * centPairWeight;
            const float dPhi = deltaPhiMinusPiToPi((float)l1.Phi(), (float)lY.Phi());
            histos.fill(HIST("deltaPhiMix"), dPhi, wSE);
            histos.fill(HIST("hCentPairTypeME"), mcacc::cent(col1), ptype, wSE);

            fillHistograms(mcacc::v0Status(t1), mcacc::v0Status(tY),
                           l1, lY, p1, pY,
                           1, meWeight, 2);
          }
        }
      }
    }
  }
  PROCESS_SWITCH(lambdaspincorrderived, processMCMEV6, "Process MC ME v6 with radius buffer", false);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<lambdaspincorrderived>(cfgc)};
}
