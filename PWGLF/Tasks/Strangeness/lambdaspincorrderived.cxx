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

#include "Common/Core/RecoDecay.h"

#include <CCDB/BasicCCDBManager.h>
#include <CCDB/CcdbApi.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <Math/GenVector/Boost.h>
#include <Math/Vector4D.h> // IWYU pragma: keep (do not replace with Math/Vector4Dfwd.h)
#include <Math/Vector4Dfwd.h>
#include <TH3.h>
#include <TMath.h>
#include <TVector3.h>

#include <algorithm>
#include <chrono>
#include <cmath> // for std::abs
#include <cstddef>
#include <cstdint>
#include <deque>
#include <iterator>
#include <random>
#include <set> // <<< CHANGED: for dedup sets
#include <string>
#include <unordered_map> // <<< CHANGED: for seenMap
#include <utility>
#include <vector>

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

    // double deta_pair = std::abs(deta1 - deta2);
    double dphi_pair = RecoDecay::constrainAngle(dphi1 - dphi2, -TMath::Pi(), harmonicDphi);
    // double deltaR = TMath::Sqrt(deta_pair * deta_pair + dphi_pair * dphi_pair);
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
      mixpairweight = 1.0;
      histos.fill(HIST("hPtYSame"), particle1.Pt(), particle1.Rapidity(), mixpairweight);
      if (tag1 == 0 && tag2 == 0) {
        if (!userapidity) {
          histos.fill(HIST("SE_LL"), dphi1, deta1, pt1, mixpairweight);
          histos.fill(HIST("SE_LL2"), dphi2, deta2, pt2, mixpairweight);
        } else {
          histos.fill(HIST("SE_LL"), dphi1, particle1.Rapidity(), pt1, mixpairweight);
          histos.fill(HIST("SE_LL2"), dphi2, particle2.Rapidity(), pt2, mixpairweight);
        }
        histos.fill(HIST("hSparseLambdaLambda"), particle1.M(), particle2.M(), cosThetaDiff, deltaR, mixpairweight);
        histos.fill(HIST("hSparseLambdaLambdaAnalysis"), particle1.M(), particle2.M(), cosThetaDiff, deltaR, deltaRap, std::abs(dphi_pair), mixpairweight);
        if (useAdditionalHisto) {
          histos.fill(HIST("hSparseRapLambdaLambda"), particle1.M(), particle2.M(), cosThetaDiff, deltaRap, mixpairweight);
          histos.fill(HIST("hSparsePhiLambdaLambda"), particle1.M(), particle2.M(), cosThetaDiff, dphi_pair, mixpairweight);
          histos.fill(HIST("hSparsePairMassLambdaLambda"), particle1.M(), particle2.M(), cosThetaDiff, pairDummy.M(), mixpairweight);
        }
      } else if (tag1 == 0 && tag2 == 1) {
        if (!userapidity) {
          histos.fill(HIST("SE_LAL"), dphi1, deta1, pt1, mixpairweight);
          histos.fill(HIST("SE_LAL2"), dphi2, deta2, pt2, mixpairweight);
        } else {
          histos.fill(HIST("SE_LAL"), dphi1, particle1.Rapidity(), pt1, mixpairweight);
          histos.fill(HIST("SE_LAL2"), dphi2, particle2.Rapidity(), pt2, mixpairweight);
        }
        histos.fill(HIST("hSparseLambdaAntiLambda"), particle1.M(), particle2.M(), cosThetaDiff, deltaR, mixpairweight);
        histos.fill(HIST("hSparseLambdaAntiLambdaAnalysis"), particle1.M(), particle2.M(), cosThetaDiff, deltaR, deltaRap, std::abs(dphi_pair), mixpairweight);
        if (useAdditionalHisto) {
          histos.fill(HIST("hSparseRapLambdaAntiLambda"), particle1.M(), particle2.M(), cosThetaDiff, deltaRap, mixpairweight);
          histos.fill(HIST("hSparsePhiLambdaAntiLambda"), particle1.M(), particle2.M(), cosThetaDiff, dphi_pair, mixpairweight);
          histos.fill(HIST("hSparsePairMassLambdaAntiLambda"), particle1.M(), particle2.M(), cosThetaDiff, pairDummy.M(), mixpairweight);
        }
      } else if (tag1 == 1 && tag2 == 0) {
        histos.fill(HIST("hSparseAntiLambdaLambda"), particle1.M(), particle2.M(), cosThetaDiff, deltaR, mixpairweight);
        histos.fill(HIST("hSparseAntiLambdaLambdaAnalysis"), particle1.M(), particle2.M(), cosThetaDiff, deltaR, deltaRap, std::abs(dphi_pair), mixpairweight);
        if (!userapidity) {
          histos.fill(HIST("SE_ALL"), dphi1, deta1, pt1, mixpairweight);
          histos.fill(HIST("SE_ALL2"), dphi2, deta2, pt2, mixpairweight);
        } else {
          histos.fill(HIST("SE_ALL"), dphi1, particle1.Rapidity(), pt1, mixpairweight);
          histos.fill(HIST("SE_ALL2"), dphi2, particle2.Rapidity(), pt2, mixpairweight);
        }
        if (useAdditionalHisto) {
          histos.fill(HIST("hSparseRapAntiLambdaLambda"), particle1.M(), particle2.M(), cosThetaDiff, deltaRap, mixpairweight);
          histos.fill(HIST("hSparsePhiAntiLambdaLambda"), particle1.M(), particle2.M(), cosThetaDiff, dphi_pair, mixpairweight);
          histos.fill(HIST("hSparsePairMassAntiLambdaLambda"), particle1.M(), particle2.M(), cosThetaDiff, pairDummy.M(), mixpairweight);
        }
      } else if (tag1 == 1 && tag2 == 1) {
        histos.fill(HIST("hSparseAntiLambdaAntiLambda"), particle1.M(), particle2.M(), cosThetaDiff, deltaR, mixpairweight);
        histos.fill(HIST("hSparseAntiLambdaAntiLambdaAnalysis"), particle1.M(), particle2.M(), cosThetaDiff, deltaR, deltaRap, std::abs(dphi_pair), mixpairweight);
        if (!userapidity) {
          histos.fill(HIST("SE_ALAL"), dphi1, deta1, pt1, mixpairweight);
          histos.fill(HIST("SE_ALAL2"), dphi2, deta2, pt2, mixpairweight);
        } else {
          histos.fill(HIST("SE_ALAL"), dphi1, particle1.Rapidity(), pt1, mixpairweight);
          histos.fill(HIST("SE_ALAL2"), dphi2, particle2.Rapidity(), pt2, mixpairweight);
        }
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
        if (!userapidity) {
          histos.fill(HIST("ME_LL"), dphi1, deta1, pt1, mixpairweight);
          histos.fill(HIST("ME_LL2"), dphi2, deta2, pt2, mixpairweight);
        } else {
          histos.fill(HIST("ME_LL"), dphi1, particle1.Rapidity(), pt1, mixpairweight);
          histos.fill(HIST("ME_LL2"), dphi2, particle2.Rapidity(), pt2, mixpairweight);
        }
        histos.fill(HIST("hSparseLambdaLambdaMixed"), particle1.M(), particle2.M(), cosThetaDiff, deltaR, weight);
        histos.fill(HIST("hSparseLambdaLambdaMixedAnalysis"), particle1.M(), particle2.M(), cosThetaDiff, deltaR, deltaRap, std::abs(dphi_pair), weight);
        if (useAdditionalHisto) {
          histos.fill(HIST("hSparseRapLambdaLambdaMixed"), particle1.M(), particle2.M(), cosThetaDiff, deltaRap, weight);
          histos.fill(HIST("hSparsePhiLambdaLambdaMixed"), particle1.M(), particle2.M(), cosThetaDiff, dphi_pair, weight);
          histos.fill(HIST("hSparsePairMassLambdaLambdaMixed"), particle1.M(), particle2.M(), cosThetaDiff, pairDummy.M(), weight);
        }
      } else if (tag1 == 0 && tag2 == 1) {
        if (!userapidity) {
          histos.fill(HIST("ME_LAL"), dphi1, deta1, pt1, mixpairweight);
          histos.fill(HIST("ME_LAL2"), dphi2, deta2, pt2, mixpairweight);
        } else {
          histos.fill(HIST("ME_LAL"), dphi1, particle1.Rapidity(), pt1, mixpairweight);
          histos.fill(HIST("ME_LAL2"), dphi2, particle2.Rapidity(), pt2, mixpairweight);
        }
        histos.fill(HIST("hSparseLambdaAntiLambdaMixed"), particle1.M(), particle2.M(), cosThetaDiff, deltaR, weight);
        histos.fill(HIST("hSparseLambdaAntiLambdaMixedAnalysis"), particle1.M(), particle2.M(), cosThetaDiff, deltaR, deltaRap, std::abs(dphi_pair), weight);
        if (useAdditionalHisto) {
          histos.fill(HIST("hSparseRapLambdaAntiLambdaMixed"), particle1.M(), particle2.M(), cosThetaDiff, deltaRap, weight);
          histos.fill(HIST("hSparsePhiLambdaAntiLambdaMixed"), particle1.M(), particle2.M(), cosThetaDiff, dphi_pair, weight);
          histos.fill(HIST("hSparsePairMassLambdaAntiLambdaMixed"), particle1.M(), particle2.M(), cosThetaDiff, pairDummy.M(), weight);
        }
      } else if (tag1 == 1 && tag2 == 0) {
        if (!userapidity) {
          histos.fill(HIST("ME_ALL"), dphi1, deta1, pt1, mixpairweight);
          histos.fill(HIST("ME_ALL2"), dphi2, deta2, pt2, mixpairweight);
        } else {
          histos.fill(HIST("ME_ALL"), dphi1, particle1.Rapidity(), pt1, mixpairweight);
          histos.fill(HIST("ME_ALL2"), dphi2, particle2.Rapidity(), pt2, mixpairweight);
        }
        histos.fill(HIST("hSparseAntiLambdaLambdaMixed"), particle1.M(), particle2.M(), cosThetaDiff, deltaR, weight);
        histos.fill(HIST("hSparseAntiLambdaLambdaMixedAnalysis"), particle1.M(), particle2.M(), cosThetaDiff, deltaR, deltaRap, std::abs(dphi_pair), weight);
        if (useAdditionalHisto) {
          histos.fill(HIST("hSparseRapAntiLambdaLambdaMixed"), particle1.M(), particle2.M(), cosThetaDiff, deltaRap, weight);
          histos.fill(HIST("hSparsePhiAntiLambdaLambdaMixed"), particle1.M(), particle2.M(), cosThetaDiff, dphi_pair, weight);
          histos.fill(HIST("hSparsePairMassAntiLambdaLambdaMixed"), particle1.M(), particle2.M(), cosThetaDiff, pairDummy.M(), weight);
        }
      } else if (tag1 == 1 && tag2 == 1) {
        if (!userapidity) {
          histos.fill(HIST("ME_ALAL"), dphi1, deta1, pt1, mixpairweight);
          histos.fill(HIST("ME_ALAL2"), dphi2, deta2, pt2, mixpairweight);
        } else {
          histos.fill(HIST("ME_ALAL"), dphi1, particle1.Rapidity(), pt1, mixpairweight);
          histos.fill(HIST("ME_ALAL2"), dphi2, particle2.Rapidity(), pt2, mixpairweight);
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
          histos.fill(HIST("deltaPhiMix"), RecoDecay::constrainAngle(t3.lambdaPhi() - t2.lambdaPhi(), -TMath::Pi(), harmonicDphi));
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

            float dPhi = RecoDecay::constrainAngle(RecoDecay::constrainAngle(lambda.Phi(), 0.0F, harmonic) - RecoDecay::constrainAngle(lambda2.Phi(), 0.0, harmonic), -TMath::Pi(), harmonicDphi);
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

            const float dPhi = RecoDecay::constrainAngle(RecoDecay::constrainAngle(lambda.Phi(), 0.0F, harmonic) - RecoDecay::constrainAngle(lambda2.Phi(), 0.0F, harmonic), -TMath::Pi(), harmonicDphi);
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
    // constructed from the task's configurables; φ is assumed already constrained upstream
    float ptMin, ptMax, ptStep;
    float etaMin, etaMax, etaStep;
    float phiMin, phiMax, phiStep;

    // configurable fixed mass-binning for mixing buffer
    float mMin, mMax, mStep;
    int nM_;

    int nPt_, nEta_, nPhi_;

    MixBinner(float ptMin_, float ptMax_, float ptStep_,
              float etaAbsMax, float etaStep_,
              float phiStep_,
              float mMin_, float mMax_, int nMassBins_)
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
        nPt_(0),
        nEta_(0),
        nPhi_(0)
    // If you want phi in [0, 2pi), use:
    // : ... phiMin(0.f), phiMax(static_cast<float>(2.0 * TMath::Pi())), ...
    {
      ptStep = (ptStep > 0.f ? ptStep : 0.1f);
      etaStep = (etaStep > 0.f ? etaStep : 0.1f);
      phiStep = (phiStep > 0.f ? phiStep : 0.1f);

      if (!(mMax > mMin)) {
        mMin = 1.09f;
        mMax = 1.14f;
      }
      mStep = (mMax - mMin) / static_cast<float>(nM_);

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
      if (!std::isfinite(v) || !std::isfinite(vmin) || !std::isfinite(step) || step <= 0.f || nBins <= 0) {
        return -1;
      }
      const float x = (v - vmin) / step;
      int b = static_cast<int>(std::floor(x + 1e-6f));
      if (b < 0) {
        return -1;
      }
      if (b >= nBins) {
        b = nBins - 1; // clamp exact-top edge
      }
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

  static inline void collectPhiNeighborBins(int phiB, int nPhi, int nNeighbor, std::vector<int>& out)
  {
    out.clear();
    out.reserve(2 * nNeighbor + 1);
    for (int d = -nNeighbor; d <= nNeighbor; ++d) {
      int b = phiB + d;
      // wrap into [0, nPhi-1]
      b %= nPhi;
      if (b < 0)
        b += nPhi;
      out.push_back(b);
    }
    // optional: unique (in case nNeighbor >= nPhi)
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
  }

  // ===================== Main mixing (with mass-bin + random unique sampling) =====================
  void processMEV4(EventCandidates const& collisions, AllTrackCandidates const& V0s)
  {
    MixBinner mb{
      ptMin.value, ptMax.value, ptMix.value,
      v0etaMixBuffer.value, etaMix.value,
      phiMix.value,
      MassMin.value, MassMax.value, cfgV5MassBins.value};

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
        const int phiB = mb.phiBin(RecoDecay::constrainAngle(t.lambdaPhi(), -TMath::Pi(), harmonic));
        const int mB = mb.massBin(t.lambdaMass());
        if (ptB < 0 || etaB < 0 || phiB < 0 || mB < 0)
          continue;

        const size_t key = linearKey(colBin, status, ptB, etaB, phiB, mB,
                                     nStat, nPt, nEta, nPhi, nM);

        buffer[key].push_back(BufferCand{
          .collisionIdx = static_cast<int64_t>(col.index()),
          .rowIndex = static_cast<int64_t>(t.globalIndex()),
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

        // Bin of t1 defines where to search (exact bin, but handle φ wrap at edges)
        const int ptB = mb.ptBin(t1.lambdaPt());
        const int etaB = mb.etaBin(t1.lambdaEta());
        const int phiB = mb.phiBin(RecoDecay::constrainAngle(t1.lambdaPhi(), -TMath::Pi(), harmonic));
        const int mB = mb.massBin(t1.lambdaMass());
        if (ptB < 0 || etaB < 0 || phiB < 0 || mB < 0)
          continue;

        // Collect partners from nominal key, plus wrapped neighbor only for φ-edge bins
        std::vector<MatchRef> matches;
        matches.reserve(128); // or keep binVec.size() if you prefer
        const int64_t curColIdx = static_cast<int64_t>(collision1.index());

        auto collectFrom = [&](int phiBinUse) {
          const size_t keyUse = linearKey(colBin, status, ptB, etaB, phiBinUse, mB,
                                          nStat, nPt, nEta, nPhi, nM);
          auto const& vec = buffer[keyUse];
          for (const auto& bc : vec) {
            if (bc.collisionIdx == curColIdx) {
              continue; // must be from different event
            }
            auto tX = V0s.iteratorAt(static_cast<uint64_t>(bc.rowIndex));
            if (!selectionV0(tX)) {
              continue;
            }
            if (!checkKinematics(t1, tX)) {
              continue;
            }
            if (!checkPairKinematics(t1, t2, tX)) {
              continue;
            }
            matches.push_back(MatchRef{bc.collisionIdx, bc.rowIndex});
          }
        };
        // 1) nominal φ-bin
        collectFrom(phiB);

        // 2) wrap only at boundaries: 0 <-> nPhi-1
        if (phiB == 0) {
          collectFrom(nPhi - 1);
        } else if (phiB == nPhi - 1) {
          collectFrom(0);
        }

        if (matches.empty()) {
          continue;
        }

        // Optional safety: dedupe exact same (collision,row) just in case
        std::sort(matches.begin(), matches.end(),
                  [](auto& a, auto& b) { return std::tie(a.collisionIdx, a.rowIndex) < std::tie(b.collisionIdx, b.rowIndex); });
        matches.erase(std::unique(matches.begin(), matches.end(),
                                  [](auto& a, auto& b) { return a.collisionIdx == b.collisionIdx && a.rowIndex == b.rowIndex; }),
                      matches.end());
        if (matches.empty()) {
          continue;
        }
        const float wBase = 1.0f / static_cast<float>(matches.size());
        for (const auto& m : matches) {
          auto tX = V0s.iteratorAt(static_cast<uint64_t>(m.rowIndex));

          auto proton = ROOT::Math::PtEtaPhiMVector(tX.protonPt(), tX.protonEta(), tX.protonPhi(), o2::constants::physics::MassProton);
          auto lambda = ROOT::Math::PtEtaPhiMVector(tX.lambdaPt(), tX.lambdaEta(), tX.lambdaPhi(), tX.lambdaMass());
          auto proton2 = ROOT::Math::PtEtaPhiMVector(t2.protonPt(), t2.protonEta(), t2.protonPhi(), o2::constants::physics::MassProton);
          auto lambda2 = ROOT::Math::PtEtaPhiMVector(t2.lambdaPt(), t2.lambdaEta(), t2.lambdaPhi(), t2.lambdaMass());

          const float dPhi = RecoDecay::constrainAngle(
            RecoDecay::constrainAngle(lambda.Phi(), 0.0F, harmonic) - RecoDecay::constrainAngle(lambda2.Phi(), 0.0F, harmonic),
            -TMath::Pi(), harmonicDphi);

          histos.fill(HIST("deltaPhiMix"), dPhi, wBase);
          fillHistograms(tX.v0Status(), t2.v0Status(), lambda, lambda2, proton, proton2, 1, wBase);
        }
      }
    }
  }
  PROCESS_SWITCH(lambdaspincorrderived, processMEV4, "Process data ME (5d buffer)", false);

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
      const int bin = colBinning.getBin(std::make_tuple(collision1.poszmc(), collision1.centmc()));

      // if pool empty, push and continue
      if (eventPools[bin].empty()) {
        auto sliced = V0sMC.sliceBy(tracksPerCollisionV0mc, collision1.index());
        eventPools[bin].emplace_back(collision1.index(), std::move(sliced));
        if ((int)eventPools[bin].size() > nEvtMixing) {
          eventPools[bin].pop_front();
        }
        continue;
      }

      // current event slice
      auto poolA = V0sMC.sliceBy(tracksPerCollisionV0mc, collision1.index());

      // loop over SE unordered pairs (t1,t2)
      for (auto& [t1, t2] : soa::combinations(o2::soa::CombinationsFullIndexPolicy(poolA, poolA))) {

        // ---- selections ----
        if (!selectionV0MC(t1) || !selectionV0MC(t2)) {
          continue;
        }
        if (t2.index() <= t1.index()) {
          continue;
        }

        // no shared daughters (use global indices stored in your MC table)
        if (t1.protonIndexmc() == t2.protonIndexmc())
          continue;
        if (t1.pionIndexmc() == t2.pionIndexmc())
          continue;
        if (t1.protonIndexmc() == t2.pionIndexmc())
          continue;
        if (t1.pionIndexmc() == t2.protonIndexmc())
          continue;

        // scan prior events for replacements for t1
        struct PV {
          AllTrackCandidatesMC* pool;
          int nRepl;
        };
        std::vector<PV> usable;
        int totalRepl = 0;

        int mixes = 0;
        for (auto it = eventPools[bin].rbegin();
             it != eventPools[bin].rend() && mixes < nEvtMixing; ++it, ++mixes) {

          const int collision2idx = it->first;
          auto& poolB = it->second;
          if (collision2idx == collision1.index()) {
            continue;
          }

          int nRepl = 0;
          for (auto& tX : poolB) {
            if (!selectionV0MC(tX))
              continue;
            if (checkKinematicsMC(t1, tX))
              ++nRepl;
          }
          if (nRepl > 0) {
            usable.push_back(PV{&poolB, nRepl});
            totalRepl += nRepl;
          }
        }

        if (totalRepl == 0) {
          continue;
        }
        const float wBase = 1.0f / static_cast<float>(totalRepl);

        // emit mixed pairs: tX replaces t1; t2 stays
        for (auto& pv : usable) {
          auto& poolB = *pv.pool;
          for (auto& tX : poolB) {
            if (!selectionV0MC(tX))
              continue;
            if (!checkKinematicsMC(t1, tX))
              continue;

            // build 4-vectors
            auto proton = ROOT::Math::PtEtaPhiMVector(tX.protonPtmc(), tX.protonEtamc(), tX.protonPhimc(), o2::constants::physics::MassProton);
            auto lambda = ROOT::Math::PtEtaPhiMVector(tX.lambdaPtmc(), tX.lambdaEtamc(), tX.lambdaPhimc(), tX.lambdaMassmc());
            auto proton2 = ROOT::Math::PtEtaPhiMVector(t2.protonPtmc(), t2.protonEtamc(), t2.protonPhimc(), o2::constants::physics::MassProton);
            auto lambda2 = ROOT::Math::PtEtaPhiMVector(t2.lambdaPtmc(), t2.lambdaEtamc(), t2.lambdaPhimc(), t2.lambdaMassmc());

            const float dPhi = RecoDecay::constrainAngle(
              RecoDecay::constrainAngle(lambda.Phi(), 0.0F, harmonic) -
                RecoDecay::constrainAngle(lambda2.Phi(), 0.0F, harmonic),
              -TMath::Pi(), harmonicDphi);

            histos.fill(HIST("deltaPhiMix"), dPhi, wBase);
            fillHistograms(tX.v0Statusmc(), t2.v0Statusmc(), lambda, lambda2, proton, proton2, 1, wBase);
          }
        }
      } // end SE pair loop

      // push current event into pool
      auto sliced = V0sMC.sliceBy(tracksPerCollisionV0mc, collision1.index());
      eventPools[bin].emplace_back(collision1.index(), std::move(sliced));
      if ((int)eventPools[bin].size() > nEvtMixing) {
        eventPools[bin].pop_front();
      }
    } // end events
  }

  // enable it
  PROCESS_SWITCH(lambdaspincorrderived, processMCMEV3, "Process MC ME (MEV3)", false);

  // -----------------------------------------------------
  // 5) MC Event Mixing using your MEV4 6D-buffer approach
  // -----------------------------------------------------
  void processMCMEV4(EventCandidatesMC const& collisions, AllTrackCandidatesMC const& V0sMC)
  {
    MixBinner mb{
      ptMin.value, ptMax.value, ptMix.value,
      v0etaMixBuffer.value, etaMix.value,
      phiMix.value,
      MassMin.value, MassMax.value, cfgV5MassBins.value};

    const int nCol = colBinning.getAllBinsCount();
    const int nStat = N_STATUS;
    const int nPt = mb.nPt();
    const int nEta = mb.nEta();
    const int nPhi = mb.nPhi();
    const int nM = mb.nM();

    const size_t nKeys = static_cast<size_t>(nCol) * nStat * nPt * nEta * nPhi * nM;
    std::vector<std::vector<BufferCand>> buffer(nKeys);

    // ---- PASS 1: fill 6D buffer from MC tables ----
    for (auto const& col : collisions) {
      const int colBin = colBinning.getBin(std::make_tuple(mcacc::posz(col), mcacc::cent(col)));
      auto slice = V0sMC.sliceBy(tracksPerCollisionV0mc, col.index());

      for (auto const& t : slice) {
        if (!selectionV0MC(t)) {
          continue;
        }

        const int status = mcacc::v0Status(t);
        if (status < 0 || status >= nStat) {
          continue;
        }

        const int ptB = mb.ptBin(mcacc::lamPt(t));
        const int etaB = mb.etaBin(mcacc::lamEta(t));
        const int phiB = mb.phiBin(RecoDecay::constrainAngle(mcacc::lamPhi(t), -TMath::Pi(), harmonic));
        const int mB = mb.massBin(mcacc::lamMass(t));
        if (ptB < 0 || etaB < 0 || phiB < 0 || mB < 0) {
          continue;
        }

        const size_t key = linearKey(colBin, status, ptB, etaB, phiB, mB,
                                     nStat, nPt, nEta, nPhi, nM);

        buffer[key].push_back(BufferCand{
          .collisionIdx = static_cast<int64_t>(col.index()),
          .rowIndex = static_cast<int64_t>(t.globalIndex()),
          .v0Status = static_cast<uint8_t>(status),
          .ptBin = static_cast<uint16_t>(ptB),
          .etaBin = static_cast<uint16_t>(etaB),
          .phiBin = static_cast<uint16_t>(phiB),
          .mBin = static_cast<uint16_t>(mB)});
      }
    }

    // ---- PASS 2: build mixed pairs for each same-event pair (t1,t2) ----
    for (auto const& collision1 : collisions) {
      const int colBin = colBinning.getBin(std::make_tuple(mcacc::posz(collision1), mcacc::cent(collision1)));
      auto poolA = V0sMC.sliceBy(tracksPerCollisionV0mc, collision1.index());

      for (auto const& [t1, t2] :
           soa::combinations(o2::soa::CombinationsFullIndexPolicy(poolA, poolA))) {

        if (!selectionV0MC(t1) || !selectionV0MC(t2)) {
          continue;
        }
        if (t2.index() <= t1.index()) {
          continue;
        }

        // no shared daughters
        if (mcacc::prIdx(t1) == mcacc::prIdx(t2))
          continue;
        if (mcacc::piIdx(t1) == mcacc::piIdx(t2))
          continue;
        if (mcacc::prIdx(t1) == mcacc::piIdx(t2))
          continue;
        if (mcacc::piIdx(t1) == mcacc::prIdx(t2))
          continue;

        const int status = mcacc::v0Status(t1);
        if (status < 0 || status >= nStat) {
          continue;
        }

        const int ptB = mb.ptBin(mcacc::lamPt(t1));
        const int etaB = mb.etaBin(mcacc::lamEta(t1));
        const int phiB = mb.phiBin(RecoDecay::constrainAngle(mcacc::lamPhi(t1), -TMath::Pi(), harmonic));
        const int mB = mb.massBin(mcacc::lamMass(t1));
        if (ptB < 0 || etaB < 0 || phiB < 0 || mB < 0) {
          continue;
        }
        std::vector<MatchRef> matches;
        matches.reserve(128);
        const int64_t curColIdx = static_cast<int64_t>(collision1.index());
        auto collectFrom = [&](int phiBinUse) {
          const size_t keyUse = linearKey(colBin, status, ptB, etaB, phiBinUse, mB,
                                          nStat, nPt, nEta, nPhi, nM);
          auto const& vec = buffer[keyUse];
          for (const auto& bc : vec) {
            if (bc.collisionIdx == curColIdx) {
              continue; // different event
            }
            auto tX = V0sMC.iteratorAt(static_cast<uint64_t>(bc.rowIndex));
            if (!selectionV0MC(tX)) {
              continue;
            }
            if (!checkKinematicsMC(t1, tX)) {
              continue;
            }
            if (!checkPairKinematicsMC(t1, t2, tX)) {
              continue;
            }
            matches.push_back(MatchRef{bc.collisionIdx, bc.rowIndex});
          }
        };

        // nominal φ-bin + wrap neighbors only at edges
        collectFrom(phiB);
        if (phiB == 0) {
          collectFrom(nPhi - 1);
        } else if (phiB == nPhi - 1) {
          collectFrom(0);
        }

        if (matches.empty()) {
          continue;
        }

        // dedupe identical (collision,row)
        std::sort(matches.begin(), matches.end(),
                  [](auto& a, auto& b) { return std::tie(a.collisionIdx, a.rowIndex) < std::tie(b.collisionIdx, b.rowIndex); });
        matches.erase(std::unique(matches.begin(), matches.end(),
                                  [](auto& a, auto& b) { return a.collisionIdx == b.collisionIdx && a.rowIndex == b.rowIndex; }),
                      matches.end());
        if (matches.empty()) {
          continue;
        }

        const float wBase = 1.0f / static_cast<float>(matches.size());

        for (const auto& m : matches) {
          auto tX = V0sMC.iteratorAt(static_cast<uint64_t>(m.rowIndex));

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

          // datatype=1 (mixed event)
          fillHistograms(mcacc::v0Status(tX), mcacc::v0Status(t2),
                         lX, l2, pX, p2,
                         /*datatype=*/1, /*mixpairweight=*/wBase);
        }
      }
    }
  }
  PROCESS_SWITCH(lambdaspincorrderived, processMCMEV4, "Process MC ME (5d buffer)", false);

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

  // symmetric neighbors for continuous mixing (pt/eta): include bin, ±1, ±2..., edge-safe
  static inline void collectNeighborBins1D(int b, int nBins, int nNeighbor, std::vector<int>& out)
  {
    out.clear();
    out.reserve(2 * nNeighbor + 1);
    for (int d = -nNeighbor; d <= nNeighbor; ++d) {
      const int bb = b + d;
      if (bb < 0 || bb >= nBins) {
        continue;
      }
      out.push_back(bb);
    }
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
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

  static inline void collectPhiBinsWithEdgeWrap(int phiB, int nPhi, std::vector<int>& out)
  {
    out.clear();
    out.reserve(2);
    out.push_back(phiB);
    if (nPhi <= 1) {
      return;
    }
    if (phiB == 0) {
      out.push_back(nPhi - 1);
    } else if (phiB == nPhi - 1) {
      out.push_back(0);
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

  static inline uint64_t splitmixmc64(uint64_t x)
  {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
  }

  void processMEV5(EventCandidates const& collisions, AllTrackCandidates const& V0s)
  {
    MixBinner mb{
      ptMin.value, ptMax.value, ptMix.value,
      v0etaMixBuffer.value, etaMix.value,
      phiMix.value,
      MassMin.value, MassMax.value, cfgV5MassBins.value};

    const int nCol = colBinning.getAllBinsCount();
    const int nStat = N_STATUS;
    const int nPt = mb.nPt();
    const int nEta = mb.nEta(); // logical "nY" if userapidity=true
    const int nPhi = mb.nPhi();
    const int nM = mb.nM();

    const size_t nKeys = static_cast<size_t>(nCol) * nStat * nPt * nEta * nPhi * nM;
    std::vector<std::vector<BufferCand>> buffer(nKeys);

    // -------- PASS 1: fill buffer --------
    for (auto const& col : collisions) {
      const int colBin = colBinning.getBin(std::make_tuple(col.posz(), col.cent()));
      if (colBin < 0) {
        continue;
      }

      auto slice = V0s.sliceBy(tracksPerCollisionV0, col.index());

      for (auto const& t : slice) {
        if (!selectionV0(t)) {
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
          etaB = mb.etaBin(lv.Rapidity()); // treat "eta axis" as rapidity axis
        }

        const int phiB = mb.phiBin(RecoDecay::constrainAngle(t.lambdaPhi(), -TMath::Pi(), harmonic));
        const int mB = mb.massBin(t.lambdaMass());

        if (ptB < 0 || etaB < 0 || phiB < 0 || mB < 0) {
          continue;
        }

        const size_t key = linearKey(colBin, status, ptB, etaB, phiB, mB,
                                     nStat, nPt, nEta, nPhi, nM);

        buffer[key].push_back(BufferCand{
          .collisionIdx = static_cast<int64_t>(col.index()),
          .rowIndex = static_cast<int64_t>(t.globalIndex()),
          .v0Status = static_cast<uint8_t>(status),
          .ptBin = static_cast<uint16_t>(ptB),
          .etaBin = static_cast<uint16_t>(etaB),
          .phiBin = static_cast<uint16_t>(phiB),
          .mBin = static_cast<uint16_t>(mB)});
      }
    }

    // Neighbor policy from configurables
    const int nN_pt = std::max(0, cfgV5NeighborPt.value);
    const int nN_eta = std::max(0, cfgV5NeighborEta.value);
    const int nN_phi = std::max(0, cfgV5NeighborPhi.value);

    std::vector<int> ptBins, etaBins, phiBins;
    std::vector<MatchRef> matches;
    matches.reserve(256);

    // -------- PASS 2: mix (replace t1 by tX, keep t2 from same event) --------
    for (auto const& col1 : collisions) {
      const int colBin = colBinning.getBin(std::make_tuple(col1.posz(), col1.cent()));
      if (colBin < 0) {
        continue;
      }

      const int64_t curColIdx = static_cast<int64_t>(col1.index());
      auto poolA = V0s.sliceBy(tracksPerCollisionV0, col1.index());

      for (auto const& [t1, t2] :
           soa::combinations(o2::soa::CombinationsFullIndexPolicy(poolA, poolA))) {

        if (!selectionV0(t1) || !selectionV0(t2)) {
          continue;
        }
        if (t2.index() <= t1.index()) {
          continue; // same-event ordering
        }

        // no shared daughters (same-event)
        if (t1.protonIndex() == t2.protonIndex())
          continue;
        if (t1.pionIndex() == t2.pionIndex())
          continue;
        if (t1.protonIndex() == t2.pionIndex())
          continue;
        if (t1.pionIndex() == t2.protonIndex())
          continue;

        const int status = static_cast<int>(t1.v0Status());
        if (status < 0 || status >= nStat) {
          continue;
        }

        const int ptB = mb.ptBin(t1.lambdaPt());

        int etaB = mb.etaBin(t1.lambdaEta());
        if (userapidity) {
          const auto lv1 = ROOT::Math::PtEtaPhiMVector(t1.lambdaPt(), t1.lambdaEta(), t1.lambdaPhi(), t1.lambdaMass());
          etaB = mb.etaBin(lv1.Rapidity());
        }

        const int phiB = mb.phiBin(RecoDecay::constrainAngle(t1.lambdaPhi(), -TMath::Pi(), harmonic));
        const int mB = mb.massBin(t1.lambdaMass());

        if (ptB < 0 || etaB < 0 || phiB < 0 || mB < 0) {
          continue;
        }

        collectNeighborBinsClamp(ptB, nPt, nN_pt, ptBins);
        collectNeighborBinsClamp(etaB, nEta, nN_eta, etaBins);
        collectNeighborBinsPhi(phiB, nPhi, nN_phi, phiBins);

        matches.clear();

        for (int ptUse : ptBins) {
          for (int etaUse : etaBins) {
            for (int phiUse : phiBins) {
              const size_t keyUse = linearKey(colBin, status, ptUse, etaUse, phiUse, mB,
                                              nStat, nPt, nEta, nPhi, nM);
              auto const& vec = buffer[keyUse];

              for (auto const& bc : vec) {
                if (bc.collisionIdx == curColIdx) {
                  continue; // enforce different event
                }

                auto tX = V0s.iteratorAt(static_cast<uint64_t>(bc.rowIndex));
                if (!selectionV0(tX)) {
                  continue;
                }

                if (!checkKinematics(t1, tX)) {
                  continue;
                }

                if (tX.globalIndex() == t1.globalIndex())
                  continue;
                if (tX.globalIndex() == t2.globalIndex())
                  continue;

                matches.push_back(MatchRef{bc.collisionIdx, bc.rowIndex});
              }
            }
          }
        }

        if (matches.empty()) {
          continue;
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
        if (matches.empty()) {
          continue;
        }

        if (cfgV5MaxMatches.value > 0 && (int)matches.size() > cfgV5MaxMatches.value) {
          uint64_t seed = cfgMixSeed.value;
          seed ^= splitmix64((uint64_t)t1.globalIndex());
          seed ^= splitmix64((uint64_t)t2.globalIndex() + 0x1234567ULL);
          seed ^= splitmix64((uint64_t)curColIdx + 0x9abcULL);

          const int K = cfgV5MaxMatches.value;
          for (int i = 0; i < K; ++i) {
            seed = splitmix64(seed);
            const int j = i + (int)(seed % (uint64_t)(matches.size() - i));
            std::swap(matches[i], matches[j]);
          }
          matches.resize(K);
        }

        const float wBase = 1.0f / static_cast<float>(matches.size());

        for (auto const& m : matches) {
          auto tX = V0s.iteratorAt(static_cast<uint64_t>(m.rowIndex));

          auto proton = ROOT::Math::PtEtaPhiMVector(tX.protonPt(), tX.protonEta(), tX.protonPhi(),
                                                    o2::constants::physics::MassProton);
          auto lambda = ROOT::Math::PtEtaPhiMVector(tX.lambdaPt(), tX.lambdaEta(), tX.lambdaPhi(),
                                                    tX.lambdaMass());

          auto proton2 = ROOT::Math::PtEtaPhiMVector(t2.protonPt(), t2.protonEta(), t2.protonPhi(),
                                                     o2::constants::physics::MassProton);
          auto lambda2 = ROOT::Math::PtEtaPhiMVector(t2.lambdaPt(), t2.lambdaEta(), t2.lambdaPhi(),
                                                     t2.lambdaMass());

          const int ptype = pairTypeCode(tX.v0Status(), t2.v0Status());
          double centPairWeight = 1.0;
          if (hweightCentPair) {
            const int bin = hweightCentPair->FindBin(col1.cent(), ptype);
            centPairWeight = hweightCentPair->GetBinContent(bin);
            if (centPairWeight <= 0.0) {
              centPairWeight = 1.0;
            }
          }
          const float meWeight = wBase * centPairWeight;
          const float dPhi = deltaPhiMinusPiToPi((float)lambda.Phi(), (float)lambda2.Phi());
          histos.fill(HIST("deltaPhiMix"), dPhi, wBase);
          histos.fill(HIST("hCentPairTypeME"), col1.cent(), ptype, wBase);
          fillHistograms(tX.v0Status(), t2.v0Status(),
                         lambda, lambda2, proton, proton2,
                         /*datatype=*/1, /*mixpairweight=*/meWeight);
        }
      }
    }
  }
  PROCESS_SWITCH(lambdaspincorrderived, processMEV5, "Process data ME v5", false);

  void processMCMEV5(EventCandidatesMC const& collisions, AllTrackCandidatesMC const& V0sMC)
  {
    MixBinner mb{
      ptMin.value, ptMax.value, ptMix.value,
      v0etaMixBuffer.value, etaMix.value,
      phiMix.value,
      MassMin.value, MassMax.value, cfgV5MassBins.value};

    const int nCol = colBinning.getAllBinsCount();
    const int nStat = N_STATUS;
    const int nPt = mb.nPt();
    const int nEta = mb.nEta(); // logical "nY" if userapidity=true
    const int nPhi = mb.nPhi();
    const int nM = mb.nM();

    const size_t nKeys = static_cast<size_t>(nCol) * nStat * nPt * nEta * nPhi * nM;
    std::vector<std::vector<BufferCand>> buffer(nKeys);

    // -------- PASS 1: fill buffer --------
    for (auto const& col : collisions) {
      const int colBin = colBinning.getBin(std::make_tuple(mcacc::posz(col), mcacc::cent(col)));
      if (colBin < 0) {
        continue;
      }

      auto slice = V0sMC.sliceBy(tracksPerCollisionV0mc, col.index());

      for (auto const& t : slice) {
        if (!selectionV0MC(t)) {
          continue;
        }

        const int status = mcacc::v0Status(t);
        if (status < 0 || status >= nStat) {
          continue;
        }

        const int ptB = mb.ptBin(mcacc::lamPt(t));

        int etaB = mb.etaBin(mcacc::lamEta(t));
        if (userapidity) {
          const auto lv = ROOT::Math::PtEtaPhiMVector(mcacc::lamPt(t), mcacc::lamEta(t),
                                                      mcacc::lamPhi(t), mcacc::lamMass(t));
          etaB = mb.etaBin(lv.Rapidity());
        }

        const int phiB = mb.phiBin(RecoDecay::constrainAngle(mcacc::lamPhi(t), -TMath::Pi(), harmonic));
        const int mB = mb.massBin(mcacc::lamMass(t));
        if (ptB < 0 || etaB < 0 || phiB < 0 || mB < 0) {
          continue;
        }

        buffer[linearKey(colBin, status, ptB, etaB, phiB, mB,
                         nStat, nPt, nEta, nPhi, nM)]
          .push_back(BufferCand{
            .collisionIdx = static_cast<int64_t>(col.index()),
            .rowIndex = static_cast<int64_t>(t.globalIndex()),
            .v0Status = static_cast<uint8_t>(status),
            .ptBin = static_cast<uint16_t>(ptB),
            .etaBin = static_cast<uint16_t>(etaB),
            .phiBin = static_cast<uint16_t>(phiB),
            .mBin = static_cast<uint16_t>(mB)});
      }
    }

    const int nN_pt = std::max(0, cfgV5NeighborPt.value);
    const int nN_eta = std::max(0, cfgV5NeighborEta.value);
    const int nN_phi = std::max(0, cfgV5NeighborPhi.value);

    std::vector<int> ptBins, etaBins, phiBins;
    std::vector<MatchRef> matches;
    matches.reserve(256);

    // -------- PASS 2: build ME --------
    for (auto const& col1 : collisions) {
      const int colBin = colBinning.getBin(std::make_tuple(mcacc::posz(col1), mcacc::cent(col1)));
      if (colBin < 0) {
        continue;
      }

      const int64_t curColIdx = static_cast<int64_t>(col1.index());
      auto poolA = V0sMC.sliceBy(tracksPerCollisionV0mc, col1.index());

      for (auto const& [t1, t2] :
           soa::combinations(o2::soa::CombinationsFullIndexPolicy(poolA, poolA))) {

        if (!selectionV0MC(t1) || !selectionV0MC(t2)) {
          continue;
        }
        if (t2.index() <= t1.index()) {
          continue;
        }

        if (mcacc::prIdx(t1) == mcacc::prIdx(t2))
          continue;
        if (mcacc::piIdx(t1) == mcacc::piIdx(t2))
          continue;
        if (mcacc::prIdx(t1) == mcacc::piIdx(t2))
          continue;
        if (mcacc::piIdx(t1) == mcacc::prIdx(t2))
          continue;

        const int status = mcacc::v0Status(t1);
        if (status < 0 || status >= nStat) {
          continue;
        }

        const int ptB = mb.ptBin(mcacc::lamPt(t1));

        int etaB = mb.etaBin(mcacc::lamEta(t1));
        if (userapidity) {
          const auto lv1 = ROOT::Math::PtEtaPhiMVector(mcacc::lamPt(t1), mcacc::lamEta(t1),
                                                       mcacc::lamPhi(t1), mcacc::lamMass(t1));
          etaB = mb.etaBin(lv1.Rapidity());
        }
        const int phiB = mb.phiBin(RecoDecay::constrainAngle(mcacc::lamPhi(t1), -TMath::Pi(), harmonic));
        const int mB = mb.massBin(mcacc::lamMass(t1));
        if (ptB < 0 || etaB < 0 || phiB < 0 || mB < 0) {
          continue;
        }
        collectNeighborBinsClamp(ptB, nPt, nN_pt, ptBins);
        collectNeighborBinsClamp(etaB, nEta, nN_eta, etaBins);
        collectNeighborBinsPhi(phiB, nPhi, nN_phi, phiBins);

        matches.clear();

        for (int ptUse : ptBins) {
          for (int etaUse : etaBins) {
            for (int phiUse : phiBins) {
              auto const& vec = buffer[linearKey(colBin, status, ptUse, etaUse, phiUse, mB,
                                                 nStat, nPt, nEta, nPhi, nM)];
              for (auto const& bc : vec) {
                if (bc.collisionIdx == curColIdx) {
                  continue;
                }

                auto tX = V0sMC.iteratorAt(static_cast<uint64_t>(bc.rowIndex));
                if (!selectionV0MC(tX)) {
                  continue;
                }
                if (!checkKinematicsMC(t1, tX)) {
                  continue;
                }

                if (tX.globalIndex() == t1.globalIndex())
                  continue;
                if (tX.globalIndex() == t2.globalIndex())
                  continue;

                matches.push_back(MatchRef{bc.collisionIdx, bc.rowIndex});
              }
            }
          }
        }

        if (matches.empty()) {
          continue;
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
        if (matches.empty()) {
          continue;
        }

        if (cfgV5MaxMatches.value > 0 && (int)matches.size() > cfgV5MaxMatches.value) {
          uint64_t seed = cfgMixSeed.value;
          seed ^= splitmix64((uint64_t)t1.globalIndex());
          seed ^= splitmix64((uint64_t)t2.globalIndex() + 0x1234567ULL);
          seed ^= splitmix64((uint64_t)curColIdx + 0x9abcULL);

          const int K = cfgV5MaxMatches.value;
          for (int i = 0; i < K; ++i) {
            seed = splitmix64(seed);
            const int j = i + (int)(seed % (uint64_t)(matches.size() - i));
            std::swap(matches[i], matches[j]);
          }
          matches.resize(K);
        }

        const float wBase = 1.0f / static_cast<float>(matches.size());

        for (auto const& m : matches) {
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
          const float meWeight = wBase * centPairWeight;
          const float dPhi = deltaPhiMinusPiToPi((float)lX.Phi(), (float)l2.Phi());
          histos.fill(HIST("deltaPhiMix"), dPhi, wBase);
          histos.fill(HIST("hCentPairTypeME"), mcacc::cent(col1), ptype, wBase);
          fillHistograms(mcacc::v0Status(tX), mcacc::v0Status(t2),
                         lX, l2, pX, p2,
                         /*datatype=*/1, /*mixpairweight=*/meWeight);
        }
      }
    }
  }
  PROCESS_SWITCH(lambdaspincorrderived, processMCMEV5, "Process MC ME v5 (paper-style)", false);

  void processMEV6(EventCandidates const& collisions, AllTrackCandidates const& V0s)
  {
    MixBinner mb{
      ptMin.value, ptMax.value, ptMix.value,
      v0etaMixBuffer.value, etaMix.value,
      phiMix.value,
      MassMin.value, MassMax.value, cfgV5MassBins.value};

    const int nCol = colBinning.getAllBinsCount();
    const int nStat = N_STATUS;
    const int nPt = mb.nPt();
    const int nEta = mb.nEta();
    const int nPhi = mb.nPhi();
    const int nM = mb.nM();

    const size_t nKeys = static_cast<size_t>(nCol) * nStat * nPt * nEta * nPhi * nM;
    std::vector<std::vector<BufferCand>> buffer(nKeys);

    // -------- PASS 1: fill buffer --------
    for (auto const& col : collisions) {
      const int colBin = colBinning.getBin(std::make_tuple(col.posz(), col.cent()));
      if (colBin < 0) {
        continue;
      }

      auto slice = V0s.sliceBy(tracksPerCollisionV0, col.index());

      for (auto const& t : slice) {
        if (!selectionV0(t)) {
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

        if (ptB < 0 || etaB < 0 || phiB < 0 || mB < 0) {
          continue;
        }

        const size_t key = linearKey(colBin, status, ptB, etaB, phiB, mB,
                                     nStat, nPt, nEta, nPhi, nM);

        buffer[key].push_back(BufferCand{
          .collisionIdx = static_cast<int64_t>(col.index()),
          .rowIndex = static_cast<int64_t>(t.globalIndex()),
          .v0Status = static_cast<uint8_t>(status),
          .ptBin = static_cast<uint16_t>(ptB),
          .etaBin = static_cast<uint16_t>(etaB),
          .phiBin = static_cast<uint16_t>(phiB),
          .mBin = static_cast<uint16_t>(mB)});
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

      if (ptB < 0 || etaB < 0 || phiB < 0 || mB < 0) {
        return;
      }

      collectNeighborBinsClamp(ptB, nPt, nN_pt, ptBins);
      collectNeighborBinsClamp(etaB, nEta, nN_eta, etaBins);
      collectNeighborBinsPhi(phiB, nPhi, nN_phi, phiBins);

      for (int ptUse : ptBins) {
        for (int etaUse : etaBins) {
          for (int phiUse : phiBins) {
            const auto& vec = buffer[linearKey(colBin, status, ptUse, etaUse, phiUse, mB,
                                               nStat, nPt, nEta, nPhi, nM)];

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

    // -------- PASS 2: two-leg mixing --------
    for (auto const& col1 : collisions) {
      const int colBin = colBinning.getBin(std::make_tuple(col1.posz(), col1.cent()));
      if (colBin < 0) {
        continue;
      }

      const int64_t curColIdx = static_cast<int64_t>(col1.index());
      auto poolA = V0s.sliceBy(tracksPerCollisionV0, col1.index());

      for (auto const& [t1, t2] :
           soa::combinations(o2::soa::CombinationsFullIndexPolicy(poolA, poolA))) {

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

        // leg 1 replaced: (t1,t2) -> (tX,t2)
        collectMatchesForReplacedLeg(t1, t2, colBin, curColIdx, matches1);

        // leg 2 replaced: (t1,t2) -> (t1,tY)
        collectMatchesForReplacedLeg(t2, t1, colBin, curColIdx, matches2);

        downsampleMatches(matches1,
                          (uint64_t)t1.globalIndex() ^ (splitmix64((uint64_t)t2.globalIndex()) + 0x111ULL) ^ splitmix64((uint64_t)curColIdx));
        downsampleMatches(matches2,
                          (uint64_t)t2.globalIndex() ^ (splitmix64((uint64_t)t1.globalIndex()) + 0x222ULL) ^ splitmix64((uint64_t)curColIdx));

        const int nReuse = static_cast<int>(matches1.size() + matches2.size());
        if (nReuse <= 0) {
          continue;
        }

        const float wSE = 1.0f / static_cast<float>(nReuse);

        // replace t1 -> tX, keep t2
        for (auto const& m : matches1) {
          auto tX = V0s.iteratorAt(static_cast<uint64_t>(m.rowIndex));

          auto proton = ROOT::Math::PtEtaPhiMVector(tX.protonPt(), tX.protonEta(), tX.protonPhi(),
                                                    o2::constants::physics::MassProton);
          auto lambda = ROOT::Math::PtEtaPhiMVector(tX.lambdaPt(), tX.lambdaEta(), tX.lambdaPhi(),
                                                    tX.lambdaMass());

          auto proton2 = ROOT::Math::PtEtaPhiMVector(t2.protonPt(), t2.protonEta(), t2.protonPhi(),
                                                     o2::constants::physics::MassProton);
          auto lambda2 = ROOT::Math::PtEtaPhiMVector(t2.lambdaPt(), t2.lambdaEta(), t2.lambdaPhi(),
                                                     t2.lambdaMass());

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

          fillHistograms(tX.v0Status(), t2.v0Status(),
                         lambda, lambda2, proton, proton2,
                         1, meWeight);
        }

        // replace t2 -> tY, keep t1
        for (auto const& m : matches2) {
          auto tY = V0s.iteratorAt(static_cast<uint64_t>(m.rowIndex));

          auto proton = ROOT::Math::PtEtaPhiMVector(t1.protonPt(), t1.protonEta(), t1.protonPhi(),
                                                    o2::constants::physics::MassProton);
          auto lambda = ROOT::Math::PtEtaPhiMVector(t1.lambdaPt(), t1.lambdaEta(), t1.lambdaPhi(),
                                                    t1.lambdaMass());

          auto proton2 = ROOT::Math::PtEtaPhiMVector(tY.protonPt(), tY.protonEta(), tY.protonPhi(),
                                                     o2::constants::physics::MassProton);
          auto lambda2 = ROOT::Math::PtEtaPhiMVector(tY.lambdaPt(), tY.lambdaEta(), tY.lambdaPhi(),
                                                     tY.lambdaMass());

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

          fillHistograms(t1.v0Status(), tY.v0Status(),
                         lambda, lambda2, proton, proton2,
                         1, meWeight);
        }
      }
    }
  }
  PROCESS_SWITCH(lambdaspincorrderived, processMEV6, "Process data ME v6 two-leg", false);

  void processMCMEV6(EventCandidatesMC const& collisions, AllTrackCandidatesMC const& V0sMC)
  {
    MixBinner mb{
      ptMin.value, ptMax.value, ptMix.value,
      v0etaMixBuffer.value, etaMix.value,
      phiMix.value,
      MassMin.value, MassMax.value, cfgV5MassBins.value};

    const int nCol = colBinning.getAllBinsCount();
    const int nStat = N_STATUS;
    const int nPt = mb.nPt();
    const int nEta = mb.nEta();
    const int nPhi = mb.nPhi();
    const int nM = mb.nM();

    const size_t nKeys = static_cast<size_t>(nCol) * nStat * nPt * nEta * nPhi * nM;
    std::vector<std::vector<BufferCand>> buffer(nKeys);

    // -------- PASS 1: fill buffer --------
    for (auto const& col : collisions) {
      const int colBin = colBinning.getBin(std::make_tuple(mcacc::posz(col), mcacc::cent(col)));
      if (colBin < 0) {
        continue;
      }

      auto slice = V0sMC.sliceBy(tracksPerCollisionV0mc, col.index());

      for (auto const& t : slice) {
        if (!selectionV0MC(t)) {
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

        if (ptB < 0 || etaB < 0 || phiB < 0 || mB < 0) {
          continue;
        }

        const size_t key = linearKey(colBin, status, ptB, etaB, phiB, mB,
                                     nStat, nPt, nEta, nPhi, nM);

        buffer[key].push_back(BufferCand{
          .collisionIdx = static_cast<int64_t>(col.index()),
          .rowIndex = static_cast<int64_t>(t.globalIndex()),
          .v0Status = static_cast<uint8_t>(status),
          .ptBin = static_cast<uint16_t>(ptB),
          .etaBin = static_cast<uint16_t>(etaB),
          .phiBin = static_cast<uint16_t>(phiB),
          .mBin = static_cast<uint16_t>(mB)});
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

      if (ptB < 0 || etaB < 0 || phiB < 0 || mB < 0) {
        return;
      }

      collectNeighborBinsClamp(ptB, nPt, nN_pt, ptBins);
      collectNeighborBinsClamp(etaB, nEta, nN_eta, etaBins);
      collectNeighborBinsPhi(phiB, nPhi, nN_phi, phiBins);

      for (int ptUse : ptBins) {
        for (int etaUse : etaBins) {
          for (int phiUse : phiBins) {
            const auto& vec = buffer[linearKey(colBin, status, ptUse, etaUse, phiUse, mB,
                                               nStat, nPt, nEta, nPhi, nM)];

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

    // -------- PASS 2: two-leg mixing --------
    for (auto const& col1 : collisions) {
      const int colBin = colBinning.getBin(std::make_tuple(mcacc::posz(col1), mcacc::cent(col1)));
      if (colBin < 0) {
        continue;
      }

      const int64_t curColIdx = static_cast<int64_t>(col1.index());
      auto poolA = V0sMC.sliceBy(tracksPerCollisionV0mc, col1.index());

      for (auto const& [t1, t2] :
           soa::combinations(o2::soa::CombinationsFullIndexPolicy(poolA, poolA))) {

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

        collectMatchesForReplacedLeg(t1, t2, colBin, curColIdx, matches1);
        collectMatchesForReplacedLeg(t2, t1, colBin, curColIdx, matches2);

        downsampleMatches(matches1,
                          (uint64_t)t1.globalIndex() ^ (splitmix64((uint64_t)t2.globalIndex()) + 0x111ULL) ^ splitmix64((uint64_t)curColIdx));
        downsampleMatches(matches2,
                          (uint64_t)t2.globalIndex() ^ (splitmix64((uint64_t)t1.globalIndex()) + 0x222ULL) ^ splitmix64((uint64_t)curColIdx));

        const int nReuse = static_cast<int>(matches1.size() + matches2.size());
        if (nReuse <= 0) {
          continue;
        }

        const float wSE = 1.0f / static_cast<float>(nReuse);

        // replace t1 -> tX, keep t2
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
                         1, meWeight);
        }

        // replace t2 -> tY, keep t1
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
                         1, meWeight);
        }
      }
    }
  }
  PROCESS_SWITCH(lambdaspincorrderived, processMCMEV6, "Process MC ME v6 two-leg", false);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<lambdaspincorrderived>(cfgc)};
}
