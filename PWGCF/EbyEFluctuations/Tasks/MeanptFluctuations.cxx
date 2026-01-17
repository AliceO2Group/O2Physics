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

/// \file MeanptFluctuations.cxx
/// \brief Task for analyzing <pT> fluctuation upto fourth order of inclusive hadrons
/// \author Swati Saha

#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CommonConstants/MathConstants.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/runDataProcessing.h"
#include <CCDB/BasicCCDBManager.h>

#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TList.h"
#include "TMath.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TRandom3.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <numeric>
#include <string>
#include <vector>

namespace o2::aod
{
namespace pt_qn
{
DECLARE_SOA_COLUMN(Q1, q1, float);                 //! sum of pT of tracks in an event
DECLARE_SOA_COLUMN(Q2, q2, float);                 //! sum of (pT)^2 of tracks in an event
DECLARE_SOA_COLUMN(Q3, q3, float);                 //! sum of (pT)^3 of tracks in an event
DECLARE_SOA_COLUMN(Q4, q4, float);                 //! sum of (pT)^4 of tracks in an event
DECLARE_SOA_COLUMN(Nch, nch, float);               //! no of charged particles/multiplicity in an event
DECLARE_SOA_COLUMN(Centrality, centrality, float); //! Centrality of event
} // namespace pt_qn
DECLARE_SOA_TABLE(MultPtQn, "AOD", "PTQN", pt_qn::Q1, pt_qn::Q2, pt_qn::Q3, pt_qn::Q4, pt_qn::Nch, pt_qn::Centrality); //! table to store e-by-e sum of pT, (pT)^2, (pT)^3, (pT)^4 of tracks, multiplicity and centrality
} // namespace o2::aod

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

struct MeanptFluctuationsQAQnTable {

  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCutPreSelEta{"cfgCutPreSelEta", 0.8f, "|eta|<x; choose x value"};
  Configurable<float> cfgCutPreSelPt{"cfgCutPreSelPt", 5.0f, "Maximum allowed pT"};
  Configurable<float> cfgCutPtLower{"cfgCutPtLower", 0.2f, "Lower pT cut"};
  Configurable<float> cfgCutPtUpper{"cfgCutPtUpper", 3.0f, "Higher pT cut"};
  Configurable<float> cfgCutTpcChi2NCl{"cfgCutTpcChi2NCl", 2.5f, "Maximum TPCchi2NCl"};
  Configurable<float> cfgCutItsChi2NCl{"cfgCutItsChi2NCl", 36.0f, "Maximum ITSchi2NCl"};
  Configurable<float> cfgCutTrackDcaZ{"cfgCutTrackDcaZ", 2.0f, "Maximum DcaZ"};
  Configurable<int> cfgITScluster{"cfgITScluster", 1, "Minimum Number of ITS cluster"};
  Configurable<int> cfgTPCcluster{"cfgTPCcluster", 80, "Minimum Number of TPC cluster"};
  Configurable<int> cfgTPCnCrossedRows{"cfgTPCnCrossedRows", 70, "Minimum Number of TPC crossed-rows"};
  ConfigurableAxis nchAxis{"nchAxis", {500, 0.5, 500.5}, "Axis for multiplicity of GlobalTracks/PVTracks"};
  ConfigurableAxis nchAxis2{"nchAxis2", {1000, 0.5, 30000.5}, "Axis for multiplicity of FT0A/FT0C/FV0A"};
  ConfigurableAxis nchAxis3{"nchAxis3", {1000, 0.5, 100000.5}, "Axis for multiplicity of FT0A/FT0C/FV0A"};
  ConfigurableAxis centAxis{"centAxis", {90, 0., 90.0}, ""};
  Configurable<bool> cfgEvSelkNoSameBunchPileup{"cfgEvSelkNoSameBunchPileup", true, "Pileup removal"};
  Configurable<bool> cfgUseGoodITSLayerAllCut{"cfgUseGoodITSLayerAllCut", true, "Remove time interval with dead ITS zone"};
  Configurable<bool> cfgEvSelkNoITSROFrameBorder{"cfgEvSelkNoITSROFrameBorder", true, "ITSROFrame border event selection cut"};
  Configurable<bool> cfgEvSelkNoTimeFrameBorder{"cfgEvSelkNoTimeFrameBorder", true, "TimeFrame border event selection cut"};
  Configurable<bool> cfgEvSelUseGoodZvtxFT0vsPV{"cfgEvSelUseGoodZvtxFT0vsPV", true, "GoodZvertex and FT0 vs PV cut"};
  Configurable<int> cfgCentralityEstimator{"cfgCentralityEstimator", 1, "Centrlaity estimatore choice: 1-->FT0C, 2-->FT0A; 3-->FT0M, 4-->FV0A"};

  O2_DEFINE_CONFIGURABLE(cfgEvSelMultCorrelation, bool, true, "Multiplicity correlation cut")
  O2_DEFINE_CONFIGURABLE(cfgEvSelV0AT0ACut, bool, true, "V0A T0A 5 sigma cut")
  struct : ConfigurableGroup {
    O2_DEFINE_CONFIGURABLE(cfgMultCentHighCutFunction, std::string, "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x + 10.*([5] + [6]*x + [7]*x*x + [8]*x*x*x + [9]*x*x*x*x)", "Functional for multiplicity correlation cut");
    O2_DEFINE_CONFIGURABLE(cfgMultCentLowCutFunction, std::string, "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x - 3.*([5] + [6]*x + [7]*x*x + [8]*x*x*x + [9]*x*x*x*x)", "Functional for multiplicity correlation cut");
    O2_DEFINE_CONFIGURABLE(cfgMultT0CCutEnabled, bool, false, "Enable Global multiplicity vs T0C centrality cut")
    Configurable<std::vector<double>> cfgMultT0CCutPars{"cfgMultT0CCutPars", std::vector<double>{143.04, -4.58368, 0.0766055, -0.000727796, 2.86153e-06, 23.3108, -0.36304, 0.00437706, -4.717e-05, 1.98332e-07}, "Global multiplicity vs T0C centrality cut parameter values"};
    O2_DEFINE_CONFIGURABLE(cfgMultPVT0CCutEnabled, bool, false, "Enable PV multiplicity vs T0C centrality cut")
    Configurable<std::vector<double>> cfgMultPVT0CCutPars{"cfgMultPVT0CCutPars", std::vector<double>{195.357, -6.15194, 0.101313, -0.000955828, 3.74793e-06, 30.0326, -0.43322, 0.00476265, -5.11206e-05, 2.13613e-07}, "PV multiplicity vs T0C centrality cut parameter values"};

    O2_DEFINE_CONFIGURABLE(cfgMultMultPVHighCutFunction, std::string, "[0]+[1]*x + 5.*([2]+[3]*x)", "Functional for multiplicity correlation cut");
    O2_DEFINE_CONFIGURABLE(cfgMultMultPVLowCutFunction, std::string, "[0]+[1]*x - 5.*([2]+[3]*x)", "Functional for multiplicity correlation cut");
    O2_DEFINE_CONFIGURABLE(cfgMultGlobalPVCutEnabled, bool, false, "Enable global multiplicity vs PV multiplicity cut")
    Configurable<std::vector<double>> cfgMultGlobalPVCutPars{"cfgMultGlobalPVCutPars", std::vector<double>{-0.140809, 0.734344, 2.77495, 0.0165935}, "PV multiplicity vs T0C centrality cut parameter values"};

    O2_DEFINE_CONFIGURABLE(cfgMultMultV0AHighCutFunction, std::string, "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x + 4.*([5] + [6]*x + [7]*x*x + [8]*x*x*x + [9]*x*x*x*x)", "Functional for multiplicity correlation cut");
    O2_DEFINE_CONFIGURABLE(cfgMultMultV0ALowCutFunction, std::string, "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x - 3.*([5] + [6]*x + [7]*x*x + [8]*x*x*x + [9]*x*x*x*x)", "Functional for multiplicity correlation cut");
    O2_DEFINE_CONFIGURABLE(cfgMultMultV0ACutEnabled, bool, false, "Enable global multiplicity vs V0A multiplicity cut")
    Configurable<std::vector<double>> cfgMultMultV0ACutPars{"cfgMultMultV0ACutPars", std::vector<double>{534.893, 184.344, 0.423539, -0.00331436, 5.34622e-06, 871.239, 53.3735, -0.203528, 0.000122758, 5.41027e-07}, "Global multiplicity vs V0A multiplicity cut parameter values"};

    std::vector<double> multT0CCutPars;
    std::vector<double> multPVT0CCutPars;
    std::vector<double> multGlobalPVCutPars;
    std::vector<double> multMultV0ACutPars;
    TF1* fMultPVT0CCutLow = nullptr;
    TF1* fMultPVT0CCutHigh = nullptr;
    TF1* fMultT0CCutLow = nullptr;
    TF1* fMultT0CCutHigh = nullptr;
    TF1* fMultGlobalPVCutLow = nullptr;
    TF1* fMultGlobalPVCutHigh = nullptr;
    TF1* fMultMultV0ACutLow = nullptr;
    TF1* fMultMultV0ACutHigh = nullptr;
    TF1* fT0AV0AMean = nullptr;
    TF1* fT0AV0ASigma = nullptr;

  } cfgFuncParas;

  O2_DEFINE_CONFIGURABLE(cfgUse22sEventCut, bool, true, "Use 22s event cut on mult correlations")
  O2_DEFINE_CONFIGURABLE(cfgUseSmallIonAdditionalEventCut, bool, true, "Use additional event cut on mult correlations for small ions")

  // Filter command***********
  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter trackFilter = (nabs(aod::track::eta) < cfgCutPreSelEta) && (aod::track::pt > cfgCutPtLower) && (aod::track::pt < cfgCutPreSelPt) && (requireGlobalTrackInFilter()) && (aod::track::tpcChi2NCl < cfgCutTpcChi2NCl) && (aod::track::itsChi2NCl < cfgCutItsChi2NCl) && (nabs(aod::track::dcaZ) < cfgCutTrackDcaZ);

  // Connect to ccdb
  Service<ccdb::BasicCCDBManager> ccdb;
  Configurable<int64_t> ccdbnolaterthan{"ccdbnolaterthan", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};
  Configurable<std::string> ccdburl{"ccdburl", "http://ccdb-test.cern.ch:8080", "url of the ccdb repository"};

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // filtering collisions and tracks***********
  using AodCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::CentFT0Ms, aod::CentFT0As, aod::CentFV0As, aod::Mults>>;
  // using aodCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>;
  using AodTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA>>;

  // Event selection cuts - Alex
  TF1* fMultPVCutLow = nullptr;
  TF1* fMultPVCutHigh = nullptr;
  TF1* fMultCutLow = nullptr;
  TF1* fMultCutHigh = nullptr;
  TF1* fMultMultPVCut = nullptr;

  // Equivalent of the AliRoot task UserCreateOutputObjects
  void init(o2::framework::InitContext&)
  {
    // Define your axes
    // Constant bin width axis
    AxisSpec vtxZAxis = {100, -20, 20};
    // Variable bin width axis
    std::vector<double> ptBinning = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.8, 3.2, 3.6, 4.};
    AxisSpec ptAxis = {ptBinning, "#it{p}_{T} (GeV/#it{c})"};

    // Add histograms to histogram manager (as in the output object of in AliPhysics)
    histos.add("hZvtx_after_sel", ";Z (cm)", kTH1F, {vtxZAxis});
    histos.add("hP", ";#it{p} (GeV/#it{c})", kTH1F, {{35, 0.2, 4.}});
    histos.add("hPt", ";#it{p}_{T} (GeV/#it{c})", kTH1F, {ptAxis});
    histos.add("hPhi", ";#phi", kTH1F, {{100, 0., o2::constants::math::TwoPI}});
    histos.add("hEta", ";#eta", kTH1F, {{100, -2.01, 2.01}});
    histos.add("hCentrality", ";centrality (%)", kTH1F, {{90, 0, 90}});
    histos.add("hDcaXY", ";#it{dca}_{XY}", kTH1F, {{1000, -5, 5}});
    histos.add("hDcaZ", ";#it{dca}_{Z}", kTH1F, {{1000, -5, 5}});
    histos.add("hMeanPt", "", kTProfile, {centAxis});
    histos.add("Hist2D_globalTracks_PVTracks", "", {HistType::kTH2D, {nchAxis, nchAxis}});
    histos.add("Hist2D_cent_nch", "", {HistType::kTH2D, {nchAxis, centAxis}});
    // before selection
    histos.add("MultCorrelationPlots/BeforeSelection/His2D_globalTracks_PVTracks_beforeSel", "", {HistType::kTH2D, {nchAxis, nchAxis}});
    histos.add("MultCorrelationPlots/BeforeSelection/His2D_globalTracks_centFT0C_beforeSel", "", {HistType::kTH2D, {centAxis, nchAxis}});
    histos.add("MultCorrelationPlots/BeforeSelection/His2D_PVTracks_centFT0C_beforeSel", "", {HistType::kTH2D, {centAxis, nchAxis}});
    histos.add("MultCorrelationPlots/BeforeSelection/His2D_globalTracks_V0ATracks_beforeSel", "", {HistType::kTH2D, {nchAxis3, nchAxis}});
    histos.add("MultCorrelationPlots/BeforeSelection/His2D_globalTracks_T0ATracks_beforeSel", "", {HistType::kTH2D, {nchAxis2, nchAxis}});
    histos.add("MultCorrelationPlots/BeforeSelection/His2D_V0ATracks_T0CTracks_beforeSel", "", {HistType::kTH2D, {nchAxis2, nchAxis3}});
    // after selection
    if (cfgUseSmallIonAdditionalEventCut) {
      histos.add("MultCorrelationPlots/AfterSelection/His2D_globalTracks_PVTracks_afterSel", "", {HistType::kTH2D, {nchAxis, nchAxis}});
      histos.add("MultCorrelationPlots/AfterSelection/His2D_globalTracks_centFT0C_afterSel", "", {HistType::kTH2D, {centAxis, nchAxis}});
      histos.add("MultCorrelationPlots/AfterSelection/His2D_PVTracks_centFT0C_afterSel", "", {HistType::kTH2D, {centAxis, nchAxis}});
      histos.add("MultCorrelationPlots/AfterSelection/His2D_globalTracks_V0ATracks_afterSel", "", {HistType::kTH2D, {nchAxis3, nchAxis}});
      histos.add("MultCorrelationPlots/AfterSelection/His2D_globalTracks_T0ATracks_afterSel", "", {HistType::kTH2D, {nchAxis2, nchAxis}});
      histos.add("MultCorrelationPlots/AfterSelection/His2D_V0ATracks_T0CTracks_afterSel", "", {HistType::kTH2D, {nchAxis2, nchAxis3}});
    }

    // Event selection - Alex
    if (cfgUse22sEventCut) {
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

    if (cfgEvSelMultCorrelation) {
      cfgFuncParas.multT0CCutPars = cfgFuncParas.cfgMultT0CCutPars;
      cfgFuncParas.multPVT0CCutPars = cfgFuncParas.cfgMultPVT0CCutPars;
      cfgFuncParas.multGlobalPVCutPars = cfgFuncParas.cfgMultGlobalPVCutPars;
      cfgFuncParas.multMultV0ACutPars = cfgFuncParas.cfgMultMultV0ACutPars;
      cfgFuncParas.fMultPVT0CCutLow = new TF1("fMultPVT0CCutLow", cfgFuncParas.cfgMultCentLowCutFunction->c_str(), 0, 100);
      cfgFuncParas.fMultPVT0CCutLow->SetParameters(&(cfgFuncParas.multPVT0CCutPars[0]));
      cfgFuncParas.fMultPVT0CCutHigh = new TF1("fMultPVT0CCutHigh", cfgFuncParas.cfgMultCentHighCutFunction->c_str(), 0, 100);
      cfgFuncParas.fMultPVT0CCutHigh->SetParameters(&(cfgFuncParas.multPVT0CCutPars[0]));
      cfgFuncParas.fMultT0CCutLow = new TF1("fMultT0CCutLow", cfgFuncParas.cfgMultCentLowCutFunction->c_str(), 0, 100);
      cfgFuncParas.fMultT0CCutLow->SetParameters(&(cfgFuncParas.multT0CCutPars[0]));
      cfgFuncParas.fMultT0CCutHigh = new TF1("fMultT0CCutHigh", cfgFuncParas.cfgMultCentHighCutFunction->c_str(), 0, 100);
      cfgFuncParas.fMultT0CCutHigh->SetParameters(&(cfgFuncParas.multT0CCutPars[0]));
      cfgFuncParas.fMultGlobalPVCutLow = new TF1("fMultGlobalPVCutLow", cfgFuncParas.cfgMultMultPVLowCutFunction->c_str(), 0, 4000);
      cfgFuncParas.fMultGlobalPVCutLow->SetParameters(&(cfgFuncParas.multGlobalPVCutPars[0]));
      cfgFuncParas.fMultGlobalPVCutHigh = new TF1("fMultGlobalPVCutHigh", cfgFuncParas.cfgMultMultPVHighCutFunction->c_str(), 0, 4000);
      cfgFuncParas.fMultGlobalPVCutHigh->SetParameters(&(cfgFuncParas.multGlobalPVCutPars[0]));
      cfgFuncParas.fMultMultV0ACutLow = new TF1("fMultMultV0ACutLow", cfgFuncParas.cfgMultMultV0ALowCutFunction->c_str(), 0, 4000);
      cfgFuncParas.fMultMultV0ACutLow->SetParameters(&(cfgFuncParas.multMultV0ACutPars[0]));
      cfgFuncParas.fMultMultV0ACutHigh = new TF1("fMultMultV0ACutHigh", cfgFuncParas.cfgMultMultV0AHighCutFunction->c_str(), 0, 4000);
      cfgFuncParas.fMultMultV0ACutHigh->SetParameters(&(cfgFuncParas.multMultV0ACutPars[0]));
      cfgFuncParas.fT0AV0AMean = new TF1("fT0AV0AMean", "[0]+[1]*x", 0, 200000);
      cfgFuncParas.fT0AV0AMean->SetParameters(-1601.0581, 9.417652e-01);
      cfgFuncParas.fT0AV0ASigma = new TF1("fT0AV0ASigma", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x", 0, 200000);
      cfgFuncParas.fT0AV0ASigma->SetParameters(463.4144, 6.796509e-02, -9.097136e-07, 7.971088e-12, -2.600581e-17);
    }

  } //! end init function

  template <typename TCollision>
  bool eventSelected(TCollision collision, const int& multTrk, const float& centrality)
  {
    if (collision.alias_bit(kTVXinTRD)) {
      // TRD triggered
      return 0;
    }
    float vtxz = -999;
    if (collision.numContrib() > 1) {
      vtxz = collision.posZ();
      float zRes = std::sqrt(collision.covZZ());
      float zResMax = 0.25;
      float numContribMin = 20;
      if (zRes > zResMax && collision.numContrib() < numContribMin)
        vtxz = -999;
    }
    auto multNTracksPV = collision.multNTracksPV();

    if ((vtxz > cfgCutVertex) || (vtxz < -1.0 * cfgCutVertex))
      return 0;
    if (multNTracksPV < fMultPVCutLow->Eval(centrality))
      return 0;
    if (multNTracksPV > fMultPVCutHigh->Eval(centrality))
      return 0;
    if (multTrk < fMultCutLow->Eval(centrality))
      return 0;
    if (multTrk > fMultCutHigh->Eval(centrality))
      return 0;
    if (multTrk > fMultMultPVCut->Eval(multNTracksPV))
      return 0;

    return 1;
  }

  template <typename TCollision>
  bool eventSelectedSmallion(TCollision collision, const int multTrk, const float centrality)
  {
    auto multNTracksPV = collision.multNTracksPV();

    if (cfgEvSelMultCorrelation) {
      if (cfgFuncParas.cfgMultPVT0CCutEnabled) {
        if (multNTracksPV < cfgFuncParas.fMultPVT0CCutLow->Eval(centrality))
          return 0;
        if (multNTracksPV > cfgFuncParas.fMultPVT0CCutHigh->Eval(centrality))
          return 0;
      }

      if (cfgFuncParas.cfgMultT0CCutEnabled) {
        if (multTrk < cfgFuncParas.fMultT0CCutLow->Eval(centrality))
          return 0;
        if (multTrk > cfgFuncParas.fMultT0CCutHigh->Eval(centrality))
          return 0;
      }

      if (cfgFuncParas.cfgMultGlobalPVCutEnabled) {
        if (multTrk < cfgFuncParas.fMultGlobalPVCutLow->Eval(multNTracksPV))
          return 0;
        if (multTrk > cfgFuncParas.fMultGlobalPVCutHigh->Eval(multNTracksPV))
          return 0;
      }

      if (cfgFuncParas.cfgMultMultV0ACutEnabled) {
        if (collision.multFV0A() < cfgFuncParas.fMultMultV0ACutLow->Eval(multTrk))
          return 0;
        if (collision.multFV0A() > cfgFuncParas.fMultMultV0ACutHigh->Eval(multTrk))
          return 0;
      }
    }

    float sigma = 5.0;
    if (cfgEvSelV0AT0ACut && (std::fabs(collision.multFV0A() - cfgFuncParas.fT0AV0AMean->Eval(collision.multFT0A())) > sigma * cfgFuncParas.fT0AV0ASigma->Eval(collision.multFT0A())))
      return 0;

    return 1;
  }

  Produces<aod::MultPtQn> multPtQn;

  // void process(aod::Collision const& coll, aod::Tracks const& inputTracks)
  void process(AodCollisions::iterator const& coll, aod::BCsWithTimestamps const&, AodTracks const& inputTracks)
  {
    if (!coll.sel8()) {
      return;
    }
    if (cfgUseGoodITSLayerAllCut && !(coll.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll))) {
      return;
    }
    if (cfgEvSelkNoSameBunchPileup && !(coll.selection_bit(o2::aod::evsel::kNoSameBunchPileup))) {
      return;
    }
    if (cfgEvSelkNoITSROFrameBorder && !(coll.selection_bit(o2::aod::evsel::kNoITSROFrameBorder))) {
      return;
    }
    if (cfgEvSelkNoTimeFrameBorder && !(coll.selection_bit(o2::aod::evsel::kNoTimeFrameBorder))) {
      return;
    }
    if (cfgEvSelUseGoodZvtxFT0vsPV && !(coll.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV))) {
      return;
    }

    histos.fill(HIST("MultCorrelationPlots/BeforeSelection/His2D_globalTracks_PVTracks_beforeSel"), coll.multNTracksPV(), inputTracks.size());
    histos.fill(HIST("MultCorrelationPlots/BeforeSelection/His2D_globalTracks_centFT0C_beforeSel"), coll.centFT0C(), inputTracks.size());
    histos.fill(HIST("MultCorrelationPlots/BeforeSelection/His2D_PVTracks_centFT0C_beforeSel"), coll.centFT0C(), coll.multNTracksPV());
    histos.fill(HIST("MultCorrelationPlots/BeforeSelection/His2D_globalTracks_V0ATracks_beforeSel"), coll.multFV0A(), inputTracks.size());
    histos.fill(HIST("MultCorrelationPlots/BeforeSelection/His2D_globalTracks_T0ATracks_beforeSel"), coll.multFT0A(), inputTracks.size());
    histos.fill(HIST("MultCorrelationPlots/BeforeSelection/His2D_V0ATracks_T0CTracks_beforeSel"), coll.multFT0C(), coll.multFV0A());

    const auto centralityFT0C = coll.centFT0C();
    if (cfgUse22sEventCut && !eventSelected(coll, inputTracks.size(), centralityFT0C))
      return;
    if (cfgUseSmallIonAdditionalEventCut && !eventSelectedSmallion(coll, inputTracks.size(), centralityFT0C))
      return;

    if (cfgUseSmallIonAdditionalEventCut) {
      histos.fill(HIST("MultCorrelationPlots/AfterSelection/His2D_globalTracks_PVTracks_afterSel"), coll.multNTracksPV(), inputTracks.size());
      histos.fill(HIST("MultCorrelationPlots/AfterSelection/His2D_globalTracks_centFT0C_afterSel"), coll.centFT0C(), inputTracks.size());
      histos.fill(HIST("MultCorrelationPlots/AfterSelection/His2D_PVTracks_centFT0C_afterSel"), coll.centFT0C(), coll.multNTracksPV());
      histos.fill(HIST("MultCorrelationPlots/AfterSelection/His2D_globalTracks_V0ATracks_afterSel"), coll.multFV0A(), inputTracks.size());
      histos.fill(HIST("MultCorrelationPlots/AfterSelection/His2D_globalTracks_T0ATracks_afterSel"), coll.multFT0A(), inputTracks.size());
      histos.fill(HIST("MultCorrelationPlots/AfterSelection/His2D_V0ATracks_T0CTracks_afterSel"), coll.multFT0C(), coll.multFV0A());
    }

    histos.fill(HIST("hZvtx_after_sel"), coll.posZ());

    double cent = 0.0;
    int centChoiceFT0C = 1;
    int centChoiceFT0A = 2;
    int centChoiceFT0M = 3;
    int centChoiceFV0A = 4;
    if (cfgCentralityEstimator == centChoiceFT0C)
      cent = coll.centFT0C();
    else if (cfgCentralityEstimator == centChoiceFT0A)
      cent = coll.centFT0A();
    else if (cfgCentralityEstimator == centChoiceFT0M)
      cent = coll.centFT0M();
    else if (cfgCentralityEstimator == centChoiceFV0A)
      cent = coll.centFV0A();

    histos.fill(HIST("hCentrality"), cent);

    histos.fill(HIST("Hist2D_globalTracks_PVTracks"), coll.multNTracksPV(), inputTracks.size());
    histos.fill(HIST("Hist2D_cent_nch"), inputTracks.size(), centralityFT0C);

    // variables
    double pTsum = 0.0;
    double nN = 0.0;

    float q1 = 0.0;
    float q2 = 0.0;
    float q3 = 0.0;
    float q4 = 0.0;
    float nCh = 0.0;

    for (const auto& track : inputTracks) { // Loop over tracks

      if (!track.has_collision()) {
        continue;
      }

      if (!track.isPVContributor()) {
        continue;
      }

      if (!(track.itsNCls() > cfgITScluster) || !(track.tpcNClsFound() >= cfgTPCcluster) || !(track.tpcNClsCrossedRows() >= cfgTPCnCrossedRows)) {
        continue;
      }

      histos.fill(HIST("hP"), track.p());
      histos.fill(HIST("hPt"), track.pt());
      histos.fill(HIST("hEta"), track.eta());
      histos.fill(HIST("hPhi"), track.phi());
      histos.fill(HIST("hDcaXY"), track.dcaXY());
      histos.fill(HIST("hDcaZ"), track.dcaZ());

      pTsum += track.pt();
      nN += 1.0;

      float pT = track.pt();
      // calculating Q1, Q2, Q3, Q4. Nch
      if (track.pt() > cfgCutPtLower && track.pt() < cfgCutPtUpper && track.sign() != 0) {
        q1 = q1 + std::pow(pT, 1.0);
        q2 = q2 + std::pow(pT, 2.0);
        q3 = q3 + std::pow(pT, 3.0);
        q4 = q4 + std::pow(pT, 4.0);
        nCh = nCh + 1;
      }
    }
    multPtQn(q1, q2, q3, q4, nCh, cent);
    // MeanPt
    if (nN > 0.0f)
      histos.fill(HIST("hMeanPt"), cent, pTsum / nN);
  }
};

struct MeanptFluctuationsAnalysis {

  Configurable<int> cfgNsubSample{"cfgNsubSample", 10, "Number of subsamples"};
  ConfigurableAxis centAxis{"centAxis", {90, 0, 90}, ""};
  ConfigurableAxis multAxis{"multAxis", {5000, 0.5, 5000.5}, ""};
  ConfigurableAxis meanpTAxis{"meanpTAxis", {500, 0, 5.0}, ""};

  float minNch = 3.0f;
  expressions::Filter nchFilter = aod::pt_qn::nch > minNch;
  using FilteredMultPtQn = soa::Filtered<aod::MultPtQn>;

  // Connect to ccdb
  Service<ccdb::BasicCCDBManager> ccdb;
  Configurable<int64_t> ccdbnolaterthan{"ccdbnolaterthan", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};
  Configurable<std::string> ccdburl{"ccdburl", "http://ccdb-test.cern.ch:8080", "url of the ccdb repository"};

  // Define output
  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};
  std::vector<std::vector<std::shared_ptr<TProfile2D>>> subSample;
  TRandom3* fRndm = new TRandom3(0);

  void init(o2::framework::InitContext&)
  {
    // AxisSpec centAxis = {90, 0, 90, "centrality (%)"};
    // AxisSpec multAxis = {5000, 0.5, 5000.5, "#it{N}_{ch,acc}"};

    registry.add("Prof_mean_t1", "", {HistType::kTProfile2D, {centAxis, multAxis}});
    registry.add("Prof_var_t1", "", {HistType::kTProfile2D, {centAxis, multAxis}});
    registry.add("Prof_skew_t1", "", {HistType::kTProfile2D, {centAxis, multAxis}});
    registry.add("Prof_kurt_t1", "", {HistType::kTProfile2D, {centAxis, multAxis}});
    registry.add("Hist2D_Nch_centrality", "", {HistType::kTH2D, {centAxis, multAxis}});
    registry.add("Hist2D_meanpt_centrality", "", {HistType::kTH2D, {centAxis, meanpTAxis}});

    // initial array
    subSample.resize(cfgNsubSample);
    for (int i = 0; i < cfgNsubSample; i++) {
      subSample[i].resize(4);
    }
    for (int i = 0; i < cfgNsubSample; i++) {
      subSample[i][0] = std::get<std::shared_ptr<TProfile2D>>(registry.add(Form("subSample_%d/Prof_mean_t1", i), "", {HistType::kTProfile2D, {centAxis, multAxis}}));
      subSample[i][1] = std::get<std::shared_ptr<TProfile2D>>(registry.add(Form("subSample_%d/Prof_var_t1", i), "", {HistType::kTProfile2D, {centAxis, multAxis}}));
      subSample[i][2] = std::get<std::shared_ptr<TProfile2D>>(registry.add(Form("subSample_%d/Prof_skew_t1", i), "", {HistType::kTProfile2D, {centAxis, multAxis}}));
      subSample[i][3] = std::get<std::shared_ptr<TProfile2D>>(registry.add(Form("subSample_%d/Prof_kurt_t1", i), "", {HistType::kTProfile2D, {centAxis, multAxis}}));
    }
  }

  float meanTerm1;
  float varianceTerm1;
  float skewnessTerm1;
  float kurtosisTerm1;

  // void process(aod::MultPtQn::iterator const& event_ptqn)
  void process(FilteredMultPtQn::iterator const& event_ptqn)
  {
    // LOGF(info, "Centrality= %f Nch= %f Q1= %f Q2= %f", event_ptqn.centrality(), event_ptqn.nch(), event_ptqn.q1(), event_ptqn.q2());

    // calculating observables
    meanTerm1 = event_ptqn.q1() / event_ptqn.nch();
    varianceTerm1 = (std::pow(event_ptqn.q1(), 2.0f) - event_ptqn.q2()) / (event_ptqn.nch() * (event_ptqn.nch() - 1.0f));
    skewnessTerm1 = (std::pow(event_ptqn.q1(), 3.0f) - 3.0f * event_ptqn.q2() * event_ptqn.q1() + 2.0f * event_ptqn.q3()) / (event_ptqn.nch() * (event_ptqn.nch() - 1.0f) * (event_ptqn.nch() - 2.0f));
    kurtosisTerm1 = (std::pow(event_ptqn.q1(), 4.0f) - (6.0f * event_ptqn.q4()) + (8.0f * event_ptqn.q1() * event_ptqn.q3()) - (6.0f * std::pow(event_ptqn.q1(), 2.0f) * event_ptqn.q2()) + (3.0f * std::pow(event_ptqn.q2(), 2.0f))) / (event_ptqn.nch() * (event_ptqn.nch() - 1.0f) * (event_ptqn.nch() - 2.0f) * (event_ptqn.nch() - 3.0f));

    // filling profiles and histograms for central values
    registry.get<TProfile2D>(HIST("Prof_mean_t1"))->Fill(event_ptqn.centrality(), event_ptqn.nch(), meanTerm1);
    registry.get<TProfile2D>(HIST("Prof_var_t1"))->Fill(event_ptqn.centrality(), event_ptqn.nch(), varianceTerm1);
    registry.get<TProfile2D>(HIST("Prof_skew_t1"))->Fill(event_ptqn.centrality(), event_ptqn.nch(), skewnessTerm1);
    registry.get<TProfile2D>(HIST("Prof_kurt_t1"))->Fill(event_ptqn.centrality(), event_ptqn.nch(), kurtosisTerm1);
    registry.fill(HIST("Hist2D_Nch_centrality"), event_ptqn.centrality(), event_ptqn.nch());
    registry.fill(HIST("Hist2D_meanpt_centrality"), event_ptqn.centrality(), meanTerm1);

    // selecting subsample and filling profiles
    float lRandom = fRndm->Rndm();
    int sampleIndex = static_cast<int>(cfgNsubSample * lRandom);
    subSample[sampleIndex][0]->Fill(event_ptqn.centrality(), event_ptqn.nch(), meanTerm1);
    subSample[sampleIndex][1]->Fill(event_ptqn.centrality(), event_ptqn.nch(), varianceTerm1);
    subSample[sampleIndex][2]->Fill(event_ptqn.centrality(), event_ptqn.nch(), skewnessTerm1);
    subSample[sampleIndex][3]->Fill(event_ptqn.centrality(), event_ptqn.nch(), kurtosisTerm1);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  // Equivalent to the AddTask in AliPhysics
  return WorkflowSpec{
    adaptAnalysisTask<MeanptFluctuationsQAQnTable>(cfgc),
    adaptAnalysisTask<MeanptFluctuationsAnalysis>(cfgc),
  };
}
