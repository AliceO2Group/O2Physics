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

#ifndef PWGLF_UTILS_SIGMA0BUILDERHELPER_H_
#define PWGLF_UTILS_SIGMA0BUILDERHELPER_H_

#include <cstdlib>
#include <cmath>
#include <array>
#include "Framework/AnalysisDataModel.h"
#include "ReconstructionDataFormats/Track.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Common/Core/trackUtilities.h"

namespace o2
{
namespace pwglf
{
namespace sigma0
{

// event selection configurables
struct evselConfigurables : o2::framework::ConfigurableGroup {
  std::string prefix = "evselOpts";

  o2::framework::Configurable<bool> doPPAnalysis{"doPPAnalysis", true, "if in pp, set to true"};
  o2::framework::Configurable<bool> requireSel8{"requireSel8", true, "require sel8 event selection"};
  o2::framework::Configurable<bool> requireTriggerTVX{"requireTriggerTVX", true, "require FT0 vertex (acceptable FT0C-FT0A time difference) at trigger level"};
  o2::framework::Configurable<bool> rejectITSROFBorder{"rejectITSROFBorder", true, "reject events at ITS ROF border"};
  o2::framework::Configurable<bool> rejectTFBorder{"rejectTFBorder", true, "reject events at TF border"};
  o2::framework::Configurable<bool> requireIsVertexITSTPC{"requireIsVertexITSTPC", true, "require events with at least one ITS-TPC track"};
  o2::framework::Configurable<bool> requireIsGoodZvtxFT0VsPV{"requireIsGoodZvtxFT0VsPV", true, "require events with PV position along z consistent (within 1 cm) between PV reconstructed using tracks and PV using FT0 A-C time difference"};
  o2::framework::Configurable<bool> requireIsVertexTOFmatched{"requireIsVertexTOFmatched", false, "require events with at least one of vertex contributors matched to TOF"};
  o2::framework::Configurable<bool> requireIsVertexTRDmatched{"requireIsVertexTRDmatched", false, "require events with at least one of vertex contributors matched to TRD"};
  o2::framework::Configurable<bool> rejectSameBunchPileup{"rejectSameBunchPileup", false, "reject collisions in case of pileup with another collision in the same foundBC"};
  o2::framework::Configurable<bool> requireNoCollInTimeRangeStd{"requireNoCollInTimeRangeStd", false, "reject collisions corrupted by the cannibalism, with other collisions within +/- 2 microseconds or mult above a certain threshold in -4 - -2 microseconds"};
  o2::framework::Configurable<bool> requireNoCollInTimeRangeStrict{"requireNoCollInTimeRangeStrict", false, "reject collisions corrupted by the cannibalism, with other collisions within +/- 10 microseconds"};
  o2::framework::Configurable<bool> requireNoCollInTimeRangeNarrow{"requireNoCollInTimeRangeNarrow", false, "reject collisions corrupted by the cannibalism, with other collisions within +/- 2 microseconds"};
  o2::framework::Configurable<bool> requireNoCollInTimeRangeVzDep{"requireNoCollInTimeRangeVzDep", false, "reject collisions corrupted by the cannibalism, with other collisions with pvZ of drifting TPC tracks from past/future collisions within 2.5 cm the current pvZ"};
  o2::framework::Configurable<bool> requireNoCollInROFStd{"requireNoCollInROFStd", false, "reject collisions corrupted by the cannibalism, with other collisions within the same ITS ROF with mult. above a certain threshold"};
  o2::framework::Configurable<bool> requireNoCollInROFStrict{"requireNoCollInROFStrict", false, "reject collisions corrupted by the cannibalism, with other collisions within the same ITS ROF"};
  o2::framework::Configurable<bool> requireINEL0{"requireINEL0", false, "require INEL>0 event selection"};
  o2::framework::Configurable<bool> requireINEL1{"requireINEL1", false, "require INEL>1 event selection"};
  o2::framework::Configurable<float> maxZVtxPosition{"maxZVtxPosition", 10., "max Z vtx position"};
  o2::framework::Configurable<bool> useEvtSelInDenomEff{"useEvtSelInDenomEff", false, "Consider event selections in the recoed <-> gen collision association for the denominator (or numerator) of the acc. x eff. (or signal loss)?"};
  o2::framework::Configurable<bool> applyZVtxSelOnMCPV{"applyZVtxSelOnMCPV", false, "Apply Z-vtx cut on the PV of the generated collision?"};
  o2::framework::Configurable<bool> useFT0CbasedOccupancy{"useFT0CbasedOccupancy", false, "Use sum of FT0-C amplitudes for estimating occupancy? (if not, use track-based definition)"};
  // fast check on occupancy
  o2::framework::Configurable<float> minOccupancy{"minOccupancy", -1, "minimum occupancy from neighbouring collisions"};
  o2::framework::Configurable<float> maxOccupancy{"maxOccupancy", -1, "maximum occupancy from neighbouring collisions"};

  // fast check on interaction rate
  o2::framework::Configurable<bool> fGetIR{"fGetIR", false, "Flag to retrieve the IR info."};
  o2::framework::Configurable<bool> fIRCrashOnNull{"fIRCrashOnNull", false, "Flag to avoid CTP RateFetcher crash."};
  o2::framework::Configurable<std::string> irSource{"irSource", "T0VTX", "Estimator of the interaction rate (Recommended: pp --> T0VTX, Pb-Pb --> ZNC hadronic)"};
  o2::framework::Configurable<float> minIR{"minIR", -1, "minimum IR collisions"};
  o2::framework::Configurable<float> maxIR{"maxIR", -1, "maximum IR collisions"};

};

// Lambda criteria:
struct lambdaselConfigurables : o2::framework::ConfigurableGroup {
  std::string prefix = "lambdaselOpts";

  o2::framework::Configurable<bool> useMLScores{"useMLScores", false, "use ML scores to select candidates"};
  o2::framework::Configurable<float> Lambda_MLThreshold{"Lambda_MLThreshold", 0.1, "Decision Threshold value to select lambdas"};
  o2::framework::Configurable<float> AntiLambda_MLThreshold{"AntiLambda_MLThreshold", 0.1, "Decision Threshold value to select antilambdas"};
  o2::framework::Configurable<float> LambdaRapidity{"LambdaRapidity", 0.8, "v0 rapidity"};
  o2::framework::Configurable<float> LambdaDauPseudoRap{"LambdaDauPseudoRap", 1.5, "Max pseudorapidity of daughter tracks"};
  o2::framework::Configurable<float> LambdaMinDCANegToPv{"LambdaMinDCANegToPv", 0.0, "min DCA Neg To PV (cm)"};
  o2::framework::Configurable<float> LambdaMinDCAPosToPv{"LambdaMinDCAPosToPv", 0.0, "min DCA Pos To PV (cm)"};
  o2::framework::Configurable<float> LambdaMaxDCAV0Dau{"LambdaMaxDCAV0Dau", 3.5, "Max DCA V0 Daughters (cm)"};
  o2::framework::Configurable<float> LambdaMinv0radius{"LambdaMinv0radius", 0.0, "Min V0 radius (cm)"};
  o2::framework::Configurable<float> LambdaMaxv0radius{"LambdaMaxv0radius", 60, "Max V0 radius (cm)"};
  o2::framework::Configurable<float> LambdaWindow{"LambdaWindow", 0.05, "Mass window around expected (in GeV/c2)"};
};

struct photonselConfigurables : o2::framework::ConfigurableGroup {
  std::string prefix = "photonselOpts";

  o2::framework::Configurable<bool> useMLScores{"useMLScores", false, "use ML scores to select candidates"};
  o2::framework::Configurable<float> Gamma_MLThreshold{"Gamma_MLThreshold", 0.1, "Decision Threshold value to select gammas"};
  o2::framework::Configurable<float> PhotonRapidity{"PhotonRapidity", 0.8, "v0 rapidity"};
  o2::framework::Configurable<float> PhotonMaxDauPseudoRap{"PhotonMaxDauPseudoRap", 1.5, "Max pseudorapidity of daughter tracks"};
  o2::framework::Configurable<float> PhotonMinDCAToPv{"PhotonMinDCAToPv", 0.0, "Min DCA daughter To PV (cm)"};
  o2::framework::Configurable<float> PhotonMaxDCAV0Dau{"PhotonMaxDCAV0Dau", 3.5, "Max DCA V0 Daughters (cm)"};
  o2::framework::Configurable<float> PhotonMinRadius{"PhotonMinRadius", 0.0, "Min photon conversion radius (cm)"};
  o2::framework::Configurable<float> PhotonMaxRadius{"PhotonMaxRadius", 240, "Max photon conversion radius (cm)"};
  o2::framework::Configurable<float> PhotonMaxMass{"PhotonMaxMass", 0.3, "Max photon mass (GeV/c^{2})"};
};

struct sigma0selConfigurables : o2::framework::ConfigurableGroup {
  std::string prefix = "sigma0selOpts";

  o2::framework::Configurable<float> Sigma0Window{"Sigma0Window", 0.1, "Mass window around expected (in GeV/c2)"};
  o2::framework::Configurable<float> SigmaMaxRap{"SigmaMaxRap", 0.8, "Max sigma0 rapidity"};

};

struct pi0selConfigurables : o2::framework::ConfigurableGroup {
  std::string prefix = "pi0selOpts";

  o2::framework::Configurable<bool> doPi0QA{"doPi0QA", true, "Flag to fill QA histos for pi0 rejection study."};  
  o2::framework::Configurable<float> Pi0PhotonMinDCADauToPv{"Pi0PhotonMinDCADauToPv", 0.0, "Min DCA daughter To PV (cm)"};
  o2::framework::Configurable<float> Pi0PhotonMaxDCAV0Dau{"Pi0PhotonMaxDCAV0Dau", 3.5, "Max DCA V0 Daughters (cm)"};
  o2::framework::Configurable<int> Pi0PhotonMinTPCCrossedRows{"Pi0PhotonMinTPCCrossedRows", 0, "Min daughter TPC Crossed Rows"};
  o2::framework::Configurable<int> Pi0PhotonMaxTPCNSigmas{"Pi0PhotonMaxTPCNSigmas", 7, "Max TPC NSigmas for daughters"};
  o2::framework::Configurable<float> Pi0PhotonMaxEta{"Pi0PhotonMaxEta", 0.8, "Max photon rapidity"};
  o2::framework::Configurable<float> Pi0PhotonMinRadius{"Pi0PhotonMinRadius", 3.0, "Min photon conversion radius (cm)"};
  o2::framework::Configurable<float> Pi0PhotonMaxRadius{"Pi0PhotonMaxRadius", 115, "Max photon conversion radius (cm)"};
  o2::framework::Configurable<float> Pi0PhotonMaxQt{"Pi0PhotonMaxQt", 0.05, "Max photon qt value (AP plot) (GeV/c)"};
  o2::framework::Configurable<float> Pi0PhotonMaxAlpha{"Pi0PhotonMaxAlpha", 0.95, "Max photon alpha absolute value (AP plot)"};
  o2::framework::Configurable<float> Pi0PhotonMinV0cospa{"Pi0PhotonMinV0cospa", 0.80, "Min V0 CosPA"};
  o2::framework::Configurable<float> Pi0PhotonMaxMass{"Pi0PhotonMaxMass", 0.10, "Max photon mass (GeV/c^{2})"};
};

struct axisConfigurables : o2::framework::ConfigurableGroup {
  std::string prefix = "axisOpts";

  // base properties
  o2::framework::ConfigurableAxis axisPt{"axisPt", {o2::framework::VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "pt axis for analysis"};
  o2::framework::ConfigurableAxis axisCentrality{"axisCentrality", {o2::framework::VARIABLE_WIDTH, 0.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.0f, 110.0f}, "Centrality"};

  // Invariant Mass
  o2::framework::ConfigurableAxis axisSigmaMass{"axisSigmaMass", {500, 1.10f, 1.30f}, "M_{#Sigma^{0}} (GeV/c^{2})"};
  o2::framework::ConfigurableAxis axisLambdaMass{"axisLambdaMass", {200, 1.05f, 1.151f}, "M_{#Lambda} (GeV/c^{2})"};
  o2::framework::ConfigurableAxis axisPhotonMass{"axisPhotonMass", {200, -0.1f, 0.5f}, "M_{#Gamma}"};
  o2::framework::ConfigurableAxis axisPi0Mass{"axisPi0Mass", {200, 0.08f, 0.18f}, "M_{#Pi^{0}}"};

  // AP plot axes
  o2::framework::ConfigurableAxis axisAPAlpha{"axisAPAlpha", {220, -1.1f, 1.1f}, "V0 AP alpha"};
  o2::framework::ConfigurableAxis axisAPQt{"axisAPQt", {220, 0.0f, 0.5f}, "V0 AP alpha"};

  // Track quality axes
  o2::framework::ConfigurableAxis axisTPCrows{"axisTPCrows", {160, 0.0f, 160.0f}, "N TPC rows"};

  // topological variable QA axes
  o2::framework::ConfigurableAxis axisDCAtoPV{"axisDCAtoPV", {500, 0.0f, 50.0f}, "DCA (cm)"};
  o2::framework::ConfigurableAxis axisXY{"axisXY", {120, -120.0f, 120.0f}, "XY axis"};
  o2::framework::ConfigurableAxis axisDCAdau{"axisDCAdau", {50, 0.0f, 5.0f}, "DCA (cm)"};
  o2::framework::ConfigurableAxis axisRadius{"axisRadius", {240, 0.0f, 120.0f}, "V0 radius (cm)"};
  o2::framework::ConfigurableAxis axisPA{"axisPA", {100, 0.0f, 1}, "Pointing angle"};
  o2::framework::ConfigurableAxis axisRapidity{"axisRapidity", {100, -2.0f, 2.0f}, "Rapidity"};
  o2::framework::ConfigurableAxis axisCandSel{"axisCandSel", {7, 0.5f, +7.5f}, "Candidate Selection"};
  o2::framework::ConfigurableAxis axisNch{"axisNch", {300, 0.0f, 3000.0f}, "N_{ch}"};
  o2::framework::ConfigurableAxis axisIRBinning{"axisIRBinning", {150, 0, 1500}, "Binning for the interaction rate (kHz)"};
};


class Sigma0BuilderModule
{
 public:
  Sigma0BuilderModule()
  {
    // constructor
  }

  template <typename TCollision, typename THistoRegistry, typename TCCDB, typename TRateFetcher>
  bool IsEventAccepted(TCollision collision, THistoRegistry& histos, TCCDB ccdb, TRateFetcher rateFetcher)
  // check whether the collision passes our collision selections
  {
    if (evselOpts.requireSel8 && !collision.sel8()) {
      return false;
    }    
    if (evselOpts.requireTriggerTVX && !collision.selection_bit(aod::evsel::kIsTriggerTVX)) {
      return false;
    }    
    if (evselOpts.rejectITSROFBorder && !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      return false;
    }    
    if (evselOpts.rejectTFBorder && !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
      return false;
    }    
    if (std::abs(collision.posZ()) > evselOpts.maxZVtxPosition) {
      return false;
    }  
    if (evselOpts.requireIsVertexITSTPC && !collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
      return false;
    }    
    if (evselOpts.requireIsGoodZvtxFT0VsPV && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }    
    if (evselOpts.requireIsVertexTOFmatched && !collision.selection_bit(o2::aod::evsel::kIsVertexTOFmatched)) {
      return false;
    }    
    if (evselOpts.requireIsVertexTRDmatched && !collision.selection_bit(o2::aod::evsel::kIsVertexTRDmatched)) {
      return false;
    }    
    if (evselOpts.rejectSameBunchPileup && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return false;
    }    
    if (evselOpts.requireNoCollInTimeRangeStd && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      return false;
    }    
    if (evselOpts.requireNoCollInTimeRangeStrict && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStrict)) {
      return false;
    }    
    if (evselOpts.requireNoCollInTimeRangeNarrow && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeNarrow)) {
      return false;
    }    
    if (evselOpts.requireNoCollInROFStd && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) {
      return false;
    }    
    if (evselOpts.requireNoCollInROFStrict && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStrict)) {
      return false;
    }    
    if (evselOpts.doPPAnalysis) { // we are in pp
      if (evselOpts.requireINEL0 && collision.multNTracksPVeta1() < 1) {
        return false;
      }      
      if (evselOpts.requireINEL1 && collision.multNTracksPVeta1() < 2) {
        return false;
      }      
    } else { // we are in Pb-Pb
      float collisionOccupancy = evselOpts.useFT0CbasedOccupancy ? collision.ft0cOccupancyInTimeRange() : collision.trackOccupancyInTimeRange();
      if (evselOpts.minOccupancy >= 0 && collisionOccupancy < evselOpts.minOccupancy) {
        return false;
      }      
      if (evselOpts.maxOccupancy >= 0 && collisionOccupancy > evselOpts.maxOccupancy) {
        return false;
      }      
    }
    // Fetch interaction rate only if required (in order to limit ccdb calls)
    double interactionRate = (evselOpts.minIR >= 0 || evselOpts.maxIR >= 0) ? rateFetcher.fetch(ccdb.service, collision.timestamp(), collision.runNumber(), evselOpts.irSource, evselOpts.fIRCrashOnNull) * 1.e-3 : -1;
    if (evselOpts.minIR >= 0 && interactionRate < evselOpts.minIR) {
      return false;
    }    
    if (evselOpts.maxIR >= 0 && interactionRate > evselOpts.maxIR) {
      return false;
    }    
    float centrality = evselOpts.doPPAnalysis ? collision.centFT0M() : collision.centFT0C();

    histos.template get<TH1>(HIST("hEventCentrality"))->Fill(centrality);
    histos.template get<TH1>(HIST("hInteractionRate"))->Fill(interactionRate);
    histos.template get<TH2>(HIST("hCentralityVsInteractionRate"))->Fill(centrality, interactionRate);

    return true;
  }

  template <typename TV0Object, typename TCollision, typename THistoRegistry>
  void runPi0QA(TV0Object const& gamma1, TV0Object const& gamma2, TCollision collision, THistoRegistry& histos)
  {
    // Check if both V0s are made of the same tracks
    if (gamma1.posTrackExtraId() == gamma2.posTrackExtraId() ||
        gamma1.negTrackExtraId() == gamma2.negTrackExtraId()) {
      return;
    }

    // Calculate pi0 properties
    std::array<float, 3> pVecGamma1{gamma1.px(), gamma1.py(), gamma1.pz()};
    std::array<float, 3> pVecGamma2{gamma2.px(), gamma2.py(), gamma2.pz()};
    std::array arrpi0{pVecGamma1, pVecGamma2};
    float pi0Mass = RecoDecay::m(arrpi0, std::array{o2::constants::physics::MassPhoton, o2::constants::physics::MassPhoton});
    float pi0Pt = RecoDecay::pt(std::array{gamma1.px() + gamma2.px(), gamma1.py() + gamma2.py()});
    float pi0Y = RecoDecay::y(std::array{gamma1.px() + gamma2.px(), gamma1.py() + gamma2.py(), gamma1.pz() + gamma2.pz()}, o2::constants::physics::MassPi0);
    float centrality = evselOpts.doPPAnalysis ? collision.centFT0M() : collision.centFT0C();

    histos.template get<TH3>(HIST("Pi0QA/h3dMassPi0BeforeSel_Candidates"))->Fill(centrality, pi0Pt, pi0Mass);

    // Photon-specific selections
    auto posTrackGamma1 = gamma1.template posTrackExtra_as<soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>>();
    auto negTrackGamma1 = gamma1.template negTrackExtra_as<soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>>();
    auto posTrackGamma2 = gamma2.template posTrackExtra_as<soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>>();
    auto negTrackGamma2 = gamma2.template negTrackExtra_as<soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>>();

    // Gamma1 Selection
    bool passedTPCGamma1 = (TMath::Abs(posTrackGamma1.tpcNSigmaEl()) < pi0selOpts.Pi0PhotonMaxTPCNSigmas) ||
                            (TMath::Abs(negTrackGamma1.tpcNSigmaEl()) < pi0selOpts.Pi0PhotonMaxTPCNSigmas);

    if (TMath::Abs(gamma1.mGamma()) > pi0selOpts.Pi0PhotonMaxMass ||
        gamma1.qtarm() >= pi0selOpts.Pi0PhotonMaxQt ||
        TMath::Abs(gamma1.alpha()) >= pi0selOpts.Pi0PhotonMaxAlpha ||
        TMath::Abs(gamma1.dcapostopv()) < pi0selOpts.Pi0PhotonMinDCADauToPv ||
        TMath::Abs(gamma1.dcanegtopv()) < pi0selOpts.Pi0PhotonMinDCADauToPv ||
        TMath::Abs(gamma1.dcaV0daughters()) > pi0selOpts.Pi0PhotonMaxDCAV0Dau ||
        TMath::Abs(gamma1.negativeeta()) >= pi0selOpts.Pi0PhotonMaxEta ||
        TMath::Abs(gamma1.positiveeta()) >= pi0selOpts.Pi0PhotonMaxEta ||
        gamma1.v0cosPA() <= pi0selOpts.Pi0PhotonMinV0cospa ||
        gamma1.v0radius() <= pi0selOpts.Pi0PhotonMinRadius ||
        gamma1.v0radius() >= pi0selOpts.Pi0PhotonMaxRadius ||
        posTrackGamma1.tpcCrossedRows() < pi0selOpts.Pi0PhotonMinTPCCrossedRows ||
        negTrackGamma1.tpcCrossedRows() < pi0selOpts.Pi0PhotonMinTPCCrossedRows ||
        !passedTPCGamma1) {
      return;
    }

    // Gamma2 Selection
    bool passedTPCGamma2 = (TMath::Abs(posTrackGamma2.tpcNSigmaEl()) < pi0selOpts.Pi0PhotonMaxTPCNSigmas) ||
                            (TMath::Abs(negTrackGamma2.tpcNSigmaEl()) < pi0selOpts.Pi0PhotonMaxTPCNSigmas);

    if (TMath::Abs(gamma2.mGamma()) > pi0selOpts.Pi0PhotonMaxMass ||
        gamma2.qtarm() >= pi0selOpts.Pi0PhotonMaxQt ||
        TMath::Abs(gamma2.alpha()) >= pi0selOpts.Pi0PhotonMaxAlpha ||
        TMath::Abs(gamma2.dcapostopv()) < pi0selOpts.Pi0PhotonMinDCADauToPv ||
        TMath::Abs(gamma2.dcanegtopv()) < pi0selOpts.Pi0PhotonMinDCADauToPv ||
        TMath::Abs(gamma2.dcaV0daughters()) > pi0selOpts.Pi0PhotonMaxDCAV0Dau ||
        TMath::Abs(gamma2.negativeeta()) >= pi0selOpts.Pi0PhotonMaxEta ||
        TMath::Abs(gamma2.positiveeta()) >= pi0selOpts.Pi0PhotonMaxEta ||
        gamma2.v0cosPA() <= pi0selOpts.Pi0PhotonMinV0cospa ||
        gamma2.v0radius() <= pi0selOpts.Pi0PhotonMinRadius ||
        gamma2.v0radius() >= pi0selOpts.Pi0PhotonMaxRadius ||
        posTrackGamma2.tpcCrossedRows() < pi0selOpts.Pi0PhotonMinTPCCrossedRows ||
        negTrackGamma2.tpcCrossedRows() < pi0selOpts.Pi0PhotonMinTPCCrossedRows ||
        !passedTPCGamma2) {
      return;
    }

    // Pi0-specific selections:
    if (TMath::Abs(pi0Y) > 0.5) {
      return;
    }

    // Fill histograms    
    histos.template get<TH3>(HIST("Pi0QA/h3dMassPi0AfterSel_Candidates"))->Fill(centrality, pi0Pt, pi0Mass);
  }

  // Process photon candidate
  template <typename TV0Object, typename TCollision, typename THistoRegistry>
  bool processPhotonCandidate(TV0Object const& gamma, TCollision collision, THistoRegistry& histos)
  {
    if (gamma.v0Type() == 0)
      return false;

    if (photonselOpts.useMLScores) {
      // Gamma selection:
      if (gamma.gammaBDTScore() <= photonselOpts.Gamma_MLThreshold)
        return false;

    } else {
      // Standard selection
      // Gamma basic selection criteria:      
      if ((gamma.mGamma() < 0) || (gamma.mGamma() > photonselOpts.PhotonMaxMass))
        return false;
      float PhotonY = RecoDecay::y(std::array{gamma.px(), gamma.py(), gamma.pz()}, o2::constants::physics::MassGamma);      
      if ((TMath::Abs(PhotonY) > photonselOpts.PhotonRapidity) || (TMath::Abs(gamma.negativeeta()) > photonselOpts.PhotonMaxDauPseudoRap) || (TMath::Abs(gamma.positiveeta()) > photonselOpts.PhotonMaxDauPseudoRap))
        return false;      
      if ((TMath::Abs(gamma.dcapostopv()) < photonselOpts.PhotonMinDCAToPv) || (TMath::Abs(gamma.dcanegtopv()) < photonselOpts.PhotonMinDCAToPv))
        return false;      
      if (TMath::Abs(gamma.dcaV0daughters()) > photonselOpts.PhotonMaxDCAV0Dau)
        return false;      
      if ((gamma.v0radius() < photonselOpts.PhotonMinRadius) || (gamma.v0radius() > photonselOpts.PhotonMaxRadius))
        return false;      
    }
    float centrality = evselOpts.doPPAnalysis ? collision.centFT0M() : collision.centFT0C();
    histos.template get<TH3>(HIST("PhotonSel/h3dPhotonMass"))->Fill(centrality, gamma.pt(), gamma.mGamma());

    return true;
  }

  // Process photon candidate
  template <typename TV0Object, typename TCollision, typename THistoRegistry>
  bool processLambdaCandidate(TV0Object const& lambda, TCollision collision, THistoRegistry& histos)
  {
    if (lambda.v0Type() != 1)
      return false;

    if (lambdaselOpts.useMLScores) {
      if ((lambda.lambdaBDTScore() <= lambdaselOpts.Lambda_MLThreshold) && (lambda.antiLambdaBDTScore() <= lambdaselOpts.AntiLambda_MLThreshold))
        return false;

    } else {
      // Lambda basic selection criteria:      
      if ((TMath::Abs(lambda.mLambda() - o2::constants::physics::MassLambda0) > lambdaselOpts.LambdaWindow) && (TMath::Abs(lambda.mAntiLambda() - o2::constants::physics::MassLambda0) > lambdaselOpts.LambdaWindow))
        return false;      
      if ((TMath::Abs(lambda.yLambda()) > lambdaselOpts.LambdaRapidity) || (TMath::Abs(lambda.negativeeta()) > lambdaselOpts.LambdaDauPseudoRap) || (TMath::Abs(lambda.positiveeta()) > lambdaselOpts.LambdaDauPseudoRap))
        return false;      
      if ((TMath::Abs(lambda.dcapostopv()) < lambdaselOpts.LambdaMinDCAPosToPv) || (TMath::Abs(lambda.dcanegtopv()) < lambdaselOpts.LambdaMinDCANegToPv))
        return false;      
      if ((lambda.v0radius() < lambdaselOpts.LambdaMinv0radius) || (lambda.v0radius() > lambdaselOpts.LambdaMaxv0radius))
        return false;      
      if (TMath::Abs(lambda.dcaV0daughters()) > lambdaselOpts.LambdaMaxDCAV0Dau)
        return false;      
    }

    float centrality = evselOpts.doPPAnalysis ? collision.centFT0M() : collision.centFT0C();
    histos.template get<TH3>(HIST("LambdaSel/h3dLambdaMass"))->Fill(centrality, lambda.pt(), lambda.mLambda());
    histos.template get<TH3>(HIST("LambdaSel/h3dALambdaMass"))->Fill(centrality, lambda.pt(), lambda.mAntiLambda());

    return true;
  }
  ///////////
  // Process sigma candidate and store properties in object
  template <typename TV0Object, typename TCollision, typename THistoRegistry>
  bool buildSigma0(TV0Object const& lambda, TV0Object const& gamma, TCollision collision, THistoRegistry& histos)
  {
    // Checking if both V0s are made of the very same tracks
    if (gamma.posTrackExtraId() == lambda.posTrackExtraId() ||
        gamma.negTrackExtraId() == lambda.negTrackExtraId()) {
      return false;
    }

    // Sigma0 candidate properties
    std::array<float, 3> pVecPhotons{gamma.px(), gamma.py(), gamma.pz()};
    std::array<float, 3> pVecLambda{lambda.px(), lambda.py(), lambda.pz()};
    auto arrMom = std::array{pVecPhotons, pVecLambda};
    float sigmaMass = RecoDecay::m(arrMom, std::array{o2::constants::physics::MassPhoton, o2::constants::physics::MassLambda0});
    float sigmaY = RecoDecay::y(std::array{gamma.px() + lambda.px(), gamma.py() + lambda.py(), gamma.pz() + lambda.pz()}, o2::constants::physics::MassSigma0);
    float SigmapT = RecoDecay::pt(std::array{gamma.px() + lambda.px(), gamma.py() + lambda.py()});
    float centrality = evselOpts.doPPAnalysis ? collision.centFT0M() : collision.centFT0C();

    // Before any selection
    histos.template get<TH3>(HIST("SigmaSel/h3dMassSigma0BeforeSel"))->Fill(centrality, SigmapT, sigmaMass);
    histos.template get<TH1>(HIST("SigmaSel/hSigmaMass"))->Fill(sigmaMass);
    histos.template get<TH1>(HIST("SigmaSel/hSigmaMassWindow"))->Fill(sigmaMass - 1.192642);

    if (TMath::Abs(sigmaMass - 1.192642) > sigma0selOpts.Sigma0Window)
      return false;
    
    histos.template get<TH1>(HIST("SigmaSel/hSigmaY"))->Fill(sigmaY);

    if (TMath::Abs(sigmaY) > sigma0selOpts.SigmaMaxRap)
      return false;

    histos.template get<TH1>(HIST("SigmaSel/hSigmaMassSelected"))->Fill(sigmaMass);
    histos.template get<TH3>(HIST("SigmaSel/h3dMassSigma0AfterSel"))->Fill(centrality, SigmapT, sigmaMass);

    return true;
  }

  // declaration of structs here
  // (N.B.: will be invisible to the outside, create your own copies)
  o2::pwglf::sigma0::evselConfigurables evselOpts;
  o2::pwglf::sigma0::lambdaselConfigurables lambdaselOpts;
  o2::pwglf::sigma0::photonselConfigurables photonselOpts;
  o2::pwglf::sigma0::sigma0selConfigurables sigma0selOpts;
  o2::pwglf::sigma0::pi0selConfigurables pi0selOpts;
  o2::pwglf::sigma0::axisConfigurables axisOpts;

  template <typename THistoRegistry, typename TEvSelOpt, typename TLambdaSelOpt, typename TPhotonSelOpt, typename TSigma0SelOpt, typename TPi0Opt, typename TAxisOpt>
  void init(THistoRegistry& histos, 
            TEvSelOpt const& external_evselopts,
            TLambdaSelOpt const& external_lambdaselopts,
            TPhotonSelOpt const& external_photonselopts,
            TSigma0SelOpt const& external_sigma0selopts,
            TPi0Opt const& external_pi0selopts,
            TAxisOpt const& external_axisopts)
  {
    // read in configurations from the task where it's used
    evselOpts = external_evselopts;
    lambdaselOpts = external_lambdaselopts;
    photonselOpts = external_photonselopts;
    sigma0selOpts = external_sigma0selopts;
    pi0selOpts = external_pi0selopts;
    axisOpts = external_axisopts;

    histos.add("hEventCentrality", "hEventCentrality", framework::kTH1D, {axisOpts.axisCentrality});    
    histos.add("hInteractionRate", "hInteractionRate", framework::kTH1F, {axisOpts.axisIRBinning});
    histos.add("hCentralityVsInteractionRate", "hCentralityVsInteractionRate", o2::framework::kTH2F, {axisOpts.axisCentrality, axisOpts.axisIRBinning});

    // For selection:
    histos.add("PhotonSel/h3dPhotonMass", "h3dPhotonMass", framework::kTH3D, {axisOpts.axisCentrality, axisOpts.axisPt, axisOpts.axisPhotonMass});
    histos.add("LambdaSel/h3dLambdaMass", "h3dLambdaMass", framework::kTH3D, {axisOpts.axisCentrality, axisOpts.axisPt, axisOpts.axisLambdaMass});
    histos.add("LambdaSel/h3dALambdaMass", "h3dALambdaMass", framework::kTH3D, {axisOpts.axisCentrality, axisOpts.axisPt, axisOpts.axisLambdaMass});

    histos.add("SigmaSel/h3dMassSigma0BeforeSel", "h3dMassSigma0BeforeSel", framework::kTH3F, {axisOpts.axisCentrality, axisOpts.axisPt, axisOpts.axisSigmaMass});
    histos.add("SigmaSel/hSigmaMass", "hSigmaMass", framework::kTH1F, {axisOpts.axisSigmaMass});
    histos.add("SigmaSel/hSigmaMassWindow", "hSigmaMassWindow", framework::kTH1F, {{200, -0.09f, 0.11f}});
    histos.add("SigmaSel/hSigmaY", "hSigmaY", framework::kTH1F, {axisOpts.axisRapidity});
    histos.add("SigmaSel/hSigmaMassSelected", "hSigmaMassSelected", framework::kTH1F, {axisOpts.axisSigmaMass});
    histos.add("SigmaSel/h3dMassSigma0AfterSel", "h3dMassSigma0AfterSel", framework::kTH3D, {axisOpts.axisCentrality, axisOpts.axisPt, axisOpts.axisSigmaMass});
  
    // For Pi0 QA
    histos.add("Pi0QA/h3dMassPi0BeforeSel_Candidates", "h3dMassPi0BeforeSel_Candidates", framework::kTH3D, {axisOpts.axisCentrality, axisOpts.axisPt, axisOpts.axisPi0Mass});
    histos.add("Pi0QA/h3dMassPi0AfterSel_Candidates", "h3dMassPi0AfterSel_Candidates", framework::kTH3D, {axisOpts.axisCentrality, axisOpts.axisPt, axisOpts.axisPi0Mass});
  }

  // ______________________________________________________
  // Real data processing - no MC subscription
  template <typename TCollision, typename TV0s, typename THistoRegistry, typename TSlicecache, typename TCCDB, typename TRateFetcher>
  void process(TCollision const& collisions, TV0s const& fullV0s, THistoRegistry& histos, TSlicecache& cache, TCCDB const& ccdb, TRateFetcher& rateFetcher)
  {
    uint64_t CollIDBuffer = 0;
    for (const auto& coll : collisions) {
      if (!IsEventAccepted(coll, histos, ccdb, rateFetcher))
        continue;

      const uint64_t collIdx = coll.globalIndex();      
      if (collIdx < CollIDBuffer) 
        LOGF(fatal, "Collision table unsorted! Previous index: %i, current index: %i", CollIDBuffer, collIdx); 

      CollIDBuffer = collIdx;

      //_______________________________________________
      // V0s loop
      std::vector<int> bestGammasArray;
      std::vector<int> bestLambdasArray;

      auto V0s = fullV0s.sliceByCached(o2::aod::v0data::straCollisionId, collIdx, cache);      
      for (auto& v0 : V0s) {
        if (processPhotonCandidate(v0, coll, histos))          // selecting photons
          bestGammasArray.push_back(v0.globalIndex()); // Save indices of best gamma candidates

        if (processLambdaCandidate(v0, coll, histos))           // selecting lambdas
          bestLambdasArray.push_back(v0.globalIndex()); // Save indices of best lambda candidates
      }
    
       
      //_______________________________________________
      // Pi0 optional loop
      if (pi0selOpts.doPi0QA) {
        for (size_t i = 0; i < bestGammasArray.size(); ++i) {
          auto gamma1 = fullV0s.rawIteratorAt(bestGammasArray[i]);
          for (size_t j = i + 1; j < bestGammasArray.size(); ++j) {
            auto gamma2 = fullV0s.rawIteratorAt(bestGammasArray[j]);
            runPi0QA(gamma1, gamma2, coll, histos);
          }
        }
      }

      //_______________________________________________
      // Sigma0 nested loop
      for (size_t i = 0; i < bestGammasArray.size(); ++i) {
        auto gamma = fullV0s.rawIteratorAt(bestGammasArray[i]);

        for (size_t j = 0; j < bestLambdasArray.size(); ++j) {
          auto lambda = fullV0s.rawIteratorAt(bestLambdasArray[j]);

          // Building sigma0 candidate
          if (!buildSigma0(lambda, gamma, coll, histos))
            continue;          
        }
      }
    }
  } // end process 
}; // end Sigma0BuilderModule

} // namespace sigma0
} // namespace pwglf
} // namespace o2

#endif // PWGLF_UTILS_SIGMA0BUILDERHELPER_H_
