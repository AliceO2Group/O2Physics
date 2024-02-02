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
//
// V0 analysis task
// ================
//
// This code loops over a V0Cores table and produces some
// standard analysis output. It is meant to be run over
// derived data.
//
//    Comments, questions, complaints, suggestions?
//    Please write to:
//    romain.schotter@cern.ch
//    david.dobrigkeit.chinellato@cern.ch
//
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Common/Core/trackUtilities.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"

#include <TFile.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TLorentzVector.h>
#include <Math/Vector4D.h>
#include <TPDGCode.h>
#include <TDatabasePDG.h>
#include <cmath>
#include <array>
#include <cstdlib>
#include "Framework/ASoAHelpers.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

using dauTracks = soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>;
using v0Candidates = soa::Join<aod::V0CollRefs, aod::V0Cores, aod::V0Extras>;
using v0MCCandidates = soa::Join<aod::V0CollRefs, aod::V0Cores, aod::V0MCCores, aod::V0Extras>;

// simple checkers
#define bitset(var, nbit) ((var) |= (1 << (nbit)))
#define bitcheck(var, nbit) ((var) & (1 << (nbit)))

struct derivedlambdakzeroanalysis {
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // master analysis switches
  Configurable<bool> analyseK0Short{"analyseK0Short", true, "process K0Short-like candidates"};
  Configurable<bool> analyseLambda{"analyseLambda", true, "process Lambda-like candidates"};
  Configurable<bool> analyseAntiLambda{"analyseAntiLambda", true, "process AntiLambda-like candidates"};

  // Selection criteria: acceptance
  Configurable<float> rapidityCut{"rapidityCut", 0.5, "rapidity"};
  Configurable<float> daughterEtaCut{"daughterEtaCut", 0.8, "max eta for daughters"};

  // Standard 5 topological criteria
  Configurable<float> v0cospa{"v0cospa", 0.97, "min V0 CosPA"};
  Configurable<float> dcav0dau{"dcav0dau", 1.0, "max DCA V0 Daughters (cm)"};
  Configurable<float> dcanegtopv{"dcanegtopv", .05, "min DCA Neg To PV (cm)"};
  Configurable<float> dcapostopv{"dcapostopv", .05, "min DCA Pos To PV (cm)"};
  Configurable<float> v0radius{"v0radius", 1.2, "minimum V0 radius (cm)"};

  // Additional selection on the AP plot (exclusive for K0Short)
  // original equation: lArmPt*5>TMath::Abs(lArmAlpha)
  Configurable<float> armPodCut{"armPodCut", 5.0f, "pT * (cut) > |alpha|, AP cut. Negative: no cut"};

  // Track quality
  Configurable<int> minTPCrows{"minTPCrows", 70, "minimum TPC crossed rows"};
  Configurable<bool> requirePosITSonly{"requirePosITSonly", false, "require that positive track is ITSonly (overrides TPC quality)"};
  Configurable<bool> requireNegITSonly{"requireNegITSonly", false, "require that negative track is ITSonly (overrides TPC quality)"};

  // PID (TPC)
  Configurable<float> TpcPidNsigmaCut{"TpcPidNsigmaCut", 5, "TpcPidNsigmaCut"};

  Configurable<bool> doQA{"doQA", true, "do topological variable QA histograms"};
  Configurable<float> qaMinPt{"qaMinPt", 0.0f, "minimum pT for QA plots"};
  Configurable<float> qaMaxPt{"qaMaxPt", 0.0f, "maximum pT for QA plots"};

  // for MC
  Configurable<bool> doMCAssociation{"doMCAssociation", true, "if MC, do MC association"};

  static constexpr float defaultLifetimeCuts[1][2] = {{30., 20.}};
  Configurable<LabeledArray<float>> lifetimecut{"lifetimecut", {defaultLifetimeCuts[0], 2, {"lifetimecutLambda", "lifetimecutK0S"}}, "lifetimecut"};

  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "pt axis for analysis"};
  ConfigurableAxis axisPtCoarse{"axisPtCoarse", {VARIABLE_WIDTH, 0.0f, 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 7.0f, 10.0f, 15.0f}, "pt axis for QA"};
  ConfigurableAxis axisK0Mass{"axisK0Mass", {200, 0.4f, 0.6f}, ""};
  ConfigurableAxis axisLambdaMass{"axisLambdaMass", {200, 1.101f, 1.131f}, ""};
  ConfigurableAxis axisCentrality{"axisCentrality", {VARIABLE_WIDTH, 0.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f}, "Centrality"};

  // topological variable QA axes
  ConfigurableAxis axisDCAtoPV{"axisDCAtoPV", {20, 0.0f, 1.0f}, "DCA (cm)"};
  ConfigurableAxis axisDCAdau{"axisDCAdau", {20, 0.0f, 2.0f}, "DCA (cm)"};
  ConfigurableAxis axisPointingAngle{"axisPointingAngle", {20, 0.0f, 2.0f}, "pointing angle (rad)"};
  ConfigurableAxis axisV0Radius{"axisV0Radius", {20, 0.0f, 60.0f}, "V0 2D radius (cm)"};

  // AP plot axes
  ConfigurableAxis axisAPAlpha{"axisAPAlpha", {220, -1.1f, 1.1f}, "V0 AP alpha"};
  ConfigurableAxis axisAPQt{"axisAPQt", {220, 0.0f, 0.5f}, "V0 AP alpha"};

  enum species { spK0Short = 0,
                 spLambda,
                 spAntiLambda };

  enum selection { selCosPA = 0,
                   selRadius,
                   selDCANegToPV,
                   selDCAPosToPV,
                   selDCAV0Dau,
                   selK0ShortRapidity,
                   selLambdaRapidity,
                   selK0ShortTPC,
                   selLambdaTPC,
                   selAntiLambdaTPC,
                   selK0ShortCTau,
                   selLambdaCTau,
                   selK0ShortArmenteros,
                   selPosGoodTPCTrack,
                   selNegGoodTPCTrack,
                   selPosItsOnly,
                   selNegItsOnly,
                   selConsiderK0Short,   // for mc tagging
                   selConsiderLambda,    // for mc tagging
                   selConsiderAntiLambda // for mc tagging
  };

  uint32_t maskTopological;
  uint32_t maskTopoNoV0Radius;
  uint32_t maskTopoNoDCANegToPV;
  uint32_t maskTopoNoDCAPosToPV;
  uint32_t maskTopoNoCosPA;
  uint32_t maskTopoNoDCAV0Dau;
  uint32_t maskTrackTypes;

  uint32_t maskK0ShortSpecific;
  uint32_t maskLambdaSpecific;
  uint32_t maskAntiLambdaSpecific;

  uint32_t maskSelectionK0Short;
  uint32_t maskSelectionLambda;
  uint32_t maskSelectionAntiLambda;

  void init(InitContext const&)
  {
    // initialise bit masks
    maskTopological = (1 << selCosPA) | (1 << selRadius) | (1 << selDCANegToPV) | (1 << selDCAPosToPV) | (1 << selDCAV0Dau);
    maskTopoNoV0Radius = (1 << selCosPA) | (1 << selDCANegToPV) | (1 << selDCAPosToPV) | (1 << selDCAV0Dau);
    maskTopoNoDCANegToPV = (1 << selCosPA) | (1 << selRadius) | (1 << selDCAPosToPV) | (1 << selDCAV0Dau);
    maskTopoNoDCAPosToPV = (1 << selCosPA) | (1 << selRadius) | (1 << selDCANegToPV) | (1 << selDCAV0Dau);
    maskTopoNoCosPA = (1 << selRadius) | (1 << selDCANegToPV) | (1 << selDCAPosToPV) | (1 << selDCAV0Dau);
    maskTopoNoDCAV0Dau = (1 << selCosPA) | (1 << selRadius) | (1 << selDCANegToPV) | (1 << selDCAPosToPV);

    maskTrackTypes = 0;
    if (requirePosITSonly) {
      maskTrackTypes = (1 << selPosItsOnly);
    } else {
      maskTrackTypes = (1 << selPosGoodTPCTrack);
    }
    if (requireNegITSonly) {
      maskTrackTypes = (1 << selNegItsOnly);
    } else {
      maskTrackTypes = (1 << selNegGoodTPCTrack);
    }

    maskK0ShortSpecific = (1 << selK0ShortRapidity) | (1 << selK0ShortTPC) | (1 << selK0ShortCTau) | (1 << selK0ShortArmenteros) | (1 << selConsiderK0Short);
    maskLambdaSpecific = (1 << selLambdaRapidity) | (1 << selLambdaTPC) | (1 << selLambdaCTau) | (1 << selConsiderLambda);
    maskAntiLambdaSpecific = (1 << selLambdaRapidity) | (1 << selAntiLambdaTPC) | (1 << selLambdaCTau) | (1 << selConsiderAntiLambda);

    maskSelectionK0Short = maskTopological | maskTrackTypes | maskK0ShortSpecific;
    maskSelectionLambda = maskTopological | maskTrackTypes | maskLambdaSpecific;
    maskSelectionAntiLambda = maskTopological | maskTrackTypes | maskAntiLambdaSpecific;

    // Event Counters
    histos.add("hEventSelection", "hEventSelection", kTH1F, {{3, -0.5f, +2.5f}});
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(1, "All collisions");
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(2, "sel8 cut");
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(3, "posZ cut");

    histos.add("hEventCentrality", "hEventCentrality", kTH1F, {{100, 0.0f, +100.0f}});

    // histograms versus mass
    if (analyseK0Short)
      histos.add("h3dMassK0Short", "h3dMassK0Short", kTH3F, {axisCentrality, axisPt, axisK0Mass});
    if (analyseLambda)
      histos.add("h3dMassLambda", "h3dMassLambda", kTH3F, {axisCentrality, axisPt, axisLambdaMass});
    if (analyseAntiLambda)
      histos.add("h3dMassAntiLambda", "h3dMassAntiLambda", kTH3F, {axisCentrality, axisPt, axisLambdaMass});

    // demo // fast
    histos.add("hMassK0Short", "hMassK0Short", kTH1F, {axisK0Mass});

    // QA histograms if requested
    if (doQA) {
      // initialize for K0short...
      if (analyseK0Short) {
        histos.add("K0Short/h4dPosDCAToPV", "h4dPosDCAToPV", kTHnF, {axisCentrality, axisPtCoarse, axisK0Mass, axisDCAtoPV});
        histos.add("K0Short/h4dNegDCAToPV", "h4dNegDCAToPV", kTHnF, {axisCentrality, axisPtCoarse, axisK0Mass, axisDCAtoPV});
        histos.add("K0Short/h4dDCADaughters", "h4dDCADaughters", kTHnF, {axisCentrality, axisPtCoarse, axisK0Mass, axisDCAdau});
        histos.add("K0Short/h4dPointingAngle", "h4dPointingAngle", kTHnF, {axisCentrality, axisPtCoarse, axisK0Mass, axisPointingAngle});
        histos.add("K0Short/h4dV0Radius", "h4dV0Radius", kTHnF, {axisCentrality, axisPtCoarse, axisK0Mass, axisV0Radius});
      }
      if (analyseLambda) {
        histos.add("Lambda/h4dPosDCAToPV", "h4dPosDCAToPV", kTHnF, {axisCentrality, axisPtCoarse, axisLambdaMass, axisDCAtoPV});
        histos.add("Lambda/h4dNegDCAToPV", "h4dNegDCAToPV", kTHnF, {axisCentrality, axisPtCoarse, axisLambdaMass, axisDCAtoPV});
        histos.add("Lambda/h4dDCADaughters", "h4dDCADaughters", kTHnF, {axisCentrality, axisPtCoarse, axisLambdaMass, axisDCAdau});
        histos.add("Lambda/h4dPointingAngle", "h4dPointingAngle", kTHnF, {axisCentrality, axisPtCoarse, axisLambdaMass, axisPointingAngle});
        histos.add("Lambda/h4dV0Radius", "h4dV0Radius", kTHnF, {axisCentrality, axisPtCoarse, axisLambdaMass, axisV0Radius});
      }
      if (analyseAntiLambda) {
        histos.add("AntiLambda/h4dPosDCAToPV", "h4dPosDCAToPV", kTHnF, {axisCentrality, axisPtCoarse, axisLambdaMass, axisDCAtoPV});
        histos.add("AntiLambda/h4dNegDCAToPV", "h4dNegDCAToPV", kTHnF, {axisCentrality, axisPtCoarse, axisLambdaMass, axisDCAtoPV});
        histos.add("AntiLambda/h4dDCADaughters", "h4dDCADaughters", kTHnF, {axisCentrality, axisPtCoarse, axisLambdaMass, axisDCAdau});
        histos.add("AntiLambda/h4dPointingAngle", "h4dPointingAngle", kTHnF, {axisCentrality, axisPtCoarse, axisLambdaMass, axisPointingAngle});
        histos.add("AntiLambda/h4dV0Radius", "h4dV0Radius", kTHnF, {axisCentrality, axisPtCoarse, axisLambdaMass, axisV0Radius});
      }

      // Check if doing the right thing in AP space please
      histos.add("GeneralQA/h2dArmenterosAll", "h2dArmenterosAll", kTH2F, {axisAPAlpha, axisAPQt});
      histos.add("GeneralQA/h2dArmenterosSelected", "h2dArmenterosSelected", kTH2F, {axisAPAlpha, axisAPQt});
    }
  }

  template <typename TV0>
  bool compatibleTPC(TV0 v0, int sp)
  {
    float pidPos = TMath::Abs(v0.template posTrackExtra_as<dauTracks>().tpcNSigmaPi());
    float pidNeg = TMath::Abs(v0.template negTrackExtra_as<dauTracks>().tpcNSigmaPi());

    if (sp == spLambda)
      pidPos = TMath::Abs(v0.template posTrackExtra_as<dauTracks>().tpcNSigmaPr());
    if (sp == spAntiLambda)
      pidNeg = TMath::Abs(v0.template negTrackExtra_as<dauTracks>().tpcNSigmaPr());

    if (pidPos < TpcPidNsigmaCut && pidNeg < TpcPidNsigmaCut)
      return true;

    // if not, then not
    return false;
  }

  template <typename TV0, typename TCollision>
  uint32_t computeReconstructionBitmap(TV0 v0, TCollision collision)
  // precalculate this information so that a check is one mask operation, not many
  {
    uint32_t bitMap = 0;
    // Base topological variables
    if (v0.v0radius() > v0radius)
      bitset(bitMap, selRadius);
    if (TMath::Abs(v0.dcapostopv()) > dcapostopv)
      bitset(bitMap, selDCAPosToPV);
    if (TMath::Abs(v0.dcanegtopv()) > dcanegtopv)
      bitset(bitMap, selDCANegToPV);
    if (v0.v0cosPA() > v0cospa)
      bitset(bitMap, selCosPA);
    if (v0.dcaV0daughters() < dcav0dau)
      bitset(bitMap, selDCAV0Dau);

    // rapidity
    if (TMath::Abs(v0.yLambda()) < rapidityCut)
      bitset(bitMap, selLambdaRapidity);
    if (TMath::Abs(v0.yK0Short()) < rapidityCut)
      bitset(bitMap, selK0ShortRapidity);

    // TPC PID
    if (compatibleTPC(v0, spK0Short))
      bitset(bitMap, selK0ShortTPC);
    if (compatibleTPC(v0, spLambda))
      bitset(bitMap, selLambdaTPC);
    if (compatibleTPC(v0, spAntiLambda))
      bitset(bitMap, selAntiLambdaTPC);

    auto posTrackExtra = v0.template posTrackExtra_as<dauTracks>();
    auto negTrackExtra = v0.template negTrackExtra_as<dauTracks>();

    // TPC quality
    if (posTrackExtra.tpcCrossedRows() >= minTPCrows)
      bitset(bitMap, selPosGoodTPCTrack);
    if (negTrackExtra.tpcCrossedRows() >= minTPCrows)
      bitset(bitMap, selNegGoodTPCTrack);

    // ITS only tag
    if (posTrackExtra.tpcCrossedRows() < 1)
      bitset(bitMap, selPosItsOnly);
    if (negTrackExtra.tpcCrossedRows() < 1)
      bitset(bitMap, selNegItsOnly);

    // proper lifetime
    if (v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0 < lifetimecut->get("lifetimecutLambda"))
      bitset(bitMap, selLambdaCTau);
    if (v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassK0Short < lifetimecut->get("lifetimecutK0S"))
      bitset(bitMap, selK0ShortCTau);

    // armenteros
    if (v0.qtarm() * armPodCut > TMath::Abs(v0.alpha()) || armPodCut < 1e-4)
      bitset(bitMap, selK0ShortArmenteros);

    return bitMap;
  }

  template <typename TV0>
  uint32_t computeMCAssociation(TV0 v0)
  // precalculate this information so that a check is one mask operation, not many
  {
    uint32_t bitMap = 0;
    // check for specific particle species

    if (v0.pdgCode() == 310 && v0.pdgCodePositive() == 211 && v0.pdgCodeNegative() == -211 && v0.isPhysicalPrimary()) {
      bitset(bitMap, selConsiderK0Short);
    }
    if (v0.pdgCode() == 3122 && v0.pdgCodePositive() == 2212 && v0.pdgCodeNegative() == -211 && v0.isPhysicalPrimary()) {
      bitset(bitMap, selConsiderLambda);
    }
    if (v0.pdgCode() == -3122 && v0.pdgCodePositive() == 211 && v0.pdgCodeNegative() == -2212 && v0.isPhysicalPrimary()) {
      bitset(bitMap, selConsiderAntiLambda);
    }
    return bitMap;
  }

  bool verifyMask(uint32_t bitmap, uint32_t mask)
  {
    return (bitmap & mask) == mask;
  }

  template <typename TV0, typename TCollision>
  void analyseCandidate(TV0 v0, TCollision collision, uint32_t selMap)
  // precalculate this information so that a check is one mask operation, not many
  {
    // __________________________________________
    // main analysis
    if (verifyMask(selMap, maskSelectionK0Short) && analyseK0Short) {
      histos.fill(HIST("GeneralQA/h2dArmenterosSelected"), v0.alpha(), v0.qtarm()); // cross-check
      histos.fill(HIST("h3dMassK0Short"), collision.centFT0C(), v0.pt(), v0.mK0Short());
      histos.fill(HIST("hMassK0Short"), v0.mK0Short());
    }
    if (verifyMask(selMap, maskSelectionLambda) && analyseLambda) {
      histos.fill(HIST("h3dMassLambda"), collision.centFT0C(), v0.pt(), v0.mLambda());
    }
    if (verifyMask(selMap, maskSelectionAntiLambda) && analyseAntiLambda) {
      histos.fill(HIST("h3dMassAntiLambda"), collision.centFT0C(), v0.pt(), v0.mAntiLambda());
    }

    // __________________________________________
    // do systematics / qa plots
    if (doQA) {
      if (analyseK0Short) {
        if (verifyMask(selMap, maskTopoNoV0Radius | maskK0ShortSpecific))
          histos.fill(HIST("K0Short/h4dV0Radius"), collision.centFT0C(), v0.pt(), v0.mK0Short(), v0.v0radius());
        if (verifyMask(selMap, maskTopoNoDCAPosToPV | maskK0ShortSpecific))
          histos.fill(HIST("K0Short/h4dPosDCAToPV"), collision.centFT0C(), v0.pt(), v0.mK0Short(), TMath::Abs(v0.dcapostopv()));
        if (verifyMask(selMap, maskTopoNoDCANegToPV | maskK0ShortSpecific))
          histos.fill(HIST("K0Short/h4dNegDCAToPV"), collision.centFT0C(), v0.pt(), v0.mK0Short(), TMath::Abs(v0.dcanegtopv()));
        if (verifyMask(selMap, maskTopoNoCosPA | maskK0ShortSpecific))
          histos.fill(HIST("K0Short/h4dPointingAngle"), collision.centFT0C(), v0.pt(), v0.mK0Short(), TMath::ACos(v0.v0cosPA()));
        if (verifyMask(selMap, maskTopoNoDCAV0Dau | maskK0ShortSpecific))
          histos.fill(HIST("K0Short/h4dDCADaughters"), collision.centFT0C(), v0.pt(), v0.mK0Short(), v0.dcaV0daughters());
      }

      if (analyseLambda) {
        if (verifyMask(selMap, maskTopoNoV0Radius | maskLambdaSpecific))
          histos.fill(HIST("Lambda/h4dV0Radius"), collision.centFT0C(), v0.pt(), v0.mLambda(), v0.v0radius());
        if (verifyMask(selMap, maskTopoNoDCAPosToPV | maskLambdaSpecific))
          histos.fill(HIST("Lambda/h4dPosDCAToPV"), collision.centFT0C(), v0.pt(), v0.mLambda(), TMath::Abs(v0.dcapostopv()));
        if (verifyMask(selMap, maskTopoNoDCANegToPV | maskLambdaSpecific))
          histos.fill(HIST("Lambda/h4dNegDCAToPV"), collision.centFT0C(), v0.pt(), v0.mLambda(), TMath::Abs(v0.dcanegtopv()));
        if (verifyMask(selMap, maskTopoNoCosPA | maskLambdaSpecific))
          histos.fill(HIST("Lambda/h4dPointingAngle"), collision.centFT0C(), v0.pt(), v0.mLambda(), TMath::ACos(v0.v0cosPA()));
        if (verifyMask(selMap, maskTopoNoDCAV0Dau | maskLambdaSpecific))
          histos.fill(HIST("Lambda/h4dDCADaughters"), collision.centFT0C(), v0.pt(), v0.mLambda(), v0.dcaV0daughters());
      }
      if (analyseAntiLambda) {
        if (verifyMask(selMap, maskTopoNoV0Radius | maskAntiLambdaSpecific))
          histos.fill(HIST("AntiLambda/h4dV0Radius"), collision.centFT0C(), v0.pt(), v0.mAntiLambda(), v0.v0radius());
        if (verifyMask(selMap, maskTopoNoDCAPosToPV | maskAntiLambdaSpecific))
          histos.fill(HIST("AntiLambda/h4dPosDCAToPV"), collision.centFT0C(), v0.pt(), v0.mAntiLambda(), TMath::Abs(v0.dcapostopv()));
        if (verifyMask(selMap, maskTopoNoDCANegToPV | maskAntiLambdaSpecific))
          histos.fill(HIST("AntiLambda/h4dNegDCAToPV"), collision.centFT0C(), v0.pt(), v0.mAntiLambda(), TMath::Abs(v0.dcanegtopv()));
        if (verifyMask(selMap, maskTopoNoCosPA | maskAntiLambdaSpecific))
          histos.fill(HIST("AntiLambda/h4dPointingAngle"), collision.centFT0C(), v0.pt(), v0.mAntiLambda(), TMath::ACos(v0.v0cosPA()));
        if (verifyMask(selMap, maskTopoNoDCAV0Dau | maskAntiLambdaSpecific))
          histos.fill(HIST("AntiLambda/h4dDCADaughters"), collision.centFT0C(), v0.pt(), v0.mAntiLambda(), v0.dcaV0daughters());
      }
    } // end systematics / qa
  }

  // ______________________________________________________
  // Real data processing - no MC subscription
  void processRealData(soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels>::iterator const& collision, v0Candidates const& fullV0s, dauTracks const&)
  {
    histos.fill(HIST("hEventSelection"), 0. /* all collisions */);
    if (!collision.sel8()) {
      return;
    }
    histos.fill(HIST("hEventSelection"), 1 /* sel8 collisions */);

    if (std::abs(collision.posZ()) > 10.f) {
      return;
    }
    histos.fill(HIST("hEventSelection"), 2 /* vertex-Z selected */);
    histos.fill(HIST("hEventCentrality"), collision.centFT0C());

    // __________________________________________
    // perform main analysis
    for (auto& v0 : fullV0s) {
      if (std::abs(v0.negativeeta()) > daughterEtaCut || std::abs(v0.positiveeta()) > daughterEtaCut)
        continue; // remove acceptance that's badly reproduced by MC / superfluous in future

      // fill AP plot for all V0s
      histos.fill(HIST("GeneralQA/h2dArmenterosAll"), v0.alpha(), v0.qtarm());

      uint32_t selMap = computeReconstructionBitmap(v0, collision);

      // consider for histograms for all species
      selMap = selMap | (1 << selConsiderK0Short) | (1 << selConsiderLambda) | (1 << selConsiderAntiLambda);

      analyseCandidate(v0, collision, selMap);
    } // end v0 loop
  }

  // ______________________________________________________
  // Simulated processing (subscribes to MC information too)
  void processMonteCarlo(soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels>::iterator const& collision, v0MCCandidates const& fullV0s, dauTracks const&)
  {
    histos.fill(HIST("hEventSelection"), 0. /* all collisions */);
    if (!collision.sel8()) {
      return;
    }
    histos.fill(HIST("hEventSelection"), 1 /* sel8 collisions */);

    if (std::abs(collision.posZ()) > 10.f) {
      return;
    }
    histos.fill(HIST("hEventSelection"), 2 /* vertex-Z selected */);
    histos.fill(HIST("hEventCentrality"), collision.centFT0C());

    // __________________________________________
    // perform main analysis
    for (auto& v0 : fullV0s) {
      if (std::abs(v0.negativeeta()) > daughterEtaCut || std::abs(v0.positiveeta()) > daughterEtaCut)
        continue; // remove acceptance that's badly reproduced by MC / superfluous in future

      // fill AP plot for all V0s
      histos.fill(HIST("GeneralQA/h2dArmenterosAll"), v0.alpha(), v0.qtarm());

      uint32_t selMap = computeReconstructionBitmap(v0, collision);

      // consider only associated candidates if asked to do so
      if (doMCAssociation) {
        selMap = selMap | computeMCAssociation(v0);
      } else {
        selMap = selMap | (1 << selConsiderK0Short) | (1 << selConsiderLambda) | (1 << selConsiderAntiLambda);
      }

      analyseCandidate(v0, collision, selMap);
    } // end v0 loop
  }

  PROCESS_SWITCH(derivedlambdakzeroanalysis, processRealData, "process as if real data", true);
  PROCESS_SWITCH(derivedlambdakzeroanalysis, processMonteCarlo, "process as if MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<derivedlambdakzeroanalysis>(cfgc)};
}
