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
///
/// \file alice3Strangeness.cxx
///
/// \brief This task produces invariant mass distributions for strange hadrons
///
/// \author Lucia Anna Tarasovičová, Pavol Jozef Šafárik University (SK)
/// \since  November 20, 2025
///

#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "ALICE3/DataModel/OTFCollision.h"
#include "ALICE3/DataModel/OTFStrangeness.h"
#include "ALICE3/DataModel/tracksAlice3.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/StaticFor.h>
#include <Framework/runDataProcessing.h>

#include <Rtypes.h>

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <string>
#include <string_view>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::constants::math;

using Alice3Tracks = soa::Join<aod::Tracks, aod::TracksCov, aod::McTrackLabels, aod::TracksDCA, aod::TracksExtraA3>;
using FullV0Candidates = soa::Join<aod::V0CandidateIndices, aod::V0CandidateCores>;
using FullCascadeCandidates = soa::Join<aod::StoredCascCores, aod::CascIndices, aod::A3CascadeMcLabels>;
using FullCollisions = soa::Join<aod::OTFLUTConfigId, aod::Collisions>;

struct Alice3Strangeness {
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<int> idGeometry{"idGeometry", 0, "geometry ID used for propagation"};
  struct : ConfigurableGroup {
    ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "pt axis for QA histograms"};
    ConfigurableAxis axisK0Mass{"axisK0Mass", {200, 0.4f, 0.6f}, ""};
    ConfigurableAxis axisLambdaMass{"axisLambdaMass", {200, 1.08f, 1.2f}, ""};
    ConfigurableAxis axisXiMass{"axisXiMass", {200, 1.22f, 1.42f}, ""};
    ConfigurableAxis axisOmegaMass{"axisOmegaMass", {200, 1.57f, 1.77f}, ""};
    ConfigurableAxis axisVertexZ{"axisVertexZ", {40, -20, 20}, "vertex Z (cm)"};
    ConfigurableAxis axisDCA{"axisDCA", {200, 0, 5}, "DCA (cm)"};
    ConfigurableAxis axisRadius{"axisRadius", {50, 0.0, 100}, "V0 radius (cm)"};
    ConfigurableAxis axisDCAV0Daughters{"axisDCAV0Daughters", {20, 0, 5}, "DCA V0 daughters"};
    ConfigurableAxis axisPointingAngle{"axisPointingAngle", {40, 0.0f, 0.4f}, "pointing angle "};
    ConfigurableAxis axisCosPA{"axisCosPA", {100, 0.91f, 1.01f}, "cos(pointing angle) "};
    ConfigurableAxis axisProperLifeTime{"axisProperLifeTime", {100, 0.0f, 100.0f}, "proper lifetime (cm)"};
    ConfigurableAxis axisNormalizedDecayLength{"axisNormalizedDecayLength", {100, 0.0f, 20}, "Normalized decay length"};
    ConfigurableAxis axisEta{"axisEta", {100, -5.0f, 5.0f}, "eta"};
  } histAxes;

  struct : ConfigurableGroup {
    std::string prefix = "v0SelectionFlags";
    Configurable<bool> applyRapiditySelection{"applyRapiditySelection", true, "apply rapidity selection"};
    Configurable<bool> applyDCAdaughterSelection{"applyDCAdaughterSelection", true, "apply DCA daughter selection"};
    Configurable<bool> applyCosOfPAngleSelection{"applyCosOfPAngleSelection", true, "apply cosine of pointing angle selection"};
    Configurable<bool> applyDCADaughtersToPVSelection{"applyDCADaughtersToPVSelection", true, "apply DCA daughters to primary vertex selection"};
    Configurable<bool> applyV0RadiusSelection{"applyV0RadiusSelection", true, "apply V0 radius selection"};
    Configurable<bool> applyArmenterosSelection{"applyArmenterosSelection", true, "apply Armenteros-Podolanski selection"};
    Configurable<bool> applyCompetingMassRejection{"applyCompetingMassRejection", true, "apply competing mass rejection"};
    Configurable<bool> applyLifetimeSelection{"applyLifetimeSelection", true, "apply lifetime selection"};
    Configurable<bool> applyEtaDaughterSelection{"applyEtaDaughterSelection", true, "apply eta daughter selection"};
    Configurable<bool> doQAforSelectionVariables{"doQAforSelectionVariables", false, "enable QA plots"};
    Configurable<bool> analyseOnlyTrueV0s{"analyseOnlyTrueV0s", false, "analyse only true V0s from MC"};
  } v0SelectionFlags;

  struct : ConfigurableGroup {
    std::string prefix = "v0SelectionValues";
    Configurable<float> yK0Selection{"yK0Selection", 0.5f, "rapidity selection for K0"};
    Configurable<float> yLambdaSelection{"yLambdaSelection", 0.5f, "rapidity selection for Lambda"};
    Configurable<float> dcaDaughterSelection{"dcaDaughterSelection", 1.0f, "DCA daughter selection"};
    Configurable<float> cosPAngleSelection{"cosPAngleSelection", 0.97f, "cosine of pointing angle selection"};
    Configurable<float> dcaDaughtersToPVSelection{"dcaDaughtersToPVSelection", 0.05f, "DCA daughters to primary vertex selection"};
    Configurable<float> v0RadiusSelection{"v0RadiusSelection", 0.5f, "V0 radius selection"};
    Configurable<float> armenterosSelection{"armenterosSelection", 0.2f, "Armenteros-Podolanski selection in [0]* alpha"};
    Configurable<float> competingMassRejectionK0{"competingMassRejectionK0", 0.005f, "competing mass rejection for K0"};
    Configurable<float> competingMassRejectionLambda{"competingMassRejectionLambda", 0.01f, "competing mass rejection for Lambda"};
    Configurable<float> lifetimecutambda{"lifetimecutambda", 30, "lifetime cut for Lambda in cm"};
    Configurable<float> lifetimecutak0{"lifetimecutak0", 20, "lifetime cut for K0 in cm"};
    Configurable<float> etaDaughterSelection{"etaDaughterSelection", 0.8f, "eta daughter selection"};
    Configurable<float> acceptedLambdaMassWindow{"acceptedLambdaMassWindow", 0.2f, "accepted Lambda mass window around PDG mass"};
    Configurable<float> acceptedK0MassWindow{"acceptedK0MassWindow", 0.3f, "accepted K0 mass window around PDG mass"};
  } v0SelectionValues;

  struct : ConfigurableGroup {
    std::string prefix = "cascadeSelectionValues";
    Configurable<float> posMinConstDCAxy{"posMinConstDCAxy", 0.5f, "[0] in |DCAxy| > [0]+[1]/pT"};
    Configurable<float> posMinPtDepDCAxy{"posMinPtDepDCAxy", 0.f, "[1] in |DCAxy| > [0]+[1]/pT"};
    Configurable<float> posMinConstDCAz{"posMinConstDCAz", -1.f, "[1] in |DCAz| > [0]+[1]/pT"};
    Configurable<float> posMinPtDepDCAz{"posMinPtDepDCAz", 0.f, "[0] in |DCAz| > [0]+[1]/pT"};
    Configurable<float> negMinConstDCAxy{"negMinConstDCAxy", 0.5f, "[0] in |DCAxy| > [0]+[1]/pT"};
    Configurable<float> negMinPtDepDCAxy{"negMinPtDepDCAxy", 0.f, "[1] in |DCAxy| > [0]+[1]/pT"};
    Configurable<float> negMinConstDCAz{"negMinConstDCAz", -1.f, "[1] in |DCAz| > [0]+[1]/pT"};
    Configurable<float> negMinPtDepDCAz{"negMinPtDepDCAz", 0.f, "[0] in |DCAz| > [0]+[1]/pT"};
    Configurable<float> bachMinConstDCAxy{"bachMinConstDCAxy", 0.5f, "[0] in |DCAxy| > [0]+[1]/pT"};
    Configurable<float> bachMinPtDepDCAxy{"bachMinPtDepDCAxy", 0.f, "[1] in |DCAxy| > [0]+[1]/pT"};
    Configurable<float> bachMinConstDCAz{"bachMinConstDCAz", -1.f, "[1] in |DCAz| > [0]+[1]/pT"};
    Configurable<float> bachMinPtDepDCAz{"bachMinPtDepDCAz", 0.f, "[0] in |DCAz| > [0]+[1]/pT"};
    Configurable<float> laMaxDauDCA{"laMaxDauDCA", 0.5f, "DCA (cm) between lambda daughters"};
    Configurable<float> laMinDecayRadius{"laMinDecayRadius", 0.5f, "Minimum lambda radius"};
    Configurable<float> laMassWindow{"laMassWindow", 0.5f, "accepted la mass window around PDG mass"};
    Configurable<float> laMinCosPA{"laMinCosPA", 0.98, "Min Lambda CosPA"};
    Configurable<float> cascMaxDauDCA{"cascMaxDauDCA", 0.5f, "DCA (cm) between cascade daughters"};
    Configurable<float> cascMaxNormalizedDecayLength{"cascMaxNormalizedDecayLength", 3, "Max cascade nomralized decay length (ctau/<ctau>)"};
    Configurable<float> cascMinCosPA{"cascMinCosPA", 0.98, "Minimum cascade CosPA"};
    Configurable<float> cascMinDecayRadius{"cascMinDecayRadius", 0.5f, "Minimum cascade decay radius"};
    Configurable<float> cascMassWindow{"cascMassWindow", 0.5f, "accepted cascade mass window around PDG mass"};
    Configurable<float> competingMassRejection{"competingMassRejection", 0.5f, "competing mass rejection"};
  } cascadeSelectionValues;

  struct : ConfigurableGroup {
    std::string prefix = "cascadeFlags";
    Configurable<int> analyseCascade{"analyseCascade", 0, "0: Xi, 1: AntiXi, 2: Omega, 3: AntiOmega"};
    Configurable<bool> analyseOnlyTrueCascades{"analyseOnlyTrueCascades", false, "analyse only true cascades from MC"};
    Configurable<bool> posDCAxy{"posDCAxy", true, "enable posDCAxy selection"};
    Configurable<bool> posDCAz{"posDCAz", false, "enable posDCAz selection"};
    Configurable<bool> negDCAxy{"negDCAxy", true, "enable negDCAxy selection"};
    Configurable<bool> negDCAz{"negDCAz", false, "enable negDCAz selection"};
    Configurable<bool> bachDCAxy{"bachDCAxy", true, "enable bachDCAxy selection"};
    Configurable<bool> bachDCAz{"bachDCAz", false, "enable bachDCAz selection"};
    Configurable<bool> laMaxDauDCA{"laMaxDauDCA", true, "enable laMaxDauDCA selection"};
    Configurable<bool> laMinDecayRadius{"laMinDecayRadius", true, "enable laMinDecayRadius selection"};
    Configurable<bool> laMassWindow{"laMassWindow", true, "enable laMassWindow selection"};
    Configurable<bool> laMinCosPA{"laMinCosPA", true, "enable laMinCosPA selection"};
    Configurable<bool> cascMaxDauDCA{"cascMaxDauDCA", true, "enable cascMaxDauDCA selection"};
    Configurable<bool> cascMinDecayRadius{"cascMinDecayRadius", true, "enable cascMinDecayRadius selection"};
    Configurable<bool> cascMaxNormalizedDecayLength{"cascMaxNormalizedDecayLength", true, "enable cascMaxNormalizedDecayLength selection"};
    Configurable<bool> cascMinCosPA{"cascMinCosPA", true, "enable cascMinCosPA selection"};
    Configurable<bool> competingMassRejection{"competingMassRejection", false, "enable competingMassRejection selection"};
  } cascadeFlags;

  uint16_t appliedSelectionCheckMask;
  double selectionCheck;
  double selectionCheckPos;
  const int posDaugDCAselIDx = 3;
  static constexpr std::string_view KSelectionNames[] = {"DCAV0Daughters", "PointingAngle", "DCAtoPVNegDaughter", "DCAtoPVPosDaughter", "V0Radius", "ProperLifeTime"};

  static constexpr float ToMicrons = 1e+4f;
  static constexpr float CtauXi = 4.91f;
  static constexpr float CtauOmega = 2.461f;

  struct Cascade {
    enum Type { Xi = 0,
                AntiXi,
                Omega,
                AntiOmega };
    Type type{};

    float pdgMass{}, pdgCompetingMass{}, ctau{};
    int sign{}, pdgCode{};

    void setCascadeType(Type newType)
    {
      type = newType;
      if (type == Xi || type == AntiXi) {
        ctau = CtauXi;
        pdgMass = o2::constants::physics::MassXiMinus;
        pdgCompetingMass = o2::constants::physics::MassOmegaMinus;
        pdgCode = (type == Xi) ? PDG_t::kXiMinus : kXiPlusBar;
        sign = (type == Xi) ? -1 : 1;
      }

      if (type == Omega || type == AntiOmega) {
        ctau = CtauOmega;
        pdgMass = o2::constants::physics::MassOmegaMinus;
        pdgCompetingMass = o2::constants::physics::MassXiMinus;
        pdgCode = (type == Omega) ? PDG_t::kOmegaMinus : kOmegaPlusBar;
        sign = (type == Omega) ? -1 : 1;
      }
    }
  } analysedCascade;

  void init(InitContext&)
  {
    histos.add("K0/hMassAllCandidates", "", kTH2D, {histAxes.axisK0Mass, histAxes.axisPt});
    histos.add("K0/hMassSelected", "", kTH2D, {histAxes.axisK0Mass, histAxes.axisPt});
    histos.add("K0/hSelections", "", kTH1D, {{10, 0, 10}});
    histos.add("K0/hDCANegDaughter", "", kTH1D, {{200, -5, 5}});
    histos.add("K0/hDCAPosDaughter", "", kTH1D, {{200, -5, 5}});
    histos.add("Xi/hMassAllCandidates", "hMassAllCandidates", kTH1D, {histAxes.axisXiMass});
    histos.add("Xi/hMassSelected", "hMassSelected", kTH1D, {histAxes.axisXiMass});

    histos.add("hPVz", "hPVz", kTH1F, {histAxes.axisVertexZ});
    histos.add("hV0CandidateCounter", "hV0CandidateCounter", kTH1F, {{11, 0, 11}});
    histos.add("reconstructedCandidates/hEtaDaughters", "hEtaDaughters", kTH1F, {histAxes.axisEta});
    histos.add("reconstructedCandidates/K0/hMass", "hMass", kTH3D, {histAxes.axisK0Mass, histAxes.axisPt, histAxes.axisEta});
    histos.add("reconstructedCandidates/K0/hMass1D", "hMass1D", kTH1D, {histAxes.axisK0Mass});
    histos.add("reconstructedCandidates/Lambda/hMass", "hMass", kTH3D, {histAxes.axisLambdaMass, histAxes.axisPt, histAxes.axisEta});
    histos.add("reconstructedCandidates/Lambda/hMass1D", "hMass1D", kTH1D, {histAxes.axisLambdaMass});
    histos.add("reconstructedCandidates/hArmeterosBeforeAllSelections", "hArmeterosBeforeAllSelections", kTH2D, {{100, -1.0f, 1.0f}, {200, 0.0f, 0.5f}});
    histos.add("reconstructedCandidates/hArmeterosAfterAllSelections", "hArmeterosAfterAllSelections", kTH2D, {{100, -1.0f, 1.0f}, {200, 0.0f, 0.5f}});

    if (doprocessFoundCascadeCandidates) {
      analysedCascade.setCascadeType(static_cast<Cascade::Type>(cascadeFlags.analyseCascade.value));
      histos.add("reconstructedCandidates/Cascade/hMassAllXiCandidates", "hMassAllXiCandidates", kTH1D, {histAxes.axisXiMass});
      histos.add("reconstructedCandidates/Cascade/hMassAllAntiXiCandidates", "hMassAllAntiXiCandidates", kTH1D, {histAxes.axisXiMass});
      histos.add("reconstructedCandidates/Cascade/hMassAllOmegaCandidates", "hMassAllOmegaCandidates", kTH1D, {histAxes.axisOmegaMass});
      histos.add("reconstructedCandidates/Cascade/hMassAllAntiOmegaCandidates", "hMassAllAntiOmegaCandidates", kTH1D, {histAxes.axisOmegaMass});

      histos.add("reconstructedCandidates/Cascade/hMassSelectedXiCandidates", "hMassSelectedXiCandidates", kTH1D, {histAxes.axisXiMass});
      histos.add("reconstructedCandidates/Cascade/hMassSelectedAntiXiCandidates", "hMassSelectedAntiXiCandidates", kTH1D, {histAxes.axisXiMass});
      histos.add("reconstructedCandidates/Cascade/hMassSelectedOmegaCandidates", "hMassSelectedOmegaCandidates", kTH1D, {histAxes.axisOmegaMass});
      histos.add("reconstructedCandidates/Cascade/hMassSelectedAntiOmegaCandidates", "hMassSelectedAntiOmegaCandidates", kTH1D, {histAxes.axisOmegaMass});

      histos.add("reconstructedCandidates/Cascade/h3dXiCandidates", "h3dXiCandidates", kTH3D, {histAxes.axisPt, histAxes.axisEta, histAxes.axisXiMass});
      histos.add("reconstructedCandidates/Cascade/h3dAntiXiCandidates", "h3dAntiXiCandidates", kTH3D, {histAxes.axisPt, histAxes.axisEta, histAxes.axisXiMass});
      histos.add("reconstructedCandidates/Cascade/h3dOmegaCandidates", "h3dOmegaCandidates", kTH3D, {histAxes.axisPt, histAxes.axisEta, histAxes.axisOmegaMass});
      histos.add("reconstructedCandidates/Cascade/h3dAntiOmegaCandidates", "h3dAntiOmegaCandidates", kTH3D, {histAxes.axisPt, histAxes.axisEta, histAxes.axisOmegaMass});

      histos.add("reconstructedCandidates/Cascade/hSelectionQa", "hSelectionQa", kTH1D, {{20, 0.5, 20.5}});
      histos.get<TH1>(HIST("reconstructedCandidates/Cascade/hSelectionQa"))->GetXaxis()->SetBinLabel(1, "all candidates");
      histos.get<TH1>(HIST("reconstructedCandidates/Cascade/hSelectionQa"))->GetXaxis()->SetBinLabel(2, "pos dcaXY");
      histos.get<TH1>(HIST("reconstructedCandidates/Cascade/hSelectionQa"))->GetXaxis()->SetBinLabel(3, "neg dcaXY");
      histos.get<TH1>(HIST("reconstructedCandidates/Cascade/hSelectionQa"))->GetXaxis()->SetBinLabel(4, "bach dcaXY");
      histos.get<TH1>(HIST("reconstructedCandidates/Cascade/hSelectionQa"))->GetXaxis()->SetBinLabel(5, "pos dcaZ");
      histos.get<TH1>(HIST("reconstructedCandidates/Cascade/hSelectionQa"))->GetXaxis()->SetBinLabel(6, "neg dcaZ");
      histos.get<TH1>(HIST("reconstructedCandidates/Cascade/hSelectionQa"))->GetXaxis()->SetBinLabel(7, "bach dcaZ");
      histos.get<TH1>(HIST("reconstructedCandidates/Cascade/hSelectionQa"))->GetXaxis()->SetBinLabel(8, "la dau dca");
      histos.get<TH1>(HIST("reconstructedCandidates/Cascade/hSelectionQa"))->GetXaxis()->SetBinLabel(9, "la decay radius");
      histos.get<TH1>(HIST("reconstructedCandidates/Cascade/hSelectionQa"))->GetXaxis()->SetBinLabel(10, "la mass window");
      histos.get<TH1>(HIST("reconstructedCandidates/Cascade/hSelectionQa"))->GetXaxis()->SetBinLabel(11, "la cosPA");
      histos.get<TH1>(HIST("reconstructedCandidates/Cascade/hSelectionQa"))->GetXaxis()->SetBinLabel(12, "casc dau dca");
      histos.get<TH1>(HIST("reconstructedCandidates/Cascade/hSelectionQa"))->GetXaxis()->SetBinLabel(13, "casc norm decay length");
      histos.get<TH1>(HIST("reconstructedCandidates/Cascade/hSelectionQa"))->GetXaxis()->SetBinLabel(14, "casc decay radius");
      histos.get<TH1>(HIST("reconstructedCandidates/Cascade/hSelectionQa"))->GetXaxis()->SetBinLabel(15, "casc cosPA");
      histos.get<TH1>(HIST("reconstructedCandidates/Cascade/hSelectionQa"))->GetXaxis()->SetBinLabel(16, "competing mass");

      histos.add("reconstructedCandidates/Cascade/BeforeSelection/hPosDCAxy", "hPosDCAxy", kTH1D, {histAxes.axisDCA});
      histos.add("reconstructedCandidates/Cascade/BeforeSelection/hNegDCAxy", "hNegDCAxy", kTH1D, {histAxes.axisDCA});
      histos.add("reconstructedCandidates/Cascade/BeforeSelection/hBachDCAxy", "hBachDCAxy", kTH1D, {histAxes.axisDCA});
      histos.add("reconstructedCandidates/Cascade/BeforeSelection/hPosDCAz", "hPosDCAz", kTH1D, {histAxes.axisDCA});
      histos.add("reconstructedCandidates/Cascade/BeforeSelection/hNegDCAz", "hNegDCAz", kTH1D, {histAxes.axisDCA});
      histos.add("reconstructedCandidates/Cascade/BeforeSelection/hBachDCAz", "hBachDCAz", kTH1D, {histAxes.axisDCA});
      histos.add("reconstructedCandidates/Cascade/BeforeSelection/hLaDauDCA", "hLaDauDCA", kTH1D, {histAxes.axisDCA});
      histos.add("reconstructedCandidates/Cascade/BeforeSelection/hLaDecayRadius", "hLaDecayRadius", kTH1D, {histAxes.axisRadius});
      histos.add("reconstructedCandidates/Cascade/BeforeSelection/hLaMassWindow", "hLaMassWindow", kTH1D, {histAxes.axisLambdaMass});
      histos.add("reconstructedCandidates/Cascade/BeforeSelection/hLaCosPA", "hLaCosPA", kTH1D, {histAxes.axisCosPA});
      histos.add("reconstructedCandidates/Cascade/BeforeSelection/hCascDauDCA", "hCascDauDCA", kTH1D, {histAxes.axisDCA});
      histos.add("reconstructedCandidates/Cascade/BeforeSelection/hCascDecayLength", "hCascDecayLength", kTH1D, {histAxes.axisNormalizedDecayLength});
      histos.add("reconstructedCandidates/Cascade/BeforeSelection/hCascDecayRadius", "hCascDecayRadius", kTH1D, {histAxes.axisRadius});
      histos.add("reconstructedCandidates/Cascade/BeforeSelection/hCascCosPA", "hCascCosPA", kTH1D, {histAxes.axisCosPA});
      histos.addClone("reconstructedCandidates/Cascade/BeforeSelection/", "reconstructedCandidates/Cascade/AfterSelection/");
    }

    if (v0SelectionFlags.doQAforSelectionVariables) {
      if (!v0SelectionFlags.applyDCADaughtersToPVSelection) {
        histos.add("reconstructedCandidates/K0/hDCAtoPVNegDaughter", "hDCAtoPVNegDaughter", kTH3D, {histAxes.axisK0Mass, histAxes.axisPt, histAxes.axisDCA});
        histos.add("reconstructedCandidates/K0/hDCAtoPVPosDaughter", "hDCAtoPVPosDaughter", kTH3D, {histAxes.axisK0Mass, histAxes.axisPt, histAxes.axisDCA});
        histos.add("reconstructedCandidates/Lambda/hDCAtoPVNegDaughter", "hDCAtoPVNegDaughter", kTH3D, {histAxes.axisLambdaMass, histAxes.axisPt, histAxes.axisDCA});
        histos.add("reconstructedCandidates/Lambda/hDCAtoPVPosDaughter", "hDCAtoPVPosDaughter", kTH3D, {histAxes.axisLambdaMass, histAxes.axisPt, histAxes.axisDCA});
      }
      if (!v0SelectionFlags.applyV0RadiusSelection) {
        histos.add("reconstructedCandidates/K0/hV0Radius", "hV0Radius", kTH3D, {histAxes.axisK0Mass, histAxes.axisPt, histAxes.axisRadius});
        histos.add("reconstructedCandidates/Lambda/hV0Radius", "hV0Radius", kTH3D, {histAxes.axisLambdaMass, histAxes.axisPt, histAxes.axisRadius});
      }
      if (!v0SelectionFlags.applyDCAdaughterSelection) {
        histos.add("reconstructedCandidates/K0/hDCAV0Daughters", "hDCAV0Daughters", kTH3D, {histAxes.axisK0Mass, histAxes.axisPt, histAxes.axisDCAV0Daughters});
        histos.add("reconstructedCandidates/Lambda/hDCAV0Daughters", "hDCAV0Daughters", kTH3D, {histAxes.axisLambdaMass, histAxes.axisPt, histAxes.axisDCAV0Daughters});
      }
      if (!v0SelectionFlags.applyCosOfPAngleSelection) {
        histos.add("reconstructedCandidates/K0/hPointingAngle", "hPointingAngle", kTH3D, {histAxes.axisK0Mass, histAxes.axisPt, histAxes.axisPointingAngle});
        histos.add("reconstructedCandidates/Lambda/hPointingAngle", "hPointingAngle", kTH3D, {histAxes.axisLambdaMass, histAxes.axisPt, histAxes.axisPointingAngle});
      }
      if (!v0SelectionFlags.applyLifetimeSelection) {
        histos.add("reconstructedCandidates/K0/hProperLifeTime", "hProperLifeTime", kTH3D, {histAxes.axisK0Mass, histAxes.axisPt, histAxes.axisProperLifeTime});
        histos.add("reconstructedCandidates/Lambda/hProperLifeTime", "hProperLifeTime", kTH3D, {histAxes.axisLambdaMass, histAxes.axisPt, histAxes.axisProperLifeTime});
      }
    }
    histos.addClone("reconstructedCandidates/Lambda/", "reconstructedCandidates/AntiLambda/");

    appliedSelectionCheckMask = 0;
    if (!v0SelectionFlags.applyDCAdaughterSelection) {
      SETBIT(appliedSelectionCheckMask, 0);
    }
    if (!v0SelectionFlags.applyCosOfPAngleSelection) {
      SETBIT(appliedSelectionCheckMask, 1);
    }
    if (!v0SelectionFlags.applyDCADaughtersToPVSelection) {
      SETBIT(appliedSelectionCheckMask, 2);
      SETBIT(appliedSelectionCheckMask, 3);
    }
    if (!v0SelectionFlags.applyV0RadiusSelection) {
      SETBIT(appliedSelectionCheckMask, 4);
    }
    if (!v0SelectionFlags.applyLifetimeSelection) {
      SETBIT(appliedSelectionCheckMask, 5);
    }
    histos.print();
  }

  void processAllFindableCandidates(aod::Collisions const& collisions, aod::McCollisions const&, aod::UpgradeV0s const& v0Recos, aod::UpgradeCascades const& cascRecos, Alice3Tracks const&)
  {
    for (const auto& collision : collisions) {
      float collisionZ = collision.posZ();
      // std::cout << "______ process V0_______" <<  collision.size() << std::endl;
      histos.fill(HIST("hPVz"), collisionZ);
    }

    for (const auto& v0Cand : v0Recos) {
      auto negV0Daughter = v0Cand.negTrack_as<Alice3Tracks>(); // de-reference negative track
      auto posV0Daughter = v0Cand.posTrack_as<Alice3Tracks>(); // de-reference positive track

      bool isK0 = v0Cand.mK0() > 0;
      if (isK0) {
        histos.fill(HIST("K0/hMassAllCandidates"), v0Cand.mK0(), v0Cand.pt());
        histos.fill(HIST("K0/hSelections"), 0); // all candidates
        histos.fill(HIST("K0/hDCANegDaughter"), negV0Daughter.dcaXY());
        histos.fill(HIST("K0/hDCAPosDaughter"), posV0Daughter.dcaXY());
        if (std::abs(negV0Daughter.dcaXY()) < v0SelectionValues.dcaDaughtersToPVSelection)
          continue;
        histos.fill(HIST("K0/hSelections"), 1); // dcaXY cut
        if (std::abs(posV0Daughter.dcaXY()) < v0SelectionValues.dcaDaughtersToPVSelection)
          continue;
        histos.fill(HIST("K0/hSelections"), 2); // dcaXY cut
        if (v0Cand.dcaV0Daughters() > v0SelectionValues.dcaDaughterSelection)
          continue;
        histos.fill(HIST("K0/hSelections"), 3); // dca between daughters
        if (v0Cand.v0Radius() < v0SelectionValues.v0RadiusSelection)
          continue;
        histos.fill(HIST("K0/hSelections"), 4); // radius cut
        if (std::abs(negV0Daughter.eta()) > v0SelectionValues.etaDaughterSelection || std::abs(posV0Daughter.eta()) > v0SelectionValues.etaDaughterSelection)
          continue;
        histos.fill(HIST("K0/hSelections"), 5); // eta cut
        histos.fill(HIST("K0/hMassSelected"), v0Cand.mK0(), v0Cand.pt());
      }
    }

    for (const auto& cascCand : cascRecos) {
      // auto bachelor = cascCand.bachTrack_as<Alice3Tracks>(); // de-reference bachelor track
      // auto negative = cascCand.negTrack_as<Alice3Tracks>();   // de-reference negative track
      // auto positive = cascCand.posTrack_as<Alice3Tracks>();   // de-reference positive track

      // Only XiMinus in the tracker for now
      histos.fill(HIST("Xi/hMassAllCandidates"), cascCand.mXi());

      // TODO Add selections
      histos.fill(HIST("Xi/hMassSelected"), cascCand.mXi());
    }
  }

  void processFoundV0Candidates(aod::Collision const& collision, FullV0Candidates const& v0Candidates, Alice3Tracks const&, aod::McParticles const&)
  {
    // if(collision.lutConfigId()!=idGeometry)
    // return;
    float collisionZ = collision.posZ();
    histos.fill(HIST("hPVz"), collisionZ);
    for (auto const& v0 : v0Candidates) {
      bool isK0 = std::abs(v0.mK0Short() - o2::constants::physics::MassK0Short) < v0SelectionValues.acceptedK0MassWindow;
      bool isLambda = std::abs(v0.mLambda() - o2::constants::physics::MassLambda0) < v0SelectionValues.acceptedLambdaMassWindow;
      bool isAntiLambda = std::abs(v0.mAntiLambda() - o2::constants::physics::MassLambda0) < v0SelectionValues.acceptedLambdaMassWindow;

      histos.fill(HIST("reconstructedCandidates/hArmeterosBeforeAllSelections"), v0.alpha(), v0.qtArm());
      histos.fill(HIST("hV0CandidateCounter"), 0.5);
      if (v0SelectionFlags.applyRapiditySelection) {
        if (isK0 && std::abs(v0.yK0Short()) > v0SelectionValues.yK0Selection)
          continue;
        if ((isLambda || isAntiLambda) && std::abs(v0.yLambda()) > v0SelectionValues.yLambdaSelection)
          continue;
      }
      histos.fill(HIST("hV0CandidateCounter"), 1.5);
      if (v0SelectionFlags.applyDCAdaughterSelection) {
        if (std::abs(v0.dcaV0Daughters()) > v0SelectionValues.dcaDaughterSelection)
          continue;
      } else {
        selectionCheck = v0.dcaV0Daughters();
      }
      histos.fill(HIST("hV0CandidateCounter"), 2.5);
      if (v0SelectionFlags.applyCosOfPAngleSelection) {
        if (v0.cosPA() < v0SelectionValues.cosPAngleSelection)
          continue;
      } else {
        selectionCheck = std::acos(v0.cosPA());
      }
      histos.fill(HIST("hV0CandidateCounter"), 3.5);
      if (v0SelectionFlags.applyDCADaughtersToPVSelection) {
        if ((std::abs(v0.dcaNegToPV()) < v0SelectionValues.dcaDaughtersToPVSelection) ||
            (std::abs(v0.dcaPosToPV()) < v0SelectionValues.dcaDaughtersToPVSelection))
          continue;
      } else {
        selectionCheckPos = std::abs(v0.dcaPosToPV());
        selectionCheck = std::abs(v0.dcaNegToPV());
      }
      histos.fill(HIST("hV0CandidateCounter"), 4.5);
      if (v0SelectionFlags.applyV0RadiusSelection) {
        if (v0.v0radius() < v0SelectionValues.v0RadiusSelection)
          continue;
      } else {
        selectionCheck = v0.v0radius();
      }
      histos.fill(HIST("hV0CandidateCounter"), 5.5);
      if (isK0) {
        if (v0SelectionFlags.applyArmenterosSelection) {
          if (v0.qtArm() < v0SelectionValues.armenterosSelection * std::abs(v0.alpha()))
            continue;
        }
      }
      histos.fill(HIST("hV0CandidateCounter"), 6.5);
      if (isK0 && v0SelectionFlags.applyCompetingMassRejection) {
        if (std::abs(v0.mLambda() - o2::constants::physics::MassLambda0) < v0SelectionValues.competingMassRejectionK0)
          continue;
        if (std::abs(v0.mAntiLambda() - o2::constants::physics::MassLambda0) < v0SelectionValues.competingMassRejectionK0)
          continue;
      }
      if ((isLambda || isAntiLambda) && v0SelectionFlags.applyCompetingMassRejection) {
        if (std::abs(v0.mK0Short() - o2::constants::physics::MassK0Short) < v0SelectionValues.competingMassRejectionLambda)
          continue;
      }
      histos.fill(HIST("hV0CandidateCounter"), 7.5);
      if (v0SelectionFlags.applyLifetimeSelection) {
        if (isK0 && v0.distOverTotMom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassK0Short > v0SelectionValues.lifetimecutak0)
          continue;
        if ((isLambda || isAntiLambda) && v0.distOverTotMom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0 > v0SelectionValues.lifetimecutambda)
          continue;
      } else {
        if (isK0)
          selectionCheck = v0.distOverTotMom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassK0Short;
        else
          selectionCheck = v0.distOverTotMom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0;
      }
      histos.fill(HIST("hV0CandidateCounter"), 8.5);
      auto posTrack = v0.template posTrack_as<Alice3Tracks>();
      auto negTrack = v0.template negTrack_as<Alice3Tracks>();
      if (v0SelectionFlags.applyEtaDaughterSelection) {
        if (std::abs(posTrack.eta()) > v0SelectionValues.etaDaughterSelection || std::abs(negTrack.eta()) > v0SelectionValues.etaDaughterSelection)
          continue;
      }
      histos.fill(HIST("reconstructedCandidates/hEtaDaughters"), posTrack.eta());
      histos.fill(HIST("reconstructedCandidates/hEtaDaughters"), negTrack.eta());
      histos.fill(HIST("hV0CandidateCounter"), 9.5);

      histos.fill(HIST("reconstructedCandidates/hArmeterosAfterAllSelections"), v0.alpha(), v0.qtArm());
      if (v0SelectionFlags.doQAforSelectionVariables) {
        static_for<0, 5>([&](auto i) {
          constexpr int In = i.value;
          if (TESTBIT(appliedSelectionCheckMask, In)) {
            if (In == posDaugDCAselIDx)
              selectionCheck = selectionCheckPos;
            if (isK0)
              histos.fill(HIST("reconstructedCandidates/K0/h") + HIST(KSelectionNames[In]), v0.mK0Short(), v0.pt(), selectionCheck);
            if (isLambda)
              histos.fill(HIST("reconstructedCandidates/Lambda/h") + HIST(KSelectionNames[In]), v0.mLambda(), v0.pt(), selectionCheck);
            if (isAntiLambda)
              histos.fill(HIST("reconstructedCandidates/AntiLambda/h") + HIST(KSelectionNames[In]), v0.mAntiLambda(), v0.pt(), selectionCheck);
          }
        });
      }
      if (isK0) {
        histos.fill(HIST("reconstructedCandidates/K0/hMass"), v0.mK0Short(), v0.pt(), v0.eta());
        histos.fill(HIST("reconstructedCandidates/K0/hMass1D"), v0.mK0Short());
      }
      if (isLambda) {
        histos.fill(HIST("reconstructedCandidates/Lambda/hMass"), v0.mLambda(), v0.pt(), v0.eta());
        histos.fill(HIST("reconstructedCandidates/Lambda/hMass1D"), v0.mLambda());
      }
      if (isAntiLambda) {
        histos.fill(HIST("reconstructedCandidates/AntiLambda/hMass"), v0.mAntiLambda(), v0.pt(), v0.eta());
        histos.fill(HIST("reconstructedCandidates/AntiLambda/hMass1D"), v0.mAntiLambda());
      }
    }
  }

  void processFoundCascadeCandidates(aod::Collision const& collision, FullCascadeCandidates const& cascadeCandidates, Alice3Tracks const&, aod::McParticles const&)
  {
    for (const auto& cascade : cascadeCandidates) {
      if (cascade.sign() < 0) {
        histos.fill(HIST("reconstructedCandidates/Cascade/hMassAllXiCandidates"), cascade.mXi());
        histos.fill(HIST("reconstructedCandidates/Cascade/hMassAllOmegaCandidates"), cascade.mOmega());
      } else {
        histos.fill(HIST("reconstructedCandidates/Cascade/hMassAllAntiXiCandidates"), cascade.mXi());
        histos.fill(HIST("reconstructedCandidates/Cascade/hMassAllAntiOmegaCandidates"), cascade.mOmega());
      }

      auto positive = cascade.template posTrack_as<Alice3Tracks>();
      auto negative = cascade.template negTrack_as<Alice3Tracks>();
      auto bachelor = cascade.template bachelor_as<Alice3Tracks>();

      const float posDCAxyCut = cascadeSelectionValues.posMinConstDCAxy + cascadeSelectionValues.posMinPtDepDCAxy / positive.pt();
      const float negDCAxyCut = cascadeSelectionValues.negMinConstDCAxy + cascadeSelectionValues.negMinPtDepDCAxy / negative.pt();
      const float bachDCAxyCut = cascadeSelectionValues.bachMinConstDCAxy + cascadeSelectionValues.bachMinPtDepDCAxy / bachelor.pt();
      const float posDCAzCut = cascadeSelectionValues.posMinConstDCAz + cascadeSelectionValues.posMinPtDepDCAz / positive.pt();
      const float negDCAzCut = cascadeSelectionValues.negMinConstDCAz + cascadeSelectionValues.negMinPtDepDCAz / negative.pt();
      const float bachDCAzCut = cascadeSelectionValues.bachMinConstDCAz + cascadeSelectionValues.bachMinPtDepDCAz / bachelor.pt();
      const float distanceFromPV = std::hypot(cascade.x() - collision.posX(), cascade.y() - collision.posY(), cascade.z() - collision.posZ());
      const float normalizedDecayLength = analysedCascade.pdgMass * distanceFromPV / (cascade.p() * analysedCascade.ctau);

      const float cascadeCandidateMass = (analysedCascade.type == Cascade::Xi || analysedCascade.type == Cascade::AntiXi) ? cascade.mXi() : cascade.mOmega();
      const float cascadeCompetingCandidateMass = (analysedCascade.type == Cascade::Xi || analysedCascade.type == Cascade::AntiXi) ? cascade.mOmega() : cascade.mXi();
      bool isAnalysedCascade = (std::abs(cascadeCandidateMass - analysedCascade.pdgMass) < cascadeSelectionValues.cascMassWindow) && (analysedCascade.sign * cascade.sign() > 0);
      if (!isAnalysedCascade) {
        continue;
      }

      histos.fill(HIST("reconstructedCandidates/Cascade/BeforeSelection/hPosDCAxy"), positive.dcaXY() * ToMicrons);
      histos.fill(HIST("reconstructedCandidates/Cascade/BeforeSelection/hNegDCAxy"), negative.dcaXY() * ToMicrons);
      histos.fill(HIST("reconstructedCandidates/Cascade/BeforeSelection/hBachDCAxy"), bachelor.dcaXY() * ToMicrons);
      histos.fill(HIST("reconstructedCandidates/Cascade/BeforeSelection/hPosDCAz"), positive.dcaZ() * ToMicrons);
      histos.fill(HIST("reconstructedCandidates/Cascade/BeforeSelection/hNegDCAz"), negative.dcaZ() * ToMicrons);
      histos.fill(HIST("reconstructedCandidates/Cascade/BeforeSelection/hBachDCAz"), bachelor.dcaZ() * ToMicrons);
      histos.fill(HIST("reconstructedCandidates/Cascade/BeforeSelection/hLaDauDCA"), cascade.dcaV0daughters() * ToMicrons);
      histos.fill(HIST("reconstructedCandidates/Cascade/BeforeSelection/hLaDecayRadius"), cascade.v0radius());
      histos.fill(HIST("reconstructedCandidates/Cascade/BeforeSelection/hLaMassWindow"), cascade.mLambda());
      histos.fill(HIST("reconstructedCandidates/Cascade/BeforeSelection/hLaCosPA"), cascade.v0cosPA(collision.posX(), collision.posY(), collision.posZ()));
      histos.fill(HIST("reconstructedCandidates/Cascade/BeforeSelection/hCascDauDCA"), cascade.dcacascdaughters() * ToMicrons);
      histos.fill(HIST("reconstructedCandidates/Cascade/BeforeSelection/hCascDecayLength"), normalizedDecayLength);
      histos.fill(HIST("reconstructedCandidates/Cascade/BeforeSelection/hCascDecayRadius"), cascade.cascradius());
      histos.fill(HIST("reconstructedCandidates/Cascade/BeforeSelection/hCascCosPA"), cascade.casccosPA(collision.posX(), collision.posY(), collision.posZ()));

      if (cascadeFlags.analyseOnlyTrueCascades) {
        if (!cascade.has_mcParticle()) {
          continue;
        }

        auto mcCasc = cascade.template mcParticle_as<aod::McParticles>();
        if (mcCasc.pdgCode() != analysedCascade.pdgCode) {
          continue;
        }
      }

      histos.fill(HIST("reconstructedCandidates/Cascade/hSelectionQa"), 1 /* All candidates*/);
      if (cascadeFlags.posDCAxy && std::abs(positive.dcaXY()) < posDCAxyCut) {
        continue;
      }

      histos.fill(HIST("reconstructedCandidates/Cascade/hSelectionQa"), 2 /* Pass pos dcaXY*/);
      if (cascadeFlags.negDCAxy && std::abs(negative.dcaXY()) < negDCAxyCut) {
        continue;
      }

      histos.fill(HIST("reconstructedCandidates/Cascade/hSelectionQa"), 3 /* Pass neg dcaXY*/);
      if (cascadeFlags.bachDCAxy && std::abs(bachelor.dcaXY()) < bachDCAxyCut) {
        continue;
      }

      histos.fill(HIST("reconstructedCandidates/Cascade/hSelectionQa"), 4 /* Pass bach dcaXY*/);
      if (cascadeFlags.posDCAz && std::abs(positive.dcaZ()) < posDCAzCut) {
        continue;
      }

      histos.fill(HIST("reconstructedCandidates/Cascade/hSelectionQa"), 5 /* Pass pos dcaZ*/);
      if (cascadeFlags.negDCAz && std::abs(negative.dcaZ()) < negDCAzCut) {
        continue;
      }

      histos.fill(HIST("reconstructedCandidates/Cascade/hSelectionQa"), 6 /* Pass neg dcaZ*/);
      if (cascadeFlags.bachDCAz && std::abs(bachelor.dcaZ()) < bachDCAzCut) {
        continue;
      }

      histos.fill(HIST("reconstructedCandidates/Cascade/hSelectionQa"), 7 /* Pass bach dcaZ*/);
      if (cascadeFlags.laMaxDauDCA && cascade.dcaV0daughters() > cascadeSelectionValues.laMaxDauDCA) {
        continue;
      }

      histos.fill(HIST("reconstructedCandidates/Cascade/hSelectionQa"), 8 /* Pass la dau dca*/);
      if (cascadeFlags.laMinDecayRadius && cascade.v0radius() < cascadeSelectionValues.laMinDecayRadius) {
        continue;
      }

      histos.fill(HIST("reconstructedCandidates/Cascade/hSelectionQa"), 9 /* Pass la decay radius*/);
      if (cascadeFlags.laMassWindow && std::abs(cascade.mLambda() - o2::constants::physics::MassLambda0) > cascadeSelectionValues.laMassWindow) {
        continue;
      }

      histos.fill(HIST("reconstructedCandidates/Cascade/hSelectionQa"), 10 /* Pass  la mass window*/);
      if (cascadeFlags.laMinCosPA && cascade.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) < cascadeSelectionValues.laMinCosPA) {
        continue;
      }

      histos.fill(HIST("reconstructedCandidates/Cascade/hSelectionQa"), 11 /* Pass la cosPA*/);
      if (cascadeFlags.cascMaxDauDCA && cascade.dcacascdaughters() > cascadeSelectionValues.cascMaxDauDCA) {
        continue;
      }

      histos.fill(HIST("reconstructedCandidates/Cascade/hSelectionQa"), 12 /* Pass casc dau dca*/);
      if (cascadeFlags.cascMaxNormalizedDecayLength && normalizedDecayLength > cascadeSelectionValues.cascMaxNormalizedDecayLength) {
        continue;
      }

      histos.fill(HIST("reconstructedCandidates/Cascade/hSelectionQa"), 13 /* Pass casc norm decay length*/);
      if (cascadeFlags.cascMinDecayRadius && cascade.cascradius() < cascadeSelectionValues.cascMinDecayRadius) {
        continue;
      }

      histos.fill(HIST("reconstructedCandidates/Cascade/hSelectionQa"), 14 /* Pass casc decay radius*/);
      if (cascadeFlags.cascMinCosPA && cascade.casccosPA(collision.posX(), collision.posY(), collision.posZ()) < cascadeSelectionValues.cascMinCosPA) {
        continue;
      }

      histos.fill(HIST("reconstructedCandidates/Cascade/hSelectionQa"), 15 /* Pass casc cosPA*/);
      if (cascadeFlags.competingMassRejection && std::abs(cascadeCompetingCandidateMass - analysedCascade.pdgCompetingMass) < cascadeSelectionValues.competingMassRejection) {
        continue;
      }

      histos.fill(HIST("reconstructedCandidates/Cascade/hSelectionQa"), 16 /* Pass casc competing mass rej*/);
      histos.fill(HIST("reconstructedCandidates/Cascade/AfterSelection/hPosDCAxy"), positive.dcaXY() * ToMicrons);
      histos.fill(HIST("reconstructedCandidates/Cascade/AfterSelection/hNegDCAxy"), negative.dcaXY() * ToMicrons);
      histos.fill(HIST("reconstructedCandidates/Cascade/AfterSelection/hBachDCAxy"), bachelor.dcaXY() * ToMicrons);
      histos.fill(HIST("reconstructedCandidates/Cascade/AfterSelection/hPosDCAz"), positive.dcaZ() * ToMicrons);
      histos.fill(HIST("reconstructedCandidates/Cascade/AfterSelection/hNegDCAz"), negative.dcaZ() * ToMicrons);
      histos.fill(HIST("reconstructedCandidates/Cascade/AfterSelection/hBachDCAz"), bachelor.dcaZ() * ToMicrons);
      histos.fill(HIST("reconstructedCandidates/Cascade/AfterSelection/hLaDauDCA"), cascade.dcaV0daughters() * ToMicrons);
      histos.fill(HIST("reconstructedCandidates/Cascade/AfterSelection/hLaDecayRadius"), cascade.v0radius());
      histos.fill(HIST("reconstructedCandidates/Cascade/AfterSelection/hLaMassWindow"), cascade.mLambda());
      histos.fill(HIST("reconstructedCandidates/Cascade/AfterSelection/hLaCosPA"), cascade.v0cosPA(collision.posX(), collision.posY(), collision.posZ()));
      histos.fill(HIST("reconstructedCandidates/Cascade/AfterSelection/hCascDauDCA"), cascade.dcacascdaughters() * ToMicrons);
      histos.fill(HIST("reconstructedCandidates/Cascade/AfterSelection/hCascDecayLength"), normalizedDecayLength);
      histos.fill(HIST("reconstructedCandidates/Cascade/AfterSelection/hCascDecayRadius"), cascade.cascradius());
      histos.fill(HIST("reconstructedCandidates/Cascade/AfterSelection/hCascCosPA"), cascade.casccosPA(collision.posX(), collision.posY(), collision.posZ()));
      if (cascade.sign() < 0) {
        histos.fill(HIST("reconstructedCandidates/Cascade/hMassSelectedXiCandidates"), cascade.mXi());
        histos.fill(HIST("reconstructedCandidates/Cascade/hMassSelectedOmegaCandidates"), cascade.mOmega());
        histos.fill(HIST("reconstructedCandidates/Cascade/h3dXiCandidates"), cascade.pt(), cascade.eta(), cascade.mXi());
        histos.fill(HIST("reconstructedCandidates/Cascade/h3dOmegaCandidates"), cascade.pt(), cascade.eta(), cascade.mOmega());
      } else {
        histos.fill(HIST("reconstructedCandidates/Cascade/hMassSelectedAntiXiCandidates"), cascade.mXi());
        histos.fill(HIST("reconstructedCandidates/Cascade/hMassSelectedAntiOmegaCandidates"), cascade.mOmega());
        histos.fill(HIST("reconstructedCandidates/Cascade/h3dAntiXiCandidates"), cascade.pt(), cascade.eta(), cascade.mXi());
        histos.fill(HIST("reconstructedCandidates/Cascade/h3dAntiOmegaCandidates"), cascade.pt(), cascade.eta(), cascade.mOmega());
      }
    }
  }

  PROCESS_SWITCH(Alice3Strangeness, processAllFindableCandidates, "", false);
  PROCESS_SWITCH(Alice3Strangeness, processFoundV0Candidates, "", true);
  PROCESS_SWITCH(Alice3Strangeness, processFoundCascadeCandidates, "", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<Alice3Strangeness>(cfgc)};
}
