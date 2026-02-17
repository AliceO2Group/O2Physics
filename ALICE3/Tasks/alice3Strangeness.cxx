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

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <DCAFitter/DCAFitterN.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DetectorsBase/Propagator.h>
#include <DetectorsVertexing/PVertexer.h>
#include <DetectorsVertexing/PVertexerHelpers.h>
#include <Field/MagneticField.h>
#include <Framework/ASoA.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/StaticFor.h>
#include <ReconstructionDataFormats/DCA.h>
#include <SimulationDataFormat/InteractionSampler.h>

#include <TGenPhaseSpace.h>
#include <TGeoGlobalMagField.h>
#include <TPDGCode.h>
#include <TRandom3.h>

#include <RtypesCore.h>

#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::constants::math;

using Alice3Tracks = soa::Join<aod::Tracks, aod::TracksCov, aod::McTrackLabels, aod::TracksDCA, aod::TracksExtraA3>;
using FullV0Candidates = soa::Join<aod::V0CandidateIndices, aod::V0CandidateCores>;
using FullCascadeCandidates = soa::Join<aod::StoredCascCores, aod::CascIndices>;
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
    ConfigurableAxis axisV0Radius{"axisV0Radius", {50, 0.0, 100}, "V0 radius (cm)"};
    ConfigurableAxis axisDCAV0Daughters{"axisDCAV0Daughters", {20, 0, 5}, "DCA V0 daughters"};
    ConfigurableAxis axisPointingAngle{"axisPointingAngle", {40, 0.0f, 0.4f}, "pointing angle "};
    ConfigurableAxis axisProperLifeTime{"axisProperLifeTime", {100, 0.0f, 100.0f}, "proper lifetime (cm)"};
    ConfigurableAxis axisEta{"axisEta", {100, -5.0f, 5.0f}, "eta"};
  } histAxes;

  struct : ConfigurableGroup {
    std::string prefix = "selectionFlags";
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
  } selectionFlags;

  struct : ConfigurableGroup {
    std::string prefix = "selectionValues";
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
    Configurable<float> acceptedXiMassWindow{"acceptedXiMassWindow", 0.5f, "accepted Xi mass window around PDG mass"};
    Configurable<float> acceptedOmegaMassWindow{"acceptedOmegaMassWindow", 0.5f, "accepted Omega mass window around PDG mass"};
  } selectionValues;

  uint16_t appliedSelectionCheckMask;
  double selectionCheck;
  double selectionCheckPos;
  const int posDaugDCAselIDx = 3;
  static constexpr std::string_view KSelectionNames[] = {"DCAV0Daughters", "PointingAngle", "DCAtoPVNegDaughter", "DCAtoPVPosDaughter", "V0Radius", "ProperLifeTime"};

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
      histos.add("reconstructedCandidates/Xi/hMassAllCandidates", "hMassAllCandidates", kTH1D, {histAxes.axisXiMass});
      histos.add("reconstructedCandidates/Xi/hMassSelected", "hMassSelected", kTH1D, {histAxes.axisXiMass});
      histos.add("reconstructedCandidates/Omega/hMassAllCandidates", "hMassAllCandidates", kTH1D, {histAxes.axisOmegaMass});
      histos.add("reconstructedCandidates/Omega/hMassSelected", "hMassSelected", kTH1D, {histAxes.axisOmegaMass});

      histos.addClone("reconstructedCandidates/Xi/", "reconstructedCandidates/AntiXi/");
      histos.addClone("reconstructedCandidates/Omega/", "reconstructedCandidates/AntiOmega/");
    }

    if (selectionFlags.doQAforSelectionVariables) {
      if (!selectionFlags.applyDCADaughtersToPVSelection) {
        histos.add("reconstructedCandidates/K0/hDCAtoPVNegDaughter", "hDCAtoPVNegDaughter", kTH3D, {histAxes.axisK0Mass, histAxes.axisPt, histAxes.axisDCA});
        histos.add("reconstructedCandidates/K0/hDCAtoPVPosDaughter", "hDCAtoPVPosDaughter", kTH3D, {histAxes.axisK0Mass, histAxes.axisPt, histAxes.axisDCA});
        histos.add("reconstructedCandidates/Lambda/hDCAtoPVNegDaughter", "hDCAtoPVNegDaughter", kTH3D, {histAxes.axisLambdaMass, histAxes.axisPt, histAxes.axisDCA});
        histos.add("reconstructedCandidates/Lambda/hDCAtoPVPosDaughter", "hDCAtoPVPosDaughter", kTH3D, {histAxes.axisLambdaMass, histAxes.axisPt, histAxes.axisDCA});
      }
      if (!selectionFlags.applyV0RadiusSelection) {
        histos.add("reconstructedCandidates/K0/hV0Radius", "hV0Radius", kTH3D, {histAxes.axisK0Mass, histAxes.axisPt, histAxes.axisV0Radius});
        histos.add("reconstructedCandidates/Lambda/hV0Radius", "hV0Radius", kTH3D, {histAxes.axisLambdaMass, histAxes.axisPt, histAxes.axisV0Radius});
      }
      if (!selectionFlags.applyDCAdaughterSelection) {
        histos.add("reconstructedCandidates/K0/hDCAV0Daughters", "hDCAV0Daughters", kTH3D, {histAxes.axisK0Mass, histAxes.axisPt, histAxes.axisDCAV0Daughters});
        histos.add("reconstructedCandidates/Lambda/hDCAV0Daughters", "hDCAV0Daughters", kTH3D, {histAxes.axisLambdaMass, histAxes.axisPt, histAxes.axisDCAV0Daughters});
      }
      if (!selectionFlags.applyCosOfPAngleSelection) {
        histos.add("reconstructedCandidates/K0/hPointingAngle", "hPointingAngle", kTH3D, {histAxes.axisK0Mass, histAxes.axisPt, histAxes.axisPointingAngle});
        histos.add("reconstructedCandidates/Lambda/hPointingAngle", "hPointingAngle", kTH3D, {histAxes.axisLambdaMass, histAxes.axisPt, histAxes.axisPointingAngle});
      }
      if (!selectionFlags.applyLifetimeSelection) {
        histos.add("reconstructedCandidates/K0/hProperLifeTime", "hProperLifeTime", kTH3D, {histAxes.axisK0Mass, histAxes.axisPt, histAxes.axisProperLifeTime});
        histos.add("reconstructedCandidates/Lambda/hProperLifeTime", "hProperLifeTime", kTH3D, {histAxes.axisLambdaMass, histAxes.axisPt, histAxes.axisProperLifeTime});
      }
    }
    histos.addClone("reconstructedCandidates/Lambda/", "reconstructedCandidates/AntiLambda/");

    appliedSelectionCheckMask = 0;
    if (!selectionFlags.applyDCAdaughterSelection)
      SETBIT(appliedSelectionCheckMask, 0);
    if (!selectionFlags.applyCosOfPAngleSelection)
      SETBIT(appliedSelectionCheckMask, 1);
    if (!selectionFlags.applyDCADaughtersToPVSelection) {
      SETBIT(appliedSelectionCheckMask, 2);
      SETBIT(appliedSelectionCheckMask, 3);
    }
    if (!selectionFlags.applyV0RadiusSelection)
      SETBIT(appliedSelectionCheckMask, 4);
    if (!selectionFlags.applyLifetimeSelection)
      SETBIT(appliedSelectionCheckMask, 5);
  }
  void processAllFindableCandidates(aod::Collisions const& collisions, aod::McCollisions const&, aod::UpgradeV0s const& v0Recos, aod::UpgradeCascades const& cascRecos, Alice3Tracks const&)
  {
    for (const auto& collision : collisions) {
      float collisionZ = collision.posZ();
      // std::cout << "______ process V0_______" <<  collision.size() << std::endl;
      histos.fill(HIST("hPVz"), collisionZ);
    }

    for (const auto& v0Cand : v0Recos) {
      auto negV0Daughter = v0Cand.negTrack_as<Alice3Tracks>(); // de-reference neg track
      auto posV0Daughter = v0Cand.posTrack_as<Alice3Tracks>(); // de-reference pos track

      bool isK0 = v0Cand.mK0() > 0;
      if (isK0) {
        histos.fill(HIST("K0/hMassAllCandidates"), v0Cand.mK0(), v0Cand.pt());
        histos.fill(HIST("K0/hSelections"), 0); // all candidates
        histos.fill(HIST("K0/hDCANegDaughter"), negV0Daughter.dcaXY());
        histos.fill(HIST("K0/hDCAPosDaughter"), posV0Daughter.dcaXY());
        if (std::abs(negV0Daughter.dcaXY()) < selectionValues.dcaDaughtersToPVSelection)
          continue;
        histos.fill(HIST("K0/hSelections"), 1); // dcaXY cut
        if (std::abs(posV0Daughter.dcaXY()) < selectionValues.dcaDaughtersToPVSelection)
          continue;
        histos.fill(HIST("K0/hSelections"), 2); // dcaXY cut
        if (v0Cand.dcaV0Daughters() > selectionValues.dcaDaughterSelection)
          continue;
        histos.fill(HIST("K0/hSelections"), 3); // dca between daughters
        if (v0Cand.v0Radius() < selectionValues.v0RadiusSelection)
          continue;
        histos.fill(HIST("K0/hSelections"), 4); // radius cut
        if (std::abs(negV0Daughter.eta()) > selectionValues.etaDaughterSelection || std::abs(posV0Daughter.eta()) > selectionValues.etaDaughterSelection)
          continue;
        histos.fill(HIST("K0/hSelections"), 5); // eta cut
        histos.fill(HIST("K0/hMassSelected"), v0Cand.mK0(), v0Cand.pt());
      }
    }

    for (const auto& cascCand : cascRecos) {
      // auto bach = cascCand.bachTrack_as<Alice3Tracks>(); // de-reference bach track
      // auto neg = cascCand.negTrack_as<Alice3Tracks>();   // de-reference neg track
      // auto pos = cascCand.posTrack_as<Alice3Tracks>();   // de-reference pos track

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
    for (auto const& v0 : v0Candidates) {
      bool isK0 = (v0.mK0Short() - o2::constants::physics::MassK0Short) < selectionValues.acceptedK0MassWindow;
      bool isLambda = (v0.mLambda() - o2::constants::physics::MassLambda0) < selectionValues.acceptedLambdaMassWindow;
      bool isAntiLambda = (v0.mAntiLambda() - o2::constants::physics::MassLambda0) < selectionValues.acceptedLambdaMassWindow;

      histos.fill(HIST("reconstructedCandidates/hArmeterosBeforeAllSelections"), v0.alpha(), v0.qtArm());
      histos.fill(HIST("hV0CandidateCounter"), 0.5);
      if (selectionFlags.applyRapiditySelection) {
        if (isK0 && std::abs(v0.yK0Short()) > selectionValues.yK0Selection)
          continue;
        if ((isLambda || isAntiLambda) && std::abs(v0.yLambda()) > selectionValues.yLambdaSelection)
          continue;
      }
      histos.fill(HIST("hV0CandidateCounter"), 1.5);
      if (selectionFlags.applyDCAdaughterSelection) {
        if (std::abs(v0.dcaV0Daughters()) > selectionValues.dcaDaughterSelection)
          continue;
      } else {
        selectionCheck = v0.dcaV0Daughters();
      }
      histos.fill(HIST("hV0CandidateCounter"), 2.5);
      if (selectionFlags.applyCosOfPAngleSelection) {
        if (v0.cosPA() < selectionValues.cosPAngleSelection)
          continue;
      } else {
        selectionCheck = std::acos(v0.cosPA());
      }
      histos.fill(HIST("hV0CandidateCounter"), 3.5);
      if (selectionFlags.applyDCADaughtersToPVSelection) {
        if ((std::abs(v0.dcaNegToPV()) < selectionValues.dcaDaughtersToPVSelection) ||
            (std::abs(v0.dcaPosToPV()) < selectionValues.dcaDaughtersToPVSelection))
          continue;
      } else {
        selectionCheckPos = std::abs(v0.dcaPosToPV());
        selectionCheck = std::abs(v0.dcaNegToPV());
      }
      histos.fill(HIST("hV0CandidateCounter"), 4.5);
      if (selectionFlags.applyV0RadiusSelection) {
        if (v0.v0radius() < selectionValues.v0RadiusSelection)
          continue;
      } else {
        selectionCheck = v0.v0radius();
      }
      histos.fill(HIST("hV0CandidateCounter"), 5.5);
      if (isK0) {
        if (selectionFlags.applyArmenterosSelection) {
          if (v0.qtArm() < selectionValues.armenterosSelection * std::abs(v0.alpha()))
            continue;
        }
      }
      histos.fill(HIST("hV0CandidateCounter"), 6.5);
      if (isK0 && selectionFlags.applyCompetingMassRejection) {
        if (std::abs(v0.mLambda() - o2::constants::physics::MassLambda0) < selectionValues.competingMassRejectionK0)
          continue;
        if (std::abs(v0.mAntiLambda() - o2::constants::physics::MassLambda0) < selectionValues.competingMassRejectionK0)
          continue;
      }
      if ((isLambda || isAntiLambda) && selectionFlags.applyCompetingMassRejection) {
        if (std::abs(v0.mK0Short() - o2::constants::physics::MassK0Short) < selectionValues.competingMassRejectionLambda)
          continue;
      }
      histos.fill(HIST("hV0CandidateCounter"), 7.5);
      if (selectionFlags.applyLifetimeSelection) {
        if (isK0 && v0.distOverTotMom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassK0Short > selectionValues.lifetimecutak0)
          continue;
        if ((isLambda || isAntiLambda) && v0.distOverTotMom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0 > selectionValues.lifetimecutambda)
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
      if (selectionFlags.applyEtaDaughterSelection) {
        if (std::abs(posTrack.eta()) > selectionValues.etaDaughterSelection || std::abs(negTrack.eta()) > selectionValues.etaDaughterSelection)
          continue;
      }
      histos.fill(HIST("reconstructedCandidates/hEtaDaughters"), posTrack.eta());
      histos.fill(HIST("reconstructedCandidates/hEtaDaughters"), negTrack.eta());
      histos.fill(HIST("hV0CandidateCounter"), 9.5);

      histos.fill(HIST("reconstructedCandidates/hArmeterosAfterAllSelections"), v0.alpha(), v0.qtArm());
      if (selectionFlags.doQAforSelectionVariables) {
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

  void processFoundCascadeCandidates(aod::Collision const&, FullCascadeCandidates const& cascadeCandidates, Alice3Tracks const&, aod::McParticles const&)
  {
    for (const auto& casc : cascadeCandidates) {
      const bool isXi = (std::abs(casc.mXi() - o2::constants::physics::MassXiMinus) < selectionValues.acceptedXiMassWindow) && casc.sign() > 0;
      const bool isOm = (std::abs(casc.mOmega() - o2::constants::physics::MassOmegaMinus) < selectionValues.acceptedOmegaMassWindow) && casc.sign() > 0;
      const bool isAntiXi = (std::abs(casc.mXi() - o2::constants::physics::MassXiMinus) < selectionValues.acceptedXiMassWindow) && casc.sign() < 0;
      const bool isAntiOm = (std::abs(casc.mOmega() - o2::constants::physics::MassOmegaMinus) < selectionValues.acceptedOmegaMassWindow) && casc.sign() < 0;

      if (isXi) {
        histos.fill(HIST("reconstructedCandidates/Xi/hMassAllCandidates"), casc.mXi());
      }

      if (isOm) {
        histos.fill(HIST("reconstructedCandidates/Omega/hMassAllCandidates"), casc.mOmega());
      }

      if (isAntiXi) {
        histos.fill(HIST("reconstructedCandidates/AntiXi/hMassAllCandidates"), casc.mXi());
      }

      if (isAntiOm) {
        histos.fill(HIST("reconstructedCandidates/AntiOmega/hMassAllCandidates"), casc.mOmega());
      }

      // TODO Add selections
      if (isXi) {
        histos.fill(HIST("reconstructedCandidates/Xi/hMassSelected"), casc.mXi());
      }

      if (isOm) {
        histos.fill(HIST("reconstructedCandidates/Omega/hMassSelected"), casc.mOmega());
      }

      if (isAntiXi) {
        histos.fill(HIST("reconstructedCandidates/AntiXi/hMassSelected"), casc.mXi());
      }

      if (isAntiOm) {
        histos.fill(HIST("reconstructedCandidates/AntiOmega/hMassSelected"), casc.mOmega());
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
