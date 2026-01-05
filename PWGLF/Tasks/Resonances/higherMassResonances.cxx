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

/// \file higherMassResonances.cxx
/// \brief glueball resonance
/// \author Sawan <sawan.sawan@cern.ch>

// #include <TDatabasePDG.h>
#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h" //

#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h" //
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h" //
#include "Common/DataModel/PIDResponseTPC.h" //
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h" //
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/StepTHn.h"
#include "Framework/runDataProcessing.h" //
#include "ReconstructionDataFormats/Track.h"

#include "Math/GenVector/Boost.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "TF1.h"
#include "TRandom3.h"
#include <TDirectory.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <THn.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TPDGCode.h>
#include <TVector2.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <memory>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::aod::rctsel;
// using namespace o2::constants::physics;
using std::array;

struct HigherMassResonances {
  SliceCache cache;
  HistogramRegistry rEventSelection{"eventSelection", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rKzeroShort{"kzeroShort", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry hglue{"hglueball", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry hMChists{"hMChists", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  struct : ConfigurableGroup {
    Configurable<bool> requireRCTFlagChecker{"requireRCTFlagChecker", true, "Check event quality in run condition table"};
    Configurable<std::string> cfgEvtRCTFlagCheckerLabel{"cfgEvtRCTFlagCheckerLabel", "CBT_hadronPID", "Evt sel: RCT flag checker label"};
    Configurable<bool> cfgEvtRCTFlagCheckerZDCCheck{"cfgEvtRCTFlagCheckerZDCCheck", false, "Evt sel: RCT flag checker ZDC check"};
    Configurable<bool> cfgEvtRCTFlagCheckerLimitAcceptAsBad{"cfgEvtRCTFlagCheckerLimitAcceptAsBad", true, "Evt sel: RCT flag checker treat Limited Acceptance As Bad"};
  } rctCut;
  RCTFlagsChecker rctChecker;

  enum MultEstimator {
    kFT0M,
    kFT0A,
    kFT0C,
    kFV0A,
    kFV0C,
    kFV0M,
    kNEstimators
  };

  struct : ConfigurableGroup {
    // PID and QA
    Configurable<bool> qAv0{"qAv0", false, "qAv0"};
    Configurable<bool> qAPID{"qAPID", true, "qAPID"};
    Configurable<bool> qAevents{"qAevents", false, "QA of events"};
    Configurable<bool> qAcorrelation2Dhist{"qAcorrelation2Dhist", true, "Lamda K0 mass correlation"};
    Configurable<bool> qAOptimisation{"qAOptimisation", false, "QA for optimisation with multiple THnSparse Axes"};
    Configurable<bool> isApplyDCAv0topv{"isApplyDCAv0topv", false, "DCA V0 to PV"};
    Configurable<bool> hasTPC{"hasTPC", false, "TPC"};
    Configurable<bool> isselectTWOKsOnly{"isselectTWOKsOnly", true, "Select only events with two K0s"};
    Configurable<bool> isapplyPairRapidityMC{"isapplyPairRapidityMC", false, "Apply pair rapidity cut on reconstructed mother (after already applying rapidity cut on generated mother)"};
    Configurable<int> cSelectMultEstimator{"cSelectMultEstimator", 0, "Select multiplicity estimator: 0 - FT0M, 1 - FT0A, 2 - FT0C"};
    // Configurable<int> configOccCut{"configOccCut", 1000, "Occupancy cut"};
    // Configurable<bool> isVertexTOFMatched{"isVertexTOFMatched", false, "Vertex TOF Matched"};
    // Configurable<bool> isNoCollInTimeRangeStandard{"isNoCollInTimeRangeStandard", false, "No collision in time range standard"};
    // Configurable<bool> isSel8{"isSel8", false, "Event Selection 8"};

    // Configurables for event selection
    // Configurable<bool> isINELgt0{"isINELgt0", true, "INEL>0 selection"};
    Configurable<bool> isTriggerTVX{"isTriggerTVX", false, "TriggerTVX"};
    // Configurable<bool> isGoodZvtxFT0vsPV{"isGoodZvtxFT0vsPV", false, "IsGoodZvtxFT0vsPV"};
    // Configurable<bool> isApplyOccCut{"isApplyOccCut", true, "Apply occupancy cut"};
    Configurable<float> cutzvertex{"cutzvertex", 10.0f, "Accepted z-vertex range (cm)"};
    Configurable<bool> timFrameEvsel{"timFrameEvsel", true, "TPC Time frame boundary cut"};
    // Configurable<bool> isNoSameBunchPileup{"isNoSameBunchPileup", true, "kNoSameBunchPileup"};
    Configurable<bool> isAllLayersGoodITS{"isAllLayersGoodITS", true, "Require all ITS layers to be good"};
    Configurable<bool> isNoTimeFrameBorder{"isNoTimeFrameBorder", true, "kNoTimeFrameBorder"};
    Configurable<bool> isNoITSROFrameBorder{"isNoITSROFrameBorder", true, "kNoITSROFrameBorder"};

    // Configurable parameters for V0 selection
    Configurable<float> confV0PtMin{"confV0PtMin", 0.f, "Minimum transverse momentum of V0"};
    Configurable<float> confV0PtMax{"confV0PtMax", 100.f, "Maximum transverse momentum of V0"};
    Configurable<float> confPiPtMin{"confPiPtMin", 0.1f, "Minimum transverse momentum of pion daughter"};
    Configurable<float> confPiPtMax{"confPiPtMax", 100.f, "Maximum transverse momentum of pion daughter"};
    Configurable<float> cMaxDeltaM{"cMaxDeltaM", 0.01f, "Sqrt((m1-mPDG)^2 + (m2-mPDG)^2) < cMaxDeltaM)"};
    Configurable<float> confV0DCADaughMax{"confV0DCADaughMax", 1.0f, "DCA b/w V0 daughters"};
    Configurable<float> v0DCApostoPV{"v0DCApostoPV", 0.06, "DCA Pos To PV"};
    Configurable<float> v0DCAnegtoPV{"v0DCAnegtoPV", 0.06, "DCA Neg To PV"};
    Configurable<double> cMaxV0DCA{"cMaxV0DCA", 0.5, "DCA V0 to PV"};
    Configurable<float> confV0CPAMin{"confV0CPAMin", 0.97f, "Minimum CPA of V0"};
    Configurable<float> confV0TranRadV0Min{"confV0TranRadV0Min", 0.5f, "Minimum transverse radius"};
    // Configurable<float> confV0TranRadV0Max{"confV0TranRadV0Max", 200.f, "Maximum transverse radius"};
    Configurable<double> cMaxV0LifeTime{"cMaxV0LifeTime", 15, "Maximum V0 life time"};
    Configurable<double> cSigmaMassKs0{"cSigmaMassKs0", 4, "n Sigma cut on Ks0 mass (Mass (Ks) - cSigmaMassKs0*cWidthKs0)"};
    Configurable<double> cWidthKs0{"cWidthKs0", 0.005, "Width of KS0"};
    Configurable<float> confDaughEta{"confDaughEta", 0.8f, "V0 Daugh sel: max eta"};
    Configurable<float> confDaughTPCnclsMin{"confDaughTPCnclsMin", 70.f, "V0 Daugh sel: Min. nCls TPC"};
    Configurable<float> confDaughPIDCutTPC{"confDaughPIDCutTPC", 5, "PID selections for KS0 daughters"};
    Configurable<float> confDaughPIDCutTOF{"confDaughPIDCutTOF", 5, "PID selections for KS0 daughters in TOF"};
    Configurable<float> confKsrapidity{"confKsrapidity", 0.5f, "Rapidity cut on K0s"};
    // Configurable<bool> isStandardV0{"isStandardV0", false, "Standard V0 selection"};
    Configurable<bool> isApplyEtaCutK0s{"isApplyEtaCutK0s", false, "Apply eta cut on K0s daughters"};
    Configurable<float> cfgETAcut{"cfgETAcut", 0.8f, "Track ETA cut"};
    Configurable<float> deltaRDaugherCut{"deltaRDaugherCut", 0.001f, "DeltaR cut on V0 daughters"};
    Configurable<bool> deltaRK0sCut{"deltaRK0sCut", false, "Apply deltaR cut between two K0s"};

    // Configurable for track selection and multiplicity
    Configurable<float> cfgPTcut{"cfgPTcut", 0.2f, "Track PT cut"};
    Configurable<int> cfgNmixedEvents{"cfgNmixedEvents", 5, "Number of mixed events"};

    // Configurable for MC
    // Configurable<bool> isMC{"isMC", false, "Is MC"};
    Configurable<bool> isallGenCollisions{"isallGenCollisions", true, "To fill all generated collisions for the signal loss calculations"};
    Configurable<bool> isavoidsplitrackMC{"isavoidsplitrackMC", false, "avoid split track in MC"};
    Configurable<bool> isapplyRapidityMC{"isapplyRapidityMC", true, "Apply rapidity cut on generated and reconstructed particles"};
    Configurable<int> selectMCparticles{"selectMCparticles", 1, "0: f0(1710), 1: f2(1525), 2: a2(1320), 3: f0(1370), 4: f0(1500), 5: f2(1270)"};
    std::vector<int> pdgCodes = {10331, 335, 115, 10221, 9030221, 225};

    // output THnSparses
    Configurable<bool> activateHelicityFrame{"activateHelicityFrame", false, "Activate the THnSparse with cosThStar w.r.t. helicity axis"};
    Configurable<bool> activateCollinsSoperFrame{"activateCollinsSoperFrame", false, "Activate the THnSparse with cosThStar w.r.t. Collins soper axis"};
    Configurable<bool> activateProductionFrame{"activateProductionFrame", false, "Activate the THnSparse with cosThStar w.r.t. production axis"};
    Configurable<bool> activateBeamAxisFrame{"activateBeamAxisFrame", true, "Activate the THnSparse with cosThStar w.r.t. beam axis (Gottified jackson frame)"};
    Configurable<bool> activateRandomFrame{"activateRandomFrame", false, "Activate the THnSparse with cosThStar w.r.t. random axis"};
    Configurable<int> cRotations{"cRotations", 3, "Number of random rotations in the rotational background"};

    // Other cuts on Ks and glueball
    Configurable<bool> isapplyCompetingcut{"isapplyCompetingcut", false, "Competing cascade rejection cut"};
    Configurable<float> competingcascrejlambda{"competingcascrejlambda", 0.005, "rejecting competing cascade lambda"};
    Configurable<int> tpcCrossedrows{"tpcCrossedrows", 70, "TPC crossed rows"};
    Configurable<float> tpcCrossedrowsOverfcls{"tpcCrossedrowsOverfcls", 0.8, "TPC crossed rows over findable clusters"};
    Configurable<int> rotationalCut{"rotationalCut", 10, "Cut value (Rotation angle pi - pi/cut and pi + pi/cut)"};

    // // Mass and pT axis as configurables
    ConfigurableAxis binsCent{"binsCent", {VARIABLE_WIDTH, 0., 5., 10., 30., 50., 70., 100., 110., 150.}, "Binning of the centrality axis"};
    ConfigurableAxis configThnAxisPOL{"configThnAxisPOL", {20, -1.0, 1.0}, "Costheta axis"};
    ConfigurableAxis configThnAxisPhi{"configThnAxisPhi", {70, 0.0f, 7.0f}, "Phi axis"}; // 0 to 2pi
    ConfigurableAxis ksMassBins{"ksMassBins", {200, 0.45f, 0.55f}, "K0s invariant mass axis"};
    ConfigurableAxis cGlueMassBins{"cGlueMassBins", {200, 0.9f, 3.0f}, "Glueball invariant mass axis"};
    ConfigurableAxis cPtBins{"cPtBins", {200, 0.0f, 20.0f}, "Glueball pT axis"};
    ConfigurableAxis configAxisDeltaM{"configAxisDeltaM", {80, 0.0, 0.08}, "#it{M} (GeV/#it{c}^{2})"};
    ConfigurableAxis configAxisAngleSep{"configAxisAngleSep", {200, 0.0, 2.0}, "Angular separation between V0s"};
    ConfigurableAxis configAxisPtCorr{"configAxisPtCorr", {1000, 0.0, 100.0}, "Pt correlation between two K0s"};

    // fixed variables
    float rapidityMotherData = 0.5;
    float beamEnergy = 13600.0;
    double beamMomentum = std::sqrt(beamEnergy * beamEnergy / 4 - o2::constants::physics::MassProton * o2::constants::physics::MassProton); // GeV
    int noOfDaughters = 2;
  } config;

  // Service<o2::framework::O2DatabasePDG> PDGdatabase;
  TRandom* rn = new TRandom();

  // variables declaration
  float multiplicity = 0.0f;
  float theta2;
  ROOT::Math::PxPyPzMVector daughter1, daughter2, daughterRot, daughterRotCM, mother, motherRot, fourVecDauCM, fourVecDauCM1;
  ROOT::Math::PxPyPzEVector mother1;
  ROOT::Math::XYZVector randomVec, beamVec, normalVec;
  ROOT::Math::XYZVectorF v1CM, zaxisHE, yaxisHE, xaxisHE;
  // ROOT::Math::XYZVector threeVecDauCM, helicityVec, randomVec, beamVec, normalVec;
  ROOT::Math::XYZVector zBeam; // ẑ: beam direction in lab frame
  ROOT::Math::PxPyPzEVector beam1{0., 0., -config.beamMomentum, 13600. / 2.};
  ROOT::Math::PxPyPzEVector beam2{0., 0., config.beamMomentum, 13600. / 2.};
  ROOT::Math::XYZVectorF beam1CM, beam2CM, zAxisCS, yAxisCS, xAxisCS;

  // const double massK0s = o2::constants::physics::MassK0Short;
  bool isMix = false;

  void init(InitContext const&)
  {
    rctChecker.init(rctCut.cfgEvtRCTFlagCheckerLabel, rctCut.cfgEvtRCTFlagCheckerZDCCheck, rctCut.cfgEvtRCTFlagCheckerLimitAcceptAsBad);

    // Axes
    AxisSpec k0ShortMassAxis = {config.ksMassBins, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    AxisSpec glueballMassAxis = {config.cGlueMassBins, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    AxisSpec vertexZAxis = {60, -15.f, 15.f, "vrtx_{Z} [cm]"}; // for histogram
    AxisSpec ptAxis = {config.cPtBins, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec multiplicityAxis = {config.binsCent, "Multiplicity Axis"};
    AxisSpec thnAxisPOL{config.configThnAxisPOL, "Configurabel theta axis"};
    AxisSpec thnAxisPhi = {config.configThnAxisPhi, "Configurabel phi axis"}; // 0 to 2pi
    AxisSpec deltaMAxis = {config.configAxisDeltaM, "#Delta M  (GeV/#it{c}^{2})"};
    AxisSpec angleSepAxis = {config.configAxisAngleSep, "Angular separation between V0s"};
    AxisSpec ptCorrAxis = {config.configAxisPtCorr, "Pt correlation between two K0s"};

    //  THnSparses
    std::array<int, 5> sparses = {config.activateHelicityFrame, config.activateCollinsSoperFrame, config.activateProductionFrame, config.activateBeamAxisFrame, config.activateRandomFrame};

    if (std::accumulate(sparses.begin(), sparses.end(), 0) == 0) {
      LOGP(fatal, "No output THnSparses enabled");
    } else {
      if (config.activateHelicityFrame) {
        LOGP(info, "THnSparse with cosThStar w.r.t. helicity axis active.");
      }
      if (config.activateCollinsSoperFrame) {
        LOGP(info, "THnSparse with cosThStar w.r.t. Collins Soper axis active.");
      }
      if (config.activateProductionFrame) {
        LOGP(info, "THnSparse with cosThStar w.r.t. production axis active.");
      }
      if (config.activateBeamAxisFrame) {
        LOGP(info, "THnSparse with cosThStar w.r.t. beam axis active. (Gottified jackson frame)");
      }
      if (config.activateRandomFrame) {
        LOGP(info, "THnSparse with cosThStar w.r.t. random axis active.");
      }
    }

    // Event selection
    if (config.qAevents) {
      rEventSelection.add("hVertexZRec", "hVertexZRec", {HistType::kTH1F, {vertexZAxis}});
      rEventSelection.add("hmultiplicity", "multiplicity percentile distribution", {HistType::kTH1F, {{150, 0.0f, 150.0f}}});
      rEventSelection.add("htrackscheck_v0", "htrackscheck_v0", kTH1I, {{15, 0, 15}});
      rEventSelection.add("htrackscheck_v0_daughters", "htrackscheck_v0_daughters", kTH1I, {{15, 0, 15}});
      hMChists.add("events_check", "No. of events in the generated MC", kTH1I, {{20, 0, 20}});
      hMChists.add("events_checkrec", "No. of events in the reconstructed MC", kTH1I, {{20, 0, 20}});

      rEventSelection.add("hEventCut", "No. of event after cuts", kTH1I, {{20, 0, 20}});
      std::shared_ptr<TH1> hCutFlow = rEventSelection.get<TH1>(HIST("hEventCut"));
      hCutFlow->GetXaxis()->SetBinLabel(1, "All Events");
      hCutFlow->GetXaxis()->SetBinLabel(2, "|Vz| < cut");
      hCutFlow->GetXaxis()->SetBinLabel(3, "sel8");
      hCutFlow->GetXaxis()->SetBinLabel(4, "kNoTimeFrameBorder");
      hCutFlow->GetXaxis()->SetBinLabel(5, "kNoITSROFrameBorder");
      hCutFlow->GetXaxis()->SetBinLabel(6, "kNoSameBunchPileup");
      hCutFlow->GetXaxis()->SetBinLabel(7, "kIsGoodITSLayersAll");
      hCutFlow->GetXaxis()->SetBinLabel(8, "Occupancy Cut");
      hCutFlow->GetXaxis()->SetBinLabel(9, "rctChecker");
      hCutFlow->GetXaxis()->SetBinLabel(10, "kIsTriggerTVX");
      // hCutFlow->GetXaxis()->SetBinLabel(11, "kIsGoodZvtxFT0vsPV");
      // hCutFlow->GetXaxis()->SetBinLabel(12, "IsINELgt0");
      // hCutFlow->GetXaxis()->SetBinLabel(13, "isVertexITSTPC");
      // hCutFlow->GetXaxis()->SetBinLabel(14, "isVertexTOFMatched");

      std::shared_ptr<TH1> hv0label = rEventSelection.get<TH1>(HIST("htrackscheck_v0"));
      hv0label->GetXaxis()->SetBinLabel(1, "All Tracks");
      hv0label->GetXaxis()->SetBinLabel(2, "DCA V0 to PV");
      hv0label->GetXaxis()->SetBinLabel(3, "y K0s");
      hv0label->GetXaxis()->SetBinLabel(4, "V0 pT cut");
      hv0label->GetXaxis()->SetBinLabel(5, "Daughter DCA");
      hv0label->GetXaxis()->SetBinLabel(6, "CosPA");
      hv0label->GetXaxis()->SetBinLabel(7, "Decay Radius");
      hv0label->GetXaxis()->SetBinLabel(8, "Lifetime");
      hv0label->GetXaxis()->SetBinLabel(9, "CompetingCascade");
      hv0label->GetXaxis()->SetBinLabel(10, "Standard V0");
      hv0label->GetXaxis()->SetBinLabel(11, "Mass Tolerance");

      std::shared_ptr<TH1> hv0DauLabel = rEventSelection.get<TH1>(HIST("htrackscheck_v0_daughters"));
      hv0DauLabel->GetXaxis()->SetBinLabel(1, "AllDau Tracks");
      hv0DauLabel->GetXaxis()->SetBinLabel(2, "has TPC");
      hv0DauLabel->GetXaxis()->SetBinLabel(3, "TPC CrossedRows");
      hv0DauLabel->GetXaxis()->SetBinLabel(4, "TPC CRFC");
      hv0DauLabel->GetXaxis()->SetBinLabel(5, "TPC Chi2NCL");
      hv0DauLabel->GetXaxis()->SetBinLabel(6, "Charge");
      hv0DauLabel->GetXaxis()->SetBinLabel(7, "Charge");
      hv0DauLabel->GetXaxis()->SetBinLabel(8, "Eta");
      hv0DauLabel->GetXaxis()->SetBinLabel(9, "PID TPC");
      hv0DauLabel->GetXaxis()->SetBinLabel(10, "PID TOF");
      hv0DauLabel->GetXaxis()->SetBinLabel(11, "Pt cut");

      std::shared_ptr<TH1> hv0labelmcrec = hMChists.get<TH1>(HIST("events_checkrec"));
      hv0labelmcrec->GetXaxis()->SetBinLabel(1, "All Tracks");
      hv0labelmcrec->GetXaxis()->SetBinLabel(2, "V0Daughter Sel.");
      hv0labelmcrec->GetXaxis()->SetBinLabel(3, "V0 Sel.");
      hv0labelmcrec->GetXaxis()->SetBinLabel(4, "V0 PDG");
      hv0labelmcrec->GetXaxis()->SetBinLabel(5, "All Mothers");
      hv0labelmcrec->GetXaxis()->SetBinLabel(6, "Mother PDG");
      hv0labelmcrec->GetXaxis()->SetBinLabel(7, "Same Mother");
      hv0labelmcrec->GetXaxis()->SetBinLabel(8, "Split Track");
      hv0labelmcrec->GetXaxis()->SetBinLabel(9, "Global Index");
      hv0labelmcrec->GetXaxis()->SetBinLabel(10, "Generator");
      hv0labelmcrec->GetXaxis()->SetBinLabel(11, "Rapidity");

      std::shared_ptr<TH1> hv0labelmcgen = hMChists.get<TH1>(HIST("events_check"));
      hv0labelmcgen->GetXaxis()->SetBinLabel(1, "All Events");
      hv0labelmcgen->GetXaxis()->SetBinLabel(2, "Event Sel.");
      hv0labelmcgen->GetXaxis()->SetBinLabel(3, "Event reconstructed");
      hv0labelmcgen->GetXaxis()->SetBinLabel(4, "PDG check");
      hv0labelmcgen->GetXaxis()->SetBinLabel(5, "Rapidity");
      hv0labelmcgen->GetXaxis()->SetBinLabel(6, "Daughters2");
      hv0labelmcgen->GetXaxis()->SetBinLabel(7, "PhysicalPrimary");
      hv0labelmcgen->GetXaxis()->SetBinLabel(8, "Daughters K0s");
    }

    if (!config.qAOptimisation) {
      hglue.add("h3glueInvMassDS", "h3glueInvMassDS", kTHnSparseF, {multiplicityAxis, ptAxis, glueballMassAxis, thnAxisPOL, thnAxisPhi}, true);
      hglue.add("h3glueInvMassME", "h3glueInvMassME", kTHnSparseF, {multiplicityAxis, ptAxis, glueballMassAxis, thnAxisPOL, thnAxisPhi}, true);
      hglue.add("h3glueInvMassRot", "h3glueInvMassRot", kTHnSparseF, {multiplicityAxis, ptAxis, glueballMassAxis, thnAxisPOL, thnAxisPhi}, true);
    } else {
      hglue.add("h3glueInvMassDS", "h3glueInvMassDS", kTHnSparseF, {multiplicityAxis, ptAxis, glueballMassAxis, deltaMAxis, angleSepAxis, ptCorrAxis}, true);
      hglue.add("h3glueInvMassME", "h3glueInvMassME", kTHnSparseF, {multiplicityAxis, ptAxis, glueballMassAxis, deltaMAxis, angleSepAxis, ptCorrAxis}, true);
      hglue.add("h3glueInvMassRot", "h3glueInvMassRot", kTHnSparseF, {multiplicityAxis, ptAxis, glueballMassAxis, deltaMAxis, angleSepAxis, ptCorrAxis}, true);
    }

    // K0s topological/PID cuts
    if (config.qAcorrelation2Dhist) {
      rKzeroShort.add("mass_lambda_kshort_before", "mass under lambda hypotheses and Kshort mass", kTH2F, {{100, 0.2, 0.8}, {100, 0.9, 1.5}});
      rKzeroShort.add("mass_lambda_kshort_after10", "mass under lambda hypotheses and Kshort mass", kTH2F, {{100, 0.2, 0.8}, {100, 0.9, 1.5}});
    }
    if (config.qAv0) {
      // Invariant Mass
      rKzeroShort.add("hMassK0Shortbefore", "hMassK0Shortbefore", kTHnSparseF, {k0ShortMassAxis, ptAxis});
      rKzeroShort.add("hK0ShortMassCorr", "hK0ShortMassCorr", kTHnSparseF, {k0ShortMassAxis, k0ShortMassAxis, deltaMAxis});
      // rKzeroShort.add("hK0ShortMassCorrAfterCut", "hK0ShortMassCorrAfterCut", kTH2F, {k0ShortMassAxis, k0ShortMassAxis});
      rKzeroShort.add("hK0sPtCorrelation", "hK0sPtCorrelation", kTH1F, {{1000, 0.0f, 100.0f}});
      rKzeroShort.add("hMassK0ShortSelected", "hMassK0ShortSelected", kTHnSparseF, {k0ShortMassAxis, ptAxis});
      // Topological histograms (after the selection)
      rKzeroShort.add("hDCAV0Daughters", "DCA between v0 daughters", {HistType::kTH1F, {{60, -3.0f, 3.0f}}});
      rKzeroShort.add("hV0CosPA", "hV0CosPA", {HistType::kTH1F, {{100, 0.96f, 1.1f}}});
      rKzeroShort.add("hLT", "hLT", {HistType::kTH1F, {{100, 0.0f, 50.0f}}});
      rKzeroShort.add("angularSeparation", "Angular distribution between two K0s vs pT", {HistType::kTH1F, {{200, 0.0f, 4.0f}}});
      rKzeroShort.add("hDauDeltaR", "Delta R of positive and negative daughers", {HistType::kTHnSparseF, {angleSepAxis, angleSepAxis}});
    }
    rKzeroShort.add("NksProduced", "Number of K0s produced", kTH1I, {{15, -0.5, 14.5}});

    if (config.qAPID) {
      rKzeroShort.add("hNSigmaPosPionK0s_before", "hNSigmaPosPionK0s_before", {HistType::kTH2F, {{ptAxis}, {100, -5.f, 5.f}}});
      rKzeroShort.add("hNSigmaPosPionK0s_after", "hNSigmaPosPionK0s_after", {HistType::kTH2F, {{ptAxis}, {100, -5.f, 5.f}}});
      rKzeroShort.add("hNSigmaNegPionK0s_before", "hNSigmaNegPionK0s_before", {HistType::kTH2F, {{ptAxis}, {100, -5.f, 5.f}}});
      rKzeroShort.add("hNSigmaNegPionK0s_after", "hNSigmaNegPionK0s_after", {HistType::kTH2F, {{ptAxis}, {100, -5.f, 5.f}}});
      // rKzeroShort.add("dE_by_dx_TPC", "dE/dx signal in the TPC as a function of pT", kTH2F, {config.axisPtfordEbydx, config.axisdEdx});
    }

    // For MC
    if (doprocessGen || doprocessRec) {
      hMChists.add("Genf1710", "Gen f_{0}(1710)", kTHnSparseF, {multiplicityAxis, ptAxis, thnAxisPOL});
      hMChists.add("Genf17102", "Gen f_{0}(1710)", kTHnSparseF, {multiplicityAxis, ptAxis, thnAxisPOL});
      hMChists.add("Recf1710_pt1", "Rec f_{0}(1710) p_{T}", kTHnSparseF, {multiplicityAxis, ptAxis, glueballMassAxis, thnAxisPOL});
      hMChists.add("Recf1710_pt2", "Rec f_{0}(1710) p_{T}", kTHnSparseF, {multiplicityAxis, ptAxis, glueballMassAxis, thnAxisPOL});
      hMChists.add("h1Recsplit", "Rec p_{T}2", kTH1F, {ptAxis});
      hMChists.add("Genf1710_mass", "Gen f_{0}(1710) mass", kTH1F, {glueballMassAxis});
      hMChists.add("Genf1710_mass2", "Gen f_{0}(1710) mass", kTH1F, {glueballMassAxis});
      hMChists.add("GenPhi", "Gen Phi", kTH1F, {{70, 0.0, 7.0f}});
      hMChists.add("GenPhi2", "Gen Phi", kTH1F, {{70, 0.0, 7.0f}});
      hMChists.add("GenEta", "Gen Eta", kTHnSparseF, {{150, -1.5f, 1.5f}});
      hMChists.add("GenEta2", "Gen Eta", kTHnSparseF, {{150, -1.5f, 1.5f}});
      hMChists.add("GenRapidity", "Gen Rapidity", kTHnSparseF, {{100, -1.0f, 1.0f}});
      hMChists.add("GenRapidity2", "Gen Rapidity", kTHnSparseF, {{100, -1.0f, 1.0f}});
      hMChists.add("RecEta", "Rec Eta", kTH1F, {{150, -1.5f, 1.5f}});
      hMChists.add("RecEta2", "Rec Eta", kTH1F, {{150, -1.5f, 1.5f}});
      hMChists.add("RecPhi", "Rec Phi", kTH1F, {{70, 0.0f, 7.0f}});
      hMChists.add("RecPhi2", "Rec Phi", kTH1F, {{70, 0.0f, 7.0f}});
      hMChists.add("RecRapidity", "Rec Rapidity", kTH1F, {{100, -1.0f, 1.0f}});
      hMChists.add("RecRapidity2", "Rec Rapidity", kTH1F, {{100, -1.0f, 1.0f}});
      hMChists.add("Rec_Multiplicity", "Multiplicity in MC", kTH1F, {multiplicityAxis});
      hMChists.add("MC_mult_after_event_sel", "Multiplicity in MC", kTH1F, {multiplicityAxis});
    }
  }

  template <typename Coll>
  bool selectionEvent(const Coll& collision, bool fillHist = true)
  {
    if (fillHist)
      rEventSelection.fill(HIST("hEventCut"), 0);

    if (std::abs(collision.posZ()) > config.cutzvertex)
      return false;
    if (fillHist)
      rEventSelection.fill(HIST("hEventCut"), 1);

    // if (config.isSel8 && !collision.sel8())
    //   return false;
    if (fillHist)
      rEventSelection.fill(HIST("hEventCut"), 2);

    if (config.isNoTimeFrameBorder && !collision.selection_bit(aod::evsel::kNoTimeFrameBorder))
      return false;
    if (fillHist)
      rEventSelection.fill(HIST("hEventCut"), 3);

    if (config.isNoITSROFrameBorder && !collision.selection_bit(aod::evsel::kNoITSROFrameBorder))
      return false;
    if (fillHist)
      rEventSelection.fill(HIST("hEventCut"), 4);

    // if (config.isNoSameBunchPileup && (!collision.selection_bit(aod::evsel::kNoSameBunchPileup)))
    //   return false;
    if (fillHist)
      rEventSelection.fill(HIST("hEventCut"), 5);

    if (config.isAllLayersGoodITS && !collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll))
      return false;
    if (fillHist)
      rEventSelection.fill(HIST("hEventCut"), 6);

    // if (config.isNoCollInTimeRangeStandard && (!collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)))
    //   return false;

    // if (config.isApplyOccCut && (std::abs(collision.trackOccupancyInTimeRange()) > config.configOccCut))
    //   return false;
    if (fillHist)
      rEventSelection.fill(HIST("hEventCut"), 7);

    if (rctCut.requireRCTFlagChecker && !rctChecker(collision))
      return false;
    if (fillHist)
      rEventSelection.fill(HIST("hEventCut"), 8);

    if (config.isTriggerTVX && !collision.selection_bit(aod::evsel::kIsTriggerTVX))
      return false;
    if (fillHist)
      rEventSelection.fill(HIST("hEventCut"), 9);

    // if (config.isGoodZvtxFT0vsPV && !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV))
    //   return false;
    // if (fillHist)
    // rEventSelection.fill(HIST("hEventCut"), 10);

    // if (config.isINELgt0 && !collision.isInelGt0()) {
    //   return false;
    // }
    // if (fillHist)
    // rEventSelection.fill(HIST("hEventCut"), 11);

    // if (config.isVertexITSTPC && !collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
    //   return false;
    // }
    // if (fillHist)
    //   rEventSelection.fill(HIST("hEventCut"), 12);

    // if (config.isVertexTOFMatched && !collision.selection_bit(aod::evsel::kIsVertexTOFmatched)) {
    //   return false;
    // }
    // if (fillHist)
    //   rEventSelection.fill(HIST("hEventCut"), 13);

    return true;
  }

  template <typename Collision, typename V0>
  bool selectionV0(Collision const& collision, V0 const& candidate, float /*multiplicity*/)
  {
    // const float qtarm = candidate.qtarm();
    // const float alph = candidate.alpha();
    // float arm = qtarm / alph;
    const float pT = candidate.pt();
    const float tranRad = candidate.v0radius();
    const float dcaDaughv0 = candidate.dcaV0daughters();
    const float cpav0 = candidate.v0cosPA();

    float ctauK0s = candidate.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassK0Short;
    float lowmasscutks0 = o2::constants::physics::MassKPlus - config.cWidthKs0 * config.cSigmaMassKs0;
    float highmasscutks0 = o2::constants::physics::MassKPlus + config.cWidthKs0 * config.cSigmaMassKs0;

    if (config.qAv0) {
      rKzeroShort.fill(HIST("hMassK0Shortbefore"), candidate.mK0Short(), candidate.pt());
      rKzeroShort.fill(HIST("hLT"), ctauK0s);
      rKzeroShort.fill(HIST("hDCAV0Daughters"), candidate.dcaV0daughters());
      rKzeroShort.fill(HIST("hV0CosPA"), candidate.v0cosPA());
    }
    if (config.qAcorrelation2Dhist)
      rKzeroShort.fill(HIST("mass_lambda_kshort_before"), candidate.mK0Short(), candidate.mLambda());

    rEventSelection.fill(HIST("htrackscheck_v0"), 0.5);

    if (config.isApplyDCAv0topv && std::fabs(candidate.dcav0topv()) > config.cMaxV0DCA) {
      return false;
    }
    rEventSelection.fill(HIST("htrackscheck_v0"), 1.5);

    if (std::abs(candidate.rapidity(0)) >= config.confKsrapidity) {
      return false;
    }
    rEventSelection.fill(HIST("htrackscheck_v0"), 2.5);

    if (pT < config.confV0PtMin || pT > config.confV0PtMax) {
      return false;
    }
    rEventSelection.fill(HIST("htrackscheck_v0"), 3.5);

    if (dcaDaughv0 > config.confV0DCADaughMax) {
      return false;
    }
    rEventSelection.fill(HIST("htrackscheck_v0"), 4.5);

    if (cpav0 < config.confV0CPAMin) {
      return false;
    }
    rEventSelection.fill(HIST("htrackscheck_v0"), 5.5);

    if (tranRad < config.confV0TranRadV0Min) {
      return false;
    }
    rEventSelection.fill(HIST("htrackscheck_v0"), 6.5);

    // if (tranRad > config.confV0TranRadV0Max) {
    //   return false;
    // }
    rEventSelection.fill(HIST("htrackscheck_v0"), 7.5);

    if (std::fabs(ctauK0s) > config.cMaxV0LifeTime) {
      return false;
    }
    rEventSelection.fill(HIST("htrackscheck_v0"), 8.5);

    if (config.isapplyCompetingcut && (std::abs(candidate.mLambda() - o2::constants::physics::MassLambda0) <= config.competingcascrejlambda || std::abs(candidate.mAntiLambda() - o2::constants::physics::MassLambda0) <= config.competingcascrejlambda)) {
      return false;
    }
    rEventSelection.fill(HIST("htrackscheck_v0"), 9.5);

    if (config.qAcorrelation2Dhist)
      rKzeroShort.fill(HIST("mass_lambda_kshort_after10"), candidate.mK0Short(), candidate.mLambda());

    if (config.qAv0) {
      rKzeroShort.fill(HIST("hMassK0ShortSelected"), candidate.mK0Short(), candidate.pt());
    }

    // if (config.isStandardV0 && candidate.v0Type() != 1) {
    //   return false; // Only standard V0s are selected
    // }
    rEventSelection.fill(HIST("htrackscheck_v0"), 10.5);

    if (candidate.mK0Short() < lowmasscutks0 || candidate.mK0Short() > highmasscutks0) {
      return false;
    }
    rEventSelection.fill(HIST("htrackscheck_v0"), 11.5);

    return true;
  }

  template <typename T, typename V0s>
  bool isSelectedV0Daughter(T const& track, float charge, double nsigmaV0DaughterTPC, V0s const& v0candidate)
  {
    if (config.qAPID) {
      // Filling the PID of the V0 daughters in the region of the K0 peak.
      (charge == 1) ? rKzeroShort.fill(HIST("hNSigmaPosPionK0s_before"), track.tpcInnerParam(), track.tpcNSigmaPi()) : rKzeroShort.fill(HIST("hNSigmaNegPionK0s_before"), track.tpcInnerParam(), track.tpcNSigmaPi());
    }
    const auto eta = track.eta();
    const auto tpcNClsF = track.tpcNClsFound();
    const auto sign = track.sign();

    rEventSelection.fill(HIST("htrackscheck_v0_daughters"), 0.5);

    if (config.hasTPC && !track.hasTPC())
      return false;
    rEventSelection.fill(HIST("htrackscheck_v0_daughters"), 1.5);

    if (track.tpcNClsCrossedRows() < config.tpcCrossedrows)
      return false;
    rEventSelection.fill(HIST("htrackscheck_v0_daughters"), 2.5);

    if (track.tpcCrossedRowsOverFindableCls() < config.tpcCrossedrowsOverfcls)
      return false;
    rEventSelection.fill(HIST("htrackscheck_v0_daughters"), 3.5);

    if (tpcNClsF < config.confDaughTPCnclsMin) {
      return false;
    }
    rEventSelection.fill(HIST("htrackscheck_v0_daughters"), 4.5);

    if (charge < 0 && sign > 0) {
      return false;
    }
    rEventSelection.fill(HIST("htrackscheck_v0_daughters"), 5.5);

    if (charge > 0 && sign < 0) {
      return false;
    }
    rEventSelection.fill(HIST("htrackscheck_v0_daughters"), 6.5);

    if (std::abs(eta) > config.confDaughEta) {
      return false;
    }
    rEventSelection.fill(HIST("htrackscheck_v0_daughters"), 7.5);

    if (std::abs(nsigmaV0DaughterTPC) > config.confDaughPIDCutTPC) {
      return false;
    }
    rEventSelection.fill(HIST("htrackscheck_v0_daughters"), 8.5);

    if (std::abs(v0candidate.tofNSigmaK0PiPlus()) > config.confDaughPIDCutTOF && v0candidate.positiveHasTOF()) {
      return false;
    }

    if (std::abs(v0candidate.tofNSigmaK0PiMinus()) > config.confDaughPIDCutTOF && v0candidate.negativeHasTOF()) {
      return false;
    }
    rEventSelection.fill(HIST("htrackscheck_v0_daughters"), 9.5);

    if (track.pt() < config.confPiPtMin || track.pt() > config.confPiPtMax) {
      return false;
    }
    rEventSelection.fill(HIST("htrackscheck_v0_daughters"), 10.5);

    if (config.qAPID) {
      (charge == 1) ? rKzeroShort.fill(HIST("hNSigmaPosPionK0s_after"), track.tpcInnerParam(), track.tpcNSigmaPi()) : rKzeroShort.fill(HIST("hNSigmaNegPionK0s_after"), track.tpcInnerParam(), track.tpcNSigmaPi());
    }

    return true;
  }

  double deltaM(double m1, double m2)
  {
    const double d1 = m1 - o2::constants::physics::MassK0Short;
    const double d2 = m2 - o2::constants::physics::MassK0Short;
    return std::sqrt(d1 * d1 + d2 * d2);
  }

  using EventCandidatesDerivedData = soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraStamps>;
  using V0CandidatesDerivedData = soa::Join<aod::V0CollRefs, aod::V0Cores, aod::V0Extras, aod::V0TOFPIDs, aod::V0TOFNSigmas>;
  // using DauTracks = soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs, aod::DauTrackTOFPIDs>;
  using DauTracks = soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>;

  template <typename TV0>
  bool isSelectedK0sDaughtersDerived(TV0 const& v0)
  {
    // Fpr derived dataset

    // de-ref track extras
    auto posTrackExtra = v0.template posTrackExtra_as<DauTracks>();
    auto negTrackExtra = v0.template negTrackExtra_as<DauTracks>();

    if (std::abs(v0.positiveeta()) > config.confDaughEta || std::abs(v0.negativeeta()) > config.confDaughEta) {
      return false;
    }

    if (posTrackExtra.tpcNClsCrossedRows() < config.tpcCrossedrows || negTrackExtra.tpcNClsCrossedRows() < config.tpcCrossedrows) {
      return false;
    }

    if (posTrackExtra.tpcNClsFound() < config.confDaughTPCnclsMin || negTrackExtra.tpcNClsFound() < config.confDaughTPCnclsMin) {
      return false;
    }
    if (posTrackExtra.tpcCrossedRowsOverFindableCls() < config.tpcCrossedrowsOverfcls || negTrackExtra.tpcCrossedRowsOverFindableCls() < config.tpcCrossedrowsOverfcls) {
      return false;
    }

    // check TPC PID
    if (((std::abs(posTrackExtra.tpcNSigmaPi()) > config.confDaughPIDCutTPC) || (std::abs(negTrackExtra.tpcNSigmaPi()) > config.confDaughPIDCutTPC))) {
      return false;
    }

    // // check TOF PID if TOF exists

    if (config.isApplyDCAv0topv && (std::abs(v0.dcapostopv()) < config.cMaxV0DCA || std::abs(v0.dcanegtopv()) < config.cMaxV0DCA)) {
      return false;
    }

    if (std::abs(v0.tofNSigmaK0PiPlus()) > config.confDaughPIDCutTOF && v0.positiveHasTOF()) {
      return false;
    }

    if (std::abs(v0.tofNSigmaK0PiMinus()) > config.confDaughPIDCutTOF && v0.negativeHasTOF()) {
      return false;
    }

    double deltaRDaugherPos = std::sqrt(TVector2::Phi_mpi_pi(v0.positivephi() - v0.negativephi()) * TVector2::Phi_mpi_pi(v0.positivephi() - v0.negativephi()) + (v0.positiveeta() - v0.negativeeta()) * (v0.positiveeta() - v0.negativeeta()));
    double deltaRDaugherNeg = std::sqrt(TVector2::Phi_mpi_pi(v0.positivephi() - v0.negativephi()) * TVector2::Phi_mpi_pi(v0.positivephi() - v0.negativephi()) + (v0.positiveeta() - v0.negativeeta()) * (v0.positiveeta() - v0.negativeeta()));

    if (config.qAv0) {
      rKzeroShort.fill(HIST("hDauDeltaR"), deltaRDaugherPos, deltaRDaugherNeg);
    }

    if (deltaRDaugherPos < config.deltaRDaugherCut || deltaRDaugherNeg < config.deltaRDaugherCut) {
      return false;
    }

    return true;
  }

  // Defining filters for events (event selection)
  // Filter eventFilter = (o2::aod::evsel::sel8 == true);
  Filter posZFilter = (nabs(o2::aod::collision::posZ) < config.cutzvertex);
  Filter acceptenceFilter = (nabs(aod::track::eta) < config.cfgETAcut && nabs(aod::track::pt) > config.cfgPTcut);

  // Filters on V0s
  Filter preFilterV0 = (nabs(aod::v0data::dcapostopv) > config.v0DCApostoPV && nabs(aod::v0data::dcanegtopv) > config.v0DCAnegtoPV);

  // Defining the type of the daughter tracks
  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::FV0Mults, aod::MultZeqs, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::CentFV0As, aod::Mults, aod::PVMults>>;
  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA, aod::pidTPCFullPi, aod::pidTOFFullPi>>;
  using V0TrackCandidate = soa::Join<aod::V0Datas, aod::V0TOFPIDs, aod::V0TOFNSigmas>;
  // For Monte Carlo
  using EventCandidatesMC = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels, aod::CentFT0Cs, aod::CentFT0Ms, aod::CentFT0As, aod::CentFV0As, aod::PVMults>;
  using TrackCandidatesMC = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::McTrackLabels>>;
  using V0TrackCandidatesMC = soa::Filtered<soa::Join<aod::V0Datas, aod::V0TOFPIDs, aod::V0TOFNSigmas, aod::McV0Labels>>;
  // zBeam direction in lab frame

  template <typename T>
  void fillInvMass(const T& mother, float multiplicity, const T& daughter1, const T& daughter2, bool isMix)
  {

    // //polarization calculations
    // zBeam = ROOT::Math::XYZVector(0.f, 0.f, 1.f); // ẑ: beam direction in lab frame

    ROOT::Math::Boost boost{mother.BoostToCM()}; // define the boost to the center of mass frame
    fourVecDauCM = boost(daughter1);             // boost the frame of daughter to the center of mass frame
    // threeVecDauCM = fourVecDauCM.Vect();         // get the 3 vector of daughter in the frame of mother

    beam1CM = ROOT::Math::XYZVectorF((boost(beam1).Vect()).Unit());
    beam2CM = ROOT::Math::XYZVectorF((boost(beam2).Vect()).Unit());

    //========================Helicity and Production frame calculation==========================
    // define y = zBeam x z: Normal to the production plane
    // ẑ: mother direction in lab, boosted into mother's rest frame

    // auto motherLabDirection = ROOT::Math::XYZVector(0, 0, mother.Vect().Z()); // ẑ axis in lab frame

    // // ŷ = zBeam × ẑ
    // auto y_axis = zBeam.Cross(motherLabDirection).Unit();

    // // x̂ = ŷ × ẑ
    // auto x_axis = y_axis.Cross(motherLabDirection).Unit();

    // // Project daughter momentum onto x–y plane
    // auto p_proj_x = threeVecDauCM.Dot(x_axis);
    // auto p_proj_y = threeVecDauCM.Dot(y_axis);

    // // Calculate φ in [-π, π]
    // auto anglePhi = std::atan2(p_proj_y, p_proj_x); // φ in radians
    //=============================================================================================

    v1CM = ROOT::Math::XYZVectorF(boost(daughter1).Vect()).Unit();
    // ROOT::Math::XYZVectorF v2_CM{(boost(daughter1).Vect()).Unit()};
    // using positive sign convention for the first track
    // ROOT::Math::XYZVectorF v_CM = (t1.sign() > 0 ? v1CM : v2_CM); // here selected decay daughter momentum is intested. here you can choose one decay daughter no need to check both case as it is neutral particle for our case
    // Helicity Frame
    zaxisHE = ROOT::Math::XYZVectorF(mother.Vect()).Unit();
    yaxisHE = ROOT::Math::XYZVectorF(beam1CM.Cross(beam2CM)).Unit();
    xaxisHE = ROOT::Math::XYZVectorF(yaxisHE.Cross(zaxisHE)).Unit();

    // CosThetaHE = zaxisHE.Dot(v_CM);

    auto anglePhi = std::atan2(yaxisHE.Dot(v1CM), xaxisHE.Dot(v1CM));
    anglePhi = RecoDecay::constrainAngle(anglePhi, 0.0);
    // if (anglePhi < 0) {
    //   anglePhi += o2::constants::math::TwoPI; // ensure phi is in [0, 2pi]
    // }

    // CS Frame
    zAxisCS = ROOT::Math::XYZVectorF((beam1CM.Unit() - beam2CM.Unit())).Unit();
    yAxisCS = ROOT::Math::XYZVectorF(beam1CM.Cross(beam2CM)).Unit();
    xAxisCS = ROOT::Math::XYZVectorF(yAxisCS.Cross(zAxisCS)).Unit();
    double cosThetaStarCS = zAxisCS.Dot(v1CM);
    auto phiCS = std::atan2(yAxisCS.Dot(v1CM), xAxisCS.Dot(v1CM));
    phiCS = RecoDecay::constrainAngle(phiCS, 0.0);

    // if (std::abs(mother.Rapidity()) < config.rapidityMotherData) {
    if (config.activateHelicityFrame) {
      // helicityVec = mother.Vect(); // 3 vector of mother in COM frame
      // auto cosThetaStarHelicity = helicityVec.Dot(threeVecDauCM) / (std::sqrt(threeVecDauCM.Mag2()) * std::sqrt(helicityVec.Mag2()));
      auto cosThetaStarHelicity = mother.Vect().Dot(fourVecDauCM.Vect()) / (std::sqrt(fourVecDauCM.Vect().Mag2()) * std::sqrt(mother.Vect().Mag2()));
      if (!isMix) {
        if (std::abs(mother.Rapidity()) < config.rapidityMotherData) {
          hglue.fill(HIST("h3glueInvMassDS"), multiplicity, mother.Pt(), mother.M(), cosThetaStarHelicity, anglePhi);
        }

        for (int i = 0; i < config.cRotations; i++) {
          theta2 = rn->Uniform(o2::constants::math::PI - o2::constants::math::PI / config.rotationalCut, o2::constants::math::PI + o2::constants::math::PI / config.rotationalCut);

          daughterRot = ROOT::Math::PxPyPzMVector(daughter1.Px() * std::cos(theta2) - daughter1.Py() * std::sin(theta2), daughter1.Px() * std::sin(theta2) + daughter1.Py() * std::cos(theta2), daughter1.Pz(), daughter1.M());

          motherRot = daughterRot + daughter2;

          ROOT::Math::Boost boost2{motherRot.BoostToCM()};
          daughterRotCM = boost2(daughterRot);

          auto cosThetaStarHelicityRot = motherRot.Vect().Dot(daughterRotCM.Vect()) / (std::sqrt(daughterRotCM.Vect().Mag2()) * std::sqrt(motherRot.Vect().Mag2()));
          auto phiHelicityRot = std::atan2(yaxisHE.Dot(daughterRotCM.Vect().Unit()), xaxisHE.Dot(daughterRotCM.Vect().Unit()));
          phiHelicityRot = RecoDecay::constrainAngle(phiHelicityRot, 0.0);
          if (motherRot.Rapidity() < config.rapidityMotherData)
            hglue.fill(HIST("h3glueInvMassRot"), multiplicity, motherRot.Pt(), motherRot.M(), cosThetaStarHelicityRot, phiHelicityRot);
        }
      } else {
        if (std::abs(mother.Rapidity()) < config.rapidityMotherData) {
          hglue.fill(HIST("h3glueInvMassME"), multiplicity, mother.Pt(), mother.M(), cosThetaStarHelicity, anglePhi);
        }
      }
    } else if (config.activateCollinsSoperFrame) {
      if (!isMix) {
        if (std::abs(mother.Rapidity()) < config.rapidityMotherData) {
          hglue.fill(HIST("h3glueInvMassDS"), multiplicity, mother.Pt(), mother.M(), cosThetaStarCS, phiCS);
        }

        for (int i = 0; i < config.cRotations; i++) {
          theta2 = rn->Uniform(o2::constants::math::PI - o2::constants::math::PI / config.rotationalCut, o2::constants::math::PI + o2::constants::math::PI / config.rotationalCut);

          daughterRot = ROOT::Math::PxPyPzMVector(daughter1.Px() * std::cos(theta2) - daughter1.Py() * std::sin(theta2), daughter1.Px() * std::sin(theta2) + daughter1.Py() * std::cos(theta2), daughter1.Pz(), daughter1.M());

          motherRot = daughterRot + daughter2;

          ROOT::Math::Boost boost2{motherRot.BoostToCM()};
          daughterRotCM = boost2(daughterRot);

          auto cosThetaStarCSrot = zAxisCS.Dot(daughterRotCM.Vect()) / std::sqrt(daughterRotCM.Vect().Mag2());
          auto phiCSrot = std::atan2(yAxisCS.Dot(daughterRotCM.Vect().Unit()), xAxisCS.Dot(daughterRotCM.Vect().Unit()));
          phiCSrot = RecoDecay::constrainAngle(phiCSrot, 0.0);

          if (motherRot.Rapidity() < config.rapidityMotherData)
            hglue.fill(HIST("h3glueInvMassRot"), multiplicity, motherRot.Pt(), motherRot.M(), cosThetaStarCSrot, phiCSrot);
        }
      } else {
        if (std::abs(mother.Rapidity()) < config.rapidityMotherData) {
          hglue.fill(HIST("h3glueInvMassME"), multiplicity, mother.Pt(), mother.M(), cosThetaStarCS, phiCS);
        }
      }
    } else if (config.activateProductionFrame) {
      normalVec = ROOT::Math::XYZVector(mother.Py(), -mother.Px(), 0.f);
      auto cosThetaProduction = normalVec.Dot(fourVecDauCM.Vect()) / (std::sqrt(fourVecDauCM.Vect().Mag2()) * std::sqrt(normalVec.Mag2()));
      if (!isMix) {
        if (std::abs(mother.Rapidity()) < config.rapidityMotherData) {
          hglue.fill(HIST("h3glueInvMassDS"), multiplicity, mother.Pt(), mother.M(), cosThetaProduction, anglePhi);
        }
        for (int i = 0; i < config.cRotations; i++) {
          theta2 = rn->Uniform(o2::constants::math::PI - o2::constants::math::PI / config.rotationalCut, o2::constants::math::PI + o2::constants::math::PI / config.rotationalCut);
          motherRot = ROOT::Math::PxPyPzMVector(mother.Px() * std::cos(theta2) - mother.Py() * std::sin(theta2), mother.Px() * std::sin(theta2) + mother.Py() * std::cos(theta2), mother.Pz(), mother.M());
          if (std::abs(motherRot.Rapidity()) < config.rapidityMotherData) {
            hglue.fill(HIST("h3glueInvMassRot"), multiplicity, motherRot.Pt(), motherRot.M(), cosThetaProduction, anglePhi);
          }
        }
      } else {
        if (std::abs(mother.Rapidity()) < config.rapidityMotherData) {
          hglue.fill(HIST("h3glueInvMassME"), multiplicity, mother.Pt(), mother.M(), cosThetaProduction, anglePhi);
        }
      }
    } else if (config.activateBeamAxisFrame) {
      beamVec = ROOT::Math::XYZVector(0.f, 0.f, 1.f);
      auto cosThetaStarBeam = beamVec.Dot(fourVecDauCM.Vect()) / std::sqrt(fourVecDauCM.Vect().Mag2());
      if (!isMix) {
        if (std::abs(mother.Rapidity()) < config.rapidityMotherData) {
          hglue.fill(HIST("h3glueInvMassDS"), multiplicity, mother.Pt(), mother.M(), cosThetaStarBeam, anglePhi);
        }
        for (int i = 0; i < config.cRotations; i++) {
          theta2 = rn->Uniform(o2::constants::math::PI - o2::constants::math::PI / config.rotationalCut, o2::constants::math::PI + o2::constants::math::PI / config.rotationalCut);
          motherRot = ROOT::Math::PxPyPzMVector(mother.Px() * std::cos(theta2) - mother.Py() * std::sin(theta2), mother.Px() * std::sin(theta2) + mother.Py() * std::cos(theta2), mother.Pz(), mother.M());
          if (std::abs(motherRot.Rapidity()) < config.rapidityMotherData) {
            hglue.fill(HIST("h3glueInvMassRot"), multiplicity, motherRot.Pt(), motherRot.M(), cosThetaStarBeam, anglePhi);
          }
        }
      } else {
        if (std::abs(mother.Rapidity()) < config.rapidityMotherData) {
          hglue.fill(HIST("h3glueInvMassME"), multiplicity, mother.Pt(), mother.M(), cosThetaStarBeam, anglePhi);
        }
      }
    } else if (config.activateRandomFrame) {
      auto phiRandom = gRandom->Uniform(0.f, constants::math::TwoPI);
      auto thetaRandom = gRandom->Uniform(0.f, constants::math::PI);

      randomVec = ROOT::Math::XYZVector(std::sin(thetaRandom) * std::cos(phiRandom), std::sin(thetaRandom) * std::sin(phiRandom), std::cos(thetaRandom));
      auto cosThetaStarRandom = randomVec.Dot(fourVecDauCM.Vect()) / std::sqrt(fourVecDauCM.Vect().Mag2());
      if (!isMix) {
        if (std::abs(mother.Rapidity()) < config.rapidityMotherData) {
          hglue.fill(HIST("h3glueInvMassDS"), multiplicity, mother.Pt(), mother.M(), cosThetaStarRandom, phiRandom);
        }
        for (int i = 0; i < config.cRotations; i++) {
          theta2 = rn->Uniform(o2::constants::math::PI - o2::constants::math::PI / config.rotationalCut, o2::constants::math::PI + o2::constants::math::PI / config.rotationalCut);
          motherRot = ROOT::Math::PxPyPzMVector(mother.Px() * std::cos(theta2) - mother.Py() * std::sin(theta2), mother.Px() * std::sin(theta2) + mother.Py() * std::cos(theta2), mother.Pz(), mother.M());
          if (std::abs(motherRot.Rapidity()) < config.rapidityMotherData) {
            hglue.fill(HIST("h3glueInvMassRot"), multiplicity, motherRot.Pt(), motherRot.M(), cosThetaStarRandom, phiRandom);
          }
        }
      } else {
        if (std::abs(mother.Rapidity()) < config.rapidityMotherData) {
          hglue.fill(HIST("h3glueInvMassME"), multiplicity, mother.Pt(), mother.M(), cosThetaStarRandom, phiRandom);
        }
      }
    }
    // }
  }

  void processSE(EventCandidates::iterator const& collision, TrackCandidates const& /*tracks*/, V0TrackCandidate const& V0s)
  {
    multiplicity = 0.0;

    if (config.cSelectMultEstimator == kFT0M) {
      multiplicity = collision.centFT0M();
    } else if (config.cSelectMultEstimator == kFT0A) {
      multiplicity = collision.centFT0A();
    } else if (config.cSelectMultEstimator == kFT0C) {
      multiplicity = collision.centFT0C();
    } else if (config.cSelectMultEstimator == kFV0A) {
      multiplicity = collision.centFV0A();
    } else {
      multiplicity = collision.centFT0M(); // default
    }

    if (!selectionEvent(collision, true)) {
      return;
    }
    // if (rctCut.requireRCTFlagChecker && !rctCut.rctChecker(collision)) {
    //   return;
    // }

    // auto occupancyNumber = collision.trackOccupancyInTimeRange();
    // if (applyOccupancyCut && occupancyNumber < occupancyCut) {
    //   return;
    // }

    if (config.qAevents) {
      rEventSelection.fill(HIST("hVertexZRec"), collision.posZ());
      rEventSelection.fill(HIST("hmultiplicity"), multiplicity);
      // rEventSelection.fill(HIST("multdist_FT0M"), collision.multFT0M());
      // rEventSelection.fill(HIST("multdist_FT0A"), collision.multFT0A());
      // rEventSelection.fill(HIST("multdist_FT0C"), collision.multFT0C());
      // rEventSelection.fill(HIST("hNcontributor"), collision.numContrib());
    }

    std::vector<int> v0indexes;
    bool allConditionsMet = 0;

    for (const auto& [v1, v2] : combinations(CombinationsFullIndexPolicy(V0s, V0s))) {

      if (v1.size() == 0 || v2.size() == 0) {
        continue;
      }

      if (!selectionV0(collision, v1, multiplicity)) {
        continue;
      }
      if (!selectionV0(collision, v2, multiplicity)) {
        continue;
      }

      auto postrack1 = v1.template posTrack_as<TrackCandidates>();
      auto negtrack1 = v1.template negTrack_as<TrackCandidates>();
      auto postrack2 = v2.template posTrack_as<TrackCandidates>();
      auto negtrack2 = v2.template negTrack_as<TrackCandidates>();

      double nTPCSigmaPos1{postrack1.tpcNSigmaPi()};
      double nTPCSigmaNeg1{negtrack1.tpcNSigmaPi()};
      double nTPCSigmaPos2{postrack2.tpcNSigmaPi()};
      double nTPCSigmaNeg2{negtrack2.tpcNSigmaPi()};

      if (!(isSelectedV0Daughter(negtrack1, -1, nTPCSigmaNeg1, v1) && isSelectedV0Daughter(postrack1, 1, nTPCSigmaPos1, v1))) {
        continue;
      }
      if (!(isSelectedV0Daughter(postrack2, 1, nTPCSigmaPos2, v2) && isSelectedV0Daughter(negtrack2, -1, nTPCSigmaNeg2, v2))) {
        continue;
      }

      if (std::find(v0indexes.begin(), v0indexes.end(), v1.globalIndex()) == v0indexes.end()) {
        v0indexes.push_back(v1.globalIndex());
      }

      if (v2.globalIndex() <= v1.globalIndex()) {
        continue;
      }

      if (postrack1.globalIndex() == postrack2.globalIndex()) {
        continue;
      }
      if (negtrack1.globalIndex() == negtrack2.globalIndex()) {
        continue;
      }

      double deltaRDaugherPos = std::sqrt(TVector2::Phi_mpi_pi(postrack1.phi() - negtrack1.phi()) * TVector2::Phi_mpi_pi(postrack1.phi() - negtrack1.phi()) + (postrack1.eta() - negtrack1.eta()) * (postrack1.eta() - negtrack1.eta()));
      double deltaRDaugherNeg = std::sqrt(TVector2::Phi_mpi_pi(postrack2.phi() - negtrack2.phi()) * TVector2::Phi_mpi_pi(postrack2.phi() - negtrack2.phi()) + (postrack2.eta() - negtrack2.eta()) * (postrack2.eta() - negtrack2.eta()));

      if (config.qAv0) {
        rKzeroShort.fill(HIST("hDauDeltaR"), deltaRDaugherPos, deltaRDaugherNeg);
      }

      if (deltaRDaugherPos < config.deltaRDaugherCut || deltaRDaugherNeg < config.deltaRDaugherCut) {
        continue;
      }

      if (config.isApplyEtaCutK0s && (v1.eta() < config.confDaughEta || v2.eta() < config.confDaughEta)) {
        continue;
      }

      allConditionsMet = 1;
      daughter1 = ROOT::Math::PxPyPzMVector(v1.px(), v1.py(), v1.pz(), o2::constants::physics::MassK0Short); // Kshort
      daughter2 = ROOT::Math::PxPyPzMVector(v2.px(), v2.py(), v2.pz(), o2::constants::physics::MassK0Short); // Kshort

      mother = daughter1 + daughter2; // invariant mass of Kshort pair
      isMix = false;

      const double deltaMass = deltaM(v1.mK0Short(), v2.mK0Short());
      if (config.qAv0) {
        rKzeroShort.fill(HIST("hK0ShortMassCorr"), v1.mK0Short(), v2.mK0Short(), deltaMass);
      }

      if (!config.qAOptimisation) {
        if (deltaMass > config.cMaxDeltaM) {
          continue;
        }
      }

      // if (config.qAv0) {
      //   rKzeroShort.fill(HIST("hK0ShortMassCorrAfterCut"), v1.mK0Short(), v2.mK0Short());
      // }

      const double ptCorr = (mother.Pt() - daughter1.Pt() != 0.) ? daughter1.Pt() / (mother.Pt() - daughter1.Pt()) : 0.;
      if (config.qAv0) {
        rKzeroShort.fill(HIST("hK0sPtCorrelation"), ptCorr);
      }

      double deltaRvalue = std::sqrt(TVector2::Phi_mpi_pi(v1.phi() - v2.phi()) * TVector2::Phi_mpi_pi(v1.phi() - v2.phi()) + (v1.eta() - v2.eta()) * (v1.eta() - v2.eta()));

      if (!config.qAOptimisation) {
        if (deltaRvalue < config.deltaRK0sCut) {
          continue;
        }
      }

      if (!config.isselectTWOKsOnly && !config.qAOptimisation) {
        fillInvMass(mother, multiplicity, daughter1, daughter2, isMix);
      }

      if (!config.isselectTWOKsOnly && config.qAOptimisation) {

        if (std::abs(mother.Rapidity()) < config.rapidityMotherData) {
          hglue.fill(HIST("h3glueInvMassDS"), multiplicity, mother.Pt(), mother.M(), deltaMass, deltaRvalue, ptCorr);
        }

        for (int i = 0; i < config.cRotations; i++) {
          double theta2 = rn->Uniform(o2::constants::math::PI - o2::constants::math::PI / config.rotationalCut, o2::constants::math::PI + o2::constants::math::PI / config.rotationalCut);

          daughterRot = ROOT::Math::PxPyPzMVector(daughter1.Px() * std::cos(theta2) - daughter1.Py() * std::sin(theta2), daughter1.Px() * std::sin(theta2) + daughter1.Py() * std::cos(theta2), daughter1.Pz(), daughter1.M());

          motherRot = daughterRot + daughter2;

          double pTcorrRot = std::abs(daughterRot.Pt() + daughter2.Pt()) / motherRot.Pt();

          if (motherRot.Rapidity() < config.rapidityMotherData)
            hglue.fill(HIST("h3glueInvMassRot"), multiplicity, motherRot.Pt(), motherRot.M(), deltaMass, deltaRvalue, pTcorrRot);
        }
      }
    }
    int sizeofv0indexes = v0indexes.size();
    rKzeroShort.fill(HIST("NksProduced"), sizeofv0indexes);
    if (config.isselectTWOKsOnly && sizeofv0indexes == config.noOfDaughters && allConditionsMet) {
      fillInvMass(mother, multiplicity, daughter1, daughter2, false);
    }
    v0indexes.clear();
  }
  PROCESS_SWITCH(HigherMassResonances, processSE, "same event process", true);

  ConfigurableAxis axisVertex{"axisVertex", {20, -10, 10}, "vertex axis for ME mixing"};
  // ConfigurableAxis axisMultiplicityClass{"axisMultiplicityClass", {10, 0, 100}, "multiplicity percentile for ME mixing"};
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {2000, 0, 10000}, "TPC multiplicity axis for ME mixing"};

  // using BinningTypeTPCMultiplicity = ColumnBinningPolicy<aod::collision::PosZ, aod::mult::MultTPC>;
  using BinningTypeFT0M = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0M>;
  using BinningTypeFT0A = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0A>;
  using BinningTypeFT0C = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C>;
  using BinningTypeFV0A = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFV0A>;

  BinningTypeFT0M binningOnFT0M{{axisVertex, axisMultiplicity}, true};
  BinningTypeFT0A binningOnFT0A{{axisVertex, axisMultiplicity}, true};
  BinningTypeFT0C binningOnFT0C{{axisVertex, axisMultiplicity}, true};
  BinningTypeFV0A binningOnFV0A{{axisVertex, axisMultiplicity}, true};

  using BinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0M>;
  BinningType colBinning{{axisVertex, axisMultiplicity}, true};                                   // for derived data only
  Preslice<V0CandidatesDerivedData> tracksPerCollisionV0Mixed = o2::aod::v0data::straCollisionId; // for derived data only

  void processME(EventCandidates const& collisions, TrackCandidates const& /*tracks*/, V0TrackCandidate const& v0s)
  {
    auto tracksTuple = std::make_tuple(v0s);
    SameKindPair<EventCandidates, V0TrackCandidate, BinningTypeFT0M> pair1{binningOnFT0M, config.cfgNmixedEvents, -1, collisions, tracksTuple, &cache};
    SameKindPair<EventCandidates, V0TrackCandidate, BinningTypeFT0A> pair2{binningOnFT0A, config.cfgNmixedEvents, -1, collisions, tracksTuple, &cache};
    SameKindPair<EventCandidates, V0TrackCandidate, BinningTypeFT0C> pair3{binningOnFT0C, config.cfgNmixedEvents, -1, collisions, tracksTuple, &cache};
    SameKindPair<EventCandidates, V0TrackCandidate, BinningTypeFV0A> pair4{binningOnFV0A, config.cfgNmixedEvents, -1, collisions, tracksTuple, &cache};

    auto runMixing = [&](auto& pair, auto multiplicityGetter) {
      for (const auto& [c1, tracks1, c2, tracks2] : pair) {

        multiplicity = multiplicityGetter(c1);

        if (!selectionEvent(c1, false) || !selectionEvent(c2, false)) {
          continue;
        }

        for (const auto& [t1, t2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {

          if (t1.size() == 0 || t2.size() == 0) {
            continue;
          }

          if (!selectionV0(c1, t1, multiplicity))
            continue;
          if (!selectionV0(c2, t2, multiplicity))
            continue;

          auto postrack1 = t1.template posTrack_as<TrackCandidates>();
          auto negtrack1 = t1.template negTrack_as<TrackCandidates>();
          auto postrack2 = t2.template posTrack_as<TrackCandidates>();
          auto negtrack2 = t2.template negTrack_as<TrackCandidates>();
          if (postrack1.globalIndex() == postrack2.globalIndex()) {
            continue;
          }
          if (negtrack1.globalIndex() == negtrack2.globalIndex()) {
            continue;
          }
          double nTPCSigmaPos1{postrack1.tpcNSigmaPi()};
          double nTPCSigmaNeg1{negtrack1.tpcNSigmaPi()};
          double nTPCSigmaPos2{postrack2.tpcNSigmaPi()};
          double nTPCSigmaNeg2{negtrack2.tpcNSigmaPi()};

          if (!isSelectedV0Daughter(postrack1, 1, nTPCSigmaPos1, t1)) {
            continue;
          }
          if (!isSelectedV0Daughter(postrack2, 1, nTPCSigmaPos2, t2)) {
            continue;
          }
          if (!isSelectedV0Daughter(negtrack1, -1, nTPCSigmaNeg1, t1)) {
            continue;
          }
          if (!isSelectedV0Daughter(negtrack2, -1, nTPCSigmaNeg2, t2)) {
            continue;
          }

          double deltaRDaugherPos = std::sqrt(TVector2::Phi_mpi_pi(postrack1.phi() - negtrack1.phi()) * TVector2::Phi_mpi_pi(postrack1.phi() - negtrack1.phi()) + (postrack1.eta() - negtrack1.eta()) * (postrack1.eta() - negtrack1.eta()));
          double deltaRDaugherNeg = std::sqrt(TVector2::Phi_mpi_pi(postrack2.phi() - negtrack2.phi()) * TVector2::Phi_mpi_pi(postrack2.phi() - negtrack2.phi()) + (postrack2.eta() - negtrack2.eta()) * (postrack2.eta() - negtrack2.eta()));

          if (deltaRDaugherPos < config.deltaRDaugherCut || deltaRDaugherNeg < config.deltaRDaugherCut) {
            continue;
          }

          if (config.isApplyEtaCutK0s && (t1.eta() < config.confDaughEta || t2.eta() < config.confDaughEta)) {
            continue;
          }

          daughter1 = ROOT::Math::PxPyPzMVector(t1.px(), t1.py(), t1.pz(), o2::constants::physics::MassK0Short); // Kshort
          daughter2 = ROOT::Math::PxPyPzMVector(t2.px(), t2.py(), t2.pz(), o2::constants::physics::MassK0Short); // Kshort

          mother = daughter1 + daughter2; // invariant mass of Kshort pair
          const double deltaMass = deltaM(t1.mK0Short(), t2.mK0Short());

          if (!config.qAOptimisation) {
            if (deltaMass > config.cMaxDeltaM) {
              continue;
            }
          }

          isMix = true;
          if (!config.qAOptimisation)
            fillInvMass(mother, multiplicity, daughter1, daughter2, isMix);

          if (config.qAOptimisation) {
            double deltaRvalue = std::sqrt(TVector2::Phi_mpi_pi(daughter1.phi() - daughter2.phi()) * TVector2::Phi_mpi_pi(daughter1.phi() - daughter2.phi()) + (daughter1.eta() - daughter2.eta()) * (daughter1.eta() - daughter2.eta()));
            const double deltaMass = deltaM(t1.mK0Short(), t2.mK0Short());
            const double ptCorr = std::abs(daughter1.Pt() + daughter2.Pt()) / mother.Pt();
            if (std::abs(mother.Rapidity()) < config.rapidityMotherData) {
              hglue.fill(HIST("h3glueInvMassME"), multiplicity, mother.Pt(), mother.M(), deltaMass, deltaRvalue, ptCorr);
            }
          }
        }
      }
    };
    // Call mixing based on selected estimator
    if (config.cSelectMultEstimator == kFT0M) {
      runMixing(pair1, [](const auto& c) { return c.centFT0M(); });
    } else if (config.cSelectMultEstimator == kFT0A) {
      runMixing(pair2, [](const auto& c) { return c.centFT0A(); });
    } else if (config.cSelectMultEstimator == kFT0C) {
      runMixing(pair3, [](const auto& c) { return c.centFT0C(); });
    } else if (config.cSelectMultEstimator == kFV0A) {
      runMixing(pair4, [](const auto& c) { return c.centFV0A(); });
    }
  }
  PROCESS_SWITCH(HigherMassResonances, processME, "mixed event process", true);

  int counter = 0;
  float multiplicityGen = 0.0;
  std::vector<bool> passKs;
  ROOT::Math::PxPyPzMVector lResonanceGen1;
  ROOT::Math::PxPyPzEVector lResonanceGen;

  void processGen(aod::McCollision const& mcCollision, aod::McParticles const& mcParticles, const soa::SmallGroups<EventCandidatesMC>& collisions)
  {
    // if (config.isMC == false) {
    //   return;
    // }
    hMChists.fill(HIST("events_check"), 0.5);

    std::vector<int64_t> selectedEvents(collisions.size());
    int nevts = 0;
    multiplicityGen = -999.0;
    for (const auto& collision : collisions) {

      // multiplicityGen = collision.centFT0M();
      if (config.cSelectMultEstimator == kFT0M) {
        multiplicityGen = collision.centFT0M();
      } else if (config.cSelectMultEstimator == kFT0A) {
        multiplicityGen = collision.centFT0A();
      } else if (config.cSelectMultEstimator == kFT0C) {
        multiplicityGen = collision.centFT0C();
      } else if (config.cSelectMultEstimator == kFV0A) {
        multiplicityGen = collision.centFV0A();
      } else {
        multiplicityGen = collision.centFT0M(); // default
      }

      if (!selectionEvent(collision, true)) {
        continue;
      }

      selectedEvents[nevts++] = collision.mcCollision_as<aod::McCollisions>().globalIndex();
    }
    selectedEvents.resize(nevts);
    hMChists.fill(HIST("events_check"), 1.5);
    const auto evtReconstructedAndSelected = std::find(selectedEvents.begin(), selectedEvents.end(), mcCollision.globalIndex()) != selectedEvents.end();

    if (!config.isallGenCollisions && !evtReconstructedAndSelected) { // Check that the event is reconstructed and that the reconstructed events pass the selection
      return;
    }
    hMChists.fill(HIST("events_check"), 2.5);
    for (const auto& mcParticle : mcParticles) {

      if (std::abs(mcParticle.pdgCode()) != config.pdgCodes[config.selectMCparticles]) // f2(1525), f0(1710)
      {
        continue;
      }
      hMChists.fill(HIST("events_check"), 3.5);

      if (config.isapplyRapidityMC && std::abs(mcParticle.y()) >= config.rapidityMotherData) {
        continue;
      }
      hMChists.fill(HIST("events_check"), 4.5);

      auto kDaughters = mcParticle.daughters_as<aod::McParticles>();
      if (kDaughters.size() != config.noOfDaughters) {
        continue;
      }
      hMChists.fill(HIST("events_check"), 5.5);

      for (const auto& kCurrentDaughter : kDaughters) {
        // int daupdg = std::abs(kCurrentDaughter.pdgCode());

        if (!kCurrentDaughter.isPhysicalPrimary()) {
          continue;
        }
        hMChists.fill(HIST("events_check"), 6.5);
        if (std::abs(kCurrentDaughter.pdgCode()) == PDG_t::kK0Short) {
          passKs.push_back(true);
          hMChists.fill(HIST("events_check"), 7.5);
          if (passKs.size() == 1) {
            daughter1 = ROOT::Math::PxPyPzMVector(kCurrentDaughter.px(), kCurrentDaughter.py(), kCurrentDaughter.pz(), o2::constants::physics::MassK0Short);
          } else if (static_cast<int>(passKs.size()) == config.noOfDaughters) {
            daughter2 = ROOT::Math::PxPyPzMVector(kCurrentDaughter.px(), kCurrentDaughter.py(), kCurrentDaughter.pz(), o2::constants::physics::MassK0Short);
          }
        }
      }
      if (static_cast<int>(passKs.size()) == config.noOfDaughters) {
        lResonanceGen = ROOT::Math::PxPyPzEVector(mcParticle.pt(), mcParticle.eta(), mcParticle.phi(), mcParticle.e());
        lResonanceGen1 = daughter1 + daughter2;

        ROOT::Math::Boost boost{lResonanceGen.BoostToCM()};
        ROOT::Math::Boost boost1{lResonanceGen1.BoostToCM()};

        fourVecDauCM = boost(daughter1);
        fourVecDauCM1 = boost1(daughter1);

        auto helicityGen = lResonanceGen.Vect().Dot(fourVecDauCM.Vect()) / (std::sqrt(fourVecDauCM.Vect().Mag2()) * std::sqrt(lResonanceGen.Vect().Mag2()));
        auto helicityGen1 = lResonanceGen1.Vect().Dot(fourVecDauCM1.Vect()) / (std::sqrt(fourVecDauCM1.Vect().Mag2()) * std::sqrt(lResonanceGen1.Vect().Mag2()));

        hMChists.fill(HIST("Genf1710"), multiplicityGen, lResonanceGen.pt(), helicityGen);
        hMChists.fill(HIST("Genf1710_mass"), lResonanceGen.M());
        hMChists.fill(HIST("GenRapidity"), mcParticle.y());
        hMChists.fill(HIST("GenEta"), mcParticle.eta());
        hMChists.fill(HIST("GenPhi"), mcParticle.phi());

        if (config.isapplyPairRapidityMC && std::abs(lResonanceGen1.Rapidity()) >= config.rapidityMotherData) {
          continue;
        }

        hMChists.fill(HIST("Genf17102"), multiplicityGen, lResonanceGen1.pt(), helicityGen1);
        hMChists.fill(HIST("Genf1710_mass2"), lResonanceGen1.M());
        hMChists.fill(HIST("GenRapidity2"), lResonanceGen1.Rapidity());
        hMChists.fill(HIST("GenEta2"), lResonanceGen1.Eta());
        hMChists.fill(HIST("GenPhi2"), lResonanceGen1.Phi());
      }
      passKs.clear(); // clear the vector for the next iteration
    }
  }
  PROCESS_SWITCH(HigherMassResonances, processGen, "Process Generated", false);

  int eventCounter = 0;
  std::vector<int> gindex1, gindex2;
  void processRec(EventCandidatesMC::iterator const& collision, TrackCandidatesMC const&, V0TrackCandidatesMC const& V0s, aod::McParticles const&, aod::McCollisions const& /*mcCollisions*/)
  {
    // if (config.isMC == false) {
    //   return;
    // }

    auto multiplicity = -999.0;
    if (config.cSelectMultEstimator == kFT0M) {
      multiplicity = collision.centFT0M();
    } else if (config.cSelectMultEstimator == kFT0A) {
      multiplicity = collision.centFT0A();
    } else if (config.cSelectMultEstimator == kFT0C) {
      multiplicity = collision.centFT0C();
    } else if (config.cSelectMultEstimator == kFV0A) {
      multiplicity = collision.centFV0A();
    } else {
      multiplicity = collision.centFT0M(); // default
    }

    if (!selectionEvent(collision, false)) {
      return;
    }
    hMChists.fill(HIST("Rec_Multiplicity"), multiplicity);

    if (!collision.has_mcCollision()) {
      return;
    }

    hMChists.fill(HIST("MC_mult_after_event_sel"), multiplicity);
    eventCounter++;

    for (const auto& v01 : V0s) {

      for (const auto& v02 : V0s) {

        if (v02.index() <= v01.index()) {
          continue;
        }

        if (!v01.has_mcParticle() || !v02.has_mcParticle()) {
          continue;
        }
        hMChists.fill(HIST("events_checkrec"), 0.5);

        auto postrack1 = v01.template posTrack_as<TrackCandidatesMC>();
        auto negtrack1 = v01.template negTrack_as<TrackCandidatesMC>();

        auto postrack2 = v02.template posTrack_as<TrackCandidatesMC>();
        auto negtrack2 = v02.template negTrack_as<TrackCandidatesMC>();

        if (!postrack1.has_mcParticle() || !postrack2.has_mcParticle())
          continue; // Checking that the daughter tracks come from particles and are not fake

        if (!negtrack1.has_mcParticle() || !negtrack2.has_mcParticle())
          continue;

        double nTPCSigmaPos1[1]{postrack1.tpcNSigmaPi()};
        double nTPCSigmaNeg1[1]{negtrack1.tpcNSigmaPi()};
        double nTPCSigmaPos2[1]{postrack2.tpcNSigmaPi()};
        double nTPCSigmaNeg2[1]{negtrack2.tpcNSigmaPi()};

        if (!isSelectedV0Daughter(postrack1, 1, nTPCSigmaPos1[0], v01) || !isSelectedV0Daughter(postrack2, 1, nTPCSigmaPos2[0], v02)) {
          continue;
        }

        if (!isSelectedV0Daughter(negtrack1, -1, nTPCSigmaNeg1[0], v01) || !isSelectedV0Daughter(negtrack2, -1, nTPCSigmaNeg2[0], v02)) {
          continue;
        }
        hMChists.fill(HIST("events_checkrec"), 1.5);

        if (!selectionV0(collision, v01, multiplicity) || !selectionV0(collision, v02, multiplicity)) {
          continue;
        }
        hMChists.fill(HIST("events_checkrec"), 2.5);

        auto mctrackv01 = v01.mcParticle();
        auto mctrackv02 = v02.mcParticle();

        int trackv0PDG1 = std::abs(mctrackv01.pdgCode());
        int trackv0PDG2 = std::abs(mctrackv02.pdgCode());

        if (std::abs(trackv0PDG1) != PDG_t::kK0Short || std::abs(trackv0PDG2) != PDG_t::kK0Short) {
          continue;
        }
        hMChists.fill(HIST("events_checkrec"), 3.5);

        for (const auto& mothertrack1 : mctrackv01.mothers_as<aod::McParticles>()) {

          // int motpdgs = std::abs(mothertrack1.pdgCode());
          gindex1.push_back(mothertrack1.globalIndex());
          if (gindex1.size() > 1) {
            if (std::find(gindex1.begin(), gindex1.end(), mothertrack1.globalIndex()) != gindex1.end()) {
              continue;
            }
          }

          for (const auto& mothertrack2 : mctrackv02.mothers_as<aod::McParticles>()) {
            hMChists.fill(HIST("events_checkrec"), 4.5);

            if (mothertrack1.pdgCode() != config.pdgCodes[config.selectMCparticles]) {
              continue;
            }
            hMChists.fill(HIST("events_checkrec"), 5.5);

            if (mothertrack1.pdgCode() != mothertrack2.pdgCode()) {
              continue;
            }
            hMChists.fill(HIST("events_checkrec"), 6.5);

            gindex2.push_back(mothertrack2.globalIndex());
            if (gindex2.size() > 1) {
              if (std::find(gindex2.begin(), gindex2.end(), mothertrack2.globalIndex()) != gindex2.end()) {
                continue;
              }
            }
            hMChists.fill(HIST("events_checkrec"), 7.5);

            if (mothertrack1.globalIndex() != mothertrack2.globalIndex()) {
              continue;
            }
            hMChists.fill(HIST("events_checkrec"), 8.5);

            if (!mothertrack1.producedByGenerator()) {
              continue;
            }
            hMChists.fill(HIST("events_checkrec"), 9.5);

            if (config.isapplyRapidityMC && std::abs(mothertrack1.y()) >= config.rapidityMotherData) {
              continue;
            }
            hMChists.fill(HIST("events_checkrec"), 10.5);

            // if (config.isavoidsplitrackMC && oldindex == mothertrack1.globalIndex()) {
            //   hMChists.fill(HIST("h1Recsplit"), mothertrack1.pt());
            //   continue;
            // }
            // hMChists.fill(HIST("events_checkrec"), 11.5);
            // oldindex = mothertrack1.globalIndex(); // split tracks is already handled using gindex1 and gindex2

            daughter1 = ROOT::Math::PxPyPzMVector(v01.px(), v01.py(), v01.pz(), o2::constants::physics::MassK0Short);
            daughter2 = ROOT::Math::PxPyPzMVector(v02.px(), v02.py(), v02.pz(), o2::constants::physics::MassK0Short);
            mother = daughter1 + daughter2;
            mother1 = ROOT::Math::PxPyPzEVector(mothertrack1.px(), mothertrack1.py(), mothertrack1.pz(), mothertrack1.e());

            ROOT::Math::Boost boost{mother.BoostToCM()};
            ROOT::Math::Boost boost1{mother1.BoostToCM()};

            fourVecDauCM = boost(daughter1);
            fourVecDauCM1 = boost1(daughter1);

            auto helicityRec = mother.Vect().Dot(fourVecDauCM.Vect()) / (std::sqrt(fourVecDauCM.Vect().Mag2()) * std::sqrt(mother.Vect().Mag2()));

            auto helicityRec2 = mother1.Vect().Dot(fourVecDauCM1.Vect()) / (std::sqrt(fourVecDauCM1.Vect().Mag2()) * std::sqrt(mother1.Vect().Mag2()));

            // const double deltaMassRec = deltaM(mctrackv01.mK0Short(), mctrackv02.mK0Short());
            // if (deltaMassRec > config.cMaxDeltaM) {
            //   continue;
            // }

            hMChists.fill(HIST("Recf1710_pt1"), multiplicity, mothertrack1.pt(), mother1.M(), helicityRec2);
            hMChists.fill(HIST("RecRapidity"), mothertrack1.y());
            hMChists.fill(HIST("RecPhi"), mothertrack1.phi());
            hMChists.fill(HIST("RecEta"), mothertrack1.eta());

            if (config.isapplyPairRapidityMC && std::abs(mother.Rapidity()) >= config.rapidityMotherData) {
              continue;
            }

            hMChists.fill(HIST("Recf1710_pt2"), multiplicity, mother.Pt(), mother.M(), helicityRec);
            hMChists.fill(HIST("RecRapidity2"), mother.Rapidity());
            hMChists.fill(HIST("RecPhi2"), mother.Phi());
            hMChists.fill(HIST("RecEta2"), mother.Eta());
          }
          gindex2.clear();
        }
        gindex1.clear();
      }
    }
  }
  PROCESS_SWITCH(HigherMassResonances, processRec, "Process Reconstructed", false);

  void processSEderived(EventCandidatesDerivedData::iterator const& collision, V0CandidatesDerivedData const& V0s, DauTracks const&)
  {
    multiplicity = 0.0;
    if (config.cSelectMultEstimator == kFT0M) {
      multiplicity = collision.centFT0M();
    } else if (config.cSelectMultEstimator == kFT0A) {
      multiplicity = collision.centFT0A();
    } else if (config.cSelectMultEstimator == kFT0C) {
      multiplicity = collision.centFT0C();
    } else if (config.cSelectMultEstimator == kFV0A) {
      multiplicity = collision.centFV0A();
    } else {
      multiplicity = collision.centFT0M(); // default
    }

    if (!selectionEvent(collision, true)) {
      return;
    }

    if (config.qAevents) {
      rEventSelection.fill(HIST("hVertexZRec"), collision.posZ());
      rEventSelection.fill(HIST("hmultiplicity"), multiplicity);
    }

    std::vector<int> v0indexes;
    bool allConditionsMet = 0;

    for (const auto& [v1, v2] : combinations(CombinationsFullIndexPolicy(V0s, V0s))) {

      if (v1.size() == 0 || v2.size() == 0) {
        continue;
      }

      if (!selectionV0(collision, v1, multiplicity)) {
        continue;
      }
      if (!selectionV0(collision, v2, multiplicity)) {
        continue;
      }

      if (!isSelectedK0sDaughtersDerived(v1) || !isSelectedK0sDaughtersDerived(v2)) {
        continue;
      }

      if (std::find(v0indexes.begin(), v0indexes.end(), v1.globalIndex()) == v0indexes.end()) {
        v0indexes.push_back(v1.globalIndex());
      }

      if (v2.globalIndex() <= v1.globalIndex()) {
        continue;
      }

      allConditionsMet = 1;
      daughter1 = ROOT::Math::PxPyPzMVector(v1.px(), v1.py(), v1.pz(), o2::constants::physics::MassK0Short); // Kshort
      daughter2 = ROOT::Math::PxPyPzMVector(v2.px(), v2.py(), v2.pz(), o2::constants::physics::MassK0Short); // Kshort

      mother = daughter1 + daughter2; // invariant mass of Kshort pair
      isMix = false;

      const double deltaMass = deltaM(v1.mK0Short(), v2.mK0Short());
      if (config.qAv0) {
        rKzeroShort.fill(HIST("hK0ShortMassCorr"), v1.mK0Short(), v2.mK0Short(), deltaMass);
      }

      if (deltaMass > config.cMaxDeltaM) {
        continue;
      }

      // if (config.qAv0) {
      //   rKzeroShort.fill(HIST("hK0ShortMassCorrAfterCut"), v1.mK0Short(), v2.mK0Short());
      // }

      const double ptCorr = (mother.Pt() - daughter1.Pt() != 0.) ? daughter1.Pt() / (mother.Pt() - daughter1.Pt()) : 0.;
      if (config.qAv0) {
        rKzeroShort.fill(HIST("hK0sPtCorrelation"), ptCorr);
      }

      double deltaRvalue = std::sqrt(TVector2::Phi_mpi_pi(v1.phi() - v2.phi()) * TVector2::Phi_mpi_pi(v1.phi() - v2.phi()) + (v1.eta() - v2.eta()) * (v1.eta() - v2.eta()));

      if (!config.qAOptimisation) {
        if (deltaRvalue < config.deltaRK0sCut) {
          continue;
        }
      }

      if (!config.isselectTWOKsOnly && !config.qAOptimisation)
        fillInvMass(mother, multiplicity, daughter1, daughter2, isMix);

      if (!config.isselectTWOKsOnly && config.qAOptimisation) {

        if (std::abs(mother.Rapidity()) < config.rapidityMotherData) {
          hglue.fill(HIST("h3glueInvMassDS"), multiplicity, mother.Pt(), mother.M(), deltaMass, deltaRvalue, ptCorr);
        }

        for (int i = 0; i < config.cRotations; i++) {
          double theta2 = rn->Uniform(o2::constants::math::PI - o2::constants::math::PI / config.rotationalCut, o2::constants::math::PI + o2::constants::math::PI / config.rotationalCut);

          daughterRot = ROOT::Math::PxPyPzMVector(daughter1.Px() * std::cos(theta2) - daughter1.Py() * std::sin(theta2), daughter1.Px() * std::sin(theta2) + daughter1.Py() * std::cos(theta2), daughter1.Pz(), daughter1.M());

          motherRot = daughterRot + daughter2;
          double pTcorrRot = std::abs(daughterRot.Pt() + daughter2.Pt()) / motherRot.Pt();
          if (motherRot.Rapidity() < config.rapidityMotherData)
            hglue.fill(HIST("h3glueInvMassRot"), multiplicity, motherRot.Pt(), motherRot.M(), deltaMass, deltaRvalue, pTcorrRot);
        }
      }
    }
    int sizeofv0indexes = v0indexes.size();
    rKzeroShort.fill(HIST("NksProduced"), sizeofv0indexes);
    if (config.isselectTWOKsOnly && sizeofv0indexes == config.noOfDaughters && allConditionsMet) {
      fillInvMass(mother, multiplicity, daughter1, daughter2, false);
    }
    v0indexes.clear();
  }
  PROCESS_SWITCH(HigherMassResonances, processSEderived, "same event process in strangeness derived data", false);

  void processMEderived(EventCandidatesDerivedData const& collisions, V0CandidatesDerivedData const& v0s, DauTracks const&)
  {

    for (const auto& [c1, c2] : selfCombinations(colBinning, config.cfgNmixedEvents, -1, collisions, collisions)) // two different centrality c1 and c2 and tracks corresponding to them
    {

      multiplicity = 0.0;
      multiplicity = c1.centFT0M();

      if (!selectionEvent(c1, false) || !selectionEvent(c2, false)) {
        continue;
      }

      auto groupV01 = v0s.sliceBy(tracksPerCollisionV0Mixed, c1.index());
      auto groupV02 = v0s.sliceBy(tracksPerCollisionV0Mixed, c2.index());
      for (const auto& [t1, t2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(groupV01, groupV02))) {

        if (t1.size() == 0 || t2.size() == 0) {
          continue;
        }

        if (!selectionV0(c1, t1, multiplicity))
          continue;
        if (!selectionV0(c2, t2, multiplicity))
          continue;

        if (!isSelectedK0sDaughtersDerived(t1) || !isSelectedK0sDaughtersDerived(t2)) {
          continue;
        }

        if (config.isApplyEtaCutK0s && (t1.eta() < config.confDaughEta || t2.eta() < config.confDaughEta)) {
          continue;
        }

        daughter1 = ROOT::Math::PxPyPzMVector(t1.px(), t1.py(), t1.pz(), o2::constants::physics::MassK0Short); // Kshort
        daughter2 = ROOT::Math::PxPyPzMVector(t2.px(), t2.py(), t2.pz(), o2::constants::physics::MassK0Short); // Kshort

        mother = daughter1 + daughter2; // invariant mass of Kshort pair
        const double deltaMass = deltaM(t1.mK0Short(), t2.mK0Short());

        if (!config.qAOptimisation) {
          if (deltaMass > config.cMaxDeltaM) {
            continue;
          }
        }

        isMix = true;
        if (!config.qAOptimisation)
          fillInvMass(mother, multiplicity, daughter1, daughter2, isMix);

        if (config.qAOptimisation) {
          double deltaRvalue = std::sqrt(TVector2::Phi_mpi_pi(daughter1.phi() - daughter2.phi()) * TVector2::Phi_mpi_pi(daughter1.phi() - daughter2.phi()) + (daughter1.eta() - daughter2.eta()) * (daughter1.eta() - daughter2.eta()));
          const double deltaMass = deltaM(t1.mK0Short(), t2.mK0Short());
          const double ptCorr = std::abs(daughter1.Pt() + daughter2.Pt()) / mother.Pt();
          if (std::abs(mother.Rapidity()) < config.rapidityMotherData) {
            hglue.fill(HIST("h3glueInvMassME"), multiplicity, mother.Pt(), mother.M(), deltaMass, deltaRvalue, ptCorr);
          }
        }
      }
    }
  }
  PROCESS_SWITCH(HigherMassResonances, processMEderived, "mixed event process in derived data", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HigherMassResonances>(cfgc)};
}
