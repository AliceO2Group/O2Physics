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

/// \file xi1820Analysis.cxx
/// \brief Invariant Mass Reconstruction of Xi(1820) Resonance
///
/// \author Adel Bensaoula <adel.bensaoula@cern.ch>, Bong-Hwi Lim <bong-hwi.lim@cern.ch>, Minjae Kim <minjae.kim@cern.ch>

#include "PWGLF/DataModel/LFResonanceTables.h"

#include "Common/Core/RecoDecay.h"

#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Configurable.h>
#include <Framework/GroupedCombinations.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/SliceCache.h>
#include <Framework/runDataProcessing.h>

#include <Math/GenVector/RotationZ.h>
#include <Math/Vector4D.h> // IWYU pragma: keep (do not replace with Math/Vector4Dfwd.h)
#include <Math/Vector4Dfwd.h>
#include <TPDGCode.h>

#include <arrow/array/concatenate.h>
#include <arrow/chunked_array.h>
#include <arrow/table.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <deque>
#include <map>
#include <memory>
#include <tuple>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::constants::physics;
using LorentzVectorSetXYZM = ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<float>>;

struct Xi1820Analysis {
  // Module-initializer tables; full tracks and V0s retain their schema.
  using ResoCollisions = aod::ResoCollisions_001;
  // Original track IDs are used as numbers; no parent AO2D dereference is needed.
  using ResoTracks = soa::Join<aod::ResoTracks, aod::ResoTrackTracks>;
  using ResoMicroTracks = aod::ResoMicroTracks_001;
  using ResoMCCols = soa::Join<ResoCollisions, aod::ResoMCCollisions_001>;

  SliceCache cache;
  Preslice<aod::ResoV0s> perResoCollisionV0 = aod::resodaughter::resoCollisionId;
  Preslice<aod::ResoTracks> perResoCollisionTrack = aod::resodaughter::resoCollisionId;
  Preslice<ResoMicroTracks> perResoCollisionMicroTrack = aod::resodaughter::resoCollisionId;
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  struct MixingCollision {
    float x, y, z, centrality;
    bool recINELgt0;
    [[nodiscard]] float posX() const { return x; }
    [[nodiscard]] float posY() const { return y; }
    [[nodiscard]] float posZ() const { return z; }
    [[nodiscard]] float cent() const { return centrality; }
  };
  struct MixingEvent {
    MixingCollision collision{};
    std::shared_ptr<arrow::Table> daughters;
  };
  std::map<int, std::deque<MixingEvent>> mixingPools;
  float mixingBField = 0.f;
  std::map<int, std::deque<MixingEvent>> neutralMixingPools;
  float neutralMixingBField = 0.f;

  template <typename Iterator>
  struct SelectedLambda {
    Iterator candidate;
    bool isLambda = false, isAntiLambda = false;
  };

  Configurable<bool> cfgFillQA{"cfgFillQA", true, "Fill reconstructed daughter QA"};
  Configurable<bool> cfgFillEventQA{"cfgFillEventQA", true, "Fill event and rotation QA"};
  Configurable<bool> cfgFillTruthQA{"cfgFillTruthQA", true, "Fill MC truth-matched daughter QA"};

  // Constants
  static constexpr float SmallMomentumDenominator = 1e-10f; // Small value to avoid division by zero
  static constexpr int PdgChargedXi1820 = 123314;           // o2-linter: disable=pdg/explicit-code (Xi(1820) PDG code not available in PDG_t or o2::constants::physics::Pdg)
  static constexpr int PdgXi1820Zero = 123324;              // o2-linter: disable=pdg/explicit-code (Xi(1820) PDG code not available in PDG_t or o2::constants::physics::Pdg)
  static constexpr int ExpectedDaughters = 2;               // Expected number of daughters for two-body decay

  // Axes
  ConfigurableAxis binsPt{"binsPt", {VARIABLE_WIDTH, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0}, "pT"};
  ConfigurableAxis binsPtQA{"binsPtQA", {VARIABLE_WIDTH, 0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0}, "pT (QA)"};
  ConfigurableAxis binsCent{"binsCent", {VARIABLE_WIDTH, 0., 1., 5., 10., 30., 50., 70., 100., 110.}, "Centrality"};

  // Invariant mass range for Xi(1820) to Λ + K
  Configurable<float> cInvMassStart{"cInvMassStart", 1.6, "Invariant mass start (GeV/c^2)"};
  Configurable<float> cInvMassEnd{"cInvMassEnd", 2.2, "Invariant mass end (GeV/c^2)"};
  Configurable<int> cInvMassBins{"cInvMassBins", 600, "Invariant mass bins"};

  // Basic pre-selections
  Configurable<float> cMinPtcut{"cMinPtcut", 0.15, "Minimum pT for candidates"};
  Configurable<float> cMaxEtaCut{"cMaxEtaCut", 0.8, "Maximum |eta|"};

  // Kaon track selections
  Configurable<float> cKaonPtMin{"cKaonPtMin", 0.15, "Minimum kaon pT"};
  Configurable<float> cKaonEtaMax{"cKaonEtaMax", 0.8, "Maximum kaon |eta|"};
  Configurable<float> cKaonDCAxyMax{"cKaonDCAxyMax", 0.1, "Maximum kaon DCAxy to PV"};
  Configurable<float> cKaonDCAzMax{"cKaonDCAzMax", 0.1, "Maximum kaon DCAz to PV"};
  Configurable<int> cKaonTPCNClusMin{"cKaonTPCNClusMin", 70, "Minimum TPC clusters for kaon"};
  Configurable<int> cKaonITSNClusMin{"cKaonITSNClusMin", 2, "Minimum ITS clusters for kaon"};

  // Kaon PID selections
  Configurable<float> cKaonTPCNSigmaMin{"cKaonTPCNSigmaMin", -3.5f, "Strict lower TPC (nSigma - mean) bound for kaons"};
  Configurable<float> cKaonTPCNSigmaMean{"cKaonTPCNSigmaMean", 0.f, "TPC nSigma mean subtracted for kaons"};
  Configurable<float> cKaonTPCNSigmaMax{"cKaonTPCNSigmaMax", 3.5, "Maximum TPC NSigma for kaon (if not using pT-dependent)"};
  Configurable<float> cKaonTOFNSigmaMin{"cKaonTOFNSigmaMin", -999.f, "Strict lower TOF (nSigma - mean) bound for kaons"};
  Configurable<float> cKaonTOFNSigmaMean{"cKaonTOFNSigmaMean", 0.f, "TOF nSigma mean subtracted for kaons"};
  Configurable<bool> cKaonBypassTOF{"cKaonBypassTOF", false, "Bypass TOF PID even when TOF information is present"};
  Configurable<float> cKaonTOFNSigmaMax{"cKaonTOFNSigmaMax", 999., "Maximum TOF NSigma for kaon (if not using pT-dependent)"};
  Configurable<bool> cKaonUsePtDepPID{"cKaonUsePtDepPID", false, "Use pT-dependent PID cuts"};
  Configurable<std::vector<float>> cKaonPIDPtBins{"cKaonPIDPtBins", {0.0f, 0.5f, 0.8f, 2.0f, 999.0f}, "pT bin edges for PID cuts (N+1 values for N bins)"};
  Configurable<std::vector<float>> cKaonTPCNSigmaCuts{"cKaonTPCNSigmaCuts", {3.0f, 3.0f, 2.0f, 2.0f}, "TPC NSigma cuts per pT bin (N values)"};
  Configurable<std::vector<float>> cKaonTOFNSigmaCuts{"cKaonTOFNSigmaCuts", {3.0f, 3.0f, 3.0f, 3.0f}, "TOF NSigma cuts per pT bin (N values)"};
  Configurable<std::vector<float>> cKaonTPCNSigmaMinCuts{"cKaonTPCNSigmaMinCuts", {-3.f, -3.f, -2.f, -2.f}, "Strict lower centered TPC bounds per pT bin"};
  Configurable<std::vector<float>> cKaonTOFNSigmaMinCuts{"cKaonTOFNSigmaMinCuts", {-3.f, -3.f, -3.f, -3.f}, "Strict lower centered TOF bounds per pT bin"};

  // V0 (Lambda) selections
  Configurable<double> cV0MinCosPA{"cV0MinCosPA", 0.995, "V0 minimum pointing angle cosine"};
  Configurable<double> cV0MaxDaughDCA{"cV0MaxDaughDCA", 0.5, "V0 daughter DCA Maximum"};
  Configurable<double> cV0MassWindow{"cV0MassWindow", 0.01, "Mass window for Lambda selection (GeV/c^2)"};
  Configurable<double> cMaxV0Etacut{"cMaxV0Etacut", 0.8, "V0 maximum eta cut"};
  Configurable<float> cV0RadiusMin{"cV0RadiusMin", 0.5, "V0 decay radius min"};
  Configurable<float> cV0RadiusMax{"cV0RadiusMax", 200.0, "V0 decay radius max"};
  Configurable<float> cV0DauPosDCAtoPVMin{"cV0DauPosDCAtoPVMin", 0.05, "V0 positive daughter DCA to PV min"};
  Configurable<float> cV0DauNegDCAtoPVMin{"cV0DauNegDCAtoPVMin", 0.05, "V0 negative daughter DCA to PV min"};
  Configurable<float> cV0ProperLifetimeMax{"cV0ProperLifetimeMax", 30.0, "Lambda proper lifetime max (cm)"};
  Configurable<bool> cV0sCrossMassRejection{"cV0sCrossMassRejection", true, "Enable K0s mass rejection for Lambda"};
  Configurable<float> cV0sCrossMassRejectionWindow{"cV0sCrossMassRejectionWindow", 0.005, "K0s mass rejection window for Lambda (GeV/c^2)"};
  Configurable<float> cLambdaDaughterPiTPCNSigmaMin{"cLambdaDaughterPiTPCNSigmaMin", -5.f, "Strict lower centered TPC PID bound"};
  Configurable<float> cLambdaDaughterPiTPCNSigmaMean{"cLambdaDaughterPiTPCNSigmaMean", 0.f, "TPC nSigma mean subtracted for the V0 daughter"};
  Configurable<float> cLambdaDaughterPiTPCNSigmaMax{"cLambdaDaughterPiTPCNSigmaMax", 5.0, "Maximum TPC NSigma for Lambda daughter pions"};
  Configurable<float> cLambdaDaughterPrTPCNSigmaMin{"cLambdaDaughterPrTPCNSigmaMin", -5.f, "Strict lower centered TPC PID bound"};
  Configurable<float> cLambdaDaughterPrTPCNSigmaMean{"cLambdaDaughterPrTPCNSigmaMean", 0.f, "TPC nSigma mean subtracted for the V0 daughter"};
  Configurable<float> cLambdaDaughterPrTPCNSigmaMax{"cLambdaDaughterPrTPCNSigmaMax", 5.0, "Maximum TPC NSigma for Lambda daughter protons"};
  Configurable<int> cLambdaPosDaughterMinCrossedRows{"cLambdaPosDaughterMinCrossedRows", 50, "Minimum TPC crossed rows for Lambda positive daughter"};
  Configurable<int> cLambdaNegDaughterMinCrossedRows{"cLambdaNegDaughterMinCrossedRows", 50, "Minimum TPC crossed rows for Lambda negative daughter"};

  // K0s selections
  Configurable<double> cK0sMinCosPA{"cK0sMinCosPA", 0.98, "K0s minimum pointing angle cosine"};
  Configurable<double> cK0sMaxDaughDCA{"cK0sMaxDaughDCA", 0.5, "K0s daughter DCA Maximum"};
  Configurable<double> cK0sMassWindow{"cK0sMassWindow", 0.025, "Mass window for K0s selection (GeV/c^2)"};
  Configurable<float> cK0sProperLifetimeMax{"cK0sProperLifetimeMax", 20.0, "K0s proper lifetime max (cm)"};
  Configurable<float> cK0sArmenterosQtMin{"cK0sArmenterosQtMin", 0.0, "K0s Armenteros qt min"};
  Configurable<float> cK0sArmenterosAlphaCoeff{"cK0sArmenterosAlphaCoeff", 0.2, "K0s Armenteros alpha max"};
  Configurable<float> cK0sDauPosDCAtoPVMin{"cK0sDauPosDCAtoPVMin", 0.05, "K0s positive daughter DCA to PV min"};
  Configurable<float> cK0sDauNegDCAtoPVMin{"cK0sDauNegDCAtoPVMin", 0.05, "K0s negative daughter DCA to PV min"};
  Configurable<float> cK0sRadiusMin{"cK0sRadiusMin", 0.5, "K0s decay radius min"};
  Configurable<float> cK0sRadiusMax{"cK0sRadiusMax", 200.0, "K0s decay radius max"};
  Configurable<bool> cK0sCrossMassRejection{"cK0sCrossMassRejection", true, "Enable Lambda mass rejection for K0s"};
  Configurable<float> cK0sCrossMassRejectionWindow{"cK0sCrossMassRejectionWindow", 0.01, "Lambda mass rejection window for K0s (GeV/c^2)"};
  Configurable<float> cK0sDaughterPiTPCNSigmaMin{"cK0sDaughterPiTPCNSigmaMin", -5.f, "Strict lower centered TPC PID bound"};
  Configurable<float> cK0sDaughterPiTPCNSigmaMean{"cK0sDaughterPiTPCNSigmaMean", 0.f, "TPC nSigma mean subtracted for the V0 daughter"};
  Configurable<float> cK0sDaughterPiTPCNSigmaMax{"cK0sDaughterPiTPCNSigmaMax", 5.0, "Maximum TPC NSigma for K0s daughter pions"};
  Configurable<int> cK0sPosDaughterMinCrossedRows{"cK0sPosDaughterMinCrossedRows", 50, "Minimum TPC crossed rows for K0s positive daughter"};
  Configurable<int> cK0sNegDaughterMinCrossedRows{"cK0sNegDaughterMinCrossedRows", 50, "Minimum TPC crossed rows for K0s negative daughter"};

  // Event Mixing
  Configurable<int> nEvtMixing{"nEvtMixing", 10, "Number of events to mix"};
  ConfigurableAxis cfgVtxBins{"cfgVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis cfgMultBins{"cfgMultBins", {VARIABLE_WIDTH, 0.0f, 1.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.0f, 110.0f}, "Mixing bins - centrality"};

  // Additional QA and configurations
  struct : ConfigurableGroup {
    Configurable<bool> cRecoINELgt0{"cRecoINELgt0", true, "Apply Reco INEL>0 event selection"};
    Configurable<bool> cMCINELgt0{"cMCINELgt0", true, "Require generator INEL>0 in reconstructed MC and generated-parent processes"};
    Configurable<bool> cMCVtxIn10{"cMCVtxIn10", true, "Require generator |vertex z| < 10 cm using isVtxIn10 in reconstructed MC and generated-parent processes"};
    Configurable<bool> cUseTruthRapidity{"cUseTruthRapidity", false, "Use truth rapidity for MC generated target"};

    Configurable<bool> cUsePtDepDCAForKaons{"cUsePtDepDCAForKaons", true, "Use pT dependent DCA cuts for kaon tracks"};
    Configurable<float> cDCAToPVByPtFirstP0{"cDCAToPVByPtFirstP0", 0.004, "pT dependent DCA cut first parameter (cm)"};
    Configurable<float> cDCAToPVByPtFirstExp{"cDCAToPVByPtFirstExp", 0.013, "Coefficient in kaon DCA cut = P0 + coefficient / pT^power (legacy key)"};
    Configurable<float> cDCAToPVByPtFirstPower{"cDCAToPVByPtFirstPower", 1.f, "Power in kaon DCA cut = P0 + coefficient / pT^power"};
    Configurable<float> cMaxDcaToPVV0Lambda{"cMaxDcaToPVV0Lambda", 1.0, "Maximum DCA to PV for Lambda candidates (cm)"};
    Configurable<float> cMaxDcaToPVV0K0s{"cMaxDcaToPVV0K0s", 1.0, "Maximum DCA to PV for K0s candidates (cm)"};

    Configurable<float> cfgRapidityCut{"cfgRapidityCut", 0.5, "Rapidity cut"};
    ConfigurableAxis multNTracksAxis{"multNTracksAxis", {500, 0, 500}, "N_{tracks}"};

    Configurable<bool> cfgFillRotBkg{"cfgFillRotBkg", true, "Fill rotated background"};
    Configurable<float> cfgMinRot{"cfgMinRot", o2::constants::math::PIHalf + o2::constants::math::PIThird, "Minimum of rotation"};
    Configurable<float> cfgMaxRot{"cfgMaxRot", o2::constants::math::TwoPI - o2::constants::math::PIHalf - o2::constants::math::PIThird, "Maximum of rotation"};
    Configurable<bool> cfgRotKaon{"cfgRotKaon", true, "Rotate Kaon"};
    Configurable<int> cfgNrotBkg{"cfgNrotBkg", 10, "Number of rotated copies (background) per each original candidate"};

    Configurable<bool> cfgDelCheck{"cfgDelCheck", true, "Check Delta R distribution of candidates"};
    Configurable<bool> cfgPVContributor{"cfgPVContributor", true, "Require PV contributor tracks"};
    Configurable<bool> cfgPrimaryTrack{"cfgPrimaryTrack", true, "Require primary tracks flags (just for check)"};

  } additionalConfig;

  using BinningTypeVertexContributor = ColumnBinningPolicy<aod::collision::PosZ, aod::resocollision::Cent>;

  void init(InitContext&)
  {
    const int sameEventModes = static_cast<int>(doprocessDataWithTracks) + static_cast<int>(doprocessDataWithMicroTracks) +
                               static_cast<int>(doprocessMCWithTracks) + static_cast<int>(doprocessMCWithMicroTracks);
    const int mixedEventModes = static_cast<int>(doprocessMixedEventWithTracks) + static_cast<int>(doprocessMixedEventWithMicroTracks) + static_cast<int>(doprocessMEDF);
    if (sameEventModes > 1 || mixedEventModes > 1 || (doprocessK0sLambda && doprocessMCK0sLambda) ||
        (doprocessK0sLambdaMixedEvent && doprocessK0sLambdaMEDF)) {
      LOG(fatal) << "Enable at most one same-event mode and one mixing mode per charged/neutral channel";
    }
    AxisSpec centAxis = {binsCent, "V0M (%)"};
    AxisSpec ptAxis = {binsPt, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec ptAxisQA = {binsPtQA, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec invMassAxis = {cInvMassBins, cInvMassStart, cInvMassEnd, "Invariant Mass (GeV/#it{c}^{2})"};
    AxisSpec lambdaMassAxis = {200, 1.08, 1.16, "#Lambda mass (GeV/#it{c}^{2})"};
    AxisSpec dcaAxis = {200, 0., 2.0, "DCA (cm)"};
    AxisSpec dcaxyAxis = {400, -0.2, 0.2, "DCA_{xy} (cm)"};
    AxisSpec dcazAxis = {400, -0.2, 0.2, "DCA_{z} (cm)"};
    AxisSpec cosPAAxis = {1000, 0.95, 1.0, "cos(PA)"};
    AxisSpec radiusAxis = {200, 0, 200, "Radius (cm)"};
    AxisSpec lifetimeAxis = {200, 0, 50, "Proper lifetime (cm)"};
    AxisSpec nsigmaAxis = {100, -5.0, 5.0, "N#sigma"};
    AxisSpec armenterosAlphaAxis = {200, -1.0, 1.0, "Armenteros alpha"};
    AxisSpec armenterosQtAxis = {500, 0.0, 0.5, "Armenteros qt (GeV/c)"};
    AxisSpec axisRot = {360, -o2::constants::math::PI, o2::constants::math::PI};
    AxisSpec axisDel = {400, -0.5, 7.0, "#DeltaR"};

    // Event QA histograms
    if (cfgFillEventQA) {
      histos.add("Event/posZ", "Event vertex Z position", kTH1F, {{200, -20., 20., "V_{z} (cm)"}});
      histos.add("Event/centrality", "Event centrality distribution", kTH1D, {centAxis});
      histos.add("Event/posZvsCent", "Vertex Z vs Centrality", kTH2F, {{200, -20., 20., "V_{z} (cm)"}, centAxis});
      histos.add("Event/nV0s", "Number of V0s per event", kTH1F, {{200, 0., 200., "N_{V0s}"}});
      histos.add("Event/nKaons", "Number of kaons per event", kTH1F, {{200, 0., 200., "N_{kaon}"}});
      histos.add("Event/nLambdasAfterCuts", "Number of Lambdas per event after cuts", kTH1F, {{100, 0., 100., "N_{Lambda}"}});
      histos.add("Event/nKaonsAfterCuts", "Number of charged kaons per event after cuts", kTH1F, {{100, 0., 100., "N_{Kaon}"}});
      if (doprocessK0sLambda || doprocessMCK0sLambda) {
        histos.add("Event/nK0s", "Number of K0s candidates per event before cuts", kTH1F, {{200, 0., 200., "N_{K^{0}_{S}}"}});
        histos.add("Event/nK0sAfterCuts", "Number of K0s per event after cuts", kTH1F, {{100, 0., 100., "N_{K^{0}_{S}}"}});
      }
      histos.add("Event/hRotBkg", "Rotated angle of rotated background", HistType::kTH1F, {axisRot});
    }

    if (cfgFillEventQA && doprocessK0sLambdaMEDF) {
      const int mixingDepth = std::max(0, nEvtMixing.value);
      const AxisSpec mixingCountAxis{mixingDepth + 1, -0.5, mixingDepth + 0.5, "N_{events}"};
      const AxisSpec mixingVertexAxis{200, -20., 20., "Current-event vertex z (cm)"};
      const AxisSpec mixingCentAxis{110, 0., 110., "Current-event centrality (%)"};
      const AxisSpec mixingV0Axis{201, -0.5, 200.5, "Current-event N_{V0s} before candidate cuts"};
      const AxisSpec previousVertexAxis{200, -20., 20., "Previous-event vertex z (cm)"};
      const AxisSpec previousCentAxis{110, 0., 110., "Previous-event centrality (%)"};
      const AxisSpec previousV0Axis{201, -0.5, 200.5, "Previous-event N_{V0s} before candidate cuts"};
      histos.add("MixingQA/K0sLambdaMEDF/hPosZvsCent", "Current mixing events;Vertex z (cm);Centrality (%)", kTH2F, {mixingVertexAxis, mixingCentAxis});
      histos.add("MixingQA/K0sLambdaMEDF/hPoolSize", "Buffered events before mixing, including INEL-rejected slots (one entry per accepted current event)", kTH1F, {mixingCountAxis});
      histos.add("MixingQA/K0sLambdaMEDF/hNPartners", "Accepted event partners before V0 cuts (one entry per accepted current event)", kTH1F, {mixingCountAxis});
      histos.add("MixingQA/K0sLambdaMEDF/hPosZPair", "Mixed event vertices", kTH2F, {previousVertexAxis, mixingVertexAxis});
      histos.add("MixingQA/K0sLambdaMEDF/hCentPair", "Mixed event centralities", kTH2F, {previousCentAxis, mixingCentAxis});
      histos.add("MixingQA/K0sLambdaMEDF/hV0CountsPair", "Mixed event V0 counts before candidate cuts", kTH2F, {previousV0Axis, mixingV0Axis});
    }

    if (cfgFillQA && (doprocessDataWithTracks || doprocessDataWithMicroTracks || doprocessMCWithTracks || doprocessMCWithMicroTracks || doprocessK0sLambda || doprocessMCK0sLambda)) {
      // Lambda QA histograms
      histos.add("QAbefore/lambdaDCAtoPV", "Lambda DCA to PV before cuts", kTH2F, {ptAxisQA, dcaAxis});
      histos.add("QAbefore/lambdaMass", "Lambda mass before cuts", kTH1F, {lambdaMassAxis});
      histos.add("QAbefore/lambdaMassAnti", "Anti-Lambda mass before cuts", kTH1F, {lambdaMassAxis});
      histos.add("QAbefore/lambdaPt", "Lambda pT before cuts", kTH1F, {ptAxisQA});
      histos.add("QAbefore/lambdaEta", "Lambda eta before cuts", kTH1F, {{100, -2.0, 2.0, "#eta"}});
      histos.add("QAbefore/lambdaCosPA", "Lambda CosPA before cuts", kTH2F, {ptAxisQA, cosPAAxis});
      histos.add("QAbefore/lambdaRadius", "Lambda radius before cuts", kTH2F, {ptAxisQA, radiusAxis});
      histos.add("QAbefore/lambdaDauDCA", "Lambda daughter DCA before cuts", kTH2F, {ptAxisQA, dcaAxis});
      histos.add("QAbefore/lambdaProperLifetime", "Lambda proper lifetime before cuts", kTH2F, {ptAxisQA, lifetimeAxis});
      histos.add("QAbefore/lambdaDauPosDCA", "Lambda positive daughter DCA before cuts", kTH2F, {ptAxisQA, dcaAxis});
      histos.add("QAbefore/lambdaDauNegDCA", "Lambda negative daughter DCA before cuts", kTH2F, {ptAxisQA, dcaAxis});
      histos.add("QAbefore/lambdaDaughterTPCNSigmaPosPr", "Lambda daughter proton TPC NSigma before cuts", kTH2F, {ptAxisQA, nsigmaAxis});
      histos.add("QAbefore/lambdaDaughterTPCNSigmaNegPi", "Lambda daughter pi- TPC NSigma before cuts", kTH2F, {ptAxisQA, nsigmaAxis});
      histos.add("QAbefore/lambdaAntiDaughterTPCNSigmaPosPi", "Anti-Lambda daughter pion TPC NSigma before cuts", kTH2F, {ptAxisQA, nsigmaAxis});
      histos.add("QAbefore/lambdaAntiDaughterTPCNSigmaNegPr", "Anti-Lambda daughter proton TPC NSigma before cuts", kTH2F, {ptAxisQA, nsigmaAxis});
      histos.add("QAbefore/lambdaNCrossedRowsPos", "Lambda positive daughter crossed rows before cuts", kTH2F, {ptAxisQA, {160, 0, 160, "TPC crossed rows"}});
      histos.add("QAbefore/lambdaNCrossedRowsNeg", "Lambda negative daughter crossed rows before cuts", kTH2F, {ptAxisQA, {160, 0, 160, "TPC crossed rows"}});
      histos.add("QAbefore/lambdaArmenterosPodolanski", "Lambda candidate Armenteros-Podolanski before cuts", kTH3F, {armenterosAlphaAxis, armenterosQtAxis, ptAxisQA});

      histos.add("QAafter/lambdaDCAtoPV", "Lambda DCA to PV after cuts", kTH2F, {ptAxisQA, dcaAxis});
      histos.add("QAafter/lambdaMass", "Lambda mass after cuts", kTH1F, {lambdaMassAxis});
      histos.add("QAafter/lambdaMassAnti", "Anti-Lambda mass after cuts", kTH1F, {lambdaMassAxis});
      histos.add("QAafter/lambdaPt", "Lambda pT after cuts", kTH1F, {ptAxisQA});
      histos.add("QAafter/lambdaEta", "Lambda eta after cuts", kTH1F, {{100, -2.0, 2.0, "#eta"}});
      histos.add("QAafter/lambdaCosPA", "Lambda CosPA after cuts", kTH2F, {ptAxisQA, cosPAAxis});
      histos.add("QAafter/lambdaRadius", "Lambda radius after cuts", kTH2F, {ptAxisQA, radiusAxis});
      histos.add("QAafter/lambdaDauDCA", "Lambda daughter DCA after cuts", kTH2F, {ptAxisQA, dcaAxis});
      histos.add("QAafter/lambdaProperLifetime", "Lambda proper lifetime after cuts", kTH2F, {ptAxisQA, lifetimeAxis});
      histos.add("QAafter/lambdaDauPosDCA", "Lambda positive daughter DCA after cuts", kTH2F, {ptAxisQA, dcaAxis});
      histos.add("QAafter/lambdaDauNegDCA", "Lambda negative daughter DCA after cuts", kTH2F, {ptAxisQA, dcaAxis});
      histos.add("QAafter/lambdaDaughterTPCNSigmaPosPr", "Lambda daughter proton TPC NSigma after cuts", kTH2F, {ptAxisQA, nsigmaAxis});
      histos.add("QAafter/lambdaDaughterTPCNSigmaNegPi", "Lambda daughter pi- TPC NSigma after cuts", kTH2F, {ptAxisQA, nsigmaAxis});
      histos.add("QAafter/lambdaAntiDaughterTPCNSigmaPosPi", "Anti-Lambda daughter pion TPC NSigma after cuts", kTH2F, {ptAxisQA, nsigmaAxis});
      histos.add("QAafter/lambdaAntiDaughterTPCNSigmaNegPr", "Anti-Lambda daughter proton TPC NSigma after cuts", kTH2F, {ptAxisQA, nsigmaAxis});
      histos.add("QAafter/lambdaNCrossedRowsPos", "Lambda positive daughter crossed rows after cuts", kTH2F, {ptAxisQA, {160, 0, 160, "TPC crossed rows"}});
      histos.add("QAafter/lambdaNCrossedRowsNeg", "Lambda negative daughter crossed rows after cuts", kTH2F, {ptAxisQA, {160, 0, 160, "TPC crossed rows"}});
      histos.add("QAafter/lambdaArmenterosPodolanski", "Lambda candidate Armenteros-Podolanski after cuts", kTH3F, {armenterosAlphaAxis, armenterosQtAxis, ptAxisQA});
    }

    if (cfgFillQA && (doprocessDataWithTracks || doprocessDataWithMicroTracks || doprocessMCWithTracks || doprocessMCWithMicroTracks)) {
      // Kaon QA histograms
      histos.add("QAbefore/kaonPt", "Kaon pT before cuts", kTH1F, {ptAxisQA});
      histos.add("QAbefore/kaonEta", "Kaon eta before cuts", kTH1F, {{100, -2.0, 2.0, "#eta"}});
      histos.add("QAbefore/kaonDCAxy", "Kaon DCAxy before cuts", kTH2F, {ptAxisQA, dcaxyAxis});
      histos.add("QAbefore/kaonDCAz", "Kaon DCAz before cuts", kTH2F, {ptAxisQA, dcazAxis});
      histos.add("QAbefore/kaonTPCNcls", "Kaon TPC clusters before cuts", kTH1F, {{160, 0, 160, "N_{TPC clusters}"}});
      histos.add("QAbefore/kaonITSNcls", "Kaon ITS clusters before cuts", kTH1F, {{10, 0, 10, "N_{ITS clusters}"}});
      histos.add("QAbefore/kaonTPCNSigma", "Kaon TPC NSigma before cuts", kTH2F, {ptAxisQA, nsigmaAxis});
      histos.add("QAbefore/kaonTOFNSigma", "Kaon TOF NSigma before cuts", kTH2F, {ptAxisQA, nsigmaAxis});

      histos.add("QAafter/kaonPt", "Kaon pT after cuts", kTH1F, {ptAxisQA});
      histos.add("QAafter/kaonEta", "Kaon eta after cuts", kTH1F, {{100, -2.0, 2.0, "#eta"}});
      histos.add("QAafter/kaonDCAxy", "Kaon DCAxy after cuts", kTH2F, {ptAxisQA, dcaxyAxis});
      histos.add("QAafter/kaonDCAz", "Kaon DCAz after cuts", kTH2F, {ptAxisQA, dcazAxis});
      histos.add("QAafter/kaonTPCNcls", "Kaon TPC clusters after cuts", kTH1F, {{160, 0, 160, "N_{TPC clusters}"}});
      histos.add("QAafter/kaonITSNcls", "Kaon ITS clusters after cuts", kTH1F, {{10, 0, 10, "N_{ITS clusters}"}});
      histos.add("QAafter/kaonTPCNSigma", "Kaon TPC NSigma after cuts", kTH2F, {ptAxisQA, nsigmaAxis});
      histos.add("QAafter/kaonTOFNSigma", "Kaon TOF NSigma after cuts", kTH2F, {ptAxisQA, nsigmaAxis});
    }

    // Resonance histograms - 4 combinations
    // K+ Lambda
    if (doprocessDataWithTracks || doprocessDataWithMicroTracks || doprocessMixedEventWithTracks || doprocessMixedEventWithMicroTracks || doprocessMEDF || doprocessMCWithTracks || doprocessMCWithMicroTracks) {

      histos.add("xi1820/kplus_lambda/hInvMassKplusLambda", "Invariant mass of K^{+} + #Lambda", kTH1F, {invMassAxis});
      histos.add("xi1820/kplus_lambda/hInvMassKplusLambda_Mix", "Mixed event Invariant mass of K^{+} + #Lambda", kTH1F, {invMassAxis});

      // K+ Anti-Lambda
      histos.add("xi1820/kplus_antilambda/hInvMassKplusAntiLambda", "Invariant mass of K^{+} + #bar{#Lambda}", kTH1F, {invMassAxis});
      histos.add("xi1820/kplus_antilambda/hInvMassKplusAntiLambda_Mix", "Mixed event Invariant mass of K^{+} + #bar{#Lambda}", kTH1F, {invMassAxis});

      // K- Lambda
      histos.add("xi1820/kminus_lambda/hInvMassKminusLambda", "Invariant mass of K^{-} + #Lambda", kTH1F, {invMassAxis});
      histos.add("xi1820/kminus_lambda/hInvMassKminusLambda_Mix", "Mixed event Invariant mass of K^{-} + #Lambda", kTH1F, {invMassAxis});

      // K- Anti-Lambda
      histos.add("xi1820/kminus_antilambda/hInvMassKminusAntiLambda", "Invariant mass of K^{-} + #bar{#Lambda}", kTH1F, {invMassAxis});
      histos.add("xi1820/kminus_antilambda/hInvMassKminusAntiLambda_Mix", "Mixed event Invariant mass of K^{-} + #bar{#Lambda}", kTH1F, {invMassAxis});

      if (!(additionalConfig.cfgDelCheck)) {
        histos.add("xi1820/kplus_lambda/hMassPtCentKplusLambda", "K^{+} + #Lambda mass vs pT vs cent", kTH3D, {invMassAxis, ptAxis, centAxis});
        histos.add("xi1820/kplus_lambda/hMassPtCentKplusLambda_Mix", "Mixed event K^{+} + #Lambda mass vs pT vs cent", kTH3D, {invMassAxis, ptAxis, centAxis});

        histos.add("xi1820/kplus_antilambda/hMassPtCentKplusAntiLambda", "K^{+} + #bar{#Lambda} mass vs pT vs cent", kTH3D, {invMassAxis, ptAxis, centAxis});
        histos.add("xi1820/kplus_antilambda/hMassPtCentKplusAntiLambda_Mix", "Mixed event K^{+} + #bar{#Lambda} mass vs pT vs cent", kTH3D, {invMassAxis, ptAxis, centAxis});
        histos.add("xi1820/kplus_antilambda/hMassPtCentKplusAntiLambda_Rot", "Rotated background K^{+} + #bar{#Lambda} mass vs pT vs cent", kTH3D, {invMassAxis, ptAxis, centAxis});

        histos.add("xi1820/kminus_lambda/hMassPtCentKminusLambda", "K^{-} + #Lambda mass vs pT vs cent", kTH3D, {invMassAxis, ptAxis, centAxis});
        histos.add("xi1820/kminus_lambda/hMassPtCentKminusLambda_Mix", "Mixed event K^{-} + #Lambda mass vs pT vs cent", kTH3D, {invMassAxis, ptAxis, centAxis});
        histos.add("xi1820/kminus_lambda/hMassPtCentKminusLambda_Rot", "Rotated background K^{-} + #Lambda mass vs pT vs cent", kTH3D, {invMassAxis, ptAxis, centAxis});

        histos.add("xi1820/kminus_antilambda/hMassPtCentKminusAntiLambda", "K^{-} + #bar{#Lambda} mass vs pT vs cent", kTH3D, {invMassAxis, ptAxis, centAxis});
        histos.add("xi1820/kminus_antilambda/hMassPtCentKminusAntiLambda_Mix", "Mixed event K^{-} + #bar{#Lambda} mass vs pT vs cent", kTH3D, {invMassAxis, ptAxis, centAxis});
      } else {
        histos.add("xi1820/kplus_lambda/hMassPtCentDelKplusLambda", "K^{+} + #Lambda mass vs pT vs cent vs #DeltaR", kTHnSparseD, {invMassAxis, ptAxis, centAxis, axisDel});
        histos.add("xi1820/kplus_lambda/hMassPtCentDelKplusLambda_Mix", "Mixed event K^{+} + #Lambda mass vs pT vs cent vs #DeltaR", kTHnSparseD, {invMassAxis, ptAxis, centAxis, axisDel});

        // K+ Anti-Lambda
        histos.add("xi1820/kplus_antilambda/hMassPtCentDelKplusAntiLambda", "K^{+} + #bar{#Lambda} mass vs pT vs cent vs #DeltaR", kTHnSparseD, {invMassAxis, ptAxis, centAxis, axisDel});
        histos.add("xi1820/kplus_antilambda/hMassPtCentDelKplusAntiLambda_Mix", "Mixed event K^{+} + #bar{#Lambda} mass vs pT vs cent vs #DeltaR", kTHnSparseD, {invMassAxis, ptAxis, centAxis, axisDel});
        histos.add("xi1820/kplus_antilambda/hMassPtCentDelKplusAntiLambda_Rot", "Rotated background K^{+} + #bar{#Lambda} mass vs pT vs cent vs #DeltaR", kTHnSparseD, {invMassAxis, ptAxis, centAxis, axisDel});

        // K- Lambda
        histos.add("xi1820/kminus_lambda/hMassPtCentDelKminusLambda", "K^{-} + #Lambda mass vs pT vs cent vs #DeltaR", kTHnSparseD, {invMassAxis, ptAxis, centAxis, axisDel});
        histos.add("xi1820/kminus_lambda/hMassPtCentDelKminusLambda_Mix", "Mixed event K^{-} + #Lambda mass vs pT vs cent vs #DeltaR", kTHnSparseD, {invMassAxis, ptAxis, centAxis, axisDel});
        histos.add("xi1820/kminus_lambda/hMassPtCentDelKminusLambda_Rot", "Rotated background K^{-} + #Lambda mass vs pT vs cent vs #DeltaR", kTHnSparseD, {invMassAxis, ptAxis, centAxis, axisDel});

        // K- Anti-Lambda
        histos.add("xi1820/kminus_antilambda/hMassPtCentDelKminusAntiLambda", "K^{-} + #bar{#Lambda} mass vs pT vs cent vs #DeltaR", kTHnSparseD, {invMassAxis, ptAxis, centAxis, axisDel});
        histos.add("xi1820/kminus_antilambda/hMassPtCentDelKminusAntiLambda_Mix", "Mixed event K^{-} + #bar{#Lambda} mass vs pT vs cent vs #DeltaR", kTHnSparseD, {invMassAxis, ptAxis, centAxis, axisDel});
      }
    }

    // MC Reco histograms for charged K + Lambda channel
    if (doprocessMCWithTracks || doprocessMCWithMicroTracks) {
      histos.add("MC/kplus_antilambda/hMCRecoInvMassKplusAntiLambda", "Invariant mass of Xi(1820) to K^{-} + #Lambda (MC Reco)", kTH1F, {invMassAxis});
      histos.add("MC/kplus_antilambda/hMCRecoMassPtCentKplusAntiLambda", "Xi(1820) mass vs pT vs cent (K^{-} + #Lambda) (MC Reco)", kTHnSparseD, {invMassAxis, ptAxis, centAxis, ptAxis});

      histos.add("MC/kminus_lambda/hMCRecoInvMassKminusLambda", "Invariant mass of Xi(1820) to K^{-} + #Lambda (MC Reco)", kTH1F, {invMassAxis});
      histos.add("MC/kminus_lambda/hMCRecoMassPtCentKminusLambda", "Xi(1820) mass vs pT vs cent (K^{-} + #Lambda) (MC Reco)", kTHnSparseD, {invMassAxis, ptAxis, centAxis, ptAxis});
    }

    // K0s QA histograms
    if (cfgFillQA && (doprocessK0sLambda || doprocessMCK0sLambda)) {
      histos.add("QAbefore/k0sDCAtoPV", "K0s DCA to PV before cuts", kTH2F, {ptAxisQA, dcaAxis});
      histos.add("QAbefore/k0sMass", "K0s mass before cuts", kTH1F, {{100, 0.4, 0.6, "K^{0}_{S} mass (GeV/#it{c}^{2})"}});
      histos.add("QAbefore/k0sPt", "K0s pT before cuts", kTH1F, {ptAxisQA});
      histos.add("QAbefore/k0sEta", "K0s eta before cuts", kTH1F, {{100, -2.0, 2.0, "#eta"}});
      histos.add("QAbefore/k0sCosPA", "K0s CosPA before cuts", kTH2F, {ptAxisQA, cosPAAxis});
      histos.add("QAbefore/k0sRadius", "K0s radius before cuts", kTH2F, {ptAxisQA, radiusAxis});
      histos.add("QAbefore/k0sDauDCA", "K0s daughter DCA before cuts", kTH2F, {ptAxisQA, dcaAxis});
      histos.add("QAbefore/k0sProperLifetime", "K0s proper lifetime before cuts", kTH2F, {ptAxisQA, lifetimeAxis});
      histos.add("QAbefore/k0sDauTPCNsigmaPosPi", "K0s daughter pion TPC NSigma before cuts", kTH2F, {ptAxisQA, nsigmaAxis});
      histos.add("QAbefore/k0sDauTPCNsigmaNegPi", "K0s daughter pion TPC NSigma before cuts", kTH2F, {ptAxisQA, nsigmaAxis});
      histos.add("QAbefore/k0sNCrossedRowsPos", "K0s positive daughter crossed rows before cuts", kTH2F, {ptAxisQA, {160, 0, 160, "TPC crossed rows"}});
      histos.add("QAbefore/k0sNCrossedRowsNeg", "K0s negative daughter crossed rows before cuts", kTH2F, {ptAxisQA, {160, 0, 160, "TPC crossed rows"}});
      histos.add("QAbefore/k0sArmenterosPodolanski", "K0s candidate Armenteros-Podolanski before cuts", kTH3F, {armenterosAlphaAxis, armenterosQtAxis, ptAxisQA});

      histos.add("QAafter/k0sDCAtoPV", "K0s DCA to PV after cuts", kTH2F, {ptAxisQA, dcaAxis});
      histos.add("QAafter/k0sMass", "K0s mass after cuts", kTH1F, {{100, 0.4, 0.6, "K^{0}_{S} mass (GeV/#it{c}^{2})"}});
      histos.add("QAafter/k0sPt", "K0s pT after cuts", kTH1F, {ptAxisQA});
      histos.add("QAafter/k0sEta", "K0s eta after cuts", kTH1F, {{100, -2.0, 2.0, "#eta"}});
      histos.add("QAafter/k0sCosPA", "K0s CosPA after cuts", kTH2F, {ptAxisQA, cosPAAxis});
      histos.add("QAafter/k0sRadius", "K0s radius after cuts", kTH2F, {ptAxisQA, radiusAxis});
      histos.add("QAafter/k0sDauDCA", "K0s daughter DCA after cuts", kTH2F, {ptAxisQA, dcaAxis});
      histos.add("QAafter/k0sProperLifetime", "K0s proper lifetime after cuts", kTH2F, {ptAxisQA, lifetimeAxis});
      histos.add("QAafter/k0sDauTPCNsigmaPosPi", "K0s daughter pion TPC NSigma after cuts", kTH2F, {ptAxisQA, nsigmaAxis});
      histos.add("QAafter/k0sDauTPCNsigmaNegPi", "K0s daughter pion TPC NSigma after cuts", kTH2F, {ptAxisQA, nsigmaAxis});
      histos.add("QAafter/k0sNCrossedRowsPos", "K0s positive daughter crossed rows after cuts", kTH2F, {ptAxisQA, {160, 0, 160, "TPC crossed rows"}});
      histos.add("QAafter/k0sNCrossedRowsNeg", "K0s negative daughter crossed rows after cuts", kTH2F, {ptAxisQA, {160, 0, 160, "TPC crossed rows"}});
      histos.add("QAafter/k0sArmenterosPodolanski", "K0s candidate Armenteros-Podolanski after cuts", kTH3F, {armenterosAlphaAxis, armenterosQtAxis, ptAxisQA});
    }

    // K0s + Lambda
    if (doprocessK0sLambda || doprocessK0sLambdaMixedEvent || doprocessK0sLambdaMEDF || doprocessMCK0sLambda) {

      histos.add("xi1820/k0s_lambda/hInvMassK0sLambda", "Invariant mass of Xi(1820) to K^{0}_{S} + #Lambda", kTH1F, {invMassAxis});
      histos.add("xi1820/k0s_lambda/hInvMassK0sLambda_Mix", "Mixed event Invariant mass of Xi(1820) to K^{0}_{S} + #Lambda", kTH1F, {invMassAxis});
      histos.add("xi1820/k0s_antilambda/hInvMassK0sAntiLambda", "Invariant mass of Xi(1820) to K^{0}_{S} + #bar{#Lambda}", kTH1F, {invMassAxis});
      histos.add("xi1820/k0s_antilambda/hInvMassK0sAntiLambda_Mix", "Mixed event Invariant mass of Xi(1820) to K^{0}_{S} + #bar{#Lambda}", kTH1F, {invMassAxis});

      if (!(additionalConfig.cfgDelCheck)) {
        histos.add("xi1820/k0s_lambda/hMassPtCentK0sLambda", "Xi(1820) mass vs pT vs cent (K^{0}_{S}#Lambda)", kTH3D, {invMassAxis, ptAxis, centAxis});
        histos.add("xi1820/k0s_lambda/hMassPtCentK0sLambda_Mix", "Mixed event Xi(1820) mass vs pT vs cent (K^{0}_{S}#Lambda)", kTH3D, {invMassAxis, ptAxis, centAxis});
        histos.add("xi1820/k0s_lambda/hMassPtCentK0sLambda_Rot", "Rotated background Xi(1820) mass vs pT vs cent (K^{0}_{S}#Lambda)", kTH3D, {invMassAxis, ptAxis, centAxis});

        // K0s + Anti-Lambda
        histos.add("xi1820/k0s_antilambda/hMassPtCentK0sAntiLambda", "Xi(1820) mass vs pT vs cent (K^{0}_{S}#bar{#Lambda})", kTH3D, {invMassAxis, ptAxis, centAxis});
        histos.add("xi1820/k0s_antilambda/hMassPtCentK0sAntiLambda_Mix", "Mixed event Xi(1820) mass vs pT vs cent (K^{0}_{S}#bar{#Lambda})", kTH3D, {invMassAxis, ptAxis, centAxis});
        histos.add("xi1820/k0s_antilambda/hMassPtCentK0sAntiLambda_Rot", "Rotated background Xi(1820) mass vs pT vs cent (K^{0}_{S}#bar{#Lambda})", kTH3D, {invMassAxis, ptAxis, centAxis});
      } else {
        histos.add("xi1820/k0s_lambda/hMassPtCentDelK0sLambda", "Xi(1820) mass vs pT vs cent vs #DeltaR (K^{0}_{S}#Lambda)", kTHnSparseD, {{invMassAxis, ptAxis, centAxis, axisDel}});
        histos.add("xi1820/k0s_lambda/hMassPtCentDelK0sLambda_Mix", "Mixed event Xi(1820) mass vs pT vs cent vs #DeltaR (K^{0}_{S}#Lambda)", kTHnSparseD, {{invMassAxis, ptAxis, centAxis, axisDel}});
        histos.add("xi1820/k0s_lambda/hMassPtCentDelK0sLambda_Rot", "Rotated background Xi(1820) mass vs pT vs cent vs #DeltaR (K^{0}_{S}#Lambda)", kTHnSparseD, {{invMassAxis, ptAxis, centAxis, axisDel}});

        // K0s + Anti-Lambda
        histos.add("xi1820/k0s_antilambda/hMassPtCentDelK0sAntiLambda", "Xi(1820) mass vs pT vs cent vs #DeltaR (K^{0}_{S}#bar{#Lambda})", kTHnSparseD, {{invMassAxis, ptAxis, centAxis, axisDel}});
        histos.add("xi1820/k0s_antilambda/hMassPtCentDelK0sAntiLambda_Mix", "Mixed event Xi(1820) mass vs pT vs cent vs #DeltaR (K^{0}_{S}#bar{#Lambda})", kTHnSparseD, {{invMassAxis, ptAxis, centAxis, axisDel}});
        histos.add("xi1820/k0s_antilambda/hMassPtCentDelK0sAntiLambda_Rot", "Rotated background Xi(1820) mass vs pT vs cent vs #DeltaR (K^{0}_{S}#bar{#Lambda})", kTHnSparseD, {{invMassAxis, ptAxis, centAxis, axisDel}});
      }
    }

    if (doprocessMCK0sLambda) {
      histos.add("MC/k0s_lambda/hMCRecoInvMassK0sLambda", "Invariant mass of Xi(1820) to K^{0}_{S} + #Lambda (MC Reco)", kTH1F, {invMassAxis});
      histos.add("MC/k0s_lambda/hMCRecoMassPtCentK0sLambda", "Xi(1820) mass vs pT vs cent (K^{0}_{S}#Lambda) (MC Reco)", kTHnSparseD, {invMassAxis, ptAxis, centAxis, ptAxis});

      histos.add("MC/k0s_antilambda/hMCRecoInvMassK0sAntiLambda", "Invariant mass of Xi(1820) to K^{0}_{S} + #bar{#Lambda} (MC Reco)", kTH1F, {invMassAxis});
      histos.add("MC/k0s_antilambda/hMCRecoMassPtCentK0sAntiLambda", "Xi(1820) mass vs pT vs cent (K^{0}_{S}#bar{#Lambda}) (MC Reco)", kTHnSparseD, {invMassAxis, ptAxis, centAxis, ptAxis});
    }

    if (cfgFillTruthQA && (doprocessMCWithTracks || doprocessMCWithMicroTracks || doprocessMCK0sLambda)) {
      histos.add("QAMCTrue/lambdaPt", "Truth-matched Lambda pT per signal pair", kTH1F, {ptAxisQA});
      histos.add("QAMCTrue/lambdaDaughterTPCNSigmaPi", "Truth-matched Lambda pion TPC PID", kTH2F, {ptAxisQA, nsigmaAxis});
      histos.add("QAMCTrue/lambdaDaughterTPCNSigmaPr", "Truth-matched Lambda proton TPC PID", kTH2F, {ptAxisQA, nsigmaAxis});
      if (doprocessMCWithTracks || doprocessMCWithMicroTracks) {
        histos.add("QAMCTrue/kaonPt", "Truth-matched kaon pT per signal pair", kTH1F, {ptAxisQA});
        histos.add("QAMCTrue/kaonDCAxy", "Truth-matched kaon DCAxy", kTH2F, {ptAxisQA, dcaxyAxis});
        histos.add("QAMCTrue/kaonDCAz", "Truth-matched kaon DCAz", kTH2F, {ptAxisQA, dcazAxis});
        histos.add("QAMCTrue/kaonTPCNSigma", "Truth-matched kaon TPC PID", kTH2F, {ptAxisQA, nsigmaAxis});
        histos.add("QAMCTrue/kaonTOFNSigma", "Truth-matched kaon TOF PID", kTH2F, {ptAxisQA, nsigmaAxis});
      }
      if (doprocessMCK0sLambda) {
        histos.add("QAMCTrue/k0sPt", "Truth-matched K0s pT per signal pair", kTH1F, {ptAxisQA});
      }
    }

    if (doprocessMCGen) {
      histos.add("multQA/h2MultCentMC", "Multiplicity vs Centrality MC", HistType::kTH2D, {centAxis, additionalConfig.multNTracksAxis});
      // MC truth invariant mass vs pT (2D)
      histos.add("MC/hMCGenPtCentMultKminusLambda", "MC Truth Mass vs pT K^{-}#Lambda", kTH3D, {ptAxis, centAxis, additionalConfig.multNTracksAxis});
      histos.add("MC/hMCGenPtCentMultKplusAntiLambda", "MC Truth Mass vs pT K^{+}#bar{#Lambda}", kTH3D, {ptAxis, centAxis, additionalConfig.multNTracksAxis});
      histos.add("MC/hMCGenPtCentMultK0sLambda", "MC Truth Mass vs pT K^{0}_{S}#Lambda", kTH3D, {ptAxis, centAxis, additionalConfig.multNTracksAxis});
      histos.add("MC/hMCGenPtCentMultK0sAntiLambda", "MC Truth Mass vs pT K^{0}_{S}#bar{#Lambda}", kTH3D, {ptAxis, centAxis, additionalConfig.multNTracksAxis});
    }

    // MC truth histograms
    AxisSpec etaAxis = {100, -2.0, 2.0, "#eta"};
    AxisSpec rapidityAxis = {100, -2.0, 2.0, "y"};

    if (doprocessMCTruth) {
      histos.add("MC/hMCTruthXi1820Pt", "MC Generated Xi(1820) pT", kTH1D, {ptAxis});
      histos.add("MC/hMCTruthXi1820PtEta", "MC Generated Xi(1820) pT vs eta", kTH2D, {ptAxis, etaAxis});
      histos.add("MC/hMCTruthXi1820Y", "MC Generated Xi(1820) rapidity", kTH1D, {rapidityAxis});

      // MC truth invariant mass (from MC particles)
      histos.add("MC/hMCTruthInvMassKminusLambda", "MC Truth Inv Mass K^{-}#Lambda", kTH1D, {invMassAxis});
      histos.add("MC/hMCTruthInvMassKplusAntiLambda", "MC Truth Inv Mass K^{+}#bar{#Lambda}", kTH1D, {invMassAxis});
      histos.add("MC/hMCTruthInvMassK0sLambda", "MC Truth Inv Mass K^{0}_{S}#Lambda", kTH1D, {invMassAxis});
      histos.add("MC/hMCTruthInvMassK0sAntiLambda", "MC Truth Inv Mass K^{0}_{S}#bar{#Lambda}", kTH1D, {invMassAxis});

      // MC truth invariant mass vs pT (2D)
      histos.add("MC/hMCTruthMassPtKminusLambda", "MC Truth Mass vs pT K^{-}#Lambda", kTH2D, {invMassAxis, ptAxis});
      histos.add("MC/hMCTruthMassPtKplusAntiLambda", "MC Truth Mass vs pT K^{+}#bar{#Lambda}", kTH2D, {invMassAxis, ptAxis});
      histos.add("MC/hMCTruthMassPtK0sLambda", "MC Truth Mass vs pT K^{0}_{S}#Lambda", kTH2D, {invMassAxis, ptAxis});
      histos.add("MC/hMCTruthMassPtK0sAntiLambda", "MC Truth Mass vs pT K^{0}_{S}#bar{#Lambda}", kTH2D, {invMassAxis, ptAxis});
    }
  }
  template <typename V0Type>
  static std::array<int, 2> v0DaughterIds(const V0Type& v0)
  {
    const auto& indices = v0.indices();
    return {indices[0], indices[1]};
  }

  template <std::size_t N, std::size_t M>
  static bool sharesAnyDaughterId(const std::array<int, N>& first, const std::array<int, M>& second)
  {
    for (const auto& firstId : first) {
      for (const auto& secondId : second) {
        if (firstId == secondId) {
          return true;
        }
      }
    }
    return false;
  }

  template <typename FirstVectorType, typename SecondVectorType>
  static float deltaR(const FirstVectorType& first, const SecondVectorType& second)
  {
    const auto dPhi = RecoDecay::constrainAngle(first.Phi() - second.Phi(), -constants::math::PI);
    const auto dEta = first.Eta() - second.Eta();
    return std::sqrt(dPhi * dPhi + dEta * dEta);
  }

  bool passesRapidity(float rapidity) const
  {
    return -additionalConfig.cfgRapidityCut.value < rapidity && rapidity < additionalConfig.cfgRapidityCut.value;
  }

  template <typename V0Type>
  void fillTruthLambdaQA(const V0Type& lambda, bool isLambda)
  {
    if (!cfgFillTruthQA) {
      return;
    }
    histos.fill(HIST("QAMCTrue/lambdaPt"), lambda.pt());
    histos.fill(HIST("QAMCTrue/lambdaDaughterTPCNSigmaPi"), lambda.pt(), isLambda ? lambda.daughterTPCNSigmaNegPi() : lambda.daughterTPCNSigmaPosPi());
    histos.fill(HIST("QAMCTrue/lambdaDaughterTPCNSigmaPr"), lambda.pt(), isLambda ? lambda.daughterTPCNSigmaPosPr() : lambda.daughterTPCNSigmaNegPr());
  }

  template <typename TrackType>
  void fillTruthKaonQA(const TrackType& track)
  {
    if (!cfgFillTruthQA) {
      return;
    }
    histos.fill(HIST("QAMCTrue/kaonPt"), track.pt());
    histos.fill(HIST("QAMCTrue/kaonDCAxy"), track.pt(), track.dcaXY());
    histos.fill(HIST("QAMCTrue/kaonDCAz"), track.pt(), track.dcaZ());
    histos.fill(HIST("QAMCTrue/kaonTPCNSigma"), track.pt(), track.tpcNSigmaKa());
    if (track.hasTOF() && !std::isnan(track.tofNSigmaKa())) {
      histos.fill(HIST("QAMCTrue/kaonTOFNSigma"), track.pt(), track.tofNSigmaKa());
    }
  }

  // Lambda/Anti-Lambda selection
  template <typename CollisionType, typename V0Type>
  bool v0Cut(const CollisionType& collision, const V0Type& v0, bool isLambda)
  {
    // Basic kinematic cuts
    if (std::abs(v0.eta()) >= cMaxV0Etacut) {
      return false;
    }
    if (v0.pt() <= cMinPtcut) {
      return false;
    }

    // DCA to PV
    if (std::abs(v0.dcav0topv()) >= additionalConfig.cMaxDcaToPVV0Lambda) {
      return false;
    }

    // Topological cuts
    if (v0.v0CosPA() <= cV0MinCosPA) {
      return false;
    }
    if (v0.daughDCA() >= cV0MaxDaughDCA) {
      return false;
    }

    // Daughter DCA to PV cuts
    if (std::abs(v0.dcapostopv()) <= cV0DauPosDCAtoPVMin) {
      return false;
    }
    if (std::abs(v0.dcanegtopv()) <= cV0DauNegDCAtoPVMin) {
      return false;
    }

    // Radius cuts
    float radius = v0.transRadius();
    if (radius <= cV0RadiusMin || radius >= cV0RadiusMax) {
      return false;
    }

    // Proper lifetime cut
    float dx = v0.decayVtxX() - collision.posX();
    float dy = v0.decayVtxY() - collision.posY();
    float dz = v0.decayVtxZ() - collision.posZ();
    float l = std::sqrt(dx * dx + dy * dy + dz * dz);
    float p = std::sqrt(v0.px() * v0.px() + v0.py() * v0.py() + v0.pz() * v0.pz());
    auto properLifetime = (l / (p + SmallMomentumDenominator)) * MassLambda;
    if (properLifetime >= cV0ProperLifetimeMax) {
      return false;
    }

    // Mass window
    if (isLambda) {
      if (std::abs(v0.mLambda() - MassLambda) >= cV0MassWindow) {
        return false;
      }
      // TPC Nsigma cuts for daughter tracks
      if (!passesPIDWindow(v0.daughterTPCNSigmaNegPi(), cLambdaDaughterPiTPCNSigmaMin, cLambdaDaughterPiTPCNSigmaMax, cLambdaDaughterPiTPCNSigmaMean)) {
        return false;
      }
      if (!passesPIDWindow(v0.daughterTPCNSigmaPosPr(), cLambdaDaughterPrTPCNSigmaMin, cLambdaDaughterPrTPCNSigmaMax, cLambdaDaughterPrTPCNSigmaMean)) {
        return false;
      }
      // Note: TOF beta cut that calibrated by primary vertex is not applied for V0 daughters in this moment
    } else {
      if (std::abs(v0.mAntiLambda() - MassLambda) >= cV0MassWindow) {
        return false;
      }
      if (!passesPIDWindow(v0.daughterTPCNSigmaPosPi(), cLambdaDaughterPiTPCNSigmaMin, cLambdaDaughterPiTPCNSigmaMax, cLambdaDaughterPiTPCNSigmaMean)) {
        return false;
      }
      if (!passesPIDWindow(v0.daughterTPCNSigmaNegPr(), cLambdaDaughterPrTPCNSigmaMin, cLambdaDaughterPrTPCNSigmaMax, cLambdaDaughterPrTPCNSigmaMean)) {
        return false;
      }
      // Note: TOF beta cut that calibrated by primary vertex is not applied for V0 daughters in this moment
    }

    if (cV0sCrossMassRejection) {
      if (std::abs(v0.mK0Short() - MassK0Short) <= cV0sCrossMassRejectionWindow) {
        return false;
      }
    }

    if (v0.nCrossedRowsPos() <= cLambdaPosDaughterMinCrossedRows) {
      return false;
    }
    if (v0.nCrossedRowsNeg() <= cLambdaNegDaughterMinCrossedRows) {
      return false;
    }

    if (v0.qtarm() >= cK0sArmenterosAlphaCoeff * std::fabs(v0.alpha())) {
      return false;
    }
    return true;
  }

  // K0s selection
  template <typename CollisionType, typename V0Type>
  bool k0sCut(const CollisionType& collision, const V0Type& v0)
  {
    // Basic kinematic cuts
    if (std::abs(v0.eta()) >= cMaxV0Etacut) {
      return false;
    }
    if (v0.pt() <= cMinPtcut) {
      return false;
    }

    // Topological cuts
    if (v0.v0CosPA() <= cK0sMinCosPA) {
      return false;
    }
    if (v0.daughDCA() >= cK0sMaxDaughDCA) {
      return false;
    }

    // Daughter DCA to PV cuts
    if (std::abs(v0.dcapostopv()) <= cK0sDauPosDCAtoPVMin) {
      return false;
    }
    if (std::abs(v0.dcanegtopv()) <= cK0sDauNegDCAtoPVMin) {
      return false;
    }

    // Radius cuts
    float radius = v0.transRadius();
    if (radius <= cK0sRadiusMin || radius >= cK0sRadiusMax) {
      return false;
    }

    // DCA to PV
    if (std::abs(v0.dcav0topv()) >= additionalConfig.cMaxDcaToPVV0K0s) {
      return false;
    }

    // Proper lifetime cut
    float dx = v0.decayVtxX() - collision.posX();
    float dy = v0.decayVtxY() - collision.posY();
    float dz = v0.decayVtxZ() - collision.posZ();
    float l = std::sqrt(dx * dx + dy * dy + dz * dz);
    float p = std::sqrt(v0.px() * v0.px() + v0.py() * v0.py() + v0.pz() * v0.pz());
    auto properLifetime = (l / (p + SmallMomentumDenominator)) * MassK0Short;
    if (properLifetime >= cK0sProperLifetimeMax) {
      return false;
    }

    // Mass window
    if (std::abs(v0.mK0Short() - MassK0Short) >= cK0sMassWindow) {
      return false;
    }

    // Competing V0 rejection: remove (Anti)Λ
    if (cK0sCrossMassRejection) {
      if (std::abs(v0.mLambda() - MassLambda) <= cK0sCrossMassRejectionWindow) {
        return false;
      }
      if (std::abs(v0.mAntiLambda() - MassLambda) <= cK0sCrossMassRejectionWindow) {
        return false;
      }
    }

    if (!passesPIDWindow(v0.daughterTPCNSigmaPosPi(), cK0sDaughterPiTPCNSigmaMin, cK0sDaughterPiTPCNSigmaMax, cK0sDaughterPiTPCNSigmaMean)) {
      return false;
    }
    if (!passesPIDWindow(v0.daughterTPCNSigmaNegPi(), cK0sDaughterPiTPCNSigmaMin, cK0sDaughterPiTPCNSigmaMax, cK0sDaughterPiTPCNSigmaMean)) {
      return false;
    }
    // Note: TOF beta cut that calibrated by primary vertex is not applied for V0 daughters in this moment

    if (v0.nCrossedRowsPos() <= cK0sPosDaughterMinCrossedRows) {
      return false;
    }
    if (v0.nCrossedRowsNeg() <= cK0sNegDaughterMinCrossedRows) {
      return false;
    }

    if (v0.qtarm() <= cK0sArmenterosAlphaCoeff * std::fabs(v0.alpha())) {
      return false;
    }

    return true;
  }

  // Strict signed window, retaining NaN rejection without negated compound comparisons.
  static bool passesPIDWindow(float nSigma, float minimum, float maximum, float mean)
  {
    const float centered = nSigma - mean;
    return minimum < centered && centered < maximum;
  }

  // Preserve pT-bin membership [low, high); PID window edges themselves are strict.
  int getPtBinIndex(float pt)
  {
    const auto& ptBins = cKaonPIDPtBins.value;
    for (std::size_t i = 1; i < ptBins.size(); ++i) {
      if (pt >= ptBins[i - 1] && pt < ptBins[i]) {
        return static_cast<int>(i - 1);
      }
    }
    return -1;
  }

  template <typename TrackType>
  bool kaonPidCut(const TrackType& track)
  {
    const bool useTOF = !cKaonBypassTOF && track.hasTOF();
    float tpcMin = cKaonTPCNSigmaMin;
    float tpcMax = cKaonTPCNSigmaMax;
    float tofMin = cKaonTOFNSigmaMin;
    float tofMax = cKaonTOFNSigmaMax;
    if (cKaonUsePtDepPID) {
      const int ptBin = getPtBinIndex(track.pt());
      if (ptBin < 0) {
        return false;
      }
      const auto bin = static_cast<std::size_t>(ptBin);
      if (bin >= cKaonTPCNSigmaCuts.value.size() || bin >= cKaonTPCNSigmaMinCuts.value.size()) {
        return false;
      }
      tpcMin = cKaonTPCNSigmaMinCuts.value[bin];
      tpcMax = cKaonTPCNSigmaCuts.value[bin];
      if (useTOF) {
        if (bin >= cKaonTOFNSigmaCuts.value.size() || bin >= cKaonTOFNSigmaMinCuts.value.size()) {
          return false;
        }
        tofMin = cKaonTOFNSigmaMinCuts.value[bin];
        tofMax = cKaonTOFNSigmaCuts.value[bin];
      }
    }
    if (!passesPIDWindow(track.tpcNSigmaKa(), tpcMin, tpcMax, cKaonTPCNSigmaMean)) {
      return false;
    }
    // No TOF-presence veto: missing or explicitly bypassed TOF uses TPC alone.
    return !useTOF ||
           passesPIDWindow(track.tofNSigmaKa(), tofMin, tofMax, cKaonTOFNSigmaMean);
  }

  // Kaon track selection (for both ResoTracks and ResoMicroTracks)
  template <bool IsResoMicrotrack, typename TrackType>
  bool kaonCut(const TrackType& track)
  {
    const float candPt = track.pt();
    if (!std::isfinite(candPt) || !std::isfinite(track.eta()) ||
        !std::isfinite(track.dcaXY()) || !std::isfinite(track.dcaZ())) {
      return false;
    }
    if (candPt <= cKaonPtMin || std::abs(track.eta()) >= cKaonEtaMax) {
      return false;
    }
    const double dcaPtCut = additionalConfig.cUsePtDepDCAForKaons
                              ? additionalConfig.cDCAToPVByPtFirstP0.value + additionalConfig.cDCAToPVByPtFirstExp.value * std::pow(candPt, -static_cast<double>(additionalConfig.cDCAToPVByPtFirstPower.value))
                              : 0.;
    const double dcaXYCut = additionalConfig.cUsePtDepDCAForKaons ? dcaPtCut : cKaonDCAxyMax.value;
    const double dcaZCut = additionalConfig.cUsePtDepDCAForKaons ? dcaPtCut : cKaonDCAzMax.value;
    if (std::abs(track.dcaXY()) >= dcaXYCut || std::abs(track.dcaZ()) >= dcaZCut) {
      return false;
    }

    // Track quality cuts - check if fields are available (only for ResoTracks)
    if constexpr (!IsResoMicrotrack) {
      if constexpr (requires { track.tpcNClsFound(); }) {
        if (track.tpcNClsFound() <= cKaonTPCNClusMin) {
          return false;
        }
      }
      if constexpr (requires { track.itsNCls(); }) {
        if (track.itsNCls() <= cKaonITSNClusMin) {
          return false;
        }
      }
    }

    // PID selection
    if (!kaonPidCut(track)) {
      return false;
    }

    // Flag selections for Primary track selection
    if (additionalConfig.cfgPVContributor && !track.isPVContributor()) {
      return false;
    }
    if (additionalConfig.cfgPrimaryTrack && !track.isPrimaryTrack()) {
      return false;
    }

    return true;
  }

  bool hasChargedSameEventProcess()
  {
    return doprocessDataWithTracks || doprocessDataWithMicroTracks || doprocessMCWithTracks || doprocessMCWithMicroTracks;
  }

  // Select each V0 once; QA and event counts must not be weighted by kaon multiplicity.
  template <bool IsMix, typename CollisionT, typename V0sT>
  auto selectLambdas(const CollisionT& collision, const V0sT& v0s, bool fillSharedQA = true)
  {
    std::vector<SelectedLambda<typename V0sT::iterator>> selected;
    selected.reserve(v0s.size());
    for (const auto& v0 : v0s) {
      // Lambda QA before cuts
      if (!IsMix && cfgFillQA && fillSharedQA) {
        histos.fill(HIST("QAbefore/lambdaDCAtoPV"), v0.pt(), v0.dcav0topv());
        histos.fill(HIST("QAbefore/lambdaMass"), v0.mLambda());
        histos.fill(HIST("QAbefore/lambdaMassAnti"), v0.mAntiLambda());
        histos.fill(HIST("QAbefore/lambdaPt"), v0.pt());
        histos.fill(HIST("QAbefore/lambdaEta"), v0.eta());
        histos.fill(HIST("QAbefore/lambdaCosPA"), v0.pt(), v0.v0CosPA());
        histos.fill(HIST("QAbefore/lambdaRadius"), v0.pt(), v0.transRadius());
        histos.fill(HIST("QAbefore/lambdaDauDCA"), v0.pt(), v0.daughDCA());
        histos.fill(HIST("QAbefore/lambdaDauPosDCA"), v0.pt(), std::abs(v0.dcapostopv()));
        histos.fill(HIST("QAbefore/lambdaDauNegDCA"), v0.pt(), std::abs(v0.dcanegtopv()));
        histos.fill(HIST("QAbefore/lambdaDaughterTPCNSigmaPosPr"), v0.pt(), v0.daughterTPCNSigmaPosPr());
        histos.fill(HIST("QAbefore/lambdaDaughterTPCNSigmaNegPi"), v0.pt(), v0.daughterTPCNSigmaNegPi());
        histos.fill(HIST("QAbefore/lambdaAntiDaughterTPCNSigmaPosPi"), v0.pt(), v0.daughterTPCNSigmaPosPi());
        histos.fill(HIST("QAbefore/lambdaAntiDaughterTPCNSigmaNegPr"), v0.pt(), v0.daughterTPCNSigmaNegPr());
        histos.fill(HIST("QAbefore/lambdaNCrossedRowsPos"), v0.pt(), v0.nCrossedRowsPos());
        histos.fill(HIST("QAbefore/lambdaNCrossedRowsNeg"), v0.pt(), v0.nCrossedRowsNeg());

        // Calculate proper lifetime manually
        float dx = v0.decayVtxX() - collision.posX();
        float dy = v0.decayVtxY() - collision.posY();
        float dz = v0.decayVtxZ() - collision.posZ();
        float l = std::sqrt(dx * dx + dy * dy + dz * dz);
        float p = std::sqrt(v0.px() * v0.px() + v0.py() * v0.py() + v0.pz() * v0.pz());
        auto properLifetime = (l / (p + SmallMomentumDenominator)) * MassLambda;
        histos.fill(HIST("QAbefore/lambdaProperLifetime"), v0.pt(), properLifetime);
        histos.fill(HIST("QAbefore/lambdaArmenterosPodolanski"), v0.alpha(), v0.qtarm(), v0.pt());
      }

      // Try Lambda
      bool isLambda = v0Cut(collision, v0, true);
      // Try Anti-Lambda
      bool isAntiLambda = v0Cut(collision, v0, false);

      if (!isLambda && !isAntiLambda) {
        continue;
      }

      if (!IsMix && cfgFillQA && fillSharedQA) {
        // QA after cuts (fill for whichever passes)
        if (isLambda) {
          histos.fill(HIST("QAafter/lambdaMass"), v0.mLambda());
          histos.fill(HIST("QAafter/lambdaDaughterTPCNSigmaPosPr"), v0.pt(), v0.daughterTPCNSigmaPosPr());
          histos.fill(HIST("QAafter/lambdaDaughterTPCNSigmaNegPi"), v0.pt(), v0.daughterTPCNSigmaNegPi());
        }
        if (isAntiLambda) {
          histos.fill(HIST("QAafter/lambdaMassAnti"), v0.mAntiLambda());
          histos.fill(HIST("QAafter/lambdaAntiDaughterTPCNSigmaPosPi"), v0.pt(), v0.daughterTPCNSigmaPosPi());
          histos.fill(HIST("QAafter/lambdaAntiDaughterTPCNSigmaNegPr"), v0.pt(), v0.daughterTPCNSigmaNegPr());
        }
        histos.fill(HIST("QAafter/lambdaDCAtoPV"), v0.pt(), v0.dcav0topv());
        histos.fill(HIST("QAafter/lambdaPt"), v0.pt());
        histos.fill(HIST("QAafter/lambdaEta"), v0.eta());
        histos.fill(HIST("QAafter/lambdaCosPA"), v0.pt(), v0.v0CosPA());
        histos.fill(HIST("QAafter/lambdaRadius"), v0.pt(), v0.transRadius());
        histos.fill(HIST("QAafter/lambdaDauDCA"), v0.pt(), v0.daughDCA());
        histos.fill(HIST("QAafter/lambdaDauPosDCA"), v0.pt(), std::abs(v0.dcapostopv()));
        histos.fill(HIST("QAafter/lambdaDauNegDCA"), v0.pt(), std::abs(v0.dcanegtopv()));
        histos.fill(HIST("QAafter/lambdaNCrossedRowsPos"), v0.pt(), v0.nCrossedRowsPos());
        histos.fill(HIST("QAafter/lambdaNCrossedRowsNeg"), v0.pt(), v0.nCrossedRowsNeg());

        float dx = v0.decayVtxX() - collision.posX();
        float dy = v0.decayVtxY() - collision.posY();
        float dz = v0.decayVtxZ() - collision.posZ();
        float l = std::sqrt(dx * dx + dy * dy + dz * dz);
        float p = std::sqrt(v0.px() * v0.px() + v0.py() * v0.py() + v0.pz() * v0.pz());
        auto properLifetime = (l / (p + SmallMomentumDenominator)) * MassLambda;
        histos.fill(HIST("QAafter/lambdaProperLifetime"), v0.pt(), properLifetime);
        histos.fill(HIST("QAafter/lambdaArmenterosPodolanski"), v0.alpha(), v0.qtarm(), v0.pt());
      }

      selected.push_back({v0, isLambda, isAntiLambda});
    }
    return selected;
  }

  template <bool IsMix, bool IsResoMicrotrack, bool IsMC, typename CollisionT, typename V0sT, typename TracksT, typename LambdaCollisionT>
  void fillChargedKLambda(const CollisionT& collision, const V0sT& v0s, const TracksT& tracks, const LambdaCollisionT& lambdaCollision) // Xi(1820) analysis: charged K + Lambda channel
  {
    auto cent = collision.cent();

    // Fill event QA histograms (only for same-event)
    if (!IsMix && cfgFillEventQA) {
      histos.fill(HIST("Event/posZ"), collision.posZ());
      histos.fill(HIST("Event/centrality"), cent);
      histos.fill(HIST("Event/posZvsCent"), collision.posZ(), cent);
      histos.fill(HIST("Event/nV0s"), v0s.size());
      histos.fill(HIST("Event/nKaons"), tracks.size());
    }

    std::vector<typename TracksT::iterator> selectedKaons;
    selectedKaons.reserve(tracks.size());

    // Build 4 combinations
    ROOT::Math::PxPyPzEVector pKaon, pLambda, pRes, lDaughterRot, lResonanceRot;

    // Loop over kaon candidates
    for (const auto& kaon : tracks) {
      // QA before cuts
      if (!IsMix && cfgFillQA) {
        histos.fill(HIST("QAbefore/kaonPt"), kaon.pt());
        histos.fill(HIST("QAbefore/kaonEta"), kaon.eta());
        histos.fill(HIST("QAbefore/kaonDCAxy"), kaon.pt(), kaon.dcaXY());
        histos.fill(HIST("QAbefore/kaonDCAz"), kaon.pt(), kaon.dcaZ());
        histos.fill(HIST("QAbefore/kaonTPCNSigma"), kaon.pt(), kaon.tpcNSigmaKa());
        if (kaon.hasTOF() && !std::isnan(kaon.tofNSigmaKa())) {
          histos.fill(HIST("QAbefore/kaonTOFNSigma"), kaon.pt(), kaon.tofNSigmaKa());
        }
        if constexpr (!IsResoMicrotrack) {
          if constexpr (requires { kaon.tpcNClsFound(); }) {
            histos.fill(HIST("QAbefore/kaonTPCNcls"), kaon.tpcNClsFound());
          }
          if constexpr (requires { kaon.itsNCls(); }) {
            histos.fill(HIST("QAbefore/kaonITSNcls"), kaon.itsNCls());
          }
        }
      }

      if (!kaonCut<IsResoMicrotrack>(kaon)) {
        continue;
      }

      if (!IsMix && cfgFillQA) {
        // QA after cuts
        histos.fill(HIST("QAafter/kaonPt"), kaon.pt());
        histos.fill(HIST("QAafter/kaonEta"), kaon.eta());
        histos.fill(HIST("QAafter/kaonDCAxy"), kaon.pt(), kaon.dcaXY());
        histos.fill(HIST("QAafter/kaonDCAz"), kaon.pt(), kaon.dcaZ());
        histos.fill(HIST("QAafter/kaonTPCNSigma"), kaon.pt(), kaon.tpcNSigmaKa());
        if (kaon.hasTOF() && !std::isnan(kaon.tofNSigmaKa())) {
          histos.fill(HIST("QAafter/kaonTOFNSigma"), kaon.pt(), kaon.tofNSigmaKa());
        }
        if constexpr (!IsResoMicrotrack) {
          if constexpr (requires { kaon.tpcNClsFound(); }) {
            histos.fill(HIST("QAafter/kaonTPCNcls"), kaon.tpcNClsFound());
          }
          if constexpr (requires { kaon.itsNCls(); }) {
            histos.fill(HIST("QAafter/kaonITSNcls"), kaon.itsNCls());
          }
        }
      }

      selectedKaons.push_back(kaon);
    }

    const auto selectedLambdas = selectLambdas<IsMix>(lambdaCollision, v0s);
    for (const auto& kaon : selectedKaons) {
      const int kaonCharge = kaon.sign();
      for (const auto& selected : selectedLambdas) {
        const auto& v0 = selected.candidate;
        if constexpr (!IsMix) {
          const auto daughterIds = v0DaughterIds(v0);
          if (kaon.trackId() == daughterIds[0] || kaon.trackId() == daughterIds[1]) {
            continue;
          }
        }

        pKaon = ROOT::Math::PxPyPzEVector(ROOT::Math::PtEtaPhiMVector(kaon.pt(), kaon.eta(), kaon.phi(), MassKaonCharged));

        // A failed rapidity or truth match skips only this mass hypothesis.
        for (const auto& isLambda : std::array{true, false}) {
          if (isLambda ? !selected.isLambda : !selected.isAntiLambda) {
            continue;
          }

          // K+ + Lambda -> Bkg channel for charged Xi(1820)
          if (kaonCharge > 0 && isLambda) {
            pLambda = ROOT::Math::PxPyPzEVector(ROOT::Math::PtEtaPhiMVector(v0.pt(), v0.eta(), v0.phi(), v0.mLambda()));
            pRes = pKaon + pLambda;
            auto pCandRapidity = pRes.Rapidity();
            if (!passesRapidity(pCandRapidity)) { // skip candidate if reconstructed rapidity is outside of cut
              continue;
            }
            if constexpr (!IsMix) {
              histos.fill(HIST("xi1820/kplus_lambda/hInvMassKplusLambda"), pRes.M());
              if (!(additionalConfig.cfgDelCheck)) {
                histos.fill(HIST("xi1820/kplus_lambda/hMassPtCentKplusLambda"), pRes.M(), pRes.Pt(), cent);
              } else {
                auto pCandDel = deltaR(pLambda, pKaon);
                histos.fill(HIST("xi1820/kplus_lambda/hMassPtCentDelKplusLambda"), pRes.M(), pRes.Pt(), cent, pCandDel);
              }
            } else {
              histos.fill(HIST("xi1820/kplus_lambda/hInvMassKplusLambda_Mix"), pRes.M());
              if (!(additionalConfig.cfgDelCheck)) {
                histos.fill(HIST("xi1820/kplus_lambda/hMassPtCentKplusLambda_Mix"), pRes.M(), pRes.Pt(), cent);
              } else {
                auto pCandDel = deltaR(pLambda, pKaon);
                histos.fill(HIST("xi1820/kplus_lambda/hMassPtCentDelKplusLambda_Mix"), pRes.M(), pRes.Pt(), cent, pCandDel);
              }
            }
          } // End of Bkg

          // K+ + Anti-Lambda -> Signal channel for Anti-charged Xi(1820)
          if (kaonCharge > 0 && !isLambda) {
            pLambda = ROOT::Math::PxPyPzEVector(ROOT::Math::PtEtaPhiMVector(v0.pt(), v0.eta(), v0.phi(), v0.mAntiLambda()));
            pRes = pKaon + pLambda;
            auto pCandRapidity = pRes.Rapidity();
            if (!passesRapidity(pCandRapidity)) { // skip candidate if reconstructed rapidity is outside of cut
              continue;
            }
            if constexpr (!IsMix) {
              histos.fill(HIST("xi1820/kplus_antilambda/hInvMassKplusAntiLambda"), pRes.M());
              if (!(additionalConfig.cfgDelCheck)) {
                histos.fill(HIST("xi1820/kplus_antilambda/hMassPtCentKplusAntiLambda"), pRes.M(), pRes.Pt(), cent);
              } else {
                auto pCandDel = deltaR(pLambda, pKaon);
                histos.fill(HIST("xi1820/kplus_antilambda/hMassPtCentDelKplusAntiLambda"), pRes.M(), pRes.Pt(), cent, pCandDel);
              }
              if (additionalConfig.cfgFillRotBkg) {
                for (int i = 0; i < additionalConfig.cfgNrotBkg; i++) {
                  const auto lRotAngle = additionalConfig.cfgNrotBkg == 1
                                           ? 0.5f * (additionalConfig.cfgMinRot + additionalConfig.cfgMaxRot)
                                           : additionalConfig.cfgMinRot + i * ((additionalConfig.cfgMaxRot - additionalConfig.cfgMinRot) / (additionalConfig.cfgNrotBkg - 1));
                  if (cfgFillEventQA) {
                    histos.fill(HIST("Event/hRotBkg"), lRotAngle - o2::constants::math::PI);
                  }
                  if (additionalConfig.cfgRotKaon) {
                    lDaughterRot = pKaon;
                    ROOT::Math::RotationZ rot(lRotAngle);
                    auto p3 = rot * lDaughterRot.Vect();
                    lDaughterRot = LorentzVectorSetXYZM(p3.X(), p3.Y(), p3.Z(), lDaughterRot.M());
                    lResonanceRot = lDaughterRot + pLambda;
                  } else {
                    lDaughterRot = pLambda;
                    ROOT::Math::RotationZ rot(lRotAngle);
                    auto p3 = rot * lDaughterRot.Vect();
                    lDaughterRot = LorentzVectorSetXYZM(p3.X(), p3.Y(), p3.Z(), lDaughterRot.M());
                    lResonanceRot = pKaon + lDaughterRot;
                  }
                  if (!(additionalConfig.cfgDelCheck)) {
                    histos.fill(HIST("xi1820/kplus_antilambda/hMassPtCentKplusAntiLambda_Rot"), lResonanceRot.M(), lResonanceRot.Pt(), cent);
                  } else {
                    auto pRotDel = additionalConfig.cfgRotKaon ? deltaR(lDaughterRot, pLambda) : deltaR(lDaughterRot, pKaon);
                    histos.fill(HIST("xi1820/kplus_antilambda/hMassPtCentDelKplusAntiLambda_Rot"), lResonanceRot.M(), lResonanceRot.Pt(), cent, pRotDel);
                  }
                }
              }

              if constexpr (IsMC) { // Calculate Acceptance x efficiency for "the particle" channel
                if (v0.motherPDG() != -PdgChargedXi1820 || kaon.motherPDG() != v0.motherPDG()) {
                  continue;
                }
                if (kaon.pdgCode() != PDG_t::kKPlus || v0.pdgCode() != PDG_t::kLambda0Bar) {
                  continue;
                }
                if (v0.motherId() < 0 || kaon.motherId() != v0.motherId()) {
                  continue;
                }
                auto pMCPt = v0.motherPt();                                                  // Check particle's pT resolution
                if (additionalConfig.cUseTruthRapidity && !passesRapidity(v0.motherRap())) { // skip candidate if True rapidity of mother particle is outside of cut
                  continue;
                }
                histos.fill(HIST("MC/kplus_antilambda/hMCRecoInvMassKplusAntiLambda"), pRes.M());
                histos.fill(HIST("MC/kplus_antilambda/hMCRecoMassPtCentKplusAntiLambda"), pRes.M(), pRes.Pt(), cent, pMCPt);

                fillTruthKaonQA(kaon);
                fillTruthLambdaQA(v0, false);
              }
            } else {
              histos.fill(HIST("xi1820/kplus_antilambda/hInvMassKplusAntiLambda_Mix"), pRes.M());
              if (!(additionalConfig.cfgDelCheck)) {
                histos.fill(HIST("xi1820/kplus_antilambda/hMassPtCentKplusAntiLambda_Mix"), pRes.M(), pRes.Pt(), cent);
              } else {
                auto pCandDel = deltaR(pLambda, pKaon);
                histos.fill(HIST("xi1820/kplus_antilambda/hMassPtCentDelKplusAntiLambda_Mix"), pRes.M(), pRes.Pt(), cent, pCandDel);
              }
            }
          } // End of signal

          // K- + Lambda -> Signal channel for Xi(1820)-
          if (kaonCharge < 0 && isLambda) {
            pLambda = ROOT::Math::PxPyPzEVector(ROOT::Math::PtEtaPhiMVector(v0.pt(), v0.eta(), v0.phi(), v0.mLambda()));
            pRes = pKaon + pLambda;
            auto pCandRapidity = pRes.Rapidity();
            if (!passesRapidity(pCandRapidity)) { // skip candidate if reconstructed rapidity is outside of cut
              continue;
            }
            if constexpr (!IsMix) {
              histos.fill(HIST("xi1820/kminus_lambda/hInvMassKminusLambda"), pRes.M());
              if (!(additionalConfig.cfgDelCheck)) {
                histos.fill(HIST("xi1820/kminus_lambda/hMassPtCentKminusLambda"), pRes.M(), pRes.Pt(), cent);
              } else {
                auto pCandDel = deltaR(pLambda, pKaon);
                histos.fill(HIST("xi1820/kminus_lambda/hMassPtCentDelKminusLambda"), pRes.M(), pRes.Pt(), cent, pCandDel);
              }

              if (additionalConfig.cfgFillRotBkg) {
                for (int i = 0; i < additionalConfig.cfgNrotBkg; i++) {
                  const auto lRotAngle = additionalConfig.cfgNrotBkg == 1
                                           ? 0.5f * (additionalConfig.cfgMinRot + additionalConfig.cfgMaxRot)
                                           : additionalConfig.cfgMinRot + i * ((additionalConfig.cfgMaxRot - additionalConfig.cfgMinRot) / (additionalConfig.cfgNrotBkg - 1));
                  if (cfgFillEventQA) {
                    histos.fill(HIST("Event/hRotBkg"), lRotAngle - o2::constants::math::PI);
                  }
                  if (additionalConfig.cfgRotKaon) {
                    lDaughterRot = pKaon;
                    ROOT::Math::RotationZ rot(lRotAngle);
                    auto p3 = rot * lDaughterRot.Vect();
                    lDaughterRot = LorentzVectorSetXYZM(p3.X(), p3.Y(), p3.Z(), lDaughterRot.M());
                    lResonanceRot = lDaughterRot + pLambda;
                  } else {
                    lDaughterRot = pLambda;
                    ROOT::Math::RotationZ rot(lRotAngle);
                    auto p3 = rot * lDaughterRot.Vect();
                    lDaughterRot = LorentzVectorSetXYZM(p3.X(), p3.Y(), p3.Z(), lDaughterRot.M());
                    lResonanceRot = pKaon + lDaughterRot;
                  }
                  if (!(additionalConfig.cfgDelCheck)) {
                    histos.fill(HIST("xi1820/kminus_lambda/hMassPtCentKminusLambda_Rot"), lResonanceRot.M(), lResonanceRot.Pt(), cent);
                  } else {
                    auto pCandDel = additionalConfig.cfgRotKaon ? deltaR(lDaughterRot, pLambda) : deltaR(lDaughterRot, pKaon);
                    histos.fill(HIST("xi1820/kminus_lambda/hMassPtCentDelKminusLambda_Rot"), lResonanceRot.M(), lResonanceRot.Pt(), cent, pCandDel);
                  }
                }
              }
              if constexpr (IsMC) { // Calculate Acceptance x efficiency for "the particle" channel
                if (v0.motherPDG() != PdgChargedXi1820 || kaon.motherPDG() != v0.motherPDG()) {
                  continue;
                }
                if (kaon.pdgCode() != PDG_t::kKMinus || v0.pdgCode() != PDG_t::kLambda0) {
                  continue;
                }
                if (v0.motherId() < 0 || kaon.motherId() != v0.motherId()) {
                  continue;
                }
                auto pMCPt = v0.motherPt();                                                  // Check particle's pT resolution
                if (additionalConfig.cUseTruthRapidity && !passesRapidity(v0.motherRap())) { // skip candidate if True rapidity of mother particle is outside of cut
                  continue;
                }
                histos.fill(HIST("MC/kminus_lambda/hMCRecoInvMassKminusLambda"), pRes.M());
                histos.fill(HIST("MC/kminus_lambda/hMCRecoMassPtCentKminusLambda"), pRes.M(), pRes.Pt(), cent, pMCPt);

                fillTruthKaonQA(kaon);
                fillTruthLambdaQA(v0, true);
              }
            } else {
              histos.fill(HIST("xi1820/kminus_lambda/hInvMassKminusLambda_Mix"), pRes.M());
              if (!(additionalConfig.cfgDelCheck)) {
                histos.fill(HIST("xi1820/kminus_lambda/hMassPtCentKminusLambda_Mix"), pRes.M(), pRes.Pt(), cent);
              } else {
                auto pCandDel = deltaR(pLambda, pKaon);
                histos.fill(HIST("xi1820/kminus_lambda/hMassPtCentDelKminusLambda_Mix"), pRes.M(), pRes.Pt(), cent, pCandDel);
              }
            }
          } // End of Signal

          // K- + Anti-Lambda -> Bkg channel for charged Xi(1820)
          if (kaonCharge < 0 && !isLambda) {
            pLambda = ROOT::Math::PxPyPzEVector(ROOT::Math::PtEtaPhiMVector(v0.pt(), v0.eta(), v0.phi(), v0.mAntiLambda()));
            pRes = pKaon + pLambda;
            auto pCandRapidity = pRes.Rapidity();
            if (!passesRapidity(pCandRapidity)) { // skip candidate if reconstructed rapidity is outside of cut
              continue;
            }
            if constexpr (!IsMix) {
              histos.fill(HIST("xi1820/kminus_antilambda/hInvMassKminusAntiLambda"), pRes.M());
              if (!(additionalConfig.cfgDelCheck)) {
                histos.fill(HIST("xi1820/kminus_antilambda/hMassPtCentKminusAntiLambda"), pRes.M(), pRes.Pt(), cent);
              } else {
                auto pCandDel = deltaR(pLambda, pKaon);
                histos.fill(HIST("xi1820/kminus_antilambda/hMassPtCentDelKminusAntiLambda"), pRes.M(), pRes.Pt(), cent, pCandDel);
              }
            } else {
              histos.fill(HIST("xi1820/kminus_antilambda/hInvMassKminusAntiLambda_Mix"), pRes.M());
              if (!(additionalConfig.cfgDelCheck)) {
                histos.fill(HIST("xi1820/kminus_antilambda/hMassPtCentKminusAntiLambda_Mix"), pRes.M(), pRes.Pt(), cent);
              } else {
                auto pCandDel = deltaR(pLambda, pKaon);
                histos.fill(HIST("xi1820/kminus_antilambda/hMassPtCentDelKminusAntiLambda_Mix"), pRes.M(), pRes.Pt(), cent, pCandDel);
              }
            }
          } // End of Bkg
        } // End of mass-hypothesis loop
      } // End of V0 loop
    }

    // Fill event QA for after-cuts counters (only for same-event)
    if (!IsMix && cfgFillEventQA) {
      histos.fill(HIST("Event/nLambdasAfterCuts"), selectedLambdas.size());
      histos.fill(HIST("Event/nKaonsAfterCuts"), selectedKaons.size());
    }
  }

  template <bool IsMix, bool IsMC, typename CollisionT, typename V0sT, typename LambdaCollisionT>
  void fillK0sLambda(const CollisionT& collision, const V0sT& k0sCands, const V0sT& lambdaCands, const LambdaCollisionT& lambdaCollision) // Xi(1820) analysis: K0s + Lambda channel, No need to MicroTrack!
  {
    auto cent = collision.cent();
    // The charged same-event callback owns common event/Lambda QA when enabled.
    const bool fillSharedQA = !hasChargedSameEventProcess();

    // Fill event QA histograms
    if (!IsMix && cfgFillEventQA) {
      if (fillSharedQA) {
        histos.fill(HIST("Event/posZ"), collision.posZ());
        histos.fill(HIST("Event/centrality"), cent);
        histos.fill(HIST("Event/posZvsCent"), collision.posZ(), cent);
        histos.fill(HIST("Event/nV0s"), lambdaCands.size());
      }
      histos.fill(HIST("Event/nK0s"), k0sCands.size());
    }

    std::vector<typename V0sT::iterator> selectedK0s;
    selectedK0s.reserve(k0sCands.size());

    // Loop over V0s for K0s
    for (const auto& k0s : k0sCands) {
      // K0s QA before cuts
      if (!IsMix && cfgFillQA) {
        histos.fill(HIST("QAbefore/k0sDCAtoPV"), k0s.pt(), k0s.dcav0topv());
        histos.fill(HIST("QAbefore/k0sMass"), k0s.mK0Short());
        histos.fill(HIST("QAbefore/k0sPt"), k0s.pt());
        histos.fill(HIST("QAbefore/k0sEta"), k0s.eta());
        histos.fill(HIST("QAbefore/k0sCosPA"), k0s.pt(), k0s.v0CosPA());
        histos.fill(HIST("QAbefore/k0sRadius"), k0s.pt(), k0s.transRadius());
        histos.fill(HIST("QAbefore/k0sDauDCA"), k0s.pt(), k0s.daughDCA());
        histos.fill(HIST("QAbefore/k0sDauTPCNsigmaPosPi"), k0s.pt(), k0s.daughterTPCNSigmaPosPi());
        histos.fill(HIST("QAbefore/k0sDauTPCNsigmaNegPi"), k0s.pt(), k0s.daughterTPCNSigmaNegPi());
        histos.fill(HIST("QAbefore/k0sNCrossedRowsPos"), k0s.pt(), k0s.nCrossedRowsPos());
        histos.fill(HIST("QAbefore/k0sNCrossedRowsNeg"), k0s.pt(), k0s.nCrossedRowsNeg());

        float dx = k0s.decayVtxX() - collision.posX();
        float dy = k0s.decayVtxY() - collision.posY();
        float dz = k0s.decayVtxZ() - collision.posZ();
        float l = std::sqrt(dx * dx + dy * dy + dz * dz);
        float p = std::sqrt(k0s.px() * k0s.px() + k0s.py() * k0s.py() + k0s.pz() * k0s.pz());
        auto k0sProperLifetime = (l / (p + 1e-10)) * MassK0Short;
        histos.fill(HIST("QAbefore/k0sProperLifetime"), k0s.pt(), k0sProperLifetime);
        histos.fill(HIST("QAbefore/k0sArmenterosPodolanski"), k0s.alpha(), k0s.qtarm(), k0s.pt());
      }

      if (!k0sCut(collision, k0s)) {
        continue;
      }

      if (!IsMix && cfgFillQA) {
        // K0s QA after cuts
        histos.fill(HIST("QAafter/k0sDCAtoPV"), k0s.pt(), k0s.dcav0topv());
        histos.fill(HIST("QAafter/k0sMass"), k0s.mK0Short());
        histos.fill(HIST("QAafter/k0sPt"), k0s.pt());
        histos.fill(HIST("QAafter/k0sEta"), k0s.eta());
        histos.fill(HIST("QAafter/k0sCosPA"), k0s.pt(), k0s.v0CosPA());
        histos.fill(HIST("QAafter/k0sRadius"), k0s.pt(), k0s.transRadius());
        histos.fill(HIST("QAafter/k0sDauDCA"), k0s.pt(), k0s.daughDCA());
        histos.fill(HIST("QAafter/k0sDauTPCNsigmaPosPi"), k0s.pt(), k0s.daughterTPCNSigmaPosPi());
        histos.fill(HIST("QAafter/k0sDauTPCNsigmaNegPi"), k0s.pt(), k0s.daughterTPCNSigmaNegPi());
        histos.fill(HIST("QAafter/k0sNCrossedRowsPos"), k0s.pt(), k0s.nCrossedRowsPos());
        histos.fill(HIST("QAafter/k0sNCrossedRowsNeg"), k0s.pt(), k0s.nCrossedRowsNeg());

        float dx = k0s.decayVtxX() - collision.posX();
        float dy = k0s.decayVtxY() - collision.posY();
        float dz = k0s.decayVtxZ() - collision.posZ();
        float l = std::sqrt(dx * dx + dy * dy + dz * dz);
        float p = std::sqrt(k0s.px() * k0s.px() + k0s.py() * k0s.py() + k0s.pz() * k0s.pz());
        auto k0sProperLifetime = (l / (p + 1e-10)) * MassK0Short;
        histos.fill(HIST("QAafter/k0sProperLifetime"), k0s.pt(), k0sProperLifetime);
        histos.fill(HIST("QAafter/k0sArmenterosPodolanski"), k0s.alpha(), k0s.qtarm(), k0s.pt());
      }

      selectedK0s.push_back(k0s);
    }

    const auto selectedLambdas = selectLambdas<IsMix>(lambdaCollision, lambdaCands, fillSharedQA);
    for (const auto& k0s : selectedK0s) {
      for (const auto& selected : selectedLambdas) {
        const auto& lambda = selected.candidate;
        if constexpr (!IsMix) {
          if (k0s.globalIndex() == lambda.globalIndex() ||
              sharesAnyDaughterId(v0DaughterIds(k0s), v0DaughterIds(lambda))) {
            continue;
          }
        }

        // 4-vectors
        ROOT::Math::PxPyPzEVector pK0s, pLambda, pRes, lDaughterRot, lResonanceRot;
        pK0s = ROOT::Math::PxPyPzEVector(ROOT::Math::PtEtaPhiMVector(k0s.pt(), k0s.eta(), k0s.phi(), MassK0Short));

        // Lambda and anti-Lambda may both pass; evaluate them independently.
        for (const auto& isLambda : std::array{true, false}) {
          if (isLambda ? !selected.isLambda : !selected.isAntiLambda) {
            continue;
          }

          if (isLambda) {
            pLambda = ROOT::Math::PxPyPzEVector(ROOT::Math::PtEtaPhiMVector(lambda.pt(), lambda.eta(), lambda.phi(), lambda.mLambda()));
            pRes = pK0s + pLambda;
            auto pCandRapidity = pRes.Rapidity();
            if (!passesRapidity(pCandRapidity)) { // skip candidate if reconstructed rapidity is outside of cut
              continue;
            }
            if constexpr (!IsMix) {
              histos.fill(HIST("xi1820/k0s_lambda/hInvMassK0sLambda"), pRes.M());
              if (!(additionalConfig.cfgDelCheck)) {
                histos.fill(HIST("xi1820/k0s_lambda/hMassPtCentK0sLambda"), pRes.M(), pRes.Pt(), cent);
              } else {
                auto pCandDel = deltaR(pLambda, pK0s);
                histos.fill(HIST("xi1820/k0s_lambda/hMassPtCentDelK0sLambda"), pRes.M(), pRes.Pt(), cent, pCandDel);
              }

              if (additionalConfig.cfgFillRotBkg) {
                for (int i = 0; i < additionalConfig.cfgNrotBkg; i++) {
                  const auto lRotAngle = additionalConfig.cfgNrotBkg == 1
                                           ? 0.5f * (additionalConfig.cfgMinRot + additionalConfig.cfgMaxRot)
                                           : additionalConfig.cfgMinRot + i * ((additionalConfig.cfgMaxRot - additionalConfig.cfgMinRot) / (additionalConfig.cfgNrotBkg - 1));
                  if (cfgFillEventQA) {
                    histos.fill(HIST("Event/hRotBkg"), lRotAngle - o2::constants::math::PI);
                  }
                  lDaughterRot = pK0s;
                  ROOT::Math::RotationZ rot(lRotAngle);
                  auto p3 = rot * lDaughterRot.Vect();
                  lDaughterRot = LorentzVectorSetXYZM(p3.X(), p3.Y(), p3.Z(), lDaughterRot.M());
                  lResonanceRot = lDaughterRot + pLambda;
                  if (!(additionalConfig.cfgDelCheck)) {
                    histos.fill(HIST("xi1820/k0s_lambda/hMassPtCentK0sLambda_Rot"), lResonanceRot.M(), lResonanceRot.Pt(), cent);
                  } else {
                    auto lRotDel = deltaR(lDaughterRot, pLambda);
                    histos.fill(HIST("xi1820/k0s_lambda/hMassPtCentDelK0sLambda_Rot"), lResonanceRot.M(), lResonanceRot.Pt(), cent, lRotDel);
                  }
                }
              }

              if constexpr (IsMC) { // Calculate Acceptance x efficiency
                if (lambda.motherPDG() != PdgXi1820Zero || k0s.motherPDG() != lambda.motherPDG()) {
                  continue;
                }
                if (std::abs(k0s.pdgCode()) != PDG_t::kK0Short || lambda.pdgCode() != PDG_t::kLambda0) {
                  continue;
                }
                if (lambda.motherId() < 0 || k0s.motherId() != lambda.motherId()) {
                  continue;
                }
                auto pMCPt = lambda.motherPt();                                                  // Check particle's pT resolution
                if (additionalConfig.cUseTruthRapidity && !passesRapidity(lambda.motherRap())) { // skip candidate if True rapidity of mother particle is outside of cut
                  continue;
                }
                histos.fill(HIST("MC/k0s_lambda/hMCRecoInvMassK0sLambda"), pRes.M());
                histos.fill(HIST("MC/k0s_lambda/hMCRecoMassPtCentK0sLambda"), pRes.M(), pRes.Pt(), cent, pMCPt);
                fillTruthLambdaQA(lambda, true);
                if (cfgFillTruthQA) {
                  histos.fill(HIST("QAMCTrue/k0sPt"), k0s.pt());
                }
              }
            } else {
              histos.fill(HIST("xi1820/k0s_lambda/hInvMassK0sLambda_Mix"), pRes.M());
              if (!(additionalConfig.cfgDelCheck)) {
                histos.fill(HIST("xi1820/k0s_lambda/hMassPtCentK0sLambda_Mix"), pRes.M(), pRes.Pt(), cent);
              } else {
                auto pCandDel = deltaR(pLambda, pK0s);
                histos.fill(HIST("xi1820/k0s_lambda/hMassPtCentDelK0sLambda_Mix"), pRes.M(), pRes.Pt(), cent, pCandDel);
              }
            }
          } // End of Lambda

          if (!isLambda) {
            pLambda = ROOT::Math::PxPyPzEVector(ROOT::Math::PtEtaPhiMVector(lambda.pt(), lambda.eta(), lambda.phi(), lambda.mAntiLambda()));
            pRes = pK0s + pLambda;
            auto pCandRapidity = pRes.Rapidity();
            if (!passesRapidity(pCandRapidity)) { // skip candidate if reconstructed rapidity is outside of cut
              continue;
            }
            if constexpr (!IsMix) {
              histos.fill(HIST("xi1820/k0s_antilambda/hInvMassK0sAntiLambda"), pRes.M());
              if (!(additionalConfig.cfgDelCheck)) {
                histos.fill(HIST("xi1820/k0s_antilambda/hMassPtCentK0sAntiLambda"), pRes.M(), pRes.Pt(), cent);
              } else {
                auto pCandDel = deltaR(pLambda, pK0s);
                histos.fill(HIST("xi1820/k0s_antilambda/hMassPtCentDelK0sAntiLambda"), pRes.M(), pRes.Pt(), cent, pCandDel);
              }

              if (additionalConfig.cfgFillRotBkg) {
                for (int i = 0; i < additionalConfig.cfgNrotBkg; i++) {
                  const auto lRotAngle = additionalConfig.cfgNrotBkg == 1
                                           ? 0.5f * (additionalConfig.cfgMinRot + additionalConfig.cfgMaxRot)
                                           : additionalConfig.cfgMinRot + i * ((additionalConfig.cfgMaxRot - additionalConfig.cfgMinRot) / (additionalConfig.cfgNrotBkg - 1));
                  if (cfgFillEventQA) {
                    histos.fill(HIST("Event/hRotBkg"), lRotAngle - o2::constants::math::PI);
                  }
                  lDaughterRot = pK0s;
                  ROOT::Math::RotationZ rot(lRotAngle);
                  auto p3 = rot * lDaughterRot.Vect();
                  lDaughterRot = LorentzVectorSetXYZM(p3.X(), p3.Y(), p3.Z(), lDaughterRot.M());
                  lResonanceRot = lDaughterRot + pLambda;
                  if (!(additionalConfig.cfgDelCheck)) {
                    histos.fill(HIST("xi1820/k0s_antilambda/hMassPtCentK0sAntiLambda_Rot"), lResonanceRot.M(), lResonanceRot.Pt(), cent);
                  } else {
                    auto lRotDel = deltaR(lDaughterRot, pLambda);
                    histos.fill(HIST("xi1820/k0s_antilambda/hMassPtCentDelK0sAntiLambda_Rot"), lResonanceRot.M(), lResonanceRot.Pt(), cent, lRotDel);
                  }
                }
              }

              if constexpr (IsMC) { // Calculate Acceptance x efficiency
                if (lambda.motherPDG() != -PdgXi1820Zero || k0s.motherPDG() != lambda.motherPDG()) {
                  continue;
                }
                if (std::abs(k0s.pdgCode()) != PDG_t::kK0Short || lambda.pdgCode() != PDG_t::kLambda0Bar) {
                  continue;
                }
                if (lambda.motherId() < 0 || k0s.motherId() != lambda.motherId()) {
                  continue;
                }
                auto pMCPt = lambda.motherPt();                                                  // Check particle's pT resolution
                if (additionalConfig.cUseTruthRapidity && !passesRapidity(lambda.motherRap())) { // skip candidate if True rapidity of mother particle is outside of cut
                  continue;
                }
                histos.fill(HIST("MC/k0s_antilambda/hMCRecoInvMassK0sAntiLambda"), pRes.M());
                histos.fill(HIST("MC/k0s_antilambda/hMCRecoMassPtCentK0sAntiLambda"), pRes.M(), pRes.Pt(), cent, pMCPt);
                fillTruthLambdaQA(lambda, false);
                if (cfgFillTruthQA) {
                  histos.fill(HIST("QAMCTrue/k0sPt"), k0s.pt());
                }
              }
            } else {
              histos.fill(HIST("xi1820/k0s_antilambda/hInvMassK0sAntiLambda_Mix"), pRes.M());
              if (!(additionalConfig.cfgDelCheck)) {
                histos.fill(HIST("xi1820/k0s_antilambda/hMassPtCentK0sAntiLambda_Mix"), pRes.M(), pRes.Pt(), cent);
              } else {
                auto pCandDel = deltaR(pLambda, pK0s);
                histos.fill(HIST("xi1820/k0s_antilambda/hMassPtCentDelK0sAntiLambda_Mix"), pRes.M(), pRes.Pt(), cent, pCandDel);
              }
            }
          }
        } // End of mass-hypothesis loop
      } // End of loop over Lambda candidates
    } // End of loop over K0s candidates

    // Fill event QA for after-cuts counters (only for same-event)
    if (!IsMix && cfgFillEventQA) {
      if (fillSharedQA) {
        histos.fill(HIST("Event/nLambdasAfterCuts"), selectedLambdas.size());
      }
      histos.fill(HIST("Event/nK0sAfterCuts"), selectedK0s.size());
    }
  }

  void processDummy(const aod::ResoCollision& /*collision*/)
  {
    // Dummy function to satisfy the compiler
  }
  PROCESS_SWITCH(Xi1820Analysis, processDummy, "Process Dummy", true);

  void processDataWithTracks(ResoCollisions::iterator const& resoCollision,
                             aod::ResoV0s const& resoV0s,
                             ResoTracks const& resoTracks)
  {
    if (additionalConfig.cRecoINELgt0 && !resoCollision.isRecINELgt0()) {
      return; // skip event if RecoINEL>0 selection is enabled and event does not pass it
    }
    fillChargedKLambda<false, false, false>(resoCollision, resoV0s, resoTracks, resoCollision);
  }
  PROCESS_SWITCH(Xi1820Analysis, processDataWithTracks, "Process Event with ResoTracks", false);

  void processDataWithMicroTracks(ResoCollisions::iterator const& resoCollision,
                                  aod::ResoV0s const& resoV0s,
                                  ResoMicroTracks const& resoMicroTracks)
  {
    if (additionalConfig.cRecoINELgt0 && !resoCollision.isRecINELgt0()) {
      return; // skip event if RecoINEL>0 selection is enabled and event does not pass it
    }
    fillChargedKLambda<false, true, false>(resoCollision, resoV0s, resoMicroTracks, resoCollision);
  }
  PROCESS_SWITCH(Xi1820Analysis, processDataWithMicroTracks, "Process Event with ResoMicroTracks", false);

  template <typename FirstCollision, typename SecondCollision>
  bool acceptsMixedCollisions(const FirstCollision& first, const SecondCollision& second) const
  {
    if (additionalConfig.cRecoINELgt0.value && (!first.isRecINELgt0() || !second.isRecINELgt0())) {
      return false;
    }
    return std::isfinite(first.bMagField()) && first.bMagField() == second.bMagField();
  }

  void processMixedEventWithTracks(ResoCollisions const& resoCollisions,
                                   aod::ResoV0s const& resoV0s,
                                   ResoTracks const& resoTracks)
  {
    if (nEvtMixing <= 0) {
      return;
    }
    auto tracksTuple = std::make_tuple(resoTracks, resoV0s);
    BinningTypeVertexContributor binning{{cfgVtxBins, cfgMultBins}, true};
    Pair<ResoCollisions, ResoTracks, aod::ResoV0s, BinningTypeVertexContributor> pairs{binning, nEvtMixing, -1, resoCollisions, tracksTuple, &cache};
    for (auto& [collision1, tracks1, collision2, v0s2] : pairs) { // o2-linter: disable=const-ref-in-for-loop (structured bindings from Pair iterator cannot be const)
      if (!acceptsMixedCollisions(collision1, collision2)) {
        continue;
      }
      fillChargedKLambda<true, false, false>(collision1, v0s2, tracks1, collision2);
    }
  }
  PROCESS_SWITCH(Xi1820Analysis, processMixedEventWithTracks, "Process Mixed Event with ResoTracks", false);

  void processMixedEventWithMicroTracks(ResoCollisions const& resoCollisions,
                                        aod::ResoV0s const& resoV0s,
                                        ResoMicroTracks const& resoMicroTracks)
  {
    if (nEvtMixing <= 0) {
      return;
    }
    auto tracksTuple = std::make_tuple(resoMicroTracks, resoV0s);
    BinningTypeVertexContributor binning{{cfgVtxBins, cfgMultBins}, true};
    Pair<ResoCollisions, ResoMicroTracks, aod::ResoV0s, BinningTypeVertexContributor> pairs{binning, nEvtMixing, -1, resoCollisions, tracksTuple, &cache};
    for (auto& [collision1, tracks1, collision2, v0s2] : pairs) { // o2-linter: disable=const-ref-in-for-loop (structured bindings from Pair iterator cannot be const)
      if (!acceptsMixedCollisions(collision1, collision2)) {
        continue;
      }
      fillChargedKLambda<true, true, false>(collision1, v0s2, tracks1, collision2);
    }
  }
  PROCESS_SWITCH(Xi1820Analysis, processMixedEventWithMicroTracks, "Process Mixed Event with ResoMicroTracks 001", false);

  template <typename DaughterTable>
  static std::shared_ptr<arrow::Table> copyMixingRows(DaughterTable const& daughters)
  {
    // Copy only this collision's rows, without retaining the source DF buffers.
    const auto source = daughters.asArrowTableConstrained();
    if (source->num_rows() == 0) {
      return arrow::Table::MakeEmpty(source->schema()).ValueOrDie();
    }
    std::vector<std::shared_ptr<arrow::ChunkedArray>> columns;
    columns.reserve(source->num_columns());
    for (const auto& column : source->columns()) {
      columns.push_back(std::make_shared<arrow::ChunkedArray>(arrow::Concatenate(column->chunks()).ValueOrDie()));
    }
    return arrow::Table::Make(source->schema(), columns);
  }

  // FIFO history is task-local and retains only the copied MicroTrack rows.
  // Use a single run/condition set; a field change clears all stored event pools.
  void processMEDF(ResoCollisions::iterator const& collision,
                   ResoMicroTracks const& tracks,
                   aod::ResoV0s const& v0s)
  {
    if (nEvtMixing <= 0 || !std::isfinite(collision.bMagField())) {
      return;
    }
    if (mixingBField != collision.bMagField()) {
      mixingPools.clear();
      mixingBField = collision.bMagField();
    }
    BinningTypeVertexContributor binning{{cfgVtxBins, cfgMultBins}, true};
    const int bin = binning.getBin(std::make_tuple(collision.posZ(), collision.cent()));
    if (bin < 0) {
      return;
    }
    const MixingCollision current{.x = collision.posX(), .y = collision.posY(), .z = collision.posZ(), .centrality = collision.cent(), .recINELgt0 = collision.isRecINELgt0()};
    auto& pool = mixingPools[bin];
    for (const auto& previous : pool) {
      if (additionalConfig.cRecoINELgt0 && (!previous.collision.recINELgt0 || !current.recINELgt0)) {
        continue;
      }
      fillChargedKLambda<true, true, false>(previous.collision, v0s, ResoMicroTracks{previous.daughters}, current);
    }
    // Insert after pairing; DF-local IDs may repeat across frames.
    pool.push_back({current, copyMixingRows(tracks)});
    while (pool.size() > static_cast<std::size_t>(nEvtMixing.value)) {
      pool.pop_front();
    }
  }
  PROCESS_SWITCH(Xi1820Analysis, processMEDF, "Process cross-DF mixing with MicroTracks 001", false);

  // K0s + Lambda analysis
  void processK0sLambda(ResoCollisions::iterator const& resoCollision,
                        aod::ResoV0s const& resoV0s)
  {
    if (additionalConfig.cRecoINELgt0 && !resoCollision.isRecINELgt0()) {
      return; // skip event if RecoINEL>0 selection is enabled and event does not pass it
    }
    fillK0sLambda<false, false>(resoCollision, resoV0s, resoV0s, resoCollision);
  }
  PROCESS_SWITCH(Xi1820Analysis, processK0sLambda, "Process K0s + Lambda", false);

  // K0s + Lambda mixed event analysis
  void processK0sLambdaMixedEvent(ResoCollisions const& resoCollisions,
                                  aod::ResoV0s const& resoV0s)
  {

    if (nEvtMixing <= 0) {
      return;
    }
    auto v0sV0sTuple = std::make_tuple(resoV0s, resoV0s);
    BinningTypeVertexContributor colBinning{{cfgVtxBins, cfgMultBins}, true};
    Pair<ResoCollisions, aod::ResoV0s, aod::ResoV0s, BinningTypeVertexContributor> pairs{colBinning, nEvtMixing, -1, resoCollisions, v0sV0sTuple, &cache};

    for (auto& [collision1, k0s1, collision2, lambda2] : pairs) { // o2-linter: disable=const-ref-in-for-loop (structured bindings from Pair iterator cannot be const)
      if (!acceptsMixedCollisions(collision1, collision2)) {
        continue;
      }
      fillK0sLambda<true, false>(collision1, k0s1, lambda2, collision2);
    }
  }
  PROCESS_SWITCH(Xi1820Analysis, processK0sLambdaMixedEvent, "Process K0s + Lambda Mixed Event", false);

  // Cross-DF neutral mixing: old-event K0s with current-event Lambda/anti-Lambda.
  // Keep a separate bounded V0 pool; use inputs from one run/condition set.
  void processK0sLambdaMEDF(ResoCollisions::iterator const& collision,
                            aod::ResoV0s const& v0s)
  {
    if (nEvtMixing <= 0 || !std::isfinite(collision.bMagField())) {
      return;
    }
    if (neutralMixingBField != collision.bMagField()) {
      neutralMixingPools.clear();
      neutralMixingBField = collision.bMagField();
    }
    BinningTypeVertexContributor binning{{cfgVtxBins, cfgMultBins}, true};
    const int bin = binning.getBin(std::make_tuple(collision.posZ(), collision.cent()));
    if (bin < 0) {
      return;
    }
    const MixingCollision current{.x = collision.posX(), .y = collision.posY(), .z = collision.posZ(), .centrality = collision.cent(), .recINELgt0 = collision.isRecINELgt0()};
    auto& pool = neutralMixingPools[bin];
    const bool fillCurrentQA = cfgFillEventQA && (!additionalConfig.cRecoINELgt0 || current.recINELgt0);
    if (fillCurrentQA) {
      histos.fill(HIST("MixingQA/K0sLambdaMEDF/hPosZvsCent"), current.posZ(), current.cent());
      histos.fill(HIST("MixingQA/K0sLambdaMEDF/hPoolSize"), pool.size());
    }
    int nPartners = 0;
    for (const auto& previous : pool) {
      if (additionalConfig.cRecoINELgt0 && (!previous.collision.recINELgt0 || !current.recINELgt0)) {
        continue;
      }
      if (cfgFillEventQA) {
        ++nPartners;
        histos.fill(HIST("MixingQA/K0sLambdaMEDF/hPosZPair"), previous.collision.posZ(), current.posZ());
        histos.fill(HIST("MixingQA/K0sLambdaMEDF/hCentPair"), previous.collision.cent(), current.cent());
        histos.fill(HIST("MixingQA/K0sLambdaMEDF/hV0CountsPair"), previous.daughters->num_rows(), v0s.size());
      }
      // Each V0 uses its own event's PV; original row IDs are not compared across events.
      fillK0sLambda<true, false>(previous.collision, aod::ResoV0s{previous.daughters}, v0s, current);
    }
    if (fillCurrentQA) {
      histos.fill(HIST("MixingQA/K0sLambdaMEDF/hNPartners"), nPartners);
    }
    // Insert after pairing to avoid self-mixing, even when DF-local IDs repeat.
    pool.push_back({current, copyMixingRows(v0s)});
    while (pool.size() > static_cast<std::size_t>(nEvtMixing.value)) {
      pool.pop_front();
    }
  }
  PROCESS_SWITCH(Xi1820Analysis, processK0sLambdaMEDF, "Process cross-DF K0s + Lambda mixing", false);

  // MC processes for charged K + Lambda analysis
  void processMCWithTracks(ResoMCCols::iterator const& resoMCcollision,
                           soa::Join<aod::ResoV0s, aod::ResoMCV0s> const& resoMCV0s,
                           soa::Join<aod::ResoTracks, aod::ResoTrackTracks, aod::ResoMCTracks> const& resoMCTracks)
  {
    if ((additionalConfig.cRecoINELgt0 && !resoMCcollision.isRecINELgt0()) ||
        (additionalConfig.cMCINELgt0 && !resoMCcollision.isINELgt0()) ||
        (additionalConfig.cMCVtxIn10 && !resoMCcollision.isVtxIn10())) {
      return;
    }
    fillChargedKLambda<false, false, true>(resoMCcollision, resoMCV0s, resoMCTracks, resoMCcollision);
  }
  PROCESS_SWITCH(Xi1820Analysis, processMCWithTracks, "Process MC for charged K + Lambda", false);

  void processMCK0sLambda(ResoMCCols::iterator const& resoMCCollision,
                          soa::Join<aod::ResoV0s, aod::ResoMCV0s> const& resoMCV0s)
  {
    if ((additionalConfig.cRecoINELgt0 && !resoMCCollision.isRecINELgt0()) ||
        (additionalConfig.cMCINELgt0 && !resoMCCollision.isINELgt0()) ||
        (additionalConfig.cMCVtxIn10 && !resoMCCollision.isVtxIn10())) {
      return;
    }
    fillK0sLambda<false, true>(resoMCCollision, resoMCV0s, resoMCV0s, resoMCCollision);
  }
  PROCESS_SWITCH(Xi1820Analysis, processMCK0sLambda, "Process MC K0s + Lambda", false);

  void processMCWithMicroTracks(ResoMCCols::iterator const& resoCollision,
                                soa::Join<aod::ResoV0s, aod::ResoMCV0s> const& resoV0s,
                                soa::Join<ResoMicroTracks, aod::ResoMCMicroTracks_001> const& resoMicroTracks)
  {
    if ((additionalConfig.cRecoINELgt0 && !resoCollision.isRecINELgt0()) ||
        (additionalConfig.cMCINELgt0 && !resoCollision.isINELgt0()) ||
        (additionalConfig.cMCVtxIn10 && !resoCollision.isVtxIn10())) {
      return;
    }
    // Module 001 contains only selected reconstructed collisions.
    fillChargedKLambda<false, true, true>(resoCollision, resoV0s, resoMicroTracks, resoCollision);
  }
  PROCESS_SWITCH(Xi1820Analysis, processMCWithMicroTracks, "Process MC with ResoMicroTracks 001", false);

  void processMCGen(ResoMCCols::iterator const& resoCollision, // Calculate denominator for the acceptance x efficiency and a part of Event-factor (for selected evennts)
                    aod::ResoMCParents_001 const& resoParents)
  {
    auto multiplicity = resoCollision.mcMultiplicity();
    auto inCent = resoCollision.cent();
    if ((additionalConfig.cRecoINELgt0 && !resoCollision.isRecINELgt0()) ||
        (additionalConfig.cMCINELgt0 && !resoCollision.isINELgt0()) ||
        (additionalConfig.cMCVtxIn10 && !resoCollision.isVtxIn10())) {
      return;
    }
    histos.fill(HIST("multQA/h2MultCentMC"), inCent, multiplicity);
    for (const auto& part : resoParents) { // loop over all pre-filtered Gen particle on selected events
      auto pdgMother = part.pdgCode();
      if (std::abs(pdgMother) != PdgChargedXi1820 && std::abs(pdgMother) != PdgXi1820Zero) {
        continue;
      }
      if (!passesRapidity(part.y())) {
        continue; // skip if rapidity of the particle is outside of cut
      }
      auto motherPt = part.pt();
      auto daughter1PDG = part.daughterPDG1();
      auto daughter2PDG = part.daughterPDG2();

      if (std::abs(pdgMother) == PdgChargedXi1820) { // Explicity check for the safety.
        // K- + Anti-Lambda,  K+ + Anti-Lambda
        if ((daughter1PDG == PDG_t::kKMinus && daughter2PDG == PDG_t::kLambda0) ||
            (daughter1PDG == PDG_t::kLambda0 && daughter2PDG == PDG_t::kKMinus)) {
          histos.fill(HIST("MC/hMCGenPtCentMultKminusLambda"), motherPt, inCent, multiplicity);
        } else if ((daughter1PDG == PDG_t::kKPlus && daughter2PDG == PDG_t::kLambda0Bar) ||
                   (daughter1PDG == PDG_t::kLambda0Bar && daughter2PDG == PDG_t::kKPlus)) {
          histos.fill(HIST("MC/hMCGenPtCentMultKplusAntiLambda"), motherPt, inCent, multiplicity);
        }
      } else {
        // K0s + Lambda, K0s + Anti-Lambda
        if ((std::abs(daughter1PDG) == PDG_t::kK0Short && daughter2PDG == PDG_t::kLambda0) ||
            (daughter1PDG == PDG_t::kLambda0 && std::abs(daughter2PDG) == PDG_t::kK0Short)) {
          histos.fill(HIST("MC/hMCGenPtCentMultK0sLambda"), motherPt, inCent, multiplicity);
        } else if ((std::abs(daughter1PDG) == PDG_t::kK0Short && daughter2PDG == PDG_t::kLambda0Bar) ||
                   (daughter1PDG == PDG_t::kLambda0Bar && std::abs(daughter2PDG) == PDG_t::kK0Short)) {
          histos.fill(HIST("MC/hMCGenPtCentMultK0sAntiLambda"), motherPt, inCent, multiplicity);
        }
      }
    }
  }
  PROCESS_SWITCH(Xi1820Analysis, processMCGen, "Process Event for MC (Generated at selected events)", false);

  void processMCTruth(aod::McParticles const& mcParticles) // ->Let's keep it and use for injected MC QA...!
  {
    // Process MC generated particles (no reconstruction requirement)
    // Charged and neutral Xi(1820) states have distinct PDG codes.

    for (const auto& mcParticle : mcParticles) {
      // Look for Xi(1820) - PDG code can vary, check for resonance mass ~1820 MeV
      int pdg = mcParticle.pdgCode();

      // Select only the configured charged and neutral Xi(1820) states.
      if (std::abs(pdg) != PdgChargedXi1820 && std::abs(pdg) != PdgXi1820Zero) {
        continue;
      }

      // Fill generated level histograms
      auto pt = mcParticle.pt();
      auto eta = mcParticle.eta();
      auto y = mcParticle.y();

      histos.fill(HIST("MC/hMCTruthXi1820Pt"), pt);
      histos.fill(HIST("MC/hMCTruthXi1820PtEta"), pt, eta);
      histos.fill(HIST("MC/hMCTruthXi1820Y"), y);

      // Get daughters
      auto daughters = mcParticle.daughters_as<aod::McParticles>();
      if (daughters.size() != ExpectedDaughters) {
        continue;
      }

      int daughter1PDG = 0, daughter2PDG = 0;
      ROOT::Math::PxPyPzEVector p1, p2, pMother;

      int iDaughter = 0;
      for (const auto& daughter : daughters) {
        if (iDaughter == 0) {
          daughter1PDG = daughter.pdgCode();
          p1.SetPxPyPzE(daughter.px(), daughter.py(), daughter.pz(), daughter.e());
        } else {
          daughter2PDG = daughter.pdgCode();
          p2.SetPxPyPzE(daughter.px(), daughter.py(), daughter.pz(), daughter.e());
        }
        iDaughter++;
      }
      pMother = p1 + p2;

      // Check decay channels
      auto motherPt = pMother.Pt();
      auto motherM = pMother.M();

      // K- + Lambda
      if ((daughter1PDG == PDG_t::kKMinus && daughter2PDG == PDG_t::kLambda0) ||
          (daughter1PDG == PDG_t::kLambda0 && daughter2PDG == PDG_t::kKMinus)) {
        histos.fill(HIST("MC/hMCTruthInvMassKminusLambda"), motherM);
        histos.fill(HIST("MC/hMCTruthMassPtKminusLambda"), motherM, motherPt);
      }

      // K+ + Anti-Lambda
      if ((daughter1PDG == PDG_t::kKPlus && daughter2PDG == PDG_t::kLambda0Bar) ||
          (daughter1PDG == PDG_t::kLambda0Bar && daughter2PDG == PDG_t::kKPlus)) {
        histos.fill(HIST("MC/hMCTruthInvMassKplusAntiLambda"), motherM);
        histos.fill(HIST("MC/hMCTruthMassPtKplusAntiLambda"), motherM, motherPt);
      }

      // K0s + Lambda
      if ((std::abs(daughter1PDG) == PDG_t::kK0Short && daughter2PDG == PDG_t::kLambda0) ||
          (daughter1PDG == PDG_t::kLambda0 && std::abs(daughter2PDG) == PDG_t::kK0Short)) {
        histos.fill(HIST("MC/hMCTruthInvMassK0sLambda"), motherM);
        histos.fill(HIST("MC/hMCTruthMassPtK0sLambda"), motherM, motherPt);
      }

      // K0s + Anti-Lambda
      if ((std::abs(daughter1PDG) == PDG_t::kK0Short && daughter2PDG == PDG_t::kLambda0Bar) ||
          (daughter1PDG == PDG_t::kLambda0Bar && std::abs(daughter2PDG) == PDG_t::kK0Short)) {
        histos.fill(HIST("MC/hMCTruthInvMassK0sAntiLambda"), motherM);
        histos.fill(HIST("MC/hMCTruthMassPtK0sAntiLambda"), motherM, motherPt);
      }
    }
  }
  PROCESS_SWITCH(Xi1820Analysis, processMCTruth, "Process MC Truth particles", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& context)
{
  return WorkflowSpec{adaptAnalysisTask<Xi1820Analysis>(context)};
}
