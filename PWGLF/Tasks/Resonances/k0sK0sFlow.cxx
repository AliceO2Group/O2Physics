// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file k0sK0sFlow.cxx
/// \brief K0S and K0S-K0S spectra, angular distributions, SP and EP observables
/// \author Satoshi Yano

#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/Qvectors.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/Configurable.h"
#include "Framework/GroupedCombinations.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include <DCAFitter/DCAFitterN.h>

#include <Math/GenVector/Boost.h>
#include <Math/Vector3D.h>
#include <Math/Vector4D.h>
#include <TVector2.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

struct K0sK0sFlow {
  enum class PionPidMode : int {
    TpcWithTofIfAvailable = 0,
    TpcOnly = 1,
    TpcAndTof = 2,
    TofRequiredAboveMomentum = 4
  };

  enum class EventSelectionStep : int {
    All = 1,
    Vertex,
    Sel8,
    NoSameBunchPileup,
    NoTimeFrameBorder,
    NoItsRofBorder,
    Centrality,
    NSteps
  };

  enum class TrackSelectionStep : int {
    All = 1,
    Kinematics,
    TpcClusters,
    TrackQuality,
    TpcPid,
    TofPid,
    NSteps
  };

  enum class K0sSelectionStep : int {
    All = 1,
    Daughters,
    Kinematics,
    Topology,
    DaughterDca,
    DecayGeometry,
    Armenteros,
    LambdaRejection,
    NSteps
  };

  enum class PairSelectionStep : int {
    All = 1,
    NoSharedDaughters,
    SelectedK0s,
    ClosePair,
    Kinematics,
    SignalMass,
    NSteps
  };

  enum class LikeSignCategory : int {
    PositivePositive,
    NegativeNegative,
    PositiveNegative
  };

  enum class SidebandRegion : int {
    Lower = -1,
    Outside = 0,
    Upper = 1
  };

  static constexpr double VectorTolerance = 1.e-12;
  static constexpr int QVectorComponents = 4;

  template <typename E>
  static constexpr auto enumValue(E value)
  {
    return static_cast<std::underlying_type_t<E>>(value);
  }

  HistogramRegistry registry{"registry"};
  SliceCache cache;
  o2::vertexing::DCAFitterN<2> lsFitter;

  struct : ConfigurableGroup {
    Configurable<float> maxAbsZvtx{"maxAbsZvtx", 10.f, "maximum |z vertex| (cm)"};
    Configurable<bool> requireSel8{"requireSel8", true, "require sel8"};
    Configurable<bool> requireNoSameBunchPileup{"requireNoSameBunchPileup", true, "reject same-bunch pileup"};
    Configurable<bool> requireNoTimeFrameBorder{"requireNoTimeFrameBorder", true, "reject TF borders"};
    Configurable<bool> requireNoITSROFrameBorder{"requireNoITSROFrameBorder", true, "reject ITS ROF borders"};
    Configurable<float> minCentrality{"minCentrality", 0.f, "minimum FT0C centrality"};
    Configurable<float> maxCentrality{"maxCentrality", 100.f, "maximum FT0C centrality"};
  } eventCuts;

  struct : ConfigurableGroup {
    Configurable<float> daughterMinPt{"daughterMinPt", 0.1f, "daughter minimum pT"};
    Configurable<float> maxAbsEta{"maxAbsEta", 0.8f, "daughter maximum |eta|"};
    Configurable<int> minTPCRows{"minTPCRows", 70, "minimum TPC crossed rows"};
    Configurable<float> minTPCRowsOverFindable{"minTPCRowsOverFindable", 0.8f, "minimum crossed/findable ratio"};
    Configurable<float> maxTPCChi2NCl{"maxTPCChi2NCl", 4.f, "maximum TPC chi2/cluster"};
    Configurable<float> maxITSChi2NCl{"maxITSChi2NCl", 36.f, "maximum ITS chi2/cluster"};
    Configurable<int> minITSClusters{"minITSClusters", 0, "minimum ITS clusters; zero disables"};
    Configurable<int> pidMode{"pidMode", static_cast<int>(PionPidMode::TpcWithTofIfAvailable), "0 TPC and TOF if present; 1 TPC only; 2 require TOF; 4 require TOF above momentum threshold"};
    Configurable<float> maxTPCNSigmaPi{"maxTPCNSigmaPi", 5.f, "maximum pion TPC |n sigma|"};
    Configurable<float> maxTOFNSigmaPi{"maxTOFNSigmaPi", 3.f, "maximum pion TOF |n sigma|"};
    Configurable<float> tofRequiredAboveP{"tofRequiredAboveP", 0.6f, "momentum threshold for PID mode 4"};
  } daughterCuts;

  struct : ConfigurableGroup {
    Configurable<float> k0sMinPt{"k0sMinPt", 0.f, "K0S minimum pT"};
    Configurable<float> k0sMaxPt{"k0sMaxPt", 20.f, "K0S maximum pT"};
    Configurable<float> k0sMaxAbsY{"k0sMaxAbsY", 0.5f, "K0S maximum |y|"};
    Configurable<float> minCosPA{"minCosPA", 0.98f, "minimum cos(PA)"};
    Configurable<float> maxDCADaughters{"maxDCADaughters", 0.3f, "maximum daughter DCA (cm)"};
    Configurable<float> minDCAPosToPV{"minDCAPosToPV", 0.05f, "minimum positive daughter |DCA to PV| (cm)"};
    Configurable<float> minDCANegToPV{"minDCANegToPV", 0.05f, "minimum negative daughter |DCA to PV| (cm)"};
    Configurable<float> minRadius{"minRadius", 0.5f, "minimum decay radius (cm)"};
    Configurable<float> maxRadius{"maxRadius", 200.f, "maximum decay radius (cm)"};
    Configurable<float> maxCtau{"maxCtau", 20.f, "maximum c tau (cm)"};
    Configurable<bool> useArmenterosCut{"useArmenterosCut", true, "apply K0S Armenteros cut"};
    Configurable<float> armenterosSlope{"armenterosSlope", 0.2f, "require qT > slope |alpha|"};
    Configurable<bool> rejectLambda{"rejectLambda", true, "reject Lambda and anti-Lambda"};
    Configurable<float> lambdaRejectionWindow{"lambdaRejectionWindow", 0.01f, "Lambda rejection half-width"};
    Configurable<float> k0sMassMin{"k0sMassMin", 0.47f, "K0S pair-building mass minimum"};
    Configurable<float> k0sMassMax{"k0sMassMax", 0.525f, "K0S pair-building mass maximum"};
    Configurable<float> k0sLowerSidebandMin{"k0sLowerSidebandMin", 0.42f, "K0S lower sideband minimum"};
    Configurable<float> k0sLowerSidebandMax{"k0sLowerSidebandMax", 0.46f, "K0S lower sideband maximum"};
    Configurable<float> k0sUpperSidebandMin{"k0sUpperSidebandMin", 0.535f, "K0S upper sideband minimum"};
    Configurable<float> k0sUpperSidebandMax{"k0sUpperSidebandMax", 0.575f, "K0S upper sideband maximum"};
  } k0Cuts;

  struct : ConfigurableGroup {
    Configurable<float> pairMinPt{"pairMinPt", 0.f, "pair minimum pT"};
    Configurable<float> pairMaxPt{"pairMaxPt", 20.f, "pair maximum pT"};
    Configurable<float> pairMaxAbsY{"pairMaxAbsY", 0.5f, "pair maximum |y|"};
    Configurable<float> pairMassMin{"pairMassMin", 0.99f, "pair mass minimum"};
    Configurable<float> pairMassMax{"pairMassMax", 3.f, "pair mass maximum"};
    Configurable<bool> enableClosePairRejection{"enableClosePairRejection", false, "enable minimum DeltaR cut"};
    Configurable<float> minDeltaR{"minDeltaR", 0.0f, "minimum K0S-K0S DeltaR"};
    Configurable<int> mixingDepth{"mixingDepth", 5, "number of mixed events"};
  } pairCuts;

  struct : ConfigurableGroup {
    Configurable<bool> enableLikeSign{"enableLikeSign", true, "build like-sign V0 and V0-V0 background candidates"};
    Configurable<float> likeSignBz{"likeSignBz", 5.f, "magnetic field Bz used by the like-sign DCA fitter (kG)"};
    Configurable<float> likeSignFitterMaxR{"likeSignFitterMaxR", 200.f, "maximum radius used by the like-sign DCA fitter (cm)"};
    Configurable<float> likeSignFitterMaxChi2{"likeSignFitterMaxChi2", 1.e9f, "maximum fitter chi2"};
  } likeSign;

  struct : ConfigurableGroup {
    Configurable<int> harmonic{"harmonic", 2, "flow harmonic"};
    Configurable<int> correctionStep{"correctionStep", 3, "Q-vector correction step (0-3)"};
    Configurable<int> numberOfDetectors{"numberOfDetectors", 7, "detectors stored in Qvectors"};
    Configurable<int> ft0cId{"ft0cId", 0, "FT0C detector index"};
    Configurable<int> ft0aId{"ft0aId", 1, "FT0A detector index"};
    Configurable<int> tpcPosId{"tpcPosId", 4, "TPC positive-eta detector index"};
    Configurable<float> minQAmplitude{"minQAmplitude", 1.e-4f, "minimum Q-vector amplitude"};
  } flow;

  ConfigurableAxis axisCent{"axisCent", {VARIABLE_WIDTH, 0., 5., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100.}, "centrality"};
  ConfigurableAxis axisK0Pt{"axisK0Pt", {VARIABLE_WIDTH, 0., 0.5, 1., 1.5, 2., 2.5, 3., 4., 6., 8., 10., 13., 20.}, "K0S pT"};
  ConfigurableAxis axisPairPt{"axisPairPt", {VARIABLE_WIDTH, 0., 0.5, 1., 1.5, 2., 2.5, 3., 4., 6., 8., 10., 13., 20.}, "pair pT"};
  ConfigurableAxis axisK0Mass{"axisK0Mass", {200, 0.4, 0.6}, "pi-pi mass"};
  ConfigurableAxis axisPairMass{"axisPairMass", {402, 0.99, 3.0}, "K0S-K0S mass"};
  ConfigurableAxis axisAbsCosTheta{"axisAbsCosTheta", {10, 0., 1.}, "|cos(theta*)|"};
  ConfigurableAxis axisPhiStar{"axisPhiStar", {12, -3.141592653589793, 3.141592653589793}, "phi*"};
  ConfigurableAxis axisMixZvtx{"axisMixZvtx", {VARIABLE_WIDTH, -10., -5., 0., 5., 10.}, "event-mixing vertex-z bins"};
  ConfigurableAxis axisMixCent{"axisMixCent", {VARIABLE_WIDTH, 0., 5., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100.}, "event-mixing centrality bins"};

  using Collisions = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::Qvectors>;
  using DaughterTracks = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTOFbeta>;
  Preslice<aod::V0Datas> perCollisionV0 = aod::v0data::collisionId;
  Preslice<DaughterTracks> perCollisionTrack = aod::track::collisionId;
  using MixingBinning = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C>;

  void init(InitContext const&)
  {
    lsFitter.setPropagateToPCA(true);
    lsFitter.setMaxR(likeSign.likeSignFitterMaxR);
    lsFitter.setMinParamChange(1.e-3f);
    lsFitter.setMinRelChi2Change(0.9f);
    lsFitter.setMaxDZIni(1.e9f);
    lsFitter.setMaxChi2(likeSign.likeSignFitterMaxChi2);
    lsFitter.setUseAbsDCA(true);
    lsFitter.setBz(likeSign.likeSignBz);
    AxisSpec cent{axisCent, "FT0C centrality (%)"}, k0pt{axisK0Pt, "p_{T}(K^{0}_{S})"}, pairpt{axisPairPt, "p_{T}(K^{0}_{S}K^{0}_{S})"};
    AxisSpec k0mass{axisK0Mass, "m_{#pi#pi}"}, pairmass{axisPairMass, "m_{K^{0}_{S}K^{0}_{S}}"};
    AxisSpec absct{axisAbsCosTheta, "|cos#theta*|"}, phistar{axisPhiStar, "#phi*"};
    constexpr int NEventSelectionSteps = enumValue(EventSelectionStep::NSteps) - 1;
    constexpr int NTrackSelectionSteps = enumValue(TrackSelectionStep::NSteps) - 1;
    constexpr int NK0sSelectionSteps = enumValue(K0sSelectionStep::NSteps) - 1;
    constexpr int NPairSelectionSteps = enumValue(PairSelectionStep::NSteps) - 1;
    registry.add("Event/hEventCutFlow", "event cut flow", HistType::kTH1F, {{NEventSelectionSteps, 0.5, NEventSelectionSteps + 0.5}});
    registry.add("Event/hVertexZBefore", "vertex z before", HistType::kTH1F, {{120, -30., 30.}});
    registry.add("Event/hVertexZAfter", "vertex z after", HistType::kTH1F, {{120, -30., 30.}});
    registry.add("Event/hCentralityBefore", "centrality before", HistType::kTH1F, {{100, 0., 100.}});
    registry.add("Event/hCentralityAfter", "centrality after", HistType::kTH1F, {{100, 0., 100.}});
    registry.add("QA/Tracks/hCutFlow", "daughter cut flow", HistType::kTH1F, {{NTrackSelectionSteps, 0.5, NTrackSelectionSteps + 0.5}});
    registry.add("QA/Tracks/BeforeSelection/hPt", "daughter pT before", HistType::kTH1F, {{200, 0., 20.}});
    registry.add("QA/Tracks/BeforeSelection/hEta", "daughter eta before", HistType::kTH1F, {{100, -1., 1.}});
    registry.add("QA/Tracks/BeforeSelection/hDCAxyVsPt", "DCAxy vs pT before", HistType::kTH2F, {{200, 0., 20.}, {200, -5., 5.}});
    registry.add("QA/Tracks/BeforeSelection/hDCAzVsPt", "DCAz vs pT before", HistType::kTH2F, {{200, 0., 20.}, {200, -5., 5.}});
    registry.add("QA/Tracks/BeforeSelection/hTPCCrossedRowsVsPt", "TPC rows vs pT before", HistType::kTH2F, {{200, 0., 20.}, {170, 0., 170.}});
    registry.add("QA/Tracks/BeforeSelection/hTPCSignalVsP", "TPC signal before", HistType::kTH2F, {{200, 0., 10.}, {300, 0., 300.}});
    registry.add("QA/Tracks/BeforeSelection/hTPCNSigmaPiVsP", "TPC pion PID before", HistType::kTH2F, {{200, 0., 10.}, {120, -6., 6.}});
    registry.add("QA/Tracks/BeforeSelection/hTOFBetaVsP", "TOF beta before", HistType::kTH2F, {{200, 0., 10.}, {140, 0., 1.4}});
    registry.add("QA/Tracks/BeforeSelection/hTOFNSigmaPiVsP", "TOF pion PID before", HistType::kTH2F, {{200, 0., 10.}, {120, -6., 6.}});
    registry.add("QA/Tracks/AfterSelection/hPt", "daughter pT", HistType::kTH1F, {{200, 0., 20.}});
    registry.add("QA/Tracks/AfterSelection/hEta", "daughter eta", HistType::kTH1F, {{100, -1., 1.}});
    registry.add("QA/Tracks/AfterSelection/hDCAxyVsPt", "DCAxy vs pT", HistType::kTH2F, {{200, 0., 20.}, {200, -5., 5.}});
    registry.add("QA/Tracks/AfterSelection/hDCAzVsPt", "DCAz vs pT", HistType::kTH2F, {{200, 0., 20.}, {200, -5., 5.}});
    registry.add("QA/Tracks/AfterSelection/hTPCCrossedRowsVsPt", "TPC rows vs pT", HistType::kTH2F, {{200, 0., 20.}, {170, 0., 170.}});
    registry.add("QA/Tracks/AfterSelection/hTPCSignalVsP", "TPC signal", HistType::kTH2F, {{200, 0., 10.}, {300, 0., 300.}});
    registry.add("QA/Tracks/AfterSelection/hTPCNSigmaPiVsP", "TPC pion PID", HistType::kTH2F, {{200, 0., 10.}, {120, -6., 6.}});
    registry.add("QA/Tracks/AfterSelection/hTOFBetaVsP", "TOF beta", HistType::kTH2F, {{200, 0., 10.}, {140, 0., 1.4}});
    registry.add("QA/Tracks/AfterSelection/hTOFNSigmaPiVsP", "TOF pion PID", HistType::kTH2F, {{200, 0., 10.}, {120, -6., 6.}});
    registry.add("QA/K0s/hCutFlow", "K0S cut flow", HistType::kTH1F, {{NK0sSelectionSteps, 0.5, NK0sSelectionSteps + 0.5}});
    registry.add("QA/K0s/BeforeSelection/hMassVsPt", "K0S before", HistType::kTH2F, {k0mass, k0pt});
    registry.add("QA/K0s/BeforeSelection/hArmenteros", "Armenteros before", HistType::kTH2F, {{100, -1., 1.}, {100, 0., 0.3}});
    registry.add("QA/K0s/BeforeSelection/hCosPA", "cosPA before", HistType::kTH1F, {{200, 0.9, 1.}});
    registry.add("QA/K0s/BeforeSelection/hDCADaughters", "DCA daughters before", HistType::kTH1F, {{200, 0., 2.}});
    registry.add("QA/K0s/BeforeSelection/hCtau", "c tau before", HistType::kTH1F, {{200, 0., 50.}});
    registry.add("QA/K0s/AfterSelection/hMassVsPt", "K0S after (no mass window)", HistType::kTH2F, {k0mass, k0pt});
    registry.add("QA/K0s/AfterSelection/hArmenteros", "Armenteros", HistType::kTH2F, {{100, -1., 1.}, {100, 0., 0.3}});
    registry.add("QA/K0s/AfterSelection/hCosPA", "cosPA", HistType::kTH1F, {{200, 0.9, 1.}});
    registry.add("QA/K0s/AfterSelection/hDCADaughters", "DCA daughters", HistType::kTH1F, {{200, 0., 2.}});
    registry.add("QA/K0s/AfterSelection/hCtau", "c tau", HistType::kTH1F, {{200, 0., 50.}});
    registry.add("QA/K0sK0s/hPairCutFlow", "pair cut flow", HistType::kTH1F, {{NPairSelectionSteps, 0.5, NPairSelectionSteps + 0.5}});
    registry.add("QA/K0sK0s/hPairCutFlowVsCent", "pair cut flow vs centrality", HistType::kTH2F, {{NPairSelectionSteps, 0.5, NPairSelectionSteps + 0.5}, cent});
    registry.add("QA/K0sK0s/SameEvent/hMass", "same-event pair mass", HistType::kTH1F, {pairmass});
    registry.add("QA/K0sK0s/SameEvent/hDeltaR", "same-event DeltaR", HistType::kTH1F, {{200, 0., 5.}});
    registry.add("QA/K0sK0s/SameEvent/hPt1VsPt2", "K0S pT correlation", HistType::kTH2F, {k0pt, k0pt});
    registry.add("QA/K0sK0s/SameEvent/hK0Mass1VsMass2", "K0S mass correlation", HistType::kTH2F, {k0mass, k0mass});
    registry.add("QA/K0sK0s/MixedEvent/hMass", "mixed-event pair mass", HistType::kTH1F, {pairmass});
    registry.add("QA/Flow/QVectors/hQxQyFT0C", "corrected FT0C q2", HistType::kTH2F, {{200, -1., 1.}, {200, -1., 1.}});
    registry.add("QA/Flow/QVectors/hQxQyFT0A", "corrected FT0A q2", HistType::kTH2F, {{200, -1., 1.}, {200, -1., 1.}});
    registry.add("QA/Flow/QVectors/hQxQyTPCpos", "corrected TPCpos q2", HistType::kTH2F, {{200, -1., 1.}, {200, -1., 1.}});
    registry.add("QA/Flow/EventPlane/hPsi2VsCentFT0C", "FT0C event plane", HistType::kTH2F, {cent, {120, -1.5708, 1.5708}});
    registry.add("QA/Flow/Calibration/hCalibrationFlag", "Q-vector calibration flag", HistType::kTH1F, {{5, -0.5, 4.5}});
    registry.add("K0s/hCounts", "K0S counts", HistType::kTHnSparseF, {k0mass, k0pt, cent});
    registry.add("K0s/hCos2Weighted", "K0S EP numerator", HistType::kTHnSparseF, {k0mass, k0pt, cent}, true);
    registry.add("K0s/hSPWeightedFT0C", "K0S SP numerator FT0C", HistType::kTHnSparseF, {k0mass, k0pt, cent}, true);
    registry.add("K0s/hSPWeightedFT0A", "K0S SP numerator FT0A", HistType::kTHnSparseF, {k0mass, k0pt, cent}, true);
    registry.add("K0s/LikeSign/PP/hCounts", "like-sign pi+pi+ candidates", HistType::kTHnSparseF, {k0mass, k0pt, cent});
    registry.add("K0s/LikeSign/MM/hCounts", "like-sign pi-pi- candidates", HistType::kTHnSparseF, {k0mass, k0pt, cent});
    for (auto const& base : {std::string("K0sK0s/SameEvent/"), std::string("K0sK0s/MixedEvent/"), std::string("K0sK0s/Sideband/Combined/"), std::string("K0sK0s/Sideband/LowerLower/"), std::string("K0sK0s/Sideband/LowerUpper/"), std::string("K0sK0s/Sideband/UpperUpper/")}) {
      registry.add((base + "hCounts").c_str(), "K0S-K0S counts", HistType::kTHnSparseF, {pairmass, pairpt, cent, absct, phistar});
      registry.add((base + "hP2Weighted").c_str(), "pair angular moment P2(cos theta*)", HistType::kTHnSparseF, {pairmass, pairpt, cent, absct, phistar}, true);
      registry.add((base + "hP4Weighted").c_str(), "pair angular moment P4(cos theta*)", HistType::kTHnSparseF, {pairmass, pairpt, cent, absct, phistar}, true);
    }
    registry.add("K0sK0s/SameEvent/hCos2Weighted", "pair EP numerator", HistType::kTHnSparseF, {pairmass, pairpt, cent, absct, phistar}, true);
    registry.add("K0sK0s/SameEvent/hSPWeightedFT0C", "pair SP numerator FT0C", HistType::kTHnSparseF, {pairmass, pairpt, cent, absct, phistar}, true);
    registry.add("K0sK0s/SameEvent/hSPWeightedFT0A", "pair SP numerator FT0A", HistType::kTHnSparseF, {pairmass, pairpt, cent, absct, phistar}, true);
    for (auto const& base : {std::string("K0sK0s/LikeSign/PPPP/"), std::string("K0sK0s/LikeSign/MMMM/"), std::string("K0sK0s/LikeSign/PPMM/"), std::string("K0sK0s/LikeSign/Combined/")}) {
      registry.add((base + "hCounts").c_str(), "like-sign K0S-K0S mimic counts", HistType::kTHnSparseF, {pairmass, pairpt, cent, absct, phistar});
      registry.add((base + "hP2Weighted").c_str(), "like-sign pair angular moment P2", HistType::kTHnSparseF, {pairmass, pairpt, cent, absct, phistar}, true);
      registry.add((base + "hP4Weighted").c_str(), "like-sign pair angular moment P4", HistType::kTHnSparseF, {pairmass, pairpt, cent, absct, phistar}, true);
    }
    registry.add("Resolution/SP/pFT0CFT0AVsCent", "FT0C-FT0A SP correlation", HistType::kTProfile, {cent});
    registry.add("Resolution/SP/pFT0CTPCposVsCent", "FT0C-TPCpos SP correlation", HistType::kTProfile, {cent});
    registry.add("Resolution/SP/pFT0ATPCposVsCent", "FT0A-TPCpos SP correlation", HistType::kTProfile, {cent});
    registry.add("Resolution/EP/pCos2PsiFT0CFT0AVsCent", "FT0C-FT0A EP correlation", HistType::kTProfile, {cent});
    registry.add("Resolution/EP/pCos2PsiFT0CTPCposVsCent", "FT0C-TPCpos EP correlation", HistType::kTProfile, {cent});
    registry.add("Resolution/EP/pCos2PsiFT0ATPCposVsCent", "FT0A-TPCpos EP correlation", HistType::kTProfile, {cent});

    constexpr std::array<const char*, NEventSelectionSteps> EventSelectionLabels{"all", "vertex", "sel8", "no same-bunch pileup", "no TF border", "no ITS ROF border", "centrality"};
    constexpr std::array<const char*, NTrackSelectionSteps> TrackSelectionLabels{"all", "kinematics", "TPC clusters", "track quality", "TPC PID", "TOF PID"};
    constexpr std::array<const char*, NK0sSelectionSteps> K0sSelectionLabels{"all", "daughters", "kinematics", "topology", "daughter DCA", "decay geometry", "Armenteros", "Lambda rejection"};
    constexpr std::array<const char*, NPairSelectionSteps> PairSelectionLabels{"all", "no shared daughters", "selected K0S", "close pair", "kinematics", "signal mass"};
    auto eventCutFlow = registry.get<TH1>(HIST("Event/hEventCutFlow"));
    auto trackCutFlow = registry.get<TH1>(HIST("QA/Tracks/hCutFlow"));
    auto k0sCutFlow = registry.get<TH1>(HIST("QA/K0s/hCutFlow"));
    auto pairCutFlow = registry.get<TH1>(HIST("QA/K0sK0s/hPairCutFlow"));
    auto pairCutFlowVsCent = registry.get<TH2>(HIST("QA/K0sK0s/hPairCutFlowVsCent"));
    for (size_t index = 0; index < EventSelectionLabels.size(); ++index)
      eventCutFlow->GetXaxis()->SetBinLabel(index + 1, EventSelectionLabels[index]);
    for (size_t index = 0; index < TrackSelectionLabels.size(); ++index)
      trackCutFlow->GetXaxis()->SetBinLabel(index + 1, TrackSelectionLabels[index]);
    for (size_t index = 0; index < K0sSelectionLabels.size(); ++index)
      k0sCutFlow->GetXaxis()->SetBinLabel(index + 1, K0sSelectionLabels[index]);
    for (size_t index = 0; index < PairSelectionLabels.size(); ++index) {
      pairCutFlow->GetXaxis()->SetBinLabel(index + 1, PairSelectionLabels[index]);
      pairCutFlowVsCent->GetXaxis()->SetBinLabel(index + 1, PairSelectionLabels[index]);
    }
  }

  SidebandRegion sidebandRegion(float mass) const
  {
    if (mass >= k0Cuts.k0sLowerSidebandMin.value && mass < k0Cuts.k0sLowerSidebandMax.value)
      return SidebandRegion::Lower;
    if (mass >= k0Cuts.k0sUpperSidebandMin.value && mass < k0Cuts.k0sUpperSidebandMax.value)
      return SidebandRegion::Upper;
    return SidebandRegion::Outside;
  }

  void fillSidebandPair(ROOT::Math::PxPyPzMVector const& d1, ROOT::Math::PxPyPzMVector const& m, float cent, SidebandRegion side1, SidebandRegion side2)
  {
    auto ang = decayAngles(d1, m);
    const double cosTheta2 = ang[0] * ang[0];
    const double p2Weight = 0.5 * (3. * cosTheta2 - 1.);
    const double p4Weight = 0.125 * (35. * cosTheta2 * cosTheta2 - 30. * cosTheta2 + 3.);
    registry.fill(HIST("K0sK0s/Sideband/Combined/hCounts"), m.M(), m.Pt(), cent, ang[0], ang[1]);
    registry.fill(HIST("K0sK0s/Sideband/Combined/hP2Weighted"), m.M(), m.Pt(), cent, ang[0], ang[1], p2Weight);
    registry.fill(HIST("K0sK0s/Sideband/Combined/hP4Weighted"), m.M(), m.Pt(), cent, ang[0], ang[1], p4Weight);
    if (side1 == SidebandRegion::Lower && side2 == SidebandRegion::Lower) {
      registry.fill(HIST("K0sK0s/Sideband/LowerLower/hCounts"), m.M(), m.Pt(), cent, ang[0], ang[1]);
      registry.fill(HIST("K0sK0s/Sideband/LowerLower/hP2Weighted"), m.M(), m.Pt(), cent, ang[0], ang[1], p2Weight);
      registry.fill(HIST("K0sK0s/Sideband/LowerLower/hP4Weighted"), m.M(), m.Pt(), cent, ang[0], ang[1], p4Weight);
    } else if (side1 == SidebandRegion::Upper && side2 == SidebandRegion::Upper) {
      registry.fill(HIST("K0sK0s/Sideband/UpperUpper/hCounts"), m.M(), m.Pt(), cent, ang[0], ang[1]);
      registry.fill(HIST("K0sK0s/Sideband/UpperUpper/hP2Weighted"), m.M(), m.Pt(), cent, ang[0], ang[1], p2Weight);
      registry.fill(HIST("K0sK0s/Sideband/UpperUpper/hP4Weighted"), m.M(), m.Pt(), cent, ang[0], ang[1], p4Weight);
    } else {
      registry.fill(HIST("K0sK0s/Sideband/LowerUpper/hCounts"), m.M(), m.Pt(), cent, ang[0], ang[1]);
      registry.fill(HIST("K0sK0s/Sideband/LowerUpper/hP2Weighted"), m.M(), m.Pt(), cent, ang[0], ang[1], p2Weight);
      registry.fill(HIST("K0sK0s/Sideband/LowerUpper/hP4Weighted"), m.M(), m.Pt(), cent, ang[0], ang[1], p4Weight);
    }
  }

  template <typename C>
  bool selectEvent(C const& c, bool fill)
  {
    float cent = c.centFT0C();
    if (fill) {
      registry.fill(HIST("Event/hVertexZBefore"), c.posZ());
      registry.fill(HIST("Event/hCentralityBefore"), cent);
      registry.fill(HIST("Event/hEventCutFlow"), enumValue(EventSelectionStep::All));
    }
    if (std::abs(c.posZ()) > eventCuts.maxAbsZvtx)
      return false;
    if (fill)
      registry.fill(HIST("Event/hEventCutFlow"), enumValue(EventSelectionStep::Vertex));
    if (eventCuts.requireSel8 && !c.sel8())
      return false;
    if (fill)
      registry.fill(HIST("Event/hEventCutFlow"), enumValue(EventSelectionStep::Sel8));
    if (eventCuts.requireNoSameBunchPileup && !c.selection_bit(aod::evsel::kNoSameBunchPileup))
      return false;
    if (fill)
      registry.fill(HIST("Event/hEventCutFlow"), enumValue(EventSelectionStep::NoSameBunchPileup));
    if (eventCuts.requireNoTimeFrameBorder && !c.selection_bit(aod::evsel::kNoTimeFrameBorder))
      return false;
    if (fill)
      registry.fill(HIST("Event/hEventCutFlow"), enumValue(EventSelectionStep::NoTimeFrameBorder));
    if (eventCuts.requireNoITSROFrameBorder && !c.selection_bit(aod::evsel::kNoITSROFrameBorder))
      return false;
    if (fill)
      registry.fill(HIST("Event/hEventCutFlow"), enumValue(EventSelectionStep::NoItsRofBorder));
    if (cent < eventCuts.minCentrality || cent >= eventCuts.maxCentrality)
      return false;
    if (fill) {
      registry.fill(HIST("Event/hEventCutFlow"), enumValue(EventSelectionStep::Centrality));
      registry.fill(HIST("Event/hVertexZAfter"), c.posZ());
      registry.fill(HIST("Event/hCentralityAfter"), cent);
    }
    return true;
  }

  template <typename T>
  bool selectDaughter(T const& t, bool qa)
  {
    bool hasTOF = t.hasTOF();
    if (qa) {
      registry.fill(HIST("QA/Tracks/hCutFlow"), enumValue(TrackSelectionStep::All));
      registry.fill(HIST("QA/Tracks/BeforeSelection/hPt"), t.pt());
      registry.fill(HIST("QA/Tracks/BeforeSelection/hEta"), t.eta());
      registry.fill(HIST("QA/Tracks/BeforeSelection/hDCAxyVsPt"), t.pt(), t.dcaXY());
      registry.fill(HIST("QA/Tracks/BeforeSelection/hDCAzVsPt"), t.pt(), t.dcaZ());
      registry.fill(HIST("QA/Tracks/BeforeSelection/hTPCCrossedRowsVsPt"), t.pt(), t.tpcNClsCrossedRows());
      registry.fill(HIST("QA/Tracks/BeforeSelection/hTPCSignalVsP"), t.p(), t.tpcSignal());
      registry.fill(HIST("QA/Tracks/BeforeSelection/hTPCNSigmaPiVsP"), t.p(), t.tpcNSigmaPi());
      if (hasTOF) {
        registry.fill(HIST("QA/Tracks/BeforeSelection/hTOFBetaVsP"), t.p(), t.beta());
        registry.fill(HIST("QA/Tracks/BeforeSelection/hTOFNSigmaPiVsP"), t.p(), t.tofNSigmaPi());
      }
    }
    if (t.pt() < daughterCuts.daughterMinPt || std::abs(t.eta()) > daughterCuts.maxAbsEta)
      return false;
    if (qa)
      registry.fill(HIST("QA/Tracks/hCutFlow"), enumValue(TrackSelectionStep::Kinematics));
    if (t.tpcNClsCrossedRows() < daughterCuts.minTPCRows || t.tpcCrossedRowsOverFindableCls() < daughterCuts.minTPCRowsOverFindable)
      return false;
    if (qa)
      registry.fill(HIST("QA/Tracks/hCutFlow"), enumValue(TrackSelectionStep::TpcClusters));
    if (t.tpcChi2NCl() > daughterCuts.maxTPCChi2NCl || t.itsChi2NCl() > daughterCuts.maxITSChi2NCl)
      return false;
    if (daughterCuts.minITSClusters > 0 && t.itsNCls() < daughterCuts.minITSClusters)
      return false;
    if (qa)
      registry.fill(HIST("QA/Tracks/hCutFlow"), enumValue(TrackSelectionStep::TrackQuality));
    if (std::abs(t.tpcNSigmaPi()) > daughterCuts.maxTPCNSigmaPi)
      return false;
    if (qa)
      registry.fill(HIST("QA/Tracks/hCutFlow"), enumValue(TrackSelectionStep::TpcPid));
    const auto pidMode = static_cast<PionPidMode>(daughterCuts.pidMode.value);
    bool requireTOF = pidMode == PionPidMode::TpcAndTof || (pidMode == PionPidMode::TofRequiredAboveMomentum && t.p() >= daughterCuts.tofRequiredAboveP);
    bool useTOF = requireTOF || (pidMode == PionPidMode::TpcWithTofIfAvailable && hasTOF);
    if (requireTOF && !hasTOF)
      return false;
    if (useTOF && std::abs(t.tofNSigmaPi()) > daughterCuts.maxTOFNSigmaPi)
      return false;
    if (qa)
      registry.fill(HIST("QA/Tracks/hCutFlow"), enumValue(TrackSelectionStep::TofPid));
    if (qa) {
      registry.fill(HIST("QA/Tracks/AfterSelection/hPt"), t.pt());
      registry.fill(HIST("QA/Tracks/AfterSelection/hEta"), t.eta());
      registry.fill(HIST("QA/Tracks/AfterSelection/hDCAxyVsPt"), t.pt(), t.dcaXY());
      registry.fill(HIST("QA/Tracks/AfterSelection/hDCAzVsPt"), t.pt(), t.dcaZ());
      registry.fill(HIST("QA/Tracks/AfterSelection/hTPCCrossedRowsVsPt"), t.pt(), t.tpcNClsCrossedRows());
      registry.fill(HIST("QA/Tracks/AfterSelection/hTPCSignalVsP"), t.p(), t.tpcSignal());
      registry.fill(HIST("QA/Tracks/AfterSelection/hTPCNSigmaPiVsP"), t.p(), t.tpcNSigmaPi());
      if (hasTOF) {
        registry.fill(HIST("QA/Tracks/AfterSelection/hTOFBetaVsP"), t.p(), t.beta());
        registry.fill(HIST("QA/Tracks/AfterSelection/hTOFNSigmaPiVsP"), t.p(), t.tofNSigmaPi());
      }
    }
    return true;
  }

  template <typename C, typename V>
  bool selectK0(C const& c, V const& v, DaughterTracks const&, bool qa)
  {
    float ctau = v.distovertotmom(c.posX(), c.posY(), c.posZ()) * o2::constants::physics::MassK0Short;
    if (qa) {
      registry.fill(HIST("QA/K0s/hCutFlow"), enumValue(K0sSelectionStep::All));
      registry.fill(HIST("QA/K0s/BeforeSelection/hMassVsPt"), v.mK0Short(), v.pt());
      registry.fill(HIST("QA/K0s/BeforeSelection/hArmenteros"), v.alpha(), v.qtarm());
      registry.fill(HIST("QA/K0s/BeforeSelection/hCosPA"), v.v0cosPA());
      registry.fill(HIST("QA/K0s/BeforeSelection/hDCADaughters"), v.dcaV0daughters());
      registry.fill(HIST("QA/K0s/BeforeSelection/hCtau"), ctau);
    }
    auto pos = v.template posTrack_as<DaughterTracks>();
    auto neg = v.template negTrack_as<DaughterTracks>();
    if (!selectDaughter(pos, qa))
      return false;
    if (!selectDaughter(neg, qa))
      return false;
    if (qa)
      registry.fill(HIST("QA/K0s/hCutFlow"), enumValue(K0sSelectionStep::Daughters));
    if (v.pt() < k0Cuts.k0sMinPt || v.pt() >= k0Cuts.k0sMaxPt || std::abs(v.rapidity(0)) >= k0Cuts.k0sMaxAbsY)
      return false;
    if (qa)
      registry.fill(HIST("QA/K0s/hCutFlow"), enumValue(K0sSelectionStep::Kinematics));
    if (v.v0cosPA() < k0Cuts.minCosPA || v.dcaV0daughters() > k0Cuts.maxDCADaughters)
      return false;
    if (qa)
      registry.fill(HIST("QA/K0s/hCutFlow"), enumValue(K0sSelectionStep::Topology));
    if (std::abs(v.dcapostopv()) < k0Cuts.minDCAPosToPV || std::abs(v.dcanegtopv()) < k0Cuts.minDCANegToPV)
      return false;
    if (qa)
      registry.fill(HIST("QA/K0s/hCutFlow"), enumValue(K0sSelectionStep::DaughterDca));
    if (v.v0radius() < k0Cuts.minRadius || v.v0radius() >= k0Cuts.maxRadius)
      return false;
    if (std::abs(ctau) > k0Cuts.maxCtau)
      return false;
    if (qa)
      registry.fill(HIST("QA/K0s/hCutFlow"), enumValue(K0sSelectionStep::DecayGeometry));
    if (k0Cuts.useArmenterosCut && v.qtarm() <= k0Cuts.armenterosSlope * std::abs(v.alpha()))
      return false;
    if (qa)
      registry.fill(HIST("QA/K0s/hCutFlow"), enumValue(K0sSelectionStep::Armenteros));
    if (k0Cuts.rejectLambda && (std::abs(v.mLambda() - o2::constants::physics::MassLambda0) < k0Cuts.lambdaRejectionWindow || std::abs(v.mAntiLambda() - o2::constants::physics::MassLambda0) < k0Cuts.lambdaRejectionWindow))
      return false;
    if (qa) {
      registry.fill(HIST("QA/K0s/hCutFlow"), enumValue(K0sSelectionStep::LambdaRejection));
      registry.fill(HIST("QA/K0s/AfterSelection/hMassVsPt"), v.mK0Short(), v.pt());
      registry.fill(HIST("QA/K0s/AfterSelection/hArmenteros"), v.alpha(), v.qtarm());
      registry.fill(HIST("QA/K0s/AfterSelection/hCosPA"), v.v0cosPA());
      registry.fill(HIST("QA/K0s/AfterSelection/hDCADaughters"), v.dcaV0daughters());
      registry.fill(HIST("QA/K0s/AfterSelection/hCtau"), ctau);
    }
    return true;
  }

  template <typename C>
  bool getQ(C const& c, int det, float& qx, float& qy) const
  {
    const int rawIndex = det * QVectorComponents + flow.correctionStep.value + (flow.harmonic.value - 2) * flow.numberOfDetectors.value * QVectorComponents;
    if (rawIndex < 0 || det < 0)
      return false;
    const auto index = static_cast<size_t>(rawIndex);
    const auto detectorIndex = static_cast<size_t>(det);
    if (index >= c.qvecRe().size() || index >= c.qvecIm().size() || detectorIndex >= c.qvecAmp().size() || c.qvecAmp()[detectorIndex] < flow.minQAmplitude.value)
      return false;
    qx = c.qvecRe()[index];
    qy = c.qvecIm()[index];
    return std::isfinite(qx) && std::isfinite(qy);
  }

  std::array<double, 2> decayAngles(ROOT::Math::PxPyPzMVector const& d, ROOT::Math::PxPyPzMVector const& mother) const
  {
    ROOT::Math::Boost boost{mother.BoostToCM()};
    auto rest = boost(d);
    ROOT::Math::XYZVector z(mother.Px(), mother.Py(), mother.Pz());
    if (z.R() == 0.)
      z.SetZ(1.);
    z = z.Unit();
    ROOT::Math::XYZVector beam(0., 0., 1.);
    ROOT::Math::XYZVector y = beam.Cross(z);
    if (y.R() < VectorTolerance)
      y.SetY(1.);
    y = y.Unit();
    ROOT::Math::XYZVector x = y.Cross(z).Unit();
    ROOT::Math::XYZVector k(rest.Px(), rest.Py(), rest.Pz());
    k = k.Unit();
    return {std::abs(k.Dot(z)), std::atan2(k.Dot(y), k.Dot(x))};
  }

  template <typename C>
  void fillResolution(C const& c, float cent, float cX, float cY, float aX, float aY, float tX, float tY)
  {
    auto dot = [](float x1, float y1, float x2, float y2) { return x1 * x2 + y1 * y2; };
    double psiC = std::atan2(cY, cX) / flow.harmonic, psiA = std::atan2(aY, aX) / flow.harmonic, psiT = std::atan2(tY, tX) / flow.harmonic;
    registry.fill(HIST("Resolution/SP/pFT0CFT0AVsCent"), cent, dot(cX, cY, aX, aY));
    registry.fill(HIST("Resolution/SP/pFT0CTPCposVsCent"), cent, dot(cX, cY, tX, tY));
    registry.fill(HIST("Resolution/SP/pFT0ATPCposVsCent"), cent, dot(aX, aY, tX, tY));
    registry.fill(HIST("Resolution/EP/pCos2PsiFT0CFT0AVsCent"), cent, std::cos(flow.harmonic * (psiC - psiA)));
    registry.fill(HIST("Resolution/EP/pCos2PsiFT0CTPCposVsCent"), cent, std::cos(flow.harmonic * (psiC - psiT)));
    registry.fill(HIST("Resolution/EP/pCos2PsiFT0ATPCposVsCent"), cent, std::cos(flow.harmonic * (psiA - psiT)));
    registry.fill(HIST("QA/Flow/QVectors/hQxQyFT0C"), cX, cY);
    registry.fill(HIST("QA/Flow/QVectors/hQxQyFT0A"), aX, aY);
    registry.fill(HIST("QA/Flow/QVectors/hQxQyTPCpos"), tX, tY);
    registry.fill(HIST("QA/Flow/EventPlane/hPsi2VsCentFT0C"), cent, psiC);
    registry.fill(HIST("QA/Flow/Calibration/hCalibrationFlag"), c.isCalibrated());
  }

  struct LikeSignCandidate {
    ROOT::Math::PxPyPzMVector p4{};
    uint64_t firstTrackId{};
    uint64_t secondTrackId{};
  };

  template <typename C, typename T>
  bool buildLikeSignCandidate(C const& c, T const& t1, T const& t2, LikeSignCandidate& candidate)
  {
    auto track1 = getTrackParCov(t1);
    auto track2 = getTrackParCov(t2);
    int nCandidates = 0;
    try {
      nCandidates = lsFitter.process(track1, track2);
    } catch (std::runtime_error const&) {
      return false;
    }
    if (nCandidates == 0)
      return false;
    auto const& vertex = lsFitter.getPCACandidate();
    std::array<float, 3> p1{};
    std::array<float, 3> p2{};
    lsFitter.getTrack(0).getPxPyPzGlo(p1);
    lsFitter.getTrack(1).getPxPyPzGlo(p2);
    ROOT::Math::PxPyPzMVector pion1(p1[0], p1[1], p1[2], o2::constants::physics::MassPionCharged);
    ROOT::Math::PxPyPzMVector pion2(p2[0], p2[1], p2[2], o2::constants::physics::MassPionCharged);
    auto pair = pion1 + pion2;
    if (pair.Pt() < k0Cuts.k0sMinPt || pair.Pt() >= k0Cuts.k0sMaxPt || std::abs(pair.Rapidity()) >= k0Cuts.k0sMaxAbsY)
      return false;
    const std::array<float, 3> pv{c.posX(), c.posY(), c.posZ()};
    const std::array<float, 3> sv{static_cast<float>(vertex[0]), static_cast<float>(vertex[1]), static_cast<float>(vertex[2])};
    const std::array<float, 3> momentum{static_cast<float>(pair.Px()), static_cast<float>(pair.Py()), static_cast<float>(pair.Pz())};
    const float cosPA = RecoDecay::cpa(pv, sv, momentum);
    const float dcaDaughters = lsFitter.getChi2AtPCACandidate();
    if (cosPA < k0Cuts.minCosPA || dcaDaughters > k0Cuts.maxDCADaughters)
      return false;
    if (std::abs(t1.dcaXY()) < k0Cuts.minDCAPosToPV || std::abs(t2.dcaXY()) < k0Cuts.minDCANegToPV)
      return false;
    const float radius = std::hypot(vertex[0], vertex[1]);
    if (radius < k0Cuts.minRadius || radius >= k0Cuts.maxRadius)
      return false;
    const float decayLength = std::sqrt((vertex[0] - c.posX()) * (vertex[0] - c.posX()) + (vertex[1] - c.posY()) * (vertex[1] - c.posY()) + (vertex[2] - c.posZ()) * (vertex[2] - c.posZ()));
    const float ctau = decayLength * o2::constants::physics::MassK0Short / pair.P();
    if (std::abs(ctau) > k0Cuts.maxCtau)
      return false;
    ROOT::Math::XYZVector total(pair.Px(), pair.Py(), pair.Pz());
    if (total.R() <= 0.)
      return false;
    total = total.Unit();
    ROOT::Math::XYZVector daughter1(p1[0], p1[1], p1[2]);
    ROOT::Math::XYZVector daughter2(p2[0], p2[1], p2[2]);
    const double longitudinal1 = daughter1.Dot(total);
    const double longitudinal2 = daughter2.Dot(total);
    const double denominator = longitudinal1 + longitudinal2;
    if (std::abs(denominator) < VectorTolerance)
      return false;
    const double alpha = (longitudinal1 - longitudinal2) / denominator;
    const double qt = (daughter1 - longitudinal1 * total).R();
    if (k0Cuts.useArmenterosCut && qt <= k0Cuts.armenterosSlope * std::abs(alpha))
      return false;
    candidate = {pair, static_cast<uint64_t>(t1.globalIndex()), static_cast<uint64_t>(t2.globalIndex())};
    return true;
  }

  void fillLikeSignPair(LikeSignCandidate const& first, LikeSignCandidate const& second, float cent, LikeSignCategory category)
  {
    if (first.firstTrackId == second.firstTrackId || first.firstTrackId == second.secondTrackId || first.secondTrackId == second.firstTrackId || first.secondTrackId == second.secondTrackId)
      return;
    ROOT::Math::PxPyPzMVector d1(first.p4.Px(), first.p4.Py(), first.p4.Pz(), o2::constants::physics::MassK0Short);
    ROOT::Math::PxPyPzMVector d2(second.p4.Px(), second.p4.Py(), second.p4.Pz(), o2::constants::physics::MassK0Short);
    auto pair = d1 + d2;
    if (std::abs(pair.Rapidity()) >= pairCuts.pairMaxAbsY || pair.Pt() < pairCuts.pairMinPt || pair.Pt() >= pairCuts.pairMaxPt || pair.M() < pairCuts.pairMassMin || pair.M() >= pairCuts.pairMassMax)
      return;
    const auto angles = decayAngles(d1, pair);
    const double cosTheta2 = angles[0] * angles[0];
    const double p2Weight = 0.5 * (3. * cosTheta2 - 1.);
    const double p4Weight = 0.125 * (35. * cosTheta2 * cosTheta2 - 30. * cosTheta2 + 3.);
    registry.fill(HIST("K0sK0s/LikeSign/Combined/hCounts"), pair.M(), pair.Pt(), cent, angles[0], angles[1]);
    registry.fill(HIST("K0sK0s/LikeSign/Combined/hP2Weighted"), pair.M(), pair.Pt(), cent, angles[0], angles[1], p2Weight);
    registry.fill(HIST("K0sK0s/LikeSign/Combined/hP4Weighted"), pair.M(), pair.Pt(), cent, angles[0], angles[1], p4Weight);
    if (category == LikeSignCategory::PositivePositive) {
      registry.fill(HIST("K0sK0s/LikeSign/PPPP/hCounts"), pair.M(), pair.Pt(), cent, angles[0], angles[1]);
      registry.fill(HIST("K0sK0s/LikeSign/PPPP/hP2Weighted"), pair.M(), pair.Pt(), cent, angles[0], angles[1], p2Weight);
      registry.fill(HIST("K0sK0s/LikeSign/PPPP/hP4Weighted"), pair.M(), pair.Pt(), cent, angles[0], angles[1], p4Weight);
    } else if (category == LikeSignCategory::NegativeNegative) {
      registry.fill(HIST("K0sK0s/LikeSign/MMMM/hCounts"), pair.M(), pair.Pt(), cent, angles[0], angles[1]);
      registry.fill(HIST("K0sK0s/LikeSign/MMMM/hP2Weighted"), pair.M(), pair.Pt(), cent, angles[0], angles[1], p2Weight);
      registry.fill(HIST("K0sK0s/LikeSign/MMMM/hP4Weighted"), pair.M(), pair.Pt(), cent, angles[0], angles[1], p4Weight);
    } else {
      registry.fill(HIST("K0sK0s/LikeSign/PPMM/hCounts"), pair.M(), pair.Pt(), cent, angles[0], angles[1]);
      registry.fill(HIST("K0sK0s/LikeSign/PPMM/hP2Weighted"), pair.M(), pair.Pt(), cent, angles[0], angles[1], p2Weight);
      registry.fill(HIST("K0sK0s/LikeSign/PPMM/hP4Weighted"), pair.M(), pair.Pt(), cent, angles[0], angles[1], p4Weight);
    }
  }

  template <typename C>
  void fillLikeSign(C const& c, DaughterTracks const& tracks)
  {
    auto tracksInCollision = tracks.sliceBy(perCollisionTrack, c.globalIndex());
    std::vector<int64_t> positiveTrackRows;
    std::vector<int64_t> negativeTrackRows;
    positiveTrackRows.reserve(tracksInCollision.size());
    negativeTrackRows.reserve(tracksInCollision.size());

    // Apply all inexpensive daughter-quality and PID selections before forming
    // pairs.  Running the DCA fitter over every track pair scales as N^2 and
    // dominated the processing time in high-multiplicity events.
    for (auto const& track : tracksInCollision) {
      if (!selectDaughter(track, false))
        continue;
      const auto localRow = static_cast<int64_t>(track.globalIndex()) - tracks.offset();
      (track.signed1Pt() > 0.f ? positiveTrackRows : negativeTrackRows).push_back(localRow);
    }

    std::vector<LikeSignCandidate> positiveCandidates;
    std::vector<LikeSignCandidate> negativeCandidates;
    auto buildCandidates = [&](std::vector<int64_t> const& selectedRows, std::vector<LikeSignCandidate>& candidates, bool positive) {
      for (size_t i = 0; i < selectedRows.size(); ++i) {
        auto track1 = tracks.iteratorAt(selectedRows[i]);
        for (size_t j = i + 1; j < selectedRows.size(); ++j) {
          auto track2 = tracks.iteratorAt(selectedRows[j]);
          LikeSignCandidate candidate;
          if (!buildLikeSignCandidate(c, track1, track2, candidate))
            continue;
          if (positive)
            registry.fill(HIST("K0s/LikeSign/PP/hCounts"), candidate.p4.M(), candidate.p4.Pt(), c.centFT0C());
          else
            registry.fill(HIST("K0s/LikeSign/MM/hCounts"), candidate.p4.M(), candidate.p4.Pt(), c.centFT0C());
          if (candidate.p4.M() >= k0Cuts.k0sMassMin && candidate.p4.M() < k0Cuts.k0sMassMax)
            candidates.push_back(candidate);
        }
      }
    };
    buildCandidates(positiveTrackRows, positiveCandidates, true);
    buildCandidates(negativeTrackRows, negativeCandidates, false);
    for (size_t i = 0; i < positiveCandidates.size(); ++i)
      for (size_t j = i + 1; j < positiveCandidates.size(); ++j)
        fillLikeSignPair(positiveCandidates[i], positiveCandidates[j], c.centFT0C(), LikeSignCategory::PositivePositive);
    for (size_t i = 0; i < negativeCandidates.size(); ++i)
      for (size_t j = i + 1; j < negativeCandidates.size(); ++j)
        fillLikeSignPair(negativeCandidates[i], negativeCandidates[j], c.centFT0C(), LikeSignCategory::NegativeNegative);
    for (auto const& positive : positiveCandidates)
      for (auto const& negative : negativeCandidates)
        fillLikeSignPair(positive, negative, c.centFT0C(), LikeSignCategory::PositiveNegative);
  }

  template <typename C, typename V>
  void fillSame(C const& c, V const& v0s, DaughterTracks const& tracks)
  {
    float cent = c.centFT0C(), cX, cY, aX, aY, tX, tY;
    bool qok = getQ(c, flow.ft0cId, cX, cY) && getQ(c, flow.ft0aId, aX, aY) && getQ(c, flow.tpcPosId, tX, tY);
    if (qok)
      fillResolution(c, cent, cX, cY, aX, aY, tX, tY);
    for (auto const& v : v0s) {
      if (!selectK0(c, v, tracks, true))
        continue;
      double phi = v.phi(), psi = qok ? std::atan2(cY, cX) / flow.harmonic : 0.;
      double ux = std::cos(flow.harmonic * phi), uy = std::sin(flow.harmonic * phi);
      registry.fill(HIST("K0s/hCounts"), v.mK0Short(), v.pt(), cent);
      if (qok) {
        registry.fill(HIST("K0s/hCos2Weighted"), v.mK0Short(), v.pt(), cent, std::cos(flow.harmonic * (phi - psi)));
        registry.fill(HIST("K0s/hSPWeightedFT0C"), v.mK0Short(), v.pt(), cent, ux * cX + uy * cY);
        registry.fill(HIST("K0s/hSPWeightedFT0A"), v.mK0Short(), v.pt(), cent, ux * aX + uy * aY);
      }
    }
    for (auto const& [v1, v2] : combinations(CombinationsStrictlyUpperIndexPolicy(v0s, v0s))) {
      registry.fill(HIST("QA/K0sK0s/hPairCutFlow"), enumValue(PairSelectionStep::All));
      registry.fill(HIST("QA/K0sK0s/hPairCutFlowVsCent"), enumValue(PairSelectionStep::All), cent);
      if (v1.posTrackId() == v2.posTrackId() || v1.posTrackId() == v2.negTrackId() || v1.negTrackId() == v2.posTrackId() || v1.negTrackId() == v2.negTrackId())
        continue;
      registry.fill(HIST("QA/K0sK0s/hPairCutFlow"), enumValue(PairSelectionStep::NoSharedDaughters));
      registry.fill(HIST("QA/K0sK0s/hPairCutFlowVsCent"), enumValue(PairSelectionStep::NoSharedDaughters), cent);
      if (!selectK0(c, v1, tracks, false) || !selectK0(c, v2, tracks, false))
        continue;
      registry.fill(HIST("QA/K0sK0s/hPairCutFlow"), enumValue(PairSelectionStep::SelectedK0s));
      registry.fill(HIST("QA/K0sK0s/hPairCutFlowVsCent"), enumValue(PairSelectionStep::SelectedK0s), cent);
      ROOT::Math::PxPyPzMVector d1(v1.px(), v1.py(), v1.pz(), o2::constants::physics::MassK0Short), d2(v2.px(), v2.py(), v2.pz(), o2::constants::physics::MassK0Short), m = d1 + d2;
      double dr = std::hypot(v1.eta() - v2.eta(), TVector2::Phi_mpi_pi(v1.phi() - v2.phi()));
      if (pairCuts.enableClosePairRejection && dr < pairCuts.minDeltaR)
        continue;
      registry.fill(HIST("QA/K0sK0s/hPairCutFlow"), enumValue(PairSelectionStep::ClosePair));
      registry.fill(HIST("QA/K0sK0s/hPairCutFlowVsCent"), enumValue(PairSelectionStep::ClosePair), cent);
      if (std::abs(m.Rapidity()) >= pairCuts.pairMaxAbsY || m.Pt() < pairCuts.pairMinPt || m.Pt() >= pairCuts.pairMaxPt || m.M() < pairCuts.pairMassMin || m.M() >= pairCuts.pairMassMax)
        continue;
      registry.fill(HIST("QA/K0sK0s/hPairCutFlow"), enumValue(PairSelectionStep::Kinematics));
      registry.fill(HIST("QA/K0sK0s/hPairCutFlowVsCent"), enumValue(PairSelectionStep::Kinematics), cent);
      const bool signalPair = v1.mK0Short() >= k0Cuts.k0sMassMin && v1.mK0Short() < k0Cuts.k0sMassMax && v2.mK0Short() >= k0Cuts.k0sMassMin && v2.mK0Short() < k0Cuts.k0sMassMax;
      const auto side1 = sidebandRegion(v1.mK0Short());
      const auto side2 = sidebandRegion(v2.mK0Short());
      if (side1 != SidebandRegion::Outside && side2 != SidebandRegion::Outside)
        fillSidebandPair(d1, m, cent, side1, side2);
      if (!signalPair)
        continue;
      registry.fill(HIST("QA/K0sK0s/hPairCutFlow"), enumValue(PairSelectionStep::SignalMass));
      registry.fill(HIST("QA/K0sK0s/hPairCutFlowVsCent"), enumValue(PairSelectionStep::SignalMass), cent);
      registry.fill(HIST("QA/K0sK0s/SameEvent/hMass"), m.M());
      registry.fill(HIST("QA/K0sK0s/SameEvent/hDeltaR"), dr);
      registry.fill(HIST("QA/K0sK0s/SameEvent/hPt1VsPt2"), v1.pt(), v2.pt());
      registry.fill(HIST("QA/K0sK0s/SameEvent/hK0Mass1VsMass2"), v1.mK0Short(), v2.mK0Short());
      auto ang = decayAngles(d1, m);
      registry.fill(HIST("K0sK0s/SameEvent/hCounts"), m.M(), m.Pt(), cent, ang[0], ang[1]);
      const double cosTheta = ang[0];
      const double cosTheta2 = cosTheta * cosTheta;
      const double p2Weight = 0.5 * (3. * cosTheta2 - 1.);
      const double p4Weight = 0.125 * (35. * cosTheta2 * cosTheta2 - 30. * cosTheta2 + 3.);
      registry.fill(HIST("K0sK0s/SameEvent/hP2Weighted"), m.M(), m.Pt(), cent, ang[0], ang[1], p2Weight);
      registry.fill(HIST("K0sK0s/SameEvent/hP4Weighted"), m.M(), m.Pt(), cent, ang[0], ang[1], p4Weight);
      if (qok) {
        double ux = std::cos(flow.harmonic * m.Phi()), uy = std::sin(flow.harmonic * m.Phi()), psi = std::atan2(cY, cX) / flow.harmonic;
        registry.fill(HIST("K0sK0s/SameEvent/hCos2Weighted"), m.M(), m.Pt(), cent, ang[0], ang[1], std::cos(flow.harmonic * (m.Phi() - psi)));
        registry.fill(HIST("K0sK0s/SameEvent/hSPWeightedFT0C"), m.M(), m.Pt(), cent, ang[0], ang[1], ux * cX + uy * cY);
        registry.fill(HIST("K0sK0s/SameEvent/hSPWeightedFT0A"), m.M(), m.Pt(), cent, ang[0], ang[1], ux * aX + uy * aY);
      }
    }
    if (likeSign.enableLikeSign) {
      fillLikeSign(c, tracks);
    }
  }

  void processSame(Collisions::iterator const& collision, aod::V0Datas const& v0s, DaughterTracks const& tracks)
  {
    if (!selectEvent(collision, true))
      return;
    fillSame(collision, v0s, tracks);
  }
  PROCESS_SWITCH(K0sK0sFlow, processSame, "same-event K0S-K0S and flow", true);

  void processMixed(Collisions const& collisions, aod::V0Datas const& v0s, DaughterTracks const& tracks)
  {
    MixingBinning binning{{axisMixZvtx, axisMixCent}, true};
    for (auto const& [c1, c2] : selfCombinations(binning, pairCuts.mixingDepth, -1, collisions, collisions)) {
      if (c1.index() == c2.index())
        continue;
      if (!selectEvent(c1, false) || !selectEvent(c2, false))
        continue;
      auto s1 = v0s.sliceBy(perCollisionV0, c1.globalIndex());
      auto s2 = v0s.sliceBy(perCollisionV0, c2.globalIndex());
      if (s1.size() == 0 || s2.size() == 0)
        continue;
      for (auto const& [v1, v2] : combinations(CombinationsFullIndexPolicy(s1, s2))) {
        if (!selectK0(c1, v1, tracks, false) || !selectK0(c2, v2, tracks, false) || v1.mK0Short() < k0Cuts.k0sMassMin || v1.mK0Short() >= k0Cuts.k0sMassMax || v2.mK0Short() < k0Cuts.k0sMassMin || v2.mK0Short() >= k0Cuts.k0sMassMax)
          continue;
        ROOT::Math::PxPyPzMVector d1(v1.px(), v1.py(), v1.pz(), o2::constants::physics::MassK0Short), d2(v2.px(), v2.py(), v2.pz(), o2::constants::physics::MassK0Short), m = d1 + d2;
        if (std::abs(m.Rapidity()) >= pairCuts.pairMaxAbsY || m.Pt() < pairCuts.pairMinPt || m.Pt() >= pairCuts.pairMaxPt || m.M() < pairCuts.pairMassMin || m.M() >= pairCuts.pairMassMax)
          continue;
        auto ang = decayAngles(d1, m);
        registry.fill(HIST("K0sK0s/MixedEvent/hCounts"), m.M(), m.Pt(), c1.centFT0C(), ang[0], ang[1]);
        const double cosTheta = ang[0];
        const double cosTheta2 = cosTheta * cosTheta;
        const double p2Weight = 0.5 * (3. * cosTheta2 - 1.);
        const double p4Weight = 0.125 * (35. * cosTheta2 * cosTheta2 - 30. * cosTheta2 + 3.);
        registry.fill(HIST("K0sK0s/MixedEvent/hP2Weighted"), m.M(), m.Pt(), c1.centFT0C(), ang[0], ang[1], p2Weight);
        registry.fill(HIST("K0sK0s/MixedEvent/hP4Weighted"), m.M(), m.Pt(), c1.centFT0C(), ang[0], ang[1], p4Weight);
        registry.fill(HIST("QA/K0sK0s/MixedEvent/hMass"), m.M());
      }
    }
  }
  PROCESS_SWITCH(K0sK0sFlow, processMixed, "mixed-event K0S-K0S (optional cross-check)", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<K0sK0sFlow>(cfgc)};
}
