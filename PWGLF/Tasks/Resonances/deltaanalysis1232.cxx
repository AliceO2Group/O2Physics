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

/// \file deltaAnalysis1232.cxx
/// \brief Task for Delta(1232) resonance reconstruction via proton-pion invariant mass analysis.
/// \author Durgesh Bhatt <durgesh.bhatt@cern.ch>

#include "PWGLF/DataModel/LFResonanceTables.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/EventSelection.h"

#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/ASoAHelpers.h>
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
#include <Framework/runDataProcessing.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

namespace delta_analysis
{
static constexpr int PdgProton{2212};
static constexpr int PdgPion{211};
static constexpr int PdgDeltaPlusPlus{2224};
static constexpr int PdgDeltaZero{2114};
static constexpr float MassProton{o2::constants::physics::MassProton};
static constexpr float MassPion{o2::constants::physics::MassPionCharged};
static constexpr float MassDelta{1.232f};
} // namespace delta_analysis

struct DeltaAnalysis1232 {
  SliceCache sliceCache;

  Preslice<aod::ResoTracks> tracksPerResonanceCollision = aod::resodaughter::resoCollisionId;
  Preslice<aod::ResoTrackDFs> tracksPerMergedDFCollision = aod::resodaughter::resoCollisionDFId;

  aod::ResoMCParents const* mcResonanceParentTable = nullptr;

  // ── Event-level configurables ────────────────────────────────────────────
  Configurable<bool> applyOccupancyInTimeRangeCut{"applyOccupancyInTimeRangeCut", false, "If true, apply a cut on the number of tracks in a time window around the collision"};
  Configurable<int> occupancyMin{"occupancyMin", 0, "Minimum track occupancy in time range"};
  Configurable<int> occupancyMax{"occupancyMax", 9999, "Maximum track occupancy in time range"};

  // ── Histogram binning configurables ─────────────────────────────────────
  Configurable<int> numberOfPtBins{"numberOfPtBins", 100, "Number of bins along the transverse momentum (pT) axis"};
  Configurable<int> numberOfInvMassBins{"numberOfInvMassBins", 120, "Number of bins along the invariant mass axis"};

  // ── Physics configurables ────────────────────────────────────────────────
  Configurable<bool> enableRotationalBackground{"enableRotationalBackground", true, "If true, compute rotational background (pion phi rotated by ~pi)"};

  // ── Track quality cuts ───────────────────────────────────────────────────
  Configurable<float> minTrackPt{"minTrackPt", 0.15f, "Minimum transverse momentum pT of a track [GeV/c]"};
  Configurable<float> minProtonMomentum{"minProtonMomentum", 0.f, "Minimum total momentum p for proton PID [GeV/c]"};
  Configurable<float> minPionMomentum{"minPionMomentum", 0.f, "Minimum total momentum p for pion PID [GeV/c]"};
  Configurable<float> minPseudorapidity{"minPseudorapidity", -0.8f, "Minimum pseudorapidity eta (detector acceptance)"};
  Configurable<float> maxPseudorapidity{"maxPseudorapidity", 0.8f, "Maximum pseudorapidity eta (detector acceptance)"};
  Configurable<float> maxDCAz{"maxDCAz", 1.0f, "Maximum allowed distance of closest approach along z-axis (DCAz) [cm]"};
  Configurable<float> minPairRapidity{"minPairRapidity", -0.5f, "Minimum rapidity y of the reconstructed pair"};
  Configurable<float> maxPairRapidity{"maxPairRapidity", 0.5f, "Maximum rapidity y of the reconstructed pair"};

  // ── TPC cluster quality cuts ─────────────────────────────────────────────
  Configurable<int> minTPCCrossedRows{"minTPCCrossedRows", 70, "Minimum number of TPC crossed pad rows"};
  Configurable<bool> applyCrossedRowsCut{"applyCrossedRowsCut", false, "If true, require at least minTPCCrossedRows"};
  Configurable<int> minTPCClustersFound{"minTPCClustersFound", 70, "Minimum number of TPC clusters found on the track"};
  Configurable<bool> applyTPCClustersCut{"applyTPCClustersCut", false, "If true, require at least minTPCClustersFound"};

  // ── pT-dependent DCAxy cuts for Protons ─────────────────────────────────
  Configurable<std::vector<float>> protonDCAPtBinEdges{"protonDCAPtBinEdges", {0.0f, 0.5f, 1.0f, 2.0f, 3.0f, 5.0f, 1000.0f}, "Proton pT bin edges [GeV/c] for DCAxy"};
  Configurable<std::vector<float>> protonMaxDCAxyPerPtBin{"protonMaxDCAxyPerPtBin", {0.10f, 0.08f, 0.05f, 0.03f, 0.02f, 0.015f}, "Maximum |DCAxy| [cm] for protons"};

  // ── pT-dependent DCAxy cuts for Pions ───────────────────────────────────
  Configurable<std::vector<float>> pionDCAPtBinEdges{"pionDCAPtBinEdges", {0.0f, 0.5f, 1.0f, 2.0f, 3.0f, 5.0f, 1000.0f}, "Pion pT bin edges [GeV/c] for DCAxy"};
  Configurable<std::vector<float>> pionMaxDCAxyPerPtBin{"pionMaxDCAxyPerPtBin", {0.20f, 0.15f, 0.10f, 0.08f, 0.05f, 0.03f}, "Maximum |DCAxy| [cm] for pions"};

  // ── Analysis mode switches ───────────────────────────────────────────────
  Configurable<bool> runQualityChecksOnly{"runQualityChecksOnly", false, "If true, only fill QA histograms and skip invariant mass computation"};
  Configurable<bool> applyDeepAngleCut{"applyDeepAngleCut", false, "Reject proton-pion pairs with very small opening angle"};
  Configurable<double> deepAngleCutValue{"deepAngleCutValue", 0.04, "Minimum allowed opening angle [rad] between proton and pion"};

  // ── Global track selection flags ─────────────────────────────────────────
  Configurable<bool> requirePrimaryTrack{"requirePrimaryTrack", true, "Require track to pass the 'isPrimaryTrack' flag"};
  Configurable<bool> requireGlobalTrackNoDCA{"requireGlobalTrackNoDCA", true, "Require track to pass 'isGlobalTrackWoDCA'"};
  Configurable<bool> requirePVContributor{"requirePVContributor", true, "Require track to be a PV contributor"};

  // ── PID configurables ────────────────────────────────────────────────────
  Configurable<bool> requireTOFForProton{"requireTOFForProton", false, "If true, only accept proton candidates that have a TOF measurement"};
  Configurable<bool> requireTOFForPion{"requireTOFForPion", false, "If true, only accept pion candidates that have a TOF measurement"};
  Configurable<bool> useTPCOnlyPID{"useTPCOnlyPID", false, "If true, use only TPC nSigma for PID (ignore TOF even if available)"};

  Configurable<float> tpcNSigmaVetoOtherSpecies{"tpcNSigmaVetoOtherSpecies", 3.0f, "Reject track if its TPC nSigma for a different species is below this threshold"};
  Configurable<float> minTPCNSigmaPion{"minTPCNSigmaPion", -6.0f, "Minimum (most negative) TPC nSigma for pion"};
  Configurable<float> minTPCNSigmaProton{"minTPCNSigmaProton", -6.0f, "Minimum (most negative) TPC nSigma for proton"};
  Configurable<float> minTOFNSigmaPion{"minTOFNSigmaPion", -6.0f, "Minimum (most negative) TOF nSigma for pion"};
  Configurable<float> minTOFNSigmaProton{"minTOFNSigmaProton", -6.0f, "Minimum (most negative) TOF nSigma for proton"};
  Configurable<float> minCombinedTPCTOFNSigmaPion{"minCombinedTPCTOFNSigmaPion", -6.0f, "Minimum combined TPC+TOF nSigma for pion"};
  Configurable<float> minCombinedTPCTOFNSigmaProton{"minCombinedTPCTOFNSigmaProton", -6.0f, "Minimum combined TPC+TOF nSigma for proton"};
  Configurable<float> tpcNSigmaVetoThreshold{"tpcNSigmaVetoThreshold", 3.0f, "General TPC nSigma veto cut to reject misidentified particles"};
  Configurable<float> tofNSigmaVetoThreshold{"tofNSigmaVetoThreshold", 3.0f, "General TOF nSigma veto cut to reject misidentified particles"};

  Configurable<double> maxTPCNSigmaProton{"maxTPCNSigmaProton", 3.0, "Maximum |TPC nSigma| for proton identification"};
  Configurable<double> combinedNSigmaCutProton{"combinedNSigmaCutProton", 3.0, "Cut on sqrt(nSigmaTPC^2 + nSigmaTOF^2) for proton. Negative value switches to asymmetric mode."};
  Configurable<std::vector<float>> protonTPCPIDMomentumBins{"protonTPCPIDMomentumBins", {0, 0.5, 0.7, 0.8}, "Momentum p bin edges [GeV/c] for TPC PID cuts on protons"};
  Configurable<std::vector<float>> protonTPCNSigmaCutPerBin{"protonTPCNSigmaCutPerBin", {5., 3.5, 2.5}, "Maximum TPC nSigma for proton in each momentum bin"};
  Configurable<std::vector<float>> protonTOFPIDMomentumBins{"protonTOFPIDMomentumBins", {0., 999.}, "Momentum p bin edges [GeV/c] for TOF PID cuts on protons"};
  Configurable<std::vector<float>> protonTOFNSigmaCutPerBin{"protonTOFNSigmaCutPerBin", {3.0}, "Maximum TOF nSigma for proton in each momentum bin"};

  Configurable<double> maxTPCNSigmaPion{"maxTPCNSigmaPion", 3.0, "Maximum |TPC nSigma| for pion identification"};
  Configurable<double> combinedNSigmaCutPion{"combinedNSigmaCutPion", 3.0, "Cut on sqrt(nSigmaTPC^2 + nSigmaTOF^2) for pion. Negative value switches to asymmetric mode."};
  Configurable<std::vector<float>> pionTPCPIDMomentumBins{"pionTPCPIDMomentumBins", {0., 0.25, 0.4, 0.5}, "Momentum p bin edges [GeV/c] for TPC PID cuts on pions"};
  Configurable<std::vector<float>> pionTPCNSigmaCutPerBin{"pionTPCNSigmaCutPerBin", {5., 3.5, 2.5}, "Maximum TPC nSigma for pion in each momentum bin"};
  Configurable<std::vector<float>> pionTOFPIDMomentumBins{"pionTOFPIDMomentumBins", {0., 999.}, "Momentum p bin edges [GeV/c] for TOF PID cuts on pions"};
  Configurable<std::vector<float>> pionTOFNSigmaCutPerBin{"pionTOFNSigmaCutPerBin", {3.0}, "Maximum TOF nSigma for pion in each momentum bin"};

  // ── Event mixing configurables ───────────────────────────────────────────
  Configurable<int> numberOfEventsToMix{"numberOfEventsToMix", 20, "Number of events to mix with each signal event"};
  ConfigurableAxis dcaZMixingBins{"dcaZMixingBins", {VARIABLE_WIDTH, -1.2f, -1.0f, -0.9f, -0.8f, -0.7f, -0.6f, -0.5f, -0.4f, -0.3f, -0.2f, -0.1f, 0.f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.2f}, "DCAz vertex bins used for mixing"};
  ConfigurableAxis vertexZMixingBins{"vertexZMixingBins", {VARIABLE_WIDTH, -12.f, -10.f, -9.f, -8.f, -7.f, -6.f, -5.f, -4.f, -3.f, -2.f, -1.f, 0.f, 1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f, 10.f, 12.f}, "Primary vertex z bins used for mixing"};
  ConfigurableAxis centralityMixingBins{"centralityMixingBins", {VARIABLE_WIDTH, 0.f, 10.f, 20.f, 30.f, 40.f, 50.f, 60.f, 70.f, 80.f, 90.f, 100.f, 200.f}, "Centrality bins used for mixing"};
  ConfigurableAxis eventPlaneMixingBins{"eventPlaneMixingBins", {VARIABLE_WIDTH, -1.5708f, -1.25664f, -0.942478f, -0.628319f, 0.f, 0.628319f, 0.942478f, 1.25664f, 1.5708f}, "Event plane angle bins used for mixing"};
  ConfigurableAxis occupancyBins{"occupancyBins", {VARIABLE_WIDTH, 0.0, 100, 500, 600, 1000, 1100, 1500, 1600, 2000, 2100, 2500, 2600, 3000, 3100, 3500, 3600, 4000, 4100, 4500, 4600, 5000, 5100, 9999}, "Track occupancy bins"};

  // ── Rotational background configurables ─────────────────────────────────
  Configurable<int> numberOfRotations{"numberOfRotations", 10, "How many times to rotate the pion phi for rotational background"};
  Configurable<float> rotationAngleWindow{"rotationAngleWindow", 6.f, "Pion rotated by angles near PI, within a window of PI/rotationAngleWindow"};

  // ── MC event selection flags ─────────────────────────────────────────────
  Configurable<bool> mcRequireAfterAllCuts{"mcRequireAfterAllCuts", false, "MC event selection: require isInAfterAllCuts flag"};
  Configurable<bool> mcRequireINELgt0{"mcRequireINELgt0", false, "MC event selection: require at least 1 charged particle in |eta|<1 (INEL>0)"};
  Configurable<bool> mcRequireSel8{"mcRequireSel8", false, "MC event selection: require the standard Sel8 event selection"};
  Configurable<bool> mcRequireVtxWithin10cm{"mcRequireVtxWithin10cm", false, "MC event selection: require primary vertex |z| < 10 cm"};
  Configurable<bool> mcRequireTriggerTVX{"mcRequireTriggerTVX", false, "MC event selection: require the TVX (T0 vertex) trigger"};
  Configurable<bool> cEvtRecINELgt0{"cEvtRecINELgt0", false, "MC event selection: require reconstructed INEL>0 condition"};

  // ── Histogram registry ───────────────────────────────────────────────────
  HistogramRegistry allHistograms{"allHistograms", {}, OutputObjHandlingPolicy::AnalysisObject};

  std::vector<float> precomputedRotations;

  // Framework-safe candidate mapping (prevents ASoA table indexing bugs)
  struct ValidCandidate {
    float px, py, pz, pt, p, phi, eta;
    int sign;
    int motherId;
    int motherPDG;
  };

  void init(InitContext const&)
  {
    if (enableRotationalBackground && numberOfRotations > 0) {
      float rotWindowHalfWidth = o2::constants::math::PI / rotationAngleWindow;
      for (int iRot = 0; iRot < numberOfRotations; iRot++) {
        if (numberOfRotations == 1) {
          precomputedRotations.push_back(o2::constants::math::PI);
        } else {
          precomputedRotations.push_back((o2::constants::math::PI - rotWindowHalfWidth) + iRot * (2.f * rotWindowHalfWidth / (numberOfRotations - 1)));
        }
      }
    }

    const AxisSpec axisCentralityPercent(110, 0, 110, "FT0 centrality (%)");
    const AxisSpec axisMomentumForPID(200, 0., 10., "p (GeV/c)");
    const AxisSpec axisPtForPID(200, 0., 10., "p_{T} (GeV/c)");
    const AxisSpec axisPt(numberOfPtBins, 0., 10., "p_{T} (GeV/c)");
    const AxisSpec axisRapidity(20, -0.7, 0.7, "y");
    const AxisSpec axisPseudorapidity(40, -1, 1, "#eta");
    const AxisSpec axisDCAxy(240, -0.12, 0.12, "DCA_{xy} (cm)");
    const AxisSpec axisTPCNClusters(200, 0, 200, "TPC N_{clusters}");
    const AxisSpec axisTPCNSigma(401, -10.025, 10.025, "n#sigma^{TPC}");
    const AxisSpec axisTOFNSigma(401, -10.025, 10.025, "n#sigma^{TOF}");
    const AxisSpec axisTPCdEdx(380, 10, 200, "#frac{dE}{dx}");
    const AxisSpec axisVertexZ(120, -12, 12, "v_{z} (cm)");
    const AxisSpec axisEventPlaneAngle(120, -3.14, 3.14, "#theta (rad)");
    const AxisSpec axisInvariantMass(numberOfInvMassBins, 0.8, 1.8, "M_{inv} (GeV/c^{2})");
    const AxisSpec axisOpenAngle(100, 0., o2::constants::math::PI, "Opening angle (rad)");

    AxisSpec axisOccupancy = {occupancyBins, "Track occupancy [-40,100]"};
    AxisSpec axisDCAz = {dcaZMixingBins, "DCA_{z} (cm)"};

    allHistograms.add("Event/centralityVsOccupancy", "Collision centrality vs track occupancy", kTH2F, {axisCentralityPercent, axisOccupancy});
    allHistograms.add("Event/primaryVertexZ", "Primary vertex z-position distribution", kTH1F, {{100, -15., 15.}});
    if (doprocessMix || doprocessMixDF || doprocessMixepDF) {
      allHistograms.add("Event/mixingBins_centralityVsVtxZVsEventPlane", "Event mixing bin occupancy", kTH3F, {axisCentralityPercent, axisVertexZ, axisEventPlaneAngle});
    }

    allHistograms.add("QAbefore/trackEta", "Track pseudorapidity (before cuts)", kTH1F, {{50, -1.0, 1.0}});
    allHistograms.add("QAbefore/trackPt", "Track pT (before cuts)", kTH1F, {axisMomentumForPID});
    allHistograms.add("QAbefore/trackPhi", "Track azimuthal angle (before cuts)", kTH1F, {{72, 0, 6.2832}});
    allHistograms.add("QAbefore/trackEtaVsPhi", "Track eta vs phi (before cuts)", kTH2F, {axisPseudorapidity, {72, 0, 6.2832}});

    allHistograms.add("QAbefore/Proton/tpcNSigmaVsMomentum", "TPC nSigma proton vs momentum (before PID cuts)", kTH2F, {axisMomentumForPID, axisTPCNSigma});
    allHistograms.add("QAbefore/Proton/tofNSigmaVsMomentum", "TOF nSigma proton vs momentum (before PID cuts)", kTH2F, {axisMomentumForPID, axisTOFNSigma});
    allHistograms.add("QAbefore/Proton/tofNSigmaVsTPCNSigma", "TOF nSigma vs TPC nSigma proton (before PID cuts)", kTH2F, {axisTPCNSigma, axisTOFNSigma});
    allHistograms.add("QAbefore/Proton/tpcDedxVsMomentum", "TPC dE/dx proton vs momentum (before cuts)", kTH2F, {axisMomentumForPID, axisTPCdEdx});

    allHistograms.add("QAbefore/Pion/tpcNSigmaVsMomentum", "TPC nSigma pion vs momentum (before PID cuts)", kTH2F, {axisMomentumForPID, axisTPCNSigma});
    allHistograms.add("QAbefore/Pion/tofNSigmaVsMomentum", "TOF nSigma pion vs momentum (before PID cuts)", kTH2F, {axisMomentumForPID, axisTOFNSigma});
    allHistograms.add("QAbefore/Pion/tofNSigmaVsTPCNSigma", "TOF nSigma vs TPC nSigma pion (before PID cuts)", kTH2F, {axisTPCNSigma, axisTOFNSigma});
    allHistograms.add("QAbefore/Pion/tpcDedxVsMomentum", "TPC dE/dx pion vs momentum (before cuts)", kTH2F, {axisMomentumForPID, axisTPCdEdx});

    allHistograms.add("QAafter/Proton/ptVsCentrality", "Proton pT vs centrality (after cuts)", kTH2F, {axisPtForPID, axisCentralityPercent});
    allHistograms.add("QAafter/Proton/dcaZVsPt", "Proton DCAz vs pT (after cuts)", kTH2F, {axisPtForPID, axisDCAz});
    allHistograms.add("QAafter/Proton/dcaXYVsPt", "Proton DCAxy vs pT (after cuts)", kTH2F, {axisPtForPID, axisDCAxy});
    allHistograms.add("QAafter/Proton/tpcDedxVsMomentum", "Proton TPC dE/dx signal vs momentum (after cuts)", kTH2F, {axisMomentumForPID, axisTPCdEdx});
    allHistograms.add("QAafter/Proton/tpcNSigmaVsPt", "Proton TPC nSigma vs pT (after cuts)", kTH2F, {axisPtForPID, axisTPCNSigma});
    allHistograms.add("QAafter/Proton/tpcNSigmaPionContamVsPt", "Proton track: TPC nSigma pion contamination check", kTH2F, {axisPtForPID, axisTPCNSigma});
    allHistograms.add("QAafter/Proton/tpcNSigmaKaonContamVsPt", "Proton track: TPC nSigma kaon contamination check", kTH2F, {axisPtForPID, axisTPCNSigma});
    allHistograms.add("QAafter/Proton/tpcNSigmaVsMomentum", "Proton TPC nSigma vs total momentum (after cuts)", kTH2F, {axisMomentumForPID, axisTPCNSigma});
    allHistograms.add("QAafter/Proton/tpcNSigmaVsCentrality", "Proton TPC nSigma vs centrality", kTH2F, {axisCentralityPercent, axisTPCNSigma});
    allHistograms.add("QAafter/Proton/tofNSigmaVsPt", "Proton TOF nSigma vs pT (after cuts)", kTH2F, {axisPtForPID, axisTOFNSigma});
    allHistograms.add("QAafter/Proton/tofNSigmaVsMomentum", "Proton TOF nSigma vs total momentum (after cuts)", kTH2F, {axisMomentumForPID, axisTOFNSigma});
    allHistograms.add("QAafter/Proton/tofNSigmaVsCentrality", "Proton TOF nSigma vs centrality", kTH2F, {axisCentralityPercent, axisTOFNSigma});
    allHistograms.add("QAafter/Proton/tofNSigmaVsTPCNSigma", "Proton TOF nSigma vs TPC nSigma (after cuts)", kTH2F, {axisTPCNSigma, axisTOFNSigma});
    allHistograms.add("QAafter/Proton/tpcCrossedRowsVsPt", "Proton TPC crossed rows vs pT", kTH2F, {axisPtForPID, {200, 0, 200}});
    allHistograms.add("QAafter/Proton/tpcClustersFoundVsPt", "Proton TPC clusters found vs pT", kTH2F, {axisPtForPID, {200, 0, 200}});

    allHistograms.add("QAafter/Pion/ptVsCentrality", "Pion pT vs centrality (after cuts)", kTH2F, {axisPtForPID, axisCentralityPercent});
    allHistograms.add("QAafter/Pion/dcaZVsPt", "Pion DCAz vs pT (after cuts)", kTH2F, {axisPtForPID, axisDCAz});
    allHistograms.add("QAafter/Pion/dcaXYVsPt", "Pion DCAxy vs pT (after cuts)", kTH2F, {axisPtForPID, axisDCAxy});
    allHistograms.add("QAafter/Pion/tpcDedxVsMomentum", "Pion TPC dE/dx signal vs momentum (after cuts)", kTH2F, {axisMomentumForPID, axisTPCdEdx});
    allHistograms.add("QAafter/Pion/tpcNSigmaProtonContamVsMomentum", "Pion track: TPC nSigma proton contamination check", kTH2F, {axisMomentumForPID, axisTPCNSigma});
    allHistograms.add("QAafter/Pion/tpcNSigmaKaonContamVsMomentum", "Pion track: TPC nSigma kaon contamination check", kTH2F, {axisMomentumForPID, axisTPCNSigma});
    allHistograms.add("QAafter/Pion/tpcNSigmaVsPt", "Pion TPC nSigma vs pT (after cuts)", kTH2F, {axisPtForPID, axisTPCNSigma});
    allHistograms.add("QAafter/Pion/tpcNSigmaVsMomentum", "Pion TPC nSigma vs total momentum (after cuts)", kTH2F, {axisMomentumForPID, axisTPCNSigma});
    allHistograms.add("QAafter/Pion/tpcNSigmaVsCentrality", "Pion TPC nSigma vs centrality", kTH2F, {axisCentralityPercent, axisTPCNSigma});
    allHistograms.add("QAafter/Pion/tofNSigmaVsPt", "Pion TOF nSigma vs pT (after cuts)", kTH2F, {axisPtForPID, axisTOFNSigma});
    allHistograms.add("QAafter/Pion/tofNSigmaVsMomentum", "Pion TOF nSigma vs total momentum (after cuts)", kTH2F, {axisMomentumForPID, axisTOFNSigma});
    allHistograms.add("QAafter/Pion/tofNSigmaVsCentrality", "Pion TOF nSigma vs centrality", kTH2F, {axisCentralityPercent, axisTOFNSigma});
    allHistograms.add("QAafter/Pion/tofNSigmaVsTPCNSigma", "Pion TOF nSigma vs TPC nSigma (after cuts)", kTH2F, {axisTPCNSigma, axisTOFNSigma});
    allHistograms.add("QAafter/Pion/tpcCrossedRowsVsPt", "Pion TPC crossed rows vs pT", kTH2F, {axisPtForPID, {200, 0, 200}});
    allHistograms.add("QAafter/Pion/tpcClustersFoundVsPt", "Pion TPC clusters found vs pT", kTH2F, {axisPtForPID, {200, 0, 200}});

    allHistograms.add("PairQA/hOpeningAngle", "Opening angle of accepted pairs", kTH1F, {axisOpenAngle});
    allHistograms.add("PairQA/hPairRapidity", "Rapidity of accepted pairs", kTH1F, {axisRapidity});
    allHistograms.add("PairQA/hPairPt", "p_{T} of accepted pairs", kTH1F, {axisPt});

    if (!doprocessMC) {
      allHistograms.add("THnSparse/invMass_DeltaPlusPlus", "Invariant mass: #Delta^{++}", kTHnSparseF, {axisInvariantMass, axisPt, axisCentralityPercent, axisRapidity, axisOccupancy});
      allHistograms.add("THnSparse/invMass_AntiDeltaPlusPlus", "Invariant mass: #bar{#Delta}^{++}", kTHnSparseF, {axisInvariantMass, axisPt, axisCentralityPercent, axisRapidity, axisOccupancy});
      allHistograms.add("THnSparse/invMass_DeltaZero", "Invariant mass: #Delta^{0}", kTHnSparseF, {axisInvariantMass, axisPt, axisCentralityPercent, axisRapidity, axisOccupancy});
      allHistograms.add("THnSparse/invMass_AntiDeltaZero", "Invariant mass: #bar{#Delta}^{0}", kTHnSparseF, {axisInvariantMass, axisPt, axisCentralityPercent, axisRapidity, axisOccupancy});

      allHistograms.add("THnSparse/invMass_Rotated_DeltaPlusPlus", "Rotational background: #Delta^{++}", kTHnSparseF, {axisInvariantMass, axisPt, axisCentralityPercent, axisRapidity, axisOccupancy});
      allHistograms.add("THnSparse/invMass_Rotated_AntiDeltaPlusPlus", "Rotational background: #bar{#Delta}^{++}", kTHnSparseF, {axisInvariantMass, axisPt, axisCentralityPercent, axisRapidity, axisOccupancy});
      allHistograms.add("THnSparse/invMass_Rotated_DeltaZero", "Rotational background: #Delta^{0}", kTHnSparseF, {axisInvariantMass, axisPt, axisCentralityPercent, axisRapidity, axisOccupancy});
      allHistograms.add("THnSparse/invMass_Rotated_AntiDeltaZero", "Rotational background: #bar{#Delta}^{0}", kTHnSparseF, {axisInvariantMass, axisPt, axisCentralityPercent, axisRapidity, axisOccupancy});

      if (doprocessMix || doprocessMixDF || doprocessMixepDF) {
        allHistograms.add("THnSparse/invMass_Mixed_DeltaPlusPlus", "Mixed background: #Delta^{++}", kTHnSparseF, {axisInvariantMass, axisPt, axisCentralityPercent, axisRapidity, axisOccupancy});
        allHistograms.add("THnSparse/invMass_Mixed_AntiDeltaPlusPlus", "Mixed background: #bar{#Delta}^{++}", kTHnSparseF, {axisInvariantMass, axisPt, axisCentralityPercent, axisRapidity, axisOccupancy});
        allHistograms.add("THnSparse/invMass_Mixed_DeltaZero", "Mixed background: #Delta^{0}", kTHnSparseF, {axisInvariantMass, axisPt, axisCentralityPercent, axisRapidity, axisOccupancy});
        allHistograms.add("THnSparse/invMass_Mixed_AntiDeltaZero", "Mixed background: #bar{#Delta}^{0}", kTHnSparseF, {axisInvariantMass, axisPt, axisCentralityPercent, axisRapidity, axisOccupancy});
      }
    }

    if (doprocessMC) {
      allHistograms.add("Event/mcEventSelectionCutflow", "MC event selection cutflow", kTH1F, {{7, 0, 7}});
      allHistograms.add("QAChecks/protonRecoLevelPt", "Reconstructed proton pT (after PID cuts) - MC", kTH1F, {axisPtForPID});
      allHistograms.add("QAChecks/pionRecoLevelPt", "Reconstructed pion pT (after PID cuts) - MC", kTH1F, {axisPtForPID});
      allHistograms.add("QAChecks/protonGenLevelPt", "Generated proton pT (truth level) - MC", kTH1F, {axisPtForPID});
      allHistograms.add("QAChecks/pionGenLevelPt", "Generated pion pT (truth level) - MC", kTH1F, {axisPtForPID});

      allHistograms.add("Analysis/mcGenerated_DeltaPlusPlus", "Generated Delta++ mass vs pT vs centrality", kTHnSparseF, {axisInvariantMass, axisPt, axisCentralityPercent});
      allHistograms.add("Analysis/mcGenerated_AntiDeltaPlusPlus", "Generated anti-Delta++ mass vs pT vs centrality", kTHnSparseF, {axisInvariantMass, axisPt, axisCentralityPercent});
      allHistograms.add("Analysis/mcGenerated_DeltaZero", "Generated Delta0 mass vs pT vs centrality", kTHnSparseF, {axisInvariantMass, axisPt, axisCentralityPercent});
      allHistograms.add("Analysis/mcGenerated_AntiDeltaZero", "Generated anti-Delta0 mass vs pT vs centrality", kTHnSparseF, {axisInvariantMass, axisPt, axisCentralityPercent});

      allHistograms.add("Analysis/mcReconstructed_DeltaPlusPlus", "Reconstructed Delta++ mass vs pT vs cent - MC", kTHnSparseF, {axisInvariantMass, axisPt, axisCentralityPercent});
      allHistograms.add("Analysis/mcReconstructed_AntiDeltaPlusPlus", "Reconstructed anti-Delta++ mass vs pT vs cent - MC", kTHnSparseF, {axisInvariantMass, axisPt, axisCentralityPercent});
      allHistograms.add("Analysis/mcReconstructed_DeltaZero", "Reconstructed Delta0 mass vs pT vs cent - MC", kTHnSparseF, {axisInvariantMass, axisPt, axisCentralityPercent});
      allHistograms.add("Analysis/mcReconstructed_AntiDeltaZero", "Reconstructed anti-Delta0 mass vs pT vs cent - MC", kTHnSparseF, {axisInvariantMass, axisPt, axisCentralityPercent});

      allHistograms.add("Analysis/mcMassResolution_DeltaPlusPlus", "Mass resolution vs pT - Delta++", kTHnSparseF, {{200, -0.05, 0.05}, axisPt, axisCentralityPercent});
      allHistograms.add("Analysis/mcMassResolution_AntiDeltaPlusPlus", "Mass resolution vs pT - anti-Delta++", kTHnSparseF, {{200, -0.05, 0.05}, axisPt, axisCentralityPercent});
      allHistograms.add("Analysis/mcMassResolution_DeltaZero", "Mass resolution vs pT - Delta0", kTHnSparseF, {{200, -0.05, 0.05}, axisPt, axisCentralityPercent});
      allHistograms.add("Analysis/mcMassResolution_AntiDeltaZero", "Mass resolution vs pT - anti-Delta0", kTHnSparseF, {{200, -0.05, 0.05}, axisPt, axisCentralityPercent});
    }

    if (doprocessMCGen) {
      allHistograms.add("SignalLoss/mcEventSelectionCutflow", "MC event selection cutflow (signal loss study)", kTH1F, {{7, 0, 7}});
      allHistograms.add("SignalLoss/mTScaled_fromDeltaPlusPlus", "mT-scaled Delta++ pT", kTHnSparseF, {axisPt, axisCentralityPercent});
      allHistograms.add("SignalLoss/mTScaled_fromAntiDeltaPlusPlus", "mT-scaled anti-Delta++ pT", kTHnSparseF, {axisPt, axisCentralityPercent});
      allHistograms.add("SignalLoss/mTScaled_fromDeltaZero", "mT-scaled Delta0 pT", kTHnSparseF, {axisPt, axisCentralityPercent});
      allHistograms.add("SignalLoss/mTScaled_fromAntiDeltaZero", "mT-scaled anti-Delta0 pT", kTHnSparseF, {axisPt, axisCentralityPercent});
    }
  }

  template <typename TrackType>
  bool passesBasicTrackSelection(TrackType const& track)
  {
    if (track.pt() < minTrackPt)
      return false;
    if (track.eta() < minPseudorapidity || track.eta() > maxPseudorapidity)
      return false;
    if (requirePrimaryTrack && !track.isPrimaryTrack())
      return false;
    if (requireGlobalTrackNoDCA && !track.isGlobalTrackWoDCA())
      return false;
    if (requirePVContributor && !track.isPVContributor())
      return false;
    if (applyCrossedRowsCut && track.tpcNClsCrossedRows() < minTPCCrossedRows)
      return false;
    if (applyTPCClustersCut && track.tpcNClsFound() < minTPCClustersFound)
      return false;
    return true;
  }

  template <typename TrackType>
  bool passesProtonDCASelection(TrackType const& track, float totalMomentum)
  {
    auto ptBinEdges = static_cast<std::vector<float>>(protonDCAPtBinEdges);
    auto maxDCAxyValues = static_cast<std::vector<float>>(protonMaxDCAxyPerPtBin);
    int numberOfBins = static_cast<int>(ptBinEdges.size()) - 1;
    bool dcaXYPassed = false;
    for (int iBin = 0; iBin < numberOfBins; iBin++) {
      if (totalMomentum >= ptBinEdges[iBin] && totalMomentum < ptBinEdges[iBin + 1] &&
          std::abs(track.dcaXY()) < maxDCAxyValues[iBin])
        dcaXYPassed = true;
    }
    if (!dcaXYPassed)
      return false;
    if (std::abs(track.dcaZ()) > maxDCAz)
      return false;
    return true;
  }

  template <typename TrackType>
  bool passesPionDCASelection(TrackType const& track, float totalMomentum)
  {
    auto ptBinEdges = static_cast<std::vector<float>>(pionDCAPtBinEdges);
    auto maxDCAxyValues = static_cast<std::vector<float>>(pionMaxDCAxyPerPtBin);
    int numberOfBins = static_cast<int>(ptBinEdges.size()) - 1;
    bool dcaXYPassed = false;
    for (int iBin = 0; iBin < numberOfBins; iBin++) {
      if (totalMomentum >= ptBinEdges[iBin] && totalMomentum < ptBinEdges[iBin + 1] &&
          std::abs(track.dcaXY()) < maxDCAxyValues[iBin])
        dcaXYPassed = true;
    }
    if (!dcaXYPassed)
      return false;
    if (std::abs(track.dcaZ()) > maxDCAz)
      return false;
    return true;
  }

  template <typename TrackType>
  bool passesProtonPID(TrackType const& track, float totalMomentum)
  {
    bool tpcPassed{false}, tofPassed{false};
    auto tpcMomBins = static_cast<std::vector<float>>(protonTPCPIDMomentumBins);
    auto tpcNSigCuts = static_cast<std::vector<float>>(protonTPCNSigmaCutPerBin);
    auto tofMomBins = static_cast<std::vector<float>>(protonTOFPIDMomentumBins);
    auto tofNSigCuts = static_cast<std::vector<float>>(protonTOFNSigmaCutPerBin);
    int nTPCBins = static_cast<int>(tpcMomBins.size());
    int nTOFBins = static_cast<int>(tofMomBins.size());

    float tpcNSigPi = std::abs(track.tpcNSigmaPi());
    float tpcNSigKaon = std::abs(track.tpcNSigmaKa());
    float tpcNSigPr = std::abs(track.tpcNSigmaPr());
    float tofNSigPi = std::abs(track.tofNSigmaPi());
    float tofNSigKaon = std::abs(track.tofNSigmaKa());
    float tofNSigPr = std::abs(track.tofNSigmaPr());

    float combinedNSigPi = tpcNSigPi * tpcNSigPi + tofNSigPi * tofNSigPi;
    float combinedNSigPr = tpcNSigPr * tpcNSigPr + tofNSigPr * tofNSigPr;
    float circularCutSq = static_cast<float>(combinedNSigmaCutProton * combinedNSigmaCutProton);
    float circularVetoCutSq = tpcNSigmaVetoThreshold * tpcNSigmaVetoThreshold + tofNSigmaVetoThreshold * tofNSigmaVetoThreshold;

    if (!useTPCOnlyPID && track.hasTOF()) {
      if (track.tofNSigmaPr() < minTOFNSigmaProton)
        return false;
      if (combinedNSigmaCutProton < 0 && totalMomentum >= minProtonMomentum) {
        for (int i = 0; i < nTOFBins - 1; ++i) {
          if (totalMomentum >= tofMomBins[i] && totalMomentum < tofMomBins[i + 1] &&
              (tofNSigPr < tofNSigCuts[i] && tofNSigPi > tofNSigmaVetoThreshold && tofNSigKaon > tofNSigmaVetoThreshold))
            tofPassed = true;
        }
        if (track.tpcNSigmaPr() < minCombinedTPCTOFNSigmaProton)
          return false;
        if (tpcNSigPr < static_cast<float>(maxTPCNSigmaProton) && tpcNSigPi > tpcNSigmaVetoThreshold && tpcNSigKaon > tpcNSigmaVetoThreshold)
          tpcPassed = true;
      }
      if ((combinedNSigmaCutProton > 0) && totalMomentum >= minProtonMomentum &&
          (combinedNSigPr < circularCutSq && combinedNSigPi > circularVetoCutSq)) {
        tofPassed = true;
        tpcPassed = true;
      }
      if (totalMomentum < minProtonMomentum && tpcNSigPr < static_cast<float>(maxTPCNSigmaProton)) {
        tofPassed = true;
        tpcPassed = true;
      }
    } else {
      tofPassed = true;
      if (track.tpcNSigmaPr() < minTPCNSigmaProton)
        return false;
      for (int i = 0; i < nTPCBins - 1; ++i) {
        if (totalMomentum >= tpcMomBins[i] && totalMomentum < tpcMomBins[i + 1] &&
            (tpcNSigPr < tpcNSigCuts[i] && tpcNSigPi > tpcNSigmaVetoOtherSpecies))
          tpcPassed = true;
      }
    }
    return (tpcPassed && tofPassed);
  }

  template <typename TrackType>
  bool passesPionPID(TrackType const& track, float totalMomentum)
  {
    bool tpcPassed{false}, tofPassed{false};
    auto tpcMomBins = static_cast<std::vector<float>>(pionTPCPIDMomentumBins);
    auto tpcNSigCuts = static_cast<std::vector<float>>(pionTPCNSigmaCutPerBin);
    auto tofMomBins = static_cast<std::vector<float>>(pionTOFPIDMomentumBins);
    auto tofNSigCuts = static_cast<std::vector<float>>(pionTOFNSigmaCutPerBin);
    int nTPCBins = static_cast<int>(tpcMomBins.size());
    int nTOFBins = static_cast<int>(tofMomBins.size());

    float tpcNSigPi = std::abs(track.tpcNSigmaPi());
    float tpcNSigKaon = std::abs(track.tpcNSigmaKa());
    float tpcNSigPr = std::abs(track.tpcNSigmaPr());
    float tofNSigPi = std::abs(track.tofNSigmaPi());
    float tofNSigKaon = std::abs(track.tofNSigmaKa());
    float tofNSigPr = std::abs(track.tofNSigmaPr());

    float combinedNSigPi = tpcNSigPi * tpcNSigPi + tofNSigPi * tofNSigPi;
    float combinedNSigPr = tpcNSigPr * tpcNSigPr + tofNSigPr * tofNSigPr;
    float circularCutSq = static_cast<float>(combinedNSigmaCutPion * combinedNSigmaCutPion);
    float circularVetoCutSq = tpcNSigmaVetoOtherSpecies * tofNSigmaVetoThreshold;

    if (!useTPCOnlyPID && track.hasTOF()) {
      if (track.tofNSigmaPi() < minTOFNSigmaPion)
        return false;
      if (combinedNSigmaCutPion < 0 && totalMomentum >= minPionMomentum) {
        for (int i = 0; i < nTOFBins - 1; ++i) {
          if (totalMomentum >= tofMomBins[i] && totalMomentum < tofMomBins[i + 1] &&
              (tofNSigPi < tofNSigCuts[i] && tofNSigPr > tofNSigmaVetoThreshold && tofNSigKaon > tofNSigmaVetoThreshold))
            tofPassed = true;
        }
        if (track.tpcNSigmaPi() < minCombinedTPCTOFNSigmaPion)
          return false;
        if (tpcNSigPi < static_cast<float>(maxTPCNSigmaPion) && tpcNSigPr > tpcNSigmaVetoThreshold && tpcNSigKaon > tpcNSigmaVetoThreshold)
          tpcPassed = true;
      }
      if ((combinedNSigmaCutPion > 0) && totalMomentum >= minPionMomentum &&
          (combinedNSigPi < circularCutSq && combinedNSigPr > circularVetoCutSq)) {
        tofPassed = true;
        tpcPassed = true;
      }
      if (totalMomentum < minPionMomentum && tpcNSigPi < static_cast<float>(maxTPCNSigmaPion)) {
        tofPassed = true;
        tpcPassed = true;
      }
    } else {
      tofPassed = true;
      if (track.tpcNSigmaPi() < minTPCNSigmaPion)
        return false;
      for (int i = 0; i < nTPCBins - 1; ++i) {
        if (totalMomentum >= tpcMomBins[i] && totalMomentum < tpcMomBins[i + 1] &&
            (tpcNSigPi < tpcNSigCuts[i] && tpcNSigPr > tpcNSigmaVetoOtherSpecies))
          tpcPassed = true;
      }
    }
    return (tpcPassed && tofPassed);
  }

  template <bool isMixedEvent, bool isMCAnalysis, typename TrackCollectionType>
  void fillInvariantMassHistograms(TrackCollectionType const& tracks1,
                                   TrackCollectionType const& tracks2,
                                   float centralityPercent,
                                   int occupancyValue = 100)
  {
    std::vector<ValidCandidate> validProtons;
    std::vector<ValidCandidate> validPions;
    validProtons.reserve(50);
    validPions.reserve(50);

    // Single-pass track pre-selection and QA filling logic
    auto processTracks = [&](const TrackCollectionType& trackCol, bool isCol1) {
      for (const auto& track : trackCol) {
        if (!passesBasicTrackSelection(track))
          continue;

        float mom = RecoDecay::p(track.px(), track.py(), track.pz());

        if (!isMixedEvent && isCol1) {
          allHistograms.fill(HIST("QAbefore/trackEta"), track.eta());
          allHistograms.fill(HIST("QAbefore/trackPt"), track.pt());
          allHistograms.fill(HIST("QAbefore/trackPhi"), track.phi());
          allHistograms.fill(HIST("QAbefore/trackEtaVsPhi"), track.eta(), track.phi());

          auto tpcNSigPr = track.tpcNSigmaPr();
          allHistograms.fill(HIST("QAbefore/Proton/tpcNSigmaVsMomentum"), mom, tpcNSigPr);
          allHistograms.fill(HIST("QAbefore/Proton/tpcDedxVsMomentum"), mom, track.tpcSignal());
          if (track.hasTOF()) {
            auto tofNSigPr = track.tofNSigmaPr();
            allHistograms.fill(HIST("QAbefore/Proton/tofNSigmaVsMomentum"), mom, tofNSigPr);
            allHistograms.fill(HIST("QAbefore/Proton/tofNSigmaVsTPCNSigma"), tpcNSigPr, tofNSigPr);
          }

          auto tpcNSigPi = track.tpcNSigmaPi();
          allHistograms.fill(HIST("QAbefore/Pion/tpcNSigmaVsMomentum"), mom, tpcNSigPi);
          allHistograms.fill(HIST("QAbefore/Pion/tpcDedxVsMomentum"), mom, track.tpcSignal());
          if (track.hasTOF()) {
            auto tofNSigPi = track.tofNSigmaPi();
            allHistograms.fill(HIST("QAbefore/Pion/tofNSigmaVsMomentum"), mom, tofNSigPi);
            allHistograms.fill(HIST("QAbefore/Pion/tofNSigmaVsTPCNSigma"), tpcNSigPi, tofNSigPi);
          }
        }

        bool isProton = passesProtonPID(track, mom) && passesProtonDCASelection(track, mom) && (!requireTOFForProton || track.hasTOF());
        bool isPion = passesPionPID(track, mom) && passesPionDCASelection(track, mom) && (!requireTOFForPion || track.hasTOF());

        if (isProton || isPion) {
          int motherId = -1, motherPDG = 0;
          if constexpr (isMCAnalysis) {
            motherId = track.motherId();
            motherPDG = track.motherPDG();
          }

          ValidCandidate cand{track.px(), track.py(), track.pz(), track.pt(), mom, track.phi(), track.eta(), track.sign(), motherId, motherPDG};

          if (isProton) {
            validProtons.push_back(cand);
            if (!isMixedEvent && isCol1) {
              allHistograms.fill(HIST("QAafter/Proton/ptVsCentrality"), track.pt(), centralityPercent);
              allHistograms.fill(HIST("QAafter/Proton/dcaZVsPt"), track.pt(), track.dcaZ());
              allHistograms.fill(HIST("QAafter/Proton/dcaXYVsPt"), track.pt(), track.dcaXY());
              allHistograms.fill(HIST("QAafter/Proton/tpcDedxVsMomentum"), mom, track.tpcSignal());
              allHistograms.fill(HIST("QAafter/Proton/tpcNSigmaPionContamVsPt"), track.pt(), track.tpcNSigmaPi());
              allHistograms.fill(HIST("QAafter/Proton/tpcNSigmaKaonContamVsPt"), track.pt(), track.tpcNSigmaKa());
              allHistograms.fill(HIST("QAafter/Proton/tpcNSigmaVsPt"), track.pt(), track.tpcNSigmaPr());
              allHistograms.fill(HIST("QAafter/Proton/tpcNSigmaVsMomentum"), mom, track.tpcNSigmaPr());
              allHistograms.fill(HIST("QAafter/Proton/tpcNSigmaVsCentrality"), centralityPercent, track.tpcNSigmaPr());
              allHistograms.fill(HIST("QAafter/Proton/tpcCrossedRowsVsPt"), track.pt(), track.tpcNClsCrossedRows());
              allHistograms.fill(HIST("QAafter/Proton/tpcClustersFoundVsPt"), track.pt(), track.tpcNClsFound());
              if (!useTPCOnlyPID && track.hasTOF()) {
                allHistograms.fill(HIST("QAafter/Proton/tofNSigmaVsPt"), track.pt(), track.tofNSigmaPr());
                allHistograms.fill(HIST("QAafter/Proton/tofNSigmaVsMomentum"), mom, track.tofNSigmaPr());
                allHistograms.fill(HIST("QAafter/Proton/tofNSigmaVsCentrality"), centralityPercent, track.tofNSigmaPr());
                allHistograms.fill(HIST("QAafter/Proton/tofNSigmaVsTPCNSigma"), track.tpcNSigmaPr(), track.tofNSigmaPr());
              }
            }
          }
          if (isPion) {
            validPions.push_back(cand);
            if (!isMixedEvent && isCol1) {
              allHistograms.fill(HIST("QAafter/Pion/ptVsCentrality"), track.pt(), centralityPercent);
              allHistograms.fill(HIST("QAafter/Pion/dcaZVsPt"), track.pt(), track.dcaZ());
              allHistograms.fill(HIST("QAafter/Pion/dcaXYVsPt"), track.pt(), track.dcaXY());
              allHistograms.fill(HIST("QAafter/Pion/tpcDedxVsMomentum"), mom, track.tpcSignal());
              allHistograms.fill(HIST("QAafter/Pion/tpcNSigmaProtonContamVsMomentum"), mom, track.tpcNSigmaPr());
              allHistograms.fill(HIST("QAafter/Pion/tpcNSigmaKaonContamVsMomentum"), mom, track.tpcNSigmaKa());
              allHistograms.fill(HIST("QAafter/Pion/tpcNSigmaVsPt"), track.pt(), track.tpcNSigmaPi());
              allHistograms.fill(HIST("QAafter/Pion/tpcNSigmaVsMomentum"), mom, track.tpcNSigmaPi());
              allHistograms.fill(HIST("QAafter/Pion/tpcNSigmaVsCentrality"), centralityPercent, track.tpcNSigmaPi());
              allHistograms.fill(HIST("QAafter/Pion/tpcCrossedRowsVsPt"), track.pt(), track.tpcNClsCrossedRows());
              allHistograms.fill(HIST("QAafter/Pion/tpcClustersFoundVsPt"), track.pt(), track.tpcNClsFound());
              if (!useTPCOnlyPID && track.hasTOF()) {
                allHistograms.fill(HIST("QAafter/Pion/tofNSigmaVsPt"), track.pt(), track.tofNSigmaPi());
                allHistograms.fill(HIST("QAafter/Pion/tofNSigmaVsMomentum"), mom, track.tofNSigmaPi());
                allHistograms.fill(HIST("QAafter/Pion/tofNSigmaVsCentrality"), centralityPercent, track.tofNSigmaPi());
                allHistograms.fill(HIST("QAafter/Pion/tofNSigmaVsTPCNSigma"), track.tpcNSigmaPi(), track.tofNSigmaPi());
              }
            }
          }
        }
      }
    };

    processTracks(tracks1, true);
    if constexpr (isMixedEvent) {
      processTracks(tracks2, false);
    }

    if (runQualityChecksOnly)
      return;

    for (size_t i = 0; i < validProtons.size(); ++i) {
      const auto& prCand = validProtons[i];
      for (size_t j = 0; j < validPions.size(); ++j) {
        const auto& piCand = validPions[j];

        // Safely avoid duplicate indexing evaluation dynamically
        if constexpr (!isMixedEvent) {
          if (&prCand == &piCand)
            continue;
        }

        float cosAngle = std::clamp((prCand.px * piCand.px + prCand.py * piCand.py + prCand.pz * piCand.pz) / (prCand.p * piCand.p), -1.f, 1.f);
        float openAngle = std::acos(cosAngle);

        if (applyDeepAngleCut && openAngle < deepAngleCutValue)
          continue;

        std::array<std::array<float, 3>, 2> momenta = {
          std::array<float, 3>{prCand.px, prCand.py, prCand.pz},
          std::array<float, 3>{piCand.px, piCand.py, piCand.pz}};

        float pairMass = RecoDecay::m(momenta, std::array{delta_analysis::MassProton, delta_analysis::MassPion});
        float pairPt = RecoDecay::pt(std::array{prCand.px + piCand.px, prCand.py + piCand.py});
        float pairY = RecoDecay::y(std::array{prCand.px + piCand.px, prCand.py + piCand.py, prCand.pz + piCand.pz}, pairMass);

        if (pairY < minPairRapidity || pairY > maxPairRapidity)
          continue;

        if constexpr (!isMixedEvent && !isMCAnalysis) {
          allHistograms.fill(HIST("PairQA/hOpeningAngle"), openAngle);
          allHistograms.fill(HIST("PairQA/hPairRapidity"), pairY);
          allHistograms.fill(HIST("PairQA/hPairPt"), pairPt);

          if (prCand.sign > 0 && piCand.sign > 0) {
            allHistograms.fill(HIST("THnSparse/invMass_DeltaPlusPlus"), pairMass, pairPt, centralityPercent, pairY, occupancyValue);
          } else if (prCand.sign < 0 && piCand.sign < 0) {
            allHistograms.fill(HIST("THnSparse/invMass_AntiDeltaPlusPlus"), pairMass, pairPt, centralityPercent, pairY, occupancyValue);
          } else if (prCand.sign > 0 && piCand.sign < 0) {
            allHistograms.fill(HIST("THnSparse/invMass_DeltaZero"), pairMass, pairPt, centralityPercent, pairY, occupancyValue);
          } else if (prCand.sign < 0 && piCand.sign > 0) {
            allHistograms.fill(HIST("THnSparse/invMass_AntiDeltaZero"), pairMass, pairPt, centralityPercent, pairY, occupancyValue);
          }

          if (enableRotationalBackground) {
            for (const auto& rotOffset : precomputedRotations) {
              float newPionPhi = RecoDecay::constrainAngle(piCand.phi + rotOffset, 0.f);
              float pxPionRotated = piCand.pt * std::cos(newPionPhi);
              float pyPionRotated = piCand.pt * std::sin(newPionPhi);

              std::array<std::array<float, 3>, 2> rotatedMomenta = {
                std::array<float, 3>{prCand.px, prCand.py, prCand.pz},
                std::array<float, 3>{pxPionRotated, pyPionRotated, piCand.pz}};

              float rotatedPairInvMass = RecoDecay::m(rotatedMomenta, std::array{delta_analysis::MassProton, delta_analysis::MassPion});
              float rotatedPairPt = RecoDecay::pt(std::array{prCand.px + pxPionRotated, prCand.py + pyPionRotated});
              float rotatedPairRapidity = RecoDecay::y(std::array{prCand.px + pxPionRotated, prCand.py + pyPionRotated, prCand.pz + piCand.pz}, rotatedPairInvMass);

              if (rotatedPairRapidity < minPairRapidity || rotatedPairRapidity > maxPairRapidity)
                continue;

              if (prCand.sign > 0 && piCand.sign > 0) {
                allHistograms.fill(HIST("THnSparse/invMass_Rotated_DeltaPlusPlus"), rotatedPairInvMass, rotatedPairPt, centralityPercent, rotatedPairRapidity, occupancyValue);
              } else if (prCand.sign < 0 && piCand.sign < 0) {
                allHistograms.fill(HIST("THnSparse/invMass_Rotated_AntiDeltaPlusPlus"), rotatedPairInvMass, rotatedPairPt, centralityPercent, rotatedPairRapidity, occupancyValue);
              } else if (prCand.sign > 0 && piCand.sign < 0) {
                allHistograms.fill(HIST("THnSparse/invMass_Rotated_DeltaZero"), rotatedPairInvMass, rotatedPairPt, centralityPercent, rotatedPairRapidity, occupancyValue);
              } else if (prCand.sign < 0 && piCand.sign > 0) {
                allHistograms.fill(HIST("THnSparse/invMass_Rotated_AntiDeltaZero"), rotatedPairInvMass, rotatedPairPt, centralityPercent, rotatedPairRapidity, occupancyValue);
              }
            }
          }
        }

        if constexpr (isMixedEvent) {
          if (prCand.sign > 0 && piCand.sign > 0) {
            allHistograms.fill(HIST("THnSparse/invMass_Mixed_DeltaPlusPlus"), pairMass, pairPt, centralityPercent, pairY, occupancyValue);
          } else if (prCand.sign < 0 && piCand.sign < 0) {
            allHistograms.fill(HIST("THnSparse/invMass_Mixed_AntiDeltaPlusPlus"), pairMass, pairPt, centralityPercent, pairY, occupancyValue);
          } else if (prCand.sign > 0 && piCand.sign < 0) {
            allHistograms.fill(HIST("THnSparse/invMass_Mixed_DeltaZero"), pairMass, pairPt, centralityPercent, pairY, occupancyValue);
          } else if (prCand.sign < 0 && piCand.sign > 0) {
            allHistograms.fill(HIST("THnSparse/invMass_Mixed_AntiDeltaZero"), pairMass, pairPt, centralityPercent, pairY, occupancyValue);
          }
        }

        if constexpr (isMCAnalysis) {
          if (prCand.motherId != piCand.motherId || prCand.motherId == -1)
            continue;
          if (prCand.motherPDG != piCand.motherPDG || prCand.motherPDG == -1)
            continue;

          int mPDG = std::abs(prCand.motherPDG);
          if (mPDG != delta_analysis::PdgDeltaPlusPlus && mPDG != delta_analysis::PdgDeltaZero)
            continue;

          float trueMassOfParent = 0.f;
          for (auto const& parentParticle : *mcResonanceParentTable) {
            if (parentParticle.mcParticleId() == prCand.motherId) {
              float tE = parentParticle.e(), tPx = parentParticle.px(), tPy = parentParticle.py(), tPz = parentParticle.pz();
              float massSq = tE * tE - tPx * tPx - tPy * tPy - tPz * tPz;
              trueMassOfParent = massSq > 0.f ? std::sqrt(massSq) : 0.f;
              break;
            }
          }
          float massResidual = pairMass - trueMassOfParent;

          if (prCand.sign > 0 && piCand.sign > 0) {
            allHistograms.fill(HIST("Analysis/mcReconstructed_DeltaPlusPlus"), pairMass, pairPt, centralityPercent);
            allHistograms.fill(HIST("Analysis/mcMassResolution_DeltaPlusPlus"), massResidual, pairPt, centralityPercent);
          } else if (prCand.sign < 0 && piCand.sign < 0) {
            allHistograms.fill(HIST("Analysis/mcReconstructed_AntiDeltaPlusPlus"), pairMass, pairPt, centralityPercent);
            allHistograms.fill(HIST("Analysis/mcMassResolution_AntiDeltaPlusPlus"), massResidual, pairPt, centralityPercent);
          } else if (prCand.sign > 0 && piCand.sign < 0) {
            allHistograms.fill(HIST("Analysis/mcReconstructed_DeltaZero"), pairMass, pairPt, centralityPercent);
            allHistograms.fill(HIST("Analysis/mcMassResolution_DeltaZero"), massResidual, pairPt, centralityPercent);
          } else if (prCand.sign < 0 && piCand.sign > 0) {
            allHistograms.fill(HIST("Analysis/mcReconstructed_AntiDeltaZero"), pairMass, pairPt, centralityPercent);
            allHistograms.fill(HIST("Analysis/mcMassResolution_AntiDeltaZero"), massResidual, pairPt, centralityPercent);
          }
        }
      }
    }
  }

  // Proper joins natively compatible with O2 evsels.
  using ResonanceCollisionsWithEP = soa::Join<aod::ResoCollisions, aod::ResoEvtPlCollisions>;
  using ResonanceMCCollisions = soa::Join<aod::ResoCollisions, aod::ResoMCCollisions>;
  using MergedDFCollisions = aod::ResoCollisionDFs;
  using ResonanceTrackTable = aod::ResoTracks;

  void processData(ResonanceCollisionsWithEP::iterator const& collision, ResonanceTrackTable const& tracks)
  {
    if (cEvtRecINELgt0 && !collision.isRecINELgt0())
      return;

    allHistograms.fill(HIST("Event/centralityVsOccupancy"), collision.cent(), 100);
    allHistograms.fill(HIST("Event/primaryVertexZ"), collision.posZ());

    auto colTracks = tracks.sliceBy(tracksPerResonanceCollision, collision.globalIndex());
    fillInvariantMassHistograms<false, false>(colTracks, colTracks, collision.cent());
  }
  PROCESS_SWITCH(DeltaAnalysis1232, processData, "Process real collision data (same-event analysis)", true);

  using MixingBinningVtxZAndCentralityAndOcc = ColumnBinningPolicy<aod::collision::PosZ, aod::resocollision::Cent>;

  void processMix(ResonanceCollisionsWithEP const& collisions, ResonanceTrackTable const& tracks)
  {
    MixingBinningVtxZAndCentralityAndOcc mixingBins{{vertexZMixingBins, centralityMixingBins}, true};
    auto trackPool = std::make_tuple(tracks);
    SameKindPair<ResonanceCollisionsWithEP, ResonanceTrackTable, MixingBinningVtxZAndCentralityAndOcc>
      eventPairs{mixingBins, numberOfEventsToMix, -1, collisions, trackPool, &sliceCache};

    for (const auto& [col1, tracks1, col2, tracks2] : eventPairs) {
      if (cEvtRecINELgt0 && !col1.isRecINELgt0())
        return;

      allHistograms.fill(HIST("Event/mixingBins_centralityVsVtxZVsEventPlane"), col1.cent(), col1.posZ(), col1.evtPl());
      fillInvariantMassHistograms<true, false>(tracks1, tracks2, col1.cent());
    }
  }
  PROCESS_SWITCH(DeltaAnalysis1232, processMix, "Event mixing for background estimation", true);

  using MergedDFTracks = aod::ResoTrackDFs;

  void processDatadf(MergedDFCollisions::iterator const& collision, MergedDFTracks const& tracks)
  {
    if (doprocessData)
      LOG(error) << "Disable processData() first when using processDatadf()!";
    if (cEvtRecINELgt0 && !collision.isRecINELgt0())
      return;
    auto occupancyValue = 100;
    if (applyOccupancyInTimeRangeCut)
      occupancyValue = collision.trackOccupancyInTimeRange();
    if (occupancyValue < occupancyMin || occupancyValue > occupancyMax)
      return;

    allHistograms.fill(HIST("Event/centralityVsOccupancy"), collision.cent(), occupancyValue);
    allHistograms.fill(HIST("Event/primaryVertexZ"), collision.posZ());

    auto colTracks = tracks.sliceBy(tracksPerMergedDFCollision, collision.globalIndex());
    fillInvariantMassHistograms<false, false>(colTracks, colTracks, collision.cent(), occupancyValue);
  }
  PROCESS_SWITCH(DeltaAnalysis1232, processDatadf, "Process real data in merged derived-data (DF) format", false);

  using MixingBinningDF = ColumnBinningPolicy<aod::collision::PosZ, aod::resocollision::Cent>;

  void processMixDF(MergedDFCollisions const& collisions, MergedDFTracks const& tracks)
  {
    if (doprocessMix)
      LOG(fatal) << "Disable processMix() first when using processMixDF()!";
    MixingBinningDF mixingBins{{vertexZMixingBins, centralityMixingBins}, true};
    auto trackPool = std::make_tuple(tracks);
    SameKindPair<MergedDFCollisions, MergedDFTracks, MixingBinningDF>
      eventPairs{mixingBins, numberOfEventsToMix, -1, collisions, trackPool, &sliceCache};

    for (const auto& [col1, tracks1, col2, tracks2] : eventPairs) {
      if (cEvtRecINELgt0 && !col1.isRecINELgt0())
        return;
      auto occupancyValue = 100;
      if (applyOccupancyInTimeRangeCut)
        occupancyValue = col1.trackOccupancyInTimeRange();
      allHistograms.fill(HIST("Event/mixingBins_centralityVsVtxZVsEventPlane"), col1.cent(), col1.posZ(), col1.evtPl());
      fillInvariantMassHistograms<true, false>(tracks1, tracks2, col1.cent(), occupancyValue);
    }
  }
  PROCESS_SWITCH(DeltaAnalysis1232, processMixDF, "Event mixing for DF-format data", false);

  using MixingBinningWithEventPlane = ColumnBinningPolicy<aod::collision::PosZ, aod::resocollision::Cent, aod::resocollision::EvtPl>;

  void processMixepDF(MergedDFCollisions const& collisions, MergedDFTracks const& tracks)
  {
    if (doprocessMix || doprocessMixDF)
      LOG(fatal) << "Disable processMix() or processMixDF() first!";
    MixingBinningWithEventPlane mixingBins{{vertexZMixingBins, centralityMixingBins, eventPlaneMixingBins}, true};
    auto trackPool = std::make_tuple(tracks);
    SameKindPair<MergedDFCollisions, MergedDFTracks, MixingBinningWithEventPlane>
      eventPairs{mixingBins, numberOfEventsToMix, -1, collisions, trackPool, &sliceCache};

    for (const auto& [col1, tracks1, col2, tracks2] : eventPairs) {
      allHistograms.fill(HIST("Event/mixingBins_centralityVsVtxZVsEventPlane"), col1.cent(), col1.posZ(), col1.evtPl());
      fillInvariantMassHistograms<true, false>(tracks1, tracks2, col1.cent());
    }
  }
  PROCESS_SWITCH(DeltaAnalysis1232, processMixepDF, "Event-plane dependent event mixing for DF-format data", false);

  void processMC(ResonanceMCCollisions::iterator const& collision, soa::Join<aod::ResoTracks, aod::ResoMCTracks> const& tracks, aod::ResoMCParents const& mcParents)
  {
    allHistograms.fill(HIST("Event/mcEventSelectionCutflow"), 0);
    if (mcRequireTriggerTVX && !collision.isTriggerTVX())
      return;
    allHistograms.fill(HIST("Event/mcEventSelectionCutflow"), 1);
    if (mcRequireVtxWithin10cm && !collision.isVtxIn10())
      return;
    allHistograms.fill(HIST("Event/mcEventSelectionCutflow"), 2);
    if (mcRequireINELgt0 && !collision.isINELgt0())
      return;
    allHistograms.fill(HIST("Event/mcEventSelectionCutflow"), 3);
    if (mcRequireSel8 && !collision.isInSel8())
      return;
    allHistograms.fill(HIST("Event/mcEventSelectionCutflow"), 4);
    if (cEvtRecINELgt0 && !collision.isRecINELgt0())
      return;
    allHistograms.fill(HIST("Event/mcEventSelectionCutflow"), 5);
    if (mcRequireAfterAllCuts && !collision.isInAfterAllCuts())
      return;
    allHistograms.fill(HIST("Event/mcEventSelectionCutflow"), 6);

    float centralityPercent = collision.cent();
    allHistograms.fill(HIST("Event/centralityVsOccupancy"), centralityPercent, 100);
    allHistograms.fill(HIST("Event/primaryVertexZ"), collision.posZ());

    auto colTracks = tracks.sliceBy(tracksPerResonanceCollision, collision.globalIndex());
    mcResonanceParentTable = &mcParents;

    fillInvariantMassHistograms<false, true>(colTracks, colTracks, centralityPercent);

    for (const auto& t0 : colTracks) {
      if (!passesBasicTrackSelection(t0))
        continue;
      float mom = RecoDecay::p(t0.px(), t0.py(), t0.pz());

      if (passesProtonPID(t0, mom) && std::abs(t0.pdgCode()) == delta_analysis::PdgProton) {
        allHistograms.fill(HIST("QAChecks/protonRecoLevelPt"), t0.pt());
      }
      if (passesPionPID(t0, mom) && std::abs(t0.pdgCode()) == delta_analysis::PdgPion) {
        allHistograms.fill(HIST("QAChecks/pionRecoLevelPt"), t0.pt());
      }
    }

    for (const auto& parentParticle : mcParents) {
      int pdg = std::abs(parentParticle.pdgCode());
      if (pdg != delta_analysis::PdgDeltaPlusPlus && pdg != delta_analysis::PdgDeltaZero)
        continue;
      if (parentParticle.y() < minPairRapidity || parentParticle.y() > maxPairRapidity)
        continue;

      bool hasPr = (std::abs(parentParticle.daughterPDG1()) == delta_analysis::PdgProton || std::abs(parentParticle.daughterPDG2()) == delta_analysis::PdgProton);
      bool hasPi = (std::abs(parentParticle.daughterPDG1()) == delta_analysis::PdgPion || std::abs(parentParticle.daughterPDG2()) == delta_analysis::PdgPion);

      if (!hasPr || !hasPi)
        continue;

      float tE = parentParticle.e(), tPx = parentParticle.px(), tPy = parentParticle.py(), tPz = parentParticle.pz();
      float massSq = tE * tE - tPx * tPx - tPy * tPy - tPz * tPz;
      float parentMass = massSq > 0.f ? std::sqrt(massSq) : 0.f;

      if (parentParticle.pdgCode() == delta_analysis::PdgDeltaPlusPlus) {
        allHistograms.fill(HIST("Analysis/mcGenerated_DeltaPlusPlus"), parentMass, parentParticle.pt(), centralityPercent);
      } else if (parentParticle.pdgCode() == -delta_analysis::PdgDeltaPlusPlus) {
        allHistograms.fill(HIST("Analysis/mcGenerated_AntiDeltaPlusPlus"), parentMass, parentParticle.pt(), centralityPercent);
      } else if (parentParticle.pdgCode() == delta_analysis::PdgDeltaZero) {
        allHistograms.fill(HIST("Analysis/mcGenerated_DeltaZero"), parentMass, parentParticle.pt(), centralityPercent);
      } else if (parentParticle.pdgCode() == -delta_analysis::PdgDeltaZero) {
        allHistograms.fill(HIST("Analysis/mcGenerated_AntiDeltaZero"), parentMass, parentParticle.pt(), centralityPercent);
      }
    }
  }
  PROCESS_SWITCH(DeltaAnalysis1232, processMC, "Process Monte Carlo simulated events", false);

  void processMCGen(ResonanceMCCollisions::iterator const& collision, aod::ResoMCParents const& mcParents)
  {
    allHistograms.fill(HIST("SignalLoss/mcEventSelectionCutflow"), 0);
    if (mcRequireTriggerTVX && !collision.isTriggerTVX())
      return;
    allHistograms.fill(HIST("SignalLoss/mcEventSelectionCutflow"), 1);
    if (mcRequireVtxWithin10cm && !collision.isVtxIn10())
      return;
    allHistograms.fill(HIST("SignalLoss/mcEventSelectionCutflow"), 2);
    if (mcRequireINELgt0 && !collision.isINELgt0())
      return;
    allHistograms.fill(HIST("SignalLoss/mcEventSelectionCutflow"), 3);
    if (mcRequireSel8 && !collision.isInSel8())
      return;
    allHistograms.fill(HIST("SignalLoss/mcEventSelectionCutflow"), 4);
    if (cEvtRecINELgt0 && !collision.isRecINELgt0())
      return;
    allHistograms.fill(HIST("SignalLoss/mcEventSelectionCutflow"), 5);
    if (mcRequireAfterAllCuts && !collision.isInAfterAllCuts())
      return;
    allHistograms.fill(HIST("SignalLoss/mcEventSelectionCutflow"), 6);

    float centralityPercent = collision.cent();

    for (const auto& parentParticle : mcParents) {
      if (parentParticle.y() < minPairRapidity || parentParticle.y() > maxPairRapidity)
        continue;

      int pdgCode = parentParticle.pdgCode();
      float parentPt = parentParticle.pt();
      float tE = parentParticle.e(), tPx = parentParticle.px(), tPy = parentParticle.py(), tPz = parentParticle.pz();
      float massSq = tE * tE - tPx * tPx - tPy * tPy - tPz * tPz;
      float parentMass = massSq > 0.f ? std::sqrt(massSq) : 0.f;

      auto computeMtScaledPtSquared = [&]() -> double {
        return (parentPt * parentPt) + (parentMass * parentMass) - (delta_analysis::MassDelta * delta_analysis::MassDelta);
      };

      double mTScaledPtSquared = -1.0;

      if (pdgCode == delta_analysis::PdgDeltaPlusPlus) {
        mTScaledPtSquared = computeMtScaledPtSquared();
        if (mTScaledPtSquared > 0)
          allHistograms.fill(HIST("SignalLoss/mTScaled_fromDeltaPlusPlus"), std::sqrt(mTScaledPtSquared), centralityPercent);
      } else if (pdgCode == -delta_analysis::PdgDeltaPlusPlus) {
        mTScaledPtSquared = computeMtScaledPtSquared();
        if (mTScaledPtSquared > 0)
          allHistograms.fill(HIST("SignalLoss/mTScaled_fromAntiDeltaPlusPlus"), std::sqrt(mTScaledPtSquared), centralityPercent);
      } else if (pdgCode == delta_analysis::PdgDeltaZero) {
        mTScaledPtSquared = computeMtScaledPtSquared();
        if (mTScaledPtSquared > 0)
          allHistograms.fill(HIST("SignalLoss/mTScaled_fromDeltaZero"), std::sqrt(mTScaledPtSquared), centralityPercent);
      } else if (pdgCode == -delta_analysis::PdgDeltaZero) {
        mTScaledPtSquared = computeMtScaledPtSquared();
        if (mTScaledPtSquared > 0)
          allHistograms.fill(HIST("SignalLoss/mTScaled_fromAntiDeltaZero"), std::sqrt(mTScaledPtSquared), centralityPercent);
      }
    }
  }
  PROCESS_SWITCH(DeltaAnalysis1232, processMCGen, "Generator-level MC signal loss study (mT scaling)", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<DeltaAnalysis1232>(cfgc)};
}
